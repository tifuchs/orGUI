"""Where a mapped frame's wall time goes, inside the real scheduler.

``benchmark_reconstruction_pipeline.py`` reports one number for the whole
mapping phase. This runs the same pipeline and splits that number by stage
-- HDF5 read, whole-frame correction, per-tile assembly, the native call,
the cross-tile merge, the checkpoint route -- and reports, for each of
them, *wall* time and *on-core* time separately.

The difference between the two is the point. A stage's wall time is what
it costs the thread running it; its on-core time is what it costs the
machine. Wall minus on-core is time the thread existed but was not
executing: waiting for the disk, waiting for the GIL, or preempted. At
depth 0 the mapping pipeline uses 7-8 of 24 cores, and which of those
three explains the gap decides whether there is anything to win.

On-core time comes from ``QueryThreadCycleTime``, not ``time.thread_time``:
the latter is ``GetThreadTimes`` on Windows and quantises to the ~15.6 ms
scheduler tick, which is larger than most of the stages measured here.
Cycles are converted with a rate calibrated by spinning at start-up.

Every stage is measured by wrapping an existing call site from outside --
no instrumentation is compiled into the pipeline, and the arithmetic is
untouched. Stages nest (the kernel, merge and route happen inside
``_map_frame_group``), so the table reports the nesting and derives the
per-tile assembly cost, which is the one stage with no call of its own,
as its parent's residual.

**Run this with ``PYTHONPATH`` pointing at the checkout**, for the reason
the pipeline benchmark gives; the resolved extension path is printed.
"""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import replace
import json
from pathlib import Path
import shutil
import tempfile
import threading
from time import perf_counter

import numpy as np

from benchmark_reconstruction_pipeline import (
    _summarize_checkpoints,
    _window_ranges,
)

from orgui.datautils.xrayutils import reconstruction
from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointRouter,
    _native_module,
)
import orgui.reconstruction_job as reconstruction_job
from orgui.reconstruction_job import (
    _base_provenance,
    _correction_pipeline,
    _execution_layout,
    _load_assets,
    _map_pending_ranges,
    read_job,
    resolve_work_block_pixels,
    split_memory_budget,
)


_kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
_kernel32.QueryThreadCycleTime.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_ulonglong),
]
_kernel32.QueryThreadCycleTime.restype = ctypes.c_int
_PSEUDO_CURRENT_THREAD = ctypes.c_void_p(_kernel32.GetCurrentThread())


def thread_cycles():
    """Cycles this thread has executed, including stall cycles.

    Unlike ``time.thread_time()`` (``GetThreadTimes``, quantised to the
    scheduler tick) this is exact enough to attribute a single millisecond
    stage, and unlike wall time it does not count a thread that is waiting
    for the GIL or for the disk.

    :rtype: int
    """
    value = ctypes.c_ulonglong()
    if not _kernel32.QueryThreadCycleTime(
        _PSEUDO_CURRENT_THREAD, ctypes.byref(value)
    ):
        raise OSError(ctypes.get_last_error(), "QueryThreadCycleTime failed")
    return value.value


def calibrate_cycle_rate(bursts=3, seconds=0.05):
    """Cycles per second of on-core time, measured by spinning.

    Taken as the maximum over several bursts: a burst that happened to be
    preempted measures too low, never too high.

    :rtype: float
    """
    best = 0.0
    for _burst in range(bursts):
        start_wall = perf_counter()
        start_cycles = thread_cycles()
        while perf_counter() - start_wall < seconds:
            pass
        elapsed = perf_counter() - start_wall
        best = max(best, (thread_cycles() - start_cycles) / elapsed)
    return best


_registry: list[dict] = []
_registry_lock = threading.Lock()
_local = threading.local()


def _bucket(name):
    """This thread's ``[wall, cycles, calls]`` accumulator for a stage.

    Per thread and merged at the end, so the hot path takes no lock: a
    shared counter would contend on exactly the resource being measured.
    """
    buckets = getattr(_local, "buckets", None)
    if buckets is None:
        buckets = _local.buckets = {}
        with _registry_lock:
            _registry.append(buckets)
    bucket = buckets.get(name)
    if bucket is None:
        bucket = buckets[name] = [0.0, 0, 0]
    return bucket


class timed:
    """Context manager accumulating wall time and cycles for one stage.

    A class rather than ``@contextmanager`` to keep the per-entry cost at
    two attribute lookups; the stages measured here are milliseconds, but
    ``np.stack`` is called several times per tile.
    """

    __slots__ = ("_bucket", "_wall", "_cycles")

    def __init__(self, name):
        self._bucket = _bucket(name)

    def __enter__(self):
        self._wall = perf_counter()
        self._cycles = thread_cycles()
        return self

    def __exit__(self, *exception):
        bucket = self._bucket
        bucket[0] += perf_counter() - self._wall
        bucket[1] += thread_cycles() - self._cycles
        bucket[2] += 1
        return False


def merged_stages():
    """Every thread's buckets, summed per stage.

    :rtype: dict[str, dict]
    """
    totals: dict[str, list] = {}
    with _registry_lock:
        snapshots = [dict(buckets) for buckets in _registry]
    for buckets in snapshots:
        for name, (wall, cycles, calls) in buckets.items():
            total = totals.setdefault(name, [0.0, 0, 0])
            total[0] += wall
            total[1] += cycles
            total[2] += calls
    return totals


def minimal_correction(correct, static_mask):
    """A correction that keeps the mask and drops the scaling arithmetic.

    A ceiling probe, not a proposal. Fusing correction into the native
    extension would remove the full-detector NumPy passes that apply the
    static factor, the exposure and the monitors -- it would not remove
    the raw copy, the repair, or the finiteness check, and it cannot make
    the stage cost less than nothing. So the wall-time difference between
    this and the real correction bounds from above what fusing could
    possibly buy, and it can be measured without writing the fused kernel.

    Which voxels a pixel reaches depends on geometry and the mask, never
    on the recorded intensity, so this arm routes the identical records
    and keeps the checkpoint fingerprint comparable; only the intensities
    and variances it carries are wrong, which is why it is a probe and
    nothing else.

    :param correct:
        The real pipeline, for its repair plan.
    :param static_mask:
        The job's static mask, or ``None``.
    """
    repair_plan = getattr(correct, "repair_plan", None)

    def correct_frame(payload, raw, frame_index):
        source_variance = getattr(payload, "variance", None)
        if source_variance is None:
            variance = np.maximum(raw, 0.0)
        else:
            variance = np.asarray(source_variance, dtype=np.float64).copy()
        intensity = np.asarray(raw, dtype=np.float64).copy()
        if static_mask is None:
            mask = np.zeros(intensity.shape, dtype=bool)
        else:
            mask = static_mask.copy()
        if repair_plan is not None:
            mask, _ = repair_plan.apply_inplace(intensity, variance)
        mask |= ~np.isfinite(intensity) | ~np.isfinite(variance)
        return intensity, variance, mask

    def pipeline(payload, raw, frame_index, tile):
        row_start, row_stop, column_start, column_stop = tile
        full = np.asarray(getattr(payload, "img", raw))
        intensity, variance, mask = correct_frame(payload, full, frame_index)
        selection = np.s_[row_start:row_stop, column_start:column_stop]
        return intensity[selection], variance[selection], mask[selection]

    pipeline.correct_frame = correct_frame
    pipeline.repair_plan = repair_plan
    return pipeline


class GilProbe(threading.Thread):
    """What share of its wall time a Python thread spends waiting for the GIL.

    Every stage's ``blocked`` column mixes GIL waiting with disk waiting
    and preemption, and no arithmetic separates them. This does: it runs
    bursts of pure Python -- no I/O, no NumPy, nothing that could release
    the GIL voluntarily -- and compares each burst's wall time against the
    cycles it actually executed. A burst is deliberately several times the
    interpreter's switch interval, so the GIL is dropped and recontended
    for during it rather than held straight through.

    Blocked fraction 0 means the GIL was there whenever this thread wanted
    it. Near 1 means a Python thread wanting to run mostly cannot, which
    is what "the pipeline is GIL-bound" would look like.

    The probe competes for the same GIL it measures, so it runs at a low
    duty cycle and its own cost is reported alongside its answer.
    """

    def __init__(self, cycle_rate, *, burst_seconds=0.02, period_seconds=0.25):
        super().__init__(name="orgui-gil-probe", daemon=True)
        self._cycle_rate = cycle_rate
        self._burst_cycles = int(burst_seconds * cycle_rate)
        self._period = period_seconds
        self._halt = threading.Event()
        self.wall = 0.0
        self.cycles = 0
        self.bursts = 0

    def halt(self):
        self._halt.set()

    def run(self):
        while not self._halt.is_set():
            start_wall = perf_counter()
            start_cycles = thread_cycles()
            value = 0
            while thread_cycles() - start_cycles < self._burst_cycles:
                for _iteration in range(2000):
                    value += 1
            self.wall += perf_counter() - start_wall
            self.cycles += thread_cycles() - start_cycles
            self.bursts += 1
            self._halt.wait(self._period)

    def summary(self):
        """:rtype: dict"""
        on_core = self.cycles / self._cycle_rate if self._cycle_rate else 0.0
        return {
            "bursts": self.bursts,
            "wall_seconds": self.wall,
            "on_core_seconds": on_core,
            "blocked_fraction": (
                (self.wall - on_core) / self.wall if self.wall else None
            ),
        }


def uncontended(scan, correct, frame_indices, cycle_rate):
    """Load and correct a few frames serially, before the pipeline starts.

    The reference the in-pipeline numbers are read against: one thread, no
    GIL competition, cold cache. A stage whose in-pipeline on-core time
    matches this is contending for the GIL or the disk but not for memory
    bandwidth; one that inflates is doing the same arithmetic more slowly
    because 24 threads are sharing the same cache and controller.

    :param frame_indices:
        Frames outside the measured window -- reading a frame here would
        otherwise leave it warm in the page cache for the run.
    :rtype: dict
    """
    correct_frame = correct.correct_frame
    totals = {"load": [0.0, 0], "correct": [0.0, 0]}
    for frame_index in frame_indices:
        wall = perf_counter()
        cycles = thread_cycles()
        payload = scan.get_raw_img(frame_index)
        totals["load"][0] += perf_counter() - wall
        totals["load"][1] += thread_cycles() - cycles
        raw = np.asarray(payload.img)
        wall = perf_counter()
        cycles = thread_cycles()
        correct_frame(payload, raw, frame_index)
        totals["correct"][0] += perf_counter() - wall
        totals["correct"][1] += thread_cycles() - cycles
    count = max(1, len(frame_indices))
    return {
        name: {
            "wall_ms_per_frame": wall / count * 1e3,
            "on_core_ms_per_frame": cycles / cycle_rate / count * 1e3,
        }
        for name, (wall, cycles) in totals.items()
    }


class _TimedKernel:
    """Forwarding proxy timing a kernel's accumulate calls."""

    def __init__(self, kernel):
        self._kernel = kernel

    def __getattr__(self, name):
        return getattr(self._kernel, name)

    def accumulate_group(self, *arguments):
        with timed("kernel"):
            return self._kernel.accumulate_group(*arguments)

    def accumulate_group_tile(self, *arguments):
        with timed("kernel"):
            return self._kernel.accumulate_group_tile(*arguments)

    def accumulate(self, *arguments):
        with timed("kernel"):
            return self._kernel.accumulate(*arguments)


def _install(scan, correct, router, *, stack_threshold_bytes=1 << 20):
    """Wrap every stage boundary in place; returns the patched pipeline.

    Patches module globals rather than editing the pipeline, so the code
    being measured is exactly the code that ships. ``np.stack`` is patched
    globally but attributed only above ``stack_threshold_bytes``, which
    separates ``_map_frame_group``'s per-tile buffers from the small
    stacks the rest of the process makes.
    """
    scan_read = scan.get_raw_img

    def get_raw_img(frame_index, *arguments, **keywords):
        with timed("load"):
            return scan_read(frame_index, *arguments, **keywords)

    scan.get_raw_img = get_raw_img

    correct_frame = correct.correct_frame

    def timed_correct_frame(payload, raw, frame_index):
        with timed("correct"):
            return correct_frame(payload, raw, frame_index)

    def timed_correct(payload, raw, frame_index, tile):
        with timed("correct_tile"):
            return correct(payload, raw, frame_index, tile)

    timed_correct.correct_frame = timed_correct_frame
    timed_correct.repair_plan = correct.repair_plan

    build_kernels = reconstruction_job._build_kernels

    def timed_build_kernels(*arguments, **keywords):
        return {
            name: _TimedKernel(kernel)
            for name, kernel in build_kernels(*arguments, **keywords).items()
        }

    reconstruction_job._build_kernels = timed_build_kernels

    reduce_batches = reconstruction._reduce_batches

    def timed_reduce_batches(batches):
        with timed("merge_tiles"):
            return reduce_batches(batches)

    reconstruction._reduce_batches = timed_reduce_batches

    map_frame_group = reconstruction_job._map_frame_group

    def timed_map_frame_group(*arguments, **keywords):
        with timed("map_frame_group"):
            return map_frame_group(*arguments, **keywords)

    reconstruction_job._map_frame_group = timed_map_frame_group

    numpy_stack = np.stack

    def timed_stack(arrays, *arguments, **keywords):
        # Sized from the inputs, so a small stack anywhere else in the
        # process is never timed rather than timed and then discounted.
        first = arrays[0] if isinstance(arrays, list | tuple) and arrays else None
        if getattr(first, "nbytes", 0) * len(arrays) < stack_threshold_bytes:
            return numpy_stack(arrays, *arguments, **keywords)
        with timed("tile_stack"):
            return numpy_stack(arrays, *arguments, **keywords)

    np.stack = timed_stack

    route = router.route
    routed = {"records": 0, "calls": 0}

    def timed_route(grid_name, frame_index, batch, **keywords):
        routed["records"] += int(batch["chunk_id"].size)
        routed["calls"] += 1
        with timed("route"):
            return route(grid_name, frame_index, batch, **keywords)

    router.route = timed_route
    return timed_correct, routed


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--start", type=int, default=1638)
    parser.add_argument("--count", type=int, default=234)
    parser.add_argument("--ranges-per-checkpoint", type=int, default=2)
    parser.add_argument("--group", type=int, default=None)
    parser.add_argument("--threads-per-image", type=int, default=None)
    parser.add_argument("--memory-gib", type=float, default=None)
    parser.add_argument("--tile-rows", type=int, default=None)
    parser.add_argument("--depth", type=int, default=None)
    parser.add_argument(
        "--baseline-frames",
        type=int,
        default=3,
        help=(
            "Frames to load and correct serially before the run, as the "
            "uncontended reference. Taken from just before the window."
        ),
    )
    parser.add_argument(
        "--correction",
        choices=("full", "minimal"),
        default="full",
        help=(
            "'minimal' keeps the mask and drops the scaling arithmetic, as "
            "an upper bound on what fusing correction could buy. It routes "
            "the same records and carries wrong intensities."
        ),
    )
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main():
    """Run the real mapping pipeline and report where its time went."""
    arguments = _arguments()
    cycle_rate = calibrate_cycle_rate()

    job = read_job(arguments.job)
    config = job.config_data
    scan = job.scan
    spec = job.internal_spec()
    if arguments.depth is not None:
        spec = replace(
            spec,
            max_depth=arguments.depth,
            work_block_pixels=resolve_work_block_pixels(
                job.work_block_pixels,
                arguments.depth,
                spec.memory_budget_bytes,
                spec.threads,
            ),
        )
    if arguments.group is not None:
        spec = replace(spec, frames_per_group=arguments.group)
    bounds = scan.exposure_angle_bounds(config, fallback=job.angle_fallback)
    ranges, tiles = _execution_layout(job, scan, config)
    if arguments.tile_rows:
        detector_rows, detector_columns = config.detector.detector.shape
        tiles = [
            (row, min(row + arguments.tile_rows, detector_rows), 0, detector_columns)
            for row in range(0, detector_rows, arguments.tile_rows)
        ]
    kept, boundaries = _window_ranges(
        ranges,
        arguments.start,
        arguments.start + arguments.count,
        arguments.ranges_per_checkpoint,
    )
    grid_names = [grid.grid_name for grid in spec.grids]
    plan = {grid_name: list(boundaries) for grid_name in grid_names}

    effective_memory = (
        int(arguments.memory_gib * 1024**3)
        if arguments.memory_gib
        else job.memory_override_bytes or job.runtime_memory_bytes
    )
    active_budget_bytes, effective_memory = split_memory_budget(
        effective_memory, max(1, len(spec.grids))
    )

    scratch = Path(
        arguments.scratch or tempfile.mkdtemp(prefix="orgui-rsmap-stages-")
    )
    checkpoint_dir = scratch / "checkpoints"
    if checkpoint_dir.exists():
        shutil.rmtree(checkpoint_dir)
    router = _CheckpointRouter(
        plan,
        spec_digest=job.digest,
        checkpoint_dir=checkpoint_dir,
        active_budget_bytes=active_budget_bytes,
    )
    assets = _load_assets(job)
    correct = _correction_pipeline(
        config, scan, assets, _base_provenance(job, config)
    )
    if arguments.correction == "minimal":
        correct = minimal_correction(
            correct,
            np.ascontiguousarray(assets["mask"], dtype=bool)
            if "mask" in assets
            else None,
        )
    baseline = (
        uncontended(
            scan,
            correct,
            [
                frame
                for frame in range(
                    kept[0][0] - arguments.baseline_frames, kept[0][0]
                )
                if frame >= 0
            ],
            cycle_rate,
        )
        if arguments.baseline_frames
        else None
    )
    timed_correct, routed = _install(scan, correct, router)
    total_images = sum(stop - start for start, stop in kept)

    grid = spec.grids[0]
    header = {
        "job": str(arguments.job),
        "frames": [kept[0][0], kept[-1][1]],
        "total_images": total_images,
        "scheduling_ranges": len(kept),
        "checkpoints": len(boundaries),
        "tiles": len(tiles),
        "grid_shape": list(grid.shape),
        "grid_step": [float(step) for step in grid.step],
        "frames_per_group": getattr(spec, "frames_per_group", None),
        "max_depth": spec.max_depth,
        "work_block_pixels": spec.work_block_pixels,
        "threads": spec.threads,
        "angle_fallback": job.angle_fallback,
        "effective_memory_bytes": int(effective_memory),
        "correction": arguments.correction,
        "cycle_rate_ghz": round(cycle_rate / 1e9, 4),
        "native": _native_module().__file__,
        "scratch": str(scratch),
    }
    print(json.dumps(header, sort_keys=True), flush=True)

    messages = []

    def capture(value, maximum, message):
        messages.append(message)

    process = None
    try:
        import psutil

        process = psutil.Process()
        cpu_before = process.cpu_times()
    except ImportError:
        cpu_before = None

    idle_probe = GilProbe(cycle_rate)
    idle_probe.start()
    idle_probe._halt.wait(1.0)
    idle_probe.halt()
    idle_probe.join(timeout=5.0)
    probe = GilProbe(cycle_rate)
    probe.start()

    started = perf_counter()
    _map_pending_ranges(
        spec,
        scan,
        config,
        bounds,
        tiles,
        kept,
        router,
        correction_pipeline=timed_correct,
        effective_memory=effective_memory,
        threads_per_image=(
            arguments.threads_per_image
            if arguments.threads_per_image is not None
            else job.threads_per_image
        ),
        accumulation_budget_bytes=job.accumulation_budget_bytes,
        total_images=total_images,
        completed_images=0,
        progress=capture,
    )
    elapsed = perf_counter() - started
    probe.halt()
    probe.join(timeout=5.0)
    process_cpu = None
    if process is not None:
        cpu_after = process.cpu_times()
        process_cpu = (
            cpu_after.user - cpu_before.user + cpu_after.system - cpu_before.system
        )

    stages = {}
    for name, (wall, cycles, calls) in merged_stages().items():
        on_core = cycles / cycle_rate if cycle_rate else float("nan")
        stages[name] = {
            "calls": calls,
            "wall_ms_per_frame": wall / total_images * 1e3,
            "on_core_ms_per_frame": on_core / total_images * 1e3,
            "blocked_ms_per_frame": (wall - on_core) / total_images * 1e3,
        }
    # The one stage with no call of its own: what _map_frame_group spends
    # outside the kernel, the merge and the route is the per-tile slice
    # and stack, plus the Python around them.
    parent = stages.get("map_frame_group")
    if parent is not None:
        # Correction is deliberately absent: on the grouped scheduler it
        # runs in the prepare pool, outside this call entirely.
        children = ("kernel", "merge_tiles", "route", "correct_tile")
        assembly = {
            key: parent[key]
            - sum(stages.get(child, {}).get(key, 0.0) for child in children)
            for key in (
                "wall_ms_per_frame",
                "on_core_ms_per_frame",
                "blocked_ms_per_frame",
            )
        }
        assembly["calls"] = parent["calls"]
        stages["tile_assembly (derived)"] = assembly

    # Every stage's on-core time is the *calling* thread's. The kernel
    # spawns its own native worker pool, whose cycles no Python thread
    # ever accounts for, so what the process spent minus what the Python
    # threads spent is (to within the instrument's reach) the native
    # pool's own CPU.
    python_cpu_ms_per_frame = sum(
        values["on_core_ms_per_frame"]
        for name, values in stages.items()
        if name in {"load", "correct", "correct_tile", "map_frame_group"}
    )
    result = {
        "seconds": elapsed,
        "seconds_per_frame": elapsed / total_images,
        "ms_per_frame": elapsed / total_images * 1e3,
        "routed_records_per_frame": routed["records"] / total_images,
        "process_cpu_seconds": process_cpu,
        "process_cpu_ms_per_frame": (
            process_cpu / total_images * 1e3 if process_cpu else None
        ),
        "python_thread_cpu_ms_per_frame": python_cpu_ms_per_frame,
        "cores_in_use": (process_cpu / elapsed) if process_cpu else None,
        "final_layout": messages[-1] if messages else None,
        "gil_probe": probe.summary(),
        "gil_probe_idle": idle_probe.summary(),
        "uncontended": baseline,
        "stages": stages,
        "grids": _summarize_checkpoints(checkpoint_dir, grid_names),
    }
    print(json.dumps(result, sort_keys=True), flush=True)

    print()
    print(
        f"{total_images} frames in {elapsed:.1f} s "
        f"({elapsed / total_images * 1e3:.1f} ms/frame)"
    )
    if process_cpu:
        print(
            f"  {process_cpu:.1f} s of CPU -> "
            f"{process_cpu / elapsed:.2f} cores in use "
            f"({process_cpu / total_images * 1e3:.0f} ms of CPU per frame, "
            f"{python_cpu_ms_per_frame:.0f} of it on Python threads)"
        )
    during = probe.summary()
    idle = idle_probe.summary()
    if during["blocked_fraction"] is not None:
        print(
            f"  GIL probe: a Python thread waited "
            f"{during['blocked_fraction'] * 100:.0f}% of its wall time during "
            f"the run ({idle['blocked_fraction'] * 100:.0f}% idle), over "
            f"{during['bursts']} bursts"
        )
    if baseline:
        print(
            "  uncontended, one thread: "
            + ", ".join(
                f"{name} {values['wall_ms_per_frame']:.1f} ms wall / "
                f"{values['on_core_ms_per_frame']:.1f} on-core"
                for name, values in baseline.items()
            )
        )
    if messages:
        print(f"  {messages[-1]}")
    print()
    print(
        f"  {'stage':<26}{'calls':>8}{'wall':>10}{'on-core':>10}"
        f"{'blocked':>10}   (ms per frame)"
    )
    order = [
        "load",
        "correct",
        "correct_tile",
        "map_frame_group",
        "tile_assembly (derived)",
        "tile_stack",
        "kernel",
        "merge_tiles",
        "route",
    ]
    for name in order:
        values = stages.get(name)
        if values is None:
            continue
        indent = "    " if name in {
            "tile_assembly (derived)",
            "tile_stack",
            "kernel",
            "merge_tiles",
            "route",
        } else "  "
        print(
            f"  {indent}{name:<{26 - len(indent) + 2}}{values['calls']:>8}"
            f"{values['wall_ms_per_frame']:>10.1f}"
            f"{values['on_core_ms_per_frame']:>10.1f}"
            f"{values['blocked_ms_per_frame']:>10.1f}"
        )
    for grid_name, values in result["grids"].items():
        print(
            f"  {grid_name}: {values['rows']:,} rows, "
            f"{values['contributors']:,} contributors, "
            f"fingerprint {values['voxel_fingerprint']}"
        )

    if arguments.output:
        arguments.output.write_text(
            json.dumps({"header": header, "result": result}, indent=2),
            encoding="utf-8",
        )
    if arguments.scratch is None:
        shutil.rmtree(scratch, ignore_errors=True)


if __name__ == "__main__":
    main()

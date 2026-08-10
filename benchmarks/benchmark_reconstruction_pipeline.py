"""End-to-end reciprocal-space mapping throughput on a prepared job.

``benchmark_reconstruction_group.py`` times one kernel call. This times the
whole mapping phase around it: the real ``_map_pending_ranges`` pipeline --
reader pool, whole-frame correction, kernel, checkpoint router, checkpoint
writes -- over a bounded window of frames. That is the thing frame grouping
has to actually improve, and it is the thing a faster kernel can fail to
improve, because the GIL-held Python work per frame does not get cheaper
when the native call does.

Checkpoint boundaries are synthesised over the measured window rather than
taken from the job's own plan: a partial run against the real plan would
leave every checkpoint short of its frame count, so nothing would flush and
there would be no output to compare. Each synthesised checkpoint is a whole
number of the job's own scheduling ranges, so it stays contiguous and
exclusion-free exactly as a planned one is.

The written checkpoint parts are read back and summarised. Frame grouping is
deliberately *not* bit-for-bit -- contributions from several frames merge in
the block map instead of in the tree accumulator -- so the summary reports
what must hold instead: the same voxel keys, the same contributor counts,
and totals equal to rounding.

**Run this with ``PYTHONPATH`` pointing at the checkout.** ``python
benchmarks/<script>.py`` puts ``benchmarks/`` on ``sys.path``, not the
repository root, so a pip-installed ``orgui`` shadows the working tree and
the run silently measures a binary you did not just build. The resolved
package path is reported for that reason.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import shutil
import tempfile
from time import perf_counter

import h5py
import numpy as np

from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointRouter,
    _native_module,
    _PARTIAL_COLUMNS,
)
import orgui.reconstruction_job as reconstruction_job
from orgui.reconstruction_job import (
    _base_provenance,
    _correction_pipeline,
    _execution_layout,
    _frame_parallelism,
    _group_pipeline_layout,
    _load_assets,
    _map_pending_ranges,
    read_job,
    split_memory_budget,
)


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--start", type=int, default=1800)
    parser.add_argument("--count", type=int, default=160)
    parser.add_argument(
        "--ranges-per-checkpoint",
        type=int,
        default=2,
        help="Scheduling ranges per synthesised checkpoint boundary.",
    )
    parser.add_argument(
        "--group",
        type=int,
        default=1,
        help="Frames per kernel call (spec.frames_per_group).",
    )
    parser.add_argument(
        "--threads-per-image",
        type=int,
        default=None,
        help="Fixed kernel_threads. Omitted means the job's own setting.",
    )
    parser.add_argument(
        "--memory-gib",
        type=float,
        default=None,
        help="Override the job's own memory budget, for both arms alike.",
    )
    parser.add_argument(
        "--tile-rows",
        type=int,
        default=None,
        help=(
            "Override the detector band height. Smaller bands cut a "
            "worker's native working set, so they buy image_workers back "
            "at the same total budget rather than by raising it."
        ),
    )
    parser.add_argument(
        "--compute-workers",
        type=int,
        default=None,
        help=(
            "Override the grouped scheduler's concurrent group calls. "
            "Native threads are split evenly between them."
        ),
    )
    parser.add_argument(
        "--depth",
        type=int,
        default=None,
        help="Override spec.max_depth (the job's accuracy setting).",
    )
    parser.add_argument("--scratch", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def _window_ranges(ranges, start, stop, ranges_per_checkpoint):
    """Scheduling ranges inside ``[start, stop)``, and the checkpoint
    boundaries covering exactly those frames.

    A boundary may only merge scheduling ranges that are adjacent: a gap
    between them is an excluded frame, and a checkpoint spanning one would
    wait forever for a frame that is never routed.
    """
    kept = [
        frame_range
        for frame_range in ranges
        if start <= frame_range[0] and frame_range[1] <= stop
    ]
    if not kept:
        raise SystemExit(
            f"No whole scheduling range lies inside frames [{start}, {stop})"
        )
    boundaries = []
    run = [kept[0]]
    for frame_range in kept[1:]:
        adjacent = frame_range[0] == run[-1][1]
        if adjacent and len(run) < max(1, ranges_per_checkpoint):
            run.append(frame_range)
            continue
        boundaries.append((run[0][0], run[-1][1]))
        run = [frame_range]
    boundaries.append((run[0][0], run[-1][1]))
    return kept, boundaries


_SUMMARY_CHUNK_ROWS = 1 << 22


def _summarize_checkpoints(checkpoint_dir, grid_names):
    """Contributor counts, totals and a key fingerprint over every part.

    These are the quantities that must agree across group sizes, in place
    of a digest: grouping changes floating-point association, so sums move
    in the last bits, while which voxels were reached and how many samples
    reached each of them must not move at all.

    ``voxel_fingerprint`` is ``sum over rows of h(key) * contributors``,
    accumulated in ``uint64`` and allowed to wrap. Because it is a sum, it
    does not care what order rows arrive in or how they were split across
    parts -- it equals the same sum taken over distinct voxels with their
    total contributor counts, which is exactly the invariant being
    claimed. Two runs reaching different voxels, or apportioning samples
    differently between them, disagree here even when every total below
    matches.

    Streamed in row chunks. A whole window's parts concatenated in memory
    is tens of gigabytes on a real job -- and would be paid on top of the
    mapping pipeline's own peak, which is already sized to the job's full
    memory budget.
    """
    summary = {}
    for grid_name in grid_names:
        parts = sorted((Path(checkpoint_dir) / grid_name).glob("ckpt*.h5"))
        frames_covered = 0
        rows = 0
        contributors = 0
        totals = dict.fromkeys(
            ("weight", "weighted_intensity", "weighted_variance"), 0.0
        )
        fingerprint = np.uint64(0)
        for part in parts:
            with h5py.File(part, "r") as handle:
                frames_covered += int(handle.attrs["frames_covered"])
                part_rows = int(handle["chunk_id"].shape[0])
                rows += part_rows
                for begin in range(0, part_rows, _SUMMARY_CHUNK_ROWS):
                    end = min(begin + _SUMMARY_CHUNK_ROWS, part_rows)
                    block = {
                        name: handle[name][begin:end]
                        for name in _PARTIAL_COLUMNS
                    }
                    contributors += int(block["contributors"].sum())
                    for name in totals:
                        totals[name] += float(block[name].sum())
                    keys = (
                        block["chunk_id"].astype(np.uint64)
                        * np.uint64(1000003)
                        + block["local_voxel_id"].astype(np.uint64)
                    )
                    with np.errstate(over="ignore"):
                        fingerprint = np.uint64(
                            fingerprint
                            + (
                                keys * block["contributors"].astype(np.uint64)
                            ).sum()
                        )
        summary[grid_name] = {
            "parts": len(parts),
            "frames_covered": frames_covered,
            "rows": rows,
            "contributors": contributors,
            "voxel_fingerprint": int(fingerprint),
            **totals,
        }
    return summary


def main():
    """Time the real mapping pipeline over a bounded window of frames."""
    arguments = _arguments()
    job = read_job(arguments.job)
    config = job.config_data
    scan = job.scan
    spec = job.internal_spec()
    if arguments.depth is not None:
        spec = replace(spec, max_depth=arguments.depth)
    if arguments.group != 1:
        if not hasattr(spec, "frames_per_group"):
            raise SystemExit(
                "This checkout has no spec.frames_per_group; --group only "
                "works once frame grouping is wired in"
            )
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
        arguments.scratch
        or tempfile.mkdtemp(prefix="orgui-rsmap-bench-")
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
    # Records *entering* the checkpoint layer, which is what frame
    # grouping is meant to reduce and what every stage downstream of the
    # kernel scales with. Rows on disk are the wrong measure: the
    # checkpoint accumulator already merges across frames, so it collapses
    # the same redundancy later and more expensively, and the saving never
    # shows up there.
    routed = {"records": 0, "calls": 0}
    original_route = router.route

    def counting_route(grid_name, frame_index, batch, **kwargs):
        routed["records"] += int(batch["chunk_id"].size)
        routed["calls"] += 1
        return original_route(grid_name, frame_index, batch, **kwargs)

    router.route = counting_route

    correct = _correction_pipeline(
        config, scan, _load_assets(job), _base_provenance(job, config)
    )
    total_images = sum(stop - start for start, stop in kept)
    seed_workers, seed_threads, worker_bytes, _accumulation = (
        _frame_parallelism(
            spec,
            tiles,
            effective_memory,
            threads_per_image=1,
            accumulation_budget_bytes=job.accumulation_budget_bytes,
        )
    )

    header = {
        "job": str(arguments.job),
        "frames": [kept[0][0], kept[-1][1]],
        "total_images": total_images,
        "scheduling_ranges": len(kept),
        "checkpoints": len(boundaries),
        "tiles": len(tiles),
        "largest_tile_pixels": max(
            (row_stop - row_start) * (column_stop - column_start)
            for row_start, row_stop, column_start, column_stop in tiles
        ),
        "frames_per_group": getattr(spec, "frames_per_group", 1),
        "max_depth": spec.max_depth,
        "work_block_pixels": spec.work_block_pixels,
        "threads": spec.threads,
        "angle_fallback": job.angle_fallback,
        "effective_memory_bytes": int(effective_memory),
        "active_budget_bytes": int(active_budget_bytes),
        # What the memory budget affords. Above one frame per group the
        # scheduler uses _group_pipeline_layout instead, which the final
        # progress line reports; these are the per-frame pipeline's.
        "seed_image_workers": int(seed_workers),
        "seed_kernel_threads": int(seed_threads),
        "worker_memory_mb": round(worker_bytes / 1e6, 1),
        "group_layout": (
            list(_group_pipeline_layout(spec, tiles, effective_memory, 4))
            if getattr(spec, "frames_per_group", 1) > 1
            else None
        ),
        "native": _native_module().__file__,
        "scratch": str(scratch),
    }
    print(json.dumps(header, sort_keys=True), flush=True)

    if arguments.compute_workers:
        # A benchmark-only hook. The layout is a memory decision inside
        # the scheduler, and the point of sweeping it is to find out what
        # that decision should be aiming at.
        original_layout = reconstruction_job._group_pipeline_layout

        def fixed_layout(spec, tiles, memory_bytes, prepare_workers):
            _workers, _threads, depth = original_layout(
                spec, tiles, memory_bytes, prepare_workers
            )
            workers = max(1, arguments.compute_workers)
            return (
                workers,
                max(1, spec.threads // workers),
                workers + max(1, depth - _workers),
            )

        reconstruction_job._group_pipeline_layout = fixed_layout

    messages = []

    def capture(value, maximum, message):
        messages.append(message)

    started = perf_counter()
    _map_pending_ranges(
        spec,
        scan,
        config,
        bounds,
        tiles,
        kept,
        router,
        correction_pipeline=correct,
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

    result = {
        "seconds": elapsed,
        "frames_per_second": total_images / elapsed,
        "seconds_per_frame": elapsed / total_images,
        "routed_records": routed["records"],
        "routed_records_per_frame": routed["records"] / total_images,
        "route_calls": routed["calls"],
        "final_layout": messages[-1] if messages else None,
        "checkpoints_written": len(router.written),
        "grids": _summarize_checkpoints(checkpoint_dir, grid_names),
    }
    print(json.dumps(result, sort_keys=True), flush=True)

    print()
    print(f"{total_images} frames in {elapsed:.1f} s "
          f"({total_images / elapsed:.3f} frames/s, "
          f"{elapsed / total_images * 1000:.1f} ms/frame)")
    if messages:
        print(f"  {messages[-1]}")
    print(f"  routed {routed['records']:,} records in {routed['calls']} calls "
          f"({routed['records'] / total_images:,.0f} per frame); "
          f"per-frame layout {seed_workers} workers x {seed_threads} threads")
    for grid_name, values in result["grids"].items():
        print(
            f"  {grid_name}: {values['rows']:,} rows, "
            f"{values['contributors']:,} contributors, "
            f"weight {values['weight']:.6e}"
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

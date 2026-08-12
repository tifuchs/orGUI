"""End-to-end tests for the centralized reconstruction job."""

import dataclasses
from hashlib import sha256
import json
import threading
import time

import h5py
import numpy as np
import pytest

from orgui.app.config_data import ConfigData
from orgui.app.database import config_data_to_json
from orgui.backend.scans import ScanReference, SimulationScan
from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg
from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointRouter,
    _GridSpec,
    _ReconstructionSpec,
)
import orgui.reconstruction_job as reconstruction_job_module
from orgui.reconstruction_job import (
    MAX_CONCURRENT_ACTIVE_CHECKPOINTS,
    WORK_BLOCK_PRESETS,
    _FRAME_RECORD_BYTES_PER_PIXEL,
    _PYTHON_CORRECTION_BYTES_PER_PIXEL,
    _GROUP_SIZE_CANDIDATES,
    _RESERVED_RECORDS_PER_PIXEL,
    ReconstructionGrid,
    ReconstructionJob,
    _angles_advance_monotonically,
    _frame_groups,
    _choose_frames_per_group,
    _frame_parallelism,
    _group_pipeline_layout,
    _node_checkpoint_plan,
    _node_excluded_frames,
    job_status,
    reconstruction_execution_settings,
    resolve_work_block_pixels,
    split_memory_budget,
    work_block_memory_cap,
    run_cluster_finalize,
    run_cluster_map_task,
    run_job,
    write_job,
)


pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


def _config():
    cell = CTRcalc.UnitCell([3.0, 3.0, 3.0], [90.0, 90.0, 90.0])
    cell.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    ub = HKLVlieg.UBCalculator(cell, 70.0)
    ub.defaultU()
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(
        729.0, 0.5, 0.5, pixelX=172.0, pixelY=172.0
    )
    detector.set_wavelength(ub.getLambda() * 1e-10)
    detector.detector.shape = (2, 2)
    detector.detector.max_shape = (2, 2)
    return ConfigData(detector, cell, ub)


def _two_frame_job(tmp_path, output_name):
    scan = SimulationScan((2, 2), 0.0, 1.0, 2)
    assets = tmp_path / "scratch" / "job-assets.nxs"
    assets.parent.mkdir()
    with h5py.File(assets, "w") as h5file:
        h5file.attrs["orgui_job_assets"] = 1
    grid = ReconstructionGrid(
        minimum=(-20.0, -20.0, -20.0),
        maximum=(20.0, 20.0, 20.0),
        step=(20.0, 20.0, 20.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    scan_reference = ScanReference.from_scan(scan).to_dict()
    job = ReconstructionJob(
        config=config_data_to_json(_config()),
        scan_reference=scan_reference,
        grids=[grid.__dict__],
        scratch_path=str(assets.parent),
        output_path=str(tmp_path / output_name),
        compression="Raw",
        assets_path=str(assets),
        assets_sha256=sha256(assets.read_bytes()).hexdigest(),
        source_fingerprint_sha256=sha256(
            json.dumps(
                scan_reference, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest(),
        threads_per_image=1,
        accumulation_budget_bytes=None,
        checkpoint_count=2,
        # One checkpoint per frame -- the grid name mirrors _GridSpec's
        # default naming for an unnamed "lab" frame grid ("q_lab").
        checkpoint_plan={"q_lab": [[0, 1], [1, 2]]},
        runtime_threads=2,
        runtime_memory_bytes=64 * 1024 * 1024,
        tile_shape=(1, 1),
        frame_batch=1,
    )
    return scan, job


def test_central_job_runs_resumes_and_cleans_verified_scratch(
    tmp_path, monkeypatch
):
    scan, job = _two_frame_job(tmp_path, "result.h5")
    job_path = tmp_path / "job.json"
    write_job(job, job_path)

    active_loads = 0
    maximum_active_loads = 0
    load_lock = threading.Lock()
    original_get_raw_img = SimulationScan.get_raw_img

    def observed_get_raw_img(self, index):
        nonlocal active_loads, maximum_active_loads
        with load_lock:
            active_loads += 1
            maximum_active_loads = max(maximum_active_loads, active_loads)
        try:
            time.sleep(0.02)
            return original_get_raw_img(self, index)
        finally:
            with load_lock:
                active_loads -= 1

    monkeypatch.setattr(
        SimulationScan, "get_raw_img", observed_get_raw_img
    )

    settings = reconstruction_execution_settings(
        job, scan=scan, config=job.config_data
    )
    assert settings["thread_budget"] == 2
    assert settings["native_threads_per_image"] == 1
    assert settings["detector_tile_shape"] == (1, 1)
    assert settings["map_tasks"] == 2
    assert settings["checkpoint_count"] == 2
    layout = settings["parallel_layouts"][0]
    assert layout["exposure"] == "stationary"
    assert layout["concurrent_image_workers"] == 2
    assert layout["native_threads_per_image"] == 1
    assert layout["tiles_per_image"] == 4
    assert layout["memory_per_image_MiB"] == pytest.approx(32.0, abs=0.01)
    assert layout["accumulation_MiB_per_worker"] == pytest.approx(
        10.33, abs=0.02
    )

    job.threads_per_image = 2
    combined = reconstruction_execution_settings(
        job, scan=scan, config=job.config_data
    )
    assert combined["native_threads_per_image"] == 2
    assert combined["parallel_layouts"][0]["concurrent_image_workers"] == 1
    assert combined["parallel_layouts"][0]["native_threads_per_image"] == 2
    job.threads_per_image = 1

    progress_updates = []
    first = run_job(
        job_path,
        progress=lambda value, maximum, message: progress_updates.append(
            (value, maximum, message)
        ),
    )
    resumed = run_job(job_path)

    assert maximum_active_loads == 2
    assert any(
        "Mapping images 1/2" in message
        for _value, _maximum, message in progress_updates
    )
    assert any(
        "Mapping images 2/2" in message
        for _value, _maximum, message in progress_updates
    )
    assert first == resumed
    assert first["status"] == "complete"
    assert not (tmp_path / "scratch" / "job-assets.nxs").exists()
    with h5py.File(first["output_path"], "r") as h5file:
        assert "entry/configuration" in h5file
        group = h5file["entry/reconstruction/results/q_lab"]
        assert np.nanmax(group["intensity"]) == 10.0
        assert np.nanmax(group["contributors"]) == 8
    assert job_status(job_path)["output_sha256"] == first["output_sha256"]
    assert job_status(job_path)["checkpoints"] == {
        "completed": 2,
        "pending": 0,
        "total": 2,
    }


def _uniform_tiles(total_side, tile_side):
    return [
        (row, row + tile_side, column, column + tile_side)
        for row in range(0, total_side, tile_side)
        for column in range(0, total_side, tile_side)
    ]


def test_memory_budget_is_split_between_accumulators_and_pipeline():
    """The checkpoint accumulators and the frame pipeline used to be sized
    against the whole budget independently, so their sum could reach twice
    what the user asked for. One split has to divide it between them.
    """
    budget = 10 * 1024**3

    per_checkpoint, pipeline = split_memory_budget(budget, 1)

    assert per_checkpoint * MAX_CONCURRENT_ACTIVE_CHECKPOINTS + pipeline <= budget
    # Every grid gets its own concurrent checkpoints, so more grids must
    # not multiply the accumulators' total share.
    two_grids, pipeline_two = split_memory_budget(budget, 2)
    assert two_grids * MAX_CONCURRENT_ACTIVE_CHECKPOINTS * 2 + pipeline_two <= budget


def test_frame_parallelism_accounts_for_a_frames_own_records():
    """A worker holds its frame's mapped records -- every tile's output
    batch and the transient copy merging them makes -- not just the image
    and correction buffers. Sizing the pool without that term let measured
    peaks run about a third above what the pool believed it had claimed.
    """
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0, threads=8)
    tiles = [(0, 1024, 0, 1024)]

    _workers, _threads, per_worker, _accumulation = _frame_parallelism(
        spec, tiles, 8 * 1024**3, threads_per_image=1
    )

    detector_pixels = 1024 * 1024
    assert per_worker >= detector_pixels * (
        _PYTHON_CORRECTION_BYTES_PER_PIXEL + _FRAME_RECORD_BYTES_PER_PIXEL
    )


def test_work_block_preset_halves_with_adaptive_depth():
    """A preset names a target working set, not a raw pixel count: deeper
    subdivision adds per-pixel state competing for the same cache, so the
    pixel count halves per depth to hold that working set roughly fixed.
    """
    memory, threads = 64 * 1024**3, 8

    assert (
        resolve_work_block_pixels("medium", 0, memory, threads)
        == WORK_BLOCK_PRESETS["medium"]
    )
    assert (
        resolve_work_block_pixels("medium", 2, memory, threads)
        == WORK_BLOCK_PRESETS["medium"] // 4
    )
    # No setting at all behaves as the default preset does.
    assert resolve_work_block_pixels(None, 2, memory, threads) == (
        resolve_work_block_pixels("medium", 2, memory, threads)
    )


def test_work_block_explicit_count_is_taken_literally():
    """A pinned number means that number -- it is not a scale to rescale."""
    memory, threads = 64 * 1024**3, 8

    assert resolve_work_block_pixels(4096, 0, memory, threads) == 4096
    assert resolve_work_block_pixels(4096, 3, memory, threads) == 4096


def test_work_block_is_capped_by_the_memory_budget():
    """The kernel reserves an arena per worker thread proportional to the
    block, so no setting -- including a pinned one -- may push that past
    the job's budget."""
    threads = 16
    tight = 64 * 1024**2

    cap = work_block_memory_cap(0, tight, threads)
    assert resolve_work_block_pixels(1_000_000_000, 0, tight, threads) == cap
    assert resolve_work_block_pixels("maximum", 0, tight, threads) <= cap
    # A generous budget must not cap a sensible request.
    assert resolve_work_block_pixels("medium", 0, 64 * 1024**3, threads) == (
        WORK_BLOCK_PRESETS["medium"]
    )


def test_work_block_rejects_an_unknown_preset():
    with pytest.raises(ValueError, match="Unknown native work block preset"):
        resolve_work_block_pixels("enormous", 0, 64 * 1024**3, 8)


def test_frame_parallelism_scopes_native_memory_to_largest_tile_not_detector_sum():
    """A worker processes its detector tiles sequentially (_map_one_frame),
    one native accumulate() call at a time -- the native-kernel share of
    its peak memory is bounded by its single largest tile, never the sum
    of every tile across the whole detector. Splitting the same total
    detector area into more, smaller tiles must not increase (and here,
    since each tile shrinks, must decrease) the native contribution --
    summing silently starved image_workers on any multi-tile job (a real
    Pilatus-6M-scale job saw concurrent_image_workers capped at 6 instead
    of the full 24-thread budget purely from this overcounting)."""
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    # depth > 0 so the native (leaf-count-driven) share is large enough
    # to dominate the fixed per-pixel Python-side buffer cost, making the
    # native-only effect being tested clearly visible.
    spec = _ReconstructionSpec(grids=(grid,), max_depth=3, threads=24)
    memory_bytes = 10_000 * 1024**2
    one_big_tile = _uniform_tiles(2048, 2048)
    four_small_tiles = _uniform_tiles(2048, 1024)
    assert len(one_big_tile) == 1
    assert len(four_small_tiles) == 4
    # Same total detector area either way -- same Python-side correction
    # buffer cost -- only the native per-call tile size differs.
    assert sum(
        (r1 - r0) * (c1 - c0) for r0, r1, c0, c1 in one_big_tile
    ) == sum((r1 - r0) * (c1 - c0) for r0, r1, c0, c1 in four_small_tiles)

    single_result = _frame_parallelism(
        spec, one_big_tile, memory_bytes, threads_per_image=1
    )
    split_result = _frame_parallelism(
        spec, four_small_tiles, memory_bytes, threads_per_image=1
    )

    # The split (4x smaller largest tile) version's per-worker memory
    # must be strictly smaller -- proving the estimate tracks the single
    # largest tile, not a detector-wide total that would be identical
    # (and therefore give an identical result) in both cases.
    assert split_result[2] < single_result[2]
    assert split_result[0] >= single_result[0]


def test_frame_parallelism_native_estimate_does_not_grow_with_depth():
    """The per-pixel native estimate is bounded by records, not leaves.

    Sized by the worst-case leaf count it grows as ``4**depth`` (or
    ``8**depth`` when the exposure rotates), which claimed 5248 bytes at
    depth 3 for a pixel that really costs about 106 and so forced any
    realistic detector onto tiny tiles. Because one worker's estimate
    divides into the budget, that also capped how many frames could be in
    flight. Leaves are transient; records are what stays resident.
    """
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    tiles = _uniform_tiles(2048, 2048)
    memory_bytes = 10_000 * 1024**2
    results = {
        depth: _frame_parallelism(
            _ReconstructionSpec(grids=(grid,), max_depth=depth, threads=24),
            tiles,
            memory_bytes,
            threads_per_image=1,
        )
        for depth in (0, 1, 2, 3, 5, 8)
    }
    # Depth 0 is exact -- one leaf can only reach one voxel -- so it stays
    # below the saturated bound; every deeper setting shares one value.
    assert results[0][2] < results[1][2]
    deep = {depth: value[2] for depth, value in results.items() if depth}
    assert len(set(deep.values())) == 1, deep
    # And the concurrency that estimate funds must not collapse with depth.
    assert results[8][0] == results[1][0]
    assert results[8][0] > 1


def test_frame_parallelism_native_estimate_mirrors_the_kernel_precheck():
    """The Python estimate and the kernel's own guard must agree.

    They are two statements of one bound, and the kernel throws where this
    side merely sizes a pool, so a tile this side considers affordable must
    not be one the kernel rejects.
    """
    native = pytest.importorskip(
        "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
    )
    rows = columns = 512
    pixels = rows * columns
    rays = np.zeros((rows + 1, columns + 1, 3), dtype=np.float64)
    rays[..., 1] = 1.0
    arguments = (
        np.ones((rows, columns)),
        np.ones((rows, columns)),
        np.zeros((rows, columns), dtype=bool),
        rays,
        np.zeros(4),
        np.zeros(4),
    )

    def kernel_for(max_depth, memory_budget_bytes):
        return native.ReconstructionKernel(
            np.array([-1.0, -1.0, -1.0]),
            np.array([0.01, 0.01, 0.01]),
            np.array([256, 256, 256], dtype=np.int64),
            np.array([16, 16, 16], dtype=np.int64),
            "lab",
            1.0,
            np.eye(3),
            np.eye(3),
            max_depth,
            1,
            64,
            memory_budget_bytes,
        )

    for max_depth in (0, 1, 3, 5, 8):
        expected = (
            128 + 2 * _RESERVED_RECORDS_PER_PIXEL * 40
            if max_depth
            else 128 + 2 * 40
        )
        # Exactly what the Python estimate asks for must be accepted. This
        # is the direction that pins the two together: a leaf-sized bound
        # demands 4**depth times more and rejects this outright from depth
        # 2 up.
        kernel_for(max_depth, pixels * expected).accumulate(*arguments)
        # One byte per pixel short of it must not be.
        with pytest.raises(ValueError, match="native memory budget"):
            kernel_for(max_depth, pixels * (expected - 1)).accumulate(*arguments)


def test_frame_parallelism_accounts_for_python_correction_buffers_and_prefetch():
    """The native-kernel-only estimate misses two real, non-native memory
    costs (design doc Sec7): _correction_pipeline.correct_frame's
    full-detector-sized Python buffers (held once per frame, not per
    tile) and the prefetch pipeline's own read-ahead queue. A large,
    single-tile, depth-0 job (negligible native share) must still be
    bounded well below "infinite" image_workers by these alone."""
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0, threads=1000)
    # A Pilatus-6M-scale detector split into many small tiles: at depth=0
    # the native share (bounded by one small tile) is negligible, so only
    # the Python-side/prefetch accounting can be limiting concurrency
    # here.
    detector_side = 2500
    tiles = _uniform_tiles(detector_side, 250)
    memory_bytes = 10_000 * 1024**2

    image_workers, _kernel_threads, per_worker_memory, _accumulation = (
        _frame_parallelism(
            spec, tiles, memory_bytes, threads_per_image=1
        )
    )

    # With a 1000-thread budget and negligible native memory, an
    # accounting that only sees the native kernel would let this run
    # essentially unbounded (unrealistic -- it would try to hold far more
    # than 10 GiB of raw/corrected frames in flight at once). The
    # Python-side buffer + prefetch-queue reservation must keep it well
    # below the thread budget instead.
    assert image_workers < 1000
    assert per_worker_memory > detector_side * detector_side * 8


def test_threads_per_image_none_round_trips_and_reports_automatic_mode(
    tmp_path,
):
    """Sec7 Phase 4b: threads_per_image=None ('automatic') must survive a
    JSON round-trip under the bumped schema version and be reported
    distinctly from a pinned int, without requiring a live run."""
    scan, job = _two_frame_job(tmp_path, "result.h5")
    job.threads_per_image = None

    assert job.schema_version == reconstruction_job_module.JOB_SCHEMA_VERSION
    restored = ReconstructionJob.from_dict(job.to_dict())
    assert restored.threads_per_image is None

    settings = reconstruction_execution_settings(
        job, scan=scan, config=job.config_data
    )
    assert settings["threads_per_image_mode"] == "automatic"
    # The automatic-mode seed matches _map_pending_ranges' own I/O-optimistic
    # starting point (kernel_threads=1) -- not a claim about the pair that
    # will actually be chosen once a run has real timing data.
    assert settings["native_threads_per_image"] == 1

    job.threads_per_image = 2
    pinned_settings = reconstruction_execution_settings(
        job, scan=scan, config=job.config_data
    )
    assert pinned_settings["threads_per_image_mode"] == "pinned"
    assert pinned_settings["native_threads_per_image"] == 2


def test_job_resumes_partial_checkpoints_without_remapping_completed_ones(
    tmp_path, monkeypatch
):
    """Interrupting after the first checkpoint, then resuming, must not
    reload frames whose checkpoint already fully covers them."""
    scan, job = _two_frame_job(tmp_path, "resumed-result.h5")
    job_path = tmp_path / "job.json"
    write_job(job, job_path)

    loaded_frames = []
    decision_lock = threading.Lock()
    first_claimed = threading.Event()
    original_get_raw_img = SimulationScan.get_raw_img

    def failing_after_first(self, index):
        with decision_lock:
            loaded_frames.append(index)
            allowed = not first_claimed.is_set()
            if allowed:
                first_claimed.set()
        if not allowed:
            raise RuntimeError("simulated crash after first checkpoint")
        return original_get_raw_img(self, index)

    monkeypatch.setattr(SimulationScan, "get_raw_img", failing_after_first)
    with pytest.raises(RuntimeError, match="simulated crash"):
        run_job(job_path)

    assert job_status(job_path)["status"] != "complete"
    # Both frames were attempted (one races the other for "first"), but
    # exactly one completed and got a durable checkpoint file before the
    # other's simulated crash aborted the run.
    assert len(loaded_frames) == 2
    loaded_frames.clear()

    def counting_get_raw_img(self, index):
        loaded_frames.append(index)
        return original_get_raw_img(self, index)

    monkeypatch.setattr(
        SimulationScan, "get_raw_img", counting_get_raw_img
    )

    result = run_job(job_path)

    assert result["status"] == "complete"
    # Only the one frame whose checkpoint was not yet resumable was
    # reloaded -- not both frames from scratch.
    assert len(loaded_frames) == 1


def _multi_frame_job(
    tmp_path, output_name, frame_count, *, cluster=False, array_task_count=None
):
    scan = SimulationScan((2, 2), 0.0, 1.0, frame_count)
    assets = tmp_path / "scratch" / "job-assets.nxs"
    assets.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(assets, "w") as h5file:
        h5file.attrs["orgui_job_assets"] = 1
    grid = ReconstructionGrid(
        minimum=(-20.0, -20.0, -20.0),
        maximum=(20.0, 20.0, 20.0),
        step=(20.0, 20.0, 20.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    scan_reference = ScanReference.from_scan(scan).to_dict()
    job = ReconstructionJob(
        config=config_data_to_json(_config()),
        scan_reference=scan_reference,
        grids=[grid.__dict__],
        scratch_path=str(assets.parent),
        output_path=str(tmp_path / output_name),
        compression="Raw",
        assets_path=str(assets),
        assets_sha256=sha256(assets.read_bytes()).hexdigest(),
        source_fingerprint_sha256=sha256(
            json.dumps(
                scan_reference, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest(),
        threads_per_image=1,
        accumulation_budget_bytes=None,
        checkpoint_count=frame_count,
        # Cluster jobs never populate checkpoint_plan (design doc Sec13:
        # per-node plans live in scratch-local sidecars, not the job
        # JSON); single-node jobs get one checkpoint per frame, mirroring
        # _two_frame_job's convention.
        checkpoint_plan=(
            {}
            if cluster
            else {"q_lab": [[i, i + 1] for i in range(frame_count)]}
        ),
        runtime_threads=2,
        runtime_memory_bytes=64 * 1024 * 1024,
        tile_shape=(1, 1),
        frame_batch=1,
        cluster_settings=(
            {"array_task_count": array_task_count} if array_task_count else {}
        ),
    )
    return scan, job


@pytest.mark.parametrize(
    "scan_length,excluded,total_tasks",
    [
        (10, set(), 3),
        (10, {3, 4, 5}, 3),
        (7, {0, 6}, 4),
        (1, set(), 5),
        (20, {2, 5, 17}, 6),
    ],
)
def test_node_excluded_frames_partitions_included_frames_exactly(
    scan_length, excluded, total_tasks
):
    """Every included frame must be owned by exactly one task_index; no
    frame is owned by zero or more than one (design doc Sec13: disjoint,
    position-based node slicing)."""
    included = {index for index in range(scan_length) if index not in excluded}
    owners = {}
    for task_index in range(total_tasks):
        node_excluded = _node_excluded_frames(
            scan_length, excluded, total_tasks, task_index
        )
        owned = set(range(scan_length)) - node_excluded
        assert owned <= included
        for frame in owned:
            assert frame not in owners, f"frame {frame} owned by multiple nodes"
            owners[frame] = task_index
    assert set(owners) == included


def test_node_checkpoint_plan_computes_once_and_caches(tmp_path, monkeypatch):
    """A node's checkpoint plan is computed once and persisted (design doc
    Sec11, extended to per-node scope); a second call for the same
    task_index must read the sidecar rather than re-running the
    calibration probe."""
    scan, job = _multi_frame_job(
        tmp_path, "sidecar.h5", 4, cluster=True, array_task_count=2
    )
    config = job.config_data
    calls = 0
    original = reconstruction_job_module.estimate_checkpoint_plan

    def counting(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(
        reconstruction_job_module, "estimate_checkpoint_plan", counting
    )

    first = _node_checkpoint_plan(
        job, scan, config, None, total_tasks=2, task_index=0
    )
    assert calls == 1
    second = _node_checkpoint_plan(
        job, scan, config, None, total_tasks=2, task_index=0
    )
    assert calls == 1
    assert first == second


@pytest.mark.parametrize("total_tasks", [1, 2, 3])
def test_cluster_run_matches_single_node_run(tmp_path, total_tasks):
    """Splitting the same scan across a cluster array must produce the
    same output as running it on a single node -- the central scientific
    claim of the checkpoint architecture's cluster extension (design doc
    Sec13): frame-range splitting must not change results."""
    frame_count = 6

    single_root = tmp_path / "single"
    single_root.mkdir()
    _single_scan, single_job = _multi_frame_job(
        single_root, "single-result.h5", frame_count
    )
    single_path = single_root / "job.json"
    write_job(single_job, single_path)
    single_result = run_job(single_path)
    assert single_result["status"] == "complete"

    cluster_root = tmp_path / f"cluster-{total_tasks}"
    cluster_root.mkdir()
    _cluster_scan, cluster_job = _multi_frame_job(
        cluster_root,
        "cluster-result.h5",
        frame_count,
        cluster=True,
        array_task_count=total_tasks,
    )
    cluster_path = cluster_root / "job.json"
    write_job(cluster_job, cluster_path)
    for task_index in range(total_tasks):
        run_cluster_map_task(
            cluster_path,
            task_index,
            total_tasks=total_tasks,
            cpus=1,
            memory_bytes=64 * 1024 * 1024,
        )
    cluster_result = run_cluster_finalize(
        cluster_path,
        total_tasks=total_tasks,
        cpus=1,
        memory_bytes=64 * 1024 * 1024,
    )
    assert cluster_result["status"] == "complete"

    with h5py.File(single_result["output_path"], "r") as single_h5:
        with h5py.File(cluster_result["output_path"], "r") as cluster_h5:
            single_group = single_h5["entry/reconstruction/results/q_lab"]
            cluster_group = cluster_h5["entry/reconstruction/results/q_lab"]
            for dataset in ("intensity", "variance", "weight", "contributors"):
                np.testing.assert_array_equal(
                    np.nan_to_num(single_group[dataset][()]),
                    np.nan_to_num(cluster_group[dataset][()]),
                    err_msg=dataset,
                )


def test_cluster_node_resumes_without_remapping_completed_frames(
    tmp_path, monkeypatch
):
    """Interrupting one node mid-map, then re-running just that node, must
    not reload frames whose checkpoint already fully covers them -- the
    same invariant as single-node resume, scoped to one node's slice."""
    scan, job = _multi_frame_job(
        tmp_path, "cluster-resume.h5", 4, cluster=True, array_task_count=2
    )
    job_path = tmp_path / "job.json"
    write_job(job, job_path)

    loaded_frames = []
    decision_lock = threading.Lock()
    first_claimed = threading.Event()
    original_get_raw_img = SimulationScan.get_raw_img

    def failing_after_first(self, index):
        with decision_lock:
            loaded_frames.append(index)
            allowed = not first_claimed.is_set()
            if allowed:
                first_claimed.set()
        if not allowed:
            raise RuntimeError("simulated crash mid-node")
        return original_get_raw_img(self, index)

    monkeypatch.setattr(SimulationScan, "get_raw_img", failing_after_first)
    with pytest.raises(RuntimeError, match="simulated crash"):
        run_cluster_map_task(
            job_path,
            0,
            total_tasks=2,
            cpus=1,
            memory_bytes=64 * 1024 * 1024,
        )

    # Node 0 owns 2 of the 4 frames; both were attempted (racing for
    # "first"), exactly one completed before the other's crash aborted
    # this node's run.
    assert len(loaded_frames) == 2
    loaded_frames.clear()

    def counting_get_raw_img(self, index):
        loaded_frames.append(index)
        return original_get_raw_img(self, index)

    monkeypatch.setattr(SimulationScan, "get_raw_img", counting_get_raw_img)

    result = run_cluster_map_task(
        job_path, 0, total_tasks=2, cpus=1, memory_bytes=64 * 1024 * 1024
    )

    assert result["status"] == "complete"
    # Only the one frame whose checkpoint was not yet resumable was
    # reloaded -- not both of node 0's frames from scratch.
    assert len(loaded_frames) == 1


def test_cluster_finalize_raises_when_a_node_is_incomplete(tmp_path):
    """Finalize must refuse to produce output while any node's planned
    checkpoints are missing, and must not mutate job.status."""
    scan, job = _multi_frame_job(
        tmp_path, "incomplete.h5", 4, cluster=True, array_task_count=2
    )
    job_path = tmp_path / "job.json"
    write_job(job, job_path)

    # Only node 0 ever runs; node 1 never does.
    run_cluster_map_task(
        job_path, 0, total_tasks=2, cpus=1, memory_bytes=64 * 1024 * 1024
    )

    with pytest.raises(RuntimeError, match="incomplete"):
        run_cluster_finalize(
            job_path, total_tasks=2, cpus=1, memory_bytes=64 * 1024 * 1024
        )

    assert job_status(job_path)["status"] != "complete"


def test_cluster_array_task_count_never_mutates_shared_job_json(tmp_path):
    """Array elements must not mutate the shared job JSON, so they may
    run fully concurrently without coordination."""
    _scan, job = _multi_frame_job(
        tmp_path, "no-mutate.h5", 2, cluster=True, array_task_count=2
    )
    job_path = tmp_path / "job.json"
    write_job(job, job_path)
    frozen_job = job_path.read_bytes()

    run_cluster_map_task(
        job_path, 0, total_tasks=2, cpus=1, memory_bytes=64 * 1024 * 1024
    )
    run_cluster_map_task(
        job_path, 1, total_tasks=2, cpus=1, memory_bytes=64 * 1024 * 1024
    )

    assert job_path.read_bytes() == frozen_job


def test_frame_groups_cut_at_the_group_size_and_at_checkpoint_boundaries():
    """A group is bounded by two things at once.

    Its own size, and the planned checkpoint it belongs to: the frames of
    a group merge inside the kernel and cannot be separated afterwards, so
    a group spanning two checkpoints could not be routed to either.
    """
    router = _CheckpointRouter(
        {"hkl": [(0, 5), (5, 12)]},
        spec_digest="test",
        checkpoint_dir="unused",
        active_budget_bytes=1024**2,
    )
    groups = _frame_groups((0, 12), router, ["hkl"], 4)

    assert groups == [[0, 1, 2, 3], [4], [5, 6, 7, 8], [9, 10, 11]]
    assert [frame for group in groups for frame in group] == list(range(12))
    assert all(len(group) <= 4 for group in groups)


def test_frame_groups_respect_every_grid_not_only_the_first():
    """Grids may be checkpointed on different boundaries, and a group has
    to be routable to all of them."""
    router = _CheckpointRouter(
        {"hkl": [(0, 8)], "q_lab": [(0, 3), (3, 8)]},
        spec_digest="test",
        checkpoint_dir="unused",
        active_budget_bytes=1024**2,
    )
    groups = _frame_groups((0, 8), router, ["hkl", "q_lab"], 8)

    assert groups == [[0, 1, 2], [3, 4, 5, 6, 7]]


def test_frame_groups_of_one_reproduce_the_per_frame_schedule():
    router = _CheckpointRouter(
        {"hkl": [(0, 4)]},
        spec_digest="test",
        checkpoint_dir="unused",
        active_budget_bytes=1024**2,
    )

    assert _frame_groups((0, 4), router, ["hkl"], 1) == [[0], [1], [2], [3]]


def test_angles_advance_monotonically_rejects_an_interlaced_order():
    """Grouping assumes frames adjacent in index are adjacent in angle.

    An interlaced scan (orgui/backend/interlacedScanLoader.py) breaks
    that: every other frame first, then the ones between them.
    """
    bounds = np.zeros((8, 2, 4), dtype=np.float64)
    sequential = np.deg2rad(0.1) * np.arange(8)
    bounds[:, 0, 1] = sequential
    bounds[:, 1, 1] = sequential

    assert _angles_advance_monotonically(bounds, list(range(8)))

    interlaced = np.concatenate([sequential[::2], sequential[1::2]])
    bounds[:, 0, 1] = interlaced
    bounds[:, 1, 1] = interlaced

    assert not _angles_advance_monotonically(bounds, list(range(8)))


def test_angles_advance_monotonically_ignores_axes_that_never_move():
    """A scan rotating one motor is judged on that motor alone; the three
    stationary angles carry floating-point noise, not direction."""
    bounds = np.zeros((6, 2, 4), dtype=np.float64)
    bounds[:, 0, 1] = np.deg2rad(0.1) * np.arange(6)
    bounds[:, 1, 1] = bounds[:, 0, 1]
    bounds[:, :, 3] = 1e-17 * np.array([1.0, -1.0, 1.0, -1.0, 1.0, -1.0])[:, None]

    assert _angles_advance_monotonically(bounds, list(range(6)))


def test_router_counts_a_group_as_all_of_its_frames(tmp_path):
    """The checkpoint countdown and ``frames_covered`` count frames, not
    route() calls -- a group that declared one frame would leave its
    checkpoint permanently short and unresumable."""
    router = _CheckpointRouter(
        {"hkl": [(0, 6)]},
        spec_digest="test",
        checkpoint_dir=tmp_path / "checkpoints",
        active_budget_bytes=1024**3,
    )
    batch = {
        "chunk_id": np.array([0], dtype=np.uint32),
        "local_voxel_id": np.array([0], dtype=np.uint32),
        "weighted_intensity": np.array([1.0]),
        "weighted_variance": np.array([1.0]),
        "weight": np.array([1.0]),
        "contributors": np.array([1], dtype=np.uint32),
    }

    router.route("hkl", 0, dict(batch), frames=3)
    assert not router.written, "half the checkpoint's frames are still missing"

    router.route("hkl", 3, dict(batch), frames=3)
    assert len(router.written) == 1

    with h5py.File(router.written[0], "r") as handle:
        assert int(handle.attrs["frames_covered"]) == 6


def test_router_rejects_a_group_straddling_a_checkpoint_boundary(tmp_path):
    router = _CheckpointRouter(
        {"hkl": [(0, 4), (4, 8)]},
        spec_digest="test",
        checkpoint_dir=tmp_path / "checkpoints",
        active_budget_bytes=1024**3,
    )
    batch = {name: np.zeros(0) for name in ("weighted_intensity",)}
    batch["chunk_id"] = np.zeros(0, dtype=np.uint32)

    with pytest.raises(ValueError, match="straddles"):
        router.route("hkl", 2, batch, frames=4)


def test_frame_parallelism_charges_a_group_for_its_resident_frames():
    """A group holds every one of its corrected frames at once, and its
    native call covers frames x tile samples -- but it emits one merged
    record batch, not one per frame, so only the correction term
    multiplies."""
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(1.0, 1.0, 1.0),
        step=(0.1, 0.1, 0.1),
        frame="lab",
        chunk_shape=(8, 8, 8),
    )
    tiles = [(0, 512, 0, 512)]
    budget = 32 * 1024**3

    def worker_bytes(frames_per_group):
        spec = _ReconstructionSpec(
            grids=(grid,),
            max_depth=0,
            threads=8,
            frames_per_group=frames_per_group,
        )
        _workers, _threads, per_worker, _accumulation = _frame_parallelism(
            spec, tiles, budget, threads_per_image=1
        )
        return per_worker

    pixels = 512 * 512
    native_per_pixel = 128 + 2 * 40
    record_term = pixels * _FRAME_RECORD_BYTES_PER_PIXEL
    for frames_per_group in (1, 2, 8):
        expected = (
            pixels * frames_per_group * native_per_pixel
            + pixels * frames_per_group * _PYTHON_CORRECTION_BYTES_PER_PIXEL
            + record_term
        )
        assert worker_bytes(frames_per_group) >= expected


def _map_with_group_size(tmp_path, frames_per_group, frame_count=8):
    """Map one range end to end at a given group size; return its output.

    Goes through _map_pending_ranges, so above one frame per group this
    exercises the streamed scheduler in full -- dispatcher, prepare pool,
    compute pool and gate -- rather than _map_frame_group alone.
    """
    root = tmp_path / f"group{frames_per_group}"
    root.mkdir()
    scan, job = _multi_frame_job(root, "result.h5", frame_count)
    # One checkpoint over every frame, so groups are cut by their own size
    # rather than by a boundary between every pair of frames.
    job.checkpoint_plan = {"q_lab": [[0, frame_count]]}
    config = job.config_data
    spec = dataclasses.replace(
        job.internal_spec(), frames_per_group=frames_per_group
    )
    boundaries = {"q_lab": [(0, frame_count)]}
    checkpoint_dir = root / "checkpoints"
    router = _CheckpointRouter(
        boundaries,
        spec_digest=job.digest,
        checkpoint_dir=checkpoint_dir,
        active_budget_bytes=64 * 1024**2,
    )
    reconstruction_job_module._map_pending_ranges(
        spec,
        scan,
        config,
        scan.exposure_angle_bounds(config, fallback=job.angle_fallback),
        [(0, 2, 0, 2)],
        [(0, frame_count)],
        router,
        correction_pipeline=reconstruction_job_module._correction_pipeline(
            config, scan, reconstruction_job_module._load_assets(job), {}
        ),
        effective_memory=64 * 1024**2,
        threads_per_image=1,
        accumulation_budget_bytes=None,
        total_images=frame_count,
        completed_images=0,
        progress=None,
    )
    contributors: dict[tuple[int, int], int] = {}
    frames_covered = 0
    for path in sorted(checkpoint_dir.glob("q_lab/ckpt*.h5")):
        with h5py.File(path, "r") as handle:
            frames_covered += int(handle.attrs["frames_covered"])
            for chunk, voxel, count in zip(
                handle["chunk_id"][:],
                handle["local_voxel_id"][:],
                handle["contributors"][:],
            ):
                key = (int(chunk), int(voxel))
                contributors[key] = contributors.get(key, 0) + int(count)
    return frames_covered, contributors


def test_streamed_group_scheduler_matches_the_per_frame_pipeline(tmp_path):
    """The grouped scheduler is a different pipeline -- one all-threads
    call at a time, correction hoisted into the prepare pool -- so it has
    to be shown to produce the per-frame pipeline's answer, not just to
    run."""
    frames_covered, reference = _map_with_group_size(tmp_path, 1)
    assert frames_covered == 8
    assert reference, "the fixture must reach at least one voxel"

    for frames_per_group in (2, 4):
        covered, contributors = _map_with_group_size(
            tmp_path, frames_per_group
        )
        # Every frame is still counted exactly once, which is what makes
        # the checkpoint resumable.
        assert covered == 8
        assert contributors == reference


def test_group_pipeline_layout_trades_concurrency_for_group_size():
    """A group buffer and a group's native working set both scale with the
    group size, so a larger group buys read-ahead out of concurrency. The
    thread budget is then split between whatever concurrency survives."""
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(1.0, 1.0, 1.0),
        step=(0.1, 0.1, 0.1),
        frame="lab",
        chunk_shape=(8, 8, 8),
    )
    tiles = [(0, 256, 0, 1024), (256, 512, 0, 1024)]
    budget = 8 * 1024**3
    previous_workers = None
    for frames_per_group in (1, 2, 4, 8, 16):
        spec = _ReconstructionSpec(
            grids=(grid,),
            max_depth=0,
            threads=24,
            frames_per_group=frames_per_group,
        )
        workers, threads, depth = _group_pipeline_layout(spec, tiles, budget, 4)
        assert workers >= 1
        assert threads >= 1
        assert workers * threads <= 24
        # Always at least one group of read-ahead beyond the calls in
        # flight, or there is no double-buffering and every call stalls at
        # its group boundary.
        assert depth > workers
        if previous_workers is not None:
            assert workers <= previous_workers
        previous_workers = workers


def test_group_pipeline_layout_never_returns_zero_workers_on_a_tiny_budget():
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(1.0, 1.0, 1.0),
        step=(0.1, 0.1, 0.1),
        frame="lab",
        chunk_shape=(8, 8, 8),
    )
    spec = _ReconstructionSpec(
        grids=(grid,), max_depth=0, threads=8, frames_per_group=16
    )

    workers, threads, depth = _group_pipeline_layout(
        spec, [(0, 2048, 0, 2048)], 1024**2, 4
    )

    assert (workers, threads) == (1, 8)
    assert depth > workers


def _layout_spec(frames_per_group=None, threads=24):
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(1.0, 1.0, 1.0),
        step=(0.1, 0.1, 0.1),
        frame="lab",
        chunk_shape=(8, 8, 8),
    )
    return _ReconstructionSpec(
        grids=(grid,),
        max_depth=0,
        threads=threads,
        frames_per_group=frames_per_group,
    )


def test_choose_frames_per_group_stops_where_frames_stop_overlapping():
    """Grouping merges frames that reach the same voxels. Once a frame
    advances a whole voxel they tile instead, and the group buffer costs
    concurrency for nothing."""
    tiles = [(0, 256, 0, 1024), (256, 512, 0, 1024)]
    budget = 8 * 1024**3

    assert (
        _choose_frames_per_group(_layout_spec(), tiles, budget, 0.7, 4) > 1
    )
    assert _choose_frames_per_group(_layout_spec(), tiles, budget, 1.0, 4) == 1
    assert _choose_frames_per_group(_layout_spec(), tiles, budget, 3.2, 4) == 1


def test_choose_frames_per_group_keeps_concurrency_rather_than_density():
    """Records fall monotonically with the group size, so density alone
    would pick the largest group memory allows -- which measured *slowest*
    on the reference job, because the group crowds out concurrent calls.
    The choice must stop where concurrency does.

    Sized so the floor actually binds inside the candidate list, the way
    it does on a real detector: a tile small enough for every candidate to
    fit would test nothing.
    """
    tiles = [(0, 1024, 0, 2048)]
    spec = _layout_spec()
    budget = 8 * 1024**3

    chosen = _choose_frames_per_group(spec, tiles, budget, 0.7, 4)

    assert 1 < chosen < max(_GROUP_SIZE_CANDIDATES)
    workers, _threads, _depth = _group_pipeline_layout(
        dataclasses.replace(spec, frames_per_group=chosen), tiles, budget, 4
    )
    assert workers >= 3
    # The next size up must be the one that broke the floor, or the choice
    # stopped early and left throughput on the table.
    larger = _group_pipeline_layout(
        dataclasses.replace(spec, frames_per_group=chosen * 2), tiles, budget, 4
    )[0]
    assert larger < 3


def test_choose_frames_per_group_falls_back_to_one_on_a_tight_budget():
    """A budget that cannot afford concurrent group calls must not group:
    one call with every thread is the configuration that measured worst."""
    tiles = [(0, 2048, 0, 2048)]

    assert (
        _choose_frames_per_group(_layout_spec(), tiles, 512 * 1024**2, 0.5, 4)
        == 1
    )


def test_frame_advance_voxels_is_zero_when_the_grid_does_not_rotate():
    """A lab-frame grid is the degenerate case worth pinning: Q in the lab
    frame does not depend on the sample angles, so every frame reaches the
    same voxels however far the sample turns. Zero travel is the correct
    answer, and it is the case grouping helps most."""
    grid = ReconstructionGrid(
        minimum=(-20.0, -20.0, -20.0),
        maximum=(20.0, 20.0, 20.0),
        step=(0.05, 0.05, 0.05),
        frame="lab",
        chunk_shape=(8, 8, 8),
    )
    config = _config()
    spec = _ReconstructionSpec(grids=(_GridSpec(**grid.__dict__),), max_depth=0)
    kernel = reconstruction_job_module._kernel_for_grid(
        spec, spec.grids[0], config.ub_calculator, threads=1
    )
    rays = reconstruction_job_module._detector_corner_rays(
        config.detector, (0, 2, 0, 2)
    )
    bounds = np.zeros((3, 2, 4), dtype=np.float64)
    bounds[:, 0, 1] = np.deg2rad(5.0) * np.arange(3)
    bounds[:, 1, 1] = bounds[:, 0, 1]

    assert reconstruction_job_module._frame_advance_voxels(
        kernel, spec.grids[0], rays, bounds, [(0, 1)]
    ) == 0.0


def test_frame_advance_voxels_scales_with_the_angular_step():
    """The probe has to measure the job, so a scan that turns twice as far
    per frame must read as twice the travel."""
    grid = ReconstructionGrid(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(0.002, 0.002, 0.002),
        frame="hkl",
        chunk_shape=(8, 8, 8),
    )
    config = _config()
    spec = _ReconstructionSpec(grids=(_GridSpec(**grid.__dict__),), max_depth=0)
    kernel = reconstruction_job_module._kernel_for_grid(
        spec, spec.grids[0], config.ub_calculator, threads=1
    )
    rays = reconstruction_job_module._detector_corner_rays(
        config.detector, (0, 2, 0, 2)
    )

    def advance(step_degrees):
        bounds = np.zeros((3, 2, 4), dtype=np.float64)
        angles = np.deg2rad(step_degrees) * np.arange(3)
        bounds[:, 0, 1] = angles
        bounds[:, 1, 1] = angles
        return reconstruction_job_module._frame_advance_voxels(
            kernel, spec.grids[0], rays, bounds, [(0, 1)]
        )

    small = advance(0.1)
    large = advance(0.2)
    assert small > 0.0
    assert large == pytest.approx(2.0 * small, rel=0.05)


def _large_detector_job(tmp_path, depth):
    scan = SimulationScan((2, 2), 0.0, 1.0, 8)
    assets = tmp_path / f"scratch{depth}" / "job-assets.nxs"
    assets.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(assets, "w") as h5file:
        h5file.attrs["orgui_job_assets"] = 1
    config = _config()
    config.detector.detector.shape = (2527, 2463)
    config.detector.detector.max_shape = (2527, 2463)
    grid = ReconstructionGrid(
        minimum=(-20.0, -20.0, -20.0),
        maximum=(20.0, 20.0, 20.0),
        step=(20.0, 20.0, 20.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    scan_reference = ScanReference.from_scan(scan).to_dict()
    job = ReconstructionJob(
        config=config_data_to_json(config),
        scan_reference=scan_reference,
        grids=[grid.__dict__],
        scratch_path=str(assets.parent),
        output_path=str(tmp_path / f"result{depth}.h5"),
        compression="Raw",
        assets_path=str(assets),
        assets_sha256=sha256(assets.read_bytes()).hexdigest(),
        source_fingerprint_sha256="0" * 64,
        threads_per_image=None,
        accumulation_budget_bytes=None,
        advanced_depth=depth,
        runtime_threads=24,
        runtime_memory_bytes=10 * 1000**3,
        frame_batch=4,
    )
    return scan, job, config


def test_detector_bands_do_not_collapse_at_high_adaptive_depth(tmp_path):
    """Band height is bounded by the record ceiling, not by the worst-case
    leaf count.

    ``8**depth`` bytes per pixel is what this site was left with when the
    two memory prechecks were corrected to the record ceiling. It is
    enormously conservative: at depth 5 it asks for 2.6 MB per pixel,
    which bands a Pilatus 6M into one row per band -- 2527 native calls
    per frame instead of six. Thin bands are not free, and one row is not
    a band at all.
    """
    band_counts = {}
    for depth in (0, 2, 5):
        scan, job, config = _large_detector_job(tmp_path, depth)
        _ranges, tiles = reconstruction_job_module._execution_layout(
            job, scan, config
        )
        band_counts[depth] = len(tiles)
        # Full detector width, so a tile is a contiguous slice of the frame.
        assert all(
            (column_start, column_stop) == (0, 2463)
            for _rs, _re, column_start, column_stop in tiles
        )

    assert band_counts[0] == band_counts[2] == band_counts[5]
    assert band_counts[5] < 20

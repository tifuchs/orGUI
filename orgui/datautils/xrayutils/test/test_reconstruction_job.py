"""End-to-end tests for the centralized reconstruction job."""

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
import orgui.reconstruction_job as reconstruction_job_module
from orgui.reconstruction_job import (
    ReconstructionGrid,
    ReconstructionJob,
    _node_checkpoint_plan,
    _node_excluded_frames,
    job_status,
    reconstruction_execution_settings,
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

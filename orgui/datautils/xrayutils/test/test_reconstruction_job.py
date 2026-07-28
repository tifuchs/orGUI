"""End-to-end tests for the centralized reconstruction job."""

from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import threading
import time

import h5py
import numpy as np
import pytest

from orgui.app.config_data import ConfigData
from orgui.app.database import config_data_to_json
from orgui.backend.scans import ScanReference, SimulationScan
from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg
from orgui.reconstruction_job import (
    ReconstructionGrid,
    ReconstructionJob,
    job_status,
    reconstruction_execution_settings,
    run_cluster_finalize,
    run_cluster_map_task,
    run_job,
    write_job,
)


pytest.importorskip("pyarrow")
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


def test_central_job_runs_resumes_and_cleans_verified_scratch(
    tmp_path, monkeypatch
):
    scan = SimulationScan((2, 2), 0.0, 1.0, 2)
    assets = tmp_path / "scratch" / "job-assets.nxs"
    assets.parent.mkdir()
    with h5py.File(assets, "w") as h5file:
        h5file.attrs["orgui_job_assets"] = 1
    assets_digest = sha256(assets.read_bytes()).hexdigest()
    grid = ReconstructionGrid(
        minimum=(-20.0, -20.0, -20.0),
        maximum=(20.0, 20.0, 20.0),
        step=(20.0, 20.0, 20.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    scan_reference = ScanReference.from_scan(scan).to_dict()
    source_digest = sha256(
        json.dumps(
            scan_reference, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    job = ReconstructionJob(
        config=config_data_to_json(_config()),
        scan_reference=scan_reference,
        grids=[grid.__dict__],
        scratch_path=str(assets.parent),
        output_path=str(tmp_path / "result.h5"),
        compression="Raw",
        assets_path=str(assets),
        assets_sha256=assets_digest,
        source_fingerprint_sha256=source_digest,
        threads_per_image=1,
        accumulation_budget_bytes=None,
        runtime_threads=2,
        runtime_memory_bytes=64 * 1024 * 1024,
        tile_shape=(1, 1),
        frame_batch=1,
        expected_map_tasks=2,
    )
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
    assert not assets.exists()
    with h5py.File(first["output_path"], "r") as h5file:
        assert "entry/configuration" in h5file
        group = h5file["entry/reconstruction/results/q_lab"]
        assert np.nanmax(group["intensity"]) == 10.0
        assert np.nanmax(group["contributors"]) == 8
    assert job_status(job_path)["output_sha256"] == first["output_sha256"]


def test_cluster_array_tasks_do_not_mutate_job_and_finalizer_collects_them(
    tmp_path,
):
    """Independent array tasks publish manifests before one finalizer writes."""
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
        output_path=str(tmp_path / "cluster-result.h5"),
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
        runtime_threads=2,
        runtime_memory_bytes=64 * 1024 * 1024,
        tile_shape=(1, 1),
        frame_batch=1,
        expected_map_tasks=2,
    )
    job_path = tmp_path / "cluster-job.json"
    write_job(job, job_path)
    frozen_job = job_path.read_bytes()

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(
                run_cluster_map_task,
                job_path,
                task_index,
                cpus=1,
                memory_bytes=64 * 1024 * 1024,
            )
            for task_index in range(2)
        ]
        first, second = (future.result() for future in futures)
    assert job_path.read_bytes() == frozen_job
    assert first["frame_range"] == [0, 1]
    assert second["frame_range"] == [1, 2]
    assert job_status(job_path)["map_tasks"] == {
        "completed": 2,
        "pending": 0,
        "total": 2,
    }
    assert job_path.read_bytes() == frozen_job

    result = run_cluster_finalize(
        job_path,
        cpus=2,
        memory_bytes=128 * 1024 * 1024,
    )

    assert result["status"] == "complete"
    assert Path(result["output_path"]).exists()
    assert job_status(job_path)["map_tasks"]["completed"] == 2

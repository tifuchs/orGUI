"""Tests for the job-level checkpoint-plan calibration probe (design doc
Sec5/Sec6): :func:`orgui.reconstruction_job.estimate_checkpoint_plan`.
"""

from __future__ import annotations

import numpy as np
import pytest

from orgui.reconstruction_job import ReconstructionGrid, estimate_checkpoint_plan


pytest.importorskip("orgui.datautils.xrayutils._reciprocal_reconstruction_cpp")


class _FakeCalibrationDetector:
    """A small, non-degenerate detector double -- large enough for the
    probe's stratified tile sampling to be meaningful, unlike the 1x1
    double used elsewhere in this test package for unrelated geometry
    tests."""

    detector = type("Detector", (), {"shape": (32, 28)})()

    def primBeamPoints(self, rows, columns):
        return np.zeros_like(rows), np.zeros_like(columns)


class _FakeCalibrationUB:
    def getUB(self):
        return np.eye(3)

    def getU(self):
        return np.eye(3)

    def getK(self):
        return 1.0


class _FakeCorrections:
    def __init__(self, excluded_frames=()):
        self.excluded_frames = tuple(excluded_frames)


class _FakeCalibrationConfig:
    def __init__(self, excluded_frames=()):
        self.detector = _FakeCalibrationDetector()
        self.ub_calculator = _FakeCalibrationUB()
        self.corrections = _FakeCorrections(excluded_frames)


class _FakeCalibrationScan:
    def __init__(self, frame_count):
        self._frame_count = frame_count

    def __len__(self):
        return self._frame_count

    def exposure_angle_bounds(self, config, fallback="stationary"):
        return np.zeros((len(self), 2, 4), dtype=np.float64)


def _grid(name="hkl"):
    return ReconstructionGrid(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(0.1, 0.1, 0.1),
        frame="hkl",
        name=name,
    )


def test_estimate_checkpoint_plan_returns_one_estimate_per_grid():
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=10)

    result = estimate_checkpoint_plan(
        config,
        scan,
        [_grid("hkl"), _grid("q_lab")],
        max_depth=0,
        threads=1,
        ram_budget_bytes=1024**3,
        checkpoint_count=10,
        budget_seconds=0.02,
    )

    assert set(result["per_grid"]) == {"hkl", "q_lab"}
    for estimate in result["per_grid"].values():
        assert estimate["job_data_bytes_estimate"] >= 0.0
        assert estimate["files_per_job"] >= 10
    assert result["files_total"] == sum(
        estimate["files_per_job"] for estimate in result["per_grid"].values()
    )


def test_estimate_checkpoint_plan_respects_checkpoint_floor():
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=5)

    result = estimate_checkpoint_plan(
        config,
        scan,
        [_grid()],
        max_depth=0,
        threads=1,
        ram_budget_bytes=1024**4,  # deliberately huge: data volume never wins
        checkpoint_count=7,
        budget_seconds=0.02,
    )

    assert result["per_grid"]["hkl"]["files_per_job"] == 7


def test_estimate_checkpoint_plan_is_deterministic_for_a_fixed_frame_choice():
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=8)

    kwargs = dict(
        max_depth=0,
        threads=1,
        ram_budget_bytes=1024**3,
        checkpoint_count=10,
        budget_seconds=0.02,
    )
    first = estimate_checkpoint_plan(config, scan, [_grid()], **kwargs)
    second = estimate_checkpoint_plan(config, scan, [_grid()], **kwargs)

    # The plan itself -- what the job is actually scheduled around -- has to
    # be reproducible. The record-volume estimate underneath it is measured
    # against a wall-clock budget, so how many pixels the probe manages to
    # sample varies slightly between calls even with the sample positions
    # pinned; requiring that figure to be bit-identical would be requiring
    # the machine to be idle, not the algorithm to be deterministic.
    assert first["files_total"] == second["files_total"]
    assert (
        first["per_grid"]["hkl"]["files_per_job"]
        == second["per_grid"]["hkl"]["files_per_job"]
    )
    assert first["per_grid"]["hkl"]["job_data_bytes_estimate"] == pytest.approx(
        second["per_grid"]["hkl"]["job_data_bytes_estimate"], rel=0.1
    )


def test_estimate_checkpoint_plan_rejects_all_frames_excluded():
    config = _FakeCalibrationConfig(excluded_frames=range(4))
    scan = _FakeCalibrationScan(frame_count=4)

    with pytest.raises(ValueError, match="No included frames"):
        estimate_checkpoint_plan(
            config,
            scan,
            [_grid()],
            max_depth=0,
            threads=1,
            ram_budget_bytes=1024**3,
            checkpoint_count=10,
        )


def test_estimate_checkpoint_plan_rejects_mask_shape_mismatch():
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=4)

    with pytest.raises(ValueError, match="mask shape"):
        estimate_checkpoint_plan(
            config,
            scan,
            [_grid()],
            max_depth=0,
            threads=1,
            ram_budget_bytes=1024**3,
            checkpoint_count=10,
            mask=np.zeros((5, 5), dtype=bool),
        )

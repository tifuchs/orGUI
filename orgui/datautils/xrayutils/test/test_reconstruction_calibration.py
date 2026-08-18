"""Tests for the job-level checkpoint-plan calibration probe (design doc
Sec5/Sec6): :func:`orgui.reconstruction_job.estimate_checkpoint_plan`.
"""

from __future__ import annotations

import numpy as np
import pytest

from orgui.datautils.xrayutils.reconstruction import (
    _MINIMUM_CONVERGENCE_TILES,
    _calibration_probe,
)
from orgui.reconstruction_job import ReconstructionGrid, estimate_checkpoint_plan


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


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


def _unreachable_kernel():
    """A kernel whose grid sits far from anything the rays reach."""
    return native.ReconstructionKernel(
        np.full(3, 500.0),
        np.full(3, 0.01),
        np.full(3, 8, dtype=np.int64),
        np.full(3, 8, dtype=np.int64),
        "hkl",
        1.0,
        np.ascontiguousarray(np.eye(3)),
        np.ascontiguousarray(np.eye(3)),
        2,
        1,
        4096,
        # Generous: the kernel's own per-call memory precheck bounds how
        # large a sample tile may be, and a tight budget would cap the
        # tiles well below the sample ceiling these tests are about.
        4 * 1024**3,
    )


def _constant_rays(rows, columns):
    rays = np.zeros((rows + 1, columns + 1, 3), dtype=np.float64)
    rays[..., 1] = 1.0
    return np.ascontiguousarray(rays)


def _flat_coverage(scan, points):
    """A frame-tagged coverage sample shaped like the real sampler's."""
    frames = np.arange(len(scan), dtype=np.int64)
    return frames, np.zeros((frames.size, points, 3), dtype=np.float64)


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


def test_estimate_checkpoint_plan_accepts_a_fractional_memory_budget():
    """A memory budget that is not a whole number of bytes still estimates.

    The GUI's memory setting is in MiB and is a float, so the live checkpoint
    estimate reaches the native kernel with ``maxMemory * 1024 * 1024`` --
    also a float, which the kernel's integer constructor argument rejects
    outright. ``prepare_job`` rounds and so always worked; the dialog did
    not, and reported "incompatible constructor arguments" instead of a file
    count.
    """
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=6)

    result = estimate_checkpoint_plan(
        config,
        scan,
        [_grid()],
        max_depth=2,
        threads=1,
        ram_budget_bytes=22272.498827934265 * 1024 * 1024,
        checkpoint_count=4,
        budget_seconds=0.02,
    )

    assert result["per_grid"]["hkl"]["files_per_job"] >= 4


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


def test_probe_does_not_treat_finding_nothing_as_convergence():
    """Empty tiles agree perfectly; that is not a converged measurement.

    The early stop compares the sampled densities' relative standard error
    against a target, and a run of tiles that found no record has a spread
    of exactly zero. That read as "converged" and stopped the probe after
    the minimum number of tiles, on precisely the grids it had learned
    nothing about -- so a small output volume the sample had not yet found
    was reported as containing no data at all.
    """
    # Large enough that the convergence check is reachable: on a small
    # frame the probe hits the pixel ceiling first and the early stop never
    # gets a chance to fire, which would make this test pass either way.
    rows = columns = 512
    mask = np.zeros((rows, columns), dtype=bool)
    kernel = _unreachable_kernel()

    result = _calibration_probe(
        kernel,
        mask,
        _constant_rays(rows, columns),
        np.zeros(4),
        np.zeros(4),
        budget_seconds=0.2,
        rng=np.random.default_rng(0),
    )

    assert result["records_per_pixel"] == 0.0
    assert result["sampled_tiles"] >= _MINIMUM_CONVERGENCE_TILES
    # Having found nothing, it keeps looking until the budget or the
    # detector runs out, rather than stopping at the convergence minimum
    # with a few hundred pixels sampled. The ceiling is the whole frame --
    # the only honest one for a volume reached by a thin locus of pixels.
    assert result["sampled_pixels"] >= rows * columns


def test_probe_sample_ceiling_is_the_detector_not_a_constant():
    """The ceiling scales with the frame rather than sitting at 2 megapixels.

    A small output volume is reached by a thin locus of pixels, so a sample
    stopped at a fixed fraction of a large frame finds it or misses it
    depending on where the tiles land -- which made the estimate depend on
    the draw rather than on the job. Deliberately uses a frame larger than
    the old constant, since that is the only size at which it was binding.
    """
    rows = columns = 1500
    assert rows * columns > 2_000_000
    mask = np.zeros((rows, columns), dtype=bool)

    result = _calibration_probe(
        _unreachable_kernel(),
        mask,
        _constant_rays(rows, columns),
        np.zeros(4),
        np.zeros(4),
        budget_seconds=2.0,
        rng=np.random.default_rng(0),
    )

    assert result["sampled_pixels"] > 2_000_000


def test_estimate_scales_a_grid_by_the_frames_that_reach_it():
    """A volume few frames reach is not credited with the whole scan.

    Automatic volume selection produces grids around single features, and a
    crystal truncation rod at a fixed ``(H, K)`` is crossed over a narrow
    range of sample rotations. Multiplying a reaching frame's record density
    by every frame in the scan would overstate such a grid by the reciprocal
    of that range.
    """
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=64)
    grids = [_grid()]
    common = dict(
        max_depth=0,
        threads=1,
        ram_budget_bytes=1024**2,
        checkpoint_count=1,
        budget_seconds=0.02,
    )
    frames, cloud = _flat_coverage(scan, points=8)

    inside = np.zeros_like(cloud)  # every sample inside the grid
    every_frame = estimate_checkpoint_plan(
        config, scan, grids, coverage=(frames, inside), **common
    )
    outside = np.full_like(cloud, 500.0)  # no sample inside, except one frame
    outside[0] = 0.0
    one_frame = estimate_checkpoint_plan(
        config, scan, grids, coverage=(frames, outside), **common
    )

    assert every_frame["per_grid"]["hkl"]["reaching_frame_fraction"] == 1.0
    assert one_frame["per_grid"]["hkl"]["reaching_frame_fraction"] == pytest.approx(
        1.0 / frames.size
    )
    assert (
        one_frame["per_grid"]["hkl"]["job_data_bytes_estimate"]
        < every_frame["per_grid"]["hkl"]["job_data_bytes_estimate"]
    )


def test_estimate_never_claims_a_grid_is_reached_by_no_frame():
    """A cloud that misses a volume leaves it unresolved, not empty.

    The coverage sample is a cloud of points, so a volume it never lands in
    may still be grazed by frames between the samples. Crediting such a grid
    with zero frames would zero its whole estimate on the strength of a
    sample that was never fine enough to say so.
    """
    config = _FakeCalibrationConfig()
    scan = _FakeCalibrationScan(frame_count=32)
    frames, cloud = _flat_coverage(scan, points=8)

    result = estimate_checkpoint_plan(
        config,
        scan,
        [_grid()],
        max_depth=0,
        threads=1,
        ram_budget_bytes=1024**2,
        checkpoint_count=1,
        budget_seconds=0.02,
        coverage=(frames, np.full_like(cloud, 500.0)),
    )

    fraction = result["per_grid"]["hkl"]["reaching_frame_fraction"]
    assert 0.0 < fraction <= 1.0 / frames.size


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

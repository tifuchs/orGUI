"""Tests for geometry-matched reciprocal-space grid sampling."""

from types import SimpleNamespace

import numpy as np

import orgui.datautils.xrayutils.reconstruction as reconstruction_module
from orgui.reconstruction_job import estimate_geometry_steps


class _Detector:
    detector = SimpleNamespace(shape=(4, 5))

    @staticmethod
    def primBeamPoints(rows, columns):
        return np.zeros_like(rows), np.zeros_like(columns)


class _UBCalculator:
    @staticmethod
    def getK():
        return 1.0

    @staticmethod
    def getUB():
        return np.eye(3)

    @staticmethod
    def getU():
        return np.eye(3)


class _Scan:
    def __len__(self):
        return 2

    @staticmethod
    def exposure_angle_bounds(config, fallback):
        bounds = np.zeros((2, 2, 4), dtype=np.float64)
        bounds[1, :, 0] = 1.0
        return bounds


class _Kernel:
    def __init__(self, *args):
        pass

    @staticmethod
    def coordinate(rays, start, end, row, column, u, v, t):
        angle = start[0] + t * (end[0] - start[0])
        return np.asarray((2.0 * u, 3.0 * v, 4.0 * angle))


def test_local_jacobian_propagates_uniform_pixel_and_scan_widths(monkeypatch):
    """Jacobian columns propagate to one-sigma axis projections."""
    monkeypatch.setattr(
        reconstruction_module,
        "_native_module",
        lambda: SimpleNamespace(ReconstructionKernel=_Kernel),
    )
    config = SimpleNamespace(
        detector=_Detector(),
        ub_calculator=_UBCalculator(),
        corrections=SimpleNamespace(excluded_frames=()),
    )

    steps = estimate_geometry_steps(
        config,
        _Scan(),
        percentile=37.0,
        detector_samples=2,
        frame_samples=2,
    )

    np.testing.assert_allclose(
        steps,
        np.asarray((2.0, 3.0, 4.0)) / np.sqrt(12.0),
    )

"""Regression tests for the per-pixel detector correction factors.

:mod:`orgui.datautils.xrayutils.corrections.detector` owns the solid-angle and
polarization arrays and the region reduction of the solid-angle part that a
structure factor has to divide back out. See
``doc/design/ctr_structure_factor_scale.md`` finding F6.
"""

import numpy as np
import pytest

pyFAI = pytest.importorskip("pyFAI")

from orgui.datautils.xrayutils import DetectorCalibration  # noqa: E402
from orgui.datautils.xrayutils.corrections import detector as detector_corrections  # noqa: E402

#: Pixel size and distance of the calibrated test detector, in meter.
PIXEL, DIST = 172e-6, 0.3
SHAPE = (619, 487)


def _calibrated_detector():
    """A real, calibrated area detector, close in so obliquity is visible."""
    det = DetectorCalibration.Detector2D_SXRD()
    det.detector = pyFAI.detectors.Detector(
        pixel1=PIXEL, pixel2=PIXEL, max_shape=SHAPE
    )
    det.poni1 = SHAPE[0] * PIXEL / 2.0
    det.poni2 = SHAPE[1] * PIXEL / 2.0
    det.rot1 = det.rot2 = det.rot3 = 0.0
    det.dist = DIST
    det.set_energy(17.7)
    det.setAzimuthalReference(np.deg2rad(90.0))
    det.setPolarization(0.0, 1.0)
    det.reset()
    det._cached_array = {}
    return det


def test_the_region_mean_is_one_at_normal_incidence():
    """At the point of normal incidence the correction is nearly neutral.

    ``solidAngleArray`` is normalized to that point, so a region centred
    there cannot introduce a scale of its own beyond its own curvature: the
    off-centre pixels of a 20-pixel region at 0.3 m still reach
    :math:`\\theta = 0.33^\\circ`, worth :math:`4\\times10^{-5}`. That
    residual is the factor being real rather than an artefact, so it is
    bounded here rather than asserted away.
    """
    det = _calibrated_detector()

    got = detector_corrections.roi_mean_inverse_solid_angle(
        det, SHAPE[0] / 2.0, SHAPE[1] / 2.0, 20.0, 20.0
    )

    np.testing.assert_allclose(got, 1.0, rtol=1e-4)
    assert got > 1.0, "the mean of 1/cos^3 over a region is never below one"


def test_the_region_mean_follows_the_obliquity():
    """It is the mean of ``1/cos^3(theta)`` over the region's pixels.

    Written out from the geometry rather than by calling the module back: the
    normalized solid angle of a flat detector is ``cos^3`` of the incidence
    angle on its face, which is what makes this factor grow away from the
    beam centre.
    """
    det = _calibrated_detector()
    row, column, rows, columns = 560.0, 240.0, 60.0, 40.0

    got = detector_corrections.roi_mean_inverse_solid_angle(
        det, row, column, rows, columns
    )

    r = np.arange(row - rows / 2, row + rows / 2)
    c = np.arange(column - columns / 2, column + columns / 2)
    dr = (r[:, None] + 0.5) * PIXEL - det.poni1
    dc = (c[None, :] + 0.5) * PIXEL - det.poni2
    theta = np.arctan(np.hypot(dr, dc) / DIST)
    expected = np.mean(1.0 / np.cos(theta) ** 3)

    np.testing.assert_allclose(got, expected, rtol=1e-6)
    assert got > 1.02, "the test region must be oblique enough to matter"


def test_the_region_mean_is_vectorized_over_a_scan():
    """orGUI resizes regions along a scan, so this is normally an array."""
    det = _calibrated_detector()
    row = np.array([310.0, 450.0, 560.0])
    size = np.array([20.0, 40.0, 60.0])

    got = detector_corrections.roi_mean_inverse_solid_angle(
        det, row, 240.0, size, size
    )

    assert got.shape == (3,)
    # Monotonic away from the beam centre, which is what makes it a shape
    # error rather than a scale error when it is left in.
    assert got[0] < got[1] < got[2]


def test_a_region_off_the_detector_is_neutral():
    """No pixels means no counts; the factor must not scale or poison them."""
    det = _calibrated_detector()

    got = detector_corrections.roi_mean_inverse_solid_angle(
        det, -500.0, -500.0, 20.0, 20.0
    )

    np.testing.assert_allclose(got, 1.0, rtol=0.0)


def test_a_region_partly_off_the_detector_uses_the_pixels_it_has():
    """Clipped at the detector edge rather than padded or rejected."""
    det = _calibrated_detector()

    edge = detector_corrections.roi_mean_inverse_solid_angle(
        det, 5.0, 240.0, 40.0, 40.0
    )

    assert np.isfinite(edge)
    assert edge > 1.0


def test_a_zero_or_negative_region_is_rejected():
    """An empty region would divide a structure factor by a meaningless mean."""
    det = _calibrated_detector()

    with pytest.raises(ValueError, match="positive height"):
        detector_corrections.roi_mean_inverse_solid_angle(
            det, 300.0, 240.0, 0.0, 20.0
        )
    with pytest.raises(ValueError, match="positive width"):
        detector_corrections.roi_mean_inverse_solid_angle(
            det, 300.0, 240.0, 20.0, -5.0
        )

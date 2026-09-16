"""Regression tests for the out-of-plane detector acceptance.

:mod:`orgui.datautils.xrayutils.corrections.acceptance` estimates the
:math:`\\Delta\\gamma` that a rocking-scan integrated intensity is
proportional to (E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532, equations 20
and 42). See ``doc/design/ctr_structure_factor_scale.md`` finding F3.
"""

import numpy as np
import pytest

pyFAI = pytest.importorskip("pyFAI")

from orgui.datautils.xrayutils import DetectorCalibration  # noqa: E402
from orgui.datautils.xrayutils.corrections import acceptance  # noqa: E402

#: Pixel size and distance of the calibrated test detector, in meter.
PIXEL, DIST = 172e-6, 1.0
SHAPE = (619, 487)


class _LinearDetector:
    """A detector whose exit angle is an exact linear function of the row.

    Lets every expectation below be written in closed form, independently of
    the pyFAI geometry, so a failure points at this module rather than at the
    detector calibration.
    """

    #: Radian per row, and per column to make the two directions separable.
    PER_ROW = 1e-4
    PER_COLUMN = 3e-5

    def surfaceAnglesPoint(self, x, y, alpha_i, gamma_arm=None, delta_arm=None):
        # First argument is pyFAI dimension 1, the row.
        arm = 0.0 if gamma_arm is None else gamma_arm
        gamma = self.PER_ROW * np.asarray(x) + self.PER_COLUMN * np.asarray(y)
        return gamma + arm - np.asarray(alpha_i), np.zeros_like(gamma)


def _calibrated_detector(rot3=0.0):
    """A real, calibrated area detector one meter from the sample."""
    det = DetectorCalibration.Detector2D_SXRD()
    det.detector = pyFAI.detectors.Detector(
        pixel1=PIXEL, pixel2=PIXEL, max_shape=SHAPE
    )
    det.poni1 = SHAPE[0] * PIXEL / 2.0
    det.poni2 = SHAPE[1] * PIXEL / 2.0
    det.rot1 = det.rot2 = 0.0
    det.rot3 = rot3
    det.dist = DIST
    det.set_energy(17.7)
    det.setAzimuthalReference(np.deg2rad(90.0))
    det.setPolarization(0.0, 1.0)
    det.reset()
    det._cached_array = {}
    return det


def test_acceptance_is_measured_edge_to_edge_not_centre_to_centre():
    """A region of ``n`` rows accepts ``n`` pixels' worth, not ``n - 1``.

    The aperture is the outer pixel boundaries. Measuring between the centres
    of the first and last row would understate every rocking acceptance by
    one pixel, which is a 1 % error on a 100-row region and 10 % on a 10-row
    one.
    """
    det = _LinearDetector()

    for rows in (1.0, 4.0, 64.0):
        got = acceptance.out_of_plane_acceptance(det, 300.0, 240.0, rows, 0.0)
        np.testing.assert_allclose(got, rows * det.PER_ROW, rtol=1e-12)


def test_acceptance_is_evaluated_at_the_centre_column():
    """Where the region sits in-plane must not change its rod acceptance.

    The rod crosses the aperture at the region's centre column, so the column
    the acceptance is read at is the centre one. On a detector whose gamma
    also varies along a row, reading it anywhere else would leak the in-plane
    position into an out-of-plane quantity.
    """
    det = _LinearDetector()

    at_left = acceptance.out_of_plane_acceptance(det, 300.0, 10.0, 40.0, 0.0)
    at_right = acceptance.out_of_plane_acceptance(det, 300.0, 470.0, 40.0, 0.0)

    np.testing.assert_allclose(at_left, at_right, rtol=1e-12)
    np.testing.assert_allclose(at_left, 40.0 * det.PER_ROW, rtol=1e-12)


def test_acceptance_is_positive_and_vectorized_over_a_scan():
    """orGUI resizes regions along a scan, so this is normally an array.

    Both the centre and the size vary per frame, and the result must follow
    both without a Python loop and without picking up a sign from the
    direction the rows run in.
    """
    det = _LinearDetector()
    row = np.array([100.0, 300.0, 500.0])
    size = np.array([20.0, 45.0, 80.0])

    got = acceptance.out_of_plane_acceptance(det, row, 240.0, size, 0.0)

    assert got.shape == (3,)
    assert np.all(got > 0)
    np.testing.assert_allclose(got, size * det.PER_ROW, rtol=1e-12)


def test_acceptance_pairs_each_region_with_its_actual_arm_position():
    """Vector arm metadata is pairwise, not an arm-by-region product."""
    det = _calibrated_detector(rot3=np.deg2rad(30.0))
    row = np.array([200.0, 350.0, 500.0])
    column = np.array([100.0, 240.0, 440.0])
    height = np.array([20.0, 40.0, 60.0])
    alpha = np.deg2rad([1.0, 5.0, 10.0])
    gamma_arm = np.deg2rad([0.0, 15.0, 30.0])
    delta_arm = np.deg2rad([0.0, 20.0, 40.0])

    got = acceptance.out_of_plane_acceptance(
        det, row, column, height, alpha, gamma_arm, delta_arm
    )
    expected = np.array(
        [
            acceptance.out_of_plane_acceptance(
                det,
                row[index],
                column[index],
                height[index],
                alpha[index],
                gamma_arm[index],
                delta_arm[index],
            )
            for index in range(row.size)
        ]
    ).reshape(-1)

    assert got.shape == row.shape
    np.testing.assert_allclose(got, expected, rtol=1e-12)


def test_a_zero_or_negative_region_is_rejected():
    """An empty region accepts no rod, and silently returning 0 would divide."""
    det = _LinearDetector()

    with pytest.raises(ValueError, match="positive height in pixels"):
        acceptance.out_of_plane_acceptance(det, 300.0, 240.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="positive height in pixels"):
        acceptance.out_of_plane_acceptance(det, 300.0, 240.0, -5.0, 0.0)


def test_pixel_acceptance_is_the_one_row_case():
    """The resolution limit, and the unit a region's acceptance is read in."""
    det = _LinearDetector()

    one = acceptance.pixel_acceptance(det, 300.0, 240.0, 0.0)
    many = acceptance.out_of_plane_acceptance(det, 300.0, 240.0, 32.0, 0.0)

    np.testing.assert_allclose(one, det.PER_ROW, rtol=1e-12)
    np.testing.assert_allclose(many / one, 32.0, rtol=1e-12)


def test_gamma_range_matches_the_column_span_when_rows_are_iso_gamma():
    """With no roll about the beam, the two estimates agree.

    ``gamma_range`` covers the whole rectangle; ``out_of_plane_acceptance``
    only its centre column. They coincide exactly when lines of constant
    gamma run along the detector rows, which is what makes their ratio a
    usable diagnostic.
    """
    det = _calibrated_detector(rot3=0.0)
    row, column, height, width = 300.0, 240.0, 60.0, 40.0
    alpha = np.deg2rad(0.6)

    column_span = acceptance.out_of_plane_acceptance(
        det, row, column, height, alpha
    )
    full_span = acceptance.gamma_range(
        det, row, column, height, width, alpha
    )

    np.testing.assert_allclose(full_span, column_span, rtol=1e-4)


def test_gamma_range_exceeds_the_column_span_on_a_rolled_detector():
    """Rolling the detector tilts the iso-gamma lines off the rows.

    Then the corners reach further in gamma than the centre column does, and
    the centre-column estimate is the one to trust for a rod-centred region
    while the ratio warns that the geometry is no longer aligned.
    """
    det = _calibrated_detector(rot3=np.deg2rad(20.0))
    row, column, height, width = 300.0, 240.0, 60.0, 120.0
    alpha = np.deg2rad(0.6)

    column_span = acceptance.out_of_plane_acceptance(
        det, row, column, height, alpha
    )
    full_span = acceptance.gamma_range(
        det, row, column, height, width, alpha
    )

    assert full_span > column_span * 1.2


def test_acceptance_scales_with_region_height_on_a_real_detector():
    """A calibrated detector one meter away, checked against small-angle optics.

    Near the beam centre a flat detector subtends ``pixel / dist`` per row to
    first order, so a region of ``n`` rows accepts about ``n * pixel / dist``.
    This pins the units -- radian, not degrees or pixels -- against a number
    written down from the geometry rather than from the code.
    """
    det = _calibrated_detector()
    alpha = np.deg2rad(0.6)
    heights = np.array([10.0, 50.0, 200.0])

    got = acceptance.out_of_plane_acceptance(
        det, SHAPE[0] / 2.0, SHAPE[1] / 2.0, heights, alpha
    )

    np.testing.assert_allclose(got, heights * PIXEL / DIST, rtol=2e-3)
    # Linear in the height to the same order.
    np.testing.assert_allclose(got / heights, got[0] / heights[0], rtol=2e-3)


def test_acceptance_shrinks_where_the_detector_is_oblique():
    """The reason Delta_gamma is not constant along a scan.

    Away from the detector normal the same number of rows subtends less
    angle. A region that keeps its pixel size therefore accepts a shorter
    piece of rod as the reflection moves up the detector, and a rocking scan
    that ignores this reports a rod whose shape follows the detector.
    """
    det = _calibrated_detector()
    alpha = np.deg2rad(0.6)
    centre = acceptance.out_of_plane_acceptance(
        det, SHAPE[0] / 2.0, SHAPE[1] / 2.0, 60.0, alpha
    )
    edge = acceptance.out_of_plane_acceptance(det, 40.0, SHAPE[1] / 2.0, 60.0, alpha)

    assert edge < centre
    assert 0.95 < edge / centre < 0.999


def test_a_moving_arm_barely_changes_the_acceptance():
    """Delta_gamma survives a moving arm; the polarization factor does not.

    Driving the gamma arm moves the region to a completely different exit
    angle -- gamma at the region center shifts by the full arm angle -- but
    the *span* the region subtends is unchanged, because that rotation is
    about the very axis gamma is measured around. A delta-arm rotation and an
    off-center column break the invariance only in the fifth digit.

    This is worth pinning because it is the opposite of the polarization
    factor, where evaluating at the calibrated position instead of the real
    one is a 10 % error at a scattering angle of 18 degrees and 33 % at 30
    (see finding F5). An acceptance estimated with the arm left out is fine;
    a polarization is not.
    """
    det = _calibrated_detector()
    alpha = np.deg2rad(0.6)
    row, column = SHAPE[0] / 2.0, SHAPE[1] / 2.0
    arm = np.deg2rad(25.0)

    at_home = det.surfaceAnglesPoint(np.array([row]), np.array([column]), alpha)[0]
    at_arm = det.surfaceAnglesPoint(
        np.array([row]), np.array([column]), alpha, arm, 0.0
    )[0]
    np.testing.assert_allclose(at_arm - at_home, arm, atol=1e-9)

    home = acceptance.out_of_plane_acceptance(det, row, column, 60.0, alpha)
    moved = acceptance.out_of_plane_acceptance(
        det, row, column, 60.0, alpha, gamma_arm=arm, delta_arm=0.0
    )
    np.testing.assert_allclose(moved, home, rtol=1e-9)

    tilted = acceptance.out_of_plane_acceptance(
        det, row, column, 60.0, alpha, gamma_arm=0.0, delta_arm=np.deg2rad(40.0)
    )
    np.testing.assert_allclose(tilted, home, rtol=1e-4)
    assert not np.isclose(tilted, home, rtol=1e-9)


def test_one_arm_angle_alone_is_rejected():
    """The calibration reference is a rotation, not two independent offsets."""
    det = _calibrated_detector()

    with pytest.raises(ValueError, match="both gamma_arm and delta_arm"):
        acceptance.out_of_plane_acceptance(
            det,
            SHAPE[0] / 2.0,
            SHAPE[1] / 2.0,
            60.0,
            np.deg2rad(0.6),
            gamma_arm=np.deg2rad(25.0),
        )


def test_the_solid_angle_correction_and_the_acceptance_carry_one_obliquity():
    """Which acceptance goes with which sum, and why F6 keeps the raw one.

    A rocking scan's omega-integrated counts in pixel row ``j`` go as that
    row's gamma height ``dgamma_j``. A raw region sum -- which is what a
    region-of-interest integration produces once F6 stops applying the
    solid-angle correction -- therefore pairs with ``sum_j dgamma_j``, which
    is what :func:`out_of_plane_acceptance` returns. A solid-angle corrected
    sum would instead pair with ``sum_j dgamma_j / Omega~_j``, larger by the
    region-mean correction that was applied to the counts.

    The two are pinned together because the second is the trap that made F6
    look optional: the correction cancels against the acceptance in a rocking
    scan, so leaving it in *looks* harmless there. It is not, because a
    stationary sum has no acceptance for it to cancel against.

    The third assertion is the other trap: a *constant* nominal ``n * pixel /
    dist`` is neither partner. On a flat detector
    ``dgamma = nominal cos^2(theta)`` while ``Omega~ = cos^3(theta)``, so the
    corrected-sum partner is ``nominal / cos(theta)`` -- above the nominal,
    where the acceptance itself is below it.
    """
    det = _calibrated_detector()
    alpha = np.deg2rad(0.6)
    column, rows = SHAPE[1] // 2, 60
    solid_angle = np.asarray(det.solidAngleArray(SHAPE), dtype=np.float64)

    # Off the beam centre, where the detector is oblique enough to separate
    # the three candidates.
    for row in (480, 560):
        edges = np.arange(row - rows / 2.0, row + rows / 2.0 + 1.0)
        gamma = np.concatenate(
            [acceptance._surface_gamma(det, e, column, alpha).ravel() for e in edges]
        )
        per_row = np.abs(np.diff(gamma))
        omega = solid_angle[row - rows // 2:row + rows // 2, column]

        edge_to_edge = acceptance.out_of_plane_acceptance(
            det, float(row), float(column), float(rows), alpha
        )
        corrected_partner = np.sum(per_row / omega)

        # The raw-sum partner is the edge-to-edge span.
        np.testing.assert_allclose(np.sum(per_row), edge_to_edge, rtol=1e-9)
        # The corrected-sum partner is that, divided by the mean correction.
        np.testing.assert_allclose(
            corrected_partner, edge_to_edge / omega.mean(), rtol=1e-4
        )
        # Neither equals the nominal, and they straddle it.
        nominal = rows * PIXEL / DIST
        assert edge_to_edge < nominal < corrected_partner


def test_the_acceptance_is_what_the_reduction_asks_for():
    """It plugs straight into the rocking angular factor.

    ``angular_factor`` refuses a rocking scan without an acceptance and
    demands radian; this is the function that supplies it, so the two are
    checked together rather than each against its own convention.
    """
    from orgui.datautils.xrayutils.corrections import measurement

    det = _calibrated_detector()
    alpha = np.deg2rad(0.6)
    delta, gamma = np.deg2rad(16.8), np.deg2rad(5.3)

    d_gamma = acceptance.out_of_plane_acceptance(
        det, SHAPE[0] / 2.0, SHAPE[1] / 2.0, 60.0, alpha
    )
    eta = measurement.angular_factor(
        measurement.ROCKING,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=d_gamma,
    )

    components = measurement.mode_components(
        measurement.ROCKING, alpha=alpha, delta=delta, gamma=gamma
    )
    expected = components["C_Lorentz"] * components["C_rod"] * d_gamma
    np.testing.assert_allclose(eta, expected, rtol=1e-12)

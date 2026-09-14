"""Regression tests for the per-pixel detector correction factors.

:mod:`orgui.datautils.xrayutils.corrections.detector` owns the solid-angle and
polarization arrays, the region reduction of the solid-angle part that a
structure factor has to divide back out (finding F6), and the factor that
moves a polarization correction from the calibrated arm position onto the arm
position a frame was actually measured at (finding F5). See
``doc/design/ctr_structure_factor_scale.md``.
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


def test_the_arm_correction_is_one_at_the_calibrated_position():
    """A fixed-arm scan must come out bit-identical.

    The whole point of expressing finding F5 as a ratio is that it changes
    nothing for a detector that does not move: the correction already applied
    is the right one there.
    """
    det = _calibrated_detector()

    got = detector_corrections.polarization_arm_correction(
        det, 300.0, 240.0, 60.0, 60.0, np.deg2rad(0.6), None, None
    )

    assert got == 1.0


def test_the_arm_correction_recovers_the_documented_errors():
    """The size of finding F5, written out from the two evaluations.

    The calibrated-position polarization is too small a correction for a
    moving arm, so the factor is above one and grows with the scattering
    angle. The reference is the ratio of the polarization at the two arm
    positions at the region centre, which is what the region means reduce to
    for a small region.
    """
    det = _calibrated_detector()
    row, column = 310.0, 244.0
    alpha = np.deg2rad(0.6)

    for two_theta, expected_percent in ((10.0, 3.2), (18.0, 10.7), (30.0, 33.6)):
        arm = np.deg2rad(two_theta)
        got = detector_corrections.polarization_arm_correction(
            det, row, column, 4.0, 4.0, alpha, arm, 0.0
        )

        p_home = det.polarizationAtPoints(
            np.array([row]), np.array([column]), alpha
        )[0]
        p_arm = det.polarizationAtPoints(
            np.array([row]), np.array([column]), alpha, arm, 0.0
        )[0]
        np.testing.assert_allclose(got, p_home / p_arm, rtol=1e-3)
        np.testing.assert_allclose(100.0 * (got - 1.0), expected_percent, rtol=5e-2)


def test_the_arm_correction_averages_over_the_region():
    """At a large scattering angle the polarization is not flat over a region.

    Taking the ratio at the region centre instead of between two region means
    is a fraction of a percent off for a large region at
    :math:`2\\theta = 30` degrees, which is why the means are used.
    """
    det = _calibrated_detector()
    row, column, alpha = 310.0, 244.0, np.deg2rad(0.6)
    arm = np.deg2rad(30.0)

    small = detector_corrections.polarization_arm_correction(
        det, row, column, 2.0, 2.0, alpha, arm, 0.0
    )
    large = detector_corrections.polarization_arm_correction(
        det, row, column, 300.0, 300.0, alpha, arm, 0.0
    )

    assert not np.isclose(small, large, rtol=1e-4)
    # Both still describe the same correction, so they stay close.
    np.testing.assert_allclose(large, small, rtol=2e-2)


def test_the_region_mean_polarization_follows_the_arm():
    """The underlying quantity, which the ratio is built from."""
    det = _calibrated_detector()
    args = (det, 310.0, 244.0, 40.0, 40.0, np.deg2rad(0.6))

    home = detector_corrections.roi_mean_inverse_polarization(*args)
    moved = detector_corrections.roi_mean_inverse_polarization(
        *args, gamma_arm=np.deg2rad(30.0), delta_arm=0.0
    )

    assert home < moved, "a larger scattering angle needs a larger correction"
    np.testing.assert_allclose(home, 1.0, rtol=1e-3)


def test_a_non_finite_region_position_is_neutral_not_fatal():
    """Real scans have frames on which the rod never reaches the detector.

    Those arrive with ``inf`` or ``nan`` region centres -- and with no counts,
    since no pixel was valid. Found on FeReO4 scan 39, where it aborted a
    stationary integration with ``OverflowError: cannot convert float
    infinity to integer`` after the frame loop had already finished.
    """
    det = _calibrated_detector()
    row = np.array([300.0, np.inf, np.nan, 400.0])
    column = np.array([240.0, 240.0, np.nan, np.inf])

    got = detector_corrections.roi_mean_inverse_solid_angle(
        det, row, column, 20.0, 20.0
    )

    assert np.all(np.isfinite(got))
    np.testing.assert_allclose(got[[1, 2, 3]], 1.0, rtol=0.0)
    assert got[0] > 1.0, "the finite frame must still be corrected"


def test_a_non_finite_region_size_is_neutral_too():
    """The size can be degenerate on an off-detector frame as well."""
    det = _calibrated_detector()

    got = detector_corrections.roi_mean_inverse_solid_angle(
        det, 300.0, 240.0, np.array([20.0, np.nan]), 20.0
    )

    assert np.all(np.isfinite(got))
    np.testing.assert_allclose(got[1], 1.0, rtol=0.0)


def test_a_finite_but_empty_region_is_still_rejected():
    """Guarding non-finite sizes must not swallow a genuinely bad one."""
    det = _calibrated_detector()

    with pytest.raises(ValueError, match="positive height"):
        detector_corrections.roi_mean_inverse_solid_angle(
            det, 300.0, 240.0, np.array([20.0, 0.0]), 20.0
        )


def test_the_arm_correction_is_neutral_for_a_non_finite_region():
    """Same guard on the F5 factor, which is fed the same coordinates."""
    det = _calibrated_detector()

    got = detector_corrections.polarization_arm_correction(
        det, np.inf, 240.0, 20.0, 20.0, np.deg2rad(0.6), np.deg2rad(30.0), 0.0
    )

    assert got == 1.0


def test_the_batched_arm_correction_matches_the_per_frame_one():
    """The whole point: one region, many frames, same numbers.

    A rocking or reflectivity scan tracks one region across a curve, with
    the incidence angle and the arm different on every frame -- exactly
    what previously forced ``orGUI._polarizationArmFactor`` into a Python
    loop calling :func:`polarization_arm_correction` once per frame, slow
    enough to make a mu scan with thousands of points unusable. This must
    reproduce that loop exactly, just batched.
    """
    det = _calibrated_detector()
    row, column, row_size, column_size = 310.0, 244.0, 40.0, 40.0
    alpha = np.deg2rad(np.linspace(0.3, 3.0, 9))
    gamma_arm = 2.0 * alpha
    delta_arm = np.deg2rad(np.linspace(-1.0, 1.0, 9))

    got = detector_corrections.polarization_arm_correction_frames(
        det, row, column, row_size, column_size, alpha, gamma_arm, delta_arm
    )
    expected = np.array([
        detector_corrections.polarization_arm_correction(
            det, row, column, row_size, column_size,
            float(alpha[i]), float(gamma_arm[i]), float(delta_arm[i]),
        )
        for i in range(alpha.size)
    ])

    np.testing.assert_allclose(got, expected, rtol=1e-10)
    assert got[0] != 1.0, "the arm does move here, or the test proves nothing"


def test_the_batched_arm_correction_is_one_at_the_calibrated_position():
    """A frame whose arm sits at the calibrated reference stays untouched."""
    det = _calibrated_detector()
    alpha = np.deg2rad([0.6, 0.6, 0.6])
    gamma_arm = np.array([0.0, np.deg2rad(10.0), 0.0])
    delta_arm = np.zeros(3)

    got = detector_corrections.polarization_arm_correction_frames(
        det, 300.0, 240.0, 60.0, 60.0, alpha, gamma_arm, delta_arm
    )

    assert got[0] == 1.0
    assert got[2] == 1.0
    assert got[1] != 1.0


def test_the_batched_arm_correction_is_neutral_for_a_non_finite_region():
    """The whole-region guard: a region with no finite position corrects nothing."""
    det = _calibrated_detector()
    alpha = np.deg2rad([0.6, 0.6])
    gamma_arm = np.deg2rad([10.0, 20.0])
    delta_arm = np.zeros(2)

    got = detector_corrections.polarization_arm_correction_frames(
        det, np.inf, 240.0, 20.0, 20.0, alpha, gamma_arm, delta_arm
    )

    np.testing.assert_array_equal(got, [1.0, 1.0])


def test_the_batched_arm_correction_skips_non_finite_frames_only(recwarn):
    """A frame on which the rod never reached the detector must not corrupt
    or warn about the others -- it is skipped, not computed and masked.

    Found by profiling the mu-scan slowdown: a naive fix that computed
    every frame and masked afterwards emitted a ``RuntimeWarning`` for
    every off-detector frame, since NaN/inf angles were fed to trigonometry
    regardless. ``recwarn`` catches that regression directly.
    """
    det = _calibrated_detector()
    alpha = np.deg2rad(np.linspace(0.3, 3.0, 6))
    gamma_arm = 2.0 * alpha
    delta_arm = np.zeros(6)

    bad_alpha = alpha.copy()
    bad_alpha[2] = np.nan
    bad_gamma_arm = gamma_arm.copy()
    bad_gamma_arm[4] = np.inf

    got = detector_corrections.polarization_arm_correction_frames(
        det, 310.0, 244.0, 40.0, 40.0, bad_alpha, bad_gamma_arm, delta_arm
    )
    expected = detector_corrections.polarization_arm_correction_frames(
        det, 310.0, 244.0, 40.0, 40.0, alpha, gamma_arm, delta_arm
    )

    assert got[2] == 1.0
    assert got[4] == 1.0
    other = [0, 1, 3, 5]
    np.testing.assert_allclose(got[other], expected[other], rtol=1e-10)
    assert not [w for w in recwarn.list if issubclass(w.category, RuntimeWarning)]

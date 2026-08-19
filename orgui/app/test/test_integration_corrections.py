"""Regression tests for the stationary-scan correction factors.

:mod:`orgui.app.integration_corrections` assembles the divisors a stationary
scan is corrected with, in the same convention the rocking-scan integration
uses: the footprint and normalization factors are divided out of the stored
intensity, and the Lorentz factor and rod interception are divided out again
to form the structure factor.
"""

import numpy as np
import pytest

from orgui.app import integration_corrections as ic
from orgui.datautils.xrayutils import geometrycorrections as gc
from orgui.datautils.xrayutils.beamprofile import gaussian_profile, top_hat_profile

roi_sum = pytest.importorskip("orgui.app._roi_sum_accel")

#: Sample size along the beam, in meters.
L = 5e-3

ALPHA = np.deg2rad(np.array([0.2, 0.5, 1.0]))
DELTA = np.deg2rad(np.array([20.0, 21.0, 22.0]))
GAMMA = np.deg2rad(np.array([2.0, 3.0, 4.0]))


class _Scan:
    """Minimal stand-in for a loaded scan."""

    auxillary_counters = ("mon", "diode", "twod")

    def __init__(self, size=3, exposure=0.5, mon=None):
        self._size = size
        if exposure is not None:
            self.exposure_time = exposure
        self.mon = np.asarray(mon if mon is not None else [100.0, 200.0, 400.0])
        self.diode = np.asarray([1.0, 2.0, 3.0])
        self.twod = np.zeros((3, 2))

    def __len__(self):
        return self._size


def test_normalization_divides_by_exposure_and_monitors():
    """Exposure and every monitor multiply into one divisor."""
    divisor, applied = ic.normalization_divisor(_Scan(), True, ("mon",), 3)

    np.testing.assert_allclose(divisor, 0.5 * np.array([100.0, 200.0, 400.0]))
    assert applied == ["exposure", "monitor:mon"]


def test_normalization_can_be_switched_off_per_part():
    """Exposure and monitors are independent switches."""
    monitors_only, applied = ic.normalization_divisor(_Scan(), False, ("mon",), 3)
    np.testing.assert_allclose(monitors_only, [100.0, 200.0, 400.0])
    assert applied == ["monitor:mon"]

    none, applied = ic.normalization_divisor(_Scan(), False, (), 3)
    np.testing.assert_allclose(none, np.ones(3))
    assert applied == []


def test_missing_exposure_is_skipped_not_fatal():
    """A backend without an exposure time still integrates.

    This matches the reciprocal-space reconstruction, which records the
    normalization as unavailable rather than failing the job.
    """
    divisor, applied = ic.normalization_divisor(_Scan(exposure=None), True, (), 3)

    np.testing.assert_allclose(divisor, np.ones(3))
    assert applied == []


def test_normalization_rejects_unusable_counters():
    """A missing or zero counter must fail loudly, not scale by infinity."""
    with pytest.raises(ValueError, match="no monitor counter named"):
        ic.normalization_divisor(_Scan(), False, ("absent",), 3)
    with pytest.raises(ValueError, match="must be finite and nonzero"):
        ic.normalization_divisor(_Scan(mon=[1.0, 0.0, 1.0]), False, ("mon",), 3)
    with pytest.raises(ValueError, match="finite and positive"):
        ic.normalization_divisor(_Scan(exposure=0.0), True, (), 3)


def test_scalar_counters_are_broadcast_over_the_scan():
    """One exposure time for the whole scan applies to every image."""
    divisor, _ = ic.normalization_divisor(_Scan(exposure=0.25), True, (), 3)

    np.testing.assert_allclose(divisor, np.full(3, 0.25))


def test_monitor_candidates_skip_multidimensional_counters():
    """Only per-image one-dimensional counters can normalize."""
    assert ic.monitor_counter_candidates(_Scan()) == ["diode", "mon"]


def test_stationary_factors_use_the_stationary_lorentz():
    """The stationary mode must not use the rocking-scan Lorentz factor."""
    factors = ic.stationary_correction_factors(ALPHA, DELTA, GAMMA, use_lorentz=True)

    np.testing.assert_allclose(
        factors["C_Lorentz"], gc.lorentz_stationary(GAMMA), rtol=1e-12
    )
    np.testing.assert_allclose(factors["C_rod"], gc.rod_interception(GAMMA), rtol=1e-12)
    assert not np.allclose(
        factors["C_Lorentz"], gc.lorentz_rocking_scan(DELTA, ALPHA, GAMMA)
    )
    assert factors.applied == ("lorentz",)


def test_footprint_factors_depend_on_the_incidence_angle():
    """The footprint corrections are evaluated at alpha, the mu circle."""
    profile = gaussian_profile(200e-6)
    factors = ic.stationary_correction_factors(
        ALPHA, DELTA, GAMMA, use_footprint=True, beam_profile=profile, sample_size=L
    )

    flux, area = profile.corrections(ALPHA, L)
    np.testing.assert_allclose(factors["C_flux_on_sample"], flux, rtol=1e-12)
    np.testing.assert_allclose(factors["C_illum_area"], area, rtol=1e-12)
    # Grazing incidence intercepts less of the beam than steeper incidence.
    assert factors["C_flux_on_sample"][0] < factors["C_flux_on_sample"][-1]


def test_footprint_without_a_profile_or_size_is_an_error():
    """Enabling the correction without its inputs must not pass silently."""
    with pytest.raises(ValueError, match="needs a beam profile"):
        ic.stationary_correction_factors(ALPHA, DELTA, GAMMA, use_footprint=True)
    with pytest.raises(ValueError, match="positive sample size"):
        ic.stationary_correction_factors(
            ALPHA,
            DELTA,
            GAMMA,
            use_footprint=True,
            beam_profile=gaussian_profile(200e-6),
            sample_size=0.0,
        )


def test_corrections_are_divisors_applied_in_the_documented_order():
    """Intensity is divided by normalization and footprint; F2 by Lorentz."""
    divisor, _ = ic.normalization_divisor(_Scan(), True, ("mon",), 3)
    profile = top_hat_profile(300e-6)
    factors = ic.stationary_correction_factors(
        ALPHA,
        DELTA,
        GAMMA,
        use_lorentz=True,
        use_footprint=True,
        beam_profile=profile,
        sample_size=L,
        normalization=divisor,
    )

    intensity = np.full(3, 1000.0)
    errors = np.full(3, 25.0)
    corrected, corrected_errors = ic.apply_stationary_corrections(
        intensity, errors, factors
    )

    expected = intensity / (
        factors["C_norm"] * factors["C_flux_on_sample"] * factors["C_illum_area"]
    )
    np.testing.assert_allclose(corrected, expected, rtol=1e-12)
    # The Lorentz factor belongs to the structure factor, not the intensity.
    assert not np.allclose(corrected, expected / factors["C_Lorentz"])

    f2, f2_errors = ic.structure_factor(corrected, corrected_errors, factors)
    np.testing.assert_allclose(
        f2, corrected / (factors["C_Lorentz"] * factors["C_rod"]), rtol=1e-12
    )
    # Dividing by a deterministic factor leaves the relative error alone.
    np.testing.assert_allclose(corrected_errors / corrected, errors / intensity)
    np.testing.assert_allclose(f2_errors / f2, errors / intensity)


def test_disabled_corrections_leave_the_intensity_untouched():
    """With nothing enabled the stored intensity is the raw ROI sum."""
    factors = ic.stationary_correction_factors(ALPHA, DELTA, GAMMA)

    intensity = np.array([10.0, 20.0, 30.0])
    corrected, errors = ic.apply_stationary_corrections(intensity, np.ones(3), factors)

    assert factors.applied == ()
    np.testing.assert_array_equal(corrected, intensity)
    np.testing.assert_array_equal(errors, np.ones(3))


def test_factors_broadcast_to_the_shape_of_alpha():
    """A scalar normalization still yields one factor per image."""
    factors = ic.stationary_correction_factors(
        ALPHA, DELTA, GAMMA, use_lorentz=True, normalization=2.0
    )

    for name in ("C_norm", "C_Lorentz", "C_rod"):
        assert factors[name].shape == ALPHA.shape
    np.testing.assert_allclose(factors["C_norm"], np.full(3, 2.0))


def test_roi_mean_correction_is_independent_of_roi_size():
    """The applied correction is the mean, not a sum scaled by the ROI area.

    orGUI multiplied the ROI-summed correction by the nominal ROI area, on
    top of an intensity that already carried that area. Because the
    projected ROI size varies over the detector, the spurious factor scaled
    two measurements of one rod taken on opposite sides of the beam apart by
    the ratio of their ROI sizes -- about 1.4 on the data this pins.
    """
    # A uniform correction of 1.25 over ROIs of very different sizes.
    pixel_count = np.array([100.0, 357.0, 500.0])
    correction_sum = 1.25 * pixel_count

    factor = ic.roi_mean_correction(correction_sum, pixel_count)

    np.testing.assert_allclose(factor, 1.25)
    # The previous expression, for contrast: it grows with the ROI.
    roi_size = np.array([168.0, 399.0, 500.0])
    previous = correction_sum * (roi_size / pixel_count)
    assert previous.max() / previous.min() > 2.0


def test_roi_mean_correction_handles_fully_masked_rois():
    """An ROI with no valid pixel yields 0 rather than a division warning."""
    factor = ic.roi_mean_correction(np.array([0.0, 12.0]), np.array([0.0, 10.0]))

    np.testing.assert_array_equal(factor, np.array([0.0, 1.2]))


def test_roi_mean_correction_recovers_a_varying_correction():
    """A per-pixel correction averages over the pixels that contributed."""
    rng = np.random.default_rng(3)
    correction = rng.uniform(0.5, 1.5, (7, 40))
    valid = rng.random((7, 40)) > 0.2

    factor = ic.roi_mean_correction((correction * valid).sum(axis=1), valid.sum(axis=1))

    for row in range(7):
        assert factor[row] == pytest.approx(correction[row][valid[row]].mean())


def _synthetic_roi_image():
    """A linear background plus a 1600-count peak inside the center ROI."""
    height = width = 60
    rows, columns = np.mgrid[0:height, 0:width]
    image = 10.0 + 0.5 * columns + 0.25 * rows
    image[28:32, 28:32] += 100.0  # 16 px * 100 counts
    center = np.array([[[25, 35], [25, 35]]], dtype=np.int64)
    backgrounds = [
        np.array([[[a, b], [c, d]]], dtype=np.int64)
        for a, b, c, d in (
            (10, 20, 25, 35),
            (40, 50, 25, 35),
            (25, 35, 10, 20),
            (25, 35, 40, 50),
        )
    ]
    return image, center, backgrounds, (rows, columns)


def _integrate(image, center, backgrounds, correction, roi_size=100.0):
    """Reproduce the ROI arithmetic of the integration paths.

    Mirrors ``orGUI.integrateROI``: sum the center and background ROIs, form
    the background-subtracted intensity scaled to the nominal ROI area, then
    multiply by the mean of the correction array over the center ROI.
    """
    counters = np.zeros((1, 4))
    correction_counters = np.zeros((1, 4))
    mask = np.zeros(image.shape, dtype=bool)
    roi_sum.processImage_Carr(
        image, mask, correction, center, *backgrounds, counters, correction_counters
    )
    croi, cpixel, bgroi, bgpixel = counters[0]
    croibg = (croi - (cpixel / bgpixel) * bgroi) * (roi_size / cpixel)
    factor = ic.roi_mean_correction(
        np.array([correction_counters[0][0]]), np.array([correction_counters[0][1]])
    )[0]
    return croibg * factor


def test_roi_integration_recovers_the_peak_counts():
    """Background subtraction and mask scaling return the true peak counts."""
    image, center, backgrounds, _ = _synthetic_roi_image()

    result = _integrate(image, center, backgrounds, np.ones(image.shape))

    assert result == pytest.approx(1600.0)


def test_roi_correction_scales_the_intensity_not_the_roi_area():
    """A uniform correction multiplies the intensity by exactly that factor.

    orGUI previously multiplied by the ROI-summed correction rescaled to the
    ROI area, which inflated this by the 100 px area of the center ROI. That
    is the regression this pins: the ROI area must not appear in the factor.
    """
    image, center, backgrounds, _ = _synthetic_roi_image()

    result = _integrate(image, center, backgrounds, np.full(image.shape, 1.25))

    assert result == pytest.approx(2000.0)
    # The previous expression, for contrast.
    counters = np.zeros((1, 4))
    correction_counters = np.zeros((1, 4))
    roi_sum.processImage_Carr(
        image,
        np.zeros(image.shape, dtype=bool),
        np.full(image.shape, 1.25),
        center,
        *backgrounds,
        counters,
        correction_counters,
    )
    previous = correction_counters[0][0] * (100.0 / correction_counters[0][1])
    assert previous == pytest.approx(125.0)  # 1.25 * the 100 px ROI area


def test_roi_correction_follows_a_smoothly_varying_correction():
    """A solid-angle-like correction is applied as its mean over the ROI."""
    image, center, backgrounds, (rows, columns) = _synthetic_roi_image()
    correction = 1.0 + 0.002 * columns + 0.001 * rows

    result = _integrate(image, center, backgrounds, correction)

    peak = np.where(
        (rows >= 28) & (rows < 32) & (columns >= 28) & (columns < 32), 100.0, 0.0
    )
    assert result == pytest.approx(float((peak * correction).sum()), rel=1e-6)

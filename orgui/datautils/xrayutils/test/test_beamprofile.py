"""Regression tests for the incident-beam footprint corrections.

These tests pin the two related quantities defined in
:mod:`orgui.datautils.xrayutils.beamprofile`: ``C_flux_on_sample`` is the
diagnostic beam/sample overlap, and ``C_illum_area`` is the numerical active
surface-area divisor that already contains that overlap.

The central requirement is that the numerical
:class:`~orgui.datautils.xrayutils.beamprofile.MeasuredBeamProfile`
evaluates the *same* definitions as the closed-form
:class:`~orgui.datautils.xrayutils.beamprofile.GaussianBeamProfile` orGUI
used before it existed: feeding a sampled Gaussian to the numerical path
must reproduce the analytical corrections, so switching an existing
Gaussian analysis to the numerical path changes nothing, and any difference
seen with a real beam is the shape of that beam rather than a change of
convention.
"""

import numpy as np
import pytest
from scipy import stats

from orgui.datautils.xrayutils.beamprofile import (
    DistributionBeamProfile,
    GaussianBeamProfile,
    MeasuredBeamProfile,
    _trapz_impl,
    generalized_normal_profile,
    profile_from_height_scan,
    read_profile_file,
    skew_normal_profile,
    smoothed_top_hat_profile,
    top_hat_profile,
    trapezoid_profile,
    triangular_profile,
    trim_to_illuminated_edge,
)

#: Angles spanning the grazing-incidence regime where the corrections matter.
ALPHAS = np.deg2rad(np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 90.0]))

#: Sample size along the beam, in meters (the dialog default of 5 mm).
L = 5e-3


def _sampled_gaussian(fwhm, half_width=8.0, n=20001):
    """Tabulate a normalized Gaussian of ``fwhm`` over +- ``half_width`` sigma."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    z = np.linspace(-half_width * sigma, half_width * sigma, n)
    return z, np.exp(-0.5 * (z / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)


def test_measured_profile_reproduces_analytical_gaussian():
    """A sampled Gaussian gives the analytical corrections back.

    This is the compatibility guarantee of the numerical path: same
    definitions, only a different way of evaluating the integrals.
    """
    fwhm = 20e-6
    analytical = GaussianBeamProfile(fwhm)
    numerical = MeasuredBeamProfile(*_sampled_gaussian(fwhm))

    flux_ana, area_ana = analytical.corrections(ALPHAS, L)
    flux_num, area_num = numerical.corrections(ALPHAS, L)

    np.testing.assert_allclose(flux_num, flux_ana, rtol=1e-6)
    np.testing.assert_allclose(area_num, area_ana, rtol=1e-6)


def test_measured_profile_recovers_gaussian_width():
    """FWHM, rms width and the three center definitions of a Gaussian."""
    fwhm = 264.7e-6
    profile = MeasuredBeamProfile(*_sampled_gaussian(fwhm))

    assert profile.fwhm == pytest.approx(fwhm, rel=1e-4)
    assert profile.rms_width == pytest.approx(
        fwhm / (2 * np.sqrt(2 * np.log(2))), rel=1e-4
    )
    # A symmetric profile has all three reference points on top of each other.
    assert profile.centroid == pytest.approx(0.0, abs=1e-9)
    assert profile.peak_position == pytest.approx(0.0, abs=1e-9)
    assert profile.median == pytest.approx(0.0, abs=1e-9)


def test_corrections_preserve_the_shape_of_alpha():
    """``alpha`` arrives as a ``(n_s, n_pts)`` array from the integrator."""
    profile = MeasuredBeamProfile(*_sampled_gaussian(20e-6))
    alpha = np.deg2rad(np.linspace(0.1, 2.0, 12).reshape(3, 4))

    flux, area = profile.corrections(alpha, L)

    assert flux.shape == alpha.shape
    assert area.shape == alpha.shape


def test_limiting_cases_of_both_profile_quantities():
    """Each profile quantity tends to 1 in its corresponding limit."""
    profile = MeasuredBeamProfile(*_sampled_gaussian(20e-6))

    # Sample much larger than the beam: it intercepts the whole flux.
    assert profile.flux_on_sample(np.deg2rad(30.0), L) == pytest.approx(1.0, rel=1e-9)
    # Footprint much smaller than the beam: all of it sees the peak intensity.
    assert profile.illuminated_area_fraction(np.deg2rad(1e-4), L) == pytest.approx(
        1.0, rel=1e-6
    )
    # alpha -> 0 is the removable singularity of C_illum_area, not a nan.
    assert profile.illuminated_area_fraction(0.0, L) == pytest.approx(1.0, rel=1e-9)


def test_top_hat_profile_has_closed_form_corrections():
    """A top hat checks the definitions independently of the Gaussian.

    For a beam of uniform density over a width ``w``, a footprint ``h``
    inside the beam intercepts the fraction ``h / w`` of the flux and is
    illuminated over its whole length.
    """
    w = 100e-6
    z = np.array([-w / 2, -w / 2 + 1e-12, w / 2 - 1e-12, w / 2])
    profile = MeasuredBeamProfile(z, np.array([0.0, 1.0, 1.0, 0.0]))

    # h = L sin(alpha) chosen well inside the top hat.
    alpha = np.arcsin(0.5 * w / L)
    assert profile.flux_on_sample(alpha, L) == pytest.approx(0.5, rel=1e-6)
    assert profile.illuminated_area_fraction(alpha, L) == pytest.approx(1.0, rel=1e-6)

    # A footprint twice the beam width takes all the flux and is half lit.
    alpha_wide = np.arcsin(2.0 * w / L)
    assert profile.flux_on_sample(alpha_wide, L) == pytest.approx(1.0, rel=1e-6)
    assert profile.illuminated_area_fraction(alpha_wide, L) == pytest.approx(
        0.5, rel=1e-6
    )


def test_asymmetric_profile_depends_on_the_centering_choice():
    """Centroid, peak and median differ once the profile is asymmetric."""
    z = np.linspace(-200e-6, 200e-6, 401)
    # Two unequal Gaussians: a strong narrow one plus a weak broad satellite.
    intensity = np.exp(-0.5 * ((z - 20e-6) / 30e-6) ** 2) + 0.4 * np.exp(
        -0.5 * ((z + 80e-6) / 60e-6) ** 2
    )

    centroid = MeasuredBeamProfile(z, intensity, center="centroid")
    peak = MeasuredBeamProfile(z, intensity, center="peak")
    median = MeasuredBeamProfile(z, intensity, center="median")

    assert centroid.sample_center < median.sample_center < peak.sample_center
    # At grazing incidence the footprint only samples the beam at the
    # centering point, so C_illum_area is the local intensity ratio, and
    # centering on the peak is the only choice that reaches 1.
    tiny = np.deg2rad(1e-4)
    # The residual below 1 is the curvature of the 1 micron profile sampling
    # across the footprint, not a property of the correction.
    assert peak.illuminated_area_fraction(tiny, L) == pytest.approx(1.0, rel=1e-4)
    assert centroid.illuminated_area_fraction(tiny, L) < 0.999
    # The intercepted flux is largest when the sample sits on the peak.
    assert peak.flux_on_sample(np.deg2rad(0.5), L) > centroid.flux_on_sample(
        np.deg2rad(0.5), L
    )


def test_offset_shifts_the_sample_center():
    """``offset`` displaces the sample relative to the chosen reference."""
    z, p = _sampled_gaussian(100e-6)
    reference = MeasuredBeamProfile(z, p)
    shifted = MeasuredBeamProfile(z, p, offset=50e-6)

    assert shifted.sample_center == pytest.approx(reference.sample_center + 50e-6)
    # Moving the sample off the beam center can only reduce the flux it takes.
    assert shifted.flux_on_sample(np.deg2rad(0.5), L) < reference.flux_on_sample(
        np.deg2rad(0.5), L
    )


def test_profile_from_height_scan_inverts_an_erf_edge():
    """Differentiating a modeled height scan returns the beam profile.

    The transmitted intensity of an edge cutting into a Gaussian beam is a
    complementary error function; its negated derivative must be the
    Gaussian back again, including a constant transmission offset that the
    derivative removes.
    """
    from scipy import special

    fwhm = 200e-6
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    z = np.linspace(-600e-6, 600e-6, 2001)
    transmitted = 1000.0 * 0.5 * special.erfc(z / (np.sqrt(2) * sigma)) + 37.0

    z_out, profile = profile_from_height_scan(z, transmitted)

    assert np.array_equal(z_out, z)
    recovered = MeasuredBeamProfile(z_out, profile)
    assert recovered.fwhm == pytest.approx(fwhm, rel=1e-3)
    assert recovered.centroid == pytest.approx(0.0, abs=1e-8)


def test_height_scan_monitor_normalisation_removes_drift():
    """Dividing by the monitor removes a smooth incident-flux drift."""
    from scipy import special

    fwhm = 100e-6
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    z = np.linspace(-300e-6, 300e-6, 601)
    clean = 1000.0 * 0.5 * special.erfc(z / (np.sqrt(2) * sigma))
    # Incident flux ramps by 50% over the scan; both the transmitted signal
    # and the monitor follow it.
    monitor = 500.0 * (1.0 + 0.5 * (z - z[0]) / (z[-1] - z[0]))
    measured = clean * monitor / float(np.mean(monitor))

    _, reference = profile_from_height_scan(z, clean)
    _, corrected = profile_from_height_scan(z, measured, monitor)
    _, uncorrected = profile_from_height_scan(z, measured)

    np.testing.assert_allclose(
        corrected, reference, rtol=1e-6, atol=1e-9 * np.max(reference)
    )
    # The drift is a real distortion, so the test would pass trivially
    # without the monitor only if it were negligible.
    assert np.max(np.abs(uncorrected - reference)) > 0.05 * np.max(reference)


def test_trim_to_illuminated_edge_cuts_the_second_edge():
    """A scan passing fully through the beam keeps only the cutting edge."""
    z = np.linspace(0.0, 1.0, 101)
    # Positive bump, then a sign change, then the mirrored leaving edge.
    profile = np.where(z < 0.5, np.sin(np.pi * z / 0.5), -1.0)

    z_cut, profile_cut = trim_to_illuminated_edge(z, profile)

    assert np.all(profile_cut > 0)
    assert z_cut[0] >= z[0]
    assert z_cut[-1] < 0.5


def test_trim_rejects_a_wrong_sign_derivative():
    """An all-negative profile means the scan was differentiated wrongly."""
    z = np.linspace(0.0, 1.0, 11)
    with pytest.raises(ValueError, match="wrong sign"):
        trim_to_illuminated_edge(z, -np.ones_like(z))


def test_measured_profile_validates_its_input():
    """Malformed profiles are rejected with an explicit message."""
    z = np.linspace(0.0, 1.0, 5)
    with pytest.raises(ValueError, match="same length"):
        MeasuredBeamProfile(z, np.ones(4))
    with pytest.raises(ValueError, match="at least two points"):
        MeasuredBeamProfile([0.0], [1.0])
    with pytest.raises(ValueError, match="strictly monotonic"):
        MeasuredBeamProfile([0.0, 1.0, 0.5, 2.0], np.ones(4))
    with pytest.raises(ValueError, match="positive integral"):
        MeasuredBeamProfile(z, -np.ones_like(z))
    with pytest.raises(ValueError, match="centroid"):
        MeasuredBeamProfile(z, np.ones_like(z), center="middle")


def test_descending_profile_is_accepted_unchanged():
    """A profile tabulated from high to low z gives the same corrections."""
    z, p = _sampled_gaussian(50e-6)
    ascending = MeasuredBeamProfile(z, p)
    descending = MeasuredBeamProfile(z[::-1], p[::-1])

    np.testing.assert_allclose(
        descending.corrections(ALPHAS, L), ascending.corrections(ALPHAS, L), rtol=1e-12
    )


def test_read_profile_file_round_trip(tmp_path):
    """Reading back a written profile reproduces the corrections."""
    z, p = _sampled_gaussian(120e-6, n=501)
    path = tmp_path / "profile.dat"
    # Written in millimeters, the unit of a diffractometer height scan.
    np.savetxt(path, np.column_stack([z * 1e3, p]), header="z/mm  profile")

    z_read, p_read = read_profile_file(str(path))

    np.testing.assert_allclose(z_read, z, atol=1e-12)
    np.testing.assert_allclose(
        MeasuredBeamProfile(z_read, p_read).corrections(ALPHAS, L),
        MeasuredBeamProfile(z, p).corrections(ALPHAS, L),
        rtol=1e-9,
    )


def test_read_profile_file_differentiates_a_height_scan(tmp_path):
    """A height-scan file is differentiated, monitor-corrected and trimmed."""
    from scipy import special

    fwhm = 200e-6
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    z = np.linspace(-500e-6, 900e-6, 351)
    # Edge cuts in, then the sample leaves the beam again beyond 400 um.
    transmitted = 1000.0 * 0.5 * special.erfc(z / (np.sqrt(2) * sigma))
    leaving = z > 400e-6
    transmitted[leaving] += (
        1000.0 * 0.5 * special.erfc((700e-6 - z[leaving]) / (np.sqrt(2) * sigma))
    )
    monitor = np.full_like(z, 200.0)

    path = tmp_path / "height_scan.dat"
    np.savetxt(path, np.column_stack([z * 1e3, transmitted, monitor]))

    z_cut, profile = read_profile_file(
        str(path), height_scan=True, monitor_col=2, trim=True
    )

    assert z_cut[-1] < 700e-6  # the leaving edge has been cut away
    assert np.all(profile > 0)
    assert MeasuredBeamProfile(z_cut, profile).fwhm == pytest.approx(fwhm, rel=1e-2)

    z_full, profile_full = read_profile_file(
        str(path), height_scan=True, monitor_col=2, trim=False
    )
    assert z_full.size > z_cut.size
    assert np.any(profile_full < 0)  # the untrimmed scan keeps the second edge


def test_distribution_profile_matches_the_closed_form_gaussian():
    """A scipy normal distribution gives the analytical corrections back.

    The closed-form Gaussian is one member of the distribution family, so
    the generic wrapper has to reproduce it to machine precision rather
    than merely closely.
    """
    fwhm = 20e-6
    analytical = GaussianBeamProfile(fwhm)
    generic = DistributionBeamProfile(stats.norm(loc=0.0, scale=analytical.sigma))

    np.testing.assert_allclose(
        generic.corrections(ALPHAS, L), analytical.corrections(ALPHAS, L), rtol=1e-12
    )
    assert generic.fwhm == pytest.approx(fwhm, rel=1e-5)


def test_top_hat_is_the_linear_slit_correction():
    """A uniform beam reproduces Gibaud, Vignaud & Sinha equation (12).

    ``C_flux = L sin(alpha) / width`` while part of the beam misses the
    sample, and the footprint is fully lit over the same range because the
    beam is flat.
    """
    width = 200e-6
    profile = top_hat_profile(width)

    alpha = np.deg2rad(np.array([0.05, 0.2, 0.5, 1.0]))
    expected = L * np.sin(alpha) / width
    assert np.all(expected < 1.0)

    np.testing.assert_allclose(profile.flux_on_sample(alpha, L), expected, rtol=1e-9)
    np.testing.assert_allclose(
        profile.illuminated_area_fraction(alpha, L), 1.0, rtol=1e-9
    )
    # Once the sample is longer than the beam it takes all of the flux.
    assert profile.flux_on_sample(np.deg2rad(10.0), L) == pytest.approx(1.0, rel=1e-12)


def test_top_hat_peak_resolves_to_the_plateau_midpoint():
    """A flat top has no single argmax; the mode must be its middle."""
    profile = top_hat_profile(200e-6, center="peak")

    assert profile.peak_position == pytest.approx(0.0, abs=1e-9)
    assert profile.sample_center == pytest.approx(0.0, abs=1e-9)


def test_trapezoid_and_triangle_widths():
    """The trapezoid is parameterized by its base and flat-top widths."""
    trapezoid = trapezoid_profile(300e-6, 100e-6)
    # Half height of a symmetric trapezoid is the mean of the two widths.
    assert trapezoid.fwhm == pytest.approx(200e-6, rel=1e-3)

    triangle = triangular_profile(300e-6)
    assert triangle.fwhm == pytest.approx(150e-6, rel=1e-3)
    # A trapezoid whose flat top is the full base is a top hat.
    np.testing.assert_allclose(
        trapezoid_profile(300e-6, 300e-6).corrections(ALPHAS, L),
        top_hat_profile(300e-6).corrections(ALPHAS, L),
        rtol=1e-9,
    )


def test_trapezoid_rejects_an_impossible_flat_top():
    """The flat top cannot be wider than the base."""
    with pytest.raises(ValueError, match="flat_width"):
        trapezoid_profile(100e-6, 200e-6)


def test_smoothed_top_hat_cdf_is_the_integral_of_its_pdf():
    """The closed-form CDF of the blurred top hat matches its density."""
    profile = smoothed_top_hat_profile(200e-6, 20e-6)
    dist = profile._dist

    z = np.linspace(-600e-6, 600e-6, 200001)
    density = dist.pdf(z)
    assert _trapz_impl(density, z) == pytest.approx(1.0, rel=1e-9)
    assert dist.cdf(0.0) == pytest.approx(0.5, abs=1e-12)
    # Compare the analytical CDF with a numerical integral of the density.
    numeric = np.concatenate(
        ([0.0], np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(z)))
    )
    np.testing.assert_allclose(dist.cdf(z), numeric, atol=1e-7)
    # Blurring symmetric edges leaves the width at half maximum unchanged.
    assert profile.fwhm == pytest.approx(200e-6, rel=1e-3)


def test_generalized_normal_interpolates_gaussian_and_top_hat():
    """``flatness`` sweeps the shape between a cusp and a flat beam."""
    fwhm = 200e-6
    alpha = np.deg2rad(0.5)

    gaussian = generalized_normal_profile(fwhm, 2.0)
    np.testing.assert_allclose(
        gaussian.corrections(ALPHAS, L),
        GaussianBeamProfile(fwhm).corrections(ALPHAS, L),
        rtol=1e-6,
    )

    peaked = generalized_normal_profile(fwhm, 1.0).flux_on_sample(alpha, L)
    flat = generalized_normal_profile(fwhm, 40.0).flux_on_sample(alpha, L)
    hat = top_hat_profile(fwhm).flux_on_sample(alpha, L)
    assert peaked < gaussian.flux_on_sample(alpha, L) < flat
    assert flat == pytest.approx(hat, rel=0.02)

    # Every flatness is calibrated to the requested half maximum.
    for flatness in (0.5, 1.0, 2.0, 8.0):
        assert generalized_normal_profile(fwhm, flatness).fwhm == pytest.approx(
            fwhm, rel=1e-3
        )


def test_skew_normal_is_calibrated_and_asymmetric():
    """The skew-normal keeps the requested FWHM and shifts its peak."""
    fwhm = 200e-6
    for skew in (-5.0, 0.0, 3.0):
        assert skew_normal_profile(fwhm, skew).fwhm == pytest.approx(fwhm, rel=1e-3)

    symmetric = skew_normal_profile(fwhm, 0.0)
    assert symmetric.peak_position == pytest.approx(symmetric.centroid, abs=1e-9)

    # A positive skew puts the tail toward larger z, so the peak sits below
    # the center of mass.
    skewed = skew_normal_profile(fwhm, 4.0)
    assert skewed.peak_position < skewed.centroid


def test_heavy_tailed_profile_refuses_centroid_centering():
    """A Cauchy beam has no center of mass to align the sample to."""
    with pytest.raises(ValueError, match="not finite"):
        DistributionBeamProfile(stats.cauchy(scale=100e-6))

    profile = DistributionBeamProfile(stats.cauchy(scale=100e-6), center="median")
    assert profile.fwhm == pytest.approx(200e-6, rel=1e-3)
    assert np.isnan(profile.centroid_position)
    # The tails are real flux that never reaches the sample, so the
    # overspill correction does not saturate even at normal incidence.
    assert profile.flux_on_sample(np.deg2rad(89.0), L) < 0.98


def test_distribution_profile_rejects_unusable_input():
    """Non-distributions and unbounded densities are refused."""
    with pytest.raises(TypeError, match="frozen continuous distribution"):
        DistributionBeamProfile(object())
    # An unfrozen family answers with its standard form, which would be a
    # beam one meter wide rather than an error.
    with pytest.raises(TypeError, match="must be frozen"):
        DistributionBeamProfile(stats.norm)
    # A density that diverges has no peak to reference the active area to.
    with pytest.raises(ValueError, match="diverges"):
        DistributionBeamProfile(stats.gamma(0.5, scale=1e-4), center="median")


def test_distribution_corrections_preserve_the_shape_of_alpha():
    """``alpha`` arrives as a ``(n_s, n_pts)`` array from the integrator."""
    profile = top_hat_profile(50e-6)
    alpha = np.deg2rad(np.linspace(0.1, 2.0, 12).reshape(3, 4))

    flux, area = profile.corrections(alpha, L)

    assert flux.shape == alpha.shape
    assert area.shape == alpha.shape


def test_distribution_offset_and_centering():
    """Centering and offset behave as they do for a measured profile."""
    profile = generalized_normal_profile(200e-6, 2.0)
    shifted = generalized_normal_profile(200e-6, 2.0, offset=50e-6)

    assert shifted.sample_center == pytest.approx(profile.sample_center + 50e-6)
    alpha = np.deg2rad(0.5)
    assert shifted.flux_on_sample(alpha, L) < profile.flux_on_sample(alpha, L)


def test_profile_curve_is_normalized_and_centered():
    """Every profile can be sampled for display around the sample center."""
    profiles = [
        GaussianBeamProfile(200e-6),
        top_hat_profile(200e-6),
        trapezoid_profile(300e-6, 100e-6),
        smoothed_top_hat_profile(200e-6, 20e-6),
        skew_normal_profile(200e-6, 3.0),
        MeasuredBeamProfile(*_sampled_gaussian(200e-6)),
    ]
    for profile in profiles:
        z, p = profile.profile_curve()

        assert z.size == p.size
        assert z[0] < 0 < z[-1]
        assert _trapz_impl(p, z) == pytest.approx(1.0, rel=1e-3)
        assert profile.centroid_position == pytest.approx(0.0, abs=1e-9)

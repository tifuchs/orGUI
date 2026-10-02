"""Independent probability/geometry references for flat sample interception."""

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.special import erf

from orgui.datautils.xrayutils.corrections import beamprofile as bp
from orgui.datautils.xrayutils.corrections.sample_interception import (
    SampleShape,
    overlap,
)


def test_gaussian_and_wide_beam_references():
    """A diamond needs the chord integral, not the bounding length."""
    shape = SampleShape("rectangle", (0.01, 0.01))
    vertical = bp.gaussian_profile(160e-6)
    alpha = np.deg2rad(0.36)
    wide = overlap(shape, vertical, bp.top_hat_profile(0.02), alpha, [0, np.pi / 4])
    np.testing.assert_allclose(
        wide.fraction, [0.178090128539, 0.178155252900], atol=1e-11
    )
    narrow = overlap(shape, vertical, bp.top_hat_profile(1e-7), alpha, [0, np.pi / 4])
    assert narrow.fraction[0] == pytest.approx(0.356180257078, rel=1e-9)
    assert narrow.fraction[1] == pytest.approx(0.486812529405, rel=2e-5)


def test_rectangle_product_and_grazing_limit():
    """Both displaced intervals and H at very small incidence remain correct."""
    vertical = bp.gaussian_profile(160e-6, offset=20e-6)
    horizontal = bp.gaussian_profile(4e-3, offset=-0.5e-3)
    alpha = np.deg2rad([0.05, 0.3, 5])
    shape = SampleShape("rectangle", (0.01, 0.006))
    result = overlap(shape, vertical, horizontal, alpha, 0)
    expected = vertical.flux_over_sine(alpha, 0.01) * horizontal.interval_mass(
        -0.003, 0.003
    )
    np.testing.assert_allclose(result.illumination, expected, rtol=1e-9)
    tiny = overlap(shape, vertical, horizontal, 1e-20, 0)
    assert tiny.illumination == pytest.approx(
        0.01 * vertical.density_at(0) * horizontal.interval_mass(-0.003, 0.003),
        rel=1e-9,
    )


@pytest.mark.parametrize("direction", [0.0, np.inf])
def test_near_right_angle_rectangle_retains_product_reference(direction):
    """Roundoff-sized edge sections preserve overlap and its error estimate."""
    shape = SampleShape("rectangle", (0.01, 0.006))
    vertical = bp.gaussian_profile(160e-6)
    horizontal = bp.top_hat_profile(0.02)
    alpha = np.deg2rad(0.36)
    angle = np.nextafter(np.pi / 2, direction)
    result = overlap(shape, vertical, horizontal, alpha, angle)
    # A quarter-turn exchanges length and width; this reference uses beam
    # interval probabilities rather than the polygon section integrator.
    expected = (
        vertical.interval_mass(-0.003 * np.sin(alpha), 0.003 * np.sin(alpha))
        / np.sin(alpha)
        * horizontal.interval_mass(-0.005, 0.005)
    )
    assert result.illumination == pytest.approx(expected, rel=1e-12)
    assert abs(result.illumination - expected) <= result.error
    assert result.error < abs(expected) * 1e-9


def _square_reference(alpha, azimuth, offset, n=100):
    """Independent tensor integration in sample coordinates, not chord sections."""
    points, weights = leggauss(n)
    x, y = np.meshgrid(
        points * 0.005 + offset[0], points * 0.005 + offset[1], indexing="ij"
    )
    u = np.cos(azimuth) * x - np.sin(azimuth) * y
    v = np.sin(azimuth) * x + np.cos(azimuth) * y
    sigma_z = 160e-6 / np.sqrt(8 * np.log(2))
    sigma_h = 0.0012
    density = np.exp(-0.5 * ((u * np.sin(alpha) / sigma_z) ** 2 + (v / sigma_h) ** 2))
    return (
        0.005**2
        * np.sin(alpha)
        / (2 * np.pi * sigma_z * sigma_h)
        * np.einsum("i,j,ij", weights, weights, density)
    )


def test_offset_rotation_and_square_symmetry():
    """Offset belongs to the rotating sample; it is not a fixed lab translation."""
    shape = SampleShape("rectangle", (0.01, 0.01))
    vertical = bp.gaussian_profile(160e-6)
    horizontal = bp.gaussian_profile(0.0012 * np.sqrt(8 * np.log(2)))
    alpha = np.deg2rad(0.36)
    angles = np.deg2rad([0, 13, 45, 77])
    result = overlap(shape, vertical, horizontal, alpha, angles, offset=(0.002, 0.001))
    expected = [_square_reference(alpha, p, (0.002, 0.001)) for p in angles]
    refined = [_square_reference(alpha, p, (0.002, 0.001), 160) for p in angles]
    np.testing.assert_allclose(expected, refined, atol=1e-12)
    np.testing.assert_allclose(result.fraction, expected, rtol=1e-9)
    base = overlap(shape, vertical, horizontal, alpha, angles)
    rotated = overlap(shape, vertical, horizontal, alpha, angles + np.pi / 2)
    np.testing.assert_allclose(base.fraction, rotated.fraction, atol=1e-12)


def test_circle_polygon_containment_and_concavity():
    """Circle symmetry and polygon winding/concavity share no shape approximation."""
    vertical, horizontal = bp.top_hat_profile(0.001), bp.top_hat_profile(0.001)
    circle = overlap(
        SampleShape("circle", (0.01,)), vertical, horizontal, 0.2, [0, 0.3, 1.5]
    )
    np.testing.assert_allclose(circle.fraction, 1, atol=1e-9)
    square = ((-0.005, -0.005), (0.005, -0.005), (0.005, 0.005), (-0.005, 0.005))
    rect = overlap(
        SampleShape("rectangle", (0.01, 0.01)), vertical, horizontal, 0.2, 0.3
    )
    polygon = overlap(
        SampleShape("polygon", vertices=square[::-1]), vertical, horizontal, 0.2, 0.3
    )
    np.testing.assert_allclose(rect.fraction, polygon.fraction, atol=1e-12)
    # L shape under a uniform beam: area/projection exactly known. Some x
    # slices of this shape have different lengths; winding must not matter.
    vertices = (
        (0, 0),
        (0.01, 0),
        (0.01, 0.004),
        (0.004, 0.004),
        (0.004, 0.01),
        (0, 0.01),
    )
    shape = SampleShape("polygon", vertices=vertices)
    result = overlap(shape, bp.top_hat_profile(0.1), bp.top_hat_profile(0.1), 0.2, 0)
    assert result.fraction == pytest.approx(shape.area * np.sin(0.2) / 0.1**2, rel=1e-9)
    # U outline produces two disjoint transverse intervals in some x slices.
    u_shape = SampleShape(
        "polygon",
        vertices=(
            (0, 0),
            (0.01, 0),
            (0.01, 0.003),
            (0.003, 0.003),
            (0.003, 0.007),
            (0.01, 0.007),
            (0.01, 0.01),
            (0, 0.01),
            (0, 0),
        ),
    )
    split = overlap(u_shape, bp.top_hat_profile(0.1), bp.top_hat_profile(0.1), 0.2, 0)
    assert split.fraction == pytest.approx(
        u_shape.area * np.sin(0.2) / 0.1**2, rel=1e-9
    )
    empty = overlap(
        SampleShape("circle", (0.001,)), vertical, horizontal, 0.2, 0, offset=(0.1, 0.1)
    )
    assert empty.fraction == 0
    assert empty.illumination == 0


@pytest.mark.parametrize("alpha", [1e-12, 0.03, 1.2])
def test_circle_gaussian_matches_radial_probability(alpha):
    """Isotropic surface Gaussian integrates to the exact disk probability."""
    radius, sigma = 0.005, 0.002
    fwhm_factor = np.sqrt(8 * np.log(2))
    result = overlap(
        SampleShape("circle", (2 * radius,)),
        bp.gaussian_profile(sigma * np.sin(alpha) * fwhm_factor),
        bp.gaussian_profile(sigma * fwhm_factor),
        alpha, np.deg2rad([0, 17, 45, 89]),
    )
    probability = -np.expm1(-radius**2 / (2 * sigma**2))
    np.testing.assert_allclose(result.fraction, probability, rtol=1e-10)
    np.testing.assert_allclose(result.illumination, probability / np.sin(alpha),
                               rtol=1e-10)


@pytest.mark.parametrize("width", [1e-7, 0.003, 0.01])
def test_circle_top_hat_matches_disk_strip_area(width):
    """Beam edges cutting a circular surface retain the exact stripe area."""
    radius, vertical_width, alpha = 0.005, 0.1, 0.03
    half = width / 2
    strip_area = 2 * (
        half * np.sqrt(max(0, radius**2 - half**2))
        + radius**2 * np.arcsin(half / radius)
    )
    result = overlap(
        SampleShape("circle", (2 * radius,)),
        bp.top_hat_profile(vertical_width), bp.top_hat_profile(width), alpha, 0,
    )
    expected_h = strip_area / (vertical_width * width)
    assert result.illumination == pytest.approx(expected_h, rel=1e-9)
    assert result.fraction == pytest.approx(np.sin(alpha) * expected_h, rel=1e-9)


def test_displaced_circle_matches_independent_polar_surface_integral():
    """Rotation and laboratory alignment retain the exact circular surface."""
    radius, sigma_z, sigma_y = 0.005, 160e-6 / np.sqrt(8 * np.log(2)), 0.0012
    offset = (0.002, -0.001)
    zv, yh = 20e-6, -0.0003
    alpha = np.array([np.deg2rad(0.36), 0.04, 0.04])
    azimuth = np.array([0.17, 0.8, 1.7])

    def _reference(order):
        nodes, weights = leggauss(order)
        radial = (nodes[:, None] + 1) / 2
        theta = np.pi * (nodes[None, :] + 1)
        area_weights = weights[:, None] * weights[None, :] * np.pi / 2
        values = []
        for incidence, angle in zip(alpha, azimuth):
            cx = offset[0] * np.cos(angle) - offset[1] * np.sin(angle)
            cy = offset[0] * np.sin(angle) + offset[1] * np.cos(angle)
            x = cx + radius * radial * np.cos(theta)
            y = cy + radius * radial * np.sin(theta)
            density = np.exp(-0.5 * (
                ((x * np.sin(incidence) + zv + 20e-6) / sigma_z)**2
                + ((y + yh - 0.0005) / sigma_y)**2
            )) / (2 * np.pi * sigma_z * sigma_y)
            values.append(np.sum(area_weights * radius**2 * radial * density))
        return np.asarray(values)

    expected = _reference(80)
    np.testing.assert_allclose(expected, _reference(120), rtol=1e-11)
    result = overlap(
        SampleShape("circle", (2 * radius,)),
        bp.gaussian_profile(160e-6, offset=20e-6),
        bp.gaussian_profile(sigma_y * np.sqrt(8 * np.log(2)), offset=-0.0005),
        alpha, azimuth, offset=offset, orientation=0.4,
        axis_vertical=zv, axis_horizontal=yh,
    )
    np.testing.assert_allclose(result.illumination, expected, rtol=1e-10)
    np.testing.assert_allclose(result.fraction, np.sin(alpha) * expected, rtol=1e-10)


@pytest.mark.parametrize(
    "vertices",
    [
        ((0, 0), (1, 1), (1, 0), (0, 1)),
        ((0, 0), (1, 0), (2, 0)),
        ((0, 0), (1, 0), (np.nan, 1)),
        ((0, 0), (1, 0), (0, 0), (0, 1)),
        ((0, 0), (2, 0), (1, 0), (1, 1), (0, 1)),
    ],
)
def test_invalid_polygons(vertices):
    """Invalid outlines are refused before numerical integration."""
    with pytest.raises(ValueError):
        SampleShape("polygon", vertices=vertices)


def test_measured_support_probability_and_scale():
    """New density obeys zero support, leaving the old private helper unchanged."""
    profile = bp.MeasuredBeamProfile([-0.001, 0, 0.001], [1, 2, 1])
    np.testing.assert_equal(profile.density_at([-0.002, 0.002]), [0, 0])
    assert profile.interval_mass(-np.inf, np.inf) == 1
    bad = bp.MeasuredBeamProfile([-0.001, 0, 0.001], [-1, 4, 1])
    with pytest.raises(ValueError, match="non-negative"):
        bad.density_at(0)
    rescaled = bp.MeasuredBeamProfile([-0.001, 0, 0.001], [3, 6, 3])
    np.testing.assert_allclose(
        profile.density_at(profile.z), rescaled.density_at(rescaled.z)
    )
    gaussian = bp.GaussianBeamProfile(0.001)
    assert gaussian.interval_mass(0.003, 0.004) > 0
    assert gaussian.interval_mass(0, 1e-20) == pytest.approx(
        1e-20 * gaussian.peak_density
    )


def test_profile_equivalence_and_area_scale():
    """Tabulation converges to the analytical model; area uses peak density."""
    z = np.linspace(-0.001, 0.001, 1001)
    gaussian = bp.gaussian_profile(160e-6)
    measured = bp.MeasuredBeamProfile(z, gaussian.density_at(z))
    shape = SampleShape("rectangle", (0.01, 0.006))
    horizontal = bp.gaussian_profile(0.004)
    exact = overlap(shape, gaussian, horizontal, 0.02, 0.3)
    table = overlap(shape, measured, horizontal, 0.02, 0.3)
    assert table.fraction == pytest.approx(exact.fraction, rel=1e-4)
    peak_density = 3e17
    flux = peak_density / (gaussian.peak_density * horizontal.peak_density)
    assert flux * exact.illumination == pytest.approx(
        peak_density * exact.effective_area
    )
    assert gaussian.interval_mass(-0.001, 0.001) == pytest.approx(
        erf(0.001 / (np.sqrt(2) * gaussian.rms_width))
    )


@pytest.mark.parametrize(
    "factory,args",
    [
        (bp.gaussian_profile, (0.001,)),
        (bp.top_hat_profile, (0.001,)),
        (bp.trapezoid_profile, (0.001, 0.0005)),
        (bp.triangular_profile, (0.001,)),
        (bp.smoothed_top_hat_profile, (0.001, 0.0001)),
        (bp.generalized_normal_profile, (0.001, 4)),
        (bp.skew_normal_profile, (0.001, 3)),
    ],
)
def test_analytical_families(factory, args):
    """Every selectable family supplies normalized density and clipped mass."""
    from scipy.integrate import quad

    profile = factory(*args, offset=0.0003)
    points = profile.integration_points
    lower, upper = -0.00025, 0.0005
    cuts = sorted([lower, upper, *points[(points > lower) & (points < upper)]])
    independent = sum(
        quad(profile.density_at, a, b, epsabs=1e-10)[0]
        for a, b in zip(cuts[:-1], cuts[1:])
    )
    assert profile.interval_mass(lower, upper) == pytest.approx(independent, abs=1e-9)
    assert profile.interval_mass(-np.inf, np.inf) == pytest.approx(1)
    result = overlap(
        SampleShape("rectangle", (0.01, 0.01)),
        profile,
        bp.top_hat_profile(0.02),
        0.3,
        0.37,
    )
    assert 0 < result.fraction <= 1


def test_heavy_tails_custom_compatibility_and_cancellation():
    """Undefined centroids stay explicit; old profiles remain instantiable."""
    from scipy.stats import cauchy

    profile = bp.DistributionBeamProfile(cauchy(scale=0.001), center="peak")
    assert np.isnan(profile.centroid_position)
    result = overlap(
        SampleShape("rectangle", (0.01, 0.01)),
        profile,
        bp.top_hat_profile(0.02),
        0.3,
        0,
    )
    expected = (2 / np.pi * np.arctan(0.005 * np.sin(0.3) / 0.001)) * 0.5
    assert result.fraction == pytest.approx(expected)

    class OldProfile(bp.BeamProfile):
        """Historical subclass implementing only the original abstract API."""

        def flux_on_sample(self, alpha, length):
            """Return an old-style diagnostic without new probability methods."""
            return np.ones_like(alpha)

        illuminated_area_fraction = flux_on_sample

        def profile_curve(self, n=512):
            """Return legacy display samples."""
            return np.array([0]), np.array([1])

        @property
        def centroid_position(self):
            """Return the legacy centred origin."""
            return 0

    old = OldProfile()
    assert old.corrections(0.3, 0.01) == (1, 1)
    with pytest.raises(NotImplementedError, match="2D interception"):
        overlap(SampleShape("circle", (0.01,)), old, profile, 0.3, 0)
    with pytest.raises(InterruptedError):
        overlap(
            SampleShape("circle", (0.01,)),
            profile,
            profile,
            0.3,
            [0, 0.1],
            progress=lambda completed, total: completed < 1,
        )

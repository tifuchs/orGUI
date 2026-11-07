"""Independent probability and conservation checks for recentered roughness."""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.stats import poisson

from ..CTRdistributions import PoissonProfile, RecenteredSkellamProfile


@pytest.mark.parametrize("up,down", [(0, 0), (0, .37), (.37, 0),
                                     (.4, .4), (.1, .8), (2.7, .2),
                                     (2000, 1500)])
@pytest.mark.parametrize("offset", [0, -.35, 2.4])
def test_moments_coverage_and_tail_bounds(up, down, offset):
    """The integer PMF conserves coverage and has the predicted variance."""
    profile = RecenteredSkellamProfile(2, up, down, offset)
    low, high = profile.support()
    layers = np.arange(low, high + 1)
    heights = layers + 1
    masses = profile.probability(heights)
    assert np.all(masses >= 0)
    assert_allclose(masses.sum(), 1, atol=2e-12)
    assert_allclose(masses @ heights, offset, atol=profile.moment_tolerance, rtol=0)
    assert_allclose(masses @ (heights - offset)**2, profile.variance,
                    rtol=2e-10, atol=2e-9)
    assert_allclose(profile.correction(layers).sum(), offset,
                    atol=profile.moment_tolerance, rtol=0)
    occupancy = profile.occupancy(layers)
    assert np.all((0 <= occupancy) & (occupancy <= 1))
    assert np.all(np.diff(occupancy) <= 2e-15)
    assert_allclose(profile.surface_occupancy(layers), masses, atol=1e-15)
    probability, moment = profile.omitted_bounds(low + 1, high + 1)
    assert probability <= profile.tail_probability
    assert moment <= profile.moment_tolerance


@pytest.mark.parametrize("rate", [0, .3, 1, 2.8, 30])
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("offset", [0, -.8, 1.7])
def test_exact_poisson_limits(rate, sign, offset):
    """Each one-sided limit equals the existing compensated Poisson model."""
    profile = RecenteredSkellamProfile(
        1, rate if sign > 0 else 0, rate if sign < 0 else 0, offset)
    reference = PoissonProfile(sign * rate, 1, offset - sign * rate)
    n = np.arange(-100, 101)
    assert_allclose(profile.occupancy(n), reference.occupancy(n), atol=5e-15)
    assert_allclose(profile.probability(n), reference.probability(n), atol=5e-15)


def test_independent_poisson_convolution_and_reflection():
    """A direct double-Poisson sum verifies signs and fractional recentering."""
    profile = RecenteredSkellamProfile(1.3, 2.1, .4)
    other = RecenteredSkellamProfile(1.3, .4, 2.1)
    counts = np.arange(60)
    differences = counts[:, None] - counts[None, :]
    weights = (poisson.pmf(counts, 1.3 * 2.1)[:, None]
               * poisson.pmf(counts, 1.3 * .4)[None, :])
    shift = -(1.3 * 2.1 - 1.3 * .4)
    lower = np.floor(shift)
    fraction = shift - lower
    n = np.arange(-30, 31)
    expected = np.array([
        (weights * ((1 - fraction) * (differences + lower == h)
                    + fraction * (differences + lower + 1 == h))).sum()
        for h in n])
    assert_allclose(profile.probability(n), expected, atol=2e-15)
    assert_allclose(profile.probability(n), other.probability(-n), atol=2e-15)


def test_centering_is_continuous_and_time_zero_is_flat():
    """Integer crossings are continuous even though derivatives can differ."""
    n = np.arange(-30, 31)
    left = RecenteredSkellamProfile(1, 2 - 1e-9, 1)
    right = RecenteredSkellamProfile(1, 2 + 1e-9, 1)
    assert_allclose(left.probability(n), right.probability(n), atol=1e-9)
    flat = RecenteredSkellamProfile(0, 4, 7)
    assert_allclose(flat.occupancy(n), n < 0)
    assert_allclose(flat.probability(n), n == 0)
    assert flat.probability(.5) == 0


def test_tighter_support_preserves_inventory_scope():
    """Zero layer coverage does not assert equal elemental inventories."""
    profile = RecenteredSkellamProfile(1, 2.1, .4)
    tighter = RecenteredSkellamProfile(
        1, 2.1, .4, tail_probability=1e-14, moment_tolerance=1e-13)
    low, high = tighter.support()
    n = np.arange(low, high + 1)
    correction = tighter.correction(n)
    assert abs(correction.sum()) < 1e-12
    # Alternating inventories explicitly invalidate an elemental-mass claim.
    assert abs(correction @ np.where(n % 2, 8., 6.)) > .01
    assert low <= profile.support()[0]
    assert high >= profile.support()[1]


@pytest.mark.parametrize("kwargs", [
    {"t": -1}, {"k_up": -1}, {"k_down": -1}, {"offset": np.nan},
    {"k_up": np.inf}, {"tail_probability": 0}, {"tail_probability": 1},
    {"moment_tolerance": 0}, {"t": 1e308, "k_up": 1e308},
])
def test_invalid_inputs(kwargs):
    """Reject unphysical or numerically undefined parameters."""
    with pytest.raises(ValueError):
        RecenteredSkellamProfile(**kwargs)

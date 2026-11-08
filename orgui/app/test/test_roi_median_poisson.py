"""Counting moments of repaired medians under the flat Poisson donor model."""

import itertools
import math

import numpy as np
import pytest

native = pytest.importorskip("orgui.app._roi_sum_cpp")
pytestmark = pytest.mark.skipif(
    not hasattr(native, "_median_poisson_moments"),
    reason="requires freshly rebuilt model-A native extension",
)


@pytest.mark.parametrize("n,m,shared", [(1, 1, 1), (2, 2, 1), (3, 3, 2), (3, 4, 2)])
def test_median_moments_match_direct_poisson_enumeration(n, m, shared):
    """Enumerate counts, including ties and both odd/even sample medians."""
    rate = 0.1
    counts = np.array(list(itertools.product(range(8), repeat=n + m - shared)))
    pmf = np.array([math.exp(-rate) * rate**k / math.factorial(k) for k in range(8)])
    probability = pmf[counts].prod(axis=1)
    left = np.median(counts[:, :n], axis=1)
    right = np.median(
        np.concatenate((counts[:, :shared], counts[:, n:]), axis=1), axis=1
    )
    mean = probability @ left
    variance = probability @ (left * left) - mean**2
    donor_covariance = probability @ (left * counts[:, 0]) - mean * rate
    covariance = probability @ (left * right) - mean * (probability @ right)
    actual = native._median_poisson_moments(n, rate)
    # The omitted Poisson tail above seven counts contributes <3e-10 here.
    np.testing.assert_allclose(actual, (mean, variance, donor_covariance), atol=3e-10)
    assert native._median_poisson_covariance(n, m, shared, rate) == pytest.approx(
        covariance, abs=3e-10
    )


@pytest.mark.parametrize("rate", [0.0, 0.01, 1.0, 100.0, 1.e6])
def test_one_and_two_donor_medians_are_linear_means(rate):
    """The limiting cases have exact linear Poisson counting moments."""
    for n in (1, 2):
        np.testing.assert_equal(
            native._median_poisson_moments(n, rate), (rate, rate / n, rate / n)
        )
    assert native._median_poisson_covariance(2, 2, 1, rate) == rate / 4
    assert native._median_poisson_covariance(6, 12, 0, rate) == 0


@pytest.fixture(scope="module")
def moment_tables():
    return {n: native._MedianPoissonMomentsTable(n) for n in (6, 7, 12)}


@pytest.fixture(scope="module")
def pair_tables():
    return {
        key: native._MedianPoissonTable(*key)
        for key in ((6, 7, 3), (6, 12, 1), (7, 12, 5))
    }


@pytest.mark.parametrize("rate", [0.0001, 0.03, 0.3, 2.3, 10.3, 100.3, 131072.3])
def test_prepared_tables_resolve_low_rates_and_high_count_tail(
    rate, moment_tables, pair_tables
):
    """Interpolation has a small quantified numerical error, without simulation."""
    for n, table in moment_tables.items():
        np.testing.assert_allclose(
            table.moments(rate), native._median_poisson_moments(n, rate),
            rtol=1.e-3, atol=2.e-7,
        )
    for key, table in pair_tables.items():
        assert table.covariance(rate) == pytest.approx(
            native._median_poisson_covariance(*key, rate), rel=1.e-3, abs=2.e-7
        )


@pytest.mark.parametrize("rate", [0.0001, 0.03, 0.3, 2.3, 18.5, 100.3])
def test_common_rate_grid_preserves_joint_covariance_psd(
    rate, moment_tables, pair_tables
):
    """Include original pixels and three partially overlapping repair medians."""
    donors = [np.arange(6), np.arange(3, 10), np.arange(5, 17)]
    matrix = np.zeros((20, 20))
    matrix[:17, :17] = np.eye(17) * rate
    for i, indices in enumerate(donors):
        _, variance, donor_covariance = moment_tables[len(indices)].moments(rate)
        matrix[17 + i, 17 + i] = variance
        matrix[17 + i, indices] = donor_covariance
        matrix[indices, 17 + i] = donor_covariance
        for j in range(i):
            shared = len(np.intersect1d(donors[j], indices))
            key = (len(donors[j]), len(indices), shared)
            value = pair_tables[key].covariance(rate)
            matrix[17 + i, 17 + j] = matrix[17 + j, 17 + i] = value
    assert np.linalg.eigvalsh(matrix)[0] >= -1.e-10 * max(1., rate)


@pytest.mark.parametrize("n,m,shared,rate", [(6, 7, 3, 0.3), (11, 12, 6, 10.)])
def test_known_rate_moments_match_independent_poisson_simulation(n, m, shared, rate):
    """Compare numerical predictions with empirical batch covariance errors."""
    rng = np.random.default_rng(728314)
    batch_moments = []
    for _ in range(20):
        counts = rng.poisson(rate, size=(10000, n + m - shared))
        left = np.median(counts[:, :n], axis=1)
        right = np.median(
            np.concatenate((counts[:, :shared], counts[:, n:]), axis=1), axis=1
        )
        cov = np.cov(np.stack((left, right, counts[:, 0])))
        batch_moments.append((cov[0, 0], cov[0, 1], cov[0, 2]))
    batches = np.asarray(batch_moments)
    _, variance, donor_covariance = native._median_poisson_moments(n, rate)
    predicted = np.array([
        variance,
        native._median_poisson_covariance(n, m, shared, rate),
        donor_covariance,
    ])
    standard_error = batches.std(axis=0, ddof=1) / np.sqrt(len(batches))
    assert np.all(np.abs(batches.mean(axis=0) - predicted) < 5 * standard_error)


def test_unsupported_moment_arguments_do_not_silently_approximate():
    """Unsupported domains fail rather than report falsely exact uncertainty."""
    with pytest.raises((ValueError, RuntimeError), match="donor"):
        native._median_poisson_moments(33, 10.)
    with pytest.raises((ValueError, RuntimeError), match="rate"):
        native._median_poisson_moments(12, -1.)
    with pytest.raises((ValueError, RuntimeError), match="shared"):
        native._median_poisson_covariance(6, 12, 7, 10.)

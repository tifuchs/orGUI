"""Polynomial background predictions and their original-pixel covariance."""

import concurrent.futures

import numpy as np
import pytest

from orgui.app import _roi_sum_accel as roi_sum


def _rois():
    return [np.array([roi], dtype=np.int64) for roi in (
        [[12, 20], [10, 18]], [[5, 12], [10, 18]],
        [[20, 27], [10, 18]], [[12, 20], [3, 10]],
        [[12, 20], [18, 25]],
    )]


def _basis(x, y, order):
    return np.array([x ** (degree - power) * y ** power
                     for degree in range(order + 1)
                     for power in range(degree + 1)]).T


def _pixel_indices(shape, roi):
    (x0, x1), (y0, y1) = roi
    return np.arange(np.prod(shape)).reshape(shape)[y0:y1, x0:x1].ravel()


def _oracle(image, mask, rois, order):
    """Use a dense least-squares influence vector, with duplicate observations."""
    raw = image.ravel()
    center = _pixel_indices(image.shape, rois[0][0])
    center = center[~mask.ravel()[center]]
    samples = np.concatenate([_pixel_indices(image.shape, roi[0])
                              for roi in rois[1:]])
    samples = samples[~mask.ravel()[samples] & ~np.isnan(raw[samples])]
    sy, sx = np.unravel_index(samples, image.shape)
    cy, cx = np.unravel_index(center, image.shape)
    origin = np.mean(rois[0][0], axis=1) - .5
    x = (sx - origin[0]) / 20.
    y = (sy - origin[1]) / 20.
    design = _basis(x, y, order)
    prediction = _basis((cx - origin[0]) / 20.,
                        (cy - origin[1]) / 20., order).sum(axis=0)
    weights = np.zeros(image.size)
    # Solve X.T @ w = prediction for the minimum-norm influence weights.
    per_sample = np.linalg.lstsq(design.T, prediction, rcond=None)[0]
    per_sample *= len(samples) / len(center)
    np.add.at(weights, samples, per_sample)
    signal = np.zeros(image.size)
    signal[center] = 1.
    variance = np.nan_to_num(raw, nan=0.)
    vc = np.sum(signal * signal * variance)
    vb = np.sum(weights * weights * variance)
    covariance = np.sum(signal * weights * variance)
    ratio = len(center) / len(samples)
    return np.array([vc, vb, covariance,
                     vc + ratio * ratio * vb - 2 * ratio * covariance])


@pytest.mark.parametrize("order", range(5))
def test_polynomial_degree_predicts_signal_background_and_uncertainty(order):
    """All total-degree terms through quartic recover a known smooth background."""
    yy, xx = np.mgrid[:30, :34]
    basis = _basis((xx.ravel() - 16) / 20., (yy.ravel() - 14) / 20., order)
    coefficients = np.linspace(.1, 1., basis.shape[1])
    coefficients[0] = 100.
    background = (basis @ coefficients).reshape(xx.shape)
    image = background.copy()
    image[10:18, 12:20] += 70.
    mask = np.zeros_like(image, dtype=bool)
    mask[10, 12] = True
    mask[5, 16] = True
    image[5, 17] = np.nan
    correction = np.ones_like(image)
    correction[mask] = np.nan
    rois = _rois()
    counters, corrections, moments = [np.zeros((1, 4)) for _ in range(3)]
    roi_sum.processImage_polybg_Carr(image, mask, correction, *rois,
                                   counters, corrections, order, moments)
    valid = ~mask[10:18, 12:20]
    expected = np.sum(background[10:18, 12:20][valid])
    assert counters[0, 2] * counters[0, 1] / counters[0, 3] == pytest.approx(
        expected, rel=1e-10)
    np.testing.assert_allclose(moments[0], _oracle(image, mask, rois, order),
                               rtol=1e-9, atol=1e-8)
    patch, stats = roi_sum.interpolate_polybg_croi(image, mask, *rois, order)
    assert stats["order"] == order
    np.testing.assert_allclose(patch, background[10:18, 12:20], rtol=1e-10)
    plan = roi_sum.PolynomialBackgroundPlan(mask, *rois, order)
    planned, planned_corr, planned_moments = [np.zeros((1, 4)) for _ in range(3)]
    plan.process(image, correction, planned, planned_corr, planned_moments)
    np.testing.assert_allclose(planned, counters, rtol=1e-10)
    np.testing.assert_allclose(planned_corr, corrections, rtol=1e-12)
    np.testing.assert_allclose(planned_moments, moments, rtol=1e-9, atol=1e-8)


@pytest.mark.parametrize("order", [0, 1, 2, 3, 4])
def test_fit_covariance_combines_overlapping_samples_before_squaring(order):
    """Repeated strip pixels and direct signal pixels remain correlated."""
    rng = np.random.default_rng(436)
    image = rng.poisson(20., size=(30, 34)).astype(float)
    mask = rng.random(image.shape) < .03
    image[8, 11] = np.nan
    rois = _rois()
    rois[1][0] = [[5, 15], [8, 20]]
    rois[3][0] = [[8, 22], [5, 13]]
    counters, corrections, moments = [np.zeros((1, 4)) for _ in range(3)]
    roi_sum.processImage_polybg_Carr(image, mask, np.ones_like(image), *rois,
                                   counters, corrections, order, moments)
    np.testing.assert_allclose(moments[0], _oracle(image, mask, rois, order),
                               rtol=1e-9, atol=1e-8)
    assert moments[0, 2] != 0.
    plan = roi_sum.PolynomialBackgroundPlan(mask, *rois, order)
    planned = np.zeros_like(moments)
    plan.process(image, np.ones_like(image), counters, corrections, planned)
    np.testing.assert_allclose(planned, moments, rtol=1e-9, atol=1e-8)


def test_quartic_normalization_handles_large_detector_coordinates():
    """Changing the absolute pixel origin does not change a local quartic fit."""
    yy, xx = np.mgrid[:30, :2000]
    x, y = (xx - 1500.) / 20., (yy - 14.) / 20.
    background = 100 + .3*x**4 + .7*x**3*y - .2*x*y**3 + .5*y**4
    image = background.copy()
    image[10:18, 1496:1504] += 50.
    rois = _rois()
    for roi in rois:
        roi[:, 0] += 1484
    mask = np.zeros_like(image, dtype=bool)
    patch, stats = roi_sum.interpolate_polybg_croi(image, mask, *rois, 4)
    assert stats["order"] == 4
    np.testing.assert_allclose(patch, background[10:18, 1496:1504], atol=1e-9)


def test_rank_deficient_fit_downgrades_and_empty_support_is_reported():
    """A one-dimensional sample cannot determine a 2D plane, even at degree4."""
    image = np.full((30, 34), 7.)
    mask = np.ones_like(image, dtype=bool)
    mask[10:18, 12:20] = False
    mask[8, 12:20] = False
    rois = _rois()
    patch, stats = roi_sum.interpolate_polybg_croi(image, mask, *rois, 4)
    assert stats["order"] == 0
    np.testing.assert_allclose(patch, 7.)
    mask[8, 12:20] = True
    patch, stats = roi_sum.interpolate_polybg_croi(image, mask, *rois, 4)
    assert not stats["success"]
    assert stats["order"] == -1
    assert np.all(np.isnan(patch))
    counters, corrections, moments = [np.zeros((1, 4)) for _ in range(3)]
    roi_sum.processImage_polybg_Carr(image, mask, np.ones_like(image), *rois,
                                   counters, corrections, 4, moments)
    np.testing.assert_array_equal(counters[0], [448., 64., 0., 0.])
    np.testing.assert_array_equal(moments[0], [448., 0., 0., 448.])


def test_repaired_prediction_footprint_does_not_add_fit_samples():
    """A repaired signal pixel is predicted using only original background data."""
    yy, xx = np.mgrid[:30, :34]
    background = 20. + .1 * xx + .2 * yy
    image = background.copy()
    image[10:18, 12:20] += 30.
    mask = np.zeros_like(image, dtype=bool)
    mask[12, 15] = True
    repaired = np.zeros_like(mask)
    repaired[12, 15] = True
    counters, corrections = [np.zeros((1, 4)) for _ in range(2)]
    roi_sum.processImage_polybg_Carr(
        image, mask, np.ones_like(image), *_rois(), counters, corrections,
        4, None, repaired)
    assert counters[0, 1] == 64.
    assert counters[0, 3] == 224.
    assert counters[0, 0] == pytest.approx(image[10:18, 12:20].sum())
    assert counters[0, 2] * 64. / 224. == pytest.approx(
        background[10:18, 12:20].sum())
    with pytest.raises(ValueError, match="repair provenance"):
        roi_sum.processImage_polybg_Carr(
            image, mask, np.ones_like(image), *_rois(), counters, corrections,
            4, np.zeros((1, 4)), repaired)


def test_empty_prediction_keeps_frame_specific_background_sample_count():
    """A fully masked center does not bypass source-NaN handling in the plan."""
    image = np.ones((30, 34))
    mask = np.zeros_like(image, dtype=bool)
    mask[10:18, 12:20] = True
    plan = roi_sum.PolynomialBackgroundPlan(mask, *_rois(), 4)
    image[5, 16] = np.nan
    actual = [np.zeros((1, 4)) for _ in range(3)]
    expected = [np.zeros((1, 4)) for _ in range(3)]
    plan.process(image, np.ones_like(image), *actual)
    roi_sum.processImage_polybg_Carr(image, mask, np.ones_like(image), *_rois(),
                                   *expected[:2], 4, expected[2])
    np.testing.assert_allclose(actual, expected)


def test_plan_owns_geometry_and_is_safe_for_concurrent_frame_workers():
    """Caller mutation and frame-specific NaNs do not alter a shared fit plan."""
    rois = _rois()
    mask = np.zeros((30, 34), dtype=bool)
    plan = roi_sum.PolynomialBackgroundPlan(mask, *rois, 4)
    original_mask = mask.copy()
    original_rois = [roi.copy() for roi in rois]
    mask[:] = True
    rois[0][:] = 0
    rng = np.random.default_rng(493)
    frames = rng.poisson(20., size=(8, 30, 34)).astype(float)
    frames[::2, 5, 15] = np.nan

    def integrate(image):
        result = [np.zeros((1, 4)) for _ in range(3)]
        plan.process(image, np.ones_like(image), *result)
        return result

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        actual = list(executor.map(integrate, frames))
    for image, result in zip(frames, actual):
        expected = [np.zeros((1, 4)) for _ in range(3)]
        roi_sum.processImage_polybg_Carr(image, original_mask, np.ones_like(image),
                                       *original_rois, *expected[:2], 4, expected[2])
        np.testing.assert_allclose(result, expected, rtol=1e-9, atol=1e-8)


@pytest.mark.parametrize("order", [-1, 5])
def test_polynomial_degree_outside_supported_range_is_rejected(order):
    """Invalid programmatic settings cannot overrun the polynomial basis."""
    image = np.ones((30, 34))
    mask = np.zeros_like(image, dtype=bool)
    with pytest.raises(ValueError, match="between 0 and 4"):
        roi_sum.interpolate_polybg_croi(image, mask, *_rois(), order)
    with pytest.raises(ValueError, match="between 0 and 4"):
        roi_sum.PolynomialBackgroundPlan(mask, *_rois(), order)

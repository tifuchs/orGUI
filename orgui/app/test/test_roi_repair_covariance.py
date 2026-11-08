"""Dense original-source covariance oracles for sparse ROI repair."""

import numpy as np
import pytest

from orgui.app import _roi_sum_accel as roi_sum
from orgui.app import _roi_sum_cpp as native


def _provenance(mask, image, radius=2):
    remaining = set(zip(*np.nonzero(mask)))
    records = []
    while remaining:
        pending = [min(remaining)]
        remaining.remove(pending[0])
        targets = []
        while pending:
            y, x = pending.pop()
            targets.append((y, x))
            for point in ((y-1, x), (y+1, x), (y, x-1), (y, x+1)):
                if point in remaining:
                    remaining.remove(point)
                    pending.append(point)
        distances = np.min([
            (np.indices(mask.shape)[0]-y)**2 + (np.indices(mask.shape)[1]-x)**2
            for y, x in targets
        ], axis=0)
        donor_mask = (~mask & ~np.isnan(image) & (distances > 0)
                      & (distances <= radius**2))
        donors = np.flatnonzero(donor_mask)
        weights = 1. / distances.ravel()[donors]
        weights /= weights.sum()
        records.append((np.ravel_multi_index(np.array(targets).T, mask.shape),
                        donors, None if len(targets) == 1 else weights))
    return records


def _dense_oracle(image, mask, rois, order):
    """Build joint covariance of original counts and nonlinear repair values."""
    records = _provenance(mask, image)
    size = image.size
    rate = float(image[~mask & ~np.isnan(image)][0])
    covariance = np.zeros((size+len(records), size+len(records)))
    covariance[np.arange(size), np.arange(size)] = rate
    for i, (_, donors, weights) in enumerate(records):
        if weights is None:
            _, variance, cross = native._median_poisson_moments(len(donors), rate)
            covariance[size+i, size+i] = variance
            covariance[size+i, donors] = cross
        else:
            covariance[size+i, size+i] = rate * np.dot(weights, weights)
            covariance[size+i, donors] = rate * weights
        covariance[:size, size+i] = covariance[size+i, :size]
        for j, (_, other, other_weights) in enumerate(records[:i]):
            shared, first, second = np.intersect1d(donors, other, return_indices=True)
            if weights is None and other_weights is None:
                cross = native._median_poisson_covariance(
                    len(donors), len(other), len(shared), rate)
            elif weights is None:
                cross = native._median_poisson_moments(len(donors), rate)[2]
                cross *= other_weights[second].sum()
            elif other_weights is None:
                cross = native._median_poisson_moments(len(other), rate)[2]
                cross *= weights[first].sum()
            else:
                cross = rate * np.dot(weights[first], other_weights[second])
            covariance[size+i, size+j] = covariance[size+j, size+i] = cross

    indices = np.arange(size).reshape(image.shape)

    def pixels(roi):
        (x0, x1), (y0, y1) = roi[0]
        return indices[y0:y1, x0:x1].ravel()

    center = pixels(rois[0])
    background = np.concatenate([pixels(roi) for roi in rois[1:]])
    unmasked_background_pixels = np.count_nonzero(~mask.ravel()[background])
    background = background[~mask.ravel()[background]
                            & ~np.isnan(image.ravel()[background])]
    coefficients = np.zeros((2, size+len(records)))
    coefficients[0, center[~mask.ravel()[center]]] = 1.
    for i, (targets, _, _) in enumerate(records):
        coefficients[0, size+i] = len(np.intersect1d(center, targets))
    if order is None:
        np.add.at(coefficients[1], background, 1.)
    else:
        cy, cx = np.unravel_index(center, image.shape)
        by, bx = np.unravel_index(background, image.shape)

        def basis(x, y):
            return np.array([x**(degree-j) * y**j for degree in range(order+1)
                             for j in range(degree+1)]).T

        weights = np.linalg.lstsq(basis(bx, by).T,
                                  basis(cx, cy).sum(axis=0), rcond=None)[0]
        weights *= len(background) / len(center)
        np.add.at(coefficients[1], background, weights)
    result = coefficients @ covariance @ coefficients.T
    ratio = len(center) / (len(background) if order is not None
                           else unmasked_background_pixels)
    return [result[0, 0], result[1, 1], result[0, 1],
            result[0, 0]+ratio**2*result[1, 1]-2*ratio*result[0, 1]]


@pytest.mark.parametrize("order", [None, 2, 4])
@pytest.mark.parametrize("missing", [False, True])
def test_shared_median_and_component_match_dense_covariance(order, missing):
    """Include shared medians, a weighted component, and overlapping strips."""
    image = np.full((11, 15), 10.)
    mask = np.zeros_like(image, dtype=bool)
    mask[5, [5, 7, 10, 11]] = True
    if missing:
        image[3, 5] = np.nan  # Changes a median's parity and donor graph.
    rois = [np.array([roi], dtype=np.int64) for roi in (
        [[4, 12], [4, 7]], [[0, 5], [1, 10]], [[11, 15], [1, 10]],
        [[2, 13], [0, 5]], [[2, 13], [6, 11]],
    )]
    gaps = np.empty((0, 2), dtype=np.int32)
    plan = roi_sum.PixelRepairPlan(mask, gaps, gaps, 4, 3, 2, 3, roi_covariance=True)
    counts, correction, moments = [np.empty((1, 4)) for _ in range(3)]
    if order is None:
        roi_sum.processImage_repair_Carr(
            image.copy(), mask, np.ones_like(image), *rois, gaps, gaps,
            counts, correction, 4, 3, 2, 3, moments=moments, plan=plan)
    else:
        roi_sum.processImage_polybg_Carr(
            image.copy(), mask, np.ones_like(image), *rois, counts, correction,
            order, moments=moments, repair_plan=plan)
    np.testing.assert_allclose(moments[0], _dense_oracle(image, mask, rois, order),
                               rtol=2.e-8, atol=1.e-8)
    assert counts[0, 0] == 240.
    if order is not None:
        fit_plan = roi_sum.PolynomialBackgroundPlan(
            mask, *rois, order, repair_plan=plan)
        cached_counts, cached_correction, cached_moments = [
            np.empty((1, 4)) for _ in range(3)
        ]
        fit_plan.process(image.copy(), np.ones_like(image), cached_counts,
                         cached_correction, cached_moments)
        np.testing.assert_allclose(cached_counts, counts, rtol=1.e-10)
        np.testing.assert_allclose(cached_moments, moments, rtol=1.e-10)


def test_weighted_repair_preserves_intensity_and_exact_variance():
    """One two-pixel component shares its weighted donors with signal and BG."""
    image = np.arange(81., dtype=float).reshape(9, 9)
    mask = np.zeros((9, 9), dtype=bool)
    mask[4, 4:6] = True
    rois = [np.array([roi], dtype=np.int64) for roi in (
        [[2, 7], [2, 7]], [[0, 3], [0, 9]], [[6, 9], [0, 9]],
        [[0, 9], [0, 3]], [[0, 9], [6, 9]],
    )]
    gaps = np.empty((0, 2), dtype=np.int32)
    outputs = [np.empty((1, 4)) for _ in range(5)]
    roi_sum.processImage_repair_Carr(image.copy(), mask, np.ones_like(image), *rois,
                                    gaps, gaps, outputs[0], outputs[1])
    roi_sum.processImage_repair_Carr(image.copy(), mask, np.ones_like(image), *rois,
                                    gaps, gaps, outputs[2], outputs[3],
                                    moments=outputs[4])
    np.testing.assert_array_equal(outputs[0], outputs[2])
    np.testing.assert_array_equal(outputs[1], outputs[3])
    targets, donors, weights = _provenance(mask, image)[0]
    a = np.zeros(image.size)
    a.reshape(image.shape)[2:7, 2:7] = 1.
    a[targets] = 0.
    a[donors] += 2 * weights
    assert outputs[4][0, 0] == pytest.approx(np.dot(a*a, image.ravel()))


def test_unsupported_median_preserves_counts_and_marks_uncertainty_unavailable():
    """Large selectable radii must not receive an unvalidated approximation."""
    image = np.full((15, 15), 7.)
    mask = np.zeros_like(image, dtype=bool)
    mask[7, 7] = True
    center = np.array([[[7, 8], [7, 8]]], dtype=np.int64)
    empty_roi = np.zeros((1, 2, 2), dtype=np.int64)
    gaps = np.empty((0, 2), dtype=np.int32)
    plan = roi_sum.PixelRepairPlan(mask, gaps, gaps, 4, 3, 4, 6, roi_covariance=True)
    assert plan.configuration()["unsupported_median_components"] == 1
    counts, correction, moments = [np.empty((1, 4)) for _ in range(3)]
    roi_sum.processImage_repair_Carr(
        image, mask, np.ones_like(image), center, *([empty_roi]*4), gaps, gaps,
        counts, correction, 4, 3, 4, 6, moments=moments, plan=plan)
    assert counts[0, 0] == 7.
    assert np.isnan(moments[0, 0])
    assert np.isnan(moments[0, 3])


def test_plan_snapshot_and_failed_repair_do_not_reuse_prediction_footprint():
    """A frame losing enough donors excludes the unrepaired target from a fit."""
    image = np.full((9, 9), 5.)
    mask = np.zeros_like(image, dtype=bool)
    mask[4, 4] = True
    gaps = np.empty((0, 2), dtype=np.int32)
    plan = roi_sum.PixelRepairPlan(mask, gaps, gaps, 4, 3, 1, 4, roi_covariance=True)
    rois = [np.array([roi], dtype=np.int64) for roi in (
        [[3, 6], [3, 6]], [[0, 3], [0, 9]], [[6, 9], [0, 9]],
        [[3, 6], [0, 3]], [[3, 6], [6, 9]],
    )]
    fit = roi_sum.PolynomialBackgroundPlan(mask, *rois, 2, repair_plan=plan)
    outputs = [np.empty((1, 4)) for _ in range(3)]
    fit.process(image.copy(), np.ones_like(image), *outputs)
    assert outputs[0][0, 1] == 9.
    image[4, 3] = np.nan
    fit.process(image.copy(), np.ones_like(image), *outputs)
    assert outputs[0][0, 1] == 8.
    direct = [np.empty((1, 4)) for _ in range(3)]
    roi_sum.processImage_polybg_Carr(
        image.copy(), mask, np.ones_like(image), *rois, direct[0], direct[1], 2,
        moments=direct[2], repair_plan=plan)
    np.testing.assert_allclose(outputs[0], direct[0])
    np.testing.assert_allclose(outputs[2], direct[2])


def test_table_preparation_scope_keeps_unprepared_geometry_correct():
    """ROI-limited setup still has a correct fallback for a new donor type."""
    image = np.full((11, 11), 10.)
    mask = np.zeros_like(image, dtype=bool)
    mask[0, 0] = mask[5, 5] = True
    region = np.zeros_like(mask)
    region[5, 5] = True
    gaps = np.empty((0, 2), dtype=np.int32)
    plan = roi_sum.PixelRepairPlan(mask, gaps, gaps, 4, 3, 2, 3,
                                  roi_covariance=True, covariance_region=region)
    assert plan.configuration()["median_moment_tables"] == 1
    center = np.array([[[0, 1], [0, 1]]], dtype=np.int64)
    empty_roi = np.zeros((1, 2, 2), dtype=np.int64)
    counts, correction, moments = [np.empty((1, 4)) for _ in range(3)]
    roi_sum.processImage_repair_Carr(
        image, mask, np.ones_like(image), center, *([empty_roi]*4), gaps, gaps,
        counts, correction, 4, 3, 2, 3, moments=moments, plan=plan)
    assert counts[0, 0] == 10.
    assert moments[0, 0] == pytest.approx(native._median_poisson_moments(5, 10.)[1])

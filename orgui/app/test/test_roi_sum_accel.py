import numpy as np
import os
import pytest
import subprocess
import sys

roi_sum = pytest.importorskip("orgui.app._roi_sum_accel")


def _roi(x0, x1, y0, y1):
    return np.array([[x0, x1], [y0, y1]], dtype=np.int64)


def _roi_stack(*rois):
    return np.ascontiguousarray(np.stack(rois), dtype=np.int64)


def _nansum_region(data, roi):
    return np.nansum(data[roi[1, 0] : roi[1, 1], roi[0, 0] : roi[0, 1]])


def _count_unmasked(mask, roi):
    return np.sum(~mask[roi[1, 0] : roi[1, 1], roi[0, 0] : roi[0, 1]])


def test_default_import_uses_cpp_without_loading_numba():
    script = (
        "import sys; "
        "import orgui.app._roi_sum_accel as roi_sum; "
        "print(roi_sum.ROI_ACCEL_BACKEND); "
        "print('numba' in sys.modules); "
        "print(hasattr(roi_sum, 'processImage_polybg_Carr'))"
    )
    env = os.environ.copy()
    env.pop("ORGUI_ACCEL_BACKEND", None)
    result = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        env=env,
        text=True,
    )
    assert result.stdout.splitlines() == ["cpp", "False", "True"]


def test_numba_backend_is_explicit_opt_in():
    script = (
        "import sys; "
        "import orgui.app._roi_sum_accel as roi_sum; "
        "print(roi_sum.ROI_ACCEL_BACKEND); "
        "print('numba' in sys.modules); "
        "print(hasattr(roi_sum, 'processImage_polybg_Carr'))"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        env={**os.environ, "ORGUI_ACCEL_BACKEND": "numba"},
        text=True,
    )
    assert result.stdout.splitlines() == ["numba", "True", "True"]


def test_numpy_backend_is_explicit_selectable():
    script = (
        "import sys; "
        "import orgui.app._roi_sum_accel as roi_sum; "
        "print(roi_sum.ROI_ACCEL_BACKEND); "
        "print('numba' in sys.modules); "
        "print(roi_sum.HAS_ACCEL_BACKEND); "
        "print(hasattr(roi_sum, 'processImage_polybg_Carr'))"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        env={**os.environ, "ORGUI_ACCEL_BACKEND": "numpy"},
        text=True,
    )
    assert result.stdout.splitlines() == ["numpy", "False", "False", "True"]


def test_set_accel_backend_switches_process_backend():
    original_backend = roi_sum.ROI_ACCEL_BACKEND
    try:
        roi_sum.set_accel_backend("numpy")
        assert roi_sum.ROI_ACCEL_BACKEND == "numpy"
        assert roi_sum.HAS_ACCEL_BACKEND is False
        roi_sum.set_accel_backend(original_backend)
        assert roi_sum.ROI_ACCEL_BACKEND == original_backend
    finally:
        roi_sum.set_accel_backend(original_backend)


def _expected(data, mask, correction, rois):
    out = np.zeros((rois[0].shape[0], 4), dtype=np.float64)
    corr = np.zeros_like(out)
    for i in range(rois[0].shape[0]):
        out[i, 0] = _nansum_region(data, rois[0][i])
        corr[i, 0] = _nansum_region(correction, rois[0][i])
        out[i, 1] = corr[i, 1] = _count_unmasked(mask, rois[0][i])
        for roi_group in rois[1:]:
            out[i, 2] += _nansum_region(data, roi_group[i])
            corr[i, 2] += _nansum_region(correction, roi_group[i])
            out[i, 3] += _count_unmasked(mask, roi_group[i])
        corr[i, 3] = out[i, 3]
    return out, corr


def test_process_image_carr_matches_masked_nansum_semantics():
    image = np.arange(36, dtype=np.float64).reshape(6, 6)
    image[1, 3] = np.nan
    original_image = image.copy()
    mask = np.zeros_like(image, dtype=bool)
    mask[0, 0] = True
    mask[2, 2] = True
    correction = np.full_like(image, 2.0)
    correction[mask] = np.nan
    rois = [
        _roi_stack(_roi(0, 3, 0, 3), _roi(2, 6, 2, 5)),
        _roi_stack(_roi(0, 1, 0, 3), _roi(1, 2, 2, 5)),
        _roi_stack(_roi(3, 4, 0, 3), _roi(5, 6, 2, 5)),
        _roi_stack(_roi(0, 3, 0, 1), _roi(2, 6, 1, 2)),
        _roi_stack(_roi(0, 3, 3, 4), _roi(2, 6, 5, 6)),
    ]
    expected_image = image.copy()
    expected_image[mask] = np.nan
    expected, expected_corr = _expected(expected_image, mask, correction, rois)
    all_counters = np.zeros((2, 4), dtype=np.float64)
    correction_counters = np.zeros((2, 4), dtype=np.float64)

    roi_sum.processImage_Carr(
        image,
        mask,
        correction,
        *rois,
        all_counters,
        correction_counters,
    )

    np.testing.assert_allclose(all_counters, expected)
    np.testing.assert_allclose(correction_counters, expected_corr)
    np.testing.assert_allclose(image, original_image, equal_nan=True)


def test_process_image_bg_carr_fills_background_counters():
    image = np.arange(25, dtype=np.float64).reshape(5, 5)
    mask = np.zeros_like(image, dtype=bool)
    mask[1, 1] = True
    correction = np.ones_like(image)
    correction[mask] = np.nan
    background = np.arange(100, 125, dtype=np.float64).reshape(5, 5)
    background[mask] = np.nan
    rois = [
        _roi_stack(_roi(0, 3, 0, 3)),
        _roi_stack(_roi(0, 1, 0, 3)),
        _roi_stack(_roi(3, 4, 0, 3)),
        _roi_stack(_roi(0, 3, 0, 1)),
        _roi_stack(_roi(0, 3, 3, 4)),
    ]
    expected_bg, _ = _expected(background, mask, correction, rois)
    all_counters = np.zeros((1, 4), dtype=np.float64)
    correction_counters = np.zeros((1, 4), dtype=np.float64)
    background_counters = np.zeros((1, 4), dtype=np.float64)

    roi_sum.processImage_bg_Carr(
        image,
        background,
        mask,
        correction,
        *rois,
        all_counters,
        correction_counters,
        background_counters,
    )

    np.testing.assert_allclose(background_counters, expected_bg)


def test_process_image_polybg_carr_fits_local_plane_under_center_roi():
    yy, xx = np.mgrid[:12, :14]
    background = 7.0 + 0.5 * xx - 0.25 * yy
    image = background.copy()
    center = _roi(5, 9, 4, 8)
    image[center[1, 0] : center[1, 1], center[0, 0] : center[0, 1]] += 100.0
    mask = np.zeros_like(image, dtype=bool)
    mask[4, 5] = True
    correction = np.ones_like(image)
    correction[mask] = np.nan
    rois = [
        _roi_stack(center),
        _roi_stack(_roi(3, 5, 4, 8)),
        _roi_stack(_roi(9, 11, 4, 8)),
        _roi_stack(_roi(5, 9, 2, 4)),
        _roi_stack(_roi(5, 9, 8, 10)),
    ]
    all_counters = np.zeros((1, 4), dtype=np.float64)
    correction_counters = np.zeros((1, 4), dtype=np.float64)

    roi_sum.processImage_polybg_Carr(
        image,
        mask,
        correction,
        *rois,
        all_counters,
        correction_counters,
        1,
    )

    expected_raw = np.nansum(
        np.where(mask, np.nan, image)[
            center[1, 0] : center[1, 1],
            center[0, 0] : center[0, 1],
        ]
    )
    expected_background = np.nansum(
        np.where(mask, np.nan, background)[
            center[1, 0] : center[1, 1],
            center[0, 0] : center[0, 1],
        ]
    )
    assert all_counters[0, 0] == pytest.approx(expected_raw)
    assert all_counters[0, 1] == 15.0
    assert all_counters[0, 3] == 32.0
    assert all_counters[0, 2] == pytest.approx(expected_background * 32.0 / 15.0)
    scaled_background = all_counters[0, 2] * all_counters[0, 1] / all_counters[0, 3]
    assert scaled_background == pytest.approx(expected_background)
    assert correction_counters[0, 0] == pytest.approx(15.0)
    assert correction_counters[0, 1] == 15.0


def test_process_image_polybg_carr_fits_local_quadratic_under_center_roi():
    yy, xx = np.mgrid[:18, :18]
    background = (
        3.0 + 0.2 * xx - 0.4 * yy + 0.03 * xx**2 + 0.01 * xx * yy - 0.02 * yy**2
    )
    image = background.copy()
    center = _roi(7, 11, 7, 11)
    image[center[1, 0] : center[1, 1], center[0, 0] : center[0, 1]] += 50.0
    mask = np.zeros_like(image, dtype=bool)
    correction = np.ones_like(image)
    rois = [
        _roi_stack(center),
        _roi_stack(_roi(4, 7, 7, 11)),
        _roi_stack(_roi(11, 14, 7, 11)),
        _roi_stack(_roi(7, 11, 4, 7)),
        _roi_stack(_roi(7, 11, 11, 14)),
    ]
    all_counters = np.zeros((1, 4), dtype=np.float64)
    correction_counters = np.zeros((1, 4), dtype=np.float64)

    roi_sum.processImage_polybg_Carr(
        image,
        mask,
        correction,
        *rois,
        all_counters,
        correction_counters,
        2,
    )

    expected_background = np.sum(
        background[
            center[1, 0] : center[1, 1],
            center[0, 0] : center[0, 1],
        ]
    )
    assert all_counters[0, 3] == 48.0
    assert all_counters[0, 2] == pytest.approx(expected_background * 48.0 / 16.0)
    scaled_background = all_counters[0, 2] * all_counters[0, 1] / all_counters[0, 3]
    assert scaled_background == pytest.approx(expected_background)


def test_interpolate_polybg_croi_returns_patch_and_fit_stats():
    yy, xx = np.mgrid[:12, :14]
    background = 4.0 + 0.3 * xx - 0.2 * yy
    image = background.copy()
    center = _roi(5, 9, 4, 8)
    image[center[1, 0] : center[1, 1], center[0, 0] : center[0, 1]] += 80.0
    mask = np.zeros_like(image, dtype=bool)
    rois = [
        _roi_stack(center),
        _roi_stack(_roi(3, 5, 4, 8)),
        _roi_stack(_roi(9, 11, 4, 8)),
        _roi_stack(_roi(5, 9, 2, 4)),
        _roi_stack(_roi(5, 9, 8, 10)),
    ]

    patch, stats = roi_sum.interpolate_polybg_croi(image, mask, *rois, 1)

    expected_patch = background[
        center[1, 0] : center[1, 1],
        center[0, 0] : center[0, 1],
    ]
    np.testing.assert_allclose(patch, expected_patch)
    assert stats["success"] is True
    assert stats["order"] == 1
    assert stats["sample_count"] == 32.0
    assert stats["rmse"] == pytest.approx(0.0, abs=1e-12)


def test_image_accumulation_helpers_match_numpy_operations():
    image = np.array([[1.0, 4.0], [np.nan, 3.0]])
    background = np.array([[0.5, 5.0], [1.0, 1.0]])
    sum_image = np.ones_like(image)
    max_image = np.array([[0.0, 6.0], [2.0, 1.0]])

    roi_sum.calcMaxSum(image, sum_image, max_image)

    np.testing.assert_allclose(sum_image, np.ones_like(image) + image, equal_nan=True)
    np.testing.assert_allclose(
        max_image,
        np.maximum(image, [[0.0, 6.0], [2.0, 1.0]]),
        equal_nan=True,
    )

    roi_sum.calcBgSub(image, background)
    np.testing.assert_allclose(image, [[0.5, -1.0], [np.nan, 2.0]], equal_nan=True)

    sum_bg = np.zeros((2, 2), dtype=np.float64)
    max_bg = np.zeros((2, 2), dtype=np.float64)
    roi_sum.calcMaxSum_bg(image, sum_bg, max_bg, background)
    expected = image - background
    np.testing.assert_allclose(sum_bg, expected, equal_nan=True)
    np.testing.assert_allclose(max_bg, np.maximum(expected, 0.0), equal_nan=True)

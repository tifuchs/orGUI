import numpy as np

from orgui.app import _peak_search_accel


def test_masked_maximum_subtracts_background_and_skips_near_mask():
    image = np.array(
        [
            [1.0, 2.0, 3.0],
            [4.0, 100.0, 6.0],
            [7.0, 8.0, 90.0],
        ]
    )
    background = np.zeros_like(image)
    background[2, 2] = 50.0
    mask = np.zeros_like(image, dtype=bool)
    mask[1, 0] = True

    result = _peak_search_accel.masked_maximum(
        image,
        mask,
        background,
        mask_distance=1,
    )

    assert result.valid is True
    assert result.value == 40.0
    assert result.x == 2.0
    assert result.y == 2.0


def test_masked_maximum_reports_invalid_when_no_finite_pixels():
    image = np.full((2, 2), np.nan)

    result = _peak_search_accel.masked_maximum(image)

    assert result.valid is False


def test_masked_roi_sum_subtracts_background_and_mask():
    image = np.arange(16, dtype=float).reshape(4, 4)
    background = np.ones_like(image)
    mask = np.zeros_like(image, dtype=bool)
    mask[1, 1] = True

    result = _peak_search_accel.masked_roi_sum(
        image,
        mask,
        background,
        1,
        3,
        1,
        3,
    )

    expected = image - background
    expected[mask] = np.nan
    assert result == np.nansum(expected[1:3, 1:3])

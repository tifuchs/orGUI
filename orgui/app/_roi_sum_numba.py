"""Optional Numba ROI summing kernels.

This backend mirrors the non-polynomial C++ ROI summing API. Polynomial
background fitting remains C++ only.
"""

from __future__ import annotations

from numba import njit
import numpy as np


@njit(inline="always")
def _nansum_unmasked(data, mask, roi, n):
    total = 0.0
    x0 = roi[n, 0, 0]
    x1 = roi[n, 0, 1]
    y0 = roi[n, 1, 0]
    y1 = roi[n, 1, 1]
    for y in range(y0, y1):
        for x in range(x0, x1):
            if not mask[y, x]:
                value = data[y, x]
                if not np.isnan(value):
                    total += value
    return total


@njit(inline="always")
def _count_unmasked(mask, roi, n):
    total = 0.0
    x0 = roi[n, 0, 0]
    x1 = roi[n, 0, 1]
    y0 = roi[n, 1, 0]
    y1 = roi[n, 1, 1]
    for y in range(y0, y1):
        for x in range(x0, x1):
            if not mask[y, x]:
                total += 1.0
    return total


@njit(nogil=True, cache=True)
def processImage_Carr(
    image,
    mask,
    correction,
    center,
    left,
    right,
    top,
    bottom,
    all_counters,
    correction_counters,
):
    for i in range(center.shape[0]):
        all_counters[i, 0] = _nansum_unmasked(image, mask, center, i)
        all_counters[i, 2] = (
            _nansum_unmasked(image, mask, left, i)
            + _nansum_unmasked(image, mask, right, i)
            + _nansum_unmasked(image, mask, top, i)
            + _nansum_unmasked(image, mask, bottom, i)
        )
        correction_counters[i, 0] = _nansum_unmasked(
            correction, mask, center, i
        )
        correction_counters[i, 2] = (
            _nansum_unmasked(correction, mask, left, i)
            + _nansum_unmasked(correction, mask, right, i)
            + _nansum_unmasked(correction, mask, top, i)
            + _nansum_unmasked(correction, mask, bottom, i)
        )
        signal_pixels = _count_unmasked(mask, center, i)
        background_pixels = (
            _count_unmasked(mask, left, i)
            + _count_unmasked(mask, right, i)
            + _count_unmasked(mask, top, i)
            + _count_unmasked(mask, bottom, i)
        )
        all_counters[i, 1] = signal_pixels
        all_counters[i, 3] = background_pixels
        correction_counters[i, 1] = signal_pixels
        correction_counters[i, 3] = background_pixels


@njit(nogil=True, cache=True)
def processImage_bg_Carr(
    image,
    background,
    mask,
    correction,
    center,
    left,
    right,
    top,
    bottom,
    all_counters,
    correction_counters,
    background_counters,
):
    processImage_Carr(
        image,
        mask,
        correction,
        center,
        left,
        right,
        top,
        bottom,
        all_counters,
        correction_counters,
    )
    for i in range(center.shape[0]):
        background_counters[i, 0] = _nansum_unmasked(
            background, mask, center, i
        )
        background_counters[i, 2] = (
            _nansum_unmasked(background, mask, left, i)
            + _nansum_unmasked(background, mask, right, i)
            + _nansum_unmasked(background, mask, top, i)
            + _nansum_unmasked(background, mask, bottom, i)
        )
        background_counters[i, 1] = all_counters[i, 1]
        background_counters[i, 3] = all_counters[i, 3]


@njit(nogil=True, cache=True)
def calcMaxSum(image, sum_image, max_image):
    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            value = image[y, x]
            sum_image[y, x] += value
            if np.isnan(value) or value > max_image[y, x]:
                max_image[y, x] = value


@njit(nogil=True, cache=True)
def calcBgSub(image, background):
    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            image[y, x] -= background[y, x]


@njit(nogil=True, cache=True)
def calcMaxSum_bg(image, sum_image, max_image, background):
    for y in range(image.shape[0]):
        for x in range(image.shape[1]):
            value = image[y, x] - background[y, x]
            sum_image[y, x] += value
            if np.isnan(value) or value > max_image[y, x]:
                max_image[y, x] = value

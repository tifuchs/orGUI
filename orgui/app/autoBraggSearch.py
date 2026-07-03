# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
"""Automatic Bragg peak search and UB seed ranking."""

from __future__ import annotations

import copy
from dataclasses import dataclass
import logging

import numpy as np

from ..datautils.xrayutils import HKLVlieg, ReciprocalNavigation as rn
from . import imagePeakFinder

logger = logging.getLogger(__name__)


@dataclass
class ImageMaximum:
    """Highest usable detector pixel in one scan image.

    :param int imageno:
        Image index in the active scan.
    :param numpy.ndarray xy:
        Detector coordinate ``(x, y)`` in pixels.
    :param float value:
        Pixel intensity in detector counts or backend-native intensity units.
    :param float sharpness:
        Robust increase score relative to preceding image maxima.
    :param float derivative_sharpness:
        Robust score of the current max-trace step relative to preceding
        first differences.
    :param float prominence:
        Candidate prominence relative to the rolling baseline, in
        backend-native intensity units.
    """

    imageno: int
    xy: np.ndarray
    value: float
    sharpness: float = np.nan
    derivative_sharpness: float = np.nan
    prominence: float = np.nan


@dataclass
class RefinedPeak:
    """Refined detector and scan-axis position of a candidate peak.

    :param ImageMaximum candidate:
        Initial single-image maximum.
    :param numpy.ndarray xy:
        Refined detector coordinate ``(x, y)`` in pixels.
    :param int imageno:
        Image index closest to the refined rocking-curve center.
    :param float axis_value:
        Scan-axis value in the scan's stored unit, usually deg.
    :param float value:
        Initial maximum pixel intensity.
    """

    candidate: ImageMaximum
    xy: np.ndarray
    imageno: int
    axis_value: float
    value: float
    rocking_curve: np.ndarray | None = None
    rocking_axis: np.ndarray | None = None


@dataclass
class UBSeed:
    """Candidate UB seed derived from one indexed Bragg peak.

    :param numpy.ndarray U:
        Orientation matrix.
    :param numpy.ndarray hkl:
        Candidate Miller indices in r.l.u.
    :param RefinedPeak peak:
        Refined peak used to construct the seed.
    :param float score:
        Lower is better. The score is the median reciprocal-space indexing
        error across available observed peak candidates.
    :param float norm_mismatch:
        Absolute Q-norm mismatch in Angstrom^-1 for the seed peak.
    :param float position_mismatch:
        Detector-position mismatch in pixels for the seed peak, if this seed
        came from predicted detector-position matching.
    """

    U: np.ndarray
    hkl: np.ndarray
    peak: RefinedPeak
    score: float
    norm_mismatch: float
    position_mismatch: float = np.nan


def _masked_maximum(img, mask=None):
    """Return max value and ``(x, y)`` for an image while ignoring masked pixels.

    :param numpy.ndarray img:
        Detector image.
    :param numpy.ndarray mask:
        Optional boolean mask where ``True`` pixels are invalid.
    :returns:
        Tuple ``(value, xy)``.
    :rtype: tuple[float, numpy.ndarray]
    """
    arr = np.asarray(img)
    if mask is not None:
        arr = np.where(mask, -np.inf, arr)
    flat = int(np.nanargmax(arr))
    y, x = np.unravel_index(flat, arr.shape)
    return float(arr[y, x]), np.array([x, y], dtype=float)


def far_from_mask(mask, xy, distance):
    """Return whether ``xy`` is at least ``distance`` pixels from masked pixels.

    :param numpy.ndarray mask:
        Boolean mask where ``True`` pixels are invalid.
    :param numpy.ndarray xy:
        Detector coordinate ``(x, y)`` in pixels.
    :param int distance:
        Rejection distance in pixels.
    :rtype: bool
    """
    if mask is None or distance <= 0:
        return True
    x, y = np.rint(xy).astype(int)
    if y < 0 or x < 0 or y >= mask.shape[0] or x >= mask.shape[1]:
        return False
    y0 = max(0, y - distance)
    y1 = min(mask.shape[0], y + distance + 1)
    x0 = max(0, x - distance)
    x1 = min(mask.shape[1], x + distance + 1)
    return not np.any(mask[y0:y1, x0:x1])


def _robust_location_scale(values):
    """Return median and robust standard-deviation estimate for ``values``."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    return median, 1.4826 * mad


def _robust_z(value, values):
    """Return robust z-score for ``value`` against previous ``values``."""
    median, scale = _robust_location_scale(values)
    if not np.isfinite(median):
        return np.nan
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = np.std(values)
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = np.finfo(float).eps
    return float((value - median) / scale)


def iter_sharp_peak_candidates(
    fscan,
    mask=None,
    excluded_images=(),
    burn_in=5,
    history=20,
    min_history=5,
    level_z=6.0,
    derivative_z=6.0,
    lookahead=2,
    min_prominence_z=4.0,
    refractory=3,
    mask_distance=3,
):
    """Yield sharp max-trace peak candidates while reading images sequentially.

    This is an online detector for the 1-D trace
    ``max(image_i)``. It never needs the full scan. Each image is read at most
    once by this detector, and candidates are emitted as soon as a rolling
    robust level/derivative test plus a short lookahead confirmation passes.
    If the caller rejects a yielded candidate, iteration can continue from the
    next unread image.

    :param orgui.backend.scans.Scan fscan:
        Active scan object.
    :param numpy.ndarray mask:
        Optional boolean mask where ``True`` pixels are invalid.
    :param excluded_images:
        Iterable of image indices to skip.
    :param int burn_in:
        Number of initial images ignored for candidate triggering.
    :param int history:
        Number of preceding max values used for rolling baseline estimation.
    :param int min_history:
        Minimum finite baseline points required before triggering.
    :param float level_z:
        Minimum robust z-score of the current max value.
    :param float derivative_z:
        Minimum robust z-score of the current first difference.
    :param int lookahead:
        Number of following images read to confirm a local maximum.
    :param float min_prominence_z:
        Minimum candidate prominence relative to rolling baseline scale.
    :param int refractory:
        Number of images skipped after yielding a candidate to avoid repeated
        candidates from one broad peak.
    :param int mask_distance:
        Reject maxima within this many pixels of a masked pixel.
    :yields:
        :class:`ImageMaximum` candidates ordered by acquisition.
    """
    excluded = {int(i) for i in np.atleast_1d(excluded_images)}
    values = {}
    maxima = {}
    next_allowed = int(burn_in)
    imageno = 0

    def read_maximum(index):
        if index in values:
            return maxima.get(index)
        if index in excluded:
            values[index] = np.nan
            return None
        image = fscan.get_raw_img(index).img
        try:
            value, xy = _masked_maximum(image, mask)
        except ValueError:
            logger.warning("Skipping image %s without finite pixels.", index)
            values[index] = np.nan
            return None
        if not np.isfinite(value) or not far_from_mask(mask, xy, mask_distance):
            values[index] = np.nan
            return None
        maximum = ImageMaximum(index, xy, value)
        values[index] = value
        maxima[index] = maximum
        return maximum

    while imageno < len(fscan):
        maximum = read_maximum(imageno)
        if maximum is None:
            imageno += 1
            continue

        history_start = max(0, imageno - history)
        baseline = np.array(
            [values.get(i, np.nan) for i in range(history_start, imageno)]
        )
        baseline = baseline[np.isfinite(baseline)]
        if imageno < next_allowed or baseline.size < min_history:
            imageno += 1
            continue

        diffs = np.diff(baseline)
        previous = [
            values.get(i, np.nan)
            for i in range(history_start, imageno)
            if np.isfinite(values.get(i, np.nan))
        ]
        current_diff = maximum.value - previous[-1]
        level_score = _robust_z(maximum.value, baseline)
        diff_score = _robust_z(current_diff, diffs)
        maximum.sharpness = level_score
        maximum.derivative_sharpness = diff_score

        if (
            not np.isfinite(level_score)
            or not np.isfinite(diff_score)
            or level_score < level_z
            or diff_score < derivative_z
        ):
            imageno += 1
            continue

        window = [maximum]
        for look_idx in range(imageno + 1, min(len(fscan), imageno + lookahead + 1)):
            look_maximum = read_maximum(look_idx)
            if look_maximum is not None:
                window.append(look_maximum)

        if not window:
            imageno += 1
            continue
        candidate = max(window, key=lambda item: item.value)
        median, scale = _robust_location_scale(baseline)
        if not np.isfinite(scale) or scale <= np.finfo(float).eps:
            scale = np.std(baseline)
        if not np.isfinite(scale) or scale <= np.finfo(float).eps:
            scale = np.finfo(float).eps
        prominence = candidate.value - median
        candidate.prominence = float(prominence)
        if prominence / scale >= min_prominence_z:
            candidate.sharpness = _robust_z(candidate.value, baseline)
            prev_values = [
                values.get(i, np.nan)
                for i in range(history_start, candidate.imageno)
                if np.isfinite(values.get(i, np.nan))
            ]
            prev_value = prev_values[-1]
            candidate.derivative_sharpness = _robust_z(
                candidate.value - prev_value, diffs
            )
            yield candidate
            next_allowed = candidate.imageno + refractory + 1

        imageno = max(imageno + 1, candidate.imageno + 1)


def scan_image_maxima(
    fscan,
    mask=None,
    excluded_images=(),
    burn_in=5,
    history=20,
    mask_distance=3,
):
    """Read scan images once and calculate the max-pixel trace.

    The returned candidates are sorted by robust sharp-increase score. The
    image data are accessed exactly once in this pass.

    :param orgui.backend.scans.Scan fscan:
        Active scan object.
    :param numpy.ndarray mask:
        Optional boolean mask where ``True`` pixels are invalid.
    :param excluded_images:
        Iterable of image indices to skip.
    :param int burn_in:
        Number of initial images excluded from candidate selection.
    :param int history:
        Number of preceding max values used for robust baseline estimation.
    :param int mask_distance:
        Reject maxima within this many pixels of a masked pixel.
    :returns:
        Sorted maxima.
    :rtype: list[ImageMaximum]
    """
    excluded = {int(i) for i in np.atleast_1d(excluded_images)}
    maxima = []
    values = np.full(len(fscan), np.nan, dtype=float)
    for imageno in range(len(fscan)):
        if imageno in excluded:
            continue
        image = fscan.get_raw_img(imageno).img
        try:
            value, xy = _masked_maximum(image, mask)
        except ValueError:
            logger.warning("Skipping image %s without finite pixels.", imageno)
            continue
        if not np.isfinite(value) or not far_from_mask(mask, xy, mask_distance):
            continue
        values[imageno] = value
        maxima.append(ImageMaximum(imageno, xy, value))

    for maximum in maxima:
        if maximum.imageno < burn_in:
            continue
        start = max(0, maximum.imageno - history)
        baseline_values = values[start : maximum.imageno]
        baseline_values = baseline_values[np.isfinite(baseline_values)]
        if baseline_values.size < max(3, min(history, 5)):
            continue
        median = np.median(baseline_values)
        mad = np.median(np.abs(baseline_values - median))
        maximum.sharpness = (maximum.value - median) / (
            1.4826 * mad + np.finfo(float).eps
        )

    return sorted(
        maxima,
        key=lambda item: (
            np.nan_to_num(item.sharpness, nan=-np.inf),
            item.value,
        ),
        reverse=True,
    )


def refine_peak_3d(
    fscan,
    candidate,
    axis_half_width=1.0,
    roi_size=(80, 80),
    fine_axis_half_width=0.4,
    fine_roi_size=(40, 40),
    mask=None,
    excluded_images=(),
    max_workers=1,
):
    """Refine one Bragg-peak candidate by detector COM and rocking COM.

    :param orgui.backend.scans.Scan fscan:
        Active scan object.
    :param ImageMaximum candidate:
        Initial maximum candidate.
    :param float axis_half_width:
        First-pass scan-axis half-width around the candidate image, in the
        scan's stored unit, usually deg.
    :param tuple roi_size:
        First-pass ``(vertical, horizontal)`` ROI size in pixels.
    :param float fine_axis_half_width:
        Second-pass scan-axis half-width around the first-pass rocking COM, in
        the scan's stored unit, usually deg.
    :param tuple fine_roi_size:
        Second-pass ``(vertical, horizontal)`` ROI size in pixels.
    :param numpy.ndarray mask:
        Optional boolean mask where ``True`` pixels are invalid.
    :param excluded_images:
        Iterable of image indices to skip.
    :param int max_workers:
        Worker count used by the existing image-range peak search.
    :returns:
        Refined peak.
    :rtype: RefinedPeak
    """
    kwargs = {
        "excluded_images": excluded_images,
        "max_workers": max_workers,
    }
    if mask is not None:
        kwargs["mask"] = mask
    axis0 = float(fscan.axis[candidate.imageno])
    result = imagePeakFinder.find_COM_Image(
        candidate.xy,
        roi_size[1],
        roi_size[0],
        fscan,
        [axis0 - axis_half_width, axis0 + axis_half_width],
        **kwargs,
    )
    result = imagePeakFinder.find_COM_Image(
        result["xy"],
        fine_roi_size[1],
        fine_roi_size[0],
        fscan,
        [
            result["axis_com"] - fine_axis_half_width,
            result["axis_com"] + fine_axis_half_width,
        ],
        return_rocking_curve=True,
        **kwargs,
    )
    imageno = int(np.argmin(np.abs(fscan.axis - result["axis_com"])))
    return RefinedPeak(
        candidate=candidate,
        xy=np.asarray(result["xy"], dtype=float),
        imageno=imageno,
        axis_value=float(result["axis_com"]),
        value=candidate.value,
        rocking_curve=result.get("rocking_curve"),
        rocking_axis=result.get("rocking_axis"),
    )


def estimate_rocking_intensity(peak, peak_half_width=1):
    """Estimate background-subtracted integrated rocking intensity.

    :param RefinedPeak peak:
        Refined peak with ``rocking_curve`` and ``rocking_axis`` populated by
        :func:`refine_peak_3d`.
    :param int peak_half_width:
        Half-width in curve samples around the rocking maximum used as the
        signal window.
    :returns:
        Dictionary with ``intensity``, ``background``, and window metadata, or
        ``None`` if no finite rocking curve is available.
    :rtype: dict or None
    """
    if peak.rocking_curve is None:
        return None
    curve = np.asarray(peak.rocking_curve, dtype=float)
    finite = np.isfinite(curve)
    if not np.any(finite):
        return None
    peak_index = int(np.nanargmax(curve))
    start = max(0, peak_index - peak_half_width)
    stop = min(curve.size, peak_index + peak_half_width + 1)
    signal_mask = np.zeros(curve.size, dtype=bool)
    signal_mask[start:stop] = True
    background_mask = finite & ~signal_mask
    signal = curve[signal_mask & finite]
    if signal.size == 0:
        return None
    if np.any(background_mask):
        background_values = curve[background_mask]
        background = float(np.nanmedian(background_values))
        mad = float(np.nanmedian(np.abs(background_values - background)))
        background_scale = 1.4826 * mad
    else:
        background = 0.0
        background_scale = np.nan
    if not np.isfinite(background_scale) or background_scale <= np.finfo(float).eps:
        background_scale = float(np.nanstd(curve[background_mask]))
    if not np.isfinite(background_scale) or background_scale <= np.finfo(float).eps:
        background_scale = np.finfo(float).eps
    intensity = float(np.nansum(signal - background))
    prominence = float(curve[peak_index] - background)
    return {
        "intensity": intensity,
        "background": background,
        "background_scale": background_scale,
        "prominence": prominence,
        "prominence_z": prominence / background_scale,
        "peak_index": peak_index,
        "start": start,
        "stop": stop,
        "points": int(signal.size),
    }


def q_phi_from_peak(peak, detector, angle_resolver, chi, phi, k):
    """Calculate the measured Q vector for a refined peak.

    :param RefinedPeak peak:
        Refined detector and image position.
    :param DetectorCalibration.Detector2D_SXRD detector:
        Detector geometry.
    :param callable angle_resolver:
        Callable returning ``(mu, omega)`` in rad for an image index.
    :param float chi:
        Diffractometer chi angle in rad.
    :param float phi:
        Diffractometer phi angle in rad.
    :param float k:
        Incident wave-vector magnitude in Angstrom^-1.
    :returns:
        Momentum transfer in the phi frame in Angstrom^-1.
    :rtype: numpy.ndarray
    """
    mu, omega = angle_resolver(peak.imageno)
    gamma, delta = detector.surfaceAnglesPoint(
        np.array([peak.xy[1]]), np.array([peak.xy[0]]), mu
    )
    pos = [mu, float(delta[0]), float(gamma[0]), omega, chi, phi]
    return HKLVlieg.calculate_q_phi(pos, k).flatten()


def _deduplicate_structure_equivalent(hkls, structure_factor, norms, decimals=6):
    """Drop sign-inverted reflections with equal structure factor and Q norm.

    This is intentionally conservative: only reflections with the same
    absolute Miller indices, rounded complex structure factor, and
    reciprocal-vector norm are merged. Reflections with equal structure factor
    but different Q, or with unrelated HKL directions, remain separate.
    The first occurrence is kept, so callers should sort before deduplication
    when they want a specific representative.
    """
    if len(hkls) == 0:
        return hkls, structure_factor, norms
    keys = np.column_stack(
        [
            np.round(np.abs(hkls), decimals),
            np.round(np.real(structure_factor), decimals),
            np.round(np.imag(structure_factor), decimals),
            np.round(norms, decimals),
        ]
    )
    _unique, idx = np.unique(keys, axis=0, return_index=True)
    idx = np.sort(idx)
    return hkls[idx], structure_factor[idx], norms[idx]


def allowed_bragg_by_qnorm(xtal, qnorm, tolerance, max_q=None, qnorm_scale=1.0):
    """Find allowed Bragg reflections with reciprocal-vector norm near ``qnorm``.

    :param CTRcalc.UnitCell or CTRcalc.SXRDCrystal xtal:
        Crystal used to calculate allowed reflections.
    :param float qnorm:
        Observed momentum-transfer norm in Angstrom^-1.
    :param float tolerance:
        Maximum absolute norm mismatch in Angstrom^-1.
    :param float max_q:
        Optional maximum reflection norm in Angstrom^-1.
    :param float qnorm_scale:
        Multiplicative scale applied to predicted reciprocal-vector norms
        before comparing them with the measured Q norm. Values larger than
        one correspond to direct lattice constants smaller than configured.
    :returns:
        Tuple ``(hkls, norms)`` sorted by norm mismatch.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    search_q = max_q if max_q is not None else qnorm + tolerance
    hkls, (structure_factor, norms) = rn.allowedReflections_G(
        xtal, maxQ=search_q, negative=True, returnF=True
    )
    finite = np.isfinite(norms)
    nonzero = np.linalg.norm(hkls, axis=1) > 0
    scaled_norms = norms * qnorm_scale
    close = np.abs(scaled_norms - qnorm) <= tolerance
    mask = finite & nonzero & close
    hkls = hkls[mask]
    structure_factor = structure_factor[mask]
    norms = norms[mask]
    order = np.argsort(np.abs(scaled_norms[mask] - qnorm))
    hkls = hkls[order]
    structure_factor = structure_factor[order]
    norms = norms[order]
    hkls, _structure_factor, norms = _deduplicate_structure_equivalent(
        hkls, structure_factor, norms
    )
    return hkls, norms


def allowed_bragg_in_q_region(xtal, qnorm, tolerance, max_q=None, qnorm_scale=1.0):
    """Return allowed Bragg reflections in a small Q-norm shell.

    :param CTRcalc.UnitCell or CTRcalc.SXRDCrystal xtal:
        Crystal used to calculate allowed reflections and structure factors.
    :param float qnorm:
        Center Q norm in Angstrom^-1.
    :param float tolerance:
        Half-width of the Q-norm shell in Angstrom^-1.
    :param float max_q:
        Optional global maximum Q norm in Angstrom^-1.
    :param float qnorm_scale:
        Multiplicative scale applied to predicted reciprocal-vector norms
        before comparing them with the measured Q norm. Values larger than
        one correspond to direct lattice constants smaller than configured.
    :returns:
        ``(hkls, intensity, norms)`` sorted by decreasing ``abs(F)**2`` and
        with likely symmetry-equivalent reflections removed.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    search_q = max_q if max_q is not None else qnorm + tolerance
    hkls, (structure_factor, norms) = rn.allowedReflections_G(
        xtal, maxQ=search_q, negative=True, returnF=True
    )
    nonzero = np.linalg.norm(hkls, axis=1) > 0
    scaled_norms = norms * qnorm_scale
    close = np.abs(scaled_norms - qnorm) <= tolerance
    mask = nonzero & close
    hkls = hkls[mask]
    structure_factor = structure_factor[mask]
    norms = norms[mask]
    intensity = np.abs(structure_factor) ** 2
    order = np.argsort(intensity)[::-1]
    hkls = hkls[order]
    structure_factor = structure_factor[order]
    norms = norms[order]
    hkls, structure_factor, norms = _deduplicate_structure_equivalent(
        hkls, structure_factor, norms
    )
    intensity = np.abs(structure_factor) ** 2
    order = np.argsort(intensity)[::-1]
    return hkls[order], intensity[order], norms[order]


def allowed_bragg_same_q(xtal, qnorm, tolerance=1e-6, max_q=None):
    """Return all allowed Bragg reflections with the same configured Q norm.

    :param CTRcalc.UnitCell or CTRcalc.SXRDCrystal xtal:
        Crystal used to calculate allowed reflections and structure factors.
    :param float qnorm:
        Configured reciprocal-vector norm in Angstrom^-1.
    :param float tolerance:
        Absolute norm tolerance in Angstrom^-1 for grouping same-Q
        reflections.
    :param float max_q:
        Optional global maximum Q norm in Angstrom^-1.
    :returns:
        ``(hkls, intensity, norms)`` sorted by decreasing ``abs(F)**2`` and
        with likely symmetry-equivalent reflections removed.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    search_q = max_q if max_q is not None else qnorm + tolerance
    hkls, (structure_factor, norms) = rn.allowedReflections_G(
        xtal, maxQ=search_q, negative=True, returnF=True
    )
    nonzero = np.linalg.norm(hkls, axis=1) > 0
    close = np.abs(norms - qnorm) <= tolerance
    mask = nonzero & close
    hkls = hkls[mask]
    structure_factor = structure_factor[mask]
    norms = norms[mask]
    intensity = np.abs(structure_factor) ** 2
    order = np.argsort(intensity)[::-1]
    hkls = hkls[order]
    structure_factor = structure_factor[order]
    norms = norms[order]
    hkls, structure_factor, norms = _deduplicate_structure_equivalent(
        hkls, structure_factor, norms
    )
    intensity = np.abs(structure_factor) ** 2
    order = np.argsort(intensity)[::-1]
    return hkls[order], intensity[order], norms[order]


def allowed_bragg_with_intensity(xtal, max_q):
    """Return allowed Bragg reflections sorted by predicted intensity.

    :param CTRcalc.UnitCell or CTRcalc.SXRDCrystal xtal:
        Crystal used to calculate allowed reflections and structure factors.
    :param float max_q:
        Maximum reciprocal-vector norm in Angstrom^-1.
    :returns:
        ``(hkls, intensity, norms)`` sorted by decreasing ``abs(F)**2``.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    hkls, (structure_factor, norms) = rn.allowedReflections_G(
        xtal, maxQ=max_q, negative=True, returnF=True
    )
    nonzero = np.linalg.norm(hkls, axis=1) > 0
    hkls = hkls[nonzero]
    structure_factor = structure_factor[nonzero]
    norms = norms[nonzero]
    intensity = np.abs(structure_factor) ** 2
    order = np.argsort(intensity)[::-1]
    hkls = hkls[order]
    structure_factor = structure_factor[order]
    norms = norms[order]
    hkls, structure_factor, norms = _deduplicate_structure_equivalent(
        hkls, structure_factor, norms
    )
    intensity = np.abs(structure_factor) ** 2
    order = np.argsort(intensity)[::-1]
    return hkls[order], intensity[order], norms[order]


def estimate_qnorm_scale(
    peaks,
    detector,
    angle_resolver,
    ub_calculator,
    xtal,
    chi,
    phi,
    tolerance,
    max_q=None,
    hkl_filter=None,
):
    """Estimate systematic reciprocal-length scale from unmatched peaks.

    :param list[RefinedPeak] peaks:
        Refined candidate peaks that did not produce a confirmed seed.
    :param DetectorCalibration.Detector2D_SXRD detector:
        Detector geometry.
    :param callable angle_resolver:
        Callable returning ``(mu, omega)`` in rad for an image index.
    :param HKLVlieg.UBCalculator ub_calculator:
        UB calculator providing lattice, energy, and wave-vector magnitude.
    :param CTRcalc.UnitCell or CTRcalc.SXRDCrystal xtal:
        Crystal used to search allowed Bragg reflections.
    :param float chi:
        Diffractometer chi angle in rad.
    :param float phi:
        Diffractometer phi angle in rad.
    :param float tolerance:
        Broad Q-norm half-width in Angstrom^-1 used for tentative matches.
    :param float max_q:
        Optional maximum allowed reflection norm in Angstrom^-1.
    :param callable hkl_filter:
        Optional vectorized predicate ``hkl_filter(peak, hkls, norms)`` that
        returns a boolean mask. If provided, only accepted reflections are used
        for the closest-Q tentative match.
    :returns:
        Dictionary with ``qnorm_scale`` and diagnostic match information, or
        ``None`` if too few finite tentative matches are available.
    :rtype: dict or None
    """
    ratios = []
    matches = []
    for peak in peaks:
        q_phi = q_phi_from_peak(
            peak, detector, angle_resolver, chi, phi, ub_calculator.getK()
        )
        qnorm = float(np.linalg.norm(q_phi))
        hkls, _intensity, norms = allowed_bragg_in_q_region(
            xtal, qnorm, tolerance, max_q=max_q
        )
        if hkl_filter is not None and len(hkls) > 0:
            keep = np.asarray(hkl_filter(peak, hkls, norms), dtype=bool)
            if keep.shape != norms.shape:
                raise ValueError("hkl_filter must return one mask value per HKL")
            hkls = hkls[keep]
            norms = norms[keep]
        if len(norms) == 0:
            continue
        index = int(np.argmin(np.abs(norms - qnorm)))
        norm = float(norms[index])
        if norm <= np.finfo(float).eps:
            continue
        ratio = qnorm / norm
        ratios.append(ratio)
        matches.append(
            {
                "peak": peak,
                "qnorm": qnorm,
                "hkl": np.asarray(hkls[index], dtype=float),
                "configured_norm": norm,
                "qnorm_scale": ratio,
            }
        )
    if len(ratios) < 2:
        return None
    ratios = np.asarray(ratios, dtype=float)
    median = float(np.median(ratios))
    scatter = float(1.4826 * np.median(np.abs(ratios - median)))
    return {
        "qnorm_scale": median,
        "direct_lattice_scale": 1.0 / median,
        "scatter": scatter,
        "matches": matches,
    }


def seed_u_from_single_reflection(ub_calculator, pos, hkl):
    """Build a UB seed using the existing one-reflection z-mode solver.

    :param HKLVlieg.UBCalculator ub_calculator:
        UB calculator initialized with the desired default geometry.
    :param numpy.ndarray pos:
        Six-circle angles ``[alpha, delta, gamma, omega, chi, phi]`` in rad.
    :param numpy.ndarray hkl:
        Candidate Miller indices in r.l.u.
    :returns:
        Candidate orientation matrix from
        :meth:`HKLVlieg.UBCalculator.zmodeUSingleRefl`.
    :rtype: numpy.ndarray
    """
    trial = copy.deepcopy(ub_calculator)
    trial.zmodeUSingleRefl(pos, hkl)
    return np.asarray(trial.getU(), dtype=float)


def score_u_against_observed_q(U, lattice, observed_q, norm_tolerance):
    """Score a U matrix by how closely observed Q vectors index as allowed hkls.

    :param numpy.ndarray U:
        Candidate orientation matrix.
    :param HKLVlieg.Lattice lattice:
        Lattice defining the B matrix.
    :param numpy.ndarray observed_q:
        Observed Q vectors in the phi frame, in Angstrom^-1, shaped ``(N, 3)``.
    :param float norm_tolerance:
        Q-norm tolerance in Angstrom^-1 used to scale the score.
    :returns:
        Median indexing residual in reciprocal-lattice units.
    :rtype: float
    """
    observed_q = np.atleast_2d(np.asarray(observed_q, dtype=float))
    hkl_float = np.linalg.solve(U @ lattice.B_mat, observed_q.T).T
    hkl_round = np.rint(hkl_float)
    q_round = (U @ lattice.B_mat @ hkl_round.T).T
    q_error = np.linalg.norm(q_round - observed_q, axis=1)
    scaled_q_error = q_error / max(norm_tolerance, np.finfo(float).eps)
    hkl_error = np.linalg.norm(hkl_float - hkl_round, axis=1)
    return float(np.median(hkl_error + scaled_q_error))


def rank_ub_seeds(
    peaks,
    detector,
    angle_resolver,
    ub_calculator,
    xtal,
    chi,
    phi,
    qnorm_tolerance=0.05,
    max_q=None,
    qnorm_scale=1.0,
):
    """Rank automatic UB seeds from refined Bragg-peak candidates.

    :param list[RefinedPeak] peaks:
        Refined peak candidates.
    :param DetectorCalibration.Detector2D_SXRD detector:
        Detector geometry.
    :param callable angle_resolver:
        Callable returning ``(mu, omega)`` in rad for an image index.
    :param HKLVlieg.UBCalculator ub_calculator:
        UB calculator providing lattice, energy, and the current/default U.
    :param CTRcalc.UnitCell or CTRcalc.SXRDCrystal xtal:
        Crystal used to search allowed Bragg reflections.
    :param float chi:
        Diffractometer chi angle in rad.
    :param float phi:
        Diffractometer phi angle in rad.
    :param float qnorm_tolerance:
        Q-norm matching tolerance in Angstrom^-1.
    :param float max_q:
        Optional maximum allowed reflection norm in Angstrom^-1.
    :param float qnorm_scale:
        Multiplicative scale applied to predicted reciprocal-vector norms
        before Q-shell matching.
    :returns:
        Candidate UB seeds sorted by increasing score.
    :rtype: list[UBSeed]
    """
    if not peaks:
        return []
    positions = []
    observed_q = []
    for peak in peaks:
        mu, omega = angle_resolver(peak.imageno)
        gamma, delta = detector.surfaceAnglesPoint(
            np.array([peak.xy[1]]), np.array([peak.xy[0]]), mu
        )
        pos = np.array([mu, float(delta[0]), float(gamma[0]), omega, chi, phi])
        positions.append(pos)
        observed_q.append(HKLVlieg.calculate_q_phi(pos, ub_calculator.getK()).flatten())
    observed_q = np.asarray(observed_q)
    seeds = []
    for peak, pos, q_phi in zip(peaks, positions, observed_q):
        qnorm = float(np.linalg.norm(q_phi))
        hkls, norms = allowed_bragg_by_qnorm(
            xtal,
            qnorm,
            qnorm_tolerance,
            max_q=max_q,
            qnorm_scale=qnorm_scale,
        )
        for hkl, norm in zip(hkls, norms):
            try:
                U = seed_u_from_single_reflection(ub_calculator, pos, hkl)
            except Exception:
                logger.debug(
                    "Cannot calculate one-reflection UB seed for hkl=%s",
                    hkl,
                    exc_info=True,
                )
                continue
            score = score_u_against_observed_q(
                U, ub_calculator.lattice, observed_q, qnorm_tolerance
            )
            seeds.append(
                UBSeed(
                    U=U,
                    hkl=np.asarray(hkl, dtype=float),
                    peak=peak,
                    score=score,
                    norm_mismatch=abs(float(norm) * qnorm_scale - qnorm),
                )
            )
    return sorted(seeds, key=lambda seed: (seed.score, seed.norm_mismatch))

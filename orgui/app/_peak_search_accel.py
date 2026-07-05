"""Compatibility wrapper for auto-Bragg peak-search acceleration."""

from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path

import numpy as np


def _import_cpp_backend():
    """Import the packaged extension or an in-tree Meson build artifact."""
    try:
        return importlib.import_module("orgui.app._peak_search_cpp")
    except ModuleNotFoundError as package_error:
        repo_root = Path(__file__).resolve().parents[2]
        candidates = sorted((repo_root / "build").glob("cp*/_peak_search_cpp*.so"))
        if not candidates:
            raise package_error
        extension_path = candidates[-1]
        spec = importlib.util.spec_from_file_location(
            "_peak_search_cpp",
            extension_path,
        )
        if spec is None or spec.loader is None:
            raise package_error
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module


try:
    _CPP_BACKEND = _import_cpp_backend()
except Exception:
    _CPP_BACKEND = None

HAS_CPP_BACKEND = _CPP_BACKEND is not None


class MaximumResult:
    """Peak maximum result returned by the fallback implementation."""

    def __init__(self, valid=False, value=np.nan, x=np.nan, y=np.nan):
        self.valid = valid
        self.value = value
        self.x = x
        self.y = y


class Candidate:
    """Peak candidate result returned by the fallback detector."""

    def __init__(
        self,
        index,
        value,
        x,
        y,
        sharpness=np.nan,
        derivative_sharpness=np.nan,
        prominence=np.nan,
    ):
        self.index = index
        self.value = value
        self.x = x
        self.y = y
        self.sharpness = sharpness
        self.derivative_sharpness = derivative_sharpness
        self.prominence = prominence


def _robust_location_scale(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan, np.nan
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    return median, 1.4826 * mad


def _robust_z(value, values):
    median, scale = _robust_location_scale(values)
    if not np.isfinite(median):
        return np.nan
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = np.std(values)
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        scale = np.finfo(float).eps
    return float((value - median) / scale)


def _near_mask(mask, y, x, distance):
    if mask is None or distance <= 0:
        return False
    y0 = max(0, y - distance)
    y1 = min(mask.shape[0], y + distance + 1)
    x0 = max(0, x - distance)
    x1 = min(mask.shape[1], x + distance + 1)
    return bool(np.any(mask[y0:y1, x0:x1]))


def _masked_maximum_fallback(image, mask=None, background=None, mask_distance=0):
    arr = np.asarray(image, dtype=float)
    if background is not None and np.shape(background) == arr.shape:
        arr = arr - np.asarray(background, dtype=float)
    valid = np.isfinite(arr)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        valid &= ~mask
        if mask_distance > 0:
            rows, cols = np.nonzero(valid)
            allowed = np.zeros_like(valid, dtype=bool)
            for y, x in zip(rows, cols):
                if not _near_mask(mask, int(y), int(x), mask_distance):
                    allowed[y, x] = True
            valid = allowed
    if not np.any(valid):
        return MaximumResult()
    values = np.where(valid, arr, -np.inf)
    flat = int(np.nanargmax(values))
    y, x = np.unravel_index(flat, values.shape)
    return MaximumResult(True, float(values[y, x]), float(x), float(y))


def _masked_roi_sum_fallback(
    image,
    mask=None,
    background=None,
    y0=0,
    y1=0,
    x0=0,
    x1=0,
):
    arr = np.asarray(image, dtype=float)
    if background is not None and np.shape(background) == arr.shape:
        arr = arr - np.asarray(background, dtype=float)
    sub = arr[y0:y1, x0:x1]
    if mask is not None:
        submask = np.asarray(mask, dtype=bool)[y0:y1, x0:x1]
        sub = np.where(submask, np.nan, sub)
    return float(np.nansum(sub))


class PeakCandidateDetector:
    """Pure Python fallback for the C++ streaming candidate detector."""

    def __init__(
        self,
        n_frames,
        burn_in,
        history,
        min_history,
        level_z,
        derivative_z,
        lookahead,
        min_prominence_z,
        refractory,
    ):
        self.n_frames = int(n_frames)
        self.burn_in = int(burn_in)
        self.history = int(history)
        self.min_history = int(min_history)
        self.level_z = float(level_z)
        self.derivative_z = float(derivative_z)
        self.lookahead = int(lookahead)
        self.min_prominence_z = float(min_prominence_z)
        self.refractory = int(refractory)
        self.next_allowed = self.burn_in
        self.eval_index = 0
        self.highest_seen = -1
        self.frames = [None] * self.n_frames

    def push(self, index, valid, value, x, y, finish=False):
        index = int(index)
        if valid:
            self.frames[index] = {
                "valid": True,
                "value": float(value),
                "x": float(x),
                "y": float(y),
                "sharpness": np.nan,
                "derivative_sharpness": np.nan,
                "prominence": np.nan,
            }
        else:
            self.frames[index] = {"valid": False, "value": np.nan}
        self.highest_seen = max(self.highest_seen, index)
        return self._drain(finish)

    def finish(self):
        return self._drain(True)

    def _ready(self, index, finish):
        if index >= self.n_frames or self.frames[index] is None:
            return False
        if finish:
            return True
        required = min(self.n_frames - 1, index + self.lookahead)
        return self.highest_seen >= required

    def _baseline(self, index):
        start = max(0, index - self.history)
        values = [
            self.frames[i]["value"]
            for i in range(start, index)
            if self.frames[i] is not None
            and self.frames[i]["valid"]
            and np.isfinite(self.frames[i]["value"])
        ]
        return np.asarray(values, dtype=float)

    def _drain(self, finish):
        candidates = []
        while self._ready(self.eval_index, finish):
            maximum = self.frames[self.eval_index]
            if not maximum["valid"] or not np.isfinite(maximum["value"]):
                self.eval_index += 1
                continue
            baseline = self._baseline(self.eval_index)
            if self.eval_index < self.next_allowed or baseline.size < self.min_history:
                self.eval_index += 1
                continue
            diffs = np.diff(baseline)
            current_diff = maximum["value"] - baseline[-1]
            level_score = _robust_z(maximum["value"], baseline)
            diff_score = _robust_z(current_diff, diffs)
            maximum["sharpness"] = level_score
            maximum["derivative_sharpness"] = diff_score
            if (
                not np.isfinite(level_score)
                or not np.isfinite(diff_score)
                or level_score < self.level_z
                or diff_score < self.derivative_z
            ):
                self.eval_index += 1
                continue
            candidate_index = self.eval_index
            candidate_frame = maximum
            for look_idx in range(
                self.eval_index + 1,
                min(self.n_frames, self.eval_index + self.lookahead + 1),
            ):
                look = self.frames[look_idx]
                if (
                    look is not None
                    and look["valid"]
                    and look["value"] > candidate_frame["value"]
                ):
                    candidate_index = look_idx
                    candidate_frame = look
            median, scale = _robust_location_scale(baseline)
            if not np.isfinite(scale) or scale <= np.finfo(float).eps:
                scale = np.std(baseline)
            if not np.isfinite(scale) or scale <= np.finfo(float).eps:
                scale = np.finfo(float).eps
            prominence = candidate_frame["value"] - median
            candidate_frame["prominence"] = float(prominence)
            if prominence / scale >= self.min_prominence_z:
                candidate_frame["sharpness"] = _robust_z(
                    candidate_frame["value"], baseline
                )
                candidate_baseline = self._baseline(candidate_index)
                prev_value = candidate_baseline[-1]
                candidate_frame["derivative_sharpness"] = _robust_z(
                    candidate_frame["value"] - prev_value, diffs
                )
                candidates.append(
                    Candidate(
                        candidate_index,
                        candidate_frame["value"],
                        candidate_frame["x"],
                        candidate_frame["y"],
                        candidate_frame["sharpness"],
                        candidate_frame["derivative_sharpness"],
                        candidate_frame["prominence"],
                    )
                )
                self.next_allowed = candidate_index + self.refractory + 1
            self.eval_index = max(self.eval_index + 1, candidate_index + 1)
        return candidates


if _CPP_BACKEND is not None:
    masked_maximum = _CPP_BACKEND.masked_maximum
    masked_roi_sum = _CPP_BACKEND.masked_roi_sum
    PeakCandidateDetector = _CPP_BACKEND.PeakCandidateDetector
else:
    masked_maximum = _masked_maximum_fallback
    masked_roi_sum = _masked_roi_sum_fallback

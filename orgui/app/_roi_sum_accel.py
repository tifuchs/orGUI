"""Compatibility wrapper for compiled ROI summing kernels."""

from __future__ import annotations

import importlib
import importlib.util
import os
from pathlib import Path

_CPP_BACKEND = None
_NUMBA_BACKEND = None
_NUMBA_IMPORT_ATTEMPTED = False
_NUMBA_IMPORT_ERROR = None
_ACCEL_BACKEND_ENV_VAR = "ORGUI_ACCEL_BACKEND"


class _NoAccelerationBackend:
    @staticmethod
    def _unavailable(*args, **kwargs):
        raise RuntimeError("ROI acceleration is disabled")

    calcBgSub = _unavailable
    calcMaxSum = _unavailable
    calcMaxSum_bg = _unavailable
    processImage_Carr = _unavailable
    processImage_bg_Carr = _unavailable


def _import_cpp_backend():
    """Import the packaged extension or an in-tree Meson build artifact."""
    try:
        return importlib.import_module("orgui.app._roi_sum_cpp")
    except ModuleNotFoundError as package_error:
        repo_root = Path(__file__).resolve().parents[2]
        candidates = sorted((repo_root / "build").glob("cp*/_roi_sum_cpp*.so"))
        if not candidates:
            raise package_error
        extension_path = candidates[-1]
        spec = importlib.util.spec_from_file_location(
            "_roi_sum_cpp",
            extension_path,
        )
        if spec is None or spec.loader is None:
            raise package_error
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module


def _load_cpp_backend():
    global _CPP_BACKEND

    if _CPP_BACKEND is None:
        _CPP_BACKEND = _import_cpp_backend()
    return _CPP_BACKEND


def _load_numba_backend():
    global _NUMBA_BACKEND
    global _NUMBA_IMPORT_ATTEMPTED
    global _NUMBA_IMPORT_ERROR

    if not _NUMBA_IMPORT_ATTEMPTED:
        _NUMBA_IMPORT_ATTEMPTED = True
        try:
            _NUMBA_BACKEND = importlib.import_module("orgui.app._roi_sum_numba")
            _NUMBA_IMPORT_ERROR = None
        except Exception as exc:
            _NUMBA_BACKEND = None
            _NUMBA_IMPORT_ERROR = exc
    return _NUMBA_BACKEND


def _default_backend():
    try:
        _load_cpp_backend()
    except Exception:
        return "numpy", _NoAccelerationBackend
    return "cpp", _CPP_BACKEND


def _select_counter_backend():
    requested = os.environ.get(_ACCEL_BACKEND_ENV_VAR)
    if requested is None:
        return _default_backend()
    backend = requested.strip().lower()
    if backend == "cpp":
        return "cpp", _load_cpp_backend()
    if backend == "numpy":
        return "numpy", _NoAccelerationBackend
    if backend == "numba":
        numba_backend = _load_numba_backend()
        if numba_backend is None:
            message = (
                f"ROI Numba acceleration backend is not available. Set "
                f"{_ACCEL_BACKEND_ENV_VAR}=cpp or install numba."
            )
            if _NUMBA_IMPORT_ERROR is not None:
                message += f" Import failed: {_NUMBA_IMPORT_ERROR}"
            raise RuntimeError(message)
        return "numba", numba_backend
    raise ValueError(
        f"{_ACCEL_BACKEND_ENV_VAR} must be 'cpp', 'numba', or 'numpy', "
        f"got {backend!r}"
    )


def _bind_backend(name, backend):
    global ROI_ACCEL_BACKEND
    global HAS_NUMBA_BACKEND
    global HAS_ACCEL_BACKEND
    global calcBgSub
    global calcMaxSum
    global calcMaxSum_bg
    global processImage_Carr
    global processImage_bg_Carr

    ROI_ACCEL_BACKEND = name
    HAS_NUMBA_BACKEND = name == "numba"
    HAS_ACCEL_BACKEND = name != "numpy"
    calcBgSub = backend.calcBgSub
    calcMaxSum = backend.calcMaxSum
    calcMaxSum_bg = backend.calcMaxSum_bg
    processImage_Carr = backend.processImage_Carr
    processImage_bg_Carr = backend.processImage_bg_Carr


def set_accel_backend(backend):
    """Select the process-global ROI acceleration backend."""
    backend = backend.strip().lower()
    if backend == "cpp":
        selected = _load_cpp_backend()
    elif backend == "numba":
        selected = _load_numba_backend()
        if selected is None:
            message = "ROI Numba acceleration backend is not available."
            if _NUMBA_IMPORT_ERROR is not None:
                message += f" Import failed: {_NUMBA_IMPORT_ERROR}"
            raise RuntimeError(message)
    elif backend == "numpy":
        selected = _NoAccelerationBackend
    else:
        raise ValueError("ROI acceleration backend must be 'cpp', 'numba', or 'numpy'")
    _bind_backend(backend, selected)


ROI_ACCEL_BACKEND, _counter_backend = _select_counter_backend()

try:
    _cpp_backend = _load_cpp_backend()
except Exception:
    _cpp_backend = None

HAS_CPP_BACKEND = _cpp_backend is not None
HAS_NUMBA_BACKEND = False
HAS_ACCEL_BACKEND = ROI_ACCEL_BACKEND != "numpy"
_bind_backend(ROI_ACCEL_BACKEND, _counter_backend)

if _cpp_backend is not None:
    interpolate_polybg_croi = _cpp_backend.interpolate_polybg_croi
    processImage_polybg_Carr = _cpp_backend.processImage_polybg_Carr

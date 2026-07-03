"""Compatibility wrapper for compiled ROI summing kernels."""

from __future__ import annotations

import importlib
import importlib.util
from pathlib import Path


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


_cpp_backend = _import_cpp_backend()

calcBgSub = _cpp_backend.calcBgSub
calcMaxSum = _cpp_backend.calcMaxSum
calcMaxSum_bg = _cpp_backend.calcMaxSum_bg
interpolate_polybg_croi = _cpp_backend.interpolate_polybg_croi
processImage_Carr = _cpp_backend.processImage_Carr
processImage_bg_Carr = _cpp_backend.processImage_bg_Carr
processImage_polybg_Carr = _cpp_backend.processImage_polybg_Carr

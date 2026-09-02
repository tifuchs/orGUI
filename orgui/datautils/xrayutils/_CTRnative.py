"""Load the optional native CTR extension without importing model classes."""

import importlib
import importlib.util
from pathlib import Path


def _import_cpp_accel():
    """Import the installed or in-tree native CTR extension."""
    try:
        return importlib.import_module("orgui.datautils.xrayutils._CTRcalc_cpp")
    except ModuleNotFoundError as package_error:
        repo_root = Path(__file__).resolve().parents[3]
        candidates = sorted(
            (
                *(repo_root / "build").glob("cp*/_CTRcalc_cpp*.so"),
                *(repo_root / "build").glob("cp*/_CTRcalc_cpp*.pyd"),
            )
        )
        if not candidates:
            raise package_error
        extension_path = candidates[-1]
        spec = importlib.util.spec_from_file_location(
            "_CTRcalc_cpp",
            extension_path,
        )
        if spec is None or spec.loader is None:
            raise package_error
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module


try:
    _CTRcalc_cpp = _import_cpp_accel()
    HAS_CPP_ACCEL = True
except Exception:
    _CTRcalc_cpp = None
    HAS_CPP_ACCEL = False

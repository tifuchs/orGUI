"""Runtime detection of the compiled-equation strided-coordinate fault."""

import logging

import numpy as np
import pytest

from .. import _pyfai_compat as compat


def test_probe_detects_ignored_coordinate_strides(monkeypatch):
    """A simulated compiled evaluator reads adjacent values, not components."""
    def _broken(x, y, z, wavelength):
        assert x.strides[-1] == 3 * x.itemsize
        flat = x.base.ravel()
        return np.arctan2(np.hypot(flat[:8], flat[1:9]), flat[2:10]).reshape(2, 4)

    monkeypatch.setattr(compat.units.TTH_RAD, "equation", _broken)
    assert not compat._strided_angles_are_safe()


def test_probe_accepts_correct_equations(monkeypatch):
    """Both angle references use independent NumPy operations, in radians."""
    monkeypatch.setattr(
        compat.units.TTH_RAD, "equation",
        lambda x, y, z, wavelength: np.arctan2(np.hypot(x, y), z),
    )
    monkeypatch.setattr(
        compat.units.CHI_RAD, "equation",
        lambda x, y, z, wavelength: np.arctan2(y, x),
    )
    assert compat._strided_angles_are_safe()


def test_probe_checks_azimuth_even_when_two_theta_is_correct(monkeypatch):
    """The polarization also needs an independently valid azimuth array."""
    monkeypatch.setattr(
        compat.units.TTH_RAD, "equation",
        lambda x, y, z, wavelength: np.arctan2(np.hypot(x, y), z),
    )
    monkeypatch.setattr(
        compat.units.CHI_RAD, "equation",
        lambda x, y, z, wavelength: np.full_like(x, np.nan),
    )
    assert not compat._strided_angles_are_safe()


@pytest.mark.parametrize("safe", [False, True])
def test_startup_probe_runs_once_and_warns_only_on_failure(monkeypatch, caplog, safe):
    """GUI, CLI and worker imports share the once-per-process result."""
    calls = []

    def _probe():
        calls.append(True)
        return safe

    monkeypatch.setattr(compat, "_strided_angles_are_safe", _probe)
    compat._ensure_pyfai_array_safety.cache_clear()
    try:
        with caplog.at_level(logging.WARNING, logger=compat.__name__):
            assert compat._ensure_pyfai_array_safety() is safe
            assert compat._ensure_pyfai_array_safety() is safe
        assert len(calls) == 1
        assert ("Environment error" in caplog.text) is not safe
        if not safe:
            record = next(r for r in caplog.records if r.name == compat.__name__)
            assert record.show_dialog is True
            assert record.levelno == logging.WARNING
            assert "Try reinstalling numexpr" in record.getMessage()
            assert "best-effort workaround" in record.getMessage()
            assert "not guaranteed" in record.getMessage()
    finally:
        compat._ensure_pyfai_array_safety.cache_clear()

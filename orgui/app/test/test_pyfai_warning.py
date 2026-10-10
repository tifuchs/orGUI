"""Route numerical environment warnings through the GUI/CLI log handlers."""

from types import SimpleNamespace

import pytest

pytest.importorskip("silx.gui.qt")

from orgui import logger_utils
from orgui.datautils.xrayutils import _pyfai_compat as compat


@pytest.mark.parametrize("mode", ["gui", "cli"])
def test_pyfai_environment_warning_uses_mode_handler(monkeypatch, mode):
    """GUI receives a popup payload; CLI warnings do not raise or show one."""
    popups = []
    monkeypatch.setattr(
        logger_utils, "_get_message_box_dispatcher",
        lambda: SimpleNamespace(show=popups.append),
    )
    monkeypatch.setattr(compat, "_strided_angles_are_safe", lambda: False)
    handler = (
        logger_utils.MessageBoxHandler() if mode == "gui"
        else logger_utils.CLIExceptionHandler()
    )
    compat.logger.addHandler(handler)
    compat._ensure_pyfai_array_safety.cache_clear()
    try:
        assert compat._ensure_pyfai_array_safety() is False
        assert compat._ensure_pyfai_array_safety() is False
        if mode == "gui":
            assert len(popups) == 1
            assert popups[0]["title"] == "Numerical environment error"
            assert "Try reinstalling numexpr" in popups[0]["message"]
            assert "not guaranteed" in popups[0]["message"]
        else:
            assert popups == []
    finally:
        compat.logger.removeHandler(handler)
        handler.close()
        compat._ensure_pyfai_array_safety.cache_clear()

"""The Q-plot must be wired to the validated conversion in :mod:`qconversion`.

The numerical equivalence of the Q-plot with orGUI's own QAlpha calculation is
established in ``orgui/app/test/test_q_conversion.py``. These tests only cover
the GUI wiring: which frames are offered, which one is preselected, and how the
selection is read back.
"""

import pytest

from orgui.app import qconversion

orGUI = pytest.importorskip("orgui.app.orGUI", reason="Qt bindings are required")


class _StubComboBox:
    """Minimal stand-in for the frame selector combo box."""

    def __init__(self, data):
        self._data = data

    def currentData(self):
        return self._data


class _StubOrGUI:
    """The frame accessor only reads the combo box."""

    def __init__(self, data):
        self.qFrameSelector = _StubComboBox(data)


@pytest.mark.parametrize("frame", qconversion.FRAMES)
def test_selected_frame_is_read_from_the_combo_box(frame):
    assert orGUI.orGUI._selectedQFrame(_StubOrGUI(frame)) == frame


@pytest.mark.parametrize("data", [None, "Q_nonsense", ""])
def test_unknown_selection_falls_back_to_the_default_frame(data):
    selected = orGUI.orGUI._selectedQFrame(_StubOrGUI(data))
    assert selected == qconversion.DEFAULT_FRAME


def test_default_frame_is_the_alpha_frame():
    """QAlpha is what the rest of orGUI uses, so it is the natural default."""
    assert qconversion.DEFAULT_FRAME == "Q_alpha"
    assert qconversion.DEFAULT_FRAME in qconversion.FRAMES


def test_every_frame_has_a_label_for_the_selector():
    for frame in qconversion.FRAMES:
        assert qconversion.FRAME_LABELS[frame]


def test_fiber_integrator_flag_is_taken_from_qconversion():
    """The action is only usable when FiberIntegrator can be imported."""
    assert orGUI.HAS_FIBER_INTEGRATOR == qconversion.HAS_FIBER_INTEGRATOR


class _StubAction:
    def __init__(self, checked):
        self._checked = checked

    def isChecked(self):
        return self._checked


class _RefreshStub:
    """Records the toggles the refresh performs."""

    def __init__(self, checked, reentrant=False):
        self.plotAgainstQAct = _StubAction(checked)
        self._qPlotRefreshing = False
        self._reentrant = reentrant
        self.calls = []

    def _convertImagetoQ(self, value):
        self.calls.append(value)
        if self._reentrant:
            # emulate a signal arriving while the plot is being rebuilt
            orGUI.orGUI._refreshQPlot(self)


def test_refresh_does_nothing_while_the_q_plot_is_off():
    stub = _RefreshStub(checked=False)

    orGUI.orGUI._refreshQPlot(stub)

    assert stub.calls == []


def test_refresh_rebuilds_in_place_and_stays_in_q_mode():
    """Changing the image must not drop out of the reciprocal-space view."""
    stub = _RefreshStub(checked=True)

    orGUI.orGUI._refreshQPlot(stub)

    # a single conversion replaces the previous one; the action stays checked
    assert stub.calls == [True]
    assert stub.plotAgainstQAct.isChecked()


def test_refresh_accepts_the_replot_signal_argument():
    """sigReplotRequest carries a bool, which the slot must tolerate."""
    stub = _RefreshStub(checked=True)

    orGUI.orGUI._refreshQPlot(stub, True)

    assert stub.calls == [True]


def test_refresh_is_not_reentrant():
    stub = _RefreshStub(checked=True, reentrant=True)

    orGUI.orGUI._refreshQPlot(stub)

    assert stub.calls == [True]
    assert stub._qPlotRefreshing is False


def test_suspended_q_plot_batches_into_one_rebuild():
    """Untoggling max/sum before replotting must not rebuild several times."""
    stub = _RefreshStub(checked=True)

    with orGUI.orGUI._suspendedQPlot(stub):
        orGUI.orGUI._refreshQPlot(stub)
        orGUI.orGUI._refreshQPlot(stub)
        assert stub.calls == []
    orGUI.orGUI._refreshQPlot(stub)

    assert stub.calls == [True]
    assert stub._qPlotRefreshing is False

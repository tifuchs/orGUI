"""Regression coverage for the integration-curve action dialog."""

import pytest
from silx.gui import qt

from orgui.app.orGUI import QPlotDeleteWindow


@pytest.fixture(scope="module")
def _application():
    application = qt.QApplication.instance() or qt.QApplication([])
    yield application


@pytest.mark.parametrize(
    "label, action, accepted",
    [("Delete curves", "delete", True), ("Hide curves", "hide", True),
     (None, None, False)],
)
def test_curve_buttons_select_the_correct_action(_application, label, action,
                                                accepted):
    """Each button finishes the dialog with its intended result and action."""
    dialog = QPlotDeleteWindow(["visible", "hidden"], [False, True])
    try:
        assert [box.isChecked() for box in dialog.boxes] == [False, True]
        if label is None:
            button = dialog.buttons.button(qt.QDialogButtonBox.Cancel)
        else:
            button = next(button for button in dialog.buttons.buttons()
                          if button.text() == label)
        button.click()
        assert dialog.action == action
        expected = qt.QDialog.Accepted if accepted else qt.QDialog.Rejected
        assert dialog.result() == expected
    finally:
        dialog.close()
        dialog.deleteLater()
        _application.processEvents()

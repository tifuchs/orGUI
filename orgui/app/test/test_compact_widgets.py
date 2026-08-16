"""Tests for the space saving widgets used by the scan selector panel."""

from silx.gui import qt

from orgui.app.qutils import CompactDoubleSpinBox, CompactSpinBox


_APPLICATION = None


def _application():
    """Return the running application, creating it on demand.

    The instance is kept alive in a module attribute: a garbage collected
    ``QApplication`` takes the widgets of the following test down with it.
    """
    global _APPLICATION
    if _APPLICATION is None:
        _APPLICATION = qt.QApplication.instance() or qt.QApplication([])
    return _APPLICATION


def test_compact_spin_box_keeps_the_full_value_range():
    """Compressing the widget must not change the accepted values."""
    _application()
    compact = CompactSpinBox()
    compact.setRange(1, 2147483647)
    compact.setValue(2147483647)
    assert compact.value() == 2147483647

    spin = CompactDoubleSpinBox()
    spin.setRange(-20000, 20000)
    spin.setDecimals(4)
    spin.setValue(-19999.1234)
    assert spin.value() == -19999.1234

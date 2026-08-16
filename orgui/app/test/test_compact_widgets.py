"""Tests for the space saving widgets used by the scan selector panel."""

from silx.gui import qt

from orgui.app.qutils import CompactDoubleSpinBox, CompactSpinBox, FlowLayout


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


def test_compact_double_spin_box_is_narrower_than_its_range():
    """A wide scientific range must not dictate the minimum width."""
    _application()
    reference = qt.QDoubleSpinBox()
    reference.setRange(-20000, 20000)
    reference.setDecimals(4)

    compact = CompactDoubleSpinBox()
    compact.setRange(-20000, 20000)
    compact.setDecimals(4)

    assert compact.minimumSizeHint().width() < reference.minimumSizeHint().width()
    assert compact.sizeHint().width() <= reference.sizeHint().width()
    assert compact.minimumSizeHint().width() <= compact.sizeHint().width()
    # the layout must be allowed to use the reduced minimum width
    assert compact.sizePolicy().horizontalPolicy() == qt.QSizePolicy.Preferred


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


def test_compact_spin_box_width_follows_the_font():
    """Widths are derived from font metrics, so they scale with the font."""
    _application()
    compact = CompactDoubleSpinBox()
    compact.setRange(-20000, 20000)
    compact.setDecimals(4)
    small = compact.minimumSizeHint().width()

    font = qt.QFont(compact.font())
    font.setPointSizeF(font.pointSizeF() * 2)
    compact.setFont(font)

    assert compact.minimumSizeHint().width() > small


def test_compact_spin_box_visible_chars_control_the_width():
    """Fewer visible characters must result in a narrower minimum width."""
    _application()
    compact = CompactDoubleSpinBox()
    compact.setRange(-20000, 20000)
    compact.setDecimals(4)
    wide = compact.minimumSizeHint().width()
    compact.setVisibleChars(compact.getVisibleChars() - 2)
    assert compact.minimumSizeHint().width() < wide


def test_compact_spin_box_preferred_chars_do_not_raise_the_minimum():
    """The preferred width may be tuned without pinning the minimum width."""
    _application()
    compact = CompactDoubleSpinBox()
    compact.setRange(-100, 100)
    compact.setDecimals(3)
    compact.setVisibleChars(5)
    compact.setPreferredChars(6)
    narrow = compact.minimumSizeHint().width()
    assert narrow < compact.sizeHint().width()

    compact.setPreferredChars(9)
    assert compact.minimumSizeHint().width() == narrow
    assert compact.sizeHint().width() > narrow


def test_flow_layout_minimum_width_is_the_widest_item():
    """The items are allowed to wrap instead of forcing one wide row."""
    _application()
    widget = qt.QWidget()
    layout = FlowLayout()
    boxes = [
        qt.QCheckBox("Use pixel mask"),
        qt.QCheckBox("Solid angle correction"),
        qt.QCheckBox("Lorentz correction"),
        qt.QCheckBox("Polarization correction"),
    ]
    for box in boxes:
        layout.addWidget(box)
    widget.setLayout(layout)

    widest = max(box.sizeHint().width() for box in boxes)
    total = sum(box.sizeHint().width() for box in boxes)
    minimum = layout.minimumSize().width()
    assert widest <= minimum < total

    # a narrow widget needs more rows than a wide one
    assert layout.heightForWidth(widest) > layout.heightForWidth(total)

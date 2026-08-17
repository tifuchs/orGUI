# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import logging

from silx.gui import qt
from silx.gui import icons
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
import numpy as np

logger = logging.getLogger(__name__)


class RobustHdf5TreeModel(Hdf5TreeModel):
    """HDF5 tree model which tolerates files that cannot be closed.

    :class:`Hdf5TreeModel` closes the files it owns when the corresponding
    item is removed from the model. h5py raises if the file cannot be written
    on close, e.g. if the drive holding the file was disconnected. Such an
    error aborts the removal of the item and propagates into the Qt event
    loop, which terminates orGUI (issue #23).

    This model logs the error instead and removes the item, so that the user
    can continue to work with the remaining files.
    """

    def _closeFileIfOwned(self, node):
        try:
            super()._closeFileIfOwned(node)
        except Exception:
            filename = getattr(node.obj, "filename", "")
            # logged as a warning on purpose: in cli context, error records
            # re-raise the exception, see orgui.logger_utils.
            logger.warning(
                "Closing of file %s failed. The file might be corrupted!",
                filename,
                exc_info=True,
            )


def messagebox_detailed_message(
    parent, title, text, detailed_text, icon, buttons=qt.QMessageBox.Ok
):
    diag = qt.QMessageBox(icon, title, text, buttons, parent)
    diag.setDetailedText(detailed_text)
    return diag.exec()


def critical_detailed_message(
    parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok
):
    return messagebox_detailed_message(
        parent, title, text, detailed_text, qt.QMessageBox.Critical, buttons=buttons
    )


def warning_detailed_message(
    parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok
):
    return messagebox_detailed_message(
        parent, title, text, detailed_text, qt.QMessageBox.Warning, buttons=buttons
    )


def information_detailed_message(
    parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok
):
    return messagebox_detailed_message(
        parent, title, text, detailed_text, qt.QMessageBox.Information, buttons=buttons
    )


def compactTabWidget(parent=None):
    """Return a :class:`QTabWidget` which does not force a wide panel.

    By default the minimum width of a tab widget contains every tab label
    at full length, which is expensive in the Fusion and Windows 11 styles
    because of their generous tab padding. Eliding the labels and enabling
    the scroll buttons lets the tab bar shrink; the document mode drops the
    additional frame drawn around the tab pages.

    :param QWidget parent: Optional parent widget.
    :rtype: QTabWidget
    """
    tabwidget = qt.QTabWidget(parent)
    tabwidget.setDocumentMode(True)
    tabwidget.setElideMode(qt.Qt.ElideRight)
    tabwidget.setUsesScrollButtons(True)
    tabwidget.tabBar().setExpanding(False)
    return tabwidget


def wrappingLabel(text, parent=None):
    """Return a :class:`QLabel` which is allowed to wrap its text.

    A non wrapping label reports the full text width as its minimum width
    and therefore forces the whole panel to stay at least that wide.
    Wrapping labels only require the width of their longest word, which
    lets the surrounding dock be resized freely.

    :param str text: Label text.
    :param QWidget parent: Optional parent widget.
    :rtype: QLabel
    """
    label = qt.QLabel(text, parent)
    label.setWordWrap(True)
    return label


def labelWithToolTip(text, tooltip, parent=None):
    """Return a :class:`QLabel` with a short text and an explaining tool tip.

    Short row labels keep a panel narrow; the tool tip carries the full
    wording for the user.

    :param str text: Label text.
    :param str tooltip: Full description shown on hover.
    :param QWidget parent: Optional parent widget.
    :rtype: QLabel
    """
    label = qt.QLabel(text, parent)
    label.setToolTip(tooltip)
    return label


class FlowLayout(qt.QLayout):
    """Layout which arranges its items in rows and wraps to the next row.

    It behaves like a horizontal box layout as long as the items fit next to
    each other, and falls back to fewer columns when the available width
    shrinks. The minimum width is therefore the width of the widest single
    item instead of the sum of all items, which allows panels holding several
    check boxes or small editors to be resized freely.

    Spacings default to the layout spacing of the current style, so they follow
    the style and the screen scaling.
    """

    def __init__(self, parent=None, margin=-1, hSpacing=-1, vSpacing=-1):
        super().__init__(parent)
        self.__items = []
        self.__hSpacing = hSpacing
        self.__vSpacing = vSpacing
        if margin >= 0:
            self.setContentsMargins(margin, margin, margin, margin)

    def __del__(self):
        self.__items.clear()

    def addItem(self, item):
        self.__items.append(item)

    def count(self):
        return len(self.__items)

    def itemAt(self, index):
        if 0 <= index < len(self.__items):
            return self.__items[index]
        return None

    def takeAt(self, index):
        if 0 <= index < len(self.__items):
            return self.__items.pop(index)
        return None

    def horizontalSpacing(self):
        """Return the horizontal gap between two items, in pixels.

        :rtype: int
        """
        if self.__hSpacing >= 0:
            return self.__hSpacing
        return self.__styleSpacing(qt.QStyle.PM_LayoutHorizontalSpacing)

    def verticalSpacing(self):
        """Return the vertical gap between two rows, in pixels.

        :rtype: int
        """
        if self.__vSpacing >= 0:
            return self.__vSpacing
        return self.__styleSpacing(qt.QStyle.PM_LayoutVerticalSpacing)

    def __styleSpacing(self, pixelMetric):
        parent = self.parent()
        if parent is None:
            return 0
        if parent.isWidgetType():
            return parent.style().pixelMetric(pixelMetric, None, parent)
        return parent.spacing()

    def expandingDirections(self):
        # Qt5 bindings expose the flags type as Qt.Orientations, which no longer
        # exists in Qt6. Qt.Orientation(0) yields the empty flag set on both.
        return qt.Qt.Orientation(0)

    def hasHeightForWidth(self):
        return True

    def heightForWidth(self, width):
        return self.__layout(qt.QRect(0, 0, width, 0), test=True)

    def setGeometry(self, rect):
        super().setGeometry(rect)
        self.__layout(rect, test=False)

    def sizeHint(self):
        return self.minimumSize()

    def minimumSize(self):
        size = qt.QSize()
        for item in self.__items:
            size = size.expandedTo(item.minimumSize())
        margins = self.contentsMargins()
        return size + qt.QSize(
            margins.left() + margins.right(), margins.top() + margins.bottom()
        )

    def __layout(self, rect, test):
        """Place the items and return the required height.

        :param QRect rect: Area available to the layout.
        :param bool test: If True, only compute the height.
        :rtype: int
        """
        margins = self.contentsMargins()
        effective = rect.adjusted(
            margins.left(), margins.top(), -margins.right(), -margins.bottom()
        )
        x = effective.x()
        y = effective.y()
        rowHeight = 0
        for item in self.__items:
            hint = item.sizeHint()
            nextX = x + hint.width() + self.horizontalSpacing()
            if nextX - self.horizontalSpacing() > effective.right() and rowHeight > 0:
                x = effective.x()
                y = y + rowHeight + self.verticalSpacing()
                nextX = x + hint.width() + self.horizontalSpacing()
                rowHeight = 0
            if not test:
                item.setGeometry(qt.QRect(qt.QPoint(x, y), hint))
            x = nextX
            rowHeight = max(rowHeight, hint.height())
        return y + rowHeight - rect.y() + margins.bottom()


class _CompactSpinBoxMixin:
    """Allow a spin box to shrink below the width of its widest value.

    Qt reserves enough horizontal space in the minimum size of a spin box to
    display the widest value of the allowed range. Scientific ranges (for
    example ``+-20000`` with four decimals) therefore reserve far more space
    than the values that are typically entered, which inflates the minimum
    width of every panel containing such a spin box.

    This mixin caps the preferred width (``sizeHint``) at ``preferredChars``
    digits and the minimum width (``minimumSizeHint``) at ``visibleChars``
    digits, in both cases plus prefix and suffix. Values wider than that are
    still accepted and can be scrolled inside the editor, and the spin box
    still grows when the layout has room to spare. All widths are derived from
    the current font metrics and the current style, so the result follows the
    font size, the widget style and high DPI scaling.

    The horizontal size policy is relaxed to ``QSizePolicy.Preferred``, since
    the default ``QSizePolicy.Minimum`` of spin boxes would prevent a layout
    from using the reduced minimum width at all.
    """

    DEFAULT_VISIBLE_CHARS = 6
    DEFAULT_PREFERRED_CHARS = 9

    def __init__(self, parent=None, visibleChars=None, preferredChars=None):
        super().__init__(parent)
        self.__visibleChars = (
            self.DEFAULT_VISIBLE_CHARS if visibleChars is None else int(visibleChars)
        )
        self.__preferredChars = (
            self.DEFAULT_PREFERRED_CHARS
            if preferredChars is None
            else int(preferredChars)
        )
        policy = self.sizePolicy()
        policy.setHorizontalPolicy(qt.QSizePolicy.Preferred)
        self.setSizePolicy(policy)

    def getVisibleChars(self):
        """Return the number of value characters kept visible at minimum width.

        :rtype: int
        """
        return self.__visibleChars

    def setVisibleChars(self, visibleChars):
        """Set the number of value characters kept visible at minimum width.

        :param int visibleChars:
            Number of digits (excluding prefix and suffix) that must remain
            readable when the spin box is compressed to its minimum width.
        """
        visibleChars = int(visibleChars)
        if visibleChars != self.__visibleChars:
            self.__visibleChars = visibleChars
            self.updateGeometry()

    def getPreferredChars(self):
        """Return the number of value characters used for the preferred width.

        :rtype: int
        """
        return self.__preferredChars

    def setPreferredChars(self, preferredChars):
        """Set the number of value characters used for the preferred width.

        :param int preferredChars:
            Number of digits (excluding prefix and suffix) the spin box asks
            for when the layout has enough room.
        """
        preferredChars = int(preferredChars)
        if preferredChars != self.__preferredChars:
            self.__preferredChars = preferredChars
            self.updateGeometry()

    def __cappedWidth(self, hint, chars):
        """Return the width of `hint` reduced to `chars` value characters.

        The width Qt reserved for the widest value of the range is replaced by
        the width of `chars` digits, keeping every style dependent extra
        (frame, button box, spacing) that `hint` already contains.

        :param QSize hint: Size computed by the base spin box implementation.
        :param int chars: Number of value characters that must fit.
        :rtype: int
        """
        metrics = self.fontMetrics()

        def advance(text):
            return metrics.horizontalAdvance(text)

        widest = max(
            advance(self.prefix() + self.textFromValue(self.minimum()) + self.suffix()),
            advance(self.prefix() + self.textFromValue(self.maximum()) + self.suffix()),
        )
        # "0" is not necessarily the widest digit, use the widest one instead
        digit = max(advance(str(d)) for d in range(10))
        compact = advance(self.prefix() + self.suffix()) + chars * digit
        return hint.width() - max(0, widest - compact)

    def sizeHint(self):
        hint = super().sizeHint()
        chars = max(self.__preferredChars, self.__visibleChars)
        return qt.QSize(self.__cappedWidth(hint, chars), hint.height())

    def minimumSizeHint(self):
        hint = super().minimumSizeHint()
        return qt.QSize(self.__cappedWidth(hint, self.__visibleChars), hint.height())


class CompactDoubleSpinBox(_CompactSpinBoxMixin, qt.QDoubleSpinBox):
    """:class:`QDoubleSpinBox` which can shrink below its full range width.

    See :class:`_CompactSpinBoxMixin` for the sizing behaviour.
    """


class CompactSpinBox(_CompactSpinBoxMixin, qt.QSpinBox):
    """:class:`QSpinBox` which can shrink below its full range width.

    See :class:`_CompactSpinBoxMixin` for the sizing behaviour.
    """


class AspectRatioPixmapLabel(qt.QLabel):
    def __init__(self, parent=None):
        qt.QLabel.__init__(self, parent)
        self.setMinimumSize(1, 1)
        self.setScaledContents(False)
        self.pix = None

    def setPixmap(self, p):
        self.pix = p
        super().setPixmap(self.scaledPixmap())

    def scaledPixmap(self):
        return self.pix.scaled(
            self.size(), qt.Qt.KeepAspectRatio, qt.Qt.SmoothTransformation
        )

    def heightForWidth(self, width):
        if self.pix is None:
            return self.height()
        else:
            return int((self.pix.height() * width) / self.pix.width())

    def sizeHint(self):
        app = qt.QApplication.instance()
        screenGeometry = app.primaryScreen().availableGeometry()
        w = int(screenGeometry.width() / 3)
        w_s = self.width()
        return qt.QSize(max(w, w_s), self.heightForWidth(w))

    def resizeEvent(self, e):
        if self.pix is not None:
            super().setPixmap(self.scaledPixmap())


class DataRangeSlider(qt.QWidget):
    """Slider widget, with 4 buttons/icons and a line edit to provide
    a way of selecting a
    """

    sigValueChanged = qt.pyqtSignal(object)

    def __init__(self, parent=None, data=None, unit=None):
        qt.QWidget.__init__(self, parent)
        self._data = None

        # Use the font size as the icon size to avoid to create bigger buttons
        fontMetric = self.fontMetrics()
        iconSize = qt.QSize(fontMetric.height(), fontMetric.height())

        self.mainLayout = qt.QHBoxLayout(self)
        self.mainLayout.setContentsMargins(0, 0, 0, 0)
        self.mainLayout.setSpacing(0)

        self.slider = qt.QSlider()
        self.slider.setOrientation(qt.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)

        self.firstButton = qt.QPushButton(self)
        self.firstButton.setIcon(icons.getQIcon("first"))
        self.firstButton.setIconSize(iconSize)
        self.previousButton = qt.QPushButton(self)
        self.previousButton.setIcon(icons.getQIcon("previous"))
        self.previousButton.setIconSize(iconSize)
        self._lineEdit = qt.QLineEdit(self)

        self._label = qt.QLabel(self)
        self.nextButton = qt.QPushButton(self)
        self.nextButton.setIcon(icons.getQIcon("next"))
        self.nextButton.setIconSize(iconSize)
        self.lastButton = qt.QPushButton(self)
        self.lastButton.setIcon(icons.getQIcon("last"))
        self.lastButton.setIconSize(iconSize)

        self.mainLayout.addWidget(self.slider)
        self.mainLayout.addWidget(self.firstButton)
        self.mainLayout.addWidget(self.previousButton)
        self.mainLayout.addWidget(self._lineEdit)
        self.mainLayout.addWidget(self._label)
        self.mainLayout.addWidget(self.nextButton)
        self.mainLayout.addWidget(self.lastButton)

        if data is None:
            first = qt.QSlider().minimum()
            last = qt.QSlider().maximum()
        else:
            first, last = data[0], data[-1]
            self.slider.setMaximum(data.size)

        self._lineEdit.setFixedWidth(
            self._lineEdit.fontMetrics().boundingRect(f"{last:.5f}").width()
        )
        validator = AxisValidator(self._lineEdit)
        self._lineEdit.setValidator(validator)
        txt = self._lineEdit.validator().fixup(first)
        self._lineEdit.setText(txt)
        if unit is None:
            self._label.setText(f"of {last:f}")
        else:
            self._label.setText(f"of {last:.5f} {unit}")

        self._lineTxt = self._lineEdit.text()

        """0-based index"""

        self.firstButton.clicked.connect(self._firstClicked)
        self.previousButton.clicked.connect(self._previousClicked)
        self.nextButton.clicked.connect(self._nextClicked)
        self.lastButton.clicked.connect(self._lastClicked)
        self._lineEdit.editingFinished.connect(self._textChangedSlot)
        self.slider.valueChanged.connect(self._sliderChangedSlot)

    def setAxis(self, data, unit=None):
        self._data = np.copy(data)
        first, last = self._data[0], self._data[-1]
        self.slider.setMaximum(self._data.size - 1)

        self._lineEdit.validator().setData(self._data)
        txt = self._lineEdit.validator().fixup(first)
        self._lineEdit.setText(txt)
        self._lineTxt = self._lineEdit.text()

        if unit is None:
            self._label.setText(f"of {last:.5f}")
        else:
            self._label.setText(f"of {last:.5f} {unit}")

        self._textChangedSlot()

    def lineEdit(self):
        """Returns the line edit provided by this widget.

        :rtype: qt.QLineEdit
        """
        return self._lineEdit

    def limitWidget(self):
        """Returns the widget displaying axes limits.

        :rtype: qt.QLabel
        """
        return self._label

    def _firstClicked(self):
        """Select first/lowest frame number"""
        if self._data is None:
            return
        self.setIndex(self.getRange()[0])

    def _previousClicked(self):
        """Select previous frame number"""
        if self._data is None:
            return
        self.setIndex(self.getIndex() - 1)

    def _nextClicked(self):
        """Select next frame number"""
        if self._data is None:
            return
        self.setIndex(self.getIndex() + 1)

    def _lastClicked(self):
        """Select last/highest frame number"""
        if self._data is None:
            return
        self.setIndex(self.getRange()[1] - 1)

    def _textChangedSlot(self):
        """Select frame number typed in the line edit widget"""
        if self._data is None:
            return
        txt = self._lineEdit.text()
        if txt == self._lineTxt:
            return
        idx = self.getIndex()
        with qt.QSignalBlocker(self.slider):
            self.slider.setValue(idx)
        ddict = {
            "event": "indexChanged",
            "oldtxt": self._lineTxt,
            "newtxt": txt,
            "idx": idx,
            "value": self._data[idx],
            "id": id(self),
        }
        self._lineTxt = txt
        self.sigValueChanged.emit(ddict)

    def _sliderChangedSlot(self):
        if self._data is None:
            return
        idx = self.slider.value()
        value = self._data[idx]

        txt = self._lineEdit.validator().fixup(value)
        if txt == self._lineTxt:
            return
        self._lineEdit.setText(txt)
        ddict = {
            "event": "indexChanged",
            "oldtxt": self._lineTxt,
            "newtxt": txt,
            "idx": idx,
            "value": self._data[idx],
            "id": id(self),
        }
        self._lineTxt = txt
        self.sigValueChanged.emit(ddict)

    def getRange(self):
        """ """
        if self._data is None:
            return 0, 0
        else:
            return 0, self._data.size - 1

    def getIndex(self):
        if self._data is None:
            raise ValueError("Data is not set.")
        return self._lineEdit.validator().getIndex(self._lineEdit.text())

    def setValue(self, value):
        if self._lineEdit.validator().validate(value, 0) == qt.QValidator.Invalid:
            raise ValueError(f"Invalid value: {value}")
        txt = self._lineEdit.validator().fixup(value)
        self._lineEdit.setText(txt)
        self._textChangedSlot()

    def getValue(self):
        """Return current frame index"""
        if self._data is None:
            raise ValueError("Data is not set.")
        idx = self.getIndex()
        value = self._data[idx]
        return value

    def setIndex(self, idx):
        """Set 0-based frame index

        Value is clipped to current range.
        """
        if self._data is None:
            return
        bottom, top = self.getRange()
        idx = int(idx)

        if idx < bottom:
            idx = bottom
        elif idx > top:
            idx = top

        value = self._data[idx]
        txt = self._lineEdit.validator().fixup(value)
        self._lineEdit.setText(txt)
        self._textChangedSlot()


class AxisValidator(
    qt.QValidator
):  # each widget needs its own validator due to API issue.
    def __init__(self, lineedit, parent=None, data=None):
        qt.QValidator.__init__(self, parent)
        self.lineedit = lineedit
        self.setData(data)

    def validate(self, input_, pos):
        try:
            val = float(input_)
        except Exception:
            return qt.QValidator.Invalid, input_, pos

        if self._data is not None:
            if self._min <= val <= self._max:
                idx = np.argmin(np.abs(self._data - val))
                exactval = self._data[idx]
                exactrepr = f"{exactval:.5f}"
                if f"{val:.5f}" == exactrepr:
                    return qt.QValidator.Acceptable, input_, pos
            return qt.QValidator.Intermediate, input_, pos
        else:
            return qt.QValidator.Acceptable, input_, pos

    def getIndex(self, input_):
        if self.validate(input_, 0) == qt.QValidator.Invalid:
            raise ValueError(f"Invalid input: {input_}")
        val = float(input_)
        idx = np.argmin(np.abs(self._data - val))
        return idx

    def getData(self):
        return self._data

    def setData(self, data):
        if data is not None:
            self._data = np.copy(data)
            self._max = np.amax(self._data)
            self._min = np.amin(self._data)
        else:
            self._data = None

    def fixup(self, input_):  # broken PyQt API - have to use workaround
        """Overrides LineEdit.text!!!!"""
        val = float(input_)
        if self._data is not None:
            idx = np.argmin(np.abs(self._data - val))
            exactval = self._data[idx]
        else:
            exactval = val
        exactrepr = f"{exactval:.5f}"
        # self.lineedit.setText(exactrepr) # API workaround - this should not be required!  # noqa: E501
        return exactrepr

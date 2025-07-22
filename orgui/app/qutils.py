# -*- coding: utf-8 -*-
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

from silx.gui import qt
from silx.gui import icons
import numpy as np


def messagebox_detailed_message(parent, title, text, detailed_text, icon, buttons=qt.QMessageBox.Ok):
    diag = qt.QMessageBox(icon, title, text, buttons, parent)
    diag.setDetailedText(detailed_text)
    return diag.exec()

def critical_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Critical, buttons=buttons)
    
def warning_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Warning, buttons=buttons)

def information_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Information, buttons=buttons)


class AspectRatioPixmapLabel(qt.QLabel):
    def __init__(self, parent=None):
        qt.QLabel.__init__(self, parent)
        self.setMinimumSize(1,1)
        self.setScaledContents(False)
        self.pix = None
        
    def setPixmap(self, p):
        self.pix = p
        super().setPixmap(self.scaledPixmap())
        
    def scaledPixmap(self):
        return self.pix.scaled(self.size(), qt.Qt.KeepAspectRatio, qt.Qt.SmoothTransformation)

    def heightForWidth(self, width):
        if self.pix is None:
            return self.height()
        else:
            return int(( self.pix.height()* width) /self.pix.width())

    def sizeHint(self):
        app = qt.QApplication.instance()
        screenGeometry = app.primaryScreen().availableGeometry()
        w = int(screenGeometry.width()/3)
        w_s = self.width()
        return qt.QSize( max(w, w_s), self.heightForWidth(w))
        
    def resizeEvent(self,e):
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
            self._lineEdit.fontMetrics().boundingRect("%.5f" % last).width()
        )
        validator = AxisValidator(self._lineEdit)
        self._lineEdit.setValidator(validator)
        txt = self._lineEdit.validator().fixup(first)
        self._lineEdit.setText(txt)
        if unit is None:
            self._label.setText("of %f" % last)
        else:
            self._label.setText("of %.5f %s" % (last, unit))
        
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
        self.slider.setMaximum(self._data.size-1)
        
        self._lineEdit.validator().setData(self._data)
        txt = self._lineEdit.validator().fixup(first)
        self._lineEdit.setText(txt)
        self._lineTxt = self._lineEdit.text()
        
        if unit is None:
            self._label.setText("of %.5f" % last)
        else:
            self._label.setText("of %.5f %s" % (last, unit))
            
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
        self.setIndex(self.getRange()[1] -1)

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
            "idx" : idx,
            "value" : self._data[idx],
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
            "idx" : idx,
            "value" : self._data[idx],
            "id": id(self),
        }
        self._lineTxt = txt
        self.sigValueChanged.emit(ddict)
        

    def getRange(self):
        """
        """
        if self._data is None:
            return 0, 0
        else:
            return 0, self._data.size -1

    def getIndex(self):
        if self._data is None:
            raise ValueError("Data is not set.")
        return self._lineEdit.validator().getIndex(self._lineEdit.text())
        
    def setValue(self, value):
        if self._lineEdit.validator().validate(value, 0) == qt.QValidator.Invalid:
            raise ValueError("Invalid value: %s" % value)
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
        
        
        
        

class AxisValidator(qt.QValidator): # each widget needs its own validator due to API issue.
    def __init__(self, lineedit, parent=None, data=None):
        qt.QValidator.__init__(self, parent)
        self.lineedit = lineedit
        self.setData(data)
            
    def validate(self, input_, pos):
        try:
            val = float(input_)
        except:
            return qt.QValidator.Invalid, input_, pos
        
        if self._data is not None:
            if self._min <= val <= self._max:
                idx = np.argmin(np.abs(self._data - val))
                exactval = self._data[idx]
                exactrepr = "%.5f" % exactval
                if "%.5f" % val == exactrepr:
                    return qt.QValidator.Acceptable, input_, pos
            return qt.QValidator.Intermediate, input_, pos
        else:
            return qt.QValidator.Acceptable, input_, pos
            
    def getIndex(self, input_):
        if self.validate(input_, 0) == qt.QValidator.Invalid:
            raise ValueError("Invalid input: %s" % input_)
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
                
    def fixup(self, input_): # broken PyQt API - have to use workaround
        """ Overrides LineEdit.text!!!!
        
        """
        val = float(input_)
        if self._data is not None:
            idx = np.argmin(np.abs(self._data - val))
            exactval = self._data[idx]
        else:
            exactval = val
        exactrepr = "%.5f" % exactval
        #self.lineedit.setText(exactrepr) # API workaround - this should not be required!
        return exactrepr
        
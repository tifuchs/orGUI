from silx.gui import qt
from . import qutils
import numpy as np
from contextlib import contextmanager

@contextmanager
def blockSignals(qobjects):
    try:
        for obj in qobjects:
            obj.blockSignals(True)
        yield
        for obj in qobjects:
            obj.blockSignals(False)
    except TypeError:
        qobject.blockSignals(True)
        yield
        qobject.blockSignals(False)

class HKLDialog(qt.QDialog):
    sigHKLchanged = qt.pyqtSignal(np.ndarray)
    
    def __init__(self, msg='', title='' ,parent=None):
        qt.QDialog.__init__(self, parent)
        self.savedParams = np.array([0.,0.,0.])
        layout = qt.QVBoxLayout()
        self.setWindowTitle(title)
        
        hkl_layout = qt.QHBoxLayout()
        
        self.hkl_editors = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.hkl_editors]
        [h.setDecimals(3) for h in self.hkl_editors]
        [h.setValue(0.) for h in self.hkl_editors]
        for editor, lbll in zip(self.hkl_editors, ["H:", "K:", "L:"]):
            hkl_layout.addWidget(qt.QLabel(lbll))
            hkl_layout.addWidget(editor)
            editor.valueChanged.connect(self._onAnyHKLChanged)

        
        self.label = qt.QLabel()
        self.label.setText(msg)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel | qt.QDialogButtonBox.Reset)
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.onCancel)
        buttons.button(qt.QDialogButtonBox.Reset).clicked.connect(self.resetParameters)
        
        layout.addWidget(self.label)
        layout.addLayout(hkl_layout)
        layout.addWidget(buttons)
        self.setLayout(layout)
        

    def _onAnyHKLChanged(self):
        self.sigHKLchanged.emit(self.get_hkl())
        
    def get_hkl(self):
        return np.array([h.value() for h in self.hkl_editors])
        
    def set_hkl(self, hkl):
        assert len(hkl) == 3
        with blockSignals(self.hkl_editors):
            for h, e in zip(hkl, self.hkl_editors):
                e.setValue(h)
        self.sigHKLchanged.emit(self.get_hkl())
        
    def showEvent(self, event):
        if event.spontaneous():
            super().showEvent(event)
        else:
            self.savedParams = self.get_hkl()
            super().showEvent(event)
            
    #def hideEvent(self, event):
    #    self.sigHide.emit()
    #    super().hideEvent(event)
            
    def resetParameters(self):
        self.set_hkl(self.savedParams)
            
    def onCancel(self):
        self.resetParameters()
        self.reject()
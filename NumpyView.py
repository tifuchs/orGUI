from PyMca5.PyMcaGui import PyMcaQt as qt
from PyMca5.PyMcaGui.misc import NumpyArrayTableWidget, NumpyArrayTableModel
from PyMca5.PyMcaGui.misc import FrameBrowser
from PyMca5.PyMcaGui.misc.NumpyArrayTableView import HorizontalHeader, VerticalHeader
import sys



class NumpyArrayEditTableModel(NumpyArrayTableModel.NumpyArrayTableModel):
    
    def flags(self,index):
        fl = super(NumpyArrayEditTableModel,self).flags(index)
        return fl | qt.Qt.ItemIsEditable
        
    
    def _data1D(self, index, role=qt.Qt.DisplayRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                # row = 0
                col = index.column()
                return NumpyArrayTableModel.MyQVariant(self._format % self._array[col])
        return NumpyArrayTableModel.MyQVariant()
    
    def _data2D(self, index, role=qt.Qt.DisplayRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                return NumpyArrayTableModel.MyQVariant(self._format % self._array[row, col])
        return NumpyArrayTableModel.MyQVariant()
    
    def _dataND(self, index, role=qt.Qt.DisplayRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                actualSelection = tuple(self._index + [row, col])
                return NumpyArrayTableModel.MyQVariant(self._format % self._array[actualSelection])
        return NumpyArrayTableModel.MyQVariant()
    
    def _data3DIndex0(self, index, role=qt.Qt.DisplayRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                return NumpyArrayTableModel.MyQVariant(self._format % self._array[self._index,
                                                                  row,
                                                                  col])
        return NumpyArrayTableModel.MyQVariant()

    def _data3DIndex1(self, index, role=qt.Qt.DisplayRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                return NumpyArrayTableModel.MyQVariant(self._format % self._array[row,
                                                                  self._index,
                                                                  col])
        return NumpyArrayTableModel.MyQVariant()

    def _data3DIndex2(self, index, role=qt.Qt.DisplayRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                return NumpyArrayTableModel.MyQVariant(self._format % self._array[row,
                                                                  col,
                                                                  self._index])
        return NumpyArrayTableModel.MyQVariant()
    
    def setData(self,index, value,role=qt.Qt.EditRole):
        if (self._setData(index, value,role)):
            self.dataChanged.emit(self.index(0,0),self.index(self.rowCount()-1,self.columnCount()-1))
            return True
        else:
            return False
        
        
    def _setData1D(self, index, value, role=qt.Qt.EditRole):
        if index.isValid():
            if role == qt.Qt.EditRole:
                # row = 0
                col = index.column()
                self._array[col] = float(value)
                return True
        return False
                
    def _setData2D(self, index, value, role=qt.Qt.EditRole):
        if index.isValid():
            if role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                self._array[row, col] = float(value)
                return True
        return False
    
    def _setDataND(self, index, value, role=qt.Qt.EditRole):
        if index.isValid():
            if role == qt.Qt.DisplayRole or role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                actualSelection = tuple(self._index + [row, col])
                self._array[actualSelection] = float(value)
                return True
        return False
    
    
    def _setData3DIndex0(self, index, value, role=qt.Qt.EditRole):
        if index.isValid():
            if role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                self._array[self._index,row,col] = float(value)
                return True
        return False

    def _setData3DIndex1(self, index, value, role=qt.Qt.EditRole):
        if index.isValid():
            if role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                self._array[row,self._index,col] = float(value)
                return True
        return False

    def _setData3DIndex2(self, index, value, role=qt.Qt.EditRole):
        if index.isValid():
            if role == qt.Qt.EditRole:
                row = index.row()
                col = index.column()
                self._array[row,col,self._index] = float(value)
                return True
        return False
        
    def assignDataFunction(self, dimension):
        shape = self._array.shape
        if len(shape) == 2:
            self._rowCount = self._rowCount2D
            self._columnCount = self._columnCount2D
            self._data = self._data2D
            self._setData = self._setData2D
        elif len(shape) == 1:
            self._rowCount = self._rowCount1D
            self._columnCount = self._columnCount1D
            self._data = self._data1D
            self._setData = self._setData1D
        elif len(shape) > 3:
            # only C order array of images supported
            self._rowCount = self._rowCountND
            self._columnCount = self._columnCountND
            self._data = self._dataND
            self._setData = self._setDataND
        else:
            if dimension == 1:
                self._rowCount = self._rowCount3DIndex1
                self._columnCount = self._columnCount3DIndex1
                self._data = self._data3DIndex1
                self._setData = self._setData3DIndex1
            elif dimension == 2:
                self._rowCount = self._rowCount3DIndex2
                self._columnCount = self._columnCount3DIndex2
                self._data = self._data3DIndex1
                self._setData = self._setData3DIndex1
            else:
                self._rowCount = self._rowCount3DIndex0
                self._columnCount = self._columnCount3DIndex0
                self._data = self._data3DIndex0
                self._setData = self._setData3DIndex0
            self._dimension = dimension
    
    def data(self, index, role=qt.Qt.DisplayRole):
        return self._data(index, role)
    
    def getArrayData(self):
        return self._array


class NumpyArrayEditTableView(qt.QTableView):
    def __init__(self, parent=None):
        qt.QTableView.__init__(self, parent)
        self._model = NumpyArrayEditTableModel(self)
        self.setModel(self._model)
        self._horizontalHeaderModel = HorizontalHeader(self._model)
        self._verticalHeaderModel = VerticalHeader(self._model)
        self.horizontalHeader().setModel(self._horizontalHeaderModel)
        self.verticalHeader().setModel(self._verticalHeaderModel)

    def setArrayData(self, data):
        t = "%s" % data.dtype
        if '|' in t:
            fmt = "%s"
        else:
            fmt = "%g"
        self._model.setFormat(fmt)
        self._model.setArrayData(data)
        #some linux distributions need this call
        self.setModel(self._model)
        if sys.platform not in ['win32']:
            self._horizontalHeaderModel = HorizontalHeader(self._model)
            self._verticalHeaderModel = VerticalHeader(self._model)
        self.horizontalHeader().setModel(self._horizontalHeaderModel)
        self.verticalHeader().setModel(self._verticalHeaderModel)

    def setCurrentArrayIndex(self, index):
        return self._model.setCurrentArrayIndex(index)
    
    def getArrayData(self):
        return self._model.getArrayData()
    
class BrowserContainer(qt.QWidget):
    def __init__(self, parent=None):
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setContentsMargins(0, 0, 0, 0)
        self.mainLayout.setSpacing(0)

class NumpyArrayEditTableWidget(qt.QTableWidget):
    def __init__(self, parent=None):
        qt.QTableWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setContentsMargins(0, 0, 0, 0)
        self.mainLayout.setSpacing(0)
        self.browserContainer = BrowserContainer(self)
        self._widgetList = []
        for i in range(4):
            browser = FrameBrowser.HorizontalSliderWithBrowser(self.browserContainer)
            self.browserContainer.mainLayout.addWidget(browser)
            self._widgetList.append(browser)
            browser.valueChanged.connect(self.browserSlot)
            if i == 0:
                browser.setEnabled(False)
                browser.hide()
        self.view = NumpyArrayEditTableView(self)
        self.mainLayout.addWidget(self.browserContainer)
        self.mainLayout.addWidget(self.view)

    def setArrayData(self, data):
        self._array = data
        nWidgets = len(self._widgetList)
        nDimensions = len(self._array.shape) 
        if nWidgets > (nDimensions - 2):
            for i in range((nDimensions - 2), nWidgets):
                browser = self._widgetList[i]
                self._widgetList[i].setEnabled(False)
                self._widgetList[i].hide()
        else:
            for i in range(nWidgets, nDimensions - 2):
                browser = FrameBrowser.HorizontalSliderWithBrowser(self.browserContainer)
                self.browserContainer.mainLayout.addWidget(browser)
                self._widgetList.append(browser)
                browser.valueChanged.connect(self.browserSlot)
                browser.setEnabled(False)
                browser.hide()
        for i in range(nWidgets):
            browser = self._widgetList[i]
            if (i + 2 ) < nDimensions:
                browser.setEnabled(True)
                if browser.isHidden():
                    browser.show()
                browser.setRange(1, self._array.shape[i])
            else:
                browser.setEnabled(False)
                browser.hide()
        self.view.setArrayData(self._array)
        
    def getArrayData(self):
        return self.view.getArrayData()

    def browserSlot(self, value):
        if len(self._array.shape) == 3:
            self.view.setCurrentArrayIndex(value - 1)
            self.view.reset()
        else:
            index = []
            for browser in self._widgetList:
                if browser.isEnabled():
                    index.append(browser.value() - 1)
            self.view.setCurrentArrayIndex(index)
            self.view.reset()
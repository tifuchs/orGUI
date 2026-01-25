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

import numpy as np
import os
import traceback
from silx.gui import qt
from silx.gui.data import ArrayTableWidget, ArrayTableModel
from silx.gui.widgets.TableWidget import TableView
from silx.gui.dialog.DataFileDialog import DataFileDialog
from silx.gui import icons

from silx.io.utils import save1D, savespec, NEXUS_HDF5_EXT
from silx.io.nxdata import save_NXdata

from orgui import resources

import logging
logger = logging.getLogger(__name__)


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


class DeleteRowAction(qt.QAction):
    """QAction to delete rows in a ArrayTableWidget

    :param table: :class:`QTableView` to which this action belongs.
    """
    def __init__(self, table, arrayTableWidget):
        if not isinstance(table, qt.QTableView):
            raise ValueError('DeleteRowAction must be initialised ' +
                             'with a QTableWidget.')
        super(DeleteRowAction, self).__init__(table)
        self.table = table
        self.arrayTableWidget = arrayTableWidget
        self.setText("Delete row")
        self.setShortcut(qt.QKeySequence.Delete)
        self.setShortcutContext(qt.Qt.WidgetShortcut)
        self.setToolTip("Delete row of the array.")
        self.triggered.connect(self.deleteSelectedRows)

    def deleteSelectedRows(self):
        """Paste text from clipboard into the table.

        :return: *True* in case of success, *False* if failed.
        """
        selected_idx = self.table.selectedIndexes()
        if len(selected_idx) < 1:
            #msgBox = qt.QMessageBox(parent=self.table)
            #msgBox.setText("A single cell must be selected to delete data")
            #msgBox.exec()
            return False
        
        
        data_model = self.table.model()
        data = data_model.getData()
        
        selected_row = np.unique([idx.row() for idx in selected_idx])
        mask = np.ones(data.shape[0], dtype=np.bool_)
        mask[selected_row] = False

        self.arrayTableWidget.updateArrayData(data[mask])
        self.arrayTableWidget.sigRowsDeleted.emit(selected_row)
        
        
class AddRowAction(qt.QAction):
    """QAction to add a row in a ArrayTableWidget

    :param table: :class:`QTableView` to which this action belongs.
    """
    def __init__(self, table, arrayTableWidget, fill_value=-1):
        if not isinstance(table, qt.QTableView):
            raise ValueError('AddRowAction must be initialised ' +
                             'with a QTableWidget.')
        super(AddRowAction, self).__init__(table)
        self.fill_value = fill_value
        self.table = table
        self.arrayTableWidget = arrayTableWidget
        self.setText("Add row")
        self.setShortcut(qt.QKeySequence.New)
        self.setShortcutContext(qt.Qt.WidgetShortcut)
        self.setToolTip("Add row to the array.")
        self.triggered.connect(self.addRow)

    def addRow(self):
        """Paste text from clipboard into the table.

        :return: *True* in case of success, *False* if failed.
        """
        selected_idx = self.table.selectedIndexes()
        if len(selected_idx) > 1:
            # msgBox = qt.QMessageBox(parent=self.table)
            # msgBox.setText("A single cell must be selected to add a row")
            # msgBox.exec()
            logger.error("Cannot add row to table: A single cell must be selected to add a row", 
                 extra={'title' : 'Cannot add row to table',
                        'show_dialog' : True,
                        'parent' : self.table})
            return False
        
        data_model = self.table.model()
        data = data_model.getData()
        
        if len(selected_idx) == 0:
            selected_row = data.shape[0]
        else:
            selected_row = selected_idx[0].row()

        newdata = np.insert(data, selected_row, np.full(data.shape[-1],self.fill_value) ,0)
        self.arrayTableWidget.updateArrayData(newdata)
        self.arrayTableWidget.sigRowAdded.emit(selected_row)
        
        
class ArrayTableHeaderModel(ArrayTableModel.ArrayTableModel):
    def headerData(self, section, orientation, role=qt.Qt.DisplayRole):
        """QAbstractTableModel method
        Return the 0-based row or column index, for display in the
        horizontal and vertical headers"""
        data = super().headerData(section, orientation, role)

        if role == qt.Qt.DisplayRole and orientation == qt.Qt.Horizontal and\
                   is_number(data):
            if hasattr(self, "header"):
                return self.header[section] if self.header is not None else data
            else:
                return data
        else:
            return data

class ArrayEditWidget(ArrayTableWidget.ArrayTableWidget):
    """Widget to display and edit an numpy array. It currently only supports 
    1d and 2d arrays. 
    
    It allows editong of entries and also the adding and removing of rows in the
    array. I.e. the shape of the numpy array can be changed.
    This however means that the array reference is changed and this widget cannot be used
    like the :class:`ArrayTableWidget`, where the reference can remain unchanged.
    Instead, :meth:`~ArrayEditWidget.getData` must be called every time!
    
    A ToolBar will be created above the array view, which includes action buttons for 
    saving and loading of the array.
    
    """
    sigRowAdded = qt.pyqtSignal(int)
    sigRowsDeleted = qt.pyqtSignal(np.ndarray)
    sigDataLoaded = qt.pyqtSignal()
    
    def __init__(self, saveact: bool = True, openact: int = -1, rowActions=True, parent = None):
        """Creates a new :class:`ArrayEditWidget`. 
        
        If ``openact`` is set, only arrays which have the same size
        of this particular axis will be opened by the ToolBar buttons. 
        Set it to ``-1`` to diable the filtering. 
        
        :param saveact: Display the save array buttons.
        :type saveact: bool
        :param openact: Display the load array buttons. Set it to ``-1`` to diable the buttons.
        :type openact: int
        :param parent: parent QWidget
        :param labels: list of labels for each dimension of the array
        """
        qt.QWidget.__init__(self, parent)
        self.mainLayout = qt.QVBoxLayout(self)
        self.mainLayout.setContentsMargins(0, 0, 0, 0)
        self.mainLayout.setSpacing(0)

        self.browserContainer = qt.QWidget(self)
        self.browserLayout = qt.QGridLayout(self.browserContainer)
        self.browserLayout.setContentsMargins(0, 0, 0, 0)
        self.browserLayout.setSpacing(0)

        self._dimensionLabelsText = []
        """List of text labels sorted in the increasing order of the dimension
        they apply to."""
        self._browserLabels = []
        """List of QLabel widgets."""
        self._browserWidgets = []
        """List of HorizontalSliderWithBrowser widgets."""

        self.axesSelector = ArrayTableWidget.AxesSelector(self)
        self.axesSelector.setVisible(False)

        self.view = TableView(self)
        
        self.toolbar = qt.QToolBar("Array modifier", self)
        self.filedialogdir = os.getcwd()
        if openact >= 0:
            self.enableOpenAction()
            
        if saveact:
            self.enableSaveAction()
        
        self.mainLayout.addWidget(self.toolbar)
        #self.mainLayout.addWidget(self.browserContainer)
        #self.mainLayout.addWidget(self.axesSelector)
        self.mainLayout.addWidget(self.view)

        self.model = ArrayTableHeaderModel(self)
        self.view.setModel(self.model)
        if rowActions:
            self.view.deleteRowAction = DeleteRowAction(self.view, self)
            self.view.addAction(self.view.deleteRowAction)
            
            self.view.addRowAction = AddRowAction(self.view, self)
            self.view.addAction(self.view.addRowAction)
        
        self.retain_axis = openact
    
    def getData(self, copy=True):
        return np.atleast_1d(np.squeeze(super().getData(copy)))
        
    def enableOpenAction(self):
        if not hasattr(self, 'openAct'):
            self.openAct = self.toolbar.addAction(icons.getQIcon('document-open'), "open ascii file")
            self.openAct.triggered.connect(lambda x: self.openLoadtxt(self.retain_axis))
        if not hasattr(self, 'openNXAct'):
            self.openNXAct = self.toolbar.addAction(resources.getQicon('document-nx-open'), "open NEXUS-like")
            self.openNXAct.triggered.connect(lambda x: self.openNXlike(self.retain_axis))
            
    def enableSaveAction(self):
        if not hasattr(self, 'saveAct'):
            self.saveAct = self.toolbar.addAction(icons.getQIcon('document-save'), "save as file")
            self.saveAct.triggered.connect(self.savetxt)
            
    def openLoadtxt(self, retain_axis=1):
        """Displays a dialog to open an array using the `:meth:~numpy.loadtxt` function.
        
        If ``retain_axis`` is set, only arrays which have the same size
        of this particular axis will be accepted. 
        Set it to ``-1`` to diable the filtering. The default will only accept 
        arrays with the same number of columns as the original array. 
        
        :param retain_axis: Axis size filtering. Set to -1 to disable.
        :type retain_axis: int
        """
        fileTypeDict = {'dat Files (*.dat)': '.dat', 'txt Files (*.txt)': '.txt', 'All files (*)': '', }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
        #print(retain_axis)    
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open file",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.filedialogdir = os.path.splitext(filename)[0]
        try:
            array_data = np.squeeze(np.loadtxt(filename))
            if retain_axis >= 0:
                if len(array_data.shape) == 1:
                    if len(self.getData().shape) != 1:
                        raise ValueError("Array must be a 1d array. Is: %s d" % len(self.getData().shape))
                    else:
                        self.updateArrayData(array_data)
                        self.sigDataLoaded.emit()
                        return True
                else:
                    if array_data.shape[retain_axis] != self.getData().shape[retain_axis]:
                        raise ValueError("The number of data columns must match. Required: %s, Is: %s" % (self.getData().shape[retain_axis], array_data.shape[retain_axis]))
                    else:
                        self.updateArrayData(array_data)
                        self.sigDataLoaded.emit()
                        return True
            else:
                self.updateArrayData(array_data, header=None, labels=None)
                self.sigDataLoaded.emit()
                return True
            
        except Exception as e:
            logger.exception("Error during loading of array", 
                             extra={'title' : 'Error during loading of array',
                                    'show_dialog' : True,
                                    'description' : str(e),
                                    'parent' : self})
            # qt.QMessageBox.critical(self,"Error during loading of array","Error during loading of array.\n%s" % traceback.format_exc())
            
    def openNXlike(self, retain_axis=1):
        """Displays a dialog to open an array from any NEXUS like data source.
        
        If ``retain_axis`` is set, only arrays which have the same size
        of this particular axis will be accepted. 
        Set it to ``-1`` to diable the filtering. The default will only accept 
        arrays with the same number of columns as the original array. 
        
        :param retain_axis: Axis size filtering. Set to -1 to disable.
        :type retain_axis: int
        """
        dialog = DataFileDialog(self)
        dialog.setFilterMode(DataFileDialog.FilterMode.ExistingDataset)
        
        def customFilter(obj):
            if len(obj.shape) == 0:
                return False
            if retain_axis >= 0:
                if len(obj.shape) == 1:
                    if len(self.getData().shape) != 1:
                        return False
                    else:
                        return True
                else:
                    if obj.shape[retain_axis] != self.getData().shape[retain_axis]:
                        return False
                    else:
                        return True
            else:
                return True
        dialog.setFilterCallback(customFilter)
        result = dialog.exec()
        if result:
            data = dialog.selectedData()
            if retain_axis >= 0:
                self.updateArrayData(data)
            else:
                self.updateArrayData(data, header=None, labels=None)
            self.sigDataLoaded.emit()
                
    def savetxt(self):
        fmt = "%.7g"
        csvdelim=";"
        fileTypeDictSave1D = {"Plain ascii file (*.dat)" : "dat", "CSV file (*.csv)" : "csv",  "NumPy format (*.npy)" : "ndarray"}
        
        fileTypeDictSPEC = {'SPEC file (*.spec)': 'spec'}
        
        fileTypeDict = {**fileTypeDictSave1D, **fileTypeDictSPEC}
        
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Save array to file",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.filedialogdir = os.path.splitext(filename)[0]

        if filetype in fileTypeDictSave1D:
            if self.header is not None:
                header = " ".join(self.header) 
            elif isinstance(self.header, str):
                header = self.header
            else:
                header = " "
            fileext = fileTypeDictSave1D[filetype]
            data = self.getData()
            if fileext == "dat":
                np.savetxt(filename, data, header=header, fmt=fmt)
            elif fileext == "csv":
                np.savetxt(filename, data, header=header, fmt=fmt, delimiter=csvdelim)
            elif fileext == "ndarray":
                np.save(filename, data)
            else:
                raise Exception("No supported file type %s" % fileext)
        
    def updateArrayData(self, data, **kwargs):
        editable = kwargs.get('editable', self.model._editable)
        labels = kwargs.get('labels', self._dimensionLabelsText)
        header = kwargs.get('header', self.header)
        self.setArrayData(data, labels=labels, editable=editable, header=self.header)
    
    def setArrayData(self, data, labels=None, copy=True, editable=False, header=None):
        if len(data.shape) == 1:
            data = data[:,np.newaxis]
        elif len(data.shape) > 2:
            raise ValueError("Only supports 1D or 2D arrays")
        if header is not None:
            if len(header) == data.shape[1]:
                self.model.header = header
                self.header = header
            else:
                raise ValueError("Header size does not match data shape.")
        else:
            self.model.header = None
            self.header = None
        super().setArrayData(data, labels, copy, editable)
        self.view.resizeColumnsToContents()
        self.view.resizeRowsToContents()
        

class ArrayTableDialog(qt.QDialog):
    
    def __init__(self, saveact: bool = True, openact: int = -1, parent = None):
        qt.QDialog.__init__(self, parent)

        self.arrayWidget = ArrayEditWidget(saveact, openact, self)
        
        self.setArrayData = self.arrayWidget.setArrayData
        self.updateArrayData = self.arrayWidget.updateArrayData
        self.getData = self.arrayWidget.getData

        layout = qt.QVBoxLayout(self)
        
        layout.addWidget(self.arrayWidget)
        
        self.setLayout(layout)
        
if __name__ == "__main__":
    
    
    app = qt.QApplication([])
    
    diag = ArrayTableDialog(True, 1)
    
    array = np.arange(10*1).reshape((10,1))#[:,np.newaxis]
    #array = np.arange(10)
    diag.setArrayData(array,editable=True, header= ['h'])
    
    diag.show()
    app.exec()
    

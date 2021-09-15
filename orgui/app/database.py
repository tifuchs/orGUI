# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2020-2021 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2021, Timo Fuchs, Olaf Magnussen all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import sys
import numpy as np
import silx.gui.hdf5
from silx.gui import icons
from silx.gui.data import DataViewerFrame
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
from silx.gui.data import DataViewerFrame
from silx.gui import qt
import tempfile
import silx.io.h5py_utils
from silx.io import nxdata
import silx.io.utils
import h5py

import datetime
import os
import traceback

from silx.io.dictdump import dicttonx ,nxtodict


class DataBase(qt.QMainWindow):
    
    def __init__(self, plot, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.nxfile = None
        self.filedialogdir = os.getcwd()
        self.plot = plot
        
        self.view = silx.gui.hdf5.Hdf5TreeView()
        self.view.setSortingEnabled(True)
        
        #self.view.addContextMenuCallback(self.nexus_treeview_callback)
        #self.hdf5model = self.view.findHdf5TreeModel()
        
        customModel = Hdf5TreeModel
        customModel.NAME_COLUMN = 0
        customModel.DESCRIPTION_COLUMN = 1
        customModel.VALUE_COLUMN = 2
        customModel.SHAPE_COLUMN = 3
        customModel.TYPE_COLUMN = 4
        customModel.NODE_COLUMN = 5
        customModel.LINK_COLUMN = 6
        
        self.hdf5model = customModel(self.view)
        self.view.setModel(self.hdf5model)

        self.dataviewer = DataViewerFrame.DataViewerFrame()
        self.dataviewerDialog = qt.QDialog(self)
        dvlayout = qt.QVBoxLayout()
        dvlayout.addWidget(self.dataviewer)
        self.dataviewerDialog.setLayout(dvlayout)
        self.dataviewerDialog.setModal(False)
        
        self.view.setExpandsOnDoubleClick(False)
        self.hdf5model.setFileMoveEnabled(True)
        #self.__treeModelSorted = silx.gui.hdf5.NexusSortFilterProxyModel(self.view)
        #self.__treeModelSorted.setSourceModel(self.hdf5model)
        #self.__treeModelSorted.sort(0, qt.Qt.AscendingOrder)
        #self.__treeModelSorted.setSortCaseSensitivity(qt.Qt.CaseInsensitive)

        #self.hdfTreeView.setModel(self.__treeModelSorted)
        
        self.view.doubleClicked.connect(self._onNEXUSDoubleClicked)
        self.view.addContextMenuCallback(self.nexus_treeview_callback)
        
        self.temp_directory = tempfile.TemporaryDirectory(dir=os.path.join(os.getcwd()))
        tempfilepath = os.path.join(self.temp_directory.name,"orgui_database.h5")
        self.createNewDBFile(tempfilepath)
        #self.add_nxdict(gauss)
        #self.add_nxdict(gauss)
        self.setCentralWidget(self.view)
        
        toolbar = qt.QToolBar("Database toolbar",self)
        
        loadDatabaseAct = toolbar.addAction(icons.getQIcon("document-open"),"Open orgui database")
        loadDatabaseAct.triggered.connect(self.onOpenDatabase)
        
        savenewact = toolbar.addAction(icons.getQIcon("layer-nx"),"Select orgui database location")
        savenewact.triggered.connect(self.onSaveNewDBFile)
                
        saveact = toolbar.addAction(icons.getQIcon("document-save"),"Save orgui database")
        saveact.triggered.connect(self.onSaveDBFile)
        
        self.addToolBar(toolbar)
        
    def _onNEXUSDoubleClicked(self,index): # ToDo add try except with popup message!
        nodes = list(self.view.selectedH5Nodes())
        if len(nodes) > 0:
            obj = nodes[0]
            if obj.ntype is h5py.Dataset:
                roi_node = self.get_roinode(obj.h5py_object)
                if roi_node is not None:
                    self.plot_signal_callback(roi_node, obj)
                    return
            self.plot_default(obj.h5py_object)
                
    def plot_default(self, h5py_object):
        try:
            nxdat = nxdata.get_default(h5py_object)
            #print(nxdat.axes, nxdat.signal, nxdat.title)
            #print(nxdat.signal_name, nxdat.errors, nxdat.axes_names)
            if nxdat is not None and len(nxdat.axes) == 1:
                self.plot.addCurve(nxdat.axes[0], nxdat.signal, legend=nxdat.title, xlabel=nxdat.axes_names[0], ylabel=nxdat.signal_name, yerror=nxdat.errors)
        except Exception as e:
            traceback.print_exc()
            print("Cannot plot data: %s" % e)

    def get_roinode(self, obj):
        if silx.io.utils.get_h5_class(obj) is None:
            return None
            
        while(obj.name != '/'):
            meta = obj.attrs.get('orgui_meta', False)
            if meta and meta == 'roi':
                return obj
            obj = obj.parent
            
    def get_scannode(self, obj):
        if silx.io.utils.get_h5_class(obj) is None:
            return None
            
        while(obj.name != '/'):
            meta = obj.attrs.get('orgui_meta', False)
            if meta and meta == 'scan':
                return obj
            obj = obj.parent
            
    def nexus_treeview_callback(self,event):
        objects = list(event.source().selectedH5Nodes())
        obj = objects[0]  # for single selection
        menu = event.menu()
        action = qt.QAction("Refresh", menu)
        action.triggered.connect(lambda:  self.hdf5model.synchronizeH5pyObject(obj))
        menu.addAction(action)
        action = qt.QAction("details", menu)
        action.triggered.connect(lambda:  self.view_data_callback(obj))
        menu.addAction(action)
        
        if obj.ntype is h5py.Dataset:
            roi_node = self.get_roinode(obj.h5py_object)
            if roi_node is not None:
                action = qt.QAction("plot %s" % obj.name, menu)
                action.triggered.connect(lambda:  self.plot_signal_callback(roi_node, obj))
                menu.addAction(action)
        
        if obj.ntype is h5py.Group:
            meta = obj.h5py_object.attrs.get('orgui_meta', False)
            if meta and meta == 'roi':
                action = qt.QAction("rename", menu)
                action.triggered.connect(lambda:  self.onRenameNode(obj.h5py_object))
                menu.addAction(action)
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda:  self.delete_node(obj.h5py_object))
                menu.addAction(action)
                
            if meta and meta == 'scan':
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda:  self.onDeleteScan(obj.h5py_object))
                menu.addAction(action)
            
        """    
        if obj.ntype is h5py.File:
            action = qt.QAction("remove", menu)
            action.triggered.connect(lambda:  self._onCloseFile())
            menu.addAction(action)
        """
    def plot_signal_callback(self, roi_node, dataset):
        try:
            nxdat = nxdata.get_default(roi_node)
            data = dataset.h5py_object[()]
            #print(nxdat.axes, nxdat.signal, nxdat.title)
            #print(nxdat.signal_name, nxdat.errors, nxdat.axes_names)
            if nxdat is not None and len(nxdat.axes) == 1:
                self.plot.addCurve(nxdat.axes[0], data, legend=nxdat.title + "_" + dataset.name, xlabel=nxdat.axes_names[0], ylabel=dataset.name)
        except Exception as e:
            traceback.print_exc()
            print("Cannot plot data: %s" % e)
    
        
        
    def view_data_callback(self,obj):
        self.dataviewer.setData(obj)
        self.dataviewerDialog.open()
        
    def onSaveNewDBFile(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Select database location",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        self.saveNewDBFile(filename)
        
    def onSaveDBFile(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Save database",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        self.saveDBFile(filename)
        
    def onOpenDatabase(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Save database",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        try:
            self.openDBFile(filename)
        except:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot open db file', 
            'Cannot open db file %s.' % filename,
            qt.QMessageBox.Ok, self)
            msgbox.setDetailedText(traceback.format_exc())
            clickedbutton = msgbox.exec()
        
        
    def saveNewDBFile(self, filename):
        alldata = nxtodict(self.nxfile)
        self.createNewDBFile(filename, alldata)
        
    def saveDBFile(self, filename):
        alldata = nxtodict(self.nxfile)
        dicttonx(alldata, filename)
        
    def createNewDBFile(self, filename, datadict=None):
        if self.nxfile is not None:
            self.hdf5model.removeH5pyObject(self.nxfile)
            self.nxfile.close()
            self.nxfile = None
            if hasattr(self, "temp_directory"):
                del self.temp_directory

        fileattrs = {"@NX_class": u"NXroot",
                     "@creator": u"orGUI version %s" % __version__,
                     "@file_name": str(os.path.basename(filename)),
                     "@file_time": datetime.datetime.utcnow().isoformat()}
        if datadict is None:
            datadict = fileattrs
        else:
            datadict.update(fileattrs)
        try:
            dicttonx(datadict, filename)
            self.openDBFile(filename)
        except:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot create db file', 
            'Cannot create file %s.' % filename,
            qt.QMessageBox.Ok, self)
            msgbox.setDetailedText(traceback.format_exc())
            clickedbutton = msgbox.exec()
        
    def openDBFile(self, filename):
        if self.nxfile is not None:
            self.hdf5model.removeH5pyObject(self.nxfile)
            self.nxfile.close()
            if hasattr(self, "temp_directory"):
                del self.temp_directory
        self.nxfile = silx.io.h5py_utils.File(filename,'a')
        self._filepath = filename
        self.hdf5model.insertH5pyObject(self.nxfile)
        self.view.expandToDepth(0)

    def add_nxdict(self, nxentry):
        if self.nxfile is None:
            raise Exception("No database file open.")
        dicttonx(nxentry, self.nxfile, update_mode='add')
        self.hdf5model.synchronizeH5pyObject(self.nxfile)
        self.view.expandToDepth(0)
        
    def onDeleteScan(self, obj):
        btn = qt.QMessageBox.question(self,"Delete scan?","Are you sure that you want to delete %s from the orgui database?" % obj.name.split("/")[-1])
        if btn == qt.QMessageBox.Yes:
            self.delete_node(obj)
        
    def delete_node(self, obj):
        basename = obj.name.split("/")[-1]
        objpar = obj.parent
        del objpar[basename]
        self.hdf5model.synchronizeH5pyObject(self.nxfile)
        self.view.expandToDepth(0)
        
    def onRenameNode(self, obj):
        basename = obj.name.split("/")[-1]
        newname, success = qt.QInputDialog.getText(self,"Rename NEXUS node",
                               "New name:", qt.QLineEdit.EchoMode.Normal,
                               basename)
        if success and newname != '':
            self.rename_node(obj, newname)
                               
    def rename_node(self, obj, newname):
        basename = obj.name.split("/")[-1]
        objpar = obj.parent
        objpar.move(basename, newname)
        self.hdf5model.synchronizeH5pyObject(self.nxfile)
        self.view.expandToDepth(0)
        
        
    def close(self):
        if self.nxfile is not None:
            self.hdf5model.removeH5pyObject(self.nxfile)
            self.nxfile.close()
            self.nxfile = None
            if hasattr(self, "temp_directory"):
                del self.temp_directory
        
    
    
    

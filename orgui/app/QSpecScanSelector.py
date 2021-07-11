# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2020 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020, Timo Fuchs, Olaf Magnussen all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"



import sys
import os
import shutil



from silx.gui import qt

from  silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
from silx.gui.dialog import ImageFileDialog
from silx.gui import icons

#from silx.gui.widgets import HorizontalSliderWithBrowser

import silx.gui.hdf5
from silx.gui.data import DataViewerFrame
import h5py

import warnings


#class ThBrowser(HorizontalSliderWithBrowser):
# todo    
    
        


class QSpecScanSelector(qt.QMainWindow):
    sigScanChanged = qt.pyqtSignal(object)
    sigImagePathChanged = qt.pyqtSignal(object)
    sigImageNoChanged = qt.pyqtSignal(object)
    sigROIChanged = qt.pyqtSignal()
    sigROIintegrate = qt.pyqtSignal()
    def __init__(self,parentmainwindow , parent=None):
        qt.QMainWindow.__init__(self ,parent=None)
        self.parentmainwindow = parentmainwindow
        
        self.mainwidget = qt.QWidget()
        self.mainLayout = qt.QVBoxLayout()
        maintab = qt.QTabWidget()
        
        
        self.openFileAction = qt.QAction(icons.getQIcon('document-open'),"Open file",self)
        self.openFileAction.triggered.connect(self._onOpenFile)
        
        self.refreshFileAction = qt.QAction(icons.getQIcon('view-refresh'),"Refresh file",self)
        self.refreshFileAction.triggered.connect(self._onRefreshFile)
        
        self.closeFileAction = qt.QAction(icons.getQIcon('close'),"Close file",self)
        self.closeFileAction.triggered.connect(self._onCloseFile)
        
        
        self.hdfTreeView = silx.gui.hdf5.Hdf5TreeView(self)
        self.hdfTreeView.setSortingEnabled(True)
        self.hdfTreeView.addContextMenuCallback(self.nexus_treeview_callback)
        self.hdf5model = Hdf5TreeModel(self.hdfTreeView,ownFiles=True)
        self.hdfTreeView.setModel(self.hdf5model)
        self.hdfTreeView.setExpandsOnDoubleClick(False)
        self.hdf5model.setFileMoveEnabled(True)
        self.__treeModelSorted = silx.gui.hdf5.NexusSortFilterProxyModel(self.hdfTreeView)
        self.__treeModelSorted.setSourceModel(self.hdf5model)
        self.__treeModelSorted.sort(0, qt.Qt.AscendingOrder)
        self.__treeModelSorted.setSortCaseSensitivity(qt.Qt.CaseInsensitive)

        self.hdfTreeView.setModel(self.__treeModelSorted)
        
        self.hdfTreeView.doubleClicked.connect(self._onNEXUSChange)
        
        self.dataviewer = DataViewerFrame.DataViewerFrame()
        self.dataviewerDialog = qt.QDialog(self)
        dvlayout = qt.QVBoxLayout()
        dvlayout.addWidget(self.dataviewer)
        self.dataviewerDialog.setLayout(dvlayout)
        self.dataviewerDialog.setModal(False)
        
        maintab.addTab(self.hdfTreeView,"NEXUS")
        
        pathSelector = qt.QSplitter(self)
        pathSelector.setOrientation(qt.Qt.Horizontal)
        qt.QLabel("File path:",pathSelector)

        self.pathedit = qt.QLineEdit(pathSelector)
        
        openButton = qt.QPushButton(pathSelector)
        openButton.setIcon(icons.getQIcon('document-open'))
        openButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        openButton.setToolTip("Choose the NEXUS file")
        openButton.clicked.connect(self._onSelectFilePath)
        
        scannoSelector = qt.QSplitter(self)
        scannoSelector.setOrientation(qt.Qt.Horizontal)
        qt.QLabel("Scan #:",scannoSelector)

        self.scannoBox = qt.QSpinBox(scannoSelector)
        self.scannoBox.setRange(1,2147483647)
        self.scannoBox.setValue(1)
        
        openScanButton = qt.QPushButton(scannoSelector)
        openScanButton.setIcon(icons.getQIcon('selected'))
        openScanButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        openScanButton.setToolTip("Open scan")
        openScanButton.clicked.connect(self._onLoadScan)

        
        self.slider = qt.QSlider()
        self.slider.setOrientation(qt.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)

        
        self.mainLayout.addWidget(maintab)
        self.mainLayout.addWidget(pathSelector)
        self.mainLayout.addWidget(scannoSelector)
        #self.mainLayout.addWidget(self.slider)
        
        
        
        self.mainwidget.setLayout(self.mainLayout)
        self.setCentralWidget(self.mainwidget)

        self.toolbar = qt.QToolBar("Image selector",self)
        
        
        #righticon = qt.QIcon(qt.QPixmap(icons.rotate_right))
        #lefticon = qt.QIcon(qt.QPixmap(icons.rotate_left))
        #increaseImageNo = qt.QAction(righticon,"next image")
        #decreaseImageNo = qt.QAction(lefticon,"previous image")
        
        
        imglabel = qt.QLabel("No:")
        self.noSelector = qt.QSpinBox()
        self.noSelector.setRange(0,0)
        
        thlabel = qt.QLabel("th:")
        self.thSelector = qt.QDoubleSpinBox()
        self.thSelector.setRange(-1000,1000)
        self.thSelector.setReadOnly(True)
        self.thSelector.setSuffix(" °")
        
        
        
        self.toolbar.addWidget(imglabel)
        self.toolbar.addWidget(self.noSelector)
        self.toolbar.addWidget(thlabel)
        self.toolbar.addWidget(self.thSelector)
        
        self.toolbar.addWidget(self.slider)
        decreaseImageNo = self.toolbar.addAction(icons.getQIcon("previous"),"previous image")
        increaseImageNo = self.toolbar.addAction(icons.getQIcon("next"),"next image")
        
        
        increaseImageNo.setShortcut(qt.QKeySequence( qt.Qt.Key_Plus))
        decreaseImageNo.setShortcut(qt.QKeySequence( qt.Qt.Key_Minus))
        
        increaseImageNo.triggered.connect(self._onIncreaseImageNo)
        decreaseImageNo.triggered.connect(self._onDecreaseImageNo)
        
        
        self.slider.valueChanged.connect(self._onSliderChanged)
        self.noSelector.valueChanged.connect(self.slider.setValue)
        self.slider.valueChanged.connect(self.noSelector.setValue)
        
        
        #self.addToolBar(qt.Qt.BottomToolBarArea,self.toolbar)
        self.th = None            
            
            
        
        self.selectedScan = None
        
        ### SCAN DISPLAY
        
        self.scanDisplayTab = qt.QWidget()
        self.scanDisplayTabLayout = qt.QVBoxLayout()
        
        ## Scan display
        
        loadScanGroup = qt.QGroupBox("load scan required to display max/sum")
        self.loadScanGroupLayout = qt.QVBoxLayout()
        
        loadScanGroupButtonlayout = qt.QHBoxLayout()
        self.loadallButton = qt.QPushButton("load scan")
        self.loadallButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.loadallButton.setToolTip("open all images and calculate sum and maximum")
        #self.loadallButton.clicked.connect(self._onLoadAll)
        
        self.showMaxButton = qt.QPushButton("show max")
        self.showMaxButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.showMaxButton.setToolTip("plot maximum of the scan")
        self.showMaxButton.setCheckable(True)
        #showMaxButton.clicked.connect(self._onLoadAll)
        
        self.showSumButton = qt.QPushButton("show sum")
        self.showSumButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.showSumButton.setToolTip("plot sum of the scan")
        self.showSumButton.setCheckable(True)
        
        loadScanGroupButtonlayout.addWidget(self.loadallButton)
        loadScanGroupButtonlayout.addWidget(self.showMaxButton)
        loadScanGroupButtonlayout.addWidget(self.showSumButton)
        
        self.loadScanGroupLayout.addLayout(loadScanGroupButtonlayout)
        
        
        
        loadScanGroup.setLayout(self.loadScanGroupLayout)
        self.scanDisplayTabLayout.addWidget(loadScanGroup)
        

        
        self.scanDisplayTab.setLayout(self.scanDisplayTabLayout)
        maintab.addTab(self.scanDisplayTab,"scan display")
        
        
        ## ROI
                
        self.roiIntegrateTab = qt.QWidget()
        self.roiIntegrateTabLayout = qt.QVBoxLayout()
        
        roiGroup = qt.QGroupBox("ROI definition (in pixel)")
        roiGroupLayout = qt.QGridLayout()
        
        self.hsize = qt.QDoubleSpinBox()
        self.vsize = qt.QDoubleSpinBox()
        self.left = qt.QDoubleSpinBox()
        self.right = qt.QDoubleSpinBox()
        self.top = qt.QDoubleSpinBox()
        self.bottom = qt.QDoubleSpinBox()
        
        self.hsize.setRange(1,20000)
        self.hsize.setDecimals(1)
        self.hsize.setSuffix(" px")
        self.hsize.setValue(6.)
        
        self.vsize.setRange(1,20000)
        self.vsize.setDecimals(1)
        self.vsize.setSuffix(" px")
        self.vsize.setValue(6.)
        
        self.left.setRange(0,20000)
        self.left.setDecimals(1)
        self.left.setSuffix(" px")
        self.left.setValue(6.)
        
        self.right.setRange(0,20000)
        self.right.setDecimals(1)
        self.right.setSuffix(" px")
        self.right.setValue(6.)

        self.top.setRange(0,20000)
        self.top.setDecimals(1)
        self.top.setSuffix(" px")
        self.top.setValue(0.)
        
        self.bottom.setRange(0,20000)
        self.bottom.setDecimals(1)
        self.bottom.setSuffix(" px")
        self.bottom.setValue(0.)

        roiGroupLayout.addWidget(qt.QLabel('center roi (h x v):'),0,0)
        roiGroupLayout.addWidget(self.hsize,0,1)
        roiGroupLayout.addWidget(self.vsize,0,2)

        roiGroupLayout.addWidget(qt.QLabel('bg roi (left, right):'),1,0)
        roiGroupLayout.addWidget(self.left,1,1)
        roiGroupLayout.addWidget(self.right,1,2)
        
        roiGroupLayout.addWidget(qt.QLabel('bg roi (top, bottom):'),2,0)
        roiGroupLayout.addWidget(self.top,2,1)
        roiGroupLayout.addWidget(self.bottom,2,2)

        self.hsize.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.vsize.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.left.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.right.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.top.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.bottom.valueChanged.connect(lambda : self.sigROIChanged.emit())
                
        roiGroup.setLayout(roiGroupLayout)
        
        self.H_0 = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.H_0]
        [h.setDecimals(2) for h in self.H_0]
        self.H_0[0].setValue(1.)
        self.H_0[1].setValue(0.)
        self.H_0[2].setValue(0.)
        [h.valueChanged.connect(lambda : self.sigROIChanged.emit()) for h in self.H_0]
        
        
        self.H_1 = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.H_1]
        [h.setDecimals(2) for h in self.H_1]
        self.H_1[0].setValue(0.)
        self.H_1[1].setValue(0.)
        self.H_1[2].setValue(1.)
        [h.valueChanged.connect(lambda : self.sigROIChanged.emit()) for h in self.H_1]
        
        directionGroup = qt.QGroupBox("Direction H_1 (hkl)")
        directionGroupLayout = qt.QHBoxLayout()
        [directionGroupLayout.addWidget(h) for h in self.H_1]
        directionGroup.setLayout(directionGroupLayout)
        
        locationGroup = qt.QGroupBox("Location vector H_0 (hkl)")
        locationGroupLayout = qt.QHBoxLayout()
        [locationGroupLayout.addWidget(h) for h in self.H_0]
        locationGroup.setLayout(locationGroupLayout)
        
        self.roiIntegrateTabLayout.addWidget(roiGroup)
        self.roiIntegrateTabLayout.addWidget(qt.QLabel("Integrate along a line given by:\nH(s) = H_1 * s + H_0"))
        self.roiIntegrateTabLayout.addWidget(directionGroup)
        self.roiIntegrateTabLayout.addWidget(locationGroup)
        
        self.integrateROIBtn = qt.QPushButton("ROI integrate scan")
        self.integrateROIBtn.clicked.connect(lambda : self.sigROIintegrate.emit())
        self.roiIntegrateTabLayout.addWidget(self.integrateROIBtn)
        
        #self.showROICheckBox = qt.QCheckBox("Show ROI")
        #self.roiIntegrateTabLayout.addWidget(self.showROICheckBox)
        
        self.roiIntegrateTab.setLayout(self.roiIntegrateTabLayout)
        
        
        
        maintab.addTab(self.roiIntegrateTab,"ROI integration")


        
    #def setImageNo(self,imageno):
    #    self.slider.setValue(imageno)
        
    #def reload_nx_callback(self,obj):
    #    self.hdf5model.synchronizeH5pyObject(obj)
        
    def view_data_callback(self,obj):
        self.dataviewer.setData(obj)
        self.dataviewerDialog.open()
        
    def _onLoadScan(self):
        ddict = dict()
        ddict['event'] = "loadScan"
        ddict['file'] = self.pathedit.text()
        ddict['scanno'] = self.scannoBox.value()
        self.sigScanChanged.emit(ddict)
        
    def _onOpenFile(self):
        fileTypeDict = {'NEXUS files (*.h5 *.hdf5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open NEXUS file",
                                                  self.parentmainwindow.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.parentmainwindow.filedialogdir = os.path.splitext(filename)[0]
        self.hdf5model.appendFile(filename)
        
    def _onSelectFilePath(self):
        fileTypeDict = {'NEXUS files (*.h5 *.hdf5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open NEXUS file",
                                                  self.parentmainwindow.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.parentmainwindow.filedialogdir = os.path.splitext(filename)[0]
        self.pathedit.setText(filename)
        
        
    def _onCloseFile(self):
        objects = list(self.hdfTreeView.selectedH5Nodes())
        if len(objects) > 0:
            obj = objects[0]
            self.hdf5model.removeH5pyObject(obj.file)

    def nexus_treeview_callback(self,event):
        objects = list(event.source().selectedH5Nodes())
        obj = objects[0]  # for single selection
        menu = event.menu()
        action = qt.QAction("Refresh", menu)
        action.triggered.connect(lambda:  self.hdf5model.synchronizeH5pyObject(obj))
        menu.addAction(action)
        if obj.ntype is h5py.Dataset:
            action = qt.QAction("display data", menu)
            action.triggered.connect(lambda:  self.view_data_callback(obj))
            menu.addAction(action)
        if obj.ntype is h5py.File:
            action = qt.QAction("remove", menu)
            action.triggered.connect(lambda:  self._onCloseFile())
            menu.addAction(action)
            
     
    def _onRefreshFile(self):
        objects = list(self.hdfTreeView.selectedH5Nodes())
        if len(objects) > 0:
            obj = objects[0]
            self.hdf5model.synchronizeH5pyObject(obj)
        
        
        #if 'NX_class' in obj.ntype.attrs:
        #    if obj.ntype.attrs['NX_class'] == 'NX_':
              
    
    
    def getScanToolbar(self):
        return self.toolbar
        
    def showToolBar(self):
        self.addToolBar(qt.Qt.BottomToolBarArea,self.toolbar)
    
    def _onSliderChanged(self,scanno):
        #self.noSelector.setValue(self.slider.value())
        self.thSelector.setValue(self.th[self.slider.value()])
        self.sigImageNoChanged.emit(self.slider.value())
    
    def _onNoSelectorChanged(self):
        self.slider.setValue(self.noSelector.value())
        self.sigImageNoChanged.emit(self.noSelector.value())
    
    def _onDecreaseImageNo(self):
        self.slider.setValue(self.slider.value() -1)
    
    def _onIncreaseImageNo(self):
        self.slider.setValue(self.slider.value() +1)
        
    def setRange(self,minimum,maximum):
        self.slider.setRange(minimum,maximum)
        self.noSelector.setRange(minimum,maximum)
        
    def setTh(self,th):
        self.th = th
        self.setRange(0,th.size-1)
    """    
    def _onSourceSelected(self,ddict):
        event = ddict['event']
        if (event == 'SourceSelected' or event == 'NewSourceSelected' or event == 'SourceReloaded'):
            dataname = ddict['sourcelist'][0]
            
            datatype = QDataSource.getSourceType(dataname)
            datapath,_ = os.path.split(dataname)
            if datatype == NexusDataSource.SOURCE_TYPE:
                #shutil.copy(dataname,datapath + '/temp.h5')
                #ds = QDataSource.QDataSource(datapath + '/temp.h5',datatype)
                if event == 'SourceReloaded':
                    pass
                else:
                    self.hdf5model.appendFile(dataname)
            elif datatype == SpecFileDataSource.SOURCE_TYPE:
                try:
                    ds = QDataSource.QDataSource(dataname,datatype)
                    #self.specfileWidget.setDataSource(ds)
                except Exception:
                    warnings.warn("didn't find good datasource in %s , try to read as P212 CrudeScan" % dataname)
                    self.sigScanChanged.emit([{'SourceName' : dataname}])
                    self.selectedScan = None
                    return 
            else:
                warnings.warn("not implemented data source: %s" % dataname)
            self.selectedScan = None
            self.sigScanChanged.emit([])
            
    """

    
    def _onReplaceSelection(self,sel_list):
        print(sel_list)
        self.sigScanChanged.emit(sel_list)
        
    def _onNEXUSChange(self,index):
        nodes = list(self.hdfTreeView.selectedH5Nodes())
        if len(nodes) > 0:
            obj = nodes[0]
            if 'NX_class' in obj.attrs:
                try:
                    nxcls = obj.attrs['NX_class'].decode("utf-8")
                    ch5523bliss = True
                except AttributeError:
                    nxcls = obj.attrs['NX_class']
                    ch5523bliss = False

                if nxcls == 'NXentry':
                    if 'instrument' in obj.h5py_target and 'fiofile' in obj.h5py_target['instrument']:
                        p212H5 = True
                        scanname = obj.basename
                        scanno, subscanno = scanname.split('.')
                    elif ch5523bliss:
                        p212H5 = False
                        scanname = obj.local_name
                        scanno = scanname.split('_')[-1]
                    else:
                        p212H5 = False
                        scanname = obj.local_name
                        scansuffix = scanname.split('_')[-1]
                        scanname_nosuffix = '_'.join(scanname.split('_')[:-1])
                        scanno, subscanno = scansuffix.split('.')
                        
                    ddict = dict()
                    ddict['event'] = "itemDoubleClicked"
                    ddict['file'] = obj.local_filename
                    ddict['name'] = obj.local_name
                    ddict['node'] = obj
                    ddict['scanno'] = int(scanno)
                    ddict['ch5523'] = ch5523bliss
                    ddict['p212H5'] = p212H5
                    self.pathedit.setText(obj.local_filename)
                    self.scannoBox.setValue(int(scanno))
                    self.sigScanChanged.emit(ddict)
    """    
    def _onSelectImageFolder(self):
        
        folder = PyMcaFileDialogs.getExistingDirectory(self,"select directory containing the images for the current scan")
        if(len(folder)):
            self.pathedit.setText(folder)
            self._onAcceptImagePath()
    """    
        
    def _onAcceptImagePath(self):
        #print(self.pathedit.text())
        self.sigImagePathChanged.emit(self.pathedit.text())

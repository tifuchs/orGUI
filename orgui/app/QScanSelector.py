# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
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
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"


import sys
import os
import shutil
from dateutil import parser as dateparser
from datetime import datetime

from silx.gui import qt

from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
from silx.gui.dialog import ImageFileDialog
from silx.gui import icons
from silx.gui.plot.AlphaSlider import NamedImageAlphaSlider
#from silx.gui.widgets import HorizontalSliderWithBrowser

import silx.gui.hdf5
from silx.gui.data import DataViewerFrame
import h5py
import traceback

import warnings
from .. import resources
from ..backend import backends, scans
from . import qutils
from .QReflectionSelector import QReflectionAnglesDialog
import runpy


class QScanSelector(qt.QMainWindow):
    sigScanChanged = qt.pyqtSignal(object)
    sigImagePathChanged = qt.pyqtSignal(object)
    sigImageNoChanged = qt.pyqtSignal(object)
    sigROIChanged = qt.pyqtSignal()
    sigROIintegrate = qt.pyqtSignal()
    sigSearchHKL = qt.pyqtSignal(list)
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
        
        self.loadBackendAction = qt.QAction("Load backend file",self)
        self.loadBackendAction.triggered.connect(self._onLoadBackend)
        
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
        
        self.hdfTreeView.doubleClicked.connect(self._onNEXUSDoubleClicked)
        
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
        openButton.setToolTip("Choose the NEXUS/SPEC/log file")
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

        btidsplit = qt.QSplitter(self)
        qt.QLabel("Backend:",btidsplit)
        self.btid = qt.QComboBox(btidsplit)
        [self.btid.addItem(bt) for bt in backends.fscans]
        self.btid.setCurrentText("id31_default")
        
        self._selectBackendBtn = qt.QPushButton("...", btidsplit)
        width = self._selectBackendBtn.fontMetrics().boundingRect("  ...  ").width() + 7
        height = self._selectBackendBtn.fontMetrics().boundingRect("  M  ").height() + 7
        self._selectBackendBtn.setMaximumWidth(width)
        self._selectBackendBtn.setMaximumHeight(height)
        self._selectBackendBtn.clicked.connect(self.loadBackendAction.trigger)
        
        self.bt_autodetect_enable = qt.QCheckBox("auto detect", btidsplit)
        self.bt_autodetect_enable.toggled.connect(lambda s : self.btid.setEnabled(not s))
        self.bt_autodetect_enable.setChecked(True)
        
        
        self.slider = qt.QSlider()
        self.slider.setOrientation(qt.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)

        
        self.mainLayout.addWidget(maintab)
        self.mainLayout.addWidget(pathSelector)
        self.mainLayout.addWidget(scannoSelector)
        self.mainLayout.addWidget(btidsplit)
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
        
        self.axislabel = qt.QLabel("axis:")
        self.axisSelector = qt.QDoubleSpinBox()
        self.axisSelector.setRange(-1000,1000)
        self.axisSelector.setReadOnly(True)
        self.axisSelector.setSuffix(" °")
        
        self.showMaxAct = qt.QAction(resources.getQicon("max_image2"), "plot maximum of the scan")
        self.showMaxAct.setCheckable(True)
        
        self.showSumAct = qt.QAction(resources.getQicon("sum_image2"), "plot sum of the scan")
        self.showSumAct.setCheckable(True)
        
        self.excludeImageAct = qt.QAction(resources.getQicon("disable-image"), "exclude image from max/sum image")
        self.excludeImageAct.setCheckable(True)
        
        self.alphaslider = NamedImageAlphaSlider(self,self.parentmainwindow.centralPlot,self.parentmainwindow.currentAddImageLabel)
        self.alphaslider.setOrientation(qt.Qt.Horizontal)
        self.alphaslider.setEnabled(True)
        
        self.alpha_menu = qt.QMenu()
        
        self.alphasliderwidget = qt.QWidgetAction(self.alpha_menu)
        self.alphasliderwidget.setDefaultWidget(self.alphaslider)

        self.alpha_menu.addAction(self.alphasliderwidget)
        
        #self.alpha_btn = qt.QToolButton(resources.getQicon("sum_image.png"),"slider")
        self.alpha_btn = qt.QToolButton()
        self.alpha_btn.setIcon(resources.getQicon("alpha"))
        self.alpha_btn.setToolTip("Transparency of max / sum image")
        self.alpha_btn.setPopupMode(qt.QToolButton.InstantPopup)
        self.alpha_btn.setMenu(self.alpha_menu)
        
        self.alpha_btn_act = qt.QWidgetAction(self)
        self.alpha_btn_act.setDefaultWidget(self.alpha_btn)
        

        self.toolbar.addAction(self.showMaxAct)
        self.toolbar.addAction(self.showSumAct)
        self.toolbar.addAction(self.alpha_btn_act)
        self.toolbar.addSeparator()
        
        self.toolbar.addAction(self.excludeImageAct)
        
        self.toolbar.addSeparator()
        
        self.toolbar.addWidget(imglabel)
        self.toolbar.addWidget(self.noSelector)
        self.toolbar.addWidget(self.axislabel)
        self.toolbar.addWidget(self.axisSelector)
        
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
        self.axis = None            
            
            
        
        self.selectedScan = None
        
        
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
        self.offsetx = qt.QDoubleSpinBox()
        self.offsety = qt.QDoubleSpinBox()
        
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
        
        self.offsetx.setRange(-20000,20000)
        self.offsetx.setDecimals(1)
        self.offsetx.setSuffix(" px")
        self.offsetx.setValue(0.)
        
        self.offsety.setRange(-20000,20000)
        self.offsety.setDecimals(1)
        self.offsety.setSuffix(" px")
        self.offsety.setValue(0.)

        roiGroupLayout.addWidget(qt.QLabel('center roi (h x v):'),0,0)
        roiGroupLayout.addWidget(self.hsize,0,1)
        roiGroupLayout.addWidget(self.vsize,0,2)

        roiGroupLayout.addWidget(qt.QLabel('bg roi (left, right):'),1,0)
        roiGroupLayout.addWidget(self.left,1,1)
        roiGroupLayout.addWidget(self.right,1,2)
        
        roiGroupLayout.addWidget(qt.QLabel('bg roi (top, bottom):'),2,0)
        roiGroupLayout.addWidget(self.top,2,1)
        roiGroupLayout.addWidget(self.bottom,2,2)
        
        roiGroupLayout.addWidget(qt.QLabel('roi loc offset (x, y):'),3,0)
        roiGroupLayout.addWidget(self.offsetx,3,1)
        roiGroupLayout.addWidget(self.offsety,3,2)

        self.hsize.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.vsize.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.left.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.right.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.top.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.bottom.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.offsetx.valueChanged.connect(lambda : self.sigROIChanged.emit())
        self.offsety.valueChanged.connect(lambda : self.sigROIChanged.emit())
                
        roiGroup.setLayout(roiGroupLayout)
        

        
        
        # hkl scan
        
        self.H_0 = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.H_0]
        [h.setDecimals(4) for h in self.H_0]
        self.H_0[0].setValue(1.)
        self.H_0[1].setValue(0.)
        self.H_0[2].setValue(0.)
        [h.valueChanged.connect(lambda : self.sigROIChanged.emit()) for h in self.H_0]
        
        
        self.H_1 = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.H_1]
        [h.setDecimals(4) for h in self.H_1]
        self.H_1[0].setValue(0.)
        self.H_1[1].setValue(0.)
        self.H_1[2].setValue(1.)
        [h.valueChanged.connect(lambda : self.sigROIChanged.emit()) for h in self.H_1]
        
        directionGroup = qt.QGroupBox(u"Direction H₁ (hkl)")
        directionGroupLayout = qt.QHBoxLayout()
        [directionGroupLayout.addWidget(h) for h in self.H_1]
        directionGroup.setLayout(directionGroupLayout)
        
        locationGroup = qt.QGroupBox(u"Location vector H₀ (hkl)")
        locationGroupLayout = qt.QHBoxLayout()
        [locationGroupLayout.addWidget(h) for h in self.H_0]
        locationGroup.setLayout(locationGroupLayout)
        
        
        hklscanwidget = qt.QWidget()
        hklscanwidgetlayout = qt.QVBoxLayout()
        
        hklscanwidgetlayout.addWidget(qt.QLabel(u"Integrate along:\nH(s) = H₁ 🞄 s + H₀"))
        hklscanwidgetlayout.addWidget(directionGroup)
        hklscanwidgetlayout.addWidget(locationGroup)
        
        hklscanwidget.setLayout(hklscanwidgetlayout)
        # static roi scan
        
        self.xy_static = [qt.QDoubleSpinBox() for i in range(2)]
        
        [h.setRange(-20000,20000) for h in self.xy_static]
        [h.setDecimals(3) for h in self.xy_static]
        self.xy_static[0].setValue(10.)
        self.xy_static[1].setValue(10.)
        [h.valueChanged.connect(lambda : self.sigROIChanged.emit()) for h in self.xy_static]
        
        self.hkl_static = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.hkl_static]
        [h.setDecimals(3) for h in self.hkl_static]
        [h.setValue(0.) for h in self.hkl_static]
        
        
        setroi_btn = qt.QToolButton()
        self.select_roi_action = qt.QAction(icons.getQIcon("crosshair"), "Select roi location by double clicking", self)
        self.select_roi_action.setCheckable(True)

        #setroi_btn.setIcon()
        #setroi_btn.setToolTip("Select roi location by double clicking")
        setroi_btn.setDefaultAction(self.select_roi_action)

        static_loc_Group = qt.QGroupBox(u"Static ROI location")
        static_loc_GroupMainVLayout = qt.QVBoxLayout()
        
        static_pix_xy_box = qt.QGroupBox(u"Pixel coordinates")
        static_loc_GroupLayout = qt.QHBoxLayout()
        static_loc_GroupLayout.addWidget(setroi_btn)
        for spinbox, lbl in zip(self.xy_static, ["x:", "y:"]):
            static_loc_GroupLayout.addWidget(qt.QLabel(lbl))
            static_loc_GroupLayout.addWidget(spinbox)
        static_pix_xy_box.setLayout(static_loc_GroupLayout)
        
        static_hkl_box = qt.QGroupBox(u"hkl coordinates")
        static_loc_HKLLayout = qt.QHBoxLayout()
        for spinbox, lbl in zip(self.hkl_static, ["H:", "K:", "L:"]):
            static_loc_HKLLayout.addWidget(qt.QLabel(lbl))
            static_loc_HKLLayout.addWidget(spinbox)
        calc_HKL_roi_btn = qt.QToolButton()
        self.roi_fromHKL_action = qt.QAction(resources.getQicon("search"), "Set ROI position to calculated position of reflection", self)
        self.roi_fromHKL_action.triggered.connect(self.search_hkl_static_loc)
        calc_HKL_roi_btn.setDefaultAction(self.roi_fromHKL_action)
        static_loc_HKLLayout.addWidget(calc_HKL_roi_btn)
        static_hkl_box.setLayout(static_loc_HKLLayout)
        
        
        static_loc_GroupMainVLayout.addWidget(static_pix_xy_box)
        static_loc_GroupMainVLayout.addWidget(static_hkl_box)
        
        static_loc_Group.setLayout(static_loc_GroupMainVLayout)



        # rocking scan
        
        self.hkl_static1 = [qt.QDoubleSpinBox() for i in range(3)]
        [h.setRange(-20000,20000) for h in self.hkl_static1]
        [h.setDecimals(3) for h in self.hkl_static1]
        [h.setValue(0.) for h in self.hkl_static1]
        self.hkl_static1[0].setValue(1.) # default H
        self.hkl_static1[2].setValue(1.) # default L

        self.roscanMaxL = qt.QDoubleSpinBox()
        self.roscanMaxL.setRange(-20000,20000)
        self.roscanMaxL.setDecimals(3)
        self.roscanMaxL.setValue(2.)

        self.roscanDeltaL = qt.QDoubleSpinBox()
        self.roscanDeltaL.setRange(-20000,20000)
        self.roscanDeltaL.setDecimals(3)
        self.roscanDeltaL.setValue(0.1)

        rocking_integration_group = qt.QGroupBox(u"Integration coordinates")
        rocking_integration_groupMainVLayout = qt.QVBoxLayout()
        
        rocking_integration_group_hkl_box = qt.QGroupBox(u"Reciprocal start coordinates")
        rocking_integration_group_loc_HKLLayout = qt.QHBoxLayout()
        for pinbox, lbll in zip(self.hkl_static1, ["H:", "K:", "L:"]):
            rocking_integration_group_loc_HKLLayout.addWidget(qt.QLabel(lbll))
            rocking_integration_group_loc_HKLLayout.addWidget(pinbox)
        calc_HKL_roi_btn2 = qt.QToolButton()
        self.roi_fromHKL_action2 = qt.QAction(resources.getQicon("search"), "Set ROI position to calculated position of reflection", self)
        self.roi_fromHKL_action2.triggered.connect(self.search_hkl_static_loc)
        calc_HKL_roi_btn2.setDefaultAction(self.roi_fromHKL_action2)
        rocking_integration_group_loc_HKLLayout.addWidget(calc_HKL_roi_btn2)
        rocking_integration_group_hkl_box.setLayout(rocking_integration_group_loc_HKLLayout)
        
        rocking_integration_group_step_box = qt.QGroupBox(u"Integration path")
        rocking_integration_group_step_Layout = qt.QHBoxLayout()
        rocking_integration_group_step_Layout.addWidget(qt.QLabel(u"L max"))
        rocking_integration_group_step_Layout.addWidget(self.roscanMaxL)

        rocking_integration_group_step_Layout.addWidget(qt.QLabel(u"L step"))
        rocking_integration_group_step_Layout.addWidget(self.roscanDeltaL)
        self.autoROIVsize = qt.QCheckBox("auto ROI-v")
        rocking_integration_group_step_Layout.addWidget(self.autoROIVsize)

        rocking_integration_group_step_box.setLayout(rocking_integration_group_step_Layout)
        
        rocking_integration_groupMainVLayout.addWidget(rocking_integration_group_hkl_box)
        rocking_integration_groupMainVLayout.addWidget(rocking_integration_group_step_box)
        
        rocking_integration_group.setLayout(rocking_integration_groupMainVLayout)
        
        
        #  roi scan tab
        
        self.scanstab = qt.QTabWidget()
        self.scanstab.addTab(hklscanwidget, "hklscan")
        self.scanstab.addTab(static_loc_Group, "fixed roi loc")
        self.scanstab.addTab(rocking_integration_group, "rocking scan integration")
        self.scanstab.currentChanged.connect(lambda : self.sigROIChanged.emit())

        
        # options group
        
        optionsGroup = qt.QGroupBox("Integration options")
        optionsGroupLayout = qt.QGridLayout()
        self.useMaskBox = qt.QCheckBox("Use pixel mask")
        self.useLorentzBox = qt.QCheckBox("Lorentz correction")
        self.useLorentzBox.setEnabled(False)
        self.useSolidAngleBox = qt.QCheckBox("Solid angle correction")
        self.usePolarizationBox = qt.QCheckBox("Polarization correction")
        
        optionsGroupLayout.addWidget(self.useMaskBox,0,0)
        optionsGroupLayout.addWidget(self.useLorentzBox,1,0)
        optionsGroupLayout.addWidget(self.useSolidAngleBox,0,1)
        optionsGroupLayout.addWidget(self.usePolarizationBox,1,1)
        
        optionsGroup.setLayout(optionsGroupLayout)
        
        self.roiIntegrateTabLayout.addWidget(roiGroup)
        self.roiIntegrateTabLayout.addWidget(self.scanstab)
        self.roiIntegrateTabLayout.addWidget(optionsGroup)
        
        self.integrateROIBtn = qt.QPushButton("ROI integrate scan")
        self.integrateROIBtn.clicked.connect(lambda : self.sigROIintegrate.emit())
        self.roiIntegrateTabLayout.addWidget(self.integrateROIBtn)
        
        #self.showROICheckBox = qt.QCheckBox("Show ROI")
        #self.roiIntegrateTabLayout.addWidget(self.showROICheckBox)
        
        self.roiIntegrateTab.setLayout(self.roiIntegrateTabLayout)
        
        
        
        maintab.addTab(self.roiIntegrateTab,"ROI integration")
    
    def set_xy_static_loc(self, x, y):
        [h.blockSignals(True) for h in self.xy_static]
        self.xy_static[0].setValue(x)
        self.xy_static[1].setValue(y)
        [h.blockSignals(False) for h in self.xy_static]
        self.sigROIChanged.emit()
        
    def search_hkl_static_loc(self):
        if self.scanstab.currentIndex() == 2:
            hkl = [h.value() for h in self.hkl_static1]
            self.sigSearchHKL.emit(hkl)    
        else:
            hkl = [h.value() for h in self.hkl_static]
            self.sigSearchHKL.emit(hkl)        
        
    def view_data_callback(self,obj):
        self.dataviewer.setData(obj)
        self.dataviewerDialog.open()
        
    def _onLoadScan(self):
        if self.bt_autodetect_enable.isChecked():
            msgbox = qt.QMessageBox(qt.QMessageBox.Warning,'Cannot auto detect backend', 
                        'Cannot auto detect beamtime id and corresponding backend in minimal mode.\nPlease first deselect the beamtime auto detection and then chose the correct beamtime or default backend.',
                        qt.QMessageBox.Ok, self)
            clickedbutton = msgbox.exec()
            return
        ddict = dict()
        ddict['event'] = "loadScan"
        ddict['file'] = self.pathedit.text()
        ddict['scanno'] = self.scannoBox.value()
        self.sigScanChanged.emit(ddict)
        
    def _onLoadBackend(self):
        fileTypeDict = {'Python backend files (*.py)': '.py', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open Backend file",
                                                  self.parentmainwindow.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.parentmainwindow.filedialogdir = os.path.splitext(filename)[0]
        
        try:
            self.loadBackendFile(filename)
        except:
            qutils.warning_detailed_message(self, "Cannot load backend", "Cannot load backend", traceback.format_exc())

    def loadBackendFile(self, filename):
        backend_file = runpy.run_path(filename)
        found_backends = []
        for e in backend_file:
            try:
                if issubclass(backend_file[e], scans.Scan) and backend_file[e] != scans.Scan:
                    found_backends.append((e, backend_file[e]))
            except:
                pass
                #traceback.print_exc()
        if not found_backends:
            raise ValueError("Found no backend in file %s" % filename)
        if len(found_backends) > 1:
            raise ValueError("Found more than one Scan class in backend file %s. Only one is permitted" % filename)
        name, scancls = found_backends[0]
        self.btid.addItem(name)
        backends.fscans[name] = scancls
        self.btid.setCurrentText(name)
        self.bt_autodetect_enable.setChecked(False)
        
    def _onOpenFile(self):
        fileTypeDict = {'NEXUS files (*.h5 *.hdf5)': '.h5', "SPEC files (*.spec *.spc)": '.spec', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open NEXUS file",
                                                  self.parentmainwindow.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.parentmainwindow.filedialogdir = os.path.splitext(filename)[0]
        try:
            self.hdf5model.appendFile(filename)
        except:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot open file', 
                        'Cannot open file %s.' % filename,
                        qt.QMessageBox.Ok, self)
            msgbox.setDetailedText(traceback.format_exc())
            clickedbutton = msgbox.exec()
            
        
    def _onSelectFilePath(self):
        fileTypeDict = {'NEXUS files (*.h5 *.hdf5)': '.h5', "SPEC files (*.spec *.spc)": '.spec', 'log files (*.log)' : '.log','All files (*)': '' }
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
        if len(objects) > 0:
            obj = objects[0]  # for single selection
            menu = event.menu()
            action = qt.QAction("Refresh", menu)
            action.triggered.connect(lambda:  self.hdf5model.synchronizeH5pyObject(obj))
            menu.addAction(action)
            #if obj.ntype is h5py.Dataset:
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
        self.axisSelector.setValue(self.axis[self.slider.value()])
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
        
    def setAxis(self,axis,label='th'):
        self.axislabel.setText(label)
        self.axis = axis
        self.setRange(0,axis.size-1)
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
        
    def _onNEXUSDoubleClicked(self,index): # ToDo add try except with popup message!
        nodes = list(self.hdfTreeView.selectedH5Nodes())
        if len(nodes) > 0:
            obj = nodes[0]
            if 'NX_class' in obj.attrs:
                try:
                    nxcls = obj.attrs['NX_class'].decode("utf-8")
                except AttributeError:
                    nxcls = obj.attrs['NX_class']

                if nxcls == 'NXentry':
                    if self.bt_autodetect_enable.isChecked():
                        try:    
                            dt = dateparser.parse(obj.h5py_target['start_time'][()])
                        except Exception as e:
                            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot open scan', 
                                'Cannot parse start time of the scan. Please manually select the beamtime id', qt.QMessageBox.Ok, self)
                            msgbox.setDetailedText(traceback.format_exc())
                            clickedbutton = msgbox.exec()
                            return
                            
                        try:
                            btid = backends.getBeamtimeId(dt)
                            self.btid.setCurrentText(btid)
                        except Exception as e:
                            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot open scan', 
                                'Cannot find matching beamtime id to the date %s:\n%s' % (dt,str(e)), qt.QMessageBox.Ok, self)
                            msgbox.setDetailedText(traceback.format_exc())
                            clickedbutton = msgbox.exec()
                            return
                    else:
                        btid = self.btid.currentText()
                    try:
                        ddict = backends.fscans[btid].parse_h5_node(obj)
                    except Exception as e:
                        msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot open scan', 
                            'Cannot parse scan number: %s' % str(e), qt.QMessageBox.Ok, self)
                        msgbox.setDetailedText(traceback.format_exc())
                        clickedbutton = msgbox.exec()
                        return
                    ddict['event'] = "itemDoubleClicked"
                    ddict['file'] = obj.local_filename
                    ddict['node'] = obj
                    ddict['beamtime'] = btid
                    self.pathedit.setText(obj.local_filename)
                    self.scannoBox.setValue(ddict['scanno'])
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

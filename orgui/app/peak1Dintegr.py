# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
# 
# class CurvesROIWidget:
# Copyright (c) 2004-2024 European Synchrotron Radiation Facility
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
__credits__ = []
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
from .. import __version__
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"


import sys
import os
from silx.gui import qt
import warnings
from silx.gui import icons

from io import StringIO
import queue
import threading
import weakref
#from IPython import embed
import silx.gui.plot
from silx.gui.plot import items


import silx
from silx.gui.plot.Profile import ProfileToolBar
from silx.gui.dialog import ImageFileDialog, DataFileDialog
from silx.gui.plot.actions import control as control_actions
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
#from silx.gui.widgets.TableWidget import TableWidget
from silx.gui.plot.CurvesROIWidget import ROITable, ROI, _FloatItem, _RoiMarkerHandler
from silx.utils.weakref import WeakMethodProxy

import traceback

from . import qutils
from .database import DataBase
from .. import resources

import numpy as np
from scipy import special
from ..datautils.xrayutils import HKLVlieg, CTRcalc
from ..datautils.xrayutils import ReciprocalNavigation as rn

from functools import partial

from silx.io.dictdump import dicttonx ,nxtodict, h5todict

import sys

QTVERSION = qt.qVersion()
DEBUG = 0

MAX_ROIS_DISPLAY = 100


class RockingPeakIntegrator(qt.QMainWindow):
    def __init__(self, database, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.database = database
        fontMetric = self.fontMetrics()
        iconSize = qt.QSize(fontMetric.height(), fontMetric.height())
        self.filedialogdir = '.'
        self._currentRoInfo = {}
        self._idx = 0
        
        dbdockwidget = qt.QDockWidget("Integrated data")
        
        self.dbside = qt.QMainWindow()
        
        self._addRoiDialog = ROICreatorDialog(np.array([0,0]), "unknown")
        self._addAllRoiDialog =  IntegrationEstimator(np.array([0,0]), "unknown")
        
        self.integrationCorrection = IntegrationCorrectionsDialog()
        
        
        self.dbview = silx.gui.hdf5.Hdf5TreeView()
        self.dbview.setSortingEnabled(True)
        self.dbview.setModel(self.database.hdf5model)
        self.dbview.setExpandsOnDoubleClick(False)
        self.dbview.addContextMenuCallback(self.database.nexus_treeview_callback)
        
        self.dbside.setCentralWidget(self.dbview)
        
        dbdockwidget.setWidget(self.dbside)
        self.addDockWidget(qt.Qt.RightDockWidgetArea,dbdockwidget)
        
        self.plotROIselect = silx.gui.plot.PlotWindow(
            parent=self,
            backend=None,
            resetzoom=True,
            autoScale=True,
            logScale=True,
            grid=True,
            curveStyle=True,
            colormap=False,
            aspectRatio=False,
            yInverted=False,
            copy=True,
            save=True,
            print_=True,
            control=False,
            position=True,
            roi=False,
            mask=False,
            fit=True
        )
        
        toolbar = qt.QToolBar()
        toolbar.addAction(control_actions.OpenGLAction(parent=toolbar, plot=self.plotROIselect))
        self.plotROIselect.addToolBar(toolbar)
        
        self.curveSlider = qutils.DataRangeSlider()
        
        self.roiwidget = CurvesROIWidget(self, "ROIs", self.plotROIselect)
        self.roiwidget.roiTable.clear()
        
        #self.roiwidget.sigROISignal
        
        #self.plotROIselect.removeToolBar(self.plotROIselect.toolBar())
        #self.plotROIselect.addToolBar(qt.Qt.LeftToolBarArea ,self.plotROIselect.toolBar())
        
        #self.plotROintegrated = silx.gui.plot.Plot1D(self)
        
        self.layout = qt.QVBoxLayout()
        self.mainwidget = qt.QWidget()
        

        
        self.zoomslider = qt.QSlider()
        self.zoomslider.setOrientation(qt.Qt.Horizontal)
        self.zoomslider.setEnabled(True)
        self.zoomslider.setMinimum(0)
        self.zoomslider.setMaximum(100)
        self.zoomslider.setValue(10)
        
        self.zoom_menu = qt.QMenu()
        
        self.zoomwidget = qt.QWidgetAction(self.zoom_menu)
        self.zoomwidget.setDefaultWidget(self.zoomslider)

        self.zoom_menu.addAction(self.zoomwidget)
        
        #self.alpha_btn = qt.QToolButton(resources.getQicon("sum_image.png"),"slider")
        self.zoom_btn = qt.QToolButton()
        self.zoom_btn.setIcon(resources.getQicon("search"))
        self.zoom_btn.setToolTip("set automatic zoom factor")
        self.zoom_btn.setPopupMode(qt.QToolButton.InstantPopup)
        self.zoom_btn.setMenu(self.zoom_menu)
        
        plotselecttools_bar = qt.QHBoxLayout()
        plotselecttools_bar.addWidget(self.curveSlider)
        plotselecttools_bar.addWidget(self.zoom_btn)
        self.autozoom_checkbox = qt.QCheckBox("auto zoom")
        self.autozoom_checkbox.setChecked(True)
        plotselecttools_bar.addWidget(self.autozoom_checkbox)
        
        _plotlayout = qt.QVBoxLayout()
        _plotlayout.addWidget(self.plotROIselect)
        _plotlayout.addLayout(plotselecttools_bar)
        
        self.layout.addLayout(_plotlayout, 2)
        #self.layout.addWidget(self.curveSlider, 1)
        
        
        # ROIS
        
        bottom_layout = qt.QHBoxLayout()
        bottom_layout.addWidget(self.roiwidget)
        
        roi_edit_layout = qt.QVBoxLayout()
        
        anchorROIsGroup = qt.QGroupBox("Anchor ROIs")
        anchorROIsGroupLayout = qt.QHBoxLayout()
        
        #self.anchorROIButton = qt.QPushButton()
        #self.anchorROIButton.setIcon(resources.getQicon("anchor-ROI"))
        #self.anchorROIButton.setIconSize(iconSize)
        #self.anchorROIButton.setCheckable(True)
        #
        #self.anchorROIButton.toggled.connect(self.onAnchorBtnToggled)
        
        self.anchorLoadRoiBtn = qt.QPushButton()
        self.anchorLoadRoiBtn.setIcon(icons.getQIcon("document-open"))
        self.anchorLoadRoiBtn.setIconSize(iconSize)
        self.anchorLoadRoiBtn.clicked.connect(self.onAnchorLoadRoi)
        
        
        self.anchorSaveRoiBtn = qt.QPushButton()
        self.anchorSaveRoiBtn.setIcon(icons.getQIcon("document-save"))
        self.anchorSaveRoiBtn.setIconSize(iconSize)
        self.anchorSaveRoiBtn.clicked.connect(self.onAnchorSaveRoi)
        
        
        self.previousROIButton = qt.QPushButton()
        self.previousROIButton.setIcon(icons.getQIcon("previous"))
        self.previousROIButton.setIconSize(iconSize)
        self.previousROIButton.clicked.connect(self.onToPreviousAnchor)
        self.nextROIButton = qt.QPushButton()
        self.nextROIButton.setIcon(icons.getQIcon("next"))
        self.nextROIButton.setIconSize(iconSize)
        self.nextROIButton.clicked.connect(self.onToNextAnchor)
        
        self.fitAnchorsButton = qt.QPushButton("Fit between anchors")
        self.fitAnchorsButton.clicked.connect(self.fit_anchors_along_rod)
        
        #anchorROIsGroupLayout.addWidget(self.anchorROIButton)
        anchorROIsGroupLayout.addWidget(self.fitAnchorsButton)
        anchorROIsGroupLayout.addStretch()
        anchorROIsGroupLayout.addWidget(self.anchorLoadRoiBtn)
        anchorROIsGroupLayout.addWidget(self.anchorSaveRoiBtn)
        anchorROIsGroupLayout.addWidget(self.previousROIButton)
        anchorROIsGroupLayout.addWidget(self.nextROIButton)
        
        anchorROIsGroup.setLayout(anchorROIsGroupLayout)
        
        roi_edit_layout.addWidget(anchorROIsGroup)
        
        modifyROIsGroup = qt.QGroupBox("Modify ROIs")
        modifyROIsGroupLayout = qt.QHBoxLayout()
        
        self.addAllROIButton = qt.QPushButton("Add All")
        self.addAllROIButton.clicked.connect(self.onAddAllROI)
        
        #self.estimateROIButton.setIcon(icons.getQIcon("previous"))
        #self.estimateROIButton.setIconSize(iconSize)
        
        self.addROIButton = qt.QPushButton("Add")
        self.addROIButton.clicked.connect(self.onAddROI)
        #self.addROIButton.setIcon(icons.getQIcon("previous"))
        #self.addROIButton.setIconSize(iconSize)
        
        self.deleteROIButton = qt.QPushButton("Delete")
        self.deleteROIButton.clicked.connect(self.onDeleteROI)
        #self.deleteROIButton.setIcon(icons.getQIcon("previous"))
        #self.deleteROIButton.setIconSize(iconSize)
        self.deleteAllROIButton = qt.QPushButton("Delete All")
        self.deleteAllROIButton.clicked.connect(self.onDeleteAllROI)
        
        
        modifyROIsGroupLayout.addWidget(self.addAllROIButton)
        modifyROIsGroupLayout.addWidget(self.addROIButton)
        modifyROIsGroupLayout.addWidget(self.deleteROIButton)
        modifyROIsGroupLayout.addWidget(self.deleteAllROIButton)
        
        modifyROIsGroup.setLayout(modifyROIsGroupLayout)
        
        roi_edit_layout.addWidget(modifyROIsGroup)
        
        integrateOptionsGroup = qt.QGroupBox("Integrate options")
        integrateOptionsGroupLayout = qt.QHBoxLayout()
        
        self.lorentzButton = qt.QCheckBox("Lorentz")
        self.footprintButton = qt.QCheckBox("footprint")
        self.footprintOptionsButton = qt.QPushButton("options")
        self.footprintOptionsButton.clicked.connect(self.integrationCorrection.exec)
        
        
        integrateOptionsGroupLayout.addWidget(self.lorentzButton)
        integrateOptionsGroupLayout.addWidget(self.footprintButton)
        integrateOptionsGroupLayout.addWidget(self.footprintOptionsButton)
        
        integrateOptionsGroup.setLayout(integrateOptionsGroupLayout)
        
        roi_edit_layout.addWidget(integrateOptionsGroup)
        
        self.integrateButton = qt.QPushButton("integrate")
        self.integrateButton.clicked.connect(self.integrate)
        
        roi_edit_layout.addWidget(self.integrateButton)
        
        bottom_layout.addLayout(roi_edit_layout)
        
        
        self.layout.addLayout(bottom_layout, 1)
        

        
        #self.zoom_btn_act = qt.QWidgetAction(self)
        #self.zoom_btn_act.setDefaultWidget(self.zoom_btn)
        

        #self.toolbar.addAction(self.showMaxAct)
        #self.toolbar.addAction(self.showSumAct)
        #self.toolbar.addAction(self.alpha_btn_act)
        
        
        
        #self.layout.addWidget(self.plotROintegrated, 1)
        
        self.mainwidget.setLayout(self.layout)
        self.setCentralWidget(self.mainwidget)
        
        toolbar = qt.QToolBar("Database toolbar",self.dbside)
        
        loadDatabaseAct = toolbar.addAction(icons.getQIcon("document-open"),"Open orgui database")
        loadDatabaseAct.triggered.connect(self.database.onOpenDatabase)
        
        savenewact = toolbar.addAction(icons.getQIcon("layer-nx"),"Select orgui database location")
        savenewact.triggered.connect(self.database.onSaveNewDBFile)
                
        saveact = toolbar.addAction(icons.getQIcon("document-save"),"Save orgui database")
        saveact.triggered.connect(self.database.onSaveDBFile)
        
        self.database.sigChangeRockingScan.connect(self.onChangeRockingScan)
        
        self.zoomslider.valueChanged.connect(self.resetXZoomScaled)
        
        self.dbside.addToolBar(toolbar)
        
        #self.curveSlider.setAxis(np.arange(100) / 10., "deg")
        
        self.curveSlider.sigValueChanged.connect(self.onSliderValueChanged)
        #self.roiwidget.sigROISignal.connect(lambda d: print(d))

        
    def set_roscan(self, name):
        ro_info = self.get_rocking_scan_info(name)
        #if ro_info['name'] + '/integration' not in self.database.nxfile:
        #    roi1D_info = self.estimate_roi1D_info(ro_info)
        #    self.database.add_nxdict(roi1D_info, update_mode='modify', h5path=ro_info['name'] + '/integration')
        
        self._currentRoInfo = ro_info
        self._idx = 0

        self.curveSlider.setAxis(self._currentRoInfo['s'], "s")
        self.curveSlider.setIndex(0)
        self.plotRoCurve(0)
        self.plotROIselect.resetZoom()
        if self.autozoom_checkbox.isChecked():
            self.resetXZoomScaled(self.zoomslider.value())
            
    def onAnchorSaveRoi(self):
        if self._currentRoInfo and 'name' in self._currentRoInfo:
            if self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
                if self._currentRoInfo['name'] + '/integration/peakpos' not in self.database.nxfile:
                    qt.QMessageBox.warning(self,'Cannot save roi locations', 'Cannot save ROI locations:\nNo ROI information availabe')
                    return
            else:
                qt.QMessageBox.warning(self,'Cannot save roi locations', 'Cannot save ROI locations:\nNo ROI information availabe')
                return
        else:
            qt.QMessageBox.warning(self,'Cannot save roi locations', 'Cannot save ROI locations:\nNo ROI information availabe')
            return
            

        fileTypeDictSave1D = {"Plain ascii file (*.dat)" : "dat", "CSV file (*.csv)" : "csv",  "NumPy format (*.npy)" : "ndarray"}
        
        fileTypeDictSPEC = {'SPEC file (*.spec)': 'spec'}
        
        fileTypeDict = {**fileTypeDictSave1D, **fileTypeDictSPEC}
        
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Save ROI to file",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.filedialogdir = os.path.splitext(filename)[0]
        
        if filetype in fileTypeDictSave1D:
            fileext = fileTypeDictSave1D[filetype]
            try:
                self.saveROI(filename, fileext)
            except Exception as e:
                qutils.critical_detailed_message(self, 'Cannot save ROIs','Cannot save ROI locations: %s' % e, traceback.format_exc())
                
    def saveROI(self, filename, fileext="dat"):
        
        roi_info = h5todict(self.database.nxfile, self._currentRoInfo['name'] + '/integration/')
        
        fmt = "%.7g"
        csvdelim=";"
            
        data = []
        header = []
        
        header.append("peakpos")
        data.append(roi_info["peakpos"])
        
        for k in list(roi_info.keys()):
            if k.startswith('sig') or k.startswith('bg'):
                header.append(k + '_from')
                data.append(roi_info[k]["from"])
                header.append(k + '_to')
                data.append(roi_info[k]["to"])
                header.append(k + '_anchor')
                data.append(roi_info[k]["anchor"].astype(float))
        
        data = np.vstack(data).T
        header = " ".join(header) 
        
        if fileext == "dat":
            np.savetxt(filename, data, header=header, fmt=fmt)
        elif fileext == "csv":
            np.savetxt(filename, data, header=header, fmt=fmt, delimiter=csvdelim)
        elif fileext == "ndarray":
            np.save(filename, data)
        else:
            raise Exception("No supported file type %s" % fileext)
        
        
    
    #def onAnchorLoadRoi(self):
        #datasetDialog = DataFileDialog.DataFileDialog(self)
        #datasetDialog.setFilterMode(DataFileDialog.DataFileDialog.FilterMode.ExistingDataset)
        #
        #def customFilter(obj):
        #    print(obj.basename)
        #    return True
        #    
        #    if "NX_class" in obj.attrs:
        #        if 'orgui_meta' in obj.attrs and obj.attrs['orgui_meta'] == 'rocking':
        #            if 'integration' in obj:
        #                return True
        #    return False
        #
        #if datasetDialog.exec():
        #    print(datasetDialog.selectedUrl())
            
    def onAnchorLoadRoi(self):
        if self._currentRoInfo:
            fileTypeDictSave1D = {"Plain ascii file (*.dat)" : "dat", "CSV file (*.csv)" : "csv"}
            
            fileTypeDictSPEC = {'SPEC file (*.spec)': 'spec'}
            
            fileTypeDict = {**fileTypeDictSave1D, **fileTypeDictSPEC}
            
            fileTypeFilter = ""
            for f in fileTypeDict:
                fileTypeFilter += f + ";;"

            filename, filetype = qt.QFileDialog.getOpenFileName(self,"Load ROI from file",
                                                      self.filedialogdir,
                                                      fileTypeFilter[:-2])
            if filename == '':
                return
            self.filedialogdir = os.path.splitext(filename)[0]
            
            try:
                roi_info = self.loadROIinfo(filename)
            except Exception as e:
                qutils.critical_detailed_message(self, 'Cannot load ROIs','Cannot load ROI locations: %s' % e, traceback.format_exc())
                return
                
            try:
                self.setROIinfo(roi_info)
            except Exception as e:
                qutils.critical_detailed_message(self, 'Cannot load ROIs','Cannot apply ROI locations to rocking scan: %s' % e, traceback.format_exc())
                return
        
    def setROIinfo(self, roi_info, interpolate=True, pktol=0.1):
        if self._currentRoInfo:
            ro_info = self.get_rocking_scan_info(self._currentRoInfo['name'])
            sc_h5 = self.database.nxfile[ro_info['name']]['rois']
            if ro_info['axisname'] == 'mu':
                pk_pos_exact = sc_h5['alpha_pk'][()]
            elif ro_info['axisname'] == 'th':
                pk_pos_exact = sc_h5['theta_pk'][()]
            elif ro_info['axisname'] == 'chi':
                pk_pos_exact = sc_h5['chi_pk'][()]
            elif ro_info['axisname'] == 'phi':
                pk_pos_exact = sc_h5['phi_pk'][()]
            else:
                raise ValueError("Cannot estimate peak position: unknown scan axis %s" % ro_info['axisname'])
            
            def _set_roi_info(roi_info):
                if self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
                    del self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
                
                for k in roi_info:
                    if k.startswith(('sig', 'bg')):
                        roi_info[k]["@NX_class"] = u"NXcollection"
                        roi_info[k]['from'] = roi_info[k]['from'].astype(np.float64)
                        roi_info[k]['to'] = roi_info[k]['to'].astype(np.float64)
                        roi_info[k]['anchor'] = roi_info[k]['anchor'].astype(bool)
                        
                roi_info["peakpos"] = roi_info["peakpos"].astype(np.float64)
                roi_info["@NX_class"] = u"NXcollection"
                roi_info["@info"] = u"ROI information for the rocking scan integration"
                
                self.database.add_nxdict(roi_info, update_mode='modify', h5path=ro_info['name'] + '/integration')
                self.plotRoCurve(self._idx)
                
            pk_pos = roi_info["peakpos"]
            
            if pk_pos.size != pk_pos_exact.size:
                if not interpolate:
                    raise ValueError('Scan length mismatch: requires interpolation is: %s, supplied roi_info: %s' % (pk_pos_exact.size, pk_pos.size))
            elif not np.all(np.isclose(pk_pos, pk_pos_exact)):
                if not interpolate:
                    raise ValueError('Peak position mismatch: requires interpolation')
            else:
                _set_roi_info(roi_info)
                return
                
            
            if np.abs(np.amax(pk_pos) - np.amax(pk_pos_exact)) >= pktol:
                raise ValueError('peak position mismatch tolerance exceeded at max')
            
            if np.abs(np.amin(pk_pos) - np.amin(pk_pos_exact)) >= pktol:
                raise ValueError('peak position mismatch tolerance exceeded at min')
            
            roi_info_interp = dict()
            
            roi_info_interp["peakpos"] = pk_pos_exact
            for k in roi_info:
                if k.startswith(('sig', 'bg')):
                    roi_info_interp[k] = dict()
                    interfrom = interp1d(roi_info["peakpos"],roi_info[k]['from'],fill_value="extrapolate")
                    roi_info_interp[k]['from'] = interfrom(pk_pos_exact)
                    
                    interto = interp1d(roi_info["peakpos"],roi_info[k]['to'],fill_value="extrapolate")
                    roi_info_interp[k]['to'] = interto(pk_pos_exact)
                    
                    interanchor = interp1d(roi_info["peakpos"],roi_info[k]['anchor'], kind='nearest',fill_value="extrapolate")
                    roi_info_interp[k]['anchor'] = interanchor(pk_pos_exact).astype(bool)
            
            _set_roi_info(roi_info_interp)
            
            
    def loadROIinfo(self, filename):
        fname, fileext = os.path.splitext(filename)
        if fileext == ".dat":
            data = np.genfromtxt(filename, names=True)
            roi_names = []
            for k in data.dtype.names:
                if k.startswith(('sig', 'bg')):
                    k_sp = k.split('_')
                    no = int(k_sp[1])
                    rname = k_sp[0] + '_' + k_sp[1]
                    if rname not in roi_names:
                        roi_names.append(rname)
            
            roi_info = dict()
            roi_info["peakpos"] = data["peakpos"]
            for k in roi_names:
                roi_info[k] = dict()
                roi_info[k]['from'] = data[k + '_from']
                roi_info[k]['to'] = data[k + '_to']
                roi_info[k]['anchor'] = data[k + '_anchor'].astype(bool)
        else:
            raise Exception("Not supported file type %s" % fileext)
        return roi_info
            
            
        
            

    #ddict = {
    #        "event": "indexChanged",
    #        "oldtxt": self._lineTxt,
    #        "newtxt": txt,
    #        "idx" : idx,
    #        "value" : self._data[idx],
    #        "id": id(self),
    #    }
    
    def onToNextAnchor(self):
        if self._currentRoInfo and self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
            roi_info = h5todict(self.database.nxfile, self._currentRoInfo['name'] + '/integration/')
            anchors_all = []
            for roikey in roi_info:
                if roikey.startswith('sig') or roikey.startswith('bg'):
                    rdict = roi_info[roikey]
                    anchors = np.nonzero(rdict['anchor'])[0]
                    anchors_all.append(anchors)
            anchors_all = np.sort(np.unique(np.concatenate(anchors_all)))
            
            anchors_next = anchors_all[anchors_all > self._idx]
            if anchors_next.size > 0:
                self.plotRoCurve(anchors_next[0])
            
    def onToPreviousAnchor(self):
        if self._currentRoInfo and self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
            roi_info = h5todict(self.database.nxfile, self._currentRoInfo['name'] + '/integration/')
            anchors_all = []
            for roikey in roi_info:
                if roikey.startswith('sig') or roikey.startswith('bg'):
                    rdict = roi_info[roikey]
                    anchors = np.nonzero(rdict['anchor'])[0]
                    anchors_all.append(anchors)
            anchors_all = np.sort(np.unique(np.concatenate(anchors_all)))
            
            anchors_prev = anchors_all[anchors_all < self._idx]
            if anchors_prev.size > 0:
                self.plotRoCurve(anchors_prev[-1])
        
    
    
    #def get_anchors(self):
        
    
    def fit_anchors_along_rod(self):
        if self._currentRoInfo and self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
            roih5grp = self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
            roi_info = h5todict(self.database.nxfile, self._currentRoInfo['name'] + '/integration/')
            pk_pos_exact = roi_info["peakpos"]
            axis = self._currentRoInfo['axis']
            for roikey in roi_info:
                if roikey.startswith('sig') or roikey.startswith('bg'):
                    rdict = roi_info[roikey]
                    anchors = np.nonzero(rdict['anchor'])[0]
                    if anchors.size == 0:
                        continue
                    from_ar_anchors = rdict['from'][anchors] - pk_pos_exact[anchors]
                    to_ar_anchors = rdict['to'][anchors] - pk_pos_exact[anchors]
                    
                    from_ar = np.zeros(rdict['from'].size)
                    to_ar = np.zeros(rdict['to'].size)
                    
                    # constant to first anchor 
                    idx_prev = 0
                    from_ar[:anchors[idx_prev]] = pk_pos_exact[:anchors[idx_prev]] + from_ar_anchors[idx_prev]
                    to_ar[:anchors[idx_prev]] = pk_pos_exact[:anchors[idx_prev]] + to_ar_anchors[idx_prev]
                    if anchors.size > 0:
                        for idx in range(1,anchors.size):
                            slope = (from_ar_anchors[idx] - from_ar_anchors[idx_prev]) / ( axis[anchors[idx]] - axis[anchors[idx_prev]] )
                            from_tmp = slope * (axis[anchors[idx_prev]: anchors[idx]] - axis[anchors[idx_prev]] ) 
                            from_ar[anchors[idx_prev]: anchors[idx]] = (from_tmp + from_ar_anchors[idx_prev]) + pk_pos_exact[anchors[idx_prev]: anchors[idx]]
                            
                            slope = (to_ar_anchors[idx] - to_ar_anchors[idx_prev]) / ( axis[anchors[idx]] - axis[anchors[idx_prev]] )
                            to_tmp = slope * (axis[anchors[idx_prev]: anchors[idx]] - axis[anchors[idx_prev]] ) 
                            to_ar[anchors[idx_prev]: anchors[idx]] = (to_tmp + to_ar_anchors[idx_prev]) + pk_pos_exact[anchors[idx_prev]: anchors[idx]]
                            
                            idx_prev = idx
                    
                    #constant from last anchor
                    from_ar[anchors[idx_prev]:] = pk_pos_exact[anchors[idx_prev]:] + from_ar_anchors[idx_prev]
                    to_ar[anchors[idx_prev]:] = pk_pos_exact[anchors[idx_prev]:] + to_ar_anchors[idx_prev]
                    
                    roih5grp[roikey]['from'][:] = from_ar
                    roih5grp[roikey]['to'][:] = to_ar
            self.plotRoCurve(self._idx)
        
    
    def integrate(self):
        if not self._currentRoInfo:
            return
        curves = self.get_all_ro_curves()
        name = self._currentRoInfo['name']
        h5_obj = self.database.nxfile[name]
        cnters = h5_obj["rois"]
        
        s_array = cnters['s'][()]
        alpha = np.deg2rad(cnters['alpha'][()])
        #alpha_pk = cnters['alpha_pk'][()]
        delta = np.deg2rad(cnters['delta'][()])
        gamma = np.deg2rad(cnters['gamma'][()])
        
        axis = curves['axis']
        if self.lorentzButton.isChecked():
            if curves['axisname'] == 'mu':
                C_Lor = 1/np.sin(2*alpha)
                C_rod = np.cos(gamma)
            elif curves['axisname'] == 'th':
                C_Lor = 1 /(np.sin(delta) * np.cos(alpha) * np.cos(gamma))
                C_rod = np.cos(gamma)
            else:
                raise NotImplementedError()
        else:
            C_Lor = 1.
            C_rod = 1.
        
        if self.footprintButton.isChecked():
            
            def total_flux_sample(alpha_i, L, sigma):
                arg = ((L * np.sin(alpha_i)) / (np.sqrt(2) * sigma)) * 0.5
                return (1 / 2) * (special.erf(arg) - special.erf(-arg))
                
            L = self.integrationCorrection.L.value() * 1e-3 # sample size (mm)
            beamsize = self.integrationCorrection.beam.value() * 1e-6 # FWHM (mum)
            
            sigma_beam = beamsize / (2*np.sqrt(2 * np.log(2)))
            
            C_flux_on_sample = total_flux_sample(alpha, L, sigma_beam)
            
            C_illum_area = (np.sqrt(2 * np.pi)* sigma_beam * C_flux_on_sample) / (L * np.sin(alpha))
            
        else:
            C_flux_on_sample = 1.
            C_illum_area = 1.
            
            
        
            

        roih5grp = self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
        roi_info = h5todict(self.database.nxfile, self._currentRoInfo['name'] + '/integration/')
        
        int_data = {}
        for roikey in roi_info:
            if roikey.startswith('sig') or roikey.startswith('bg'):
                int_data[roikey] = {
                    'cnts' : [],
                    'cnts_errors' : [],
                    'raw_cnts' : [],
                    'raw_cnts_errors' : [],
                    'int_interval' : [],
                    'C_Lor' : [],
                    'C_rod' : [],
                    'C_flux_on_sample' : [],
                    'C_illum_area' : []
                }
                
        deltaaxis = np.gradient(axis)
        
        progress = qt.QProgressDialog("Integrating rocking scans","abort",0,s_array.size,self)
        progress.setWindowModality(qt.Qt.WindowModal)
        
        for i, s in enumerate(s_array):
            croibg = curves['croibg'][i]
            croibg_errors = curves['croibg_errors'][i]
            
            for roikey in int_data:
                roi = roi_info[roikey]
                
                idx_from = np.argmin(np.abs(axis - roi['from'][i]))
                idx_to = np.argmin(np.abs(axis - roi['to'][i]))
                if idx_from > idx_to:
                    idx_from, idx_to = idx_to, idx_from
                int_interval = abs(axis[idx_to] - axis[idx_from])
                int_data[roikey]['int_interval'].append(int_interval)
                
                cnts = croibg[idx_from:idx_to]
                cnts_errors = croibg_errors[idx_from:idx_to]
                
                C_corr = np.ones(cnts.size,dtype=float)
                if self.lorentzButton.isChecked():
                    int_data[roikey]['C_Lor'].append( np.mean(C_Lor[i][idx_from:idx_to]) )
                    int_data[roikey]['C_rod'].append( np.mean(C_rod[i][idx_from:idx_to]) )
                
                if self.footprintButton.isChecked():
                    int_data[roikey]['C_flux_on_sample'].append( np.mean(C_flux_on_sample[i][idx_from:idx_to]) )
                    int_data[roikey]['C_illum_area'].append( np.mean(C_illum_area[i][idx_from:idx_to]) )
                    C_corr *= C_flux_on_sample[i][idx_from:idx_to] * C_illum_area[i][idx_from:idx_to]
                    
                I_raw = np.trapz(cnts, axis[idx_from:idx_to])
                I_corr = np.trapz(cnts / C_corr , axis[idx_from:idx_to])
                
                I_raw_error = np.sum(cnts_errors * deltaaxis[idx_from:idx_to])
                I_corr_error = np.sum((cnts_errors / C_corr) * deltaaxis[idx_from:idx_to])
                
                int_data[roikey]['raw_cnts'].append(I_raw)
                int_data[roikey]['raw_cnts_errors'].append(I_raw_error)
                
                int_data[roikey]['cnts'].append(I_corr)
                int_data[roikey]['cnts_errors'].append(I_corr_error)
            
            progress.setValue(i)
            if progress.wasCanceled():
                break
        progress.setValue(s_array.size)
        
        for roikey in int_data:
            for d in list(int_data[roikey].keys()):
                int_data[roikey][d] = np.array(int_data[roikey][d])
        
        # signals:
        
        croi = np.zeros(s_array.size,dtype=float)
        croi_errors = np.zeros(s_array.size,dtype=float)
        raw_croi = np.zeros(s_array.size,dtype=float)
        raw_croi_errors = np.zeros(s_array.size,dtype=float)
        sig_interval = np.zeros(s_array.size,dtype=float)
        C_Lorentz = np.zeros(s_array.size,dtype=float)
        C_rod_intersect = np.zeros(s_array.size,dtype=float)
        
        for roikey in int_data:
            if roikey.startswith('sig'):
                sig_interval += int_data[roikey]['int_interval']
                
        for roikey in int_data:
            if roikey.startswith('sig'):
                croi += int_data[roikey]['cnts']
                croi_errors += int_data[roikey]['cnts_errors']**2
                raw_croi += int_data[roikey]['raw_cnts']
                raw_croi_errors += int_data[roikey]['raw_cnts_errors']**2
                if self.lorentzButton.isChecked():
                    C_Lorentz += int_data[roikey]['C_Lor'] * (int_data[roikey]['int_interval'] / sig_interval)
                    C_rod_intersect += int_data[roikey]['C_rod'] * (int_data[roikey]['int_interval'] / sig_interval)
        
        raw_croi_errors = np.sqrt(raw_croi_errors)
        croi_errors = np.sqrt(croi_errors)
        
        
        bgroi = np.zeros(s_array.size,dtype=float)
        bgroi_errors = np.zeros(s_array.size,dtype=float)
        raw_bgroi = np.zeros(s_array.size,dtype=float)
        raw_bgroi_errors = np.zeros(s_array.size,dtype=float)
        bg_interval = np.zeros(s_array.size,dtype=float)
        
        for roikey in int_data:
            if roikey.startswith('bg'):
                bg_interval += int_data[roikey]['int_interval']
                
        for roikey in int_data:
            if roikey.startswith('bg'):
                ratio = (sig_interval / bg_interval) * (int_data[roikey]['int_interval'] / bg_interval)
                bgroi += int_data[roikey]['cnts'] * ratio
                bgroi_errors += (int_data[roikey]['cnts_errors'] * ratio)**2 # should improve error propagation here!
                raw_bgroi += int_data[roikey]['raw_cnts'] * ratio
                raw_bgroi_errors += (int_data[roikey]['raw_cnts_errors'] * ratio)**2
        
        raw_bgroi_errors = np.sqrt(raw_croi_errors)
        bgroi_errors = np.sqrt(croi_errors)
        
        croibg = croi - bgroi # already normalized bg
        croibg_errors = np.sqrt(croi_errors**2 + bgroi_errors**2)
                
        raw_croibg = raw_croi - raw_bgroi # already normalized bg
        raw_croibg_errors = np.sqrt(raw_croi_errors**2 + raw_bgroi_errors**2)
        
        if self.lorentzButton.isChecked():
            F2_hkl = croibg / (C_Lorentz * C_rod_intersect)
            F2_hkl_errors = raw_croibg_errors / (C_Lorentz * C_rod_intersect)
        
        
        int_data["@NX_class"] = u"NXdetector"
        
        H_1 = cnters['H_1'][0][()]
        H_0 = cnters['H_0'][0][()]

        traj1 = {
            "@NX_class": u"NXcollection",
            "@direction" : u" Integrated rocking scan along H_1*s + H_0 in reciprocal space",
            "H_1"  : H_1,
            "H_0" : H_0,
            "s" : s_array
        }

        suffix = ''
        i = 0
        name1 = str(H_1) + "*s1+" + str(H_0)
        while(self._currentRoInfo['name'] + "/measurement/" + name1 + suffix in self.database.nxfile):
            suffix = "_%s" % i
            i += 1
        availname1 = name1 + suffix
                    
        datas1 = {
            "@NX_class": u"NXdata",
            "sixc_angles": {
                "@NX_class": u"NXpositioner",
                "alpha" : cnters['alpha_pk'][()],
                "omega" :  cnters['omega_pk'][()],
                "theta" :  cnters['theta_pk'][()],
                "delta" : cnters['delta_pk'][()],
                "gamma" :  cnters['gamma_pk'][()],
                "chi" :  cnters['chi_pk'][()],
                "phi" :  cnters['phi_pk'][()],
                "@unit" : u"deg"
            },
            "hkl": {
                "@NX_class": u"NXcollection",
                "h" :  cnters['HKL_pk'][:, 0][()],
                "k" :  cnters['HKL_pk'][:, 1][()],
                "l" : cnters['HKL_pk'][:, 2][()]
            },
            "counters":{
                "@NX_class": u"NXdetector",
                "croibg"  : croibg,
                "croibg_errors" :  croibg_errors,
                "croi" :  croi,
                "croi_errors" :  croi_errors,
                "bgroi"  : bgroi,
                "bgroi_errors"  : bgroi_errors,
                "raw_croibg"  : raw_croibg,
                "raw_croibg_errors" :  raw_croibg_errors,
                "raw_croi" :  raw_croi,
                "raw_croi_errors" :  raw_croi_errors,
                "raw_bgroi"  : raw_bgroi,
                "raw_bgroi_errors"  : raw_bgroi_errors,
                "integrated": int_data
            },
            "pixelcoord": {
                "@NX_class": u"NXdetector",
                "x" : cnters['x'][:, 0][()],
                "y"  : cnters['x'][:, 1][()]
            },
            "trajectory" : traj1,
            "@signal" : u"counters/croibg",
            "@axes": u"trajectory/s",
            "@title": self._currentRoInfo['name'] + "_" + availname1,
            "@orgui_meta": u"roi"
        }
        
        measurement =  {
            "@NX_class": u"NXentry",
            "@default": availname1 ,
            availname1 : datas1
        }
        
        if self.lorentzButton.isChecked():
            measurement[availname1]['counters']['F2_hkl'] = F2_hkl
            measurement[availname1]['counters']['F2_hkl_errors'] = F2_hkl_errors
            measurement[availname1]['@signal'] = u"counters/F2_hkl"

        self.database.add_nxdict(measurement, update_mode='modify', h5path=self._currentRoInfo['name'] + "/measurement")
        
    
    def onSliderValueChanged(self, ddict):
        self.plotRoCurve(ddict['idx'])
        
    def onRoiChanged(self, roi):
        #roi_dict = self.get_roi1D_info(self._idx)
        roi_t = self.roiwidget.roiTable.roidict[roi]
        new_roi_dict = roi_t.toDict()
        if self._currentRoInfo:
            #from IPython import embed; embed()
            if self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
                h5grp = self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
                if roi in h5grp:
                    if h5grp[roi]['from'][self._idx] != new_roi_dict['from'] or h5grp[roi]['to'][self._idx] != new_roi_dict['to']:
                        h5grp[roi]['from'][self._idx] = new_roi_dict['from']
                        h5grp[roi]['to'][self._idx] = new_roi_dict['to']
                        h5grp[roi]['anchor'][self._idx] = True
                        self.roiwidget.roiTable.isModelSetting = True
                        roi_t.setAnchor(True)
                        self.roiwidget.roiTable._updateRoiInfo(roi_t.getID())
                        self.roiwidget.roiTable.isModelSetting = False
                        #print('set true: %s -> %s' % (h5grp[roi]['anchor'][self._idx], roi_t.isAnchor()))
                    else:
                        if h5grp[roi]['anchor'][self._idx] != roi_t.isAnchor(): # toggled anchor 
                            #print('toggled: %s -> %s' % (h5grp[roi]['anchor'][self._idx], roi_t.isAnchor()))
                            h5grp[roi]['anchor'][self._idx] = roi_t.isAnchor()
                    #else:
                    #    print('New Anchor set: True')
                    #    h5grp[roi]['anchor'][self._idx] = True
                    #    roi_t.setAnchor(True)
                        #roi_t.anchor_updating = False
                    
        #with qt.QSignalBlocker(self.anchorROIButton):
        #    self.anchorROIButton.setChecked(True)
        
    def get_ro_curve(self, idx):
        name = self._currentRoInfo['name']
        h5_obj = self.database.nxfile[name]
        cnters = h5_obj["rois"]
        curve = {
            'axisname' : self._currentRoInfo['axisname'],
            'axis' : self._currentRoInfo['axis'],
            'croibg' : cnters['croibg'][idx][()],
            'croibg_errors' : cnters['croibg_errors'][idx][()]
        }
        return curve
        
    def get_all_ro_curves(self):
        name = self._currentRoInfo['name']
        h5_obj = self.database.nxfile[name]
        cnters = h5_obj["rois"]
        curve = {
            'axisname' : self._currentRoInfo['axisname'],
            'axis' : self._currentRoInfo['axis'],
            'croibg' : cnters['croibg'][()],
            'croibg_errors' : cnters['croibg_errors'][()]
        }
        return curve
        
    #def onAnchorBtnToggled(self, state):
    #    if self._currentRoInfo:
    #        if self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
    #            h5grp = self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
    #            for k in h5grp:
    #                if k.startswith('sig') or k.startswith('bg'):
    #                    h5grp[k]['anchor'][self._idx] = state
    #                
    #    print(state)
        
        
    def onAddAllROI(self):
        if self._currentRoInfo:
            try:
                self._addAllRoiDialog.set_axis(self._currentRoInfo['axis'], self._currentRoInfo['axisname'])
                if self._addAllRoiDialog.exec() == qt.QDialog.Accepted:
                    roidict = self._addAllRoiDialog.get_ranges()
                    
                    if self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
                        del self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
                        
                    ro_info = self.get_rocking_scan_info(self._currentRoInfo['name'])
                    sc_h5 = self.database.nxfile[ro_info['name']]['rois']
                    if ro_info['axisname'] == 'mu':
                        pk_pos_exact = sc_h5['alpha_pk'][()]
                    elif ro_info['axisname'] == 'th':
                        pk_pos_exact = sc_h5['theta_pk'][()]
                    elif ro_info['axisname'] == 'chi':
                        pk_pos_exact = sc_h5['chi_pk'][()]
                    elif ro_info['axisname'] == 'phi':
                        pk_pos_exact = sc_h5['phi_pk'][()]
                    else:
                        raise ValueError("Cannot estimate peak position: unknown scan axis %s" % ro_info['axisname'])
                
                        
                    roi1Ddict = {
                        "@NX_class": u"NXcollection",
                        "@info" : u"ROI information for the rocking scan integration",
                        "peakpos" : pk_pos_exact,
                        'sig_1' : {
                                "@NX_class": u"NXcollection",
                                'from' :  pk_pos_exact + roidict['sig_1']['from'],
                                'to' : pk_pos_exact + roidict['sig_1']['to'],
                                'anchor' : np.zeros_like(pk_pos_exact, dtype=bool)
                            },
                        'bg_1' : {
                                "@NX_class": u"NXcollection",
                                'from' :  pk_pos_exact + roidict['bg_1']['from'],
                                'to' : pk_pos_exact + roidict['bg_1']['to'],
                                'anchor' : np.zeros_like(pk_pos_exact, dtype=bool)
                            },
                        'bg_2' : {
                                "@NX_class": u"NXcollection",
                                'from' :  pk_pos_exact + roidict['bg_2']['from'],
                                'to' : pk_pos_exact + roidict['bg_2']['to'],
                                'anchor' : np.zeros_like(pk_pos_exact, dtype=bool)
                            }
                    }
                    
                    self.database.add_nxdict(roi1Ddict, update_mode='modify', h5path=ro_info['name'] + '/integration')
                    self.plotRoCurve(self._idx)
            except Exception as e:
                qutils.warning_detailed_message(self, "Cannot add ROI", "Cannot add ROI: %s" % e, traceback.format_exc())
        
    def onAddROI(self):
        if self._currentRoInfo:
            try:
                self._addRoiDialog.set_axis(self._currentRoInfo['axis'], self._currentRoInfo['axisname'])
                if self._addRoiDialog.exec() == qt.QDialog.Accepted:
                    roidict = self._addRoiDialog.get_roi()
                    ro_info = self.get_rocking_scan_info(self._currentRoInfo['name'])
                    sc_h5 = self.database.nxfile[ro_info['name']]['rois']
                    if ro_info['axisname'] == 'mu':
                        pk_pos_exact = sc_h5['alpha_pk'][()]
                    elif ro_info['axisname'] == 'th':
                        pk_pos_exact = sc_h5['theta_pk'][()]
                    elif ro_info['axisname'] == 'chi':
                        pk_pos_exact = sc_h5['chi_pk'][()]
                    elif ro_info['axisname'] == 'phi':
                        pk_pos_exact = sc_h5['phi_pk'][()]
                    else:
                        raise ValueError("Cannot estimate peak position: unknown scan axis %s" % ro_info['axisname'])
                
                    roi_dict = self.get_roi1D_info(self._idx) # check existing ROIS
                    
                    prefix = 'sig' if roidict['signal'] else 'bg'
                    no = 0
                    for k in roi_dict:
                        if k.startswith(prefix):
                            no = max(no, int(k.split('_')[1]))
                            
                    no += 1
                            
                    roi1Ddict = {
                        "@NX_class": u"NXcollection",
                        "@info" : u"ROI information for the rocking scan integration",
                        "peakpos" : pk_pos_exact,
                        '%s_%s' % (prefix, no) : {
                                "@NX_class": u"NXcollection",
                                'from' :  pk_pos_exact + roidict['from'],
                                'to' : pk_pos_exact + roidict['to'],
                                'anchor' : np.zeros_like(pk_pos_exact, dtype=bool)
                        }
                    }
                    self.database.add_nxdict(roi1Ddict, update_mode='modify', h5path=ro_info['name'] + '/integration')
                    self.plotRoCurve(self._idx)
            except Exception as e:
                qutils.warning_detailed_message(self, "Cannot add ROI", "Cannot add ROI: %s" % e, traceback.format_exc())
    
    def onDeleteROI(self):
        if self._currentRoInfo:
            roi_name = self.roiwidget.currentRoi.getName()
            if qt.QMessageBox.Yes == qt.QMessageBox.question(self, "Delete ROI", "Are you sure you want to delete ROI %s?" % roi_name):
                if self._currentRoInfo['name'] + '/integration/' + roi_name in self.database.nxfile:
                    del self.database.nxfile[self._currentRoInfo['name'] + '/integration/' + roi_name]
                    self.plotRoCurve(self._idx)
                else:
                    qutils.warning_detailed_message(self, "Cannot delete ROI", "Cannot delete ROI: no such ROI in database", "")
    
    def onDeleteAllROI(self):
        if self._currentRoInfo:
            if qt.QMessageBox.Yes == qt.QMessageBox.question(self, "Delete All ROIs", "Are you sure you want to delete all ROIs?"):
                if self._currentRoInfo['name'] + '/integration/'in self.database.nxfile:
                    del self.database.nxfile[self._currentRoInfo['name'] + '/integration/']
                    self.plotRoCurve(self._idx)
                else:
                    qutils.warning_detailed_message(self, "Cannot delete ROIs", "Cannot delete ROI: no ROIs in database", "")            
            
        
    def plotRoCurve(self, idx):
        ro_curve = self.get_ro_curve(idx)
        s = self._currentRoInfo['s'][idx]
        hkl = self._currentRoInfo['H_1'] * s + self._currentRoInfo['H_0']
        title = "Rocking scan at s = %s, HKL = [%.2f %.2f %.2f]" % (s, *hkl)
        #print('before table clear')
        self.roiwidget.roiTable.clear()
        #print('after table clear')
        
        #print('before plot clear')
        self.plotROIselect.clear()
        #print('after plot clear')
        
        lbl = self.plotROIselect.addCurve(ro_curve['axis'], ro_curve['croibg'], legend=title, 
                           xlabel="%s / deg" % self._currentRoInfo['axisname'], 
                           ylabel='center ROI background subtracted / cnts', 
                           yerror=ro_curve['croibg_errors'], resetzoom=False)
        
        self.plotROIselect.setActiveCurve(lbl)
        self._idx = idx
        
        with qt.QSignalBlocker(self.curveSlider):
            self.curveSlider.setIndex(idx)

        
        roi_dict = self.get_roi1D_info(idx)
        #with qt.QSignalBlocker(self.anchorROIButton):
        #    for k in roi_dict:
        #        dk = roi_dict[k]
        #        if dk['anchor']:
        #            self.anchorROIButton.setChecked(True)
        #            break
        #    else:
        #        self.anchorROIButton.setChecked(False)
            
        minfrom = np.inf
        maxto = -np.inf
        for key in roi_dict:
            if key.startswith('sig'):
                roi_d = roi_dict[key]
                roi = AnchorROI(key, fromdata=roi_d['from'], todata=roi_d['to'], type_=str(roi_d['type']), anchor=roi_d['anchor'])
                minfrom = min(minfrom, roi_d['from'])
                maxto = max(maxto, roi_d['to'])
                self.roiwidget.roiTable.isModelSetting = True # workaround for a signal loop, which resets anchor status to False
                # took me 2 days to concede, do not try to fix this bug. Workaround is enough
                self.roiwidget.roiTable.addRoi(roi)
                self.roiwidget.roiTable._updateRoiInfo(roi.getID())
                roi.sigChanged.connect(partial(self.onRoiChanged, key))
                self.roiwidget.roiTable.isModelSetting = False
            elif key.startswith('bg'): # color?
                roi_d = roi_dict[key]
                roi = AnchorROI(key, fromdata=roi_d['from'], todata=roi_d['to'], type_=str(roi_d['type']), anchor=roi_d['anchor'])
                minfrom = min(minfrom, roi_d['from'])
                maxto = max(maxto, roi_d['to'])
                self.roiwidget.roiTable.isModelSetting = True # workaround for a signal loop, which resets anchor status to False
                # took me 2 days to concede, do not try to fix this bug. Workaround is enough
                self.roiwidget.roiTable.addRoi(roi)
                self.roiwidget.roiTable._updateRoiInfo(roi.getID())
                roi.sigChanged.connect(partial(self.onRoiChanged, key))
                self.roiwidget.roiTable.isModelSetting = False
        
        # toDo: set from GUI
        if self.autozoom_checkbox.isChecked():
            if np.all(np.isfinite([minfrom, maxto])):
                self.resetXZoomScaled(self.zoomslider.value())
    
    def resetXZoomScaled(self, value):
        try:
            roi_dict = self.get_roi1D_info(self._idx)
            minfrom = np.inf
            maxto = -np.inf
            for key in roi_dict:
                if key.startswith('sig'):
                    roi_d = roi_dict[key]
                    minfrom = min(minfrom, roi_d['from'])
                    maxto = max(maxto, roi_d['to'])
                elif key.startswith('bg'): # color?
                    roi_d = roi_dict[key]
                    minfrom = min(minfrom, roi_d['from'])
                    maxto = max(maxto, roi_d['to'])
            if np.all(np.isfinite([minfrom, maxto])):
                add_rnge = (maxto - minfrom) * ((value**1.5)/100)
                self.plotROIselect.setGraphXLimits(minfrom - add_rnge, maxto + add_rnge)
        except Exception as e:
            pass
            #print("Cannot zoom:" , e)
            # can fail, add better error handling later!
        
                
                
    
    def get_roi1D_info(self, idx):
        name = self._currentRoInfo['name']
        if 'integration' in self.database.nxfile[name]:
            h5_obj = self.database.nxfile[name]['integration']
            
            roi1Ddict = {}
            for k in h5_obj:
                if k.startswith('sig') or k.startswith('bg') :
                    roi1Ddict[k] = {
                        'name' : k,
                        'from' : h5_obj[k]['from'][idx][()],
                        'to' : h5_obj[k]['to'][idx][()],
                        'type' : 'signal',
                        'anchor' : h5_obj[k]['anchor'][idx][()]
                    }
        else:
            roi1Ddict = {}
        return roi1Ddict

    #def estimate_roi1D_info(self, ro_info):
    #    # ToDo make settings available from GUI 
    #    
    #    # estimate peak pos
    #    axis_pk_pos = [] 
    #    axis_pk_pos_idx = [] 
    #
    #    sc_h5 = self.database.nxfile[ro_info['name']]['rois']
    #    if ro_info['axisname'] == 'mu':
    #        pk_pos_exact = sc_h5['alpha_pk'][()]
    #    elif ro_info['axisname'] == 'th':
    #        pk_pos_exact = sc_h5['theta_pk'][()]
    #    elif ro_info['axisname'] == 'chi':
    #        pk_pos_exact = sc_h5['chi_pk'][()]
    #    elif ro_info['axisname'] == 'phi':
    #        pk_pos_exact = sc_h5['phi_pk'][()]
    #    else:
    #        raise ValueError("Cannot estimate peak position: unknown scan axis %s" % ro_info['axisname'])
    #    #idx_pk = np.argmin(np.abs(sc_h5['axis'][()] - pk_pos_exact), axis=1)
    #
    #    # ToDo make settings available from GUI 
    #    sig_width = 0.1 # deg
    #    bg_1 = -0.2 , -0.1
    #    bg_2 = 0.1 , 0.2
    #    
    #    roi1Ddict = {
    #        "@NX_class": u"NXcollection",
    #        "@info" : u"ROI information for the rocking scan integration",
    #        "peakpos" : pk_pos_exact,
    #        'sig_1' : {
    #                "@NX_class": u"NXcollection",
    #                'from' :  pk_pos_exact - 0.5*sig_width,
    #                'to' : pk_pos_exact + 0.5*sig_width,
    #                'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
    #            },
    #        'bg_1' : {
    #                "@NX_class": u"NXcollection",
    #                'from' :  pk_pos_exact + bg_1[0],
    #                'to' : pk_pos_exact + bg_1[1],
    #                'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
    #            },
    #        'bg_2' : {
    #                "@NX_class": u"NXcollection",
    #                'from' :  pk_pos_exact + bg_2[0],
    #                'to' : pk_pos_exact + bg_2[1],
    #                'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
    #            }
    #    }
    #    
    #    return roi1Ddict
        
        
        
    #def set_
        
    def onChangeRockingScan(self, name):
        try:
            self.set_roscan(name)
        except Exception as e:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical, 'Invalid rocking scan info', 
            'Invalid or missig rocking scan info: %s.\n%s' % (name, e),
            qt.QMessageBox.Ok, self)
            msgbox.setDetailedText(traceback.format_exc())
            clickedbutton = msgbox.exec()
            return
        
        
    def get_rocking_scan_info(self, name):
        if name not in self.database.nxfile:
            raise ValueError("scan %s is not in the database" % name)
        
        h5_obj = self.database.nxfile[name]
        if 'orgui_meta' not in h5_obj.attrs or h5_obj.attrs['orgui_meta'] != 'rocking':
            raise ValueError("scan %s is not a valid rocking scan" % name)

        
        scangroup = h5_obj.parent.parent
        positioners = scangroup['instrument/positioners']
        if len(positioners) > 1:
            raise NotImplementedError("Multiple positioner changes are not yet implemented")
        
        axisname = list(positioners.keys())[0]
        axis = positioners[axisname][()]
        
        H_1 = h5_obj["rois"]['H_1'][()]
        H_0 = h5_obj["rois"]['H_0'][()]
        
        if not ( np.allclose(H_1.T[0], H_1.T[0][0]) and np.allclose(H_1.T[1], H_1.T[1][0]) and np.allclose(H_1.T[2], H_1.T[2][0]) ):
            raise ValueError("Rocking scan H_1 mismatch: Are these multiple scans?")
        
        if not ( np.allclose(H_0.T[0], H_0.T[0][0]) and np.allclose(H_0.T[1], H_0.T[1][0]) and np.allclose(H_0.T[2], H_0.T[2][0])):
            raise ValueError("Rocking scan H_0 mismatch: Are these multiple scans?")
        
        s_array = h5_obj["rois"]['s'][()]
        
        ddict = {
            'name' : name,
            'axisname' : axisname,
            'axis' : axis,
            's' : s_array,
            'H_0' : H_0[0],
            'H_1' : H_1[0]
        }
        return ddict
    
    
class IntegrationCorrectionsDialog(qt.QDialog):
    def __init__(self, parent=None):
        qt.QDialog.__init__(self, parent)
        verticalLayout = qt.QVBoxLayout(self)
        verticalLayout.setContentsMargins(0, 0, 0, 0)
        img = qutils.AspectRatioPixmapLabel(self)
        pixmp = qt.QPixmap(resources.getPath("incident_corrections.png"))
        img.setPixmap(pixmp)
        
        verticalLayout.addWidget(img)
        
        layout = qt.QGridLayout()
        
        self._L_save = 5
        self._beam_save = 20
        
        layout.addWidget(qt.QLabel("Sample size L:"),0, 0)

        self.L = qt.QDoubleSpinBox()
        self.L.setRange(0.00001, 1000000)
        self.L.setDecimals(4)
        self.L.setSuffix(u" mm")
        self.L.setValue(5)
        layout.addWidget(self.L, 0,1)
        
        layout.addWidget(qt.QLabel("beam size (FWHM):"),1, 0)
        self.beam = qt.QDoubleSpinBox()
        self.beam.setRange(0.000001, 1000000)
        self.beam.setDecimals(4)
        self.beam.setSuffix(u" microns")
        self.beam.setValue(20)
        layout.addWidget(self.beam, 1,1)
        

        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,2,0,1,-1)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.onCancel)
        
        verticalLayout.addLayout(layout)
        
        self.setLayout(verticalLayout)
        
    def onOk(self):
        self._L_save = self.L.value()
        self._beam_save = self.beam.value()
        self.accept()
        
    def onCancel(self):
        self.L.setValue(self._L_save)
        self.beam.setValue(self._beam_save)
        self.reject()
        


class IntegrationEstimator(qt.QDialog):

    
    def __init__(self,axis, axisname, parent=None):
        qt.QDialog.__init__(self, parent)
        
        layout = qt.QGridLayout()
        
        self._scan_label = qt.QLabel("%s-scan: from %s to %s\nroi ranges will be set relative to calculated peak position" % (axisname, axis[0], axis[-1]))
        layout.addWidget(self._scan_label,0,0, 1, -1)
        layout.addWidget(qt.QLabel("center roi:"),1, 0)
        layout.addWidget(qt.QLabel("from (rel):"),1, 1)
        layout.addWidget(qt.QLabel("to (rel):"),1, 3)
        
        layout.addWidget(qt.QLabel("bg roi 1:"),2,0)
        layout.addWidget(qt.QLabel("from (rel):"),2, 1)
        layout.addWidget(qt.QLabel("to (rel):"),2, 3)
        
        layout.addWidget(qt.QLabel("bg roi 2:"),3,0)
        layout.addWidget(qt.QLabel("from (rel):"),3, 1)
        layout.addWidget(qt.QLabel("to (rel):"),3, 3)

        self.croi_from = qt.QDoubleSpinBox()
        self.croi_from.setRange(-1000000, 1000000)
        self.croi_from.setDecimals(4)
        self.croi_from.setSuffix(u" °")
        self.croi_from.setValue(-0.1)
        layout.addWidget(self.croi_from, 1,2)
        
        self.croi_to = qt.QDoubleSpinBox()
        self.croi_to.setRange(-1000000, 1000000)
        self.croi_to.setDecimals(4)
        self.croi_to.setSuffix(u" °")
        self.croi_to.setValue(0.1)
        layout.addWidget(self.croi_to, 1,4)
        
        self.bgroi1_from = qt.QDoubleSpinBox()
        self.bgroi1_from.setRange(-1000000, 1000000)
        self.bgroi1_from.setDecimals(4)
        self.bgroi1_from.setSuffix(u" °")
        self.bgroi1_from.setValue(-0.3)
        layout.addWidget(self.bgroi1_from, 2,2)
        
        self.bgroi1_to = qt.QDoubleSpinBox()
        self.bgroi1_to.setRange(-1000000, 1000000)
        self.bgroi1_to.setDecimals(4)
        self.bgroi1_to.setSuffix(u" °")
        self.bgroi1_to.setValue(-0.12)
        layout.addWidget(self.bgroi1_to, 2,4)
        
        self.bgroi2_from = qt.QDoubleSpinBox()
        self.bgroi2_from.setRange(-1000000, 1000000)
        self.bgroi2_from.setDecimals(4)
        self.bgroi2_from.setSuffix(u" °")
        self.bgroi2_from.setValue(0.12)
        layout.addWidget(self.bgroi2_from, 3,2)
        
        self.bgroi2_to = qt.QDoubleSpinBox()
        self.bgroi2_to.setRange(-1000000, 1000000)
        self.bgroi2_to.setDecimals(4)
        self.bgroi2_to.setSuffix(u" °")
        self.bgroi2_to.setValue(0.3)
        layout.addWidget(self.bgroi2_to, 3,4)
        
        
        layout.addWidget(qt.QLabel("Attention: This will reset and override all existing ROIs!"),4,0, 1, -1)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,5,0, 1,-1)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def set_axis(self,axis, axisname):
        self._scan_label.setText("%s-scan: from %s to %s, roi ranges will be set relative to calculated peak position" % (axisname, axis[0], axis[-1]))
        
    def _verify_ranges(self):
        if self.croi_from.value() > self.croi_to.value():
            raise ValueError("Invalid input of center roi: from %s > to %s" % (self.croi_from.value(),self.croi_to.value()) )
        if self.bgroi1_from.value() > self.bgroi1_to.value():
            raise ValueError("Invalid input of bg roi 1: from %s > to %s" % (self.bgroi1_from.value(),self.bgroi1_to.value()) )
        if self.bgroi2_from.value() > self.bgroi2_to.value():
            raise ValueError("Invalid input of bg roi 2: from %s > to %s" % (self.bgroi2_from.value(),self.bgroi2_to.value()) )
        return True
        
    def get_ranges(self):
        self._verify_ranges()
        ddict = {
            'sig_1' : {
                'name' : 'sig_1',
                'from' : self.croi_from.value(),
                'to' : self.croi_to.value()
            },
            'bg_1' : {
                'name' : 'bg_1',
                'from' : self.bgroi1_from.value(),
                'to' : self.bgroi1_to.value()
            },
            'bg_2' : {
                'name' : 'bg_2',
                'from' : self.bgroi2_from.value(),
                'to' : self.bgroi2_to.value()
            }
        }
        return ddict
        
        
    def onOk(self):
        try:
            self._verify_ranges()
        except Exception as e:
            qt.QMessageBox.warning(self, "Invalid input", str(e))
            return
        self.accept()
        
        
class ROICreatorDialog(qt.QDialog):

    
    def __init__(self,axis, axisname, parent=None):
        qt.QDialog.__init__(self, parent)
        
        layout = qt.QGridLayout()
        
        self._scan_label = qt.QLabel("%s-scan: from %s to %s, roi ranges will be set relative to calculated peak position" % (axisname, axis[0], axis[-1]))
        layout.addWidget(self._scan_label,0,0, 1, -1)
        
        
        
        layout.addWidget(qt.QLabel("roi:"),1, 0)
        layout.addWidget(qt.QLabel("from (rel):"),1, 1)
        layout.addWidget(qt.QLabel("to (rel):"),1, 3)

        self.roi_from = qt.QDoubleSpinBox()
        self.roi_from.setRange(-1000000, 1000000)
        self.roi_from.setDecimals(4)
        self.roi_from.setSuffix(u" °")
        self.roi_from.setValue(-0.1)
        layout.addWidget(self.roi_from, 1,2)
        
        self.roi_to = qt.QDoubleSpinBox()
        self.roi_to.setRange(-1000000, 1000000)
        self.roi_to.setDecimals(4)
        self.roi_to.setSuffix(u" °")
        self.roi_to.setValue(0.1)
        layout.addWidget(self.roi_to, 1,4)
        
        checkboxlayout = qt.QHBoxLayout()
        
        
        self.signalCheckbox = qt.QCheckBox("Signal")
        self.backgroundCheckbox = qt.QCheckBox("Background")
        
        self.btngroup = qt.QButtonGroup()
        self.btngroup.addButton(self.signalCheckbox)
        self.btngroup.addButton(self.backgroundCheckbox)
        self.btngroup.setExclusive(True)
        
        self.signalCheckbox.setChecked(True)
        
        checkboxlayout.addWidget(qt.QLabel("ROI type:"))
        checkboxlayout.addWidget(self.signalCheckbox)
        checkboxlayout.addWidget(self.backgroundCheckbox)
        
        layout.addLayout(checkboxlayout, 2, 0, 1, -1)
        
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,3,0,1,-1)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def set_axis(self,axis, axisname):
        self._scan_label.setText("%s-scan: from %s to %s, roi ranges will be set relative to calculated peak position" % (axisname, axis[0], axis[-1]))
        
    def _verify_ranges(self):
        if self.roi_from.value() > self.roi_to.value():
            raise ValueError("Invalid input of roi: from %s > to %s" % (self.roi_from.value(),self.roi_to.value()) )
        return True
        
    def get_roi(self):
        self._verify_ranges()
        signal = self.signalCheckbox.isChecked()
        
        ddict = {
            'from' : self.roi_from.value(),
            'to' : self.roi_to.value(),
            'signal' : signal
        }
        return ddict
        
        
    def onOk(self):
        try:
            self._verify_ranges()
        except Exception as e:
            qt.QMessageBox.warning(self, "Invalid input", str(e))
            return
        self.accept()
        
        
class AnchorROI(ROI):
    
    #sigAnchorChanged = qt.Signal()
    
    def __init__(self, name, fromdata=None, todata=None, type_=None, anchor=False):
        ROI.__init__(self, name, fromdata, todata, type_)
        self._anchor = bool(anchor)
        self.anchor_updating = False
    
    def isAnchor(self):
        return self._anchor
        
    def setAnchor(self, anchor):
        #print(self.getName(), self._anchor, anchor)
        if self._anchor != bool(anchor):
            self._anchor = bool(anchor)
            self.sigChanged.emit()
        
        
    def toDict(self):
        """

        :return: dict containing the roi parameters
        """
        ddict = super().toDict()
        ddict['anchor'] = self.isAnchor()
        return ddict

    @staticmethod
    def _fromDict(dic):
        assert "name" in dic
        roi = AnchorROI(name=dic["name"])
        roi._extraInfo = {}
        for key in dic:
            if key == "from":
                roi.setFrom(dic["from"])
            elif key == "to":
                roi.setTo(dic["to"])
            elif key == "type":
                roi.setType(dic["type"])
            elif key == "anchor":
                roi.setAnchor(dic["anchor"])
            else:
                roi._extraInfo[key] = dic[key]
        return roi
        
        

        
        
        
class SelectableROITable(ROITable):
    COLUMNS_INDEX = dict(
        [
            ("Anchor", 0),
            ("ID", 1),
            ("ROI", 2),
            ("Type", 3),
            ("From", 4),
            ("To", 5),
            ("Raw Counts", 6),
            ("Net Counts", 7),
            ("Raw Area", 8),
            ("Net Area", 9)
        ]
    )
    
    COLUMNS = list(COLUMNS_INDEX.keys())
    
    
    
    def __init__(self, parent=None, plot=None, rois=None):
        super().__init__(parent, plot, rois)
        self.isModelSetting = False
    
    def _updateRoiInfo(self, roiID):
        if self._userIsEditingRoi is True:
            return
        if roiID not in self._roiDict:
            return
        super()._updateRoiInfo(roiID)
        
        roi = self._roiDict[roiID]
        itemID = self._getItem(name="ID", roi=roi, row=None)
        
        itemAnchor = self._getItem(name="Anchor", row=itemID.row(), roi=roi)
        #print('update roi info: ', roi.getName(), roi.isAnchor())
        self.setAnchorState(roiID, roi.isAnchor())
            
    def setAnchorState(self, roiID , anchor):
        if roiID not in self._roiDict:
            return
        roi = self._roiDict[roiID]
        itemID = self._getItem(name="ID", roi=roi, row=None)
        
            
        itemAnchor = self._getItem(name="Anchor", row=itemID.row(), roi=roi)
        #with qt.QSignalBlocker(self):
        if anchor:
            itemAnchor.setCheckState(qt.Qt.Checked)
        else:
            itemAnchor.setCheckState(qt.Qt.Unchecked)
    
    def addRoi(self, roi):
        """

        :param :class:`ROI` roi: roi to add to the table
        """
        assert isinstance(roi, AnchorROI)
        self._getItem(name="ID", row=None, roi=roi)
        self._roiDict[roi.getID()] = roi
        self._markersHandler.add(roi, _ColorRoiMarkerHandler(roi, self.plot))
        self._updateRoiInfo(roi.getID())
        callback = partial(WeakMethodProxy(self._updateRoiInfo), roi.getID())
        roi.sigChanged.connect(callback)
        # set it as the active one
        self.setActiveRoi(roi)
        
        
    def setRois(self, rois, order=None):
        """Set the ROIs by providing a dictionary of ROI information.

        The dictionary keys are the ROI names.
        Each value is a sub-dictionary of ROI info with the following fields:

        - ``"from"``: x coordinate of the left limit, as a float
        - ``"to"``: x coordinate of the right limit, as a float
        - ``"type"``: type of ROI, as a string (e.g "channels", "energy")


        :param roidict: Dictionary of ROIs
        :param str order: Field used for ordering the ROIs.
             One of "from", "to", "type".
             None (default) for no ordering, or same order as specified
             in parameter ``rois`` if provided as a dict.
        """
        assert order in [None, "from", "to", "type"]
        self.clear()

        # backward compatibility since 0.10.0
        if isinstance(rois, dict):
            for roiName, roi in rois.items():
                if isinstance(roi, AnchorROI):
                    _roi = roi
                else:
                    roi["name"] = roiName
                    _roi = AnchorROI._fromDict(roi)
                self.addRoi(_roi)
        else:
            for roi in rois:
                assert isinstance(roi, AnchorROI)
                self.addRoi(roi)
        self._updateMarkers()
    
    def load(self, filename):
        """
        Load ROI widget information from a file storing a dict of ROI.

        :param str filename: The file from which to load ROI
        """
        roisDict = dictdump.load(filename)
        rois = []

        # Remove rawcounts and netcounts from ROIs
        for roiDict in roisDict["ROI"]["roidict"].values():
            roiDict.pop("rawcounts", None)
            roiDict.pop("netcounts", None)
            rois.append(AnchorROI._fromDict(roiDict))

        self.setRois(rois)
        
        
    def _getItem(self, name, row, roi):
        if row:
            item = self.item(row, self.COLUMNS_INDEX[name])
        else:
            item = None
        if item:
            return item
        else:
            if name == "ID":
                assert roi
                if roi.getID() in self._roiToItems:
                    return self._roiToItems[roi.getID()]
                else:
                    # create a new row
                    row = self.rowCount()
                    self.setRowCount(self.rowCount() + 1)
                    item = qt.QTableWidgetItem(
                        str(roi.getID()), type=qt.QTableWidgetItem.Type
                    )
                    self._roiToItems[roi.getID()] = item
            elif name == "ROI":
                item = qt.QTableWidgetItem(
                    roi.getName() if roi else "", type=qt.QTableWidgetItem.Type
                )
                if roi.getName().upper() in ("ICR", "DEFAULT"):
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                else:
                    item.setFlags(
                        qt.Qt.ItemIsSelectable
                        | qt.Qt.ItemIsEnabled
                        | qt.Qt.ItemIsEditable
                    )
            elif name == "Type":
                item = qt.QTableWidgetItem(type=qt.QTableWidgetItem.Type)
                item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
            elif name in ("To", "From"):
                item = _FloatItem()
                if roi.getName().upper() in ("ICR", "DEFAULT"):
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                else:
                    item.setFlags(
                        qt.Qt.ItemIsSelectable
                        | qt.Qt.ItemIsEnabled
                        | qt.Qt.ItemIsEditable
                    )
            elif name in ("Raw Counts", "Net Counts", "Raw Area", "Net Area"):
                item = _FloatItem()
                item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
            elif name == "Anchor":
                item = qt.QTableWidgetItem(type=qt.QTableWidgetItem.Type)
                if roi.getName().upper() in ("ICR", "DEFAULT"):
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                else:
                    item.setFlags(
                        qt.Qt.ItemIsSelectable
                        | qt.Qt.ItemIsEnabled
                        | qt.Qt.ItemIsUserCheckable
                    )
                    item.setIcon(resources.getQicon("anchor-ROI"))
                    item.setCheckState(qt.Qt.Unchecked)
            else:
                raise ValueError("item type not recognized")

            self.setItem(row, self.COLUMNS_INDEX[name], item)
            return item
            
    def _itemChanged(self, item):
        def getRoi():
            IDItem = self.item(item.row(), self.COLUMNS_INDEX["ID"])
            assert IDItem
            id = int(IDItem.text())
            assert id in self._roiDict
            roi = self._roiDict[id]
            return roi

        def signalChanged(roi):
            if self.activeRoi and roi.getID() == self.activeRoi.getID():
                self.activeROIChanged.emit()
        
        super()._itemChanged(item)

        
        self._userIsEditingRoi = True
        if not self.isModelSetting:
            if item.column() == self.COLUMNS_INDEX["Anchor"]:
                roi = getRoi()
                anchor = item.checkState() == qt.Qt.Checked
                if anchor != roi.isAnchor():
                    #print('anchor changed:', getRoi().getName() ,getRoi().getID())
                    roi.setAnchor(anchor)

        self._userIsEditingRoi = False
    
        
class CurvesROIWidget(qt.QWidget):
    """
    Widget displaying a table of ROI information.

    Implements also the following behavior:

    * if the roiTable has no ROI when showing create the default ICR one

    :param parent: See :class:`QWidget`
    :param str name: The title of this widget
    """

    sigROIWidgetSignal = qt.Signal(object)
    """Signal of ROIs modifications.

    Modification information if given as a dict with an 'event' key
    providing the type of events.

    Type of events:

    - AddROI, DelROI, LoadROI and ResetROI with keys: 'roilist', 'roidict'
    - selectionChanged with keys: 'row', 'col' 'roi', 'key', 'colheader',
      'rowheader'
    """

    sigROISignal = qt.Signal(object)

    def __init__(self, parent=None, name=None, plot=None):
        super(CurvesROIWidget, self).__init__(parent)
        if name is not None:
            self.setWindowTitle(name)
        self.__lastSigROISignal = None
        """Store the last value emitted for the sigRoiSignal. In the case the
        active curve change we need to add this extra step in order to make
        sure we won't send twice the sigROISignal.
        This come from the fact sigROISignal is connected to the 
        activeROIChanged signal which is emitted when raw and net counts
        values are changing but are not embed in the sigROISignal.
        """
        assert plot is not None
        self._plotRef = weakref.ref(plot)
        self._showAllMarkers = False
        self.currentROI = None

        layout = qt.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.roiTable = SelectableROITable(self, plot=plot)
        rheight = self.roiTable.horizontalHeader().sizeHint().height()
        self.roiTable.setMinimumHeight(4 * rheight)
        layout.addWidget(self.roiTable)
        self._roiFileDir = qt.QDir.home().absolutePath()
        #self._showAllCheckBox.toggled.connect(self.roiTable.showAllMarkers)
        self.roiTable.showAllMarkers(True)


        self._isConnected = False  # True if connected to plot signals
        self._isInit = False

    def getROIListAndDict(self):
        return self.roiTable.getROIListAndDict()

    def getPlotWidget(self):
        """Returns the associated PlotWidget or None

        :rtype: Union[~silx.gui.plot.PlotWidget,None]
        """
        return None if self._plotRef is None else self._plotRef()

    def showEvent(self, event):
        self._visibilityChangedHandler(visible=True)
        qt.QWidget.showEvent(self, event)

    @property
    def roiFileDir(self):
        """The directory from which to load/save ROI from/to files."""
        if not os.path.isdir(self._roiFileDir):
            self._roiFileDir = qt.QDir.home().absolutePath()
        return self._roiFileDir

    @roiFileDir.setter
    def roiFileDir(self, roiFileDir):
        self._roiFileDir = str(roiFileDir)

    def setRois(self, rois, order=None):
        return self.roiTable.setRois(rois, order)

    def getRois(self, order=None):
        return self.roiTable.getRois(order)

    def setMiddleROIMarkerFlag(self, flag=True):
        return self.roiTable.setMiddleROIMarkerFlag(flag)

    def _add(self):
        """Add button clicked handler"""

        def getNextRoiName():
            rois = self.roiTable.getRois(order=None)
            roisNames = []
            [roisNames.append(roiName) for roiName in rois]
            nrois = len(rois)
            if nrois == 0:
                return "ICR"
            else:
                i = 1
                newroi = "newroi %d" % i
                while newroi in roisNames:
                    i += 1
                    newroi = "newroi %d" % i
                return newroi

        roi = AnchorROI(name=getNextRoiName())

        if roi.getName() == "ICR":
            roi.setType("Default")
        else:
            roi.setType(self.getPlotWidget().getXAxis().getLabel())

        xmin, xmax = self.getPlotWidget().getXAxis().getLimits()
        fromdata = xmin + 0.25 * (xmax - xmin)
        todata = xmin + 0.75 * (xmax - xmin)
        if roi.isICR():
            fromdata, dummy0, todata, dummy1 = self._getAllLimits()
        roi.setFrom(fromdata)
        roi.setTo(todata)
        self.roiTable.addRoi(roi)

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "AddROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def _del(self):
        """Delete button clicked handler"""
        self.roiTable.deleteActiveRoi()

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "DelROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def _reset(self):
        """Reset button clicked handler"""
        self.roiTable.clear()
        old = self.blockSignals(True)  # avoid several sigROISignal emission
        self._add()
        self.blockSignals(old)

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "ResetROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def _load(self):
        """Load button clicked handler"""
        dialog = qt.QFileDialog(self)
        dialog.setNameFilters(["INI File  *.ini", "JSON File *.json", "All *.*"])
        dialog.setFileMode(qt.QFileDialog.ExistingFile)
        dialog.setDirectory(self.roiFileDir)
        if not dialog.exec():
            dialog.close()
            return

        # pyflakes bug http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=666494
        outputFile = dialog.selectedFiles()[0]
        dialog.close()

        self.roiFileDir = os.path.dirname(outputFile)
        self.roiTable.load(outputFile)

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "LoadROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def load(self, filename):
        """Load ROI widget information from a file storing a dict of ROI.

        :param str filename: The file from which to load ROI
        """
        self.roiTable.load(filename)

    def _save(self):
        """Save button clicked handler"""
        dialog = qt.QFileDialog(self)
        dialog.setNameFilters(["INI File  *.ini", "JSON File *.json"])
        dialog.setFileMode(qt.QFileDialog.AnyFile)
        dialog.setAcceptMode(qt.QFileDialog.AcceptSave)
        dialog.setDirectory(self.roiFileDir)
        if not dialog.exec():
            dialog.close()
            return

        outputFile = dialog.selectedFiles()[0]
        extension = "." + dialog.selectedNameFilter().split(".")[-1]
        dialog.close()

        if not outputFile.endswith(extension):
            outputFile += extension

        if os.path.exists(outputFile):
            try:
                os.remove(outputFile)
            except IOError:
                msg = qt.QMessageBox(self)
                msg.setIcon(qt.QMessageBox.Critical)
                msg.setText("Input Output Error: %s" % (sys.exc_info()[1]))
                msg.exec()
                return
        self.roiFileDir = os.path.dirname(outputFile)
        self.save(outputFile)

    def save(self, filename):
        """Save current ROIs of the widget as a dict of ROI to a file.

        :param str filename: The file to which to save the ROIs
        """
        self.roiTable.save(filename)

    def calculateRois(self, roiList=None, roiDict=None):
        """Compute ROI information"""
        return self.roiTable.calculateRois()

    def showAllMarkers(self, _show=True):
        self.roiTable.showAllMarkers(_show)

    def _getAllLimits(self):
        """Retrieve the limits based on the curves."""
        plot = self.getPlotWidget()
        curves = () if plot is None else plot.getAllCurves()
        if not curves:
            return 1.0, 1.0, 100.0, 100.0

        xmin, ymin = None, None
        xmax, ymax = None, None

        for curve in curves:
            x = curve.getXData(copy=False)
            y = curve.getYData(copy=False)
            if xmin is None:
                xmin = x.min()
            else:
                xmin = min(xmin, x.min())
            if xmax is None:
                xmax = x.max()
            else:
                xmax = max(xmax, x.max())
            if ymin is None:
                ymin = y.min()
            else:
                ymin = min(ymin, y.min())
            if ymax is None:
                ymax = y.max()
            else:
                ymax = max(ymax, y.max())

        return xmin, ymin, xmax, ymax

    def showEvent(self, event):
        self._visibilityChangedHandler(visible=True)
        qt.QWidget.showEvent(self, event)

    def hideEvent(self, event):
        self._visibilityChangedHandler(visible=False)
        qt.QWidget.hideEvent(self, event)

    def _visibilityChangedHandler(self, visible):
        """Handle widget's visibility updates.

        It is connected to plot signals only when visible.
        """
        if visible:
            # if no ROI existing yet, add the default one
            if self.roiTable.rowCount() == 0:
                old = self.blockSignals(True)  # avoid several sigROISignal emission
                self._add()
                self.blockSignals(old)
                self.calculateRois()

    def fillFromROIDict(self, *args, **kwargs):
        self.roiTable.fillFromROIDict(*args, **kwargs)

    def _emitCurrentROISignal(self):
        ddict = {}
        ddict["event"] = "currentROISignal"
        if self.roiTable.activeRoi is not None:
            ddict["ROI"] = self.roiTable.activeRoi.toDict()
            ddict["current"] = self.roiTable.activeRoi.getName()
        else:
            ddict["current"] = None

        if self.__lastSigROISignal != ddict:
            self.__lastSigROISignal = ddict
            self.sigROISignal.emit(ddict)

            
    @property
    def currentRoi(self):
        return self.roiTable.activeRoi

class _ColorRoiMarkerHandler(_RoiMarkerHandler):
    
    def __init__(self, roi, plot):
        super().__init__(roi, plot)
        if roi.getName().startswith('sig'):
            self._color = 'red'
        elif roi.getName().startswith('bg'):
            self._color = 'blue'
        else:
            self._color = 'black'
        
    



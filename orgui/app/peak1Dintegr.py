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
from silx.gui.dialog import ImageFileDialog
from silx.gui.plot.actions import control as control_actions
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
#from silx.gui.widgets.TableWidget import TableWidget
from silx.gui.plot.CurvesROIWidget import ROITable, ROI

import traceback

from . import qutils
from .database import DataBase
from .. import resources

import numpy as np
from ..datautils.xrayutils import HKLVlieg, CTRcalc
from ..datautils.xrayutils import ReciprocalNavigation as rn

from silx.io.dictdump import dicttonx ,nxtodict

import sys

QTVERSION = qt.qVersion()
DEBUG = 0

MAX_ROIS_DISPLAY = 100


class RockingPeakIntegrator(qt.QMainWindow):
    def __init__(self, database, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.database = database
        
        self._currentRoInfo = {}
        self._idx = 0
        
        dbdockwidget = qt.QDockWidget("Integrated data")
        
        self.dbside = qt.QMainWindow()
        
        
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
        self.layout.addWidget(self.roiwidget, 1)
        

        
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
        
    def set_roscan(self, name):
        ro_info = self.get_rocking_scan_info(name)
        if ro_info['name'] + '/integration' not in self.database.nxfile:
            roi1D_info = self.estimate_roi1D_info(ro_info)
            self.database.add_nxdict(roi1D_info, update_mode='modify', h5path=ro_info['name'] + '/integration')
        
        self._currentRoInfo = ro_info
        self._idx = 0

        self.curveSlider.setAxis(self._currentRoInfo['s'], "s")
        self.curveSlider.setIndex(0)
        self.plotRoCurve(0)
        self.plotROIselect.resetZoom()
        if self.autozoom_checkbox.isChecked():
            self.resetXZoomScaled(self.zoomslider.value())


    #ddict = {
    #        "event": "indexChanged",
    #        "oldtxt": self._lineTxt,
    #        "newtxt": txt,
    #        "idx" : idx,
    #        "value" : self._data[idx],
    #        "id": id(self),
    #    }
    def onSliderValueChanged(self, ddict):
        self.plotRoCurve(ddict['idx'])
        
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
        
    def plotRoCurve(self, idx):
        ro_curve = self.get_ro_curve(idx)
        s = self._currentRoInfo['s'][idx]
        hkl = self._currentRoInfo['H_1'] * s + self._currentRoInfo['H_0']
        title = "Rocking scan at s = %s, HKL = [%.2f %.2f %.2f]" % (s, *hkl)
        
        self.plotROIselect.clear()
        
        lbl = self.plotROIselect.addCurve(ro_curve['axis'], ro_curve['croibg'], legend=title, 
                           xlabel="%s / deg" % self._currentRoInfo['axisname'], 
                           ylabel='center ROI background subtracted / cnts', 
                           yerror=ro_curve['croibg_errors'], resetzoom=False)
        
        self.plotROIselect.setActiveCurve(lbl)
        self._idx = idx
                           
        self.roiwidget.roiTable.clear()
        
        roi_dict = self.get_roi1D_info(idx)
        
        minfrom = np.inf
        maxto = -np.inf
        for key in roi_dict:
            if key.startswith('sig'):
                roi_d = roi_dict[key]
                roi = ROI(key, fromdata=roi_d['from'], todata=roi_d['to'], type_=str(roi_d['type']))
                minfrom = min(minfrom, roi_d['from'])
                maxto = max(maxto, roi_d['to'])
                self.roiwidget.roiTable.addRoi(roi)
            elif key.startswith('bg'): # color?
                roi_d = roi_dict[key]
                roi = ROI(key, fromdata=roi_d['from'], todata=roi_d['to'], type_=str(roi_d['type']))
                minfrom = min(minfrom, roi_d['from'])
                maxto = max(maxto, roi_d['to'])
                self.roiwidget.roiTable.addRoi(roi)
        
        # toDo: set from GUI
        if self.autozoom_checkbox.isChecked():
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
            print("Cannot zoom:" , e)
            # can fail, add better error handling later!
        
                
                
    
    def get_roi1D_info(self, idx):
        name = self._currentRoInfo['name']
        h5_obj = self.database.nxfile[name]['integration']
        
        roi1Ddict = {}
        for k in h5_obj:
            if k.startswith('sig') or k.startswith('bg') :
                roi1Ddict[k] = {
                    'name' : k,
                    'from' : h5_obj[k]['from'][idx][()],
                    'to' : h5_obj[k]['to'][idx][()],
                    'type' : 'signal',
                    'fixed' : h5_obj[k]['fixed'][idx][()]
                }
        return roi1Ddict

    def estimate_roi1D_info(self, ro_info):
        # ToDo make settings available from GUI 
        
        # estimate peak pos
        axis_pk_pos = [] 
        axis_pk_pos_idx = [] 

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
        #idx_pk = np.argmin(np.abs(sc_h5['axis'][()] - pk_pos_exact), axis=1)

        # ToDo make settings available from GUI 
        sig_width = 0.1 # deg
        bg_1 = -0.2 , -0.1
        bg_2 = 0.1 , 0.2
        
        roi1Ddict = {
            "@NX_class": u"NXcollection",
            "@info" : u"ROI information for the rocking scan integration",
            "peakpos" : pk_pos_exact,
            'sig_1' : {
                    "@NX_class": u"NXcollection",
                    'from' :  pk_pos_exact - 0.5*sig_width,
                    'to' : pk_pos_exact + 0.5*sig_width,
                    'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
                },
            'bg_1' : {
                    "@NX_class": u"NXcollection",
                    'from' :  pk_pos_exact + bg_1[0],
                    'to' : pk_pos_exact + bg_1[1],
                    'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
                },
            'bg_2' : {
                    "@NX_class": u"NXcollection",
                    'from' :  pk_pos_exact + bg_2[0],
                    'to' : pk_pos_exact + bg_2[1],
                    'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
                }
        }
        
        return roi1Ddict
        
        
        
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
    
    
    def __estimate_roi1D_info_old_manygroups(self, ro_info):
        # ToDo make settings available from GUI 
        
        # estimate peak pos
        axis_pk_pos = [] 
        axis_pk_pos_idx = [] 
        for sc in ro_info['scanpath']:
            sc_h5 = self.database.nxfile[sc]
            sixc_angles_pk = sc_h5['trajectory/HKL_sixc_angles']
            if ro_info['axisname'] == 'mu':
                pk_pos_exact = sixc_angles_pk['alpha'][()]
            elif ro_info['axisname'] == 'th':
                pk_pos_exact = sixc_angles_pk['theta'][()]
            else:
                raise ValueError("Cannot estimate peak position: unknown scan axis %s" % ro_info['axisname'])
            idx_pk = np.argmin(np.abs(ro_info['axis'] - pk_pos_exact))
            
            axis_pk_pos.append(pk_pos_exact)
            axis_pk_pos_idx.append(idx_pk)

        axis_pk_pos = np.array(axis_pk_pos)
        axis_pk_pos_idx = np.array(axis_pk_pos_idx)
        
        # ToDo make settings available from GUI 
        sig_width = 0.1 # deg
        bg_1 = -0.2 , -0.1
        bg_2 = 0.1 , 0.2
        
        roi1Ddict = {}
        for i, sc in enumerate(ro_info['scans']):
            roi1Ddict[sc] = {
                "@NX_class": u"NXcollection",
                "@info" : u"ROI information for the rocking scan integration",
                "peak_position" : axis_pk_pos[i]
            }
            roi1Ddict[sc]['sigROI_1'] = {
                'name' : 'sigROI_1',
                'from' : axis_pk_pos[i] - 0.5*sig_width,
                'to' : axis_pk_pos[i] + 0.5*sig_width,
                'type' : 'signal',
                "@NX_class": u"NXcollection"
            }
            roi1Ddict[sc]['bgROI_1'] = {
                'name' : 'bgROI_1',
                'from' : axis_pk_pos[i] + bg_1[0],
                'to' : axis_pk_pos[i] + bg_1[1],
                'type' : 'signal',
                "@NX_class": u"NXcollection"
            }
            roi1Ddict[sc]['bgROI_2'] = {
                'name' : 'bgROI_2',
                'from' : axis_pk_pos[i] + bg_2[0],
                'to' : axis_pk_pos[i] + bg_2[1],
                'type' : 'signal',
                "@NX_class": u"NXcollection"
            }
        
        return roi1Ddict
        
        
    def __get_rocking_scan_info_old_manygroups(self, name):
        if name not in self.database.nxfile:
            raise ValueError("scan %s is not in the database" % name)
        
        h5_obj = self.database.nxfile[name]
        if 'orgui_meta' not in h5_obj.attrs or h5_obj.attrs['orgui_meta'] != 'rocking':
            raise ValueError("scan %s is not a valid rocking scan" % name)
        
        scans = np.array(list(h5_obj["rois"].keys()))
        
        scangroup = h5_obj.parent.parent
        positioners = scangroup['instrument/positioners']
        if len(positioners) > 1:
            raise NotImplementedError("Multiple positioner changes are not yet implemented")
        
        axisname = list(positioners.keys())[0]
        axis = positioners[axisname][()]
        
        H_1 = h5_obj["rois"][scans[0]]['trajectory/H_1'][()]
        H_0 = h5_obj["rois"][scans[0]]['trajectory/H_0'][()]
        
        s_array = []
        for sc in scans:
            s = h5_obj["rois"][sc]['trajectory/s'][()]
            s_array.append(s)
            if not np.all(np.isclose(h5_obj["rois"][sc]['trajectory/H_1'][()], H_1)):
                raise ValueError("Rocking scan H_1 mismatch: is: %s, expected: %s. Are these multiple scans?" % 
                    (h5_obj["rois"][sc]['trajectory/H_1'][()], H_1) )
            if not np.all(np.isclose(h5_obj["rois"][sc]['trajectory/H_0'][()], H_0)):
                raise ValueError("Rocking scan H_0 mismatch: is: %s, expected: %s. Are these multiple scans?" % 
                    (h5_obj["rois"][sc]['trajectory/H_0'][()], H_0) )
        
        s_array = np.array(s_array)
        order = np.argsort(s_array)
        
        scans = scans[order]
        s_array = s_array[order]
        
        fullscanspath = []
        for sc in scans:
            fullscanspath.append(name + "/rois/" + sc)
        
        ddict = {
            'name' : name,
            'axisname' : axisname,
            'axis' : axis,
            'scans' : scans,
            'scanpath' : fullscanspath,
            's' : s_array,
            'H_0' : H_0,
            'H_1' : H_1
        }
        return ddict
        
        
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

        #widgetAllCheckbox = qt.QWidget(parent=self)
        #self._showAllCheckBox = qt.QCheckBox("show all ROI", parent=widgetAllCheckbox)
        #widgetAllCheckbox.setLayout(qt.QHBoxLayout())
        #spacer = qt.QWidget(parent=widgetAllCheckbox)
        #spacer.setSizePolicy(qt.QSizePolicy.Expanding, qt.QSizePolicy.Fixed)
        #widgetAllCheckbox.layout().addWidget(spacer)
        #widgetAllCheckbox.layout().addWidget(self._showAllCheckBox)
        #layout.addWidget(widgetAllCheckbox)

        self.roiTable = ROITable(self, plot=plot)
        rheight = self.roiTable.horizontalHeader().sizeHint().height()
        self.roiTable.setMinimumHeight(4 * rheight)
        layout.addWidget(self.roiTable)
        self._roiFileDir = qt.QDir.home().absolutePath()
        #self._showAllCheckBox.toggled.connect(self.roiTable.showAllMarkers)
        self.roiTable.showAllMarkers(True)

        hbox = qt.QWidget(self)
        hboxlayout = qt.QHBoxLayout(hbox)
        hboxlayout.setContentsMargins(0, 0, 0, 0)
        hboxlayout.setSpacing(0)

        hboxlayout.addStretch(0)

        self.addButton = qt.QPushButton(hbox)
        self.addButton.setText("Add ROI")
        self.addButton.setToolTip("Create a new ROI")
        self.delButton = qt.QPushButton(hbox)
        self.delButton.setText("Delete ROI")
        self.addButton.setToolTip("Remove the selected ROI")
        self.resetButton = qt.QPushButton(hbox)
        self.resetButton.setText("Reset")
        self.addButton.setToolTip(
            "Clear all created ROIs. We only let the " "default ROI"
        )

        hboxlayout.addWidget(self.addButton)
        hboxlayout.addWidget(self.delButton)
        hboxlayout.addWidget(self.resetButton)

        hboxlayout.addStretch(0)

        self.loadButton = qt.QPushButton(hbox)
        self.loadButton.setText("Load")
        self.loadButton.setToolTip("Load ROIs from a .ini file")
        self.saveButton = qt.QPushButton(hbox)
        self.saveButton.setText("Save")
        self.loadButton.setToolTip("Save ROIs to a .ini file")
        hboxlayout.addWidget(self.loadButton)
        hboxlayout.addWidget(self.saveButton)
        layout.setStretchFactor(self.roiTable, 1)
        layout.setStretchFactor(hbox, 0)

        layout.addWidget(hbox)

        # Signal / Slot connections
        self.addButton.clicked.connect(self._add)
        self.delButton.clicked.connect(self._del)
        self.resetButton.clicked.connect(self._reset)

        self.loadButton.clicked.connect(self._load)
        self.saveButton.clicked.connect(self._save)

        self.roiTable.activeROIChanged.connect(self._emitCurrentROISignal)

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

        roi = ROI(name=getNextRoiName())

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





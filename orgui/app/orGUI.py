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
__credits__ = ['Finn Schroeter']
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
from .. import __version__
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import gc
import sys
import os
from silx.gui import qt
import warnings

from io import StringIO
import concurrent.futures
import queue
import threading

#from IPython import embed
import silx.gui.plot
from silx.gui.plot import items
from silx.gui.colors import Colormap
import weakref

#from silx import sx

import silx
from silx.utils.weakref import WeakMethodProxy
from silx.gui.plot.Profile import ProfileToolBar
from silx.gui.plot.AlphaSlider import NamedImageAlphaSlider
from silx.gui.dialog import ImageFileDialog
from silx.gui.plot.tools.roi import RegionOfInterestManager
from silx.gui.plot.tools.roi import RegionOfInterestTableWidget
from silx.gui.plot.items.roi import RectangleROI, PolygonROI, ArcROI
from silx.gui.plot.actions import control as control_actions

try:
    from silx.gui import console
except:
    console = False

import traceback

from . import qutils
from .QScanSelector import QScanSelector
from .QReflectionSelector import QReflectionSelector, QReflectionAnglesDialog
from .QUBCalculator import QUBCalculator
from .peak1Dintegr import RockingPeakIntegrator
from .ArrayTableDialog import ArrayTableDialog
from .bgroi import RectangleBgROI
from .database import DataBase
from ..backend.scans import SimulationScan
from ..backend import backends
from ..backend import universalScanLoader
from .. import resources

try:
    from . import _roi_sum_accel
    HAS_ACCEL = True
except:
    print(traceback.format_exc())
    HAS_ACCEL = False

import numpy as np
from ..datautils.xrayutils import HKLVlieg, CTRcalc
from ..datautils.xrayutils import ReciprocalNavigation as rn
import pyFAI.detectors

import sys
#legacy import:
from ..backend.beamline.id31_tools import BlissScan_EBS, Fastscan, BlissScan

QTVERSION = qt.qVersion()
DEBUG = 0

MAX_ROIS_DISPLAY = 100

silx.config.DEFAULT_PLOT_SYMBOL = '.'


class orGUI(qt.QMainWindow):
    def __init__(self,configfile,parent=None):
        qt.QMainWindow.__init__(self, parent)
        #self.setWindowIcon(resources.getQicon("orguiicon"))
        self.h5database = None # must be a h5py file-like, by default not opened to avoid reading issues at beamtimes!
        self.images_loaded = False
        self.resetZoom = True
        #icon = resources.getQicon("sum_image.svg")
        self.fscan = None
        self.activescanname = "scan"
        self.numberthreads = int(min(os.cpu_count(), 16)) if os.cpu_count() is not None else 1 
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            self.numberthreads = int(os.environ['SLURM_CPUS_ON_NODE'])
        
        self.filedialogdir = os.getcwd()
        
        self.excludedImagesDialog = ArrayTableDialog(True, 1)
        self.excludedImagesDialog.setArrayData(np.array([-1]),editable=True, header= ['image no'])
        
        
        self.centralPlot = Plot2DHKL(self.newXyHKLConverter(),parent=self)
        self.centralPlot.setDefaultColormap(Colormap(name='jet',normalization='log'))
        self.centralPlot.setCallback(self._graphCallback)
        toolbar = qt.QToolBar()
        toolbar.addAction(control_actions.OpenGLAction(parent=toolbar, plot=self.centralPlot))
        self.centralPlot.addToolBar(toolbar)
        
        self.currentImageLabel = None
        self.currentAddImageLabel = None
        
        selectorDock = qt.QDockWidget("Scan data")
        selectorDock.setAllowedAreas(qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea)
        self.scanSelector = QScanSelector(self)
        selectorDock.setWidget(self.scanSelector)
        self.addDockWidget(qt.Qt.LeftDockWidgetArea,selectorDock)
    
        self.imagepath = ''
        self.imageno = 0
        
    
        ubWidget = qt.QSplitter(qt.Qt.Vertical)
        ubWidget.setChildrenCollapsible(False)
        #ubLayout = qt.QVBoxLayout()
        self.ubcalc = QUBCalculator(configfile, self)
        self.ubcalc.sigNewReflection.connect(self._onNewReflection)
        
        
        
        maincentralwidget = qt.QTabWidget()
        
        self.integrdataPlot = silx.gui.plot.Plot1D(self)
        legendwidget = self.integrdataPlot.getLegendsDockWidget()
        
        toolbar = qt.QToolBar()
        toolbar.addAction(control_actions.OpenGLAction(parent=toolbar, plot=self.integrdataPlot))
        self.integrdataPlot.addToolBar(toolbar)
        
        self.integrdataPlot.addDockWidget(qt.Qt.RightDockWidgetArea,legendwidget)
        legendwidget.show()
        self.database = DataBase(self.integrdataPlot)
        dbdockwidget = qt.QDockWidget("Integrated data")
        dbdockwidget.setWidget(self.database)
        
        self.integrdataPlot.addDockWidget(qt.Qt.RightDockWidgetArea,dbdockwidget)
        
        self.roPkIntegrTab = RockingPeakIntegrator(self.database)
        
        

        self.scanSelector.sigImageNoChanged.connect(self._onSliderValueChanged)

        self.scanSelector.sigImagePathChanged.connect(self._onImagePathChanged)
        self.scanSelector.sigScanChanged.connect(self._onScanChanged)

        self.scanSelector.showMaxAct.toggled.connect(self._onMaxToggled)
        self.scanSelector.showSumAct.toggled.connect(self._onSumToggled)
        
        self.scanSelector.sigROIChanged.connect(self.updateROI)
        self.scanSelector.sigROIintegrate.connect(self.integrateROI)
        self.scanSelector.sigSearchHKL.connect(self.onSearchHKLforStaticROI)
        
        self.scanSelector.excludeImageAct.toggled.connect(self._onToggleExcludeImage)
        

        toolbar = self.scanSelector.getScanToolbar()

        
        self.centralPlot.addToolBar(qt.Qt.BottomToolBarArea,toolbar)
        
        maincentralwidget.addTab(self.centralPlot,"Scan Image browser")
        maincentralwidget.addTab(self.integrdataPlot,"ROI integrated data")
        maincentralwidget.addTab(self.roPkIntegrTab, "Rocking scan integrate")
        
        self.setCentralWidget(maincentralwidget)
        
        
        # Create the object controlling the ROIs and set it up
        self.roiManager = RegionOfInterestManager(self.centralPlot)
        self.roiManager.setColor('pink')  # Set the color of ROI
        
        #self.roiTable = RegionOfInterestTableWidget()
        #self.roiTable.setRegionOfInterestManager(self.roiManager)
        
        self.rocking_rois = []
        
        self.roiS1 = RectangleBgROI()
        self.roiS1.setLineWidth(2)
        self.roiS1.setLineStyle('-')
        self.roiS1.setColor('red')
        self.roiS1.setBgStyle('pink', '-', 2.)
        self.roiS1.setVisible(True)
        self.roiS1.setGeometry(origin=(0, 0), size=(0, 0))
        # S1 is also fixed roi, editing enabled when static ROI:
        self.roiS1.sigEditingFinished.connect(self._onStaticROIedited)
        
        self.roiS2 = RectangleBgROI()
        self.roiS2.setLineWidth(2)
        self.roiS2.setLineStyle('-')
        self.roiS2.setColor('red')
        self.roiS2.setBgStyle('pink', '-', 2.)
        self.roiS2.setVisible(True)
        self.roiS2.setGeometry(origin=(0, 0), size=(0, 0))
        self.roiManager.addRoi(self.roiS1,useManagerColor=False)
        self.roiManager.addRoi(self.roiS2,useManagerColor=False)
        
        #self.reflTable.view._model.dataChanged.connect(printmodel)
        #self.reflTable.setArrayData(np.array([0,0,0,0,10,10],dtype=np.float64))
        ubDock = qt.QDockWidget("Reciprocal space navigation")
        ubDock.setAllowedAreas(qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea | qt.Qt.BottomDockWidgetArea)
        
        self.reflectionSel = QReflectionSelector(self.centralPlot, self.ubcalc, self)
        self.reflectionSel.sigQueryImageChange.connect(self._onChangeImage)
        
        self.ubcalc.setReflectionHandler(self.getReflections)
        
        self.ubcalc.sigPlottableMachineParamsChanged.connect(self._onPlotMachineParams)
        self.ubcalc.sigReplotRequest.connect(self.updatePlotItems)
        self.allimgsum = None
        self.allimgmax = None

        ubWidget.addWidget(self.reflectionSel)
        ubWidget.addWidget(self.ubcalc)
        
        
        #ubWidget.setLayout(ubLayout)
        ubDock.setWidget(ubWidget)
        self.centralPlot.addDockWidget(qt.Qt.RightDockWidgetArea,ubDock)
        
        
        
        menu_bar = qt.QMenuBar() 
        file = menu_bar.addMenu("&File")
        file.addAction(self.scanSelector.openFileAction)
        file.addAction(self.scanSelector.refreshFileAction)
        file.addAction(self.scanSelector.closeFileAction)
        file.addSeparator()
        
        self.folderToScan = file.addAction("Generate scan from images")
        self.folderToScan.triggered.connect(self._onLoadScanFromImages)
        file.addSeparator()

        self.loadImagesAct = file.addAction("reload images")
        self.loadImagesAct.triggered.connect(self._onLoadAll)
        
        config_menu =  menu_bar.addMenu("&Config")
        loadConfigAct = qt.QAction("Load config",self) # connected with UBCalculator below
        #loadXtalAct = qt.QAction("Load Crystal file",self)
        machineParamsAct = qt.QAction("Machine parameters",self)
        machineParamsAct.setCheckable(True)
        xtalParamsAct = qt.QAction("Crystal parameters",self)
        xtalParamsAct.setCheckable(True)
        cpucountAct = qt.QAction("Set CPU count",self)
        
        
        loadConfigAct.triggered.connect(self.ubcalc._onLoadConfig)        
        machineParamsAct.toggled.connect(lambda checked: self.ubcalc.machineDialog.setVisible(checked))
        self.ubcalc.machineDialog.sigHide.connect(lambda : machineParamsAct.setChecked(False))
        
        xtalParamsAct.toggled.connect(lambda checked: self.ubcalc.xtalDialog.setVisible(checked))
        self.ubcalc.xtalDialog.sigHide.connect(lambda : xtalParamsAct.setChecked(False))
        
        cpucountAct.triggered.connect(self._onSelectCPUcount)
        
        self.autoLoadAct = qt.QAction("Auto load scans",self)
        self.autoLoadAct.setCheckable(True)
        self.autoLoadAct.setChecked(True)
        
        self.showExcludedImagesAct = qt.QAction("Excluded images",self)
        self.showExcludedImagesAct.setCheckable(True)
        self.showExcludedImagesAct.toggled.connect(lambda visible : self.excludedImagesDialog.setVisible(visible))
        

        
        config_menu.addAction(loadConfigAct)
        #config_menu.addAction(loadXtalAct)
        config_menu.addSeparator()
        config_menu.addAction(machineParamsAct)
        config_menu.addAction(xtalParamsAct)
        config_menu.addSeparator()
        config_menu.addAction(cpucountAct)
        config_menu.addAction(self.autoLoadAct)
        config_menu.addAction(self.showExcludedImagesAct)
        
        view_menu = menu_bar.addMenu("&View")
        showRefReflectionsAct = view_menu.addAction("reference reflections")
        showRefReflectionsAct.setCheckable(True)
        showRefReflectionsAct.setChecked(True)
        showRefReflectionsAct.toggled.connect(lambda checked: self.reflectionSel.setReferenceReflectionsVisible(checked))
        
        showBraggAct = view_menu.addAction("allowed Bragg reflections")
        showBraggAct.setCheckable(True)
        showBraggAct.setChecked(False)
        showBraggAct.toggled.connect(self.onShowBragg)
        
        showROIAct = view_menu.addAction("show ROI")
        showROIAct.setCheckable(True)
        showROIAct.setChecked(False)
        showROIAct.toggled.connect(self.onShowROI)
        self.roivisible = False
        #self.scanSelector.showROICheckBox.addAction(showROIAct)
        
        self.showCTRreflAct = view_menu.addAction("CTR reflections")
        self.showCTRreflAct.setCheckable(True)
        self.showCTRreflAct.setChecked(False)
        self.showCTRreflAct.toggled.connect(self.onShowCTRreflections)
        self.reflectionsVisible = False
        
        self.showMachineParamsAct = view_menu.addAction("machine parameters")
        self.showMachineParamsAct.setCheckable(True)
        self.showMachineParamsAct.setChecked(False)
        self.showMachineParamsAct.toggled.connect(self._onPlotMachineParams)
        
        view_menu.addSeparator()
        
        saveBraggAct = view_menu.addAction("save allowed Bragg reflections")
        saveBraggAct.setCheckable(False)
        saveBraggAct.triggered.connect(self.saveBraggRefl)
        if console:
            view_menu.addSeparator()
            
            custom_banner = f"""orGUI v. {__version__} console 
Available variables:
orgui : top level gui
ub : gui for UB matrix and angle calculations 
"""
            
            self.console_dockwidget = console.IPythonDockWidget(self, {'orgui': self, 'ub': self.ubcalc}, custom_banner, "orGUI console")
            
            self.console_dockwidget.setAllowedAreas(qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea | qt.Qt.BottomDockWidgetArea)
            self.tabifyDockWidget(selectorDock,self.console_dockwidget)
            #self.addDockWidget(qt.Qt.LeftDockWidgetArea,self.console_dockwidget)
            self.console_dockwidget.setVisible(False)
            consoleViewAct = self.console_dockwidget.toggleViewAction()
            view_menu.addAction(consoleViewAct)
        
        
        ##############################
        
        editUAct = qt.QAction("Edit orientation matrix",self)
        editUAct.setCheckable(True)
        editUAct.toggled.connect(lambda checked: self.ubcalc.ueditDialog.setVisible(checked))
        self.ubcalc.ueditDialog.sigHide.connect(lambda : editUAct.setChecked(False))
        
        calcCTRsAvailableAct = qt.QAction("Calculate available CTRs",self)
        calcCTRsAvailableAct.triggered.connect(self._onCalcAvailableCTR)
        rs = menu_bar.addMenu("&Reciprocal space")
        rs.addAction(calcCTRsAvailableAct)
        rs.addAction(editUAct)
        
        simul = menu_bar.addMenu("&Simulation")
        
        createScanAct = simul.addAction("Create dummy scan")
        createScanAct.triggered.connect(self._onCreateScan)
    
        
        helpmenu = menu_bar.addMenu("&Help")
        
        diffractAct = helpmenu.addAction("Diffraction geometry")
        diffractAct.triggered.connect(self._onShowDiffractionGeometry)
        
        helpmenu.addSeparator()
        
        aboutAct = helpmenu.addAction("About")
        aboutAct.triggered.connect(self._onShowAbout)
        

        aboutQtAct = helpmenu.addAction("About Qt")
        aboutQtAct.triggered.connect(lambda : qt.QMessageBox.aboutQt(self))
        
        self.setMenuBar(menu_bar)

    def get_rocking_coordinates(self, H_0=None, H_1=None, maxValue=None,step_width=None, **kwargs):
        # going back to the more universal integration along H_0 + s*H_1 positions 
        if H_0 is None:
            H_0 = self.scanSelector.ro_H_0_dialog.get_hkl()
        if H_1 is None:
            H_1 = self.scanSelector.ro_H_1_dialog.get_hkl() # default to CTR integration
        if step_width is None:
            step_width = self.scanSelector.roscanDeltaS.value()
        if maxValue is None:
            maxValue = self.scanSelector.roscanMaxS.value()
            
        dc = self.ubcalc.detectorCal
        xoffset = kwargs.get('xoffset', self.scanSelector.offsetx.value())
        yoffset = kwargs.get('yoffset', self.scanSelector.offsety.value())

        step_nr = round(maxValue/step_width) + 1
        s_points = np.linspace(0,maxValue,step_nr) 
        
        hkl_desired = np.outer(H_1, s_points).T + H_0 # F contiguous is faster in anglesToHKL
        
        refldict = self.ubcalc.calcReflection(hkl_desired) # F contiguous is faster

        ymask1 = np.logical_and(refldict['xy_1'][...,1] >= 0, refldict['xy_1'][...,1] < dc.detector.shape[0])
        xmask1 = np.logical_and(refldict['xy_1'][...,0] >= 0, refldict['xy_1'][...,0] < dc.detector.shape[1])
        yxmask1 = np.logical_and(xmask1,ymask1)
    
        ymask2 = np.logical_and(refldict['xy_2'][...,1] >= 0, refldict['xy_2'][...,1] < dc.detector.shape[0])
        xmask2 = np.logical_and(refldict['xy_2'][...,0] >= 0, refldict['xy_2'][...,0] < dc.detector.shape[1])
        yxmask2 = np.logical_and(xmask2,ymask2)
        
        refldict['mask_1'] = yxmask1
        refldict['mask_2'] = yxmask2
        refldict['s'] = s_points
        
        if xoffset != 0. or yoffset != 0.:
            #warnings.warn("Nonzero pixel offset selected. Experimental feature! Angles and hkl are incorrect!!!")
            refldict['xy_1'][..., 0] += xoffset
            refldict['xy_1'][..., 0] += xoffset
            refldict['xy_2'][..., 1] += yoffset
            refldict['xy_2'][..., 1] += yoffset
        refldict['H_0'] = H_0
        refldict['H_1'] = H_1
        return refldict
        
    def intkeys_rocking(self, refldict, **kwargs):
        vsize = kwargs.get('vsize' ,int(self.scanSelector.vsize.value()))
        hsize = kwargs.get('hsize' ,int(self.scanSelector.hsize.value()))
        apply_mask = kwargs.get('mask' ,True)
        autoROIVsize = kwargs.get('autovsize', self.scanSelector.autoROIVsize.isChecked())
        autoROIHsize = kwargs.get('autohsize', self.scanSelector.autoROIHsize.isChecked())
        if self.scanSelector.intersS1Act.isChecked():
            intersect = 1
        elif self.scanSelector.intersS2Act.isChecked():
            intersect = 2
        else:
            intersect = 1 # default
        intersect = kwargs.get('intersect', intersect)
        
        xy = refldict['xy_%s' % int(intersect)]
        if apply_mask:
            xy = xy[refldict['mask_%s' % int(intersect)]]
        
        step_nr = xy.shape[0]
        if step_nr == 0:
            return {'center' : [], 'vsize' : vsize, 'hsize': hsize}
        if step_nr > 1:
            if autoROIVsize:
                #dist_in_pixels = np.abs(xy[0][1] - xy[-1][1])
                dist_in_pixels = np.median(np.abs(np.diff(xy[:,1])))
                #roi_vlength = np.ceil(dist_in_pixels/step_nr)
                vsize = int(np.ceil(dist_in_pixels))
                
            if autoROIHsize:
                #dist_in_pixels = np.abs(xy[0][1] - xy[-1][1])
                dist_in_pixels = np.median(np.abs(np.diff(xy[:,0])))
                #roi_vlength = np.ceil(dist_in_pixels/step_nr)
                hsize = int(np.ceil(dist_in_pixels))

        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape
        
        coord_restr = np.clip( xy, [0,0], [dethsize, detvsize])

        vhalfsize = vsize // 2
        hhalfsize = hsize // 2
        fromcoords = np.round(coord_restr - np.array([hhalfsize, vhalfsize]))
        tocoords = np.round(coord_restr + np.array([hhalfsize, vhalfsize]))
        
        if hsize % 2:
            remainder_mask = coord_restr[:,0] % 1 < 0.5
            tocoords[remainder_mask, 0] += 1
            fromcoords[~remainder_mask, 0] -= 1
        if vsize % 2:
            remainder_mask = coord_restr[:,1] % 1 < 0.5
            tocoords[remainder_mask, 1] += 1
            fromcoords[~remainder_mask, 1] -= 1
                
        fromcoords = np.clip( fromcoords, [0,0], [dethsize, detvsize])
        tocoords = np.clip( tocoords, [0,0], [dethsize, detvsize])
        
        locations = []
        for roifrom, toroi in zip(fromcoords, tocoords): # any way to do this with ndarray operations?
            locations.append(tuple(slice(int(fromcoord), int(tocoord)) for fromcoord, tocoord in zip(roifrom,toroi)))
        return {'center' : locations, 'vsize' : vsize, 'hsize': hsize}
        
    def intbkgkeys_rocking(self, refldict, **kwargs):
        left = kwargs.get( 'left' ,int(self.scanSelector.left.value()))
        right = kwargs.get( 'right' ,int(self.scanSelector.right.value()))
        top = kwargs.get( 'top' ,int(self.scanSelector.top.value()))
        bottom = kwargs.get( 'bottom' ,int(self.scanSelector.bottom.value()))
        
        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape
        
        roi_dict = self.intkeys_rocking(refldict, **kwargs)
        crois = roi_dict['center']
        leftrois = []
        rightrois = []
        toprois = []
        bottomrois = []
        
        for croi in crois:
            leftrois.append((slice(int(np.clip(croi[0].start - left, 0, dethsize)), croi[0].start), croi[1]))
            rightrois.append((slice(croi[0].stop,int(np.clip(croi[0].stop + right, 0, dethsize))), croi[1]))
            toprois.append((croi[0], slice(int(np.clip(croi[1].start - top, 0, detvsize)), croi[1].start)))
            bottomrois.append((croi[0], slice(croi[1].stop, int(np.clip(croi[1].stop + bottom,0,detvsize)) )))
        roi_dict['left'] = leftrois
        roi_dict['right'] = rightrois
        roi_dict['top'] = toprois
        roi_dict['bottom'] = bottomrois
        return roi_dict

    def rocking_extraction(self):

        if self.fscan is None: #or isinstance(self.fscan, SimulationScan):
            qt.QMessageBox.warning(self, "No scan loaded", "Cannot integrate scan: No scan loaded.")
            return
        refldict = self.get_rocking_coordinates()

        if self.scanSelector.intersS1Act.isChecked():
            intersect = 1
        elif self.scanSelector.intersS2Act.isChecked():
            intersect = 2
        else:
            intersect = 1 # default
        mask = refldict['mask_%s' % intersect]
        xy = refldict['xy_%s' % intersect][mask]
        
        refldict['angles'] = refldict['angles_%s' % intersect][mask]
        
        roi_keys = self.intbkgkeys_rocking(refldict)
        hkl_del_gam = self.getStaticROIparams(xy)

        self.rocking_integrate(xy, roi_keys, hkl_del_gam, refldict)
            


    def rocking_integrate(self,xylist, rois, hkl_del_gam, refldict):
        try:
            image = self.fscan.get_raw_img(0)
        except Exception as e:
            print("no images found! %s" % e)
            return
        if self.database.nxfile is None:
            print("No database available")
            return
        dc = self.ubcalc.detectorCal
        
        imgmask = None
        
        if self.scanSelector.useMaskBox.isChecked():
            if self.centralPlot.getMaskToolsDockWidget().getSelectionMask() is None:
                btn = qt.QMessageBox.question(self,"No mask available","""No mask was selected with the masking tool.
    Do you want to continue without mask?""")
                if btn != qt.QMessageBox.Yes:
                    return
            else:
                imgmask = self.centralPlot.getMaskToolsDockWidget().getSelectionMask() > 0.
        
        corr = self.scanSelector.useSolidAngleBox.isChecked() or\
            self.scanSelector.usePolarizationBox.isChecked()
        
        if corr:
            C_arr = np.ones(dc.detector.shape,dtype=np.float64)
            if self.scanSelector.useSolidAngleBox.isChecked():
                C_arr /= dc.solidAngleArray()
            if self.scanSelector.usePolarizationBox.isChecked():
                C_arr /= dc.polarization(factor=dc._polFactor,axis_offset=dc._polAxis)


        def fill_counters(image,pixelavail, key, bkgkey):
            
            cimg = image[key[::-1]]
            
            # !!!!!!!!!! add mask here  !!!!!!!!!
            croi = np.nansum(cimg)
            cpixel = np.nansum(pixelavail[key[::-1]])
            bgroi = 0.
            bgpixel = 0.
            for bg in bkgkey:
                bgimg = image[bg[::-1]]
                bgroi += np.nansum(bgimg)
                bgpixel += np.nansum(pixelavail[bg[::-1]])


            return (croi, cpixel, bgroi, bgpixel)
        
        hkl_del_gam_1 = hkl_del_gam[0] # needed to initialize integration 

        # initialize 1d np arrays for storing roi integration counters for all images
        croi1_a = np.zeros_like(hkl_del_gam_1.shape[0],dtype=np.float64)
        cpixel1_a = np.zeros_like(hkl_del_gam_1.shape[0],dtype=np.float64)
        bgroi1_a = np.zeros_like(hkl_del_gam_1.shape[0],dtype=np.float64)
        bgpixel1_a = np.zeros_like(hkl_del_gam_1.shape[0],dtype=np.float64)

        # initialize 2d np array to store roi integration counters together for all images/ROIs
        croi1_all = np.zeros((hkl_del_gam_1.shape[0],) + (xylist.shape[0],))
        cpixel1_all = np.zeros((hkl_del_gam_1.shape[0],) + (xylist.shape[0],))
        bgroi1_all = np.zeros((hkl_del_gam_1.shape[0],) + (xylist.shape[0],))
        bgpixel1_all = np.zeros((hkl_del_gam_1.shape[0],) + (xylist.shape[0],))
        
        progress = qt.QProgressDialog("Integrating images","abort",0,len(self.fscan),self)
        progress.setWindowModality(qt.Qt.WindowModal)
        
        if HAS_ACCEL:
            if imgmask is not None:
                mask = np.ascontiguousarray(imgmask, dtype=bool)
            else:
                mask = np.zeros(image.img.shape, dtype=bool)
            if corr:
                C_arr = np.ascontiguousarray(C_arr, dtype=np.float64)
            else:
                C_arr = np.ones(image.img.shape, dtype=bool)
                
            roi_lists_numba = []
            for roiname in ['center', 'left', 'right', 'top', 'bottom']:
                roi_list = [] 
                for r in rois[roiname]:
                    roi_list.append(np.array([[r[0].start , r[0].stop], [r[1].start , r[1].stop]]))
                roi_list = np.ascontiguousarray(np.stack(roi_list), dtype=np.int64)
                roi_lists_numba.append(roi_list)
                
            def sumImage(i):
                image = np.ascontiguousarray(self.fscan.get_raw_img(i).img, dtype=np.float64) # unlocks gil during file read

                all_counters = np.zeros((roi_lists_numba[0].shape[0],) + (4,)) # need gil for python object creation
                _roi_sum_accel.processImage(image, mask, C_arr, *roi_lists_numba, all_counters) # numba nopython and nogil mode
                return all_counters
            
        else:
            

            def sumImage(i):
                image = self.fscan.get_raw_img(i).img.astype(np.float64)
                if imgmask is not None:
                    image[imgmask] = np.nan
                    pixelavail = (~imgmask).astype(np.float64)
                else:
                    pixelavail = np.ones_like(image)
                if corr:
                    image *= C_arr

                all_counters1 = np.zeros((xylist.shape[0],) + (4,))

                for crnr in range(xylist.shape[0]):

                    # set ROI (moved to rocking-function)

                    # get roi
                    key = rois['center'][crnr]
                    bgkey = [rois['left'][crnr], rois['right'][crnr], rois['top'][crnr], rois['bottom'][crnr]]
                    # fill counters
                    counters1 = fill_counters(image,pixelavail,key, bgkey)

                    all_counters1[crnr] = counters1


                return all_counters1
            

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.numberthreads) as executor: # speedup only for the file reads 
            futures = {}
            for i in range(len(self.fscan)):
                futures[executor.submit(sumImage, i)] = i
            
            status = "no error"
            for f in concurrent.futures.as_completed(futures): # iteration over jobs
                try:
                    for j in range(len(f.result())): # iteration over ROIs
                        (croi1, cpixel1, bgroi1, bgpixel1) = f.result()[j]
                        i = futures[f]
                        croi1_all[i][j] = croi1
                        cpixel1_all[i][j] = cpixel1
                        bgroi1_all[i][j] = bgroi1
                        bgpixel1_all[i][j] = bgpixel1
                    progress.setValue(futures[f])
                    del f
                except concurrent.futures.CancelledError:
                    del f
                except Exception as e:
                    print("Cannot read image:\n%s" % traceback.format_exc())
                    print("Cancel to avoid memory leak")
                    [f.cancel() for f in futures]
                    status = 'error'
                    del f
                if progress.wasCanceled():
                    status = 'error'
                    [f.cancel() for f in futures]
                    
                    break
                #gc.collect()

        progress.setValue(len(self.fscan))
        if status == 'error':
            return

        currentPlotCount = len(self.integrdataPlot.getAllCurves())
        numberOfNewPlots = xylist.shape[0]
        maxAmountOfPlots = 30
        plotOnlyNth = (numberOfNewPlots // max((maxAmountOfPlots-currentPlotCount),1)) + 1

        #print('Number of integration curves: ' + str(numberOfPlots))
        #print('We can plot every ' + str(plotOnlyNth) + '-th curve.' )
        
        ro_name = "rocking_[%.2f %.2f %.2f]_H0_[%.2f %.2f %.2f]_H1" % (*refldict['H_0'], *refldict['H_1'])
        suffix = ''
        i = 0
        while(self.activescanname + "/measurement/" + ro_name + suffix in self.database.nxfile):
            suffix = "_%s" % i
            i += 1
        ro_name = ro_name + suffix
        
        auxcounters = {"@NX_class": u"NXcollection"}
        for auxname in self.fscan.auxillary_counters:
            if hasattr(self.fscan, auxname):
                cntr = getattr(self.fscan, auxname)
                if cntr is not None:
                    auxcounters[auxname] = cntr
            
        if hasattr(self.fscan, "title"):
            title = str(self.fscan.title)
        else:
            title = u"%s-scan" % self.fscan.axisname
        
        mu, om = self.getMuOm()
        if len(np.asarray(om).shape) == 0:
            om = np.full_like(mu,om)
        if len(np.asarray(mu).shape) == 0:
            mu = np.full_like(om,mu)
                    
        data = {self.activescanname:{
                    "instrument": {
                        "@NX_class": u"NXinstrument",
                        "positioners": {
                            "@NX_class": u"NXcollection",
                            self.fscan.axisname: self.fscan.axis
                        }
                    },
                    "auxillary" : auxcounters,
                    "measurement": {
                        "@NX_class": u"NXentry",
                        "@default": ro_name ,
                        ro_name : {
                            "@NX_class": u"NXentry",
                            "@default": "rois",
                            "@orgui_meta": u"rocking",
                            "rois" : {
                                "@NX_class": u"NXcollection",
                                "@default": None,
                                "@orgui_meta": u"roi rocking",
                            }
                        }
                    },
                    "title":u"%s" % title,
                    "@NX_class": u"NXentry",
                    "@default": u"measurement/%s" % ro_name,
                    "@orgui_meta": u"scan"
                }
            }

        #plot and save data in database
        for d in range(croi1_all.shape[1]):

            hkl_del_gam_1 = hkl_del_gam[d]

            croi1_a = croi1_all[...,d]
            cpixel1_a = cpixel1_all[...,d]
            bgroi1_a = bgroi1_all[...,d]
            bgpixel1_a = bgpixel1_all[...,d]
            
        
            if np.any(bgpixel1_a):
                croibg1_a = croi1_a - (cpixel1_a/bgpixel1_a) * bgroi1_a
                croibg1_err_a = np.sqrt(croi1_a + ((cpixel1_a/bgpixel1_a)**2)  * bgroi1_a)
            else:
                croibg1_a = croi1_a
                croibg1_err_a = np.sqrt(croi1_a)
                

            rod_mask1 = np.isfinite(croibg1_a)

            
            s1_masked = hkl_del_gam_1[:,5][rod_mask1]
            
            croibg1_a_masked = croibg1_a[rod_mask1]

            croibg1_err_a_masked = croibg1_err_a[rod_mask1]


            # save data
            
            x, y = xylist[d] 
            name1 = "rocking_%.5fs_[%.2f %.2f %.2f]" % (refldict['s'][d], *(refldict['H_1']*refldict['s'][d] + refldict['H_0']))
            
            alpha1, delta1, gamma1, omega1, chi1, phi1 = refldict['angles'][d]
            sixc_angles_hkl = {
                    "@NX_class": u"NXpositioner",
                    "alpha" : np.rad2deg(alpha1),
                    "omega" :  np.rad2deg(omega1),
                    "theta" :  np.rad2deg(-1*omega1),
                    "delta" : np.rad2deg(delta1),
                    "gamma" :  np.rad2deg(gamma1),
                    "chi" :  np.rad2deg(chi1),
                    "phi" :  np.rad2deg(phi1),
                    "@unit" : u"deg"
            }
            
            traj1 = {
                "@direction" : u"Rocking scan at fixed pixel location along H_1*s + H_0 in reciprocal space",
                "@NX_class": u"NXcollection",
                "axis" : hkl_del_gam_1[:,5],
                "s" : refldict['s'][d], 
                "H_1" : refldict['H_1'],
                "H_0" : refldict['H_0'],
                "HKL" : refldict['H_1']*refldict['s'][d] + refldict['H_0'],
                "HKL_sixc_angles" : sixc_angles_hkl
            }
            
            suffix = ''
            i = 0

            while(self.activescanname + "/measurement/" + ro_name + "/" + name1 + suffix in self.database.nxfile):
                suffix = "_%s" % i
                i += 1
                
            availname1 = name1 + suffix

            x_coord1_a = xylist[:,0]
            y_coord1_a = xylist[:,1]
                        
            datas1 = {
                "@NX_class": u"NXdata",
                "sixc_angles": {
                    "@NX_class": u"NXpositioner",
                    "alpha" : np.rad2deg(mu),
                    "omega" :  np.rad2deg(om),
                    "theta" :  np.rad2deg(-1*om),
                    "delta" : np.rad2deg(hkl_del_gam_1[:,3]),
                    "gamma" :  np.rad2deg(hkl_del_gam_1[:,4]),
                    "chi" :  np.rad2deg(self.ubcalc.chi),
                    "phi" :  np.rad2deg(self.ubcalc.phi),
                    "@unit" : u"deg"
                },
                "hkl": {
                    "@NX_class": u"NXcollection",
                    "h" :  hkl_del_gam_1[:,0],
                    "k" :  hkl_del_gam_1[:,1],
                    "l" : hkl_del_gam_1[:,2]
                },
                "counters":{
                    "@NX_class": u"NXdetector",
                    "croibg"  : croibg1_a,
                    "croibg_errors" :  croibg1_err_a,
                    "croi" :  croi1_a,
                    "bgroi"  : bgroi1_a,
                    "croi_pix"  : cpixel1_a,
                    "bgroi_pix" : bgpixel1_a
                },
                "pixelcoord": {
                    "@NX_class": u"NXdetector",
                    "x" : x_coord1_a,
                    "y"  : y_coord1_a
                },
                "trajectory" : traj1,
                "@signal" : u"counters/croibg",
                "@axes": u"trajectory/axis",
                "@title": self.activescanname + "_" + availname1,
                "@orgui_meta": u"roi rocking"
            }
                
            data[self.activescanname]["measurement"][ro_name]["rois"]["@default"] = availname1
            if np.any(cpixel1_a > 0.):
                data[self.activescanname]["measurement"][ro_name]["rois"][availname1] = datas1
                if d % plotOnlyNth == 0:
                    self.integrdataPlot.addCurve(s1_masked,croibg1_a_masked,legend=self.activescanname + "_" + availname1,
                                                    xlabel="trajectory/s", ylabel="counters/croibg", yerror=croibg1_err_a_masked)

        
        
        # lets keep legacy data structure for now
            
        data_2d_structured = {self.activescanname:{
                    "instrument": {
                        "@NX_class": u"NXinstrument",
                        "positioners": {
                            "@NX_class": u"NXcollection",
                            self.fscan.axisname: self.fscan.axis
                        }
                    },
                    "auxillary" : auxcounters,
                    "measurement": {
                        "@NX_class": u"NXentry",
                        "@default": ro_name ,
                        ro_name : {
                            "@NX_class": u"NXentry",
                            "@default": "rois",
                            "@orgui_meta": u"rocking"
                        }
                    },
                    "title":u"%s" % title,
                    "@NX_class": u"NXentry",
                    "@default": u"measurement/%s" % ro_name,
                    "@orgui_meta": u"scan"
                }
            }
        alpha = []; theta = []; delta = []; gamma = []; chi = []; phi = []; omega = []; 
        alpha_pk = []; theta_pk = []; delta_pk = []; gamma_pk = []; chi_pk = []; phi_pk = []; omega_pk = []; 
        x = []; y = []; h = []; k = []; l = []; 
        croibg = []; croibg_errors = []; croi = []; bgroi = []; croi_pix = []; bgroi_pix = []; 
        axis = []; s = []; H_0 = []; H_1 = []; HKL = []
        
        #from IPython import embed; embed()

        for sc in data[self.activescanname]["measurement"][ro_name]["rois"]:
            if sc.startswith('@'):
                continue
            try:
                dsc = data[self.activescanname]["measurement"][ro_name]["rois"][sc]
                
                
                
                # 2D arrays
                alpha.append(dsc["sixc_angles"]["alpha"])
                theta.append(dsc["sixc_angles"]["theta"])
                delta.append(dsc["sixc_angles"]["delta"])
                gamma.append(dsc["sixc_angles"]["gamma"])
                chi.append(dsc["sixc_angles"]["chi"])
                phi.append(dsc["sixc_angles"]["phi"])
                omega.append(dsc["sixc_angles"]["omega"])
                
                # 2D arrays
                h.append(dsc["hkl"]["h"])
                k.append(dsc["hkl"]["k"])
                l.append(dsc["hkl"]["l"])
                
                # 2D arrays
                croibg.append(dsc["counters"]["croibg"])
                croibg_errors.append(dsc["counters"]["croibg_errors"])
                croi.append(dsc["counters"]["croi"])
                bgroi.append(dsc["counters"]["bgroi"])
                croi_pix.append(dsc["counters"]["croi_pix"])
                bgroi_pix.append(dsc["counters"]["bgroi_pix"])
                
                # 2D arrays
                x.append(dsc["pixelcoord"]["x"])
                y.append(dsc["pixelcoord"]["y"])
                
                axis.append(dsc["trajectory"]["axis"])
                s.append(dsc["trajectory"]["s"]) # 1D array
                H_1.append(dsc["trajectory"]["H_1"])
                H_0.append(dsc["trajectory"]["H_0"])
                HKL.append(dsc["trajectory"]["HKL"])
                
                # 1d Array
                alpha_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["alpha"])
                theta_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["theta"])
                delta_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["delta"])
                gamma_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["gamma"])
                chi_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["chi"])
                phi_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["phi"])
                omega_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["omega"])
            except Exception as e:
                from IPython import embed; embed()
                sys.exit(0)

        
        rois = {
                "@NX_class": u"NXcollection",
                "@default": "croibg",
                "@orgui_meta": u"roi rocking",
                "alpha" : np.vstack(alpha),
                "theta" : np.vstack(theta),
                "delta" : np.vstack(delta),
                "gamma" : np.vstack(gamma),
                "chi" : np.vstack(chi),
                "phi" : np.vstack(phi),
                "omega" : np.vstack(omega),
                "h" : np.vstack(h),
                "k" : np.vstack(k),
                "l" : np.vstack(l),
                "croibg" : np.vstack(croibg),
                "croibg_errors" : np.vstack(croibg_errors),
                "croi" : np.vstack(croi),
                "bgroi" : np.vstack(bgroi),
                "croi_pix" : np.vstack(croi_pix),
                "bgroi_pix" : np.vstack(bgroi_pix),
                "x" : np.vstack(x),
                "y" : np.vstack(y),
                "axis" : np.vstack(axis),
                "H_1" : np.vstack(H_1),
                "H_0" : np.vstack(H_0),
                "HKL_pk" : np.vstack(HKL),
                "s" : np.array(s),
                "alpha_pk" : np.array(alpha_pk),
                "theta_pk" : np.array(theta_pk),
                "delta_pk" : np.array(delta_pk),
                "gamma_pk" : np.array(gamma_pk),
                "chi_pk" : np.array(chi_pk),
                "phi_pk" : np.array(phi_pk),
                "omega_pk" : np.array(omega_pk)
            }
        scsize = np.array(s).shape[0]
        for t in rois:
            if t.startswith('@'):
                continue
            if rois[t].shape[0] != scsize:
                from IPython import embed; embed()
                sys.exit(0)
            
        data_2d_structured[self.activescanname]["measurement"][ro_name]["rois"] = rois
                    
        self.database.add_nxdict(data_2d_structured)
        
                
    def updatePlotItems(self, recalculate=True):
        if self.roivisible:
            try:
                self.updateROI()
            except Exception:
                pass
        
        if self.reflectionsVisible:
            if recalculate:
                try: 
                    hkm = self.calculateAvailableCTR()
                    hk = np.unique(hkm[:,:2],axis=0)
                    H_0 = np.hstack((hk, np.zeros((hk.shape[0],1))))
                    H_1 = np.array([0,0,1])
                    self.reflectionsToDisplay = H_0, H_1
                except Exception:
                    self.showCTRreflAct.setChecked(False)
                    self.reflectionsVisible = False
            self.updateReflections()

        if self.reflectionSel.showBraggReflections:
            try:
                self.calcBraggRefl()
            except:
                pass
        
        
    def onShowBragg(self,visible):
        try:
            self.reflectionSel.setBraggReflectionsVisible(visible)
            self.calcBraggRefl()
        except Exception:
            qutils.warning_detailed_message(self, "Cannot show show Bragg reflections", "Cannot show Bragg reflections", traceback.format_exc())
            #qt.QMessageBox.critical(self,"Cannot show show Bragg reflections", "Cannot Cannot show Bragg reflections:\n%s" % traceback.format_exc())
        
    def onShowROI(self,visible):
        self.roivisible = visible
        try:
            self.updateROI()
        except Exception:
            qutils.warning_detailed_message(self, "Cannot show ROI", "Cannot show ROI", traceback.format_exc())
            #qt.QMessageBox.critical(self,"Cannot show ROI", "Cannot Cannot show ROI:\n%s" % traceback.format_exc())
            
    def onShowCTRreflections(self,visible):
        self.reflectionsVisible = visible
        if self.reflectionsVisible:
            try: 
                hkm = self.calculateAvailableCTR()
            except Exception:
                qutils.warning_detailed_message(self, "Cannot calculate CTR locations", "Cannot calculate CTR locatons", traceback.format_exc())
                #qt.QMessageBox.critical(self,"Cannot calculate CTR locatons", "Cannot calculate CTR locatons:\n%s" % traceback.format_exc())
                self.showCTRreflAct.setChecked(False)
                self.reflectionsVisible = False
                return
            hk = np.unique(hkm[:,:2],axis=0)
            H_0 = np.hstack((hk, np.zeros((hk.shape[0],1))))
            H_1 = np.array([0,0,1])
            
            self.reflectionsToDisplay = H_0, H_1
        self.updateReflections()
        
    def _onShowAbout(self):
        dial = AboutDialog(self, __version__)
        dial.exec()
#        messageStr = """Copyright (c) 2020-2024 Timo Fuchs, published under MIT License
#        <br> <br>
#orGUI: Orientation and Integration with 2D detectors (1.0.0).<br>
#Zenodo. <a href=\"https://doi.org/10.5281/zenodo.12592485\">https://doi.org/10.5281/zenodo.12592485</a> <br> <br> 
#New software updates will be published under <a href=\"https://doi.org/10.5281/zenodo.12592485\">Zenodo</a>.
#<br> <br>
#Help requests can be send via Email to Timo Fuchs. 
#<br> <br>
#"orGUI" was developed during the PhD work of Timo Fuchs,
#within the group of Olaf Magnussen.
#"""
#        msg0 = qt.QMessageBox(self)
#        msg0.setWindowTitle("About orGUI")
#        msg0.setText(messageStr)
#        msg0.setTextInteractionFlags(qt.Qt.TextBrowserInteraction)
#        msg0.setTextFormat(qt.Qt.RichText)
#        msg0.exec()
        
    def _onShowDiffractionGeometry(self):
        if hasattr(self, 'diffractometerdialog'):
            self.diffractometerdialog.show()
        else:
            self.diffractometerdialog = QDiffractometerImageDialog()
            self.diffractometerdialog.show()

    def _onToggleExcludeImage(self, exclude):
        currentimgno = self.scanSelector.slider.value()
        data = self.excludedImagesDialog.getData()
        imgno_in_excudearray = currentimgno in data
        
        if imgno_in_excudearray and exclude:
            return
        if not imgno_in_excudearray and exclude:
            data = np.hstack([data, currentimgno])
            self.excludedImagesDialog.updateArrayData(data)
        else:
            data = data[data != currentimgno]
            self.excludedImagesDialog.updateArrayData(data)        

    def _onSelectCPUcount(self):
        maxavail = os.cpu_count() if os.cpu_count() is not None else 1 
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            maxavail = int(os.environ['SLURM_CPUS_ON_NODE'])
        
        cpus, success = qt.QInputDialog.getInt(self,"CPU count",
                               "CPU count (detected: %s)" % maxavail,
                               self.numberthreads,1)
        if success:
            self.numberthreads = cpus
        
        
    def calcBraggRefl(self):
        if self.fscan is not None and self.reflectionSel.showBraggReflections:
            if self.fscan.axisname != 'th':
                raise NotImplementedError("Calculation of available Bragg reflections is not implemented for %s - scans" % self.fscan.axisname)
            try:
                xtal = self.ubcalc.crystal
                ommin = np.deg2rad(np.amin(self.fscan.omega))
                ommax = np.deg2rad(np.amax(self.fscan.omega))
                dc = self.ubcalc.detectorCal
                mu = self.ubcalc.mu
                ub = self.ubcalc.ubCal
                chi = self.ubcalc.chi
                phi = self.ubcalc.phi
                xtal.setEnergy(ub.getEnergy()*1e3)
                hkls, yx, angles = rn.thscanBragg(xtal,ub,mu,dc,(ommin,ommax), chi=chi, phi=phi)
                self.reflectionSel.setBraggReflections(hkls, yx, angles)
            except Exception:
                qutils.warning_detailed_message(self, "Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections", traceback.format_exc())
                #qt.QMessageBox.critical(self,"Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections:\n%s" % traceback.format_exc())
        

    def saveBraggRefl(self):
        try:
            hkls, yx, angles = self.reflectionSel.getBraggReflections()
        except ValueError:
            if self.fscan is not None:
                try:
                    xtal = self.ubcalc.crystal
                    ommin = np.deg2rad(np.amin(self.fscan.omega))
                    ommax = np.deg2rad(np.amax(self.fscan.omega))
                    dc = self.ubcalc.detectorCal
                    mu = self.ubcalc.mu
                    ub = self.ubcalc.ubCal
                    chi = self.ubcalc.chi
                    phi = self.ubcalc.phi
                    xtal.setEnergy(ub.getEnergy()*1e3)
                    hkls, yx, angles = rn.thscanBragg(xtal,ub,mu,dc,(ommin,ommax), chi=chi, phi=phi)
                    
                    
                    #self.reflectionSel.setBraggReflections(hkls, yx, angles)
                except Exception:
                    qutils.warning_detailed_message(self, "Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections", traceback.format_exc())
                    #qt.QMessageBox.critical(self,"Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections:\n%s" % traceback.format_exc())
                    return
            else:
                qt.QMessageBox.critical(self,"Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections:\nNo scan loaded.")
                return

        hkm = np.concatenate((hkls, yx[:,::-1], np.rad2deg(angles)), axis=1)

        sio = StringIO()
        np.savetxt(sio,hkm,fmt="%.3f", delimiter='\t',header="H K L x y alpha delta gamma omega chi phi")

        #Question dialog for saving the possible CTR locations     
        msgbox = qt.QMessageBox(qt.QMessageBox.Question,'Saving Bragg reflection ...', 
                                'Found possible Bragg reflections. Do you want to save the following positions?',
                                qt.QMessageBox.Yes | qt.QMessageBox.No, self)
        
        msgbox.setDetailedText(sio.getvalue())
        
        clickedbutton = msgbox.exec()
        #Question dialog for saving the possible CTR locations        
        #clickedbutton=qt.QMessageBox.question(self, 'Saving CTR locations...', 'Do you want to save the following positions: \n' + hkstring +"?");
        
        if clickedbutton==qt.QMessageBox.Yes:
            #File saving
            fileTypeDict = {'dat Files (*.dat)': '.dat', 'txt Files (*.txt)': '.txt', 'All files (*)': '', }
            fileTypeFilter = ""
            for f in fileTypeDict:
                fileTypeFilter += f + ";;"
                
            filename, filetype = qt.QFileDialog.getSaveFileName(self,"Save reflections",
                                                      self.filedialogdir,
                                                      fileTypeFilter[:-2])
            if filename == '':
                return
            
            self.filedialogdir = os.path.splitext(filename)[0]
            filename += fileTypeDict[filetype]
            np.savetxt(filename,hkm,fmt="%.3f",header="H K L x y alpha delta gamma omega chi phi")

    def calculateAvailableCTR(self):
        if self.fscan is None:
            raise Exception("No scan selected!")
        if self.fscan.axisname != 'th':
            raise NotImplementedError("Calculation of available CTRs is not implemented for %s - scans" % self.fscan.axisname)
        xtal = self.ubcalc.crystal
        ommin = np.deg2rad(np.amin(self.fscan.omega))
        ommax = np.deg2rad(np.amax(self.fscan.omega))
        dc = self.ubcalc.detectorCal
        mu = self.ubcalc.mu
        chi = self.ubcalc.chi
        phi = self.ubcalc.phi
        ub = self.ubcalc.ubCal
        xtal.setEnergy(ub.getEnergy()*1e3)
        hk, xmirror = rn.thscanCTRs(xtal,ub,mu,dc,(ommin,ommax), chi=chi, phi=phi)
        xmirror = np.array(xmirror).astype(np.float64)
        #making the hk list of arrays into a reasonable string
        hkm = np.concatenate((np.array(hk), xmirror.reshape((1,xmirror.size)).T), axis=1)
        return hkm
        
    def _onCalcAvailableCTR(self):
        try: 
            hkm = self.calculateAvailableCTR()
        except Exception:
            qutils.warning_detailed_message(self, "Cannot calculate CTR locatons", "Cannot calculate CTR locatons", traceback.format_exc())
            #qt.QMessageBox.critical(self,"Cannot calculate CTR locatons", "Cannot calculate CTR locatons:\n%s" % traceback.format_exc())
            return
        sio = StringIO()
        np.savetxt(sio,hkm,fmt="%.3f", delimiter='\t',header="H K detectorRight")

        #Question dialog for saving the possible CTR locations     
        msgbox = qt.QMessageBox(qt.QMessageBox.Question,'Saving CTR locations...', 
                                'Found CTRs. Do you want to save the following positions?',
                                qt.QMessageBox.Yes | qt.QMessageBox.No, self)
        msgbox.setDetailedText(sio.getvalue())
        
        clickedbutton = msgbox.exec()
        #Question dialog for saving the possible CTR locations        
        #clickedbutton=qt.QMessageBox.question(self, 'Saving CTR locations...', 'Do you want to save the following positions: \n' + hkstring +"?");
        
        if clickedbutton==qt.QMessageBox.Yes:
            #File saving
            fileTypeDict = {'dat Files (*.dat)': '.dat', 'txt Files (*.txt)': '.txt', 'All files (*)': '', }
            fileTypeFilter = ""
            for f in fileTypeDict:
                fileTypeFilter += f + ";;"
                
            filename, filetype = qt.QFileDialog.getSaveFileName(self,"Save reflections",
                                                      self.filedialogdir,
                                                      fileTypeFilter[:-2])
            if filename == '':
                return
            
            self.filedialogdir = os.path.splitext(filename)[0]
            filename += fileTypeDict[filetype]
            np.savetxt(filename,hkm,fmt="%.3f",header="H K mirror")
            
        
    def getReflections(self):
        hkls = []
        angles = []
        for refl in self.reflectionSel.reflections:
            #print(refl.xy)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([refl.xy[1]]),np.array([refl.xy[0]]),self.ubcalc.mu)
            delta = float(delta); gamma = float(gamma)
            try:
                pos = np.array([self.ubcalc.mu,delta,gamma,self.imageNoToOmega(refl.imageno),self.ubcalc.chi,self.ubcalc.phi])
            except:
                from IPython import embed; embed()
                raise
            #print(pos)
            hkls.append(refl.hkl)
            angles.append(pos)
        return np.array(hkls), np.array(angles)


    def _onPlotMachineParams(self, enable=None):
        #[cp,azimxy,polax] = paramslist
        if enable is None:
            enable = self.showMachineParamsAct.isChecked()
        if enable:
            fit2DCal = self.ubcalc.detectorCal.getFit2D()
            cp = fit2DCal['centerX'], fit2DCal['centerY']
            gam_p,_ = self.ubcalc.detectorCal.rangegamdel_p
            azimy,azimx = self.ubcalc.detectorCal.pixelsPrimeBeam(gam_p[1]/5, 0 )[0]
            self.centralPlot.addMarker(cp[0],cp[1],legend="CentralPixel",text="CP",color='yellow',symbol='+')
            self.centralPlot.addMarker(azimx,azimy,legend="azimuth",text="Azim",color='yellow',symbol='+')
        else:
            self.centralPlot.removeMarker("CentralPixel")
            self.centralPlot.removeMarker("azimuth")
            
    def searchPixelCoordHKL(self, hkl):
        refldict = self.ubcalc.calcReflection(hkl)
        axisname = self.fscan.axisname
        dc = self.ubcalc.detectorCal
        
        if self.fscan.axisname == 'mu':
            angle_idx = 0
            sign = 1.
        elif self.fscan.axisname == 'th':
            angle_idx = 3
            sign = -1.
        else:
            qt.QMessageBox.warning(self,"Cannot calculate reflection","Cannot calculate reflection.\n%s is no supported scan axis." % self.fscan.axisname)
            return 
        try:
            imageno1 = self.axisToImageNo(np.rad2deg(refldict['angles_1'][angle_idx]) * sign)
            refldict['imageno_1'] = imageno1
        except:
            imageno1 = None
        xy = refldict['xy_1']
        onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and \
                     (xy[1] >= 0 and xy[1] < dc.detector.shape[0])
        if onDetector:
            refldict['selectable_1'] = True
        else:
            refldict['selectable_1'] = False
            
        try:
            imageno2 = self.axisToImageNo(np.rad2deg(refldict['angles_2'][angle_idx]) * sign)
            refldict['imageno_2'] = imageno2
        except:
            imageno2 = None
        xy = refldict['xy_2']
        onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and \
                     (xy[1] >= 0 and xy[1] < dc.detector.shape[0])
        if onDetector:
            refldict['selectable_2'] = True
        else:
            refldict['selectable_2'] = False
        return refldict
            
    def onSearchHKLforStaticROI(self, hkl):
        try:
            refldict = self.searchPixelCoordHKL(hkl)
        except Exception as e:
            qutils.warning_detailed_message(self, "Cannot calculate location of reflection", "Cannot calculate position of reflection:\n%s" % e, traceback.format_exc())
            return
        refl_dialog = QReflectionAnglesDialog(refldict,"Select reflection location", self)
        if qt.QDialog.Accepted == refl_dialog.exec():
            for i, cb in enumerate(refl_dialog.checkboxes,1):
                if cb.isChecked():
                    xy = refldict['xy_%s' % i]
                    self.scanSelector.set_xy_static_loc(xy[0], xy[1])
                    return      
    
    def _onStaticROIedited(self):
        xy = self.roiS1.getCenter()
        hsize, vsize = np.round(self.roiS1.getSize())
        self.scanSelector.hsize.blockSignals(True)
        self.scanSelector.vsize.blockSignals(True)
        self.scanSelector.hsize.setValue(hsize)
        self.scanSelector.vsize.setValue(vsize)
        self.scanSelector.hsize.blockSignals(False)
        self.scanSelector.vsize.blockSignals(False)
        self.scanSelector.set_xy_static_loc(xy[0], xy[1])
        

    def _onNewReflection(self,refldict):
        axisname = self.fscan.axisname
        dc = self.ubcalc.detectorCal
        
        if self.fscan.axisname == 'mu':
            angle_idx = 0
            sign = 1.
        elif self.fscan.axisname == 'th':
            angle_idx = 3
            sign = -1.
        else:
            qt.QMessageBox.warning(self,"Cannot calculate reflection","Cannot calculate reflection.\n%s is no supported scan axis." % self.fscan.axisname)
            return 
        try:
            imageno1 = self.axisToImageNo(np.rad2deg(refldict['angles_1'][angle_idx]) * sign)
            refldict['imageno_1'] = imageno1
            xy = refldict['xy_1']
            onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and \
                         (xy[1] >= 0 and xy[1] < dc.detector.shape[0])
            if onDetector:
                refldict['selectable_1'] = True
            else:
                refldict['selectable_1'] = False
        except:
            imageno1 = None
            refldict['selectable_1'] = False
        try:
            imageno2 = self.axisToImageNo(np.rad2deg(refldict['angles_2'][angle_idx]) * sign)
            refldict['imageno_2'] = imageno2
            xy = refldict['xy_2']
            onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and \
                         (xy[1] >= 0 and xy[1] < dc.detector.shape[0])
            if onDetector:
                refldict['selectable_2'] = True
            else:
                refldict['selectable_2'] = False
        except:
            imageno2 = None
            refldict['selectable_2'] = False
            
        refl_dialog = QReflectionAnglesDialog(refldict,"Select reflections to add into list of reference reflections", self)
        if qt.QDialog.Accepted == refl_dialog.exec():
            for i, cb in enumerate(refl_dialog.checkboxes,1):
                if cb.isChecked():
                    xy = refldict['xy_%s' % i]
                    eventdict = {'x' : xy[0], 'y': xy[1]}
                    self.reflectionSel.addReflection(eventdict,refldict['imageno_%s' % i],refldict['hkl'])
            
    def newXyHKLConverter(self):
        def xyToHKL(x,y):
            #print("xytoHKL:")
            #print("x,y = %s, %s" % (x,y))
            if self.fscan is None:
                return np.array([np.nan,np.nan,np.nan, np.nan, np.nan])
            mu, om = self.getMuOm(self.imageno)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([y]),np.array([x]), mu)
            #print(self.ubcalc.detectorCal)
            #print(x,y)
            #print(self.ubcalc.detectorCal.tth(np.array([y]),np.array([x])))
            pos = [mu,delta[0],gamma[0],om,self.ubcalc.chi,self.ubcalc.phi]
            pos = HKLVlieg.crystalAngles(pos,self.ubcalc.n)
            hkl = np.concatenate((np.array(self.ubcalc.angles.anglesToHkl(*pos)),np.rad2deg([delta[0],gamma[0]])))
            return hkl
        return xyToHKL
        
    def getMuOm(self, imageno=None):
        if imageno is not None:
            if self.fscan.axisname == 'th':
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.imageNoToAxis(imageno))
            elif self.fscan.axisname == 'mu':
                mu = np.deg2rad(self.imageNoToAxis(imageno))
                om = -1 * np.deg2rad(self.fscan.th)
                if len(np.asarray(om).shape) > 0:
                    om = om[0]
            else:
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.fscan.th)
                if len(np.asarray(om).shape) > 0:
                    om = om[0]
            return mu, om
        else:
            if self.fscan.axisname == 'th':
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.fscan.axis)
            elif self.fscan.axisname == 'mu':
                mu = np.deg2rad(self.fscan.axis)
                om = -1 * np.deg2rad(self.fscan.th)
            else:
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.fscan.th)
            return mu, om
        
    def omegaToImageNo(self,omega):
        if self.fscan is not None:
            omrad = np.deg2rad(self.fscan.omega)
            ommax = np.amax(omrad)
            ommin = np.amin(omrad)
            #print(ommin,omega,ommax)
            if omega < ommin or omega > ommax:
                omdeg =  np.rad2deg([ommin,omega,ommax])
                raise Exception("omega not in range: %s < %s < %s" % tuple(omdeg))
            return np.argmin(np.abs(omrad -omega))
        else:
            raise Exception("No Scan selected")
            
    def imageNoToOmega(self,imageno):
        if self.fscan is not None:
            return np.deg2rad(self.fscan.omega[imageno])
        else:
            return 0.

    def imageNoToAxis(self,imageno):
        if self.fscan is not None:
            return self.fscan.axis[imageno]
        else:
            return 0.

    def axisToImageNo(self,axisval):
        if self.fscan is not None:
            #axis = np.deg2rad(self.fscan.axis)
            axismax = np.amax(self.fscan.axis)
            axismin = np.amin(self.fscan.axis)
            #print(ommin,omega,ommax)
            if axisval < axismin or axisval > axismax:
                axisrange = [axismin,axisval,axismax]
                raise Exception("Value of scan axis \"%s\" not in range: %s < %s < %s" % tuple([self.fscan.axisname]+axisrange))
            return np.argmin(np.abs(self.fscan.axis - axisval))
        else:
            raise Exception("No Scan selected")
            
    def _onCreateScan(self):
        try:
            mu, om = self.getMuOm(self.imageno)
        except:
            mu = self.ubcalc.mu
            om = 0.
        th = om*-1.
        muTh = np.rad2deg([mu,th]) #defaults if fixed 
        diag = QScanCreator(muTh)
        if diag.exec() == qt.QDialog.Accepted:
            shape = self.ubcalc.detectorCal.detector.shape
            try:
                axis = diag.scanaxis.currentText()
                if axis == 'theta':
                    axis = 'th'
                elif axis == 'mu':
                    pass
                fscan = SimulationScan(shape, diag.omstart.value(),
                                        diag.omend.value(),
                                        diag.no.value(),
                                        axis, diag.fixedAngle.value())
                self._onScanChanged(fscan)
            except MemoryError:
                qutils.warning_detailed_message(self, "Can not create simulation scan","Can not create simualtion scan. Memory is insufficient for the scan size. See details for further information.", traceback.format_exc())
        
    def _onLoadScanFromImages(self):
        # generates a scan from a selected folder containing raw detector images
        
        # generate file source selection GUI

        # create filter of scan image formats (following code is copied from silx view)
        extensions = {}
        for description, ext in silx.io.supported_extensions().items():
            extensions[description] = " ".join(sorted(list(ext)))

        extensions["NeXus layout from EDF files"] = "*.edf"
        extensions["NeXus layout from TIFF image files"] = "*.tif *.tiff"
        extensions["NeXus layout from CBF files"] = "*.cbf"
        extensions["NeXus layout from MarCCD image files"] = "*.mccd"

        all_supported_extensions = set()
        for name, exts in extensions.items():
            exts = exts.split(" ")
            all_supported_extensions.update(exts)
        all_supported_extensions = sorted(list(all_supported_extensions))

        filters = []
        filters.append("All supported files (%s)" % " ".join(all_supported_extensions))
        for name, extension in extensions.items():
            filters.append("%s (%s)" % (name, extension))
        filters.append("All files (*)")

        fileTypeFilter = ""
        for f in filters:
            fileTypeFilter += f + ";;"

        # call dialog
        filename,_ = qt.QFileDialog.getOpenFileName(self,"Open data source",'',fileTypeFilter[:-2])

        # Qt dialog returns '' if cancelled
        if filename == '':
            qt.QMessageBox.warning(self,"Error - Open data source","No data source selected")
            return

        # search files using ImportImagesScan backend
        importedscan = universalScanLoader.ImportImagesScan(filename)

        if importedscan.inpath == None:
            qt.QMessageBox.critical(self,
                                    "Images could not be imported",
                                    "The selected data source is not suitable\n"\
                                    "It is necessary to select a file containing raw detector image(s)!")
            return
        
        if importedscan.shape != self.ubcalc.detectorCal.detector.shape:
            qt.QMessageBox.critical(self,
                                    "Detector data mismatch",
                                    "The selected image data shape does not match to the detector data shape:\n"\
                                    "Detector Size %sx%s\n"\
                                    "Data size %sx%s\n"\
                                    "Please first adjust the detector configuration to load this data" % 
                                    (*self.ubcalc.detectorCal.detector.shape, *importedscan.shape))
            return
        
        [imagePrefix, found_scanfiles] = importedscan.inpath

        # generate dialog with list of files and frames
        nrofFilesfound = len(found_scanfiles)
        messageStr = 'Found ' + str(nrofFilesfound) + ' files in selected directory'

        if importedscan.FramesPerFile > 1:

            if nrofFilesfound == 0:
                messageStr = 'No images found!!!'
                fullStr = messageStr
            elif 0 < nrofFilesfound < 4:
                messageStr += ':\n'
                for i in range(nrofFilesfound-1):
                    messageStr += imagePrefix + found_scanfiles[i] + ': ' + str(importedscan.FramesPerFile) + ' frames' ### mark expected nr when file not actually loaded
                    if i > 0:
                        messageStr += ' (expected)'
                    messageStr += '\n'
                messageStr += imagePrefix + found_scanfiles[nrofFilesfound-1] + ': ' + str(importedscan.FramesLastFile) + ' frames\n' + str(importedscan.nopoints) + ' frames in total.'
                fullStr = messageStr
            else:
                messageStr += ':\n'
                for i in range(0,3):
                    messageStr += imagePrefix + found_scanfiles[i] + ': ' + str(importedscan.FramesPerFile) + ' frames'
                    if i > 0:
                        messageStr += ' (expected)'
                    messageStr += '\n'
                fullStr = messageStr
                for i in range(3,nrofFilesfound-1):
                    fullStr += imagePrefix + found_scanfiles[i] + ': ' + str(importedscan.FramesPerFile) + ' frames (expected) \n'

                messageStr += '...' + '\n' + imagePrefix + found_scanfiles[nrofFilesfound-1] + ': ' + str(importedscan.FramesLastFile) + ' frames\n' + str(importedscan.nopoints) + ' frames in total.'
                fullStr += imagePrefix + found_scanfiles[nrofFilesfound-1] + ': ' + str(importedscan.FramesLastFile) + ' frames\n' + str(importedscan.nopoints) + ' frames in total.'


        else:

            if nrofFilesfound == 0:
                messageStr = 'No images found!!!'
                fullStr = messageStr
            elif 0 < nrofFilesfound < 4:
                messageStr += ':\n'
                for i in range(nrofFilesfound-1):
                    messageStr += imagePrefix + found_scanfiles[i] + '\n'
                messageStr += imagePrefix + found_scanfiles[nrofFilesfound-1] 
                fullStr = messageStr
            else:
                messageStr += ':\n'
                for i in range(0,3):
                    messageStr += imagePrefix + found_scanfiles[i] + '\n'
                fullStr = messageStr
                for i in range(3,nrofFilesfound):
                    fullStr += imagePrefix + found_scanfiles[i] + '\n'
                messageStr += '...' + '\n' + imagePrefix + found_scanfiles[nrofFilesfound-1]
                fullStr += '\n' + str(importedscan.nopoints) + ' frames in total.'
            
        

        msg0 = qt.QMessageBox(self)
        msg0.setWindowTitle("Manual scan import")
        msg0.setText(messageStr)
        msg0.setDetailedText(fullStr)
        msg0.exec()

        # angle conversions
        try:
            mu, om = self.getMuOm(self.imageno)
        except:
            mu = self.ubcalc.mu
            om = 0.
        th = om*-1.
        muTh = np.rad2deg([mu,th]) #defaults if fixed 
        
        # open scan creator GUI to let the user insert missing scan angles
        diag = QImportScanCreator(muTh)        
        # detector pixel nr and frame nr is adapted from opened image file
        diag.no.setValue(importedscan.nopoints)
        
        if diag.exec() == qt.QDialog.Accepted:
            try:
                axis = diag.scanaxis.currentText()
                if axis == 'theta':
                    axis = 'th'
                elif axis == 'mu':
                    pass
                # pass inserted angles to scan object
                importedscan.set_axis(diag.omstart.value(),diag.omend.value(),axis,diag.fixedAngle.value())
                self._onScanChanged(importedscan)
            except MemoryError:
                qutils.warning_detailed_message(self, "Can not create scan","Can not create scan. Memory is insufficient for the scan size. See details for further information.", traceback.format_exc())
        
    def _onScanChanged(self,sel_list):
        self.resetZoom = True
        #print(sel_list)
        self.activescanname = "scan"
        if isinstance(sel_list,list): 
            self.sel_list = sel_list
            if len(sel_list):
                self.specfile = sel_list[0]['SourceName']
                try:
                    self.scanno = int(float(sel_list[0]['Key']))-1
                    self.fscan = Fastscan(self.specfile,self.scanno)
                    self.imageno = 0
                except Exception:
                    self.scanno = 0
                    self.fscan = FioFastsweep(self.specfile)
                    self.imageno = 0
                self.reflectionSel.setImage(self.imageno)
                if self.imagepath != '':
                    self.fscan.set_image_folder(self.imagepath)
                    self.plotImage()
                    self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
                    
            else:
                self.scanSelector.setRange(0,0)
                self.imageno = 0
                self.reflectionSel.setImage(self.imageno)
                #print(self.centralPlot._callback)

        elif isinstance(sel_list,universalScanLoader.ImportImagesScan):
            self.scanno = 1
            self.fscan = sel_list
            self.imageno = 0
            self.plotImage()
            self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
            self.activescanname = "%s-rawImport %s-%s" % (self.fscan.axisname, np.amin(self.fscan.axis),np.amax(self.fscan.axis))

            self.images_loaded = False
            if self.fscan is not None and self.autoLoadAct.isChecked():
                self.loadAll()
                self.scanSelector.showMaxAct.setChecked(False)
                self.scanSelector.showMaxAct.setChecked(True)

        elif isinstance(sel_list,SimulationScan):
            self.scanno = 1
            self.fscan = sel_list
            self.imageno = 0
            self.plotImage()
            self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
            self.activescanname = "%s-sim %s-%s" % (self.fscan.axisname, np.amin(self.fscan.axis),np.amax(self.fscan.axis))
        else:
            if 'name' in sel_list:
                self.activescanname = sel_list['name']
            else:
                self.activescanname = "scan"
            self.hdffile = sel_list['file']
            #self.scanname = sel_list['name'].strip("/")
            try:
                msg = qt.QMessageBox(self)
                msg.setWindowTitle("Loading Scan")
                msg.setText("Loading Scan. This might take a while...")
                msg.setStandardButtons(qt.QMessageBox.Cancel)
                msg.setModal(True)
                msg.show()
                if 'beamtime' in sel_list:
                    self.fscan = backends.openScan(sel_list['beamtime'], sel_list)
                else:
                    self.fscan = backends.openScan(self.scanSelector.btid.currentText(), sel_list)
                    
                self.plotImage()
                self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
                msg.hide()
                self.images_loaded = False
                if self.fscan is not None and self.autoLoadAct.isChecked():
                    self.loadAll()
                    self.scanSelector.showMaxAct.setChecked(False)
                    self.scanSelector.showMaxAct.setChecked(True)
            except Exception:
                msg.hide()
                qutils.warning_detailed_message(self, "Cannot open scan", "Cannot open scan" , traceback.format_exc())
                #qt.QMessageBox.critical(self,"Cannot open scan", "Cannot open scan:\n%s" % traceback.format_exc())
        if hasattr(self.fscan, 'name'):
            self.activescanname = self.fscan.name
                
            
            
    def _onImagePathChanged(self,path):
        #print("newpath %s" % path)
        self.imagepath = path
        if self.fscan is not None:
            self.fscan.set_image_folder(self.imagepath)
            self.plotImage()
            self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
            #self.scanSelector.slider.setMinimum(0)
            #self.scanSelector.slider.setMaximum(self.fscan.nopoints-1)
        else:
            self.scanSelector.setRange(0,0)
            self.imageno = 0
            self.reflectionSel.setImage(self.imageno)
            #print(self.centralPlot._callback)
        
    def _onChangeImage(self,imageno):
        if self.fscan is not None:
            self.scanSelector.slider.setValue(imageno)
            self.plotImage(self.scanSelector.slider.value())
        
    def _onSliderValueChanged(self,value):
        if self.fscan is not None: 
            self.plotImage(value)
        #print(self.centralPlot._callback)
            
        
    def _onLoadAll(self):
        self.images_loaded = False
        if self.fscan is not None:
            self.loadAll()
            self.scanSelector.showMaxAct.setChecked(False)
            self.scanSelector.showMaxAct.setChecked(True)
            
    def loadAll(self):
        try:
            image = self.fscan.get_raw_img(0)
        except Exception as e:
            print("no images found! %s" % e)
            return
        self.allimgsum = np.zeros_like(image.img)
        self.allimgmax = np.zeros_like(image.img)
        progress = qt.QProgressDialog("Reading images","abort",0,len(self.fscan),self)
        progress.setWindowModality(qt.Qt.WindowModal)
        #tasks = queue.Queue()
        #[tasks.put(i) for i in range(len(self.fscan))]
        lock = threading.Lock()
        self.images_loaded = True
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.numberthreads) as executor: # speedup only for the file reads 
            futures = {}
            def readfile_max(imgno):
                if imgno in self.excludedImagesDialog.getData(): # skip if excluded
                    return imgno
                image = self.fscan.get_raw_img(imgno) # here speedup during file read
                with lock:
                    self.allimgsum += image.img
                    self.allimgmax = np.maximum(self.allimgmax,image.img)
                return imgno
            for i in range(len(self.fscan)):
                futures[executor.submit(readfile_max, i)] = i
            
            for f in concurrent.futures.as_completed(futures):
                try:
                    imgno = f.result()
                    progress.setValue(imgno)
                except concurrent.futures.CancelledError:
                    pass
                except Exception as e:
                    print("Cannot read image:\n%s" % traceback.format_exc())

                if progress.wasCanceled():
                    [f.cancel() for f in futures]
                    self.images_loaded = False
                    break

        progress.setValue(len(self.fscan))
        
    def _onMaxToggled(self,value):
        if self.scanSelector.showSumAct.isChecked():
            self.scanSelector.showSumAct.setChecked(False)
        if value:
            if not self.images_loaded and self.fscan is not None:
                btn = qt.QMessageBox.question(self,"Incomplete sum / max image", "Sum/Max image was not loaded completely. Displayed maximum image will be incomplete! Do you want to load all images?",qt.QMessageBox.Yes | qt.QMessageBox.No | qt.QMessageBox.Cancel)
                if btn == qt.QMessageBox.Yes: 
                    self.loadAll()
                elif btn == qt.QMessageBox.Cancel:
                    self.scanSelector.showMaxAct.setChecked(False)
                    return
            if self.allimgmax is not None:
                self.currentAddImageLabel = self.centralPlot.addImage(self.allimgmax,legend="special",
                                                               replace=False,resetzoom=False,copy=True,z=1)
                self.centralPlot.setActiveImage(self.currentAddImageLabel)
                self.scanSelector.alphaslider.setLegend(self.currentAddImageLabel)
            else:
                self.scanSelector.showMaxAct.setChecked(False)
        else:
            if self.currentAddImageLabel is not None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
                self.centralPlot.removeImage(self.currentAddImageLabel)
                self.currentAddImageLabel = None

        
    def _onSumToggled(self,value):
        if self.scanSelector.showMaxAct.isChecked():
            self.scanSelector.showMaxAct.setChecked(False)
        if value:
            if not self.images_loaded and self.fscan is not None:
                btn = qt.QMessageBox.question(self,"Incomplete sum / max image", "Sum/Max image was not loaded completely. Displayed sum image will be incomplete! Do you want to load all images?",qt.QMessageBox.Yes | qt.QMessageBox.No | qt.QMessageBox.Cancel)
                if btn == qt.QMessageBox.Yes: 
                    self.loadAll()
                elif btn == qt.QMessageBox.Cancel:
                    self.scanSelector.showSumAct.setChecked(False)
                    return
            if self.allimgsum is not None:
                self.currentAddImageLabel = self.centralPlot.addImage(self.allimgsum,legend="special",
                                                               replace=False,resetzoom=False,copy=True,z=1)
                self.centralPlot.setActiveImage(self.currentAddImageLabel)
                self.scanSelector.alphaslider.setLegend(self.currentAddImageLabel)
            else:
                self.scanSelector.showSumAct.setChecked(False)
        else:
            if self.currentAddImageLabel is not None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
                self.centralPlot.removeImage(self.currentAddImageLabel)
                self.currentAddImageLabel = None

        
    
        
    def plotImage(self,key=0):
        try:
            image = self.fscan.get_raw_img(key)
            #if self.currentImageLabel is not None:
            #    self.centralPlot.removeImage(self.currentImageLabel)

            self.currentImageLabel = self.centralPlot.addImage(image.img,legend="scan_image",
                                                               replace=False,resetzoom=self.resetZoom,copy=True)
            if self.currentAddImageLabel is None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
            self.resetZoom = False
            self.imageno = key
            self.reflectionSel.setImage(self.imageno)
            self.updateROI(image_changed=True)
            self.updateReflections()
            
            mu, om = self.getMuOm(self.imageno)
            self.ubcalc.uedit.setAngles(mu, self.ubcalc.chi, self.ubcalc.phi, om)
            
            self.scanSelector.excludeImageAct.blockSignals(True)
            self.scanSelector.excludeImageAct.setChecked(key in self.excludedImagesDialog.getData())
            self.scanSelector.excludeImageAct.blockSignals(False)
                        
        except Exception as e:
            print(traceback.format_exc())
            print("no image %s" % e)


    def updateReflections(self):
        if not self.reflectionsVisible:
            self.centralPlot.removeCurve('all_image_reflections')
            return
        H_0, H_1 = self.reflectionsToDisplay
        #H_0 = np.array([[1,0,0], [1,1,0]])
        #H_1 = np.array([[0,0,1], [0,0,1]])
        
        hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(self.imageno, H_0, H_1, intersect=True)
        
        mask1 = hkl_del_gam_1[:,-1].nonzero()
        mask2 = hkl_del_gam_2[:,-1].nonzero()
        
        masked_hkl_del_gam = np.vstack((hkl_del_gam_1[mask1],hkl_del_gam_2[mask2]))
        
        self.centralPlot.addCurve(masked_hkl_del_gam[:,-3],masked_hkl_del_gam[:,-2],legend='all_image_reflections',
                                  linestyle=' ', symbol='.', color='y',resetzoom=False)
        
    
        
    def updateROI(self, **kwargs):
        if not self.roivisible:
            #for roi in self.rois:
            self.roiS1.setVisible(False)
            self.roiS2.setVisible(False)
            if self.rocking_rois:
                for roi in self.rocking_rois:
                    roi.setVisible(False)
            self.roiManager._roisUpdated()
            return
            #self.centralPlot.removeMarker('main_croi_loc')
        
        current_mode = self.scanSelector.scanstab.currentIndex()
        if (current_mode == 0 or current_mode == 1):
            if self.rocking_rois:
                for roi in self.rocking_rois:
                    roi.setVisible(False)
            try:
                hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(self.imageno)
            except:
                #for roi in self.rois:
                #    roi.setVisible(False)
                self.roiS1.setVisible(False)
                self.roiS2.setVisible(False)
                #self.centralPlot.removeMarker('main_croi_loc')
                return


            if hkl_del_gam_1[0,-1]:
                if self.roivisible:
                    self.plotROI(hkl_del_gam_1[0,6:8], self.roiS1)
                for i, spinbox in enumerate(self.scanSelector.hkl_static):
                    spinbox.setValue(hkl_del_gam_1[0,i])
            else:
                self.roiS1.setVisible(False)

            if hkl_del_gam_2[0,-1]:
                if self.roivisible:
                    self.plotROI(hkl_del_gam_2[0,6:8], self.roiS2)
                for i, spinbox in enumerate(self.scanSelector.hkl_static):
                    spinbox.setValue(hkl_del_gam_1[0,i])
            else:
                self.roiS2.setVisible(False)
                
            if current_mode == 1:
                self.roiS1.setEditable(True)
            else:
                self.roiS1.setEditable(False)
                
        elif (current_mode == 2 and not kwargs.get('image_changed', False)):
            self.roiS1.setVisible(False)
            self.roiS2.setVisible(False)
            try:
                refldict = self.get_rocking_coordinates()
            except:
                #print(traceback.format_exc())
                if self.rocking_rois:
                    for roi in self.rocking_rois:
                        roi.setVisible(False)
                return
            roi_keys = self.intbkgkeys_rocking(refldict)
            self.scanSelector.autoSize_label.setText("%s x %s" % (roi_keys['hsize'], roi_keys['vsize']))
            
            number_rois = len(roi_keys['center'])
            divider = 1
            if number_rois > MAX_ROIS_DISPLAY:
                divider = np.ceil(number_rois / MAX_ROIS_DISPLAY)
            no_rois_to_display = int(np.floor(number_rois / divider))
            
            # lazy create ROIs 
            if len(self.rocking_rois) < no_rois_to_display:
                for i in range(no_rois_to_display - len(self.rocking_rois)):
                    roi = RectangleBgROI()
                    roi.setLineWidth(1)
                    roi.setLineStyle('-')
                    roi.setColor('red')
                    roi.setBgStyle('pink', '-', 1.)
                    roi.setVisible(False)
                    roi.setGeometry(origin=(0, 0), size=(0, 0))
                    self.rocking_rois.append(roi)
                    self.roiManager.addRoi(roi,useManagerColor=False)
            
            for roino, i in enumerate(np.arange(no_rois_to_display)*divider):
                ckey = roi_keys['center'][int(i)]
                leftkey = roi_keys['left'][int(i)]
                rightkey = roi_keys['right'][int(i)]
                topkey = roi_keys['top'][int(i)]
                bottomkey = roi_keys['bottom'][int(i)]

                origin =(ckey[0].start, ckey[1].start)
                size = (ckey[0].stop - ckey[0].start, ckey[1].stop - ckey[1].start)
                #leftkey, rightkey, topkey, bottomkey
                left = leftkey[0].stop - leftkey[0].start
                right = rightkey[0].stop - rightkey[0].start
                top = topkey[1].stop - topkey[1].start
                bottom = bottomkey[1].stop - bottomkey[1].start
                self.rocking_rois[roino].setGeometry(origin=origin, size=size, left=left, right=right, top=top, bottom=bottom)
                self.rocking_rois[roino].setVisible(True)
            for roi in self.rocking_rois[no_rois_to_display:]:
                roi.setVisible(False)
        
        self.roiManager._roisUpdated()
                #self.centralPlot.removeMarker('main_croi_loc')
            

    def getStaticROIparams(self, xy, **kwargs):
        if self.fscan is None:
            raise Exception("No scan loaded!")
        mu, om = self.getMuOm()
        #mu_cryst = HKLVlieg.crystalAngles_singleArray(mu, self.ubcalc.n)
        dc = self.ubcalc.detectorCal
        angles = self.ubcalc.angles
        
        if 'mask' in kwargs:
            mask = kwargs['mask']
            xy = xy[mask]

        if len(np.asarray(om).shape) == 0:
            om = np.full(len(self.fscan),om)
            
        hkl_del_gam = np.empty((xy.shape[0],len(self.fscan), 6), dtype=np.float64)
        for i, xy_i in enumerate(xy):
            x = np.full(len(self.fscan),xy_i[0])
            y = np.full(len(self.fscan),xy_i[1])
            gamma, delta, alpha = self.ubcalc.detectorCal.crystalAnglesPoint(y, x, mu, self.ubcalc.n)

            if len(np.asarray(alpha).shape) == 0:
                alpha = np.full(len(self.fscan),alpha)

            hkl = self.ubcalc.angles.anglesToHkl(alpha, delta, gamma, om, self.ubcalc.chi, self.ubcalc.phi)
            #for i in range(len(self.fscan)):

            hkl_del_gam[i,:,:3] = np.array(hkl).T
            hkl_del_gam[i,:, 3] = delta
            hkl_del_gam[i,:, 4] = gamma
            hkl_del_gam[i,:, 5] = self.fscan.axis
        return hkl_del_gam
        
        

    def getROIloc(self, imageno=None, H_0=None, H_1=None, **kwargs):
        if self.fscan is None:
            raise Exception("No scan loaded!")
        

        mu, om = self.getMuOm(imageno)
        mu_cryst = HKLVlieg.crystalAngles_singleArray(mu, self.ubcalc.n)
        dc = self.ubcalc.detectorCal
        #mu = self.ubcalc.mu
        angles = self.ubcalc.angles

        if (self.scanSelector.scanstab.currentIndex() == 1 and not kwargs.get('intersect', False)):
            if imageno is None:
                hkl_del_gam_1 = np.ones((len(self.fscan),6),dtype=np.float64)
                x = np.full(len(self.fscan),self.scanSelector.xy_static[0].value())
                y = np.full(len(self.fscan),self.scanSelector.xy_static[1].value())
                
                gamma, delta, alpha = self.ubcalc.detectorCal.crystalAnglesPoint(y, x, mu,  self.ubcalc.n)
                s = np.arange(len(self.fscan))
                
                if len(np.asarray(om).shape) == 0:
                    om = np.full(len(self.fscan),om)
                
                if len(np.asarray(alpha).shape) == 0:
                    alpha = np.full(len(self.fscan),alpha)

                yx1 = np.vstack((y,x)).T
                yx2 = np.full_like(yx1, np.inf)
                
                for i in range(len(self.fscan)):
                    
                    pos = [alpha[i],delta[i],gamma[i],
                            om[i],
                            self.ubcalc.chi,
                            self.ubcalc.phi]
                    #pos = np.vstack(pos).T
                    hkl = np.array(self.ubcalc.angles.anglesToHkl(*pos))
                    hkl_del_gam_1[i, :3] = hkl
                hkl_del_gam_1[:, 3] = delta
                hkl_del_gam_1[:, 4] = gamma
                hkl_del_gam_1[:, 5] = self.fscan.axis
                hkl_del_gam_2 = np.full_like(hkl_del_gam_1, -1)
            else:
                hkl_del_gam_1 = np.ones(6,dtype=np.float64)
                x = self.scanSelector.xy_static[0].value()
                y = self.scanSelector.xy_static[1].value()
                yx1 = np.zeros((1,2))
                yx1[0][0] = y
                yx1[0][1] = x
                yx2 = np.full_like(yx1, np.inf)
                
                if len(np.asarray(om).shape) > 0:
                    om = om[imageno]
                if len(np.asarray(mu).shape) > 0:
                    mu = mu[imageno]
                gamma, delta, alpha = self.ubcalc.detectorCal.crystalAnglesPoint(np.array([y]),np.array([x]), mu,  self.ubcalc.n)
                pos = [alpha,delta,gamma,om,self.ubcalc.chi,self.ubcalc.phi]
                hkl_del_gam_1[:3] = np.concatenate(self.ubcalc.angles.anglesToHkl(*pos))
                hkl_del_gam_1[3] = delta
                hkl_del_gam_1[4] = gamma
                hkl_del_gam_1[5] = self.fscan.axis[imageno]
                hkl_del_gam_2 = np.full_like(hkl_del_gam_1, -1)
        
        else:
            if H_0 is None or H_1 is None:
                H_1 = np.array([h.value() for h in self.scanSelector.H_1])
                H_0 = np.array([h.value() for h in self.scanSelector.H_0])
            
            hkl_del_gam_1, hkl_del_gam_2, Qa_1, Qa_2 = angles.anglesIntersectLineEwald(H_0, H_1, mu_cryst, om, self.ubcalc.phi,self.ubcalc.chi, Qalpha=True)
            # H, K, L ,delta_1, gamma_1, HKL_Q1[-1]=s
            
            delta1 = hkl_del_gam_1[...,3]
            delta2 = hkl_del_gam_2[...,3]
            gam1 = hkl_del_gam_1[...,4]
            gam2 = hkl_del_gam_2[...,4]
            
            Qmin, Qmax = dc.Qrange
            Qa_1_n = np.linalg.norm(Qa_1, axis=-1)
            Qa_2_n = np.linalg.norm(Qa_2, axis=-1)
            
            
            mask1 = np.logical_and(Qmin <= Qa_1_n , Qmax >= Qa_1_n)
            mask2 = np.logical_and(Qmin <= Qa_2_n , Qmax >= Qa_2_n)
            
            yx1 = dc.pixelsCrystalAngles(gam1, delta1, mu, self.ubcalc.n)
            yx2 = dc.pixelsCrystalAngles(gam2, delta2, mu, self.ubcalc.n)
            yx1[~mask1] = np.inf
            yx2[~mask2] = np.inf
        
        ymask1 = np.logical_and(yx1[...,0] >= 0, yx1[...,0] < dc.detector.shape[0])
        xmask1 = np.logical_and(yx1[...,1] >= 0, yx1[...,1] < dc.detector.shape[1])
        yxmask1 = np.logical_and(xmask1,ymask1)
    
        ymask2 = np.logical_and(yx2[...,0] >= 0, yx2[...,0] < dc.detector.shape[0])
        xmask2 = np.logical_and(yx2[...,1] >= 0, yx2[...,1] < dc.detector.shape[1])
        yxmask2 = np.logical_and(xmask2,ymask2)
        
        xy1 = yx1[...,::-1]
        xy2 = yx2[...,::-1]
        
        if not kwargs.get('intersect', False):
            xoffset = self.scanSelector.offsetx.value()
            yoffset = self.scanSelector.offsety.value()
            
            if xoffset != 0. or yoffset != 0.:
                warnings.warn("Nonzero pixel offset selected. Experimental feature! Angles and hkl are incorrect!!!")
                xy1[..., 0] += xoffset
                xy2[..., 0] += xoffset
                xy1[..., 1] += yoffset
                xy2[..., 1] += yoffset
        
        return np.concatenate((np.atleast_2d(hkl_del_gam_1), xy1, yxmask1[...,np.newaxis]),axis=-1),\
               np.concatenate((np.atleast_2d(hkl_del_gam_2), xy2, yxmask2[...,np.newaxis]),axis=-1)

    def plotROI(self, loc, roi):

        key = self.intkey(loc)
        leftkey, rightkey, topkey, bottomkey = self.bkgkeys(loc)
        
        #print([(roi, roi.isEditable()) for roi in self.rois])
        
        #croi:
        origin =(key[0].start, key[1].start)
        size = (key[0].stop - key[0].start, key[1].stop - key[1].start)
        #leftkey, rightkey, topkey, bottomkey
        left = leftkey[0].stop - leftkey[0].start
        right = rightkey[0].stop - rightkey[0].start
        top = topkey[1].stop - topkey[1].start
        bottom = bottomkey[1].stop - bottomkey[1].start
        roi.setGeometry(origin=origin, size=size, left=left, right=right, top=top, bottom=bottom)
        roi.setVisible(True)
        #self.roiManager._roisUpdated()
        
    def integrateROI(self):

        if self.scanSelector.scanstab.currentIndex() == 2:
            self.rocking_extraction()
            return
            
        try:
            image = self.fscan.get_raw_img(0)
        except Exception as e:
            print("no images found! %s" % e)
            return
        if self.database.nxfile is None:
            print("No database available")
            return
        dc = self.ubcalc.detectorCal
        #mu = self.ubcalc.mu
        angles = self.ubcalc.angles

        H_1 = np.array([h.value() for h in self.scanSelector.H_1])
        H_0 = np.array([h.value() for h in self.scanSelector.H_0])
        
        imgmask = None
        
        if self.scanSelector.useMaskBox.isChecked():
            if self.centralPlot.getMaskToolsDockWidget().getSelectionMask() is None:
                btn = qt.QMessageBox.question(self,"No mask available","""No mask was selected with the masking tool.
Do you want to continue without mask?""")
                if btn != qt.QMessageBox.Yes:
                    return
            else:
                imgmask = self.centralPlot.getMaskToolsDockWidget().getSelectionMask() > 0.
        
        corr = self.scanSelector.useSolidAngleBox.isChecked() or\
            self.scanSelector.usePolarizationBox.isChecked()
        
        if corr:
            C_arr = np.ones(dc.detector.shape,dtype=np.float64)
            if self.scanSelector.useSolidAngleBox.isChecked():
                C_arr /= dc.solidAngleArray()
            if self.scanSelector.usePolarizationBox.isChecked():
                C_arr /= dc.polarization(factor=dc._polFactor,axis_offset=dc._polAxis)

        
        hkl_del_gam_s1, hkl_del_gam_s2 = self.getROIloc()
        
        nodatapoints = len(self.fscan)
        #print(hkl_del_gam_1s.shape)
        
        if hkl_del_gam_s1.shape[0] == 1:
            hkl_del_gam_1 = np.zeros((nodatapoints,hkl_del_gam_s1.shape[1]), dtype=np.float64)
            hkl_del_gam_2 = np.zeros((nodatapoints,hkl_del_gam_s1.shape[1]), dtype=np.float64)
            hkl_del_gam_1[:] = hkl_del_gam_s1[0]
            hkl_del_gam_2[:] = hkl_del_gam_s2[0]
        else:
            hkl_del_gam_1, hkl_del_gam_2 = hkl_del_gam_s1, hkl_del_gam_s2

        
        dataavail = np.logical_or(hkl_del_gam_1[:,-1],hkl_del_gam_2[:,-1])

        croi1_a = np.zeros_like(dataavail,dtype=np.float64)
        cpixel1_a = np.zeros_like(dataavail,dtype=np.float64)
        bgroi1_a = np.zeros_like(dataavail,dtype=np.float64)
        bgpixel1_a = np.zeros_like(dataavail,dtype=np.float64)
        x_coord1_a = hkl_del_gam_1[:,6]
        y_coord1_a = hkl_del_gam_1[:,7]
        
        croi2_a = np.zeros_like(dataavail,dtype=np.float64)
        cpixel2_a = np.zeros_like(dataavail,dtype=np.float64)
        bgroi2_a = np.zeros_like(dataavail,dtype=np.float64)
        bgpixel2_a = np.zeros_like(dataavail,dtype=np.float64)
        x_coord2_a = hkl_del_gam_2[:,6]
        y_coord2_a = hkl_del_gam_2[:,7]

        progress = qt.QProgressDialog("Integrating images","abort",0,len(self.fscan),self)
        progress.setWindowModality(qt.Qt.WindowModal)
        
        def sumImage(i):
            if not dataavail[i]:
                croi1 = np.nan; croi2 = np.nan
                cpixel1 = np.nan; cpixel2 = np.nan
                bgroi1 = np.nan; bgroi2 = np.nan
                bgpixel1 = np.nan; bgpixel2 = np.nan
            else:
                image = self.fscan.get_raw_img(i).img.astype(np.float64)
                if imgmask is not None:
                    image[imgmask] = np.nan
                    pixelavail = (~imgmask).astype(np.float64)
                else:
                    pixelavail = np.ones_like(image)
                if corr:
                    image *= C_arr
                    
                if hkl_del_gam_1[i,-1]:
                    key = self.intkey(hkl_del_gam_1[i,6:8])
                    bkgkey = self.bkgkeys(hkl_del_gam_1[i,6:8])
                    
                    cimg = image[key[::-1]]
                    
                    # !!!!!!!!!! add mask here  !!!!!!!!!
                    croi1 = np.nansum(cimg)
                    cpixel1 = np.nansum(pixelavail[key[::-1]])
                    bgroi1 = 0.
                    bgpixel1 = 0.
                    for bg in bkgkey:
                        bgimg = image[bg[::-1]]
                        bgroi1 += np.nansum(bgimg)
                        bgpixel1 += np.nansum(pixelavail[bg[::-1]])
                else:
                    croi1 = np.nan
                    cpixel1 = np.nan
                    bgroi1 = np.nan
                    bgpixel1 = np.nan
                
                if hkl_del_gam_2[i,-1]:
                    key = self.intkey(hkl_del_gam_2[i,6:8])
                    bkgkey = self.bkgkeys(hkl_del_gam_2[i,6:8])
                    
                    cimg = image[key[::-1]]
                    # !!!!!!!!!! add mask here  !!!!!!!!!
                    croi2 = np.nansum(cimg)
                    cpixel2 = np.nansum(pixelavail[key[::-1]])
                    bgroi2 = 0.
                    bgpixel2 = 0.
                    for bg in bkgkey:
                        bgimg = image[bg[::-1]]
                        bgroi2 += np.nansum(bgimg)
                        bgpixel2 += np.nansum(pixelavail[bg[::-1]])

                else:
                    croi2 = np.nan
                    cpixel2 = np.nan
                    bgroi2 = np.nan
                    bgpixel2 = np.nan

            return (croi1, cpixel1, bgroi1, bgpixel1), (croi2, cpixel2, bgroi2, bgpixel2)
                

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.numberthreads) as executor: # speedup only for the file reads 
            futures = {}
            for i in range(len(self.fscan)):
                futures[executor.submit(sumImage, i)] = i
            
            for f in concurrent.futures.as_completed(futures):
                try:
                    (croi1, cpixel1, bgroi1, bgpixel1), (croi2, cpixel2, bgroi2, bgpixel2) = f.result()
                    i = futures[f]
                    croi1_a[i] = croi1
                    cpixel1_a[i] = cpixel1
                    bgroi1_a[i] = bgroi1
                    bgpixel1_a[i] = bgpixel1
                    croi2_a[i] = croi2
                    cpixel2_a[i] = cpixel2
                    bgroi2_a[i] = bgroi2
                    bgpixel2_a[i] = bgpixel2
                    
                    progress.setValue(futures[f])
                except concurrent.futures.CancelledError:
                    pass
                except Exception as e:
                    print("Cannot read image:\n%s" % traceback.format_exc())

                if progress.wasCanceled():
                    [f.cancel() for f in futures]
                    break
            
            
        progress.setValue(len(self.fscan))
        
        if np.any(bgpixel1_a):
            croibg1_a = croi1_a - (cpixel1_a/bgpixel1_a) * bgroi1_a
            croibg1_err_a = np.sqrt(croi1_a + ((cpixel1_a/bgpixel1_a)**2)  * bgroi1_a)
        else:
            croibg1_a = croi1_a
            croibg1_err_a = np.sqrt(croi1_a)
            
        if np.any(bgpixel2_a):
            croibg2_a = croi2_a - (cpixel2_a/bgpixel2_a) * bgroi2_a
            croibg2_err_a = np.sqrt(croi2_a + ((cpixel2_a/bgpixel2_a)**2) * bgroi2_a)
        else:
            croibg2_a = croi2_a
            croibg2_err_a = np.sqrt(croi2_a)

        rod_mask1 = np.isfinite(croibg1_a)
        rod_mask2 = np.isfinite(croibg2_a)
        
        s1_masked = hkl_del_gam_1[:,5][rod_mask1]
        s2_masked = hkl_del_gam_2[:,5][rod_mask2]
        
        croibg1_a_masked = croibg1_a[rod_mask1]
        croibg2_a_masked = croibg2_a[rod_mask2]
        
        croibg1_err_a_masked = croibg1_err_a[rod_mask1]
        croibg2_err_a_masked = croibg2_err_a[rod_mask2]
        
        #name = str(H_1) + "*s+" + str(H_0)
        if self.scanSelector.scanstab.currentIndex() == 1:
            x = self.scanSelector.xy_static[0].value()
            y = self.scanSelector.xy_static[1].value()
            name1 = "pixloc[%.2f %.2f]" % (x,y)
            name2 = "pixloc[%.2f %.2f]_2" % (x,y) # does not exist, Just for compatibility
            traj1 = {
                "@NX_class": u"NXcollection",
                "@direction" : u"Fixed pixel coordinates",
                "s" : hkl_del_gam_1[:,5]
            }
            traj2 = {
                "@NX_class": u"NXcollection",
                "@direction" : u"Fixed pixel coordinates",
                "s" : hkl_del_gam_2[:,5]
            }
        else:
            name1 = str(H_1) + "*s1+" + str(H_0)
            name2 = str(H_1) + "*s2+" + str(H_0)
            traj1 = {
                "@NX_class": u"NXcollection",
                "@direction" : u"Intergrated along H_1*s + H_0 in reciprocal space",
                "H_1"  : H_1,
                "H_0" : H_0,
                "s" : hkl_del_gam_1[:,5]
            }
            traj2 = {
                "@NX_class": u"NXcollection",
                "@direction" : u"Intergrated along H_1*s + H_0 in reciprocal space",
                "H_1"  : H_1,
                "H_0" : H_0,
                "s" : hkl_del_gam_2[:,5]
            }

        
        defaultS1 = croibg1_a_masked.size > croibg2_a_masked.size
        
        if hasattr(self.fscan, "title"):
            title = str(self.fscan.title)
        else:
            title = u"%s-scan" % self.fscan.axisname
        
        mu, om = self.getMuOm()
        if len(np.asarray(om).shape) == 0:
            om = np.full_like(mu,om)
        if len(np.asarray(mu).shape) == 0:
            mu = np.full_like(om,mu)
        
        suffix = ''
        i = 0

        while(self.activescanname + "/measurement/" + name1 + suffix in self.database.nxfile):
            suffix = "_%s" % i
            i += 1
        availname1 = name1 + suffix
        
        suffix = ''
        i = 0
        while(self.activescanname + "/measurement/" + name2 + suffix in self.database.nxfile):
            suffix = "_%s" % i
            i += 1
        
        availname2 = name2 + suffix
                                         
        auxcounters = {"@NX_class": u"NXcollection"}
        for auxname in self.fscan.auxillary_counters:
            if hasattr(self.fscan, auxname):
                cntr = getattr(self.fscan, auxname)
                if cntr is not None:
                    auxcounters[auxname] = cntr
                    
                    
        datas1 = {
            "@NX_class": u"NXdata",
            "sixc_angles": {
                "@NX_class": u"NXpositioner",
                "alpha" : np.rad2deg(mu),
                "omega" :  np.rad2deg(om),
                "theta" :  np.rad2deg(-1*om),
                "delta" : np.rad2deg(hkl_del_gam_1[:,3]),
                "gamma" :  np.rad2deg(hkl_del_gam_1[:,4]),
                "chi" :  np.rad2deg(self.ubcalc.chi),
                "phi" :  np.rad2deg(self.ubcalc.phi),
                "@unit" : u"deg"
            },
            "hkl": {
                "@NX_class": u"NXcollection",
                "h" :  hkl_del_gam_1[:,0],
                "k" :  hkl_del_gam_1[:,1],
                "l" : hkl_del_gam_1[:,2]
            },
            "counters":{
                "@NX_class": u"NXdetector",
                "croibg"  : croibg1_a,
                "croibg_errors" :  croibg1_err_a,
                "croi" :  croi1_a,
                "bgroi"  : bgroi1_a,
                "croi_pix"  : cpixel1_a,
                "bgroi_pix" : bgpixel1_a
            },
            "pixelcoord": {
                "@NX_class": u"NXdetector",
                "x" : x_coord1_a,
                "y"  : y_coord1_a
            },
            "trajectory" : traj1,
            "@signal" : u"counters/croibg",
            "@axes": u"trajectory/s",
            "@title": self.activescanname + "_" + availname1,
            "@orgui_meta": u"roi"
        }
        
        datas2 = {
            "@NX_class": u"NXdata",
            "sixc_angles": {
                "@NX_class": u"NXpositioner",
                "alpha" : np.rad2deg(mu),
                "omega" :  np.rad2deg(om),
                "theta" :  np.rad2deg(-1*om),
                "delta" : np.rad2deg(hkl_del_gam_2[:,3]),
                "gamma" :  np.rad2deg(hkl_del_gam_2[:,4]),
                "chi" :  np.rad2deg(self.ubcalc.chi),
                "phi" :  np.rad2deg(self.ubcalc.phi),
                "@unit" : u"deg"
            },
            "hkl": {
                "@NX_class": u"NXcollection",
                "h" :  hkl_del_gam_2[:,0],
                "k" :  hkl_del_gam_2[:,1],
                "l" : hkl_del_gam_2[:,2]
            },
            "counters":{
                "@NX_class": u"NXdetector",
                "croibg"  : croibg2_a,
                "croibg_errors" :  croibg2_err_a,
                "croi" :  croi2_a,
                "bgroi"  : bgroi2_a,
                "croi_pix"  : cpixel2_a,
                "bgroi_pix" : bgpixel2_a
            },
            "pixelcoord": {
                "@NX_class": u"NXdetector",
                "x" : x_coord2_a,
                "y"  : y_coord2_a
            },
            "trajectory" : traj2,
            "@signal" : u"counters/croibg",
            "@axes": u"trajectory/s",
            "@title": self.activescanname + "_" + availname2,
            "@orgui_meta": u"roi"
        }
            
        data = {self.activescanname:{
                    "instrument": {
                        "@NX_class": u"NXinstrument",
                        "positioners": {
                            "@NX_class": u"NXcollection",
                            self.fscan.axisname: self.fscan.axis
                        }
                    },
                    "auxillary" : auxcounters,
                    "measurement": {
                        "@NX_class": u"NXentry",
                        "@default": availname1 if defaultS1 else availname2,
                    },
                    "title":u"%s" % title,
                    "@NX_class": u"NXentry",
                    "@default": u"measurement/%s" % (availname1 if defaultS1 else availname2),
                    "@orgui_meta": u"scan"
                }
            }
            
        if np.any(cpixel1_a > 0.):
            
            self.integrdataPlot.addCurve(s1_masked,croibg1_a_masked,legend=self.activescanname + "_" + availname1,
                                         xlabel="trajectory/s", ylabel="counters/croibg", yerror=croibg1_err_a_masked)
            
            data[self.activescanname]["measurement"][availname1] = datas1
        if np.any(cpixel2_a > 0.):
            
            self.integrdataPlot.addCurve(s2_masked,croibg2_a_masked,legend=self.activescanname + "_" + availname2,
                                         xlabel="trajectory/s", ylabel="counters/croibg", yerror=croibg2_err_a_masked)
            
            data[self.activescanname]["measurement"][availname2] = datas2
            
        self.database.add_nxdict(data)
        
        
            
        
    def _graphCallback(self,eventdict):
        #print(eventdict)
        if eventdict['event'] == 'mouseDoubleClicked':
            #newReflection = np.array([1,1,1,self.imageno,eventdict['x'],eventdict['y']])
            if self.scanSelector.select_roi_action.isChecked():
                self.scanSelector.set_xy_static_loc(eventdict['x'], eventdict['y'])
                self.scanSelector.select_roi_action.setChecked(False)
            else:
                hkl = self.centralPlot.xyHKLConverter(eventdict['x'],eventdict['y'])[:3]
                self.reflectionSel.addReflection(eventdict,self.imageno,hkl)
            
        if eventdict['event'] == 'markerMoved':
            if eventdict['label'].startswith('__'):
                return
            self.reflectionSel.moveReflection(eventdict['label'],[eventdict['x'],eventdict['y']])
            self.reflectionSel.setReflectionActive(eventdict['label'])
        if eventdict['event'] == 'markerClicked':
            if eventdict['label'].startswith('__'):
                return
            self.reflectionSel.setReflectionActive(eventdict['label'])
        
    def intkey(self, coords):

        vsize = int(self.scanSelector.vsize.value())
        hsize = int(self.scanSelector.hsize.value())
        
        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape
        
        coord_restr = np.clip( np.asarray(coords), [0,0], [dethsize, detvsize])

        
        vhalfsize = vsize // 2
        hhalfsize = hsize // 2
        fromcoords = np.round(np.asarray(coord_restr) - np.array([hhalfsize, vhalfsize]))
        tocoords = np.round(np.asarray(coord_restr) + np.array([hhalfsize, vhalfsize]))
        
        if hsize % 2:
            if coord_restr[0] % 1 < 0.5:
                tocoords[0] += 1
            else:
                fromcoords[0] -= 1
        if vsize % 2:
            if coord_restr[1] % 1 < 0.5:
                tocoords[1] += 1
            else:
                fromcoords[1] -= 1
                
        fromcoords = np.clip( np.asarray(fromcoords), [0,0], [dethsize, detvsize])
        tocoords = np.clip( np.asarray(tocoords), [0,0], [dethsize, detvsize])

        loc = tuple(slice(int(fromcoord), int(tocoord)) for fromcoord, tocoord in zip(fromcoords,tocoords))
        
        #from IPython import embed; embed()

        return loc

    def bkgkeys(self, coords):

        left = int(self.scanSelector.left.value())
        right = int(self.scanSelector.right.value())
        top = int(self.scanSelector.top.value())
        bottom = int(self.scanSelector.bottom.value())
        
        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape
        
        croi = self.intkey(coords)
        hcroikey = croi[0]
        vcroikey = croi[1]
        
        leftkey = (slice(int(np.clip(croi[0].start - left, 0, dethsize)), croi[0].start), croi[1]) 
        rightkey = (slice(croi[0].stop,int(np.clip(croi[0].stop + right, 0, dethsize))), croi[1])
        
        topkey = (croi[0], slice(int(np.clip(croi[1].start - top, 0, detvsize)), croi[1].start))
        bottomkey = (croi[0], slice(croi[1].stop, int(np.clip(croi[1].stop + bottom,0,detvsize)) ))
        return leftkey, rightkey, topkey, bottomkey
        
    def closeEvent(self,event):
        self.database.close()
        super().closeEvent(event)
        
        
class Plot2DHKL(silx.gui.plot.PlotWindow):
    sigKeyPressDelete = qt.pyqtSignal()

    def __init__(self,xyHKLConverter,parent=None,backend=None):
        self.xyHKLConverter = xyHKLConverter
        
        
        posInfo = [
            ('X', lambda x, y: x),
            ('Y', lambda x, y: y),
            ('H', lambda x, y: self.xyHKLConverter(x,y)[0]),
            ('K', lambda x, y: self.xyHKLConverter(x,y)[1]),
            ('L', lambda x, y: self.xyHKLConverter(x,y)[2]),
            ('del', lambda x, y: self.xyHKLConverter(x,y)[3]),
            ('gam', lambda x, y: self.xyHKLConverter(x,y)[4]),
            ('Data', WeakMethodProxy(self._getImageValue))]
        
        super(Plot2DHKL, self).__init__(parent=parent, backend=backend,
                             resetzoom=True, autoScale=False,
                             logScale=False, grid=False,
                             curveStyle=False, colormap=True,
                             aspectRatio=True, yInverted=True,
                             copy=True, save=True, print_=True,
                             control=True, position=posInfo,
                             roi=False, mask=True)
        
        if parent is None:
            self.setWindowTitle('Plot2D')
        self.getXAxis().setLabel('Columns')
        self.getYAxis().setLabel('Rows')

        #if silx.config.DEFAULT_PLOT_IMAGE_Y_AXIS_ORIENTATION == 'downward':
        self.getYAxis().setInverted(True)
        

        self.profile = ProfileToolBar(plot=self)
        self.addToolBar(self.profile)

        self.colorbarAction.setVisible(True)
        self.getColorBarWidget().setVisible(True)

        # Put colorbar action after colormap action
        actions = self.toolBar().actions()
        for action in actions:
            if action is self.getColormapAction():
                break

        self.sigActiveImageChanged.connect(self.__activeImageChanged)
        
    def keyPressEvent(self, event):
        key = event.key()
        if key == qt.Qt.Key_Delete and not event.isAutoRepeat():
            self.sigKeyPressDelete.emit()
        super(Plot2DHKL, self).keyPressEvent(event)
        
    def setXyHKLconverter(self,xyHKLConverter):
        self.xyHKLConverter = xyHKLConverter
        
    def __activeImageChanged(self, previous, legend):
        """Handle change of active image

        :param Union[str,None] previous: Legend of previous active image
        :param Union[str,None] legend: Legend of current active image
        """
        if previous is not None:
            item = self.getImage(previous)
            if item is not None:
                item.sigItemChanged.disconnect(self.__imageChanged)

        if legend is not None:
            item = self.getImage(legend)
            item.sigItemChanged.connect(self.__imageChanged)

        positionInfo = self.getPositionInfoWidget()
        if positionInfo is not None:
            positionInfo.updateInfo()

    def __imageChanged(self, event):
        """Handle update of active image item

        :param event: Type of changed event
        """
        if event == items.ItemChangedType.DATA:
            positionInfo = self.getPositionInfoWidget()
            if positionInfo is not None:
                positionInfo.updateInfo()

    def _getImageValue(self, x, y):
        """Get status bar value of top most image at position (x, y)

        :param float x: X position in plot coordinates
        :param float y: Y position in plot coordinates
        :return: The value at that point or '-'
        """
        value = '-'
        valueZ = -float('inf')
        mask = 0
        maskZ = -float('inf')

        for image in self.getAllImages():
            data = image.getData(copy=False)
            isMask = isinstance(image, items.MaskImageData)
            if isMask:
                zIndex = maskZ
            else:
                zIndex = valueZ
            if image.getZValue() >= zIndex:
                # This image is over the previous one
                ox, oy = image.getOrigin()
                sx, sy = image.getScale()
                row, col = (y - oy) / sy, (x - ox) / sx
                if row >= 0 and col >= 0:
                    # Test positive before cast otherwise issue with int(-0.5) = 0
                    row, col = int(row), int(col)
                    if (row < data.shape[0] and col < data.shape[1]):
                        v, z = data[row, col], image.getZValue()
                        if not isMask:
                            value = v
                            valueZ = z
                        else:
                            mask = v
                            maskZ = z
        if maskZ > valueZ and mask > 0:
            return value, "Masked"
        return value

    def _getImageDims(self, *args):
        activeImage = self.getActiveImage()
        if (activeImage is not None and
                    activeImage.getData(copy=False) is not None):
            dims = activeImage.getData(copy=False).shape[1::-1]
            return 'x'.join(str(dim) for dim in dims)
        else:
            return '-'

    def getProfileToolbar(self):
        """Profile tools attached to this plot

        See :class:`silx.gui.plot.Profile.ProfileToolBar`
        """
        return self.profile

    #@deprecated(replacement="getProfilePlot", since_version="0.5.0")
    def getProfileWindow(self):
        return self.getProfilePlot()

    def getProfilePlot(self):
        """Return plot window used to display profile curve.

        :return: :class:`Plot1D`
        """
        return self.profile.getProfilePlot()

class QImportScanCreator(qt.QDialog):
    
    def __init__(self,defaultMuTh, parent=None):
        qt.QDialog.__init__(self, parent)
        self.defaultMuTh = defaultMuTh
        
        layout = qt.QGridLayout()
        
        layout.addWidget(qt.QLabel("scan axis:"),0,0)
        layout.addWidget(qt.QLabel("axis start:"),1,0)
        layout.addWidget(qt.QLabel("axis end:"),2,0)
        self.fixed_label = qt.QLabel("mu (fixed):")
        layout.addWidget(self.fixed_label,3,0)
        layout.addWidget(qt.QLabel("no frames:"),5,0)

        self.scanaxis = qt.QComboBox()
        self.scanaxis.addItem("theta")
        self.scanaxis.addItem("mu")
        self.scanaxis.setCurrentIndex(0)
        self.scanaxis.currentIndexChanged.connect(self.onScanAxisChanged)

        self.omstart = qt.QDoubleSpinBox()
        self.omstart.setRange(-180,180)
        self.omstart.setDecimals(4)
        self.omstart.setSuffix(u" °")
        self.omstart.setValue(-90)
        
        self.omend = qt.QDoubleSpinBox()
        self.omend.setRange(-180,180)
        self.omend.setDecimals(4)
        self.omend.setSuffix(u" °")
        self.omend.setValue(90)
        
        self.no = qt.QSpinBox()
        self.no.setReadOnly(True)
        self.no.setRange(1,1000000000)
        self.no.setValue(180)
        
        self.fixedAngle = qt.QDoubleSpinBox()
        self.fixedAngle.setRange(-180,180)
        self.fixedAngle.setValue(self.defaultMuTh[0])
        
        layout.addWidget(self.scanaxis,0,1)
        layout.addWidget(self.omstart,1,1)
        layout.addWidget(self.omend,2,1)
        layout.addWidget(self.no,5,1)
        layout.addWidget(self.fixedAngle,3,1)
        
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,6,0,-1,-1)

        test = qt.QLabel("Parameters determined from loaded scan:")
        layout.addWidget(test,4,0,1,2)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def onScanAxisChanged(self, index):
        if index == 0:
            self.omstart.setValue(-90.)
            self.omend.setValue(90.)
            self.fixed_label.setText("mu (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[0])
            
        elif index == 1:
            self.omstart.setValue(0.)
            self.omend.setValue(15.)
            self.fixed_label.setText("theta (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[1])

class QScanCreator(qt.QDialog):
    
    def __init__(self,defaultMuTh, parent=None):
        qt.QDialog.__init__(self, parent)
        self.defaultMuTh = defaultMuTh
        
        layout = qt.QGridLayout()
        
        layout.addWidget(qt.QLabel("scan axis:"),0,0)
        layout.addWidget(qt.QLabel("axis start:"),1,0)
        layout.addWidget(qt.QLabel("axis end:"),2,0)
        layout.addWidget(qt.QLabel("no points:"),3,0)
        self.fixed_label = qt.QLabel("mu (fixed):")
        layout.addWidget(self.fixed_label,4,0)

        self.scanaxis = qt.QComboBox()
        self.scanaxis.addItem("theta")
        self.scanaxis.addItem("mu")
        self.scanaxis.setCurrentIndex(0)
        self.scanaxis.currentIndexChanged.connect(self.onScanAxisChanged)

        self.omstart = qt.QDoubleSpinBox()
        self.omstart.setRange(-180,180)
        self.omstart.setDecimals(4)
        self.omstart.setSuffix(u" °")
        self.omstart.setValue(-90)
        
        self.omend = qt.QDoubleSpinBox()
        self.omend.setRange(-180,180)
        self.omend.setDecimals(4)
        self.omend.setSuffix(u" °")
        self.omend.setValue(90)
        
        self.no = qt.QSpinBox()
        self.no.setRange(1,1000000000)
        self.no.setValue(180)
        
        self.fixedAngle = qt.QDoubleSpinBox()
        self.fixedAngle.setRange(-180,180)
        self.fixedAngle.setValue(self.defaultMuTh[0])
        
        layout.addWidget(self.scanaxis,0,1)
        layout.addWidget(self.omstart,1,1)
        layout.addWidget(self.omend,2,1)
        layout.addWidget(self.no,3,1)
        layout.addWidget(self.fixedAngle,4,1)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,5,0,-1,-1)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def onScanAxisChanged(self, index):
        if index == 0:
            self.omstart.setValue(-90.)
            self.omend.setValue(90.)
            self.fixed_label.setText("mu (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[0])
            
        elif index == 1:
            self.omstart.setValue(0.)
            self.omend.setValue(15.)
            self.fixed_label.setText("theta (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[1])


class QDiffractometerImageDialog(qt.QDialog):
    def __init__(self, parent=None):
        qt.QDialog.__init__(self, parent)
        verticalLayout = qt.QVBoxLayout(self)
        verticalLayout.setContentsMargins(0, 0, 0, 0)
        img = qutils.AspectRatioPixmapLabel(self)
        pixmp = qt.QPixmap(resources.getDiffractometerPath())
        #img.setScaledContents(False)
        img.setPixmap(pixmp)
        
        verticalLayout.addWidget(img)

        #reader = qt.QImageReader()
        #img = reader.read()
        #view = qt.QGraphicsView(self)
        #svg = qt.QSvgWidget(self)
        #svg.load(resources.getDiffractometerPath())

        #verticalLayout.addWidget(svg)
        self.setLayout(verticalLayout)


            
class AboutDialog(qt.QDialog):
    def __init__(self,version, msg='' ,parent=None):
        qt.QDialog.__init__(self, parent)
        layout = qt.QVBoxLayout()
        self.setWindowTitle("About orGUI")
        
        pixmap = resources.getSplashScreen(str(version))
        self.logo = qt.QLabel()
        app = qt.QApplication.instance()
        screenGeometry = app.primaryScreen().availableGeometry()
        splashpm = pixmap.scaledToHeight(int(screenGeometry.height()/5), qt.Qt.SmoothTransformation)
        self.logo.setPixmap(splashpm)
        
        messageStr = "orGUI version %s" % version
        messageStr += msg
        messageStr += """<br> <br>
Copyright (c) 2020-2024 Timo Fuchs, published under MIT License
<br> <br>
orGUI: Orientation and Integration with 2D detectors.<br>
Zenodo. <a href=\"https://doi.org/10.5281/zenodo.12592485\">https://doi.org/10.5281/zenodo.12592485</a> <br> <br> 
New software updates will be published under <a href=\"https://doi.org/10.5281/zenodo.12592485\">Zenodo</a>.
<br> <br>
Help requests can be send via Email to Timo Fuchs. 
<br> <br>
"orGUI" was developed during the PhD work of Timo Fuchs,<br>
within the group of Olaf Magnussen.
""" 
        self.label = qt.QLabel()
        self.label.setText(messageStr)
        self.label.setTextInteractionFlags(qt.Qt.TextBrowserInteraction)
        self.label.setTextFormat(qt.Qt.RichText)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok)
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        layout.addWidget(self.logo)
        layout.addWidget(self.label)
        layout.addWidget(buttons)
        self.setLayout(layout)
            
        
class UncaughtHook(qt.QObject):
    #_exception_caught = qt.Signal(object)
 
    def __init__(self, *args, **kwargs):
        super(UncaughtHook, self).__init__(*args, **kwargs)

        # this registers the exception_hook() function as hook with the Python interpreter
        sys.excepthook = self.exception_hook
        
        self.orgui = None

        # connect signal to execute the message box function always on main thread
        #self._exception_caught.connect(show_exception_box)

    def set_orgui(self,orgui):
        self.orgui = orgui
 
    def exception_hook(self, exc_type, exc_value, exc_traceback):
        """Function handling uncaught exceptions.
        It is triggered each time an uncaught exception occurs. 
        """
        if issubclass(exc_type, KeyboardInterrupt):
            # ignore keyboard interrupt to support console applications
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
        else:
            exc_info = (exc_type, exc_value, exc_traceback)
            log_msg = '\n'.join([''.join(traceback.format_tb(exc_traceback)),
                                 '{0}: {1}'.format(exc_type.__name__, exc_value)])
            print("Uncaught exception:\n {0}".format(log_msg))# exc_info=exc_info)

            # trigger message box show
            #self._exception_caught.emit(log_msg)

            if qt.QApplication.instance() is not None:
                if self.orgui is None:
                    errorbox = qt.QMessageBox(qt.QMessageBox.Critical, 
                                              "Uncaught Exception",
                                              "An unexpected error occured. The program will terminate now:\n{0}".format(log_msg),
                                              qt.QMessageBox.Ok)
                    errorbox.exec()
                    sys.exit(1)
                else:
                    resBtn = qutils.critical_detailed_message(self.orgui, "Uncaught Exception", "An unexpected error has occured.\norGUI will terminate now.\nDo you want to try to save the database before terminating?" ,log_msg, qt.QMessageBox.Save | qt.QMessageBox.Discard)
                    #errorbox = qt.QMessageBox(qt.QMessageBox.Critical, 
                    #                          "Uncaught Exception",
                    #                          "An unexpected error occured:\n{0}\nDo you want to try to save the data before terminating?".format(log_msg),
                    #                           qt.QMessageBox.Save | qt.QMessageBox.Discard)
                                              
                    #resBtn = errorbox.exec()
                    
                    if resBtn == qt.QMessageBox.Save:
                        try:
                            self.orgui.database.onSaveDBFile()
                            
                        except Exception:
                            print("Fatal error: Cannot save database:\n%s" % traceback.format_exc())
                            qutils.critical_detailed_message(self.orgui, "Fatal error", "Cannot save database." ,traceback.format_exc())
                    self.orgui.database.close()
            else:
                print("No QApplication instance available.")
            sys.exit(1)
        
def main(configfile):

    a = qt.QApplication(['orGUI'])
    
    qt_exception_hook = UncaughtHook()
    
    mainWindow = orGUI(configfile)
    qt_exception_hook.set_orgui(mainWindow)
    mainWindow.show()
    #a.lastWindowClosed.connect(a.quit)
    return a.exec_()
    
            
if __name__ == '__main__':
    main("./config")

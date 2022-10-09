# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2020 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module description

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
from .ArrayTableDialog import ArrayTableDialog
from .database import DataBase
from ..backend.scans import SimulationScan
from ..backend import backends
from .. import resources

import numpy as np
from datautils.xrayutils import HKLVlieg, CTRcalc
from datautils.xrayutils import ReciprocalNavigation as rn

import sys
from datautils.xrayutils.id31_tools import BlissScan_EBS, Fastscan, BlissScan

QTVERSION = qt.qVersion()
DEBUG = 0

silx.config.DEFAULT_PLOT_SYMBOL = '.'


class orGUI(qt.QMainWindow):
    def __init__(self,configfile,parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.h5database = None # must be a h5py file-like, by default not opened to avoid reading issues at beamtimes!
        self.images_loaded = False
        self.resetZoom = True
        #icon = resources.getQicon("sum_image.svg")
        self.fscan = None
        self.activescanname = "scan"
        self.numberthreads = int(min(os.cpu_count(), 8)) if os.cpu_count() is not None else 1 
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
        
        self.setCentralWidget(maincentralwidget)
        
        
        # Create the object controlling the ROIs and set it up
        self.roiManager = RegionOfInterestManager(self.centralPlot)
        self.roiManager.setColor('pink')  # Set the color of ROI
        
        #self.roiTable = RegionOfInterestTableWidget()
        #self.roiTable.setRegionOfInterestManager(self.roiManager)
        
        #roi order: left, top, right, bottom,  croi
        self.rois = []
        for i in range(5):
            
            roi = RectangleROI()
            roi.setGeometry(origin=(0, 0), size=(0, 0))
            if i == 4:
                roi.setColor('red') # center color
            else:
                roi.setColor('pink') # bg color

                #roi.setEditable(False)
            roi.setLineWidth(2)
            roi.setLineStyle('-')
            roi.setVisible(True)

            self.roiManager.addRoi(roi,useManagerColor=False)
            self.rois.append(roi)
        
        #self.reflTable.view._model.dataChanged.connect(printmodel)
        #self.reflTable.setArrayData(np.array([0,0,0,0,10,10],dtype=np.float))
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
        
        self.loadImagesAct = file.addAction("load scan images")
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
        aboutAct = helpmenu.addAction("About")
        aboutAct.triggered.connect(self._onShowAbout)
        
        #hackAct = helpmenu.addAction("Hack")
        #hackAct.triggered.connect(self.hack_rocking_extraction)
        
        aboutQtAct = helpmenu.addAction("About Qt")
        aboutQtAct.triggered.connect(lambda : qt.QMessageBox.aboutQt(self))
        
        self.setMenuBar(menu_bar)
        
    def hack_rocking_extraction(self):
        # To delete when not needed anymore!
        raise Exception("FOOOO!")
        scans = [(j,'ascan_%s' % j) for j in range(1,15)]
        lar = np.loadtxt("C:/Timo_loc/CH5523/ch5523/CTR/XRR_orgui/ctr_straj.txt").T[0]
        for no, sc in scans:
            ddict = {'name' : sc, 'scanno' : no}
            ddict['event'] = "itemDoubleClicked"
            ddict['file'] = self.hdffile
            #ddict['node'] = obj
            ddict['beamtime'] = 'ch5523'
            self._onScanChanged(ddict)
            #lar = np.linspace(0.28, 3.88,181)
            for l in lar:
                self.scanSelector.H_0[2].setValue(l)
                self.integrateROI()
    
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
        qt.QMessageBox.about(self, "About orGUI", 
        """Copyright (c) 2021 Timo Fuchs, Olaf Magnussen all rights reserved

For internal use only. 
Do not redistribute!

"orGUI" and its dependeny "datautils" were developed during the PhD work of Timo Fuchs,
within the group of Olaf Magnussen. Usage within the group is hereby granted.
""")

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
        xmirror = np.array(xmirror).astype(np.float)
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
            pos = [self.ubcalc.mu,delta,gamma,self.imageNoToOmega(refl.imageno),self.ubcalc.chi,self.ubcalc.phi]
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
            
    def onSearchHKLforStaticROI(self, hkl):
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
        refl_dialog = QReflectionAnglesDialog(refldict,"Select reflection location", self)
        if qt.QDialog.Accepted == refl_dialog.exec():
            for i, cb in enumerate(refl_dialog.checkboxes,1):
                if cb.isChecked():
                    xy = refldict['xy_%s' % i]
                    self.scanSelector.set_xy_static_loc(xy[0], xy[1])
                    return
        

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
            
        """
        try:
            imageno = self.omegaToImageNo(omega)
        except:
            warnings.warn("Not xmirrored: Didn't find the corresponding image")
            [hkl,x,y,omega] = self.ubcalc.calcReflection(hkl,True)
            mirrored = [hkl,x,y,omega]
            try:
                imageno = self.omegaToImageNo(omega)
            except Exception as e:
                errormsg = "[hkl, x, y, om]\nnot mirrored: %s\nmirrored: %s" % (notmirrored,mirrored)
                qutils.warning_detailed_message(self, "Could not find reflection","Didn't find the corresponding reflection on any image.\nError: %s\nShould be at location%s" % (str(e),errormsg), traceback.format_exc())
                #qt.QMessageBox.warning(self,"Could not find reflection","Didn't find the corresponding reflection on any image.\nError: %s\nShould be at location%s" % (str(e),errormsg))
                return
        eventdict = {'x' : x, 'y': y}
        self.reflectionSel.addReflection(eventdict,imageno,hkl)
        """
            
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
        diag.shape1.setValue(self.ubcalc.detectorCal.detector.shape[0])
        diag.shape2.setValue(self.ubcalc.detectorCal.detector.shape[1])
        if diag.exec() == qt.QDialog.Accepted:
            shape = (diag.shape1.value(), diag.shape2.value())
            self.ubcalc.detectorCal.detector.shape = shape
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
        
        
    def _onScanChanged(self,sel_list):
        self.resetZoom = True
        print(sel_list)
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
            self.updateROI()
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
        
    
        
    def updateROI(self):
        #dc = self.ubcalc.detectorCal
        #mu = self.ubcalc.mu
        #angles = self.ubcalc.angles
        
        #H_1 = np.array([h.value() for h in self.scanSelector.H_1])
        #H_0 = np.array([h.value() for h in self.scanSelector.H_0])
        try:
            hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(self.imageno)
        except:
            for roi in self.rois:
                roi.setVisible(False)
            self.centralPlot.removeMarker('main_croi_loc')
            return

        
        """
        hkl_del_gam_1, hkl_del_gam_2 = angles.anglesIntersectLineEwald(H_0, H_1, mu,self.imageNoToOmega(self.imageno),self.ubcalc.phi,self.ubcalc.chi)
        
        delta1 = hkl_del_gam_1[...,3]
        delta2 = hkl_del_gam_2[...,3]
        gam1 = hkl_del_gam_1[...,4]
        gam2 = hkl_del_gam_2[...,4]
        
        yx1 = dc.pixelsSurfaceAngles(gam1, delta1, mu)
        yx2 = dc.pixelsSurfaceAngles(gam2, delta2, mu)
        
        ymask1 = np.logical_and(yx1[...,0] >= 0, yx1[...,0] < dc.detector.shape[0])
        xmask1 = np.logical_and(yx1[...,1] >= 0, yx1[...,1] < dc.detector.shape[1])
        yxmask1 = np.logical_and(xmask1,ymask1)
    
        ymask2 = np.logical_and(yx2[...,0] >= 0, yx2[...,0] < dc.detector.shape[0])
        xmask2 = np.logical_and(yx2[...,1] >= 0, yx2[...,1] < dc.detector.shape[1])
        yxmask2 = np.logical_and(xmask2,ymask2)
        """
        if not self.roivisible:
            for roi in self.rois:
                roi.setVisible(False)
            self.roiManager._roisUpdated()
            self.centralPlot.removeMarker('main_croi_loc')
            
            
        if hkl_del_gam_1[0,-1] or hkl_del_gam_2[0,-1]:
            if hkl_del_gam_1[0,-1]:
                if self.roivisible:
                    self.plotROI(hkl_del_gam_1[0,6:8])
                for i, spinbox in enumerate(self.scanSelector.hkl_static):
                    spinbox.setValue(hkl_del_gam_1[0,i])

            if hkl_del_gam_2[0,-1]:
                if self.roivisible:
                    self.plotROI(hkl_del_gam_2[0,6:8])
                for i, spinbox in enumerate(self.scanSelector.hkl_static):
                    spinbox.setValue(hkl_del_gam_1[0,i])
        else:
            for roi in self.rois:
                roi.setVisible(False)
            self.centralPlot.removeMarker('main_croi_loc')
        
    def getROIloc(self, imageno=None, H_0=None, H_1=None, **kwargs):
        if self.fscan is None:
            raise Exception("No scan loaded!")
        

        mu, om = self.getMuOm(imageno)
            
        dc = self.ubcalc.detectorCal
        #mu = self.ubcalc.mu
        angles = self.ubcalc.angles
        
        """
                    gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([y]),np.array([x]), mu)
            #print(self.ubcalc.detectorCal)
            #print(x,y)
            #print(self.ubcalc.detectorCal.tth(np.array([y]),np.array([x])))
            pos = [mu,delta[0],gamma[0],om,self.ubcalc.chi,self.ubcalc.phi]
            pos = HKLVlieg.crystalAngles(pos,self.ubcalc.n)
            hkl = np.concatenate((np.array(self.ubcalc.angles.anglesToHkl(*pos)),np.rad2deg([delta[0],gamma[0]])))
        """
        #print(self.scanSelector.scanstab.currentIndex())
        if self.scanSelector.scanstab.currentIndex() == 1 and not kwargs.get('intersect', False):
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
            
            hkl_del_gam_1, hkl_del_gam_2, Qa_1, Qa_2 = angles.anglesIntersectLineEwald(H_0, H_1, mu, om, self.ubcalc.phi,self.ubcalc.chi, Qalpha=True)
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

    def plotROI(self, loc):

        key = self.intkey(loc)
        bkgkey = self.bkgkeys(loc)
        for roi in self.rois:
           roi.setVisible(True)
        
        #print([(roi, roi.isEditable()) for roi in self.rois])
        
        #croi:
        self.rois[4].setGeometry(origin=(key[0].start, key[1].start), size=(key[0].stop - key[0].start, key[1].stop - key[1].start))
        #self.rois[4].setVisible(True)
        for i,k in enumerate(bkgkey,0):
            self.rois[i].setGeometry(origin=(k[0].start, k[1].start), size=( k[0].stop - k[0].start, k[1].stop - k[1].start))

        self.roiManager._roisUpdated()
        
        #print([str(r) for r in self.roiManager.getRois()])
        self.centralPlot.addMarker(loc[0],loc[1],legend='main_croi_loc')
        
    def integrateROI(self):
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
        for auxname in backends.auxillary_counters:
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
            self.reflectionSel.moveReflection(eventdict['label'],[eventdict['x'],eventdict['y']])
            self.reflectionSel.setReflectionActive(eventdict['label'])
        if eventdict['event'] == 'markerClicked':
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
        layout.addWidget(qt.QLabel("Det pixel1:"),5,0)
        layout.addWidget(qt.QLabel("Det pixel2:"),6,0)

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
        
        self.shape1 = qt.QSpinBox()
        self.shape1.setRange(1,1000000000)
        self.shape1.setValue(1)
        
        self.shape2 = qt.QSpinBox()
        self.shape2.setRange(1,1000000000)
        self.shape2.setValue(1)
        
        layout.addWidget(self.scanaxis,0,1)
        layout.addWidget(self.omstart,1,1)
        layout.addWidget(self.omend,2,1)
        layout.addWidget(self.no,3,1)
        layout.addWidget(self.fixedAngle,4,1)
        
        layout.addWidget(self.shape1,5,1)
        layout.addWidget(self.shape2,6,1)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,7,0,-1,-1)
        
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

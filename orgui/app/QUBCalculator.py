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

from io import StringIO
from silx.gui import qt
import pyFAI
import numpy as np

import traceback
from datautils.xrayutils import HKLVlieg, CTRcalc
from datautils.xrayutils import DetectorCalibration
from datautils.xrayutils import unitcells
import warnings
import configparser
import os
import copy
from contextlib import contextmanager

from .ArrayTableDialog import ArrayEditWidget

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

# reflectionhandler must implement the method getReflections

class QUBCalculator(qt.QTabWidget):
    sigNewReflection = qt.pyqtSignal(list)
    sigPlottableMachineParamsChanged = qt.pyqtSignal(list)
    #sigQueryImageChange = qt.pyqtSignal(int)
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self,configfile=None, parent=None):
        qt.QTabWidget.__init__(self, parent=None)
        #self.setOrientation(qt.Qt.Vertical)
        
        self.configdir = os.getcwd()
        #self.mainLayout = qt.QVBoxLayout()
        
        umatrixwidget = qt.QSplitter()
        umatrixwidget.setOrientation(qt.Qt.Vertical)
        
        
        self.reflectionWidget = qt.QSplitter()
        self.reflectionWidget.setOrientation(qt.Qt.Horizontal)
        qt.QLabel("H:",self.reflectionWidget)
        self.Hbox = qt.QDoubleSpinBox(self.reflectionWidget)
        self.Hbox.setRange(-100,100)
        self.Hbox.setDecimals(2)
        qt.QLabel("K:",self.reflectionWidget)
        self.Kbox = qt.QDoubleSpinBox(self.reflectionWidget)
        self.Kbox.setRange(-100,100)
        self.Kbox.setDecimals(2)
        qt.QLabel("L:",self.reflectionWidget)
        self.Lbox = qt.QDoubleSpinBox(self.reflectionWidget)
        self.Lbox.setRange(-100,100)
        self.Lbox.setDecimals(2)
        estimateButton = qt.QPushButton("calculate",self.reflectionWidget)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        estimateButton.setToolTip("calculate position of the reflection with given HKL")
        
        estimateButton.clicked.connect(self._onCalcReflection)
        
        umatrixwidget.addWidget(self.reflectionWidget)
        
        umatrixsplitter = qt.QSplitter()
        umatrixsplitter.setOrientation(qt.Qt.Horizontal)
        
        self.Ueditor = qt.QTextEdit("")
        umatrixsplitter.addWidget(self.Ueditor)
        self.calUButton = qt.QPushButton("calculate U")
        self.calUButton.setToolTip("calculate orientation matrix based on the given reflections")
                
        vertCalUSplitter = qt.QSplitter()
        vertCalUSplitter.setOrientation(qt.Qt.Vertical)
        
        fitUbox = qt.QGroupBox("fit options (only available with enough reflections)")
        
        self.latnofit = qt.QRadioButton("don't fit lattice")
        self.latnofit.setChecked(True)
        self.latscale = qt.QRadioButton("fit scale of lattice")
        self.latfitall = qt.QRadioButton("fit all lattice parameters")
        
        fitUboxlayout = qt.QVBoxLayout()
        fitUboxlayout.addWidget(self.latnofit)
        fitUboxlayout.addWidget(self.latscale)
        fitUboxlayout.addWidget(self.latfitall)
        fitUboxlayout.addStretch(1)
        
        fitUbox.setLayout(fitUboxlayout)
        
        
        vertCalUSplitter.addWidget(fitUbox)
        
        vertCalUSplitter.addWidget(self.calUButton)
        
        
        
        
        umatrixsplitter.addWidget(vertCalUSplitter)
        
        umatrixwidget.addWidget(umatrixsplitter)
        

        
        
        
        self.calUButton.clicked.connect(self._onCalcU)
        self.addTab(umatrixwidget,"U Matrix")
        
        self.crystalparams = QCrystalParameter()
        #self.crystalparams.setValues(self.crystal,self.n)
        #self.addTab(self.crystalparams,"Crystal")
        
        
        self.crystalparams.sigCrystalParamsChanged.connect(self._onCrystalParamsChanged)
        
        #paramsSplitter.setOrientation(qt.Qt.Horizontal)
        
        self.machineParams = QMachineParameters()
        self.machineParams.sigMachineParamsChanged.connect(self._onMachineParamsChanged)
        self.machineParams.loadConfigButton.clicked.connect(self._onLoadConfig)
        #self.addTab(self.machineParams,"Machine")
        
        self.uedit = QUEdit()
        self.ueditDialog = QUEditDialog(self.uedit)
        
        self.machineDialog = QMachineParametersDialog(self.machineParams)
        self.xtalDialog = QCrystalParameterDialog(self.crystalparams)
        
                
        if configfile is not None:
            if not self.readConfig(configfile):
                self.toFallbackConfig()
        else:
            self.toFallbackConfig()
        
        """
        
        editorSplitter = qt.QSplitter()
        editorSplitter.setOrientation(qt.Qt.Horizontal)

        self.refleditor = qt.QTextEdit("H\tK\tL\tx\ty\timageno\n")
        editorSplitter.addWidget(self.refleditor)
        
        
        fromEditorButton = qt.QPushButton("from editor",self.reflectionWidget)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        fromEditorButton.setToolTip("take refelctions from editor")
        fromEditorButton.clicked.connect(self.reflectionsFromEditor)
        editorSplitter.addWidget(fromEditorButton)
        
        
        self.addWidget(editorSplitter)
        """
        
    def calcReflection(self,hkl,mirrorx=False):
        pos = self.angles.anglesZmode(hkl,self.mu,'in',self.chi,self.phi,mirrorx=mirrorx)
        alpha, delta, gamma, omega, chi, phi = HKLVlieg.vacAngles(pos,self.n)
        print(np.rad2deg([alpha, delta, gamma, omega, chi, phi]))
        y,x = self.detectorCal.pixelsSurfacePoint([gamma],[delta],alpha)[0]
        return [hkl,x,y,omega]
        
    def _onCalcReflection(self):
        hkl = [self.Hbox.value(),self.Kbox.value(),self.Lbox.value()]
        self.sigNewReflection.emit(self.calcReflection(hkl))
        
        
    def _onCrystalParamsChanged(self,crystal,n):
        #a,alpha,_,_ = crystal.getLatticeParameters()
        #self.crystal.setLattice(a,np.rad2deg(alpha))
        
        self.crystal = crystal
        self.n = n
        self.ubCal.setCrystal(self.crystal)
        #self.ubCal.defaultU()
    
    def _onMachineParamsChanged(self,params):
        [E,mu,sdd,pixsize,cp,polax,polf,azim,chi,phi] = params
        self.mu = mu
        self.chi = chi
        self.phi = phi
        self.ubCal.setEnergy(E)
        fit2DCal = self.detectorCal.getFit2D()
        fit2DCal['centerX'] = cp[0]
        fit2DCal['centerY'] = cp[1]
        fit2DCal['directDist'] = sdd*1e3
        f2d = [fit2DCal['directDist'],
               fit2DCal['centerX'],
               fit2DCal['centerY'],
               fit2DCal['tilt'],
               fit2DCal['tiltPlanRotation'],
               pixsize*1e6,
               pixsize*1e6,
               fit2DCal['splineFile']]
        self.detectorCal.setFit2D(*f2d)
        self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
        self.detectorCal.setAzimuthalReference(azim)
        self.azimuth = azim
        self.detectorCal.setPolarization(polax,polf)
        self.crystal.setEnergy(E*1e3)
        try:
            gam_p,_ = self.detectorCal.rangegamdel_p
            azimy,azimx = self.detectorCal.pixelsPrimeBeam(gam_p[1]/5, 0 )[0]
            self.sigPlottableMachineParamsChanged.emit([cp,[azimx,azimy],polax])
        except Exception as e:
            # here is a bug with the init of the detector cal
            print(traceback.format_exc())
            #pass
        #print(self.detectorCal.get_wavelength())
        #print(self.detectorCal.getFit2D())
        
    def _onLoadConfig(self):
        fileTypeDict = {'config files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open config file",
                                                  self.configdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.configdir = os.path.splitext(filename)[0]
        self.readConfig(filename)

    def setReflectionHandler(self,refls):
        self.reflections = refls
        
    def readConfig(self,configfile):
        config = configparser.ConfigParser()
        try:
            if os.path.isfile(configfile):
                config.read(configfile)
            else:
                raise Exception("File does not exist")
        except Exception as e:
            qt.QMessageBox.warning(self,"Can not read config","Can not read config file:\nException occured during read of configfile %s,\nException:\n %s" % (configfile,e))
            return False
        try:
            machine = config['Machine']
            lattice = config['Lattice']
            diffrac = config['Diffractometer']
            
            self.azimuth = np.deg2rad(diffrac.getfloat('azimuthal_reference',0))
            self.polaxis = np.deg2rad(diffrac.getfloat('polarization_axis',0))
            self.polfactor = diffrac.getfloat('polarization_factor',0)
            
            sdd = machine.getfloat('SDD',0.729) #m
            E =  machine.getfloat('E',78.0) #keV
            pixelsize = machine.getfloat('pixelsize',172e-6) #m
            cpx = machine.getfloat('cpx',731)
            cpy = machine.getfloat('cpy',1587)
            cp = [cpx,cpy]
            
            self.mu = np.deg2rad(diffrac.getfloat('mu',0.05))
            self.chi = np.deg2rad(diffrac.getfloat('chi',0.0))
            self.phi = np.deg2rad(diffrac.getfloat('phi',0.0))
            

            a1 = lattice.getfloat('a1',-1)
            a2 = lattice.getfloat('a2',-1)
            a3 = lattice.getfloat('a3',-1)
            alpha1 = lattice.getfloat('alpha1',-1)
            alpha2 = lattice.getfloat('alpha2',-1)
            alpha3 = lattice.getfloat('alpha3',-1)
            self.n = 1 - lattice.getfloat('refractionindex',0.0)
            
            lat = np.array([a1,a2,a3])
            
            latticeoverride = True
            latangle = np.array([alpha1,alpha2,alpha3])
            if np.any(lat < 0.) or np.any(latangle < 0):
                latticeoverride = False
                a1 = a2 = a3 = 1.
                alpha1 = alpha2 = alpha3 = 90.
                print("Fallback lattice vectors")
                
            
            self.crystal = CTRcalc.UnitCell([a1,a2,a3],[alpha1,alpha2,alpha3])
            self.crystal.addAtom('Pt',[0.,0.,0.],0.1,0.1,1.)
            self.crystal.setEnergy(E*1e3)
            
            self.ubCal = HKLVlieg.UBCalculator(self.crystal,E)
            self.ubCal.defaultU()
            self.uedit.setU(self.ubCal.getU())
            self.angles = HKLVlieg.VliegAngles(self.ubCal)
            
            if 'crystal' in lattice:
                idx = self.crystalparams.crystalComboBox.findText(lattice['crystal'],qt.Qt.MatchFixedString)
                if idx == -1:
                    qt.QMessageBox.warning(self,"Did not find crystal","Can not find crystal <%s> \nException occured during read of configfile %s,\nException:\n%s" % (unitcellsconfigfile,traceback.format_exc()))
                else:
                    self.crystalparams.crystalComboBox.setCurrentIndex(idx)
                    self.crystalparams.onSwitchCrystal(idx)
                    #self.crystal = self.crystalparams.getCrystal()

            if latticeoverride:
                print('foo')
                self.crystal.setLattice([a1,a2,a3],[alpha1,alpha2,alpha3])
            
            self.detectorCal = DetectorCalibration.Detector2D_SXRD()
            if 'poni' in machine:
                if machine['poni']:
                    self.detectorCal.load(machine['poni'])
                    self.ubCal.setLambda(self.detectorCal.get_wavelength()*1e10)
                    
                else:
                    self.detectorCal.setFit2D(sdd*1e3,cpx,cpy,pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
                    self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
                    self.detectorCal.detector.shape = (2880,2880) # Perkin 
                    
            else:
                self.detectorCal.setFit2D(sdd*1e3,cpx,cpy,pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
                self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
                self.detectorCal.detector.shape = (2880,2880)
                
            self.detectorCal.setAzimuthalReference(self.azimuth)
            self.detectorCal.setPolarization(self.polaxis,self.polfactor)
            
            fit2dCal = self.detectorCal.getFit2D()

            paramlist = [self.ubCal.getEnergy(),self.mu,fit2dCal['directDist']/1e3,
                         fit2dCal['pixelX']*1e-6,[fit2dCal['centerX'],fit2dCal['centerY']],self.polaxis
                         ,self.polfactor,self.azimuth,self.chi,self.phi]
            self.crystalparams.setValues(self.crystal,self.n)
            self.machineParams.setValues(paramlist)
            return True
            
        except Exception as e:
            qt.QMessageBox.warning(self,"Can not parse config","Can not parse config file:\nException occured during parsing of configfile %s,\nException:\n %s" % (configfile,traceback.format_exc()))
            return False
        
    def toFallbackConfig(self):
        sdd = 0.729 #m
        E = 78.
        pixelsize = 172e-6
        cp = [731.0,1587.856]
        self.mu = np.deg2rad(0.05)
        self.chi = 0.
        self.phi = 0.
        self.n = 1 - 1.1415e-06
        self.crystal = CTRcalc.UnitCell([3.9242,3.9242,3.9242],[90.,90.,90.])
        self.crystal.addAtom('Pt',[0.,0.,0.],0.1,0.1,1.)
        self.crystal.setEnergy(E*1e3)
        
        self.ubCal = HKLVlieg.UBCalculator(self.crystal,E)
        self.ubCal.defaultU()
        self.uedit.setU(self.ubCal.getU())
        self.polaxis = 0
        self.polfactor = 0
        self.azimuth = 0
        self.detectorCal = DetectorCalibration.Detector2D_SXRD()
        self.detectorCal.detector = pyFAI.detector_factory("Pilatus2m")
        self.detectorCal.setFit2D(sdd*1e3,cp[0],cp[1],pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
        self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
        
        fit2dCal = self.detectorCal.getFit2D()
        self.angles = HKLVlieg.VliegAngles(self.ubCal)
        paramlist = [self.ubCal.getEnergy(),self.mu,fit2dCal['directDist']/1e3,
                     fit2dCal['pixelX']*1e-6,[fit2dCal['centerX'],fit2dCal['centerY']],self.polaxis,
                     self.polfactor,self.azimuth,self.chi,self.phi]
        self.crystalparams.setValues(self.crystal,self.n)
        self.machineParams.setValues(paramlist)
        
        
        
    def _onCalcU(self):
        hkls,angles = self.reflections()
        if len(hkls) < 1:
            qt.QMessageBox.warning(self,"Not enough reflections","You must select at least 1 reflection to calculate an orientation matrix")
            return
        elif len(hkls) == 1:
            try:
                self.ubCal.zmodeUSingleRefl(angles[0],hkls[0])
            except Exception as e:
                qt.QMessageBox.critical(self,"Cannot calculate UB matrix","Error during UB matrix calculation:\n%s" % traceback.format_exc())
                return
        else:
            try:
                self.ubCal.setPrimaryReflection(angles[0],hkls[0])
                self.ubCal.setSecondayReflection(angles[1],hkls[1])
                self.ubCal.calculateU()
            except Exception as e:
                qt.QMessageBox.critical(self,"Cannot calculate UB matrix","Error during UB matrix calculation:\n%s" % traceback.format_exc())
                return
        
        if len(hkls) > 2:
            if self.latnofit.isChecked():
                self.ubCal.refineU(hkls,angles)
                #print(self.ubCal.getU())
                
            if self.latscale.isChecked():
                if len(hkls) > 3:
                    self.ubCal.refineULattice(hkls,angles,'scale')
                else:
                    qt.QMessageBox.warning(self,"Not enough reflections","You must select at least 4 reflections to fit lattice scale and U")
            if self.latfitall.isChecked():
                if len(hkls) > 5:
                    self.ubCal.refineULattice(hkls,angles,'lat')
                else:
                    qt.QMessageBox.warning(self,"Not enough reflections","You must select at least 6 reflections to fit lattice and U")
            #print(self.ubCal.getU())
            self.crystalparams.setValues(self.crystal,self.n)
        #print(self.ubCal.getU())
        self.uedit.setU(self.ubCal.getU())
        #self.Ueditor.setPlainText(str(self.ubCal.getU()))
            
        
        
class QCrystalParameter(qt.QSplitter):
    sigCrystalParamsChanged = qt.pyqtSignal(HKLVlieg.Crystal,float)
    def __init__(self,parent=None):
        qt.QSplitter.__init__(self, parent=None)
        self.setOrientation(qt.Qt.Vertical)
        
        latticeParamsGroup = qt.QGroupBox("Lattice parameters")
        latticeParamsLayout = qt.QGridLayout()
        
        self.filedialogdir = '.'
        
        self._uc = CTRcalc.UnitCell([3.9242,3.9242,3.9242],[90.,90.,90.])
        self._uc.addAtom('Pt',[0.,0.,0.],0.1,0.1,1.)
        self._uc.setEnergy(70000.)
        
        self._n = 1.
        
        latticeParamsLayout.addWidget(qt.QLabel("a1:"),0,0)
        latticeParamsLayout.addWidget(qt.QLabel("a2:"),1,0)
        latticeParamsLayout.addWidget(qt.QLabel("a3:"),2,0)
        self.a1box = qt.QDoubleSpinBox()
        self.a1box.setRange(0,100)
        self.a1box.setDecimals(4)
        self.a1box.setSuffix(u" \u212B")
        
        latticeParamsLayout.addWidget(self.a1box,0,1)
        
        self.a2box = qt.QDoubleSpinBox()
        self.a2box.setRange(0,100)
        self.a2box.setDecimals(4)
        self.a2box.setSuffix(u" \u212B")

        latticeParamsLayout.addWidget(self.a2box,1,1)
        
        self.a3box = qt.QDoubleSpinBox()
        self.a3box.setRange(0,100)
        self.a3box.setDecimals(4)
        self.a3box.setSuffix(u" \u212B")
        
        latticeParamsLayout.addWidget(self.a3box,2,1)

        latticeParamsLayout.addWidget(qt.QLabel("alpha1:"),0,2)
        latticeParamsLayout.addWidget(qt.QLabel("alpha2:"),1,2)
        latticeParamsLayout.addWidget(qt.QLabel("alpha3:"),2,2)
        self.alpha1box = qt.QDoubleSpinBox()
        self.alpha1box.setRange(0,180)
        self.alpha1box.setDecimals(2)
        self.alpha1box.setSuffix(" °")
        
        latticeParamsLayout.addWidget(self.alpha1box,0,3)
        
        self.alpha2box = qt.QDoubleSpinBox()
        self.alpha2box.setRange(0,180)
        self.alpha2box.setDecimals(2)
        self.alpha2box.setSuffix(" °")
       
        latticeParamsLayout.addWidget(self.alpha2box,1,3)
        
        self.alpha3box = qt.QDoubleSpinBox()
        self.alpha3box.setRange(0,180)
        self.alpha3box.setDecimals(2)
        self.alpha3box.setSuffix(" °")
        
        latticeParamsLayout.addWidget(self.alpha3box,2,3)
        
        
        latticeParamsGroup.setLayout(latticeParamsLayout)
        
        
        
        refractionindexGroup = qt.QGroupBox("refraction index")
        refractionindexLayout = qt.QHBoxLayout()
        #refractionindexLayout.setOrientation(qt.Qt.Horizontal)
        
        
        refractionindexLayout.addWidget(qt.QLabel("delta / 1e-6:"))
        self.refractionIndexBox = qt.QDoubleSpinBox()
        
        self.refractionIndexBox.setRange(0,1000)
        self.refractionIndexBox.setDecimals(3)
        
        refractionindexLayout.addWidget(self.refractionIndexBox)
        
        refractionindexGroup.setLayout(refractionindexLayout)
        #self.setValues(crystal,n)
        
        
        
        crystalParamsGroup = qt.QGroupBox("Crystal")
        crystalParamsLayout = qt.QGridLayout()
        
        crystalParamsLayout.addWidget(qt.QLabel("Crystal:"),0,0)
        
        self.crystalComboBox = qt.QComboBox()
        crystalParamsLayout.addWidget(self.crystalComboBox,0,1)
        
        for uc in unitcells.availablebulk:
            self.crystalComboBox.addItem(uc, uc)
            
        loadxtalbtn = qt.QPushButton("Load bulk file")
        crystalParamsLayout.addWidget(loadxtalbtn,1,1)
        loadxtalbtn.clicked.connect(self.onLoadXtal)
        
        self.crystalComboBox.activated.connect(self.onSwitchCrystal)
        
        
        crystalParamsGroup.setLayout(crystalParamsLayout)
        
        self.addWidget(crystalParamsGroup)
        self.addWidget(latticeParamsGroup)
        self.addWidget(refractionindexGroup)
        
        
        
        self.a1box.valueChanged.connect(self._onAnyValueChanged)
        self.a2box.valueChanged.connect(self._onAnyValueChanged)
        self.a3box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha1box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha2box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha3box.valueChanged.connect(self._onAnyValueChanged)
        self.refractionIndexBox.valueChanged.connect(self._onAnyValueChanged)
    
    def onLoadXtal(self):
        fileTypeDict = {'bulk files (*.bul)': '.bul', 'Crystal Files (*.xtal)': '.txt', 'All files (*)': '', }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open crystal file with atom locations",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        
        if filetype == 'bulk files (*.bul)':
            try:
                uc_bulk = CTRcalc.UnitCell.fromBULfile(filename)
            except Exception:
                qt.QMessageBox.critical(self,"Cannot open xtal", "Cannot open:\n%s" % traceback.format_exc())
                return
        else:
            qt.QMessageBox.critical(self,"Cannot open xtal", "File extension not understood")
            return
            
        uc_name = os.path.splitext(os.path.basename(filename) )[0]
        self.crystalComboBox.addItem(uc_name, uc_bulk)
        idx = self.crystalComboBox.findText(uc_name)
        self.crystalComboBox.setCurrentIndex(idx)
        self.onSwitchCrystal(idx)
        
    def onSwitchCrystal(self, index):
        selectiondata = self.crystalComboBox.itemData(index)
        if isinstance(selectiondata, str):
            uc = unitcells.unitcell(selectiondata)
            self.crystalComboBox.setItemData(index,uc)
            uc = copy.deepcopy(uc)
        else:
            uc = copy.deepcopy(selectiondata)
        self.setValues(uc,self.getRefractionIndex())
        
    def setValues(self,crystal,n):
        [a1,a2,a3],alpha,_,_ = crystal.getLatticeParameters()
        [alpha1,alpha2,alpha3] = np.rad2deg(alpha)
        signList = [self.a1box, self.a2box, self.a3box,
                    self.alpha1box, self.alpha2box, self.alpha3box,
                    self.refractionIndexBox]
        with blockSignals(signList):
            self.a1box.setValue(a1)
            self.a2box.setValue(a2)
            self.a3box.setValue(a3)
            self.alpha1box.setValue(alpha1)
            self.alpha2box.setValue(alpha2)
            self.alpha3box.setValue(alpha3)
            self.refractionIndexBox.setValue((1.-n)*1e6)
        #self.blockSignals(False)
        self._uc = crystal
        self._onAnyValueChanged()
        
    def getCrystal(self):
        a = np.array([self.a1box.value(),self.a2box.value(),self.a3box.value()])
        alpha = np.array([self.alpha1box.value(),self.alpha2box.value(),self.alpha3box.value()])
        if np.any(a == 0) or np.any(alpha == 0):
            raise Exception("No crystal set")
        self._uc.setLattice(a,alpha)
        if len(self._uc.basis) < 1.:
            self._uc.addAtom('Pt',[0.,0.,0.],0.1,0.1,1.)
        return self._uc
        
    def getRefractionIndex(self):
        return 1 - self.refractionIndexBox.value()*1e-6
        
    def _onAnyValueChanged(self):
        try:
            newCrystal = self.getCrystal()
            n = self.getRefractionIndex()
            self.sigCrystalParamsChanged.emit(newCrystal,n)
        except Exception as e:
            qt.QMessageBox.warning(self,"Can not calculate B Matrix","The B Matrix can not be calculated\nError: %s" % str(e))
                

class QUEdit(qt.QWidget):
    sigUChanged = qt.pyqtSignal(np.ndarray)
    
    def __init__(self,parent=None):
        qt.QWidget.__init__(self, parent=None)
        
        alignGroup = qt.QGroupBox("Align hkl onto xyz")
        orientationLayout = qt.QGridLayout()
        
        self.hkl = []
        for i, index in enumerate(['H', 'K', 'L']):
            orientationLayout.addWidget(qt.QLabel("%s:" % index),i,0)
            milleredit = qt.QDoubleSpinBox()
            milleredit.setRange(0.001,1000)
            milleredit.setDecimals(4)
            orientationLayout.addWidget(milleredit,i,1)
            milleredit.setValue(1. if i == 2 else 0.)
            self.hkl.append(milleredit)
            
        self.override_angles = qt.QCheckBox("Override angles")
        orientationLayout.addWidget(self.override_angles, 3, 0, 1, 2)
            
        self.angles = []
        for i, index in enumerate(['alpha', 'chi', 'phi', 'theta']):
            orientationLayout.addWidget(qt.QLabel("%s:" % index),i,2)
            edit = qt.QDoubleSpinBox()
            edit.setRange(-360,360)
            edit.setDecimals(4)
            edit.setSuffix(u" °")
            edit.setValue(0.)
            orientationLayout.addWidget(edit,i,3)
            self.angles.append(edit)
        
        self.xyz = []
        for i, index in enumerate(['x', 'y', 'z']):
            orientationLayout.addWidget(qt.QLabel("%s:" % index),i,4)
            edit = qt.QDoubleSpinBox()
            edit.setRange(-1000,1000)
            edit.setDecimals(4)
            orientationLayout.addWidget(edit,i,5)
            edit.setValue(1. if i == 2 else 0.)
            self.xyz.append(edit)
            
        alignGroup.setLayout(orientationLayout)
        
        referenceGroup = qt.QGroupBox("Reference frame")
        referenceLayout = qt.QVBoxLayout()
        self.alphaFrame = qt.QRadioButton("surface")
        self.alphaFrame.setChecked(True)
        self.labFrame = qt.QRadioButton("laboratory")
        
        referenceLayout.addWidget(self.alphaFrame)
        referenceLayout.addWidget(self.labFrame)
        referenceLayout.addStretch(1)
        referenceGroup.setLayout(referenceLayout)
        
        alignbtnlayout = qt.QVBoxLayout()
        alignbtnlayout.addWidget(referenceGroup)
        
        self.alignBtn = qt.QPushButton("Align")
        alignbtnlayout.addWidget(self.alignBtn)
        
        editLayout = qt.QHBoxLayout()
        
        editLayout.addWidget(alignGroup)
        editLayout.addLayout(alignbtnlayout)
        
        mainLayout = qt.QVBoxLayout()
        mainLayout.addLayout(editLayout)
        
        uGroup = qt.QGroupBox("Orientation matrix")
        self.uview = ArrayEditWidget(True, 1, False)
        self.uview.model.dataChanged.connect(self.onUChanged)
        la = qt.QVBoxLayout()
        la.addWidget(self.uview)
        uGroup.setLayout(la)
        mainLayout.addWidget(uGroup)
        
        self.setLayout(mainLayout)
        
    def setU(self, U):
        self.uview.setArrayData(U, None, True, True)
        
    def getU(self):
        return self.uview.getData()
        
    def onUChanged(self):
        print("changed")
        self.sigUChanged.emit(self.getU())

        
class QUEditDialog(qt.QDialog):
    sigHide = qt.pyqtSignal()

    def __init__(self,uedit,parent=None):
        qt.QDialog.__init__(self, parent=None)
        self.uedit = uedit
        layout = qt.QVBoxLayout()
        layout.addWidget(uedit)
        
        self.savedU = self.uedit.getU()
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel | qt.QDialogButtonBox.Reset,
                                      qt.Qt.Horizontal)
        layout.addWidget(buttons)
        
        okbtn = buttons.button(qt.QDialogButtonBox.Ok)
        cancelbtn = buttons.button(qt.QDialogButtonBox.Cancel)
        resetbtn = buttons.button(qt.QDialogButtonBox.Reset)
        
        okbtn.clicked.connect(self.hide)
        cancelbtn.clicked.connect(self.onCancel)
        resetbtn.clicked.connect(self.resetParameters)
        
        self.setLayout(layout)
        
    def showEvent(self, event):
        if event.spontaneous():
            super().showEvent(event)
        else:
            self.savedU = self.uedit.getU()
            super().showEvent(event)
            
    def hideEvent(self, event):
        self.sigHide.emit()
        super().hideEvent(event)
            
    def resetParameters(self):
        self.uedit.setU(self.savedU)
        self.uedit.onUChanged()
            
    def onCancel(self):
        self.resetParameters()
        self.hide()
    
        

        
class QMachineParameters(qt.QWidget):
    sigMachineParamsChanged = qt.pyqtSignal(list)
    def __init__(self,parent=None):
        qt.QWidget.__init__(self, parent=None)
        
        #[E,mu,sdd,pixsize,cp,chi,phi] = params
        
        mainLayout = qt.QGridLayout()
        
        mainLayout.addWidget(qt.QLabel("Energy:"),0,0)
        mainLayout.addWidget(qt.QLabel("SDD:"),1,0)
        mainLayout.addWidget(qt.QLabel("mu:"),2,0)
        mainLayout.addWidget(qt.QLabel("chi:"),3,0)
        mainLayout.addWidget(qt.QLabel("azimuth:"),4,0)
        mainLayout.addWidget(qt.QLabel("Pol factor:"),5,0)
        
        mainLayout.addWidget(qt.QLabel("Centr X:"),0,2)
        mainLayout.addWidget(qt.QLabel("Centr Y:"),1,2)
        mainLayout.addWidget(qt.QLabel("pixel size:"),2,2)
        mainLayout.addWidget(qt.QLabel("phi:"),3,2)
        mainLayout.addWidget(qt.QLabel("polarization axis:"),4,2)
        
        self.loadConfigButton = qt.QPushButton("load config",self) 
        self.loadConfigButton.setToolTip("load machine and crystal configuration from configfile,\naccepts poni file from pyFAI")
        mainLayout.addWidget(self.loadConfigButton,5,2)
         
        self.Ebox = qt.QDoubleSpinBox()
        self.Ebox.setRange(0.1,1000)
        self.Ebox.setDecimals(4)
        self.Ebox.setSuffix(u" keV")
        
        mainLayout.addWidget(self.Ebox,0,1)
        
        self.SDDbox = qt.QDoubleSpinBox()
        self.SDDbox.setRange(0.001,100)
        self.SDDbox.setDecimals(4)
        self.SDDbox.setSuffix(u" m")
        mainLayout.addWidget(self.SDDbox,1,1)
        
        self.mubox = qt.QDoubleSpinBox()
        self.mubox.setRange(0,90)
        self.mubox.setDecimals(4)
        self.mubox.setSuffix(u" °")
        
        mainLayout.addWidget(self.mubox,2,1)
        
        self.chibox = qt.QDoubleSpinBox()
        self.chibox.setRange(0,90)
        self.chibox.setDecimals(4)
        self.chibox.setSuffix(u" °")
        
        mainLayout.addWidget(self.chibox,3,1)
        
        self.azimbox = qt.QDoubleSpinBox()
        self.azimbox.setRange(0,360)
        self.azimbox.setDecimals(4)
        self.azimbox.setSuffix(u" °")
        
        mainLayout.addWidget(self.azimbox,4,1)        

        self.polfbox = qt.QDoubleSpinBox()
        self.polfbox.setRange(-1,1)
        self.polfbox.setDecimals(4)
        #self.polfbox.setSuffix(u" °")
        
        mainLayout.addWidget(self.polfbox,5,1)
        
        
        self.cpXbox = qt.QDoubleSpinBox()
        self.cpXbox.setRange(-100000,100000)
        self.cpXbox.setDecimals(4)
        #self.cpXbox.setSuffix(u" keV")
        
        mainLayout.addWidget(self.cpXbox,0,3)
        
        self.cpYbox = qt.QDoubleSpinBox()
        self.cpYbox.setRange(-100000,100000)
        self.cpYbox.setDecimals(4)
        #self.cpYbox.setSuffix(u" m")
        
        mainLayout.addWidget(self.cpYbox,1,3)
        
        self.pixsizebox = qt.QDoubleSpinBox()
        self.pixsizebox.setRange(0.000001,10)
        self.pixsizebox.setDecimals(5)
        self.pixsizebox.setSuffix(u" mm")
        
        mainLayout.addWidget(self.pixsizebox,2,3)
        
        self.phibox = qt.QDoubleSpinBox()
        self.phibox.setRange(0,90)
        self.phibox.setDecimals(4)
        self.phibox.setSuffix(u" °")
        
        mainLayout.addWidget(self.phibox,3,3)
        
        self.polaxbox = qt.QDoubleSpinBox()
        self.polaxbox.setRange(0,360)
        self.polaxbox.setDecimals(4)
        self.polaxbox.setSuffix(u" °")
        
        mainLayout.addWidget(self.polaxbox,4,3)
        
        #self.setValues(params)
        
        self.Ebox.valueChanged.connect(self._onAnyValueChanged)
        self.mubox.valueChanged.connect(self._onAnyValueChanged)
        self.SDDbox.valueChanged.connect(self._onAnyValueChanged)
        self.cpXbox.valueChanged.connect(self._onAnyValueChanged)
        self.cpYbox.valueChanged.connect(self._onAnyValueChanged)
        self.chibox.valueChanged.connect(self._onAnyValueChanged)
        self.phibox.valueChanged.connect(self._onAnyValueChanged)
        self.pixsizebox.valueChanged.connect(self._onAnyValueChanged)
        self.polaxbox.valueChanged.connect(self._onAnyValueChanged)
        self.azimbox.valueChanged.connect(self._onAnyValueChanged)
        self.polfbox.valueChanged.connect(self._onAnyValueChanged)
        
        self.setLayout(mainLayout)
        
    def setValues(self,params):
        [E,mu,sdd,pixsize,cp,polax,polf,azim,chi,phi] = params
        signList = [self.Ebox, self.SDDbox, self.mubox,
                    self.chibox, self.cpXbox, self.cpYbox,
                    self.pixsizebox, self.phibox, self.polaxbox,
                    self.azimbox, self.polfbox]
        with blockSignals(signList):
            self.Ebox.setValue(E)
            self.SDDbox.setValue(sdd)
            self.mubox.setValue(np.rad2deg(mu))
            self.chibox.setValue(np.rad2deg(chi))
            self.cpXbox.setValue(cp[0])
            self.cpYbox.setValue(cp[1])
            self.pixsizebox.setValue(pixsize*1e3)
            self.phibox.setValue(np.rad2deg(phi))
            self.polaxbox.setValue(np.rad2deg(polax))
            self.azimbox.setValue(np.rad2deg(azim))
            self.polfbox.setValue(polf)
        self._onAnyValueChanged()
        
    def getParameters(self):
        E = self.Ebox.value()
        mu = np.deg2rad(self.mubox.value())
        sdd = self.SDDbox.value()
        pixsize = self.pixsizebox.value()*1e-3
        cp = [self.cpXbox.value(),self.cpYbox.value()]
        chi = np.deg2rad(self.chibox.value())
        phi = np.deg2rad(self.phibox.value())
        azim = np.deg2rad(self.azimbox.value())
        polax = np.deg2rad(self.polaxbox.value())
        polf = self.polfbox.value()
        return [E,mu,sdd,pixsize,cp,polax,polf,azim,chi,phi]
        
    def _onAnyValueChanged(self):
        self.sigMachineParamsChanged.emit(self.getParameters())
        
        
        
class QMachineParametersDialog(qt.QDialog):
    sigHide = qt.pyqtSignal()

    def __init__(self,machineparams,parent=None):
        qt.QDialog.__init__(self, parent=None)
        self.machineparams = machineparams
        layout = qt.QVBoxLayout()
        layout.addWidget(machineparams)
        
        self.savedParams = self.machineparams.getParameters()
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel | qt.QDialogButtonBox.Reset,
                                      qt.Qt.Horizontal)
        layout.addWidget(buttons)
        
        okbtn = buttons.button(qt.QDialogButtonBox.Ok)
        cancelbtn = buttons.button(qt.QDialogButtonBox.Cancel)
        resetbtn = buttons.button(qt.QDialogButtonBox.Reset)
        
        okbtn.clicked.connect(self.hide)
        cancelbtn.clicked.connect(self.onCancel)
        resetbtn.clicked.connect(self.resetParameters)
        
        self.setLayout(layout)
        
    def showEvent(self, event):
        if event.spontaneous():
            super().showEvent(event)
        else:
            self.savedParams = self.machineparams.getParameters()
            super().showEvent(event)
            
    def hideEvent(self, event):
        self.sigHide.emit()
        super().hideEvent(event)
            
    def resetParameters(self):
        self.machineparams.setValues(self.savedParams)
        #self.machineparams._onAnyValueChanged()
            
    def onCancel(self):
        self.resetParameters()
        self.hide()
    
        
class QCrystalParameterDialog(qt.QDialog):
    sigHide = qt.pyqtSignal()

    def __init__(self,crystalparams,parent=None):
        qt.QDialog.__init__(self, parent=None)
        self.crystalparams = crystalparams
        layout = qt.QVBoxLayout()
        layout.addWidget(crystalparams)
        
        self.savedParams = None
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel | qt.QDialogButtonBox.Reset,
                                      qt.Qt.Horizontal)
        layout.addWidget(buttons)
        
        okbtn = buttons.button(qt.QDialogButtonBox.Ok)
        cancelbtn = buttons.button(qt.QDialogButtonBox.Cancel)
        resetbtn = buttons.button(qt.QDialogButtonBox.Reset)
        
        okbtn.clicked.connect(self.hide)
        cancelbtn.clicked.connect(self.onCancel)
        resetbtn.clicked.connect(self.resetParameters)
        
        self.setLayout(layout)
        
    def showEvent(self, event):
        if event.spontaneous():
            super().showEvent(event)
        else:
            try:
                self.savedParams = self.crystalparams.getCrystal(), self.crystalparams.getRefractionIndex()
            except Exception:
                self.savedParams = None
            super().showEvent(event)
            
    def hideEvent(self, event):
        self.sigHide.emit()
        super().hideEvent(event)
            
    def resetParameters(self):
        if self.savedParams is not None:
            self.crystalparams.setValues(*self.savedParams)
            #self.crystalparams._onAnyValueChanged()
            
    def onCancel(self):
        self.resetParameters()
        self.hide()

        

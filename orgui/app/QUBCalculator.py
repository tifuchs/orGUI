# -*- coding: utf-8 -*-
from io import StringIO
from silx.gui import qt

#from PyMca5.PyMcaGui import PyMca_Icons as icons
#from PyMca5.PyMcaGui.io import PyMcaFileDialogs
import numpy as np

import traceback

from datautils.xrayutils import HKLVlieg
from datautils.xrayutils import DetectorCalibration
import warnings
import configparser
import os 

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
        a,alpha,_,_ = crystal.getLatticeParameters()
        self.crystal.setLattice(a,np.rad2deg(alpha))
        self.n = n
        #self.ubCal.defaultU()
        #print("New %s,\nyou have to recalculate the UB matrix to apply the changes" % self.crystal) 
        #print(crystal)
        #print(n)

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
        self.detectorCal.setPolarization(polax,polf)
        try:
            azimy,azimx = self.detectorCal.pixelsPrimPoint([0.05],[0])[0]
            self.sigPlottableMachineParamsChanged.emit([cp,[azimx,azimy],polax])
        except Exception as e:
            # here is a bug with the init of the detector cal
            pass
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
            

            a1 = lattice.getfloat('a1')
            a2 = lattice.getfloat('a2')
            a3 = lattice.getfloat('a3')
            alpha1 = lattice.getfloat('alpha1')
            alpha2 = lattice.getfloat('alpha2')
            alpha3 = lattice.getfloat('alpha3')
            self.n = 1 - lattice.getfloat('refractionindex',0.0)
            
            self.crystal = HKLVlieg.Crystal([a1,a2,a3],[alpha1,alpha2,alpha3])
            self.ubCal = HKLVlieg.UBCalculator(self.crystal,E)
            self.ubCal.defaultU()
            self.angles = HKLVlieg.VliegAngles(self.ubCal)
            
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
            qt.QMessageBox.warning(self,"Can not parse config","Can not parse config file:\nException occured during parsing of configfile %s,\nException:\n %s" % (configfile,e))
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
        self.crystal = HKLVlieg.Crystal([3.9242,3.9242,3.9242],[90.,90.,90.])
        self.ubCal = HKLVlieg.UBCalculator(self.crystal,E)
        self.ubCal.defaultU()
        self.polaxis = 0
        self.polfactor = 0
        self.azimuth = 0
        self.detectorCal = DetectorCalibration.Detector2D_SXRD()
        
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
            self.crystalparams.setValues(self.crystal,self.n)
        self.Ueditor.setPlainText(str(self.ubCal.getU()))
            
        
        
class QCrystalParameter(qt.QSplitter):
    sigCrystalParamsChanged = qt.pyqtSignal(HKLVlieg.Crystal,float)
    def __init__(self,parent=None):
        qt.QSplitter.__init__(self, parent=None)
        self.setOrientation(qt.Qt.Vertical)
        
        crystalParamsGroup = qt.QGroupBox("Lattice parameters")
        crystalParamsLayout = qt.QGridLayout()
        
        
        
        crystalParamsLayout.addWidget(qt.QLabel("a1:"),0,0)
        crystalParamsLayout.addWidget(qt.QLabel("a2:"),1,0)
        crystalParamsLayout.addWidget(qt.QLabel("a3:"),2,0)
        self.a1box = qt.QDoubleSpinBox()
        self.a1box.setRange(0,100)
        self.a1box.setDecimals(4)
        self.a1box.setSuffix(u" \u212B")
        
        crystalParamsLayout.addWidget(self.a1box,0,1)
        
        self.a2box = qt.QDoubleSpinBox()
        self.a2box.setRange(0,100)
        self.a2box.setDecimals(4)
        self.a2box.setSuffix(u" \u212B")

        crystalParamsLayout.addWidget(self.a2box,1,1)
        
        self.a3box = qt.QDoubleSpinBox()
        self.a3box.setRange(0,100)
        self.a3box.setDecimals(4)
        self.a3box.setSuffix(u" \u212B")
        
        crystalParamsLayout.addWidget(self.a3box,2,1)

        crystalParamsLayout.addWidget(qt.QLabel("alpha1:"),0,2)
        crystalParamsLayout.addWidget(qt.QLabel("alpha2:"),1,2)
        crystalParamsLayout.addWidget(qt.QLabel("alpha3:"),2,2)
        self.alpha1box = qt.QDoubleSpinBox()
        self.alpha1box.setRange(0,180)
        self.alpha1box.setDecimals(2)
        self.alpha1box.setSuffix(" °")
        
        crystalParamsLayout.addWidget(self.alpha1box,0,3)
        
        self.alpha2box = qt.QDoubleSpinBox()
        self.alpha2box.setRange(0,180)
        self.alpha2box.setDecimals(2)
        self.alpha2box.setSuffix(" °")
       
        crystalParamsLayout.addWidget(self.alpha2box,1,3)
        
        self.alpha3box = qt.QDoubleSpinBox()
        self.alpha3box.setRange(0,180)
        self.alpha3box.setDecimals(2)
        self.alpha3box.setSuffix(" °")
        
        crystalParamsLayout.addWidget(self.alpha3box,2,3)
        
        
        crystalParamsGroup.setLayout(crystalParamsLayout)
        self.addWidget(crystalParamsGroup)
        
        
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
        
        self.addWidget(refractionindexGroup)
        self.a1box.valueChanged.connect(self._onAnyValueChanged)
        self.a2box.valueChanged.connect(self._onAnyValueChanged)
        self.a3box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha1box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha2box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha3box.valueChanged.connect(self._onAnyValueChanged)
        self.refractionIndexBox.valueChanged.connect(self._onAnyValueChanged)
    
    def setValues(self,crystal,n):
        [a1,a2,a3],alpha,_,_ = crystal.getLatticeParameters()
        [alpha1,alpha2,alpha3] = np.rad2deg(alpha)
        self.a1box.setValue(a1)
        self.a2box.setValue(a2)    
        self.a3box.setValue(a3)
        self.alpha1box.setValue(alpha1)
        self.alpha2box.setValue(alpha2)
        self.alpha3box.setValue(alpha3)
        self.refractionIndexBox.setValue((1.-n)*1e6)
        
    def getCrystal(self):
        a = np.array([self.a1box.value(),self.a2box.value(),self.a3box.value()])
        alpha = np.array([self.alpha1box.value(),self.alpha2box.value(),self.alpha3box.value()])
        if np.any(a == 0) or np.any(alpha == 0):
            raise Exception("No crystal set")
        return HKLVlieg.Crystal(a,alpha)
        
    def getRefractionIndex(self):
        return 1 - self.refractionIndexBox.value()*1e-6
        
    def _onAnyValueChanged(self):
        try:
            newCrystal = self.getCrystal()
            n = self.getRefractionIndex()
            self.sigCrystalParamsChanged.emit(newCrystal,n)
        except Exception as e:
            qt.QMessageBox.warning(self,"Can not calculate B Matrix","The B Matrix can not be calculated\nError: %s" % str(e))
                
        
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
        self.blockSignals(True)
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
        self.blockSignals(False)
        
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
        self.machineparams._onAnyValueChanged()
            
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
            self.crystalparams._onAnyValueChanged()
            
    def onCancel(self):
        self.resetParameters()
        self.hide()

        

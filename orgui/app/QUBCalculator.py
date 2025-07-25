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

from io import StringIO
from silx.gui import qt
from silx.gui import icons
import pyFAI, pyFAI.detectors
import numpy as np
from scipy.spatial import transform
try:
    import ase.io
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    

import traceback
from ..datautils.xrayutils import HKLVlieg, CTRcalc
from ..datautils.xrayutils import DetectorCalibration
from ..datautils.xrayutils import unitcells
import warnings
import configparser
import os
import copy
from enum import Enum, auto
from contextlib import contextmanager

from .ArrayTableDialog import ArrayEditWidget

from ..backend import udefaults, backends
from .. import resources
from . import qutils

from pyFAI.gui.dialog import DetectorSelectorDialog
from pyFAI.gui.widgets import GeometryTabs

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

@contextmanager
def disconnectTemporarily(signal, reciever):
    signal.disconnect(reciever)
    yield
    signal.connect(reciever)
    
class LatIndex(Enum):
    A1 = auto()
    A2 = auto()
    A3 = auto()


# reflectionhandler must implement the method getReflections

class QUBCalculator(qt.QSplitter):
    sigNewReflection = qt.pyqtSignal(dict)
    sigPlottableMachineParamsChanged = qt.pyqtSignal()
    sigReplotRequest = qt.pyqtSignal(bool)
    #sigQueryImageChange = qt.pyqtSignal(int)
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self,configfile, parent):
        qt.QSplitter.__init__(self, parent)
        self.setOrientation(qt.Qt.Vertical)
        
        self.mainGui = parent
        
        self.configdir = os.getcwd()
        #self.mainLayout = qt.QVBoxLayout()

        self.setChildrenCollapsible(False)
        
        self.reflectionWidget = qt.QToolBar()
        self.reflectionWidget.setFloatable(False)
        self.reflectionWidget.setMovable(False)
        #self.reflectionWidget.setChildrenCollapsible(False)
        self.reflectionWidget.setOrientation(qt.Qt.Horizontal)
        label = qt.QLabel("H:")
        self.reflectionWidget.addWidget(label)
        self.Hbox = qt.QDoubleSpinBox()
        self.Hbox.setRange(-100,100)
        self.Hbox.setDecimals(3)
        self.reflectionWidget.addWidget(self.Hbox)
        label = qt.QLabel("K:")
        self.reflectionWidget.addWidget(label)
        self.Kbox = qt.QDoubleSpinBox()
        self.Kbox.setRange(-100,100)
        self.Kbox.setDecimals(3)
        self.reflectionWidget.addWidget(self.Kbox)
        label = qt.QLabel("L:")
        self.reflectionWidget.addWidget(label)
        self.Lbox = qt.QDoubleSpinBox()
        self.Lbox.setRange(-100,100)
        self.Lbox.setDecimals(3)
        self.reflectionWidget.addWidget(self.Lbox)
        searchReflAct = qt.QAction(resources.getQicon('search'), "search reflection", self)
        self.reflectionWidget.addAction(searchReflAct)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        searchReflAct.setToolTip("calculate position of the reflection with given HKL")
        
        searchReflAct.triggered.connect(self._onCalcReflection)
        
        self.addWidget(self.reflectionWidget)
        
        umatrixsplitter = qt.QSplitter()
        umatrixsplitter.setOrientation(qt.Qt.Horizontal)
        umatrixsplitter.setChildrenCollapsible(False)
        
        #self.Ueditor = qt.QTextEdit("")
        #umatrixsplitter.addWidget(self.Ueditor)
        self.calUButton = qt.QPushButton("calculate U")
        self.calUButton.setToolTip("calculate orientation matrix based on the given reflections")
                
        vertCalUSplitter = qt.QSplitter()
        vertCalUSplitter.setOrientation(qt.Qt.Vertical)
        vertCalUSplitter.setChildrenCollapsible(False)
        fitUbox = qt.QGroupBox("fit options")
        fitUbox.setToolTip("only available with enough reflections")
        
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
        
        self.calMiscutButton = qt.QPushButton("calculate miscut")
        self.calMiscutButton.setToolTip("calculate miscut based on deviation from ideal orientation matrix")
        
        vertCalUSplitter.addWidget(self.calMiscutButton)
        
        
        
        
        umatrixsplitter.addWidget(vertCalUSplitter)
        
        self.addWidget(umatrixsplitter)
        

        
        
        
        self.calUButton.clicked.connect(self._onCalcU)
        self.calMiscutButton.clicked.connect(self._onCalMiscut)
        
        self.crystalparams = QCrystalParameter()
        
        
        self.crystalparams.sigCrystalParamsChanged.connect(self._onCrystalParamsChanged)
        self.machineParams = QMachineParameters()
        self.machineParams.sigMachineParamsChanged.connect(self._onMachineParamsChanged)
        self.machineParams.loadConfigButton.clicked.connect(self._onLoadConfig)

        self.uedit = QUEdit()
        self.ueditDialog = QUEditDialog(self.uedit)
        
        self.uedit.sigResetRequest.connect(self._onResetU)
        self.uedit.sigAlignRequest.connect(self._onAlignU)
        

        self.machineDialog = QMachineParametersDialog(self.machineParams)
        self.xtalDialog = QCrystalParameterDialog(self.crystalparams)
        
                
        if configfile is not None:
            if not self.readConfig(configfile):
                self.toFallbackConfig()
        else:
            self.toFallbackConfig()
        
        self.uedit.sigUChanged.connect(self._onUchanged)
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
        
        
    def calcReflection(self,hkl, axisname=None):
        if self.mainGui.fscan is not None:
            if axisname is None:
                axisname = self.mainGui.fscan.axisname
            mu, om = self.mainGui.getMuOm()
        else:
            mu = self.mu; om = 0.0
        angle_factors = np.unique((om+np.pi) // (2*np.pi)) # possible angle offset factors
        angle_factors = np.sort(angle_factors)[::-1]  # prefer solutions close to offset factor of 0, first positive elements
        srt = np.argsort(np.abs(angle_factors))
        angle_factors = angle_factors[srt] # prefer solutions close to offset factor of 0
        
        ommax = np.amax(om)
        ommin = np.amin(om)
        
        mu_cryst = HKLVlieg.crystalAngles_singleArray(mu,self.n)
        hkl = np.asarray(hkl)
        if len(hkl.shape) > 1:
            hkl = hkl.T # for anglesZmode, is a bit inconsistent
            
        if axisname == 'th':
            pos1 = self.angles.anglesZmode(hkl,mu_cryst,'in',self.chi,self.phi,mirrorx=False)
            pos2 = self.angles.anglesZmode(hkl,mu_cryst,'in',self.chi,self.phi,mirrorx=True)
        elif axisname == 'mu':
            pos1 = self.angles.anglesZmode(hkl,mu_cryst,'eq',self.chi,self.phi,mirrorx=False)
            pos2 = self.angles.anglesZmode(hkl,mu_cryst,'eq',self.chi,self.phi,mirrorx=True)
        else:
            raise ValueError("No scan axis given or no scan loaded.")

        pos1_refr = HKLVlieg.vacAngles(pos1,self.n)
        pos2_refr = HKLVlieg.vacAngles(pos2,self.n)
        
        def _adjust_omega_array(omega, ommin, ommax, angle_factors):
            """ Original function, numpy parralelized with AI:
            It finds a 'best guess' for the om range used by the experiment
            The order of the tested angle_factors determines the priority.
            I.e. it should usually start (highest priority) with zero and increase.
            
            minfactor = 0
            metric = np.inf
            for factor in angle_factors:
                if ommin <= omega1 + factor*np.pi <= ommax:
                    omega1 = omega1 + factor*np.pi
                    break
                else:
                    fmetric = min( abs(omega1 + factor*np.pi  - ommin) , abs(omega1 + factor*np.pi  - ommax))
                    if fmetric < metric:
                        minfactor = factor
                        metric = fmetric
            else:
                omega1 = omega1 + minfactor*np.pi
            """
            omega = np.atleast_1d(omega)
            angle_factors = np.atleast_1d(angle_factors)

            
            # Compute all candidate shifts: shape = (N, M)
            # N = len(omega), M = len(angle_factors)
            candidates = omega[:, None] + angle_factors[None, :] * 2 *np.pi

            # Determine which candidates lie within [ommin, ommax]
            inside = (candidates >= ommin) & (candidates <= ommax)

            # For each ω₁, find the first factor that yields an “inside” result
            first_inside_idx = np.argmax(inside, axis=1)
            has_inside = inside.any(axis=1)

            # Compute distance to the nearest bound for all candidates
            dist_to_min = np.abs(candidates - ommin)
            dist_to_max = np.abs(candidates - ommax)
            nearest_dist = np.minimum(dist_to_min, dist_to_max)

            # For those with no “inside” candidate, pick the factor with minimal distance
            closest_idx = np.argmin(nearest_dist, axis=1)

            # Choose per-row index: first_inside_idx if in‐range exists, else closest_idx
            choice_idx = np.where(has_inside, first_inside_idx, closest_idx)

            # Extract the adjusted ω₁ values
            adjusted = candidates[np.arange(omega.shape[0]), choice_idx]
            return adjusted

        if len(hkl.shape) > 1:
            alpha1, delta1, gamma1, omega1, chi1, phi1 = pos1_refr.T
            alpha2, delta2, gamma2, omega2, chi2, phi2 = pos2_refr.T
            hkl = hkl.T
            omega1 = _adjust_omega_array(omega1, ommin, ommax, angle_factors)
            omega2 = _adjust_omega_array(omega2, ommin, ommax, angle_factors)
            pos1_refr.T[3] = omega1
            pos2_refr.T[3] = omega2
            
        else:
            alpha1, delta1, gamma1, omega1, chi1, phi1 = pos1_refr
            alpha2, delta2, gamma2, omega2, chi2, phi2 = pos2_refr
            omega1 = float(np.squeeze(_adjust_omega_array(omega1, ommin, ommax, angle_factors)))
            omega2 = float(np.squeeze(_adjust_omega_array(omega2, ommin, ommax, angle_factors)))
            pos1_refr[3] = omega1
            pos2_refr[3] = omega2
                
        
        xy1 = self.detectorCal.pixelsSurfaceAngles(gamma1,delta1,alpha1)[:,::-1]
        xy2 = self.detectorCal.pixelsSurfaceAngles(gamma2,delta2,alpha2)[:,::-1]
        
        di = {
           'hkl' : hkl,
           'xy_1' : np.squeeze(xy1),
           'xy_2' : np.squeeze(xy2),
           'angles_1' : pos1_refr,
           'angles_2' : pos2_refr
        }
        return di
            
        
    def _onCalcReflection(self):
        hkl = [self.Hbox.value(),self.Kbox.value(),self.Lbox.value()]
        try:
            refl = self.calcReflection(hkl)
        except Exception as e:
            qutils.warning_detailed_message(self, "Cannot calculate reflection", 
                                            "Cannot calculate reflection:\n%s" % e,
                                            traceback.format_exc())
            return
        self.sigNewReflection.emit(refl)
        
    def _onUchanged(self, U):
        self.ubCal.setU(U)
        self.sigReplotRequest.emit(True)
        
    def _onResetU(self, func):
        func(self.ubCal)
        self.uedit.setU(self.ubCal.getU())
        self.sigReplotRequest.emit(True)

    def _onAlignU(self, ddict):
        angles = ddict['angles']  # ['alpha', 'chi', 'phi', 'theta']
        pos = [angles[0], None, None, angles[3], angles[1], angles[2]]
        if ddict['frame'] == 'surface':
            self.ubCal.alignU_alpha(ddict['hkl'], pos, ddict['xyz'])
        elif ddict['frame'] == 'lab':
            self.ubCal.alignU_lab(ddict['hkl'], pos, ddict['xyz'])
        self.uedit.setU(self.ubCal.getU())
        self.sigReplotRequest.emit(True)
        
    def _onCrystalParamsChanged(self,crystal,n):
        #a,alpha,_,_ = crystal.getLatticeParameters()
        #self.crystal.setLattice(a,np.rad2deg(alpha))
        
        self.crystal = crystal
        self.n = n
        self.ubCal.setLattice(self.crystal)
        #self.ubCal.defaultU()
        self.sigReplotRequest.emit(True)
    
    def _onMachineParamsChanged(self,params):
        diffrac = params['diffractometer']
        self.mu = diffrac['mu']
        self.chi = diffrac['chi']
        self.phi = diffrac['phi']
        self.ubCal.setEnergy(params['source']['E'])
        detCal = params['SXRD_geometry']
        self.detectorCal.set_config(detCal.get_config())
        
        azim = detCal.getAzimuthalReference()
        polax, polf = detCal.getPolarization()
        self.detectorCal.setAzimuthalReference(azim)
        self.azimuth = azim
        self.detectorCal.setPolarization(polax,polf)
        self.crystal.setEnergy(params['source']['E']*1e3)
        
        angles_u = self.uedit.cached_angles

        self.uedit.setAngles(self.mu, self.chi, self.phi, angles_u[-1])
        
        try:
            self.sigPlottableMachineParamsChanged.emit()
            self.sigReplotRequest.emit(True)
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
            
            azimuth = np.deg2rad(diffrac.getfloat('azimuthal_reference',0))
            polaxis = np.deg2rad(diffrac.getfloat('polarization_axis',0))
            polfactor = diffrac.getfloat('polarization_factor',0)
            
            sdd = machine.getfloat('SDD',0.729) #m
            E =  machine.getfloat('E',78.0) #keV
            pixelsize = machine.getfloat('pixelsize',172e-6) #m
            cpx = machine.getfloat('cpx',731)
            cpy = machine.getfloat('cpy',1587)
            cp = [cpx,cpy]
            det_sizex =  machine.getfloat('sizex',3000)
            det_sizey =  machine.getfloat('sizey',3000)
            det_shape = (det_sizey, det_sizex)
            
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
                #print("Fallback lattice vectors")
                
            
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
                    try:
                        if os.path.isabs(lattice['crystal']):
                            xtalpath = lattice['crystal']
                        else:
                            p = os.path.abspath(configfile)
                            xtalpath = os.path.join(os.path.dirname(p), lattice['crystal'])
                        self.crystalparams.loadUnitCell(xtalpath)
                    except Exception:
                        qt.QMessageBox.warning(self,"Did not find crystal","Can not find crystal <%s> \nException occured during read of configfile %s,\nException:\n%s" % (lattice['crystal'],traceback.format_exc()))
                else:
                    self.crystalparams.crystalComboBox.setCurrentIndex(idx)
                    self.crystalparams.onSwitchCrystal(idx)
                    #self.crystal = self.crystalparams.getCrystal()

            if latticeoverride:
                self.crystal.setLattice([a1,a2,a3],[alpha1,alpha2,alpha3])
            
            self.detectorCal = DetectorCalibration.Detector2D_SXRD()
            if 'poni' in machine:
                if machine['poni']:
                    if os.path.isabs(machine['poni']):
                        ponipath = machine['poni']
                    else:
                        p = os.path.abspath(configfile)
                        ponipath = os.path.join(os.path.dirname(p), machine['poni'])
                    self.detectorCal.load(ponipath)
                    self.ubCal.setLambda(self.detectorCal.get_wavelength()*1e10)
                    self.crystal.setEnergy(self.detectorCal.get_energy()*1e3)
                else:
                    self.detectorCal.setFit2D(sdd*1e3,cpx,cpy,pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
                    self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
                    self.detectorCal.detector.shape = det_shape # Perkin 
                    self.detectorCal.detector.max_shape = det_shape # Perkin det_shape
                    
            else:
                self.detectorCal.setFit2D(sdd*1e3,cpx,cpy,pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
                self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
                self.detectorCal.detector.shape = det_shape
                self.detectorCal.detector.max_shape = det_shape
                
            self.detectorCal.setAzimuthalReference(azimuth)
            self.detectorCal.setPolarization(polaxis,polfactor)
            
            #fit2dCal = self.detectorCal.getFit2D()
            settings = {'diffractometer' : {
                    'mu' : self.mu,
                    'phi' : self.phi,
                    'chi' : self.chi
                },
                'source' :  {
                    'E' : self.detectorCal.get_energy()
                },
                'SXRD_geometry' : self.detectorCal
            }

            self.crystalparams.setValues(self.crystal,self.n)
            self.machineParams.setValues(settings)
            
            if 'backend' in config:
                if 'file' in config['backend']:
                    if os.path.isabs(config['backend']['file']):
                        backendpath = config['backend']['file']
                    else:
                        p = os.path.abspath(configfile)
                        backendpath = os.path.join(os.path.dirname(p), config['backend']['file'])
                    self.mainGui.scanSelector.loadBackendFile(backendpath)
                elif 'beamtime' in config['backend']:
                    beamtime = config['backend']['beamtime']
                    if beamtime in backends.fscans:
                        self.mainGui.scanSelector.btid.setCurrentText(beamtime)
                        self.mainGui.scanSelector.bt_autodetect_enable.setChecked(False)
                    else:
                        raise ValueError("Cannot find beamtime %s in the list of available backends" % beamtime)
            #self.machineParams.set_detector(self.detectorCal.detector)
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
        self.detectorCal.setAzimuthalReference(np.deg2rad(90.))
        self.detectorCal.setPolarization(0.,0.)

        self.angles = HKLVlieg.VliegAngles(self.ubCal)
        settings = {'diffractometer' : {
                'mu' : self.mu,
                'phi' : self.phi,
                'chi' : self.chi
            },
            'source' :  {
                'E' : self.detectorCal.get_energy()
            },
            'SXRD_geometry' : self.detectorCal
        }
        self.crystalparams.setValues(self.crystal,self.n)
        self.machineParams.setValues(settings)
        
        
        
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
        self.sigReplotRequest.emit(False)
        #self.Ueditor.setPlainText(str(self.ubCal.getU()))
        
    def _onCalMiscut(self):
        om, chi, phi = self.angles.anglesOrientationAlpha([0,0,1], [0,0,1])
        chi, phi = float(chi), float(phi)
        chi_displ = float(np.rad2deg(chi))
        phi_displ = float(np.rad2deg(phi))
        btn = qt.QMessageBox.question(self,"Miscut calcuated","The miscut is along\nphi = %.4f\nchi = %.4f\nReset chi and phi to negative these values?" % (phi_displ,chi_displ),qt.QMessageBox.Yes | qt.QMessageBox.No, qt.QMessageBox.No)
        if btn == qt.QMessageBox.Yes:
            self.chi = -chi
            self.phi = -phi
            angles_u = self.uedit.cached_angles
            self.uedit.setAngles(angles_u[0], -chi, -phi, angles_u[-1])

        
class QCrystalParameter(qt.QWidget):
    sigCrystalParamsChanged = qt.pyqtSignal(HKLVlieg.Lattice,float)
    def __init__(self,parent=None):
        qt.QWidget.__init__(self, parent=None)
        
        mainLayout = qt.QVBoxLayout()
        
        self.toolbar = qt.QToolBar("Crystal parameter tools")
        self.toolbar.setMovable(False)
        self.toolbar.setFloatable(False)
        
        
        self.loadxtalAct = qt.QAction(icons.getQIcon('document-open'), "Load unit cell")
        self.toolbar.addAction(self.loadxtalAct)
        toolbarbtn = self.toolbar.widgetForAction(self.loadxtalAct)
        toolbarbtn.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
        self.loadxtalAct.triggered.connect(self.onLoadXtal)
        
        self.showxtalAct = qt.QAction(resources.getQicon("lattice-view"), "show unit cell")
        self.toolbar.addAction(self.showxtalAct)
        toolbarbtn_2 = self.toolbar.widgetForAction(self.showxtalAct)
        toolbarbtn_2.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
        self.showxtalAct.triggered.connect(self.onShowXtal)
        
        spacer_widget = qt.QWidget()
        spacer_widget.setSizePolicy(qt.QSizePolicy.Expanding, qt.QSizePolicy.Preferred)
        self.toolbar.addWidget(spacer_widget)
        
        self.latticelinkgrp = qt.QActionGroup(self)
        self.latticelinkgrp.setExclusive(True)
        self.notLinkedAct = self.latticelinkgrp.addAction(resources.getQicon("lattice-no-link"), "scale a1, a2, a3 individually")
        self.horizLinkedAct = self.latticelinkgrp.addAction(resources.getQicon("lattice-horizontal-link"), "scale a1, a2 with a common factor")
        self.allLinkedAct = self.latticelinkgrp.addAction(resources.getQicon("lattice-all-link"), "scale a1, a2, a3 with a common factor")
        
        self.notLinkedAct.setCheckable(True)
        self.horizLinkedAct.setCheckable(True)
        self.allLinkedAct.setCheckable(True)

        self.link_menu = qt.QMenu()
        self.link_menu.addAction(self.notLinkedAct)
        self.link_menu.addAction(self.horizLinkedAct)
        self.link_menu.addAction(self.allLinkedAct)
        
        self.link_btn = qt.QToolButton()
        self.link_btn.setIcon(resources.getQicon("alpha"))
        self.link_btn.setToolTip("Set linking of lattice parameters")
        self.link_btn.setPopupMode(qt.QToolButton.InstantPopup)
        self.link_btn.setMenu(self.link_menu)
        
        self.notLinkedAct.triggered.connect(self._onLinkLatticeChanged)
        self.horizLinkedAct.triggered.connect(self._onLinkLatticeChanged)
        self.allLinkedAct.triggered.connect(self._onLinkLatticeChanged)
        
        self.notLinkedAct.trigger()
        
        self.toolbar.addWidget(self.link_btn)
        
        
        #crystalParamsLayout.addWidget(loadxtalbtn,1,1)
        
        
        
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
            
        
        self.crystalComboBox.activated.connect(self.onSwitchCrystal)
        
        
        crystalParamsGroup.setLayout(crystalParamsLayout)
        
        mainLayout.addWidget(self.toolbar)
        mainLayout.addWidget(crystalParamsGroup)
        mainLayout.addWidget(latticeParamsGroup)
        mainLayout.addWidget(refractionindexGroup)
        self.setLayout(mainLayout)
        
        
        self.a1box.valueChanged.connect(lambda : self._onLatticeParamsChanged(LatIndex.A1))
        self.a2box.valueChanged.connect(lambda : self._onLatticeParamsChanged(LatIndex.A2))
        self.a3box.valueChanged.connect(lambda : self._onLatticeParamsChanged(LatIndex.A3))
        self.alpha1box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha2box.valueChanged.connect(self._onAnyValueChanged)
        self.alpha3box.valueChanged.connect(self._onAnyValueChanged)
        self.refractionIndexBox.valueChanged.connect(self._onAnyValueChanged)
        
    def _onLinkLatticeChanged(self, checked):
        self.link_btn.setIcon(self.latticelinkgrp.checkedAction().icon())

    def onShowXtal(self):
        repr_xtal = str(self.getCrystal())
        qutils.information_detailed_message(self, 
                                            "Unit cell parameters", 
                                            "Show the detailed text to display the coordinates of the atoms in the chosen unit cell.",
                                            repr_xtal)
        
    def _onLatticeParamsChanged(self, which):
        act = self.latticelinkgrp.checkedAction()
        if act is self.notLinkedAct:
            self._onAnyValueChanged()
            return
        elif act is self.horizLinkedAct:
            if which == LatIndex.A1:
                ratio = self.a1box.value() / self._uc.a[0]
                with blockSignals([self.a2box]):
                    self.a2box.setValue(self._uc.a[1] * ratio)
            elif which == LatIndex.A2:
                ratio = self.a2box.value() / self._uc.a[1]
                with blockSignals([self.a1box]):
                    self.a1box.setValue(self._uc.a[0] * ratio)
            self._onAnyValueChanged()
            return
        elif act is self.allLinkedAct:
            if which == LatIndex.A1:
                ratio = self.a1box.value() / self._uc.a[0]
                with blockSignals([self.a2box, self.a3box]):
                    self.a2box.setValue(self._uc.a[1] * ratio)
                    self.a3box.setValue(self._uc.a[2] * ratio)
            elif which == LatIndex.A2:
                ratio = self.a2box.value() / self._uc.a[1]
                with blockSignals([self.a1box, self.a3box]):
                    self.a1box.setValue(self._uc.a[0] * ratio)
                    self.a3box.setValue(self._uc.a[2] * ratio)
            elif which == LatIndex.A3:
                ratio = self.a3box.value() / self._uc.a[2]
                with blockSignals([self.a1box, self.a2box]):
                    self.a1box.setValue(self._uc.a[0] * ratio)
                    self.a2box.setValue(self._uc.a[1] * ratio)
            self._onAnyValueChanged()
            return
        else:
            raise ValueError("Invalid lattice parameter changed: %s" % which)
            
    def loadUnitCell(self, filename):
        ext = os.path.splitext(filename)[1]
        if ext in ['.xtal', '.h5', '.xpr']:
            xtal = CTRcalc.SXRDCrystal.fromFile(filename)
            uc_bulk = xtal['bulk']
        else:
            uc_bulk = CTRcalc.UnitCell.fromFile(filename)
        uc_name = os.path.splitext(os.path.basename(filename) )[0]
        self.crystalComboBox.addItem(uc_name, uc_bulk)
        idx = self.crystalComboBox.findText(uc_name)
        self.crystalComboBox.setCurrentIndex(idx)
        self.onSwitchCrystal(idx)
        
    def onLoadXtal(self):
        fileTypeDict = {'ANA ROD files (*.bul *.sur)': '.bul', 'Crystal Files (*.xtal *.xpr)': '.xtal .xpr'}
        if HAS_ASE:
            fmt = ase.io.formats.ioformats
            r_ext = list(filter(lambda d : fmt[d].can_read, fmt))
            r_ext = [ '.' + r for r in r_ext]
            ase_extensions = ' '.join(r_ext)
            fileTypeDict['ASE supported files (*.cif *.vasp *.xyz *.abinit-in *.abinit-out *.espresso-in *.espresso-out *.mol *)'] = ase_extensions
        fileTypeDict['All files (*)'] = ''

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
        try:
            self.loadUnitCell(filename)
        except Exception:
            qt.QMessageBox.critical(self,"Cannot open unit cell file", "Cannot open:\n%s" % traceback.format_exc())
            return
        
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
    sigResetRequest = qt.pyqtSignal(object)
    sigAlignRequest = qt.pyqtSignal(object)
    
    def __init__(self,parent=None):
        qt.QWidget.__init__(self, parent=None)
        
        alignGroup = qt.QGroupBox("Align hkl onto xyz")
        orientationLayout = qt.QGridLayout()
        
        self.hkl = []
        for i, index in enumerate(['H', 'K', 'L']):
            orientationLayout.addWidget(qt.QLabel("%s:" % index),i,0)
            milleredit = qt.QDoubleSpinBox()
            milleredit.setRange(-1000,1000)
            milleredit.setDecimals(4)
            orientationLayout.addWidget(milleredit,i,1)
            milleredit.setValue(1. if i == 2 else 0.)
            self.hkl.append(milleredit)
            
        self.override_angles = qt.QCheckBox("Override angles")
        self.override_angles.toggled.connect(self._onOverride)
        orientationLayout.addWidget(self.override_angles, 3, 0, 1, 2)
        
        self.cached_angles = [0. ,0., 0., 0.] 
        self.angles = []
        for i, index in enumerate(['alpha', 'chi', 'phi', 'theta']):
            orientationLayout.addWidget(qt.QLabel("%s:" % index),i,2)
            edit = qt.QDoubleSpinBox()
            edit.setRange(-360,360)
            edit.setDecimals(4)
            edit.setSuffix(u" °")
            edit.setValue(0.)
            edit.setEnabled(False)
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
        self.alignBtn.clicked.connect(self.onAlignU)
        alignbtnlayout.addWidget(self.alignBtn)

        
        editLayout = qt.QHBoxLayout()
        
        editLayout.addWidget(alignGroup)
        editLayout.addLayout(alignbtnlayout)
        
        
        
        
        bottomLayout = qt.QHBoxLayout()
        
        uGroup = qt.QGroupBox("Orientation matrix")
        self.uview = ArrayEditWidget(True, 1, False)
        self.uview.model.dataChanged.connect(self.onUChanged)
        self.uview.sigDataLoaded.connect(self.onUChanged)
        la = qt.QVBoxLayout()
        la.addWidget(self.uview)
        uGroup.setLayout(la)
        
        bottomLayout.addWidget(uGroup)
        
        defaultGroup = qt.QGroupBox("Default geometries")
        
        self.uDefaults = qt.QComboBox()
        for i, geometry in enumerate(udefaults.u_defaults):
            self.uDefaults.addItem(geometry, udefaults.u_defaults[geometry])
            self.uDefaults.setItemData(i, udefaults.u_defaults[geometry].__doc__, qt.Qt.ToolTipRole)
        udef = qt.QVBoxLayout()
        udef.addWidget(self.uDefaults)
        resetUbtn = qt.QPushButton("reset")
        udef.addWidget(resetUbtn)
        resetUbtn.clicked.connect(self.onResetU)
        
        defaultGroup.setLayout(udef)

        bottomLayout.addWidget(defaultGroup)
        
        rotateUGroup = qt.QGroupBox("Manual rotate U")
        rotateUlayout = qt.QHBoxLayout()
        
        self.euler_xyz = []
        for i, index in enumerate(['x', 'y', 'z']):
            rotateUlayout.addWidget(qt.QLabel("%s:" % index))
            edit = qt.QDoubleSpinBox()
            edit.setRange(-360,360)
            edit.setDecimals(5)
            edit.setSingleStep(0.01)
            rotateUlayout.addWidget(edit)
            edit.setValue(0)
            edit.valueChanged.connect(self.onEulerChanged)
            self.euler_xyz.append(edit)
        rotateUGroup.setLayout(rotateUlayout)
        
        mainLayout = qt.QVBoxLayout()
        mainLayout.addLayout(editLayout)
        mainLayout.addWidget(rotateUGroup)
        mainLayout.addLayout(bottomLayout)
        
        self.setLayout(mainLayout)
        
    def _onOverride(self, override):
        if override:
            for edit, val in zip(self.angles, self.cached_angles):
                edit.setEnabled(True)
        else:
            for i, (edit, val) in enumerate(zip(self.angles, self.cached_angles)):
                if i == 3:
                    edit.setValue(-np.rad2deg(val))
                else:
                    edit.setValue(np.rad2deg(val))
                edit.setEnabled(False)
    
    def onAlignU(self):
        hkl = np.array([edit.value() for edit in self.hkl])
        angles = self.getAngles()
        xyz = np.array([edit.value() for edit in self.xyz])
        if self.alphaFrame.isChecked():
            frame = 'surface'
        elif self.labFrame.isChecked():
            frame = 'lab'
        ddict = {'hkl' : hkl,
                 'xyz' : xyz,
                 'angles' : angles,
                 'frame' : frame}
        self.sigAlignRequest.emit(ddict)
        
    def onEulerChanged(self):
        xyz_angles = np.deg2rad(np.array([edit.value() for edit in self.euler_xyz]))
        newU = transform.Rotation.from_euler("xyz", xyz_angles).as_matrix()
        self.setU(newU)
        self.onUChanged()
        
    def setU(self, U):
        with blockSignals(self.euler_xyz):
            self.uview.setArrayData(U, None, True, True)
            new_euler_xyz_val = np.rad2deg(transform.Rotation.from_matrix(U).as_euler('xyz'))
            [edit.setValue(v) for edit, v in zip(self.euler_xyz, new_euler_xyz_val)]
        
        
    def getU(self):
        return self.uview.getData()
        
    def setAngles(self, alpha, chi, phi, omega):
        self.cached_angles = [alpha, chi, phi, omega] 
        if not self.override_angles.isChecked():
            for i, (edit, val) in enumerate(zip(self.angles, self.cached_angles)):
                if i == 3:
                    edit.setValue(-np.rad2deg(val))
                else:
                    edit.setValue(np.rad2deg(val))
                    
    def getAngles(self):
        ang = []
        for i, edit in enumerate(self.angles):
            if i == 3:
                ang.append(-np.deg2rad(edit.value()))
            else:
                ang.append(np.deg2rad(edit.value()))
        return np.array(ang)
    
    def onResetU(self):
        self.sigResetRequest.emit(self.uDefaults.currentData())
        
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
    sigMachineParamsChanged = qt.pyqtSignal(dict)
    def __init__(self,parent=None):
        qt.QWidget.__init__(self, parent=None)
        
        self._detector = None
        self._detectorDialog = DetectorSelectorDialog.DetectorSelectorDialog(self)
        
        #[E,mu,sdd,pixsize,cp,chi,phi] = params
        
        vertical_layout = qt.QVBoxLayout()        
        
        diffractometer_box = qt.QGroupBox("Diffractometer settings", self)
        source_box = qt.QGroupBox("Source settings", self)
        
        sourceLayout = qt.QGridLayout()
        sourceLayout.addWidget(qt.QLabel("Energy:"),0,0)
        sourceLayout.addWidget(qt.QLabel("Wavelength:"),1,0)
        
        self.Ebox = qt.QDoubleSpinBox()
        self.Ebox.setRange(0.01,1000)
        self.Ebox.setDecimals(4)
        self.Ebox.setSuffix(u" keV")
        
        sourceLayout.addWidget(self.Ebox,0,1)
        
        self.wavelengthBox = qt.QDoubleSpinBox()
        self.wavelengthBox.setRange(0.001,1000)
        self.wavelengthBox.setDecimals(4)
        self.wavelengthBox.setSuffix(u" \u212B")
        
        sourceLayout.addWidget(self.wavelengthBox,1,1)
        
        source_box.setLayout(sourceLayout)
        
        mainLayout = qt.QGridLayout()
        
        diffractometerLayout = qt.QGridLayout()
        
        
        #diffractometerLayout.addWidget(qt.QLabel("SDD:"),1,0)
        diffractometerLayout.addWidget(qt.QLabel("mu:"),0,0)
        diffractometerLayout.addWidget(qt.QLabel("chi:"),1,0)
        diffractometerLayout.addWidget(qt.QLabel("phi:"),2,0)
        diffractometerLayout.addWidget(qt.QLabel("azimuth:"),3,0)
        diffractometerLayout.addWidget(qt.QLabel("Pol factor:"),4,0)
        
        
        diffractometerLayout.addWidget(qt.QLabel("polarization axis:"),5,0)
        
        
        self.mubox = qt.QDoubleSpinBox()
        self.mubox.setRange(0,90)
        self.mubox.setSingleStep(0.01)
        self.mubox.setDecimals(4)
        self.mubox.setSuffix(u" °")
        
        diffractometerLayout.addWidget(self.mubox,0,1)
        
        self.chibox = qt.QDoubleSpinBox()
        self.chibox.setRange(0,90)
        self.chibox.setDecimals(4)
        self.chibox.setSuffix(u" °")
        
        diffractometerLayout.addWidget(self.chibox,1,1)
        
        self.phibox = qt.QDoubleSpinBox()
        self.phibox.setRange(0,90)
        self.phibox.setDecimals(4)
        self.phibox.setSuffix(u" °")
        
        diffractometerLayout.addWidget(self.phibox,2,1)
        
        
        
        self.azimbox = qt.QDoubleSpinBox()
        self.azimbox.setRange(0,360)
        self.azimbox.setDecimals(4)
        self.azimbox.setSuffix(u" °")
        
        diffractometerLayout.addWidget(self.azimbox,3,1)        

        self.polfbox = qt.QDoubleSpinBox()
        self.polfbox.setRange(-1,1)
        self.polfbox.setSingleStep(0.1)
        self.polfbox.setDecimals(4)
        #self.polfbox.setSuffix(u" °")
        
        diffractometerLayout.addWidget(self.polfbox,4,1)
        
        self.polaxbox = qt.QDoubleSpinBox()
        self.polaxbox.setRange(0,360)
        self.polaxbox.setDecimals(4)
        self.polaxbox.setSuffix(u" °")
        
        diffractometerLayout.addWidget(self.polaxbox,5,1)
        
        #self.setValues(params)
        
        self.Ebox.valueChanged.connect(lambda val: self._onSourceChanged({'E' : val}) )
        self.wavelengthBox.valueChanged.connect(lambda val: self._onSourceChanged({'wavelength' : val}) )
        
        self.mubox.valueChanged.connect(self._onAnyValueChanged)
        self.chibox.valueChanged.connect(self._onAnyValueChanged)
        self.phibox.valueChanged.connect(self._onAnyValueChanged)
        self.polaxbox.valueChanged.connect(self._onAnyValueChanged)
        self.azimbox.valueChanged.connect(self._onAnyValueChanged)
        self.polfbox.valueChanged.connect(self._onAnyValueChanged)
        
        diffractometer_box.setLayout(diffractometerLayout)
        
        detector_box = qt.QGroupBox("Detector", self)
        
        self._selectDetectorBtn = qt.QPushButton("...")
        width = self._selectDetectorBtn.fontMetrics().boundingRect("  ...  ").width() + 7
        height = self._selectDetectorBtn.fontMetrics().boundingRect("  M  ").height() + 7
        self._selectDetectorBtn.setMaximumWidth(width)
        self._selectDetectorBtn.setMaximumHeight(height)
        
        self._selectDetectorBtn.clicked.connect(self._onSelectDetector)
        
        self._detectorLabel = qt.QLabel("No detector")
        self._detectorSize = qt.QLabel("")
        self._detectorPixelSize = qt.QLabel("")
        self._detectorFileDescription = qt.QLabel("")
        self._detectorFileDescriptionTitle = qt.QLabel("")
        self._detectorSizeUnit = qt.QLabel("px")
        self._detectorSizeUnit.setVisible(False)
        
        self._detectorLabel.setStyleSheet("QLabel { color: red }")
        self._detectorFileDescription.setVisible(False)
        self._detectorFileDescriptionTitle.setVisible(False)
        
        detector_panel_layout = qt.QGridLayout()
        
        detector_panel_layout.addWidget(qt.QLabel("Name:"),0,0)
        detector_panel_layout.addWidget(self._detectorLabel,0,1)
        detector_panel_layout.addWidget(self._selectDetectorBtn, 0, 2)
        
        detector_panel_layout.addWidget(qt.QLabel("Size (hxw):"),1,0)
        detector_panel_layout.addWidget(self._detectorSize, 1, 1)
        detector_panel_layout.addWidget(self._detectorSizeUnit,1,2)
        
        detector_panel_layout.addWidget(qt.QLabel("Pixel Size (hxw):"),2,0)
        detector_panel_layout.addWidget(self._detectorPixelSize, 2, 1)
        detector_panel_layout.addWidget(qt.QLabel(u"\u03BCm"), 2, 2)

        detector_panel_layout.addWidget(self._detectorFileDescription, 3, 0)
        detector_panel_layout.addWidget(self._detectorFileDescriptionTitle, 4, 0)
        
        
        detector_box.setLayout(detector_panel_layout)
        
        geometry_box = qt.QGroupBox("Detector Geometry", self)
        geometry_box_layout = qt.QVBoxLayout()
        
       
        self.geometryTabs = GeometryTabs.GeometryTabs()
        self.geometryTabs.geometryModel().changed.connect(self._onAnyValueChanged)
        geometry_box_layout.addWidget(self.geometryTabs)
        geometry_box.setLayout(geometry_box_layout)
        
        
        mainLayout.addWidget(detector_box, 0,0)
        mainLayout.addWidget(diffractometer_box, 1, 0)
        mainLayout.addWidget(source_box, 0, 1)
        mainLayout.addWidget(geometry_box, 1, 1)
        
        
        vertical_layout.addLayout(mainLayout)
        
        horizonal_buttons_layout = qt.QHBoxLayout()
        self.loadConfigButton = qt.QPushButton("load config",self) 
        self.loadConfigButton.setToolTip("load machine and crystal configuration from configfile,\naccepts poni file from pyFAI")
        horizonal_buttons_layout.addWidget(self.loadConfigButton)
        
        self.loadPoniButton = qt.QPushButton("load poni",self) 
        self.loadPoniButton.clicked.connect(self._onLoadPoni)
        horizonal_buttons_layout.addWidget(self.loadPoniButton)
        
        vertical_layout.addLayout(horizonal_buttons_layout)
        
        self.setLayout(vertical_layout)
        
    def _onLoadPoni(self):
        f,_ = qt.QFileDialog.getOpenFileName(self,"Open PyFAI calibration file","","PyFAI poni file (*.poni), All files (*)")
        if f != '':
            try:
                az = pyFAI.load(f)
                #self._detectorDialog.selectDetector(az.detector)
                self.set_Xray_source({'wavelength' : az.get_wavelength()*1e10})
                model = self.geometryTabs.geometryModel()
                model.lockSignals()
                model.distance().setValue(az.get_dist())
                model.poni1().setValue(az.get_poni1())
                model.poni2().setValue(az.get_poni2())
                model.rotation1().setValue(az.get_rot1())
                model.rotation2().setValue(az.get_rot2())
                model.rotation3().setValue(az.get_rot3())
                #model.wavelength().setValue(az.get_wavelength())
                model.unlockSignals()
                self.set_detector(az.detector)
                self._onAnyValueChanged()
                
            except Exception:
                qt.QMessageBox.warning(self,"Cannot load calibration","Cannot load poni file:\n%s" % str(traceback.format_exc()))
    
    def set_Xray_source(self, source_config):
        if 'E' in source_config:
            E = source_config['E'] # keV
            wavelength = 12.398419843320026 / E # Angstrom
        elif 'wavelength' in source_config:
            wavelength = source_config['wavelength'] # Angstrom
            E = 12.398419843320026 / wavelength # keV
        else:
            raise ValueError("misformed X-ray source settings: requries either E or wavelength, is: %s" % source_config)
        
        with blockSignals([self.Ebox, self.wavelengthBox]):
            self.Ebox.setValue(E)
            self.wavelengthBox.setValue(wavelength)
        model = self.geometryTabs.geometryModel()
        with disconnectTemporarily(self.geometryTabs.geometryModel().changed, self._onAnyValueChanged):
            model.wavelength().setValue(wavelength*1e-10)
        
    
    def _onSourceChanged(self, source_config):
        self.set_Xray_source(source_config)
        self._onAnyValueChanged()


    def set_detector(self, detector):
        with disconnectTemporarily(self.geometryTabs.geometryModel().changed, self._onAnyValueChanged):
            self.geometryTabs.setDetector(detector)
        self._detectorSizeUnit.setVisible(detector is not None)
        if detector is None:
            self._detectorLabel.setStyleSheet("QLabel { color: red }")
            self._detectorLabel.setText("No detector")
            self._detectorSize.setText("")
            self._detectorPixelSize.setText("")
            self._detectorFileDescription.setVisible(False)
            self._detectorFileDescriptionTitle.setVisible(False)
        else:
            self._detectorLabel.setStyleSheet("QLabel { }")
            self._detectorLabel.setText(detector.name)
            text = [str(s) for s in detector.max_shape]
            text = u" × ".join(text)
            self._detectorSize.setText(text)
            try:
                text = ["%0.1f" % (s * 10 ** 6) for s in [detector.pixel1, detector.pixel2]]
                text = u" × ".join(text)
            except Exception as e:
                # Is heterogeneous detectors have pixel size?
                #_logger.debug(e, exc_info=True)
                text = "N.A."
            self._detectorPixelSize.setText(text)

            if detector.HAVE_TAPER or detector.__class__ == pyFAI.detectors.Detector:
                fileDescription = detector.get_splineFile()
            elif isinstance(detector, pyFAI.detectors.NexusDetector):
                fileDescription = detector.filename
            else:
                fileDescription = None
            if fileDescription is not None:
                fileDescription = fileDescription.strip()
            if fileDescription == "":
                fileDescription = None

            self._detectorFileDescription.setVisible(fileDescription is not None)
            self._detectorFileDescriptionTitle.setVisible(fileDescription is not None)
            self._detectorFileDescription.setText(fileDescription if fileDescription else "")
        self._detector = detector
        
    def get_SXRD_geometry(self):
        detectorCal = DetectorCalibration.Detector2D_SXRD()
        model = self.geometryTabs.geometryModel()
        dist = model.distance().value()
        poni1 = model.poni1().value()
        poni2 = model.poni2().value()
        rot1 = model.rotation1().value()
        rot2 = model.rotation2().value()
        rot3 = model.rotation3().value()
        wavelength = model.wavelength().value()
        detectorCal.setPyFAI(dist=dist,
                          poni1=poni1,
                          poni2=poni2,
                          rot1=rot1,
                          rot2=rot2,
                          rot3=rot3,
                          detector=self.get_detector(),
                          wavelength=wavelength)
        azim = np.deg2rad(self.azimbox.value())
        polax = np.deg2rad(self.polaxbox.value())
        polf = self.polfbox.value()
        detectorCal.setAzimuthalReference(azim)
        detectorCal.setPolarization(polax,polf)
        return detectorCal
        
    def get_detector(self):
        return self._detector
        
    def _onSelectDetector(self):
        detector = self.get_detector()
        self._detectorDialog.selectDetector(detector)
        if self._detectorDialog.exec():
            newdetector = self._detectorDialog.selectedDetector()
            self.set_detector(newdetector)
            self._onAnyValueChanged()
            
        
    def setValues(self,params):
        signList = [self.mubox, self.chibox,
                     self.phibox, self.polaxbox,
                    self.azimbox, self.polfbox]
                    
        self.set_Xray_source(params['source'])
        detCal = params['SXRD_geometry']
        diffrac = params['diffractometer']
                    
        with blockSignals(signList):
            self.mubox.setValue(np.rad2deg(diffrac['mu']))
            self.chibox.setValue(np.rad2deg(diffrac['chi']))
            self.phibox.setValue(np.rad2deg(diffrac['phi']))
            self.azimbox.setValue(np.rad2deg(detCal.getAzimuthalReference()))
            polax, polfac = detCal.getPolarization()
            self.polaxbox.setValue(np.rad2deg(polax))
            self.polfbox.setValue(polfac)
        model = self.geometryTabs.geometryModel()
        with disconnectTemporarily(self.geometryTabs.geometryModel().changed, self._onAnyValueChanged):
            with model.lockContext():
                model.distance().setValue(detCal.dist)
                model.poni1().setValue(detCal.poni1)
                model.poni2().setValue(detCal.poni2)
                model.rotation1().setValue(detCal.rot1)
                model.rotation2().setValue(detCal.rot2)
                model.rotation3().setValue(detCal.rot3)
                model.wavelength().setValue(detCal.wavelength)


        self.set_detector(detCal.detector)
        #self._onAnyValueChanged()
        
    def getParameters(self):
        E = self.Ebox.value()
        mu = np.deg2rad(self.mubox.value())
        chi = np.deg2rad(self.chibox.value())
        phi = np.deg2rad(self.phibox.value())
        settings = {'diffractometer' : {
                'mu' : mu,
                'phi' : phi,
                'chi' : chi
            },
            'source' :  {
                'E' : E
            },
            'SXRD_geometry' : self.get_SXRD_geometry()
        }
        return settings
        
    def _onAnyValueChanged(self):
        self.sigMachineParamsChanged.emit(self.getParameters())
        
        
        
class QMachineParametersDialog(qt.QDialog):
    sigHide = qt.pyqtSignal()

    def __init__(self,machineparams,parent=None):
        qt.QDialog.__init__(self, parent=None)
        self.machineparams = machineparams
        layout = qt.QVBoxLayout()
        layout.addWidget(machineparams)
        
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

        

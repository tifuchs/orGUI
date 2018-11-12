
from io import StringIO
from PyMca5.PyMcaGui import PyMcaQt as qt

from PyMca5.PyMcaGui import PyMca_Icons as icons
from PyMca5.PyMcaGui.io import PyMcaFileDialogs
import numpy as np

from datautils.xrayutils import HKLVlieg
from datautils.xrayutils import DetectorCalibration

# reflectionhandler must implement the method getReflections

class QUBCalculator(qt.QTabWidget):
    sigNewReflection = qt.pyqtSignal(list)
    #sigQueryImageChange = qt.pyqtSignal(int)
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self, parent=None):
        qt.QTabWidget.__init__(self, parent=None)
        sdd = 0.729
        E = 68.
        pixelsize = 172e-6
        cp = [731.0,1587.856]
        #self.reflections = reflectionhandler
        self.crystal = HKLVlieg.Crystal([3.9242,3.9242,3.9242],[90.,90.,90.])
        self.detectorCal = DetectorCalibration.DetectorCalibration(E,self.crystal,pixelsize)
        self.detectorCal.setCalibration(cp,sdd)
        self.ubCal = HKLVlieg.UBCalculator(self.crystal,E)
        self.ubCal.defaultU()
        self.mu = np.deg2rad(0.05)
        self.chi = 0.
        self.phi = 0.
        self.angles = HKLVlieg.VliegAngles(self.ubCal)
        self.n = 1 - 1.1415e-06
        
        paramlist = [E,self.mu,sdd,pixelsize,cp,0.,0.]
        
        #self.setOrientation(qt.Qt.Vertical)
        
        
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
        
        self.Ueditor = qt.QTextEdit(str(self.ubCal.getU()))
        umatrixsplitter.addWidget(self.Ueditor)
        self.calUButton = qt.QPushButton("calculate U")
        self.calUButton.setToolTip("calculate orientation matrix based on the given reflections")
        umatrixsplitter.addWidget(self.calUButton)
        umatrixwidget.addWidget(umatrixsplitter)
        
        
        self.calUButton.clicked.connect(self._onCalcU)
        self.addTab(umatrixwidget,"U Matrix")
        
        self.crystalparams = QCrystalParameter(self.crystal,self.n)
        self.addTab(self.crystalparams,"Crystal")
        
        self.crystalparams.sigCrystalParamsChanged.connect(self._onCrystalParamsChanged)
        
        #paramsSplitter.setOrientation(qt.Qt.Horizontal)
        
        self.machineParams = QMachineParameters(paramlist)
        self.machineParams.sigMachineParamsChanged.connect(self._onMachineParamsChanged)
        self.addTab(self.machineParams,"Machine")
        
        
        
        
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
        gamma += self.mu
        x,y = self.detectorCal.delGamToxy(delta,gamma)
        return [hkl,x,y,omega]
        
    def _onCalcReflection(self):
        hkl = [self.Hbox.value(),self.Kbox.value(),self.Lbox.value()]
        self.sigNewReflection.emit(self.calcReflection(hkl))
        
        
    def _onCrystalParamsChanged(self,crystal,n):
        a,alpha,_,_ = crystal.getLatticeParameters()
        self.crystal.setLattice(a,alpha)
        self.n = n
        #self.ubCal.defaultU()
        print("New %s,\nyou have to recalculate the UB matrix to apply the changes" % self.crystal) 
        #print(crystal)
        #print(n)

    def _onMachineParamsChanged(self,params):
        [E,mu,sdd,pixsize,cp,chi,phi] = params
        self.mu = mu
        self.chi = chi
        self.phi = phi
        self.detectorCal.setE(E)
        self.detectorCal.setCalibration(cp,sdd)
        self.ubCal.setEnergy(E)
        # pixsize missing!
        

    def setReflectionHandler(self,refls):
        self.reflections = refls
        
    def _onCalcU(self):
        hkls,angles = self.reflections()
        if len(hkls) < 2:
            qt.QMessageBox.warning(self,"Not enough reflections","You must select at least 2 reflections to calculate a orientation matrix")
            return
        self.ubCal.setPrimaryReflection(angles[0],hkls[0])
        self.ubCal.setSecondayReflection(angles[1],hkls[1])
        self.ubCal.calculateU()
        if len(hkls) > 2:
            self.ubCal.refineU(hkls,angles)
        self.Ueditor.setPlainText(str(self.ubCal.getU()))
            
        
        
class QCrystalParameter(qt.QSplitter):
    sigCrystalParamsChanged = qt.pyqtSignal(HKLVlieg.Crystal,float)
    def __init__(self,crystal, n,parent=None):
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
        self.setValues(crystal,n)
        
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
        
        
    def _onAnyValueChanged(self):
        a = np.array([self.a1box.value(),self.a2box.value(),self.a3box.value()])
        alpha = np.array([self.alpha1box.value(),self.alpha2box.value(),self.alpha3box.value()])
        if np.any(a == 0) or np.any(alpha == 0):
            return
        try:
            newCrystal = HKLVlieg.Crystal(a,alpha)
            n = 1 - self.refractionIndexBox.value()*1e-6
            self.sigCrystalParamsChanged.emit(newCrystal,n)
        except Exception as e:
            qt.QMessageBox.warning(self,"Can not calculate B Matrix","The B Matrix can not be calculated\nError: %s" % str(e))
        
        
        
        
        
class QMachineParameters(qt.QWidget):
    sigMachineParamsChanged = qt.pyqtSignal(list)
    def __init__(self,params,parent=None):
        qt.QWidget.__init__(self, parent=None)
        
        [E,mu,sdd,pixsize,cp,chi,phi] = params
        
        mainLayout = qt.QGridLayout()
        
        mainLayout.addWidget(qt.QLabel("Energy:"),0,0)
        mainLayout.addWidget(qt.QLabel("SDD:"),1,0)
        mainLayout.addWidget(qt.QLabel("mu:"),2,0)
        mainLayout.addWidget(qt.QLabel("chi:"),3,0)
        
        mainLayout.addWidget(qt.QLabel("Centr X:"),0,2)
        mainLayout.addWidget(qt.QLabel("Centr Y:"),1,2)
        mainLayout.addWidget(qt.QLabel("pixel size:"),2,2)
        mainLayout.addWidget(qt.QLabel("phi:"),3,2)
        
        
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
        
        self.setValues(params)
        
        self.Ebox.valueChanged.connect(self._onAnyValueChanged)
        self.mubox.valueChanged.connect(self._onAnyValueChanged)
        self.SDDbox.valueChanged.connect(self._onAnyValueChanged)
        self.cpXbox.valueChanged.connect(self._onAnyValueChanged)
        self.cpYbox.valueChanged.connect(self._onAnyValueChanged)
        self.chibox.valueChanged.connect(self._onAnyValueChanged)
        self.phibox.valueChanged.connect(self._onAnyValueChanged)
        self.pixsizebox.valueChanged.connect(self._onAnyValueChanged)
        
        
        self.setLayout(mainLayout)
        
    def setValues(self,params):
        [E,mu,sdd,pixsize,cp,chi,phi] = params
        self.Ebox.setValue(E)
        self.SDDbox.setValue(sdd)
        self.mubox.setValue(np.rad2deg(mu))
        self.chibox.setValue(np.rad2deg(chi))
        self.cpXbox.setValue(cp[0])
        self.cpYbox.setValue(cp[1])
        self.pixsizebox.setValue(pixsize*1e3)
        self.phibox.setValue(np.rad2deg(phi))
        
        
        
    def _onAnyValueChanged(self):
        E = self.Ebox.value()
        mu = np.deg2rad(self.mubox.value())
        sdd = self.SDDbox.value()
        pixsize = self.pixsizebox.value()*1e-3
        cp = [self.cpXbox.value(),self.cpYbox.value()]
        chi = np.deg2rad(self.chibox.value())
        phi = np.deg2rad(self.phibox.value())
        sigList = [E,mu,sdd,pixsize,cp,chi,phi]
        self.sigMachineParamsChanged.emit(sigList)
        
        
        
        
        
        
        

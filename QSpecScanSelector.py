# -*- coding: utf-8 -*-
import sys
import os
import shutil
from PyMca5.PyMcaGui import PyMcaQt as qt

from PyMca5.PyMcaGui.io.QSpecFileWidget import QSpecFileWidget
from PyMca5.PyMcaGui.pymca import QDataSource
from PyMca5.PyMcaGui.io.QSourceSelector import QSourceSelector
from PyMca5.PyMcaGui import PyMca_Icons as icons
from PyMca5.PyMcaGui.io import PyMcaFileDialogs
from PyMca5.PyMcaGui.pymca.PyMcaNexusWidget import PyMcaNexusWidget
from PyMca5.PyMcaGui.io.hdf5 import QNexusWidget
from PyMca5.PyMcaCore import SpecFileDataSource
import warnings

from PyMca5.PyMcaCore import NexusDataSource


class QSpecScanSelector(qt.QMainWindow):
    sigScanChanged = qt.pyqtSignal(object)
    sigImagePathChanged = qt.pyqtSignal(object)
    sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self, parent=None):
        qt.QMainWindow.__init__(self, parent=None)
        
        self.mainwidget = qt.QWidget()
        self.mainLayout = qt.QVBoxLayout()
        maintab = qt.QTabWidget()
        
        
        self.specfileWidget = QSpecFileWidget(self,True)
        self.sourceSelector = QSourceSelector(self,["Spec Files (*spec)",
                                                    "Spec Files (*dat)",
                                                    "SPE Files (*SPE)",
                                                    "hdf Files (*.h5)",
                                                    "All Files (*)"])
        self.nexusWidget = PyMcaNexusWidget()
        
    
        self.specfileWidget.autoOffBox.setEnabled(False)
        self.specfileWidget.autoAddBox.setEnabled(False)
        self.specfileWidget.autoReplaceBox.setEnabled(False)
        self.specfileWidget.meshBox.setEnabled(False)
        self.specfileWidget.forceMcaBox.setEnabled(False)
        self.specfileWidget.object3DBox.setEnabled(False)
        self.specfileWidget.mainTab.setTabEnabled(1,False)
        
        #print(self.specfileWidget.mainLayout.count())
        
        # ugly!!!!
        self.specfileWidget.mainLayout.itemAt(self.specfileWidget.mainLayout.count()-1).widget().setVisible(False)
        self.specfileWidget.mainLayout.itemAt(self.specfileWidget.mainLayout.count()-2).widget().setVisible(False)
    
        self.sourceSelector.sigSourceSelectorSignal.connect(self._onSourceSelected)
        
        self.specfileWidget.sigReplaceSelection.connect(self._onReplaceSelection)
        self.nexusWidget.hdf5Widget.sigHDF5WidgetSignal.connect(self._onNEXUSChange)        
        
        
        pathSelector = qt.QSplitter(self)
        pathSelector.setOrientation(qt.Qt.Horizontal)
        qt.QLabel("image path:",pathSelector)
        
        
        
        self.pathedit = qt.QLineEdit(pathSelector)
        
        acceptButton= qt.QPushButton(pathSelector)
        acceptButton.setIcon(qt.QIcon(qt.QPixmap(icons.selected)))
        acceptButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        acceptButton.setToolTip("Accept path to imagefiles")
        acceptButton.clicked.connect(self._onAcceptImagePath)
        
        openButton = qt.QPushButton(pathSelector)
        openButton.setIcon(qt.QIcon(qt.QPixmap(icons.fileopen)))
        openButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        openButton.setToolTip("Choose the directory to the images")
        openButton.clicked.connect(self._onSelectImageFolder)
        
        self.slider = qt.QSlider()
        self.slider.setOrientation(qt.Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)
        
        
        self.mainLayout.addWidget(self.sourceSelector)
        self.mainLayout.addWidget(pathSelector)

        maintab.addTab(self.specfileWidget,"Spec file")
        maintab.addTab(self.nexusWidget,"NEXUS")
        
        self.mainLayout.addWidget(maintab)
        
        #self.mainLayout.addWidget(self.slider)
        
        
        
        self.mainwidget.setLayout(self.mainLayout)
        self.setCentralWidget(self.mainwidget)

        self.toolbar = qt.QToolBar("Image selector",self)
        
        
        righticon = qt.QIcon(qt.QPixmap(icons.rotate_right))
        lefticon = qt.QIcon(qt.QPixmap(icons.rotate_left))
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
        decreaseImageNo = self.toolbar.addAction(lefticon,"previous image")
        increaseImageNo = self.toolbar.addAction(righticon,"next image")
        
        
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
        
        
        self.integrateTab = qt.QWidget()
        
        self.integrateTabLayout = qt.QGridLayout()
        
        self.loadallButton = qt.QPushButton("load all")
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
        
        self.integrateTabLayout.addWidget(self.loadallButton,0,0)
        self.integrateTabLayout.addWidget(self.showMaxButton,1,0)
        self.integrateTabLayout.addWidget(self.showSumButton,2,0)
        
        self.integrateTabLayout.addWidget(qt.QLabel("hier kommt integration basierend\n auf dem masking widget"),0,1)
        
        self.integrateTab.setLayout(self.integrateTabLayout)
        
        self.specfileWidget.mainTab.addTab(self.integrateTab,"Integrate")
        self.specfileWidget.mainTab.setCurrentWidget(self.integrateTab)
        
        #QDataSource.QDataSource()
        
    #def setImageNo(self,imageno):
    #    self.slider.setValue(imageno)
    
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
        
    def _onSourceSelected(self,ddict):
        event = ddict['event']
        if (event == 'SourceSelected' or event == 'NewSourceSelected' or event == 'SourceReloaded'):
            dataname = ddict['sourcelist'][0]
            datatype = QDataSource.getSourceType(dataname)
            datapath,_ = os.path.split(dataname)
            
            
            if datatype == NexusDataSource.SOURCE_TYPE:
                shutil.copy(dataname,datapath + '/temp.h5')
                ds = QDataSource.QDataSource(datapath + '/temp.h5',datatype)
                self.nexusWidget.setDataSource(ds)
            elif datatype == SpecFileDataSource.SOURCE_TYPE:
                ds = QDataSource.QDataSource(dataname,datatype)
                self.specfileWidget.setDataSource(ds)
            else:
                warnings.warn("not implemented data source: %s" % dataname)
            self.selectedScan = None
            self.sigScanChanged.emit([])

    
    def _onReplaceSelection(self,sel_list):
        print(sel_list)
        self.sigScanChanged.emit(sel_list)
        
    def _onNEXUSChange(self,ddict):
        if ddict['event'] == "itemDoubleClicked":
            self.sigScanChanged.emit(ddict)
        
    def _onSelectImageFolder(self):
        folder = PyMcaFileDialogs.getExistingDirectory(self,"select directory containing the images for the current scan")
        if(len(folder)):
            self.pathedit.setText(folder)
            self._onAcceptImagePath()
        
        
    def _onAcceptImagePath(self):
        #print(self.pathedit.text())
        self.sigImagePathChanged.emit(self.pathedit.text())

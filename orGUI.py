# -*- coding: utf-8 -*-
import sys
import os
from PyMca5.PyMcaGui import PyMcaQt as qt
import warnings

from PyMca5.PyMcaGui.io.QSpecFileWidget import QSpecFileWidget
from PyMca5.PyMcaGui.pymca import QDataSource
from PyMca5.PyMcaGui.io.QSourceSelector import QSourceSelector
from PyMca5.PyMcaGui import PyMca_Icons as icons
from PyMca5.PyMcaGui.misc import NumpyArrayTableWidget, NumpyArrayTableModel
from PyMca5.PyMcaGui.misc import FrameBrowser
from PyMca5.PyMcaGui.misc.NumpyArrayTableView import HorizontalHeader, VerticalHeader
from PyMca5.PyMcaGui.io import PyMcaFileDialogs
#from SilxMaskImageWidget import SilxMaskImageWidget

import silx.gui.plot
from silx.gui.plot import items
from silx.gui.colors import Colormap
import weakref

import silx
from silx.utils.weakref import WeakMethodProxy
from silx.gui.plot.Profile import ProfileToolBar
from silx.gui.plot.AlphaSlider import NamedImageAlphaSlider

from NumpyView import NumpyArrayEditTableWidget
from QSpecScanSelector import QSpecScanSelector
from QReflectionSelector import QReflectionSelector
from QUBCalculator import QUBCalculator
import numpy as np
from datautils.xrayutils import HKLVlieg
#if __name__ == "__main__":
#    os.chdir("..")
    
from datautils.xrayutils.id31_tools_5 import Fastscan, BlissScan
from datautils.xrayutils.P212_tools import CrudeThScan

QTVERSION = qt.qVersion()
DEBUG = 0


class orGUI(qt.QMainWindow):
    def __init__(self,parent=None):
        qt.QMainWindow.__init__(self, parent)
        
        self.resetZoom = True
        
        self.fscan = None
        
        self.currentImageLabel = None
        self.currentAddImageLabel = None
        
        selectorDock = qt.QDockWidget()
        selectorDock.setAllowedAreas(qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea)
        self.scanSelector = QSpecScanSelector()
        selectorDock.setWidget(self.scanSelector)
        self.addDockWidget(qt.Qt.LeftDockWidgetArea,selectorDock)
    
        self.imagepath = ''
        self.imageno = 0
        
        #self.reflections = np.array([])
    
        ubWidget = qt.QWidget()
        ubLayout = qt.QHBoxLayout()
        self.ubcalc = QUBCalculator("./config")
        self.ubcalc.sigNewReflection.connect(self._onNewReflection)
        
        

        self.centralPlot = Plot2DHKL(self.newXyHKLConverter(),parent=self)
        self.centralPlot.setDefaultColormap(Colormap(name='jet',normalization='log'))
        self.centralPlot.setCallback(self._graphCallback)

        self.scanSelector.sigImageNoChanged.connect(self._onSliderValueChanged)
        
        
        self.scanSelector.sigImagePathChanged.connect(self._onImagePathChanged)
        self.scanSelector.sigScanChanged.connect(self._onScanChanged)
        
        self.scanSelector.loadallButton.clicked.connect(self._onLoadAll)
        self.scanSelector.showMaxButton.toggled.connect(self._onMaxToggled)
        self.scanSelector.showSumButton.toggled.connect(self._onSumToggled)
        
        self.alphaslider = NamedImageAlphaSlider(self,self.centralPlot,self.currentAddImageLabel)
        self.alphaslider.setOrientation(qt.Qt.Horizontal)
        self.alphaslider.setEnabled(True)
        
        self.scanSelector.integrateTabLayout.addWidget(self.alphaslider,2,1)
        


        
        toolbar = self.scanSelector.getScanToolbar()
        self.centralPlot.addToolBar(qt.Qt.BottomToolBarArea,toolbar)
        

        self.setCentralWidget(self.centralPlot)
        
        
        
        
        #self.reflTable.view._model.dataChanged.connect(printmodel)
        #self.reflTable.setArrayData(np.array([0,0,0,0,10,10],dtype=np.float))
        ubDock = qt.QDockWidget()
        ubDock.setAllowedAreas(qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea | qt.Qt.BottomDockWidgetArea)
        
        self.reflectionSel = QReflectionSelector(self.centralPlot)
        self.reflectionSel.sigQueryImageChange.connect(self._onChangeImage)
        self.reflectionSel.sigQuerySaveReflections.connect(self._onSaveReflections)
        self.reflectionSel.sigQueryLoadReflections.connect(self._onLoadReflections)
        
        
        self.ubcalc.setReflectionHandler(self.getReflections)
        
        self.allimgsum = None
        self.allimgmax = None

        
        ubLayout.addWidget(self.ubcalc)
        ubLayout.addWidget(self.reflectionSel)
        
        ubWidget.setLayout(ubLayout)
        ubDock.setWidget(ubWidget)
        self.addDockWidget(qt.Qt.BottomDockWidgetArea,ubDock)
        
    def getReflections(self):
        hkls = []
        angles = []
        for i in self.reflectionSel.reflections:
            refl = self.reflectionSel.reflections[i]
            #print(refl.xy)
            delta, gamma = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([refl.xy[0]]),np.array([refl.xy[1]]),self.ubcalc.mu)
            delta = float(delta); gamma = float(gamma)
            pos = [self.ubcalc.mu,delta,gamma,self.imageNoToOmega(refl.imageno),self.ubcalc.chi,self.ubcalc.phi]
            #print(pos)
            hkls.append(refl.hkl)
            angles.append(pos)
        return np.array(hkls), np.array(angles)
        
    def _onSaveReflections(self):
        fileTypeList = ['dat Files (*.dat)', 'txt Files (*.txt)','All files (*)']
        newfile, filetype = PyMcaFileDialogs.getFileList(parent=self,
                                               filetypelist=fileTypeList,
                                               message="please select a file",
                                               mode="SAVE",
                                               getfilter=True,
                                               single=True)
        if not newfile:
            return
        filename = newfile[0]
        hkls, angles = self.getReflections()
        angles = np.rad2deg(angles)
        angles[:,3] *= -1 # om to th
        hklangles = np.concatenate([hkls,angles],axis=1)
        np.savetxt(filename,hklangles,header="H K L mu del gam th chi phi",fmt='%.5f')
        
        #print(hkls,angles)
        
    def _onLoadReflections(self):
        fileTypeList = ['dat Files (*.dat)', 'txt Files (*.txt)','All files (*)']
        newfile, filetype = PyMcaFileDialogs.getFileList(parent=self,
                                               filetypelist=fileTypeList,
                                               message="please select a file",
                                               mode="LOAD",
                                               getfilter=True,
                                               single=True)
        if not newfile:
            return
        filename = newfile[0]
        hklangles = np.loadtxt(filename,skiprows=1)
        hkls, angles = hklangles[:,:3], np.deg2rad(hklangles[:,3:])
        self.getReflections()
        angles = np.rad2deg(angles)
        angles[:,3] *= -1 # om to th
        hklangles = np.concatenate([hkls,angles],axis=1)
        
        #print(hkls,angles)
        

    def _onNewReflection(self,refl):
        [hkl,x,y,omega] = refl
        try:
            imageno = self.omegaToImageNo(omega)
        except:
            warnings.warn("Not xmirrored: Didn't find the corresponding image")
            [hkl,x,y,omega] = self.ubcalc.calcReflection(hkl,True)
            try:
                imageno = self.omegaToImageNo(omega)
            except Exception as e:
                qt.QMessageBox.warning(self,"Could not find reflection","Didn't find the corresponding reflection on any image.\n%s" % str(e))
                return
        eventdict = {'x' : x, 'y': y}
        self.reflectionSel.addReflection(eventdict,imageno,hkl)
            
    def newXyHKLConverter(self):
        def xyToHKL(x,y):
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([x]),np.array([y]),self.ubcalc.mu)
            pos = [self.ubcalc.mu,delta,gamma,self.imageNoToOmega(self.imageno),self.ubcalc.chi,self.ubcalc.phi]
            pos = HKLVlieg.crystalAngles(pos,self.ubcalc.n)
            hkl = self.ubcalc.angles.anglesToHkl(pos)
            #print(self.ubcalc.crystal)
            return hkl.A1
        return xyToHKL
        
    def omegaToImageNo(self,omega):
        if self.fscan is not None:
            omrad = np.deg2rad(self.fscan.omega)
            ommax = np.amax(omrad)
            ommin = np.amin(omrad)
            #print(ommin,omega,ommax)
            if omega < ommin or omega > ommax:
                raise Exception("omega not in range")
            return np.argmin(np.abs(omrad -omega))
        else:
            raise Exception("No Scan selected")
            
    def imageNoToOmega(self,imageno):
        if self.fscan is not None:
            return np.deg2rad(self.fscan.omega[imageno])
        else:
            return 0.
        
        
    def _onScanChanged(self,sel_list):
        self.resetZoom = True
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
                    self.fscan = CrudeThScan(self.specfile,'PE1',r"C:\Timo_loc\P21_2_comissioning\Pt111_HClO4_0.4\PE1\dark00001.tif.gz")
                    self.imageno = 0
                    self.imagepath = self.fscan.path + "/" + 'PE1'
                self.reflectionSel.setImage(self.imageno)
                if self.imagepath != '':
                    self.fscan.set_image_folder(self.imagepath)
                    self.plotImage()
                    
                    #self.scanSelector.slider.setMinimum(0)
                    #self.scanSelector.slider.setMaximum(self.fscan.nopoints-1)
                    self.scanSelector.setTh(self.fscan.th)
                    #print(self.fscan.nopoints)
                    #self.readAllImages()
                    #print(self.centralPlot._callback)
                    
            else:
                self.scanSelector.setRange(0,0)
                self.imageno = 0
                self.reflectionSel.setImage(self.imageno)
                #print(self.centralPlot._callback)
        else:
            self.hdffile = sel_list['file']
            self.scanname = sel_list['name'].strip("/")
            self.fscan = BlissScan(self.hdffile,self.scanname)
            
                
            
            
    def _onImagePathChanged(self,path):
        #print("newpath %s" % path)
        self.imagepath = path
        if self.fscan is not None:
            self.fscan.set_image_folder(self.imagepath)
            self.plotImage()
            self.scanSelector.setTh(self.fscan.th)
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
        if self.fscan is not None and self.imagepath != '':
            self.loadAll()
            
    def loadAll(self):
        try:
            image = self.fscan.get_raw_PE_img(0)
        except Exception as e:
            print("no images found! %s" % e)
            return
        self.allimgsum = np.zeros_like(image.img)
        self.allimgmax = np.zeros_like(image.img)
        progress = qt.QProgressDialog("Reading images","abort",0,len(self.fscan),self)
        progress.setWindowModality(qt.Qt.WindowModal)
        for i in range(len(self.fscan)):
            image = self.fscan.get_raw_PE_img(i)
            self.allimgsum += image.img
            self.allimgmax = np.maximum(self.allimgmax,image.img)
            progress.setValue(i)
            if progress.wasCanceled():
                #self.allimgsum = None
                #self.allimgmax = None
                break
        progress.setValue(len(self.fscan))
        
    def _onMaxToggled(self,value):
        if self.scanSelector.showSumButton.isChecked():
            self.scanSelector.showSumButton.setChecked(False)
        if value:
            if self.allimgmax is not None:
                self.currentAddImageLabel = self.centralPlot.addImage(self.allimgmax,legend="special",
                                                               replace=False,resetzoom=False,copy=True,z=1)
                self.centralPlot.setActiveImage(self.currentAddImageLabel)
                self.alphaslider.setLegend(self.currentAddImageLabel)
            else:
                self.scanSelector.showMaxButton.setChecked(False)
        else:
            if self.currentAddImageLabel is not None:
                self.centralPlot.removeImage(self.currentAddImageLabel)
                self.currentAddImageLabel = None
                self.centralPlot.setActiveImage(self.currentImageLabel)
        
    def _onSumToggled(self,value):
        if self.scanSelector.showMaxButton.isChecked():
            self.scanSelector.showMaxButton.setChecked(False)
        if value:
            if self.allimgsum is not None:
                self.currentAddImageLabel = self.centralPlot.addImage(self.allimgsum,legend="special",
                                                               replace=False,resetzoom=False,copy=True,z=1)
                self.centralPlot.setActiveImage(self.currentAddImageLabel)
                self.alphaslider.setLegend(self.currentAddImageLabel)
            else:
                self.scanSelector.showSumButton.setChecked(False)
        else:
            if self.currentAddImageLabel is not None:
                self.centralPlot.removeImage(self.currentAddImageLabel)
                self.currentAddImageLabel = None
                self.centralPlot.setActiveImage(self.currentImageLabel)
        
        
        
    def plotImage(self,key=0):
        try:
            image = self.fscan.get_raw_PE_img(key)
            if self.currentImageLabel is not None:
                self.centralPlot.removeImage(self.currentImageLabel)

            self.currentImageLabel = self.centralPlot.addImage(image.img,legend="No %i" % key,
                                                               replace=False,resetzoom=self.resetZoom,copy=True)
            if self.currentAddImageLabel is None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
            self.resetZoom = False
            self.imageno = key
            self.reflectionSel.setImage(self.imageno)
        except Exception as e:
            print("no image %s" % e)
        
        
        
    def _graphCallback(self,eventdict):
        #print(eventdict)
        if eventdict['event'] == 'mouseDoubleClicked':
            #newReflection = np.array([1,1,1,self.imageno,eventdict['x'],eventdict['y']])
            hkl = self.centralPlot.xyHKLConverter(eventdict['x'],eventdict['y'])
            self.reflectionSel.addReflection(eventdict,self.imageno,hkl)
            
        if eventdict['event'] == 'markerMoved':
            self.reflectionSel.moveReflection(eventdict['label'],[eventdict['x'],eventdict['y']])
            self.reflectionSel.setReflectionActive(eventdict['label'])
        if eventdict['event'] == 'markerClicked':
            self.reflectionSel.setReflectionActive(eventdict['label'])
        
        
class Plot2DHKL(silx.gui.plot.PlotWindow):
    def __init__(self,xyHKLConverter,parent=None,backend=None):
        self.xyHKLConverter = xyHKLConverter
        
        
        posInfo = [
            ('X', lambda x, y: x),
            ('Y', lambda x, y: y),
            ('H', lambda x, y: self.xyHKLConverter(x,y)[0]),
            ('K', lambda x, y: self.xyHKLConverter(x,y)[1]),
            ('L', lambda x, y: self.xyHKLConverter(x,y)[2]),
            ('Data', WeakMethodProxy(self._getImageValue))]
        
        super(Plot2DHKL, self).__init__(parent=parent, backend=backend,
                             resetzoom=True, autoScale=False,
                             logScale=False, grid=False,
                             curveStyle=False, colormap=True,
                             aspectRatio=True, yInverted=True,
                             copy=True, save=True, print_=True,
                             control=False, position=posInfo,
                             roi=False, mask=True)
        
        if parent is None:
            self.setWindowTitle('Plot2D')
        self.getXAxis().setLabel('Columns')
        self.getYAxis().setLabel('Rows')

        if silx.config.DEFAULT_PLOT_IMAGE_Y_AXIS_ORIENTATION == 'downward':
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

    def getProfileToolbar(self):
        """Profile tools attached to this plot

        See :class:`silx.gui.plot.Profile.ProfileToolBar`
        """
        return self.profile


    def getProfilePlot(self):
        """Return plot window used to display profile curve.

        :return: :class:`Plot1D`
        """
        return self.profile.getProfilePlot()
        
    


a = qt.QApplication(sys.argv)

mainWindow = orGUI()
mainWindow.show()



a.lastWindowClosed.connect(a.quit)

a.exec_()


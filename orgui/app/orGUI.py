# -*- coding: utf-8 -*-
import sys
import os
from silx.gui import qt
import warnings


import silx.gui.plot
from silx.gui.plot import items
from silx.gui.colors import Colormap
import weakref

import silx
from silx.utils.weakref import WeakMethodProxy
from silx.gui.plot.Profile import ProfileToolBar
from silx.gui.plot.AlphaSlider import NamedImageAlphaSlider
from silx.gui.dialog import ImageFileDialog
import traceback


from .QSpecScanSelector import QSpecScanSelector
from .QReflectionSelector import QReflectionSelector
from .QUBCalculator import QUBCalculator
import numpy as np
from datautils.xrayutils import HKLVlieg, CTRcalc
from datautils.xrayutils import ReciprocalNavigation as rn


#if __name__ == "__main__":
#    os.chdir("..")

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import sys
#sys.path.append('/home/fuchstim/repos/datautils/datautils/xrayutils')
from datautils.xrayutils.id31_tools_5 import BlissScan_EBS, Fastscan, BlissScan
#from datautils.xrayutils.P212_tools import CrudeThScan, FioFastsweep
#from P212_tools import BlissScan_EBS#CrudeThScan, FioFastsweep

QTVERSION = qt.qVersion()
DEBUG = 0


class orGUI(qt.QMainWindow):
    def __init__(self,configfile,parent=None):
        qt.QMainWindow.__init__(self, parent)
        
        self.resetZoom = True
        
        self.fscan = None
        
        self.filedialogdir = os.getcwd()
        
        self.currentImageLabel = None
        self.currentAddImageLabel = None
        
        selectorDock = qt.QDockWidget()
        selectorDock.setAllowedAreas(qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea)
        self.scanSelector = QSpecScanSelector(self)
        selectorDock.setWidget(self.scanSelector)
        self.addDockWidget(qt.Qt.LeftDockWidgetArea,selectorDock)
    
        self.imagepath = ''
        self.imageno = 0
        
        menu_bar = qt.QMenuBar() 
        file = menu_bar.addMenu("&File")
        file.addAction(self.scanSelector.openFileAction)
        file.addAction(self.scanSelector.refreshFileAction)
        file.addAction(self.scanSelector.closeFileAction)
        
        calcCTRsAvailableAct = qt.QAction("Calculate available CTRs",self)
        calcCTRsAvailableAct.triggered.connect(self._onCalcAvailableCTR)
        rs = menu_bar.addMenu("&Reciprocal space")
        rs.addAction(calcCTRsAvailableAct)
        
        
        self.setMenuBar(menu_bar)
        
        
        #self.reflections = np.array([])
    
        ubWidget = qt.QWidget()
        ubLayout = qt.QHBoxLayout()
        self.ubcalc = QUBCalculator(configfile)
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
        
        self.ubcalc.sigPlottableMachineParamsChanged.connect(self._onPlotMachineParams)
        self.allimgsum = None
        self.allimgmax = None

        
        ubLayout.addWidget(self.ubcalc)
        ubLayout.addWidget(self.reflectionSel)
        
        ubWidget.setLayout(ubLayout)
        ubDock.setWidget(ubWidget)
        self.addDockWidget(qt.Qt.BottomDockWidgetArea,ubDock)
        
    def _onCalcAvailableCTR(self):
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
                xtal = CTRcalc.UnitCell.fromBULfile(filename)
            except Exception:
                qt.QMessageBox.critical(self,"Cannot open xtal", "Cannot open:\n%s" % traceback.format_exc())
                return
        else:
            qt.QMessageBox.critical(self,"Cannot open xtal", "File extension not understood")
            return
        
        try:
            ommin = np.amin(self.fscan.omega)
            ommax = np.amax(self.fscan.omega)
            dc = self.ubcalc.detectorCal
            mu = np.rad2deg(self.ubcalc.mu)
            ub = self.ubcalc.ubCal
            xtal.setEnergy(ub.getEnergy()*1e3)
            hk = rn.thscanCTRs(xtal,ub,mu,dc,(ommin,ommax))
        except Exception:
            qt.QMessageBox.critical(self,"Cannot calculate CTR locatons", "Cannot calculate CTR locatons:\n%s" % traceback.format_exc())
            return
        
        #making the hk list of arrays into a reasonable string
        hkstring="";
        rowcount=0;
        for i in hk:
            if rowcount!=5:
                hkstring=hkstring+"%s "%i;
                rowcount=rowcount+1;
            else:
                hkstring=hkstring+"%s \n"%i;
                rowcount=0;
                
        #Question dialog for saving the possible CTR locations        
        clickedbutton=qt.QMessageBox.question(self, 'Saving CTR locations...', 'Do you want to save the following positions: \n' + hkstring +"?");
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
            np.savetxt(filename,hk,header="H K",fmt='%d')
            
        
    def getReflections(self):
        hkls = []
        angles = []
        for i in self.reflectionSel.reflections:
            refl = self.reflectionSel.reflections[i]
            #print(refl.xy)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([refl.xy[1]]),np.array([refl.xy[0]]),self.ubcalc.mu)
            delta = float(delta); gamma = float(gamma)
            pos = [self.ubcalc.mu,delta,gamma,self.imageNoToOmega(refl.imageno),self.ubcalc.chi,self.ubcalc.phi]
            #print(pos)
            hkls.append(refl.hkl)
            angles.append(pos)
        return np.array(hkls), np.array(angles)
        
    def _onSaveReflections(self):
        
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
        hkls, angles = self.getReflections()
        if angles.size > 0:
            angles = np.rad2deg(angles)
            angles[:,3] *= -1 # om to th
            hklangles = np.concatenate([hkls,angles],axis=1)
            np.savetxt(filename,hklangles,header="H K L mu del gam th chi phi",fmt='%.5f')
        
        #print(hkls,angles)
        
    def _onLoadReflections(self):
        fileTypeDict = {'dat Files (*.dat)': '.dat', 'txt Files (*.txt)': '.txt', 'All files (*)': '', }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Open file with reflections",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.filedialogdir = os.path.splitext(filename)[0]
        try:
            hklangles = np.loadtxt(filename,skiprows=1)
            hkls, angles = hklangles[:,:3], np.deg2rad(hklangles[:,3:])
            
            #self.getReflections()
            angles = np.rad2deg(angles)
            angles[:,3] *= -1 # om to th
            hklangles = np.concatenate([hkls,angles],axis=1)
            for refl in hklangles:
                self._onNewReflection(refl)
        except Exception:
            qt.QMessageBox.critical(self,"Error during loading of reflections","Error during loading of reflections.\n%s" % traceback.format_exc())
        
        #print(hkls,angles)
        
    def _onPlotMachineParams(self,paramslist):
        [cp,azimxy,polax] = paramslist
        self.centralPlot.addMarker(cp[0],cp[1],legend="CentralPixel",text="CP",color='yellow',symbol='+')
        self.centralPlot.addMarker(azimxy[0],azimxy[1],legend="azimuth",text="Azim",color='yellow',symbol='+')

    def _onNewReflection(self,refl):
        [hkl,x,y,omega] = refl
        notmirrored = [hkl,x,y,omega]
        
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
                
                qt.QMessageBox.warning(self,"Could not find reflection","Didn't find the corresponding reflection on any image.\nError: %s\nShould be at location%s" % (str(e),errormsg))
                return
        eventdict = {'x' : x, 'y': y}
        self.reflectionSel.addReflection(eventdict,imageno,hkl)
            
    def newXyHKLConverter(self):
        def xyToHKL(x,y):
            #print("xytoHKL:")
            #print("x,y = %s, %s" % (x,y))
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([y]),np.array([x]),self.ubcalc.mu)
            #print(self.ubcalc.detectorCal)
            #print(x,y)
            #print(self.ubcalc.detectorCal.tth(np.array([y]),np.array([x])))
            
            pos = [self.ubcalc.mu,delta[0],gamma[0],self.imageNoToOmega(self.imageno),self.ubcalc.chi,self.ubcalc.phi]
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
        
        
    def _onScanChanged(self,sel_list):
        self.resetZoom = True
        print(sel_list)

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
                    #self.fscan = CrudeThScan(self.specfile,'PE1',r"C:\Timo_loc\P21_2_comissioning\Pt111_HClO4_0.4\PE1\dark00001.tif.gz")
                    self.fscan = FioFastsweep(self.specfile)
                    self.imageno = 0
                    #self.imagepath = self.fscan.path + "/" + 'PE1'
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
            #self.scanname = sel_list['name'].strip("/")
            try:
                msg = qt.QMessageBox(self)
                msg.setWindowTitle("Loading Scan")
                msg.setText("Loading Scan. This might take a while...")
                msg.setStandardButtons(qt.QMessageBox.Cancel)
                msg.setModal(True)
                msg.show()
                ch5523 = sel_list.get('ch5523',False)
                if ch5523:
                    self.fscan = BlissScan(self.hdffile,sel_list['name'].strip('/'))
                else:
                    if 'node' in sel_list:
                        self.fscan = BlissScan_EBS(sel_list['node'],sel_list['scanno'])
                    else:
                        self.fscan = BlissScan_EBS(self.hdffile,sel_list['scanno'])
                self.plotImage()
                self.scanSelector.setTh(self.fscan.th)
                msg.hide()
                self._onLoadAll()
                self.scanSelector.showMaxButton.setChecked(False)
                self.scanSelector.showMaxButton.setChecked(True)
            except Exception:
                msg.hide()
                qt.QMessageBox.critical(self,"Cannot open scan", "Cannot open scan:\n%s" % traceback.format_exc())
                
            
            
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
        if self.fscan is not None:
            self.loadAll()
            
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
        for i in range(len(self.fscan)):
            image = self.fscan.get_raw_img(i)
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
            image = self.fscan.get_raw_img(key)
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
        
def main(configfile):
	a = qt.QApplication(['orGUI'])
	mainWindow = orGUI(configfile)
	mainWindow.show()
	a.lastWindowClosed.connect(a.quit)
	return a.exec_()
	
            
if __name__ == '__main__':
	main("./config")

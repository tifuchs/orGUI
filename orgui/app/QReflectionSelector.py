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

import sys
import traceback
import os
from io import StringIO
from silx.gui import qt
from silx.gui import icons
from .ArrayTableDialog import ArrayEditWidget
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import ClassVar
import copy

import numpy as np
from ..datautils.xrayutils import HKLVlieg

from .. import resources
from . import qutils
from . import imagePeakFinder

@dataclass
class HKLReflection:
    xy: np.ndarray
    hkl: np.ndarray
    imageno: int
    identifier : str
        
    @classmethod
    def fromStr(cls,string):
        h,k,l,x,y,imageno = string.split()
        refl = cls([x,y],[h,k,l],int(imageno))
        return refl
    
    @staticmethod
    def fromArray(arr, identifier):
        reflist = []
        if arr.ndim == 1:
            h,k,l,x,y,imageno = arr
            refl = HKLReflection([x,y],[h,k,l],int(imageno), identifier)
            return refl
        if len(identifier) != len(arr):
            raise ValueError("Mismatch between number of reflections in array and identifiers.")
        for row, ident in zip(arr, identifier):
            h,k,l,x,y,imageno = row
            reflist.append(HKLReflection([x,y],[h,k,l],int(imageno), ident))
        return reflist
        
    def __array__(self):
        return np.concatenate([self.hkl,self.xy,[self.imageno]])

class QReflectionSelector(qt.QWidget):
    sigQueryImageChange = qt.pyqtSignal(int)
    sigQuerySaveReflections = qt.pyqtSignal()
    sigQueryLoadReflections = qt.pyqtSignal()
    sigQueryCenterPlot = qt.pyqtSignal(object)
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self,plot, ubcalc, parent=None):
        qt.QWidget.__init__(self, parent=None)
        layout = qt.QVBoxLayout()
        #self.setOrientation(qt.Qt.Vertical)
        self.plot = plot
        self.ubcalc = ubcalc
        self.orparent = parent
        self.pkImgDiag = PeakImgRangeDialog([0, 360], 'NoAxis')
        self.peakReflDialog = PeakReflDialog(self)
        
        self.reflections = []
        self.reflBragg = []
        
        self.nextNo = 0
        self.imageno = 0
        
        self._showReferenceReflections = True
        self.showBraggReflections = False
        
        self.activeReflection = None
        
        self.plot.sigKeyPressDelete.connect(self._onDelete)

        editorTabWidget = qt.QTabWidget()

        self.refleditor = ArrayEditWidget(True, 1)
        header = ['H', 'K', 'L', 'x', 'y', 'imageno']
        self.refleditor.setArrayData(np.array([]).reshape((-1,6)), editable=True, header=header)
        #self.updateEditor()
        self.refleditor.model.dataChanged.connect(self._onDataChanged)
        self.refleditor.sigRowAdded.connect(self._onRowAdded)
        self.refleditor.sigDataLoaded.connect(self.reflectionsFromEditor)
        self.refleditor.sigRowsDeleted.connect(self._onRowsDeleted)
        self.refleditor.view.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.refleditor.view.selectionModel().currentRowChanged.connect(self._onSelectedRowChanged)

        editorTabWidget.addTab(self.refleditor, "intrinsic coord")

        self.refleditor_angles = ArrayEditWidget(True, -1, False)
        header = ['H', 'K', 'L', 'alpha', 'delta', 'gamma', 'omega', 'chi', 'phi']
        self.refleditor_angles.setArrayData(np.array([]).reshape((-1,9)), editable=False, header=header)
        self.refleditor_angles.view.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.refleditor_angles.view.selectionModel().currentRowChanged.connect(self._onSelectedRowChanged)
        editorTabWidget.addTab(self.refleditor_angles, "SIXC angles")
        
        layout.addWidget(editorTabWidget)
        #self.addWidget(editorTabWidget)
        
        self.controlToolbar = qt.QToolBar()
        self.controlToolbar.setMovable(False)
        self.controlToolbar.setOrientation(qt.Qt.Horizontal)
        
        spacer_widget = qt.QWidget()
        spacer_widget.setSizePolicy(qt.QSizePolicy.Expanding, qt.QSizePolicy.Preferred)
        
        self.setImageAct = qt.QAction(resources.getQicon('select-image'), "select Image", self)
        self.setImageAct.setToolTip("set this image for the chosen reflection")
        self.setImageAct.triggered.connect(self._onSetImage)
        
        self.gotoImageAct = qt.QAction(resources.getQicon('search-image'), "search image",self)
        self.gotoImageAct.setToolTip("diplay image with the chosen reflection")
        self.gotoImageAct.triggered.connect(self._onGotoImage)
        
        self.peakSearchInImage = qt.QAction(resources.getQicon('search-peak'), "2D peak search",self) 
        self.peakSearchInImage.setToolTip("perform peak search in the image")
        self.peakSearchInImage.triggered.connect(self._onSinglePeakSearch)
        
        self.addBraggReflAct = qt.QAction( "add Bragg reflection",self) 
        self.addBraggReflAct.setToolTip("add the next Bragg reflection, which improves the U matrix the most")
        self.addBraggReflAct.triggered.connect(self._onAddBraggReflection)
        
        self.controlToolbar.addAction(self.gotoImageAct)
        self.controlToolbar.addSeparator()
        self.controlToolbar.addAction(self.setImageAct)
        self.controlToolbar.addAction(self.peakSearchInImage)
        self.controlToolbar.addWidget(spacer_widget)
        self.controlToolbar.addAction(self.addBraggReflAct)
        
        
        layout.addWidget(self.controlToolbar)
        
        self.setLayout(layout)
        
    def _onSinglePeakSearch(self):
        
        if self.activeReflection is not None:
            try:
                refl = self.reflections[self.indexReflection(self.activeReflection)]
                if self.orparent.fscan is not None:
                    self.pkImgDiag.set_axis(self.orparent.fscan.axis, self.orparent.fscan.axisname)
                    if self.pkImgDiag.exec() == qt.QDialog.Accepted:
                        roidict = self.pkImgDiag.get_roi()
                        #mu, om = self.orparent.getMuOm(refl.imageno)
                        axis = self.orparent.fscan.axis
                        ax_min = roidict['from'] + axis[refl.imageno]
                        ax_max = roidict['to'] + axis[refl.imageno]
                        xy_start = refl.xy
                        com_d = self.searchPeak(xy_start, 
                                                roidict['vsize'], 
                                                roidict['hsize'],
                                                [ax_min, ax_max])
                        ax_min = roidict['from_2'] + com_d['axis_com']
                        ax_max = roidict['to_2'] + com_d['axis_com']                                     
                        xy_start = com_d['xy']
                        com_d = self.searchPeak(xy_start, 
                                                roidict['vsize_2'], 
                                                roidict['hsize_2'],
                                                [ax_min, ax_max])
                        
                        refl.xy = com_d['xy']
                        refl.imageno = np.argmin(np.abs(axis - com_d['axis_com']))
                        self.updateEditor()
                        self.sigQueryImageChange.emit(refl.imageno)
                        
            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot calculte peak position:\n%s" % e, traceback.format_exc())
        
    def searchPeak(self, xy_start, vsize, hsize, ax_range):
        kwargs = {
            'excluded_images' : self.orparent.excludedImagesDialog.getData(),
            'max_workers' : self.orparent.numberthreads
        }
        if self.orparent.centralPlot.getMaskToolsDockWidget().getSelectionMask() is not None:
            kwargs['mask'] = self.orparent.centralPlot.getMaskToolsDockWidget().getSelectionMask() > 0.0
        com_d = imagePeakFinder.find_COM_Image(xy_start, 
                                                vsize, 
                                                hsize,
                                                self.orparent.fscan,
                                                ax_range,
                                                **kwargs)
        return com_d
                                                
        
        
    def getBraggCandidates(self, recalculate=False):
        if recalculate:
            self.orparent.calcBraggRefl()
        if not self.reflBragg:
            raise ValueError("No Bragg reflections calculated or available")
        
        refl_new = []
        refl_hkl_current = np.vstack([r.hkl for r in self.reflections])
        
        for refl in self.reflBragg:
            for cref in self.reflections:
                if np.allclose(refl.hkl, cref.hkl):
                    break
            else:
                refl_new.append(refl)
                
        refl_new_hkl_Bragg = np.vstack([r.hkl for r in refl_new])
        
        Q_current = self.orparent.ubcalc.ubCal.lattice.getB() @ refl_hkl_current.T
        Q_candidates = self.orparent.ubcalc.ubCal.lattice.getB() @ refl_new_hkl_Bragg.T
        
        refl_new_norm = []
        
        for Qc in Q_candidates.T:
            norm = np.linalg.norm(Q_current.T - Qc, axis=1)
            refl_new_norm.append(np.sum(norm))
            
        refl_new_norm = np.array(refl_new_norm)
        order = np.argsort(refl_new_norm)[::-1]
        
        ordered_refl = []
        for o in order:
            ordered_refl.append(refl_new[o])
        return ordered_refl
        
    def _onAddBraggReflection(self):
        if self.peakReflDialog.isVisible():
            self.peakReflDialog.activateWindow()
            return
        if not self.reflBragg:
            qutils.warning_detailed_message(self, "Error","No Bragg reflection positions calculated.\nUse view->show Bragg reflections to calculate them", '')
            return
        self.peakReflDialog.set_candidates(self.getBraggCandidates())
        self.peakReflDialog.show()
        
        
    def _onDataChanged(self, *dat):
        data = self.refleditor.getData()
        idx_no = dat[0].row()
        identifier = self.reflections[idx_no].identifier
        if data.ndim > 1:
            refl = HKLReflection.fromArray(data[idx_no], identifier)
        else:
            refl = HKLReflection.fromArray(data, identifier)
        self.reflections[idx_no] = refl
        if self.activeReflection != identifier:
            self.setReflectionActive(identifier)
        else:
            self.redrawActiveReflection()
            
    def anglesReflection(self, refl):
        if self.orparent.fscan is None:
            return np.array([np.nan,np.nan,np.nan, np.nan, np.nan, np.nan])
        print(refl.imageno)
        mu, om = self.orparent.getMuOm(refl.imageno)
        
        gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([refl.xy[1]]),np.array([refl.xy[0]]), mu)
        pos = [mu,delta[0],gamma[0],om,self.ubcalc.chi,self.ubcalc.phi]
        pos = HKLVlieg.crystalAngles(pos,self.ubcalc.n)
        return np.array(pos)
        
            
    def _onRowAdded(self, row):
        identifier = 'ref_'+str(self.nextNo)
        if self.refleditor.getData().ndim > 1:
            refl = HKLReflection.fromArray(self.refleditor.getData()[row], identifier)
        else:
            refl = HKLReflection.fromArray(self.refleditor.getData(), identifier)
        self.nextNo += 1
        angles = np.concatenate([refl.hkl,np.rad2deg(self.anglesReflection(refl))])
        angledata = self.refleditor_angles.getData()
        if angledata.ndim > 1:
            angledata = np.insert(angledata, row, angles, axis=0)
        elif angledata.ndim == 1:
            if row == 0:
                angledata = np.vstack([angles, angledata])
            elif row == 1:
                angledata = np.vstack([angledata, angles])
        else:
            angledata = angles[np.newaxis,:]
        self.refleditor_angles.updateArrayData(angledata)
        
        if self._showReferenceReflections:
            self.plot.addMarker(*refl.xy,legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        self.reflections.insert(row, refl)
        self.setReflectionActive(refl.identifier)
        
    def _onRowsDeleted(self, rows):
        idents = []
        for r in rows:
            idents.append(self.reflections[r].identifier)
        for ident in idents:
            self.deleteReflection(ident)
        """
        print(self.refleditor_angles.getData())
        angledata = np.atleast_2d(self.refleditor_angles.getData())
        print(angledata)
        mask = np.ones(len(angledata),dtype=bool)
        print(mask)
        mask[np.array(rows)] = False
        angledata = angledata[mask]
        print(angledata)
        print(angledata.ndim )
        if angledata.ndim == 1:
            angledata = angledata[np.newaxis,:]
        elif angledata.ndim == 0:
            angledata = np.array([]).reshape((-1,9))
        self.refleditor_angles.updateArrayData(angledata)
        """
            
    def _onSelectedRowChanged(self, selected, deselected):
        row = selected.row()
        refl = self.reflections[row]
        if self._showReferenceReflections:
            if self.activeReflection is not None:
                self.setReflectionInactive(self.activeReflection)
            self.plot.removeMarker(refl.identifier)
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
            self.activeReflection = refl.identifier
        else:
            self.activeReflection = refl.identifier
        
    def setBraggReflectionsVisible(self, visible):
        if visible:
            self.showBraggReflections = True
            self.redrawBraggReflections()
        else:
            self.showBraggReflections = False
            self.clearPlotBraggReflections()
        
    def setBraggReflections(self, hkls, yx, angles):
        self.__cached_bragg_refl = hkls, yx, angles
        if self.showBraggReflections:
            self.clearPlotBraggReflections()
        self.reflBragg.clear()
        for i, (hkl, yxr, pos) in enumerate(zip(hkls, yx, angles)):
            identifier = 'bragg_'+str(i)
            if self.orparent.fscan.axisname == 'mu':
                angle_idx = 0
                sign = 1.
            elif self.orparent.fscan.axisname == 'th':
                angle_idx = 3
                sign = -1.
            else:
                qt.QMessageBox.warning(self,"Cannot calculate reflection","Cannot calculate reflection.\n%s is no supported scan axis." % self.orparent.fscan.axisname)
                return 
            try:
                imageno1 = self.orparent.axisToImageNo(np.rad2deg(pos[angle_idx]) * sign)
            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot calculate axis position:\n%s" % e, traceback.format_exc())

            self.reflBragg.append(HKLReflection(yxr[::-1], hkl, imageno1, identifier))
        if self.showBraggReflections:
            self.redrawBraggReflections()
        else:
            self.clearPlotBraggReflections()
            
    def getBraggReflections(self):
        if hasattr(self, "__cached_bragg_refl") and self.showBraggReflections:
            return self.__cached_bragg_refl
        else:
            raise ValueError("Bragg reflections are not calculated or are not shown!")
        
    def setReferenceReflectionsVisible(self, visible):
        if visible:
            self._showReferenceReflections = True
            self.redrawRefReflections()
        else:
            self._showReferenceReflections = False
            self.clearPlotRefReflections()
        
    def setImage(self,imageno):
        self.imageno = imageno
        
    def _onSaveReflections(self):
        self.sigQuerySaveReflections.emit()
        
    def _onLoadReflections(self):
        self.sigQueryLoadReflections.emit()
        
        
    def updateEditor(self):
        if self.reflections:
            array = np.vstack(self.reflections)
            self.refleditor.updateArrayData(array)
            anglearray = []
            for refl in self.reflections:
                anglearray.append(np.concatenate([refl.hkl,np.rad2deg(self.anglesReflection(refl))]))
            anglearray = np.vstack(anglearray)
            self.refleditor_angles.updateArrayData(anglearray)
        else:
            self.refleditor.updateArrayData(np.array([]).reshape((-1,6)))
            self.refleditor_angles.updateArrayData(np.array([]).reshape((-1,9)))
        try:
            idx = self.indexReflection(self.activeReflection)
            self.refleditor.view.selectRow(idx)
            self.refleditor_angles.view.selectRow(idx)
        except:
            if self.reflections:
                self.setReflectionActive(self.reflections[0].identifier)
            else:
                self.activeReflection = None
                return
        
        
    def reflectionsFromEditor(self):
        data = self.refleditor.getData()
        idents = ['ref_'+str(i) for i in range(len(data))]
        self.nextNo = len(data)
        try:
            if len(data) > 1:
                reflist = HKLReflection.fromArray(data, idents)
                self.setReflections(reflist)
            else:
                self.setReflections([])
        except Exception:
            qt.QMessageBox.critical(self,"Cannot read reflections", "Cannot read reflections %s" % traceback.format_exc())
            return
            
    def clearPlotRefReflections(self):
        for refl in self.reflections:
            self.plot.removeMarker(refl.identifier)
            
    def clearPlotBraggReflections(self):
        for refl in self.reflBragg:
            self.plot.removeMarker(refl.identifier)
            
    def redrawBraggReflections(self):
        self.clearPlotBraggReflections()
        if self.showBraggReflections:
            for refl in self.reflBragg:
                #refl = self.reflBragg[identifier]
                #text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl)
                self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,color='g',selectable=False,draggable=False,symbol='x', text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl))
        
            
    def redrawRefReflections(self):
        self.setReflections(copy.deepcopy(self.reflections))
        
    def setReflections(self,refls):
        self.clearPlotRefReflections()
        self.activeReflection = None
        self.nextNo = 0
        self.reflections.clear()
        for refl in refls:
            refl.identifier = 'ref_'+str(self.nextNo)
            if self._showReferenceReflections:
                self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
            
            self.reflections.append(refl)
            self.nextNo += 1
        if self.reflections:
            newactive = self.reflections[0].identifier
            self.setReflectionActive(newactive)
        
    def addReflection(self,eventdict,imageno,hkl=np.array([np.nan,np.nan,np.nan])):
        identifier = 'ref_'+str(self.nextNo)
        refl = HKLReflection([eventdict['x'],eventdict['y']],hkl,imageno,identifier)
        if self._showReferenceReflections:
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        self.reflections.append(refl)
        self.nextNo += 1
        self.updateEditor()
        self.setReflectionActive(identifier)
        return refl

        
    def indexReflection(self, identifier):
        return [ref.identifier for ref in self.reflections].index(identifier)
        
    def setReflectionActive(self,identifier):
        idx = self.indexReflection(identifier)
        refl = self.reflections[idx]
        if self._showReferenceReflections:
            if self.activeReflection is not None:
                self.setReflectionInactive(self.activeReflection)
            self.plot.removeMarker(identifier)
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(identifier),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
            self.activeReflection = identifier
        else:
            self.activeReflection = identifier
        selectModel = self.refleditor.view.selectionModel()
        if idx != selectModel.currentIndex().row():
            self.refleditor.view.selectRow(idx)
            self.refleditor_angles.view.selectRow(idx)
    
    def setReflectionInactive(self,identifier):
        if self._showReferenceReflections:
            self.plot.removeMarker(identifier)
            refl = self.reflections[self.indexReflection(identifier)]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(identifier),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        
    def moveReflection(self,identifier,xy):
        self.reflections[self.indexReflection(identifier)].xy = xy
        self.updateEditor()
        
    def redrawActiveReflection(self):
        if self._showReferenceReflections:
            self.plot.removeMarker(self.activeReflection)
            refl = self.reflections[self.indexReflection(self.activeReflection)]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(self.activeReflection),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
        
    def _onGotoImage(self):
        if self.activeReflection is not None:
            self.sigQueryImageChange.emit(self.reflections[self.indexReflection(self.activeReflection)].imageno)
            self.sigQueryCenterPlot.emit(self.reflections[self.indexReflection(self.activeReflection)].xy)
            
    def _onSetImage(self):
        if self.activeReflection is not None:
            self.reflections[self.indexReflection(self.activeReflection)].imageno = self.imageno
            self.updateEditor()
    
    def _onDelete(self):
        try:
            idx = self.indexReflection(self.activeReflection)
        except:
            return
        refl = self.reflections[idx]
        btn = qt.QMessageBox.question(self, "Delete reflection?", 
                                "Do you want to delete the reflection?\n%s" % refl.hkl,
                                qt.QMessageBox.Yes | qt.QMessageBox.No, qt.QMessageBox.Yes)
        if btn == qt.QMessageBox.Yes:
            self.deleteReflection(self.activeReflection)
    
    def deleteReflection(self,identifier):
        try:
            del self.reflections[self.indexReflection(identifier)]
        except:
            return
        
        if self._showReferenceReflections:
            self.plot.removeMarker(identifier)
        if self.activeReflection == identifier:
            self.activeReflection = None
            if self.reflections:
                newactive = self.reflections[0]
                self.setReflectionActive(newactive.identifier)
        self.updateEditor()
        
        
        
        
        
class QReflectionAnglesDialog(qt.QDialog):
    """Displays the angles for a reflection. Allows for multiple solutions.
    
    """
    
    class QSelectableLabel(qt.QLabel):
        def __init__(self, *args, **kwags):
            qt.QLabel.__init__(self,*args, **kwags)
            self.setTextInteractionFlags(qt.Qt.TextSelectableByKeyboard | qt.Qt.TextSelectableByMouse)
    
    def __init__(self,reflection, message, parent=None):
        qt.QDialog.__init__(self, parent=None)
        
        layout = qt.QVBoxLayout()
        
        layout.addWidget(qt.QLabel(message))
        
        layout.addWidget(QReflectionAnglesDialog.QSelectableLabel("Reflection: %s" % reflection['hkl']))
        
        reflectionLayout = qt.QGridLayout()
        
        available_data = [k.split('_')[0] for k in reflection.keys()]
        self.header = ""
        self.idx_dict = {}
        i = 1
        if 'xy' in available_data:
            for lbl in ['x', 'y']:
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel(lbl), 0, i)
                self.idx_dict[lbl] = i
                i += 1
                self.header += lbl + '\t'
        if 'angles' in available_data:
            for lbl in ['mu', 'delta', 'gamma', 'theta', 'chi', 'phi']:
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel(lbl), 0, i)
                self.idx_dict[lbl] = i
                i += 1
                self.header += lbl + '\t'
        if 'imageno' in available_data:
            for lbl in ['imageno']:
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel(lbl), 0, i)
                self.idx_dict[lbl] = i
                i += 1
                self.header += lbl + '\t'
        
        self.checkboxes = []
        self.datastr = ""
        number_refl = max([int(k.split('_')[-1]) for k in reflection.keys() if '_' in k])
        for i in range(1,number_refl+1):
            box = qt.QCheckBox()
            reflectionLayout.addWidget(box, i, 0)
            
            if not reflection.get('selectable_%s' % i, False):
                box.setEnabled(False)
                
            self.checkboxes.append(box)
            if 'xy_%s' % i in reflection.keys():
                xy = reflection['xy_%s' % i]
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % xy[0]), i, self.idx_dict['x'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % xy[1]), i, self.idx_dict['y'])
                self.datastr += "%.5f\t" % xy[0]
                self.datastr += "%.5f\t" % xy[1]
            elif 'xy' in available_data:
                self.datastr += "\t\t"
                
            if 'angles_%s' % i in reflection.keys():
                angles = np.rad2deg(reflection['angles_%s' % i])
                angles[3] *= -1. # to theta
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[0]), i, self.idx_dict['mu'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[1]), i, self.idx_dict['delta'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[2]), i, self.idx_dict['gamma'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[3]), i, self.idx_dict['theta'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[4]), i, self.idx_dict['chi'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[5]), i, self.idx_dict['phi'])
                for ang in angles:
                    self.datastr += "%.5f\t" % ang
            elif 'angles' in available_data:
                self.datastr += "\t\t\t\t\t\t"
                
            if 'imageno_%s' % i in reflection.keys():
                imageno = reflection['imageno_%s' % i]
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%s" % imageno), i, self.idx_dict['imageno'])
                self.datastr += "%s" % imageno
            elif 'imageno' in available_data:
                self.datastr += "\t"
            self.datastr += "\n"
        layout.addLayout(reflectionLayout)
                
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel,
                                      qt.Qt.Horizontal)
        copybutton = buttons.addButton("Copy to clipboad",qt.QDialogButtonBox.HelpRole)
        copybutton.clicked.connect(self.copy_to_clipboard)
        
        layout.addWidget(buttons)
        
        okbtn = buttons.button(qt.QDialogButtonBox.Ok)
        cancelbtn = buttons.button(qt.QDialogButtonBox.Cancel)
        
        okbtn.clicked.connect(self.accept)
        cancelbtn.clicked.connect(self.reject)

        self.setLayout(layout)
        
    def copy_to_clipboard(self):
        app = qt.QApplication.instance()
        app.clipboard().setText(self.header + "\n" + self.datastr)
        
        
        
class PeakReflDialog(qt.QDialog):

    
    def __init__(self, refl_selector):
        qt.QDialog.__init__(self, refl_selector)
        self.refl_selector = refl_selector
        self.candidate_reflections = []
        self.currentRefl = None
        self.idx = 0
        
        layout = qt.QVBoxLayout()
        
        toplayout = qt.QHBoxLayout()
        hkl = np.array([np.nan, np.nan, np.nan])
        self._hkl_label = qt.QLabel("New reflection: %s (%s/%s)" % (hkl, 0,0))
        toplayout.addWidget(self._hkl_label)
        
        
        self.refreshBtn = qt.QPushButton()
        self.refreshBtn.setIcon(icons.getQIcon('view-refresh'))
        self.refreshBtn.setToolTip("Refresh and reorder Bragg reflection list")
        self.refreshBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.refreshBtn.clicked.connect(self._onRefreshBragg)
        toplayout.addWidget(self.refreshBtn)
        
        self.gotoImageBtn = qt.QPushButton("search image")
        self.gotoImageBtn.setIcon(resources.getQicon('search-image'))
        self.gotoImageBtn.setToolTip("diplay image with the chosen reflection")
        self.gotoImageBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.gotoImageBtn.clicked.connect(self._onGotoImage)
        toplayout.addWidget(self.gotoImageBtn)
        
        self.next_btn = qt.QPushButton('next')
        self.next_btn.setIcon(icons.getQIcon('next'))
        self.next_btn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.next_btn.setToolTip("Goto next Bragg reflection")
        self.next_btn.clicked.connect(self._onNextBragg)
        toplayout.addWidget(self.next_btn)
        
        self.prev_btn = qt.QPushButton('previous')
        self.prev_btn.setIcon(icons.getQIcon('previous'))
        self.prev_btn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.prev_btn.setToolTip("Goto previous Bragg reflection")
        self.prev_btn.clicked.connect(self._onPrevBragg)
        toplayout.addWidget(self.prev_btn)
        
        layout.addLayout(toplayout)
        
        editlayout = qt.QHBoxLayout()
        
        editlayout.addWidget(qt.QLabel("Edit reflection: "))
        
        self.setImageBtn = qt.QPushButton("select Image")
        self.setImageBtn.setIcon(resources.getQicon('select-image'))
        self.setImageBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.setImageBtn.setToolTip("set this image for the chosen reflection")
        self.setImageBtn.clicked.connect(self._onSelectImg)
        editlayout.addWidget(self.setImageBtn)
        
        self.peakSearchBtn = qt.QPushButton("2D peak search")
        self.peakSearchBtn.setIcon(resources.getQicon('search-peak'))
        self.peakSearchBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.peakSearchBtn.setToolTip("perform peak search in the image")
        self.peakSearchBtn.clicked.connect(self._onPeakSearch)
        editlayout.addWidget(self.peakSearchBtn)
        
        editlayout.addWidget(qt.QLabel("Accept reflection:"))
        

        self.acceptReflBtn = qt.QPushButton('Add reflection')
        self.acceptReflBtn.setIcon(icons.getQIcon('selected'))
        self.acceptReflBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.acceptReflBtn.setToolTip("add reflection to the list of reference reflections")
        self.acceptReflBtn.clicked.connect(self._onAcceptRefl)
        editlayout.addWidget(self.acceptReflBtn)
        
        layout.addLayout(editlayout)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.onCancel)
        
        layout.addWidget(buttons)
        
        self.setLayout(layout)
        
    def _onGotoImage(self):
        if self.currentRefl is not None:
            try:
                self.refl_selector.setReflectionActive(self.currentRefl.identifier)
                self.refl_selector.sigQueryImageChange.emit(self.currentRefl.imageno)
                self.refl_selector.sigQueryCenterPlot.emit(self.currentRefl.xy)
            except Exception as e:
                print("Cannot goto reflection %s: %s" % (self.currentRefl, e))
                
        
    def _onNextBragg(self):
        if self.currentRefl is not None and self.candidate_reflections:
            if self.idx < len(self.candidate_reflections)-1:
                self.set_current_candidate(self.idx + 1)
            else:
                qt.QMessageBox.warning(self, "Max Bragg reflections reached", "The maximum of available Bragg reflections has been reached.")
    
    def _onPrevBragg(self):
        if self.currentRefl is not None and self.candidate_reflections:
            if self.idx > 0:
                self.set_current_candidate(self.idx - 1)
            else:
                qt.QMessageBox.warning(self, "Min Bragg reflections reached", "The lowest of available Bragg reflections has been reached.")

    def _onSelectImg(self):
        if self.currentRefl is not None and self.candidate_reflections:
            try:
                self.refl_selector.setReflectionActive(self.currentRefl.identifier)
                self.refl_selector._onSetImage()
            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot select image number:\n%s" % e, traceback.format_exc())
        
    def _onRefreshBragg(self):
        try:
            refl_list = self.refl_selector.getBraggCandidates(True)
            self.set_candidates(refl_list)
        except Exception as e:
            qutils.warning_detailed_message(self, "Error","Cannot refresh Bragg reflections:\n%s" % e, traceback.format_exc())
            
    def _onPeakSearch(self):
        if self.currentRefl is not None and self.candidate_reflections:
            try:
                self.refl_selector.setReflectionActive(self.currentRefl.identifier)
                self.refl_selector._onSinglePeakSearch()
            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot perform peak search:\n%s" % e, traceback.format_exc())
        
    def _onAcceptRefl(self):
        if self.currentRefl is not None and self.candidate_reflections:
            del self.candidate_reflections[self.idx]
            self.currentRefl = None
            if self.idx < len(self.candidate_reflections) -1:
                self.set_current_candidate(self.idx)
            elif self.candidate_reflections:
                self.set_current_candidate(self.idx-1)
            else:
                self.idx = 0
                hkl = np.array([np.nan, np.nan, np.nan])
                self._hkl_label.setText("New reflection: %s (%s/%s)" % (hkl, 0,0))
                
        
    def set_current_candidate(self, idx):
        refl = self.candidate_reflections[idx]
        if self.currentRefl is not None:
            try:
                self.refl_selector.deleteReflection(self.currentRefl.identifier)
            except Exception as e:
                print("Cannot remove reflection %s: %s" % (self.currentRefl, e))
        dd = {'x' : refl.xy[0], 'y' : refl.xy[1]}
        refl = self.refl_selector.addReflection(dd,refl.imageno,hkl=refl.hkl)
        self.candidate_reflections[idx] = refl
        self.currentRefl = refl
        self.idx = idx
        self._hkl_label.setText("New reflection: %s (%s/%s)" % (refl.hkl, idx+1, len(self.candidate_reflections)))
        self.refl_selector.sigQueryImageChange.emit(refl.imageno)
        self.refl_selector.sigQueryCenterPlot.emit(refl.xy)
        
        
    def set_candidates(self,refl_list):
        if refl_list:
            self.candidate_reflections = copy.deepcopy(refl_list)
            self.set_current_candidate(0)
        else:
            if self.currentRefl is not None:
                try:
                    self.refl_selector.deleteReflection(self.currentRefl.identifier)
                except Exception as e:
                    print("Cannot remove reflection %s: %s" % (self.currentRefl, e))
            self.idx = 0
            self.currentRefl = None
            self.candidate_reflections = []
            hkl = np.array([np.nan, np.nan, np.nan])
            self._hkl_label.setText("New reflection: %s (%s/%s)" % (hkl, 0,0))
            
        
    def onOk(self):
        self.set_candidates([])
        self.accept()

    def onCancel(self):
        self.set_candidates([])
        self.reject()
        
class PeakImgRangeDialog(qt.QDialog):

    
    def __init__(self,axis, axisname, parent=None):
        qt.QDialog.__init__(self, parent)
        
        layout = qt.QGridLayout()
        
        self._scan_label = qt.QLabel("%s-scan: from %s to %s, max image will be calculated relative to set peak position" % (axisname, axis[0], axis[-1]))
        layout.addWidget(self._scan_label,0,0, 1, -1)
        
        layout.addWidget(qt.QLabel("1st pass (rough alignment):"),1,0, 1, -1)
        
        layout.addWidget(qt.QLabel("max img range:"),2, 0)
        layout.addWidget(qt.QLabel("from (rel):"),2, 1)
        layout.addWidget(qt.QLabel("to (rel):"),2, 3)

        self.roi_from = qt.QDoubleSpinBox()
        self.roi_from.setRange(-1000000, 1000000)
        self.roi_from.setDecimals(2)
        self.roi_from.setSuffix(u" °")
        self.roi_from.setValue(-2.0)
        layout.addWidget(self.roi_from, 2,2)
        
        self.roi_to = qt.QDoubleSpinBox()
        self.roi_to.setRange(-1000000, 1000000)
        self.roi_to.setDecimals(2)
        self.roi_to.setSuffix(u" °")
        self.roi_to.setValue(2.0)
        layout.addWidget(self.roi_to, 2,4)
        
        
        layout.addWidget(qt.QLabel("ROI size:"),3, 0)
        layout.addWidget(qt.QLabel("size v:"),3, 1)
        layout.addWidget(qt.QLabel("size h:"),3, 3)

        self.vsize = qt.QSpinBox()
        self.vsize.setRange(1, 1000000)
        self.vsize.setSuffix(u" pix")
        self.vsize.setValue(100)
        layout.addWidget(self.vsize, 3,2)
        
        self.hsize = qt.QSpinBox()
        self.hsize.setRange(1, 1000000)
        self.hsize.setSuffix(u" pix")
        self.hsize.setValue(100)
        layout.addWidget(self.hsize, 3,4)
        
        
        
        layout.addWidget(qt.QLabel("2nd pass (fine alignment):"),4,0, 1, -1)
        
        layout.addWidget(qt.QLabel("max img range:"),5, 0)
        layout.addWidget(qt.QLabel("from (rel):"),5, 1)
        layout.addWidget(qt.QLabel("to (rel):"),5, 3)

        self.roi_from_2 = qt.QDoubleSpinBox()
        self.roi_from_2.setRange(-1000000, 1000000)
        self.roi_from_2.setDecimals(2)
        self.roi_from_2.setSuffix(u" °")
        self.roi_from_2.setValue(-1.0)
        layout.addWidget(self.roi_from_2, 5,2)
        
        self.roi_to_2 = qt.QDoubleSpinBox()
        self.roi_to_2.setRange(-1000000, 1000000)
        self.roi_to_2.setDecimals(2)
        self.roi_to_2.setSuffix(u" °")
        self.roi_to_2.setValue(1.0)
        layout.addWidget(self.roi_to_2, 5,4)
        
        
        layout.addWidget(qt.QLabel("ROI size:"),6, 0)
        layout.addWidget(qt.QLabel("size v:"),6, 1)
        layout.addWidget(qt.QLabel("size h:"),6, 3)

        self.vsize_2 = qt.QSpinBox()
        self.vsize_2.setRange(1, 1000000)
        self.vsize_2.setSuffix(u" pix")
        self.vsize_2.setValue(60)
        layout.addWidget(self.vsize_2, 6,2)
        
        self.hsize_2 = qt.QSpinBox()
        self.hsize_2.setRange(1, 1000000)
        self.hsize_2.setSuffix(u" pix")
        self.hsize_2.setValue(60)
        layout.addWidget(self.hsize_2, 6,4)
        
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,7,0,1,-1)
        
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)
        
        self.setLayout(layout)
        
    def set_axis(self,axis, axisname):
        self._scan_label.setText("%s-scan: from %s to %s, max image will be calculated relative to set peak position" % (axisname, axis[0], axis[-1]))
        
    def _verify_ranges(self):
        if self.roi_from.value() > self.roi_to.value():
            raise ValueError("Invalid input of max img range: from %s > to %s" % (self.roi_from.value(),self.roi_to.value()) )
        if self.roi_from_2.value() > self.roi_to_2.value():
            raise ValueError("Invalid input of max img range: from %s > to %s" % (self.roi_from.value(),self.roi_to.value()) )
        return True
        
    def get_roi(self):
        self._verify_ranges()
        ddict = {
            'from' : self.roi_from.value(),
            'to' : self.roi_to.value(),
            'vsize' : self.vsize.value(),
            'hsize' : self.hsize.value(),
            'from_2' : self.roi_from_2.value(),
            'to_2' : self.roi_to_2.value(),
            'vsize_2' : self.vsize_2.value(),
            'hsize_2' : self.hsize_2.value()
        }
        return ddict
        
        
    def onOk(self):
        try:
            self._verify_ranges()
        except Exception as e:
            qt.QMessageBox.warning(self, "Invalid input", str(e))
            return
        self.accept()
        
        
      

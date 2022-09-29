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

import sys
import traceback
import os
from io import StringIO
from silx.gui import qt
from .ArrayTableDialog import ArrayEditWidget
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import ClassVar

import numpy as np
from datautils.xrayutils import HKLVlieg

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

class QReflectionSelector(qt.QSplitter):
    sigQueryImageChange = qt.pyqtSignal(int)
    sigQuerySaveReflections = qt.pyqtSignal()
    sigQueryLoadReflections = qt.pyqtSignal()
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self,plot, ubcalc, parent=None):
        qt.QSplitter.__init__(self, parent=None)
        self.setOrientation(qt.Qt.Vertical)
        self.plot = plot
        self.ubcalc = ubcalc
        self.orparent = parent
        
        self.reflections = []
        self.reflBragg = []
        
        self.nextNo = 0
        self.imageno = 0
        
        self._showReferenceReflections = True
        self.showBraggReflections = False
        
        self.activeReflection = None
        
        self.plot.sigKeyPressDelete.connect(self._onDelete)
        
        #self.mainLayout = qt.QVBoxLayout()
        self.reflectionWidget = qt.QSplitter(self)
        self.reflectionWidget.setOrientation(qt.Qt.Horizontal)
        
        qt.QLabel("No:",self.reflectionWidget)
        self.markerNo = qt.QSpinBox(self.reflectionWidget)
        self.markerNo.setRange(0,0)
        
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

        self.actionbuttons = qt.QSplitter(self)
        self.actionbuttons.setOrientation(qt.Qt.Horizontal)
        
        applyButton = qt.QPushButton("Apply",self.actionbuttons)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        applyButton.setToolTip("Accept changes")
        applyButton.clicked.connect(self._onApplyChange)
        
        gotoImageButton = qt.QPushButton("goto Image",self.actionbuttons)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        gotoImageButton.setToolTip("diplay image with this reflection")
        gotoImageButton.clicked.connect(self._onGotoImage)
        
        setImageButton = qt.QPushButton("set Image",self.actionbuttons)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        setImageButton.setToolTip("set this image for this reflection")
        setImageButton.clicked.connect(self._onSetImage)
        
        deleteButton = qt.QPushButton("Delete",self.actionbuttons)
        #applyButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        deleteButton.setToolTip("Delete Reflection")
        deleteButton.clicked.connect(self._onDelete)
        
        self.addWidget(self.reflectionWidget)
        self.addWidget(self.actionbuttons)
        
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
        
        self.addWidget(editorTabWidget)
        
        
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
        self.markerNo.setMaximum(len(self.reflections))
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
            self._setactiveHKLVals(refl.hkl)
            self.markerNo.setValue(int(refl.identifier.split('_')[-1]))
        else:
            self.activeReflection = refl.identifier
            self._setactiveHKLVals(refl.hkl)
            self.markerNo.setValue(int(refl.identifier.split('_')[-1]))
        
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
            self.reflBragg.append(HKLReflection(yxr[::-1], hkl,1, identifier))
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
            qt.QMessageBox.critical(self,"Cannot read reflections", "Cannot read reeflections %s" % traceback.format_exc())
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
        self.setReflections(self.reflections)
        
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
        if self.nextNo > 0:
            self.markerNo.setMaximum(self.nextNo-1)
        else:
            self.markerNo.setMaximum(0)
        if self.reflections:
            newactive = self.reflections[0].identifier
            self.setReflectionActive(newactive)
        
    def addReflection(self,eventdict,imageno,hkl=np.array([np.nan,np.nan,np.nan])):
        identifier = 'ref_'+str(self.nextNo)
        refl = HKLReflection([eventdict['x'],eventdict['y']],hkl,imageno,identifier)
        if self._showReferenceReflections:
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        self.reflections.append(refl)
        self.markerNo.setMaximum(self.nextNo)
        self.nextNo += 1
        self.updateEditor()
        self.setReflectionActive(identifier)

        
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
            self._setactiveHKLVals(refl.hkl)
            self.markerNo.setValue(int(identifier.split('_')[-1]))
        else:
            self.activeReflection = identifier
            self._setactiveHKLVals(refl.hkl)
            self.markerNo.setValue(int(identifier.split('_')[-1]))
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
        
    def _setactiveHKLVals(self,hkl):
        nanto0 = lambda x : x if not np.isnan(x) else 0 
        self.Hbox.setValue(nanto0(hkl[0]))
        self.Kbox.setValue(nanto0(hkl[1]))
        self.Lbox.setValue(nanto0(hkl[2]))
        
    def _onApplyChange(self):
        if self.activeReflection is not None:
            hkl = np.empty(3)
            hkl[0] = self.Hbox.value()
            hkl[1] = self.Kbox.value()
            hkl[2] = self.Lbox.value()
            self.reflections[self.indexReflection(self.activeReflection)].hkl = hkl
            self.redrawActiveReflection()
            self.updateEditor()
            
    def redrawActiveReflection(self):
        if self._showReferenceReflections:
            self.plot.removeMarker(self.activeReflection)
            refl = self.reflections[self.indexReflection(self.activeReflection)]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(self.activeReflection),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
        
    def _onGotoImage(self):
        if self.activeReflection is not None:
            self.sigQueryImageChange.emit(self.reflections[self.indexReflection(self.activeReflection)].imageno)
            
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

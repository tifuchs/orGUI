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

import numpy as np

class HKLReflection(object):
    def __init__(self,xy,hkl,imageno):
        self.xy = np.array(xy,dtype=np.float)
        self.hkl = np.array(hkl,dtype=np.float)
        self.imageno = int(imageno)
        
    def __str__(self):
        x, y = self.xy
        h, k, l = self.hkl
        return "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f" % (h,k,l,x,y,self.imageno)
        
    @classmethod
    def fromStr(cls,string):
        h,k,l,x,y,imageno = string.split()
        refl = cls([x,y],[h,k,l],imageno)
        return refl
    
    @staticmethod
    def fromArray(arr):
        reflist = []
        if arr.ndim == 1:
            h,k,l,x,y,imageno = arr
            reflist.append(HKLReflection([x,y],[h,k,l],imageno))
            return reflist
        for row in arr:
            h,k,l,x,y,imageno = row
            reflist.append(HKLReflection([x,y],[h,k,l],imageno))
        return reflist
    
        

class QReflectionSelector(qt.QSplitter):
    sigQueryImageChange = qt.pyqtSignal(int)
    sigQuerySaveReflections = qt.pyqtSignal()
    sigQueryLoadReflections = qt.pyqtSignal()
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self,plot, parent=None):
        qt.QSplitter.__init__(self, parent=None)
        self.setOrientation(qt.Qt.Vertical)
        self.plot = plot
        
        self.reflections = {}
        self.reflBragg = {}
        
        self.nextNo = 0
        self.imageno = 0
        
        self._showReferenceReflections = True
        self.showBraggReflections = False
        
        self.activeReflection = None
        
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
        
        editorSplitter = qt.QSplitter()
        editorSplitter.setOrientation(qt.Qt.Horizontal)

        self.refleditor = qt.QTextEdit()
        self.updateEditor()
        editorSplitter.addWidget(self.refleditor)
        
        
        operations = qt.QSplitter()
        operations.setOrientation(qt.Qt.Vertical)
        
        
        fromEditorButton = qt.QPushButton("from editor")
        fromEditorButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Minimum, qt.QSizePolicy.Minimum))
        fromEditorButton.setToolTip("take refelctions from editor")
        fromEditorButton.clicked.connect(self.reflectionsFromEditor)
        
        saveReflButton = qt.QPushButton("save reflections")
        saveReflButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Minimum, qt.QSizePolicy.Minimum))
        saveReflButton.setToolTip("save all reflections")
        saveReflButton.clicked.connect(self._onSaveReflections)
        
        loadReflButton = qt.QPushButton("load reflections")
        loadReflButton.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Minimum, qt.QSizePolicy.Minimum))
        loadReflButton.setToolTip("load reflections from file")
        loadReflButton.clicked.connect(self._onLoadReflections)
        
        
        operations.addWidget(fromEditorButton)
        operations.addWidget(saveReflButton)
        operations.addWidget(loadReflButton)
        
        editorSplitter.addWidget(operations)
        
        
        self.addWidget(editorSplitter)
        
    def setBraggReflectionsVisible(self, visible):
        if visible:
            self.showBraggReflections = True
            self.redrawBraggReflections()
        else:
            self.showBraggReflections = False
            self.clearPlotBraggReflections()
        
    def setBraggReflections(self, hkls, yx, angles):
        for i, (hkl, yxr, pos) in enumerate(zip(hkls, yx, angles)):
            identifier = 'bragg_'+str(i)
            self.reflBragg[identifier] = HKLReflection(yxr[::-1], hkl,1)
        if self.showBraggReflections:
            self.redrawBraggReflections()
        else:
            self.clearPlotBraggReflections()
            
        
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
        s = "H\tK\tL\tx\ty\timageno\n"
        for k in self.reflections.keys():
            s += str(self.reflections[k]) + "\n"
        self.refleditor.setPlainText(s)
        
    def reflectionsFromEditor(self):
        plaintxt = self.refleditor.toPlainText()
        try:
            if len(plaintxt.splitlines()) > 1:
                txt = StringIO(self.refleditor.toPlainText())
                asarr = np.genfromtxt(txt,skip_header=1)
                if asarr.size == 0:
                    return
                reflist = HKLReflection.fromArray(asarr)
                self.setReflections(reflist)
            else:
                self.setReflections([])
        except Exception:
            qt.QMessageBox.critical(self,"Cannot read reflections", "Cannot read reeflections %s" % traceback.format_exc())
            return
            
    def clearPlotRefReflections(self):
        for refl in self.reflections:
            self.plot.removeMarker(refl)
            
    def clearPlotBraggReflections(self):
        for refl in self.reflBragg:
            self.plot.removeMarker(refl)
            
    def redrawBraggReflections(self):
        self.clearPlotBraggReflections()
        if self.showBraggReflections:
            for identifier in self.reflBragg:
                refl = self.reflBragg[identifier]
                #text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl)
                self.plot.addMarker(refl.xy[0],refl.xy[1],legend=identifier,color='g',selectable=False,draggable=False,symbol='x', text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl))
        
            
    def redrawRefReflections(self):
        self.setReflections(self.reflections.values())
        
    def setReflections(self,refls):
        self.clearPlotRefReflections()
        self.activeReflection = None
        self.nextNo = 0
        self.reflections = {}
        for refl in refls:
            identifier = 'ref_'+str(self.nextNo)
            if self._showReferenceReflections:
                self.plot.addMarker(refl.xy[0],refl.xy[1],legend=identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
            self.reflections[identifier] = refl
            self.nextNo += 1
        if self.nextNo > 0:
            self.markerNo.setMaximum(self.nextNo-1)
        else:
            self.markerNo.setMaximum(0)
        if self.reflections:
            newactive = next(iter(self.reflections.keys()))
            self.setReflectionActive(newactive)
        
    def addReflection(self,eventdict,imageno,hkl=np.array([np.nan,np.nan,np.nan])):
        identifier = 'ref_'+str(self.nextNo)
        if self._showReferenceReflections:
            self.plot.addMarker(eventdict['x'],eventdict['y'],legend=identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(hkl),color='blue',selectable=True,draggable=True,symbol='.')
        self.reflections[identifier] = HKLReflection([eventdict['x'],eventdict['y']],hkl,imageno)
        self.markerNo.setMaximum(self.nextNo)
        self.setReflectionActive(identifier)
        self.nextNo += 1
        self.updateEditor()
        
    def setReflectionActive(self,identifier):
        if self._showReferenceReflections:
            if self.activeReflection is not None:
                self.setReflectionInactive(self.activeReflection)
            self.plot.removeMarker(identifier)
            refl = self.reflections[identifier]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(identifier),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
            self.activeReflection = identifier
            self._setactiveHKLVals(refl.hkl)
            self.markerNo.setValue(int(identifier.split('_')[-1]))
        else:
            self.activeReflection = identifier
            refl = self.reflections[identifier]
            self._setactiveHKLVals(refl.hkl)
            self.markerNo.setValue(int(identifier.split('_')[-1]))
    
    def setReflectionInactive(self,identifier):
        if self._showReferenceReflections:
            self.plot.removeMarker(identifier)
            refl = self.reflections[identifier]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(identifier),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        
    def moveReflection(self,identifier,xy):
        self.reflections[identifier].xy = xy
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
            self.reflections[self.activeReflection].hkl = hkl
            self.redrawActiveReflection()
            self.updateEditor()
            
    def redrawActiveReflection(self):
        if self._showReferenceReflections:
            self.plot.removeMarker(self.activeReflection)
            refl = self.reflections[self.activeReflection]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(self.activeReflection),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
        
    def _onGotoImage(self):
        if self.activeReflection is not None:
            self.sigQueryImageChange.emit(self.reflections[self.activeReflection].imageno)
            
    def _onSetImage(self):
        if self.activeReflection is not None:
            self.reflections[self.activeReflection].imageno = self.imageno
            self.updateEditor()
        
    def getActiveReflection(self):
        pass
    
    def _onDelete(self):
        self.deleteReflection(self.activeReflection)
    
    def deleteReflection(self,identifier):
        if identifier not in self.reflections:
            return
        del self.reflections[identifier]
        if self._showReferenceReflections:
            self.plot.removeMarker(identifier)
        if self.activeReflection == identifier:
            self.activeReflection = None
            if self.reflections:
                newactive = next(iter(self.reflections.keys()))
                self.setReflectionActive(newactive)
        self.updateEditor()

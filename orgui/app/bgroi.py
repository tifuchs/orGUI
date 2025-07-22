# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2025 Timo Fuchs
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
__copyright__ = "Copyright 2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import logging
import numpy

from silx.gui import utils
from silx.gui.plot import items
from silx.gui.colors import rgba
from silx.image.shapes import Polygon
from silx.image._boundingbox import _BoundingBox
from silx.utils.proxy import docstring
from silx.gui.plot.utils.intersections import segments_intersection
from silx.gui.plot.items._roi_base import _RegionOfInterestBase

# He following imports have to be exposed by this module
from silx.gui.plot.items._roi_base import RegionOfInterest
from silx.gui.plot.items._roi_base import HandleBasedROI
from silx.gui.plot.items._arc_roi import ArcROI  # noqa
from silx.gui.plot.items._roi_base import InteractionModeMixIn  # noqa
from silx.gui.plot.items._roi_base import RoiInteractionMode  # noqa


logger = logging.getLogger(__name__)




class RectangleBgROI(HandleBasedROI, items.LineMixIn):
    """A ROI identifying a rectangle in a 2D plot.

    This ROI provides 1 anchor for each corner, plus an anchor in the
    center to translate the full ROI.
    """

    ICON = 'add-shape-rectangle'
    NAME = 'rectangle with bg ROI'
    SHORT_NAME = "rectangleBg"
    """Metadata for this kind of ROI"""

    _plotShape = "rectangle"
    """Plot shape which is used for the first interaction"""

    def __init__(self, parent=None):
        HandleBasedROI.__init__(self, parent=parent)
        items.LineMixIn.__init__(self)
        self._handleTopLeft = self.addHandle()
        self._handleTopRight = self.addHandle()
        self._handleBottomLeft = self.addHandle()
        self._handleBottomRight = self.addHandle()
        self._handleCenter = self.addTranslateHandle()
        self._handleLabel = self.addLabelHandle()
        
        self.__bgshapes = {}
        self.__bgsize = {}
        for n in ['left', 'right', 'top', 'bottom']:
            shape = items.Shape("rectangle")
            shape.setPoints([[0, 0], [0, 0]])
            shape.setFill(False)
            shape.setOverlay(True)
            shape.setLineStyle(self.getLineStyle())
            shape.setLineWidth(self.getLineWidth())
            shape.setColor(rgba(self.getColor()))
            self.__bgshapes[n] = shape
            self.__bgsize[n] = 0
            self.addItem(shape)
            
        self.centerMarker = items.Marker() #setPosition
        self.centerMarker.setPosition(0,0)
        self.addItem(self.centerMarker)
        
        shape = items.Shape("rectangle")
        shape.setPoints([[0, 0], [0, 0]])
        shape.setFill(False)
        shape.setOverlay(True)
        shape.setLineStyle(self.getLineStyle())
        shape.setLineWidth(self.getLineWidth())
        shape.setColor(rgba(self.getColor()))
        self.__shape = shape
        self.addItem(shape)


    def _updated(self, event=None, checkVisibility=True):
        if event in [items.ItemChangedType.VISIBLE]:
            self._updateItemProperty(event, self, [self.__shape, *self.__bgshapes.values(), self.centerMarker])
            #self._updateItemProperty(event, self, self.__shape)
        super(RectangleBgROI, self)._updated(event, checkVisibility)

    def _updatedStyle(self, event, style):
        super(RectangleBgROI, self)._updatedStyle(event, style)
        self.__shape.setColor(style.getColor())
        self.__shape.setLineStyle(style.getLineStyle())
        self.__shape.setLineWidth(style.getLineWidth())
        for n in self.__bgshapes:
            self.__bgshapes[n].setColor(style.getColor())
            self.__bgshapes[n].setLineStyle(style.getLineStyle())
            self.__bgshapes[n].setLineWidth(style.getLineWidth())

    def setFirstShapePoints(self, points):
        assert len(points) == 2
        self._setBound(points)

    def _setBound(self, points):
        """Initialize the rectangle from a bunch of points"""
        top = max(points[:, 1])
        bottom = min(points[:, 1])
        left = min(points[:, 0])
        right = max(points[:, 0])
        size = right - left, top - bottom
        self._updateGeometry(origin=(left, bottom), size=size)

    def _updateText(self, text):
        self._handleLabel.setText(text)
        
    def setBgSize(self, left=None, right=None, top=None, bottom=None):
        if left is not None:
            self.__bgsize['left'] = left
        if right is not None:
            self.__bgsize['right'] = right
        if top is not None:
            self.__bgsize['top'] = top
        if bottom is not None:
            self.__bgsize['bottom'] = bottom

        origin = self.getOrigin()
        size = self.getSize()
        self._updateGeometry(origin=origin, size=size)
        
    def _setBgSize(self, left=None, right=None, top=None, bottom=None):
        if left is not None:
            self.__bgsize['left'] = left
        if right is not None:
            self.__bgsize['right'] = right
        if top is not None:
            self.__bgsize['top'] = top
        if bottom is not None:
            self.__bgsize['bottom'] = bottom
        
    def setBgStyle(self, color, linestyle, linewidth):
        for n in self.__bgshapes:
            self.__bgshapes[n].setColor(color)
            self.__bgshapes[n].setLineStyle(linestyle)
            self.__bgshapes[n].setLineWidth(linewidth)
            
    def getCenter(self):
        """Returns the central point of this rectangle

        :rtype: numpy.ndarray([float,float])
        """
        pos = self._handleCenter.getPosition()
        return numpy.array(pos)

    def getOrigin(self):
        """Returns the corner point with the smaller coordinates

        :rtype: numpy.ndarray([float,float])
        """
        pos = self._handleBottomLeft.getPosition()
        return numpy.array(pos)

    def getSize(self):
        """Returns the size of this rectangle

        :rtype: numpy.ndarray([float,float])
        """
        vmin = self._handleBottomLeft.getPosition()
        vmax = self._handleTopRight.getPosition()
        vmin, vmax = numpy.array(vmin), numpy.array(vmax)
        return vmax - vmin

    def setOrigin(self, position):
        """Set the origin position of this ROI

        :param numpy.ndarray position: Location of the smaller corner of the ROI
        """
        size = self.getSize()
        self.setGeometry(origin=position, size=size)

    def setSize(self, size):
        """Set the size of this ROI

        :param numpy.ndarray size: Size of the center of the ROI
        """
        origin = self.getOrigin()
        self.setGeometry(origin=origin, size=size)

    def setCenter(self, position):
        """Set the size of this ROI

        :param numpy.ndarray position: Location of the center of the ROI
        """
        size = self.getSize()
        self.setGeometry(center=position, size=size)

    def setGeometry(self, origin=None, size=None, center=None, left=None, right=None, top=None, bottom=None):
        """Set the geometry of the ROI
        """
        if ((origin is None or numpy.array_equal(origin, self.getOrigin())) and
                (center is None or numpy.array_equal(center, self.getCenter())) and
                numpy.array_equal(size, self.getSize())):
                    
            if left is None and right is None and top is None and bottom is None:
                return  # Nothing has changed
        
        self._setBgSize(left, right, top, bottom)
        self._updateGeometry(origin, size, center)

    def _updateGeometry(self, origin=None, size=None, center=None):
        """Forced update of the geometry of the ROI"""
        if origin is not None:
            origin = numpy.array(origin)
            size = numpy.array(size)
            points = numpy.array([origin, origin + size])
            center = origin + size * 0.5
        elif center is not None:
            center = numpy.array(center)
            size = numpy.array(size)
            points = numpy.array([center - size * 0.5, center + size * 0.5])
        else:
            raise ValueError("Origin or center expected")

        with utils.blockSignals(self._handleBottomLeft):
            self._handleBottomLeft.setPosition(points[0, 0], points[0, 1])
        with utils.blockSignals(self._handleBottomRight):
            self._handleBottomRight.setPosition(points[1, 0], points[0, 1])
        with utils.blockSignals(self._handleTopLeft):
            self._handleTopLeft.setPosition(points[0, 0], points[1, 1])
        with utils.blockSignals(self._handleTopRight):
            self._handleTopRight.setPosition(points[1, 0], points[1, 1])
        with utils.blockSignals(self._handleCenter):
            self._handleCenter.setPosition(center[0], center[1])
        with utils.blockSignals(self._handleLabel):
            self._handleLabel.setPosition(points[0, 0], points[0, 1])
        
        # left:
        leftorigin = [points[0, 0] - self.__bgsize['left'], points[0, 1]]
        leftcornerur = [points[0, 0], points[1, 1]]
        self.__bgshapes['left'].setPoints(numpy.array([leftorigin, leftcornerur]))
        #right:
        rightorigin = [points[1, 0] , points[0, 1]]
        rightcornerur = [points[1, 0] + self.__bgsize['right'], points[1, 1]]
        self.__bgshapes['right'].setPoints(numpy.array([rightorigin, rightcornerur]))
        #top:
        toporigin = [points[0, 0] , points[1, 1] ]
        topcornerur = [points[1, 0] , points[1, 1] + self.__bgsize['top']]
        self.__bgshapes['top'].setPoints(numpy.array([toporigin, topcornerur]))
        #bottom:
        bottomorigin = [points[0, 0] , points[0, 1] - self.__bgsize['bottom'] ]
        bottomcornerur = [points[1, 0] , points[0, 1] ]
        self.__bgshapes['bottom'].setPoints(numpy.array([bottomorigin, bottomcornerur]))
        self.__shape.setPoints(points)
        
        self.centerMarker.setPosition(*center)
        
        self.sigRegionChanged.emit()

    @docstring(HandleBasedROI)
    def contains(self, position):
        assert isinstance(position, (tuple, list, numpy.array))
        points = self.__shape.getPoints()
        bb1 = _BoundingBox.from_points(points)
        return bb1.contains(position)

    def handleDragUpdated(self, handle, origin, previous, current):
        if handle is self._handleCenter:
            # It is the center anchor
            size = self.getSize()
            self._updateGeometry(center=current, size=size)
        else:
            opposed = {
                self._handleBottomLeft: self._handleTopRight,
                self._handleTopRight: self._handleBottomLeft,
                self._handleBottomRight: self._handleTopLeft,
                self._handleTopLeft: self._handleBottomRight,
            }
            handle2 = opposed[handle]
            current2 = handle2.getPosition()
            points = numpy.array([current, current2])

            # Switch handles if they were crossed by interaction
            if self._handleBottomLeft.getXPosition() > self._handleBottomRight.getXPosition():
                self._handleBottomLeft, self._handleBottomRight = self._handleBottomRight, self._handleBottomLeft

            if self._handleTopLeft.getXPosition() > self._handleTopRight.getXPosition():
                self._handleTopLeft, self._handleTopRight = self._handleTopRight, self._handleTopLeft

            if self._handleBottomLeft.getYPosition() > self._handleTopLeft.getYPosition():
                self._handleBottomLeft, self._handleTopLeft = self._handleTopLeft, self._handleBottomLeft

            if self._handleBottomRight.getYPosition() > self._handleTopRight.getYPosition():
                self._handleBottomRight, self._handleTopRight = self._handleTopRight, self._handleBottomRight

            self._setBound(points)

    def __str__(self):
        origin = self.getOrigin()
        w, h = self.getSize()
        params = origin[0], origin[1], w, h
        params = 'origin: %f %f; width: %f; height: %f' % params
        return "%s(%s)" % (self.__class__.__name__, params)
        
        
if __name__ == '__main__':
    from silx.gui import qt
    from silx.gui.plot import Plot2D
    from silx.gui.plot.tools.roi import RegionOfInterestManager
    app = qt.QApplication([])
    import numpy as np
    img = np.arange(100*100).reshape((100,100))
    pt = Plot2D()
    pt.addImage(img)
    
    roiManager = RegionOfInterestManager(pt)
    roiManager.setColor('pink')  # Set the color of ROI
    
    #self.roiTable = RegionOfInterestTableWidget()
    #self.roiTable.setRegionOfInterestManager(self.roiManager)
    
    #roi order: left, top, right, bottom,  croi
    rois = []

    roi = RectangleBgROI()
    
    roi.setGeometry(origin=(30, 30), size=(10, 10), left=5, right=6, top=7, bottom=8)

    roi.setColor('red') # bg color
    

    roi.setLineWidth(2)
    roi.setLineStyle('-')
    roi.setBgStyle('pink', '-', 2.)
    roi.setVisible(True)
    
    roi.setEditable(True)

    roiManager.addRoi(roi,useManagerColor=False)
    rois.append(roi)
    
    pt.show()
    app.exec()
    
    
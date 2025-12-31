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
__credits__ = []
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
from .. import __version__
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import copy
import sys
import os
import numpy as np


def cos_incidence_12(xy, sxrddetector):
    
    xy = np.asarray(xy)
    if len(xy.shape) > 1:
        x = xy.T[0]
        y = xy.T[1]
    else:
        x = np.asarray(xy[0]) 
        y = np.asarray(xy[1]) 
        
    p1, p2, p3 = sxrddetector._calc_cartesian_positions(x, y)
    if p3 is not None:
        raise NotImplementedError('Non planar detectors are not supported.')
        
    #R = np.array([
    #    [np.cos(sxrddetector._deltaChi), -np.sin(sxrddetector._deltaChi)],
    #    [np.sin(sxrddetector._deltaChi),  np.cos(sxrddetector._deltaChi)]
    #])
    #
    #p12 = np.vstack([p1, p2])
    #p12_rot = R @ p12
        
    cosa1 = sxrddetector.dist / np.sqrt(sxrddetector.dist * sxrddetector.dist + (p1 * p1 ))
    cosa2 = sxrddetector.dist / np.sqrt(sxrddetector.dist * sxrddetector.dist + (p2 * p2 ))
    
    return cosa1, cosa2
    
def parallax_correction(xy, sxrddetector, sizeBeamX, sizeBeamY):
    cosa1, cosa2 = cos_incidence_12(xy, sxrddetector)
    projected_X = sizeBeamX / cosa1
    projected_Y = sizeBeamY / cosa2
    return projected_X, projected_Y
    
    
def effective_sample_size(sxrddetector, sizeX, sizeY, sizeZ):
    """Rotate by azimuth and get max outline along x,z
    
    """
    
    points = np.array([
        [0, sizeZ, 0, sizeZ],
        [0, 0, sizeX, sizeX]
    ])
    
    R = np.array([
        [np.cos(sxrddetector._deltaChi), -np.sin(sxrddetector._deltaChi)],
        [np.sin(sxrddetector._deltaChi),  np.cos(sxrddetector._deltaChi)]
    ])
    
    points_rot = np.abs(R @ points)
    
    return np.amax(points_rot[0]), sizeY, np.amax(points_rot[1])
    
def projected_beamsize(xy, sxrddetector, sizeX, sizeY, sizeZ):
    esizeX, esizeY, esizeZ = effective_sample_size(sxrddetector, sizeX, sizeY, sizeZ)
    
    xy = np.asarray(xy)
    if len(xy.shape) > 1:
        x = xy.T[0]
        y = xy.T[1]
    else:
        x = np.asarray(xy[0]) 
        y = np.asarray(xy[1]) 
    
    # calculate delta and gamma in detector frame
    azimuth = sxrddetector.chi(x,y)
    tth = sxrddetector.tth(x,y)
    sintth = np.sin(tth); costth = np.cos(tth)
    delta_p = np.arctan2(sintth*np.sin(azimuth), costth)
    gamma_p = np.arctan2(sintth*np.cos(azimuth)*np.cos(delta_p), costth)
    
    # beam X Y refer to the detector (horiz, vert), esizeX, esizeY, esizeZ are in diffractometer (lab frame) coordinates
    beamX = np.abs(esizeY * np.sin(delta_p)) + np.abs(esizeX * np.cos(delta_p))
    beamY = np.abs(esizeY * np.sin(gamma_p)) + np.abs(esizeZ * np.cos(gamma_p))
    
    return beamX, beamY
    
    
def calc_corrections(xy, sxrddetector, 
                    roisize0=np.array([0,0]),
                    samplesize=None, 
                    parallax=True,
                    factor=1.):
    """samplesize must be a dict, sizes in m
    
    """
    if np.all(roisize0 == np.array([0,0])) and samplesize is None:
        raise ValueError("You must either provide the sample size or a minimum roisize")
    
    roisize_X_real = roisize0[0] * sxrddetector.detector.pixel2
    roisize_Y_real = roisize0[1] * sxrddetector.detector.pixel1
    if samplesize is not None:
        sizeX = samplesize['sizeX']
        sizeY = samplesize['sizeY']
        sizeZ = samplesize['sizeZ']
        
        beamX, beamY = projected_beamsize(xy, sxrddetector, sizeX, sizeY, sizeZ)
        
        beamX = np.maximum(beamX*factor, roisize_X_real)
        beamY = np.maximum(beamY*factor, roisize_Y_real)
    else:
        beamX = roisize_X_real
        beamY = roisize_Y_real
        
    if parallax:
        beamX, beamY = parallax_correction(xy, sxrddetector, beamX*factor, beamY*factor)
    
    roi_X = np.round(beamX / sxrddetector.detector.pixel2)
    roi_Y = np.round(beamY / sxrddetector.detector.pixel1)
    
    return np.vstack([roi_X, roi_Y]).T
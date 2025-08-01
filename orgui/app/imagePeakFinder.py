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
__copyright__ = "Copyright 2025 Timo Fuchs"
__license__ = "MIT License"
from .. import __version__
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import sys
import os
from silx.gui import qt
import warnings

import concurrent.futures
import queue
import threading
import numpy as np
import traceback
from scipy import ndimage


import numpy as np


def center_of_mass(img, mask=()):
    M = img[mask].sum()

    # y_indices is the row index (0,1,2...), x_indices is the column index (0,1,2...)
    y_indices, x_indices = np.indices(img.shape)

    # Compute weighted sums along each axis
    x_center = (img[mask] * x_indices[mask]).sum() / M
    y_center = (img[mask] * y_indices[mask]).sum() / M
    
    return x_center, y_center
    
    
def center_of_mass_1d(arr, mask=()):
    M = arr[mask].sum()
    indices = np.squeeze(np.indices(arr.shape))

    # Compute weighted sums
    com = (arr[mask] * indices[mask]).sum() / M
    
    return com


def roiCoords(yxcenter,yxsize,det_shape):
    ydetsize,xdetsize = det_shape
    ycenter, xcenter = yxcenter
    ysize, xsize = yxsize
    
    rangex = np.array([xcenter-xsize/2,xcenter+xsize/2])
    rangey = np.array([ycenter-ysize/2,ycenter+ysize/2])
    rangex = rangex.round(0)
    rangey = rangey.round(0)
    
    rangex = np.clip(rangex,0,xdetsize).astype(int)
    rangey = np.clip(rangey,0,ydetsize).astype(int)
    
    return rangey, rangex
    
    
def calc_image_range(fscan, axis_range, **kargs):
    max_workers = kargs.get('max_workers', 1)
    excluded_images = kargs.get('excluded_images', [])
    
    idx_1 = np.argmin(np.abs(fscan.axis - axis_range[0]))
    idx_2 = np.argmin(np.abs(fscan.axis - axis_range[1]))
    if idx_1 < idx_2:
        idx_min = idx_1; idx_max = idx_2
    else:
        idx_min = idx_2; idx_max = idx_1
    idx_min = max(0, idx_min)
    idx_max = min(idx_max, len(fscan)-1)
    
    image = fscan.get_raw_img(0)
    #global imgsum
    #global imgmax
    imgsum = np.zeros_like(image.img)
    imgmax = np.zeros_like(image.img)
    rocking_curve = np.full(len(list(range(idx_min, idx_max+1))), np.nan,dtype=float)
    rocking_axis = fscan.axis[idx_min: idx_max+1]
    lock = threading.Lock()
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor: # speedup only for the file reads 
        futures = {}
        def readfile_max(idx, imgno):
            if imgno in excluded_images: # skip if excluded
                return imgno
            image = fscan.get_raw_img(imgno) # here speedup during file read
            with lock:
                imgsum[:] += image.img
                np.maximum(imgmax, image.img, out=imgmax)
                if 'rocking_curve' in kargs:
                    ro_params = kargs['rocking_curve']
                    yxsizestart = np.array([ro_params['ysize'], ro_params['xsize']])
                    yxcenter = np.array(ro_params['xy_start'])[::-1]
                    yxcoords = roiCoords(yxcenter,yxsizestart,imgmax.shape)
                    subimage = np.copy(image.img[slice(*yxcoords[0]),slice(*yxcoords[1])])
                    if 'mask' in kargs:
                        submask = kargs['mask'][slice(*yxcoords[0]),slice(*yxcoords[1])]
                        subimage[submask] = 0
                    rocking_curve[idx] = np.nansum(subimage)
                
            return imgno
        for j, i in enumerate(range(idx_min, idx_max+1)):
            futures[executor.submit(readfile_max, j, i)] = i
        
        for f in concurrent.futures.as_completed(futures):
            try:
                imgno = f.result()
            except Exception as e:
                print("Cannot read image:\n%s" % traceback.format_exc())
                
    return {'max' : imgmax, 'sum' : imgsum, 'rocking_curve' : rocking_curve, 'rocking_axis' : rocking_axis}

    
def find_COM_Image(xy_start, xsize, ysize, fscan, axis_range, **kargs):
    kargs['rocking_curve'] = {
        'xy_start' : xy_start, 'xsize' : xsize, 'ysize' : ysize
    }
    images = calc_image_range(fscan, axis_range, **kargs)
    
    img_max = images['max']
    
    yxcenter = np.array(xy_start)[::-1]
    yxsizestart = np.array([ysize, xsize])
    
    deltapix = 1
    counter = 0
    
    #optimize peak position
    while(deltapix > 0.1 and counter < 10):
        yxcenter_old = np.copy(yxcenter)
        yxcoords = roiCoords(yxcenter,yxsizestart,img_max.shape)
        if 'mask' in kargs:
            submask = np.logical_not(kargs['mask'][slice(*yxcoords[0]),slice(*yxcoords[1])])
        else:
            submask = ()
        subimage = img_max[slice(*yxcoords[0]),slice(*yxcoords[1])]
        com = np.array(center_of_mass(subimage, submask))[::-1]
        
        yxcenter = com + np.array([yxcoords[0][0],yxcoords[1][0]])
        deltapix = np.linalg.norm(yxcenter - yxcenter_old) 
        counter += 1
    
    # find COM and peak of the rocking curve
    ro_curve = images['rocking_curve']
    ro_mask = np.logical_not(np.isnan(ro_curve))
    
    com_ro = center_of_mass_1d(ro_curve, ro_mask)
    com_idx = round(com_ro)
    axis_com = images['rocking_axis'][com_idx]
    
    peak_idx = np.nanargmax(ro_curve)
    axis_peak = images['rocking_axis'][peak_idx]
    
    return {'xy' : yxcenter[::-1], 'axis_com' : axis_com, 'axis_peak' : axis_peak} 
    
    
    
    
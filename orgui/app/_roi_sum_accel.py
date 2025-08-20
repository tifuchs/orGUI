# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
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


from numba import njit, config

import numpy as np

config.THREADING_LAYER = 'threadsafe'

@njit('f8[:,::1](f8[:,::1], b1[:, ::1], f8[:,::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], f8[:,::1])', nogil=True, cache=True)
def processImage(image,     mask,       C_corr,    croi,      leftroi,   rightroi,  toproi,   bottomroi, all_counters):

    for i in range(image.shape[0]): # this is faster than image[mask] = np.nan for some reason
        for j in range(image.shape[1]):
            if mask[i, j]:
                image[i, j] = np.nan
            else:
                image[i, j] *= C_corr[i, j]
    
    invmask = np.logical_not(mask)
    
    for i in range(croi.shape[0]):
        ckey = croi[i]
        leftkey = leftroi[i]
        rightkey = rightroi[i]
        topkey = toproi[i]
        bottomkey = bottomroi[i]
        
        # signal
        all_counters[i,0] = np.nansum(image[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        all_counters[i,1] = np.sum(invmask[ckey[1, 0]: ckey[1, 1], ckey[0, 0]: ckey[0, 1]])
        
        # background
        all_counters[i,2] = 0.
        all_counters[i,3] = 0.
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_counters[i, 2] += np.nansum(image[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
            all_counters[i, 3] += np.sum(invmask[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
    
    return all_counters

@njit('void(f8[:,::1], f8[:,::1], b1[:, ::1], f8[:,::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], f8[:,::1], f8[:,::1], f8[:,::1])', nogil=True, cache=True)
def processImage_bg_Carr(image,  bg,        mask,       C_corr,    croi,      leftroi,   rightroi,  toproi,   bottomroi, all_counters, all_Carr, all_Bgimg):
    invmask = np.logical_not(mask)
    """bgimage and C_corr mask must already be applied!!!
    
    
    """
    for i in range(image.shape[0]): # this is faster than image[mask] = np.nan for some reason
        for j in range(image.shape[1]):
            if mask[i, j]:
                image[i, j] = np.nan
                
    for i in range(croi.shape[0]):
        ckey = croi[i]
        leftkey = leftroi[i]
        rightkey = rightroi[i]
        topkey = toproi[i]
        bottomkey = bottomroi[i]
        
        # image signal
        all_counters[i,0] = np.nansum(image[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        
        # image background
        all_counters[i,2] = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_counters[i, 2] += np.nansum(image[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
        
        # Correction image
        all_Carr[i,0] = np.nansum(C_corr[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        all_Carr[i,2] = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_Carr[i, 2] += np.nansum(C_corr[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
        
        # Background image
        all_Bgimg[i,0] = np.nansum(bg[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        all_Bgimg[i,2] = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_Bgimg[i, 2] += np.nansum(bg[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
            
        # number of available pixels
        sig_pix = np.sum(invmask[ckey[1, 0]: ckey[1, 1], ckey[0, 0]: ckey[0, 1]])
        bg_pix = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            bg_pix += np.sum(invmask[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
            
        all_counters[i,1] = sig_pix
        all_counters[i, 3] = bg_pix
        all_Carr[i,1] = sig_pix
        all_Carr[i, 3] = bg_pix        
        all_Bgimg[i,1] = sig_pix
        all_Bgimg[i, 3] = bg_pix

@njit('void(f8[:,::1],  b1[:, ::1], f8[:,::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], f8[:,::1], f8[:,::1])', nogil=True, cache=True)
def processImage_Carr(image,        mask,       C_corr,    croi,      leftroi,   rightroi,  toproi,   bottomroi, all_counters, all_Carr):
    invmask = np.logical_not(mask)
    """bgimage and C_corr mask must already be applied!!!
    
    
    """
    for i in range(image.shape[0]): # this is faster than image[mask] = np.nan for some reason
        for j in range(image.shape[1]):
            if mask[i, j]:
                image[i, j] = np.nan
                
    for i in range(croi.shape[0]):
        ckey = croi[i]
        leftkey = leftroi[i]
        rightkey = rightroi[i]
        topkey = toproi[i]
        bottomkey = bottomroi[i]
        
        # image signal
        all_counters[i,0] = np.nansum(image[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        
        # image background
        all_counters[i,2] = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_counters[i, 2] += np.nansum(image[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
        
        # Correction image
        all_Carr[i,0] = np.nansum(C_corr[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        all_Carr[i,2] = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_Carr[i, 2] += np.nansum(C_corr[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
            
        # number of available pixels
        sig_pix = np.sum(invmask[ckey[1, 0]: ckey[1, 1], ckey[0, 0]: ckey[0, 1]])
        bg_pix = 0.0
        for key in [leftkey, rightkey, topkey, bottomkey]:
            bg_pix += np.sum(invmask[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
            
        all_counters[i,1] = sig_pix
        all_counters[i, 3] = bg_pix
        all_Carr[i,1] = sig_pix
        all_Carr[i, 3] = bg_pix        

@njit('f8[:,::1](f8[:,::1], f8[:,::1], b1[:, ::1], f8[:,::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], i8[:,:, ::1], f8[:,::1])', nogil=True, cache=True)
def processImage_bg(image,  bg,        mask,       C_corr,    croi,      leftroi,   rightroi,  toproi,   bottomroi, all_counters):
    image -= bg
    for i in range(image.shape[0]): # this is faster than image[mask] = np.nan for some reason
        for j in range(image.shape[1]):
            if mask[i, j]:
                image[i, j] = np.nan
            else:
                image[i, j] *= C_corr[i, j]
    
    invmask = np.logical_not(mask)
    
    for i in range(croi.shape[0]):
        ckey = croi[i]
        leftkey = leftroi[i]
        rightkey = rightroi[i]
        topkey = toproi[i]
        bottomkey = bottomroi[i]
        
        # signal
        all_counters[i,0] = np.nansum(image[ckey[1, 0]: ckey[1, 1] , ckey[0, 0]: ckey[0, 1]])
        all_counters[i,1] = np.sum(invmask[ckey[1, 0]: ckey[1, 1], ckey[0, 0]: ckey[0, 1]])
        
        # background
        all_counters[i,2] = 0.
        all_counters[i,3] = 0.
        for key in [leftkey, rightkey, topkey, bottomkey]:
            all_counters[i, 2] += np.nansum(image[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
            all_counters[i, 3] += np.sum(invmask[key[1, 0]:key[1, 1], key[0, 0]:key[0, 1]])
    
    return all_counters

@njit('void(f8[:,::1], f8[:,::1], f8[:,::1])', nogil=True, cache=True)
def calcMaxSum(image,     sumimg,    maximg):
    sumimg += image
    maximg[:] = np.maximum(image, maximg)
    
@njit('void(f8[:,::1], f8[:,::1])', nogil=True, cache=True)
def calcBgSub(image,     bg):
    image -= bg
    
@njit('void(f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1])', nogil=True, cache=True)
def calcMaxSum_bg(image,     sumimg,    maximg,  bg):
    imgbgsub = image - bg
    sumimg += imgbgsub
    maximg[:] = np.maximum(imgbgsub, maximg)
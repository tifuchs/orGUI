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


import numpy as np
import numpy.linalg as LA
from .. import util
import math as m
import scipy.optimize as opt
import pyFAI
from pyFAI import geometry
from pyFAI.utils import binning
#import pyFAI.azimuthalIntegrator
from pyFAI.ext import invert_geometry
import copy
import warnings

from .HKLVlieg import crystalAngles_singleArray, vacAngles_singleArray, primBeamAngles, vliegDiffracAngles

def load(ponifile):
    det = Detector2D_SXRD()
    det.load(ponifile)
    return det

def loadNXdict(nxdict):
    det = Detector2D_SXRD()
    det.fromNXdict(nxdict)
    return det
    

"""

Attention!!!!

There is something wierd with pyFAI

x,y seem to be inverted

"""
class Detector2D_SXRD(geometry.Geometry):
    
    def __init__(self,*args,**keyargs):
        super(Detector2D_SXRD,self).__init__(*args,**keyargs)
        self.setAzimuthalReference(0)
        self.setPolarization(0,0)
    
    def toNXdict(self):
        """To be used with silx.io.dictdump.dicttonx
        to save and load the data from any nexus file.
        
        """
        pyFAI_dict = self.getPyFAI()
        
        if hasattr(self, '_roi'):
            pyFAI_dict['poni1'] = self.poni1_max
            pyFAI_dict['poni2'] = self.poni2_max
            
        nxdict = {
            "detector_SXRD":
                {
                    "config" : {**pyFAI_dict},
                    "title" : u"pyFAI detector calibration",
                    "@NX_class": u"NXcollection",
                    "azimuth" : self._deltaChi,
                    "polarization" : self._polFactor,
                    "polarization_axis" : self._polAxis,
                    "binning" : np.array(self.detector.binning)
                },
            '@NX_class': u'NXcollection',
            "@creator" : "datautils v %s" % __version__
        }

        if hasattr(self, '_roi'):
            nxdict['detector_SXRD']["roi"] = np.array(self._roi)
        return nxdict

    def fromNXdict(self, nxdict):
        """Reads detector config from nexus dict written by toNXdict().
        Use silx.io.dictdump.nxtodict to read these nexus files.
        
        First searches for a "detector_SXRD" entry in the passed dict, otherwise
        assumes that the config is in the current group.
        """
        if "detector_SXRD" in nxdict:
            detdict = nxdict["detector_SXRD"]
        else:
            detdict = nxdict
        for entry in list(detdict['config']):
            if entry.startswith("@"):
                del detdict['config'][entry]
        detdict['config']['detector'] = str(detdict['config']['detector'])
        detdict['config']['max_shape'] = tuple(detdict['config']['max_shape'])
        self.setPyFAI(**detdict['config'])
        self.setAzimuthalReference(detdict['azimuth'])
        self.setPolarization(detdict['polarization_axis'], detdict['polarization'])
        if 'binning' in detdict:
            self.detector.set_binning(tuple(detdict['binning']))
        if 'roi' in detdict:
            roi = detdict['roi']
            self.set_roi(list(roi[0]), list(roi[1]))
        return self
        
                
    def copy_bin_roi_applied(self):
        """Returns a new Detector2D_SXRD instance, where the current binning and shape 
        have been set to the source image file. This is intended for further binning 
        of images which already have been processed with ``pply_bin_roi_to_image` once.
        """
        new_det = Detector2D_SXRD()
        ddict = self.toNXdict()
        new_det.fromNXdict(ddict)
        binns = new_det.detector.binning
        new_det.detector.max_shape = new_det.detector.shape
        new_det.poni1_max = new_det.poni1
        new_det.poni2_max = new_det.poni2
        new_det.detector.pixel1 *= binns[0]
        new_det.detector.pixel2 *= binns[1]
        new_det.detector.set_binning((1,1))
        new_det.set_roi([0,new_det.detector.shape[0]], [0,new_det.detector.shape[1]])
        new_det.reset()
        return new_det
        
    def __bin_img(self, dat):
        dim1, dim2 = self.detector.binning
        if dim1 == 1 and dim2 == 1:
            return dat
        t1 = dat.shape[0] % dim1
        t2 = dat.shape[1] % dim2
        if t1 == 0:
            t1 = None
        if t2 == 0:
            t2 = None
        dat_binned = binning(dat[:-t1, :-t2], (dim1, dim2))
        return dat_binned

    def set_roi(self, range1,  range2):
        self._roi = range1,  range2
        offset1 = range1[0]
        offset2 = range2[0]
        self.roi_range = range1, range2
        if not hasattr(self, 'poni1_max'):
            self.poni1_max = self.poni1
        if not hasattr(self, 'poni2_max'):
            self.poni2_max = self.poni2
        self.poni1 = self.poni1_max - offset1 * self.detector.pixel1
        self.poni2 = self.poni2_max - offset2 * self.detector.pixel2
        binns = self.detector.binning
        self.detector.set_binning((1,1))
        self.detector.set_binning(binns) # reset shape
        self.detector.shape = (range1[1] - range1[0], range2[1] - range2[0])
        self.reset()
        
    def apply_bin_roi_to_image(self, img):
        mask = np.isfinite(img).astype(np.int32)
        img_binned = self.__bin_img(np.nan_to_num(img))
        img_cnts = self.__bin_img(mask).astype(np.float64)
        img_cnts[~np.isfinite(img_cnts)] = np.nan
        img_binned *= (np.prod(self.detector.binning) / img_cnts)
        if hasattr(self, 'roi_range'):
            range1, range2 = self.roi_range
            return img_binned[range1[0]: range1[1], range2[0]: range2[1]]
        else:
            return img_binned

        
    def getBornAgain(self):
        """Only for surface facing upwards (along y), i.e. azimuth = 90°.
        """
        ba_config = {}
        
        ba_config['NbinsX'] = self.detector.shape[1]
        ba_config['Width'] = (self.detector.shape[1] * self.detector.pixel2) *1e3 # mm
        ba_config['NbinsY'] = self.detector.shape[0]
        ba_config['Height'] = (self.detector.shape[0] * self.detector.pixel1) *1e3 # mm
        ba_config['wavelength'] = self.get_wavelength()*1e9 # nm
        _cal = self.getFit2D()
        ba_config['Distance'] = _cal['directDist'] # mm
        ba_config['v0(vertBeamPos)'] = (self.detector.shape[0] - _cal['centerY']) * self.detector.pixel1 *1e3 # mm
        ba_config['u0(horizontBeamPos)'] = _cal['centerX'] * self.detector.pixel2 *1e3 # mm
        return ba_config
    
    #in rad!!!
    def setAzimuthalReference(self,deltaChi):
        self._deltaChi = copy.copy(deltaChi)
        self._cached_array = {}
        
    def getAzimuthalReference(self):
        return self._deltaChi
        
    def setPolarization(self,angle,factor):
        self._polAxis = angle
        self._polFactor = factor
        self._cached_array = {}
        
    def getPolarization(self):
        return self._polAxis, self._polFactor
        
    def primBeamAngles(self,shape=None):
        """gives angles in laboratory reference frame.
        
        """
        if self._cached_array.get("gamma_p") is None or self._cached_array.get("delta_p") is None:
            azimuth = self.chiArray(shape=shape) + self._deltaChi
            tth = self.array_from_unit(unit=pyFAI.units.TTH_RAD)
            sintth = np.sin(tth); costth = np.cos(tth)
            
            delta_p = np.arctan2(sintth*np.sin(azimuth), costth)
            gamma_p = np.arctan2(sintth*np.cos(azimuth)*np.cos(delta_p), costth)

            self._cached_array['gamma_p'] = gamma_p
            self._cached_array['delta_p'] = delta_p
        
        return self._cached_array['gamma_p'], self._cached_array['delta_p']
    
    #numpy array!!!
    def primBeamPoints(self,x,y):
        azimuth = self.chi(x,y) + self._deltaChi 
        tth = self.tth(x,y)
        sintth = np.sin(tth); costth = np.cos(tth)
        delta_p = np.arctan2(sintth*np.sin(azimuth), costth)
        gamma_p = np.arctan2(sintth*np.cos(azimuth)*np.cos(delta_p), costth)

        return gamma_p, delta_p
    
    def surfaceAngles(self,alpha_i,shape=None):
        """Angles in the reference frame, where the crystal is tilted by alpha_i.
        
        
        """
        if hasattr(self,'_alpha_i'):
            if self._alpha_i == alpha_i:
                if self._cached_array.get("gamma") is not None and self._cached_array.get("delta") is not None:
                    return self._cached_array.get("gamma"), self._cached_array.get("delta")
        self._alpha_i = alpha_i           
        gamma_p, delta_p = self.primBeamAngles(shape)
        gamma = np.arcsin( np.cos(alpha_i)*np.sin(gamma_p) - np.sin(alpha_i)*np.cos(delta_p)*np.cos(gamma_p) )
        delta = np.arcsin( (np.sin(delta_p)*np.cos(gamma_p))/np.cos(gamma) ) # evil, revise!!!
        self._cached_array["gamma"] = gamma
        self._cached_array["delta"] = delta
        return self._cached_array.get("gamma"), self._cached_array.get("delta")
    
    def crystalAngles(self,alpha_i,refraction_index,shape=None):
        if hasattr(self,'_n_refr'):
            if self._n_refr == refraction_index:
                if self._cached_array.get("gamma_cry") is not None and hasattr(self,'_alpha_i_ref'):
                    return self._cached_array.get("gamma_cry"), self._cached_array.get("delta"), self._alpha_i_ref
        self._n_refr = refraction_index
        gamma, delta = self.surfaceAngles(alpha_i,shape)
        self._cached_array["gamma_cry"] = crystalAngles_singleArray(gamma,self._n_refr)
        self._alpha_i_ref = crystalAngles_singleArray(self._alpha_i,self._n_refr)
        return self._cached_array.get("gamma_cry"), self._cached_array.get("delta"), self._alpha_i_ref
        
        
    def surfaceAnglesPoint(self,x,y,alpha_i):
        gamma_p, delta_p = self.primBeamPoints(x,y)
        gamma = np.arcsin( np.cos(alpha_i)*np.sin(gamma_p) - np.sin(alpha_i)*np.cos(delta_p)*np.cos(gamma_p) )
        delta = np.arcsin( (np.sin(delta_p)*np.cos(gamma_p))/np.cos(gamma) ) # evil, revise!!!
        return gamma, delta
    
    def crystalAnglesPoint(self,x,y,alpha_i,refraction_index):
        gamma, delta = self.surfaceAnglesPoint(x,y,alpha_i)
        gamma_cry = crystalAngles_singleArray(gamma,refraction_index)
        alpha_i_cry = crystalAngles_singleArray(alpha_i,refraction_index)
        return gamma_cry, delta, alpha_i_cry
        
    def pixelsTthChi(self,tth,chi):
        tth = np.atleast_1d(np.asarray(tth))
        chi = np.atleast_1d(np.asarray(chi))
        shape = tth.shape
        assert tth.shape == chi.shape
        
        # from here everything flat
        #tanchi = np.tan(chi.flatten())
        #nu = np.tan(tth.flatten()) * np.sin(chi.flatten())
        
        R = self.rotation_matrix() # detector rotation
        #M = np.empty((np.prod(shape), 2, 2))
        #A = np.empty((np.prod(shape), 2))
        chi = chi.flatten()
        tth = tth.flatten() 
        
        sinchi = np.sin(chi)
        coschi = np.cos(chi)
        
        sintth = np.sin(tth)
        costth = np.cos(tth)
        
        nu = sintth * sinchi
        
        R = self.rotation_matrix() # detector rotation
        
        a = R[0,0] * coschi - R[1,0] * sinchi
        b = R[0,1] * coschi - R[1,1] * sinchi
        c = R[0,0] * costth - R[2,0] * nu
        d = R[0,1] * costth - R[2,1] * nu
        
        A1 = (R[1,2] * sinchi - R[0,2] * coschi) * self.dist
        A2 = (R[2,2] * nu - R[0,2] * costth) * self.dist
        
        """
        The rest is a fast way to solve
        M @ ptilde = A
        with

        M = np.empty((np.prod(shape), 2, 2))
        A = np.empty((np.prod(shape), 2))
        
        A[:, 0] = A1
        A[:, 1] = A2
        
        M[:, 0, 0] = a
        M[:, 0, 1] = b
        M[:, 1, 0] = c
        M[:, 1, 1] = d
        
        linalg.solve doesn't work since singular matrices are present
        
        linalg.lstsq - solution: extremely slow!!!
        (> 10 s for Pilatus 2M pixels ~1600*1400 )
        
        ptilde = np.empty((np.prod(shape), 2))
        for i in range(a.size):
            ptilde[i] = LA.lstsq(M[i],A[i], rcond=None)[0]
        
        # Direct SVD:
        (about 10 s for Pilatus 2M pixels ~1600*1400)
        M_pinv = LA.pinv(M)
        ptilde_svd = np.empty((np.prod(shape), 2))
        for i in range(a.size):
            ptilde_svd[i] = M_pinv[i] @ A[i]
            
        so have to do it explicitly and find failing solutions manually:
        (about 1 s for Pilatus 2M pixels ~1600*1400)
        """
        determ = b*c - a*d
        
        ptilde1 = (b*A2 - d*A1) / determ # this raises divide by 0 warnings, not sure how to suppress them
        ptilde2 = (c*A1 - a*A2) / determ # this raises divide by 0 warnings, not sure how to suppress them
        
        mask = np.isclose(chi, 0., atol=1e-12) # solution fails at chi = 0
        # Special case chi = 0
        if mask.any():
            A = np.zeros((3,np.count_nonzero(mask)))
            A[1,:] = sintth[mask]
            A[2,:] = costth[mask]
            
            X = R.T @ A
            
            s = X[2] / self.dist
            
            ptilde1[mask] = X[0] / s
            ptilde2[mask] = X[1] / s
        
        
        p = (np.column_stack((ptilde1,ptilde2)) + np.array([self.poni1, self.poni2])) / np.array([self.pixel1, self.pixel2])
        
        return p.reshape((*shape, 2))
    

    
    @property
    def rangegamdel_p(self):
        """Is not identical to rangegamdel_p_full_det, where each pixel is calculated, but close within 1e-5
        
        Not sure why... 
        """
        xx = np.array([0.,self.detector.shape[0],0,self.detector.shape[0]]) #+ 0.5 # corner pixel center coordinates
        yy = np.array([0.,0., self.detector.shape[1],self.detector.shape[1]]) #+ 0.5
        gamma_p, delta_p = self.primBeamPoints(xx,yy)
        delrange = np.amin(delta_p), np.amax(delta_p)
        gamrange = np.amin(gamma_p), np.amax(gamma_p)
        return gamrange, delrange
    
    @property        
    def Qmax(self):
        """in A-1
        """
        xx = np.array([0.,self.detector.shape[0],0,self.detector.shape[0]]) #+ 0.5 # corner pixel center coordinates
        yy = np.array([0.,0., self.detector.shape[1],self.detector.shape[1]]) #+ 0.5
        Q = self.qFunction(xx,yy) / 10.
        return np.amax(Q)
        
    @property 
    def Qmin(self):
        """in A-1
        """
        qmin, qmax = self.Qrange
        return qmin
    
    @property    
    def Qrange(self):
        """in A-1
        """
        xx = np.array([0.,self.detector.shape[0],0,self.detector.shape[0]]) #+ 0.5 # corner pixel center coordinates
        yy = np.array([0.,0., self.detector.shape[1],self.detector.shape[1]]) #+ 0.5
        Q = self.qFunction(xx,yy) / 10.
        f2d_cal = self.getFit2D()
        # beam on detector ?
        if 0 <= f2d_cal['centerX'] <= self.detector.shape[1] and 0 <= f2d_cal['centerY'] <= self.detector.shape[0]:
            return 0., np.amax(Q)
        else: # fallback to numerical search
            edge_1_x = np.arange(self.detector.shape[0])
            edge_1_y = np.zeros(self.detector.shape[0])
            edge_2_x = edge_1_x
            edge_2_y = np.full(self.detector.shape[0], self.detector.shape[1] -1)
            edge_3_x = np.zeros(self.detector.shape[1]) 
            edge_3_y = np.arange(self.detector.shape[1])
            edge_4_x = np.full(self.detector.shape[1], self.detector.shape[0]-1)
            edge_4_y = edge_3_y
            #                   left edge right edge top edge  bottom edge
            xx = np.concatenate([edge_1_x, edge_2_x, edge_3_x, edge_4_x])
            yy = np.concatenate([edge_1_y, edge_2_y, edge_3_y, edge_4_y])
            Q = self.qFunction(xx,yy) / 10.
            return np.amin(Q), np.amax(Q) 
    
    def rangegamdel(self, alpha_i):
        """Is not identical to _rangegamdel_full_det, where each pixel is calculated, but close within 1e-5
        
        Not sure why... 
        """
        xx = np.array([0.,self.detector.shape[0],0,self.detector.shape[0]]) #+ 0.5 # corner pixel center coordinates
        yy = np.array([0.,0., self.detector.shape[1],self.detector.shape[1]]) #+ 0.5
        gamma, delta = self.surfaceAnglesPoint(xx,yy,alpha_i)
        delrange = np.amin(delta), np.amax(delta)
        gamrange = np.amin(gamma), np.amax(gamma)
        return gamrange, delrange
    

    @property
    def _rangegamdel_p_full_det(self):
        gamma_p, delta_p = self.primBeamAngles()
        delrange = np.amin(delta_p), np.amax(delta_p)
        gamrange = np.amin(gamma_p), np.amax(gamma_p)
        return gamrange, delrange
        
    def _rangegamdel_full_det(self, alpha_i):
        gamma_p, delta_p = self.surfaceAngles(alpha_i)
        delrange = np.amin(delta_p), np.amax(delta_p)
        gamrange = np.amin(gamma_p), np.amax(gamma_p)
        return gamrange, delrange

 
    
    def pixelsPrimeBeam(self,gamma_p,delta_p):
        gamma_p = np.atleast_1d(np.asarray(gamma_p))
        delta_p = np.atleast_1d(np.asarray(delta_p))
        shape = gamma_p.shape
        assert gamma_p.shape == delta_p.shape
        
        singam = np.sin(gamma_p); cosgam = np.cos(gamma_p)
        
        azimuth = np.arctan2(np.sin(delta_p) * cosgam, singam )
        tth = np.arctan2(singam, np.cos(delta_p) * cosgam * np.cos(azimuth))

        chi = azimuth - self._deltaChi
        
        return self.pixelsTthChi(tth, chi)
    
    def pixelsSurfaceAngles(self,gamma,delta,alpha_i):
        gamma_p = np.arcsin( np.cos(alpha_i)*np.sin(gamma) + np.sin(alpha_i)*np.cos(delta)*np.cos(gamma) )
        delta_p = np.arcsin( (np.sin(delta)*np.cos(gamma))/np.cos(gamma_p) )
        return self.pixelsPrimeBeam(gamma_p,delta_p)
        
    def pixelsCrystalAngles(self,gamma_cry,delta,alpha_i_cry, refraction_index):
        gamma = vacAngles_singleArray(gamma_cry,refraction_index)
        alpha_i = vacAngles_singleArray(alpha_i_cry,refraction_index)
        return self.pixelsSurfaceAngles(gamma,delta,alpha_i)
        
        
    
    # for now numerical solution, only iterables!
    # in rad
    def pixelsPrimPoint(self,gamma_p,delta_p,shape=None):
        warnings.warn("pixelsPrimPoint is deprecated and will be removed in the future. Use pixelsPrimeBeam instead.", FutureWarning)
        """
        gamma_p_det, delta_p_det = self.primBeamAngles(shape)
        ig = invert_geometry.InvertGeometry(gamma_p_det,delta_p_det)
        xy = []
        for gam_p, del_p in zip(gamma_p,delta_p):
            xy.append(ig(gam_p,del_p))
        return xy
        """
        return self.pixelsPrimeBeam(np.asarray(gamma_p),np.asarray(delta_p))
        
    
    def pixelsSurfacePoint(self,gamma,delta,alpha_i,shape=None):
        warnings.warn("pixelsSurfacePoint is deprecated and will be removed in the future. Use pixelsSurfaceAngles instead.", FutureWarning)
        """
        gamma_det, delta_det = self.surfaceAngles(alpha_i,shape)
        ig = invert_geometry.InvertGeometry(gamma_det,delta_det)
        xy = []
        for gam, delt in zip(gamma,delta):
            xy.append(ig(gam,delt))
        return xy
        """
        return self.pixelsSurfaceAngles(np.asarray(gamma),np.asarray(delta),np.asarray(alpha_i))
    
    def correctionArray(self,shape=None):
        if shape is None:
            shape = self.get_shape()
        if self._cached_array.get("corrarr") is not None:
            return self._cached_array.get("corrarr")
        else:
            self._cached_array["corrarr"] = self.solidAngleArray(shape)*self.polarization(shape=shape,factor=self._polFactor,axis_offset=self._polAxis)
        return self._cached_array.get("corrarr")
    
    
class DetectorCalibration():
    
    def __init__(self,E,calibrationCrystal,pixelsize):
        self._energy = E
        self._crystal = calibrationCrystal
        self._pixelsize = pixelsize

    def setE(self,E):
        self._energy = E

    # shapes: hkl = [[h1,k1,l1], [h2,k2,l2], ...] , xy = [[x1,y1], [x2,y2], ...]
    def calibrateFromReflections(self,hkl,xy):
        x_centralPixel = 730.
        xy = xy.T
        xyz = np.zeros((3,xy.shape[1]))
        xyz[:2] = xy
        
        twoTheta = self._crystal.get2ThetaFromHKL(hkl,self._energy)
        
        """
        p[0]: th_x, p[1]: th_y, p[2]: d_detector, p[3]: y_centralPixel
        """
        def Chi2(p):
            centralPixel = np.array([x_centralPixel,p[3],0.])
            _xyz = np.matrix((xyz.T - centralPixel).T)
            # uncomment to enable detector tilt correction:
            #_xyz = util.x_rotation(p[0])*util.y_rotation(p[1])*_xyz
            _xyz[2] = 0
            d = LA.norm(_xyz,axis=0)*self._pixelsize
            print(list((twoTheta - np.arctan(d/p[2])) / twoTheta ))
            return LA.norm((twoTheta - np.arctan(d/p[2]))/twoTheta)
        
        fitbounds = [(-0.1,0.1),(-0.1,0.1),(0.5,1.),(1400,1430)]
        res = opt.minimize(Chi2,np.array([0.0,0.0,0.72,1415.]),bounds=fitbounds,method='TNC')
        print(res)
        res = opt.minimize(Chi2,res.x)
        print(res)
        self.setCalibration([x_centralPixel,res.x[3]],res.x[2],[res.x[0],res.x[1]])
        self._rotz = None
        
    def setCalibration(self,centralPixel,distance,rotz=None):
        self._centralPixel = centralPixel
        self._distance = distance
        self._rotz = rotz
        
    def xyToDelGam(self,x,y):
        x = (x - self._centralPixel[0])*self._pixelsize
        y = (y - self._centralPixel[1])*-self._pixelsize
        #print(x)
        xy = np.vstack([x,y])
        if self._rotz is not None:
            mat = np.matrix([[np.cos(self._rotz),-np.sin(self._rotz)],
                             [np.sin(self._rotz),np.cos(self._rotz)]])
            #print(xy)
            xy = mat*xy
            xy = xy.A
            #print(xy)
            #raise NotImplementedError("Detector tilt is not implemented")
        
        delta = np.arctan(xy[0]/self._distance)
        gamma = np.arctan(xy[1]/self._distance)
        
        #print(delta)
        return delta , gamma
    
    def delGamToxy(self,delta,gamma):
        x = np.tan(delta)*self._distance
        y = np.tan(gamma)*self._distance
        
        x = (x/self._pixelsize) + self._centralPixel[0]
        y = -(y/self._pixelsize) + self._centralPixel[1]
        
        return [x,y]
    
    def xyToDelGam_corr(self,x,y):
        x = (x - self._centralPixel[0])*self._pixelsize
        y = (y - self._centralPixel[1])*-self._pixelsize
        #print(x)
        #xy = np.vstack([x,y])
        """
        if self._rotz is not None:
            mat = np.matrix([[np.cos(self._rotz),-np.sin(self._rotz)],
                             [np.sin(self._rotz),np.cos(self._rotz)]])
            #print(xy)
            xy = mat*xy
            xy = xy.A
            #print(xy)
            #raise NotImplementedError("Detector tilt is not implemented")
        """
        delta = np.arctan(x/self._distance)
        gamma = np.arctan(y/self._distance)
        
        deldelta = np.arctan((x+self._pixelsize)*(1/self._distance)) - delta
        delgamma = np.arctan((y+self._pixelsize)*(1/self._distance)) - gamma
        
        dd,dg = np.meshgrid(deldelta,delgamma)
                
        correction = dd*dg
        correction /= np.amax(correction)
        #print(delta)
        return delta , gamma, correction
    
    def xyToDelGam_grid(self,x,y):
        x = (np.array(x,dtype=np.float64) - self._centralPixel[0])*self._pixelsize
        y = (np.array(y,dtype=np.float64) - self._centralPixel[1])*-self._pixelsize
        #print(x)
            
        xx, yy = np.meshgrid(x,y)
        
        if self._rotz is not None:
            xx = np.cos(self._rotz)*xx - np.sin(self._rotz)*yy
            yy = np.sin(self._rotz)*xx + np.cos(self._rotz)*yy
        
        
        
        delta = np.arctan(xx*(1/self._distance))
        gamma = np.arctan(yy*(1/self._distance))
        
        deldelta = np.arctan((xx+self._pixelsize)*(1/self._distance)) - delta
        delgamma = np.arctan((yy+self._pixelsize)*(1/self._distance)) - gamma
        
        correction = deldelta*delgamma
        correction /= np.amax(correction)
        
        return delta , gamma, correction
    
    # **** Depricated  ****
    # Doesn't work good with just 4 points, but nice idea
    # P1 and P2 are one pair, P3, P4 the other, lower theta first
    def powderCalibration(self,P1,P2,P3,P4):
        hkl1 , pos1 = P1 
        hkl2 , pos2 = P2
        hkl3 , pos3 = P3
        hkl4 , pos4 = P4
        
        twoTheta1 = self._crystal.get2ThetaFromHKL(hkl1,self._energy)
        twoTheta2 = self._crystal.get2ThetaFromHKL(hkl2,self._energy)
        gamma = twoTheta2 -twoTheta1
        
        l1 = LA.norm(pos1 - pos2)*pixelsize
        l2 = LA.norm(pos3 - pos4)*pixelsize
        
        delta = util.solveTrigEquation(twoTheta1,twoTheta2,l1,l2)
        
        b1 = (l1/m.sin(gamma))*m.cos(delta - twoTheta2)
        
        A = (m.cos(delta - twoTheta1) / m.sin(twoTheta1) )**2 - 1.
        B = -2. * b1 *m.cos(delta - twoTheta1) * (1/m.tan(twoTheta1))
        print(b1)
        
        ld = - (B/(2*A)) + m.sqrt((B/(2*A))**2 - (b1**2)/A)
        print(m.sqrt((B/(2*A))**2 - (b1**2)/A))
        
        print("angle : %s, distance: %s pixel" % (np.rad2deg(delta),ld/pixelsize ))
        
        direction = (pos1 - pos2)/LA.norm(pos1 - pos2)
        print ("Central pixel: %s" % (direction*ld/pixelsize + pos1))
        #print("prev: %s" % LA.norm(pos1 - centralPixel))

    

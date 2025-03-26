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
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import unittest

from .. import DetectorCalibration
import pyFAI

import numpy as np
import os
import warnings

try:
    from silx.io import dictdump
    silx_avail = True
except:
    silx_avail = False


class TestRWDetector2D_SXRD(unittest.TestCase):
    
    def setUp(self):
        self.sxrddet = DetectorCalibration.Detector2D_SXRD()
        self.sxrddet.detector = pyFAI.detector_factory('Pilatus 2M CdTe')
        self.sxrddet.poni1 = 0.1
        self.sxrddet.poni2 = 0.1
        self.sxrddet.rot1 = np.pi/20
        self.sxrddet.rot2 = -np.pi/20
        self.sxrddet.rot3 = 0
        self.sxrddet.dist = 1.5
        self.sxrddet.set_energy(15.)
        self.sxrddet.setAzimuthalReference(np.deg2rad(90.))
        self.sxrddet.setPolarization(np.deg2rad(90.), 0.75)

    def test_from_to_dict(self):
        nxdict = self.sxrddet.toNXdict()
        othersxrddet = DetectorCalibration.Detector2D_SXRD()
        othersxrddet.fromNXdict(nxdict)
        self.check_sxrd_equal(othersxrddet)
        
    def check_sxrd_equal(self, other):
        self.assertDictEqual(self.sxrddet.getPyFAI(), other.getPyFAI())
        self.assertEqual(self.sxrddet._polFactor, other._polFactor)
        self.assertEqual(self.sxrddet._polAxis, other._polAxis)
        self.assertEqual(self.sxrddet._deltaChi, other._deltaChi)
        
    @unittest.skipUnless(silx_avail, "silx not available")
    def test_write_nx(self):
        self._nxfilename = "./detcal_test.nx"
        self.addCleanup(self._destruct_file, self._nxfilename)
        
        nxdict = self.sxrddet.toNXdict()
        dictdump.dicttonx(nxdict, self._nxfilename)
        self.assertTrue(os.path.isfile(self._nxfilename))
        os.remove(self._nxfilename)

    @unittest.skipUnless(silx_avail, "silx not available")
    def test_read_write_nx(self):
        self._nxfilename = "./detcal_test.nx"
        self.addCleanup(self._destruct_file, self._nxfilename)
        
        nxdict = self.sxrddet.toNXdict()
        dictdump.dicttonx(nxdict, self._nxfilename)

        read_dict = dictdump.nxtodict(self._nxfilename)
        othersxrddet = DetectorCalibration.Detector2D_SXRD()
        othersxrddet.fromNXdict(read_dict)
        self.check_sxrd_equal(othersxrddet)
        
        othersxrddet = DetectorCalibration.loadNXdict(read_dict)
        self.check_sxrd_equal(othersxrddet)
        
        os.remove(self._nxfilename)
    
        
    def _destruct_file(self, filename):
        if os.path.exists(filename):
            os.remove(self._nxfilename)
        
    
class TestAnglePixelConversion(unittest.TestCase):

    def setUp(self):
        self.sxrddet = DetectorCalibration.Detector2D_SXRD()
        self.sxrddet.detector = pyFAI.detector_factory('Pilatus 2M CdTe')
        self.sxrddet.poni1 = 0.1
        self.sxrddet.poni2 = 0.1
        self.sxrddet.rot1 = np.pi/20
        self.sxrddet.rot2 = -np.pi/20
        self.sxrddet.rot3 = 0
        self.sxrddet.dist = 1.5
        self.sxrddet.set_energy(15.)
        self.sxrddet.setAzimuthalReference(np.deg2rad(90.))
        self.sxrddet.setPolarization(np.deg2rad(90.), 0.75)
        
        self.p1 = np.arange(self.sxrddet.detector.shape[1] ) + 0.5 # pixel center
        self.p2 = np.arange(self.sxrddet.detector.shape[0] ) + 0.5
        self.p12 = np.moveaxis(np.array(np.meshgrid(self.p1,self.p2)),0, -1)[:,:,::-1] 
        
        self.mu = np.deg2rad(0.1) # should be differed, but probably ok
        
    
    def assertPixelErrorSurfaceAnglesLessThan(self, abserr=1e-3, msg=""):
        gamma, delta = self.sxrddet.surfaceAngles(self.mu)
        p12_conv = self.sxrddet.pixelsSurfaceAngles(gamma, delta, self.mu)
        maxerror = np.nanmax(np.abs(self.p12 - p12_conv))
        if maxerror > abserr:
            warnings.warn("too large error: %.5f pixel coord from surface angles, %s" % (maxerror, msg))
        #self.assertLessEqual(maxerror, abserr,  "too large error pixel coord from surface angles, %s" % msg)
        
    def assertPixelErrorTthChiLessThan(self, abserr=1e-3, msg=""):
        tth = self.sxrddet.twoThetaArray()
        chi = self.sxrddet.chiArray()
        p12_conv = self.sxrddet.pixelsTthChi(tth, chi)
        maxerror = np.nanmax(np.abs(self.p12 - p12_conv))
        if maxerror > abserr:
            warnings.warn("too large error: %.5f pixel coord from tth and chi, %s" % (maxerror, msg))
        #self.assertLessEqual(maxerror, abserr, "too large error pixel coord from tth and chi, %s" % msg)
        
    def assertGamDelRangeErrorLessThan(self, abserr=1e-3, msg=""):
        exact = self.sxrddet._rangegamdel_p_full_det
        corner = self.sxrddet.rangegamdel_p

        diff =  np.array(exact) - np.array(corner)
        max_rel_diff = np.amax(np.abs(diff))
        if max_rel_diff > abserr:
            warnings.warn("too large error: %.5f approx det corners, %s" % (max_rel_diff, msg))
        #self.assertLessEqual(max_rel_diff, abserr, "too large error approx det corners, %s" % msg)
    
    def assertQrangeValid(self):
        Q = self.sxrddet.qArray() / 10.
        Qmin = np.amin(Q)
        Qmax = np.amax(Q)
        Qmin_fast , Qmax_fast = self.sxrddet.Qrange
        f2d_cal = self.sxrddet.getFit2D()
        # beam on detector ?
        if 0 <= f2d_cal['centerX'] <= self.sxrddet.detector.shape[1] and 0 <= f2d_cal['centerY'] <= self.sxrddet.detector.shape[0]:
            Qmin = 0.
        
        self.assertTrue(np.allclose(Qmin, Qmin_fast, 1e-5))
        self.assertTrue(np.allclose(Qmax, Qmax_fast, 1e-5))

    def test_poni1(self):
        for p1 in np.linspace(-3,3,5):
            self.sxrddet.poni1 = p1
            msg="poni1 = %s" % p1
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.poni1 = 0.1
    
    def test_poni2(self):
        for p1 in np.linspace(-3,3,5):
            self.sxrddet.poni2 = p1
            msg="poni2 = %s" % p1
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.poni2 = 0.1
    
    def test_rot1(self):
        for p1 in np.linspace(0,np.pi,8):
            self.sxrddet.rot1 = p1
            msg="rot1 = %s" % p1
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.rot1 = np.pi/20
    
    def test_rot2(self):
        for p1 in np.linspace(0,np.pi,8):
            self.sxrddet.rot2 = p1
            msg="rot2 = %s" % p1
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.rot2 = -np.pi/20
        
    def test_rot3(self):
        for p1 in np.linspace(0,np.pi,8):
            self.sxrddet.rot3 = p1
            msg="rot3 = %s" % p1
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.rot3 = 0
        
    def test_dist(self):
        for d in np.linspace(0.01,10,5):
            self.sxrddet.dist = d
            msg="dist = %s" % d
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.dist = 1.5
        

        
        
"""
def test_del_gam_range():
    sxrddet = Detector2D_SXRD()
    sxrddet.detector = pyFAI.detector_factory('Pilatus 2M CdTe')
    sxrddet.poni1 = 0.1
    sxrddet.poni2 = 0.1
    sxrddet.rot1 = np.pi/4
    sxrddet.rot2 = np.pi/5
    sxrddet.rot3 = np.pi/6
    sxrddet.dist = 1.5
    sxrddet.setAzimuthalReference(np.deg2rad(90.))
    

    def checkRange():
        exact = sxrddet._rangegamdel_p_full_det
        corner = sxrddet.rangegamdel_p
        #print("Exact: ", exact)
        #print("Corner: ", corner)
        diff =  np.array(exact) - np.array(corner)
        max_rel_diff = np.amax(np.abs(diff))
        return max_rel_diff, diff
        #print("Difference:", np.array(exact) - np.array(corner))
    
    sxrddet.poni2 = 0.0
    for p1 in np.linspace(-5,5,20):
        sxrddet.poni1 = p1
        maxerr, diff = checkRange()
        print(f"Max numerical error at poni1 = {p1} m: {maxerr}")
        
    sxrddet.poni1 = 0.0
    for p2 in np.linspace(-5,5,20):
        sxrddet.poni2 = p2
        maxerr, diff = checkRange()
        print(f"Max numerical error at poni2 = {p2} m: {maxerr}")

    # test rot1:
    for r1 in np.linspace(0,np.pi,8):
        sxrddet.rot1 = r1
        maxerr, diff = checkRange()
        print(f"Max numerical error at rot1 = {np.rad2deg(r1)} deg: {maxerr}")
    sxrddet.rot1 = np.pi/4
    for r2 in np.linspace(0,np.pi,8):
        sxrddet.rot2 = r2
        maxerr, diff = checkRange()
        print(f"Max numerical error at rot2 = {np.rad2deg(r2)} deg: {maxerr}")
    
    sxrddet.rot2 = np.pi/4
    for r3 in np.linspace(0,np.pi,8):
        sxrddet.rot3 = r3
        maxerr, diff = checkRange()
        print(f"Max numerical error at rot3 = {np.rad2deg(r3)} deg: {maxerr}")
    sxrddet.rot3 = np.pi/4
    
    for d in np.linspace(0.001,10,5):
        sxrddet.dist = d
        maxerr, diff = checkRange()
        print(f"Max numerical error at dist = {d} m: {maxerr}")
    
    
    def testPixelConversion():
    sxrddet = Detector2D_SXRD()
    sxrddet.detector = pyFAI.detector_factory('Pilatus 2M CdTe')
    sxrddet.poni1 = 0.1
    sxrddet.poni2 = 0.1
    sxrddet.rot1 = np.pi/4
    sxrddet.rot2 = np.pi/5
    sxrddet.rot3 = np.pi/6
    sxrddet.dist = 1.5
    sxrddet.setAzimuthalReference(np.deg2rad(90.))
    
    # pixel coordinates:
    p1 = np.arange(sxrddet.detector.shape[1] ) + 0.5 # pixel center
    p2 = np.arange(sxrddet.detector.shape[0] ) + 0.5
    p12 = np.moveaxis(np.array(np.meshgrid(p1,p2)),0, -1)[:,:,::-1] 
    # this seems to be overcomplicated... is there a better method?
    
    
    def checkPixel_sxrd():
        gamma, delta = sxrddet.surfaceAngles(np.deg2rad(0.1))
        p12_conv = sxrddet.pixelsSurfaceAngles(gamma, delta, np.deg2rad(0.1))
        return np.nanmax(np.abs(p12 - p12_conv))
        #return np.allclose(p12, p12_conv, atol=1e-5)
        
    def checkPixel_tth():
        tth = sxrddet.twoThetaArray()
        chi = sxrddet.chiArray()
        p12_conv = sxrddet.pixelsTthChi(tth, chi)
        return np.nanmax(np.abs(p12 - p12_conv))
        
    checkPixel = checkPixel_sxrd
    
    sxrddet.poni2 = 0.0
    for p1 in np.linspace(-5,5,20):
        sxrddet.poni1 = p1
        maxerr = checkPixel()
        print(f"Max numerical error at poni1 = {p1} m: {maxerr}")
    sxrddet.poni1 = 0.0
    for p2 in np.linspace(-5,5,20):
        sxrddet.poni2 = p2
        maxerr = checkPixel()
        print(f"Max numerical error at poni2 = {p2} m: {maxerr}")

    # test rot1:
    for r1 in np.linspace(0,np.pi,8):
        sxrddet.rot1 = r1
        maxerr = checkPixel()
        print(f"Max numerical error at rot1 = {np.rad2deg(r1)} deg: {maxerr}")
    sxrddet.rot1 = np.pi/4
    for r2 in np.linspace(0,np.pi,8):
        sxrddet.rot2 = r2
        maxerr = checkPixel()
        print(f"Max numerical error at rot2 = {np.rad2deg(r2)} deg: {maxerr}")
    
    sxrddet.rot2 = np.pi/4
    for r3 in np.linspace(0,np.pi,8):
        sxrddet.rot3 = r3
        maxerr = checkPixel()
        print(f"Max numerical error at rot3 = {np.rad2deg(r3)} deg: {maxerr}")
    sxrddet.rot3 = np.pi/4
    
    for d in np.linspace(0.001,10,5):
        sxrddet.dist = d
        maxerr = checkPixel()
        print(f"Max numerical error at dist = {d} m: {maxerr}")
           
    
"""

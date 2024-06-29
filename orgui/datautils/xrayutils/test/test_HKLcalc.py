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

from .. import HKLVlieg
from ... import util

import numpy as np
import os


class TestLattice(unittest.TestCase):

    def testCreateLattice(self):
        lat = HKLVlieg.Lattice([3.9242, 3.9242, 3.9242] ,[90.0000 ,90.0000, 120.0000])
        self.assertTrue(np.allclose([3.9242, 3.9242, 3.9242], lat.a))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 120.0000]), lat.alpha))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 60.0000]), lat.beta))
        self.assertTrue(np.allclose([1.84883478, 1.84883478, 1.60113789], lat.b))
        
        lat.setLattice([3.9242, 3.9242, 3.9242] ,[90.0000 ,90.0000, 120.0000])
        self.assertTrue(np.allclose([3.9242, 3.9242, 3.9242], lat.a))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 120.0000]), lat.alpha))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 60.0000]), lat.beta))
        self.assertTrue(np.allclose([1.84883478, 1.84883478, 1.60113789], lat.b))
        
        lat.a = [3.9242, 3.9242, 3.9242]
        self.assertTrue(np.allclose([3.9242, 3.9242, 3.9242], lat.a))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 120.0000]), lat.alpha))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 60.0000]), lat.beta))
        self.assertTrue(np.allclose([1.84883478, 1.84883478, 1.60113789], lat.b))
        
        lat.alpha = np.deg2rad([90.0000 ,90.0000, 120.0000])
        self.assertTrue(np.allclose([3.9242, 3.9242, 3.9242], lat.a))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 120.0000]), lat.alpha))
        self.assertTrue(np.allclose(np.deg2rad([90.0000 ,90.0000, 60.0000]), lat.beta))
        self.assertTrue(np.allclose([1.84883478, 1.84883478, 1.60113789], lat.b))
        
        
    def testIsReciprocal(self):
        
        lat = HKLVlieg.Lattice([3.9242, 6.9242, 2.9242] ,[100.0000 ,90.0000, 120.0000])
        
        xyz_rel = (np.arange(100*3).reshape((3,100)) - 150) / 10.
        hkl = (np.arange(100*3).reshape((3,100)) - 150) / 10.
        
        phase = 2*np.pi * np.sum(xyz_rel*hkl,axis=0)
        
        xyz_cart = lat.directVectorCart(xyz_rel.T) 
        self.assertTrue(np.allclose(xyz_cart, lat.R_mat @ xyz_rel))

        Q_cart = lat.reciprocalVectorCart(hkl.T) 
        self.assertTrue(np.allclose(Q_cart, lat.B_mat @ hkl))
        
        phase_lat = np.sum(xyz_cart*Q_cart,axis=0)
        
        self.assertTrue(np.allclose(phase_lat, phase))



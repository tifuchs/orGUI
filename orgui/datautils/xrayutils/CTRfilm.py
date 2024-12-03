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
__version__ = "1.2.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

from .HKLVlieg import Lattice, UBCalculator, VliegAngles, GeometryCorrection
import numpy as np
import numpy.linalg as LA
from silx.io import dictdump
from .. import util
import xraydb
import warnings
import random
import json
import math
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.ndimage import gaussian_filter1d
import os
import re
from scipy.special import erf
#random.seed(45)
from collections import OrderedDict
import dataclasses
from dataclasses import dataclass ,field
import typing
import enum
import copy
import errno
import glob
from scipy.stats import skellam
from .element_data import cov_radii_array, rgb_array


from .CTRutil import (special_elementcolors, ParameterType, Parameter,
                     _ensure_contiguous, next_skip_comment, DWtoDisorder,
                     readWaasmaier, readDispersion, atomic_number,
                     estimateDispersionCompound)

from .CTRuc import UnitCell, HAS_NUMBA_ACCEL

class EpitaxyInterface():
    parameterOrder = "Width/cells Skew/cells"
    
    parameterLookup = {'W' : 1, 'S' : 2}
    
    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))
    
    def __init__(self,uc_layers_top, uc_layers_bottom, type='skellam' ,**kwargs):
        """
        
        uc_layers_top and uc_layers_bottom must be a list of UnitCell.
        Each UnitCell in uc_layers_top must have the same lattice constants
        and may be used to describe different layers in the unit cell. 
        Inteface cells are created taking into account the layers.

        
        """
        
        self._dirty = True
        self.type = type
        self.sigma_calc = kwargs.get('sigma_calc', 4)
        if isinstance(uc_layers_top, UnitCell):
            self.uc_layers_top = [uc_layers_top]
            self.uc_layers_bottom = [uc_layers_bottom]
        else:
            self.uc_layers_top = uc_layers_top
            self.uc_layers_bottom = uc_layers_bottom
        if len(self.uc_layers_top) != len(self.uc_layers_bottom):
            raise ValueError("Number of top and bottom layers must be equal.")
        a_top = self.uc_layers_top[0].a
        a_bot = self.uc_layers_bottom[0].a
        for uc in self.uc_layers_top:
            if np.any(uc.a != a_top):
                raise ValueError("Lattice constants of top layers must be equal.")
        for uc in self.uc_layers_bottom:
            if np.any(uc.a != a_bot):
                raise ValueError("Lattice constants of bottom layers must be equal.")
        self.reference_uc = self.uc_layers_bottom
        for uc in self.uc_layers_top:
            uc.setReferenceUnitCell(self.uc_layers_bottom[0], kwargs.get('rot', np.identity(3)))

        self.basis = np.array([0., 0.])
        self.parameters = {
            'absolute' : [],
            'relative' : []
        }
        
        self.basis_0 = np.array([0., 0.])
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = "unnamed"
            
    def setReferenceUnitCell(self,uc,rotMatrix=np.identity(3)):
        """   
        set reference unit cell. 
        When any F_hkl is called, first hkl will be transformed from the reference 
        lattice into the lattice of this unit cell,
        using equation
        B' * H' = O * B * H, where a ' indicates the entities in the frame of this
        lattice. O is an optional rotation matrix, which describes the rotation of
        this lattice with respect to the reference lattice.
        
        Equally the same applies to the relative coordinates of the respective 
        lattices:
        R' * x' = O * R * x
        
        !!! This has to be checked !!!
        """
        for uc in self.uc_layers_top:
            uc.setReferenceUnitCell(uc, rotMatrix)
        for uc in self.uc_layers_bottom:
            uc.setReferenceUnitCell(uc, rotMatrix)
        self.reference_uc = uc
        #self._dirty = True
    
    def setEnergy(self,E):
        self.E = E
        for uc in self.uc_layers_top:
            uc.setEnergy(E)
        for uc in self.uc_layers_bottom:
            uc.setEnergy(E)
            
    @property
    def loc_absolute(self):
        if self._dirty:
            self.createInterfaceCells()
        return self._loc_absolute
        
    @property
    def height_absolute(self):
        if self._dirty:
            self.createInterfaceCells()
        H = (self.uc_layers_top[-1].coherentDomainMatrix[-1][2,3] +1)*self.uc_layers_top[1-].a[2]
        return H
        
    def createInterfaceCells(self):
        sigma = self.basis[0] * len(self.uc_layers_top)
        skew = self.basis[1] / sigma #(len(self.uc_layers_top))**2
        
        
        if self.type == 'skellam':
            if abs(sigma * skew) > 1:
                raise ValueError("abs(sigma * skew) must be smaller than one.")
            mu1 = 0.5* sigma**2 * (1. + skew*sigma)
            mu2 = sigma**2 - mu1
            
            loc = mu1 - mu2
            loc_int = int(round(loc/len(self.uc_layers_top),0))*len(self.uc_layers_top)
            
            
            uc_number = (int(np.ceil((self.sigma_calc*sigma) / len(self.uc_layers_top))) + 1) * len(self.uc_layers_top)
            uc_number = int(uc_number)
            unitcells = np.arange(-uc_number, uc_number) + loc_int
            assert unitcells.size % len(self.uc_layers_top) == 0
            
            probability_top = skellam.cdf(unitcells, mu1, mu2).reshape((-1, len(self.uc_layers_top)))
            probability_bottom = 1. - probability_top
            

            a3_top = self.uc_layers_top[0].a[2]
            a3_bottom = self.uc_layers_bottom[0].a[2]
            
            
            for i, (uc_t, uc_b) in enumerate(zip(self.uc_layers_top, self.uc_layers_bottom)):
                uc_t.coherentDomainMatrix = []
                uc_t.coherentDomainOccupancy = probability_top.T[i]
                uc_b.coherentDomainMatrix = []
                uc_b.coherentDomainOccupancy = probability_bottom.T[i]
                
            mat_0 = np.vstack((np.identity(3).T,np.array([0,0,0]))).T
            
            ratio_top = (a3_bottom/a3_top)
            ratio_bottom = 1/ratio_top
            h_top = 0.
            h_bottom = 0.
            h = 0.
            
            for p_t, p_b in zip(probability_top, probability_bottom):
                
                for i, (uc_t, uc_b) in enumerate(zip(self.uc_layers_top, self.uc_layers_bottom)):
                    mat_top_i = np.copy(mat_0)
                    top_strain_and_h = 1.*p_t[i] + ratio_top*p_b[i]
                    mat_top_i[2,2] = top_strain_and_h
                    mat_top_i[2,3] = h / a3_top
                    
                    uc_t.coherentDomainMatrix.append(mat_top_i)
                    
                    mat_bottom_i = np.copy(mat_0)
                    bottom_strain_and_h = ratio_bottom*p_t[i] + 1.*p_b[i]
                    mat_bottom_i[2,2] = bottom_strain_and_h
                    mat_bottom_i[2,3] = h / a3_bottom
                    #h_bottom += bottom_strain_and_h
                    uc_b.coherentDomainMatrix.append(mat_bottom_i)
                h += a3_top * top_strain_and_h
                #h_bottom += bottom_strain_and_h
                #h_top += top_strain_and_h
            
            loc_rescaled = loc - loc_int + uc_number
            uc_no_loc = int(np.floor(loc_rescaled)) // len(self.uc_layers_top)
            layer_no_loc = int(np.floor(loc_rescaled)) % len(self.uc_layers_top)
            loc_remainder = (loc_rescaled % len(self.uc_layers_top)) % 1
            loc_mat = self.uc_layers_top[layer_no_loc].coherentDomainMatrix[uc_no_loc]
            self._loc_absolute = loc_mat[2,3]*a3_top + loc_remainder*loc_mat[2,2]*a3_top
            
            
            
            self._dirty = False
        else:
            raise NotImplementedError("%s is not a valid interface model" % self.type)
    
    def F_uc(self,h,k,l):
        if self._dirty:
            self.createInterfaceCells()
        if HAS_NUMBA_ACCEL:
            h,k,l = _ensure_contiguous(h,k,l, testOnly=False, astype=np.float64)
        F = np.zeros_like(l, dtype=np.complex128)
        for uc_t, uc_b in zip(self.uc_layers_top, self.uc_layers_bottom):
            F += uc_t.F_uc(h,k,l)
            F += uc_b.F_uc(h,k,l)
        return F
        
    def zDensity_G(self,z,h,k):
        if self._dirty:
            self.createInterfaceCells()
        rho = np.zeros_like(z, dtype=np.complex128)
        for uc_t, uc_b in zip(self.uc_layers_top, self.uc_layers_bottom):
            rho += uc_t.zDensity_G(z,h,k)
            rho += uc_b.zDensity_G(z,h,k)
        return rho
        
        

        


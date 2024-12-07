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
                     estimateDispersionCompound, LinearFitFunctions)

from .CTRuc import UnitCell, HAS_NUMBA_ACCEL




class EpitaxyInterface(LinearFitFunctions):
    parameterOrder = "Width/cells Skew/cells"
    
    parameterLookup = {'W' : 0, 'S' : 1}
    
    avail_types = ['skellam']
    
    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))
    
    def __init__(self,uc_top, uc_bottom, type='skellam' ,**kwargs):
        """
        sigma_calc : automatic selection of number of unitcells for calculation of interface
                     this will shift the location of the interface.
        fixed_ucs : create interface with fixed number of unitcells. Location of interface
                    will only change if skew != 0
                    total number of uc will be fixed_ucs*2 +1
        

        
        """
        super(EpitaxyInterface,self).__init__()
        self.type = type
        self.sigma_calc = kwargs.get('sigma_calc', 4)
        self.fixed_ucs = kwargs.get('fixed_ucs', False)
        self.set_ucs(uc_top, uc_bottom, **kwargs)
        self.basis = np.array([0., 0.])
        self._basis_created = np.array([np.nan, np.nan])
        self.basis_0 = np.array([0., 0.])
        self.errors = None
        if 'name' in kwargs:
            self.name = kwargs['name']
        else:
            self.name = "unnamed"
            
    def set_ucs(self, uc_top, uc_bottom, **kwargs):
        if not np.all(uc_top.layers == uc_bottom.layers):
            raise ValueError("Top and bottom layer numbers must be equal. " 
            "Is Top: %s, Bottom: %s" % (uc_top.layers, uc_bottom.layers))
        self.uc_top = uc_top
        self.uc_bottom = uc_bottom
        self.reference_uc = self.uc_bottom

        self.uc_top.setReferenceUnitCell(self.uc_bottom, kwargs.get('rot', np.identity(3)))
        
        self.uc_layers_top = self.uc_top.split_in_layers()
        self.uc_layers_bottom = self.uc_bottom.split_in_layers()
        
        self.top_layers = np.array(list(self.uc_layers_top.values()))
        self.bottom_layers = np.array(list(self.uc_layers_bottom.values()))
        
        # sort for height:
        avg_height_top = []
        for uc_l in self.uc_layers_top:
            avg_height_top.append(np.mean(self.uc_layers_top[uc_l].basis[:,3]))
        avg_height_top = np.array(avg_height_top)
        
        avg_height_bottom = []
        for uc_l in self.uc_layers_bottom:
            avg_height_bottom.append(np.mean(self.uc_layers_bottom[uc_l].basis[:,3]))
        avg_height_bottom = np.array(avg_height_bottom)
        
        sort_top = np.argsort(avg_height_top)
        sort_bottom = np.argsort(avg_height_bottom)
        
        self.top_layers = self.top_layers[sort_top]
        self.bottom_layers = self.bottom_layers[sort_bottom]
        
            
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
        self.uc_top.setReferenceUnitCell(uc, rotMatrix)
        self.uc_bottom.setReferenceUnitCell(uc, rotMatrix)
        self.reference_uc = uc
        #self._dirty = True
    
    def setEnergy(self,E):
        self.E = E
        self.uc_top.setEnergy(E)
        self.uc_bottom.setEnergy(E)
        
            
    @property
    def loc_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        return self._loc_absolute
        
    @property
    def height_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        H = (self.top_layers[-1].coherentDomainMatrix[-1][2,3] +1)*self.top_layers[-1].a[2]
        return H
    
    @property
    def pos_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        H = self.top_layers[0].coherentDomainMatrix[0][2,3]*self.top_layers[0].a[2]
        return H
    
    @pos_absolute.setter
    def pos_absolute(self, pos):
        for l in self.top_layers:
            l.pos_absolute = pos
        for l in self.bottom_layers:
            l.pos_absolute = pos
            
    
    @property
    def end_layer_number(self):
        return self.top_layers[-1].basis[0,7]
        
        
    def createInterfaceCells(self):
        n_layers = len(self.uc_top.layers)
        sigma = self.basis[0] * n_layers
        if abs(self.basis[1]) > 1:
            raise ValueError("skew must be between -1 and 1.")
        if abs(self.basis[1]) == 1.:
            skew = (self.basis[1] - np.sign(self.basis[1]) * 1e-6) / sigma
        else:
            skew = self.basis[1] / sigma #(len(self.top_layers))**2
        
        if self.type == 'skellam':
            if abs(sigma * skew) > 1:
                raise ValueError("abs(sigma * skew) must be smaller than one.")
            mu1 = 0.5* sigma**2 * (1. + skew*sigma)
            mu2 = sigma**2 - mu1
            
            loc = mu1 - mu2
            loc_int = int(round(loc/n_layers,0))*n_layers
            
            if self.fixed_ucs:
                uc_number = self.fixed_ucs
            else:
                uc_number = (int(np.ceil((self.sigma_calc*sigma) / n_layers)) + 1) * n_layers
                uc_number = int(uc_number)
            unitcells = np.arange(-uc_number, uc_number) + loc_int
            assert unitcells.size % len(self.top_layers) == 0
            
            probability_top = skellam.cdf(unitcells, mu1, mu2).reshape((-1, n_layers))
            probability_bottom = 1. - probability_top
            

            a3_top = self.top_layers[0].a[2]
            a3_bottom = self.bottom_layers[0].a[2]
            
            
            for i, (uc_t, uc_b) in enumerate(zip(self.top_layers, self.bottom_layers)):
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
                
                for i, (uc_t, uc_b) in enumerate(zip(self.top_layers, self.bottom_layers)):
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
            uc_no_loc = int(np.floor(loc_rescaled)) // n_layers
            layer_no_loc = int(np.floor(loc_rescaled)) % n_layers
            loc_remainder = (loc_rescaled % n_layers) % 1
            loc_mat = self.top_layers[layer_no_loc].coherentDomainMatrix[uc_no_loc]
            self._loc_absolute = loc_mat[2,3]*a3_top + loc_remainder*loc_mat[2,2]*a3_top
            
            
            
            self._basis_created = np.copy(self.basis)
        else:
            raise NotImplementedError("%s is not a valid interface model" % self.type)
    
    def F_uc(self,h,k,l):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        if HAS_NUMBA_ACCEL:
            h,k,l = _ensure_contiguous(h,k,l, testOnly=False, astype=np.float64)
        F = np.zeros_like(l, dtype=np.complex128)
        for uc_t, uc_b in zip(self.top_layers, self.bottom_layers):
            F += uc_t.F_uc(h,k,l)
            F += uc_b.F_uc(h,k,l)
        return F
        
    def zDensity_G(self,z,h,k):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        rho = np.zeros_like(z, dtype=np.complex128)
        for uc_t, uc_b in zip(self.top_layers, self.bottom_layers):
            rho += uc_t.zDensity_G(z,h,k)
            rho += uc_b.zDensity_G(z,h,k)
        return rho

    def addFitParameter(self,indexarray,limits=(-np.inf,np.inf),**kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`   
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addFitParameter(self,indexarray,limits,**kwarg)
        if 'unitcell' not in kwarg:
            raise ValueError("Missing unit cell name. Provide unit cell name as kwarg \'unitcell\'")
        if isinstance(kwarg['unitcell'], list):
            fp = []
            for ucn in kwarg['unitcell']:
                fp.append(self[ucn].addFitParameter(indexarray,limits,**kwarg))
            return fp
        else:
            return self[kwarg['unitcell']].addFitParameter(indexarray,limits,**kwarg)

    def addRelParameter(self,indexarray,factors,limits=(-np.inf,np.inf),**kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`   
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addRelParameter(self,indexarray,factors,limits,**kwarg)
        if 'unitcell' not in kwarg:
            raise ValueError("Missing unit cell name. Provide unit cell name as kwarg \'unitcell\'")
        if isinstance(kwarg['unitcell'], list):
            fp = []
            for ucn in kwarg['unitcell']:
                fp.append(self[ucn].addRelParameter(indexarray,factors,limits,**kwarg))
            return fp
        else:
            return self[kwarg['unitcell']].addRelParameter(indexarray,factors,limits,**kwarg)

        
    def getStartParamAndLimits(self, force_recalculate=False):
        #if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        x0, lower, upper = super().getStartParamAndLimits(force_recalculate) # absolute and relative
        top_x0, top_lower, top_upper = self.uc_top.getStartParamAndLimits(force_recalculate)
        bottom_x0, bottom_lower, bottom_upper = self.uc_bottom.getStartParamAndLimits(force_recalculate)
        return (np.concatenate([x0, top_x0, bottom_x0]),
               np.concatenate([lower, top_lower, bottom_lower]),
               np.concatenate([upper, top_upper, bottom_upper]))
    
    def setFitParameters(self,x):
        abs_rel_no = len(self.parameters['absolute']) + len(self.parameters['relative'])
        fp_top_no = len(self.uc_top.fitparnames)
        fp_bottom_no = len(self.uc_bottom.fitparnames)
        super().setFitParameters(x[:abs_rel_no])
        self.uc_top.setFitParameters(x[abs_rel_no: abs_rel_no+fp_top_no])
        self.uc_bottom.setFitParameters(x[abs_rel_no+fp_top_no: abs_rel_no+fp_top_no+fp_bottom_no])
        
    def setLimits(self,lim):
        abs_rel_no = len(self.parameters['absolute']) + len(self.parameters['relative'])
        fp_top_no = len(self.uc_top.fitparnames)
        fp_bottom_no = len(self.uc_bottom.fitparnames)
        super().setLimits(lim[:abs_rel_no])
        self.uc_top.setLimits(lim[abs_rel_no: abs_rel_no+fp_top_no])
        self.uc_bottom.setLimits(lim[abs_rel_no+fp_top_no: abs_rel_no+fp_top_no+fp_bottom_no])
        
    def setFitErrors(self,errors):
        abs_rel_no = len(self.parameters['absolute']) + len(self.parameters['relative'])
        fp_top_no = len(self.uc_top.fitparnames)
        fp_bottom_no = len(self.uc_bottom.fitparnames)
        super().setFitErrors(errors[:abs_rel_no])
        self.uc_top.setFitErrors(errors[abs_rel_no: abs_rel_no+fp_top_no])
        self.uc_bottom.setFitErrors(errors[abs_rel_no+fp_top_no: abs_rel_no+fp_top_no+fp_bottom_no])
        
    def getFitErrors(self):
        err = super().getFitErrors()
        err_t = self.uc_top.getFitErrors()
        err_b = self.uc_bottom.getFitErrors()
        return np.concatenate([err, err_t, err_b])
        
    @property
    def fitparnames(self):
        return super().fitparnames + self.uc_top.fitparnames + self.uc_bottom.fitparnames

    @property
    def priors(self):
        return super().priors + self.uc_top.priors + self.uc_bottom.priors
        
    def parametersToDict(self):
        d = super().parametersToDict()
        d['unitcells'] = {}
        d['unitcells']['top'] = self.uc_top.parametersToDict()
        d['unitcells']['bottom'] = self.uc_bottom.parametersToDict()
        return d
    
    def clearParameters(self):
        super().clearParameters()
        self.uc_top.clearParameters()
        self.uc_bottom.clearParameters()
    
    def parametersFromDict(self, d, override_values=True):
        self.uc_top.parametersFromDict(d['unitcells']['top'])
        self.uc_bottom.parametersFromDict(d['unitcells']['bottom'])
        super().parametersFromDict(d)
        
    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters 
        """
        self.uc_top.updateFromParameters()
        self.uc_bottom.updateFromParameters()
        super().updateFromParameters()

    def __getitem__(self,uc_name_or_index):
        if isinstance(uc_name_or_index,str):
            if uc_name_or_index.lower() in ['top', 't', 'upper', self.uc_top.name]:
                return self.uc_top
            elif uc_name_or_index.lower() in ['bottom', 'b', 'lower', self.uc_bottom.name]:
                return self.uc_bottom
            else:
                raise KeyError("No unit cell %s in EpitaxyInterface %s" % (uc_name_or_index, self.name))
        else:
            raise ValueError("must be str, not {}".format(type(uc_name_or_index)) )
            
    def parameter_list(self):
        return super().parameter_list() + self.uc_top.parameter_list() + self.uc_bottom.parameter_list()
       
    @classmethod
    def fromStr(cls, string):
        xprfile = False
        with util.StringIO(string) as f:
            # parse header 
            line = next_skip_comment(f).split()
            if line[0].lower() != 'type':
                raise ValueError("You must specify a epitaxy type in line 1."
                " Available are %s" % EpitaxyInterface.avail_types)
            if line[1].lower() not in EpitaxyInterface.avail_types:
                raise ValueError("Expitaxy type %s is not valid."
                " Must be one of %s" % (line[1], EpitaxyInterface.avail_types))
            ep_type = line[1].lower()
            
            statistics = dict()
            line = next_skip_comment(f)
            while('Width' in line or '=' in line): # parameter header or statistics line
                if '=' in line:
                    try:
                        splitted = [n.split(',') for n in line.split('=')]
                        splitted = [item for sublist in splitted for item in sublist]
                        for i in range(0,len(splitted),2):
                            statistics[splitted[i].strip()] = float(splitted[i+1])
                    except Exception:
                        print("Cannot read statistics string: %s" % line)
                line = next_skip_comment(f)
            # epitaxy parameters
            sline = line.split()
            if '+-' in sline:
                params = re.findall(r'\(([^)]+)',line)
                params_array = np.array([np.array(p.split('+-'),dtype=np.float64) for p in params]).T
                basis = params_array[0]
                errors = params_array[1]
            else:
                basis = np.array(sline,dtype=np.float64)
                errors = None
        
        # very explicit searching for the lines containing TopUnitCell and BottomUnitCell:
        sp_str = string.splitlines()
        top_pos = -1
        bottom_pos = -1
        for i, l in enumerate(sp_str):
            if top_pos == -1:
                if 'TopUnitCell' in l:
                    top_pos = i # found it, and save line number
            if bottom_pos == -1:
                if 'BottomUnitCell' in l:
                    bottom_pos = i # found it, and save line number
            if top_pos != -1 and bottom_pos != -1:
                break
        else:
            msg = "Cannot create EpitaxyInterface. "
            if top_pos < 0:
                msg += "No TopUnitCell provided. "
            if bottom_pos < 0:
                msg += "No BottomUnitCell provided."
            raise ValueError(msg)
        
        classname_top,top_name = sp_str[top_pos].split(maxsplit=1)
        classname_bottom,bottom_name = sp_str[bottom_pos].split(maxsplit=1)
        
        assert classname_top == 'TopUnitCell'
        assert classname_bottom == 'BottomUnitCell'
        
        if top_pos < bottom_pos:
            top_uc_str = '\n'.join(sp_str[top_pos+1:bottom_pos])
            bottom_uc_str = '\n'.join(sp_str[bottom_pos+1:])
        else:
            top_uc_str = '\n'.join(sp_str[top_pos+1:])
            bottom_uc_str = '\n'.join(sp_str[bottom_pos+1:top_pos])
        
        uc_top = UnitCell.fromStr(top_uc_str)
        uc_bottom = UnitCell.fromStr(bottom_uc_str)
        
        uc_top.name = top_name
        uc_bottom.name = bottom_name
        
        epit = cls(uc_top, uc_bottom, ep_type)
        epit.statistics = statistics
        epit.basis = basis
        epit.basis_0 = np.copy(basis)
        epit.errors = errors
        return epit
        
        
    def toStr(self):
        s = "type %s" % self.type
        s += "\n" + EpitaxyInterface.parameterOrder + "\n" + self.epitToStr()
        s += "\n\n"
        s += "TopUnitCell %s\n" % self.uc_top.name
        s += self.uc_top.toStr() + "\n\n"
        s += "BottomUnitCell %s\n" % self.uc_bottom.name
        s += self.uc_bottom.toStr() + "\n"
        return s
        
    def __repr__(self):
        return self.toStr()
        
    def epitToStr(self,showErrors=True):
        param = self.basis
        if (self.errors is not None) and showErrors:
            err = self.errors
            l = []
            for p, err in zip(param,err):
                l.append("(%.5f +- %.5f)" % (p, err))
            return "   ".join(l)
        else:
            l = []
            for p in param:
                l.append("%.5f " % p)
            return "   ".join(l)
        
            
    
        
            
            
        

        


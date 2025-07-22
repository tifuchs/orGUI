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

from .CTRutil import (special_elementcolors, ParameterType, Parameter,
                     _ensure_contiguous, next_skip_comment, DWtoDisorder,
                     readWaasmaier, readDispersion, atomic_number,
                     estimateDispersionCompound)

from .CTRuc import WaterModel, UnitCell
from .CTRfilm import EpitaxyInterface, Film

class SXRDCrystal(object):

    def __init__(self,uc_bulk,*uc_surface,**keyargs):
        self.uc_bulk = uc_bulk
        self.uc_surface_list = list(uc_surface)
        self.enable_uc_stacking = keyargs.get('enable_stacking', False)
        if not self.enable_uc_stacking:
            for uc in self.uc_surface_list:
                if isinstance(uc, (Film, EpitaxyInterface) ):
                    self.enable_uc_stacking = True
                    break
        
        self.uc_stacking = keyargs.get('stacking', np.arange(len(self.uc_surface_list))[::-1])
        order = np.argsort(self.uc_stacking)
        self.uc_surface_list_ordered = np.array(self.uc_surface_list)[order]
        self.uc_stacking_ordered = self.uc_stacking[order]
        self.domains = [[(np.identity(3,dtype=np.float64),1)] for i in self.uc_surface_list]
        self.atten = keyargs.get('atten', 0.01)
        self.specularRes = keyargs.get('spec_res', 0.0)
        self.name = keyargs.get('name', 'xtal')
        
        self.weights = np.ones(len(self.uc_surface_list),dtype=np.float64)/len(self.uc_surface_list)
        self.weights_0 = np.copy(self.weights)
        self._weights_parvalues = None
        self.werrors = None
        self._werrors_parvalues = None
        self.parameters = {
            'coupled' : [],
            'weight' : [],
            'domain' : []
        }
        self._parIdNo = 0
        self.fit_metadata_cache = None
        
    def parametersToDict(self):
        d = dict()
        d['weights_0'] = self.weights_0
        d['_parIdNo'] = self._parIdNo
        
        for p in self.parameters:
            d[p] = dict()
            for i, param in enumerate(self.parameters[p]):
                d[p]["%s_%s" % (i,param.name)] = param.asdict()
        
        d['unitcells'] = dict()
        for uc in self.uc_surface_list:
            d['unitcells'][uc.name] = uc.parametersToDict()
        d['unitcells']['bulk'] = self.uc_bulk.parametersToDict()
        return d
        
    def clearParameters(self):
        for p in self.parameters:
            self.parameters[p] = []
        self._parIdNo = 0
        self.fit_metadata_cache = None
        self._weights_parvalues = None
        self._werrors_parvalues = None
        for uc in self.uc_surface_list:
            uc.clearParameters()
        self.uc_bulk.clearParameters()
    
    def parametersFromDict(self, d, override_values=True):
        for p in self.parameters:
            self.parameters[p] = []
            
        self.weights_0 = d['weights_0'].astype(np.float64)
        self._parIdNo = int(d['_parIdNo']) 

        for p in self.parameters:
            if p in d:
                for dkey in sorted(d[p].keys()):
                    self.parameters[p].append(Parameter(**d[p][dkey]))
        
        for uc in d['unitcells']:
            self[uc].parametersFromDict(d['unitcells'][uc], override_values)

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn("Can not apply weight fit parameter values to this SXRDCrystal: %s" % e) 
        else:
            warnings.warn("Weight values not updated from Parameter values. This can cause a mismatch between weight and fitparameter values!")
        
    def updateFromParameters(self):
        self.weights = np.copy(self.weights_0)
        self.werrors = np.full_like(self.weights,np.nan)
        no_errors = False
        for par in self.parameters['weight']:
            if par.value is not None:
                self.weights[par.indices] += par.factors*par.value
                self.werrors[par.indices] = np.nan_to_num(self.werrors[par.indices],nan=0.0)
            else:
                raise ValueError("Can not set weight values from parameters. Value of Parameter %s is None." % par.name)
            if par.error is not None:
                self.werrors[par.indices] += par.factors*par.error
            else:   
                no_errors = True
        self._weights_parvalues = np.copy(self.weights)
        if not no_errors:
            self._werrors_parvalues = np.copy(self.werrors)
        
        
    def setGlobalReferenceUnitCell(self,uc_reference,rotMatrix=np.identity(3)):
        """   
        sets a global reference unit cell. 
        When any F_hkl is called, first hkl will be transformed from the reference 
        lattice into the lattice of this unit cell,
        using equation
        B' * H' = O * B * H, where a ' indicates the entities in the frame of this
        lattice. O is an optional rotation matrix, which describes the rotation of
        this lattice with respect to the reference lattice.
        
        Equally the same applies to the relative coordinates of the respective 
        lattices:
        R' * x' = O * R * x
        
        """
        self.uc_bulk.setReferenceUnitCell(uc_reference,rotMatrix)
        [uc.setReferenceUnitCell(uc_reference,rotMatrix) for uc in self.uc_surface_list]
    
    def F_surf(self,harray,karray,Larray):
        F = np.empty_like(Larray,dtype=np.complex128)
        for uc, weight in zip(self.uc_surface_list,self.weights):
            F += weight*uc.F_uc(harray,karray,Larray)
        return F
    
    def apply_stacking(self):
        self.enable_uc_stacking = True
        hnew = 0.
        i = 0
        layer_number_new = -1
        locnew = 0
        for l in np.unique(self.uc_stacking_ordered):
            h = hnew
            layer_number = layer_number_new
            loc = locnew
            while(self.uc_stacking_ordered[i] == l):
                uc = self.uc_surface_list_ordered[i]
                if isinstance(uc, Film):
                    uc.set_below(loc, h)
                uc.start_layer_number = layer_number
                uc.pos_absolute = h
                hnew = uc.height_absolute
                locnew = uc.loc_absolute
                layer_number_new = uc.end_layer_number
                i += 1
                if i == len(self.uc_stacking_ordered):
                    return
                
    
    def F(self,harray,karray,Larray):
        F = self.uc_bulk.F_bulk(harray,karray,Larray,self.atten)
        hkl = np.vstack((harray,karray,Larray))
        if self.enable_uc_stacking:
            self.apply_stacking()

        for uc, weight, domains in zip(self.uc_surface_list,self.weights,self.domains):
            for matrix, occup in domains:
                hkl_n = np.dot(matrix,hkl)
                F += occup*weight*uc.F_uc(hkl_n[0],hkl_n[1],hkl_n[2])
        return F
    

    def setDomain(self,uc_no,domains):
        """
        should be a list of tuples with:
        (rotmatrix, occupancy)
        """
        self.domains[uc_no] = domains
    
    def addRelParameter(self,index0_or_name,indexneg_or_name=None,initial=None,limits=(-np.inf,np.inf),prior=None,**keyargs):
        """
        sets a new fit parameter for the weight w of a unit cell (Legacy API!, consider using addWeightParameter)
        
        F_xtal = w * F_uc1
        
        if indexneg_or_name is not None the structure factor will be calculated
        to
        
        F_xtal = w * F_uc1 +  (1 - w) * F_uc2
        
        Parameters
        ----------
        index0_or_name : str or int
            index or name of the unit cell 1 (with F_uc1)
        indexneg_or_name : default: None, str or int
            index or name of the unit cell 2 (with F_uc2)
        initial: float
            initial value of w
        limits: list
            fit limits of w 
        """
        
        index0 = self.getUcIndex(index0_or_name)
        
        if initial is None:
            initial = self.weights[index0]
        
        self.weights_0[index0] = 0.
        self.weights[index0] = initial
        if indexneg_or_name is not None:
            indexneg = self.getUcIndex(indexneg_or_name)
            self.weights_0[indexneg] = 1.
            self.weights[indexneg] = 1. - initial
            names = np.array([index0, indexneg],dtype=np.intp)
            factors = np.array([1., -1.])
        else:
            names = np.array([index0])
            factors = np.array([1.])
            
        keyargs['limits'] = limits
        keyargs['prior'] = prior
        
        return self.addWeightParameter(names,factors,**keyargs)
        
    def addWeightParameter(self,indices_or_names,factors,**keyargs):
        """
        sets a new fit parameter for the weight w of a unit cell
        
        Parameters
        ----------
        index0_or_name : str or int
            index or name of the unit cell 1 (with F_uc1)
        indexneg_or_name : default: None, str or int
            index or name of the unit cell 2 (with F_uc2)
        initial: float
            initial value of w
        limits: list
            fit limits of w 
        """
        
        prior = keyargs.get('prior', None)
        name = keyargs.get('name', "%s_weight_%s" % (self.name, self._parIdNo))
        limits = keyargs.get('limits', (-np.inf,np.inf))
        initial = keyargs.get('initial', None)
        basis_init = keyargs.get('basis_init', None)
        
        indexarray = np.array([self.getUcIndex(ucn) for ucn in indices_or_names],dtype=np.intp)
        
        factors = np.array(factors,dtype=np.float64)
        
        if basis_init is not None:
            self.weights_0[indexarray] = basis_init
            
        if initial is not None:
            self.weights[indexarray] = self.weights_0[indexarray] + factors*initial
        par = Parameter(name,indexarray, ParameterType.RELATIVE, limits, prior, keyargs, factors)
        self.parameters['weight'].append(par)

        self._parIdNo += 1
        
        self.fit_metadata_cache = None
        
        return par
    
    
    def addFitParameter(self,parameter,limits=(-np.inf,np.inf),**keyargs):
        """
        API for fitparameters:
        
        par = {
            'ucname1':
            {
                'atoms' : (1, 2, 3), # atom index
                'par' : ('z', 'oDW', 'iDW'), # structure parameter
                'factors' : (1, -1, 1) # optional, if present, this will be a "relative" fit parameter
                'settings' : optional, dict with keyargs passed to the unit cell fitparameter constructor 
            }
        }
        
        
        """
        prior = keyargs.get('prior', None)
        name = keyargs.get('name', "%s_par_%s" % (self.name, self._parIdNo))
        keyargs['name'] = name
        keyargs['prior'] = prior
        
        
        ucnames = np.array(list(parameter.keys()))
        ucindices = np.array([self.getUcIndex(ucn) for ucn in ucnames])
        sidx = np.argsort(ucindices)
        
        ucindices = ucindices[sidx]
        ucnames = ucnames[sidx]
        
        if not isinstance(parameter[ucnames[0]], dict): #  legacy API fix!
            for uc_identifier in ucnames:
                atoms, par = parameter[uc_identifier]
                parameter[uc_identifier] = {
                    'atoms' : atoms,
                    'par' : par
                }
        
        #children = []
        
        partype = ParameterType(0)
        for uc_identifier in ucnames:
            par = parameter[uc_identifier]
            settings = par.get('settings', dict())
            settings.update(keyargs)
            if 'factors' in par: # relative parameters
                partype |= ParameterType.RELATIVE
                fp = self[uc_identifier].addRelParameter((par['atoms'], par['par']), par['factors'] ,limits,**settings)
                #children.append(fp)
            else:
                partype |= ParameterType.ABSOLUTE
                fp = self[uc_identifier].addFitParameter((par['atoms'], par['par']),limits,**settings)
                #children.append(fp)

        par = Parameter(name, ucindices, partype, limits, prior, keyargs)
        #par.children = children

        self.parameters['coupled'].append(par)
        self._parIdNo += 1
        
        self.fit_metadata_cache = None

        return par
        

    def getSurfaceBasis(self):
        return np.concatenate([uc.basis for uc in self if isinstance(uc,UnitCell)])
        
    def getSurfaceDWConstraintEnable(self):
        return np.concatenate([uc.dw_increase_constraint for uc in self if isinstance(uc,UnitCell)])

    def fitparameterList(self, verbosity=2):
        s = "## Parameters\n"
        for ucn in self.getUcNames():
            uc = self[ucn]
            s += "# %s %s\n" % (uc.__class__.__name__ ,uc.name)
            s += uc.fitparameterList(verbosity) + "\n"
            
        s += "## UnitCell Weights\n"
        s += "{:10} {:10} {:10} {:10} {:10}\n".format("uc1","occup","uc2","occup","Limits")
        for [index0,indexneg,value], lim in zip(self.weightparameters, self.weightparameterLimits):
            uc1name = self.uc_surface_list[index0].name
            uc1value = value
            if indexneg is not None:
                uc2name = self.uc_surface_list[indexneg].name
                uc2value = 1. - value
            else:
                uc2name = ""
                uc2value = ""
            s += "{:10} {:10} {:10} {:10} {:10}\n".format(str(uc1name),str(uc1value),str(uc2name),str(uc2value),str(lim))
            
        return s
        
    def __uc_fitparnames(self):
        names = []
        names += self.uc_bulk.fitparnames 
        for uc in self.uc_surface_list:
            names += uc.fitparnames
        return names
                    
    @property
    def fitparnames(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        names = [par.name for par in self.parameters['coupled'] + self.parameters['weight'] + self.parameters['domain']]
        names += self.__uc_fitparnames()
        names = list(np.array(names)[self.fit_metadata_cache['par_mask']])
        return names
    
    @property
    def priors(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        priorlist = []
        for par in self.parameters['coupled'] + self.parameters['weight'] + self.parameters['domain']:
            if par.prior is None:
                if tuple(par.limits) == (-np.inf,np.inf):
                    raise Exception(f"Infinite range given as prior for parameter {par.name}. You must provide an explicit prior.")
                else:
                    priorlist.append(tuple(par.limits)) #has to be converted to correct distribution in CTROptimizer
            else:
                priorlist.append(par.prior) # real prior distribution
        
        priorlist += self.uc_bulk.priors
        for uc in self.uc_surface_list:
            priorlist += uc.priors
        priorlist = list(np.array(priorlist)[self.fit_metadata_cache['par_mask']])
        return priorlist
    
    def parameter_list(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        l = self.parameters['coupled'] + self.parameters['weight'] + self.parameters['domain']
        l += self.uc_bulk.parameter_list() 
        for uc in self.uc_surface_list:
            l += uc.parameter_list()
        l = list(np.array(l)[self.fit_metadata_cache['par_mask']])
        return l
        
    
    def getInitialParameters(self, force_recalculate=False):
        p, _, _ = self.getStartParamAndLimits(force_recalculate)
        return p    
        
    
    def getStartParamAndLimits(self, force_recalculate=False):
        x0 = []; lower = []; higher = []
        try: 
            mismatch = not np.allclose(self._weights_parvalues, self.weights, equal_nan=True)
        except Exception:
            mismatch = True
        recalculate = force_recalculate or mismatch
        
        number_xtalpar = sum(len(p) for p in self.parameters.values()) 
        
        par_names = []    
        for par in self.parameters['weight']:
            if recalculate or par.value is None:  
                x0.append(np.mean( (self.weights[par.indices] - self.weights_0[par.indices]) / par.factors))
            else:
                x0.append(par.value)
            lower.append(par.limits[0]); higher.append(par.limits[1])
            par_names.append(par.name)
        
        for par in self.parameters['domain']:
            raise NotImplementedError("Domain fitting not yet implemented")
            lower.append(par.limits[0]); higher.append(par.limits[1])
            par_names.append(par.name)
        
        uc_numberpars = []
        p, low, high = self.uc_bulk.getStartParamAndLimits()
        x0.extend(p); lower.extend(low); higher.extend(high)
        
        uc_numberpars.append(len(p))
        
        for uc in self.uc_surface_list:
            p, low, high = uc.getStartParamAndLimits()
            uc_numberpars.append(len(p))
            x0.extend(p); lower.extend(low); higher.extend(high)
        
        par_names += self.__uc_fitparnames()
        
        number_coupled = len(self.parameters['coupled'])
        
        par_names = np.array(par_names)
        uc_par_idx = np.arange(par_names.size) + number_coupled
        uc_par_mask = np.ones(par_names.size,dtype=np.bool_)
        
        x0_coupled = []; lower_coupled = []; higher_coupled = []
        idx = 0
        for par in self.parameters['coupled']:
            paridx = par_names == par.name
            x0_coupled.append(np.mean(np.array(x0)[paridx]))
            uc_par_idx[paridx] = idx
            uc_par_mask[paridx] = False
            
            lower_coupled.append(par.limits[0]); higher_coupled.append(par.limits[1])
            idx += 1
            
        uc_par_idx_trunc = uc_par_idx[uc_par_mask]
        uc_par_idx_trunc_idx = np.arange(uc_par_idx_trunc.size) + number_coupled
        uc_par_idx[uc_par_mask] = uc_par_idx_trunc_idx
        
        par_idx_sortarray = np.concatenate((np.arange(number_coupled),uc_par_idx)).astype(np.intp)
        
        self.fit_metadata_cache = {
            'par_idx_sortarray' : par_idx_sortarray,
            'par_mask' : np.concatenate((np.ones(number_coupled,dtype=np.bool_),uc_par_mask)),
            'number_xtalpar' : number_xtalpar,
            'number_coupled' : number_coupled,
            'uc_numberpars' : uc_numberpars
        }
        
        x0 = np.array(x0)[uc_par_mask]
        lower = np.array(lower)[uc_par_mask]
        higher = np.array(higher)[uc_par_mask]

        return np.concatenate((x0_coupled,x0)), np.concatenate((lower_coupled,lower)), np.concatenate((higher_coupled,higher))
        
        
        
    """
    def _readFitparValues(self):

        updates fitp_values from fitparameters
        
        reads parameter values in self.fitparameters and writes 
        the mean value to the corresponding entry in fitp_values

        Returns
        -------
        None.


        values = np.zeros(len(self.fitp_values))
        count = np.zeros(len(self.fitp_values))
        for uc_pars in self.fitparameters:
            for p in uc_pars:
                values[int(p[0])] += p[2]
                count[int(p[0])] += 1.
            
        self.fitp_values = values/count        
    """
    def setParameters(self,x):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        x = np.asarray(x)
        x_r = x[self.fit_metadata_cache['par_idx_sortarray']] # reorder fitpars, this handles coupled parameters
        x_ucs = x_r[self.fit_metadata_cache['number_xtalpar']:]
        uc_numberpars = self.fit_metadata_cache['uc_numberpars']
        
        lower = 0
        for uc,noparam in zip([self.uc_bulk] + self.uc_surface_list,uc_numberpars):
            uc.setFitParameters(x_ucs[lower:noparam + lower])
            lower = noparam + lower
        
        for i, par in enumerate(self.parameters['coupled']):
            par.value = x[i]
        
        self.weights = np.copy(self.weights_0)
        idx = self.fit_metadata_cache['number_coupled']
        for par in self.parameters['weight']:
            self.weights[par.indices] += par.factors*x_r[idx]
            par.value = x_r[idx]
            idx += 1
            
        
        self._weights_parvalues = np.copy(self.weights)
        
        
    def setLimits(self,lim):
        """fit bounds.
        
        lim shape: (n, 2)
        
        """
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        lim = np.asarray(lim)
        x_r = lim[self.fit_metadata_cache['par_idx_sortarray']] # reorder fitpars, this handles coupled parameters
        x_ucs = x_r[self.fit_metadata_cache['number_xtalpar']:]
        uc_numberpars = self.fit_metadata_cache['uc_numberpars']
        
        lower = 0
        for uc,noparam in zip([self.uc_bulk] + self.uc_surface_list,uc_numberpars):
            uc.setLimits(x_ucs[lower:noparam + lower])
            lower = noparam + lower
        
        for i, par in enumerate(self.parameters['coupled']):
            if lim[i][0] > lim[i][1]:
                raise ValueError("Upper fit bound of coupled par must be larger than the lower bound.")
            par.limits = (lim[i][0],lim[i][1])
            
        idx = self.fit_metadata_cache['number_coupled']
        for par in self.parameters['weight']:
            if x_r[idx][0] > x_r[idx][1]:
                raise ValueError("Upper fit bound of weight pars must be larger than the lower bound.")
            par.limits = (x_r[idx][0],x_r[idx][1])
            idx += 1
        
    def setFitErrors(self,errors):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        errors = np.asarray(errors)
        errors_r = errors[self.fit_metadata_cache['par_idx_sortarray']] # reorder fitpars
        errors_ucs = errors_r[self.fit_metadata_cache['number_xtalpar']:]
        uc_numberpars = self.fit_metadata_cache['uc_numberpars']
        
        lower = 0
        for uc,noparam in zip([self.uc_bulk] + self.uc_surface_list,uc_numberpars):
            uc.setFitErrors(errors_ucs[lower:noparam + lower])
            lower = noparam + lower
        
        for i, par in enumerate(self.parameters['coupled']):
            par.error = errors[i]
        
        self.werrors = np.full_like(self.weights,np.nan)
        idx = 0            
        for par in self.parameters['weight']:
            self.werrors[par.indices] = np.nan_to_num(self.werrors[par.indices],nan=0.0)
            self.werrors[par.indices] += par.factors*errors_r[idx]
            par.error = errors_r[idx]
            idx += 1
        self._werrors_parvalues = np.copy(self.werrors)
    
    def getFitErrors(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits() # will generate metadata, if not yet done so.
        if self.werrors is None:
            raise ValueError("No errors for SXRDCrystal UnitCell weights have been set.")
        try: 
            mismatch = not np.allclose(self._werrors_parvalues, self.werrors, equal_nan=True)
        except Exception:
            mismatch = True
        err0 = []

        number_xtalpar = sum(len(p) for p in self.parameters.values()) 
        
        par_names = []    
        for par in self.parameters['weight']:
            if mismatch:
                err0.append(np.mean(self.werrors[par.indices] / np.abs(par.factors)))
            else:
                err0.append(par.error)
            par_names.append(par.name)
        
        for par in self.parameters['domain']:
            raise NotImplementedError("Domain fitting not yet implemented")
            #lower.append(par.limits[0]); higher.append(par.limits[1])
            #par_names.append(par.name)
        
        #uc_numberpars = []
        p = self.uc_bulk.getFitErrors()
        err0.extend(p)
        
        #uc_numberpars.append(len(p))
        
        for uc in self.uc_surface_list:
            p = uc.getFitErrors()
            #uc_numberpars.append(len(p))
            err0.extend(p)
        
        par_names += self.__uc_fitparnames()
        
        par_names = np.array(par_names)

        err0_coupled = []
        for par in self.parameters['coupled']:
            paridx = par_names == par.name
            err0_coupled.append(np.mean(np.array(err0)[paridx]))

        err = np.concatenate((err0_coupled,err0))[self.fit_metadata_cache['par_mask']]

        return err

    
    def setEnergy(self,E):
        for uc in self.uc_surface_list:
            uc.setEnergy(E)
        self.uc_bulk.setEnergy(E)
        
        
    def zDensity_G(self,z,h,k):
        if self.enable_uc_stacking:
            self.apply_stacking()
        rho = self.uc_bulk.zDensity_G_asbulk(z,h,k)
        for uc,w in zip(self.uc_surface_list,self.weights):
            rho += uc.zDensity_G(z,h,k)*w
        return rho
        
        
    def toRODStr(self):
        s = "E = %.5f keV\n" % (self.uc_bulk._E*1e-3)
        for i,w,uc in zip(self.uc_stacking, self.weights,self.uc_surface_list):
            s += "# %s %s\n" % (uc.__class__.__name__ ,uc.name)
            s += "%04i %.5f\n" % (i,w)
            s += uc.toRODStr() + "\n"
            
        s+= "# UnitCell bulk\n" + self.uc_bulk.toRODStr()
        return s
    
    def toStr(self,showErrors=True):
        s = "E = %.5f keV\n" % (self.uc_bulk._E*1e-3)
        for i ,(no,w,uc) in enumerate(zip(self.uc_stacking,self.weights,self.uc_surface_list)):
            s += "# %s %s\n" % (uc.__class__.__name__ ,uc.name)
            if showErrors and self.werrors is not None:
                s += "%04i occupancy = %.5f +- %.5f\n" % (no,w,self.werrors[i])
            else:
                s += "%04i occupancy = %.5f\n" % (i,w)
            s += uc.toStr() + "\n"
            
        s+= "# UnitCell bulk\n" + self.uc_bulk.toStr()
        return s
    
    @staticmethod
    def fromStr(string,atten=0.01):
        weights = []
        werrors = []
        uc_suface = []
        uc_stacking = []
        
        strio = util.StringIO(string)
        Estr = next_skip_comment(strio).split('=')[1].split('keV')[0]
        E = float(Estr)*1e3
        string = strio.read()
        errors = False
        splstring = string.split('#')
        unitCells = splstring[1:]
        for i in range(len(splstring)-1,0, -1): # ugly fix for commenting of line with #
            lineprefix = splstring[i-1][splstring[i-1].rfind('\n'):] # get contents of line before '#'
            if '//' in lineprefix:
                del unitCells[i-1] # index of unitCells is already shifted by 1 in "unitCells = splstring[1:]"
        
        for ucstr in unitCells[:-1]:
            strio = util.StringIO(ucstr)
            line1 = next_skip_comment(strio)
            classname,name = line1.split(maxsplit=1)
            name = name.strip()
            line2 = next_skip_comment(strio).split()
            uc_stacking.append(int(line2[0]))
            if '+-' in line2:
                errors = True
                weights.append(float(line2[line2.index('+-')-1]))
                werrors.append(float(line2[line2.index('+-')+1]))
            elif '=' in line2:
                weights.append(float(line2[line2.index('=')+1]))
                werrors.append(float('nan'))
            else:
                weights.append(float(line2[1]))
                werrors.append(float('nan'))
            if classname == 'UnitCell':
                uc = UnitCell.fromStr(strio.read())
            elif classname == 'WaterModel':
                uc = WaterModel.fromStr(strio.read())
            elif classname == 'EpitaxyInterface':
                uc = EpitaxyInterface.fromStr(strio.read())
            elif classname == 'Film':
                uc = Film.fromStr(strio.read())
            else:
                raise NotImplementedError("class name not understood: %s" % classname)
            uc.name = name
            uc.setEnergy(E)
            uc_suface.append(uc)
        bulkstrio = util.StringIO(unitCells[-1])
        next_skip_comment(bulkstrio)
        uc_bulk = UnitCell.fromStr(bulkstrio.read())
        uc_bulk.setEnergy(E)
        uc_stacking = np.array(uc_stacking)
        xtal = SXRDCrystal(uc_bulk,*uc_suface,atten=atten,stacking=uc_stacking)
        xtal.weights = np.array(weights)
        xtal.weights_0 = np.copy(weights)
        xtal.werrors = np.array(werrors) if errors else None
        xtal.setGlobalReferenceUnitCell(xtal['bulk'])
        return xtal
        
    def __getitem__(self,uc_name_or_index):
        if isinstance(uc_name_or_index,int):
            if uc_name_or_index == -1:
                return self.uc_bulk
            return self.uc_surface_list[uc_name_or_index]
        elif isinstance(uc_name_or_index,str):
            if uc_name_or_index == 'bulk':
                return self.uc_bulk
            namelist = [uc.name for uc in self.uc_surface_list]
            idx = namelist.index(uc_name_or_index)
            return self.uc_surface_list[idx]
        else:
            raise ValueError("must be str or int, not {}".format(type(uc_name_or_index)) )
            
    def getUcIndex(self,uc_name_or_index):
        if isinstance(uc_name_or_index,(int,np.integer)):
            return int(uc_name_or_index)
        elif isinstance(uc_name_or_index,str):
            if uc_name_or_index == 'bulk':
                return -1
            namelist = [uc.name for uc in self.uc_surface_list]
            idx = namelist.index(uc_name_or_index)
            return idx
        else:
            raise ValueError("must be str or int, not {}".format(type(uc_name_or_index)) )
    
    def plot3d(self,ucx=2,ucy=2,ucz=-3,dwon=False,occuon=False,figure=None,translate=np.array([0.,0.,0.]), **keyargs):
        try:
            from mayavi import mlab
        except ImportError:
            warnings.warn("can not import mayavi: 3D plotting not supported")
            return
        
        if figure is None:
            figure = mlab.figure()
            
        for mat in self.uc_bulk.coherentDomainMatrix:
            self.uc_bulk.plot3d(ucx,ucy,ucz,dwon,occuon,figure,mat,**keyargs)
            
        for uc in self.uc_surface_list:
            for mat in uc.coherentDomainMatrix:
                uc.plot3d(ucx,ucy,1,dwon,occuon,figure,mat,**keyargs)
        
        
            
    def getUcNames(self):
        return [uc.name for uc in self.uc_surface_list] + ['bulk']
        
    def _ipython_key_completions_(self):
        return self.getUcNames()

    def toFile(self,filename):
        with open(filename,'w') as f:
            f.write(self.toStr())
        
    def toRODFile(self,filename):
        with open(filename,'w') as f:
            f.write(self.toRODStr())
    
    @staticmethod
    def fromFile(filename):
        f, ext = os.path.splitext(filename)
        if ext == '':
            for fext in [".xpr", '.xtal']:
                if os.path.isfile(filename + fext):
                    with open(filename + fext,'r') as f:
                        fstr = f.read()
                    xtal = SXRDCrystal.fromStr(fstr)
                    xtal.name = os.path.basename(filename + fext)
                    break
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
            if os.path.isfile(filename + ".h5"):
                try:
                    from silx.io import dictdump
                    fitdict = dictdump.h5todict(filename + ".h5")
                    xtal.parametersFromDict(fitdict)
                except Exception as e:
                    warnings.warn("Fitsettings for SXRDCrystal %s seem to exist, but cannot be read: %s" % (os.path.basename(filename), e))
        else:
            with open(filename,'r') as f:
                fstr = f.read()
            xtal = SXRDCrystal.fromStr(fstr)
            xtal.name = os.path.basename(filename)
        return xtal
    
    def __repr__(self):
        return self.toStr()
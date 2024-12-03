from .HKLVlieg import Lattice
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

from .element_data import cov_radii_array, rgb_array

from .CTRutil import (special_elementcolors, ParameterType, Parameter,
                     _ensure_contiguous, next_skip_comment, DWtoDisorder,
                     readWaasmaier, readDispersion, atomic_number,
                     estimateDispersionCompound)

try:
    from . import _CTRcalc_accel
    HAS_NUMBA_ACCEL = True
except:
    HAS_NUMBA_ACCEL = False

class WaterModel(Lattice):
    
    #path = os.path.split(__file__)[0]

    #f_water = np.loadtxt("%s/water_scattering.dat" % path).T
    #f_water[0] *= 4*np.pi/10.
    #f_water[1] *= np.sqrt(18.015)
    #f_water[2] *= np.sqrt(18.015)
    
    # 2 for measured bulk liquid water
    # 1 for calculated "free atom model" 
    #f0_water_fun = interp1d(f_water[0],f_water[2])
    supportedTypes = ['step', 'layered','1layer','1layer_O','layered_O','1layer_s_O']
    
    parameterOrder = "zpos/frac layer_spacing/frac sigma_0/frac  sigma_bar/frac rho_0_r"
    
    parameterLookup = {'z' : 0, 'spacing' : 1,
                       'sigma_0' : 2, 'sigma_bar' : 3, 'rho_0' : 4}
    
    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))
    
    def __init__(self,a,alpha,watertype='step',**keyargs):
        super(WaterModel,self).__init__(a,alpha)
        self.basis = np.array([0.,2.8,0.5,1.,0.])
        self.pw = 0.0334
        self.fitparameters = []
        self.relfitparam = []
        self.fitparlimits = []
        self.fitparameters_name = []
        self.relfitparam_name = []
        self.relfitlimits = []
        self.relfitparam_prior = []
        self.fitparameter_prior = []
        self.parameters = {
            'absolute' : [],
            'relative' : []
        }
        self.basis_0 = None
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None
        if 'name' in keyargs:
            self.name = keyargs['name']
        else:
            self.name = "unnamed water model"
        if watertype not in WaterModel.supportedTypes:
            raise Exception("Water model %s is not supported, should be one of %s" % (watertype,WaterModel.supportedTypes))
        self.type = watertype
            
        self.refHKLTransform = np.identity(3)
        self.refRealTransform = np.identity(3)
        
    def parametersToDict(self):
        d = dict()
        d['basis_0'] = self.basis_0
        d['dw_increase_constraint'] = self.dw_increase_constraint
        for p in self.parameters:
            d[p] = dict()
            for i, param in enumerate(self.parameters[p]):
                d[p]["%s_%s" % (i,param.name)] = param.asdict()
        return d
    
    def clearParameters(self):
        for p in self.parameters:
            self.parameters[p] = []
        self._basis_parvalues = None
        self._errors_parvalues = None
    
    def parametersFromDict(self, d, override_values=True):
        for p in self.parameters:
            self.parameters[p] = []
            
        self.basis_0 = d['basis_0'].astype(np.float64)
        self.dw_increase_constraint = d['dw_increase_constraint'].astype(np.bool_)
        for p in self.parameters:
            for dkey in sorted(d[p].keys()):
                self.parameters[p].append(Parameter(**d[p][dkey]))

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn("Can not apply fit parameter values to this unitcell: %s" % e) 
        else:
            warnings.warn("Basis values not updated from Parameter values. This can cause a mismatch between basis and fitparameter values!")

    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters 
        """
        no_errors = False
        self.basis = np.copy(self.basis_0)
        self.errors = np.full_like(self.basis,np.nan)
        for par in self.parameters['absolute']:
            if par.value is not None:
                self.basis[par.indices] = par.value
            else:
                raise ValueError("Can not set basis values from parameters. Value of Parameter %s is None." % par.name)
            if par.error is not None:
                self.errors[par.indices] = par.error
            else:
                no_errors = True
        for val, par in zip(x_r,self.parameters['relative']):
            if par.value is not None:
                self.basis[par.indices] += par.factors*par.value
            else:
                raise ValueError("Can not set basis values from parameters. Value of Parameter %s is None." % par.name)
            if par.error is not None:
                self.errors[par.indices] = np.nan_to_num(self.errors[par.indices],nan=0.0)
                self.errors[par.indices] += par.factors*par.error
            else:
                no_errors = True

        self._basis_parvalues = np.copy(self.basis)
        if not no_errors:
            self._errors_parvalues = np.copy(self.errors)
    
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
        self.refRealTransform = self.R_mat_inv @ rotMatrix @ uc.R_mat 
        self.refHKLTransform = self.B_mat_inv @ rotMatrix @ uc.B_mat

    def F_uc(self,h,k,l):
        mask = np.logical_and(np.isclose(h,0), np.isclose(k,0))
        if not np.any(mask):
            return np.zeros_like(l,dtype=np.complex128)
        a,alpha,b,beta = self.getLatticeParameters()
        #Layered water model '''simga_0 is rms of water profile simga_bar is successive broadening where: simga_j = #sqrt(simga_0**2 + j*sigma_bar**2) spacing is the layer spacing.'''
        #self.pw = 0.0334
        #density of water
        #f = np.zeros(h.size,dtype=np.complex128)
        #hkl = (self.refHKLTransform * np.vstack((h[mask],k[mask],l[mask]))).A
        #l = hkl[2]
        #Q_cart2 = ((self._BMatrix * hkl).A)**2
        #Q2 = np.sum(Q_cart2,axis=0) #squared!!!
        
        #for i in range(2):
        #    f[:] +=  self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
        #    for j in range(5):
        #        f += self.f[i][j]*np.exp(- self.f[i][j+5]*Q2)
        #zwater = WaterModel.f0_water_fun(Q) + self.f #scattering factor of water 
        #sigma_bar = self.basis[3]
        #spacing = self.basis[1]
        l_masked = l[mask]
        zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
        
        
        
        #f = np.empty_like(l_masked,dtype=np.complex128)
        
        if self.type == 'step':
            raise NotImplementedError("This water type uses data, that is not distributed with this version of datautils.")
            #fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            #F_wat = 1j*self.pw*fwater*self.volume*np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(2*np.pi*l_masked)
        elif self.type == 'layered':
            raise NotImplementedError("This water type uses data, that is not distributed with this version of datautils.")
            #fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            #f_w_layer = fwater*self.pw*self.volume*d_layering
            #if A0 == 1.:
            #    sf = np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(1. - np.exp(-2*(np.pi * sigma_bar * l_masked)**2)*np.exp(2j*np.pi * l_masked * d_layering))
            #    F_wat = f_w_layer*sf
            #else:
            #    sigma_1 = np.sqrt(sigma_0**2 + sigma_bar**2)
            #    sf = np.exp(-2*(np.pi * sigma_1 * l_masked)**2)/(1. - np.exp(-2*(np.pi * sigma_bar * l_masked)**2)*np.exp(2j*np.pi * l_masked * d_layering))
            #    F_wat = f_w_layer*( A0*np.exp(-0.5*(sigma_0*l_masked*b[2])**2) + sf*np.exp(2j*np.pi * l_masked * d_layering))
        elif self.type == 'layered_O':
            #fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            
            f = np.empty_like(l_masked,dtype=np.complex128)
            f[:] =  self.f[0][10] + self.f[0][11] + 1j*self.f[0][12]
            for j in range(5):
                f += self.f[0][j]*np.exp(- self.f[0][j+5]*(l_masked*b[2])**2)
            #f *= (10/8)
            
            f_w_layer = f*self.pw*self.volume*d_layering
            
            if A0 == 1.:
                sf = np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(1. - np.exp(-2*(np.pi * sigma_bar * l_masked)**2)*np.exp(2j*np.pi * l_masked * d_layering))
                F_wat = f_w_layer*sf
            else:
                sigma_1 = np.sqrt(sigma_0**2 + sigma_bar**2)
                sf = np.exp(-2*(np.pi * sigma_1 * l_masked)**2)/(1. - np.exp(-2*(np.pi * sigma_bar * l_masked)**2)*np.exp(2j*np.pi * l_masked * d_layering))
                F_wat = f_w_layer*( A0*np.exp(-0.5*(sigma_0*l_masked*b[2])**2) + sf*np.exp(2j*np.pi * l_masked * d_layering))
            
        elif self.type == '1layer':
            raise NotImplementedError("This water type uses data, that is not distributed with this version of datautils.")
            #fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            #rho, mu_offset = self._1layer_firstGauss()
            #F_wat = 1j*self.pw*fwater*self.volume*np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(2*np.pi*l_masked)
            #F_wat = F_wat + A0*fwater*np.exp(-0.5*(sigma_0*l_masked*b[2])**2)*np.exp(2j*np.pi * l_masked * mu_offset)
        
        elif self.type == '1layer_O':
            rho, mu_offset = self._1layer_firstGauss(0.00)
            f = np.empty_like(l_masked,dtype=np.complex128)
            f[:] =  self.f[0][10] + self.f[0][11] + 1j*self.f[0][12]
            for j in range(5):
                f += self.f[0][j]*np.exp(- self.f[0][j+5]*(l_masked*b[2])**2)
                
            f *= (10/8)
            F_wat = 1j*self.pw*f*self.volume*np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(2*np.pi*l_masked)
                            
            F_wat = F_wat + A0*f*np.exp(-0.5*(sigma_0*l_masked*b[2])**2)*np.exp(2j*np.pi * l_masked * mu_offset)
            
        elif self.type == '1layer_s_O':
            rho, mu_offset = self._1layer_firstGauss(0.0)
            f = np.empty_like(l_masked,dtype=np.complex128)
            f[:] =  self.f[2][10] + self.f[2][11] + 1j*self.f[2][12] # self.f[2] is O2-
            for j in range(5):
                f += self.f[2][j]*np.exp(- self.f[2][j+5]*(l_masked*b[2])**2)
                
            F_wat = 1j*self.pw*f*self.volume*np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(2*np.pi*l_masked)
                            
            F_wat = F_wat + A0*f*np.exp(-0.5*(sigma_0*l_masked*b[2])**2)
            
        else:
            warnings.warn("unknown water type %s !" % self.type)
            F_wat = np.zeros_like(l_masked)

        F_wat *= np.exp(2j*np.pi * l_masked * zpos) # position of water
        
        #temp = np.exp(-0.5*(Q*sigma_bar)**2)*np.exp(1j*Q*spacing) 
        #flayer = spacing*(zwater)*pw*self.uc_area*np.exp(-0.5*(Q*sigma_0)**2)
        #F = (np.exp(2j*np.pi*l[mask]*zpos)*(flayer/(1-temp)))/self.volume
        
        F_water = np.zeros_like(l,dtype=np.complex128)
        F_water[mask] = F_wat
        return F_water/self.volume
    
    def setWaterParameters(self,zpos,layer_spacing,sigma_0,sigma_bar,A0):
        self.basis = np.array([zpos,layer_spacing,sigma_0,sigma_bar,A0])

    def setEnergy(self,E):
        self._E = E
        self.lookupScatteringFactors(E)

    def lookupScatteringFactors(self,E):        
        self.f = np.empty((3,13),dtype=np.float64)
        for i,name in enumerate(['O','H','O2-']):
            self.f[i,:11] = readWaasmaier(name)
            self.f[i,11:] = readDispersion(name,E)
        f1,f2 = UnitCell.special_formfactors['H2O'][1](E)
        self.wat_dispersion = f1 + 1j*f2

    def addFitParameter(self,indexarray,limits=(-np.inf,np.inf),**keyargs):
        """
        

        Parameters
        ----------
        indexarray : TYPE
            DESCRIPTION.
        limits : list, optional
            list with [lower, upper] bounds for the fit. The default is [-np.inf,np.inf].

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        int
            internal id of the parameter.

        """
        
        par = indexarray
        try:
            parindexes = np.array([p if isinstance(p, (np.integer,int)) else WaterModel.parameterLookup[p] for p in np.atleast_1d(par)],dtype=np.intp)
        except Exception as e:
            raise ValueError("Invalid atom parameter name %s." % par) from e

        if 'name' in keyargs:
            name = keyargs['name']
        else:
            parameternames = tuple([WaterModel.parameterLookup_inv[n] for n in parindexes])
            name = self.name + " " + " ".join([p for p in parameternames])
        
        if name in self.fitparnames:
            raise ValueError("Absolute fit parameter %s already exists." % name)
        
        indexarray = parindexes
        try:
            curr_val = self.basis[indexarray]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e
        
        if not limits[0] <= np.mean(curr_val) <= limits[1]:
            raise ValueError("start parameter is not within limits! Is not: %s <= %s <= %s" % (limits[0],np.mean(self.basis[indexarray]), limits[1]))
        prior = keyargs.get('prior', None)
        keyargs['name'] = name
        keyargs['prior'] = prior
        
        par = Parameter(name,indexarray, ParameterType.ABSOLUTE, limits, prior, keyargs)
        self.parameters['absolute'].append(par)
        return par
    
    def addRelParameter(self,indexarray,factors,limits=(-np.inf,np.inf),**keyargs):
        par = indexarray
        try:
            parindexes = np.array([p if isinstance(p, (np.integer,int)) else UnitCell.parameterLookup[p] for p in np.atleast_1d(par)],dtype=np.intp)
        except Exception as e:
            raise ValueError("Invalid atom parameter name %s." % par) from e
        
        if 'name' in keyargs:
            name = keyargs['name']
        else:
            parameternames = tuple([WaterModel.parameterLookup_inv[n] + '_r' for n in parindexes])
            name = self.name + " " + " ".join([p for p in parameternames])
        
        if name in self.fitparnames:
            raise ValueError("Relative fit parameter %s already exists." % name)
        
        prior = keyargs.get('prior', None)
        keyargs['name'] = name
        keyargs['prior'] = prior
        
        factors = np.atleast_1d(factors)
        try:
            curr_val = self.basis[parindexes]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e
        
        if curr_val.shape != factors.shape:
            raise ValueError("Number of basis parameters does not match number of factors.")

        par = Parameter(name,parindexes, ParameterType.RELATIVE, limits, prior, keyargs, factors)
        self.parameters['relative'].append(par)
        return par
    
    def showFitparameters(self):
        print(self.fitparameterList())
        
    def parameter_list(self):
        return self.parameters['absolute'] + self.parameters['relative']
    
    @property
    def fitparnames(self):
        #pars = [self.fitparToStr(i,True) for i in range(len(self.fitparameters))] + [self.fitparToStr(i + len(self.fitparameters),True) for i in range(len(self.relfitparam))]
        #pars = 
        return [p.name for p in self.parameters['absolute']] + [p.name for p in self.parameters['relative']]
    
    @property
    def priors(self):
        priorlist = []
        for par in (self.parameters['absolute'] + self.parameters['relative']):
            if par.prior is None:
                if tuple(par.limits) == (-np.inf,np.inf):
                    raise Exception(f"Infinite range given as prior for parameter {par.name}. You must provide an explicit prior.")
                else:
                    priorlist.append(tuple(par.limits)) #has to be converted to correct distribution in CTROptimizer
            else:
                priorlist.append(par.prior) # real prior distribution
        return priorlist
    
    def getInitialParameters(self, force_recalculate=False):
        x0, _, _ = self.getStartParamAndLimits(force_recalculate)
        return x0         
        
    def getStartParamAndLimits(self, force_recalculate=False):
        #if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        try: 
            mismatch = not np.allclose(self._basis_parvalues, self.basis, equal_nan=True)
        except Exception:
            mismatch = True
        recalculate = force_recalculate or mismatch
        abspar = self.parameters['absolute']
        relpar = self.parameters['relative']
        
        x0 = []
        lower = []
        upper = []
        for par in abspar:
            if recalculate or par.value is None:
                x0.append(np.mean(self.basis[par.indices]))
            else:
                x0.append(par.value) 
            lower.append(par.limits[0])
            upper.append(par.limits[1])
        
        for par in relpar:
            if recalculate or par.value is None:    
                x0.append(np.mean( (self.basis[par.indices] - self.basis_0[par.indices]) / par.factors))
            else:
                x0.append(par.value) 
            lower.append(par.limits[0])
            upper.append(par.limits[1])
        
        return np.array(x0), np.array(lower), np.array(upper)
    
    
    def fitparameterList(self,verbosity=3):
        if verbosity > 1:
            outstr = "{}".format("# Direct fit parameters:\n")
        else:
            outstr = ""
            
        if len(self.parameters['absolute']) > 0:
            outstr += "{:10} {:20} {:10} {:10}".format("Id","parameters","Value","Limits")
            
            for par in self.parameters['absolute']:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No direct fit parameters"
        
        outstr += "\n{}".format("# Relative fit parameters:\n")
        
        if len(self.parameters['relative']) > 0:
            
            outstr += "{:10} {:20} {:10} {:10} {:10}".format("Id","parameters","Value","Factors","Limits")
            
            for par in self.parameters['relative']:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No rel. parameters"
        
        return outstr
    
    def fitparToStr(self, par):
        if self.basis_0 is None:
            basis_0 = np.copy(self.basis)
        else:
            basis_0 = self.basis_0

        if par.kind & ParameterType.ABSOLUTE: 
            val = np.mean(self.basis[par.indices])
            parameternames = tuple([WaterModel.parameterLookup_inv[n] for n in par.indices])
            return "{:10} {:15} {:20} {:10} {:10}".format(par.name,str(parameternames),val,str(par.limits))
        elif par.kind & ParameterType.RELATIVE:
            val = np.mean( (self.basis[par.indices] - basis_0[par.indices]) / par.factors)
            parameternames = tuple([WaterModel.parameterLookup_inv[n] + '_r' for n in par.indices])


            return "{:10} {:20} {:10} {:10} {:10}".format(par.name,str(parameternames),val, str(par.factors) ,str(par.limits))
        else:
            raise ValueError("Unvalid parameter type %s for WaterModel" % par.kind)

    def setFitParameters(self,x):
        x_0 = x[:len(self.parameters['absolute'])]
        x_r = x[len(self.parameters['absolute']):]
        self.basis[:] = self.basis_0
        for val, par in zip(x_0,self.parameters['absolute']):
            self.basis[par.indices] = val
            par.value = val
        for val, par in zip(x_r,self.parameters['relative']):
            self.basis[par.indices] += par.factors*val
            par.value = val
        self._basis_parvalues = np.copy(self.basis)
        
    def setLimits(self, lim):
        x_0 = lim[:len(self.parameters['absolute'])]
        x_r = lim[len(self.parameters['absolute']):]
        for val, par in zip(x_0,self.parameters['absolute']):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0],val[1])
        for val, par in zip(x_r,self.parameters['relative']):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0],val[1])
            
    def setFitErrors(self,errors):
        self.errors = np.full_like(self.basis,np.nan)
        err_0 = errors[:len(self.parameters['absolute'])]
        err_r = errors[len(self.parameters['absolute']):]
        for val, par in zip(err_0,self.parameters['absolute']):
            self.errors[par.indices] = val
            par.error = val
        for val, par in zip(err_r,self.parameters['relative']):
            self.errors[par.indices] = np.nan_to_num(self.errors[par.indices],nan=0.0)
            self.errors[par.indices] += par.factors*val
            par.error = val
        self._errors_parvalues = np.copy(self._errors_parvalues)

    def getFitErrors(self):
        if self.errors is None: 
            raise ValueError("No errors have been set.")
        try: 
            mismatch = not np.allclose(self._errors_parvalues, self.errors, equal_nan=True)
        except Exception:
            mismatch = True
        abspar = self.parameters['absolute']
        relpar = self.parameters['relative']
        err0 = []
        for par in abspar:
            if mismatch:
                err0.append(np.mean(self.errors[par.indices]))
            else:
                err0.append(par.error)
                
        for par in relpar:
            if mismatch:
                err0.append(np.mean( self.errors[par.indices]  / np.abs(par.factors)))
            else:
                err0.append(par.error)
        
        err0 = np.array(err0)
        if np.any(~np.isfinite(err0)):
            raise ValueError("Some errors are non-finite or not set.")
        return err0
            
    def _1layer_firstGauss(self,biaswidth=0.2547):
        p = [ 0.85886939,  1.03841354, -1.35438352]
        zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
        #a,alpha,_,_ = self.getLatticeParameters()
        #uc_area = a[0]*a[1]*np.sin(alpha[2])
        
        rho_0 = A0/(self.volume*np.sqrt(2*np.pi*(sigma_0**2 + (biaswidth/self._a[2])**2)))
        
        A = rho_0/self.pw
        
        mu_rel = np.sqrt(2*np.log(2)) + util.stepdown(np.log10(A),p[0],p[1],p[2])
        mu = mu_rel*np.sqrt(sigma_0**2 + (biaswidth/self._a[2])**2)
        
        return rho_0*(10.+self.wat_dispersion), mu
    
        
        
        
    def zDensity_G(self, z, h,k):
        """
        calculates h,k-th Fourier component of the electron density
        Water models only contribute to the 0,0-th component
        Will return an array with zeros for h, k != 0 
        
        The density is normalized to the surface area of the unit cell

        Parameters
        ----------
        z : 1-d array
            z coordinates in Angstrom, should be equidestant 
            and monotonally increasing to avoid numerical issues with convolutions
        h : float
            h-th component index
        k : float
            k-th component index

        Returns
        -------
        1d- array complex128
            complex h,k-th Fourier component of the electron density 
            in electrons/Angstrom**3

        """    
        
        if abs(h) > 0.001 or  abs(k) > 0.001:
            return np.zeros_like(z)
        zstep = np.diff(z)
        zstep_mean = np.mean(zstep)
        if not np.allclose(zstep,zstep_mean):
            warnings.warn("zDensity: z stepsize is not equal in given z array."
                          "This will result in numerical errors in electron density calculation!")
        zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
        if self.type == 'step':
            rho = 0.5*(10.+self.wat_dispersion)*self.pw*(1. + erf((z - zpos*self._a[2])/(np.sqrt(2)* sigma_0 * self._a[2])))
            rho = gaussian_filter1d(np.abs(rho),0.2547 / zstep_mean).astype(np.complex128) #estimation of water molecular form factor
        elif self.type == 'layered':
            zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
            zrange_waterstructure = np.amax(z)/self._a[2] - zpos # lattice units
            nogausseans_inrange = np.ceil(zrange_waterstructure/d_layering)
            nogausseans = int(nogausseans_inrange + 20) # some additional gausseans to reduce edge effects
            layer_density = (10.+self.wat_dispersion)*self.pw*d_layering*self._a[2]
            rho = np.zeros_like(z,dtype=np.complex128)
            for i in range(nogausseans):
                sigma2_rel = sigma_0**2 + i*sigma_bar**2 + (0.2547/self._a[2])**2
                if i == 0:
                    rho += A0*np.exp(- ((z/self._a[2] - zpos - i*d_layering)**2)/(2 * sigma2_rel) ) / (np.sqrt(2*np.pi * sigma2_rel) * self._a[2]) 
                else:
                    rho += np.exp(- ((z/self._a[2] - zpos - i*d_layering)**2)/(2 * sigma2_rel) ) / (np.sqrt(2*np.pi * sigma2_rel) * self._a[2]) 
            rho *= layer_density
        elif self.type == 'layered_O':
            zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
            zrange_waterstructure = np.amax(z)/self._a[2] - zpos # lattice units
            nogausseans_inrange = np.ceil(zrange_waterstructure/d_layering)
            nogausseans = int(nogausseans_inrange + 5) # some additional gausseans to reduce edge effects
            layer_density = self.pw*d_layering*self._a[2]
            rho = np.zeros_like(z,dtype=np.complex128)
            
            for i in range(nogausseans):
                sigma2_rel = sigma_0**2 + i*sigma_bar**2
                deltaZ2i = sigma2_rel*self._a[2]**2
                
                z_i = (zpos + i*d_layering)*self._a[2]
                
                rho_i = ((self.f[0][10] + self.f[0][11] + 1j*self.f[0][12]) / (np.sqrt(2*np.pi*deltaZ2i))) *np.exp( -0.5 *(  ((z- z_i)**2)/deltaZ2i ))
                for j in range(5):
                    exp_dz = self.f[0][j+5] + 0.5*deltaZ2i
                    rho_i += (self.f[0][j]/(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- (((z- z_i)**2)/(4*exp_dz))  )
                #rho_i *= np.exp(-2j*np.pi*((z- z_i)/self._a[2]))
                if i == 0:
                    rho += A0*rho_i*layer_density
                else:
                    rho += rho_i*layer_density
            
        elif self.type == '1layer':
            sigma2_rel = sigma_0**2 + (0.2547/self._a[2])**2
            rho = 0.5*(10.+self.wat_dispersion)*self.pw*(1. + erf((z - zpos*self._a[2])/(np.sqrt(2)* sigma_0 * self._a[2])))
            rho = gaussian_filter1d(np.abs(rho),0.2547 / zstep_mean).astype(np.complex128) #estimation of water molecular form factor
            rho_layer, mu_offset = self._1layer_firstGauss()
            rho += rho_layer*np.exp(- ((z/self._a[2] - zpos - mu_offset)**2)/(2 * sigma2_rel) )
        elif self.type == '1layer_O':

            

            #rho = 0.5*(10.+self.wat_dispersion)*self.pw*(1. + erf((z - zpos*self._a[2])/(np.sqrt(2)* sigma_0 * self._a[2])))
            #rho = gaussian_filter1d(np.abs(rho),0.2547 / zstep_mean).astype(np.complex128) #estimation of water molecular form factor
            rho_layer, mu_offset = self._1layer_firstGauss(0.00)
            
            deltaZ2i = (sigma_0*self._a[2])**2
            z_i = (zpos + mu_offset)*self._a[2]
            
            rho_i = ((self.f[0][10] + self.f[0][11] + 1j*self.f[0][12]) / (np.sqrt(2*np.pi*deltaZ2i))) *np.exp( -0.5 *(  ((z- z_i)**2)/deltaZ2i ))
            rho_step = util.stepup(z,self.f[0][10] + self.f[0][11] + 1j*self.f[0][12],sigma_0*self._a[2],zpos*self._a[2])
            for j in range(5):
                #exp_dpara = self.f[i][j+5] + 0.5*deltaPara2i
                exp_dz = self.f[0][j+5] + 0.5*deltaZ2i
                eff_sigma = np.sqrt(exp_dz/2)
                rho_step += util.stepup(z,self.f[0][j],eff_sigma,zpos*self._a[2])
                rho_i += (self.f[0][j]/(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- (((z- z_i)**2)/(4*exp_dz))  )
            rho_i *= A0*np.exp(-2j*np.pi*((z- z_i)/self._a[2]))
            
            rho_step *= (10/8)*self.pw
            
            rho = rho_i/self.uc_area + rho_step
        elif self.type == '1layer_s_O':
            rho_layer, mu_offset = self._1layer_firstGauss(0.0)
            
            deltaZ2i = (sigma_0*self._a[2])**2
            z_i = (zpos)*self._a[2]
            
            rho_i = ((self.f[2][10] + self.f[2][11] + 1j*self.f[2][12]) / (np.sqrt(2*np.pi*deltaZ2i))) *np.exp( -0.5 *(  ((z- z_i)**2)/deltaZ2i ))
            rho_step = util.stepup(z,self.f[2][10] + self.f[2][11] + 1j*self.f[2][12],sigma_0*self._a[2],zpos*self._a[2])
            for j in range(5):
                #exp_dpara = self.f[i][j+5] + 0.5*deltaPara2i
                exp_dz = self.f[2][j+5] + 0.5*deltaZ2i
                eff_sigma = np.sqrt(exp_dz/2)
                rho_step += util.stepup(z,self.f[2][j],eff_sigma,zpos*self._a[2])
                rho_i += (self.f[2][j]/(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- (((z- z_i)**2)/(4*exp_dz))  )
            rho_i *= A0*np.exp(-2j*np.pi*((z- z_i)/self._a[2]))
            
            rho_step *= self.pw
            
            rho = rho_i/self.uc_area + rho_step
            
            
        return rho

    def __repr__(self):
        repr_super = super(WaterModel,self).__repr__()
        return "%s\n water model" % repr_super
        
    def waterModelToStr(self,showErrors=True):
        param = self.basis
        if (self.errors is not None) and showErrors:
            err = self.errors
            l = []
            for t in zip(param,err):
                [l.append(ti) for ti in t]
            return "Water layer  (%.5f +- %.5f)  (%.5f +- %.5f)  (%.5f +- %.5f)  (%.4f +- %.4f)  (%.4f +- %.4f)" % tuple(l)
        else:
            return "%.5f     %.5f     %.5f  %.5f  %.5f" % tuple(param)

    def latticeRODStr(self):
        a,alpha,_,_ = self.getLatticeParameters()
        return "return\n%s\n%.4f %.4f %.4f %.4f %.4f %.4f" % (self.type,*a,*np.rad2deg(alpha))
    
    
    def parameterStr(self):
        return self.waterModelToStr() + "\n"
    
    def parameterStrRod(self):
        return self.waterModelToStr(False) + "\n"
            
    def __str__(self):
        st = repr(self) + "\n"
        st += "id  " + WaterModel.parameterOrder + "\n"
        return st + self.parameterStr()
            
    def writeWATfile(self,filename):
        with open(filename,'w') as f:
            st = self.latticeRODStr() + "\n" + self.parameterStrRod()
            f.write(st)
            
    def toRODStr(self):
        return self.latticeRODStr() + "\n" + self.parameterStrRod()
    
    def toStr(self):
        return self.latticeRODStr() + "\n" + WaterModel.parameterOrder + "\n" + self.parameterStr()
    
    def pos(self):
        return self.basis[0]*self._a[2]
        
    def pos_cart_error(self):
        z = self.pos()
        error_zeros = np.nan_to_num(self.errors,0)
        z_error = error_zeros[0]*self._a[2]
        return z, z_error
    
    def width_error(self):
        ds0 = self.basis[2]*self._a[2]
        ds1 = self.basis[3]*self._a[2]
        errds0 = self.errors[2]*self._a[2]
        errds1 = self.errors[3]*self._a[2]
        return ds0, ds1, errds0, errds1
    
    @staticmethod
    def fromWATfile(watfile):
        with open(watfile,'r') as f:
            next_skip_comment(f) # return line
            nline = next_skip_comment(f).split()
            try:
                latticeparams = np.array(nline,dtype=np.float64)
                watertype = 'step' # default behaviour
            except ValueError:
                watertype = nline[0]
                latticeparams = np.array(next_skip_comment(f).split(),dtype=np.float64)
            basis = np.loadtxt(f)
        uc = WaterModel(latticeparams[:3],latticeparams[3:],watertype)
        if basis.size < 5: # compatibility to old versions
            basis = np.concatenate((basis,[1.]))
        uc.basis = basis
        uc.basis_0 = np.copy(uc.basis)
        return uc
    
    @staticmethod
    def fromStr(string):
        with util.StringIO(string) as f:
            #next_skip_comment(f) # return line
            nline = next_skip_comment(f).split()
            try:
                latticeparams = np.array(nline,dtype=np.float64)
                watertype = 'step' # default behaviour
            except ValueError:
                watertype = nline[0]
                latticeparams = np.array(next_skip_comment(f).split(),dtype=np.float64)
            try:
                while True:
                    nline = next_skip_comment(f)
                    sline = nline.split()
                    if sline and not 'zpos/frac' in sline:
                        if '+-' in sline:
                            params = re.findall(r'\(([^)]+)',nline)
                            params_array = np.array([np.array(p.split('+-'),dtype=np.float64) for p in params]).T
                            basis = params_array[0]
                            errors = params_array[1]
                        else:
                            basis = np.array(sline,dtype=np.float64)
                            errors = None
                        if basis.size < 5: # compatibility to old versions
                            basis = np.concatenate((basis,[1.]))
                            if errors is not None:
                                errors = np.concatenate((errors,[np.nan]))
            except StopIteration:
                pass
        uc = WaterModel(latticeparams[:3],latticeparams[3:],watertype)
        uc.basis = basis
        uc.errors = errors
        uc.basis_0 = np.copy(uc.basis)
        return uc


class UnitCell(Lattice):
    
    parameterOrder = "Name   x/frac     y/frac     z/frac     iDW     oDW      occup"
    
    parameterLookup = {'x' : 1, 'y' : 2, 'z' : 3,
                       'iDW' : 4, 'oDW' : 5,
                       'occ' : 6}
    
    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))
    
    # lookuptable for special materials
    # give two functions: 1. for f0 (Q) , 2. lookup for f1 and f2 (E)
    #special_formfactors = OrderedDict([('H2O' , (f0_water_fun, lambda E : estimateDispersionCompound('H2O', E) ) )])
    special_formfactors = OrderedDict()
    #special_eDensity = OrderedDict([('H2O' , eDens_water_fun )])
    special_eDensity = OrderedDict()
    __db = xraydb.get_xraydb()
    if hasattr(__db, 'atomic_symbols'): # XrayDB <= 4.5
        special_onset = len(__db.atomic_symbols)
    else:
        __elems = __db.get_cache('elements')
        special_onset = len([e.element for e in __elems])
    __db.close()
        
    special_numbers = OrderedDict(zip(special_formfactors.keys(),np.arange(len(special_formfactors))+special_onset))
    
    def __init__(self,a,alpha,**keyargs):
        super(UnitCell,self).__init__(a,alpha)
        self.basis = np.array([])
        self.names = []
        
        """
        self.fitparameters = []
        self.fitparameters_name = []
        self.relfitparam = []
        self.relfitparam_name = []
        self.fitparlimits = []
        self.relfitlimits = []
        self.relfitparam_prior = []
        self.fitparameter_prior = []
        """
        
        self.parameters = {
            'absolute' : [],
            'relative' : []
        }
        
        self.basis_0 = np.array([])
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None
        if 'name' in keyargs:
            self.name = keyargs['name']
        else:
            self.name = "unnamed"
        self.refHKLTransform = np.identity(3,dtype=np.float64)
        self.refRealTransform = np.identity(3,dtype=np.float64)
        self.coherentDomainMatrix = [np.vstack((np.identity(3).T,np.array([0,0,0]))).T]
        self.coherentDomainOccupancy = [1.]
        self.dw_increase_constraint = np.array([],dtype=np.bool_)
        self._special_formfactors_present = False
        
    def parametersToDict(self):
        d = dict()
        d['basis_0'] = self.basis_0
        d['dw_increase_constraint'] = self.dw_increase_constraint
        for p in self.parameters:
            d[p] = dict()
            for i, param in enumerate(self.parameters[p]):
                d[p]["%s_%s" % (i,param.name)] = param.asdict()
        return d
    
    def parametersFromDict(self, d, override_values=True):
        for p in self.parameters:
            self.parameters[p] = []
            
        self.basis_0[:] = d['basis_0'].astype(np.float64)
        self.dw_increase_constraint = d['dw_increase_constraint'].astype(np.bool_)
        for p in self.parameters:
            for dkey in sorted(d[p].keys()):
                self.parameters[p].append(Parameter(**d[p][dkey]))

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn("Can not apply fit parameter values to this unitcell: %s" % e) 
        else:
            warnings.warn("Basis values not updated from Parameter values. This can cause a mismatch between basis and fitparameter values!")


    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters 
        """
        no_errors = False
        self.basis[:] = self.basis_0
        self.errors = np.full_like(self.basis,np.nan)
        for par in self.parameters['absolute']:
            if par.value is not None:
                self.basis[par.indices] = par.value
            else:
                raise ValueError("Can not set basis values from parameters. Value of Parameter %s is None." % par.name)
            if par.error is not None:
                self.errors[par.indices] = par.error
            else:
                no_errors = True
        for par in self.parameters['relative']:
            if par.value is not None:
                self.basis[par.indices] += par.factors*par.value
            else:
                raise ValueError("Can not set basis values from parameters. Value of Parameter %s is None." % par.name)
            if par.error is not None:
                self.errors[par.indices] = np.nan_to_num(self.errors[par.indices],nan=0.0)
                self.errors[par.indices] += par.factors*par.error
            else:
                no_errors = True

        self._basis_parvalues = np.copy(self.basis)
        if not no_errors:
            self._errors_parvalues = np.copy(self.errors)
    
    def clearParameters(self):
        for p in self.parameters:
            self.parameters[p] = []
        self._basis_parvalues = None
        self._errors_parvalues = None

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
        self.refRealTransform = self.R_mat_inv @ rotMatrix @ uc.R_mat 
        self.refHKLTransform = self.B_mat_inv @ rotMatrix @ uc.B_mat
        
    def _test_special_formfactors(self):
        for name in self.names:
            if name in UnitCell.special_formfactors:
                self._special_formfactors_present = True
                return True
        else:
            self._special_formfactors_present = False
            return False
    
    def addAtom(self,element_or_param, xyz_rel=None, iDW=None, oDW=None, occu=None):
        if xyz_rel is not None:
            if not isinstance(element_or_param, str):
                raise ValueError("You must provide a atom symbol. Got: %s" % element_or_param)
            atom_no = atomic_number(element_or_param)

            parameters = np.array([atom_no,
                                   *xyz_rel,
                                   iDW,
                                   oDW,
                                   occu],dtype=np.float64)
            self.names.append(element_or_param) 
        else:
            if not isinstance(element_or_param[0], str):
                raise ValueError("You must provide a atom symbol. Got: %s" % element_or_param[0])
            self.names.append(element_or_param[0])
            atom_no = atomic_number(element_or_param[0])
            element_or_param[0] = atom_no
            parameters = np.array(element_or_param,dtype=np.float64)
            
        if self.basis.size > 0:
            self.basis = np.vstack([self.basis,parameters])
        else:
            self.basis = np.array([parameters])
        if self.basis_0.size > 0:
            self.basis_0 = np.vstack([self.basis_0,parameters])
        else:
            self.basis_0 = np.array([parameters])
        self.errors = None
        self.dw_increase_constraint = np.insert(self.dw_increase_constraint, len(self.dw_increase_constraint), True)
        self._test_special_formfactors()
        return
    
    def insertAtom(self,index,element_or_param, xyz_rel=None, iDW=None, oDW=None, occu=None):
        if xyz_rel is not None:
            if not isinstance(element_or_param, str):
                raise ValueError("You must provide a atom symbol. Got: %s" % element_or_param)
            atom_no = atomic_number(element_or_param)

            parameters = np.array([atom_no,
                                   *xyz_rel,
                                   iDW,
                                   oDW,
                                   occu],dtype=np.float64)
            self.names.insert(element_or_param,index) 
        else:
            if not isinstance(element_or_param[0], str):
                raise ValueError("You must provide a atom symbol. Got: %s" % element_or_param[0])
            self.names.insert(element_or_param[0],index)
            atom_no = atomic_number(element_or_param[0])
            element_or_param[0] = atom_no
            parameters = np.array(element_or_param,dtype=np.float64)
            
        if self.basis.size > 0:
            self.basis = np.insert(self.basis,index,parameters,axis=0)
        else:
            self.basis = np.array([parameters])
        if self.basis_0.size > 0:
            self.basis_0 = np.insert(self.basis_0,index,parameters,axis=0)
        else:
            self.basis_0 = np.array([parameters])
        self.dw_increase_constraint = np.insert(self.dw_increase_constraint, index, True)
        self.errors = None
        self._test_special_formfactors()
        #newfp = []
        """
        for atomindex, parindex in self.fitparameters:
            try:
                if atomindex >= index:
                    atomindex += 1
                newfp.append((atomindex,parindex))
            except TypeError:
                idxarray = []
                for idx in atomindex:
                    if idx >= index:
                        idx += 1
                    idxarray.append(idx)
                newfp.append((tuple(idxarray),parindex))
        self.fitparameters = newfp
        """
    # in eV
    def setEnergy(self,E):
        """
        for dispersion and absorption correction
        in eV
        """
        self._E = E
        self.lookupScatteringFactors(E)
        
    def addFitParameter(self,indexarray,limits=(-np.inf,np.inf),**keyargs):
        """
        

        Parameters
        ----------
        indexarray : TYPE
            DESCRIPTION.
        limits : list, optional
            list with [lower, upper] bounds for the fit. The default is [-np.inf,np.inf].

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        int
            internal id of the parameter.

        """
        
        atoms, par = indexarray
        try:
            parindexes = np.array([p if isinstance(p, (np.integer,int)) else UnitCell.parameterLookup[p] for p in np.atleast_1d(par)],dtype=np.intp)
        except Exception as e:
            raise ValueError("Invalid atom parameter name %s." % par) from e
        atoms = np.atleast_1d(atoms).astype(np.intp)
        
        if 'name' in keyargs:
            name = keyargs['name']
        else:
            parameternames = tuple([UnitCell.parameterLookup_inv[n] for n in parindexes])
            atomnames = tuple([f"{n}_{self.names[n]}" for n in atoms])
            name = self.name + " " + " ".join(["_".join((n,p)) for p, n in zip(parameternames,atomnames)])
        
        if name in self.fitparnames:
            raise ValueError("Absolute fit parameter %s already exists." % name)
        
        indexarray = atoms, parindexes
        try:
            curr_val = self.basis[indexarray]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e
        
        if not limits[0] <= np.mean(curr_val) <= limits[1]:
            raise ValueError("start parameter is not within limits! Is not: %s <= %s <= %s" % (limits[0],np.mean(self.basis[indexarray]), limits[1]))
        prior = keyargs.get('prior', None)
        keyargs['name'] = name
        keyargs['prior'] = prior
        
        par = Parameter(name,indexarray, ParameterType.ABSOLUTE, limits, prior, keyargs)
        self.parameters['absolute'].append(par)
        return par
    
    def addRelParameter(self,indexarray,factors,limits=(-np.inf,np.inf),**keyargs):
        #if not limits[0] <= np.mean(self.basis[indexarray]) <= limits[1]:
        #    raise ValueError("start parameter is not within limits!")
        
        atoms, par = indexarray
        try:
            parindexes = np.array([p if isinstance(p, (np.integer,int)) else UnitCell.parameterLookup[p] for p in np.atleast_1d(par)],dtype=np.intp)
        except Exception as e:
            raise ValueError("Invalid atom parameter name %s." % par) from e
        atoms = np.atleast_1d(atoms).astype(np.intp)
        
        if 'name' in keyargs:
            name = keyargs['name']
        else:
            parameternames = tuple([UnitCell.parameterLookup_inv[n] + '_r' for n in parindexes])
            atomnames = tuple([f"{n}_{self.names[n]}" for n in atoms])
            name = self.name + " " + " ".join(["_".join((n,p)) for p, n in zip(parameternames,atomnames)])
        
        if name in self.fitparnames:
            raise ValueError("Relative fit parameter %s already exists." % name)
        
        prior = keyargs.get('prior', None)
        keyargs['name'] = name
        keyargs['prior'] = prior
        
        factors = np.atleast_1d(factors)
        try:
            curr_val = self.basis[(atoms, parindexes)]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e
        
        if curr_val.shape != factors.shape:
            raise ValueError("Number of basis parameters does not match number of factors.")

        par = Parameter(name,(atoms, parindexes), ParameterType.RELATIVE, limits, prior, keyargs, factors)
        self.parameters['relative'].append(par)
        return par
    
    def showFitparameters(self):
        print(self.fitparameterList())
        
    def parameter_list(self):
        return self.parameters['absolute'] + self.parameters['relative']
    
    @property
    def fitparnames(self):
        #pars = [self.fitparToStr(i,True) for i in range(len(self.fitparameters))] + [self.fitparToStr(i + len(self.fitparameters),True) for i in range(len(self.relfitparam))]
        #pars = 
        return [p.name for p in self.parameters['absolute']] + [p.name for p in self.parameters['relative']]
    
    @property
    def priors(self):
        priorlist = []
        for par in (self.parameters['absolute'] + self.parameters['relative']):
            if par.prior is None:
                if tuple(par.limits) == (-np.inf,np.inf):
                    raise Exception(f"Infinite range given as prior for parameter {par.name}. You must provide an explicit prior.")
                else:
                    priorlist.append(tuple(par.limits)) #has to be converted to correct distribution in CTROptimizer
            else:
                priorlist.append(par.prior) # real prior distribution
        return priorlist
    
    def fitparameterList(self,verbosity=3):
        if verbosity > 1:
            outstr = "{}".format("# Direct fit parameters:\n")
        else:
            outstr = ""
            
        if len(self.parameters['absolute']) > 0:
            outstr += "{:10} {:15} {:20} {:10} {:10}".format("Id","atoms","parameters","Value","Limits")
            
            for par in self.parameters['absolute']:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No direct fit parameters"
        
        outstr += "\n{}".format("# Relative fit parameters:\n")
        if len(self.parameters['relative']) > 0:
            outstr += "{:10} {:15} {:20} {:10} {:10} {:10}".format("Id","atoms","parameters","Value","Factors","Limits")
            
            for par in self.parameters['relative']:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No rel. parameters"
        
        return outstr
    
    def fitparToStr(self, par):
        if self.basis_0 is None:
            basis_0 = np.copy(self.basis)
        else:
            basis_0 = self.basis_0
        
        if par.kind & ParameterType.ABSOLUTE: 
            val = np.mean(self.basis[par.indices])
            #if isinstance(par.indices[1], (np.integer,int)):
            #    parameternames = (UnitCell.parameterLookup_inv[par.indices[1]],)
            #    atoms = (f"{par.indices[0]}_{self.names[par.indices[0]]}",)
            #else:
            parameternames = tuple([UnitCell.parameterLookup_inv[n] for n in par.indices[1]])
            atoms = tuple([f"{n}_{self.names[n]}" for n in par.indices[0]])
            return "{:10} {:15} {:20} {:10} {:10}".format(par.name,str(atoms),str(parameternames),val,str(par.limits))
        
        
        elif par.kind & ParameterType.RELATIVE:
            val = np.mean( (self.basis[par.indices] - basis_0[par.indices]) / par.factors)
            #if isinstance(par.indices[1], (np.integer,int)):
            #    parameternames = (UnitCell.parameterLookup_inv[par.indices[1]] + '_r',)
            #    atoms = (f"{par.indices[0]}_{self.names[par.indices[0]]}",)
            #else:
            parameternames = tuple([UnitCell.parameterLookup_inv[n] + '_r' for n in par.indices[1]])
            atoms = tuple([f"{n}_{self.names[n]}" for n in par.indices[0]])
            return "{:10} {:15} {:20} {:10} {:10} {:10}".format(par.name ,str(atoms),str(parameternames),val, str(par.factors) ,str(par.limits))
            
        else:
            raise ValueError("Unvalid parameter type %s for UnitCell" % par.kind)

    #def fixInitialBasis(self):
    #    self.basis_0 = np.copy(self.basis)

    def getInitialParameters(self, force_recalculate=False):
        x0, _, _ = self.getStartParamAndLimits(force_recalculate)
        return x0   
        
    def getStartParamAndLimits(self, force_recalculate=False):
        #if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        try: 
            mismatch = not np.allclose(self._basis_parvalues, self.basis, equal_nan=True)
        except Exception:
            mismatch = True
        recalculate = force_recalculate or mismatch
        abspar = self.parameters['absolute']
        relpar = self.parameters['relative']
        
        x0 = []
        lower = []
        upper = []
        for par in abspar:
            if recalculate or par.value is None:
                x0.append(np.mean(self.basis[par.indices]))
            else:
                x0.append(par.value) 
            lower.append(par.limits[0])
            upper.append(par.limits[1])
        
        for par in relpar:
            if recalculate or par.value is None:    
                x0.append(np.mean( (self.basis[par.indices] - self.basis_0[par.indices]) / par.factors))
            else:
                x0.append(par.value) 
            lower.append(par.limits[0])
            upper.append(par.limits[1])
        
        return np.array(x0), np.array(lower), np.array(upper)
        
    def setFitParameters(self,x):
        x_0 = x[:len(self.parameters['absolute'])]
        x_r = x[len(self.parameters['absolute']):]
        self.basis[:] = self.basis_0
        for val, par in zip(x_0,self.parameters['absolute']):
            self.basis[par.indices] = val
            par.value = val
        for val, par in zip(x_r,self.parameters['relative']):
            self.basis[par.indices] += par.factors*val
            par.value = val
        self._basis_parvalues = np.copy(self.basis)
        
    def setLimits(self, lim):
        x_0 = lim[:len(self.parameters['absolute'])]
        x_r = lim[len(self.parameters['absolute']):]
        for val, par in zip(x_0,self.parameters['absolute']):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0],val[1])
        for val, par in zip(x_r,self.parameters['relative']):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0],val[1])

    def setFitErrors(self,errors):
        self.errors = np.full_like(self.basis,np.nan)
        err_0 = errors[:len(self.parameters['absolute'])]
        err_r = errors[len(self.parameters['absolute']):]
        for val, par in zip(err_0,self.parameters['absolute']):
            self.errors[par.indices] = val
            par.error = val
        for val, par in zip(err_r,self.parameters['relative']):
            self.errors[par.indices] = np.nan_to_num(self.errors[par.indices],nan=0.0)
            self.errors[par.indices] += np.abs(par.factors)*val
            par.error = val
        self._errors_parvalues = np.copy(self.errors)

    def getFitErrors(self):
        if self.errors is None: 
            raise ValueError("No errors have been set.")
        try: 
            mismatch = not np.allclose(self._errors_parvalues, self.errors, equal_nan=True)
        except Exception:
            mismatch = True
        abspar = self.parameters['absolute']
        relpar = self.parameters['relative']
        err0 = []
        for par in abspar:
            if mismatch:
                err0.append(np.mean(self.errors[par.indices]))
            else:
                err0.append(par.error)
                
        for par in relpar:
            if mismatch:
                err0.append(np.mean( self.errors[par.indices]  / np.abs(par.factors)))
            else:
                err0.append(par.error)
        
        err0 = np.array(err0)
        if np.any(~np.isfinite(err0)):
            raise ValueError("Some errors are non-finite or not set.")
        return err0
        

    # returns the structure factor of the unit cell
    # h,k,l have to be 1d arrays 
    def F_uc_bulk(self,h,k,l,atten=0):
        if HAS_NUMBA_ACCEL:
            h,k,l = _ensure_contiguous(h,k,l, testOnly=False, astype=np.float64)
            F = _CTRcalc_accel.unitcell_F_uc_bulk(h,
                    k,
                    l,
                    atten,
                    self.basis,
                    self.f,
                    self.refHKLTransform,
                    self.B_mat,
                    self.R_mat,
                    self.R_mat_inv,
                    np.asarray(self.coherentDomainMatrix),
                    np.asarray(self.coherentDomainOccupancy),
                    self.volume)
            return F
        else:
            F = np.zeros(h.size,dtype=np.complex128)
            f = np.zeros(h.size,dtype=np.complex128)
            hkl = self.refHKLTransform @ np.vstack((h,k,l))
            Q_cart2 = (self.B_mat @ hkl)**2
            Q_para2 = np.sum(Q_cart2[:2],axis=0) #squared!!!
            Q_perp2 = Q_cart2[2] #squared!!!
            Q2 = Q_para2 + Q_perp2 #squared!!!
            
            domainmatrix = [self.R_mat_inv @ mat[:,:-1] @ self.R_mat for mat in self.coherentDomainMatrix]
            
            for i in range(len(self.basis)):
                f[:] =  self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
                for j in range(5):
                    f += self.f[i][j]*np.exp(- self.f[i][j+5]*Q2)
                f *= np.exp(- (self.basis[i][4] * Q_para2 + self.basis[i][5] * Q_perp2)/ (16*np.pi**2))
                f *= self.basis[i][6]
                for mat, weight, eff_mat in zip(self.coherentDomainMatrix,self.coherentDomainOccupancy, domainmatrix):
                    xyz_rel = np.dot(eff_mat,self.basis[i][1:4]) + mat[:,-1]
                    F += weight * f * np.exp(2j*np.pi * np.sum(hkl.T * xyz_rel,axis=1) ) * math.exp(atten*xyz_rel[2])
            return F/self.volume
    
    # returns the structure factor of the unit cell
    # h,k,l have to be 1d arrays 
    def F_uc_bulk_direct(self,h,k,l,atten=0):
        F = np.zeros(h.size,dtype=np.complex128)
        f = np.zeros(h.size,dtype=np.complex128)
        hkl =  np.vstack((h,k,l))
        Q_cart2 = (self.B_mat @ hkl)**2
        Q_para2 = np.sum(Q_cart2[:2],axis=0) #squared!!!
        Q_perp2 = Q_cart2[2] #squared!!!
        Q2 = Q_para2 + Q_perp2 #squared!!!
        
        domainmatrix = [self.R_mat_inv @ mat[:,:-1] @ self.R_mat for mat in self.coherentDomainMatrix]
        
        for i in range(len(self.basis)):
            f[:] =  self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
            for j in range(5):
                f += self.f[i][j]*np.exp(- self.f[i][j+5]*Q2)
            f *= np.exp(- (self.basis[i][4] * Q_para2 + self.basis[i][5] * Q_perp2)/ (16*np.pi**2))
            f *= self.basis[i][6]
            for mat, weight, eff_mat in zip(self.coherentDomainMatrix,self.coherentDomainOccupancy,domainmatrix):
                xyz_rel = np.dot(eff_mat,self.basis[i][1:4]) + mat[:,-1]
                F += weight * f * np.exp(2j*np.pi * np.sum(hkl.T * xyz_rel,axis=1) ) * math.exp(atten*xyz_rel[2])
        return F/self.volume
    
    def F_uc(self,h,k,l):
        if HAS_NUMBA_ACCEL and not self._special_formfactors_present:
            h,k,l = _ensure_contiguous(h,k,l, testOnly=False, astype=np.float64)
            F = _CTRcalc_accel.unitcell_F_uc(h,
                    k,
                    l,
                    self.basis,
                    self.f,
                    self.refHKLTransform,
                    self.B_mat,
                    self.R_mat,
                    self.R_mat_inv,
                    np.asarray(self.coherentDomainMatrix),
                    np.asarray(self.coherentDomainOccupancy),
                    self.volume)
            return F
        else:
            F = np.zeros(h.size,dtype=np.complex128)
            f = np.zeros(h.size,dtype=np.complex128)
            hkl = self.refHKLTransform @ np.vstack((h,k,l))
            Q_cart2 = (self.B_mat @ hkl)**2
            Q_para2 = np.sum(Q_cart2[:2],axis=0) #squared!!!
            Q_perp2 = Q_cart2[2] #squared!!!
            Q2 = Q_para2 + Q_perp2 #squared!!!
            
            domainmatrix = [self.R_mat_inv @ mat[:,:-1] @ self.R_mat for mat in self.coherentDomainMatrix]
                
            for i,name in zip(range(len(self.basis)),self.names):
                if name in UnitCell.special_formfactors:
                    f[:] = UnitCell.special_formfactors[name][0](np.sqrt(Q2)) + self.f[i][11] + 1j*self.f[i][12]
                else:
                    f[:] =  self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
                    for j in range(5):
                        f += self.f[i][j]*np.exp(- self.f[i][j+5]*Q2)
                f *= np.exp(- (self.basis[i][4] * Q_para2 + self.basis[i][5] * Q_perp2)/ (16*np.pi**2))
                f *= self.basis[i][6]
                for mat, weight, eff_mat in zip(self.coherentDomainMatrix,self.coherentDomainOccupancy, domainmatrix):
                    xyz_rel = np.dot(eff_mat,self.basis[i][1:4]) + mat[:,-1]
                    F += weight * f * np.exp(2j*np.pi * np.sum(hkl.T * xyz_rel,axis=1) )
            return F/self.volume
    
    def F_bulk(self,h,k,l,atten=0):
        if HAS_NUMBA_ACCEL:
            h,k,l = _ensure_contiguous(h,k,l, testOnly=False, astype=np.float64)
            F = _CTRcalc_accel.unitcell_F_bulk(h,
                    k,
                    l,
                    atten,
                    self.basis,
                    self.f,
                    self.refHKLTransform,
                    self.B_mat,
                    self.R_mat,
                    self.R_mat_inv,
                    np.asarray(self.coherentDomainMatrix),
                    np.asarray(self.coherentDomainOccupancy),
                    self.volume)
            return F
        else:
            hkl = self.refHKLTransform @ np.vstack((h,k,l))
            Fuc = self.F_uc_bulk_direct(*hkl,atten)
            return Fuc/(1- np.exp(- 2j*np.pi * l - atten ))
    
    SQRT2pi = np.sqrt(2*np.pi)
    def zDensity_G(self, z, h,k):
        """
        calculates h,k-th Fourier component of the electron density of the unit cell 
        i.e. 0,0-th component is the commonly used z-projected electron density
        
        The density is normalized to the surface area of the unit cell

        Parameters
        ----------
        z : 1-d array
            z coordinates in Angstrom, should be equidestant 
            and monotonally increasing to avoid numerical issues with convolutions
        h : float
            h-th component index
        k : float
            k-th component index

        Returns
        -------
        1d- array complex128
            complex h,k-th Fourier component of the electron density 
            in electrons/Angstrom**3
            calculate the absolute value to get the electron density

        """
        hkl = (self.refHKLTransform @ np.array([h,k,0.])).flatten()
        Qpara2 = np.sum(np.dot(self.B_mat , hkl)**2)
        a,alpha,_,_ = self.getLatticeParameters()
        #uc_area = a[0]*a[1]*np.sin(alpha[2])
        #dispersion = self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
        zstep = np.diff(z)
        zstep_mean = np.mean(zstep)
        if not np.allclose(zstep,zstep_mean):
            warnings.warn("zDensity: z stepsize is not equal in given z array."
                          "This will result in numerical errors in electron density calculation!")
        
        rho = np.zeros_like(z,dtype=np.complex128)
        rho_i = np.empty_like(z,dtype=np.complex128)
        
        domainmatrix = [self.R_mat_inv @ mat[:,:-1] @ self.R_mat for mat in self.coherentDomainMatrix]
        
        for i,name in zip(range(len(self.basis)),self.names):
            deltaZ2i = self.basis[i][5]/(8*(np.pi)**2)
            deltaPara2i = self.basis[i][4]/(8*(np.pi)**2)
            for mat, weight, eff_mat in zip(self.coherentDomainMatrix,self.coherentDomainOccupancy, domainmatrix):
                xyz_rel = eff_mat @ self.basis[i][1:4] + mat[:,-1]
                z_i = xyz_rel[2]*self._a[2]
                y_i_frac = xyz_rel[1]
                x_i_frac = xyz_rel[0]
                if name in UnitCell.special_eDensity:
                    rho_i[:] = (self.f[i][11] + 1j*self.f[i][12])/ (np.sqrt(2*np.pi*deltaZ2i))
                    rho_i *= np.exp( -0.5 *( deltaPara2i * Qpara2 +  ((z- z_i)**2)/deltaZ2i ))
                    rho_i += gaussian_filter1d(UnitCell.special_eDensity[name](z- z_i)*np.exp( -0.5 *( deltaPara2i * Qpara2)),np.sqrt(deltaZ2i)/ zstep_mean)
                else:                
                    rho_i[:] = (self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]) / (np.sqrt(2*np.pi*deltaZ2i))
                    rho_i *= np.exp( -0.5 *( deltaPara2i * Qpara2 +  ((z- z_i)**2)/deltaZ2i ))
                    for j in range(5):
                        exp_dpara = self.f[i][j+5] + 0.5*deltaPara2i
                        exp_dz = self.f[i][j+5] + 0.5*deltaZ2i
                        rho_i += (self.f[i][j]/(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- exp_dpara * Qpara2  - (((z- z_i)**2)/(4*exp_dz))  )
                        
                rho_i *= np.exp(-2j*np.pi*(h*x_i_frac + k*y_i_frac + (z- z_i)/self._a[2]))
                rho += rho_i*self.basis[i][6] * weight
        
        return rho/self.uc_area
    
    
    def zDensity_G_asbulk(self, z, h,k, noUC=30):
        hkl = (self.refHKLTransform @ np.array([h,k,0.])).flatten()
        Qpara2 = np.sum(np.dot(self.B_mat , hkl)**2)
        a,alpha,_,_ = self.getLatticeParameters()
        #uc_area = a[0]*a[1]*np.sin(alpha[2])
        #dispersion = self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
        rho = np.zeros_like(z,dtype=np.complex128)
        rho_i = np.empty_like(z,dtype=np.complex128)
        for no in range(noUC):
            noA = no*self._a[2]
            for i in range(len(self.basis)):
                deltaZ2i = self.basis[i][5]/(8*(np.pi)**2)
                deltaPara2i = self.basis[i][4]/(8*(np.pi)**2)
                z_i = self.basis[i][3]*self._a[2]
                y_i_frac = self.basis[i][2]
                x_i_frac = self.basis[i][1]
                
                rho_i[:] = (self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]) / (np.sqrt(2*np.pi*deltaZ2i))
                rho_i *= np.exp( -0.5 *( deltaPara2i * Qpara2 +  ((z- z_i + noA)**2)/deltaZ2i ))
                        
                
                for j in range(5):
                    exp_dpara = self.f[i][j+5] + 0.5*deltaPara2i
                    exp_dz = self.f[i][j+5] + 0.5*deltaZ2i
                    rho_i += (self.f[i][j]/(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- exp_dpara * Qpara2  - (((z- z_i + noA)**2)/(4*exp_dz))  )
                    
                rho_i *= np.exp(-2j*np.pi*(h*x_i_frac + k*y_i_frac + (z- z_i + noA)/self._a[2]))
                rho += rho_i*self.basis[i][6]
        
        return rho/self.uc_area
    
    def lookupScatteringFactors(self,E):
        self.f = np.empty((self.basis.shape[0],13),dtype=np.float64)
        for i,name in enumerate(self.names):
            if name in UnitCell.special_formfactors:
                self.f[i,:11] = 0.
                self.f[i,11:] = UnitCell.special_formfactors[name][1](E)
            else:
                self.f[i,:11] = readWaasmaier(name)
                self.f[i,11:] = readDispersion(name,E)
        return
    
    def plot3d(self,ucx=1,ucy=1,ucz=1,dwon=False,occuon=False,figure=None,translate=np.array([0.,0.,0.]), **keyargs):
        try:
            from mayavi import mlab
        except ImportError:
            warnings.warn("can not import mayavi: 3D plotting not supported")
            return

        if figure is None:
            figure = mlab.figure()
        
        if ucx == 0 or ucy == 0 or ucz == 0:
            raise ValueError("One of ucx,ucy or ucz is zero. Must plot at least one unit cell.")
        
        def color_generator(name):
            return special_elementcolors.get(name, tuple(rgb_array[atomic_number(name)-1]))
            
        elcolors = keyargs.get('color')
        if elcolors is None:
            elcolors = [color_generator(name) for name in self.names]
        elif isinstance(elcolors ,tuple):
            elcolors = [elcolors for params in self.basis]
            
        resolution = keyargs.get('resolution')
        if resolution is None:
            resolution = 25

        for i,params in enumerate(self.basis):
            
            radius = cov_radii_array[atomic_number(self.names[i])-1][2]*2
            #elcolor_c = keyargs.get('color')
            #if elcolor_c is None:
            #    elcolor_c = elements.rgb(int(params[0]))
            if occuon:
                occup = params[6]
            else:
                occup = 1
            x,y,z = params[1:4]
            positions = np.empty((np.abs(ucx*ucy*ucz),3))
            iDW, oDW = DWtoDisorder(params[4:6])
            signx = np.sign(ucx)
            signy = np.sign(ucy)
            signz = np.sign(ucz)
            no = 0
            for xno in range(1,abs(ucx)+1):
                for yno in range(1,abs(ucy)+1):
                    for zno in range(1,abs(ucz)+1):
                        if len(translate.shape) > 1: 
                            x_t, y_t, z_t = np.array([x,y,z]) + translate[:,-1].flatten()
                            positions[no] = np.array([x_t + signx*(xno-1),y_t + signy*(yno-1),z_t + signz*(zno-1)])
                        else:
                            positions[no] = np.array([x + signx*(xno-1),y + signy*(yno-1),z+ signz*(zno-1)]) + translate
                        if random.random() > occup:
                            positions[no] = np.nan
                        no += 1
            position_cart = self.RealspaceMatrix @ positions.T
            if len(translate.shape) > 1:
                position_cart = np.asarray(translate[:,:-1] @ position_cart)

            sigmas = np.empty_like(position_cart)
            sigmas[:2] = iDW
            sigmas[2] = oDW
            
            if dwon:
                position_cart = np.random.normal(position_cart,sigmas)
            mlab.points3d(*position_cart,scale_factor=radius,color=elcolors[i],resolution=resolution,figure=figure)
            
        atomlist = keyargs.get('atomlist')
        if atomlist is not None:
            atomlist.append(self.pos_cart_all(ucx,ucy,ucz,translate))
                
        return figure
    
    
    def pos_cart_all(self,ucx=1,ucy=1,ucz=1,translate=np.array([0.,0.,0.])):
        dt = np.dtype([('name', 'U6'),('x', np.float64), ('y', np.float64), ('z', np.float64)])
        xyz_array = np.empty(np.abs(ucx*ucy*ucz)*len(self.names),dtype=dt)
        position_cart = []
        names = []
        for i,params in enumerate(self.basis):
            x,y,z = params[1:4]
            positions = np.empty((np.abs(ucx*ucy*ucz),3))
            signx = np.sign(ucx)
            signy = np.sign(ucy)
            signz = np.sign(ucz)
            no = 0
            for xno in range(1,abs(ucx)+1):
                for yno in range(1,abs(ucy)+1):
                    for zno in range(1,abs(ucz)+1):
                        if len(translate.shape) > 1: 
                            x_t, y_t, z_t = np.array([x,y,z]) + translate[:,-1].flatten()
                            positions[no] = np.array([x_t + signx*(xno-1),y_t + signy*(yno-1),z_t + signz*(zno-1)])
                        else:
                            positions[no] = np.array([x + signx*(xno-1),y + signy*(yno-1),z+ signz*(zno-1)]) + translate
                        no += 1
            positions_c = self.RealspaceMatrix @ positions.T
            
            if len(translate.shape) > 1:
                positions_c = np.asarray(translate[:,:-1] @ positions_c)
            position_cart.append(positions_c)
            names.append([self.names[i] for k in range(no)])
            
        positions_all = np.concatenate(position_cart,axis=1)
        xyz_array['x'] = positions_all[0]
        xyz_array['y'] = positions_all[1]
        xyz_array['z'] = positions_all[2]
        xyz_array['name'] = np.concatenate(names)
        return xyz_array

    
    def pos_cart(self,atomNo):
        return self.RealspaceMatrix @ self.basis[atomNo][1:4].T
        
    def pos_cart_error(self,atomNo):
        xyz = self.pos_cart(atomNo)
        realMatrix = self.RealspaceMatrix**2
        error_zeros = np.nan_to_num(self.errors,0)
        xyz_error = np.sqrt(realMatrix @ (error_zeros[atomNo][1:4]**2).T)
        return xyz, xyz_error
        
        
    def distanceVectorAtom(self,atomNo1,atomNo2):
        return self.pos_cart(atomNo1) - self.pos_cart(atomNo2)
    
    def bondingLength(self,xyz_atom):
        pos_atom_cart = self.RealspaceMatrix @ xyz_atom.T
        blength = np.empty(self.basis.shape[0],dtype=np.float64)
        for i, params in enumerate(self.basis):
            #occup = params[6]
            xyz = params[1:4]
            pos_cart = self.RealspaceMatrix @ xyz
            #iDW, oDW = DWtoDisorder(params[4:6])
            blength[i] = LA.norm(pos_atom_cart - pos_cart)
        return blength
            
    def disorder_error(self,atomNo):
        dr = np.sqrt(self.basis[atomNo][4]/(8*np.pi**2))
        dz = np.sqrt(self.basis[atomNo][5]/(8*np.pi**2))
        errdr = 0.5*dr*(self.errors[atomNo][4]/self.basis[atomNo][4])
        errdz = 0.5*dz*(self.errors[atomNo][5]/self.basis[atomNo][5])
        return dr, dz, errdr, errdz
        
    
    def __repr__(self):
        repr_super = super(UnitCell,self).__repr__()
        return "%s\n%i atoms in unit cell" % (repr_super,len(self.names))
        
    def atomToStr(self,no,showErrors=True):
        param = self.basis[no][1:]
        name = self.names[no]
        if (self.errors is not None) and showErrors:
            err = self.errors[no][1:]
            l = []
            for t in zip(param,err):
                [l.append(ti) for ti in t]
            return "%s  (%.5f +- %.5f)  (%.5f +- %.5f)  (%.5f +- %.5f)  (%.4f +- %.4f)  (%.4f +- %.4f)  (%.4f +- %.4f)" % (name,*l)
        else:
            return "%s     %.5f     %.5f     %.5f  %.4f  %.4f  %.4f" % (name,*param)
        
    def domainsToStr(self):
        s = ""
        for mat, occu in zip(self.coherentDomainMatrix,self.coherentDomainOccupancy):
            matT = np.vstack((mat[:,:-1],mat[:,-1]))
            s += ("Coherent %.5f   " % occu) + np.array2string(np.array(matT).flatten(), formatter={'float_kind':lambda x: "%.5f" % x},max_line_width=100000)[1:-1] + "\n"
        return s
            
    def latticeRODStr(self):
        a,alpha,_,_ = self.getLatticeParameters()
        return "%.4f %.4f %.4f %.4f %.4f %.4f" % (*a,*np.rad2deg(alpha))
    
    
    def parameterStr(self):
        st = ""
        for i in range(len(self.names)):
            st += str(i).zfill(2) + "  " + self.atomToStr(i) + "\n"
        return st
    
    def parameterStrRod(self):
        st = ""
        for i in range(len(self.names)):
            st += self.atomToStr(i,False) + "\n"
        return st
            
    def __str__(self):
        st = repr(self) + "\n"
        st += "id  " + UnitCell.parameterOrder + "\n"
        return st + self.parameterStr()
            
    def writeSURfile(self,filename):
        with open(filename,'w') as f:
            st = "return\n" + self.domainsToStr() + self.latticeRODStr() + "\n" + self.parameterStrRod()
            f.write(st)
            
    def toRODStr(self):
        return "return\n" + self.domainsToStr() + self.latticeRODStr() + "\n" + self.parameterStrRod()
    
    def toStr(self):
        return "return\n" + self.domainsToStr() + self.latticeRODStr() + "\n" + UnitCell.parameterOrder + "\n" + self.parameterStr()
    
    def toXYZfile(self,xyzfile, ucx=1,ucy=1,ucz=1,translate=np.array([0.,0.,0.])):
        xyz_array = self.pos_cart_all(ucx,ucy,ucz,translate)
        with open(xyzfile,'w') as f:
            noatoms = len(xyz_array)
            f.write("{}\n".format(noatoms))
            f.write("{}\n".format(self.latticeRODStr()))
            np.savetxt(f,xyz_array,fmt=['%s','%.6f','%.6f','%.6f'])
            
    @staticmethod
    def fromFile(filename):
        f, ext = os.path.splitext(filename)
        if ext == '':
            files = glob.glob(filename + ".*")
            if files:
                extensions = [os.path.splitext(fi)[1] for fi in files]
                for fext in extensions:
                    if fext in ['.sur', '.bul']:
                        ext = os.path.splitext(filename)[1]
                        break
                else:
                    if len(extensions) == 1:
                        ext = extensions[0]
                    else:
                        ext = extensions
            else:
                if os.path.isfile(filename):
                    raise IOError("File %s has no file extension." % filename)
                else:
                    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
        elif ext not in ['.sur', '.bul']:
            if not os.path.isfile(filename):
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), filename)
                
        if ext == '.sur':
            return UnitCell.fromSURfile(f + '.sur')
        elif ext == '.bul':
            return UnitCell.fromBULfile(f + '.bul')
        else:
            try:
                import ase.io
            except ImportError:
                raise IOError("File %s has no known file extension. Try to install ASE for more supported file types." % filename)
            if type(ext) == list:
                for e in ext:
                    try:
                        atoms = ase.io.read(f + e)
                        break
                    except:
                        pass
                else:
                    raise IOError("Cannot read file %s." % filename)
            else:
                atoms = ase.io.read(f + ext)
            if not atoms.cell:
                raise IOError("File %s contains no valid crystal lattice parameters" % filename)
            uc = UnitCell(atoms.cell.lengths(), atoms.cell.angles())
            coord = atoms.get_scaled_positions()
            symb = atoms.get_chemical_symbols()
            for sym, xyz in zip(symb, coord):
                uc.addAtom(sym, xyz, 0., 0., 1.)
            return uc
    
    @classmethod
    def fromSURfile(cls, surfile):
        with open(surfile,'r') as f:
            #next(f) # return line
            line = next_skip_comment(f).split()
            domainoccu = []
            domainmatrix = []
            while(line[0] == 'Coherent'):
                domainoccu.append(float(line[1]))
                mat = np.array(line[2:],dtype=np.float64)
                if mat.size > 10:
                    domainmatrix.append(np.ascontiguousarray(np.vstack( (mat[:-3].reshape((3,3)).T,mat[-3:].T)).T.astype(np.float64)))
                else:
                    domainmatrix.append(np.ascontiguousarray(np.vstack((mat.reshape((3,3)).T, np.array([0,0,0]))).T.astype(np.float64)))
                line = next_skip_comment(f).split()
            if not domainoccu:
                domainoccu.append(1.)
                domainmatrix.append(np.ascontiguousarray(np.vstack((np.identity(3),np.array([0,0,0]))).T.astype(np.float64)))
            latticeparams = np.array(line,dtype=np.float64)
            names = []
            basis = []
            for l in f:
                line = l.split()
                if line:
                    names.append(line[0])
                    basis.append(np.concatenate(([atomic_number(line[0])], np.array(line[1:],dtype=np.float64) )))
        uc = cls(latticeparams[:3],latticeparams[3:])
        uc.coherentDomainMatrix = domainmatrix
        uc.coherentDomainOccupancy = domainoccu
        if len(basis) >= 2:
            basis = np.vstack(basis)
            if basis.shape[1] == 7:
                uc.basis = basis
                uc.basis_0 = np.copy(uc.basis)
                uc.names = names
                uc.dw_increase_constraint = np.ones(uc.basis.shape[0],dtype=np.bool_)
                uc._test_special_formfactors()
                return uc
            else:
                raise ValueError("wrong number of atomic parameters. read basis is: {}".format(basis))
        elif basis[0].size == 7:
            uc.basis = np.array([basis])
            uc.dw_increase_constraint = np.ones(basis.shape[0],dtype=np.bool_)
            uc.names = names
            uc.basis_0 = np.copy(uc.basis)
            uc._test_special_formfactors()
            return uc
        else:
            raise ValueError("wrong number of atomic parameters. read basis is: {}".format(basis))
    
    @classmethod
    def fromStr(cls, string):
        xprfile = False
        with util.StringIO(string) as f:
            line = next_skip_comment(f).split()
            domainoccu = []
            domainmatrix = []
            while(line[0] == 'Coherent'):
                domainoccu.append(float(line[1]))
                mat = np.array(line[2:],dtype=np.float64)
                if mat.size > 10:
                    domainmatrix.append(np.ascontiguousarray(np.vstack( (mat[:-3].reshape((3,3)).T,mat[-3:].T)).T.astype(np.float64)))
                else:
                    domainmatrix.append(np.ascontiguousarray(np.vstack((mat.reshape((3,3)).T, np.array([0,0,0]))).T.astype(np.float64)))
                line = next_skip_comment(f).split()
            if not domainoccu:
                domainoccu.append(1.)
                domainmatrix.append(np.ascontiguousarray(np.vstack((np.identity(3),np.array([0,0,0]))).T.astype(np.float64)))
            latticeparams = np.array(line,dtype=np.float64)
            names = []
            basis = []
            errors = []
            indices = []
            statistics = dict()
            for l in f:
                line = l.rsplit('//')[0]
                line = line.split()
                if line:
                    if 'Name' in line or '=' in line: # parameter or statistics line
                        if '=' in line:
                            try:
                                splitted = [n.split(',') for n in l.split('=')]
                                splitted = [item for sublist in splitted for item in sublist]
                                for i in range(0,len(splitted),2):
                                    statistics[splitted[i].strip()] = float(splitted[i+1])
                            except Exception:
                                print("Cannot read statistics string: %s" % l)
                        continue
                    if '+-' in line:
                        xprfile = True
                        indices.append(int(line[0]))
                        names.append(line[1])
                        params = re.findall(r'\(([^)]+)',l)
                        params_array = np.array([np.array(p.split('+-'),dtype=np.float64) for p in params]).T
                        basis.append(np.concatenate(([int(line[0])], params_array[0] )))
                        errors.append(np.concatenate(([int(line[0])], params_array[1] )))
                    else:
                        if line[0].isnumeric():
                            
                            names.append(line[1])
                            basis.append(np.concatenate(([int(line[0])], np.array(line[2:],dtype=np.float64) )))
                            indices.append(int(line[0]))
                        else:
                            names.append(line[0])
                            indices.append(np.nan)
                            basis.append(np.concatenate(([np.nan], np.array(line[1:],dtype=np.float64) )))
            
            basis_save = np.vstack(basis).astype(np.float64)
            basis = np.empty_like(basis_save)
            
            indices = np.array(indices)
            names_save = np.array(names)
            names = np.copy(names_save)
            
            mask_read = ~np.isnan(indices) # atoms with explicitly given index
            indices = indices.astype(np.intp) # convert to int to allow indexing.
            indices_explicit = indices[mask_read] # only the entries, where an index is given
            
            mask_write = np.zeros_like(mask_read, dtype=np.bool_)
            mask_write[indices_explicit] = True # mask where to store the atoms (if array is sorted!)
            sortidx = np.argsort(indices_explicit) # get order of atoms with explicit indices
            
            if indices_explicit.size != np.unique(indices_explicit).size:
                raise ValueError("Atom indices are not unique. This may be caused by assigning two atoms the same index.")
            
            basis[mask_write] = basis_save[mask_read][sortidx] # store atoms with explicit index
            names[mask_write] = names_save[mask_read][sortidx]
            
            basis[~mask_write] = basis_save[~mask_read] # all other just retain the order of the file
            names[~mask_write] = names_save[~mask_read]
            
            names = list(names)

            if xprfile:
                errors_save = np.vstack(errors).astype(np.float64)
                errors = np.copy(errors_save)
                errors[mask_write] = errors_save[mask_read][sortidx]
                errors[~mask_write] = errors_save[~mask_read]
                
        uc = cls(latticeparams[:3],latticeparams[3:])
        uc.coherentDomainMatrix = domainmatrix
        uc.coherentDomainOccupancy = domainoccu
        uc.statistics = statistics
        uc.dw_increase_constraint = np.ones(uc.basis.shape[0],dtype=np.bool_)
        if len(basis.shape) >= 2:
            basis = np.vstack(basis)
            if basis.shape[1] == 7:
                uc.basis = basis
                if xprfile:
                    uc.errors = errors
                uc.names = names
                uc.basis_0 = np.copy(uc.basis)
                uc.dw_increase_constraint = np.ones(uc.basis.shape[0],dtype=np.bool_)
                uc._test_special_formfactors()
                return uc
            else:
                raise ValueError("wrong number of atomic parameters. read basis is: {}".format(basis))
        elif basis[0].size == 7:
            uc.basis = np.array([basis])
            if xprfile:
                uc.errors = np.array([errors])
            uc.names = names
            uc.basis_0 = np.copy(uc.basis)
            uc.dw_increase_constraint = np.ones(uc.basis.shape[0],dtype=np.bool_)
            uc._test_special_formfactors()
            return uc
        else:
            raise ValueError("wrong number of atomic parameters. read basis is: {}".format(basis))

    
    #ugly!!!! please redo this for faster file reads!
    @classmethod
    def fromBULfile(cls, bulfile,DW=0.):
        with open(bulfile,'r') as f:
            #next_skip_comment(f) # return line
            line = next_skip_comment(f).split()
            domainoccu = []
            domainmatrix = []
            while(line[0] == 'Coherent'):
                domainoccu.append(float(line[1]))
                mat = np.matrix(line[2:],dtype=np.float64)
                if mat.size > 10:
                    domainmatrix.append(np.ascontiguousarray(np.vstack( (mat[:-3].reshape((3,3)).T,mat[-3:].T)).T.astype(np.float64)))
                else:
                    domainmatrix.append(np.ascontiguousarray(np.vstack((mat.reshape((3,3)).T, np.array([0,0,0]))).T.astype(np.float64)))
                line = next_skip_comment(f).split()
            if not domainoccu:
                domainoccu.append(1.)
                domainmatrix.append(np.ascontiguousarray(np.vstack((np.identity(3),np.array([0,0,0]))).T.astype(np.float64)))
            latticeparams = np.array(line,dtype=np.float64)
            #basis = np.loadtxt(f,converters={0 : lambda x :  atomic_number(x.decode("utf-8")) })
            #f.seek(basisstartline)
            names = []
            basis = []
            for l in f:
                line = l.split()
                if line:
                    names.append(line[0])
                    basis.append(np.concatenate(([atomic_number(line[0])], np.array(line[1:],dtype=np.float64) )))
            basis = np.vstack(basis).astype(np.float64)
        uc = cls(latticeparams[:3],latticeparams[3:])
        uc.coherentDomainMatrix = domainmatrix
        uc.coherentDomainOccupancy = domainoccu
        if len(basis.shape) >= 2:
            basis = np.vstack(basis)
            if basis.shape[1] == 7:
                for i,b in enumerate(basis):
                    uc.addAtom(names[i],b[1:4],DW,DW,1)
                uc._test_special_formfactors()
                return uc
            else:
                raise ValueError("wrong number of atomic parameters. read basis is: {}".format(basis))
        elif basis[0].size == 7:
            uc.basis = np.array([basis]).astype(np.float64)
            for i,b in enumerate(basis):
                uc.addAtom(names[i],b[1:4],DW,DW,1)
            uc._test_special_formfactors()
            return uc
        else:
            raise ValueError("wrong number of atomic parameters. read basis is: {}".format(basis))
   


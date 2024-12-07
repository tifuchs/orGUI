import numpy as np
import numpy.linalg as LA
from silx.io import dictdump
import xraydb
from .. import util
import warnings
import random
import json
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
from abc import ABC

from .element_data import cov_radii_array, rgb_array


special_elementcolors = {
'Pt2+' : (1.0, 0.75, 0.),
'Pt4+' : (0., 0., 1.)
}

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


class ParameterType(enum.IntFlag):
    UNDEFINED = enum.auto()
    ABSOLUTE = enum.auto()
    RELATIVE = enum.auto()
    

@dataclass
class Parameter:
    name: str
    indices: np.ndarray
    kind: ParameterType
    limits: tuple = (-np.inf, np.inf)
    prior: typing.Any = None
    settings: dict = field(default_factory=dict)
    factors: np.ndarray = None
    #children: list = field(default_factory=list)
    value: float = None
    error: float = None
    
    def __post_init__(self): # ensures correct data types
        for field in dataclasses.fields(self):
            value = getattr(self, field.name)
            if field.name == 'indices':
                idx = np.array(value, dtype=np.intp)
                if len(idx.shape) == 2:
                    idx = tuple(idx)
                setattr(self, field.name, idx)
            elif field.name == 'factors':
                if value is not None:
                    factors = np.array(value, dtype=np.float64)
                    setattr(self, field.name, factors)
            elif field.name in ['prior']:
                pass # ToDo convert to distribution here!
            else:
                if value is not None:
                    setattr(self, field.name, field.type(value))
    
    def asdict(self):
        d = dataclasses.asdict(self)
        prior_a = np.array(d['prior'])
        if not np.issubdtype(prior_a.dtype, np.number):
            d['prior'] = None        
        #d['children'] = [c.asdict() for c in d['children']]
        for p in ['indices_t', 'factors_t']:
            d.pop(p, None)
        return d
        
class LinearFitFunctions(ABC):
    parameterOrder = ""
    parameterLookup = {}
    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))
    
    def __init__(self):
        self.basis = np.array([])
        self.parameters = {
            'absolute' : [],
            'relative' : []
        }
        self.basis_0 = np.array([])
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None
        
    def parametersToDict(self):
        d = dict()
        d['basis_0'] = self.basis_0
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
        for p in self.parameters:
            for dkey in sorted(d[p].keys()):
                self.parameters[p].append(Parameter(**d[p][dkey]))

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn("Can not apply fit parameter values to this %s: %s" % (self.__class__.__name__,e)) 
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
        
    def addFitParameter(self,indexarray,limits=(-np.inf,np.inf),**keyargs):
        par = indexarray
        try:
            parindexes = np.array([p if isinstance(p, (np.integer,int)) else self.parameterLookup[p] for p in np.atleast_1d(par)],dtype=np.intp)
        except Exception as e:
            raise ValueError("Invalid parameter name %s., Parameters for class %s are %s" % (par, self.__class__.__name__, list(self.parameterLookup.keys()))) from e

        if 'name' in keyargs:
            name = keyargs['name']
        else:
            parameternames = tuple([self.parameterLookup_inv[n] for n in parindexes])
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
            parindexes = np.array([p if isinstance(p, (np.integer,int)) else self.parameterLookup[p] for p in np.atleast_1d(par)],dtype=np.intp)
        except Exception as e:
            raise ValueError("Invalid parameter name %s., Parameters for class %s are %s" % (par, self.__class__.__name__, list(self.parameterLookup.keys()))) from e
        
        if 'name' in keyargs:
            name = keyargs['name']
        else:
            parameternames = tuple([self.parameterLookup_inv[n] + '_r' for n in parindexes])
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
            parameternames = tuple([self.parameterLookup_inv[n] for n in par.indices])
            return "{:10} {:15} {:20} {:10} {:10}".format(par.name,str(parameternames),val,str(par.limits))
        elif par.kind & ParameterType.RELATIVE:
            val = np.mean( (self.basis[par.indices] - basis_0[par.indices]) / par.factors)
            parameternames = tuple([self.parameterLookup_inv[n] + '_r' for n in par.indices])


            return "{:10} {:20} {:10} {:10} {:10}".format(par.name,str(parameternames),val, str(par.factors) ,str(par.limits))
        else:
            raise ValueError("Unvalid parameter type %s for class %s" % (par.kind, self.__class__.__name__))

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
        
    def parameter_list(self):
        return self.parameters['absolute'] + self.parameters['relative']
        
    @property
    def fitparnames(self):
        #pars = [self.fitparToStr(i,True) for i in range(len(self.fitparameters))] + [self.fitparToStr(i + len(self.fitparameters),True) for i in range(len(self.relfitparam))]
        pars = []
        for p_kind in self.parameters:
            pars += [p.name for p in self.parameters[p_kind]]
        return pars
    
    @property
    def priors(self):
        priorlist = []
        pars = []
        for p_kind in self.parameters:
            pars += [p for p in self.parameters[p_kind]]
        for par in pars:
            if par.prior is None:
                if tuple(par.limits) == (-np.inf,np.inf):
                    raise Exception(f"Infinite range given as prior for parameter {par.name}. You must provide an explicit prior.")
                else:
                    priorlist.append(tuple(par.limits)) #has to be converted to correct distribution in CTROptimizer
            else:
                priorlist.append(par.prior) # real prior distribution
        return priorlist
        
        
def _ensure_contiguous(*arrays, testOnly=False, astype=None):
    if testOnly:
        for a in arrays:
            if not a.flags['C_CONTIGUOUS']:
                raise ValueError("ndArray is not contiguous. Use np.ascontiguousarray() to convert the array")
            if astype is not None and not a.dtype == astype:
                raise ValueError("Wrong data type of ndArray. Should be %s, but is %s" % (astype,a.dtype))
        return arrays
    else:
        a_c = []
        for a in arrays:
            if not a.flags['C_CONTIGUOUS']:
                warnings.warn("h,k,l is not contiguous. This causes calculation speed loss. Use np.ascontiguousarray() to convert the arrays")
                a_c.append(np.ascontiguousarray(a, dtype=astype))
            elif astype is not None and not a.dtype == astype:
                warnings.warn("Not matching data type of ndArray. Should be %s, but is %s. This causes calculation speed loss." % (astype,a.dtype) )
                a_c.append(np.ascontiguousarray(a, dtype=astype))
            else:
                a_c.append(a)
        return a_c
    
            
   
def next_skip_comment(it, comment=('//','return')):
    while(True):
        line = next(it)
        if not line.startswith(comment):
            line = line.rsplit(comment[0])[0].strip()
            if line:
                break
    return line
        
    
SQRT2pi = np.sqrt(2*np.pi)
def _eDensity_old(z ,element, deltaZ, Qpara, deltaPara, E):
    atom_parameters = xraydb.read_Waasmaier(element)
    atom_parameters['exponents'] /= (4*np.pi)**2
    dispersion = xraydb.f1_chantler(element,E) + 1j*xraydb.f2_chantler(element,E)
    
    rho_c = ((atom_parameters['c'] + dispersion) / (SQRT2pi*deltaZ)) * \
            np.exp( -0.5 *( (deltaPara**2) * (Qpara**2) +  (z**2)/(deltaZ**2) ))
            
    rho_0 = np.zeros_like(z)
    for a,b in zip(atom_parameters['scale'],atom_parameters['exponents']):
        exp_dpara = b + 0.5*(deltaPara**2)
        exp_dz = b + 0.5*(deltaZ**2)
        rho_0 += (a/(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- exp_dpara * (Qpara**2)  - ((z**2)/(4*exp_dz))  )
        #rho_0 += (1./(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- ((z**2)/(4*exp_dz))  )
    
    return rho_c + rho_0

def DWtoDisorder(dw):
    return np.sqrt(dw/(8*np.pi**2))
        
        
# returns a,b,c for Q = 4pi/lambda * sin th (instead of s)
def readWaasmaier(element):
    xraydb_t = xraydb.get_xraydb()
    wtab = xraydb_t.tables['Waasmaier']

    row = xraydb_t.query(wtab)
    if isinstance(element, int):
        row = row.filter(wtab.c.atomic_number==element).one()
    else:
        row = row.filter(wtab.c.ion==element.title()).one()
    #if len(row) > 0:
    #    row = row[0]

    c = row.offset
    a = json.loads(row.scale); b = np.array(json.loads(row.exponents)) / (4*np.pi)**2;
    xraydb_t.close()
    return np.concatenate((a,b,[c]))
        
# incorrect for ions!!
def readDispersion(element,E):
    return xraydb.f1_chantler(atomic_number(element),E), xraydb.f2_chantler(atomic_number(element),E)
        
def atomic_number(elementname):
    if elementname in special_numbers:
        return special_numbers[elementname]
    elif elementname[-1] == '+' or elementname[-1] == '-':
        return xraydb.atomic_number(elementname[:-2])
    else:
        return xraydb.atomic_number(elementname)
    
def estimateDispersionCompound(name,E):
    f1 = 0.; f2 = 0.
    for atom, frac in xraydb.chemparse(name).items():
        f1 += xraydb.f1_chantler(atom, E)*frac
        f2 += xraydb.f2_chantler(atom, E)*frac
    return f1, f2
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
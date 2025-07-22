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
import copy
import warnings
from functools import cached_property
#from functools import partial
from .. import util
from .CTRplotutil import ctrfigure, CTR, CTRCollection

from .CTRcalc import UnitCell, SXRDCrystal

from scipy import stats

from typing import Callable

class CTROptimizer():
    
    def __init__(self,xtal,CTRs):
        self.CTRs = copy.deepcopy(CTRs)
        self.CTRs.sort(key=lambda x: abs(x.hk[0]) + abs(x.hk[1]))
        self.xtal = copy.deepcopy(xtal)
        self.scaling = util.get_scale_chi2
        
    def prepareFit(self):
        self.startp, self.lower_bounds, self.higher_bounds = self.xtal.getStartParamAndLimits()
        self.bounds = (self.lower_bounds, self.higher_bounds)
        for ctr in self.CTRs:
            ctr.invrelerrsqrd_weight = ctr.weight * ctr.err**-2
        
    def get_bounds(self):
        return self.bounds
    
    def get_parameters(self):
        return self.xtal.getInitialParameters()
    
    def weighted_residues2(self,x=None):
        if x is not None:    
            self.xtal.setParameters(x)
        residues = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
            residues.append( (ctr.invrelerrsqrd_weight/scale**2) * ((ctr.sfI*scale - F_theo )**2))
        return np.concatenate(residues)
    
    def residues(self,x=None):
        if x is not None:    
            self.xtal.setParameters(x)
        residues = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
            residues.append(ctr.sfI*scale - F_theo )
        return np.concatenate(residues)
    
    def flat_data(self,specular=True):
        dat = []
        err = []
        for i,ctr in enumerate(filter(lambda x: specular or x.hk != (0,0),self.CTRs)):
            dat.append(ctr.sfI)
            err.append(ctr.err)
        return np.concatenate(dat), np.concatenate(err)
    
    @property    
    def nopoints(self):
        F,err = self.flat_data()
        return F.size
    
    def flat_Fcalc(self, x=None):
        if x is not None:    
            self.xtal.setParameters(x)
        F = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
            F.append(F_theo/scale)
        return np.concatenate(F)
            
    def Rfactor(self,x=None):
        if x is not None:    
            self.xtal.setParameters(x)
        residues = []
        Fobs = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
            residues.append(np.abs(ctr.sfI*scale - F_theo) )
            Fobs.append(np.abs(ctr.sfI*scale))
        residues = np.sum(np.concatenate(residues))
        Fobs = np.sum(np.concatenate(Fobs))
        return residues/Fobs
        
    
    def weighted_residues(self,x=None):
        if x is not None:    
            self.xtal.setParameters(x)
        residues = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
            residues.append( np.sqrt(ctr.weight) * ( (ctr.sfI*scale - F_theo ) / (ctr.err*scale) ))
        return np.concatenate(residues)
    
    def fitness(self,x):
        self.xtal.setParameters(x)
        sumchi2 = 0.
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
            sumchi2 += np.sum( (ctr.invrelerrsqrd_weight/scale**2) * ((ctr.sfI*scale - F_theo )**2) )
        return [sumchi2]
    
    def statistics(self,x):
        self.xtal.setParameters(x)
        residues2 = self.weighted_residues2()
        
        Rfactor = self.Rfactor()
        
        stat = dict()
        
        #variance = np.concatenate([ctr.err**2 for ctr in self.CTRs])
        #varmat_i = np.diag(1/variance)
        
        chi2_result = np.sum(residues2)
        pvalue = 1 - stats.chi2.cdf(chi2_result,residues2.size - x.size)
        chi2_red = chi2_result/(residues2.size - x.size)
        
        pcov = util.leastsq_covariance(self.weighted_residues,x)

        self.errors = np.sqrt(np.diag(pcov)*chi2_red)
        self.xtal.setFitErrors(self.errors)
        self.xtal.setParameters(x)
        
        stat['Chisqr'] = chi2_result
        stat['nodatapoints'] = residues2.size
        stat['Chisqr_red'] = chi2_red
        stat['noparameters'] = x.size
        stat['pvalue'] = pvalue
        stat['Rfactor'] = Rfactor
        stat['covariance'] = pcov
        
        return stat
        
    
    def evaluateStatistics(self,x):
        warnings.warn("usage of evaluateStatistics is deprecated, use CTROptimizer.statistics instead!",DeprecationWarning)
        
        self.xtal.setParameters(x)
        residues2 = self.weighted_residues2()
        
        #variance = np.concatenate([ctr.err**2 for ctr in self.CTRs])
        #varmat_i = np.diag(1/variance)
        
        chi2_result = np.sum(residues2)
    
        pvalue = 1 - stats.chi2.cdf(chi2_result,residues2.size - x.size)
        
        chi2_red = chi2_result/(residues2.size - x.size)
        
        pcov = util.leastsq_covariance(self.residues,x)

        errors = np.sqrt(np.diag(pcov)*chi2_red)
        self.xtal.setFitErrors(errors)
        self.xtal.setParameters(x)
        
        return chi2_result, chi2_red, pvalue, residues2.size
    
    def printStatistics(self,x):
        #chi2_result, chi2_red , pvalue, nodatapoints = self.evaluateStatistics(x)
        stat = self.statistics(x)
        print("Chisqr = %.4f, Chisqr_red = %.4f, R-factor = %.4f ,p-value = %.6f, n_refl = %s" % (stat['Chisqr'],stat['Chisqr_red'],stat['Rfactor'],stat['pvalue'],stat['noparameters']))
        
    
    def get_name(self):
        return "CTR optimizer"
    
    def setCTRPlotSettings(self,lrange,plotsize,**settings):
        self.lrange = lrange
        self.plotsize = plotsize
        self.settings = settings
    
class FitCallback(object):
    global_counter = 1
    
    def __init__(self, function: Callable, 
                       bounds_low: list,
                       bounds_high: list,
                       init: list, **kwargs):
        """Wrapper for a fit callback

        function : [[SXRDCrystal, list[float]], None]
        The function takes the SXRDCrystal and an array of parameters as arguments
        
        """
        self.name = kwargs.get('name', "default-%s" % FitCallback.global_counter)
        FitCallback.global_counter += 1
        
        self.inital = np.asarray(init)
        self.n_pars = self.inital.size
        
        self.current_values = np.copy(self.inital)
        
        up_bnds = np.asarray(bounds_low)
        if up_bnds.size != self.n_pars:
            raise ValueError("Number of upper bounds does not match number of initial parameters.")
        low_bnds = np.asarray(bounds_high)
        if low_bnds.size != self.n_pars:
            raise ValueError("Number of lower bounds does not match number of initial parameters.")
            
        self.bounds = (low_bnds, up_bnds)
        
        self.function = function
        
    def __call__(self, xtal: SXRDCrystal, x: list) -> None:
        self.current_values = np.copy(x)
        return self.function(xtal, x)
        
    def get_parameters(self, xtal: SXRDCrystal = None):
        return self.current_values
        
    def set_parameters(self, xtal: SXRDCrystal, x: list):
        self.current_values = np.copy(x)
        return self.function(xtal, x)
        
    def set_errors(self, xerror):
        self.errors = xerror
        
    def get_bounds(self):
        return self.bounds
        
    def __repr__(self):
        return "<%s : %s (%s pars)>" % (type(self).__name__, self.name, self.n_pars)
    

class CTROptAngleCorrection(CTROptimizer):
    
    def __init__(self,*args,**kwargs):
        super().__init__(*args, **kwargs)
        self.scaleindividual = False
        self.useAnglecorr = False
        self.dw_zconstraints = False
        self.phasevelocity = 1.
        self.nic = 0
        self.callbacks = []
    
    def register_fit_callback(self, function: Callable, 
                       bounds_low: list,
                       bounds_high: list,
                       init: list,
                       **kwargs):
        callback = FitCallback(function, bounds_low, bounds_high, init, **kwargs)
        self.callbacks.insert(0,callback)
        return callback.name
    
    @property
    def callback_names(self):
        names = [n.name for n in self.callbacks]
        return names
        
    def unregister_fit_callback(self, name: str):
        try:
            idx = self.callback_names.index(name)
        except ValueError as e:
            raise ValueError("%s is not a registered callback") from e
        del self.callbacks[idx]
    
    def prepareFit(self, phaselim=[0,2*np.pi],amplim=[0,.75], start=[0.,0.]):
        
        self.startp, self.lower_bounds, self.higher_bounds = self.xtal.getStartParamAndLimits()
        if self.useAnglecorr:
            self.bounds = (np.concatenate(([phaselim[0],amplim[0]],self.lower_bounds)), np.concatenate(([phaselim[1],amplim[1]],self.higher_bounds)))
        else:
            self.bounds = (self.lower_bounds, self.higher_bounds)
        for ctr in self.CTRs:
            ctr.invrelerrsqrd_weight = np.sqrt(ctr.weight) / ctr.err
        self.phase, self.amp = start
        
        for cb in reversed(self.callbacks):
            self.bounds = (np.concatenate((cb.bounds[0],self.bounds[0])), np.concatenate((cb.bounds[1],self.bounds[1])))
        
        if self.dw_zconstraints:
            constr = self.get_inequalconstraints()
            self.nic = constr.size
        self.fitparnames = self.xtal.fitparnames
        self.priors = self.xtal.priors

            
    def get_inequalconstraints(self):
        constraints = []
        if self.dw_zconstraints:
            dwc_enable = self.xtal.getSurfaceDWConstraintEnable()
            sur_basis = self.xtal.getSurfaceBasis()[dwc_enable]
            sur_basis = sur_basis[np.argsort(sur_basis[:,3])] # order in z direction
            #iDW_bulk = self.xtal.uc_bulk.basis[0,4]
            #oDW_bulk = self.xtal.uc_bulk.basis[0,5]
            iDWconstr = np.diff(self.xtal.uc_bulk.basis[0,4] - sur_basis[:,4], prepend=0)
            oDWconstr = np.diff(self.xtal.uc_bulk.basis[0,5] - sur_basis[:,5], prepend=0)
            dw_constraints = np.concatenate((iDWconstr,oDWconstr))
            constraints.append(dw_constraints)
        return np.concatenate(constraints)
        
    def applyCorrections(self):
        F_obs = []
        F_t = []
        F_err = []
        if self.useAnglecorr:
            for i,ctr in enumerate(self.CTRs):
                F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
                if hasattr(ctr,'angles'):
                    anglecorr = self.get_anglecorrection(ctr.angles['omega'])
                else:
                    anglecorr = 1.
                ctr *= anglecorr
                if self.scaleindividual or ctr.hk == (0,0):
                    scale = self.scaling(F_theo, ctr.sfI, ctr.err) # scale CTR
                    ctr *= scale
                else:
                    F_obs.append(ctr.sfI)
                    F_t.append(F_theo)
                    F_err.append(ctr.err)
            if not self.scaleindividual:
                scale = self.scaling(np.concatenate(F_t), np.concatenate(F_obs), np.concatenate(F_err))
                for i,ctr in enumerate(filter(lambda x: x.hk != (0,0),self.CTRs)):
                    ctr *= scale
                    ctr.invrelerrsqrd_weight = np.sqrt(ctr.weight) / ctr.err
            
            self.amp = 0.
        else:
            warnings.warn("Angle correction was not enabled. Skip applyCorrections.")
    
    def get_nic(self):
        return self.nic

    def fitness(self,x):
        if self.dw_zconstraints:
            return np.concatenate(([np.sum(self.weighted_residues2(x))], self.get_inequalconstraints()))
        else:
            return [np.sum(self.weighted_residues2(x))]
            
    def log_prob(self, x):
        resid, err = self.weighted_residues_errors(x)
        return -0.5 * np.sum( resid**2 + np.log(2*np.pi*err**2) )
    
    def get_parameters(self):
        if self.useAnglecorr:
            pars = np.concatenate(([self.phase, self.amp],self.xtal.getInitialParameters()))
        else:
            pars = self.xtal.getInitialParameters()
        
        for cb in reversed(self.callbacks):
            pars = np.concatenate((cb.get_parameters(self.xtal), pars))
        return pars
        
    def set_parameters(self, x):
        counter = 0
        for cb in self.callbacks:
            cb.set_parameters(self.xtal, x[counter:counter+cb.n_pars])
            counter += cb.n_pars

        x = x[counter:]
        
        if self.useAnglecorr:
            self.phase, self.amp = x[:2]
            self.xtal.setParameters(x[2:])
        else:
            self.xtal.setParameters(x)
            
    def set_errors(self,xerror):
        self.errors = xerror
        
        counter = 0
        for cb in self.callbacks:
            cb.set_errors(self.xtal, xerror[counter:counter+cb.n_pars])
            counter += cb.n_pars
        xerror = xerror[counter:]
        
        if self.useAnglecorr:
            self.xtal.setFitErrors(xerror[2:])
        else:
            self.xtal.setFitErrors(xerror)
        
    def get_anglecorrection_(self, omega, x=None): # old, unused
        if x is not None:
            self.phase, self.amp = x[:2]
        return 1 + self.amp*np.sin(self.phasevelocity*(omega + self.phase))
        
    def get_anglecorrection(self, omega, x=None): # improved the fit a bit, but not significantly...
        if self.useAnglecorr:
            if x is not None:
                self.phase, self.amp = x[:2]
            return np.sqrt(1 + self.amp*np.sin(self.phasevelocity*(omega + self.phase)))
        else:
            return 1.
        
    def weighted_residues2(self,x):
        return self.weighted_residues(x)**2
    
    def residues(self,x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            if hasattr(ctr,'angles'):
                anglecorr = self.get_anglecorrection(ctr.angles['omega'])
            else:
                anglecorr = 1.
            F_obs_corr = ctr.sfI*anglecorr
            F_err_corr = ctr.err*anglecorr
            if self.scaleindividual or ctr.hk == (0,0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr) # scale CTR
                residues.append(F_obs_corr*scale - F_theo  )
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            return np.concatenate(residues)
        else:
            scale = self.scaling(np.concatenate(F_t), np.concatenate(F_obs), np.concatenate(F_err))
            for i,ctr in enumerate(filter(lambda x: x.hk != (0,0),self.CTRs)):
                residues.append( F_obs[i]*scale - F_t[i] )
            return np.concatenate(residues)
        return np.concatenate(residues)
        
        
    def Rfactor(self,x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            if hasattr(ctr,'angles'):
                anglecorr = self.get_anglecorrection(ctr.angles['omega'])
            else:
                anglecorr = 1.
            F_obs_corr = ctr.sfI*anglecorr
            F_err_corr = ctr.err*anglecorr
            if self.scaleindividual or ctr.hk == (0,0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr) # scale CTR
                residues.append(F_obs_corr*scale - F_theo  )
                F_obs.append(F_obs_corr*scale)
                F_err.append(F_err_corr)
                F_t.append(F_theo)
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            residues = np.concatenate(residues)
            F_obs = np.concatenate(F_obs)
        else:
            scale = self.scaling(np.concatenate(F_t), np.concatenate(F_obs), np.concatenate(F_err))
            for i,ctr in enumerate(filter(lambda x: x.hk != (0,0),self.CTRs)):
                residues.append( F_obs[i]*scale - F_t[i] )
                F_obs[i] *= scale
            residues = np.concatenate(residues)
            F_obs = np.concatenate(F_obs)
        residues = np.sum(np.abs(residues))
        return residues/np.sum(np.asarray(F_obs))
        
    
    def weighted_residues(self,x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            if hasattr(ctr,'angles'):
                anglecorr = self.get_anglecorrection(ctr.angles['omega'])
            else:
                anglecorr = 1.
            F_obs_corr = ctr.sfI*anglecorr
            F_err_corr = ctr.err*anglecorr
            if self.scaleindividual or ctr.hk == (0,0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr) # scale CTR
                residues.append((ctr.weight/(F_err_corr*scale)) * (F_obs_corr*scale - F_theo ) )
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            return np.concatenate(residues)
        else:
            _, err = self.flat_data(False)
            scale = self.scaling(np.concatenate(F_t), np.concatenate(F_obs), np.asarray(err))
            for i,ctr in enumerate(filter(lambda x: x.hk != (0,0),self.CTRs)):
                residues.append( (ctr.weight/(scale*F_err[i])) * (F_obs[i]*scale - F_t[i] ) )
            return np.concatenate(residues)
            
    def weighted_residues_errors(self,x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        scaled_errors = []
        for i,ctr in enumerate(self.CTRs):
            F_theo = np.abs(self.xtal.F(ctr.harr,ctr.karr,ctr.l))
            if hasattr(ctr,'angles'):
                anglecorr = self.get_anglecorrection(ctr.angles['omega'])
            else:
                anglecorr = 1.
            F_obs_corr = ctr.sfI*anglecorr
            F_err_corr = ctr.err*anglecorr
            if self.scaleindividual or ctr.hk == (0,0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr) # scale CTR
                residues.append((ctr.weight/(F_err_corr*scale)) * (F_obs_corr*scale - F_theo ) )
                scaled_errors.append(F_err_corr*scale)
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            return np.concatenate(residues), np.concatenate(scaled_errors)
        else:
            _, err = self.flat_data(False)
            scale = self.scaling(np.concatenate(F_t), np.concatenate(F_obs), np.asarray(err))
            for i,ctr in enumerate(filter(lambda x: x.hk != (0,0),self.CTRs)):
                residues.append( (ctr.weight/(scale*F_err[i])) * (F_obs[i]*scale - F_t[i] ) )
                scaled_errors.append(F_err[i]*scale)
            return np.concatenate(residues), np.concatenate(scaled_errors)
        
    def statistics(self,x=None):
        
        if x is None:
            x = self.get_parameters() 
        
        #self.xtal.setParameters(x)
        residues2 = self.weighted_residues2(x)
        
        Rfactor = self.Rfactor(x)
        
        stat = dict()
        
        #variance = np.concatenate([ctr.err**2 for ctr in self.CTRs])
        #varmat_i = np.diag(1/variance)
        
        chi2_result = np.sum(residues2)
        pvalue = 1 - stats.chi2.cdf(chi2_result,residues2.size - x.size)
        chi2_red = chi2_result/(residues2.size - x.size)
        
        pcov = util.leastsq_covariance(self.weighted_residues,x)

        errors = np.sqrt(np.diag(pcov)*chi2_red)
        
        self.set_errors(errors)
        self.set_parameters(x)

        
        stat['Chisqr'] = chi2_result
        stat['nodatapoints'] = residues2.size
        stat['Chisqr_red'] = chi2_red
        stat['noparameters'] = x.size
        stat['pvalue'] = pvalue
        stat['Rfactor'] = Rfactor
        stat['covariance'] = pcov
        
        return stat
    
                    
    def set_archi_result(self, archi):
        islandid = int(np.argmin([f[0] for f in archi.get_champions_f()]))
        minisland = archi[islandid]
        pop_min = minisland.get_population()
        
        res = pop_min.champion_x
        chi2_res = pop_min.champion_f
        
        stat = self.statistics(res)
    
        popsize = archi[0].get_population().get_f().shape[0]

        params = { name: np.empty((len(archi), popsize))  for name in self.fitparnames}
        params['chisqr'] = np.empty((len(archi), popsize))
        #params['logpdf'] = np.empty((len(archi), popsize))

        for i, isl in enumerate(archi):
            pop = isl.get_population()
            for j, p in enumerate(self.fitparnames):
                params[p][i] = pop.get_x()[:, j]
            params['chisqr'][i] = pop.get_f()[:,0]
            #params['logpdf'][i] = stats.chi2.logpdf(pop.get_f()[:,0],self._cached_flat_data[0].size - len(self.fitparnames)) # correct???
        
        cov = stat.pop('covariance',np.array([])) # cannot save it as netcdf!
        
        import arviz as az
        fittrace = az.from_dict(params, attrs=stat)
        return fittrace
    
    """
    def defaultCTRplotsettings(self):
        self.lrange = [0.,9.]
        self.plotsize = (19,12)
        self.settings = {linestyle='', marker='.',color='black',zorder=2,elinewidth=0.5,capsize=1.,'markersize' : 2.}
    
    def plotParametersetCTR(self,x,lrange=[0.,9.],plotsize=(19,12)):
        self.xtal.setParameters(x)
        fitCTRs = self.CTRs.generateCollectionFromXtal(self.xtal,1000,lrange)
        
        ctroverviewfig = ctrfigure(figsize=plotsize)
        
        self.CTRs.setPlotSettings(linestyle='', marker='.',color='black',zorder=2,elinewidth=0.5,capsize=1.,markersize=2.)
        self.CTRs.setAllToDefaultID()
        ctroverviewfig.addCollection(self.CTRs)
        #ideal_rods.setPlotSettings(linestyle='-', marker='',color=(0.7,0.7,0.7,1.),zorder=1,linewidth=4)
        #ideal_rods.setAllToDefaultID()
        #ctroverviewfig.addCollection(ideal_rods)
        fitCTRs.setPlotSettings(linestyle='-', marker='',color='red',zorder=3)
        fitCTRs.setAllToDefaultID()
        ctroverviewfig.addCollection(fitCTRs)
        
        ctroverviewfig.settings(wspace=0.02,hspace=0.05,ylabels='|$F$| / a.u.',ylim=[1e-1,1e1],xlim=[-0.05,8.05],xlabels="$L$ / r.l.u.") 
        ctroverviewfig.generateCTRplot(2)
        
    """
        

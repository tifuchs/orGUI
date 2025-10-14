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


import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.transforms as mtransforms
import scipy.interpolate as interp
from scipy.stats.mstats import gmean
import numpy as np 
from matplotlib import rc
import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import math
import json
import warnings
from collections import OrderedDict
from .CTRcalc import SXRDCrystal
from .HKLVlieg import VliegAngles
from .. import util

# don't init CTRFigure directly, use ctrfigure instead
class CTRFigure(mplfig.Figure):
    
    def __init__(self,**figargs):
        super(CTRFigure, self).__init__(**figargs)
        self.data = OrderedDict()
        self.xlabels = "L / r.l.u."
        self.ylabels = "Structure factor / arb. units"
        self.xlim = None #[0,5]
        self.ylim = None #[3e1,2e3]
        self.wspace = 0
        self.hspace = 0.05
        
        
    def settings(self,**settings):
        self.__dict__.update(settings)
        
        
    def addCTR(self,ctr,*plotargs,**keyargs):
        if ctr.ctr_id in self.data:
            self.data[ctr.ctr_id].append([ctr,plotargs,keyargs])
        else:
            self.data[ctr.ctr_id] = [[ctr,plotargs,keyargs]]
            
    def addCollection(self,collection):
        for rod in collection:
            self.addCTR(rod,*collection.plotsett,**collection.plotkeyargs)
        
            
    def generateCTRplot(self,cols=2,maxrelerr=None,**keyargs):
        rows = math.ceil(float(len(self.data))/float(cols))
        ax1 = self.add_subplot( rows,cols,1 )
        [self.add_subplot( rows,cols,i ) for i in range(2,rows*cols+1)]
        self.axes_ctr_id = []
        self.axes_ctr_hk = []
        share = keyargs.get('share_axes', True)
        self.to_sharex = np.array([share for i in range(len(self.axes))])
        self.to_sharey = np.array([share for i in range(len(self.axes))])
        if keyargs.get('sort_hk', False):
            data = sorted(self.data,key=lambda x : np.sum(np.abs(x[0])))
        else:
            data = self.data
        offset = 0
        for i,ctrkey in enumerate(data):
            i += offset
            while(i in keyargs.get('skip_panel', [])):
                i += 1
                offset += 1
                self.axes_ctr_id.append(i)
                self.axes_ctr_hk.append(i)
            ax = self.axes[i]
            ctrname, ctrid = ctrkey
            self.axes_ctr_id.append(ctrkey)
            self.axes_ctr_hk.append(ctrname)
            
            def formatMillerName(h):
                h = round(h, 3)
                if not h % 1:
                    h = int(h)
                if h < 0:
                    return r"\overline{%s}" % abs(h)
                else:
                    return r"%s" % abs(h)
            
            if 'rodLabelL' in keyargs:
                ctrstr = r"$( \, %s \, %s \, \ell \, )$" % (formatMillerName(ctrname[0]),formatMillerName(ctrname[1]))
            else:
                ctrstr = r"$( \,%s \, %s \,)$" % (formatMillerName(ctrname[0]),formatMillerName(ctrname[1]) )
            
            if hasattr(self,"rodlabelsize"):
                size = self.rodlabelsize
            else:
                size = None
            
            if hasattr(self,"rodlabelweight"):
                weight = self.rodlabelweight
            else:
                weight = 'bold'

            if 'rodlabel' in keyargs:
                if keyargs['rodlabel'] == 'bottom':
                    txt = ax.text(keyargs.get('labelxpos',0.02),keyargs.get('labelypos',0.02),ctrstr,horizontalalignment='left',
                        verticalalignment='bottom' ,size=size,transform=ax.transAxes, weight=weight)
                if keyargs['rodlabel'] == 'topright':
                    txt = ax.text(keyargs.get('labelxpos',0.95),keyargs.get('labelypos',0.95),ctrstr,horizontalalignment='right',
                        verticalalignment='top' ,size=size,transform=ax.transAxes, weight=weight)
            else:
                txt = ax.text(keyargs.get('labelxpos',0.05),keyargs.get('labelypos',0.95),ctrstr,horizontalalignment='left',
                        verticalalignment='top' ,size=size,transform=ax.transAxes, weight=weight)
            #txt.set_bbox(dict(facecolor='white', alpha=1., edgecolor=))
            xlim = None
            ylim = None
            for ctr,plotargs,keyargs_c in self.data[ctrkey]:
                xlim = keyargs_c.pop('xlim', xlim)
                ylim = keyargs_c.pop('ylim', ylim)
                if ctr.isWithError:
                    if maxrelerr is not None:
                        relerr = ctr.err / ctr.sfI
                        mask = relerr > maxrelerr
                        ax.errorbar(ctr.l[~mask],ctr.sfI[~mask],yerr=ctr.err[~mask],**keyargs_c)                        
                        if 'plotcaps' in keyargs:
                            errbar = ax.errorbar(ctr.l[mask],ctr.sfI[mask],yerr=ctr.err[mask],**keyargs_c)
                            _,_,barlinecols = errbar.lines
                            barlinecols[0].set_visible(False)
                    else:
                        ax.errorbar(ctr.l,ctr.sfI,yerr=ctr.err,**keyargs_c) # ,fmt='.',elinewidth=0.5,capsize=1. ,errorevery=1,color='k',zorder=1
                else:
                    ax.plot(ctr.l,ctr.sfI,*plotargs,**keyargs_c)
            if self.data[ctrkey][0][0].difference:    
                ax.set_yscale('linear')
            else:
                ax.set_yscale('log')
            if xlim is not None:
                self.to_sharex[i] = False
                ax.set_xlim(xlim)
            elif self.xlim is not None:
                ax.set_xlim(self.xlim)
            
            if ylim is not None:
                self.to_sharey[i] = False
                ax.set_ylim(ylim)
            elif self.xlim is not None:
                ax.set_ylim(self.ylim)

        if np.any(self.to_sharey):
            shareyaxes = np.array(self.axes)[self.to_sharey]
            for axy in shareyaxes:
                if axy is shareyaxes[0]:
                    continue
                axy.sharey(shareyaxes[0])
        
        if np.any(self.to_sharex):
            sharexaxes = np.array(self.axes)[self.to_sharex]
            for axx in sharexaxes:
                if axx is sharexaxes[0]:
                    continue
                axx.sharex(sharexaxes[0])


        axes = np.reshape(self.axes,(rows,cols))
        [a.tick_params(axis='x', labelbottom=False) for a in axes[:-1, :].flat]
        [a.tick_params(axis='y', labelleft=False) for a in axes[:, 1:].flat]
        if isinstance(self.xlabels, list):
            [a.set_ylabel(ylbl) for ylbl, a in zip(self.ylabels, axes[:, 0])]
            [a.set_xlabel(xlbl) for xlbl,a in zip(self.xlabels, axes[-1, :])]
        elif keyargs.get('sharexyLabels', False):
            if keyargs.get('sharexyLabels') == 'x':
                axes[-1,0].set_xlabel(self.xlabels)
                [a.set_ylabel(self.ylabels) for a in axes[:, 0]]
                xy_lablel = 'x'
            elif keyargs.get('sharexyLabels') == 'y':
                axes[0,0].set_ylabel(self.ylabels)
                [a.set_xlabel(self.xlabels) for a in axes[-1, :]]
                xy_lablel = 'y'
            else:
                axes[0,0].set_ylabel(self.ylabels) # create labels for the tight_layout() call below
                axes[-1,0].set_xlabel(self.xlabels)
                xy_lablel = 'xy'
        else:
            [a.set_ylabel(self.ylabels) for a in axes[:, 0]]
            [a.set_xlabel(self.xlabels) for a in axes[-1, :]]

        self.tight_layout()
        self.subplots_adjust(wspace=self.wspace,hspace=self.hspace)
        
        if keyargs.get('sharexyLabels', False):
            if 'x' in xy_lablel:
                avepos_width = 0.5 * (self.subplotpars.left + self.subplotpars.right)
                transform_xlabel = mtransforms.blended_transform_factory(self.transFigure, mtransforms.IdentityTransform())
                axes[-1,0].xaxis.label.set_transform(transform_xlabel)
                axes[-1,0].xaxis.label.set_x(avepos_width)
            if 'y' in xy_lablel:
                avepos_height = 0.5 * (self.subplotpars.bottom + self.subplotpars.top)
                transform_ylabel = mtransforms.blended_transform_factory(mtransforms.IdentityTransform(), self.transFigure) 
                axes[0,0].yaxis.label.set_transform(transform_ylabel)
                axes[0,0].yaxis.label.set_y(avepos_height)
            
        
    def get_ctr_ax(self, hk, ident=None):
        hk = (float(hk[0]), float(hk[1]))
        if ident is None:
            idx = self.axes_ctr_hk.index(hk)
        else:
            idx = self.axes_ctr_id.index((hk,ident))
        return self.axes[idx]
        
    def set_ctr_xlim(self, hk, xlim, ident=None):
        hk = (float(hk[0]), float(hk[1]))
        if self.axes:
            ax = self.get_ctr_ax(hk, ident)
            if len(ax.get_shared_x_axes().get_siblings(ax)) > 1: 
                raise RuntimeError("CTR axes limits must be set before plot creation.")
            else:
                ax.set_xlim(xlim)
        else:
            if ident is None:
                for ids in self.data:
                    if ids[0] == hk:
                        self.data[ids][0][2]['xlim'] = xlim
                        break
                else:
                    raise KeyError("CTR %s not found in figure data" % str(hk))
            else:
                self.data[(hk, ident)][0][2]['xlim'] = xlim
                
                
    def set_ctr_ylim(self, hk, ylim, ident=None):
        hk = (float(hk[0]), float(hk[1]))
        if self.axes:
            ax = self.get_ctr_ax(hk, ident)
            if len(ax.get_shared_y_axes().get_siblings(ax)) > 1: 
                raise RuntimeError("CTR axes limits must be set before plot creation.")
            else:
                ax.set_ylim(ylim)
        else:
            if ident is None:
                for ids in self.data:
                    if ids[0] == hk:
                        self.data[ids][0][2]['ylim'] = ylim
                        break
                else:
                    raise KeyError("CTR %s not found in figure data" % str(hk))
            else:
                self.data[(hk, ident)][0][2]['ylim'] = ylim

                
class CTR(object):
    
    ctrtype = 100
    
    optional_counters = ['bgI', 'ctrI', 'croi_pix', 'bgroi_pix', 'weight']
    
    def __init__(self,hk,l=None,sfI=None,err=None,phi=None,**keyargs):
        self.hk = tuple(hk)
        h,k = hk
        self.l = np.ascontiguousarray(l)
        self.sfI = np.ascontiguousarray(sfI)
        if err is not None:
            self.err = np.ascontiguousarray(err)
        else:
            self.err = None
        self.ctrtype = CTR.ctrtype
        CTR.ctrtype += 1
        self.withErr = True
        self.phi = phi
        self.difference = False
        self.weight = 1
        if 'name' in keyargs:
            self.name = keyargs['name']
        else:
            self.name = 'default'
        self.harr = np.full_like(self.l,h)
        self.karr = np.full_like(self.l,k)
        
    def toNXdict(self):
        nxdict = {
            "@NX_class": u"NXdata",
            "sixc_angles": {
                "@NX_class": u"NXpositioner",
                "@unit" : u"deg"
            },
            "hkl": {
                "@NX_class": u"NXcollection",
                "h" :  self.harr,
                "k" :  self.harr,
                "l" :  self.l,
                "@unit" : u"r.l.u."
            },
            "counters":{
                "@NX_class": u"NXdetector",
                "structurefactor" : self.sfI
            },
            "@signal" : u"counters/structurefactor",
            "@axes": u"hkl/l",
            "@title" : repr(self),
            "@name" : self.name,
            "@difference" : self.difference
        }
        for cnter in CTR.optional_counters:
            if hasattr(self, cnter):
                nxdict["counters"][cnter] = getattr(self, cnter)
                
        if self.err is not None:
            nxdict["counters"]["structurefactor_errors"] = self.err
        
        if self.phi is not None:
            nxdict["counters"]["phase"] = self.phi
            
        if hasattr(self, 'angles'):
            for ang in self.angles.dtype.fields:
                nxdict["sixc_angles"][ang] = self.angles[ang]
        
        return nxdict
    
    @classmethod
    def fromNXdict(cls, nxdict):
        h = nxdict['hkl']['h']
        k = nxdict['hkl']['k']
        l = nxdict['hkl']['l']
        sfI = nxdict['counters']['structurefactor']
        err = nxdict['counters'].get('structurefactor_errors', None)
        phi = nxdict['counters'].get('phase', None)
        name = nxdict.get('@name', 'default')
        
        ctr = cls((h[0], k[0]), l, sfI, err, phi, name=name)
        ctr.harr = h
        ctr.karr = k
        
        if '@difference' in nxdict:
            ctr.difference = nxdict['@difference']
        
        for cnter in CTR.optional_counters:
            if hasattr(nxdict['counters'], cnter):
                setattr(ctr, cnter, nxdict["counters"][cnter])
        
        angles = []
        angles_names = []
        for ang in nxdict['sixc_angles']:
            if not ang.startswith("@"):
                angles.append(nxdict['sixc_angles'][ang])
                angles_names.append(ang)
        if angles:
            dt = np.dtype([(ang, 'f8') for ang in angles_names])
            angles = np.vstack(angles).T
            angles = np.core.records.fromarrays(angles,dtype=dt)
            ctr.angles = angles
        return ctr

        
    def getPlotLabel(self):
        if self.plotlabel in self:
            return self.plotlabel
        else:
            return self.name
    
    def convertToF(self,excludeInvalid=True):
        self.sfI = np.sqrt(self.sfI)
        if excludeInvalid:
            mask = ~np.isnan(self.sfI)
        else:
            mask = np.ones_like(self.sfI,dtype=np.bool_)
        
        if self.isWithError:
            self.err = 0.5*(self.err/self.sfI)[mask]
        if hasattr(self,'bgI'):
            self.bgI = np.sqrt(self.bgI)[mask]
        if hasattr(self,'ctrI'):
            self.ctrI = np.sqrt(self.ctrI)[mask]
        self.sfI = self.sfI[mask]
        self.l = self.l[mask]
        self.harr = self.harr[mask]
        self.karr = self.karr[mask]
        if self.isWithPhase:
            self.phi = self.phi[mask]
        
    # in degrees
    def setPhase(self,phi):
        self.phi = phi
    
    def setWithError(self,err):
        self.withErr = err
    
    @property
    def isWithPhase(self):
        if not isinstance(self.phi,np.ndarray):
            return False
        else:
            return True
        
    @property
    def isWithError(self):
        if not isinstance(self.err,np.ndarray):
            return False
        return self.withErr
    
    def getComplexSF(self):
        if not self.isWithPhase:
            raise Exception("%s:\nNo phase informaion available." % repr(self))
        return self.sfI*np.exp(1j*np.deg2rad(self.phi))
    
    @property
    def ctr_id(self):
        return tuple(np.around(self.hk,2)),self.ctrtype
    
    def setToDefaultID(self):
        self.ctrtype = 0
    
    def generateDifference(self,other):
        otherinter = interp.interp1d(other.l,other.sfI)
        self.sfI -= otherinter(self.l)
        self.difference = True
        return self
    
    def meanSF(self,lowerL,upperL):
        upper = np.nanargmin(np.abs(self.l - upperL))
        lower = np.nanargmin(np.abs(self.l - lowerL))
        return np.nanmean(self.sfI[lower:upper])
        
    def __imul__(self,valOrArray):
        self.sfI *= valOrArray
        if isinstance(self.err,np.ndarray):
            self.err *= valOrArray
        return self
    
    def __iadd__(self,valOrArray):
        self.sfI += valOrArray
        #raise NotImplementedError()
        #if isinstance(self.err,np.ndarray):
        #    self.err += valOrArray
        return self
    
    def cut(self,lower,upper, invert=False):
        """Restricts the CTR to the selected lower and upper index
        
        If invert is True, will remove the data within the selected range.
        
        If invert is set to 'insertNAN', the sfI of the lower index will be set to nan.
        This is useful for CTR plotting to interrupt the lines at this point.
        But should never be used for any CTR that is supposed to be used for 
        further computations!
        
        """
        if invert:
            mask = np.ones_like(self.l,dtype=np.bool_)
            if invert=='insertNAN':
                mask[lower+1:upper] = False
                self.sfI[lower] = np.nan
            else:
                mask[lower:upper] = False
                
        else:
            mask = slice(lower,upper)
        self.l = self.l[mask]
        self.harr = self.harr[mask]
        self.karr = self.karr[mask]
        self.sfI = self.sfI[mask]
        if self.isWithError:
            self.err = self.err[mask]
        if self.isWithPhase:
            self.phi = self.phi[mask]
        if hasattr(self,'bgI'):
            self.bgI = self.bgI[mask]
        if hasattr(self,'ctrI'):
            self.ctrI = self.ctrI[mask]
            
    def cutToL(self,lowerL,upperL, invert=False):
        """Restricts the CTR to the selected lowerL and upperL.
        
        If invert is True, will remove the data within the selected range.
        
        If invert is set to 'insertNAN', the sfI of the lower index will be set to nan.
        This is useful for CTR plotting to interrupt the lines at this point.
        But should never be used for any CTR that is supposed to be used for 
        further computations!
        """
        upper = np.nanargmin(np.abs(self.l - upperL))
        lower = np.nanargmin(np.abs(self.l - lowerL))
        self.cut(lower,upper,invert)
        
        
    def cutToROIfile(self,jsonfile):
        with open(jsonfile,'r') as f:
            udict = json.load(f)
        
        rois = udict['ROI']['roidict']
        del rois['ICR']
        mask = np.zeros_like(self.l,dtype=np.bool_)
        
        for roikey in rois:
            fr = rois[roikey]['from']
            to = rois[roikey]['to']
            lower = np.nanargmin(np.abs(self.l - fr))
            upper = np.nanargmin(np.abs(self.l - to))
            mask[lower:upper] = 1.
        self.l = self.l[mask]
        self.harr = self.harr[mask]
        self.karr = self.karr[mask]
        self.sfI = self.sfI[mask]
        
        if self.isWithError:
            self.err = self.err[mask]
        if self.isWithPhase:
            self.phi = self.phi[mask]
        if hasattr(self,'bgI'):
            self.bgI = self.bgI[mask]
        if hasattr(self,'ctrI'):
            self.ctrI = self.ctrI[mask]
            
    def get_scale(self, xtal, omitErrors=False ,lognorm=False):
        if not hasattr(self,'err') or omitErrors:
            err = None
        else:
            err = self.err
        F_cryst = np.abs( xtal.F(self.harr,self.karr,self.l))
        if lognorm:
            return util.get_scale_logchi2(F_cryst, self.sfI)
        else:
            return util.get_scale_chi2(F_cryst, self.sfI, err)
        
    def scaleToXtal(self,xtal, omitErrors=False, lognorm=False):
        self.__imul__(self.get_scale(xtal,omitErrors,lognorm))

        
    def calcAnglesZmode(self,vliegangles,fixedangle=np.deg2rad(0.1),
                        fixed='in', chi=0.,phi=0., **keyargs):
        l = self.l
        h = self.harr
        k = self.karr
        hkl = np.vstack((h,k,l))
    
        pos = vliegangles.anglesZmode(hkl,fixedangle,fixed='in',chi=0,phi=0,**keyargs)
        dt = np.dtype([('alpha', 'f8'), ('delta','f8'), ('gamma','f8'), ('omega','f8'), ('chi','f8'), ('phi','f8')])
        self.angles = np.core.records.fromarrays(pos.T,dtype=dt)
        return self.angles
        
        
    def toArray(self,mode=None):
        if mode is not None:
            data = np.empty((6,self.l.size))
            data[5] = mode
        else:
            data = np.empty((5,self.l.size))
        data[0] = self.harr
        data[1] = self.karr
        data[2] = self.l
        data[3] = self.sfI
        if self.err is not None: 
            data[4] = self.err
        elif self.phi is not None:
            data[4] = self.phi
        else:
            data[4] = np.nan
            
        return data.T

    #returns a list of CTRs!!!
    @staticmethod
    def fromANAROD(filenameOrArray,RODexport=False):
        warnings.warn("CTR.fromANAROD is deprecated, use CTRCollection.fromANAROD instead!", DeprecationWarning)
        if not isinstance(filenameOrArray,np.ndarray):
            filenameOrArray = np.loadtxt(filenameOrArray,skiprows=1)
        rods = np.unique(filenameOrArray[:,:2],axis=0)
        CTRs = []
        for hk in rods:
            rodmask = np.logical_and(filenameOrArray[:,0] == hk[0],filenameOrArray[:,1] == hk[1])
            rod = filenameOrArray[rodmask]
            l = rod[:,2]
            sfI = rod[:,3]
            if RODexport:
                ctr = CTR(tuple(hk),l,sfI)
                ctr.setPhase(rod[:,4])
                CTRs.append(ctr)
            else:
                err = rod[:,4] if rod.shape[1] > 4 else None
                CTRs.append(CTR(tuple(hk),l,sfI,err))
            
        return CTRCollection(CTRs)
    
    @classmethod
    def fromArray(cls, array,RODexport=False):
        h = array[:,0][0]
        k = array[:,1][0]
        l = array[:,2]
        sfI = array[:,3]
        if RODexport:
            err = None
            phase = array[:,4]
        else:
            err = array[:,4] if array.shape[1] > 4 else None
            phase = None
        return cls([h,k],l,sfI,err,phase)
    
    def millerIdentifier(self):
        h,k = self.hk
        idstr = "%i_%02i_%i_%02i" % (round(h), round((h % 1) * 100),round(k), round((k % 1) * 100))
        return idstr
    
    def rollbackHK(self):
        idstrsp = self.name.split('_')
        h = float(idstrsp[0] + '.' + idstrsp[1])
        k = float(idstrsp[2] + '.' + idstrsp[3])
        self.hk = h,k
        
    def generateAverage(self, step_size=None, **kwargs):
        """Creates a new averaged CTR.
        
        args: step_size:  step size along l
              nbins : number of bins along l
              overlap (float from 0 to 1): define how much the first and last bins exceed the data range.
                                           This can solve issues with edge data points               
        provide either step_size or nbins
        
        """
        overlap = kwargs.get('overlap', 0.25)
        
        lmax = np.amax(self.l)
        lmin = np.amin(self.l)
        size_max = self.l.size
        
        l_range = lmax - lmin
        
        if 'nbins' in kwargs:
            nbins = int(kwargs.get('nbins'))
            #
            l_full_range = l_range / (1. - (2*overlap) / nbins)
            step = l_full_range / nbins
            l_first_bin = lmin - step*overlap
            bin_edges = l_first_bin + step*np.arange(nbins+1)
            
        elif step_size is not None:
            nbins = int(np.floor(l_range / abs(step_size))) + 1
            l_full_range = l_range / (1. - (2*overlap) / nbins)
            
            nbins = int(np.floor(l_full_range / abs(step_size))) + 1
            
            l_first_bin = lmin - step_size*overlap
            bin_edges = l_first_bin + step_size*np.arange(nbins+1)
        else:
            nbins = size_max + 1

            l_full_range = l_range / (1. - (2*overlap) / nbins)
            step = l_full_range / nbins
            
            l_first_bin = lmin - step*overlap
            bin_edges = l_first_bin + step*np.arange(nbins+1)
            
        l_cntr = np.zeros(nbins)
        
        indexes = np.digitize(self.l, bin_edges)
        
        if np.any(indexes == 0) or np.any(indexes == nbins+1):
            raise Exception("bin edges were chosen incorrectly. This is probably a bug.")
            
        indexes -= 1
        
        weights = np.zeros_like(l_cntr)
        np.add.at(weights,indexes,1.)
        
        np.add.at(l_cntr, indexes, self.l)
        l_cntr /= weights
        
        I = np.zeros_like(l_cntr)
        np.add.at(I, indexes, self.sfI)
        I /= weights
        
        if hasattr(self,'err') and self.err is not None:
            Ierr = np.zeros_like(l_cntr)
            np.add.at(Ierr, indexes, self.err**2)
            Ierr = np.sqrt(Ierr)
            Ierr /= weights
        else:
            Ierr = None
        
        mask = np.logical_and(weights != 0, np.isfinite(I))
        
        l_masked = l_cntr[mask]
        I_masked = I[mask]
        if Ierr is not None:
            Ierr_masked = Ierr[mask]
        else:
            Ierr_masked = None
        weights_masked = weights[mask]

        newctr = CTR(self.hk,l_masked,I_masked,final_error)
        newctr.contributions = weights_masked
        return newctr
    
    def __repr__(self):
        #return "<CTR %s ctrtype %s at %016X>" % (tuple(np.around(self.hk,2)), self.ctrtype , id(self))
        return "<CTR<%s> %s ctrtype %s>" % (self.name,tuple(np.around(self.hk,2)), self.ctrtype)
        

class CTRCollection(list):
    def __init__(self,iterableCTRs=[], **kwargs):
        super(CTRCollection, self).__init__(iterableCTRs)
        self.plotsett = ()
        self.plotkeyargs = {'linestyle':'-', 'marker':'.','color':'k','zorder':1}
        #self._updaterodtypes()
        self.name = kwargs.get('name', 'CTRs')
        
    def setPlotSettings(self,*settings,**keyargs):
        self.plotsett = settings
        self.plotkeyargs = keyargs
    
    def getReprList(self):
        return [repr(rod) for rod in self]
        
    def getHKList(self):
        return [rod.hk for rod in self]
    
    def setAllToDefaultID(self):
        for rod in self:
            rod.setToDefaultID()
    
    def deleteRod(self,key):
        if isinstance(key,tuple):
            for rod in self:
                if rod.hk == key:
                    self.remove(rod)
            
    def __and__(self,other):
        coll = CTRCollection()
        for rod in self:
            if repr(rod) in other.getReprList():
                coll.append(rod)
        #coll._updaterodtypes()
        return coll
        
    def convertToF(self):
        for ctr in self:
            ctr.convertToF()
            
    def generateAverage(self, step_size=None, **kwargs):
        """Creates a new CTRCollection with CTRs, which were individually
        averaged.
        
        args: step_size:  step size along l
              nbins : number of bins along l
              overlap (float from 0 to 1): define how much the first and last bins exceed the data range.
                                           This can solve issues with edge data points               
        provide either step_size or nbins
        
        """
        coll = CTRCollection(name='AVE: ' + self.name)
        for ctr in self:
            coll.append(ctr.generateAverage(step_size, **kwargs))
        return coll
    
    def generateDifferenceCollection(self,other,sortby='repr'):
        coll = CTRCollection()
        if isinstance(other,CTRCollection):
            for rod in self:
                if sortby == 'repr':
                    try:
                        idx = other.getReprList().index(repr(rod))
                    except ValueError:
                        continue
                else:
                    try:
                        idx = other.getHKList().index(rod.hk)
                    except ValueError:
                        continue
                coll.append(rod.generateDifference(other[idx]))
            return coll
        elif isinstance(other,SXRDCrystal):
            for ctr in self:
                l = np.copy(ctr.l)
                h,k = ctr.hk
                harr = np.full_like(l,h)
                karr = np.full_like(l,k)
                F_theo = np.abs(other.F(harr,karr,l))
                scale = np.prod(ctr.sfI/F_theo)**(1/ctr.l.size)
                F_theo *= scale
                diffCTR = CTR((h,k), l, ctr.sfI - F_theo,name=ctr.name)
                diffCTR.difference = True
                coll.append(diffCTR)
            return coll
        else:
            raise NotImplementedError("Can not generate difference collection for type %s" % type(other))
        
    def addScalar(self,scalar):
        for rod in self:
            rod += scalar
        return self

    def generateCollectionFromXtal(self,xtal,samples=None,lrange=None):
        CTRs = CTRCollection()
        if samples is not None or lrange is not None:
            if lrange is not None:
                for ctr in self:
                    h,k = ctr.hk
                    l = np.linspace(lrange[0],lrange[1],samples)
                    h = np.full_like(l,h)
                    k = np.full_like(l,k)
                    F =  xtal.F(h,k,l)
                    CTRs.append( CTR( ctr.hk, l, np.abs(F),phi=np.angle(F)))
                return CTRs
            else:
                for ctr in self:
                    h,k = ctr.hk
                    l = np.linspace(np.amin(ctr.l),np.amax(ctr.l),samples)
                    h = np.full_like(l,h)
                    k = np.full_like(l,k)
                    F =  xtal.F(h,k,l)
                    CTRs.append( CTR(ctr.hk, l, np.abs(F),phi=np.angle(F)))
                return CTRs
        else:
            for ctr in self:
                F =  xtal.F(ctr.harr,ctr.karr,ctr.l)
                CTRs.append( CTR(ctr.hk, ctr.l, np.abs(F),phi=np.angle(F)))
            return CTRs
        
    def toANAROD(self,filename, mode=-3):
        if self.__getitem__(0).isWithError:
            header = "H  K  L  F_HKL  errorF  mode"
        else:
            header = "H  K  L  F_HKL  phi"
            
        data_combined = np.vstack([ctr.toArray(mode) for ctr in self])
        
        np.savetxt(filename,data_combined,header=header,fmt='%.5f')
        
                    
    @staticmethod
    def fromANAROD(filenameOrArray,RODexport=False, **kwargs):
        if not isinstance(filenameOrArray,np.ndarray):
            name = kwargs.get('name', os.path.basename(filenameOrArray))
            filenameOrArray = np.loadtxt(filenameOrArray,skiprows=1)
        else:
            name = kwargs.get('name', 'Array_CTR_import')
        rods = np.unique(filenameOrArray[:,:2],axis=0)
        CTRs = []
        for hk in rods:
            rodmask = np.logical_and(filenameOrArray[:,0] == hk[0],filenameOrArray[:,1] == hk[1])
            rod = filenameOrArray[rodmask]
            l = rod[:,2]
            sfI = rod[:,3]
            if RODexport:
                ctr = CTR(tuple(hk),l,sfI)
                ctr.setPhase(rod[:,4])
                CTRs.append(ctr)
            else:
                err = rod[:,4] if rod.shape[1] > 4 else None
                CTRs.append(CTR(tuple(hk),l,sfI,err))
            
        return CTRCollection(CTRs, name=name)
        
    def toNXdict(self):

        nxdict = {
            "@title":u"%s" % self.name,
            "@NX_class": u"NXentry",
            "@name" : u"%s" % self.name
        }
        for ctr in self:
            d = ctr.toNXdict()
            nxdict[d['@title']] = d
            nxdict["@default"] = d['@title']
            
        return nxdict
    
    @classmethod
    def fromNXdict(cls, nxdict):
        ctrs = []
        for dt in nxdict:
            if not dt.startswith("@"):
                ctr = CTR.fromNXdict(nxdict[dt])
                ctrs.append(ctr)
        name = nxdict["@name"]
        
        CTRs = cls(ctrs, name=name)
        return CTRs
        
        
    def get_flat(self):
        h = []; k = []; l = []; F = []
        for ctr in self:
            h.append(ctr.harr)
            k.append(ctr.karr)
            l.append(ctr.l)
            F.append(ctr.sfI)
        return (np.concatenate(h),np.concatenate(k),np.concatenate(l)),np.concatenate(F)
        
    def get_err_flat(self):
        err = []
        for ctr in self:
            err.append(ctr.err)
        return np.concatenate(err)
        
    def get_scale(self, xtal, omitErrors=False ,lognorm=False):
        hkl, F = self.get_flat()
        F_cryst = np.abs( xtal.F(*hkl))
        if omitErrors:
            err = None
        else:
            err = self.get_err_flat()
        if lognorm:
            return util.get_scale_logchi2(F_cryst, F)
        else:
            return util.get_scale_chi2(F_cryst, F, err)
        
        
    def scaleToXtal(self,xtal, individual=True, omitErrors=False, lognorm=False):
        if individual:
            for ctr in self:
                ctr.scaleToXtal(xtal,lognorm, omitErrors)
        else:
            self.__imul__(self.get_scale(xtal,omitErrors,lognorm))
    
    def __imul__(self,valOrArray):
        for rod in self:
            rod *= valOrArray
        return self

    def __getitem__(self,key):
        if isinstance(key,tuple):
            for rod in self:
                if rod.hk == key:
                    return rod
            else:
                raise KeyError("<%s: %s> CTR indices %s not found." % (type(self).__name__, self.name, str(key)))
        return super(CTRCollection, self).__getitem__(key)
        
    def __repr__(self):
        s = "<%s: %s\n" % (type(self).__name__, self.name)
        s += super().__repr__()
        s += ">"
        return s

     
def ctrfigure(**figargs):
    figargs['FigureClass'] = CTRFigure
    fig = plt.figure(**figargs)
    return fig
    
if __name__ == "__main__":
    fig = ctrfigure(figsize=(12,8))
    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V17/CH4977_0V17_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r,linestyle='', marker='.',color='k',zorder=1) for r in ctrlist]
    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V47/CH4977_0V47_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r,linestyle='', marker='.',color='g',zorder=1) for r in ctrlist]
    
    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V72/CH4977_0V72_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r,linestyle='', marker='.',color='b',zorder=1) for r in ctrlist]
    
    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V02/CH4977_0V02_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r,linestyle='', marker='.',color='y',zorder=1) for r in ctrlist]
    
    fig.settings(wspace=0.05,hspace=0,ylabels='|F| / a.u.',ylim=[1e1,1e3]) 
    fig.generateCTRplot()
    fig.show()
    














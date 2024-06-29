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

from . import CTRcalc, HKLVlieg, DetectorCalibration
import numpy as np

import numpy.linalg as LA
from typing import Union


def allowedReflections(xtal: Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell],order: int = 5,**keyargs):
    """
    returns Bragg-reflections of the unit cell or SXRDCrystal
    brute force algorithm... has the advantage that it always works ;)
    
    params:
        xtal: CTRcalc.UnitCell or CTRcalc.SXRDCrystal
        
        order: calculate list of reflections up 
        to max{ h, k, l} = order - 1 (default: 5)
        (i.e. maximum order**3 reflections)
        
    keyargs:
        
        negative: bool take negative values of h k l into account (default: False)
        
        returnF : bool returns structure factor of the unit cell at the Bragg-reflection
        
        hklrange : (h,k,l) : tuple with arrays with the values of h,k,l to compute
        
    returns:
        ndarray shape (n,3) with all n valid reflections 
        sorted via absolute value of momentum transfer
        
        if returnF: tuple of hkls, (F, Q)
        F: structure factor (complex)
        Q: momentum transfer (in inverse Angstroms)
    
    """
    if 'hklrange' in keyargs:
        h,k,l = keyargs['hklrange']
        hkls = np.array(np.meshgrid(h,k,l)).T.reshape(-1,3)
    else:
        if 'negative' in keyargs and keyargs['negative']:
            miller= np.arange(-order+1,order)
        else:
            miller= np.arange(order)
        hkls = np.array(np.meshgrid(miller,miller,miller)).T.reshape(-1,3)


    if isinstance(xtal,CTRcalc.SXRDCrystal):
        F = xtal.uc_bulk.F_uc(*hkls.T)
        G = LA.norm(xtal.uc_bulk.reciprocalVectorCart(hkls).T,axis=1)
    else:
        F = xtal.F_uc(*hkls.T)
        G = LA.norm(xtal.reciprocalVectorCart(hkls).T,axis=1)
    
    indices_valid = np.abs(F) > 1e-3
    
    hkls = hkls[indices_valid]
    F = F[indices_valid]
    G = G[indices_valid]
    
    sortarg = np.argsort(G)
    if 'returnF' in keyargs and keyargs['returnF']:
        return hkls[sortarg], (F[sortarg],G[sortarg])
    else:
        return hkls[sortarg]
    
    
def allowedReflections_G(xtal: Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell],maxQ: float = 12.0,**keyargs):
    if isinstance(xtal,CTRcalc.SXRDCrystal):
        b = xtal.uc_bulk.b
        #_,_,b,_ = xtal.uc_bulk.getLatticeParameters()
    else:
        #_,_,b,_ = xtal.getLatticeParameters()
        b = xtal.b
    hklmax = np.floor(maxQ/b).astype(np.int64)
    if 'negative' in keyargs and keyargs['negative']:
        h = np.arange(-hklmax[0]+1,hklmax[0])
        k = np.arange(-hklmax[1]+1,hklmax[1])
        l = np.arange(-hklmax[2]+1,hklmax[2])
    else:
        h = np.arange(hklmax[0])
        k = np.arange(hklmax[1])
        l = np.arange(hklmax[2])
        
    keyargs['hklrange'] = (h,k,l)
    
    if 'returnF' in keyargs:
        returnF = keyargs['returnF']
    else:
        returnF = False
        
    keyargs['returnF'] = True
    hkls, (F,G) = allowedReflections(xtal,  1, **keyargs)
    if not keyargs.get('singleAxis',False):
        Qinrange = G < maxQ
        hkls = hkls[Qinrange]
        F = F[Qinrange]
        G = G[Qinrange]
    if returnF:
        return hkls, (F,G)
    else:
        return hkls
    
    
    
    
def allowedCTRs(xtal: Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell],order: int = 5,**keyargs):
    if 'returnF' in keyargs and keyargs['returnF']:
        hkls, (F,G) = allowedReflections(xtal,  order, **keyargs)
        hk, indx = np.unique(hkls[:,:2],axis=0,return_index=True)
        F = F[indx]
        G = G[indx]
        return hk, (F,G)
    else:
        return np.unique(allowedReflections(xtal,order, **keyargs)[:,:2],axis=0)
    
def allowedCTRs_G(xtal: Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell],maxQ: float = 12.0,**keyargs):
    keyargs['singleAxis'] = True
    if 'returnF' in keyargs and keyargs['returnF']:
        hkls, (F,G) = allowedReflections_G(xtal,  maxQ, **keyargs)
        hk, indx = np.unique(hkls[:,:2],axis=0,return_index=True)
        hkls = hkls[indx]
        F = F[indx]
        hkls[:,2] = 0.
        if isinstance(xtal,CTRcalc.SXRDCrystal):
            G = LA.norm(xtal.uc_bulk.reciprocalVectorCart(hkls).T,axis=1)
        else:
            G = LA.norm(xtal.reciprocalVectorCart(hkls).T,axis=1)
        
        Qinrange = G < maxQ
        hk = hk[Qinrange]
        F = F[Qinrange]
        G = G[Qinrange]
        return hk, (F,G)
    else:
        hkls = allowedReflections_G(xtal,  maxQ, **keyargs)
        hk, indx = np.unique(hkls[:,:2],axis=0,return_index=True)
        hkls = hkls[indx]
        hkls[:,2] = 0.
        if isinstance(xtal,CTRcalc.SXRDCrystal):
            G = LA.norm(xtal.uc_bulk.reciprocalVectorCart(hkls).T,axis=1)
        else:
            G = LA.norm(xtal.reciprocalVectorCart(hkls).T,axis=1)
        Qinrange = G < maxQ
        hk = hk[Qinrange]
        return hk
        
# pos = [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] (angles)
def thscanBragg(xtal: Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell], 
               ub : HKLVlieg.UBCalculator, 
               fixedangle : float, 
               sxrddetector : DetectorCalibration.Detector2D_SXRD, 
               omrange : tuple, **keyargs):
    """
    don't forget to set the correct azimuthal reference in Detector2D_SXRD!

    Parameters
    ----------
    xtal : Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell]
        DESCRIPTION.
    ub : HKLVlieg.UBCalculator
        DESCRIPTION.
    fixedangle : float
        in rad!!!
    sxrddetector : DetectorCalibration.Detector2D_SXRD
        DESCRIPTION.
    omrange : tuple(float,float)
        in rad!!!
    **keyargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    ommin, ommax = omrange[0], omrange[1]
    chi = keyargs.get('chi',0.)
    phi = keyargs.get('phi',0.)
    fixed = keyargs.get('fixed','in')
    if fixed != 'in':
        raise NotImplementedError("Angle of incidence needed... Have to calculate this here in the future")
    alpha = fixedangle
    keyargs['ub'] = ub
    keyargs['mirror'] = 'both'
    keyargs['negative'] = True

    Qmax = sxrddetector.Qmax # np.abs(sxrddetector.qArray()).max() / 10. # to Angstrom-1

    hkls, angles = anglesAllowedReflections_G(xtal,None, alpha,Qmax,**keyargs)
    
    ommask = np.logical_and(angles[:,3] > ommin ,  angles[:,3] < ommax)
    
    #Lmask = hkls[:,2] > 0.1
    
    #mask = np.logical_and(Lmask,ommask)
    
    mask = ommask
    
    hkls = hkls[mask]
    angles = angles[mask]
    
    delta = angles[:,1]
    gamma = angles[:,2]
    
    yx = sxrddetector.pixelsSurfaceAngles(gamma, delta, alpha)
    
    ymask = np.logical_and(yx[:,0] >= 0, yx[:,0] < sxrddetector.detector.shape[0])
    xmask = np.logical_and(yx[:,1] >= 0, yx[:,1] < sxrddetector.detector.shape[1])
    
    mask = np.logical_and(xmask,ymask)
    
    yx = yx[mask]
    angles = angles[mask]
    hkls = hkls[mask]
            
    return hkls, yx, angles
    
        
def thscanCTRs(xtal: Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell], 
               ub : HKLVlieg.UBCalculator, 
               fixedangle : float, 
               sxrddetector : DetectorCalibration.Detector2D_SXRD, 
               omrange : tuple, **keyargs):
    """
    
    don't forget to set the correct azimuthal reference in Detector2D_SXRD!

    Parameters
    ----------
    xtal : Union[CTRcalc.SXRDCrystal,CTRcalc.UnitCell]
        DESCRIPTION.
    ub : HKLVlieg.UBCalculator
        DESCRIPTION.
    fixedangle : float
        in rad!!!
    sxrddetector : DetectorCalibration.Detector2D_SXRD
        DESCRIPTION.
    omrange : tuple(float,float)
        in rad!!!
    **keyargs : TYPE
        om : omega positions where pixel positions should be calculated
        hk : calculate 

    Returns
    -------
    None.

    """
    
    chi = keyargs.pop('chi',0.)
    phi = keyargs.pop('phi',0.)
    fixed = keyargs.pop('fixed','in')
    if fixed != 'in':
        raise NotImplementedError("Angle of incidence needed... Have to calculate this here in the future")
    alpha = fixedangle
    gamma_range, delta_range = sxrddetector.rangegamdel(alpha)
    gam_min, gam_max = gamma_range
    del_min, del_max = delta_range
    delabsmax = max(abs(del_min),abs(del_max))
    
    Qmaxinplane = ((4*np.pi)/ub.getLambda())*np.sin(delabsmax)
    
    possible_hk = allowedCTRs_G(xtal,Qmaxinplane,negative=True)
    angle_calc = HKLVlieg.VliegAngles(ub)
    
    found_hk = []
    xmirror = []
    
    degrange = np.rad2deg(np.abs(omrange[0] - omrange[1]))
    testom = np.linspace(omrange[0],omrange[1],int(np.ceil(degrange)))

    for hk in possible_hk:
        (L1,gam1,delta1), (L2,gam2,delta2) = angle_calc.hkIntersect(hk,alpha,testom,phi,chi)
        mask1 = ~np.isnan(L1)
        mask2 = ~np.isnan(L2)
        if np.all(~mask1) and np.all(~mask2):
            continue

        yx1 = sxrddetector.pixelsSurfaceAngles(gam1[mask1], delta1[mask1], alpha)
        yx2 = sxrddetector.pixelsSurfaceAngles(gam2[mask2], delta2[mask2], alpha)
        
        ymask1 = np.logical_and(yx1[:,0] >= 0, yx1[:,0] < sxrddetector.detector.shape[0])
        xmask1 = np.logical_and(yx1[:,1] >= 0, yx1[:,1] < sxrddetector.detector.shape[1])
        yxmask1 = np.logical_and(xmask1,ymask1)
        
        ymask2 = np.logical_and(yx2[:,0] >= 0, yx2[:,0] < sxrddetector.detector.shape[0])
        xmask2 = np.logical_and(yx2[:,1] >= 0, yx2[:,1] < sxrddetector.detector.shape[1])
        yxmask2 = np.logical_and(xmask2,ymask2)
        
        if np.any(yxmask1):
            found_hk.append(hk)
            xmirror.append(np.sum(delta1[mask1] < 0.) > delta1[mask1].size/2)
            
        if np.any(yxmask2):
            found_hk.append(hk)
            xmirror.append(np.sum(delta2[mask2] < 0.) > delta2[mask2].size/2)
            
    return found_hk, xmirror
    
# pos = [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] (angles)

def anglesAllowedReflections(xtal,energy, fixedangle,order = 5,**keyargs):
    if 'ub' in keyargs:
        ub = keyargs['ub']
    elif isinstance(xtal,CTRcalc.SXRDCrystal):
        ub = HKLVlieg.UBCalculator(xtal.uc_bulk, energy)
        ub.defaultU() # Grazing incidence, for TSD you could use either chi or phi = 90deg
    else:
        ub = HKLVlieg.UBCalculator(xtal, energy)
        ub.defaultU() # Grazing incidence, for TSD you could use either chi or phi = 90deg
    
    
    angle_calc = HKLVlieg.VliegAngles(ub)
    
    hkls = allowedReflections(xtal,order,**keyargs)
    
    chi = keyargs.pop('chi',0.)
    phi = keyargs.pop('phi',0.)
    
    fixed = keyargs.pop('fixed','in')
    
    mirror = keyargs.pop('mirror',False)
    
    if mirror == 'both':
        rightangles = angle_calc.anglesZmode(hkls.T, fixedangle, fixed,chi,phi,mirrorx=False)
        leftangles = angle_calc.anglesZmode(hkls.T, fixedangle, fixed,chi,phi,mirrorx=True)
        hkls = np.hstack((hkls,hkls)).reshape((-1,3))
        angles = np.hstack((rightangles,leftangles)).reshape((-1,6))
    else:
        angles = angle_calc.anglesZmode(hkls.T, fixedangle, fixed,chi,phi,mirrorx=mirror)
        
    return hkls, angles
    
    
def anglesAllowedReflections_G(xtal,energy, fixedangle,maxQ: float = 12.0,**keyargs):
    if 'ub' in keyargs:
        ub = keyargs['ub']
    elif isinstance(xtal,CTRcalc.SXRDCrystal):
        ub = HKLVlieg.UBCalculator(xtal.uc_bulk, energy)
        ub.defaultU() # Grazing incidence, for TSD you could use either chi or phi = 90deg
    else:
        ub = HKLVlieg.UBCalculator(xtal, energy)
        ub.defaultU() # Grazing incidence, for TSD you could use either chi or phi = 90deg
    
    angle_calc = HKLVlieg.VliegAngles(ub)
    
    hkls = allowedReflections_G(xtal,maxQ,**keyargs)
    
    chi = keyargs.pop('chi',0.)
    phi = keyargs.pop('phi',0.)
    
    fixed = keyargs.pop('fixed','in')
    
    mirror = keyargs.pop('mirror',False)
    
    if mirror == 'both':
        rightangles = angle_calc.anglesZmode(hkls.T, fixedangle, fixed,chi,phi,mirrorx=False)
        leftangles = angle_calc.anglesZmode(hkls.T, fixedangle, fixed,chi,phi,mirrorx=True)
        hkls = np.hstack((hkls,hkls)).reshape((-1,3))
        angles = np.hstack((rightangles,leftangles)).reshape((-1,6))
    else:
        angles = angle_calc.anglesZmode(hkls.T, fixedangle, fixed,chi,phi,mirrorx=mirror)
        
    return hkls, angles
        





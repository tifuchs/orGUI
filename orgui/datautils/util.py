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
import numpy.linalg as LA
from scipy.linalg import svd
import scipy.optimize as opt
from scipy import special
from scipy.optimize._numdiff import approx_derivative
from scipy.stats.mstats import gmean
from scipy.spatial.transform import Rotation

from .xrayutils import HKLVlieg
import warnings

# numdifftools sets TLS to static?
# causes OSError while hdf5 read 
#import numdifftools as nd
from math import *
import os
import copy
import matplotlib.pyplot as plt
from matplotlib import colors as colors
import configparser
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
    
    
def atomlist_to_xyzfile(filename,atomlist):
    """
    

    Parameters
    ----------
    filename : str
        file path, .xyz will not be appended.
    atomlist : structured array
        fields: [name, x ,y ,z]

    Returns
    -------
    None.

    """
    
    with open(filename,'w') as f:
        noatoms = atomlist.shape[0]
        f.write("{}\n".format(noatoms))
        f.write(" \n")
        np.savetxt(f,atomlist,fmt=['%s','%.6f','%.6f','%.6f'])
    
    
    
def formatMeasurement_LatexFloat(x,dx,pad=1,concise=False):
    roundedx, roundederr, prec = roundMeasurement(x,dx)
    if prec < 0:
        prec = 0
    
    if concise:
        valstr = "${:{pad}.{prec}f}".format(float(roundedx),prec=prec,pad=pad)
        errstr = "({:d})$".format(int(roundederr*10**prec))
        
    else:
        valstr = r"${:{pad}.{prec}f} \pm ".format(float(roundedx),prec=prec,pad=pad)    
        errstr = r"{:{pad}.{prec}f}$".format(float(roundederr),prec=prec,pad=pad)
        
    return valstr + errstr


def roundMeasurement(x,dx, fallbackdigits=4):
    if np.abs(dx / x) < 1e-8:
        warnings.warn("Relative error smaller than 1e-8. Use fallback digits!")
        precision = fallbackdigits
        errprec = fallbackdigits
        roundedx = np.round(x/errprec)*errprec
        roundederr = np.ceil(dx/errprec)*errprec
        return roundedx, roundederr, precision
        
    firstdigit = 10**np.floor(np.log10(dx)) #mask first digit != 0

    smaller3 = np.array(dx < 3*firstdigit,dtype=np.float64)
    
    errprec = 10**(np.floor(np.log10(dx)) - smaller3)
    
    roundederr = np.ceil(dx/errprec)*errprec
    
    roundedx = np.round(x/errprec)*errprec
    precision = int(-1*np.log10(errprec))
    
    return roundedx, roundederr, precision

def F_stat(chi1, chi2, p1, p2, n):
    nom = (chi1 - chi2)/(p2-p1)
    denom = chi2/(n - p2)
    return nom/denom
    
def leastsq_covariance(fun,x0):
    J = approx_derivative(fun, x0, method='3-point')
    
    if J.ndim != 2:
        J = np.atleast_2d(J)
    _, s, VT = LA.svd(J, full_matrices=False)
    threshold = np.finfo(float).eps * max(J.shape) * s[0]
    s = s[s > threshold]
    VT = VT[:s.size]
    pcov = np.dot(VT.T / s**2, VT)
    return pcov

def leastsq_covariance_unstable(fun,x0):
    Jfun = approx_derivative(fun, x0, method='3-point')
    J =  Jfun(x0)
    if J.ndim != 2:
        J = np.atleast_2d(J)
    pcov = LA.pinv(np.dot(J.T ,J))
    
    return pcov
    
def leastsq_covariance_lm(fun,x0):
    
    res = opt.least_squares(fun,x0,jac='3-point')
    J = res.jac
    if J.ndim != 2:
        J = np.atleast_2d(J)
    _, s, VT = LA.svd(J, full_matrices=False)
    threshold = np.finfo(float).eps * max(J.shape) * s[0]
    s = s[s > threshold]
    VT = VT[:s.size]
    pcov = np.dot(VT.T / s**2, VT)
    return pcov
    
    
def get_scale_chi2(obs, theo, err=None):
    """returns scale ``s`` so that chi2 is minimized.
    
    chi2 = sum_i ((obs_i - s * theo_i)/err_i)**2    
    """
    
    if err is None:
        err = 1.
    err2 = err**2
    scale = np.sum((obs*theo)/ err2) / np.sum((theo**2)/ err2)
    return scale
    
def get_scale_logchi2(obs, theo):
    """returns scale ``s`` so that chi2 is minimized.
    
    chi2 = sum_i ((log(obs_i) - log(s * theo_i)))**2    
    """
    return gmean(obs/theo)


def getScalefactor_fitted(ctr,xtal):
    h,k = ctr.hk
    harr = np.full_like(ctr.l,h)
    karr = np.full_like(ctr.l,k)
    logI = np.log( xtal.F(harr,karr,ctr.l)**2 )
    
    def minfun(x):
        if x[0] < 0:
            return np.inf
        return np.sum((np.log( (ctr.sfI*x[0])**2 ) - logI )**2)
    res = opt.minimize(minfun,[1.])
    return res.x

def getScalefactor(ctr,xtal):
    h,k = ctr.hk
    harr = np.full_like(ctr.l,h)
    karr = np.full_like(ctr.l,k)
    F_cryst = np.abs( xtal.F(harr,karr,ctr.l))
    return np.prod(F_cryst/ctr.sfI)**(1/ctr.l.size)

def getScalefactor_collection(coll,xtal):
    scale = []
    for ctr in coll:
        scale.append(getScalefactor(ctr,xtal))
    return scale
    
def stepup(z,A,sigma,mu):
    return A*0.5*(1 + special.erf(np.sqrt(0.5)*((z-mu)/(sigma))))

def stepdown(z,A,sigma,mu):
    return A*0.5*(1 - special.erf(np.sqrt(0.5)*((z-mu)/(sigma))))
    
    
def averageCTRs_fcc111(CTRcoll, cutoff=2,nosymmetry_factor=3, pclip=0.3, **keyargs):
    fcc111_lattice = HKLVlieg.Crystal([1./np.sqrt(2),1./np.sqrt(2),np.sqrt(3)],[90.,90.,120.])
    B = np.asarray(fcc111_lattice.getB())
    rods = [B @ np.array([*ctr.hk,0.]) for ctr in CTRcoll]
    rotaxis = np.array([0,0,1])
    
    Bi = np.linalg.inv(B)
    
    equiv = equivalentReflectionsRotation(3,rotaxis,rods)
    
    equiv_hkl = [(Bi @ e.T).T for e in equiv]
    
    equivalent = [[] for _ in range(len(equiv))]
    
    decimals = keyargs.get('decimals', 3)
    
    for ctr in CTRcoll:
        ctrhkl = np.array([*ctr.hk,0.])
        for i, ehkl in enumerate(equiv_hkl):
            if np.any(np.isclose(np.around(ehkl,decimals), np.around(ctrhkl,decimals)).all(axis=1)):
                equivalent[i].append(ctr)
                break
        else:
            raise Exception("Did not find the group of symmetry equivalent rods for %s." % ctr)
    
    return averageCTRs(equivalent, cutoff,nosymmetry_factor,pclip)
                
    

def averageCTRs_fcc100(CTRcoll, cutoff=2,nosymmetry_factor=3, pclip=0.3):
    #temporary fix!!
    #from datautils.xrayutils import CTRplotutil
    #newColl = CTRplotutil.CTRCollection()
    rodnames = []
    for ctr in CTRcoll:
        h, k = ctr.hk
        if abs(h) > abs(k):
            hk = abs(h), abs(k)
        else:
            hk = abs(k), abs(h)
        rodnames.append(hk)
        
    #CTRcoll_copy = copy.deepcopy(CTRcoll)
    rodnames = list(sorted(set(rodnames)))
    equivalent = []
    for n in rodnames:
        rodtpe = [c for c in CTRcoll if np.sum(np.array(n)**2) == np.sum(np.array(c.hk)**2) ]
        equivalent.append(rodtpe)

    return averageCTRs(equivalent, cutoff,nosymmetry_factor,pclip)


# for fcc100
def averageCTRs_Pt100(CTRcoll, cutoff=2,nosymmetry_factor=3, pclip=0.3):
    return averageCTRs_fcc100(CTRcoll, cutoff,nosymmetry_factor,pclip)


def equivalentReflectionsRotation(numberrot, rotaxis, reflections, decimals=3):
    """Returns all unique vectors which are generated by rotation
    of the vectors ``reflections`` around the axis ``rotaxis``.
    
    The vectors are rotated by ``0``, ``(2pi / numberrot) * 2``,
    ``(2pi / numberrot) * 2``... ``(2pi / numberrot) * (numberrot-1)``. 
    
    All vectors (i.e. reflections) are sorted into groups of
    vectors which can be generated by the symmetry operation.
    The vectors can be provided as a list of 2d arrays. The vectors in
    the 2d arrays will remain part of the same group. This allows
    concatenation of multiple symmetry operations.
    
    :param int numberrot: number of rotations (2*pi / n)
    :param np.ndarray rotaxis: rotation axis (shape (3,))
    :param list reflections: vectors to rotate. (List of arrays of shape (n,3))
    :param int decimals: Decimals to round before comparing vectors. 
    """
    rotvec_norm = rotaxis / np.linalg.norm(rotaxis)
    rotation_vectors = [2*np.pi * (i / numberrot) * rotvec_norm for i in range(numberrot)]
    rotations = Rotation.from_rotvec(rotation_vectors)
    
    reflections = [np.atleast_2d(r) for r in reflections]
    
    all_symmetry = []

    for refl in reflections:
        symmetry_reflections = np.vstack([rotations.apply(r) for r in refl])
        _, idx = np.unique(np.around(symmetry_reflections,decimals),return_index=True,axis=0)
        all_symmetry.append(symmetry_reflections[idx])
    
    reduced_symmety = []
    
    # very slow and super ugly code :( !!!
    for i in range(len(all_symmetry)):
        for sym_refl in all_symmetry[i]:
            found = False
            for j in range(i+1,len(all_symmetry)):
                if np.any(np.isclose(np.around(all_symmetry[j],decimals), np.around(sym_refl,decimals)).all(axis=1)):
                    symmstack = np.vstack([all_symmetry[j], sym_refl])
                    _, idx = np.unique(np.around(symmstack,decimals),return_index=True,axis=0)
                    all_symmetry[j] = symmstack[idx]
                    found = True
                    break
            else:
                continue
            if found:
                break
        else:
            reduced_symmety.append(all_symmetry[i])
    return reduced_symmety
    
    
    
def averageCTRs(equivalent_list, cutoff=2,nosymmetry_factor=3, pclip=0.3):
    """Averages symmetry equivalent reflections.
    
    :param list equivalent_list: List with lists of symmetry equivalent CTR  
    :param float cutoff: reflections with ``sf*cutoff > error`` are merked as\
    sufficient quality 
    :param float nosymmetry_factor: If no symmetry equivalent reflection exists\
    will set systematic error to ``meanagreement*sf*nosymmetry_factor``,\
    where meanagreement is either calculated from this rod or uses the\
    agreement of the whole dataset.
    :param float pclip: Sets the minimum systematic error to the ``pclip``\
    -percentile of the agreement of the whole dataset.
    
    """
    
    #temporary fix!!
    from datautils.xrayutils import CTRplotutil
    newColl = CTRplotutil.CTRCollection()
    for e in equivalent_list:
        e_curr = list(sorted(e,key=lambda x : x.sfI.size))
        #print(e_curr)
        l = e_curr[0].l
        err = e_curr[0].err**2
        sfI = e_curr[0].sfI
        count = np.ones_like(l)
        #print(e_curr)
        """
        first calculate the mean structure factor 
        and the resulting statistical error
        """
        for ctrno in range(1,len(e_curr)):
            #print(l)
            ctr = copy.deepcopy(e_curr[ctrno])
            val,i1,i2 = np.intersect1d(l,ctr.l,return_indices=True)
            mask1 = np.zeros_like(l, dtype=np.bool_)
            mask1[i1] = True
            mask2 = np.zeros_like(ctr.l, dtype=np.bool_)
            mask2[i2] = True
            err[mask1] += ctr.err[mask2]**2
            sfI[mask1] += ctr.sfI[mask2]
            count[mask1] += 1
            l = np.concatenate((l,ctr.l[~mask2]))
            err = np.concatenate((err,ctr.err[~mask2]**2))
            sfI = np.concatenate((sfI,ctr.sfI[~mask2]))
            count = np.concatenate((count,np.ones(np.sum(~mask2))))
            
            sortarg = np.argsort(l)
            
            l = l[sortarg]
            err = err[sortarg]
            sfI = sfI[sortarg]
            count = count[sortarg]
            
        sfI /= count # mean structure factor of symmetry equivalent rods
        err = np.sqrt(err) / count # statistical error
        """
        now calculate the variance of symmetry equivalent reflections
        """
        varerr = np.zeros_like(err)
        
        for ctrno in range(1,len(e_curr)):
            #print(l)
            ctr = copy.deepcopy(e_curr[ctrno])
            val,i1,i2 = np.intersect1d(l,ctr.l,return_indices=True)
            mask1 = np.zeros_like(l, dtype=np.bool_)
            mask1[i1] = True
            mask2 = np.zeros_like(ctr.l, dtype=np.bool_)
            mask2[i2] = True
            varerr[mask1] += (sfI[mask1] - ctr.sfI[mask2])**2
            
        varerr = np.sqrt(varerr/count) # systematic error, only averaged, hence the count in the sqrt !
        agreement_factor = varerr/sfI
        no_symmetry_equivalent = varerr == 0.
        if not no_symmetry_equivalent.all():
            meanagreement = np.nanmean((varerr/sfI)[~no_symmetry_equivalent])
            #agreement_factor[:] = meanagreement
            varerr[no_symmetry_equivalent] = meanagreement*sfI[no_symmetry_equivalent]*nosymmetry_factor # no symmetry equivalent reflections available
        else:
            meanagreement = None
        
        #sufficient_qualtity_indicator = sfI*cutoff > err
        
        final_error = np.sqrt(err**2 + varerr**2)
        
        hk = list(reversed(sorted([ctr.hk for ctr in e])))[0]
        newctr = CTRplotutil.CTR(hk,l,sfI,final_error)
        newctr.contributions = count
        newctr.systematicerr = varerr
        newctr.staterr = err
        #newctr.quality_reflection = sufficient_qualtity_indicator
        newctr.agreement_factor = agreement_factor
        newctr.no_symmetry_equivalent = no_symmetry_equivalent
        newctr.meanagreement = meanagreement
        newColl.append(newctr)
    
    agreements = np.concatenate([ctr.agreement_factor for ctr in newColl])
    globalagreement = np.nanmean(agreements)
    clipagreement = np.nanpercentile(agreements,pclip*100)
    for ctr in newColl:
        if ctr.no_symmetry_equivalent.all():
            ctr.systematicerr = globalagreement*ctr.sfI*nosymmetry_factor
            ctr.agreement_factor[:] = globalagreement*nosymmetry_factor
            #ctr.err = np.sqrt(ctr.staterr**2 + ctr.systematicerr**2)
        ctr.effagreement = np.clip(ctr.agreement_factor,clipagreement,None)
        ctr.err = np.sqrt((ctr.sfI**2)*(ctr.effagreement)**2 + ctr.staterr**2)
        ctr.quality_reflection = ctr.sfI*cutoff > ctr.err

    newColl.globalagreement = globalagreement
    newColl.clipagreement = clipagreement
        
    return newColl

def openRodfile(name):
    with open(name,'r') as ctr:
        names = next(ctr).split('\t')[:-1]
        names[0] = names[0][2:]
        formats = [np.float64 for name in names]
        return np.loadtxt(ctr,dtype={'names': names,'formats': formats} )


def makeSurePathExists(path):
    if not os.path.exists(path):
        os.makedirs(path)

def as_ndarray(obj):
    """make sure a float, int, list of floats or ints,
    or tuple of floats or ints, acts as a numpy array
    """
    if isinstance(obj, (float, int)):
        return np.array([obj])
    return np.asarray(obj)
    
def readNDarrayConfig(filename):
    
    header = StringIO()
    
    with open(filename,'r') as f:
        for line in f:
            if not line.startswith('#'):
                break
            header.write(line[2:])
    config = configparser.ConfigParser()
    config.read_string(header.getvalue())
    
    data = np.loadtxt(filename)
    return data.T, config

def rotation_matrix_from_vectors(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix
    

def sumRoi(image,roi):
    xlo,xhi = roi[0]
    ylo,yhi = roi[1]
    image_tmp = image[xlo:xhi,ylo:yhi]
    I = np.nansum(image_tmp)
    pixel = np.count_nonzero(~np.isnan(image_tmp))
    return I, pixel

def calcIntersect(g1,g2):
    g1P1 , g1P2 = g1
    g2P1 , g2P2 = g2
    
    a = g1P1 - g1P2
    b = g2P1 - g2P2
    
    res = LA.solve(np.array([a,-b]),g2P1-g1P1)
    
    return a*res[0] + g1P1
"""
gives delta in:
    
"""

def solveTrigEquation(theta1,theta2,l1,l2):
    
    def Chi2(delta):
        res = function(theta1,theta2,l1,l2,delta)
        #print("%s : %s" % (delta,res))
        return res
    
    res = opt.minimize(Chi2,np.array([0.0]),bounds=[(-np.pi,np.pi)],method='TNC')
    print(res)
    print(np.rad2deg(res.x))
    return res.x

def function(theta1,theta2,l1,l2,delta):
    t1 = cos(theta1 - delta[0])**2 * l1**2 
    t2 = cos(theta1 + delta[0])**2 * l2**2
    t3 = 2 * l1 * l2 * cos(theta1 + delta[0]) * cos(theta1 - delta[0]) * cos(2*theta2)
    t4 = ( (l2**2) * cos(theta1 + delta[0])**2 * sin(2*theta2)**2 )/( cos(delta[0] - theta2)**2 )
    return fabs(t1 + t2 - t3 - t4)

#solveTrigEquation(np.deg2rad(2.),np.deg2rad(76.),1.,1.1)

"""
only single value of th:
"""
def x_rotation(th):
    return np.array(((1., 0., 0.), (0., np.cos(th), -np.sin(th)), (0., np.sin(th), np.cos(th))))


def y_rotation(th):
    return np.array(((np.cos(th), 0., np.sin(th)), (0, 1., 0.), (-np.sin(th), 0., np.cos(th))))


def z_rotation(th):
    return np.array(((np.cos(th), -np.sin(th), 0.), (np.sin(th), np.cos(th), 0.), (0., 0., 1.)))


def rot_trans_matrix(xrot, yrot, zrot, translate):
    rotmat = x_rotation(xrot) @ y_rotation(yrot) @ z_rotation(zrot)
    return np.vstack((rotmat.T,translate)).T

"""
th is now a 1d-array, returns an array of rotation matrices
"""
def x_rotationArray(th):
    matrices = np.zeros((th.size,3,3))
    matrices[:,0,0] = 1
    matrices[:,1,1] = np.cos(th)
    matrices[:,1,2] = -np.sin(th)
    matrices[:,2,1] = np.sin(th)
    matrices[:,2,2] = np.cos(th)
    return matrices
    

def y_rotationArray(th):
    matrices = np.zeros((th.size,3,3))
    matrices[:,1,1] = 1
    matrices[:,0,0] = np.cos(th)
    matrices[:,0,2] = np.sin(th)
    matrices[:,2,0] = -np.sin(th)
    matrices[:,2,2] = np.cos(th)
    return matrices


def z_rotationArray(th):
    matrices = np.zeros((th.size,3,3))
    matrices[:,2,2] = 1
    matrices[:,0,0] = np.cos(th)
    matrices[:,0,1] = -np.sin(th)
    matrices[:,1,0] = np.sin(th)
    matrices[:,1,1] = np.cos(th)
    return matrices

def orthogonal(matrix):
    matrix = np.array(matrix)
    
    SMALL = 1e-4
    
    def normalise(m):
        d = LA.norm(m)
        if d < SMALL:
            raise Exception("Error: can't make matrix orthogonal")
        return m / d
    #print(LA.norm(v1))
    v1 = normalise(matrix[:,0])
    print(matrix[:,0])
    print(LA.norm(v1))
    v2 = normalise(matrix[:,1])
    v3 = normalise(matrix[:,2])
    
    return np.hstack([v1, v2, v3]).A


def calcHighPixel(image,threshold):
    return np.argwhere(image > threshold)
    
    #old:
    highintensity = image > threshold
    highpixel = []

    for x,xrow in enumerate(highintensity):
        for y in range(xrow.size):
            if highintensity[x][y]:
                highpixel.append([y,x])

    highpixel = np.array(highpixel).T
    return highpixel
    
def lines_concatenated(x, m, x0, b0, **kwargs):
    """Series of concatenated lines
    
    y_1 = m_1 * x + b0
    
    or (if f1_x_intersect == True):
    y_1 = m_1 * (x - b0)

    y_n = m_n * x + b_n
    y_n+1 = m_n+1*x + (m_n - m_n+1)*x_n + b_n
    
    :param x np.ndarray: datapoints where to calcualte the lines
    :param m np.ndarray: gradients of lines, shape n 
    :param x0 np.ndarray: switchover point between lines, shape n-1
    :param b0 float: constant offset or x intersect of first line
    :param kwargs: optional: when f1_x_intersect is set to True, the b0 parameter will be the x intersect of the first line. I.e.
    """
    y = np.zeros_like(x, dtype=np.float64)
    
    line_switch_idx = np.argmin(np.abs(x - x0[0]))
    if x[line_switch_idx] - x0[0] < 0.:
        line_switch_idx += 1
    if kwargs.get('f1_x_intersect', False):
        y[:line_switch_idx] = m[0] * (x[:line_switch_idx] - b0)
        b_n = -b0*m[0]
    else:
        y[:line_switch_idx] = m[0] * x[:line_switch_idx] + b0
        b_n = b0
        
    last_idx = line_switch_idx
    m_n = m[0]
    x_n = x0[0]
    
    for i, x_np in enumerate(x0[1:], start=1):
        line_switch_idx = np.argmin(np.abs(x - x_np))
        if x[line_switch_idx] - x_np < 0.:
            line_switch_idx += 1
        m_np = m[i]
        b_np = (m_n - m_np) * x_n + b_n
        y[last_idx:line_switch_idx] = m_np * x[last_idx:line_switch_idx] + b_np
        
        
        #print(i, m[i], last_idx, line_switch_idx)
        m_n, b_n, x_n = m_np, b_np, x_np
        last_idx = line_switch_idx
    
    m_np = m[-1]
    b_np = (m_n - m_np) * x_n + b_n
    y[last_idx:] = m_np * x[last_idx:] + b_np
        
    return y

def plotP3Image(image,vmin=0,vmax=None,cmap='jet',thresholdMarker=None,**keyargs):
    if 'figure' in keyargs:
        fig = keyargs['figure']
        ax = fig.add_subplot(111)
    elif 'axis' in keyargs:
        ax = keyargs['axis']
        fig = ax.get_figure()
    else:
        fig = plt.figure(figsize=(12,14))
        ax = fig.add_subplot(111)

    if not vmax:
        ax.imshow(image,interpolation='none',cmap=plt.get_cmap(cmap),norm=colors.SymLogNorm(linthresh=1,linscale=1,vmin=vmin))
    else:
        ax.imshow(image,interpolation='none',cmap=plt.get_cmap(cmap),norm=colors.SymLogNorm(linthresh=1,linscale=1,vmin=vmin,vmax=vmax))
    
    if thresholdMarker is not None:
        highpixel = calcHighPixel(thresholdMarker)
        if highpixel.size > 0:
            ax.plot(highpixel[0],highpixel[1],'ro')
    ax.set_ylim([1700,0])
    ax.set_xlim([0,1475])
    numrows, numcols = image.shape
    def format_coord(x, y):
        col = int(x + 0.5)
        row = int(y + 0.5)
        if col >= 0 and col < numcols and row >= 0 and row < numrows:
            z = image[row, col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f' % (x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f' % (x, y)
    ax.format_coord = format_coord
    return fig, ax

#def deltaGamma(drr,gamma):
#    return np.arctan((1 + drr) / np.abs(np.cos(gamma))) - gamma

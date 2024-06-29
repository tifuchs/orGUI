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


from numba import njit, prange

import numpy as np


# returns the structure factor of the unit cell
# h,k,l have to be 1d arrays

@njit('c16[:](f8[::1], f8[::1], f8[::1], f8, f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,:,::1], f8[::1] , f8)', nogil=True, cache=True)
def unitcell_F_uc_bulk(h,k,l,atten,
                        basis,
                        f_factors,
                        refHKLTransform,
                        B_mat,
                        R_mat,
                        R_mat_inv,
                        coherentDomainMatrix,
                        coherentDomainOccupancy,
                        volume
                       ):
    F = np.zeros(h.size,dtype=np.complex128)
    f = np.zeros(h.size,dtype=np.complex128)
    hkl = refHKLTransform @ np.vstack((h,k,l))
    Q_cart2 = (B_mat @ hkl)**2
    Q_para2 = np.sum(Q_cart2[:2],axis=0) #squared!!!
    Q_perp2 = Q_cart2[2] #squared!!!
    Q2 = Q_para2 + Q_perp2 #squared!!!
    
    domainmatrix = np.empty((coherentDomainMatrix.shape[0],3,3))
    for i in range(coherentDomainMatrix.shape[0]):
        domainmatrix[i] = R_mat_inv @ coherentDomainMatrix[i,:,:-1] @ R_mat
    
    #domainmatrix = [R_mat_inv @ mat[:,:-1] @ R_mat for mat in coherentDomainMatrix]
    
    for i in range(basis.shape[0]):
        f[:] =  f_factors[i][10] + f_factors[i][11] + 1j*f_factors[i][12]
        for j in range(5):
            f += f_factors[i][j]*np.exp(- f_factors[i][j+5]*Q2)
        f *= np.exp(- (basis[i][4] * Q_para2 + basis[i][5] * Q_perp2)/ (16*np.pi**2))
        f *= basis[i][6]
        for mat, weight, eff_mat in zip(coherentDomainMatrix,coherentDomainOccupancy, domainmatrix):
            xyz_rel = eff_mat @ basis[i][1:4] + mat[:,-1]
            F += weight * f * np.exp(2j*np.pi * np.sum(hkl.T * xyz_rel,axis=1) ) * np.exp(atten*xyz_rel[2])
    return F/volume

# returns the structure factor of the unit cell
# h,k,l have to be 1d arrays 
@njit('c16[:](f8[::1], f8[::1], f8[::1], f8, f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,:,::1], f8[::1] , f8)', nogil=True, cache=True)
def unitcell_F_uc_bulk_direct(h,k,l,atten,
                        basis,
                        f_factors,
                        refHKLTransform,
                        B_mat,
                        R_mat,
                        R_mat_inv,
                        coherentDomainMatrix,
                        coherentDomainOccupancy,
                        volume
                       ):
    F = np.zeros(h.size,dtype=np.complex128)
    f = np.zeros(h.size,dtype=np.complex128)
    hkl =  np.vstack((h,k,l))
    Q_cart2 = (B_mat @ hkl)**2
    Q_para2 = np.sum(Q_cart2[:2],axis=0) #squared!!!
    Q_perp2 = Q_cart2[2] #squared!!!
    Q2 = Q_para2 + Q_perp2 #squared!!!
    
    domainmatrix = np.empty((coherentDomainMatrix.shape[0],3,3))
    for i in range(coherentDomainMatrix.shape[0]):
        domainmatrix[i] = R_mat_inv @ coherentDomainMatrix[i,:,:-1] @ R_mat
    
    #domainmatrix = [R_mat_inv @ mat[:,:-1] @ R_mat for mat in coherentDomainMatrix]
    
    for i in range(basis.shape[0]):
        f[:] =  f_factors[i][10] + f_factors[i][11] + 1j*f_factors[i][12]
        for j in range(5):
            f += f_factors[i][j]*np.exp(- f_factors[i][j+5]*Q2)
        f *= np.exp(- (basis[i][4] * Q_para2 + basis[i][5] * Q_perp2)/ (16*np.pi**2))
        f *= basis[i][6]
        for mat, weight, eff_mat in zip(coherentDomainMatrix,coherentDomainOccupancy, domainmatrix):
            xyz_rel = eff_mat @ basis[i][1:4] + mat[:,-1]
            F += weight * f * np.exp(2j*np.pi * np.sum(hkl.T * xyz_rel,axis=1) ) * np.exp(atten*xyz_rel[2])
    return F/volume

@njit('c16[:](f8[::1], f8[::1], f8[::1], f8, f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,:,::1], f8[::1] , f8)', nogil=True, cache=True)
def unitcell_F_bulk(h,k,l,atten,
                        basis,
                        f_factors,
                        refHKLTransform,
                        B_mat,
                        R_mat,
                        R_mat_inv,
                        coherentDomainMatrix,
                        coherentDomainOccupancy,
                        volume
                       ):
    hkl = refHKLTransform @ np.vstack((h,k,l))
    Fuc = unitcell_F_uc_bulk_direct(hkl[0], hkl[1], hkl[2],atten,
                        basis,
                        f_factors,
                        refHKLTransform,
                        B_mat,
                        R_mat,
                        R_mat_inv,
                        coherentDomainMatrix,
                        coherentDomainOccupancy,
                        volume)
    return Fuc/(1- np.exp(- 2j*np.pi * l - atten ))

@njit('c16[:](f8[::1], f8[::1], f8[::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,::1], f8[:,:,::1], f8[::1] , f8)', nogil=True, cache=True)
def unitcell_F_uc(h,k,l,
                        basis,
                        f_factors,
                        refHKLTransform,
                        B_mat,
                        R_mat,
                        R_mat_inv,
                        coherentDomainMatrix,
                        coherentDomainOccupancy,
                        volume
                       ):
    F = np.zeros(h.size,dtype=np.complex128)
    f = np.zeros(h.size,dtype=np.complex128)
    hkl = refHKLTransform @ np.vstack((h,k,l))
    Q_cart2 = (B_mat @ hkl)**2
    Q_para2 = np.sum(Q_cart2[:2],axis=0) #squared!!!
    Q_perp2 = Q_cart2[2] #squared!!!
    Q2 = Q_para2 + Q_perp2 #squared!!!
    
    domainmatrix = np.empty((coherentDomainMatrix.shape[0],3,3))
    for i in range(coherentDomainMatrix.shape[0]):
        domainmatrix[i] = R_mat_inv @ coherentDomainMatrix[i,:,:-1] @ R_mat
    
    #domainmatrix = [R_mat_inv @ mat[:,:-1] @ R_mat for mat in coherentDomainMatrix]
    
    for i in range(basis.shape[0]):
        f[:] =  f_factors[i][10] + f_factors[i][11] + 1j*f_factors[i][12]
        for j in range(5):
            f += f_factors[i][j]*np.exp(- f_factors[i][j+5]*Q2)
        f *= np.exp(- (basis[i][4] * Q_para2 + basis[i][5] * Q_perp2)/ (16*np.pi**2))
        f *= basis[i][6]
        for mat, weight, eff_mat in zip(coherentDomainMatrix,coherentDomainOccupancy, domainmatrix):
            xyz_rel = eff_mat @ basis[i][1:4] + mat[:,-1]
            F += weight * f * np.exp(2j*np.pi * np.sum(hkl.T * xyz_rel,axis=1) )
    return F/volume


    
    
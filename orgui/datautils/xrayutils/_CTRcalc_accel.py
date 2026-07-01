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


import math

from numba import njit

import numpy as np


@njit(nogil=True, cache=True)
def _unitcell_F_core(
    h,
    k,
    l,
    atten,
    apply_reference_transform,
    apply_attenuation,
    apply_bulk_lattice_sum,
    basis,
    f_factors,
    refHKLTransform,
    B_mat,
    R_mat,
    R_mat_inv,
    coherentDomainMatrix,
    coherentDomainOccupancy,
):
    F = np.empty(h.size, dtype=np.complex128)

    domainmatrix = np.empty((coherentDomainMatrix.shape[0], 3, 3))
    for d in range(coherentDomainMatrix.shape[0]):
        for row in range(3):
            for col in range(3):
                value = 0.0
                for left in range(3):
                    for right in range(3):
                        value += (
                            R_mat_inv[row, left]
                            * coherentDomainMatrix[d, left, right]
                            * R_mat[right, col]
                        )
                domainmatrix[d, row, col] = value

    two_pi = 2.0 * math.pi
    debye_waller_denominator = 16.0 * math.pi**2
    for p in range(h.size):
        if apply_reference_transform:
            hh = (
                refHKLTransform[0, 0] * h[p]
                + refHKLTransform[0, 1] * k[p]
                + refHKLTransform[0, 2] * l[p]
            )
            kk = (
                refHKLTransform[1, 0] * h[p]
                + refHKLTransform[1, 1] * k[p]
                + refHKLTransform[1, 2] * l[p]
            )
            ll = (
                refHKLTransform[2, 0] * h[p]
                + refHKLTransform[2, 1] * k[p]
                + refHKLTransform[2, 2] * l[p]
            )
        else:
            hh = h[p]
            kk = k[p]
            ll = l[p]

        qx = B_mat[0, 0] * hh + B_mat[0, 1] * kk + B_mat[0, 2] * ll
        qy = B_mat[1, 0] * hh + B_mat[1, 1] * kk + B_mat[1, 2] * ll
        qz = B_mat[2, 0] * hh + B_mat[2, 1] * kk + B_mat[2, 2] * ll
        q_para2 = qx * qx + qy * qy
        q_perp2 = qz * qz
        q2 = q_para2 + q_perp2

        amplitude = 0.0 + 0.0j
        for i in range(basis.shape[0]):
            form_factor = (
                f_factors[i, 10]
                + f_factors[i, 11]
                + 1j * f_factors[i, 12]
            )
            for j in range(5):
                form_factor += f_factors[i, j] * math.exp(
                    -f_factors[i, j + 5] * q2
                )
            form_factor *= math.exp(
                -(basis[i, 4] * q_para2 + basis[i, 5] * q_perp2)
                / debye_waller_denominator
            )
            form_factor *= basis[i, 6]

            x = basis[i, 1]
            y = basis[i, 2]
            z = basis[i, 3]
            for d in range(coherentDomainMatrix.shape[0]):
                x_rel = (
                    domainmatrix[d, 0, 0] * x
                    + domainmatrix[d, 0, 1] * y
                    + domainmatrix[d, 0, 2] * z
                    + coherentDomainMatrix[d, 0, 3]
                )
                y_rel = (
                    domainmatrix[d, 1, 0] * x
                    + domainmatrix[d, 1, 1] * y
                    + domainmatrix[d, 1, 2] * z
                    + coherentDomainMatrix[d, 1, 3]
                )
                z_rel = (
                    domainmatrix[d, 2, 0] * x
                    + domainmatrix[d, 2, 1] * y
                    + domainmatrix[d, 2, 2] * z
                    + coherentDomainMatrix[d, 2, 3]
                )
                phase = two_pi * (hh * x_rel + kk * y_rel + ll * z_rel)
                domain_factor = coherentDomainOccupancy[d] * (
                    math.cos(phase) + 1j * math.sin(phase)
                )
                if apply_attenuation:
                    domain_factor *= math.exp(atten * z_rel)
                amplitude += domain_factor * form_factor
        if apply_bulk_lattice_sum:
            denominator_phase = -two_pi * l[p]
            denominator = 1.0 - math.exp(-atten) * (
                math.cos(denominator_phase)
                + 1j * math.sin(denominator_phase)
            )
            amplitude /= denominator
        F[p] = amplitude
    return F


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
                        uc_area
                       ):
    """Return one attenuated bulk-cell amplitude in electrons."""
    return _unitcell_F_core(
        h,
        k,
        l,
        atten,
        True,
        True,
        False,
        basis,
        f_factors,
        refHKLTransform,
        B_mat,
        R_mat,
        R_mat_inv,
        coherentDomainMatrix,
        coherentDomainOccupancy,
    )

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
                        uc_area
                       ):
    """Return one bulk-cell amplitude in electrons without frame conversion."""
    return _unitcell_F_core(
        h,
        k,
        l,
        atten,
        False,
        True,
        False,
        basis,
        f_factors,
        refHKLTransform,
        B_mat,
        R_mat,
        R_mat_inv,
        coherentDomainMatrix,
        coherentDomainOccupancy,
    )

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
                        uc_area
                       ):
    """Return a semi-infinite bulk amplitude in electrons."""
    return _unitcell_F_core(
        h,
        k,
        l,
        atten,
        True,
        True,
        True,
        basis,
        f_factors,
        refHKLTransform,
        B_mat,
        R_mat,
        R_mat_inv,
        coherentDomainMatrix,
        coherentDomainOccupancy,
    )

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
                        uc_area
                       ):
    """Return an unnormalized unit-cell amplitude in electrons."""
    return _unitcell_F_core(
        h,
        k,
        l,
        0.0,
        True,
        False,
        False,
        basis,
        f_factors,
        refHKLTransform,
        B_mat,
        R_mat,
        R_mat_inv,
        coherentDomainMatrix,
        coherentDomainOccupancy,
    )

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


from .HKLVlieg import Lattice
import numpy as np
import numpy.linalg as LA
from .. import util
import xraydb
import warnings
import random
import math
import copy
import dataclasses
from scipy.ndimage import gaussian_filter1d
import os
import re
import importlib
import importlib.util
from pathlib import Path
from scipy.special import erf

# random.seed(45)
from collections import OrderedDict
import errno
import glob

from .element_data import cov_radii_array, rgb_array

from .CTRutil import (
    special_elementcolors,
    ParameterType,
    Parameter,
    _ensure_contiguous,
    next_skip_comment,
    DWtoDisorder,
    readWaasmaier,
    readDispersion,
    atomic_number,
    LinearFitFunctions,
)
from .CTRstacking import (
    LayerCycle,
    LayerState,
    LayerTransition,
    resolve_upper_start,
)

def _import_cpp_accel():
    try:
        return importlib.import_module("orgui.datautils.xrayutils._CTRcalc_cpp")
    except ModuleNotFoundError as package_error:
        repo_root = Path(__file__).resolve().parents[3]
        candidates = sorted((repo_root / "build").glob("cp*/_CTRcalc_cpp*.so"))
        if not candidates:
            raise package_error
        extension_path = candidates[-1]
        spec = importlib.util.spec_from_file_location(
            "_CTRcalc_cpp",
            extension_path,
        )
        if spec is None or spec.loader is None:
            raise package_error
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module


try:
    _CTRcalc_cpp = _import_cpp_accel()
    HAS_CPP_ACCEL = True
except Exception:
    _CTRcalc_cpp = None
    HAS_CPP_ACCEL = False

_CTRcalc_accel = None
_ACCEL_BACKEND_ENV_VAR = "ORGUI_ACCEL_BACKEND"
_NUMBA_ACCEL_IMPORT_ATTEMPTED = False
_NUMBA_ACCEL_IMPORT_ERROR = None
HAS_NUMBA_ACCEL = False

if HAS_CPP_ACCEL:
    CTR_ACCEL_BACKEND = "cpp"
else:
    CTR_ACCEL_BACKEND = "numpy"


def _load_numba_accel():
    global _CTRcalc_accel
    global _NUMBA_ACCEL_IMPORT_ATTEMPTED
    global _NUMBA_ACCEL_IMPORT_ERROR
    global HAS_NUMBA_ACCEL

    if not _NUMBA_ACCEL_IMPORT_ATTEMPTED:
        _NUMBA_ACCEL_IMPORT_ATTEMPTED = True
        try:
            _CTRcalc_accel = importlib.import_module(
                "._CTRcalc_accel", package=__package__
            )
            _NUMBA_ACCEL_IMPORT_ERROR = None
            HAS_NUMBA_ACCEL = True
        except Exception as exc:
            _CTRcalc_accel = None
            _NUMBA_ACCEL_IMPORT_ERROR = exc
            HAS_NUMBA_ACCEL = False
    return _CTRcalc_accel


def ctr_numba_accel_available():
    """Return whether the optional CTR Numba backend can be loaded.

    Calling this function imports the optional Numba backend when it has not
    already been probed.

    :returns:
        ``True`` when the environment opt-in is set and the backend imports.
    :rtype: bool
    """
    return _load_numba_accel() is not None


def set_accel_backend(backend):
    """Select the CTR structure-factor acceleration backend.

    :param str backend:
        One of ``"cpp"``, ``"numba"``, or ``"numpy"``.
    :raises ValueError:
        If the requested backend is unknown or unavailable.
    """
    if backend not in {"cpp", "numba", "numpy"}:
        raise ValueError("CTR acceleration backend must be 'cpp', 'numba', or 'numpy'")
    if backend == "cpp" and not HAS_CPP_ACCEL:
        raise ValueError("CTR C++ acceleration backend is not available")
    if backend == "numba" and _load_numba_accel() is None:
        message = (
            "CTR Numba acceleration backend is not available."
        )
        if _NUMBA_ACCEL_IMPORT_ERROR is not None:
            message += f" Import failed: {_NUMBA_ACCEL_IMPORT_ERROR}"
        raise ValueError(message)
    global CTR_ACCEL_BACKEND
    CTR_ACCEL_BACKEND = backend


_requested_backend = os.environ.get(_ACCEL_BACKEND_ENV_VAR)
if _requested_backend is not None:
    set_accel_backend(_requested_backend.strip().lower())


def ctr_accel_enabled():
    """Return whether CTR calculations use an accelerated backend."""
    return CTR_ACCEL_BACKEND != "numpy"


def _ctr_accel_module():
    if CTR_ACCEL_BACKEND == "numpy":
        return None
    if CTR_ACCEL_BACKEND == "cpp":
        return _CTRcalc_cpp
    if CTR_ACCEL_BACKEND == "numba":
        return _load_numba_accel()
    raise ValueError(
        "CTR acceleration backend must be 'cpp', 'numba', or 'numpy'"
    )


class WaterModel(Lattice, LinearFitFunctions):
    # path = os.path.split(__file__)[0]

    # f_water = np.loadtxt("%s/water_scattering.dat" % path).T
    # f_water[0] *= 4*np.pi/10.
    # f_water[1] *= np.sqrt(18.015)
    # f_water[2] *= np.sqrt(18.015)

    # 2 for measured bulk liquid water
    # 1 for calculated "free atom model"
    # f0_water_fun = interp1d(f_water[0],f_water[2])
    supportedTypes = [
        "step",
        "layered",
        "1layer",
        "1layer_O",
        "layered_O",
        "1layer_s_O",
    ]

    parameterOrder = "zpos/frac layer_spacing/frac sigma_0/frac  sigma_bar/frac rho_0_r"

    parameterLookup = {"z": 0, "spacing": 1, "sigma_0": 2, "sigma_bar": 3, "rho_0": 4}

    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))

    def __init__(self, a, alpha, watertype="step", **keyargs):
        super().__init__(a, alpha)
        LinearFitFunctions.__init__(self)
        self.basis = np.array([0.0, 2.8, 0.5, 1.0, 0.0])
        self.pw = 0.0334
        self.fitparameters = []
        self.relfitparam = []
        self.fitparlimits = []
        self.fitparameters_name = []
        self.relfitparam_name = []
        self.relfitlimits = []
        self.relfitparam_prior = []
        self.fitparameter_prior = []
        self.parameters = {"absolute": [], "relative": []}
        self.basis_0 = None
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None
        if "name" in keyargs:
            self.name = keyargs["name"]
        else:
            self.name = "unnamed water model"
        if watertype not in WaterModel.supportedTypes:
            raise Exception(
                f"Water model {watertype} is not supported, should be one of {WaterModel.supportedTypes}"  # noqa: E501
            )
        self.type = watertype

        self.refHKLTransform = np.identity(3)
        self.refRealTransform = np.identity(3)
        self._pos_absolute = 0.0

    def setReferenceUnitCell(self, uc, rotMatrix=np.identity(3)):
        """Set the reciprocal- and real-space reference coordinate system.

        Input ``h``, ``k``, and ``l`` values are expressed in reciprocal
        lattice units of ``uc`` and transformed into this object's lattice
        before its structure factor is evaluated.

        :param UnitCell uc:
            Reference unit cell.
        :param numpy.ndarray rotMatrix:
            Optional 3-by-3 rotation from the reference crystal frame into
            this object's crystal frame.
        """
        self.refRealTransform[:] = self.R_mat_inv @ rotMatrix @ uc.R_mat
        self.refHKLTransform[:] = self.B_mat_inv @ rotMatrix @ uc.B_mat

    @property
    def pos_absolute(self):
        return self._pos_absolute

    @pos_absolute.setter
    def pos_absolute(self, pos):
        self._pos_absolute = pos

    @property
    def height_absolute(self):
        """Return the water onset height in Angstrom.

        Water is continuous towards positive z and therefore does not advance
        the structural stacking cursor beyond its onset.
        """
        return self.pos_absolute

    @property
    def loc_absolute(self):
        """Return the water onset location in Angstrom."""
        return self.pos_absolute

    @property
    def stacking_height_absolute(self):
        """Return the height passed to components above the water model."""
        return self.pos_absolute

    @property
    def stacking_loc_absolute(self):
        """Return the location passed to components above the water model."""
        return self.pos_absolute

    @property
    def end_layer_number(self):
        """Preserve the structural layer number below continuous water."""
        return getattr(self, "start_layer_number", -1)

    def stack_on(self, below_loc, below_height, below_layer=-1, below_state=None):
        """Place the water onset at the top of the object below it.

        :param float below_loc:
            Reference location of the object below in Angstrom.
        :param float below_height:
            Top height of the object below in Angstrom.
        :param float below_layer:
            Top structural layer identifier of the object below.
        :param below_state:
            Accepted for compatibility with other stackable components.
        """
        self.start_layer_number = below_layer
        self.pos_absolute = below_height

    def F_uc(self, h, k, l):  # noqa: E741
        """Return the water-model structure factor in electrons.

        The result is the unnormalized scattering amplitude associated with
        this model's lateral unit cell. Area conversion between heterogeneous
        components is performed by
        :meth:`~orgui.datautils.xrayutils.CTRcalc.SXRDCrystal.F`.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex structure-factor amplitude in electrons.
        :rtype: numpy.ndarray
        """
        mask = np.logical_and(np.isclose(h, 0), np.isclose(k, 0))
        if not np.any(mask):
            return np.zeros_like(l, dtype=np.complex128)
        a, alpha, b, beta = self.getLatticeParameters()
        # Layered water model '''simga_0 is rms of water profile simga_bar is successive broadening where: simga_j = #sqrt(simga_0**2 + j*sigma_bar**2) spacing is the layer spacing.'''  # noqa: E501
        # self.pw = 0.0334
        # density of water
        # f = np.zeros(h.size,dtype=np.complex128)
        # hkl = (self.refHKLTransform * np.vstack((h[mask],k[mask],l[mask]))).A
        # l = hkl[2]
        # Q_cart2 = ((self._BMatrix * hkl).A)**2
        # Q2 = np.sum(Q_cart2,axis=0) #squared!!!

        # for i in range(2):
        #    f[:] +=  self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
        #    for j in range(5):
        #        f += self.f[i][j]*np.exp(- self.f[i][j+5]*Q2)
        # zwater = WaterModel.f0_water_fun(Q) + self.f #scattering factor of water
        # sigma_bar = self.basis[3]
        # spacing = self.basis[1]
        l_masked = l[mask]
        zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis

        zpos += self.pos_absolute / self.a[2]

        # f = np.empty_like(l_masked,dtype=np.complex128)

        if self.type == "step":
            raise NotImplementedError(
                "This water type uses data, that is not distributed with this version of datautils."  # noqa: E501
            )
            # fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            # F_wat = 1j*self.pw*fwater*self.volume*np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(2*np.pi*l_masked)  # noqa: E501
        elif self.type == "layered":
            raise NotImplementedError(
                "This water type uses data, that is not distributed with this version of datautils."  # noqa: E501
            )
            # fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            # f_w_layer = fwater*self.pw*self.volume*d_layering
            # if A0 == 1.:
            #    sf = np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(1. - np.exp(-2*(np.pi * sigma_bar * l_masked)**2)*np.exp(2j*np.pi * l_masked * d_layering))  # noqa: E501
            #    F_wat = f_w_layer*sf
            # else:
            #    sigma_1 = np.sqrt(sigma_0**2 + sigma_bar**2)
            #    sf = np.exp(-2*(np.pi * sigma_1 * l_masked)**2)/(1. - np.exp(-2*(np.pi * sigma_bar * l_masked)**2)*np.exp(2j*np.pi * l_masked * d_layering))  # noqa: E501
            #    F_wat = f_w_layer*( A0*np.exp(-0.5*(sigma_0*l_masked*b[2])**2) + sf*np.exp(2j*np.pi * l_masked * d_layering))  # noqa: E501
        elif self.type == "layered_O":
            # fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion

            f = np.empty_like(l_masked, dtype=np.complex128)
            f[:] = self.f[0][10] + self.f[0][11] + 1j * self.f[0][12]
            for j in range(5):
                f += self.f[0][j] * np.exp(-self.f[0][j + 5] * (l_masked * b[2]) ** 2)
            # f *= (10/8)

            f_w_layer = f * self.pw * self.volume * d_layering

            if A0 == 1.0:
                sf = np.exp(-2 * (np.pi * sigma_0 * l_masked) ** 2) / (
                    1.0
                    - np.exp(-2 * (np.pi * sigma_bar * l_masked) ** 2)
                    * np.exp(2j * np.pi * l_masked * d_layering)
                )
                F_wat = f_w_layer * sf
            else:
                sigma_1 = np.sqrt(sigma_0**2 + sigma_bar**2)
                sf = np.exp(-2 * (np.pi * sigma_1 * l_masked) ** 2) / (
                    1.0
                    - np.exp(-2 * (np.pi * sigma_bar * l_masked) ** 2)
                    * np.exp(2j * np.pi * l_masked * d_layering)
                )
                F_wat = f_w_layer * (
                    A0 * np.exp(-0.5 * (sigma_0 * l_masked * b[2]) ** 2)
                    + sf * np.exp(2j * np.pi * l_masked * d_layering)
                )

        elif self.type == "1layer":
            raise NotImplementedError(
                "This water type uses data, that is not distributed with this version of datautils."  # noqa: E501
            )
            # fwater = f0_water_fun(l_masked*b[2]) + self.wat_dispersion
            # rho, mu_offset = self._1layer_firstGauss()
            # F_wat = 1j*self.pw*fwater*self.volume*np.exp(-2*(np.pi * sigma_0 * l_masked)**2)/(2*np.pi*l_masked)  # noqa: E501
            # F_wat = F_wat + A0*fwater*np.exp(-0.5*(sigma_0*l_masked*b[2])**2)*np.exp(2j*np.pi * l_masked * mu_offset)  # noqa: E501

        elif self.type == "1layer_O":
            rho, mu_offset = self._1layer_firstGauss(0.00)
            f = np.empty_like(l_masked, dtype=np.complex128)
            f[:] = self.f[0][10] + self.f[0][11] + 1j * self.f[0][12]
            for j in range(5):
                f += self.f[0][j] * np.exp(-self.f[0][j + 5] * (l_masked * b[2]) ** 2)

            f *= 10 / 8
            F_wat = (
                1j
                * self.pw
                * f
                * self.volume
                * np.exp(-2 * (np.pi * sigma_0 * l_masked) ** 2)
                / (2 * np.pi * l_masked)
            )

            F_wat = F_wat + A0 * f * np.exp(
                -0.5 * (sigma_0 * l_masked * b[2]) ** 2
            ) * np.exp(2j * np.pi * l_masked * mu_offset)

        elif self.type == "1layer_s_O":
            rho, mu_offset = self._1layer_firstGauss(0.0)
            f = np.empty_like(l_masked, dtype=np.complex128)
            f[:] = (
                self.f[2][10] + self.f[2][11] + 1j * self.f[2][12]
            )  # self.f[2] is O2-
            for j in range(5):
                f += self.f[2][j] * np.exp(-self.f[2][j + 5] * (l_masked * b[2]) ** 2)

            F_wat = (
                1j
                * self.pw
                * f
                * self.volume
                * np.exp(-2 * (np.pi * sigma_0 * l_masked) ** 2)
                / (2 * np.pi * l_masked)
            )

            F_wat = F_wat + A0 * f * np.exp(-0.5 * (sigma_0 * l_masked * b[2]) ** 2)

        else:
            warnings.warn(f"unknown water type {self.type} !")
            F_wat = np.zeros_like(l_masked)

        F_wat *= np.exp(2j * np.pi * l_masked * zpos)  # position of water

        # temp = np.exp(-0.5*(Q*sigma_bar)**2)*np.exp(1j*Q*spacing)
        # flayer = spacing*(zwater)*pw*self.uc_area*np.exp(-0.5*(Q*sigma_0)**2)
        # F = (np.exp(2j*np.pi*l[mask]*zpos)*(flayer/(1-temp)))/self.volume

        F_water = np.zeros_like(l, dtype=np.complex128)
        F_water[mask] = F_wat
        return F_water

    def setWaterParameters(self, zpos, layer_spacing, sigma_0, sigma_bar, A0):
        self.basis = np.array([zpos, layer_spacing, sigma_0, sigma_bar, A0])

    def setEnergy(self, E):
        self._E = E
        self.lookupScatteringFactors(E)

    def lookupScatteringFactors(self, E):
        self.f = np.empty((3, 13), dtype=np.float64)
        for i, name in enumerate(["O", "H", "O2-"]):
            self.f[i, :11] = readWaasmaier(name)
            self.f[i, 11:] = readDispersion(name, E)
        self.wat_dispersion = (
            self.f[0, 11]
            + 2.0 * self.f[1, 11]
            + 1j * (self.f[0, 12] + 2.0 * self.f[1, 12])
        )

    def optical_profile(self, noUC=30, z_step=None, z_origin=None):
        """Return continuous water optical constants sampled towards positive z.

        :param int noUC:
            Length of the sampled water region in water-model unit cells.
        :param float z_step:
            Uniform z sampling interval in Angstrom. Defaults to this model's
            lattice height.
        :param float z_origin:
            Optional atomistic-grid origin in Angstrom used to align samples.
        :returns:
            ``(N, 3)`` float64 array containing z in Angstrom, delta, and
            beta.
        :rtype: numpy.ndarray
        """
        from .CTRoptics import water_optical_profile

        return water_optical_profile(
            self, noUC=noUC, z_step=z_step, z_origin=z_origin
        )
        # f1,f2 = UnitCell.special_formfactors['H2O'][1](E)
        # self.wat_dispersion = f1 + 1j*f2

    def _1layer_firstGauss(self, biaswidth=0.2547):
        p = [0.85886939, 1.03841354, -1.35438352]
        zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
        # a,alpha,_,_ = self.getLatticeParameters()
        # uc_area = a[0]*a[1]*np.sin(alpha[2])

        rho_0 = A0 / (
            self.volume
            * np.sqrt(2 * np.pi * (sigma_0**2 + (biaswidth / self._a[2]) ** 2))
        )

        A = rho_0 / self.pw

        mu_rel = np.sqrt(2 * np.log(2)) + util.stepdown(np.log10(A), p[0], p[1], p[2])
        mu = mu_rel * np.sqrt(sigma_0**2 + (biaswidth / self._a[2]) ** 2)

        return rho_0 * (10.0 + self.wat_dispersion), mu

    def zDensity_G(self, z, h, k):
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

        if abs(h) > 0.001 or abs(k) > 0.001:
            return np.zeros_like(z)
        zstep = np.diff(z)
        zstep_mean = np.mean(zstep)
        if not np.allclose(zstep, zstep_mean):
            warnings.warn(
                "zDensity: z stepsize is not equal in given z array."
                "This will result in numerical errors in electron density calculation!"
            )
        zpos, d_layering, sigma_0, sigma_bar, A0 = self.basis
        zpos += self.pos_absolute / self.a[2]
        if self.type == "step":
            rho = (
                0.5
                * (10.0 + self.wat_dispersion)
                * self.pw
                * (
                    1.0
                    + erf((z - zpos * self._a[2]) / (np.sqrt(2) * sigma_0 * self._a[2]))
                )
            )
            rho_real = gaussian_filter1d(rho.real, 0.2547 / zstep_mean)
            rho_imag = gaussian_filter1d(rho.imag, 0.2547 / zstep_mean)
            rho = rho_real + 1j * rho_imag  # estimation of water molecular form factor
        elif self.type == "layered":
            zrange_waterstructure = np.amax(z) / self._a[2] - zpos  # lattice units
            nogausseans_inrange = np.ceil(zrange_waterstructure / d_layering)
            nogausseans = int(
                nogausseans_inrange + 20
            )  # some additional gausseans to reduce edge effects
            layer_density = (
                (10.0 + self.wat_dispersion) * self.pw * d_layering * self._a[2]
            )
            rho = np.zeros_like(z, dtype=np.complex128)
            for i in range(nogausseans):
                sigma2_rel = sigma_0**2 + i * sigma_bar**2 + (0.2547 / self._a[2]) ** 2
                if i == 0:
                    rho += (
                        A0
                        * np.exp(
                            -((z / self._a[2] - zpos - i * d_layering) ** 2)
                            / (2 * sigma2_rel)
                        )
                        / (np.sqrt(2 * np.pi * sigma2_rel) * self._a[2])
                    )
                else:
                    rho += np.exp(
                        -((z / self._a[2] - zpos - i * d_layering) ** 2)
                        / (2 * sigma2_rel)
                    ) / (np.sqrt(2 * np.pi * sigma2_rel) * self._a[2])
            rho *= layer_density
        elif self.type == "layered_O":
            zrange_waterstructure = np.amax(z) / self._a[2] - zpos  # lattice units
            nogausseans_inrange = np.ceil(zrange_waterstructure / d_layering)
            nogausseans = int(
                nogausseans_inrange + 5
            )  # some additional gausseans to reduce edge effects
            layer_density = self.pw * d_layering * self._a[2]
            rho = np.zeros_like(z, dtype=np.complex128)

            for i in range(nogausseans):
                sigma2_rel = sigma_0**2 + i * sigma_bar**2
                deltaZ2i = sigma2_rel * self._a[2] ** 2

                z_i = (zpos + i * d_layering) * self._a[2]

                rho_i = (
                    (self.f[0][10] + self.f[0][11] + 1j * self.f[0][12])
                    / (np.sqrt(2 * np.pi * deltaZ2i))
                ) * np.exp(-0.5 * (((z - z_i) ** 2) / deltaZ2i))
                for j in range(5):
                    exp_dz = self.f[0][j + 5] + 0.5 * deltaZ2i
                    rho_i += (self.f[0][j] / (np.sqrt(4.0 * np.pi * exp_dz))) * np.exp(
                        -(((z - z_i) ** 2) / (4 * exp_dz))
                    )
                # rho_i *= np.exp(-2j*np.pi*((z- z_i)/self._a[2]))
                if i == 0:
                    rho += A0 * rho_i * layer_density
                else:
                    rho += rho_i * layer_density

        elif self.type == "1layer":
            sigma2_rel = sigma_0**2 + (0.2547 / self._a[2]) ** 2
            rho = (
                0.5
                * (10.0 + self.wat_dispersion)
                * self.pw
                * (
                    1.0
                    + erf((z - zpos * self._a[2]) / (np.sqrt(2) * sigma_0 * self._a[2]))
                )
            )
            rho_real = gaussian_filter1d(rho.real, 0.2547 / zstep_mean)
            rho_imag = gaussian_filter1d(rho.imag, 0.2547 / zstep_mean)
            rho = rho_real + 1j * rho_imag  # estimation of water molecular form factor
            rho_layer, mu_offset = self._1layer_firstGauss()
            rho += rho_layer * np.exp(
                -((z / self._a[2] - zpos - mu_offset) ** 2) / (2 * sigma2_rel)
            )
        elif self.type == "1layer_O":
            # rho = 0.5*(10.+self.wat_dispersion)*self.pw*(1. + erf((z - zpos*self._a[2])/(np.sqrt(2)* sigma_0 * self._a[2])))  # noqa: E501
            # rho = gaussian_filter1d(np.abs(rho),0.2547 / zstep_mean).astype(np.complex128) #estimation of water molecular form factor  # noqa: E501
            rho_layer, mu_offset = self._1layer_firstGauss(0.00)

            deltaZ2i = (sigma_0 * self._a[2]) ** 2
            z_i = (zpos + mu_offset) * self._a[2]

            rho_i = (
                (self.f[0][10] + self.f[0][11] + 1j * self.f[0][12])
                / (np.sqrt(2 * np.pi * deltaZ2i))
            ) * np.exp(-0.5 * (((z - z_i) ** 2) / deltaZ2i))
            rho_step = util.stepup(
                z,
                self.f[0][10] + self.f[0][11] + 1j * self.f[0][12],
                sigma_0 * self._a[2],
                zpos * self._a[2],
            )
            for j in range(5):
                # exp_dpara = self.f[i][j+5] + 0.5*deltaPara2i
                exp_dz = self.f[0][j + 5] + 0.5 * deltaZ2i
                eff_sigma = np.sqrt(exp_dz / 2)
                rho_step += util.stepup(z, self.f[0][j], eff_sigma, zpos * self._a[2])
                rho_i += (self.f[0][j] / (np.sqrt(4.0 * np.pi * exp_dz))) * np.exp(
                    -(((z - z_i) ** 2) / (4 * exp_dz))
                )
            rho_i *= A0 * np.exp(-2j * np.pi * ((z - z_i) / self._a[2]))

            rho_step *= (10 / 8) * self.pw

            rho = rho_i / self.uc_area + rho_step
        elif self.type == "1layer_s_O":
            rho_layer, mu_offset = self._1layer_firstGauss(0.0)

            deltaZ2i = (sigma_0 * self._a[2]) ** 2
            z_i = (zpos) * self._a[2]

            rho_i = (
                (self.f[2][10] + self.f[2][11] + 1j * self.f[2][12])
                / (np.sqrt(2 * np.pi * deltaZ2i))
            ) * np.exp(-0.5 * (((z - z_i) ** 2) / deltaZ2i))
            rho_step = util.stepup(
                z,
                self.f[2][10] + self.f[2][11] + 1j * self.f[2][12],
                sigma_0 * self._a[2],
                zpos * self._a[2],
            )
            for j in range(5):
                # exp_dpara = self.f[i][j+5] + 0.5*deltaPara2i
                exp_dz = self.f[2][j + 5] + 0.5 * deltaZ2i
                eff_sigma = np.sqrt(exp_dz / 2)
                rho_step += util.stepup(z, self.f[2][j], eff_sigma, zpos * self._a[2])
                rho_i += (self.f[2][j] / (np.sqrt(4.0 * np.pi * exp_dz))) * np.exp(
                    -(((z - z_i) ** 2) / (4 * exp_dz))
                )
            rho_i *= A0 * np.exp(-2j * np.pi * ((z - z_i) / self._a[2]))

            rho_step *= self.pw

            rho = rho_i / self.uc_area + rho_step

        return rho

    def __repr__(self):
        repr_super = super().__repr__()
        return f"{repr_super}\n water model"

    def waterModelToStr(self, showErrors=True):
        param = self.basis
        if (self.errors is not None) and showErrors:
            err = self.errors
            l = []  # noqa: E741
            for t in zip(param, err):
                [l.append(ti) for ti in t]
            return (
                "Water layer  ({:.5f} +- {:.5f})  ({:.5f} +- {:.5f})  ({:.5f} +- {:.5f})  ({:.4f} +- {:.4f})  ({:.4f} +- {:.4f})".format(  # noqa: E501
                    *tuple(l)
                )  # noqa: E501
            )
        else:
            return "{:.5f}     {:.5f}     {:.5f}  {:.5f}  {:.5f}".format(*tuple(param))

    def latticeRODStr(self):
        a, alpha, _, _ = self.getLatticeParameters()
        return "return\n{}\n{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(
            self.type,
            *a,
            *np.rad2deg(alpha),
        )

    def parameterStr(self, showErrors=True):
        """Return water-model parameters as plain text.

        :param bool showErrors:
            Include propagated fit errors when available.
        :returns:
            Serialized parameter line.
        :rtype: str
        """
        return self.waterModelToStr(showErrors=showErrors) + "\n"

    def parameterStrRod(self):
        return self.waterModelToStr(False) + "\n"

    def __str__(self):
        st = repr(self) + "\n"
        st += "id  " + WaterModel.parameterOrder + "\n"
        return st + self.parameterStr()

    def writeWATfile(self, filename):
        with open(filename, "w") as f:
            st = self.latticeRODStr() + "\n" + self.parameterStrRod()
            f.write(st)

    def toRODStr(self):
        return self.latticeRODStr() + "\n" + self.parameterStrRod()

    def toStr(self, showErrors=True):
        """Serialize the water model as plain text.

        :param bool showErrors:
            Include propagated fit errors when available.
        :returns:
            Plain-text water-model representation.
        :rtype: str
        """
        return (
            self.latticeRODStr()
            + "\n"
            + WaterModel.parameterOrder
            + "\n"
            + self.parameterStr(showErrors=showErrors)
        )

    def pos(self):
        return self.basis[0] * self._a[2]

    def pos_cart_error(self):
        z = self.pos()
        error_zeros = np.nan_to_num(self.errors, 0)
        z_error = error_zeros[0] * self._a[2]
        return z, z_error

    def width_error(self):
        ds0 = self.basis[2] * self._a[2]
        ds1 = self.basis[3] * self._a[2]
        errds0 = self.errors[2] * self._a[2]
        errds1 = self.errors[3] * self._a[2]
        return ds0, ds1, errds0, errds1

    @staticmethod
    def fromWATfile(watfile):
        with open(watfile) as f:
            next_skip_comment(f)  # return line
            nline = next_skip_comment(f).split()
            try:
                latticeparams = np.array(nline, dtype=np.float64)
                watertype = "step"  # default behaviour
            except ValueError:
                watertype = nline[0]
                latticeparams = np.array(next_skip_comment(f).split(), dtype=np.float64)
            basis = np.loadtxt(f)
        uc = WaterModel(latticeparams[:3], latticeparams[3:], watertype)
        if basis.size < 5:  # compatibility to old versions
            basis = np.concatenate((basis, [1.0]))
        uc.basis = basis
        uc.basis_0 = np.copy(uc.basis)
        return uc

    @staticmethod
    def fromStr(string):
        with util.StringIO(string) as f:
            # next_skip_comment(f) # return line
            nline = next_skip_comment(f).split()
            try:
                latticeparams = np.array(nline, dtype=np.float64)
                watertype = "step"  # default behaviour
            except ValueError:
                watertype = nline[0]
                latticeparams = np.array(next_skip_comment(f).split(), dtype=np.float64)
            try:
                while True:
                    nline = next_skip_comment(f)
                    sline = nline.split()
                    if sline and "zpos/frac" not in sline:
                        if "+-" in sline:
                            params = re.findall(r"\(([^)]+)", nline)
                            params_array = np.array(
                                [
                                    np.array(p.split("+-"), dtype=np.float64)
                                    for p in params
                                ]
                            ).T
                            basis = params_array[0]
                            errors = params_array[1]
                        else:
                            basis = np.array(sline, dtype=np.float64)
                            errors = None
                        if basis.size < 5:  # compatibility to old versions
                            basis = np.concatenate((basis, [1.0]))
                            if errors is not None:
                                errors = np.concatenate((errors, [np.nan]))
            except StopIteration:
                pass
        uc = WaterModel(latticeparams[:3], latticeparams[3:], watertype)
        uc.basis = basis
        uc.errors = errors
        uc.basis_0 = np.copy(uc.basis)
        return uc


class UnitCell(Lattice):
    _atom_column_names = (
        "Name",
        "x/frac",
        "y/frac",
        "z/frac",
        "iDW",
        "oDW",
        "occup",
        "layerIdx",
    )
    _atom_column_widths = (6, 12, 12, 12, 10, 10, 10, 10)
    _atom_error_column_widths = (6, 24, 24, 24, 22, 22, 22, 14)
    parameterOrder = "".join(
        f"{name:<{width}}" if index == 0 else f"{name:>{width}}"
        for index, (name, width) in enumerate(
            zip(_atom_column_names, _atom_column_widths)
        )
    ).rstrip()

    parameterLookup = {"x": 1, "y": 2, "z": 3, "iDW": 4, "oDW": 5, "occ": 6, "layer": 7}

    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))

    # lookuptable for special materials
    # give two functions: 1. for f0 (Q) , 2. lookup for f1 and f2 (E)
    # special_formfactors = OrderedDict([('H2O' , (f0_water_fun, lambda E : estimateDispersionCompound('H2O', E) ) )])  # noqa: E501
    special_formfactors = OrderedDict()
    # special_eDensity = OrderedDict([('H2O' , eDens_water_fun )])
    special_eDensity = OrderedDict()
    __db = xraydb.get_xraydb()
    if hasattr(__db, "atomic_symbols"):  # XrayDB <= 4.5
        special_onset = len(__db.atomic_symbols)
    else:
        __elems = __db.get_cache("elements")
        special_onset = len([e.element for e in __elems])
    __db.close()

    special_numbers = OrderedDict(
        zip(
            special_formfactors.keys(),
            np.arange(len(special_formfactors)) + special_onset,
        )
    )

    def __init__(self, a, alpha, **keyargs):
        super().__init__(a, alpha)
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

        self.parameters = {"absolute": [], "relative": []}
        self.layer_behaviour = keyargs.get(
            "layer_behavior",
            keyargs.get("layer_behaviour", "ignore"),
        )
        self._explicit_layer_cycle = keyargs.get("layer_cycle")
        transition = keyargs.get("layer_transition")
        if transition is not None and not isinstance(transition, LayerTransition):
            transition = LayerTransition(transition)
        self.layer_transition = transition
        self._start_layer = -1
        self.layerpos = {0.0: 0.0}
        self.basis_0 = np.array([])
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None
        if "name" in keyargs:
            self.name = keyargs["name"]
        else:
            self.name = "unnamed"
        self.refHKLTransform = np.identity(3, dtype=np.float64)
        self.refRealTransform = np.identity(3, dtype=np.float64)
        self.coherentDomainMatrix = [
            np.vstack((np.identity(3).T, np.array([0, 0, 0]))).T
        ]
        self.coherentDomainOccupancy = [1.0]
        self.dw_increase_constraint = np.array([], dtype=np.bool_)
        self._special_formfactors_present = False
        self.symmetry_metadata = keyargs.get("symmetry_metadata", None)

    def parametersToDict(self):
        d = dict()
        d["basis_0"] = self.basis_0
        d["dw_increase_constraint"] = self.dw_increase_constraint
        for p in self.parameters:
            d[p] = dict()
            for i, param in enumerate(self.parameters[p]):
                d[p][f"{i}_{param.name}"] = param.asdict()
        return d

    def parametersFromDict(self, d, override_values=True):
        for p in self.parameters:
            self.parameters[p] = []
        basis_0 = d["basis_0"].astype(np.float64)
        if basis_0.shape[1] == 7:  # add layer parameter
            basis_0 = np.insert(basis_0, basis_0.shape[1], 0, axis=1)
        self.basis_0[:] = basis_0
        self.dw_increase_constraint = d["dw_increase_constraint"].astype(np.bool_)
        for p in self.parameters:
            for dkey in sorted(d[p].keys()):
                self.parameters[p].append(Parameter(**d[p][dkey]))

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn(
                    f"Can not apply fit parameter values to this unitcell: {e}"
                )
        else:
            warnings.warn(
                "Basis values not updated from Parameter values. This can cause a mismatch between basis and fitparameter values!"  # noqa: E501
            )

    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters"""
        no_errors = False
        self.basis[:] = self.basis_0
        self.errors = np.full_like(self.basis, np.nan)
        for par in self.parameters["absolute"]:
            if par.value is not None:
                self.basis[par.indices] = par.value
            else:
                raise ValueError(
                    f"Can not set basis values from parameters. Value of Parameter {par.name} is None."  # noqa: E501
                )
            if par.error is not None:
                self.errors[par.indices] = par.error
            else:
                no_errors = True
        for par in self.parameters["relative"]:
            if par.value is not None:
                value = self._relative_internal_value(par, par.value)
                self.basis[par.indices] += par.factors * value
            else:
                raise ValueError(
                    f"Can not set basis values from parameters. Value of Parameter {par.name} is None."  # noqa: E501
                )
            if par.error is not None:
                self.errors[par.indices] = np.nan_to_num(
                    self.errors[par.indices], nan=0.0
                )
                self.errors[par.indices] += par.factors * par.error
            else:
                no_errors = True

        self._basis_parvalues = np.copy(self.basis)
        if not no_errors:
            self._errors_parvalues = np.copy(self.errors)

    @staticmethod
    def _relative_internal_value(parameter, value):
        """Convert an exposed parameter value to its internal relative delta."""
        wyckoff = parameter.settings.get("wyckoff", {})
        if wyckoff.get("value_kind") == "absolute":
            return value - wyckoff["reference_value"]
        return value

    def clearParameters(self):
        for p in self.parameters:
            self.parameters[p] = []
        self._basis_parvalues = None
        self._errors_parvalues = None

    def setReferenceUnitCell(self, uc, rotMatrix=np.identity(3)):
        """Set the reciprocal- and real-space reference coordinate system.

        :param UnitCell uc:
            Unit cell defining input reciprocal lattice units.
        :param numpy.ndarray rotMatrix:
            Optional 3-by-3 rotation from the reference crystal frame into
            this unit cell's crystal frame.
        """
        self.refRealTransform[:] = self.R_mat_inv @ rotMatrix @ uc.R_mat
        self.refHKLTransform[:] = self.B_mat_inv @ rotMatrix @ uc.B_mat

    def _test_special_formfactors(self):
        for name in self.names:
            if name in UnitCell.special_formfactors:
                self._special_formfactors_present = True
                return True
        else:
            self._special_formfactors_present = False
            return False

    def addAtom(
        self, element_or_param, xyz_rel=None, iDW=None, oDW=None, occu=None, layer=0
    ):
        if xyz_rel is not None:
            if not isinstance(element_or_param, str):
                raise ValueError(
                    f"You must provide a atom symbol. Got: {element_or_param}"
                )
            atom_no = atomic_number(element_or_param)

            parameters = np.array(
                [atom_no, *xyz_rel, iDW, oDW, occu, layer], dtype=np.float64
            )
            self.names.append(element_or_param)
        else:
            if not isinstance(element_or_param[0], str):
                raise ValueError(
                    f"You must provide a atom symbol. Got: {element_or_param[0]}"
                )
            self.names.append(element_or_param[0])
            atom_no = atomic_number(element_or_param[0])
            element_or_param[0] = atom_no
            parameters = np.array(element_or_param, dtype=np.float64)

        if self.basis.size > 0:
            self.basis = np.vstack([self.basis, parameters])
        else:
            self.basis = np.array([parameters])
        if self.basis_0.size > 0:
            self.basis_0 = np.vstack([self.basis_0, parameters])
        else:
            self.basis_0 = np.array([parameters])
        if float(layer) not in self.layerpos:
            self.layerpos[float(layer)] = 0.0  # default to 0.
        self.errors = None
        self.dw_increase_constraint = np.insert(
            self.dw_increase_constraint, len(self.dw_increase_constraint), True
        )
        self._test_special_formfactors()
        return

    def insertAtom(
        self,
        index,
        element_or_param,
        xyz_rel=None,
        iDW=None,
        oDW=None,
        occu=None,
        layer=0,
    ):
        if xyz_rel is not None:
            if not isinstance(element_or_param, str):
                raise ValueError(
                    f"You must provide a atom symbol. Got: {element_or_param}"
                )
            atom_no = atomic_number(element_or_param)

            parameters = np.array(
                [atom_no, *xyz_rel, iDW, oDW, occu, layer], dtype=np.float64
            )
            self.names.insert(element_or_param, index)
        else:
            if not isinstance(element_or_param[0], str):
                raise ValueError(
                    f"You must provide a atom symbol. Got: {element_or_param[0]}"
                )
            self.names.insert(element_or_param[0], index)
            atom_no = atomic_number(element_or_param[0])
            element_or_param[0] = atom_no
            parameters = np.array(element_or_param, dtype=np.float64)

        if self.basis.size > 0:
            self.basis = np.insert(self.basis, index, parameters, axis=0)
        else:
            self.basis = np.array([parameters])
        if self.basis_0.size > 0:
            self.basis_0 = np.insert(self.basis_0, index, parameters, axis=0)
        else:
            self.basis_0 = np.array([parameters])
        if float(layer) not in self.layerpos:
            self.layerpos[float(layer)] = 0.0  # default to 0.
        self.dw_increase_constraint = np.insert(
            self.dw_increase_constraint, index, True
        )
        self.errors = None
        self._test_special_formfactors()
        # newfp = []
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

    def split_in_layers(self, ordered=True):
        layer_numbers = np.sort(np.unique(self.basis[:, 7]))
        layers = OrderedDict()
        if not hasattr(self, "f"):
            self.setEnergy(
                10000.0
            )  # populate f with some values to enable in-place modification

        for l in layer_numbers:  # noqa: E741
            where = (self.basis[:, 7] == l).nonzero()[0]
            idx_low = where[0]
            idx_high = where[-1] + 1
            uc = UnitCell(self.a, np.rad2deg(self.alpha), name=self.name + f"_layer{l}")
            uc.basis = self.basis[idx_low:idx_high]
            if self.errors is not None:
                uc.errors = self.errors[idx_low:idx_high]
            uc.names = self.names[idx_low:idx_high]
            uc.basis_0 = self.basis_0[idx_low:idx_high]
            uc.dw_increase_constraint = self.dw_increase_constraint[idx_low:idx_high]
            uc._test_special_formfactors()
            uc.f = self.f[idx_low:idx_high]
            if hasattr(self, "_E"):
                uc._E = self._E
            uc.refRealTransform = self.refRealTransform
            uc.refHKLTransform = self.refHKLTransform
            layers[l] = uc
        if self._explicit_layer_cycle is not None:
            cycle = tuple(self._explicit_layer_cycle)
            missing = set(layers) - set(cycle)
            if missing:
                raise ValueError(
                    f"layer_cycle omits unit-cell layers {sorted(missing)!r}"
                )
            layers = OrderedDict((n, layers[n]) for n in cycle)
        elif ordered and len(layer_numbers) > 1:
            avg_height = []
            uc_nms = []
            for uc_l in layers:
                avg_height.append(np.mean(layers[uc_l].basis[:, 3]))
                uc_nms.append(uc_l)
            avg_height = np.array(avg_height)

            sort_uc = np.argsort(avg_height)
            uc_nms = np.array(uc_nms)[sort_uc]

            layers = OrderedDict([(n, layers[n]) for n in uc_nms])

        ordered_layers = tuple(layers)
        origins = []
        for layer in ordered_layers:
            origin = self.layerpos.get(float(layer))
            if origin is None:
                atom_mask = self.basis[:, 7] == layer
                origin = float(np.mean(self.basis[atom_mask, 3]))
            origins.append(float(origin))
        for index, layer in enumerate(ordered_layers):
            spacing = origins[(index + 1) % len(ordered_layers)] - origins[index]
            while spacing <= 0.0:
                spacing += 1.0
            layers[layer]._optical_layer_origin = origins[index]
            layers[layer]._optical_layer_thickness_fraction = spacing

        return layers

    def optical_profile(self):
        """Return homogeneous optical constants for every layer and domain.

        The returned ``(N, 3)`` float64 array has columns ``z`` in Angstrom,
        ``delta``, and ``beta``. Each row is a domain-transformed layer
        contribution with coherent-domain occupancy already applied.

        :returns:
            One row for every structural layer and coherent domain.
        :rtype: numpy.ndarray
        :raises ValueError:
            If :meth:`setEnergy` has not populated the anomalous scattering
            factors.
        """
        from .CTRoptics import optical_profile

        return optical_profile(self)

    def optical_profile_asbulk(self, noUC=30):
        """Return a finite representation of the semi-infinite bulk profile.

        The unit cell is repeated towards negative z exactly as in
        :meth:`zDensity_G_asbulk`.

        :param int noUC:
            Number of unit cells to represent below the termination.
        :returns:
            ``(N, 3)`` float64 array containing z in Angstrom, delta, and
            beta.
        :rtype: numpy.ndarray
        """
        from .CTRoptics import optical_profile_asbulk

        return optical_profile_asbulk(self, noUC=noUC)

    @property
    def layers(self):
        return np.sort(np.unique(self.basis[:, 7]))

    @property
    def layer_cycle(self):
        """Return the ordered local structural-layer cycle."""
        if self._explicit_layer_cycle is not None:
            return LayerCycle(self._explicit_layer_cycle)
        ordered_layers = self.split_in_layers(ordered=True)
        return LayerCycle(tuple(ordered_layers))

    @layer_cycle.setter
    def layer_cycle(self, layers):
        self._explicit_layer_cycle = tuple(layers)

    @staticmethod
    def _translation_from_affine_input(translation):
        translation = np.asarray(translation, dtype=np.float64)
        if translation.shape == (3,):
            return translation
        if translation.shape == (3, 4):
            linear = translation[:, :3]
            vector = translation[:, 3]
        elif translation.shape == (4, 4):
            if not np.allclose(translation[3], [0.0, 0.0, 0.0, 1.0]):
                raise ValueError(
                    "Homogeneous affine transforms must end with "
                    "[0, 0, 0, 1]."
                )
            linear = translation[:3, :3]
            vector = translation[:3, 3]
        else:
            raise ValueError(
                "translation must be a vector with shape (3,), "
                "a 3-by-4 matrix, or a homogeneous 4-by-4 matrix."
            )
        if not np.allclose(linear, np.identity(3)):
            raise ValueError(
                "UnitCell.translate_layered currently supports "
                "translations only; the affine linear part must be identity."
            )
        return vector

    @staticmethod
    def _remap_parameter_atom_indices(parameter, old_to_new):
        if not isinstance(parameter.indices, tuple):
            return
        atoms, parindexes = parameter.indices
        atoms = np.asarray(atoms, dtype=np.intp)
        parameter.indices = (old_to_new[atoms], parindexes)

    @staticmethod
    def _update_absolute_parameter_values_after_basis_transform(
        transformed,
        original_basis,
        old_to_new,
    ):
        new_to_old = np.empty_like(old_to_new)
        new_to_old[old_to_new] = np.arange(old_to_new.size, dtype=np.intp)
        for parameter in transformed.parameters["absolute"]:
            if parameter.value is None or not isinstance(parameter.indices, tuple):
                continue
            atoms, parindexes = parameter.indices
            atoms = np.asarray(atoms, dtype=np.intp)
            old_atoms = new_to_old[atoms]
            old_values = original_basis[(old_atoms, parindexes)]
            new_values = transformed.basis[(atoms, parindexes)]
            shifts = np.asarray(new_values - old_values, dtype=np.float64)
            if np.allclose(shifts, shifts.flat[0]):
                parameter.value += float(shifts.flat[0])
            else:
                parameter.value = None

    @staticmethod
    def _transform_symmetry_metadata(
        metadata,
        old_to_new,
        layer_map,
        xy_translation,
        z_shift_by_layer,
    ):
        if metadata is None:
            return None
        metadata = copy.deepcopy(metadata)
        coordinate_index = {"x": 0, "y": 1, "z": 2}
        atoms = []
        site_parent_shifts = {}
        for atom in metadata.atoms:
            z_shift = z_shift_by_layer[atom.layer]
            coordinate_shift = np.asarray(
                [xy_translation[0], xy_translation[1], z_shift],
                dtype=np.float64,
            )
            couplings = []
            for coupling in atom.couplings:
                couplings.append(
                    dataclasses.replace(
                        coupling,
                        atom_index=int(old_to_new[coupling.atom_index]),
                        constant=(
                            coupling.constant
                            + coordinate_shift[coordinate_index[coupling.coordinate]]
                        ),
                    )
                )
            site_couplings = []
            for coupling in atom.site_couplings:
                site_couplings.append(
                    dataclasses.replace(
                        coupling,
                        atom_index=int(old_to_new[coupling.atom_index]),
                    )
                )
            parent_shift = metadata.surface_spec.transform @ coordinate_shift
            site_parent_shifts.setdefault(atom.site_id, parent_shift)
            atoms.append(
                dataclasses.replace(
                    atom,
                    atom_index=int(old_to_new[atom.atom_index]),
                    parent_fractional=(
                        np.asarray(atom.parent_fractional, dtype=np.float64)
                        + parent_shift
                    ),
                    surface_fractional=(
                        np.asarray(atom.surface_fractional, dtype=np.float64)
                        + coordinate_shift
                    ),
                    layer=layer_map[atom.layer],
                    couplings=tuple(couplings),
                    site_couplings=tuple(site_couplings),
                )
            )
        sites = []
        for site in metadata.sites:
            representative = site.representative_parent_fractional
            if representative is not None and site.site_id in site_parent_shifts:
                representative = tuple(
                    np.asarray(representative, dtype=np.float64)
                    + site_parent_shifts[site.site_id]
                )
            sites.append(
                dataclasses.replace(
                    site,
                    representative_parent_fractional=representative,
                )
            )
        metadata.sites = tuple(sites)
        metadata.atoms = sorted(atoms, key=lambda atom: atom.atom_index)
        return metadata

    def _layer_positions_for_cycle(self, cycle):
        layer_positions = {
            layer: float(self.layerpos.get(float(layer), np.nan)) for layer in cycle
        }
        for layer, position in layer_positions.items():
            if np.isnan(position):
                where = self.basis[:, 7] == layer
                layer_positions[layer] = float(np.mean(self.basis[where, 3]))
        return layer_positions

    def _validated_layered_translation(self, translation):
        translation = self._translation_from_affine_input(translation)
        xy_translation = translation[:2]
        if not np.allclose(xy_translation, np.rint(xy_translation)):
            raise ValueError("x and y affine translations must be integers.")
        if not np.isclose(translation[2], np.rint(translation[2])):
            raise ValueError("z affine translation must be an integer layer step.")
        return np.rint(xy_translation).astype(np.float64), int(
            np.rint(translation[2])
        )

    @staticmethod
    def _shift_layered_basis_values(
        values,
        cycle,
        xy_translation,
        z_shift_by_layer,
        layer_map=None,
    ):
        if values is None or values.size == 0:
            return values
        values_new = np.array(values, copy=True)
        values_new[:, 1:3] += xy_translation
        for layer in cycle:
            where = values[:, 7] == layer
            if not np.any(where):
                continue
            values_new[where, 3] += z_shift_by_layer[layer]
            if layer_map is not None:
                values_new[where, 7] = layer_map[layer]
        return values_new

    def _finish_layered_transform(
        self,
        transformed,
        old_to_new,
        layer_map,
        xy_translation,
        z_shift_by_layer,
        test_special_formfactors=False,
    ):
        for parameters in transformed.parameters.values():
            for parameter in parameters:
                self._remap_parameter_atom_indices(parameter, old_to_new)
        self._update_absolute_parameter_values_after_basis_transform(
            transformed,
            self.basis,
            old_to_new,
        )
        transformed.symmetry_metadata = self._transform_symmetry_metadata(
            self.symmetry_metadata,
            old_to_new,
            layer_map,
            xy_translation,
            z_shift_by_layer,
        )
        if test_special_formfactors:
            transformed._test_special_formfactors()
        return transformed

    def _apply_reordered_layered_transform(
        self,
        transformed,
        cycle,
        xy_translation,
        layer_map,
        z_shift_by_layer,
    ):
        basis = self._shift_layered_basis_values(
            self.basis,
            cycle,
            xy_translation,
            z_shift_by_layer,
            layer_map,
        )
        basis_0 = self._shift_layered_basis_values(
            self.basis_0,
            cycle,
            xy_translation,
            z_shift_by_layer,
            layer_map,
        )
        basis_parvalues = self._shift_layered_basis_values(
            self._basis_parvalues,
            cycle,
            xy_translation,
            z_shift_by_layer,
            layer_map,
        )

        order_lookup = {layer: index for index, layer in enumerate(cycle)}
        row_order = np.array(
            sorted(range(basis.shape[0]), key=lambda i: (order_lookup[basis[i, 7]], i)),
            dtype=np.intp,
        )
        old_to_new = np.empty_like(row_order)
        old_to_new[row_order] = np.arange(row_order.size)

        transformed.basis = basis[row_order]
        transformed.basis_0 = basis_0[row_order]
        if basis_parvalues is not None:
            transformed._basis_parvalues = basis_parvalues[row_order]
        if self.errors is not None:
            transformed.errors = np.array(self.errors, copy=True)[row_order]
        if self._errors_parvalues is not None:
            transformed._errors_parvalues = np.array(
                self._errors_parvalues,
                copy=True,
            )[row_order]
        transformed.names = [self.names[i] for i in row_order]
        transformed.dw_increase_constraint = self.dw_increase_constraint[row_order]
        if hasattr(self, "f"):
            transformed.f = self.f[row_order]

        return self._finish_layered_transform(
            transformed,
            old_to_new,
            layer_map,
            xy_translation,
            z_shift_by_layer,
            test_special_formfactors=True,
        )

    def translate_layered(self, translation, name=None):
        """Return a translated copy with cyclic z-layer wrapping.

        The translation is expressed in unit-cell fractional coordinates. The
        ``x`` and ``y`` translations must be integer unit-cell offsets, and
        ``z`` is an integer number of structural-layer steps. Atom fractional
        z coordinates and layer origins are shifted by ``z / len(layer_cycle)``
        while layer identifiers are reassigned through the ordered
        :class:`CTRstacking.LayerCycle`. A positive z step moves lower layers
        upward and wraps the former top layer to the bottom.

        Symmetry metadata is updated by remapping atom indices, atom layers,
        surface coordinates, parent coordinates, and Wyckoff coupling
        constants.

        :param translation:
            Translation vector with shape ``(3,)``, affine matrix with shape
            ``(3, 4)``, or homogeneous affine matrix with shape ``(4, 4)``.
        :param str name:
            Optional name for the returned unit cell. Defaults to this unit
            cell's name.
        :returns:
            Transformed unit-cell copy.
        :rtype:
            UnitCell
        :raises ValueError:
            If the affine linear part is not identity or any requested
            translation is not an integer.
        """
        xy_translation, z_steps = self._validated_layered_translation(translation)
        transformed = copy.deepcopy(self)
        if name is not None:
            transformed.name = name

        if self.basis.size == 0:
            if z_steps:
                raise ValueError("Cannot cyclically shift layers in an empty UnitCell.")
            return transformed

        cycle = self.layer_cycle.layers
        layer_count = len(cycle)
        z_steps_effective = z_steps % layer_count
        z_translation = z_steps / layer_count

        if z_steps_effective == 0 and np.all(xy_translation == 0):
            if z_steps == 0:
                return transformed
        if z_steps_effective == 0:
            z_shift_by_layer = {layer: z_translation for layer in cycle}
            transformed.basis = self._shift_layered_basis_values(
                transformed.basis,
                cycle,
                xy_translation,
                z_shift_by_layer,
            )
            transformed.basis_0 = self._shift_layered_basis_values(
                transformed.basis_0,
                cycle,
                xy_translation,
                z_shift_by_layer,
            )
            transformed._basis_parvalues = self._shift_layered_basis_values(
                transformed._basis_parvalues,
                cycle,
                xy_translation,
                z_shift_by_layer,
            )
            for layer in cycle:
                layer_key = float(layer)
                if layer_key not in transformed.layerpos:
                    where = self.basis[:, 7] == layer
                    transformed.layerpos[layer_key] = float(
                        np.mean(self.basis[where, 3])
                    )
                transformed.layerpos[layer_key] += z_translation
            old_to_new = np.arange(self.basis.shape[0], dtype=np.intp)
            layer_map = {layer: layer for layer in cycle}
            return self._finish_layered_transform(
                transformed,
                old_to_new,
                layer_map,
                xy_translation,
                z_shift_by_layer,
            )

        layer_positions = self._layer_positions_for_cycle(cycle)
        layer_map = {
            layer: cycle[(index + z_steps_effective) % layer_count]
            for index, layer in enumerate(cycle)
        }
        z_shift_by_layer = {layer: z_translation for layer in cycle}
        transformed = self._apply_reordered_layered_transform(
            transformed,
            cycle,
            xy_translation,
            layer_map,
            z_shift_by_layer,
        )

        transformed.layerpos = copy.deepcopy(self.layerpos)
        for layer in cycle:
            new_layer = layer_map[layer]
            transformed.layerpos[float(new_layer)] = (
                layer_positions[layer] + z_translation
            )
        if transformed._explicit_layer_cycle is not None:
            transformed._explicit_layer_cycle = tuple(cycle)
        return transformed

    def affine_layer_transform(self, translation, name=None):
        """Return a canonical layered translation copy.

        The translation input follows :meth:`translate_layered`, but z layer
        steps are wrapped back into this unit cell's structural-layer
        positions. This preserves a ``0 <= z < 1`` unit-cell representation
        for cells whose layer origins are inside that interval.

        Symmetry metadata is updated by remapping atom indices, atom layers,
        surface coordinates, parent coordinates, and Wyckoff coupling
        constants.

        :param translation:
            Translation vector with shape ``(3,)``, affine matrix with shape
            ``(3, 4)``, or homogeneous affine matrix with shape ``(4, 4)``.
        :param str name:
            Optional name for the returned unit cell. Defaults to this unit
            cell's name.
        :returns:
            Transformed unit-cell copy with wrapped layer coordinates.
        :rtype:
            UnitCell
        :raises ValueError:
            If the affine linear part is not identity or any requested
            translation is not an integer.
        """
        xy_translation, z_steps = self._validated_layered_translation(translation)
        transformed = copy.deepcopy(self)
        if name is not None:
            transformed.name = name

        if self.basis.size == 0:
            if z_steps:
                raise ValueError("Cannot cyclically shift layers in an empty UnitCell.")
            return transformed

        cycle = self.layer_cycle.layers
        layer_count = len(cycle)
        z_steps_effective = z_steps % layer_count
        if z_steps_effective == 0:
            z_shift_by_layer = {layer: 0.0 for layer in cycle}
            transformed.basis = self._shift_layered_basis_values(
                transformed.basis,
                cycle,
                xy_translation,
                z_shift_by_layer,
            )
            transformed.basis_0 = self._shift_layered_basis_values(
                transformed.basis_0,
                cycle,
                xy_translation,
                z_shift_by_layer,
            )
            transformed._basis_parvalues = self._shift_layered_basis_values(
                transformed._basis_parvalues,
                cycle,
                xy_translation,
                z_shift_by_layer,
            )
            old_to_new = np.arange(self.basis.shape[0], dtype=np.intp)
            layer_map = {layer: layer for layer in cycle}
            return self._finish_layered_transform(
                transformed,
                old_to_new,
                layer_map,
                xy_translation,
                z_shift_by_layer,
            )

        layer_positions = self._layer_positions_for_cycle(cycle)
        layer_map = {
            layer: cycle[(index + z_steps_effective) % layer_count]
            for index, layer in enumerate(cycle)
        }
        z_shift_by_layer = {
            layer: layer_positions[layer_map[layer]] - layer_positions[layer]
            for layer in cycle
        }
        transformed = self._apply_reordered_layered_transform(
            transformed,
            cycle,
            xy_translation,
            layer_map,
            z_shift_by_layer,
        )

        transformed.layerpos = copy.deepcopy(self.layerpos)
        for layer in cycle:
            transformed.layerpos[float(layer)] = layer_positions[layer]
        if transformed._explicit_layer_cycle is not None:
            transformed._explicit_layer_cycle = tuple(cycle)
        return transformed

    def as_surface_termination(self, layer, name=None, origin=None):
        """Return a copy whose complete structure is one surface termination.

        The atomic ``layer`` column is used only as a stacking/termination
        selector in the returned cell.  All atoms therefore receive the same
        layer identifier, while their fractional coordinates and fit
        parameters remain unchanged.  This permits a multi-layer relaxed slab
        to be selected as one termination with ``layer_behavior='select'``.

        :param float layer:
            Layer identifier in the primitive Film stacking cycle.
        :param str name:
            Optional name for the returned unit cell.
        :param float origin:
            Fractional z coordinate in this cell that is placed at the exposed
            terrace height.  By default the highest structural-layer origin is
            used.
        :returns:
            Independent termination-cell copy.
        :rtype: UnitCell
        """
        transformed = copy.deepcopy(self)
        if name is not None:
            transformed.name = name
        layer = float(layer)
        if origin is None:
            origins = list(self.layerpos.values())
            if origins:
                origin = max(origins)
            elif self.basis.size:
                origin = float(np.max(self.basis[:, 3]))
            else:
                origin = 0.0
        origin = float(origin)

        for values in (
            transformed.basis,
            transformed.basis_0,
            transformed._basis_parvalues,
        ):
            if values is not None and values.size:
                values[:, 7] = layer
        if transformed.symmetry_metadata is not None:
            transformed.symmetry_metadata.atoms = tuple(
                dataclasses.replace(atom, layer=layer)
                for atom in transformed.symmetry_metadata.atoms
            )
        transformed.layerpos = {layer: origin}
        transformed._explicit_layer_cycle = (layer,)
        transformed.layer_behavior = "select"
        transformed._start_layer = layer
        return transformed

    def supercell(self, repeats, symmetry="preserve", name=None):
        """Return a repeated unit-cell copy.

        :param repeats:
            Integer repeat counts along the ``a``, ``b``, and ``c`` lattice
            directions.
        :param str symmetry:
            ``"preserve"`` keeps generated atoms on the original Wyckoff
            sites, so a later Wyckoff fit parameter is shared across all
            repeated copies. ``"independent"`` creates one copied Wyckoff site
            per generated unit-cell repeat, so each copy can be fitted
            independently.
        :param str name:
            Optional name for the returned unit cell. Defaults to this unit
            cell's name.
        :returns:
            Repeated unit-cell copy.
        :rtype:
            UnitCell
        :raises ValueError:
            If repeat counts are not positive integers or ``symmetry`` is
            unknown.
        """
        repeat_values = np.asarray(repeats, dtype=np.float64)
        if (
            repeat_values.shape != (3,)
            or not np.all(np.isfinite(repeat_values))
            or np.any(repeat_values <= 0)
            or not np.all(repeat_values == np.rint(repeat_values))
        ):
            raise ValueError("repeats must contain three positive integers.")
        repeats = repeat_values.astype(np.intp)
        if symmetry not in {"preserve", "independent"}:
            raise ValueError("symmetry must be 'preserve' or 'independent'.")

        new_a = self.a * repeats
        transformed = UnitCell(
            new_a,
            np.rad2deg(self.alpha),
            name=self.name if name is None else name,
            layer_behavior=self.layer_behavior,
        )
        transformed.layer_transition = copy.deepcopy(self.layer_transition)
        transformed.refHKLTransform = np.copy(self.refHKLTransform)
        transformed.refRealTransform = np.copy(self.refRealTransform)
        transformed.coherentDomainMatrix = []
        repeats_float = repeats.astype(np.float64)
        for matrix in self.coherentDomainMatrix:
            matrix_new = np.array(matrix, copy=True)
            matrix_new[:, -1] = matrix_new[:, -1] / repeats_float
            transformed.coherentDomainMatrix.append(matrix_new)
        transformed.coherentDomainOccupancy = copy.deepcopy(
            self.coherentDomainOccupancy
        )
        transformed._special_formfactors_present = self._special_formfactors_present
        if hasattr(self, "_E"):
            transformed._E = self._E

        if self.basis.size == 0:
            return transformed

        cycle = self.layer_cycle.layers
        layer_positions = self._layer_positions_for_cycle(cycle)
        layer_ids = self._supercell_layer_ids(cycle, int(repeats[2]))

        rows = []
        for iz in range(int(repeats[2])):
            for iy in range(int(repeats[1])):
                for ix in range(int(repeats[0])):
                    cell_offset = np.asarray([ix, iy, iz], dtype=np.float64)
                    copy_index = (
                        iz * int(repeats[0]) * int(repeats[1])
                        + iy * int(repeats[0])
                        + ix
                    )
                    for atom_index, row in enumerate(self.basis):
                        layer = row[7]
                        layer_index = cycle.index(layer)
                        new_row = np.array(row, copy=True)
                        new_row[1:4] = (row[1:4] + cell_offset) / repeats
                        new_row[7] = layer_ids[(iz, layer)]
                        rows.append(
                            {
                                "basis": new_row,
                                "basis_0": np.array(
                                    self.basis_0[atom_index], copy=True
                                ),
                                "name": self.names[atom_index],
                                "constraint": self.dw_increase_constraint[atom_index],
                                "old_atom": atom_index,
                                "copy_index": copy_index,
                                "cell_offset": cell_offset,
                                "layer_order": iz * len(cycle) + layer_index,
                            }
                        )
                        rows[-1]["basis_0"][1:4] = (
                            self.basis_0[atom_index, 1:4] + cell_offset
                        ) / repeats
                        rows[-1]["basis_0"][7] = layer_ids[(iz, layer)]
                        if self.errors is not None:
                            rows[-1]["errors"] = np.array(
                                self.errors[atom_index], copy=True
                            )
                        if self._errors_parvalues is not None:
                            rows[-1]["errors_parvalues"] = np.array(
                                self._errors_parvalues[atom_index], copy=True
                            )
                        if self._basis_parvalues is not None:
                            rows[-1]["basis_parvalues"] = np.array(
                                self._basis_parvalues[atom_index], copy=True
                            )
                            rows[-1]["basis_parvalues"][1:4] = (
                                self._basis_parvalues[atom_index, 1:4] + cell_offset
                            ) / repeats
                            rows[-1]["basis_parvalues"][7] = layer_ids[(iz, layer)]
                        if hasattr(self, "f"):
                            rows[-1]["f"] = np.array(self.f[atom_index], copy=True)

        rows = sorted(
            enumerate(rows),
            key=lambda item: (item[1]["layer_order"], item[0]),
        )
        rows = [row for _, row in rows]

        transformed.basis = np.vstack([row["basis"] for row in rows])
        transformed.basis_0 = np.vstack([row["basis_0"] for row in rows])
        transformed.names = [row["name"] for row in rows]
        transformed.dw_increase_constraint = np.asarray(
            [row["constraint"] for row in rows],
            dtype=np.bool_,
        )
        if self.errors is not None:
            transformed.errors = np.vstack([row["errors"] for row in rows])
        if self._errors_parvalues is not None:
            transformed._errors_parvalues = np.vstack(
                [row["errors_parvalues"] for row in rows]
            )
        if self._basis_parvalues is not None:
            transformed._basis_parvalues = np.vstack(
                [row["basis_parvalues"] for row in rows]
            )
        if hasattr(self, "f"):
            transformed.f = np.vstack([row["f"] for row in rows])

        transformed.layerpos = {}
        for iz in range(int(repeats[2])):
            for layer in cycle:
                transformed.layerpos[float(layer_ids[(iz, layer)])] = (
                    layer_positions[layer] + iz
                ) / repeats[2]
        transformed._explicit_layer_cycle = tuple(
            layer_ids[(iz, layer)]
            for iz in range(int(repeats[2]))
            for layer in cycle
        )
        transformed.parameters = {"absolute": [], "relative": []}
        transformed.basis_0 = np.copy(transformed.basis)
        transformed._basis_parvalues = np.copy(transformed.basis)
        transformed.symmetry_metadata = self._supercell_symmetry_metadata(
            rows,
            repeats.astype(np.float64),
            layer_ids,
            symmetry,
        )
        transformed._test_special_formfactors()
        return transformed

    def _supercell_layer_ids(self, cycle, z_repeats):
        if z_repeats == 1:
            return {(0, layer): layer for layer in cycle}
        numeric = np.asarray(cycle, dtype=np.float64)
        expected = np.arange(1, len(cycle) + 1, dtype=np.float64)
        if np.allclose(numeric, expected):
            return {
                (iz, layer): float(layer + iz * len(cycle))
                for iz in range(z_repeats)
                for layer in cycle
            }
        return {
            (iz, layer): float(iz * len(cycle) + index)
            for iz in range(z_repeats)
            for index, layer in enumerate(cycle)
        }

    def _supercell_symmetry_metadata(self, rows, repeats, layer_ids, symmetry):
        if self.symmetry_metadata is None:
            return None
        metadata = copy.deepcopy(self.symmetry_metadata)
        site_by_copy = {}
        if symmetry == "preserve":
            sites = tuple(copy.deepcopy(metadata.sites))
        else:
            sites = []
            for row in rows:
                atom = metadata.atom_wyckoff_metadata(row["old_atom"])
                if atom is None:
                    continue
                key = (row["copy_index"], atom.site_id)
                if key in site_by_copy:
                    continue
                site = next(
                    site for site in metadata.sites if site.site_id == atom.site_id
                )
                site_id = f"{site.site_id}_copy{row['copy_index']}"
                site_by_copy[key] = site_id
                representative = site.representative_parent_fractional
                if representative is not None:
                    representative = tuple(
                        np.asarray(representative, dtype=np.float64)
                        + metadata.surface_spec.transform @ row["cell_offset"]
                    )
                sites.append(
                    dataclasses.replace(
                        site,
                        site_id=site_id,
                        representative_parent_fractional=representative,
                    )
                )
            sites = tuple(sites)

        atoms = []
        coordinate_index = {"x": 0, "y": 1, "z": 2}
        for new_index, row in enumerate(rows):
            atom = metadata.atom_wyckoff_metadata(row["old_atom"])
            if atom is None:
                continue
            site_id = atom.site_id
            if symmetry == "independent":
                site_id = site_by_copy[(row["copy_index"], atom.site_id)]
            surface_fractional = (
                np.asarray(atom.surface_fractional, dtype=np.float64)
                + row["cell_offset"]
            ) / repeats
            parent_fractional = (
                np.asarray(atom.parent_fractional, dtype=np.float64)
                + metadata.surface_spec.transform @ row["cell_offset"]
            )
            couplings = []
            for coupling in atom.couplings:
                new_site_id = site_id
                couplings.append(
                    dataclasses.replace(
                        coupling,
                        atom_index=new_index,
                        site_id=new_site_id,
                        constant=(
                            coupling.constant
                            + row["cell_offset"][
                                coordinate_index[coupling.coordinate]
                            ]
                        )
                        / repeats[coordinate_index[coupling.coordinate]],
                        factor=(
                            coupling.factor
                            / repeats[coordinate_index[coupling.coordinate]]
                        ),
                    )
                )
            site_couplings = []
            for coupling in atom.site_couplings:
                site_couplings.append(
                    dataclasses.replace(
                        coupling,
                        atom_index=new_index,
                        site_id=site_id,
                        factor=(
                            coupling.factor
                            / repeats[coordinate_index[coupling.coordinate]]
                        ),
                    )
                )
            atoms.append(
                dataclasses.replace(
                    atom,
                    atom_index=new_index,
                    site_id=site_id,
                    parent_fractional=parent_fractional,
                    surface_fractional=surface_fractional,
                    layer=layer_ids[(int(row["cell_offset"][2]), atom.layer)],
                    couplings=tuple(couplings),
                    site_couplings=tuple(site_couplings),
                )
            )

        metadata.surface_spec = dataclasses.replace(
            metadata.surface_spec,
            transform=metadata.surface_spec.transform @ np.diag(repeats),
            layer_origins=tuple(
                (self._layer_positions_for_cycle(self.layer_cycle.layers)[layer] + iz)
                / repeats[2]
                for iz in range(int(repeats[2]))
                for layer in self.layer_cycle.layers
            ),
        )
        metadata.sites = sites
        metadata.atoms = atoms
        return metadata

    @property
    def layer_behavior(self):
        """Return how layer selection is handled during crystal stacking."""
        return self.layer_behaviour

    @layer_behavior.setter
    def layer_behavior(self, behavior):
        self.layer_behaviour = behavior

    @property
    def height_absolute(self):
        H = (self.coherentDomainMatrix[-1][2, 3] + 1) * self.a[2]
        return H

    @property
    def end_layer_number(self):
        if self.layer_behavior == "select" and self.start_layer_number != -1:
            return self.start_layer_number
        idxmax = np.argmax(self.basis[:, 3])
        return self.basis[idxmax][7]

    @property
    def pos_absolute(self):
        return self.coherentDomainMatrix[0][2, 3] * self.a[2]

    @pos_absolute.setter
    def pos_absolute(self, pos):
        current_pos = self.pos_absolute
        for mat in self.coherentDomainMatrix:
            mat[2, 3] += (pos - current_pos) / self.a[2]

    @property
    def loc_absolute(self):
        return self.pos_absolute

    def stack_on(self, below_loc, below_height, below_layer=-1, below_state=None):
        """Place this unit cell on the object below it.

        :param float below_loc:
            Absolute reference location of the object below in Angstrom.
            Unit cells do not otherwise use this value.
        :param float below_height:
            Absolute top height of the object below in Angstrom.
        :param float below_layer:
            Top cyclic layer identifier of the object below.
        """
        if below_state is None:
            below_state = LayerState(self.layer_cycle, below_layer)
        if self.layer_behavior == "select":
            self._start_layer = resolve_upper_start(
                below_state, self.layer_cycle, self.layer_transition
            )
        else:
            self.start_layer_number = below_layer
        self.pos_absolute = below_height

    @property
    def layer_state(self):
        """Return the top structural-layer state of this unit cell."""
        return LayerState(self.layer_cycle, self.end_layer_number)

    @property
    def stacking_height_absolute(self):
        """Return the nominal height passed to the object above."""
        return self.height_absolute

    @property
    def stacking_loc_absolute(self):
        """Return the nominal reference location passed upward."""
        return self.loc_absolute

    @property
    def start_layer_number(self):
        return self._start_layer

    @start_layer_number.setter
    def start_layer_number(self, ln):
        if self.layer_behavior != "select":
            return
        if ln == -1:
            self._start_layer = self.layers[0]
            return
        matches = np.flatnonzero(self.layers == ln)
        if matches.size == 0:
            raise ValueError(f"Layer {ln} does not exist in UnitCell {self.name}.")
        self._start_layer = self.layers[(matches[0] + 1) % len(self.layers)]

    # in eV
    def setEnergy(self, E):
        """
        for dispersion and absorption correction
        in eV
        """
        self._E = E
        self.lookupScatteringFactors(E)

    def addFitParameter(self, indexarray, limits=(-np.inf, np.inf), **keyargs):
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

        """  # noqa: E501

        atoms, par = indexarray
        try:
            parindexes = np.array(
                [
                    p
                    if isinstance(p, np.integer | int)
                    else UnitCell.parameterLookup[p]
                    for p in np.atleast_1d(par)
                ],
                dtype=np.intp,
            )
        except Exception as e:
            raise ValueError(f"Invalid atom parameter name {par}.") from e
        atoms = np.atleast_1d(atoms).astype(np.intp)

        if "name" in keyargs:
            name = keyargs["name"]
        else:
            parameternames = tuple(
                [UnitCell.parameterLookup_inv[n] for n in parindexes]
            )
            atomnames = tuple([f"{n}_{self.names[n]}" for n in atoms])
            name = (
                self.name
                + " "
                + " ".join(
                    ["_".join((n, p)) for p, n in zip(parameternames, atomnames)]
                )
            )

        if name in self.fitparnames:
            raise ValueError(f"Absolute fit parameter {name} already exists.")

        indexarray = atoms, parindexes
        try:
            curr_val = self.basis[indexarray]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e

        if not limits[0] <= np.mean(curr_val) <= limits[1]:
            raise ValueError(
                f"start parameter is not within limits! Is not: {limits[0]} <= {np.mean(self.basis[indexarray])} <= {limits[1]}"  # noqa: E501
            )
        prior = keyargs.get("prior", None)
        keyargs["name"] = name
        keyargs["prior"] = prior

        par = Parameter(
            name, indexarray, ParameterType.ABSOLUTE, limits, prior, keyargs
        )
        self.parameters["absolute"].append(par)
        return par

    def addRelParameter(self, indexarray, factors, limits=(-np.inf, np.inf), **keyargs):
        # if not limits[0] <= np.mean(self.basis[indexarray]) <= limits[1]:
        #    raise ValueError("start parameter is not within limits!")

        atoms, par = indexarray
        try:
            parindexes = np.array(
                [
                    p
                    if isinstance(p, np.integer | int)
                    else UnitCell.parameterLookup[p]
                    for p in np.atleast_1d(par)
                ],
                dtype=np.intp,
            )
        except Exception as e:
            raise ValueError(f"Invalid atom parameter name {par}.") from e
        atoms = np.atleast_1d(atoms).astype(np.intp)

        if "name" in keyargs:
            name = keyargs["name"]
        else:
            parameternames = tuple(
                [UnitCell.parameterLookup_inv[n] + "_r" for n in parindexes]
            )
            atomnames = tuple([f"{n}_{self.names[n]}" for n in atoms])
            name = (
                self.name
                + " "
                + " ".join(
                    ["_".join((n, p)) for p, n in zip(parameternames, atomnames)]
                )
            )

        if name in self.fitparnames:
            raise ValueError(f"Relative fit parameter {name} already exists.")

        prior = keyargs.get("prior", None)
        keyargs["name"] = name
        keyargs["prior"] = prior

        factors = np.atleast_1d(factors)
        try:
            curr_val = self.basis[(atoms, parindexes)]
        except IndexError as e:
            raise ValueError(
                f"UnitCell {self.name}: Invalid parameter indices ({atoms}, {parindexes}) ."  # noqa: E501
            ) from e

        if curr_val.shape != factors.shape:
            raise ValueError(
                "Number of basis parameters does not match number of factors."
            )

        par = Parameter(
            name,
            (atoms, parindexes),
            ParameterType.RELATIVE,
            limits,
            prior,
            keyargs,
            factors,
        )
        self.parameters["relative"].append(par)
        return par

    def wyckoff_sites(self):
        """Return Wyckoff site metadata for this unit cell.

        :returns:
            List of site dictionaries. The list is empty when no symmetry
            metadata is attached.
        :rtype:
            list
        """
        if self.symmetry_metadata is None:
            return []
        return self.symmetry_metadata.wyckoff_sites(self.parameters)

    def wyckoff(self, site_id):
        """Return metadata for one Wyckoff site.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :returns:
            The matching site metadata dictionary.
        :rtype:
            dict
        :raises ValueError:
            If this unit cell has no symmetry metadata or the site is unknown.
        """
        return self._wyckoff_site(site_id)

    def set_wyckoff_atom_parameter(self, site_id, parameter, value):
        """Set one field on every generated atom in a Wyckoff site.

        Coordinates are fractional coordinates in the surface unit cell.
        Coordinate values may be a scalar, applied to every atom, or one value
        per atom in the order reported by
        ``wyckoff(site_id)["atom_indices"]``. The site-wide ``iDW``, ``oDW``,
        and ``occ`` fields require a scalar. Both the active basis and its
        unfitted baseline are updated.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :param str parameter:
            One of ``"x"``, ``"y"``, ``"z"``, ``"iDW"``, ``"oDW"``,
            or ``"occ"``.
        :param value:
            Scalar value or one value per generated atom.
        :raises ValueError:
            If the parameter is invalid or the value cannot be broadcast to
            the atoms in the site.
        """
        allowed = {"x", "y", "z", "iDW", "oDW", "occ"}
        if parameter not in allowed:
            raise ValueError(
                "parameter must be one of 'x', 'y', 'z', 'iDW', 'oDW', or 'occ'."
            )
        site = self.wyckoff(site_id)
        atoms = np.asarray(site["atom_indices"], dtype=np.intp)
        if not len(atoms):
            raise ValueError(f"Wyckoff site {site_id} has no generated atoms.")
        value_array = np.asarray(value, dtype=np.float64)
        if parameter in {"iDW", "oDW", "occ"} and value_array.ndim:
            raise ValueError(f"{parameter} must be a scalar site-wide value.")
        try:
            values = np.broadcast_to(
                value_array,
                atoms.shape,
            )
        except ValueError as error:
            raise ValueError(
                f"value must be scalar or contain one value for each of the "
                f"{len(atoms)} atoms in Wyckoff site {site_id}."
            ) from error

        column = self.parameterLookup[parameter]
        self.basis_0[atoms, column] = values
        self.basis[atoms, column] = values
        self._update_wyckoff_metadata_from_atom_values(
            site_id,
            parameter,
            atoms,
            values,
        )
        self._refresh_after_baseline_change()

    def set_wyckoff_site_parameter(self, site_id, parameter, value):
        """Set one representative coordinate or site-wide physical parameter.

        An absolute parent conventional fractional coordinate is propagated to
        all generated surface-cell atoms through their stored space-group and
        surface-transform couplings. ``iDW``, ``oDW``, and ``occ`` are applied
        as scalar values to every atom in the site. Both the active basis and
        its unfitted baseline are updated.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :param str parameter:
            Parent conventional fractional axis ``"x"``, ``"y"``, or
            ``"z"``; or site-wide ``"iDW"``, ``"oDW"``, or ``"occ"``.
        :param float value:
            New absolute coordinate in parent fractional units or scalar
            physical parameter value.
        :raises ValueError:
            If the parameter is invalid or required symmetry metadata is absent.
        """
        allowed = {"x", "y", "z", "iDW", "oDW", "occ"}
        if parameter not in allowed:
            raise ValueError(
                "parameter must be one of 'x', 'y', 'z', 'iDW', 'oDW', or 'occ'."
            )
        if parameter in {"iDW", "oDW", "occ"}:
            self.set_wyckoff_atom_parameter(site_id, parameter, value)
            return

        site = self.wyckoff(site_id)
        representative = site.get("representative_parent_fractional")
        if representative is None:
            raise ValueError(
                f"Wyckoff site {site_id} has no representative parent coordinate."
            )

        axis_index = {"x": 0, "y": 1, "z": 2}[parameter]
        value = float(value)
        delta = value - representative[axis_index]
        couplings = [
            coupling
            for coupling in self.wyckoff_site_couplings(site_id)
            if coupling.axis == parameter
        ]
        if not couplings:
            raise ValueError(
                f"No site-displacement couplings found for Wyckoff site "
                f"{site_id} and parent axis {parameter}."
            )

        surface_deltas = {}
        for coupling in couplings:
            column = self.parameterLookup[coupling.coordinate]
            change = delta * coupling.factor
            self.basis_0[coupling.atom_index, column] += change
            self.basis[coupling.atom_index, column] += change
            surface_deltas.setdefault(
                coupling.atom_index,
                np.zeros(3, dtype=np.float64),
            )[{"x": 0, "y": 1, "z": 2}[coupling.coordinate]] += change

        self._update_wyckoff_metadata_from_site_value(
            site_id,
            axis_index,
            value,
            surface_deltas,
        )
        self._refresh_after_baseline_change()

    def _refresh_after_baseline_change(self):
        self.errors = None
        self._basis_parvalues = None
        self._errors_parvalues = None

    def _update_wyckoff_metadata_from_atom_values(
        self,
        site_id,
        parameter,
        atoms,
        values,
    ):
        model = self.symmetry_metadata
        if parameter in {"iDW", "oDW", "occ"}:
            model.sites = tuple(
                dataclasses.replace(site, **{parameter: float(values[0])})
                if site.site_id == site_id
                else site
                for site in model.sites
            )
            return

        coordinate = {"x": 0, "y": 1, "z": 2}[parameter]
        values_by_atom = dict(zip(atoms.tolist(), values.tolist()))
        transform = model.surface_spec.transform
        updated_atoms = []
        for atom in model.atoms:
            if atom.atom_index not in values_by_atom:
                updated_atoms.append(atom)
                continue
            surface = np.array(atom.surface_fractional, copy=True)
            surface_delta = values_by_atom[atom.atom_index] - surface[coordinate]
            surface[coordinate] = values_by_atom[atom.atom_index]
            parent = np.array(atom.parent_fractional, copy=True)
            parent += transform[:, coordinate] * surface_delta
            updated_atoms.append(
                dataclasses.replace(
                    atom,
                    surface_fractional=surface,
                    parent_fractional=parent,
                )
            )
        model.atoms = updated_atoms

    def _update_wyckoff_metadata_from_site_value(
        self,
        site_id,
        axis_index,
        value,
        surface_deltas,
    ):
        model = self.symmetry_metadata
        model.sites = tuple(
            dataclasses.replace(
                site,
                representative_parent_fractional=tuple(
                    value if index == axis_index else coordinate
                    for index, coordinate in enumerate(
                        site.representative_parent_fractional
                    )
                ),
            )
            if site.site_id == site_id
            else site
            for site in model.sites
        )
        transform = model.surface_spec.transform
        model.atoms = [
            dataclasses.replace(
                atom,
                surface_fractional=np.asarray(atom.surface_fractional)
                + surface_deltas[atom.atom_index],
                parent_fractional=np.asarray(atom.parent_fractional)
                + transform @ surface_deltas[atom.atom_index],
            )
            if atom.atom_index in surface_deltas
            else atom
            for atom in model.atoms
        ]

    def wyckoff_couplings(self, site_id=None):
        """Return symmetry couplings for generated Wyckoff atom coordinates.

        :param str site_id:
            Optional site identifier. If omitted, all couplings are returned.
        :returns:
            List of coordinate coupling metadata objects.
        :rtype:
            list
        """
        if self.symmetry_metadata is None:
            return []
        return self.symmetry_metadata.wyckoff_couplings(site_id)

    def wyckoff_site_couplings(self, site_id=None):
        """Return couplings for representative Wyckoff-site displacement.

        :param str site_id:
            Optional site identifier. If omitted, all couplings are returned.
        :returns:
            List of site-displacement coupling metadata objects.
        :rtype:
            list
        """
        if self.symmetry_metadata is None:
            return []
        return self.symmetry_metadata.wyckoff_site_couplings(site_id)

    def atom_wyckoff_metadata(self, atom_index):
        """Return symmetry metadata for one atom.

        :param int atom_index:
            Index in ``basis``.
        :returns:
            Atom metadata or ``None`` when no metadata is available.
        :rtype:
            object or None
        """
        if self.symmetry_metadata is None:
            return None
        return self.symmetry_metadata.atom_wyckoff_metadata(atom_index)

    @staticmethod
    def _limit_for(limit_spec, parameter_name):
        if isinstance(limit_spec, dict):
            return limit_spec.get(parameter_name, (-np.inf, np.inf))
        return limit_spec

    @staticmethod
    def _delta_limits(limits, absolute_limits, reference_value):
        if absolute_limits is None:
            return limits
        if limits != (-np.inf, np.inf):
            raise ValueError("Use either limits or absolute_limits, not both.")
        return (
            absolute_limits[0] - reference_value,
            absolute_limits[1] - reference_value,
        )

    def _wyckoff_site(self, site_id):
        if self.symmetry_metadata is None:
            raise ValueError(f"UnitCell {self.name} has no symmetry metadata.")
        site = next(
            (
                site
                for site in self.wyckoff_sites()
                if site["site_id"] == site_id
            ),
            None,
        )
        if site is None:
            raise ValueError(f"Unknown Wyckoff site {site_id}.")
        return site

    @staticmethod
    def _coupling_parameter_arrays(couplings):
        atoms = np.asarray(
            [coupling.atom_index for coupling in couplings],
            dtype=np.intp,
        )
        coordinates = tuple(coupling.coordinate for coupling in couplings)
        factors = np.asarray(
            [coupling.factor for coupling in couplings],
            dtype=np.float64,
        )
        return atoms, coordinates, factors

    def _add_wyckoff_relative_parameter(
        self,
        site_id,
        selector_name,
        selector_value,
        kind,
        couplings,
        limits,
        default_name,
        empty_message,
        keyargs,
    ):
        if not couplings:
            raise ValueError(empty_message)
        atoms, coordinates, factors = self._coupling_parameter_arrays(couplings)
        keyargs.setdefault("name", default_name)
        settings = keyargs.setdefault("wyckoff", {})
        settings.update(
            {
                "site_id": site_id,
                "kind": kind,
                selector_name: selector_value,
                "value_kind": "delta",
            }
        )
        return self.addRelParameter(
            (atoms, coordinates),
            factors,
            limits=limits,
            **keyargs,
        )

    def addWyckoffParameter(
        self,
        site_id,
        variable,
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **keyargs,
    ):
        """Add an absolute symmetry-preserving Wyckoff variable parameter.

        The exposed fit value is the absolute Wyckoff variable. Internally,
        its difference from the site's reference value is propagated through
        every affine coupling, preserving symmetry-equivalent coordinates.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :param str variable:
            Wyckoff variable name, for example ``"u"``.
        :param tuple limits:
            Absolute variable limits in parent fractional units.
        :param tuple absolute_limits:
            Deprecated alias for ``limits`` retained for compatibility.
        :returns:
            Parameter exposing the absolute Wyckoff variable value.
        :rtype:
            CTRutil.Parameter
        :raises ValueError:
            If no matching affine couplings exist.
        """
        site = self._wyckoff_site(site_id)
        if not site["variables"]:
            raise ValueError(
                f"Wyckoff site {site_id} has no positional coordinate variables."
            )
        if variable not in site["variables"]:
            raise ValueError(f"Wyckoff site {site_id} has no variable {variable}.")

        if absolute_limits is not None:
            if limits != (-np.inf, np.inf):
                raise ValueError("Use either limits or absolute_limits, not both.")
            limits = absolute_limits
        reference_value = site["variables"][variable]
        couplings = [
            coupling
            for coupling in self.wyckoff_couplings(site_id)
            if coupling.variable == variable
        ]
        parameter = self._add_wyckoff_relative_parameter(
            site_id,
            "variable",
            variable,
            "coordinate",
            couplings,
            limits,
            f"{self.name} {site_id}_{variable}_wyckoff_coordinate",
            (
                f"No affine couplings found for Wyckoff site {site_id} "
                f"and variable {variable}."
            ),
            keyargs,
        )
        parameter.settings["wyckoff"].update(
            {
                "value_kind": "absolute",
                "reference_value": reference_value,
            }
        )
        return parameter

    def addWyckoffParameters(
        self,
        site_id,
        variables=None,
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **keyargs,
    ):
        """Add absolute symmetry-preserving Wyckoff variable parameters.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :param variables:
            Iterable of coordinate variables to fit. If ``None``, all free
            variables on the site are fitted.
        :type variables:
            iterable or None
        :param limits:
            One absolute ``(lower, upper)`` tuple for every variable or a
            dictionary mapping variable names to absolute limits.
        :param absolute_limits:
            Deprecated alias for absolute ``limits`` retained for compatibility.
        :returns:
            Created absolute variable parameters in variable order.
        :rtype:
            list
        :raises ValueError:
            If the site has no free coordinate variables.
        """
        site_variables = tuple(self._wyckoff_site(site_id)["variables"])
        if not site_variables:
            raise ValueError(
                f"Wyckoff site {site_id} has no positional coordinate variables."
            )
        if variables is None:
            variables = site_variables
        variables = tuple(variables)

        parameters = []
        for variable in variables:
            parameter_keyargs = copy.deepcopy(keyargs)
            parameter_limits = self._limit_for(limits, variable)
            parameter_absolute_limits = (
                None
                if absolute_limits is None
                else self._limit_for(absolute_limits, variable)
            )
            parameters.append(
                self.addWyckoffParameter(
                    site_id,
                    variable,
                    limits=parameter_limits,
                    absolute_limits=parameter_absolute_limits,
                    **parameter_keyargs,
                )
            )
        return parameters

    def addWyckoffShift(
        self,
        site_id,
        axis,
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **keyargs,
    ):
        """Add a relative shift from a symmetry site along a parent direction.

        ``axis`` is a parent conventional-cell direction. The exposed value is
        a relative displacement from the representative symmetry-site
        coordinate, not an absolute parent coordinate. Generated atoms move
        through the stored space-group operation and surface-cell transform.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :param str axis:
            Parent representative coordinate, one of ``"x"``, ``"y"``, or
            ``"z"``.
        :param tuple limits:
            Fit limits for the coordinate change in parent fractional units.
        :param tuple absolute_limits:
            Optional absolute parent-coordinate bounds converted to relative
            shift bounds. The fitted value remains a relative displacement.
        :returns:
            Created relative fit parameter.
        :rtype:
            CTRutil.Parameter
        :raises ValueError:
            If no matching site-displacement couplings exist.
        """
        if axis not in {"x", "y", "z"}:
            raise ValueError("axis must be one of 'x', 'y', or 'z'.")
        site = self._wyckoff_site(site_id)
        representative = site.get("representative_parent_fractional")
        if representative is None:
            raise ValueError(
                f"Wyckoff site {site_id} has no representative parent coordinate."
            )
        axis_index = {"x": 0, "y": 1, "z": 2}[axis]
        limits = self._delta_limits(
            limits,
            absolute_limits,
            representative[axis_index],
        )

        couplings = [
            coupling
            for coupling in self.wyckoff_site_couplings(site_id)
            if coupling.axis == axis
        ]
        return self._add_wyckoff_relative_parameter(
            site_id,
            "axis",
            axis,
            "site_displacement",
            couplings,
            limits,
            f"{self.name} {site_id}_{axis}_wyckoff_site",
            (
                f"No site-displacement couplings found for Wyckoff site "
                f"{site_id} and parent axis {axis}."
            ),
            keyargs,
        )

    def addWyckoffShifts(
        self,
        site_id,
        axes=("x", "y", "z"),
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **keyargs,
    ):
        """Add relative symmetry-site shifts along parent-cell directions.

        Every exposed value is a displacement from the representative
        symmetry-site coordinate along the selected parent conventional-cell
        direction; it is not an absolute coordinate.

        :param str site_id:
            Site identifier returned by :meth:`wyckoff_sites`.
        :param axes:
            Parent representative coordinate axes to fit.
        :type axes:
            iterable
        :param limits:
            Either one ``(lower, upper)`` tuple applied to every axis or a
            dictionary mapping axis names to delta limits.
        :param absolute_limits:
            Optional absolute parent-coordinate bounds converted to relative
            shift bounds. The fitted values remain relative displacements.
        :returns:
            Created relative fit parameters in axis order.
        :rtype:
            list
        """

        parameters = []
        for axis in axes:
            parameter_keyargs = copy.deepcopy(keyargs)
            parameter_limits = self._limit_for(limits, axis)
            parameter_absolute_limits = (
                None
                if absolute_limits is None
                else self._limit_for(absolute_limits, axis)
            )
            parameters.append(
                self.addWyckoffShift(
                    site_id,
                    axis,
                    limits=parameter_limits,
                    absolute_limits=parameter_absolute_limits,
                    **parameter_keyargs,
                )
            )
        return parameters

    def showFitparameters(self):
        print(self.fitparameterList())

    def parameter_list(self):
        return self.parameters["absolute"] + self.parameters["relative"]

    @property
    def fitparnames(self):
        # pars = [self.fitparToStr(i,True) for i in range(len(self.fitparameters))] + [self.fitparToStr(i + len(self.fitparameters),True) for i in range(len(self.relfitparam))]  # noqa: E501
        # pars =
        return [p.name for p in self.parameters["absolute"]] + [
            p.name for p in self.parameters["relative"]
        ]

    @property
    def priors(self):
        priorlist = []
        for par in self.parameters["absolute"] + self.parameters["relative"]:
            if par.prior is None:
                if tuple(par.limits) == (-np.inf, np.inf):
                    raise Exception(
                        f"Infinite range given as prior for parameter {par.name}. You must provide an explicit prior."  # noqa: E501
                    )
                else:
                    priorlist.append(
                        tuple(par.limits)
                    )  # has to be converted to correct distribution in CTROptimizer
            else:
                priorlist.append(par.prior)  # real prior distribution
        return priorlist

    def fitparameterList(self, verbosity=3):
        if verbosity > 1:
            outstr = "{}".format("# Direct fit parameters:\n")
        else:
            outstr = ""

        if len(self.parameters["absolute"]) > 0:
            outstr += "{:10} {:15} {:20} {:10} {:10}".format(
                "Id", "atoms", "parameters", "Value", "Limits"
            )

            for par in self.parameters["absolute"]:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No direct fit parameters"

        outstr += "\n{}".format("# Relative fit parameters:\n")
        if len(self.parameters["relative"]) > 0:
            outstr += "{:10} {:15} {:20} {:10} {:10} {:10}".format(
                "Id", "atoms", "parameters", "Value", "Factors", "Limits"
            )

            for par in self.parameters["relative"]:
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
            # if isinstance(par.indices[1], (np.integer,int)):
            #    parameternames = (UnitCell.parameterLookup_inv[par.indices[1]],)
            #    atoms = (f"{par.indices[0]}_{self.names[par.indices[0]]}",)
            # else:
            parameternames = tuple(
                [UnitCell.parameterLookup_inv[n] for n in par.indices[1]]
            )
            atoms = tuple([f"{n}_{self.names[n]}" for n in par.indices[0]])
            return f"{par.name:10} {str(atoms):15} {str(parameternames):20} {val:10} {str(par.limits):10}"  # noqa: E501

        elif par.kind & ParameterType.RELATIVE:
            val = np.mean(
                (self.basis[par.indices] - basis_0[par.indices]) / par.factors
            )
            wyckoff = par.settings.get("wyckoff", {})
            if wyckoff.get("value_kind") == "absolute":
                val += wyckoff["reference_value"]
            # if isinstance(par.indices[1], (np.integer,int)):
            #    parameternames = (UnitCell.parameterLookup_inv[par.indices[1]] + '_r',)
            #    atoms = (f"{par.indices[0]}_{self.names[par.indices[0]]}",)
            # else:
            parameternames = tuple(
                [UnitCell.parameterLookup_inv[n] + "_r" for n in par.indices[1]]
            )
            atoms = tuple([f"{n}_{self.names[n]}" for n in par.indices[0]])
            return f"{par.name:10} {str(atoms):15} {str(parameternames):20} {val:10} {str(par.factors):10} {str(par.limits):10}"  # noqa: E501

        else:
            raise ValueError(f"Unvalid parameter type {par.kind} for UnitCell")

    # def fixInitialBasis(self):
    #    self.basis_0 = np.copy(self.basis)

    def getInitialParameters(self, force_recalculate=False):
        x0, _, _ = self.getStartParamAndLimits(force_recalculate)
        return x0

    def getStartParamAndLimits(self, force_recalculate=False):
        # if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        try:
            mismatch = not np.allclose(
                self._basis_parvalues, self.basis, equal_nan=True
            )
        except Exception:
            mismatch = True
        recalculate = force_recalculate or mismatch
        abspar = self.parameters["absolute"]
        relpar = self.parameters["relative"]

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
                value = np.mean(
                    (self.basis[par.indices] - self.basis_0[par.indices])
                    / par.factors
                )
                wyckoff = par.settings.get("wyckoff", {})
                if wyckoff.get("value_kind") == "absolute":
                    value += wyckoff["reference_value"]
                x0.append(value)
            else:
                x0.append(par.value)
            lower.append(par.limits[0])
            upper.append(par.limits[1])

        return np.array(x0), np.array(lower), np.array(upper)

    def setFitParameters(self, x):
        x_0 = x[: len(self.parameters["absolute"])]
        x_r = x[len(self.parameters["absolute"]) :]
        self.basis[:] = self.basis_0
        for val, par in zip(x_0, self.parameters["absolute"]):
            self.basis[par.indices] = val
            par.value = val
        for val, par in zip(x_r, self.parameters["relative"]):
            internal_value = self._relative_internal_value(par, val)
            self.basis[par.indices] += par.factors * internal_value
            par.value = val
        self._basis_parvalues = np.copy(self.basis)

    def setLimits(self, lim):
        x_0 = lim[: len(self.parameters["absolute"])]
        x_r = lim[len(self.parameters["absolute"]) :]
        for val, par in zip(x_0, self.parameters["absolute"]):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0], val[1])
        for val, par in zip(x_r, self.parameters["relative"]):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0], val[1])

    def setFitErrors(self, errors):
        self.errors = np.full_like(self.basis, np.nan)
        err_0 = errors[: len(self.parameters["absolute"])]
        err_r = errors[len(self.parameters["absolute"]) :]
        for val, par in zip(err_0, self.parameters["absolute"]):
            self.errors[par.indices] = val
            par.error = val
        for val, par in zip(err_r, self.parameters["relative"]):
            self.errors[par.indices] = np.nan_to_num(self.errors[par.indices], nan=0.0)
            self.errors[par.indices] += np.abs(par.factors) * val
            par.error = val
        self._errors_parvalues = np.copy(self.errors)

    def getFitErrors(self):
        if self.errors is None:
            raise ValueError("No errors have been set.")
        try:
            mismatch = not np.allclose(
                self._errors_parvalues, self.errors, equal_nan=True
            )
        except Exception:
            mismatch = True
        abspar = self.parameters["absolute"]
        relpar = self.parameters["relative"]
        err0 = []
        for par in abspar:
            if mismatch:
                err0.append(np.mean(self.errors[par.indices]))
            else:
                err0.append(par.error)

        for par in relpar:
            if mismatch:
                err0.append(np.mean(self.errors[par.indices] / np.abs(par.factors)))
            else:
                err0.append(par.error)

        err0 = np.array(err0)
        if np.any(~np.isfinite(err0)):
            raise ValueError("Some errors are non-finite or not set.")
        return err0

    def build_selected_basis(self):
        if self.layer_behaviour == "select":
            if self.start_layer_number == -1.0:
                warnings.warn(
                    "Layer behaviour is >select<, but start number is -1 "
                    "(undefined). Proceed with the first layer."
                )
                ln = self.layers[0]
            else:
                ln = self.start_layer_number

            mask_layer = self.basis[:, 7] == ln
            if not np.any(mask_layer):
                raise ValueError(f"Layer {ln} does not exist in UnitCell.")
            return (
                self.basis[mask_layer],
                self.f[mask_layer],
                np.asarray(self.names)[mask_layer],
            )

        elif self.layer_behaviour == "ignore":
            return self.basis, self.f, self.names
        else:
            warnings.warn(
                f"unknown layer behaviour {self.layer_behaviour}: proceed calculating F while ignoring layer behaviour."  # noqa: E501
            )
            return self.basis, self.f, self.names

    # returns the structure factor of the unit cell
    # h,k,l have to be 1d arrays
    def F_uc_bulk(self, h, k, l, atten=0):  # noqa: E741
        """Return one attenuated bulk unit-cell amplitude in electrons.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :param float atten:
            Dimensionless attenuation exponent per unit-cell translation.
        :returns:
            Complex structure-factor amplitude in electrons.
        :rtype: numpy.ndarray
        """
        basis, formf, names = self.build_selected_basis()
        if ctr_accel_enabled():
            h, k, l = _ensure_contiguous(h, k, l, testOnly=False, astype=np.float64)  # noqa: E741
            accel = _ctr_accel_module()
            F = accel.unitcell_F_uc_bulk(
                h,
                k,
                l,
                atten,
                basis,
                formf,
                self.refHKLTransform,
                self.B_mat,
                self.R_mat,
                self.R_mat_inv,
                np.asarray(self.coherentDomainMatrix),
                np.asarray(self.coherentDomainOccupancy),
                self.uc_area,
            )
            return F
        else:
            F = np.zeros(h.size, dtype=np.complex128)
            f = np.zeros(h.size, dtype=np.complex128)
            hkl = self.refHKLTransform @ np.vstack((h, k, l))
            Q_cart2 = (self.B_mat @ hkl) ** 2
            Q_para2 = np.sum(Q_cart2[:2], axis=0)  # squared!!!
            Q_perp2 = Q_cart2[2]  # squared!!!
            Q2 = Q_para2 + Q_perp2  # squared!!!

            domainmatrix = [
                self.R_mat_inv @ mat[:, :-1] @ self.R_mat
                for mat in self.coherentDomainMatrix
            ]

            for i in range(len(basis)):
                f[:] = formf[i][10] + formf[i][11] + 1j * formf[i][12]
                for j in range(5):
                    f += formf[i][j] * np.exp(-formf[i][j + 5] * Q2)
                f *= np.exp(
                    -(basis[i][4] * Q_para2 + basis[i][5] * Q_perp2) / (16 * np.pi**2)
                )
                f *= basis[i][6]
                for mat, weight, eff_mat in zip(
                    self.coherentDomainMatrix,
                    self.coherentDomainOccupancy,
                    domainmatrix,
                ):
                    xyz_rel = np.dot(eff_mat, basis[i][1:4]) + mat[:, -1]
                    F += (
                        weight
                        * f
                        * np.exp(2j * np.pi * np.sum(hkl.T * xyz_rel, axis=1))
                        * math.exp(atten * xyz_rel[2])
                    )
            return F

    # returns the structure factor of the unit cell
    # h,k,l have to be 1d arrays
    def F_uc_bulk_direct(self, h, k, l, atten=0):  # noqa: E741
        """Return one bulk-cell amplitude without reference conversion.

        The result is unnormalized and has units of electrons.

        :param numpy.ndarray h:
            Reciprocal coordinate in this unit cell's r.l.u.
        :param numpy.ndarray k:
            Reciprocal coordinate in this unit cell's r.l.u.
        :param numpy.ndarray l:
            Reciprocal coordinate in this unit cell's r.l.u.
        :param float atten:
            Dimensionless attenuation exponent.
        :returns:
            Complex structure-factor amplitude in electrons.
        :rtype: numpy.ndarray
        """
        basis, formf, names = self.build_selected_basis()
        F = np.zeros(h.size, dtype=np.complex128)
        f = np.zeros(h.size, dtype=np.complex128)
        hkl = np.vstack((h, k, l))
        Q_cart2 = (self.B_mat @ hkl) ** 2
        Q_para2 = np.sum(Q_cart2[:2], axis=0)  # squared!!!
        Q_perp2 = Q_cart2[2]  # squared!!!
        Q2 = Q_para2 + Q_perp2  # squared!!!

        domainmatrix = [
            self.R_mat_inv @ mat[:, :-1] @ self.R_mat
            for mat in self.coherentDomainMatrix
        ]
        for i in range(len(basis)):
            f[:] = formf[i][10] + formf[i][11] + 1j * formf[i][12]
            for j in range(5):
                f += formf[i][j] * np.exp(-formf[i][j + 5] * Q2)
            f *= np.exp(
                -(basis[i][4] * Q_para2 + basis[i][5] * Q_perp2) / (16 * np.pi**2)
            )
            f *= basis[i][6]
            for mat, weight, eff_mat in zip(
                self.coherentDomainMatrix, self.coherentDomainOccupancy, domainmatrix
            ):
                xyz_rel = np.dot(eff_mat, basis[i][1:4]) + mat[:, -1]
                F += (
                    weight
                    * f
                    * np.exp(2j * np.pi * np.sum(hkl.T * xyz_rel, axis=1))
                    * math.exp(atten * xyz_rel[2])
                )
        return F

    def F_uc(self, h, k, l):  # noqa: E741
        """Return the canonical unit-cell structure factor in electrons.

        No unit-cell area or volume normalization is applied. Input reciprocal
        coordinates are interpreted in the configured reference unit cell and
        transformed into this unit cell before evaluation.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex structure-factor amplitude in electrons.
        :rtype: numpy.ndarray
        """
        basis, formf, names = self.build_selected_basis()
        if ctr_accel_enabled() and not self._special_formfactors_present:
            h, k, l = _ensure_contiguous(h, k, l, testOnly=False, astype=np.float64)  # noqa: E741
            accel = _ctr_accel_module()
            F = accel.unitcell_F_uc(
                h,
                k,
                l,
                basis,
                formf,
                self.refHKLTransform,
                self.B_mat,
                self.R_mat,
                self.R_mat_inv,
                np.asarray(self.coherentDomainMatrix),
                np.asarray(self.coherentDomainOccupancy),
                self.uc_area,
            )
            return F
        else:
            F = np.zeros(h.size, dtype=np.complex128)
            f = np.zeros(h.size, dtype=np.complex128)
            hkl = self.refHKLTransform @ np.vstack((h, k, l))
            Q_cart2 = (self.B_mat @ hkl) ** 2
            Q_para2 = np.sum(Q_cart2[:2], axis=0)  # squared!!!
            Q_perp2 = Q_cart2[2]  # squared!!!
            Q2 = Q_para2 + Q_perp2  # squared!!!

            domainmatrix = [
                self.R_mat_inv @ mat[:, :-1] @ self.R_mat
                for mat in self.coherentDomainMatrix
            ]
            for i, name in zip(range(len(basis)), names):
                if name in UnitCell.special_formfactors:
                    f[:] = (
                        UnitCell.special_formfactors[name][0](np.sqrt(Q2))
                        + formf[i][11]
                        + 1j * formf[i][12]
                    )
                else:
                    f[:] = formf[i][10] + formf[i][11] + 1j * formf[i][12]
                    for j in range(5):
                        f += formf[i][j] * np.exp(-formf[i][j + 5] * Q2)
                f *= np.exp(
                    -(basis[i][4] * Q_para2 + basis[i][5] * Q_perp2) / (16 * np.pi**2)
                )
                f *= basis[i][6]
                for mat, weight, eff_mat in zip(
                    self.coherentDomainMatrix,
                    self.coherentDomainOccupancy,
                    domainmatrix,
                ):
                    xyz_rel = np.dot(eff_mat, basis[i][1:4]) + mat[:, -1]
                    F += (
                        weight
                        * f
                        * np.exp(2j * np.pi * np.sum(hkl.T * xyz_rel, axis=1))
                    )
            return F

    def F_bulk(self, h, k, l, atten=0):  # noqa: E741
        """Return the semi-infinite bulk structure factor in electrons.

        The amplitude represents one lateral bulk unit cell. The geometric
        lattice sum is applied only along the out-of-plane direction.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :param float atten:
            Dimensionless attenuation exponent per bulk unit cell.
        :returns:
            Complex bulk amplitude in electrons per lateral bulk cell.
        :rtype: numpy.ndarray
        """
        basis, formf, names = self.build_selected_basis()
        if ctr_accel_enabled():
            h, k, l = _ensure_contiguous(h, k, l, testOnly=False, astype=np.float64)  # noqa: E741
            accel = _ctr_accel_module()
            F = accel.unitcell_F_bulk(
                h,
                k,
                l,
                atten,
                basis,
                formf,
                self.refHKLTransform,
                self.B_mat,
                self.R_mat,
                self.R_mat_inv,
                np.asarray(self.coherentDomainMatrix),
                np.asarray(self.coherentDomainOccupancy),
                self.uc_area,
            )
            return F
        else:
            hkl = self.refHKLTransform @ np.vstack((h, k, l))
            Fuc = self.F_uc_bulk_direct(*hkl, atten)
            return Fuc / (1 - np.exp(-2j * np.pi * l - atten))

    SQRT2pi = np.sqrt(2 * np.pi)

    def zDensity_G(self, z, h, k):
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
        basis, formf, names = self.build_selected_basis()

        hkl = (self.refHKLTransform @ np.array([h, k, 0.0])).flatten()
        Qpara2 = np.sum(np.dot(self.B_mat, hkl) ** 2)
        a, alpha, _, _ = self.getLatticeParameters()
        # uc_area = a[0]*a[1]*np.sin(alpha[2])
        # dispersion = self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
        zstep = np.diff(z)
        zstep_mean = np.mean(zstep)
        if not np.allclose(zstep, zstep_mean):
            warnings.warn(
                "zDensity: z stepsize is not equal in given z array."
                "This will result in numerical errors in electron density calculation!"
            )

        rho = np.zeros_like(z, dtype=np.complex128)
        rho_i = np.empty_like(z, dtype=np.complex128)

        domainmatrix = [
            self.R_mat_inv @ mat[:, :-1] @ self.R_mat
            for mat in self.coherentDomainMatrix
        ]

        for i, name in zip(range(len(basis)), names):
            deltaZ2i = basis[i][5] / (8 * (np.pi) ** 2)
            deltaPara2i = basis[i][4] / (8 * (np.pi) ** 2)
            for mat, weight, eff_mat in zip(
                self.coherentDomainMatrix, self.coherentDomainOccupancy, domainmatrix
            ):
                xyz_rel = eff_mat @ basis[i][1:4] + mat[:, -1]
                z_i = xyz_rel[2] * self._a[2]
                y_i_frac = xyz_rel[1]
                x_i_frac = xyz_rel[0]
                if name in UnitCell.special_eDensity:
                    rho_i[:] = (formf[i][11] + 1j * formf[i][12]) / (
                        np.sqrt(2 * np.pi * deltaZ2i)
                    )
                    rho_i *= np.exp(
                        -0.5 * (deltaPara2i * Qpara2 + ((z - z_i) ** 2) / deltaZ2i)
                    )
                    rho_i += gaussian_filter1d(
                        UnitCell.special_eDensity[name](z - z_i)
                        * np.exp(-0.5 * (deltaPara2i * Qpara2)),
                        np.sqrt(deltaZ2i) / zstep_mean,
                    )
                else:
                    rho_i[:] = (formf[i][10] + formf[i][11] + 1j * formf[i][12]) / (
                        np.sqrt(2 * np.pi * deltaZ2i)
                    )
                    rho_i *= np.exp(
                        -0.5 * (deltaPara2i * Qpara2 + ((z - z_i) ** 2) / deltaZ2i)
                    )
                    for j in range(5):
                        exp_dpara = formf[i][j + 5] + 0.5 * deltaPara2i
                        exp_dz = formf[i][j + 5] + 0.5 * deltaZ2i
                        rho_i += (
                            formf[i][j] / (np.sqrt(4.0 * np.pi * exp_dz))
                        ) * np.exp(
                            -exp_dpara * Qpara2 - (((z - z_i) ** 2) / (4 * exp_dz))
                        )

                rho_i *= np.exp(-2j * np.pi * (h * x_i_frac + k * y_i_frac))
                rho += rho_i * basis[i][6] * weight

        return rho / self.uc_area

    def zDensity_G_asbulk(self, z, h, k, noUC=30):
        basis, formf, names = self.build_selected_basis()

        hkl = (self.refHKLTransform @ np.array([h, k, 0.0])).flatten()
        Qpara2 = np.sum(np.dot(self.B_mat, hkl) ** 2)
        a, alpha, _, _ = self.getLatticeParameters()
        # uc_area = a[0]*a[1]*np.sin(alpha[2])
        # dispersion = self.f[i][10] + self.f[i][11] + 1j*self.f[i][12]
        rho = np.zeros_like(z, dtype=np.complex128)
        rho_i = np.empty_like(z, dtype=np.complex128)
        for no in range(noUC):
            noA = no * self._a[2]
            for i in range(len(basis)):
                deltaZ2i = basis[i][5] / (8 * (np.pi) ** 2)
                deltaPara2i = basis[i][4] / (8 * (np.pi) ** 2)
                z_i = basis[i][3] * self._a[2]
                y_i_frac = basis[i][2]
                x_i_frac = basis[i][1]

                rho_i[:] = (formf[i][10] + formf[i][11] + 1j * formf[i][12]) / (
                    np.sqrt(2 * np.pi * deltaZ2i)
                )
                rho_i *= np.exp(
                    -0.5 * (deltaPara2i * Qpara2 + ((z - z_i + noA) ** 2) / deltaZ2i)
                )

                for j in range(5):
                    exp_dpara = formf[i][j + 5] + 0.5 * deltaPara2i
                    exp_dz = formf[i][j + 5] + 0.5 * deltaZ2i
                    rho_i += (formf[i][j] / (np.sqrt(4.0 * np.pi * exp_dz))) * np.exp(
                        -exp_dpara * Qpara2 - (((z - z_i + noA) ** 2) / (4 * exp_dz))
                    )

                rho_i *= np.exp(-2j * np.pi * (h * x_i_frac + k * y_i_frac))
                rho += rho_i * basis[i][6]

        return rho / self.uc_area

    def lookupScatteringFactors(self, E):
        if hasattr(self, "f"):
            if self.f.shape != (self.basis.shape[0], 13):
                self.f = np.empty((self.basis.shape[0], 13), dtype=np.float64)
        else:
            self.f = np.empty((self.basis.shape[0], 13), dtype=np.float64)

        for i, name in enumerate(self.names):
            if name in UnitCell.special_formfactors:
                self.f[i, :11] = 0.0
                self.f[i, 11:] = UnitCell.special_formfactors[name][1](E)
            else:
                self.f[i, :11] = readWaasmaier(name)
                self.f[i, 11:] = readDispersion(name, E)
        return

    def plot3d(
        self,
        ucx=1,
        ucy=1,
        ucz=1,
        dwon=False,
        occuon=False,
        figure=None,
        translate=np.array([0.0, 0.0, 0.0]),
        domain=0,
        **keyargs,
    ):
        try:
            from mayavi import mlab
        except ImportError:
            warnings.warn("can not import mayavi: 3D plotting not supported")
            return

        if figure is None:
            figure = mlab.figure()
        if keyargs.get("useSelected", True):
            basis, formf, names = self.build_selected_basis()
        else:
            basis, _formf, names = self.basis, self.f, self.names
        if ucx == 0 or ucy == 0 or ucz == 0:
            raise ValueError(
                "One of ucx,ucy or ucz is zero. Must plot at least one unit cell."
            )

        def color_generator(name):
            return special_elementcolors.get(
                name, tuple(rgb_array[atomic_number(name) - 1])
            )

        elcolors = keyargs.get("color")
        if elcolors is None:
            elcolors = [color_generator(name) for name in names]
        elif isinstance(elcolors, tuple):
            elcolors = [elcolors for params in basis]

        resolution = keyargs.get("resolution")
        if resolution is None:
            resolution = 25

        mat = self.coherentDomainMatrix[domain]
        domainmatrix = self.R_mat_inv @ mat[:, :-1] @ self.R_mat
        for i, params in enumerate(basis):
            radius = cov_radii_array[atomic_number(names[i]) - 1][2] * 2
            # elcolor_c = keyargs.get('color')
            # if elcolor_c is None:
            #    elcolor_c = elements.rgb(int(params[0]))
            if occuon:
                occup = params[6]
            else:
                occup = 1
            x, y, z = params[1:4]
            positions = np.empty((np.abs(ucx * ucy * ucz), 3))
            iDW, oDW = DWtoDisorder(params[4:6])
            signx = np.sign(ucx)
            signy = np.sign(ucy)
            signz = np.sign(ucz)
            no = 0
            for xno in range(1, abs(ucx) + 1):
                for yno in range(1, abs(ucy) + 1):
                    for zno in range(1, abs(ucz) + 1):
                        if len(translate.shape) > 1:
                            x_t, y_t, z_t = (
                                np.array([x, y, z]) + translate[:, -1].flatten()
                            )
                            positions[no] = np.array(
                                [
                                    x_t + signx * (xno - 1),
                                    y_t + signy * (yno - 1),
                                    z_t + signz * (zno - 1),
                                ]
                            )
                        else:
                            positions[no] = (
                                np.array(
                                    [
                                        x + signx * (xno - 1),
                                        y + signy * (yno - 1),
                                        z + signz * (zno - 1),
                                    ]
                                )
                                + translate
                            )
                        if random.random() > occup:
                            positions[no] = np.nan
                        no += 1
            position_cart = self.R_mat @ (
                domainmatrix @ positions.T + np.atleast_2d(mat[:, -1]).T
            )
            # position_cart = self.RealspaceMatrix @ positions.T
            if len(translate.shape) > 1:
                position_cart = np.asarray(translate[:, :-1] @ position_cart)

            sigmas = np.empty_like(position_cart)
            sigmas[:2] = iDW
            sigmas[2] = oDW

            if dwon:
                position_cart = np.random.default_rng().normal(position_cart, sigmas)
            mlab.points3d(
                *position_cart,
                scale_factor=radius,
                color=elcolors[i],
                resolution=resolution,
                figure=figure,
            )

        atomlist = keyargs.get("atomlist")
        if atomlist is not None:
            atomlist.append(self.pos_cart_all(ucx, ucy, ucz, translate))

        return figure

    def pos_cart_all(
        self, ucx=1, ucy=1, ucz=1, translate=np.array([0.0, 0.0, 0.0]), domain=0
    ):
        dt = np.dtype(
            [("name", "U6"), ("x", np.float64), ("y", np.float64), ("z", np.float64)]
        )
        xyz_array = np.empty(np.abs(ucx * ucy * ucz) * len(self.names), dtype=dt)
        position_cart = []
        names = []
        mat = self.coherentDomainMatrix[domain]
        domainmatrix = self.R_mat_inv @ mat[:, :-1] @ self.R_mat

        for i, params in enumerate(self.basis):
            x, y, z = params[1:4]
            positions = np.empty((np.abs(ucx * ucy * ucz), 3))
            signx = np.sign(ucx)
            signy = np.sign(ucy)
            signz = np.sign(ucz)
            no = 0
            for xno in range(1, abs(ucx) + 1):
                for yno in range(1, abs(ucy) + 1):
                    for zno in range(1, abs(ucz) + 1):
                        if len(translate.shape) > 1:
                            x_t, y_t, z_t = (
                                np.array([x, y, z]) + translate[:, -1].flatten()
                            )
                            positions[no] = np.array(
                                [
                                    x_t + signx * (xno - 1),
                                    y_t + signy * (yno - 1),
                                    z_t + signz * (zno - 1),
                                ]
                            )
                        else:
                            positions[no] = (
                                np.array(
                                    [
                                        x + signx * (xno - 1),
                                        y + signy * (yno - 1),
                                        z + signz * (zno - 1),
                                    ]
                                )
                                + translate
                            )
                        no += 1
            positions_c = self.R_mat @ (
                domainmatrix @ positions.T + np.atleast_2d(mat[:, -1]).T
            )

            if len(translate.shape) > 1:
                positions_c = np.asarray(translate[:, :-1] @ positions_c)
            position_cart.append(positions_c)
            names.append([self.names[i] for k in range(no)])

        positions_all = np.concatenate(position_cart, axis=1)
        xyz_array["x"] = positions_all[0]
        xyz_array["y"] = positions_all[1]
        xyz_array["z"] = positions_all[2]
        xyz_array["name"] = np.concatenate(names)
        return xyz_array

    def pos_cart(self, atomNo, domain=0):
        mat = self.coherentDomainMatrix[domain]
        domainmatrix = self.R_mat_inv @ mat[:, :-1] @ self.R_mat
        xyz_rel = domainmatrix @ self.basis[atomNo][1:4] + mat[:, -1]
        return self.R_mat @ xyz_rel

    def pos_cart_error(self, atomNo, domain=0):
        xyz = self.pos_cart(atomNo, domain)
        realMatrix = self.RealspaceMatrix**2
        error_zeros = np.nan_to_num(self.errors, 0)
        xyz_error = np.sqrt(realMatrix @ (error_zeros[atomNo][1:4] ** 2).T)
        return xyz, xyz_error

    def distanceVectorAtom(self, atomNo1, atomNo2, domain=0):
        return self.pos_cart(atomNo1, domain) - self.pos_cart(atomNo2, domain)

    def bondingLength(self, xyz_atom, domain=0):
        pos_atom_cart = self.RealspaceMatrix @ xyz_atom.T
        blength = np.empty(self.basis.shape[0], dtype=np.float64)
        for i, params in enumerate(self.basis):
            # occup = params[6]
            xyz = params[1:4]
            pos_cart = self.RealspaceMatrix @ xyz
            # iDW, oDW = DWtoDisorder(params[4:6])
            blength[i] = LA.norm(pos_atom_cart - pos_cart)
        return blength

    def disorder_error(self, atomNo):
        dr = np.sqrt(self.basis[atomNo][4] / (8 * np.pi**2))
        dz = np.sqrt(self.basis[atomNo][5] / (8 * np.pi**2))
        errdr = 0.5 * dr * (self.errors[atomNo][4] / self.basis[atomNo][4])
        errdz = 0.5 * dz * (self.errors[atomNo][5] / self.basis[atomNo][5])
        return dr, dz, errdr, errdz

    def __repr__(self):
        repr_super = super().__repr__()
        return f"{repr_super}\n{len(self.names):d} atoms in unit cell"

    def atomToStr(self, no, showErrors=True):
        param = self.basis[no][1:]
        name = self.names[no]
        if (self.errors is not None) and showErrors:
            err = self.errors[no][1:]
            values = [
                f"({param[index]:.5f} +- {err[index]:.5f})"
                for index in range(3)
            ]
            values.extend(
                f"({param[index]:.4f} +- {err[index]:.4f})"
                for index in range(3, 6)
            )
            values.append(f"({param[6]:.0f} +- {err[6]:.0f})")
            return f"{name:<{self._atom_error_column_widths[0]}}" + "".join(
                f"{value:>{width}}"
                for value, width in zip(
                    values,
                    self._atom_error_column_widths[1:],
                )
            )

        formats = (".5f", ".5f", ".5f", ".4f", ".4f", ".4f", ".0f")
        values = [format(value, spec) for value, spec in zip(param, formats)]
        return f"{name:<{self._atom_column_widths[0]}}" + "".join(
            f"{value:>{width}}"
            for value, width in zip(values, self._atom_column_widths[1:])
        )

    def _parameter_header(self, showErrors=True):
        widths = (
            self._atom_error_column_widths
            if self.errors is not None and showErrors
            else self._atom_column_widths
        )
        return "".join(
            f"{name:<{width}}" if index == 0 else f"{name:>{width}}"
            for index, (name, width) in enumerate(
                zip(self._atom_column_names, widths)
            )
        ).rstrip()

    def domainsToStr(self):
        s = ""
        for mat, occu in zip(self.coherentDomainMatrix, self.coherentDomainOccupancy):
            matT = np.vstack((mat[:, :-1], mat[:, -1]))
            s += (
                (f"Coherent {occu:.5f}   ")
                + np.array2string(
                    np.array(matT).flatten(),
                    formatter={"float_kind": lambda x: f"{x:.5f}"},
                    max_line_width=100000,
                )[1:-1]
                + "\n"
            )
        return s

    def latticeRODStr(self):
        a, alpha, _, _ = self.getLatticeParameters()
        return "{:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}".format(
            *a, *np.rad2deg(alpha)
        )  # noqa: E501

    def parameterStr(self, showErrors=True):
        """Return atom and layer parameters as plain text.

        :param bool showErrors:
            Include propagated fit errors when available.
        :returns:
            Serialized unit-cell parameter block.
        :rtype: str
        """
        st = ""
        for i in range(len(self.names)):
            st += str(i).zfill(2) + "  " + self.atomToStr(i, showErrors) + "\n"
        layer_positions = ", ".join(f"{p} = {self.layerpos[p]}" for p in self.layerpos)
        st += f"layerpos: {layer_positions}\n"
        st += f"layer_behaviour: {self.layer_behaviour}\n"
        if self._explicit_layer_cycle is not None:
            st += "layer_cycle: {}\n".format(
                ", ".join(str(layer) for layer in self._explicit_layer_cycle)
            )
        if self.symmetry_metadata is not None:
            from .CTRsymmetry import symmetry_metadata_to_lines

            st += "\n".join(symmetry_metadata_to_lines(self.symmetry_metadata)) + "\n"
        return st

    def parameterStrRod(self):
        st = ""
        for i in range(len(self.names)):
            st += self.atomToStr(i, False) + "\n"
        return st

    def __str__(self):
        st = repr(self) + "\n"
        st += "id  " + self._parameter_header() + "\n"
        return st + self.parameterStr()

    def writeSURfile(self, filename):
        with open(filename, "w") as f:
            st = (
                "return\n"
                + self.domainsToStr()
                + self.latticeRODStr()
                + "\n"
                + self.parameterStrRod()
            )
            f.write(st)

    def toRODStr(self):
        return (
            "return\n"
            + self.domainsToStr()
            + self.latticeRODStr()
            + "\n"
            + self.parameterStrRod()
        )

    def toStr(self, showErrors=True):
        """Serialize the unit cell as plain text.

        :param bool showErrors:
            Include propagated fit errors when available.
        :returns:
            Plain-text unit-cell representation.
        :rtype: str
        """
        return (
            "return\n"
            + self.domainsToStr()
            + self.latticeRODStr()
            + "\n"
            + self._parameter_header(showErrors)
            + "\n"
            + self.parameterStr(showErrors=showErrors)
        )

    def toXYZfile(
        self, xyzfile, ucx=1, ucy=1, ucz=1, translate=np.array([0.0, 0.0, 0.0])
    ):
        xyz_array = self.pos_cart_all(ucx, ucy, ucz, translate)
        with open(xyzfile, "w") as f:
            noatoms = len(xyz_array)
            f.write(f"{noatoms}\n")
            f.write(f"{self.latticeRODStr()}\n")
            np.savetxt(f, xyz_array, fmt=["%s", "%.6f", "%.6f", "%.6f"])

    @staticmethod
    def fromFile(filename):
        f, ext = os.path.splitext(filename)
        if ext == "":
            files = glob.glob(filename + ".*")
            if files:
                extensions = [os.path.splitext(fi)[1] for fi in files]
                for fext in extensions:
                    if fext in [".sur", ".bul"]:
                        ext = os.path.splitext(filename)[1]
                        break
                else:
                    if len(extensions) == 1:
                        ext = extensions[0]
                    else:
                        ext = extensions
            else:
                if os.path.isfile(filename):
                    raise OSError(f"File {filename} has no file extension.")
                else:
                    raise FileNotFoundError(
                        errno.ENOENT, os.strerror(errno.ENOENT), filename
                    )
        elif ext not in [".sur", ".bul"]:
            if not os.path.isfile(filename):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), filename
                )

        if ext == ".sur":
            return UnitCell.fromSURfile(f + ".sur")
        elif ext == ".bul":
            return UnitCell.fromBULfile(f + ".bul")
        else:
            try:
                import ase.io
            except ImportError:
                raise OSError(
                    f"File {filename} has no known file extension. Try to install ASE for more supported file types."  # noqa: E501
                )
            if isinstance(ext, list):
                for e in ext:
                    try:
                        atoms = ase.io.read(f + e)
                        break
                    except Exception:
                        pass
                else:
                    raise OSError(f"Cannot read file {filename}.")
            else:
                atoms = ase.io.read(f + ext)
            if not atoms.cell:
                raise OSError(
                    f"File {filename} contains no valid crystal lattice parameters"
                )
            uc = UnitCell(atoms.cell.lengths(), atoms.cell.angles())
            coord = atoms.get_scaled_positions()
            symb = atoms.get_chemical_symbols()
            for sym, xyz in zip(symb, coord):
                uc.addAtom(sym, xyz, 0.0, 0.0, 1.0)
            return uc

    @classmethod
    def fromSURfile(cls, surfile):
        with open(surfile) as f:
            # next(f) # return line
            line = next_skip_comment(f).split()
            domainoccu = []
            domainmatrix = []
            while line[0] == "Coherent":
                domainoccu.append(float(line[1]))
                mat = np.array(line[2:], dtype=np.float64)
                if mat.size > 10:
                    domainmatrix.append(
                        np.ascontiguousarray(
                            np.vstack(
                                (mat[:-3].reshape((3, 3)).T, mat[-3:].T)
                            ).T.astype(np.float64)
                        )
                    )
                else:
                    domainmatrix.append(
                        np.ascontiguousarray(
                            np.vstack(
                                (mat.reshape((3, 3)).T, np.array([0, 0, 0]))
                            ).T.astype(np.float64)
                        )
                    )
                line = next_skip_comment(f).split()
            if not domainoccu:
                domainoccu.append(1.0)
                domainmatrix.append(
                    np.ascontiguousarray(
                        np.vstack((np.identity(3), np.array([0, 0, 0]))).T.astype(
                            np.float64
                        )
                    )
                )
            latticeparams = np.array(line, dtype=np.float64)
            names = []
            basis = []
            for l in f:  # noqa: E741
                line = l.split()
                if line:
                    names.append(line[0])
                    basis.append(
                        np.concatenate(
                            (
                                [atomic_number(line[0])],
                                np.array(line[1:], dtype=np.float64),
                            )
                        )
                    )
        uc = cls(latticeparams[:3], latticeparams[3:])
        uc.coherentDomainMatrix = domainmatrix
        uc.coherentDomainOccupancy = domainoccu
        if len(basis) >= 2:
            basis = np.vstack(basis)
            if basis.shape[1] == 7:
                basis = np.insert(basis, basis.shape[1], 0, axis=1)  # layer number
            elif basis.shape[1] == 8:
                pass  # all good, layer parameter provided
            else:
                raise ValueError(
                    f"wrong number of atomic parameters. read basis is: {basis}"
                )
            uc.basis = basis
            uc.basis_0 = np.copy(uc.basis)
            uc.names = names
            uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
            uc._test_special_formfactors()
            return uc

        elif len(basis) == 1:
            if basis[0].size == 7:
                basis[0] = np.concatenate((basis[0], [0]))
            elif basis.shape[1] == 8:
                pass  # all good, layer parameter provided
            else:
                raise ValueError(
                    f"wrong number of atomic parameters. read basis is: {basis}"
                )
            uc.basis = np.array([basis])
            uc.dw_increase_constraint = np.ones(basis.shape[0], dtype=np.bool_)
            uc.names = names
            uc.basis_0 = np.copy(uc.basis)
            uc._test_special_formfactors()
            return uc
        else:
            return uc  # no atoms

    @classmethod
    def fromStr(cls, string):
        xprfile = False
        with util.StringIO(string) as f:
            line = next_skip_comment(f).split()
            domainoccu = []
            domainmatrix = []
            while line[0] == "Coherent":
                domainoccu.append(float(line[1]))
                mat = np.array(line[2:], dtype=np.float64)
                if mat.size > 10:
                    domainmatrix.append(
                        np.ascontiguousarray(
                            np.vstack(
                                (mat[:-3].reshape((3, 3)).T, mat[-3:].T)
                            ).T.astype(np.float64)
                        )
                    )
                else:
                    domainmatrix.append(
                        np.ascontiguousarray(
                            np.vstack(
                                (mat.reshape((3, 3)).T, np.array([0, 0, 0]))
                            ).T.astype(np.float64)
                        )
                    )
                line = next_skip_comment(f).split()
            if not domainoccu:
                domainoccu.append(1.0)
                domainmatrix.append(
                    np.ascontiguousarray(
                        np.vstack((np.identity(3), np.array([0, 0, 0]))).T.astype(
                            np.float64
                        )
                    )
                )
            latticeparams = np.array(line, dtype=np.float64)
            names = []
            basis = []
            errors = []
            indices = []
            statistics = dict()
            layerpos = dict()
            layer_behaviour = "ignore"
            layer_cycle = None
            symmetry_lines = []
            reading_symmetry = False
            for l in f:  # noqa: E741
                line = l.rsplit("//")[0]
                stripped = line.strip()
                if reading_symmetry or stripped.startswith(
                    (
                        "spacegroup:",
                        "parent_a:",
                        "parent_alpha:",
                        "surface_transform:",
                        "surface_origin:",
                        "wyckoff_sites:",
                        "wyckoff_atoms:",
                        "wyckoff_couplings:",
                    )
                ):
                    reading_symmetry = True
                    if stripped:
                        symmetry_lines.append(stripped)
                    continue
                if line.startswith("layerpos:"):
                    if "=" in line:
                        try:
                            splitted = [
                                n.split(",")
                                for n in line[len("layerpos:") :].split("=")
                            ]
                            splitted = [
                                item for sublist in splitted for item in sublist
                            ]
                            for i in range(0, len(splitted), 2):
                                layerpos[float(splitted[i].strip())] = float(
                                    splitted[i + 1]
                                )
                        except Exception:
                            print(f"Cannot read layerpos string: {l}")
                    continue
                if line.startswith("layerbehaviour:"):
                    layer_behaviour = line[len("layerbehaviour:") :].strip()
                    continue
                if line.startswith("layer_behaviour:"):
                    layer_behaviour = line[len("layer_behaviour:") :].strip()
                    continue
                if line.startswith("layer_cycle:"):
                    layer_cycle = tuple(
                        float(value) for value in line.split(":", 1)[1].split(",")
                    )
                    continue

                try:
                    if line.strip():
                        if (
                            "Name" in line or "=" in line
                        ):  # parameter or statistics line
                            if "=" in line:
                                try:
                                    splitted = [n.split(",") for n in line.split("=")]
                                    splitted = [
                                        item for sublist in splitted for item in sublist
                                    ]
                                    for i in range(0, len(splitted), 2):
                                        statistics[splitted[i].strip()] = float(
                                            splitted[i + 1]
                                        )
                                except Exception:
                                    print(f"Cannot read statistics string: {l}")
                            continue
                        if "+-" in line:
                            line_sp = line.split()
                            xprfile = True
                            indices.append(int(line_sp[0]))
                            names.append(line_sp[1])
                            params = re.findall(r"\(([^)]+)", line)
                            params_array = np.array(
                                [
                                    np.array(p.split("+-"), dtype=np.float64)
                                    for p in params
                                ]
                            ).T
                            basis.append(
                                np.concatenate(([int(line_sp[0])], params_array[0]))
                            )
                            errors.append(
                                np.concatenate(([int(line_sp[0])], params_array[1]))
                            )
                        else:
                            line_sp = line.split()
                            if line_sp[0].isnumeric():
                                names.append(line_sp[1])
                                basis.append(
                                    np.concatenate(
                                        (
                                            [int(line_sp[0])],
                                            np.array(line_sp[2:], dtype=np.float64),
                                        )
                                    )
                                )
                                indices.append(int(line_sp[0]))
                            else:
                                names.append(line_sp[0])
                                indices.append(np.nan)
                                basis.append(
                                    np.concatenate(
                                        (
                                            [np.nan],
                                            np.array(line_sp[1:], dtype=np.float64),
                                        )
                                    )
                                )

                except Exception as e:
                    raise OSError(
                        f'Cannot read line "{line}" describing atom parameters of UnitCell'  # noqa: E501
                    ) from e
            basis_save = np.vstack(basis).astype(np.float64)
            basis = np.empty_like(basis_save)

            indices = np.array(indices)
            names_save = np.array(names)
            names = np.copy(names_save)

            mask_read = ~np.isnan(indices)  # atoms with explicitly given index
            indices = indices.astype(np.intp)  # convert to int to allow indexing.
            indices_explicit = indices[
                mask_read
            ]  # only the entries, where an index is given

            mask_write = np.zeros_like(mask_read, dtype=np.bool_)
            mask_write[indices_explicit] = (
                True  # mask where to store the atoms (if array is sorted!)
            )
            sortidx = np.argsort(
                indices_explicit
            )  # get order of atoms with explicit indices

            if indices_explicit.size != np.unique(indices_explicit).size:
                raise ValueError(
                    "Atom indices are not unique. This may be caused by assigning two atoms the same index."  # noqa: E501
                )

            basis[mask_write] = basis_save[mask_read][
                sortidx
            ]  # store atoms with explicit index
            names[mask_write] = names_save[mask_read][sortidx]

            basis[~mask_write] = basis_save[
                ~mask_read
            ]  # all other just retain the order of the file
            names[~mask_write] = names_save[~mask_read]

            names = list(names)

            if xprfile:
                errors_save = np.vstack(errors).astype(np.float64)
                errors = np.copy(errors_save)
                errors[mask_write] = errors_save[mask_read][sortidx]
                errors[~mask_write] = errors_save[~mask_read]

        uc = cls(latticeparams[:3], latticeparams[3:])
        uc.coherentDomainMatrix = domainmatrix
        uc.coherentDomainOccupancy = domainoccu
        uc.statistics = statistics
        uc.layerpos = layerpos
        uc.layer_behaviour = layer_behaviour
        uc._explicit_layer_cycle = layer_cycle
        # uc.dw_increase_constraint = np.ones(uc.basis.shape[0],dtype=np.bool_)
        if len(basis) >= 2:
            basis = np.vstack(basis)
            if basis.shape[1] == 7:
                basis = np.insert(basis, basis.shape[1], 0, axis=1)  # layer number
                if xprfile:
                    errors = np.insert(errors, errors.shape[1], np.nan, axis=1)
            elif basis.shape[1] == 8:
                pass  # all good, layer parameter provided
            else:
                raise ValueError(
                    f"wrong number of atomic parameters. read basis is: {basis}"
                )
            if xprfile:
                uc.errors = errors
            uc.names = names
            uc.basis = basis
            uc.basis_0 = np.copy(uc.basis)
            uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
            uc._test_special_formfactors()
            if symmetry_lines:
                from .CTRsymmetry import symmetry_metadata_from_lines

                uc.symmetry_metadata = symmetry_metadata_from_lines(
                    symmetry_lines,
                    uc,
                )
            return uc
        elif len(basis) == 1:
            if (
                basis[0].size == 7
            ):  # basis is already ndarray with ndim=2 for some reason
                basis = np.insert(basis, basis.shape[1], 0, axis=1)
                if xprfile:
                    errors = np.insert(errors, errors.shape[1], np.nan, axis=1)
            elif basis.shape[1] == 8:
                pass  # all good, layer parameter provided
            else:
                raise ValueError(
                    f"wrong number of atomic parameters. read basis is: {basis}"
                )
            uc.basis = basis
            if xprfile:
                uc.errors = errors
            uc.names = names
            uc.basis_0 = np.copy(uc.basis)
            uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
            uc._test_special_formfactors()
            if symmetry_lines:
                from .CTRsymmetry import symmetry_metadata_from_lines

                uc.symmetry_metadata = symmetry_metadata_from_lines(
                    symmetry_lines,
                    uc,
                )
            return uc
        else:
            raise ValueError(
                f"wrong number of atomic parameters. read basis is: {basis}"
            )

    # ugly!!!! please redo this for faster file reads!
    @classmethod
    def fromBULfile(cls, bulfile, DW=0.0):
        with open(bulfile) as f:
            # next_skip_comment(f) # return line
            line = next_skip_comment(f).split()
            domainoccu = []
            domainmatrix = []
            while line[0] == "Coherent":
                domainoccu.append(float(line[1]))
                mat = np.matrix(line[2:], dtype=np.float64)
                if mat.size > 10:
                    domainmatrix.append(
                        np.ascontiguousarray(
                            np.vstack(
                                (mat[:-3].reshape((3, 3)).T, mat[-3:].T)
                            ).T.astype(np.float64)
                        )
                    )
                else:
                    domainmatrix.append(
                        np.ascontiguousarray(
                            np.vstack(
                                (mat.reshape((3, 3)).T, np.array([0, 0, 0]))
                            ).T.astype(np.float64)
                        )
                    )
                line = next_skip_comment(f).split()
            if not domainoccu:
                domainoccu.append(1.0)
                domainmatrix.append(
                    np.ascontiguousarray(
                        np.vstack((np.identity(3), np.array([0, 0, 0]))).T.astype(
                            np.float64
                        )
                    )
                )
            latticeparams = np.array(line, dtype=np.float64)
            # basis = np.loadtxt(f,converters={0 : lambda x :  atomic_number(x.decode("utf-8")) })  # noqa: E501
            # f.seek(basisstartline)
            names = []
            basis = []
            for l in f:  # noqa: E741
                line = l.split()
                if line:
                    names.append(line[0])
                    basis.append(
                        np.concatenate(
                            (
                                [atomic_number(line[0])],
                                np.array(line[1:], dtype=np.float64),
                            )
                        )
                    )
            basis = np.vstack(basis).astype(np.float64)
        uc = cls(latticeparams[:3], latticeparams[3:])
        uc.coherentDomainMatrix = domainmatrix
        uc.coherentDomainOccupancy = domainoccu
        if len(basis) >= 2:
            if basis.shape[1] < 7:
                for i, b in enumerate(basis):
                    uc.addAtom(names[i], b[1:4], DW, DW, 1)
                uc._test_special_formfactors()
                return uc
            elif basis.shape[1] == 7:
                basis = np.insert(basis, basis.shape[1], 0, axis=1)  # layer number
                uc.basis = basis
                uc.names = names
                uc.basis_0 = np.copy(uc.basis)
                uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
                uc._test_special_formfactors()
                return uc
            elif basis.shape[1] == 8:
                uc.basis = basis
                uc.names = names
                uc.basis_0 = np.copy(uc.basis)
                uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
                uc._test_special_formfactors()
                return uc
            else:
                raise ValueError(
                    f"wrong number of atomic parameters. read basis is: {basis}"
                )
        elif len(basis) == 1:
            if basis.shape[1] < 7:
                uc.addAtom(names[0], basis[0][1:4], DW, DW, 1)
                uc._test_special_formfactors()
                return uc
            elif basis.shape[1] == 7:
                basis = np.insert(basis, basis.shape[1], 0, axis=1)  # layer number
                uc.basis = basis
                uc.names = names
                uc.basis_0 = np.copy(uc.basis)
                uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
                uc._test_special_formfactors()
                return uc
            elif basis.shape[1] == 8:
                uc.basis = basis
                uc.names = names
                uc.basis_0 = np.copy(uc.basis)
                uc.dw_increase_constraint = np.ones(uc.basis.shape[0], dtype=np.bool_)
                uc._test_special_formfactors()
                return uc
            else:
                raise ValueError(
                    f"wrong number of atomic parameters. read basis is: {basis}"
                )

        else:
            return uc

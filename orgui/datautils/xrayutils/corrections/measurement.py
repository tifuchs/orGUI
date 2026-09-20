# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
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
r"""Integrated intensities and structure factors on one common scale.

E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532 writes the integrated intensity
of a surface reflection as one prefactor times :math:`|F_{hkl}|^2` times a
product of correction factors. Rocking scans, stationary area-detector
measurements and reflectivity differ only in *which* correction factors
appear. Once every factor is divided out, the three modes land on the same
scale and can be refined as a single data set -- the point of Vlieg's section
5 and of section 3.4 of J. Drnec *et al.*, *J. Appl. Cryst.* **47** (2014)
365.

This module implements that reduction as pure functions, in one place, so
that a rocking scan and a stationary scan of the same rod are guaranteed to
be normalized identically. It is deliberately independent of how the counts
were obtained: :mod:`orgui.app.peak1Dintegr` and :mod:`orgui.app.orGUI`
produce the integrated counts, this module turns them into
:math:`|F_{hkl}|^2`.

The master equation
-------------------

With :math:`I` the *normalized* integrated intensity of
:func:`normalized_intensity` -- counts per second per monitor unit, already
divided by the polarization factor :math:`P` --

.. math::

    I = \Phi_0 \frac{r_e^2 A \lambda^2}{A_u^2}\,
        |F_{hkl}|^2 \, \eta \, C_\mathrm{det}

where

``Phi_0``
    incident flux *density*, photons per second and square meter, so that
    :math:`\Phi_0 A_0` is the total flux on the sample.
``r_e``
    classical electron radius, :data:`CLASSICAL_ELECTRON_RADIUS`, in meter.
``A``
    illuminated *active* surface area, in square meter. This is Vlieg's
    :math:`A = A_0 C_\mathrm{area} C_\mathrm{beam}` (equation 41): the beam
    cross-section corrected for the projection onto the surface and for the
    beam profile over a finite sample. See :ref:`active-area` below.
``lambda``
    wavelength, in Angstrom.
``A_u``
    area of the surface unit cell, in square Angstrom
    (:attr:`orgui.datautils.xrayutils.HKLVlieg.Lattice.uc_area`).
``eta``
    the angle-dependent factor of the measurement mode,
    :func:`angular_factor`.
``C_det``
    in-plane detector acceptance (Vlieg equations 24, 27, 28); ``1`` when the
    region of interest is wide enough to contain the whole in-plane peak
    profile. Not modelled here -- pass a measured or fitted value.

:math:`|F_{hkl}|^2` then comes out in electron units squared, the same scale
a :class:`~orgui.datautils.xrayutils.CTRcalc.SXRDCrystal` calculates in.

The mode-dependent factor
-------------------------

==============================  ============================================
Rocking scan (``th``/omega)     :math:`\eta = L_\varphi\, C_\mathrm{rod}\,
                                \Delta\gamma`
Reflectivity rocking (``mu``)   :math:`\eta = L_r\, C_\mathrm{rod}\,
                                \Delta\gamma`
Stationary area detector        :math:`\eta = L_s`
Specular reflectivity           :math:`\eta = L_s` with
                                :math:`\gamma = \alpha`
==============================  ============================================

with :math:`L_\varphi`, :math:`L_r`, :math:`L_s` and :math:`C_\mathrm{rod}`
the z-axis entries implemented in :mod:`~.geometry`. This table is the single
place the choice is made: :func:`mode_components` names the factors a mode
applies, and both the rocking and the stationary integration ask it rather
than selecting a Lorentz factor of their own.

Two factors of the rocking expression are easy to lose and are the reason a
rocking scan and a stationary scan of the same rod do not otherwise agree:

* the rocking angle must be integrated in **radian**. Integrating in degrees
  scales every rocking structure factor by :math:`180/\pi`.
* :math:`\Delta\gamma`, the **out-of-plane angular acceptance** of the
  region of interest, in radian. A rocking scan intercepts a slice of rod
  whose length is proportional to :math:`\Delta\gamma` (Vlieg equation 20),
  so its integrated intensity is too; a stationary measurement intercepts the
  whole rod cross-section and carries no such factor. Because orGUI sizes
  regions of interest per detector position
  (:func:`orgui.app.ROIutils.calc_corrections`), :math:`\Delta\gamma` is not
  even constant within one scan.

.. _active-area:

The active area
---------------

:math:`A` is whatever surface area actually contributes to the measured
counts. Which of two limits applies is an experimental question, and
:mod:`~.activearea` builds both: :func:`~.activearea.slit_limited_area` when
post-sample slits cut the footprint, :func:`~.activearea.beam_limited_area`
for the open-slit area-detector case orGUI usually runs in.

Either way :math:`A` is the same for a rocking and for a stationary
measurement at the same incidence angle, so it cancels from their ratio and
cannot be the reason the two disagree. It does set the absolute scale.

The explicit total-flux entry points use the algebraically equivalent
:math:`Y=N_\mathrm{net}/(QH)` and
:math:`K=r_e^2\lambda^2/A_u^2`. They do not reinterpret ``flux_density`` or
the legacy active-area APIs.
"""

import numpy as np

from . import geometry

__all__ = [
    "CLASSICAL_ELECTRON_RADIUS",
    "REFLECTIVITY_ROCKING",
    "ROCKING",
    "SPECULAR",
    "STATIONARY",
    "angular_factor",
    "integrated_intensity",
    "mode_components",
    "normalized_intensity",
    "photon_yield_from_structure_factor",
    "reflectivity_from_structure_factor",
    "scale_factor",
    "structure_factor_from_reflectivity",
    "structure_factor_squared",
    "structure_factor_squared_from_photon_yield",
    "total_flux_prefactor",
]

#: Classical electron radius in meter (CODATA 2018).
CLASSICAL_ELECTRON_RADIUS = 2.8179403262e-15

#: Rocking scan about the sample rotation; see :func:`angular_factor`.
ROCKING = geometry.ROCKING

#: Rocking scan about the incidence angle, how a reflectivity curve is
#: measured; see :func:`angular_factor`.
REFLECTIVITY_ROCKING = geometry.REFLECTIVITY_ROCKING

#: Stationary measurement on an area detector; see :func:`angular_factor`.
STATIONARY = geometry.STATIONARY

#: Specular reflectivity: the stationary factor at ``gamma = alpha``.
SPECULAR = "specular"

#: Modes whose integrated intensity is proportional to the out-of-plane
#: acceptance of the region of interest, because they intercept a slice of
#: rod rather than its whole cross-section.
_ROCKING_MODES = (ROCKING, REFLECTIVITY_ROCKING)

_ANGSTROM = 1e-10


def _reflectivity_prefactor(wavelength, unitcell_area):
    r""":math:`r_e^2\lambda^2/A_u^2`, dimensionless.

    The same combination :func:`scale_factor` builds, without the flux and
    the area, which is what Vlieg's reflectivity expression (equation 63)
    needs. Kept separate so the reflectivity functions do not depend on
    :func:`scale_factor`'s unit defaults being ``1``.

    :param wavelength: X-ray wavelength, in Angstrom.
    :param unitcell_area: Surface unit-cell area, in square Angstrom.
    :rtype: numpy.ndarray
    :raises ValueError: If either argument is not positive.
    """
    wavelength = np.asarray(wavelength, dtype=np.float64)
    unitcell_area = np.asarray(unitcell_area, dtype=np.float64)
    if np.any(wavelength <= 0):
        raise ValueError("wavelength must be positive, in Angstrom")
    if np.any(unitcell_area <= 0):
        raise ValueError("unitcell_area must be positive, in square Angstrom")
    lam = wavelength * _ANGSTROM
    a_u = unitcell_area * _ANGSTROM**2
    return CLASSICAL_ELECTRON_RADIUS**2 * lam**2 / a_u**2


def normalized_intensity(counts, exposure_time=1.0, monitor=1.0, angle_unit=None):
    r"""Counts brought to per second, per monitor unit, per radian.

    Both Vlieg's rocking-scan expression (equation 42) and his stationary one
    (equation 54) contain the incident flux and the counting time, so neither
    an integrated rocking curve nor a stationary frame means anything until
    it is divided by the counting time and by whatever monitor the flux is
    tracked with. Doing it in one place is what lets the two modes share a
    scale.

    :param counts: Background-subtracted integrated counts. For a stationary
        measurement, the sum over the region of interest of one frame. For a
        rocking scan, the integral of the per-frame counts over the rocking
        angle, in counts times ``angle_unit``.
    :param exposure_time: Counting time of a single frame, in seconds. A
        rocking scan uses the *per-frame* time, not the sum over the scan:
        the rocking angle is the integration variable, the time is not.
    :param monitor: Monitor counter value the flux is normalized against.
        ``1.0`` leaves the intensity in counts per second.
    :param str angle_unit: ``None`` for a stationary sum, ``'rad'`` or
        ``'deg'`` for a rocking-scan integral. ``'deg'`` converts to radian,
        which is what :func:`angular_factor` and the published equations
        assume.
    :returns: The normalized intensity, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If ``angle_unit`` is not ``None``, ``'rad'`` or
        ``'deg'``, or if a divisor is not finite and non-zero.
    """
    counts = np.asarray(counts, dtype=np.float64)
    exposure_time = np.asarray(exposure_time, dtype=np.float64)
    monitor = np.asarray(monitor, dtype=np.float64)
    if np.any(exposure_time <= 0) or not np.all(np.isfinite(exposure_time)):
        raise ValueError("exposure_time must be finite and positive")
    if np.any(monitor == 0) or not np.all(np.isfinite(monitor)):
        raise ValueError("monitor must be finite and non-zero")

    if angle_unit is None:
        scale = 1.0
    elif angle_unit == "rad":
        scale = 1.0
    elif angle_unit == "deg":
        scale = np.deg2rad(1.0)
    else:
        raise ValueError(
            f"unknown angle_unit {angle_unit!r}; expected None, 'rad' or 'deg'"
        )
    return counts * scale / (exposure_time * monitor)


def _resolve_specular(mode, alpha, gamma):
    """Rewrite :data:`SPECULAR` as the stationary case at ``gamma = alpha``.

    :returns: ``(mode, gamma)`` with the specular mode replaced.
    :rtype: tuple
    :raises ValueError: If ``alpha`` is missing, or ``gamma`` contradicts the
        specular condition.
    """
    if mode != SPECULAR:
        return mode, gamma
    if alpha is None:
        raise ValueError("the specular factor needs alpha")
    if gamma is not None and not np.allclose(gamma, alpha):
        raise ValueError(
            "the specular factor is the stationary one at gamma = alpha; "
            "use STATIONARY for a non-specular exit angle"
        )
    return STATIONARY, alpha


def mode_components(mode, alpha=None, delta=None, gamma=None):
    r"""Named angular factors a measurement mode applies.

    The single place the choice of correction factors is tied to the kind of
    scan. Both the rocking and the stationary integration ask this rather
    than selecting a Lorentz factor themselves, and they store the returned
    names beside ``F2_hkl``, so the saved data records which mode produced
    it.

    Returned keys:

    ``C_Lorentz``
        :math:`L_\varphi = 1/(\sin\delta\cos\alpha\cos\gamma)` for
        :data:`ROCKING`, :math:`L_r = 1/\sin 2\alpha` for
        :data:`REFLECTIVITY_ROCKING`, :math:`L_s = 1/\sin\gamma` for
        :data:`STATIONARY` and :data:`SPECULAR`.
    ``C_rod``
        :math:`\cos\gamma`, the rod interception, for the two rocking modes
        only. A stationary measurement integrates across the whole rod and
        has none (Vlieg 1997, before equation 54).

    :param str mode: :data:`ROCKING`, :data:`REFLECTIVITY_ROCKING`,
        :data:`STATIONARY` or :data:`SPECULAR`.
    :param alpha: Incidence angle, in radian.
    :param delta: In-plane detector angle, in radian.
    :param gamma: Out-of-plane detector angle -- the exit angle in the z-axis
        geometry -- in radian.
    :returns: Mapping of factor name to value, broadcast over the inputs.
    :rtype: dict
    :raises ValueError: If ``mode`` is unknown or a required angle is
        missing.
    """
    mode, gamma = _resolve_specular(mode, alpha, gamma)
    if mode == ROCKING:
        return {
            "C_Lorentz": geometry.lorentz_factor(
                geometry.ROCKING, alpha=alpha, delta=delta, gamma=gamma
            ),
            "C_rod": geometry.rod_interception(gamma),
        }
    if mode == REFLECTIVITY_ROCKING:
        return {
            "C_Lorentz": geometry.lorentz_factor(
                geometry.REFLECTIVITY_ROCKING, alpha=alpha
            ),
            "C_rod": geometry.rod_interception(gamma),
        }
    if mode == STATIONARY:
        return {
            "C_Lorentz": geometry.lorentz_factor(geometry.STATIONARY, gamma=gamma)
        }
    raise ValueError(
        f"unknown mode {mode!r}; expected one of {ROCKING!r}, "
        f"{REFLECTIVITY_ROCKING!r}, {STATIONARY!r} or {SPECULAR!r}"
    )


def angular_factor(
    mode, alpha=None, delta=None, gamma=None, detector_acceptance=None
):
    r"""Angle-dependent factor :math:`\eta` of the integrated intensity.

    The product of the :func:`mode_components` of the mode and, for the two
    rocking modes, the out-of-plane acceptance: :data:`ROCKING` gives
    :math:`L_\varphi\,C_\mathrm{rod}\,\Delta\gamma` (Vlieg equations 16, 20
    and 23, Drnec equation 3), :data:`STATIONARY` gives
    :math:`L_s = 1/\sin\gamma` (Vlieg equation 53, Drnec equation 5), and
    :data:`SPECULAR` is the stationary factor at :math:`\gamma = \alpha`,
    the case Vlieg treats in his section 3.2.

    :param str mode: :data:`ROCKING`, :data:`REFLECTIVITY_ROCKING`,
        :data:`STATIONARY` or :data:`SPECULAR`.
    :param alpha: Incidence angle, in radian. Required for the rocking modes
        and for :data:`SPECULAR`.
    :param delta: In-plane detector angle, in radian. Required for
        :data:`ROCKING`.
    :param gamma: Out-of-plane detector angle, in radian. Required for
        :data:`ROCKING`, :data:`REFLECTIVITY_ROCKING` and
        :data:`STATIONARY`.
    :param detector_acceptance: Out-of-plane angular acceptance
        :math:`\Delta\gamma` of the region of interest, in **radian**.
        Required for the rocking modes, rejected otherwise: a stationary
        measurement intercepts the whole rod and does not depend on it.
    :returns: :math:`\eta`, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If ``mode`` is unknown, a required angle is missing,
        or ``detector_acceptance`` is given for a mode that has none.
    """
    resolved, _ = _resolve_specular(mode, alpha, gamma)
    if resolved in _ROCKING_MODES:
        if detector_acceptance is None:
            raise ValueError(
                "a rocking scan needs the out-of-plane acceptance of its "
                "region of interest, in radian; its integrated intensity is "
                "proportional to it (Vlieg 1997, equations 20 and 42)"
            )
        acceptance = np.asarray(detector_acceptance, dtype=np.float64)
        if np.any(acceptance <= 0) or not np.all(np.isfinite(acceptance)):
            raise ValueError("detector_acceptance must be finite and positive")
    elif detector_acceptance is not None:
        raise ValueError(
            "a stationary measurement intercepts the whole rod cross "
            "section and has no detector-acceptance factor "
            "(Vlieg 1997, equation 54)"
        )
    else:
        acceptance = 1.0

    factor = acceptance
    for value in mode_components(
        mode, alpha=alpha, delta=delta, gamma=gamma
    ).values():
        factor = factor * value
    return factor


def scale_factor(wavelength, unitcell_area, active_area=1.0, flux_density=1.0):
    r"""Mode-independent prefactor
    :math:`\Phi_0 r_e^2 A \lambda^2 / A_u^2`, in 1/s.

    Leaving ``flux_density`` and ``active_area`` at ``1`` gives the relative
    scale factor, which is all that is needed to put different scan modes of
    one experiment on a *common* scale. Supplying the measured flux density
    and the illuminated area additionally puts them on the *absolute* one, so
    that :func:`structure_factor_squared` returns electron units.

    :param wavelength: X-ray wavelength, in Angstrom.
    :param unitcell_area: Area :math:`A_u` of the surface unit cell, in
        square Angstrom.
    :param active_area: Illuminated active surface area :math:`A`, in square
        meter; see :mod:`~.activearea`.
    :param flux_density: Incident flux density :math:`\Phi_0`, in photons per
        second and square meter, so that ``flux_density * A_0`` is the flux
        on the sample.
    :returns: The prefactor, in 1/s, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If the wavelength or the unit-cell area is not
        positive.
    """
    wavelength = np.asarray(wavelength, dtype=np.float64)
    unitcell_area = np.asarray(unitcell_area, dtype=np.float64)
    if np.any(wavelength <= 0):
        raise ValueError("wavelength must be positive, in Angstrom")
    if np.any(unitcell_area <= 0):
        raise ValueError("unitcell_area must be positive, in square Angstrom")
    lam = wavelength * _ANGSTROM
    a_u = unitcell_area * _ANGSTROM**2
    return (
        np.asarray(flux_density, dtype=np.float64)
        * np.asarray(active_area, dtype=np.float64)
        * CLASSICAL_ELECTRON_RADIUS**2
        * lam**2
        / a_u**2
    )


def total_flux_prefactor(wavelength, unitcell_area):
    r"""Total-flux CTR prefactor :math:`K=r_e^2\lambda^2/A_u^2`.

    This is the prefactor paired with photon-normalized yield
    :math:`Y=N_\mathrm{net}/(QH)`. It contains neither a flux density nor an
    active area: the total incident frame fluence :math:`Q` and the
    dimensionless illumination divisor :math:`H` replace that legacy product.

    :param wavelength: X-ray wavelength, in Angstrom.
    :param unitcell_area: Surface unit-cell area, in square Angstrom.
    :returns: Dimensionless prefactor, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If either physical input is not finite and positive.
    """
    wavelength = np.asarray(wavelength, dtype=np.float64)
    unitcell_area = np.asarray(unitcell_area, dtype=np.float64)
    if not np.all(np.isfinite(wavelength)):
        raise ValueError("wavelength must be finite and positive, in Angstrom")
    if not np.all(np.isfinite(unitcell_area)):
        raise ValueError(
            "unitcell_area must be finite and positive, in square Angstrom"
        )
    return _reflectivity_prefactor(wavelength, unitcell_area)


def _total_flux_efficiency(detector_efficiency):
    efficiency = np.asarray(detector_efficiency, dtype=np.float64)
    if np.any(efficiency <= 0) or not np.all(np.isfinite(efficiency)):
        raise ValueError("detector_efficiency must be finite and positive")
    return efficiency


def structure_factor_squared_from_photon_yield(
    photon_yield,
    mode,
    alpha=None,
    delta=None,
    gamma=None,
    detector_acceptance=None,
    wavelength=None,
    unitcell_area=None,
    detector_efficiency=1.0,
):
    r"""Reduce :math:`Y=N_\mathrm{net}/(QH)` to :math:`|F_{hkl}|^2`.

    For a rocking scan ``photon_yield`` is the angular integral of the
    framewise yield in radians. For a stationary frame it is the scalar yield
    itself. The mode-dependent factor remains :func:`angular_factor`; only the
    illumination/fluence convention differs from
    :func:`structure_factor_squared`.

    :param photon_yield: Counts divided framewise by calibrated incident
        photons ``Q`` and dimensionless illumination ``H``.
    :param str mode: Measurement mode accepted by :func:`angular_factor`.
    :param alpha: Incidence angle, radian.
    :param delta: In-plane detector angle, radian.
    :param gamma: Out-of-plane detector angle, radian.
    :param detector_acceptance: Rocking out-of-plane acceptance, radian.
    :param wavelength: X-ray wavelength, Angstrom.
    :param unitcell_area: Surface unit-cell area, square Angstrom.
    :param detector_efficiency: Further multiplicative detector response.
    :returns: Structure factor squared, in electron units squared for a
        calibrated ``Q`` and stated detector response.
    :rtype: numpy.ndarray
    :raises ValueError: If scale inputs or detector efficiency are invalid.
    """
    if wavelength is None or unitcell_area is None:
        raise ValueError(
            "wavelength (Angstrom) and unitcell_area (square Angstrom) set the "
            "scale and must both be given"
        )
    eta = angular_factor(
        mode,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=detector_acceptance,
    )
    return np.asarray(photon_yield, dtype=np.float64) / (
        total_flux_prefactor(wavelength, unitcell_area)
        * eta
        * _total_flux_efficiency(detector_efficiency)
    )


def photon_yield_from_structure_factor(
    f2,
    mode,
    alpha=None,
    delta=None,
    gamma=None,
    detector_acceptance=None,
    wavelength=None,
    unitcell_area=None,
    detector_efficiency=1.0,
):
    r"""Forward total-flux model, inverse of the yield reduction.

    Arguments and units match
    :func:`structure_factor_squared_from_photon_yield`; ``f2`` is in electron
    units squared and the result is the expected counts per ``Q H``.

    :returns: Photon-normalized detector yield.
    :rtype: numpy.ndarray
    """
    if wavelength is None or unitcell_area is None:
        raise ValueError(
            "wavelength (Angstrom) and unitcell_area (square Angstrom) set the "
            "scale and must both be given"
        )
    eta = angular_factor(
        mode,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=detector_acceptance,
    )
    return (
        np.asarray(f2, dtype=np.float64)
        * total_flux_prefactor(wavelength, unitcell_area)
        * eta
        * _total_flux_efficiency(detector_efficiency)
    )


def structure_factor_squared(
    intensity,
    mode,
    alpha=None,
    delta=None,
    gamma=None,
    detector_acceptance=None,
    wavelength=None,
    unitcell_area=None,
    active_area=1.0,
    flux_density=1.0,
    detector_efficiency=1.0,
):
    r"""Turn a normalized integrated intensity into :math:`|F_{hkl}|^2`.

    Inverts the master equation of this module. ``intensity`` must already be
    the output of :func:`normalized_intensity` and must already be divided by
    the polarization factor :math:`P`; everything else is divided out here.

    :param intensity: Normalized integrated intensity, from
        :func:`normalized_intensity`.
    :param str mode: :data:`ROCKING`, :data:`STATIONARY` or :data:`SPECULAR`.
    :param alpha: Incidence angle, in radian.
    :param delta: In-plane detector angle, in radian.
    :param gamma: Out-of-plane detector angle, in radian.
    :param detector_acceptance: Out-of-plane acceptance :math:`\Delta\gamma`
        of the region of interest, in radian; rocking scans only.
    :param wavelength: X-ray wavelength, in Angstrom.
    :param unitcell_area: Surface unit-cell area, in square Angstrom.
    :param active_area: Illuminated active area, in square meter.
    :param flux_density: Incident flux density, in photons per second and
        square meter.
    :param detector_efficiency: In-plane acceptance :math:`C_\mathrm{det}`,
        or any further multiplicative correction of the *intensity*.
    :returns: :math:`|F_{hkl}|^2`, in electron units squared when the flux
        density and the active area are the measured ones, and on a common
        but arbitrary scale otherwise.
    :rtype: numpy.ndarray
    :raises ValueError: If the mode or the angles are inconsistent, or if the
        wavelength or unit-cell area is missing.
    """
    if wavelength is None or unitcell_area is None:
        raise ValueError(
            "wavelength (Angstrom) and unitcell_area (square Angstrom) set the "
            "scale and must both be given"
        )
    eta = angular_factor(
        mode,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=detector_acceptance,
    )
    scale = scale_factor(wavelength, unitcell_area, active_area, flux_density)
    return np.asarray(intensity, dtype=np.float64) / (
        scale * eta * np.asarray(detector_efficiency, dtype=np.float64)
    )


def integrated_intensity(
    f2,
    mode,
    alpha=None,
    delta=None,
    gamma=None,
    detector_acceptance=None,
    wavelength=None,
    unitcell_area=None,
    active_area=1.0,
    flux_density=1.0,
    detector_efficiency=1.0,
):
    r"""Forward model: the normalized intensity a given :math:`|F_{hkl}|^2`
    produces.

    The exact inverse of :func:`structure_factor_squared`, kept as a public
    function because it is how a simulated measurement -- and therefore a
    regression test of the whole reduction -- is built.

    Arguments and units are those of :func:`structure_factor_squared`, with
    ``f2`` in electron units squared.

    :returns: The normalized integrated intensity, in the units
        :func:`normalized_intensity` produces.
    :rtype: numpy.ndarray
    """
    eta = angular_factor(
        mode,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=detector_acceptance,
    )
    scale = scale_factor(wavelength, unitcell_area, active_area, flux_density)
    return (
        np.asarray(f2, dtype=np.float64)
        * scale
        * eta
        * np.asarray(detector_efficiency, dtype=np.float64)
    )


def reflectivity_from_structure_factor(
    f2, wavelength, unitcell_area, alpha, beta_out=None, polarization=1.0
):
    r"""Absolute reflectivity of a rod, from its structure factor.

    Vlieg equation 63,

    .. math::

        R = \frac{r_e^2 \lambda^2 P_r}{A_u^2 \sin\alpha \sin\beta_\mathrm{out}}
            |F_{hkl}|^2 ,

    the fraction of the incident flux scattered into the rod. For the
    specular rod :math:`\beta_\mathrm{out} = \alpha` and this is the familiar
    :math:`1/\sin^2\alpha`: one power from the beam footprint
    :math:`A_r = A_0/\sin\alpha`, one from the stationary Lorentz factor
    :math:`1/\sin\beta_\mathrm{out}`.

    Nothing beyond the structure factor is needed, so an absolutely scaled
    :math:`|F_{hkl}|^2` from any scan mode yields an absolute reflectivity
    without a separate measurement. This is the kinematic result: it does not
    hold near a bulk Bragg peak, and below the critical angle refraction and
    multiple scattering are not described by it.

    Off the specular rod, :math:`R` is the fraction of the incident flux
    scattered into that rod; it is a well-defined number for a truncation rod
    integrated across its cross-section, but not for diffuse scattering,
    where only a differential cross-section is meaningful.

    :param f2: :math:`|F_{hkl}|^2` in electron units squared, on an absolute
        scale.
    :param wavelength: X-ray wavelength, in Angstrom.
    :param unitcell_area: Surface unit-cell area, in square Angstrom.
    :param alpha: Incidence angle, in radian.
    :param beta_out: Exit angle, in radian; ``None`` means the specular
        condition ``beta_out = alpha``.
    :param polarization: Polarization factor :math:`P_r`; for the specular
        condition Vlieg equation 59 gives
        ``p_h * cos(2 alpha)**2 + (1 - p_h)``.
    :returns: The reflectivity, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If the wavelength or the unit-cell area is not
        positive.
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    beta_out = alpha if beta_out is None else np.asarray(beta_out, dtype=np.float64)
    return (
        np.asarray(f2, dtype=np.float64)
        * _reflectivity_prefactor(wavelength, unitcell_area)
        * np.asarray(polarization, dtype=np.float64)
        / (np.sin(alpha) * np.sin(beta_out))
    )


def structure_factor_from_reflectivity(
    reflectivity, wavelength, unitcell_area, alpha, beta_out=None, polarization=1.0
):
    r"""Inverse of :func:`reflectivity_from_structure_factor`.

    Puts a measured reflectivity curve on the same :math:`|F_{hkl}|^2` scale
    as the crystal truncation rods, so that a reflectivity and a set of rods
    can be refined together.

    Arguments and units are those of
    :func:`reflectivity_from_structure_factor`.

    :returns: :math:`|F_{hkl}|^2` in electron units squared.
    :rtype: numpy.ndarray
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    beta_out = alpha if beta_out is None else np.asarray(beta_out, dtype=np.float64)
    return (
        np.asarray(reflectivity, dtype=np.float64)
        * np.sin(alpha)
        * np.sin(beta_out)
        / (
            _reflectivity_prefactor(wavelength, unitcell_area)
            * np.asarray(polarization, dtype=np.float64)
        )
    )

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
r"""Adapter between the integration workflow and the correction physics.

Every correction factor itself lives in
:mod:`orgui.datautils.xrayutils.corrections`, which takes numbers and returns
numbers and knows nothing about scans, configurations or widgets. What is
left here is the translation in the other direction: turning *which
corrections were switched on* and *what a loaded scan calls its counters*
into arguments for those functions, and assembling the result into the bundle
that is stored beside the integrated intensity.

The normalization and numerical active-area factor are intensity divisors:

.. math::

    I_\mathrm{corr} = \frac{I}{C_\mathrm{norm}\,C_\mathrm{area}}
    \qquad
    F^2_{hkl} = \frac{I_\mathrm{corr}}{L_\mathrm{stationary}}

``C_flux_on_sample`` is retained as a diagnostic component of
``C_illum_area``. It must not be divided out separately: the numerical active
area already contains the same beam/sample overlap integral. Which angular
factors a stationary measurement applies -- and that it has no
rod-interception factor, unlike a rocking scan -- follows the z-axis table of
:mod:`~orgui.datautils.xrayutils.corrections.geometry`.

The exposure and monitor normalization mirrors the reciprocal-space
reconstruction (:mod:`orgui.reconstruction_job`), so a stationary integration
and a reconstruction of the same scan are normalized identically.

This module holds no Qt state and reads only public scan attributes, so it is
safe in CLI and batch use.
"""

import numpy as np

from ..datautils.xrayutils.corrections import geometry
from ..datautils.xrayutils.corrections.normalization import (
    normalization_divisor as _divisor_from_counters,
)
from ..datautils.xrayutils.corrections.roi import (  # noqa: F401
    CorrectionFactors,
    roi_mean_correction,
)

__all__ = [
    "CorrectionFactors",
    "apply_stationary_corrections",
    "monitor_counter_candidates",
    "normalization_divisor",
    "roi_mean_correction",
    "stationary_correction_factors",
    "structure_factor",
]


def monitor_counter_candidates(scan):
    """Counter names of ``scan`` usable as a monitor normalization.

    A counter qualifies when it is a one-dimensional array with one value per
    image, which is what a divisive normalization needs.

    :param scan: A loaded scan object.
    :returns: Sorted counter names.
    :rtype: list[str]
    """
    names = []
    candidates = set(getattr(scan, "auxillary_counters", ()) or ())
    try:
        length = len(scan)
    except TypeError:
        length = None
    for name in candidates:
        value = getattr(scan, name, None)
        if value is None:
            continue
        array = np.atleast_1d(np.asarray(value))
        if array.ndim != 1:
            continue
        if length is not None and array.size not in (1, length):
            continue
        names.append(name)
    return sorted(names)


def normalization_divisor(scan, normalize_exposure, monitor_corrections, size):
    """Exposure and monitor divisor for every image of a scan.

    Pulls the counters off the scan object and hands them to
    :func:`~orgui.datautils.xrayutils.corrections.normalization.normalization_divisor`,
    which owns the arithmetic and the validation. Mirrors the reciprocal-space
    reconstruction, which divides each frame by its exposure time and by every
    configured monitor counter.

    :param scan: Scan object providing ``exposure_time`` and the monitor
        counters by attribute.
    :param bool normalize_exposure: Divide by ``scan.exposure_time`` when the
        scan provides it.
    :param monitor_corrections: Iterable of monitor counter names.
    :param int size: Number of images, used to broadcast scalar counters.
    :returns: ``(divisor, applied)`` -- an array of shape ``(size,)`` and the
        names of the normalizations that contributed.
    :rtype: tuple
    :raises ValueError: If a counter is missing, or has a non-positive or
        non-finite value that would make the normalization undefined.
    """
    exposure = None
    if normalize_exposure:
        # A scan backend that reports no exposure time is not an error; the
        # reconstruction records the normalization as unavailable rather than
        # failing the job, and this follows it.
        exposure = getattr(scan, "exposure_time", None)

    monitors = {}
    for name in monitor_corrections:
        if not hasattr(scan, name):
            raise ValueError(f"Active scan has no monitor counter named {name!r}")
        monitors[name] = getattr(scan, name)

    return _divisor_from_counters(
        size, exposure_time=exposure, monitors=monitors
    )


def stationary_correction_factors(
    alpha,
    delta,
    gamma,
    use_lorentz=False,
    use_footprint=False,
    beam_profile=None,
    sample_size=None,
    normalization=None,
):
    r"""Correction divisors for one stationary-scan trajectory.

    :param alpha: Incidence angle per image, in radian. This is the ``mu``
        circle, and the angle the footprint corrections depend on.
    :param delta: In-plane detector angle per image, in radian.
    :param gamma: Out-of-plane detector angle per image, in radian.
    :param bool use_lorentz: Add ``C_Lorentz`` -- the *stationary-mode*
        factor :math:`1/\sin\gamma`, not the rocking-scan one.
        Stationary integration has no rod-interception factor.
    :param bool use_footprint: Add the ``C_illum_area`` divisor and its
        diagnostic numerator ``C_flux_on_sample``.
    :param beam_profile: A
        :class:`~orgui.datautils.xrayutils.corrections.beamprofile.BeamProfile`,
        required when ``use_footprint`` is set.
    :param float sample_size: Sample size along the beam in meters, required
        when ``use_footprint`` is set.
    :param normalization: Optional per-image exposure and monitor divisor
        from :func:`normalization_divisor`, stored as ``C_norm``.
    :returns: The factors, each broadcast to the shape of ``alpha``.
    :rtype: CorrectionFactors
    :raises ValueError: If the footprint correction is requested without a
        beam profile or a sample size.
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    factors = {}
    applied = []

    if normalization is not None:
        factors["C_norm"] = np.broadcast_to(
            np.asarray(normalization, dtype=np.float64), alpha.shape
        ).copy()
        applied.append("normalization")

    if use_footprint:
        if beam_profile is None:
            raise ValueError(
                "the footprint correction needs a beam profile; configure one "
                "in the integration corrections dialog"
            )
        if sample_size is None or not sample_size > 0:
            raise ValueError(
                f"the footprint correction needs a positive sample size, got "
                f"{sample_size!r}"
            )
        flux, area = beam_profile.corrections(alpha, sample_size)
        factors["C_flux_on_sample"] = np.broadcast_to(flux, alpha.shape).copy()
        factors["C_illum_area"] = np.broadcast_to(area, alpha.shape).copy()
        applied.append("footprint")

    if use_lorentz:
        factors["C_Lorentz"] = np.broadcast_to(
            geometry.lorentz_factor(geometry.STATIONARY, gamma=gamma), alpha.shape
        ).copy()
        applied.append("lorentz")

    return CorrectionFactors(factors, applied)


def apply_stationary_corrections(intensity, errors, factors):
    """Divide an intensity and its errors by the intensity-level factors.

    ``C_flux_on_sample`` is not an additional divisor. The numerical active
    area ``C_illum_area`` already contains that overlap integral. The Lorentz
    factor is applied separately by :func:`structure_factor`.

    :param intensity: Integrated intensity per image.
    :param errors: 1-sigma errors of ``intensity``.
    :param CorrectionFactors factors: Divisors to apply.
    :returns: ``(intensity, errors)`` corrected.
    :rtype: tuple of numpy.ndarray
    """
    divisor = factors.divisor("C_norm", "C_illum_area")
    return np.asarray(intensity) / divisor, np.asarray(errors) / divisor


def structure_factor(intensity, errors, factors):
    r"""Form :math:`F^2_{hkl}` from an already corrected intensity.

    :math:`F^2 = I_\mathrm{corr} / L_\mathrm{stationary}`. Unlike a rocking
    scan, stationary area-detector integration has no rod-interception
    factor.

    :param intensity: Corrected intensity per image.
    :param errors: 1-sigma errors of ``intensity``.
    :param CorrectionFactors factors: Must contain ``C_Lorentz``.
    :returns: ``(F2_hkl, F2_hkl_errors)``.
    :rtype: tuple of numpy.ndarray
    :raises KeyError: If the Lorentz factors are absent.
    """
    divisor = factors["C_Lorentz"]
    return np.asarray(intensity) / divisor, np.asarray(errors) / divisor

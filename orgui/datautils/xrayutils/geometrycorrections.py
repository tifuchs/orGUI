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
r"""Geometrical correction factors of the z-axis diffractometer.

orGUI applies the correction factors tabulated for the **z-axis geometry** in
Appendix A of the ANA/ROD manual (E. Vlieg, *ANA -- program for the analysis
of surface X-ray diffraction data*), which reproduces
E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532. The z-axis, 5-circle and
6-circle geometries are special cases of one another; only the z-axis column
is implemented here, and it is the one orGUI's angle convention matches.

Angles follow :mod:`orgui.datautils.xrayutils.HKLVlieg`:

``alpha``
    incidence angle of the beam on the surface, in radian. This is the ``mu``
    circle of the diffractometer, and the angle the footprint corrections of
    :mod:`orgui.datautils.xrayutils.beamprofile` depend on.
``delta``
    in-plane detector angle, in radian.
``gamma``
    out-of-plane detector angle, in radian. In the z-axis geometry this is
    the exit angle of the diffracted beam from the surface.

The table entries used are, verbatim from the manual:

=========================================  ================================
Lorentz factor rocking scan                :math:`1/(\sin\delta\,
                                           \cos\alpha\,\cos\gamma)`
Lorentz factor stationary mode             :math:`1/\sin\gamma`
Lorentz factor reflectivity rocking scan   :math:`1/\sin 2\alpha`
Rod interception                           :math:`\cos\gamma`
Area correction (ignoring footprint and    :math:`1/\sin\delta`
sample size)
Beam profile and finite sample size        calculated numerically, see
                                           :mod:`~.beamprofile`
=========================================  ================================

Structure factors follow orGUI's established convention, in which the
tabulated Lorentz factor and the rod interception are *divided out* of the
integrated intensity:

.. math:: F^2_{hkl} = \frac{I}{L \cdot C_\mathrm{rod}}

The three Lorentz factors are alternatives, selected by how the intensity was
measured, not factors to be combined:

* a rocking scan about the sample rotation (``th``/``omega``) uses the
  rocking-scan factor,
* a rocking scan about the incidence angle (``mu``), which is how a
  reflectivity curve is measured, uses the reflectivity factor,
* a scan with the sample stationary, integrated across the rod on an area
  detector, uses the stationary factor.

The area correction :math:`1/\sin\delta` is listed for completeness and is
**not** applied by orGUI; the numerically evaluated beam-profile and finite
sample-size corrections of :mod:`~.beamprofile` are used instead.
"""

import numpy as np

__all__ = [
    "AREA",
    "REFLECTIVITY_ROCKING",
    "ROCKING",
    "STATIONARY",
    "area_correction",
    "lorentz_factor",
    "lorentz_reflectivity_rocking_scan",
    "lorentz_rocking_scan",
    "lorentz_stationary",
    "rod_interception",
]

#: Lorentz factor of a rocking scan about the sample rotation.
ROCKING = "rocking"

#: Lorentz factor of a rocking scan about the incidence angle.
REFLECTIVITY_ROCKING = "reflectivity_rocking"

#: Lorentz factor of a scan measured with the sample stationary.
STATIONARY = "stationary"

#: Identifier of the area correction, for provenance strings.
AREA = "area"


def lorentz_rocking_scan(delta, alpha, gamma):
    r"""Lorentz factor of a rocking scan about the sample rotation.

    :math:`L = 1/(\sin\delta\,\cos\alpha\,\cos\gamma)`, the z-axis
    "Lorentz factor rocking scan" entry. The absolute value is taken so that
    integrated intensities stay positive whichever way the scan runs.

    :param delta: In-plane detector angle, in radian.
    :param alpha: Incidence angle, in radian.
    :param gamma: Out-of-plane detector angle, in radian.
    :returns: The Lorentz factor, broadcast over the inputs.
    :rtype: numpy.ndarray
    """
    return np.abs(1.0 / (np.sin(delta) * np.cos(alpha) * np.cos(gamma)))


def lorentz_reflectivity_rocking_scan(alpha):
    r"""Lorentz factor of a reflectivity scan rocked in the incidence angle.

    :math:`L = 1/\sin 2\alpha`, the z-axis "Lorentz factor reflectivity
    rocking scan" entry.

    :param alpha: Incidence angle, in radian.
    :returns: The Lorentz factor, broadcast over the input.
    :rtype: numpy.ndarray
    """
    return np.abs(1.0 / np.sin(2.0 * np.asarray(alpha, dtype=float)))


def lorentz_stationary(gamma):
    r"""Lorentz factor of a measurement with the sample stationary.

    :math:`L = 1/\sin\gamma`, the z-axis "Lorentz factor stationary mode"
    entry. In this geometry :math:`\gamma` is the exit angle of the
    diffracted beam from the surface, so the factor diverges as the rod
    approaches the surface plane, where a stationary measurement carries no
    information about the rod profile.

    :param gamma: Out-of-plane detector angle, in radian.
    :returns: The Lorentz factor, broadcast over the input.
    :rtype: numpy.ndarray
    """
    return np.abs(1.0 / np.sin(np.asarray(gamma, dtype=float)))


def rod_interception(gamma):
    r"""Rod interception factor, :math:`\cos\gamma`.

    :param gamma: Out-of-plane detector angle, in radian.
    :returns: The rod interception factor, broadcast over the input.
    :rtype: numpy.ndarray
    """
    return np.cos(np.asarray(gamma, dtype=float))


def area_correction(delta):
    r"""Area correction ignoring footprint and sample size, :math:`1/\sin\delta`.

    Listed in the manual for completeness. orGUI does not apply it: the
    footprint corrections of :mod:`~.beamprofile` evaluate the beam profile
    and finite sample size numerically instead, which is the row the manual
    marks as "calculated numerically in ANA".

    :param delta: In-plane detector angle, in radian.
    :returns: The area correction, broadcast over the input.
    :rtype: numpy.ndarray
    """
    return np.abs(1.0 / np.sin(np.asarray(delta, dtype=float)))


def lorentz_factor(mode, alpha=None, delta=None, gamma=None):
    r"""Lorentz factor for the given measurement mode.

    Dispatches to the three alternatives of the z-axis table. They are
    alternatives, not factors to be combined: which one applies is decided by
    how the intensity was measured.

    :param str mode: :data:`ROCKING`, :data:`REFLECTIVITY_ROCKING` or
        :data:`STATIONARY`.
    :param alpha: Incidence angle in radian; required except for
        :data:`STATIONARY`.
    :param delta: In-plane detector angle in radian; required for
        :data:`ROCKING`.
    :param gamma: Out-of-plane detector angle in radian; required except for
        :data:`REFLECTIVITY_ROCKING`.
    :returns: The Lorentz factor, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If ``mode`` is unknown or a required angle is
        missing.
    """
    required = {
        ROCKING: ("delta", "alpha", "gamma"),
        REFLECTIVITY_ROCKING: ("alpha",),
        STATIONARY: ("gamma",),
    }
    if mode not in required:
        raise ValueError(
            f"unknown Lorentz mode {mode!r}; expected one of "
            f"{ROCKING!r}, {REFLECTIVITY_ROCKING!r} or {STATIONARY!r}"
        )
    given = {"alpha": alpha, "delta": delta, "gamma": gamma}
    missing = [name for name in required[mode] if given[name] is None]
    if missing:
        raise ValueError(
            f"the {mode!r} Lorentz factor needs {', '.join(required[mode])}; "
            f"missing {', '.join(missing)}"
        )
    if mode == ROCKING:
        return lorentz_rocking_scan(delta, alpha, gamma)
    if mode == REFLECTIVITY_ROCKING:
        return lorentz_reflectivity_rocking_scan(alpha)
    return lorentz_stationary(gamma)

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
r"""Per-pixel correction factors of a detector image.

The factors that depend on where a photon landed on the detector rather than
on the scan: the solid angle each pixel subtends and the polarization factor
at its scattering angle. Both are properties of the calibrated geometry, so
they are evaluated once and reused for every frame that shares it.

They are returned as a **multiplicative** factor, the reciprocal of the
physical quantities, because that is how a corrected intensity is formed:

.. math::

    I_\mathrm{corr} = I \cdot \frac{1}{\Omega\,P}

The rocking integration, the stationary integration and the reciprocal-space
reconstruction each built this array themselves; this is the one definition
they share.

.. warning::

    :meth:`~.DetectorCalibration.Detector2D_SXRD.polarizationArray` evaluates
    the polarization at the **calibrated** detector position. That is correct
    for a detector whose arm does not move -- every pixel already carries its
    own scattering angle -- but not for a scan that drives the arm, where it
    understates the correction badly: 3 % at a scattering angle of 10
    degrees, 10 % at 18 and 33 % at 30.
    :meth:`~.DetectorCalibration.Detector2D_SXRD.polarizationAtPoints`
    follows the arm and is what such a scan needs. This function reproduces
    the historical, arm-blind behavior; see
    ``doc/design/ctr_structure_factor_scale.md`` finding F5.
"""

import numpy as np

__all__ = ["pixel_factors"]


def pixel_factors(detector, solid_angle=False, polarization=False, shape=None):
    r"""Per-pixel multiplicative correction factor of a detector image.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`
        with its calibrated geometry and polarization set.
    :param bool solid_angle: Divide by the solid angle each pixel subtends.
    :param bool polarization: Divide by the polarization factor.
    :param shape: Detector shape; taken from the detector when omitted.
    :returns: The factor for every pixel, or ``None`` when neither
        correction is enabled -- which lets a caller skip the multiplication
        entirely rather than multiply by an array of ones.
    :rtype: numpy.ndarray or None
    """
    factor = None
    if solid_angle:
        factor = 1.0 / np.asarray(
            detector.solidAngleArray(shape) if shape is not None
            else detector.solidAngleArray(),
            dtype=np.float64,
        )
    if polarization:
        values = np.asarray(
            detector.polarizationArray(shape) if shape is not None
            else detector.polarizationArray(),
            dtype=np.float64,
        )
        factor = 1.0 / values if factor is None else factor / values
    return factor

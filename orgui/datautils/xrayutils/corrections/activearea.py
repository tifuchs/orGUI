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
r"""The illuminated active surface area :math:`A`, in square meter.

:math:`A` is the surface area that actually contributes to the measured
counts. It multiplies the integrated intensity of every scan mode -- Vlieg
equation 42 for a rocking scan, 54 for a stationary one, 62 for reflectivity
-- so it sets the absolute scale of :math:`|F_{hkl}|^2` but cancels from the
ratio of two modes measured at the same incidence angle.

Which of two limits applies is an experimental question, not a preference:

**Slit limited.** Post-sample slits narrower than the footprint cut the
visible length along the surface, and the area is Vlieg equation 37,
:math:`A = s_1 s_2 / (\sin\delta\,\cos(\alpha-\beta_\mathrm{in}))`. It varies
with the in-plane detector angle and therefore along a rod.
:func:`slit_limited_area`.

**Beam and sample limited.** With an area detector and open slits -- orGUI's
usual configuration -- the whole illuminated footprint is seen and the area
is set by the beam profile and the sample, with no :math:`\delta` dependence.
:func:`beam_limited_area`.

Only the second is two-dimensional in the sense that matters: the beam's
*vertical* extent is what the projection onto the surface stretches, and it
is described in full by a
:class:`~orgui.datautils.xrayutils.corrections.beamprofile.BeamProfile`. The
horizontal extent is not projected and enters as a plain width, which is the
one number needed to turn the dimensionless factor a beam profile returns
into an absolute area:

.. math::

    A(\alpha) = w \, L \, C_\mathrm{illum}(\alpha, L)

with :math:`w` the horizontal beam width, :math:`L` the sample length along
the beam and :math:`C_\mathrm{illum}` the profile's
:meth:`~.beamprofile.BeamProfile.illuminated_area_fraction`, which is already
the mean of :math:`p(z)/p_\mathrm{max}` over the projected sample footprint.

.. note::

    The familiar closed form :math:`A = w \min(L, h/\sin\alpha)` is exactly
    this expression for a
    :func:`~.beamprofile.top_hat_profile` of full width :math:`h`, and only
    for that profile. A Gaussian of the same width differs from it by up to
    19 % at grazing incidence and by 6.5 % once the beam is fully on the
    sample, so the closed form is not a general shortcut and this module does
    not offer one.

For the alternative total-incident-flux convention, this module also returns
the intercepted fraction :math:`f_\mathrm{hit}=f_zf_x` and the dimensionless
illumination divisor :math:`H=f_\mathrm{hit}/\sin\alpha`. These replace the
legacy density-times-area term; they are not additional corrections. The
horizontal fraction is required explicitly because a vertical profile cannot
determine it.
"""

import numpy as np

from . import geometry

__all__ = [
    "beam_limited_area",
    "footprint_length",
    "illumination_divisor",
    "intercepted_fraction",
    "slit_limited_area",
]


def _total_flux_inputs(alpha, sample_length, horizontal_fraction, profile):
    """Validate and return inputs shared by total-flux overlap helpers."""
    if profile is None:
        raise ValueError("total-flux illumination needs a beam profile")
    if not float(sample_length) > 0 or not np.isfinite(sample_length):
        raise ValueError("sample_length must be finite and positive, in meter")
    alpha = np.asarray(alpha, dtype=np.float64)
    if (
        np.any(alpha <= 0)
        or np.any(alpha > np.pi / 2)
        or not np.all(np.isfinite(alpha))
    ):
        raise ValueError(
            "alpha must be finite and in the physical interval (0, pi/2] rad"
        )
    horizontal = np.asarray(horizontal_fraction, dtype=np.float64)
    if (
        np.any(horizontal <= 0)
        or np.any(horizontal > 1)
        or not np.all(np.isfinite(horizontal))
    ):
        raise ValueError("horizontal_fraction must be finite and in (0, 1]")
    return alpha, horizontal


def intercepted_fraction(
    alpha,
    sample_length,
    profile,
    *,
    horizontal_fraction,
):
    r"""Fraction :math:`f_\mathrm{hit}=f_z f_x` of total flux on sample.

    The vertical fraction ``f_z`` is the integral of ``profile`` over the
    sample's projected interval. ``horizontal_fraction`` supplies the
    independently known horizontal interception; a one-dimensional vertical
    profile cannot infer it.

    :param alpha: Incidence angle(s), in radian and in ``(0, pi/2]``.
    :param float sample_length: Sample length along the beam, in meter.
    :param profile: A :class:`~.beamprofile.BeamProfile`.
    :param horizontal_fraction: Known horizontal intercepted fraction
        :math:`f_x`, in ``(0, 1]``.
    :returns: Dimensionless total intercepted fraction.
    :rtype: numpy.ndarray
    :raises ValueError: If the geometry, profile or fraction is invalid.
    """
    alpha, horizontal = _total_flux_inputs(
        alpha, sample_length, horizontal_fraction, profile
    )
    vertical = np.asarray(
        profile.flux_on_sample(alpha, float(sample_length)), dtype=np.float64
    )
    if (
        np.any(vertical < 0)
        or np.any(vertical > 1.0 + 1e-12)
        or not np.all(np.isfinite(vertical))
    ):
        raise ValueError("beam profile returned an invalid intercepted fraction")
    return vertical * horizontal


def illumination_divisor(
    alpha,
    sample_length,
    profile,
    *,
    horizontal_fraction,
):
    r"""Total-flux illumination divisor :math:`H=f_\mathrm{hit}/\sin\alpha`.

    This is the dimensionless geometrical divisor paired with a total incident
    frame fluence :math:`Q_f`. It is a replacement convention for the legacy
    density-times-active-area product, not an additional footprint factor.
    The profile evaluates ``f_z / sin(alpha)`` directly so the finite
    grazing-incidence limit is retained without silently clipping the angle.

    :param alpha: Incidence angle(s), in radian and in ``(0, pi/2]``.
    :param float sample_length: Sample length along the beam, in meter.
    :param profile: A :class:`~.beamprofile.BeamProfile`.
    :param horizontal_fraction: Known horizontal intercepted fraction
        :math:`f_x`, in ``(0, 1]``.
    :returns: Dimensionless :math:`H`, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If the geometry or overlap is nonphysical.
    """
    alpha, horizontal = _total_flux_inputs(
        alpha, sample_length, horizontal_fraction, profile
    )
    ratio = np.asarray(
        profile.flux_over_sine(alpha, float(sample_length)), dtype=np.float64
    )
    result = ratio * horizontal
    if np.any(result <= 0) or not np.all(np.isfinite(result)):
        raise ValueError("sample/beam overlap must be finite and positive")
    return result


def footprint_length(alpha, beam_height, sample_length):
    r"""Illuminated length along the beam, in meter.

    :math:`\min(L, h/\sin\alpha)`, the geometric footprint of a beam of
    uniform vertical size. It is the length :func:`beam_limited_area` would
    use for a top-hat beam, and is kept as a separate function because it is
    the quantity to report when describing an experiment, not because
    :func:`beam_limited_area` needs it.

    :param alpha: Incidence angle, in radian.
    :param float beam_height: Vertical beam size, in meter.
    :param float sample_length: Sample length along the beam, in meter.
    :returns: The footprint length in meter, broadcast over ``alpha``.
    :rtype: numpy.ndarray
    :raises ValueError: If a size is not positive.
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    for name, value in (
        ("beam_height", beam_height),
        ("sample_length", sample_length),
    ):
        if not float(value) > 0:
            raise ValueError(f"{name} must be positive, in meter, got {value!r}")
    sin_alpha = np.sin(alpha)
    projected = np.divide(
        float(beam_height),
        sin_alpha,
        out=np.full(sin_alpha.shape, np.inf),
        where=sin_alpha > 0,
    )
    return np.minimum(float(sample_length), projected)


def beam_limited_area(alpha, beam_width, sample_length, profile):
    r"""Active area of an open-slit measurement, in square meter.

    :math:`A = w\,L\,C_\mathrm{illum}(\alpha, L)`. The vertical direction
    comes from the beam profile, which already integrates :math:`p(z)` over
    the projected sample and normalizes to the peak density; the horizontal
    direction is the plain beam width, since it is not projected.

    :param alpha: Incidence angle, in radian.
    :param float beam_width: Horizontal beam size, in meter.
    :param float sample_length: Sample length along the beam, in meter.
    :param profile: A
        :class:`~.beamprofile.BeamProfile` describing the vertical beam
        profile.
    :returns: The active area in square meter, broadcast over ``alpha``.
    :rtype: numpy.ndarray
    :raises ValueError: If a size is not positive, or no profile is given.
    """
    if profile is None:
        raise ValueError(
            "the beam-limited area needs a beam profile; a uniform beam of "
            "full width h is beamprofile.top_hat_profile(h)"
        )
    for name, value in (
        ("beam_width", beam_width),
        ("sample_length", sample_length),
    ):
        if not float(value) > 0:
            raise ValueError(f"{name} must be positive, in meter, got {value!r}")
    fraction = profile.illuminated_area_fraction(alpha, float(sample_length))
    return float(beam_width) * float(sample_length) * np.asarray(fraction)


def slit_limited_area(delta, slit_width, slit_height, alpha=None, beta_in=None):
    r"""Active area seen through post-sample slits, in square meter.

    Vlieg equation 37 with equation 38,
    :math:`A = s_1 s_2 / (\sin\delta \cos(\alpha - \beta_\mathrm{in}))`. The
    :math:`1/\sin\delta` is :func:`~.geometry.area_correction`, the row the
    ANA/ROD z-axis table lists; the remaining cosine is one in the z-axis
    mode, where the incidence angle is the angle to the surface
    (:math:`\alpha = \beta_\mathrm{in}`), which is why the table omits it.

    :param delta: In-plane detector angle, in radian.
    :param float slit_width: Beam size across the scattering plane
        (Vlieg's :math:`s_1`), in meter.
    :param float slit_height: Detector slit opening projected onto the
        surface (Vlieg's :math:`s_2`), in meter.
    :param alpha: Incidence angle, in radian. Only needed together with
        ``beta_in``.
    :param beta_in: Angle of the incident beam to the surface, in radian.
        ``None`` means the z-axis mode, ``beta_in = alpha``, and the cosine
        drops out.
    :returns: The active area in square meter, broadcast over the inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If a slit size is not positive, or ``beta_in`` is
        given without ``alpha``.
    """
    for name, value in (
        ("slit_width", slit_width),
        ("slit_height", slit_height),
    ):
        if not float(value) > 0:
            raise ValueError(f"{name} must be positive, in meter, got {value!r}")
    tilt = 1.0
    if beta_in is not None:
        if alpha is None:
            raise ValueError("beta_in needs alpha to form alpha - beta_in")
        tilt = np.cos(
            np.asarray(alpha, dtype=np.float64)
            - np.asarray(beta_in, dtype=np.float64)
        )
    projected = geometry.area_correction(delta) / tilt
    return float(slit_width) * float(slit_height) * projected

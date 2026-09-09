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
r"""Incident-beam profiles for numerical footprint corrections.

Two related quantities are evaluated at incidence angle :math:`\alpha` on a
sample of length :math:`L` along the beam:

``C_flux_on_sample``
    the fraction of the total incident flux that actually strikes the sample.
    This is the overlap integral used below and is stored as a diagnostic; it
    is not an additional intensity divisor.

``C_illum_area``
    the numerical active surface area, referenced to a sample fully bathed in
    a beam of the profile's peak intensity. This is the footprint divisor
    applied to an integrated intensity.

Both follow from the normalized vertical beam profile :math:`p(z)`, with
:math:`\int p(z)\,\mathrm{d}z = 1`, where :math:`z` runs perpendicular to the
beam in the scattering plane. With the projected sample size

.. math:: h(\alpha) = L \sin\alpha

and the sample centered at :math:`z_0` in the beam,

.. math::

    C_\mathrm{flux} = \int_{z_0 - h/2}^{z_0 + h/2} p(z)\,\mathrm{d}z
    \qquad
    C_\mathrm{area} = \frac{C_\mathrm{flux}}{p_\mathrm{max}\, h}

:math:`C_\mathrm{area}` is the mean of :math:`p(z)/p_\mathrm{max}` over the
projected sample footprint. It is the one-dimensional form of Vlieg's
numerically illuminated area, so it already contains the overlap integral
:math:`C_\mathrm{flux}`. Multiplying the two would count beam overspill twice.
It tends to 1 when the sample is small compared with a centered beam and falls
off as :math:`1/\sin\alpha` once the beam is fully on the sample.

For a Gaussian :math:`p` both integrals have the closed form orGUI used before
this module existed, reproduced exactly by :class:`GaussianBeamProfile`.
:class:`MeasuredBeamProfile` evaluates the same two definitions by numerical
integration of a tabulated profile, so an asymmetric or multiply-peaked beam
measured at the beamline can be used instead of the Gaussian idealization.

All lengths in this module are in **meters**, all angles in **radians**.
"""

from abc import ABC, abstractmethod

import numpy as np
from scipy import optimize, special, stats

__all__ = [
    "BeamProfile",
    "DistributionBeamProfile",
    "GaussianBeamProfile",
    "MeasuredBeamProfile",
    "gaussian_profile",
    "generalized_normal_profile",
    "profile_from_height_scan",
    "read_profile_file",
    "skew_normal_profile",
    "smoothed_top_hat_profile",
    "top_hat_profile",
    "trapezoid_profile",
    "triangular_profile",
    "trim_to_illuminated_edge",
]

#: Quantile at which a distribution's plotted and searched range is cut.
_TAIL = 1e-6

#: Half-width of that range in interquartile ranges, used as a second bound.
#: A heavy-tailed distribution puts its ``_TAIL`` quantile absurdly far out
#: -- a Cauchy profile of 0.3 mm FWHM reaches its at about 40 m -- which
#: would leave the peak unresolved by any practical number of samples. The
#: interquartile range stays finite for every distribution, and 4 of them
#: cover a Gaussian to 5.4 sigma, so the quantile bound still wins for
#: light-tailed and bounded profiles and nothing changes for them.
_RANGE_IQR = 4.0

#: Number of samples used to locate maxima and half-maximum crossings.
_SCAN_POINTS = 4001

if hasattr(np, "trapezoid"):  # ToDo remove for orGUI release >1.5
    _trapz_impl = np.trapezoid  # numpy >= 2.0
else:
    _trapz_impl = np.trapz  # noqa: NPY201  # numpy < 2.0

#: Conversion of a Gaussian FWHM to its standard deviation.
_FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


class BeamProfile(ABC):
    """Vertical intensity profile of the incident beam.

    Subclasses provide the applied active-area factor and its intercepted-flux
    diagnostic for arbitrarily shaped arrays of incidence angles.
    """

    @abstractmethod
    def flux_on_sample(self, alpha, L):
        """Fraction of the incident flux intercepted by the sample.

        :param alpha: Incidence angle(s) in radian, any array shape.
            Expected in ``(0, pi/2]``.
        :param float L: Sample size along the beam, in meters.
        :returns: ``C_flux_on_sample``, broadcast to the shape of ``alpha``.
        :rtype: numpy.ndarray
        """

    @abstractmethod
    def illuminated_area_fraction(self, alpha, L):
        """Illuminated fraction of the projected sample footprint.

        :param alpha: Incidence angle(s) in radian, any array shape.
        :param float L: Sample size along the beam, in meters.
        :returns: ``C_illum_area``, broadcast to the shape of ``alpha``.
        :rtype: numpy.ndarray
        """

    def corrections(self, alpha, L):
        """Return the intercepted flux and numerical active area.

        :param alpha: Incidence angle(s) in radian, any array shape.
        :param float L: Sample size along the beam, in meters.
        :returns: ``(C_flux_on_sample, C_illum_area)``. Only
            ``C_illum_area`` is an integrated-intensity divisor;
            ``C_flux_on_sample`` is its diagnostic numerator.
        :rtype: tuple
        """
        return self.flux_on_sample(alpha, L), self.illuminated_area_fraction(alpha, L)

    @abstractmethod
    def profile_curve(self, n=512):
        """Sample the profile for display.

        :param int n: Requested number of samples. A tabulated profile
            returns its own points and ignores this.
        :returns: ``(z, p)`` with ``z`` relative to the sample center in
            meters and ``p`` the normalized profile in 1/meter.
        :rtype: tuple of numpy.ndarray
        """

    @property
    @abstractmethod
    def centroid_position(self):
        """Center of mass relative to the sample center, in meters.

        ``0.0`` when the sample is centered on the centroid, and ``nan``
        for a profile whose first moment does not converge.

        :rtype: float
        """


def _half_maximum_width(z, p, pmax):
    """Width between the outermost crossings of ``pmax / 2``.

    Uses the outermost crossings so that a profile with several maxima
    reports its full extent rather than the width of a single sub-peak.

    :param z: Sample positions, increasing.
    :param p: Profile values at ``z``.
    :param float pmax: Peak value the half maximum is taken from.
    :returns: The full width at half maximum, in the units of ``z``.
    :rtype: float
    """
    half = pmax / 2.0
    above = np.flatnonzero(p >= half)
    if above.size == 0:
        return float("nan")
    lo, hi = int(above[0]), int(above[-1])
    if lo == 0:
        left = z[lo]
    else:
        left = np.interp(half, [p[lo - 1], p[lo]], [z[lo - 1], z[lo]])
    if hi == z.size - 1:
        right = z[hi]
    else:
        right = np.interp(half, [p[hi + 1], p[hi]], [z[hi + 1], z[hi]])
    return float(right - left)


class _CenteredProfile(BeamProfile):
    """Shared centering and correction evaluation for profiles with a CDF.

    Both corrections follow from the cumulative integral of the profile and
    from its peak density, so a subclass only has to supply
    ``_cumulative(x)`` and ``_density_at(x)`` -- both in coordinates relative
    to the sample center -- together with ``_pmax`` and the ``_tiny``
    threshold below which the ``alpha -> 0`` limit is used.
    """

    def _set_center(self, center, offset, centroid, peak_position, median):
        """Resolve the requested centering into :attr:`sample_center`.

        :param str center: ``"centroid"``, ``"peak"`` or ``"median"``.
        :param float offset: Extra displacement of the sample center.
        :param float centroid: Center of mass of the profile.
        :param float peak_position: Position of the maximum density.
        :param float median: Position at which half the flux has passed.
        :raises ValueError: If ``center`` is unknown, or names a reference
            point this profile does not have (a heavy-tailed distribution
            has no finite centroid).
        """
        try:
            reference = {
                "centroid": centroid,
                "peak": peak_position,
                "median": median,
            }[center]
        except KeyError:
            raise ValueError(
                f"center must be 'centroid', 'peak' or 'median', got {center!r}"
            ) from None
        if not np.isfinite(reference):
            raise ValueError(
                f"the {center} of this beam profile is not finite, so the "
                "sample cannot be centered on it. Heavy-tailed profiles have "
                "no center of mass; center on 'median' or 'peak' instead."
            )
        self.centroid = float(centroid)
        self.peak_position = float(peak_position)
        self.median = float(median)
        self.center = center
        self.offset = float(offset)
        #: Position of the sample center in the profile's own coordinate.
        self.sample_center = float(reference) + self.offset

    @property
    def centroid_position(self):
        """Center of mass relative to the sample center, in meters."""
        return float(self.centroid - self.sample_center)

    def _flux(self, h):
        """``C_flux_on_sample`` for a projected sample size ``h``."""
        return self._cumulative(h / 2.0) - self._cumulative(-h / 2.0)

    def flux_on_sample(self, alpha, L):
        """Fraction of the incident flux intercepted by the sample.

        Integrates the profile over the projected sample size
        ``h = L sin(alpha)``, centered on the sample center.
        """
        return self._flux(L * np.sin(np.asarray(alpha, dtype=float)))

    def illuminated_area_fraction(self, alpha, L):
        """Illuminated fraction of the projected sample footprint.

        ``C_flux_on_sample / (p_max * h)``, continued to its limit
        ``p(0) / p_max`` as ``h -> 0``.
        """
        h = L * np.sin(np.asarray(alpha, dtype=float))
        flux = self._flux(h)
        limit = self._density_at(0.0) / self._pmax
        denom = self._pmax * h
        return np.where(
            np.abs(h) > self._tiny,
            flux / np.where(denom == 0.0, 1.0, denom),
            limit,
        )


class GaussianBeamProfile(BeamProfile):
    """Analytical Gaussian beam profile.

    Reproduces orGUI's original closed-form footprint correction exactly; it
    is the reference against which :class:`MeasuredBeamProfile` is validated.

    :param float fwhm: Full width at half maximum of the vertical beam
        profile, in meters.
    :raises ValueError: If ``fwhm`` is not positive.
    """

    def __init__(self, fwhm):
        fwhm = float(fwhm)
        if not fwhm > 0:
            raise ValueError(f"beam FWHM must be positive, got {fwhm:g}")
        self.fwhm = fwhm
        self.sigma = fwhm * _FWHM_TO_SIGMA

    def flux_on_sample(self, alpha, L):
        """Fraction of the incident flux intercepted by the sample.

        Evaluates ``erf(h / (2 sqrt(2) sigma))`` with ``h = L sin(alpha)``.
        """
        arg = ((L * np.sin(alpha)) / (np.sqrt(2) * self.sigma)) * 0.5
        return (1 / 2) * (special.erf(arg) - special.erf(-arg))

    def illuminated_area_fraction(self, alpha, L):
        """Illuminated fraction of the projected sample footprint."""
        return (np.sqrt(2 * np.pi) * self.sigma * self.flux_on_sample(alpha, L)) / (
            L * np.sin(alpha)
        )

    def profile_curve(self, n=512):
        """Sample the Gaussian over +- 4 sigma around the sample center."""
        z = np.linspace(-4.0 * self.sigma, 4.0 * self.sigma, int(n))
        return z, np.exp(-0.5 * (z / self.sigma) ** 2) / (
            np.sqrt(2 * np.pi) * self.sigma
        )

    @property
    def centroid_position(self):
        """A Gaussian is symmetric, so its centroid is the sample center."""
        return 0.0

    def __repr__(self):
        return f"GaussianBeamProfile(fwhm={self.fwhm:g} m)"


class MeasuredBeamProfile(_CenteredProfile):
    """Tabulated beam profile integrated numerically.

    The profile is treated as piecewise linear between the supplied sample
    points and as exactly zero outside their range, so its cumulative
    integral -- and with it ``C_flux_on_sample`` -- is evaluated in closed
    form on each interval instead of being re-sampled onto an auxiliary grid.

    A measured profile carries no absolute position information, so it is
    re-referenced to the point of the beam that the center of the sample is
    aligned to; see ``center``. For a symmetric profile every choice
    coincides, and the results then reduce to :class:`GaussianBeamProfile`.

    :param z: Positions perpendicular to the beam in the scattering plane,
        in meters. Must be strictly monotonic; the direction is preserved
        (increasing ``z`` upward, as in a sample height scan).
    :param intensity: Beam intensity at ``z``, in arbitrary units. It is
        normalized internally, so any scale or monitor normalization is fine.
    :param str center: Which point of the profile the center of the sample
        sits on: ``"centroid"`` (center of mass, the default), ``"peak"``
        (maximum intensity), or ``"median"`` (the half-cut position an
        edge-scan alignment converges to).
    :param float offset: Additional displacement of the sample center from
        the ``center`` reference, in meters. Positive moves the sample
        center toward larger ``z``.
    :raises ValueError: If ``z`` is not strictly monotonic, if fewer than
        two points are given, if any value is not finite, or if the profile
        does not have a positive integral.
    """

    def __init__(self, z, intensity, center="centroid", offset=0.0):
        z = np.asarray(z, dtype=float).ravel()
        intensity = np.asarray(intensity, dtype=float).ravel()
        if z.size != intensity.size:
            raise ValueError(
                "z and intensity must have the same length, got "
                f"{z.size:d} and {intensity.size:d}"
            )
        if z.size < 2:
            raise ValueError("a beam profile needs at least two points")
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(intensity)):
            raise ValueError("beam profile contains non-finite values")

        dz = np.diff(z)
        if np.all(dz < 0):
            z, intensity = z[::-1], intensity[::-1]
            dz = np.diff(z)
        if not np.all(dz > 0):
            raise ValueError("z must be strictly monotonic")

        norm = _trapz_impl(intensity, z)
        if not norm > 0:
            raise ValueError(
                f"beam profile must have a positive integral, got {norm:g}. A "
                "height scan differentiated with the wrong sign gives a "
                "negative one."
            )

        p = intensity / norm
        # Exact cumulative integral of the piecewise-linear profile at the
        # sample points; interior positions are handled in _cumulative().
        cum = np.concatenate(([0.0], np.cumsum(0.5 * (p[:-1] + p[1:]) * dz)))
        cum[-1] = 1.0

        self._z_raw = z
        self._p = p
        self._cum = cum
        self._pmax = float(p.max())
        # Guard for the alpha -> 0 limit, scaled to the profile's support so
        # that it stays meaningful whatever length scale the caller works on.
        self._tiny = 1e-9 * float(z[-1] - z[0])

        centroid = float(_trapz_impl(z * p, z))
        self._set_center(
            center,
            offset,
            centroid,
            float(z[np.argmax(p)]),
            float(np.interp(0.5, cum, z)),
        )
        self.rms_width = float(
            np.sqrt(max(float(_trapz_impl((z - centroid) ** 2 * p, z)), 0.0))
        )

    def profile_curve(self, n=512):
        """Return the tabulated points; ``n`` is ignored."""
        return self.z, self._p

    @property
    def z(self):
        """Profile positions relative to the sample center, in meters."""
        return self._z_raw - self.sample_center

    @property
    def density(self):
        """Normalized profile ``p(z)`` in 1/meter, matching :attr:`z`."""
        return self._p

    @property
    def support(self):
        """``(z_min, z_max)`` of the profile relative to the sample center."""
        return (
            float(self._z_raw[0] - self.sample_center),
            float(self._z_raw[-1] - self.sample_center),
        )

    @property
    def fwhm(self):
        """Full width at half maximum of the profile, in meters.

        Measured between the outermost crossings of half the peak value, so
        a profile with several maxima reports its full extent rather than
        the width of a single sub-peak.
        """
        return _half_maximum_width(self._z_raw, self._p, self._pmax)

    def _density_at(self, x):
        """Normalized profile at ``x``, relative to the sample center."""
        return np.interp(
            np.asarray(x, dtype=float) + self.sample_center, self._z_raw, self._p
        )

    def _cumulative(self, x):
        """Integral of the normalized profile from ``-inf`` up to ``x``.

        ``x`` is given relative to the sample center. Outside the tabulated
        range the profile is taken to be zero, so the result saturates at 0
        and 1.
        """
        z, p, cum = self._z_raw, self._p, self._cum
        x = np.clip(np.asarray(x, dtype=float) + self.sample_center, z[0], z[-1])
        i = np.clip(np.searchsorted(z, x, side="right") - 1, 0, z.size - 2)
        d = z[i + 1] - z[i]
        t = (x - z[i]) / d
        # Exact integral of one linear segment over [z[i], x].
        return cum[i] + d * t * (p[i] + 0.5 * t * (p[i + 1] - p[i]))

    def __repr__(self):
        return (
            f"MeasuredBeamProfile({self._z_raw.size:d} points, "
            f"fwhm={self.fwhm:g} m, center={self.center!r})"
        )


class DistributionBeamProfile(_CenteredProfile):
    """Beam profile given by an analytical probability distribution.

    Both corrections need only the cumulative distribution and the peak
    density, so any continuous distribution with a ``cdf``, a ``pdf`` and a
    ``ppf`` can describe the beam -- in particular any frozen
    :mod:`scipy.stats` distribution. The named constructors in this module
    (:func:`top_hat_profile`, :func:`trapezoid_profile`,
    :func:`smoothed_top_hat_profile`, :func:`generalized_normal_profile`,
    :func:`skew_normal_profile`) wrap the shapes that describe real beams,
    in physical parameters rather than raw distribution arguments.

    :class:`GaussianBeamProfile` is the same model evaluated in closed form;
    passing ``scipy.stats.norm`` here reproduces it to machine precision.

    ``scipy.stats`` exposes no mode, so the peak is located once at
    construction by scanning the distribution's central range and refining
    the best sample. A plateau, as in a top hat, resolves to its midpoint.

    :param dist: A frozen continuous distribution, e.g.
        ``scipy.stats.norm(scale=1e-5)``.
    :param str center: Which point of the profile the center of the sample
        sits on: ``"centroid"`` (mean, the default), ``"peak"`` (mode) or
        ``"median"``. All three coincide for a symmetric distribution.
    :param float offset: Additional displacement of the sample center from
        the ``center`` reference, in meters.
    :raises TypeError: If ``dist`` is not a frozen distribution.
    :raises ValueError: If the distribution has no finite central range, or
        no finite positive maximum density -- a beam profile with an
        unbounded peak has no active-area reference and cannot be used.
    """

    def __init__(self, dist, center="centroid", offset=0.0):
        for name in ("cdf", "pdf", "ppf"):
            if not callable(getattr(dist, name, None)):
                raise TypeError(
                    "dist must be a frozen continuous distribution with cdf, "
                    f"pdf and ppf methods, got {dist!r}. Freeze a scipy.stats "
                    "distribution first, e.g. scipy.stats.norm(scale=1e-5)."
                )
        if isinstance(dist, stats.rv_continuous):
            # An unfrozen distribution answers every call with its standard
            # form, which here would silently describe a beam one meter wide.
            raise TypeError(
                f"dist must be frozen, got the distribution family {dist!r}. "
                "Supply its parameters, e.g. scipy.stats.norm(scale=1e-5) "
                "rather than scipy.stats.norm."
            )
        lo, hi = self._central_range(dist)

        self._dist = dist
        self._lo = lo
        self._hi = hi
        peak_position, pmax = self._locate_peak(dist, lo, hi)
        if not np.isfinite(pmax) or pmax <= 0:
            raise ValueError(
                f"the distribution has no finite positive maximum density "
                f"(got {pmax:g}); the illuminated-area correction has no "
                "reference intensity for such a profile."
            )
        self._reject_unbounded_density(dist, pmax)
        self._pmax = pmax
        self._tiny = 1e-9 * (hi - lo)

        self._set_center(
            center, offset, self._finite_moment(dist.mean), peak_position, dist.ppf(0.5)
        )
        self.rms_width = self._finite_moment(dist.std)

    @staticmethod
    def _central_range(dist):
        """Range over which the profile is scanned and displayed.

        The narrower of the ``_TAIL`` quantile range and ``_RANGE_IQR``
        interquartile ranges around the median. This bounds only where the
        peak is looked for, the width is measured and the preview is drawn;
        the corrections always integrate the full distribution, so a
        heavy-tailed profile still reports the flux that misses the sample.

        :returns: ``(lo, hi)``.
        :rtype: tuple of float
        :raises ValueError: If the distribution is not localized.
        """
        lo = float(dist.ppf(_TAIL))
        hi = float(dist.ppf(1.0 - _TAIL))
        median = float(dist.ppf(0.5))
        iqr = float(dist.ppf(0.75)) - float(dist.ppf(0.25))
        if np.isfinite(iqr) and iqr > 0:
            lo = max(lo, median - _RANGE_IQR * iqr)
            hi = min(hi, median + _RANGE_IQR * iqr)
        if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
            raise ValueError(
                f"the distribution has no finite central range, got {lo:g} to "
                f"{hi:g}. A beam profile must be localized."
            )
        return lo, hi

    @staticmethod
    def _reject_unbounded_density(dist, pmax):
        """Refuse a distribution whose density diverges at an end.

        The peak is searched between finite quantiles, so a density that
        runs to infinity -- as a gamma distribution with shape below one
        does at zero -- still yields a large but finite maximum there, and
        that value would depend on nothing but the chosen cut. Probing a
        much deeper quantile catches it: for a bounded density the value
        stops growing, for a divergent one it does not.

        :raises ValueError: If the density is still growing at either end.
        """
        for quantile in (_TAIL * 1e-3, 1.0 - _TAIL * 1e-3):
            edge = float(dist.ppf(quantile))
            if not np.isfinite(edge):
                continue
            density = float(dist.pdf(edge))
            if not np.isfinite(density) or density > pmax * 1.01:
                raise ValueError(
                    "the density of this distribution diverges, so its "
                    "maximum depends only on where the tail is cut and the "
                    "illuminated-area correction has no reference intensity. "
                    "Use a distribution with a bounded peak."
                )

    @staticmethod
    def _finite_moment(func):
        """Evaluate a distribution moment, returning ``nan`` if it diverges."""
        try:
            value = float(func())
        except (ValueError, TypeError, ZeroDivisionError):
            return float("nan")
        return value if np.isfinite(value) else float("nan")

    @staticmethod
    def _locate_peak(dist, lo, hi):
        """Find the position and value of the maximum density.

        Scans the central range first so that a plateau or a second maximum
        cannot be missed, then refines a single interior maximum. The
        plateau of a top hat resolves to its midpoint.

        :returns: ``(peak_position, peak_density)``.
        :rtype: tuple of float
        """
        z = np.linspace(lo, hi, _SCAN_POINTS)
        p = np.asarray(dist.pdf(z), dtype=float)
        p = np.where(np.isfinite(p), p, -np.inf)
        pmax = float(p.max())
        if not np.isfinite(pmax) or pmax <= 0:
            return float("nan"), pmax
        # A flat top has no single argmax; take the middle of the plateau.
        flat = np.flatnonzero(p >= pmax * (1.0 - 1e-12))
        if flat.size > 1:
            return float(0.5 * (z[flat[0]] + z[flat[-1]])), pmax
        best = int(flat[0])
        if best == 0 or best == z.size - 1:
            return float(z[best]), pmax
        refined = optimize.minimize_scalar(
            lambda x: -float(dist.pdf(x)),
            bounds=(z[best - 1], z[best + 1]),
            method="bounded",
            options={"xatol": (hi - lo) * 1e-12},
        )
        position = float(refined.x)
        density = float(dist.pdf(position))
        return (position, density) if density >= pmax else (float(z[best]), pmax)

    def _cumulative(self, x):
        """Integral of the profile up to ``x``, relative to the sample center."""
        return self._dist.cdf(np.asarray(x, dtype=float) + self.sample_center)

    def _density_at(self, x):
        """Profile density at ``x``, relative to the sample center."""
        return self._dist.pdf(np.asarray(x, dtype=float) + self.sample_center)

    def profile_curve(self, n=512):
        """Sample the distribution over its central range."""
        z = np.linspace(self._lo, self._hi, int(n))
        return z - self.sample_center, np.asarray(self._dist.pdf(z), dtype=float)

    @property
    def support(self):
        """``(z_min, z_max)`` of the plotted range, relative to the center."""
        return (
            float(self._lo - self.sample_center),
            float(self._hi - self.sample_center),
        )

    @property
    def fwhm(self):
        """Full width at half maximum of the distribution, in meters."""
        z = np.linspace(self._lo, self._hi, _SCAN_POINTS)
        return _half_maximum_width(
            z, np.asarray(self._dist.pdf(z), dtype=float), self._pmax
        )

    def __repr__(self):
        return (
            f"DistributionBeamProfile({self._dist!r}, "
            f"fwhm={self.fwhm:g} m, center={self.center!r})"
        )


class _SmoothedTopHat:
    """Top hat of full ``width`` convolved with a Gaussian of ``sigma``.

    The frozen-distribution interface :class:`DistributionBeamProfile`
    needs, implemented in closed form. This is the beam a pair of defining
    slits produces once the finite source size blurs the slit edges, so the
    profile is flat in the middle with error-function flanks.

    :param float width: Full width of the top hat, in meters.
    :param float sigma: Standard deviation of the edge blur, in meters.
    """

    def __init__(self, width, sigma):
        width = float(width)
        sigma = float(sigma)
        if not width > 0:
            raise ValueError(f"width must be positive, got {width:g}")
        if not sigma > 0:
            raise ValueError(f"edge sigma must be positive, got {sigma:g}")
        self.width = width
        self.sigma = sigma
        self._a = -0.5 * width
        self._b = 0.5 * width

    def pdf(self, x):
        """Density: the difference of two shifted normal CDFs."""
        x = np.asarray(x, dtype=float)
        return (
            special.ndtr((x - self._a) / self.sigma)
            - special.ndtr((x - self._b) / self.sigma)
        ) / self.width

    def cdf(self, x):
        """Cumulative distribution, the antiderivative of :meth:`pdf`."""
        x = np.asarray(x, dtype=float)

        def term(edge):
            u = (x - edge) / self.sigma
            return (x - edge) * special.ndtr(u) + self.sigma * np.exp(
                -0.5 * u**2
            ) / np.sqrt(2 * np.pi)

        return (term(self._a) - term(self._b)) / self.width

    def ppf(self, q):
        """Quantile function, inverted numerically from :meth:`cdf`."""
        q = float(q)
        if q <= 0.0 or q >= 1.0:
            raise ValueError(f"quantile must be in (0, 1), got {q:g}")
        span = 0.5 * self.width + 12.0 * self.sigma
        return float(optimize.brentq(lambda x: self.cdf(x) - q, -span, span))

    def mean(self):
        """The profile is symmetric about zero."""
        return 0.0

    def std(self):
        """Quadrature sum of the top-hat and blur widths."""
        return float(np.sqrt(self.width**2 / 12.0 + self.sigma**2))

    def __repr__(self):
        return f"_SmoothedTopHat(width={self.width:g}, sigma={self.sigma:g})"


def gaussian_profile(fwhm, center="centroid", offset=0.0):
    """Gaussian beam profile that supports centering and an offset.

    Numerically the same model as :class:`GaussianBeamProfile`, which
    evaluates it in closed form but always sits centered on the sample. Use
    this one when the sample is displaced from the beam center.

    :param float fwhm: Full width at half maximum, in meters.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    :raises ValueError: If ``fwhm`` is not positive.
    """
    fwhm = float(fwhm)
    if not fwhm > 0:
        raise ValueError(f"fwhm must be positive, got {fwhm:g}")
    return DistributionBeamProfile(
        stats.norm(loc=0.0, scale=fwhm * _FWHM_TO_SIGMA), center=center, offset=offset
    )


def top_hat_profile(width, center="centroid", offset=0.0):
    """Beam profile of uniform intensity over ``width``.

    The slit-limited beam. With this profile ``C_flux_on_sample`` is
    ``L sin(alpha) / width`` until the sample intercepts the whole beam,
    the linear correction of Gibaud, Vignaud & Sinha (1993),
    *Acta Cryst.* A49, 642, equation (12), and ``C_illum_area`` stays at 1
    over the same range.

    :param float width: Full width of the beam, in meters.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    :raises ValueError: If ``width`` is not positive.
    """
    width = float(width)
    if not width > 0:
        raise ValueError(f"width must be positive, got {width:g}")
    return DistributionBeamProfile(
        stats.uniform(loc=-0.5 * width, scale=width), center=center, offset=offset
    )


def trapezoid_profile(full_width, flat_width, center="centroid", offset=0.0):
    """Beam profile of a trapezoid with symmetric flanks.

    Two slits of different apertures convolve to a trapezoid whose flat top
    is the smaller aperture and whose base is the larger one; matched slits
    give the triangle of :func:`triangular_profile`.

    :param float full_width: Width at the base, in meters.
    :param float flat_width: Width of the flat top, in meters, between 0 and
        ``full_width``.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    :raises ValueError: If the widths are not positive and ordered.
    """
    full_width = float(full_width)
    flat_width = float(flat_width)
    if not full_width > 0:
        raise ValueError(f"full_width must be positive, got {full_width:g}")
    if not 0.0 <= flat_width <= full_width:
        raise ValueError(
            f"flat_width must be between 0 and full_width, got {flat_width:g} "
            f"with full_width {full_width:g}"
        )
    ramp = 0.5 * (1.0 - flat_width / full_width)
    return DistributionBeamProfile(
        stats.trapezoid(ramp, 1.0 - ramp, loc=-0.5 * full_width, scale=full_width),
        center=center,
        offset=offset,
    )


def triangular_profile(full_width, center="centroid", offset=0.0):
    """Beam profile of a symmetric triangle, from two matched slits.

    :param float full_width: Width at the base, in meters.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    """
    return trapezoid_profile(full_width, 0.0, center=center, offset=offset)


def smoothed_top_hat_profile(width, edge_sigma, center="centroid", offset=0.0):
    """Beam profile of a top hat with error-function flanks.

    A slit-defined beam whose edges are blurred by the finite source size:
    flat in the middle, with flanks of standard deviation ``edge_sigma``.

    :param float width: Full width of the flat part before blurring, in
        meters.
    :param float edge_sigma: Standard deviation of the edge blur, in meters.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    """
    return DistributionBeamProfile(
        _SmoothedTopHat(width, edge_sigma), center=center, offset=offset
    )


def generalized_normal_profile(fwhm, flatness, center="centroid", offset=0.0):
    """Beam profile interpolating between a peaked, Gaussian and flat beam.

    The generalized normal density is proportional to
    ``exp(-|z / scale| ** flatness)``: ``flatness = 1`` gives an exponential
    cusp, ``flatness = 2`` is exactly Gaussian, and large values approach a
    top hat. It is the one-parameter family for a focused beam that is
    neither Gaussian nor flat.

    :param float fwhm: Full width at half maximum, in meters.
    :param float flatness: Shape exponent, positive.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    :raises ValueError: If ``fwhm`` or ``flatness`` is not positive.
    """
    fwhm = float(fwhm)
    flatness = float(flatness)
    if not fwhm > 0:
        raise ValueError(f"fwhm must be positive, got {fwhm:g}")
    if not flatness > 0:
        raise ValueError(f"flatness must be positive, got {flatness:g}")
    # exp(-(w/2 / scale)**flatness) = 1/2 at the half maximum.
    scale = 0.5 * fwhm / np.log(2.0) ** (1.0 / flatness)
    return DistributionBeamProfile(
        stats.gennorm(flatness, loc=0.0, scale=scale), center=center, offset=offset
    )


def skew_normal_profile(fwhm, skew, center="centroid", offset=0.0):
    """Asymmetric beam profile with a Gaussian core.

    The skew-normal density, scaled so that its full width at half maximum
    is ``fwhm``. ``skew = 0`` is Gaussian; positive values put the tail
    toward larger ``z``.

    :param float fwhm: Full width at half maximum, in meters.
    :param float skew: Shape parameter; 0 is symmetric.
    :param str center: Centering, see :class:`DistributionBeamProfile`.
    :param float offset: Sample-center displacement, in meters.
    :rtype: DistributionBeamProfile
    :raises ValueError: If ``fwhm`` is not positive.
    """
    fwhm = float(fwhm)
    if not fwhm > 0:
        raise ValueError(f"fwhm must be positive, got {fwhm:g}")
    skew = float(skew)
    # The skew normal has no closed-form FWHM; measure it once at unit
    # scale and rescale, since the family is a location-scale family.
    unit = stats.skewnorm(skew)
    z = np.linspace(unit.ppf(_TAIL), unit.ppf(1.0 - _TAIL), _SCAN_POINTS)
    p = unit.pdf(z)
    unit_fwhm = _half_maximum_width(z, p, float(p.max()))
    return DistributionBeamProfile(
        stats.skewnorm(skew, loc=0.0, scale=fwhm / unit_fwhm),
        center=center,
        offset=offset,
    )


def profile_from_height_scan(z, intensity, monitor=None):
    r"""Differentiate a sample height scan into a beam profile.

    In a height (``samz``) scan the sample edge cuts progressively further
    into the beam, so the transmitted intensity is the beam profile
    integrated over the part of the beam still passing the sample,
    :math:`I(z) = I_0 \int_z^\infty p(z')\,\mathrm{d}z'`. The profile is
    therefore the *negated* derivative of the measured curve,

    .. math:: p(z) \propto -\frac{\mathrm{d}I}{\mathrm{d}z}

    which also removes any constant transmission offset. ``z`` keeps its
    direction, so the returned profile is the beam profile in the same
    (upward) coordinate as the scan.

    Only the range over which the edge cuts into the beam is meaningful. A
    scan that continues until the sample leaves the beam again has a
    negative derivative there and must be sliced before use, otherwise the
    result is a profile followed by its mirror image.

    :param z: Height positions of the scan, in meters, strictly monotonic.
    :param intensity: Transmitted intensity at each ``z``.
    :param monitor: Optional incident-beam monitor at each ``z``. When
        given, ``intensity`` is divided by it before differentiation, which
        removes storage-ring current drift over the scan. The mean monitor
        value is multiplied back in, so the profile keeps the scale of the
        raw counts.
    :returns: ``(z, profile)``, the profile in arbitrary units and of the
        same length as ``z``.
    :rtype: tuple of numpy.ndarray
    :raises ValueError: If the inputs have different lengths, if fewer than
        two points are given, or if ``monitor`` contains zeros.
    """
    z = np.asarray(z, dtype=float).ravel()
    intensity = np.asarray(intensity, dtype=float).ravel()
    if z.size != intensity.size:
        raise ValueError(
            "z and intensity must have the same length, got "
            f"{z.size:d} and {intensity.size:d}"
        )
    if z.size < 2:
        raise ValueError("a height scan needs at least two points")
    if monitor is not None:
        monitor = np.asarray(monitor, dtype=float).ravel()
        if monitor.size != z.size:
            raise ValueError(
                "monitor must have the same length as z, got "
                f"{monitor.size:d} and {z.size:d}"
            )
        if np.any(monitor == 0):
            raise ValueError("monitor contains zeros")
        intensity = intensity / monitor * float(np.mean(monitor))
    return z, -np.gradient(intensity, z)


def trim_to_illuminated_edge(z, profile):
    """Restrict a differentiated height scan to its single cutting edge.

    A height scan long enough to move the sample fully through the beam
    contains two edges: the sample cutting in, whose negated derivative is
    the beam profile, and -- if the beam clears the sample again -- the
    sample leaving, which appears as a negative excursion. Only the first is
    a beam profile.

    The returned slice is the run of samples around the maximum over which
    the profile stays positive, i.e. it is cut at the sign changes flanking
    the peak. A profile that is positive everywhere is returned unchanged.

    :param z: Positions, strictly monotonic and increasing.
    :param profile: Profile values at ``z``, as returned by
        :func:`profile_from_height_scan`.
    :returns: ``(z, profile)`` restricted to the edge.
    :rtype: tuple of numpy.ndarray
    :raises ValueError: If no positive sample exists at all, which means the
        scan was differentiated with the wrong sign.
    """
    z = np.asarray(z, dtype=float).ravel()
    profile = np.asarray(profile, dtype=float).ravel()
    if z.size != profile.size:
        raise ValueError(
            "z and profile must have the same length, got "
            f"{z.size:d} and {profile.size:d}"
        )
    if not np.any(profile > 0):
        raise ValueError(
            "the differentiated height scan has no positive values; it is "
            "most likely differentiated with the wrong sign"
        )
    peak = int(np.argmax(profile))
    nonpositive = np.flatnonzero(profile <= 0)
    before = nonpositive[nonpositive < peak]
    after = nonpositive[nonpositive > peak]
    start = int(before[-1]) + 1 if before.size else 0
    stop = int(after[0]) if after.size else z.size
    return z[start:stop], profile[start:stop]


def read_profile_file(
    path,
    z_scale=1e-3,
    height_scan=False,
    usecols=(0, 1),
    monitor_col=None,
    trim=True,
):
    """Read a two-column text file into a beam profile.

    The file is read with :func:`numpy.loadtxt`, so ``#`` comments and any
    common whitespace separator are accepted.

    :param str path: File to read.
    :param float z_scale: Factor converting the file's position column to
        meters. The default ``1e-3`` reads millimeters, the usual unit of a
        diffractometer height scan; use ``1e-6`` for micrometers.
    :param bool height_scan: If ``True``, the intensity column is a measured
        height scan and is differentiated by
        :func:`profile_from_height_scan`. If ``False`` (default), it already
        is a beam profile.
    :param usecols: ``(position, intensity)`` column indices.
    :param monitor_col: Optional column index of a monitor counter, only
        used together with ``height_scan``.
    :param bool trim: Restrict a differentiated height scan to its cutting
        edge with :func:`trim_to_illuminated_edge`. Ignored unless
        ``height_scan`` is set.
    :returns: ``(z, profile)`` with ``z`` in meters.
    :rtype: tuple of numpy.ndarray
    """
    cols = tuple(usecols) if monitor_col is None else tuple(usecols) + (monitor_col,)
    data = np.loadtxt(path, usecols=cols, unpack=True, ndmin=2)
    z = np.asarray(data[0], dtype=float) * float(z_scale)
    intensity = np.asarray(data[1], dtype=float)
    monitor = np.asarray(data[2], dtype=float) if monitor_col is not None else None
    if not height_scan:
        return z, intensity
    if z.size > 1 and z[1] < z[0]:  # keep the edge search on an increasing axis
        z, intensity = z[::-1], intensity[::-1]
        monitor = None if monitor is None else monitor[::-1]
    z, profile = profile_from_height_scan(z, intensity, monitor)
    if trim:
        z, profile = trim_to_illuminated_edge(z, profile)
    return z, profile

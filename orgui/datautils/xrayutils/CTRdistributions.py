"""Statistical occupancy profiles for CTR interfaces and surfaces."""

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy.stats import poisson, skellam


DEFAULT_TAIL_PROBABILITY = 1e-10
DEFAULT_MOMENT_TOLERANCE = 1e-10
SURFACE_OCCUPANCY_TOLERANCE = 1e-12


def _validate_tail_probability(tail_probability):
    tail_probability = float(tail_probability)
    if not 0.0 < tail_probability < 1.0:
        raise ValueError("tail_probability must be between zero and one")
    return tail_probability


class SurfaceProfile(ABC):
    """Contract for discrete surface-height occupancy profiles.

    Offsets passed to :meth:`surface_occupancy` are consecutive structural
    layers ordered from bottom to top.  :meth:`occupancy` returns the total
    material fraction at each offset.  The exposed fraction is the decrease
    in material occupancy toward the next layer, with the highest represented
    layer treated as completely exposed.
    """

    @abstractmethod
    def support(self):
        """Return the inclusive structural-layer support offsets."""

    @abstractmethod
    def occupancy(self, offsets):
        """Return cumulative material occupancy at structural-layer offsets."""

    def correction(self, offsets):
        """Return material occupancy relative to a sharp surface at zero."""
        offsets = np.asarray(offsets, dtype=np.float64)
        sharp = (offsets < 0).astype(np.float64)
        return self.occupancy(offsets) - sharp

    def surface_occupancy(self, offsets):
        """Return the truly exposed fraction of each represented layer.

        :param numpy.ndarray offsets:
            Consecutive integer structural-layer offsets ordered bottom to
            top.
        :returns:
            Exposed surface occupancy for every supplied offset.
        :rtype: numpy.ndarray
        :raises ValueError:
            If offsets are not consecutive or the cumulative occupancy is
            non-finite, outside zero to one, or increases toward the surface.
        """
        offsets = np.asarray(offsets, dtype=np.float64)
        if offsets.ndim != 1:
            raise ValueError("surface offsets must be a one-dimensional array")
        if offsets.size == 0:
            return np.array([], dtype=np.float64)
        if not np.all(np.isclose(offsets, np.rint(offsets))):
            raise ValueError("surface offsets must be integer structural layers")
        if offsets.size > 1 and not np.all(np.diff(offsets) == 1):
            raise ValueError("surface offsets must be consecutive and bottom-to-top")

        occupancy = np.asarray(self.occupancy(offsets), dtype=np.float64)
        if occupancy.shape != offsets.shape:
            raise ValueError("surface occupancy must match the offsets shape")
        if not np.all(np.isfinite(occupancy)):
            raise ValueError("surface occupancy must be finite")
        tolerance = SURFACE_OCCUPANCY_TOLERANCE
        if np.any(occupancy < -tolerance) or np.any(occupancy > 1.0 + tolerance):
            raise ValueError("surface occupancy must be between zero and one")
        occupancy = np.clip(occupancy, 0.0, 1.0)

        exposed = np.empty_like(occupancy)
        exposed[:-1] = occupancy[:-1] - occupancy[1:]
        exposed[-1] = occupancy[-1]
        if np.any(exposed < -tolerance):
            raise ValueError("surface occupancy must not increase toward the surface")
        return np.clip(exposed, 0.0, 1.0)


class PoissonProfile(SurfaceProfile):
    """Describe signed step--Poisson growth or etching of a surface.

    The random surface-height change is the convolution

    ``Step(offset + (1 - alpha) * mean_change)``

    ``+ sign(mean_change) * Poisson(alpha * abs(mean_change))``.

    For ``s = sign(s) * (floor(abs(s)) + f)``, ``Step(s)`` is
    ``sign(s) * (floor(abs(s)) + Bernoulli(f))``. This two-layer
    interpolation keeps a fractional deterministic offset continuous on the
    structural-layer grid. Consequently, ``alpha = 0`` gives ideal
    layer-by-layer growth or dissolution, and ``alpha = 1`` combines the
    deterministic structural cursor with a Poisson process.

    :param float mean_change:
        Signed mean process width in structural layers. Positive values model
        growth and negative values model dissolution or etching.
    :param float alpha:
        Fraction of the absolute mean assigned to the Poisson component.
        Must lie between zero and one.
    :param float offset:
        Deterministic surface-height offset in structural layers.
    :param float tail_probability:
        Maximum omitted probability mass in the unrepresented Poisson tail.
    :param float mean:
        Legacy alias for ``mean_change``.
    """

    def __init__(
        self,
        mean_change=None,
        alpha=1.0,
        offset=0.0,
        tail_probability=DEFAULT_TAIL_PROBABILITY,
        *,
        mean=None,
    ):
        if mean_change is None:
            if mean is None:
                raise TypeError("mean_change must be provided")
            mean_change = mean
        elif mean is not None:
            raise TypeError("provide mean_change or mean, not both")
        self.mean_change = float(mean_change)
        self.alpha = float(alpha)
        if not 0.0 <= self.alpha <= 1.0:
            raise ValueError("alpha must be between zero and one")
        self.offset = float(offset)
        self.tail_probability = _validate_tail_probability(tail_probability)

    @property
    def mean(self):
        """Return the legacy name for the signed mean height change."""
        return self.mean_change

    @property
    def rate(self):
        """Return the Poisson component mean in structural layers."""
        return self.alpha * abs(self.mean_change)

    @property
    def step_mean(self):
        """Return the unsigned process step mean in structural layers."""
        return (1.0 - self.alpha) * abs(self.mean_change)

    @property
    def structural_mean(self):
        """Return the signed deterministic cursor mean in structural layers."""
        return self.offset + (1.0 - self.alpha) * self.mean_change

    @property
    def expected_height_change(self):
        """Return the expected surface-height change in structural layers."""
        return self.offset + self.mean_change

    def _step_parameters(self):
        step_sign = -1 if self.structural_mean < 0 else 1
        magnitude = abs(self.structural_mean)
        step_floor = int(np.floor(magnitude))
        step_fraction = magnitude - step_floor
        step_low = step_sign * step_floor
        step_high = step_low + step_sign
        return step_low, step_high, step_fraction

    def _height_cdf(self, values):
        """Return the CDF of the signed convolved surface-height change."""
        values = np.asarray(values, dtype=np.float64)
        step_low, step_high, step_fraction = self._step_parameters()
        if self.rate == 0:
            lower = (values >= step_low).astype(np.float64)
            upper = (values >= step_high).astype(np.float64)
        elif self.mean_change > 0:
            lower = poisson.cdf(np.floor(values - step_low), self.rate)
            upper = poisson.cdf(np.floor(values - step_high), self.rate)
        else:
            lower = poisson.sf(np.ceil(step_low - values) - 1, self.rate)
            upper = poisson.sf(np.ceil(step_high - values) - 1, self.rate)
        return (1.0 - step_fraction) * lower + step_fraction * upper

    def probability(self, changes):
        """Return probability masses for signed integer height changes."""
        changes = np.asarray(changes, dtype=np.float64)
        integer = np.isclose(changes, np.rint(changes))
        changes = np.rint(changes)
        probability = self._height_cdf(changes) - self._height_cdf(changes - 1)
        return np.where(integer, probability, 0.0)

    def support(self):
        """Return the inclusive structural-layer support offsets."""
        step_low, step_high, step_fraction = self._step_parameters()
        if step_fraction == 0.0:
            step_high = step_low
        if self.rate == 0:
            process_low = min(step_low, step_high)
            process_high = max(step_low, step_high)
        else:
            quantile = int(poisson.ppf(1.0 - self.tail_probability, self.rate))
            if self.mean_change > 0:
                process_low = min(step_low, step_high)
                process_high = max(step_low, step_high) + quantile
            else:
                process_low = min(step_low, step_high) - quantile
                process_high = max(step_low, step_high)
        lower = int(np.floor(min(0.0, process_low))) - 1
        upper = int(np.ceil(max(0.0, process_high)))
        return lower, upper

    def occupancy(self, offsets):
        """Return material occupancy at structural-layer offsets."""
        offsets = np.asarray(offsets, dtype=np.float64)
        return 1.0 - self._height_cdf(offsets)

@dataclass(frozen=True)
class RecenteredSkellamProfile(SurfaceProfile):
    """Recenter a Skellam height distribution on integer structural layers.

    Heights are ``Poisson(k_up*t) - Poisson(k_down*t) + Step(s)``, where
    ``s = offset - (k_up-k_down)*t`` and ``Step(s)`` mixes the two adjacent
    integers with mean ``s``. Thus the mean height is always ``offset``.
    Zero offset conserves layer coverage; physical mass conservation also
    requires equivalent material inventories in the transported layers.
    This recentered family is not a constant-rate kinetic semigroup.

    :param float t: Elapsed time in seconds, nonnegative.
    :param float k_up: Upward single-layer event rate in inverse seconds.
    :param float k_down: Downward single-layer event rate in inverse seconds.
    :param float offset: Target mean displacement in structural layers.
    :param float tail_probability: Maximum probability outside the support.
    :param float moment_tolerance:
        Maximum omitted absolute first moment, in structural layers.
    """

    t: float = 0.0
    k_up: float = 0.0
    k_down: float = 0.0
    offset: float = 0.0
    tail_probability: float = DEFAULT_TAIL_PROBABILITY
    moment_tolerance: float = DEFAULT_MOMENT_TOLERANCE

    def __post_init__(self):
        for name in ("t", "k_up", "k_down", "offset",
                     "tail_probability", "moment_tolerance"):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite")
            object.__setattr__(self, name, value)
        if min(self.t, self.k_up, self.k_down) < 0:
            raise ValueError("t, k_up and k_down must be nonnegative")
        _validate_tail_probability(self.tail_probability)
        if self.moment_tolerance <= 0:
            raise ValueError("moment_tolerance must be positive")
        if not np.isfinite(self.mu_up + self.mu_down):
            raise ValueError("rate-time products and their sum must be finite")
        if max(abs(self.offset),
               abs(self.offset - self.mu_up + self.mu_down)) >= 2**52:
            raise ValueError("centering shift exceeds float64 layer resolution")

    @property
    def mu_up(self):
        """Return the dimensionless expected upward event count."""
        return self.k_up * self.t

    @property
    def mu_down(self):
        """Return the dimensionless expected downward event count."""
        return self.k_down * self.t

    @property
    def expected_height_change(self):
        """Return the target mean displacement in structural layers."""
        return self.offset

    @property
    def mean(self):
        """Return the actual mean height change in structural layers."""
        return self.offset

    @property
    def variance(self):
        """Return height variance in layers squared, including interpolation."""
        _, fraction = self._shift()
        return self.mu_up + self.mu_down + fraction * (1 - fraction)

    def _shift(self):
        shift = self.offset - (self.mu_up - self.mu_down)
        lower = int(np.floor(shift))
        return lower, shift - lower

    def _raw_cdf(self, values, up, down):
        values = np.floor(np.asarray(values, dtype=np.float64))
        if up == 0 and down == 0:
            result = (values >= 0).astype(np.float64)
        elif down == 0:
            result = poisson.cdf(values, up)
        elif up == 0:
            result = poisson.sf(-values - 1, down)
        else:
            result = skellam.cdf(values, up, down)
        if not np.all(np.isfinite(result)):
            raise ValueError("Skellam CDF could not be evaluated accurately")
        return result

    def _height_cdf(self, values):
        lower, fraction = self._shift()
        return ((1 - fraction) * self._raw_cdf(
            np.asarray(values) - lower, self.mu_up, self.mu_down)
            + fraction * self._raw_cdf(
                np.asarray(values) - lower - 1, self.mu_up, self.mu_down))

    def occupancy(self, offsets):
        """Return ``P(H > n)`` at structural-layer offsets ``n``."""
        lower, fraction = self._shift()
        values = np.floor(np.asarray(offsets, dtype=np.float64))
        # Reflect X to evaluate the small upper tail directly, rather than
        # subtracting a CDF close to one (SciPy's generic discrete sf).
        return ((1 - fraction) * self._raw_cdf(
            lower - values - 1, self.mu_down, self.mu_up)
            + fraction * self._raw_cdf(
                lower - values, self.mu_down, self.mu_up))

    def correction(self, offsets):
        """Return the stable signed material correction to a sharp Film."""
        offsets = np.asarray(offsets, dtype=np.float64)
        return np.where(offsets < 0, -self._height_cdf(offsets),
                        self.occupancy(offsets))

    def probability(self, changes):
        """Return masses at integer height changes; zero off the layer grid."""
        changes = np.asarray(changes, dtype=np.float64)
        lower, fraction = self._shift()

        def pmf(values):
            if self.mu_up == 0 and self.mu_down == 0:
                return (values == 0).astype(np.float64)
            if self.mu_down == 0:
                return poisson.pmf(values, self.mu_up)
            if self.mu_up == 0:
                return poisson.pmf(-values, self.mu_down)
            return skellam.pmf(values, self.mu_up, self.mu_down)

        result = ((1 - fraction) * pmf(changes - lower)
                  + fraction * pmf(changes - lower - 1))
        if not np.all(np.isfinite(result)):
            raise ValueError("Skellam PMF could not be evaluated accurately")
        return np.where(changes == np.floor(changes), result, 0.0)

    def surface_occupancy(self, offsets):
        """Return exact exposed masses ``P(H=n+1)`` without folding tails."""
        super().surface_occupancy(offsets)  # Validate the layer grid.
        return self.probability(np.asarray(offsets) + 1)

    def omitted_bounds(self, lower, upper):
        """Bound probability and absolute first moment outside height bounds.

        :param int lower: Inclusive retained height, at or below zero.
        :param int upper: Inclusive retained height, at or above zero.
        :returns: Probability bound and absolute first-moment bound in layers.

        Chernoff bounds use the centered cumulant generating function;
        Cauchy--Schwarz then bounds the omitted first moment. These conservative
        bounds avoid cancellation between large signed raw Poisson moments.
        """
        if lower > 0 or upper < 0 or lower > upper:
            raise ValueError("height bounds must contain zero")
        if self.mu_up + self.mu_down == 0:
            shift, fraction = self._shift()
            heights = np.array([shift, shift + 1])
            masses = np.array([1 - fraction, fraction])
            omitted = (heights < lower) | (heights > upper)
            return (float(np.sum(masses[omitted])),
                    float(np.sum(np.abs(heights[omitted]) * masses[omitted])))
        shift, fraction = self._shift()

        def log_bound(sign, threshold):
            if sign < 0 and self.mu_down == 0 and -threshold < shift:
                return -np.inf
            if (sign > 0 and self.mu_up == 0
                    and threshold > shift + int(fraction > 0)):
                return -np.inf
            distance = threshold - sign * self.offset
            if distance <= 0:
                return 0.0
            tilt = min(20.0, np.log1p(distance / (self.variance + 1)))
            q = sign * tilt
            # Each expm1 remainder is nonnegative, including near zero.
            cgf = (self.mu_up * (np.expm1(q) - q)
                   + self.mu_down * (np.expm1(-q) + q)
                   + q * self.offset
                   + np.log1p(fraction * np.expm1(q)) - fraction * q)
            return min(0.0, cgf - tilt * threshold)

        log_probability = min(0.0, np.logaddexp(
            log_bound(-1, 1 - lower), log_bound(1, upper + 1)))
        second_moment = self.variance + self.offset**2
        moment = np.exp(0.5 * (np.log(second_moment) + log_probability))
        return float(np.exp(log_probability)), float(moment)

    def support(self):
        """Return layer bounds satisfying probability and first-moment limits."""
        shift, fraction = self._shift()
        if self.mu_up + self.mu_down == 0:
            return min(-1, shift - 1), max(0, shift + int(fraction > 0))
        radius = max(2.0, 8 * np.sqrt(self.variance))
        for _ in range(64):
            lower = min(-1, int(np.floor(self.offset - radius)) - 1)
            upper = max(0, int(np.ceil(self.offset + radius)))
            if self.mu_down == 0:
                lower = max(lower, min(-1, shift - 1))
            if self.mu_up == 0:
                upper = min(upper, max(0, shift + int(fraction > 0)))
            if max(abs(lower), abs(upper)) >= 2**52:
                raise ValueError("support exceeds float64 layer resolution")
            probability, moment = self.omitted_bounds(lower + 1, upper + 1)
            # Reserve margin for renormalization of a selected domain ensemble.
            if (probability <= self.tail_probability / 4
                    and moment <= self.moment_tolerance / 4
                    and abs(self.offset) * probability
                    <= self.moment_tolerance / 4):
                return lower, upper
            radius *= 1.5
        raise ValueError("Could not bound the recentered Skellam support")

    def tight_support(self):
        """Shrink the legacy interval using the same bounds and safety margins.

        This opt-in numerical diagnostic leaves :meth:`support` and all legacy
        surface behavior unchanged. Returned bounds are structural-layer
        offsets; retained heights are one index greater at both endpoints.
        """
        lower, upper = self.support()

        def acceptable(lo, hi):
            probability, moment = self.omitted_bounds(lo+1, hi+1)
            return (probability <= self.tail_probability/4
                    and moment <= self.moment_tolerance/4
                    and abs(self.offset)*probability <= self.moment_tolerance/4)

        while lower < -1 and acceptable(lower+1, upper):
            lower += 1
        while upper > 0 and acceptable(lower, upper-1):
            upper -= 1
        return lower, upper


@dataclass(frozen=True)
class SkellamProfile:
    """Describe a two-sided Skellam interface profile.

    ``width`` and ``asymmetry`` retain the existing CTR interface convention:
    width is in unit cells and asymmetry is dimensionless.

    :param float width:
        Interface width in unit cells.
    :param float asymmetry:
        Dimensionless interface asymmetry between -1 and 1.
    :param float tail_probability:
        Maximum total omitted probability mass.
    """

    width: float
    asymmetry: float = 0.0
    tail_probability: float = DEFAULT_TAIL_PROBABILITY

    def __post_init__(self):
        if self.width < 0:
            raise ValueError("Skellam width must be greater than or equal to zero")
        if abs(self.asymmetry) > 1:
            raise ValueError("Skellam asymmetry must be between -1 and 1")
        _validate_tail_probability(self.tail_probability)

    def parameters(self, layer_count):
        """Return ``mu1`` and ``mu2`` for a unit cell with ``layer_count`` layers."""
        sigma = self.width * int(layer_count)
        if sigma == 0:
            return 0.0, 0.0
        asymmetry = self.asymmetry
        if abs(asymmetry) == 1.0:
            asymmetry -= np.sign(asymmetry) * 1e-6
        skew = asymmetry / sigma
        if abs(sigma * skew) > 1:
            raise ValueError("abs(sigma * skew) must be smaller than one")
        mu1 = 0.5 * sigma**2 * (1.0 + skew * sigma)
        return mu1, sigma**2 - mu1

    def support(self, layer_count):
        """Return inclusive structural-layer support offsets."""
        mu1, mu2 = self.parameters(layer_count)
        if mu1 == 0 and mu2 == 0:
            return 0, 0
        tail = self.tail_probability / 2.0
        lower = int(skellam.ppf(tail, mu1, mu2))
        upper = int(skellam.ppf(1.0 - tail, mu1, mu2))
        return lower, upper

    def occupancy(self, offsets, layer_count):
        """Return upper-material probabilities at structural-layer offsets."""
        mu1, mu2 = self.parameters(layer_count)
        offsets = np.asarray(offsets, dtype=np.int64)
        if mu1 == 0 and mu2 == 0:
            return (offsets >= 0).astype(np.float64)
        return skellam.cdf(offsets, mu1, mu2)


PROFILE_TYPES = {
    "poisson": PoissonProfile,
    "skellam": SkellamProfile,
    "recentered_skellam": RecenteredSkellamProfile,
}


@dataclass(frozen=True)
class SupportCertificate:
    """Numerical upper-tail bounds, separate from physical exhaustion.

    First and second moments are in parent index units and index squared.
    ``normalization_factor`` applies only to non-exhausted retained states.
    The physical exhausted probability is never renormalized.
    """

    upper: int
    omitted_probability: float
    omitted_first_moment: float
    omitted_second_moment: float
    retained_probability: float
    normalization_factor: float


def capped_skellam_states(profile, minimum, *, second_moment_tolerance=1e-8,
                          max_states=4096):
    """Return states of ``max(H, minimum)`` and an upper-tail certificate.

    :param RecenteredSkellamProfile profile: Unbounded parent distribution.
    :param int minimum: Physical exhausted state in structural indices.
    :param float second_moment_tolerance: Omitted index-squared budget.
    :param int max_states: Explicit resource limit; never a clipping policy.
    :returns: Read-only indices, probabilities, and a SupportCertificate.
    :raises ValueError: For invalid settings or insufficient resources.
    """
    if not isinstance(profile, RecenteredSkellamProfile):
        raise TypeError("A RecenteredSkellamProfile is required")
    if int(minimum) != minimum or minimum > 0:
        raise ValueError("minimum must be an integer at or below zero")
    if (not np.isfinite(second_moment_tolerance)
            or second_moment_tolerance <= 0):
        raise ValueError("second_moment_tolerance must be finite and positive")
    if int(max_states) != max_states or max_states < 1:
        raise ValueError("max_states must be a positive integer")
    minimum = int(minimum)
    shift, fraction = profile._shift()
    variance = profile.variance
    mean = profile.mean
    bernoulli = fraction * (1 - fraction)
    third = profile.mu_up - profile.mu_down + bernoulli * (1 - 2*fraction)
    fourth_cumulant = (profile.mu_up + profile.mu_down
                       + bernoulli * (1 - 6*bernoulli))
    fourth = max(0., fourth_cumulant + 4*mean*third + 3*variance**2
                 + 6*mean**2*variance + mean**4)

    def bounds(upper):
        if profile.mu_up == 0 and upper >= shift + int(fraction > 0):
            return 0., 0., 0.
        distance = upper + 1 - mean
        if distance <= 0:
            return 1., np.sqrt(variance + mean**2), np.sqrt(fourth)
        tilt = min(20., np.log1p(distance/(variance + 1)))
        cgf = (profile.mu_up * (np.expm1(tilt) - tilt)
               + profile.mu_down * (np.expm1(-tilt) + tilt)
               + tilt*mean + np.log1p(fraction*np.expm1(tilt))
               - fraction*tilt)
        logp = min(0., cgf - tilt*(upper + 1))
        # Cauchy--Schwarz bounds E[|H|;tail] and E[H²;tail].
        rootp = np.exp(logp/2)
        return (np.exp(logp), np.sqrt(variance + mean**2)*rootp,
                np.sqrt(fourth)*rootp)

    upper = max(0, int(np.ceil(mean + 8*np.sqrt(variance))))
    if profile.mu_up == 0:
        upper = max(minimum, shift + int(fraction > 0))
    while True:
        p_bound, m_bound, s_bound = bounds(upper)
        if (p_bound <= profile.tail_probability/4
                and m_bound <= profile.moment_tolerance/4
                and s_bound <= second_moment_tolerance/4):
            break
        upper += max(1, (upper - minimum)//4)
        if upper - minimum + 1 > max_states:
            raise ValueError("max_states cannot meet the upper-tail budgets")
    if upper - minimum + 1 > max_states:
        raise ValueError("max_states cannot represent the physical state range")
    indices = np.arange(minimum, upper + 1, dtype=np.int64)
    probabilities = profile.probability(indices)
    probabilities[0] = float(profile._height_cdf(minimum))
    retained = float(probabilities.sum())
    above = float(profile.occupancy(minimum))
    retained_above = float(probabilities[1:].sum())
    factor = above/retained_above if retained_above else 1.
    if above and not retained_above:
        raise ValueError("No non-exhausted state represents positive mass")
    probabilities[1:] *= factor
    if (not np.all(np.isfinite(probabilities))
            or np.any(probabilities < 0)
            or abs(probabilities.sum() - 1) > 1e-12):
        raise ValueError("Could not normalize bounded Skellam states accurately")
    indices.flags.writeable = probabilities.flags.writeable = False
    return indices, probabilities, SupportCertificate(
        upper, float(p_bound), float(m_bound), float(s_bound), retained,
        float(factor))

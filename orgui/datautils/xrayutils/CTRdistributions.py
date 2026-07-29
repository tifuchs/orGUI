"""Statistical occupancy profiles for CTR interfaces and surfaces."""

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from scipy.stats import poisson, skellam


DEFAULT_TAIL_PROBABILITY = 1e-10
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
}

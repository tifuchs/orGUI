"""Statistical occupancy profiles for CTR interfaces and surfaces."""

from dataclasses import dataclass

import numpy as np
from scipy.stats import poisson, skellam


DEFAULT_TAIL_PROBABILITY = 1e-10


def _validate_tail_probability(tail_probability):
    tail_probability = float(tail_probability)
    if not 0.0 < tail_probability < 1.0:
        raise ValueError("tail_probability must be between zero and one")
    return tail_probability


class PoissonProfile:
    """Describe signed Poisson growth or etching of a surface.

    The random surface-height change in structural layers is
    ``offset + sign(mean_change) * Poisson(abs(mean_change))``.

    :param float mean_change:
        Signed mean process width in structural layers. Positive values model
        growth and negative values model dissolution or etching.
    :param float offset:
        Deterministic surface-height offset in structural layers. Set this to
        ``-mean_change`` for mean-preserving roughness.
    :param float tail_probability:
        Maximum omitted probability mass in the unrepresented Poisson tail.
    :param float mean:
        Legacy alias for ``mean_change``.
    """

    def __init__(
        self,
        mean_change=None,
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
        self.offset = float(offset)
        self.tail_probability = _validate_tail_probability(
            tail_probability
        )

    @property
    def mean(self):
        """Return the legacy name for the signed mean height change."""
        return self.mean_change

    @property
    def rate(self):
        """Return the non-negative Poisson rate in structural layers."""
        return abs(self.mean_change)

    @property
    def expected_height_change(self):
        """Return the expected surface-height change in structural layers."""
        return self.offset + self.mean_change

    def support(self):
        """Return the inclusive structural-layer support offsets."""
        if self.rate == 0:
            process_low = process_high = self.offset
        else:
            quantile = int(
                poisson.ppf(
                    1.0 - self.tail_probability, self.rate
                )
            )
            if self.mean_change > 0:
                process_low = self.offset
                process_high = self.offset + quantile
            else:
                process_low = self.offset - quantile
                process_high = self.offset
        lower = int(np.floor(min(0.0, process_low))) - 1
        upper = int(np.ceil(max(0.0, process_high)))
        return lower, upper

    def occupancy(self, offsets):
        """Return material occupancy at structural-layer offsets."""
        offsets = np.asarray(offsets, dtype=np.float64)
        if self.rate == 0:
            return (offsets < self.offset).astype(np.float64)
        if self.mean_change > 0:
            return poisson.sf(offsets - self.offset, self.rate)
        threshold = self.offset - offsets
        return poisson.cdf(np.ceil(threshold) - 1, self.rate)

    def correction(self, offsets):
        """Return occupancy relative to a sharp surface at zero."""
        offsets = np.asarray(offsets, dtype=np.float64)
        sharp = (offsets < 0).astype(np.float64)
        return self.occupancy(offsets) - sharp


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

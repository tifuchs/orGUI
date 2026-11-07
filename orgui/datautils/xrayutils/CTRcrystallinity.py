"""Independent recovery of coherent order with depth below a local boundary."""

import copy
from dataclasses import dataclass

import numpy as np

from .CTRdistributions import PoissonProfile, SkellamProfile
from .CTRutil import LinearFitFunctions


class CrystallinityProfile(LinearFitFunctions):
    """Coherent-order recovery C(d) from a distribution of disorder depths.

    ``C(d) = minimum + (1-minimum)*P(D <= d-offset)`` multiplies the
    material's cell/layer crystallinity. Depth zero is the first occupied
    structural layer adjacent to the local boundary. Deep material recovers
    its source crystallinity. This object never changes material occupancy.

    :param profile: Independent :class:`PoissonProfile` or :class:`SkellamProfile`.
        Poisson parameters and ``offset`` use structural layers. Skellam width
        retains its unit-cell convention, converted using the host layer count.
    :param float offset: Additional recovery-depth shift in structural layers.
    :param float minimum: Minimum dimensionless amplitude fraction in [0, 1].
    :param str name: Name used by the fit-parameter machinery.
    """

    def __init__(self, profile, *, offset=0.0, minimum=0.0, name="order"):
        super().__init__()
        if isinstance(profile, PoissonProfile):
            self.kind = "poisson"
            keys = ["mean_change", "alpha", "offset", "depth_offset", "minimum"]
            values = [
                profile.mean_change,
                profile.alpha,
                profile.offset,
                offset,
                minimum,
            ]
        elif isinstance(profile, SkellamProfile):
            self.kind = "skellam"
            keys = ["width", "asymmetry", "depth_offset", "minimum"]
            values = [profile.width, profile.asymmetry, offset, minimum]
        else:
            raise TypeError(
                "CrystallinityProfile requires PoissonProfile or SkellamProfile"
            )
        self.name = name
        self.tail_probability = profile.tail_probability
        self.parameterLookup = {key: index for index, key in enumerate(keys)}
        self.parameterLookup_inv = dict(enumerate(keys))
        self.parameterOrder = " ".join(keys)
        self.basis = np.asarray(values, dtype=float)
        self.basis_0 = self.basis.copy()
        self._validate()

    def _validate(self):
        if not np.all(np.isfinite(self.basis)) or not 0 <= self.basis[-1] <= 1:
            raise ValueError(
                "Crystallinity profile values must be finite; minimum in [0, 1]."
            )
        self._distribution()

    def _distribution(self):
        if self.kind == "poisson":
            return PoissonProfile(
                *self.basis[:3], tail_probability=self.tail_probability
            )
        return SkellamProfile(*self.basis[:2], tail_probability=self.tail_probability)

    def evaluate(self, depth, layer_count=1):
        """Return dimensionless coherent order at local depths in structural layers."""
        self._validate()
        depth = np.asarray(depth, dtype=float)
        if not np.all(np.isfinite(depth)) or np.any(depth < 0):
            raise ValueError(
                "Local crystallinity depths must be finite and non-negative."
            )
        if not isinstance(layer_count, (int, np.integer)) or layer_count < 1:
            raise ValueError("layer_count must be a positive integer.")
        coordinate = depth - self.basis[-2]
        distribution = self._distribution()
        if self.kind == "poisson":
            cumulative = 1.0 - distribution.occupancy(coordinate)
        else:
            cumulative = distribution.occupancy(np.floor(coordinate), layer_count)
        return self.basis[-1] + (1.0 - self.basis[-1]) * cumulative

    def setFitParameters(self, values):
        """Apply independent order parameters and validate their physical bounds."""
        super().setFitParameters(values)
        self._validate()

    def updateFromParameters(self):
        """Restore fitted order values without changing any occupancy profile."""
        super().updateFromParameters()
        self._validate()

    def to_settings(self, include_fit=False):
        """Return JSON/HDF5-compatible live profile settings and optional fit state."""
        self._validate()
        result = dict(
            kind=self.kind,
            name=self.name,
            values=self.basis.tolist(),
            tail_probability=self.tail_probability,
        )
        if include_fit:
            for parameter, value in zip(
                self.parameter_list(), self.getInitialParameters()
            ):
                parameter.value = value
            result["fit"] = super().parametersToDict()
        return result

    @classmethod
    def from_settings(cls, settings):
        """Restore a profile written by :meth:`to_settings`."""
        values = np.asarray(settings["values"], dtype=float)
        kind = settings["kind"]
        if isinstance(kind, bytes):
            kind = kind.decode()
        tail = float(settings["tail_probability"])
        if kind == "poisson" and values.shape == (5,):
            distribution = PoissonProfile(*values[:3], tail_probability=tail)
        elif kind == "skellam" and values.shape == (4,):
            distribution = SkellamProfile(*values[:2], tail_probability=tail)
        else:
            raise ValueError("Invalid crystallinity profile family or parameter count.")
        name = settings["name"]
        if isinstance(name, bytes):
            name = name.decode()
        result = cls(distribution, offset=values[-2], minimum=values[-1], name=name)
        if "fit" in settings:
            result.parametersFromDict(settings["fit"])
            result.basis[:] = values
        return result


@dataclass
class _OrderBinding:
    owner: object
    role: str
    offset: int
    depths: dict

    def factor(self, layer):
        """Return the local dimensionless order multiplier for one source layer."""
        return self.owner._order_factor(self.role, self.offset, self.depths[layer])


def _depths(cell):
    labels = cell._crystallinity_labels()
    layers = sorted(
        cell.crystallinity_layers,
        key=lambda layer: np.mean(cell.basis[labels == layer, 3]),
        reverse=True,
    )
    return {layer: depth for depth, layer in enumerate(layers)}


def _bound_cell(cell, matrix, occupancy, owner, role, offset):
    """Make one occurrence with live order, retaining its physical material weight."""
    result = copy.copy(cell)
    result._crystallinity_source = cell
    result._layer_crystallinity = {}
    result.coherentDomainMatrix = [matrix]
    result.coherentDomainOccupancy = [occupancy]
    result._crystallinity_modifier = _OrderBinding(
        owner, role, int(offset), _depths(cell)
    )
    return result

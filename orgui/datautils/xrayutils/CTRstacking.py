"""Layer-cycle state and transition helpers for CTR crystal stacks."""

from dataclasses import dataclass


@dataclass(frozen=True)
class LayerCycle:
    """Ordered local structural-layer identifiers for one crystal structure."""

    layers: tuple

    def __init__(self, layers):
        layers = tuple(layers)
        if not layers:
            raise ValueError("A layer cycle must contain at least one layer")
        if len(set(layers)) != len(layers):
            raise ValueError("Layer identifiers in a cycle must be unique")
        object.__setattr__(self, "layers", layers)

    def successor(self, layer):
        """Return the cyclic successor of ``layer``."""
        try:
            index = self.layers.index(layer)
        except ValueError as exc:
            raise ValueError(
                f"Layer {layer!r} is not present in cycle {self.layers!r}"
            ) from exc
        return self.layers[(index + 1) % len(self.layers)]

    def is_rotation_of(self, other):
        """Return whether ``other`` has the same cyclic order."""
        if len(self.layers) != len(other.layers):
            return False
        doubled = self.layers + self.layers
        size = len(self.layers)
        return any(
            doubled[start:start + size] == other.layers
            for start in range(size)
        )


@dataclass(frozen=True)
class LayerState:
    """Top structural-layer state propagated during stack layout."""

    cycle: LayerCycle
    layer: float


@dataclass(frozen=True)
class LayerTransition:
    """Map lower terminating layers to upper starting layers."""

    mapping: dict

    def __init__(self, mapping):
        object.__setattr__(self, "mapping", dict(mapping))

    def upper_start(self, lower_layer):
        """Return the upper starting layer for ``lower_layer``."""
        try:
            return self.mapping[lower_layer]
        except KeyError as exc:
            raise ValueError(
                "No layer transition is defined for lower layer "
                f"{lower_layer!r}"
            ) from exc


def resolve_upper_start(lower_state, upper_cycle, transition=None):
    """Resolve the first upper layer from the lower stack state."""
    if lower_state is None or lower_state.layer == -1:
        return upper_cycle.layers[0]
    if transition is not None:
        start = transition.upper_start(lower_state.layer)
        if start not in upper_cycle.layers:
            raise ValueError(
                f"Transition selects layer {start!r} outside upper cycle "
                f"{upper_cycle.layers!r}"
            )
        return start
    if not lower_state.cycle.is_rotation_of(upper_cycle):
        raise ValueError(
            f"Layer cycles {lower_state.cycle.layers!r} and "
            f"{upper_cycle.layers!r} differ; provide layer_transition"
        )
    return upper_cycle.successor(lower_state.layer)

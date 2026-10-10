"""Runtime numerical checks for pyFAI's detector-array calculations."""

from functools import lru_cache
import logging

import numpy as np
from pyFAI import units


logger = logging.getLogger(__name__)


def _strided_angles_are_safe():
    """Check radians on the interleaved coordinate layout used by pyFAI.

    Direct compiled NumExpr calls in some pyFAI/NumExpr combinations ignore
    the coordinate views' strides. A contiguous-input test misses this fault.
    The reference uses NumPy, independently of pyFAI's unit equations.
    """
    positions = np.arange(1.0, 25.0).reshape(2, 4, 3)
    x, y, z = positions[..., 2], positions[..., 1], positions[..., 0]
    expected = (
        np.arctan2(np.hypot(x, y), z),
        np.arctan2(y, x),
    )
    for unit, reference in zip((units.TTH_RAD, units.CHI_RAD), expected):
        result = np.asarray(unit.equation(x, y, z, 1e-10))
        if (result.shape != reference.shape
                or not np.all(np.isfinite(result))
                or np.max(np.abs(result - reference)) > 1e-12):
            return False
    return True


@lru_cache(maxsize=1)
def _ensure_pyfai_array_safety():
    """Probe once per process, warning when contiguous inputs are needed."""
    safe = _strided_angles_are_safe()
    if not safe:
        logger.warning(
            "Environment error: orGUI detected numerical artifacts in "
            "pyFAI's strided angle-array check. Some numexpr versions cause "
            "incorrect geometry and polarization corrections. Try reinstalling "
            "numexpr in the Python environment used by orGUI, then restart "
            "orGUI. Reinstalling may fix the problem. orGUI is attempting a "
            "best-effort workaround using contiguous coordinate inputs, but "
            "correct results are not guaranteed.",
            extra={"show_dialog": True, "title": "Numerical environment error"},
        )
    return safe


def _contiguous_equation(equation, x, y, z=None, *args, **kwargs):
    """Keep coordinates/units intact, copying only noncontiguous arrays."""
    coordinates = tuple(
        np.ascontiguousarray(value)
        if isinstance(value, np.ndarray) and not value.flags.c_contiguous
        else value
        for value in (x, y, z)
    )
    return equation(*coordinates, *args, **kwargs)

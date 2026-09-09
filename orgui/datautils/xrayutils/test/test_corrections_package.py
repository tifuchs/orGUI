"""The corrections package: shared definitions and its released aliases.

:mod:`orgui.datautils.xrayutils.corrections` exists so that the rocking
integration, the stationary integration and the reciprocal-space
reconstruction correct their data with the same code rather than three
copies. These tests pin the two properties that makes true: that the shared
functions really are what the callers use, and that the module paths released
before the package existed still work.
"""

import importlib

import numpy as np
import pytest

from orgui.datautils.xrayutils import corrections
from orgui.datautils.xrayutils.corrections import (
    detector,
    normalization,
    roi,
)


class _Detector:
    """Minimal stand-in for a calibrated detector."""

    def __init__(self, shape=(4, 3)):
        self._shape = shape

    def solidAngleArray(self):
        return np.linspace(0.9, 1.0, int(np.prod(self._shape))).reshape(self._shape)

    def polarizationArray(self):
        return np.linspace(0.8, 1.0, int(np.prod(self._shape))).reshape(self._shape)


def test_the_package_exposes_every_correction_module():
    """A maintainer looking for a correction factor finds them in one place."""
    assert set(corrections.__all__) == {
        "beamprofile",
        "detector",
        "geometry",
        "normalization",
        "roi",
    }
    for name in corrections.__all__:
        assert getattr(corrections, name) is not None


@pytest.mark.parametrize(
    "legacy, moved",
    [
        (
            "orgui.datautils.xrayutils.geometrycorrections",
            "orgui.datautils.xrayutils.corrections.geometry",
        ),
        (
            "orgui.datautils.xrayutils.beamprofile",
            "orgui.datautils.xrayutils.corrections.beamprofile",
        ),
    ],
)
def test_released_module_paths_still_import(legacy, moved):
    """Both names were released, so both must keep working.

    Not just importable: the objects have to be the same ones, or a caller
    holding the old path would be testing a different implementation from the
    one the integration uses.
    """
    old = importlib.import_module(legacy)
    new = importlib.import_module(moved)

    assert old.__all__ == new.__all__
    for name in old.__all__:
        assert getattr(old, name) is getattr(new, name)


def test_pixel_factors_is_the_reciprocal_of_the_enabled_arrays():
    """The one definition of the per-pixel divisors.

    ``None`` for nothing enabled is part of the contract: the reciprocal-space
    reconstruction skips its multiplication entirely on that, rather than
    walking a full detector of ones.
    """
    det = _Detector()

    assert detector.pixel_factors(det) is None
    np.testing.assert_allclose(
        detector.pixel_factors(det, solid_angle=True),
        1.0 / det.solidAngleArray(),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        detector.pixel_factors(det, polarization=True),
        1.0 / det.polarizationArray(),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        detector.pixel_factors(det, solid_angle=True, polarization=True),
        1.0 / (det.solidAngleArray() * det.polarizationArray()),
        rtol=1e-12,
    )


def test_normalization_takes_values_not_a_scan():
    """The physics layer never learns a beamline's counter names.

    Pulling counters off a scan object is the application's job; this module
    only knows what a usable value is.
    """
    divisor, applied = normalization.normalization_divisor(
        3, exposure_time=0.5, monitors={"mon": [100.0, 200.0, 400.0]}
    )

    np.testing.assert_allclose(divisor, 0.5 * np.array([100.0, 200.0, 400.0]))
    assert applied == ["exposure", "monitor:mon"]

    none, applied = normalization.normalization_divisor(3)
    np.testing.assert_allclose(none, np.ones(3))
    assert applied == []

    with pytest.raises(ValueError, match="finite and positive"):
        normalization.normalization_divisor(3, exposure_time=0.0)
    with pytest.raises(ValueError, match="finite and nonzero"):
        normalization.normalization_divisor(3, monitors={"mon": [1.0, 0.0, 1.0]})
    with pytest.raises(ValueError, match="has 2 values for 3 images"):
        normalization.normalization_divisor(3, monitors={"mon": [1.0, 2.0]})


def test_correction_factors_bundle_multiplies_the_named_subset():
    """The container the applied factors are carried and stored in."""
    factors = roi.CorrectionFactors(
        {"C_norm": np.full(3, 2.0), "C_Lorentz": np.full(3, 4.0)},
        ["normalization", "lorentz"],
    )

    np.testing.assert_allclose(factors.divisor("C_norm"), np.full(3, 2.0))
    np.testing.assert_allclose(
        factors.divisor("C_norm", "C_Lorentz"), np.full(3, 8.0)
    )
    # A correction that was not enabled is skipped, not an error.
    np.testing.assert_allclose(factors.divisor("C_illum_area"), 1.0)
    assert factors.applied == ("normalization", "lorentz")

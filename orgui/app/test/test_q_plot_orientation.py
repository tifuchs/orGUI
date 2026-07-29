"""The Q-plot action must use the sample orientations validated against QAlpha.

The numerical equivalence of the Q-plot conversion with orGUI's own QAlpha
calculation is established in
``orgui/datautils/xrayutils/test/test_q_conversion.py`` for a fixed mapping of
the azimuthal reference onto pyFAI's EXIF sample orientation. This test keeps
the GUI implementation tied to that validated mapping.
"""

import numpy as np
import pytest

from orgui.datautils.xrayutils.test.test_q_conversion import AZIMUTH_ORIENTATION

orGUI = pytest.importorskip("orgui.app.orGUI", reason="Qt bindings are required")

# Which plot axis has to be flipped so that the converted image keeps the same
# handedness as the pixel image, see orGUI._qConversionSampleOrientation.
EXPECTED_FLIP = {0.0: None, 90.0: "y", 180.0: "x", 270.0: None}


class _StubDetectorCal:
    def __init__(self, azimuth):
        self._azimuth = azimuth

    def getAzimuthalReference(self):
        return self._azimuth


class _StubUBCalculator:
    def __init__(self, azimuth):
        self.detectorCal = _StubDetectorCal(azimuth)


class _StubOrGUI:
    """Minimal stand-in: the method only reads the azimuthal reference."""

    def __init__(self, azimuth):
        self.ubcalc = _StubUBCalculator(azimuth)


@pytest.mark.parametrize("azimuth_deg, orientation", AZIMUTH_ORIENTATION)
def test_orientation_matches_validated_mapping(azimuth_deg, orientation):
    stub = _StubOrGUI(np.deg2rad(azimuth_deg))

    result, flipaxis = orGUI.orGUI._qConversionSampleOrientation(stub)

    assert result == orientation
    assert flipaxis == EXPECTED_FLIP[azimuth_deg]


def test_fiber_integrator_flag_matches_pyfai_version():
    """The Q-plot action is only offered when FiberIntegrator can be imported."""
    if orGUI.HAS_FIBER_INTEGRATOR:
        assert hasattr(orGUI, "FiberIntegrator")
    else:
        assert not hasattr(orGUI, "FiberIntegrator")

"""The plot status line reports HKL and the selected reciprocal-space frame.

The readout shows ``HKL`` and ``Q[<frame>]`` as bracketed triplets, and the
``Q`` field is relabelled when a different frame is selected.
"""

import numpy as np
import pytest

from orgui.app import qconversion
from orgui.datautils.xrayutils import HKLVlieg

orGUI = pytest.importorskip("orgui.app.orGUI", reason="Qt bindings are required")


ENERGY = 78.0  # keV

# refracted [alpha, delta, gamma, omega, chi, phi] in rad
POSITION = [
    np.deg2rad(0.6),
    np.deg2rad(12.0),
    np.deg2rad(4.0),
    np.deg2rad(23.0),
    np.deg2rad(3.0),
    np.deg2rad(-5.0),
]

EXPECTED_FIELD_NAMES = {
    "Q_alpha": "Q[alpha]",
    "Q_lab": "Q[lab]",
    "Q_omega": "Q[omega]",
    "Q_chi": "Q[chi]",
    "Q_phi": "Q[phi]",
    "Q_cryst": "Q[cryst]",
}


def _ub_calculator(with_orientation=True):
    lattice = HKLVlieg.Lattice(np.array([4.0, 4.0, 4.0]), np.array([90.0, 90.0, 90.0]))
    ub = HKLVlieg.UBCalculator(lattice, ENERGY)
    if with_orientation:
        ub.defaultU()
    return ub


class _StubUBCalculatorWidget:
    def __init__(self, ub):
        self.ubCal = ub
        self.angles = HKLVlieg.VliegAngles(ub)


class _StubOrGUI:
    """The readout helper only needs the frame and the UB calculator."""

    def __init__(self, frame, ub):
        self._frame = frame
        self.ubcalc = _StubUBCalculatorWidget(ub)

    def _selectedQFrame(self):
        return self._frame


def test_triplet_is_formatted_as_a_bracketed_vector():
    """Five decimals are the minimum useful precision for hkl and for Q."""
    formatted = orGUI.formatPositionTriplet([1.23456, -0.5, 2.0])

    assert formatted == "[1.23456 -0.50000 2.00000]"


def test_triplet_template_is_wide_enough_for_the_format():
    """The size hint must cover the widest value the format can produce."""
    widest = orGUI.formatPositionTriplet([-99.999994, -99.999994, -99.999994])

    assert len(widest) <= len(orGUI.POSITION_TRIPLET_TEMPLATE)


@pytest.mark.parametrize("values", [[np.nan, 1.0, 2.0], [np.inf, 1.0, 2.0]])
def test_triplet_without_a_scan_shows_a_placeholder(values):
    assert orGUI.formatPositionTriplet(values) == "------"


def test_every_frame_has_an_expected_field_name():
    assert set(EXPECTED_FIELD_NAMES) == set(qconversion.FRAMES)


@pytest.mark.parametrize("frame, expected", sorted(EXPECTED_FIELD_NAMES.items()))
def test_field_name_follows_the_frame(frame, expected):
    assert orGUI.qFieldName(frame) == expected


def test_alpha_frame_reports_qalpha():
    stub = _StubOrGUI("Q_alpha", _ub_calculator())

    q = orGUI.orGUI._qInSelectedFrame(stub, POSITION)

    expected = stub.ubcalc.angles.QAlpha(*POSITION[:3])
    np.testing.assert_allclose(q, expected, atol=1e-12)


def test_crystal_frame_reports_b_times_hkl():
    """Q_cryst is the crystal Cartesian frame, so it must equal B times hkl."""
    ub = _ub_calculator()
    stub = _StubOrGUI("Q_cryst", ub)

    q = orGUI.orGUI._qInSelectedFrame(stub, POSITION)

    hkl = np.array(stub.ubcalc.angles.anglesToHkl(*POSITION))
    np.testing.assert_allclose(q, ub.lattice.B_mat @ hkl, atol=1e-9)


def test_crystal_frame_without_an_orientation_matrix_is_unavailable():
    stub = _StubOrGUI("Q_cryst", _ub_calculator(with_orientation=False))

    q = orGUI.orGUI._qInSelectedFrame(stub, POSITION)

    assert np.all(np.isnan(q))
    assert orGUI.formatPositionTriplet(q) == "------"


class _StubTitle:
    """Stand-in for the bold field title that silx creates once."""

    def __init__(self, text):
        self._text = text

    def text(self):
        return self._text

    def setText(self, text):
        self._text = text


class _StubPositionInfo:
    def __init__(self, titles):
        self._titles = titles
        self.updates = 0

    def findChildren(self, cls):
        del cls
        return self._titles

    def updateInfo(self):
        self.updates += 1


class _StubPlot:
    def __init__(self):
        self._qFieldName = "Q[alpha]"
        self.titles = [_StubTitle("<b>HKL:</b>"), _StubTitle("<b>Q[alpha]:</b>")]
        self.positionInfo = _StubPositionInfo(self.titles)

    def getPositionInfoWidget(self):
        return self.positionInfo


@pytest.mark.parametrize("frame, name", sorted(EXPECTED_FIELD_NAMES.items()))
def test_selecting_a_frame_relabels_only_the_q_field(frame, name):
    plot = _StubPlot()

    orGUI.Plot2DHKL.setQFrame(plot, frame)

    assert plot.titles[1].text() == f"<b>{name}:</b>"
    assert plot.titles[0].text() == "<b>HKL:</b>"  # untouched
    assert plot._qFieldName == name


def test_reselecting_the_same_frame_is_a_noop():
    plot = _StubPlot()

    orGUI.Plot2DHKL.setQFrame(plot, "Q_alpha")

    assert plot.positionInfo.updates == 0
    assert plot.titles[1].text() == "<b>Q[alpha]:</b>"


def test_relabelling_refreshes_the_displayed_values():
    plot = _StubPlot()

    orGUI.Plot2DHKL.setQFrame(plot, "Q_phi")

    assert plot.positionInfo.updates == 1

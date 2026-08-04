"""The Q-plot must be wired to the validated conversion in :mod:`qconversion`.

The numerical equivalence of the Q-plot with orGUI's own QAlpha calculation is
established in ``orgui/app/test/test_q_conversion.py``. These tests only cover
the GUI wiring: which frames are offered, which one is preselected, and how the
selection is read back.
"""

import os

import numpy as np
import pyFAI
import pytest

from orgui.app import qconversion
from orgui.datautils.xrayutils import DetectorCalibration, HKLVlieg

orGUI = pytest.importorskip("orgui.app.orGUI", reason="Qt bindings are required")


class _StubComboBox:
    """Minimal stand-in for the frame selector combo box."""

    def __init__(self, data):
        self._data = data

    def currentData(self):
        return self._data


class _StubOrGUI:
    """The frame accessor only reads the combo box."""

    def __init__(self, data):
        self.qFrameSelector = _StubComboBox(data)


@pytest.mark.parametrize("frame", qconversion.FRAMES)
def test_selected_frame_is_read_from_the_combo_box(frame):
    assert orGUI.orGUI._selectedQFrame(_StubOrGUI(frame)) == frame


@pytest.mark.parametrize("data", [None, "Q_nonsense", ""])
def test_unknown_selection_falls_back_to_the_default_frame(data):
    selected = orGUI.orGUI._selectedQFrame(_StubOrGUI(data))
    assert selected == qconversion.DEFAULT_FRAME


def test_default_frame_is_the_alpha_frame():
    """QAlpha is what the rest of orGUI uses, so it is the natural default."""
    assert qconversion.DEFAULT_FRAME == "Q_alpha"
    assert qconversion.DEFAULT_FRAME in qconversion.FRAMES


def test_every_frame_has_a_label_for_the_selector():
    for frame in qconversion.FRAMES:
        assert qconversion.FRAME_LABELS[frame]


def test_fiber_integrator_flag_is_taken_from_qconversion():
    """The action is only usable when FiberIntegrator can be imported."""
    assert orGUI.HAS_FIBER_INTEGRATOR == qconversion.HAS_FIBER_INTEGRATOR


@pytest.mark.parametrize("suffix", [".svg", ".png"])
def test_the_toolbar_icon_ships_with_the_package(suffix):
    """The action loads its icon by name, which silently yields a blank icon.

    ``QIcon`` never reports a missing file, so a packaging mistake would only
    show up as an empty toolbar button.
    """
    from orgui import resources

    assert os.path.isfile(resources.getPath("q-plot" + suffix))


class _StubAction:
    def __init__(self, checked):
        self._checked = checked

    def isChecked(self):
        return self._checked


class _RefreshStub:
    """Records the toggles the refresh performs."""

    def __init__(self, checked, reentrant=False):
        self.plotAgainstQAct = _StubAction(checked)
        self._qPlotRefreshing = False
        self._reentrant = reentrant
        self.calls = []

    def _convertImagetoQ(self, value):
        self.calls.append(value)
        if self._reentrant:
            # emulate a signal arriving while the plot is being rebuilt
            orGUI.orGUI._refreshQPlot(self)


def test_refresh_does_nothing_while_the_q_plot_is_off():
    stub = _RefreshStub(checked=False)

    orGUI.orGUI._refreshQPlot(stub)

    assert stub.calls == []


def test_refresh_rebuilds_in_place_and_stays_in_q_mode():
    """Changing the image must not drop out of the reciprocal-space view."""
    stub = _RefreshStub(checked=True)

    orGUI.orGUI._refreshQPlot(stub)

    # a single conversion replaces the previous one; the action stays checked
    assert stub.calls == [True]
    assert stub.plotAgainstQAct.isChecked()


def test_refresh_accepts_the_replot_signal_argument():
    """sigReplotRequest carries a bool, which the slot must tolerate."""
    stub = _RefreshStub(checked=True)

    orGUI.orGUI._refreshQPlot(stub, True)

    assert stub.calls == [True]


def test_refresh_is_not_reentrant():
    stub = _RefreshStub(checked=True, reentrant=True)

    orGUI.orGUI._refreshQPlot(stub)

    assert stub.calls == [True]
    assert stub._qPlotRefreshing is False


def test_suspended_q_plot_batches_into_one_rebuild():
    """Untoggling max/sum before replotting must not rebuild several times."""
    stub = _RefreshStub(checked=True)

    with orGUI.orGUI._suspendedQPlot(stub):
        orGUI.orGUI._refreshQPlot(stub)
        orGUI.orGUI._refreshQPlot(stub)
        assert stub.calls == []
    orGUI.orGUI._refreshQPlot(stub)

    assert stub.calls == [True]
    assert stub._qPlotRefreshing is False


# --- the position readout while the Q-plot is on -------------------------
#
# The plot axes are momentum transfer there, not pixels, so the readout has to
# invert the conversion. These tests drive the real methods against a fully
# specified geometry.

ENERGY = 78.0  # keV
ALPHA = np.deg2rad(0.6)
OMEGA = np.deg2rad(23.0)


def _geometry(azimuth_deg=45.0):
    det = DetectorCalibration.Detector2D_SXRD()
    det.detector = pyFAI.detector_factory("Pilatus2m")
    det.wavelength = (12.39842 / ENERGY) * 1e-10
    det.dist = 0.729
    det.poni1 = 0.21
    det.poni2 = 0.14
    det.setAzimuthalReference(np.deg2rad(azimuth_deg))
    return det


class _StubUBCalculatorWidget:
    def __init__(self, det):
        lattice = HKLVlieg.Lattice(
            np.array([4.0, 4.0, 4.0]), np.array([90.0, 90.0, 90.0])
        )
        self.ubCal = HKLVlieg.UBCalculator(lattice, ENERGY)
        self.ubCal.defaultU()
        self.angles = HKLVlieg.VliegAngles(self.ubCal)
        self.detectorCal = det
        self.chi = 0.0
        self.phi = 0.0
        self.n = 1.0


class _ReadoutStub:
    """Everything the Q-plot readout inversion reads."""

    def __init__(self, frame="Q_alpha"):
        self._frame = frame
        self.ubcalc = _StubUBCalculatorWidget(_geometry())

    def _selectedQFrame(self):
        return self._frame

    def _qPlotAngles(self):
        return ALPHA, OMEGA

    _qPlotOrientationMatrix = orGUI.orGUI._qPlotOrientationMatrix
    _qPlotToPosition = orGUI.orGUI._qPlotToPosition
    _isOnDetector = orGUI.orGUI._isOnDetector


def _plot_coordinates(stub, row, col):
    """Q-plot coordinates of one detector pixel, in inverse Angstrom."""
    det = stub.ubcalc.detectorCal
    gamma, delta = det.surfaceAnglesPoint(
        np.array([float(row)]), np.array([float(col)]), ALPHA
    )
    q_alpha = stub.ubcalc.angles.QAlpha(ALPHA, delta, gamma)[0]
    rotation = qconversion.frameMatrix(
        stub._selectedQFrame(),
        alpha=ALPHA,
        omega=OMEGA,
        chi=stub.ubcalc.chi,
        phi=stub.ubcalc.phi,
        U=stub.ubcalc.ubCal.getU(),
    )
    Q = rotation @ q_alpha
    return np.copysign(np.hypot(Q[0], Q[1]), Q[0]), Q[2], delta[0], gamma[0]


@pytest.mark.parametrize("frame", qconversion.FRAMES)
@pytest.mark.parametrize("pixel", [(900, 700), (600, 1000), (150, 1400)])
def test_readout_inverts_the_conversion_back_to_the_same_pixel(frame, pixel):
    """Hovering a Q-plot point must report the angles of the pixel it came from."""
    stub = _ReadoutStub(frame)
    q_ip, q_oop, delta, gamma = _plot_coordinates(stub, *pixel)

    resolved = orGUI.orGUI._qPlotToPosition(stub, q_ip, q_oop)

    assert resolved is not None
    alpha, omega, recovered_delta, recovered_gamma, _ = resolved
    assert (alpha, omega) == (ALPHA, OMEGA)
    np.testing.assert_allclose(
        [recovered_delta, recovered_gamma], [delta, gamma], atol=1e-9
    )


@pytest.mark.parametrize("frame", qconversion.FRAMES)
@pytest.mark.parametrize("pixel", [(900, 700), (600, 1000), (150, 1400)])
def test_reported_q_reproduces_the_cursor_coordinates(frame, pixel):
    """The reported Q must be the point the cursor is on, not a corrected one.

    The Q-plot axes carry no refraction correction, so applying one to the
    reported momentum transfer would contradict the displayed position.
    """
    stub = _ReadoutStub(frame)
    q_ip, q_oop, _, _ = _plot_coordinates(stub, *pixel)

    _, _, _, _, Q = orGUI.orGUI._qPlotToPosition(stub, q_ip, q_oop)

    np.testing.assert_allclose(np.copysign(np.hypot(Q[0], Q[1]), Q[0]), q_ip, atol=1e-9)
    np.testing.assert_allclose(Q[2], q_oop, atol=1e-9)


def test_readout_reports_nothing_outside_the_detector():
    """The rebinned grid is rectangular and reaches past the detector."""
    stub = _ReadoutStub()

    assert orGUI.orGUI._qPlotToPosition(stub, 0.0, -5.0) is None


def test_readout_reports_nothing_beyond_the_ewald_sphere():
    stub = _ReadoutStub()
    k = stub.ubcalc.ubCal.getK()

    assert orGUI.orGUI._qPlotToPosition(stub, 10.0 * k, 0.0) is None


# --- which way the Q-plot ordinate runs ----------------------------------
#
# The detector image is drawn with the row index increasing downwards. The
# Q-plot should keep that visual orientation, which means its ordinate has to
# follow the geometry rather than always pointing upwards.


def _out_of_plane_at_top_and_bottom(det, alpha):
    """``q_perp`` at the top-centre and bottom-centre pixel, in A^-1."""
    lattice = HKLVlieg.Lattice(np.array([4.0, 4.0, 4.0]), np.array([90.0, 90.0, 90.0]))
    ub = HKLVlieg.UBCalculator(lattice, ENERGY)
    ub.defaultU()
    shape = det.detector.shape
    rows = np.array([0.0, float(shape[0] - 1)])
    cols = np.full(2, 0.5 * (shape[1] - 1))
    gamma, delta = det.surfaceAnglesPoint(rows, cols, alpha)
    return HKLVlieg.VliegAngles(ub).QAlpha(alpha, delta, gamma)[..., 2]


@pytest.mark.parametrize("azimuth_deg", [0.0, 45.0, 90.0, 135.0, 180.0, 270.0])
@pytest.mark.parametrize("alpha_deg", [0.6, -0.6])
def test_ordinate_direction_matches_the_detector_image(azimuth_deg, alpha_deg):
    """The flag must say where q_perp really grows on the displayed image."""
    det = _geometry(azimuth_deg)
    alpha = np.deg2rad(alpha_deg)
    top, bottom = _out_of_plane_at_top_and_bottom(det, alpha)

    flipped = qconversion.outOfPlaneIncreasesWithRow(det, alpha)

    assert flipped == bool(bottom > top)


def test_upward_scattering_keeps_the_ordinate_pointing_up():
    """The ordinary geometry must not be inverted."""
    det = _geometry(90.0)

    assert not qconversion.outOfPlaneIncreasesWithRow(det, np.deg2rad(0.6))


def test_inverted_geometry_inverts_the_ordinate():
    """ID31's downward scattering setup, reached by rotating the azimuth."""
    det = _geometry(270.0)

    assert qconversion.outOfPlaneIncreasesWithRow(det, np.deg2rad(0.6))


@pytest.mark.parametrize("frame", qconversion.FRAMES)
def test_ordinate_direction_is_defined_for_every_frame(frame):
    lattice = HKLVlieg.Lattice(np.array([4.0, 4.0, 4.0]), np.array([90.0, 90.0, 90.0]))
    ub = HKLVlieg.UBCalculator(lattice, ENERGY)
    ub.defaultU()
    det = _geometry(270.0)

    flipped = qconversion.outOfPlaneIncreasesWithRow(
        det,
        np.deg2rad(0.6),
        frame=frame,
        omega=OMEGA,
        chi=0.05,
        phi=-0.1,
        U=ub.getU(),
    )

    assert isinstance(flipped, bool)


# --- restoring the axis orientation when the Q-plot is switched off -------


class _AxisPlotStub:
    """Records the axis orientation the Q-plot toggle applies."""

    def __init__(self, xInverted, yInverted):
        self.xInverted = xInverted
        self.yInverted = yInverted
        self.zoomResets = 0

    def isXAxisInverted(self):
        return self.xInverted

    def isYAxisInverted(self):
        return self.yInverted

    def setXAxisInverted(self, flag):
        self.xInverted = flag

    def setYAxisInverted(self, flag):
        self.yInverted = flag

    def resetZoom(self):
        self.zoomResets += 1

    # not exercised by the restore path
    def setQPlotActive(self, active):
        pass

    def getAllImages(self):
        return []

    def setActiveImage(self, legend):
        pass


class _RestoreStub:
    """The part of orGUI the axis restore touches."""

    def __init__(self, xInverted, yInverted):
        self.centralPlot = _AxisPlotStub(xInverted, yInverted)
        self._qPlotPreviousXInverted = None
        self._qPlotPreviousYInverted = None
        self.currentAddImageLabel = None
        self.currentImageLabel = "scan_image"
        self.allimgmax = None
        self.allimgsum = None

    def _enter(self, flipOrdinate=False):
        """The axis part of switching the Q-plot on."""
        if self._qPlotPreviousYInverted is None:
            self._qPlotPreviousYInverted = self.centralPlot.isYAxisInverted()
        if self._qPlotPreviousXInverted is None:
            self._qPlotPreviousXInverted = self.centralPlot.isXAxisInverted()
        self.centralPlot.setYAxisInverted(flipOrdinate)
        self.centralPlot.setXAxisInverted(False)


@pytest.mark.parametrize("xInverted", [False, True])
@pytest.mark.parametrize("yInverted", [False, True])
def test_leaving_the_q_plot_restores_both_axes(xInverted, yInverted):
    """Both origins must come back, whatever they were set to."""
    stub = _RestoreStub(xInverted, yInverted)
    stub._enter()

    orGUI.orGUI._convertImagetoQ(stub, False)

    assert stub.centralPlot.xInverted == xInverted
    assert stub.centralPlot.yInverted == yInverted
    assert stub._qPlotPreviousXInverted is None
    assert stub._qPlotPreviousYInverted is None


def test_restore_does_not_depend_on_the_max_sum_image():
    """A conversion that failed after flipping the axes must still roll back.

    The restore used to sit inside the block that rebuilds the maximum or sum
    image, so it was skipped whenever no additional image was displayed.
    """
    stub = _RestoreStub(xInverted=True, yInverted=True)
    stub._enter(flipOrdinate=True)
    assert stub.currentAddImageLabel is None  # nothing was ever added

    orGUI.orGUI._convertImagetoQ(stub, False)

    assert stub.centralPlot.xInverted is True
    assert stub.centralPlot.yInverted is True
    assert stub.centralPlot.zoomResets == 1


def test_leaving_without_ever_entering_touches_nothing():
    stub = _RestoreStub(xInverted=True, yInverted=False)

    orGUI.orGUI._convertImagetoQ(stub, False)

    assert stub.centralPlot.xInverted is True
    assert stub.centralPlot.yInverted is False
    assert stub.centralPlot.zoomResets == 0


def test_rebuilding_in_place_keeps_the_recorded_orientation():
    """Stepping through images re-enters the Q-plot; the record must survive."""
    stub = _RestoreStub(xInverted=False, yInverted=True)
    stub._enter()
    stub._enter()  # a rebuild, not a fresh entry

    orGUI.orGUI._convertImagetoQ(stub, False)

    assert stub.centralPlot.yInverted is True

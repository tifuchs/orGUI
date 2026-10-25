"""Mu-scan overlays must follow the same fixed-azimuth Ewald intersection."""

from types import SimpleNamespace

import numpy as np
import pytest
from silx.gui import qt

from orgui.app.orGUI import orGUI
from orgui.app.QUBCalculator import QMachineParameters, QUBCalculator
from orgui.datautils.xrayutils import DetectorCalibration, HKLVlieg


@pytest.mark.parametrize("chi,phi", [(0, 0.06), (0.03, 0.06), (-0.03, -0.06)])
@pytest.mark.parametrize("omega", [0.0, 23.0])
@pytest.mark.parametrize("n", [1.0, 1.0 - 1.1415e-6])
def test_mu_scan_00l_overlay_matches_stationary(chi, phi, omega, n):
    """The same tilted rod must land on the same pixel in both ROI modes."""
    mu, omega, chi, phi = np.deg2rad([0.5, omega, chi, phi])
    ub = HKLVlieg.UBCalculator(HKLVlieg.Lattice([4, 4, 6.8], [90, 90, 90]), 77.0)
    ub.setU(np.eye(3))
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(1000, 1500, 1500, pixelX=100, pixelY=100)
    detector.detector.shape = (3000, 3000)
    detector.wavelength = ub.getLambda() * 1e-10
    gui = SimpleNamespace(
        fscan=[0, 1, 2],
        imageno=0,
        getMuOm=lambda imageno=None: (
            (mu, omega) if imageno is not None
            else (mu + np.array([0, 0.001, 0.002]), np.full(3, omega))
        ),
        getArmAngles=lambda imageno=None: (None, None),
        _qRangePerFrame=lambda shape, imageno: (0, np.inf),
        scanSelector=SimpleNamespace(
            scanstab=SimpleNamespace(currentIndex=lambda: 0),
            roioptions=SimpleNamespace(get_offsets=lambda: (0, 0)),
        ),
    )
    calculator = SimpleNamespace(
        mainGui=gui,
        mu=mu,
        chi=chi,
        phi=phi,
        n=n,
        angles=HKLVlieg.VliegAngles(ub),
        detectorCal=detector,
    )
    calculator.calcReflection = lambda hkl: QUBCalculator.calcReflection(
        calculator, hkl, axisname="mu"
    )
    gui.ubcalc = calculator
    gui._pixelsCrystalAnglesPerFrame = lambda *args: orGUI._pixelsCrystalAnglesPerFrame(
        gui, *args
    )
    origin = np.array([0.0, 0.0, 0.0])
    direction = np.array([0.0, 0.0, 1.0])
    stationary = orGUI.getROIloc(gui, 0, origin, direction)
    # The other intersection is the direct beam at L=0.
    point = max((branch[0] for branch in stationary), key=lambda row: row[2])
    assert point[-1]
    rocking = orGUI.get_rocking_coordinates(
        gui, origin, direction, maxValue=point[2], step_width=point[2]
    )
    for branch in (1, 2):
        np.testing.assert_allclose(rocking[f"xy_{branch}"][1], point[6:8], atol=1e-7)
        np.testing.assert_allclose(
            rocking[f"angles_{branch}"][1, [0, 3, 4, 5]],
            [mu, omega, chi, phi],
            atol=1e-10,
        )
    single = calculator.calcReflection(point[:3])
    np.testing.assert_allclose(single["xy_1"], point[6:8], atol=1e-7)


@pytest.mark.parametrize("chi,phi", [(-0.03, -0.06), (-120, 180), (270, -270)])
def test_machine_angles_round_trip_signed_degrees(chi, phi):
    """Loading and editing signed angles must preserve their radian values."""
    app = qt.QApplication.instance() or qt.QApplication([])
    widget = QMachineParameters()
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(1000, 1500, 1500, pixelX=100, pixelY=100)
    detector.detector.max_shape = detector.detector.shape = (3000, 3000)
    detector.wavelength = 0.16e-10
    try:
        widget.setValues(
            {
                "diffractometer": dict(
                    zip(("mu", "chi", "phi"), np.deg2rad([0.5, chi, phi]))
                ),
                "source": {"E": 77.0},
                "SXRD_geometry": detector,
            }
        )
        params = widget.getParameters()["diffractometer"]
        np.testing.assert_allclose(
            [params["chi"], params["phi"]], np.deg2rad([chi, phi])
        )
        emitted = []
        widget.sigMachineParamsChanged.connect(emitted.append)
        widget.phibox.setValue(phi + 0.01)
        assert emitted[-1]["diffractometer"]["phi"] == pytest.approx(
            np.deg2rad(phi + 0.01)
        )
    finally:
        widget.close()
        widget.deleteLater()
        app.processEvents()

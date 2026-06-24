import logging
import copy
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app.orGUI import orGUI
from orgui.app.QReflectionSelector import QReflectionSelector
from orgui.app.QUBCalculator import (
    QUBCalculator,
    QUBFitDialog,
    _DeprecatedFitOption,
)
from orgui.datautils.xrayutils import DetectorCalibration, HKLVlieg


class FakeScan:
    axisname = "th"
    axis = np.array([0.0, 1.0, 2.0])
    omega = np.array([0.0, 1.0, 2.0])

    def __len__(self):
        return len(self.axis)


class FakeDetectorCalibration:
    def surfaceAnglesPoint(self, y, x, mu):
        return np.array([0.5]), np.array([0.25])


class FakeBraggParent:
    fscan = SimpleNamespace(axisname="th")

    def axisToImageNo(self, axisval):
        if 0.0 <= axisval <= 2.0:
            return int(axisval)
        raise ValueError("axis outside scan")


def make_gui(reflections=None):
    gui = orGUI.__new__(orGUI)
    gui.fscan = FakeScan()
    gui.ubcalc = SimpleNamespace(
        detectorCal=FakeDetectorCalibration(),
        mu=0.1,
        chi=0.2,
        phi=0.3,
    )
    gui.reflectionSel = SimpleNamespace(reflections=reflections or [])
    return gui


def test_image_number_converters_reject_stale_indices():
    gui = make_gui()

    assert gui.imageNoToAxis(1) == pytest.approx(1.0)
    assert gui.imageNoToOmega(1) == pytest.approx(np.deg2rad(1.0))

    with pytest.raises(IndexError):
        gui.imageNoToAxis(3)
    with pytest.raises(IndexError):
        gui.imageNoToOmega(3)
    with pytest.raises(IndexError):
        gui.getMuOm(3)


def test_get_reflections_skips_stale_image_numbers(caplog):
    reflections = [
        SimpleNamespace(
            hkl=np.array([1.0, 0.0, 0.0]),
            xy=np.array([10.0, 20.0]),
            imageno=1,
        ),
        SimpleNamespace(
            hkl=np.array([0.0, 1.0, 0.0]),
            xy=np.array([30.0, 40.0]),
            imageno=5,
        ),
    ]
    gui = make_gui(reflections)

    with caplog.at_level(logging.WARNING, logger="orgui.app.orGUI"):
        hkls, angles = gui.getReflections()

    assert hkls.tolist() == [[1.0, 0.0, 0.0]]
    assert angles.shape == (1, 6)
    assert angles[0, 3] == pytest.approx(np.deg2rad(1.0))
    assert "Skipping stale reflection" in caplog.text


def test_change_image_skips_stale_image_numbers(caplog):
    gui = make_gui()
    gui.scanSelector = SimpleNamespace(
        slider=SimpleNamespace(
            setValue=lambda imageno: pytest.fail("stale image should not be selected"),
            value=lambda: pytest.fail("stale image should not be plotted"),
        )
    )
    gui.plotImage = lambda imageno: pytest.fail("stale image should not be plotted")

    with caplog.at_level(logging.WARNING, logger="orgui.app.orGUI"):
        gui._onChangeImage(4)

    assert "Skipping request to display stale image number" in caplog.text


def test_bragg_reflection_list_skips_stale_image_numbers(caplog):
    selector = QReflectionSelector.__new__(QReflectionSelector)
    selector.showBraggReflections = False
    selector.reflBragg = []
    selector.plot = SimpleNamespace(removeMarker=lambda identifier: None)
    selector.orparent = FakeBraggParent()

    hkls = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    yx = np.array([[20.0, 10.0], [40.0, 30.0]])
    angles = np.zeros((2, 6))
    angles[0, 3] = np.deg2rad(-1.0)
    angles[1, 3] = np.deg2rad(-5.0)

    with caplog.at_level(logging.WARNING, logger="orgui.app.QReflectionSelector"):
        selector.setBraggReflections(hkls, yx, angles)

    assert len(selector.reflBragg) == 1
    assert selector.reflBragg[0].hkl.tolist() == [1.0, 0.0, 0.0]
    assert selector.reflBragg[0].imageno == 1
    assert "Skipping Bragg reflection" in caplog.text


def test_reflection_mismatch_tolerates_table_update_in_progress():
    class FakeModel:
        def __init__(self, shape):
            self.data = np.empty(shape)
            self.colors = None

        def getData(self):
            return self.data

        def setArrayColors(self, bgcolors=None, fgcolors=None):
            self.colors = bgcolors

    class FakeViewport:
        def update(self):
            pass

    def fake_editor(shape):
        return SimpleNamespace(
            model=FakeModel(shape),
            view=SimpleNamespace(viewport=lambda: FakeViewport()),
        )

    selector = QReflectionSelector.__new__(QReflectionSelector)
    selector.reflections = [object(), object()]
    selector.refleditor = fake_editor((2, 6))
    selector.refleditor_angles = fake_editor((0, 9))
    label = SimpleNamespace(text=None)
    selector.mismatchLabel = SimpleNamespace(
        setText=lambda text: setattr(label, "text", text)
    )
    mismatch = {
        "angle_mismatch": np.array([0.01, 0.02]),
        "norm_mismatch": np.array([0.1, 0.2]),
        "relative_norm_mismatch": np.array([0.03, 0.04]),
    }

    selector.setReflectionMismatch(mismatch)

    assert selector.refleditor.model.colors.shape == (2, 6, 3)
    assert selector.refleditor_angles.model.colors is None
    assert label.text.startswith("Mismatch:")


def test_set_reflections_synchronizes_both_tables():
    selector = QReflectionSelector.__new__(QReflectionSelector)
    selector.reflections = []
    selector.activeReflection = None
    selector.nextNo = 0
    selector._showReferenceReflections = False
    selector.plot = SimpleNamespace(removeMarker=lambda identifier: None)
    update_calls = []
    selector.updateEditor = lambda: update_calls.append(True)
    reflections = [
        SimpleNamespace(identifier="old-0"),
        SimpleNamespace(identifier="old-1"),
    ]

    selector.setReflections(reflections)

    assert update_calls == [True]
    assert [refl.identifier for refl in selector.reflections] == [
        "ref_0",
        "ref_1",
    ]


def test_add_reflection_refreshes_mismatch_through_editor_update():
    selector = QReflectionSelector.__new__(QReflectionSelector)
    selector.reflections = []
    selector.nextNo = 0
    selector._showReferenceReflections = False
    refresh_calls = []
    selector.updateEditor = lambda: refresh_calls.append(True)
    selector.setReflectionActive = lambda identifier: None

    selector.addReflection(
        {"x": 10.0, "y": 20.0},
        imageno=3,
        hkl=np.array([1.0, 0.0, 0.0]),
    )

    assert refresh_calls == [True]
    assert len(selector.reflections) == 1


def test_u_change_refreshes_mismatch():
    calls = []
    calculator = SimpleNamespace(
        ubCal=SimpleNamespace(setU=lambda value: calls.append(("setU", value))),
        updateReflectionMismatch=lambda: calls.append(("refresh", None)),
        sigReplotRequest=SimpleNamespace(
            emit=lambda value: calls.append(("replot", value))
        ),
    )
    u_matrix = np.eye(3)

    QUBCalculator._onUchanged(calculator, u_matrix)

    assert calls[0][0] == "setU"
    assert np.array_equal(calls[0][1], u_matrix)
    assert calls[1:] == [("refresh", None), ("replot", True)]


def test_lattice_change_refreshes_mismatch():
    calls = []
    lattice = object()
    calculator = SimpleNamespace(
        ubCal=SimpleNamespace(
            setLattice=lambda value: calls.append(("setLattice", value))
        ),
        updateReflectionMismatch=lambda: calls.append(("refresh", None)),
        sigReplotRequest=SimpleNamespace(
            emit=lambda value: calls.append(("replot", value))
        ),
    )

    QUBCalculator._onCrystalParamsChanged(calculator, lattice, 0.999)

    assert calculator.crystal is lattice
    assert calculator.n == 0.999
    assert calls == [
        ("setLattice", lattice),
        ("refresh", None),
        ("replot", True),
    ]


def test_setting_reflection_handler_refreshes_mismatch():
    calls = []
    calculator = SimpleNamespace(
        updateReflectionMismatch=lambda: calls.append("refresh")
    )

    def reflection_handler():
        return np.empty((0, 3)), np.empty((0, 6))

    QUBCalculator.setReflectionHandler(calculator, reflection_handler)

    assert calculator.reflections is reflection_handler
    assert calls == ["refresh"]


def test_experiment_fit_recovers_detector_distance():
    lattice = HKLVlieg.Lattice([3.9, 4.1, 5.2], [90.0, 90.0, 90.0])
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(
        1000.0, 500.0, 600.0, pixelX=100.0, pixelY=100.0
    )
    detector.wavelength = 12.39842 / 70.0 * 1e-10
    xy = np.array(
        [[300.0, 400.0], [700.0, 450.0],
         [350.0, 800.0], [800.0, 900.0]]
    )
    expected_u = HKLVlieg.Rotation.from_euler(
        "xyz", [0.1, -0.2, 0.3]
    ).as_matrix()
    wavevector = 2 * np.pi / (detector.wavelength * 1e10)
    angles = []
    q_phi = []
    for index, (x, y) in enumerate(xy):
        gamma, delta = detector.surfaceAnglesPoint(
            np.array([y]), np.array([x]), 0.01
        )
        position = np.array(
            [0.01, delta[0], gamma[0],
             -0.3 + 0.15 * index, 0.02, -0.04]
        )
        angles.append(position)
        q_phi.append(HKLVlieg.calculate_q_phi(position, wavevector))
    angles = np.asarray(angles)
    q_phi = np.asarray(q_phi)
    hkls = np.linalg.solve(
        expected_u @ lattice.B_mat, q_phi.T
    ).T

    starting_detector = copy.deepcopy(detector)
    starting_detector.dist = 0.85
    starting_detector.reset()
    reset_calls = []
    original_reset = starting_detector.reset

    def counted_reset():
        reset_calls.append(True)
        original_reset()

    starting_detector.reset = counted_reset
    ub_calculator = HKLVlieg.UBCalculator(lattice, 70.0)
    ub_calculator.setU(np.eye(3))
    calculator = SimpleNamespace(
        mainGui=SimpleNamespace(
            getReflectionFitData=lambda: (hkls, angles, xy),
            reflectionSel=SimpleNamespace(updateEditor=lambda: None),
        ),
        detectorCal=starting_detector,
        crystal=lattice,
        ubCal=ub_calculator,
        n=1.0,
        mu=0.01,
        chi=0.02,
        phi=-0.04,
        crystalparams=SimpleNamespace(setValues=lambda *args: None),
        machineParams=SimpleNamespace(setValues=lambda *args: None),
        uedit=SimpleNamespace(setU=lambda *args: None),
        updateReflectionMismatch=lambda: None,
        sigReplotRequest=SimpleNamespace(emit=lambda *args: None),
    )

    rms = QUBCalculator.fitExperiment(calculator, {"dist"})

    assert calculator.detectorCal.dist == pytest.approx(1.0)
    assert calculator.ubCal.getU() == pytest.approx(expected_u)
    assert rms < 1e-10
    assert reset_calls == [True]


def test_deprecated_lattice_fit_option_redirects_to_dialog():
    modes = []
    dialog = SimpleNamespace(
        setLegacyLatticeMode=lambda mode: modes.append(mode),
        legacyLatticeMode=lambda: modes[-1] if modes else "none",
    )
    calculator = SimpleNamespace(fitDialog=dialog)
    option = _DeprecatedFitOption(calculator, "scale", "latscale")

    with pytest.warns(DeprecationWarning):
        option.setChecked(True)
    with pytest.warns(DeprecationWarning):
        assert option.isChecked()

    assert modes == ["scale"]


def test_legacy_scale_mode_routes_calc_u_to_new_fit_api():
    fit_calls = []
    calculator = SimpleNamespace(
        reflections=lambda: (
            np.ones((4, 3)),
            np.ones((4, 6)),
        ),
        ubCal=SimpleNamespace(
            calculateUFromReflections=lambda hkls, angles: None,
            getU=lambda: np.eye(3),
        ),
        uedit=SimpleNamespace(setU=lambda value: None),
        updateReflectionMismatch=lambda: None,
        sigReplotRequest=SimpleNamespace(emit=lambda value: None),
        fitDialog=SimpleNamespace(legacyLatticeMode=lambda: "scale"),
        fitExperiment=lambda parameters: fit_calls.append(parameters),
    )

    QUBCalculator._onCalcU(calculator)

    assert fit_calls == [{"lattice_scale"}]


def test_lattice_only_fit_does_not_copy_trial_detector(monkeypatch):
    true_lattice = HKLVlieg.Lattice(
        [3.9, 4.1, 5.2], [90.0, 90.0, 90.0]
    )
    starting_lattice = HKLVlieg.Lattice(
        true_lattice.a * 1.2, [90.0, 90.0, 90.0]
    )
    expected_u = HKLVlieg.Rotation.from_euler(
        "xyz", [0.1, -0.2, 0.3]
    ).as_matrix()
    angles = np.deg2rad(
        [
            [0.2, 4.0, 2.0, 10.0, 1.0, -3.0],
            [0.2, 7.0, -1.0, 15.0, 1.0, -3.0],
            [0.2, 3.0, 5.0, 25.0, 1.0, -3.0],
            [0.2, -2.0, 6.0, 35.0, 1.0, -3.0],
        ]
    )
    energy = 70.0
    wavevector = 2 * np.pi / (12.39842 / energy)
    q_phi = np.array([
        HKLVlieg.calculate_q_phi(position, wavevector)
        for position in angles
    ])
    hkls = np.linalg.solve(
        expected_u @ true_lattice.B_mat, q_phi.T
    ).T

    detector = SimpleNamespace(
        dist=1.0,
        poni1=0.0,
        poni2=0.0,
        rot1=0.0,
        rot2=0.0,
        rot3=0.0,
        energy=energy,
        reset=lambda: None,
    )
    ub_calculator = HKLVlieg.UBCalculator(starting_lattice, energy)
    ub_calculator.setU(np.eye(3))
    calculator = SimpleNamespace(
        mainGui=SimpleNamespace(
            getReflectionFitData=lambda: (
                hkls, angles, np.zeros((len(hkls), 2))
            ),
            reflectionSel=SimpleNamespace(updateEditor=lambda: None),
        ),
        detectorCal=detector,
        crystal=starting_lattice,
        ubCal=ub_calculator,
        n=1.0,
        mu=0.01,
        chi=0.0,
        phi=0.0,
        crystalparams=SimpleNamespace(setValues=lambda *args: None),
        machineParams=SimpleNamespace(setValues=lambda *args: None),
        uedit=SimpleNamespace(setU=lambda *args: None),
        sigReplotRequest=SimpleNamespace(emit=lambda *args: None),
    )
    monkeypatch.setattr(
        "orgui.app.QUBCalculator.copy.deepcopy",
        lambda value: pytest.fail("lattice-only fit copied detector geometry"),
    )

    QUBCalculator.fitExperiment(calculator, {"lattice_scale"})

    assert starting_lattice.a == pytest.approx(true_lattice.a)


def test_fit_dialog_updates_mean_discrepancies():
    class FakeLabel:
        def __init__(self):
            self.text = None

        def setText(self, text):
            self.text = text

    dialog = SimpleNamespace(
        calculator=SimpleNamespace(
            getReflectionMismatch=lambda: {
                "angle_mismatch": np.deg2rad([1.0, 3.0]),
                "norm_mismatch": np.array([0.1, 0.3]),
            }
        ),
        meanQLabel=FakeLabel(),
        meanAngleLabel=FakeLabel(),
    )

    QUBFitDialog.updateDiscrepancy(dialog)

    assert dialog.meanQLabel.text == (
        "Mean Q discrepancy: 0.2 Angstrom^-1"
    )
    assert dialog.meanAngleLabel.text == (
        "Mean angle discrepancy: 2 deg"
    )


def test_fit_dialog_reset_restores_opening_state():
    calls = []
    state = {"marker": "opening"}
    dialog = SimpleNamespace(
        _savedState=state,
        _restoreState=lambda value: calls.append(("restore", value)),
        setCurrentValues=lambda: calls.append(("values", None)),
        updateDiscrepancy=lambda: calls.append(("discrepancy", None)),
        rmsLabel=SimpleNamespace(
            setText=lambda text: calls.append(("rms", text))
        ),
    )

    QUBFitDialog.resetParameters(dialog)

    assert calls == [
        ("restore", state),
        ("values", None),
        ("discrepancy", None),
        ("rms", "RMS relative Q discrepancy: unavailable"),
    ]

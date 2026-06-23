import logging
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app.orGUI import orGUI
from orgui.app.QReflectionSelector import QReflectionSelector
from orgui.app.QUBCalculator import QUBCalculator


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
    mismatch = {
        "angle_mismatch": np.array([0.01, 0.02]),
        "relative_norm_mismatch": np.array([0.03, 0.04]),
    }

    selector.setReflectionMismatch(mismatch)

    assert selector.refleditor.model.colors.shape == (2, 6, 3)
    assert selector.refleditor_angles.model.colors is None


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

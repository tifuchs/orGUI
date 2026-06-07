import logging
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app.orGUI import orGUI
from orgui.app.QReflectionSelector import QReflectionSelector


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

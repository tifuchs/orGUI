"""The integration settings must survive a save and a reload, in full.

Before the typed NeXus layout, reloading a configuration restored three of
the seven correction switches and none of the region-of-interest settings, so
a saved dataset could not be reproduced from what was stored with it. These
tests drive the real :class:`~orgui.app.QScanSelector.QScanSelector` widgets
rather than a stand-in, because the gap was in the wiring between the widgets
and :class:`~orgui.app.config_data.ConfigData`, not in either alone.
"""

import numpy as np
import pytest
from silx.gui import qt

from orgui.app.config_data import ROIState, roi_from_nxdict, roi_to_nxdict
from orgui.app.QScanSelector import QScanSelector


@pytest.fixture(scope="module")
def qapp():
    application = qt.QApplication.instance()
    if application is None:
        application = qt.QApplication([])
    return application


class _StubMainWindow(qt.QMainWindow):
    """The smallest parent a selector will build against.

    ``QScanSelector.__init__`` reads two attributes off its parent -- the
    plot the alpha slider drives and the current image legend -- and also
    passes it to child dialogs as a Qt parent, so it has to be a real
    widget rather than a namespace.
    """

    def __init__(self):
        super().__init__()
        from silx.gui.plot import Plot2D

        self.centralPlot = Plot2D(parent=self)
        self.currentAddImageLabel = "image"


def _make_selector():
    parent = _StubMainWindow()
    selector = QScanSelector(parent)
    # Keep the parent alive for as long as the selector needs it.
    selector._test_parent = parent
    return selector


@pytest.fixture
def selector(qapp):
    widget = _make_selector()
    yield widget
    widget.deleteLater()


#: Every switch off its default, so a dropped one is visible.
SWITCHES = {
    "mask": True,
    "solid_angle": True,
    "polarization": True,
    "lorentz": True,
    "footprint": True,
    "normalization": True,
}

REGION = {
    "hsize": 21.0,
    "vsize": 7.0,
    "left": 19.0,
    "right": 18.0,
    "top": 4.0,
    "bottom": 3.0,
    "auto_hsize": True,
    "auto_vsize": True,
}

ADVANCED = {
    "detector_inclination": True,
    "project_sample_size": True,
    "offset_x": 1.5,
    "offset_y": -2.5,
    "sample_size_x": 5e-4,
    "sample_size_y": 7e-3,
    "sample_size_z": 5e-3,
    "factor": 1.25,
    "fitted_background": True,
    "fitted_background_order": 2,
}

ROCKING = {"delta_s": 0.0019, "max_s": 5.0}


def test_every_setting_round_trips_through_the_widgets(selector):
    """Set everything, read it back, get the same thing."""
    selector.set_integration_options({
        **SWITCHES,
        "region": REGION,
        "advanced": ADVANCED,
        "rocking_scan": ROCKING,
    })

    options = selector.get_integration_options()

    for name, value in SWITCHES.items():
        assert options[name] == value, name
    assert options["region"] == pytest.approx(REGION)
    assert options["advanced"] == pytest.approx(ADVANCED)
    assert options["rocking_scan"] == pytest.approx(ROCKING)


def test_restoring_delta_s_does_not_re_clip_it(selector):
    """The stored value is already the resolution-clipped one.

    Assigning it normally triggers ``onRoSChanged``, which solves for the
    step that keeps rod points one pixel apart and writes the result back.
    Re-running that against whatever scan is loaded would not restore what
    was saved, so the setter blocks it.
    """
    selector.set_integration_options({"rocking_scan": {"delta_s": 0.0019}})

    assert selector.roscanDeltaS.value() == pytest.approx(0.0019)


def test_the_settings_survive_a_nexus_group(selector, tmp_path):
    """Widgets -> ROIState -> NeXus -> ROIState -> widgets is the identity."""
    from silx.io.dictdump import dicttonx, nxtodict

    selector.set_integration_options({
        "region": REGION, "advanced": ADVANCED, "rocking_scan": ROCKING,
    })
    options = selector.get_integration_options()
    state = ROIState(
        region=dict(options["region"]),
        advanced=dict(options["advanced"]),
        rocking_scan=dict(options["rocking_scan"]),
    )

    path = tmp_path / "roi.h5"
    dicttonx({"roi": roi_to_nxdict(state)}, path)
    restored = roi_from_nxdict(nxtodict(path)["roi"])

    fresh = _make_selector()
    try:
        fresh.set_integration_options({
            "region": restored.region,
            "advanced": restored.advanced,
            "rocking_scan": restored.rocking_scan,
        })
        after = fresh.get_integration_options()
    finally:
        fresh.deleteLater()

    assert after["region"] == pytest.approx(options["region"])
    assert after["advanced"] == pytest.approx(options["advanced"])
    assert after["rocking_scan"] == pytest.approx(options["rocking_scan"])


def test_a_legacy_option_dictionary_is_still_accepted(selector):
    """``startup_setup.py``-style scripts keep working for a cycle."""
    with pytest.deprecated_call():
        selector.set_integration_options({
            "solidAngle": True,
            "advanced": {
                "DetectorInclination": True,
                "ProjectSampleSize": True,
                "xoffset": 0.0,
                "yoffset": 0.0,
                "sizeX": 5e-4,
                "sizeY": 7e-3,
                "sizeZ": 5e-3,
                "factor": 1.0,
            },
        })

    options = selector.get_integration_options()
    assert options["solid_angle"] is True
    assert options["advanced"]["detector_inclination"] is True
    assert options["advanced"]["sample_size_y"] == pytest.approx(7e-3)
    # and the legacy spelling still reads off the result, with a warning
    with pytest.deprecated_call():
        assert options["solidAngle"] is True


def test_an_unknown_option_is_ignored(selector):
    """What lets an older orGUI open a newer configuration."""
    selector.set_integration_options({"mask": True, "from_the_future": 1})

    assert selector.get_integration_options()["mask"] is True


def test_sample_sizes_are_stored_in_meter(selector):
    """The widgets show micrometer; the dictionary and the file use meter."""
    selector.set_integration_options({"advanced": {**ADVANCED,
                                                   "sample_size_y": 7e-3}})

    assert selector.get_integration_options()["advanced"][
        "sample_size_y"
    ] == pytest.approx(7e-3)
    assert np.isclose(selector.roioptions._sizeYsample.value(), 7000.0)

"""Live footprint diagnostics, angle sources and foreground-only updates."""

import numpy as np
import pytest
from silx.gui import qt

from orgui.app.peak1Dintegr import IntegrationCorrectionsDialog
from orgui.app import sample_interception_dialog as interception_dialog


class _Selector(qt.QObject):
    sigScanChanged = qt.Signal(object)
    sigImageNoChanged = qt.Signal(int)


class _Scan:
    axisname = "th"
    phi = np.array([0., 45.])

    def __len__(self):
        return len(self.phi)


@pytest.fixture
def preview():
    """Create an initially hidden editor with a replaceable scan source."""
    app = qt.QApplication.instance() or qt.QApplication([])
    host = qt.QWidget()
    host.setAttribute(qt.Qt.WA_DontShowOnScreen, True)
    selector = _Selector(host)
    main = host
    main.fscan = _Scan()
    main.imageno = 0
    main.scanSelector = selector
    main.getMuOm = lambda: (np.deg2rad([0.36, 0.5]), 0.)
    host._mainWindow = lambda: main
    dialog = IntegrationCorrectionsDialog(host)
    dialog.setAttribute(qt.Qt.WA_DontShowOnScreen, True)
    dialog.setSettings({
        "analytical": True, "shape": "Gaussian", "shape_values": [160],
        "sample_interception": {
            "version": 1, "enabled": True,
            "shape": {"kind": "rectangle", "dimensions_m": [0.01, 0.006]},
            "azimuth_source": "phi", "azimuth_unit": "deg",
            "offset_m": [0.002, 0.001],
            "horizontal": {
                "analytical": True, "shape": "Gaussian", "shape_values": [2800],
            },
        },
    })
    editor = dialog.sampleEditor
    editor._live_timer.setInterval(0)
    yield app, dialog, editor, main
    dialog.hide()
    dialog.deleteLater()
    host.deleteLater()
    app.sendPostedEvents(None, qt.QEvent.DeferredDelete)


def _events(app):
    for _ in range(4):
        app.processEvents()


def _show(app, dialog):
    dialog.show()
    app.setActiveWindow(dialog)
    _events(app)


def _image(editor):
    return editor.footprint_plot.getImage("Beam density (1/m²)").getData(copy=True)


def test_live_edits_update_right_panel_without_quadrature(preview, monkeypatch):
    """Both profiles and geometry update without calculating a scan curve."""
    app, dialog, editor, main = preview
    monkeypatch.setattr(interception_dialog, "shape_frame_factors",
                        lambda *a, **k: pytest.fail("live preview ran quadrature"))
    assert editor.footprint_plot is None
    _show(app, dialog)
    assert editor._livePreviewActive()
    assert editor.live_panel.isVisibleTo(dialog)
    assert editor.splitter.widget(1) is editor.live_panel
    assert editor.overlap_plot is None
    before = _image(editor)
    dialog.profileOffset.setValue(40)
    _events(app)
    assert not np.array_equal(before, _image(editor))
    before = _image(editor)
    editor.horizontal.profileOffset.setValue(500)
    _events(app)
    assert not np.array_equal(before, _image(editor))
    before = editor.footprint_plot.getCurve("Sample").getXData(copy=True)
    dialog.L.setValue(14)
    _events(app)
    after = editor.footprint_plot.getCurve("Sample").getXData()
    assert not np.array_equal(before, after)


def test_scan_changes_and_independent_manual_angles(preview):
    """Frame/scan signals refresh readouts; a manual angle stays independent."""
    app, dialog, editor, main = preview
    _show(app, dialog)
    assert editor.live_readouts["incidence"].text() == "0.36 deg"
    main.imageno = 1
    main.scanSelector.sigImageNoChanged.emit(1)
    _events(app)
    assert editor.frame_index.value() == 1
    assert editor.live_readouts["incidence"].text() == "0.5 deg"
    assert editor.live_readouts["azimuth"].text() == "45 deg"
    saved = dialog.settings()
    editor.live_angle_modes["incidence"].setCurrentIndex(1)
    editor.preview_alpha.setValue(0.8)
    _events(app)
    assert "Incidence 0.8 deg" in editor.live_status.text()
    assert dialog.settings() == saved
    main.fscan = _Scan()
    main.fscan.phi = np.array([-12.])
    main.getMuOm = lambda: (np.deg2rad(0.4), 0.)
    main.imageno = 0
    main.scanSelector.sigScanChanged.emit(main.fscan)
    _events(app)
    assert editor.frame_index.maximum() == 0
    assert editor.live_readouts["azimuth"].text() == "-12 deg"
    assert "Incidence 0.8 deg" in editor.live_status.text()
    editor.live_angle_modes["azimuth"].setCurrentIndex(1)
    editor.preview_azimuth.setValue(27)
    main.fscan = None
    main.scanSelector.sigScanChanged.emit(None)
    _events(app)
    assert "azimuth readback 27 deg" in editor.live_status.text()
    assert editor.footprint_plot.getImage("Beam density (1/m²)") is not None


@pytest.mark.parametrize("change", ["image", "scan"])
def test_main_window_navigation_refreshes_visible_preview(preview, change):
    """Selecting images in the main window must not require clicking the dialog."""
    app, dialog, editor, main = preview
    main.show()
    _show(app, dialog)
    before = _image(editor)
    app.setActiveWindow(main)
    assert main.isActiveWindow()
    assert not dialog.isActiveWindow()
    main.imageno = 1
    if change == "scan":
        main.fscan = _Scan()
        main.scanSelector.sigScanChanged.emit(main.fscan)
    else:
        main.scanSelector.sigImageNoChanged.emit(1)
    _events(app)
    assert editor.frame_index.value() == 1
    assert editor.live_readouts["incidence"].text() == "0.5 deg"
    assert editor.live_readouts["azimuth"].text() == "45 deg"
    assert not np.array_equal(before, _image(editor))
    app.setActiveWindow(None)
    main.imageno = 0
    main.scanSelector.sigImageNoChanged.emit(0)
    _events(app)
    assert editor.frame_index.value() == 1
    app.setActiveWindow(main)
    _events(app)
    assert editor.frame_index.value() == 0
    dialog.hide()
    main.imageno = 1
    main.scanSelector.sigImageNoChanged.emit(1)
    _events(app)
    assert editor.frame_index.value() == 0


def test_live_preview_pauses_hidden_and_unfocused_and_catches_up(preview, monkeypatch):
    """No image work is done offscreen; activation renders the newest inputs."""
    app, dialog, editor, main = preview
    calls = []
    draw = editor._drawFootprint

    def _record(*args, **kwargs):
        calls.append(args)
        return draw(*args, **kwargs)

    monkeypatch.setattr(editor, "_drawFootprint", _record)
    dialog.L.setValue(12)
    _events(app)
    assert calls == []
    _show(app, dialog)
    assert len(calls) == 1
    app.setActiveWindow(None)
    dialog.L.setValue(13)
    main.imageno = 1
    main.scanSelector.sigImageNoChanged.emit(1)
    _events(app)
    assert len(calls) == 1
    assert not editor._live_timer.isActive()
    app.setActiveWindow(dialog)
    _events(app)
    assert len(calls) == 2
    assert editor.frame_index.value() == 1
    dialog.hide()
    dialog.L.setValue(15)
    _events(app)
    assert len(calls) == 2
    _show(app, dialog)
    assert len(calls) == 3
    assert calls[-1][0]["shape"]["dimensions_m"][0] == pytest.approx(0.015)


def test_invalid_geometry_clears_stale_preview_and_recovers(preview):
    """Invalid edits report inline and recover without a modal error dialog."""
    app, dialog, editor, main = preview
    _show(app, dialog)
    editor.kind.setCurrentText("polygon")
    editor.vertices.setPlainText("0 0\n1 1")
    _events(app)
    assert editor.footprint_plot.getImage("Beam density (1/m²)") is None
    assert "polygon" in editor.live_status.text()
    editor.vertices.setPlainText("-5 -3\n5 -3\n5 3\n-5 3")
    _events(app)
    assert editor.footprint_plot.getImage("Beam density (1/m²)") is not None
    main.fscan = None
    main.scanSelector.sigScanChanged.emit(None)
    _events(app)
    assert "Load a scan" in editor.live_status.text()
    assert editor.footprint_plot.getImage("Beam density (1/m²)") is None


def test_manual_preview_does_not_change_full_scan_alignment(preview):
    """An illustration override must not change the full-scan diagnostic curve."""
    app, dialog, editor, main = preview
    _show(app, dialog)
    editor.scan_preview.setChecked(True)
    editor._preview()
    expected = editor.overlap_plot.getCurve("H").getYData(copy=True)
    editor.live_angle_modes["incidence"].setCurrentIndex(1)
    editor.preview_alpha.setValue(1.2)
    editor._preview()
    np.testing.assert_array_equal(editor.overlap_plot.getCurve("H").getYData(),
                                  expected)

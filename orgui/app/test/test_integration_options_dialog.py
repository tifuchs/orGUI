"""Lifecycle regressions for the integration-correction options dialog."""

from types import SimpleNamespace

import numpy as np
import pytest
import silx
from silx.gui import qt
from silx.gui.plot import Plot2D
from silx.gui.plot.backends.BackendMatplotlib import BackendMatplotlib
from silx.io.dictdump import dicttonx, nxtodict

from orgui.app.config_data import ConfigData, CorrectionState
from orgui.app.QScanSelector import IntegrationOptionsDialog


@pytest.fixture(scope="session")
def qapp():
    """Keep a Qt application alive for the dialog tests."""
    application = qt.QApplication.instance()
    if application is None:
        application = qt.QApplication([])
    return application


class _Scan:
    """Minimal scan exposing one rate-like beam monitor."""

    auxillary_counters = ("ic2", "temperature")
    ic2 = np.array([10.0, 11.0, 12.0])
    temperature = np.array([295.0, 295.1, 295.2])
    exposure_time = np.array([1.0, 2.0, 3.0])

    def __len__(self):
        return 3


def _selector(main=None, *, rocking=False):
    """Return the correction state needed by the options dialog."""
    selector = SimpleNamespace(
        parentmainwindow=main,
        useMaskBox=qt.QCheckBox("Use pixel mask"),
        useSolidAngleBox=qt.QCheckBox("Solid angle correction"),
        usePolarizationBox=qt.QCheckBox("Polarization correction"),
        useLorentzBox=qt.QCheckBox("Lorentz correction"),
        useFootprintBox=qt.QCheckBox("Beam footprint"),
        useNormalizationBox=qt.QCheckBox("Normalize integrated intensities"),
    )
    if rocking:
        selector.scanstab = qt.QTabWidget()
        for name in ("stationary", "other", "rocking", "reflectivity"):
            selector.scanstab.addTab(qt.QWidget(), name)
        selector.scanstab.setCurrentIndex(2)
    return selector


def test_close_and_reopen_preserves_contents_and_settings(qapp):
    """Closing the persistent dialog only hides its live widget tree."""
    main_window = qt.QMainWindow()
    selector = _selector()
    dialog = IntegrationOptionsDialog(selector, parent=main_window)
    try:
        dialog.show()
        qapp.processEvents()
        groups = tuple(dialog.findChildren(qt.QGroupBox))
        selector.useSolidAngleBox.setChecked(True)

        assert not dialog.testAttribute(qt.Qt.WA_DeleteOnClose)
        assert {group.title() for group in groups} == {
            "Detector signal",
            "Incident beam",
            "CTR result",
        }

        dialog.close()
        qapp.processEvents()
        assert not dialog.isVisible()

        dialog.show()
        qapp.processEvents()
        assert dialog.isVisible()
        assert tuple(dialog.findChildren(qt.QGroupBox)) == groups
        assert all(group.isVisibleTo(dialog) for group in groups)
        assert selector.useSolidAngleBox.isChecked()
    finally:
        dialog.deleteLater()
        main_window.deleteLater()


@pytest.mark.parametrize("backend", ["matplotlib", "opengl"])
def test_hdf5_load_and_preview_reopen(qapp, tmp_path, monkeypatch, backend):
    """A preview retains its backend without changing its parent's surface."""
    if backend == "opengl":
        if qapp.platformName() in {"offscreen", "minimal"}:
            pytest.skip("OpenGL widget rendering requires a native Qt platform")
        context = qt.QOpenGLContext()
        if not context.create():
            pytest.skip("No OpenGL context available on this Qt platform")
        from silx.gui.plot.backends.BackendOpenGL import BackendOpenGL
        backend_class = BackendOpenGL
    else:
        backend_class = BackendMatplotlib
    monkeypatch.setattr(silx.config, "DEFAULT_PLOT_BACKEND", backend)
    # Installed-package tests cannot resolve examples relative to __file__.
    config_path = tmp_path / "config.ini"
    config_path.write_text(
        "[Machine]\nE = 77\nSDD = 0.78\npixelsize = 172e-6\n"
        "sizex = 1475\nsizey = 1679\ncpx = 731\ncpy = 1587.856\n"
        "[Lattice]\na1 = 2.774\na2 = 2.774\na3 = 6.796\n"
        "alpha1 = 90\nalpha2 = 90\nalpha3 = 120\n"
        "refractionindex = 1.1415e-6\n"
        "[Diffractometer]\nazimuthal_reference = 90\n"
        "polarization_axis = 0\npolarization_factor = 1\nmu = 0.1\n",
        encoding="utf-8",
    )
    config = ConfigData.from_ini(config_path)
    config.corrections = CorrectionState(
        sample_length_m=0.01, sample_width_m=0.01,
        beam_shape_analytical=True, beam_shape_name="Gaussian",
        beam_shape_values=(160.0,), use_footprint=True,
        sample_interception={
            "version": 1, "enabled": True,
            "shape": {"kind": "rectangle", "dimensions_m": [0.01, 0.01]},
            "azimuth_source": "fixed", "fixed_azimuth_deg": 0.0,
            "normal_rotation_confirmed": True,
            "reference_incidence_deg": 0.36,
            "horizontal": {
                "analytical": True, "shape": "Top hat", "shape_values": [20000],
            },
        },
    )
    filename = tmp_path / "database.h5"
    dicttonx({"configuration": config.to_nxdict(role="scan")}, filename)
    loaded = ConfigData.from_nxdict(nxtodict(filename)["configuration"])
    main = qt.QMainWindow()
    main.setAttribute(qt.Qt.WA_DontShowOnScreen, True)
    main.setCentralWidget(Plot2D(main))
    main.ubcalc = SimpleNamespace()
    selector = _selector(main)
    selector.set_integration_options = lambda values: (
        selector.useFootprintBox.setChecked(values["footprint"])
    )
    main.scanSelector = selector
    dialog = selector.correctionsDialog = IntegrationOptionsDialog(
        selector, parent=main
    )
    dialog.setAttribute(qt.Qt.WA_DontShowOnScreen, True)
    messages = []
    previous_handler = qt.qInstallMessageHandler(
        lambda kind, context, message: messages.append(message)
    )
    try:
        main.show()
        qapp.processEvents()
        dialog.show()
        qapp.processEvents()
        groups = tuple(dialog.findChildren(qt.QGroupBox))
        dialog.close()
        loaded.apply_to_gui(main)
        beam = dialog.footprintOptions_shared()
        assert beam.testAttribute(qt.Qt.WA_NativeWindow)
        assert not beam.sampleEditor.horizontal.testAttribute(qt.Qt.WA_NativeWindow)
        beam.setAttribute(qt.Qt.WA_DontShowOnScreen, True)
        beam.show()
        for index in range(beam.tabs.count()):
            beam.tabs.setCurrentIndex(index)
            qapp.processEvents()
        beam.sampleEditor._preview()
        qapp.processEvents()
        assert len(beam.sampleEditor.overlap_plot.getAllCurves()) == 2
        image = beam.sampleEditor.footprint_plot.getImage("Beam density (1/m²)")
        assert image is not None
        for plot in (
            beam.profilePlot, beam.sampleEditor.horizontal.profilePlot,
            beam.sampleEditor.overlap_plot, beam.sampleEditor.footprint_plot,
        ):
            assert isinstance(plot.getBackend(), backend_class)
            if backend == "opengl":
                gl_widgets = plot.findChildren(qt.QOpenGLWidget)
                assert gl_widgets and all(widget.isValid() for widget in gl_widgets[:1])
        colormap_action = beam.sampleEditor.footprint_plot.getColormapAction()
        colormap_dialog = colormap_action.getColormapDialog()
        colormap_dialog.setAttribute(qt.Qt.WA_DontShowOnScreen, True)
        colormap_action.trigger()
        qapp.processEvents()
        assert colormap_dialog.isVisible()
        assert colormap_dialog.getColormap() is image.getColormap()
        colormap_dialog.hide()
        beam.close()
        qapp.processEvents()
        for _ in range(2):
            dialog.show()
            qapp.processEvents()
            assert tuple(dialog.findChildren(qt.QGroupBox))[:3] == groups
            assert all(group.isVisibleTo(dialog) for group in groups)
            assert selector.useFootprintBox.isChecked()
            assert beam.L.value() == pytest.approx(10.0)
            assert beam.W.value() == pytest.approx(10.0)
            assert silx.config.DEFAULT_PLOT_BACKEND == backend
            if hasattr(qt, "QSurface"):
                assert dialog.windowHandle().surfaceType() == qt.QSurface.RasterSurface
            dialog.close()
        assert not any(
            "non-opengl surface" in message or "Failed to create QRhi" in message
            or "Failed to make context current" in message
            for message in messages
        )
    finally:
        qt.qInstallMessageHandler(previous_handler)
        dialog.deleteLater()
        main.deleteLater()
        qapp.sendPostedEvents(None, qt.QEvent.DeferredDelete)
        qapp.processEvents()


def test_new_session_uses_relative_total_flux_and_requires_horizontal_choice(qapp):
    """A new session is explicit about its relative and incomplete scale."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        assert dialog.scaleModeCombo.currentData() == "relative"
        assert main.ctr_correction_state.total_flux_calibrated is False
        assert not selector.useNormalizationBox.isChecked()
        assert dialog.scaleModeCombo.count() == 2
        assert dialog.scaleModeCombo.findData("legacy") == -1
        assert dialog.legacyNormalization.isHidden()
        assert not dialog.legacyNormalization.isEnabled()
        assert "Relative scale pending" in dialog.scaleStatus.text()
        assert "horizontal interception" in dialog.scaleStatus.text()
        assert "No primary monitor" == dialog.primaryMonitorCombo.currentText()
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_calibrated_rate_monitor_is_one_counter_with_exposure_once(qapp):
    """A rate-like monitor such as QM2 ic2 multiplies one frame exposure."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        selector.useNormalizationBox.setChecked(True)
        selector.useFootprintBox.setChecked(True)
        dialog.scaleModeCombo.setCurrentIndex(
            dialog.scaleModeCombo.findData("calibrated")
        )
        dialog.totalFluxEdit.setText("2.5e10")
        dialog.primaryMonitorCombo.setCurrentIndex(
            dialog.primaryMonitorCombo.findData("ic2")
        )
        dialog.monitorKindCombo.setCurrentIndex(
            dialog.monitorKindCombo.findData("rate")
        )
        dialog.monitorUnitEdit.setText("counts/s")
        dialog.referenceMonitorEdit.setText("1000")
        footprint = dialog.footprintOptions_shared()
        footprint.setHorizontalInterception("full")
        dialog._onTotalFluxChanged()

        state = main.ctr_correction_state
        assert state.total_incident_flux == pytest.approx(2.5e10)
        assert state.primary_monitor == "ic2"
        assert state.primary_monitor_kind == "rate"
        assert state.monitor_reference_exposure_s is None
        assert "exposure × M/Mref" in dialog.normalizationFormula.text()
        assert "Calibrated F² scale" in dialog.scaleStatus.text()
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_integrated_monitor_formula_does_not_apply_frame_exposure_twice(qapp):
    """Integrated counters use their own accumulated frame exposure."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        dialog.scaleModeCombo.setCurrentIndex(
            dialog.scaleModeCombo.findData("calibrated")
        )
        dialog.primaryMonitorCombo.setCurrentIndex(
            dialog.primaryMonitorCombo.findData("ic2")
        )
        dialog.monitorKindCombo.setCurrentIndex(
            dialog.monitorKindCombo.findData("integrated")
        )
        dialog._onTotalFluxChanged()

        assert "no second frame exposure" in dialog.normalizationFormula.text()
        assert dialog.referenceExposureEdit.isEnabled()
        assert "positive reference exposure" in dialog.scaleStatus.text()
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_legacy_multi_monitor_settings_are_preserved_but_hidden(qapp):
    """Opening an old configuration does not reinterpret its monitor product."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    main.ctr_correction_state = CorrectionState(
        normalize_exposure=True,
        monitor_corrections=("ic2", "temperature"),
    )
    main.reconstruction_normalize_exposure = True
    main.reconstruction_monitor_corrections = ("ic2", "temperature")
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        assert dialog.scaleModeCombo.currentIndex() == -1
        assert dialog.scaleModeCombo.findData("legacy") == -1
        assert dialog.legacyNormalization.isHidden()
        assert not dialog.legacyNormalization.isEnabled()
        assert dialog.monitorEdit.text() == "ic2, temperature"
        assert "Deprecated normalization preserved" in dialog.scaleStatus.text()
        assert dialog.primaryMonitorCombo.currentData() == ""
        dialog._onTotalFluxChanged()
        dialog.refresh()
        state = main.ctr_correction_state
        assert state.total_flux_calibrated is None
        assert state.monitor_corrections == ("ic2", "temperature")
        assert state.normalize_exposure is True
        footprint = dialog.footprintOptions_shared()
        assert footprint.beamFlux.isHidden()
        assert not footprint.beamFlux.isEnabled()
        assert footprint.legacyDensity.isHidden()
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_hidden_legacy_editor_cannot_change_settings_or_select_legacy(qapp):
    """Deprecated controls retain state without allowing GUI edits."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        assert not dialog.legacyNormalization.isEnabled()
        assert dialog.legacyNormalization.isHidden()
        dialog.normalizeExposureBox.setChecked(False)
        dialog.monitorEdit.setText("ic2, temperature")
        dialog._onNormalizationChanged()

        state = main.ctr_correction_state
        assert state.normalize_exposure is True
        assert state.monitor_corrections == ()
        assert not hasattr(main, "reconstruction_monitor_corrections")
        dialog.monitorEdit.setText("stale")
        dialog.refresh()
        assert dialog.monitorEdit.text() == ""
        assert dialog.scaleModeCombo.currentData() == "relative"

        dialog.scaleModeCombo.setCurrentIndex(
            dialog.scaleModeCombo.findData("calibrated")
        )
        dialog.scaleModeCombo.setCurrentIndex(
            dialog.scaleModeCombo.findData("relative")
        )
        assert state.total_flux_calibrated is False
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_imported_legacy_settings_can_explicitly_switch_to_total_flux(qapp):
    """An explicit mode selection migrates without reusing legacy monitors."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    main.ctr_correction_state = CorrectionState(
        normalize_exposure=False, monitor_corrections=("ic2", "temperature")
    )
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        dialog.scaleModeCombo.setCurrentIndex(
            dialog.scaleModeCombo.findData("relative")
        )
        state = main.ctr_correction_state
        assert state.total_flux_calibrated is False
        assert state.primary_monitor is None
        assert state.monitor_corrections == ("ic2", "temperature")
        assert state.normalize_exposure is False
        assert not selector.useNormalizationBox.isChecked()
        assert not dialog.legacyNormalization.isEnabled()
        assert "Deprecated" not in dialog.scaleStatus.text()
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_rocking_scan_reports_that_ctr_is_calculated_during_reduction(qapp):
    """The extraction dialog identifies the later rocking-reducer boundary."""
    main = qt.QMainWindow()
    selector = _selector(main, rocking=True)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        assert "calculated after rocking-curve integration" in (
            dialog.measurementModeInfo.text()
        )
        assert "angular acceptance" in dialog.measurementModeInfo.text()
    finally:
        dialog.deleteLater()
        main.deleteLater()

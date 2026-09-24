"""Lifecycle regressions for the integration-correction options dialog."""

from types import SimpleNamespace

import numpy as np
import pytest
from silx.gui import qt

from orgui.app.config_data import CorrectionState
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


def test_new_session_uses_relative_total_flux_and_requires_horizontal_choice(qapp):
    """A new session is explicit about its relative and incomplete scale."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        assert dialog.scaleModeCombo.currentData() == "relative"
        assert main.ctr_correction_state.total_flux_calibrated is False
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


def test_legacy_multi_monitor_settings_remain_separate_and_labeled(qapp):
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
        assert dialog.scaleModeCombo.currentData() == "legacy"
        assert dialog.monitorEdit.text() == "ic2, temperature"
        assert "Legacy density/area convention" in dialog.scaleStatus.text()
        assert dialog.primaryMonitorCombo.currentData() == ""
    finally:
        dialog.deleteLater()
        main.deleteLater()


def test_monitor_editor_updates_the_shared_correction_state(qapp):
    """The one editor owns the monitor selection captured by mapping."""
    main = qt.QMainWindow()
    main.fscan = _Scan()
    selector = _selector(main)
    dialog = IntegrationOptionsDialog(selector, parent=main)
    try:
        assert not dialog.legacyNormalization.isEnabled()
        dialog.scaleModeCombo.setCurrentIndex(
            dialog.scaleModeCombo.findData("legacy")
        )
        assert dialog.legacyNormalization.isEnabled()
        dialog.normalizeExposureBox.setChecked(False)
        dialog.monitorEdit.setText("ic2, temperature")
        dialog._onNormalizationChanged()

        state = main.ctr_correction_state
        assert state.normalize_exposure is False
        assert state.monitor_corrections == ("ic2", "temperature")
        assert not hasattr(main, "reconstruction_monitor_corrections")
        dialog.monitorEdit.setText("stale")
        dialog.refresh()
        assert dialog.monitorEdit.text() == "ic2, temperature"
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

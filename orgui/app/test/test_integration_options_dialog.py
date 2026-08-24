"""Lifecycle regressions for the integration-correction options dialog."""

from types import SimpleNamespace

import pytest
from silx.gui import qt

from orgui.app.QScanSelector import IntegrationOptionsDialog


@pytest.fixture(scope="session")
def qapp():
    """Keep a Qt application alive for the dialog tests."""
    application = qt.QApplication.instance()
    if application is None:
        application = qt.QApplication([])
    return application


def _selector():
    """Return the correction state needed by the options dialog."""
    return SimpleNamespace(
        parentmainwindow=None,
        useMaskBox=qt.QCheckBox("Use pixel mask"),
        useSolidAngleBox=qt.QCheckBox("Solid angle correction"),
        usePolarizationBox=qt.QCheckBox("Polarization correction"),
        useLorentzBox=qt.QCheckBox("Lorentz and rod interception"),
        useFootprintBox=qt.QCheckBox("Beam footprint"),
        useNormalizationBox=qt.QCheckBox("Normalize integrated intensities"),
    )


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
            "Detector corrections",
            "Geometrical corrections",
            "Exposure and monitor normalization",
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

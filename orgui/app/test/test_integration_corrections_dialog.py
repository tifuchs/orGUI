"""Tests for the footprint-correction options dialog.

:class:`orgui.app.peak1Dintegr.IntegrationCorrectionsDialog` owns the user
side of the numerical active-area correction: it converts the
millimeter and micrometer values shown to the user into the meters
:mod:`orgui.datautils.xrayutils.beamprofile` works in, and decides whether an
integration runs against an analytical beam shape or a measured profile.

These tests pin that boundary. They construct the dialog directly, without
showing it, so no user interaction is involved and no plot widget is built.
"""

import numpy as np
import pytest
from silx.gui import qt

from orgui.app.peak1Dintegr import BEAM_SHAPES, IntegrationCorrectionsDialog
from orgui.datautils.xrayutils.beamprofile import (
    DistributionBeamProfile,
    GaussianBeamProfile,
    MeasuredBeamProfile,
)

#: Sample size along the beam, in meters, matching the dialog default of 5 mm.
L = 5e-3

#: Angles at which corrections are compared, in radian.
ALPHAS = np.deg2rad(np.array([0.1, 0.5, 1.0, 5.0]))


@pytest.fixture(scope="session")
def qapp():
    """The Qt application, kept referenced for the whole test session.

    A ``QApplication`` that is not held on to is garbage-collected, and
    creating a widget afterwards aborts the interpreter.
    """
    application = qt.QApplication.instance()
    if application is None:
        application = qt.QApplication([])
    return application


@pytest.fixture
def dialog(qapp):
    """An unshown corrections dialog on a running Qt application."""
    widget = IntegrationCorrectionsDialog()
    yield widget
    widget.deleteLater()


def _write_gaussian_profile(path, fwhm, n=2001):
    """Write a Gaussian profile file with its position column in mm."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    z = np.linspace(-8 * sigma, 8 * sigma, n)
    intensity = np.exp(-0.5 * (z / sigma) ** 2)
    np.savetxt(path, np.column_stack([z * 1e3, intensity]), header="z/mm  profile")
    return z, intensity


def _select_shape(dialog, name):
    """Select an analytical shape by name and return its parameter spin boxes."""
    dialog.shapeSelector.setCurrentIndex(dialog.shapeSelector.findText(name))
    shape = dialog.currentShape()
    assert shape.name == name
    return dialog.shapeParameters[: len(shape.parameters)]


def test_default_reproduces_the_original_gaussian_correction(dialog):
    """An untouched dialog gives the correction orGUI has always applied.

    The default path now runs through the generic distribution wrapper, so
    what has to hold is that its numbers still match the closed form, not
    that it is the same class.
    """
    assert dialog.analyticalButton.isChecked()
    assert dialog.currentShape().name == "Gaussian"

    profile = dialog.beamProfile()
    # The spin box shows micrometers, the profile is built in meters.
    reference = GaussianBeamProfile(dialog.shapeParameters[0].value() * 1e-6)

    np.testing.assert_allclose(
        profile.corrections(ALPHAS, L), reference.corrections(ALPHAS, L), rtol=1e-12
    )


def test_every_shape_builds_a_usable_profile(dialog):
    """Each offered shape produces finite corrections from its defaults."""
    for shape in BEAM_SHAPES:
        _select_shape(dialog, shape.name)

        flux, area = dialog.beamProfile().corrections(ALPHAS, L)

        assert np.all(np.isfinite(flux)), shape.name
        assert np.all(np.isfinite(area)), shape.name
        assert np.all(flux > 0) and np.all(flux <= 1.0 + 1e-12), shape.name
        assert np.all(area > 0) and np.all(area <= 1.0 + 1e-12), shape.name


def test_only_the_selected_shapes_parameters_are_shown(dialog):
    """The numeric controls follow the selected shape."""
    for shape in BEAM_SHAPES:
        _select_shape(dialog, shape.name)

        for index, spin in enumerate(dialog.shapeParameters):
            label = dialog.shapeParameterLabels[index]
            used = index < len(shape.parameters)
            assert spin.isVisibleTo(dialog) is used, (shape.name, index)
            assert label.isVisibleTo(dialog) is used, (shape.name, index)
            if used:
                assert label.text() == shape.parameters[index].label


def test_top_hat_matches_the_analytical_slit_correction(dialog):
    """A uniform beam gives ``C_flux = L sin(alpha) / width``.

    This is the linear geometrical correction of Gibaud, Vignaud & Sinha
    (1993), *Acta Cryst.* A49, 642, equation (12), and it holds until the
    sample intercepts the whole beam.
    """
    width_um = 400.0
    (width,) = _select_shape(dialog, "Top hat")
    width.setValue(width_um)

    profile = dialog.beamProfile()
    alpha = np.deg2rad(np.array([0.1, 0.5, 1.0, 2.0]))
    expected = np.minimum(L * np.sin(alpha) / (width_um * 1e-6), 1.0)

    np.testing.assert_allclose(profile.flux_on_sample(alpha, L), expected, rtol=1e-9)
    # The whole footprint sits inside a flat beam, so it is fully lit.
    np.testing.assert_allclose(
        profile.illuminated_area_fraction(alpha[expected < 1.0], L), 1.0, rtol=1e-9
    )


def test_shape_parameters_are_micrometers(dialog):
    """Width controls shown in micrometers reach the profile as meters."""
    (fwhm,) = _select_shape(dialog, "Gaussian")
    fwhm.setValue(150.0)

    assert dialog.beamProfile().fwhm == pytest.approx(150e-6, rel=1e-4)


def test_shape_values_survive_switching_shapes(dialog):
    """Editing one shape, visiting another and coming back keeps the value."""
    (width, _flat) = _select_shape(dialog, "Trapezoid")
    width.setValue(77.0)

    _select_shape(dialog, "Gaussian")
    (width, _flat) = _select_shape(dialog, "Trapezoid")

    assert width.value() == pytest.approx(77.0)


def test_offset_applies_to_analytical_shapes(dialog):
    """Displacing the sample reduces the flux a symmetric beam delivers."""
    (fwhm,) = _select_shape(dialog, "Gaussian")
    fwhm.setValue(200.0)
    alpha = np.deg2rad(0.5)

    centered = dialog.beamProfile().flux_on_sample(alpha, L)
    dialog.profileOffset.setValue(60.0)
    displaced = dialog.beamProfile().flux_on_sample(alpha, L)

    assert displaced < centered
    assert dialog.beamProfile().sample_center == pytest.approx(60e-6)


def test_measured_profile_selected_without_a_file_is_an_error(dialog):
    """Integrating without a loaded profile must fail loudly, not silently."""
    dialog.measuredButton.setChecked(True)

    with pytest.raises(ValueError, match="No beam profile loaded"):
        dialog.beamProfile()


def test_unreadable_profile_file_is_reported_and_leaves_no_profile(dialog, tmp_path):
    """A malformed file is rejected while the dialog is open."""
    broken = tmp_path / "broken.dat"
    broken.write_text("this is not a profile\n")
    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(broken))

    assert dialog.loadProfile() is False
    with pytest.raises(ValueError, match="No beam profile loaded"):
        dialog.beamProfile()


def test_stale_profile_is_dropped_when_a_later_load_fails(dialog, tmp_path):
    """A failed reload must not leave the previous profile in place."""
    good = tmp_path / "good.dat"
    _write_gaussian_profile(good, 100e-6)
    broken = tmp_path / "broken.dat"
    broken.write_text("nonsense\n")

    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(good))
    assert dialog.loadProfile() is True

    dialog.profileFileEdit.setText(str(broken))
    assert dialog.loadProfile() is False
    with pytest.raises(ValueError, match="No beam profile loaded"):
        dialog.beamProfile()


def test_loaded_profile_reproduces_the_gaussian_corrections(dialog, tmp_path):
    """A Gaussian read through the dialog matches the analytical path.

    This is the unit-conversion check of the whole GUI boundary: file in
    millimeters, spin boxes in micrometers, corrections in meters.
    """
    fwhm = 120e-6
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, fwhm)

    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    assert dialog.loadProfile() is True

    measured = dialog.beamProfile()
    assert isinstance(measured, MeasuredBeamProfile)
    assert measured.fwhm == pytest.approx(fwhm, rel=1e-3)

    np.testing.assert_allclose(
        measured.corrections(ALPHAS, L),
        GaussianBeamProfile(fwhm).corrections(ALPHAS, L),
        rtol=1e-5,
    )


def test_position_unit_selection_rescales_the_profile(dialog, tmp_path):
    """The position column is read in the unit chosen in the dialog."""
    fwhm = 120e-6
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, fwhm)

    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    dialog.profileUnit.setCurrentIndex(dialog.profileUnit.findText("mm"))
    assert dialog.loadProfile() is True
    assert dialog.beamProfile().fwhm == pytest.approx(fwhm, rel=1e-3)

    # The same file read as micrometers describes a beam 1000x narrower.
    dialog.profileUnit.setCurrentIndex(dialog.profileUnit.findText("microns"))
    assert dialog.loadProfile() is True
    assert dialog.beamProfile().fwhm == pytest.approx(fwhm * 1e-3, rel=1e-3)


def test_offset_spin_box_is_micrometers(dialog, tmp_path):
    """The sample offset shown in micrometers arrives in meters."""
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, 100e-6)
    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    assert dialog.loadProfile() is True

    unshifted = dialog.beamProfile().sample_center
    dialog.profileOffset.setValue(25.0)

    assert dialog.beamProfile().sample_center == pytest.approx(unshifted + 25e-6)


def test_mode_switch_enables_the_matching_widgets(dialog):
    """Only the widgets of the selected beam model are editable."""
    assert dialog.shapeGroup.isEnabled()
    assert not dialog.profileGroup.isEnabled()

    dialog.measuredButton.setChecked(True)

    assert not dialog.shapeGroup.isEnabled()
    assert dialog.profileGroup.isEnabled()


def test_preview_summary_reports_width_and_centroid(dialog, tmp_path):
    """The summary line describes whichever profile is selected."""
    (fwhm,) = _select_shape(dialog, "Gaussian")
    fwhm.setValue(250.0)
    assert "FWHM 250.0 microns" in dialog.profileInfo.text()

    # An asymmetric measured profile puts its centroid off the peak, so
    # centering on the peak has to show a non-zero centroid.
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, 100e-6)
    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    assert dialog.loadProfile() is True

    assert "FWHM 100.0 microns" in dialog.profileInfo.text()
    assert "centroid at" in dialog.profileInfo.text()


def test_preview_summary_reports_a_missing_profile(dialog):
    """Selecting a measured profile without a file says so in the preview."""
    dialog.measuredButton.setChecked(True)

    assert "No beam profile loaded" in dialog.profileInfo.text()


def test_centroid_position_follows_the_centering_choice(dialog, tmp_path):
    """The centroid the preview marks moves with the centering reference."""
    z = np.linspace(-200e-6, 200e-6, 401)
    intensity = np.exp(-0.5 * ((z - 20e-6) / 30e-6) ** 2) + 0.4 * np.exp(
        -0.5 * ((z + 80e-6) / 60e-6) ** 2
    )
    path = tmp_path / "skewed.dat"
    np.savetxt(path, np.column_stack([z * 1e3, intensity]))

    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    assert dialog.loadProfile() is True

    dialog.profileCenter.setCurrentIndex(dialog.profileCenter.findText("centroid"))
    assert dialog.beamProfile().centroid_position == pytest.approx(0.0, abs=1e-12)

    dialog.profileCenter.setCurrentIndex(dialog.profileCenter.findText("peak"))
    # Centering on the peak leaves the centroid on the weak satellite side.
    assert dialog.beamProfile().centroid_position < -1e-6


def test_cancel_restores_the_last_accepted_settings(dialog, tmp_path):
    """Cancelling reverts every field, as the dialog did for L and beam."""
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, 100e-6)
    dialog.L.setValue(3.0)
    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    dialog.loadProfile()
    dialog.onOk()

    dialog.L.setValue(7.5)
    dialog.profileOffset.setValue(40.0)
    dialog.profileCenter.setCurrentIndex(dialog.profileCenter.findText("peak"))
    _select_shape(dialog, "Top hat")
    dialog.analyticalButton.setChecked(True)
    dialog.onCancel()

    assert dialog.L.value() == pytest.approx(3.0)
    assert dialog.profileOffset.value() == pytest.approx(0.0)
    assert dialog.profileCenter.currentText() == "centroid"
    assert dialog.currentShape().name == "Gaussian"
    assert dialog.measuredButton.isChecked()
    # The restored state is still usable, not just cosmetically reset.
    assert isinstance(dialog.beamProfile(), MeasuredBeamProfile)


def test_settings_round_trip(dialog, tmp_path):
    """:meth:`settings` and :meth:`setSettings` are inverse."""
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, 100e-6)
    dialog.L.setValue(2.25)
    dialog.measuredButton.setChecked(True)
    dialog.profileFileEdit.setText(str(path))
    dialog.loadProfile()
    (base, flat) = _select_shape(dialog, "Trapezoid")
    base.setValue(90.0)
    flat.setValue(30.0)
    dialog.profileCenter.setCurrentIndex(dialog.profileCenter.findText("median"))
    dialog.profileOffset.setValue(-12.5)
    saved = dialog.settings()

    restored = IntegrationCorrectionsDialog()
    try:
        restored.setSettings(saved)
        assert restored.settings() == saved
        assert restored.beamProfile().sample_center == pytest.approx(
            dialog.beamProfile().sample_center
        )
    finally:
        restored.deleteLater()


def test_choosing_a_profile_file_selects_it_as_the_beam(dialog, tmp_path):
    """A loaded profile must not sit behind an analytical shape.

    Loading a profile while the analytical shape stayed selected left the
    default 20 micron Gaussian correcting the data. Against a 10 mm sample
    that is a wholly different correction -- it falls monotonically as
    1/sin(alpha) instead of peaking where the beam fills the sample -- with
    nothing in the result to show which beam was used.
    """
    path = tmp_path / "gauss.dat"
    _write_gaussian_profile(path, 250e-6)
    assert dialog.analyticalButton.isChecked()

    dialog.profileFileEdit.setText(str(path))
    dialog._onProfileFileChosen()

    assert dialog.measuredButton.isChecked()
    assert isinstance(dialog.beamProfile(), MeasuredBeamProfile)
    assert dialog.beamProfile().fwhm == pytest.approx(250e-6, rel=1e-3)


def test_an_unreadable_profile_does_not_switch_the_beam(dialog, tmp_path):
    """A file that fails to load leaves the analytical shape in charge."""
    broken = tmp_path / "broken.dat"
    broken.write_text("not a profile\n")

    dialog.profileFileEdit.setText(str(broken))
    dialog._onProfileFileChosen()

    assert dialog.analyticalButton.isChecked()
    assert isinstance(dialog.beamProfile(), DistributionBeamProfile)

"""Unit boundaries, provenance, frame policy and safe shape replacement."""

import ast
import json
import subprocess
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from silx.gui import qt
from silx.io.dictdump import dicttonx

from orgui.app.config_data import (
    CURVE_CORRECTIONS_GROUP,
    ConfigData,
    CorrectionState,
    CurveCorrectionRecord,
    corrections_from_nxdict,
    corrections_to_nxdict,
    curve_correction_record_from_nxdict,
    curve_correction_record_to_nxdict,
)
from orgui.app.integration_corrections import (
    corrected_curve_from_record,
    frame_correction_policy,
)
from orgui.app.peak1Dintegr import IntegrationCorrectionsDialog
from orgui.app.sample_interception_config import (
    embed_profile,
    profile_from_settings,
    replacement_shape_divisor,
    shape_frame_factors,
    sample_angle_inputs,
    parallel_omega_reading,
    orientation_for_parallel_reading,
)
from orgui.backend.scans import sample_azimuth
from orgui.datautils.xrayutils.corrections import beamprofile as bp


@pytest.fixture(scope="module")
def qapp():
    """Keep the shared Qt application alive for widget regression checks."""
    return qt.QApplication.instance() or qt.QApplication([])


def _settings():
    return {
        "version": 1,
        "enabled": True,
        "shape": {"kind": "rectangle", "dimensions_m": [0.01, 0.01]},
        "offset_m": [0, 0],
        "reference_incidence_deg": 0.36,
        "normal_rotation_confirmed": True,
        "azimuth_source": "phi",
        "azimuth_sign": 1,
        "horizontal": {"analytical": True, "shape": "Top hat", "shape_values": [20000]},
    }


@pytest.mark.parametrize("total_flux", [False, True])
def test_shape_policy_reconstructs_known_photon_curve(total_flux):
    """An independently specified square reference recovers signal and variance."""
    settings = _settings()
    state = CorrectionState(sample_interception=settings)
    if total_flux:
        state.total_flux_calibrated = True
        state.total_incident_flux = 3e10
    scan = SimpleNamespace(phi=np.array([0, 45]), exposure_time=np.array([1.0, 2.0]))
    alpha = np.deg2rad([0.36, 0.36])
    profile = bp.gaussian_profile(160e-6)
    policy = frame_correction_policy(
        scan,
        state,
        2,
        use_normalization=True,
        use_illumination=True,
        alpha=alpha,
        beam_profile=profile,
    )
    reference_h = np.array([0.178090128539, 0.178155252900]) / np.sin(alpha)
    if not total_flux:
        reference_h /= profile.peak_density * 50 * 0.01**2
    np.testing.assert_allclose(policy.illumination_divisor, reference_h, rtol=1e-10)
    assert policy.vertical_intercepted_fraction is None
    assert policy.horizontal_intercepted_fraction is None
    q = policy.normalization_divisor
    base = q * reference_h * 7
    record = CurveCorrectionRecord(
        algorithm="shape_interception_total_flux_v1"
        if total_flux
        else "shape_interception_legacy_v1",
        output_quantity="curve",
        scale_convention=policy.scale_convention,
        normalization_status="applied",
        normalization_divisor=q,
        illumination_status="applied",
        illumination_divisor=policy.illumination_divisor,
        illumination_convention=policy.illumination_convention,
        base_croibg=base,
        base_croibg_variance=(q * reference_h * 0.3) ** 2,
        profile_provenance=policy.interception_provenance,
        alpha=alpha,
    )
    signal, errors, _, _ = corrected_curve_from_record(record)
    np.testing.assert_allclose(signal, 7, rtol=1e-10)
    np.testing.assert_allclose(errors, 0.3, rtol=1e-10)
    settings["orientation_deg"] = 45
    replacement = replacement_shape_divisor(settings, profile, record)
    one = corrected_curve_from_record(
        record,
        footprint_action="apply",
        replacement_illumination=replacement,
        replacement_convention=record.illumination_convention,
    )[0]
    two = corrected_curve_from_record(
        record,
        footprint_action="apply",
        replacement_illumination=replacement,
        replacement_convention=record.illumination_convention,
    )[0]
    np.testing.assert_array_equal(one, two)
    settings["azimuth_source"] = "omega"
    with pytest.raises(ValueError, match="re-extraction"):
        replacement_shape_divisor(settings, profile, record)


def test_angle_mapping_missing_profile_and_frame_exclusion():
    """Declared acquisition units/order are retained; invalid frames are gaps."""
    settings = _settings()
    scan = SimpleNamespace(phi=[350.0, 360.0, 370.0])
    np.testing.assert_allclose(sample_azimuth(scan, settings, 3), np.deg2rad(scan.phi))
    with pytest.raises(ValueError, match="source"):
        sample_azimuth(scan, {"azimuth_source": "missing"}, 3)
    with pytest.raises(ValueError, match="per frame"):
        sample_azimuth(scan, settings, 4)
    settings["azimuth_sign"] = -1
    h, fraction, _, valid, _ = shape_frame_factors(
        settings,
        bp.gaussian_profile(160e-6),
        np.deg2rad([0, 0.36]),
        np.deg2rad([0, 45]),
    )
    np.testing.assert_array_equal(valid, [False, True])
    assert np.isnan(h[0]) and np.isnan(fraction[0])
    settings["horizontal"] = {}
    with pytest.raises(ValueError, match="horizontal"):
        shape_frame_factors(settings, bp.gaussian_profile(160e-6), 0.01, 0)


def test_embedded_profile_preserves_alignment_and_does_not_need_file(tmp_path):
    """Recentring normalized embedded coordinates must not erase a sample offset."""
    positions = np.array([-0.002, -0.001, 0, 0.001])
    raw = bp.MeasuredBeamProfile(
        positions, [1, 2, 4, 1], center="median", offset=0.0004
    )
    settings = embed_profile(
        {
            "analytical": False,
            "profile_center": "median",
            "profile_offset": 400,
            "profile_file": "missing.dat",
        },
        raw,
    )
    rebuilt = profile_from_settings(settings)
    np.testing.assert_allclose(
        raw.density_at(raw.z), rebuilt.density_at(raw.z), rtol=1e-13
    )
    model = _settings()
    model["horizontal"] = settings
    state = CorrectionState(sample_interception=model)
    path = tmp_path / "shape.nxs"
    dicttonx({"corrections": corrections_to_nxdict(state)}, path)
    with h5py.File(path) as file:
        loaded = corrections_from_nxdict(file["corrections"])
    assert loaded.sample_interception == model
    result = shape_frame_factors(
        loaded.sample_interception, bp.gaussian_profile(160e-6), 0.01, 0.3
    )
    assert np.isfinite(result[0])


def test_old_reader_refuses_shape_record_and_new_reader_keeps_lazy_arrays(tmp_path):
    """Exercise the actual prior reducer guard, not a simulated version comparison."""
    source = subprocess.check_output(
        ["git", "show", "HEAD:orgui/app/peak1Dintegr.py"], text=True
    )
    tree = ast.parse(source)
    cls = next(
        n
        for n in tree.body
        if isinstance(n, ast.ClassDef) and n.name == "RockingPeakIntegrator"
    )
    function = next(
        n
        for n in cls.body
        if isinstance(n, ast.FunctionDef) and n.name == "_legacy_rocking_curve_group"
    )
    function.decorator_list = []
    namespace = {
        "CURVE_CORRECTIONS_GROUP": CURVE_CORRECTIONS_GROUP,
        "CURVE_CORRECTIONS_SCHEMA_VERSION": 3,
    }
    exec(
        compile(
            ast.Module(body=[function], type_ignores=[]), "old_reducer_guard", "exec"
        ),
        namespace,
    )
    record = CurveCorrectionRecord(
        algorithm="shape_interception_total_flux_v1",
        output_quantity="curve",
        scale_convention="total_flux_relative",
        base_croibg=np.ones((2, 3)),
        base_croibg_variance=np.ones((2, 3)),
        profile_provenance={"sample_interception_json": json.dumps(_settings())},
    )
    path = tmp_path / "curve.nxs"
    dicttonx(
        {
            "curve": {
                CURVE_CORRECTIONS_GROUP: curve_correction_record_to_nxdict(record),
                "rois": {"croibg": np.ones((2, 3))},
            }
        },
        path,
    )
    with h5py.File(path) as file:
        with pytest.raises(ValueError, match="must not reinterpret"):
            namespace[function.name](file["curve"])
        loaded = curve_correction_record_from_nxdict(
            file["curve"][CURVE_CORRECTIONS_GROUP]
        )
        assert isinstance(loaded.base_croibg, h5py.Dataset)
        assert loaded.profile_provenance == record.profile_provenance


def test_shape_editor_units_and_embedded_horizontal_restore():
    """UI dimensions cross to SI once; saved profiles work after relocation."""
    app = qt.QApplication.instance() or qt.QApplication([])
    parent = IntegrationCorrectionsDialog()
    parent.setSettings({"sample_interception": _settings()})
    dialog = parent.sampleEditor
    dialog.values["offset_x"].setValue(2)
    model = dialog.settings()
    assert model["offset_m"] == [0.002, 0]
    assert model["shape"]["dimensions_m"] == [0.01, 0.01]
    assert model["horizontal"]["shape_values"] == [20000]
    measured = bp.MeasuredBeamProfile(
        [-0.002, 0, 0.002], [1, 4, 1], center="peak", offset=0.0002
    )
    horizontal = embed_profile(
        {
            "analytical": False,
            "profile_center": "peak",
            "profile_offset": 200,
            "profile_file": "missing.dat",
        },
        measured,
    )
    dialog.horizontal.setSettings(horizontal)
    np.testing.assert_allclose(dialog.horizontal.measuredProfile().z, measured.z)
    np.testing.assert_allclose(
        dialog.horizontal.measuredProfile().density, measured.density
    )
    # Changing the offset moves the sample, rather than re-centring saved data.
    dialog.horizontal.profileOffset.setValue(300)
    np.testing.assert_allclose(
        dialog.horizontal.measuredProfile().z, measured.z - 0.0001
    )
    dialog.enabled.setChecked(False)
    assert dialog.settings()["horizontal"]["profile_offset"] == 300
    assert dialog.settings()["horizontal"]["positions_m"]
    dialog.close()
    parent.close()
    assert app is not None


def test_ini_json_old_defaults_and_headless_profile(tmp_path):
    """Optional shape INI settings and JSON snapshots preserve the same model."""
    from pathlib import Path

    template = Path("examples/config_minimal").read_text()
    model = _settings()
    config_path = tmp_path / "shape.ini"
    config_path.write_text(
        template + "\n[SampleInterception]\nsettings = " + json.dumps(model)
    )
    config = ConfigData.from_ini(config_path)
    state = CorrectionState.from_dict(
        json.loads(json.dumps(config.corrections.to_dict()))
    )
    assert state.sample_interception == model
    old_path = tmp_path / "old.ini"
    old_path.write_text(template)
    old = ConfigData.from_ini(old_path)
    assert old.corrections.sample_interception == {}
    assert corrections_to_nxdict(old.corrections)["@orgui_schema_version"] == 3
    state.beam_shape_analytical = True
    state.beam_shape_name = "Gaussian"
    state.beam_shape_values = (160,)
    policy = frame_correction_policy(
        SimpleNamespace(phi=[0, 45]),
        state,
        2,
        use_normalization=False,
        use_illumination=True,
        alpha=np.deg2rad(0.36),
    )
    np.testing.assert_allclose(
        policy.intercepted_fraction, [0.178090128539, 0.178155252900], atol=1e-11
    )


def test_source_frame_preview_uses_acquisition_angles():
    """Source-frame diagnostics use declared readbacks and selected frame incidence."""
    app = qt.QApplication.instance() or qt.QApplication([])
    host = qt.QWidget()
    host._mainWindow = lambda: SimpleNamespace(
        fscan=SimpleScan(),
        getMuOm=lambda: (np.deg2rad([0.36, 0.5]), 0),
    )
    parent = IntegrationCorrectionsDialog(host)
    parent.setSettings({
        "analytical": True, "shape": "Gaussian", "shape_values": [160],
        "sample_interception": _settings(),
    })
    dialog = parent.sampleEditor
    dialog.scan_preview.setChecked(True)
    dialog.frame_index.setValue(1)
    dialog._preview()
    assert "0.5 deg" in dialog.status.text()
    assert "45 deg" in dialog.status.text()
    assert dialog.overlap_plot.getGraphXLabel() == "Source frame index"
    assert dialog.tabs.isEnabled()
    dialog.close()
    parent.close()
    host.close()
    assert app is not None


class SimpleScan:
    """Minimal ordered angle source for the diagnostics widget."""

    phi = np.array([0, 45])

    def __len__(self):
        return 2


@pytest.mark.parametrize("axisname", ["mu", "th"])
def test_automatic_angles_follow_acquisition_convention(axisname):
    """Mu varies incidence; theta varies omega with the established minus sign."""
    settings = _settings()
    settings.pop("azimuth_source")
    settings.pop("reference_incidence_deg")
    settings.pop("normal_rotation_confirmed")
    scan = SimpleNamespace(axisname=axisname, axis=np.array([1., 2.]), th=30.)
    mu = np.deg2rad(scan.axis if axisname == "mu" else 0.36)
    alpha, omega, resolved = sample_angle_inputs(scan, settings, mu, count=2)
    np.testing.assert_allclose(alpha, mu if axisname == "mu" else [mu, mu])
    expected = -np.deg2rad(scan.axis if axisname == "th" else [30., 30.])
    np.testing.assert_allclose(omega, expected)
    assert resolved["reference_incidence_deg"] == pytest.approx(
        1. if axisname == "mu" else 0.36
    )
    assert "reference_incidence_deg" not in settings
    # Supplied getMuOm values take precedence over backend fallback metadata.
    _, omega, _ = sample_angle_inputs(scan, settings, mu, count=2, omega=[0.1, 0.2])
    np.testing.assert_allclose(omega, [0.1, 0.2])


def test_automatic_policy_preserves_rocking_curve_angle_arrays():
    """Multiple ROI curves retain their per-frame incidence and overlap arrays."""
    settings = _settings()
    settings["azimuth_source"] = "auto"
    alpha = np.deg2rad([[0.36, 0.5], [0.36, 0.5]])
    policy = frame_correction_policy(
        SimpleNamespace(), CorrectionState(sample_interception=settings), 2,
        use_normalization=False, use_illumination=True,
        alpha=alpha, omega=np.deg2rad([0, 45]),
        beam_profile=bp.gaussian_profile(160e-6),
    )
    assert policy.illumination_divisor.shape == (2, 2)
    np.testing.assert_array_equal(policy.interception_provenance["sample_incidence_rad"],
                                  alpha)


@pytest.mark.parametrize("unit", ["deg", "rad"])
def test_manual_counter_overrides_and_saved_incidence_replacement(unit):
    """Footprint overrides are saved separately from diffraction incidence."""
    settings = _settings()
    settings.pop("reference_incidence_deg")
    settings["incidence_source"] = "measured_alpha"
    settings["incidence_unit"] = unit
    settings["azimuth_unit"] = unit
    settings["offset_m"] = [0.001, 0.002]
    incidence = np.deg2rad([0.5, 0.8])
    azimuth = np.deg2rad([5, 40])
    scan = SimpleNamespace(
        measured_alpha=np.rad2deg(incidence) if unit == "deg" else incidence,
        phi=np.rad2deg(azimuth) if unit == "deg" else azimuth,
    )
    original_alpha = np.deg2rad([0.36, 0.36])
    vertical = bp.gaussian_profile(160e-6)
    policy = frame_correction_policy(
        scan, CorrectionState(sample_interception=settings), 2,
        use_normalization=False, use_illumination=True,
        alpha=original_alpha, beam_profile=vertical,
    )
    np.testing.assert_allclose(policy.interception_provenance["sample_incidence_rad"],
                               incidence)
    np.testing.assert_allclose(policy.interception_provenance["sample_azimuth_rad"],
                               azimuth)
    resolved = json.loads(policy.interception_provenance["sample_interception_json"])
    assert resolved["reference_incidence_deg"] == pytest.approx(0.5)
    record = CurveCorrectionRecord(
        algorithm="shape_interception_legacy_v1", output_quantity="curve",
        alpha=original_alpha, scale_convention=policy.scale_convention,
        profile_provenance=policy.interception_provenance,
    )
    np.testing.assert_allclose(replacement_shape_divisor(settings, vertical, record),
                               policy.illumination_divisor)
    settings["incidence_source"] = "other_counter"
    with pytest.raises(ValueError, match="incidence source requires re-extraction"):
        replacement_shape_divisor(settings, vertical, record)
    settings["incidence_source"] = "fixed"
    settings["fixed_incidence_deg"] = 1
    assert np.all(np.isfinite(replacement_shape_divisor(settings, vertical, record)))


def test_fixed_overrides_do_not_require_scan_counters():
    """Fixed incidence and azimuth are degrees regardless of counter units."""
    settings = _settings()
    settings.update(incidence_source="fixed", fixed_incidence_deg=0.7,
                    azimuth_source="fixed", fixed_azimuth_deg=32,
                    incidence_unit="rad", azimuth_unit="rad")
    alpha, azimuth, _ = sample_angle_inputs(SimpleNamespace(), settings, 0.01, count=3)
    np.testing.assert_allclose(alpha, np.deg2rad([0.7] * 3))
    np.testing.assert_allclose(azimuth, np.deg2rad([32] * 3))


def test_reference_and_rotation_offset_preserve_distinct_off_axis_geometry():
    """An angular mounting offset does not rotate the sample centre's orbit."""
    from orgui.app.sample_interception_config import shape_from_settings

    settings = _settings()
    settings.update(reference_azimuth_deg=20, orientation_deg=30,
                    offset_m=[0.002, 0.001])
    shape = shape_from_settings(settings)
    reference = shape.outline(0, settings["offset_m"], np.deg2rad(30))
    np.testing.assert_allclose(reference.mean(axis=0), settings["offset_m"])
    # Centre rotates by the readback difference (90 deg), not by mounting offset.
    moved = shape.outline(np.deg2rad(110 - 20), settings["offset_m"], np.deg2rad(30))
    np.testing.assert_allclose(moved.mean(axis=0), [-0.001, 0.002], atol=1e-18)
    settings["normal_rotation_confirmed"] = False  # Obsolete UI acknowledgement.
    vertical = bp.gaussian_profile(160e-6)
    before = shape_frame_factors(settings, vertical, 0.01, np.deg2rad(50))[0]
    settings["normal_rotation_confirmed"] = True
    np.testing.assert_array_equal(
        before, shape_frame_factors(settings, vertical, 0.01, np.deg2rad(50))[0]
    )


def test_editor_auto_defaults_overrides_and_cancel(qapp):
    """Angle overrides are optional and Cancel restores both angle modes."""
    parent = IntegrationCorrectionsDialog()
    try:
        editor = parent.sampleEditor
        editor.enabled.setChecked(True)
        assert editor.settings()["azimuth_source"] == "auto"
        assert editor.settings()["incidence_source"] == "auto"
        assert "reference_incidence_deg" not in editor.settings()
        assert editor.override_panel.isHidden()
        assert editor.off_centre_panel.isHidden()
        assert not editor.parallel_reading.isHidden()
        parent.onOk()
        saved = parent.settings()
        parent.show()
        editor.override_button.click()
        editor.angle_modes["incidence"].setCurrentText("Counter")
        editor.angle_sources["incidence"].setText("alpha_readback")
        editor.angle_modes["azimuth"].setCurrentText("Fixed")
        editor.fixed_angles["azimuth"].setValue(17)
        assert not editor.override_panel.isHidden()
        assert editor.settings()["incidence_source"] == "alpha_readback"
        assert editor.settings()["fixed_azimuth_deg"] == 17
        parent.reject()
        assert parent.settings() == saved
        assert editor.override_panel.isHidden()
    finally:
        parent.close()


def test_automatic_source_frame_preview_matches_policy(qapp):
    """Loaded-frame diagnostics resolve automatic and overridden angles alike."""
    host = qt.QWidget()
    main = SimpleNamespace(
        fscan=SimpleScan(),
        getMuOm=lambda: (np.deg2rad([0.36, 0.5]), np.deg2rad([10, 20])),
    )
    host._mainWindow = lambda: main
    parent = IntegrationCorrectionsDialog(host)
    settings = _settings()
    settings.pop("reference_incidence_deg")
    settings["azimuth_source"] = "auto"
    try:
        parent.setSettings({"sample_interception": settings, "shape_values": [160]})
        editor = parent.sampleEditor
        editor.scan_preview.setChecked(True)
        editor.frame_index.setValue(1)
        editor._preview()
        assert "0.5 deg" in editor.status.text()
        assert "20 deg" in editor.status.text()
        editor.angle_modes["incidence"].setCurrentText("Fixed")
        editor.fixed_angles["incidence"].setValue(0.8)
        editor._preview()
        assert "0.8 deg" in editor.status.text()
    finally:
        parent.close()
        host.close()


@pytest.mark.parametrize("sign", [1, -1])
@pytest.mark.parametrize("kind", ["rectangle", "polygon"])
def test_centred_alignment_conversion_preserves_old_angles_and_overlap(
    qapp, sign, kind,
):
    """Unwrapped parallel readings preserve either sign and asymmetric shapes."""
    settings = _settings()
    settings.update(azimuth_sign=sign, reference_azimuth_deg=713.123456789123,
                    orientation_deg=-431.234567891234)
    if kind == "polygon":
        settings["shape"] = {
            "kind": kind,
            "vertices_m": [[-0.004, -0.003], [0.005, -0.003],
                           [0.005, 0.003], [-0.004, 0.002]],
        }
    else:
        settings["shape"]["dimensions_m"] = [0.01, 0.006]
    expected = settings["reference_azimuth_deg"] - sign * settings["orientation_deg"]
    assert parallel_omega_reading(settings) == expected
    assert orientation_for_parallel_reading(settings, expected) == pytest.approx(
        settings["orientation_deg"]
    )
    canonical = dict(settings, reference_azimuth_deg=expected, orientation_deg=0)
    # Include near alignment after cancellation of large unwrapped angles;
    # section-local quadrature must preserve the equivalent representation.
    omega = np.deg2rad(expected + np.array([0, 1e-6, 7, 30, 73]))
    alpha = np.deg2rad([0.36, 0.36, 0.36, 0.5, 0.7])
    vertical = bp.gaussian_profile(160e-6)
    old_factors = shape_frame_factors(settings, vertical, alpha, omega)
    centred_factors = shape_frame_factors(canonical, vertical, alpha, omega)
    for old, centred in zip(old_factors[:-1], centred_factors[:-1]):
        np.testing.assert_allclose(old, centred, rtol=1e-10, atol=1e-13)
    # Equivalent angles can yield different roundoff-sized corner pieces and
    # hence different conservative error estimates, especially at exact zero.
    for factors in (old_factors, centred_factors):
        assert np.all(factors[-1] <= np.abs(factors[0]) * 1e-9)
    parent = IntegrationCorrectionsDialog()
    try:
        parent.setSettings({"sample_interception": settings})
        editor = parent.sampleEditor
        assert editor.parallel_reading.value() == pytest.approx(expected, abs=1e-10)
        assert not editor.off_centre_button.isChecked()
        assert editor.off_centre_panel.isHidden()
        assert not editor.parallel_reading.isHidden()
        assert ("polygon's x axis" if kind == "polygon" else "long edge") in (
            editor.parallel_label.text()
        )
        assert "0°" in editor.alignment_note.text() and "90°" in (
            editor.alignment_note.text()
        )
        saved = editor.settings()
        assert saved["reference_azimuth_deg"] == settings["reference_azimuth_deg"]
        assert saved["orientation_deg"] == settings["orientation_deg"]
        restored_state = corrections_from_nxdict(corrections_to_nxdict(
            CorrectionState(sample_interception=saved)
        ))
        assert restored_state.sample_interception == saved
        for old, restored in zip(old_factors,
                                 shape_frame_factors(saved, vertical, alpha, omega)):
            np.testing.assert_array_equal(old, restored)
        parent.setSettings({"sample_interception": saved})
        assert editor.settings()["orientation_deg"] == settings["orientation_deg"]
    finally:
        parent.close()


@pytest.mark.parametrize("sign", [1, -1])
def test_centred_alignment_edits_and_sign_keep_placement_reference(qapp, sign):
    """Editing alignment changes orientation without redefining future offsets."""
    settings = _settings()
    settings.update(azimuth_sign=sign, reference_azimuth_deg=37, orientation_deg=23)
    parent = IntegrationCorrectionsDialog()
    try:
        parent.setSettings({"sample_interception": settings})
        editor = parent.sampleEditor
        editor.parallel_reading.setValue(12)
        saved = editor.settings()
        assert saved["reference_azimuth_deg"] == 37
        assert saved["orientation_deg"] == sign * 25
        for reading, expected_edge_angle in [(12, 0), (12 + sign * 90, 90)]:
            edge_angle = sign * (reading - saved["reference_azimuth_deg"])
            edge_angle += saved["orientation_deg"]
            assert edge_angle == expected_edge_angle
        # Reversing direction retains the centred field's physical anchor.
        editor.sign.setCurrentText(f"{-sign:+d}")
        assert editor.parallel_reading.value() == 12
        assert editor.settings()["orientation_deg"] == -sign * 25
        editor.off_centre_button.click()
        editor.values["offset_x"].setValue(2)
        off_centre = editor.settings()
        assert editor.parallel_reading.isHidden()
        assert off_centre["reference_azimuth_deg"] == 37
        assert off_centre["orientation_deg"] == -sign * 25
        assert off_centre["offset_m"] == [0.002, 0]
        # Off-centre direction changes preserve both placement parameters.
        editor.sign.setCurrentText(f"{sign:+d}")
        assert editor.settings()["orientation_deg"] == off_centre["orientation_deg"]
    finally:
        parent.close()


@pytest.mark.parametrize("sign", [1, -1])
@pytest.mark.parametrize("kind", ["rectangle", "polygon", "circle"])
def test_off_centre_collapse_keeps_saved_placement_and_results(qapp, sign, kind):
    """Section visibility cannot canonicalize the offset orbit or the mounting."""
    settings = _settings()
    settings.update(azimuth_sign=sign, reference_azimuth_deg=21.234567891234,
                    orientation_deg=31.345678912345,
                    reference_incidence_deg=0.364567891234567,
                    offset_m=[0.001234567891234, -0.000987654321098])
    if kind == "circle":
        settings["shape"] = {"kind": kind, "dimensions_m": [0.01]}
    elif kind == "polygon":
        settings["shape"] = {
            "kind": kind, "vertices_m": [[-0.004, -0.003], [0.005, -0.003],
                                         [0.005, 0.003], [-0.004, 0.002]],
        }
    parent = IntegrationCorrectionsDialog()
    try:
        parent.setSettings({"sample_interception": settings})
        editor = parent.sampleEditor
        assert editor.off_centre_button.isChecked()
        assert not editor.off_centre_panel.isHidden()
        assert editor.parallel_reading.isHidden()
        before = editor.settings()
        for key in ("reference_azimuth_deg", "orientation_deg",
                    "reference_incidence_deg", "offset_m"):
            assert before[key] == settings[key]
        alpha, omega = np.deg2rad([0.36, 0.5]), np.deg2rad([21, 61])
        vertical = bp.gaussian_profile(160e-6)
        factors = shape_frame_factors(settings, vertical, alpha, omega)
        for old, loaded in zip(factors, shape_frame_factors(
            before, vertical, alpha, omega
        )):
            np.testing.assert_array_equal(old, loaded)
        editor.off_centre_button.click()
        assert editor.off_centre_panel.isHidden()
        assert editor.parallel_reading.isHidden()
        assert editor.settings() == before
        for old, collapsed in zip(factors, shape_frame_factors(
            editor.settings(), vertical, alpha, omega
        )):
            np.testing.assert_array_equal(old, collapsed)
        editor.off_centre_button.click()
        assert editor.settings() == before
        parent.setSettings({"sample_interception": before})
        assert editor.off_centre_button.isChecked()
        assert editor.settings() == before
    finally:
        parent.close()


def test_centred_circle_hides_orientation_and_restores_other_shapes(qapp):
    """Circle symmetry hides alignment without discarding the stored mounting."""
    settings = _settings()
    settings.update(reference_azimuth_deg=12, orientation_deg=5)
    parent = IntegrationCorrectionsDialog()
    try:
        parent.setSettings({"sample_interception": settings})
        editor = parent.sampleEditor
        editor.kind.setCurrentText("circle")
        assert editor.parallel_reading.isHidden()
        assert "needs no angular alignment" in editor.alignment_note.text()
        editor.off_centre_button.click()
        assert editor.values["orientation_deg"].isHidden()
        assert editor.settings()["orientation_deg"] == 5
        editor.kind.setCurrentText("rectangle")
        assert not editor.parallel_reading.isHidden()
        assert editor.parallel_reading.value() == 7
    finally:
        parent.close()


@pytest.mark.parametrize("sign", [1, -1])
def test_alignment_uses_actual_long_edge_when_rectangle_width_is_larger(qapp, sign):
    """A wide rectangle's long edge is local y, leaving old x placement intact."""
    settings = _settings()
    settings["shape"]["dimensions_m"] = [0.006, 0.01]
    settings.update(azimuth_sign=sign, reference_azimuth_deg=20, orientation_deg=30)
    expected = 20 - sign * (30 + 90)
    assert parallel_omega_reading(settings) == expected
    assert orientation_for_parallel_reading(settings, expected) == 30
    parent = IntegrationCorrectionsDialog()
    try:
        parent.setSettings({"sample_interception": settings})
        editor = parent.sampleEditor
        assert editor.parallel_reading.value() == expected
        assert editor.settings()["orientation_deg"] == 30
        editor.parallel_reading.setValue(15)
        assert editor.settings()["orientation_deg"] == sign * (20 - 15) - 90
        editor.values["length"].setValue(12)
        # Changing which dimension is longer updates only the displayed reading.
        assert editor.parallel_reading.value() == 15 + sign * 90
        assert editor.settings()["orientation_deg"] == sign * 5 - 90
    finally:
        parent.close()


@pytest.mark.parametrize("axis", ["mu", "th"])
def test_default_motor_labels_refresh_without_changing_saved_sources(qapp, axis):
    """Scan-dependent motor labels must not change mode selection or persistence."""
    host = qt.QWidget()
    scan = SimpleScan()
    scan.axisname = axis
    main = SimpleNamespace(fscan=scan, getMuOm=lambda: (
        np.deg2rad([0.4, 0.5]) if scan.axisname == "mu" else np.deg2rad(0.4),
        np.deg2rad(-30) if scan.axisname == "mu" else np.deg2rad([-30, -31]),
    ))
    host._mainWindow = lambda: main
    parent = IntegrationCorrectionsDialog(host)
    try:
        model = _settings()
        model.update(incidence_source="auto", azimuth_source="auto")
        parent.setSettings({"sample_interception": model})
        editor = parent.sampleEditor
        before = editor.settings()
        for scan_axis in (axis, "th" if axis == "mu" else "mu"):
            scan.axisname = scan_axis
            editor.refreshSourceFrames()
            incidence = ("mu (scan axis)" if scan_axis == "mu" else
                         "mu (configuration)")
            azimuth = ("omega = -th (scan axis)" if scan_axis == "th" else
                       "omega = -th (scan readback)")
            assert editor.angle_modes["incidence"].currentText() == incidence
            assert editor.angle_modes["azimuth"].currentText() == azimuth
            assert incidence in editor.angle_summary.text()
            assert azimuth in editor.angle_summary.text()
            assert editor.settings() == before
            # Restoring auto uses its stable ID, not the scan-dependent label.
            editor.setSettings(before)
            assert editor.angle_modes["incidence"].currentText() == incidence
            assert editor.angle_modes["azimuth"].currentText() == azimuth
        editor.angle_modes["incidence"].setCurrentText("Counter")
        editor.angle_sources["incidence"].setText("alpha_counter")
        editor.angle_modes["azimuth"].setCurrentText("Fixed")
        editor.fixed_angles["azimuth"].setValue(17)
        overridden = editor.settings()
        scan.axisname = axis
        editor.refreshSourceFrames()
        assert editor.angle_modes["incidence"].currentText() == "Counter"
        assert editor.angle_modes["azimuth"].currentText() == "Fixed"
        assert editor.settings() == overridden
    finally:
        parent.close()
        host.close()

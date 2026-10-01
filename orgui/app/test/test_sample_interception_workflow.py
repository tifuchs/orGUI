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
)
from orgui.backend.scans import sample_azimuth
from orgui.datautils.xrayutils.corrections import beamprofile as bp


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
    from orgui.app.sample_interception_dialog import SampleInterceptionDialog

    app = qt.QApplication.instance() or qt.QApplication([])
    parent = IntegrationCorrectionsDialog()
    dialog = SampleInterceptionDialog(parent, _settings())
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
    from orgui.app.sample_interception_dialog import SampleInterceptionDialog

    app = qt.QApplication.instance() or qt.QApplication([])
    host = qt.QWidget()
    host._mainWindow = lambda: SimpleNamespace(
        fscan=SimpleScan(),
        getMuOm=lambda: (np.deg2rad([0.36, 0.5]), 0),
    )
    parent = IntegrationCorrectionsDialog(host)
    parent.setSettings({"analytical": True, "shape": "Gaussian", "shape_values": [160]})
    dialog = SampleInterceptionDialog(parent, _settings())
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

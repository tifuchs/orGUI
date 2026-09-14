"""Round-trip tests for the typed NeXus layout of the stored settings.

The integration corrections used to be persisted as one opaque JSON string,
and the region-of-interest settings -- the advanced options, the ``s``
sampling and the automatic-sizing switches -- were not persisted at all. A
saved dataset was therefore not reproducible: reloading its configuration
restored three of the seven correction switches and none of the region
settings.

These tests pin the replacement: every field of :class:`CorrectionState` and
every region-of-interest setting survives
``state -> NeXus dict -> HDF5 file -> NeXus dict -> state``, and a
configuration written before the change still loads.
"""

import dataclasses
import json

import numpy as np
import pytest
from silx.io.dictdump import dicttonx, nxtodict

from orgui.app.config_data import (
    CORRECTIONS_SCHEMA_VERSION,
    CorrectionState,
    ROIState,
    corrections_from_nxdict,
    corrections_to_nxdict,
    roi_from_nxdict,
    roi_to_nxdict,
)


def _through_file(nxdict, tmp_path, name="group"):
    """Write one group to HDF5 and read it back, as a database would."""
    path = tmp_path / f"{name}.h5"
    dicttonx({name: nxdict}, path)
    return nxtodict(path)[name]


def _populated_corrections():
    """A correction state with every field away from its default."""
    return CorrectionState(
        use_mask=True,
        use_background=True,
        use_solid_angle=True,
        use_polarization=True,
        use_lorentz=True,
        use_footprint=True,
        use_normalization=False,
        repair_masked_pixels=True,
        repair_max_component_pixels=4,
        repair_max_span=3,
        repair_radius=2,
        repair_min_valid_neighbors=6,
        repair_use_pyfai_gaps=False,
        repair_gap_size_px=6,
        normalize_exposure=False,
        monitor_corrections=("mondio", "ic1"),
        excluded_frames=(3, 7, 11),
        mask_asset="assets/mask",
        background_asset="assets/bg",
        background_variance_asset="assets/bgvar",
        uncertainty_provenance={"background": "measured"},
        sample_length_m=5e-3,
        sample_width_m=3e-3,
        beam_flux_density=1.2e16,
        beam_shape_analytical=False,
        beam_shape_name="Trapezoid",
        beam_shape_values=(90.0, 30.0),
        beam_profile_file="profiles/beam.dat",
        beam_profile_content="height scan (-dI/dz)",
        beam_profile_unit="mm",
        beam_profile_center="median",
        beam_profile_offset_um=-12.5,
    )


def test_every_correction_field_is_covered_by_the_layout():
    """Adding a field to the dataclass must not silently stop being saved.

    This is the guard that makes the round-trip test below meaningful: it
    fails when a new field is added to :class:`CorrectionState` and left out
    of the test fixture, rather than the field quietly never being written.
    """
    populated = _populated_corrections()
    default = CorrectionState()
    unexercised = [
        entry.name
        for entry in dataclasses.fields(CorrectionState)
        if getattr(populated, entry.name) == getattr(default, entry.name)
    ]
    assert not unexercised, (
        f"these CorrectionState fields are still at their default in the "
        f"round-trip fixture, so the test would not notice them being "
        f"dropped: {unexercised}"
    )


def test_a_fully_populated_correction_state_round_trips(tmp_path):
    """Every field survives state -> NeXus -> file -> NeXus -> state."""
    state = _populated_corrections()

    loaded = corrections_from_nxdict(
        _through_file(corrections_to_nxdict(state), tmp_path, "corrections")
    )

    assert loaded == state


def test_a_default_correction_state_round_trips(tmp_path):
    """The all-defaults case, where most datasets are absent entirely.

    ``None`` and empty sequences are encoded by omission, so this is the
    path where the reader has to supply the defaults itself.
    """
    state = CorrectionState()
    nxdict = corrections_to_nxdict(state)

    # The optional values must genuinely not be in the file.
    assert "excluded_frames" not in nxdict
    assert "monitor_corrections" not in nxdict["normalization"]
    assert "max_span" not in nxdict["pixel_repair"]
    assert "use_lorentz" not in nxdict["switches"]
    assert nxdict["assets"] == {"@NX_class": "NXcollection"}

    loaded = corrections_from_nxdict(
        _through_file(nxdict, tmp_path, "corrections")
    )
    assert loaded == state


def test_the_unrecorded_switches_stay_none(tmp_path):
    """``None`` must not collapse to ``False`` on the way through a file.

    It is the difference between "this configuration says the Lorentz
    correction was off" and "this configuration predates the switch being
    stored", and only the first may change a loaded GUI.
    """
    loaded = corrections_from_nxdict(
        _through_file(corrections_to_nxdict(CorrectionState()), tmp_path)
    )

    assert loaded.use_lorentz is None
    assert loaded.use_footprint is None
    assert loaded.use_normalization is None


def test_a_nested_uncertainty_provenance_is_refused():
    """Rather than silently reintroducing a serialized blob."""
    state = CorrectionState(uncertainty_provenance={"a": {"b": 1}})

    with pytest.raises(ValueError, match="flat mapping"):
        corrections_to_nxdict(state)


def test_unknown_datasets_are_ignored(tmp_path):
    """A configuration from a newer orGUI still loads."""
    nxdict = corrections_to_nxdict(CorrectionState(use_mask=True))
    nxdict["switches"]["use_something_new"] = True
    nxdict["a_whole_new_group"] = {"@NX_class": "NXcollection", "value": 1}

    loaded = corrections_from_nxdict(_through_file(nxdict, tmp_path))

    assert loaded.use_mask is True


def test_the_layout_is_browsable_and_versioned(tmp_path):
    """The point of replacing the JSON string: values are real datasets."""
    nxdict = _through_file(
        corrections_to_nxdict(_populated_corrections()), tmp_path
    )

    assert nxdict["@orgui_schema_version"] == CORRECTIONS_SCHEMA_VERSION
    assert "json" not in nxdict
    assert bool(nxdict["switches"]["use_solid_angle"]) is True
    assert int(nxdict["pixel_repair"]["radius"]) == 2
    assert np.array_equal(np.asarray(nxdict["excluded_frames"]), [3, 7, 11])


def test_a_legacy_json_configuration_still_reads():
    """Existing databases must keep loading; only writing changed."""
    legacy = CorrectionState(
        use_mask=True, use_solid_angle=True, monitor_corrections=("mondio",)
    )
    group = {
        "@NX_class": "NXcollection",
        "json": json.dumps(legacy.to_dict(), sort_keys=True),
    }

    # The reader branch ConfigData.from_nxdict takes for such a file.
    loaded = CorrectionState.from_dict(json.loads(group["json"]))

    assert loaded.use_mask is True
    assert loaded.use_solid_angle is True
    assert loaded.monitor_corrections == ("mondio",)
    assert loaded.use_lorentz is None


def test_the_region_of_interest_settings_round_trip(tmp_path):
    """Advanced options, region sizes and s sampling all survive."""
    state = ROIState(
        region={
            "hsize": 20.0,
            "vsize": 6.0,
            "left": 20.0,
            "right": 20.0,
            "top": 5.0,
            "bottom": 5.0,
            "auto_hsize": True,
            "auto_vsize": False,
        },
        advanced={
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
        },
        rocking_scan={"delta_s": 0.0019001086, "max_s": 5.0},
    )

    loaded = roi_from_nxdict(_through_file(roi_to_nxdict(state), tmp_path, "roi"))

    assert loaded.region == pytest.approx(state.region)
    assert loaded.advanced == pytest.approx(state.advanced)
    assert loaded.rocking_scan == pytest.approx(state.rocking_scan)
    # delta_s is the value the auto-pinning settled on, to full precision.
    assert loaded.rocking_scan["delta_s"] == pytest.approx(
        0.0019001086, rel=0, abs=1e-12
    )


def test_the_region_of_interest_units_are_recorded(tmp_path):
    """A reader should not have to guess metres from pixels."""
    state = ROIState(
        region={"hsize": 20.0},
        advanced={"sample_size_x": 5e-4, "offset_x": 1.0},
        rocking_scan={"delta_s": 0.002},
    )

    nxdict = roi_to_nxdict(state)

    assert nxdict["region"]["@unit"] == "px"
    assert nxdict["advanced"]["@sample_size_unit"] == "m"
    assert nxdict["advanced"]["@offset_unit"] == "px"
    assert nxdict["rocking_scan"]["@unit"] == "rlu"


def test_an_absent_region_group_is_empty_not_fatal():
    """Every configuration written before this change has no such group."""
    state = roi_from_nxdict(None)

    assert state.is_empty()
    assert state.region == {}
    assert roi_to_nxdict(state) == {
        "@NX_class": "NXcollection",
        "@orgui_schema_version": 1,
    }


def _fake_gui(options, captured):
    """A stand-in exposing only what ConfigData reads and writes."""
    from types import SimpleNamespace

    from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg

    cell = CTRcalc.UnitCell([3.0, 3.0, 5.0], [90.0, 90.0, 90.0], name="bulk")
    ub = HKLVlieg.UBCalculator(cell, 15.0)
    ub.defaultU_GID()
    return SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=DetectorCalibration.Detector2D_SXRD(),
            crystal=cell,
            ubCal=ub,
            mu=0.0,
            chi=0.0,
            phi=0.0,
            n=1.0,
        ),
        scanSelector=SimpleNamespace(
            get_integration_options=lambda: dict(options),
            set_integration_options=captured.update,
        ),
    )


def test_the_gui_switches_and_region_settings_survive_from_gui_to_apply():
    """Closes the gap: reloading restored 3 of 7 switches and no region.

    ``from_gui`` must collect every switch and the region settings, and
    ``apply_to_gui`` must put them all back.
    """
    from orgui.app.config_data import ConfigData

    options = {
        "mask": True,
        "solid_angle": True,
        "polarization": True,
        "lorentz": True,
        "footprint": True,
        "normalization": True,
        "region": {"hsize": 21.0, "auto_vsize": True},
        "advanced": {"detector_inclination": True, "sample_size_x": 5e-4},
        "rocking_scan": {"delta_s": 0.0019, "max_s": 5.0},
    }
    captured = {}
    config = ConfigData.from_gui(_fake_gui(options, captured))

    assert config.corrections.use_lorentz is True
    assert config.corrections.use_footprint is True
    assert config.corrections.use_normalization is True
    assert config.roi.region["hsize"] == 21.0
    assert config.roi.rocking_scan["delta_s"] == pytest.approx(0.0019)

    config.apply_to_gui(_fake_gui(options, captured))

    for name in ("mask", "solid_angle", "polarization",
                 "lorentz", "footprint", "normalization"):
        assert captured[name] is True, name
    assert captured["region"]["hsize"] == 21.0
    assert captured["advanced"]["sample_size_x"] == pytest.approx(5e-4)
    assert captured["rocking_scan"]["max_s"] == 5.0


def test_an_old_configuration_does_not_touch_the_unrecorded_switches():
    """A file predating the switches must leave the GUI as the user has it."""
    from orgui.app.config_data import ConfigData

    captured = {}
    config = ConfigData.from_gui(_fake_gui({"mask": True}, captured))
    config.corrections.use_lorentz = None
    config.corrections.use_footprint = None
    config.corrections.use_normalization = None
    config.roi = ROIState()

    captured.clear()
    config.apply_to_gui(_fake_gui({}, captured))

    assert "lorentz" not in captured
    assert "footprint" not in captured
    assert "normalization" not in captured
    assert "region" not in captured

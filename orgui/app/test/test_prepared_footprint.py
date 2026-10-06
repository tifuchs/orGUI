"""Prepared footprint inputs remain reproducible after reopening and relocation."""

import copy
import json
import shutil
from hashlib import sha256
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app.config_data import ConfigData, CorrectionState
from orgui.app.sample_interception_config import (
    profile_from_settings, vertical_settings,
)
from orgui.backend.scans import Scan, SimulationScan, h5_Image, sample_azimuth
from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg
import orgui.reconstruction_job as jobs


class _GuardedScan(SimulationScan):
    @property
    def derived_azimuth(self):
        """Reject readback unless a derivation was supplied on this live scan."""
        if not hasattr(self, "_derived"):
            raise ValueError("Derived angle has not been calculated")
        return self._derived


class _SlicedScan(Scan):
    def __init__(self, source):
        self.filename = source
        self.axisname = "th"
        self.axis = np.arange(10, 16, dtype=float)
        self.th = self.axis
        self.omega = -self.th
        self.mu = 0.5
        self.offsetindex = 0
        self.scandatapoints = 6

    def __len__(self):
        return len(self.axis)

    def get_raw_img(self, index):
        """Return constant synthetic detector counts."""
        return h5_Image(np.full((2, 2), 10.0))

    @classmethod
    def parse_h5_node(cls, node):
        """Return empty constructor metadata for this synthetic backend."""
        return {}

    def slice(self, start, stop):
        """Select absolute frame indices using the backend's slice convention."""
        selected = type(self)(self.filename)
        selected.axis = selected.th = self.axis[start:stop]
        selected.omega = -selected.th
        selected.offsetindex = start
        return selected


def _config():
    cell = CTRcalc.UnitCell([3, 3, 3], [90, 90, 90])
    ub = HKLVlieg.UBCalculator(cell, 70)
    ub.defaultU()
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(729, 0.5, 0.5, pixelX=172, pixelY=172)
    detector.set_wavelength(ub.getLambda() * 1e-10)
    detector.detector.shape = detector.detector.max_shape = (2, 2)
    config = ConfigData(detector, cell, ub)
    config.corrections = CorrectionState(
        shared_frame_normalization=True, use_normalization=False,
        shared_frame_illumination=True, use_footprint=True,
        total_flux_calibrated=False, horizontal_interception="full",
        sample_length_m=0.01, beam_shape_analytical=True,
        beam_shape_name="Top hat", beam_shape_values=(160,),
    )
    return config


def _prepare(monkeypatch, root, config, scan):
    root.mkdir(parents=True, exist_ok=True)
    gui = SimpleNamespace(
        fscan=scan, ubcalc=SimpleNamespace(detectorCal=config.detector),
        database=SimpleNamespace(compression=None), numberthreads=1, maxMemory=64,
    )
    monkeypatch.setattr(jobs.ConfigData, "from_gui", lambda gui: config)
    # No mapping/probe is needed to exercise real snapshot/assets serialization.
    monkeypatch.setattr(jobs, "estimate_checkpoint_plan", lambda *args, **kwargs: {
        "per_grid": {"q_lab": {"files_per_job": 3}},
    })
    return jobs.prepare_job(
        gui, root / "job.json",
        grids=[dict(minimum=(-20, -20, -20), maximum=(20, 20, 20),
                    step=(20, 20, 20), frame="lab")],
        scratch_path=root / "scratch", output_path=root / "map.h5",
        compression_override="Raw", thread_override=1, threads_per_image=1,
    )


@pytest.mark.parametrize("unit", ["deg", "rad"])
@pytest.mark.parametrize("guarded", [False, True])
def test_preparation_freezes_custom_source_values_and_divisors(
    tmp_path, monkeypatch, unit, guarded,
):
    """Raw selected angles survive reopening, even when derivation is guarded."""
    config = _config()
    scan = (_GuardedScan if guarded else SimulationScan)(
        (2, 2), 10, 12, 3, fixed=0.5
    )
    values = np.array([10.01, 11.02, 12.03])
    if unit == "rad":
        values = np.deg2rad(values)
    if guarded:
        scan._derived = values
    else:
        scan.derived_azimuth = values
    config.corrections.sample_interception = {
        "version": 1, "enabled": True,
        "shape": {"kind": "rectangle", "dimensions_m": [0.01, 0.005]},
        "offset_m": [0, 0], "reference_incidence_deg": 0.5,
        "azimuth_source": "derived_azimuth", "azimuth_unit": unit,
        "azimuth_derivation": "calibrated theta plus frame-dependent offset",
        "horizontal": {"analytical": True, "shape": "Top hat",
                       "shape_values": [20000]},
    }
    before = jobs._reconstruction_frame_policy(config, scan)
    prepared = _prepare(monkeypatch, tmp_path, config, scan)
    saved = jobs.read_job(tmp_path / "job.json")
    jobs.verify_job(saved)
    reopened = saved.scan
    after = jobs._reconstruction_frame_policy(saved.config_data, reopened)
    snapshot = saved.config_data.corrections.sample_interception["azimuth_snapshot"]
    assert snapshot["source"] == "derived_azimuth"
    assert snapshot["unit"] == unit
    assert snapshot["frame_indices"] == [0, 1, 2]
    assert snapshot["frame_count"] == 3
    assert snapshot["derivation"].startswith("calibrated theta")
    np.testing.assert_array_equal(snapshot["values"], values)
    np.testing.assert_array_equal(
        before.interception_provenance["sample_azimuth_rad"],
        after.interception_provenance["sample_azimuth_rad"],
    )
    np.testing.assert_array_equal(
        before.illumination_divisor, after.illumination_divisor,
    )
    # A changed setting changes checkpoint identity, even on the same grid.
    changed = copy.deepcopy(saved)
    changed.config["orgui"]["integration_corrections"]["footprint"][
        "shared_frame_illumination"
    ] = False
    assert changed.digest != prepared.digest
    # Re-preparation cannot reuse a prior job's angles for another loaded scan.
    with pytest.raises(ValueError, match="prepare a new job|has not been calculated"):
        _prepare(monkeypatch, tmp_path / "new", saved.config_data, reopened)
    assert not (tmp_path / "new" / "scratch").exists()
    if guarded and unit == "deg":
        pytest.importorskip("orgui.datautils.xrayutils._reciprocal_reconstruction_cpp")
        result = jobs.run_job(tmp_path / "job.json")
        assert result["status"] == "complete"
        applied = jobs.read_job(tmp_path / "job.json").correction_provenance[
            "illumination"
        ]
        np.testing.assert_array_equal(applied["divisor"], before.illumination_divisor)
        np.testing.assert_array_equal(
            applied["sample_azimuth_rad"],
            before.interception_provenance["sample_azimuth_rad"],
        )


def test_snapshot_checks_slice_order_source_unit_and_reference():
    """A sliced source keeps its offset and rejects reordered or foreign angles."""
    scan = SimpleNamespace(offsetindex=4, _orgui_scan_reference_sha256="scan-id")
    settings = {
        "azimuth_source": "derived", "azimuth_unit": "deg",
        "azimuth_snapshot": {
            "version": 1, "source": "derived", "unit": "deg",
            "frame_count": 3, "frame_indices": [4, 5, 6],
            "values": [10, None, 12], "scan_reference_sha256": "scan-id",
        },
    }
    np.testing.assert_allclose(
        sample_azimuth(scan, settings, 3), np.deg2rad([10, np.nan, 12]), equal_nan=True,
    )
    for key, value in (
        ("source", "other"), ("unit", "rad"), ("frame_count", 2),
        ("frame_indices", [6, 5, 4]), ("scan_reference_sha256", "other-id"),
        ("values", [10]),
    ):
        wrong = copy.deepcopy(settings)
        wrong["azimuth_snapshot"][key] = value
        with pytest.raises(ValueError, match="Saved sample azimuth"):
            sample_azimuth(scan, wrong, 3)


def test_prepared_sliced_scan_reopens_matching_indices(tmp_path, monkeypatch):
    """A real ScanReference slice round trip keeps the matching derived values."""
    source = tmp_path / "scan.dat"
    source.touch()
    scan = _SlicedScan(str(source)).slice(2, 5)
    scan.derived = np.array([12.01, 13.02, 14.03])
    config = _config()
    config.corrections.sample_interception = {
        "version": 1, "enabled": True, "azimuth_source": "derived",
        "shape": {"kind": "rectangle", "dimensions_m": [0.01, 0.005]},
        "horizontal": {"analytical": True, "shape": "Top hat",
                       "shape_values": [20000]},
    }
    before = jobs._reconstruction_frame_policy(config, scan)
    _prepare(monkeypatch, tmp_path / "job", config, scan)
    saved = jobs.read_job(tmp_path / "job" / "job.json")
    reopened = saved.scan
    assert reopened.offsetindex == 2 and len(reopened) == 3
    shape = saved.config_data.corrections.sample_interception
    assert shape["azimuth_snapshot"]["frame_indices"] == [2, 3, 4]
    np.testing.assert_array_equal(sample_azimuth(reopened, shape, 3),
                                  np.deg2rad(scan.derived))
    after = jobs._reconstruction_frame_policy(saved.config_data, reopened)
    np.testing.assert_array_equal(
        before.illumination_divisor, after.illumination_divisor,
    )


@pytest.mark.parametrize("storage", [None, "file"])
def test_measured_profile_job_relocation(tmp_path, monkeypatch, storage):
    """Embedding is self-contained; explicit files resolve beside the job and hash."""
    root = tmp_path / "original"
    fixture = root / "fixture" / "profile.dat"
    fixture.parent.mkdir(parents=True)
    np.savetxt(fixture, [[-0.12, 1], [-0.04, 3], [0.03, 2], [0.1, 1]])
    config = _config()
    state = config.corrections
    state.beam_shape_analytical = False
    state.beam_profile_file = "fixture/profile.dat"
    state.beam_profile_base = "."
    state.beam_profile_storage = storage
    state.beam_profile_unit = "mm"
    state.beam_profile_center = "median"
    state.beam_profile_offset_um = 23
    scan = SimulationScan((2, 2), 10, 12, 3, fixed=0.5)
    prepared = _prepare(monkeypatch, root, config, scan)
    before_profile = profile_from_settings(vertical_settings(state), base_path=root)
    before = jobs._reconstruction_frame_policy(prepared.config_data, scan)
    copied = tmp_path / "copied"
    shutil.copytree(root, copied)
    fixture.unlink()
    saved = jobs.read_job(copied / "job.json")
    saved_state = saved.config_data.corrections
    assert saved_state.beam_profile_unit == "mm"
    assert saved_state.beam_profile_center == "median"
    assert saved_state.beam_profile_offset_um == 23
    after_profile = profile_from_settings(
        vertical_settings(saved_state), base_path=copied,
    )
    after = jobs._reconstruction_frame_policy(saved.config_data, saved.scan)
    np.testing.assert_allclose(
        before_profile.profile_curve(), after_profile.profile_curve(),
    )
    np.testing.assert_array_equal(
        before.illumination_divisor, after.illumination_divisor,
    )
    if storage == "file":
        assert not saved_state.beam_profile_positions_m
        assert saved_state.beam_profile_sha256 == sha256(
            (copied / "fixture" / "profile.dat").read_bytes()
        ).hexdigest()
        (copied / "fixture" / "profile.dat").write_text("-1 1\n1 1\n")
        with pytest.raises(ValueError, match="content changed"):
            jobs._reconstruction_frame_policy(saved.config_data, saved.scan)
    else:
        assert saved_state.beam_profile_positions_m
    (copied / "fixture" / "profile.dat").unlink()
    if storage == "file":
        with pytest.raises(ValueError, match="Cannot read measured beam profile"):
            jobs._reconstruction_frame_policy(saved.config_data, saved.scan)
    else:
        policy = jobs._reconstruction_frame_policy(saved.config_data, saved.scan)
        np.testing.assert_array_equal(
            before.illumination_divisor, policy.illumination_divisor,
        )


def test_legacy_status_and_recorded_status_preserve_saved_jobs(tmp_path, monkeypatch):
    """Compatibility and recorded application are visible without rewriting jobs."""
    config = _config()
    scan = SimulationScan((2, 2), 10, 12, 3, fixed=0.5)
    job = _prepare(monkeypatch, tmp_path, config, scan)
    path = tmp_path / "job.json"
    new = jobs.job_status(path)["illumination"]
    assert new["effective"] and new["resolved_status"] == "applied"
    assert new["applied_status"] == "not_recorded"
    assert new["convention"] == "total_flux_H"
    values = job.to_dict()
    del values["config"]["orgui"]["integration_corrections"]["footprint"][
        "shared_frame_illumination"
    ]
    path.write_text(json.dumps(values))
    original = path.read_bytes()
    legacy_digest = jobs.read_job(path).digest
    legacy = jobs.job_status(path)["illumination"]
    assert legacy["footprint_requested"] and not legacy["effective"]
    assert "compatibility" in legacy["reason"]
    assert path.read_bytes() == original
    assert jobs.read_job(path).digest == legacy_digest
    job.status = "complete"
    job.correction_provenance = {"illumination": {
        "status": "applied", "convention": "total_flux_H", "divisor": [2, 3, 4],
    }}
    jobs.write_job(job, path)
    monkeypatch.setattr(
        jobs.ScanReference, "open", lambda self: pytest.fail("reopened"),
    )
    recorded = jobs.job_status(path)["illumination"]
    assert recorded["applied_status"] == "applied"
    assert recorded["resolved_status"] == "recorded"
    assert recorded["divisor"] == [2, 3, 4]


def test_profile_reference_fields_round_trip_without_runtime_base():
    """JSON/NeXus retain explicit references without persisting a runtime root."""
    from orgui.app.config_data import corrections_from_nxdict, corrections_to_nxdict

    state = CorrectionState(
        beam_profile_base=".", beam_profile_sha256="abcd", beam_profile_storage="file",
    )
    state._profile_root = "/temporary/runtime/root"
    assert "_profile_root" not in state.to_dict()
    assert CorrectionState.from_dict(state.to_dict()) == state
    assert corrections_from_nxdict(corrections_to_nxdict(state)) == state

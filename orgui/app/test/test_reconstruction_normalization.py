"""Frame normalization shared by ROI integration and reconstruction."""

from types import SimpleNamespace
import json

import numpy as np
import pytest

from orgui.app.config_data import (
    CorrectionState, corrections_from_nxdict, corrections_to_nxdict,
)
from orgui.app.integration_corrections import frame_correction_policy
from orgui.reconstruction_job import (
    _correction_pipeline, _reconstruction_frame_policy, prepare_job,
)
import orgui.reconstruction_job as job_module


class _Scan:
    exposure_time = np.array([2.0])
    exposure_time_variance = np.array([0.04])
    mon = np.array([5.0])
    mon_variance = np.array([1.0])

    def __len__(self):
        return 1


@pytest.mark.parametrize(
    ("kind", "expected_divisor", "expected_variance"),
    [("rate", 10.0, 6.0), ("integrated", 5.0, 20.0)],
)
def test_reconstruction_uses_integration_primary_monitor_policy(
    kind, expected_divisor, expected_variance
):
    state = CorrectionState(
        use_normalization=True,
        shared_frame_normalization=True,
        total_flux_calibrated=False,
        primary_monitor="mon",
        primary_monitor_kind=kind,
        monitor_corrections=("unused_legacy_counter",),
    )
    scan = _Scan()
    policy = frame_correction_policy(
        scan, state, 1, use_normalization=True, use_illumination=False
    )
    provenance = {}
    pipeline = _correction_pipeline(
        SimpleNamespace(corrections=state, detector=object()),
        scan, {}, provenance,
    )

    intensity, variance, mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )

    assert policy.normalization_divisor[0] == expected_divisor
    assert intensity[0, 0] == pytest.approx(100.0 / expected_divisor)
    assert variance[0, 0] == pytest.approx(expected_variance)
    assert not mask[0, 0]
    assert provenance["factor_uncertainty"]["normalization"] == "propagated"


def test_reconstruction_respects_shared_normalization_switch():
    state = CorrectionState(
        use_normalization=False,
        shared_frame_normalization=True,
        total_flux_calibrated=False,
        primary_monitor="mon",
        primary_monitor_kind="rate",
    )
    pipeline = _correction_pipeline(
        SimpleNamespace(corrections=state, detector=object()),
        _Scan(), {}, {},
    )

    intensity, variance, _mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )

    assert intensity[0, 0] == 100.0
    assert variance[0, 0] == 100.0


def test_reconstruction_uses_calibrated_primary_monitor_reference():
    state = CorrectionState(
        use_normalization=True,
        shared_frame_normalization=True,
        total_flux_calibrated=True,
        total_incident_flux=100.0,
        primary_monitor="mon",
        primary_monitor_kind="rate",
        monitor_reference_reading=10.0,
    )
    pipeline = _correction_pipeline(
        SimpleNamespace(corrections=state, detector=object()),
        _Scan(), {}, {},
    )

    intensity, _variance, _mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )

    # 100 photons/s * 2 s * (5 / 10) = 100 incident photons.
    assert intensity[0, 0] == pytest.approx(1.0)


def test_legacy_job_keeps_exposure_times_monitor_product():
    state = CorrectionState(
        use_normalization=False,
        normalize_exposure=True,
        monitor_corrections=("mon",),
        total_flux_calibrated=False,
        primary_monitor="mon",
        primary_monitor_kind="integrated",
    )
    pipeline = _correction_pipeline(
        SimpleNamespace(corrections=state, detector=object()),
        _Scan(), {}, {},
    )

    intensity, _variance, _mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )

    assert intensity[0, 0] == pytest.approx(10.0)


def _footprint_state(**kwargs):
    values = dict(
        shared_frame_illumination=True, use_footprint=True,
        shared_frame_normalization=True, use_normalization=True,
        total_flux_calibrated=False, primary_monitor="mon",
        primary_monitor_kind="rate", horizontal_interception="full",
        sample_length_m=0.01, beam_shape_analytical=True,
        beam_shape_name="Top hat", beam_shape_values=(160.0,),
    )
    values.update(kwargs)
    return CorrectionState(**values)


@pytest.mark.parametrize("native", [False, True])
@pytest.mark.parametrize("kind,normalization,variance", [
    ("rate", 10.0, 6.0), ("integrated", 5.0, 20.0),
])
def test_map_divides_signal_and_monitor_variance_by_illumination(
    native, kind, normalization, variance, monkeypatch,
):
    """Top-hat interception has an independent analytic H = L / height."""
    if native:
        pytest.importorskip(
            "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
        )
    else:
        monkeypatch.setattr(job_module, "_correction_extension", lambda: None)
    state = _footprint_state(primary_monitor_kind=kind)
    provenance = {}
    config = SimpleNamespace(corrections=state, detector=object(), mu=0.006)
    pipeline = _correction_pipeline(config, _Scan(), {}, provenance)
    intensity, result_variance, mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )
    h = 0.01 / 160e-6
    assert intensity[0, 0] == pytest.approx(100.0 / normalization / h)
    assert result_variance[0, 0] == pytest.approx(variance / h**2)
    assert not mask[0, 0]
    assert provenance["illumination"]["divisor"] == pytest.approx([h])
    assert provenance["illumination"]["convention"] == "total_flux_H"
    assert provenance["factor_uncertainty"]["illumination"] == (
        "deterministic-no-uncertainty"
    )
    json.dumps(provenance)


@pytest.mark.parametrize("total_flux", [False, True])
def test_map_footprint_without_normalization_preserves_convention(total_flux):
    """Legacy area and total-flux H retain their different scales."""
    state = _footprint_state(use_normalization=False)
    if not total_flux:
        state.primary_monitor = state.primary_monitor_kind = None
        state.total_flux_calibrated = state.horizontal_interception = None
    scan = _Scan()
    scan.mu = 2.0  # Scan readbacks are degrees and override fixed config mu.
    config = SimpleNamespace(corrections=state, detector=object(), mu=0.006)
    provenance = {}
    pipeline = _correction_pipeline(config, scan, {}, provenance)
    intensity, variance, _mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )
    h = (1 / np.sin(np.deg2rad(2.0)) if total_flux
         else 160e-6 / (0.01 * np.sin(np.deg2rad(2.0))))
    assert intensity[0, 0] == pytest.approx(100.0 / h)
    assert variance[0, 0] == pytest.approx(100.0 / h**2)
    assert provenance["illumination"]["convention"] == (
        "total_flux_H" if total_flux else "legacy_C_illum_area"
    )


def test_map_uses_sample_shape_acquisition_angles_and_overrides():
    """Square/Gaussian interception matches an independent numerical reference."""
    shape = {
        "version": 1, "enabled": True,
        "shape": {"kind": "rectangle", "dimensions_m": [0.01, 0.01]},
        "offset_m": [0, 0], "azimuth_source": "auto",
        "horizontal": {
            "analytical": True, "shape": "Top hat", "shape_values": [20000],
        },
    }
    state = _footprint_state(
        use_normalization=False, beam_shape_name="Gaussian",
        sample_interception=shape,
    )

    class Scan:
        th = np.array([0.0, -45.0])
        mu = np.array([0.36, 0.36])

        def __len__(self):
            return 2

    config = SimpleNamespace(corrections=state, detector=object(), mu=0.0)
    provenance = {}
    pipeline = _correction_pipeline(config, Scan(), {}, provenance)
    h = np.array([0.178090128539, 0.178155252900]) / np.sin(np.deg2rad(0.36))
    for frame in range(2):
        intensity, variance, _mask = pipeline.correct_frame(
            SimpleNamespace(), np.array([[100.0]]), frame
        )
        assert intensity[0, 0] == pytest.approx(100.0 / h[frame], rel=1e-10)
        assert variance[0, 0] == pytest.approx(100.0 / h[frame]**2, rel=1e-10)
    saved = provenance["illumination"]
    np.testing.assert_allclose(saved["sample_azimuth_rad"], np.deg2rad([0, 45]))
    assert saved["interception_method"] == "surface_chord_quadrature_v1"
    json.dumps(provenance)
    shape.update(incidence_source="fixed", fixed_incidence_deg=2.0,
                 azimuth_source="fixed", fixed_azimuth_deg=10.0)
    policy = _reconstruction_frame_policy(config, Scan())
    np.testing.assert_allclose(policy.interception_provenance["sample_incidence_rad"],
                               np.deg2rad([2.0, 2.0]))
    np.testing.assert_allclose(policy.interception_provenance["sample_azimuth_rad"],
                               np.deg2rad([10.0, 10.0]))


@pytest.mark.parametrize("native", [False, True])
def test_map_masks_frames_with_undefined_illumination(native, monkeypatch):
    """Invalid incidence contributes no pixels through either scaling path."""
    if native:
        pytest.importorskip(
            "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
        )
    else:
        monkeypatch.setattr(job_module, "_correction_extension", lambda: None)

    class Scan:
        mu = np.array([0.0, 0.36, -1.0])

        def __len__(self):
            return 3

    config = SimpleNamespace(
        corrections=_footprint_state(use_normalization=False), detector=object(),
    )
    pipeline = _correction_pipeline(config, Scan(), {}, {})
    for frame in range(3):
        intensity, variance, mask = pipeline.correct_frame(
            SimpleNamespace(), np.ones((2, 2)), frame
        )
        assert np.all(mask) == (frame != 1)
        assert np.all(np.isfinite(intensity)) == (frame == 1)
        assert np.all(np.isfinite(variance)) == (frame == 1)


@pytest.mark.parametrize("old_job,enabled", [(True, True), (False, False)])
def test_old_jobs_and_disabled_footprint_skip_unconfigured_inputs(old_job, enabled):
    """Opt-in and checkbox are both required; unused profiles are never read."""
    state = _footprint_state(
        shared_frame_illumination=not old_job, use_footprint=enabled,
        beam_shape_name=None,
    )
    provenance = {}
    pipeline = _correction_pipeline(
        SimpleNamespace(corrections=state, detector=object()), _Scan(), {}, provenance,
    )
    intensity, variance, _mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )
    assert intensity[0, 0] == pytest.approx(10.0)
    assert variance[0, 0] == pytest.approx(6.0)
    assert "illumination" not in provenance


def test_prepare_rejects_missing_footprint_before_creating_assets(
    monkeypatch, tmp_path,
):
    """Preparation validates requested illumination without starting a run."""
    config = SimpleNamespace(
        corrections=_footprint_state(beam_shape_name=None), mu=0.006,
    )
    monkeypatch.setattr(job_module.ConfigData, "from_gui", lambda gui: config)
    scratch = tmp_path / "scratch"
    gui = SimpleNamespace(fscan=_Scan())
    with pytest.raises(ValueError, match="unknown analytical beam shape"):
        prepare_job(
            gui, tmp_path / "job.json",
            grids=[dict(minimum=(0, 0, 0), maximum=(1, 1, 1), step=(1, 1, 1),
                        frame="lab")],
            scratch_path=scratch, output_path=tmp_path / "map.h5",
        )
    assert not scratch.exists()


def test_illumination_opt_in_round_trips_and_defaults_off_for_old_jobs():
    """Saved correction snapshots explicitly distinguish corrected map jobs."""
    state = _footprint_state()
    assert CorrectionState.from_dict(state.to_dict()) == state
    assert corrections_from_nxdict(corrections_to_nxdict(state)) == state
    old = state.to_dict()
    del old["shared_frame_illumination"]
    assert not CorrectionState.from_dict(old).shared_frame_illumination
    saved = corrections_to_nxdict(state)
    del saved["footprint"]["shared_frame_illumination"]
    assert not corrections_from_nxdict(saved).shared_frame_illumination


def test_map_uses_embedded_measured_profile_without_original_file():
    """Embedded, centred profile coordinates remain independent of source paths."""
    state = _footprint_state(
        use_normalization=False, beam_shape_analytical=False,
        beam_profile_file="missing-profile.dat", beam_profile_offset_um=123.0,
        beam_profile_positions_m=(-80e-6, 80e-6),
        beam_profile_density_per_m=(6250.0, 6250.0),
    )
    config = SimpleNamespace(corrections=state, detector=object(), mu=0.006)
    pipeline = _correction_pipeline(config, _Scan(), {}, {})
    intensity, variance, _mask = pipeline.correct_frame(
        SimpleNamespace(), np.array([[100.0]]), 0
    )
    h = 0.01 / 160e-6
    assert intensity[0, 0] == pytest.approx(100 / h)
    assert variance[0, 0] == pytest.approx(100 / h**2)


def test_map_rejects_scan_without_any_physical_incidence():
    """A requested correction cannot silently become a unity factor."""
    config = SimpleNamespace(corrections=_footprint_state(), mu=0.0)
    with pytest.raises(ValueError, match="at least one frame"):
        _reconstruction_frame_policy(config, _Scan())

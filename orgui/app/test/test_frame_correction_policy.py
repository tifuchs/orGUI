"""Regression tests for the versioned framewise CTR correction policy."""

from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app import integration_corrections as ic
from orgui.app.config_data import CorrectionState, CurveCorrectionRecord
from orgui.datautils.xrayutils.corrections import measurement
from orgui.datautils.xrayutils.corrections.beamprofile import top_hat_profile


def _state(**updates):
    values = {
        "total_incident_flux": 2.0e10,
        "total_flux_calibrated": True,
        "primary_monitor": "ic2",
        "primary_monitor_kind": "rate",
        "primary_monitor_unit": "count/s",
        "monitor_reference_reading": 100.0,
        "horizontal_interception": "full",
    }
    values.update(updates)
    return CorrectionState(**values)


def test_calibrated_rate_monitor_counts_exposure_exactly_once():
    """A rate-like ion chamber is converted to frame photon fluence."""
    exposure = np.array([0.5, 1.0, 2.0])
    monitor = np.array([50.0, 100.0, 200.0])
    scan = SimpleNamespace(exposure_time=exposure, ic2=monitor)

    policy = ic.frame_correction_policy(
        scan,
        _state(),
        3,
        use_normalization=True,
        use_illumination=False,
    )

    expected = 2.0e10 * exposure * monitor / 100.0
    np.testing.assert_allclose(policy.normalization_divisor, expected)
    assert policy.normalization_unit == "photons"
    assert policy.normalization_components == (
        "exposure",
        "primary_monitor:ic2:rate",
    )


def test_integrated_monitor_does_not_multiply_frame_exposure_again():
    """An integrated monitor already contains the acquisition duration."""
    scan = SimpleNamespace(
        exposure_time=np.array([0.5, 1.0, 2.0]),
        ic2=np.array([50.0, 100.0, 200.0]),
    )
    state = _state(
        primary_monitor_kind="integrated",
        monitor_reference_exposure_s=0.25,
    )

    policy = ic.frame_correction_policy(
        scan,
        state,
        3,
        use_normalization=True,
        use_illumination=False,
    )

    expected = 2.0e10 * 0.25 * scan.ic2 / 100.0
    np.testing.assert_allclose(policy.normalization_divisor, expected)
    assert policy.normalization_components == (
        "primary_monitor:ic2:integrated",
    )


def test_total_flux_illumination_requires_horizontal_provenance():
    """The unmeasured beam direction cannot silently be called complete."""
    scan = SimpleNamespace(exposure_time=np.ones(3), ic2=np.ones(3) * 100.0)
    with pytest.raises(ValueError, match="horizontal_interception"):
        ic.frame_correction_policy(
            scan,
            _state(horizontal_interception=None),
            3,
            use_normalization=True,
            use_illumination=True,
            alpha=np.deg2rad([0.2, 0.5, 1.0]),
            beam_profile=top_hat_profile(300e-6),
            sample_length=5e-3,
        )


def test_total_flux_illumination_excludes_frames_at_the_horizon():
    """A scan starting at alpha = 0 drops that frame instead of failing.

    Regression: ``flyscan th 0 12`` on the CHESS QM2 backend (th -> alpha)
    rejected the whole extraction because of its first frame.
    """
    scan = SimpleNamespace(exposure_time=np.ones(4), ic2=np.ones(4) * 100.0)
    alpha = np.deg2rad([[0.0, 0.5, 1.0, -0.1]] * 2)
    profile = top_hat_profile(300e-6)

    policy = ic.frame_correction_policy(
        scan,
        _state(),
        4,
        use_normalization=True,
        use_illumination=True,
        alpha=alpha,
        beam_profile=profile,
        sample_length=5e-3,
    )

    valid = np.array([[False, True, True, False]] * 2)
    expected = ic.activearea.illumination_divisor(
        alpha[valid], 5e-3, profile, horizontal_fraction=1.0
    )
    for values in (
        policy.illumination_divisor,
        policy.intercepted_fraction,
        policy.vertical_intercepted_fraction,
    ):
        assert values.shape == alpha.shape
        np.testing.assert_array_equal(np.isfinite(values), valid)
    np.testing.assert_allclose(policy.illumination_divisor[valid], expected)
    assert policy.illumination_status == "applied"


def test_total_flux_illumination_needs_one_physical_frame():
    """Excluding frames must not turn an all-invalid scan into all-NaN data."""
    with pytest.raises(ValueError, match="at least one frame"):
        ic.framewise_illumination_divisor(
            np.deg2rad([0.0, -1.0]),
            5e-3,
            top_hat_profile(300e-6),
            horizontal_fraction=1.0,
        )


def test_legacy_policy_preserves_counter_product_and_active_area():
    """Absent version-3 fields retain the established numerical convention."""
    exposure = np.array([0.5, 1.0, 2.0])
    monitor = np.array([10.0, 20.0, 40.0])
    alpha = np.deg2rad([0.2, 0.5, 1.0])
    profile = top_hat_profile(300e-6)
    scan = SimpleNamespace(exposure_time=exposure, mon=monitor)
    state = CorrectionState(
        normalize_exposure=True, monitor_corrections=("mon",)
    )

    policy = ic.frame_correction_policy(
        scan,
        state,
        3,
        use_normalization=True,
        use_illumination=True,
        alpha=alpha,
        beam_profile=profile,
        sample_length=5e-3,
    )

    vertical, area = profile.corrections(alpha, 5e-3)
    np.testing.assert_allclose(policy.normalization_divisor, exposure * monitor)
    np.testing.assert_allclose(policy.vertical_intercepted_fraction, vertical)
    np.testing.assert_allclose(policy.illumination_divisor, area)
    assert policy.new_contract is False
    assert policy.scale_convention == "legacy_relative"


def _record(**updates):
    values = {
        "algorithm": "framewise_ctr_total_flux_v1",
        "output_quantity": "rocking_ctr_photon_curve",
        "scale_convention": "total_flux_calibrated",
        "normalization_status": "applied",
        "normalization_divisor": np.array([2.0, 4.0]),
        "illumination_status": "applied",
        "illumination_divisor": np.array([0.5, 0.25]),
        "illumination_convention": "total_flux_H",
        "base_croibg": np.array([80.0, 96.0]),
        "base_croibg_variance": np.array([16.0, 36.0]),
    }
    values.update(updates)
    return CurveCorrectionRecord(**values)


def test_footprint_actions_always_restart_from_the_stored_base():
    """Keep, remove, and replace never compound a previous transformation."""
    record = _record()

    kept, kept_errors, status, convention = ic.corrected_curve_from_record(record)
    removed, removed_errors, removed_status, _ = ic.corrected_curve_from_record(
        record, footprint_action=ic.FOOTPRINT_REMOVE
    )
    replaced, replaced_errors, replaced_status, _ = ic.corrected_curve_from_record(
        record,
        footprint_action=ic.FOOTPRINT_APPLY,
        replacement_illumination=np.array([0.25, 0.5]),
        replacement_convention="total_flux_H",
    )

    np.testing.assert_allclose(kept, [80.0, 96.0] / np.array([2.0, 4.0]) / [0.5, 0.25])
    np.testing.assert_allclose(removed, [80.0, 96.0] / np.array([2.0, 4.0]))
    np.testing.assert_allclose(
        replaced, [80.0, 96.0] / np.array([2.0, 4.0]) / [0.25, 0.5]
    )
    np.testing.assert_allclose(
        kept_errors, [4.0, 6.0] / np.array([2.0, 4.0]) / [0.5, 0.25]
    )
    np.testing.assert_allclose(removed_errors, [4.0, 6.0] / np.array([2.0, 4.0]))
    np.testing.assert_allclose(
        replaced_errors, [4.0, 6.0] / np.array([2.0, 4.0]) / [0.25, 0.5]
    )
    assert (status, convention) == ("applied", "total_flux_H")
    assert removed_status == "removed"
    assert replaced_status == "replaced"


def test_excluded_illumination_frame_is_nan_only_in_that_frame():
    """A NaN H marks an excluded frame; a nonpositive H is still an error."""
    curve, errors, _status, _convention = ic.corrected_curve_from_record(
        _record(illumination_divisor=np.array([np.nan, 0.25]))
    )
    assert np.isnan(curve[0]) and np.isnan(errors[0])
    np.testing.assert_allclose(curve[1], 96.0 / 4.0 / 0.25)
    np.testing.assert_allclose(errors[1], 6.0 / 4.0 / 0.25)

    for bad in ([0.0, 0.25], [np.inf, 0.25], [np.nan, np.nan]):
        with pytest.raises(ValueError, match="positive where defined"):
            ic.corrected_curve_from_record(
                _record(illumination_divisor=np.array(bad))
            )


def test_convention_change_and_unknown_provenance_are_explicit_errors():
    """Unsafe migrations are never inferred from a curve's numeric values."""
    with pytest.raises(ValueError, match="changing illumination convention"):
        ic.corrected_curve_from_record(
            _record(),
            footprint_action=ic.FOOTPRINT_APPLY,
            replacement_illumination=np.ones(2),
            replacement_convention="legacy_C_illum_area",
        )

    with pytest.raises(ValueError, match="normalization provenance"):
        ic.corrected_curve_from_record(_record(normalization_status="unknown"))


def test_calibrated_stationary_structure_factor_uses_total_flux_prefactor():
    """The total-flux scale omits the legacy flux-density/area product."""
    alpha = np.deg2rad(np.array([0.5]))
    delta = np.deg2rad(np.array([20.0]))
    gamma = np.deg2rad(np.array([3.0]))
    factors = ic.stationary_correction_factors(
        alpha, delta, gamma, use_lorentz=True
    )
    policy = ic.FrameCorrectionPolicy(
        new_contract=True,
        scale_convention="total_flux_calibrated",
        normalization_status="applied",
        illumination_status="applied",
        calibrated=True,
    )

    result, errors = ic.structure_factor_from_policy(
        np.array([3.0]),
        np.array([0.3]),
        factors,
        policy,
        wavelength=1.0,
        unitcell_area=16.0,
    )

    prefactor = measurement.total_flux_prefactor(1.0, 16.0)
    expected = 3.0 / factors["C_Lorentz"] / prefactor
    np.testing.assert_allclose(result, expected)
    np.testing.assert_allclose(errors / result, 0.1)

"""Frame normalization shared by ROI integration and reconstruction."""

from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app.config_data import CorrectionState
from orgui.app.integration_corrections import frame_correction_policy
from orgui.reconstruction_job import _correction_pipeline


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

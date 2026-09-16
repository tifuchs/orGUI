"""Tests for calibrated and relative per-frame incident fluence."""

import numpy as np
import pytest

from orgui.datautils.xrayutils.corrections import normalization


def test_constant_total_flux_uses_each_frame_exposure():
    exposure = np.array([0.1, 0.25, 0.8])

    fluence = normalization.frame_fluence(4.0e11, exposure_time=exposure)

    np.testing.assert_allclose(fluence, 4.0e11 * exposure)


def test_rate_and_integrated_monitors_give_the_same_photons():
    """An integrated reading contains exposure; a rate reading does not."""
    reference_flux = 5.0e11
    reference_rate = 2.0e6
    reference_exposure = 0.4
    exposure = np.array([0.1, 0.35, 1.2])
    flux = np.array([4.0e11, 6.0e11, 2.5e11])
    rate = reference_rate * flux / reference_flux
    integrated = rate * exposure

    from_rate = normalization.frame_fluence(
        reference_flux,
        exposure_time=exposure,
        monitor=rate,
        monitor_kind="rate",
        reference_monitor=reference_rate,
    )
    from_integrated = normalization.frame_fluence(
        reference_flux,
        exposure_time=exposure,
        monitor=integrated,
        monitor_kind="integrated",
        reference_monitor=reference_rate * reference_exposure,
        reference_exposure=reference_exposure,
    )

    expected = flux * exposure
    np.testing.assert_allclose(from_rate, expected, rtol=1e-14)
    np.testing.assert_allclose(from_integrated, expected, rtol=1e-14)


def test_integrated_monitor_does_not_apply_frame_exposure_twice():
    integrated = np.array([10.0, 20.0, 40.0])

    first = normalization.frame_fluence(
        1.0e10,
        exposure_time=np.array([0.1, 0.2, 0.4]),
        monitor=integrated,
        monitor_kind="integrated",
        reference_monitor=20.0,
        reference_exposure=0.5,
    )
    second = normalization.frame_fluence(
        1.0e10,
        exposure_time=np.array([7.0, 8.0, 9.0]),
        monitor=integrated,
        monitor_kind="integrated",
        reference_monitor=20.0,
        reference_exposure=0.5,
    )

    np.testing.assert_array_equal(first, second)


def test_uncalibrated_fluence_preserves_relative_frame_scale():
    exposure = np.array([0.2, 0.5, 1.0])
    rate = np.array([4.0, 3.0, 2.0])
    integrated = rate * exposure

    np.testing.assert_allclose(
        normalization.relative_frame_fluence(
            exposure, monitor=rate, monitor_kind="rate"
        ),
        exposure * rate,
    )
    np.testing.assert_allclose(
        normalization.relative_frame_fluence(
            exposure, monitor=integrated, monitor_kind="integrated"
        ),
        integrated,
    )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"total_flux": 0.0, "exposure_time": 1.0}, "total_flux"),
        ({"total_flux": 1.0}, "exposure_time"),
        (
            {"total_flux": 1.0, "exposure_time": 1.0, "monitor": 2.0},
            "monitor_kind",
        ),
        (
            {
                "total_flux": 1.0,
                "exposure_time": 1.0,
                "monitor": 2.0,
                "monitor_kind": "rate",
            },
            "reference_monitor",
        ),
        (
            {
                "total_flux": 1.0,
                "monitor": 2.0,
                "monitor_kind": "integrated",
                "reference_monitor": 2.0,
            },
            "reference_exposure",
        ),
        (
            {
                "total_flux": 1.0,
                "monitor": -2.0,
                "monitor_kind": "integrated",
                "reference_monitor": 2.0,
                "reference_exposure": 1.0,
            },
            "monitor must be finite and positive",
        ),
    ],
)
def test_calibrated_fluence_rejects_missing_or_nonphysical_inputs(
    kwargs, message
):
    with pytest.raises(ValueError, match=message):
        normalization.frame_fluence(**kwargs)

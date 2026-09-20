"""Independent forward validation of the calibrated CTR reduction path.

The simulated counts below are assembled directly on calibrated detector rays
from the measurement equation.  In particular, the forward calculation does
not call the production structure-factor, normalization, illumination or
angular-factor helpers.  The production reduction is then required to recover
the explicitly quadrature-averaged input ``|F|^2``.
"""

from dataclasses import dataclass

import numpy as np
import pytest

pyFAI = pytest.importorskip("pyFAI")

from orgui.app.config_data import CorrectionState, CurveCorrectionRecord  # noqa: E402
from orgui.app.integration_corrections import (  # noqa: E402
    apply_stationary_corrections,
    corrected_curve_from_record,
    frame_correction_policy,
    stationary_correction_factors,
    structure_factor_from_policy,
)
from orgui.app.peak1Dintegr import _compute_rocking_integration  # noqa: E402
from orgui.datautils.xrayutils import DetectorCalibration  # noqa: E402
from orgui.datautils.xrayutils.corrections import acceptance, beamprofile  # noqa: E402


ELECTRON_RADIUS_M = 2.8179403262e-15
WAVELENGTH_ANGSTROM = 0.7005
UNIT_CELL_AREA_ANGSTROM2 = 15.8
PIXEL_M = 172e-6
DISTANCE_M = 0.24
SHAPE = (241, 261)


@dataclass
class _Scan:
    """Small scan object exposing the counters used by the policy adapter."""

    exposure_time: np.ndarray
    monitor: np.ndarray

    @property
    def auxillary_counters(self):
        return ("monitor",)

    def __len__(self):
        return self.exposure_time.size


def _detector():
    """Return a calibrated flat detector with nontrivial arm geometry."""
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.detector = pyFAI.detectors.Detector(
        pixel1=PIXEL_M, pixel2=PIXEL_M, max_shape=SHAPE
    )
    detector.poni1 = SHAPE[0] * PIXEL_M / 2.0
    detector.poni2 = SHAPE[1] * PIXEL_M / 2.0
    detector.rot1 = detector.rot2 = detector.rot3 = 0.0
    detector.dist = DISTANCE_M
    detector.set_energy(17.7)
    detector.setAzimuthalReference(np.deg2rad(90.0))
    detector.setPolarization(0.0, 1.0)
    detector.reset()
    detector._cached_array = {}
    return detector


def _prefactor():
    """Independent evaluation of ``r_e**2 lambda**2 / A_u**2``."""
    wavelength_m = WAVELENGTH_ANGSTROM * 1e-10
    area_m2 = UNIT_CELL_AREA_ANGSTROM2 * 1e-20
    return ELECTRON_RADIUS_M**2 * wavelength_m**2 / area_m2**2


def _profile(fwhm_m, offset_m):
    """A densely sampled measured profile, including a sample offset."""
    sigma = fwhm_m / np.sqrt(8.0 * np.log(2.0))
    z = np.linspace(-8.0 * sigma, 8.0 * sigma, 4001)
    density = np.exp(-0.5 * (z / sigma) ** 2)
    return beamprofile.MeasuredBeamProfile(z, density, offset=offset_m)


def _piecewise_integral(x, y, lower, upper):
    """Integrate tabulated linear data without using the beam-profile API."""
    points = np.concatenate(
        ([lower], x[(x > lower) & (x < upper)], [upper])
    )
    values = np.interp(points, x, y, left=0.0, right=0.0)
    return np.trapezoid(values, points)


def _independent_illumination(profile, alpha, length, horizontal):
    """Evaluate ``H = f_z f_x / sin(alpha)`` from the sampled profile."""
    result = []
    for angle in np.atleast_1d(alpha):
        projected = length * np.sin(angle)
        fraction = _piecewise_integral(
            profile.z, profile.density, -projected / 2.0, projected / 2.0
        )
        result.append(fraction * horizontal / np.sin(angle))
    return np.asarray(result)


def _independent_fluence(
    flux, exposure, monitor, monitor_kind, reference, reference_exposure
):
    """Calculate photons per frame from the stated monitor convention."""
    if monitor_kind == "rate":
        return flux * exposure * monitor / reference
    return flux * reference_exposure * monitor / reference


def _state(flux, monitor_kind, reference, reference_exposure, horizontal):
    """Return an explicitly calibrated total-flux correction state."""
    return CorrectionState(
        use_normalization=True,
        use_footprint=True,
        total_incident_flux=flux,
        total_flux_calibrated=True,
        primary_monitor="monitor",
        primary_monitor_kind=monitor_kind,
        primary_monitor_unit=("counts/s" if monitor_kind == "rate" else "counts"),
        monitor_reference_reading=reference,
        monitor_reference_exposure_s=(
            reference_exposure if monitor_kind == "integrated" else None
        ),
        horizontal_interception="fraction",
        horizontal_intercepted_fraction=horizontal,
    )


def _pixel_rays(detector, center, size, alpha, gamma_arm, delta_arm):
    """Return calibrated ray angles and edge-to-edge angular pixel widths."""
    row0, column0 = center
    rows, columns = size
    row = row0 + np.arange(rows) - (rows - 1.0) / 2.0
    column = column0 + np.arange(columns) - (columns - 1.0) / 2.0
    rr, cc = np.meshgrid(row, column, indexing="ij")
    args = (alpha, gamma_arm, delta_arm)
    gamma, delta = detector.surfaceAnglesPoint(rr, cc, *args)
    gamma_low = detector.surfaceAnglesPoint(rr - 0.5, cc, *args)[0]
    gamma_high = detector.surfaceAnglesPoint(rr + 0.5, cc, *args)[0]
    delta_low = detector.surfaceAnglesPoint(rr, cc - 0.5, *args)[1]
    delta_high = detector.surfaceAnglesPoint(rr, cc + 0.5, *args)[1]
    return (
        np.asarray(gamma),
        np.asarray(delta),
        np.abs(np.asarray(gamma_high) - np.asarray(gamma_low)),
        np.abs(np.asarray(delta_high) - np.asarray(delta_low)),
    )


def _polarization(alpha, gamma, delta):
    """Horizontal-polarization factor written independently from the model."""
    projection = (
        np.sin(alpha) * np.cos(delta) * np.cos(gamma)
        + np.cos(alpha) * np.sin(gamma)
    )
    return 1.0 - projection**2


def _f2_profile(gamma, center, sigma):
    """A smooth, finite-width rod whose value varies across the aperture."""
    coordinate = (gamma - center) / sigma
    return 420.0 * (1.0 + 0.08 * coordinate + 0.035 * coordinate**2)


@pytest.mark.parametrize(
    "monitor_kind, roi_size, arms, fwhm_um, offset_um, horizontal",
    [
        ("rate", (41, 39), (0.08, 0.11), 34.0, 0.0, 1.0),
        ("integrated", (9, 11), (0.14, 0.22), 76.0, 19.0, 0.63),
    ],
)
def test_stationary_reduction_recovers_the_pixel_ray_resolution_average(
    monitor_kind, roi_size, arms, fwhm_um, offset_um, horizontal
):
    """Flux, exposure, arm position, ROI clipping and footprint all cancel."""
    detector = _detector()
    center = (145.0, 169.0)
    alpha = np.array([0.018, 0.032, 0.051])
    gamma_arm = arms[0] + np.array([-0.012, 0.0, 0.017])
    delta_arm = arms[1] + np.array([-0.018, 0.0, 0.021])
    exposure = np.array([0.37, 0.82, 1.41])
    reference = 2.4e6
    reference_exposure = 0.65
    if monitor_kind == "rate":
        monitor = reference * np.array([0.72, 1.13, 0.91])
    else:
        monitor = reference * exposure / reference_exposure
        monitor *= np.array([0.72, 1.13, 0.91])
    scan = _Scan(exposure, monitor)
    flux = 3.7e10
    length = 4.2e-3
    profile = _profile(fwhm_um * 1e-6, offset_um * 1e-6)
    state = _state(
        flux, monitor_kind, reference, reference_exposure, horizontal
    )

    policy = frame_correction_policy(
        scan,
        state,
        alpha.size,
        use_normalization=True,
        use_illumination=True,
        alpha=alpha,
        beam_profile=profile,
        sample_length=length,
    )
    q_forward = _independent_fluence(
        flux,
        exposure,
        monitor,
        monitor_kind,
        reference,
        reference_exposure,
    )
    h_forward = _independent_illumination(profile, alpha, length, horizontal)
    np.testing.assert_allclose(policy.normalization_divisor, q_forward, rtol=2e-14)
    np.testing.assert_allclose(policy.illumination_divisor, h_forward, rtol=2e-12)

    raw = []
    expected = []
    center_gamma = []
    center_delta = []
    for index in range(alpha.size):
        gamma, delta, dgamma, ddelta = _pixel_rays(
            detector,
            center,
            roi_size,
            alpha[index],
            gamma_arm[index],
            delta_arm[index],
        )
        gamma0, delta0 = detector.surfaceAnglesPoint(
            np.array([center[0]]),
            np.array([center[1]]),
            alpha[index],
            gamma_arm[index],
            delta_arm[index],
        )
        gamma0, delta0 = float(gamma0[0]), float(delta0[0])
        center_gamma.append(gamma0)
        center_delta.append(delta0)
        sigma_gamma = 3.6 * PIXEL_M / DISTANCE_M
        sigma_delta = 3.1 * PIXEL_M / DISTANCE_M
        resolution = np.exp(-0.5 * ((gamma - gamma0) / sigma_gamma) ** 2)
        resolution *= np.exp(-0.5 * ((delta - delta0) / sigma_delta) ** 2)
        resolution /= 2.0 * np.pi * sigma_gamma * sigma_delta
        polarization = _polarization(alpha[index], gamma, delta)
        f2 = _f2_profile(gamma, gamma0, sigma_gamma)
        lorentz = 1.0 / np.sin(gamma)
        quadrature = dgamma * ddelta
        photon_sum = np.sum(
            f2 * lorentz * polarization * resolution * quadrature
        )
        polarization_correction = np.mean(1.0 / polarization)
        raw.append(
            q_forward[index]
            * h_forward[index]
            * _prefactor()
            * photon_sum
            * polarization_correction
        )
        expected.append(
            photon_sum * polarization_correction * np.sin(gamma0)
        )

    factors = stationary_correction_factors(
        alpha,
        np.asarray(center_delta),
        np.asarray(center_gamma),
        use_lorentz=True,
        normalization=policy.normalization_divisor,
        illumination_divisor=policy.illumination_divisor,
    )
    corrected, errors = apply_stationary_corrections(
        np.asarray(raw), np.zeros(alpha.size), factors
    )
    actual, _ = structure_factor_from_policy(
        corrected,
        errors,
        factors,
        policy,
        wavelength=WAVELENGTH_ANGSTROM,
        unitcell_area=UNIT_CELL_AREA_ANGSTROM2,
    )

    # The forward target is a finite-resolution integral, not F2 at one ray.
    expected = np.asarray(expected)
    np.testing.assert_allclose(actual, expected, rtol=2e-10)
    if min(roi_size) < 15:
        point_f2 = _f2_profile(np.asarray(center_gamma), np.asarray(center_gamma), 1.0)
        assert np.max(np.abs(actual / point_f2 - 1.0)) > 0.05


def _rocking_result(sample_count, roi_size):
    """Forward-simulate and reduce one rocking curve."""
    detector = _detector()
    center = (154.0, 171.0)
    alpha = 0.031
    gamma_arm, delta_arm = 0.105, 0.19
    gamma, delta, dgamma, ddelta = _pixel_rays(
        detector, center, roi_size, alpha, gamma_arm, delta_arm
    )
    gamma0, delta0 = detector.surfaceAnglesPoint(
        np.array([center[0]]),
        np.array([center[1]]),
        alpha,
        gamma_arm,
        delta_arm,
    )
    gamma0, delta0 = float(gamma0[0]), float(delta0[0])
    sigma_delta = 3.0 * PIXEL_M / DISTANCE_M
    sigma_gamma = 4.0 * PIXEL_M / DISTANCE_M
    transverse = np.exp(-0.5 * ((delta - delta0) / sigma_delta) ** 2)
    transverse /= np.sqrt(2.0 * np.pi) * sigma_delta
    polarization = _polarization(alpha, gamma, delta)
    f2 = _f2_profile(gamma, gamma0, sigma_gamma)
    cross_section = np.sum(
        f2 * polarization * transverse * ddelta * dgamma
    )
    polarization_correction = np.mean(1.0 / polarization)

    sigma_omega = np.deg2rad(0.055)
    axis = np.linspace(-4.5 * sigma_omega, 4.5 * sigma_omega, sample_count)
    rocking_profile = np.exp(-0.5 * (axis / sigma_omega) ** 2)
    rocking_profile /= np.sqrt(2.0 * np.pi) * sigma_omega
    exposure = 0.25 + 0.7 * (axis - axis.min()) / np.ptp(axis)
    reference = 1.8e6
    monitor = reference * (1.0 + 0.17 * np.sin(axis / sigma_omega))
    scan = _Scan(exposure, monitor)
    flux = 8.1e9
    reference_exposure = 0.5
    horizontal = 0.71
    length = 3.8e-3
    profile = _profile(58e-6, -13e-6)
    state = _state(flux, "rate", reference, reference_exposure, horizontal)
    policy = frame_correction_policy(
        scan,
        state,
        sample_count,
        use_normalization=True,
        use_illumination=True,
        alpha=np.full(sample_count, alpha),
        beam_profile=profile,
        sample_length=length,
    )
    q_forward = _independent_fluence(
        flux, exposure, monitor, "rate", reference, reference_exposure
    )
    h_forward = _independent_illumination(
        profile, np.full(sample_count, alpha), length, horizontal
    )
    lorentz = 1.0 / (
        np.sin(delta0) * np.cos(alpha) * np.cos(gamma0)
    )
    rod = np.cos(gamma0)
    base_curve = (
        q_forward
        * h_forward
        * _prefactor()
        * lorentz
        * rod
        * rocking_profile
        * cross_section
        * polarization_correction
    )
    record = CurveCorrectionRecord(
        algorithm="framewise_ctr_total_flux_v1",
        output_quantity="ctr_photon_curve",
        scale_convention="total_flux_calibrated",
        normalization_status="applied",
        illumination_status="applied",
        pixel_correction_status="polarization_only",
        normalization_divisor=policy.normalization_divisor,
        normalization_unit="photons",
        illumination_divisor=policy.illumination_divisor,
        illumination_convention="total_flux_H",
        base_croibg=base_curve,
        base_croibg_variance=np.zeros(sample_count),
    )
    corrected, corrected_errors, _status, _convention = (
        corrected_curve_from_record(record)
    )
    detector_acceptance = float(
        np.asarray(
            acceptance.out_of_plane_acceptance(
                detector,
                center[0],
                center[1],
                roi_size[0],
                alpha,
                gamma_arm,
                delta_arm,
            )
        ).reshape(-1)[0]
    )
    result = _compute_rocking_integration(
        np.array([0.0]),
        axis,
        corrected[np.newaxis, :],
        corrected_errors[np.newaxis, :],
        {"sig_1": {"from": np.array([axis[0]]), "to": np.array([axis[-1]])}},
        {},
        True,
        False,
        C_Lor=np.full((1, sample_count), lorentz),
        C_rod=np.full((1, sample_count), rod),
        detector_acceptance=np.array([detector_acceptance]),
        angle_unit="rad",
    )
    actual = result["F2_hkl"][0] / _prefactor()
    expected = (
        np.trapezoid(rocking_profile, axis)
        * cross_section
        * polarization_correction
        / detector_acceptance
    )
    return actual, expected


@pytest.mark.parametrize("roi_size", [(41, 41), (11, 7)])
def test_rocking_reduction_recovers_full_and_clipped_pixel_ray_integrals(roi_size):
    """A clipped transverse peak is compared with its captured integral."""
    actual, expected = _rocking_result(257, roi_size)
    np.testing.assert_allclose(actual, expected, rtol=3e-11)


def test_rocking_sampling_converges_before_the_validation_tolerance_is_set():
    """The chosen grid is demonstrably inside its numerical convergence regime."""
    reference = _rocking_result(2049, (31, 31))[0]
    approximations = np.array(
        [_rocking_result(count, (31, 31))[0] for count in (9, 17, 33, 65)]
    )
    errors = np.abs(approximations - reference)
    assert np.all(np.diff(errors) < 0.0)
    assert errors[-1] / reference < 3e-7

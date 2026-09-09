"""Regression tests for the common structure-factor scale.

:mod:`orgui.datautils.xrayutils.corrections.measurement` reduces an integrated
intensity to :math:`|F_{hkl}|^2` for any scan mode. Every expectation here is
written out independently from the published equations -- E. Vlieg,
*J. Appl. Cryst.* **30** (1997) 532 and J. Drnec *et al.*,
*J. Appl. Cryst.* **47** (2014) 365 -- rather than by calling the module back,
so that the tests pin the physics and not the implementation.

See ``doc/design/ctr_structure_factor_scale.md`` for the analysis these tests
belong to, and `issue #15 <https://github.com/tifuchs/orGUI/issues/15>`_ and
`issue #82 <https://github.com/tifuchs/orGUI/issues/82>`_ for the requests.
"""

import numpy as np
import pytest

from orgui.datautils.xrayutils.corrections import measurement as ii

#: A z-axis trajectory along a rod: fixed incidence, exit angle climbing.
ALPHA = np.deg2rad(np.full(6, 0.6))
DELTA = np.deg2rad(np.array([16.77, 16.79, 16.83, 16.89, 16.93, 16.95]))
GAMMA = np.deg2rad(np.array([1.76, 3.53, 5.31, 8.88, 12.49, 17.38]))

#: Out-of-plane acceptance of a region of interest, in radian.
DGAMMA = np.deg2rad(0.35)

#: Pt(111) hexagonal surface cell at 17.7 keV.
WAVELENGTH = 0.70047
UC_AREA = 2.7748 * 2.7748 * np.sin(np.deg2rad(120.0))


def test_normalized_intensity_divides_by_time_and_monitor():
    """Both published expressions carry ``Phi_0 T``; both must come out."""
    got = ii.normalized_intensity(
        np.array([100.0, 200.0]), exposure_time=0.5, monitor=4.0
    )

    np.testing.assert_allclose(got, np.array([50.0, 100.0]))


def test_normalized_intensity_converts_the_rocking_angle_to_radian():
    """A rocking integral taken in degrees is 180/pi too large.

    :func:`~orgui.app.peak1Dintegr._compute_rocking_integration` integrates
    over the rocking axis as stored, which is degrees. Vlieg's equation 42
    and Drnec's equation 2 integrate in radian.
    """
    counts_deg = 573.0

    in_deg = ii.normalized_intensity(counts_deg, angle_unit="deg")
    in_rad = ii.normalized_intensity(counts_deg, angle_unit="rad")

    np.testing.assert_allclose(in_deg, np.deg2rad(counts_deg))
    np.testing.assert_allclose(in_rad / in_deg, np.rad2deg(1.0))


def test_normalized_intensity_rejects_unusable_divisors():
    """A zero monitor or exposure must fail loudly, not scale by infinity."""
    with pytest.raises(ValueError, match="finite and positive"):
        ii.normalized_intensity(1.0, exposure_time=0.0)
    with pytest.raises(ValueError, match="finite and non-zero"):
        ii.normalized_intensity(1.0, monitor=0.0)
    with pytest.raises(ValueError, match="unknown angle_unit"):
        ii.normalized_intensity(1.0, angle_unit="degrees")


def test_rocking_factor_is_the_published_z_axis_expression():
    """``L_phi * C_rod * Delta_gamma``, Vlieg equations 16, 20, 23 and 42."""
    expected = (
        1.0 / (np.sin(DELTA) * np.cos(ALPHA) * np.cos(GAMMA))  # L_phi, eq (16)
        * np.cos(GAMMA)  # C_rod, eq (23) in the z-axis mode
        * DGAMMA  # rod length intercepted, eq (20)
    )

    got = ii.angular_factor(
        ii.ROCKING,
        alpha=ALPHA,
        delta=DELTA,
        gamma=GAMMA,
        detector_acceptance=DGAMMA,
    )

    np.testing.assert_allclose(got, expected, rtol=1e-12)


def test_stationary_factor_is_the_published_z_axis_expression():
    """``L_s = 1/sin(gamma)``, Vlieg equation 53 and Drnec equation 5."""
    got = ii.angular_factor(ii.STATIONARY, gamma=GAMMA)

    np.testing.assert_allclose(got, 1.0 / np.sin(GAMMA), rtol=1e-12)


def test_rocking_intensity_is_linear_in_the_detector_acceptance():
    """Twice the out-of-plane acceptance intercepts twice as much rod.

    This is the factor that makes a rocking scan disagree with a stationary
    measurement of the same rod, and it is not constant over a scan whenever
    the region of interest is resized per detector position.
    """
    single = ii.angular_factor(
        ii.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA, detector_acceptance=DGAMMA
    )
    double = ii.angular_factor(
        ii.ROCKING,
        alpha=ALPHA,
        delta=DELTA,
        gamma=GAMMA,
        detector_acceptance=2.0 * DGAMMA,
    )

    np.testing.assert_allclose(double, 2.0 * single, rtol=1e-12)


def test_the_acceptance_belongs_to_rocking_scans_only():
    """Neither mode may silently accept the other's arguments."""
    with pytest.raises(ValueError, match="out-of-plane acceptance"):
        ii.angular_factor(ii.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA)
    with pytest.raises(ValueError, match="whole rod cross"):
        ii.angular_factor(ii.STATIONARY, gamma=GAMMA, detector_acceptance=DGAMMA)
    with pytest.raises(ValueError, match="unknown mode"):
        ii.angular_factor("rsm", gamma=GAMMA)


def test_specular_is_the_stationary_factor_at_the_specular_condition():
    """Vlieg section 3.2: reflectivity is the stationary case, gamma = alpha."""
    alpha = np.deg2rad(np.array([0.4, 1.0, 2.5]))

    specular = ii.angular_factor(ii.SPECULAR, alpha=alpha)

    np.testing.assert_allclose(specular, 1.0 / np.sin(alpha), rtol=1e-12)
    np.testing.assert_allclose(
        specular, ii.angular_factor(ii.STATIONARY, gamma=alpha), rtol=1e-12
    )
    with pytest.raises(ValueError, match="non-specular exit angle"):
        ii.angular_factor(ii.SPECULAR, alpha=alpha, gamma=alpha + 0.1)


def test_scale_factor_carries_its_documented_units():
    """``Phi_0 r_e^2 A lambda^2 / A_u^2``, with lengths mixed as documented.

    The wavelength and the unit-cell area are in Angstrom, the area and the
    flux density in meter; getting that conversion wrong is a factor of
    1e20 and is exactly the kind of error the module exists to prevent.
    """
    r_e = 2.8179403262e-15
    expected = (
        1e16 * r_e**2 * 5.73e-7 * (WAVELENGTH * 1e-10) ** 2 / (UC_AREA * 1e-20) ** 2
    )

    got = ii.scale_factor(
        WAVELENGTH, UC_AREA, active_area=5.73e-7, flux_density=1e16
    )

    np.testing.assert_allclose(got, expected, rtol=1e-12)
    # The prefactor is a rate: counts per second for |F|^2 = 1.
    assert 1e-4 < got < 1e2


def test_scale_factor_rejects_nonphysical_lengths():
    """A negative wavelength or unit cell is a units mistake, not a value."""
    with pytest.raises(ValueError, match="wavelength must be positive"):
        ii.scale_factor(-1.0, UC_AREA)
    with pytest.raises(ValueError, match="unitcell_area must be positive"):
        ii.scale_factor(WAVELENGTH, 0.0)


def test_structure_factor_round_trips_through_the_forward_model():
    """The reduction is the exact inverse of the forward model, per mode."""
    f2 = np.array([273.6, 148.3, 108.2, 98.0, 865.9, 5000.0])
    common = dict(
        wavelength=WAVELENGTH,
        unitcell_area=UC_AREA,
        active_area=5.73e-7,
        flux_density=1e16,
        detector_efficiency=0.93,
    )

    for kwargs in (
        dict(mode=ii.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA,
             detector_acceptance=DGAMMA),
        dict(mode=ii.STATIONARY, gamma=GAMMA),
        dict(mode=ii.SPECULAR, alpha=GAMMA),
    ):
        intensity = ii.integrated_intensity(f2, **kwargs, **common)
        back = ii.structure_factor_squared(intensity, **kwargs, **common)
        np.testing.assert_allclose(back, f2, rtol=1e-12)


def test_structure_factor_needs_the_scale_to_be_stated():
    """The wavelength and unit-cell area are not optional defaults."""
    with pytest.raises(ValueError, match="set the scale"):
        ii.structure_factor_squared(1.0, ii.STATIONARY, gamma=GAMMA)


def test_the_two_modes_agree_on_one_structure_factor():
    """The point of issue #82: one rod, two scan modes, one ``|F|^2``.

    A rocking scan and a stationary measurement of the same reflection are
    forward-simulated with completely different counting times, monitor
    values and detector acceptances, then reduced. Only if every
    mode-dependent factor is right do the two land on the same number.
    """
    f2 = np.array([273.6, 148.3, 108.2, 98.0, 865.9, 5000.0])
    scale = dict(
        wavelength=WAVELENGTH,
        unitcell_area=UC_AREA,
        active_area=5.73e-7,
        flux_density=1e16,
    )
    # Different acceptance at every point, as a resized region of interest
    # gives; if it were not divided out this test could not pass.
    acceptance = DGAMMA * np.linspace(0.7, 1.6, GAMMA.size)

    rocking_counts = (
        ii.integrated_intensity(
            f2, ii.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA,
            detector_acceptance=acceptance, **scale
        )
        * 0.4  # per-frame counting time, s
        * 3.0  # monitor
        * np.rad2deg(1.0)  # the rocking axis was integrated in degrees
    )
    stationary_counts = (
        ii.integrated_intensity(f2, ii.STATIONARY, gamma=GAMMA, **scale)
        * 2.0  # counting time, s
        * 7.0  # monitor
    )

    f2_rocking = ii.structure_factor_squared(
        ii.normalized_intensity(
            rocking_counts, exposure_time=0.4, monitor=3.0, angle_unit="deg"
        ),
        ii.ROCKING,
        alpha=ALPHA,
        delta=DELTA,
        gamma=GAMMA,
        detector_acceptance=acceptance,
        **scale,
    )
    f2_stationary = ii.structure_factor_squared(
        ii.normalized_intensity(stationary_counts, exposure_time=2.0, monitor=7.0),
        ii.STATIONARY,
        gamma=GAMMA,
        **scale,
    )

    np.testing.assert_allclose(f2_rocking, f2, rtol=1e-12)
    np.testing.assert_allclose(f2_stationary, f2, rtol=1e-12)


def test_mode_ratio_matches_vlieg_equation_65():
    """The published ratio of the two integrated intensities.

    Vlieg equation 65 for the z-axis mode,
    ``I_s / I_phi = T omega_0 sin(delta) cos(beta_in) / (Delta_gamma sin(gamma))``.
    Written here for intensities already divided by ``T`` and by the rotation
    speed, the ``T omega_0`` drops out and what remains is a statement about
    the two angular factors alone.
    """
    rocking = ii.angular_factor(
        ii.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA, detector_acceptance=DGAMMA
    )
    stationary = ii.angular_factor(ii.STATIONARY, gamma=GAMMA)

    expected = np.sin(DELTA) * np.cos(ALPHA) / (np.sin(GAMMA) * DGAMMA)

    np.testing.assert_allclose(stationary / rocking, expected, rtol=1e-12)


def test_reflectivity_reproduces_the_fresnel_asymptote():
    """Vlieg equation 63 against the textbook far-field Fresnel limit.

    For a semi-infinite substrate far from a bulk Bragg peak the truncation
    rod amplitude is ``|F| = rho_e A_u / q_z``, and the kinematic
    reflectivity of such a substrate is ``(q_c / 2 q_z)^4`` with
    ``q_c^2 = 16 pi r_e rho_e``. Any error in the wavelength, unit-cell area
    or ``sin`` powers of equation 63 breaks this identity, so it validates
    the absolute scale end to end and not just an algebraic rearrangement.
    """
    r_e_ang = ii.CLASSICAL_ELECTRON_RADIUS / 1e-10  # Angstrom
    rho_e = 5.16  # electrons per cubic Angstrom, close to Pt
    alpha = np.deg2rad(np.array([0.5, 0.8, 1.2, 2.0, 3.0]))
    q_z = 4.0 * np.pi * np.sin(alpha) / WAVELENGTH  # 1/Angstrom
    q_c = np.sqrt(16.0 * np.pi * r_e_ang * rho_e)

    f2 = (rho_e * UC_AREA / q_z) ** 2
    reflectivity = ii.reflectivity_from_structure_factor(
        f2, WAVELENGTH, UC_AREA, alpha
    )

    np.testing.assert_allclose(reflectivity, (q_c / (2.0 * q_z)) ** 4, rtol=1e-12)
    # Well above the critical angle, where the kinematic result is valid.
    assert np.all(np.rad2deg(np.arcsin(q_c * WAVELENGTH / (4.0 * np.pi))) < 0.5)


def test_reflectivity_round_trips_and_generalizes_off_specular():
    """The inverse recovers ``|F|^2``, and a non-specular exit angle works."""
    alpha = np.deg2rad(0.6)
    beta_out = np.deg2rad(np.array([0.6, 2.0, 8.0]))
    f2 = np.array([1500.0, 900.0, 120.0])

    reflectivity = ii.reflectivity_from_structure_factor(
        f2, WAVELENGTH, UC_AREA, alpha, beta_out=beta_out, polarization=0.98
    )
    back = ii.structure_factor_from_reflectivity(
        reflectivity, WAVELENGTH, UC_AREA, alpha, beta_out=beta_out, polarization=0.98
    )

    np.testing.assert_allclose(back, f2, rtol=1e-12)
    # At the specular condition the two sines collapse to the familiar
    # 1/sin^2(alpha) of Vlieg equation 63.
    specular = ii.reflectivity_from_structure_factor(
        f2[0], WAVELENGTH, UC_AREA, alpha
    )
    prefactor = (
        ii.CLASSICAL_ELECTRON_RADIUS**2
        * (WAVELENGTH * 1e-10) ** 2
        / (UC_AREA * 1e-20) ** 2
    )
    np.testing.assert_allclose(
        specular, prefactor * f2[0] / np.sin(alpha) ** 2, rtol=1e-12
    )


def test_footprint_area_cancels_between_the_modes():
    """The active area cannot be why two scan modes disagree.

    It enters both expressions identically, so it drops out of their ratio
    even when it is wrong. That is why it is listed under the absolute scale
    (issue #15) and not under the mode equivalence (issue #82).
    """
    kwargs = dict(wavelength=WAVELENGTH, unitcell_area=UC_AREA, flux_density=1e16)
    ratios = []
    for area in (1e-7, 5.73e-7, 2e-6):
        rocking = ii.integrated_intensity(
            1.0, ii.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA,
            detector_acceptance=DGAMMA, active_area=area, **kwargs
        )
        stationary = ii.integrated_intensity(
            1.0, ii.STATIONARY, gamma=GAMMA, active_area=area, **kwargs
        )
        ratios.append(rocking / stationary)

    for ratio in ratios[1:]:
        np.testing.assert_allclose(ratio, ratios[0], rtol=1e-12)

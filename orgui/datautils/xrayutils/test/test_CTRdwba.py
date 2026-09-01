"""Focused validation of the native bulk DWBA path."""

from unittest import mock

import numpy as np
import pytest

from .. import CTRcalc, CTRuc, CTRutil, HKLVlieg
from .. import CTRoptics
from ..CTRoptics import HC_KEV_ANGSTROM


def _kinematics(cell, energy_eV, alpha_i, alpha_f, azimuth):
    wavelength = HC_KEV_ANGSTROM / (energy_eV * 1e-3)
    k0 = 2.0 * np.pi / wavelength
    incident_parallel = k0 * np.cos(alpha_i)
    exit_parallel = k0 * np.cos(alpha_f)
    Q_parallel = np.sqrt(
        incident_parallel**2
        + exit_parallel**2
        - 2.0 * incident_parallel * exit_parallel * np.cos(azimuth)
    )
    Qz = k0 * (np.sin(alpha_i) + np.sin(alpha_f))
    h = Q_parallel * cell.a[0] / (2.0 * np.pi)
    l_value = Qz * cell.a[2] / (2.0 * np.pi)
    return h, l_value, k0


def _one_atom_cell(occupancy=1.0):
    cell = CTRuc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0])
    cell.addAtom("C", [0.17, 0.23, 0.31], 0.18, 0.37, occupancy)
    cell.setEnergy(10000.0)
    cell.f[0, 11] = 0.27
    cell.f[0, 12] = 0.19
    return cell


@pytest.mark.parametrize(
    ("polarization_i", "polarization_f", "contraction"),
    [
        ("s", "s", lambda ai, af, az: np.cos(az)),
        ("s", "p", lambda ai, af, az: -np.sin(af) * np.sin(az)),
        ("p", "s", lambda ai, af, az: -np.sin(ai) * np.sin(az)),
        (
            "p",
            "p",
            lambda ai, af, az: (
                np.cos(ai) * np.cos(af)
                - np.sin(ai) * np.sin(af) * np.cos(az)
            ),
        ),
    ],
)
def test_vacuum_born_limit_matches_renaud_contractions(
    polarization_i, polarization_f, contraction
):
    cell = _one_atom_cell(occupancy=0.0)
    crystal = CTRcalc.SXRDCrystal(cell)
    alpha_i = 0.21
    alpha_f = 0.34
    azimuth = 0.43
    h, l_value, _ = _kinematics(
        cell, cell._E, alpha_i, alpha_f, azimuth
    )
    prepared = crystal.prepare_DWBA(
        h,
        0.0,
        l_value,
        alpha_i,
        alpha_f,
        polarization_i=polarization_i,
        polarization_f=polarization_f,
    )
    np.testing.assert_array_equal(prepared.field_i.n, np.ones(2))
    np.testing.assert_array_equal(prepared.field_i.A_minus[-1], 0.0)
    np.testing.assert_array_equal(prepared.field_f.A_minus[-1], 0.0)

    cell.basis[0, 6] = 1.0
    actual = crystal.F_DWBA_prepared(prepared, bulk_mode="unit_cell")
    ordinary = cell.F_uc(
        np.array([h]), np.array([0.0]), np.array([l_value])
    )[0]
    expected = ordinary * contraction(alpha_i, alpha_f, azimuth)
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
    assert isinstance(actual, complex)
    attenuation = 1.0e-10
    bulk = crystal.F_DWBA_prepared(
        prepared,
        bulk_mode="semi_infinite",
        bulk_attenuation=attenuation,
    )
    ordinary_bulk = cell.F_bulk(
        np.array([h]),
        np.array([0.0]),
        np.array([l_value]),
        attenuation,
    )[0]
    np.testing.assert_allclose(
        bulk,
        ordinary_bulk * contraction(alpha_i, alpha_f, azimuth),
        rtol=1e-9,
        atol=2e-13,
    )


def test_complex_channel_oracle_domains_anomalous_and_anisotropic_dw():
    cell = _one_atom_cell()
    identity = np.hstack((np.identity(3), np.zeros((3, 1))))
    translated = identity.copy()
    translated[:, 3] = [0.11, -0.07, 0.03]
    cell.coherentDomainMatrix = [identity, translated]
    cell.coherentDomainOccupancy = [0.65, 0.35]
    crystal = CTRcalc.SXRDCrystal(cell)
    alpha_i = 0.035
    alpha_f = 0.052
    azimuth = 0.28
    h, l_value, _ = _kinematics(
        cell, cell._E, alpha_i, alpha_f, azimuth
    )
    prepared = crystal.prepare_DWBA(
        h, 0.0, l_value, alpha_i, alpha_f
    )

    field_i = prepared.field_i
    field_f = prepared.field_f
    kz_i = field_i.kz[-1, 0]
    kz_f = field_f.kz[-1, 0]
    Qz = -(kz_i + kz_f)
    Qx, Qy = prepared.Q_parallel[0]
    Q_parallel_squared = Qx**2 + Qy**2
    Q_squared = Q_parallel_squared + Qz**2
    assert abs(Q_squared.imag) > 0.0
    row = cell.f[0]
    form_factor = row[10] + row[11] + 1j * row[12]
    for term in range(5):
        form_factor += row[term] * np.exp(-row[term + 5] * Q_squared)
    basis = cell.basis[0]
    form_factor *= np.exp(
        -(
            basis[4] * Q_parallel_squared + basis[5] * Qz**2
        )
        / (16.0 * np.pi**2)
    )
    form_factor *= basis[6]
    domain_sum = 0.0j
    attenuated_domain_sum = 0.0j
    attenuation = 3.0e-7
    for matrix, occupancy in zip(
        cell.coherentDomainMatrix, cell.coherentDomainOccupancy
    ):
        effective = cell.R_mat_inv @ matrix[:, :3] @ cell.R_mat
        fractional_position = effective @ basis[1:4] + matrix[:, 3]
        position = cell.R_mat @ fractional_position
        phase = np.exp(
            1j * (Qx * position[0] + Qy * position[1] + Qz * position[2])
        )
        domain_sum += occupancy * phase
        attenuated_domain_sum += (
            occupancy
            * np.exp(attenuation * fractional_position[2])
            * phase
        )
    coefficient = (
        field_i.A_plus[-1, 0]
        * field_f.A_plus[-1, 0]
        * prepared.cos_azimuth[0]
        * np.exp(
            1j
            * (
                kz_i * field_i.z_reference[-1]
                + kz_f * field_f.z_reference[-1]
            )
        )
    )
    expected_unit = coefficient * form_factor * domain_sum

    CTRuc._CTRcalc_cpp.clear_form_factor_cache()
    CTRuc._CTRcalc_cpp.reset_form_factor_cache_stats()
    actual_unit = crystal.F_DWBA_prepared(
        prepared, bulk_mode="unit_cell"
    )
    first_stats = CTRuc._CTRcalc_cpp.form_factor_cache_stats()
    repeated = crystal.F_DWBA_prepared(prepared, bulk_mode="unit_cell")
    second_stats = CTRuc._CTRcalc_cpp.form_factor_cache_stats()
    np.testing.assert_allclose(actual_unit, expected_unit, rtol=2e-13)
    np.testing.assert_allclose(repeated, actual_unit, rtol=0.0, atol=0.0)
    assert second_stats["hits"] > first_stats["hits"]

    repeat_phase = (
        Qx * cell.R_mat[0, 2]
        + Qy * cell.R_mat[1, 2]
        + Qz * cell.R_mat[2, 2]
    )
    expected_bulk = (
        coefficient * form_factor * attenuated_domain_sum
    ) / (
        -np.expm1(-1j * repeat_phase - attenuation)
    )
    actual_bulk = crystal.F_DWBA_prepared(
        prepared,
        bulk_mode="semi_infinite",
        bulk_attenuation=attenuation,
    )
    np.testing.assert_allclose(actual_bulk, expected_bulk, rtol=2e-13)
    one_call = crystal.F_DWBA(
        h,
        0.0,
        l_value,
        alpha_i,
        alpha_f,
        bulk_mode="unit_cell",
    )
    np.testing.assert_allclose(one_call, actual_unit, rtol=0.0, atol=0.0)


def test_empirical_attenuation_is_continuous_at_the_bulk_repeat():
    cell = CTRuc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0])
    cell.addAtom("C", [0.0, 0.0, -1.0], 0.0, 0.0, 0.0)
    cell.addAtom("C", [0.0, 0.0, -0.5], 0.0, 0.0, 0.0)
    cell.setEnergy(10000.0)
    crystal = CTRcalc.SXRDCrystal(cell)

    wavelength = HC_KEV_ANGSTROM / (cell._E * 1e-3)
    k0 = 2.0 * np.pi / wavelength
    alpha = np.arcsin(np.pi / (cell.a[2] * k0))
    prepared = crystal.prepare_DWBA(0.0, 0.0, 1.0, alpha, alpha)
    cell.basis[:, 6] = 1.0

    attenuation = 0.01
    actual = crystal.F_DWBA_prepared(
        prepared,
        bulk_attenuation=attenuation,
    )
    expected = cell.F_bulk(
        np.array([0.0]),
        np.array([0.0]),
        np.array([1.0]),
        attenuation,
    )[0]
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
    assert abs(actual) > 1.0


def test_dwba_empirical_attenuation_defaults_to_zero():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    alpha_i = 0.18
    alpha_f = 0.27
    h, l_value, _ = _kinematics(cell, cell._E, alpha_i, alpha_f, 0.2)
    prepared = crystal.prepare_DWBA(
        h, 0.0, l_value, alpha_i, alpha_f
    )

    default = crystal.F_DWBA_prepared(prepared)
    explicit_zero = crystal.F_DWBA_prepared(
        prepared, bulk_attenuation=0.0
    )
    np.testing.assert_array_equal(default, explicit_zero)


def test_dwba_equivalent_kinematical_attenuation_and_setter():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell, atten=0.0)
    alpha_i = 0.035
    alpha_f = np.array([0.02, 0.08, 0.3])

    actual = CTRutil.attenuation_from_dwba(
        crystal, alpha_i, alpha_f
    )
    field_i = crystal.wavefield(alpha_i)
    field_f = crystal.wavefield(alpha_f)
    repeat_scale = abs(cell.refHKLTransform[2, 2])
    expected = abs(cell.R_mat[2, 2]) * (
        np.imag(field_i.kz[-1]) + np.imag(field_f.kz[-1])
    ) / repeat_scale
    np.testing.assert_allclose(actual, expected, rtol=2e-15)
    assert actual.shape == alpha_f.shape

    incident_only = CTRutil.set_atten_from_dwba(crystal, alpha_i)
    expected_incident = (
        abs(cell.R_mat[2, 2])
        * np.imag(field_i.kz[-1])
        / repeat_scale
    )
    assert incident_only == pytest.approx(expected_incident)
    assert crystal.atten == incident_only
    assert np.all(actual > incident_only)

    with pytest.raises(ValueError, match="requires scalar angles"):
        CTRutil.set_atten_from_dwba(crystal, alpha_i, alpha_f)


def test_broadcasting_preparation_reuse_and_frozen_geometry():
    cell = _one_atom_cell(occupancy=0.0)
    crystal = CTRcalc.SXRDCrystal(cell)
    alpha_i = np.array([[0.18], [0.24]])
    alpha_f = np.array([[0.27, 0.31, 0.35]])
    azimuth = 0.2
    h, l_value, _ = _kinematics(
        cell, cell._E, alpha_i, alpha_f, azimuth
    )
    with mock.patch.object(
        CTRcalc,
        "solve_electric_field",
        wraps=CTRcalc.solve_electric_field,
    ) as solve:
        prepared = crystal.prepare_DWBA(
            h, 0.0, l_value, alpha_i, alpha_f
        )
        assert solve.call_count == 2
        cell.basis[0, 6] = 1.0
        result = crystal.F_DWBA_prepared(
            prepared, bulk_mode="unit_cell"
        )
        assert solve.call_count == 2
    assert result.shape == (2, 3)
    with pytest.raises(ValueError, match="read-only"):
        prepared.alpha_i[0] = 0.5
    cell.setEnergy(11000.0)
    with pytest.raises(ValueError, match="energy changed"):
        crystal.F_DWBA_prepared(prepared)


def test_geometry_validation_and_exact_zero_attenuation_pole():
    cell = _one_atom_cell(occupancy=0.0)
    crystal = CTRcalc.SXRDCrystal(cell)
    with pytest.raises(ValueError, match="normal geometry"):
        crystal.prepare_DWBA(0.0, 0.0, 0.2, 0.1, 0.1)

    wavelength = HC_KEV_ANGSTROM / (cell._E * 1e-3)
    k0 = 2.0 * np.pi / wavelength
    alpha = np.arcsin(np.pi / (cell.a[2] * k0))
    prepared = crystal.prepare_DWBA(0.0, 0.0, 1.0, alpha, alpha)
    cell.basis[0, 6] = 1.0
    with pytest.raises(ValueError, match="exact Bragg pole"):
        crystal.F_DWBA_prepared(
            prepared,
            bulk_mode="semi_infinite",
            bulk_attenuation=0.0,
        )


def test_bulk_scope_and_missing_native_extension_errors():
    bulk = _one_atom_cell()
    surface = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(bulk, surface)
    with pytest.raises(NotImplementedError, match="only the bulk"):
        crystal.prepare_DWBA(0.0, 0.0, 0.1, 0.1, 0.1)

    bulk_only = CTRcalc.SXRDCrystal(bulk)
    with mock.patch.object(CTRcalc, "HAS_CPP_ACCEL", False):
        with pytest.raises(RuntimeError, match="native CTR extension"):
            bulk_only.prepare_DWBA(0.0, 0.0, 0.1, 0.1, 0.1)
    profile = np.array([[-1.0, 1.0e-5, 0.0], [1.0, 0.0, 0.0]])
    with mock.patch.object(CTRoptics, "HAS_CPP_ACCEL", False):
        with pytest.raises(RuntimeError, match="native CTR extension"):
            CTRoptics.solve_electric_field(profile, 10000.0, 0.1)


def _reference_zmode_angles(cell, hkl, fixedangle, fixed, **keyargs):
    """Independent z-mode solution for the same crystal and wavelength."""
    ub_calculator = HKLVlieg.UBCalculator(cell, cell._E * 1e-3)
    ub_calculator.setLambda(HC_KEV_ANGSTROM / (cell._E * 1e-3))
    ub_calculator.setU(np.identity(3))
    return HKLVlieg.VliegAngles(ub_calculator).anglesZmode(
        np.asarray(hkl, dtype=np.float64), fixedangle, fixed=fixed, **keyargs
    )


def test_zmode_matches_explicit_angles_and_ewald_normal_geometry():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    L = np.array([0.7, 1.0, 1.3])
    hkl = np.vstack((np.ones_like(L), np.zeros_like(L), L))
    fixedangle = 0.05

    prepared, angles = crystal.prepare_DWBA_Zmode(
        1.0, 0.0, L, fixedangle, fixed="in", return_angles=True
    )
    reference = _reference_zmode_angles(cell, hkl, fixedangle, "in")
    assert angles.shape == (3, 6)
    np.testing.assert_allclose(angles, reference, rtol=0.0, atol=1e-14)

    alpha_i = angles[:, 0]
    alpha_f = angles[:, 2]
    np.testing.assert_allclose(alpha_i, fixedangle, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(prepared.alpha_i, alpha_i)
    np.testing.assert_allclose(prepared.alpha_f, alpha_f)

    # The identity orientation matrix must place the surface normal on the
    # omega axis, so the z-mode angles reproduce the crystal-frame Qz.
    k0 = 2.0 * np.pi / (HC_KEV_ANGSTROM / (cell._E * 1e-3))
    Qz = (cell.B_mat @ (cell.refHKLTransform @ hkl))[2]
    np.testing.assert_allclose(
        k0 * (np.sin(alpha_i) + np.sin(alpha_f)), Qz, rtol=1e-12
    )

    np.testing.assert_allclose(
        crystal.F_DWBA_prepared(prepared),
        crystal.F_DWBA(1.0, 0.0, L, alpha_i, alpha_f),
        rtol=0.0,
        atol=0.0,
    )


def test_zmode_equal_angle_constraint_is_specular():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    prepared, angles = crystal.prepare_DWBA_Zmode(
        0.0, 0.0, 1.0, 0.0, fixed="eq", return_angles=True
    )
    alpha, delta, gamma, omega, chi, phi = angles
    assert isinstance(angles, tuple) and len(angles) == 6
    assert alpha == pytest.approx(gamma)
    assert delta == pytest.approx(0.0, abs=1e-12)
    assert (chi, phi) == (0.0, 0.0)
    assert prepared.alpha_i.shape == (1,)
    assert prepared.scalar

    # The fixed exit-angle constraint must mirror the fixed incidence one.
    _, angles_in = crystal.prepare_DWBA_Zmode(
        1.0, 0.0, 1.0, 0.05, fixed="in", return_angles=True
    )
    _, angles_out = crystal.prepare_DWBA_Zmode(
        1.0, 0.0, 1.0, angles_in[2], fixed="out", return_angles=True
    )
    assert angles_out[0] == pytest.approx(angles_in[0])
    assert angles_out[2] == pytest.approx(angles_in[2])


def test_zmode_broadcasting_and_single_batched_field_solve():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    h = np.array([[1.0], [1.2]])
    L = np.array([[0.9, 1.0, 1.1]])
    with mock.patch.object(
        CTRcalc,
        "solve_electric_field",
        wraps=CTRcalc.solve_electric_field,
    ) as solve:
        prepared, angles = crystal.prepare_DWBA_Zmode(
            h, 0.0, L, 0.05, fixed="in", return_angles=True
        )
        assert solve.call_count == 2
        result = crystal.F_DWBA_prepared(prepared)
        assert solve.call_count == 2
    assert result.shape == (2, 3)
    assert angles.shape == (2, 3, 6)
    np.testing.assert_allclose(angles[..., 0], 0.05, rtol=0.0, atol=1e-15)

    mirrored = crystal.prepare_DWBA_Zmode(
        h, 0.0, L, 0.05, fixed="in", mirrorx=True, return_angles=True
    )[1]
    np.testing.assert_allclose(mirrored[..., 1], -angles[..., 1])
    np.testing.assert_allclose(mirrored[..., 0], angles[..., 0])
    np.testing.assert_allclose(mirrored[..., 2], angles[..., 2])


def test_zmode_rejects_unsupported_constraints_and_unreachable_angles():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    with pytest.raises(ValueError, match="fixed must be one of"):
        crystal.prepare_DWBA_Zmode(0.0, 0.0, 1.0, 0.0, fixed="azimuth")
    with pytest.raises(ValueError, match="chi=phi=0"):
        crystal.prepare_DWBA_Zmode(0.0, 0.0, 1.0, 0.05, chi=0.01)
    with pytest.raises(ValueError, match="chi=phi=0"):
        crystal.prepare_DWBA_Zmode(0.0, 0.0, 1.0, 0.05, phi=0.01)
    # Out of the Ewald sphere for this cell and energy.
    with pytest.raises(ValueError, match="no solution"):
        crystal.prepare_DWBA_Zmode(0.0, 0.0, 8.0, 0.0, fixed="eq")
    # Reachable glancing angles, but no detector solution for delta.
    with pytest.raises(ValueError, match="no solution"):
        crystal.prepare_DWBA_Zmode(1.0, 0.0, 2.8, 0.05, fixed="in")
    # The exit wave would leave along the surface or below it.
    with pytest.raises(ValueError, match=r"\(0, pi/2\]"):
        crystal.prepare_DWBA_Zmode(1.0, 0.0, 0.1, 0.05, fixed="in")
    with pytest.raises(NotImplementedError, match="only the bulk"):
        CTRcalc.SXRDCrystal(cell, _one_atom_cell()).prepare_DWBA_Zmode(
            0.0, 0.0, 1.0, 0.0, fixed="eq"
        )


def test_l_from_glancing_angles_satisfies_the_ewald_condition():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    h, k = 1.0, 0.4
    alpha_i = 0.05
    alpha_f = np.array([0.02, 0.09, 0.21])

    l_value = crystal.l_from_glancing_angles(h, k, alpha_i, alpha_f)
    assert l_value.shape == alpha_f.shape

    k0 = 2.0 * np.pi / (HC_KEV_ANGSTROM / (cell._E * 1e-3))
    hkl = cell.refHKLTransform @ np.vstack(
        (np.full_like(l_value, h), np.full_like(l_value, k), l_value)
    )
    np.testing.assert_allclose(
        (cell.B_mat @ hkl)[2],
        k0 * (np.sin(alpha_i) + np.sin(alpha_f)),
        rtol=1e-14,
    )
    assert isinstance(
        crystal.l_from_glancing_angles(h, k, alpha_i, 0.09), float
    )


def test_prepare_DWBA_from_angles_matches_an_explicit_l():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    h, k = 1.0, 0.0
    alpha_i = 0.05
    alpha_f = np.array([0.03, 0.11, 0.19])

    prepared, l_value = crystal.prepare_DWBA_from_angles(
        h, k, alpha_i, alpha_f, return_l=True
    )
    np.testing.assert_allclose(
        l_value, crystal.l_from_glancing_angles(h, k, alpha_i, alpha_f)
    )
    np.testing.assert_allclose(prepared.alpha_i, alpha_i)
    np.testing.assert_allclose(prepared.alpha_f, alpha_f)
    np.testing.assert_array_equal(
        crystal.F_DWBA_prepared(prepared),
        crystal.F_DWBA(h, k, l_value, alpha_i, alpha_f),
    )

    scalar = crystal.prepare_DWBA_from_angles(h, k, alpha_i, 0.11)
    assert scalar.scalar
    assert isinstance(crystal.F_DWBA_prepared(scalar), complex)


def test_prepare_DWBA_from_angles_inverts_the_zmode_constraint():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    h, k = 1.0, 0.0
    L = np.array([0.8, 1.0, 1.2])
    fixedangle = 0.05

    prepared_z, angles = crystal.prepare_DWBA_Zmode(
        h, k, L, fixedangle, fixed="in", return_angles=True
    )
    prepared_a, l_value = crystal.prepare_DWBA_from_angles(
        h, k, angles[:, 0], angles[:, 2], return_l=True
    )
    np.testing.assert_allclose(l_value, L, rtol=1e-13)
    np.testing.assert_array_equal(
        crystal.F_DWBA_prepared(prepared_a),
        crystal.F_DWBA_prepared(prepared_z),
    )


def test_from_angles_rejects_a_degenerate_reference_axis_and_bad_input():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    with pytest.raises(ValueError, match="finite and nonempty"):
        crystal.l_from_glancing_angles(0.0, 0.0, np.nan, 0.05)
    with pytest.raises(ValueError, match=r"\(0, pi/2\]"):
        crystal.prepare_DWBA_from_angles(0.0, 0.0, 0.05, -0.05)

    # An out-of-plane reference axis lying in the surface plane leaves Q_z
    # independent of l. SXRDCrystal resets refHKLTransform, so reindex after
    # constructing it.
    crystal.uc_bulk.refHKLTransform = np.array(
        [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]
    )
    with pytest.raises(ValueError, match="no surface-normal component"):
        crystal.l_from_glancing_angles(0.0, 0.0, 0.05, 0.05)
    with pytest.raises(ValueError, match="no surface-normal component"):
        crystal.prepare_DWBA_from_angles(0.0, 0.0, 0.05, 0.05)

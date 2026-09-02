"""Focused validation of the public DWBA state and native bulk kernel."""

from dataclasses import FrozenInstanceError
from unittest import mock

import numpy as np
import pytest

from .. import CTRcalc, CTRdwba, CTRuc
from ..CTRoptics import HC_KEV_ANGSTROM


def _one_atom_cell(occupancy=1.0):
    cell = CTRuc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0])
    cell.addAtom("C", [0.17, 0.23, 0.31], 0.18, 0.37, occupancy)
    cell.setEnergy(10000.0)
    cell.f[0, 11] = 0.27
    cell.f[0, 12] = 0.19
    return cell


def _kinematics(cell, alpha_i, alpha_f, azimuth):
    wavelength = HC_KEV_ANGSTROM / (cell._E * 1e-3)
    k0 = 2.0 * np.pi / wavelength
    ki = k0 * np.cos(alpha_i)
    kf = k0 * np.cos(alpha_f)
    q_parallel = np.sqrt(ki**2 + kf**2 - 2.0 * ki * kf * np.cos(azimuth))
    qz = k0 * (np.sin(alpha_i) + np.sin(alpha_f))
    h = q_parallel * cell.a[0] / (2.0 * np.pi)
    l_value = qz * cell.a[2] / (2.0 * np.pi)
    return h, l_value


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
    alpha_i, alpha_f, azimuth = 0.21, 0.34, 0.43
    h, l_value = _kinematics(cell, alpha_i, alpha_f, azimuth)
    prepared = crystal.dwba.prepare_from_glancing(
        h,
        0.0,
        alpha_i,
        alpha_f,
        polarization_i=polarization_i,
        polarization_f=polarization_f,
    )
    np.testing.assert_array_equal(prepared.field_i.n, np.ones(2))
    cell.basis[0, 6] = 1.0
    actual = crystal.dwba.evaluate_prepared(
        prepared, bulk_mode="unit_cell"
    ).F_contrast
    expected = (
        cell.F_uc(np.array([h]), np.array([0.0]), np.array([l_value]))[0]
        * contraction(alpha_i, alpha_f, azimuth)
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)
    assert isinstance(actual, complex)


@pytest.mark.parametrize(
    ("polarization_i", "polarization_f"),
    [("s", "s"), ("s", "p"), ("p", "s"), ("p", "p")],
)
def test_native_four_channel_coefficients_use_renaud_vectors(
    polarization_i, polarization_f
):
    alpha_i, alpha_f, azimuth, k0 = 0.23, 0.31, 0.41, 5.3
    kz_i = np.array([[-k0 * np.sin(alpha_i)]], dtype=np.complex128)
    kz_f = np.array([[-k0 * np.sin(alpha_f)]], dtype=np.complex128)
    incident = {1: 0.73 + 0.11j, -1: -0.19 + 0.07j}
    exit_ = {1: 0.81 - 0.05j, -1: 0.16 + 0.09j}
    Aip = np.array([[incident[1]]])
    Aim = np.array([[incident[-1]]])
    Afp = np.array([[exit_[1]]])
    Afm = np.array([[exit_[-1]]])
    factors = np.zeros((1, 13))
    factors[0, 10] = 1.0
    basis = np.array([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
    domains = np.zeros((1, 3, 4))
    domains[0, :, :3] = np.identity(3)
    actual = CTRuc._CTRcalc_cpp.unitcell_F_DWBA_bulk(
        np.array([[0.17, 0.09]]),
        np.array([0], dtype=np.int64),
        np.array([1.0 + 0.0j]),
        np.array([0], dtype=np.int64),
        np.array([0], dtype=np.int64),
        kz_i,
        Aip,
        Aim,
        Aip / kz_i,
        Aim / kz_i,
        kz_f,
        Afp,
        Afm,
        Afp / kz_f,
        Afm / kz_f,
        np.array([0.0]),
        np.array([0.0]),
        np.array([alpha_i]),
        np.array([alpha_f]),
        np.array([np.cos(azimuth)]),
        np.array([np.sin(azimuth)]),
        k0,
        polarization_i,
        polarization_f,
        np.array([0], dtype=np.int64),
        False,
        0.0,
        basis,
        factors,
        np.identity(3),
        np.identity(3),
        domains,
        np.array([1.0]),
    )[0]
    expected = 0.0j
    for sigma_i in (1, -1):
        for sigma_f in (1, -1):
            if polarization_i == polarization_f == "s":
                factor = np.cos(azimuth)
            elif polarization_i == "s":
                factor = -sigma_f * np.sin(alpha_f) * np.sin(azimuth)
            elif polarization_f == "s":
                factor = -sigma_i * np.sin(alpha_i) * np.sin(azimuth)
            else:
                factor = (
                    np.cos(alpha_i) * np.cos(alpha_f)
                    - sigma_i
                    * sigma_f
                    * np.sin(alpha_i)
                    * np.sin(alpha_f)
                    * np.cos(azimuth)
                )
            expected += incident[sigma_i] * exit_[sigma_f] * factor
    np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=2e-15)


def test_lazy_state_and_legacy_api_removal():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    serialized = crystal.toStr()
    assert "_dwba_state" not in crystal.__dict__
    assert crystal.dwba is crystal.dwba
    assert isinstance(crystal.dwba, CTRdwba.DWBAState)
    assert crystal.toStr() == serialized
    for name in (
        "prepare_DWBA",
        "prepare_DWBA_Zmode",
        "F_DWBA",
        "F_DWBA_prepared",
        "specular_DWBA_reflectivity",
    ):
        assert not hasattr(crystal, name)


def test_fixed_incident_exit_equal_and_automatic_specular_geometry():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    crystal.dwba.set_ctr_geometry(alpha_i=0.05)
    incident = crystal.dwba.prepare(1.0, 0.0, 0.8)
    np.testing.assert_allclose(incident.alpha_i, 0.05, atol=1e-15)

    crystal.dwba.set_ctr_geometry(alpha_f=0.09, rods=[(2, 0)])
    exit_ = crystal.dwba.prepare(2.0, 0.0, 0.8)
    np.testing.assert_allclose(exit_.alpha_f, 0.09, atol=1e-15)

    specular = crystal.dwba.prepare(0.0, 0.0, np.array([0.4, 0.8]))
    np.testing.assert_allclose(specular.alpha_i, specular.alpha_f)
    assert np.all(specular.is_specular)

    crystal.dwba.set_ctr_geometry(equal_angles=True, rods=[(1, 0)])
    equal = crystal.dwba.prepare(1.0, 0.0, 0.8)
    np.testing.assert_allclose(equal.alpha_i, equal.alpha_f)


def test_geometry_rule_validation_and_conflicting_mixed_rods():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    with pytest.raises(ValueError, match="exactly one"):
        crystal.dwba.set_ctr_geometry(alpha_i=0.05, alpha_f=0.06)
    with pytest.raises(ValueError, match="exactly one"):
        crystal.dwba.set_ctr_geometry()
    crystal.dwba.set_ctr_geometry(alpha_i=0.05)
    crystal.dwba.set_ctr_geometry(alpha_f=0.08, rods=[(2, 0)])
    with pytest.raises(ValueError, match="split by rod"):
        crystal.dwba.prepare([1.0, 2.0], 0.0, [0.8, 0.8])


def test_global_mixed_geometry_groups_specular_and_offspecular_points():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    crystal.dwba.set_ctr_geometry(alpha_i=0.05)
    prepared = crystal.dwba.prepare([0.0, 1.0], 0.0, [0.6, 0.8])
    assert prepared.is_specular.tolist() == [1, 0]
    assert prepared.alpha_i[0] == pytest.approx(prepared.alpha_f[0])
    assert prepared.alpha_i[1] == pytest.approx(0.05)
    result = crystal.dwba.evaluate_prepared(prepared)
    assert result.unperturbed_amplitude[0] != 0.0
    assert result.unperturbed_amplitude[1] == 0.0
    with pytest.raises(ValueError, match="only for specular"):
        _ = result.reflectivity


def test_glancing_inversion_and_vlieg_round_trip_with_orientation():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    theta = 0.17
    U = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    crystal.dwba.set_orientation(U)
    crystal.dwba.set_ctr_geometry(alpha_i=0.05)
    prepared = crystal.dwba.prepare(1.0, 0.0, np.array([0.7, 1.0]))
    measured = np.moveaxis(prepared.vlieg_angles, -1, 0)
    round_trip = crystal.dwba.prepare_from_vlieg(*measured)
    np.testing.assert_allclose(round_trip.hkl, prepared.hkl, atol=2e-14)
    np.testing.assert_allclose(round_trip.alpha_i, prepared.alpha_i)
    np.testing.assert_allclose(round_trip.alpha_f, prepared.alpha_f)

    direct = crystal.dwba.prepare_from_glancing(
        1.0, 0.0, prepared.alpha_i, prepared.alpha_f
    )
    np.testing.assert_allclose(direct.hkl[2], prepared.hkl[2], rtol=1e-13)


def test_orientation_validation_and_stale_preparation_rejection():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    with pytest.raises(ValueError, match="3-by-3"):
        crystal.dwba.set_orientation(np.ones((2, 2)))
    with pytest.raises(ValueError, match="nonsingular"):
        crystal.dwba.set_orientation(np.zeros((3, 3)))
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    crystal.dwba.set_orientation(np.diag([-1.0, -1.0, 1.0]))
    with pytest.raises(ValueError, match="orientation changed"):
        crystal.dwba.evaluate_prepared(prepared)


def test_unique_angle_tables_and_preparation_cache_avoid_field_solves():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    alpha_i = np.full(3, 0.05)
    alpha_f = np.array([0.08, 0.10, 0.08])
    h = np.array(
        [_kinematics(cell, ai, af, 0.25)[0] for ai, af in zip(alpha_i, alpha_f)]
    )
    with mock.patch.object(
        CTRdwba,
        "solve_electric_field",
        wraps=CTRdwba.solve_electric_field,
    ) as solve:
        prepared = crystal.dwba.prepare_from_glancing(h, 0.0, alpha_i, alpha_f)
        assert solve.call_count == 2
        assert prepared.field_i.angle_shape == (1,)
        assert prepared.field_f.angle_shape == (2,)
        np.testing.assert_array_equal(prepared.field_i_index, [0, 0, 0])
        np.testing.assert_array_equal(prepared.field_f_index, [0, 1, 0])
        repeated = crystal.dwba.prepare_from_glancing(h, 0.0, alpha_i, alpha_f)
        assert repeated is prepared
        crystal.dwba.evaluate_prepared(prepared)
        assert solve.call_count == 2
    info = crystal.dwba.cache_info()
    assert info["prepared_hits"] == 1
    assert info["field_misses"] == 3


def test_atomic_changes_reuse_fields_but_live_reference_changes_rebuild_them():
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    info = crystal.dwba.cache_info()
    cell.basis[0, 4] *= 1.5
    changed_dw = crystal.dwba.evaluate_prepared(prepared)
    assert np.isfinite(changed_dw.structure_factor_squared)
    assert crystal.dwba.cache_info()["field_misses"] == info["field_misses"]

    frozen_n = prepared.reference.n.copy()
    cell.basis[0, 6] *= 0.8
    live = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    assert not np.array_equal(live.reference.n, frozen_n)
    np.testing.assert_array_equal(prepared.reference.n, frozen_n)
    frozen_result = crystal.dwba.evaluate_prepared(prepared)
    assert np.isfinite(frozen_result.structure_factor_squared)


def test_dwba_wavefield_and_reflectivity_share_the_exact_bulk_reference():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    public_field = crystal.wavefield(0.04)
    np.testing.assert_array_equal(prepared.reference.n, public_field.n)
    np.testing.assert_array_equal(prepared.field_i.kz[:, 0], public_field.kz)
    np.testing.assert_allclose(
        np.abs(prepared.field_i.r_S[0]) ** 2,
        crystal.specular_reflectivity(0.04),
        rtol=0.0,
        atol=0.0,
    )


def test_result_observables_and_unpolarized_reflectivity():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    angles = np.array([0.001, 0.006, 0.04])
    result = crystal.dwba.evaluate_from_glancing(0.0, 0.0, angles, angles)
    np.testing.assert_allclose(
        result.total_amplitude,
        result.unperturbed_amplitude + result.scattered_amplitude,
    )
    np.testing.assert_allclose(
        result.structure_factor_squared, np.abs(result.F_contrast) ** 2
    )
    np.testing.assert_allclose(
        result.scattered_amplitude_squared,
        np.abs(result.scattered_amplitude) ** 2,
    )
    np.testing.assert_allclose(
        result.differential_cross_section_kernel,
        2.8179403262e-5**2 * np.abs(result.F_contrast) ** 2,
    )
    np.testing.assert_allclose(result.reflectivity, np.abs(result.total_amplitude) ** 2)
    expected_first = np.abs(result.unperturbed_amplitude) ** 2 + 2.0 * np.real(
        np.conj(result.unperturbed_amplitude) * result.scattered_amplitude
    )
    np.testing.assert_allclose(result.first_order_reflectivity, expected_first)

    s = crystal.dwba.reflectivity(0.0, 0.0, result.prepared.hkl[2], polarization="s")
    p = crystal.dwba.reflectivity(0.0, 0.0, result.prepared.hkl[2], polarization="p")
    unpolarized = crystal.dwba.reflectivity(
        0.0, 0.0, result.prepared.hkl[2], polarization="unpolarized"
    )
    np.testing.assert_allclose(unpolarized, 0.5 * (s + p))


def test_scalar_results_and_immutable_preparations():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    result = crystal.dwba.evaluate_prepared(prepared)
    assert isinstance(result.F_contrast, complex)
    assert isinstance(result.reflectivity, float)
    with pytest.raises(ValueError, match="read-only"):
        prepared.alpha_i[0] = 0.2
    with pytest.raises(FrozenInstanceError):
        result.bulk_mode = "unit_cell"


def test_zero_attenuation_pole_and_unit_cell_reflectivity_rejection():
    cell = _one_atom_cell(occupancy=0.0)
    crystal = CTRcalc.SXRDCrystal(cell)
    k0 = 2.0 * np.pi / (HC_KEV_ANGSTROM / (cell._E * 1e-3))
    alpha = np.arcsin(np.pi / (cell.a[2] * k0))
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, alpha, alpha)
    cell.basis[0, 6] = 1.0
    with pytest.raises(ValueError, match="exact Bragg pole"):
        crystal.dwba.evaluate_prepared(prepared, bulk_attenuation=0.0)
    unit = crystal.dwba.evaluate_prepared(prepared, bulk_mode="unit_cell")
    with pytest.raises(ValueError, match="semi-infinite"):
        _ = unit.reflectivity


def test_energy_and_lattice_changes_reject_retained_preparation():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    crystal.setEnergy(11000.0)
    with pytest.raises(ValueError, match="energy changed"):
        crystal.dwba.evaluate_prepared(prepared)

    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    crystal.uc_bulk.B_mat = crystal.uc_bulk.B_mat.copy()
    crystal.uc_bulk.B_mat[0, 0] *= 1.001
    with pytest.raises(ValueError, match="B_mat changed"):
        crystal.dwba.evaluate_prepared(prepared)


def test_bulk_scope_and_missing_native_extension_errors():
    bulk = _one_atom_cell()
    surface = _one_atom_cell()
    with pytest.raises(NotImplementedError, match="only the bulk"):
        CTRcalc.SXRDCrystal(bulk, surface).dwba.prepare_from_glancing(
            0.0, 0.0, 0.04, 0.04
        )
    crystal = CTRcalc.SXRDCrystal(bulk)
    with mock.patch.object(CTRdwba, "HAS_CPP_ACCEL", False):
        with pytest.raises(RuntimeError, match="native CTR extension"):
            crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)


def test_cache_clear_and_info_contract():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    before = crystal.dwba.cache_info()
    assert before["prepared_size"] == 1
    assert before["reference_size"] == 1
    crystal.dwba.clear_cache()
    after = crystal.dwba.cache_info()
    assert after["prepared_size"] == 0
    assert after["reference_size"] == 0
    assert after["field_size"] == 0


def test_preparation_lru_is_bounded_but_retained_handles_remain_usable():
    crystal = CTRcalc.SXRDCrystal(_one_atom_cell())
    retained = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.02, 0.02)
    for angle in np.linspace(0.021, 0.08, 36):
        crystal.dwba.prepare_from_glancing(0.0, 0.0, angle, angle)
    info = crystal.dwba.cache_info()
    assert info["prepared_size"] == info["prepared_capacity"] == 32
    result = crystal.dwba.evaluate_prepared(retained)
    assert np.isfinite(result.reflectivity)

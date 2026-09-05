"""Focused validation of the public DWBA state and native record kernel."""

from dataclasses import FrozenInstanceError
from unittest import mock

import numpy as np
import pytest

from .. import CTRcalc, CTRdwba, CTRfilm, CTRuc
from ..CTRfilm import PoissonProfile, SkellamProfile
from ..CTRoptics import HC_KEV_ANGSTROM


def _one_atom_cell(occupancy=1.0):
    cell = CTRuc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0])
    cell.addAtom("C", [0.17, 0.23, 0.31], 0.18, 0.37, occupancy)
    cell.setEnergy(10000.0)
    cell.f[0, 11] = 0.27
    cell.f[0, 12] = 0.19
    return cell


def _layered_cell(name="layered"):
    cell = CTRuc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0], name=name)
    cell.addAtom("C", [0.13, 0.17, 0.0], 0.11, 0.19, 1.0, layer=1)
    cell.addAtom("O", [0.31, 0.29, 0.5], 0.14, 0.23, 1.0, layer=2)
    cell.layerpos = {1.0: 0.0, 2.0: 0.5}
    cell.setEnergy(10000.0)
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


def _synthetic_record_call(
    *,
    q_parallel,
    kz_i,
    kz_f,
    A_i_plus,
    A_i_minus,
    A_f_plus,
    A_f_minus,
    alpha_i=0.23,
    alpha_f=0.31,
    cos_azimuth=np.cos(0.41),
    sin_azimuth=np.sin(0.41),
    polarization_i="s",
    polarization_f="s",
    is_specular=False,
    interfaces=(0.0,),
    reference_shares=None,
    basis=None,
    factors=None,
    finite_positions=None,
    finite_atoms=None,
    finite_records=None,
    finite_weights=None,
    reference_area=1.0,
):
    """Call the packed native kernel with one synthetic scan point."""
    kz_i = np.ascontiguousarray(kz_i, dtype=np.complex128).reshape(-1, 1)
    kz_f = np.ascontiguousarray(kz_f, dtype=np.complex128).reshape(-1, 1)
    media = len(kz_i)
    A_i_plus = np.ascontiguousarray(A_i_plus, dtype=np.complex128).reshape(
        media, 1
    )
    A_i_minus = np.ascontiguousarray(A_i_minus, dtype=np.complex128).reshape(
        media, 1
    )
    A_f_plus = np.ascontiguousarray(A_f_plus, dtype=np.complex128).reshape(
        media, 1
    )
    A_f_minus = np.ascontiguousarray(A_f_minus, dtype=np.complex128).reshape(
        media, 1
    )
    if reference_shares is None:
        reference_shares = np.zeros((2, media, 2), dtype=np.float64)
    if basis is None:
        basis = np.empty((0, 7), dtype=np.float64)
    if factors is None:
        factors = np.empty((0, 13), dtype=np.float64)
    if finite_positions is None:
        finite_positions = np.empty((0, 3), dtype=np.float64)
    count = len(finite_positions)
    if finite_atoms is None:
        finite_atoms = np.zeros(count, dtype=np.int64)
    if finite_records is None:
        finite_records = np.empty((0,), dtype=np.int64)
    if finite_weights is None:
        finite_weights = np.ones(count, dtype=np.float64)
    return CTRuc._CTRcalc_cpp.unitcell_F_DWBA_records(
        np.asarray([q_parallel], dtype=np.float64),
        np.asarray([is_specular], dtype=np.int64),
        np.ones(media, dtype=np.complex128),
        np.zeros(1, dtype=np.int64),
        np.zeros(1, dtype=np.int64),
        kz_i,
        A_i_plus,
        A_i_minus,
        A_i_plus / kz_i,
        A_i_minus / kz_i,
        kz_f,
        A_f_plus,
        A_f_minus,
        A_f_plus / kz_f,
        A_f_minus / kz_f,
        np.zeros(media),
        np.zeros(media),
        np.asarray(interfaces, dtype=np.float64),
        np.asarray([alpha_i]),
        np.asarray([alpha_f]),
        np.asarray([cos_azimuth]),
        np.asarray([sin_azimuth]),
        5.3,
        polarization_i,
        polarization_f,
        np.ascontiguousarray(reference_shares, dtype=np.float64),
        reference_area,
        True,
        0.0,
        np.ascontiguousarray(basis, dtype=np.float64),
        np.ascontiguousarray(factors, dtype=np.float64),
        np.ascontiguousarray(finite_positions, dtype=np.float64).reshape(-1, 3),
        np.ascontiguousarray(finite_atoms, dtype=np.int64),
        np.ascontiguousarray(finite_records, dtype=np.int64),
        np.ascontiguousarray(finite_weights, dtype=np.float64),
        np.empty((0, 3), dtype=np.float64),
        np.empty((0,), dtype=np.int64),
        np.empty((0,), dtype=np.float64),
        np.empty((0,), dtype=np.float64),
        np.asarray([0.0, 0.0, 4.0]),
    )


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


@pytest.mark.parametrize(
    ("polarization_i", "polarization_f"),
    [("s", "s"), ("s", "p"), ("p", "s"), ("p", "p")],
)
def test_packed_record_kernel_matches_direct_four_channel_atom_sum(
    polarization_i, polarization_f
):
    q_parallel = np.array([0.37, 0.19])
    kz_i = np.array([-1.21 + 0.04j, -1.07 + 0.09j])
    kz_f = np.array([-1.46 + 0.03j, -1.18 + 0.07j])
    A_i_plus = np.array([0.81 + 0.06j, 0.73 - 0.04j])
    A_i_minus = np.array([-0.17 + 0.08j, 0.21 + 0.03j])
    A_f_plus = np.array([0.77 - 0.05j, 0.69 + 0.02j])
    A_f_minus = np.array([0.14 + 0.07j, -0.12 + 0.09j])
    basis = np.array([[6.0, 0.0, 0.0, 0.0, 0.23, 0.31, 0.67]])
    factors = np.zeros((1, 13))
    factors[0, 0] = 1.7
    factors[0, 5] = 0.42
    factors[0, 10:13] = [0.8, 0.13, 0.21]
    position = np.array([[0.27, -0.18, -0.43]])
    actual = _synthetic_record_call(
        q_parallel=q_parallel,
        kz_i=kz_i,
        kz_f=kz_f,
        A_i_plus=A_i_plus,
        A_i_minus=A_i_minus,
        A_f_plus=A_f_plus,
        A_f_minus=A_f_minus,
        polarization_i=polarization_i,
        polarization_f=polarization_f,
        basis=basis,
        factors=factors,
        finite_positions=position,
        finite_records=np.array([1]),
        finite_weights=np.array([-0.72]),
    )[0][1, 0]

    medium = 1
    alpha_i, alpha_f, azimuth, k0 = 0.23, 0.31, 0.41, 5.3
    expected = 0.0j
    for sigma_i, incident, incident_over_kz in (
        (1, A_i_plus[medium], A_i_plus[medium] / kz_i[medium]),
        (-1, A_i_minus[medium], A_i_minus[medium] / kz_i[medium]),
    ):
        for sigma_f, exit_, exit_over_kz in (
            (1, A_f_plus[medium], A_f_plus[medium] / kz_f[medium]),
            (-1, A_f_minus[medium], A_f_minus[medium] / kz_f[medium]),
        ):
            qz = -(sigma_i * kz_i[medium] + sigma_f * kz_f[medium])
            q_squared = q_parallel @ q_parallel + qz**2
            atom_factor = (
                1.7 * np.exp(-0.42 * q_squared) + 0.8 + 0.13 + 0.21j
            )
            atom_factor *= np.exp(
                -(0.23 * (q_parallel @ q_parallel) + 0.31 * qz**2)
                / (16.0 * np.pi**2)
            )
            atom_factor *= 0.67
            if polarization_i == polarization_f == "s":
                vector = incident * exit_ * np.cos(azimuth)
            elif polarization_i == "s":
                vector = (
                    -sigma_f
                    * incident
                    * exit_
                    * np.sin(alpha_f)
                    * np.sin(azimuth)
                )
            elif polarization_f == "s":
                vector = (
                    -sigma_i
                    * incident
                    * exit_
                    * np.sin(alpha_i)
                    * np.sin(azimuth)
                )
            else:
                tangent = (
                    -sigma_i
                    * sigma_f
                    * incident
                    * exit_
                    * np.sin(alpha_i)
                    * np.sin(alpha_f)
                    * np.cos(azimuth)
                )
                incident_normal = (
                    -k0
                    * np.cos(alpha_i)
                    * np.sin(alpha_i)
                    * incident_over_kz
                )
                exit_normal = (
                    -k0
                    * np.cos(alpha_f)
                    * np.sin(alpha_f)
                    * exit_over_kz
                )
                vector = tangent + incident_normal * exit_normal
            expected += (
                -0.72
                * vector
                * atom_factor
                * np.exp(1j * (q_parallel @ position[0, :2] + qz * position[0, 2]))
            )
    np.testing.assert_allclose(actual, expected, rtol=3e-14, atol=3e-14)


def test_packed_record_kernel_keeps_distinct_debye_waller_states():
    common = {
        "q_parallel": np.array([0.37, 0.19]),
        "kz_i": np.array([-1.21 + 0.04j, -1.07 + 0.09j]),
        "kz_f": np.array([-1.46 + 0.03j, -1.18 + 0.07j]),
        "A_i_plus": np.array([0.81 + 0.06j, 0.73 - 0.04j]),
        "A_i_minus": np.array([-0.17 + 0.08j, 0.21 + 0.03j]),
        "A_f_plus": np.array([0.77 - 0.05j, 0.69 + 0.02j]),
        "A_f_minus": np.array([0.14 + 0.07j, -0.12 + 0.09j]),
    }
    basis = np.array(
        [
            [6.0, 0.0, 0.0, 0.0, 0.08, 0.21, 0.67],
            [6.0, 0.0, 0.0, 0.0, 0.39, 0.57, -0.42],
        ]
    )
    factor = np.zeros(13)
    factor[0] = 1.7
    factor[5] = 0.42
    factor[10:13] = [0.8, 0.13, 0.21]
    factors = np.repeat(factor[None, :], 2, axis=0)
    positions = np.array(
        [
            [0.27, -0.18, -0.43],
            [-0.11, 0.22, -0.67],
            [0.08, 0.31, -0.91],
            [-0.24, -0.07, -1.12],
        ]
    )
    atom_indices = np.array([0, 1, 0, 1])
    weights = np.array([-0.72, 0.31, 0.44, -0.19])

    combined = _synthetic_record_call(
        **common,
        basis=basis,
        factors=factors,
        finite_positions=positions,
        finite_atoms=atom_indices,
        finite_records=np.ones(4, dtype=np.int64),
        finite_weights=weights,
    )[0][1, 0]
    separate = 0.0j
    for position, atom, weight in zip(positions, atom_indices, weights):
        separate += _synthetic_record_call(
            **common,
            basis=basis[atom : atom + 1],
            factors=factors[atom : atom + 1],
            finite_positions=position[None, :],
            finite_atoms=np.zeros(1, dtype=np.int64),
            finite_records=np.ones(1, dtype=np.int64),
            finite_weights=np.array([weight]),
        )[0][1, 0]
    np.testing.assert_allclose(combined, separate, rtol=3e-14, atol=3e-14)


def test_packed_record_medium_assignment_is_half_open_at_an_interface():
    q_parallel = np.array([0.17, 0.08])
    kz_i = np.array([-1.0 + 0.02j, -1.3 + 0.05j])
    kz_f = np.array([-1.1 + 0.03j, -1.4 + 0.06j])
    basis = np.array([[6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
    factors = np.zeros((1, 13))
    factors[0, 10] = 1.0
    positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -1e-12]])
    atomic = _synthetic_record_call(
        q_parallel=q_parallel,
        kz_i=kz_i,
        kz_f=kz_f,
        A_i_plus=np.array([1.0, 2.0]),
        A_i_minus=np.zeros(2),
        A_f_plus=np.array([1.0, 3.0]),
        A_f_minus=np.zeros(2),
        basis=basis,
        factors=factors,
        finite_positions=positions,
        finite_records=np.array([1, 2]),
        reference_shares=np.zeros((3, 2, 2)),
    )[0]
    qz_upper = -(kz_i[0] + kz_f[0])
    qz_lower = -(kz_i[1] + kz_f[1])
    np.testing.assert_allclose(atomic[1, 0], np.cos(0.41))
    np.testing.assert_allclose(
        atomic[2, 0],
        6.0 * np.cos(0.41) * np.exp(1j * qz_lower * -1e-12),
    )
    assert qz_upper != qz_lower


@pytest.mark.parametrize("q", [0.73 - 0.19j, 1e-12 + 0.0j])
def test_packed_record_finite_reference_uses_stable_phi(q):
    kzi = np.full(3, -0.5 * q, dtype=np.complex128)
    shares = np.zeros((2, 3, 2))
    shares[1, 1, 0] = (
        2.0 * np.pi * 2.8179403262e-5 / 5.3**2
    )
    reference = _synthetic_record_call(
        q_parallel=np.zeros(2),
        kz_i=kzi,
        kz_f=kzi,
        A_i_plus=np.ones(3),
        A_i_minus=np.zeros(3),
        A_f_plus=np.ones(3),
        A_f_minus=np.zeros(3),
        cos_azimuth=1.0,
        sin_azimuth=0.0,
        is_specular=True,
        interfaces=(2.0, 0.0),
        reference_shares=shares,
    )[1]
    expected = np.exp(1j * q * 0.0) * np.expm1(1j * q * 2.0) / (1j * q)
    np.testing.assert_allclose(reference[1, 0], expected, rtol=3e-14, atol=3e-14)


def test_packed_record_terminal_reference_uses_half_space_integral():
    q = 0.73 - 0.19j
    kzi = np.full(3, -0.5 * q, dtype=np.complex128)
    shares = np.zeros((2, 3, 2))
    shares[1, 2, 0] = (
        2.0 * np.pi * 2.8179403262e-5 / 5.3**2
    )
    reference = _synthetic_record_call(
        q_parallel=np.zeros(2),
        kz_i=kzi,
        kz_f=kzi,
        A_i_plus=np.ones(3),
        A_i_minus=np.zeros(3),
        A_f_plus=np.ones(3),
        A_f_minus=np.zeros(3),
        cos_azimuth=1.0,
        sin_azimuth=0.0,
        is_specular=True,
        interfaces=(2.0, 0.0),
        reference_shares=shares,
    )[1]
    np.testing.assert_allclose(reference[1, 0], 1.0 / (1j * q), rtol=3e-14)


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


@pytest.mark.parametrize(
    ("bulk_mode", "attenuation"),
    [("unit_cell", 0.0), ("semi_infinite", 0.073)],
)
def test_surface_free_packed_kernel_preserves_the_bulk_kernel(
    bulk_mode, attenuation
):
    cell = _one_atom_cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    alpha_i = np.array([0.04, 0.055])
    alpha_f = np.array([0.04, 0.08])
    h = np.array(
        [_kinematics(cell, ai, af, 0.27)[0] for ai, af in zip(alpha_i, alpha_f)]
    )
    h[0] = 0.0
    prepared = crystal.dwba.prepare_from_glancing(h, 0.0, alpha_i, alpha_f)
    actual = crystal.dwba.evaluate_prepared(
        prepared,
        bulk_mode=bulk_mode,
        bulk_attenuation=attenuation,
    ).F_contrast
    basis, factors, _ = cell.build_selected_basis()
    domains, occupancies = CTRuc._coherent_domain_arrays(
        cell.coherentDomainMatrix,
        cell.coherentDomainOccupancy,
    )
    media = len(prepared.reference.n)

    def field_array(value):
        return np.ascontiguousarray(np.asarray(value).reshape(media, -1))

    expected = CTRuc._CTRcalc_cpp.unitcell_F_DWBA_bulk(
        prepared.Q_parallel,
        prepared.is_specular,
        prepared.reference.n,
        prepared.field_i_index,
        prepared.field_f_index,
        field_array(prepared.field_i.kz),
        field_array(prepared.field_i.A_plus),
        field_array(prepared.field_i.A_minus),
        field_array(prepared.field_i._A_plus_over_kz),
        field_array(prepared.field_i._A_minus_over_kz),
        field_array(prepared.field_f.kz),
        field_array(prepared.field_f.A_plus),
        field_array(prepared.field_f.A_minus),
        field_array(prepared.field_f._A_plus_over_kz),
        field_array(prepared.field_f._A_minus_over_kz),
        prepared.field_i.z_reference,
        prepared.field_f.z_reference,
        prepared.alpha_i,
        prepared.alpha_f,
        prepared.cos_azimuth,
        prepared.sin_azimuth,
        prepared.k0,
        prepared.polarization_i,
        prepared.polarization_f,
        np.full(len(basis), media - 1, dtype=np.int64),
        bulk_mode == "semi_infinite",
        attenuation,
        np.ascontiguousarray(basis),
        np.ascontiguousarray(factors),
        np.ascontiguousarray(cell.R_mat),
        np.ascontiguousarray(cell.R_mat_inv),
        np.ascontiguousarray(domains),
        np.ascontiguousarray(occupancies),
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=3e-13)


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


def test_unit_cell_surface_support_and_missing_native_extension_errors():
    bulk = _one_atom_cell()
    surface = _one_atom_cell()
    surface.name = "surface"
    crystal = CTRcalc.SXRDCrystal(bulk, surface)
    result = crystal.dwba.evaluate_from_glancing(
        0.0, 0.0, 0.04, 0.04, bulk_mode="unit_cell"
    )
    assert [item.role for item in result.contributions] == [
        "bulk",
        "unit_cell_layer",
    ]
    assert result.contributions[1].component_name == "surface"
    assert result.F_reference == 0.0j

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


def test_contributions_sum_and_reference_is_specular_only():
    bulk = _layered_cell("bulk")
    film = CTRfilm.Film(_layered_cell("film_cell"), name="film")
    film.basis[0] = 3.0
    crystal = CTRcalc.SXRDCrystal(bulk, film, stacking=np.array([1]))
    crystal.dwba.set_ctr_geometry(alpha_i=0.05)
    prepared = crystal.dwba.prepare([0.0, 1.0], 0.0, [0.7, 0.8])
    result = crystal.dwba.evaluate_prepared(prepared)

    assert [item.role for item in result.contributions] == [
        "bulk",
        "film_layer",
        "film_layer",
    ]
    np.testing.assert_allclose(
        result.F_atomic,
        sum(item.F_atomic for item in result.contributions),
    )
    np.testing.assert_allclose(
        result.F_reference,
        sum(item.F_reference for item in result.contributions),
    )
    np.testing.assert_array_equal(
        result.F_contrast,
        sum(item.F_contrast for item in result.contributions),
    )
    assert result.F_reference[0] != 0.0j
    assert result.F_reference[1] == 0.0j
    for item in result.contributions:
        assert item.F_reference[1] == 0.0j
        with pytest.raises(ValueError, match="read-only"):
            item.F_atomic[0] = 0.0j

    diagnostic = crystal.dwba.evaluate_prepared(
        prepared, bulk_mode="unit_cell"
    )
    np.testing.assert_array_equal(diagnostic.F_reference, 0.0j)
    np.testing.assert_array_equal(diagnostic.unperturbed_amplitude, 0.0j)


def test_commensurate_surface_cells_preserve_reference_area_amplitude():
    bulk = _one_atom_cell(occupancy=0.0)
    primitive = CTRuc.UnitCell(
        [3.0, 3.0, 4.0], [90.0, 90.0, 90.0], name="primitive"
    )
    primitive.addAtom("C", [0.1, 0.2, 0.3], 0.18, 0.27, 1.0, layer=1)
    primitive.setEnergy(10000.0)
    supercell = CTRuc.UnitCell(
        [6.0, 3.0, 4.0], [90.0, 90.0, 90.0], name="supercell"
    )
    supercell.addAtom("C", [0.05, 0.2, 0.3], 0.18, 0.27, 1.0, layer=1)
    supercell.addAtom("C", [0.55, 0.2, 0.3], 0.18, 0.27, 1.0, layer=1)
    supercell.setEnergy(10000.0)

    amplitudes = []
    for surface in (primitive, supercell):
        crystal = CTRcalc.SXRDCrystal(bulk, surface)
        crystal.dwba.set_ctr_geometry(alpha_i=0.05)
        result = crystal.dwba.evaluate(
            1.0, 0.0, 0.8, bulk_mode="unit_cell"
        )
        amplitudes.append(result.contributions[1].F_atomic)
    np.testing.assert_allclose(amplitudes[0], amplitudes[1], rtol=3e-13)


@pytest.mark.parametrize(
    "model_factory, expected_roles",
    [
        (
            lambda: CTRfilm.EpitaxyInterface(
                _layered_cell("top"),
                _layered_cell("bottom"),
                profile=SkellamProfile(0.35, 0.0),
                name="interface",
            ),
            {"interface_top", "interface_bottom"},
        ),
        (
            lambda: CTRfilm.PoissonSurface(
                _layered_cell("termination"),
                profile=PoissonProfile(1.5, 0.5),
                name="surface",
            ),
            {
                "surface_termination",
                "covered_film",
                "sharp_film_correction",
            },
        ),
    ],
)
def test_generated_model_records_preserve_the_prepared_optical_profile(
    model_factory, expected_roles
):
    bulk = _layered_cell("bulk")
    model = model_factory()
    components = [model]
    if isinstance(model, CTRfilm.PoissonSurface):
        film = CTRfilm.Film(_layered_cell("film"), name="film")
        film.basis[0] = 3.0
        components.insert(0, film)
    crystal = CTRcalc.SXRDCrystal(
        bulk,
        *components,
        stacking=np.arange(1, len(components) + 1),
    )
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    shares = prepared._atomic_model.reference_shares

    np.testing.assert_allclose(
        shares.sum(axis=0),
        prepared.reference.profile.values[::-1, 1:],
        rtol=2e-14,
        atol=1e-18,
    )
    roles = {item[5] for item in prepared._atomic_model.descriptors}
    assert expected_roles <= roles
    result = crystal.dwba.evaluate_prepared(prepared)
    assert np.isfinite(result.structure_factor_squared)


def test_zero_width_epitaxy_bottom_correction_cancels_in_dwba():
    """A sharp interface adds and removes the same lower-material planes."""
    bulk = _layered_cell("bulk")
    interface = CTRfilm.EpitaxyInterface(
        _layered_cell("top"),
        _layered_cell("bottom"),
        profile=SkellamProfile(0.0, 0.0),
        fixed_ucs=2,
        name="interface",
    )
    crystal = CTRcalc.SXRDCrystal(
        bulk, interface, stacking=np.array([1])
    )

    # The two points exercise the sharp-interface correction at Q_parallel = 0
    # and Q_parallel != 0 under the same prepared optical reference.
    prepared = crystal.dwba.prepare_from_glancing(
        np.array([0.0, 0.6]),
        0.0,
        np.array([0.04, 0.04]),
        np.array([0.04, 0.05]),
    )
    result = crystal.dwba.evaluate_prepared(prepared)
    records = crystal.dwba._live_records()
    np.testing.assert_array_equal(prepared.is_specular, [1, 0])

    bottom_records = [
        (record, contribution, prepared._atomic_model.reference_shares[index])
        for index, (record, contribution) in enumerate(
            zip(records, result.contributions)
        )
        if record.role == "interface_bottom"
    ]
    assert len(bottom_records) == 2

    for record, contribution, reference_shares in bottom_records:
        occupancy = np.asarray(record.cell.coherentDomainOccupancy)
        positive = np.flatnonzero(occupancy > 0.0)
        negative = np.flatnonzero(occupancy < 0.0)
        np.testing.assert_array_equal(occupancy[positive], [1.0])
        np.testing.assert_array_equal(occupancy[negative], [-1.0])
        np.testing.assert_array_equal(
            record.cell.coherentDomainMatrix[positive[0]],
            record.cell.coherentDomainMatrix[negative[0]],
        )

        # The optical profile retains the same signed pair before the DWBA
        # stratification combines coincident planes into a zero reference share.
        optical_profile = record.cell.optical_profile()
        positive = optical_profile[optical_profile[:, 1] > 0.0]
        negative = optical_profile[optical_profile[:, 1] < 0.0]
        assert positive.shape == negative.shape == (1, 3)
        np.testing.assert_array_equal(positive[:, 0], negative[:, 0])
        np.testing.assert_allclose(positive[:, 1:], -negative[:, 1:])
        np.testing.assert_array_equal(reference_shares, 0.0)

        # Each record is zero at both the specular and off-specular point;
        # checking atomic and reference terms separately prevents contrast-only
        # cancellation from concealing a decomposition error.
        np.testing.assert_allclose(contribution.F_atomic, 0.0j, atol=2e-13)
        np.testing.assert_array_equal(contribution.F_reference, 0.0j)
        np.testing.assert_allclose(contribution.F_contrast, 0.0j, atol=2e-13)

    # Isolate the two signed halves under the same frozen wavefield. This
    # rules out a vacuous zero caused by the native kernel skipping the lower
    # interface records: each half scatters, and their amplitudes are opposite.
    original_occupancies = [
        record.cell.coherentDomainOccupancy.copy()
        for record, _, _ in bottom_records
    ]
    stacking_enabled = crystal.enable_uc_stacking
    try:
        # Prevent evaluate_prepared() from regenerating the interface and
        # restoring both signs while this diagnostic isolates one at a time.
        crystal.enable_uc_stacking = False
        for (record, _, _), occupancy in zip(
            bottom_records, original_occupancies
        ):
            record.cell.coherentDomainOccupancy = np.where(
                occupancy > 0.0, occupancy, 0.0
            )
        positive_result = crystal.dwba.evaluate_prepared(prepared)

        for (record, _, _), occupancy in zip(
            bottom_records, original_occupancies
        ):
            record.cell.coherentDomainOccupancy = np.where(
                occupancy < 0.0, occupancy, 0.0
            )
        negative_result = crystal.dwba.evaluate_prepared(prepared)
    finally:
        for (record, _, _), occupancy in zip(
            bottom_records, original_occupancies
        ):
            record.cell.coherentDomainOccupancy = occupancy
        crystal.enable_uc_stacking = stacking_enabled

    positive_atomic = [
        item.F_atomic
        for item in positive_result.contributions
        if item.role == "interface_bottom"
    ]
    negative_atomic = [
        item.F_atomic
        for item in negative_result.contributions
        if item.role == "interface_bottom"
    ]
    for positive, negative in zip(positive_atomic, negative_atomic):
        assert np.all(np.abs(positive) > 1e-12)
        np.testing.assert_allclose(positive, -negative, rtol=2e-13, atol=2e-13)


@pytest.mark.parametrize("mean_change", [-1.5, 0.0, 1.5])
def test_signed_poisson_records_are_coherent_and_zero_width_cancels(mean_change):
    bulk = _layered_cell("bulk")
    film = CTRfilm.Film(_layered_cell("film"), name="film")
    film.basis[0] = 4.0
    surface = CTRfilm.PoissonSurface(
        _layered_cell("termination"),
        profile=PoissonProfile(mean_change, 0.5),
        name="surface",
    )
    crystal = CTRcalc.SXRDCrystal(
        bulk, film, surface, stacking=np.array([1, 2])
    )
    result = crystal.dwba.evaluate_from_glancing(
        0.0, 0.0, 0.04, 0.04, bulk_mode="unit_cell"
    )
    surface_records = [
        item for item in result.contributions if item.component_name == "surface"
    ]
    assert surface_records
    np.testing.assert_allclose(
        sum(item.F_contrast for item in surface_records),
        sum(item.F_atomic for item in surface_records),
    )
    if mean_change == 0.0:
        np.testing.assert_allclose(
            sum(item.F_atomic for item in surface_records),
            0.0j,
            atol=2e-13,
        )


def test_retained_surface_preparation_allows_atomic_changes_but_not_topology():
    bulk = _layered_cell("bulk")
    film = CTRfilm.Film(_layered_cell("film"), name="film")
    film.basis[0] = 2.0
    crystal = CTRcalc.SXRDCrystal(bulk, film, stacking=np.array([1]))
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)

    film.layer_ucs[0].basis[0, 1] += 0.03
    live = crystal.dwba.evaluate_prepared(prepared)
    assert np.isfinite(live.F_contrast)
    live_film = sum(
        item.F_atomic for item in live.contributions if item.component_index == 0
    )
    frozen_reference = sum(
        item.F_reference for item in live.contributions if item.component_index == 0
    )

    crystal.weights[0] = 0.4
    reweighted = crystal.dwba.evaluate_prepared(prepared)
    np.testing.assert_allclose(
        sum(
            item.F_atomic
            for item in reweighted.contributions
            if item.component_index == 0
        ),
        0.4 * live_film,
    )
    np.testing.assert_array_equal(
        sum(
            item.F_reference
            for item in reweighted.contributions
            if item.component_index == 0
        ),
        frozen_reference,
    )
    rebuilt = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    assert rebuilt is not prepared
    assert not np.array_equal(rebuilt.reference.n, prepared.reference.n)

    film.basis[0] = 4.0
    with pytest.raises(ValueError, match="record geometry changed"):
        crystal.dwba.evaluate_prepared(prepared)


def test_water_special_factors_and_unsupported_transforms_fail_clearly(
    monkeypatch,
):
    bulk = _layered_cell("bulk")
    water = CTRuc.WaterModel(bulk.a, bulk.alpha, "step")
    water.setEnergy(10000.0)
    with pytest.raises(NotImplementedError, match="WaterModel"):
        CTRcalc.SXRDCrystal(bulk, water).dwba.prepare_from_glancing(
            0.0, 0.0, 0.04, 0.04
        )

    surface = _layered_cell("surface")
    surface.coherentDomainMatrix[0][0, 0] = 1.1
    crystal = CTRcalc.SXRDCrystal(bulk, surface)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    with pytest.raises(NotImplementedError, match="lateral-area"):
        crystal.dwba.evaluate_prepared(prepared)

    special = _one_atom_cell()
    monkeypatch.setitem(
        CTRuc.UnitCell.special_formfactors,
        "C",
        (lambda q: np.ones_like(q), lambda energy: (0.0, 0.0)),
    )
    special._test_special_formfactors()
    crystal = CTRcalc.SXRDCrystal(special)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.04, 0.04)
    with pytest.raises(NotImplementedError, match="special form-factor"):
        crystal.dwba.evaluate_prepared(prepared)

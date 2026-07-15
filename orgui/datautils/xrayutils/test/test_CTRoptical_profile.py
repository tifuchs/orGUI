import unittest

import numpy as np
import xraydb

from .. import CTRcalc, CTRfilm, CTRoptics, unitcells


AVOGADRO = 6.02214076e23


def density_from_cell(formula_mass, formula_units, volume):
    """Return density in g / cm**3 from an Angstrom**3 cell volume."""
    return formula_mass * formula_units / (volume * 1e-24 * AVOGADRO)


class TestUnitCellOpticalProfile(unittest.TestCase):
    energy_eV = 10000.0

    def test_requires_energy(self):
        cell = CTRcalc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0])
        cell.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)

        with self.assertRaisesRegex(ValueError, "Set the UnitCell energy"):
            cell.optical_profile()

    def test_pt100_matches_xraydb(self):
        cell = unitcells.unitcell("Pt100")
        cell.setEnergy(self.energy_eV)
        profile = cell.optical_profile()
        np.testing.assert_array_equal(profile, CTRoptics.optical_profile(cell))
        density = density_from_cell(
            xraydb.atomic_mass("Pt"), 4, cell.volume
        )
        delta, beta, _ = xraydb.xray_delta_beta("Pt", density, self.energy_eV)

        self.assertEqual(profile.shape, (1, 3))
        self.assertEqual(profile.dtype, np.float64)
        self.assertTrue(profile.flags.c_contiguous)
        self.assertEqual(profile[0, 0], 0.0)
        np.testing.assert_allclose(profile[:, 1], delta, rtol=1e-7)
        np.testing.assert_allclose(profile[:, 2], beta, rtol=1e-7)

    def test_tio2_100_has_two_homogeneous_layers(self):
        cell = unitcells.crystal("TiO2(100)").uc_bulk
        cell.setEnergy(self.energy_eV)
        profile = cell.optical_profile()
        formula_mass = xraydb.atomic_mass("Ti") + 2.0 * xraydb.atomic_mass("O")
        density = density_from_cell(formula_mass, 4, cell.volume)
        delta, beta, _ = xraydb.xray_delta_beta("TiO2", density, self.energy_eV)

        self.assertEqual(profile.shape, (2, 3))
        np.testing.assert_allclose(profile[:, 0], [-6.5807, -3.29035])
        np.testing.assert_allclose(profile[:, 1], delta, rtol=1e-7)
        np.testing.assert_allclose(profile[:, 2], beta, rtol=1e-7)

    def test_domain_translation_and_occupancy_weight_optical_contribution(self):
        cell = unitcells.unitcell("Pt100")
        cell.setEnergy(self.energy_eV)
        translated = np.eye(3, 4)
        translated[2, 3] = 0.25
        cell.coherentDomainMatrix = [np.eye(3, 4), translated]
        cell.coherentDomainOccupancy = [0.75, 0.25]

        profile = cell.optical_profile()

        np.testing.assert_allclose(profile[:, 0], [0.0, 0.25 * cell.a[2]])
        np.testing.assert_allclose(profile[1, 1], profile[0, 1] / 3.0)
        np.testing.assert_allclose(profile[1, 2], profile[0, 2] / 3.0)

    def test_film_and_crystal_combine_layer_profiles(self):
        bulk = unitcells.crystal("TiO2(100)").uc_bulk
        bulk.setEnergy(self.energy_eV)
        film = CTRfilm.Film(bulk)
        film.basis[0] = 2.0

        film_profile = film.optical_profile()
        crystal_profile = CTRcalc.SXRDCrystal(bulk, film).optical_profile()

        self.assertEqual(film_profile.shape, (2, 3))
        np.testing.assert_allclose(film_profile[:, 0], [0.0, 3.29035])
        self.assertEqual(crystal_profile.shape, (63, 3))
        np.testing.assert_allclose(
            crystal_profile[-5:, 0],
            [-6.5807, -3.29035, 0.0, 3.29035, 6.5807],
        )
        np.testing.assert_allclose(
            crystal_profile[:-1, 1:],
            np.tile(film_profile[0, 1:], (62, 1)),
        )
        np.testing.assert_array_equal(crystal_profile[-1, 1:], [0.0, 0.0])

    def test_combine_profiles_sums_coincident_positions(self):
        combined = CTRoptics.combine_profiles(
            np.array([[1.0, 2.0, 3.0]]),
            np.array([[1.0, 4.0, 5.0], [2.0, 6.0, 7.0]]),
        )

        np.testing.assert_allclose(combined, [[1.0, 6.0, 8.0], [2.0, 6.0, 7.0]])
        self.assertTrue(combined.flags.c_contiguous)

    def test_combine_profiles_merges_nearby_layer_origins(self):
        combined = CTRoptics.combine_profiles(
            np.array([[3.2728, 2.0, 3.0]]),
            np.array([[3.2731576542, 4.0, 5.0], [4.0, 6.0, 7.0]]),
        )

        np.testing.assert_allclose(combined, [[3.2728, 6.0, 8.0], [4.0, 6.0, 7.0]])

    def test_structural_layers_are_added_to_continuous_samples(self):
        combined = CTRoptics.add_structural_to_sampled_profile(
            np.array([[1.1, 2.0, 3.0], [4.0, 5.0, 6.0]]),
            np.array([[1.0, 10.0, 20.0], [1.6, 10.0, 20.0], [2.2, 10.0, 20.0]]),
        )

        np.testing.assert_allclose(
            combined,
            [[1.0, 12.0, 23.0], [1.6, 10.0, 20.0], [2.2, 10.0, 20.0], [4.0, 5.0, 6.0]],
        )

    def test_bulk_profile_repeats_towards_negative_z(self):
        bulk = unitcells.unitcell("Pt100")
        bulk.setEnergy(self.energy_eV)

        profile = bulk.optical_profile_asbulk(noUC=3)

        self.assertEqual(profile.shape, (3, 3))
        np.testing.assert_allclose(profile[:, 0], [0.0, -bulk.a[2], -2 * bulk.a[2]])
        np.testing.assert_allclose(profile[:, 1:], np.tile(profile[0, 1:], (3, 1)))

    def test_simplify_profile_conserves_finite_optical_thickness(self):
        profile = np.array(
            [
                [0.0, 2.0e-6, 2.0e-8],
                [1.0, 1.0e-6, 1.0e-8],
                [3.0, 1.1e-6, 1.1e-8],
                [6.0, 1.2e-6, 1.2e-8],
                [10.0, 5.0e-6, 5.0e-8],
                [12.0, 0.0, 0.0],
            ]
        )

        stratified = CTRoptics.stratify_profile(
            profile, delta_tolerance=0.25e-6, beta_tolerance=0.25e-8
        )
        simplified = stratified.values

        np.testing.assert_allclose(
            simplified[:, 0], [0.0, 4.25, 9.5, 12.0]
        )
        np.testing.assert_allclose(stratified.boundaries, [0.5, 8.0, 11.0])
        original_boundaries = CTRoptics.profile_boundaries(profile)
        original_integral = np.sum(
            profile[1:-1, 1:]
            * np.diff(original_boundaries)[:, None],
            axis=0,
        )
        simplified_integral = np.sum(
            simplified[1:-1, 1:]
            * np.diff(stratified.boundaries)[:, None],
            axis=0,
        )
        np.testing.assert_allclose(simplified_integral, original_integral)
        np.testing.assert_array_equal(simplified[[0, -1]], profile[[0, -1]])

    def test_simplify_profile_respects_beta_and_group_range(self):
        profile = np.array(
            [
                [0.0, 2.0e-6, 2.0e-8],
                [1.0, 1.00e-6, 1.0e-8],
                [2.0, 1.09e-6, 1.0e-8],
                [3.0, 1.18e-6, 4.0e-8],
                [4.0, 0.0, 0.0],
            ]
        )

        simplified = CTRoptics.simplify_profile(
            profile, delta_tolerance=0.15e-6, beta_tolerance=0.5e-8
        )

        np.testing.assert_allclose(simplified[:, 0], [0.0, 1.5, 3.0, 4.0])

    def test_profile_boundaries_are_centered_between_samples(self):
        profile = np.array(
            [
                [-4.0, 2.0e-6, 0.0],
                [2.0, 1.0e-6, 0.0],
                [10.0, 0.0, 0.0],
            ]
        )

        stratified = CTRoptics.stratify_profile(profile)

        np.testing.assert_array_equal(stratified.values, profile)
        np.testing.assert_allclose(stratified.boundaries, [-1.0, 6.0])

    def test_water_profile_matches_bulk_water_and_stacks_on_film(self):
        bulk = unitcells.crystal("TiO2(100)").uc_bulk
        bulk.setEnergy(self.energy_eV)
        film = CTRfilm.Film(bulk)
        film.basis[0] = 2.0
        water = CTRcalc.WaterModel(bulk.a, bulk.alpha, "step")
        water.setEnergy(self.energy_eV)
        crystal = CTRcalc.SXRDCrystal(
            bulk, film, water, stacking=np.array([1, 2])
        )

        crystal.apply_stacking()
        self.assertAlmostEqual(water.pos_absolute, film.stacking_height_absolute)
        initial_water_position = water.pos_absolute
        dz = bulk.a[2] / 2.0
        profile = water.optical_profile(noUC=10, z_step=dz)
        self.assertTrue(np.allclose(np.diff(profile[:, 0]), dz))
        density = water.pw * 18.01528 * 1.66053906660
        expected = xraydb.xray_delta_beta("H2O", density, self.energy_eV)
        np.testing.assert_allclose(profile[-1, 1:], expected[:2], rtol=5e-5)

        crystal_profile = crystal.optical_profile()
        np.testing.assert_allclose(crystal_profile[-1, 1:], expected[:2], rtol=5e-5)
        water_region = crystal_profile[crystal_profile[:, 0] >= profile[0, 0]]
        np.testing.assert_allclose(np.diff(water_region[:, 0]), dz)

        film.basis[0] = 4.0
        crystal.apply_stacking()
        self.assertAlmostEqual(water.pos_absolute, film.stacking_height_absolute)
        self.assertGreater(water.pos_absolute, initial_water_position)


class TestLayeredWavefield(unittest.TestCase):
    energy_eV = 10000.0

    @staticmethod
    def _q(n, energy_eV, alpha):
        wavelength = CTRoptics.HC_KEV_ANGSTROM / (energy_eV * 1e-3)
        k0 = 2.0 * np.pi / wavelength
        if n == 1.0:
            return k0 * np.sin(np.deg2rad(alpha))
        q = k0 * np.sqrt(
            n**2 - np.cos(np.deg2rad(alpha)) ** 2
        )
        if q.imag > 0.0:
            q = -q
        return q

    def test_single_interface_matches_fresnel_and_conserves_flux(self):
        profile = np.array([[-10.0, 1.0e-5, 0.0], [0.0, 0.0, 0.0]])
        alpha = 1.0

        for polarization in ("s", "p"):
            field = CTRoptics.solve_wavefield(
                profile, self.energy_eV, alpha, polarization
            )
            n = np.array([1.0 + 0.0j, 1.0 - 1.0e-5 + 0.0j])
            q = np.array([self._q(value, self.energy_eV, alpha) for value in n])
            admittance = q if polarization == "s" else n**2 / q
            expected_r = (admittance[0] - admittance[1]) / (
                admittance[0] + admittance[1]
            )
            expected_t = 2.0 * admittance[0] / (
                admittance[0] + admittance[1]
            )

            np.testing.assert_allclose(field.r_S, expected_r, rtol=1e-13)
            np.testing.assert_allclose(field.t_S, expected_t, rtol=1e-13)
            flux = abs(field.r_S) ** 2 + (
                admittance[1].real / admittance[0].real
            ) * abs(field.t_S) ** 2
            np.testing.assert_allclose(flux, 1.0, rtol=1e-12)
            np.testing.assert_allclose(field.A_minus[-1], 0.0, atol=1e-14)

    def test_non_vacuum_incident_medium_uses_its_tangential_wavevector(self):
        incident_delta = 4.6e-7
        substrate_delta = 2.05e-6
        profile = np.array(
            [
                [-10.0, substrate_delta, 0.0],
                [0.0, incident_delta, 0.0],
            ]
        )
        alpha = 0.02

        field = CTRoptics.solve_wavefield(
            profile, self.energy_eV, alpha, "s"
        )
        wavelength = CTRoptics.HC_KEV_ANGSTROM / (self.energy_eV * 1e-3)
        k0 = 2.0 * np.pi / wavelength
        n_incident = 1.0 - incident_delta
        expected_q_incident = (
            k0 * n_incident * np.sin(np.deg2rad(alpha))
        )

        np.testing.assert_allclose(-field.kz[0], expected_q_incident, rtol=1e-13)
        np.testing.assert_allclose(abs(field.r_S) ** 2, 1.0, rtol=1e-12)

    def test_amplitudes_satisfy_interface_continuity(self):
        profile = np.array(
            [
                [-12.0, 1.2e-5, 0.0],
                [0.0, 6.0e-6, 0.0],
                [20.0, 0.0, 0.0],
            ]
        )
        field = CTRoptics.solve_wavefield(profile, self.energy_eV, 1.2, "p")
        q = -field.kz
        admittance = field.n**2 / q

        for interface, z_value in enumerate(field.z_interfaces):
            depth = field.z_reference[interface] - z_value
            down = field.A_plus[interface] * np.exp(
                -1j * q[interface] * depth
            )
            up = field.A_minus[interface] * np.exp(
                1j * q[interface] * depth
            )
            lower_down = field.A_plus[interface + 1]
            lower_up = field.A_minus[interface + 1]
            np.testing.assert_allclose(
                down + up, lower_down + lower_up, rtol=1e-12, atol=1e-12
            )
            np.testing.assert_allclose(
                admittance[interface] * (down - up),
                admittance[interface + 1] * (lower_down - lower_up),
                rtol=1e-12,
                atol=1e-12,
            )

    def test_absorbing_substrate_field_decays_with_depth(self):
        profile = np.array([[-40.0, 1.0e-5, 2.0e-6], [0.0, 0.0, 0.0]])
        field = CTRoptics.solve_wavefield(profile, self.energy_eV, 1.0, "s")

        self.assertLess(abs(field.psi[0]), abs(field.t_S))
        q_substrate = -field.kz[-1]
        depth = field.z_interfaces[-1] - field.z[0]
        expected = abs(field.t_S) * np.exp(q_substrate.imag * depth)
        np.testing.assert_allclose(abs(field.psi[0]), expected, rtol=1e-12)

    def test_single_interface_critical_angle(self):
        delta = 1.0e-5
        profile = np.array([[-10.0, delta, 0.0], [0.0, 0.0, 0.0]])
        alpha_c = np.rad2deg(np.sqrt(2.0 * delta))

        below = CTRoptics.solve_wavefield(
            profile, self.energy_eV, 0.5 * alpha_c, "s"
        )
        above = CTRoptics.solve_wavefield(
            profile, self.energy_eV, 2.0 * alpha_c, "s"
        )

        np.testing.assert_allclose(abs(below.r_S) ** 2, 1.0, rtol=1e-12)
        self.assertLess(abs(above.r_S) ** 2, 1.0)
        self.assertLess((-below.kz[-1]).imag, 0.0)

    def test_one_film_matches_analytic_slab_expression(self):
        thickness = 35.0
        profile = np.array(
            [
                [-0.5 * thickness, 1.4e-5, 0.0],
                [0.5 * thickness, 5.0e-6, 0.0],
                [1.5 * thickness, 0.0, 0.0],
            ]
        )
        alpha = 1.5
        field = CTRoptics.solve_wavefield(profile, self.energy_eV, alpha, "s")
        n = np.array([1.0, 1.0 - 5.0e-6, 1.0 - 1.4e-5])
        q = np.array([self._q(value, self.energy_eV, alpha) for value in n])
        r01 = (q[0] - q[1]) / (q[0] + q[1])
        r12 = (q[1] - q[2]) / (q[1] + q[2])
        phase = np.exp(-1j * q[1] * thickness)
        expected = (r01 + r12 * phase**2) / (
            1.0 + r01 * r12 * phase**2
        )

        np.testing.assert_allclose(field.r_S, expected, rtol=1e-12)

    def test_crystal_reflectivity_polarizations(self):
        bulk = unitcells.unitcell("Pt100")
        bulk.setEnergy(self.energy_eV)
        crystal = CTRcalc.SXRDCrystal(bulk)
        angles = np.array([0.5, 1.0, 2.0])

        s = crystal.specular_reflectivity(angles, "s")
        p = crystal.specular_reflectivity(angles, "p")
        unpolarized = crystal.specular_reflectivity(angles, "unpolarized")

        np.testing.assert_allclose(unpolarized, 0.5 * (s + p))
        self.assertIsInstance(crystal.specular_reflectivity(1.0), float)
        self.assertEqual(crystal.wavefield(1.0).polarization, "s")

    def test_crystal_simplification_removes_redundant_bulk_layers(self):
        bulk = unitcells.unitcell("Pt100")
        bulk.setEnergy(self.energy_eV)
        crystal = CTRcalc.SXRDCrystal(bulk)
        full = crystal.optical_profile()
        simplified = crystal.simplified_optical_profile()
        angles = np.array([0.5, 1.0, 2.0])

        self.assertEqual(len(full), 31)
        self.assertEqual(len(simplified), 3)
        exact = crystal.specular_reflectivity(angles, delta_tolerance=None)
        reduced = crystal.specular_reflectivity(
            angles, delta_tolerance=1e-9
        )
        np.testing.assert_allclose(reduced, exact, rtol=1e-13, atol=1e-15)

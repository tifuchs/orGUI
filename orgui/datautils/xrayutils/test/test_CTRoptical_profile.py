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

    def test_pt100_uses_waasmaier_forward_factor_and_dispersion(self):
        cell = unitcells.unitcell("Pt100")
        cell.setEnergy(self.energy_eV)
        profile = cell.optical_profile()
        np.testing.assert_array_equal(profile, CTRoptics.optical_profile(cell))
        wavelength = CTRoptics.HC_KEV_ANGSTROM / (self.energy_eV * 1e-3)
        scale = 2.8179403262e-5 * wavelength**2 / (2.0 * np.pi)
        forward = np.sum(cell.f[:, :5], axis=1) + cell.f[:, 10]
        expected = scale * np.sum(
            cell.basis[:, 6]
            * (forward + cell.f[:, 11] + 1j * cell.f[:, 12])
        ) / cell.volume

        self.assertEqual(profile.shape, (1, 3))
        self.assertEqual(profile.dtype, np.float64)
        self.assertTrue(profile.flags.c_contiguous)
        self.assertEqual(profile[0, 0], 0.0)
        np.testing.assert_allclose(
            profile[0, 1:], [expected.real, expected.imag], rtol=1e-14
        )
        homogeneous = CTRoptics.homogeneous_bulk_profile(cell)
        np.testing.assert_allclose(
            homogeneous.values[0, 1:],
            [expected.real, expected.imag],
            rtol=1e-14,
        )

    def test_ionic_profile_uses_ionic_forward_scattering_factor(self):
        cell = CTRcalc.UnitCell([10.0, 10.0, 10.0], [90.0, 90.0, 90.0])
        cell.addAtom("O2-", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
        cell.setEnergy(self.energy_eV)

        profile = cell.optical_profile()

        wavelength = CTRoptics.HC_KEV_ANGSTROM / (self.energy_eV * 1e-3)
        scale = 2.8179403262e-5 * wavelength**2 / (2.0 * np.pi)
        ionic_f0 = np.sum(cell.f[0, :5]) + cell.f[0, 10]
        expected = scale * (
            ionic_f0 + cell.f[0, 11] + 1j * cell.f[0, 12]
        ) / cell.volume
        np.testing.assert_allclose(
            profile[0, 1:],
            [expected.real, expected.imag],
            rtol=1e-14,
        )
        self.assertGreater(ionic_f0, xraydb.atomic_number("O"))

    def test_tio2_100_has_two_homogeneous_layers(self):
        cell = unitcells.crystal("TiO2(100)").uc_bulk
        cell.setEnergy(self.energy_eV)
        profile = cell.optical_profile()
        formula_mass = xraydb.atomic_mass("Ti") + 2.0 * xraydb.atomic_mass("O")
        density = density_from_cell(formula_mass, 4, cell.volume)
        delta, beta, _ = xraydb.xray_delta_beta("TiO2", density, self.energy_eV)

        self.assertEqual(profile.shape, (2, 3))
        np.testing.assert_allclose(profile[:, 0], [-6.5807, -3.29035])
        # xraydb uses the exact neutral-atom electron count at Q=0, whereas
        # the CTR/DWBA path intentionally uses the Waasmaier--Kirfel fit at
        # Q=0 so that optical and atomic amplitudes have identical limits.
        np.testing.assert_allclose(profile[:, 1], delta, rtol=1e-3)
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

    def test_domain_normal_strain_preserves_areal_optical_content(self):
        cell = unitcells.unitcell("Pt100")
        cell.setEnergy(self.energy_eV)
        unstrained = cell.optical_profile()
        stretched = np.eye(3, 4)
        stretched[2, 2] = 2.0
        cell.coherentDomainMatrix = [stretched]
        cell.coherentDomainOccupancy = [1.0]

        profile = cell.optical_profile()

        np.testing.assert_allclose(profile[:, 1:], unstrained[:, 1:] / 2.0)

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

    def test_combine_profiles_bounds_full_group_span(self):
        combined = CTRoptics.combine_profiles(
            np.array(
                [
                    [0.0, 1.0, 0.0],
                    [0.5, 2.0, 0.0],
                    [1.0, 4.0, 0.0],
                ]
            ),
            z_tolerance=0.6,
        )

        np.testing.assert_allclose(
            combined,
            [[0.0, 3.0, 0.0], [1.0, 4.0, 0.0]],
        )

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
    def _kz(n, energy_eV, alpha):
        wavelength = CTRoptics.HC_KEV_ANGSTROM / (energy_eV * 1e-3)
        k0 = 2.0 * np.pi / wavelength
        if n == 1.0:
            return -k0 * np.sin(alpha)
        kz = -k0 * np.sqrt(n**2 - np.cos(alpha) ** 2)
        if kz.imag < 0.0 or (kz.imag == 0.0 and kz.real > 0.0):
            kz = -kz
        return kz

    def test_single_interface_matches_fresnel_and_conserves_flux(self):
        profile = np.array([[-10.0, 1.0e-5, 0.0], [0.0, 0.0, 0.0]])
        alpha = np.deg2rad(1.0)

        for polarization in ("s", "p"):
            field = CTRoptics.solve_electric_field(
                profile, self.energy_eV, alpha, polarization
            )
            n = np.array([1.0 + 0.0j, 1.0 - 1.0e-5 + 0.0j])
            kz = np.array(
                [self._kz(value, self.energy_eV, alpha) for value in n]
            )
            admittance = -kz if polarization == "s" else -(n**2) / kz
            tangential_r = (admittance[0] - admittance[1]) / (
                admittance[0] + admittance[1]
            )
            expected_r = -tangential_r if polarization == "p" else tangential_r
            expected_t = 2.0 * admittance[0] / (
                admittance[0] + admittance[1]
            )

            np.testing.assert_allclose(field.r_S, expected_r, rtol=1e-13)
            np.testing.assert_allclose(field.t_S, expected_t, rtol=1e-13)
            np.testing.assert_allclose(field.A_minus[0], field.r_S, rtol=0.0)
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
        alpha = np.deg2rad(0.02)

        field = CTRoptics.solve_electric_field(
            profile, self.energy_eV, alpha, "s"
        )
        wavelength = CTRoptics.HC_KEV_ANGSTROM / (self.energy_eV * 1e-3)
        k0 = 2.0 * np.pi / wavelength
        n_incident = 1.0 - incident_delta
        expected_kz_incident = -k0 * n_incident * np.sin(alpha)

        np.testing.assert_allclose(field.kz[0], expected_kz_incident, rtol=1e-13)
        np.testing.assert_allclose(abs(field.r_S) ** 2, 1.0, rtol=1e-12)

    def test_amplitudes_satisfy_interface_continuity(self):
        profile = np.array(
            [
                [-12.0, 1.2e-5, 0.0],
                [0.0, 6.0e-6, 0.0],
                [20.0, 0.0, 0.0],
            ]
        )
        field = CTRoptics.solve_electric_field(
            profile, self.energy_eV, np.deg2rad(1.2), "p"
        )
        kz = field.kz
        admittance = -(field.n**2) / kz

        for interface, z_value in enumerate(field.z_interfaces):
            depth = field.z_reference[interface] - z_value
            down = field.A_plus[interface] * np.exp(
                1j * kz[interface] * depth
            )
            up = field.A_minus[interface] * np.exp(
                -1j * kz[interface] * depth
            )
            lower_down = field.A_plus[interface + 1]
            lower_up = field.A_minus[interface + 1]
            if field.polarization == "p":
                tangential = down - up
                lower_tangential = lower_down - lower_up
                flux_field = down + up
                lower_flux_field = lower_down + lower_up
            else:
                tangential = down + up
                lower_tangential = lower_down + lower_up
                flux_field = down - up
                lower_flux_field = lower_down - lower_up
            np.testing.assert_allclose(
                tangential, lower_tangential, rtol=1e-12, atol=1e-12
            )
            np.testing.assert_allclose(
                admittance[interface] * flux_field,
                admittance[interface + 1] * lower_flux_field,
                rtol=1e-12,
                atol=1e-12,
            )

    def test_absorbing_substrate_field_decays_with_depth(self):
        profile = np.array([[-40.0, 1.0e-5, 2.0e-6], [0.0, 0.0, 0.0]])
        field = CTRoptics.solve_electric_field(
            profile, self.energy_eV, np.deg2rad(1.0), "s"
        )
        sampled = CTRoptics.sample_electric_field(field, profile[0, 0])

        self.assertLess(abs(sampled), abs(field.t_S))
        kz_substrate = field.kz[-1]
        depth = field.z_interfaces[-1] - profile[0, 0]
        expected = abs(field.t_S) * np.exp(-kz_substrate.imag * depth)
        np.testing.assert_allclose(abs(sampled), expected, rtol=1e-12)

    def test_single_interface_critical_angle(self):
        delta = 1.0e-5
        profile = np.array([[-10.0, delta, 0.0], [0.0, 0.0, 0.0]])
        alpha_c = np.sqrt(2.0 * delta)

        below = CTRoptics.solve_electric_field(
            profile, self.energy_eV, 0.5 * alpha_c, "s"
        )
        above = CTRoptics.solve_electric_field(
            profile, self.energy_eV, 2.0 * alpha_c, "s"
        )

        np.testing.assert_allclose(abs(below.r_S) ** 2, 1.0, rtol=1e-12)
        self.assertLess(abs(above.r_S) ** 2, 1.0)
        self.assertGreater(below.kz[-1].imag, 0.0)

    def test_p_polarization_is_finite_at_exact_substrate_critical_angle(self):
        delta = 1.0e-5
        profile = np.array([[-10.0, delta, 0.0], [0.0, 0.0, 0.0]])
        alpha_c = np.arccos(1.0 - delta)

        field = CTRoptics.solve_electric_field(
            profile, self.energy_eV, alpha_c, "p"
        )

        np.testing.assert_allclose(-field.kz[-1], 0.0, atol=0.0)
        np.testing.assert_allclose(field.r_S, 1.0, rtol=1e-14)
        np.testing.assert_allclose(field.t_S, 0.0, atol=0.0)
        sampled = CTRoptics.sample_electric_field(field, profile[:, 0])
        self.assertTrue(np.all(np.isfinite(sampled)))
        self.assertTrue(np.all(np.isfinite(field.A_plus)))
        self.assertTrue(np.all(np.isfinite(field.A_minus)))

    def test_p_polarization_is_continuous_at_internal_critical_angle(self):
        delta = 1.0e-5
        profile = np.array(
            [
                [-40.0, 2.0 * delta, 0.0],
                [-10.0, delta, 0.0],
                [20.0, 0.0, 0.0],
            ]
        )
        alpha_c = np.arccos(1.0 - delta)

        exact = CTRoptics.solve_electric_field(
            profile, self.energy_eV, alpha_c, "p"
        )
        below = CTRoptics.solve_electric_field(
            profile, self.energy_eV, alpha_c * (1.0 - 1.0e-8), "p"
        )
        above = CTRoptics.solve_electric_field(
            profile, self.energy_eV, alpha_c * (1.0 + 1.0e-8), "p"
        )

        np.testing.assert_allclose(-exact.kz[1], 0.0, atol=0.0)
        exact_E = CTRoptics.sample_electric_field(exact, profile[:, 0])
        below_E = CTRoptics.sample_electric_field(below, profile[:, 0])
        above_E = CTRoptics.sample_electric_field(above, profile[:, 0])
        self.assertTrue(np.all(np.isfinite(exact_E)))
        self.assertTrue(np.isfinite(exact.r_S))
        np.testing.assert_allclose(below.r_S, exact.r_S, rtol=1e-7)
        np.testing.assert_allclose(above.r_S, exact.r_S, rtol=1e-7)
        np.testing.assert_allclose(below_E, exact_E, rtol=1e-7)
        np.testing.assert_allclose(above_E, exact_E, rtol=1e-7)

    def test_vectorized_reflection_matches_scalar_wavefields(self):
        delta = 1.0e-5
        profile = np.array(
            [
                [-40.0, 2.0 * delta, 0.0],
                [-10.0, delta, 0.0],
                [20.0, 0.0, 0.0],
            ]
        )
        alpha_c = np.arccos(1.0 - delta)
        angles = np.array([np.deg2rad(0.1), alpha_c, np.deg2rad(0.5), np.deg2rad(1.0)])

        for polarization in ("s", "p"):
            actual = CTRoptics._specular_reflection(
                profile, self.energy_eV, angles, polarization
            )
            expected = np.array(
                [
                    CTRoptics.solve_electric_field(
                        profile, self.energy_eV, angle, polarization
                    ).r_S
                    for angle in angles
                ]
            )

            np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    def test_vectorized_wavefield_matches_scalar_wavefields(self):
        delta = 1.0e-5
        profile = np.array(
            [
                [-40.0, 2.0 * delta, 0.0],
                [-10.0, delta, 0.0],
                [20.0, 0.0, 0.0],
            ]
        )
        alpha_c = np.arccos(1.0 - delta)
        angles = np.array(
            [
                [np.deg2rad(0.1), alpha_c],
                [np.deg2rad(0.5), np.deg2rad(1.0)],
            ]
        )

        vector = CTRoptics.solve_electric_field(
            profile, self.energy_eV, angles, "p"
        )
        vector_E = CTRoptics.sample_electric_field(vector, profile[:, 0])

        self.assertEqual(vector_E.shape, (len(profile), *angles.shape))
        self.assertEqual(vector.kz.shape, (len(profile), *angles.shape))
        self.assertEqual(vector.A_plus.shape, (len(profile), *angles.shape))
        self.assertEqual(vector.A_minus.shape, (len(profile), *angles.shape))
        self.assertEqual(vector.r_S.shape, angles.shape)
        self.assertEqual(vector.t_S.shape, angles.shape)
        for index in np.ndindex(angles.shape):
            scalar = CTRoptics.solve_electric_field(
                profile, self.energy_eV, angles[index], "p"
            )
            scalar_E = CTRoptics.sample_electric_field(
                scalar, profile[:, 0]
            )
            np.testing.assert_allclose(
                vector_E[(slice(None), *index)],
                scalar_E,
                rtol=1e-12,
                atol=1e-12,
            )
            np.testing.assert_allclose(
                vector.A_plus[(slice(None), *index)],
                scalar.A_plus,
                rtol=1e-12,
                atol=1e-12,
            )
            np.testing.assert_allclose(
                vector.A_minus[(slice(None), *index)],
                scalar.A_minus,
                rtol=1e-12,
                atol=1e-12,
            )
            np.testing.assert_allclose(vector.r_S[index], scalar.r_S)
            np.testing.assert_allclose(vector.t_S[index], scalar.t_S)

    def test_one_film_matches_analytic_slab_expression(self):
        thickness = 35.0
        profile = np.array(
            [
                [-0.5 * thickness, 1.4e-5, 0.0],
                [0.5 * thickness, 5.0e-6, 0.0],
                [1.5 * thickness, 0.0, 0.0],
            ]
        )
        alpha = np.deg2rad(1.5)
        field = CTRoptics.solve_electric_field(
            profile, self.energy_eV, alpha, "s"
        )
        n = np.array([1.0, 1.0 - 5.0e-6, 1.0 - 1.4e-5])
        kz = np.array(
            [self._kz(value, self.energy_eV, alpha) for value in n]
        )
        r01 = (kz[0] - kz[1]) / (kz[0] + kz[1])
        r12 = (kz[1] - kz[2]) / (kz[1] + kz[2])
        phase = np.exp(1j * kz[1] * thickness)
        expected = (r01 + r12 * phase**2) / (
            1.0 + r01 * r12 * phase**2
        )

        np.testing.assert_allclose(field.r_S, expected, rtol=1e-12)

    def test_crystal_reflectivity_polarizations(self):
        bulk = unitcells.unitcell("Pt100")
        bulk.setEnergy(self.energy_eV)
        crystal = CTRcalc.SXRDCrystal(bulk)
        angles = np.deg2rad([0.5, 1.0, 2.0])

        s = crystal.specular_reflectivity(angles, "s")
        p = crystal.specular_reflectivity(angles, "p")
        unpolarized = crystal.specular_reflectivity(angles, "unpolarized")

        np.testing.assert_allclose(unpolarized, 0.5 * (s + p))
        self.assertIsInstance(
            crystal.specular_reflectivity(np.deg2rad(1.0)), float
        )
        self.assertEqual(
            crystal.wavefield(np.deg2rad(1.0)).polarization, "s"
        )
        grid = angles.reshape(1, -1)
        np.testing.assert_allclose(
            crystal.specular_reflectivity(grid, "p"),
            p.reshape(1, -1),
        )

    def test_crystal_simplification_removes_redundant_bulk_layers(self):
        bulk = unitcells.unitcell("Pt100")
        bulk.setEnergy(self.energy_eV)
        crystal = CTRcalc.SXRDCrystal(bulk)
        full = crystal.optical_profile()
        simplified = crystal.simplified_optical_profile()
        angles = np.deg2rad([0.5, 1.0, 2.0])

        self.assertEqual(len(full), 31)
        self.assertEqual(len(simplified), 3)
        exact = crystal.specular_reflectivity(angles, delta_tolerance=None)
        reduced = crystal.specular_reflectivity(
            angles, delta_tolerance=1e-9
        )
        np.testing.assert_allclose(reduced, exact, rtol=1e-13, atol=1e-15)

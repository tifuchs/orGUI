# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
import copy
import importlib.util
import subprocess
import sys
import unittest
from unittest import mock

import numpy as np

from .. import CTRsymmetry
from ..CTRuc import UnitCell


HAS_PYXTAL = importlib.util.find_spec("pyxtal") is not None
HAS_PYMATGEN = importlib.util.find_spec("pymatgen") is not None


class TestSymmetryUtilities(unittest.TestCase):
    @staticmethod
    def make_metadata_unitcell():
        site = CTRsymmetry.WyckoffSiteSpec(
            site_id="O_2a",
            element="O",
            wyckoff_label="2a",
            coordinates=(),
            representative_parent_fractional=(0.2, 0.3, 0.4),
            occ=0.9,
            iDW=0.5,
            oDW=0.6,
        )
        atoms = []
        for index, (position, x_factor) in enumerate(
            (((0.2, 0.3, 0.4), 1.0), ((0.8, 0.7, 0.6), -1.0))
        ):
            atoms.append(
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=index,
                    element="O",
                    site_id="O_2a",
                    wyckoff_label="2a",
                    parent_fractional=np.asarray(position),
                    surface_fractional=np.asarray(position),
                    layer=0,
                    site_couplings=(
                        CTRsymmetry.WyckoffSiteCoupling(
                            atom_index=index,
                            coordinate="x",
                            axis="x",
                            factor=x_factor,
                            site_id="O_2a",
                        ),
                    ),
                )
            )
        model = CTRsymmetry.SurfaceSymmetryModel(
            CTRsymmetry.SurfaceCellSpec(
                (4.0, 4.0, 4.0),
                (90.0, 90.0, 90.0),
                np.identity(3),
                translation_range=0,
            ),
            (site,),
            atoms,
        )
        return model.build_unitcell("test")

    def test_wyckoff_returns_one_site_and_reports_site_parameters(self):
        unitcell = self.make_metadata_unitcell()

        self.assertEqual(unitcell.wyckoff("O_2a"), unitcell.wyckoff_sites()[0])
        self.assertEqual(unitcell.wyckoff("O_2a")["occ"], 0.9)
        self.assertEqual(unitcell.wyckoff("O_2a")["iDW"], 0.5)
        self.assertEqual(unitcell.wyckoff("O_2a")["oDW"], 0.6)
        with self.assertRaisesRegex(ValueError, "Unknown Wyckoff site"):
            unitcell.wyckoff("missing")

    def test_set_wyckoff_atom_parameter_updates_all_site_atoms(self):
        unitcell = self.make_metadata_unitcell()

        unitcell.set_wyckoff_atom_parameter("O_2a", "occ", 0.75)
        unitcell.set_wyckoff_atom_parameter("O_2a", "z", [0.1, 0.9])

        np.testing.assert_allclose(unitcell.basis[:, 6], 0.75)
        np.testing.assert_allclose(unitcell.basis_0[:, 6], 0.75)
        np.testing.assert_allclose(unitcell.basis[:, 3], [0.1, 0.9])
        np.testing.assert_allclose(unitcell.basis_0[:, 3], [0.1, 0.9])
        self.assertEqual(unitcell.wyckoff("O_2a")["occ"], 0.75)
        np.testing.assert_allclose(
            [atom.surface_fractional[2] for atom in unitcell.symmetry_metadata.atoms],
            [0.1, 0.9],
        )

    def test_set_wyckoff_atom_parameter_allows_uninitialized_fit_parameters(self):
        unitcell = self.make_metadata_unitcell()
        unitcell.addFitParameter(
            ([0], ["z"]),
            name="manual_top_oxygen_z",
        )
        original_z = unitcell.basis[:, 3].copy()

        unitcell.set_wyckoff_atom_parameter("O_2a", "iDW", 0.35)

        np.testing.assert_allclose(unitcell.basis[:, 4], 0.35)
        np.testing.assert_allclose(unitcell.basis_0[:, 4], 0.35)
        np.testing.assert_allclose(unitcell.basis[:, 3], original_z)
        self.assertIsNone(unitcell.parameters["absolute"][0].value)

    def test_set_wyckoff_site_parameter_propagates_parent_coordinate(self):
        unitcell = self.make_metadata_unitcell()

        unitcell.set_wyckoff_site_parameter("O_2a", "x", 0.25)

        np.testing.assert_allclose(unitcell.basis[:, 1], [0.25, 0.75])
        np.testing.assert_allclose(unitcell.basis_0[:, 1], [0.25, 0.75])
        self.assertEqual(
            unitcell.wyckoff("O_2a")["representative_parent_fractional"],
            (0.25, 0.3, 0.4),
        )
        np.testing.assert_allclose(
            [atom.parent_fractional[0] for atom in unitcell.symmetry_metadata.atoms],
            [0.25, 0.75],
        )

    def test_set_wyckoff_site_parameter_sets_site_wide_physical_values(self):
        unitcell = self.make_metadata_unitcell()

        unitcell.set_wyckoff_site_parameter("O_2a", "iDW", 0.7)
        unitcell.set_wyckoff_site_parameter("O_2a", "oDW", 0.8)
        unitcell.set_wyckoff_site_parameter("O_2a", "occ", 0.65)

        np.testing.assert_allclose(unitcell.basis[:, 4], 0.7)
        np.testing.assert_allclose(unitcell.basis[:, 5], 0.8)
        np.testing.assert_allclose(unitcell.basis[:, 6], 0.65)
        self.assertEqual(unitcell.wyckoff("O_2a")["iDW"], 0.7)
        self.assertEqual(unitcell.wyckoff("O_2a")["oDW"], 0.8)
        self.assertEqual(unitcell.wyckoff("O_2a")["occ"], 0.65)

    def test_legacy_expanded_couplings_remain_readable(self):
        lines = [
            "parent_a: 4 4 4",
            "parent_alpha: 90 90 90",
            "surface_transform:",
            "1 0 0",
            "0 1 0",
            "0 0 1",
            "wyckoff_sites:",
            "site_id element wyckoff_label variables representative_x "
            "representative_y representative_z occ iDW oDW",
            "O_1a O 1a u=0.2 0.2 0 0 1 0.5 0.5",
            "wyckoff_atoms:",
            "atom_index element site_id wyckoff_label parent_x parent_y "
            "parent_z surface_x surface_y surface_z layer",
            "0 O O_1a 1a 0.2 0 0 0.2 0 0 0",
            "wyckoff_couplings:",
            "atom_index coordinate site_id variable constant factor",
            "0 x O_1a u 0 1",
            "wyckoff_site_couplings:",
            "atom_index coordinate site_id axis factor",
            "0 x O_1a x 1",
        ]

        model = CTRsymmetry.symmetry_metadata_from_lines(lines)

        self.assertEqual(len(model.wyckoff_couplings("O_1a")), 1)
        self.assertEqual(len(model.wyckoff_site_couplings("O_1a")), 1)

    def test_duplicate_wyckoff_site_ids_are_disambiguated(self):
        sites = CTRsymmetry._with_unique_site_ids(
            (
                CTRsymmetry.WyckoffSiteSpec(
                    site_id="O_4f",
                    element="O",
                    wyckoff_label="4f",
                    coordinates=(),
                ),
                CTRsymmetry.WyckoffSiteSpec(
                    site_id="O_4f",
                    element="O",
                    wyckoff_label="4f",
                    coordinates=(),
                ),
                CTRsymmetry.WyckoffSiteSpec(
                    site_id="Ru_2a",
                    element="Ru",
                    wyckoff_label="2a",
                    coordinates=(),
                ),
            )
        )

        self.assertEqual(
            [site.site_id for site in sites],
            ["O_4f_1", "O_4f_2", "Ru_2a"],
        )


@unittest.skipUnless(
    HAS_PYXTAL and HAS_PYMATGEN,
    "PyXtal symmetry tests require PyXtal",
)
class TestPyxtalRutileSurfaceSymmetry(unittest.TestCase):
    def make_rutile_seed(self, oxygen_u=0.30569):
        """Return a conventional rutile parent structure."""
        from pymatgen.core import Lattice, Structure

        u = oxygen_u
        return Structure(
            Lattice.tetragonal(4.653255, 2.9692),
            ["Ru", "Ru", "O", "O", "O", "O"],
            [
                [0.0, 0.0, 0.0],
                [0.5, 0.5, 0.5],
                [u, u, 0.0],
                [1.0 - u, 1.0 - u, 0.0],
                [0.5 - u, 0.5 + u, 0.5],
                [0.5 + u, 0.5 - u, 0.5],
            ],
        )

    def make_rutile_110_unitcell(self):
        surface_spec = CTRsymmetry.rutile_110_surface_spec(
            (4.653255, 4.653255, 2.9692),
        )
        model = CTRsymmetry.model_from_seed(
            self.make_rutile_seed(),
            surface_spec,
            tol=1e-3,
            iDW=0.5033,
            oDW=0.5033,
        )
        return model.build_unitcell("RuO2")

    def make_pyxtal_generated_unitcell(
        self,
        spacegroup,
        wyckoff_label,
        variable_values,
        lattice,
        transform=None,
        element="Si",
    ):
        """Generate a UnitCell from one PyXtal Wyckoff position."""
        from pymatgen.core import Structure
        from pyxtal.symmetry import Group

        wyckoff_position = Group(spacegroup).get_wyckoff_position(wyckoff_label)
        coordinates = wyckoff_position.get_all_positions(variable_values)
        seed = Structure(lattice, [element] * len(coordinates), coordinates)
        surface_spec = CTRsymmetry.SurfaceCellSpec(
            lattice.abc,
            lattice.angles,
            np.identity(3) if transform is None else transform,
            layer_origins=(0.0,),
            translation_range=0,
        )
        model = CTRsymmetry.model_from_seed(seed, surface_spec, tol=1e-3)
        return model.build_unitcell(element)

    def assert_variables_allclose(self, actual, expected):
        """Compare Wyckoff ``variables`` dicts, tolerating float rounding.

        pyxtal/spglib can return values that are a few ULP off an exact
        input (e.g. 0.11999999999999997 instead of 0.12), so this avoids
        exact dict equality on floats.
        """
        self.assertEqual(set(actual), set(expected))
        np.testing.assert_allclose(
            [actual[key] for key in expected], [expected[key] for key in expected]
        )

    @staticmethod
    def coupling_factor_vectors(unitcell, site_id):
        grouped = {}
        for coupling in unitcell.wyckoff_couplings(site_id):
            key = (coupling.atom_index, coupling.coordinate)
            grouped.setdefault(key, {})[coupling.variable] = coupling.factor
        return {
            tuple(round(factors.get(variable, 0.0), 8) for variable in ("u", "v", "w"))
            for factors in grouped.values()
        }

    @staticmethod
    def site_coupling_factor_vectors(unitcell, site_id):
        grouped = {}
        for coupling in unitcell.wyckoff_site_couplings(site_id):
            key = (coupling.atom_index, coupling.coordinate)
            grouped.setdefault(key, {})[coupling.axis] = coupling.factor
        return {
            tuple(round(factors.get(axis, 0.0), 8) for axis in ("x", "y", "z"))
            for factors in grouped.values()
        }

    def test_pyxtal_assigns_rutile_wyckoff_sites(self):
        sites, number, symbol = CTRsymmetry.sites_from_seed(
            self.make_rutile_seed(),
            tol=1e-3,
        )

        self.assertEqual(number, 136)
        self.assertEqual(symbol, "P42/mnm")
        self.assertEqual([site.site_id for site in sites], ["Ru_2a", "O_4f"])
        self.assertEqual(sites[1].variables, {"u": 0.30569})

    def test_possible_wyckoff_positions_exposes_group_table(self):
        positions = CTRsymmetry.possible_wyckoff_positions(136)

        labels = [position["label"] for position in positions]
        self.assertIn("2a", labels)
        self.assertIn("4f", labels)
        oxygen_position = next(
            position for position in positions if position["label"] == "4f"
        )
        self.assertEqual(oxygen_position["dof"], 1)

    def test_rutile_110_generation_reproduces_surface_basis(self):
        unitcell = self.make_rutile_110_unitcell()

        self.assertEqual(unitcell.basis.shape, (12, 8))
        self.assertEqual(unitcell.parameters["absolute"], [])
        self.assertEqual(unitcell.parameters["relative"], [])
        np.testing.assert_allclose(unitcell.a, [6.5807, 2.9692, 6.5807], atol=5e-5)
        self.assertEqual(unitcell.layerpos, {1.0: 0.0, 2.0: 0.5})

        expected = {
            ("O", 0.5, 0.0, 0.80569, 2.0),
            ("O", 0.0, 0.0, 0.69431, 2.0),
            ("Ru", 0.5, 0.0, 0.5, 2.0),
            ("Ru", 0.0, 0.5, 0.5, 2.0),
            ("O", 0.30569, 0.5, 0.5, 2.0),
            ("O", 0.69431, 0.5, 0.5, 2.0),
            ("O", 0.0, 0.0, 0.30569, 1.0),
            ("O", 0.5, 0.0, 0.19431, 1.0),
            ("Ru", 0.0, 0.0, 0.0, 1.0),
            ("Ru", 0.5, 0.5, 0.0, 1.0),
            ("O", 0.19431, 0.5, 0.0, 1.0),
            ("O", 0.80569, 0.5, 0.0, 1.0),
        }
        actual = {
            (
                name,
                round(row[1], 5),
                round(row[2], 5),
                round(row[3], 5),
                row[7],
            )
            for name, row in zip(unitcell.names, unitcell.basis)
        }
        self.assertEqual(actual, expected)

    def test_wyckoff_query_exposes_oxygen_u_couplings(self):
        unitcell = self.make_rutile_110_unitcell()

        sites = unitcell.wyckoff_sites()
        self.assertEqual([site["site_id"] for site in sites], ["Ru_2a", "O_4f"])
        self.assertEqual(sites[0]["spacegroup_number"], 136)
        self.assert_variables_allclose(sites[1]["variables"], {"u": 0.30569})
        self.assertEqual(sites[1]["status"], "metadata_only")

        couplings = unitcell.wyckoff_couplings("O_4f")
        self.assertEqual(len(couplings), 8)
        self.assertTrue(all(coupling.variable == "u" for coupling in couplings))
        expressions = {
            (round(coupling.constant, 5), round(coupling.factor, 5))
            for coupling in couplings
        }
        self.assertEqual(
            expressions,
            {(0.0, 1.0), (1.0, -1.0), (0.5, 1.0), (0.5, -1.0)},
        )

    def test_wyckoff_coordinate_parameter_preserves_rutile_u_symmetry(self):
        unitcell = self.make_rutile_110_unitcell()
        couplings = unitcell.wyckoff_couplings("O_4f")
        original = unitcell.basis.copy()

        parameter = unitcell.addWyckoffParameter(
            "O_4f",
            "u",
            absolute_limits=(0.2, 0.4),
        )

        self.assertEqual(parameter.settings["wyckoff"]["kind"], "coordinate")
        self.assertEqual(parameter.settings["wyckoff"]["value_kind"], "absolute")
        np.testing.assert_allclose(unitcell.getInitialParameters(), [0.30569])
        np.testing.assert_allclose(
            unitcell.getStartParamAndLimits()[1:],
            ([0.2], [0.4]),
        )
        self.assertEqual(unitcell.wyckoff_sites()[1]["status"], "symmetry_preserving")

        unitcell.setFitParameters([0.31569])

        for coupling in couplings:
            column = unitcell.parameterLookup[coupling.coordinate]
            expected = original[coupling.atom_index, column] + 0.01 * coupling.factor
            self.assertAlmostEqual(
                unitcell.basis[coupling.atom_index, column],
                expected,
            )

    def test_wyckoff_legacy_delta_parameter_dict_matches_absolute_result(self):
        """A pre-fix, ``value_kind="delta"`` saved Wyckoff parameter dict
        must reproduce the same basis as the current absolute-value
        convention, since old saved fits stored the delta from the site's
        reference value rather than the absolute coordinate."""
        unitcell = self.make_rutile_110_unitcell()
        unitcell.addWyckoffParameter("O_4f", "u", absolute_limits=(0.2, 0.4))
        unitcell.setFitParameters([0.31569])
        expected_basis = unitcell.basis.copy()

        saved = unitcell.parametersToDict()
        (key, param_dict), = saved["relative"].items()
        wyckoff_settings = param_dict["settings"]["wyckoff"]
        self.assertEqual(wyckoff_settings["value_kind"], "absolute")
        reference_value = wyckoff_settings["reference_value"]

        # Simulate a parameter dict saved before value_kind="absolute" was
        # introduced: no reference_value, value_kind="delta", and the stored
        # value is the delta from the reference rather than the absolute
        # coordinate.
        legacy_dict = copy.deepcopy(saved)
        legacy_settings = legacy_dict["relative"][key]["settings"]["wyckoff"]
        legacy_settings["value_kind"] = "delta"
        del legacy_settings["reference_value"]
        legacy_dict["relative"][key]["value"] = param_dict["value"] - reference_value

        restored = self.make_rutile_110_unitcell()
        restored.parametersFromDict(legacy_dict)

        np.testing.assert_allclose(restored.basis, expected_basis)

    def test_wyckoff_site_parameter_displaces_fixed_rutile_site(self):
        unitcell = self.make_rutile_110_unitcell()
        couplings = [
            coupling
            for coupling in unitcell.wyckoff_site_couplings("Ru_2a")
            if coupling.axis == "x"
        ]
        original = unitcell.basis.copy()

        parameter = unitcell.addWyckoffShift(
            "Ru_2a",
            "x",
            absolute_limits=(-0.1, 0.1),
        )

        self.assertEqual(parameter.settings["wyckoff"]["kind"], "site_displacement")
        self.assertEqual(unitcell.wyckoff_sites()[0]["status"], "site_displaced")
        np.testing.assert_allclose(
            unitcell.getStartParamAndLimits()[1:],
            ([-0.1], [0.1]),
        )

        unitcell.setFitParameters([0.02])

        for coupling in couplings:
            column = unitcell.parameterLookup[coupling.coordinate]
            expected = original[coupling.atom_index, column] + 0.02 * coupling.factor
            self.assertAlmostEqual(
                unitcell.basis[coupling.atom_index, column],
                expected,
            )

    def test_wyckoff_site_parameter_lowers_oxygen_site_symmetry(self):
        unitcell = self.make_rutile_110_unitcell()

        unitcell.addWyckoffParameter("O_4f", "u")
        self.assertEqual(unitcell.wyckoff_sites()[1]["status"], "symmetry_preserving")

        unitcell = self.make_rutile_110_unitcell()
        unitcell.addWyckoffShift("O_4f", "x")

        self.assertEqual(unitcell.wyckoff_sites()[1]["status"], "site_displaced")

    def test_wyckoff_coordinate_parameters_fit_all_site_variables(self):
        from pymatgen.core import Lattice

        unitcell = self.make_pyxtal_generated_unitcell(
            62,
            "8d",
            [0.12, 0.23, 0.34],
            Lattice.orthorhombic(4.0, 5.0, 6.0),
            element="Mg",
        )
        original = unitcell.basis.copy()

        parameters = unitcell.addWyckoffParameters("Mg_8d")

        self.assertEqual(
            [par.settings["wyckoff"]["variable"] for par in parameters],
            ["u", "v", "w"],
        )
        np.testing.assert_allclose(
            unitcell.getInitialParameters(),
            [0.12, 0.23, 0.34],
        )
        self.assertEqual(unitcell.wyckoff_sites()[0]["status"], "symmetry_preserving")

        unitcell.setFitParameters([0.13, 0.21, 0.37])

        deltas = {"u": 0.01, "v": -0.02, "w": 0.03}
        for coordinate_name in ("x", "y", "z"):
            column = unitcell.parameterLookup[coordinate_name]
            expected = original[:, column].copy()
            for coupling in unitcell.wyckoff_couplings("Mg_8d"):
                if coupling.coordinate == coordinate_name:
                    expected[coupling.atom_index] += (
                        coupling.factor * deltas[coupling.variable]
                    )
            np.testing.assert_allclose(unitcell.basis[:, column], expected)

    def test_fixed_wyckoff_site_rejects_coordinate_parameter(self):
        unitcell = self.make_rutile_110_unitcell()

        with self.assertRaisesRegex(ValueError, "no positional coordinate variables"):
            unitcell.addWyckoffParameter("Ru_2a", "u")

    def test_trigonal_general_site_exposes_u_minus_v_couplings(self):
        from pymatgen.core import Lattice

        unitcell = self.make_pyxtal_generated_unitcell(
            152,
            "6c",
            [0.12, 0.27, 0.34],
            Lattice.hexagonal(4.0, 6.0),
        )

        self.assert_variables_allclose(
            unitcell.wyckoff_sites()[0]["variables"],
            {"u": 0.12, "v": 0.27, "w": 0.34},
        )
        vectors = self.coupling_factor_vectors(unitcell, "Si_6c")

        self.assertIn((1.0, -1.0, 0.0), vectors)
        self.assertTrue(
            any(abs(vector[0]) == 1.0 and abs(vector[1]) == 1.0 for vector in vectors)
        )

    def test_surface_transform_can_mix_u_v_w_into_one_coordinate(self):
        from pymatgen.core import Lattice

        parent_to_surface = np.asarray(
            [
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        transform = np.linalg.inv(parent_to_surface)
        unitcell = self.make_pyxtal_generated_unitcell(
            62,
            "8d",
            [0.12, 0.23, 0.34],
            Lattice.orthorhombic(4.0, 5.0, 6.0),
            transform=transform,
            element="Mg",
        )

        vectors = self.coupling_factor_vectors(unitcell, "Mg_8d")

        self.assertIn((1.0, 1.0, 1.0), vectors)

    def test_general_site_displacement_exposes_operation_factor_vectors(self):
        from pymatgen.core import Lattice

        unitcell = self.make_pyxtal_generated_unitcell(
            62,
            "8d",
            [0.12, 0.23, 0.34],
            Lattice.orthorhombic(4.0, 5.0, 6.0),
            element="Mg",
        )

        vectors = self.site_coupling_factor_vectors(unitcell, "Mg_8d")

        self.assertIn((1.0, 0.0, 0.0), vectors)
        self.assertIn((-1.0, 0.0, 0.0), vectors)
        self.assertIn((0.0, 0.0, 1.0), vectors)
        self.assertIn((0.0, 0.0, -1.0), vectors)

    def test_surface_transform_can_mix_site_displacement_axes(self):
        from pymatgen.core import Lattice

        parent_to_surface = np.asarray(
            [
                [1.0, 1.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        transform = np.linalg.inv(parent_to_surface)
        unitcell = self.make_pyxtal_generated_unitcell(
            62,
            "8d",
            [0.12, 0.23, 0.34],
            Lattice.orthorhombic(4.0, 5.0, 6.0),
            transform=transform,
            element="Mg",
        )

        vectors = self.site_coupling_factor_vectors(unitcell, "Mg_8d")

        self.assertIn((1.0, 1.0, 1.0), vectors)

    def test_manual_atom_parameter_marks_wyckoff_site_as_partially_overridden(self):
        unitcell = self.make_rutile_110_unitcell()
        atom_index = unitcell.wyckoff_couplings("O_4f")[0].atom_index

        unitcell.addFitParameter(
            ([atom_index], ["z"]),
            limits=(-1.0, 2.0),
            name="manual_oxygen_coordinate",
        )

        self.assertEqual(unitcell.wyckoff_sites()[1]["status"], "partially_overridden")

    def test_symmetry_metadata_round_trips_without_pyxtal_construction(self):
        unitcell = self.make_rutile_110_unitcell()
        text = unitcell.toStr()

        self.assertIn("wyckoff_coupling_matrices:", text)
        self.assertIn("wyckoff_site_coupling_matrices:", text)
        self.assertNotIn("\nwyckoff_couplings:\n", text)
        self.assertNotIn("\nwyckoff_site_couplings:\n", text)

        restored = UnitCell.fromStr(text)

        self.assert_variables_allclose(
            restored.wyckoff_sites()[1]["variables"], {"u": 0.30569}
        )
        self.assertEqual(len(restored.wyckoff_couplings("O_4f")), 8)
        self.assertGreater(len(restored.wyckoff_site_couplings("Ru_2a")), 0)
        original_couplings = sorted(
            (
                coupling.atom_index,
                coupling.coordinate,
                coupling.variable,
                round(coupling.constant, 10),
                round(coupling.factor, 10),
            )
            for coupling in unitcell.wyckoff_couplings()
        )
        restored_couplings = sorted(
            (
                coupling.atom_index,
                coupling.coordinate,
                coupling.variable,
                round(coupling.constant, 10),
                round(coupling.factor, 10),
            )
            for coupling in restored.wyckoff_couplings()
        )
        self.assertEqual(restored_couplings, original_couplings)
        restored.addWyckoffParameter("O_4f", "u", absolute_limits=(0.2, 0.4))
        self.assertEqual(restored.wyckoff_sites()[1]["status"], "symmetry_preserving")

        restored = UnitCell.fromStr(text)
        restored.addWyckoffShift("Ru_2a", "x", absolute_limits=(-0.1, 0.1))
        self.assertEqual(restored.wyckoff_sites()[0]["status"], "site_displaced")

    def test_symmetry_metadata_round_trip_preserves_parent_lattice(self):
        unitcell = self.make_rutile_110_unitcell()

        restored = UnitCell.fromStr(unitcell.toStr())
        rebuilt = restored.symmetry_metadata.build_unitcell("rebuilt")

        np.testing.assert_allclose(
            restored.symmetry_metadata.surface_spec.parent_a,
            unitcell.symmetry_metadata.surface_spec.parent_a,
        )
        np.testing.assert_allclose(
            restored.symmetry_metadata.surface_spec.parent_alpha,
            unitcell.symmetry_metadata.surface_spec.parent_alpha,
        )
        self.assertEqual(len(rebuilt.basis), len(restored.basis))
        self.assertEqual(rebuilt.names, restored.names)
        np.testing.assert_allclose(rebuilt.basis[:, 1:4], restored.basis[:, 1:4])
        np.testing.assert_allclose(rebuilt.a, unitcell.a, atol=1e-8)
        np.testing.assert_allclose(rebuilt.alpha, unitcell.alpha, atol=1e-8)


@unittest.skipUnless(
    HAS_PYXTAL and HAS_PYMATGEN,
    "PyXtal symmetry tests require PyXtal",
)
class TestPyxtalSeedSetting(unittest.TestCase):
    """Seeds whose cell differs from the spglib standardized setting."""

    # Columns: C2/c cell vectors of a tilted perovskite -> 2x2x2 pseudocubic.
    C2C_TO_SURFACE = np.asarray([[0, 1, 1], [0, -1, 1], [1, 0, 1]], float).T
    ATOL_ANGSTROM = 1e-10

    @staticmethod
    def make_c2c_seed(o_general=(0.282, 0.001, 0.283)):
        """Tilted LaNiO3 in C2/c; spglib standardizes it as (x, -y, 2x - z)."""
        from pymatgen.core import Lattice, Structure

        a_ip, c_pc = 3.905, 3.807
        b = c = np.sqrt(2.0) * a_ip
        a = np.hypot(2.0 * c_pc, c)
        beta = 180.0 - np.degrees(np.arctan2(2.0 * c_pc, c))
        return Structure.from_spacegroup(
            15,
            Lattice.monoclinic(a, b, c, beta),
            ["La", "Ni", "O", "O"],
            [
                [0.0, 0.25, 0.25],
                [0.25, 0.25, 0.0],
                [0.0, 0.810, 0.25],
                list(o_general),
            ],
        )

    def assert_positions_match(self, elements, fractional, structure):
        """Assert both atom sets coincide modulo lattice translations."""
        lattice = structure.lattice.matrix
        self.assertEqual(len(fractional), len(structure))
        for element, position in zip(elements, fractional):
            difference = structure.frac_coords - np.asarray(position)
            difference -= np.round(difference)
            distances = np.linalg.norm(difference @ lattice, axis=1)
            same_element = np.asarray(
                [site.species.elements[0].symbol == element for site in structure]
            )
            self.assertLess(distances[same_element].min(), self.ATOL_ANGSTROM)

    def assert_sites_reproduce_seed(self, sites, seed):
        elements = []
        positions = []
        for site in sites:
            for position in site.parent_positions():
                elements.append(site.element)
                positions.append(position)
        self.assert_positions_match(elements, positions, seed)
        for site in sites:
            representative = np.asarray(site.representative_parent_fractional)
            np.testing.assert_allclose(
                site.parent_positions()[0],
                representative,
                atol=1e-12,
            )

    def make_c2c_surface_unitcell(self, seed, tol=1e-3):
        spec = CTRsymmetry.SurfaceCellSpec(
            seed.lattice.abc,
            seed.lattice.angles,
            self.C2C_TO_SURFACE,
            translation_range=2,
        )
        return CTRsymmetry.surface_unitcell_from_seed(seed, spec, "LNO", tol=tol)

    def assert_unitcell_matches_supercell(self, unitcell, seed):
        supercell = seed.copy()
        supercell.make_supercell(self.C2C_TO_SURFACE.T)
        self.assert_positions_match(
            unitcell.names,
            unitcell.basis[:, 1:4],
            supercell,
        )

    def test_c2c_sites_are_mapped_back_to_seed_setting(self):
        seed = self.make_c2c_seed()

        for tol in (1e-3, 1e-4, 1e-5):
            with self.subTest(tol=tol):
                sites, number, _ = CTRsymmetry.sites_from_seed(seed, tol=tol)

                self.assertEqual(number, 15)
                self.assertEqual(
                    [site.site_id for site in sites],
                    ["La_4e", "Ni_4c", "O_4e", "O_8f"],
                )
                self.assert_sites_reproduce_seed(sites, seed)
                general = sites[3]
                self.assertEqual(general.variables.keys(), {"u", "v", "w"})
                np.testing.assert_allclose(
                    [general.variables[name] for name in ("u", "v", "w")],
                    general.representative_parent_fractional,
                    atol=1e-12,
                )

    def test_c2c_surface_cell_matches_pymatgen_supercell(self):
        seed = self.make_c2c_seed()

        unitcell = self.make_c2c_surface_unitcell(seed)

        self.assertEqual(len(unitcell.basis), 40)
        self.assert_unitcell_matches_supercell(unitcell, seed)

    def test_origin_shifted_seed_is_reproduced(self):
        from pymatgen.core import Lattice, Structure

        u = 0.30569
        shift = np.asarray([0.1, 0.2, 0.3])
        seed = Structure(
            Lattice.tetragonal(4.653255, 2.9692),
            ["Ru", "Ru", "O", "O", "O", "O"],
            np.asarray(
                [
                    [0.0, 0.0, 0.0],
                    [0.5, 0.5, 0.5],
                    [u, u, 0.0],
                    [1.0 - u, 1.0 - u, 0.0],
                    [0.5 - u, 0.5 + u, 0.5],
                    [0.5 + u, 0.5 - u, 0.5],
                ]
            )
            + shift,
        )

        sites, number, _ = CTRsymmetry.sites_from_seed(seed, tol=1e-3)

        self.assertEqual(number, 136)
        self.assert_sites_reproduce_seed(sites, seed)
        np.testing.assert_allclose(
            sites[1].representative_parent_fractional,
            [u + 0.1, u + 0.2, 0.3],
            atol=1e-12,
        )

    def test_set_wyckoff_site_parameter_updates_variables(self):
        seed = self.make_c2c_seed()
        unitcell = self.make_c2c_surface_unitcell(seed)
        site = unitcell.wyckoff("O_8f")
        new_representative = np.asarray(site["representative_parent_fractional"])
        new_representative += (0.004, -0.002, 0.003)

        for axis, value in zip("xyz", new_representative):
            unitcell.set_wyckoff_site_parameter("O_8f", axis, value)

        site = unitcell.wyckoff("O_8f")
        np.testing.assert_allclose(
            [site["variables"][name] for name in ("u", "v", "w")],
            new_representative,
            atol=1e-12,
        )
        np.testing.assert_allclose(
            site["representative_parent_fractional"],
            new_representative,
            atol=1e-12,
        )
        # The moved cell equals a cell built from the moved seed.
        moved_seed = self.make_c2c_seed(new_representative)
        self.assert_unitcell_matches_supercell(unitcell, moved_seed)

        # Absolute Wyckoff parameters start from the updated reference values.
        basis = unitcell.basis.copy()
        unitcell.addWyckoffParameters("O_8f")
        np.testing.assert_allclose(
            unitcell.getInitialParameters(),
            new_representative,
            atol=1e-12,
        )
        unitcell.setFitParameters(new_representative)
        np.testing.assert_allclose(unitcell.basis, basis, atol=1e-12)

    def test_set_wyckoff_site_parameter_updates_special_site_variable(self):
        seed = self.make_c2c_seed()
        unitcell = self.make_c2c_surface_unitcell(seed)
        site = unitcell.wyckoff("O_4e")
        (variable,) = site["variables"]
        free_axis = "y"
        value = site["representative_parent_fractional"][1] + 0.01

        unitcell.set_wyckoff_site_parameter("O_4e", free_axis, value)

        self.assertAlmostEqual(unitcell.wyckoff("O_4e")["variables"][variable], value)
        # A displacement off the Wyckoff position leaves the variable unchanged.
        unitcell.set_wyckoff_site_parameter("O_4e", "x", 0.01)
        self.assertAlmostEqual(unitcell.wyckoff("O_4e")["variables"][variable], value)


class TestOptionalSymmetryImports(unittest.TestCase):
    def test_missing_pyxtal_reports_optional_dependency(self):
        real_import = __import__

        def import_without_pyxtal(name, *args, **kwargs):
            if name.startswith("pyxtal"):
                raise ImportError("simulated missing pyxtal")
            return real_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=import_without_pyxtal):
            with self.assertRaisesRegex(ImportError, "optional 'symmetry'"):
                CTRsymmetry.possible_wyckoff_positions(136)

    def test_importing_ctr_modules_does_not_import_heavy_symmetry_dependencies(self):
        command = [
            sys.executable,
            "-c",
            (
                "import sys; "
                "import orgui.datautils.xrayutils.CTRcalc; "
                "import orgui.datautils.xrayutils.CTRuc; "
                "import orgui.datautils.xrayutils.CTRsymmetry; "
                "print(any(name in sys.modules for name in "
                "('pyxtal', 'spglib', 'pymatgen', 'ase')))"
            ),
        ]
        completed = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            # A blocked child would otherwise hang the whole job.
            timeout=120,
        )

        self.assertEqual(completed.stdout.strip(), "False")


if __name__ == "__main__":
    unittest.main()

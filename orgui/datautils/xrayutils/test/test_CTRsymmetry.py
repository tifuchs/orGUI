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
        self.assertEqual(sites[1]["variables"], {"u": 0.30569})
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

    def test_wyckoff_parameter_preserves_rutile_u_symmetry(self):
        unitcell = self.make_rutile_110_unitcell()
        couplings = unitcell.wyckoff_couplings("O_4f")
        original = unitcell.basis.copy()

        parameter = unitcell.addWyckoffParameter(
            "O_4f",
            "u",
            absolute_limits=(0.2, 0.4),
        )

        self.assertEqual(parameter.settings["wyckoff"]["value_kind"], "delta")
        np.testing.assert_allclose(unitcell.getInitialParameters(), [0.0])
        np.testing.assert_allclose(
            unitcell.getStartParamAndLimits()[1:],
            ([-0.10569], [0.09431]),
        )
        self.assertEqual(unitcell.wyckoff_sites()[1]["status"], "symmetry_preserving")

        unitcell.setFitParameters([0.01])

        for coupling in couplings:
            column = unitcell.parameterLookup[coupling.coordinate]
            expected = original[coupling.atom_index, column] + 0.01 * coupling.factor
            self.assertAlmostEqual(
                unitcell.basis[coupling.atom_index, column],
                expected,
            )

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

        restored = UnitCell.fromStr(text)

        self.assertEqual(restored.wyckoff_sites()[1]["variables"], {"u": 0.30569})
        self.assertEqual(len(restored.wyckoff_couplings("O_4f")), 8)
        restored.addWyckoffParameter("O_4f", "u", absolute_limits=(0.2, 0.4))
        self.assertEqual(restored.wyckoff_sites()[1]["status"], "symmetry_preserving")


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
        )

        self.assertEqual(completed.stdout.strip(), "False")


if __name__ == "__main__":
    unittest.main()

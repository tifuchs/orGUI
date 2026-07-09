# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
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
# Serialized XPR fixtures intentionally preserve their exact line layout.
# ruff: noqa: E501
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import unittest
from unittest import mock
import os
import tempfile

import numpy as np

from ... import util
from .. import CTRcalc, CTRfilm, CTRplotutil, CTRsymmetry, CTRuc
from ..CTRdistributions import PoissonProfile, SkellamProfile


class TestReadSXRDCrystal(unittest.TestCase):
    xpr_file = """E = 68.00000 keV
# UnitCell relaxations
0000 occupancy = 1.00000 +- nan
return
Coherent 1.00000   1.00000 2.00000 3.00000 4.00000 5.00000 6.00000 7.00000 8.00000 9.00000 10.00000 11.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (0.00249 +- 0.00007)  (0.4871 +- nan)  (0.5663 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (0.51591 +- 0.00010)  (0.7070 +- nan)  (0.8379 +- nan)  (1.0000 +- nan)

# UnitCell adsorbates
0001 occupancy = 1.00000 +- nan
return
Coherent 0.50000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
Coherent 0.50000   0.00000 -1.00000 0.00000 1.00000 0.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  O  (0.00000 +- nan)  (0.50000 +- nan)  (1.00100 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (0.0000 +- nan)

# UnitCell bulk
return
Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (-1.00000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (-0.50000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
"""

    xpr_file_orig = """E = 68.00000 keV
# UnitCell relaxations
0000 occupancy = 1.00000 +- nan
return
Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (0.00249 +- 0.00007)  (0.4871 +- nan)  (0.5663 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (0.51591 +- 0.00010)  (0.7070 +- nan)  (0.8379 +- nan)  (1.0000 +- nan)

# UnitCell adsorbates
0001 occupancy = 1.00000 +- nan
return
Coherent 0.50000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
Coherent 0.50000   0.00000 -1.00000 0.00000 1.00000 0.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  O  (0.00000 +- nan)  (0.50000 +- nan)  (1.00100 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (0.0000 +- nan)

# UnitCell bulk
return
Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (-1.00000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (-0.50000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
"""

    def testFromStr(self):
        xtal = CTRcalc.SXRDCrystal.fromStr(TestReadSXRDCrystal.xpr_file)
        self.assertIsInstance(xtal["relaxations"], CTRcalc.UnitCell)
        with self.assertRaises(ValueError):
            xtal["notexisting"]
        serialized = str(xtal)
        restored = CTRcalc.SXRDCrystal.fromStr(serialized)
        self.assertEqual(str(restored), serialized)

        self.assertTrue(
            np.array_equal(
                xtal["relaxations"].coherentDomainMatrix[0],
                np.array(
                    [
                        [1.00000, 2.00000, 3.00000, 10.00000],
                        [4.00000, 5.00000, 6.00000, 11.00000],
                        [7.00000, 8.00000, 9.00000, 0.00000],
                    ]
                ),
            )
        )

    def testFromFile(self):
        fp = os.path.split(__file__)[0]
        xtal = CTRcalc.SXRDCrystal.fromFile(
            os.path.join(fp, "testdata/0V12_calculated.xpr")
        )
        serialized = str(xtal)
        restored = CTRcalc.SXRDCrystal.fromStr(serialized)
        self.assertEqual(str(restored), serialized)


class TestPoissonSurface(unittest.TestCase):
    def setUp(self):
        self.unitcell = CTRcalc.UnitCell(
            [3.0, 3.0, 4.0],
            [90.0, 90.0, 90.0],
            name="layered",
        )
        self.unitcell.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=0)
        self.unitcell.addAtom("C", [0.0, 0.0, 0.5], 0.1, 0.1, 1.0, layer=1)
        self.unitcell.coherentDomainOccupancy = [0.75]

    def test_serialization_round_trip(self):
        surface = CTRfilm.PoissonSurface(self.unitcell)
        surface.basis = np.array([4.0, 1.0])
        surface.basis_0 = np.copy(surface.basis)
        surface.errors = np.array([0.2, 0.1])

        restored = CTRfilm.PoissonSurface.fromStr(surface.toStr())

        np.testing.assert_allclose(restored.basis, surface.basis)
        np.testing.assert_allclose(restored.errors, surface.errors)
        np.testing.assert_allclose(
            restored.unitcell.basis[:, 1:],
            self.unitcell.basis[:, 1:],
        )
        self.assertEqual(restored.unitcell.names, self.unitcell.names)
        self.assertEqual(restored.unitcell.name, self.unitcell.name)
        self.assertEqual(restored.toStr(), surface.toStr())

    def test_create_layers_assigns_poisson_occupancies(self):
        surface = CTRfilm.PoissonSurface(self.unitcell)
        surface.basis = np.array([4.0, 1.0])

        surface.set_below(0.0, self.unitcell.a[2])

        expected = 0.75 * np.array(
            [
                1.0,
                1.0 - np.exp(-1.0),
                1.0 - 2.0 * np.exp(-1.0),
            ]
        )
        np.testing.assert_allclose(
            surface.layer_ucs[0].coherentDomainOccupancy,
            expected[[0, 2]],
        )
        np.testing.assert_allclose(
            surface.layer_ucs[1].coherentDomainOccupancy,
            expected[[1]],
        )
        self.assertEqual(
            [len(uc.coherentDomainMatrix) for uc in surface.layer_ucs],
            [2, 1],
        )

    def test_profile_support_and_serialization_are_numerical_only(self):
        profile = PoissonProfile(
            mean_change=1.0,
            offset=-1.0,
            tail_probability=1e-8,
        )
        surface = CTRfilm.PoissonSurface(self.unitcell, profile=profile)
        surface.basis_0[:] = surface.basis
        surface.set_below(7.0, 11.0)

        self.assertEqual(surface.pos_absolute, 11.0)
        self.assertEqual(surface.stacking_height_absolute, 11.0)
        self.assertGreater(surface.height_absolute, 11.0)
        self.assertLessEqual(
            surface.layer_ucs[-1].coherentDomainOccupancy[-1],
            self.unitcell.coherentDomainOccupancy[0],
        )

        restored = CTRfilm.PoissonSurface.fromStr(surface.toStr())
        self.assertFalse(restored._legacy_absolute_width)
        self.assertAlmostEqual(restored.profile.mean_change, 1.0)
        self.assertAlmostEqual(restored.profile.offset, -1.0)
        self.assertAlmostEqual(restored.profile.tail_probability, 1e-8)
        self.assertEqual(restored.toStr(), surface.toStr())

    def test_distribution_support_bounds_omitted_tail(self):
        poisson_profile = PoissonProfile(2.5, tail_probability=1e-9)
        _, upper = poisson_profile.support()
        self.assertLessEqual(poisson_profile.occupancy([upper])[0], 1e-9)

        skellam_profile = SkellamProfile(0.35, asymmetry=0.0, tail_probability=1e-8)
        lower, upper = skellam_profile.support(layer_count=2)
        self.assertLess(lower, 0)
        self.assertGreater(upper, 0)

    def test_signed_poisson_process_changes_mean_surface_height(self):
        for width, offset in (
            (2.0, 0.0),
            (-2.0, 0.0),
            (2.0, -2.0),
            (-2.0, 2.0),
        ):
            profile = PoissonProfile(
                mean_change=width,
                offset=offset,
                tail_probability=1e-12,
            )
            lower, upper = profile.support()
            layers = np.arange(lower, upper + 1)
            self.assertAlmostEqual(
                np.sum(profile.correction(layers)),
                width + offset,
                places=10,
            )

    def test_poisson_surface_growth_and_etching_are_signed(self):
        growth = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=1.0),
        )
        growth.set_below(0.0, 8.0)
        growth_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in growth.layer_ucs]
        )
        self.assertTrue(np.all(growth_occupancy > 0))
        self.assertAlmostEqual(growth.mean_height_absolute, 10.0)

        etching = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=-1.0),
        )
        etching.set_below(0.0, 8.0)
        etching_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in etching.layer_ucs]
        )
        self.assertTrue(np.all(etching_occupancy < 0))
        self.assertAlmostEqual(etching.mean_height_absolute, 6.0)

        compensated = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=1.0, offset=-1.0),
        )
        compensated.set_below(0.0, 8.0)
        compensated_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in compensated.layer_ucs]
        )
        self.assertTrue(np.any(compensated_occupancy < 0))
        self.assertTrue(np.any(compensated_occupancy > 0))
        self.assertAlmostEqual(compensated.mean_height_absolute, 8.0)

    def test_zero_width_poisson_surface_stacks_without_domains(self):
        surface = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=0.0),
        )
        crystal = CTRcalc.SXRDCrystal(self.unitcell, surface, stacking=np.array([1]))

        crystal.apply_stacking()

        self.assertEqual(surface.stacking_height_absolute, 0.0)
        self.assertTrue(
            all(not layer.coherentDomainMatrix for layer in surface.layer_ucs)
        )
        h = np.zeros(3)
        ell = np.linspace(0.1, 0.3, 3)
        np.testing.assert_allclose(
            surface.F_uc(h, h, ell),
            0.0,
        )
        self.assertTrue(np.all(np.isfinite(crystal.F(h, h, ell))))
        np.testing.assert_allclose(
            surface.zDensity_G(np.linspace(-1.0, 1.0, 5), 0, 0),
            0.0,
        )

    def test_nominal_delta_width_serialization_migrates_to_offset(self):
        transitional = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=1.0),
        ).toStr()
        transitional = transitional.replace(
            "W/layers offset/layers",
            "Width/layers deltaW/layers",
        ).replace(
            "1.00000    0.00000",
            "0.00000    1.00000",
        )

        restored = CTRfilm.PoissonSurface.fromStr(transitional)

        np.testing.assert_allclose(restored.basis, [1.0, -1.0])
        self.assertAlmostEqual(restored.profile.expected_height_change, 0.0)
        self.assertIn("W/layers offset/layers", restored.toStr())


class TestLayerStacking(unittest.TestCase):
    @staticmethod
    def make_layered_unitcell(name="layered"):
        unitcell = CTRcalc.UnitCell(
            [3.0, 3.0, 6.0],
            [90.0, 90.0, 90.0],
            name=name,
        )
        for layer, z in zip((1, 2, 3), (0.0, 1 / 3, 2 / 3)):
            unitcell.addAtom("C", [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=layer)
            unitcell.layerpos[float(layer)] = z
        return unitcell

    def test_unitcell_selects_layer_after_below_layer(self):
        unitcell = self.make_layered_unitcell()
        unitcell.layer_behavior = "select"

        unitcell.start_layer_number = 3

        self.assertEqual(unitcell.start_layer_number, 1)
        self.assertEqual(unitcell.end_layer_number, 1)
        self.assertEqual(unitcell.layer_behaviour, "select")

    def test_unitcell_translate_layered_cycles_top_to_bottom(self):
        unitcell = self.make_layered_unitcell()
        unitcell.basis[:, 1] = [0.1, 0.2, 0.3]
        unitcell.basis[:, 3] += 0.05
        unitcell.basis_0[:] = unitcell.basis

        transformed = unitcell.translate_layered([0, 0, 1])

        np.testing.assert_allclose(transformed.basis[:, 1], [0.3, 0.1, 0.2])
        np.testing.assert_allclose(
            transformed.basis[:, 3],
            [1.0 + 0.05, 1 / 3 + 0.05, 2 / 3 + 0.05],
        )
        np.testing.assert_array_equal(transformed.basis[:, 7], [1, 2, 3])
        self.assertEqual(
            transformed.layerpos,
            {0.0: 0.0, 1.0: 1.0, 2.0: 1 / 3, 3.0: 2 / 3},
        )
        self.assertTrue(transformed.layer_cycle.is_rotation_of(unitcell.layer_cycle))

    def test_unitcell_translate_layered_cycles_bottom_to_top(self):
        unitcell = self.make_layered_unitcell()
        unitcell.basis[:, 1] = [0.1, 0.2, 0.3]
        unitcell.basis[:, 3] += 0.05
        unitcell.basis_0[:] = unitcell.basis

        transformed = unitcell.translate_layered([0, 0, -1])

        np.testing.assert_allclose(transformed.basis[:, 1], [0.2, 0.3, 0.1])
        np.testing.assert_allclose(
            transformed.basis[:, 3],
            [0.05, 1 / 3 + 0.05, -1 / 3 + 0.05],
        )
        np.testing.assert_array_equal(transformed.basis[:, 7], [1, 2, 3])

    def test_unitcell_translate_layered_preserves_noop_metadata(self):
        unitcell = self.make_layered_unitcell()

        transformed = unitcell.translate_layered([0, 0, 0])

        self.assertIsNot(transformed, unitcell)
        np.testing.assert_allclose(transformed.basis, unitcell.basis)
        np.testing.assert_allclose(transformed.basis_0, unitcell.basis_0)
        self.assertEqual(transformed.layerpos, unitcell.layerpos)
        self.assertEqual(transformed.toStr(), unitcell.toStr())

    def test_unitcell_translate_layered_translates_full_two_layer_cycle(self):
        unitcell = CTRcalc.UnitCell(
            [3.0, 3.0, 6.0],
            [90.0, 90.0, 90.0],
            name="two_layer",
        )
        for layer, z in zip((1, 2), (0.0, 0.5)):
            unitcell.addAtom("C", [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=layer)
            unitcell.layerpos[float(layer)] = z

        transformed = unitcell.translate_layered([0, 0, -2])

        np.testing.assert_allclose(transformed.basis[:, 3], [-1.0, -0.5])
        np.testing.assert_allclose(transformed.basis_0[:, 3], [-1.0, -0.5])
        np.testing.assert_array_equal(transformed.basis[:, 7], unitcell.basis[:, 7])
        self.assertEqual(
            transformed.layerpos,
            {0.0: 0.0, 1.0: -1.0, 2.0: -0.5},
        )

    def test_unitcell_translate_layered_preserves_full_cycle_symmetry_metadata(self):
        unitcell = CTRcalc.UnitCell(
            [3.0, 3.0, 6.0],
            [90.0, 90.0, 90.0],
            name="two_layer_reversed",
        )
        unitcell.addAtom("C", [0.0, 0.0, 0.5], 0.1, 0.1, 1.0, layer=2)
        unitcell.layerpos[2.0] = 0.5
        unitcell.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
        unitcell.layerpos[1.0] = 0.0
        unitcell.symmetry_metadata = CTRsymmetry.SurfaceSymmetryModel(
            CTRsymmetry.SurfaceCellSpec(
                (3.0, 3.0, 6.0),
                (90.0, 90.0, 90.0),
                np.identity(3),
            ),
            (),
            atoms=[
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=0,
                    element="C",
                    site_id="C_0",
                    wyckoff_label="1a",
                    parent_fractional=np.array([0.0, 0.0, 0.5]),
                    surface_fractional=np.array([0.0, 0.0, 0.5]),
                    layer=2,
                    couplings=(
                        CTRsymmetry.WyckoffCoupling(
                            atom_index=0,
                            coordinate="z",
                            variable="u",
                            constant=0.5,
                            factor=1.0,
                            site_id="C_0",
                        ),
                    ),
                ),
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=1,
                    element="C",
                    site_id="C_1",
                    wyckoff_label="1a",
                    parent_fractional=np.array([0.0, 0.0, 0.0]),
                    surface_fractional=np.array([0.0, 0.0, 0.0]),
                    layer=1,
                ),
            ],
        )

        transformed = unitcell.translate_layered([0, 0, -2])

        np.testing.assert_array_equal(transformed.basis[:, 7], [2, 1])
        np.testing.assert_allclose(transformed.basis[:, 3], [-0.5, -1.0])
        self.assertEqual(transformed.layerpos[2.0], -0.5)
        self.assertEqual(transformed.layerpos[1.0], -1.0)
        self.assertIsNotNone(transformed.symmetry_metadata)
        atom = transformed.atom_wyckoff_metadata(0)
        self.assertEqual(atom.layer, 2)
        np.testing.assert_allclose(atom.surface_fractional, [0.0, 0.0, -0.5])
        self.assertEqual(atom.couplings[0].atom_index, 0)
        self.assertAlmostEqual(atom.couplings[0].constant, -0.5)

    def test_unitcell_translate_layered_translates_in_plane(self):
        unitcell = self.make_layered_unitcell()
        original_basis = np.copy(unitcell.basis)
        affine = np.array(
            [
                [1.0, 0.0, 0.0, 2.0],
                [0.0, 1.0, 0.0, -1.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        )
        transformed = unitcell.translate_layered(affine, name="shifted")

        np.testing.assert_allclose(transformed.basis[:, 1], original_basis[:, 1] + 2)
        np.testing.assert_allclose(transformed.basis[:, 2], original_basis[:, 2] - 1)
        np.testing.assert_array_equal(transformed.basis[:, 7], original_basis[:, 7])
        np.testing.assert_allclose(unitcell.basis, original_basis)
        self.assertEqual(transformed.name, "shifted")

    def test_unitcell_translate_layered_remaps_fit_parameter_atoms(self):
        unitcell = self.make_layered_unitcell()
        unitcell.addFitParameter(([2], "x"), name="top_x")

        transformed = unitcell.translate_layered([0, 0, 1])

        atoms, parameters = transformed.parameters["absolute"][0].indices
        np.testing.assert_array_equal(atoms, [0])
        np.testing.assert_array_equal(parameters, [1])

    def test_unitcell_translate_layered_remaps_symmetry_metadata_atoms(self):
        unitcell = self.make_layered_unitcell()
        unitcell.symmetry_metadata = CTRsymmetry.SurfaceSymmetryModel(
            CTRsymmetry.SurfaceCellSpec(
                (3.0, 3.0, 6.0),
                (90.0, 90.0, 90.0),
                np.identity(3),
            ),
            (),
            atoms=[
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=2,
                    element="C",
                    site_id="C_2",
                    wyckoff_label="1a",
                    parent_fractional=np.array([0.0, 0.0, 2 / 3]),
                    surface_fractional=np.array([0.0, 0.0, 2 / 3]),
                    layer=3,
                    couplings=(
                        CTRsymmetry.WyckoffCoupling(
                            atom_index=2,
                            coordinate="z",
                            variable="u",
                            constant=2 / 3,
                            factor=1.0,
                            site_id="C_2",
                        ),
                    ),
                ),
            ],
        )

        transformed = unitcell.translate_layered([0, 0, 1])

        atom = transformed.atom_wyckoff_metadata(0)
        self.assertIsNotNone(atom)
        self.assertEqual(atom.layer, 1)
        np.testing.assert_allclose(atom.surface_fractional, [0.0, 0.0, 1.0])
        self.assertEqual(atom.couplings[0].atom_index, 0)
        self.assertAlmostEqual(atom.couplings[0].constant, 1.0)

    def test_unitcell_translate_layered_rejects_unsupported_affines(self):
        unitcell = self.make_layered_unitcell()

        with self.assertRaisesRegex(ValueError, "x and y"):
            unitcell.translate_layered([0.5, 0, 0])
        with self.assertRaisesRegex(ValueError, "linear part"):
            unitcell.translate_layered(
                np.array(
                    [
                        [0.0, -1.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 0.0],
                    ]
                )
            )

    def test_unitcell_translate_layered_stacks_with_existing_layers(self):
        transformed = self.make_layered_unitcell().translate_layered([0, 0, 1])
        film = CTRfilm.Film(transformed)
        film.basis[0] = 2
        film.start_layer_number = 1
        film.createLayers()

        self.assertEqual(film.start_layer_number, 2)
        self.assertEqual(film.end_layer_number, 3)
        np.testing.assert_array_equal(film.layer_order, [2, 3, 1])

    def test_unitcell_affine_layer_transform_wraps_layers_to_unit_cell(self):
        unitcell = self.make_layered_unitcell()
        unitcell.basis[:, 1] = [0.1, 0.2, 0.3]
        unitcell.basis[:, 3] += 0.05
        unitcell.basis_0[:] = unitcell.basis

        transformed = unitcell.affine_layer_transform([0, 0, 1])

        np.testing.assert_allclose(transformed.basis[:, 1], [0.3, 0.1, 0.2])
        np.testing.assert_allclose(
            transformed.basis[:, 3],
            [0.05, 1 / 3 + 0.05, 2 / 3 + 0.05],
        )
        np.testing.assert_array_equal(transformed.basis[:, 7], [1, 2, 3])
        self.assertEqual(transformed.layerpos, unitcell.layerpos)

    def test_unitcell_affine_layer_transform_full_cycle_preserves_wrapped_z(self):
        unitcell = CTRcalc.UnitCell(
            [3.0, 3.0, 6.0],
            [90.0, 90.0, 90.0],
            name="two_layer",
        )
        for layer, z in zip((1, 2), (0.0, 0.5)):
            unitcell.addAtom("C", [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=layer)
            unitcell.layerpos[float(layer)] = z

        transformed = unitcell.affine_layer_transform([0, 0, -2])

        np.testing.assert_allclose(transformed.basis, unitcell.basis)
        np.testing.assert_allclose(transformed.basis_0, unitcell.basis_0)
        self.assertEqual(transformed.layerpos, unitcell.layerpos)

    def test_unitcell_affine_layer_transform_remaps_symmetry_metadata_atoms(self):
        unitcell = self.make_layered_unitcell()
        unitcell.symmetry_metadata = CTRsymmetry.SurfaceSymmetryModel(
            CTRsymmetry.SurfaceCellSpec(
                (3.0, 3.0, 6.0),
                (90.0, 90.0, 90.0),
                np.identity(3),
            ),
            (),
            atoms=[
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=2,
                    element="C",
                    site_id="C_2",
                    wyckoff_label="1a",
                    parent_fractional=np.array([0.0, 0.0, 2 / 3]),
                    surface_fractional=np.array([0.0, 0.0, 2 / 3]),
                    layer=3,
                    couplings=(
                        CTRsymmetry.WyckoffCoupling(
                            atom_index=2,
                            coordinate="z",
                            variable="u",
                            constant=2 / 3,
                            factor=1.0,
                            site_id="C_2",
                        ),
                    ),
                ),
            ],
        )

        transformed = unitcell.affine_layer_transform([0, 0, 1])

        atom = transformed.atom_wyckoff_metadata(0)
        self.assertIsNotNone(atom)
        self.assertEqual(atom.layer, 1)
        np.testing.assert_allclose(atom.surface_fractional, [0.0, 0.0, 0.0])
        np.testing.assert_allclose(atom.parent_fractional, [0.0, 0.0, 0.0])
        self.assertEqual(atom.couplings[0].atom_index, 0)
        self.assertAlmostEqual(atom.couplings[0].constant, 0.0)

    @staticmethod
    def make_single_wyckoff_cell():
        unitcell = CTRcalc.UnitCell(
            [3.0, 4.0, 5.0],
            [90.0, 90.0, 90.0],
            name="single",
        )
        unitcell.addAtom("C", [0.25, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
        unitcell.layerpos[1.0] = 0.0
        unitcell.symmetry_metadata = CTRsymmetry.SurfaceSymmetryModel(
            CTRsymmetry.SurfaceCellSpec(
                (3.0, 4.0, 5.0),
                (90.0, 90.0, 90.0),
                np.identity(3),
                layer_origins=(0.0,),
            ),
            (
                CTRsymmetry.WyckoffSiteSpec(
                    site_id="C_1",
                    element="C",
                    wyckoff_label="1a",
                    coordinates=(),
                    representative_parent_fractional=(0.25, 0.0, 0.0),
                    variables={"u": 0.25},
                    occ=1.0,
                    iDW=0.1,
                    oDW=0.1,
                ),
            ),
            atoms=[
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=0,
                    element="C",
                    site_id="C_1",
                    wyckoff_label="1a",
                    parent_fractional=np.array([0.25, 0.0, 0.0]),
                    surface_fractional=np.array([0.25, 0.0, 0.0]),
                    layer=1,
                    couplings=(
                        CTRsymmetry.WyckoffCoupling(
                            atom_index=0,
                            coordinate="x",
                            variable="u",
                            constant=0.0,
                            factor=1.0,
                            site_id="C_1",
                        ),
                    ),
                    site_couplings=(
                        CTRsymmetry.WyckoffSiteCoupling(
                            atom_index=0,
                            coordinate="x",
                            axis="x",
                            factor=1.0,
                            site_id="C_1",
                        ),
                    ),
                ),
            ],
        )
        return unitcell

    def test_unitcell_supercell_scales_lattice_and_coordinates(self):
        unitcell = self.make_single_wyckoff_cell()

        supercell = unitcell.supercell((2, 3, 1), name="super")

        np.testing.assert_allclose(supercell.a, [6.0, 12.0, 5.0])
        self.assertEqual(supercell.name, "super")
        self.assertEqual(supercell.basis.shape[0], 6)
        np.testing.assert_allclose(
            supercell.basis[:2, 1:4],
            [[0.125, 0.0, 0.0], [0.625, 0.0, 0.0]],
        )
        np.testing.assert_allclose(
            supercell.symmetry_metadata.surface_spec.transform,
            np.diag([2.0, 3.0, 1.0]),
        )
        self.assertEqual(supercell.parameters, {"absolute": [], "relative": []})

    def test_unitcell_supercell_repeats_layers_along_z(self):
        unitcell = self.make_layered_unitcell()

        supercell = unitcell.supercell((1, 1, 2))

        np.testing.assert_array_equal(supercell.layers, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        self.assertEqual(set(supercell.layerpos), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0})
        np.testing.assert_allclose(
            [supercell.layerpos[layer] for layer in supercell.layers],
            [0.0, 1 / 6, 1 / 3, 0.5, 2 / 3, 5 / 6],
        )

    def test_unitcell_supercell_preserves_shared_wyckoff_site(self):
        unitcell = self.make_single_wyckoff_cell()

        supercell = unitcell.supercell((2, 1, 1), symmetry="preserve")

        self.assertEqual(
            [site["site_id"] for site in supercell.wyckoff_sites()],
            ["C_1"],
        )
        couplings = supercell.wyckoff_couplings("C_1")
        self.assertEqual([coupling.atom_index for coupling in couplings], [0, 1])
        np.testing.assert_allclose(
            [coupling.constant for coupling in couplings],
            [0.0, 0.5],
        )
        np.testing.assert_allclose(
            [coupling.factor for coupling in couplings],
            [0.5, 0.5],
        )

        supercell.addWyckoffParameter("C_1", "u", limits=(-0.2, 0.2))
        supercell.setFitParameters([0.1])

        np.testing.assert_allclose(supercell.basis[:, 1], [0.175, 0.675])

        supercell = unitcell.supercell((2, 1, 1), symmetry="preserve")
        supercell.addWyckoffShift("C_1", "x", limits=(-0.2, 0.2))
        supercell.setFitParameters([0.1])

        np.testing.assert_allclose(supercell.basis[:, 1], [0.175, 0.675])

    def test_unitcell_supercell_independent_wyckoff_sites_fit_separately(self):
        unitcell = self.make_single_wyckoff_cell()

        supercell = unitcell.supercell((2, 1, 1), symmetry="independent")

        self.assertEqual(
            [site["site_id"] for site in supercell.wyckoff_sites()],
            ["C_1_copy0", "C_1_copy1"],
        )
        supercell.addWyckoffParameter(
            "C_1_copy0",
            "u",
            limits=(-0.2, 0.2),
        )
        supercell.setFitParameters([0.1])

        np.testing.assert_allclose(supercell.basis[:, 1], [0.175, 0.625])

        supercell = unitcell.supercell((2, 1, 1), symmetry="independent")
        supercell.addWyckoffShift("C_1_copy0", "x", limits=(-0.2, 0.2))
        supercell.setFitParameters([0.1])

        np.testing.assert_allclose(supercell.basis[:, 1], [0.175, 0.625])

    def test_film_forwards_wyckoff_parameters_to_template_unitcell(self):
        film = CTRfilm.Film(self.make_single_wyckoff_cell())

        film.addWyckoffParameter("C_1", "u", limits=(-0.2, 0.2))
        film.setFitParameters([0.1])

        self.assertEqual(film.unitcell.wyckoff_sites()[0]["status"], "symmetry_preserving")
        np.testing.assert_allclose(film.unitcell.basis[:, 1], [0.35])

    def test_poisson_surface_forwards_wyckoff_shifts_to_template_unitcell(self):
        surface = CTRfilm.PoissonSurface(self.make_single_wyckoff_cell())

        surface.addWyckoffShift("C_1", "x", limits=(-0.2, 0.2))
        surface.setFitParameters([0.1])

        self.assertEqual(surface.unitcell.wyckoff_sites()[0]["status"], "site_displaced")
        np.testing.assert_allclose(surface.unitcell.basis[:, 1], [0.35])

    def test_epitaxy_interface_uses_unitcell_selector_for_wyckoff_parameters(self):
        interface = CTRfilm.EpitaxyInterface(
            self.make_single_wyckoff_cell(),
            self.make_single_wyckoff_cell(),
        )

        with self.assertRaisesRegex(ValueError, "Missing unit cell name"):
            interface.addWyckoffParameter("C_1", "u", limits=(-0.2, 0.2))

        interface.addWyckoffShift(
            "C_1",
            "x",
            limits=(-0.2, 0.2),
            unitcell="top",
        )
        interface.setFitParameters([0.1])

        self.assertEqual(interface.uc_top.wyckoff_sites()[0]["status"], "site_displaced")
        self.assertEqual(interface.uc_bottom.wyckoff_sites()[0]["status"], "metadata_only")
        np.testing.assert_allclose(interface.uc_top.basis[:, 1], [0.35])
        np.testing.assert_allclose(interface.uc_bottom.basis[:, 1], [0.25])

    def test_epitaxy_interface_accepts_unitcell_list_for_wyckoff_parameters(self):
        interface = CTRfilm.EpitaxyInterface(
            self.make_single_wyckoff_cell(),
            self.make_single_wyckoff_cell(),
        )

        parameters = interface.addWyckoffParameter(
            "C_1",
            "u",
            limits=(-0.2, 0.2),
            unitcell=["top", "bottom"],
        )
        interface.setFitParameters([0.1, 0.1])

        self.assertEqual(len(parameters), 2)
        self.assertEqual(interface.uc_top.wyckoff_sites()[0]["status"], "symmetry_preserving")
        self.assertEqual(interface.uc_bottom.wyckoff_sites()[0]["status"], "symmetry_preserving")
        np.testing.assert_allclose(interface.uc_top.basis[:, 1], [0.35])
        np.testing.assert_allclose(interface.uc_bottom.basis[:, 1], [0.35])

    def test_sxrdcrystal_links_wyckoff_parameters_across_unitcells(self):
        bulk = self.make_single_wyckoff_cell()
        film1 = self.make_single_wyckoff_cell()
        film2 = self.make_single_wyckoff_cell()
        bulk.name = "bulk_template"
        film1.name = "film1"
        film2.name = "film2"
        crystal = CTRcalc.SXRDCrystal(bulk, film1, film2)

        parameter = crystal.addWyckoffParameter(
            {
                "film1": ("C_1", "u"),
                "film2": ("C_1", "u"),
            },
            name="shared_u",
            limits=(-0.2, 0.2),
        )
        crystal.setParameters([0.1])

        self.assertEqual(parameter.name, "shared_u")
        self.assertEqual(crystal.fitparnames, ["shared_u"])
        self.assertEqual(film1.wyckoff_sites()[0]["status"], "symmetry_preserving")
        self.assertEqual(film2.wyckoff_sites()[0]["status"], "symmetry_preserving")
        np.testing.assert_allclose(film1.basis[:, 1], [0.35])
        np.testing.assert_allclose(film2.basis[:, 1], [0.35])

    def test_sxrdcrystal_links_wyckoff_shifts_across_unitcells(self):
        bulk = self.make_single_wyckoff_cell()
        film1 = self.make_single_wyckoff_cell()
        film2 = self.make_single_wyckoff_cell()
        bulk.name = "bulk_template"
        film1.name = "film1"
        film2.name = "film2"
        crystal = CTRcalc.SXRDCrystal(bulk, film1, film2)

        crystal.addWyckoffShift(
            {
                "film1": ("C_1", "x"),
                "film2": ("C_1", "x"),
            },
            name="shared_shift",
            limits=(-0.2, 0.2),
        )
        crystal.setParameters([0.1])

        self.assertEqual(crystal.fitparnames, ["shared_shift"])
        self.assertEqual(film1.wyckoff_sites()[0]["status"], "site_displaced")
        self.assertEqual(film2.wyckoff_sites()[0]["status"], "site_displaced")
        np.testing.assert_allclose(film1.basis[:, 1], [0.35])
        np.testing.assert_allclose(film2.basis[:, 1], [0.35])

    def test_sxrdcrystal_links_wyckoff_parameters_inside_interface(self):
        bulk = self.make_single_wyckoff_cell()
        top = self.make_single_wyckoff_cell()
        bottom = self.make_single_wyckoff_cell()
        bulk.name = "bulk_template"
        interface = CTRfilm.EpitaxyInterface(top, bottom, name="interface")
        crystal = CTRcalc.SXRDCrystal(bulk, interface)

        with self.assertRaisesRegex(ValueError, "EpitaxyInterface"):
            crystal.addWyckoffParameter({"interface": ("C_1", "u")})

        crystal.addWyckoffParameter(
            {
                "interface": {
                    "site_id": "C_1",
                    "variable": "u",
                    "unitcell": ("top", "bottom"),
                }
            },
            name="shared_interface_u",
            limits=(-0.2, 0.2),
        )
        crystal.setParameters([0.1])

        self.assertEqual(crystal.fitparnames, ["shared_interface_u"])
        self.assertEqual(interface.uc_top.wyckoff_sites()[0]["status"], "symmetry_preserving")
        self.assertEqual(interface.uc_bottom.wyckoff_sites()[0]["status"], "symmetry_preserving")
        np.testing.assert_allclose(interface.uc_top.basis[:, 1], [0.35])
        np.testing.assert_allclose(interface.uc_bottom.basis[:, 1], [0.35])

    def test_sxrdcrystal_links_wyckoff_shift_to_interface_unitcell_key(self):
        bulk = self.make_single_wyckoff_cell()
        top = self.make_single_wyckoff_cell()
        bottom = self.make_single_wyckoff_cell()
        bulk.name = "bulk_template"
        interface = CTRfilm.EpitaxyInterface(top, bottom, name="interface")
        crystal = CTRcalc.SXRDCrystal(bulk, interface)

        crystal.addWyckoffShift(
            {("interface", "top"): ("C_1", "x")},
            name="top_shift",
            limits=(-0.2, 0.2),
        )
        crystal.setParameters([0.1])

        self.assertEqual(crystal.fitparnames, ["top_shift"])
        self.assertEqual(interface.uc_top.wyckoff_sites()[0]["status"], "site_displaced")
        self.assertEqual(interface.uc_bottom.wyckoff_sites()[0]["status"], "metadata_only")
        np.testing.assert_allclose(interface.uc_top.basis[:, 1], [0.35])
        np.testing.assert_allclose(interface.uc_bottom.basis[:, 1], [0.25])

    def test_film_and_poisson_rotate_cyclic_layer_order(self):
        film = CTRfilm.Film(self.make_layered_unitcell("film"))
        film.basis[0] = 2
        film.start_layer_number = 1
        film.createLayers()

        self.assertEqual(film.start_layer_number, 2)
        self.assertEqual(film.end_layer_number, 3)
        np.testing.assert_array_equal(film.layer_order, [2, 3, 1])
        self.assertEqual(
            [len(layer.coherentDomainMatrix) for layer in film.layer_ucs],
            [1, 1, 0],
        )
        self.assertAlmostEqual(film.height_absolute, 4.0)

        surface = CTRfilm.PoissonSurface(self.make_layered_unitcell("surface"))
        surface.basis[:] = [2, 0]
        surface.start_layer_number = 2
        surface.createLayers()

        self.assertEqual(surface.start_layer_number, 3)
        self.assertEqual(surface.end_layer_number, 1)
        np.testing.assert_array_equal(surface.layer_order, [3, 1, 2])
        self.assertAlmostEqual(
            surface.layer_ucs[0].coherentDomainMatrix[0][2, 3],
            -2 / 3,
        )
        self.assertAlmostEqual(
            surface.layer_ucs[1].coherentDomainMatrix[0][2, 3],
            1 / 3,
        )

    def test_epitaxy_interface_rotates_cyclic_layer_order(self):
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_unitcell("top"),
            self.make_layered_unitcell("bottom"),
            fixed_ucs=3,
        )
        interface.basis[:] = [0.5, 0.0]
        interface.stack_on(0.0, 12.0, 1)

        self.assertEqual(interface.start_layer_number, 2)
        self.assertEqual(interface.end_layer_number, 1)
        np.testing.assert_array_equal(interface.layer_order, [2, 3, 1])
        self.assertAlmostEqual(interface.pos_absolute, 12.0)
        self.assertGreater(interface.height_absolute, 12.0)
        self.assertGreater(
            interface.top_layers[-1].coherentDomainMatrix[0][2, 3],
            interface.top_layers[0].coherentDomainMatrix[0][2, 3],
        )

    def test_sxrdcrystal_passes_top_layer_to_object_above(self):
        lower = CTRfilm.Film(self.make_layered_unitcell("lower"))
        lower.name = "lower"
        lower.basis[0] = 2

        upper = CTRfilm.Film(self.make_layered_unitcell("upper"))
        upper.name = "upper"
        upper.basis[0] = 4
        upper.pos_absolute = 100.0

        surface = CTRfilm.PoissonSurface(self.make_layered_unitcell("surface"))
        surface.name = "surface"
        surface.basis[:] = [5, 0]

        crystal = CTRcalc.SXRDCrystal(
            self.make_layered_unitcell("bulk"),
            lower,
            upper,
            surface,
            stacking=np.array([1, 2, 3]),
        )

        crystal.apply_stacking()

        self.assertEqual(lower.start_layer_number, 1)
        self.assertEqual(lower.end_layer_number, 2)
        self.assertEqual(upper.start_layer_number, 3)
        self.assertEqual(upper.end_layer_number, 1)
        self.assertAlmostEqual(upper.pos_absolute, lower.height_absolute)
        self.assertAlmostEqual(
            upper.height_absolute,
            lower.height_absolute + 4.0,
        )
        self.assertAlmostEqual(surface.pos_absolute, upper.height_absolute)
        self.assertGreater(surface.height_absolute, surface.pos_absolute)

    def test_stacked_objects_default_to_full_crystal_weight(self):
        film = CTRfilm.Film(self.make_layered_unitcell("film"))
        film.basis[0] = 3
        surface = CTRfilm.PoissonSurface(
            self.make_layered_unitcell("surface"),
            profile=PoissonProfile(1.0),
        )

        stacked = CTRcalc.SXRDCrystal(
            self.make_layered_unitcell("bulk"),
            film,
            surface,
            stacking=np.array([1, 2]),
        )
        np.testing.assert_array_equal(stacked.weights, [1.0, 1.0])
        np.testing.assert_array_equal(stacked.weights_0, [1.0, 1.0])

        first = self.make_layered_unitcell("first")
        second = self.make_layered_unitcell("second")
        alternatives = CTRcalc.SXRDCrystal(
            self.make_layered_unitcell("bulk"),
            first,
            second,
        )
        np.testing.assert_allclose(alternatives.weights, [0.5, 0.5])

    def test_different_cycles_require_an_explicit_transition(self):
        bulk = CTRcalc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0], name="bulk")
        bulk.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
        bulk.addAtom("C", [0.0, 0.0, 0.5], 0.1, 0.1, 1.0, layer=2)
        bulk.layerpos = {1.0: 0.0, 2.0: 0.5}

        upper = CTRfilm.Film(self.make_layered_unitcell("upper"))
        upper.basis[0] = 2
        crystal = CTRcalc.SXRDCrystal(bulk, upper, stacking=np.array([1]))

        with self.assertRaisesRegex(ValueError, "layer_transition"):
            crystal.apply_stacking()

    def test_transition_maps_lower_termination_to_upper_start(self):
        bulk = CTRcalc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0], name="bulk")
        bulk.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
        bulk.addAtom("C", [0.0, 0.0, 0.5], 0.1, 0.1, 1.0, layer=2)
        bulk.layerpos = {1.0: 0.0, 2.0: 0.5}

        upper = CTRfilm.Film(
            self.make_layered_unitcell("upper"),
            layer_transition={1: 2, 2: 3},
        )
        upper.basis[0] = 2
        crystal = CTRcalc.SXRDCrystal(bulk, upper, stacking=np.array([1]))
        crystal.apply_stacking()

        self.assertEqual(bulk.end_layer_number, 2)
        self.assertEqual(upper.start_layer_number, 3)
        np.testing.assert_array_equal(upper.layer_order, [3, 1, 2])
        self.assertEqual(upper.end_layer_number, 1)

        restored = CTRfilm.Film.fromStr(upper.toStr())
        self.assertEqual(restored.layer_transition.mapping, {1.0: 2.0, 2.0: 3.0})

    def test_rotated_cycle_is_compatible_without_transition(self):
        lower = self.make_layered_unitcell("lower")
        lower.layer_cycle = (2, 3, 1)
        upper = CTRfilm.Film(self.make_layered_unitcell("upper"))
        upper.basis[0] = 2
        crystal = CTRcalc.SXRDCrystal(lower, upper, stacking=np.array([1]))

        crystal.apply_stacking()

        self.assertEqual(lower.end_layer_number, 3)
        self.assertEqual(upper.start_layer_number, 1)

    def test_profile_interface_does_not_advance_nominal_cursor(self):
        lower = self.make_layered_unitcell("lower")
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_unitcell("top"),
            self.make_layered_unitcell("bottom"),
            profile=SkellamProfile(0.5),
        )
        film = CTRfilm.Film(self.make_layered_unitcell("film"))
        film.basis[0] = 3
        crystal = CTRcalc.SXRDCrystal(
            lower,
            interface,
            film,
            stacking=np.array([1, 2]),
        )

        crystal.apply_stacking()

        self.assertEqual(interface.stacking_height_absolute, 0.0)
        self.assertEqual(film.pos_absolute, 0.0)
        self.assertAlmostEqual(film.height_absolute, 6.0)

        restored = CTRfilm.EpitaxyInterface.fromStr(interface.toStr())
        self.assertFalse(restored._legacy_support_cursor)
        self.assertIsInstance(restored.profile, SkellamProfile)
        np.testing.assert_allclose(restored.basis, interface.basis)

    def test_profile_interface_is_a_correction_not_extra_material(self):
        unitcell = self.make_layered_unitcell("material")
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_unitcell("top"),
            self.make_layered_unitcell("bottom"),
            profile=SkellamProfile(0.5),
        )
        interface.stack_on(
            0.0,
            0.0,
            unitcell.end_layer_number,
            below_state=unitcell.layer_state,
        )

        top_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in interface.top_layers]
        )
        bottom_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in interface.bottom_layers]
        )

        self.assertTrue(np.any(top_occupancy < 0.0))
        self.assertTrue(np.any(bottom_occupancy < 0.0))
        np.testing.assert_allclose(
            np.sort(top_occupancy),
            np.sort(-bottom_occupancy),
            atol=1e-15,
        )


class TestLegacyLayeredCTR(unittest.TestCase):
    def test_legacy_xtal_reconstructs_reference_interface(self):
        repository_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
        )
        fixture_root = os.path.join(
            repository_root, "examples", "CTR", "test"
        )
        xtal_path = os.path.join(
            fixture_root, "0001_fit_2V036_reference.xtal"
        )
        reference_path = os.path.join(fixture_root, "CTRs_reference.dat")
        for path in (xtal_path, reference_path):
            if not os.path.exists(path):
                self.skipTest(f"Missing optional CTR fixture: {path}")
        xtal = CTRcalc.SXRDCrystal.fromFile(
            xtal_path
        )

        self.assertIsInstance(xtal["RuO2"], CTRfilm.Film)
        self.assertIsInstance(xtal["TiO2toRuO2"], CTRfilm.EpitaxyInterface)
        np.testing.assert_allclose(xtal["TiO2toRuO2"].basis, [0.35, 0.0])
        self.assertTrue(xtal["TiO2toRuO2"]._legacy_support_cursor)

        xtal["RuO2"].basis[0] = 17.0
        xtal["RuO2"].basis_0[0] = 17.0
        xtal.atten = 0.01
        reference = np.loadtxt(
            reference_path,
            skiprows=1,
            max_rows=8,
        )
        l_values = np.linspace(0.0, 7.25, 2000)[:8]
        calculated = np.abs(
            xtal.F(
                np.zeros(8, dtype=np.float64),
                np.zeros(8, dtype=np.float64),
                l_values,
            )
        )
        np.testing.assert_allclose(
            calculated,
            reference[:, 3] * xtal.reference_area,
            rtol=2e-5,
            atol=2e-5,
        )


class TestStructureFactorNormalization(unittest.TestCase):
    @staticmethod
    def make_carbon_cell(a, atoms, name):
        cell = CTRcalc.UnitCell(
            [a, 3.0, 4.0],
            [90.0, 90.0, 90.0],
            name=name,
        )
        for x in atoms:
            cell.addAtom("C", [x, 0.0, 0.0], 0.1, 0.1, 1.0)
        cell.setEnergy(10000.0)
        return cell

    def test_raw_unitcell_factor_scales_with_supercell_area(self):
        primitive = self.make_carbon_cell(3.0, [0.0], "primitive")
        supercell = self.make_carbon_cell(6.0, [0.0, 0.5], "supercell")
        h = np.zeros(1)

        primitive_factor = primitive.F_uc(h, h, h)
        supercell_factor = supercell.F_uc(h, h, h)

        np.testing.assert_allclose(supercell_factor, 2.0 * primitive_factor)
        self.assertAlmostEqual(supercell.uc_area, 2.0 * primitive.uc_area)

    def test_crystal_factor_is_invariant_to_lateral_supercell(self):
        bulk = self.make_carbon_cell(3.0, [0.0], "bulk")
        primitive = self.make_carbon_cell(3.0, [0.0], "primitive")
        supercell = self.make_carbon_cell(6.0, [0.0, 0.5], "supercell")
        primitive_crystal = CTRcalc.SXRDCrystal(bulk, primitive)
        supercell_crystal = CTRcalc.SXRDCrystal(bulk, supercell)
        h = np.zeros(1)

        np.testing.assert_allclose(
            primitive_crystal.F(h, h, h),
            supercell_crystal.F(h, h, h),
            rtol=1e-2,
        )

    def test_constructor_propagates_bulk_reference_to_layers(self):
        bulk = self.make_carbon_cell(3.0, [0.0], "bulk")
        film_cell = self.make_carbon_cell(6.0, [0.0, 0.5], "film")
        film = CTRfilm.Film(film_cell)
        film.basis[0] = 1

        crystal = CTRcalc.SXRDCrystal(bulk, film, stacking=np.array([1]))

        self.assertIs(crystal.reference_uc, bulk)
        self.assertEqual(crystal.reference_area, bulk.uc_area)
        np.testing.assert_allclose(
            film.layer_ucs[0].refHKLTransform,
            film.unitcell.refHKLTransform,
        )
        self.assertFalse(
            np.allclose(
                film.layer_ucs[0].refHKLTransform,
                np.identity(3),
            )
        )

    def test_interface_uses_lower_unitcell_area(self):
        lower = self.make_carbon_cell(3.0, [0.0], "lower")
        upper = self.make_carbon_cell(6.0, [0.0, 0.5], "upper")
        interface = CTRfilm.EpitaxyInterface(
            upper,
            lower,
            profile=SkellamProfile(0.25),
        )
        interface.setReferenceUnitCell(lower)
        interface.createInterfaceCells()
        h = np.zeros(1)
        expected = np.zeros(1, dtype=np.complex128)
        for top_layer, bottom_layer in zip(
            interface.top_layers, interface.bottom_layers
        ):
            expected += lower.uc_area * top_layer.F_uc(h, h, h) / top_layer.uc_area
            expected += (
                lower.uc_area * bottom_layer.F_uc(h, h, h) / bottom_layer.uc_area
            )

        self.assertEqual(interface.uc_area, lower.uc_area)
        np.testing.assert_allclose(interface.F_uc(h, h, h), expected)


class TestCTRTextFilesAndFitParameters(unittest.TestCase):
    @staticmethod
    def make_layered_cell(name):
        cell = CTRcalc.UnitCell(
            [3.0, 3.0, 4.0],
            [90.0, 90.0, 90.0],
            name=name,
        )
        cell.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
        cell.addAtom("C", [0.0, 0.0, 0.5], 0.1, 0.1, 1.0, layer=2)
        cell.layerpos = {1.0: 0.0, 2.0: 0.5}
        cell.setEnergy(10000.0)
        return cell

    def make_crystal_with_errors(self):
        bulk = self.make_layered_cell("bulk")
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_cell("interface_top"),
            self.make_layered_cell("interface_bottom"),
            profile=SkellamProfile(0.35, 0.0),
            name="interface",
        )
        film = CTRfilm.Film(self.make_layered_cell("film_cell"), name="film")
        film.basis[:] = [6.0]
        film.basis_0[:] = film.basis
        surface = CTRfilm.PoissonSurface(
            self.make_layered_cell("surface_cell"),
            profile=PoissonProfile(mean_change=1.5, offset=-1.5),
            name="surface",
        )
        surface.basis_0[:] = surface.basis

        film.errors = np.array([0.2])
        interface.errors = np.array([0.03, 0.04])
        surface.errors = np.array([0.1, 0.15])
        film.unitcell.errors = np.full_like(film.unitcell.basis, np.nan)
        film.unitcell.errors[:, 1:7] = 0.01

        crystal = CTRcalc.SXRDCrystal(
            bulk,
            interface,
            film,
            surface,
            stacking=np.array([1, 2, 3]),
        )
        crystal.werrors = np.array([0.01, 0.02, 0.03])
        return crystal

    def test_xtal_and_xpr_plain_text_round_trip(self):
        crystal = self.make_crystal_with_errors()
        with tempfile.TemporaryDirectory() as directory:
            xtal_path = os.path.join(directory, "model.xtal")
            xpr_path = os.path.join(directory, "model.xpr")
            crystal.toFile(xtal_path)
            crystal.toFile(xpr_path)

            with open(xtal_path, encoding="utf-8") as stream:
                xtal_text = stream.read()
            with open(xpr_path, encoding="utf-8") as stream:
                xpr_text = stream.read()

            self.assertNotIn("+-", xtal_text)
            self.assertIn("+-", xpr_text)

            restored_xtal = CTRcalc.SXRDCrystal.fromFile(xtal_path)
            restored_xpr = CTRcalc.SXRDCrystal.fromFile(xpr_path)

        self.assertIsNone(restored_xtal.werrors)
        self.assertIsNone(restored_xtal["film"].errors)
        self.assertIsNone(restored_xtal["interface"].errors)
        self.assertIsNone(restored_xtal["surface"].errors)
        np.testing.assert_allclose(restored_xpr.werrors, crystal.werrors)
        np.testing.assert_allclose(restored_xpr["film"].errors, [0.2])
        np.testing.assert_allclose(restored_xpr["interface"].errors, [0.03, 0.04])
        np.testing.assert_allclose(restored_xpr["surface"].errors, [0.1, 0.15])
        np.testing.assert_allclose(
            restored_xpr["film"].unitcell.errors[:, 1:7],
            crystal["film"].unitcell.errors[:, 1:7],
        )

    def test_ctrfilm_fit_parameters_and_errors(self):
        film = CTRfilm.Film(self.make_layered_cell("film_cell"), name="film")
        film.basis[:] = [6.0]
        film.basis_0[:] = film.basis
        film.addFitParameter("W", limits=(1.0, 10.0))
        np.testing.assert_allclose(film.getStartParamAndLimits()[0], [6.0])
        film.setFitParameters([7.0])
        film.setFitErrors([0.25])
        np.testing.assert_allclose(film.basis, [7.0])
        np.testing.assert_allclose(film.getFitErrors(), [0.25])

        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_cell("top"),
            self.make_layered_cell("bottom"),
            profile=SkellamProfile(0.35, 0.0),
            name="interface",
        )
        interface.basis_0[:] = interface.basis
        interface.addRelParameter(("W", "S"), (1.0, -1.0), limits=(-0.2, 0.2))
        interface.setFitParameters([0.1])
        interface.setFitErrors([0.02])
        np.testing.assert_allclose(interface.basis, [0.45, -0.1])
        np.testing.assert_allclose(interface.errors, [0.02, 0.02])
        self.assertAlmostEqual(interface.profile.width, 0.45)
        self.assertAlmostEqual(interface.profile.asymmetry, -0.1)

        surface = CTRfilm.PoissonSurface(
            self.make_layered_cell("surface_cell"),
            profile=PoissonProfile(mean_change=2.0, offset=-2.0),
            name="surface",
        )
        surface.basis_0[:] = surface.basis
        surface.addRelParameter(
            ("W", "offset"),
            (1.0, -1.0),
            limits=(-1.0, 1.0),
        )
        surface.setFitParameters([0.5])
        surface.setFitErrors([0.1])
        np.testing.assert_allclose(surface.basis, [2.5, -2.5])
        np.testing.assert_allclose(surface.errors, [0.1, 0.1])
        self.assertAlmostEqual(surface.profile.mean_change, 2.5)
        self.assertAlmostEqual(surface.profile.offset, -2.5)
        self.assertAlmostEqual(surface.profile.expected_height_change, 0.0)

    def test_linear_fit_parameter_dictionary_round_trip(self):
        original = CTRfilm.PoissonSurface(
            self.make_layered_cell("surface_cell"),
            profile=PoissonProfile(mean_change=2.0, offset=-2.0),
            name="surface",
        )
        original.basis_0[:] = original.basis
        original.addRelParameter(
            ("W", "offset"),
            (1.0, -1.0),
            limits=(-1.0, 1.0),
            name="mean_preserving_width",
        )
        original.setFitParameters([0.4])
        original.setFitErrors([0.05])

        restored = CTRfilm.PoissonSurface(
            self.make_layered_cell("surface_cell"),
            profile=PoissonProfile(mean_change=2.0, offset=-2.0),
            name="surface",
        )
        restored.parametersFromDict(original.parametersToDict())

        np.testing.assert_allclose(restored.basis, original.basis)
        np.testing.assert_allclose(restored.errors, original.errors)
        np.testing.assert_allclose(restored.getFitErrors(), original.getFitErrors())
        self.assertEqual(restored.fitparnames, original.fitparnames)
        self.assertAlmostEqual(restored.profile.mean_change, 2.4)
        self.assertAlmostEqual(restored.profile.offset, -2.4)

    def test_crystal_fit_vector_updates_ctrfilm_components(self):
        crystal = self.make_crystal_with_errors()
        crystal["interface"].clearParameters()
        crystal["film"].clearParameters()
        crystal["surface"].clearParameters()
        crystal["interface"].addFitParameter("W", limits=(0.0, 1.0))
        crystal["film"].addFitParameter("W", limits=(1.0, 10.0))
        crystal["surface"].addFitParameter("W", limits=(-5.0, 5.0))

        np.testing.assert_allclose(crystal.getInitialParameters(), [0.35, 6.0, 1.5])
        crystal.setParameters([0.4, 7.0, 2.0])
        crystal.setFitErrors([0.01, 0.2, 0.1])

        np.testing.assert_allclose(crystal["interface"].basis, [0.4, 0.0])
        np.testing.assert_allclose(crystal["film"].basis, [7.0])
        np.testing.assert_allclose(crystal["surface"].basis, [2.0, -1.5])
        np.testing.assert_allclose(crystal.getFitErrors(), [0.01, 0.2, 0.1])
        self.assertAlmostEqual(crystal["interface"].profile.width, 0.4)
        self.assertAlmostEqual(crystal["surface"].profile.mean_change, 2.0)

    def test_example_xtal_and_xpr_describe_same_poisson_stack(self):
        repository_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
        )
        example_root = os.path.join(repository_root, "examples", "CTR")
        xtal_path = os.path.join(example_root, "RuO2_TiO2_Poisson_etching.xtal")
        xpr_path = os.path.join(example_root, "RuO2_TiO2_Poisson_etching.xpr")
        for path in (xtal_path, xpr_path):
            if not os.path.exists(path):
                self.skipTest(f"Missing optional CTR fixture: {path}")
        xtal = CTRcalc.SXRDCrystal.fromFile(xtal_path)
        xpr = CTRcalc.SXRDCrystal.fromFile(xpr_path)

        expected_types = [
            CTRfilm.PoissonSurface,
            CTRfilm.Film,
            CTRfilm.EpitaxyInterface,
        ]
        self.assertEqual(
            [type(component) for component in xtal.uc_surface_list],
            expected_types,
        )
        self.assertEqual(
            [type(component) for component in xpr.uc_surface_list],
            expected_types,
        )
        np.testing.assert_array_equal(xtal.uc_stacking, [3, 2, 1])
        np.testing.assert_array_equal(xpr.uc_stacking, [3, 2, 1])
        np.testing.assert_allclose(xtal["TiO2toRuO2"].basis, [0.35, 0.0])
        np.testing.assert_allclose(xtal["RuO2"].basis, [17.0])
        np.testing.assert_allclose(xtal["RuO2surface"].basis, [-6.0, 0.0])
        for name in ("TiO2toRuO2", "RuO2", "RuO2surface"):
            np.testing.assert_allclose(xtal[name].basis, xpr[name].basis)

        l = np.linspace(0.05, 7.0, 200)  # noqa: E741
        zeros = np.zeros_like(l)
        z = np.linspace(-20.0, 100.0, 1000)
        np.testing.assert_allclose(xtal.F(zeros, zeros, l), xpr.F(zeros, zeros, l))
        np.testing.assert_allclose(
            xtal.zDensity_G(z, 0, 0),
            xpr.zDensity_G(z, 0, 0),
        )
        self.assertIsNone(xtal["RuO2"].unitcell.errors)
        self.assertIsNotNone(xpr["RuO2"].unitcell.errors)


class StructureFactorValidationMixin:
    def assert_structure_factors_match_volume_normalized_reference(self):
        """Compare current amplitudes with legacy volume-normalized F."""
        calc_CTRs = self.CTRs.generateCollectionFromXtal(self.xtal_unitcells)

        for calc, reference in zip(calc_CTRs, self.CTRs):
            expected = reference.sfI * self.reference_scale
            self.assertTrue(np.allclose(calc.sfI, expected, rtol=1e-02))


class TestCTRAccelerationBackendSelection(unittest.TestCase):
    def setUp(self):
        self.original_backend = CTRuc.CTR_ACCEL_BACKEND
        self.addCleanup(CTRuc.set_accel_backend, self.original_backend)
        self.addCleanup(self._reset_numba_probe)
        self._reset_numba_probe()

    @staticmethod
    def _reset_numba_probe():
        CTRuc._CTRcalc_accel = None
        CTRuc._NUMBA_ACCEL_IMPORT_ATTEMPTED = False
        CTRuc._NUMBA_ACCEL_IMPORT_ERROR = None
        CTRuc.HAS_NUMBA_ACCEL = False

    def testDefaultBackendIsValid(self):
        self.assertIn(CTRuc.CTR_ACCEL_BACKEND, {"cpp", "numba", "numpy"})

    def testNumbaBackendImportsLazilyWhenSelected(self):
        numba_backend = object()
        with mock.patch.object(
            CTRuc.importlib,
            "import_module",
            return_value=numba_backend,
        ) as import_module:
            CTRuc.set_accel_backend("numba")

        import_module.assert_called_once_with(
            "._CTRcalc_accel",
            package=CTRuc.__package__,
        )
        self.assertIs(CTRuc._CTRcalc_accel, numba_backend)
        self.assertTrue(CTRuc.HAS_NUMBA_ACCEL)

    def testSetAccelBackendCanSelectNumpy(self):
        CTRuc.set_accel_backend("numpy")
        self.assertEqual(CTRuc.CTR_ACCEL_BACKEND, "numpy")


class TestCTRcalculationNumPy(StructureFactorValidationMixin, unittest.TestCase):
    def setUp(self):
        original_backend = CTRuc.CTR_ACCEL_BACKEND
        self.addCleanup(CTRuc.set_accel_backend, original_backend)
        fp = os.path.split(__file__)[0]
        self.xtal_unitcells = CTRcalc.SXRDCrystal.fromFile(
            os.path.join(fp, "testdata/0V12_calculated.xpr")
        )
        self.CTRs = CTRplotutil.CTRCollection.fromANAROD(
            os.path.join(fp, "testdata/0V12_calculated.dat"),
            RODexport=True,
        )
        pt100 = CTRcalc.UnitCell([3.9242, 3.9242, 3.9242], [90.0000, 90.0000, 90.0000])
        self.xtal_unitcells.setGlobalReferenceUnitCell(
            pt100, util.z_rotation(np.deg2rad(45.0))
        )
        # This legacy reference file was generated when amplitudes were
        # normalized by unit-cell volume. Canonical F values in electrons are
        # therefore larger by the reference-cell volume.
        self.reference_scale = pt100.volume
        CTRuc.set_accel_backend("numpy")

    def testStructureFactorEqual(self):
        self.assert_structure_factors_match_volume_normalized_reference()


class TestCTRcalculationNumba(StructureFactorValidationMixin, unittest.TestCase):
    def setUp(self):
        original_backend = CTRuc.CTR_ACCEL_BACKEND
        self.addCleanup(CTRuc.set_accel_backend, original_backend)
        fp = os.path.split(__file__)[0]
        self.xtal_unitcells = CTRcalc.SXRDCrystal.fromFile(
            os.path.join(fp, "testdata/0V12_calculated.xpr")
        )
        self.CTRs = CTRplotutil.CTRCollection.fromANAROD(
            os.path.join(fp, "testdata/0V12_calculated.dat"),
            RODexport=True,
        )
        pt100 = CTRcalc.UnitCell([3.9242, 3.9242, 3.9242], [90.0000, 90.0000, 90.0000])
        self.xtal_unitcells.setGlobalReferenceUnitCell(
            pt100, util.z_rotation(np.deg2rad(45.0))
        )
        # This legacy reference file was generated when amplitudes were
        # normalized by unit-cell volume. Canonical F values in electrons are
        # therefore larger by the reference-cell volume.
        self.reference_scale = pt100.volume
        if not CTRuc.ctr_numba_accel_available():
            self.skipTest(
                "Cannot perform Numba tests: _CTRcalc_accel library was not "
                "imported. Is Numba installed?"
            )
        CTRuc.set_accel_backend("numba")

    def testStructureFactorEqual(self):
        self.assert_structure_factors_match_volume_normalized_reference()


class TestCTRcalculationCpp(StructureFactorValidationMixin, unittest.TestCase):
    def setUp(self):
        if not CTRuc.HAS_CPP_ACCEL:
            self.skipTest(
                "Cannot perform C++ tests: _CTRcalc_cpp library was not "
                "imported. Was the C++ extension built?"
            )
        original_backend = CTRuc.CTR_ACCEL_BACKEND
        self.addCleanup(CTRuc.set_accel_backend, original_backend)
        fp = os.path.split(__file__)[0]
        self.xtal_unitcells = CTRcalc.SXRDCrystal.fromFile(
            os.path.join(fp, "testdata/0V12_calculated.xpr")
        )
        self.CTRs = CTRplotutil.CTRCollection.fromANAROD(
            os.path.join(fp, "testdata/0V12_calculated.dat"),
            RODexport=True,
        )
        pt100 = CTRcalc.UnitCell([3.9242, 3.9242, 3.9242], [90.0000, 90.0000, 90.0000])
        self.xtal_unitcells.setGlobalReferenceUnitCell(
            pt100, util.z_rotation(np.deg2rad(45.0))
        )
        # This legacy reference file was generated when amplitudes were
        # normalized by unit-cell volume. Canonical F values in electrons are
        # therefore larger by the reference-cell volume.
        self.reference_scale = pt100.volume
        CTRuc.set_accel_backend("cpp")

    def testStructureFactorEqual(self):
        self.assert_structure_factors_match_volume_normalized_reference()

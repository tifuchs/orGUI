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

import dataclasses
import copy
import unittest
from unittest import mock
import os
import tempfile

import numpy as np

from ... import util
from .. import CTRcalc, CTRfilm, CTRplotutil, CTRsymmetry, CTRuc
from ..CTRdistributions import PoissonProfile, SkellamProfile, SurfaceProfile
from ..CTRutil import generate_surface_termination_cells


class _FakePy3DmolModel:
    def __init__(self):
        self.styles = []

    def setStyle(self, selection, style):
        self.styles.append((selection, style))


class _FakePy3DmolView:
    def __init__(self, displayed=False):
        self._orgui_plot3d_backend = "py3dmol"
        self.uniqueid = "displayed" if displayed else None
        self.models = []
        self.zoom_count = 0
        self.update_count = 0

    def addModel(self, data, file_format):
        self.models.append((data, file_format, _FakePy3DmolModel()))

    def getModel(self):
        return self.models[-1][2]

    def zoomTo(self):
        self.zoom_count += 1

    def update(self):
        self.update_count += 1


class TestPlot3d(unittest.TestCase):
    def setUp(self):
        self.cell = CTRuc.UnitCell(
            [2.0, 3.0, 4.0],
            [90.0, 90.0, 90.0],
            name="plot-test",
        )
        self.cell.addAtom("C", [0.5, 0.0, 0.0], 0.1, 0.1, 1.0)
        self.cell.addAtom("O", [0.0, 0.5, 0.0], 0.1, 0.1, 1.0)
        self.cell.f = np.empty((2, 13), dtype=np.float64)

    @staticmethod
    def _xyz_positions(model_data):
        lines = model_data.splitlines()
        return np.array(
            [[float(value) for value in line.split()[1:]] for line in lines[2:]]
        )

    def test_py3dmol_models_preserve_positions_styles_and_increment(self):
        view = _FakePy3DmolView(displayed=True)

        result = self.cell.plot3d(
            2,
            1,
            1,
            figure=view,
            backend="py3dmol",
        )

        self.assertIs(result, view)
        self.assertEqual(len(view.models), 2)
        self.assertEqual(view.update_count, 1)
        self.assertEqual(view.zoom_count, 1)
        carbon_xyz, file_format, carbon_model = view.models[0]
        self.assertEqual(file_format, "xyz")
        self.assertEqual(carbon_xyz.splitlines()[0], "2")
        expected = self.cell.pos_cart_all(2, 1, 1)
        expected_carbon = np.column_stack(
            (
                expected["x"][expected["name"] == "C"],
                expected["y"][expected["name"] == "C"],
                expected["z"][expected["name"] == "C"],
            )
        )
        np.testing.assert_allclose(
            self._xyz_positions(carbon_xyz),
            expected_carbon,
            atol=1e-14,
        )
        selection, style = carbon_model.styles[0]
        self.assertEqual(selection, {})
        self.assertEqual(style["sphere"]["radius"], 0.76)
        self.assertRegex(style["sphere"]["color"], r"^#[0-9a-f]{6}$")

        second_result = self.cell.plot3d(figure=view)
        self.assertIs(second_result, view)
        self.assertEqual(len(view.models), 4)
        self.assertEqual(view.update_count, 2)

    def test_radius_scale_applies_to_py3dmol_radius(self):
        view = _FakePy3DmolView()

        self.cell.plot3d(
            figure=view,
            backend="py3dmol",
            radius_scale=0.5,
        )

        radii = [
            model.styles[0][1]["sphere"]["radius"]
            for _, _, model in view.models
        ]
        self.assertEqual(radii, [0.38, 0.33])

    def test_py3dmol_omits_unoccupied_atoms(self):
        self.cell.basis[0, 6] = 0.0
        view = _FakePy3DmolView()

        self.cell.plot3d(
            figure=view,
            backend="py3dmol",
            occuon=True,
        )

        self.assertEqual(len(view.models), 1)
        self.assertTrue(view.models[0][0].splitlines()[2].startswith("O "))

    def test_mayavi_backend_keeps_points3d_contract(self):
        figure = mock.Mock()
        figure._orgui_plot3d_backend = "mayavi"
        mlab = mock.Mock()

        with mock.patch.object(CTRuc.importlib, "import_module", return_value=mlab):
            result = self.cell.plot3d(
                figure=figure,
                backend="mayavi",
                color=(0.1, 0.2, 0.3),
                resolution=17,
            )

        self.assertIs(result, figure)
        self.assertEqual(mlab.points3d.call_count, 2)
        diameters = [
            call.kwargs["scale_factor"] for call in mlab.points3d.call_args_list
        ]
        self.assertEqual(diameters, [1.52, 1.32])
        kwargs = mlab.points3d.call_args.kwargs
        self.assertEqual(kwargs["color"], (0.1, 0.2, 0.3))
        self.assertEqual(kwargs["resolution"], 17)
        self.assertIs(kwargs["figure"], figure)

    def test_backend_validation_and_missing_dependency(self):
        for radius_scale in (0, -1, np.nan, np.inf, None, "invalid"):
            with self.subTest(radius_scale=radius_scale):
                with self.assertRaisesRegex(ValueError, "positive and finite"):
                    self.cell.plot3d(radius_scale=radius_scale)

        with self.assertRaisesRegex(ValueError, "Unknown 3D plotting backend"):
            self.cell.plot3d(backend="unknown")

        mayavi_figure = mock.Mock()
        mayavi_figure._orgui_plot3d_backend = "mayavi"
        with self.assertRaisesRegex(TypeError, "uses mayavi"):
            self.cell.plot3d(figure=mayavi_figure, backend="py3dmol")

        with mock.patch.object(
            CTRuc.importlib,
            "import_module",
            side_effect=ImportError("simulated missing dependency"),
        ):
            with self.assertRaisesRegex(ImportError, r"orGUI\[full\]"):
                self.cell.plot3d(backend="py3dmol")

    def test_auto_prefers_py3dmol_in_kernel(self):
        view = _FakePy3DmolView()
        py3dmol = mock.Mock()
        py3dmol.view.return_value = view

        with (
            mock.patch.object(CTRuc, "_in_jupyter_kernel", return_value=True),
            mock.patch.object(
                CTRuc.importlib,
                "import_module",
                return_value=py3dmol,
            ) as import_module,
        ):
            result = self.cell.plot3d(backend="auto")

        self.assertIs(result, view)
        import_module.assert_called_once_with("py3Dmol")
        self.assertEqual(len(view.models), 2)

    def test_crystal_updates_one_shared_view_once(self):
        surface = copy.deepcopy(self.cell)
        crystal = CTRcalc.SXRDCrystal(self.cell, surface)
        view = _FakePy3DmolView(displayed=True)

        result = crystal.plot3d(
            1,
            1,
            1,
            figure=view,
            backend="py3dmol",
            radius_scale=0.5,
        )

        self.assertIs(result, view)
        self.assertEqual(len(view.models), 4)
        self.assertEqual(view.update_count, 1)
        self.assertEqual(view.models[0][2].styles[0][1]["sphere"]["radius"], 0.38)


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


class TestUnitCellCifImport(unittest.TestCase):
    def test_cif_import_orders_atoms_by_fractional_coordinate(self):
        cif = """data_ordering
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt P1
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
O1 O 0.3 0.4 0.8
O2 O 0.1 0.8 0.2
C1 C 0.4 0.2 0.2
Fe1 Fe 0.9 0.1 0.2
"""
        with tempfile.NamedTemporaryFile(suffix=".cif", mode="w") as handle:
            handle.write(cif)
            handle.flush()
            unitcell = CTRcalc.UnitCell.fromFile(handle.name)

        self.assertEqual(unitcell.names, ["Fe", "C", "O", "O"])
        np.testing.assert_allclose(
            unitcell.basis[:, 1:4],
            [[0.9, 0.1, 0.2], [0.4, 0.2, 0.2], [0.1, 0.8, 0.2], [0.3, 0.4, 0.8]],
        )
        self.assertTrue(unitcell._allow_layer_order_normalization)


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

    def bind_surface(self, surface, loc=0.0, height=4.0):
        film = CTRfilm.Film(copy.deepcopy(self.unitcell), name="underlying_film")
        film.basis[0] = 2.0
        film.createLayers()
        surface.stack_on(
            loc,
            height,
            film.end_layer_number,
            below_state=film.layer_state,
            below_component=film,
        )
        return film

    def test_serialization_round_trip(self):
        surface = CTRfilm.PoissonSurface(self.unitcell)
        surface.basis = np.array([4.0, 1.0, 0.25])
        surface.basis_0 = np.copy(surface.basis)
        surface.errors = np.array([0.2, 0.1, 0.05])

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

    def test_create_layers_assigns_convolved_occupancies(self):
        surface = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=2.0, alpha=0.5),
        )

        self.bind_surface(surface)

        lower, upper = surface.profile.support()
        offsets = np.arange(lower, upper + 1)
        material = surface.profile.occupancy(offsets)
        top_stop = np.flatnonzero(
            material > surface.profile.tail_probability
        )[-1] + 1
        offsets = offsets[:top_stop]
        material = material[:top_stop]
        expected_surface = 0.75 * surface.profile.surface_occupancy(offsets)
        expected_film = 0.75 * (material - (offsets < 0))
        expected_reference = -expected_surface
        represented = (expected_surface > surface.profile.tail_probability) | (
            np.abs(expected_film) > surface.profile.tail_probability
        )
        actual_surface = np.concatenate(
            [uc.coherentDomainOccupancy for uc in surface.layer_ucs]
        )
        actual_film = np.concatenate(
            [uc.coherentDomainOccupancy for uc in surface.film_layer_ucs]
        )
        actual_reference = np.concatenate(
            [
                uc.coherentDomainOccupancy
                for uc in surface._film_termination_ucs.values()
            ]
        )
        np.testing.assert_allclose(
            np.sort(actual_surface), np.sort(expected_surface[represented])
        )
        np.testing.assert_allclose(
            np.sort(actual_film), np.sort(expected_film[represented])
        )
        np.testing.assert_allclose(
            np.sort(actual_reference), np.sort(expected_reference[represented])
        )

    def test_profile_support_and_serialization_are_numerical_only(self):
        profile = PoissonProfile(
            mean_change=1.0,
            alpha=0.25,
            offset=0.25,
            tail_probability=1e-8,
        )
        surface = CTRfilm.PoissonSurface(self.unitcell, profile=profile)
        surface.basis_0[:] = surface.basis
        self.bind_surface(surface, loc=7.0, height=11.0)

        self.assertEqual(surface.pos_absolute, 11.0)
        self.assertEqual(surface.stacking_height_absolute, 13.5)
        self.assertGreater(surface.height_absolute, 13.5)
        self.assertLessEqual(
            surface.layer_ucs[-1].coherentDomainOccupancy[-1],
            self.unitcell.coherentDomainOccupancy[0],
        )

        restored = CTRfilm.PoissonSurface.fromStr(surface.toStr())
        self.assertAlmostEqual(restored.profile.mean_change, 1.0)
        self.assertAlmostEqual(restored.profile.alpha, 0.25)
        self.assertAlmostEqual(restored.profile.offset, 0.25)
        self.assertAlmostEqual(restored.profile.tail_probability, 1e-8)
        self.assertEqual(restored.toStr(), surface.toStr())

    def test_supercell_termination_bank_selects_complete_relaxed_slabs(self):
        slab = self.unitcell.supercell((1, 1, 2), symmetry="independent")
        termination_1 = slab.affine_layer_transform([0, 0, 0]).as_surface_termination(
            1,
            name="surface_termination_1",
        )
        termination_0 = slab.affine_layer_transform([0, 0, 1]).as_surface_termination(
            0,
            name="surface_termination_0",
        )
        top_atom = np.argmax(termination_0.basis[:, 3])
        termination_0.basis[top_atom, 3] += 0.025
        termination_0.basis_0[top_atom, 3] += 0.025

        surface = CTRfilm.PoissonSurface(
            {0: termination_0, 1: termination_1},
            profile=PoissonProfile(1.0, alpha=1.0),
            name="surface",
        )
        serialized = surface.toStr()
        restored = CTRfilm.PoissonSurface.fromStr(serialized)
        self.assertEqual(set(restored.termination_cells), {0.0, 1.0})
        self.assertIn("TerminationUnitCell 0", restored.toStr())

        self.bind_surface(surface)
        self.assertEqual(
            {cell.basis.shape[0] for cell in surface.termination_cells.values()},
            {4},
        )
        for layer, cell in surface.termination_cells.items():
            np.testing.assert_array_equal(cell.basis[:, 7], layer)
            self.assertEqual(cell.layer_behavior, "select")
        self.assertTrue(
            all(
                cell.coherentDomainMatrix
                for cell in surface.layer_ucs
            )
        )
        self.assertEqual(surface.toStr(), serialized)

    def test_surface_termination_helper_accepts_film_cycle_owner(self):
        film = CTRfilm.Film(copy.deepcopy(self.unitcell))
        slab = self.unitcell.supercell((1, 1, 2), symmetry="independent")

        terminations = generate_surface_termination_cells(
            slab,
            film,
            name_template="relaxed_{layer:g}",
        )

        self.assertEqual(set(terminations), {0.0, 1.0})
        self.assertEqual(
            [terminations[layer].name for layer in (0.0, 1.0)],
            ["relaxed_0", "relaxed_1"],
        )
        for layer, cell in terminations.items():
            self.assertEqual(cell.basis.shape[0], 4)
            np.testing.assert_array_equal(cell.basis[:, 7], layer)
            self.assertEqual(cell.layer_behavior, "select")
            self.assertEqual(tuple(cell.layer_cycle.layers), (layer,))

        with self.assertRaisesRegex(ValueError, "integer multiple"):
            generate_surface_termination_cells(slab, (0, 1, 2))

    def test_termination_bank_requires_one_cell_per_film_cycle(self):
        termination = self.unitcell.as_surface_termination(0)
        surface = CTRfilm.PoissonSurface(
            {0: termination},
            profile=PoissonProfile(0.0),
        )
        film = CTRfilm.Film(copy.deepcopy(self.unitcell))
        film.basis[0] = 2.0
        crystal = CTRcalc.SXRDCrystal(
            self.unitcell,
            film,
            surface,
            stacking=np.array([1, 2]),
        )
        with self.assertRaisesRegex(ValueError, "exactly one surface unit cell"):
            crystal.apply_stacking()

    def test_distribution_support_bounds_omitted_tail(self):
        poisson_profile = PoissonProfile(2.5, tail_probability=1e-9)
        _, upper = poisson_profile.support()
        self.assertLessEqual(poisson_profile.occupancy([upper])[0], 1e-9)

        skellam_profile = SkellamProfile(0.35, asymmetry=0.0, tail_probability=1e-8)
        lower, upper = skellam_profile.support(layer_count=2)
        self.assertLess(lower, 0)
        self.assertGreater(upper, 0)

    def test_step_poisson_convolution_preserves_signed_mean(self):
        for width in (2.3, -2.3):
            for alpha in (0.0, 0.25, 1.0):
                offset = 0.4
                profile = PoissonProfile(
                    mean_change=width,
                    alpha=alpha,
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

    def test_surface_occupancy_is_neighbor_difference(self):
        for profile in (
            PoissonProfile(2.5, alpha=0.4, tail_probability=1e-12),
            PoissonProfile(-2.5, alpha=0.4, tail_probability=1e-12),
        ):
            lower, upper = profile.support()
            offsets = np.arange(lower, upper + 1)
            material = profile.occupancy(offsets)
            exposed = profile.surface_occupancy(offsets)

            np.testing.assert_allclose(exposed[:-1], material[:-1] - material[1:])
            self.assertEqual(exposed[-1], material[-1])
            self.assertTrue(np.all(exposed >= 0.0))
            self.assertAlmostEqual(np.sum(exposed), material[0])

    def test_surface_profile_rejects_nonmonotonic_occupancy(self):
        class InvalidProfile(SurfaceProfile):
            def support(self):
                return 0, 2

            def occupancy(self, offsets):
                return np.array([0.5, 0.6, 0.1])

        with self.assertRaisesRegex(ValueError, "must not increase"):
            InvalidProfile().surface_occupancy(np.arange(3))

    def test_probability_is_step_poisson_convolution(self):
        profile = PoissonProfile(mean_change=2.0, alpha=0.25)
        changes = np.arange(0, 8)
        poisson_probability = np.exp(-0.5) * 0.5 ** np.arange(7) / np.array(
            [1, 1, 2, 6, 24, 120, 720]
        )
        expected = np.zeros_like(changes, dtype=np.float64)
        expected[1:] += 0.5 * poisson_probability
        expected[2:] += 0.5 * poisson_probability[:-1]
        np.testing.assert_allclose(profile.probability(changes), expected)

        dissolution = PoissonProfile(mean_change=-2.0, alpha=0.25)
        np.testing.assert_allclose(
            dissolution.probability(-changes),
            expected,
        )

    def test_alpha_selects_layer_by_layer_and_poisson_limits(self):
        layer_by_layer = PoissonProfile(mean_change=-1.25, alpha=0.0)
        np.testing.assert_allclose(
            layer_by_layer.probability([-1, -2]),
            [0.75, 0.25],
        )

        poisson_only = PoissonProfile(mean_change=1.25, alpha=1.0)
        np.testing.assert_allclose(
            poisson_only.probability([0, 1, 2]),
            np.exp(-1.25) * np.array([1.0, 1.25, 1.25**2 / 2.0]),
        )

        with self.assertRaisesRegex(ValueError, "alpha"):
            PoissonProfile(
                mean_change=1.0,
                alpha=1.01,
            )

    def test_poisson_surface_growth_and_etching_are_signed(self):
        growth = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=1.0, alpha=0.4),
        )
        self.bind_surface(growth, height=8.0)
        growth_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in growth.layer_ucs]
        )
        self.assertTrue(np.all(growth_occupancy > 0))
        self.assertAlmostEqual(growth.mean_height_absolute, 10.0)

        etching = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=-1.0, alpha=0.4),
        )
        self.bind_surface(etching, height=8.0)
        etching_surface_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in etching.layer_ucs]
        )
        etching_film_occupancy = np.concatenate(
            [layer.coherentDomainOccupancy for layer in etching.film_layer_ucs]
        )
        self.assertTrue(np.all(etching_surface_occupancy > 0))
        self.assertTrue(np.all(etching_film_occupancy < 0))
        self.assertAlmostEqual(etching.mean_height_absolute, 6.0)

    def test_offset_is_an_independent_fit_parameter(self):
        surface = CTRfilm.PoissonSurface(
            self.unitcell,
            profile=PoissonProfile(mean_change=1.0, alpha=0.4, offset=0.25),
        )
        surface.basis_0[:] = surface.basis
        surface.addFitParameter("offset", limits=(-2.0, 2.0))

        np.testing.assert_allclose(surface.getStartParamAndLimits()[0], [0.25])
        surface.setFitParameters([-0.5])

        self.assertAlmostEqual(surface.profile.mean_change, 1.0)
        self.assertAlmostEqual(surface.profile.alpha, 0.4)
        self.assertAlmostEqual(surface.profile.offset, -0.5)
        self.assertAlmostEqual(surface.profile.expected_height_change, 0.5)

    def test_zero_width_identical_surface_is_zero_net_correction(self):
        film = CTRfilm.Film(copy.deepcopy(self.unitcell), name="film")
        film.basis[0] = 2.0
        surface = CTRfilm.PoissonSurface(
            copy.deepcopy(self.unitcell),
            profile=PoissonProfile(mean_change=0.0),
        )
        crystal = CTRcalc.SXRDCrystal(
            self.unitcell, film, surface, stacking=np.array([1, 2])
        )

        crystal.apply_stacking()

        self.assertEqual(surface.stacking_height_absolute, film.height_absolute)
        self.assertTrue(
            any(layer.coherentDomainMatrix for layer in surface.layer_ucs)
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
        np.testing.assert_allclose(surface.optical_profile()[:, 1:], 0.0)

    def test_zero_width_replaces_top_film_layer_with_surface_structure(self):
        film_cell = copy.deepcopy(self.unitcell)
        surface_cell = copy.deepcopy(self.unitcell)
        surface_cell.names[:] = ["O", "O"]
        surface_cell.basis[:, 0] = 8
        surface_cell.setEnergy(10000.0)
        film = CTRfilm.Film(film_cell, name="film")
        film.basis[0] = 2.0
        surface = CTRfilm.PoissonSurface(
            surface_cell,
            profile=PoissonProfile(mean_change=0.0),
            name="surface",
        )
        crystal = CTRcalc.SXRDCrystal(
            self.unitcell, film, surface, stacking=np.array([1, 2])
        )
        crystal.apply_stacking()

        surface_occupancies = np.concatenate(
            [layer.coherentDomainOccupancy for layer in surface.layer_ucs]
        )
        film_corrections = np.concatenate(
            [layer.coherentDomainOccupancy for layer in surface.film_layer_ucs]
        )
        reference_corrections = np.concatenate(
            [
                cell.coherentDomainOccupancy
                for cell in surface._film_termination_ucs.values()
            ]
        )
        np.testing.assert_allclose(surface_occupancies, [0.75])
        np.testing.assert_allclose(film_corrections, [0.0])
        np.testing.assert_allclose(reference_corrections, [-0.75])
        h = np.zeros(3)
        ell = np.linspace(0.1, 0.3, 3)
        self.assertTrue(np.any(np.abs(surface.F_uc(h, h, ell)) > 0.0))

    def test_surface_requires_compatible_underlying_film(self):
        surface = CTRfilm.PoissonSurface(
            copy.deepcopy(self.unitcell), profile=PoissonProfile(1.0)
        )
        crystal = CTRcalc.SXRDCrystal(
            self.unitcell, surface, stacking=np.array([1])
        )
        with self.assertRaisesRegex(ValueError, "immediately above a Film"):
            crystal.apply_stacking()

        mismatched = copy.deepcopy(self.unitcell)
        mismatched.a[2] *= 1.1
        film = CTRfilm.Film(mismatched)
        film.basis[0] = 2.0
        crystal = CTRcalc.SXRDCrystal(
            self.unitcell, film, surface, stacking=np.array([1, 2])
        )
        with self.assertRaisesRegex(ValueError, "integer multiple"):
            crystal.apply_stacking()

    def test_repeated_stacking_refreshes_underlying_film_views(self):
        film = CTRfilm.Film(copy.deepcopy(self.unitcell), name="film")
        film.basis[0] = 2.0
        surface = CTRfilm.PoissonSurface(
            copy.deepcopy(self.unitcell),
            profile=PoissonProfile(1.0, alpha=0.5),
            name="surface",
        )
        crystal = CTRcalc.SXRDCrystal(
            self.unitcell, film, surface, stacking=np.array([1, 2])
        )
        crystal.apply_stacking()
        first_counts = [
            len(layer.coherentDomainMatrix) for layer in surface.film_layer_ucs
        ]

        film.unitcell.basis[:, 6] *= 0.5
        crystal.apply_stacking()

        self.assertEqual(
            [len(layer.coherentDomainMatrix) for layer in surface.film_layer_ucs],
            first_counts,
        )
        self.assertTrue(
            all(
                np.shares_memory(layer.basis, film.unitcell.basis)
                for layer in surface.film_layer_ucs
            )
        )
        np.testing.assert_allclose(
            np.concatenate([layer.basis[:, 6] for layer in surface.film_layer_ucs]),
            0.5,
        )

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

    def test_split_in_layers_rejects_interleaved_native_layers(self):
        unitcell = CTRcalc.UnitCell([3.0, 3.0, 6.0], [90.0, 90.0, 90.0])
        for name, z, layer in (("Fe", 0.1, 1), ("Re", 0.6, 2), ("O", 0.2, 1)):
            unitcell.addAtom(name, [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=layer)

        with self.assertRaisesRegex(ValueError, "must be written in layer order"):
            unitcell.split_in_layers()

    def test_split_in_layers_normalizes_cif_metadata_and_couplings(self):
        unitcell = CTRcalc.UnitCell([3.0, 3.0, 6.0], [90.0, 90.0, 90.0])
        for name, z, layer in (("Fe", 0.1, 1), ("Re", 0.6, 2), ("O", 0.2, 1), ("O", 0.7, 2)):
            unitcell.addAtom(name, [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=layer)
        unitcell._allow_layer_order_normalization = True
        unitcell.symmetry_metadata = CTRsymmetry.SurfaceSymmetryModel(
            CTRsymmetry.SurfaceCellSpec(
                (3.0, 3.0, 6.0),
                (90.0, 90.0, 90.0),
                np.identity(3),
                layer_origins=(0.0, 0.5),
            ),
            (),
            atoms=[
                CTRsymmetry.GeneratedWyckoffAtom(
                    atom_index=2,
                    element="O",
                    site_id="O_1",
                    wyckoff_label="1a",
                    parent_fractional=np.array([0.0, 0.0, 0.2]),
                    surface_fractional=np.array([0.0, 0.0, 0.2]),
                    layer=1,
                    couplings=(
                        CTRsymmetry.WyckoffCoupling(
                            atom_index=2,
                            coordinate="z",
                            variable="u",
                            constant=0.2,
                            factor=1.0,
                            site_id="O_1",
                        ),
                    ),
                    site_couplings=(
                        CTRsymmetry.WyckoffSiteCoupling(
                            atom_index=2,
                            coordinate="z",
                            axis="z",
                            factor=1.0,
                            site_id="O_1",
                        ),
                    ),
                ),
            ],
        )

        layers = unitcell.split_in_layers()

        self.assertEqual(unitcell.names, ["Fe", "O", "Re", "O"])
        self.assertEqual(layers[1.0].names, ["Fe", "O"])
        self.assertEqual(layers[2.0].names, ["Re", "O"])
        self.assertTrue(np.shares_memory(layers[1.0].basis, unitcell.basis))
        self.assertTrue(np.shares_memory(layers[2.0].basis, unitcell.basis))
        atom = unitcell.atom_wyckoff_metadata(1)
        self.assertEqual(atom.couplings[0].atom_index, 1)
        self.assertEqual(atom.site_couplings[0].atom_index, 1)

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

    def test_unitcell_translate_layered_shifts_absolute_parameter_values(self):
        unitcell = self.make_layered_unitcell()
        unitcell.addFitParameter((0, "z"), limits=(-0.1, 2.0), name="z_abs")
        unitcell.setFitParameters([0.05])

        transformed = unitcell.translate_layered([0, 0, 3])

        np.testing.assert_allclose(transformed.getInitialParameters(), [1.05])
        transformed.setFitParameters([1.05])
        self.assertAlmostEqual(transformed.basis[0, 3], 1.05)

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

    def test_unitcell_translate_layered_updates_wyckoff_representatives(self):
        unitcell = self.make_single_wyckoff_cell()

        transformed = unitcell.translate_layered([1, -1, 0])

        site = transformed.wyckoff_sites()[0]
        np.testing.assert_allclose(
            site["representative_parent_fractional"],
            [1.25, -1.0, 0.0],
        )
        parameter = transformed.addWyckoffShift(
            "C_1",
            "x",
            absolute_limits=(1.0, 1.5),
        )
        self.assertEqual(parameter.limits, (-0.25, 0.25))

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

    def test_unitcell_supercell_rejects_fractional_repeat_counts(self):
        unitcell = self.make_single_wyckoff_cell()

        with self.assertRaisesRegex(ValueError, "positive integers"):
            unitcell.supercell((1.5, 1, 1))

    def test_unitcell_supercell_rescales_coherent_domain_translations(self):
        unitcell = CTRcalc.UnitCell(
            [2.0, 3.0, 4.0],
            [90.0, 90.0, 90.0],
            name="domain",
        )
        unitcell.addAtom("C", [0.25, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
        unitcell.layerpos[1.0] = 0.0
        domain = np.eye(3, 4)
        domain[0, 3] = 0.5
        unitcell.coherentDomainMatrix = [domain]

        supercell = unitcell.supercell((2, 1, 1))

        np.testing.assert_allclose(
            supercell.coherentDomainMatrix[0][:, -1],
            [0.25, 0.0, 0.0],
        )
        np.testing.assert_allclose(supercell.pos_cart(0), unitcell.pos_cart(0))

    def test_unitcell_supercell_uses_current_coordinates_as_new_reference(self):
        unitcell = self.make_single_wyckoff_cell()
        unitcell.addRelParameter(([0], "x"), [1.0], (-1.0, 1.0), name="dx")
        unitcell.setFitParameters([0.05])

        supercell = unitcell.supercell((2, 1, 1))

        np.testing.assert_allclose(supercell.basis[:, 1], [0.15, 0.65])
        np.testing.assert_allclose(supercell.basis_0, supercell.basis)
        supercell.addRelParameter(([0], "x"), [1.0], (-1.0, 1.0), name="dx2")
        supercell.setFitParameters([0.0])
        np.testing.assert_allclose(supercell.basis[:, 1], [0.15, 0.65])

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

        supercell.addWyckoffParameter("C_1", "u", limits=(0.1, 0.4))
        np.testing.assert_allclose(supercell.getInitialParameters(), [0.25])
        supercell.setFitParameters([0.35])

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
            limits=(0.1, 0.4),
        )
        supercell.setFitParameters([0.35])

        np.testing.assert_allclose(supercell.basis[:, 1], [0.175, 0.625])

        supercell = unitcell.supercell((2, 1, 1), symmetry="independent")
        supercell.addWyckoffShift("C_1_copy0", "x", limits=(-0.2, 0.2))
        supercell.setFitParameters([0.1])

        np.testing.assert_allclose(supercell.basis[:, 1], [0.175, 0.625])

    def test_film_forwards_wyckoff_parameters_to_template_unitcell(self):
        film = CTRfilm.Film(self.make_single_wyckoff_cell())

        film.addWyckoffParameter("C_1", "u", limits=(0.1, 0.4))
        np.testing.assert_allclose(film.getInitialParameters(), [0.25])
        film.setFitParameters([0.35])

        self.assertEqual(film.unitcell.wyckoff_sites()[0]["status"], "symmetry_preserving")
        np.testing.assert_allclose(film.unitcell.basis[:, 1], [0.35])

    def test_film_forwards_absolute_plural_wyckoff_parameters(self):
        film = CTRfilm.Film(self.make_single_wyckoff_cell())

        parameters = film.addWyckoffParameters(
            "C_1",
            limits={"u": (0.1, 0.4)},
        )

        self.assertEqual(len(parameters), 1)
        np.testing.assert_allclose(
            film.getStartParamAndLimits(),
            ([0.25], [0.1], [0.4]),
        )
        film.setFitParameters([0.3])
        np.testing.assert_allclose(film.unitcell.basis[:, 1], [0.3])

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
            limits=(0.1, 0.4),
            unitcell=["top", "bottom"],
        )
        np.testing.assert_allclose(interface.getInitialParameters(), [0.25, 0.25])
        interface.setFitParameters([0.35, 0.35])

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
            limits=(0.1, 0.4),
        )
        np.testing.assert_allclose(
            crystal.getStartParamAndLimits(),
            ([0.25], [0.1], [0.4]),
        )
        np.testing.assert_allclose(
            film1.getStartParamAndLimits()[1:],
            ([0.1], [0.4]),
        )
        crystal.setLimits([(0.2, 0.3)])
        np.testing.assert_allclose(
            film1.getStartParamAndLimits()[1:],
            ([0.2], [0.3]),
        )
        crystal.setParameters([0.35])

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
            limits=(0.1, 0.4),
        )
        np.testing.assert_allclose(crystal.getInitialParameters(), [0.25])
        crystal.setParameters([0.35])

        self.assertEqual(crystal.fitparnames, ["shared_interface_u"])
        self.assertEqual(interface.uc_top.wyckoff_sites()[0]["status"], "symmetry_preserving")
        self.assertEqual(interface.uc_bottom.wyckoff_sites()[0]["status"], "symmetry_preserving")
        np.testing.assert_allclose(interface.uc_top.basis[:, 1], [0.35])
        np.testing.assert_allclose(interface.uc_bottom.basis[:, 1], [0.35])

    def test_sxrdcrystal_distributes_absolute_wyckoff_parameter(self):
        bulk = self.make_single_wyckoff_cell()
        film1 = self.make_single_wyckoff_cell()
        film2 = self.make_single_wyckoff_cell()
        bulk.name = "bulk_template"
        film1.name = "film1"
        film2.name = "film2"
        film2.basis[:, 1] = 0.27
        film2.basis_0[:, 1] = 0.27
        film2.symmetry_metadata.sites = (
            dataclasses.replace(
                film2.symmetry_metadata.sites[0],
                variables={"u": 0.27},
            ),
        )
        crystal = CTRcalc.SXRDCrystal(bulk, film1, film2)

        crystal.addWyckoffParameter(
            {"film1": ("C_1", "u"), "film2": ("C_1", "u")},
            name="mean_u",
            limits=(0.2, 0.35),
        )

        np.testing.assert_allclose(crystal.getInitialParameters(), [0.26])
        crystal.setParameters([0.28])
        np.testing.assert_allclose(film1.basis[:, 1], [0.28])
        np.testing.assert_allclose(film2.basis[:, 1], [0.28])

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
        surface.basis[:] = [2, 0, 0]
        surface._bind_underlying_component(film)
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
            -1 / 3,
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
        self.assertAlmostEqual(interface.loc_absolute, 12.0)
        self.assertLess(interface.pos_absolute, interface.loc_absolute)
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
        surface.basis[:] = [5, 0, 0]

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

    def test_profile_interface_support_is_not_owned_by_film(self):
        lower = self.make_layered_unitcell("lower")
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_unitcell("top"),
            self.make_layered_unitcell("bottom"),
            profile=SkellamProfile(0.5),
        )
        film = CTRfilm.Film(self.make_layered_unitcell("film"))
        film.basis[0] = 18
        crystal = CTRcalc.SXRDCrystal(
            lower,
            interface,
            film,
            stacking=np.array([1, 2]),
        )

        crystal.apply_stacking()

        self.assertAlmostEqual(
            interface.stacking_height_absolute, interface.height_absolute
        )
        self.assertAlmostEqual(
            interface.stacking_loc_absolute, interface.loc_absolute
        )
        self.assertAlmostEqual(film.pos_absolute, interface.height_absolute)
        self.assertGreater(film.pos_absolute, interface.loc_absolute)
        self.assertAlmostEqual(
            film.height_absolute - interface.loc_absolute,
            film.basis[0] * film.unitcell.a[2] / len(film.layers),
        )
        layer_cycle = interface.layer_cycle.layers
        next_index = (layer_cycle.index(interface.end_layer_number) + 1) % len(
            layer_cycle
        )
        self.assertEqual(film.start_layer_number, layer_cycle[next_index])

        restored = CTRfilm.EpitaxyInterface.fromStr(interface.toStr())
        self.assertIsInstance(restored.profile, SkellamProfile)
        np.testing.assert_allclose(restored.basis, interface.basis)
        self.assertNotIn("support_cursor", interface.toStr())

    def test_profile_interface_can_consume_entire_film_width(self):
        lower = self.make_layered_unitcell("lower")
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_unitcell("top"),
            self.make_layered_unitcell("bottom"),
            profile=SkellamProfile(0.5),
        )
        interface.stack_on(
            0.0,
            0.0,
            lower.end_layer_number,
            below_state=lower.layer_state,
        )
        film = CTRfilm.Film(self.make_layered_unitcell("film"))
        film.basis[0] = (
            (interface.height_absolute - interface.loc_absolute)
            * len(film.layers) / film.unitcell.a[2]
        )
        crystal = CTRcalc.SXRDCrystal(
            lower, interface, film, stacking=np.array([1, 2])
        )

        crystal.apply_stacking()

        self.assertAlmostEqual(film.pos_absolute, interface.height_absolute)
        self.assertAlmostEqual(film.height_absolute, interface.height_absolute)
        self.assertEqual(
            sum(len(layer.coherentDomainMatrix) for layer in film.layer_ucs), 0
        )
        np.testing.assert_allclose(
            film.F_uc(
                np.array([0.0]), np.array([0.0]), np.array([1.0])
            ),
            0.0,
        )

    def test_profile_interface_rejects_film_narrower_than_support(self):
        lower = self.make_layered_unitcell("lower")
        interface = CTRfilm.EpitaxyInterface(
            self.make_layered_unitcell("top"),
            self.make_layered_unitcell("bottom"),
            profile=SkellamProfile(0.5),
        )
        film = CTRfilm.Film(self.make_layered_unitcell("film"))
        film.basis[0] = 3
        crystal = CTRcalc.SXRDCrystal(
            lower, interface, film, stacking=np.array([1, 2])
        )

        with self.assertRaisesRegex(ValueError, "shorter than lower-component"):
            crystal.apply_stacking()

    def test_profile_interface_owns_upper_support_material(self):
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

        self.assertTrue(np.all(top_occupancy >= 0.0))
        self.assertTrue(np.any(top_occupancy > 0.9))
        self.assertTrue(np.any(bottom_occupancy < 0.0))
        self.assertTrue(np.any(bottom_occupancy > 0.0))
        total_occupancy = top_occupancy + bottom_occupancy
        self.assertTrue(np.any(np.isclose(total_occupancy, 0.0)))
        self.assertTrue(np.any(np.isclose(total_occupancy, 1.0)))


class TestLegacyLayeredCTR(unittest.TestCase):
    def test_legacy_xtal_uses_corrected_interface_support(self):
        repository_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "..", "..")
        )
        fixture_root = os.path.join(
            repository_root, "examples", "CTR", "test"
        )
        xtal_path = os.path.join(
            fixture_root, "0001_fit_2V036_reference.xtal"
        )
        for path in (xtal_path,):
            if not os.path.exists(path):
                self.skipTest(f"Missing optional CTR fixture: {path}")
        xtal = CTRcalc.SXRDCrystal.fromFile(
            xtal_path
        )

        self.assertIsInstance(xtal["RuO2"], CTRfilm.Film)
        self.assertIsInstance(xtal["TiO2toRuO2"], CTRfilm.EpitaxyInterface)
        np.testing.assert_allclose(xtal["TiO2toRuO2"].basis, [0.35, 0.0])
        interface = xtal["TiO2toRuO2"]
        xtal.apply_stacking()
        self.assertAlmostEqual(
            interface.stacking_height_absolute, interface.height_absolute
        )
        self.assertAlmostEqual(
            interface.stacking_loc_absolute, interface.loc_absolute
        )
        self.assertIsInstance(interface.profile, SkellamProfile)
        self.assertNotIn("support_cursor", interface.toStr())


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

    def test_bulk_repeat_follows_reference_unit_cell(self):
        reference = CTRcalc.UnitCell(
            [3.0, 3.0, 4.0],
            [90.0, 90.0, 90.0],
            name="reference",
        )
        bulk = CTRcalc.UnitCell(
            [3.0, 3.0, 8.0],
            [90.0, 90.0, 90.0],
            name="bulk",
        )
        bulk.addAtom("C", [0.0, 0.0, 0.25], 0.1, 0.1, 1.0)
        bulk.setEnergy(10000.0)
        bulk.setReferenceUnitCell(reference)
        h = np.zeros(2)
        k = np.zeros(2)
        l = np.array([0.17, 0.41])
        atten = 0.07

        hkl_bulk = bulk.refHKLTransform @ np.vstack((h, k, l))
        repeat_atten = 2.0 * atten
        expected = bulk.F_uc_bulk_direct(*hkl_bulk, repeat_atten)
        expected /= 1 - np.exp(
            -2j * np.pi * hkl_bulk[2] - repeat_atten
        )

        with mock.patch.object(CTRuc, "CTR_ACCEL_BACKEND", "numpy"):
            actual = bulk.F_bulk(h, k, l, atten)

        np.testing.assert_allclose(actual, expected)
        np.testing.assert_allclose(hkl_bulk[2], 2.0 * l)

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
    def test_unitcell_atom_table_uses_fixed_width_columns(self):
        cell = CTRcalc.UnitCell([1.0, 1.0, 1.0], [90.0, 90.0, 90.0])
        cell.addAtom("O", [0.5, 0.0, 0.80569], 0.4119, 0.4119, 1.0, layer=2)

        lines = cell.toStr(showErrors=False).splitlines()

        self.assertIn(
            "Name        x/frac      y/frac      z/frac       iDW       oDW"
            "     occup  layerIdx",
            lines,
        )
        self.assertIn(
            "00  O          0.50000     0.00000     0.80569    0.4119"
            "    0.4119    1.0000         2",
            lines,
        )

        cell.errors = np.full_like(cell.basis, np.nan)
        cell.errors[0, 3] = 0.00007
        error_lines = cell.toStr(showErrors=True).splitlines()
        header = next(line for line in error_lines if line.startswith("Name"))
        atom = next(line for line in error_lines if line.startswith("00"))
        for column in ("x/frac", "y/frac", "z/frac", "iDW", "oDW", "occup"):
            self.assertGreater(header.index(column), 0)
        self.assertIn("(0.80569 +- 0.00007)", atom)
        self.assertEqual(len(header), len(atom) - 4)

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
            profile=PoissonProfile(mean_change=1.5, alpha=0.5),
            name="surface",
        )
        surface.basis_0[:] = surface.basis

        film.errors = np.array([0.2])
        interface.errors = np.array([0.03, 0.04])
        surface.errors = np.array([0.1, 0.15, 0.2])
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
        np.testing.assert_allclose(
            restored_xpr["surface"].errors,
            [0.1, 0.15, 0.2],
        )
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
            profile=PoissonProfile(mean_change=2.0, alpha=0.4),
            name="surface",
        )
        surface.basis_0[:] = surface.basis
        surface.addRelParameter(
            ("W", "alpha"),
            (1.0, 0.1),
            limits=(-0.4, 0.4),
        )
        surface.setFitParameters([0.25])
        surface.setFitErrors([0.1])
        np.testing.assert_allclose(surface.basis, [2.25, 0.425, 0.0])
        np.testing.assert_allclose(surface.errors[:2], [0.1, 0.01])
        self.assertTrue(np.isnan(surface.errors[2]))
        self.assertAlmostEqual(surface.profile.mean_change, 2.25)
        self.assertAlmostEqual(surface.profile.alpha, 0.425)
        self.assertAlmostEqual(surface.profile.expected_height_change, 2.25)

    def test_linear_fit_parameter_dictionary_round_trip(self):
        original = CTRfilm.PoissonSurface(
            self.make_layered_cell("surface_cell"),
            profile=PoissonProfile(mean_change=2.0, alpha=0.4),
            name="surface",
        )
        original.basis_0[:] = original.basis
        original.addRelParameter(
            ("W", "alpha"),
            (1.0, 0.1),
            limits=(-0.4, 0.4),
            name="coupled_roughness",
        )
        original.setFitParameters([0.4])
        original.setFitErrors([0.05])

        restored = CTRfilm.PoissonSurface(
            self.make_layered_cell("surface_cell"),
            profile=PoissonProfile(mean_change=2.0, alpha=0.4),
            name="surface",
        )
        restored.parametersFromDict(original.parametersToDict())

        np.testing.assert_allclose(restored.basis, original.basis)
        np.testing.assert_allclose(restored.errors, original.errors)
        np.testing.assert_allclose(restored.getFitErrors(), original.getFitErrors())
        self.assertEqual(restored.fitparnames, original.fitparnames)
        self.assertAlmostEqual(restored.profile.mean_change, 2.4)
        self.assertAlmostEqual(restored.profile.alpha, 0.44)

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
        np.testing.assert_allclose(crystal["surface"].basis, [2.0, 0.5, 0.0])
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
        np.testing.assert_allclose(xtal["RuO2surface"].basis, [-6.0, 1.0, 0.0])
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

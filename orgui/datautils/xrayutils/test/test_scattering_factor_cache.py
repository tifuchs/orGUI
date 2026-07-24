import unittest
from unittest import mock

import numpy as np
import xraydb

from .. import CTRuc, CTRutil


class TestScatteringFactorCaching(unittest.TestCase):
    def tearDown(self):
        CTRutil._readWaasmaier_cached.cache_clear()
        CTRutil._readDispersion_cached.cache_clear()

    def test_cached_values_match_xraydb_and_cannot_be_mutated(self):
        for species in ("Fe", "Re", "O", "Ti", "O2-", "Fe2+"):
            coefficients = CTRutil.readWaasmaier(species)
            baseline = CTRutil._readWaasmaier_cached(
                CTRutil._normalize_species(species)
            )
            np.testing.assert_allclose(coefficients, baseline)
            coefficients[0] += 1.0
            np.testing.assert_allclose(
                CTRutil.readWaasmaier(species), baseline
            )
            for energy in (10000.0, 20000.0):
                expected = (
                    xraydb.f1_chantler(CTRutil.atomic_number(species), energy),
                    xraydb.f2_chantler(CTRutil.atomic_number(species), energy),
                )
                np.testing.assert_allclose(
                    CTRutil.readDispersion(species, energy), expected
                )

    def test_repeated_energy_and_species_groups_do_not_repeat_lookups(self):
        cell = CTRuc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
        for index, species in enumerate(("Fe", "Re", "O") * 32):
            cell.addAtom(species, [index / 96.0, 0.0, 0.0], 0.0, 0.0, 1.0)

        with (
            mock.patch.object(
                CTRuc, "readWaasmaier", wraps=CTRuc.readWaasmaier
            ) as waasmaier,
            mock.patch.object(
                CTRuc, "readDispersion", wraps=CTRuc.readDispersion
            ) as dispersion,
        ):
            cell.setEnergy(10000.0)
            first = np.array(cell.f, copy=True)
            cell.setEnergy(10000.0)
            self.assertEqual(waasmaier.call_count, 3)
            self.assertEqual(dispersion.call_count, 3)
            cell.setEnergy(20000.0)
            self.assertEqual(waasmaier.call_count, 6)
            self.assertEqual(dispersion.call_count, 6)
        self.assertFalse(np.allclose(first[:, 11:], cell.f[:, 11:]))

    def test_layer_splitting_defers_or_copies_form_factors(self):
        cell = CTRuc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
        cell.addAtom("Fe", [0.0, 0.0, 0.0], 0.0, 0.0, 1.0, layer=0)
        cell.addAtom("O", [0.0, 0.0, 0.5], 0.0, 0.0, 1.0, layer=1)
        layers = cell.split_in_layers()
        self.assertFalse(hasattr(cell, "f"))
        self.assertTrue(all(not hasattr(layer, "f") for layer in layers.values()))

        cell.setEnergy(20000.0)
        layers = cell.split_in_layers()
        first_layer = next(iter(layers.values()))
        before = cell.f[0, 0]
        first_layer.f[0, 0] += 1.0
        self.assertEqual(cell.f[0, 0], before)
        self.assertEqual(first_layer._E, 20000.0)


if __name__ == "__main__":
    unittest.main()

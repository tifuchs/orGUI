"""Regression tests for the Poisson surface flat-height decomposition.

The reference for every state correction is an independently stacked ``Film``
from ``_poisson_oracle``, which never touches ``PoissonSurface``. These cover
increment 2 of ``doc/design/incoherent_ctr_models.md``.
"""

import copy
import unittest

import numpy as np

from .. import CTRfilm
from ..CTRdistributions import PoissonProfile
from . import _poisson_oracle


class PoissonFlatHeightMixin(unittest.TestCase):
    """Shared coordinates and builders for flat-height tests."""

    H = np.zeros(5)
    K = np.zeros(5)
    L = np.array([0.35, 0.8, 1.3, 1.85, 2.4])

    def bound_surface(self, profile, n_layers=2, **keyargs):
        """Return a stacked crystal, its surface, and the base width."""
        w_base = _poisson_oracle.minimum_base_width(profile)
        crystal, surface = _poisson_oracle.poisson_crystal(
            profile, w_base=w_base, n_layers=n_layers, **keyargs
        )
        crystal.apply_stacking()
        return crystal, surface, w_base


class TestFlatHeightCorrections(PoissonFlatHeightMixin):
    """State corrections reproduce independently built flat-height crystals."""

    def test_state_corrections_match_the_flat_height_oracle(self):
        """Adding one state correction to the base crystal gives that height.

        The base crystal is the same sharp Film boundary the correction is
        defined against, so ``common + correction`` must reproduce a crystal
        independently stacked at that flat height.
        """
        cases = {
            "integer growth": PoissonProfile(2.0, alpha=0.0),
            "fractional step": PoissonProfile(1.5, alpha=0.0),
            "negative offset": PoissonProfile(0.0, alpha=0.0, offset=-2.25),
            "poisson growth": PoissonProfile(1.5, alpha=0.6),
            "poisson etching": PoissonProfile(-1.5, alpha=0.6),
            "step and poisson": PoissonProfile(2.3, alpha=0.4, offset=0.4),
        }
        for label, profile in cases.items():
            with self.subTest(case=label):
                profile.tail_probability = 1e-14
                _, surface, w_base = self.bound_surface(profile)
                common = _poisson_oracle.flat_height_crystal(
                    -1, w_base=w_base
                ).F(self.H, self.K, self.L)
                states = surface.flat_domain_corrections(
                    self.H, self.K, self.L
                )
                self.assertGreater(states.layer_numbers.size, 0)
                for state in states.iter_states():
                    oracle = _poisson_oracle.flat_height_crystal(
                        state.layer_number, w_base=w_base
                    ).F(self.H, self.K, self.L)
                    np.testing.assert_allclose(
                        common + state.amplitude, oracle, rtol=1e-12
                    )

    def test_multi_layer_cycles_and_termination_banks(self):
        """Layer cycles and an explicit termination bank keep the agreement."""
        profile = PoissonProfile(1.5, alpha=0.5, tail_probability=1e-14)
        for n_layers in (1, 2, 3):
            with self.subTest(n_layers=n_layers):
                _, surface, w_base = self.bound_surface(
                    profile, n_layers=n_layers
                )
                common = _poisson_oracle.flat_height_crystal(
                    -1, w_base=w_base, n_layers=n_layers
                ).F(self.H, self.K, self.L)
                states = surface.flat_domain_corrections(
                    self.H, self.K, self.L
                )
                for state in states.iter_states():
                    oracle = _poisson_oracle.flat_height_crystal(
                        state.layer_number,
                        w_base=w_base,
                        n_layers=n_layers,
                    ).F(self.H, self.K, self.L)
                    np.testing.assert_allclose(
                        common + state.amplitude, oracle, rtol=1e-12
                    )

    def test_copied_surfaces_still_produce_corrections(self):
        """``copy.deepcopy`` is performed once per fit and must survive it."""
        profile = PoissonProfile(1.5, alpha=0.5, tail_probability=1e-14)
        crystal, surface, w_base = self.bound_surface(profile)
        original = surface.flat_domain_corrections(
            self.H, self.K, self.L
        ).as_array()

        copied_crystal = copy.deepcopy(crystal)
        copied_crystal.apply_stacking()
        copied = copied_crystal.uc_surface_list[1].flat_domain_corrections(
            self.H, self.K, self.L
        ).as_array()
        np.testing.assert_allclose(copied, original, rtol=1e-12)

    def test_parameter_changes_are_picked_up(self):
        """Corrections follow the basis through the ``F_uc`` sync path."""
        profile = PoissonProfile(1.5, alpha=0.0, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(profile)
        before = surface.flat_domain_corrections(self.H, self.K, self.L)

        surface.basis[0] = 3.0
        after = surface.flat_domain_corrections(self.H, self.K, self.L)
        self.assertNotEqual(
            before.layer_numbers.tolist(), after.layer_numbers.tolist()
        )


class TestFlatHeightSelection(PoissonFlatHeightMixin):
    """The ensemble selection policy and its reported masses."""

    def test_interior_bins_agree_with_the_shifted_probability(self):
        """``surface_occupancy(n)`` is ``probability(n + 1)`` except at the top.

        Layer ``n`` is the top filled layer exactly when the height change is
        ``n + 1``. The two expressions therefore agree on every interior bin,
        and deliberately disagree at the terminal one, where
        ``surface_occupancy`` folds in the whole upper tail. Both halves are
        asserted so a refactor cannot silently exchange them.
        """
        profile = PoissonProfile(1.5, alpha=0.6, tail_probability=1e-6)
        _, surface, _ = self.bound_surface(profile)
        candidates = surface._profile_candidates()

        np.testing.assert_allclose(
            candidates.surface_occupancy[:-1],
            candidates.exposed_probability[:-1],
            atol=1e-15,
        )
        self.assertGreater(
            abs(
                candidates.surface_occupancy[-1]
                - candidates.exposed_probability[-1]
            ),
            0.0,
        )

    def test_the_two_retention_policies_select_different_sets(self):
        """Coherent retention keeps Film-correction layers the ensemble drops.

        A single deterministic grown state exposes only its own layer, but the
        coherent assembly must also keep the Film layer it added underneath.
        Merging the two policies into one mask would silently give one path
        the other's cutoff.
        """
        profile = PoissonProfile(2.0, alpha=0.0, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(profile)
        candidates = surface._profile_candidates()

        coherent = candidates.layer_numbers[
            surface._coherent_retention(candidates)
        ]
        low, high, _, _, _, _ = surface._height_state_selection(candidates, 10)
        ensemble = candidates.layer_numbers[low : high + 1]

        self.assertEqual(ensemble.tolist(), [1])
        self.assertEqual(coherent.tolist(), [0, 1])

    def test_exact_layer_count_is_a_policy_not_a_cap(self):
        """Below the count every state is kept; above it the mass decides."""
        narrow = PoissonProfile(1.5, alpha=0.0, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(narrow)
        states = surface.flat_domain_corrections(
            self.H, self.K, self.L, exact_layer_count=10
        )
        self.assertEqual(states.layer_numbers.size, 2)

        wide = PoissonProfile(6.0, alpha=1.0, tail_probability=1e-6)
        _, surface, _ = self.bound_surface(wide)
        candidates = surface._profile_candidates()
        populated = int(np.sum(candidates.exposed_probability > 0.0))
        self.assertGreater(populated, 10)

        # More than ten states are retained because the probability-mass
        # criterion demands them, and asking for fewer does not cap it.
        for requested in (3, 10):
            with self.subTest(exact_layer_count=requested):
                states = surface.flat_domain_corrections(
                    self.H, self.K, self.L, exact_layer_count=requested
                )
                self.assertGreater(states.layer_numbers.size, 10)
                self.assertLessEqual(
                    states.excluded_probability, wide.tail_probability
                )

        everything = surface.flat_domain_corrections(
            self.H, self.K, self.L, exact_layer_count=populated + 5
        )
        self.assertEqual(everything.layer_numbers.size, populated)

    def test_masses_are_normalized_once_and_add_to_one(self):
        """Retained and excluded masses partition the whole distribution."""
        for tail in (1e-14, 1e-3):
            with self.subTest(tail_probability=tail):
                profile = PoissonProfile(2.0, alpha=1.0, tail_probability=tail)
                _, surface, _ = self.bound_surface(profile)
                states = surface.flat_domain_corrections(
                    self.H, self.K, self.L
                )
                self.assertAlmostEqual(
                    float(states.probabilities.sum()), 1.0, places=12
                )
                self.assertAlmostEqual(
                    states.raw_retained_probability
                    + states.excluded_lower_probability
                    + states.excluded_upper_probability,
                    1.0,
                    places=12,
                )
                self.assertAlmostEqual(
                    states.excluded_probability,
                    states.excluded_lower_probability
                    + states.excluded_upper_probability,
                    places=15,
                )

    def test_excluded_mass_is_not_an_amplitude_error_bound(self):
        """A looser tail costs more amplitude error than dimensionless mass.

        This is why the design record refuses to use the excluded probability
        directly as a reconstruction tolerance.
        """
        profile = PoissonProfile(2.0, alpha=1.0, tail_probability=1e-3)
        _, surface, _ = self.bound_surface(profile)
        direct = surface.F_uc(self.H, self.K, self.L)
        states = surface.flat_domain_corrections(self.H, self.K, self.L)
        rebuilt = sum(
            state.probability * state.amplitude
            for state in states.iter_states()
        )
        relative = np.max(np.abs(direct - rebuilt)) / np.max(np.abs(direct))
        self.assertGreater(states.excluded_probability, 0.0)
        self.assertGreater(relative, states.excluded_probability)

    def test_tight_tail_reconstructs_the_coherent_amplitude(self):
        """With a negligible tail the state sum reproduces ``F_uc``."""
        for label, profile in {
            "growth": PoissonProfile(1.5, alpha=0.6, tail_probability=1e-14),
            "etching": PoissonProfile(-1.5, alpha=0.6, tail_probability=1e-14),
        }.items():
            with self.subTest(case=label):
                _, surface, _ = self.bound_surface(profile)
                direct = surface.F_uc(self.H, self.K, self.L)
                states = surface.flat_domain_corrections(
                    self.H, self.K, self.L
                )
                rebuilt = sum(
                    state.probability * state.amplitude
                    for state in states.iter_states()
                )
                np.testing.assert_allclose(rebuilt, direct, rtol=1e-11)


class TestFlatHeightContract(PoissonFlatHeightMixin):
    """API behavior of the returned result."""

    def test_evaluation_leaves_the_coherent_assembly_untouched(self):
        """Flat states must not perturb the stored coherent domains."""
        profile = PoissonProfile(1.5, alpha=0.5, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(profile)
        before = surface.F_uc(self.H, self.K, self.L)
        snapshot = [
            (
                [np.copy(matrix) for matrix in uc.coherentDomainMatrix],
                list(uc.coherentDomainOccupancy),
            )
            for uc in surface.layer_ucs + surface.film_layer_ucs
        ]

        surface.flat_domain_corrections(
            self.H, self.K, self.L
        ).as_array()

        for uc, (matrices, occupancies) in zip(
            surface.layer_ucs + surface.film_layer_ucs, snapshot
        ):
            self.assertEqual(len(uc.coherentDomainMatrix), len(matrices))
            for stored, expected in zip(uc.coherentDomainMatrix, matrices):
                np.testing.assert_array_equal(stored, expected)
            np.testing.assert_allclose(
                uc.coherentDomainOccupancy, occupancies
            )
        np.testing.assert_array_equal(
            surface.F_uc(self.H, self.K, self.L), before
        )

    def test_as_array_matches_the_streamed_states(self):
        """The diagnostic materialization agrees with the streaming path."""
        profile = PoissonProfile(1.5, alpha=0.5, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(profile)
        states = surface.flat_domain_corrections(self.H, self.K, self.L)

        streamed = {
            state.layer_number: state.amplitude
            for state in states.iter_states()
        }
        materialized = states.as_array()
        self.assertEqual(
            materialized.shape, (states.layer_numbers.size, self.L.size)
        )
        for index, layer in enumerate(states.layer_numbers):
            np.testing.assert_array_equal(
                materialized[index], streamed[int(layer)]
            )

    def test_reported_arrays_are_read_only(self):
        """Callers must not be able to edit the reported selection."""
        profile = PoissonProfile(1.5, alpha=0.5, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(profile)
        states = surface.flat_domain_corrections(self.H, self.K, self.L)
        with self.assertRaises(ValueError):
            states.layer_numbers[0] = 0
        with self.assertRaises(ValueError):
            states.probabilities[0] = 0.0

    def test_requires_a_bound_film(self):
        """A surface which is not stacked has no boundary to correct."""
        surface = CTRfilm.PoissonSurface(
            _poisson_oracle.layered_cell(2, "film"),
            profile=PoissonProfile(1.5, alpha=0.5),
            name="surface",
        )
        with self.assertRaisesRegex(ValueError, "stacked immediately above"):
            surface.flat_domain_corrections(self.H, self.K, self.L)

    def test_rejects_an_invalid_state_count(self):
        """``exact_layer_count`` is a positive integer policy setting."""
        profile = PoissonProfile(1.5, alpha=0.5, tail_probability=1e-14)
        _, surface, _ = self.bound_surface(profile)
        for bad in (0, -1, 2.5):
            with self.subTest(exact_layer_count=bad):
                with self.assertRaisesRegex(ValueError, "positive integer"):
                    surface.flat_domain_corrections(
                        self.H, self.K, self.L, exact_layer_count=bad
                    )


if __name__ == "__main__":
    unittest.main()

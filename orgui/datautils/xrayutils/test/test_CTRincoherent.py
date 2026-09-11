"""Contract tests for the incoherent CTR wrapper.

Every model here is synthetic. Increment 3 of
``doc/design/incoherent_ctr_models.md`` requires the contract to hold without
the Poisson implementation, so nothing in this file constructs or imports a
``PoissonSurface``.
"""

import copy
import unittest

import numpy as np

from .. import CTRcalc, CTRincoherent
from ..CTRincoherent import (
    CoherentStateEnsembleModel,
    IncoherentF2Model,
    IncoherentState,
    available_incoherent_models,
    create_incoherent_model,
    create_incoherent_model_from_config,
    register_incoherent_model,
    unregister_incoherent_model,
)
from . import _poisson_oracle


def _crystal(fitted=0):
    """Return a plain coherent crystal with ``fitted`` crystal parameters."""
    film = CTRcalc.UnitCell([3.0, 3.0, 6.0], [90.0, 90.0, 90.0], name="film")
    film.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=0)
    crystal = CTRcalc.SXRDCrystal(
        _poisson_oracle.layered_cell(1, "bulk"), film, stacking=np.array([1])
    )
    # Crystal fit parameters come from the component cell, which is how a
    # real crystal acquires its tail.
    definitions = (("iDW", (0.0, 5.0)), ("occ", (0.0, 1.0)))
    for index in range(fitted):
        parameter, limits = definitions[index]
        crystal["film"].addFitParameter(
            (0, parameter), limits=limits, name=f"film {parameter}"
        )
    return crystal


class _NoLocalModel(IncoherentF2Model):
    """Wrapper with no local fit parameters."""

    model_type = "test_no_local"

    def _evaluate_F2(self, context):  # noqa: N802
        return np.abs(context.coherent.total) ** 2


class _OneLocalModel(IncoherentF2Model):
    """Wrapper with one local fit parameter."""

    model_type = "test_one_local"
    parameterLookup = {"scale": 0}
    parameterLookup_inv = {0: "scale"}

    def __init__(self, crystal, *, scale=0.5, name="one"):
        super().__init__(crystal, name=name)
        self.basis = np.array([float(scale)])
        self.basis_0 = np.array([float(scale)])
        self.errors = None

    def _config_settings(self):
        return {}

    def _evaluate_F2(self, context):  # noqa: N802
        return self.basis[0] * np.abs(context.coherent.total) ** 2


class _ThreeLocalModel(IncoherentF2Model):
    """Wrapper with three local fit parameters."""

    model_type = "test_three_local"
    parameterLookup = {"a": 0, "b": 1, "c": 2}
    parameterLookup_inv = {0: "a", 1: "b", 2: "c"}

    def __init__(self, crystal, *, name="three"):
        super().__init__(crystal, name=name)
        self.basis = np.array([0.1, 0.2, 0.3])
        self.basis_0 = np.copy(self.basis)
        self.errors = None

    def _evaluate_F2(self, context):  # noqa: N802
        return float(np.sum(self.basis)) * np.abs(context.coherent.total) ** 2


class _TwoStateSource:
    """State source standing in for a real surface decomposition."""

    def __init__(self, amplitudes, probabilities):
        self._amplitudes = amplitudes
        self._probabilities = probabilities

    def iter_states(self):
        """Yield the states one at a time."""
        for index, (amplitude, probability) in enumerate(
            zip(self._amplitudes, self._probabilities)
        ):
            yield IncoherentState(index, probability, amplitude)


class _TwoStateEnsemble(CoherentStateEnsembleModel):
    """Ensemble of two explicitly supplied complete state amplitudes."""

    model_type = "test_two_state"

    def __init__(self, crystal, *, kappa=0.0, offsets=(5.0, -5.0), name="two"):
        super().__init__(crystal, name=name)
        self.basis = np.array([float(kappa)])
        self.basis_0 = np.array([float(kappa)])
        self.errors = None
        self.offsets = tuple(offsets)
        self.probabilities = (0.5, 0.5)

    parameterLookup = {"incoherent_fraction": 0}
    parameterLookup_inv = {0: "incoherent_fraction"}

    @property
    def incoherent_fraction(self):
        """Return the mixing fraction held in the local basis."""
        return float(self.basis[0])

    def state_amplitudes(self, context):
        """Return the two complete state amplitudes for this context."""
        total = context.coherent.total
        return [total + offset for offset in self.offsets]

    def _iter_states(self, context):
        source = _TwoStateSource(
            self.state_amplitudes(context), self.probabilities
        )
        return source.iter_states()


class ContractMixin(unittest.TestCase):
    """Shared coordinates for the contract tests."""

    H = np.zeros(4)
    K = np.zeros(4)
    L = np.array([0.4, 0.9, 1.5, 2.2])


class TestParameterBlockContract(ContractMixin):
    """The composite fit API keeps one stable local-then-coherent order."""

    def models(self):
        """Yield one wrapper of each local-parameter count."""
        no_local = _NoLocalModel(_crystal(fitted=2))

        one_local = _OneLocalModel(_crystal(fitted=2))
        one_local.addFitParameter("scale", limits=(0.0, 1.0))

        three_local = _ThreeLocalModel(_crystal(fitted=1))
        for parameter in ("a", "b", "c"):
            three_local.addFitParameter(parameter, limits=(0.0, 1.0))

        return {
            "zero local": (no_local, 0, 2),
            "one local": (one_local, 1, 2),
            "three local": (three_local, 3, 1),
        }

    def test_every_contract_array_has_the_same_length_and_order(self):
        """One check catches any composite override left un-overridden.

        ``LinearFitFunctions`` implements these as local-only methods, so a
        forgotten override does not raise: it returns a vector one block
        short. Comparing every array at once is what makes that visible.
        """
        for label, (model, n_local, n_coherent) in self.models().items():
            with self.subTest(case=label):
                total = n_local + n_coherent
                start, lower, upper = model.getStartParamAndLimits()
                self.assertEqual(len(model.fitparnames), total)
                self.assertEqual(len(model.getInitialParameters()), total)
                self.assertEqual(len(model.parameter_list()), total)
                self.assertEqual(len(start), total)
                self.assertEqual(len(lower), total)
                self.assertEqual(len(upper), total)
                self.assertEqual(model.n_local_parameters, n_local)
                self.assertEqual(model.n_coherent_parameters, n_coherent)

                # The coherent block is the tail, in the crystal's own order.
                crystal_names = list(model.coherent_model.fitparnames)
                self.assertEqual(
                    model.fitparnames[n_local:], crystal_names
                )
                np.testing.assert_allclose(
                    model.getInitialParameters()[n_local:],
                    model.coherent_model.getInitialParameters(),
                )

    def test_priors_cover_both_blocks(self):
        """Priors follow the same order and length as the names."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        self.assertEqual(len(model.priors), len(model.fitparnames))
        # Prior entries may be arrays, so compare them element-wise.
        for mine, theirs in zip(
            model.priors[1:], list(model.coherent_model.priors)
        ):
            np.testing.assert_array_equal(mine, theirs)
        np.testing.assert_array_equal(model.priors[0], (0.0, 1.0))

    def test_set_parameters_splits_the_vector(self):
        """Each block receives its own segment."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))

        model.setParameters([0.75, 1.5, 2.5])
        self.assertAlmostEqual(model.basis[0], 0.75)
        np.testing.assert_allclose(
            model.coherent_model.getInitialParameters(), [1.5, 2.5]
        )
        np.testing.assert_allclose(
            model.getInitialParameters(), [0.75, 1.5, 2.5]
        )

    def test_bad_lengths_leave_the_model_untouched(self):
        """The whole vector is validated before either block is written."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        before = model.getInitialParameters().copy()

        for bad in ([0.5], [0.5, 1.0], [0.5, 1.0, 2.0, 3.0]):
            with self.subTest(length=len(bad)):
                with self.assertRaisesRegex(ValueError, "entries but this model"):
                    model.setParameters(bad)
                np.testing.assert_allclose(
                    model.getInitialParameters(), before
                )

    def test_local_only_setter_refuses(self):
        """``setFitParameters`` would silently drop the coherent tail."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        with self.assertRaisesRegex(NotImplementedError, "setParameters"):
            model.setFitParameters([0.5])

    def test_errors_round_trip_across_both_blocks(self):
        """Errors split like parameters, and clear together."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))

        model.setFitErrors([0.01, 0.02, 0.03])
        np.testing.assert_allclose(
            model.getFitErrors(), [0.01, 0.02, 0.03]
        )

        model.setFitErrors(None)
        with self.assertRaises(ValueError):
            model.getFitErrors()

    def test_empty_local_block_is_skipped_not_queried(self):
        """A block with no parameters must not raise for missing errors."""
        model = _NoLocalModel(_crystal(fitted=2))
        model.setFitErrors([0.05, 0.06])
        np.testing.assert_allclose(model.getFitErrors(), [0.05, 0.06])

    def test_bounds_reach_both_blocks(self):
        """``setLimits`` splits like the parameter vector."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        model.setLimits([(0.1, 0.9), (0.2, 4.0), (0.3, 4.5)])
        _, lower, upper = model.getStartParamAndLimits()
        np.testing.assert_allclose(lower, [0.1, 0.2, 0.3])
        np.testing.assert_allclose(upper, [0.9, 4.0, 4.5])

    def test_wrapping_requires_a_crystal(self):
        """The wrapper owns one primary coherent model."""
        with self.assertRaisesRegex(TypeError, "wraps one SXRDCrystal"):
            _NoLocalModel(object())


class TestSerializationAndRegistry(ContractMixin):
    """Copying, configuration round-trips, and registry behavior."""

    def test_deepcopy_keeps_both_blocks(self):
        """``CTROptimizer`` deep-copies the model once per fit."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        model.setParameters([0.75, 1.5, 2.5])

        copied = copy.deepcopy(model)
        np.testing.assert_allclose(
            copied.getInitialParameters(), model.getInitialParameters()
        )
        self.assertEqual(copied.fitparnames, model.fitparnames)

        copied.setParameters([0.25, 1.0, 2.0])
        np.testing.assert_allclose(
            model.getInitialParameters(), [0.75, 1.5, 2.5]
        )

    def test_parameter_dict_round_trip(self):
        """Both subtrees survive a dictionary round-trip."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        model.setParameters([0.75, 1.5, 2.5])
        stored = model.parametersToDict()
        self.assertEqual(set(stored), {"local", "coherent"})

        restored = _OneLocalModel(_crystal(fitted=2))
        restored.parametersFromDict(stored)
        self.assertEqual(restored.fitparnames, model.fitparnames)
        np.testing.assert_allclose(
            restored.getInitialParameters(), model.getInitialParameters()
        )

    def test_clear_parameters_clears_both(self):
        """Clearing is composite, like the other parameter operations."""
        model = _OneLocalModel(_crystal(fitted=2))
        model.addFitParameter("scale", limits=(0.0, 1.0))
        model.clearParameters()
        self.assertEqual(model.fitparnames, [])
        self.assertEqual(model.n_local_parameters, 0)
        self.assertEqual(model.n_coherent_parameters, 0)

    def test_config_round_trip_restores_only_local_state(self):
        """The factory wraps a supplied crystal and restores local state."""
        register_incoherent_model(_OneLocalModel)
        try:
            model = _OneLocalModel(_crystal(fitted=2))
            model.addFitParameter("scale", limits=(0.0, 1.0))
            model.setParameters([0.75, 1.5, 2.5])
            config = model.to_config()
            self.assertEqual(config["type"], "test_one_local")
            self.assertEqual(set(config), {"type", "settings", "parameters"})

            crystal = _crystal(fitted=2)
            restored = create_incoherent_model_from_config(crystal, config)
            self.assertIs(restored.coherent_model, crystal)
            self.assertEqual(
                restored.fitparnames[: restored.n_local_parameters],
                ["one scale"],
            )
            self.assertAlmostEqual(restored.basis[0], 0.75)
        finally:
            unregister_incoherent_model("test_one_local")

    def test_registry_rejects_bad_entries(self):
        """Only concrete, keyed, quantity-bearing subclasses register."""
        register_incoherent_model(_NoLocalModel)
        try:
            self.assertIn("test_no_local", available_incoherent_models())
            info = available_incoherent_models()["test_no_local"]
            self.assertIs(info.model_class, _NoLocalModel)
            self.assertEqual(info.output_quantity, "F2")
            self.assertIn("kinematical", info.supported_forward_models)

            with self.assertRaisesRegex(ValueError, "already registered"):
                register_incoherent_model(_NoLocalModel)

            with self.assertRaisesRegex(TypeError, "IncoherentModel"):
                register_incoherent_model(dict)

            with self.assertRaisesRegex(TypeError, "abstract"):
                register_incoherent_model(IncoherentF2Model)

            class _Unkeyed(_NoLocalModel):
                model_type = ""

            with self.assertRaisesRegex(ValueError, "model_type"):
                register_incoherent_model(_Unkeyed)
        finally:
            unregister_incoherent_model("test_no_local")

    def test_registry_mapping_is_read_only(self):
        """Discovery returns a view, not the live registry."""
        models = available_incoherent_models()
        with self.assertRaises(TypeError):
            models["anything"] = None

    def test_unknown_keys_fail_deterministically(self):
        """Deserialization instantiates registered keys only."""
        with self.assertRaisesRegex(ValueError, "Unknown incoherent model"):
            create_incoherent_model("not_registered", _crystal())
        with self.assertRaisesRegex(ValueError, "Unknown incoherent model"):
            create_incoherent_model_from_config(
                _crystal(), {"type": "not_registered"}
            )
        with self.assertRaisesRegex(ValueError, "Unexpected keys"):
            create_incoherent_model_from_config(
                _crystal(), {"type": "x", "surprise": 1}
            )

    def test_third_party_models_need_no_core_changes(self):
        """A model defined outside the package registers and evaluates."""

        @register_incoherent_model
        class _ThirdParty(IncoherentF2Model):
            """Externally defined model."""

            model_type = "test_third_party"

            def _evaluate_F2(self, context):  # noqa: N802
                return 2.0 * np.abs(context.coherent.total) ** 2

        try:
            model = create_incoherent_model("test_third_party", _crystal())
            np.testing.assert_allclose(
                model.F2(self.H, self.K, self.L),
                2.0 * model.coherent_model.F2(self.H, self.K, self.L),
                rtol=1e-14,
            )
        finally:
            unregister_incoherent_model("test_third_party")


class TestF2Boundary(ContractMixin):
    """The squared-structure-factor boundary and its guarantees."""

    def test_F2_is_real_and_nonnegative(self):
        """The public boundary validates what it returns."""
        model = _NoLocalModel(_crystal())
        value = model.F2(self.H, self.K, self.L)
        self.assertEqual(value.dtype, np.float64)
        self.assertTrue(np.all(value >= 0.0))
        np.testing.assert_allclose(
            value, model.coherent_model.F2(self.H, self.K, self.L), rtol=1e-14
        )

    def test_no_F_method_is_exposed(self):
        """A mixed state has no unique complex amplitude."""
        model = _NoLocalModel(_crystal())
        self.assertFalse(hasattr(model, "F"))

    def test_negative_or_nonfinite_results_are_rejected(self):
        """A broken model fails at the boundary rather than downstream."""

        class _Negative(IncoherentF2Model):
            """Deliberately invalid model."""

            model_type = "test_negative"

            def _evaluate_F2(self, context):  # noqa: N802
                return -np.ones_like(np.asarray(context.coordinates[2]))

        class _NotFinite(IncoherentF2Model):
            """Deliberately invalid model."""

            model_type = "test_not_finite"

            def _evaluate_F2(self, context):  # noqa: N802
                return np.full_like(
                    np.asarray(context.coordinates[2]), np.nan
                )

        with self.assertRaisesRegex(ValueError, "negative F2"):
            _Negative(_crystal()).F2(self.H, self.K, self.L)
        with self.assertRaisesRegex(ValueError, "non-finite F2"):
            _NotFinite(_crystal()).F2(self.H, self.K, self.L)

    def test_validate_rejects_unsupported_forward_models(self):
        """An F2 wrapper must not be reinterpreted as a reflectivity."""
        model = _NoLocalModel(_crystal())
        model.validate(forward_model="kinematical")
        with self.assertRaisesRegex(ValueError, "does not support"):
            model.validate(forward_model="dwba")


class TestEnsembleMixing(ContractMixin):
    """Endpoints and interpolation of the generic state ensemble."""

    def ensemble(self, kappa, offsets=(5.0, -5.0)):
        """Return a two-state ensemble with the given mixing fraction."""
        return _TwoStateEnsemble(_crystal(), kappa=kappa, offsets=offsets)

    def test_coherent_endpoint_matches_the_wrapped_crystal(self):
        """``kappa = 0`` is the live coherent result, exactly."""
        model = self.ensemble(0.0)
        np.testing.assert_allclose(
            model.F2(self.H, self.K, self.L),
            model.coherent_model.F2(self.H, self.K, self.L),
            rtol=1e-14,
        )

    def test_incoherent_endpoint_is_the_squared_average(self):
        """``kappa = 1`` averages the squared state amplitudes."""
        model = self.ensemble(1.0)
        context = CTRincoherent.KinematicIncoherentContext(
            model.coherent_model, self.H, self.K, self.L
        )
        amplitudes = model.state_amplitudes(context)
        expected = sum(
            0.5 * np.abs(amplitude) ** 2 for amplitude in amplitudes
        )
        np.testing.assert_allclose(
            model.F2(self.H, self.K, self.L), expected, rtol=1e-12
        )

    def test_partial_mixing_is_the_convex_interpolation(self):
        """Intermediate fractions interpolate the two endpoints linearly."""
        coherent = self.ensemble(0.0).F2(self.H, self.K, self.L)
        incoherent = self.ensemble(1.0).F2(self.H, self.K, self.L)
        for kappa in (0.25, 0.5, 0.75):
            with self.subTest(kappa=kappa):
                np.testing.assert_allclose(
                    self.ensemble(kappa).F2(self.H, self.K, self.L),
                    (1.0 - kappa) * coherent + kappa * incoherent,
                    rtol=1e-12,
                )

    def test_the_endpoints_actually_differ(self):
        """Otherwise the interpolation test would pass on a constant."""
        coherent = self.ensemble(0.0).F2(self.H, self.K, self.L)
        incoherent = self.ensemble(1.0).F2(self.H, self.K, self.L)
        self.assertGreater(
            np.max(np.abs(incoherent - coherent)) / np.max(incoherent), 0.01
        )
        # The incoherent endpoint drops a cross term, so it never falls below.
        self.assertTrue(np.all(incoherent >= coherent - 1e-9))

    def test_a_single_state_is_independent_of_kappa(self):
        """With one populated height there is no variance to expose."""

        class _SingleState(_TwoStateEnsemble):
            """One-state ensemble."""

            model_type = "test_single_state"

            def __init__(self, crystal, *, kappa=0.0):
                super().__init__(crystal, kappa=kappa, offsets=(0.0,))
                self.probabilities = (1.0,)

        values = [
            _SingleState(_crystal(), kappa=kappa).F2(self.H, self.K, self.L)
            for kappa in (0.0, 0.5, 1.0)
        ]
        for value in values[1:]:
            np.testing.assert_allclose(value, values[0], rtol=1e-12)

    def test_unnormalized_masses_are_rejected(self):
        """Probabilities must be normalized over the retained states once."""

        class _Unnormalized(_TwoStateEnsemble):
            """Ensemble whose masses do not sum to one."""

            model_type = "test_unnormalized"

            def __init__(self, crystal):
                super().__init__(crystal, kappa=1.0)
                self.probabilities = (0.5, 0.2)

        with self.assertRaisesRegex(ValueError, "sum to"):
            _Unnormalized(_crystal()).F2(self.H, self.K, self.L)

    def test_invalid_mixing_fractions_are_rejected(self):
        """``kappa`` is dimensionless and lives in ``[0, 1]``."""
        for kappa in (-0.1, 1.5, np.nan):
            with self.subTest(kappa=kappa):
                model = self.ensemble(1.0)
                model.basis[0] = kappa
                with self.assertRaisesRegex(ValueError, "incoherent_fraction"):
                    model.F2(self.H, self.K, self.L)

    def test_reconstruction_diagnostic_is_available(self):
        """The state sum is a diagnostic against the live coherent result."""
        model = self.ensemble(1.0, offsets=(2.0, -2.0))
        self.assertAlmostEqual(
            model.coherence_reconstruction_error(self.H, self.K, self.L),
            0.0,
            places=12,
        )


class TestEvaluationContext(ContractMixin):
    """The lazy decomposition shared by every model."""

    def test_the_coherent_decomposition_is_evaluated_once(self):
        """Streaming many states must not re-evaluate the common crystal."""
        crystal = _crystal()
        calls = []
        original = crystal.evaluate_kinematic

        def counting(*args, **keyargs):
            calls.append(args)
            return original(*args, **keyargs)

        crystal.evaluate_kinematic = counting
        context = CTRincoherent.KinematicIncoherentContext(
            crystal, self.H, self.K, self.L
        )
        for _ in range(5):
            context.coherent  # noqa: B018
        context.component("film")
        context.common_amplitude("film")
        self.assertEqual(len(calls), 1)

    def test_common_amplitude_removes_only_the_named_component(self):
        """The common part is everything except the target component."""
        crystal = _crystal()
        context = CTRincoherent.KinematicIncoherentContext(
            crystal, self.H, self.K, self.L
        )
        np.testing.assert_allclose(
            context.common_amplitude("film"), context.coherent.bulk, rtol=1e-14
        )

    def test_states_are_evaluated_at_transformed_coordinates(self):
        """An outer domain changes where a component state must be evaluated.

        The context must call the evaluator once per domain at that domain's
        coordinates, and combine the results with the crystal's own area,
        weight, and occupancy scaling.
        """
        crystal = _crystal()
        transforms = [
            (np.identity(3), 0.6),
            (np.diag([1.0, 1.0, 0.5]), 0.4),
        ]
        crystal.setDomain(0, transforms)
        context = CTRincoherent.KinematicIncoherentContext(
            crystal, self.H, self.K, self.L
        )

        seen = []

        def evaluator(component, h, k, l):  # noqa: E741
            seen.append(np.copy(np.asarray(l)))
            return _TwoStateSource(
                [
                    np.asarray(l, dtype=np.complex128),
                    np.zeros_like(l, dtype=np.complex128),
                ],
                [0.5, 0.5],
            )

        states = list(context.iter_component_states("film", evaluator))

        self.assertEqual(len(seen), 2)
        np.testing.assert_allclose(seen[0], self.L)
        np.testing.assert_allclose(seen[1], 0.5 * self.L)

        common = context.common_amplitude("film")
        area_weight = crystal.reference_area / crystal.uc_surface_list[0].uc_area
        expected = common + area_weight * (
            0.6 * self.L + 0.4 * (0.5 * self.L)
        )
        self.assertEqual(len(states), 2)
        np.testing.assert_allclose(states[0].amplitude, expected, rtol=1e-12)
        self.assertEqual(states[0].probability, 0.5)

    def test_inconsistent_state_metadata_is_rejected(self):
        """Only the amplitudes may depend on the domain transform."""
        crystal = _crystal()
        crystal.setDomain(0, [(np.identity(3), 0.5), (np.identity(3), 0.5)])
        context = CTRincoherent.KinematicIncoherentContext(
            crystal, self.H, self.K, self.L
        )
        masses = iter([(0.5, 0.5), (0.25, 0.75)])

        def evaluator(component, h, k, l):  # noqa: E741
            zeros = np.zeros_like(np.asarray(l), dtype=np.complex128)
            return _TwoStateSource([zeros, zeros], next(masses))

        with self.assertRaisesRegex(ValueError, "disagree on state"):
            list(context.iter_component_states("film", evaluator))


if __name__ == "__main__":
    unittest.main()

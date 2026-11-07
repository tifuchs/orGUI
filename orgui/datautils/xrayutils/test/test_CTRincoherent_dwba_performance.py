"""Numerical equivalence and deterministic work-count guards for terrace DWBA."""

import copy
from collections import Counter

import numpy as np
import pytest
from numpy.testing import assert_allclose

from .. import CTRdwba, _CTRincoherent_dwba as ensemble
from ..CTRcrystallinity import CrystallinityProfile
from ..CTRdistributions import PoissonProfile, RecenteredSkellamProfile
from .test_CTRincoherent_dwba import _hkl, _model


def _old_intensity(model, hkl, **options):
    """Reference the original complete-state algorithm, without the fast basis."""
    native = model.coherent_model.dwba.evaluate(*hkl, **options)
    candidates = model.target_surface()._profile_candidates()
    low, high, masses, *_ = model.target_surface()._height_state_selection(
        candidates, model.exact_layer_count)
    amplitude = np.zeros_like(native.reflectivity, dtype=complex)
    intensity = np.zeros_like(native.reflectivity)
    for height, mass in zip(candidates.layer_numbers[low:high + 1], masses):
        if mass == 0:
            continue
        state = ensemble._flat_crystal(model, height).dwba
        with state.batch():
            _, atomic, profile = state._atomic_model_snapshot()
            if model.optical_reference == "shared":
                prepared = ensemble._shared_preparation(
                    state, native.prepared, atomic, profile)
            else:
                prepared = state._prepare_explicit(
                    ensemble._geometry(native.prepared), state._reference(profile),
                    atomic, native.prepared.polarization_i,
                    native.prepared.polarization_f, 1e-7, 1e-8)
            result = state.evaluate_prepared(
                prepared, bulk_attenuation=options.get("bulk_attenuation", 0.0))
        phase = ensemble._phase_to(prepared, native.prepared)
        field = (native.unperturbed_amplitude + native.amplitude_prefactor * (
            result.F_atomic * phase - native.F_reference)
            if model.optical_reference == "shared"
            else result.total_amplitude * phase)
        amplitude += mass * field
        intensity += mass * abs(field)**2
    coherent = (native.reflectivity if model.optical_reference == "shared"
                else abs(amplitude)**2)
    kappa = model.incoherent_fraction
    return (1 - kappa)*coherent + kappa*intensity


def _work_counts(monkeypatch):
    counts = Counter()
    for method, counter in (("_atomic_model_snapshot", "snapshot_builds"),
                            ("_pack_atomic_records", "packing_builds")):
        original = getattr(CTRdwba.DWBAState, method)

        def _count(self, *args, _original=original, _counter=counter, **kwargs):
            before = self._cache_counters[_counter]
            result = _original(self, *args, **kwargs)
            counts[_counter] += self._cache_counters[_counter] - before
            return result

        monkeypatch.setattr(CTRdwba.DWBAState, method, _count)
    for owner, name, key in ((CTRdwba, "solve_electric_field", "field_solves"),
                             (ensemble, "_flat_crystal", "crystal_copies"),
                             (CTRdwba.DWBAState, "_evaluate_raw", "native_calls")):
        original = getattr(owner, name)

        def _count(*args, _original=original, _key=key, **kwargs):
            counts[_key] += 1
            return _original(*args, **kwargs)

        monkeypatch.setattr(owner, name, _count)
    return counts


@pytest.mark.parametrize("profile", [PoissonProfile(-.7), PoissonProfile(.7),
                                    RecenteredSkellamProfile(.2, .4, .3)])
@pytest.mark.parametrize("channels", [("s", "s"), ("p", "p"), ("s", "p")])
def test_shared_basis_matches_complete_states(profile, channels):
    """Growth, etching, Skellam, specular and Bragg points retain interference."""
    model = _model(profile)
    options = dict(polarization_i=channels[0], polarization_f=channels[1],
                   bulk_attenuation=.025)
    expected = _old_intensity(model, _hkl(), **options)
    actual = model.dwba.evaluate(*_hkl(), **options).reflectivity
    assert_allclose(actual, expected, rtol=3e-10, atol=2e-15)


def test_shared_basis_preserves_component_weights_and_rotated_domains():
    """Domain phases and unequal component weights remain in each full state."""
    model = _model(PoissonProfile(.6))
    crystal = model.coherent_model
    crystal.uc_bulk.name = "surface"  # Bulk names do not select surface components.
    for cell in (crystal.uc_surface_list[0].unitcell,
                 model.target_surface().unitcell):
        cell.basis[:, 1] = cell.basis_0[:, 1] = .13
        cell.basis[:, 2] = cell.basis_0[:, 2] = .27
    model.target_surface().setFitParameters([])
    crystal.weights[:] = [.8, .65]
    for index in (0, 1):
        crystal.setDomain(index, [(np.eye(3), .7), (np.diag([-1., -1., 1.]), .3)])
    expected = _old_intensity(model, _hkl())
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity,
                    expected, rtol=3e-10, atol=2e-15)


@pytest.mark.parametrize("kind", ["layer", "depth"])
@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_general_states_match_original_with_crystallinity(kind, policy):
    """Refining shared fields preserves state-dependent disordered density."""
    model = _model(policy=policy)
    surface = model.target_surface()
    if kind == "depth":
        surface.setCrystallinityProfile(
            CrystallinityProfile(PoissonProfile(1.3, alpha=.4)))
    else:
        surface.unitcell.setLayerCrystallinity(0, .3)
        surface.unitcell.setLayerCrystallinity(1, .8)
        surface.setFitParameters([])
    options = dict(polarization_i="p", polarization_f="p")
    expected = _old_intensity(model, _hkl(), **options)
    actual = model.dwba.evaluate(*_hkl(), **options).reflectivity
    assert_allclose(actual, expected, rtol=3e-10, atol=2e-15)


def test_partial_order_switches_safely_between_specular_and_offspecular():
    """Off-specular order weighting can use the basis in a mixed-rod batch."""
    model = _model()
    model.target_surface().unitcell.setLayerCrystallinity(0, .4)
    model.target_surface().setFitParameters([])
    rods = [(0., 0., np.array([.01, .03])),
            (1., 0., np.array([.06, .45, 1.999]))]
    expected = [_old_intensity(model, rod) for rod in rods]
    with model.dwba.batch():
        for index in (0, 1, 0):
            assert_allclose(model.dwba.evaluate(*rods[index]).reflectivity,
                            expected[index], rtol=3e-10, atol=2e-15)
        assert model.dwba._scope["basis"] is not None


def test_shared_work_does_not_scale_with_height_count(monkeypatch):
    """64 terraces need two packings and one correction kernel per channel."""
    rate = (5.15 / 3)**2 / 2
    model = _model(RecenteredSkellamProfile(1, rate, rate))
    counts = _work_counts(monkeypatch)
    with model.dwba.batch():
        for ell in (.5, .7):
            model.dwba.polarization_reflectivity(
                1., 0., ell, s_fraction=.6, outgoing="unanalysed")
        heights, _ = model.dwba._selection()
        assert len(heights) == 64
        assert counts["crystal_copies"] == 0
        assert counts["snapshot_builds"] == 1
        assert counts["packing_builds"] == 2
        assert counts["native_calls"] == 16  # two rods, four channels, two kernels
        assert counts["field_solves"] == 6  # two incident and four exit solves
    assert model.dwba._scope is None


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_general_states_reused_across_channels_rods_and_nested_batches(
    monkeypatch, policy
):
    """A batch owns one preparation per state, including the fallback path."""
    model = _model(policy=policy)
    model.target_surface().setCrystallinityProfile(
        CrystallinityProfile(PoissonProfile(1.3, alpha=.4)))
    counts = _work_counts(monkeypatch)
    retained = []
    with model.dwba.batch():
        model.dwba.polarization_reflectivity(
            1., 0., .5, s_fraction=.6, outgoing="unanalysed")
        first = counts.copy()
        retained.extend(state for state, _, _ in model.dwba._scope["states"])
        with model.dwba.batch():
            model.dwba.polarization_reflectivity(
                1., 0., .7, s_fraction=.6, outgoing="unanalysed")
        for key in ("crystal_copies", "snapshot_builds", "packing_builds"):
            assert counts[key] == first[key]
        assert counts["crystal_copies"] == 2
        assert counts["snapshot_builds"] == 3
        assert counts["packing_builds"] == 3
        assert len({id(state._batch_snapshot) for state in retained}) == 2
    assert model.dwba._scope is None
    assert all(state._batch_depth == 0 and state._batch_snapshot is None
               and not state._batch_packings for state in retained)
    assert all(state.cache_info()["field_size"] == 0
               and state.cache_info()["reference_size"] == 0
               and state.cache_info()["prepared_size"] == 0
               for state in retained)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_exception_cleanup_and_parameter_change_invalidate_batch(policy):
    """An aborted evaluation cannot freeze the next fit proposal or a copy."""
    model = _model(policy=policy)
    facade = model.dwba
    with pytest.raises(RuntimeError, match="abort"):
        with facade.batch():
            original = facade.evaluate(*_hkl()).reflectivity
            raise RuntimeError("abort")
    assert facade._scope is None
    assert model.coherent_model.dwba._batch_depth == 0
    model.target_surface().basis[2] = .8
    with facade.batch():
        actual = facade.evaluate(*_hkl()).reflectivity
    fresh = copy.deepcopy(model)
    assert fresh.dwba is not facade
    assert fresh.dwba._model is fresh
    assert fresh.dwba._scope is None
    assert_allclose(actual, fresh.dwba.evaluate(*_hkl()).reflectivity, rtol=1e-12)
    assert not np.allclose(actual, original, rtol=1e-7, atol=1e-20)


def test_terrace_field_storage_released_after_exception():
    """Retained private states cannot keep field tables alive after a batch."""
    model = _model(policy="state_specific")
    retained = []
    with pytest.raises(RuntimeError, match="abort"):
        with model.dwba.batch():
            model.dwba.evaluate(*_hkl())
            retained.extend(state for state, _, _ in model.dwba._scope["states"])
            assert all(state.cache_info()["field_size"] > 0 for state in retained)
            raise RuntimeError("abort")
    assert model.dwba._scope is None
    assert all(state.cache_info()["field_size"] == 0
               and state.cache_info()["reference_size"] == 0
               and state._batch_snapshot is None for state in retained)
    # The wrapped crystal's persistent fields still belong to the user.
    assert model.coherent_model.dwba.cache_info()["field_size"] > 0


def test_optimizer_shares_terraces_across_rods(monkeypatch):
    """The objective batches the ensemble itself, rather than just its crystal."""
    from .. import CTRopt, CTRplotutil
    from .test_CTRopt_dwba import _ctr
    model = _model(policy="state_specific")
    ctrs = CTRplotutil.CTRCollection([
        _ctr(model.coherent_model, hk=hk, s_fraction=.6, outgoing="unanalysed",
             geometry=CTRplotutil.CTRScanGeometry("in", .003))
        for hk in ((1, 0), (0, 1))])
    optimizer = CTRopt.CTROptimizer(model, ctrs)
    optimizer.set_dwba()
    optimizer.prepareFit()
    counts = _work_counts(monkeypatch)
    optimizer.flat_prediction()
    assert counts["crystal_copies"] == 2
    assert counts["snapshot_builds"] == 3
    assert counts["packing_builds"] == 3
    assert optimizer.model.dwba._scope is None

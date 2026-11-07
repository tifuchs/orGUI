"""Batch reuse must preserve physical dependencies and release owned storage."""

import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from .. import CTRuc
from .test_CTRdwba import _synthetic_record_call
from .test_CTRincoherent_dwba import _hkl, _model


@pytest.fixture
def native_arguments(monkeypatch):
    """Capture a nontrivial native problem with signed finite and bulk atoms."""
    args = []
    native = CTRuc._CTRcalc_cpp.unitcell_F_DWBA_records
    with monkeypatch.context() as patch:
        patch.setattr(CTRuc._CTRcalc_cpp, "unitcell_F_DWBA_records",
                      lambda *values: args.extend(values))
        _synthetic_record_call(
            q_parallel=[.17, .08], kz_i=[-1 + .02j, -1.3 + .05j],
            kz_f=[-1.1 + .03j, -1.4 + .06j],
            A_i_plus=[1., .8], A_i_minus=[.1, .2],
            A_f_plus=[1., .7], A_f_minus=[.2, .1],
            basis=np.array([[6., 0., 0., 0., .2, .3, .7]]),
            factors=np.array([[1.] + [0.] * 12]),
            finite_positions=np.array([[.13, .27, -.5], [.4, .3, -.5]]),
            finite_records=np.array([1, 1]), finite_weights=np.array([.8, -.3]),
            is_specular=True, reference_shares=np.full((2, 2, 2), 1e-7))
    args[35:39] = [np.array([[.2, .1, -.3]]), np.array([0], dtype=np.int64),
                   np.array([.9]), np.array([-.1])]
    return native, args


def test_native_reuses_finite_work_across_channels_and_damping(native_arguments):
    """All four polarizations and both damping controls match uncached sums."""
    native, args = native_arguments
    cache = CTRuc._CTRcalc_cpp._DWBAFiniteCache()
    native(*args, cache)
    first = CTRuc._CTRcalc_cpp._dwba_motif_work()
    assert first["finite_cache_hits"] == 0
    for attenuation in (0., .005):
        for pi, pf in (("s", "s"), ("s", "p"), ("p", "s"), ("p", "p")):
            args[23:25] = pi, pf
            args[28] = attenuation
            expected = native(*args)
            actual = native(*args, cache)
            work = CTRuc._CTRcalc_cpp._dwba_motif_work()
            assert work["finite_cache_hits"] == 1
            assert work["instances"] == first["instances"] - 2
            assert_allclose(actual, expected, rtol=1e-13, atol=1e-28)
            repeated = native(*args, cache)
            assert_array_equal(repeated, actual)


@pytest.mark.parametrize("index, offset", [(0, 0), (5, 0), (10, 0), (17, 0),
                                         (29, 4), (29, 5), (29, 6),
                                         (30, 0), (31, 2), (34, 0)])
def test_native_exact_dependencies_invalidate(native_arguments, index, offset):
    """Even one-bit changes of geometry, DWF, density or atomic factors miss."""
    native, args = native_arguments
    cache = CTRuc._CTRcalc_cpp._DWBAFiniteCache()
    native(*args, cache)
    changed = args[index].copy()
    if np.iscomplexobj(changed):
        changed.flat[offset] += 1e-12j
    else:
        changed.flat[offset] = np.nextafter(changed.flat[offset], np.inf)
    args[index] = changed
    actual = native(*args, cache)
    assert CTRuc._CTRcalc_cpp._dwba_motif_work()["finite_cache_hits"] == 0
    assert_allclose(actual, native(*args), rtol=1e-13, atol=1e-28)


def test_large_native_request_uses_complete_uncached_sum(native_arguments):
    """The one-MiB cache budget never truncates atoms or changes predictions."""
    native, args = native_arguments
    cache = CTRuc._CTRcalc_cpp._DWBAFiniteCache()
    for index in (31, 32, 33, 34):
        args[index] = np.repeat(args[index], 15000, axis=0)
    expected = native(*args)
    for _ in range(2):
        actual = native(*args, cache)
        work = CTRuc._CTRcalc_cpp._dwba_motif_work()
        assert work["instances"] == 30001
        assert work["finite_cache_hits"] == work["finite_cache_bytes"] == 0
        assert_array_equal(actual, expected)


def test_dense_scan_bypasses_retention_without_dropping_points(native_arguments):
    """The dense kernel evaluates all points, including beyond the cache cutoff."""
    native, args = native_arguments
    cache = CTRuc._CTRcalc_cpp._DWBAFiniteCache()
    for index in (0, 1, 3, 4, 18, 19, 20, 21):
        args[index] = np.repeat(args[index], 129, axis=0)
    expected = native(*args)
    actual = native(*args, cache)
    assert actual[0].shape == (2, 129)
    assert CTRuc._CTRcalc_cpp._dwba_motif_work()["finite_cache_bytes"] == 0
    assert_array_equal(actual, expected)


def test_automatic_scope_restores_outer_native_reuse():
    """A facade request cannot disable reuse in its caller's native batch."""
    model = _model(policy="state_specific")
    state = model.coherent_model.dwba
    with state.batch():
        model.dwba.evaluate(*_hkl())
        assert state._batch_reuse
        assert state._batch_native is None


def test_copying_crystal_does_not_copy_native_scratch():
    """Opaque native reuse cannot make an otherwise copyable crystal fail."""
    crystal = _model().coherent_model
    with crystal.dwba.batch():
        crystal.dwba.evaluate(*_hkl())
        assert crystal.dwba._batch_native is not None
        duplicate = copy.deepcopy(crystal)
        assert duplicate.dwba._batch_native is None
        assert duplicate.dwba._crystal is duplicate


def test_shared_basis_keeps_its_existing_reuse_path():
    """Mean/basis alternation does not build ineffective cross-channel caches."""
    model = _model(policy="shared")
    state = model.coherent_model.dwba
    with model.dwba.batch():
        model.dwba.polarization_reflectivity(*_hkl(), outgoing="unanalysed")
        assert "states" not in model.dwba._scope
        assert state._batch_profiles is None
        assert state._batch_native is None
        assert not state._batch_field_blocks
        assert state.cache_info()["profile_hits"] == 0
    assert state._batch_profiles == {} and state._batch_reuse


def test_zero_channel_is_not_reused_as_a_computed_sum(native_arguments):
    """A formerly absent optical channel must compute its finite contribution."""
    native, args = native_arguments
    cache = CTRuc._CTRcalc_cpp._DWBAFiniteCache()
    original = args[7]
    args[7] = np.zeros_like(original)
    native(*args, cache)
    args[7] = original
    actual = native(*args, cache)
    assert CTRuc._CTRcalc_cpp._dwba_motif_work()["finite_cache_hits"] == 0
    assert_allclose(actual, native(*args), rtol=1e-13, atol=1e-28)


@pytest.mark.parametrize("factor, coefficient, depth", [
    (1e200, 1e-150, -230.), (1e-200, 1e150, 230.)])
def test_extreme_range_falls_back_to_scaled_sum(
    native_arguments, factor, coefficient, depth
):
    """Separating optical and atomic factors cannot overflow or underflow."""
    native, args = native_arguments
    args[1][:] = 0
    args[5][:] = args[10][:] = -1 + 1j
    args[6][:] = args[11][:] = 0
    args[7][:] = args[12][:] = coefficient
    args[29][:, 4:6] = 0
    args[30][:] = 0
    args[30][:, 10] = factor
    args[31][:, 2] = depth
    for index in (35, 36, 37, 38):
        args[index] = args[index][:0]
    expected = native(*args)
    assert np.all(np.isfinite(expected))
    cache = CTRuc._CTRcalc_cpp._DWBAFiniteCache()
    actual = native(*args, cache)
    assert_allclose(actual, expected, rtol=1e-13, atol=0.)
    assert CTRuc._CTRcalc_cpp._dwba_motif_work()["finite_cache_bytes"] == 0


@pytest.mark.parametrize("change", ["occupancy", "factors", "energy", "origin",
                                   "domains", "planes", "cycle"])
def test_profile_reuse_uses_complete_cell_state(change):
    """Equal copies share optics; density, planes, energy and domains do not."""
    state = _model().coherent_model.dwba
    cell = state._crystal.uc_bulk
    with state.batch():
        first = state._cell_profile(cell)
        other = copy.deepcopy(cell)
        assert state._cell_profile(other) is first
        if change == "occupancy":
            other.basis[0, 6] *= .9
        elif change == "factors":
            other.f[0, 11] += .1
        elif change == "energy":
            other.setEnergy(19.)
        elif change == "origin":
            other.layerpos[0.] += .01
        elif change == "domains":
            other.coherentDomainOccupancy[0] *= .9
        elif change == "planes":
            other._optical_plane_geometry = {0.: (0., .4), 1.: (.4, .6)}
        else:
            other._explicit_layer_cycle = tuple(other.layers[::-1])
        actual = state._cell_profile(other)
        assert actual is not first
        assert_array_equal(actual, other.optical_profile())
    assert not state._batch_profiles


@pytest.mark.parametrize("abort", [False, True])
def test_ensemble_shares_profiles_and_releases_batch_storage(abort):
    """Snapshots reuse cells, and nested scopes release native/field storage."""
    model = _model(policy="state_specific")
    retained = []
    try:
        with model.dwba.batch():
            first = model.dwba.evaluate(*_hkl())
            retained = [model.coherent_model.dwba] + [
                state for state, _, _ in model.dwba._scope["states"]]
            assert len({id(state._batch_profiles) for state in retained}) == 1
            assert sum(s.cache_info()["profile_hits"] for s in retained) > 0
            with model.dwba.batch():
                again = model.dwba.evaluate(*_hkl())
                assert_array_equal(first.reflectivity, again.reflectivity)
            assert all(s._batch_native is not None for s in retained)
            assert all(s._batch_field_blocks for s in retained)
            if abort:
                raise RuntimeError("test abort")
    except RuntimeError as error:
        assert abort and str(error) == "test abort"
    assert all(not s._batch_profiles and not s._batch_field_blocks
               and s._batch_native is None for s in retained)

"""Shared terrace grids preserve optics and shares at floating-point limits."""

from dataclasses import replace

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from .. import _CTRincoherent_dwba as ensemble
from ..CTRoptics import StratifiedProfile, _validated_boundaries
from .test_CTRincoherent_dwba import _hkl, _model


@pytest.mark.parametrize("direction", [-np.inf, None, np.inf])
def test_shared_preparation_coincident_boundaries_and_reuse(direction):
    """One-ULP state offsets retain common interfaces and both share arrays."""
    model = _model()
    common = model.coherent_model.dwba.prepare(*_hkl())
    flat = ensemble._flat_crystal(model, 0)
    _, atomic, profile = flat.dwba._atomic_model_snapshot()
    boundaries = profile.boundaries.copy()
    if direction is not None:
        boundaries = np.nextafter(boundaries, direction)
    shifted = StratifiedProfile(profile.values, boundaries)
    # Both source profiles are valid; only their naive union is problematic.
    _validated_boundaries(shifted.values, shifted.boundaries)
    original_shares = atomic.reference_shares.copy()
    original_layers = atomic.layer_reference_shares.copy()
    expected = ensemble._shared_preparation(
        flat.dwba, common, atomic, profile)
    actual = ensemble._shared_preparation(flat.dwba, common, atomic, shifted)
    assert_array_equal(actual.reference.profile.boundaries,
                       common.reference.profile.boundaries)
    assert_array_equal(actual.reference.profile.values[:, 1:],
                       common.reference.profile.values[:, 1:])
    assert_array_equal(actual._atomic_model.reference_shares, original_shares)
    assert_array_equal(actual._atomic_model.layer_reference_shares,
                       original_layers)
    _validated_boundaries(actual.reference.profile.values,
                          actual.reference.profile.boundaries)
    result = flat.dwba.evaluate_prepared(actual)
    assert_allclose(result.total_amplitude,
                    flat.dwba.evaluate_prepared(expected).total_amplitude,
                    rtol=2e-13, atol=1e-18)
    again = ensemble._shared_preparation(flat.dwba, common, atomic, shifted)
    assert again.reference is actual.reference
    assert_array_equal(flat.dwba.evaluate_prepared(again).total_amplitude,
                       result.total_amplitude)
    assert_array_equal(atomic.reference_shares, original_shares)
    assert_array_equal(atomic.layer_reference_shares, original_layers)
    assert_array_equal(shifted.boundaries, boundaries)


@pytest.mark.parametrize("direction", [-np.inf, None, np.inf])
def test_shared_coincident_public_endpoints_and_batch(monkeypatch, direction):
    """Roundoff in state interfaces preserves endpoint mixing and batch reuse."""
    model = _model(fraction=1)
    # Exercise the reference-density term as well as atomic scattering.
    model.coherent_model.uc_bulk.crystallinity = .4
    expected = model.dwba.evaluate(*_hkl())
    original = ensemble.DWBAState._atomic_model_snapshot
    shifted_states = []

    def _shifted_snapshot(state, *args, **kwargs):
        result = original(state, *args, **kwargs)
        if state._crystal is model.coherent_model:
            return result
        records, atomic, profile = result
        boundaries = profile.boundaries.copy()
        if direction is not None:
            boundaries = np.nextafter(boundaries, direction)
        shifted_states.append(state)
        return records, atomic, replace(profile, boundaries=boundaries)

    monkeypatch.setattr(ensemble.DWBAState, "_atomic_model_snapshot",
                        _shifted_snapshot)
    coherent = model.coherent_model.dwba.evaluate(*_hkl()).reflectivity
    for fraction in (0., .20542522935100496, 1.):
        model.incoherent_fraction = fraction
        with model.dwba.batch():
            actual = model.dwba.evaluate(*_hkl(), bulk_attenuation=0.)
            again = model.dwba.evaluate(*_hkl(), bulk_attenuation=0.)
        assert_allclose(actual.reflectivity,
                        (1 - fraction) * coherent
                        + fraction * expected.incoherent_intensity,
                        rtol=2e-12, atol=1e-18)
        assert_array_equal(actual.reflectivity, again.reflectivity)
    assert shifted_states


def test_shared_resolvable_slivers_keep_optics_and_state_shares():
    """Even two-ULP optical slabs and independent state splits are retained."""
    model = _model()
    common = model.coherent_model.dwba.prepare(*_hkl())
    flat = ensemble._flat_crystal(model, 0)
    _, atomic, profile = flat.dwba._atomic_model_snapshot()
    # Use the reported positive height, with a real narrow optical layer.
    z = 70.99394978244919
    middle = np.nextafter(z, np.inf)
    top = np.nextafter(middle, np.inf)
    shared = StratifiedProfile(
        np.array([[z - 1, 1e-5, 1e-7], [middle, 2e-5, 2e-7],
                  [top + 1, 0., 0.]]), np.array([z, top]))
    common = replace(common, reference=flat.dwba._reference(shared))
    state = StratifiedProfile(
        np.array([[z - 2, 1e-5, 1e-7], [z - .5, 3e-5, 3e-7],
                  [top + 1, 0., 0.]]), np.array([z - 1, middle]))
    shares = np.zeros((len(atomic.reference_shares), 3, 2))
    layers = np.zeros((len(atomic.layer_reference_shares), 3, 2))
    shares[0] = state.values[::-1, 1:]
    layers[0] = state.values[::-1, 1:]
    atomic = replace(atomic, reference_shares=shares,
                     layer_reference_shares=layers)
    prepared = ensemble._shared_preparation(flat.dwba, common, atomic, state)
    result = prepared.reference.profile
    assert_array_equal(result.boundaries, [z - 1, z, top])
    assert_array_equal(result.values[:, 1:], shared.values[[0, 0, 1, 2], 1:])
    _validated_boundaries(result.values, result.boundaries)
    # The independently resolved state split survives in the substrate medium.
    assert_array_equal(prepared._atomic_model.reference_shares[:, -2, :],
                       shares[:, 1, :])
    assert_array_equal(prepared._atomic_model.layer_reference_shares[:, -2, :],
                       layers[:, 1, :])


def test_shared_grid_retains_resolvable_state_slivers():
    """Float-scale state intervals with an interior float are not coalesced."""
    z = 70.99394978244919
    ulp = np.spacing(z)
    shared = StratifiedProfile(
        np.array([[z - 1, 1e-5, 1e-7], [z + .5, 2e-5, 2e-7],
                  [z + 2, 0., 0.]]), np.array([z, z + 1]))
    # The middle state boundary is unresolvable with its left neighbor; the
    # two-ULP interval on either side must survive, even within one medium.
    boundaries = z + np.array([2., 3., 4.]) * ulp
    result = ensemble._shared_grid(shared, boundaries)
    assert_array_equal(result.boundaries, [z, z + 2 * ulp, z + 4 * ulp, z + 1])
    assert_array_equal(result.values[:, 1:], shared.values[[0, 1, 1, 1, 2], 1:])
    _validated_boundaries(result.values, result.boundaries)

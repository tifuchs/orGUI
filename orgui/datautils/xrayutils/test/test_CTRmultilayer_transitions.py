"""Explicit starting-layer transitions across renamed strict stack cycles."""

import copy
import pickle

import numpy as np
import pytest

from .. import CTRcalc, CTRfilm, CTRuc
from ..CTRstacking import LayerTransition, StackValidationError
from .test_CTRmultilayer import _parts, _hkl


def _rename(cell, ids):
    cell = copy.deepcopy(cell)
    old_cycle = cell.layer_cycle.layers
    mapping = dict(zip(old_cycle, ids))
    for basis in (cell.basis, cell.basis_0):
        original = basis[:, 7].copy()
        for old, new in mapping.items():
            basis[original == old, 7] = new
    cell.layerpos = {mapping[key]: cell.layerpos[key] for key in old_cycle}
    cell.layer_cycle = tuple(mapping[key] for key in old_cycle)
    return cell


def _next_map(lower_ids, upper_ids):
    return {lower: upper_ids[(i + 1) % len(upper_ids)]
            for i, lower in enumerate(lower_ids)}


def _models(*, widths=(40, 41), positions=(0., .5),
            shift=(0., 0., 1.), colliding=False, skew=.1, coupling=.4,
            select_first=False):
    bulk, parts = _parts(widths=widths, positions=positions, shift=shift,
                         repeats=(6., 5.2, 7.1, 6.))
    for interface in parts[::2]:
        interface.basis[1:] = [skew, coupling, .15]
        interface.basis_0[:] = interface.basis
    count = len(positions)
    original = tuple(range(1, count + 1))
    if colliding:
        ids = [tuple(np.roll(original, i)) for i in range(5)]
    else:
        ids = [tuple(100 * (i + 1) + j * 7 for j in range(count)) for i in range(5)]
    renamed_bulk = _rename(bulk, ids[0])
    renamed = []
    def transition(lower, upper):
        if select_first:
            return {layer: upper[0] for layer in lower}
        return _next_map(lower, upper)

    for i, (interface, film) in enumerate(zip(parts[::2], parts[1::2])):
        interface_ids, film_ids = ids[2 * i + 1:2 * i + 3]
        new_interface = CTRfilm.EpitaxyInterface(
            _rename(interface.uc_top, interface_ids),
            _rename(interface.uc_bottom, interface_ids),
            name=interface.name, profile=interface.profile,
            layer_transition=transition(ids[2 * i], interface_ids),
        )
        new_interface.basis[:] = new_interface.basis_0[:] = interface.basis
        new_film = CTRfilm.Film(
            _rename(film.unitcell, film_ids), name=film.name,
            width_layers=film.basis[0],
            layer_transition=transition(interface_ids, film_ids),
        )
        renamed.extend((new_interface, new_film))
    if select_first:
        for component in parts:
            component.layer_transition = LayerTransition(transition(original, original))
    baseline = CTRcalc.SXRDCrystal.from_stack(bulk, parts)
    actual = CTRcalc.SXRDCrystal.from_stack(renamed_bulk, renamed)
    for crystal in (baseline, actual):
        crystal.setEnergy(18000.)
    return baseline, actual


@pytest.mark.parametrize("colliding", [False, True])
@pytest.mark.parametrize("widths", [(40, 41), (41, 40)])
@pytest.mark.parametrize("positions,shift", [
    ((0., .5), (0., 0., 1.)),
    ((.1, .3, .8), (.5, .25, 1.)),
])
def test_relabelled_cycles_preserve_geometry_and_observables(
    colliding, widths, positions, shift,
):
    """Label namespaces cannot move the physical lattice or change its phase."""
    baseline, actual = _models(colliding=colliding, widths=widths,
                               positions=positions, shift=shift)
    np.testing.assert_allclose(actual.F(*_hkl()), baseline.F(*_hkl()),
                               rtol=1e-12, atol=1e-10)
    for left, right in zip(baseline.uc_surface_list, actual.uc_surface_list):
        np.testing.assert_allclose(
            [left.pos_absolute, left.height_absolute, left.stacking_loc_absolute],
            [right.pos_absolute, right.height_absolute, right.stacking_loc_absolute],
            rtol=0, atol=1e-11)
    z = np.linspace(-50, 350, 501)
    np.testing.assert_allclose(actual.zDensity_G(z, 0., 0.),
                               baseline.zDensity_G(z, 0., 0.), atol=1e-11)
    np.testing.assert_allclose(actual.optical_profile(), baseline.optical_profile(),
                               rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("boundary", [0, 1, 2, 3])
@pytest.mark.parametrize("invalid", ["missing", "source", "target", "geometry"])
def test_invalid_transition_or_geometry_rejected(boundary, invalid):
    """Explicit mappings do not bypass geometry checks at either component type."""
    _, crystal = _models()
    component = crystal.uc_surface_list[boundary]
    original = copy.deepcopy(component.layer_transition)
    if invalid == "missing":
        component.layer_transition = LayerTransition({})
    elif invalid == "source":
        component.layer_transition.mapping[-999] = next(iter(original.mapping.values()))
    elif invalid == "target":
        key = next(iter(original.mapping))
        component.layer_transition.mapping[key] = -999
    else:
        cell = (component.uc_bottom if isinstance(component, CTRfilm.EpitaxyInterface)
                else component.unitcell)
        key = next(iter(cell.layerpos))
        cell.layerpos[key] += .02
    with pytest.raises(StackValidationError):
        crystal.F(*_hkl())


def test_no_transition_does_not_infer_relabelling():
    """Renamed labels require explicit starting-phase selection."""
    _, crystal = _models()
    crystal["interface_0"].layer_transition = None
    with pytest.raises(StackValidationError, match="layer geometry"):
        crystal.apply_stacking()


@pytest.mark.parametrize("method", ["deepcopy", "pickle", "text"])
def test_transition_persistence_and_parameter_history(method):
    """Derived registration maps rebuild on copies, load and parameter changes."""
    _, crystal = _models(colliding=True)
    expected = crystal.F(*_hkl())
    if method == "deepcopy":
        restored = copy.deepcopy(crystal)
    elif method == "pickle":
        restored = pickle.loads(pickle.dumps(crystal))
    else:
        restored = CTRcalc.SXRDCrystal.fromStr(crystal.toStr(False))
    np.testing.assert_allclose(restored.F(*_hkl()), expected, atol=1e-10)
    for skew in (-.1, .2, .1):
        restored["interface_1"].basis[1] = skew
        fresh = _models(colliding=True)[1]
        fresh["interface_1"].basis[1] = skew
        np.testing.assert_allclose(restored.F(*_hkl()), fresh.F(*_hkl()), atol=1e-10)


@pytest.mark.parametrize("colliding", [False, True])
def test_dwba_preserves_relabelled_stack_amplitudes(colliding):
    """Signed reference records preserve total DWBA despite different labels."""
    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA requires the native CTR extension")
    baseline, actual = _models(colliding=colliding)
    results = []
    for crystal in (baseline, actual):
        prepared = crystal.dwba.prepare_from_glancing(
            np.array([0., .6]), 0., np.array([.04, .04]), np.array([.04, .05]))
        results.append(crystal.dwba.evaluate_prepared(prepared))
    for attribute in ("F_atomic", "F_reference", "F_h"):
        np.testing.assert_allclose(getattr(results[0], attribute),
                                   getattr(results[1], attribute),
                                   rtol=1e-11, atol=1e-10)


@pytest.mark.parametrize("coupling", [0., 1.])
@pytest.mark.parametrize("colliding", [False, True])
def test_many_to_one_start_selection_is_not_a_label_correspondence(
    coupling, colliding,
):
    """Different terminations may deliberately select the same starting layer."""
    baseline, actual = _models(coupling=coupling, colliding=colliding,
                               select_first=True, widths=(41, 40))
    np.testing.assert_allclose(actual.F(*_hkl()), baseline.F(*_hkl()),
                               rtol=1e-12, atol=1e-10)


def test_partial_map_revalidated_when_film_termination_changes():
    """A partial map is valid only while it covers the live lower termination."""
    _, crystal = _models(widths=(40, 40))
    interface = crystal["interface_1"]
    terminal = crystal["film_0"].layer_state.layer
    original = copy.deepcopy(interface.layer_transition)
    interface.layer_transition = LayerTransition(
        {terminal: original.mapping[terminal]})
    expected = crystal.F(*_hkl())
    crystal["film_0"].basis[0] = 41
    with pytest.raises(StackValidationError, match="No layer transition"):
        crystal.F(*_hkl())
    interface.layer_transition = original
    assert np.isfinite(crystal.F(*_hkl())).all()
    crystal["film_0"].basis[0] = 40
    np.testing.assert_array_equal(crystal.F(*_hkl()), expected)


def test_transition_cannot_join_different_layer_counts():
    """An explicit selector does not define a two-to-three-layer interface."""
    from .test_CTRmultilayer import _cell

    bulk, parts = _parts()
    bulk = _cell("bulk", positions=(0., 1 / 3, 2 / 3))
    parts[0].layer_transition = LayerTransition({1: 1, 2: 1, 3: 1})
    with pytest.raises(StackValidationError, match="unequal layer counts"):
        CTRcalc.SXRDCrystal.from_stack(bulk, parts)


def test_inactive_transition_cannot_bypass_geometry_validation():
    """Ignore mode must not use a transition that stacking will not execute."""
    _, crystal = _models()
    crystal["interface_0"].layer_behavior = "ignore"
    with pytest.raises(StackValidationError, match="layer geometry"):
        crystal.apply_stacking()

"""Regressions for fixed lower-material registration in epitaxy interfaces."""

import copy

import numpy as np
import pytest
from scipy.optimize import brentq

from ..CTRdistributions import SkellamProfile
from .test_CTRinterface_anchors import _H, _K, _L, _cell, _evaluate, _sites, _stack


@pytest.mark.parametrize("coupling", [0.0, 0.4, 1.0])
@pytest.mark.parametrize("threshold", [-2 / 3, 2 / 3])
def test_three_layer_registration_threshold_is_continuous(coupling, threshold):
    """Changing the rounded Skellam mean cannot translate the lower crystal."""
    crystal = _stack(
        width=32.25,
        positions=(0.0, 1 / 3, 2 / 3),
        ids=(1, 2, 3),
        interface_width=0.5,
        coupling=coupling,
    )
    center = _evaluate(crystal, threshold)
    counts = []
    for epsilon in (1e-4, 1e-6, 1e-8):
        left = _evaluate(crystal, threshold - epsilon)
        counts.append(crystal["film"]._layers_to_create)
        right = _evaluate(crystal, threshold + epsilon)
        counts.append(crystal["film"]._layers_to_create)
        for value in (left, right):
            for actual, expected in zip(value, center):
                np.testing.assert_allclose(
                    actual, expected, rtol=0, atol=2e3 * epsilon + 1e-9
                )
    assert len(set(counts)) == 1


def _reference_atoms(crystal):
    """Enumerate a finite lower stack independently of interface domains."""
    bulk = crystal.uc_bulk
    sites = [
        (atom[3] * bulk.a[2] - repeat * bulk.a[2], atom[6])
        for repeat in range(100)
        for atom in bulk.basis
    ]
    registry = next(
        (uc for uc in crystal.uc_surface_list if uc.name == "registry"), None
    )
    if registry is not None:
        ids = list(registry.layer_order)
        positions = registry.unitcell.layerpos
        intervals = [
            (positions[ids[(i + 1) % len(ids)]] - positions[layer]) % 1
            for i, layer in enumerate(ids)
        ]
        z = registry.below_H
        for index in range(int(registry.basis[0])):
            sites.append((z, 1.0))
            z += intervals[index % len(ids)] * registry.unitcell.a[2]
    return sites


def _oracle(crystal, skew):
    """Enumerate replacement atoms using cumulative physical intervals only."""
    interface = crystal["interface"]
    top, bottom = interface.uc_top, interface.uc_bottom
    ids = list(interface.layer_order)
    n = len(ids)
    profile = SkellamProfile(
        interface.basis[0], skew, interface.profile.tail_probability
    )
    low, high = profile.support(n)
    indices = np.arange(min((low // n) * n, 0), max((high // n + 1) * n, n))
    p = profile.occupancy(indices, n)
    intervals = [
        ((top.layerpos[ids[(i + 1) % n]] - top.layerpos[layer]) % 1)
        for i, layer in enumerate(ids)
    ]
    distances = np.tile(intervals, len(indices) // n)
    ideal = np.r_[0.0, np.cumsum(distances)]
    mean = interface.basis[0] ** 2 * n**2 * skew
    upper_anchor = np.interp(mean - indices[0], np.arange(len(ideal)), ideal)
    # Repeat zero belongs to F_bulk. Finite Film uses its actual top boundary.
    if len(crystal.uc_surface_list) == 3:
        bulk = crystal.uc_bulk
        origin = bulk.layerpos[ids[0]]
        origin += np.floor(max(bulk.layerpos[i] for i in bulk.layers) - origin) + 1
        origin *= bulk.a[2]
    else:
        registry = crystal["registry"]
        order = list(registry.layer_order)
        source = registry.unitcell
        origin = registry.below_H
        # Walk the source lattice forward to the first requested boundary
        # at/above the end of the finite registry Film.
        for index in range(int(registry.basis[0]) + len(order)):
            layer = order[index % len(order)]
            if index >= registry.basis[0] and layer == ids[0]:
                break
            next_layer = order[(index + 1) % len(order)]
            origin += (
                (source.layerpos[next_layer] - source.layerpos[layer]) % 1 * source.a[2]
            )
    lower_start = origin + indices[0] / n * bottom.a[2]
    heights = distances * (p * top.a[2] + (1 - p) * bottom.a[2])
    coupled = lower_start + np.r_[0.0, np.cumsum(heights)]
    k = interface.basis[2]
    offset = interface.basis[3] * bottom.a[2]
    offset_profile = offset * (p - p.min()) / (p.max() - p.min())
    upper_z = (
        (1 - k) * (interface.below_H + (ideal[:-1] - upper_anchor) * top.a[2])
        + k * coupled[:-1]
        + (1 - k) * offset
        + k * offset_profile
    )
    lower_z = (
        (1 - k) * (lower_start + ideal[:-1] * bottom.a[2])
        + k * coupled[:-1]
        + k * offset_profile
    )
    removed = [
        (z, -occupancy)
        for z, occupancy in _reference_atoms(crystal)
        if lower_start - 1e-8 <= z < lower_start + ideal[-1] * bottom.a[2] - 1e-8
    ]
    return (list(zip(upper_z, p)), list(zip(lower_z, 1 - p)) + removed)


def _atom_amplitude(crystal, element, sites):
    atom = _cell("oracle", element=element, positions=(0.0,), ids=(1,))
    atom.setReferenceUnitCell(crystal.uc_bulk)
    atom.setEnergy(18000.0)
    single = atom.F_uc(_H, _K, _L)
    return sum(
        occupancy * single * np.exp(2j * np.pi * _L * z / crystal.uc_bulk.a[2])
        for z, occupancy in sites
    )


@pytest.mark.parametrize("registry", [None, 2, 3])
@pytest.mark.parametrize("coupling", [0.0, 0.4, 1.0])
@pytest.mark.parametrize(
    "positions,ids,transition",
    [
        ((0.0, 0.5), (1, 2), None),
        ((0.1, 0.3, 0.8), (7, 42, 13), {7: 7, 42: 42, 13: 13}),
    ],
)
def test_signed_atom_and_amplitude_oracle(
    registry, coupling, positions, ids, transition
):
    """Remove exactly the supplied finite atoms, including registry/bulk overlap."""
    crystal = _stack(
        width=64.25,
        registry=registry,
        coupling=coupling,
        positions=positions,
        ids=ids,
        transition=transition,
        interface_width=0.5,
        mismatch=0.85,
        offset=0.2,
    )
    for skew in (-0.7, 0.0, 0.7):
        _evaluate(crystal, skew)
        top, bottom = _oracle(crystal, skew)
        interface = crystal["interface"]
        # Compare signed atom sums via complex phase at several irrational
        # spatial frequencies; zero-weight domain bookkeeping is irrelevant.
        for actual, expected in (
            (_sites(interface.top_layers), top),
            (_sites(interface.bottom_layers + interface._reference_layers), bottom),
        ):
            for sign in (-1, 1):
                actual_sites = sorted((z, w) for z, w in actual if sign * w > 1e-13)
                expected_sites = sorted((z, w) for z, w in expected if sign * w > 1e-13)
                np.testing.assert_allclose(
                    actual_sites, expected_sites, rtol=0, atol=1e-10
                )
            for q in (0.0, np.sqrt(2), np.sqrt(3), np.pi):
                np.testing.assert_allclose(
                    sum(w * np.exp(1j * q * z) for z, w in actual),
                    sum(w * np.exp(1j * q * z) for z, w in expected),
                    rtol=0,
                    atol=1e-10,
                )
        expected = _atom_amplitude(crystal, "Si", top) + _atom_amplitude(
            crystal, "C", bottom
        )
        np.testing.assert_allclose(
            interface.F_uc(_H, _K, _L), expected, rtol=0, atol=1e-9
        )


@pytest.mark.parametrize("n", [2, 3, 4])
def test_all_reachable_rounding_thresholds_and_history(n):
    """Positive/negative ties, rotations, copies, and reverse updates agree."""
    width = 1.0
    thresholds = (np.arange(-n, n) + 0.5) / n
    for start in range(n):
        ids = tuple(range(1, n + 1))
        crystal = _stack(
            width=100.25,
            interface_width=width,
            positions=tuple(np.arange(n) / n),
            ids=ids,
            transition={i: ids[start] for i in ids},
        )
        for skew in thresholds:
            center = _evaluate(crystal, skew)
            clone = copy.deepcopy(crystal)
            for candidate in (crystal, clone):
                for value in (skew + 1e-8, skew - 1e-8, skew):
                    actual = _evaluate(candidate, value)
                    for a, b in zip(actual, center):
                        np.testing.assert_allclose(
                            a, b, rtol=0, atol=1e-3 if value != skew else 1e-9
                        )


@pytest.mark.parametrize("coupling", [0.0, 0.4, 1.0])
@pytest.mark.parametrize("start", [1, 2, 3])
@pytest.mark.parametrize("registry", [None, 2, 3])
def test_support_extension_with_mismatch_and_unequal_intervals(
    coupling, start, registry
):
    """Tail extension keeps occupied lower sites fixed and only adds tail error."""
    crystals = [
        _stack(
            width=100.25,
            registry=registry,
            interface_width=0.5,
            coupling=coupling,
            positions=(0.1, 0.3, 0.8),
            ids=(1, 2, 3),
            transition={1: start, 2: start, 3: start},
            mismatch=0.85,
            offset=0.2,
            tail=tail,
        )
        for tail in (1e-7, 1e-12)
    ]
    if registry is not None:
        # Film places its first native .1 boundary at absolute zero. Put the
        # bulk on that same lattice for a meaningful tail-convergence test.
        # A physically misregistered lower stack is a different model.
        for crystal in crystals:
            bulk = crystal.uc_bulk
            bulk.basis[:, 3] -= 0.1
            bulk.basis_0[:, 3] -= 0.1
            bulk.layerpos = {key: value - 0.1 for key, value in bulk.layerpos.items()}
    results = [_evaluate(crystal, 0.35) for crystal in crystals]
    assert (
        crystals[0]["film"]._layers_to_create != crystals[1]["film"]._layers_to_create
    )
    # The coarse support loses < 1e-7 occupancy per omitted layer. Allow a
    # conservative 100-layer budget times a repeat for accumulated strain,
    # and 20 electrons per atom for amplitude plus displacement sensitivity.
    np.testing.assert_allclose(results[0][0][4:], results[1][0][4:], rtol=0, atol=1e-5)
    np.testing.assert_allclose(results[0][1], results[1][1], rtol=0, atol=2e-4)
    np.testing.assert_allclose(results[0][2], results[1][2], rtol=0, atol=2e-4)


def test_energy_refresh_and_fresh_construction_preserve_reference():
    """Reference copies refresh after source scattering factors change."""
    crystal = _stack(reverse=True, registry=3, interface_width=0.5)
    _evaluate(crystal, 0.7)
    crystal.setEnergy(12000.0)
    actual = crystal["interface"].F_uc(_H, _K, _L)
    fresh = _stack(reverse=True, registry=3, interface_width=0.5)
    _evaluate(fresh, 0.7)
    fresh.setEnergy(12000.0)
    fresh.apply_stacking()
    np.testing.assert_allclose(
        actual, fresh["interface"].F_uc(_H, _K, _L), rtol=0, atol=1e-9
    )


@pytest.mark.parametrize("coupling", [0.0, 0.4, 1.0])
def test_half_occupancy_does_not_change_subtraction(coupling):
    """Crossing CDF=0.5 cannot remove or restore an underlying atom."""
    threshold = brentq(
        lambda skew: SkellamProfile(0.5, skew).occupancy(np.array([0]), 3)[0] - 0.5,
        -0.99,
        0.99,
    )
    crystal = _stack(
        width=64.25,
        interface_width=0.5,
        coupling=coupling,
        positions=(0.0, 1 / 3, 2 / 3),
        ids=(1, 2, 3),
    )
    left, right = (_evaluate(crystal, threshold + sign * 1e-8) for sign in (-1, 1))
    for a, b in zip(left, right):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-4)


@pytest.mark.parametrize("n", [2, 3])
def test_fixed_support_uses_complete_cycles_and_fixed_registration(n):
    """A rounded support-window shift changes tails, not the lattice phase."""
    ids = tuple(range(1, n + 1))
    crystal = _stack(
        width=100.25, interface_width=0.5, positions=tuple(np.arange(n) / n), ids=ids
    )
    # An odd requested half-width must also include complete cycles.
    crystal["interface"].fixed_ucs = 23
    # A wider distribution puts the crossing inside the allowed skew range.
    crystal["interface"].basis[0] = 1.0
    crystal["interface"].basis_0[0] = 1.0
    threshold = 0.5 / n
    left, right = (_evaluate(crystal, threshold + sign * 1e-8) for sign in (-1, 1))
    for a, b in zip(left[1:], right[1:]):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-3)
    np.testing.assert_allclose(left[0][4:], right[0][4:], rtol=0, atol=1e-7)


@pytest.mark.parametrize("skew", [-1.0, 1.0])
def test_support_spans_sharp_boundary_at_large_skew(skew):
    """Fully occupied material between disjoint sharp/statistical boundaries stays."""
    crystal = _stack(width=300.25, interface_width=3.0)
    _evaluate(crystal, skew)
    top, bottom = _oracle(crystal, skew)
    interface = crystal["interface"]
    low, high = interface.profile.support(2)
    assert low > 0 or high < 0  # The original support misses the sharp boundary.
    expected = _atom_amplitude(crystal, "Si", top) + _atom_amplitude(
        crystal, "C", bottom
    )
    np.testing.assert_allclose(interface.F_uc(_H, _K, _L), expected, rtol=0, atol=1e-9)


def test_explicit_transition_between_differently_named_cycles():
    """Renaming only the bulk cycle does not move its physical lower reference."""
    base = _stack(transition={1: 1, 2: 1}, interface_width=0.5)
    renamed = _stack(transition={42: 1, 7: 1}, interface_width=0.5)
    bulk = renamed.uc_bulk
    bulk.basis[:, 7] = [42, 7]
    bulk.basis_0[:, 7] = [42, 7]
    bulk.layerpos = {42: 0.0, 7: 0.5}
    for skew in (-0.7, 0.0, 0.7):
        for a, b in zip(_evaluate(base, skew), _evaluate(renamed, skew)):
            np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)

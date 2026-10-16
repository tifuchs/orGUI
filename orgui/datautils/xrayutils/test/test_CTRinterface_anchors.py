"""Physical boundary regressions for cyclic epitaxy upper anchors.

The independent references use cumulative layer intervals and explicit atom
sites. Lower registration is covered in test_CTRinterface_registration.py.
"""

import copy
from unittest import mock

import h5py
import numpy as np
import pytest
from silx.io import dictdump

from .. import CTRcalc, CTRfilm
from ..CTRdistributions import PoissonProfile, SkellamProfile


_H = np.repeat([0.0, 1.0, 2.0], 10)
_K = np.zeros(30)
_L = np.tile(np.linspace(0.13, 2.87, 10), 3)
_Z = np.linspace(-30.0, 140.0, 1001)


def _cell(name, element="C", repeat=6.58069633,
          positions=(0.0, 0.5), ids=(1, 2)):
    uc = CTRcalc.UnitCell([3.0, 3.0, repeat], [90.0] * 3, name=name)
    for layer, z in zip(ids, positions):
        uc.addAtom(element, [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=layer)
        uc.layerpos[float(layer)] = z
    return uc


def _stack(width=32, reverse=False, coupling=1.0, offset=0.0,
           mismatch=1.0, registry=None, termination="poisson",
           positions=(0.0, 0.5), ids=(1, 2), repeat=6.58069633,
           tail=1e-10, interface_width=0.15, transition=None):
    bulk = _cell("bulk", repeat=repeat, positions=positions, ids=ids)
    if reverse:
        bulk = bulk.translate_layered([0, 0, -1])
    top = _cell("top", "Si", repeat * mismatch, positions, ids)
    interface = CTRfilm.EpitaxyInterface(
        top, _cell("bottom", repeat=repeat, positions=positions, ids=ids),
        profile=SkellamProfile(interface_width, tail_probability=tail),
        name="interface", layer_transition=transition,
    )
    interface.basis[:] = [interface_width, 0.0, coupling, offset]
    interface.basis_0[:] = interface.basis
    interface.addFitParameter("S", limits=(-1.0, 1.0))
    film = CTRfilm.Film(copy.deepcopy(top), name="film")
    film.basis[:] = film.basis_0[:] = width
    if termination == "poisson":
        surface = CTRfilm.PoissonSurface(
            copy.deepcopy(top), profile=PoissonProfile(0.0), name="surface",
        )
    else:
        surface = copy.deepcopy(top)
        surface.name = "surface"
    components = [interface, film, surface]
    if registry is not None:
        lower = CTRfilm.Film(copy.deepcopy(bulk), name="registry")
        lower.basis[:] = lower.basis_0[:] = registry
        components.insert(0, lower)
    crystal = CTRcalc.SXRDCrystal(
        bulk, *components, stacking=np.arange(1, len(components) + 1),
    )
    crystal.setEnergy(18000.0)
    return crystal


def _evaluate(crystal, skew):
    crystal.setParameters([skew])
    crystal.apply_stacking()
    interface, film = crystal["interface"], crystal["film"]
    geometry = np.array([
        interface.pos_absolute, interface.height_absolute,
        interface.stacking_loc_absolute, film.pos_absolute,
        film.height_absolute, crystal["surface"].height_absolute,
    ])
    return (geometry, crystal.F(_H, _K, _L),
            crystal.zDensity_G(_Z, 0.0, 0.0))


def _sites(layers):
    """Return sorted (z in Angstrom, signed occupancy) for generated atoms."""
    sites = []
    for layer in layers:
        for matrix, occupancy in zip(
            layer.coherentDomainMatrix, layer.coherentDomainOccupancy,
        ):
            for atom in layer.basis:
                z = (matrix[2, :3] @ atom[1:4] + matrix[2, 3]) * layer.a[2]
                sites.append((z, occupancy * atom[6]))
    return np.asarray(sorted(sites))


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("width", [17, 18, 32])
def test_zero_skew_keeps_film_layers_and_converges(width, reverse):
    """Shrinking skew must not remove a repeat from a fixed nominal Film."""
    crystal = _stack(width=width, reverse=reverse)
    zero = _evaluate(crystal, 0.0)
    expected_order = [2, 1] if reverse else [1, 2]
    np.testing.assert_array_equal(crystal["interface"].layer_order, expected_order)
    assert crystal["film"]._layers_to_create == width - 6
    initial_differences = None
    for epsilon in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14):
        left = _evaluate(crystal, -epsilon)
        assert crystal["film"]._layers_to_create == width - 6
        right = _evaluate(crystal, epsilon)
        assert crystal["film"]._layers_to_create == width - 6
        # Coupling one and matched lattices leave physical support stationary.
        for current in (left, right):
            np.testing.assert_allclose(
                current[0][[0, 1, 3, 4, 5]], zero[0][[0, 1, 3, 4, 5]],
                rtol=0, atol=1e-9,
            )
        differences = np.array([
            np.max(np.abs(a - b)) for a, b in zip(left, right)
        ])
        if initial_differences is None:
            initial_differences = differences
            # Distinct materials have a real finite-skew occupancy response.
            assert differences[1] > 1e-6
        budget = initial_differences * epsilon / 1e-4 * 1.05 + 1e-9
        assert np.all(differences <= budget)
        for current in (left, right):
            assert np.all([
                np.max(np.abs(a - b)) <= bound
                for a, b, bound in zip(current, zero, budget)
            ])


@pytest.mark.parametrize("positions,ids", [
    ((0.0, 0.5), (1, 2)),
    ((0.0, 0.2, 0.7), (7, 42, 13)),
    ((0.1, 0.3, 0.8), (11, 3, 91)),
])
def test_anchor_matches_cumulative_boundaries_in_every_rotation(positions, ids):
    """Independent physical intervals locate the anchor, regardless of labels."""
    n = len(ids)
    for start in range(n):
        interface = CTRfilm.EpitaxyInterface(
            _cell("top", positions=positions, ids=ids),
            _cell("bottom", positions=positions, ids=ids),
            profile=SkellamProfile(1.0), fixed_ucs=12,
        )
        interface._set_start_layer(ids[start])
        interface.basis[2] = 0.0  # The ideal support starts at -anchor.
        intervals = np.diff(positions, append=positions[0] + 1.0)
        intervals = np.roll(intervals, -start)
        boundaries = np.concatenate(([0.0], np.cumsum(np.tile(intervals, 24 // n))))
        for boundary in range(-n, n + 1):
            for delta in (-1e-8, 0.0, 1e-8):
                interface.basis[1] = (boundary + delta) / n**2
                profile = SkellamProfile(1.0, interface.basis[1])
                mu1, mu2 = profile.parameters(n)
                mean = mu1 - mu2
                support_first = -12 + round(mean / n) * n
                expected = np.interp(
                    mean - support_first, np.arange(25), boundaries,
                ) * interface.uc_top.a[2]
                interface.createInterfaceCells()
                np.testing.assert_allclose(
                    -interface.pos_absolute, expected, rtol=0, atol=1e-11,
                )


@pytest.mark.parametrize("coupling", [0.0, 0.4, 1.0])
@pytest.mark.parametrize("mismatch", [1.0, 0.85])
@pytest.mark.parametrize("offset", [0.0, 0.2])
def test_zero_crossing_with_strain_and_offset(coupling, mismatch, offset):
    """Both independent and coupled lattices remain continuous near zero."""
    crystal = _stack(reverse=True, coupling=coupling,
                     mismatch=mismatch, offset=offset)
    left = _evaluate(crystal, -1e-8)
    left_count = crystal["film"]._layers_to_create
    right = _evaluate(crystal, 1e-8)
    assert crystal["film"]._layers_to_create == left_count
    for a, b in zip(left, right):
        np.testing.assert_allclose(a, b, rtol=0, atol=2e-6)
    assert crystal["film"].pos_absolute == pytest.approx(
        crystal["interface"].height_absolute, abs=1e-9,
    )


@pytest.mark.parametrize("registry", [2, 3])
@pytest.mark.parametrize("termination", ["poisson", "unitcell"])
def test_native_history_copies_and_fixed_inputs(registry, termination):
    """Parameter history cannot alter fixed geometry or width references."""
    crystal = _stack(registry=registry, termination=termination)
    interface, film = crystal["interface"], crystal["film"]
    cells = [crystal.uc_bulk, interface.uc_top, interface.uc_bottom, film.unitcell]
    fixed = [(c.basis.copy(), c.basis_0.copy(), dict(c.layerpos)) for c in cells]
    for sequence in ((-1e-4, 0.0, 1e-4), (1e-4, 0.0, -1e-4)):
        for skew in sequence:
            actual = _evaluate(crystal, skew)
            fresh = _evaluate(_stack(registry=registry, termination=termination), skew)
            copied = _evaluate(copy.deepcopy(crystal), skew)
            for result in (fresh, copied):
                for a, b in zip(actual, result):
                    np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)
    for cell, (basis, basis_0, layerpos) in zip(cells, fixed):
        np.testing.assert_array_equal(cell.basis, basis)
        np.testing.assert_array_equal(cell.basis_0, basis_0)
        assert cell.layerpos == layerpos
    np.testing.assert_array_equal(film.basis, [32])
    np.testing.assert_array_equal(film.basis_0, [32])


@pytest.mark.parametrize("registry", [2, 3])
def test_upper_material_has_no_missing_or_duplicate_sites(registry):
    """Explicit uniform sites account for interface support and remaining Film."""
    crystal = _stack(registry=registry)
    _evaluate(crystal, 0.0)
    interface, film = crystal["interface"], crystal["film"]
    sites = _sites(interface.top_layers + film.layer_ucs)
    spacing = film.unitcell.a[2] / 2
    expected_z = np.arange(
        interface.pos_absolute, film.height_absolute - spacing / 2, spacing,
    )
    np.testing.assert_allclose(sites[:, 0], expected_z, rtol=0, atol=1e-9)
    first, last = interface.profile.support(2)
    support = np.arange(first // 2 * 2, (last // 2 + 1) * 2)
    expected_occupancy = np.concatenate((
        interface.profile.occupancy(support, 2),
        np.ones(film._layers_to_create),
    ))
    np.testing.assert_allclose(sites[:, 1], expected_occupancy, rtol=0, atol=1e-14)


def _legacy_anchor(self, coordinate, domains, layer_space):
    # Compatibility reference only: preserve the original, incorrect formula.
    n = len(self.layer_order)
    cycle, layer = divmod(int(np.floor(coordinate)), n)
    matrix = domains[layer][cycle]
    return self.top_layers[0].a[2] * (
        matrix[2, 3] + (coordinate % n) / n * matrix[2, 2]
    )


@pytest.mark.parametrize("coupling", [0.0, 0.4, 1.0])
@pytest.mark.parametrize("reverse", [False, True])
def test_default_compatibility_and_unchanged_lower_domains(coupling, reverse):
    """Only the inconsistent upper reference changes; lower domains stay fixed."""
    crystal = _stack(reverse=reverse, coupling=coupling,
                     mismatch=0.85, offset=0.2)
    for skew in (-0.125, 0.0, 0.125):
        actual = _evaluate(crystal, skew)
        interface = crystal["interface"]
        lower = [np.copy(layer.coherentDomainMatrix)
                 for layer in interface.bottom_layers]
        sites = _sites(interface.bottom_layers)
        with mock.patch.object(CTRfilm.EpitaxyInterface, "_ideal_top_anchor",
                               _legacy_anchor):
            legacy = _evaluate(crystal, skew)
            np.testing.assert_array_equal(_sites(interface.bottom_layers), sites)
            for layer, matrices in zip(interface.bottom_layers, lower):
                np.testing.assert_array_equal(layer.coherentDomainMatrix, matrices)
        if not reverse:
            for a, b in zip(actual, legacy):
                np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)


def test_consistent_relabelling_and_explicit_transition_preserve_predictions():
    """An explicit transition selects a rotation independently of numeric IDs."""
    base = _stack(transition={1: 1, 2: 2})
    relabelled = _stack(ids=(42, 7), transition={42: 42, 7: 7})
    for skew in (-1e-4, 0.0, 1e-4):
        expected, actual = _evaluate(base, skew), _evaluate(relabelled, skew)
        np.testing.assert_array_equal(base["interface"].layer_order, [2, 1])
        np.testing.assert_array_equal(relabelled["interface"].layer_order, [7, 42])
        for a, b in zip(expected, actual):
            np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)


def test_anchor_support_endpoints_and_affine_coordinates():
    """Native origins and affine strain reconstruct physical interval endpoints."""
    positions = (0.1, 0.3, 0.8)
    interface = CTRfilm.EpitaxyInterface(
        _cell("top", positions=positions, ids=(11, 3, 91)),
        _cell("bottom", positions=positions, ids=(11, 3, 91)),
    )
    interface._set_start_layer(3)
    intervals = np.array([0.5, 0.3, 0.2])
    boundaries = np.concatenate(([0.0], np.cumsum(np.tile(intervals, 2))))
    domains = [[], [], []]
    scale, origin = 1.3, 1.7
    for index, boundary in enumerate(boundaries[:-1]):
        layer = index % 3
        matrix = np.column_stack((np.eye(3), np.zeros(3)))
        matrix[2, 2] = scale
        matrix[2, 3] = origin + scale * (
            boundary - interface.uc_top.layerpos[interface.layer_order[layer]]
        )
        domains[layer].append(matrix)
    for coordinate in (0.0, 1e-8, 1 - 1e-8, 1.0, 1 + 1e-8,
                       3 - 1e-8, 3.0, 3 + 1e-8, 6 - 1e-8, 6.0):
        expected = (origin + scale * np.interp(
            coordinate, np.arange(7), boundaries,
        )) * interface.uc_top.a[2]
        assert interface._ideal_top_anchor(
            coordinate, domains, intervals,
        ) == pytest.approx(expected, rel=0, abs=1e-12)
    for coordinate in (-1e-8, 6 + 1e-8, np.nan, np.inf):
        with pytest.raises(ValueError, match="outside the generated support"):
            interface._ideal_top_anchor(coordinate, domains, intervals)
    for invalid in ([], [[], [], []], [domains[0], domains[1], domains[2][:1]]):
        with pytest.raises(ValueError, match="complete support cycles"):
            interface._ideal_top_anchor(0.0, invalid, intervals)


@pytest.mark.parametrize("positions", [(0.0, 1.0), (0.0, np.nan)])
def test_invalid_upper_intervals_are_rejected(positions):
    """Invalid native geometry must not divide by zero or create a false anchor."""
    interface = CTRfilm.EpitaxyInterface(_cell("top"), _cell("bottom"))
    interface.layerpos = np.array(positions)
    with pytest.raises(ValueError, match="upper layer intervals"):
        interface.createInterfaceCells()


@pytest.mark.parametrize("termination", ["poisson", "unitcell"])
@pytest.mark.parametrize("interface_width,skew", [(0.15, 0.125), (1.0, 0.75)])
def test_text_and_fit_settings_roundtrip_preserve_references(
    termination, interface_width, skew,
):
    """Text geometry plus companion fit settings regenerate identical predictions."""
    # Exactly representable at the existing text format's decimal precision.
    crystal = _stack(repeat=6.5, registry=3, termination=termination,
                     transition={1: 2, 2: 1}, interface_width=interface_width)
    expected = _evaluate(crystal, skew)
    restored = CTRcalc.SXRDCrystal.fromStr(crystal.toStr(showErrors=False))
    # Exercise the same silx HDF5 dictionary format as companion fit files,
    # using an in-memory HDF5 file rather than a temporary filesystem path.
    with h5py.File("epitaxy-fit.h5", "w", driver="core", backing_store=False) as h5:
        dictdump.dicttoh5(crystal.parametersToDict(), h5)
        restored.parametersFromDict(dictdump.h5todict(h5))
    actual = _evaluate(restored, skew)
    for a, b in zip(expected, actual):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-9)
    for name in ("registry", "interface", "film"):
        np.testing.assert_array_equal(restored[name].basis_0, crystal[name].basis_0)
        np.testing.assert_array_equal(restored[name].basis, crystal[name].basis)
        assert restored[name].layer_state == crystal[name].layer_state
    assert restored["interface"].layer_transition.mapping == {1: 2, 2: 1}
    for name in ("uc_top", "uc_bottom"):
        assert getattr(restored["interface"], name).layerpos == getattr(
            crystal["interface"], name,
        ).layerpos


def test_support_extension_only_changes_truncated_occupancy():
    """Moving the support window transfers full upper sites to/from Film once."""
    crystals = [_stack(reverse=True, tail=tail) for tail in (1e-6, 1e-12)]
    results = [_evaluate(crystal, 0.0) for crystal in crystals]
    tables = []
    counts = []
    for crystal in crystals:
        interface, film = crystal["interface"], crystal["film"]
        counts.append(film._layers_to_create)
        table = {}
        for element, layers in (
            ("Si", interface.top_layers + film.layer_ucs),
            ("C", interface.bottom_layers + interface._reference_layers),
        ):
            for z, occupancy in _sites(layers):
                site = (element, round(z / (film.unitcell.a[2] / 2)))
                table[site] = table.get(site, 0.0) + occupancy
        tables.append(table)
    assert counts[0] != counts[1]  # The test actually changes generated support.
    differences = [
        abs(tables[0].get(site, 0.0) - tables[1].get(site, 0.0))
        for site in tables[0].keys() | tables[1].keys()
    ]
    # All fully occupied sites agree. Only discarded tails differ.
    assert max(differences) < 1e-6
    assert sum(differences) < 1e-6
    # 20 electrons bounds either atom's amplitude on these sampled rods.
    assert np.max(np.abs(results[0][1] - results[1][1])) < 20 * sum(differences) + 1e-9
    np.testing.assert_allclose(results[0][0][4:], results[1][0][4:],
                               rtol=0, atol=1e-9)

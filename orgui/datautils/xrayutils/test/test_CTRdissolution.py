"""Independent material and complete-state oracles for finite dissolution."""

import copy
import json
import pickle

from ..CTRfilm import LegacySurfaceBoundary

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from ..CTRcalc import SXRDCrystal
from ..CTRdistributions import (
    RecenteredSkellamProfile, HeightSupport, bounded_height_states,
)
from ._surface_boundary_helpers import _path
from ..CTRfilm import EpitaxyInterface, Film, SkellamSurface
from ..CTRincoherent import create_incoherent_model_from_config
from ._surface_boundary_helpers import _bounded_model
from .._CTRincoherent_dwba import _phase_to, _shared_preparation
from ._poisson_oracle import layered_cell


def _crystal(width=6, interface=False, substrate_film=False):
    bulk = layered_cell(2, "substrate")
    film_cell = layered_cell(2, "material")
    film_cell.basis[:, 6] *= .7
    film_cell.basis_0[:] = film_cell.basis
    components = []
    if substrate_film:
        components.append(Film(copy.deepcopy(bulk), width_layers=40,
                               name="substrate_film"))
    if interface:
        transition = EpitaxyInterface(copy.deepcopy(film_cell),
                                      copy.deepcopy(bulk), name="interface")
        transition.basis[:] = [.1, 0., 1., .03]
        transition.basis_0[:] = transition.basis
        components.append(transition)
    film = Film(film_cell, width_layers=width, name="film")
    # Build the independent pre-dissolution reference using legacy corrections.
    surface = SkellamSurface(copy.deepcopy(film_cell),
        boundary=LegacySurfaceBoundary(), name="surface")
    components.extend((film, surface))
    crystal = SXRDCrystal(bulk, *components,
                          stacking=np.arange(1, len(components)+1))
    crystal.setEnergy(10000)
    crystal.apply_stacking()
    crystal.dwba.set_ctr_geometry(alpha_i=.03)
    return crystal


def _flat(width, bulk=None, occupancy=.7):
    bulk = layered_cell(2, "substrate") if bulk is None else copy.deepcopy(bulk)
    components = []
    if width:
        cell = layered_cell(2, "flat")
        cell.basis[:, 6] *= occupancy
        cell.basis_0[:] = cell.basis
        components.append(Film(cell, width_layers=width, name="flat"))
    crystal = SXRDCrystal(bulk, *components, stacking=np.arange(len(components)))
    crystal.setEnergy(10000)
    if components:
        crystal.apply_stacking()
    crystal.dwba.set_ctr_geometry(alpha_i=.03)
    return crystal


def _hkl():
    return np.array([0., 1., 1., 1.]), np.zeros(4), np.array([.02, .43, .9, 1.99])


@pytest.mark.parametrize("up,down,offset", [(0., 0., -8.5), (0., 4., -2.3),
                                          (2., 0., -5.1), (1., 2., -4.2)])
def test_bounded_probability_and_moments(up, down, offset):
    """Cap an independent parent PMF; exhausted mass is not conditioning."""
    parent = RecenteredSkellamProfile(1, up, down, offset)
    indices, probabilities, certificate = bounded_height_states(parent, -6)
    raw = np.arange(-300, 301)
    pmf = parent.probability(raw)
    bounded = np.maximum(raw, -6)
    assert_allclose(probabilities.sum(), 1, atol=2e-15)
    assert_allclose(probabilities[0], pmf[raw <= -6].sum(), atol=3e-15)
    assert_allclose(indices @ probabilities, bounded @ pmf, atol=1e-10)
    assert_allclose(indices**2 @ probabilities, bounded**2 @ pmf, atol=1e-8)
    outside = raw > indices[-1]
    assert pmf[outside].sum() <= certificate.omitted_probability + 1e-16
    assert (raw[outside]**2 @ pmf[outside]
            <= certificate.omitted_second_moment + 1e-16)
    assert not probabilities.flags.writeable


def test_explicit_resource_limit():
    """The resource ceiling raises instead of clipping probability."""
    with pytest.raises(ValueError, match="max_states"):
        bounded_height_states(RecenteredSkellamProfile(1, 20, 20), -6,
                              support=HeightSupport(max_states=4))


@pytest.mark.parametrize("index", [-6, -5, -2, 0, 1, 4])
def test_deterministic_states_against_plain_films(index):
    """Each structural state equals an independently built finite Film."""
    original = _crystal()
    before = original.F(*_hkl())
    path = _path(original, films=["film"], surface="surface")
    state = path.state(index)
    expected = _flat(6+index)
    assert_allclose(state.crystal.F(*_hkl()), expected.F(*_hkl()),
                    rtol=1e-12, atol=1e-11)
    assert_allclose(state.front_height_A, 3*(6+index))
    assert_allclose(sum(state.inventory["film"].values()), .7*(6+index))
    assert_array_equal(original.F(*_hkl()), before)
    if index == -6:
        assert state.inventory["film"] == {}


@pytest.mark.parametrize("substrate_film", [False, True])
def test_interface_inventory_exhaustion_and_protected_material(substrate_film):
    """Film material in the mixed interface is removed without restoring STO."""
    crystal = _crystal(width=40, interface=True, substrate_film=substrate_film)
    path = _path(crystal, films=["film"], surface="surface",
                           exposure_policy="preserve")
    assert_allclose(path.state(0).crystal.F(*_hkl()), crystal.F(*_hkl()),
                    rtol=1e-12, atol=2e-11)
    exhausted = path.state(path.minimum)
    assert exhausted.inventory["film"] == {}
    # Independent endpoint: zero the upper material and Film contribution,
    # preserving all lower/reference records and their exact coordinates.
    expected = copy.deepcopy(crystal)
    expected.enable_uc_stacking = False
    expected._stack_policy = None
    expected.weights[-2:] = 0
    for cell in expected["interface"].top_layers:
        cell.coherentDomainOccupancy = np.zeros(len(cell.coherentDomainMatrix))
    assert_allclose(exhausted.crystal.F(*_hkl()), expected.F(*_hkl()),
                    rtol=2e-12, atol=2e-11)
    protected = {k: v for k, v in path.initial_inventory.items() if k != "film"}
    for index in [path.minimum, path.minimum+1, -2, 0]:
        state = path.state(index)
        for owner, atoms in protected.items():
            assert state.inventory[owner] == atoms


def test_explicit_substrate_film_can_dissolve():
    """A finite Film of bulk material is removable without making bulk finite."""
    crystal = _crystal(width=40, interface=True, substrate_film=True)
    path = _path(crystal, films=["substrate_film", "film"],
                           surface="surface", exposure_policy="preserve")
    endpoint = path.state(path.minimum)
    assert endpoint.inventory["film"] == {}
    assert endpoint.inventory["substrate_film"] == {}
    expected = _flat(0)
    assert_allclose(endpoint.crystal.F(*_hkl()), expected.F(*_hkl()),
                    rtol=1e-12, atol=1e-11)


@pytest.mark.parametrize("fraction", [0., .205425229351, 1.])
def test_complete_kinematic_ensemble(fraction):
    """An independent two-height ensemble includes full substrate interference."""
    crystal = _crystal()
    model = _bounded_model(
        crystal, films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(offset=-5.25),
        incoherent_fraction=fraction)
    f0, f1 = (_flat(n).F(*_hkl()) for n in (0, 1))
    coherent = abs(.25*f0+.75*f1)**2
    incoherent = .25*abs(f0)**2+.75*abs(f1)**2
    assert_allclose(model.F2(*_hkl()), (1-fraction)*coherent+fraction*incoherent,
                    rtol=1e-12, atol=1e-10)
    assert model.coherence_reconstruction_error(*_hkl()) < 1e-10
    stats = model.statistics()
    assert stats["exhausted_probability"] == .25
    assert_allclose(stats["mean_front_height_A"], 2.25)
    assert_allclose(stats["rms_front_height_A"], np.sqrt(1.6875))
    for species, initial in stats["initial_inventory"]["film"].items():
        assert_allclose(initial+stats["added_inventory"]["film"][species],
                        stats["remaining_inventory"]["film"].get(species, 0)
                        + stats["removed_inventory"]["film"][species])


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("fraction", [0., .205425229351, 1.])
def test_dwba_exhausted_deterministic_limit(policy, fraction):
    """Exhaustion recovers the native protected bulk for both field policies."""
    model = _bounded_model(
        _crystal(), films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(offset=-20),
        optical_reference=policy, incoherent_fraction=fraction)
    actual = model.dwba.evaluate(*_hkl())
    expected = _flat(0).dwba.evaluate(*_hkl())
    assert_allclose(actual.reflectivity, expected.reflectivity, rtol=1e-10, atol=1e-18)


@pytest.mark.parametrize("channels", [("s", "s"), ("s", "p"),
                                     ("p", "s"), ("p", "p")])
def test_state_specific_dwba_independent_states(channels):
    """State-specific fields mix total amplitudes in a common ambient gauge."""
    model = _bounded_model(
        _crystal(), films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(offset=-5.25), incoherent_fraction=.3)
    options = dict(polarization_i=channels[0], polarization_f=channels[1])
    result = model.dwba.evaluate(*_hkl(), **options)
    amplitude, intensity = 0., 0.
    for width, probability in [(0, .25), (1, .75)]:
        state = _flat(width).dwba.evaluate(*_hkl(), **options)
        field = state.total_amplitude*_phase_to(state.prepared, result.prepared)
        amplitude += probability*field
        intensity += probability*abs(field)**2
    assert_allclose(result.coherent_intensity, abs(amplitude)**2,
                    rtol=2e-10, atol=1e-18)
    assert_allclose(result.reflectivity, .7*abs(amplitude)**2+.3*intensity,
                    rtol=2e-10, atol=1e-18)


def test_config_copy_pickle_and_parameter_update():
    """The native fit/configuration APIs retain physical boundary policies."""
    model = _bounded_model(_crystal(), films=["film"], surface="surface",
                                 profile=RecenteredSkellamProfile(offset=-5.25))
    model.target_surface().addFitParameter(
        "offset", name="parent_offset", limits=(-20, 2))
    model.setParameters([-7.])
    expected = model.F2(*_hkl())
    config = json.loads(json.dumps(model.to_config(), default=lambda a: a.tolist()))
    for other in [copy.deepcopy(model), pickle.loads(pickle.dumps(model)),
                  create_incoherent_model_from_config(
                      copy.deepcopy(model.coherent_model), config)]:
        assert_allclose(other.F2(*_hkl()), expected)
        other.setParameters([-5.5])
        assert other.statistics()["exhausted_probability"] == .5
    assert model.statistics()["exhausted_probability"] == 1.


def test_invalid_ownership_and_recovery():
    """Invalid topology or rates fail explicitly and do not corrupt the source."""
    crystal = _crystal(width=40, interface=True, substrate_film=True)
    for films in [["substrate"], ["substrate_film"], ["film", "film"]]:
        with pytest.raises(ValueError):
            _path(crystal, films=films, surface="surface")
    model = _bounded_model(crystal, films=["film"], surface="surface")
    expected = model.F2(*_hkl())
    model.target_surface().basis[1] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        model.F2(*_hkl())
    model.target_surface().basis[1] = 0
    assert_allclose(model.F2(*_hkl()), expected)


def test_native_text_and_companion_config_roundtrip(tmp_path):
    """Native text plus required ensemble metadata preserve the physical model."""
    crystal = _crystal(width=40, interface=True)
    model = _bounded_model(
        crystal, films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(1., 0., 3., -35.),
        incoherent_fraction=.3)
    path = tmp_path / "dissolution.xtal"
    crystal.toFile(str(path))
    config = json.loads(json.dumps(model.to_config(), default=lambda a: a.tolist()))
    restored_crystal = SXRDCrystal.fromFile(str(path))
    restored_crystal.setEnergy(10000.)
    restored = create_incoherent_model_from_config(restored_crystal, config)
    # Native text has finite decimal precision; use the declared 1e-9
    # complex-amplitude scale budget, rather than a bitwise state comparison.
    assert_allclose(np.sqrt(restored.F2(*_hkl())),
                    np.sqrt(model.F2(*_hkl())), rtol=1e-9)
    assert_allclose(restored.statistics()["exhausted_probability"],
                    model.statistics()["exhausted_probability"], rtol=1e-12)


@pytest.mark.parametrize("dwba", [False, True])
def test_native_optimizer_parameter_updates(dwba):
    """The ordinary optimizer evaluates live dissolution parameters natively."""
    from ..CTRopt import CTROptimizer
    from ..CTRplotutil import CTR, CTRCollection
    model = _bounded_model(
        _crystal(), films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(offset=-5.25))
    model.target_surface().addFitParameter(
        "offset", name="parent_offset", limits=(-10., 2.))
    ell = np.array([.4, .6, .8])
    data = CTRCollection([CTR((1., 0.), ell, np.ones(3), np.ones(3))])
    if dwba:
        from .test_CTRopt_dwba import _ctr
        data = CTRCollection([_ctr(model.coherent_model, hk=(1., 0.),
                                   l_values=ell, quantity="reflectivity")])
    fit = CTROptimizer(model, data, scale_policy={"F": "fixed"})
    fit.set_dwba(dwba)
    fit.prepareFit()
    fit.set_parameters([-6.])
    prediction = fit.flat_prediction()
    reference = _flat(0)
    if dwba:
        angles = data[0].angles
        result = reference.dwba.evaluate_from_vlieg(
            *(angles[name] for name in angles.dtype.names))
        expected = result.reflectivity
    else:
        expected = abs(reference.F(np.ones(3), np.zeros(3), ell))
    assert_allclose(prediction, expected, rtol=1e-10, atol=1e-10)
    assert model.target_surface()._make_profile().offset == -5.25


def test_tight_support_keeps_existing_certificates():
    """C1 support tightening is explicit and does not relax the old budgets."""
    profile = RecenteredSkellamProfile(
        1, 0, 7.320904705207793, -3.107558383490982, 1e-6, .01)
    assert profile.support() == (-26, 5)
    lo, hi = profile.tight_support()
    assert lo+1 == -22
    probability, moment = profile.omitted_bounds(lo+1, hi+1)
    assert probability <= profile.tail_probability/4
    assert moment <= profile.moment_tolerance/4
    assert profile.support() == (-26, 5)


@pytest.mark.parametrize("remove_lower", [False, True])
@pytest.mark.parametrize("lateral", [(3., 3., 3., 3.), (3., 4., 5., 6.)])
def test_material_between_interfaces(remove_lower, lateral):
    """Signed replacement records preserve owners and lateral-area conversion."""
    from .test_CTRmultilayer import _parts
    bulk, components = _parts(lateral=lateral, roughness=.1)
    components.append(SkellamSurface(
        copy.deepcopy(components[-1].unitcell),
            boundary=LegacySurfaceBoundary(), name="surface"))
    crystal = SXRDCrystal(bulk, *components,
                          stacking=np.arange(1, len(components)+1))
    crystal.setEnergy(18000.)
    crystal.apply_stacking()
    owners = ["film_0", "film_1"] if remove_lower else ["film_1"]
    path = _path(crystal, films=owners, surface="surface",
                           exposure_policy="preserve")
    assert_allclose(path.state(0).crystal.F(*_hkl()), crystal.F(*_hkl()),
                    rtol=1e-12, atol=1e-10)
    exhausted = path.state(path.minimum)
    for owner in owners:
        assert exhausted.inventory[owner] == {}
    if remove_lower:
        expected = copy.deepcopy(crystal)
        expected.enable_uc_stacking = False
        expected._stack_policy = None
        expected.weights[1:] = 0
        for cell in expected["interface_0"].top_layers:
            cell.coherentDomainOccupancy = np.zeros(len(cell.coherentDomainMatrix))
        assert_allclose(exhausted.crystal.F(*_hkl()),
                        expected.F(*_hkl()),
                        rtol=1e-12, atol=1e-10)
    else:
        assert exhausted.inventory["film_0"] == path.initial_inventory["film_0"]


def test_upper_support_convergence_and_crystallinity():
    """Tighter numerical support converges native observables at fixed physics."""
    source = _crystal()
    source["film"].unitcell.crystallinity = .45
    parent = RecenteredSkellamProfile(1., 3., 5., -3., 1e-8, 1e-8)
    model = _bounded_model(
        source, films=["film"], surface="surface", profile=parent,
        incoherent_fraction=.2, exposure_policy="preserve")
    tighter = _bounded_model(
        source, films=["film"], surface="surface", profile=parent,
        incoherent_fraction=.2, exposure_policy="preserve",
        tail_probability=1e-12, moment_tolerance=1e-12,
        second_moment_tolerance=1e-12)
    assert_allclose(model.F2(*_hkl()), tighter.F2(*_hkl()), rtol=1e-11)
    assert model.coherence_reconstruction_error(*_hkl()) < 1e-10
    assert_allclose(model.statistics()["rms_front_height_A"],
                    tighter.statistics()["rms_front_height_A"], rtol=1e-12)


def test_shared_dwba_independent_state_fields():
    """The averaged optical medium is used consistently for every full state."""
    model = _bounded_model(
        _crystal(), films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(offset=-5.25),
        incoherent_fraction=.3, optical_reference="shared")
    with model.dwba.batch():
        result = model.dwba.evaluate(*_hkl())
        common = model.dwba._scope["average"].dwba.evaluate_prepared(result.prepared)
    amplitude, intensity = 0., 0.
    for width, probability in [(0, .25), (1, .75)]:
        state = _flat(width).dwba
        _, atomic, profile = state._atomic_model_snapshot()
        prepared = _shared_preparation(state, common.prepared, atomic, profile)
        value = state.evaluate_prepared(prepared)
        field = common.unperturbed_amplitude + common.amplitude_prefactor * (
            value.F_atomic*_phase_to(prepared, common.prepared)-common.F_reference)
        amplitude += probability*field
        intensity += probability*abs(field)**2
    assert_allclose(amplitude, common.total_amplitude, rtol=1e-10, atol=1e-14)
    assert_allclose(result.reflectivity,
                    .7*abs(amplitude)**2+.3*intensity, rtol=1e-10, atol=1e-18)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("change", ["profile", "occupancy", "geometry", "domain"])
def test_live_changes_match_fresh_and_batch(policy, change):
    """No state or wavefield from an earlier standalone call survives mutation."""
    model = _bounded_model(
        _crystal(), films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(offset=-5.25), optical_reference=policy,
        incoherent_fraction=.3)
    model.dwba.evaluate(*_hkl())
    if change == "profile":
        model.target_surface().basis[3] = -4.75
    elif change == "occupancy":
        film = model.coherent_model["film"].unitcell
        film.basis[:, 6] *= .9
        film.basis_0[:] = film.basis
    elif change == "geometry":
        model.coherent_model.dwba.set_ctr_geometry(alpha_i=.04)
    else:
        rotation = np.diag([-1., -1., 1.])
        model.coherent_model.domains = [[(rotation, 1.)] for _ in range(2)]
    fresh = copy.deepcopy(model)
    expected = fresh.dwba.evaluate(*_hkl()).reflectivity
    with model.dwba.batch():
        assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity, expected,
                        rtol=1e-12, atol=1e-18)
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity, expected,
                    rtol=1e-12, atol=1e-18)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_interface_exhaustion_dwba_and_statistics(policy):
    """Interface depletion preserves protected records in optical and atomic sums."""
    source = _crystal(width=40, interface=True)
    path = _path(source, films=["film"], surface="surface")
    model = _bounded_model(
        source, films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=path.minimum-2))
    expected = copy.deepcopy(source)
    expected.enable_uc_stacking = False
    expected._stack_policy = None
    expected["interface"].top_layers = []
    expected.uc_surface_list = expected.uc_surface_list[:1]
    expected.uc_surface_list_ordered = np.array(expected.uc_surface_list)
    expected.weights = np.ones(1)
    expected.domains = expected.domains[:1]
    expected.uc_stacking = np.array([1])
    expected.dwba.set_ctr_geometry(alpha_i=.03)
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity,
                    expected.dwba.evaluate(*_hkl()).reflectivity,
                    rtol=1e-9, atol=1e-18)
    stats = model.statistics()
    assert stats["exhausted_probability"] == 1
    assert stats["exposed_owner_probabilities"]["__bulk__"] == 1
    assert stats["interface_support_audit"]["interface"][
        "missing_plane_inventory_bound"] >= 0


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("interface", [False, True])
def test_nominal_dwba_matches_original_native_crystal(policy, interface):
    """The zero-removal state preserves atomic and optical native geometry."""
    source = _crystal(width=40, interface=interface)
    model = _bounded_model(
        source, films=["film"], surface="surface", optical_reference=policy)
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity,
                    source.dwba.evaluate(*_hkl()).reflectivity,
                    rtol=1e-9, atol=1e-18)


def test_motif_follows_material_across_support_partition():
    """Kinematic motif eligibility is independent of the component partition."""
    from ..CTRdistributions import SkellamProfile
    source = _crystal(width=40, interface=True)
    source["interface"].profile = SkellamProfile(.1, tail_probability=1e-14)
    source["interface"].createInterfaceCells()
    surface = source["surface"]
    surface.unitcell.basis[:, 6] *= .5
    surface.unitcell.basis_0[:] = surface.unitcell.basis
    surface.set_ucs(surface.unitcell)
    source.apply_stacking()
    expanded = copy.deepcopy(source)
    expanded["interface"].fixed_ucs = 20
    expanded["interface"].createInterfaceCells()
    options = dict(films=["film"], surface="surface",
                   profile=RecenteredSkellamProfile(offset=-18.5),
                   incoherent_fraction=.2)
    model, wider = (_bounded_model(c, **options)
                    for c in (source, expanded))
    assert model.coherent_model._surface_state_plan().top_film._layers_to_create != (
        wider.coherent_model._surface_state_plan().top_film._layers_to_create)
    assert_allclose(wider.F2(*_hkl()), model.F2(*_hkl()), rtol=1e-10)


def test_nonuniform_spacing_and_finite_substrate_partition():
    """The front follows native intervals; a substrate split adds no atoms."""
    source = _crystal()
    for component in [source["film"].unitcell, source["surface"].unitcell]:
        component.layerpos[1.] = .3
    source["film"].set_ucs(source["film"].unitcell)
    source["surface"].set_ucs(source["surface"].unitcell)
    source.apply_stacking()
    path = _path(source, films=["film"], surface="surface",
                           exposure_policy="preserve")
    assert_allclose([path.state(i).front_height_A for i in range(-6, 1)],
                    [0., 1.8, 6., 7.8, 12., 13.8, 18.])
    # A protected reference at the lower boundary and finite material above it
    # must equal the same perfect crystal after a common origin translation.
    bulk = layered_cell(2, "bulk")
    film = Film(copy.deepcopy(bulk), width_layers=4, name="finite_substrate")
    surface = SkellamSurface(copy.deepcopy(bulk),
        boundary=LegacySurfaceBoundary(), name="surface")
    # Native F_bulk includes repeat zero. Shift its atoms down by one complete
    # cell so the finite replacement starts directly above the protected bulk.
    bulk.basis[:, 3] -= 1
    bulk.basis_0[:] = bulk.basis
    crystal = SXRDCrystal(bulk, film, surface, stacking=np.array([1, 2]), atten=0.)
    crystal.setEnergy(10000)
    crystal.apply_stacking()
    path = _path(crystal, films=["finite_substrate"], surface="surface")
    h, k, l = _hkl()  # noqa: E741
    reference = _flat(0, bulk=bulk)
    reference.atten = 0
    assert_allclose(path.state(0).crystal.F(h, k, l),
                    reference.F(h, k, l)*np.exp(2j*np.pi*l*2),
                    rtol=1e-12, atol=1e-10)

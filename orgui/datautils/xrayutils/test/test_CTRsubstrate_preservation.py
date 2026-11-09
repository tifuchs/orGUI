"""The complete original substrate profile survives film dissolution."""

import copy
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose

from ..CTRcalc import SXRDCrystal
from ..CTRdistributions import HeightSupport, SkellamProfile
from ..CTRfilm import (
    EpitaxyInterface, Film, FiniteFilmBoundary, LegacySurfaceBoundary,
    PoissonSurface, SkellamSurface,
    _FrozenMaterial,
)
from ..CTRincoherent import SurfaceHeightDomains
from ..CTRuc import UnitCell
from .test_CTRdissolution import _hkl


def _perovskite(a, names, name):
    cell = UnitCell([a, a, a], [90., 90., 90.], name=name)
    for species, xyz, layer in [
        (names[0], [0., 0., 0.], 0), ("O", [.5, .5, 0.], 0),
        (names[1], [.5, .5, .5], 1), ("O", [0., .5, .5], 1),
        ("O", [.5, 0., .5], 1),
    ]:
        cell.addAtom(species, xyz, .5, .5, 1., layer=layer)
        cell.layerpos[float(layer)] = layer / 2
    return cell


def _source(kind=SkellamSurface):
    bulk = _perovskite(3.905, ("Sr", "Ti"), "STO")
    material = _perovskite(3.8, ("La", "Ni"), "LNO")
    interface = EpitaxyInterface(copy.deepcopy(material), copy.deepcopy(bulk),
                                 name="interface", profile=SkellamProfile(.6, -.7))
    interface.basis[:] = [.6, -.7, .2, .06]
    interface.basis_0 = interface.basis.copy()
    film = Film(material, width_layers=40, name="film")
    surface = kind(copy.deepcopy(material), name="surface")
    crystal = SXRDCrystal(bulk, interface, film, surface,
                          stacking=np.array([1, 2, 3]))
    crystal.setEnergy(20000)
    crystal.apply_stacking()
    crystal.dwba.set_ctr_geometry(alpha_i=.03)
    return crystal


def _original_substrate(source):
    """Independent native STO oracle, including its oxygen and reference terms."""
    original = copy.deepcopy(source)
    original["surface"].boundary = LegacySurfaceBoundary()
    original["surface"].basis[:] = original["surface"]._default_basis
    original.apply_stacking()
    interface = original["interface"]
    sto = _FrozenMaterial([*interface.bottom_layers, *interface._reference_layers],
                          "STO_interface", original.reference_area)
    expected = SXRDCrystal(original.uc_bulk, sto)
    expected.dwba.set_ctr_geometry(alpha_i=.03)
    return expected


def _substrate_of_state(state):
    """Select all native STO cells, not just cations or z-negative atoms."""
    components = []
    for component in state.uc_surface_list:
        # Each STO structural layer in this fixture contains Sr or Ti. Keep
        # the complete cell, including O and signed substrate references.
        native_cells = (component.layer_ucs if isinstance(component, Film)
                        else component.cells)
        cells = [copy.deepcopy(cell) for cell in native_cells
                 if any(name in {"Sr", "Ti"} for name in cell.names)]
        if cells:
            components.append(_FrozenMaterial(cells, component.name,
                                               state.reference_area))
    result = SXRDCrystal(copy.deepcopy(state.uc_bulk), *components)
    result.weights[:] = 1.
    return result


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
@pytest.mark.parametrize("hk", [(0., 0.), (1., 0.)])
def test_complete_sto_profile_is_invariant_on_both_sides_of_interface(kind, hk):
    """Full STO density and amplitude remain unchanged at every removal depth."""
    source = _source(kind)
    original = _original_substrate(source)
    path = source._surface_state_plan()
    z = np.linspace(-40, 30, 351)
    rho = original.zDensity_G(z, *hk)
    if hk == (0., 0.):
        bare = original.uc_bulk.zDensity_G_asbulk(z, *hk)
        # This catches the incorrect bare-bulk-only implementation: the
        # initial interface contains real substrate material above z=0.
        assert np.max((rho - bare).real[z > 0]) > .01
    for index in [path.minimum, path.minimum + 1, -20, 0, 2]:
        actual = _substrate_of_state(path.state(index).crystal)
        assert_allclose(actual.zDensity_G(z, *hk), rho, rtol=1e-12, atol=2e-12)
        assert_allclose(actual.F(*_hkl()), original.F(*_hkl()),
                        rtol=1e-12, atol=1e-10)
    # Probability averaging must preserve exactly the same complete profile.
    surface = source["surface"]
    if kind is SkellamSurface:
        surface.basis[:] = [1., 0., 5., -25.]
    else:
        surface.basis[:] = [-5., .4, -25.]
    average = _substrate_of_state(source._surface_state_plan()._average())
    assert_allclose(average.zDensity_G(z, *hk), rho, rtol=1e-12, atol=2e-12)


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("fraction", [0., .37, 1.])
def test_exhaustion_keeps_complete_sto_in_density_kinematic_and_dwba(
        kind, policy, fraction):
    """Exhaustion leaves the original rough substrate, including its extension."""
    source = _source(kind)
    expected = _original_substrate(source)
    source["surface"].basis[source["surface"].parameterLookup["offset"]] = -200
    model = SurfaceHeightDomains(source, surface="surface",
                                 optical_reference=policy,
                                 incoherent_fraction=fraction)
    z = np.linspace(-40, 30, 351)
    assert_allclose(source.zDensity_G(z, 0, 0), expected.zDensity_G(z, 0, 0),
                    rtol=1e-12, atol=2e-12)
    assert_allclose(model.F2(*_hkl()), expected.F2(*_hkl()),
                    rtol=1e-12, atol=1e-10)
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity,
                    expected.dwba.evaluate(*_hkl()).reflectivity,
                    rtol=1e-9, atol=1e-18)


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
def test_saved_models_retain_complete_original_substrate(kind):
    """Text, dictionary and pickle loading must not strip the STO extension."""
    source = _source(kind)
    source["surface"].basis[source["surface"].parameterLookup["offset"]] = -200
    source["surface"].basis_0 = source["surface"].basis.copy()
    from_dict = _source(kind)
    from_dict.parametersFromDict(source.parametersToDict())
    original = _original_substrate(source)
    z = np.linspace(-20, 20, 201)
    for restored in [SXRDCrystal.fromStr(source.toStr()), from_dict,
                     pickle.loads(pickle.dumps(source))]:
        assert_allclose(restored.zDensity_G(z, 0, 0),
                        original.zDensity_G(z, 0, 0), rtol=1e-10, atol=1e-10)


def _material_signature(crystal):
    cells = [crystal.uc_bulk, *(cell for component in crystal.uc_surface_list
                               for cell in component.cells)]
    return tuple((tuple(cell.names), cell.basis.tobytes(), cell.a.tobytes(),
                  cell.alpha.tobytes(),
                  np.asarray(cell.coherentDomainMatrix).tobytes(),
                  np.asarray(cell.coherentDomainOccupancy).tobytes())
                 for cell in cells)


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
@pytest.mark.parametrize("change", ["distribution", "offset", "support", "film", "cap"])
def test_substrate_is_independent_of_unrelated_parameters(kind, change):
    """Unrelated changes cannot translate, duplicate or reweight protected atoms."""
    source = _source(kind)
    reference = _substrate_of_state(source._surface_state_plan().state(0).crystal)
    surface = source["surface"]
    surface.basis[:] = ([1., 0., 5., -25.] if kind is SkellamSurface
                        else [-5., .4, -25.])
    if change == "distribution":
        surface.basis[:] = ([2., 3., 4., -25.] if kind is SkellamSurface
                            else [-9., .8, -25.])
    elif change == "offset":
        surface.basis[-1] = -48.
    elif change == "support":
        surface.height_support = HeightSupport(
            tail_probability=1e-13, moment_tolerance=1e-11,
            second_moment_tolerance=1e-10)
    elif change == "film":
        source["film"].basis[0] += 4
    elif change == "cap":
        cap = Film(copy.deepcopy(source["film"].unitcell), width_layers=2,
                   name="cap")
        source = SXRDCrystal(source.uc_bulk, *source.uc_surface_list, cap,
                             stacking=np.arange(1, 5))
        source.setEnergy(20000)
    actual = _substrate_of_state(source._surface_state_plan()._average())
    assert _material_signature(actual) == _material_signature(reference)
    z = np.linspace(-25, 15, 101)
    np.testing.assert_array_equal(actual.zDensity_G(z, 0, 0),
                                  reference.zDensity_G(z, 0, 0))
    np.testing.assert_array_equal(actual.F(*_hkl()), reference.F(*_hkl()))
    np.testing.assert_array_equal(actual.optical_profile(), reference.optical_profile())


@pytest.mark.parametrize("parameter", ["width", "skew", "coupling", "offset"])
def test_substrate_interface_parameters_still_update_protected_material(parameter):
    """Protection is relative to dissolution, not a permanent geometry cache."""
    source = _source()
    before = source._surface_state_plan()
    reference = before._crystal(before._protected)
    index = ["width", "skew", "coupling", "offset"].index(parameter)
    source["interface"].basis[index] += .05
    after = source._surface_state_plan()
    actual = after._crystal(after._protected)
    assert _material_signature(actual) != _material_signature(reference)
    z = np.linspace(-25, 15, 101)
    assert np.max(abs(actual.zDensity_G(z, 0, 0)
                      - reference.zDensity_G(z, 0, 0))) > 1e-7


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
def test_finite_substrate_partition_preserves_initial_structure(kind):
    """A selectable substrate Film replaces bulk repeats without adding atoms."""
    source = _source(kind)
    n = len(source.uc_bulk.layers)
    source.uc_bulk = source.uc_bulk.translate_layered([0, 0, -n])
    source.atten = 0.  # Exact bulk/finite amplitude identity, no bulk damping.
    source.apply_stacking()
    low = source["interface"]._support_bounds[0]
    repeats = int(np.ceil(-low / source.uc_bulk.a[2])) + 1
    depth = repeats * source.uc_bulk.a[2]
    cell = source.uc_bulk.translate_layered([0, 0, n])
    film = Film(cell, width_layers=repeats*n, name="finite_STO")
    split = SXRDCrystal(copy.deepcopy(source.uc_bulk), film,
                        *copy.deepcopy(source.uc_surface_list),
                        stacking=np.arange(1, 5), atten=0.)
    split["surface"].boundary = FiniteFilmBoundary(("film",))
    original = source._surface_state_plan().state(0).crystal
    path = split._surface_state_plan()
    initial = path.state(0).crystal
    # Native split coordinates start at the bottom of the finite substrate.
    # Plotting/evaluation coordinates z+depth recover the original origin.
    z = np.linspace(-depth-5, 90, 401)
    assert_allclose(initial.zDensity_G(z+depth, 0, 0),
                    original.zDensity_G(z, 0, 0), rtol=1e-11, atol=1e-11)
    h, k, ell = _hkl()
    phase = np.exp(2j*np.pi*ell*depth/source.reference_uc.a[2])
    assert_allclose(initial.F(h, k, ell), original.F(h, k, ell)*phase,
                    rtol=1e-11, atol=1e-9)
    source["surface"].basis[-1] = split["surface"].basis[-1] = -200
    assert_allclose(split.zDensity_G(z+depth, 0, 0),
                    source.zDensity_G(z, 0, 0), rtol=1e-11, atol=1e-11)
    # The same initial structure can now explicitly permit finite STO removal.
    split["surface"].boundary = FiniteFilmBoundary()
    removable = split._surface_state_plan()
    assert_allclose(removable.state(0).crystal.zDensity_G(z+depth, 0, 0),
                    original.zDensity_G(z, 0, 0), rtol=1e-11, atol=1e-11)
    exhausted = removable.state(removable.minimum)
    assert exhausted.inventory["finite_STO"] == {}
    assert_allclose(exhausted.crystal.zDensity_G(z+depth, 0, 0),
                    split.uc_bulk.zDensity_G_asbulk(z+depth, 0, 0),
                    rtol=1e-11, atol=1e-11)

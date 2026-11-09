"""Physical closure and fixed registry of complementary interface populations."""

import copy
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose

from .test_CTRinterface_anchors import _stack, _evaluate, _cell
from ..CTRcalc import SXRDCrystal
from ..CTRdistributions import SkellamProfile
from ..CTRfilm import EpitaxyInterface, Film


@pytest.mark.parametrize("coupling", [0., .4, 1.])
@pytest.mark.parametrize("skew", [-.8, 0., .7])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("positions,ids", [
    ((0.,), (1,)), ((0., .5), (1, 2)), ((.1, .3, .8), (7, 42, 13)),
])
def test_matched_lattices_share_sites_without_registry_offset(
        coupling, skew, reverse, positions, ids):
    """Zero offset must not turn skew into relative lattice displacement."""
    crystal = _stack(coupling=coupling, interface_width=.6, reverse=reverse,
                     positions=positions, ids=ids)
    _evaluate(crystal, skew)
    interface = crystal["interface"]
    for top, bottom in zip(interface.top_layers, interface.bottom_layers):
        count = len(top._crystallinity_domain_offsets)
        for mt, mb in zip(top.coherentDomainMatrix[:count],
                          bottom.coherentDomainMatrix[:count]):
            z_top = (mt[2, :3] @ top.basis[0, 1:4] + mt[2, 3]) * top.a[2]
            z_bottom = (mb[2, :3] @ bottom.basis[0, 1:4] + mb[2, 3]) * bottom.a[2]
            assert_allclose(z_top, z_bottom, rtol=0, atol=1e-11)
        assert_allclose(top.coherentDomainOccupancy[:count]
                        + bottom.coherentDomainOccupancy[:count], 1., atol=1e-15)


@pytest.mark.parametrize("coupling", [0., .4, 1.])
def test_offset_is_the_only_residual_registration_difference(coupling):
    """Matched lattices differ by (1-kappa)*offset*c, independent of skew."""
    for skew in [-.8, 0., .7]:
        crystal = _stack(coupling=coupling, interface_width=.6, offset=.13,
                         reverse=True, positions=(0.,), ids=(1,))
        _evaluate(crystal, skew)
        interface = crystal["interface"]
        top, bottom = interface.top_layers[0], interface.bottom_layers[0]
        count = len(top._crystallinity_domain_offsets)
        actual = [(mt[2, 3] - mb[2, 3]) * top.a[2]
                  for mt, mb in zip(top.coherentDomainMatrix[:count],
                                    bottom.coherentDomainMatrix[:count])]
        assert_allclose(actual, (1-coupling)*.13*bottom.a[2], rtol=0, atol=1e-11)


@pytest.mark.parametrize("coupling", [0., .4, 1.])
@pytest.mark.parametrize("skew", [-.8, .7])
def test_identical_material_has_no_artificial_interface(coupling, skew):
    """An interface within one ideal crystal equals a plain finite Film."""
    crystal = _stack(coupling=coupling, interface_width=.6, reverse=True,
                     positions=(0.,), ids=(1,))
    interface = crystal["interface"]
    cell = copy.deepcopy(interface.uc_bottom)
    interface.set_ucs(copy.deepcopy(cell), interface.uc_bottom)
    crystal["film"].set_ucs(copy.deepcopy(cell))
    crystal["surface"].set_ucs(copy.deepcopy(cell))
    crystal.setEnergy(18000.)
    _evaluate(crystal, skew)
    count = crystal["film"].height_absolute / cell.a[2]
    assert_allclose(count, round(count), atol=1e-12)
    expected = SXRDCrystal(copy.deepcopy(crystal.uc_bulk),
                           Film(cell, width_layers=round(count)),
                           stacking=np.array([1]))
    expected.setEnergy(18000.)
    h, k, ell = np.array([0., 1., 1.]), 0., np.array([.03, .37, 1.21])
    assert_allclose(crystal.F(h, k, ell), expected.F(h, k, ell),
                    rtol=1e-12, atol=1e-10)
    z = np.linspace(-35, 45, 401)
    assert_allclose(crystal.zDensity_G(z, 0, 0), expected.zDensity_G(z, 0, 0),
                    rtol=1e-12, atol=1e-11)
    for model in (crystal, expected):
        model.dwba.set_ctr_geometry(alpha_i=.03)
    assert_allclose(crystal.dwba.evaluate(h, k, ell).reflectivity,
                    expected.dwba.evaluate(h, k, ell).reflectivity,
                    rtol=2e-10, atol=1e-18)


@pytest.mark.parametrize("coupling", [0., .4, 1.])
@pytest.mark.parametrize("positions,ids", [((0.,), (1,)), ((.1, .3, .8), (7, 42, 13))])
def test_fixed_origin_preserves_oblique_repeat_registry(coupling, positions, ids):
    """Matched oblique repeat vectors align laterally as well as vertically."""
    crystal = _stack(coupling=coupling, interface_width=.6, reverse=True,
                     positions=positions, ids=ids)
    interface = crystal["interface"]
    for cell in (crystal.uc_bulk, interface.uc_top, interface.uc_bottom,
                 crystal["film"].unitcell, crystal["surface"].unitcell):
        cell.stacking_translation = (.25, -.125, 1.)
    _evaluate(crystal, -.8)
    for top, bottom in zip(interface.top_layers, interface.bottom_layers):
        count = len(top._crystallinity_domain_offsets)
        for mt, mb in zip(top.coherentDomainMatrix[:count],
                          bottom.coherentDomainMatrix[:count]):
            assert_allclose(mt[:, :3] @ top.basis[0, 1:4] + mt[:, 3],
                            mb[:, :3] @ bottom.basis[0, 1:4] + mb[:, 3],
                            rtol=0, atol=1e-12)


def test_old_cached_interface_domains_are_rebuilt():
    """A saved old domain cache cannot retain the previous registration."""
    crystal = _stack(coupling=.4, interface_width=.6)
    _evaluate(crystal, -.8)
    interface = crystal["interface"]
    hkl = tuple(np.array([value]) for value in (1., 0., .37))
    expected = interface.F_uc(*hkl)
    old = copy.deepcopy(interface)
    old._repeat_created = old._repeat_signature()[1]
    for cell in old.top_layers:
        for matrix in cell.coherentDomainMatrix:
            matrix[2, 3] += .123
    restored = pickle.loads(pickle.dumps(old))
    assert_allclose(restored.F_uc(*hkl), expected, rtol=1e-12, atol=1e-12)




@pytest.mark.parametrize("positions,ids", [
    ((0., .5), (1, 2)), ((.1, .3, .8), (7, 42, 13)),
])
@pytest.mark.parametrize("oblique", [False, True])
def test_independent_upper_lattice_does_not_slide_with_skew(positions, ids, oblique):
    """Skew changes occupancies, not the phase of an independent lattice."""
    crystal = _stack(coupling=0., mismatch=.85, registry=3,
                     interface_width=.6, positions=positions, ids=ids)
    interface = crystal["interface"]
    if oblique:
        for cell in (crystal.uc_bulk, interface.uc_top, interface.uc_bottom,
                     crystal["registry"].unitcell, crystal["film"].unitcell,
                     crystal["surface"].unitcell):
            cell.stacking_translation = (.25, -.125, 1.)
    samples = []
    for skew in (-.8, 0., .7):
        _evaluate(crystal, skew)
        sites = {}
        for cell in interface.top_layers:
            for index, matrix in zip(cell._crystallinity_domain_offsets,
                                     cell.coherentDomainMatrix):
                sites[int(index)] = (matrix[:, :3] @ cell.basis[0, 1:4]
                                     + matrix[:, 3]) * cell.a
        samples.append(sites)
    for index in set.intersection(*(set(sites) for sites in samples)):
        for sites in samples[1:]:
            assert_allclose(sites[index], samples[0][index], rtol=0, atol=1e-11)


@pytest.mark.parametrize("positions", [(0., .5), (0., .2, .7)])
@pytest.mark.parametrize("oblique", [False, True])
@pytest.mark.parametrize("skew", [-.8, .7])
@pytest.mark.parametrize("coupling", [.4, 1.])
def test_multilayer_identical_crystal_matches_plain_native_dwba(
        positions, oblique, skew, coupling):
    """Unequal and oblique cycles close against an independent plain Film."""
    n = len(positions)
    cell = _cell("material", positions=positions, ids=tuple(range(1, n+1)))
    if oblique:
        cell.stacking_translation = (.25, -.125, 1.)
    bulk = cell.translate_layered([0, 0, -n], name="bulk")
    interface = EpitaxyInterface(copy.deepcopy(cell), copy.deepcopy(cell),
                                 profile=SkellamProfile(.6, skew), name="interface")
    interface.basis[2] = coupling
    film = Film(copy.deepcopy(cell), width_layers=48, name="film")
    crystal = SXRDCrystal(bulk, interface, film, stacking=np.array([1, 2]))
    crystal.setEnergy(18000.)
    crystal.apply_stacking()
    # Enumerate the plain lattice up to the same external top boundary. This
    # reference contains no interface domains, probabilities or signed terms.
    top = film.height_absolute
    sites = [(repeat+position)*cell.a[2] for repeat in range(100)
             for position in positions]
    count = sum(z < top-1e-9 for z in sites)
    plain = SXRDCrystal(copy.deepcopy(bulk),
                        Film(copy.deepcopy(cell), width_layers=count),
                        stacking=np.array([1]))
    plain.setEnergy(18000.)
    z = np.linspace(-30, top+5, 501)
    assert_allclose(crystal.zDensity_G(z, 0, 0), plain.zDensity_G(z, 0, 0),
                    rtol=1e-12, atol=1e-11)
    h, k, ell = np.array([0., 1., 2.]), np.zeros(3), np.array([.13, .37, 1.21])
    assert_allclose(crystal.F(h, k, ell), plain.F(h, k, ell),
                    rtol=1e-12, atol=1e-10)
    for model in (crystal, plain):
        model.dwba.set_ctr_geometry(alpha_i=.03)
    assert_allclose(crystal.dwba.evaluate(h, k, ell).reflectivity,
                    plain.dwba.evaluate(h, k, ell).reflectivity,
                    rtol=2e-10, atol=1e-18)

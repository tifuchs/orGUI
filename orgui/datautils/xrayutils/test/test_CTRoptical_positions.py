"""Empty interface bookkeeping must not displace occupied optical centers."""

import copy
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from .. import CTRcalc, CTRdwba, CTRoptics
from ..CTRdistributions import RecenteredSkellamProfile, SkellamProfile
from ._surface_boundary_helpers import _bounded_model
from .test_CTRdissolution import _crystal, _hkl


@pytest.mark.parametrize("contribution", [
    [2., 3.], [-2., -3.], [0., 1e-300], [1e-300, 0.],
])
def test_empty_origin_does_not_move_nonzero_profile(contribution):
    """Use exact nonzero contributions, including signed and tiny tails."""
    values = np.array([[3., 0., -0.], [3.18, *contribution], [6., 0., 0.]])
    before = values.copy()
    expected = values[1:].copy()
    assert_array_equal(CTRoptics.combine_profiles(values), expected)
    combined, shares = CTRdwba._combine_tagged(values, values[None, :, 1:])
    assert_array_equal(combined, expected)
    assert_array_equal(shares[0], expected[:, 1:])
    assert_array_equal(values, before)


def test_signed_cancellation_and_empty_groups_keep_their_origins():
    """Cancellation does not make contributors empty; real gaps remain."""
    values = np.array([[3., 0., 0.], [3.18, 2., 3.], [3.3, -2., -3.],
                       [6., 0., 0.], [6.2, 0., 0.], [9., 4., 5.]])
    expected = np.array([[3.18, 0., 0.], [6., 0., 0.], [9., 4., 5.]])
    assert_array_equal(CTRoptics.combine_profiles(values), expected)
    shares = np.zeros((2, len(values), 2))
    shares[0, 1] = values[1, 1:]
    shares[1, 2] = values[2, 1:]
    shares[0, -1] = values[-1, 1:]
    combined, tagged = CTRdwba._combine_tagged(values, shares)
    assert_array_equal(combined, expected)
    assert_array_equal(tagged[:, 0], [[2., 3.], [-2., -3.]])
    assert_array_equal(tagged.sum(axis=1), shares.sum(axis=1))


def test_tagged_reference_only_share_selects_coordinate():
    """Areal reference inventory remains active even without net density."""
    values = np.array([[3., 0., 0.], [3.18, 0., 0.], [6., 0., 0.]])
    shares = np.zeros((2, 3, 2))
    shares[1, 1] = [-1e-300, 0.]
    combined, tagged = CTRdwba._combine_tagged(values, shares, density_count=1)
    assert_array_equal(combined, [[3.18, 0., 0.], [6., 0., 0.]])
    assert_array_equal(tagged[:, 0], shares[:, 1])


def test_opposing_tags_at_one_origin_are_not_empty():
    """Net-zero density retains the origin and both signed reference shares."""
    values = np.array([[3., 0., 0.], [3.18, 0., 0.], [6., 0., 0.]])
    shares = np.zeros((2, 3, 2))
    shares[:, 1] = [[2., 3.], [-2., -3.]]
    combined, tagged = CTRdwba._combine_tagged(values, shares)
    assert_array_equal(combined, [[3.18, 0., 0.], [6., 0., 0.]])
    assert_array_equal(tagged[:, 0], shares[:, 1])


def test_only_marked_zero_domains_are_removed_before_grouping():
    """Preserve user empties and every nonzero signed or tiny domain."""
    cell = SimpleNamespace(
        coherentDomainOccupancy=[0., 1., -1., 1e-300, 0.],
        _optical_bookkeeping_domains=[True, True, True, True, False],
    )
    profile = np.array([[0., 0., 0.], [.5, 2., 0.], [.9, -2., 0.],
                        [2., 1e-300, 0.], [4., 0., 0.]])
    filtered = CTRoptics._without_empty_bookkeeping(cell, profile)
    assert_array_equal(filtered, profile[1:])
    assert_array_equal(CTRoptics.combine_profiles(filtered),
                       [[.5, 0., 0.], [2., 1e-300, 0.], [4., 0., 0.]])
    # Public array-only calls cannot distinguish either kind of zero.
    assert_array_equal(CTRoptics.combine_profiles(profile)[:2, 0], [.5, .9])


def test_user_zero_occupancy_cell_keeps_existing_grouping():
    """A user-supplied empty cell retains its existing optical interpretation."""
    def cell(z, occupancy, name):
        result = CTRcalc.UnitCell([3., 3., 6.], [90., 90., 90.], name=name)
        result.addAtom("C", [0., 0., 0.], .1, .1, occupancy)
        matrix = np.eye(3, 4)
        matrix[2, 3] = z / 6.
        result.coherentDomainMatrix = [matrix]
        return result

    crystal = CTRcalc.SXRDCrystal(cell(0., 1., "bulk"), cell(6., 0., "empty"),
                                  cell(6.5, 1., "first"),
                                  cell(6.9, 1., "second"))
    crystal.weights[:] = 1.
    crystal.setEnergy(10000.)
    for profile in (crystal.optical_profile(),
                    crystal.dwba._atomic_model_snapshot(raw_structure=True)):
        assert_allclose(profile[-3:-1, 0], [6.5, 6.9], atol=1e-12)


def _support_pair():
    """Same atoms, but additional zero-weight interface subtraction records."""
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
    expanded.apply_stacking()
    return source, expanded


def _offset_support_pair():
    source, expanded = _support_pair()
    for crystal in (source, expanded):
        crystal["interface"].basis[3] = .15  # 0.9 Angstrom displacement
        crystal["interface"].createInterfaceCells()
        crystal.apply_stacking()
    return source, expanded


def test_generated_subtraction_slot_cannot_split_occupied_group():
    """Real generated zero slots cannot determine occupied group membership."""
    _, expanded = _offset_support_pair()
    cell = expanded["interface"].bottom_layers[0]
    raw = cell.optical_profile()
    marked = np.asarray(cell._optical_bookkeeping_domains)
    index = np.flatnonzero(marked & (cell.coherentDomainOccupancy == 0))[-1]
    empty = raw[index:index + 1]
    active = np.array([[empty[0, 0] + .5, 1e-6, 0.],
                       [empty[0, 0] + .9, 2e-6, 0.]])
    filtered = CTRoptics._without_empty_bookkeeping(cell, raw)
    assert not np.any(filtered[:, 0] == empty[0, 0])
    # Retaining this generated zero changes the occupied grouping.
    assert len(CTRoptics.combine_profiles(empty, active)) == 2
    expected = CTRoptics.combine_profiles(active)
    assert_array_equal(expected, [[active[0, 0], 3e-6, 0.]])
    tagged, shares = CTRdwba._combine_tagged(active, active[None, :, 1:])
    assert_array_equal(tagged, expected)
    assert_array_equal(shares[0], expected[:, 1:])


def test_generated_sharp_subtraction_support_does_not_change_optics():
    """Extra zero sharp-subtraction slots do not create optical media."""
    source, expanded = _offset_support_pair()
    assert_allclose(expanded.F(*_hkl()), source.F(*_hkl()), rtol=1e-11)
    for getter in (lambda c: c.optical_profile(),
                   lambda c: c.dwba._atomic_model_snapshot(raw_structure=True)):
        assert_allclose(getter(expanded), getter(source), rtol=1e-11, atol=1e-12)
    narrow, wide = [c.dwba.evaluate(*_hkl()) for c in (source, expanded)]
    assert_allclose(wide.total_amplitude, narrow.total_amplitude, rtol=1e-9)


def test_saved_interface_reconstructs_bookkeeping_provenance():
    """Older interfaces need no serialized bookkeeping marker."""
    _, crystal = _offset_support_pair()
    expected = crystal.dwba._atomic_model_snapshot(raw_structure=True)
    for cell in crystal["interface"].bottom_layers:
        del cell._optical_bookkeeping_domains
    assert_allclose(crystal.dwba._atomic_model_snapshot(raw_structure=True),
                    expected, rtol=0, atol=0)
    assert all(hasattr(cell, "_optical_bookkeeping_domains")
               for cell in crystal["interface"].bottom_layers)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_bounded_offset_support_preserves_dwba(policy):
    """Rebuild bookkeeping provenance for averaged and exhausted states."""
    models = [_bounded_model(
        c, films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=-18.5), incoherent_fraction=.2)
        for c in _offset_support_pair()]
    narrow, wide = [m.dwba.evaluate(*_hkl()) for m in models]
    assert_allclose(wide.reflectivity, narrow.reflectivity, rtol=1e-9, atol=1e-15)


def test_interface_support_preserves_optical_centers_and_boundaries():
    """Extra empty subtraction planes preserve public and tagged profiles."""
    source, expanded = _support_pair()
    assert_allclose(expanded.F(*_hkl()), source.F(*_hkl()), rtol=1e-11)
    for getter in (lambda c: c.optical_profile(),
                   lambda c: c.dwba._atomic_model_snapshot(raw_structure=True)):
        narrow, wide = getter(source), getter(expanded)
        # Repartitioned stacking sums the same translations in a different
        # order. Allow a few ulps in Angstrom, not an optical merge tolerance.
        roundoff = 4 * np.finfo(float).eps * np.max(abs(narrow[:, 0]))
        assert_allclose(wide[:, 0], narrow[:, 0], rtol=0, atol=roundoff)
        assert_allclose(wide[:, 1:], narrow[:, 1:], rtol=1e-12, atol=1e-18)
        assert_allclose(CTRoptics.profile_boundaries(wide),
                        CTRoptics.profile_boundaries(narrow), rtol=0,
                        atol=roundoff)


@pytest.mark.parametrize("polarization", ["s", "p"])
def test_native_dwba_support_invariant(polarization):
    """Native complex amplitudes and intensities agree across support changes."""
    source, expanded = _support_pair()
    narrow, wide = [c.dwba.evaluate(*_hkl(), polarization_i=polarization,
                                  polarization_f=polarization)
                    for c in (source, expanded)]
    assert_allclose(wide.total_amplitude, narrow.total_amplitude, rtol=1e-9)
    assert_allclose(wide.reflectivity, narrow.reflectivity, rtol=1e-9)


@pytest.mark.parametrize("polarization", ["s", "p"])
@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("incoherent_fraction", [0., .2, 1.])
def test_dissolution_dwba_support_invariant(
        policy, polarization, incoherent_fraction):
    """Both field policies preserve coherent, incoherent and mixed observables."""
    models = [_bounded_model(
        c, films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=-18.5),
        incoherent_fraction=incoherent_fraction) for c in _support_pair()]
    narrow, wide = [m.dwba.evaluate(*_hkl(), polarization_i=polarization,
                                  polarization_f=polarization)
                    for m in models]
    assert_allclose(wide.reflectivity, narrow.reflectivity, rtol=1e-9)

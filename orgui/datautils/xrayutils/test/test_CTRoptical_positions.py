"""Empty interface bookkeeping must not displace occupied optical centers."""

import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from .. import CTRdwba, CTRoptics
from ..CTRdistributions import RecenteredSkellamProfile, SkellamProfile
from ..CTRincoherent import FiniteFilmDissolution
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
    models = [FiniteFilmDissolution(
        c, films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=-18.5),
        incoherent_fraction=incoherent_fraction) for c in _support_pair()]
    narrow, wide = [m.dwba.evaluate(*_hkl(), polarization_i=polarization,
                                  polarization_f=polarization)
                    for m in models]
    assert_allclose(wide.reflectivity, narrow.reflectivity, rtol=1e-9)

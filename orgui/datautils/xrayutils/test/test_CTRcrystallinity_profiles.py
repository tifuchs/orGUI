"""Local order recovery is conditioned on material existence, not counted twice."""

import copy
import pickle

import numpy as np
import pytest

from .. import CTRcalc, CTRfilm
from ..CTRcrystallinity import CrystallinityProfile
from ..CTRdistributions import PoissonProfile, SkellamProfile
from ._poisson_oracle import flat_height_crystal, height_states, poisson_crystal
from .test_CTRcrystallinity import _cell, _hkl, _prepare


def _flat_ordered_amplitude(height, profile, width=16):
    crystal = flat_height_crystal(int(height), w_base=width)
    crystal.setEnergy(10000.0)
    crystal.apply_stacking()
    film = crystal.uc_surface_list[0]
    layers = []
    for cell in film.layer_ucs:
        for index, matrix in enumerate(cell.coherentDomainMatrix):
            z = matrix[2, 2] * cell._optical_layer_origin + matrix[2, 3]
            layers.append((z, cell, index))
    for depth, (_, cell, index) in enumerate(sorted(layers, key=lambda row: -row[0])):
        cell.coherentDomainOccupancy[index] *= profile.evaluate(depth, 2)
    # Evaluate directly: SXRDCrystal.F re-stacks and regenerates domain weights.
    return crystal.uc_bulk.F_bulk(*_hkl(), crystal.atten) + film.F_uc(*_hkl())


@pytest.mark.parametrize("family", [PoissonProfile(2, alpha=0), SkellamProfile(0)])
def test_profile_depth_conventions(family):
    """A deterministic two-layer skin suppresses exactly the first two layers."""
    offset = 2 if isinstance(family, SkellamProfile) else 0
    profile = CrystallinityProfile(family, offset=offset)
    np.testing.assert_array_equal(profile.evaluate([0, 1, 2, 3], 2), [0, 0, 1, 1])
    profile.basis[-1] = 0.25
    np.testing.assert_array_equal(profile.evaluate([0, 1, 2], 2), [0.25, 0.25, 1])


@pytest.mark.parametrize(
    "height_profile",
    [
        PoissonProfile(0, alpha=0),
        PoissonProfile(1.7, alpha=0.6),
        PoissonProfile(-1.4, alpha=0.7),
        PoissonProfile(0, offset=0.35),
    ],
)
@pytest.mark.parametrize("depth", [1.0, 3.0])
def test_surface_matches_explicit_local_terrace_sum(height_profile, depth):
    """Average disordered flat terraces, including the buried film layers."""
    crystal, surface = poisson_crystal(height_profile, w_base=16)
    crystal.setEnergy(10000.0)
    optics = crystal.optical_profile().copy()
    order = CrystallinityProfile(PoissonProfile(depth, alpha=0.4))
    surface.setCrystallinityProfile(order)
    heights, masses = height_states(height_profile)
    expected = sum(
        mass * _flat_ordered_amplitude(h, order) for h, mass in zip(heights, masses)
    )
    np.testing.assert_allclose(crystal.F(*_hkl()), expected, atol=2e-8, rtol=2e-9)
    np.testing.assert_array_equal(crystal.optical_profile(), optics)


def _interface():
    interface = CTRfilm.EpitaxyInterface(
        _cell("top"),
        _cell("bottom"),
        profile=SkellamProfile(0.8, 0.2),
        fixed_ucs=32,
    )
    crystal = CTRcalc.SXRDCrystal(_cell("bulk"), interface, stacking=np.array([1]))
    return crystal, interface


@pytest.mark.parametrize("side", ["top", "bottom"])
def test_interface_order_is_conditional_on_local_boundary(side):
    """Explicit boundary masses agree with the conditional order of each occurrence."""
    crystal, interface = _interface()
    crystal.apply_stacking()
    order = CrystallinityProfile(PoissonProfile(2.3, alpha=0.7))
    interface.setCrystallinityProfile(order, side=side)
    boundary = np.arange(-70, 71)
    distribution = SkellamProfile(*interface.basis[:2])
    mass = distribution.occupancy(boundary, 2) - distribution.occupancy(boundary - 1, 2)
    for offset in range(-5, 6):
        depth = offset - boundary if side == "top" else boundary - 1 - offset
        present = depth >= 0
        expected = np.dot(mass[present], order.evaluate(depth[present], 2))
        probability = (
            distribution.occupancy([offset], 2)[0]
            if side == "top"
            else 1 - distribution.occupancy([offset], 2)[0]
        )
        actual = interface._order_factor(side, offset, 0) * probability
        np.testing.assert_allclose(actual, expected, atol=2e-10)


@pytest.mark.parametrize("kind", ["surface", "interface"])
def test_order_fit_reuses_dwba_reference_and_survives_saved_models(kind, tmp_path):
    """Fitted order changes contrast, preserving optics and source coverage."""
    if kind == "surface":
        crystal, component = poisson_crystal(PoissonProfile(0.7, alpha=0.5), w_base=16)
        film = crystal.uc_surface_list[0]
        film.basis_0[:] = film.basis
    else:
        crystal, component = _interface()
    crystal.setEnergy(10000.0)
    profile = component.setCrystallinityProfile(
        CrystallinityProfile(PoissonProfile(2.0, alpha=0.5), minimum=0.2)
    )
    profile.addFitParameter("minimum", name="skin_order", limits=(0, 1))
    crystal = copy.deepcopy(crystal)
    component = crystal.uc_surface_list[-1]
    prepared = _prepare(crystal)
    initial = crystal.dwba.evaluate_prepared(prepared)
    coverage = crystal.optical_profile().copy()
    np.testing.assert_allclose(crystal.getInitialParameters(), [0.2])
    crystal.setParameters([0.7])
    with crystal.dwba.batch():
        changed = crystal.dwba.evaluate_prepared(prepared)
    assert not np.allclose(changed.F_h, initial.F_h)
    np.testing.assert_array_equal(changed.F_reference, initial.F_reference)
    np.testing.assert_array_equal(
        changed.unperturbed_amplitude, initial.unperturbed_amplitude
    )
    np.testing.assert_allclose(
        crystal.optical_profile(), coverage, rtol=2e-14, atol=1e-18
    )
    for clone in (copy.deepcopy(crystal), pickle.loads(pickle.dumps(crystal))):
        np.testing.assert_allclose(clone.F(*_hkl()), crystal.F(*_hkl()))
    path = str(tmp_path / "order.xtal")
    crystal.toFile(path)
    restored = CTRcalc.SXRDCrystal.fromFile(path)
    restored.setEnergy(10000.0)
    np.testing.assert_allclose(restored.F(*_hkl()), crystal.F(*_hkl()), atol=1e-10)
    np.testing.assert_allclose(
        restored.dwba.evaluate_prepared(_prepare(restored)).F_h, changed.F_h, atol=1e-10
    )


def test_no_disorder_profile_recovers_legacy_amplitudes():
    """Enabling independent order at unity preserves both scattering paths."""
    crystal, surface = poisson_crystal(PoissonProfile(1.2), w_base=16)
    crystal.setEnergy(10000.0)
    before = crystal.dwba.evaluate_prepared(_prepare(crystal))
    kinematic = crystal.F(*_hkl())
    surface.setCrystallinityProfile(PoissonProfile(0))
    after = crystal.dwba.evaluate_prepared(_prepare(crystal))
    np.testing.assert_allclose(after.F_h, before.F_h, atol=1e-10)
    np.testing.assert_allclose(after.F_reference, before.F_reference, atol=1e-10)
    np.testing.assert_allclose(crystal.F(*_hkl()), kinematic, atol=1e-10)


def test_ordered_flat_states_match_explicit_terraces():
    """Coherent and incoherent consumers see the same local order on each terrace."""
    crystal, surface = poisson_crystal(PoissonProfile(1.2, alpha=0.7), w_base=16)
    crystal.setEnergy(10000.0)
    order = surface.setCrystallinityProfile(PoissonProfile(3, alpha=0.5))
    crystal.apply_stacking()
    states = surface.flat_domain_corrections(*_hkl())
    baseline = crystal.uc_bulk.F_bulk(*_hkl(), crystal.atten) + crystal.uc_surface_list[
        0
    ].F_uc(*_hkl())
    for state in states.iter_states():
        np.testing.assert_allclose(
            baseline + state.amplitude,
            _flat_ordered_amplitude(state.layer_number, order),
            atol=1e-10,
        )


def test_interface_rejects_truncated_disorder():
    """A long disorder profile must not silently stop at a finite support edge."""
    crystal, interface = _interface()
    crystal.setEnergy(10000.0)
    interface.setCrystallinityProfile(PoissonProfile(100, alpha=0))
    with pytest.raises(ValueError, match="increase EpitaxyInterface fixed_ucs"):
        crystal.F(*_hkl())


def test_skellam_surface_order_and_detachment():
    """A skewed recovery distribution agrees with terraces and detaches cleanly."""
    occupancy = PoissonProfile(-1.2, alpha=0.6)
    crystal, surface = poisson_crystal(occupancy, w_base=16)
    crystal.setEnergy(10000.0)
    original = crystal.F(*_hkl())
    supplied = CrystallinityProfile(SkellamProfile(0.6, 0.3), offset=2, minimum=0.1)
    owned = surface.setCrystallinityProfile(supplied)
    supplied.basis[-1] = 0.9
    assert owned.basis[-1] == 0.1
    heights, masses = height_states(occupancy)
    expected = sum(
        mass * _flat_ordered_amplitude(height, owned)
        for height, mass in zip(heights, masses)
    )
    np.testing.assert_allclose(crystal.F(*_hkl()), expected, atol=2e-8)
    surface.setCrystallinityProfile(None)
    np.testing.assert_allclose(crystal.F(*_hkl()), original, atol=1e-10)


def test_both_interface_profiles_keep_independent_fit_settings(tmp_path):
    """Each interface side retains its own family, values and uncertainty on save."""
    from silx.io.dictdump import dicttoh5

    crystal, interface = _interface()
    crystal.setEnergy(10000.0)
    bottom = interface.setCrystallinityProfile(
        CrystallinityProfile(SkellamProfile(0.4), offset=2, minimum=0.3), side="bottom"
    )
    top = interface.setCrystallinityProfile(
        CrystallinityProfile(PoissonProfile(1, alpha=0.5), minimum=0.2), side="top"
    )
    top.addFitParameter("minimum", name="top_order", limits=(0, 1))
    bottom.addFitParameter("minimum", name="bottom_order", limits=(0, 1))
    np.testing.assert_allclose(crystal.getInitialParameters(), [0.2, 0.3])
    crystal.setParameters([0.4, 0.6])
    interface.setFitErrors([0.01, 0.02])
    optical = crystal.optical_profile().copy()
    expected = crystal.F(*_hkl())
    path = str(tmp_path / "two_sides.xtal")
    crystal.toFile(path)
    dicttoh5(crystal.parametersToDict(), str(tmp_path / "two_sides.h5"))
    restored = CTRcalc.SXRDCrystal.fromFile(str(tmp_path / "two_sides"))
    restored.setEnergy(10000.0)
    saved = restored.uc_surface_list[-1]
    np.testing.assert_allclose(restored.getInitialParameters(), [0.4, 0.6])
    np.testing.assert_allclose(saved.getFitErrors(), [0.01, 0.02])
    np.testing.assert_allclose(restored.F(*_hkl()), expected)
    np.testing.assert_allclose(restored.optical_profile(), optical)

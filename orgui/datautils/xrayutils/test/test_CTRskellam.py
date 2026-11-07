"""Surface, persistence, and scattering regressions for recentered Skellam."""

import copy
from dataclasses import replace
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from .. import CTRcalc, CTRfilm, CTRincoherent
from ..CTRdistributions import PoissonProfile, RecenteredSkellamProfile
from ..CTRutil import generate_surface_termination_cells
from . import _poisson_oracle as oracle


def _hkl():
    return np.array([0., 0., 1., 1.]), np.zeros(4), np.array([.37, .83, .63, 1.7])


def _crystal(profile=None, n_layers=2, bank=False, width=80):
    if profile is None:
        profile = RecenteredSkellamProfile(1.3, .2, .7)
    source = oracle.layered_cell(n_layers, "film")
    film = CTRfilm.Film(source, width_layers=width, name="film")
    template = oracle.layered_cell(n_layers, "surface")
    if bank:
        template = generate_surface_termination_cells(template, template.layer_cycle)
        for cell in template.values():
            cell.basis[:, 4] = .3
            cell.basis_0[:] = cell.basis
    cls = (CTRfilm.PoissonSurface if isinstance(profile, PoissonProfile)
           else CTRfilm.SkellamSurface)
    surface = cls(template, profile=profile, name="surface")
    crystal = CTRcalc.SXRDCrystal(
        oracle.layered_cell(n_layers, "bulk"), film, surface,
        stacking=np.array([1, 2]))
    crystal.setEnergy(10000)
    crystal.apply_stacking()
    return crystal, surface


@pytest.mark.parametrize("n_layers", [1, 2, 3])
def test_coherent_and_incoherent_against_independent_films(n_layers):
    """Complete amplitudes agree with explicitly stacked flat Film states."""
    crystal, surface = _crystal(n_layers=n_layers)
    states = surface.flat_domain_corrections(*_hkl())
    common = oracle.flat_height_crystal(-1, w_base=80, n_layers=n_layers)
    common.setEnergy(10000)
    base = common.F(*_hkl())
    expected_f = np.zeros(4, dtype=complex)
    expected_f2 = np.zeros(4)
    for state in states.iter_states():
        flat = oracle.flat_height_crystal(
            state.layer_number, w_base=80, n_layers=n_layers)
        flat.setEnergy(10000)
        amplitude = flat.F(*_hkl())
        assert_allclose(base + state.amplitude, amplitude, atol=2e-10, rtol=2e-11)
        expected_f += state.probability * amplitude
        expected_f2 += state.probability * abs(amplitude)**2
    assert_allclose(crystal.F(*_hkl()), expected_f, atol=2e-10, rtol=2e-11)
    assert abs(states.probabilities @ (states.layer_numbers + 1)) < 1e-10
    for fraction in (0, .35, 1):
        wrapper = CTRincoherent.SurfaceHeightDomains(
            crystal, surface="surface", incoherent_fraction=fraction)
        expected = (1 - fraction) * abs(expected_f)**2 + fraction * expected_f2
        assert_allclose(wrapper.F2(*_hkl()), expected, atol=2e-9, rtol=2e-11)
        restored = CTRincoherent.create_incoherent_model_from_config(
            crystal, wrapper.to_config())
        assert_allclose(restored.F2(*_hkl()), expected, atol=2e-9, rtol=2e-11)
    with pytest.raises(ValueError, match="PoissonSurface"):
        CTRincoherent.PoissonHeightDomains(crystal, surface="surface")
    wrapper.validate(forward_model="dwba")


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("bank", [False, True])
def test_one_sided_surface_equals_compensated_poisson(sign, bank):
    """Film replacement and termination banks recover both Poisson limits."""
    profile = RecenteredSkellamProfile(1, .7 if sign > 0 else 0,
                                     .7 if sign < 0 else 0)
    reference = PoissonProfile(sign * .7, 1, -sign * .7,
                               tail_probability=1e-14)
    actual, _ = _crystal(profile, bank=bank)
    expected, _ = _crystal(reference, bank=bank)
    assert_allclose(actual.F(*_hkl()), expected.F(*_hkl()), atol=2e-10, rtol=2e-10)
    for crystal in (actual, expected):
        crystal.dwba.set_ctr_geometry(alpha_i=.05)
    result = actual.dwba.evaluate(*_hkl())
    reference_result = expected.dwba.evaluate(*_hkl())
    # Different retained tails place the incident-wave reference plane at
    # different heights in vacuum. Compare complex amplitudes in one gauge.
    prepared, reference_prepared = result.prepared, reference_result.prepared
    dz_i = (prepared.field_i.z_reference[0]
            - reference_prepared.field_i.z_reference[0])
    dz_f = (prepared.field_f.z_reference[0]
            - reference_prepared.field_f.z_reference[0])
    phase = np.exp(1j * prepared.k0 * (
        np.sin(prepared.alpha_i) * dz_i + np.sin(prepared.alpha_f) * dz_f))
    assert_allclose(result.F_atomic * phase, reference_result.F_atomic,
                    atol=2e-9, rtol=2e-8)
    # On specular rods the planar-reference split also depends on the slab
    # support. Its change cancels against r0 in the physical total amplitude.
    assert_allclose(result.total_amplitude * phase,
                    reference_result.total_amplitude, atol=2e-12, rtol=2e-8)
    assert_allclose(result.reflectivity, reference_result.reflectivity, rtol=2e-8)


def test_flat_limit_mean_height_registry_and_inventory():
    """Zero time is sharp; equal-layer inventory and mean registry stay fixed."""
    crystal, surface = _crystal(RecenteredSkellamProfile(0, 1, 2), n_layers=1)
    flat = oracle.flat_height_crystal(-1, w_base=80, n_layers=1)
    flat.setEnergy(10000)
    assert_allclose(crystal.F(*_hkl()), flat.F(*_hkl()), atol=1e-11)
    for time in (.3, 2, 3):
        surface.basis[0] = time
        crystal.apply_stacking()
        assert surface.mean_height_absolute == surface.below_H
        cells = [*surface.film_layer_ucs, *surface.layer_ucs,
                 *surface._film_termination_ucs.values()]
        atoms = sum(c.basis[:, 6].sum() * sum(c.coherentDomainOccupancy)
                    for c in cells)
        assert abs(atoms) < 1e-10
    # An oblique stacking repeat must use the recentered mean as well.
    surface._below_registry = (1.2, -.4, 0.)
    assert_allclose(surface.layer_state.registry_cartesian, surface._below_registry)


@pytest.mark.parametrize("extension", ["xtal", "xpr"])
def test_native_files_fit_settings_and_copy(tmp_path, extension):
    """Native files and optimizer copies preserve live rates and tolerances."""
    profile = RecenteredSkellamProfile(
        1.234567891, .223456789, .734567891,
        tail_probability=1e-13, moment_tolerance=1e-12)
    crystal, surface = _crystal(profile, bank=True)
    surface.addFitParameter("k_up", limits=(0, 1), name="up_rate")
    assert_allclose(crystal.getInitialParameters(), [.223456789])
    before = crystal.F(*_hkl())
    crystal.setParameters([.5])
    after = crystal.F(*_hkl())
    assert not np.allclose(before, after)
    assert surface.mean_height_absolute == surface.below_H
    path = str(tmp_path / f"skellam.{extension}")
    crystal.toFile(path)
    restored = CTRcalc.SXRDCrystal.fromFile(path)
    clones = restored, copy.deepcopy(crystal), pickle.loads(pickle.dumps(crystal))
    for clone in clones:
        clone.setEnergy(10000)
        assert_allclose(clone.F(*_hkl()), after, atol=2e-9)
        target = clone.uc_surface_list[-1]
        assert isinstance(target, CTRfilm.SkellamSurface)
        assert_array_equal(target.basis, surface.basis)
        assert target._profile_settings() == surface._profile_settings()
    settings = crystal.parametersToDict()
    restored.parametersFromDict(settings)
    assert_allclose(restored.getInitialParameters(), [.5])
    restored.setParameters([.3])
    crystal.setParameters([.3])
    assert_allclose(restored.F(*_hkl()), crystal.F(*_hkl()), atol=2e-9)


def test_water_order_and_dwba_cache_invalidation():
    """Water follows terraces; order preserves optics; rate edits refresh DWBA."""
    from .test_CTRwater import _water

    crystal, surface = _crystal()
    water = _water("step")
    surface.attach_water(water)
    wet = crystal.F(*_hkl())
    heights, masses, tail = surface.terrace_distribution()
    assert_allclose(masses @ heights, surface.below_H, atol=1e-9)
    assert tail <= surface.profile.tail_probability
    optics = crystal.optical_profile().copy()
    surface.setCrystallinityProfile(PoissonProfile(2, alpha=.4))
    assert not np.allclose(crystal.F(*_hkl()), wet)
    assert_allclose(crystal.optical_profile(), optics, atol=1e-17)
    crystal.dwba.set_ctr_geometry(alpha_i=.05)
    result = crystal.dwba.evaluate(*_hkl())
    surface.basis[1] += .2
    changed = crystal.dwba.evaluate(*_hkl())
    assert np.isfinite(changed.F_h).all()
    assert not np.allclose(result.F_h, changed.F_h)
    restored = CTRcalc.SXRDCrystal.fromStr(crystal.toStr())
    restored.setEnergy(10000)
    assert_allclose(restored.F(*_hkl()), crystal.F(*_hkl()), atol=2e-9)
    # Tolerances are part of geometry invalidation even with unchanged basis.
    surface.profile = replace(surface.profile, moment_tolerance=1e-15)
    assert surface._repeat_state_changed()
    crystal.F(*_hkl())
    assert not surface._repeat_state_changed()


def test_rejects_wrong_profiles_invalid_rates_and_film_overlap():
    """Invalid candidates must fail explicitly without clipping or shifting."""
    with pytest.raises(TypeError, match="RecenteredSkellamProfile"):
        CTRfilm.SkellamSurface(oracle.layered_cell(), profile=PoissonProfile(1))
    with pytest.raises(ValueError, match="below the underlying Film"):
        _crystal(width=1)
    crystal, surface = _crystal()
    surface.basis[1] = -1
    with pytest.raises(ValueError, match="nonnegative"):
        crystal.F(*_hkl())


def test_optimizer_copies_rates_and_recovers_after_invalid_candidate():
    """Coherent and domain models use the fitted profile in optimizer copies."""
    from .. import CTRopt, CTRplotutil

    crystal, surface = _crystal()
    surface.addFitParameter("k_up", name="up_rate", limits=(0, 1))
    ell = np.array([.4, .8, 1.4])
    measured = abs(crystal.F(np.ones(3), np.zeros(3), ell))
    rod = CTRplotutil.CTR((1., 0.), ell, measured, np.ones(3))
    data = CTRplotutil.CTRCollection([rod])
    for model in (crystal, CTRincoherent.SurfaceHeightDomains(
            crystal, surface="surface", incoherent_fraction=.3)):
        optimizer = CTRopt.CTROptimizer(model, data, scale_policy={"F": "fixed"})
        optimizer.prepareFit()
        initial = optimizer.flat_Fcalc([.2]).copy()
        changed = optimizer.flat_Fcalc([.5])
        assert not np.allclose(initial, changed)
        with pytest.raises(ValueError, match="nonnegative"):
            optimizer.flat_Fcalc([-1])
        assert_allclose(optimizer.flat_Fcalc([.2]), initial, atol=1e-12)
    assert surface.basis[1] == .2

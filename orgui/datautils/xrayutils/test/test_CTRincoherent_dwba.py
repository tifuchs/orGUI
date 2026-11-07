"""Complete terrace DWBA limits, physical oracles, and public fitting path."""

import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from .. import CTRincoherent, CTRopt, CTRplotutil, CTRresolution
from ..CTRdistributions import PoissonProfile, RecenteredSkellamProfile
from .._CTRincoherent_dwba import (
    _flat_crystal, _phase_to, _shared_preparation,
)
from . import _poisson_oracle as oracle
from .test_CTRopt_dwba import _ctr, _generated_angles
from .test_CTRskellam import _crystal as skellam_crystal


def _model(profile=None, policy="shared", fraction=.205):
    if profile is None:
        profile = PoissonProfile(0, offset=.35)
    if isinstance(profile, PoissonProfile):
        crystal, _ = oracle.poisson_crystal(profile, w_base=40)
        crystal.setEnergy(10000)
        crystal.apply_stacking()
    else:
        crystal, _ = skellam_crystal(profile, width=40)
    crystal.dwba.set_ctr_geometry(alpha_i=.003)
    return CTRincoherent.SurfaceHeightDomains(
        crystal, surface="surface", incoherent_fraction=fraction,
        optical_reference=policy)


def _hkl():
    # Specular critical region plus non-specular critical incidence,
    # anti-Bragg/interference minima and both sides of a Bragg feature.
    return (np.array([0, 0, 1, 1, 1, 1, 1.]), np.zeros(7),
            np.array([.01, .03, .06, .45, 1.7, 1.999, 2.001]))


@pytest.mark.parametrize("profile", [PoissonProfile(-.4),
                                    RecenteredSkellamProfile(.3, .2, .5)])
def test_shared_coherent_endpoint_and_probability_weighted_states(profile):
    """Shared fields preserve native coherence and mix complete intensities."""
    model = _model(profile, fraction=0)
    expected = model.coherent_model.dwba.evaluate(*_hkl()).reflectivity
    assert_array_equal(model.dwba.evaluate(*_hkl()).reflectivity, expected)
    model.incoherent_fraction = 1
    large = model.dwba.evaluate(*_hkl())
    model.incoherent_fraction = .205
    partial = model.dwba.evaluate(*_hkl())
    assert_allclose(partial.reflectivity, .795 * expected + .205 * large.reflectivity)
    assert_allclose(large.reflectivity, large.incoherent_intensity)
    assert np.all(np.isfinite(partial.reflectivity))


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("fraction", [0, .205, 1])
def test_single_height_limit(policy, fraction):
    """A deterministic surface is independent of coherence and reference policy."""
    model = _model(PoissonProfile(0, offset=-2), policy, fraction)
    native = model.coherent_model.dwba.evaluate(*_hkl())
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity,
                    native.reflectivity, rtol=3e-12, atol=1e-18)


@pytest.mark.parametrize("channels", [("s", "s"), ("s", "p"),
                                     ("p", "s"), ("p", "p")])
def test_state_specific_against_independent_flat_films(channels):
    """Independent complete films verify probabilities, channels and phase gauge."""
    model = _model(policy="state_specific")
    crystal = model.coherent_model
    options = dict(polarization_i=channels[0], polarization_f=channels[1])
    reference = crystal.dwba.evaluate(*_hkl(), **options).prepared
    states = model.target_surface().flat_domain_corrections(*_hkl())
    amplitude = np.zeros(7, complex)
    intensity = np.zeros(7)
    for height, probability in zip(states.layer_numbers, states.probabilities):
        flat = oracle.flat_height_crystal(height, w_base=40)
        flat.setEnergy(10000)
        flat.dwba.set_ctr_geometry(alpha_i=.003)
        result = flat.dwba.evaluate(*_hkl(), **options)
        field = result.total_amplitude * _phase_to(result.prepared, reference)
        amplitude += probability * field
        intensity += probability * abs(field)**2
    result = model.dwba.evaluate(*_hkl(), **options)
    assert_allclose(result.incoherent_intensity, intensity, rtol=3e-9, atol=1e-17)
    assert_allclose(result.coherent_intensity, abs(amplitude)**2,
                    rtol=3e-9, atol=1e-17)
    assert_allclose(result.reflectivity, .795 * abs(amplitude)**2 + .205 * intensity,
                    rtol=3e-9, atol=1e-17)


def test_shared_fields_against_independent_flat_films():
    """Flat-film atomic sums on shared fields retain full state interference."""
    model = _model(PoissonProfile(-.4), fraction=1)
    native = model.coherent_model.dwba.evaluate(*_hkl())
    states = model.target_surface().flat_domain_corrections(*_hkl())
    amplitude = np.zeros(7, complex)
    intensity = np.zeros(7)
    for height, probability in zip(states.layer_numbers, states.probabilities):
        flat = oracle.flat_height_crystal(height, w_base=40)
        flat.setEnergy(10000)
        _, atomic, profile = flat.dwba._atomic_model_snapshot()
        prepared = _shared_preparation(flat.dwba, native.prepared, atomic, profile)
        result = flat.dwba.evaluate_prepared(prepared)
        field = native.unperturbed_amplitude + native.amplitude_prefactor * (
            result.F_atomic * _phase_to(prepared, native.prepared)
            - native.F_reference)
        amplitude += probability * field
        intensity += probability * abs(field)**2
    assert_allclose(amplitude, native.total_amplitude, rtol=2e-9, atol=1e-12)
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity,
                    intensity, rtol=2e-9, atol=1e-17)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("channels", [("s", "s"), ("s", "p"),
                                     ("p", "s"), ("p", "p")])
def test_weak_optics_born_limit_matches_kinematic_ensemble(policy, channels):
    """Dilute density recovers F2 with identical attenuation and vector factors."""
    model = _model(policy=policy)
    crystal = model.coherent_model
    for cell in (crystal.uc_bulk, crystal.uc_surface_list[0].unitcell,
                 model.target_surface().unitcell):
        cell.basis[:, 6] *= 1e-7
        cell.basis_0[:, 6] *= 1e-7
    model.target_surface().setFitParameters([])
    crystal.atten = .03
    crystal.dwba.set_ctr_geometry(alpha_i=.12)
    h, k, ell = np.ones(3), np.zeros(3), np.array([.7, 1.4, 2.03])
    result = model.dwba.evaluate(h, k, ell, bulk_attenuation=crystal.atten,
                                 polarization_i=channels[0],
                                 polarization_f=channels[1])
    p = result.prepared
    ai, af, ca, sa = p.alpha_i, p.alpha_f, p.cos_azimuth, p.sin_azimuth
    contraction = {
        ("s", "s"): ca,
        ("s", "p"): -np.sin(af) * sa,
        ("p", "s"): np.sin(ai) * sa,
        ("p", "p"): np.cos(ai)*np.cos(af) - np.sin(ai)*np.sin(af)*ca,
    }[channels]
    expected = model.F2(h, k, ell) * contraction**2
    actual = result.reflectivity / abs(result.amplitude_prefactor)**2
    assert_allclose(actual, expected, rtol=2e-6, atol=1e-24)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_channels_domains_batching_and_config_round_trip(policy):
    """Outer and cell domains survive serialization and arbitrary point batches."""
    model = _model(policy=policy)
    crystal = model.coherent_model
    rotation = np.diag([-1., -1., 1.])
    crystal.setDomain(0, [(np.eye(3), .6), (rotation, .4)])
    crystal.setDomain(1, [(np.eye(3), .6), (rotation, .4)])
    options = dict(s_fraction=.63, outgoing="unanalysed")
    expected = model.dwba.polarization_reflectivity(*_hkl(), **options)
    explicit = sum(weight * model.dwba.evaluate(
        *_hkl(), polarization_i=pi, polarization_f=pf).reflectivity
        for pi, pf, weight in crystal.dwba._polarization_pairs(.63, "unanalysed"))
    assert_allclose(expected, explicit, rtol=2e-13)
    restored = CTRincoherent.create_incoherent_model_from_config(
        copy.deepcopy(crystal), model.to_config())
    assert restored.optical_reference == policy
    pieces = []
    with restored.dwba.batch():
        for indices in (slice(0, 2), slice(2, 5), slice(5, 7)):
            pieces.append(restored.dwba.polarization_reflectivity(
                *(v[indices] for v in _hkl()), **options))
    assert_allclose(np.concatenate(pieces), expected, rtol=2e-12)
    scalar = restored.dwba.evaluate(1., 0., .45).reflectivity
    assert isinstance(scalar, float)


@pytest.mark.parametrize("change", ["geometry", "domains", "roughness",
                                  "crystallinity", "energy", "occupancy"])
@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_live_mutations_equal_fresh_evaluation(change, policy):
    """No terrace snapshot outlives geometry, structural or optical changes."""
    model = _model(policy=policy)
    facade = model.dwba
    before = facade.evaluate(*_hkl()).reflectivity
    crystal = model.coherent_model
    if change == "geometry":
        crystal.dwba.set_ctr_geometry(alpha_i=.007)
    elif change == "domains":
        crystal.setDomain(1, [(np.eye(3), .7)])
    elif change == "roughness":
        model.target_surface().basis[2] = .7
    elif change == "crystallinity":
        model.target_surface().unitcell.crystallinity = .3
        model.target_surface().setFitParameters([])
    elif change == "energy":
        crystal.setEnergy(11000)
    else:
        model.target_surface().unitcell.basis[:, 6] *= .5
        model.target_surface().unitcell.basis_0[:, 6] *= .5
        model.target_surface().setFitParameters([])
    actual = facade.evaluate(*_hkl()).reflectivity
    fresh = copy.deepcopy(model)
    fresh.coherent_model.dwba.clear_cache()
    assert_allclose(actual, fresh.dwba.evaluate(*_hkl()).reflectivity, rtol=2e-12)
    assert not np.allclose(actual, before, rtol=1e-7, atol=1e-20)


@pytest.mark.parametrize("calculation", ["convolve", "sample"])
@pytest.mark.parametrize("quantity", ["reflectivity", "structure_factor"])
def test_optimizer_convolves_mixed_field_intensity_before_conversion(
    calculation, quantity
):
    """The public fit path reduces the convolved ensemble field intensity."""
    model = _model()
    ctr = _ctr(model.coherent_model, quantity=quantity, hk=(1, 0),
               s_fraction=.6, outgoing="unanalysed",
               polarization_factor=(np.full(3, .8)
                                    if quantity == "structure_factor" else None),
               geometry=CTRplotutil.CTRScanGeometry("in", .003),
               l_values=[.4, .5, .6])
    optimizer = CTRopt.CTROptimizer(model, CTRplotutil.CTRCollection([ctr]))
    optimizer.set_dwba()
    resolution = CTRresolution.GaussianResolution(.015, 0, 0)
    optimizer.set_resolution(resolution, calculation=calculation)
    optimizer.prepareFit()
    angles = _generated_angles(optimizer.xtal, ctr)

    def field_at(h, k, ell):
        """Evaluate field intensity at the resolution quadrature points."""
        local = _generated_angles(optimizer.xtal, ctr, h=h, k=k, l_values=ell)
        return optimizer.model.dwba.polarization_reflectivity_from_vlieg(
            *(local[name] for name in local.dtype.names), s_fraction=.6,
            outgoing="unanalysed")

    if calculation == "sample":
        field = CTRresolution.sample_intensity(
            ctr.harr, ctr.karr, ctr.l, field_at, resolution, angles)
    else:
        field = CTRresolution.fast_convolve_intensity(
            ctr.harr, ctr.karr, ctr.l, field_at(ctr.harr, ctr.karr, ctr.l),
            resolution, angles)
    prefactor = optimizer.xtal.dwba.evaluate_from_vlieg(
        *(angles[name] for name in angles.dtype.names)).amplitude_prefactor
    expected = (field if quantity == "reflectivity" else
                np.sqrt(field / .8) / abs(prefactor))
    assert_allclose(optimizer._dwba_prediction(ctr), expected, rtol=2e-11)


def test_flat_states_preserve_local_crystallinity_and_original_model():
    """Depth-dependent order follows each terrace, without mutating the fit."""
    from ..CTRcrystallinity import CrystallinityProfile
    model = _model()
    model.target_surface().setCrystallinityProfile(
        CrystallinityProfile(PoissonProfile(1.3, alpha=.4)))
    original = model.F2(*_hkl()).copy()
    target = model._bound_target()
    states = target.flat_domain_corrections(*_hkl())
    common = model.coherent_model.evaluate_kinematic(*_hkl())
    base = common.total - common.component("surface").amplitude
    for state in states.iter_states():
        flat = _flat_crystal(model, state.layer_number)
        assert_allclose(flat.F(*_hkl()), base + state.amplitude,
                        rtol=2e-10, atol=2e-9)
    model.dwba.evaluate(*_hkl())
    assert_array_equal(model.F2(*_hkl()), original)


def test_old_configuration_and_invalid_policy():
    """Legacy configs select shared fields; invalid policies fail explicitly."""
    model = _model()
    config = model.to_config()
    del config["settings"]["optical_reference"]
    restored = CTRincoherent.create_incoherent_model_from_config(
        model.coherent_model, config)
    assert restored.optical_reference == "shared"
    with pytest.raises(ValueError, match="optical_reference"):
        model.optical_reference = "ambiguous"
    with pytest.raises(ValueError, match="semi_infinite"):
        model.dwba.evaluate(*_hkl(), bulk_mode="unit_cell")


def test_continuous_water_requires_state_specific_fields():
    """Attached water follows the flat states and uses its scalar phase origin."""
    from .test_CTRwater import _water
    model = _model(policy="state_specific")
    model.target_surface().attach_water(_water("step_tanh"))
    result = model.dwba.evaluate_from_glancing(0, 0, [.01, .06], [.01, .06])
    assert np.all(np.isfinite(result.reflectivity))
    assert np.all(result.reflectivity >= 0)
    model.optical_reference = "shared"
    with pytest.raises(NotImplementedError, match="layered optical reference"):
        model.dwba.evaluate_from_glancing(0, 0, .01, .01)

"""Independent closed-form and no-discretization checks for water DWBA."""

from pathlib import Path

import numpy as np
import pytest

from .._CTRwater import _fourier, _parameters
from .._CTRwater_analytic import (
    _contrast_transform, _epstein_coefficients, _interface_fields, _tanh_ratio,
)
from .test_CTRwater import _water


def _terrace_model(kind, roughness, a0=1.0):
    """Build the 20 keV RuO2/TiO2 case using only installed test resources."""
    from ..CTRcalc import SXRDCrystal
    from ..CTRfilm import PoissonSurface
    from ..CTRuc import WaterModel

    path = Path(__file__).parent / "testdata" / "RuO2_TiO2_Poisson_etching.xtal"
    model = SXRDCrystal.fromFile(path)
    surface = next(c for c in model.uc_surface_list if isinstance(c, PoissonSurface))
    surface.basis[1] = roughness
    if kind != "dry":
        family, reference = kind.split("_")
        cell = model.uc_bulk
        height = cell.volume / cell.uc_area
        water = WaterModel(cell.a.copy(), np.rad2deg(cell.alpha), family,
                           optical_reference=reference, name="water")
        # Physical lengths in Angstrom, converted to fractional normal repeats.
        water.setWaterParameters(3.0/height, 2.8/height, .6/height, .4/height, a0)
        water.adsorption_position = "onset"
        surface.attach_water(water)
    model.setEnergy(20000.0)
    model.apply_stacking()
    return model


@pytest.mark.parametrize("kind", ["step", "layered", "tanh", "1layer_O"])
@pytest.mark.parametrize("sigma", [0., .6])
def test_analytic_contrast_against_existing_fourier(kind, sigma):
    """Check signs, molecular factors and phase against the kinematic sum."""
    water = _water(kind)
    q = np.array([.01, .2, 1., 2.7, -.3, .2+.03j])
    zg = water.reference_position + .21
    expected = (_fourier(water, q)/water.uc_area*np.exp(-1j*q*zg)
                - water.bulk_density*1j/q*_tanh_ratio(q, sigma))
    np.testing.assert_allclose(_contrast_transform(water, q, zg, sigma),
                               expected, rtol=3e-11, atol=2e-13)


@pytest.mark.parametrize("kind", ["step", "layered", "tanh"])
def test_equimolar_zero_channel_and_frozen_shift(kind):
    """Opposite DWBA branches have q=0; no divergent half-spaces are formed."""
    water = _water(kind)
    q = np.array([0, 1e-14, -1e-14, 1e-14j])
    zg = water.reference_position
    np.testing.assert_allclose(_contrast_transform(water, q, zg, .6), 0,
                               atol=1e-13)
    np.testing.assert_allclose(_contrast_transform(water, q, zg+.3, .6),
                               .3*water.bulk_density, atol=1e-13)


def test_hosemann_against_explicit_gaussian_layers():
    """The infinite expression equals explicit layers where the sum converges."""
    from .._CTRwater import _factors
    water = _water("layered")
    z, d, s, b, a = _parameters(water)
    q = np.array([.2+.04j, .6+.04j, 1.4+.04j])
    j = np.arange(4000)[:, None]
    weights = np.ones_like(j, float)
    weights[0] = a
    f = sum(c*np.exp(-v*q*q) for c, v in _factors(water))
    expected = water.pw*d*f*np.sum(
        weights*np.exp(1j*q*(z+j*d)-.5*(s*s+j*b*b)*q*q), axis=0)
    zg = water.reference_position
    expected = expected*np.exp(-1j*q*zg)-water.bulk_density*1j/q
    np.testing.assert_allclose(_contrast_transform(water, q, zg), expected,
                               rtol=2e-12, atol=1e-13)


@pytest.mark.parametrize("kind", ["step_sharp", "layered"])
@pytest.mark.parametrize("polarization", ["s", "p"])
@pytest.mark.parametrize("energy", [1000., 10000.])
def test_four_channels_against_constant_medium_integral(kind, polarization, energy):
    """Integrate vector fields independently, including absorbing water for p."""
    from scipy.integrate import simpson
    from ..CTRcalc import SXRDCrystal
    from .._CTRwater_analytic import _water_contrast
    from .test_CTRcrystallinity import _cell

    water = _water(kind)
    crystal = SXRDCrystal(_cell(), water)
    crystal.setEnergy(energy)
    angles = np.array([.003, .01, .08])  # glancing angles in radians
    prepared = crystal.dwba.prepare_from_glancing(
        0., 0., angles, angles,
        polarization_i=polarization, polarization_f=polarization)
    records = crystal.dwba._live_records()
    actual = _water_contrast(prepared, records[-1], len(records)-1)
    spec, field = prepared.reference.continuous, prepared.field_i
    center = spec.centers[0]
    host = len(field.n)-1-np.searchsorted(spec.boundaries, center, side="right")
    k = field.kz[host]
    plus, minus = field.A_plus[host], field.A_minus[host]
    low, high = water.transition_extent()
    expected = np.zeros(len(angles), complex)
    for left, right, rho0 in ((low, center, 0), (center, high, water.bulk_density)):
        z = np.linspace(left, right, 20001)
        contrast = water.zDensity_G(z, 0., 0.) - rho0
        phase = k*(z[:, None]-field.z_reference[host])
        down, up = plus*np.exp(-1j*phase), minus*np.exp(1j*phase)
        if polarization == "s":
            product = (down + up)**2
        else:
            # Construct Cartesian E components from transversality k.E=0.
            # The reciprocal exit field reverses the tangential direction;
            # use a bilinear product (no complex conjugation) in DWBA.
            tangent = np.sin(angles)*(down - up)
            parallel_k = prepared.k0*field.n[0]*np.cos(angles)
            normal = -parallel_k/k*np.sin(angles)*(down + up)
            product = normal**2 - tangent**2
        expected += simpson(contrast[:, None]*product, x=z, axis=0)
    expected *= prepared.reference_area
    np.testing.assert_allclose(actual, expected, rtol=2e-6, atol=2e-7)


def test_tanh_interface_against_exact_solution_and_sharp_limit():
    """Both directions use exact isolated Epstein amplitudes, including phase."""
    from ..CTRcontinuous import solve_tanh_field
    from ..CTRoptics import solve_electric_field
    values = np.array([[-1., 8e-6, 1e-7], [1., 2e-6, 2e-8]])
    angles = np.array([.001, .003, .004, .02, .1])
    hard = solve_electric_field(values, 20000., angles, boundaries=[0.])
    for sigma in [.6, 2.]:
        exact = solve_tanh_field(hard.n[1], hard.n[0], 0., sigma,
                                 20000., angles, method="analytic")
        fast = _interface_fields(values, np.array([0.]), np.array([sigma]),
                                 20000., angles, "s")
        np.testing.assert_allclose(fast.r_S, exact.r_S, rtol=2e-9, atol=2e-12)
        np.testing.assert_allclose(fast.t_S, exact.t_S, rtol=2e-9, atol=2e-12)
    r, t = _epstein_coefficients(hard.kz[0], hard.kz[1], 1e-8)
    np.testing.assert_allclose(r, hard.r_S, atol=2e-12)
    np.testing.assert_allclose(t, hard.t_S, atol=2e-12)


def test_interface_composition_sharp_limit():
    """Check amplitudes and reference planes against the native transfer solver."""
    from ..CTRoptics import solve_electric_field
    boundaries = np.array([-8., 1., 7.])
    values = np.array([[-9., 8e-6, 1e-7], [-3., 4e-6, 2e-8],
                       [4., 1e-6, 1e-8], [8., 2e-6, 2e-8]])
    angles = np.array([.001, .003, .004, .02, .1])
    hard = solve_electric_field(values, 20000., angles, boundaries=boundaries)
    fast = _interface_fields(values, boundaries, np.full(3, 1e-8),
                             20000., angles, "s")
    for name in ("r_S", "t_S", "A_plus", "A_minus"):
        np.testing.assert_allclose(getattr(fast, name), getattr(hard, name),
                                   rtol=1e-10, atol=1e-10)


def test_critical_tanh_coefficients_are_regular():
    """Gamma recurrences keep the lossless critical endpoint finite."""
    incident = np.array([-.03+0j])
    r, t = _epstein_coefficients(incident, np.array([0j]), .6)
    assert np.all(np.isfinite(r)) and np.all(np.isfinite(t))
    np.testing.assert_allclose(abs(r), 1., atol=2e-15)
    for transmitted in [-1e-14+0j, 1e-14j]:
        rn, tn = _epstein_coefficients(incident, np.array([transmitted]), .6)
        np.testing.assert_allclose(rn, r, atol=2e-12)
        np.testing.assert_allclose(tn, t, atol=2e-12)


def test_lossless_tanh_stack_conserves_flux():
    """Check reverse coefficients and propagation through several interfaces."""
    values = np.array([[-9., 8e-6, 0.], [-3., 4e-6, 0.],
                       [4., 1e-6, 0.], [8., 2e-6, 0.]])
    field = _interface_fields(
        values, np.array([-8., 1., 7.]), np.array([1.1, .6, .2]),
        20000., np.array([.001, .003, .006, .04, .1]), "s")
    flux = abs(field.r_S)**2 + abs(field.t_S)**2*field.kz[-1].real/field.kz[0].real
    np.testing.assert_allclose(flux, 1., rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize("kind", ["step_tanh", "step_sharp",
                                 "layered_sharp", "layered_tanh"])
@pytest.mark.parametrize("roughness", [0., .15])
def test_production_never_samples_density_or_slices(kind, roughness, monkeypatch):
    """Adsorption/layering cannot grow the optical stack or trigger quadrature."""
    from .. import _CTRwater, _CTRwater_continuous, _CTRwater_dwba
    from ..CTRuc import WaterModel
    from ..CTRoptics import LayeredElectricField
    import scipy.integrate

    def forbidden(*args, **kwargs):
        raise AssertionError("Production water attempted real-space integration")
    for module in (_CTRwater, _CTRwater_continuous, _CTRwater_dwba):
        monkeypatch.setattr(module, "_integration_breaks", forbidden)
    monkeypatch.setattr(_CTRwater_continuous, "_evaluate_continuous", forbidden)
    monkeypatch.setattr(scipy.integrate, "solve_ivp", forbidden)
    monkeypatch.setattr(WaterModel, "zDensity_G", forbidden)
    monkeypatch.setattr(WaterModel, "transition_extent", forbidden)
    model = _terrace_model(kind, roughness)
    surface = next(c for c in model.uc_surface_list if getattr(c, "water", None))
    water = surface.water
    prepared = model.dwba.prepare_from_glancing(
        0., 0., [.003, .01, .1], [.003, .01, .1])
    assert isinstance(prepared.field_i, LayeredElectricField)
    count = len(prepared.reference.n)
    dry = _terrace_model("dry", roughness).dwba.prepare_from_glancing(
        0., 0., [.003, .01, .1], [.003, .01, .1])
    # Only physical terrace onsets are added to the existing dry reference.
    assert count <= len(dry.reference.n) + len(water._placements[0])
    result = model.dwba.evaluate_prepared(prepared)
    assert np.all(np.isfinite(result.total_amplitude))
    water.basis[2] *= .01
    water.basis[3] *= .01
    changed = model.dwba.evaluate_prepared(prepared)
    assert np.all(np.isfinite(changed.total_amplitude))
    assert len(prepared.reference.n) == count
    assert not np.array_equal(changed.F_h, result.F_h)


@pytest.mark.parametrize("roughness", [0., .15])
@pytest.mark.parametrize("a0", [0., 1.])
def test_full_ctrs_tanh_against_erf_confirmation(roughness, a0):
    """Bound the fitting approximation for 00L/01L/11L, through L=3 r.l.u.

    This is a model-specific accuracy bound, not an identity between optical
    references. The independent numerical erf path includes the full physical
    water profile. Cover the critical region, fringes and near-Bragg points;
    allow 0.5% intensity error with a 1e-8 peak floor at deep minima.
    """
    models = [_terrace_model("step_" + ref, roughness, a0)
              for ref in ("tanh", "erf")]
    ell = np.array([.008, .015, .025, .04, .07, .13, .33, .71,
                    1.15, 1.64, 1.93, 2.07, 2.32, 2.71, 3.0])
    for h, k in ((0, 0), (0, 1), (1, 1)):
        amplitudes = []
        for model in models:
            model.dwba.set_ctr_geometry(equal_angles=True, rods=[(h, k)])
            result = model.dwba.evaluate(np.full_like(ell, h),
                                         np.full_like(ell, k), ell)
            prepared = result.prepared
            origin = float(np.atleast_1d(prepared.field_i.z_reference)[0])
            # Normalize both incident/reciprocal waves at z=0 Angstrom.
            # In an absorbing ambient a phase shift also changes the modulus.
            phase = np.exp(1j*prepared.k0*prepared.reference.n[0]*origin
                           * (np.sin(prepared.alpha_i)+np.sin(prepared.alpha_f)))
            amplitudes.append(result.F_effective*phase)
            water = next(c.F_h for c in result.contributions if c.role == "water")
            assert np.all(np.isfinite(result.total_amplitude))
            if h == k == 0 and model is models[0]:
                assert np.max(abs(water)) > 1e-3
            elif h != 0 or k != 0:
                np.testing.assert_array_equal(water, 0.)
        actual, expected = (abs(a)**2 for a in amplitudes)
        np.testing.assert_allclose(actual, expected, rtol=.005,
                                   atol=1e-8*expected.max(), err_msg=f"{h}{k}L")


def test_tanh_p_requires_explicit_reference_choice():
    """Do not silently replace an analytic calculation with an ODE for p."""
    from ..CTRcalc import SXRDCrystal
    from .test_CTRcrystallinity import _cell
    model = SXRDCrystal(_cell(), _water("step_tanh"))
    with pytest.raises(ValueError, match="s polarization only"):
        model.dwba.prepare_from_glancing(0., 0., .01, .01, polarization_i="p")

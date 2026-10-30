"""Independent density/transform and surface-following water checks."""

import copy

import numpy as np
import pytest
from scipy.integrate import simpson

from ..CTRuc import WaterModel
from .._CTRwater import _fourier


def _water(kind="step"):
    options = {}
    if kind.startswith("step_"):
        options["optical_reference"] = kind.removeprefix("step_")
        kind = "step"
    water = WaterModel(
        [3.0, 3.0, 4.0], [90.0, 90.0, 90.0], kind, name="water", **options
    )
    water.setEnergy(10000.0)
    water.setWaterParameters(0.2, 0.7, 0.15, 0.25, 0.8)
    return water


@pytest.mark.parametrize("kind", WaterModel.supportedTypes + WaterModel.legacyTypes)
def test_density_fourier_pair(kind):
    """Transform density derivative to avoid the divergent bulk half-space."""
    water = _water(kind)
    z = np.linspace(-30.0, 150.0, 36001)
    density = water.zDensity_G(z, 0.0, 0.0)
    derivative = np.gradient(density, z, edge_order=2)
    q = np.array([0.17, 0.43, 1.1])
    actual = (
        -water.uc_area
        * simpson(derivative[:, None] * np.exp(1j * z[:, None] * q), x=z, axis=0)
        / (1j * q)
    )
    np.testing.assert_allclose(actual, _fourier(water, q), rtol=1e-5, atol=2e-6)


@pytest.mark.parametrize("kind", WaterModel.supportedTypes + WaterModel.legacyTypes)
def test_reference_conserves_integrated_density(kind):
    """Equimolar reference uses all Gaussian layers and their broadening."""
    water = _water(kind)
    z = np.linspace(-30.0, 160.0, 38001)
    zg = water.reference_position
    # Integrate actual density; the reference step integral is exact.
    excess = simpson(water.zDensity_G(z, 0.0, 0.0), x=z) - water.bulk_density * (
        z[-1] - zg
    )
    assert abs(excess) < 2e-5


@pytest.mark.parametrize("kind", WaterModel.supportedTypes + WaterModel.legacyTypes)
def test_scalar_broadcast_and_text_roundtrip(kind):
    """Profiles preserve serializable parameters and physical translations."""
    water = _water(kind)
    text = water.toStr().removeprefix("return\n")
    restored = WaterModel.fromStr(text)
    restored.setEnergy(10000.0)
    np.testing.assert_allclose(restored.F_uc(0.0, 0.0, 0.3), water.F_uc(0.0, 0.0, 0.3))
    assert water.F_uc(1.0, 0.0, 0.0) == 0
    water.pos_absolute = 2.0
    q = 0.3 * 2 * np.pi / 4
    np.testing.assert_allclose(
        water.F_uc(0.0, 0.0, 0.3), restored.F_uc(0.0, 0.0, 0.3) * np.exp(2j * q)
    )
    with pytest.raises(ValueError, match="q_z=0"):
        water.F_uc(0.0, 0.0, 0.0)


def test_zero_excess_single_layer_and_nondecaying_rejection():
    """Zero excess must not evaluate log(0); periodic water is not an ambient."""
    water = _water("1layer")
    water.basis[4] = 0.0
    assert np.all(np.isfinite(water.zDensity_G(np.arange(-10.0, 10.0), 0.0, 0.0)))
    water.type = "layered"
    water.basis[3] = 0.0
    with pytest.raises(ValueError, match="sigma_bar"):
        water.transition_extent()


def test_placements_are_coherent_translations():
    """Direct density and amplitude sums independently check terrace averaging."""
    water = _water("layered")
    flat = copy.deepcopy(water)
    heights = np.array([-3.0, 0.0, 3.0])
    probabilities = np.array([0.2, 0.5, 0.3])
    water._placements = heights, probabilities
    z = np.linspace(-10.0, 20.0, 120)
    expected = sum(
        p * flat.zDensity_G(z - h, 0.0, 0.0) for h, p in zip(heights, probabilities)
    )
    np.testing.assert_allclose(water.zDensity_G(z, 0.0, 0.0), expected)
    q = 0.4 * 2 * np.pi / 4
    expected_f = sum(
        p * np.exp(1j * q * h) for h, p in zip(heights, probabilities)
    ) * flat.F_uc(0.0, 0.0, 0.4)
    np.testing.assert_allclose(water.F_uc(0.0, 0.0, 0.4), expected_f)


@pytest.mark.parametrize("mean", [-2.0, 0.0, 2.0])
def test_attached_water_and_roundtrip(mean):
    """Surface water uses physical terrace tops and survives complete model I/O."""
    from ..CTRdistributions import PoissonProfile
    from ..CTRcalc import SXRDCrystal
    from ._poisson_oracle import poisson_crystal

    crystal, surface = poisson_crystal(PoissonProfile(mean_change=mean, alpha=0.4))
    crystal.setEnergy(10000.0)
    crystal.apply_stacking()
    dry = surface.F_uc(np.array([0.0]), np.array([0.0]), np.array([0.3]))
    surface.attach_water(_water())
    water = surface._sync_water()
    heights, probabilities, tail = surface.terrace_distribution()
    assert abs(probabilities.sum() - 1) < 1e-14
    assert tail < 1e-6
    expected_height = surface.below_H + 3 * mean
    assert abs(heights @ probabilities - expected_height) < 1e-5
    wet = surface.F_uc(np.array([0.0]), np.array([0.0]), np.array([0.3]))
    np.testing.assert_allclose(
        wet - dry, water.F_uc(0.0, 0.0, 0.3) * surface.uc_area / water.uc_area
    )
    restored = SXRDCrystal.fromStr(crystal.toStr())
    np.testing.assert_allclose(restored.F(0.0, 0.0, 0.3), crystal.F(0.0, 0.0, 0.3))


@pytest.mark.parametrize("kind", ["step", "layered", "tanh", "erf"])
def test_water_dwba_ambient_and_reference(kind):
    """Water supplies a non-vacuum ambient and a localized specular contrast."""
    from ..CTRcalc import SXRDCrystal
    from .test_CTRcrystallinity import _cell

    bulk = _cell()
    water = _water(kind)
    crystal = SXRDCrystal(bulk, water)
    prepared = crystal.dwba.prepare_from_glancing(
        0.0, 0.0, np.array([0.02, 0.1]), np.array([0.02, 0.1])
    )
    assert prepared.reference.n[0].real < 1
    result = crystal.dwba.evaluate_prepared(prepared)
    assert np.all(np.isfinite(result.total_amplitude))
    assert result.contributions[-1].role == "water"
    non = crystal.dwba.evaluate_from_glancing(0.2, 0.0, 0.1, 0.1)
    assert non.contributions[-1].F_h == 0


@pytest.mark.parametrize("kind", ["step", "tanh", "erf", "layered"])
def test_attached_dwba_fit_and_cache(kind):
    """Owned water participates in fitting and frozen/automatic preparations."""
    from ..CTRdistributions import PoissonProfile
    from ._poisson_oracle import poisson_crystal

    crystal, surface = poisson_crystal(PoissonProfile(mean_change=-1.0, alpha=0.2))
    crystal.setEnergy(10000.0)
    crystal.apply_stacking()
    surface.attach_water(_water(kind))
    water = surface.water
    water.addFitParameter("sigma_0", limits=(0.05, 0.3), name="water_width")
    values, _, _ = surface.getStartParamAndLimits()
    assert surface.fitparnames[-1] == "water_width"
    values[-1] = 0.19
    surface.setFitParameters(values)
    assert water.basis[2] == 0.19
    original = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.08, 0.08)
    old = crystal.dwba.evaluate_prepared(original)
    water.basis[2] = 0.25
    frozen = crystal.dwba.evaluate_prepared(original)
    assert abs(frozen.F_h - old.F_h) > 1e-8
    fresh = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.08, 0.08)
    if kind in {"tanh", "erf"}:
        assert fresh.reference.fingerprint != original.reference.fingerprint
    if kind == "tanh":
        from ..CTRoptics import LayeredElectricField
        assert isinstance(fresh.field_i, LayeredElectricField)
    water.pw *= 1.01
    with pytest.raises(ValueError, match="ambient changed"):
        crystal.dwba.evaluate_prepared(original)


@pytest.mark.parametrize("kind", ["tanh", "erf", "layered"])
def test_water_reference_records_and_high_angle_limit(kind):
    """Continuum DWBA retains decomposition and the leading kinematic amplitude."""
    from ..CTRcalc import SXRDCrystal
    from .test_CTRcrystallinity import _cell

    crystal = SXRDCrystal(_cell(), _water(kind), atten=0.0)
    p = crystal.dwba.prepare_from_glancing(
        0.0, 0.0, np.array([0.1, 0.15]), np.array([0.1, 0.15])
    )
    if kind == "tanh":
        from ..CTRoptics import LayeredElectricField
        assert isinstance(p.field_i, LayeredElectricField)
    r = crystal.dwba.evaluate_prepared(p)
    np.testing.assert_allclose(r.F_h, sum(c.F_h for c in r.contributions))
    # At high angle the effective amplitude includes reference reflection.
    fkin = crystal.F(*p.hkl)
    np.testing.assert_allclose(abs(r.F_effective), abs(fkin), rtol=0.01)


@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("terraces", [False, True])
@pytest.mark.parametrize("kind", ["erf", "step_erf"])
def test_continuous_dwba_against_fine_slabs(pol, terraces, kind, monkeypatch):
    """A separately solved fine slab reference checks continuum matrix elements."""
    from dataclasses import replace
    from ..CTRcalc import SXRDCrystal
    from ..CTRoptics import StratifiedProfile, solve_electric_field
    from .._CTRwater_dwba import _scale
    from .test_CTRcrystallinity import _cell

    crystal = SXRDCrystal(_cell(), _water(kind), atten=0.0)
    if terraces:
        from ..CTRdistributions import PoissonProfile
        from ._poisson_oracle import poisson_crystal

        crystal, surface = poisson_crystal(PoissonProfile(mean_change=-1.3, alpha=0.2))
        crystal.setEnergy(10000.0)
        surface.attach_water(_water(kind))
    angles = np.array([0.012, 0.07])
    prepared = crystal.dwba.prepare_from_glancing(
        0.0, 0.0, angles, angles, polarization_i=pol, polarization_f=pol
    )
    spec = prepared.reference.continuous
    low, high = spec.bounds()
    # Adsorbed peaks require a finer independent slab discretization than
    # monotonic forward-factor continua; keep the same error threshold.
    nodes = 3001 if kind.startswith("step_") else 1501
    edges = np.unique(np.r_[np.linspace(low, high, nodes), spec.boundaries])
    z = np.r_[edges[0] - 1.0, (edges[:-1] + edges[1:]) / 2, edges[-1] + 1.0]
    density = spec.reference_density(z)
    n = np.sqrt(spec.epsilon(z))
    n[0], n[-1] = prepared.reference.n[-1], prepared.reference.n[0]
    values = np.column_stack((z, 1 - n.real, -n.imag))
    shares = _scale(10000.0) * np.stack((density.real, density.imag), axis=-1)
    shares[:, 0] = prepared._atomic_model.reference_shares[:, -1]
    shares[:, -1] = prepared._atomic_model.reference_shares[:, 0]
    profile = StratifiedProfile(values, edges)
    fields = solve_electric_field(values, 10000.0, angles, pol, boundaries=edges)
    ref = replace(
        prepared.reference, profile=profile, n=n[::-1].copy(), continuous=None
    )
    model = replace(
        prepared._atomic_model, reference_shares=np.ascontiguousarray(shares[:, ::-1])
    )
    sampled = replace(
        prepared,
        reference=ref,
        _atomic_model=model,
        field_i=fields,
        field_f=fields,
        field_i_index=np.arange(2),
        field_f_index=np.arange(2),
    )
    actual = crystal.dwba.evaluate_prepared(prepared)
    # The fine-slab oracle explicitly integrates continuous densities. It is
    # not the production constant-medium molecular-form-factor approximation.
    from .. import _CTRwater_analytic, _CTRwater_dwba
    with monkeypatch.context() as patch:
        patch.setattr(_CTRwater_analytic, "_water_contrast",
                      _CTRwater_dwba._water_contrast)
        oracle = crystal.dwba.evaluate_prepared(sampled)
    # Match phase origins: native fields use the highest interface; continuous
    # fields use the tanh midpoint. Translate the two-wave product accordingly.
    kappa = prepared.k0 * prepared.reference.n[0] * np.sin(angles)
    phase = np.exp(2j * kappa * (fields.z_reference[0] - prepared.field_i.z_reference))
    np.testing.assert_allclose(actual.F_h, oracle.F_h * phase, rtol=8e-4, atol=2e-4)


@pytest.mark.parametrize("mean", [-1.4, 1.4])
def test_correlated_incoherent_water_against_flat_films(mean):
    """Each independent flat Film carries water at its own structural top."""
    from ..CTRcalc import SXRDCrystal
    from ..CTRdistributions import PoissonProfile
    from ..CTRincoherent import PoissonHeightDomains
    from ._poisson_oracle import poisson_crystal, flat_height_crystal

    crystal, surface = poisson_crystal(PoissonProfile(mean_change=mean, alpha=0.0))
    crystal.setEnergy(10000.0)
    surface.attach_water(_water())
    ell = np.linspace(0.1, 1.7, 37)
    lower = int(np.floor(mean))
    frac = mean - lower
    amplitudes = []
    for change in (lower, lower + 1):
        flat = flat_height_crystal(change - 1)
        water = _water()
        flat = SXRDCrystal(
            flat.uc_bulk, *flat.uc_surface_list, water, stacking=np.array([1, 2])
        )
        flat.setEnergy(10000.0)
        amplitudes.append(flat.F(0.0, 0.0, ell))
    expected = (1 - frac) * abs(amplitudes[0]) ** 2 + frac * abs(amplitudes[1]) ** 2
    model = PoissonHeightDomains(crystal, surface="surface", incoherent_fraction=1.0)
    zeros = np.zeros_like(ell)
    np.testing.assert_allclose(model.F2(zeros, zeros, ell), expected, rtol=1e-11)
    np.testing.assert_allclose(
        crystal.F(0.0, 0.0, ell),
        (1 - frac) * amplitudes[0] + frac * amplitudes[1],
        rtol=1e-11,
    )


@pytest.mark.parametrize("kind", ["tanh", "erf", "layered"])
def test_owned_water_hdf5_and_copy(tmp_path, kind):
    """Saved fitting state and copies preserve water identity and errors."""
    from silx.io.dictdump import dicttoh5, h5todict
    from ..CTRcalc import SXRDCrystal
    from ..CTRdistributions import PoissonProfile
    from ._poisson_oracle import poisson_crystal

    crystal, surface = poisson_crystal(PoissonProfile(mean_change=0.0, alpha=0.0))
    crystal.setEnergy(10000.0)
    surface.attach_water(_water(kind))
    surface.water.pw = 0.031
    surface.water.addFitParameter("sigma_0", limits=(0.01, 0.5), name="water_width")
    surface.getStartParamAndLimits()
    surface.setFitParameters([0.2])
    surface.setFitErrors([0.03])
    path = str(tmp_path / "water.h5")
    dicttoh5(surface.parametersToDict(), path)
    restored = SXRDCrystal.fromStr(crystal.toStr())
    target = restored.uc_surface_list[-1]
    target.parametersFromDict(h5todict(path))
    assert target.water.type == ("layered" if kind == "layered" else "step")
    assert target.water._settings() == surface.water._settings()
    assert target.water.pw == 0.031
    np.testing.assert_allclose(target.getInitialParameters(), [0.2])
    np.testing.assert_allclose(target.getFitErrors(), [0.03])
    clone = copy.deepcopy(surface)
    assert clone.water is not surface.water
    assert clone.water._surface_owner is clone
    with pytest.raises(ValueError, match="already"):
        clone.attach_water(surface.water)
    clone.attach_water(None)
    assert clone.water is None


def test_narrow_water_layer_quadrature():
    """Narrow excess peaks remain resolved when smaller than a wave period."""
    from ..CTRcalc import SXRDCrystal
    from .._CTRwater_dwba import _water_contrast, _water_contrast_order
    from .test_CTRcrystallinity import _cell

    water = _water("1layer_O")
    water.basis[2] = 1e-4
    crystal = SXRDCrystal(_cell(), water)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.09, 0.09)
    records = crystal.dwba._live_records()
    actual = _water_contrast(prepared, records[-1], len(records) - 1)
    oracle = _water_contrast_order(prepared, records[-1], len(records) - 1, 96)
    np.testing.assert_allclose(actual, oracle, rtol=1e-7, atol=1e-9)


@pytest.mark.parametrize("broadening", [0.1, 0.25])
def test_layered_contrast_against_sampled_fields(broadening):
    """Independent Simpson integration checks short and long layering tails."""
    from ..CTRcalc import SXRDCrystal
    from ..CTRoptics import sample_electric_field
    from .._CTRwater_dwba import _scale, _water_contrast
    from .test_CTRcrystallinity import _cell

    water = _water("layered")
    water.basis[3] = broadening
    crystal = SXRDCrystal(_cell(), water)
    angles = np.array([0.02, 0.08])
    prepared = crystal.dwba.prepare_from_glancing(0., 0., angles, angles)
    records = crystal.dwba._live_records()
    actual = _water_contrast(prepared, records[-1], len(records) - 1)
    low, high = water.transition_extent()
    boundaries = prepared.reference.profile.boundaries
    low, high = min(low, boundaries[0]), max(high, boundaries[-1])
    edges = np.unique(np.r_[low, boundaries, high])
    shares = prepared._atomic_model.reference_shares[-1]
    oracle = np.zeros(len(angles), complex)
    for left, right in zip(edges[:-1], edges[1:]):
        # Reference density is constant on each open interval; keep its
        # one-sided value at the endpoints for the Simpson rule.
        medium = len(shares) - 1 - np.searchsorted(boundaries, (left + right)/2)
        z = np.linspace(left, right, max(3, 2*int(np.ceil((right-left)/0.01))+1))
        contrast = (
            water.zDensity_G(z, 0., 0.) - complex(*shares[medium])/_scale(water._E)
        )
        field = sample_electric_field(prepared.field_i, z)
        oracle += (
            simpson(contrast[:, None]*field**2, x=z, axis=0)*prepared.reference_area
        )
    np.testing.assert_allclose(actual, oracle, rtol=2e-6, atol=2e-7)


@pytest.mark.parametrize("kind", ["step", "tanh", "erf", "layered"])
def test_water_optimizer_local_scan_and_external_angle_rejection(kind):
    """Generated local geometry fits water; external measured angles need optics."""
    from ..CTRcalc import SXRDCrystal
    from ..CTRopt import CTROptimizer
    from ..CTRplotutil import CTRCollection
    from .test_CTRopt_dwba import _ctr
    from .test_CTRcrystallinity import _cell

    water = _water(kind)
    water.addFitParameter("sigma_0", limits=(0.05, 0.5), name="water_width")
    crystal = SXRDCrystal(_cell(), water)
    rod = _ctr(crystal, records=False)
    optimizer = CTROptimizer(crystal, CTRCollection([rod]), scale_policy={"R": "fixed"})
    optimizer.set_dwba()
    optimizer.prepareFit()
    a = optimizer.flat_prediction([0.15]).copy()
    b = optimizer.flat_prediction([0.3]).copy()
    assert not np.array_equal(a, b)
    np.testing.assert_array_equal(a, optimizer.flat_prediction([0.15]))
    with pytest.raises(ValueError, match="External Vlieg"):
        crystal.dwba.prepare_from_vlieg(0.02, 0.0, 0.02, 0.0, 0.0, 0.0)


def test_water_unpolarized_and_tanh_tail_tolerance():
    """Unpolarized water is the s/p intensity mean; tail tolerance is explicit."""
    from ..CTRcalc import SXRDCrystal
    from .test_CTRcrystallinity import _cell

    water = _water("tanh")
    low, high = water.transition_extent(tolerance=1e-12)
    np.testing.assert_allclose(
        water.zDensity_G(np.array([low, high]), 0.0, 0.0) / water.bulk_density,
        [0.0, 1.0],
        atol=1.1e-12,
    )
    water.optical_reference = "sharp"  # Analytic s and p; tanh is s-only.
    crystal = SXRDCrystal(_cell(), water)
    args = (0.0, 0.0, 0.03, 0.03)
    expected = (
        sum(
            crystal.dwba.evaluate_from_glancing(
                *args, polarization_i=pol, polarization_f=pol
            ).reflectivity
            for pol in ("s", "p")
        )
        / 2
    )
    actual = crystal.dwba.polarization_reflectivity_from_glancing(
        *args, s_fraction=0.5, outgoing="unanalysed"
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-12)


def test_continuous_integral_cache_keeps_live_atoms_and_invalidates_water(monkeypatch):
    """Frozen integrals are reused, while atomic and water changes remain live."""
    from ..CTRcalc import SXRDCrystal
    from .test_CTRcrystallinity import _cell
    from .._CTRwater_continuous import _WaterReference

    water = _water("erf")  # Explicit numerical confirmation path only.
    crystal = SXRDCrystal(_cell(), water)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.05, 0.05)
    initial = crystal.dwba.evaluate_prepared(prepared)
    original = _WaterReference.reference_density
    calls = []

    def counted(self, z):
        """Record density evaluations without changing the reference."""
        calls.append(np.size(z))
        return original(self, z)

    monkeypatch.setattr(_WaterReference, "reference_density", counted)
    crystal.uc_bulk.basis[0, 3] += 0.01
    changed = crystal.dwba.evaluate_prepared(prepared)
    assert changed.F_h != initial.F_h
    assert calls == []
    water.basis[2] *= 1.1
    crystal.dwba.evaluate_prepared(prepared)
    assert calls
    assert crystal.dwba.cache_info()["continuous_integral_size"] == 2
    crystal.dwba.clear_cache()
    assert crystal.dwba.cache_info()["continuous_integral_size"] == 0


@pytest.mark.parametrize("kind", ["layered", "step_sharp"])
@pytest.mark.parametrize("polarization", ["s", "p"])
def test_analytic_water_keeps_live_atoms_and_water(kind, polarization, monkeypatch):
    """Live fits need no integral cache; preserve frozen-ambient checks."""
    from ..CTRcalc import SXRDCrystal
    from .. import _CTRwater_dwba
    from .test_CTRcrystallinity import _cell

    water = _water(kind)
    crystal = SXRDCrystal(_cell(), water)
    angles = np.array([0.02, 0.08, 0.13])
    prepared = crystal.dwba.prepare_from_glancing(
        0., 0., angles, angles, polarization_i=polarization, polarization_f=polarization
    )
    initial = crystal.dwba.evaluate_prepared(prepared)
    calls = []
    original = _CTRwater_dwba._water_contrast_order

    def counted(*args, **kwargs):
        """Count expensive integrations while retaining their actual values."""
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(_CTRwater_dwba, "_water_contrast_order", counted)
    np.testing.assert_array_equal(
        crystal.dwba.evaluate_prepared(prepared).total_amplitude,
        initial.total_amplitude,
    )
    crystal.uc_bulk.basis[0, 3] += 0.01
    changed_atoms = crystal.dwba.evaluate_prepared(prepared)
    assert np.max(abs(changed_atoms.F_h-initial.F_h)) > 1e-8
    assert not calls
    water.basis[2] *= 1.1
    changed_water = crystal.dwba.evaluate_prepared(prepared)
    assert not calls
    assert np.max(abs(changed_water.F_h-changed_atoms.F_h)) > 1e-8
    records = crystal.dwba._live_records()
    from .._CTRwater_analytic import _water_contrast
    expected = _water_contrast(prepared, records[-1], len(records)-1)
    np.testing.assert_allclose(changed_water.contributions[-1].F_h, expected,
                               rtol=1e-12, atol=1e-12)
    for _ in range(5):
        water.basis[2] *= 1.01
        crystal.dwba.evaluate_prepared(prepared)
    assert crystal.dwba.cache_info()["water_contrast_size"] == 0
    water.pw *= 1.01
    with pytest.raises(ValueError, match="ambient changed"):
        crystal.dwba.evaluate_prepared(prepared)
    crystal.dwba.clear_cache()
    assert crystal.dwba.cache_info()["water_contrast_size"] == 0


@pytest.mark.parametrize("kind", WaterModel.supportedTypes + WaterModel.legacyTypes)
@pytest.mark.parametrize("sigma", [0., .15])
def test_step_density_scalar_vector_agree(kind, sigma):
    """ODE scalar sampling retains the vector density, including legacy factors."""
    from .._CTRwater import _step_density_function

    water = _water(kind)
    if water.type != "step":
        return
    water.basis[2] = sigma
    if sigma == 0:
        water.basis[4] = 0.  # A delta adsorption layer cannot be point-sampled.
    water._placements = (np.array([-2., 0., 3.]), np.array([.2, .5, .3]))
    density = _step_density_function(water, averaged=True)
    z = np.linspace(-4., 8., 71)
    np.testing.assert_allclose(np.array([density(v) for v in z]), density(z),
                               rtol=2e-14, atol=1e-15)


@pytest.mark.parametrize("width", [0.12, 0.65])
@pytest.mark.parametrize("kind", ["tanh", "step_tanh"])
def test_piecewise_poisson_reference_without_ode(width, kind, monkeypatch):
    """Both separated and overlapping terraces retain their physical density."""
    from .. import CTRcontinuous
    from ..CTRdistributions import PoissonProfile
    from ._poisson_oracle import poisson_crystal

    crystal, surface = poisson_crystal(PoissonProfile(mean_change=-1.3, alpha=0.2))
    crystal.setEnergy(10000.0)
    water = _water(kind)
    water.basis[2] = width
    surface.attach_water(water)

    def forbidden(*args, **kwargs):
        """Prevent a regression to numerical propagation on the fitting path."""
        raise AssertionError("Tanh s-polarized Poisson fields must remain analytic")

    monkeypatch.setattr(CTRcontinuous, "solve_ivp", forbidden)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, [0.01, 0.06], [0.01, 0.06])
    result = crystal.dwba.evaluate_prepared(prepared)
    assert np.all(np.isfinite(result.total_amplitude))
    spec = prepared.reference.continuous
    z = np.linspace(*spec.bounds(), 501)
    reference = spec.reference_density(z)
    physical = spec.atomic_density(z).sum(axis=0) + water.zDensity_G(z, 0.0, 0.0)
    residual = water.zDensity_G(z, 0.0, 0.0) - reference[spec.water_index]
    np.testing.assert_allclose(
        reference.sum(axis=0) + residual, physical, rtol=2e-14, atol=2e-14
    )


@pytest.mark.parametrize("pol", ["s", "p"])
@pytest.mark.parametrize("kind", ["erf", "step_erf"])
def test_erf_reference_is_full_profile_and_converges(pol, kind):
    """The confirmation reference solves the full averaged erf susceptibility."""
    from ..CTRdistributions import PoissonProfile
    from ._poisson_oracle import poisson_crystal

    crystal, surface = poisson_crystal(PoissonProfile(mean_change=-1.3, alpha=0.2))
    crystal.setEnergy(10000.0)
    surface.attach_water(_water(kind))
    prepared = crystal.dwba.prepare_from_glancing(
        0.0,
        0.0,
        [0.002, 0.01, 0.06],
        [0.002, 0.01, 0.06],
        polarization_i=pol,
        polarization_f=pol,
    )
    spec = prepared.reference.continuous
    z = np.linspace(*spec.bounds(), 151)
    physical = spec.atomic_density(z).sum(axis=0) + surface.water.zDensity_G(
        z, 0.0, 0.0
    )
    np.testing.assert_allclose(
        spec.epsilon(z), (1 - spec._scale * physical) ** 2, rtol=0.0, atol=5e-16
    )
    tighter = spec.solve(prepared.field_i.alpha, pol, rtol=2e-12)
    np.testing.assert_allclose(prepared.field_i.r_S, tighter.r_S, rtol=2e-6, atol=2e-9)
    np.testing.assert_allclose(
        prepared.field_i.sample(z)[0], tighter.sample(z)[0], rtol=2e-6, atol=2e-9
    )

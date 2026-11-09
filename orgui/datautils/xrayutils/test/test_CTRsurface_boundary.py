"""Native API contracts for finite material boundaries on both surface laws."""

import copy
from dataclasses import FrozenInstanceError

from ..CTRfilm import LegacySurfaceBoundary

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from ..CTRcalc import SXRDCrystal
from ..CTRdistributions import (
    HeightSupport,
    PoissonProfile,
    RecenteredSkellamProfile,
    bounded_height_states,
)
from ..CTRfilm import Film, FiniteFilmBoundary, PoissonSurface, SkellamSurface
from ..CTRincoherent import SurfaceHeightDomains, create_incoherent_model
from .test_CTRdissolution import _crystal, _flat, _hkl


def _native(kind, profile, width=6, **kwargs):
    template = _crystal(width=width)
    surface = kind(
        copy.deepcopy(template["film"].unitcell),
        name="surface",
        profile=profile,
        **kwargs,
    )
    crystal = SXRDCrystal(
        template.uc_bulk, template["film"], surface, stacking=np.array([1, 2])
    )
    crystal.setEnergy(10000)
    crystal.dwba.set_ctr_geometry(alpha_i=0.03)
    return crystal


@pytest.mark.parametrize(
    "mean,alpha,offset", [(-8, 0, 0.3), (-4, 0.7, -2.3), (4, 0.8, -5.1), (0, 1, -5.25)]
)
def test_signed_poisson_bounded_moments(mean, alpha, offset):
    """Independent native Poisson PMF agrees with the bounded representation."""
    profile = PoissonProfile(mean, alpha, offset)
    indices, probabilities, certificate = bounded_height_states(profile, -6)
    raw = np.arange(-300, 301)
    pmf = profile.probability(raw)
    bounded = np.maximum(raw, -6)
    assert_allclose(probabilities[0], pmf[raw <= -6].sum(), atol=2e-15)
    assert_allclose(probabilities.sum(), 1.0, atol=2e-15)
    assert_allclose(indices @ probabilities, bounded @ pmf, atol=1e-10)
    assert_allclose(indices**2 @ probabilities, bounded**2 @ pmf, atol=1e-8)
    assert certificate.omitted_second_moment <= 1e-8 / 4


@pytest.mark.parametrize(
    "kind,profile",
    [
        (SkellamSurface, RecenteredSkellamProfile(offset=-5.25)),
        (PoissonSurface, PoissonProfile(-5.0, alpha=0, offset=-0.25)),
    ],
)
def test_native_coherent_complete_states_and_density(kind, profile, monkeypatch):
    """Coherent native entry points include the same entire bounded stack."""
    crystal = _native(kind, profile)
    basis = crystal["surface"].basis.copy()

    def forbidden(*args, **kwargs):
        pytest.fail("Bounded binding requested the legacy distribution support")

    monkeypatch.setattr(kind, "_profile_candidates", forbidden)
    f0, f1 = (_flat(n).F(*_hkl()) for n in (0, 1))
    expected = 0.25 * f0 + 0.75 * f1
    result = crystal.evaluate_kinematic(*_hkl())
    assert [c.name for c in result.components] == ["film", "surface"]
    assert_allclose(result.total, expected, atol=1e-10)
    assert_allclose(result.bulk + sum(c.amplitude for c in result.components), expected)
    assert_allclose(crystal.F_surf(*_hkl()), expected - result.bulk)
    assert_allclose(crystal.F2(*_hkl()), abs(expected) ** 2)
    domains = SurfaceHeightDomains(crystal, surface="surface")
    assert_allclose(domains.F2(*_hkl()), crystal.F2(*_hkl()), atol=1e-10)
    z = np.linspace(-10, 22, 65)
    assert_allclose(
        crystal.zDensity_G(z, 0, 0),
        0.25 * _flat(0).zDensity_G(z, 0, 0) + 0.75 * _flat(1).zDensity_G(z, 0, 0),
    )
    stats = crystal.surface_height_statistics(surface="surface")
    assert stats["parent_mean_index"] == -5.25
    assert stats["parent_variance_index2"] == 0.1875
    assert stats["realized_mean_index"] == -5.25
    assert stats["realized_variance_index2"] == .1875
    assert stats["exhausted_probability"] == 0.25
    assert stats["mean_front_height_A"] == 2.25
    assert (
        crystal.surface_height_state(surface="surface", index=-6).inventory["film"]
        == {}
    )
    assert_array_equal(crystal["surface"].basis, basis)


@pytest.mark.parametrize(
    "kind,profile",
    [
        (SkellamSurface, RecenteredSkellamProfile(1, 0, 3, -5.25)),
        (PoissonSurface, PoissonProfile(-3, 1, -2.25)),
    ],
)
@pytest.mark.parametrize("channels", [("s", "s"), ("s", "p"), ("p", "s"), ("p", "p")])
def test_native_dwba_equals_shared_coherent_endpoint(kind, profile, channels):
    """Native DWBA retains the all-state optical grid, including exposed bulk."""
    crystal = _native(kind, profile)
    options = dict(polarization_i=channels[0], polarization_f=channels[1])
    native = crystal.dwba.evaluate(*_hkl(), **options)
    domains = SurfaceHeightDomains(
        crystal, surface="surface", optical_reference="shared"
    )
    ensemble = domains.dwba.evaluate(*_hkl(), **options)
    assert_allclose(native.reflectivity, ensemble.reflectivity, rtol=1e-10, atol=1e-28)
    assert_array_equal(
        native.prepared.reference.profile.values,
        ensemble.prepared.reference.profile.values,
    )
    assert_array_equal(
        crystal.optical_profile(), native.prepared.reference.profile.values
    )
    with crystal.dwba.batch():
        first = crystal.dwba.evaluate(*_hkl(), **options)
        second = crystal.dwba.evaluate(*_hkl(), **options)
        assert_array_equal(first.reflectivity, second.reflectivity)
    # The source still owns preparation identity and geometry invalidation.
    crystal["film"].basis[0] += 2
    with pytest.raises(ValueError):
        crystal.dwba.evaluate_prepared(native.prepared)


@pytest.mark.parametrize(
    "kind,profile",
    [
        (SkellamSurface, RecenteredSkellamProfile(offset=-5.25)),
        (PoissonSurface, PoissonProfile(-5, 0, -0.25)),
    ],
)
def test_native_persistence_guard_and_no_companion_requirement(kind, profile):
    """A standalone crystal retains the physical boundary without wrapper config."""
    crystal = _native(kind, profile)
    text = crystal.toStr(showErrors=False)
    assert "RequiresFiniteFilmBoundary" in text
    restored = SXRDCrystal.fromStr(text)
    assert restored["surface"].boundary == crystal["surface"].boundary
    assert restored["surface"].height_support == HeightSupport()
    assert_allclose(restored.F(*_hkl()), crystal.F(*_hkl()), atol=1e-10)
    settings = crystal.parametersToDict()
    restored.parametersFromDict(settings)
    assert_allclose(restored.F(*_hkl()), crystal.F(*_hkl()), atol=1e-10)
    # Old RoughSurface readers skip the first line and then parse numbers:
    # the second mandatory non-comment token must make that fail.
    surface_text = crystal["surface"].toStr()
    with pytest.raises(ValueError):
        np.array(surface_text.splitlines()[1].split(), dtype=float)
    with pytest.raises(ValueError, match="guard"):
        kind.fromStr(surface_text.replace("RequiresFiniteFilmBoundary", "Ignored"))


def test_default_native_boundary_and_unpublished_api_removal():
    """New surfaces dissolve by default; a named policy keeps legacy correction."""
    from .. import CTRincoherent

    assert not hasattr(CTRincoherent, "FiniteFilmDissolution")
    with pytest.raises((ValueError, KeyError)):
        create_incoherent_model("finite_film_dissolution", _crystal())
    crystal = _crystal()
    crystal["surface"].basis[:] = [1, 0, 10, -30]
    crystal["surface"].boundary = FiniteFilmBoundary()
    crystal["surface"].height_support = HeightSupport()
    assert np.all(np.isfinite(crystal.F(*_hkl())))
    crystal["surface"].boundary = LegacySurfaceBoundary()
    crystal["surface"].height_support = None
    with pytest.raises(ValueError, match="etching support"):
        crystal.F(*_hkl())
    with pytest.raises(FrozenInstanceError):
        FiniteFilmBoundary(("film",)).films = ("other",)
    with pytest.raises(ValueError, match="height_support"):
        PoissonSurface(crystal["film"].unitcell, boundary=LegacySurfaceBoundary(),
                       height_support=HeightSupport())


def test_native_boundary_retains_overlayer_and_protects_bulk():
    crystal = _native(SkellamSurface, RecenteredSkellamProfile(offset=-50))
    overlayer = Film(
        copy.deepcopy(crystal["film"].unitcell), width_layers=2, name="cap"
    )
    capped = SXRDCrystal(
        crystal.uc_bulk,
        *crystal.uc_surface_list,
        overlayer,
        stacking=np.array([1, 2, 3]),
    )
    capped.setEnergy(10000)
    assert_allclose(capped.F(*_hkl()), _flat(2).F(*_hkl()), atol=1e-10)
    capped.apply_stacking()
    assert capped["cap"].pos_absolute == 0
    assert_allclose(crystal.F(*_hkl()), _flat(0).F(*_hkl()), atol=1e-10)
    stats = crystal.surface_height_statistics(surface="surface")
    assert stats["parent_mean_index"] == -50
    assert stats["realized_mean_index"] == -6
    assert stats["realized_variance_index2"] == 0


@pytest.mark.parametrize(
    "kind,profile",
    [
        (SkellamSurface, RecenteredSkellamProfile(offset=-5.25)),
        (PoissonSurface, PoissonProfile(-5.0, 0.0, -0.25)),
    ],
)
def test_fitted_native_height_hdf5_and_wrapper_config(tmp_path, kind, profile):
    """Native fitted height parameters survive HDF5, independently of mixing."""
    from silx.io.dictdump import dicttoh5, h5todict
    from ..CTRincoherent import create_incoherent_model_from_config

    crystal = _native(kind, profile)
    target = crystal["surface"]
    target.addFitParameter("offset", name="offset", limits=(-10, 0))
    model = SurfaceHeightDomains(crystal, surface="surface", incoherent_fraction=0.3)
    model.setParameters([-7.0])
    path = str(tmp_path / "surface.h5")
    dicttoh5(model.parametersToDict(), path)
    other = SurfaceHeightDomains(
        SXRDCrystal.fromStr(crystal.toStr()), surface="surface"
    )
    other.parametersFromDict(h5todict(path))
    restored = create_incoherent_model_from_config(
        other.coherent_model, model.to_config()
    )
    assert_allclose(restored.F2(*_hkl()), model.F2(*_hkl()))
    assert restored.target_surface().boundary == target.boundary
    assert restored.target_surface().fitparnames == ["offset"]
    assert restored.basis.shape == (1,)


def test_stratified_profile_preserves_union_boundaries():
    """Wavefields inherit physical edges instead of reconstructing center midpoints."""
    from ..CTRoptics import solve_electric_field

    crystal = _native(SkellamSurface, RecenteredSkellamProfile(offset=-5.25))
    reference = crystal.dwba._atomic_model_snapshot()[2]
    public = crystal.stratified_optical_profile(delta_tolerance=None)
    assert_array_equal(public.boundaries, reference.boundaries)
    assert_array_equal(public.values, reference.values)
    # Unsplit and merged fields must describe the same physical medium.
    expected = solve_electric_field(
        public.values, 10000, 0.03, boundaries=public.boundaries
    )
    actual = crystal.wavefield(0.03, delta_tolerance=None)
    assert_allclose(actual.kz, expected.kz)


def test_local_surface_scattering_cannot_hide_material_removal():
    crystal = _native(SkellamSurface, RecenteredSkellamProfile(offset=-5.25))
    crystal.apply_stacking()
    with pytest.raises(ValueError, match="complete SXRDCrystal"):
        crystal["surface"].F_uc(*_hkl())


def test_saved_boundary_rejects_unknown_owners_and_versions():
    crystal = _native(SkellamSurface, RecenteredSkellamProfile(offset=-5.25))
    crystal["surface"].boundary = FiniteFilmBoundary(("film",))
    text = crystal.toStr()
    with pytest.raises(ValueError, match="bulk is protected"):
        SXRDCrystal.fromStr(text.replace('"films": ["film"]', '"films": ["bulk"]'))
    with pytest.raises(ValueError, match="Unsupported"):
        SXRDCrystal.fromStr(text.replace('"schema_version": 1', '"schema_version": 2'))
    assert_allclose(SXRDCrystal.fromStr(text.replace("\n", "\r\n")).F(*_hkl()),
                    crystal.F(*_hkl()))


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
def test_default_removes_all_finite_films_and_explicit_selection_limits_it(kind):
    """Finite substrate is removable by default; the infinite bulk never is."""
    # Native adjacent-Film widths are cumulative: 40 lower + 6 upper planes.
    template = _crystal(width=46, substrate_film=True)
    surface = kind(copy.deepcopy(template["film"].unitcell), name="surface")
    surface.basis[surface.parameterLookup["offset"]] = -100
    crystal = SXRDCrystal(
        template.uc_bulk, *template.uc_surface_list[:-1], surface,
        stacking=template.uc_stacking.copy())
    crystal.setEnergy(10000)
    crystal.dwba.set_ctr_geometry(alpha_i=.03)
    assert surface.boundary == FiniteFilmBoundary()
    assert surface.boundary.films is None
    expected = _flat(0)
    assert_allclose(crystal.F(*_hkl()), expected.F(*_hkl()), atol=1e-10)
    stats = crystal.surface_height_statistics(surface="surface")
    assert stats["realized_mean_index"] == -46
    for fraction in (0, .4, 1):
        domains = SurfaceHeightDomains(
            crystal, surface="surface", incoherent_fraction=fraction,
            optical_reference="shared")
        assert_allclose(domains.F2(*_hkl()), expected.F2(*_hkl()), atol=1e-10)
        assert_allclose(domains.dwba.evaluate(*_hkl()).reflectivity,
                        expected.dwba.evaluate(*_hkl()).reflectivity,
                        rtol=1e-10, atol=1e-28)
    # A named upper sequence protects the explicit substrate Film as well.
    surface.boundary = FiniteFilmBoundary(("film",))
    restricted = _flat(40, occupancy=1)
    assert_allclose(crystal.F(*_hkl()), restricted.F(*_hkl()), atol=1e-10)
    assert crystal.surface_height_statistics(surface="surface")[
        "realized_mean_index"] == -6


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
def test_inferred_selection_persists_without_freezing_stack(tmp_path, kind):
    """Text/HDF5 preserve inference; attaching a finite lower Film extends it."""
    from silx.io.dictdump import dicttoh5, h5todict

    profile = (PoissonProfile(-100, 0) if kind is PoissonSurface
               else RecenteredSkellamProfile(offset=-100))
    original = _native(kind, profile)
    restored = SXRDCrystal.fromStr(original.toStr())
    path = str(tmp_path / "inferred.h5")
    dicttoh5(original.parametersToDict(), path)
    restored.parametersFromDict(h5todict(path))
    assert restored["surface"].boundary == FiniteFilmBoundary()
    lower = Film(copy.deepcopy(restored.uc_bulk), width_layers=4, name="lower")
    expanded = SXRDCrystal(
        restored.uc_bulk, lower, *restored.uc_surface_list,
        stacking=np.array([1, 2, 3]))
    expanded.setEnergy(10000)
    # The same nominal total height now contains four lower and two upper planes.
    assert expanded._surface_state_plan().films == ("lower", "film")
    assert expanded.surface_height_statistics(surface="surface")[
        "realized_mean_index"] == -6
    assert_allclose(expanded.F(*_hkl()), _flat(0).F(*_hkl()), atol=1e-10)


def test_empty_explicit_selection_is_not_unrestricted():
    with pytest.raises(ValueError, match="distinct finite Film names"):
        FiniteFilmBoundary(())


@pytest.mark.parametrize("kind,profile", [
    (PoissonSurface, PoissonProfile(-50, .8)),
    (SkellamSurface, RecenteredSkellamProfile(1, 0, 50, -40)),
])
def test_default_broad_distribution_includes_mixed_interface(kind, profile):
    """Inferred ownership matches explicit all-Film selection at mixed interfaces."""
    template = _crystal(width=60, interface=True, substrate_film=True)
    surface = kind(copy.deepcopy(template["film"].unitcell), name="surface",
                   profile=profile)
    crystal = SXRDCrystal(
        template.uc_bulk, *template.uc_surface_list[:-1], surface,
        stacking=template.uc_stacking.copy())
    crystal.setEnergy(10000)
    crystal.dwba.set_ctr_geometry(alpha_i=.03)
    explicit = copy.deepcopy(crystal)
    explicit["surface"].boundary = FiniteFilmBoundary(("substrate_film", "film"))
    assert_allclose(crystal.F(*_hkl()), explicit.F(*_hkl()), atol=1e-10)
    stats = crystal.surface_height_statistics(surface="surface")
    assert stats == explicit.surface_height_statistics(surface="surface")
    assert 0 < stats["exhausted_probability"] < 1
    exhausted = crystal.surface_height_state(
        surface="surface", index=crystal._surface_state_plan().minimum)
    assert exhausted.inventory["film"] == {}
    assert exhausted.inventory["substrate_film"] == {}


@pytest.mark.parametrize("kind", [PoissonSurface, SkellamSurface])
def test_legacy_files_keep_explicit_legacy_behavior(kind, tmp_path):
    """Absent persisted metadata is legacy even though constructors changed default."""
    cell = _crystal()["film"].unitcell
    legacy = kind(cell, boundary=LegacySurfaceBoundary())
    assert legacy.toStr().startswith("LegacySurfaceBoundaryV1 ")
    assert legacy.parametersToDict()["height_boundary"] == {
        "kind": "legacy", "schema_version": 1}
    restored = kind.fromStr(legacy.toStr())
    assert restored.boundary == LegacySurfaceBoundary()
    assert restored.height_support is None
    fresh = kind(copy.deepcopy(cell))
    assert fresh.boundary == FiniteFilmBoundary()
    fresh.parametersFromDict(legacy.parametersToDict())
    assert fresh.boundary == LegacySurfaceBoundary()
    assert fresh.height_support is None
    old_text = legacy.toStr().split("\n", 1)[1]
    assert kind.fromStr(old_text).boundary == LegacySurfaceBoundary()
    old_parameters = legacy.parametersToDict()
    del old_parameters["height_boundary"]
    fresh.parametersFromDict(old_parameters)
    assert fresh.boundary == LegacySurfaceBoundary()
    from silx.io.dictdump import dicttoh5, h5todict
    path = str(tmp_path / "legacy.h5")
    dicttoh5(legacy.parametersToDict(), path)
    fresh.parametersFromDict(h5todict(path))
    assert fresh.boundary == LegacySurfaceBoundary()


@pytest.mark.parametrize("kind,profile", [
    (PoissonSurface, PoissonProfile(-7, .8)),
    (SkellamSurface, RecenteredSkellamProfile(1, 0, 7, -5)),
])
def test_none_and_omission_are_identical_in_both_scattering_paths(kind, profile):
    default = _native(kind, profile)
    explicit_none = _native(kind, profile, boundary=None)
    assert explicit_none["surface"].boundary == FiniteFilmBoundary()
    assert_array_equal(default.F(*_hkl()), explicit_none.F(*_hkl()))
    assert_array_equal(default.dwba.evaluate(*_hkl()).reflectivity,
                       explicit_none.dwba.evaluate(*_hkl()).reflectivity)
    for fraction in (0, .4, 1):
        models = [SurfaceHeightDomains(
            crystal, surface="surface", incoherent_fraction=fraction,
            optical_reference="shared") for crystal in (default, explicit_none)]
        assert_array_equal(models[0].F2(*_hkl()), models[1].F2(*_hkl()))
        assert_array_equal(models[0].dwba.evaluate(*_hkl()).reflectivity,
                           models[1].dwba.evaluate(*_hkl()).reflectivity)


@pytest.mark.parametrize("kind,profile", [
    (PoissonSurface, PoissonProfile(-1, .5)),
    (SkellamSurface, RecenteredSkellamProfile(1, 0, .5)),
])
def test_old_pickle_policy_and_explicit_none_assignment(kind, profile):
    import pickle

    legacy = _native(kind, profile, width=80, boundary=LegacySurfaceBoundary())
    expected = legacy.F(*_hkl())
    for stored_none in (False, True):
        old = copy.deepcopy(legacy)
        del old["surface"].__dict__["_boundary"]
        if stored_none:
            old["surface"].__dict__["boundary"] = None
        restored = pickle.loads(pickle.dumps(old))
        assert restored["surface"].boundary == LegacySurfaceBoundary()
        assert_array_equal(restored.F(*_hkl()), expected)
        restored["surface"].boundary = None
        assert restored["surface"].boundary == FiniteFilmBoundary()
        assert restored["surface"].height_support == HeightSupport()
        assert_allclose(restored.F(*_hkl()),
                        _native(kind, profile, width=80).F(*_hkl()), atol=1e-10)


def test_default_requires_an_explicit_finite_inventory():
    template = _crystal()
    surface = PoissonSurface(copy.deepcopy(template.uc_bulk), name="surface")
    crystal = SXRDCrystal(template.uc_bulk, surface, stacking=np.array([1]))
    with pytest.raises(ValueError, match="at least one finite Film"):
        crystal.F(*_hkl())


def test_inferred_selection_does_not_hide_unordered_interface_fronts():
    """Automatic selection preserves topology validation instead of falling back."""
    from .test_CTRinterface_anchors import _stack

    crystal = _stack(width=64.25, registry=2, coupling=0., mismatch=.85,
                     offset=.2, interface_width=.5)
    crystal.setParameters([0.])
    crystal["surface"].boundary = FiniteFilmBoundary()
    crystal["surface"].height_support = HeightSupport()
    with pytest.raises(ValueError, match="front heights are not ordered"):
        crystal.F(*_hkl())

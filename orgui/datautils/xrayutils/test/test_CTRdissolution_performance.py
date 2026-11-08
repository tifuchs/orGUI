"""Native complete-state oracles and deterministic dissolution work budgets."""

from collections import Counter
import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from ..CTRdissolution import DissolutionPath
from ..CTRdistributions import RecenteredSkellamProfile
from ..CTRincoherent import FiniteFilmDissolution
from ..CTRuc import UnitCell
from .test_CTRdissolution import _crystal, _hkl
from ._dissolution_dwba_oracle import _complete_dwba_batch, _complete_dwba_oracle


@pytest.mark.parametrize("interface,substrate", [(False, False), (True, False),
                                              (True, True)])
@pytest.mark.parametrize("exposure", ["preserve", "native_or_unreconstructed"])
@pytest.mark.parametrize("shape", ["scalar", "broadcast"])
def test_prefix_amplitudes_match_complete_native_states(interface, substrate,
                                                       exposure, shape):
    """Retain signed interfaces, outer domains, growth and exhaustion phases."""
    source = _crystal(width=40 if interface else 6, interface=interface,
                      substrate_film=substrate)
    source.atten = .023
    source["film"].unitcell.crystallinity = .63
    rotation = np.diag([-1., -1., 1.])
    source.domains = [[(np.eye(3), .7), (rotation, .3)]
                      for _ in source.uc_surface_list]
    films = ["film", "substrate_film"] if substrate else ["film"]
    path = DissolutionPath(source, films=films, surface="surface",
                           exposure_policy=exposure)
    indices = sorted({path.minimum, path.minimum+1, -3, -1, 0, 1, 4})
    hkl = ((1., -1., .43) if shape == "scalar" else
           (np.array([[0.], [1.]]), 0., np.array([.02, .43, .9, 1.99])))
    before = source.F(*hkl)
    expected = [path.state(index).crystal.F(*hkl) for index in indices]
    actual = list(path._iter_kinematic_amplitudes(indices, *hkl))
    assert [index for index, _ in actual] == indices
    for (_, value), reference in zip(actual, expected):
        assert_allclose(value, reference, rtol=2e-12, atol=2e-10)
    assert_array_equal(source.F(*hkl), before)


@pytest.mark.parametrize("fraction", [0., .2054, 1.])
def test_kinematic_work_and_mixture(monkeypatch, fraction):
    """Native fast path evaluates bulk once and constructs no full states."""
    model = FiniteFilmDissolution(
        _crystal(), films=["film"], surface="surface",
        profile=RecenteredSkellamProfile(1., 2., 3., -4.25),
        incoherent_fraction=fraction)
    path, indices, probabilities, _ = model._snapshot()
    amplitudes = np.array([path.state(index).crystal.F(*_hkl())
                           for index in indices])
    expected = ((1-fraction)*abs(probabilities @ amplitudes)**2
                + fraction*(probabilities @ abs(amplitudes)**2))
    counts = Counter()
    original = UnitCell.F_bulk

    def bulk(*args, **kwargs):
        counts["bulk"] += 1
        return original(*args, **kwargs)

    def forbidden(*args, **kwargs):
        pytest.fail("Kinematic evaluation constructed full states or DWBA optics")

    from .. import CTRdwba
    monkeypatch.setattr(UnitCell, "F_bulk", bulk)
    monkeypatch.setattr(DissolutionPath, "state", forbidden)
    monkeypatch.setattr(CTRdwba.DWBAState, "_atomic_model_snapshot", forbidden)
    monkeypatch.setattr(CTRdwba, "solve_electric_field", forbidden)
    assert_allclose(model.F2(*_hkl()), expected, rtol=2e-12, atol=1e-9)
    assert counts["bulk"] == 1


def test_custom_cell_retains_complete_state_fallback(monkeypatch):
    """Unknown evaluators are not assumed to be linear in domain selection."""
    class CustomCell(UnitCell):
        pass

    source = _crystal()
    source.uc_bulk.__class__ = CustomCell
    path = DissolutionPath(source, films=["film"], surface="surface")
    indices = [-6, -3, 1]
    expected = [path.state(index).crystal.F(*_hkl()) for index in indices]
    calls = []
    original = path.state

    def state(index):
        calls.append(index)
        return original(index)

    monkeypatch.setattr(path, "state", state)
    actual = list(path._iter_kinematic_amplitudes(indices, *_hkl()))
    assert calls == indices
    assert_array_equal([value for _, value in actual], expected)


@pytest.mark.parametrize("remove_lower", [False, True])
def test_prefix_lateral_areas_and_reference_frame(remove_lower):
    """Native area conversion remains correct across unlike finite materials."""
    from ..CTRcalc import SXRDCrystal
    from ..CTRfilm import SkellamSurface
    from .test_CTRmultilayer import _parts, _cell
    bulk, components = _parts(lateral=(3., 4., 5., 6.), roughness=.1)
    components.append(SkellamSurface(
        copy.deepcopy(components[-1].unitcell), name="surface"))
    source = SXRDCrystal(
        bulk, *components, stacking=np.arange(1, len(components)+1),
        reference_uc=_cell("reference", lateral=7.),
        reference_rotation=np.diag([-1., -1., 1.]))
    source.setEnergy(18000.)
    source.apply_stacking()
    path = DissolutionPath(
        source, films=["film_0", "film_1"] if remove_lower else ["film_1"],
        surface="surface", exposure_policy="preserve")
    indices = [path.minimum, path.minimum+1, -20, 0, 3]
    expected = [path.state(index).crystal.F(*_hkl()) for index in indices]
    actual = [value for _, value in path._iter_kinematic_amplitudes(indices, *_hkl())]
    assert_allclose(actual, expected, rtol=2e-12, atol=2e-10)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("channels", [("s", "s"), ("s", "p"), ("p", "s"),
                                     ("p", "p")])
@pytest.mark.parametrize("fraction", [0., .2054, 1.])
def test_dwba_reuse_matches_complete_states(monkeypatch, policy, channels,
                                           fraction):
    """Reused fields preserve all mixtures and polarization near exhaustion."""
    from .._CTRdissolution_dwba import DissolutionDWBA
    source = _crystal(interface=True)
    source["film"].unitcell.crystallinity = .63
    model = FiniteFilmDissolution(
        source, films=["film"], surface="surface", optical_reference=policy,
        incoherent_fraction=fraction)
    model.basis[3] = model.path().minimum+.5
    options = dict(polarization_i=channels[0], polarization_f=channels[1],
                   bulk_attenuation=.025)
    with monkeypatch.context() as context:
        context.setattr(DissolutionDWBA, "batch", _complete_dwba_batch)
        context.setattr(DissolutionDWBA, "_evaluate_in_batch", _complete_dwba_oracle)
        expected = model.dwba.evaluate(*_hkl(), **options)
    with model.dwba.batch():
        actual = model.dwba.evaluate(*_hkl(), **options)
        repeated = model.dwba.evaluate(*_hkl(), **options)
    for key in ("reflectivity", "coherent_intensity", "incoherent_intensity",
                "amplitude_prefactor"):
        assert_allclose(getattr(actual, key), getattr(expected, key),
                        rtol=1e-10, atol=1e-28)
        assert_array_equal(getattr(actual, key), getattr(repeated, key))
    assert_array_equal(actual.prepared.reference.profile.values,
                       expected.prepared.reference.profile.values)
    assert_array_equal(actual.prepared.reference.profile.boundaries,
                       expected.prepared.reference.profile.boundaries)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
@pytest.mark.parametrize("fail", [False, True])
def test_dwba_batch_work_and_cleanup(monkeypatch, policy, fail):
    """No repeat packing/field solves; all private caches leave with the batch."""
    from .. import CTRdwba
    from .. import _CTRdissolution_dwba as dissolution
    counts = Counter()
    for name, key in [("_atomic_model_snapshot", "snapshot_builds"),
                      ("_pack_atomic_records", "packing_builds")]:
        original = getattr(CTRdwba.DWBAState, name)

        def counted(self, *args, _original=original, _key=key, **kwargs):
            before = self._cache_counters[_key]
            result = _original(self, *args, **kwargs)
            counts[_key] += self._cache_counters[_key]-before
            return result

        monkeypatch.setattr(CTRdwba.DWBAState, name, counted)
    for obj, name in [(CTRdwba, "solve_electric_field"),
                      (dissolution, "_regrid_model")]:
        original = getattr(obj, name)

        def counted(*args, _original=original, _key=name, **kwargs):
            counts[_key] += 1
            return _original(*args, **kwargs)

        monkeypatch.setattr(obj, name, counted)
    model = FiniteFilmDissolution(
        _crystal(), films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=-5.25), incoherent_fraction=.3)
    retained = []
    try:
        with model.dwba.batch():
            scope = model.dwba._scope
            retained = [scope["average"].dwba, *scope["states"].values()]
            nstates = len(scope["states"])
            first = model.dwba.evaluate(*_hkl())
            assert counts["snapshot_builds"] == nstates+1
            assert counts["packing_builds"] == nstates+int(policy == "shared")
            assert counts["solve_electric_field"] <= (
                4 if policy == "shared" else 2*(nstates+1))
            assert counts["_regrid_model"] == (
                nstates+1 if policy == "shared" else 0)
            counts.clear()
            with model.dwba.batch():
                assert model.dwba._scope is scope
                repeated = model.dwba.evaluate(*_hkl())
                assert_array_equal(first.reflectivity, repeated.reflectivity)
            assert all(value == 0 for value in counts.values()), counts
            if fail:
                model.dwba.evaluate(*_hkl(), bulk_attenuation=-.1)
    except ValueError as error:
        assert fail and "bulk_attenuation" in str(error)
    else:
        assert not fail
    assert model.dwba._scope is None
    for state in retained:
        assert state._batch_depth == 0
        assert state._batch_snapshot is None
        assert not state._batch_packings
        assert not state._batch_profiles
        assert state._batch_native is None
        assert all(value == 0 for key, value in state.cache_info().items()
                   if key.endswith("_size"))
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity, first.reflectivity,
                    rtol=1e-12, atol=1e-28)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_dwba_changed_channels_and_rays_in_batch(monkeypatch, policy):
    """Common-grid reuse does not bind fields to a prior channel or ray block."""
    from .._CTRdissolution_dwba import DissolutionDWBA
    model = FiniteFilmDissolution(
        _crystal(), films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=-5.25), incoherent_fraction=.3)

    def scan():
        values = []
        with model.dwba.batch():
            for pi, pf in [("s", "s"), ("p", "s"), ("p", "p"), ("s", "p")]:
                for shift in [0., .001, 0.]:
                    h, k, ell = _hkl()
                    result = model.dwba.evaluate(
                        h, k, ell+shift, polarization_i=pi, polarization_f=pf,
                        bulk_attenuation=shift)
                    values.append([result.reflectivity, result.coherent_intensity,
                                   result.incoherent_intensity])
        return values

    with monkeypatch.context() as context:
        context.setattr(DissolutionDWBA, "batch", _complete_dwba_batch)
        context.setattr(DissolutionDWBA, "_evaluate_in_batch", _complete_dwba_oracle)
        expected = scan()
    assert_allclose(scan(), expected, rtol=1e-10, atol=1e-28)


@pytest.mark.parametrize("policy", ["shared", "state_specific"])
def test_dwba_explicit_preparation_cache_bypass(policy):
    """An explicit cache=False keeps native preparation caching disabled."""
    model = FiniteFilmDissolution(
        _crystal(), films=["film"], surface="surface", optical_reference=policy,
        profile=RecenteredSkellamProfile(offset=-5.25), incoherent_fraction=.3)
    with model.dwba.batch():
        first = model.dwba.evaluate(*_hkl(), cache=False)
        repeated = model.dwba.evaluate(*_hkl(), cache=False)
        assert_array_equal(first.reflectivity, repeated.reflectivity)
        scope = model.dwba._scope
        for state in [scope["average"].dwba, *scope["states"].values()]:
            assert state.cache_info()["prepared_size"] == 0

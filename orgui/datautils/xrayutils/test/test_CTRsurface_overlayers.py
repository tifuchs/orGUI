"""Finite dissolution retains native overlayers at a common mean front."""

import copy

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from ..CTRcalc import SXRDCrystal
from ..CTRdistributions import (
    PoissonProfile, RecenteredSkellamProfile, bounded_height_states,
)
from ..CTRfilm import (
    Film, FiniteFilmBoundary, LegacySurfaceBoundary, PoissonSurface, SkellamSurface,
)
from ..CTRincoherent import SurfaceHeightDomains
from ..CTRuc import UnitCell, WaterModel
from .test_CTRdissolution import _flat, _hkl
from .test_CTRsurface_boundary import _native


def _cap(kind):
    cell = copy.deepcopy(_flat(2)["flat"].unitcell)
    cell.name = "cap"
    if kind == "tanh":
        return WaterModel(cell.a, np.rad2deg(cell.alpha), "step", name="cap")
    if kind == "water":
        return WaterModel(cell.a, np.rad2deg(cell.alpha), "layered", name="cap")
    if kind == "amorphous":
        cell.crystallinity = 0.
    return (Film(cell, width_layers=2, name="cap")
            if kind in {"film", "amorphous"} else cell)


def _capped(kind, profile, cap_kind="film", **kwargs):
    base = _native(kind, profile, **kwargs)
    result = SXRDCrystal(base.uc_bulk, *base.uc_surface_list, _cap(cap_kind),
                         stacking=np.array([1, 2, 3]))
    result.setEnergy(10000)
    result.dwba.set_ctr_geometry(alpha_i=.03)
    return result


def _oracle(index, anchor, cap_kind):
    base = _flat(6 + int(index))
    cap = _cap(cap_kind)
    cap.stack_on(anchor, anchor, 1)
    result = SXRDCrystal(base.uc_bulk, *base.uc_surface_list, cap,
                         stacking=np.arange(len(base.uc_surface_list)+1))
    result.enable_uc_stacking = False
    result.setEnergy(10000)
    result.dwba.set_ctr_geometry(alpha_i=.03)
    return result


@pytest.mark.parametrize("kind,profile", [
    (PoissonSurface, PoissonProfile(-7, .5)),
    (SkellamSurface, RecenteredSkellamProfile(1, 0, 3, -6)),
])
@pytest.mark.parametrize("cap_kind", ["film", "unitcell", "water", "amorphous"])
def test_common_overlayer_kinematic_states_and_mean_anchor(kind, profile, cap_kind):
    crystal = _capped(kind, profile, cap_kind)
    indices, probabilities, _ = bounded_height_states(profile, -6)
    anchor = float(probabilities @ (3*(6+indices)))
    crystal.apply_stacking()
    assert_allclose(crystal["cap"].pos_absolute, anchor, atol=1e-12)
    path = crystal._surface_state_plan()
    assert path.films == ("film",)  # The cap is above the selected rough boundary.
    expected = []
    for index in indices:
        state = path.state(int(index))
        assert_allclose(state.crystal["cap"].pos_absolute, anchor, atol=1e-12)
        oracle = _oracle(index, anchor, cap_kind)
        amplitude = oracle.F(*_hkl())
        assert_allclose(state.crystal.F(*_hkl()), amplitude, rtol=1e-12, atol=1e-10)
        expected.append(amplitude)
    expected = np.array(expected)
    coherent = probabilities @ expected
    incoherent = probabilities @ abs(expected)**2
    assert_allclose(crystal.F(*_hkl()), coherent, rtol=1e-12, atol=1e-10)
    for fraction in (0, .4, 1):
        model = SurfaceHeightDomains(crystal, surface="surface",
                                     incoherent_fraction=fraction)
        assert_allclose(model.F2(*_hkl()),
                        (1-fraction)*abs(coherent)**2 + fraction*incoherent,
                        rtol=1e-12, atol=1e-9)
    # Geometry and live parameter changes must invalidate the mean anchor.
    crystal["surface"].basis[-1] += 2
    crystal.apply_stacking()
    assert crystal["cap"].pos_absolute > anchor
    restored = SXRDCrystal.fromStr(crystal.toStr(showErrors=False))
    assert_allclose(restored.F(*_hkl()), crystal.F(*_hkl()), atol=1e-9)


@pytest.mark.parametrize("kind,profile", [
    (PoissonSurface, PoissonProfile(-5, 0, -.25)),
    (SkellamSurface, RecenteredSkellamProfile(offset=-5.25)),
])
@pytest.mark.parametrize("optical_reference", ["shared", "state_specific"])
@pytest.mark.parametrize("cap_kind", ["film", "water", "amorphous"])
def test_overlayer_dwba_complete_states(kind, profile, optical_reference, cap_kind):
    from .._CTRincoherent_dwba import _phase_to, _shared_preparation

    crystal = _capped(kind, profile, cap_kind)
    # Independent two-state material oracle: 1/4 empty + 3/4 one remaining plane.
    states = [_oracle(index, 2.25, cap_kind) for index in (-6, -5)]
    if cap_kind in {"film", "amorphous"}:
        for item in (crystal, *states):
            item.weights[-1] = .37
            item.domains[-1] = [(np.diag([-1., -1., 1.]), .6), (np.eye(3), .4)]
    probabilities = (.25, .75)
    common = crystal.dwba.evaluate(*_hkl())
    if cap_kind == "water" and optical_reference == "shared":
        model = SurfaceHeightDomains(crystal, surface="surface",
                                     optical_reference="shared")
        with pytest.raises(NotImplementedError, match="layered optical reference"):
            model.dwba.evaluate(*_hkl())
        return
    amplitudes = []
    for state in states:
        if optical_reference == "shared":
            _, atomic, optical = state.dwba._atomic_model_snapshot()
            prepared = _shared_preparation(
                state.dwba, common.prepared, atomic, optical)
            result = state.dwba.evaluate_prepared(prepared)
            amplitude = (common.unperturbed_amplitude + common.amplitude_prefactor
                         * (result.F_atomic * _phase_to(prepared, common.prepared)
                            - common.F_reference))
        else:
            result = state.dwba.evaluate(*_hkl())
            amplitude = (result.total_amplitude
                         * _phase_to(result.prepared, common.prepared))
        amplitudes.append(amplitude)
    coherent = abs(sum(p*a for p, a in zip(probabilities, amplitudes)))**2
    incoherent = sum(p*abs(a)**2 for p, a in zip(probabilities, amplitudes))
    for fraction in (0, .4, 1):
        model = SurfaceHeightDomains(crystal, surface="surface",
                                     incoherent_fraction=fraction,
                                     optical_reference=optical_reference)
        with model.dwba.batch():
            first = model.dwba.evaluate(*_hkl())
            second = model.dwba.evaluate(*_hkl())
        assert_array_equal(first.reflectivity, second.reflectivity)
        assert_allclose(first.reflectivity, (1-fraction)*coherent + fraction*incoherent,
                        rtol=1e-10, atol=1e-28)
        if optical_reference == "shared" and fraction == 0:
            assert_allclose(first.reflectivity, common.reflectivity, rtol=1e-10)


@pytest.mark.parametrize("kind,profile", [
    (PoissonSurface, PoissonProfile(-2, 0, -.25)),
    (SkellamSurface, RecenteredSkellamProfile(offset=-2.25)),
])
def test_no_exhaustion_reproduces_legacy_common_overlayer(kind, profile):
    native = _capped(kind, profile)
    legacy = _capped(kind, profile, boundary=LegacySurfaceBoundary())
    assert_allclose(native.F(*_hkl()), legacy.F(*_hkl()), atol=1e-10)
    for crystal in (native, legacy):
        crystal.apply_stacking()
    assert_allclose(native["cap"].pos_absolute, legacy["cap"].pos_absolute, atol=1e-12)


def test_weighted_overlayer_is_common_without_repeating_bulk_or_optics(monkeypatch):
    crystal = _capped(SkellamSurface, RecenteredSkellamProfile(offset=-5.25))
    crystal.weights[-1] = .37
    crystal.domains[-1] = [(np.diag([-1., -1., 1.]), .6), (np.eye(3), .4)]
    expected = []
    for index in (-6, -5):
        state = _oracle(index, 2.25, "film")
        state.weights[-1] = crystal.weights[-1]
        state.domains[-1] = copy.deepcopy(crystal.domains[-1])
        expected.append(state.F(*_hkl()))
    calls = []
    original = UnitCell.F_bulk

    def counted(self, *args, **kwargs):
        calls.append(self)
        return original(self, *args, **kwargs)

    def forbidden(*args, **kwargs):
        pytest.fail("Kinematic overlayer evaluation requested optics")

    monkeypatch.setattr(UnitCell, "F_bulk", counted)
    monkeypatch.setattr(UnitCell, "optical_profile", forbidden)
    model = SurfaceHeightDomains(crystal, surface="surface", incoherent_fraction=.4)
    actual = model.F2(*_hkl())
    assert len(calls) == 1
    coherent = abs(.25*expected[0] + .75*expected[1])**2
    incoherent = .25*abs(expected[0])**2 + .75*abs(expected[1])**2
    assert_allclose(actual, .6*coherent + .4*incoherent, rtol=1e-12, atol=1e-9)


def test_multiple_overlayers_keep_native_relative_placement():
    crystal = _capped(SkellamSurface, RecenteredSkellamProfile(offset=-5.25))
    upper = _cap("unitcell")
    upper.name = "upper"
    crystal = SXRDCrystal(crystal.uc_bulk, *crystal.uc_surface_list, upper,
                         stacking=np.arange(4))
    crystal.apply_stacking()
    assert_allclose(crystal["cap"].pos_absolute, 2.25)
    assert_allclose(crystal["upper"].pos_absolute, 8.25)
    for index in (-6, -5, 0, 1):
        state = crystal.surface_height_state(surface="surface", index=index).crystal
        assert_allclose(state["cap"].pos_absolute, 2.25)
        assert_allclose(state["upper"].pos_absolute, 8.25)


def test_analytic_water_reference_survives_complete_exhaustion():
    crystal = _capped(SkellamSurface, RecenteredSkellamProfile(offset=-20), "tanh")
    expected = _oracle(-6, 0, "tanh").dwba.evaluate(*_hkl())
    actual = crystal.dwba.evaluate(*_hkl())
    assert actual.prepared.reference.continuous is not None
    assert_allclose(actual.reflectivity, expected.reflectivity, rtol=1e-11)
    model = SurfaceHeightDomains(crystal, surface="surface", incoherent_fraction=.4,
                                 optical_reference="state_specific")
    assert_allclose(model.dwba.evaluate(*_hkl()).reflectivity, expected.reflectivity,
                    rtol=1e-11)


def test_overlayer_cannot_be_selected_as_underlying_removable_material():
    crystal = _capped(SkellamSurface, RecenteredSkellamProfile(offset=-5.25))
    crystal["surface"].boundary = FiniteFilmBoundary(("film", "cap"))
    with pytest.raises(ValueError, match="below the target surface"):
        crystal.F(*_hkl())

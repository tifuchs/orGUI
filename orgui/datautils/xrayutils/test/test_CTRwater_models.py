"""Consolidated physical water families and pre-consolidation file migration."""

import copy
from contextlib import nullcontext

import numpy as np
import pytest
from scipy.integrate import simpson

from ..CTRuc import WaterModel
from ..CTRcalc import SXRDCrystal
from .test_CTRcrystallinity import _cell
from .test_CTRwater import _water


# Complex F(0,0,.37) recorded before consolidation at 10 keV, with a=3,3,4 A
# and parameters .2,.7,.15,.25,.8. Retain independent migration baselines.
_LEGACY_AMPLITUDES = {
    "step": -2.1398289926568177 + 4.24360970726116j,
    "tanh": -2.1991539147675976 + 4.361856578376978j,
    "erf": -2.197565163069214 + 4.358705408738295j,
    "layered": 0.9955682891930647 + 4.540003393186669j,
    "1layer": 2.452589388108655 + 9.983685772378703j,
    "1layer_O": 2.4616454124650353 + 10.042476830499282j,
    "layered_O": 0.7991970959253962 + 3.653341798581633j,
    "1layer_s_O": 4.3841205017615525 + 7.485622002506086j,
}


@pytest.mark.parametrize("kind", _LEGACY_AMPLITUDES)
def test_legacy_files_preserve_amplitude_and_roundtrip(kind, tmp_path):
    """Read old names, then write only canonical families with explicit settings."""
    old = f"{kind}\n3 3 4 90 90 90\n.2 .7 .15 .25 .8\n"
    warning = (
        pytest.warns(DeprecationWarning)
        if kind in WaterModel.legacyTypes
        else nullcontext()
    )
    with warning:
        water = WaterModel.fromStr(old)
    water.setEnergy(10000.0)
    expected = _LEGACY_AMPLITUDES[kind]
    np.testing.assert_allclose(water.F_uc(0.0, 0.0, 0.37), expected, rtol=2e-13)
    assert water.type in {"step", "layered"}
    path = tmp_path / "water.wat"
    water.writeWATfile(path)
    restored = WaterModel.fromWATfile(path)
    restored.setEnergy(10000.0)
    assert restored._settings() == water._settings()
    np.testing.assert_allclose(restored.F_uc(0.0, 0.0, 0.37), expected, rtol=2e-13)
    np.testing.assert_allclose(restored.pw, water.pw)


def test_untyped_four_parameter_step_keeps_adsorption_disabled():
    """Old files default to step and A0=1, without adding an adsorbed layer."""
    water = WaterModel.fromStr("3 3 4 90 90 90\n.2 .7 .15 .25\n")
    water.setEnergy(10000.0)
    assert water.type == "step" and water.basis[4] == 1.0
    assert water.optical_reference == "sharp" and not water.adsorption
    before = water.F_uc(0.0, 0.0, 0.37)
    water.basis[4] = 4.0
    assert water.F_uc(0.0, 0.0, 0.37) == before
    water.adsorption = True
    assert abs(water.F_uc(0.0, 0.0, 0.37) - before) > 1.0


def test_step_adsorption_normalization_and_one_layer_alias():
    """A0 counts molecules per lateral cell; the old one-layer physics survives."""
    water = _water("step_sharp")
    empty = copy.deepcopy(water)
    empty.basis[4] = 0.0
    z = np.linspace(-15.0, 25.0, 16001)
    excess = water.zDensity_G(z, 0.0, 0.0) - empty.zDensity_G(z, 0.0, 0.0)
    forward_factor = water.bulk_density / water.pw
    np.testing.assert_allclose(
        simpson(excess, x=z) * water.uc_area,
        water.basis[4] * forward_factor,
        rtol=2e-12,
        atol=1e-12,
    )
    with pytest.warns(DeprecationWarning):
        legacy = _water("1layer")
    np.testing.assert_allclose(water.F_uc(0.0, 0.0, 0.37), legacy.F_uc(0.0, 0.0, 0.37))
    np.testing.assert_allclose(
        water.zDensity_G(z, 0.0, 0.0), legacy.zDensity_G(z, 0.0, 0.0)
    )


@pytest.mark.parametrize("a0", [0.0, 0.8])
def test_optical_reference_does_not_change_physical_density_or_factors(a0):
    """Tanh/erf/sharp choose optics, not a different molecule or density profile."""
    water = _water("step")
    water.basis[4] = a0
    z, ell = np.linspace(-8.0, 15.0, 101), np.array([0.1, 0.3, 1.2])
    density, amplitude = water.zDensity_G(z, 0.0, 0.0), water.F_uc(0.0, 0.0, ell)
    fingerprints = []
    for reference in ("sharp", "tanh", "erf"):
        water.optical_reference = reference
        np.testing.assert_array_equal(water.zDensity_G(z, 0.0, 0.0), density)
        np.testing.assert_array_equal(water.F_uc(0.0, 0.0, ell), amplitude)
        crystal = SXRDCrystal(_cell(), water)
        p = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.03, 0.03)
        fingerprints.append(p.reference.fingerprint)
        assert np.isfinite(crystal.dwba.evaluate_prepared(p).total_amplitude).all()
    assert len(set(fingerprints)) == 3


def test_reference_and_live_density_settings_invalidate_caches():
    """Physical convention edits cannot reuse stale frozen-water integrals."""
    water = _water("step_erf")
    crystal = SXRDCrystal(_cell(), water)
    prepared = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.03, 0.03)
    old = crystal.dwba.evaluate_prepared(prepared).F_h
    water.adsorption = False
    changed = crystal.dwba.evaluate_prepared(prepared).F_h
    assert abs(changed - old) > 1e-5
    fresh = crystal.dwba.prepare_from_glancing(0.0, 0.0, 0.03, 0.03)
    assert fresh.reference.fingerprint != prepared.reference.fingerprint
    water.adsorption = True
    np.testing.assert_array_equal(crystal.dwba.evaluate_prepared(prepared).F_h, old)


def test_legacy_hdf5_settings_keep_inactive_step_a0(tmp_path):
    """Unversioned fit state must not turn ignored A0 into physical adsorption."""
    from silx.io.dictdump import dicttoh5, h5todict

    water = _water("step")
    water.addFitParameter("rho_0", limits=(0.0, 2.0), name="old_a0")
    water.getStartParamAndLimits()
    state = water.parametersToDict()
    del state["water_settings"]
    path = str(tmp_path / "legacy.h5")
    dicttoh5(state, path)
    water.parametersFromDict(h5todict(path))
    assert not water.adsorption and water.optical_reference == "sharp"
    before = water.F_uc(0.0, 0.0, 0.37)
    water.setFitParameters([1.5])
    assert water.F_uc(0.0, 0.0, 0.37) == before


def test_invalid_reference_does_not_silently_change_layered_physics():
    """Reject unsupported optics choices instead of silently ignoring them."""
    assert WaterModel.supportedTypes == ["step", "layered"]
    with pytest.raises(ValueError, match="sharp"):
        WaterModel(
            [3.0, 3.0, 4.0], [90.0, 90.0, 90.0], "layered", optical_reference="erf"
        )
    water = _water()
    settings = water._settings()
    with pytest.raises(ValueError, match="version"):
        water._restore_settings({**settings, "version": 2.5})
    with pytest.raises(ValueError, match="boolean"):
        water._restore_settings({**settings, "adsorption": "false"})


@pytest.mark.parametrize("kind", ["step", "1layer_O", "layered_O", "1layer_s_O"])
def test_versioned_fit_and_error_text_roundtrip(kind, tmp_path):
    """Canonical adsorption and imported oxygen conventions survive fit storage."""
    from silx.io.dictdump import dicttoh5, h5todict

    with (
        pytest.warns(DeprecationWarning)
        if kind in WaterModel.legacyTypes
        else nullcontext()
    ):
        water = _water(kind)
    water.addFitParameter("A0", limits=(0.0, 2.0), name="adsorbed_water")
    water.getStartParamAndLimits()
    water.setFitParameters([1.2])
    water.setFitErrors([0.03])
    water.pw = 0.031
    path = str(tmp_path / "settings.h5")
    dicttoh5(water.parametersToDict(), path)
    restored = WaterModel.fromStr(water.toStr().removeprefix("return\n"))
    restored.parametersFromDict(h5todict(path))
    restored.setEnergy(10000.0)
    assert restored._settings() == water._settings()
    assert restored.pw == water.pw
    np.testing.assert_allclose(restored.getFitErrors(), [0.03])
    np.testing.assert_allclose(
        restored.F_uc(0.0, 0.0, 0.37), water.F_uc(0.0, 0.0, 0.37)
    )

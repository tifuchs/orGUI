"""Coherent order changes contrast without changing material or optics."""

import copy
import pickle

import numpy as np
import pytest

from .. import CTRcalc, CTRfilm, CTRuc
from ..CTRutil import generate_surface_termination_cells
from ._poisson_oracle import flat_height_crystal, poisson_crystal


def _cell(name="cell", crystallinity=1.0):
    cell = CTRuc.UnitCell(
        [3.0, 3.0, 4.0], [90.0, 90.0, 90.0],
        name=name, crystallinity=crystallinity,
    )
    cell.addAtom("C", [0.13, 0.17, 0.0], 0.11, 0.19, 1.0, layer=1)
    cell.addAtom("O", [0.31, 0.29, 0.5], 0.14, 0.23, 1.0, layer=2)
    cell.layerpos = {1.0: 0.0, 2.0: 0.5}
    cell.setEnergy(10000.0)
    return cell


def _hkl():
    return np.array([0., 1., 1.]), np.zeros(3), np.array([0.27, 0.43, 0.81])


def _prepare(crystal):
    return crystal.dwba.prepare_from_glancing(
        np.array([0.0, 0.6]), 0.0,
        np.array([0.04, 0.04]), np.array([0.04, 0.05]),
    )


@pytest.mark.parametrize("value", [-0.1, 1.1, np.nan, np.inf, -np.inf, [0.5]])
def test_invalid_crystallinity(value):
    """Invalid amplitudes must never silently enter either scattering path."""
    with pytest.raises(ValueError, match="crystallinity"):
        _cell(crystallinity=value)


@pytest.mark.parametrize("backend", ["numpy", "cpp", "numba"])
@pytest.mark.parametrize("repeat", [(0., 0., 1.), (0.2, -0.1, 1.)])
def test_kinematic_amplitudes_and_raw_density(backend, repeat, monkeypatch):
    """Every motif/bulk backend scales once, leaving physical density intact."""
    if backend == "cpp" and not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("native CTR extension unavailable")
    if backend == "numba":
        pytest.importorskip("numba")
    monkeypatch.setattr(CTRuc, "CTR_ACCEL_BACKEND", backend)
    cell = _cell()
    cell.stacking_translation = repeat
    methods = (cell.F_uc, cell.F_uc_bulk, cell.F_uc_bulk_direct, cell.F_bulk)
    initial = [method(*_hkl()) for method in methods]
    profile = cell.optical_profile()
    z = np.linspace(-2., 6., 101)
    density = cell.zDensity_G(z, 0., 0.)
    for value in [0.37, 0.0, 1.0]:
        cell.crystallinity = value
        for method, amplitude in zip(methods, initial):
            np.testing.assert_allclose(method(*_hkl()), value * amplitude)
        np.testing.assert_array_equal(cell.optical_profile(), profile)
        np.testing.assert_array_equal(cell.zDensity_G(z, 0., 0.), density)


@pytest.mark.parametrize(
    "clone", [copy.deepcopy, lambda x: pickle.loads(pickle.dumps(x))],
)
def test_generated_layers_follow_copied_source(clone):
    """Copies retain internal ownership without referring back to the original."""
    source = _cell()
    film = CTRfilm.Film(source, width_layers=4, name="film")
    original = film.F_uc(*_hkl())
    copied = clone(film)
    copied.unitcell.crystallinity = 0.23
    np.testing.assert_allclose(copied.F_uc(*_hkl()), 0.23 * original)
    np.testing.assert_allclose(film.F_uc(*_hkl()), original)
    copied.unitcell.crystallinity = 0.0
    np.testing.assert_array_equal(copied.F_uc(*_hkl()), 0.0j)


def test_cell_transformations_preserve_value():
    """Independent transformed cells retain the material's initial order."""
    cell = _cell(crystallinity=0.4)
    cell.addCrystallinityParameter()
    transformed = [
        cell.translate_layered([0., 0., 1.]),
        cell.affine_layer_transform([0., 0., 1.]),
        cell.as_surface_termination(1),
        cell.supercell((1, 1, 2)),
        *generate_surface_termination_cells(cell, cell.layer_cycle).values(),
    ]
    cell.crystallinity = 0.9
    for item in transformed:
        assert item.crystallinity == 0.4
        item.getInitialParameters()


@pytest.mark.parametrize("backend", ["numpy", "cpp"])
def test_zero_order_bulk_has_no_bragg_pole(backend, monkeypatch):
    """Zero coherent contrast must not yield zero times a divergent sum."""
    monkeypatch.setattr(CTRuc, "CTR_ACCEL_BACKEND", backend)
    cell = _cell(crystallinity=0.)
    np.testing.assert_array_equal(cell.F_bulk(np.zeros(1), np.zeros(1),
                                             np.zeros(1)), 0.j)


@pytest.mark.parametrize("explicit", [False, True])
def test_full_model_text_preserves_source_order(tmp_path, explicit):
    """Film and surface values survive their nested native text sections."""
    crystal, surface = poisson_crystal(CTRfilm.PoissonProfile(-0.7, alpha=0.5))
    crystal.apply_stacking()
    if explicit:
        surface.set_ucs(copy.deepcopy(surface.termination_cells))
    crystal.uc_bulk.crystallinity = 0.9
    crystal.uc_surface_list[0].unitcell.crystallinity = 0.65
    for index, cell in enumerate(surface._owned_unitcells()):
        cell.crystallinity = 0.2 + 0.1 * index
    crystal.setEnergy(10000.)
    initial = crystal.F(*_hkl())
    path = str(tmp_path / "crystallinity.xtal")
    crystal.toFile(path)
    restored = CTRcalc.SXRDCrystal.fromFile(path)
    restored.setEnergy(10000.)
    assert restored.uc_bulk.crystallinity == 0.9
    assert restored.uc_surface_list[0].unitcell.crystallinity == 0.65
    np.testing.assert_allclose(restored.F(*_hkl()), initial, atol=1e-11)
    before = crystal.dwba.evaluate_prepared(_prepare(crystal))
    after = restored.dwba.evaluate_prepared(_prepare(restored))
    np.testing.assert_allclose(after.F_h, before.F_h, atol=1e-11)


def test_fit_parameters_and_legacy_roundtrip(tmp_path):
    """The scalar fits beside atomic parameters and survives text/HDF5 saves."""
    from silx.io.dictdump import dicttoh5, h5todict

    cell = _cell(crystallinity=0.8)
    cell.addFitParameter(([0], "x"), limits=(0., 1.), name="x")
    cell.addRelParameter(([1], "z"), [1.], limits=(-0.2, 0.2), name="z_shift")
    cell.addCrystallinityParameter(name="order", limits=(0.1, 1.))
    np.testing.assert_allclose(cell.getInitialParameters(), [0.13, 0., 0.8])
    cell.setFitParameters([0.2, 0.05, 0.35])
    cell.setFitErrors([0.01, 0.02, 0.03])
    np.testing.assert_allclose(cell.getFitErrors(), [0.01, 0.02, 0.03])
    assert cell.fitparnames == ["x", "z_shift", "order"]
    assert len(cell.priors) == 3
    assert "crystallinity=0.35" in cell.fitparameterList()
    cell.setLimits([(0., 1.), (-0.2, 0.2), (0.2, 0.9)])
    with pytest.raises(ValueError, match="Crystallinity limits"):
        cell.setLimits([(0., 1.), (-0.2, 0.2), (-0.1, 1.)])

    # A direct edit must also be saved, even after a previous fit.
    cell.crystallinity = 0.45
    settings = cell.parametersToDict()
    path = str(tmp_path / "fit.h5")
    dicttoh5(settings, path)
    restored = CTRuc.UnitCell.fromStr(cell.toStr())
    assert restored.crystallinity == 0.45
    restored.parametersFromDict(h5todict(path))
    assert restored.crystallinity == 0.45
    np.testing.assert_allclose(restored.getFitErrors(), [0.01, 0.02, 0.03])
    np.testing.assert_allclose(restored.getInitialParameters(), [0.2, 0.05, 0.45])
    assert restored.fitparnames == cell.fitparnames

    legacy_text = "\n".join(
        line for line in cell.toStr().splitlines()
        if not line.startswith("crystallinity:")
    )
    assert CTRuc.UnitCell.fromStr(legacy_text).crystallinity == 1.0
    settings.pop("crystallinity")
    settings.pop("crystallinity_value")
    restored.parametersFromDict(settings)
    assert restored.crystallinity == 1.0
    restored.clearParameters()
    assert restored.fitparnames == []
    assert restored.crystallinity == 1.0


def test_dwba_live_order_preserves_prepared_optics():
    """Both specular and non-specular contrast scale on the same preparation."""
    cell = _cell()
    crystal = CTRcalc.SXRDCrystal(cell)
    prepared = _prepare(crystal)
    before = crystal.dwba.evaluate_prepared(prepared)
    profile = crystal.optical_profile().copy()
    for value in [0.35, 0., 1.]:
        cell.crystallinity = value
        result = crystal.dwba.evaluate_prepared(prepared)
        np.testing.assert_allclose(result.F_h, value * before.F_h, atol=1e-12)
        np.testing.assert_array_equal(result.F_reference, before.F_reference)
        np.testing.assert_array_equal(result.unperturbed_amplitude,
                                      before.unperturbed_amplitude)
        np.testing.assert_allclose(result.F_atomic - result.F_reference, result.F_h)
        np.testing.assert_allclose(sum(c.F_h for c in result.contributions), result.F_h)
        np.testing.assert_array_equal(crystal.optical_profile(), profile)
        assert _prepare(crystal) is prepared
        if value == 0.:
            np.testing.assert_array_equal(result.F_h, 0.j)
            np.testing.assert_array_equal(result.total_amplitude,
                                          before.unperturbed_amplitude)


def test_interface_removal_retains_underlying_order():
    """Equal structures with independent order must remain separate records."""
    bulk, top, bottom = (_cell(name) for name in ("bulk", "top", "bottom"))
    interface = CTRfilm.EpitaxyInterface(
        top, bottom, profile=CTRfilm.SkellamProfile(0., 0.),
        fixed_ucs=2, name="interface",
    )
    crystal = CTRcalc.SXRDCrystal(bulk, interface, stacking=np.array([1]))
    prepared = _prepare(crystal)
    initial = crystal.dwba.evaluate_prepared(prepared)
    removed = [c for c in initial.contributions if c.role == "interface_reference"]
    lower = [c for c in initial.contributions if c.role == "interface_bottom"]
    assert removed and lower
    np.testing.assert_allclose(sum(c.F_h for c in removed + lower), 0., atol=1e-12)
    original_kinematic = crystal.F(*_hkl())
    lower_kinematic = sum(c.F_uc(*_hkl()) for c in interface.bottom_layers)

    bottom.crystallinity = 0.
    result = crystal.dwba.evaluate_prepared(prepared)
    np.testing.assert_allclose(result.F_h,
                               initial.F_h - sum(c.F_h for c in lower), atol=1e-12)
    for previous, current in zip(initial.contributions, result.contributions):
        if current.role == "interface_reference":
            np.testing.assert_array_equal(current.F_h, previous.F_h)
        elif current.role == "interface_bottom":
            np.testing.assert_array_equal(current.F_h, 0.j)
    np.testing.assert_allclose(crystal.F(*_hkl()),
                               original_kinematic - lower_kinematic, atol=1e-12)
    assert _prepare(crystal) is prepared


@pytest.mark.parametrize("offset", [-1., 1.])
@pytest.mark.parametrize("explicit", [False, True])
def test_poisson_material_removal_and_exposed_order(offset, explicit):
    """An unordered exposed layer still replaces the ordered film beneath it."""
    profile = CTRfilm.PoissonProfile(0., alpha=0., offset=offset)
    crystal, surface = poisson_crystal(profile)
    if explicit:
        crystal.apply_stacking()
        surface.set_ucs(copy.deepcopy(surface.termination_cells))
    crystal.setEnergy(10000.)
    # Change values only after preparation, including generated terminations.
    prepared = _prepare(crystal)
    baseline = crystal.dwba.evaluate_prepared(prepared)
    original_surface = [c for c in baseline.contributions
                        if c.role == "surface_termination"]
    film = crystal.uc_surface_list[0]
    film.unitcell.crystallinity = 0.4
    for cell in surface._owned_unitcells():
        cell.crystallinity = 0.

    # Each exposed termination is a full two-layer cell. Flat state
    # n=offset-1 has width w_base+offset; removing that entire exposed slab
    # leaves two fewer crystalline layers, not just its outermost plane.
    expected = flat_height_crystal(int(offset) - 3)
    expected.setEnergy(10000.)
    expected.uc_surface_list[0].unitcell.crystallinity = 0.4
    np.testing.assert_allclose(crystal.F(*_hkl()), expected.F(*_hkl()), atol=1e-11)

    result = crystal.dwba.evaluate_prepared(prepared)
    for previous, current in zip(baseline.contributions, result.contributions):
        factor = (1. if current.role == "bulk" else
                  0. if current.role == "surface_termination" else 0.4)
        np.testing.assert_allclose(current.F_h, factor * previous.F_h, atol=1e-12)
    assert original_surface
    assert _prepare(crystal) is prepared


def test_crystal_fit_updates_generated_dwba_records():
    """Public crystal fit updates reach layers and reuse cached wavefields."""
    source = _cell("film", crystallinity=0.8)
    source.addCrystallinityParameter(name="film_order")
    film = CTRfilm.Film(source, width_layers=4, name="film")
    crystal = CTRcalc.SXRDCrystal(_cell("bulk"), film, stacking=np.array([1]))
    crystal = copy.deepcopy(crystal)
    prepared = _prepare(crystal)
    before = crystal.dwba.evaluate_prepared(prepared)
    np.testing.assert_allclose(crystal.getInitialParameters(), [0.8])
    crystal.setParameters([0.2])
    with crystal.dwba.batch():
        result = crystal.dwba.evaluate_prepared(prepared)
        again = crystal.dwba.evaluate_prepared(prepared)
    for initial, current in zip(before.contributions, result.contributions):
        factor = 1. if current.role == "bulk" else 0.25
        np.testing.assert_allclose(current.F_h, factor * initial.F_h, atol=1e-12)
    np.testing.assert_array_equal(result.F_h, again.F_h)
    assert _prepare(crystal) is prepared

"""Cartesian geometry and scattering oracles for nonvertical crystal repeats."""

import copy
import pickle

import numpy as np
import pytest

from .. import CTRfilm, CTRsymmetry, CTRuc


def _cell(shift=(0.23, -0.17, 1)):
    cell = CTRuc.UnitCell((4, 5, 6), (90, 90, 76), name="shifted",
                          stacking_translation=shift)
    cell.addAtom("Ti", (0.12, 0.22, 0.15), 0.4, 0.4, 0.8, layer=1)
    cell.addAtom("O", (0.34, 0.17, 0.65), 0.3, 0.3, 1, layer=2)
    cell.layerpos = {1: 0.0, 2: 0.4}
    cell.setEnergy(15000)
    return cell


def _hkl():
    return np.array([[0.31, 1.2, -0.72], [0.41, -0.2, 0.6], [0.37, 1.43, 2.12]])


@pytest.mark.parametrize("value", [(1, 2, 0), (1, 2), (np.nan, 0, 1)])
def test_invalid_repeat_rejected(value):
    """Invalid repeats fail at the state boundary."""
    with pytest.raises(ValueError, match="stacking_translation"):
        _cell(value)


def test_repeat_persistence_and_layer_copies():
    """A saved/copy/split cell retains a repeat with phase-relevant precision."""
    cell = _cell((np.sqrt(2) / 5, -np.pi / 9, 1))
    for restored in [copy.deepcopy(cell), pickle.loads(pickle.dumps(cell)),
                     CTRuc.UnitCell.fromStr(cell.toStr()),
                     *cell.split_in_layers().values()]:
        assert restored.stacking_translation == cell.stacking_translation
    with pytest.raises(ValueError, match="ROD export"):
        cell.parameterStrRod()
    cell.coherentDomainMatrix[0][:, 3] = [np.sqrt(2) / 3, -np.pi / 9, 0.123456789]
    restored = CTRuc.UnitCell.fromStr(cell.toStr())
    restored.setEnergy(15000)
    np.testing.assert_allclose(restored.F_bulk(*_hkl(), atten=0.1),
                               cell.F_bulk(*_hkl(), atten=0.1), rtol=2e-12)


@pytest.mark.parametrize("backend", ["numpy", "numba", "cpp"])
@pytest.mark.parametrize("with_reference", [False, True])
def test_bulk_matches_explicit_cartesian_repeat_sum(backend, with_reference):
    """The infinite sum agrees with a converged finite physical stack."""
    old = CTRuc.CTR_ACCEL_BACKEND
    try:
        try:
            CTRuc.set_accel_backend(backend)
        except (ValueError, ImportError, RuntimeError):
            pytest.skip(f"{backend} unavailable")
        cell = _cell()
        if with_reference:
            reference = CTRuc.UnitCell((4.3, 5.2, 7), (90, 90, 83))
            cell.setReferenceUnitCell(reference)
        hkl = _hkl()
        eta = 0.25
        repeat_eta = eta * abs(cell.refHKLTransform[2, 2])
        q = cell.B_mat @ cell.refHKLTransform @ hkl
        phase = q.T @ cell.repeat_cartesian
        repeats = np.arange(180)
        expected = cell.F_uc_bulk(*hkl, repeat_eta) * np.sum(
            np.exp(-repeats[:, None] * (repeat_eta + 1j * phase)), axis=0
        )
        np.testing.assert_allclose(cell.F_bulk(*hkl, eta), expected, rtol=2e-12)
    finally:
        CTRuc.set_accel_backend(old)


def test_projected_surface_matches_original_physical_motif():
    """Projection preserves occupancy, B, parent positions and bulk phase."""
    pm = pytest.importorskip("pymatgen.core")
    pytest.importorskip("pyxtal")
    source = pm.Structure(pm.Lattice.orthorhombic(5.7233, 7.9147, 5.4543),
                          [{"Tb": 0.96}], [[0, 0, 0]],
                          site_properties={"B_iso": [0.6]})
    native = CTRsymmetry.surface_terminations(source, miller=(1, 0, 1),
                                              cell_mode="lattice")[0]
    surface = CTRsymmetry.surface_terminations(source, miller=(1, 0, 1),
                                               cell_mode="surface")[0]
    a = native.build_unitcell("native")
    b = surface.build_unitcell("surface")
    np.testing.assert_allclose(np.rad2deg(b.alpha[:2]), 90, atol=1e-10)
    np.testing.assert_allclose(b.stacking_translation, [0, 0.47594794115389494, 1],
                               atol=1e-12)
    np.testing.assert_allclose(b.a[2], 3.948446313120785)
    ap = (native.surface_spec.transform @ a.basis[:, 1:4].T).T
    bp = (surface.surface_spec.transform @ b.basis[:, 1:4].T).T
    np.testing.assert_allclose(ap, bp, atol=1e-12)
    np.testing.assert_allclose(a.basis[:, 4:7], b.basis[:, 4:7])
    a.setEnergy(15000)
    b.setEnergy(15000)
    parent_hkl = _hkl()
    np.testing.assert_allclose(
        a.F_bulk(*(native.surface_spec.transform.T @ parent_hkl), atten=0.1),
        b.F_bulk(*(surface.surface_spec.transform.T @ parent_hkl), atten=0.1),
        rtol=2e-12, atol=1e-10,
    )
    restored = CTRuc.UnitCell.fromStr(b.toStr())
    np.testing.assert_allclose(restored.symmetry_metadata.surface_spec.bulk_transform,
                               native.surface_spec.transform)
    restored.setEnergy(15000)
    np.testing.assert_allclose(restored.F_bulk(*_hkl(), atten=0.1),
                               b.F_bulk(*_hkl(), atten=0.1), rtol=2e-12)


def test_whole_repeat_translation_and_supercell():
    """Translations and supercells place copies at the physical repeat vector."""
    cell = _cell()
    shifted = cell.translate_layered((0, 0, -2))
    np.testing.assert_allclose(shifted.basis[:, 1:4],
                               cell.basis[:, 1:4] - cell.stacking_translation)
    expanded = cell.supercell((2, 1, 3))
    expected = np.concatenate([
        cell.basis[:, 1:4] + np.array([i, 0, 0]) + n*np.array(cell.stacking_translation)
        for n in range(3) for i in range(2)
    ]) @ cell.R_mat.T
    actual = expanded.basis[:, 1:4] @ expanded.R_mat.T
    # The supercell groups rows by layer; compare Cartesian sets.
    def sort(v):
        """Sort Cartesian rows for a representation-independent comparison."""
        return v[np.lexsort(v.T[::-1])]
    np.testing.assert_allclose(sort(actual), sort(expected), atol=1e-12)
    np.testing.assert_allclose(expanded.stacking_translation, (0.345, -0.51, 1))


def test_film_repeat_registry():
    """Film layer copies include lateral repeat offsets without normal strain."""
    cell = _cell()
    film = CTRfilm.Film(cell)
    film.basis[0] = 6
    film.createLayers()
    expected = cell.F_uc(*_hkl()) * sum(
        np.exp(2j*np.pi*n*(np.array(cell.stacking_translation) @ _hkl()))
        for n in range(3)
    )
    np.testing.assert_allclose(film.F_uc(*_hkl()), expected, rtol=2e-12, atol=1e-10)


def test_film_continuation_preserves_registry_and_normal_gap():
    """Splitting a stack preserves positions; an added gap is purely normal."""
    first = CTRfilm.Film(_cell(), layer_behavior="select")
    first.basis[0] = 3
    first.createLayers()
    second = CTRfilm.Film(_cell(), layer_behavior="select")
    second.basis[0] = 3
    second.stack_on(first.height_absolute, first.height_absolute,
                    below_state=first.layer_state)
    whole = CTRfilm.Film(_cell())
    whole.basis[0] = 6
    np.testing.assert_allclose(first.F_uc(*_hkl()) + second.F_uc(*_hkl()),
                               whole.F_uc(*_hkl()), atol=1e-10)
    before = second.F_uc(*_hkl())
    second.pos_absolute += 1.7
    np.testing.assert_allclose(second.F_uc(*_hkl()),
                               before * np.exp(2j*np.pi*_hkl()[2]*1.7/6), atol=1e-10)


def test_changed_repeat_invalidates_film_geometry():
    """A changed source repeat cannot leave generated film domains stale."""
    cell = _cell()
    film = CTRfilm.Film(cell)
    film.basis[0] = 4
    film.F_uc(*_hkl())
    cell.stacking_translation = (-0.13, 0.42, 1)
    expected = cell.F_uc(*_hkl()) * (1 + np.exp(
        2j*np.pi*(np.array(cell.stacking_translation) @ _hkl())))
    np.testing.assert_allclose(film.F_uc(*_hkl()), expected, atol=1e-10)


@pytest.mark.parametrize("growth", [-1.1, 0, 1.3])
def test_poisson_shifted_stack_matches_flat_height_oracle(monkeypatch, growth):
    """Height-state corrections retain off-specular lateral phase."""
    from . import _poisson_oracle as oracle
    from ..CTRdistributions import PoissonProfile

    original = oracle.layered_cell

    def shifted(*args, **kwargs):
        cell = original(*args, **kwargs)
        cell.stacking_translation = (0.23, -0.17, 1)
        cell.basis[:, 1:3] += cell.basis[:, 3, None] * [0.23, -0.17]
        cell.basis_0 = cell.basis.copy()
        return cell

    monkeypatch.setattr(oracle, "layered_cell", shifted)
    profile = PoissonProfile(growth, alpha=0.3, tail_probability=1e-12)
    width = oracle.minimum_base_width(profile)
    crystal, surface = oracle.poisson_crystal(profile, w_base=width)
    crystal.apply_stacking()
    expected = oracle.coherent_reference(profile, *_hkl(), w_base=width)
    np.testing.assert_allclose(crystal.F(*_hkl()), expected, rtol=2e-10, atol=1e-9)
    common = oracle.flat_height_crystal(-1, w_base=width).F(*_hkl())
    for state in surface.flat_domain_corrections(*_hkl()).iter_states():
        reference = oracle.flat_height_crystal(state.layer_number, w_base=width)
        np.testing.assert_allclose(common + state.amplitude, reference.F(*_hkl()),
                                   rtol=2e-10, atol=1e-9)


@pytest.mark.parametrize("coupling", [0, 0.4, 1])
def test_equal_material_interface_cancels_with_shift(coupling):
    """Replacing a material by itself cannot create off-specular scattering."""
    from .. import CTRcalc
    from ..CTRdistributions import SkellamProfile

    source = _cell()
    source.layerpos = {1: 0, 2: 0.5}
    bulk = source.translate_layered((0, 0, -2))
    interface = CTRfilm.EpitaxyInterface(copy.deepcopy(source), copy.deepcopy(source),
        profile=SkellamProfile(0.4, tail_probability=1e-12))
    interface.basis[2] = coupling
    film = CTRfilm.Film(copy.deepcopy(source))
    film.basis[0] = 16
    crystal = CTRcalc.SXRDCrystal(bulk, interface, film, stacking=np.array([1, 2]))
    crystal.setEnergy(15000)
    crystal.apply_stacking()
    reference = CTRfilm.Film(copy.deepcopy(source))
    reference.basis[0] = 16
    np.testing.assert_allclose(interface.F_uc(*_hkl()) + film.F_uc(*_hkl()),
                               reference.F_uc(*_hkl()), rtol=2e-10, atol=1e-9)


@pytest.mark.parametrize("h,k", [(0, 0), (0.31, -0.42)])
def test_density_bulk_has_repeat_fourier_phase(h, k):
    """Real-space density follows the same lateral registry as diffraction."""
    cell = _cell()
    z = np.linspace(-30, 6, 500)
    expected = sum(cell.zDensity_G(z + n*cell.a[2], h, k) * np.exp(
        2j*np.pi*n*(h*cell.stacking_translation[0] + k*cell.stacking_translation[1])
    ) for n in range(5))
    np.testing.assert_allclose(cell.zDensity_G_asbulk(z, h, k, noUC=5), expected,
                               atol=1e-10, rtol=2e-12)


def test_dwba_vacuum_bulk_limit_and_repeat_cache():
    """The DWBA vector repeat agrees with Born scattering and invalidates caches."""
    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA native backend unavailable")
    from .. import CTRcalc
    from .test_CTRdwba import _one_atom_cell, _kinematics

    cell = _one_atom_cell(occupancy=0.0)
    cell.stacking_translation = (0.23, -0.17, 1)
    crystal = CTRcalc.SXRDCrystal(cell)
    ai, af, azimuth = 0.21, 0.34, 0.43
    h, l_value = _kinematics(cell, ai, af, azimuth)
    prepared = crystal.dwba.prepare_from_glancing(h, 0.0, ai, af)
    cell.basis[0, 6] = 1.0
    actual = crystal.dwba.evaluate_prepared(prepared).F_h
    expected = cell.F_bulk(np.array([h]), np.array([0.0]), np.array([l_value]))[0]
    np.testing.assert_allclose(actual, expected * np.cos(azimuth), rtol=2e-12)
    cell.stacking_translation = (-0.1, 0.2, 1)
    with pytest.raises(ValueError):
        crystal.dwba.evaluate_prepared(prepared)


def test_defaults_and_explicit_matrix_geometry():
    """Only Miller-based construction changes its default coordinate basis."""
    pm = pytest.importorskip("pymatgen.core")
    pytest.importorskip("pyxtal")
    source = pm.Structure(pm.Lattice.orthorhombic(5.7233, 7.9147, 5.4543),
                          ["Tb"], [[0, 0, 0]], site_properties={"B_iso": [0.6]})
    automatic = CTRsymmetry.surface_terminations(source, miller=(1, 0, 1))[0]
    assert automatic.surface_spec.bulk_transform is not None
    matrix = automatic.surface_spec.bulk_transform
    explicit = CTRsymmetry.surface_terminations(source, transform=matrix)[0]
    np.testing.assert_array_equal(explicit.surface_spec.transform, matrix)
    assert explicit.surface_spec.bulk_transform is None
    projected = CTRsymmetry.surface_terminations(
        source, transform=matrix, cell_mode="surface")[0]
    np.testing.assert_allclose(projected.surface_spec.transform,
                               automatic.surface_spec.transform)


@pytest.mark.parametrize("kind", ["coordinate", "representative"])
def test_projected_fitting_survives_supercell_and_save(kind):
    """Projected fit derivatives preserve parent-frame Cartesian movements."""
    pytest.importorskip("pyxtal")
    pytest.importorskip("pymatgen.core")
    from .test_CTRsymmetry import TestPyxtalRutileSurfaceSymmetry

    source = TestPyxtalRutileSurfaceSymmetry().make_rutile_seed()
    source.add_site_property("B_iso", [0.5] * len(source))
    cells = []
    for mode in ["lattice", "surface"]:
        option = CTRsymmetry.surface_terminations(
            source, miller=(1, 0, 1), cell_mode=mode, tol=1e-3)[0]
        cell = option.build_unitcell(mode).supercell((2, 1, 2))
        cell = CTRuc.UnitCell.fromStr(cell.toStr())
        original = cell.basis[:, 1:4].copy()
        if kind == "coordinate":
            cell.addWyckoffParameter("O_4f", "u")
            cell.setFitParameters([0.31569])
        else:
            cell.addWyckoffShift("Ru_2a", "x")
            cell.setFitParameters([0.02])
        cells.append((cell, cell.basis[:, 1:4] - original))
    native, surface = cells
    np.testing.assert_allclose(
        native[1] @ native[0].symmetry_metadata.surface_spec.transform.T,
        surface[1] @ surface[0].symmetry_metadata.surface_spec.transform.T,
        atol=2e-6,  # legacy text files retain their original coordinate precision
    )


def test_registry_propagates_through_vertical_cell_and_resets():
    """A vertical cell passes incoming registry onward and clears it on restack."""
    from ..CTRstacking import LayerState

    cell = _cell((0, 0, 1))
    state = LayerState(cell.layer_cycle, -1, (1.2, -0.7, 0))
    cell.stack_on(0, 0, below_state=state)
    np.testing.assert_allclose(cell.layer_state.registry_cartesian,
                               state.registry_cartesian)
    cell.stack_on(0, 0)
    assert cell.layer_state.registry_cartesian is None
    np.testing.assert_allclose(cell.coherentDomainMatrix[0][:, 3], 0, atol=1e-15)


def test_unsupported_export_does_not_truncate_file(tmp_path):
    """Rejecting an unsupported export leaves an existing output intact."""
    target = tmp_path / "existing.sur"
    target.write_text("existing model")
    with pytest.raises(ValueError, match="ROD export"):
        _cell().writeSURfile(target)
    assert target.read_text() == "existing model"


def test_dwba_record_refresh_and_realspace_repeats():
    """DWBA packing and real-space copies both observe a changed source repeat."""
    from ..CTRdwba import _component_record_cells

    cell = _cell()
    film = CTRfilm.Film(cell)
    film.basis[0] = 4
    _component_record_cells(film)
    cell.stacking_translation = (-0.13, 0.42, 1)
    records = _component_record_cells(film)
    expected = cell.F_uc(*_hkl()) * (1 + np.exp(
        2j*np.pi*(np.array(cell.stacking_translation) @ _hkl())))
    actual = sum(layer.F_uc(*_hkl()) for _, layer in records)
    np.testing.assert_allclose(actual, expected, atol=1e-10)
    positions = cell.pos_cart_all(ucz=2)
    xyz = np.column_stack([positions[axis] for axis in ["x", "y", "z"]])
    np.testing.assert_allclose(xyz[1::2] - xyz[::2],
                               np.tile(cell.repeat_cartesian, (2, 1)), atol=1e-12)


def test_bulk_rejects_domains_that_change_repeat():
    """A common denominator cannot silently represent different domain repeats."""
    cell = _cell()
    cell.coherentDomainMatrix[0][0, 0] = -1
    with pytest.raises(ValueError, match="preserve stacking_translation"):
        cell.F_bulk(*_hkl(), atten=0.1)


def test_film_normal_strain_preserves_lateral_registry():
    """Stretching only z leaves the lateral part of every repeat unchanged."""
    cell = _cell()
    strain = 1.07
    cell.coherentDomainMatrix[0][2, 2] = strain
    film = CTRfilm.Film(cell)
    film.basis[0] = 6
    repeat = np.array(cell.stacking_translation) * [1, 1, strain]
    expected = cell.F_uc(*_hkl()) * sum(
        np.exp(2j*np.pi*n*(repeat @ _hkl())) for n in range(3))
    np.testing.assert_allclose(film.F_uc(*_hkl()), expected, atol=1e-10)

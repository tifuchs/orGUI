"""Reference-coordinate invariants for repeated unit cells."""

import copy
import pickle

import numpy as np

from .. import CTRuc


def _cell():
    cell = CTRuc.UnitCell((3.1, 4.2, 5.3), (81, 95, 73))
    cell.addAtom("C", (0.17, 0.23, 0.31), 0.1, 0.1, 1, layer=1)
    cell.layerpos = {1: 0.0}
    return cell


def test_default_supercell_remains_self_indexed():
    """Unreferenced cells keep the returned supercell's own indexing."""
    cell = _cell()
    repeated = cell.supercell((2, 3, 1))
    np.testing.assert_array_equal(repeated.refHKLTransform, np.eye(3))
    np.testing.assert_array_equal(repeated.refRealTransform, np.eye(3))
    assert not repeated._reference_unit_cell_set


def test_explicit_identity_reference_survives_copies_and_layers():
    """Explicit identity references retain their meaning through derivation."""
    cell = _cell()
    cell.setReferenceUnitCell(cell)
    np.testing.assert_allclose(cell.refHKLTransform, np.eye(3), atol=1e-14)
    # Exercise the ambiguous case exactly, independent of roundoff in the
    # lattice inverses used by setReferenceUnitCell.
    cell.refHKLTransform[:] = np.eye(3)
    cell.refRealTransform[:] = np.eye(3)
    for source in (cell, copy.deepcopy(cell), pickle.loads(pickle.dumps(cell)),
                   *cell.split_in_layers().values()):
        repeated = source.supercell((2, 3, 1))
        np.testing.assert_allclose(repeated.B_mat @ repeated.refHKLTransform,
                                   cell.B_mat, atol=1e-13)
        np.testing.assert_allclose(repeated.R_mat @ repeated.refRealTransform,
                                   cell.R_mat, atol=1e-13)
        twice = repeated.supercell((1, 1, 2))
        np.testing.assert_allclose(twice.B_mat @ twice.refHKLTransform,
                                   cell.B_mat, atol=1e-13)


def test_rotated_reference_preserves_vectors_and_factor(monkeypatch):
    """Skew supercells preserve reference vectors and the repeated amplitude."""
    cell = _cell()
    reference = CTRuc.UnitCell((4.7, 3.8, 6.1), (85, 99, 77))
    angle = 0.37
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    cell.setReferenceUnitCell(reference, rotation)
    cell.setEnergy(15000)
    repeated = cell.supercell((2, 3, 2))
    hkl = np.array([[0.31, 0.52], [0.47, -0.14], [0.29, 0.61]])
    xyz = np.array([0.12, 0.34, -0.28])
    np.testing.assert_allclose(repeated.B_mat @ repeated.refHKLTransform @ hkl,
                               rotation @ reference.B_mat @ hkl, atol=1e-13)
    np.testing.assert_allclose(repeated.R_mat @ repeated.refRealTransform @ xyz,
                               rotation @ reference.R_mat @ xyz, atol=1e-13)
    monkeypatch.setattr(CTRuc, "CTR_ACCEL_BACKEND", "numpy")
    source_hkl = cell.refHKLTransform @ hkl
    phase_sum = sum(
        np.exp(2j * np.pi * (np.array([i, j, k]) @ source_hkl))
        for i in range(2) for j in range(3) for k in range(2)
    )
    np.testing.assert_allclose(repeated.F_uc(*hkl), cell.F_uc(*hkl) * phase_sum,
                               rtol=1e-12, atol=1e-12)


def test_nonidentity_legacy_transforms_are_preserved():
    """Preserve nonidentity reference transforms from objects without a marker."""
    cell = _cell()
    cell.refHKLTransform = np.diag([1.2, 0.8, 1.1])
    cell.refRealTransform = np.diag([0.9, 1.3, 0.7])
    del cell._reference_unit_cell_set  # Older pickles lack the marker.
    repeated = cell.supercell((2, 1, 3))
    np.testing.assert_allclose(repeated.B_mat @ repeated.refHKLTransform,
                               cell.B_mat @ cell.refHKLTransform, atol=1e-13)
    np.testing.assert_allclose(repeated.R_mat @ repeated.refRealTransform,
                               cell.R_mat @ cell.refRealTransform, atol=1e-13)

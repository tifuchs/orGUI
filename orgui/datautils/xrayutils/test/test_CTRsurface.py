"""Physical invariants for CIF import and bulk termination enumeration."""

import copy

import numpy as np
import pytest

pytest.importorskip("pymatgen")
pytest.importorskip("pyxtal")

from pymatgen.core import Lattice, Structure  # noqa: E402

from .. import CTRsymmetry as symmetry  # noqa: E402
from ..CTRuc import UnitCell  # noqa: E402


def _sto():
    return Structure(
        Lattice.cubic(4),
        [{"Sr": 0.9}, "Ti", "O", "O", "O"],
        [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]],
        site_properties={"U_iso": [0.01, 0.02, 0.03, 0.03, 0.03]},
    )


def _cif(tmp_path, text=None):
    path = tmp_path / "source.cif"
    path.write_text(
        text
        or """data_metadata
_publ_section_title 'Bibliography only'
data_crystal
_cell_length_a 4
_cell_length_b 5
_cell_length_c 6
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1 21 1'
loop_
_symmetry_equiv_pos_as_xyz
'x,y,z'
'-x,y+1/2,-z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
_atom_site_adp_type
O1 O 0.13 0.21 0.37 0.75(2) 0.015(1) Uani
loop_
_atom_site_aniso_label
_atom_site_aniso_U_11
_atom_site_aniso_U_22
_atom_site_aniso_U_33
_atom_site_aniso_U_12
_atom_site_aniso_U_13
_atom_site_aniso_U_23
O1 0.01 0.02 0.03 0.002 0.003 0.004
""",
        encoding="utf-8",
    )
    return path


def test_cif_preserves_published_scalar_and_rotates_tensor(tmp_path):
    """Published Ueq remains independent of tensor-derived Ueq on every image."""
    seed = symmetry.structure_from_cif(_cif(tmp_path))
    assert len(seed) == 2
    np.testing.assert_allclose([sum(s.species.values()) for s in seed], 0.75)
    np.testing.assert_allclose(seed.site_properties["B_iso"], 8 * np.pi**2 * 0.015)
    tensors = np.asarray(seed.site_properties["U_cart"])
    np.testing.assert_allclose(np.trace(tensors, axis1=1, axis2=2), 0.06)
    np.testing.assert_allclose(sorted(tensors[:, 0, 1]), [-0.002, 0.002], atol=1e-15)
    np.testing.assert_allclose(sorted(tensors[:, 1, 2]), [-0.004, 0.004], atol=1e-15)
    np.testing.assert_allclose(tensors[:, 0, 2], 0.003, atol=1e-15)
    assert seed.site_properties["cif_block"] == ["crystal"] * 2
    assert seed.site_properties["cif_label"] == ["O1"] * 2
    with pytest.raises(ValueError, match="Anisotropic"):
        symmetry.sites_from_seed(seed)
    sites, _, _ = symmetry.sites_from_seed(seed, displacement="isotropic_equivalent")
    assert all(np.isclose(s.iDW, 8 * np.pi**2 * 0.015) for s in sites)
    assert all(s.occ == 0.75 for s in sites)


def test_cif_b_values_and_missing_values_are_not_confused(tmp_path):
    """B input is not multiplied by 8 pi squared, and unknown B stays missing."""
    path = _cif(tmp_path)
    text = path.read_text().split("loop_\n_atom_site_aniso_label")[0]
    text = text.replace("_atom_site_U_iso_or_equiv", "_atom_site_B_iso_or_equiv")
    text = text.replace("0.015(1) Uani", "0.8(1) Biso")
    seed = symmetry.structure_from_cif(_cif(tmp_path, text))
    assert seed.site_properties["B_iso"] == [0.8, 0.8]
    seed = symmetry.structure_from_cif(_cif(tmp_path, text.replace("0.8(1)", "?")))
    assert seed.site_properties["B_iso"] == [None, None]
    with pytest.raises(ValueError, match="Missing source"):
        symmetry.sites_from_seed(seed, fallback_B=None)


def test_monoclinic_cif_tensor_uses_crystallographic_metric(tmp_path):
    """Isotropic Cartesian motion has nonzero CIF U13 when beta is 120 degrees."""
    path = _cif(tmp_path)
    text = path.read_text().replace("_cell_angle_beta 90", "_cell_angle_beta 120")
    text = text.replace("'P 1 21 1'", "'P 1'").replace("'-x,y+1/2,-z'\n", "")
    text = text.replace("0.015(1) Uani", "? Uani")
    text = text.replace(
        "O1 0.01 0.02 0.03 0.002 0.003 0.004", "O1 0.02 0.02 0.02 0 0.01 0"
    )
    seed = symmetry.structure_from_cif(_cif(tmp_path, text))
    assert seed.site_properties["B_iso"] == [None]
    np.testing.assert_allclose(
        seed[0].properties["U_cart"], 0.02 * np.eye(3), atol=1e-15
    )
    sites, _, _ = symmetry.sites_from_seed(seed, displacement="isotropic_equivalent")
    assert sites[0].iDW == pytest.approx(8 * np.pi**2 * 0.02)


def test_cif_tensor_inconsistent_with_site_symmetry_fails(tmp_path):
    """Do not select an arbitrary rotation for a non-invariant special-site ADP."""
    path = _cif(tmp_path)
    text = path.read_text().replace("'P 1 21 1'", "'P 1 2 1'")
    text = text.replace("'-x,y+1/2,-z'", "'-x,y,-z'")
    text = text.replace("0.13 0.21 0.37", "0 0.21 0")
    with pytest.raises(ValueError, match="violates site symmetry"):
        symmetry.structure_from_cif(_cif(tmp_path, text))


def test_cif_multiple_blocks_require_selection(tmp_path):
    """Bibliography blocks are ignored; multiple structures require a choice."""
    path = _cif(tmp_path)
    text = path.read_text()
    second = text[text.index("data_crystal") :].replace("data_crystal", "data_second")
    path = _cif(tmp_path, text + second)
    with pytest.raises(ValueError, match="Select a CIF block"):
        symmetry.structure_from_cif(path)
    assert len(symmetry.structure_from_cif(path, block="second")) == 2
    with pytest.raises(ValueError, match="No structural"):
        symmetry.structure_from_cif(path, block="metadata")


@pytest.mark.parametrize("value", ["?", ".", "1.2", "0", "-0.1"])
def test_cif_invalid_occupancy_is_not_normalized(tmp_path, value):
    """Unknown and unsupported occupancy values fail before parser normalization."""
    path = _cif(tmp_path)
    path = _cif(tmp_path, path.read_text().replace("0.75(2)", value))
    with pytest.raises(ValueError, match="occupancies"):
        symmetry.structure_from_cif(path)


def test_seed_preserves_occupancy_and_displacements_without_mutation():
    """Geometry assignment must not fill vacancies or mutate the caller's seed."""
    seed = _sto()
    original = seed.as_dict()
    sites, _, _ = symmetry.sites_from_seed(seed)
    values = {s.element: (s.occ, s.iDW, s.oDW) for s in sites}
    np.testing.assert_allclose(
        values["Sr"], [0.9, 8 * np.pi**2 * 0.01, 8 * np.pi**2 * 0.01]
    )
    np.testing.assert_allclose(
        values["O"], [1, 8 * np.pi**2 * 0.03, 8 * np.pi**2 * 0.03]
    )
    assert seed.as_dict() == original


def test_mixed_sites_are_rejected():
    """A mixed site must not silently become its first element."""
    seed = Structure(Lattice.cubic(4), [{"Sr": 0.8, "Ca": 0.2}], [[0, 0, 0]])
    with pytest.raises(ValueError, match="Mixed-species"):
        symmetry.sites_from_seed(seed)
    with pytest.raises(ValueError, match="Mixed-species"):
        symmetry.surface_terminations(seed, miller=(0, 0, 1))


@pytest.mark.parametrize("property_name", ["occupancy", "U_iso"])
def test_inconsistent_orbits_fail_instead_of_using_one_representative(property_name):
    """All symmetry-equivalent source atoms must agree on inherited properties."""
    seed = _sto()
    if property_name == "occupancy":
        seed.replace(2, {"O": 0.8}, properties=seed[2].properties)
    else:
        seed.add_site_property("U_iso", [0.01, 0.02, 0.04, 0.03, 0.03])
    with pytest.raises(ValueError, match="varies within Wyckoff orbit"):
        symmetry.sites_from_seed(seed)
    sites, _, _ = symmetry.sites_from_seed(seed, occ=1, iDW=0.4, oDW=0.6)
    assert all((s.occ, s.iDW, s.oDW) == (1, 0.4, 0.6) for s in sites)


def test_missing_displacement_requires_explicit_fallback_for_new_api():
    """New construction rejects missing B while legacy construction warns."""
    seed = Structure(Lattice.cubic(4), ["Si"], [[0, 0, 0]])
    option = symmetry.surface_terminations(seed, miller=(0, 0, 1))[0]
    with pytest.raises(ValueError, match="Missing source"):
        option.build_unitcell("Si")
    with pytest.warns(UserWarning, match="fallback_B=0.0"):
        cell = option.build_unitcell("Si", fallback_B=0.0)
    np.testing.assert_allclose(cell.basis[:, 4:6], 0)
    with pytest.warns(UserWarning, match="fallback_B=0.5"):
        sites, _, _ = symmetry.sites_from_seed(seed)
    assert sites[0].iDW == 0.5


def test_sto_candidates_and_built_scattering_match_manual_reference():
    """Both chemical terminations preserve the full repeat and its scattering."""
    options = symmetry.surface_terminations(_sto(), miller=(0, 0, 1))
    assert len(options) == 2
    assert options[0].layer_sequence[0]["composition"] == {"O": 1.0, "Sr": 0.9}
    assert options[1].layer_sequence[0]["composition"] == {"O": 2.0, "Ti": 1.0}
    for option in options:
        cell = option.build_unitcell("STO")
        assert len(cell.basis) == 5
        assert np.all((cell.basis[:, 3] < 0) & (cell.basis[:, 3] >= -1))
        np.testing.assert_allclose(cell.basis[:, 6].sum(), 4.9)
        reference = UnitCell(cell.a, np.rad2deg(cell.alpha), name="manual")
        for layer in option.layer_sequence:
            for atom in layer["atoms"]:
                b = 8 * np.pi**2 * atom["U_iso"]
                reference.addAtom(
                    atom["element"], atom["fractional"], b, b, atom["occupancy"]
                )
        cell.lookupScatteringFactors(10000)
        reference.lookupScatteringFactors(10000)
        np.testing.assert_allclose(
            cell.F_uc(
                np.array([0.0, 1.0, 0.5]),
                np.array([0.0, 0.0, 0.2]),
                np.array([0.7, 1.3, 2.1]),
            ),
            reference.F_uc(
                np.array([0.0, 1.0, 0.5]),
                np.array([0.0, 0.0, 0.2]),
                np.array([0.7, 1.3, 2.1]),
            ),
            atol=1e-10,
        )
        second = option.build_unitcell("independent")
        cell.basis[:, 6] = 0
        assert second.basis[:, 6].sum() == pytest.approx(4.9)


def test_repeated_cell_merges_equivalent_cuts_and_records_them():
    """Normal translations in an unreduced full matrix remain discoverable."""
    matrix = np.diag([1, 1, 3])
    options = symmetry.surface_terminations(_sto(), transform=matrix)
    raw = symmetry.surface_terminations(_sto(), transform=matrix, unique=False)
    assert len(raw) == 6
    assert len(options) == 2
    assert [len(o.equivalent_cuts) for o in options] == [3, 3]
    for option in options:
        np.testing.assert_array_equal(option.surface_spec.transform, matrix)
        assert len(option.build_unitcell("repeat").basis) == 15


def test_terminated_cell_round_trip_preserves_negative_repeat_and_metadata():
    """Serialization must retain the chosen interface position and Wyckoff data."""
    cell = symmetry.surface_terminations(_sto(), transform=np.eye(3))[1].build_unitcell(
        "STO"
    )
    restored = UnitCell.fromStr(cell.toStr())
    assert restored.names == cell.names
    np.testing.assert_allclose(restored.basis[:, 1:4], cell.basis[:, 1:4], atol=1e-8)
    # Legacy atom rows serialize B values to four decimal places and replace
    # atomic numbers with row indices; names retain the element identity.
    np.testing.assert_allclose(restored.basis[:, 4:], cell.basis[:, 4:], atol=5e-5)
    assert restored.wyckoff_sites()
    np.testing.assert_allclose(
        restored.symmetry_metadata.surface_spec.transform,
        cell.symmetry_metadata.surface_spec.transform,
    )
    rebuilt = restored.symmetry_metadata.build_unitcell("rebuilt")
    assert rebuilt.names == cell.names
    np.testing.assert_allclose(rebuilt.basis[:, 1:4], cell.basis[:, 1:4], atol=1e-8)
    np.testing.assert_allclose(rebuilt.basis[:, 4:], cell.basis[:, 4:], atol=5e-5)


def test_same_top_composition_different_stacking_is_not_merged():
    """Top-layer composition alone cannot establish termination equivalence."""
    seed = Structure(
        Lattice.orthorhombic(4, 5, 8),
        ["Sr", "O", "Ti", "O"],
        [[0, 0, 0], [0, 0, 0.2], [0, 0, 0.5], [0, 0, 0.7]],
    )
    options = symmetry.surface_terminations(seed, transform=np.eye(3))
    assert len(options) == 4
    assert sum(o.layer_sequence[0]["composition"] == {"O": 1.0} for o in options) == 2


def test_normal_reversing_operations_do_not_merge_cuts():
    """A centrosymmetric two-plane motif still has two oriented terminations."""
    seed = Structure(
        Lattice.orthorhombic(4, 5, 10), ["O", "O"], [[0.1, 0.2, 0.15], [0.1, 0.2, 0.85]]
    )
    options = symmetry.surface_terminations(seed, transform=np.eye(3))
    assert len(options) == 2


@pytest.mark.parametrize(
    "property_name, values",
    [
        ("B_iso", [0.4, 0.8]),
        (
            "U_cart",
            [
                np.diag([0.01, 0.02, 0.03]).tolist(),
                np.diag([0.02, 0.02, 0.03]).tolist(),
            ],
        ),
    ],
)
def test_displacement_properties_prevent_false_equivalence(property_name, values):
    """A geometric translation is insufficient when physical properties differ."""
    seed = Structure(
        Lattice.orthorhombic(4, 5, 8), ["O", "O"], [[0, 0, 0], [0, 0, 0.5]]
    )
    assert len(symmetry.surface_terminations(seed, transform=np.eye(3))) == 1
    seed.add_site_property(property_name, values)
    assert len(symmetry.surface_terminations(seed, transform=np.eye(3))) == 2


def test_layer_tolerance_uses_normal_distance_and_handles_periodic_buckling():
    """Oblique C length must not replace the perpendicular repeat height."""
    seed = Structure(
        Lattice.cubic(4), ["Sr", "O", "Ti"], [[0, 0, 0.99], [0, 0, 0.01], [0, 0, 0.5]]
    )
    matrix = np.array([[1, 0, 8], [0, 1, 0], [0, 0, 1]])
    options = symmetry.surface_terminations(
        seed, transform=matrix, layer_tolerance=0.1, unique=False
    )
    assert len(options) == 2
    assert (
        len(
            symmetry.surface_terminations(
                seed, transform=matrix, layer_tolerance=0.01, unique=False
            )
        )
        == 3
    )
    assert all(
        sum(layer["spacing_to_next_angstrom"] for layer in o.layer_sequence)
        == pytest.approx(4)
        for o in options
    )
    with pytest.raises(ValueError, match="no open layer gaps"):
        symmetry.surface_terminations(seed, transform=matrix, layer_tolerance=5)


def test_automatic_translation_bounds_cover_large_sheared_cell():
    """Automatic generation matches an independently expanded pymatgen cell."""
    matrix = np.array([[3, 1, 4], [0, 2, 3], [0, 0, 2]])
    seed = _sto()
    original = copy.deepcopy(seed.as_dict())
    option = symmetry.surface_terminations(seed, transform=matrix)[0]
    cell = option.build_unitcell("large")
    assert len(cell.basis) == 60
    assert cell.basis[:, 6].sum() == pytest.approx(12 * 4.9)
    expanded = seed.copy()
    expanded.make_supercell(matrix.T)
    shift = np.linalg.solve(matrix, option.surface_spec.origin)
    expected = np.mod(expanded.frac_coords - shift, 1)
    for row in cell.basis:
        eligible = np.array([next(iter(s.species)).Z == int(row[0]) for s in expanded])
        distance = expanded.lattice.get_all_distances([row[1:4]], expected[eligible])
        assert distance.min() < 1e-8
    assert seed.as_dict() == original


def test_miller_direction_and_explicit_basis_are_respected():
    """Automatic basis obeys the input setting and normal sign."""
    for hkl in ((1, 0, 1), (-1, 0, -1)):
        options = symmetry.surface_terminations(_sto(), miller=hkl, in_plane=(1, 0, -1))
        matrix = options[0].surface_spec.transform
        np.testing.assert_array_equal(matrix[:, 0], [1, 0, -1])
        np.testing.assert_allclose(np.array(hkl) @ matrix[:, :2], 0)
        assert np.dot(hkl, matrix[:, 2]) > 0
        assert np.linalg.det(matrix) > 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {},
        {"miller": (0, 0, 0)},
        {"miller": (0, 0, 0.5)},
        {"miller": (0, 0, 1), "in_plane": (0, 0, 1)},
        {"miller": (0, 0, 1), "transform": np.eye(3)},
        {"transform": np.eye(2)},
        {"transform": np.diag([1, 1, -1])},
        {"transform": np.diag([1, 1, 0.5])},
        {"transform": np.eye(3), "layer_tolerance": -1},
    ],
)
def test_invalid_surface_requests_fail_clearly(kwargs):
    """Reject ambiguous orientation and unsupported noninteger surface cells."""
    with pytest.raises(ValueError):
        symmetry.surface_terminations(_sto(), **kwargs)

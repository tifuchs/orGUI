"""Bulk-truncated termination enumeration, exported through CTRsymmetry."""

from dataclasses import dataclass, field
import math

import numpy as np

from ._CTRsource import _read_source


@dataclass(frozen=True)
class SurfaceTermination:
    """One oriented bulk termination and its equivalent cuts.

    ``surface_spec`` describes the wrapped repeat (0 <= z < 1); the exposed
    boundary is its upper edge. :meth:`build_unitcell` translates the complete
    repeat below z=0, preserving the symmetry metadata. ``cut_interval`` and
    ``equivalent_cuts`` are open gaps in the unshifted oriented fractional z
    coordinate; a wraparound gap can end above 1. ``layer_sequence`` lists
    layers from the surface inward, including occupancy-weighted composition,
    mean height and spacing in Angstrom, and individual atom dictionaries.
    Heights are measured along the true surface normal relative to the cut.
    IDs are deterministic for a given input and enumeration settings, not
    persistent crystallographic identifiers across different settings.
    """

    id: str
    surface_spec: object
    cut_interval: tuple[float, float]
    equivalent_cuts: tuple[tuple[float, float], ...]
    layer_sequence: tuple[dict, ...]
    _seed: object = field(repr=False, compare=False)
    _tol: float = field(repr=False, compare=False)
    _a_tol: float = field(repr=False, compare=False)

    @property
    def outermost_atoms(self):
        """Return atom dictionaries for the outermost grouped layer."""
        return self.layer_sequence[0]["atoms"]

    def build_unitcell(
        self,
        name,
        layer_behavior="ignore",
        *,
        displacement=None,
        fallback_B=None,
        iDW=None,
        oDW=None,
        occ=None,
    ):
        """Build a CTR UnitCell with this termination at z=0, crystal below.

        :param str name: Unit-cell name.
        :param str layer_behavior: Existing CTR UnitCell layer behavior.
        :param displacement: ``"isotropic_equivalent"`` explicitly reduces
            anisotropic source displacement data; otherwise they are rejected.
        :param fallback_B: Explicit B fallback in Angstrom squared for missing
            source data. Unlike the legacy seed helpers, no fallback is assumed.
        :param iDW: Optional in-plane B override in Angstrom squared.
        :param oDW: Optional out-of-plane B override in Angstrom squared.
        :param occ: Optional site occupancy override.
        :returns: A new UnitCell with Wyckoff metadata on each call.
        """
        from .CTRsymmetry import model_from_seed

        model = model_from_seed(
            self._seed.copy(),
            self.surface_spec,
            tol=self._tol,
            a_tol=self._a_tol,
            displacement=displacement,
            fallback_B=fallback_B,
            iDW=iDW,
            oDW=oDW,
            occ=occ,
        )
        unitcell = model.build_unitcell(name, layer_behavior=layer_behavior)
        return unitcell.translate_layered((0, 0, -len(self.layer_sequence)))


def surface_terminations(
    seed,
    *,
    miller=None,
    in_plane=None,
    transform=None,
    cell_mode=None,
    layer_tolerance=0.25,
    unique=True,
    tol=1e-5,
    a_tol=5.0,
    block=None,
):
    """Enumerate distinct bulk-truncated surfaces without adding vacuum.

    :param seed: CIF path or pymatgen Structure with preserved site properties.
    :param miller: Nonzero integer (h, k, l) in the input cell setting. Defines
        the oriented normal; reversing its sign exposes the opposite side.
    :param in_plane: Optional integer direction for the first surface vector.
        Requires ``miller``. Its common integer factor is removed.
    :param transform: Alternative full integer 3-by-3 matrix, whose columns
        are surface vectors in parent fractional coordinates. It must have
        positive determinant. It is retained exactly in lattice mode.
    :param cell_mode: ``"surface"`` projects c onto the surface normal and
        retains the crystal repeat as a stacking translation. ``"lattice"``
        retains the periodic oblique cell. Defaults to surface mode for Miller
        indices and lattice mode for an explicit ``transform`` matrix.
    :param float layer_tolerance: Maximum adjacent height separation in
        Angstrom within a layer (periodic single-linkage grouping). Cuts are
        midpoints of the remaining open gaps, never inside a grouped layer.
    :param bool unique: Merge only cuts related by property-preserving bulk
        symmetry with the same oriented surface normal. False returns all gaps.
    :param float tol: Cartesian symmetry/matching tolerance in Angstrom.
    :param float a_tol: Symmetry angle tolerance in degrees.
    :param str block: CIF data block; required for multiple structural blocks.
    :returns: Tuple of :class:`SurfaceTermination` candidates. Their labels
        are descriptive; no surface energy or preferred termination is inferred.
    :raises ValueError: For invalid/ambiguous inputs, mixed-species sites, or
        a layer tolerance joining the entire repeat into one uncuttable group.

    This enumerates bulk cuts for the specified cell and grouping, not surface
    reconstructions, vacancies added at the surface, or film/substrate registry.
    Source tensors are retained for equivalence checks even when building later
    requests an isotropic approximation. Optional dependencies remain lazy.
    """
    from .CTRsymmetry import SurfaceCellSpec

    if cell_mode is None:
        cell_mode = "surface" if transform is None else "lattice"
    if cell_mode not in {"surface", "lattice"}:
        raise ValueError("cell_mode must be 'surface' or 'lattice'.")
    if not np.isfinite(layer_tolerance) or layer_tolerance < 0:
        raise ValueError("layer_tolerance must be finite and nonnegative (Angstrom).")
    if not np.isfinite(tol) or tol <= 0 or not np.isfinite(a_tol) or a_tol < 0:
        raise ValueError("Symmetry tolerances must be finite; tol > 0 and a_tol >= 0.")
    structure = _read_source(seed, block)
    if not len(structure):
        raise ValueError("Cannot enumerate an empty structure.")
    matrix = _surface_matrix(structure, miller, in_plane, transform)
    oriented = structure.copy()
    oriented.make_supercell(matrix.T)
    cart = oriented.lattice.matrix
    normal = np.cross(cart[0], cart[1])
    normal /= np.linalg.norm(normal)
    height = float(cart[2] @ normal)
    coordinate_matrix = matrix.astype(float).copy()
    if cell_mode == "surface":
        coordinate_matrix[:, 2] = np.linalg.solve(
            structure.lattice.matrix.T, height * normal
        )
    projection = np.linalg.solve(coordinate_matrix, matrix)
    z = np.mod(oriented.frac_coords[:, 2], 1)
    order = np.argsort(z, kind="stable")
    sorted_z = z[order]
    next_z = np.r_[sorted_z[1:], sorted_z[0] + 1]
    gaps = [
        (float(lo), float(hi))
        for lo, hi in zip(sorted_z, next_z)
        if (hi - lo) * height > max(layer_tolerance, 1e-9)
    ]
    if not gaps:
        raise ValueError(
            "layer_tolerance joins the entire repeat; no open layer gaps remain."
        )
    cuts = np.array([(lo + hi) / 2 % 1 for lo, hi in gaps])
    operations = _surface_symmetries(oriented, tol, a_tol) if unique else []
    groups = []
    for index, gap in enumerate(gaps):
        for group in groups:
            representative = group[0]
            if any(
                _gaps_match(
                    gap, gaps[representative], op.translation_vector[2], tol / height
                )
                for op in operations
            ):
                group.append(index)
                break
        else:
            groups.append([index])

    candidates = []
    for group in groups:
        index = group[0]
        shift = cuts[index]
        boundaries = np.sort(np.mod(cuts - shift, 1))
        boundaries[np.isclose(boundaries, 0, atol=1e-12)] = 0
        spec = SurfaceCellSpec(
            structure.lattice.abc,
            structure.lattice.angles,
            coordinate_matrix.copy(),
            origin=tuple(matrix[:, 2] * shift),
            layer_origins=tuple(boundaries),
            translation_range=None,
            bulk_transform=matrix.copy() if cell_mode == "surface" else None,
        )
        layers = _layer_sequence(oriented, z, shift, boundaries, height)
        if cell_mode == "surface":
            for layer in layers:
                for atom in layer["atoms"]:
                    atom["fractional"] = tuple(projection @ atom["fractional"])
        candidates.append(
            SurfaceTermination(
                f"cut-{index:03d}",
                spec,
                gaps[index],
                tuple(gaps[i] for i in group),
                layers,
                structure.copy(),
                tol,
                a_tol,
            )
        )
    return tuple(candidates)


def _integer_vector(value, name):
    vector = np.asarray(value, dtype=float)
    if (
        vector.shape != (3,)
        or not np.isfinite(vector).all()
        or not np.allclose(vector, np.round(vector), atol=1e-10, rtol=0)
        or not np.any(vector)
    ):
        raise ValueError(f"{name} must be a nonzero integer three-vector.")
    vector = np.round(vector).astype(int)
    return vector // math.gcd(*map(int, vector))


def _surface_matrix(structure, miller, in_plane, transform):
    if transform is not None:
        if miller is not None or in_plane is not None:
            raise ValueError("Use either transform or miller/in_plane, not both.")
        matrix = np.asarray(transform, dtype=float)
        if (
            matrix.shape != (3, 3)
            or not np.isfinite(matrix).all()
            or not np.allclose(matrix, np.round(matrix), atol=1e-10, rtol=0)
            or np.linalg.det(matrix) < 0.5
        ):
            raise ValueError(
                "transform must be an integer 3-by-3 matrix with positive determinant."
            )
        return np.round(matrix).astype(int)
    if miller is None:
        raise ValueError("Specify either miller or transform.")
    from pymatgen.core import Structure
    from pymatgen.core.surface import SlabGenerator

    hkl = _integer_vector(miller, "miller")
    # SlabGenerator's metadata analysis is geometric; preserve the source copy.
    geometry = Structure(
        structure.lattice,
        [next(iter(site.species)).symbol for site in structure],
        structure.frac_coords,
    )
    generator = SlabGenerator(
        geometry,
        tuple(hkl),
        min_slab_size=1,
        min_vacuum_size=0,
        primitive=False,
        reorient_lattice=False,
    )
    matrix = np.array(generator.slab_scale_factor.T, dtype=int)
    if in_plane is not None:
        first = _integer_vector(in_plane, "in_plane")
        if hkl @ first != 0:
            raise ValueError("in_plane must lie in the specified Miller plane.")
        choices = [v for v in matrix[:, :2].T if np.any(np.cross(first, v))]
        second = min(
            choices, key=lambda v: np.linalg.norm(v @ structure.lattice.matrix)
        )
        matrix[:, 0], matrix[:, 1] = first, second.copy()
    if hkl @ matrix[:, 2] < 0:
        matrix[:, 2] *= -1
    if np.linalg.det(matrix) < 0:
        matrix[:, 1] *= -1
    return matrix


def _gaps_match(source, target, translation, tolerance):
    moved = np.asarray(source) + translation
    delta = moved - target
    delta -= np.round(delta[0])
    return np.all(np.abs(delta) <= tolerance)


def _surface_symmetries(structure, tol, a_tol):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from scipy.optimize import linear_sum_assignment

    operations = SpacegroupAnalyzer(
        structure, symprec=tol, angle_tolerance=a_tol
    ).get_symmetry_operations()
    a = structure.lattice.matrix.T
    inverse = np.linalg.inv(a)
    result = []
    for op in operations:
        if not np.allclose(op.rotation_matrix[2], [0, 0, 1], atol=1e-8, rtol=0):
            continue  # A surface cannot exchange its crystal and vacuum sides.
        rotation = a @ op.rotation_matrix @ inverse
        moved = op.operate_multi(structure.frac_coords)
        distances = structure.lattice.get_all_distances(moved, structure.frac_coords)
        allowed = distances <= tol
        for i, source in enumerate(structure):
            for j in np.flatnonzero(allowed[i]):
                allowed[i, j] = _properties_match(source, structure[j], rotation)
        rows, columns = linear_sum_assignment(np.where(allowed, distances, 1e20))
        if np.all(allowed[rows, columns]):
            result.append(op)
    return result


def _properties_match(source, target, rotation):
    if source.species != target.species:
        return False
    for key in ("B_iso", "U_iso", "displacement_anisotropic"):
        first, second = source.properties.get(key), target.properties.get(key)
        if (first is None) != (second is None):
            return False
        if first is not None and not np.allclose(first, second, rtol=1e-6, atol=1e-8):
            return False
    first, second = source.properties.get("U_cart"), target.properties.get("U_cart")
    if (first is None) != (second is None):
        return False
    return first is None or np.allclose(
        rotation @ np.asarray(first) @ rotation.T,
        second,
        rtol=1e-6,
        atol=1e-8,
    )


def _layer_sequence(structure, z, shift, boundaries, height):
    shifted_z = np.mod(z - shift, 1)
    layer_indices = np.searchsorted(boundaries, shifted_z, side="right") - 1
    sequence = []
    for layer in reversed(range(len(boundaries))):
        indices = np.flatnonzero(layer_indices == layer)
        atoms = []
        composition = {}
        for index in indices:
            site = structure[index]
            element = next(iter(site.species)).symbol
            occupancy = float(sum(site.species.values()))
            composition[element] = composition.get(element, 0) + occupancy
            fractional = site.frac_coords.copy()
            fractional[2] = shifted_z[index] - 1
            atoms.append(
                {
                    "label": site.properties.get("cif_label", site.label),
                    "element": element,
                    "occupancy": occupancy,
                    "fractional": tuple(map(float, fractional)),
                    "height_angstrom": float(fractional[2] * height),
                    "B_iso": site.properties.get("B_iso"),
                    "U_iso": site.properties.get("U_iso"),
                    "U_cart": site.properties.get("U_cart"),
                }
            )
        sequence.append(
            {
                "composition": dict(sorted(composition.items())),
                "height_angstrom": float(np.mean(shifted_z[indices] - 1) * height),
                "atoms": tuple(atoms),
            }
        )
    for index, layer in enumerate(sequence):
        next_height = sequence[(index + 1) % len(sequence)]["height_angstrom"]
        if index == len(sequence) - 1:
            next_height -= height
        layer["spacing_to_next_angstrom"] = layer["height_angstrom"] - next_height
    return tuple(sequence)

"""Source-property preservation for the public CTRsymmetry import helpers."""

from dataclasses import replace
from pathlib import Path
import warnings

import numpy as np


def structure_from_cif(path, block=None):
    """Read a CIF without discarding occupancy or displacement information.

    :param path: CIF path, including a :class:`pathlib.Path`.
    :param str block: Data-block name. Required if several structural blocks exist.
    :returns: A pymatgen Structure in the input cell, with source labels,
        ``B_iso`` (published B or converted Ueq, in Angstrom squared), and
        ``U_cart`` (Cartesian displacement tensors, in Angstrom squared) as
        site properties. Missing displacement data remain ``None``. Tensors
        are rotated for each symmetry image; published equivalent scalars
        are retained independently. ``cif_block`` records the selected block.
    :raises ValueError: For ambiguous blocks, labels, unsupported mixed sites,
        unknown occupancies, or invalid displacement values.
    """
    from pymatgen.io.cif import CifFile, CifParser, str2float

    document = CifFile.from_file(str(path))
    structural = {
        name: value
        for name, value in document.data.items()
        if "_atom_site_fract_x" in {key.lower() for key in value.data}
    }
    if block is None:
        if len(structural) != 1:
            raise ValueError(
                "Select a CIF block explicitly; structural blocks: "
                + ", ".join(structural)
            )
        block = next(iter(structural))
    if block not in structural:
        raise ValueError(f"No structural CIF block {block!r}.")
    parser = CifParser.from_str(
        str(structural[block]),
        frac_tolerance=0,
        site_tolerance=1e-6,
    )
    data = {key.lower(): value for key, value in structural[block].data.items()}
    labels = data["_atom_site_label"]
    if len(set(labels)) != len(labels):
        raise ValueError("CIF atom-site labels must be unique.")

    def _number(value):
        if value is None or str(value).strip() in {"?", "."}:
            return None
        result = str2float(value)
        if not np.isfinite(result):
            raise ValueError(f"Nonfinite CIF value {value!r}.")
        return result

    occupancies = [
        _number(v) for v in data.get("_atom_site_occupancy", ["1"] * len(labels))
    ]
    if any(v is None or not 0 < v <= 1 for v in occupancies):
        raise ValueError(
            "CIF occupancies must be known and in (0, 1]; no normalization is applied."
        )
    structure = parser.parse_structures(primitive=False)[0]
    _validate_species(structure)
    label_indices = {label: index for index, label in enumerate(labels)}
    b_values = [None] * len(labels)
    for key, factor in (
        ("_atom_site_b_iso_or_equiv", 1.0),
        ("_atom_site_u_iso_or_equiv", 8 * np.pi**2),
    ):
        for index, value in enumerate(data.get(key, [])):
            value = _number(value)
            if value is not None:
                value *= factor
                if value < 0:
                    raise ValueError("CIF displacement parameters must be nonnegative.")
                if b_values[index] is not None and not np.isclose(
                    b_values[index], value, rtol=1e-3, atol=1e-6
                ):
                    raise ValueError(f"Inconsistent CIF U and B for {labels[index]}.")
                b_values[index] = value

    # CIF Uij use reciprocal-axis length factors, not fractional covariance.
    a_matrix = structure.lattice.matrix.T
    to_cart = a_matrix @ np.diag(
        structure.lattice.reciprocal_lattice_crystallographic.abc
    )
    tensors = {}
    for index, label in enumerate(data.get("_atom_site_aniso_label", [])):
        if label not in label_indices or label in tensors:
            raise ValueError(f"Ambiguous anisotropic CIF label {label!r}.")
        candidates = []
        for kind, factor in (("u", 1.0), ("b", 1 / (8 * np.pi**2))):
            keys = [
                f"_atom_site_aniso_{kind}_{ij}"
                for ij in ("11", "22", "33", "12", "13", "23")
            ]
            if not any(key in data for key in keys):
                continue
            values = [
                _number(data[key][index]) if key in data else None for key in keys
            ]
            if any(value is None for value in values):
                raise ValueError(
                    f"Incomplete anisotropic displacement tensor for {label}."
                )
            u11, u22, u33, u12, u13, u23 = np.asarray(values) * factor
            tensor = (
                to_cart
                @ np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
                @ to_cart.T
            )
            if np.linalg.eigvalsh(tensor).min() < -1e-8:
                raise ValueError(f"Nonphysical displacement tensor for {label}.")
            candidates.append(tensor)
        if len(candidates) == 2 and not np.allclose(
            candidates[0], candidates[1], rtol=1e-3, atol=1e-6
        ):
            raise ValueError(f"Inconsistent anisotropic U and B for {label}.")
        if candidates:
            tensors[label] = candidates[0]

    source_xyz = np.array(
        [
            [_number(data[f"_atom_site_fract_{axis}"][index]) for axis in "xyz"]
            for index in range(len(labels))
        ],
        dtype=float,
    )
    if not np.isfinite(source_xyz).all():
        raise ValueError("CIF coordinates must be finite.")
    operations = parser.get_symops(structural[block])
    inverse_a = np.linalg.inv(a_matrix)
    properties = {
        "B_iso": [],
        "U_cart": [],
        "cif_label": [],
        "cif_block": [],
        "displacement_anisotropic": [],
    }
    adp_types = data.get(
        "_atom_site_adp_type",
        data.get("_atom_site_thermal_displace_type", [""] * len(labels)),
    )
    for site in structure:
        if site.label not in label_indices:
            raise ValueError(
                f"Cannot associate parsed site {site.label!r} with a CIF row."
            )
        index = label_indices[site.label]
        if not np.isclose(
            sum(site.species.values()), occupancies[index], atol=1e-10, rtol=0
        ):
            raise ValueError(f"CIF parsing changed occupancy at {site.label}.")
        tensor = tensors.get(site.label)
        if tensor is not None:
            images = [op.operate(source_xyz[index]) for op in operations]
            distances = structure.lattice.get_all_distances([site.frac_coords], images)[
                0
            ]
            match = int(np.argmin(distances))
            if distances[match] > 1e-5:
                raise ValueError(
                    f"Cannot match CIF tensor to {site.label} symmetry image."
                )
            rotated = []
            for image in np.flatnonzero(distances <= 1e-5):
                rotation = a_matrix @ operations[image].rotation_matrix @ inverse_a
                rotated.append(rotation @ tensor @ rotation.T)
            if not np.allclose(rotated, rotated[0], rtol=1e-5, atol=1e-7):
                raise ValueError(
                    f"Displacement tensor violates site symmetry at {site.label}."
                )
            tensor = rotated[0]
        properties["B_iso"].append(b_values[index])
        properties["U_cart"].append(None if tensor is None else tensor.tolist())
        properties["cif_label"].append(site.label)
        properties["cif_block"].append(block)
        properties["displacement_anisotropic"].append(
            tensor is not None or str(adp_types[index]).lower() in {"uani", "bani"}
        )
    for key, values in properties.items():
        structure.add_site_property(key, values)
    _validate_species(structure)
    return structure


def _validate_species(structure):
    for site in structure:
        if len(site.species) != 1:
            raise ValueError(
                "Mixed-species sites are not supported by CTR symmetry models."
            )
        occupancy = float(sum(site.species.values()))
        if not np.isfinite(occupancy) or not 0 < occupancy <= 1:
            raise ValueError("Site occupancies must be finite and in (0, 1].")
        for key in ("B_iso", "U_iso"):
            value = site.properties.get(key)
            if value is not None and (
                not np.isscalar(value) or not np.isfinite(value) or value < 0
            ):
                raise ValueError(f"{key} must be finite and nonnegative.")
        tensor = site.properties.get("U_cart")
        if tensor is not None:
            tensor = np.asarray(tensor, dtype=float)
            if (
                tensor.shape != (3, 3)
                or not np.isfinite(tensor).all()
                or not np.allclose(tensor, tensor.T, atol=1e-10, rtol=0)
                or np.linalg.eigvalsh(tensor).min() < -1e-8
            ):
                raise ValueError(
                    "U_cart must be a positive semidefinite Cartesian tensor."
                )


def _read_source(seed, block=None):
    from pymatgen.core import Structure

    if isinstance(seed, (str, Path)):
        return structure_from_cif(seed, block=block)
    if not isinstance(seed, Structure):
        raise TypeError("Expected a CIF path or pymatgen Structure.")
    if block is not None:
        raise ValueError("block is only applicable to CIF paths.")
    _validate_species(seed)
    return seed.copy()


def _site_b(site, displacement=None):
    if displacement not in (None, "isotropic_equivalent"):
        raise ValueError("displacement must be None or 'isotropic_equivalent'.")
    props = site.properties
    tensor = props.get("U_cart")
    anisotropic = props.get("displacement_anisotropic", False) or tensor is not None
    if anisotropic and displacement != "isotropic_equivalent":
        raise ValueError(
            "Anisotropic source data require displacement='isotropic_equivalent' "
            "or explicit iDW and oDW overrides."
        )
    b = props.get("B_iso")
    u = props.get("U_iso")
    if u is not None:
        converted = 8 * np.pi**2 * float(u)
        if b is not None and not np.isclose(b, converted, rtol=1e-3, atol=1e-6):
            raise ValueError("Inconsistent B_iso and U_iso site properties.")
        b = converted
    if b is None and tensor is not None:
        tensor = np.asarray(tensor, dtype=float)
        if tensor.shape != (3, 3) or not np.isfinite(tensor).all():
            raise ValueError("U_cart must be a finite 3-by-3 Cartesian tensor.")
        b = 8 * np.pi**2 * np.trace(tensor) / 3
    if b is not None and (not np.isfinite(b) or b < 0):
        raise ValueError("Displacement parameters must be finite and nonnegative.")
    return None if b is None else float(b)


def _restore_site_properties(
    sites, structure, tol, occ, iDW, oDW, displacement, fallback_B
):
    """Map every source atom, refusing orbits with inconsistent properties."""
    from scipy.optimize import linear_sum_assignment

    if displacement not in (None, "isotropic_equivalent"):
        raise ValueError("displacement must be None or 'isotropic_equivalent'.")
    _validate_species(structure)
    missing = False
    result = []
    assigned = set()
    for site in sites:
        positions = site.parent_positions()
        # Conventional orbits can contain equivalent images in primitive seeds.
        # Require a bijection after deduplication, not just geometric coverage.
        distinct = []
        for position in positions:
            if not distinct or structure.lattice.get_all_distances(
                [position], distinct
            ).min() > max(2 * tol, 1e-6):
                distinct.append(position)
        positions = distinct
        distances = structure.lattice.get_all_distances(
            structure.frac_coords, positions
        )
        matches = [
            i
            for i, atom in enumerate(structure)
            if next(iter(atom.species)).symbol == site.element
            and distances[i].min() <= max(2 * tol, 1e-6)
        ]
        if not matches or np.any(distances[matches].min(axis=0) > max(2 * tol, 1e-6)):
            raise ValueError(
                f"Cannot map all source atoms for Wyckoff orbit {site.site_id}."
            )
        if len(matches) != len(positions):
            raise ValueError(f"Source atom count differs for orbit {site.site_id}.")
        rows, columns = linear_sum_assignment(distances[matches])
        if np.any(distances[matches][rows, columns] > max(2 * tol, 1e-6)):
            raise ValueError(f"No one-to-one source mapping for orbit {site.site_id}.")
        if assigned.intersection(matches):
            raise ValueError("Ambiguous source-to-Wyckoff orbit assignment.")
        assigned.update(matches)
        values = []
        for index in matches:
            atom = structure[index]
            b = _site_b(atom, displacement) if iDW is None or oDW is None else None
            if b is None and (iDW is None or oDW is None):
                if fallback_B is None:
                    raise ValueError(
                        "Missing source displacement data; specify fallback_B "
                        "in Angstrom squared."
                    )
                b = fallback_B
                missing = True
            values.append(
                (
                    sum(atom.species.values()) if occ is None else occ,
                    b if iDW is None else iDW,
                    b if oDW is None else oDW,
                )
            )
        if not np.allclose(values, values[0], atol=1e-8, rtol=1e-6):
            raise ValueError(
                "Source occupancy/displacement varies within Wyckoff orbit "
                f"{site.site_id}; "
                "use a lower-symmetry model or explicit overrides."
            )
        occupancy, in_plane, out_of_plane = map(float, values[0])
        if (
            not np.isfinite(values).all()
            or not 0 <= occupancy <= 1
            or min(in_plane, out_of_plane) < 0
        ):
            raise ValueError("Invalid occupancy or displacement override.")
        result.append(replace(site, occ=occupancy, iDW=in_plane, oDW=out_of_plane))
    if len(assigned) != len(structure):
        raise ValueError("Some source atoms were lost during Wyckoff assignment.")
    if missing:
        warnings.warn(
            "Source displacement data missing; using "
            f"fallback_B={fallback_B} Angstrom squared.",
            UserWarning,
            stacklevel=3,
        )
    return tuple(result)

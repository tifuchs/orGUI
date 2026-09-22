"""Independent reciprocal-volume reference calculations.

These routines deliberately use only a callable coordinate transform. They do
not share determinant or subdivision code with the native reconstruction
kernel, so they can serve as a bounded validation oracle for that kernel.
They are much too slow for production detector mapping.
"""

from __future__ import annotations

import numpy as np


def _coordinate_array(coordinate, u, v, t):
    value = np.asarray(coordinate(float(u), float(v), float(t)), dtype=np.float64)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError("coordinate must return one finite three-vector")
    return value


def _centred_secant_jacobian(corners):
    """Return centred-secant columns for corners indexed by ``(u, v, t)``."""
    corners = np.asarray(corners, dtype=np.float64)
    if corners.shape != (2, 2, 2, 3):
        raise ValueError("corners must have shape (2, 2, 2, 3)")
    du = np.mean(corners[1] - corners[0], axis=(0, 1))
    dv = np.mean(corners[:, 1] - corners[:, 0], axis=(0, 1))
    dt = np.mean(corners[:, :, 1] - corners[:, :, 0], axis=(0, 1))
    return np.column_stack((du, dv, dt))


def _jacobian_diagnostics(jacobian):
    """Return signed volume, conditioning, and normalized determinant."""
    jacobian = np.asarray(jacobian, dtype=np.float64)
    signed_volume = float(np.linalg.det(jacobian))
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    smallest = float(singular_values[-1])
    condition = (
        float(singular_values[0] / smallest)
        if smallest > 0.0
        else float("inf")
    )
    column_norm_product = float(np.prod(np.linalg.norm(jacobian, axis=0)))
    normalized_determinant = (
        abs(signed_volume) / column_norm_product
        if column_norm_product > 0.0
        else 0.0
    )
    return {
        "signed_volume": signed_volume,
        "volume": abs(signed_volume),
        "condition": condition,
        "normalized_determinant": normalized_determinant,
    }


def _dyadic_leaf_oracle(coordinate, depth):
    """Evaluate root and uniform dyadic-leaf centred-secant volumes.

    The returned leaf centres use the exact coordinate callable rather than a
    trilinear average of their corners. This matches the native deepest-leaf
    voxel assignment while keeping the volume calculation independent.
    """
    depth = int(depth)
    if not 0 <= depth <= 5:
        raise ValueError("validation depth must be between 0 and 5")
    count = 1 << depth
    edges = np.linspace(0.0, 1.0, count + 1)
    lattice = np.empty((count + 1, count + 1, count + 1, 3), dtype=np.float64)
    for iu, u in enumerate(edges):
        for iv, v in enumerate(edges):
            for it, t in enumerate(edges):
                lattice[iu, iv, it] = _coordinate_array(coordinate, u, v, t)

    selected = np.asarray([0, count], dtype=np.int64)
    root_corners = lattice[np.ix_(selected, selected, selected)]
    root = _jacobian_diagnostics(_centred_secant_jacobian(root_corners))
    signed = np.empty((count, count, count), dtype=np.float64)
    volumes = np.empty_like(signed)
    centres = np.empty((*signed.shape, 3), dtype=np.float64)
    for iu in range(count):
        for iv in range(count):
            for it in range(count):
                corners = lattice[iu : iu + 2, iv : iv + 2, it : it + 2]
                value = float(np.linalg.det(_centred_secant_jacobian(corners)))
                signed[iu, iv, it] = value
                volumes[iu, iv, it] = abs(value)
                centres[iu, iv, it] = _coordinate_array(
                    coordinate,
                    0.5 * (edges[iu] + edges[iu + 1]),
                    0.5 * (edges[iv] + edges[iv + 1]),
                    0.5 * (edges[it] + edges[it + 1]),
                )

    mean_volume = float(np.mean(volumes))
    root_fraction = root["volume"] / volumes.size
    scale = max(mean_volume, np.finfo(np.float64).tiny)
    nonzero_signs = np.sign(signed[np.abs(signed) > scale * 1e-13])
    return {
        "root": root,
        "leaf_signed_volumes": signed,
        "leaf_volumes": volumes,
        "leaf_centres": centres,
        "leaf_total_volume": float(np.sum(volumes)),
        "leaf_volume_cv": (
            float(np.std(volumes) / mean_volume) if mean_volume > 0.0 else 0.0
        ),
        "leaf_root_fraction_rms": float(
            np.sqrt(np.mean(((volumes - root_fraction) / scale) ** 2))
        ),
        "leaf_root_fraction_max": float(
            np.max(np.abs(volumes - root_fraction)) / scale
        ),
        "orientation_sign_change": bool(
            np.any(nonzero_signs < 0.0) and np.any(nonzero_signs > 0.0)
        ),
    }


def _gauss_legendre_volume(coordinate, order=3, difference_step=1e-5):
    """Integrate ``abs(det(dy/d(u,v,t)))`` over the unit source cell.

    Tensor-product Gauss-Legendre nodes provide the quadrature and centred
    finite differences estimate each Jacobian column. This intentionally adds
    coordinate evaluations and is suitable only for a small validation sample.
    """
    order = int(order)
    if not 1 <= order <= 8:
        raise ValueError("order must be between 1 and 8")
    difference_step = float(difference_step)
    if not np.isfinite(difference_step) or difference_step <= 0.0:
        raise ValueError("difference_step must be finite and positive")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes = 0.5 * (nodes + 1.0)
    weights = 0.5 * weights
    boundary_distance = float(np.min(np.minimum(nodes, 1.0 - nodes)))
    step = min(difference_step, 0.25 * boundary_distance)
    volume = 0.0
    signed_volume = 0.0
    signs = []
    conditions = []
    normalized_determinants = []
    for iu, u in enumerate(nodes):
        for iv, v in enumerate(nodes):
            for it, t in enumerate(nodes):
                jacobian = np.column_stack(
                    (
                        (
                            _coordinate_array(coordinate, u + step, v, t)
                            - _coordinate_array(coordinate, u - step, v, t)
                        )
                        / (2.0 * step),
                        (
                            _coordinate_array(coordinate, u, v + step, t)
                            - _coordinate_array(coordinate, u, v - step, t)
                        )
                        / (2.0 * step),
                        (
                            _coordinate_array(coordinate, u, v, t + step)
                            - _coordinate_array(coordinate, u, v, t - step)
                        )
                        / (2.0 * step),
                    )
                )
                diagnostics = _jacobian_diagnostics(jacobian)
                quadrature_weight = weights[iu] * weights[iv] * weights[it]
                volume += quadrature_weight * diagnostics["volume"]
                signed_volume += (
                    quadrature_weight * diagnostics["signed_volume"]
                )
                signs.append(np.sign(diagnostics["signed_volume"]))
                conditions.append(diagnostics["condition"])
                normalized_determinants.append(
                    diagnostics["normalized_determinant"]
                )
    signs = np.asarray(signs)
    nonzero = signs[signs != 0.0]
    return {
        "volume": float(volume),
        "signed_volume": float(signed_volume),
        "condition_max": float(np.max(conditions)),
        "normalized_determinant_min": float(
            np.min(normalized_determinants)
        ),
        "orientation_sign_change": bool(
            np.any(nonzero < 0.0) and np.any(nonzero > 0.0)
        ),
        "order": order,
        "difference_step": step,
    }

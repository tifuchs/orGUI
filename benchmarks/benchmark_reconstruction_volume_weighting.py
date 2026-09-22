"""Validate production reciprocal-volume weighting against bounded oracles.

The native root-cell mode is compared with two deliberately slower references:

* uniform dyadic leaves calculate one centred-secant determinant per leaf and
  allocate that volume at the exact leaf centre;
* tensor-product Gauss-Legendre quadrature integrates ``abs(det(J))`` for a
  bounded subset of source cells using numerical derivatives.

The default 8 x 8 synthetic detector is small because the reference paths call
the coordinate transform from Python hundreds of times per source pixel. A
prepared job can be supplied to validate a real geometry and corrected frame.

Run with ``PYTHONPATH`` pointing at the checkout so the report names and uses
the native extension that was just rebuilt.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from benchmarks.benchmark_reconstruction_subdivision import (
    _job_case,
    _synthetic_case,
)
from benchmarks.reconstruction_volume_oracle import (
    _dyadic_leaf_oracle,
    _gauss_legendre_volume,
)
from orgui.datautils.xrayutils.reconstruction import _native_module


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job", type=Path, help="prepared reconstruction job")
    parser.add_argument("--frame", type=int, default=1800)
    parser.add_argument("--rows", type=int, default=8)
    parser.add_argument("--columns", type=int, default=8)
    parser.add_argument("--row-origin", type=int, default=1000)
    parser.add_argument("--column-origin", type=int, default=1000)
    parser.add_argument("--synthetic-sweep-degrees", type=float, default=0.1)
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--quadrature-order", type=int, default=3)
    parser.add_argument("--quadrature-pixels", type=int, default=8)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def _kernel(case, depth, threads):
    return _native_module().ReconstructionKernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        case["frame"],
        case["wavevector"],
        case["ub_inverse"],
        case["u_inverse"],
        depth,
        threads,
        max(1, 16384 >> depth),
        1 << 50,
        "reciprocal_volume_average",
    )


def _record_key(chunk, local):
    return (int(chunk) << 32) | int(local)


def _native_records(result):
    records = {}
    for index in range(len(result["chunk_id"])):
        key = _record_key(
            result["chunk_id"][index], result["local_voxel_id"][index]
        )
        records[key] = np.asarray(
            [
                result["weighted_intensity"][index],
                result["weighted_variance"][index],
                result["weight"][index],
                result["contributors"][index],
            ],
            dtype=np.float64,
        )
    return records


def _voxel_key(case, coordinate):
    index = np.floor(
        (np.asarray(coordinate) - case["minimum"]) / case["step"]
    ).astype(np.int64)
    if np.any(index < 0) or np.any(index >= case["shape"]):
        return None
    chunk_shape = case["chunk_shape"]
    chunk_index, local_index = np.divmod(index, chunk_shape)
    chunk_grid = (case["shape"] + chunk_shape - 1) // chunk_shape
    chunk = int(
        (chunk_index[0] * chunk_grid[1] + chunk_index[1])
        * chunk_grid[2]
        + chunk_index[2]
    )
    local = int(
        (local_index[0] * chunk_shape[1] + local_index[1])
        * chunk_shape[2]
        + local_index[2]
    )
    return _record_key(chunk, local)


def _coordinate(kernel, case, row, column):
    def evaluate(u, v, t):
        return kernel.coordinate(
            case["rays"],
            case["angles_start"],
            np.ascontiguousarray(case["angles_start"] + case["sweep"]),
            row,
            column,
            u,
            v,
            t,
        )

    return evaluate


def _leaf_reference(kernel, case, depth):
    records = {}
    diagnostics = []
    pixels = []
    for row, column in np.ndindex(case["intensity"].shape):
        if case["mask"][row, column]:
            continue
        coordinate = _coordinate(kernel, case, row, column)
        oracle = _dyadic_leaf_oracle(coordinate, depth)
        per_voxel = {}
        for centre, volume in zip(
            oracle["leaf_centres"].reshape(-1, 3),
            oracle["leaf_volumes"].ravel(),
        ):
            key = _voxel_key(case, centre)
            if key is not None and volume > 0.0:
                per_voxel[key] = per_voxel.get(key, 0.0) + float(volume)
        for key, weight in per_voxel.items():
            record = records.setdefault(key, np.zeros(4, dtype=np.float64))
            record[0] += weight * case["intensity"][row, column]
            record[1] += weight**2 * case["variance"][row, column]
            record[2] += weight
            record[3] += 1.0
        diagnostics.append(
            {
                "row": row,
                "column": column,
                "root_volume": oracle["root"]["volume"],
                "leaf_volume": oracle["leaf_total_volume"],
                "root_condition": oracle["root"]["condition"],
                "root_normalized_determinant": oracle["root"][
                    "normalized_determinant"
                ],
                "leaf_volume_cv": oracle["leaf_volume_cv"],
                "leaf_root_fraction_rms": oracle[
                    "leaf_root_fraction_rms"
                ],
                "leaf_root_fraction_max": oracle[
                    "leaf_root_fraction_max"
                ],
                "orientation_sign_change": oracle[
                    "orientation_sign_change"
                ],
            }
        )
        pixels.append((row, column))
    return records, diagnostics, pixels


def _weighted_quantile(values, weights, fraction):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    target = fraction * cumulative[-1]
    return float(values[order[np.searchsorted(cumulative, target)]])


def _compare_records(root, leaf):
    shared = sorted(set(root) & set(leaf))
    root_only = set(root) - set(leaf)
    leaf_only = set(leaf) - set(root)
    pulls = []
    intensity_relative = []
    weights = []
    for key in shared:
        root_record = root[key]
        leaf_record = leaf[key]
        if root_record[2] <= 0.0 or leaf_record[2] <= 0.0:
            continue
        root_intensity = root_record[0] / root_record[2]
        leaf_intensity = leaf_record[0] / leaf_record[2]
        sigma = np.sqrt(leaf_record[1]) / leaf_record[2]
        difference = abs(root_intensity - leaf_intensity)
        if sigma > 0.0:
            pulls.append(difference / sigma)
            intensity_relative.append(
                difference / max(abs(leaf_intensity), np.finfo(float).tiny)
            )
            weights.append(leaf_record[2])
    root_weight = sum(record[2] for record in root.values())
    leaf_weight = sum(record[2] for record in leaf.values())
    return {
        "shared_voxels": len(shared),
        "root_only_voxels": len(root_only),
        "leaf_only_voxels": len(leaf_only),
        "root_to_leaf_total_weight": root_weight
        / max(leaf_weight, np.finfo(float).tiny),
        "pull_weighted_median": _weighted_quantile(pulls, weights, 0.5),
        "pull_weighted_p90": _weighted_quantile(pulls, weights, 0.9),
        "pull_weighted_p99": _weighted_quantile(pulls, weights, 0.99),
        "relative_intensity_weighted_median": _weighted_quantile(
            intensity_relative, weights, 0.5
        ),
        "relative_intensity_weighted_p99": _weighted_quantile(
            intensity_relative, weights, 0.99
        ),
    }


def _summary(values):
    values = np.asarray(values, dtype=np.float64)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {
            "minimum": float("nan"),
            "median": float("nan"),
            "p99": float("nan"),
            "maximum": float("nan"),
        }
    return {
        "minimum": float(np.min(finite)),
        "median": float(np.median(finite)),
        "p99": float(np.quantile(finite, 0.99)),
        "maximum": float(np.max(finite)),
    }


def _quadrature_reference(kernel, case, pixels, depth, order, maximum):
    if maximum <= 0 or not pixels:
        return {"pixels": 0}
    positions = np.linspace(
        0, len(pixels) - 1, min(maximum, len(pixels)), dtype=int
    )
    entries = []
    for position in positions:
        row, column = pixels[position]
        coordinate = _coordinate(kernel, case, row, column)
        leaf = _dyadic_leaf_oracle(coordinate, depth)
        quadrature = _gauss_legendre_volume(coordinate, order=order)
        denominator = max(quadrature["volume"], np.finfo(float).tiny)
        entries.append(
            {
                "row": row,
                "column": column,
                "root_relative_error": (
                    leaf["root"]["volume"] - quadrature["volume"]
                )
                / denominator,
                "leaf_relative_error": (
                    leaf["leaf_total_volume"] - quadrature["volume"]
                )
                / denominator,
                "quadrature_condition_max": quadrature["condition_max"],
                "quadrature_normalized_determinant_min": quadrature[
                    "normalized_determinant_min"
                ],
                "orientation_sign_change": quadrature[
                    "orientation_sign_change"
                ],
            }
        )
    return {
        "pixels": len(entries),
        "order": order,
        "root_absolute_relative_error": _summary(
            [abs(entry["root_relative_error"]) for entry in entries]
        ),
        "leaf_absolute_relative_error": _summary(
            [abs(entry["leaf_relative_error"]) for entry in entries]
        ),
        "condition_max": _summary(
            [entry["quadrature_condition_max"] for entry in entries]
        ),
        "normalized_determinant_min": _summary(
            [entry["quadrature_normalized_determinant_min"] for entry in entries]
        ),
        "orientation_sign_change_pixels": sum(
            entry["orientation_sign_change"] for entry in entries
        ),
        "entries": entries,
    }


def main():
    """Run root/leaf/quadrature reciprocal-volume validation."""
    arguments = _arguments()
    if arguments.job:
        tile = (
            arguments.row_origin,
            arguments.row_origin + arguments.rows,
            arguments.column_origin,
            arguments.column_origin + arguments.columns,
        )
        case = _job_case(arguments.job, arguments.frame, tile)
        source = {
            "job": str(arguments.job),
            "frame": arguments.frame,
            "tile": tile,
        }
    else:
        case = _synthetic_case(
            arguments.rows,
            arguments.columns,
            arguments.row_origin,
            arguments.column_origin,
            arguments.synthetic_sweep_degrees,
        )
        source = {
            "synthetic": True,
            "rows": arguments.rows,
            "columns": arguments.columns,
            "origin": [arguments.row_origin, arguments.column_origin],
            "sweep_degrees": arguments.synthetic_sweep_degrees,
        }
    kernel = _kernel(case, arguments.depth, arguments.threads)
    angles_end = np.ascontiguousarray(case["angles_start"] + case["sweep"])

    started = perf_counter()
    native_result = kernel.accumulate(
        case["intensity"],
        case["variance"],
        case["mask"],
        case["rays"],
        case["angles_start"],
        angles_end,
        False,
    )
    native_seconds = perf_counter() - started
    root_records = _native_records(native_result)

    started = perf_counter()
    leaf_records, diagnostics, pixels = _leaf_reference(
        kernel, case, arguments.depth
    )
    leaf_seconds = perf_counter() - started
    root_to_leaf = [
        entry["root_volume"]
        / max(entry["leaf_volume"], np.finfo(float).tiny)
        for entry in diagnostics
    ]

    started = perf_counter()
    quadrature = _quadrature_reference(
        kernel,
        case,
        pixels,
        arguments.depth,
        arguments.quadrature_order,
        arguments.quadrature_pixels,
    )
    quadrature_seconds = perf_counter() - started
    report = {
        "source": source,
        "settings": {
            "depth": arguments.depth,
            "threads": arguments.threads,
            "quadrature_order": arguments.quadrature_order,
            "quadrature_pixels": arguments.quadrature_pixels,
        },
        "native": {
            "module": _native_module().__file__,
            "seconds": native_seconds,
            "records": len(root_records),
        },
        "leaf_oracle": {
            "seconds": leaf_seconds,
            "records": len(leaf_records),
            "root_to_leaf_volume": _summary(root_to_leaf),
            "leaf_volume_cv": _summary(
                [entry["leaf_volume_cv"] for entry in diagnostics]
            ),
            "leaf_root_fraction_rms": _summary(
                [entry["leaf_root_fraction_rms"] for entry in diagnostics]
            ),
            "leaf_root_fraction_max": _summary(
                [entry["leaf_root_fraction_max"] for entry in diagnostics]
            ),
            "root_condition": _summary(
                [entry["root_condition"] for entry in diagnostics]
            ),
            "root_normalized_determinant": _summary(
                [entry["root_normalized_determinant"] for entry in diagnostics]
            ),
            "orientation_sign_change_pixels": sum(
                entry["orientation_sign_change"] for entry in diagnostics
            ),
        },
        "root_vs_leaf_map": _compare_records(root_records, leaf_records),
        "quadrature": {**quadrature, "seconds": quadrature_seconds},
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if arguments.output:
        arguments.output.write_text(
            json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
        )


if __name__ == "__main__":
    main()

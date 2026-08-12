"""Depth-sweep cost and output fingerprints for the pixel-subdivision path.

Two jobs in one tool:

- **Cost.** Coordinate evaluations, leaves and nanoseconds per valid pixel at
  each ``max_depth``, for both stationary and continuous exposures. This is
  what shows the shared-corner cache switching off above depth three.
- **Fingerprints.** XXH3-128 of every output array, so a restructuring of the
  subdivision can be held to producing bit-identical records. Capture a
  baseline before changing anything, then ``--compare`` against it.

Fingerprints are only comparable within one build on one machine: floating
point contraction differs between compilers, so a mismatch against a baseline
captured elsewhere means nothing. That is why this lives in ``benchmarks``
rather than the test suite, where
``test/test_reconstruction_subdivision.py`` covers portable behaviour against
an independent reference implementation instead.

Runs against the synthetic geometry by default, or a prepared job with
``--job``.

**Run this with ``PYTHONPATH`` pointing at the checkout.** ``python
benchmarks/<script>.py`` puts ``benchmarks/`` on ``sys.path``, not the
repository root, so a pip-installed ``orgui`` shadows the working tree and the
sweep silently measures a binary you did not just build. Every report records
the resolved package path, and ``--compare`` refuses baselines that were taken
against a different one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from orgui.datautils.xrayutils.reconstruction import (
    _detector_corner_rays,
    _native_module,
    _xxh3_128,
)


RECORD_FIELDS = (
    "chunk_id",
    "local_voxel_id",
    "weighted_intensity",
    "weighted_variance",
    "weight",
    "contributors",
)


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--job", type=Path, help="prepared reconstruction job")
    parser.add_argument("--frame", type=int, default=1800)
    parser.add_argument("--rows", type=int, default=256)
    parser.add_argument("--columns", type=int, default=256)
    parser.add_argument("--row-origin", type=int, default=1000)
    parser.add_argument("--column-origin", type=int, default=1000)
    parser.add_argument("--depths", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5])
    parser.add_argument(
        "--exposures",
        choices=("stationary", "continuous"),
        nargs="+",
        default=["stationary", "continuous"],
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--compare", type=Path, help="baseline JSON to check fingerprints against"
    )
    return parser.parse_args()


def _synthetic_case(rows, columns, row_origin, column_origin):
    """A flat detector off-axis, with a grid sized to the pixel footprint.

    Mirrors the regime the real job runs in -- roughly half a voxel per pixel
    step -- because a grid much coarser than the detector resolves every cell
    immediately and measures nothing.
    """
    distance_mm = 500.0
    pitch_mm = 0.172
    wavevector = 7.601
    ub = np.array(
        [
            [1.6122, -0.5878, 0.0000],
            [0.0000, 1.4899, -0.1051],
            [0.0000, 0.0000, 1.2083],
        ]
    )
    row_edges = np.arange(rows + 1) + row_origin - 0.5
    column_edges = np.arange(columns + 1) + column_origin - 0.5
    x = np.broadcast_to(column_edges[None, :] * pitch_mm, (rows + 1, columns + 1))
    z = np.broadcast_to(row_edges[:, None] * pitch_mm, (rows + 1, columns + 1))
    rays = np.stack((x, np.full_like(x, distance_mm), z), axis=-1).astype(np.float64)
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)
    rays = np.ascontiguousarray(rays)

    angles_start = np.array([0.10471975511965978, 0.7, 0.0, 0.0])
    sweep = np.array([0.0, 0.00174533, 0.0, 0.0])

    native = _native_module()
    probe = native.ReconstructionKernel(
        np.full(3, -1e12),
        np.ones(3),
        np.full(3, 2_000_000, dtype=np.int64),
        np.ones(3, dtype=np.int64),
        "hkl",
        wavevector,
        np.ascontiguousarray(np.linalg.inv(ub)),
        np.ascontiguousarray(np.eye(3)),
        0,
        1,
        1,
        1024 * 1024,
    )

    def at(row, column, t):
        return np.asarray(
            probe.coordinate(
                rays, angles_start, angles_start + sweep, row, column, 0.5, 0.5, t
            )
        )

    centre = at(0, 0, 0.5)
    pitch = np.maximum.reduce(
        [
            np.abs(at(0, 1, 0.5) - centre),
            np.abs(at(1, 0, 0.5) - centre),
            np.abs(at(0, 0, 1.0) - at(0, 0, 0.0)),
        ]
    )
    step = np.maximum(pitch, np.max(pitch) * 0.2) / 0.5
    corners = np.asarray(
        [
            at(row, column, t)
            for row in (0, rows - 1)
            for column in (0, columns - 1)
            for t in (0.0, 1.0)
        ]
    )
    minimum = corners.min(axis=0) - 4.0 * step - 0.31718 * step
    shape = np.ceil(
        (corners.max(axis=0) + 4.0 * step - minimum) / step
    ).astype(np.int64)

    rng = np.random.default_rng(0)
    return {
        "rays": rays,
        "angles_start": np.ascontiguousarray(angles_start),
        "sweep": sweep,
        "minimum": minimum,
        "step": step,
        "shape": shape,
        "chunk_shape": np.array([64, 64, 64], dtype=np.int64),
        "frame": "hkl",
        "wavevector": wavevector,
        "ub_inverse": np.ascontiguousarray(np.linalg.inv(ub)),
        "u_inverse": np.ascontiguousarray(np.eye(3)),
        "intensity": rng.uniform(1.0, 100.0, size=(rows, columns)),
        "variance": np.abs(rng.normal(4.0, 0.5, size=(rows, columns))),
        "mask": np.zeros((rows, columns), dtype=bool),
    }


def _job_case(job_path, frame, tile):
    """The same case assembled from a prepared job and a real frame."""
    from orgui.reconstruction_job import (
        _correction_pipeline,
        _load_assets,
        read_job,
    )

    job = read_job(job_path)
    config = job.config_data
    scan = job.scan
    spec = job.internal_spec()
    grid = spec.grids[0]
    bounds = np.asarray(
        scan.exposure_angle_bounds(config, fallback=job.angle_fallback),
        dtype=np.float64,
    )
    correct = _correction_pipeline(config, scan, _load_assets(job), {})
    payload = scan.get_raw_img(frame)
    image = np.asarray(payload.img)
    corrected = correct.correct_frame(payload, image, frame)
    selection = np.s_[tile[0] : tile[1], tile[2] : tile[3]]
    ub = np.asarray(config.ub_calculator.getUB(), dtype=np.float64)
    u = np.asarray(config.ub_calculator.getU(), dtype=np.float64)
    # The job's own bounds may already be swept; a stationary/continuous pair
    # is built here from its centre so both arms are always available.
    centre = 0.5 * (bounds[frame, 0] + bounds[frame, 1])
    following = 0.5 * (bounds[frame + 1, 0] + bounds[frame + 1, 1])
    return {
        "rays": _detector_corner_rays(config.detector, tile),
        "angles_start": np.ascontiguousarray(centre),
        "sweep": np.ascontiguousarray(following - centre),
        "minimum": np.asarray(grid.minimum),
        "step": np.asarray(grid.step),
        "shape": np.asarray(grid.shape, dtype=np.int64),
        "chunk_shape": np.asarray(grid.chunk_shape, dtype=np.int64),
        "frame": grid.frame,
        "wavevector": float(config.ub_calculator.getK()),
        "ub_inverse": np.ascontiguousarray(np.linalg.inv(ub)),
        "u_inverse": np.ascontiguousarray(np.linalg.inv(u)),
        "intensity": np.ascontiguousarray(corrected[0][selection], dtype=np.float64),
        "variance": np.ascontiguousarray(corrected[1][selection], dtype=np.float64),
        "mask": np.ascontiguousarray(corrected[2][selection], dtype=bool),
    }


def _run(case, max_depth, exposure, threads, repeat):
    native = _native_module()
    kernel = native.ReconstructionKernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        case["frame"],
        case["wavevector"],
        case["ub_inverse"],
        case["u_inverse"],
        max_depth,
        threads,
        # The block-size preset halves with depth; the memory ceiling is
        # lifted because the kernel's own precheck is sized from
        # children**depth, which rejects any useful tile above depth two.
        max(1, 16384 >> max_depth),
        1 << 50,
    )
    angles_start = case["angles_start"]
    angles_end = (
        angles_start
        if exposure == "stationary"
        else np.ascontiguousarray(angles_start + case["sweep"])
    )
    best = None
    for _ in range(repeat):
        started = perf_counter()
        result = kernel.accumulate(
            case["intensity"],
            case["variance"],
            case["mask"],
            case["rays"],
            angles_start,
            angles_end,
            True,
        )
        elapsed = perf_counter() - started
        if best is None or elapsed < best[0]:
            best = (elapsed, result)
    elapsed, result = best
    profile = dict(result["_profile"])
    valid = max(1, profile["valid_pixels"])
    return {
        "depth": max_depth,
        "exposure": exposure,
        "valid_pixels": profile["valid_pixels"],
        "records": int(profile["final_records"]),
        "evaluations_per_pixel": profile["coordinate_evaluations"] / valid,
        "leaves_per_pixel": profile["voxel_weights"] / valid,
        "max_leaves_one_pixel": int(profile["maximum_weights_per_pixel"]),
        "records_per_pixel": profile["final_records"] / valid,
        "map_ns_per_pixel": profile["block_mapping_cpu_seconds"] * 1e9 / valid,
        "wall_seconds": elapsed,
        "fingerprints": {
            name: _xxh3_128(result[name]) for name in RECORD_FIELDS
        },
    }


def _native_build_identity():
    """Which extension module this process actually imported."""
    return {
        "module": _native_module().__file__,
        "mtime": Path(_native_module().__file__).stat().st_mtime,
    }


def _compare(baseline, current):
    """Report fingerprint agreement and the speed-up, case by case."""
    was_module = baseline.get("native", {}).get("module")
    now_module = current["native"]["module"]
    if was_module != now_module:
        raise SystemExit(
            "Baseline was taken against a different extension module, so "
            "nothing here would be comparable:\n"
            f"  baseline: {was_module}\n"
            f"  current:  {now_module}\n"
            "Re-capture the baseline with the same PYTHONPATH."
        )
    index = {
        (entry["depth"], entry["exposure"]): entry for entry in baseline["cases"]
    }
    identical = True
    print(
        f"{'exposure':>11} {'depth':>5} {'records':>9} {'bitwise':>8} "
        f"{'evals/px':>9} {'was':>9} {'ns/px':>10} {'was':>10} {'speedup':>8}"
    )
    for entry in current["cases"]:
        was = index.get((entry["depth"], entry["exposure"]))
        if was is None:
            print(f"{entry['exposure']:>11} {entry['depth']:>5}   (no baseline)")
            continue
        same = entry["fingerprints"] == was["fingerprints"]
        identical = identical and same
        speedup = (
            was["map_ns_per_pixel"] / entry["map_ns_per_pixel"]
            if entry["map_ns_per_pixel"]
            else float("nan")
        )
        print(
            f"{entry['exposure']:>11} {entry['depth']:>5} {entry['records']:>9} "
            f"{'same' if same else 'DIFFER':>8} "
            f"{entry['evaluations_per_pixel']:>9.2f} "
            f"{was['evaluations_per_pixel']:>9.2f} "
            f"{entry['map_ns_per_pixel']:>10.1f} "
            f"{was['map_ns_per_pixel']:>10.1f} {speedup:>8.2f}x"
        )
    print(
        "\nAll fingerprints identical."
        if identical
        else "\nFINGERPRINTS DIFFER -- output is not bit-for-bit."
    )
    return identical


def main():
    """Sweep depth and exposure, reporting cost and output fingerprints."""
    arguments = _arguments()
    if arguments.job:
        tile = (
            arguments.row_origin,
            arguments.row_origin + arguments.rows,
            arguments.column_origin,
            arguments.column_origin + arguments.columns,
        )
        case = _job_case(arguments.job, arguments.frame, tile)
        source = {"job": str(arguments.job), "frame": arguments.frame, "tile": tile}
    else:
        case = _synthetic_case(
            arguments.rows,
            arguments.columns,
            arguments.row_origin,
            arguments.column_origin,
        )
        source = {
            "synthetic": True,
            "rows": arguments.rows,
            "columns": arguments.columns,
        }

    cases = []
    for exposure in arguments.exposures:
        for max_depth in arguments.depths:
            entry = _run(
                case, max_depth, exposure, arguments.threads, arguments.repeat
            )
            cases.append(entry)
            reported = {
                key: value
                for key, value in entry.items()
                if key != "fingerprints"
            }
            print(json.dumps(reported, sort_keys=True), flush=True)
    report = {
        "source": source,
        "native": _native_build_identity(),
        "threads": arguments.threads,
        "cases": cases,
    }

    if arguments.compare:
        baseline = json.loads(arguments.compare.read_text(encoding="utf-8"))
        print()
        _compare(baseline, report)
    if arguments.output:
        arguments.output.write_text(
            json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
        )


if __name__ == "__main__":
    main()

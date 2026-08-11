"""Does raising the subdivision depth change a reconstructed intensity?

``doc/design/reciprocal_space_mapping_locality_and_subdivision.md``
observes that emitted records converge at depth 3 (0.632 -> 0.650 across
depths 3-6) while cost rises 3x/6x/13x, and asks whether ``very_high`` and
``maximum`` are worth optimising at all. That is a question about
intensities, not about record counts: a voxel can gain a finer
partial-volume weight without its intensity moving anywhere a measurement
could see.

So this compares the quantity a user actually reads --
``weighted_intensity / weight``, exactly as ``_reduce_chunk`` forms it --
across depths, against that same voxel's own propagated error bar
``sqrt(weighted_variance) / weight``. The statistic is the pull

    (I_d - I_ref) / sigma_ref

which is dimensionless and already scaled by what the measurement can
resolve. Because sigma falls as 1/sqrt(contributors) while a systematic
from the subdivision does not fall at all, the pull measured here also
says at what exposure the difference *would* become visible:

    contributors_critical = contributors_now / pull**2

That last number is the honest answer to "measurably", independent of how
many frames this particular probe could afford to run. Calibrate it
against the exposure a real job reaches: this probe's 64 frames give
20-45 contributors per voxel, and a 3651-frame job gives roughly
1,000-2,500.

Unlike the other reconstruction benchmarks this one compares numbers
rather than times, so it is immune to the machine noise documented in
``reciprocal_space_performance_findings.md`` and needs no interleaving or
repeats.

**Run this with ``PYTHONPATH`` pointing at the checkout.** ``python
benchmarks/<script>.py`` puts ``benchmarks/`` on ``sys.path``, not the
repository root, so a pip-installed ``orgui`` shadows the working tree.
The resolved extension path is reported for that reason.
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
)


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--frame", type=int, default=1800)
    parser.add_argument("--frames", type=int, default=64)
    parser.add_argument("--rows", type=int, default=192)
    parser.add_argument("--columns", type=int, default=192)
    parser.add_argument("--row-origin", type=int, default=1100)
    parser.add_argument("--column-origin", type=int, default=1100)
    parser.add_argument("--depths", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5])
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def _case(job_path, frames, tile):
    """Corrected tile arrays and exposure bounds for a run of frames."""
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
    selection = np.s_[tile[0] : tile[1], tile[2] : tile[3]]
    intensity = []
    variance = []
    mask = []
    for frame in frames:
        payload = scan.get_raw_img(frame)
        corrected = correct.correct_frame(payload, np.asarray(payload.img), frame)
        intensity.append(np.ascontiguousarray(corrected[0][selection], np.float64))
        variance.append(np.ascontiguousarray(corrected[1][selection], np.float64))
        mask.append(np.ascontiguousarray(corrected[2][selection], bool))
    ub = np.asarray(config.ub_calculator.getUB(), dtype=np.float64)
    u = np.asarray(config.ub_calculator.getU(), dtype=np.float64)
    return {
        "rays": _detector_corner_rays(config.detector, tile),
        "bounds": bounds,
        "frames": list(frames),
        "intensity": intensity,
        "variance": variance,
        "mask": mask,
        "minimum": np.asarray(grid.minimum),
        "step": np.asarray(grid.step),
        "shape": np.asarray(grid.shape, dtype=np.int64),
        "chunk_shape": np.asarray(grid.chunk_shape, dtype=np.int64),
        "frame_name": grid.frame,
        "wavevector": float(config.ub_calculator.getK()),
        "ub_inverse": np.ascontiguousarray(np.linalg.inv(ub)),
        "u_inverse": np.ascontiguousarray(np.linalg.inv(u)),
        "angle_fallback": job.angle_fallback,
    }


def _accumulate(case, max_depth, threads):
    """Merge every frame's records at one depth into per-voxel sums.

    Frames are accumulated one at a time and summed in Python rather than
    through ``accumulate_group``: the question here is what depth does,
    and mixing in the grouped kernel's different summation order would
    put a second variable in the comparison.
    """
    native = _native_module()
    kernel = native.ReconstructionKernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        case["frame_name"],
        case["wavevector"],
        case["ub_inverse"],
        case["u_inverse"],
        max_depth,
        threads,
        max(1, 16384 >> max_depth),
        1 << 50,
    )
    keys = []
    columns = {name: [] for name in
               ("weighted_intensity", "weighted_variance", "weight", "contributors")}
    started = perf_counter()
    for position, frame in enumerate(case["frames"]):
        result = kernel.accumulate(
            case["intensity"][position],
            case["variance"][position],
            case["mask"][position],
            case["rays"],
            np.ascontiguousarray(case["bounds"][frame, 0]),
            np.ascontiguousarray(case["bounds"][frame, 1]),
            False,
        )
        keys.append(
            (np.asarray(result["chunk_id"]).astype(np.int64) << np.int64(32))
            | np.asarray(result["local_voxel_id"]).astype(np.int64)
        )
        for name in columns:
            columns[name].append(np.asarray(result[name]).astype(np.float64))
    elapsed = perf_counter() - started
    key = np.concatenate(keys)
    unique, inverse = np.unique(key, return_inverse=True)
    totals = {
        name: np.bincount(
            inverse, weights=np.concatenate(values), minlength=unique.size
        )
        for name, values in columns.items()
    }
    return {
        "key": unique,
        "seconds": elapsed,
        **totals,
    }


def _compare(current, reference):
    """Pull of one depth's intensities against a deeper reference."""
    shared, here, there = np.intersect1d(
        current["key"], reference["key"], return_indices=True
    )
    weight_here = current["weight"][here]
    weight_there = reference["weight"][there]
    live = (weight_here > 0) & (weight_there > 0)
    intensity_here = current["weighted_intensity"][here][live] / weight_here[live]
    intensity_there = (
        reference["weighted_intensity"][there][live] / weight_there[live]
    )
    sigma = (
        np.sqrt(reference["weighted_variance"][there][live]) / weight_there[live]
    )
    usable = sigma > 0
    pull = np.abs(intensity_here[usable] - intensity_there[usable]) / sigma[usable]
    weights = weight_there[live][usable]
    contributors = reference["contributors"][there][live][usable]
    order = np.argsort(pull)
    cumulative = np.cumsum(weights[order])
    total = cumulative[-1] if cumulative.size else 1.0

    def weighted_quantile(fraction):
        if not cumulative.size:
            return float("nan")
        return float(pull[order][np.searchsorted(cumulative, fraction * total)])

    median_pull = weighted_quantile(0.5)
    mean_contributors = float(np.average(contributors, weights=weights))
    return {
        "shared_voxels": int(shared.size),
        "voxels_here_only": int(current["key"].size - shared.size),
        "voxels_reference_only": int(reference["key"].size - shared.size),
        "weight_shared_fraction": float(
            weight_here[live].sum() / max(1e-30, current["weight"].sum())
        ),
        "mean_contributors": mean_contributors,
        "pull_median": median_pull,
        "pull_p90": weighted_quantile(0.9),
        "pull_p99": weighted_quantile(0.99),
        "pull_rms": float(
            np.sqrt(np.average(pull**2, weights=weights))
        ),
        # Contributors at which this depth's systematic would equal the
        # statistical error bar. Beyond the probe's own exposure, so it is
        # an extrapolation -- but the scaling it rests on (sigma ~
        # 1/sqrt(N), systematic flat) is not in doubt.
        "contributors_for_unit_pull": (
            float(mean_contributors / median_pull**2)
            if median_pull > 0
            else float("inf")
        ),
        "total_weight_ratio": float(
            current["weight"].sum() / max(1e-30, reference["weight"].sum())
        ),
        "total_intensity_ratio": float(
            current["weighted_intensity"].sum()
            / max(1e-30, reference["weighted_intensity"].sum())
        ),
    }


def main():
    """Compare reconstructed intensities across subdivision depths."""
    arguments = _arguments()
    tile = (
        arguments.row_origin,
        arguments.row_origin + arguments.rows,
        arguments.column_origin,
        arguments.column_origin + arguments.columns,
    )
    frames = list(range(arguments.frame, arguments.frame + arguments.frames))
    case = _case(arguments.job, frames, tile)
    print(
        json.dumps(
            {
                "job": str(arguments.job),
                "tile": tile,
                "frames": [frames[0], frames[-1] + 1],
                "angle_fallback": case["angle_fallback"],
                "native": _native_module().__file__,
                "depths": arguments.depths,
            },
            sort_keys=True,
        ),
        flush=True,
    )
    depths = sorted(arguments.depths)
    accumulated = {}
    for depth in depths:
        accumulated[depth] = _accumulate(case, depth, arguments.threads)
        print(
            f"depth {depth}: {accumulated[depth]['key'].size:,} voxels, "
            f"{accumulated[depth]['seconds']:.1f} s",
            flush=True,
        )
    reference = depths[-1]
    report = []
    for depth in depths[:-1]:
        entry = {"depth": depth, "reference": reference}
        entry.update(_compare(accumulated[depth], accumulated[reference]))
        report.append(entry)
        print(json.dumps(entry, sort_keys=True), flush=True)
    if arguments.output:
        arguments.output.write_text(
            json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
        )


if __name__ == "__main__":
    main()

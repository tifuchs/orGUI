"""Frame-group bricks against per-image blocks, on a prepared job.

``accumulate_group`` maps several adjacent images in one call, with the work
block a brick in (row, column, frame) rather than a run of one flattened
image. This measures what that is worth at each group size: records emitted
before the cross-block merge -- which is what every stage downstream of the
kernel scales with -- and wall time.

The per-image arm is the same frames through ``accumulate``, one at a time,
merged with ``_reduce_batches``, i.e. exactly what the pipeline does today.

**Run this with ``PYTHONPATH`` pointing at the checkout.** ``python
benchmarks/<script>.py`` puts ``benchmarks/`` on ``sys.path``, not the
repository root, so a pip-installed ``orgui`` shadows the working tree and
the sweep silently measures a binary you did not just build. The resolved
package path is reported for that reason.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from orgui.datautils.xrayutils.reconstruction import (
    _detector_corner_rays,
    _kernel_for_grid,
    _native_module,
    _reduce_batches,
)
from orgui.reconstruction_job import (
    _correction_pipeline,
    _load_assets,
    read_job,
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
    parser.add_argument("job", type=Path)
    parser.add_argument("--frame", type=int, default=1800)
    parser.add_argument("--groups", type=int, nargs="+", default=[1, 2, 4, 8, 16])
    parser.add_argument("--rows", type=int, default=512)
    parser.add_argument("--row-origin", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--repeat", type=int, default=2)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def _load_window(job, frames, tile):
    """Corrected intensity, variance and mask for one window of frames."""
    config = job.config_data
    scan = job.scan
    correct = _correction_pipeline(config, scan, _load_assets(job), {})
    selection = np.s_[tile[0] : tile[1], tile[2] : tile[3]]
    rows = tile[1] - tile[0]
    columns = tile[3] - tile[2]
    intensity = np.empty((len(frames), rows, columns))
    variance = np.empty_like(intensity)
    mask = np.empty(intensity.shape, dtype=bool)
    for slot, frame in enumerate(frames):
        payload = scan.get_raw_img(frame)
        corrected = correct.correct_frame(
            payload, np.asarray(payload.img), frame
        )
        intensity[slot] = corrected[0][selection]
        variance[slot] = corrected[1][selection]
        mask[slot] = corrected[2][selection]
    return intensity, variance, mask


def _best_of(repeat, call):
    best = None
    for _ in range(repeat):
        started = perf_counter()
        result = call()
        elapsed = perf_counter() - started
        if best is None or elapsed < best[0]:
            best = (elapsed, result)
    return best


def _per_frame_arm(kernel, rays, start, end, payload, repeat):
    """Today's path: every frame mapped alone, then merged."""
    intensity, variance, mask = payload
    seconds = 0.0
    block_records = 0
    batches = []
    for slot in range(intensity.shape[0]):
        elapsed, batch = _best_of(
            repeat,
            lambda slot=slot: kernel.accumulate(
                np.ascontiguousarray(intensity[slot]),
                np.ascontiguousarray(variance[slot]),
                np.ascontiguousarray(mask[slot]),
                rays,
                np.ascontiguousarray(start[slot]),
                np.ascontiguousarray(end[slot]),
                True,
            ),
        )
        seconds += elapsed
        block_records += int(batch["_profile"]["reduced_block_records"])
        batches.append({name: batch[name] for name in RECORD_FIELDS})
    merged = _reduce_batches(batches)
    return {
        "mode": "per-frame",
        "frames": intensity.shape[0],
        "block_records": block_records,
        "records_after_merge": int(merged["chunk_id"].size),
        "seconds": seconds,
    }


def _group_arm(kernel, rays, start, end, payload, group, repeat):
    """Frames mapped in windows of ``group``, covering as many as fit."""
    intensity, variance, mask = payload
    frames = intensity.shape[0]
    windows = max(1, frames // group)
    seconds = 0.0
    block_records = 0
    records = 0
    brick = None
    for origin in range(0, windows * group, group):
        window = slice(origin, origin + group)
        elapsed, batch = _best_of(
            repeat,
            lambda window=window: kernel.accumulate_group(
                np.ascontiguousarray(intensity[window]),
                np.ascontiguousarray(variance[window]),
                np.ascontiguousarray(mask[window]),
                rays,
                np.ascontiguousarray(start[window]),
                np.ascontiguousarray(end[window]),
                True,
            ),
        )
        seconds += elapsed
        block_records += int(batch["_profile"]["reduced_block_records"])
        records += int(batch["chunk_id"].size)
        brick = [
            int(batch["_profile"]["brick_rows"]),
            int(batch["_profile"]["brick_columns"]),
            int(batch["_profile"]["frames"]),
        ]
    return {
        "mode": "group",
        "group": group,
        "frames": windows * group,
        "brick": brick,
        "block_records": block_records,
        "records_after_merge": records,
        "seconds": seconds,
    }


def main():
    """Compare frame-group bricks against per-image blocks."""
    arguments = _arguments()
    job = read_job(arguments.job)
    config = job.config_data
    spec = job.internal_spec()
    grid = spec.grids[0]
    bounds = np.asarray(
        job.scan.exposure_angle_bounds(config, fallback=job.angle_fallback),
        dtype=np.float64,
    )
    _rows, columns = config.detector.detector.shape
    tile = (
        arguments.row_origin,
        arguments.row_origin + arguments.rows,
        0,
        columns,
    )
    biggest = max(arguments.groups)
    frames = list(range(arguments.frame, arguments.frame + biggest))
    payload = _load_window(job, frames, tile)
    rays = _detector_corner_rays(config.detector, tile)
    start = np.ascontiguousarray(bounds[frames, 0])
    end = np.ascontiguousarray(bounds[frames, 1])
    kernel = _kernel_for_grid(
        spec, grid, config.ub_calculator, threads=arguments.threads
    )
    valid_samples = int((~payload[2]).sum())

    header = {
        "job": str(arguments.job),
        "tile": list(tile),
        "frames": frames[0],
        "valid_samples": valid_samples,
        "max_depth": spec.max_depth,
        "work_block_pixels": spec.work_block_pixels,
        "threads": arguments.threads,
        "native": _native_module().__file__,
    }
    print(json.dumps(header, sort_keys=True), flush=True)

    baseline = _per_frame_arm(
        kernel, rays, start, end, payload, arguments.repeat
    )
    print(json.dumps(baseline, sort_keys=True), flush=True)

    results = [baseline]
    for group in arguments.groups:
        if group > biggest:
            continue
        entry = _group_arm(
            kernel, rays, start, end, payload, group, arguments.repeat
        )
        share = entry["frames"] / baseline["frames"]
        entry["relative_block_records"] = entry["block_records"] / (
            baseline["block_records"] * share
        )
        entry["relative_seconds"] = entry["seconds"] / (
            baseline["seconds"] * share
        )
        results.append(entry)
        print(json.dumps(entry, sort_keys=True), flush=True)

    print()
    print(
        f"{'group':>6} {'brick':>16} {'block records':>14} "
        f"{'vs per-frame':>13} {'time':>8}"
    )
    for entry in results[1:]:
        brick = "x".join(str(value) for value in entry["brick"])
        print(
            f"{entry['group']:>6} {brick:>16} {entry['block_records']:>14,} "
            f"{entry['relative_block_records']:>13.3f} "
            f"{entry['relative_seconds']:>8.3f}"
        )

    if arguments.output:
        arguments.output.write_text(
            json.dumps({"header": header, "results": results}, indent=2),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()

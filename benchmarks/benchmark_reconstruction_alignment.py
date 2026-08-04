"""Measure native reconstruction sensitivity to input-buffer alignment."""

from __future__ import annotations

import argparse
from dataclasses import replace
import gc
import json
from pathlib import Path
from statistics import median
from time import perf_counter

import numpy as np

from orgui.datautils.xrayutils.reconstruction import (
    _detector_corner_rays,
    _kernel_for_grid,
)
from orgui.reconstruction_job import (
    _correction_pipeline,
    _load_assets,
    read_job,
)


def _aligned_copy(values, address_mod_64):
    values = np.ascontiguousarray(values)
    storage = np.empty(values.nbytes + 63, dtype=np.uint8)
    offset = (address_mod_64 - storage.ctypes.data) % 64
    result = storage[offset : offset + values.nbytes].view(values.dtype)
    result = result.reshape(values.shape)
    np.copyto(result, values)
    if result.ctypes.data % 64 != address_mod_64:
        raise RuntimeError("Failed to construct requested input alignment")
    return result


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument("--tile-size", type=int, default=1024)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--offsets",
        type=int,
        nargs="+",
        default=tuple(range(0, 64, 8)),
    )
    return parser.parse_args()


def main():
    """Benchmark naturally aligned and deliberately offset native inputs."""
    arguments = _arguments()
    job = read_job(arguments.job)
    config = job.config_data
    scan = job.scan
    payload = scan.get_raw_img(arguments.frame)
    image = np.asarray(payload.img)
    tile = (
        0,
        min(arguments.tile_size, image.shape[0]),
        0,
        min(arguments.tile_size, image.shape[1]),
    )
    selection = np.s_[tile[0] : tile[1], tile[2] : tile[3]]
    correction = _correction_pipeline(
        config,
        scan,
        _load_assets(job),
        {},
    )
    corrected = correction.correct_frame(payload, image, arguments.frame)
    inputs = [
        np.ascontiguousarray(corrected[0][selection], dtype=np.float64),
        np.ascontiguousarray(corrected[1][selection], dtype=np.float64),
        np.ascontiguousarray(corrected[2][selection], dtype=bool),
        _detector_corner_rays(config.detector, tile),
    ]
    bounds = scan.exposure_angle_bounds(
        config,
        fallback=job.angle_fallback,
    )[arguments.frame]
    angles_start = np.ascontiguousarray(bounds[0], dtype=np.float64)
    angles_end = np.ascontiguousarray(bounds[1], dtype=np.float64)
    spec = replace(
        job.internal_spec(),
        max_depth=arguments.depth,
    )
    kernel = _kernel_for_grid(
        spec,
        spec.grids[0],
        config.ub_calculator,
        threads=arguments.threads,
        memory_budget_bytes=spec.memory_budget_bytes,
    )

    offsets = tuple(arguments.offsets)
    if any(offset < 0 or offset >= 64 for offset in offsets):
        raise ValueError("Alignment offsets must be in [0, 64)")
    aligned_inputs = {
        offset: tuple(_aligned_copy(values, offset) for values in inputs)
        for offset in offsets
    }
    kernel.accumulate(
        *aligned_inputs[0],
        angles_start,
        angles_end,
    )
    timings = {offset: [] for offset in offsets}
    gc.disable()
    try:
        for repeat in range(arguments.repeats):
            order = offsets if repeat % 2 == 0 else tuple(reversed(offsets))
            for offset in order:
                started = perf_counter()
                kernel.accumulate(
                    *aligned_inputs[offset],
                    angles_start,
                    angles_end,
                )
                timings[offset].append(perf_counter() - started)
    finally:
        gc.enable()

    baseline = median(timings[0])
    result = {
        "tile": tile,
        "threads": arguments.threads,
        "depth": arguments.depth,
        "repeats": arguments.repeats,
        "results": [
            {
                "address_mod_64": offset,
                "median_seconds": median(timings[offset]),
                "relative_to_aligned": median(timings[offset]) / baseline,
                "samples": timings[offset],
            }
            for offset in offsets
        ],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

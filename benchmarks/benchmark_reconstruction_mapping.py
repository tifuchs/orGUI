"""Benchmark reciprocal-space mapping layouts on a prepared orGUI job."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import gc
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np

from orgui.datautils.xrayutils.reconstruction import (
    _detector_corner_rays,
    _kernel_for_grid,
    _map_frame_range,
    _reduce_batches,
    _xxh3_128,
)
import orgui.datautils.xrayutils.reconstruction as reconstruction
from orgui.reconstruction_job import (
    _correction_pipeline,
    _load_assets,
    read_job,
)


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--frame", type=int, default=0)
    parser.add_argument("--threads", type=int, nargs="+", default=[4, 6])
    parser.add_argument(
        "--layouts",
        choices=("current", "equal"),
        nargs="+",
        default=["current", "equal"],
    )
    parser.add_argument("--tile-sizes", type=int, nargs="+")
    parser.add_argument("--work-blocks", type=int, nargs="+")
    parser.add_argument("--depths", type=int, nargs="+")
    parser.add_argument("--parallel-frames", type=int, default=0)
    parser.add_argument("--native-profile", action="store_true")
    parser.add_argument("--profile-load-frames", type=int, default=0)
    parser.add_argument("--profile-load-workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def _factor_grid(count, rows, columns):
    candidates = []
    for row_count in range(1, count + 1):
        if count % row_count:
            continue
        column_count = count // row_count
        tile_rows = rows / row_count
        tile_columns = columns / column_count
        distortion = abs(math.log(tile_rows / tile_columns))
        candidates.append((distortion, row_count, column_count))
    _, row_count, column_count = min(candidates)
    return row_count, column_count


def _equal_tiles(rows, columns, count):
    row_count, column_count = _factor_grid(count, rows, columns)
    row_edges = [(index * rows) // row_count for index in range(row_count + 1)]
    column_edges = [
        (index * columns) // column_count
        for index in range(column_count + 1)
    ]
    return tuple(
        (
            row_edges[row],
            row_edges[row + 1],
            column_edges[column],
            column_edges[column + 1],
        )
        for row in range(row_count)
        for column in range(column_count)
    )


def _fixed_tiles(rows, columns, tile_rows, tile_columns):
    return tuple(
        (
            row,
            min(row + tile_rows, rows),
            column,
            min(column + tile_columns, columns),
        )
        for row in range(0, rows, tile_rows)
        for column in range(0, columns, tile_columns)
    )


def _run_case(
    *,
    name,
    tiles,
    threads,
    spec,
    job,
    config,
    payload,
    image,
    angles,
    assets,
    native_profile=False,
):
    provenance = {}
    correct = _correction_pipeline(config, job.scan, assets, provenance)
    rays = {
        tile: _detector_corner_rays(config.detector, tile)
        for tile in tiles
    }
    kernel = _kernel_for_grid(
        spec,
        spec.grids[0],
        config.ub_calculator,
        threads=threads,
        memory_budget_bytes=spec.memory_budget_bytes,
    )
    batches = []
    native_profiles = []
    correction_seconds = 0.0
    native_seconds = 0.0
    started = perf_counter()
    correction_started = perf_counter()
    corrected_frame = correct.correct_frame(payload, image, 0)
    correction_seconds = perf_counter() - correction_started
    for tile in tiles:
        row_start, row_stop, column_start, column_stop = tile
        selection = np.s_[
            row_start:row_stop, column_start:column_stop
        ]
        intensity, variance, mask = (
            values[selection] for values in corrected_frame
        )
        native_started = perf_counter()
        batch = kernel.accumulate(
                np.ascontiguousarray(intensity, dtype=np.float64),
                np.ascontiguousarray(variance, dtype=np.float64),
                np.ascontiguousarray(mask, dtype=bool),
                rays[tile],
                angles[0],
                angles[1],
                profile=native_profile,
            )
        if native_profile:
            native_profiles.append(batch.pop("_profile"))
        batches.append(batch)
        native_seconds += perf_counter() - native_started
    reduction_started = perf_counter()
    reduced = _reduce_batches(batches)
    reduction_seconds = perf_counter() - reduction_started
    elapsed = perf_counter() - started
    rows, columns = image.shape
    result = {
        "name": name,
        "threads": threads,
        "tiles": len(tiles),
        "largest_tile_pixels": max(
            (tile[1] - tile[0]) * (tile[3] - tile[2])
            for tile in tiles
        ),
        "seconds": elapsed,
        "megapixels_per_second": rows * columns / elapsed / 1e6,
        "correction_seconds": correction_seconds,
        "native_seconds": native_seconds,
        "reduction_seconds": reduction_seconds,
        "records": int(reduced["chunk_id"].size),
    }
    if native_profiles:
        result["native_profile"] = {
            name: (
                max(profile[name] for profile in native_profiles)
                if name == "maximum_weights_per_pixel"
                else sum(profile[name] for profile in native_profiles)
            )
            for name in native_profiles[0]
        }
        result["result_xxh3_128"] = {
            name: _xxh3_128(values)
            for name, values in reduced.items()
        }
    del batches, reduced, kernel, rays
    gc.collect()
    return result


def main():
    """Benchmark current square tiling against equal-area thread-count tiling."""
    arguments = _arguments()
    job = read_job(arguments.job)
    config = job.config_data
    scan = job.scan
    if arguments.profile_load_frames:
        frame_count = min(len(scan), arguments.profile_load_frames)
        started = perf_counter()
        total_bytes = 0
        def load(frame):
            return np.asarray(scan.get_raw_img(frame).img).nbytes

        with ThreadPoolExecutor(
            max_workers=arguments.profile_load_workers
        ) as executor:
            total_bytes = sum(executor.map(load, range(frame_count)))
        elapsed = perf_counter() - started
        print(
            json.dumps(
                {
                    "frames": frame_count,
                    "workers": arguments.profile_load_workers,
                    "seconds": elapsed,
                    "frames_per_second": frame_count / elapsed,
                    "megabytes_per_second": total_bytes / elapsed / 1e6,
                },
                sort_keys=True,
            )
        )
        return
    payload = scan.get_raw_img(arguments.frame)
    image = np.asarray(payload.img)
    bounds = scan.exposure_angle_bounds(
        config,
        fallback=job.angle_fallback,
    )
    angles = np.ascontiguousarray(bounds[arguments.frame], dtype=np.float64)
    assets = _load_assets(job)
    spec = job.internal_spec()
    rows, columns = image.shape
    cases = []
    for threads in arguments.threads:
        tile_sizes = arguments.tile_sizes or [1024]
        work_blocks = arguments.work_blocks or [spec.work_block_pixels]
        depths = arguments.depths or [spec.max_depth]
        for tile_size in tile_sizes:
            if "current" in arguments.layouts:
                for work_block in work_blocks:
                    for depth in depths:
                        cases.append(
                            (
                                f"fixed-{tile_size}/{threads}/"
                                f"block-{work_block}/depth-{depth}",
                                _fixed_tiles(
                                    rows,
                                    columns,
                                    tile_size,
                                    tile_size,
                                ),
                                threads,
                                work_block,
                                depth,
                            )
                        )
        if "equal" in arguments.layouts and arguments.tile_sizes is None:
            for work_block in work_blocks:
                for depth in depths:
                    cases.append(
                        (
                            f"equal-{threads}/{threads}/"
                            f"block-{work_block}/depth-{depth}",
                            _equal_tiles(rows, columns, threads),
                            threads,
                            work_block,
                            depth,
                        )
                    )
    if arguments.parallel_frames:
        frame_indices = list(
            range(
                arguments.frame,
                min(len(scan), arguments.frame + arguments.parallel_frames),
            )
        )
        payloads = {
            frame: scan.get_raw_img(frame)
            for frame in frame_indices
        }
        for payload in payloads.values():
            np.asarray(payload.img)
        original_write = reconstruction._write_parquet
        original_checksum = reconstruction._uri_checksum_and_size
        reconstruction._write_parquet = lambda *args, **kwargs: None
        reconstruction._uri_checksum_and_size = lambda *args, **kwargs: ("", 0)
        try:
            parallel_results = []
            for name, tiles, threads, work_block, depth in cases:
                case_spec = replace(
                    spec,
                    threads=24,
                    work_block_pixels=work_block,
                    max_depth=depth,
                )
                rays = {
                    tile: _detector_corner_rays(config.detector, tile)
                    for tile in tiles
                }
                fingerprints = {
                    tile: _xxh3_128(values)
                    for tile, values in rays.items()
                }
                worker_count = max(1, 24 // threads)
                provenance = {}
                correct = _correction_pipeline(
                    config, scan, assets, provenance
                )
                started = perf_counter()

                def map_frame(frame):
                    return _map_frame_range(
                        case_spec,
                        scan,
                        config.detector,
                        config.ub_calculator,
                        (frame, frame + 1),
                        tiles,
                        bounds[frame : frame + 1],
                        Path("."),
                        correction_pipeline=correct,
                        job_digest=job.digest,
                        image_payloads={frame: payloads[frame]},
                        corner_rays=rays,
                        corner_rays_fingerprints=fingerprints,
                        kernel_threads=threads,
                        kernel_memory_budget_bytes=spec.memory_budget_bytes,
                        accumulation_budget_bytes=spec.memory_budget_bytes,
                    )

                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    manifests = list(executor.map(map_frame, frame_indices))
                elapsed = perf_counter() - started
                result = {
                    "name": name,
                    "threads": threads,
                    "image_workers": worker_count,
                    "tiles": len(tiles),
                    "frames": len(frame_indices),
                    "seconds": elapsed,
                    "frames_per_second": len(frame_indices) / elapsed,
                    "records": sum(
                        sum(partition.rows for partition in manifest.partitions)
                        for manifest in manifests
                    ),
                }
                parallel_results.append(result)
                print(json.dumps(result, sort_keys=True), flush=True)
                del manifests, rays
                gc.collect()
            if arguments.output:
                arguments.output.write_text(
                    json.dumps(parallel_results, indent=2, sort_keys=True),
                    encoding="utf-8",
                )
            return
        finally:
            reconstruction._write_parquet = original_write
            reconstruction._uri_checksum_and_size = original_checksum
    results = []
    for name, tiles, threads, work_block, depth in cases:
        result = _run_case(
            name=name,
            tiles=tiles,
            threads=threads,
            spec=replace(
                spec,
                threads=24,
                work_block_pixels=work_block,
                max_depth=depth,
            ),
            job=job,
            config=config,
            payload=payload,
            image=image,
            angles=angles,
            assets=assets,
            native_profile=arguments.native_profile,
        )
        results.append(result)
        print(json.dumps(result, sort_keys=True), flush=True)
    if arguments.output:
        arguments.output.write_text(
            json.dumps(results, indent=2, sort_keys=True),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()

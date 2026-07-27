"""Benchmark reconstruction map-task and Parquet storage controls."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np

from orgui.datautils.xrayutils.reconstruction import (
    _detector_corner_rays,
    _map_frame_range,
    _xxh3_128,
)
from orgui.reconstruction_job import (
    _correction_pipeline,
    _load_assets,
    read_job,
)

from .benchmark_reconstruction_mapping import _fixed_tiles


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--frames", type=int, default=16)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--tile-size", type=int, default=1024)
    parser.add_argument("--work-block", type=int, default=65536)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def _cases():
    return (
        ("baseline", 4, 1024, 2000),
        ("frame-batch-1", 1, 1024, 2000),
        ("frame-batch-16", 16, 1024, 2000),
        ("accumulation-64", 16, 64, 2000),
        ("accumulation-256", 16, 256, 2000),
        ("accumulation-2048", 16, 2048, 2000),
        ("partition-span-256", 4, 1024, 256),
        ("partition-span-1000", 4, 1024, 1000),
        ("partition-span-4000", 4, 1024, 4000),
    )


def main():
    """Run controlled map/Parquet cases while retaining 24 total threads."""
    arguments = _arguments()
    job = read_job(arguments.job)
    config = job.config_data
    scan = job.scan
    assets = _load_assets(job)
    frame_count = min(arguments.frames, len(scan))
    frame_indices = list(range(frame_count))
    payloads = {frame: scan.get_raw_img(frame) for frame in frame_indices}
    for payload in payloads.values():
        np.asarray(payload.img)
    bounds = scan.exposure_angle_bounds(
        config,
        fallback=job.angle_fallback,
    )
    rows, columns = config.detector.detector.shape
    tiles = _fixed_tiles(
        rows,
        columns,
        arguments.tile_size,
        arguments.tile_size,
    )
    rays = {
        tile: _detector_corner_rays(config.detector, tile)
        for tile in tiles
    }
    fingerprints = {
        tile: _xxh3_128(values)
        for tile, values in rays.items()
    }
    root = arguments.output or Path("build/reconstruction-storage-benchmark")
    root.mkdir(parents=True, exist_ok=True)
    results = []
    for name, frame_batch, accumulation_mib, partition_span in _cases():
        ranges = [
            (start, min(start + frame_batch, frame_count))
            for start in range(0, frame_count, frame_batch)
        ]
        spec = replace(
            job.internal_spec(),
            threads=24,
            work_block_pixels=arguments.work_block,
            partition_chunk_span=partition_span,
        )
        provenance = {}
        correct = _correction_pipeline(
            config,
            scan,
            assets,
            provenance,
        )
        started = perf_counter()
        with TemporaryDirectory(prefix=f"{name}-", dir=root) as output:

            def map_range(frame_range):
                return _map_frame_range(
                    spec,
                    scan,
                    config.detector,
                    config.ub_calculator,
                    frame_range,
                    tiles,
                    bounds[frame_range[0] : frame_range[1]],
                    output,
                    correction_pipeline=correct,
                    job_digest=job.digest,
                    image_payloads={
                        frame: payloads[frame]
                        for frame in range(*frame_range)
                    },
                    corner_rays=rays,
                    corner_rays_fingerprints=fingerprints,
                    kernel_threads=arguments.threads,
                    kernel_memory_budget_bytes=spec.memory_budget_bytes,
                    accumulation_budget_bytes=accumulation_mib * 1024**2,
                )

            workers = min(
                len(ranges),
                max(1, 24 // arguments.threads),
            )
            with ThreadPoolExecutor(max_workers=workers) as executor:
                manifests = list(executor.map(map_range, ranges))
            elapsed = perf_counter() - started
            partitions = [
                partition
                for manifest in manifests
                for partition in manifest.partitions
            ]
            result = {
                "name": name,
                "frame_batch": frame_batch,
                "accumulation_MiB": accumulation_mib,
                "partition_span": partition_span,
                "workers": workers,
                "seconds": elapsed,
                "frames_per_second": frame_count / elapsed,
                "parquet_files": len(partitions),
                "parquet_MiB": sum(
                    partition.size_bytes or 0 for partition in partitions
                )
                / 1024**2,
                "segments": sum(
                    int(manifest.metadata["accumulation_segments"])
                    for manifest in manifests
                ),
            }
            results.append(result)
            print(json.dumps(result, sort_keys=True), flush=True)
    (root / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

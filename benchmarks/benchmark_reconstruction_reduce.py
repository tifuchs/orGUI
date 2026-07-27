"""Benchmark bounded reciprocal-space reduction on an existing job."""

from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
from time import perf_counter

from orgui.datautils.xrayutils import reconstruction


def _arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("job", type=Path)
    parser.add_argument("--chunks", type=int, default=64)
    parser.add_argument("--memory-gib", type=float, default=16.0)
    return parser.parse_args()


def main():
    """Stream and reduce representative chunks from the largest map bucket."""
    arguments = _arguments()
    job = json.loads(arguments.job.read_text(encoding="utf-8"))
    manifests = [
        reconstruction._read_manifest(path)
        for path in job["map_manifests"]
    ]
    grouped = defaultdict(list)
    for manifest in manifests:
        for partition in manifest.partitions:
            grouped[(partition.grid_name, partition.bucket)].append(partition)
    group_key, partitions = max(
        grouped.items(),
        key=lambda item: sum(partition.rows for partition in item[1]),
    )
    grid_name, bucket = group_key
    spec = reconstruction._ReconstructionSpec.from_dict(manifests[0].spec)
    grid = next(grid for grid in spec.grids if grid.grid_name == grid_name)
    memory_bytes = int(arguments.memory_gib * 1024**3)
    batch_rows = max(
        4096,
        min(
            131072,
            memory_bytes // (max(1, len(partitions)) * 48 * 4),
        ),
    )
    readers = [
        reconstruction._ParquetRangeReader(
            partition.uri,
            batch_size=batch_rows,
        )
        for partition in sorted(partitions, key=lambda item: item.uri)
    ]
    chunk_grid = tuple(
        math.ceil(size / chunk)
        for size, chunk in zip(grid.shape, grid.chunk_shape)
    )
    chunk_start = bucket * spec.partition_chunk_span
    chunk_stop = min(
        math.prod(chunk_grid),
        chunk_start + spec.partition_chunk_span,
    )
    metadata_minima = []
    for reader in readers:
        for row_group in range(reader.parquet.metadata.num_row_groups):
            statistics = reader.parquet.metadata.row_group(row_group).column(
                reader.chunk_column
            ).statistics
            if statistics is not None and statistics.has_min_max:
                metadata_minima.append(int(statistics.min))
    if metadata_minima:
        chunk_start = max(chunk_start, min(metadata_minima))
    chunk_stop = min(chunk_stop, chunk_start + arguments.chunks)

    input_rows = 0
    output_rows = 0
    started = perf_counter()
    for ordinal, chunk_id in enumerate(
        range(chunk_start, chunk_stop),
        start=1,
    ):
        coordinates = reconstruction._chunk_coordinates(chunk_id, grid)
        local_shape = tuple(
            min(chunk, size - coordinate * chunk)
            for coordinate, chunk, size in zip(
                coordinates,
                grid.chunk_shape,
                grid.shape,
            )
        )
        local_stop = (
            (local_shape[0] - 1)
            * grid.chunk_shape[1]
            * grid.chunk_shape[2]
            + (local_shape[1] - 1) * grid.chunk_shape[2]
            + local_shape[2]
        )
        levels = []
        for reader in readers:
            batch = reader.read(chunk_id, 0, local_stop)
            input_rows += int(batch["chunk_id"].size)
            if not batch["chunk_id"].size:
                continue
            level = 0
            while level < len(levels) and levels[level] is not None:
                batch = reconstruction._merge_sorted_batches(
                    levels[level], batch
                )
                levels[level] = None
                level += 1
            if level == len(levels):
                levels.append(batch)
            else:
                levels[level] = batch
        reduced = reconstruction._empty_batch()
        for batch in reversed(levels):
            if batch is not None:
                reduced = reconstruction._merge_sorted_batches(reduced, batch)
        output_rows += int(reduced["chunk_id"].size)
        if ordinal % 8 == 0 or chunk_id + 1 == chunk_stop:
            elapsed = perf_counter() - started
            print(
                f"{ordinal}/{chunk_stop - chunk_start} chunks; "
                f"{input_rows / elapsed:,.0f} input rows/s; "
                f"{input_rows * 48 / elapsed / 1024**2:,.1f} MiB/s"
            )
    elapsed = perf_counter() - started
    print(
        f"bucket={bucket}, partitions={len(partitions)}, "
        f"batch_rows={batch_rows}, input_rows={input_rows}, "
        f"output_rows={output_rows}, seconds={elapsed:.3f}"
    )


if __name__ == "__main__":
    main()

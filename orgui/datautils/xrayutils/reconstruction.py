"""Out-of-core reciprocal-space reconstruction.

The numerical hot path is implemented by
``_reciprocal_reconstruction_cpp``.  This module owns the scientific Python
boundary, uncertainty propagation, immutable Parquet task products, and final
NeXus/HDF5 serialization.

Momentum-transfer grids use ``Angstrom^-1``.  HKL grids use reciprocal lattice
units.  All diffractometer angles accepted here are in radians.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field
from hashlib import sha256
import importlib
import json
import math
import os
from pathlib import Path
from queue import Queue
import threading
import time
from typing import Any

import h5py
import numpy as np


_PARTIAL_COLUMNS = (
    "chunk_id",
    "local_voxel_id",
    "weighted_intensity",
    "weighted_variance",
    "weight",
    "contributors",
)
_Q_FRAMES = {"lab", "alpha", "omega", "chi", "phi", "crystal"}
_ALL_FRAMES = _Q_FRAMES | {"hkl"}
_MIN_REDUCER_WORKER_MEMORY = 64 * 1024**2
_CHECKPOINT_BYTES_PER_ROW = 48
"""Fixed on-the-wire row size for the six ``_PARTIAL_COLUMNS`` fields
(two uint64 key columns, three float64 value columns, one uint64 count
column -- 6 * 8 bytes). Used to convert calibration-probe record counts
into byte estimates for the checkpoint file-count formula."""


def _triple(values, name, dtype=float):
    result = tuple(dtype(value) for value in values)
    if len(result) != 3:
        raise ValueError(f"{name} must contain exactly three values")
    return result


@dataclass(frozen=True)
class _GridSpec:
    """Describe a regular reciprocal-space grid.

    :param minimum:
        Lower voxel edges in r.l.u. for ``hkl`` or ``Angstrom^-1`` for Q.
    :param maximum:
        Exclusive upper grid bounds in the same units as ``minimum``.
    :param step:
        Voxel widths in the same units as ``minimum``.
    :param frame:
        ``hkl`` or one of ``lab``, ``alpha``, ``omega``, ``chi``, ``phi``,
        and ``crystal``.
    :param chunk_shape:
        HDF5 chunk shape in voxels.
    :param name:
        Optional HDF5 group name.
    """

    minimum: tuple[float, float, float]
    maximum: tuple[float, float, float]
    step: tuple[float, float, float]
    frame: str
    chunk_shape: tuple[int, int, int] = (64, 64, 64)
    name: str | None = None

    def __post_init__(self):
        object.__setattr__(self, "minimum", _triple(self.minimum, "minimum"))
        object.__setattr__(self, "maximum", _triple(self.maximum, "maximum"))
        object.__setattr__(self, "step", _triple(self.step, "step"))
        object.__setattr__(
            self, "chunk_shape", _triple(self.chunk_shape, "chunk_shape", int)
        )
        frame = self.frame.removeprefix("q_").lower()
        if frame not in _ALL_FRAMES:
            raise ValueError(f"Unsupported reciprocal-space frame: {self.frame}")
        object.__setattr__(self, "frame", frame)
        if any(not math.isfinite(value) for value in self.minimum + self.maximum):
            raise ValueError("Grid bounds must be finite")
        if any(value <= 0 or not math.isfinite(value) for value in self.step):
            raise ValueError("Grid steps must be finite and positive")
        if any(upper <= lower for lower, upper in zip(self.minimum, self.maximum)):
            raise ValueError("Each grid maximum must exceed its minimum")
        if any(value <= 0 for value in self.chunk_shape):
            raise ValueError("Chunk dimensions must be positive")
        effective_chunk = tuple(
            min(size, chunk)
            for size, chunk in zip(self.shape, self.chunk_shape)
        )
        chunk_bytes = math.prod(effective_chunk) * np.dtype(np.float64).itemsize
        if chunk_bytes >= 2**32:
            raise ValueError(
                "The effective HDF5 chunk must be smaller than 4 GiB; "
                f"{effective_chunk} requires {chunk_bytes / 1024**3:.2f} GiB "
                "per float64 dataset"
            )
        if self.name is not None and not self.name:
            raise ValueError("Grid name cannot be empty")

    @property
    def shape(self) -> tuple[int, int, int]:
        """Grid shape in voxels."""
        return tuple(
            int(math.ceil((upper - lower) / width))
            for lower, upper, width in zip(self.minimum, self.maximum, self.step)
        )

    @property
    def effective_maximum(self) -> tuple[float, float, float]:
        """Upper edge implied by the integer grid shape."""
        return tuple(
            lower + size * width
            for lower, size, width in zip(self.minimum, self.shape, self.step)
        )

    @property
    def grid_name(self) -> str:
        """Stable HDF5-safe name for this grid."""
        raw = self.name or ("hkl" if self.frame == "hkl" else f"q_{self.frame}")
        return "".join(char if char.isalnum() or char in "_-" else "_" for char in raw)


@dataclass(frozen=True)
class _ReconstructionSpec:
    """Configuration shared by mapping, reduction, and finalization."""

    grids: tuple[_GridSpec, ...]
    max_depth: int = 2
    threads: int = 1
    work_block_pixels: int = 4096
    memory_budget_bytes: int = 512 * 1024 * 1024
    partition_chunk_span: int = 256
    compression: str = "bitshuffle-lz4"
    infer_angle_bounds: bool = False
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        grids = tuple(
            grid if isinstance(grid, _GridSpec) else _GridSpec(**grid)
            for grid in self.grids
        )
        object.__setattr__(self, "grids", grids)
        if not grids:
            raise ValueError("At least one output grid is required")
        names = [grid.grid_name for grid in grids]
        if len(names) != len(set(names)):
            raise ValueError("Grid names must be unique")
        if not 0 <= self.max_depth <= 8:
            raise ValueError("max_depth must be between 0 and 8")
        if self.threads < 1:
            raise ValueError("threads must be positive")
        if self.work_block_pixels < 1:
            raise ValueError("work_block_pixels must be positive")
        if self.memory_budget_bytes < 1024 * 1024:
            raise ValueError("memory_budget_bytes must be at least 1 MiB")
        if self.partition_chunk_span < 1:
            raise ValueError("partition_chunk_span must be positive")

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable representation."""
        result = asdict(self)
        result["metadata"] = dict(self.metadata)
        return result

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]):
        """Create a specification from JSON-compatible values."""
        data = dict(values)
        data["grids"] = tuple(_GridSpec(**grid) for grid in data["grids"])
        return cls(**data)

    @property
    def digest(self) -> str:
        """SHA-256 digest identifying scientifically equivalent settings."""
        encoded = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), default=str
        ).encode("utf-8")
        return sha256(encoded).hexdigest()


@dataclass
class _PartitionFile:
    """One immutable Parquet object produced by a mapping task."""

    grid_name: str
    bucket: int
    uri: str
    rows: int
    checksum: str
    size_bytes: int | None = None


@dataclass
class _ChunkFile:
    """One reduced Parquet shard belonging to one output spatial chunk."""

    grid_name: str
    chunk_id: int
    shard_start: int
    shard_stop: int
    uri: str
    rows: int
    checksum: str
    size_bytes: int | None = None


@dataclass
class _TaskManifest:
    """Serializable, restart-safe mapping or reduction task manifest."""

    kind: str
    task_id: str
    spec_hash: str
    status: str
    spec: dict[str, Any]
    frame_range: tuple[int, int] | None = None
    detector_tile: tuple[int, int, int, int] | None = None
    partitions: list[_PartitionFile] = field(default_factory=list)
    chunks: list[_ChunkFile] = field(default_factory=list)
    source_tasks: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible manifest."""
        result = asdict(self)
        return result

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]):
        """Deserialize a task manifest."""
        data = dict(values)
        data["partitions"] = [
            value if isinstance(value, _PartitionFile) else _PartitionFile(**value)
            for value in data.get("partitions", [])
        ]
        data["chunks"] = [
            value if isinstance(value, _ChunkFile) else _ChunkFile(**value)
            for value in data.get("chunks", [])
        ]
        if data.get("frame_range") is not None:
            data["frame_range"] = tuple(data["frame_range"])
        if data.get("detector_tile") is not None:
            data["detector_tile"] = tuple(data["detector_tile"])
        return cls(**data)


def _native_module():
    try:
        return importlib.import_module(
            "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
        )
    except ImportError as exc:
        raise RuntimeError(
            "The native reciprocal-space reconstruction extension is unavailable. "
            "Install orGUI from a wheel containing its C++ extensions or rebuild it "
            "with Meson."
        ) from exc


def _xxh3_128(array) -> str:
    """Return the XXH3-128 fingerprint of one contiguous array or buffer."""
    value = np.ascontiguousarray(array)
    return _native_module().xxh3_128(value)


def _detector_corner_rays(detector, detector_tile) -> np.ndarray:
    """Calculate calibrated outgoing unit rays at detector-pixel corners.

    :param DetectorCalibration.Detector2D_SXRD detector:
        Calibrated detector geometry.
    :param detector_tile:
        ``(row_start, row_stop, column_start, column_stop)`` using exclusive
        stop indices.
    :returns:
        C-contiguous array shaped ``(rows + 1, columns + 1, 3)`` in the
        laboratory frame.
    """
    row_start, row_stop, column_start, column_stop = map(int, detector_tile)
    if row_start < 0 or column_start < 0:
        raise ValueError("Detector tile starts must be non-negative")
    if row_stop <= row_start or column_stop <= column_start:
        raise ValueError("Detector tile stops must exceed starts")
    rows = np.arange(row_start, row_stop + 1, dtype=np.float64) - 0.5
    columns = np.arange(column_start, column_stop + 1, dtype=np.float64) - 0.5
    row_grid, column_grid = np.meshgrid(rows, columns, indexing="ij")
    gamma_p, delta_p = detector.primBeamPoints(row_grid, column_grid)
    cosine_gamma = np.cos(gamma_p)
    rays = np.empty((*row_grid.shape, 3), dtype=np.float64)
    rays[..., 0] = np.sin(delta_p) * cosine_gamma
    rays[..., 1] = np.cos(delta_p) * cosine_gamma
    rays[..., 2] = np.sin(gamma_p)
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)
    return np.ascontiguousarray(rays)


def _kernel_for_grid(
    spec,
    grid,
    ub_calculator,
    *,
    threads=None,
    memory_budget_bytes=None,
):
    ub = np.asarray(ub_calculator.getUB(), dtype=np.float64)
    u = np.asarray(ub_calculator.getU(), dtype=np.float64)
    return _native_module().ReconstructionKernel(
        np.asarray(grid.minimum, dtype=np.float64),
        np.asarray(grid.step, dtype=np.float64),
        np.asarray(grid.shape, dtype=np.int64),
        np.asarray(grid.chunk_shape, dtype=np.int64),
        grid.frame,
        float(ub_calculator.getK()),
        np.ascontiguousarray(np.linalg.inv(ub)),
        np.ascontiguousarray(np.linalg.inv(u)),
        spec.max_depth,
        spec.threads if threads is None else threads,
        spec.work_block_pixels,
        (
            spec.memory_budget_bytes
            if memory_budget_bytes is None
            else memory_budget_bytes
        ),
    )


def _files_per_job(job_data_bytes, ram_budget_bytes, checkpoint_count) -> int:
    """Minimum number of resumable checkpoint files for one job.

    ``checkpoint_count`` is a floor, not a target: the result exceeds it
    whenever the estimated data volume would not otherwise fit the memory
    budget across that many files.

    :param float job_data_bytes:
        Estimated **reduced** (post in-frame-dedup) record volume for the
        job's frame-range slice. Not raw per-pixel/per-corner records.
    :param float ram_budget_bytes:
        Combined memory budget for the whole job (all threads, one node) --
        the user's own resource request, never a hardcoded value.
    :param int checkpoint_count:
        User-requested minimum number of checkpoint files.
    :returns:
        ``max(checkpoint_count, ceil(job_data_bytes / ram_budget_bytes))``.
    :rtype: int
    :raises ValueError:
        If ``checkpoint_count`` or ``ram_budget_bytes`` are not positive, or
        ``job_data_bytes`` is negative.
    """
    checkpoint_count = int(checkpoint_count)
    if checkpoint_count < 1:
        raise ValueError("checkpoint_count must be at least one")
    if ram_budget_bytes <= 0:
        raise ValueError("ram_budget_bytes must be positive")
    if job_data_bytes < 0:
        raise ValueError("job_data_bytes must not be negative")
    files_for_memory = math.ceil(job_data_bytes / ram_budget_bytes)
    return max(checkpoint_count, files_for_memory)


def _calibration_probe(
    kernel,
    mask,
    corner_rays,
    angles_start,
    angles_end,
    *,
    budget_seconds=0.1,
    kernel_threads=1,
    bootstrap_tile=16,
    sample_tiles=9,
    min_sample_pixels=256,
    max_sample_pixels=2_000_000,
    rng=None,
):
    """Estimate per-pixel mapping cost and record volume, live, in-budget.

    Runs the native kernel on small tiles scattered across the real
    detector mask with dummy (zero) intensity/variance: coordinate-
    evaluation cost and voxel occupancy depend only on geometry and the
    mask, not on pixel values, so no image data or disk I/O is needed.
    Bounded to ``budget_seconds`` wall time (a two-pass design: a tiny
    bootstrap tile first estimates the per-pixel rate, then a properly
    sized, stratified sample spends the remaining budget), so it is safe
    to call live, e.g. from the GUI, on every settings change.

    The returned rate is specific to the machine, job settings, and sample
    this call happened to draw -- it must be measured fresh for every job,
    never cached or hardcoded across jobs or machines.

    :param kernel:
        A constructed ``ReconstructionKernel`` (already encodes
        ``max_depth`` and its native thread count).
    :param np.ndarray mask:
        Boolean mask for the full detector, shape ``(rows, columns)``.
        ``True`` marks excluded pixels, matching the kernel's convention.
    :param np.ndarray corner_rays:
        Full-detector corner rays, shape ``(rows + 1, columns + 1, 3)``.
    :param np.ndarray angles_start:
        Shape ``(4,)`` diffractometer angles in radians at exposure start.
    :param np.ndarray angles_end:
        Shape ``(4,)`` diffractometer angles in radians at exposure end.
    :param float budget_seconds:
        Wall-time budget for the whole probe, including the bootstrap pass.
    :param int kernel_threads:
        Recorded in the result as bookkeeping only; does not affect
        ``kernel``'s own configured thread count.
    :param int bootstrap_tile:
        Side length in pixels of the bootstrap tile.
    :param int sample_tiles:
        Number of scattered tiles sampled in the sized pass.
    :param int min_sample_pixels, max_sample_pixels:
        Clamp on the sized pass's total target sample size.
    :param rng:
        Optional :class:`numpy.random.Generator` for deterministic tests;
        defaults to a fresh, unseeded generator.
    :returns:
        Dict with ``kernel_threads``, ``sampled_pixels``,
        ``sampled_seconds``, ``seconds_per_pixel`` (kernel CPU-seconds per
        sampled pixel), and ``records_per_pixel`` (post-dedup records per
        sampled pixel).
    :rtype: dict
    """
    if rng is None:
        rng = np.random.default_rng()
    rows, columns = mask.shape
    started = time.perf_counter()

    def sample_tile(row_start, row_stop, column_start, column_stop):
        tile_shape = (row_stop - row_start, column_stop - column_start)
        intensity = np.zeros(tile_shape, dtype=np.float64)
        variance = np.zeros(tile_shape, dtype=np.float64)
        tile_mask = np.ascontiguousarray(
            mask[row_start:row_stop, column_start:column_stop]
        )
        tile_rays = np.ascontiguousarray(
            corner_rays[row_start : row_stop + 1, column_start : column_stop + 1]
        )
        result = kernel.accumulate(
            intensity, variance, tile_mask, tile_rays,
            angles_start, angles_end, True,
        )
        profile = result["_profile"]
        pixels = tile_shape[0] * tile_shape[1]
        return (
            pixels,
            float(profile["block_mapping_cpu_seconds"]),
            int(profile["reduced_block_records"]),
        )

    bootstrap_side = max(1, min(int(bootstrap_tile), rows, columns))
    total_pixels, total_seconds, total_records = sample_tile(
        0, bootstrap_side, 0, bootstrap_side
    )
    elapsed = time.perf_counter() - started
    rate = total_seconds / max(1, total_pixels)

    remaining_budget = budget_seconds - elapsed
    if remaining_budget > 0 and rate > 0:
        target_pixels = int(
            min(max(remaining_budget / rate, min_sample_pixels), max_sample_pixels)
        )
        tile_side = max(1, int(math.sqrt(target_pixels / max(1, sample_tiles))))
        tile_side = min(tile_side, rows, columns)
        for _ in range(int(sample_tiles)):
            if time.perf_counter() - started >= budget_seconds:
                break
            row_start = int(rng.integers(0, max(1, rows - tile_side + 1)))
            column_start = int(rng.integers(0, max(1, columns - tile_side + 1)))
            pixels, seconds, records = sample_tile(
                row_start, row_start + tile_side,
                column_start, column_start + tile_side,
            )
            total_pixels += pixels
            total_seconds += seconds
            total_records += records

    return {
        "kernel_threads": int(kernel_threads),
        "sampled_pixels": total_pixels,
        "sampled_seconds": total_seconds,
        "seconds_per_pixel": total_seconds / max(1, total_pixels),
        "records_per_pixel": total_records / max(1, total_pixels),
    }


def _calibration_probe_all_grids(
    spec,
    ub_calculator,
    mask,
    corner_rays,
    angles_start,
    angles_end,
    *,
    budget_seconds=0.1,
    kernel_threads=None,
):
    """Run :func:`_calibration_probe` once per grid in ``spec.grids``.

    Reuses the same sampled-tile geometry (mask, corner rays, angles)
    across every grid -- only the kernel differs per grid, so no new
    geometry sampling is needed per grid (see the design doc's multi-grid
    calibration-probe extension). The overall ``budget_seconds`` is split
    evenly across grids.

    :returns:
        Dict keyed by ``grid.grid_name``, values are
        :func:`_calibration_probe` result dicts.
    :rtype: dict[str, dict]
    """
    per_grid_budget = budget_seconds / max(1, len(spec.grids))
    results = {}
    for grid in spec.grids:
        kernel = _kernel_for_grid(spec, grid, ub_calculator, threads=kernel_threads)
        results[grid.grid_name] = _calibration_probe(
            kernel,
            mask,
            corner_rays,
            angles_start,
            angles_end,
            budget_seconds=per_grid_budget,
            kernel_threads=(
                kernel_threads if kernel_threads is not None else spec.threads
            ),
        )
    return results


def _empty_batch() -> dict[str, np.ndarray]:
    return {
        "chunk_id": np.empty(0, dtype=np.uint64),
        "local_voxel_id": np.empty(0, dtype=np.uint64),
        "weighted_intensity": np.empty(0, dtype=np.float64),
        "weighted_variance": np.empty(0, dtype=np.float64),
        "weight": np.empty(0, dtype=np.float64),
        "contributors": np.empty(0, dtype=np.uint64),
    }


def _tree_insert(levels, batch):
    """Insert one already-reduced batch into a binary-counter merge tree.

    Mutates ``levels`` in place: descends occupied levels, merging and
    cascading exactly as this insert step always has inside
    :func:`_reduce_batches`, so calling this once per batch incrementally
    (as batches become available over time) is behaviorally identical to
    building the whole ``levels`` list from a pre-collected batch list up
    front. This is the mechanism :class:`_CheckpointAccumulator` reuses to
    merge frames into a checkpoint as they finish, rather than only once at
    the end.

    :param list levels:
        Mutable list of ``batch | None``, one slot per tree level.
    :param batch:
        A native record batch (dict of the six ``_PARTIAL_COLUMNS``
        arrays), already reduced/sorted.
    """
    if not np.asarray(batch["chunk_id"]).size:
        return
    level = 0
    while level < len(levels) and levels[level] is not None:
        batch = _merge_sorted_batches(levels[level], batch)
        levels[level] = None
        level += 1
    if level == len(levels):
        levels.append(batch)
    else:
        levels[level] = batch


def _tree_finalize(levels):
    """Fold whatever ``levels`` remain into one merged batch.

    :param list levels:
        The same ``levels`` list :func:`_tree_insert` maintains.
    :returns:
        One merged, sorted, deduplicated native record batch.
    """
    result = _empty_batch()
    for batch in reversed(levels):
        if batch is not None:
            result = _merge_sorted_batches(result, batch)
    return result


def _reduce_batches(batches: Iterable[Mapping[str, np.ndarray]]):
    levels = []
    for batch in batches:
        _tree_insert(levels, batch)
    return _tree_finalize(levels)


def _merge_sorted_batches(left, right):
    """Linearly merge two sorted, already-reduced native record batches."""
    if not left["chunk_id"].size:
        return right
    if not right["chunk_id"].size:
        return left
    arguments = [
        np.ascontiguousarray(batch[name])
        for batch in (left, right)
        for name in _PARTIAL_COLUMNS
    ]
    return _native_module().merge_sorted_batches(*arguments)


class _CheckpointAccumulator:
    """In-memory tree-merge accumulator for one checkpoint's records.

    Not grid-aware: callers key a ``dict[(checkpoint_index, grid_name),
    _CheckpointAccumulator]`` externally (design doc Sec9/Sec12) -- this
    class only tracks one grid's records for one checkpoint. Thread-safe:
    :meth:`insert` may be called concurrently by multiple worker threads
    finishing different frames; the lock protects only the amortized
    O(log N) binary-counter bookkeeping (:func:`_tree_insert`), not the
    expensive per-frame kernel work that happens before a caller ever
    calls :meth:`insert` -- this is what keeps contention low even under
    many concurrent workers (design doc Sec9).
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._levels: list = []
        self._byte_total = 0

    def insert(self, batch):
        """Merge one already-reduced batch into this checkpoint's tree.

        :param batch:
            A native record batch (dict of the six ``_PARTIAL_COLUMNS``
            arrays), already reduced/sorted -- e.g. one frame's kernel
            output.
        """
        batch_bytes = sum(np.asarray(values).nbytes for values in batch.values())
        with self._lock:
            _tree_insert(self._levels, batch)
            self._byte_total += batch_bytes

    def should_flush(self, budget_bytes) -> bool:
        """Whether the running byte total has crossed ``budget_bytes``."""
        with self._lock:
            return self._byte_total >= budget_bytes

    def finalize(self):
        """Fold remaining levels into one merged batch, then reset.

        :returns:
            The folded batch -- see :func:`_tree_finalize`.
        """
        with self._lock:
            result = _tree_finalize(self._levels)
            self._levels = []
            self._byte_total = 0
        return result


def _checkpoint_batch_digest(batch) -> str:
    """XXH3-128 over the batch's six columns, concatenated in
    ``_PARTIAL_COLUMNS`` order. Shared by :func:`_write_checkpoint` (hashed
    while the data is already in memory, not by re-reading the file
    afterward) and :func:`_verify_checkpoint`'s full re-verify path (which
    must re-derive the digest the same way to compare against it)."""
    return _xxh3_128(
        np.concatenate(
            [
                np.ascontiguousarray(batch[name]).view(np.uint8)
                for name in _PARTIAL_COLUMNS
            ]
        )
    )


def _write_checkpoint(batch, path, *, spec_digest, metadata=None, chunk_rows=65536):
    """Write one reduced record batch as a checkpoint HDF5 file.

    Columnar layout (design doc Sec8): one dataset per ``_PARTIAL_COLUMNS``
    field, ``bitshuffle-lz4`` compression, chunked at ``chunk_rows`` rows --
    the measured best-performing configuration (Sec8), used unconditionally
    for checkpoints regardless of the job's own final-output
    ``compression`` setting, since checkpoints are internal scratch data.
    Checksummed inline while the data is already in memory (Sec8), not by
    re-reading the file afterward. Written atomically -- to ``.tmp``, then
    ``.replace()`` -- matching the pattern already used by
    :func:`_write_manifest` and :func:`_finalize_reconstruction`.

    :param batch:
        A native record batch (dict of the six ``_PARTIAL_COLUMNS``
        arrays).
    :param path:
        Destination path for the checkpoint file.
    :param str spec_digest:
        The job's ``_ReconstructionSpec.digest``, stored as an HDF5
        attribute so resume can detect a stale checkpoint by comparing
        digests (design doc Sec11) without a separate manifest registry.
    :param dict metadata:
        Optional extra JSON-serializable scalar attributes to store (e.g.
        checkpoint index, grid name, frame range).
    :param int chunk_rows:
        HDF5 chunk size in rows for every column dataset.
    :returns:
        The XXH3-128 hex digest written to the file.
    :rtype: str
    """
    rows = int(np.asarray(batch["chunk_id"]).size)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    digest = _checkpoint_batch_digest(batch)
    dataset_kwargs = (
        {
            "chunks": (min(int(chunk_rows), rows),),
            **_compression_kwargs("bitshuffle-lz4"),
        }
        if rows > 0
        else {}
    )
    with h5py.File(temporary, "w") as h5file:
        h5file.attrs["spec_sha256"] = spec_digest
        h5file.attrs["xxh3_128"] = digest
        h5file.attrs["rows"] = rows
        for key, value in (metadata or {}).items():
            h5file.attrs[key] = value
        for name in _PARTIAL_COLUMNS:
            h5file.create_dataset(
                name, data=np.ascontiguousarray(batch[name]), **dataset_kwargs
            )
        h5file.flush()
    temporary.replace(path)
    return digest


def _read_checkpoint(path) -> dict[str, np.ndarray]:
    """Read a checkpoint file's record batch back into the standard shape.

    :returns:
        Dict of the six ``_PARTIAL_COLUMNS`` arrays.
    :rtype: dict[str, np.ndarray]
    """
    with h5py.File(path, "r") as h5file:
        return {name: h5file[name][()] for name in _PARTIAL_COLUMNS}


def _verify_checkpoint(path, *, spec_digest, full=False) -> bool:
    """Check a checkpoint file's identity, optionally its data integrity.

    Cheap path (default): the file exists and its ``spec_sha256`` attribute
    matches ``spec_digest`` -- this is the entire resume mechanism (design
    doc Sec11), with no separate manifest registry, since atomic
    write-then-replace already guarantees a file present at its final name
    is complete.

    :param bool full:
        If true, also re-read every column and re-verify the stored
        XXH3-128 digest against freshly hashed data -- cheap at XXH3 speed
        even for real job sizes, so safe to do routinely, unlike the
        SHA-256-by-reread pattern this replaces (design doc Sec8).
    :returns:
        ``True`` if the checkpoint is present and matches.
    :rtype: bool
    """
    path = Path(path)
    if not path.exists():
        return False
    with h5py.File(path, "r") as h5file:
        if h5file.attrs.get("spec_sha256") != spec_digest:
            return False
        if full:
            batch = {name: h5file[name][()] for name in _PARTIAL_COLUMNS}
            if _checkpoint_batch_digest(batch) != h5file.attrs.get("xxh3_128"):
                return False
    return True


def _require_pyarrow():
    try:
        import pyarrow as pa
        import pyarrow.fs as pafs
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise RuntimeError(
            "Parquet reconstruction storage requires the optional 'pyarrow' "
            "dependency."
        ) from exc
    return pa, pafs, pq


def _filesystem_path(uri):
    _, pafs, _ = _require_pyarrow()
    text = os.fspath(uri)
    if "://" not in text:
        path = str(Path(text).absolute())
        return pafs.LocalFileSystem(), path
    return pafs.FileSystem.from_uri(text)


def _write_parquet(batch, uri, metadata):
    pa, _, pq = _require_pyarrow()
    filesystem, path = _filesystem_path(uri)
    parent = str(Path(path).parent).replace("\\", "/")
    try:
        filesystem.create_dir(parent, recursive=True)
    except (NotImplementedError, OSError):
        pass
    table = pa.table({name: batch[name] for name in _PARTIAL_COLUMNS})
    encoded_metadata = {
        str(key).encode("utf-8"): json.dumps(value, sort_keys=True).encode("utf-8")
        for key, value in metadata.items()
    }
    table = table.replace_schema_metadata(encoded_metadata)
    pq.write_table(table, path, filesystem=filesystem, compression="zstd")


def _read_parquet(uri):
    _, _, pq = _require_pyarrow()
    filesystem, path = _filesystem_path(uri)
    table = pq.read_table(path, filesystem=filesystem, columns=list(_PARTIAL_COLUMNS))
    return {name: table[name].to_numpy() for name in _PARTIAL_COLUMNS}


class _ParquetRangeReader:
    """Stream monotonically increasing key ranges from sorted Parquet."""

    def __init__(self, uri, *, batch_size=131072, use_threads=True):
        _, _, pq = _require_pyarrow()
        filesystem, path = _filesystem_path(uri)
        self.parquet = pq.ParquetFile(path, filesystem=filesystem)
        names = self.parquet.schema_arrow.names
        self.chunk_column = names.index("chunk_id")
        self.local_column = names.index("local_voxel_id")
        self.batch_size = int(batch_size)
        self.use_threads = bool(use_threads)
        self.iterator = None
        self.batch = None
        self.position = 0
        self.previous_key = None

    def _start(self, chunk_id, local_start):
        row_groups = []
        for index in range(self.parquet.metadata.num_row_groups):
            metadata = self.parquet.metadata.row_group(index)
            chunk_statistics = metadata.column(
                self.chunk_column
            ).statistics
            if (
                chunk_statistics is not None
                and chunk_statistics.has_min_max
                and int(chunk_statistics.max) < chunk_id
            ):
                continue
            if (
                chunk_statistics is not None
                and chunk_statistics.has_min_max
                and int(chunk_statistics.max) == chunk_id
            ):
                local_statistics = metadata.column(
                    self.local_column
                ).statistics
                if (
                    local_statistics is not None
                    and local_statistics.has_min_max
                    and int(local_statistics.max) < local_start
                ):
                    continue
            row_groups.append(index)
        self.iterator = self.parquet.iter_batches(
            batch_size=self.batch_size,
            row_groups=row_groups,
            columns=list(_PARTIAL_COLUMNS),
            use_threads=self.use_threads,
        )

    def _advance(self):
        try:
            record_batch = next(self.iterator)
        except StopIteration:
            self.batch = None
            return False
        self.batch = {
            name: record_batch.column(name).to_numpy(zero_copy_only=False)
            for name in _PARTIAL_COLUMNS
        }
        self.position = 0
        return True

    def read(self, chunk_id, local_start, local_stop):
        """Read one exact range without revisiting earlier Parquet data."""
        key = (int(chunk_id), int(local_start))
        if self.previous_key is not None and key < self.previous_key:
            raise ValueError("Parquet ranges must be requested in sorted order")
        self.previous_key = key
        if self.iterator is None:
            self._start(chunk_id, local_start)
        parts = []
        while self.batch is not None or self._advance():
            chunks = self.batch["chunk_id"]
            position = self.position
            if position >= chunks.size:
                self.batch = None
                continue
            position += int(
                np.searchsorted(
                    chunks[position:],
                    np.uint64(chunk_id),
                    side="left",
                )
            )
            if position >= chunks.size:
                self.batch = None
                continue
            current_chunk = int(chunks[position])
            if current_chunk > chunk_id:
                self.position = position
                break
            chunk_stop = int(
                np.searchsorted(
                    chunks,
                    np.uint64(chunk_id),
                    side="right",
                    sorter=None,
                )
            )
            local = self.batch["local_voxel_id"]
            start = position + int(
                np.searchsorted(
                    local[position:chunk_stop],
                    np.uint64(local_start),
                    side="left",
                )
            )
            stop = position + int(
                np.searchsorted(
                    local[position:chunk_stop],
                    np.uint64(local_stop),
                    side="left",
                )
            )
            if stop > start:
                parts.append(
                    {
                        name: values[start:stop]
                        for name, values in self.batch.items()
                    }
                )
            self.position = stop
            if stop < chunk_stop:
                break
            self.position = chunk_stop
            if chunk_stop < chunks.size:
                break
            self.batch = None
        if not parts:
            return _empty_batch()
        if len(parts) == 1:
            return parts[0]
        return {
            name: np.concatenate([part[name] for part in parts])
            for name in _PARTIAL_COLUMNS
        }


def _uri_checksum_and_size(uri):
    filesystem, path = _filesystem_path(uri)
    digest = sha256()
    size = 0
    with filesystem.open_input_file(path) as stream:
        while True:
            block = stream.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
            size += len(block)
    return digest.hexdigest(), size


def _uri_checksum(uri):
    return _uri_checksum_and_size(uri)[0]


def _uri_size(uri):
    filesystem, path = _filesystem_path(uri)
    size = int(filesystem.get_file_info(path).size)
    if size < 0:
        raise OSError(f"Could not determine scratch-file size: {uri}")
    return size


def _mark_scratch_verified(file, verification_cache):
    if verification_cache is not None:
        verification_cache.add((file.uri, file.checksum))


def _verify_scratch_file(file, verification_cache=None):
    key = (file.uri, file.checksum)
    if verification_cache is not None and key in verification_cache:
        return True
    if file.size_bytes is not None and _uri_size(file.uri) != file.size_bytes:
        raise OSError(f"Scratch-file size mismatch: {file.uri}")
    if _uri_checksum(file.uri) != file.checksum:
        raise OSError(f"Checksum mismatch for {file.uri}")
    if verification_cache is not None:
        verification_cache.add(key)
    return True


def _join_uri(base, *parts):
    base = os.fspath(base).rstrip("/\\")
    separator = "/" if "://" in base else os.sep
    return separator.join((base, *(str(part).strip("/\\") for part in parts)))


def _jsonable(value):
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, str | int | float | bool) or value is None:
        return value
    return repr(value)


def _write_manifest(manifest: _TaskManifest, uri):
    """Write a task manifest atomically on a local filesystem."""
    path = Path(uri)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(manifest.to_dict(), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    temporary.replace(path)
    return str(path)


def _read_manifest(value) -> _TaskManifest:
    """Read a manifest path or normalize an existing manifest."""
    if isinstance(value, _TaskManifest):
        return value
    if isinstance(value, Mapping):
        return _TaskManifest.from_dict(value)
    return _TaskManifest.from_dict(
        json.loads(Path(value).read_text(encoding="utf-8"))
    )


def _map_frame_range(
    spec: _ReconstructionSpec,
    scan,
    detector,
    ub_calculator,
    frame_range,
    detector_tiles,
    angle_bounds_rad,
    output_uri,
    *,
    correction_pipeline: Callable,
    job_digest: str,
    image_payloads: Mapping[int, object] | None = None,
    corner_rays: Mapping[tuple[int, int, int, int], np.ndarray] | None = None,
    corner_rays_fingerprints: Mapping[
        tuple[int, int, int, int], str
    ]
    | None = None,
    verification_cache: set[tuple[str, str]] | None = None,
    kernel_threads: int | None = None,
    kernel_memory_budget_bytes: int | None = None,
    accumulation_budget_bytes: int = 64 * 1024 * 1024,
    image_progress: Callable[[int, int, int], None] | None = None,
) -> _TaskManifest:
    """Map a streaming frame range and write byte-bounded Parquet partials.

    ``angle_bounds_rad`` must have shape ``(frames, 2, 4)`` in the order
    ``alpha, omega, chi, phi``. The central correction object may provide a
    ``correct_frame`` method, which is called once before detector tiling.
    This preserves full-detector pixel-repair neighborhoods and avoids
    repeating static corrections for every tile. The tile callback remains
    available for correction objects without that method.

    Every image is loaded once, then all detector tiles are processed. Reduced
    native records are retained until ``accumulation_budget_bytes`` would be
    exceeded, at which point they are reduced across images and written as one
    deterministic Parquet segment.

    ``image_payloads``, ``corner_rays``, ``corner_rays_fingerprints``,
    ``verification_cache``, ``kernel_threads``,
    ``kernel_memory_budget_bytes``, and ``accumulation_budget_bytes`` are local
    execution controls. They do not alter scientific values.

    ``image_progress`` is called after each completed image with the frame
    index, retained-record bytes, and number of flushed segments.
    """
    start, stop = map(int, frame_range)
    if start < 0 or stop <= start or stop > len(scan):
        raise ValueError("Invalid frame range")
    detector_tiles = tuple(
        tuple(map(int, detector_tile)) for detector_tile in detector_tiles
    )
    if not detector_tiles:
        raise ValueError("At least one detector tile is required")
    detector_rows, detector_columns = detector.detector.shape
    for row_start, row_stop, column_start, column_stop in detector_tiles:
        if (
            row_start < 0
            or column_start < 0
            or row_stop > detector_rows
            or column_stop > detector_columns
            or row_stop <= row_start
            or column_stop <= column_start
        ):
            raise ValueError("Invalid detector tile")
    bounds = np.asarray(angle_bounds_rad, dtype=np.float64)
    if bounds.shape != (stop - start, 2, 4):
        raise ValueError("angle_bounds_rad must have shape (frame_count, 2, 4)")
    if not np.all(np.isfinite(bounds)):
        raise ValueError("angle_bounds_rad contains non-finite values")
    if not callable(correction_pipeline):
        raise TypeError("correction_pipeline must be the central job correction")
    if not job_digest:
        raise ValueError("job_digest is required")
    if accumulation_budget_bytes < 1:
        raise ValueError("accumulation_budget_bytes must be positive")
    ray_arrays = {}
    ray_fingerprints = {}
    for detector_tile in detector_tiles:
        row_start, row_stop, column_start, column_stop = detector_tile
        rays = (
            _detector_corner_rays(detector, detector_tile)
            if corner_rays is None or detector_tile not in corner_rays
            else np.ascontiguousarray(
                corner_rays[detector_tile], dtype=np.float64
            )
        )
        expected_ray_shape = (
            row_stop - row_start + 1,
            column_stop - column_start + 1,
            3,
        )
        if rays.shape != expected_ray_shape:
            raise ValueError(
                f"corner_rays must have shape {expected_ray_shape}"
            )
        ray_arrays[detector_tile] = rays
        fingerprint = (
            None
            if corner_rays_fingerprints is None
            else corner_rays_fingerprints.get(detector_tile)
        )
        if fingerprint is None:
            fingerprint = _xxh3_128(rays)
        elif not isinstance(fingerprint, str) or len(fingerprint) != 32:
            raise ValueError(
                "Corner-ray fingerprints must be 32-character "
                "XXH3-128 digests"
            )
        ray_fingerprints[detector_tile] = fingerprint
    ub = np.ascontiguousarray(ub_calculator.getUB(), dtype=np.float64)
    u = np.ascontiguousarray(ub_calculator.getU(), dtype=np.float64)
    angle_bounds_fingerprint = _xxh3_128(bounds)
    ub_fingerprint = _xxh3_128(ub)
    u_fingerprint = _xxh3_128(u)
    wavevector = float(ub_calculator.getK())
    detector_config = (
        _jsonable(detector.get_config()) if hasattr(detector, "get_config") else None
    )
    base_context = {
        "job_sha256": job_digest,
        "ub": ub.tolist(),
        "u": u.tolist(),
        "wavevector_A^-1": wavevector,
        "detector": detector_config,
    }
    base_context_hash = sha256(
        json.dumps(base_context, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    scientific_components = {
        "algorithm": "xxh3-128-components-v1",
        "job_sha256": job_digest,
        "angle_bounds_xxh3_128": angle_bounds_fingerprint,
        "corner_rays_xxh3_128": [
            {
                "detector_tile": list(detector_tile),
                "digest": ray_fingerprints[detector_tile],
            }
            for detector_tile in detector_tiles
        ],
        "ub_xxh3_128": ub_fingerprint,
        "u_xxh3_128": u_fingerprint,
        "wavevector_A^-1": wavevector,
    }
    scientific_context_hash = sha256(
        json.dumps(
            scientific_components, sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()
    scientific_context = {
        **base_context,
        "base_context_sha256": base_context_hash,
        **scientific_components,
        "scientific_context_sha256": scientific_context_hash,
        "detector_tiles": [list(detector_tile) for detector_tile in detector_tiles],
    }
    task_seed = {
        "spec": spec.digest,
        "frame_range": [start, stop],
        "detector_tiles": [
            list(detector_tile) for detector_tile in detector_tiles
        ],
        "scientific_context": scientific_context_hash,
    }
    task_id = sha256(
        json.dumps(task_seed, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:24]
    kernels = {
        grid.grid_name: _kernel_for_grid(
            spec,
            grid,
            ub_calculator,
            threads=kernel_threads,
            memory_budget_bytes=kernel_memory_budget_bytes,
        )
        for grid in spec.grids
    }
    grid_batches: dict[str, list[Mapping[str, np.ndarray]]] = {
        grid.grid_name: [] for grid in spec.grids
    }
    accumulated_bytes = 0
    segment_index = 0
    partitions = []

    def flush():
        nonlocal accumulated_bytes, segment_index
        if not any(grid_batches.values()):
            return
        for grid in spec.grids:
            batch = _reduce_batches(grid_batches[grid.grid_name])
            if batch["chunk_id"].size == 0:
                continue
            buckets = batch["chunk_id"] // spec.partition_chunk_span
            for bucket in np.unique(buckets):
                selected = buckets == bucket
                subset = {
                    name: values[selected] for name, values in batch.items()
                }
                uri = _join_uri(
                    output_uri,
                    f"grid={grid.grid_name}",
                    f"bucket={int(bucket):08d}",
                    (
                        f"task={task_id}-"
                        f"segment={segment_index:06d}.parquet"
                    ),
                )
                _write_parquet(
                    subset,
                    uri,
                    {
                        "spec_hash": spec.digest,
                        "task_id": task_id,
                        "segment": segment_index,
                        "grid_name": grid.grid_name,
                        "bucket": int(bucket),
                    },
                )
                checksum, size_bytes = _uri_checksum_and_size(uri)
                partition = _PartitionFile(
                    grid_name=grid.grid_name,
                    bucket=int(bucket),
                    uri=uri,
                    rows=int(subset["chunk_id"].size),
                    checksum=checksum,
                    size_bytes=size_bytes,
                )
                _mark_scratch_verified(partition, verification_cache)
                partitions.append(partition)
        for batches in grid_batches.values():
            batches.clear()
        accumulated_bytes = 0
        segment_index += 1

    for offset, frame_index in enumerate(range(start, stop)):
        image_payload = (
            scan.get_raw_img(frame_index)
            if image_payloads is None
            else image_payloads[frame_index]
        )
        image = np.asarray(image_payload.img)
        frame_correction = getattr(
            correction_pipeline, "correct_frame", None
        )
        corrected_frame = (
            frame_correction(image_payload, image, frame_index)
            if callable(frame_correction)
            else None
        )
        for detector_tile in detector_tiles:
            row_start, row_stop, column_start, column_stop = detector_tile
            selection = np.s_[
                row_start:row_stop, column_start:column_stop
            ]
            if corrected_frame is None:
                intensity, variance, mask = correction_pipeline(
                    image_payload,
                    image[selection],
                    frame_index,
                    detector_tile,
                )
            else:
                intensity, variance, mask = (
                    values[selection] for values in corrected_frame
                )
            intensity = np.ascontiguousarray(intensity, dtype=np.float64)
            variance = np.ascontiguousarray(variance, dtype=np.float64)
            mask = np.ascontiguousarray(mask, dtype=bool)
            for grid in spec.grids:
                batch = kernels[grid.grid_name].accumulate(
                    intensity,
                    variance,
                    mask,
                    ray_arrays[detector_tile],
                    np.ascontiguousarray(bounds[offset, 0]),
                    np.ascontiguousarray(bounds[offset, 1]),
                )
                batch_bytes = sum(
                    np.asarray(values).nbytes for values in batch.values()
                )
                if (
                    accumulated_bytes
                    and accumulated_bytes + batch_bytes
                    > accumulation_budget_bytes
                ):
                    flush()
                grid_batches[grid.grid_name].append(batch)
                accumulated_bytes += batch_bytes
                if accumulated_bytes >= accumulation_budget_bytes:
                    flush()
        if image_progress is not None:
            image_progress(
                frame_index,
                accumulated_bytes,
                segment_index,
            )
    flush()
    scientific_context["accumulation_segments"] = segment_index
    return _TaskManifest(
        kind="map",
        task_id=task_id,
        spec_hash=spec.digest,
        status="complete",
        spec=spec.to_dict(),
        frame_range=(start, stop),
        detector_tile=None,
        partitions=partitions,
        metadata=scientific_context,
    )


def _map_partition(
    spec: _ReconstructionSpec,
    scan,
    detector,
    ub_calculator,
    frame_range,
    detector_tile,
    angle_bounds_rad,
    output_uri,
    *,
    correction_pipeline: Callable,
    job_digest: str,
    image_payloads: Mapping[int, object] | None = None,
    corner_rays: np.ndarray | None = None,
    corner_rays_fingerprint: str | None = None,
    verification_cache: set[tuple[str, str]] | None = None,
    kernel_threads: int | None = None,
    kernel_memory_budget_bytes: int | None = None,
    accumulation_budget_bytes: int = 64 * 1024 * 1024,
    image_progress: Callable[[int, int, int], None] | None = None,
) -> _TaskManifest:
    """Map one detector tile through the streaming range implementation."""
    tile = tuple(map(int, detector_tile))
    result = _map_frame_range(
        spec,
        scan,
        detector,
        ub_calculator,
        frame_range,
        (tile,),
        angle_bounds_rad,
        output_uri,
        correction_pipeline=correction_pipeline,
        job_digest=job_digest,
        image_payloads=image_payloads,
        corner_rays=None if corner_rays is None else {tile: corner_rays},
        corner_rays_fingerprints=(
            None
            if corner_rays_fingerprint is None
            else {tile: corner_rays_fingerprint}
        ),
        verification_cache=verification_cache,
        kernel_threads=kernel_threads,
        kernel_memory_budget_bytes=kernel_memory_budget_bytes,
        accumulation_budget_bytes=accumulation_budget_bytes,
        image_progress=image_progress,
    )
    result.detector_tile = tile
    return result


def _reduce_partition(
    manifest_set,
    output_uri,
    *,
    verification_cache: set[tuple[str, str]] | None = None,
    memory_budget_bytes: int = 1024**3,
    workers: int = 1,
    progress: Callable[[int, int, str], None] | None = None,
    checkpoint_root=None,
) -> _TaskManifest:
    """Externally reduce mapping partitions into bounded chunk shards.

    Contiguous shard ranges run concurrently with private forward-only
    Parquet readers. The total memory budget is divided between workers.
    """
    memory_budget_bytes = max(1, int(memory_budget_bytes))
    requested_workers = max(1, int(workers))
    worker_limit = max(
        1,
        memory_budget_bytes // _MIN_REDUCER_WORKER_MEMORY,
    )
    configured_workers = min(requested_workers, worker_limit)
    worker_memory_budget = max(
        1,
        memory_budget_bytes // configured_workers,
    )
    manifests = [_read_manifest(value) for value in manifest_set]
    if not manifests:
        raise ValueError("At least one mapping manifest is required")
    if any(manifest.kind != "map" for manifest in manifests):
        raise ValueError("reduce_partition accepts mapping manifests only")
    spec_hashes = {manifest.spec_hash for manifest in manifests}
    if len(spec_hashes) != 1:
        raise ValueError("All mapping manifests must use the same specification")
    spec = _ReconstructionSpec.from_dict(manifests[0].spec)
    base_contexts = {
        manifest.metadata.get("base_context_sha256") for manifest in manifests
    }
    if len(base_contexts) != 1:
        raise ValueError(
            "Mapping manifests contain different geometry, inputs, or correction "
            "identities"
        )
    grouped: dict[tuple[str, int], list[_PartitionFile]] = {}
    for manifest in manifests:
        if manifest.status != "complete":
            raise ValueError(f"Mapping task {manifest.task_id} is not complete")
        for partition in manifest.partitions:
            grouped.setdefault((partition.grid_name, partition.bucket), []).append(
                partition
            )
    source_tasks = sorted(manifest.task_id for manifest in manifests)
    reduce_id = sha256(
        (spec.digest + "\0" + "\0".join(source_tasks)).encode("utf-8")
    ).hexdigest()[:24]
    base_metadata = {
        **{
            key: value
            for key, value in manifests[0].metadata.items()
            if key
            not in {
                "angle_bounds_sha256",
                "corner_rays_sha256",
                "angle_bounds_xxh3_128",
                "corner_rays_xxh3_128",
                "ub_xxh3_128",
                "u_xxh3_128",
                "scientific_context_sha256",
            }
        },
        "source_scientific_contexts": {
            manifest.task_id: manifest.metadata.get(
                "scientific_context_sha256"
            )
            for manifest in manifests
        },
    }
    grid_by_name = {grid.grid_name: grid for grid in spec.grids}
    local_span = max(
        256 * 1024,
        min(
            4 * 1024 * 1024,
            worker_memory_budget // (8 * 48),
        ),
    )
    plans = []
    for (grid_name, bucket), partitions in sorted(grouped.items()):
        grid = grid_by_name[grid_name]
        chunk_grid = tuple(
            math.ceil(size / chunk)
            for size, chunk in zip(grid.shape, grid.chunk_shape)
        )
        total_chunks = math.prod(chunk_grid)
        chunk_start = bucket * spec.partition_chunk_span
        chunk_stop = min(
            total_chunks,
            (bucket + 1) * spec.partition_chunk_span,
        )
        for chunk_id in range(chunk_start, chunk_stop):
            coordinates = _chunk_coordinates(chunk_id, grid)
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
            for shard_start in range(0, local_stop, local_span):
                plans.append(
                    (
                        grid_name,
                        bucket,
                        chunk_id,
                        shard_start,
                        min(shard_start + local_span, local_stop),
                        tuple(sorted(partitions, key=lambda item: item.uri)),
                    )
                )

    checkpoint_path = (
        None
        if checkpoint_root is None
        else Path(checkpoint_root) / f"{reduce_id}.json"
    )
    chunks = []
    completed_shards = set()
    if checkpoint_path is not None and checkpoint_path.exists():
        try:
            checkpoint = _read_manifest(checkpoint_path)
            if (
                checkpoint.kind != "reduce"
                or checkpoint.task_id != reduce_id
                or checkpoint.spec_hash != spec.digest
                or checkpoint.source_tasks != source_tasks
            ):
                raise ValueError("Reduction checkpoint identity mismatch")
            for chunk in checkpoint.chunks:
                if _verify_scratch_file(chunk, verification_cache):
                    chunks.append(chunk)
            completed_shards.update(
                checkpoint.metadata.get("completed_shards", ())
            )
        except (OSError, ValueError, RuntimeError, TypeError):
            chunks = []
            completed_shards.clear()

    def shard_key(grid_name, chunk_id, shard_start, shard_stop):
        return (
            f"{grid_name}:{chunk_id}:"
            f"{shard_start}:{shard_stop}"
        )

    def checkpoint(status):
        ordered_chunks = sorted(
            chunks,
            key=lambda item: (
                item.grid_name,
                item.chunk_id,
                item.shard_start,
                item.shard_stop,
                item.uri,
            ),
        )
        manifest = _TaskManifest(
            kind="reduce",
            task_id=reduce_id,
            spec_hash=spec.digest,
            status=status,
            spec=spec.to_dict(),
            chunks=ordered_chunks,
            source_tasks=source_tasks,
            metadata={
                **base_metadata,
                "local_voxel_shard_span": local_span,
                "reducer_worker_capacity": configured_workers,
                "reducer_workers_used": actual_workers,
                "reducer_memory_bytes_per_worker": worker_memory_budget,
                "completed_shards": sorted(completed_shards),
            },
        )
        if checkpoint_path is not None:
            _write_manifest(manifest, checkpoint_path)
        return manifest

    existing = {
        shard_key(
            chunk.grid_name,
            chunk.chunk_id,
            chunk.shard_start,
            chunk.shard_stop,
        )
        for chunk in chunks
    }
    checkpoint_interval = 64
    total_plans = len(plans)
    completed = 0
    active_plans = []
    for plan in plans:
        grid_name, _, chunk_id, shard_start, shard_stop, _ = plan
        key = shard_key(
            grid_name,
            chunk_id,
            shard_start,
            shard_stop,
        )
        if key in existing or key in completed_shards:
            completed += 1
            if progress is not None:
                progress(
                    completed,
                    total_plans,
                    f"Reusing reduced shard {completed}/{total_plans}",
                )
            continue
        active_plans.append(plan)

    actual_workers = min(configured_workers, len(active_plans))
    worker_memory_budget = max(
        1,
        memory_budget_bytes // max(1, actual_workers),
    )

    active_partitions = {
        (partition.uri, partition.checksum): partition
        for plan in active_plans
        for partition in plan[5]
    }
    partitions_to_verify = [
        partition
        for key, partition in active_partitions.items()
        if verification_cache is None or key not in verification_cache
    ]
    if partitions_to_verify:
        verification_workers = min(
            actual_workers,
            len(partitions_to_verify),
        )
        with ThreadPoolExecutor(
            max_workers=verification_workers
        ) as executor:
            futures = {
                executor.submit(
                    _verify_scratch_file,
                    partition,
                    None,
                ): partition
                for partition in partitions_to_verify
            }
            for verified, future in enumerate(
                as_completed(futures),
                start=1,
            ):
                partition = futures[future]
                future.result()
                _mark_scratch_verified(partition, verification_cache)
                if progress is not None:
                    progress(
                        completed,
                        total_plans,
                        (
                            "Verified mapping partition "
                            f"{verified}/{len(partitions_to_verify)} "
                            f"with {verification_workers} workers"
                        ),
                    )

    result_queue = Queue()

    def reduce_plan_batch(plan_batch):
        readers = {}
        active_group = None
        try:
            for (
                grid_name,
                bucket,
                chunk_id,
                shard_start,
                shard_stop,
                partitions,
            ) in plan_batch:
                group_key = (grid_name, bucket)
                if group_key != active_group:
                    readers.clear()
                    active_group = group_key
                for partition in partitions:
                    if partition.uri in readers:
                        continue
                    batch_rows = max(
                        4096,
                        min(
                            131072,
                            worker_memory_budget
                            // (max(1, len(partitions)) * 48 * 4),
                        ),
                    )
                    readers[partition.uri] = _ParquetRangeReader(
                        partition.uri,
                        batch_size=batch_rows,
                        use_threads=actual_workers == 1,
                    )
                levels = []
                for partition in partitions:
                    batch = readers[partition.uri].read(
                        chunk_id,
                        shard_start,
                        shard_stop,
                    )
                    if not batch["chunk_id"].size:
                        continue
                    level = 0
                    while level < len(levels) and levels[level] is not None:
                        batch = _merge_sorted_batches(levels[level], batch)
                        levels[level] = None
                        level += 1
                    if level == len(levels):
                        levels.append(batch)
                    else:
                        levels[level] = batch
                reduced = _empty_batch()
                for batch in reversed(levels):
                    if batch is not None:
                        reduced = _merge_sorted_batches(reduced, batch)
                chunk = None
                if reduced["chunk_id"].size:
                    uri = _join_uri(
                        output_uri,
                        f"grid={grid_name}",
                        (
                            f"chunk={int(chunk_id):016d}-"
                            f"shard={shard_start:016d}-{shard_stop:016d}.parquet"
                        ),
                    )
                    _write_parquet(
                        reduced,
                        uri,
                        {
                            "spec_hash": spec.digest,
                            "reduce_id": reduce_id,
                            "grid_name": grid_name,
                            "chunk_id": int(chunk_id),
                            "shard_start": shard_start,
                            "shard_stop": shard_stop,
                            "source_bucket": bucket,
                        },
                    )
                    checksum, size_bytes = _uri_checksum_and_size(uri)
                    chunk = _ChunkFile(
                        grid_name=grid_name,
                        chunk_id=int(chunk_id),
                        shard_start=shard_start,
                        shard_stop=shard_stop,
                        uri=uri,
                        rows=int(reduced["chunk_id"].size),
                        checksum=checksum,
                        size_bytes=size_bytes,
                    )
                result_queue.put(
                    (
                        shard_key(
                            grid_name,
                            chunk_id,
                            shard_start,
                            shard_stop,
                        ),
                        chunk,
                        grid_name,
                        chunk_id,
                        shard_start,
                        shard_stop,
                    )
                )
        finally:
            result_queue.put(None)

    if active_plans:
        plans_per_worker, extra = divmod(
            len(active_plans),
            actual_workers,
        )
        plan_batches = []
        start = 0
        for worker_index in range(actual_workers):
            stop = start + plans_per_worker + (worker_index < extra)
            plan_batches.append(active_plans[start:stop])
            start = stop

        with ThreadPoolExecutor(max_workers=actual_workers) as executor:
            futures = [
                executor.submit(reduce_plan_batch, batch)
                for batch in plan_batches
            ]
            finished_workers = 0
            while finished_workers < len(futures):
                outcome = result_queue.get()
                if outcome is None:
                    finished_workers += 1
                    continue
                (
                    key,
                    chunk,
                    grid_name,
                    chunk_id,
                    shard_start,
                    shard_stop,
                ) = outcome
                if chunk is not None:
                    _mark_scratch_verified(chunk, verification_cache)
                    chunks.append(chunk)
                completed_shards.add(key)
                completed += 1
                if (
                    checkpoint_path is not None
                    and (
                        completed % checkpoint_interval == 0
                        or completed == total_plans
                    )
                ):
                    checkpoint("running")
                if progress is not None:
                    progress(
                        completed,
                        total_plans,
                        (
                            f"Reduced shard {completed}/{total_plans} "
                            f"with {actual_workers} workers: "
                            f"{grid_name} chunk {chunk_id}, "
                            f"{shard_start}:{shard_stop}"
                        ),
                    )
            for future in futures:
                future.result()
    return checkpoint("complete")


def _compression_kwargs(name):
    if name.startswith("database:"):
        from ...app.database import FILTERS

        filter_name = name.partition(":")[2]
        if filter_name not in FILTERS:
            raise ValueError(f"Unknown orGUI database compression: {filter_name}")
        selected = FILTERS[filter_name]
        if selected is None:
            return {}
        if isinstance(selected, Mapping):
            return dict(selected)
        return {"compression": selected}
    normalized = name.lower()
    if normalized in {"none", "raw"}:
        return {}
    if normalized == "bitshuffle-lz4":
        try:
            import hdf5plugin

            return {"compression": hdf5plugin.Bitshuffle(cname="lz4")}
        except ImportError:
            return {"compression": "gzip", "compression_opts": 4, "shuffle": True}
    if normalized == "lzf":
        return {"compression": "lzf", "shuffle": True}
    if normalized == "gzip":
        return {"compression": "gzip", "compression_opts": 4, "shuffle": True}
    raise ValueError(f"Unsupported HDF5 compression: {name}")


def _chunk_coordinates(chunk_id, grid):
    chunk_grid = tuple(
        math.ceil(size / chunk) for size, chunk in zip(grid.shape, grid.chunk_shape)
    )
    chunk_x, remainder = divmod(chunk_id, chunk_grid[1] * chunk_grid[2])
    chunk_y, chunk_z = divmod(remainder, chunk_grid[2])
    return chunk_x, chunk_y, chunk_z


def _write_chunk(
    group,
    grid,
    chunk_files,
    verification_cache=None,
):
    chunk_files = sorted(chunk_files, key=lambda item: item.shard_start)
    if not chunk_files:
        return
    chunk_id = chunk_files[0].chunk_id
    if any(chunk.chunk_id != chunk_id for chunk in chunk_files):
        raise ValueError("Reduced shard group contains multiple chunks")
    coordinates = _chunk_coordinates(chunk_id, grid)
    starts = tuple(
        coordinate * chunk
        for coordinate, chunk in zip(coordinates, grid.chunk_shape)
    )
    stops = tuple(
        min(start + chunk, size)
        for start, chunk, size in zip(starts, grid.chunk_shape, grid.shape)
    )
    local_shape = tuple(stop - start for start, stop in zip(starts, stops))
    weighted_intensity = np.zeros(local_shape, dtype=np.float64)
    weighted_variance = np.zeros(local_shape, dtype=np.float64)
    weight = np.zeros(local_shape, dtype=np.float64)
    contributors = np.zeros(local_shape, dtype=np.uint64)
    for chunk_file in chunk_files:
        _verify_scratch_file(chunk_file, verification_cache)
        batch = _read_parquet(chunk_file.uri)
        if batch["chunk_id"].size and np.any(
            batch["chunk_id"] != np.uint64(chunk_id)
        ):
            raise ValueError(f"{chunk_file.uri} contains more than one chunk")
        local = batch["local_voxel_id"].astype(np.uint64, copy=False)
        if local.size and (
            np.any(local < np.uint64(chunk_file.shard_start))
            or np.any(local >= np.uint64(chunk_file.shard_stop))
        ):
            raise ValueError(
                f"{chunk_file.uri} contains records outside its shard"
            )
        local_x, remainder = np.divmod(
            local, np.uint64(grid.chunk_shape[1] * grid.chunk_shape[2])
        )
        local_y, local_z = np.divmod(
            remainder,
            np.uint64(grid.chunk_shape[2]),
        )
        valid = (
            (local_x < local_shape[0])
            & (local_y < local_shape[1])
            & (local_z < local_shape[2])
        )
        if not np.all(valid):
            raise ValueError(
                f"{chunk_file.uri} contains an invalid local voxel ID"
            )
        index = (
            local_x.astype(int),
            local_y.astype(int),
            local_z.astype(int),
        )
        weighted_intensity[index] = batch["weighted_intensity"]
        weighted_variance[index] = batch["weighted_variance"]
        weight[index] = batch["weight"]
        contributors[index] = batch["contributors"]
    populated = weight > 0
    weighted_intensity[~populated] = np.nan
    weighted_variance[~populated] = np.nan
    np.divide(
        weighted_intensity,
        weight,
        out=weighted_intensity,
        where=populated,
    )
    np.divide(
        weighted_variance,
        weight,
        out=weighted_variance,
        where=populated,
    )
    np.divide(
        weighted_variance,
        weight,
        out=weighted_variance,
        where=populated,
    )
    selection = tuple(slice(start, stop) for start, stop in zip(starts, stops))
    group["intensity"][selection] = weighted_intensity
    group["variance"][selection] = weighted_variance
    group["weight"][selection] = weight
    group["contributors"][selection] = contributors


def _finalize_reconstruction(
    manifest_set,
    output_path,
    *,
    provenance: Mapping[str, Any] | None = None,
    config=None,
    chunk_progress: Callable | None = None,
    verification_cache: set[tuple[str, str]] | None = None,
):
    """Create a conventional HDF5 reconstruction from reduced chunk files."""
    manifests = [_read_manifest(value) for value in manifest_set]
    if not manifests:
        raise ValueError("At least one reduction manifest is required")
    if any(manifest.kind != "reduce" for manifest in manifests):
        raise ValueError("finalize_reconstruction accepts reduction manifests only")
    spec_hashes = {manifest.spec_hash for manifest in manifests}
    if len(spec_hashes) != 1:
        raise ValueError("Reduction manifests use different specifications")
    spec = _ReconstructionSpec.from_dict(manifests[0].spec)
    grid_by_name = {grid.grid_name: grid for grid in spec.grids}
    chunk_files = {}
    for manifest in manifests:
        if manifest.status != "complete":
            raise ValueError(f"Reduction task {manifest.task_id} is not complete")
        for chunk in manifest.chunks:
            key = (chunk.grid_name, chunk.chunk_id)
            chunk_files.setdefault(key, []).append(chunk)
    for key, shards in chunk_files.items():
        ordered = sorted(shards, key=lambda item: item.shard_start)
        for previous, current in zip(ordered, ordered[1:]):
            if current.shard_start < previous.shard_stop:
                raise ValueError(f"Overlapping reduced shards for chunk {key}")
    output_path = Path(output_path).absolute()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(output_path.name + ".tmp")
    compression = _compression_kwargs(spec.compression)
    with h5py.File(temporary_path, "w") as h5file:
        entry = h5file.create_group("entry")
        entry.attrs["NX_class"] = "NXentry"
        process = entry.create_group("reconstruction")
        process.attrs["NX_class"] = "NXprocess"
        process.attrs["program"] = "orGUI"
        process.attrs["spec_sha256"] = spec.digest
        process.create_dataset(
            "configuration_json",
            data=json.dumps(spec.to_dict(), sort_keys=True),
        )
        if provenance:
            process.create_dataset(
                "provenance_json",
                data=json.dumps(dict(provenance), sort_keys=True, default=str),
            )
        process.create_dataset(
            "scientific_context_json",
            data=json.dumps(manifests[0].metadata, sort_keys=True, default=str),
        )
        process.create_dataset(
            "marginal_variance_warning",
            data=(
                "Adaptive footprint splitting and pixel repair can introduce "
                "cross-voxel or spatial covariance. Variance datasets contain "
                "marginal variances only."
            ),
        )
        if config is not None:
            from ...app.config_data import ConfigHandler

            ConfigHandler(create_dataset_args=compression).write_scan_config(
                entry, config, source="reciprocal_space_reconstruction"
            )
        results = process.create_group("results")
        groups = {}
        for grid in spec.grids:
            group = results.create_group(grid.grid_name)
            group.attrs["NX_class"] = "NXdata"
            group.attrs["signal"] = "intensity"
            group.attrs["axes"] = np.asarray(
                ["h", "k", "l"] if grid.frame == "hkl" else ["qx", "qy", "qz"],
                dtype=h5py.string_dtype(),
            )
            group.attrs["coordinate_frame"] = grid.frame
            group.attrs["units"] = "r.l.u." if grid.frame == "hkl" else "Angstrom^-1"
            axis_names = ("h", "k", "l") if grid.frame == "hkl" else ("qx", "qy", "qz")
            for axis, axis_name in enumerate(axis_names):
                values = grid.minimum[axis] + (
                    np.arange(grid.shape[axis], dtype=np.float64) + 0.5
                ) * grid.step[axis]
                dataset = group.create_dataset(axis_name, data=values)
                dataset.attrs["units"] = (
                    "r.l.u." if grid.frame == "hkl" else "Angstrom^-1"
                )
                dataset.attrs["long_name"] = f"{axis_name} voxel center"
            hdf5_chunks = tuple(
                min(size, chunk)
                for size, chunk in zip(grid.shape, grid.chunk_shape)
            )
            group.create_dataset(
                "intensity",
                shape=grid.shape,
                dtype=np.float64,
                chunks=hdf5_chunks,
                fillvalue=np.nan,
                **compression,
            )
            group.create_dataset(
                "variance",
                shape=grid.shape,
                dtype=np.float64,
                chunks=hdf5_chunks,
                fillvalue=np.nan,
                **compression,
            )
            group.create_dataset(
                "weight",
                shape=grid.shape,
                dtype=np.float64,
                chunks=hdf5_chunks,
                fillvalue=0.0,
                **compression,
            )
            group.create_dataset(
                "contributors",
                shape=grid.shape,
                dtype=np.uint64,
                chunks=hdf5_chunks,
                fillvalue=0,
                **compression,
            )
            groups[grid.grid_name] = group
        for written, ((grid_name, _), chunk_shards) in enumerate(
            sorted(chunk_files.items()), start=1
        ):
            if grid_name not in grid_by_name:
                raise ValueError(f"Unknown grid in chunk manifest: {grid_name}")
            _write_chunk(
                groups[grid_name],
                grid_by_name[grid_name],
                chunk_shards,
                verification_cache,
            )
            if chunk_progress is not None:
                chunk_progress(written, len(chunk_files))
        h5file.flush()
    temporary_path.replace(output_path)
    with h5py.File(output_path, "r") as h5file:
        results = h5file["entry/reconstruction/results"]
        for grid in spec.grids:
            group = results[grid.grid_name]
            for dataset in ("intensity", "variance", "weight", "contributors"):
                if group[dataset].shape != grid.shape:
                    raise OSError(
                        f"Final reconstruction dataset has invalid shape: "
                        f"{grid.grid_name}/{dataset}"
                    )
    digest = sha256()
    with output_path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    result = {
        "path": str(output_path),
        "sha256": digest.hexdigest(),
        "bytes": output_path.stat().st_size,
        "chunks_written": len(chunk_files),
    }
    return result


__all__ = []

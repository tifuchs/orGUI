"""Central orGUI job model for reciprocal-space reconstruction."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from functools import lru_cache
from hashlib import sha256
import importlib
import json
import math
from pathlib import Path
from queue import Empty, SimpleQueue
import shutil
import threading
import time
from typing import Any
import warnings

import h5py
import numpy as np

from .app.config_data import ConfigData
from .app.database import FILTERS, config_data_from_json, config_data_to_json
from .app.mask_config import create_pixel_repair_plan
from .backend.scans import ScanReference
from .datautils.xrayutils.reconstruction import (
    _CHECKPOINT_BYTES_PER_ROW,
    _CheckpointRouter,
    _GridSpec,
    _ReconstructionSpec,
    _build_kernels,
    _calibration_probe_all_grids,
    _detector_corner_rays,
    _discover_checkpoint_state,
    _files_per_job,
    _finalize_reconstruction,
    _kernel_for_grid,
    _kernel_threads_sweep,
    _map_frame_group,
    _tile_ray_arrays,
    _validate_mapping_setup,
)


JOB_SCHEMA_VERSION = 6
ACCURACY_DEPTHS = {
    "center": 0,
    "low": 1,
    "balanced": 2,
    "high": 3,
    "very_high": 4,
    "maximum": 5,
}
WORK_BLOCK_PRESETS = {
    "minimum": 1024,
    "tiny": 2048,
    "small": 4096,
    "low": 8192,
    "medium": 16384,
    "high": 32768,
    "maximum": 65536,
}
"""Named native work-block sizes, in detector pixels at ``center`` accuracy.

A block's working set is the pixel stream it reads -- about 41 bytes per
pixel of intensity, variance, mask and detector corner rays -- plus the
accumulator it builds, roughly 72 bytes per distinct voxel reached. So a
name selects a *cache scale* rather than a raw count: ``minimum`` is about
50-75 KiB and lands in a typical per-core L1 data cache, ``medium`` is
about 0.8-1.2 MiB and lands in a typical 1 MiB L2. ``medium`` measured
fastest across grid resolutions and thread counts and is the default.

The pixel count halves with each adaptive depth, because deeper
subdivision adds per-pixel state (the dyadic coordinate cache alone grows
as ``(2**(depth+1)+1)**3``) that competes for the same cache, so holding
the name fixed holds the working set roughly fixed.
"""

DEFAULT_WORK_BLOCK = "medium"

_MAP_NODE_BYTES = 72
"""Bytes the native accumulator spends per distinct voxel: an 8-byte key,
a 32-byte accumulator, and 32 bytes of tree links."""

_ARENA_BUDGET_SHARE = 4
"""The per-worker arena may claim at most this fraction's reciprocal of one
thread's share of the memory budget."""

_RESERVED_RECORDS_PER_PIXEL = 4
"""Records the native arena reserves per block pixel, mirroring the kernel.
Measured density stays between 0.46 and 0.87 records per pixel at every
adaptive depth, so this leaves several times the headroom actually needed."""

AUTO_MAX_FRAMES_PER_TASK = 64
ACCUMULATION_TRANSIENT_FACTOR = 3
AUTO_MAX_ACCUMULATION_BYTES = 2 * 1024**3
MAX_CONCURRENT_ACTIVE_CHECKPOINTS = 2
"""Design doc Sec10: how many checkpoints (per grid) may have live,
non-empty accumulator state at once -- small and fixed, not user-facing."""
_RAW_IMAGE_BYTES_PER_PIXEL = 8
"""Conservative per-pixel size for one decoded detector frame held in
memory -- whether still sitting in the prefetch queue or already handed
to a compute worker. Assumes float64 (the widest dtype involved) as a
safe upper bound regardless of a scan backend's actual raw dtype."""
_CORRECTED_ARRAY_BYTES_PER_PIXEL = 8 + 8 + 1
"""_correction_pipeline.correct_frame's per-pixel output: float64
intensity, float64 variance, bool mask."""
_FRAME_RECORD_BYTES_PER_PIXEL = 2 * _CHECKPOINT_BYTES_PER_ROW
"""A frame's own mapped records, held Python-side while it is in flight:
every detector tile's output batch, plus the transient copy
_reduce_batches makes while merging those batches into one. Sized at one
record per detector pixel, which bounds every density measured across
adaptive depths (0.46 to 0.87 records per pixel), and doubled for the
merge. Without this term the worker pool was sized as though a frame cost
only its image and correction buffers, and measured peaks ran about a
third above what the pool believed it had claimed."""

_PYTHON_CORRECTION_BYTES_PER_PIXEL = (
    _RAW_IMAGE_BYTES_PER_PIXEL + _CORRECTED_ARRAY_BYTES_PER_PIXEL
)
"""Sec7's Python-side frame footprint: correct_frame runs once per
frame (not per detector tile) and holds the raw image plus every
corrected array, all full-detector-sized, simultaneously. This is real
per-in-flight-frame resident memory the native kernel's own per-tile
working-set estimate (_frame_parallelism) does not see at all."""


def _build_metadata():
    from . import __version__

    try:
        from ._build_config import BUILD_CONFIG
    except ImportError:
        build_config = {"build_config": "unavailable-from-source-tree"}
    else:
        build_config = dict(BUILD_CONFIG)
    return {
        "orgui_version": __version__,
        "numpy_version": np.__version__,
        "h5py_version": h5py.__version__,
        "native": build_config,
    }


def _sha256_file(path):
    digest = sha256()
    with Path(path).open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path, values):
    path = Path(path).absolute()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(values, indent=2, sort_keys=True), encoding="utf-8"
    )
    temporary.replace(path)


def _compression_name(value):
    for name, available in FILTERS.items():
        try:
            if value is available or value == available:
                return name
        except Exception:
            continue
    raise ValueError("The active database compression is not registered")


@dataclass(frozen=True)
class ReconstructionGrid:
    """One user-selected reciprocal-space output grid."""

    minimum: tuple[float, float, float]
    maximum: tuple[float, float, float]
    step: tuple[float, float, float]
    frame: str = "hkl"
    name: str | None = None
    chunk_shape: tuple[int, int, int] = (64, 64, 64)

    def to_spec(self):
        """Return the internal native/storage grid representation."""
        return _GridSpec(**asdict(self))


@dataclass
class ReconstructionJob:
    """Immutable experiment snapshot plus mutable deterministic task state."""

    config: dict[str, Any]
    scan_reference: dict[str, Any]
    grids: list[dict[str, Any]]
    scratch_path: str
    output_path: str
    compression: str
    assets_path: str
    assets_sha256: str
    source_fingerprint_sha256: str
    threads_per_image: int | None
    accumulation_budget_bytes: int | None
    checkpoint_count: int = 10
    checkpoint_plan: dict[str, list[list[int]]] = field(default_factory=dict)
    build_metadata: dict[str, Any] = field(default_factory=dict)
    runtime_threads: int = 1
    runtime_memory_bytes: int = 6_000 * 1024 * 1024
    angle_fallback: str = "stationary"
    accuracy: str = "balanced"
    advanced_depth: int | None = None
    thread_override: int | None = None
    memory_override_bytes: int | None = None
    frame_batch: int | None = None
    tile_shape: tuple[int, int] | None = None
    #: ``None`` for the default preset, a :data:`WORK_BLOCK_PRESETS` name,
    #: or an explicit pixel count. See :func:`resolve_work_block_pixels`.
    work_block_pixels: int | str | None = None
    user_note: str = ""
    status: str = "prepared"
    output_sha256: str | None = None
    correction_provenance: dict[str, Any] = field(default_factory=dict)
    cleanup_errors: list[str] = field(default_factory=list)
    cluster_settings: dict[str, Any] = field(default_factory=dict)
    schema_version: int = JOB_SCHEMA_VERSION

    def to_dict(self):
        """Return a JSON-compatible job descriptor."""
        return asdict(self)

    @classmethod
    def from_dict(cls, values):
        """Build a reconstruction job from JSON-compatible values."""
        values = dict(values)
        if values.get("schema_version") != JOB_SCHEMA_VERSION:
            raise ValueError("Unsupported reconstruction job schema")
        if values.get("tile_shape") is not None:
            values["tile_shape"] = tuple(values["tile_shape"])
        return cls(**values)

    @property
    def digest(self):
        """Return the immutable scientific job digest.

        Also used as every checkpoint file's resume-identity fingerprint
        (design doc Sec11) -- stronger than the bare
        ``_ReconstructionSpec.digest``, since it also covers scan/config
        identity, so two different jobs that happen to share identical
        grid/depth/thread settings can never be mistaken for each other's
        checkpoints even if they were prepared into the same scratch path.
        """
        values = self.to_dict()
        for key in (
            "status",
            "checkpoint_plan",
            "output_sha256",
            "correction_provenance",
            "cleanup_errors",
        ):
            values.pop(key, None)
        encoded = json.dumps(
            values, sort_keys=True, separators=(",", ":")
        ).encode()
        return sha256(encoded).hexdigest()

    @property
    def config_data(self):
        """Return the central configuration snapshot."""
        return config_data_from_json(self.config)

    @property
    def scan(self):
        """Reopen and verify the referenced scan."""
        return ScanReference.from_dict(self.scan_reference).open()

    def internal_spec(self):
        """Build the internal native/storage specification."""
        depth = (
            self.advanced_depth
            if self.advanced_depth is not None
            else ACCURACY_DEPTHS[self.accuracy]
        )
        memory = self.memory_override_bytes or self.runtime_memory_bytes
        threads = self.thread_override or self.runtime_threads
        work_block = resolve_work_block_pixels(
            self.work_block_pixels, depth, memory, threads
        )
        return _ReconstructionSpec(
            grids=tuple(_GridSpec(**grid) for grid in self.grids),
            max_depth=depth,
            threads=threads,
            work_block_pixels=work_block,
            memory_budget_bytes=memory,
            checkpoint_count=self.checkpoint_count,
            compression=f"database:{self.compression}",
            infer_angle_bounds=self.angle_fallback == "midpoint",
        )


def read_job(path):
    """Read a reconstruction job JSON file."""
    return ReconstructionJob.from_dict(
        json.loads(Path(path).read_text(encoding="utf-8"))
    )


def write_job(job, path):
    """Atomically write a reconstruction job JSON file."""
    _atomic_json(path, job.to_dict())
    return str(Path(path).absolute())


def _snapshot_assets(gui, config, assets_path):
    assets_path = Path(assets_path).absolute()
    assets_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = assets_path.with_name(assets_path.name + ".tmp")
    correction = config.corrections
    shape = tuple(gui.ubcalc.detectorCal.detector.shape)
    mask = gui.get_detector_mask(shape) if correction.use_mask else None
    if correction.use_mask and mask is None:
        raise ValueError(
            "Pixel-mask correction is enabled, but no active mask matches "
            f"the detector image shape {shape}. Load a matching mask or "
            "disable pixel-mask correction before preparing the job."
        )
    background = (
        getattr(gui, "background_image", None)
        if correction.use_background
        else None
    )
    background_variance = getattr(gui, "background_variance", None)
    with h5py.File(temporary, "w") as h5file:
        h5file.attrs["NX_class"] = "NXroot"
        h5file.attrs["orgui_job_assets"] = JOB_SCHEMA_VERSION
        if mask is not None:
            h5file.create_dataset("mask", data=np.asarray(mask, dtype=np.bool_))
            correction.mask_asset = "/mask"
        if background is not None:
            h5file.create_dataset(
                "background", data=np.asarray(background, dtype=np.float64)
            )
            correction.background_asset = "/background"
        if background_variance is not None:
            h5file.create_dataset(
                "background_variance",
                data=np.asarray(background_variance, dtype=np.float64),
            )
            correction.background_variance_asset = "/background_variance"
        if correction.repair_masked_pixels:
            enabled, settings, rows, columns = gui._repair_config_for_image(shape)
            correction.repair_masked_pixels = bool(enabled)
            if enabled:
                repair = h5file.create_group("repair")
                repair.attrs["max_component_pixels"] = settings.max_component_pixels
                repair.attrs["max_span"] = settings.max_span
                repair.attrs["radius"] = settings.radius
                repair.attrs["min_valid_neighbors"] = settings.min_valid_neighbors
                repair.create_dataset("row_gaps", data=rows)
                repair.create_dataset("column_gaps", data=columns)
    temporary.replace(assets_path)
    return _sha256_file(assets_path)


def derive_grid(config, scan, *, frame="hkl", name=None):
    """Estimate grid coverage and spacing from current experiment geometry."""
    detector = config.detector
    rows, columns = detector.detector.shape
    detector_rows = np.array([-0.5, -0.5, rows - 0.5, rows - 0.5])
    detector_columns = np.array(
        [-0.5, columns - 0.5, -0.5, columns - 0.5]
    )
    gamma, delta = detector.primBeamPoints(
        detector_rows, detector_columns
    )
    cosine_gamma = np.cos(gamma)
    edge_rays = np.column_stack(
        (
            np.sin(delta) * cosine_gamma,
            np.cos(delta) * cosine_gamma,
            np.sin(gamma),
        )
    )
    corner_rays = np.ascontiguousarray(edge_rays.reshape(2, 2, 3))
    bounds = scan.exposure_angle_bounds(config, fallback="stationary")
    from .datautils.xrayutils.reconstruction import _native_module

    ub = config.ub_calculator
    kernel = _native_module().ReconstructionKernel(
        np.full(3, -1e12),
        np.ones(3),
        # Only kernel.coordinate() is called here, which never looks at the
        # grid shape; it just has to be large enough to stand in for
        # "unbounded" and small enough for the kernel's packed voxel
        # identifier (64 bits across the three axes).
        np.full(3, 2_000_000, dtype=np.int64),
        np.ones(3, dtype=np.int64),
        frame,
        float(ub.getK()),
        np.ascontiguousarray(np.linalg.inv(ub.getUB())),
        np.ascontiguousarray(np.linalg.inv(ub.getU())),
        0,
        1,
        1,
        1024 * 1024,
    )
    coordinates = []
    for angles_start, angles_end in bounds:
        angles_start = np.ascontiguousarray(angles_start, dtype=np.float64)
        angles_end = np.ascontiguousarray(angles_end, dtype=np.float64)
        for u in (0.0, 1.0):
            for v in (0.0, 1.0):
                for t in (0.0, 1.0):
                    coordinates.append(
                        kernel.coordinate(
                            corner_rays,
                            angles_start,
                            angles_end,
                            0,
                            0,
                            u,
                            v,
                            t,
                        )
                    )
    coordinates = np.asarray(coordinates)
    minimum = np.min(coordinates, axis=0)
    maximum = np.max(coordinates, axis=0)
    extent = np.maximum(maximum - minimum, 1e-6)
    samples = max(64, min(512, max(rows, columns, len(scan))))
    step = extent / samples
    minimum -= step
    maximum += step
    return ReconstructionGrid(
        tuple(minimum),
        tuple(maximum),
        tuple(step),
        frame=frame,
        name=name,
    )


def estimate_geometry_steps(
    config,
    scan,
    *,
    frame="hkl",
    percentile=10.0,
    detector_samples=5,
    frame_samples=32,
):
    """Estimate axis steps from local detector/scan geometry Jacobians.

    Detector-pixel rows and columns are treated as unit-width uniform
    variables. The third Jacobian direction is the exposure sweep, or the
    adjacent-frame center displacement for stationary exposures. The returned
    values are the requested percentile of the local one-sigma axis
    projections in r.l.u. for ``hkl`` or ``Angstrom^-1`` for Q frames.

    :param config:
        Central :class:`ConfigData` experiment snapshot.
    :param scan:
        Active scan backend providing exposure angle bounds in radians.
    :param str frame:
        Output coordinate frame.
    :param float percentile:
        Robust percentile of sampled local axis resolutions.
    :param int detector_samples:
        Number of representative pixels sampled along each detector axis.
    :param int frame_samples:
        Maximum number of representative scan frames.
    :returns:
        Three geometry-matched axis steps.
    :rtype: tuple[float, float, float]
    :raises ValueError:
        If the percentile or sampling counts are invalid, or the scan does not
        span sufficient reciprocal-space dimensions.
    """
    percentile = float(percentile)
    if not 0.0 < percentile <= 100.0:
        raise ValueError("Resolution percentile must be above 0 and at most 100")
    detector_samples = int(detector_samples)
    frame_samples = int(frame_samples)
    if detector_samples < 2 or frame_samples < 1:
        raise ValueError(
            "Use at least two detector samples and one frame sample"
        )

    detector = config.detector
    rows, columns = detector.detector.shape
    bounds = np.asarray(
        scan.exposure_angle_bounds(config, fallback="stationary"),
        dtype=np.float64,
    )
    if bounds.shape != (len(scan), 2, 4):
        raise ValueError(
            "Exposure angle bounds must have shape (frames, 2, 4)"
        )
    excluded = set(config.corrections.excluded_frames)
    included = np.asarray(
        [index for index in range(len(scan)) if index not in excluded],
        dtype=np.int64,
    )
    if included.size == 0:
        raise ValueError("No included frames are available for resolution sampling")

    detector_rows = np.unique(
        np.rint(
            np.linspace(0, rows - 1, min(detector_samples, rows))
        ).astype(np.int64)
    )
    detector_columns = np.unique(
        np.rint(
            np.linspace(0, columns - 1, min(detector_samples, columns))
        ).astype(np.int64)
    )
    sampled_positions = np.unique(
        np.rint(
            np.linspace(0, included.size - 1, min(frame_samples, included.size))
        ).astype(np.int64)
    )
    sampled_frames = included[sampled_positions]

    from .datautils.xrayutils.reconstruction import _native_module

    ub = config.ub_calculator
    kernel = _native_module().ReconstructionKernel(
        np.full(3, -1e12),
        np.ones(3),
        # Only kernel.coordinate() is called here, which never looks at the
        # grid shape; it just has to be large enough to stand in for
        # "unbounded" and small enough for the kernel's packed voxel
        # identifier (64 bits across the three axes).
        np.full(3, 2_000_000, dtype=np.int64),
        np.ones(3, dtype=np.int64),
        frame,
        float(ub.getK()),
        np.ascontiguousarray(np.linalg.inv(ub.getUB())),
        np.ascontiguousarray(np.linalg.inv(ub.getU())),
        0,
        1,
        1,
        1024 * 1024,
    )

    def pixel_rays(row, column):
        ray_rows = np.asarray([row - 0.5, row + 0.5])
        ray_columns = np.asarray([column - 0.5, column + 0.5])
        row_grid, column_grid = np.meshgrid(
            ray_rows, ray_columns, indexing="ij"
        )
        gamma, delta = detector.primBeamPoints(row_grid, column_grid)
        cosine_gamma = np.cos(gamma)
        rays = np.empty((2, 2, 3), dtype=np.float64)
        rays[..., 0] = np.sin(delta) * cosine_gamma
        rays[..., 1] = np.cos(delta) * cosine_gamma
        rays[..., 2] = np.sin(gamma)
        rays /= np.linalg.norm(rays, axis=-1, keepdims=True)
        return np.ascontiguousarray(rays)

    def coordinate(rays, frame_index, u, v, t):
        return kernel.coordinate(
            rays,
            np.ascontiguousarray(bounds[frame_index, 0]),
            np.ascontiguousarray(bounds[frame_index, 1]),
            0,
            0,
            u,
            v,
            t,
        )

    local_sigmas = []
    scan_spans = []
    for frame_index in sampled_frames:
        included_position = int(np.searchsorted(included, frame_index))
        neighbor_index = None
        if included.size > 1:
            neighbor_position = (
                included_position + 1
                if included_position + 1 < included.size
                else included_position - 1
            )
            neighbor_index = int(included[neighbor_position])
        for row in detector_rows:
            for column in detector_columns:
                rays = pixel_rays(int(row), int(column))
                row_vector = (
                    coordinate(rays, frame_index, 1.0, 0.5, 0.5)
                    - coordinate(rays, frame_index, 0.0, 0.5, 0.5)
                )
                column_vector = (
                    coordinate(rays, frame_index, 0.5, 1.0, 0.5)
                    - coordinate(rays, frame_index, 0.5, 0.0, 0.5)
                )
                scan_vector = (
                    coordinate(rays, frame_index, 0.5, 0.5, 1.0)
                    - coordinate(rays, frame_index, 0.5, 0.5, 0.0)
                )
                if (
                    np.linalg.norm(scan_vector)
                    <= np.finfo(np.float64).eps
                    and neighbor_index is not None
                ):
                    scan_vector = (
                        coordinate(rays, neighbor_index, 0.5, 0.5, 0.5)
                        - coordinate(rays, frame_index, 0.5, 0.5, 0.5)
                    )
                scan_spans.append(np.linalg.norm(scan_vector))
                jacobian = np.column_stack(
                    (row_vector, column_vector, scan_vector)
                )
                local_sigmas.append(
                    np.sqrt(np.sum(jacobian * jacobian, axis=1) / 12.0)
                )

    local_sigmas = np.asarray(local_sigmas, dtype=np.float64)
    if not scan_spans or max(scan_spans) <= np.finfo(np.float64).eps:
        raise ValueError(
            "The scan has no exposure sweep or adjacent-frame displacement; "
            "a three-dimensional geometry-matched step cannot be estimated"
        )
    steps = np.percentile(local_sigmas, percentile, axis=0)
    if np.any(~np.isfinite(steps)) or np.any(steps <= 0):
        raise ValueError(
            "The sampled geometry does not span all three output axes; "
            "geometry-matched 3-D steps cannot be estimated"
        )
    return tuple(float(value) for value in steps)


def estimate_checkpoint_plan(
    config,
    scan,
    grids,
    *,
    max_depth,
    threads,
    ram_budget_bytes,
    checkpoint_count,
    angle_fallback="stationary",
    mask=None,
    budget_seconds=0.1,
    memory_budget_bytes=None,
    extra_excluded_frames=(),
):
    """Estimate per-grid job data volume and checkpoint file counts.

    Mirrors :func:`derive_grid`/:func:`estimate_geometry_steps`'s sampling
    idiom -- a single representative included frame's real angle bounds and
    a real detector mask if the caller has one, with dummy pixel values,
    since mapping cost and voxel occupancy depend only on geometry and the
    mask (see the design doc's calibration-probe rationale). Bounded to a
    small wall-time budget, so it is safe to call live, before a job is
    prepared, e.g. from the GUI on every settings change.

    :param config:
        Central :class:`ConfigData` experiment snapshot.
    :param scan:
        Active scan backend providing exposure angle bounds in radians and
        ``get_raw_img`` (not called by this function -- no image data is
        loaded).
    :param grids:
        Iterable of :class:`ReconstructionGrid` or grid-spec-compatible
        mappings, same convention as :func:`prepare_job`'s ``grids``
        parameter.
    :param int max_depth:
        Adaptive footprint-splitting depth the estimate should reflect.
    :param int threads:
        Native kernel thread count the estimate should reflect.
    :param float ram_budget_bytes:
        Combined memory budget for the whole job -- the user's own
        resource request, never a hardcoded value.
    :param int checkpoint_count:
        User-requested minimum number of checkpoint files.
    :param str angle_fallback:
        ``"stationary"`` or ``"midpoint"``, passed through to
        ``scan.exposure_angle_bounds``.
    :param mask:
        Optional boolean detector mask (``True`` = excluded), shape
        matching ``config.detector.detector.shape``. If ``None``, every
        pixel is treated as valid. The estimate is only as good as the
        mask it is given -- callers with a real mask available (a live
        GUI's detector mask, or a prepared job's assets) should supply it.
    :param float budget_seconds:
        Wall-time budget for the underlying calibration probe.
    :param memory_budget_bytes:
        Optional override for the native kernel's own memory-budget
        precheck; defaults to ``ram_budget_bytes``.
    :param extra_excluded_frames:
        Additional frame indices to treat as excluded, unioned with
        ``config.corrections.excluded_frames``. Used to scope the estimate
        to one cluster node's own slice of the scan (design doc Sec13)
        without touching the job-level exclusion set; empty for the
        single-node case.
    :returns:
        Dict with ``per_grid`` (mapping ``grid_name`` ->
        ``{"job_data_bytes_estimate", "files_per_job"}``) and
        ``files_total`` (sum of ``files_per_job`` across grids).
    :rtype: dict
    :raises ValueError:
        If ``checkpoint_count``/``ram_budget_bytes`` are not positive, no
        frames are included, or ``mask``'s shape does not match the
        detector.
    """
    checkpoint_count = int(checkpoint_count)
    if checkpoint_count < 1:
        raise ValueError("checkpoint_count must be at least one")
    if ram_budget_bytes <= 0:
        raise ValueError("ram_budget_bytes must be positive")

    detector = config.detector
    rows, columns = detector.detector.shape
    if mask is None:
        mask = np.zeros((rows, columns), dtype=bool)
    else:
        mask = np.ascontiguousarray(mask, dtype=bool)
        if mask.shape != (rows, columns):
            raise ValueError("mask shape must match the detector shape")

    excluded = set(config.corrections.excluded_frames) | set(extra_excluded_frames)
    included = [index for index in range(len(scan)) if index not in excluded]
    if not included:
        raise ValueError("No included frames are available for estimation")
    representative_frame = included[len(included) // 2]

    bounds = scan.exposure_angle_bounds(config, fallback=angle_fallback)
    angles_start = np.ascontiguousarray(bounds[representative_frame, 0])
    angles_end = np.ascontiguousarray(bounds[representative_frame, 1])
    corner_rays = _detector_corner_rays(detector, (0, rows, 0, columns))

    grid_values = [
        asdict(grid) if isinstance(grid, ReconstructionGrid) else dict(grid)
        for grid in grids
    ]
    if not grid_values:
        raise ValueError("At least one reconstruction grid is required")
    grid_specs = tuple(_GridSpec(**values) for values in grid_values)
    spec = _ReconstructionSpec(
        grids=grid_specs,
        max_depth=max_depth,
        threads=threads,
        memory_budget_bytes=(
            memory_budget_bytes
            if memory_budget_bytes is not None
            else ram_budget_bytes
        ),
    )

    probe = _calibration_probe_all_grids(
        spec,
        config.ub_calculator,
        mask,
        corner_rays,
        angles_start,
        angles_end,
        budget_seconds=budget_seconds,
        kernel_threads=threads,
    )

    valid_pixels = int((~mask).sum())
    frame_count = len(included)
    per_grid = {}
    files_total = 0
    for grid_name, result in probe.items():
        job_data_bytes_estimate = (
            result["records_per_pixel"]
            * _CHECKPOINT_BYTES_PER_ROW
            * valid_pixels
            * frame_count
        )
        files_per_job = _files_per_job(
            job_data_bytes_estimate, ram_budget_bytes, checkpoint_count
        )
        per_grid[grid_name] = {
            "job_data_bytes_estimate": job_data_bytes_estimate,
            "files_per_job": files_per_job,
        }
        files_total += files_per_job

    return {"per_grid": per_grid, "files_total": files_total}


def _included_ranges(count, excluded, batch_size):
    included = [index for index in range(count) if index not in excluded]
    ranges = []
    start = None
    previous = None
    for index in included:
        if (
            start is None
            or index != previous + 1
            or index - start >= batch_size
        ):
            if start is not None:
                ranges.append((start, previous + 1))
            start = index
        previous = index
    if start is not None:
        ranges.append((start, previous + 1))
    return ranges


def _checkpoint_frame_ranges(count, excluded, files_per_job):
    """Contiguous, exclusion-respecting frame ranges for one grid's planned
    checkpoints (design doc Sec5/Sec11).

    Splits all included frames into ``files_per_job`` roughly equal
    contiguous groups by reusing :func:`_included_ranges`'s existing
    gap-respecting batching. An excluded-frame gap forces an extra
    boundary regardless of target size, which can produce a few more
    ranges than ``files_per_job`` -- consistent with ``files_per_job``
    being a floor, not an exact count (Sec4).
    """
    included_count = sum(1 for index in range(count) if index not in excluded)
    if included_count == 0:
        return []
    batch_size = max(1, math.ceil(included_count / max(1, int(files_per_job))))
    return _included_ranges(count, excluded, batch_size)


def _node_excluded_frames(scan_length, excluded, total_tasks, task_index) -> set[int]:
    """Frames outside this cluster node's disjoint share of the included
    frames, unioned with the job's own excluded frames (design doc Sec13).

    Splits the *included* frame list into ``total_tasks`` contiguous
    shares by position (``included[i*n//total_tasks : (i+1)*n//total_tasks]``),
    not by re-running :func:`_included_ranges`'s gap-driven batching --
    that split is sized to bound checkpoint/task granularity and can
    produce more pieces than ``total_tasks`` when excluded frames create
    extra gaps, which would leave a piece owned by no node. Position-based
    splitting always yields exactly ``total_tasks`` shares (some possibly
    empty), covering every included frame exactly once. A share can still
    be internally non-contiguous (excluded frames inside it); that is
    fine -- :func:`_included_ranges` already decomposes a sparse included
    set into multiple ``(start, stop)`` ranges everywhere else it is used.

    :returns:
        The set of frame indices this node does *not* own -- pass to
        ``extra_excluded_frames`` on :func:`_execution_layout` and
        :func:`estimate_checkpoint_plan`, or union directly into an
        ``excluded`` argument to :func:`_checkpoint_frame_ranges`.
    """
    included = [index for index in range(scan_length) if index not in excluded]
    total_tasks = max(1, int(total_tasks))
    task_index = int(task_index)
    share_start = (task_index * len(included)) // total_tasks
    share_stop = ((task_index + 1) * len(included)) // total_tasks
    owned = set(included[share_start:share_stop])
    return (set(range(scan_length)) - owned) | set(excluded)


def prepare_job(
    gui,
    job_path,
    *,
    grids,
    scratch_path,
    output_path,
    accuracy="balanced",
    advanced_depth=None,
    compression_override=None,
    angle_fallback="stationary",
    user_note="",
    thread_override=None,
    memory_override_bytes=None,
    frame_batch=None,
    tile_shape=None,
    work_block_pixels=None,
    threads_per_image=None,
    accumulation_budget_bytes=None,
    checkpoint_count=10,
    cluster_settings=None,
):
    """Freeze current orGUI state into an immutable reconstruction job."""
    if gui.fscan is None:
        raise RuntimeError("Load a scan before preparing reconstruction")
    if accuracy not in ACCURACY_DEPTHS:
        raise ValueError(f"Unknown accuracy preset: {accuracy}")
    if angle_fallback not in {"stationary", "midpoint"}:
        raise ValueError("Angle fallback must be stationary or midpoint")
    if advanced_depth is not None and not 0 <= advanced_depth <= 8:
        raise ValueError("Advanced split depth must be between 0 and 8")
    if compression_override is not None and compression_override not in FILTERS:
        raise ValueError(
            f"Unknown HDF5 compression override: {compression_override}"
        )
    if threads_per_image is not None and int(threads_per_image) < 1:
        raise ValueError("Threads per image must be at least one")
    if int(checkpoint_count) < 1:
        raise ValueError("checkpoint_count must be at least one")
    if (
        accumulation_budget_bytes is not None
        and int(accumulation_budget_bytes) < 1024**2
    ):
        raise ValueError(
            "Per-worker accumulation budget must be at least 1 MiB"
        )
    grid_values = [
        asdict(grid) if isinstance(grid, ReconstructionGrid) else dict(grid)
        for grid in grids
    ]
    if not grid_values:
        raise ValueError("At least one reconstruction grid is required")
    for values in grid_values:
        _GridSpec(**values)
    config = ConfigData.from_gui(gui)
    for monitor in config.corrections.monitor_corrections:
        if not hasattr(gui.fscan, monitor):
            raise ValueError(
                f"Active scan has no monitor counter named {monitor!r}"
            )
    reference = ScanReference.from_scan(gui.fscan)
    scratch = Path(scratch_path).absolute()
    for component in (scratch, *scratch.parents):
        if component.exists() and not component.is_dir():
            raise ValueError(
                "Scratch directory cannot be created because a path "
                f"component is a file: {component}"
            )
    scratch.mkdir(parents=True, exist_ok=True)
    assets = scratch / "job-assets.nxs"
    assets_sha256 = _snapshot_assets(gui, config, assets)
    fingerprint = sha256(
        json.dumps(
            reference.to_dict(), sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    job = ReconstructionJob(
        config=config_data_to_json(config),
        scan_reference=reference.to_dict(),
        grids=grid_values,
        scratch_path=str(scratch),
        output_path=str(Path(output_path).absolute()),
        compression=(
            compression_override
            if compression_override is not None
            else _compression_name(gui.database.compression)
        ),
        assets_path=str(assets),
        assets_sha256=assets_sha256,
        source_fingerprint_sha256=fingerprint,
        threads_per_image=(
            None if threads_per_image is None else int(threads_per_image)
        ),
        accumulation_budget_bytes=(
            None
            if accumulation_budget_bytes is None
            else int(accumulation_budget_bytes)
        ),
        checkpoint_count=int(checkpoint_count),
        build_metadata=_build_metadata(),
        runtime_threads=max(1, int(gui.numberthreads)),
        runtime_memory_bytes=max(1, int(gui.maxMemory * 1024 * 1024)),
        accuracy=accuracy,
        advanced_depth=advanced_depth,
        angle_fallback=angle_fallback,
        user_note=user_note,
        thread_override=thread_override,
        memory_override_bytes=memory_override_bytes,
        frame_batch=frame_batch,
        tile_shape=tile_shape,
        work_block_pixels=work_block_pixels,
        cluster_settings=dict(cluster_settings or {}),
    )
    depth = (
        job.advanced_depth
        if job.advanced_depth is not None
        else ACCURACY_DEPTHS[job.accuracy]
    )
    threads = job.thread_override or job.runtime_threads
    memory = job.memory_override_bytes or job.runtime_memory_bytes
    mask = _load_assets(job).get("mask")
    plan = estimate_checkpoint_plan(
        config,
        gui.fscan,
        grid_values,
        max_depth=depth,
        threads=threads,
        ram_budget_bytes=memory,
        checkpoint_count=job.checkpoint_count,
        angle_fallback=angle_fallback,
        mask=mask,
    )
    excluded = set(config.corrections.excluded_frames)
    job.checkpoint_plan = {
        grid_name: [
            list(item)
            for item in _checkpoint_frame_ranges(
                len(gui.fscan), excluded, estimate["files_per_job"]
            )
        ]
        for grid_name, estimate in plan["per_grid"].items()
    }
    write_job(job, job_path)
    return job


def verify_job(job):
    """Verify immutable job assets and source fingerprints."""
    ScanReference.from_dict(job.scan_reference).verify()
    reference_digest = sha256(
        json.dumps(
            job.scan_reference, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    if reference_digest != job.source_fingerprint_sha256:
        raise RuntimeError("Reconstruction scan-reference digest mismatch")
    if _sha256_file(job.assets_path) != job.assets_sha256:
        raise RuntimeError("Reconstruction job assets checksum mismatch")
    return True


def _load_assets(job):
    result = {}
    with h5py.File(job.assets_path, "r") as h5file:
        for name in ("mask", "background", "background_variance"):
            if name in h5file:
                result[name] = h5file[name][()]
        if "repair" in h5file:
            repair = h5file["repair"]
            result["repair"] = {
                "max_component_pixels": int(
                    repair.attrs["max_component_pixels"]
                ),
                "max_span": int(repair.attrs["max_span"]),
                "radius": int(repair.attrs["radius"]),
                "min_valid_neighbors": int(
                    repair.attrs["min_valid_neighbors"]
                ),
                "row_gaps": repair["row_gaps"][()],
                "column_gaps": repair["column_gaps"][()],
            }
    return result


@lru_cache(maxsize=1)
def _correction_extension():
    """The native extension, or ``None`` where it is not installed.

    Correction has always had a pure-NumPy implementation and keeps one:
    unlike mapping, it is not a place the pipeline refuses to run without
    the extension, and the two must stay numerically interchangeable.
    """
    try:
        return importlib.import_module(
            "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
        )
    except ImportError:
        return None


def _correction_pipeline(config, scan, assets, provenance):
    correction = config.corrections
    detector = config.detector
    background = (
        np.asarray(assets["background"], dtype=np.float64)
        if correction.use_background and "background" in assets
        else None
    )
    background_variance = (
        np.asarray(assets["background_variance"], dtype=np.float64)
        if background is not None and "background_variance" in assets
        else None
    )
    static_mask = (
        np.ascontiguousarray(assets["mask"], dtype=bool)
        if correction.use_mask and "mask" in assets
        else None
    )
    static_factor = None
    if correction.use_solid_angle:
        static_factor = 1.0 / np.asarray(
            detector.solidAngleArray(), dtype=np.float64
        )
        provenance.setdefault("factor_uncertainty", {})[
            "solid_angle"
        ] = "deterministic-no-uncertainty"
    if correction.use_polarization:
        polarization = np.asarray(
            detector.polarization(
                factor=detector._polFactor, axis_offset=detector._polAxis
            ),
            dtype=np.float64,
        )
        if static_factor is None:
            static_factor = 1.0 / polarization
        else:
            static_factor /= polarization
        provenance.setdefault("factor_uncertainty", {})[
            "polarization"
        ] = "deterministic-no-uncertainty"
    if static_factor is not None:
        static_factor = np.ascontiguousarray(static_factor, dtype=np.float64)
        static_factor_squared = np.square(static_factor)
    else:
        static_factor_squared = None

    repair_plan = None
    if (
        correction.repair_masked_pixels
        and static_mask is not None
        and "repair" in assets
    ):
        repair_plan = create_pixel_repair_plan(
            static_mask, **dict(assets["repair"])
        )
        provenance["repair_plan"] = repair_plan.configuration()
        provenance["repair_covariance"] = "marginal-variance-only"

    exposure = (
        np.asarray(scan.exposure_time, dtype=np.float64)
        if correction.normalize_exposure and hasattr(scan, "exposure_time")
        else None
    )
    exposure_variance = (
        np.asarray(scan.exposure_time_variance, dtype=np.float64)
        if exposure is not None and hasattr(scan, "exposure_time_variance")
        else None
    )
    monitor_values = {
        name: np.asarray(getattr(scan, name), dtype=np.float64)
        for name in correction.monitor_corrections
    }
    monitor_variances = {
        name: (
            np.asarray(getattr(scan, f"{name}_variance"), dtype=np.float64)
            if hasattr(scan, f"{name}_variance")
            else None
        )
        for name in correction.monitor_corrections
    }

    def frame_value(values, frame_index):
        return (
            float(values)
            if values.ndim == 0
            else float(values[frame_index])
        )

    def apply_factor(intensity, variance, factor, factor_variance, name):
        if factor_variance is not None:
            variance *= factor**2
            variance += intensity**2 * factor_variance
            intensity *= factor
        else:
            intensity *= factor
            variance *= factor**2

    def record_factor(name, factor_variance):
        provenance.setdefault("factor_uncertainty", {})[name] = (
            "propagated"
            if factor_variance is not None
            else "deterministic-no-uncertainty"
        )

    def scalar_factors(frame_index):
        """This frame's exposure and monitor factors, in application order.

        :returns:
            ``[(factor, factor_squared, factor_variance_or_None, name)]``.
            ``factor_squared`` is carried rather than recomputed because
            it is what the NumPy form multiplied by, and ``x ** 2`` is
            ``pow`` -- the native pass must scale by the same value, not
            by one that agrees to the last bit on most platforms.
        :rtype: list[tuple]
        """
        found = []
        if exposure is not None:
            value = frame_value(exposure, frame_index)
            if value <= 0 or not math.isfinite(value):
                raise ValueError("Exposure time must be finite and positive")
            factor = 1.0 / value
            factor_variance = None
            if exposure_variance is not None:
                value_variance = frame_value(exposure_variance, frame_index)
                factor_variance = value_variance / value**4
            found.append((factor, factor**2, factor_variance, "exposure"))
        elif correction.normalize_exposure:
            provenance["exposure_normalization"] = "unavailable"
        for name, values in monitor_values.items():
            value = frame_value(values, frame_index)
            if value == 0 or not math.isfinite(value):
                raise ValueError(f"Monitor {name} must be finite and nonzero")
            monitor_variance = monitor_variances[name]
            factor_variance = None
            if monitor_variance is not None:
                value_variance = frame_value(monitor_variance, frame_index)
                factor_variance = value_variance / value**4
            factor = 1.0 / value
            found.append((factor, factor**2, factor_variance, f"monitor:{name}"))
        return found

    def apply_scaling(intensity, variance, mask, factors):
        """Scale, propagate and mask off non-finite pixels.

        Fused into one native pass where the extension is available. It
        is the same arithmetic on the same values in the same order --
        correction is element-wise, with nothing to reassociate -- so the
        two paths are bit-for-bit, and a test pins that. What changes is
        that eight or nine full-detector NumPy passes over ~50 MB each
        become one, and that the pass releases the GIL for its whole
        duration instead of taking and dropping it between operations.
        """
        native = _correction_extension()
        if native is not None and all(
            array.flags.c_contiguous for array in (intensity, variance, mask)
        ):
            native.apply_correction_factors(
                intensity,
                variance,
                mask,
                static_factor,
                static_factor_squared,
                [factor for factor, _squared, _variance, _name in factors],
                [squared for _factor, squared, _variance, _name in factors],
                [
                    math.nan if variance_ is None else variance_
                    for _factor, _squared, variance_, _name in factors
                ],
            )
            return mask
        if static_factor is not None:
            intensity *= static_factor
            variance *= static_factor_squared
        for factor, _squared, factor_variance, name in factors:
            apply_factor(intensity, variance, factor, factor_variance, name)
        mask |= ~np.isfinite(intensity) | ~np.isfinite(variance)
        return mask

    def correct_frame(payload, raw, frame_index):
        source_variance = getattr(payload, "variance", None)
        image_provenance = getattr(payload, "processing_provenance", None)
        if image_provenance:
            provenance.setdefault("image_processing", {}).update(
                dict(image_provenance)
            )
        if source_variance is None:
            variance = np.maximum(raw, 0.0)
            provenance["variance_fallback"] = "clipped-current-image-poisson"
        else:
            variance = np.asarray(source_variance, dtype=np.float64).copy()
        intensity = np.asarray(raw, dtype=np.float64).copy()
        if background is not None:
            intensity -= background
            if background_variance is not None:
                variance += background_variance
            else:
                provenance["background_variance"] = "deterministic"
        if static_mask is None:
            mask = np.zeros(intensity.shape, dtype=bool)
        else:
            mask = static_mask.copy()
        if repair_plan is not None:
            mask, _ = repair_plan.apply_inplace(intensity, variance)
        factors = scalar_factors(frame_index)
        for _factor, _squared, factor_variance, name in factors:
            record_factor(name, factor_variance)
        mask = apply_scaling(intensity, variance, mask, factors)
        return intensity, variance, mask

    def correct(payload, raw, frame_index, tile):
        """Compatibility call shape used by internal focused tests."""
        row_start, row_stop, column_start, column_stop = tile
        full_raw = np.asarray(getattr(payload, "img", raw))
        intensity, variance, mask = correct_frame(
            payload, full_raw, frame_index
        )
        selection = np.s_[
            row_start:row_stop, column_start:column_stop
        ]
        return (
            intensity[selection],
            variance[selection],
            mask[selection],
        )

    correct.correct_frame = correct_frame
    correct.repair_plan = repair_plan
    return correct


CHECKPOINT_MEMORY_SHARE = 0.25
"""Fraction of the memory budget reserved for in-memory checkpoint
accumulators; the rest funds the frame pipeline.

The two used to be sized against the whole budget independently -- the
worker pool took as many concurrent frames as the budget allowed, while
``MAX_CONCURRENT_ACTIVE_CHECKPOINTS`` accumulators were each allowed the
budget divided by that same count -- so nothing stopped their sum
reaching twice what the user asked for, and measured peaks did reach
1.4x. Splitting the budget once, here, is what makes the two shares add
up to it instead.

Accumulating longer before flushing mostly saves scratch writes, and
measurably little: quartering the share grew a checkpoint set by 6%.
"""


def split_memory_budget(memory_bytes, grid_count):
    """Divide a memory budget between checkpoint accumulators and pipeline.

    :param int memory_bytes:
        The job's whole memory budget.
    :param int grid_count:
        Number of output grids; each has its own active checkpoints.
    :returns:
        ``(accumulator_bytes_per_checkpoint, pipeline_bytes)`` -- the
        first bounds one active checkpoint's in-memory records, the
        second is what the frame pipeline (in-flight frames, prefetch
        queue, native buffers) may use.
    :rtype: tuple[int, int]
    """
    memory_bytes = max(2 * 1024**2, int(memory_bytes))
    accumulator_total = int(memory_bytes * CHECKPOINT_MEMORY_SHARE)
    per_checkpoint = max(
        1024**2,
        accumulator_total
        // (MAX_CONCURRENT_ACTIVE_CHECKPOINTS * max(1, int(grid_count))),
    )
    return per_checkpoint, max(1024**2, memory_bytes - accumulator_total)


def work_block_memory_cap(depth, memory_bytes, threads):
    """Largest work block whose native arenas stay inside the memory budget.

    The kernel reserves one arena per worker thread, sized for the worst
    case that every leaf of every pixel in a block reaches a distinct
    voxel: ``block * children**depth * 72`` bytes. That grows by the
    subdivision factor with depth, so a block size that is unremarkable at
    ``center`` accuracy reserves tens of gigabytes per thread at the
    deepest settings. Eight children are assumed because a spec cannot
    know whether an individual frame range's exposure rotates; a
    stationary one subdivides four ways and so has headroom to spare.

    :param int depth:
        Adaptive subdivision depth.
    :param int memory_bytes:
        The job's memory budget.
    :param int threads:
        Thread budget; arenas are live one per worker thread at once.
    :returns:
        Pixel count, at least 1.
    :rtype: int
    """
    leaves = min(8 ** max(0, int(depth)), _RESERVED_RECORDS_PER_PIXEL)
    bytes_per_block_pixel = leaves * _MAP_NODE_BYTES
    per_thread = max(1, int(memory_bytes)) // max(1, int(threads))
    return max(1, per_thread // _ARENA_BUDGET_SHARE // bytes_per_block_pixel)


def resolve_work_block_pixels(setting, depth, memory_bytes, threads):
    """Resolve a work-block setting to a pixel count.

    ``None`` selects the default preset, a string names one of
    :data:`WORK_BLOCK_PRESETS`, and an integer is taken literally. Preset
    (and default) sizes halve with each adaptive depth so that the name
    keeps selecting a working-set scale rather than a raw count; an
    explicit integer is what the caller asked for and is not rescaled.

    Every route is capped by :func:`work_block_memory_cap`, so no setting
    -- including a pinned integer -- can make the kernel reserve more
    memory than the job's budget allows.

    :raises ValueError:
        If ``setting`` names an unknown preset or is not positive.
    :rtype: int
    """
    if setting is None:
        scaled = WORK_BLOCK_PRESETS[DEFAULT_WORK_BLOCK] >> min(int(depth), 16)
    elif isinstance(setting, str):
        if setting not in WORK_BLOCK_PRESETS:
            raise ValueError(
                f"Unknown native work block preset: {setting!r}; choose one of "
                + ", ".join(
                    sorted(WORK_BLOCK_PRESETS, key=WORK_BLOCK_PRESETS.get)
                )
            )
        scaled = WORK_BLOCK_PRESETS[setting] >> min(int(depth), 16)
    else:
        scaled = int(setting)
        if scaled < 1:
            raise ValueError("Native work block pixels must be positive")
    return max(
        1, min(max(1, scaled), work_block_memory_cap(depth, memory_bytes, threads))
    )


def _execution_layout(job, scan, config, *, extra_excluded_frames=()):
    """Compute the automatic scheduling ranges and detector tiles.

    :param extra_excluded_frames:
        Additional frame indices to exclude, unioned with
        ``config.corrections.excluded_frames``. Used to scope scheduling
        to one cluster node's own slice of the scan (design doc Sec13);
        empty for the single-node case.
    """
    excluded = set(config.corrections.excluded_frames) | set(extra_excluded_frames)
    rows, columns = config.detector.detector.shape
    effective_memory = (
        job.memory_override_bytes or job.runtime_memory_bytes
    )
    depth = (
        job.advanced_depth
        if job.advanced_depth is not None
        else ACCURACY_DEPTHS[job.accuracy]
    )
    # Bounded by the record ceiling the arena and both memory prechecks
    # already use, not by the worst-case leaf count. 8**depth is what this
    # site was left with when the prechecks were corrected, and it is far
    # more conservative than anything a pixel actually leaves behind: at
    # depth 2 it claims 5248 bytes per pixel against a real ~106, which
    # bands this detector into 13 strips of 194 rows, and by depth 5 it
    # asks for 2.6 MB per pixel and collapses a band to a single row.
    #
    # Thin bands are not free. Measured at depth 0 on the reference job,
    # banding finer than the budget requires costs 10-15% -- more kernel
    # calls per frame group, each with less work to spread over its
    # threads, and a brick that is shorter in the row direction than the
    # locality it is trying to exploit.
    estimated_native_bytes_per_pixel = (
        128
        + 2
        * min(8**depth, _RESERVED_RECORDS_PER_PIXEL)
        * _CHECKPOINT_BYTES_PER_ROW
    )
    tile_pixels = max(
        1,
        min(
            1024**2,
            effective_memory
            // max(
                (len(job.grids) + 2) * estimated_native_bytes_per_pixel,
                1,
            ),
        ),
    )
    tile_rows, tile_columns = job.tile_shape or (
        # Full detector width, so a tile is a contiguous slice of the
        # corrected frame: the per-tile copy the mapping loop makes is
        # then a no-op instead of a strided gather, and every tile holds
        # the same number of pixels (and so the same number of native
        # work blocks) rather than the ragged remainders a square grid
        # leaves along two edges.
        max(1, min(rows, tile_pixels // max(1, columns))),
        columns,
    )
    included_count = max(
        1,
        len(scan)
        - len({frame for frame in excluded if 0 <= frame < len(scan)}),
    )
    if job.frame_batch is None:
        spec = job.internal_spec()
        seed_threads_per_image = (
            1 if job.threads_per_image is None else job.threads_per_image
        )
        native_threads = max(
            1, min(seed_threads_per_image, spec.threads)
        )
        image_workers = max(1, spec.threads // native_threads)
        target_tasks = max(1, image_workers * 4)
        frame_batch = min(
            AUTO_MAX_FRAMES_PER_TASK,
            max(1, math.ceil(included_count / target_tasks)),
        )
    else:
        frame_batch = job.frame_batch
    ranges = _included_ranges(len(scan), excluded, frame_batch)
    if job.tile_shape is None:
        # Spread the rows evenly over the bands the budget allows, rather
        # than filling each to tile_rows and leaving a short remainder:
        # every band then carries the same work to within one row.
        band_count = max(1, math.ceil(rows / tile_rows))
        edges = [
            round(rows * index / band_count) for index in range(band_count + 1)
        ]
        tiles = [
            (edges[index], edges[index + 1], 0, columns)
            for index in range(band_count)
        ]
    else:
        tiles = [
            (
                row,
                min(row + tile_rows, rows),
                column,
                min(column + tile_columns, columns),
            )
            for row in range(0, rows, tile_rows)
            for column in range(0, columns, tile_columns)
        ]
    return ranges, tiles


def _frame_parallelism(
    spec,
    tiles,
    memory_bytes,
    *,
    frames_per_task=1,
    threads_per_image=1,
    accumulation_budget_bytes=None,
):
    """Derive bounded parallelism across images and within each image.

    Each frame worker processes that frame's detector tiles sequentially,
    one native ``accumulate()`` call at a time -- never more than one tile's
    worth of native working memory is live at once for a given worker. The
    configured CPU budget is divided between concurrent image workers and
    native threads used by each image. The worker count is also bounded by a
    depth-aware working-set estimate scoped to that single largest tile
    (mirroring the native kernel's own per-call memory precheck exactly,
    ``ReconstructionKernel::accumulate`` in
    ``reciprocal_reconstruction_cpp.cpp``) -- summing every tile across the
    whole detector would overstate one worker's actual peak memory need by
    roughly the tile count, silently starving ``image_workers`` on jobs with
    more than one detector tile.

    The native-kernel estimate alone still misses two real, non-native
    memory costs (design doc Sec7): the Python-side correction step
    (``_correction_pipeline.correct_frame``), which runs once per frame
    rather than once per tile and holds the raw image plus every
    corrected array at full-detector size for as long as that frame is
    in flight; and the prefetch pipeline's own read-ahead queue, which
    can hold up to ``_PREFETCH_QUEUE_SLACK`` additional decoded (but not
    yet claimed) frames beyond ``image_workers``. Both are folded into
    this function's own worker/reservation accounting so ``image_workers``
    reflects the job's actual peak resident memory, not just the native
    kernel's.

    ``accumulation_budget_bytes``/the returned ``accumulation`` value model a
    per-worker retained-record buffer from the retired Parquet-era mapping
    path; the checkpoint-routing path (design doc Sec9) has no equivalent
    per-worker buffer of its own; ``run_job`` computes the checkpoint
    accumulators' own memory budget separately (Sec10) and does not forward
    this value onward. Kept only as a conservative extra bound on
    ``image_workers`` when a caller supplies it explicitly.

    :param _ReconstructionSpec spec:
        Frozen reconstruction compute settings.
    :param iterable tiles:
        Detector tiles as ``(row_start, row_stop, column_start,
        column_stop)`` tuples.
    :param int memory_bytes:
        Total memory budget in bytes.
    :param int frames_per_task:
        Number of images streamed by one resumable task. Images are not
        retained together.
    :param int threads_per_image:
        Requested native reconstruction threads for each concurrent image.
    :param accumulation_budget_bytes:
        Optional requested retained-record bytes per worker.
    :returns:
        ``(image_workers, kernel_threads, per_worker_memory_bytes,
        accumulation_budget_bytes)``.
    :rtype: tuple[int, int, int, int]
    """
    tiles = list(tiles)
    if not tiles:
        minimum = 1024**2
        return 1, 1, max(minimum, memory_bytes), minimum
    tile_pixels = [
        (row_stop - row_start) * (column_stop - column_start)
        for row_start, row_stop, column_start, column_stop in tiles
    ]
    # One worker's peak native working set is bounded by its single
    # largest tile (see docstring) -- not the sum of every tile across
    # the whole detector.
    largest_tile_pixels = max(tile_pixels)
    # The whole detector's pixel count: correct_frame (unlike the native
    # kernel) runs once per frame, not once per tile, and its buffers are
    # always full-detector-sized regardless of how many native tiles the
    # frame is split into.
    detector_pixels = sum(tile_pixels)
    # Mirrors ReconstructionKernel::accumulate's own memory precheck
    # exactly (a fixed per-pixel baseline plus twice the records a pixel is
    # assumed to leave behind, times one native Record's on-the-wire size),
    # so this Python-side estimate cannot drift from what the kernel itself
    # actually enforces.
    #
    # That estimate used to be the worst-case *leaf* count, 4**depth or
    # 8**depth, which is where the two could drift: a spec cannot know
    # whether an individual frame range's exposure rotates, so this side had
    # to assume the rotating case and reserve twice what the kernel would.
    # Bounding both by the same record ceiling removes the asymmetry
    # outright -- above depth 0 the bound saturates, so the estimate no
    # longer depends on the exposure model at all, and the two sides now
    # agree by construction rather than by careful mirroring.
    native_bytes_per_pixel = (
        128 + 2 * _RESERVED_RECORDS_PER_PIXEL * _CHECKPOINT_BYTES_PER_ROW
        if spec.max_depth > 0
        else 128 + 2 * _CHECKPOINT_BYTES_PER_ROW
    )
    # A frame group is one native call over (frames, tile rows, tile
    # columns) samples, and the kernel's own precheck is written against
    # that sample count -- so a worker's native working set scales with
    # the group size exactly as the tile size, and this side must scale
    # with it too or the two would disagree in the one direction that
    # matters (this side permissive, the kernel throwing).
    frames_per_group = _frames_per_group(spec)
    native_memory = (
        largest_tile_pixels * frames_per_group * native_bytes_per_pixel
    )
    # The native estimate alone misses real, per-in-flight-frame Python
    # memory: the raw decoded image plus every corrected array
    # (_PYTHON_CORRECTION_BYTES_PER_PIXEL), full-detector-sized.
    #
    # Only the correction half multiplies by the group size. A group
    # holds all F corrected frames at once -- each tile is a row band, so
    # every frame is revisited once per band and none can be released
    # early -- but it produces one merged record batch for the whole
    # group, not F of them, and that batch is measured at 0.575x the sum
    # of the F it replaces. Charging the record term per frame would
    # reserve nearly twice what a group actually holds.
    python_memory = detector_pixels * (
        frames_per_group * _PYTHON_CORRECTION_BYTES_PER_PIXEL
        + _FRAME_RECORD_BYTES_PER_PIXEL
    )
    image_memory = native_memory + python_memory
    worker_memory = max(
        1024**2,
        image_memory,
    )
    kernel_threads = max(
        1, min(int(threads_per_image), max(1, int(spec.threads)))
    )
    cpu_workers = max(1, int(spec.threads) // kernel_threads)
    minimum_accumulation = 1024**2
    # Design doc Sec7's prefetch-queue memory constraint: up to
    # _PREFETCH_QUEUE_SLACK frames beyond image_workers can sit in the
    # ready queue, already decoded but not yet claimed by a compute
    # worker -- raw image bytes only, no correction buffers yet. Reserved
    # off the top before dividing the remaining budget among workers, the
    # same way the ray-corner cache is reserved by callers before this
    # function ever sees the budget.
    # The queue holds groups, not frames, so a slack slot now costs a
    # whole group's raw images.
    prefetch_reserve_bytes = (
        _PREFETCH_QUEUE_SLACK
        * frames_per_group
        * detector_pixels
        * _RAW_IMAGE_BYTES_PER_PIXEL
    )
    usable_memory_bytes = max(
        minimum_accumulation, memory_bytes - prefetch_reserve_bytes
    )
    if accumulation_budget_bytes is None:
        memory_workers = max(1, usable_memory_bytes // worker_memory)
    else:
        requested = max(
            minimum_accumulation, int(accumulation_budget_bytes)
        )
        required = (
            worker_memory
            + ACCUMULATION_TRANSIENT_FACTOR * requested
        )
        memory_workers = max(1, usable_memory_bytes // required)
    image_workers = max(
        1,
        min(cpu_workers, memory_workers),
    )
    safe_accumulation = max(
        minimum_accumulation,
        (
            usable_memory_bytes // image_workers - worker_memory
        )
        // ACCUMULATION_TRANSIENT_FACTOR,
    )
    if accumulation_budget_bytes is None:
        accumulation = min(
            AUTO_MAX_ACCUMULATION_BYTES, safe_accumulation
        )
    else:
        accumulation = min(
            max(minimum_accumulation, int(accumulation_budget_bytes)),
            safe_accumulation,
        )
    per_worker_memory = (
        worker_memory
        + ACCUMULATION_TRANSIENT_FACTOR * accumulation
    )
    return (
        image_workers,
        kernel_threads,
        per_worker_memory,
        accumulation,
    )


def reconstruction_execution_settings(job, scan=None, config=None):
    """Return the effective map-task and native execution layout.

    Automatic values are derived with the same functions used by
    :func:`run_job`; unset advanced job fields therefore remain automatic.

    :param ReconstructionJob job:
        Prepared reconstruction job.
    :param scan:
        Optional already-open scan. The job scan reference is opened when
        omitted.
    :param ConfigData config:
        Optional decoded configuration snapshot.
    :returns:
        JSON-compatible detected execution settings.
    :rtype: dict
    """
    scan = job.scan if scan is None else scan
    config = job.config_data if config is None else config
    spec = job.internal_spec()
    # Automatic mode (job.threads_per_image is None, Sec7 Phase 4b) has no
    # live-measured I/O rate available before a job runs, so this preview
    # reports the same I/O-optimistic seed _map_pending_ranges itself
    # starts a live-balanced run from, not a claim about the pair that
    # will actually be chosen once real timing data exists.
    effective_threads_per_image = (
        1 if job.threads_per_image is None else job.threads_per_image
    )
    ranges, tiles = _execution_layout(job, scan, config)
    bounds = scan.exposure_angle_bounds(
        config, fallback=job.angle_fallback
    )
    ray_cache_bytes = sum(
        (tile[1] - tile[0] + 1)
        * (tile[3] - tile[2] + 1)
        * 3
        * np.dtype(np.float64).itemsize
        for tile in tiles
    )
    cache_detector_rays = (
        ray_cache_bytes <= spec.memory_budget_bytes // 8
    )
    scheduler_memory = max(
        1024**2,
        spec.memory_budget_bytes
        - (ray_cache_bytes if cache_detector_rays else 0),
    )
    layouts = {}
    for start, stop in ranges:
        task_bounds = bounds[start:stop]
        stationary = bool(
            np.array_equal(task_bounds[:, 0], task_bounds[:, 1])
        )
        (
            image_workers,
            kernel_threads,
            per_worker_memory,
            accumulation_budget,
        ) = (
            _frame_parallelism(
                spec,
                tiles,
                scheduler_memory,
                frames_per_task=stop - start,
                threads_per_image=effective_threads_per_image,
                accumulation_budget_bytes=job.accumulation_budget_bytes,
            )
        )
        key = (
            "stationary" if stationary else "swept",
            image_workers,
            kernel_threads,
            per_worker_memory,
            accumulation_budget,
        )
        layouts[key] = {
            "exposure": key[0],
            "concurrent_image_workers": image_workers,
            "native_threads_per_image": kernel_threads,
            "tiles_per_image": len(tiles),
            "memory_per_image_MiB": per_worker_memory / 1024**2,
            "accumulation_MiB_per_worker": (
                accumulation_budget / 1024**2
            ),
        }
    if tiles:
        tile_shape = (
            max(row_stop - row_start for row_start, row_stop, _, _ in tiles),
            max(
                column_stop - column_start
                for _, _, column_start, column_stop in tiles
            ),
        )
    else:
        tile_shape = (0, 0)
    return {
        "thread_budget": spec.threads,
        "native_threads_per_image": max(
            1, min(effective_threads_per_image, spec.threads)
        ),
        "threads_per_image_mode": (
            "automatic" if job.threads_per_image is None else "pinned"
        ),
        "memory_budget_MiB": spec.memory_budget_bytes / 1024**2,
        "accumulation_budget_MiB_per_worker": min(
            (
                layout["accumulation_MiB_per_worker"]
                for layout in layouts.values()
            ),
            default=1.0,
        ),
        "frames_per_task": max(
            (stop - start for start, stop in ranges), default=0
        ),
        "detector_tile_shape": tile_shape,
        "native_work_block_pixels": spec.work_block_pixels,
        "checkpoint_count": spec.checkpoint_count,
        "frame_tasks": len(ranges),
        "detector_tiles": len(tiles),
        "map_tasks": len(ranges),
        "detector_ray_cache_MiB": (
            ray_cache_bytes / 1024**2 if cache_detector_rays else 0.0
        ),
        "parallel_layouts": list(layouts.values()),
    }


def _job_checkpoint_boundaries(job):
    return {
        grid_name: [tuple(item) for item in items]
        for grid_name, items in job.checkpoint_plan.items()
    }


def _cluster_status_counts(job, *, total_tasks):
    """Best-effort checkpoint completion counts for a cluster job.

    Only reads whatever node sidecars already exist on disk -- a node
    that has not run yet simply contributes 0 planned/0 completed rather
    than triggering a calibration probe just to answer a status query.
    Authoritative completeness is still :func:`run_cluster_finalize`'s
    own check.
    """
    total_checkpoints = 0
    completed_checkpoints = 0
    for task_index in range(int(total_tasks)):
        sidecar = _read_node_checkpoint_sidecar(
            job, task_index, total_tasks=total_tasks
        )
        if sidecar is None:
            continue
        boundaries = {
            grid_name: [tuple(item) for item in items]
            for grid_name, items in sidecar["boundaries"].items()
        }
        total_checkpoints += sum(len(items) for items in boundaries.values())
        node_dir = _node_checkpoint_dir(job, task_index)
        resumed, _files = _discover_checkpoint_state(
            node_dir, boundaries, job.digest
        )
        completed_checkpoints += len(resumed)
    return total_checkpoints, completed_checkpoints


def job_status(path):
    """Return verified completion state for a job JSON."""
    job = read_job(path)
    array_task_count = job.cluster_settings.get("array_task_count")
    if job.status != "complete":
        verify_job(job)
    if array_task_count and job.status != "complete":
        # Node sidecars are deleted on successful finalize (matching the
        # single-node cleanup), so a completed cluster job has nothing
        # left on disk to count -- `status`/`output_sha256` already tell
        # the real story at that point.
        total_checkpoints, completed_checkpoints = _cluster_status_counts(
            job, total_tasks=array_task_count
        )
    elif array_task_count:
        total_checkpoints = completed_checkpoints = 0
    else:
        boundaries = _job_checkpoint_boundaries(job)
        total_checkpoints = sum(len(items) for items in boundaries.values())
        if job.status != "complete":
            checkpoint_dir = Path(job.scratch_path) / "checkpoints"
            resumed, _files = _discover_checkpoint_state(
                checkpoint_dir, boundaries, job.digest
            )
            completed_checkpoints = len(resumed)
        else:
            completed_checkpoints = total_checkpoints
    result = {
        "status": job.status,
        "job_sha256": job.digest,
        "checkpoints": {
            "completed": completed_checkpoints,
            "pending": max(0, total_checkpoints - completed_checkpoints),
            "total": total_checkpoints,
        },
        "output_path": job.output_path,
        "output_sha256": job.output_sha256,
        "grids": [
            {
                "name": _GridSpec(**grid).grid_name,
                "frame": grid["frame"],
                "shape": _GridSpec(**grid).shape,
            }
            for grid in job.grids
        ],
        "cleanup_errors": job.cleanup_errors,
    }
    return result


def _base_provenance(job, config):
    return {
        **job.correction_provenance,
        "job_sha256": job.digest,
        "user_note": job.user_note,
        "source_fingerprint_sha256": job.source_fingerprint_sha256,
        "correction_state": config.corrections.to_dict(),
        "cross_voxel_covariance": "marginal variances only",
        "native_build": job.build_metadata,
    }


def _merge_provenance(target, source):
    """Deep-merge one node's correction provenance into another's.

    Every cluster node runs :func:`_correction_pipeline` in its own
    process, mutating its own provenance dict independently -- unlike the
    single-node path, where one dict is shared across every frame within
    one process, cluster provenance genuinely needs merging back together
    at finalize time (design doc Sec13).
    """
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _merge_provenance(target[key], value)
        else:
            target[key] = value


def _execution_spec(job, *, threads=None, memory_bytes=None):
    """Build the job's spec, optionally overriding threads/memory for one
    execution.

    A cluster array task's own scheduler allocation (``cpus``,
    ``memory_bytes``) can differ from the job's captured
    ``runtime_threads``/``runtime_memory_bytes`` (design doc Sec13); this
    lets :func:`run_cluster_map_task` build a spec reflecting what that
    task actually has available, without mutating the job's own recorded
    configuration. Returns ``job.internal_spec()`` unchanged when both
    overrides are ``None`` (the single-node case).
    """
    spec = job.internal_spec()
    if threads is None and memory_bytes is None:
        return spec
    overrides = {}
    if threads is not None:
        overrides["threads"] = max(1, int(threads))
    if memory_bytes is not None:
        overrides["memory_budget_bytes"] = max(1024**2, int(memory_bytes))
    return _ReconstructionSpec.from_dict({**spec.to_dict(), **overrides})


def _node_checkpoint_dir(job, task_index):
    """Scratch subdirectory owning one cluster node's checkpoint files
    (design doc Sec13) -- disjoint from every other node's, and from the
    single-node case's flat ``checkpoints/{grid_name}/...`` layout."""
    return Path(job.scratch_path) / "checkpoints" / f"node{int(task_index):04d}"


def _read_node_checkpoint_sidecar(job, task_index, *, total_tasks):
    """Read one node's persisted checkpoint-plan sidecar, if present and
    valid for the current job/total_tasks/task_index -- no scan/config
    access, no computation. Used by :func:`_node_checkpoint_plan` (which
    falls back to computing it) and by :func:`job_status` (which reports
    "not started yet" rather than paying for a calibration probe just to
    answer a status query).

    :returns:
        ``{"excluded": [...], "boundaries": {...}}`` or ``None``.
    :rtype: dict | None
    """
    sidecar_path = _node_checkpoint_dir(job, task_index) / "plan.json"
    if not sidecar_path.exists():
        return None
    try:
        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
        if (
            sidecar.get("job_sha256") == job.digest
            and sidecar.get("total_tasks") == int(total_tasks)
            and sidecar.get("task_index") == int(task_index)
        ):
            return {
                "excluded": sidecar["excluded"],
                "boundaries": sidecar["boundaries"],
            }
    except (OSError, ValueError, KeyError):
        pass
    return None


def _node_checkpoint_plan(job, scan, config, mask, *, total_tasks, task_index):
    """Load or compute-and-persist one cluster node's checkpoint plan.

    Mirrors :func:`prepare_job`'s single-node checkpoint-plan computation,
    scoped to this node's own disjoint frame share
    (:func:`_node_excluded_frames`) and persisted as a small JSON sidecar
    under this node's own scratch subdirectory instead of the job JSON
    (design doc Sec11: computed once, never recomputed on resume --
    extended to per-node scope so a scheduler retry of the same array
    task reuses the same boundaries rather than risking different ones
    from calibration-probe sampling variance).

    :returns:
        ``None`` if this node owns no frames (``task_index`` beyond the
        real share count for ``total_tasks``). Otherwise a dict with
        ``excluded`` (sorted list of frame indices this node does not
        own) and ``boundaries`` (``{grid_name: [[start, stop], ...]}``,
        this node's own planned checkpoint frame ranges).
    :rtype: dict | None
    """
    node_excluded = _node_excluded_frames(
        len(scan), set(config.corrections.excluded_frames), total_tasks, task_index
    )
    if len(node_excluded) >= len(scan):
        return None

    cached = _read_node_checkpoint_sidecar(
        job, task_index, total_tasks=total_tasks
    )
    if cached is not None:
        return cached

    depth = (
        job.advanced_depth
        if job.advanced_depth is not None
        else ACCURACY_DEPTHS[job.accuracy]
    )
    threads = job.thread_override or job.runtime_threads
    memory = job.memory_override_bytes or job.runtime_memory_bytes
    plan = estimate_checkpoint_plan(
        config,
        scan,
        job.grids,
        max_depth=depth,
        threads=threads,
        ram_budget_bytes=memory,
        checkpoint_count=job.checkpoint_count,
        angle_fallback=job.angle_fallback,
        mask=mask,
        extra_excluded_frames=node_excluded,
    )
    boundaries = {
        grid_name: [
            list(item)
            for item in _checkpoint_frame_ranges(
                len(scan), node_excluded, estimate["files_per_job"]
            )
        ]
        for grid_name, estimate in plan["per_grid"].items()
    }
    result = {"excluded": sorted(node_excluded), "boundaries": boundaries}
    sidecar_path = _node_checkpoint_dir(job, task_index) / "plan.json"
    _atomic_json(
        sidecar_path,
        {
            "job_sha256": job.digest,
            "total_tasks": int(total_tasks),
            "task_index": int(task_index),
            **result,
        },
    )
    return result


def _range_is_resumed(router, boundaries, resumed, frame_range):
    """Whether every frame in ``frame_range`` is, for every grid, already
    covered by an already-resumed checkpoint. Shared by :func:`run_job`
    and :func:`run_cluster_map_task`."""
    start, stop = frame_range
    for grid_name in boundaries:
        for frame_index in range(start, stop):
            index = router.checkpoint_index_for_frame(grid_name, frame_index)
            if (grid_name, index) not in resumed:
                return False
    return True


class _AdjustablePool:
    """A resizable pool of raw worker threads.

    ``concurrent.futures.ThreadPoolExecutor`` cannot be resized once
    created; both the prefetch pool (grown/shrunk continuously by
    blocked-fraction, design doc Sec7) and, in a later phase, the compute
    pool (resized rarely, on a joint ``kernel_threads``/``image_workers``
    rebalance) need to change worker count live, so both use this
    instead.

    ``worker_fn(retire, cancellation)`` runs one worker's entire loop; it
    must poll both events at every blocking boundary (never a bare
    blocking call) and only exit *between* work items, never mid-item, so
    a shrink or a cancellation never abandons in-flight work.
    """

    def __init__(self, worker_fn, *, initial_size, name):
        self._worker_fn = worker_fn
        self._name = name
        self._cancellation = threading.Event()
        self._lock = threading.Lock()
        self._active: list[tuple[threading.Thread, threading.Event]] = []
        self._retiring: list[threading.Thread] = []
        self._next_id = 0
        self.retarget(initial_size)

    @property
    def size(self) -> int:
        """Current target worker count (excludes stragglers still
        finishing their last item after a shrink)."""
        with self._lock:
            return len(self._active)

    def _spawn(self):
        retire = threading.Event()
        thread = threading.Thread(
            target=self._worker_fn,
            args=(retire, self._cancellation),
            name=f"{self._name}-{self._next_id}",
            daemon=True,
        )
        self._next_id += 1
        thread.start()
        return thread, retire

    def retarget(self, size) -> None:
        """Grow to ``size`` workers, or mark the excess for retirement."""
        size = max(0, int(size))
        with self._lock:
            if self._cancellation.is_set():
                return
            current = len(self._active)
            if size > current:
                self._active.extend(self._spawn() for _ in range(size - current))
            elif size < current:
                retiring = self._active[size:]
                self._active = self._active[:size]
                for thread, retire in retiring:
                    retire.set()
                    self._retiring.append(thread)

    def reap(self) -> None:
        """Drop worker threads that have finished retiring. Call this
        periodically from the coordinator loop; retired threads are not
        joined eagerly, so the caller never stalls waiting for one to
        finish its in-flight item."""
        with self._lock:
            self._retiring = [thread for thread in self._retiring if thread.is_alive()]

    def shutdown(self, *, wait=True) -> None:
        with self._lock:
            self._cancellation.set()
            threads = [thread for thread, _retire in self._active] + self._retiring
            self._active = []
            self._retiring = []
        if wait:
            for thread in threads:
                thread.join()


class _BoundedGate:
    """A semaphore-like backpressure gate whose capacity can be retargeted
    live, decoupling a bounded in-flight count from the pool objects on
    either side of it (design doc Sec7: prefetch queue depth
    ``N ~= image_workers + small_constant``) -- retargeting it never
    requires tearing down the reader or compute pools.
    """

    def __init__(self, capacity):
        self._condition = threading.Condition()
        self._capacity = max(1, int(capacity))
        self._in_flight = 0

    def acquire(self, *, poll_timeout, should_stop) -> bool:
        """Block (polling ``should_stop``) until under capacity, then
        reserve a slot. Returns ``False`` without reserving if
        ``should_stop`` fires first."""
        with self._condition:
            while self._in_flight >= self._capacity:
                if should_stop():
                    return False
                self._condition.wait(timeout=poll_timeout)
            if should_stop():
                return False
            self._in_flight += 1
            return True

    def release(self) -> None:
        with self._condition:
            self._in_flight = max(0, self._in_flight - 1)
            self._condition.notify()

    def retarget(self, capacity) -> None:
        with self._condition:
            self._capacity = max(1, int(capacity))
            self._condition.notify_all()


_PREFETCH_QUEUE_SLACK = 3
"""Design doc Sec7: prefetch queue depth N ~= image_workers + small
constant (+2 to +4). A structural default, not a measured value -- same
category as MAX_CONCURRENT_ACTIVE_CHECKPOINTS."""
_PREFETCH_POOL_INITIAL = 4
"""Sec7 point 1: a small fixed-default prefetch pool, chosen because
diminishing returns past a handful of concurrent reads is a common
pattern across storage backends -- not a claim that 4 is optimal for any
particular job's storage."""
_PREFETCH_POOL_MAX = 16
_BLOCKED_FRACTION_GROW = 0.20
_BLOCKED_FRACTION_SHRINK = 0.02
"""Sec7's own already-resolved open item: illustrative 20%/2% bands, not
tuned against real production load."""
_COORDINATOR_TICK_SECONDS = 0.3
_POLL_TIMEOUT_SECONDS = 0.2
_REBALANCE_INITIAL_SECONDS = 30
"""How soon automatic mode first re-evaluates its thread/worker pair.

A single fixed cadence cannot serve both ends of the range. At ten
minutes -- what this used to be -- a job shorter than that never
re-evaluated at all, and every job spent its first ten minutes on
whatever the seed guessed. Checking early and then backing off costs one
sweep against a sample tile, and the sweep stops itself at the first
plateau."""

_REBALANCE_MAXIMUM_SECONDS = 600
"""Ceiling the interval backs off to once the pair stops changing."""

_REBALANCE_BACKOFF = 2.0
"""Interval growth per check that leaves the pair unchanged. A check that
does change it resets to :data:`_REBALANCE_INITIAL_SECONDS`, so the
scheduler stays attentive while conditions are still moving and goes
quiet once they settle."""
"""Sec7 Phase 4b: periodic, not continuous, re-evaluation cadence for the
joint kernel_threads/image_workers rebalance. Illustrative, like the
blocked-fraction bands above -- not tuned against real production load."""
_REBALANCE_RATE_HYSTERESIS = 0.25
"""Only commit a kernel_threads rebuild if the measured delivery rate has
moved by more than this fraction since the pair currently in effect was
chosen (design doc Sec17's resolved item) -- avoids rebuilding on noise
even when a rebalance check happens to land on a slightly different
candidate."""
_KERNEL_SWEEP_TILE_PIXELS = 1_048_576
_KERNEL_SWEEP_PLATEAU_RATIO = 0.9


def _kernel_threads_candidates(total_threads, include=()):
    """Candidate ``kernel_threads`` values for the Sec7 sweep.

    Powers of two up to ``total_threads``, plus ``total_threads`` itself
    if not already a power of two, plus every value in ``include`` (e.g.
    the currently-active ``kernel_threads``, which the caller must always
    inject -- see :func:`_map_pending_ranges`). Generated from the job's
    own thread budget, never a hardcoded absolute list: the design doc's
    own ``{1,2,4,8}`` vs. ``{1,4,8,16,24+}`` examples are per-depth
    illustrations of this same construction at ``total_threads=24``, not
    values to hardcode.

    :param int total_threads:
        The job/node's own thread budget.
    :param iterable include:
        Additional values to force into the candidate set.
    :returns:
        Sorted, deduplicated list of positive ``kernel_threads`` values,
        each at most ``total_threads``.
    :rtype: list[int]
    """
    total_threads = max(1, int(total_threads))
    candidates = {1}
    value = 1
    while value < total_threads:
        value *= 2
        candidates.add(min(value, total_threads))
    for extra in include:
        extra = int(extra)
        if 1 <= extra <= total_threads:
            candidates.add(extra)
    return sorted(candidates)


_GROUP_SIZE_CANDIDATES = (1, 2, 4, 8, 16)
"""Group sizes the automatic choice considers."""
_GROUP_ADVANCE_VOXEL_LIMIT = 1.0
"""Per-frame reciprocal-space travel, in voxels, past which grouping stops
paying.

Frame grouping merges frames that reach the same voxels. Measured over a
striding sweep of the reference scan, the gain is already gone by the time
a frame advances a whole voxel -- past that, consecutive frames tile
rather than overlap, and a group buffer buys nothing while still costing
concurrency. This is a statement about grid step against angular step: a
finer grid moves a job into the regime, a coarser scan moves it out."""
_GROUP_MINIMUM_CONCURRENT_CALLS = 3
"""Concurrent native calls a group size has to leave affordable.

Purely empirical, and the reason the automatic choice is not simply "the
largest group that fits". On the reference job, four frames per group with
three concurrent calls mapped at 0.866x the per-frame pipeline; eight with
two was 1.014x, and eight with one was 1.24x. Records keep falling as the
group grows, so density alone would pick the losing configuration -- what
turns over is the concurrency the group buffer and its native working set
crowd out."""


def _frames_per_group(spec):
    """Resolved group size for ``spec``: ``None`` means one frame."""
    value = getattr(spec, "frames_per_group", None)
    return 1 if value is None else max(1, int(value))


_GROUP_PIPELINE_MINIMUM_READAHEAD = 1
"""Groups prepared *beyond* the ones being mapped.

At least one, or there is no double-buffering at all: every concurrent
call would stop dead at its group boundary while the next group is loaded
and corrected, which is the synchronisation point that makes grouping
lose."""
_GROUP_PIPELINE_MAXIMUM_READAHEAD = 3
"""More than a few groups of read-ahead buys nothing and costs the largest
resident term in the pipeline, one full group buffer each."""


class _GroupAssembly:
    """One frame group being loaded and corrected, frame by frame.

    Prepare workers cooperate on a single group rather than each taking a
    whole group of their own. A group buffer is ``F`` full-detector
    corrected frames -- the largest resident term in the pipeline -- so
    per-worker groups would multiply it by the pool size. Frame
    granularity keeps the number of live buffers down to the pipeline
    depth while still spreading a group's ``F`` loads and corrections
    across every worker in the pool.
    """

    __slots__ = ("frames", "payloads", "corrected", "_remaining", "_lock", "failed")

    def __init__(self, frames):
        self.frames = list(frames)
        self.payloads = [None] * len(self.frames)
        self.corrected = [None] * len(self.frames)
        self._remaining = len(self.frames)
        self._lock = threading.Lock()
        self.failed = False

    def complete_one(self, *, failed=False) -> bool:
        """Record one frame finished. ``True`` when the group is whole."""
        with self._lock:
            if failed:
                self.failed = True
            self._remaining -= 1
            return self._remaining <= 0


def _group_pipeline_layout(spec, tiles, memory_bytes, prepare_workers):
    """Concurrent group calls, native threads each, and read-ahead depth.

    The grouped scheduler spends its memory differently from the
    per-frame one. There, every concurrent frame carried its own native
    working set *and* its own record batch, so the budget bought
    ``image_workers`` directly. Here the budget is split three ways:
    concurrent group calls, the read-ahead buffer in front of them, and
    the prepare pool's own transients.

    One concurrent call with every thread was the design's first guess
    and it measured badly -- three calls of eight threads beat one call
    of twenty-four by 1.38x on real data at four frames per group. A
    group call has thousands of bricks, so the thread budget is not
    limited by work available to hand out; it is limited by what the call
    does around the parallel region. Concurrency therefore still has to
    be bought, exactly as it did before grouping, and this decides how
    much of it the budget affords.

    :param int prepare_workers:
        Prepare pool size. Each worker holds one raw image and one
        corrected frame while it works, independently of the group
        buffers those frames land in.
    :returns:
        ``(compute_workers, kernel_threads, pipeline_depth)``.
    :rtype: tuple[int, int, int]
    """
    tiles = list(tiles)
    tile_pixels = [
        (row_stop - row_start) * (column_stop - column_start)
        for row_start, row_stop, column_start, column_stop in tiles
    ]
    if not tile_pixels:
        return (
            1,
            max(1, int(spec.threads)),
            1 + _GROUP_PIPELINE_MINIMUM_READAHEAD,
        )
    detector_pixels = sum(tile_pixels)
    frames_per_group = _frames_per_group(spec)
    native_bytes_per_pixel = (
        128 + 2 * _RESERVED_RECORDS_PER_PIXEL * _CHECKPOINT_BYTES_PER_ROW
        if spec.max_depth > 0
        else 128 + 2 * _CHECKPOINT_BYTES_PER_ROW
    )
    # Per concurrent group call: its native working set over the largest
    # tile times the group size, plus the merged record batch it builds.
    call_bytes = (
        max(tile_pixels) * frames_per_group * native_bytes_per_pixel
        + detector_pixels * _FRAME_RECORD_BYTES_PER_PIXEL
    )
    group_buffer_bytes = max(
        1, frames_per_group * detector_pixels * _PYTHON_CORRECTION_BYTES_PER_PIXEL
    )
    prepare_bytes = max(1, int(prepare_workers)) * detector_pixels * (
        _RAW_IMAGE_BYTES_PER_PIXEL + _PYTHON_CORRECTION_BYTES_PER_PIXEL
    )
    spare = max(0, int(memory_bytes) - prepare_bytes)
    # Every concurrent call needs a group buffer to read from, and the
    # read-ahead sits on top of that -- so a worker costs both terms.
    per_worker_bytes = call_bytes + group_buffer_bytes
    compute_workers = max(
        1,
        min(
            max(1, int(spec.threads)),
            spare // max(1, per_worker_bytes),
        ),
    )
    kernel_threads = max(1, int(spec.threads) // compute_workers)
    remaining = spare - compute_workers * per_worker_bytes
    pipeline_depth = compute_workers + max(
        _GROUP_PIPELINE_MINIMUM_READAHEAD,
        min(
            _GROUP_PIPELINE_MAXIMUM_READAHEAD,
            remaining // group_buffer_bytes,
        ),
    )
    return compute_workers, kernel_threads, pipeline_depth


def _warn_if_nothing_was_routed(router, routed_before, total_images, progress):
    """Say so when a whole mapping run produced no records at all.

    Mapping every frame of a range and emitting nothing is possible
    legitimately -- a grid can cover a region this slice of the scan never
    reaches -- so this warns rather than raising. But it is also what an
    intermittent failure observed in 2026-08 looks like from the outside:
    roughly one run in twenty mapped nothing, wrote a checkpoint claiming
    its full frame count with zero rows, and exited successfully, in about
    a fifth of the usual time. Nothing downstream could tell that from a
    fast, empty-but-correct run, and on resume the empty part counts as
    done.

    Reported through ``progress`` as well as :mod:`warnings`, because a
    warning alone is invisible in the GUI, which is where an empty
    reconstruction is most expensive to discover late.
    """
    routed = getattr(router, "routed_records", None)
    if routed is None or total_images <= 0 or routed > routed_before:
        return
    message = (
        f"Mapped {total_images} frames and produced no records at all. "
        "This is expected only if the grids genuinely cover no part of "
        "the reciprocal space these frames reach; otherwise the frames "
        "were masked away or their geometry placed every sample outside "
        "every grid, and the resulting checkpoints will be empty."
    )
    warnings.warn(message, RuntimeWarning, stacklevel=2)
    if progress is not None:
        progress(total_images, total_images, message)


def _angles_advance_monotonically(bounds, frame_indices):
    """Whether frames adjacent in index are also adjacent in angle.

    Frame grouping merges several images inside one native call on the
    premise that consecutive frames sit next to each other in reciprocal
    space. An interlaced scan (``orgui/backend/interlacedScanLoader.py``)
    breaks that premise: frames adjacent in file order can be half a
    rotation apart, and grouping them would collapse nothing while still
    costing the group buffer. Rather than assume the scan is sequential,
    check the angles the job actually resolved.

    Every axis that moves at all across the sampled frames must move in
    one direction only. An axis that never moves is ignored, so a scan
    that rotates a single motor is judged on that motor alone.

    :param bounds:
        ``(len(scan), 2, 4)`` exposure angle bounds in radians.
    :param frame_indices:
        Frame indices in the order they would be grouped.
    :returns:
        ``True`` when grouping is sound for these frames.
    :rtype: bool
    """
    if len(frame_indices) < 3:
        return True
    starts = np.asarray(bounds, dtype=np.float64)[list(frame_indices), 0]
    steps = np.diff(starts, axis=0)
    for axis in range(starts.shape[1]):
        column = steps[:, axis]
        span = float(np.ptp(starts[:, axis]))
        if span <= 1e-9:
            continue
        if not (np.all(column >= -1e-12) or np.all(column <= 1e-12)):
            return False
    return True


def _frame_advance_voxels(
    kernel, grid, corner_rays, bounds, frame_pairs, *, samples_per_axis=8
):
    """Median distance a pixel travels between adjacent frames, in voxels.

    The quantity finding A rests on, measured on the job rather than
    inferred from its angular step: two frames merge usefully only while
    they land within about a voxel of each other. Geometry and the grid
    decide it, so this needs no image data and does no I/O -- it asks the
    kernel where a pixel lands on one frame and on the next, and takes the
    difference in units of the grid step.

    A median rather than a mean: the distribution across a detector is
    wide (measured p10 0.22, p90 1.48 voxels on the reference job), and
    the tails are the corners rather than the bulk of the frame.

    :param kernel:
        Constructed ``ReconstructionKernel`` for ``grid``.
    :param corner_rays:
        Corner rays for one detector tile, ``(rows + 1, columns + 1, 3)``.
        Pixels are sampled inside that tile.
    :param bounds:
        ``(len(scan), 2, 4)`` exposure angle bounds in radians.
    :param frame_pairs:
        ``(first, second)`` adjacent frame indices to sample across.
    :returns:
        Median displacement in voxels; ``0.0`` when nothing was sampled.
    :rtype: float
    """
    rows = int(corner_rays.shape[0]) - 1
    columns = int(corner_rays.shape[1]) - 1
    if rows < 1 or columns < 1:
        return 0.0
    step = np.asarray(grid.step, dtype=np.float64)
    row_samples = np.unique(
        np.linspace(0, rows - 1, min(samples_per_axis, rows)).astype(int)
    )
    column_samples = np.unique(
        np.linspace(0, columns - 1, min(samples_per_axis, columns)).astype(int)
    )
    distances = []
    for first, second in frame_pairs:
        angles = [
            (
                np.ascontiguousarray(bounds[index, 0]),
                np.ascontiguousarray(bounds[index, 1]),
            )
            for index in (first, second)
        ]
        for row in row_samples:
            for column in column_samples:
                positions = [
                    np.asarray(
                        kernel.coordinate(
                            corner_rays, start, end, int(row), int(column)
                        )
                    )
                    for start, end in angles
                ]
                distances.append(
                    float(
                        np.linalg.norm((positions[1] - positions[0]) / step)
                    )
                )
    return float(np.median(distances)) if distances else 0.0


def _choose_frames_per_group(
    spec,
    tiles,
    memory_bytes,
    advance_voxels,
    prepare_workers,
    candidates=_GROUP_SIZE_CANDIDATES,
):
    """Largest group size that is both in the regime and affordable.

    Two independent gates, and the second is the one that is easy to get
    wrong. Records emitted fall monotonically as the group grows, so
    choosing on record density alone picks the largest group that fits
    memory -- which on the reference job is the configuration that maps
    *slowest*, because the group buffer and the call's native working set
    crowd out the concurrency that was paying for itself. So the rule is
    the largest group that still leaves
    ``_GROUP_MINIMUM_CONCURRENT_CALLS`` native calls affordable, not the
    largest group that fits.

    :param float advance_voxels:
        Per-frame reciprocal-space travel from
        :func:`_frame_advance_voxels`.
    :returns:
        Frames per group, at least 1.
    :rtype: int
    """
    if advance_voxels >= _GROUP_ADVANCE_VOXEL_LIMIT:
        return 1
    chosen = 1
    for candidate in sorted(int(value) for value in candidates):
        if candidate <= 1:
            continue
        workers, _threads, _depth = _group_pipeline_layout(
            replace(spec, frames_per_group=candidate),
            tiles,
            memory_bytes,
            prepare_workers,
        )
        if workers < _GROUP_MINIMUM_CONCURRENT_CALLS:
            break
        chosen = candidate
    return chosen


def _frame_groups(frame_range, router, grid_names, frames_per_group):
    """Split one scheduling range into the groups mapped in one call each.

    A group must be contiguous in frame index and must lie inside one
    planned checkpoint range for *every* grid: its frames merge inside the
    kernel and can no longer be told apart, so they have to share a
    checkpoint. Scheduling ranges are already contiguous and free of
    excluded frames (:func:`_included_ranges` splits on both), so the only
    cut this has to make beyond the group size is at a checkpoint
    boundary -- which need not align with the scheduling ranges, since the
    two are sized independently.

    :returns:
        List of frame-index lists, each at most ``frames_per_group`` long,
        together covering the range in order.
    :rtype: list[list[int]]
    """
    start, stop = frame_range
    size = max(1, int(frames_per_group))
    groups = []
    current: list[int] = []
    previous_keys = None
    for frame in range(start, stop):
        keys = tuple(
            router.checkpoint_index_for_frame(grid_name, frame)
            for grid_name in grid_names
        )
        if current and (keys != previous_keys or len(current) >= size):
            groups.append(current)
            current = []
        previous_keys = keys
        current.append(frame)
    if current:
        groups.append(current)
    return groups


def _map_frame_groups_streamed(
    spec,
    scan,
    config,
    bounds,
    detector_tiles,
    ray_arrays,
    frame_groups,
    router,
    *,
    correction_pipeline,
    kernel_memory_budget,
    compute_workers,
    kernel_threads,
    pipeline_depth,
    prepare_workers,
    maximum_prepare_workers,
    total_images,
    completed_images,
    progress,
):
    """Map frame groups with one all-threads native call at a time.

    The scheduler the per-frame pipeline uses answers a question frame
    grouping dissolves. There, one frame's kernel call could not
    profitably use every thread, so the budget was split jointly between
    ``image_workers`` and ``kernel_threads`` and rebalanced live. A group
    call has thousands of bricks and saturates the thread budget on its
    own, so there is nothing left to split: ``kernel_threads`` is pinned
    at the budget and exactly one group is mapped at a time.

    What replaces the joint balance is a simpler question -- how many
    workers it takes to keep that one call fed. Loading and correcting a
    frame is GIL-held Python; mapping is not. So correction moves out of
    the compute worker and into the prepare pool, where group N+1's
    Python work overlaps group N's native call instead of queueing behind
    it. Without that move a group is a hard synchronisation point and the
    pipeline stalls once per group, which is what made grouping lose
    end-to-end even while it cut records by 40%.

    The pool is sized the same way the reader pool always was, by the
    blocked fraction the compute worker measures: grow eagerly when it is
    starved, shrink cautiously when it is not. That signal now steers
    load *and* correction together, which is the whole per-frame Python
    cost rather than only its I/O half.

    Mutates ``router`` in place; returns nothing.

    :param int compute_workers:
        Concurrent group calls.
    :param int kernel_threads:
        Native threads per call.
    :param int pipeline_depth:
        Groups that may be resident at once
        (:func:`_group_pipeline_layout`).
    :param int prepare_workers:
        Initial prepare pool size.
    :param int maximum_prepare_workers:
        Ceiling for the pool, one when the scan cannot be read
        concurrently.
    """
    frame_correction = getattr(correction_pipeline, "correct_frame", None)
    corrects_whole_frame = callable(frame_correction)

    progress_events = SimpleQueue()
    frame_queue = SimpleQueue()
    ready_queue = SimpleQueue()
    gate = _BoundedGate(pipeline_depth)
    dispatch_done = threading.Event()
    mapped_images = completed_images

    counter_lock = threading.Lock()
    dispatched = [0]
    completed = [0]

    exception_lock = threading.Lock()
    first_exception: list[BaseException] = []

    def record_exception(exc):
        with exception_lock:
            if not first_exception:
                first_exception.append(exc)

    def should_stop():
        return bool(first_exception)

    def all_groups_retired():
        with counter_lock:
            return dispatch_done.is_set() and completed[0] >= dispatched[0]

    def dispatch_loop(cancellation):
        """Hand out one group's frames at a time, under the gate.

        Issuing every group's frames up front would let the pool spread
        itself over all of them and hold a group buffer for each. The
        gate is acquired per group and released when that group has been
        mapped, so the number of live buffers is the pipeline depth and
        nothing else.
        """
        def stop_waiting():
            return cancellation.is_set() or should_stop()

        try:
            for group in frame_groups:
                if stop_waiting():
                    break
                if not gate.acquire(
                    poll_timeout=_POLL_TIMEOUT_SECONDS,
                    should_stop=stop_waiting,
                ):
                    break
                assembly = _GroupAssembly(group)
                with counter_lock:
                    dispatched[0] += 1
                for slot in range(len(group)):
                    frame_queue.put((assembly, slot))
        except BaseException as exc:  # noqa: BLE001 -- must reach the coordinator
            record_exception(exc)
        finally:
            dispatch_done.set()

    def prepare_loop(retire, cancellation):
        """Load and correct one frame per item, into its group's slot.

        Exits only between items, never mid-frame: a claimed slot must be
        completed one way or the other, since its group can never become
        whole otherwise and the coordinator would wait forever for a
        group that no one is still working on.
        """
        try:
            while True:
                if retire.is_set() or cancellation.is_set():
                    return
                try:
                    assembly, slot = frame_queue.get(
                        timeout=_POLL_TIMEOUT_SECONDS
                    )
                except Empty:
                    if dispatch_done.is_set() and frame_queue.empty():
                        return
                    continue
                failed = False
                if should_stop():
                    # Drain without doing the work: the group still has to
                    # complete so the compute side can retire it and
                    # release its gate slot.
                    failed = True
                else:
                    frame_index = assembly.frames[slot]
                    try:
                        payload = scan.get_raw_img(frame_index)
                        if corrects_whole_frame:
                            assembly.corrected[slot] = frame_correction(
                                payload, np.asarray(payload.img), frame_index
                            )
                        else:
                            assembly.payloads[slot] = payload
                    except BaseException as exc:  # noqa: BLE001
                        record_exception(exc)
                        failed = True
                if assembly.complete_one(failed=failed):
                    ready_queue.put(assembly)
        except BaseException as exc:  # noqa: BLE001
            record_exception(exc)

    blocked_counts: dict[int, list[int]] = {}
    blocked_counts_lock = threading.Lock()

    def compute_loop(retire, cancellation):
        counts = [0, 0]  # [blocked, total], only this thread ever writes
        with blocked_counts_lock:
            blocked_counts[id(counts)] = counts
        try:
            try:
                kernels = _build_kernels(
                    spec,
                    config.ub_calculator,
                    threads=kernel_threads,
                    memory_budget_bytes=kernel_memory_budget,
                )
                while True:
                    if retire.is_set() or cancellation.is_set():
                        return
                    try:
                        assembly = ready_queue.get(
                            timeout=_POLL_TIMEOUT_SECONDS
                        )
                    except Empty:
                        counts[0] += 1
                        counts[1] += 1
                        if all_groups_retired() and ready_queue.empty():
                            return
                        continue
                    counts[1] += 1
                    try:
                        if not assembly.failed:
                            group = assembly.frames
                            _map_frame_group(
                                spec,
                                kernels,
                                ray_arrays,
                                detector_tiles,
                                correction_pipeline,
                                None if corrects_whole_frame else assembly.payloads,
                                group,
                                np.ascontiguousarray(bounds[group, 0]),
                                np.ascontiguousarray(bounds[group, 1]),
                                router,
                                corrected_frames=(
                                    assembly.corrected
                                    if corrects_whole_frame
                                    else None
                                ),
                            )
                    except BaseException as exc:  # noqa: BLE001
                        record_exception(exc)
                    finally:
                        # Drop the group buffer before releasing its slot,
                        # so the gate bounds what is actually resident
                        # rather than what has merely been claimed.
                        assembly.corrected = []
                        assembly.payloads = []
                        with counter_lock:
                            completed[0] += 1
                        gate.release()
                    if not assembly.failed:
                        progress_events.put(len(assembly.frames))
            except BaseException as exc:  # noqa: BLE001
                record_exception(exc)
        finally:
            with blocked_counts_lock:
                blocked_counts.pop(id(counts), None)

    dispatch_cancellation = threading.Event()
    dispatcher = threading.Thread(
        target=dispatch_loop,
        args=(dispatch_cancellation,),
        name="orgui-rsmap-dispatch",
        daemon=True,
    )
    dispatcher.start()
    prepare_pool = _AdjustablePool(
        prepare_loop, initial_size=prepare_workers, name="orgui-rsmap-prepare"
    )
    compute_pool = _AdjustablePool(
        compute_loop, initial_size=compute_workers, name="orgui-rsmap-compute"
    )
    previous_blocked = 0
    previous_total = 0
    try:
        while True:
            deadline = time.monotonic() + _COORDINATOR_TICK_SECONDS
            while time.monotonic() < deadline:
                try:
                    mapped_in_group = progress_events.get(
                        timeout=max(0.0, deadline - time.monotonic())
                    )
                except Empty:
                    break
                mapped_images += mapped_in_group
                if progress is not None:
                    progress(
                        mapped_images,
                        total_images + 1,
                        (
                            f"Mapping images {mapped_images}/{total_images} "
                            f"({spec.frames_per_group} frames/call, "
                            f"{compute_pool.size} concurrent calls x "
                            f"{kernel_threads} native threads, "
                            f"{prepare_pool.size} prepare workers)"
                        ),
                    )
            with blocked_counts_lock:
                blocked = sum(counts[0] for counts in blocked_counts.values())
                total = sum(counts[1] for counts in blocked_counts.values())
            window_blocked = max(0, blocked - previous_blocked)
            window_total = max(0, total - previous_total)
            previous_blocked, previous_total = blocked, total
            if window_total > 0 and not all_groups_retired() and not first_exception:
                blocked_fraction = window_blocked / window_total
                if blocked_fraction > _BLOCKED_FRACTION_GROW:
                    prepare_pool.retarget(
                        min(maximum_prepare_workers, prepare_pool.size + 2)
                    )
                elif (
                    blocked_fraction < _BLOCKED_FRACTION_SHRINK
                    and prepare_pool.size > 1
                ):
                    prepare_pool.retarget(prepare_pool.size - 1)
            prepare_pool.reap()
            compute_pool.reap()
            with blocked_counts_lock:
                still_computing = bool(blocked_counts)
            if not still_computing and (all_groups_retired() or first_exception):
                # The second condition is not redundant. A compute worker
                # that dies -- building its kernels, or on a mapping
                # failure -- leaves dispatched groups that nothing will
                # ever complete, so waiting for the counters to agree
                # would hang here forever instead of raising.
                break
        if first_exception:
            raise first_exception[0]
    finally:
        dispatch_cancellation.set()
        dispatcher.join(timeout=_POLL_TIMEOUT_SECONDS * 10)
        prepare_pool.shutdown(wait=True)
        compute_pool.shutdown(wait=True)


def _map_pending_ranges(
    spec,
    scan,
    config,
    bounds,
    tiles,
    pending_ranges,
    router,
    *,
    correction_pipeline,
    effective_memory,
    threads_per_image,
    accumulation_budget_bytes,
    total_images,
    completed_images,
    progress,
):
    """Map every range in ``pending_ranges`` against ``router``.

    A producer/consumer pipeline (design doc Sec7): a prefetch pool of
    reader threads loads images and feeds a bounded-backpressure queue; a
    pool of compute workers drains it, correcting and accumulating each
    group of frames via :func:`_map_frame_group`. The reader pool adapts
    continuously via a blocked-fraction signal (grow eagerly when compute
    is starved for images, shrink cautiously when it isn't).

    The unit of work is ``spec.frames_per_group`` consecutive frames
    mapped in one native call -- one frame per group by default, which is
    the degenerate case of the same machinery rather than a separate
    path. Groups are cut at scheduling-range and checkpoint boundaries
    (:func:`_frame_groups`), and dropped to one frame per group entirely
    when the scan's angles do not advance monotonically
    (:func:`_angles_advance_monotonically`).

    Above one frame per group the whole scheduling question changes, and
    :func:`_map_frame_groups_streamed` answers the new one instead: a
    group call saturates the thread budget by itself, so there is no pair
    left to balance, and what matters is keeping that one call fed. See
    its docstring.

    At one frame per group the joint balance still applies, and is
    unchanged. ``threads_per_image=None`` (Sec7 Phase 4b, "automatic")
    starts from an I/O-optimistic seed (``kernel_threads=1``, all budget
    as ``image_workers``) and periodically rebalances
    ``kernel_threads``/``image_workers`` live against a measured
    frame-delivery rate (every ``_REBALANCE_INITIAL_SECONDS``, backing off
    while nothing changes), via a wall-clock ``_kernel_threads_sweep`` and
    the design doc's joint-balancing rule. A concrete
    ``threads_per_image`` int keeps the static,
    ``_frame_parallelism``-derived pair fixed for the whole run (today's
    Phase 4a behavior, unchanged) -- the explicit override escape hatch.

    Extracted from :func:`run_job`'s mapping loop so
    :func:`run_cluster_map_task` can reuse it unchanged, scoped to one
    node's own ``pending_ranges`` instead of the whole job's. Mutates
    ``router`` (and the checkpoint files it writes) in place; returns
    nothing. Not resume-aware itself -- callers decide ``pending_ranges``
    (see :func:`_range_is_resumed`).
    """
    if not pending_ranges:
        return

    automatic = threads_per_image is None
    seed_threads_per_image = 1 if automatic else threads_per_image

    detector_tiles, bounds = _validate_mapping_setup(
        scan, config.detector, tiles, bounds, correction_pipeline
    )

    frame_indices = [
        frame_index
        for start, stop in pending_ranges
        for frame_index in range(start, stop)
    ]
    # Grouping rests on frames adjacent in index being adjacent in angle.
    # Where they are not -- an interlaced scan -- the honest answer is one
    # frame per call, whatever else the job asked for. Cheap, and it gates
    # both the requested and the measured group size, so it is settled
    # here; the rest of the choice needs the ray arrays and the memory
    # split and waits for them.
    grouping_is_sound = _angles_advance_monotonically(bounds, frame_indices)

    ray_cache_bytes = sum(
        (tile[1] - tile[0] + 1)
        * (tile[3] - tile[2] + 1)
        * 3
        * np.dtype(np.float64).itemsize
        for tile in tiles
    )
    cache_detector_rays = ray_cache_bytes <= effective_memory // 8
    ray_cache = {}
    if cache_detector_rays:
        for tile in tiles:
            ray_cache[tuple(tile)] = _detector_corner_rays(config.detector, tile)
    ray_arrays = _tile_ray_arrays(
        config.detector, detector_tiles, ray_cache if cache_detector_rays else None
    )

    scheduler_memory = max(
        1024**2,
        effective_memory - (ray_cache_bytes if cache_detector_rays else 0),
    )

    max_readers = (
        _PREFETCH_POOL_MAX
        if getattr(scan, "supports_concurrent_read", True)
        else 1
    )
    seed_prepare_workers = min(_PREFETCH_POOL_INITIAL, max_readers)
    requested_frames_per_group = getattr(spec, "frames_per_group", None)
    if not grouping_is_sound:
        frames_per_group = 1
    elif requested_frames_per_group is not None:
        frames_per_group = max(1, int(requested_frames_per_group))
    elif len(frame_indices) < 2:
        frames_per_group = 1
    else:
        # Sample the travel per frame at a few places in the scan, not
        # one: a scan with a non-uniform angular step is the mild version
        # of the interlaced problem, and one pair at the start would miss
        # it entirely.
        last_pair = len(frame_indices) - 2
        probe_pairs = [
            (frame_indices[position], frame_indices[position + 1])
            for position in sorted(
                {0, min(last_pair, len(frame_indices) // 2), last_pair}
            )
            # Skip a pair that straddles an excluded frame: its travel is
            # two frames' worth and would read as out of the regime.
            if frame_indices[position] + 1 == frame_indices[position + 1]
        ]
        frames_per_group = (
            _choose_frames_per_group(
                spec,
                tiles,
                scheduler_memory,
                _frame_advance_voxels(
                    _kernel_for_grid(
                        spec, spec.grids[0], config.ub_calculator, threads=1
                    ),
                    spec.grids[0],
                    ray_arrays[detector_tiles[0]],
                    bounds,
                    probe_pairs,
                ),
                seed_prepare_workers,
            )
            if probe_pairs
            else 1
        )
    spec = replace(spec, frames_per_group=frames_per_group)

    def layout_for(threads_per_image_seed):
        """Worker count and native threads this seed would run with."""
        workers = []
        natives = []
        for frame_range in pending_ranges:
            image_limit, native_threads, _worker_memory, _accumulation_limit = (
                _frame_parallelism(
                    spec,
                    tiles,
                    scheduler_memory,
                    frames_per_task=frame_range[1] - frame_range[0],
                    threads_per_image=threads_per_image_seed,
                    accumulation_budget_bytes=accumulation_budget_bytes,
                )
            )
            workers.append(image_limit)
            natives.append(native_threads)
        return (
            min(len(pending_ranges), min(workers)) if workers else 1,
            min(natives) if natives else 1,
        )

    if automatic:
        # Seed from how many frames can actually be in flight, not from
        # one native thread each. The memory budget caps concurrent
        # frames, and giving each a single thread then leaves most of the
        # thread budget idle -- five of twenty-four at balanced accuracy
        # on a real job, measured 1.85x slower than using them all. The
        # live rebalance would eventually find this, but not before its
        # first interval, which outlasts many jobs.
        affordable_workers, _seed_threads = layout_for(1)
        seed_threads_per_image = max(
            1, math.ceil(spec.threads / max(1, affordable_workers))
        )
    image_workers, seed_kernel_threads = layout_for(seed_threads_per_image)
    # A mutable box: read by compute_loop at each worker's kernel-build
    # time (not closed over as a plain value), so an automatic-mode
    # kernel_threads change only affects newly-spawned compute workers --
    # see the rebalance block in the coordinator loop below.
    current_kernel_threads = [seed_kernel_threads]
    # The native budget is a per-call guard. The scheduler independently
    # bounds the aggregate worker working set above.
    kernel_memory_budget = effective_memory

    grid_names = [grid.grid_name for grid in spec.grids]
    frame_groups = [
        group
        for frame_range in pending_ranges
        for group in _frame_groups(
            frame_range, router, grid_names, frames_per_group
        )
    ]

    routed_before = getattr(router, "routed_records", 0)

    if frames_per_group > 1:
        group_workers, group_threads, group_depth = _group_pipeline_layout(
            spec, tiles, scheduler_memory, seed_prepare_workers
        )
        _map_frame_groups_streamed(
            spec,
            scan,
            config,
            bounds,
            detector_tiles,
            ray_arrays,
            frame_groups,
            router,
            correction_pipeline=correction_pipeline,
            kernel_memory_budget=effective_memory,
            compute_workers=group_workers,
            kernel_threads=group_threads,
            pipeline_depth=group_depth,
            prepare_workers=seed_prepare_workers,
            maximum_prepare_workers=max_readers,
            total_images=total_images,
            completed_images=completed_images,
            progress=progress,
        )
        _warn_if_nothing_was_routed(
            router, routed_before, total_images, progress
        )
        return

    if automatic:
        # A single representative tile/frame is enough for the sweep: it
        # measures kernel_threads scaling behavior (thread-count
        # plateau shape), not a scientifically exact per-job estimate --
        # unlike the checkpoint accumulation itself, this only informs a
        # scheduling decision. A trivial all-included mask is used
        # (real per-frame masks only exist after correction, which is
        # not needed here).
        sweep_tile = detector_tiles[0]
        sweep_mask = np.zeros(
            (sweep_tile[1] - sweep_tile[0], sweep_tile[3] - sweep_tile[2]),
            dtype=bool,
        )
        sweep_rays = ray_arrays[sweep_tile]
        sweep_angles_start = np.ascontiguousarray(bounds[frame_indices[0], 0])
        sweep_angles_end = np.ascontiguousarray(bounds[frame_indices[0], 1])
        # Every tile of a frame is mapped, so a frame's kernel work is the
        # whole detector's pixels -- what the sweep's single sample tile
        # has to be scaled up to.
        frame_pixels = sum(
            (row_stop - row_start) * (column_stop - column_start)
            for row_start, row_stop, column_start, column_stop in detector_tiles
        )

    reader_pool_size = min(_PREFETCH_POOL_INITIAL, max_readers)

    progress_events = SimpleQueue()
    ready_queue = SimpleQueue()
    gate = _BoundedGate(image_workers + _PREFETCH_QUEUE_SLACK)
    mapped_images = completed_images

    work_lock = threading.Lock()
    # The unit of work is a group of frames mapped in one native call --
    # one frame per group when frames_per_group is 1. The gate therefore
    # bounds in-flight *groups* rather than frames: a reader loads a whole
    # group before queueing it, so a per-frame gate would let several
    # readers each hold a partial group and none of them make progress.
    # _frame_parallelism reserves the queue slack in whole groups to match.
    work_iterator = iter(frame_groups)
    remaining = [len(frame_groups)]
    remaining_lock = threading.Lock()
    readers_done = threading.Event()
    # Successfully-read frames since the last rebalance check -- the
    # reader pool's own delivery rate, not end-to-end throughput (which
    # would be capped by compute when compute is the bottleneck, making
    # the rebalance rule unable to ever discover more image_workers would
    # help; see design doc Sec7 Phase 4b notes).
    delivered_in_window = [0]

    exception_lock = threading.Lock()
    first_exception: list[BaseException] = []

    def record_exception(exc):
        with exception_lock:
            if not first_exception:
                first_exception.append(exc)

    def should_stop():
        return bool(first_exception)

    def mark_group_delivered(group_frames, *, succeeded=True):
        """One group off the reader pool's plate.

        ``remaining`` counts groups (the work unit), but
        ``delivered_in_window`` stays in frames: it feeds the rebalance's
        Little's-law rule against a sweep that reports a per-frame time,
        and mixing the two units there would silently misprice every
        candidate thread count.
        """
        with remaining_lock:
            remaining[0] -= 1
            if succeeded:
                delivered_in_window[0] += group_frames
            if remaining[0] <= 0:
                readers_done.set()

    def reader_loop(retire, pool_cancellation):
        """A personal ``retire`` (reader-pool shrink) only stops this
        worker from claiming its *next* frame from ``work_iterator`` --
        once a frame index is claimed, it must be delivered (successfully
        or as a recorded failure) before this worker exits, matching
        ``compute_loop``'s own in-flight-completes-first contract.
        ``work_iterator`` is a one-shot generator: a claimed-then-abandoned
        frame index can never be reclaimed by another reader, which would
        permanently stall ``remaining`` above zero and hang the
        coordinator waiting for a ``readers_done`` that can never fire.
        """
        def local_should_stop():
            return retire.is_set() or pool_cancellation.is_set() or should_stop()

        # The outer try/except is a safety net independent of the
        # get_raw_img-specific one below: a bug anywhere in this loop's
        # own logic (not just in get_raw_img) must still reach the
        # coordinator via record_exception, not die silently -- Python's
        # default behavior for an uncaught exception in a plain
        # threading.Thread is to print a traceback and let the thread
        # exit with no other observable effect, which would otherwise
        # hang the coordinator waiting for completion signals that will
        # never come.
        def abandon_should_stop():
            # Deliberately excludes this reader's own `retire` -- once a
            # frame has been claimed from work_iterator below, abandoning
            # the gate wait just because *this* reader is being retired
            # (e.g. a live reader-pool shrink) would drop that frame on
            # the floor forever: nothing else can reclaim it from a
            # one-shot iterator, so remaining[0] would never reach zero
            # and readers_done would never fire, hanging the coordinator.
            # A hard pool_cancellation (final teardown) or a job-wide
            # should_stop() (some other failure) still abort it -- both
            # are handled below by explicitly marking it delivered
            # (failed), keeping the bookkeeping consistent either way.
            return pool_cancellation.is_set() or should_stop()

        try:
            while not local_should_stop():
                with work_lock:
                    group = next(work_iterator, None)
                if group is None:
                    return
                if not gate.acquire(
                    poll_timeout=_POLL_TIMEOUT_SECONDS,
                    should_stop=abandon_should_stop,
                ):
                    mark_group_delivered(len(group), succeeded=False)
                    return
                try:
                    image_payloads = [
                        scan.get_raw_img(frame_index) for frame_index in group
                    ]
                except BaseException as exc:  # noqa: BLE001
                    gate.release()
                    record_exception(exc)
                    mark_group_delivered(len(group), succeeded=False)
                    return
                ready_queue.put((group, image_payloads))
                mark_group_delivered(len(group))
        except BaseException as exc:  # noqa: BLE001 -- must reach the coordinator
            record_exception(exc)

    blocked_counts: dict[int, list[int]] = {}
    blocked_counts_lock = threading.Lock()

    def compute_loop(retire, pool_cancellation):
        """A personal ``retire`` (this specific worker being asked to step
        down, e.g. a future live pool shrink) or ``pool_cancellation``
        exits immediately -- other workers remain to drain the queue. A
        job-wide ``should_stop()`` (an exception anywhere) does *not*
        exit immediately: it stops accepting the possibility of new work
        arriving from readers (who do stop immediately), but drains
        whatever is already sitting in the ready queue first, so an
        already-loaded frame's real, already-paid-for I/O is never
        silently discarded -- resumability is about not repeating
        finished work, not just about not corrupting it.
        """
        counts = [0, 0]  # [blocked, total], only this thread ever writes
        with blocked_counts_lock:
            blocked_counts[id(counts)] = counts
        try:
            # Same catch-all rationale as reader_loop: any bug in this
            # loop (including kernel construction) must reach the
            # coordinator via record_exception, never die silently.
            try:
                kernels = _build_kernels(
                    spec,
                    config.ub_calculator,
                    threads=current_kernel_threads[0],
                    memory_budget_bytes=kernel_memory_budget,
                )
                while True:
                    if retire.is_set() or pool_cancellation.is_set():
                        return
                    try:
                        group, image_payloads = ready_queue.get(
                            timeout=_POLL_TIMEOUT_SECONDS
                        )
                    except Empty:
                        counts[0] += 1
                        counts[1] += 1
                        if (
                            readers_done.is_set() or should_stop()
                        ) and ready_queue.empty():
                            return
                        continue
                    counts[1] += 1
                    gate.release()
                    try:
                        _map_frame_group(
                            spec,
                            kernels,
                            ray_arrays,
                            detector_tiles,
                            correction_pipeline,
                            image_payloads,
                            group,
                            np.ascontiguousarray(bounds[group, 0]),
                            np.ascontiguousarray(bounds[group, 1]),
                            router,
                        )
                    except BaseException as exc:  # noqa: BLE001
                        record_exception(exc)
                        return
                    progress_events.put(len(group))
            except BaseException as exc:  # noqa: BLE001
                record_exception(exc)
        finally:
            with blocked_counts_lock:
                blocked_counts.pop(id(counts), None)

    reader_pool = _AdjustablePool(
        reader_loop, initial_size=reader_pool_size, name="orgui-rsmap-reader"
    )
    compute_pool = _AdjustablePool(
        compute_loop, initial_size=image_workers, name="orgui-rsmap-compute"
    )
    previous_blocked = 0
    previous_total = 0
    last_rebalance_monotonic = time.monotonic()
    rate_at_last_rebalance = 0.0
    rebalance_interval = _REBALANCE_INITIAL_SECONDS
    try:
        while True:
            deadline = time.monotonic() + _COORDINATOR_TICK_SECONDS
            while time.monotonic() < deadline:
                try:
                    mapped_in_group = progress_events.get(
                        timeout=max(0.0, deadline - time.monotonic())
                    )
                except Empty:
                    break
                mapped_images += mapped_in_group
                if progress is not None:
                    grouping = (
                        f" x {frames_per_group} frames"
                        if frames_per_group > 1
                        else ""
                    )
                    progress(
                        mapped_images,
                        total_images + 1,
                        (
                            f"Mapping images {mapped_images}/{total_images} "
                            f"({compute_pool.size} image workers{grouping}, "
                            f"{current_kernel_threads[0]} native threads/image, "
                            f"{reader_pool.size} prefetch readers)"
                        ),
                    )
            with blocked_counts_lock:
                blocked = sum(counts[0] for counts in blocked_counts.values())
                total = sum(counts[1] for counts in blocked_counts.values())
            window_blocked = max(0, blocked - previous_blocked)
            window_total = max(0, total - previous_total)
            previous_blocked, previous_total = blocked, total
            if window_total > 0 and not readers_done.is_set() and not first_exception:
                blocked_fraction = window_blocked / window_total
                if blocked_fraction > _BLOCKED_FRACTION_GROW:
                    reader_pool.retarget(min(max_readers, reader_pool.size + 2))
                elif (
                    blocked_fraction < _BLOCKED_FRACTION_SHRINK
                    and reader_pool.size > 1
                ):
                    reader_pool.retarget(reader_pool.size - 1)
            if automatic and not readers_done.is_set() and not first_exception:
                now = time.monotonic()
                elapsed = now - last_rebalance_monotonic
                if elapsed >= rebalance_interval:
                    with remaining_lock:
                        delivered = delivered_in_window[0]
                        delivered_in_window[0] = 0
                    last_rebalance_monotonic = now
                    rate = delivered / elapsed if elapsed > 0 else 0.0
                    candidates = _kernel_threads_candidates(
                        spec.threads, include=(current_kernel_threads[0],)
                    )
                    sweep = _kernel_threads_sweep(
                        spec,
                        spec.grids[0],
                        config.ub_calculator,
                        sweep_mask,
                        sweep_rays,
                        sweep_angles_start,
                        sweep_angles_end,
                        candidates=candidates,
                        tile_pixels=_KERNEL_SWEEP_TILE_PIXELS,
                        frame_pixels=frame_pixels,
                        memory_budget_bytes=kernel_memory_budget,
                        plateau_ratio=_KERNEL_SWEEP_PLATEAU_RATIO,
                    )
                    feasible = []
                    for candidate_threads, per_frame_time in sweep.items():
                        # Little's law: to consume frames arriving at
                        # `rate` while each takes `per_frame_time` to map,
                        # that many must be in flight at once. The sweep
                        # reports a whole frame's time, so this is a
                        # worker count rather than a fraction of one --
                        # measured against a sample tile alone it came out
                        # low by the ratio of a frame to that tile, six
                        # times on a 6.2-megapixel detector, which made
                        # every candidate look affordable.
                        needed = (
                            max(1, math.ceil(rate * per_frame_time))
                            if rate > 0
                            else 1
                        )
                        ceiling, _nt, _pm, _acc = _frame_parallelism(
                            spec,
                            tiles,
                            scheduler_memory,
                            threads_per_image=candidate_threads,
                            accumulation_budget_bytes=accumulation_budget_bytes,
                        )
                        if needed <= ceiling:
                            feasible.append(
                                (candidate_threads, min(needed, ceiling))
                            )
                    if feasible:
                        new_kernel_threads, new_image_workers = max(
                            feasible, key=lambda pair: pair[0]
                        )
                        rate_ref = max(rate_at_last_rebalance, 1e-9)
                        rate_moved = (
                            abs(rate - rate_at_last_rebalance) / rate_ref
                            > _REBALANCE_RATE_HYSTERESIS
                        )
                        if (
                            new_kernel_threads != current_kernel_threads[0]
                            and rate_moved
                        ):
                            # kernel_threads change: each compute worker's
                            # kernel is built once at worker-start, so this
                            # needs a full generation swap, not a resize.
                            # ready_queue/gate are untouched -- readers and
                            # already-queued frames are unaffected, only
                            # which pool drains the queue changes.
                            current_kernel_threads[0] = new_kernel_threads
                            gate.retarget(
                                new_image_workers + _PREFETCH_QUEUE_SLACK
                            )
                            old_compute_pool = compute_pool
                            compute_pool = _AdjustablePool(
                                compute_loop,
                                initial_size=new_image_workers,
                                name="orgui-rsmap-compute",
                            )
                            # Blocks until every straggler on the retired
                            # generation finishes its in-flight
                            # _map_frame_group call -- a deliberate, rare
                            # stall (this whole block runs at most once per
                            # rebalance interval), never abandons
                            # in-flight work.
                            old_compute_pool.shutdown(wait=True)
                            rate_at_last_rebalance = rate
                            rebalance_interval = _REBALANCE_INITIAL_SECONDS
                        elif new_image_workers != compute_pool.size:
                            gate.retarget(
                                new_image_workers + _PREFETCH_QUEUE_SLACK
                            )
                            compute_pool.retarget(new_image_workers)
                            rate_at_last_rebalance = rate
                            rebalance_interval = _REBALANCE_INITIAL_SECONDS
                        else:
                            # Nothing moved, so look again later rather
                            # than paying for a sweep at the same cadence
                            # for the rest of a long job.
                            rebalance_interval = min(
                                _REBALANCE_MAXIMUM_SECONDS,
                                rebalance_interval * _REBALANCE_BACKOFF,
                            )
            reader_pool.reap()
            compute_pool.reap()
            # Completion: every frame delivered or errored, and every
            # compute worker has drained the queue and exited on its own.
            if readers_done.is_set():
                with blocked_counts_lock:
                    still_computing = bool(blocked_counts)
                if not still_computing:
                    break
        if first_exception:
            raise first_exception[0]
    finally:
        reader_pool.shutdown(wait=True)
        compute_pool.shutdown(wait=True)
    _warn_if_nothing_was_routed(router, routed_before, total_images, progress)


def run_cluster_map_task(
    path,
    task_index,
    *,
    total_tasks,
    cpus=1,
    memory_bytes=1024**3,
    progress=None,
):
    """Execute one cluster map-array task.

    Runs the full single-node mapping pipeline (design doc Sec13) against
    this node's own disjoint frame slice, computed from
    ``task_index``/``total_tasks`` alone -- never read from a scheduler
    environment variable, since most schedulers besides Slurm never
    expose a reliable total-array-size variable to a running task (see
    ``doc/design/reciprocal_space_scratch_architecture.md`` Sec13). Never
    writes the shared job JSON: array elements must not mutate shared
    state, so they may run fully concurrently.

    :param int task_index:
        This node's 0-based position in the array (already normalized
        from the scheduler's own per-task index -- SGE's 1-based
        ``SGE_TASK_ID``, Slurm's 0-based ``SLURM_ARRAY_TASK_ID`` -- by the
        generated script).
    :param int total_tasks:
        Total array size, supplied explicitly (baked into the generated
        script by
        :func:`orgui.reconstruction_cluster.generate_cluster_scripts`),
        never inferred from scheduler environment variables.
    :param int cpus, memory_bytes:
        This array task's own scheduler resource allocation -- can differ
        from the job's captured
        ``runtime_threads``/``runtime_memory_bytes``.
    :returns:
        Status dict: ``task_index``, ``total_tasks``, ``frames`` (owned
        by this node), ``mapped`` (newly loaded this run), ``reused``
        (already-resumed frames skipped).
    :rtype: dict
    """
    path = Path(path).absolute()
    job = read_job(path)
    verify_job(job)
    scan = job.scan
    config = job.config_data
    assets = _load_assets(job)
    mask = assets.get("mask")

    plan = _node_checkpoint_plan(
        job, scan, config, mask, total_tasks=total_tasks, task_index=task_index
    )
    if plan is None:
        return {
            "status": "complete",
            "task_index": int(task_index),
            "total_tasks": int(total_tasks),
            "frames": 0,
            "mapped": 0,
            "reused": 0,
        }

    node_excluded = set(plan["excluded"])
    node_dir = _node_checkpoint_dir(job, task_index)
    boundaries = {
        grid_name: [tuple(item) for item in items]
        for grid_name, items in plan["boundaries"].items()
    }
    execution_spec = _execution_spec(job, threads=cpus, memory_bytes=memory_bytes)
    effective_memory = execution_spec.memory_budget_bytes
    number_of_grids = max(1, len(execution_spec.grids))
    active_budget_bytes, effective_memory = split_memory_budget(
        effective_memory, number_of_grids
    )
    resumed, _existing_files = _discover_checkpoint_state(
        node_dir, boundaries, job.digest, cleanup_stale=True
    )
    router = _CheckpointRouter(
        boundaries,
        spec_digest=job.digest,
        checkpoint_dir=node_dir,
        active_budget_bytes=active_budget_bytes,
        resumed=resumed,
    )

    provenance = _base_provenance(job, config)
    correct = _correction_pipeline(config, scan, assets, provenance)

    ranges, tiles = _execution_layout(
        job, scan, config, extra_excluded_frames=node_excluded
    )
    bounds = scan.exposure_angle_bounds(config, fallback=job.angle_fallback)
    pending_ranges = [
        frame_range
        for frame_range in ranges
        if not _range_is_resumed(router, boundaries, resumed, frame_range)
    ]
    total_frames = sum(stop - start for start, stop in ranges)
    reused_frames = total_frames - sum(
        stop - start for start, stop in pending_ranges
    )

    _map_pending_ranges(
        execution_spec,
        scan,
        config,
        bounds,
        tiles,
        pending_ranges,
        router,
        correction_pipeline=correct,
        effective_memory=effective_memory,
        threads_per_image=job.threads_per_image,
        accumulation_budget_bytes=job.accumulation_budget_bytes,
        total_images=total_frames,
        completed_images=reused_frames,
        progress=progress,
    )

    resumed_after, _files = _discover_checkpoint_state(
        node_dir, boundaries, job.digest, cleanup_stale=False
    )
    expected = {
        (grid_name, index)
        for grid_name, ranges_for_grid in boundaries.items()
        for index in range(len(ranges_for_grid))
    }
    missing = expected - resumed_after
    if missing:
        raise RuntimeError(
            f"Node {task_index} mapping did not produce {len(missing)} "
            f"expected checkpoint(s): {sorted(missing)[:10]}"
        )

    _atomic_json(node_dir / "provenance.json", provenance)

    return {
        "status": "complete",
        "task_index": int(task_index),
        "total_tasks": int(total_tasks),
        "frames": total_frames,
        "mapped": total_frames - reused_frames,
        "reused": reused_frames,
    }


def run_cluster_finalize(
    path,
    *,
    total_tasks,
    cpus=1,
    memory_bytes=1024**3,
    progress=None,
):
    """Verify every cluster node's checkpoints, then finalize their job.

    Reads every node's persisted checkpoint-plan sidecar
    (:func:`_node_checkpoint_plan`) and provenance, verifies each node's
    planned checkpoints are all present and digest-matching
    (:func:`_discover_checkpoint_state`, the same mechanism as the
    single-node path -- design doc Sec11), then merges every node's
    checkpoint files into one :func:`_finalize_reconstruction` call: the
    "more readers, same tree-merge" extension from design doc Sec13 --
    no change needed inside ``reconstruction.py`` for cluster support.

    :param int total_tasks:
        Must match the array size the map tasks actually ran with.
    :param int cpus:
        Accepted for CLI/script-generation symmetry with the map tasks;
        not currently used algorithmically -- finalize has no internal
        worker pool to size, matching the single-node finalize path's own
        single-threaded per-grid chunk loop.
    :param int memory_bytes:
        Accepted for CLI/script-generation symmetry; finalize already
        bounds memory per checkpoint file via
        ``_CheckpointRangeReader`` regardless of this value, matching the
        single-node path, so it is currently informational only.
    """
    path = Path(path).absolute()
    job = read_job(path)
    if job.status == "complete":
        return job_status(path)
    verify_job(job)
    if int(cpus) < 1:
        raise ValueError("cpus must be at least one")
    if int(memory_bytes) < 1:
        raise ValueError("memory_bytes must be positive")
    scan = job.scan
    config = job.config_data
    assets = _load_assets(job)
    mask = assets.get("mask")

    checkpoint_files: dict[str, list] = {}
    provenance = _base_provenance(job, config)
    incomplete = []
    for task_index in range(int(total_tasks)):
        plan = _node_checkpoint_plan(
            job,
            scan,
            config,
            mask,
            total_tasks=total_tasks,
            task_index=task_index,
        )
        if plan is None:
            continue
        node_dir = _node_checkpoint_dir(job, task_index)
        boundaries = {
            grid_name: [tuple(item) for item in items]
            for grid_name, items in plan["boundaries"].items()
        }
        resumed, files = _discover_checkpoint_state(
            node_dir, boundaries, job.digest
        )
        expected = {
            (grid_name, index)
            for grid_name, ranges_for_grid in boundaries.items()
            for index in range(len(ranges_for_grid))
        }
        if expected - resumed:
            incomplete.append(task_index)
            continue
        for grid_name, paths in files.items():
            checkpoint_files.setdefault(grid_name, []).extend(paths)
        provenance_path = node_dir / "provenance.json"
        if provenance_path.exists():
            _merge_provenance(
                provenance,
                json.loads(provenance_path.read_text(encoding="utf-8")),
            )

    if incomplete:
        raise RuntimeError(
            f"Cannot finalize: {len(incomplete)} of {total_tasks} cluster "
            f"map tasks are incomplete or missing: {incomplete[:20]}"
        )

    job.correction_provenance = provenance
    spec = job.internal_spec()
    result = _finalize_reconstruction(
        spec,
        checkpoint_files,
        job.output_path,
        provenance={
            **provenance,
            "scan_reference": job.scan_reference,
        },
        config=config,
        chunk_progress=(
            None
            if progress is None
            else lambda written, count: progress(
                written, count, f"Finalizing HDF5 chunk {written}/{count}"
            )
        ),
    )
    job.output_sha256 = result["sha256"]
    job.status = "complete"
    for task_index in range(int(total_tasks)):
        target = _node_checkpoint_dir(job, task_index)
        try:
            if target.is_dir():
                shutil.rmtree(target)
        except OSError as error:
            job.cleanup_errors.append(f"{target}: {error}")
    checkpoints_root = Path(job.scratch_path) / "checkpoints"
    try:
        if checkpoints_root.is_dir() and not any(checkpoints_root.iterdir()):
            checkpoints_root.rmdir()
    except OSError as error:
        job.cleanup_errors.append(f"{checkpoints_root}: {error}")
    try:
        Path(job.assets_path).unlink()
    except OSError as error:
        job.cleanup_errors.append(f"{job.assets_path}: {error}")
    write_job(job, path)
    return job_status(path)


def run_job(
    path,
    *,
    progress=None,
    execution_memory_bytes=None,
):
    """Run or resume one prepared reconstruction job."""
    path = Path(path).absolute()
    job = read_job(path)
    if job.status == "complete":
        return job_status(path)
    verify_job(job)
    scan = job.scan
    config = job.config_data
    assets = _load_assets(job)
    spec = job.internal_spec()
    bounds = scan.exposure_angle_bounds(
        config, fallback=job.angle_fallback
    )
    ranges, tiles = _execution_layout(job, scan, config)
    checkpoint_dir = Path(job.scratch_path) / "checkpoints"
    boundaries = _job_checkpoint_boundaries(job)

    effective_memory = (
        max(1024**2, int(execution_memory_bytes))
        if execution_memory_bytes is not None
        else job.memory_override_bytes or job.runtime_memory_bytes
    )
    number_of_grids = max(1, len(spec.grids))
    # One split of the budget: what the accumulators may hold and what the
    # frame pipeline may use, so the two cannot each claim most of it.
    active_budget_bytes, effective_memory = split_memory_budget(
        effective_memory, number_of_grids
    )

    resumed, _existing_files = _discover_checkpoint_state(
        checkpoint_dir, boundaries, job.digest, cleanup_stale=True
    )
    router = _CheckpointRouter(
        boundaries,
        spec_digest=job.digest,
        checkpoint_dir=checkpoint_dir,
        active_budget_bytes=active_budget_bytes,
        resumed=resumed,
    )
    provenance = _base_provenance(job, config)
    correct = _correction_pipeline(config, scan, assets, provenance)

    pending_ranges = [
        frame_range
        for frame_range in ranges
        if not _range_is_resumed(router, boundaries, resumed, frame_range)
    ]
    total_images = sum(stop - start for start, stop in ranges)
    completed_images = total_images - sum(
        stop - start for start, stop in pending_ranges
    )
    if progress is not None:
        progress(
            completed_images,
            total_images + 1,
            f"Mapping images {completed_images}/{total_images}",
        )

    _map_pending_ranges(
        spec,
        scan,
        config,
        bounds,
        tiles,
        pending_ranges,
        router,
        correction_pipeline=correct,
        effective_memory=effective_memory,
        threads_per_image=job.threads_per_image,
        accumulation_budget_bytes=job.accumulation_budget_bytes,
        total_images=total_images,
        completed_images=completed_images,
        progress=progress,
    )

    if progress is not None:
        progress(total_images, total_images + 1, "Verifying checkpoints")
    resumed_after, checkpoint_files = _discover_checkpoint_state(
        checkpoint_dir, boundaries, job.digest, cleanup_stale=False
    )
    expected = {
        (grid_name, index)
        for grid_name, ranges_for_grid in boundaries.items()
        for index in range(len(ranges_for_grid))
    }
    missing = expected - resumed_after
    if missing:
        raise RuntimeError(
            "Reconstruction mapping did not produce "
            f"{len(missing)} expected checkpoint(s): {sorted(missing)[:10]}"
        )
    job.correction_provenance = provenance
    job.status = "finalizing"
    write_job(job, path)
    if progress is not None:
        progress(total_images, total_images + 1, "Finalizing HDF5")
    result = _finalize_reconstruction(
        spec,
        checkpoint_files,
        job.output_path,
        provenance={
            **provenance,
            "scan_reference": job.scan_reference,
        },
        config=config,
        chunk_progress=(
            None
            if progress is None
            else lambda written, count: progress(
                total_images,
                total_images + 1,
                f"Finalizing HDF5 chunk {written}/{count}",
            )
        ),
    )
    job.output_sha256 = result["sha256"]
    job.status = "complete"
    for target in (checkpoint_dir, Path(job.assets_path)):
        try:
            if target.is_dir():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()
        except OSError as error:
            job.cleanup_errors.append(f"{target}: {error}")
    write_job(job, path)
    if progress is not None:
        progress(total_images + 1, total_images + 1, "Complete")
    return job_status(path)

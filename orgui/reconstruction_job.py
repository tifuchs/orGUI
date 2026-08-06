"""Central orGUI job model for reciprocal-space reconstruction."""

from __future__ import annotations

from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import asdict, dataclass, field
from hashlib import sha256
import json
import math
from pathlib import Path
from queue import Empty, SimpleQueue
import shutil
from threading import Event
from typing import Any

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
    _calibration_probe_all_grids,
    _detector_corner_rays,
    _discover_checkpoint_state,
    _files_per_job,
    _finalize_reconstruction,
    _map_frame_range,
)


JOB_SCHEMA_VERSION = 4
ACCURACY_DEPTHS = {
    "center": 0,
    "low": 1,
    "balanced": 2,
    "high": 3,
    "very_high": 4,
    "maximum": 5,
}
AUTO_MAX_FRAMES_PER_TASK = 64
ACCUMULATION_TRANSIENT_FACTOR = 3
AUTO_MAX_ACCUMULATION_BYTES = 2 * 1024**3
MAX_CONCURRENT_ACTIVE_CHECKPOINTS = 2
"""Design doc Sec10: how many checkpoints (per grid) may have live,
non-empty accumulator state at once -- small and fixed, not user-facing."""


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
    threads_per_image: int
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
    work_block_pixels: int | None = None
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
        work_block = self.work_block_pixels or max(
            1024, min(65536, memory // max(threads * 256, 1))
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
        np.full(3, 2_000_000_000, dtype=np.int64),
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
        np.full(3, 2_000_000_000, dtype=np.int64),
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

    excluded = set(config.corrections.excluded_frames)
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
    threads_per_image=4,
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
    if int(threads_per_image) < 1:
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
        threads_per_image=int(threads_per_image),
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
            provenance.setdefault("factor_uncertainty", {})[
                name
            ] = "propagated"
        else:
            intensity *= factor
            variance *= factor**2
            provenance.setdefault("factor_uncertainty", {})[
                name
            ] = "deterministic-no-uncertainty"

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
        if static_factor is not None:
            intensity *= static_factor
            variance *= static_factor_squared
        if exposure is not None:
            value = frame_value(exposure, frame_index)
            if value <= 0 or not math.isfinite(value):
                raise ValueError("Exposure time must be finite and positive")
            factor = 1.0 / value
            factor_variance = None
            if exposure_variance is not None:
                value_variance = frame_value(
                    exposure_variance, frame_index
                )
                factor_variance = value_variance / value**4
            apply_factor(
                intensity,
                variance,
                factor,
                factor_variance,
                "exposure",
            )
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
            apply_factor(
                intensity,
                variance,
                1.0 / value,
                factor_variance,
                f"monitor:{name}",
            )
        mask |= ~np.isfinite(intensity) | ~np.isfinite(variance)
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


def _execution_layout(job, scan, config):
    rows, columns = config.detector.detector.shape
    effective_memory = (
        job.memory_override_bytes or job.runtime_memory_bytes
    )
    depth = (
        job.advanced_depth
        if job.advanced_depth is not None
        else ACCURACY_DEPTHS[job.accuracy]
    )
    worst_leaves = 8**depth
    estimated_native_bytes_per_pixel = 128 + 2 * worst_leaves * 40
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
    tile_side = max(1, int(math.sqrt(tile_pixels)))
    tile_rows, tile_columns = job.tile_shape or (
        min(tile_side, rows),
        min(tile_side, columns),
    )
    included_count = max(
        1,
        len(scan)
        - len(
            {
                frame
                for frame in config.corrections.excluded_frames
                if 0 <= frame < len(scan)
            }
        ),
    )
    if job.frame_batch is None:
        spec = job.internal_spec()
        native_threads = max(
            1, min(job.threads_per_image, spec.threads)
        )
        image_workers = max(1, spec.threads // native_threads)
        target_tasks = max(1, image_workers * 4)
        frame_batch = min(
            AUTO_MAX_FRAMES_PER_TASK,
            max(1, math.ceil(included_count / target_tasks)),
        )
    else:
        frame_batch = job.frame_batch
    ranges = _included_ranges(
        len(scan), set(config.corrections.excluded_frames), frame_batch
    )
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
    stationary,
    frames_per_task=1,
    threads_per_image=1,
    accumulation_budget_bytes=None,
):
    """Derive bounded parallelism across images and within each image.

    Each frame worker processes that frame's detector tiles sequentially. The
    configured CPU budget is divided between concurrent image workers and
    native threads used by each image. The worker count is also bounded by a
    depth-aware working-set estimate. Native blocks are reduced before task
    records are retained, so the unattainable bound in which every footprint
    leaf survives for the complete detector must not be multiplied by the
    number of image workers.

    ``accumulation_budget_bytes``/the returned ``accumulation`` value model a
    per-worker retained-record buffer from the retired Parquet-era mapping
    path; the checkpoint-routing path (design doc Sec9) has no equivalent
    per-worker buffer of its own; ``run_job`` computes the checkpoint
    accumulators' own memory budget separately (Sec10) and does not forward
    this value onward. Kept unchanged here (this function belongs to the
    still-unimplemented Sec7 thread-allocation phase) so a caller-supplied
    ``accumulation_budget_bytes`` still conservatively bounds
    ``image_workers``.

    :param _ReconstructionSpec spec:
        Frozen reconstruction compute settings.
    :param iterable tiles:
        Detector tiles as ``(row_start, row_stop, column_start,
        column_stop)`` tuples.
    :param int memory_bytes:
        Total memory budget in bytes.
    :param bool stationary:
        Whether exposure start and end angles are identical.
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
    detector_pixels = sum(tile_pixels)
    children = 4 if stationary else 8
    depth_scale = 1 << max(0, spec.max_depth - 2)
    sweep_scale = 2 if children == 8 else 1
    bytes_per_pixel = 128 * depth_scale * sweep_scale
    image_memory = detector_pixels * bytes_per_pixel
    worker_memory = max(
        1024**2,
        image_memory,
    )
    kernel_threads = max(
        1, min(int(threads_per_image), max(1, int(spec.threads)))
    )
    cpu_workers = max(1, int(spec.threads) // kernel_threads)
    minimum_accumulation = 1024**2
    if accumulation_budget_bytes is None:
        memory_workers = max(1, memory_bytes // worker_memory)
    else:
        requested = max(
            minimum_accumulation, int(accumulation_budget_bytes)
        )
        required = (
            worker_memory
            + ACCUMULATION_TRANSIENT_FACTOR * requested
        )
        memory_workers = max(1, memory_bytes // required)
    image_workers = max(
        1,
        min(cpu_workers, memory_workers),
    )
    safe_accumulation = max(
        minimum_accumulation,
        (
            memory_bytes // image_workers - worker_memory
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
                stationary=stationary,
                frames_per_task=stop - start,
                threads_per_image=job.threads_per_image,
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
            1, min(job.threads_per_image, spec.threads)
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


def job_status(path):
    """Return verified completion state for a job JSON."""
    job = read_job(path)
    boundaries = _job_checkpoint_boundaries(job)
    total_checkpoints = sum(len(items) for items in boundaries.values())
    if job.status != "complete":
        verify_job(job)
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


def run_cluster_map_task(
    path,
    task_index,
    *,
    cpus=1,
    memory_bytes=1024**3,
    progress=None,
):
    """Execute one cluster map-array task.

    Not yet implemented: cluster execution is being reworked for the
    checkpoint architecture (design doc Sec13), where each array task is a
    full independent node running the whole single-node pipeline against
    its own disjoint frame-range slice -- a structurally different split
    from the retired Parquet-era wiring this replaces (one frame-range
    mapping task feeding a shared job-wide reduce). Run the job with
    ``run``/``resume`` on a single node until cluster support lands.
    """
    raise NotImplementedError(
        "Cluster execution is being reworked for the checkpoint "
        "architecture (design doc Sec13) and is not yet available in this "
        "build. Use 'run'/'resume' on a single node instead."
    )


def run_cluster_finalize(
    path,
    *,
    cpus=1,
    memory_bytes=1024**3,
    progress=None,
):
    """Verify all cluster map tasks, then finalize their job.

    Not yet implemented -- see :func:`run_cluster_map_task`.
    """
    raise NotImplementedError(
        "Cluster execution is being reworked for the checkpoint "
        "architecture (design doc Sec13) and is not yet available in this "
        "build. Use 'run'/'resume' on a single node instead."
    )


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
    active_budget_bytes = max(
        1024**2,
        effective_memory
        // (MAX_CONCURRENT_ACTIVE_CHECKPOINTS * number_of_grids),
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

    def range_is_resumed(frame_range):
        start, stop = frame_range
        for grid_name in boundaries:
            for frame_index in range(start, stop):
                index = router.checkpoint_index_for_frame(grid_name, frame_index)
                if (grid_name, index) not in resumed:
                    return False
        return True

    pending_ranges = [
        frame_range for frame_range in ranges if not range_is_resumed(frame_range)
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

    worker_limits = []
    kernel_thread_limits = []
    scheduler_memory = max(
        1024**2,
        effective_memory - (ray_cache_bytes if cache_detector_rays else 0),
    )
    for frame_range in pending_ranges:
        task_bounds = bounds[frame_range[0] : frame_range[1]]
        stationary = bool(
            np.array_equal(task_bounds[:, 0], task_bounds[:, 1])
        )
        image_limit, native_threads, _worker_memory, _accumulation_limit = (
            _frame_parallelism(
                spec,
                tiles,
                scheduler_memory,
                stationary=stationary,
                frames_per_task=frame_range[1] - frame_range[0],
                threads_per_image=job.threads_per_image,
                accumulation_budget_bytes=job.accumulation_budget_bytes,
            )
        )
        worker_limits.append(image_limit)
        kernel_thread_limits.append(native_threads)
    image_workers = (
        min(len(pending_ranges), min(worker_limits)) if worker_limits else 1
    )
    kernel_threads = min(kernel_thread_limits) if kernel_thread_limits else 1
    # The native budget is a per-call guard. The scheduler independently
    # bounds the aggregate worker working set above.
    kernel_memory_budget = effective_memory

    progress_events = SimpleQueue()
    cancellation = Event()
    mapped_images = completed_images

    def publish_image_progress(frame_index):
        if cancellation.is_set():
            raise RuntimeError("Reconstruction mapping cancelled")
        progress_events.put(frame_index)

    def drain_progress_events():
        nonlocal mapped_images
        while True:
            try:
                progress_events.get_nowait()
            except Empty:
                return
            mapped_images += 1
            if progress is not None:
                progress(
                    mapped_images,
                    total_images + 1,
                    (
                        f"Mapping images {mapped_images}/{total_images} "
                        f"({image_workers} image workers, "
                        f"{kernel_threads} native threads/image)"
                    ),
                )

    def map_frame(frame_range):
        task_bounds = bounds[frame_range[0] : frame_range[1]]
        _map_frame_range(
            spec,
            scan,
            config.detector,
            config.ub_calculator,
            frame_range,
            tiles,
            task_bounds,
            router,
            correction_pipeline=correct,
            corner_rays=ray_cache if cache_detector_rays else None,
            kernel_threads=kernel_threads,
            kernel_memory_budget_bytes=kernel_memory_budget,
            image_progress=publish_image_progress,
        )
        return frame_range

    executor = (
        None
        if image_workers == 1
        else ThreadPoolExecutor(
            max_workers=image_workers,
            thread_name_prefix="orgui-rsmap-image",
        )
    )
    try:
        for wave_start in range(0, len(pending_ranges), image_workers):
            wave = pending_ranges[wave_start : wave_start + image_workers]
            if executor is None:
                for frame_task in wave:
                    map_frame(frame_task)
                    drain_progress_events()
            else:
                future_set = {
                    executor.submit(map_frame, frame_task) for frame_task in wave
                }
                while future_set:
                    done, future_set = wait(
                        future_set,
                        timeout=0.1,
                        return_when=FIRST_COMPLETED,
                    )
                    drain_progress_events()
                    for future in done:
                        future.result()
                drain_progress_events()
    except BaseException:
        cancellation.set()
        raise
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

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

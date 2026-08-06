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
    _GridSpec,
    _MIN_REDUCER_WORKER_MEMORY,
    _ReconstructionSpec,
    _TaskManifest,
    _calibration_probe_all_grids,
    _detector_corner_rays,
    _files_per_job,
    _finalize_reconstruction,
    _map_frame_range,
    _read_manifest,
    _reduce_partition,
    _verify_scratch_file,
    _write_manifest,
    _xxh3_128,
)


JOB_SCHEMA_VERSION = 3
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
    partition_chunk_span: int | None = None
    user_note: str = ""
    status: str = "prepared"
    expected_map_tasks: int = 0
    map_manifests: list[str] = field(default_factory=list)
    reduction_manifest: str | None = None
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
        """Return the immutable scientific job digest."""
        values = self.to_dict()
        for key in (
            "status",
            "map_manifests",
            "reduction_manifest",
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
        chunk_span = self.partition_chunk_span or max(
            16, min(4096, memory // (64**3 * 32))
        )
        return _ReconstructionSpec(
            grids=tuple(_GridSpec(**grid) for grid in self.grids),
            max_depth=depth,
            threads=threads,
            work_block_pixels=work_block,
            memory_budget_bytes=memory,
            partition_chunk_span=chunk_span,
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
    partition_chunk_span=None,
    threads_per_image=4,
    accumulation_budget_bytes=None,
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
        partition_chunk_span=partition_chunk_span,
        cluster_settings=dict(cluster_settings or {}),
    )
    ranges, tiles = _execution_layout(job, gui.fscan, config)
    job.expected_map_tasks = len(ranges)
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
        "maximum_parallel_reducer_workers": min(
            spec.threads,
            max(
                1,
                spec.memory_budget_bytes
                // _MIN_REDUCER_WORKER_MEMORY,
            ),
        ),
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
        "parquet_chunk_span": spec.partition_chunk_span,
        "frame_tasks": len(ranges),
        "detector_tiles": len(tiles),
        "map_tasks": len(ranges),
        "detector_ray_cache_MiB": (
            ray_cache_bytes / 1024**2 if cache_detector_rays else 0.0
        ),
        "parallel_layouts": list(layouts.values()),
    }


def job_status(path):
    """Return verified completion state for a job JSON."""
    job = read_job(path)
    map_manifests = list(job.map_manifests)
    if job.status != "complete":
        verify_job(job)
        scan = job.scan
        config = job.config_data
        ranges, tiles = _execution_layout(job, scan, config)
        discovered = _discover_map_manifests(
            job,
            ranges,
            tiles,
            job.internal_spec(),
            verify_partitions=False,
        )
        map_manifests = list(discovered.values())
    result = {
        "status": job.status,
        "job_sha256": job.digest,
        "map_tasks": {
            "completed": len(map_manifests),
            "pending": max(
                0, job.expected_map_tasks - len(map_manifests)
            ),
            "total": job.expected_map_tasks,
        },
        "reduction_manifest": job.reduction_manifest,
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


def _valid_map_manifest(
    path,
    frame_range,
    detector_tiles,
    *,
    spec_hash,
    job_digest,
    verification_cache=None,
    verify_partitions=True,
):
    try:
        manifest = _read_manifest(path)
        return (
            manifest.kind == "map"
            and manifest.status == "complete"
            and manifest.spec_hash == spec_hash
            and manifest.metadata.get("job_sha256") == job_digest
            and tuple(manifest.frame_range) == tuple(frame_range)
            and manifest.detector_tile is None
            and manifest.metadata.get("detector_tiles")
            == [list(tile) for tile in detector_tiles]
            and (
                not verify_partitions
                or all(
                    _verify_scratch_file(partition, verification_cache)
                    for partition in manifest.partitions
                )
            )
        )
    except (OSError, ValueError, RuntimeError):
        return False


def _discover_map_manifests(
    job,
    ranges,
    tiles,
    spec,
    *,
    verification_cache=None,
    only_ranges=None,
    verify_partitions=True,
):
    manifest_root = Path(job.scratch_path) / "manifests"
    candidates = {
        str(Path(value).absolute()) for value in job.map_manifests
    }
    if manifest_root.exists():
        candidates.update(
            str(path.absolute()) for path in manifest_root.glob("*.json")
        )
    expected = {
        tuple(frame_range)
        for frame_range in (
            ranges if only_ranges is None else only_ranges
        )
    }
    existing = {}
    for manifest_path in sorted(candidates):
        try:
            manifest = _read_manifest(manifest_path)
            frame_range = tuple(manifest.frame_range)
            if (
                frame_range in expected
                and _valid_map_manifest(
                    manifest_path,
                    frame_range,
                    tiles,
                    spec_hash=spec.digest,
                    job_digest=job.digest,
                    verification_cache=verification_cache,
                    verify_partitions=verify_partitions,
                )
            ):
                existing[frame_range] = manifest_path
        except (OSError, TypeError, ValueError, RuntimeError):
            continue
    return existing


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
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _merge_provenance(target[key], value)
        else:
            target[key] = value


def run_cluster_map_task(
    path,
    task_index,
    *,
    cpus=1,
    memory_bytes=1024**3,
    progress=None,
):
    """Execute one deterministic cluster map-array task.

    The task writes only its immutable map manifest and Parquet partitions;
    it never mutates the shared reconstruction job JSON.
    """
    path = Path(path).absolute()
    job = read_job(path)
    verify_job(job)
    scan = job.scan
    config = job.config_data
    assets = _load_assets(job)
    spec = job.internal_spec()
    ranges, tiles = _execution_layout(job, scan, config)
    task_index = int(task_index)
    if task_index < 0 or task_index >= len(ranges):
        raise IndexError(
            f"Map task index {task_index} is outside 0..{len(ranges) - 1}"
        )
    frame_range = ranges[task_index]
    cpus = max(1, int(cpus))
    memory_bytes = max(1024**2, int(memory_bytes))
    verification_cache = set()
    existing = _discover_map_manifests(
        job,
        ranges,
        tiles,
        spec,
        verification_cache=verification_cache,
        only_ranges=(frame_range,),
    )
    if tuple(frame_range) in existing:
        return {
            "status": "complete",
            "reused": True,
            "task_index": task_index,
            "frame_range": list(frame_range),
            "manifest": existing[tuple(frame_range)],
        }
    bounds = scan.exposure_angle_bounds(
        config, fallback=job.angle_fallback
    )
    task_bounds = bounds[frame_range[0] : frame_range[1]]
    stationary = bool(
        np.array_equal(task_bounds[:, 0], task_bounds[:, 1])
    )
    execution_spec = _ReconstructionSpec.from_dict(
        {
            **spec.to_dict(),
            "threads": cpus,
            "memory_budget_bytes": memory_bytes,
        }
    )
    _, kernel_threads, _, accumulation_budget = _frame_parallelism(
        execution_spec,
        tiles,
        memory_bytes,
        stationary=stationary,
        frames_per_task=frame_range[1] - frame_range[0],
        threads_per_image=cpus,
        accumulation_budget_bytes=job.accumulation_budget_bytes,
    )
    provenance = _base_provenance(job, config)
    correct = _correction_pipeline(config, scan, assets, provenance)
    ray_cache = {}
    for tile in tiles:
        rays = _detector_corner_rays(config.detector, tile)
        ray_cache[tuple(tile)] = (rays, _xxh3_128(rays))

    completed = 0

    def image_progress(frame_index, retained_bytes, segments):
        nonlocal completed
        completed += 1
        if progress is not None:
            progress(
                completed,
                frame_range[1] - frame_range[0],
                (
                    f"Cluster map task {task_index}; frame {frame_index}; "
                    f"{retained_bytes / 1024**2:.0f} MiB retained; "
                    f"{segments} segments flushed"
                ),
            )

    manifest = _map_frame_range(
        spec,
        scan,
        config.detector,
        config.ub_calculator,
        frame_range,
        tiles,
        task_bounds,
        Path(job.scratch_path) / "map",
        correction_pipeline=correct,
        job_digest=job.digest,
        corner_rays={
            tile: values[0] for tile, values in ray_cache.items()
        },
        corner_rays_fingerprints={
            tile: values[1] for tile, values in ray_cache.items()
        },
        verification_cache=verification_cache,
        kernel_threads=kernel_threads,
        kernel_memory_budget_bytes=memory_bytes,
        accumulation_budget_bytes=accumulation_budget,
        image_progress=image_progress,
    )
    manifest.metadata["correction_provenance"] = provenance
    manifest_root = Path(job.scratch_path) / "manifests"
    manifest_root.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_root / f"{manifest.task_id}.json"
    _write_manifest(manifest, manifest_path)
    return {
        "status": "complete",
        "reused": False,
        "task_index": task_index,
        "frame_range": list(frame_range),
        "manifest": str(manifest_path),
        "partitions": len(manifest.partitions),
    }


def run_cluster_finalize(
    path,
    *,
    cpus=1,
    memory_bytes=1024**3,
    progress=None,
):
    """Verify all cluster map tasks, then reduce and finalize their job."""
    path = Path(path).absolute()
    job = read_job(path)
    if job.status == "complete":
        return job_status(path)
    verify_job(job)
    scan = job.scan
    config = job.config_data
    spec = job.internal_spec()
    ranges, tiles = _execution_layout(job, scan, config)
    verification_cache = set()
    existing = _discover_map_manifests(
        job,
        ranges,
        tiles,
        spec,
        verification_cache=verification_cache,
    )
    missing = [
        index
        for index, frame_range in enumerate(ranges)
        if tuple(frame_range) not in existing
    ]
    if missing:
        preview = ", ".join(map(str, missing[:20]))
        suffix = "..." if len(missing) > 20 else ""
        raise RuntimeError(
            f"Cannot finalize: {len(missing)} map array tasks are missing "
            f"({preview}{suffix})"
        )
    job.map_manifests = [
        existing[tuple(frame_range)] for frame_range in ranges
    ]
    provenance = _base_provenance(job, config)
    for manifest_path in job.map_manifests:
        manifest = _read_manifest(manifest_path)
        _merge_provenance(
            provenance,
            manifest.metadata.get("correction_provenance", {}),
        )
    job.correction_provenance = provenance
    write_job(job, path)
    return run_job(
        path,
        progress=progress,
        execution_threads=max(1, int(cpus)),
        execution_memory_bytes=max(1024**2, int(memory_bytes)),
    )


def run_job(
    path,
    *,
    progress=None,
    execution_threads=None,
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
    manifest_root = Path(job.scratch_path) / "manifests"
    manifest_root.mkdir(parents=True, exist_ok=True)
    verification_cache = set()
    existing = _discover_map_manifests(
        job,
        ranges,
        tiles,
        spec,
        verification_cache=verification_cache,
    )
    job.map_manifests = list(existing.values())
    provenance = _base_provenance(job, config)
    correct = _correction_pipeline(config, scan, assets, provenance)
    total_tasks = len(ranges)
    completed_tasks = sum(
        tuple(frame_range) in existing
        for frame_range in ranges
    )
    total_images = sum(stop - start for start, stop in ranges)
    completed_images = sum(
        stop - start
        for start, stop in ranges
        if (start, stop) in existing
    )
    if progress is not None:
        progress(
            completed_images,
            total_images + 2,
            (
                f"Mapping images {completed_images}/{total_images} "
                f"({completed_tasks}/{total_tasks} tasks committed)"
            ),
        )
    effective_memory = (
        max(1024**2, int(execution_memory_bytes))
        if execution_memory_bytes is not None
        else job.memory_override_bytes or job.runtime_memory_bytes
    )
    reducer_threads = (
        max(1, int(execution_threads))
        if execution_threads is not None
        else spec.threads
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
            rays = _detector_corner_rays(config.detector, tile)
            ray_cache[tuple(tile)] = (rays, _xxh3_128(rays))

    pending_frames = [
        frame_range
        for frame_range in ranges
        if tuple(frame_range) not in existing
    ]

    worker_limits = []
    kernel_thread_limits = []
    accumulation_limits = []
    scheduler_memory = max(
        1024**2,
        effective_memory - (ray_cache_bytes if cache_detector_rays else 0),
    )
    for frame_range in pending_frames:
        task_bounds = bounds[frame_range[0] : frame_range[1]]
        stationary = bool(
            np.array_equal(task_bounds[:, 0], task_bounds[:, 1])
        )
        (
            image_limit,
            native_threads,
            _worker_memory,
            accumulation_limit,
        ) = _frame_parallelism(
            spec,
            tiles,
            scheduler_memory,
            stationary=stationary,
            frames_per_task=frame_range[1] - frame_range[0],
            threads_per_image=job.threads_per_image,
            accumulation_budget_bytes=job.accumulation_budget_bytes,
        )
        worker_limits.append(image_limit)
        kernel_thread_limits.append(native_threads)
        accumulation_limits.append(accumulation_limit)
    image_workers = (
        min(len(pending_frames), min(worker_limits))
        if worker_limits
        else 1
    )
    kernel_threads = (
        min(kernel_thread_limits) if kernel_thread_limits else 1
    )
    accumulation_budget = (
        min(accumulation_limits) if accumulation_limits else 1024**2
    )
    # The native budget is a per-call guard. The scheduler independently
    # bounds the aggregate worker working set above.
    kernel_memory_budget = effective_memory
    progress_events = SimpleQueue()
    cancellation = Event()
    mapped_images = completed_images
    parallel_mapping = image_workers > 1

    def publish_image_progress(frame_index, retained_bytes, segments):
        if cancellation.is_set():
            raise RuntimeError("Reconstruction mapping cancelled")
        event = (frame_index, retained_bytes, segments)
        if parallel_mapping:
            progress_events.put(event)
        else:
            report_image_progress(event)

    def report_image_progress(event):
        nonlocal mapped_images
        frame_index, retained_bytes, segments = event
        mapped_images += 1
        if progress is not None:
            progress(
                mapped_images,
                total_images + 2,
                (
                    f"Mapping images {mapped_images}/{total_images}; "
                    f"frame {frame_index}; "
                    f"{completed_tasks}/{total_tasks} tasks committed; "
                    f"{retained_bytes / 1024**2:.0f} MiB retained; "
                    f"{segments} segments flushed"
                ),
            )

    def drain_progress_events():
        while True:
            try:
                event = progress_events.get_nowait()
            except Empty:
                return
            report_image_progress(event)

    def map_frame(frame_range):
        task_bounds = bounds[frame_range[0] : frame_range[1]]
        return frame_range, _map_frame_range(
            spec,
            scan,
            config.detector,
            config.ub_calculator,
            frame_range,
            tiles,
            task_bounds,
            Path(job.scratch_path) / "map",
            correction_pipeline=correct,
            job_digest=job.digest,
            corner_rays={
                tile: values[0] for tile, values in ray_cache.items()
            },
            corner_rays_fingerprints={
                tile: values[1] for tile, values in ray_cache.items()
            },
            verification_cache=verification_cache,
            kernel_threads=kernel_threads,
            kernel_memory_budget_bytes=kernel_memory_budget,
            accumulation_budget_bytes=accumulation_budget,
            image_progress=publish_image_progress,
        )

    executor = (
        None
        if image_workers == 1
        else ThreadPoolExecutor(
            max_workers=image_workers,
            thread_name_prefix="orgui-rsmap-image",
        )
    )
    try:
        for wave_start in range(0, len(pending_frames), image_workers):
            wave = pending_frames[
                wave_start : wave_start + image_workers
            ]
            if executor is None:
                mapped_frames = [map_frame(frame_task) for frame_task in wave]
            else:
                future_set = {
                    executor.submit(map_frame, frame_task)
                    for frame_task in wave
                }
                mapped_frames = []
                while future_set:
                    done, future_set = wait(
                        future_set,
                        timeout=0.1,
                        return_when=FIRST_COMPLETED,
                    )
                    drain_progress_events()
                    mapped_frames.extend(
                        future.result() for future in done
                    )
                drain_progress_events()
                mapped_frames.sort(key=lambda item: item[0])
            # This bounded worker wave has finished, so correction provenance
            # is stable while the coordinator updates resumable state.
            for _frame_range, manifest in mapped_frames:
                manifest.metadata["correction_provenance"] = provenance
                manifest_path = (
                    manifest_root / f"{manifest.task_id}.json"
                )
                _write_manifest(manifest, manifest_path)
                job.map_manifests.append(str(manifest_path))
                completed_tasks += 1
            job.correction_provenance = provenance
            write_job(job, path)
            if progress is not None:
                progress(
                    mapped_images,
                    total_images + 2,
                    (
                        f"Mapping images {mapped_images}/{total_images}; "
                        f"{completed_tasks}/{total_tasks} tasks committed "
                        f"({image_workers} image workers, "
                        f"{kernel_threads} native threads/image, "
                        f"{accumulation_budget / 1024**2:.0f} MiB "
                        "accumulation/worker)"
                    ),
                )
    except BaseException:
        cancellation.set()
        raise
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)
    if progress is not None:
        progress(
            total_images,
            total_images + 2,
            f"Reducing {len(job.map_manifests)} mapping tasks",
        )
    if job.map_manifests:
        reduced = _reduce_partition(
            job.map_manifests,
            Path(job.scratch_path) / "reduced",
            verification_cache=verification_cache,
            memory_budget_bytes=effective_memory,
            workers=reducer_threads,
            checkpoint_root=manifest_root,
            progress=(
                None
                if progress is None
                else lambda done, count, message: progress(
                    total_images,
                    total_images + 2,
                    f"Reducing shards {done}/{count}: {message}",
                )
            ),
        )
    else:
        reduced = _TaskManifest(
            kind="reduce",
            task_id=sha256((spec.digest + job.digest).encode()).hexdigest()[:24],
            spec_hash=spec.digest,
            status="complete",
            spec=spec.to_dict(),
            metadata={"job_sha256": job.digest, "empty_input": True},
        )
    reduced_path = manifest_root / f"{reduced.task_id}.json"
    _write_manifest(reduced, reduced_path)
    job.reduction_manifest = str(reduced_path)
    job.status = "finalizing"
    write_job(job, path)
    if progress is not None:
        progress(
            total_images + 1,
            total_images + 2,
            "Finalizing HDF5",
        )
    result = _finalize_reconstruction(
        [reduced_path],
        job.output_path,
        provenance={
            **provenance,
            "scan_reference": job.scan_reference,
        },
        config=config,
        verification_cache=verification_cache,
        chunk_progress=(
            None
            if progress is None
            else lambda written, count: progress(
                total_images + 1,
                total_images + 2,
                f"Finalizing HDF5 chunk {written}/{count}",
            )
        ),
    )
    job.output_sha256 = result["sha256"]
    job.status = "complete"
    for target in (
        Path(job.scratch_path) / "map",
        Path(job.scratch_path) / "reduced",
        manifest_root,
        Path(job.assets_path),
    ):
        try:
            if target.is_dir():
                shutil.rmtree(target)
            elif target.exists():
                target.unlink()
        except OSError as error:
            job.cleanup_errors.append(f"{target}: {error}")
    write_job(job, path)
    if progress is not None:
        progress(total_images + 2, total_images + 2, "Complete")
    return job_status(path)

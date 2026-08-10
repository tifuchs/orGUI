"""Out-of-core reciprocal-space reconstruction.

The numerical hot path is implemented by
``_reciprocal_reconstruction_cpp``.  This module owns the scientific Python
boundary, uncertainty propagation, checkpointed HDF5 scratch storage, and
final NeXus/HDF5 serialization.

Momentum-transfer grids use ``Angstrom^-1``.  HKL grids use reciprocal lattice
units.  All diffractometer angles accepted here are in radians.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from hashlib import sha256
import bisect
import importlib
import json
import math
import re
from pathlib import Path
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
_CHECKPOINT_BYTES_PER_ROW = 40
"""Fixed on-the-wire row size for the six ``_PARTIAL_COLUMNS`` fields
(two uint32 key columns, three float64 value columns, one uint32 count
column -- 2 * 4 + 3 * 8 + 4 = 40 bytes). Used to convert calibration-probe
record counts into byte estimates for the checkpoint file-count formula."""

_CHECKPOINT_PART_PATTERN = re.compile(r"^ckpt(?P<index>\d+)_p(?P<part>\d+)\.h5$")


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
        # The native kernel packs each record's within-chunk voxel index
        # and chunk index into uint32 fields (RecordKey.local/.chunk):
        # local's range is bounded by one chunk's own voxel count, chunk's
        # by the grid's total chunk count -- both computed from this same
        # (unclamped) chunk_shape. A config that could overflow either
        # bound must be rejected here, not silently corrupt indices deep
        # inside the native kernel.
        uint32_max = 2**32 - 1
        local_bound = math.prod(self.chunk_shape)
        if local_bound > uint32_max:
            raise ValueError(
                f"chunk_shape {self.chunk_shape} packs {local_bound:,} voxels per "
                f"chunk, which exceeds the native kernel's uint32 within-chunk "
                f"voxel index (max {uint32_max:,}); use a smaller chunk_shape"
            )
        chunk_bound = math.prod(
            math.ceil(size / chunk) for size, chunk in zip(self.shape, self.chunk_shape)
        )
        if chunk_bound > uint32_max:
            raise ValueError(
                f"grid shape {self.shape} with chunk_shape {self.chunk_shape} "
                f"produces {chunk_bound:,} chunks, which exceeds the native "
                f"kernel's uint32 chunk index (max {uint32_max:,}); use a "
                "larger chunk_shape or a smaller grid"
            )
        # The native kernel identifies a voxel internally by packing the
        # three axis indices into one uint64, one bit field per axis (see
        # voxel_id()/record_key()), so that the chunk/local split needs
        # shifts rather than divisions. Rounding each axis up to a whole
        # number of bits costs at most three bits in total, so this only
        # rejects grids already within a factor of eight of overflowing a
        # 64-bit voxel count.
        voxel_bits = sum(max(1, int(size - 1).bit_length()) for size in self.shape)
        if voxel_bits > 64:
            raise ValueError(
                f"grid shape {self.shape} needs {voxel_bits} bits of voxel "
                "index, which exceeds the native kernel's 64-bit packed voxel "
                "identifier; use a coarser step or a smaller grid"
            )
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
    """Configuration shared by mapping, checkpointing, and finalization."""

    grids: tuple[_GridSpec, ...]
    max_depth: int = 2
    threads: int = 1
    #: Pixels per native work block. The best value depends on ``max_depth``
    #: (see ReconstructionJob.internal_spec, which derives it); 4096 is the
    #: measured optimum for this class's own default depth of 2.
    work_block_pixels: int = 4096
    #: Frames mapped in one native call. On a rotation scan two adjacent
    #: frames are as close in reciprocal space as two adjacent pixels are
    #: (measured 0.71 voxels per frame against 0.49 and 0.80 for the next
    #: column and row), so a work block that is a brick in (row, column,
    #: frame) collapses redundancy a per-image block cannot see. One frame
    #: per call is the degenerate case of the same machinery, not a
    #: separate path. Past roughly one voxel of travel per frame the gain
    #: is gone, so this is a per-job quantity rather than a constant --
    #: ``None`` measures the job and chooses (the normal case), an
    #: explicit integer overrides that measurement.
    frames_per_group: int | None = None
    memory_budget_bytes: int = 512 * 1024 * 1024
    checkpoint_count: int = 10
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
        if self.frames_per_group is not None and self.frames_per_group < 1:
            raise ValueError("frames_per_group must be positive when set")
        if self.memory_budget_bytes < 1024 * 1024:
            raise ValueError("memory_budget_bytes must be at least 1 MiB")
        if self.checkpoint_count < 1:
            raise ValueError("checkpoint_count must be at least one")

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


def _admissible_sample_pixels(kernel, angles_start, angles_end):
    """Largest pixel count one ``accumulate()`` call will accept.

    Mirrors the kernel's own per-call precheck, which refuses a tile whose
    worst-case record footprint would exceed its memory budget. That limit
    shrinks by the subdivision factor at every adaptive depth -- a 10 GB
    budget admits about 50 million pixels at depth 0 but only 128 thousand
    at depth 5, and a rotating exposure subdivides eight ways rather than
    four -- so a sampler using a fixed pixel ceiling eventually asks for
    more than the kernel will take and gets an exception instead of an
    estimate.

    :returns:
        Pixel count that is safe to pass to ``accumulate()``.
    :rtype: int
    """
    configuration = kernel.configuration
    children = 4 if np.array_equal(angles_start, angles_end) else 8
    worst_leaves = children ** int(configuration["max_depth"])
    bytes_per_pixel = 128 + 2 * worst_leaves * _CHECKPOINT_BYTES_PER_ROW
    return max(
        1, int(configuration["memory_budget_bytes"]) // bytes_per_pixel
    )


def _representative_tile_origin(mask, side):
    """Pick a bootstrap-tile origin that reflects the whole detector.

    Anchoring the bootstrap at the detector corner measures a per-pixel
    cost well below average: a corner maps partly outside the output grid
    and is frequently masked, and on a real detector it produced half the
    records per pixel of the centre. Since the sized pass divides the
    remaining time budget by that rate, an unrepresentatively cheap
    bootstrap inflates every later sample. Candidates are tried from the
    centre outwards and the least-masked one wins, so a masked beam stop
    at the centre does not simply move the problem.

    :returns:
        ``(row_start, column_start)`` for a ``side``-by-``side`` tile.
    :rtype: tuple[int, int]
    """
    rows, columns = mask.shape
    best = None
    # Centre first, then progressively further out, so an unmasked
    # detector keeps the most representative position and a masked one
    # (a beam stop, a dead module) still finds clean pixels nearby.
    fractions = (0.5, 0.3, 0.7, 0.15, 0.85)
    for row_fraction in fractions:
        for column_fraction in fractions:
            row = min(
                max(0, int(rows * row_fraction) - side // 2), max(0, rows - side)
            )
            column = min(
                max(0, int(columns * column_fraction) - side // 2),
                max(0, columns - side),
            )
            masked = float(mask[row : row + side, column : column + side].mean())
            if masked == 0.0:
                return row, column
            if best is None or masked < best[0]:
                best = (masked, row, column)
    return best[1], best[2]


_MINIMUM_TILE_PIXELS = 1
"""Smallest sample tile the probe will ask for. Only ever reached when the
per-pixel cost is so high that a single pixel already fills a useful share
of the budget, which is the case at the deepest subdivision levels."""

_BOOTSTRAP_TIME_SHARE = 0.02
"""Fraction of the budget the first tile should take before the probe stops
enlarging it. Growing into this rather than fixing a size in advance is
what keeps the first measurement from costing more than the whole budget
when the per-pixel cost is extreme."""

_MINIMUM_CONVERGENCE_TILES = 6
"""Tiles required before the probe may stop early on its own error
estimate: a standard error computed from two or three samples is itself
too noisy to act on."""


def _sobol_tile_origins(rows, columns, count, rng=None):
    """Low-discrepancy tile origins, as fractions of the placeable range.

    The probe averages quantities that vary smoothly across the detector,
    so where its handful of tiles land decides how good the average is.
    Independent uniform draws clump: at eight tiles it is entirely
    ordinary for a quadrant to go unsampled. A scrambled Sobol sequence
    stratifies by construction while staying randomised, so repeated
    probes of the same job do not all inherit one unlucky layout.

    Requested counts are rounded up to a power of two, which is where a
    Sobol sequence's balance properties hold.

    :param int count:
        Lower bound on how many origins to generate.
    :param rng:
        Optional generator, used only to seed the scramble so tests can
        pin the sequence.
    :returns:
        List of ``(row_fraction, column_fraction)`` in ``[0, 1]``.
    :rtype: list[tuple[float, float]]
    """
    # Imported here rather than at module scope: scipy.stats is expensive
    # to import and this module is pulled in by the GUI at start-up.
    from scipy.stats import qmc

    exponent = max(1, math.ceil(math.log2(max(2, int(count)))))
    # A fixed default seed rather than an arbitrary one: the scramble is
    # there to stratify, not to randomise, and pinning it means two
    # preparations of the same job on the same machine sample the same
    # positions. Callers wanting independent draws pass their own
    # generator.
    seed = 0 if rng is None else int(rng.integers(0, 2**32))
    sampler = qmc.Sobol(d=2, scramble=True, seed=seed)
    points = sampler.random_base2(exponent)
    return [(float(point[0]), float(point[1])) for point in points]


def _relative_standard_error(values):
    """Standard error of the mean of ``values``, relative to that mean.

    Reported alongside the probe's estimate so a caller can tell a
    well-converged number from one drawn on a detector whose record
    density varies strongly with position, and size its safety margin on
    that rather than on a fixed guess.

    :rtype: float
    """
    if len(values) < 2:
        return float("inf")
    mean = sum(values) / len(values)
    if mean <= 0.0:
        return 0.0 if all(value == 0.0 for value in values) else float("inf")
    variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
    return math.sqrt(variance / len(values)) / mean


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
    target_relative_error=0.05,
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
        Side length in pixels of the first tile. Deliberately small: the
        adaptive loop grows subsequent tiles from the rate this one
        measures, so starting large only risks spending the whole budget
        before anything has been learned.
    :param float target_relative_error:
        Stop early once the sampled record density's relative standard
        error falls below this, rather than spending the whole budget on
        a detector that is already well characterised.
    :param int sample_tiles:
        Target number of tiles to spread the budget over. Rounded up to a
        power of two for the Sobol sequence; sampling stops on budget or
        convergence, so it is an upper bound rather than a fixed count.
    :param int min_sample_pixels, max_sample_pixels:
        Clamp on the sized pass's total target sample size.
    :param rng:
        Optional :class:`numpy.random.Generator` for deterministic tests;
        defaults to a fresh, unseeded generator.
    :returns:
        Dict with ``kernel_threads``, ``sampled_pixels``,
        ``sampled_seconds``, ``sampled_tiles``, ``seconds_per_pixel``
        (kernel CPU-seconds per sampled pixel), ``records_per_pixel``
        (post-dedup records per sampled pixel), and
        ``records_per_pixel_relative_error``, the standard error of that
        density relative to itself.
    :rtype: dict
    """
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

    admissible_pixels = _admissible_sample_pixels(
        kernel, angles_start, angles_end
    )
    admissible_side = max(1, int(math.isqrt(admissible_pixels)))
    # Grow the first tile until it is long enough to time, instead of
    # fixing its size up front. At the deepest subdivision a single pixel
    # already costs tens of milliseconds, so any fixed starting size
    # spends several budgets before the probe has learned anything.
    bootstrap_cap = max(
        1, min(int(bootstrap_tile), rows, columns, admissible_side)
    )
    total_pixels = 0
    total_seconds = 0.0
    total_records = 0
    # Per-tile record densities, kept so the estimate can report how well
    # it has converged rather than only what it averaged to.
    densities = []
    bootstrap_side = 1
    while True:
        bootstrap_row, bootstrap_column = _representative_tile_origin(
            mask, bootstrap_side
        )
        tile_started = time.perf_counter()
        pixels, seconds, records = sample_tile(
            bootstrap_row,
            bootstrap_row + bootstrap_side,
            bootstrap_column,
            bootstrap_column + bootstrap_side,
        )
        total_pixels += pixels
        total_seconds += seconds
        total_records += records
        densities.append(records / max(1, pixels))
        if (
            time.perf_counter() - tile_started
            >= budget_seconds * _BOOTSTRAP_TIME_SHARE
            or bootstrap_side >= bootstrap_cap
            or time.perf_counter() - started >= budget_seconds
        ):
            break
        bootstrap_side = min(bootstrap_side * 2, bootstrap_cap)
    # The ramp's earlier tiles exist to time the kernel, not to measure
    # density: a one-pixel tile yields a density of exactly zero or one and
    # would dominate the scatter the error estimate is built from. Their
    # pixels still count towards the totals, only their spread is dropped.
    densities = densities[-1:]
    # Low-discrepancy positions rather than independent uniform draws: the
    # quantities being averaged vary smoothly across the detector (solid
    # angle, incidence angle, which part of the output grid a pixel
    # reaches), so this is a quadrature problem, and at the handful of
    # tiles a 0.1 s budget allows, independent draws leave whole regions
    # unsampled often enough to matter.
    origins = _sobol_tile_origins(rows, columns, int(sample_tiles), rng)

    previous_pixels = total_pixels
    for fraction_row, fraction_column in origins:
        elapsed = time.perf_counter() - started
        remaining_budget = budget_seconds - elapsed
        if remaining_budget <= 0 or total_pixels >= max_sample_pixels:
            break
        if (
            len(densities) >= _MINIMUM_CONVERGENCE_TILES
            and total_pixels >= min_sample_pixels
            and _relative_standard_error(densities) <= target_relative_error
        ):
            break
        rate = total_seconds / max(1, total_pixels)
        if rate <= 0:
            next_pixels = max_sample_pixels
        else:
            # Spend the remaining budget over the tiles still planned,
            # re-estimated after every tile: a first tile that happened to
            # be unusually cheap costs one more small tile rather than the
            # whole budget, which is what the old size-once design did.
            share = remaining_budget / max(1, int(sample_tiles) - len(densities) + 1)
            next_pixels = share / rate
        next_pixels = min(
            next_pixels,
            max_sample_pixels - total_pixels,
            admissible_pixels,
            # Bound how fast the tile may grow, so one underestimated rate
            # cannot produce a single tile that overruns the budget.
            8 * previous_pixels,
        )
        # The floor is per tile and deliberately tiny: min_sample_pixels is
        # a target for the sample as a whole, and forcing it on every tile
        # would make a single tile cost more than the entire budget once
        # the per-pixel cost is high enough (a 256-pixel tile takes a
        # quarter of a second at the deepest subdivision).
        next_pixels = max(next_pixels, _MINIMUM_TILE_PIXELS)
        tile_side = max(1, int(math.isqrt(int(next_pixels))))
        tile_side = min(tile_side, rows, columns, admissible_side)
        row_start = int(round(fraction_row * max(0, rows - tile_side)))
        column_start = int(round(fraction_column * max(0, columns - tile_side)))
        pixels, seconds, records = sample_tile(
            row_start,
            row_start + tile_side,
            column_start,
            column_start + tile_side,
        )
        total_pixels += pixels
        total_seconds += seconds
        total_records += records
        densities.append(records / max(1, pixels))
        previous_pixels = pixels

    return {
        "kernel_threads": int(kernel_threads),
        "sampled_pixels": total_pixels,
        "sampled_seconds": total_seconds,
        "sampled_tiles": len(densities),
        "seconds_per_pixel": total_seconds / max(1, total_pixels),
        "records_per_pixel": total_records / max(1, total_pixels),
        "records_per_pixel_relative_error": _relative_standard_error(densities),
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


def _kernel_threads_sweep(
    spec,
    grid,
    ub_calculator,
    mask,
    corner_rays,
    angles_start,
    angles_end,
    *,
    candidates,
    tile_pixels=1_048_576,
    frame_pixels=None,
    memory_budget_bytes=None,
    plateau_ratio=0.9,
):
    """Measure real wall-clock ``per_call_time(kernel_threads)`` (design doc
    Sec7's live joint-balancing rule), one real kernel and one real
    ``accumulate()`` call per candidate against a single tile near
    ``tile_pixels``.

    Deliberately not :func:`_calibration_probe`'s
    ``block_mapping_cpu_seconds``: that figure is aggregate CPU-time summed
    across ``kernel_threads`` native worker threads, which hides the
    parallel-efficiency plateau by construction (it grows with thread count
    even when wall-clock time has flattened). Also deliberately one large
    tile rather than many small scattered ones -- the plateau is a property
    of one call's native block count (design doc Sec7/Sec17's confirmed
    sweep), not something many small calls would reveal.

    :param _ReconstructionSpec spec:
        Frozen reconstruction compute settings (``max_depth``,
        ``work_block_pixels`` are read from here; ``threads`` is not --
        each candidate supplies its own thread count).
    :param _GridSpec grid:
        The grid to build sweep kernels for.
    :param ub_calculator:
        Live UB calculator (same object passed to :func:`_build_kernels`).
    :param np.ndarray mask:
        Boolean mask for the full detector, shape ``(rows, columns)``.
    :param np.ndarray corner_rays:
        Full-detector corner rays, shape ``(rows + 1, columns + 1, 3)``.
    :param np.ndarray angles_start, angles_end:
        Shape ``(4,)`` diffractometer angles in radians.
    :param iterable candidates:
        Candidate ``kernel_threads`` values, sorted ascending.
    :param int tile_pixels:
        Target sample-tile pixel count; clamped to the detector's actual
        size.
    :param memory_budget_bytes:
        Forwarded to each candidate's kernel construction.
    :param float plateau_ratio:
        Stop measuring further candidates once a candidate's time is at
        least this fraction of the previous candidate's time (diminishing
        returns). Later, unmeasured candidates are not claimed to be
        faster: they simply inherit the last measured time.
    :param frame_pixels:
        Total pixels one frame maps, across every detector tile. The
        sample tile is a fraction of a frame -- and is clamped to the
        detector besides -- so a caller reasoning about how many frames
        can be in flight needs the time to map a whole one, not one
        sample tile. Given this, results are scaled to a whole frame;
        omitted, they stay per sample tile.
    :returns:
        Dict keyed by ``kernel_threads`` (every requested candidate,
        including unmeasured ones past the plateau), values are wall-clock
        seconds to map one frame when ``frame_pixels`` is given, or one
        ``tile_pixels`` call otherwise.
    :rtype: dict[int, float]
    """
    rows, columns = mask.shape
    tile_side = max(1, min(int(math.sqrt(tile_pixels)), rows, columns))
    # Scale from what was actually measured, not from the requested tile
    # size: the sample tile is clamped to the detector, so on a small
    # detector the two differ.
    sampled_pixels = tile_side * tile_side
    frame_scale = (
        1.0
        if frame_pixels is None
        else max(1.0, float(frame_pixels) / max(1, sampled_pixels))
    )
    tile_mask = np.ascontiguousarray(mask[:tile_side, :tile_side])
    tile_rays = np.ascontiguousarray(
        corner_rays[: tile_side + 1, : tile_side + 1]
    )
    intensity = np.zeros((tile_side, tile_side), dtype=np.float64)
    variance = np.zeros((tile_side, tile_side), dtype=np.float64)

    results: dict[int, float] = {}
    last_measured_time = None
    plateaued = False
    for kernel_threads in candidates:
        if plateaued:
            results[kernel_threads] = last_measured_time * frame_scale
            continue
        kernel = _kernel_for_grid(
            spec,
            grid,
            ub_calculator,
            threads=kernel_threads,
            memory_budget_bytes=memory_budget_bytes,
        )
        started = time.perf_counter()
        kernel.accumulate(
            intensity, variance, tile_mask, tile_rays,
            angles_start, angles_end, True,
        )
        elapsed = time.perf_counter() - started
        results[kernel_threads] = elapsed * frame_scale
        if (
            last_measured_time is not None
            and last_measured_time > 0
            and elapsed >= plateau_ratio * last_measured_time
        ):
            plateaued = True
        last_measured_time = elapsed
    return results


def _empty_batch() -> dict[str, np.ndarray]:
    return {
        "chunk_id": np.empty(0, dtype=np.uint32),
        "local_voxel_id": np.empty(0, dtype=np.uint32),
        "weighted_intensity": np.empty(0, dtype=np.float64),
        "weighted_variance": np.empty(0, dtype=np.float64),
        "weight": np.empty(0, dtype=np.float64),
        "contributors": np.empty(0, dtype=np.uint32),
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
    :func:`_finalize_reconstruction`.

    :param batch:
        A native record batch (dict of the six ``_PARTIAL_COLUMNS``
        arrays).
    :param path:
        Destination path for the checkpoint file.
    :param str spec_digest:
        The job's resume-identity digest, stored as an HDF5 attribute so
        resume can detect a stale checkpoint by comparing digests (design
        doc Sec11) without a separate manifest registry.
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


def _checkpoint_part_path(checkpoint_dir, grid_name, index, part):
    """Path of one checkpoint part file (design doc Sec9/Sec10/Sec11)."""
    return (
        Path(checkpoint_dir)
        / grid_name
        / f"ckpt{int(index):04d}_p{int(part):04d}.h5"
    )


class _CheckpointRouter:
    """Route each frame's per-grid reduced batch to its checkpoint file.

    Owns the state one job's mapping phase shares across every
    concurrently running frame-range worker: which planned checkpoint a
    given frame belongs to per grid (Sec5/Sec12), the live in-memory
    :class:`_CheckpointAccumulator` for each currently active checkpoint
    (Sec9), and the safety-valve mid-range split when a checkpoint's own
    memory-budget share is exceeded before its planned frame range drains
    (Sec10). :meth:`route` is safe to call concurrently from multiple
    frame workers; routing decisions for the whole job currently serialize
    behind one lock, since per-call bookkeeping cost is small relative to
    the image-load/correction/kernel work that happens before a worker
    ever calls :meth:`route` (mirroring Sec9's own rationale for why one
    lock per checkpoint is cheap enough) -- splitting this into a
    finer-grained per-checkpoint lock is a candidate throughput
    refinement for the Sec7 thread-allocation phase, not a correctness
    concern here.

    Each planned checkpoint (one contiguous frame range per grid, decided
    once at prepare time and passed in via ``boundaries``) is written as
    one or more numbered "part" files under
    ``checkpoint_dir/{grid_name}/ckpt{index:04d}_p{part:04d}.h5``: exactly
    one part in the common case, more than one only when the safety valve
    trips mid-range. Every part file records how many distinct frames
    contributed to it (``frames_covered`` attribute); a planned checkpoint
    is resumable exactly when its parts' ``frames_covered`` values sum to
    its full planned frame count (see :func:`_discover_checkpoint_state`).
    """

    def __init__(
        self,
        boundaries: Mapping[str, Iterable[tuple[int, int]]],
        *,
        spec_digest: str,
        checkpoint_dir,
        active_budget_bytes: int,
        resumed: Iterable[tuple[str, int]] = (),
    ):
        self._lock = threading.Lock()
        self._spec_digest = spec_digest
        self._checkpoint_dir = Path(checkpoint_dir)
        self._active_budget_bytes = max(1, int(active_budget_bytes))
        self._boundaries = {
            grid_name: [tuple(map(int, item)) for item in ranges]
            for grid_name, ranges in boundaries.items()
        }
        self._starts = {
            grid_name: [item[0] for item in ranges]
            for grid_name, ranges in self._boundaries.items()
        }
        self._resumed = set(resumed)
        self._accumulators: dict[tuple[str, int], _CheckpointAccumulator] = {}
        self._parts: dict[tuple[str, int], int] = {}
        self._frames_in_part: dict[tuple[str, int], int] = {}
        self._remaining: dict[tuple[str, int], int] = {}
        # Count of insert() calls currently in flight per key -- lets route()
        # release the router lock before calling the (potentially slow)
        # accumulator merge, while still knowing exactly when it's safe to
        # flush: not merely when every frame has been routed, but when every
        # one of their insert() calls has actually finished.
        self._inflight: dict[tuple[str, int], int] = {}
        for grid_name, ranges in self._boundaries.items():
            for index, (start, stop) in enumerate(ranges):
                key = (grid_name, index)
                self._remaining[key] = 0 if key in self._resumed else stop - start
        self.written: list[Path] = []

    def checkpoint_index_for_frame(self, grid_name, frame_index) -> int:
        """Return the planned checkpoint index containing ``frame_index``."""
        starts = self._starts[grid_name]
        position = bisect.bisect_right(starts, int(frame_index)) - 1
        if position < 0:
            raise ValueError(
                f"Frame {frame_index} precedes grid {grid_name!r}'s first "
                "planned checkpoint boundary"
            )
        start, stop = self._boundaries[grid_name][position]
        if not start <= frame_index < stop:
            raise ValueError(
                f"Frame {frame_index} is outside grid {grid_name!r}'s "
                "planned checkpoint boundaries"
            )
        return position

    def route(self, grid_name, frame_index, batch, *, frames=1) -> None:
        """Insert one reduced batch into its checkpoint accumulator.

        ``frames`` is how many frames the batch covers, starting at
        ``frame_index``: one for a single-frame batch, ``F`` for a frame
        group mapped in one native call, whose contributions merge inside
        the kernel and can no longer be told apart afterwards. The
        checkpoint's own remaining-frame countdown and the
        ``frames_covered`` attribute that makes a part resumable both
        count frames, not ``route()`` calls, so a group must declare its
        own size here or its checkpoint would never reach zero and could
        never be resumed from.

        A group must lie entirely inside one planned checkpoint range and
        must be contiguous in frame index; both are checked, since getting
        either wrong corrupts the resume bookkeeping silently rather than
        loudly.

        A no-op if that ``(grid, checkpoint)`` was already resumed from a
        previous run. May flush the checkpoint (or a part of it, under the
        Sec10 safety valve) to disk; the write itself happens without
        holding the router's lock.

        The merge itself (:meth:`_CheckpointAccumulator.insert`) runs
        outside the router's lock, relying on the accumulator's own lock
        for correctness -- only the short bookkeeping around it is
        serialized job-wide. This matters because ``insert()`` is not
        cheap: its binary-counter merge tree periodically cascades
        pairwise merges over progressively larger accumulated batches, and
        holding a single job-wide lock for that whole cascade would block
        every other concurrently running frame worker's own ``route()``
        call, even for frames belonging to a different checkpoint or grid
        entirely. Flushing is gated on an in-flight counter, not just the
        remaining-frame count, so that the *last* frame's ``route()`` call
        reaching ``_remaining[key] == 0`` can never flush the accumulator
        while an *earlier* frame's ``insert()`` for the same key is still
        running -- which would otherwise silently lose that frame's
        contribution to a ``finalize()``-reset accumulator.
        """
        frames = max(1, int(frames))
        index = self.checkpoint_index_for_frame(grid_name, frame_index)
        if frames > 1:
            last = self.checkpoint_index_for_frame(
                grid_name, int(frame_index) + frames - 1
            )
            if last != index:
                raise ValueError(
                    f"Frame group [{int(frame_index)}, "
                    f"{int(frame_index) + frames}) straddles grid "
                    f"{grid_name!r}'s planned checkpoint boundaries"
                )
        key = (grid_name, index)
        with self._lock:
            if key in self._resumed:
                return
            accumulator = self._accumulators.get(key)
            if accumulator is None:
                accumulator = _CheckpointAccumulator()
                self._accumulators[key] = accumulator
            self._inflight[key] = self._inflight.get(key, 0) + 1
            self._remaining[key] -= frames
        accumulator.insert(batch)
        to_flush = None
        with self._lock:
            self._inflight[key] -= 1
            self._frames_in_part[key] = (
                self._frames_in_part.get(key, 0) + frames
            )
            # Gate both flush triggers on no-other-insert-in-flight-for-this-
            # key, not just "done": an over-budget flush racing a still-
            # running sibling insert() for the same key would tear down
            # _accumulators[key]/_inflight[key] out from under it.
            no_inflight = self._inflight[key] <= 0
            done = no_inflight and self._remaining[key] <= 0
            over_budget = no_inflight and accumulator.should_flush(
                self._active_budget_bytes
            )
            if done or over_budget:
                del self._accumulators[key]
                del self._inflight[key]
                part = self._parts.get(key, 0)
                self._parts[key] = part + 1
                frames_covered = self._frames_in_part.pop(key, 0)
                to_flush = (accumulator, part, frames_covered)
        if to_flush is not None:
            accumulator, part, frames_covered = to_flush
            self._write_part(grid_name, index, part, accumulator, frames_covered)

    def _write_part(self, grid_name, index, part, accumulator, frames_covered):
        batch = accumulator.finalize()
        path = _checkpoint_part_path(self._checkpoint_dir, grid_name, index, part)
        _write_checkpoint(
            batch,
            path,
            spec_digest=self._spec_digest,
            metadata={
                "checkpoint_index": index,
                "grid_name": grid_name,
                "part": part,
                "frames_covered": frames_covered,
            },
        )
        with self._lock:
            self.written.append(path)


def _discover_checkpoint_state(
    checkpoint_dir, boundaries, spec_digest, *, cleanup_stale=False
):
    """Determine which planned checkpoints are already fully resumable.

    A planned ``(grid_name, index)`` checkpoint is resumable when its
    on-disk part files (Sec9/Sec10/Sec11) all carry a matching
    ``spec_sha256`` attribute and their ``frames_covered`` values sum to
    exactly its planned frame count -- the whole resume mechanism, with no
    separate manifest registry (design doc Sec11), extended one level to
    accommodate the Sec10 safety valve's occasional multi-part checkpoints.

    :param checkpoint_dir:
        Root scratch directory containing one subdirectory per grid.
    :param boundaries:
        ``dict[grid_name, list[(start, stop)]]`` -- the job's persisted
        checkpoint plan (design doc Sec5/Sec11: computed once, never
        recomputed on resume).
    :param str spec_digest:
        The job's resume-identity digest.
    :param bool cleanup_stale:
        If true, delete on-disk part files for any planned checkpoint
        found to be stale/partial/corrupt, so a following mapping run
        starts that checkpoint from a clean slate. Leave false when only
        inspecting status.
    :returns:
        ``(resumed, files)``: ``resumed`` is the set of fully-resumable
        ``(grid_name, index)`` pairs; ``files`` is
        ``dict[grid_name, list[Path]]`` of every valid checkpoint part
        file found, grouped by grid -- directly usable as
        :func:`_finalize_reconstruction`'s checkpoint-file input once
        every planned checkpoint is resumed.
    :rtype: tuple[set, dict]
    """
    checkpoint_dir = Path(checkpoint_dir)
    resumed = set()
    files: dict[str, list[Path]] = {}
    for grid_name, ranges in boundaries.items():
        grid_dir = checkpoint_dir / grid_name
        files[grid_name] = []
        by_index: dict[int, list[Path]] = {}
        if grid_dir.is_dir():
            for path in sorted(grid_dir.glob("ckpt*_p*.h5")):
                match = _CHECKPOINT_PART_PATTERN.match(path.name)
                if match is None:
                    continue
                by_index.setdefault(int(match.group("index")), []).append(path)
        for index, (start, stop) in enumerate(ranges):
            planned_frames = stop - start
            covered = 0
            valid_parts = []
            for path in by_index.get(index, ()):
                try:
                    with h5py.File(path, "r") as h5file:
                        if h5file.attrs.get("spec_sha256") != spec_digest:
                            continue
                        covered += int(h5file.attrs.get("frames_covered", 0))
                        valid_parts.append(path)
                except (OSError, KeyError):
                    continue
            if valid_parts and covered == planned_frames:
                resumed.add((grid_name, index))
                files[grid_name].extend(valid_parts)
            elif cleanup_stale:
                for path in by_index.get(index, ()):
                    try:
                        path.unlink()
                    except OSError:
                        pass
    return resumed, files


def _validate_mapping_setup(
    scan, detector, detector_tiles, angle_bounds_rad, correction_pipeline
) -> tuple[tuple[tuple[int, int, int, int], ...], np.ndarray]:
    """Validate a mapping run's fixed inputs once, up front.

    Design doc Sec7: the pipeline processes frames continuously through a
    shared prefetch/compute pool rather than in per-range calls, so this
    validation (previously repeated once per ``_map_frame_range`` call)
    now runs exactly once for the whole run instead of once per
    scheduling range.

    :param angle_bounds_rad:
        Shape ``(len(scan), 2, 4)`` in the order ``alpha, omega, chi,
        phi``, covering every frame in the scan (not just pending ones).
    :returns:
        ``(detector_tiles, bounds)`` -- ``detector_tiles`` normalized to a
        tuple of int tuples, ``bounds`` as a validated float array.
    """
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
    if bounds.shape != (len(scan), 2, 4):
        raise ValueError("angle_bounds_rad must have shape (len(scan), 2, 4)")
    if not np.all(np.isfinite(bounds)):
        raise ValueError("angle_bounds_rad contains non-finite values")
    if not callable(correction_pipeline):
        raise TypeError("correction_pipeline must be the central job correction")
    return detector_tiles, bounds


def _tile_ray_arrays(detector, detector_tiles, corner_rays=None):
    """Build (or reuse cached) corner rays for every detector tile, once.

    :param corner_rays:
        Optional ``{detector_tile: array}`` cache (design doc Sec7: shared
        read-only across every compute worker, built once per mapping run
        rather than once per frame or per worker).
    """
    ray_arrays = {}
    for detector_tile in detector_tiles:
        row_start, row_stop, column_start, column_stop = detector_tile
        rays = (
            _detector_corner_rays(detector, detector_tile)
            if corner_rays is None or detector_tile not in corner_rays
            else np.ascontiguousarray(corner_rays[detector_tile], dtype=np.float64)
        )
        expected_ray_shape = (
            row_stop - row_start + 1,
            column_stop - column_start + 1,
            3,
        )
        if rays.shape != expected_ray_shape:
            raise ValueError(f"corner_rays must have shape {expected_ray_shape}")
        ray_arrays[detector_tile] = rays
    return ray_arrays


def _build_kernels(spec, ub_calculator, *, threads=None, memory_budget_bytes=None):
    """One kernel per grid, for one compute worker's exclusive use.

    Each compute worker in the Sec7 pipeline builds and owns its own
    kernel set (kernels are not safe to call concurrently from multiple
    Python threads -- each call spawns its own native worker threads
    already, per ``kernel_threads``); this is the shared construction
    logic every worker calls once at start-up.
    """
    return {
        grid.grid_name: _kernel_for_grid(
            spec,
            grid,
            ub_calculator,
            threads=threads,
            memory_budget_bytes=memory_budget_bytes,
        )
        for grid in spec.grids
    }


def _map_frame_group(
    spec: _ReconstructionSpec,
    kernels: Mapping[str, object],
    ray_arrays: Mapping[tuple[int, int, int, int], np.ndarray],
    detector_tiles: Iterable[tuple[int, int, int, int]],
    correction_pipeline: Callable,
    image_payloads: Sequence[object],
    frame_indices: Sequence[int],
    angles_start: np.ndarray,
    angles_end: np.ndarray,
    router: _CheckpointRouter,
    *,
    corrected_frames: Sequence[tuple] | None = None,
) -> None:
    """Process one already-loaded group of frames: correction, per-tile-
    per-grid kernel accumulate, merge, and route.

    This is the compute-side half of the mapping pipeline (design doc
    Sec7): image loading (I/O) happens separately, in the prefetch
    pipeline, so ``image_payloads`` here are already loaded -- this
    function does no I/O of its own and is safe to call concurrently for
    different groups from multiple compute workers, provided each worker
    uses its own ``kernels`` (see :func:`_build_kernels`) and every
    worker shares the same read-only ``ray_arrays``.

    ``frame_indices`` must be contiguous and ascending, and must lie
    inside one planned checkpoint range for every grid: the group's
    frames merge inside the kernel and cannot be told apart afterwards,
    so they must all belong to the same checkpoint. The caller forms
    groups under those rules; :meth:`_CheckpointRouter.route` re-checks
    the checkpoint one.

    A single-frame group is the degenerate case of the same call, not a
    separate path -- but note it is *not* bit-for-bit with a per-image
    ``accumulate()``: a brick partitions the detector into rectangles
    where a block partitions it into runs of the flattened image, so
    pixels group into accumulators differently and their sums associate
    differently. Same voxels, same contributor counts, totals equal to
    rounding, and about 3.4% fewer records from the better block shape.

    The correction pipeline's whole-frame ``correct_frame`` step (pixel
    repair, static factors) runs here by default, in the compute worker.
    Pass ``corrected_frames`` to have it run somewhere else instead: the
    grouped scheduler corrects in its prepare pool, so that the GIL-held
    Python work of group N+1 overlaps group N's native call rather than
    queueing behind it. A pipeline with no ``correct_frame`` corrects per
    tile here either way, since there is nothing whole-frame to hoist.

    :param image_payloads:
        Loaded payloads, one per frame in the group, in frame order. May
        be ``None`` when ``corrected_frames`` supplies the arrays and the
        pipeline has a ``correct_frame`` step, since the raw images are
        then no longer needed and are better released.
    :param frame_indices:
        The group's frame indices: contiguous, ascending, all inside one
        checkpoint range per grid.
    :param corrected_frames:
        Optional pre-corrected ``(intensity, variance, mask)`` triples,
        one per frame, each full-detector sized.
    :param angles_start:
        ``(frames, 4)`` exposure-start diffractometer angles, radians.
    :param angles_end:
        ``(frames, 4)`` exposure-end angles, radians. Equal to
        ``angles_start`` for a stationary exposure; the kernel decides
        per frame, so a group may mix the two.
    """
    frame_indices = [int(index) for index in frame_indices]
    images = (
        [np.asarray(payload.img) for payload in image_payloads]
        if image_payloads is not None
        else [None] * len(frame_indices)
    )
    frame_correction = getattr(correction_pipeline, "correct_frame", None)
    # The whole group's corrected frames stay resident until every tile
    # has been mapped: a tile is a row band of the detector, so each
    # frame is revisited once per band. This is the group buffer, and it
    # is what _frame_parallelism sizes a compute worker against.
    if corrected_frames is None and callable(frame_correction):
        corrected_frames = [
            frame_correction(payload, image, frame_index)
            for payload, image, frame_index in zip(
                image_payloads, images, frame_indices
            )
        ]
    # Collect every tile's contribution per grid, then route exactly once
    # per (group, grid) -- the router's remaining-frame countdown
    # (Sec9/Sec10) counts frames, so a multi-tile group must be merged
    # down to one routed batch per grid before reaching it.
    tile_batches: dict[str, list[Mapping[str, np.ndarray]]] = {
        grid.grid_name: [] for grid in spec.grids
    }
    for detector_tile in detector_tiles:
        row_start, row_stop, column_start, column_stop = detector_tile
        selection = np.s_[row_start:row_stop, column_start:column_stop]
        if corrected_frames is None:
            tile_values = [
                correction_pipeline(
                    payload,
                    image[selection],
                    frame_index,
                    detector_tile,
                )
                for payload, image, frame_index in zip(
                    image_payloads, images, frame_indices
                )
            ]
        else:
            tile_values = [
                tuple(values[selection] for values in corrected)
                for corrected in corrected_frames
            ]
        # One contiguous (frames, rows, columns) buffer per array, which
        # is the copy the per-tile ascontiguousarray already made for a
        # single frame -- the group spends it over F frames instead of 1.
        intensity = np.ascontiguousarray(
            np.stack([values[0] for values in tile_values]), dtype=np.float64
        )
        variance = np.ascontiguousarray(
            np.stack([values[1] for values in tile_values]), dtype=np.float64
        )
        mask = np.ascontiguousarray(
            np.stack([values[2] for values in tile_values]), dtype=bool
        )
        del tile_values
        for grid in spec.grids:
            batch = kernels[grid.grid_name].accumulate_group(
                intensity,
                variance,
                mask,
                ray_arrays[detector_tile],
                angles_start,
                angles_end,
            )
            tile_batches[grid.grid_name].append(batch)
    for grid in spec.grids:
        merged = _reduce_batches(tile_batches[grid.grid_name])
        router.route(
            grid.grid_name,
            frame_indices[0],
            merged,
            frames=len(frame_indices),
        )


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


class _CheckpointRangeReader:
    """Stream monotonically increasing ``(chunk_id, local_voxel_id)`` ranges
    from one checkpoint HDF5 file, without loading it into memory at once.

    Structurally mirrors the retired Parquet-era range reader: a
    forward-only cursor over sorted rows, refilled in windows of
    ``batch_rows``. Checkpoint files are internally sorted by
    ``(chunk_id, local_voxel_id)`` (Sec8/Sec9's tree-merge guarantee), so
    no per-row-group statistics are needed to skip ahead -- a plain
    forward scan is sufficient, and :meth:`read` requires callers to
    request ranges in non-decreasing order, exactly as the old reader did.
    """

    def __init__(self, path, *, batch_rows=131072):
        self._file = h5py.File(path, "r")
        self._total_rows = int(self._file["chunk_id"].shape[0])
        self._batch_rows = max(1, int(batch_rows))
        self._cursor = 0
        self.batch = None
        self.position = 0
        self.previous_key = None

    def close(self):
        """Close the underlying HDF5 file."""
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def _advance(self):
        if self._cursor >= self._total_rows:
            self.batch = None
            return False
        stop = min(self._total_rows, self._cursor + self._batch_rows)
        self.batch = {
            name: self._file[name][self._cursor : stop] for name in _PARTIAL_COLUMNS
        }
        self._cursor = stop
        self.position = 0
        return True

    def read(self, chunk_id, local_start, local_stop):
        """Read one exact range without revisiting earlier checkpoint data."""
        key = (int(chunk_id), int(local_start))
        if self.previous_key is not None and key < self.previous_key:
            raise ValueError("Checkpoint ranges must be requested in sorted order")
        self.previous_key = key
        parts = []
        while self.batch is not None or self._advance():
            chunks = self.batch["chunk_id"]
            position = self.position
            if position >= chunks.size:
                self.batch = None
                continue
            position += int(
                np.searchsorted(chunks[position:], np.uint32(chunk_id), side="left")
            )
            if position >= chunks.size:
                self.batch = None
                continue
            current_chunk = int(chunks[position])
            if current_chunk > chunk_id:
                self.position = position
                break
            chunk_stop = int(
                np.searchsorted(chunks, np.uint32(chunk_id), side="right")
            )
            local = self.batch["local_voxel_id"]
            start = position + int(
                np.searchsorted(
                    local[position:chunk_stop], np.uint32(local_start), side="left"
                )
            )
            stop = position + int(
                np.searchsorted(
                    local[position:chunk_stop], np.uint32(local_stop), side="left"
                )
            )
            if stop > start:
                parts.append(
                    {name: values[start:stop] for name, values in self.batch.items()}
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


def _write_chunk_from_checkpoints(group, grid, readers, chunk_id):
    """Merge every checkpoint's contribution to one chunk and write it.

    Different checkpoints can legitimately contain overlapping
    ``(chunk_id, local_voxel_id)`` records -- unlike the retired
    reduce-phase Parquet shards, which were already deduplicated against
    each other before finalize ran, checkpoints are only deduplicated
    within their own frame range (Sec9). The Sec9 tree-merge is reused
    here, across checkpoints instead of across frames, to combine them
    before writing.
    """
    coordinates = _chunk_coordinates(chunk_id, grid)
    starts = tuple(
        coordinate * chunk for coordinate, chunk in zip(coordinates, grid.chunk_shape)
    )
    stops = tuple(
        min(start + chunk, size)
        for start, chunk, size in zip(starts, grid.chunk_shape, grid.shape)
    )
    local_shape = tuple(stop - start for start, stop in zip(starts, stops))
    local_stop = (
        (local_shape[0] - 1) * grid.chunk_shape[1] * grid.chunk_shape[2]
        + (local_shape[1] - 1) * grid.chunk_shape[2]
        + local_shape[2]
    )
    levels = []
    for reader in readers:
        batch = reader.read(chunk_id, 0, local_stop)
        if batch["chunk_id"].size:
            _tree_insert(levels, batch)
    reduced = _tree_finalize(levels)
    if not reduced["chunk_id"].size:
        return
    if np.any(reduced["chunk_id"] != np.uint32(chunk_id)):
        raise ValueError(f"Checkpoint data for chunk {chunk_id} is inconsistent")
    local = reduced["local_voxel_id"].astype(np.uint32, copy=False)
    local_x, remainder = np.divmod(
        local, np.uint32(grid.chunk_shape[1] * grid.chunk_shape[2])
    )
    local_y, local_z = np.divmod(remainder, np.uint32(grid.chunk_shape[2]))
    valid = (
        (local_x < local_shape[0])
        & (local_y < local_shape[1])
        & (local_z < local_shape[2])
    )
    if not np.all(valid):
        raise ValueError(f"Chunk {chunk_id} contains an invalid local voxel ID")
    index = (
        local_x.astype(int),
        local_y.astype(int),
        local_z.astype(int),
    )
    weighted_intensity = np.zeros(local_shape, dtype=np.float64)
    weighted_variance = np.zeros(local_shape, dtype=np.float64)
    weight = np.zeros(local_shape, dtype=np.float64)
    contributors = np.zeros(local_shape, dtype=np.uint32)
    weighted_intensity[index] = reduced["weighted_intensity"]
    weighted_variance[index] = reduced["weighted_variance"]
    weight[index] = reduced["weight"]
    contributors[index] = reduced["contributors"]
    populated = weight > 0
    weighted_intensity[~populated] = np.nan
    weighted_variance[~populated] = np.nan
    # Divides variance by weight twice (i.e. by weight**2), matching the
    # pre-existing reduce/finalize behavior this rewrite must not silently
    # change (AGENTS.md: scientific conventions are not this refactor's to
    # alter).
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
    spec: _ReconstructionSpec,
    checkpoint_files: Mapping[str, Iterable],
    output_path,
    *,
    provenance: Mapping[str, Any] | None = None,
    scientific_context: Mapping[str, Any] | None = None,
    config=None,
    chunk_progress: Callable | None = None,
):
    """Create the final HDF5 reconstruction directly from checkpoint files.

    Replaces the retired reduce-then-finalize two-step: each grid's output
    chunks are streamed and merged directly from its checkpoint part files
    (Sec9's tree-merge, reused across checkpoints instead of across
    frames), with no intermediate reduced-shard files. Memory-bounded the
    same way mapping was -- at most one window's worth of data from each
    reader is ever materialized at a time
    (:class:`_CheckpointRangeReader`), never a whole grid's total data
    volume at once.

    :param _ReconstructionSpec spec:
        The job's specification (grids, compression, ...). Unlike the
        retired manifest-based path, this is supplied directly by the
        caller rather than reconstructed from scratch-file metadata.
    :param checkpoint_files:
        ``dict[grid_name, Iterable[path]]`` -- every checkpoint part file
        to merge for each grid, e.g. from
        :func:`_discover_checkpoint_state`.
    """
    grid_by_name = {grid.grid_name: grid for grid in spec.grids}
    unknown = set(checkpoint_files) - set(grid_by_name)
    if unknown:
        raise ValueError(f"Unknown grid(s) in checkpoint files: {sorted(unknown)}")
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
        if scientific_context:
            process.create_dataset(
                "scientific_context_json",
                data=json.dumps(
                    dict(scientific_context), sort_keys=True, default=str
                ),
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
                dtype=np.uint32,
                chunks=hdf5_chunks,
                fillvalue=0,
                **compression,
            )
            groups[grid.grid_name] = group
        per_grid_chunks = {}
        for grid in spec.grids:
            paths = [Path(path) for path in checkpoint_files.get(grid.grid_name, ())]
            populated: set[int] = set()
            for path in paths:
                with h5py.File(path, "r") as source:
                    populated.update(
                        int(value) for value in np.unique(source["chunk_id"][()])
                    )
            per_grid_chunks[grid.grid_name] = (paths, sorted(populated))
        total_chunks = sum(len(chunks) for _paths, chunks in per_grid_chunks.values())
        written = 0
        for grid in spec.grids:
            paths, chunk_ids = per_grid_chunks[grid.grid_name]
            readers = []
            try:
                for path in paths:
                    readers.append(_CheckpointRangeReader(path))
                for chunk_id in chunk_ids:
                    _write_chunk_from_checkpoints(
                        groups[grid.grid_name], grid, readers, chunk_id
                    )
                    written += 1
                    if chunk_progress is not None:
                        chunk_progress(written, total_chunks)
            finally:
                for reader in readers:
                    reader.close()
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
    return {
        "path": str(output_path),
        "sha256": digest.hexdigest(),
        "bytes": output_path.stat().st_size,
        "chunks_written": total_chunks,
    }


__all__ = []

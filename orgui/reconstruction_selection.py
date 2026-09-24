# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
"""Automatic reciprocal-space volume selection for the reconstruction.

Instead of one large grid covering the whole measured volume, a scan is often
better reconstructed as many small grids, each holding one crystallographic
feature:

- a *crystal truncation rod* is a column along ``L`` at one allowed integer
  ``(H, K)``, and
- a *fractional-order rod* is a column along ``L`` at an explicitly selected
  fractional ``(H, K)``, and
- a *Bragg reflection* is a compact box around one allowed integer
  ``(H, K, L)``, and
- a *fractional-order Bragg peak* is a compact box around an explicitly
  selected noninteger ``(H, K, L)``.

Integer CTR and Bragg selections enumerate features from the reference unit
cell; fractional rods and peaks use explicit index rules. All keep only the
ones the active scan actually reaches and return one
:class:`~orgui.reconstruction_job.ReconstructionGrid` per feature. All grids of
one selection share a frame, a voxel step, and therefore units, so the
reconstruction extracts them together in a single pass over the images -- see
``_map_frame_group`` in
:mod:`orgui.datautils.xrayutils.reconstruction`, which maps every grid of the
job from the same corrected frame group. Reading and correcting the images is
therefore paid once, however many features were selected.

The per-grid cost is not zero, though. Each grid gets its own output group and
its own checkpoint stream, ``split_memory_budget`` divides the checkpoint
memory between the grids, and the detector tile the mapping loop uses shrinks
as the grid count grows. A few dozen features is routine; several hundred is
not, which is what :data:`DEFAULT_MAX_GRIDS` guards against.

Reciprocal-lattice coordinates ``(h, k, l)`` are in r.l.u., momentum transfer
in ``Angstrom^-1``, and all diffractometer angles handled here are in radians.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import lcm
import re

import numpy as np

from .datautils.xrayutils import ReciprocalNavigation
from .datautils.xrayutils.reconstruction import _sample_angle_bounds
from .reconstruction_job import ReconstructionGrid


#: Output frames whose relation to ``hkl`` is a single fixed matrix, and which
#: automatic selection can therefore place a feature box in. ``lab``,
#: ``alpha``, ``omega``, ``chi``, and ``phi`` rotate with the sample during the
#: scan, so a fixed ``(H, K, L)`` feature has no fixed box in them.
STATIC_FRAMES = ("hkl", "crystal")

#: Refuse to build more grids than this unless the caller raises the limit.
#: Each grid is a separate HDF5 group with its own checkpoint stream, so an
#: accidental selection of thousands of features is a job that will not finish.
DEFAULT_MAX_GRIDS = 512


def _normalized_frame(frame):
    return str(frame).removeprefix("q_").lower()


def hkl_to_frame_matrix(ub_calculator, frame="hkl"):
    """Return the fixed linear map from r.l.u. to one output frame.

    The reconstruction kernel derives ``hkl`` from a lab-frame scattering
    vector as ``UB^-1 q`` and the crystal frame as ``U^-1 q``, so
    ``q_crystal = U^-1 UB h`` is constant for the whole scan, while the
    remaining frames still carry the per-frame sample rotation.

    :param ub_calculator:
        Active :class:`~orgui.datautils.xrayutils.HKLVlieg.UBCalculator`.
    :param str frame:
        ``hkl`` or ``crystal`` (``q_crystal`` is accepted as a synonym).
    :returns:
        ``(3, 3)`` matrix mapping r.l.u. to r.l.u. (``hkl``) or to
        ``Angstrom^-1`` (``crystal``).
    :rtype: numpy.ndarray
    :raises ValueError:
        If the frame is not one of :data:`STATIC_FRAMES`.
    """
    normalized = _normalized_frame(frame)
    if normalized == "hkl":
        return np.eye(3)
    if normalized == "crystal":
        return np.ascontiguousarray(
            np.linalg.inv(ub_calculator.getU()) @ ub_calculator.getUB()
        )
    raise ValueError(
        "Automatic volume selection needs a frame that does not rotate with "
        f"the sample, one of {', '.join(STATIC_FRAMES)}; got {frame!r}"
    )


def sample_hkl_coverage(
    config,
    scan,
    *,
    detector_samples=33,
    frame_samples=128,
):
    """Sample the reciprocal-space positions the scan actually reaches.

    A flattened view of :func:`sample_hkl_coverage_by_frame`; see there for
    the sampling itself and for the per-frame form, which is what a caller
    asking *which frames* reach a feature needs.

    The returned cloud is a *sample* of the swept detector surface, not its
    exact hull: a feature touched only by a few detector pixels of a few frames
    can fall between samples and then be missed by the selection functions
    below. Raise ``detector_samples``/``frame_samples`` when that matters.

    Frames listed in ``config.corrections.excluded_frames`` are left out, and
    both exposure endpoints are sampled so that a rocking exposure contributes
    its whole sweep.

    :param config:
        Central :class:`~orgui.app.config_data.ConfigData` experiment snapshot.
    :param scan:
        Active scan backend providing exposure angle bounds in radians.
    :param int detector_samples:
        Pixel centers sampled along each detector axis; at least two.
    :param int frame_samples:
        Maximum number of included scan frames sampled; at least one.
    :returns:
        ``(n, 3)`` sampled reciprocal-lattice coordinates in r.l.u.
    :rtype: numpy.ndarray
    :raises ValueError:
        If the sampling counts are invalid or no frame is included.
    """
    _frames, coordinates = sample_hkl_coverage_by_frame(
        config,
        scan,
        detector_samples=detector_samples,
        frame_samples=frame_samples,
    )
    return coordinates.reshape(-1, 3)


def sample_hkl_coverage_by_frame(
    config,
    scan,
    *,
    detector_samples=33,
    frame_samples=128,
    frames=None,
):
    """Sample the reciprocal space the scan reaches, keeping the frame.

    Same sampling as :func:`sample_hkl_coverage`, but the samples stay
    grouped by the frame they came from, so a caller can ask which frames
    reach a given volume rather than only whether any of them does. That
    distinction is what separates a feature the whole scan sweeps through
    from one a handful of frames graze: a crystal truncation rod at a fixed
    ``(H, K)`` is crossed over a narrow range of sample rotations, while the
    specular rod is reached by every frame.

    :param config:
        Central :class:`~orgui.app.config_data.ConfigData` experiment
        snapshot.
    :param scan:
        Active scan backend providing exposure angle bounds in radians.
    :param int detector_samples:
        Pixel centers sampled along each detector axis; at least two.
    :param int frame_samples:
        Maximum number of included scan frames sampled; at least one.
    :param frames:
        Optional explicit frame indices to draw the sample from, instead of
        every frame the config includes. A caller estimating one cluster
        node's own slice needs the sample drawn from that slice: spread over
        the whole scan, a few dozen samples may land in a short range not at
        all, and a volume is then scored against frames the node will never
        map.
    :returns:
        ``(frames, coordinates)`` -- the sampled frame indices, shape
        ``(f,)``, and their reciprocal-lattice coordinates in r.l.u., shape
        ``(f, n, 3)`` with ``n`` samples per frame.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    :raises ValueError:
        If the sampling counts are invalid or no frame is included.
    """
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
    if bounds.shape not in {(len(scan), 2, 4), (len(scan), 2, 6)}:
        raise ValueError(
            "Exposure angle bounds must have shape (frames, 2, 4) or "
            "(frames, 2, 6)"
        )
    bounds = _sample_angle_bounds(bounds)
    if frames is None:
        excluded = set(config.corrections.excluded_frames)
        included = np.asarray(
            [index for index in range(len(scan)) if index not in excluded],
            dtype=np.int64,
        )
    else:
        included = np.unique(np.asarray(frames, dtype=np.int64))
    if included.size == 0:
        raise ValueError(
            "No included frames are available for coverage sampling"
        )

    sampled_rows = np.unique(
        np.rint(np.linspace(0, rows - 1, min(detector_samples, rows)))
    )
    sampled_columns = np.unique(
        np.rint(np.linspace(0, columns - 1, min(detector_samples, columns)))
    )
    sampled_frames = included[
        np.unique(
            np.rint(
                np.linspace(
                    0, included.size - 1, min(frame_samples, included.size)
                )
            ).astype(np.int64)
        )
    ]

    # One ray per sampled pixel center, built through the same detector
    # geometry path estimate_geometry_steps uses.
    row_grid, column_grid = np.meshgrid(
        sampled_rows, sampled_columns, indexing="ij"
    )
    gamma, delta = detector.primBeamPoints(row_grid, column_grid)
    cosine_gamma = np.cos(gamma)
    pixel_rays = np.empty(np.shape(gamma) + (3,), dtype=np.float64)
    pixel_rays[..., 0] = np.sin(delta) * cosine_gamma
    pixel_rays[..., 1] = np.cos(delta) * cosine_gamma
    pixel_rays[..., 2] = np.sin(gamma)
    pixel_rays /= np.linalg.norm(pixel_rays, axis=-1, keepdims=True)

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
        "hkl",
        float(ub.getK()),
        np.ascontiguousarray(np.linalg.inv(ub.getUB())),
        np.ascontiguousarray(np.linalg.inv(ub.getU())),
        0,
        1,
        1,
        1024 * 1024,
    )

    coordinates = []
    for frame_index in sampled_frames:
        angles_start = np.ascontiguousarray(bounds[frame_index, 0])
        angles_end = np.ascontiguousarray(bounds[frame_index, 1])
        per_frame = []
        for row_index in range(pixel_rays.shape[0]):
            for column_index in range(pixel_rays.shape[1]):
                # A degenerate 2x2 corner patch of the one pixel-center ray:
                # the kernel interpolates over (u, v), so repeating the ray
                # makes every (u, v) the pixel center and leaves only the
                # exposure parameter t free.
                rays = np.ascontiguousarray(
                    np.broadcast_to(
                        pixel_rays[row_index, column_index], (2, 2, 3)
                    )
                )
                for t in (0.0, 1.0):
                    per_frame.append(
                        kernel.coordinate(
                            rays, angles_start, angles_end, 0, 0, 0.0, 0.0, t
                        )
                    )
        coordinates.append(per_frame)
    return (
        np.asarray(sampled_frames, dtype=np.int64),
        np.asarray(coordinates, dtype=np.float64),
    )


def _bulk_unit_cell(xtal):
    """Return the bulk unit cell backing a crystal or unit-cell object."""
    return getattr(xtal, "uc_bulk", xtal)


def _format_index(value):
    return str(int(round(float(value)))).replace("-", "m")


def _format_fractional_index(value):
    fraction = Fraction(float(value)).limit_denominator(24)
    if fraction.denominator == 1:
        return _format_index(value)
    return f"{fraction.numerator}p{fraction.denominator}".replace("-", "m")


def _fractional_condition(source, axes):
    """Parse a half-order or integer-parity index condition."""
    source = source.strip().lower()
    if source in ("any half", "all half"):
        return (source,)
    match = re.fullmatch(
        rf"([{axes}](?:\s*\+\s*[{axes}])*)\s*(?:=\s*|\s+)(even|odd)",
        source,
    )
    if match:
        indices = tuple(axes.index(axis) for axis in match.group(1)
                        if axis in axes)
        if len(indices) != len(set(indices)):
            raise ValueError(f"Repeated index in fractional rule {source!r}")
        return (match.group(2), indices)
    raise ValueError(f"Unsupported fractional condition {source!r}")


def _fractional_pattern(source, axes):
    """Parse index phases modulo one; ``*`` matches every phase."""
    parts = [part.strip() for part in source.split(",")]
    if len(parts) != len(axes):
        raise ValueError(f"A fractional pattern needs {','.join(axes).upper()} parts")
    pattern = []
    for part in parts:
        if part == "*":
            pattern.append(None)
            continue
        try:
            phase = Fraction(part)
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError(f"Invalid fractional phase {part!r}") from error
        if not 0 <= phase < 1 or phase.denominator > 24:
            raise ValueError(
                "Fractional phases must be in [0, 1) with denominator <= 24"
            )
        pattern.append(phase)
    return tuple(pattern)


def _fractional_rules(families, exclusions, axes):
    """Parse semicolon-separated existence and exclusion descriptors."""
    include = []
    for entry in str(families or "").split(";"):
        if not entry.strip():
            continue
        parts = re.split(r"\s+where\s+", entry.strip(), maxsplit=1, flags=re.I)
        if parts[0].lower() in ("any half", "all half"):
            pattern = (None,) * len(axes)
            conditions = [_fractional_condition(parts[0], axes)]
        else:
            pattern = _fractional_pattern(parts[0], axes)
            if not any(phase not in (None, 0) for phase in pattern):
                raise ValueError("A fractional family needs a nonzero phase")
            conditions = []
        if len(parts) == 2:
            conditions.append(_fractional_condition(parts[1], axes))
        include.append((pattern, conditions))
    if not include:
        raise ValueError("Specify at least one fractional family")
    exclude = []
    for entry in str(exclusions or "").split(";"):
        if not entry.strip():
            continue
        exclude.append(
            ("pattern", _fractional_pattern(entry, axes)) if "," in entry
            else ("condition", _fractional_condition(entry, axes))
        )
    return include, exclude


def _fractional_condition_matches(numerators, condition, denominator):
    if condition[0] in ("any half", "all half"):
        half = [value % denominator == denominator // 2 for value in numerators]
        return any(half) if condition[0] == "any half" else all(half)
    index_sum = sum(numerators[axis] for axis in condition[1])
    target = 0 if condition[0] == "even" else denominator
    return index_sum % (2 * denominator) == target


def _fractional_pattern_matches(numerators, pattern, denominator):
    return all(
        phase is None
        or numerators[axis] % denominator
        == phase.numerator * (denominator // phase.denominator)
        for axis, phase in enumerate(pattern)
    )


def _fractional_indices(families, exclusions, limits, axes):
    """Return sorted noninteger index tuples matching the rules."""
    include, exclude = _fractional_rules(families, exclusions, axes)
    family_denominator = lcm(*(
        [phase.denominator for pattern, _ in include for phase in pattern
         if phase is not None]
        + [2 for _, conditions in include for condition in conditions
           if condition[0] in ("any half", "all half")]
        + [1]
    ))
    denominator = lcm(*(
        [family_denominator]
        + [phase.denominator for kind, pattern in exclude
           if kind == "pattern" for phase in pattern if phase is not None]
        + [2 for kind, condition in exclude if kind == "condition"
           and condition[0] in ("any half", "all half")]
    ))
    if family_denominator > 24 or denominator > 120:
        raise ValueError("Combined fractional denominator is too large")
    candidates = set()
    for pattern, conditions in include:
        candidate_axes = []
        for phase, (lower, upper) in zip(pattern, limits):
            residues = (
                range(family_denominator) if phase is None
                else (phase.numerator * (family_denominator // phase.denominator),)
            )
            candidate_axes.append([
                value for value in range(lower * family_denominator,
                                         upper * family_denominator + 1)
                if value % family_denominator in residues
            ])
        if len(candidates) + np.prod([len(axis) for axis in candidate_axes]) > 100_000:
            raise ValueError("Fractional rules enumerate too many candidates")
        candidates.update(
            values for values in product(*candidate_axes)
            if any(value % family_denominator for value in values)
            and all(_fractional_condition_matches(values, condition, family_denominator)
                    for condition in conditions)
        )
    selected = []
    for values in sorted(candidates):
        numerators = tuple(
            value * (denominator // family_denominator) for value in values
        )
        if any(
            _fractional_pattern_matches(numerators, rule, denominator)
            if kind == "pattern"
            else _fractional_condition_matches(numerators, rule, denominator)
            for kind, rule in exclude
        ):
            continue
        selected.append(tuple(value / family_denominator for value in values))
    return selected


def _axis_integers(limits, coverage_lower, coverage_upper):
    """Integer indices for one axis, from explicit limits or from coverage.

    A scalar ``n`` is read as the symmetric ``(-n, n)``; a pair is read as
    ``(lower, upper)``; ``None`` spans the measured coverage.
    """
    if limits is None:
        return np.arange(
            int(np.floor(coverage_lower)),
            int(np.ceil(coverage_upper)) + 1,
            dtype=np.int64,
        )
    values = np.atleast_1d(np.asarray(limits)).reshape(-1)
    if values.size == 1:
        magnitude = abs(int(values[0]))
        lower, upper = -magnitude, magnitude
    elif values.size == 2:
        lower, upper = int(values[0]), int(values[1])
    else:
        raise ValueError(
            "Axis limits must be a single symmetric magnitude or a "
            "(lower, upper) pair"
        )
    if upper < lower:
        raise ValueError("Axis limits must be ordered as (lower, upper)")
    return np.arange(lower, upper + 1, dtype=np.int64)


def _feature_box(matrix, corners_hkl, half_width):
    """Axis-aligned output-frame box around one feature.

    ``corners_hkl`` is ``(m, 3)`` in r.l.u.: one point for a Bragg reflection,
    the two segment endpoints for a rod. The box is the bounding box of those
    points in the output frame, expanded by ``half_width`` on each axis. For a
    rod in a non-orthogonal cell this follows the rod at the cost of covering a
    little more than the rod itself.
    """
    transformed = np.asarray(corners_hkl, dtype=np.float64) @ np.asarray(
        matrix, dtype=np.float64
    ).T
    return (
        np.min(transformed, axis=0) - half_width,
        np.max(transformed, axis=0) + half_width,
    )


def _feature_grid(minimum, maximum, step, frame, name, chunk_shape):
    """Turn one output-frame box into a :class:`ReconstructionGrid`."""
    step = np.asarray(step, dtype=np.float64).reshape(-1)
    if step.size != 3:
        raise ValueError("step must contain exactly three values")
    if np.any(~np.isfinite(step)) or np.any(step <= 0):
        raise ValueError("Voxel steps must be finite and positive")
    if np.any(maximum <= minimum):
        raise ValueError(
            f"Feature {name!r} has no extent on every axis; give it "
            "non-zero half-widths"
        )
    return ReconstructionGrid(
        minimum=tuple(float(value) for value in minimum),
        maximum=tuple(float(value) for value in maximum),
        step=tuple(float(value) for value in step),
        frame=_normalized_frame(frame),
        name=name,
        chunk_shape=tuple(int(value) for value in chunk_shape),
    )


def _check_grid_count(grids, max_grids, kind):
    if max_grids is not None and len(grids) > int(max_grids):
        raise ValueError(
            f"Automatic {kind} selection produced {len(grids)} grids, above "
            f"the limit of {int(max_grids)}. Narrow the H/K limits, coarsen "
            "the steps, or raise max_grids deliberately -- every grid costs "
            "its own checkpoint stream and output group."
        )
    return grids


def _prepare_coverage(config, scan, coverage, detector_samples, frame_samples):
    if coverage is None:
        coverage = sample_hkl_coverage(
            config,
            scan,
            detector_samples=detector_samples,
            frame_samples=frame_samples,
        )
    coverage = np.asarray(coverage, dtype=np.float64)
    if coverage.ndim != 2 or coverage.shape[1] != 3:
        raise ValueError("coverage must have shape (n, 3)")
    if coverage.shape[0] == 0:
        raise ValueError("The coverage sample is empty")
    return coverage


def derive_ctr_grids(
    config,
    scan,
    *,
    step,
    half_width=(0.05, 0.05, 0.0),
    frame="hkl",
    h_limits=None,
    k_limits=None,
    coverage=None,
    structure_factor_samples=21,
    detector_samples=33,
    frame_samples=128,
    chunk_shape=(64, 64, 64),
    max_grids=DEFAULT_MAX_GRIDS,
):
    """One output grid per crystal truncation rod the scan reaches.

    A rod is the column at an allowed integer ``(H, K)`` running along ``L``
    over the measured ``L`` range. Candidate ``(H, K)`` are enumerated from the
    reference unit cell, keeping those whose unit-cell structure factor is
    non-zero somewhere in the measured ``L`` range, and are then kept only if a
    sampled coverage point falls inside the rod's own box.

    :param config:
        Central :class:`~orgui.app.config_data.ConfigData` experiment snapshot.
    :param scan:
        Active scan backend providing exposure angle bounds in radians.
    :param step:
        Three voxel widths, in r.l.u. for ``hkl`` and ``Angstrom^-1`` for
        ``crystal``.
    :param half_width:
        Padding added on each side of the rod, on each output-frame axis, in
        the same units as ``step``. For ``hkl`` the first two entries are the
        ``H`` and ``K`` half-widths of the column and the third extends the
        measured ``L`` range, which the default leaves as measured. A single
        number is read as the same padding on all three axes.
    :param str frame:
        ``hkl`` or ``crystal``; see :data:`STATIC_FRAMES`.
    :param h_limits:
        Optional ``(lower, upper)`` integer ``H`` limits in r.l.u. A single
        number ``n`` is read as the symmetric ``(-n, n)``. ``None`` takes the
        limits from the measured coverage.
    :param k_limits:
        As ``h_limits``, for ``K``.
    :param coverage:
        Optional pre-computed ``(n, 3)`` r.l.u. coverage cloud from
        :func:`sample_hkl_coverage`, reused instead of sampling again.
    :param int structure_factor_samples:
        Number of ``L`` values at which the unit-cell structure factor is
        probed when deciding whether a rod is allowed.
    :param int detector_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param int frame_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param chunk_shape:
        HDF5 chunk shape in voxels, shared by every returned grid.
    :param max_grids:
        Refuse to return more than this many grids; ``None`` disables the
        guard.
    :returns:
        Rod grids ordered by ``(H, K)``, named ``ctr_<H>_<K>`` with ``m``
        standing in for a minus sign.
    :rtype: list[orgui.reconstruction_job.ReconstructionGrid]
    :raises ValueError:
        If the frame rotates with the sample, the steps are not positive, or
        the selection exceeds ``max_grids``.
    """
    matrix = hkl_to_frame_matrix(config.ub_calculator, frame)
    coverage = _prepare_coverage(
        config, scan, coverage, detector_samples, frame_samples
    )
    coverage_minimum = np.min(coverage, axis=0)
    coverage_maximum = np.max(coverage, axis=0)

    h_values = _axis_integers(
        h_limits, coverage_minimum[0], coverage_maximum[0]
    )
    k_values = _axis_integers(
        k_limits, coverage_minimum[1], coverage_maximum[1]
    )
    l_lower = float(coverage_minimum[2])
    l_upper = float(coverage_maximum[2])
    l_probe = np.linspace(
        l_lower, l_upper, max(2, int(structure_factor_samples))
    )

    allowed = np.atleast_2d(
        ReciprocalNavigation.allowedCTRs(
            _bulk_unit_cell(config.unit_cell),
            hklrange=(
                h_values.astype(np.float64),
                k_values.astype(np.float64),
                l_probe,
            ),
        )
    )
    if allowed.size == 0:
        return []
    allowed = allowed[np.lexsort((allowed[:, 1], allowed[:, 0]))]

    # The rod box lives in the output frame, so the reachability test has to
    # be made there too. The coverage cloud maps exactly, because hkl -> frame
    # is the single fixed matrix above.
    coverage_frame = coverage @ matrix.T
    half_width = np.broadcast_to(
        np.asarray(half_width, dtype=np.float64).reshape(-1), (3,)
    )

    grids = []
    for h, k in allowed:
        endpoints = np.asarray(
            [[float(h), float(k), l_lower], [float(h), float(k), l_upper]]
        )
        minimum, maximum = _feature_box(matrix, endpoints, half_width)
        if not np.any(
            np.all(
                (coverage_frame >= minimum) & (coverage_frame <= maximum),
                axis=1,
            )
        ):
            continue
        grids.append(
            _feature_grid(
                minimum,
                maximum,
                step,
                frame,
                f"ctr_{_format_index(h)}_{_format_index(k)}",
                chunk_shape,
            )
        )
    return _check_grid_count(grids, max_grids, "CTR")


def derive_fractional_rod_grids(
    config,
    scan,
    *,
    step,
    half_width=(0.05, 0.05, 0.0),
    frame="hkl",
    h_limits=None,
    k_limits=None,
    families="any half",
    exclude_rods="",
    coverage=None,
    detector_samples=33,
    frame_samples=128,
    chunk_shape=(64, 64, 64),
    max_grids=DEFAULT_MAX_GRIDS,
):
    """Select fractional ``(H, K)`` rods across the measured ``L`` range.

    Every candidate has at least one noninteger in-plane index. Fractional
    rods are explicit existence hypotheses; the bulk unit-cell structure
    factor is not used. The rod and coverage boxes use the same fixed frame
    convention as :func:`derive_ctr_grids`.

    :param config:
        Experiment snapshot with the reference UB matrix.
    :param scan:
        Active scan backend, used when ``coverage`` is not supplied.
    :param step:
        Three voxel widths in r.l.u. for ``hkl`` or ``Angstrom^-1`` for
        ``crystal``.
    :param half_width:
        Per-axis box padding in the output frame's units. The third entry
        extends the measured ``L`` range, and defaults to zero.
    :param str frame:
        ``hkl`` or ``crystal``.
    :param h_limits:
        Integer bounds for candidate ``H`` in r.l.u.; a scalar is symmetric.
        ``None`` derives bounds from measured coverage.
    :param k_limits:
        As ``h_limits``, for ``K``.
    :param str families:
        Semicolon-separated fractional-part patterns such as ``1/2,0`` or
        ``1/2,1/2 where h+k even``. ``any half`` and ``all half`` are
        shortcuts. Each rule selects a family; families are combined by OR.
    :param str exclude_rods:
        Semicolon-separated phase patterns or conditions such as
        ``1/2,*`` and ``h+k odd``. Exclusions are combined by OR and run
        after family selection.
    :param coverage:
        Optional sampled ``(n, 3)`` HKL coverage in reference r.l.u.
    :param int detector_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param int frame_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param chunk_shape:
        HDF5 chunk shape shared by returned grids.
    :param max_grids:
        Refuse more grids than this count; ``None`` disables the guard.
    :returns:
        Rod grids ordered by ``(H, K)``, named ``fractional_rod_<H>_<K>``.
        ``p`` denotes a fraction bar and ``m`` a minus sign.
    :rtype: list[orgui.reconstruction_job.ReconstructionGrid]
    :raises ValueError:
        If rules, limits, frame, steps, or grid count are invalid.
    """
    matrix = hkl_to_frame_matrix(config.ub_calculator, frame)
    coverage = _prepare_coverage(
        config, scan, coverage, detector_samples, frame_samples
    )
    coverage_minimum = np.min(coverage, axis=0)
    coverage_maximum = np.max(coverage, axis=0)
    h_values = _axis_integers(
        h_limits, coverage_minimum[0], coverage_maximum[0]
    )
    k_values = _axis_integers(
        k_limits, coverage_minimum[1], coverage_maximum[1]
    )
    indices = _fractional_indices(
        families,
        exclude_rods,
        ((int(h_values[0]), int(h_values[-1])),
         (int(k_values[0]), int(k_values[-1]))),
        "hk",
    )
    coverage_frame = coverage @ matrix.T
    half_width = np.broadcast_to(
        np.asarray(half_width, dtype=np.float64).reshape(-1), (3,)
    )
    l_lower = float(coverage_minimum[2])
    l_upper = float(coverage_maximum[2])
    grids = []
    for h, k in indices:
        endpoints = np.asarray(
            [[h, k, l_lower], [h, k, l_upper]], dtype=np.float64
        )
        minimum, maximum = _feature_box(matrix, endpoints, half_width)
        if not np.any(np.all(
            (coverage_frame >= minimum) & (coverage_frame <= maximum), axis=1
        )):
            continue
        name = (
            f"fractional_rod_{_format_fractional_index(h)}_"
            f"{_format_fractional_index(k)}"
        )
        grids.append(
            _feature_grid(minimum, maximum, step, frame, name, chunk_shape)
        )
    return _check_grid_count(grids, max_grids, "fractional rod")


def derive_bragg_grids(
    config,
    scan,
    *,
    step,
    half_width=(0.1, 0.1, 0.1),
    frame="hkl",
    h_limits=None,
    k_limits=None,
    l_limits=None,
    strain=None,
    coverage=None,
    detector_samples=33,
    frame_samples=128,
    chunk_shape=(64, 64, 64),
    max_grids=DEFAULT_MAX_GRIDS,
):
    """One box-shaped output grid per Bragg reflection the scan reaches.

    Reflections are the allowed integer ``(H, K, L)`` of the reference unit
    cell within the enumerated index limits. Each becomes a box of
    ``center +- half_width`` in the output frame, and is kept only if a sampled
    coverage point falls inside that box.

    With ``strain``, the reflections belong to a strained lattice whose
    lattice constants are ``a_i * (1 + strain_i)`` at unchanged cell angles and
    orientation -- the same convention as the rocking-scan Bragg extraction
    (``orGUI.get_Bragg_rocking_coordinates``). Each reciprocal basis vector
    then scales by ``1 / (1 + strain_i)``, so the strained reflection
    ``(H, K, L)`` is centered at reference ``(H / (1 + strain_a),
    K / (1 + strain_b), L / (1 + strain_c))`` r.l.u. Grid frame and names stay
    those of the reference lattice and the strained indices respectively.

    :param config:
        Central :class:`~orgui.app.config_data.ConfigData` experiment snapshot.
    :param scan:
        Active scan backend providing exposure angle bounds in radians.
    :param step:
        Three voxel widths, in r.l.u. for ``hkl`` and ``Angstrom^-1`` for
        ``crystal``.
    :param half_width:
        Half box size on each output-frame axis, in the same units as
        ``step``. A single number is read as the symmetric cube.
    :param str frame:
        ``hkl`` or ``crystal``; see :data:`STATIC_FRAMES`.
    :param h_limits:
        Optional ``(lower, upper)`` integer ``H`` limits in r.l.u.; a single
        number ``n`` is read as the symmetric ``(-n, n)``, and ``None`` takes
        the limits from the measured coverage.
    :param k_limits:
        As ``h_limits``, for ``K``.
    :param l_limits:
        As ``h_limits``, for ``L``. All index limits refer to the strained
        lattice when ``strain`` is given.
    :param strain:
        Optional fractional strain ``(strain_a, strain_b, strain_c)`` of the
        direct lattice constants, dimensionless (``0.01`` is 1 %). A single
        number applies to all three. Each entry must be above ``-1``.
        ``None`` or zeros select the reference lattice.
    :param coverage:
        Optional pre-computed ``(n, 3)`` r.l.u. coverage cloud from
        :func:`sample_hkl_coverage`, in the reference lattice.
    :param int detector_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param int frame_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param chunk_shape:
        HDF5 chunk shape in voxels, shared by every returned grid.
    :param max_grids:
        Refuse to return more than this many grids; ``None`` disables the
        guard.
    :returns:
        Reflection grids ordered by ``(H, K, L)``, named
        ``bragg_<H>_<K>_<L>`` with ``m`` standing in for a minus sign.
    :rtype: list[orgui.reconstruction_job.ReconstructionGrid]
    :raises ValueError:
        If the frame rotates with the sample, the steps are not positive, a
        strain is not above ``-1``, or the selection exceeds ``max_grids``.
    """
    lattice_scale = np.ones(3)
    if strain is not None:
        strain = np.broadcast_to(
            np.asarray(strain, dtype=np.float64).reshape(-1), (3,)
        )
        if np.any(~np.isfinite(strain)) or np.any(strain <= -1.0):
            raise ValueError(
                "Lattice strain must be finite and above -1 (-100 %)"
            )
        # Direct lattice constant a_i -> a_i (1 + strain_i); reference
        # r.l.u. = strained index / (1 + strain_i).
        lattice_scale = 1.0 + strain
    matrix = hkl_to_frame_matrix(config.ub_calculator, frame)
    coverage = _prepare_coverage(
        config, scan, coverage, detector_samples, frame_samples
    )
    # Enumerate in strained indices: the positive per-axis scale keeps the
    # coverage bounds ordered.
    coverage_minimum = np.min(coverage, axis=0) * lattice_scale
    coverage_maximum = np.max(coverage, axis=0) * lattice_scale

    axis_values = [
        _axis_integers(limits, coverage_minimum[axis], coverage_maximum[axis])
        for axis, limits in enumerate((h_limits, k_limits, l_limits))
    ]

    allowed = np.atleast_2d(
        ReciprocalNavigation.allowedReflections(
            _bulk_unit_cell(config.unit_cell),
            hklrange=tuple(
                values.astype(np.float64) for values in axis_values
            ),
        )
    )
    if allowed.size == 0:
        return []
    allowed = allowed[np.lexsort((allowed[:, 2], allowed[:, 1], allowed[:, 0]))]

    coverage_frame = coverage @ matrix.T
    half_width = np.broadcast_to(
        np.asarray(half_width, dtype=np.float64).reshape(-1), (3,)
    )

    grids = []
    for hkl in allowed:
        center = (
            np.asarray([[float(value) for value in hkl]]) / lattice_scale
        )
        minimum, maximum = _feature_box(matrix, center, half_width)
        if not np.any(
            np.all(
                (coverage_frame >= minimum) & (coverage_frame <= maximum),
                axis=1,
            )
        ):
            continue
        name = "bragg_" + "_".join(_format_index(value) for value in hkl)
        grids.append(
            _feature_grid(minimum, maximum, step, frame, name, chunk_shape)
        )
    return _check_grid_count(grids, max_grids, "Bragg")


def derive_fractional_bragg_grids(
    config,
    scan,
    *,
    step,
    half_width=(0.1, 0.1, 0.1),
    frame="hkl",
    h_limits=None,
    k_limits=None,
    l_limits=None,
    strain=None,
    families="any half",
    exclude_peaks="",
    coverage=None,
    detector_samples=33,
    frame_samples=128,
    chunk_shape=(64, 64, 64),
    max_grids=DEFAULT_MAX_GRIDS,
):
    """Select compact boxes around explicit fractional ``(H, K, L)`` peaks.

    Every candidate has at least one noninteger index. The reference bulk
    structure factor is not consulted. With strain, indices belong to the
    strained lattice and centers are mapped into the reference lattice as
    ``(H, K, L) / (1 + strain)`` before converting to the output frame, the
    same convention as :func:`derive_bragg_grids`.

    :param config:
        Experiment snapshot with the reference UB matrix.
    :param scan:
        Active scan backend, used when ``coverage`` is not supplied.
    :param step:
        Three voxel widths in r.l.u. for ``hkl`` or ``Angstrom^-1`` for
        ``crystal``.
    :param half_width:
        Three per-axis box half-widths in the output frame's units.
    :param str frame:
        ``hkl`` or ``crystal``.
    :param h_limits:
        Integer bounds on ``H`` in the enumerated lattice's r.l.u.; a scalar
        is symmetric. ``None`` derives bounds from measured coverage.
    :param k_limits:
        As ``h_limits``, for ``K``.
    :param l_limits:
        As ``h_limits``, for ``L``.
    :param strain:
        Optional fractional direct-lattice strain ``(a, b, c)``; each entry
        is dimensionless and greater than ``-1``. A scalar applies to all.
    :param str families:
        Semicolon-separated fractional-part patterns such as ``1/2,0,0``
        or ``1/2,1/2,0 where h+k even``. ``any half`` and ``all half`` are
        shortcuts. Families are combined by OR.
    :param str exclude_peaks:
        Semicolon-separated phase patterns or conditions such as
        ``1/2,*,*`` and ``h+k odd``, combined by OR.
    :param coverage:
        Optional sampled ``(n, 3)`` HKL coverage in reference r.l.u.
    :param int detector_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param int frame_samples:
        Forwarded to :func:`sample_hkl_coverage`.
    :param chunk_shape:
        HDF5 chunk shape shared by returned grids.
    :param max_grids:
        Refuse more grids than this count; ``None`` disables the guard.
    :returns:
        Grids ordered by ``(H, K, L)`` and named
        ``fractional_bragg_<H>_<K>_<L>``. ``p`` denotes a fraction bar and
        ``m`` a minus sign.
    :rtype: list[orgui.reconstruction_job.ReconstructionGrid]
    :raises ValueError:
        If rules, limits, strain, frame, steps, or grid count are invalid.
    """
    lattice_scale = np.ones(3)
    if strain is not None:
        strain = np.broadcast_to(
            np.asarray(strain, dtype=np.float64).reshape(-1), (3,)
        )
        if np.any(~np.isfinite(strain)) or np.any(strain <= -1.0):
            raise ValueError("Lattice strain must be finite and above -1 (-100 %)")
        lattice_scale = 1.0 + strain
    matrix = hkl_to_frame_matrix(config.ub_calculator, frame)
    coverage = _prepare_coverage(
        config, scan, coverage, detector_samples, frame_samples
    )
    coverage_minimum = np.min(coverage, axis=0) * lattice_scale
    coverage_maximum = np.max(coverage, axis=0) * lattice_scale
    axis_values = [
        _axis_integers(limits, coverage_minimum[axis], coverage_maximum[axis])
        for axis, limits in enumerate((h_limits, k_limits, l_limits))
    ]
    indices = _fractional_indices(
        families,
        exclude_peaks,
        tuple((int(values[0]), int(values[-1])) for values in axis_values),
        "hkl",
    )
    coverage_frame = coverage @ matrix.T
    half_width = np.broadcast_to(
        np.asarray(half_width, dtype=np.float64).reshape(-1), (3,)
    )
    grids = []
    for hkl in indices:
        center = np.asarray([hkl], dtype=np.float64) / lattice_scale
        minimum, maximum = _feature_box(matrix, center, half_width)
        if not np.any(np.all(
            (coverage_frame >= minimum) & (coverage_frame <= maximum), axis=1
        )):
            continue
        name = "fractional_bragg_" + "_".join(
            _format_fractional_index(value) for value in hkl
        )
        grids.append(
            _feature_grid(minimum, maximum, step, frame, name, chunk_shape)
        )
    return _check_grid_count(grids, max_grids, "fractional Bragg")

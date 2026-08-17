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
- a *Bragg reflection* is a compact box around one allowed integer
  ``(H, K, L)``.

Both selections enumerate features from the reference unit cell, keep only the
ones the active scan actually reaches, and return one
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

import numpy as np

from .datautils.xrayutils import ReciprocalNavigation
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
        raise ValueError("Exposure angle bounds must have shape (frames, 2, 4)")
    excluded = set(config.corrections.excluded_frames)
    included = np.asarray(
        [index for index in range(len(scan)) if index not in excluded],
        dtype=np.int64,
    )
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
                    coordinates.append(
                        kernel.coordinate(
                            rays, angles_start, angles_end, 0, 0, 0.0, 0.0, t
                        )
                    )
    return np.asarray(coordinates, dtype=np.float64)


def _bulk_unit_cell(xtal):
    """Return the bulk unit cell backing a crystal or unit-cell object."""
    return getattr(xtal, "uc_bulk", xtal)


def _format_index(value):
    return str(int(round(float(value)))).replace("-", "m")


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
        As ``h_limits``, for ``L``.
    :param coverage:
        Optional pre-computed ``(n, 3)`` r.l.u. coverage cloud from
        :func:`sample_hkl_coverage`.
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
        If the frame rotates with the sample, the steps are not positive, or
        the selection exceeds ``max_grids``.
    """
    matrix = hkl_to_frame_matrix(config.ub_calculator, frame)
    coverage = _prepare_coverage(
        config, scan, coverage, detector_samples, frame_samples
    )
    coverage_minimum = np.min(coverage, axis=0)
    coverage_maximum = np.max(coverage, axis=0)

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
        center = np.asarray([[float(value) for value in hkl]])
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

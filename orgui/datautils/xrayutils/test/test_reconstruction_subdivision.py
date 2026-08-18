"""Oracle and invariant tests for the adaptive pixel-subdivision path.

These pin the *behaviour* of ``split_pixel`` and its two fast paths (the
``max_depth == 0`` centre-only case and ``split_pixel_stationary_depth2``)
against an independent Python implementation of the same adaptive rule, so
that the native side can be restructured -- shared-corner cache replaced by
corner-passing recursion -- with something better than self-consistency to
check it against.

The reference deliberately obtains every voxel through the kernel's own
``coordinate()``, so it cannot drift from the kernel on frame conventions,
ray interpolation or rounding: what it independently reproduces is the
subdivision *algorithm* -- which lattice points are visited, which cells
resolve, and what weight each leaf carries.
"""

from __future__ import annotations

import numpy as np
import pytest


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


DETECTOR_DISTANCE_MM = 500.0
PIXEL_PITCH_MM = 0.172
WAVEVECTOR = 7.601
# A triclinic-ish UB, so the hkl frame applies a genuinely non-diagonal
# transform rather than a scaled identity.
UB = np.array(
    [
        [1.6122, -0.5878, 0.0000],
        [0.0000, 1.4899, -0.1051],
        [0.0000, 0.0000, 1.2083],
    ]
)


def _detector_rays(rows, columns, row_origin=0, column_origin=0):
    """Unit rays at pixel corners of a flat detector facing the beam.

    :param int rows, columns:
        Tile size in pixels; the returned array is one larger in each axis.
    :param int row_origin, column_origin:
        Offset of the tile within a detector centred on the direct beam, so
        a tile can be placed off-axis where the angular spread is largest.
    :returns:
        ``(rows + 1, columns + 1, 3)`` unit rays, beam along ``+y``.
    """
    row_edges = np.arange(rows + 1) + row_origin - 0.5
    column_edges = np.arange(columns + 1) + column_origin - 0.5
    x = column_edges[None, :] * PIXEL_PITCH_MM
    z = row_edges[:, None] * PIXEL_PITCH_MM
    x, z = np.broadcast_arrays(np.broadcast_to(x, (rows + 1, columns + 1)), z)
    rays = np.stack(
        (x, np.full_like(x, DETECTOR_DISTANCE_MM), z), axis=-1
    ).astype(np.float64)
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)
    return np.ascontiguousarray(rays)


def _probe_kernel():
    """A kernel used only through ``coordinate()``, to bound the volume."""
    return native.ReconstructionKernel(
        np.full(3, -1e12),
        np.ones(3),
        np.full(3, 2_000_000, dtype=np.int64),
        np.ones(3, dtype=np.int64),
        "hkl",
        WAVEVECTOR,
        np.ascontiguousarray(np.linalg.inv(UB)),
        np.ascontiguousarray(np.eye(3)),
        0,
        1,
        1,
        1024 * 1024,
    )


def _grid_for(rays, angles_start, angles_end, voxels_per_pixel, chunk_shape):
    """A grid whose voxel is a chosen multiple of the pixel footprint.

    Sizing the voxel against the real pixel-to-pixel step is what makes the
    subdivision actually subdivide: a grid much coarser than the detector
    resolves every cell at depth zero and exercises nothing.

    :param float voxels_per_pixel:
        Target sample pitch measured in voxels, applied to whichever of the
        row, column and exposure directions moves the coordinate furthest
        along each axis. Below one, neighbouring samples usually share a
        voxel; above one, most cells straddle a boundary and never resolve.
    :returns:
        ``(minimum, step, shape)``.
    """
    probe = _probe_kernel()
    rows = rays.shape[0] - 1
    columns = rays.shape[1] - 1
    corners = []
    for row in (0, rows - 1):
        for column in (0, columns - 1):
            for u, v, t in ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)):
                corners.append(
                    probe.coordinate(
                        rays, angles_start, angles_end, row, column, u, v, t
                    )
                )
    corners = np.asarray(corners)

    def at(row, column, t):
        return np.asarray(
            probe.coordinate(
                rays, angles_start, angles_end, row, column, 0.5, 0.5, t
            )
        )

    # Each axis is sized against the largest step any of the three sampling
    # directions takes along it. Sizing from the column direction alone
    # leaves the axis a column barely moves along thousands of voxels wide
    # per *row*, so no cell ever resolves and subdivision never converges.
    centre = at(0, 0, 0.5)
    pitch = np.maximum.reduce(
        [
            np.abs(at(0, 1, 0.5) - centre),
            np.abs(at(1, 0, 0.5) - centre),
            np.abs(at(0, 0, 1.0) - at(0, 0, 0.0)),
        ]
    )
    step = np.maximum(pitch, np.max(pitch) * 0.2) / voxels_per_pixel
    lower = corners.min(axis=0) - 4.0 * step
    upper = corners.max(axis=0) + 4.0 * step
    shape = np.ceil((upper - lower) / step).astype(np.int64)
    # An offset that is not a whole number of voxels, so no coordinate is
    # likely to sit within rounding distance of a voxel edge.
    lower = lower - 0.31718 * step
    return lower, step, shape


def _kernel(minimum, step, shape, chunk_shape, max_depth, threads, block):
    return native.ReconstructionKernel(
        np.ascontiguousarray(minimum),
        np.ascontiguousarray(step),
        np.ascontiguousarray(shape),
        np.ascontiguousarray(np.asarray(chunk_shape, dtype=np.int64)),
        "hkl",
        WAVEVECTOR,
        np.ascontiguousarray(np.linalg.inv(UB)),
        np.ascontiguousarray(np.eye(3)),
        max_depth,
        threads,
        block,
        1 << 40,
    )


def _reference_leaves(
    voxel_at, max_depth, stationary
):
    """Leaves one pixel's adaptive subdivision produces.

    An independent implementation of ``split_pixel``: a cell resolves when
    all of its corners agree on one valid voxel, otherwise it subdivides,
    and a cell still unresolved at ``max_depth`` contributes its whole
    weight at its centre. Child visit order matches the native explicit
    stack (children pushed high-to-low, popped from the back), so leaves
    come out in the same order and floating-point summation matches.

    :param voxel_at:
        ``(u, v, t) -> (voxel, valid)``.
    :returns:
        List of ``(voxel, weight)`` in emission order.
    """
    if max_depth == 0:
        # Mirrors the kernel's centre-only fast path, which skips the
        # corner test entirely and assigns weight exactly one.
        voxel, valid = voxel_at(0.5, 0.5, 0.5)
        return [(voxel, 1.0)] if valid else []

    corner_count = 4 if stationary else 8
    leaves = []
    stack = [
        (
            0.0,
            1.0,
            0.0,
            1.0,
            0.5 if stationary else 0.0,
            0.5 if stationary else 1.0,
            1.0,
            0,
        )
    ]
    while stack:
        u0, u1, v0, v1, t0, t1, weight, depth = stack.pop()
        first_voxel = None
        first_valid = False
        all_same = True
        for corner in range(corner_count):
            u = u1 if corner & 1 else u0
            v = v1 if corner & 2 else v0
            t = 0.5 if stationary else (t1 if corner & 4 else t0)
            voxel, valid = voxel_at(u, v, t)
            if corner == 0:
                first_valid = valid
                first_voxel = voxel
            elif valid != first_valid or (valid and voxel != first_voxel):
                all_same = False
        if all_same and first_valid:
            leaves.append((first_voxel, weight))
            continue
        if depth >= max_depth:
            voxel, valid = voxel_at(
                0.5 * (u0 + u1), 0.5 * (v0 + v1), 0.5 * (t0 + t1)
            )
            if valid:
                leaves.append((voxel, weight))
            continue
        um = 0.5 * (u0 + u1)
        vm = 0.5 * (v0 + v1)
        tm = 0.5 * (t0 + t1)
        child_weight = weight / corner_count
        for child in range(corner_count - 1, -1, -1):
            stack.append(
                (
                    um if child & 1 else u0,
                    u1 if child & 1 else um,
                    vm if child & 2 else v0,
                    v1 if child & 2 else vm,
                    0.5 if stationary else (tm if child & 4 else t0),
                    0.5 if stationary else (t1 if child & 4 else tm),
                    child_weight,
                    depth + 1,
                )
            )
    return leaves


def _reference_records(
    kernel,
    rays,
    intensity,
    variance,
    mask,
    angles_start,
    angles_end,
    minimum,
    step,
    shape,
    chunk_shape,
    max_depth,
):
    """Records the kernel should emit, built from the reference subdivision.

    :returns:
        ``{(chunk_id, local_voxel_id): [weighted_intensity,
        weighted_variance, weight, contributors]}``.
    """
    stationary = bool(np.array_equal(angles_start, angles_end))
    chunk_shape = np.asarray(chunk_shape, dtype=np.int64)
    chunk_grid = np.ceil(np.asarray(shape) / chunk_shape).astype(np.int64)
    rows, columns = intensity.shape
    records: dict[tuple[int, int], list] = {}

    for row in range(rows):
        for column in range(columns):
            if (
                mask[row, column]
                or not np.isfinite(intensity[row, column])
                or not np.isfinite(variance[row, column])
                or variance[row, column] < 0.0
            ):
                continue
            cache: dict[tuple[float, float, float], tuple] = {}

            def voxel_at(u, v, t, row=row, column=column, cache=cache):
                key = (u, v, t)
                hit = cache.get(key)
                if hit is not None:
                    return hit
                coordinate = np.asarray(
                    kernel.coordinate(
                        rays, angles_start, angles_end, row, column, u, v, t
                    )
                )
                index = np.floor((coordinate - minimum) / step)
                valid = bool(
                    np.all(index >= 0.0) and np.all(index < np.asarray(shape))
                )
                voxel = tuple(index.astype(np.int64)) if valid else None
                cache[key] = (voxel, valid)
                return voxel, valid

            leaves = _reference_leaves(voxel_at, max_depth, stationary)
            if not leaves:
                continue
            # The kernel sorts a pixel's leaves by packed voxel id and sums
            # the weight of each run before touching the accumulator, so a
            # pixel contributes at most once per distinct voxel.
            leaves.sort(key=lambda entry: entry[0])
            offset = 0
            while offset < len(leaves):
                voxel = leaves[offset][0]
                weight = 0.0
                while offset < len(leaves) and leaves[offset][0] == voxel:
                    weight += leaves[offset][1]
                    offset += 1
                index = np.asarray(voxel, dtype=np.int64)
                chunk_index, local_index = np.divmod(index, chunk_shape)
                key = (
                    int(
                        (chunk_index[0] * chunk_grid[1] + chunk_index[1])
                        * chunk_grid[2]
                        + chunk_index[2]
                    ),
                    int(
                        (local_index[0] * chunk_shape[1] + local_index[1])
                        * chunk_shape[2]
                        + local_index[2]
                    ),
                )
                entry = records.setdefault(key, [0.0, 0.0, 0.0, 0])
                entry[0] += weight * intensity[row, column]
                entry[1] += weight * weight * variance[row, column]
                entry[2] += weight
                entry[3] += 1
    return records


def _payload(rows, columns, seed=0):
    rng = np.random.default_rng(seed)
    intensity = rng.uniform(1.0, 100.0, size=(rows, columns))
    variance = np.abs(rng.normal(4.0, 0.5, size=(rows, columns)))
    mask = np.zeros((rows, columns), dtype=bool)
    mask[rows // 3, columns // 2] = True
    return intensity, variance, mask


def _case(rows, columns, *, moving, voxels_per_pixel=0.5, chunk_shape=(16, 16, 16)):
    """Assemble one complete geometry/payload case."""
    rays = _detector_rays(rows, columns, row_origin=300, column_origin=400)
    angles_start = np.array([0.10471975511965978, 0.7, 0.0, 0.0])
    angles_end = (
        angles_start + np.array([0.0, 0.00174533, 0.0, 0.0])
        if moving
        else angles_start.copy()
    )
    minimum, step, shape = _grid_for(
        rays, angles_start, angles_end, voxels_per_pixel, chunk_shape
    )
    intensity, variance, mask = _payload(rows, columns)
    return {
        "rays": rays,
        "angles_start": np.ascontiguousarray(angles_start),
        "angles_end": np.ascontiguousarray(angles_end),
        "minimum": minimum,
        "step": step,
        "shape": shape,
        "chunk_shape": chunk_shape,
        "intensity": intensity,
        "variance": variance,
        "mask": mask,
    }


def _pixel_footprint(case, row, column):
    """The mapped coordinates of one pixel's four corners and its centre.

    :returns:
        ``(corners, centre)``; ``corners`` is ``(4, 3)``, ``centre`` is
        ``(3,)``, both in the grid's own r.l.u. frame.
    """
    probe = _probe_kernel()

    def at(u, v, t):
        return np.asarray(
            probe.coordinate(
                case["rays"],
                case["angles_start"],
                case["angles_end"],
                row,
                column,
                u,
                v,
                t,
            )
        )

    corners = np.asarray(
        [at(u, v, 0.0) for u in (0.0, 1.0) for v in (0.0, 1.0)]
    )
    return corners, at(0.5, 0.5, 0.5)


def _with_grid(case, minimum, step, shape):
    """``case`` with its output grid replaced, payload and rays untouched."""
    replaced = dict(case)
    replaced["minimum"] = np.ascontiguousarray(np.asarray(minimum, dtype=np.float64))
    replaced["step"] = np.ascontiguousarray(np.asarray(step, dtype=np.float64))
    replaced["shape"] = np.ascontiguousarray(np.asarray(shape, dtype=np.int64))
    return replaced


def _assert_matches_reference(case, max_depth):
    """The native records for ``case`` are exactly the reference's."""
    result = _accumulate(case, max_depth)
    kernel = _kernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
        1,
        64,
    )
    expected = _reference_records(
        kernel,
        case["rays"],
        case["intensity"],
        case["variance"],
        case["mask"],
        case["angles_start"],
        case["angles_end"],
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
    )
    produced = dict(
        zip(
            zip(
                result["chunk_id"].tolist(),
                result["local_voxel_id"].tolist(),
            ),
            zip(
                result["weighted_intensity"].tolist(),
                result["weighted_variance"].tolist(),
                result["weight"].tolist(),
                result["contributors"].tolist(),
            ),
        )
    )
    assert set(produced) == set(expected)
    for key, (intensity, variance, weight, contributors) in produced.items():
        reference = expected[key]
        assert contributors == reference[3], key
        assert intensity == pytest.approx(reference[0], rel=1e-12), key
        assert variance == pytest.approx(reference[1], rel=1e-12), key
        assert weight == pytest.approx(reference[2], rel=1e-12), key
    return result


def _accumulate(case, max_depth, threads=1, block=64, profile=False):
    kernel = _kernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
        threads,
        block,
    )
    return kernel.accumulate(
        case["intensity"],
        case["variance"],
        case["mask"],
        case["rays"],
        case["angles_start"],
        case["angles_end"],
        profile,
    )


@pytest.mark.parametrize("max_depth", [0, 1, 2, 3, 4])
def test_stationary_subdivision_matches_reference(max_depth):
    """Every depth of the stationary path reproduces the reference exactly.

    Depth 0 covers the centre-only fast path and depth 2 the hand-written
    ``split_pixel_stationary_depth2``; the rest exercise generic
    ``split_pixel``, cached below depth 4 and uncached above it.
    """
    case = _case(6, 7, moving=False)
    result = _accumulate(case, max_depth)
    kernel = _kernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
        1,
        64,
    )
    expected = _reference_records(
        kernel,
        case["rays"],
        case["intensity"],
        case["variance"],
        case["mask"],
        case["angles_start"],
        case["angles_end"],
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
    )

    produced = dict(
        zip(
            zip(
                result["chunk_id"].tolist(),
                result["local_voxel_id"].tolist(),
            ),
            zip(
                result["weighted_intensity"].tolist(),
                result["weighted_variance"].tolist(),
                result["weight"].tolist(),
                result["contributors"].tolist(),
            ),
        )
    )
    assert len(expected) > 1, "geometry must reach more than one voxel"
    assert set(produced) == set(expected)
    for key, (intensity, variance, weight, contributors) in produced.items():
        reference = expected[key]
        assert contributors == reference[3], key
        assert intensity == pytest.approx(reference[0], rel=1e-12), key
        assert variance == pytest.approx(reference[1], rel=1e-12), key
        assert weight == pytest.approx(reference[2], rel=1e-12), key


@pytest.mark.parametrize("max_depth", [0, 1, 2, 3, 4])
def test_moving_subdivision_matches_reference(max_depth):
    """The eight-corner octree path reproduces the reference exactly.

    Continuous exposures are what ``angle_fallback="midpoint"`` produces and
    what a rotating scan really does. Depth 4 matters most here: it is where
    the shared-corner cache used to switch off, so it is the case the
    corner-passing recursion changed the most.
    """
    case = _case(5, 5, moving=True)
    result = _accumulate(case, max_depth)
    kernel = _kernel(
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
        1,
        64,
    )
    expected = _reference_records(
        kernel,
        case["rays"],
        case["intensity"],
        case["variance"],
        case["mask"],
        case["angles_start"],
        case["angles_end"],
        case["minimum"],
        case["step"],
        case["shape"],
        case["chunk_shape"],
        max_depth,
    )
    produced = set(
        zip(result["chunk_id"].tolist(), result["local_voxel_id"].tolist())
    )
    assert len(expected) > 1
    assert produced == set(expected)
    weights = dict(
        zip(
            zip(
                result["chunk_id"].tolist(),
                result["local_voxel_id"].tolist(),
            ),
            result["weight"].tolist(),
        )
    )
    for key, weight in weights.items():
        assert weight == pytest.approx(expected[key][2], rel=1e-12), key


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [0, 1, 2, 3, 4])
def test_subdivision_conserves_weight_and_intensity(max_depth, moving):
    """Subdivision partitions each pixel: weights sum to the pixel count.

    A pixel's leaves carry weights summing to exactly one unless part of it
    leaves the grid, so with a grid that contains the whole tile the total
    weight is the valid-pixel count and the total weighted intensity is the
    plain sum of those pixels' intensities -- at every depth. This is the
    invariant that a corner-passing rewrite could most plausibly break, by
    dropping or double-counting a child.
    """
    case = _case(6, 6, moving=moving, voxels_per_pixel=0.5)
    result = _accumulate(case, max_depth)
    valid = ~case["mask"]
    assert float(np.sum(result["weight"])) == pytest.approx(
        float(np.count_nonzero(valid)), rel=1e-12
    )
    assert float(np.sum(result["weighted_intensity"])) == pytest.approx(
        float(np.sum(case["intensity"][valid])), rel=1e-12
    )


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [0, 1, 2, 3, 4])
def test_subdivision_is_independent_of_threads_and_block_size(
    max_depth, moving
):
    """Splitting the tile differently across workers changes nothing.

    Block size decides how pixels group into per-block accumulators and
    thread count decides who runs them; neither may reach the scientific
    result. Checked on the full output arrays, not summary statistics, so
    an ordering or merge regression cannot hide inside a matching total.
    """
    case = _case(6, 7, moving=moving)
    reference = _accumulate(case, max_depth, threads=1, block=4096)
    for threads, block in ((1, 7), (4, 5), (8, 4096)):
        candidate = _accumulate(case, max_depth, threads=threads, block=block)
        for name in ("chunk_id", "local_voxel_id", "contributors"):
            np.testing.assert_array_equal(
                candidate[name], reference[name], err_msg=f"{name} {threads}/{block}"
            )
        for name in ("weighted_intensity", "weighted_variance", "weight"):
            np.testing.assert_allclose(
                candidate[name],
                reference[name],
                rtol=1e-12,
                err_msg=f"{name} {threads}/{block}",
            )


def test_output_keys_are_sorted_and_unique_at_every_depth():
    """The merge contract every downstream stage relies on."""
    for moving in (False, True):
        case = _case(6, 7, moving=moving)
        for max_depth in range(0, 5):
            result = _accumulate(case, max_depth, threads=4, block=8)
            key = (
                result["chunk_id"].astype(np.uint64) << np.uint64(32)
            ) | result["local_voxel_id"].astype(np.uint64)
            assert key.size > 1
            assert np.all(np.diff(key) > 0), (moving, max_depth)


@pytest.mark.parametrize("moving", [False, True])
def test_deeper_subdivision_resolves_more_voxels_but_converges(moving):
    """Depth buys voxels quickly and then only refines weights.

    Guards the premise the subdivision exists for. Note that the count is
    **not** monotone in depth and must not be asserted to be: a cell still
    unresolved at ``max_depth`` contributes at its own centre, so raising
    the depth replaces one centre sample by four elsewhere, and a voxel the
    coarser sample happened to land in can go unvisited. Real frames average
    that away over millions of pixels; a small tile does not.
    """
    case = _case(16, 16, moving=moving)
    counts = [
        _accumulate(case, max_depth)["chunk_id"].size for max_depth in range(0, 6)
    ]
    assert counts[1] > counts[0], counts
    assert counts[-1] < 2 * counts[1], counts
    # Converged: the last step changes the count by only a few per cent.
    assert counts[-1] - counts[-2] <= 0.05 * counts[-2], counts


# ---------------------------------------------------------------------------
# Whole-pixel rejection (sparse output volumes)
#
# Automatic reciprocal-space volume selection maps many small grids in one
# pass, so most pixels miss most grids. The subdivision resolves a cell whose
# corners share one voxel, but had no case for a cell whose corners are all
# *outside* the grid, so a missing pixel used to cost the full lattice. These
# pin both halves of the guard that fixes it: that it fires, and that it never
# fires when the footprint could still reach the grid.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [1, 2, 3])
def test_unreachable_grid_is_settled_at_the_pixel_corners(moving, max_depth):
    """A grid no pixel can reach costs only that pixel's corners.

    Four for a stationary exposure, eight for a moving one, whose reject
    reasons about the chord between the exposure's two ends and so needs
    both. Either way the whole lattice below them -- 41 evaluations per
    pixel at stationary depth two -- is what the guard skips.
    """
    case = _case(6, 7, moving=moving)
    # Far outside anything the detector maps to, in every axis at once.
    offset = 1000.0 * np.asarray(case["step"]) * np.asarray(case["shape"])
    unreachable = _with_grid(
        case, np.asarray(case["minimum"]) + offset, case["step"], case["shape"]
    )

    result = _accumulate(unreachable, max_depth, profile=True)
    profile = result["_profile"]

    assert result["chunk_id"].size == 0
    assert profile["valid_pixels"] == 6 * 7 - 1  # one masked pixel
    assert profile["skipped_pixels"] == profile["valid_pixels"]
    corners_per_pixel = 8 if moving else 4
    assert (
        profile["coordinate_evaluations"]
        == corners_per_pixel * profile["valid_pixels"]
    )


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [1, 2, 3])
def test_sparse_grid_matches_reference(moving, max_depth):
    """A grid only part of the detector reaches still maps exactly.

    The reference implements the subdivision *without* the reject, so any
    pixel wrongly proved to miss shows up here as a missing record.
    """
    case = _case(16, 16, moving=moving)
    step = np.asarray(case["step"])
    shape = np.asarray(case["shape"])
    # One eighth of the mapped volume along each axis, offset off-centre so
    # the boundary crosses the detector rather than clipping a corner.
    small_shape = np.maximum(1, shape // 8)
    minimum = np.asarray(case["minimum"]) + 0.4 * shape * step
    sparse = _with_grid(case, minimum, step, small_shape)

    result = _assert_matches_reference(sparse, max_depth)

    profile = _accumulate(sparse, max_depth, profile=True)["_profile"]
    assert result["chunk_id"].size > 0, "the grid must still be reached"
    assert profile["skipped_pixels"] > 0, "and most pixels must still miss it"


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [1, 2, 3])
def test_grid_smaller_than_a_pixel_footprint_is_not_rejected(moving, max_depth):
    """A grid that fits *between* a pixel's corners is still subdivided.

    The case the reject has to get right, and the reason it tests a bounding
    box rather than the cheaper "all four corners are invalid": a grid
    narrower than one pixel's footprint leaves every corner outside while the
    footprint's interior still crosses it.

    Whether that pixel then emits a record is a separate question and not
    asserted here -- the subdivision only records the voxel a *leaf centre*
    falls in, and a grid this small sits between the leaf centres at most
    depths. What must hold is that the pixel is still subdivided rather than
    proved away, and that the result matches the reference either way.
    """
    case = _case(6, 7, moving=moving)
    corners, centre = _pixel_footprint(case, 3, 3)
    spread = corners.max(axis=0) - corners.min(axis=0)
    assert np.all(spread > 0.0)
    # One voxel, an eighth of the pixel footprint across, centred on the
    # pixel centre -- so no corner of that pixel is inside it, and no
    # neighbouring pixel's footprint reaches it either.
    step = spread / 8.0
    minimum = centre - 0.5 * step
    tiny = _with_grid(case, minimum, step, np.ones(3, dtype=np.int64))
    inside = [
        np.all(corner >= minimum) and np.all(corner < minimum + step)
        for corner in corners
    ]
    assert not any(inside), "the test pixel's corners must all be outside"

    _assert_matches_reference(tiny, max_depth)

    profile = _accumulate(tiny, max_depth, profile=True)["_profile"]
    # Every pixel but the one the grid sits inside is provably outside it.
    assert profile["skipped_pixels"] == profile["valid_pixels"] - 1

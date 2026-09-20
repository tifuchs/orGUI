"""Tests for mapping several adjacent frames in one native call.

``accumulate_group`` exists to collapse redundancy that spans frames: on a
rotation scan two adjacent images land as close together in reciprocal space
as two adjacent pixels do, so a work block shaped as a brick in (row,
column, frame) shares voxels the per-image block never sees. What it must
not change is the answer -- a group has to produce what merging the same
frames one at a time produces, to floating-point association.
"""

from __future__ import annotations

import numpy as np
import pytest

from orgui.datautils.xrayutils.reconstruction import _reduce_batches


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


ROWS, COLUMNS = 24, 28
WAVEVECTOR = 7.601
UB = np.array(
    [
        [1.6122, -0.5878, 0.0000],
        [0.0000, 1.4899, -0.1051],
        [0.0000, 0.0000, 1.2083],
    ]
)


def _rays(rows=ROWS, columns=COLUMNS):
    """Unit rays at pixel corners of a flat detector, placed off-axis."""
    row_edges = np.arange(rows + 1) + 300 - 0.5
    column_edges = np.arange(columns + 1) + 400 - 0.5
    x = np.broadcast_to(column_edges[None, :] * 0.172, (rows + 1, columns + 1))
    z = np.broadcast_to(row_edges[:, None] * 0.172, (rows + 1, columns + 1))
    rays = np.stack((x, np.full_like(x, 500.0), z), axis=-1).astype(np.float64)
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)
    return np.ascontiguousarray(rays)


def _angles(frames, moving):
    """Per-frame exposure bounds for a 0.1 deg/frame rotation."""
    step = np.deg2rad(0.1)
    start = np.zeros((frames, 4))
    start[:, 0] = 0.10471975511965978
    start[:, 1] = 0.7 + step * np.arange(frames)
    end = start.copy()
    if moving:
        end[:, 1] += step
    return np.ascontiguousarray(start), np.ascontiguousarray(end)


def _grid(rays, start, end, voxels_per_pixel=0.5):
    """A grid whose voxel is a chosen multiple of the sample pitch."""
    probe = native.ReconstructionKernel(
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

    def at(frame, row, column, t):
        return np.asarray(
            probe.coordinate(
                rays, start[frame], end[frame], row, column, 0.5, 0.5, t
            )
        )

    centre = at(0, 0, 0, 0.5)
    pitch = np.maximum.reduce(
        [
            np.abs(at(0, 0, 1, 0.5) - centre),
            np.abs(at(0, 1, 0, 0.5) - centre),
            np.abs(at(0, 0, 0, 1.0) - at(0, 0, 0, 0.0)),
            np.abs(at(len(start) - 1, 0, 0, 0.5) - centre) / max(1, len(start) - 1),
        ]
    )
    step = np.maximum(pitch, np.max(pitch) * 0.2) / voxels_per_pixel
    corners = np.asarray(
        [
            at(frame, row, column, t)
            for frame in (0, len(start) - 1)
            for row in (0, ROWS - 1)
            for column in (0, COLUMNS - 1)
            for t in (0.0, 1.0)
        ]
    )
    minimum = corners.min(axis=0) - 4.0 * step - 0.31718 * step
    shape = np.ceil(
        (corners.max(axis=0) + 4.0 * step - minimum) / step
    ).astype(np.int64)
    return minimum, step, shape


def _kernel(minimum, step, shape, max_depth, threads, block):
    return native.ReconstructionKernel(
        np.ascontiguousarray(minimum),
        np.ascontiguousarray(step),
        np.ascontiguousarray(shape),
        np.array([16, 16, 16], dtype=np.int64),
        "hkl",
        WAVEVECTOR,
        np.ascontiguousarray(np.linalg.inv(UB)),
        np.ascontiguousarray(np.eye(3)),
        max_depth,
        threads,
        block,
        1 << 40,
    )


def _payload(frames, seed=0):
    rng = np.random.default_rng(seed)
    intensity = rng.uniform(1.0, 100.0, size=(frames, ROWS, COLUMNS))
    variance = np.abs(rng.normal(4.0, 0.5, size=(frames, ROWS, COLUMNS)))
    mask = np.zeros((frames, ROWS, COLUMNS), dtype=bool)
    mask[:, ROWS // 3, COLUMNS // 2] = True
    mask[frames // 2, 0, 0] = True
    return intensity, variance, mask


def _case(frames, moving, max_depth, threads=1, block=256):
    rays = _rays()
    start, end = _angles(frames, moving)
    minimum, step, shape = _grid(rays, start, end)
    kernel = _kernel(minimum, step, shape, max_depth, threads, block)
    intensity, variance, mask = _payload(frames)
    return kernel, rays, start, end, intensity, variance, mask


def _per_frame_reference(kernel, rays, start, end, intensity, variance, mask):
    """What the same frames produce mapped one at a time, then merged."""
    batches = [
        kernel.accumulate(
            np.ascontiguousarray(intensity[frame]),
            np.ascontiguousarray(variance[frame]),
            np.ascontiguousarray(mask[frame]),
            rays,
            np.ascontiguousarray(start[frame]),
            np.ascontiguousarray(end[frame]),
        )
        for frame in range(intensity.shape[0])
    ]
    return _reduce_batches(batches)


def _as_mapping(batch):
    return dict(
        zip(
            zip(
                batch["chunk_id"].tolist(),
                batch["local_voxel_id"].tolist(),
            ),
            zip(
                batch["weighted_intensity"].tolist(),
                batch["weighted_variance"].tolist(),
                batch["weight"].tolist(),
                batch["contributors"].tolist(),
            ),
        )
    )


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [0, 1, 2])
@pytest.mark.parametrize("frames", [1, 2, 8])
def test_group_matches_frames_mapped_one_at_a_time(frames, max_depth, moving):
    """The scientific contract: grouping changes cost, never the answer.

    Contributions from several frames now meet in the block map rather than
    in the checkpoint accumulator, so the sums associate differently and are
    compared to tolerance -- but every voxel, every contributor count and
    every total must agree.
    """
    case = _case(frames, moving, max_depth)
    kernel, rays, start, end, intensity, variance, mask = case
    grouped = _as_mapping(
        kernel.accumulate_group(intensity, variance, mask, rays, start, end)
    )
    expected = _as_mapping(_per_frame_reference(*case))

    assert len(expected) > 1, "geometry must reach more than one voxel"
    assert set(grouped) == set(expected)
    for key, (found, variance_, weight, contributors) in grouped.items():
        reference = expected[key]
        assert contributors == reference[3], key
        assert found == pytest.approx(reference[0], rel=1e-12), key
        assert variance_ == pytest.approx(reference[1], rel=1e-12), key
        assert weight == pytest.approx(reference[2], rel=1e-12), key


def test_single_frame_group_reaches_the_same_voxels_as_accumulate():
    """A group of one is the fallback for scans a brick cannot help.

    It is scientifically identical to the per-image path but **not**
    bit-for-bit: a brick partitions the detector into rectangles where a
    block partitions it into runs of the flattened image, so pixels group
    into per-block accumulators differently and their sums associate
    differently. Every voxel and every contributor count still agrees
    exactly; the totals agree to rounding.
    """
    for max_depth in (0, 1, 2, 3):
        case = _case(1, False, max_depth)
        kernel, rays, start, end, intensity, variance, mask = case
        grouped = kernel.accumulate_group(
            intensity, variance, mask, rays, start, end
        )
        single = kernel.accumulate(
            np.ascontiguousarray(intensity[0]),
            np.ascontiguousarray(variance[0]),
            np.ascontiguousarray(mask[0]),
            rays,
            np.ascontiguousarray(start[0]),
            np.ascontiguousarray(end[0]),
        )
        for name in ("chunk_id", "local_voxel_id", "contributors"):
            np.testing.assert_array_equal(grouped[name], single[name], name)
        for name in ("weighted_intensity", "weighted_variance", "weight"):
            np.testing.assert_allclose(
                grouped[name], single[name], rtol=1e-14, err_msg=name
            )


def _partial_grid_case(frames, moving, max_depth, block=64):
    """A grid covering only a corner of what the detector reaches.

    The brick reject needs both kinds of brick present to be worth
    testing: some that provably miss the grid and some that do not.
    """
    rays = _rays()
    start, end = _angles(frames, moving)
    minimum, step, shape = _grid(rays, start, end)
    # Keep the low corner of the full footprint, a quarter of it per
    # axis, so most of the detector maps outside.
    partial_shape = np.maximum(1, shape // 4)
    kernel = _kernel(minimum, step, partial_shape, max_depth, 1, block)
    intensity, variance, mask = _payload(frames)
    return kernel, rays, start, end, intensity, variance, mask


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [0, 2, 3])
def test_brick_reject_skips_bricks_without_changing_the_answer(
    moving, max_depth
):
    """Rejecting a whole brick must be invisible in the output.

    A brick is proved to miss the grid by the same test a pixel is, one
    level up, so a rejected brick is one whose pixels would each have
    been rejected in turn. Skipping it is therefore an optimization with
    no scientific content -- which only means anything if it is checked
    against the path that does not have it. ``accumulate`` partitions the
    detector into runs of the flattened image rather than bricks and has
    no brick reject, so it is an independent reference here.
    """
    case = _partial_grid_case(4, moving, max_depth)
    kernel, rays, start, end, intensity, variance, mask = case

    grouped = kernel.accumulate_group(
        intensity, variance, mask, rays, start, end, True
    )
    profile = grouped["_profile"]
    bricks_total = profile["skipped_bricks"] + 1
    assert profile["skipped_bricks"] > 0, (
        "the partial grid must leave some brick provably outside it, "
        "or this test proves nothing"
    )
    assert grouped["chunk_id"].size > 0, (
        "and some brick must still reach it"
    )
    assert bricks_total > 1

    reference = _per_frame_reference(
        kernel, rays, start, end, intensity, variance, mask
    )
    found = _as_mapping(grouped)
    expected = _as_mapping(reference)
    assert set(found) == set(expected)
    for key, (weighted, variance_, weight, contributors) in expected.items():
        got = found[key]
        assert got[3] == contributors, key
        assert got[0] == pytest.approx(weighted, rel=1e-12), key
        assert got[1] == pytest.approx(variance_, rel=1e-12), key
        assert got[2] == pytest.approx(weight, rel=1e-12), key


@pytest.mark.parametrize("moving", [False, True])
def test_frames_reach_grid_never_rejects_a_frame_that_maps_records(moving):
    """The frame filter's claim, checked against what mapping produces.

    ``frames_reach_grid`` decides from geometry alone, before any frame
    is read, so the pipeline can skip a frame's read, correction and
    kernel call outright. That is only safe because the claim is
    one-directional: a rejected frame provably maps nothing. An accepted
    frame may still map nothing, which costs only the work it would have
    cost anyway.

    A full rotation against a grid sized to one frame's footprint is the
    case the filter exists for -- most sample angles carry the detector
    nowhere near that volume.
    """
    rays = _rays()
    frames = 48
    start = np.zeros((frames, 4))
    start[:, 0] = 0.10471975511965978
    start[:, 1] = np.linspace(0.0, 2.0 * np.pi, frames, endpoint=False)
    end = start.copy()
    if moving:
        end[:, 1] += np.deg2rad(0.1)
    start = np.ascontiguousarray(start)
    end = np.ascontiguousarray(end)

    # A grid covering what the first frame reaches, and little else.
    minimum, step, shape = _grid(
        rays,
        np.ascontiguousarray(start[:1]),
        np.ascontiguousarray(end[:1]),
    )
    kernel = _kernel(minimum, step, shape, 2, 1, 256)

    reaches = kernel.frames_reach_grid(rays, start, end)
    assert reaches.shape == (frames,)
    assert reaches.any(), "the grid was built from a frame that reaches it"
    assert not reaches.all(), (
        "a full rotation against a one-frame grid must leave some frame "
        "provably outside it, or this test proves nothing"
    )

    intensity, variance, mask = _payload(1)
    mapped = 0
    for frame in range(frames):
        batch = kernel.accumulate_group(
            intensity,
            variance,
            mask,
            rays,
            np.ascontiguousarray(start[frame : frame + 1]),
            np.ascontiguousarray(end[frame : frame + 1]),
        )
        if batch["chunk_id"].size:
            mapped += 1
            assert reaches[frame], (
                f"frame {frame} mapped {batch['chunk_id'].size} records "
                "but the filter would have skipped it"
            )
    assert mapped > 0


@pytest.mark.parametrize("moving", [False, True])
def test_brick_reject_costs_a_frame_that_misses_almost_nothing(moving):
    """A frame that reaches no part of the grid must not be walked.

    This is the case the reject exists for: on a small volume swept
    through a full rotation about half the scan's frames reach it at all,
    and each of the rest used to pay four or eight coordinate evaluations
    for every one of its pixels to establish that. The cost is now one
    test per brick, so the evaluation count must fall to a small multiple
    of the brick count rather than of the pixel count.
    """
    rays = _rays()
    start, end = _angles(2, moving)
    minimum, step, shape = _grid(rays, start, end)
    # Far outside anything the detector maps to, on every axis at once.
    unreachable = np.asarray(minimum) + 1000.0 * np.asarray(
        step
    ) * np.asarray(shape)
    kernel = _kernel(unreachable, step, shape, 3, 1, 64)
    intensity, variance, mask = _payload(2)

    batch = kernel.accumulate_group(
        intensity, variance, mask, rays, start, end, True
    )
    profile = batch["_profile"]

    assert batch["chunk_id"].size == 0
    assert profile["skipped_bricks"] > 0
    assert profile["skipped_pixels"] == profile["pixels_seen"]
    # The whole point: proportional to bricks, not to pixels. A block of
    # 64 pixels over a 24x28 detector is a handful of bricks per frame,
    # so anything near the pixel count means the reject did not fire.
    assert profile["coordinate_evaluations"] < profile["pixels_seen"], (
        "a missing frame must not cost a per-pixel evaluation"
    )


@pytest.mark.parametrize("moving", [False, True])
def test_group_is_independent_of_threads_and_brick_size(moving):
    """Splitting the group across workers or bricks changes nothing."""
    reference = None
    for threads, block in ((1, 4096), (1, 64), (4, 128), (8, 32)):
        kernel, rays, start, end, intensity, variance, mask = _case(
            8, moving, 0, threads=threads, block=block
        )
        candidate = kernel.accumulate_group(
            intensity, variance, mask, rays, start, end
        )
        if reference is None:
            reference = candidate
            continue
        for name in ("chunk_id", "local_voxel_id", "contributors"):
            np.testing.assert_array_equal(
                candidate[name], reference[name], f"{name} {threads}/{block}"
            )
        for name in ("weighted_intensity", "weighted_variance", "weight"):
            np.testing.assert_allclose(
                candidate[name],
                reference[name],
                rtol=1e-12,
                err_msg=f"{name} {threads}/{block}",
            )


def test_group_emits_fewer_records_than_the_frames_do_apart():
    """The point of the exercise, on geometry that has the redundancy.

    Adjacent frames must collapse into shared voxels, so a group emits
    strictly fewer records than the same frames mapped separately -- while
    reaching exactly the same set of voxels, which the equivalence test
    above already pins.
    """
    case = _case(8, False, 0)
    kernel, rays, start, end, intensity, variance, mask = case
    grouped = kernel.accumulate_group(
        intensity, variance, mask, rays, start, end, True
    )
    apart = sum(
        kernel.accumulate(
            np.ascontiguousarray(intensity[frame]),
            np.ascontiguousarray(variance[frame]),
            np.ascontiguousarray(mask[frame]),
            rays,
            np.ascontiguousarray(start[frame]),
            np.ascontiguousarray(end[frame]),
            True,
        )["_profile"]["reduced_block_records"]
        for frame in range(intensity.shape[0])
    )
    brick_records = grouped["_profile"]["reduced_block_records"]
    assert brick_records < apart
    # The pre-merge saving is what reaches the checkpoint layer, so it is
    # worth pinning that it is a substantial fraction rather than a token.
    assert brick_records < 0.85 * apart


def _tile_arguments(rays, intensity, variance, mask, tile):
    """The same group, as whole frames plus a rectangle."""
    row_start, row_stop, column_start, column_stop = tile
    return (
        [np.ascontiguousarray(frame) for frame in intensity],
        [np.ascontiguousarray(frame) for frame in variance],
        [np.ascontiguousarray(frame) for frame in mask],
        np.ascontiguousarray(
            rays[row_start : row_stop + 1, column_start : column_stop + 1]
        ),
    )


@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize("max_depth", [0, 1, 2])
@pytest.mark.parametrize(
    "tile", [(0, ROWS, 0, COLUMNS), (5, 19, 3, 20), (0, 7, 21, COLUMNS)]
)
def test_tile_view_is_bit_for_bit_with_the_copied_tile(tile, max_depth, moving):
    """Reading a tile in place must be the identical computation.

    The mapping pipeline hands the kernel whole corrected frames and the
    rectangle to map, rather than gathering each tile into its own
    ``(frames, rows, columns)`` buffer first. That is a pure optimisation:
    the same values, in the same order, into the same accumulators. Not
    "agrees to rounding" -- bit for bit, because nothing about the
    arithmetic or its association has changed.
    """
    frames = 4
    case = _case(frames, moving, max_depth)
    kernel, rays, start, end, intensity, variance, mask = case
    row_start, row_stop, column_start, column_stop = tile
    selection = np.s_[:, row_start:row_stop, column_start:column_stop]
    tile_rays = np.ascontiguousarray(
        rays[row_start : row_stop + 1, column_start : column_stop + 1]
    )
    copied = kernel.accumulate_group(
        np.ascontiguousarray(intensity[selection]),
        np.ascontiguousarray(variance[selection]),
        np.ascontiguousarray(mask[selection]),
        tile_rays,
        start,
        end,
    )
    frame_intensity, frame_variance, frame_mask, _rays_again = _tile_arguments(
        rays, intensity, variance, mask, tile
    )
    viewed = kernel.accumulate_group_tile(
        frame_intensity,
        frame_variance,
        frame_mask,
        tile_rays,
        start,
        end,
        row_start,
        row_stop,
        column_start,
        column_stop,
    )
    assert copied["chunk_id"].size > 1
    for name in (
        "chunk_id",
        "local_voxel_id",
        "contributors",
        "weighted_intensity",
        "weighted_variance",
        "weight",
    ):
        np.testing.assert_array_equal(viewed[name], copied[name], name)


def test_tiles_partition_the_detector():
    """Row bands mapped separately and merged reach the whole frame's voxels.

    The tiling the mapping pipeline uses is a set of row bands, and their
    merged batches have to reproduce the ungrouped whole-frame call --
    which is what makes reading a band in place safe at the band's edges.
    """
    case = _case(4, True, 1)
    kernel, rays, start, end, intensity, variance, mask = case
    whole = kernel.accumulate_group(intensity, variance, mask, rays, start, end)
    bands = []
    for row_start in range(0, ROWS, 9):
        row_stop = min(row_start + 9, ROWS)
        frame_intensity, frame_variance, frame_mask, band_rays = _tile_arguments(
            rays, intensity, variance, mask, (row_start, row_stop, 0, COLUMNS)
        )
        bands.append(
            kernel.accumulate_group_tile(
                frame_intensity,
                frame_variance,
                frame_mask,
                band_rays,
                start,
                end,
                row_start,
                row_stop,
                0,
                COLUMNS,
            )
        )
    merged = _reduce_batches(bands)
    assert len(bands) > 1
    assert set(_as_mapping(merged)) == set(_as_mapping(whole))
    for name in ("contributors",):
        assert merged[name].sum() == whole[name].sum()
    for name in ("weighted_intensity", "weight"):
        np.testing.assert_allclose(
            merged[name].sum(), whole[name].sum(), rtol=1e-12, err_msg=name
        )


def test_tile_view_rejects_impossible_requests():
    kernel, rays, start, end, intensity, variance, mask = _case(4, False, 0)
    frames = [np.ascontiguousarray(frame) for frame in intensity]
    variances = [np.ascontiguousarray(frame) for frame in variance]
    masks = [np.ascontiguousarray(frame) for frame in mask]
    with pytest.raises(ValueError, match="outside the frame"):
        kernel.accumulate_group_tile(
            frames, variances, masks, rays, start, end, 0, ROWS + 1, 0, COLUMNS
        )
    with pytest.raises(ValueError, match="non-empty"):
        kernel.accumulate_group_tile(
            frames, variances, masks, rays, start, end, 4, 4, 0, COLUMNS
        )
    with pytest.raises(ValueError, match="one frame each"):
        kernel.accumulate_group_tile(
            frames, variances[:2], masks, rays, start, end, 0, ROWS, 0, COLUMNS
        )
    with pytest.raises(ValueError, match=r"\(rows \+ 1, columns \+ 1, 3\)"):
        kernel.accumulate_group_tile(
            frames, variances, masks, rays, start, end, 0, 6, 0, COLUMNS
        )


def test_group_rejects_mismatched_shapes():
    kernel, rays, start, end, intensity, variance, mask = _case(4, False, 0)
    with pytest.raises(ValueError, match="three-dimensional"):
        kernel.accumulate_group(
            intensity[0], variance[0], mask[0], rays, start, end
        )
    with pytest.raises(ValueError, match="match the intensity shape"):
        kernel.accumulate_group(
            intensity, variance[:2], mask, rays, start, end
        )
    with pytest.raises(ValueError, match=r"\(frames, 4\)"):
        kernel.accumulate_group(
            intensity, variance, mask, rays, start[:2], end
        )

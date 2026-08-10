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

"""Regression tests for native reciprocal-space reconstruction."""

from __future__ import annotations

import os
import sys
import time

import numpy as np
import pytest

from orgui.backend.scans import h5_Image
from orgui.app.config_data import CorrectionState
from orgui.app.mask_config import repair_intensity_variance
from orgui.datautils.xrayutils import HKLVlieg
from orgui.reconstruction_job import _correction_pipeline, derive_grid
from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointRouter,
    _GridSpec,
    _ReconstructionSpec,
    _build_kernels,
    _admissible_sample_pixels,
    _calibration_probe,
    _relative_standard_error,
    _representative_tile_origin,
    _sobol_tile_origins,
    _files_per_job,
    _kernel_threads_sweep,
    _map_frame_group,
    _merge_sorted_batches,
    _reduce_batches,
    _tile_ray_arrays,
)


native = pytest.importorskip("orgui.datautils.xrayutils._reciprocal_reconstruction_cpp")

#: The GitHub-hosted macOS runner overshoots short ``time.sleep`` calls badly
#: -- a 20 ms sleep was measured at 178 ms there -- which breaks timing
#: assertions that compare two independently measured sweeps. Only that runner
#: substitutes a fake clock; every other environment, local runs included,
#: keeps measuring real elapsed time.
_FAKE_SWEEP_CLOCK = (
    sys.platform == "darwin" and os.environ.get("GITHUB_ACTIONS") == "true"
)


def test_native_merge_sorted_batches():
    """Native merging preserves order and reduces matching voxel keys."""

    def batch(local, intensity, contributors):
        size = len(local)
        return {
            "chunk_id": np.zeros(size, dtype=np.uint64),
            "local_voxel_id": np.asarray(local, dtype=np.uint64),
            "weighted_intensity": np.asarray(intensity, dtype=np.float64),
            "weighted_variance": np.ones(size, dtype=np.float64),
            "weight": np.ones(size, dtype=np.float64),
            "contributors": np.asarray(contributors, dtype=np.uint64),
        }

    merged = _merge_sorted_batches(
        batch([1, 3], [10.0, 30.0], [1, 2]),
        batch([2, 3], [20.0, 4.0], [3, 4]),
    )

    np.testing.assert_array_equal(merged["local_voxel_id"], [1, 2, 3])
    np.testing.assert_array_equal(merged["weighted_intensity"], [10.0, 20.0, 34.0])
    np.testing.assert_array_equal(merged["weighted_variance"], [1.0, 1.0, 2.0])
    np.testing.assert_array_equal(merged["weight"], [1.0, 1.0, 2.0])
    np.testing.assert_array_equal(merged["contributors"], [1, 3, 6])


def test_grid_validates_effective_hdf5_chunk_bytes():
    """Requested chunks are clamped to grid shape before HDF5 validation."""
    valid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(514.0, 514.0, 514.0),
        step=(1.0, 1.0, 1.0),
        frame="hkl",
        chunk_shape=(1024, 1024, 1024),
    )
    assert valid.shape == (514, 514, 514)

    with pytest.raises(ValueError, match="smaller than 4 GiB"):
        _GridSpec(
            minimum=(0.0, 0.0, 0.0),
            maximum=(1024.0, 1024.0, 1024.0),
            step=(1.0, 1.0, 1.0),
            frame="hkl",
            chunk_shape=(1024, 1024, 1024),
        )


def test_grid_rejects_chunk_shape_that_overflows_uint32_local_index():
    """RecordKey.local is uint32 in the native kernel -- chunk_shape's own
    voxel product must stay under 2**32-1, or per-record indices would
    silently wrap. 1626**3 = 4,298,942,376 exceeds the bound; 1625**3 =
    4,291,015,625 does not.
    """
    with pytest.raises(ValueError, match="within-chunk voxel index"):
        _GridSpec(
            minimum=(0.0, 0.0, 0.0),
            maximum=(5000.0, 5000.0, 5000.0),
            step=(1.0, 1.0, 1.0),
            frame="hkl",
            chunk_shape=(1626, 1626, 1626),
        )


def test_grid_rejects_config_that_overflows_uint32_chunk_index():
    """RecordKey.chunk is uint32 too -- an unrealistically fine chunk_shape
    relative to a large grid can overflow the total chunk count even
    though each individual chunk is tiny (so the unrelated 4 GiB
    per-chunk check does not catch it). ceil(5000/3)**3 = 4,632,407,963
    exceeds the bound.
    """
    with pytest.raises(ValueError, match="uint32 chunk index"):
        _GridSpec(
            minimum=(0.0, 0.0, 0.0),
            maximum=(5000.0, 5000.0, 5000.0),
            step=(1.0, 1.0, 1.0),
            frame="hkl",
            chunk_shape=(3, 3, 3),
        )


def test_grid_accepts_default_chunk_shape_at_a_large_grid():
    """The default chunk_shape has wide uint32 headroom even at a
    5000**3-voxel grid (493,039 chunks, 262,144 voxels/chunk -- both far
    under 2**32-1)."""
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(5000.0, 5000.0, 5000.0),
        step=(1.0, 1.0, 1.0),
        frame="hkl",
    )
    assert grid.shape == (5000, 5000, 5000)
    assert grid.chunk_shape == (64, 64, 64)


def test_grid_rejects_shape_needing_more_than_64_bits_of_voxel_index():
    """The native kernel packs the three axis indices into one uint64 voxel
    identifier, one bit field per axis. 2,640,625 voxels per axis needs 22
    bits each, 66 in total -- while staying inside both uint32 RecordKey
    bounds, so only this check can catch it.
    """
    with pytest.raises(ValueError, match="64-bit packed voxel identifier"):
        _GridSpec(
            minimum=(0.0, 0.0, 0.0),
            maximum=(2640625.0, 2640625.0, 2640625.0),
            step=(1.0, 1.0, 1.0),
            frame="hkl",
            chunk_shape=(1625, 1625, 1625),
        )


def test_grid_accepts_a_realistic_grids_voxel_index_width():
    """A 5000**3 grid needs 13 bits per axis, 39 of the 64 available."""
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(5000.0, 5000.0, 5000.0),
        step=(1.0, 1.0, 1.0),
        frame="hkl",
    )
    assert sum(int(size - 1).bit_length() for size in grid.shape) == 39


@pytest.mark.parametrize(
    "chunk_shape",
    [
        pytest.param((8, 4, 16), id="power-of-two-chunks"),
        pytest.param((5, 7, 9), id="arbitrary-chunks"),
    ],
)
def test_native_chunk_local_ids_match_independently_computed_voxel_indices(
    chunk_shape,
):
    """End-to-end check of voxel_id()'s packed identifier and record_key()'s
    decomposition of it: every emitted (chunk_id, local_voxel_id) pair must
    decode back to the voxel the pixel's own coordinate falls in, computed
    here from ``coordinate()`` alone. Deliberately uses axes of differing
    extent (so the packed bit fields differ in width) and covers both the
    shift/mask chunk split and the division fallback for a chunk_shape that
    is not a power of two.
    """
    minimum = np.array([-0.6, -0.35, -0.45])
    step = np.array([0.02, 0.015, 0.025])
    shape = np.array([37, 61, 23], dtype=np.int64)
    kernel = native.ReconstructionKernel(
        minimum,
        step,
        shape,
        np.array(chunk_shape, dtype=np.int64),
        "lab",
        1.0,
        np.eye(3),
        np.eye(3),
        0,
        1,
        16,
        1024 * 1024,
    )
    rows = columns = 12
    rays = np.zeros((rows + 1, columns + 1, 3), dtype=np.float64)
    rays[..., 0] = np.linspace(-0.35, 0.35, columns + 1)[None, :]
    rays[..., 2] = np.linspace(-0.3, 0.3, rows + 1)[:, None]
    rays[..., 1] = 1.0
    rays /= np.linalg.norm(rays, axis=2, keepdims=True)
    angles = np.zeros(4)

    expected = set()
    chunk_grid = np.ceil(shape / np.array(chunk_shape)).astype(np.int64)
    for row in range(rows):
        for column in range(columns):
            coordinate = np.asarray(
                kernel.coordinate(rays, angles, angles, row, column)
            )
            index = np.floor((coordinate - minimum) / step).astype(np.int64)
            if np.any(index < 0) or np.any(index >= shape):
                continue
            chunk_index, local_index = np.divmod(index, np.array(chunk_shape))
            expected.add(
                (
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
            )

    result = kernel.accumulate(
        np.ones((rows, columns)),
        np.ones((rows, columns)),
        np.zeros((rows, columns), dtype=bool),
        rays,
        angles,
        angles,
    )

    assert expected
    assert len(expected) > 1
    assert (
        set(
            zip(
                result["chunk_id"].tolist(),
                result["local_voxel_id"].tolist(),
            )
        )
        == expected
    )


def test_native_center_accumulation_survives_revisiting_voxels_out_of_order():
    """The mapping loop remembers the previous pixel's accumulator, which is
    only correct if a key revisited after an intervening different key still
    lands on its own accumulator. This tile alternates between voxels so
    that cache hit and miss alternate too, and checks every record's summed
    values against the same sums accumulated independently here.
    """
    minimum = np.array([-0.4, -0.4, -0.4])
    step = np.array([0.05, 0.05, 0.05])
    shape = np.array([16, 16, 16], dtype=np.int64)
    kernel = native.ReconstructionKernel(
        minimum,
        step,
        shape,
        np.array([8, 8, 8], dtype=np.int64),
        "lab",
        1.0,
        np.eye(3),
        np.eye(3),
        0,
        1,
        4096,
        1024 * 1024,
    )
    rows = 4
    columns = 12
    rays = np.zeros((rows + 1, columns + 1, 3), dtype=np.float64)
    # Alternating corner offsets, so consecutive pixel centres keep moving
    # back and forth between a small set of voxels.
    offsets = np.array([0.0, 0.0, 0.22, 0.22] * ((columns + 4) // 4))
    rays[..., 0] = offsets[: columns + 1][None, :]
    rays[..., 2] = np.linspace(-0.05, 0.05, rows + 1)[:, None]
    rays[..., 1] = 1.0
    rays /= np.linalg.norm(rays, axis=2, keepdims=True)
    rng = np.random.default_rng(7)
    intensity = rng.random((rows, columns))
    variance = rng.random((rows, columns)) + 0.5
    mask = np.zeros((rows, columns), dtype=bool)
    angles = np.zeros(4)

    expected: dict[tuple[int, int], list[float]] = {}
    chunk_shape = np.array([8, 8, 8])
    chunk_grid = np.ceil(shape / chunk_shape).astype(np.int64)
    for row in range(rows):
        for column in range(columns):
            coordinate = np.asarray(
                kernel.coordinate(rays, angles, angles, row, column)
            )
            index = np.floor((coordinate - minimum) / step).astype(np.int64)
            if np.any(index < 0) or np.any(index >= shape):
                continue
            chunk_index, local_index = np.divmod(index, chunk_shape)
            key = (
                int(
                    (chunk_index[0] * chunk_grid[1] + chunk_index[1]) * chunk_grid[2]
                    + chunk_index[2]
                ),
                int(
                    (local_index[0] * chunk_shape[1] + local_index[1]) * chunk_shape[2]
                    + local_index[2]
                ),
            )
            totals = expected.setdefault(key, [0.0, 0.0, 0.0, 0])
            totals[0] += intensity[row, column]
            totals[1] += variance[row, column]
            totals[2] += 1.0
            totals[3] += 1

    result = kernel.accumulate(intensity, variance, mask, rays, angles, angles)

    # More voxels than one, revisited more often than they are distinct --
    # otherwise the cache would never be exercised in both directions.
    assert len(expected) > 2
    assert result["chunk_id"].size == len(expected)
    assert sum(totals[3] for totals in expected.values()) > 2 * len(expected)
    for position in range(result["chunk_id"].size):
        key = (
            int(result["chunk_id"][position]),
            int(result["local_voxel_id"][position]),
        )
        totals = expected[key]
        assert result["weighted_intensity"][position] == totals[0]
        assert result["weighted_variance"][position] == totals[1]
        assert result["weight"][position] == totals[2]
        assert int(result["contributors"][position]) == totals[3]


def test_native_accumulate_produces_uint32_chunk_and_local_ids_across_multiple_chunks():
    """A grid split into several chunks along one axis exercises the real
    record_key() chunk/local decomposition (not just the trivial
    single-chunk case), confirming uint32 output end to end."""
    kernel = native.ReconstructionKernel(
        np.array([-20.0, -20.0, -20.0]),
        np.array([0.01, 0.01, 0.01]),
        np.array([4000, 4000, 4000], dtype=np.int64),
        np.array([16, 4000, 4000], dtype=np.int64),
        "lab",
        1.0,
        np.linalg.inv(np.eye(3)),
        np.eye(3),
        0,
        1,
        16,
        1024 * 1024,
    )
    rays = np.zeros((4, 4, 3), dtype=np.float64)
    rays[..., 0] = np.linspace(-1.5, 1.5, 4)[:, None]
    rays[..., 1] = 1.0
    result = kernel.accumulate(
        np.ones((3, 3)),
        np.ones((3, 3)),
        np.zeros((3, 3), dtype=bool),
        rays,
        np.zeros(4),
        np.zeros(4),
    )
    assert result["chunk_id"].dtype == np.uint32
    assert result["local_voxel_id"].dtype == np.uint32
    assert result["contributors"].dtype == np.uint32
    assert result["chunk_id"].size > 0
    # This tile spans multiple 16-voxel-wide chunks along axis 0 -- confirm
    # record_key() actually exercised more than one chunk, not just the
    # trivial chunk=0 case every other test in this file uses.
    assert np.unique(result["chunk_id"]).size > 1


def test_native_accumulate_reused_arena_reduce_is_thread_count_and_depth_independent():
    """accumulate_block's per-worker-thread reused pmr::map arena replaces
    the old per-block vector+stable_sort reduce -- this must not change
    the scientific result. Cross-checks identical (record count, summed
    weighted_intensity, summed contributors) across thread counts (1, 4,
    8 -- exercising a single worker reusing its arena across many blocks
    as well as several workers each with their own arena) and across
    max_depth values (0 through 4, spanning the max_depth==0 fast path,
    the depth==2 stationary fast path, and the general adaptive
    split_pixel path -- all three feed the same shared reduce)."""
    rows, columns = 64, 64
    rng_intensity = np.random.default_rng(0)
    rng_variance = np.random.default_rng(1)
    intensity = rng_intensity.uniform(1.0, 100.0, size=(rows, columns))
    variance = np.abs(rng_variance.normal(1.0, 0.1, size=(rows, columns)))
    mask = np.zeros((rows, columns), dtype=bool)
    rays = _constant_rays(rows, columns)
    rays[..., 0] = np.linspace(-3.0, 3.0, columns + 1)[None, :]

    for max_depth in (0, 1, 2, 3, 4):
        reference = None
        for threads in (1, 4, 8):
            kernel = native.ReconstructionKernel(
                np.array([-20.0, -20.0, -20.0]),
                np.array([0.02, 0.02, 0.02]),
                np.array([2000, 2000, 2000], dtype=np.int64),
                np.array([32, 32, 32], dtype=np.int64),
                "lab",
                1.0,
                np.linalg.inv(np.eye(3)),
                np.eye(3),
                max_depth,
                threads,
                64,
                512 * 1024 * 1024,
            )
            result = kernel.accumulate(
                intensity, variance, mask, rays, np.zeros(4), np.zeros(4)
            )
            chunk = result["chunk_id"]
            local = result["local_voxel_id"]
            combined_key = chunk.astype(np.uint64) << 32 | local.astype(np.uint64)
            assert np.all(np.diff(combined_key) > 0), (
                f"max_depth={max_depth} threads={threads}: output must be "
                "strictly sorted with no duplicate keys"
            )
            summary = (
                chunk.size,
                float(np.sum(result["weighted_intensity"])),
                int(np.sum(result["contributors"])),
            )
            if reference is None:
                reference = summary
            else:
                assert summary[0] == reference[0], (max_depth, threads)
                assert summary[1] == pytest.approx(reference[1], rel=1e-9), (
                    max_depth,
                    threads,
                )
                assert summary[2] == reference[2], (max_depth, threads)


def test_native_loser_tree_merges_duplicate_voxels_across_many_blocks():
    """A hand-constructed scenario where every valid pixel maps to the
    SAME single voxel, forced across many separate blocks (small
    work_block_pixels relative to the tile) -- exercises exactly the
    cross-block duplicate-key case the loser-tree merge must combine
    correctly (accumulate_block's own per-block reduce cannot produce
    this duplication on its own; only the Level-2 merge across blocks
    can). The expected combined result is exactly predictable: one
    output record, contributors equal to the valid pixel count,
    weighted_intensity equal to the sum of every valid pixel's own
    intensity (weight is 1.0 per pixel at max_depth=0)."""
    rows, columns = 20, 20
    intensity = np.arange(1, rows * columns + 1, dtype=np.float64).reshape(
        rows, columns
    )
    variance = np.ones((rows, columns), dtype=np.float64)
    mask = np.zeros((rows, columns), dtype=bool)
    mask[0, 0] = True  # one excluded pixel, to also confirm it is excluded
    rays = _constant_rays(rows, columns)  # every corner ray identical -> one voxel

    kernel = native.ReconstructionKernel(
        np.array([-20.0, -20.0, -20.0]),
        np.array([0.01, 0.01, 0.01]),
        np.array([4000, 4000, 4000], dtype=np.int64),
        np.array([64, 64, 64], dtype=np.int64),  # same default as elsewhere here
        "lab",
        1.0,
        np.linalg.inv(np.eye(3)),
        np.eye(3),
        0,  # max_depth
        4,  # threads
        8,  # work_block_pixels -- forces (400 valid pixels)/8 = 50 blocks
        512 * 1024 * 1024,
    )
    result = kernel.accumulate(
        intensity, variance, mask, rays, np.zeros(4), np.zeros(4)
    )

    valid_pixels = intensity.size - int(np.sum(mask))
    expected_intensity = float(np.sum(intensity)) - float(intensity[0, 0])

    assert result["chunk_id"].size == 1
    assert int(result["contributors"][0]) == valid_pixels
    assert result["weighted_intensity"][0] == pytest.approx(expected_intensity)
    # weight is 1.0 per contributing pixel at max_depth=0.
    assert result["weight"][0] == pytest.approx(float(valid_pixels))


def test_files_per_job_formula_floors_at_checkpoint_count():
    """The checkpoint-count floor wins when data comfortably fits budget."""
    assert _files_per_job(0, 24e9, checkpoint_count=10) == 10
    assert _files_per_job(6e9, 24e9, checkpoint_count=10) == 10
    # Right at the boundary: ceil(24e9 / 24e9) == 1, still below the floor.
    assert _files_per_job(24e9, 24e9, checkpoint_count=10) == 10


def test_files_per_job_formula_scales_with_data_volume():
    """The data-volume term wins once it exceeds the checkpoint-count floor."""
    assert _files_per_job(129.7e9, 4e9, checkpoint_count=10) == 33
    assert _files_per_job(129.7e9, 8e9, checkpoint_count=10) == 17
    assert _files_per_job(129.7e9, 13e9, checkpoint_count=10) == 10
    # Exact multiples must not be rounded up unnecessarily.
    assert _files_per_job(100.0, 10.0, checkpoint_count=1) == 10


def test_files_per_job_formula_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="checkpoint_count"):
        _files_per_job(1.0, 1.0, checkpoint_count=0)
    with pytest.raises(ValueError, match="ram_budget_bytes"):
        _files_per_job(1.0, 0.0, checkpoint_count=1)
    with pytest.raises(ValueError, match="job_data_bytes"):
        _files_per_job(-1.0, 1.0, checkpoint_count=1)


def _kernel(frame="lab", *, threads=1, max_depth=2, ub=None):
    if ub is None:
        ub = np.eye(3)
    return native.ReconstructionKernel(
        np.array([-20.0, -20.0, -20.0]),
        np.array([0.01, 0.01, 0.01]),
        np.array([4000, 4000, 4000], dtype=np.int64),
        np.array([64, 64, 64], dtype=np.int64),
        frame,
        1.0,
        np.linalg.inv(ub),
        np.eye(3),
        max_depth,
        threads,
        16,
        1024 * 1024,
    )


def _constant_rays(rows, columns):
    rays = np.zeros((rows + 1, columns + 1, 3), dtype=np.float64)
    rays[..., 1] = 1.0
    return rays


def test_native_average_variance_for_two_pixels():
    kernel = _kernel()
    result = kernel.accumulate(
        np.array([[10.0, 20.0]]),
        np.array([[10.0, 20.0]]),
        np.zeros((1, 2), dtype=bool),
        _constant_rays(1, 2),
        np.zeros(4),
        np.zeros(4),
    )

    assert result["chunk_id"].size == 1
    assert result["weighted_intensity"][0] == 30.0
    assert result["weighted_variance"][0] == 30.0
    assert result["weight"][0] == 2.0
    assert result["contributors"][0] == 2
    assert result["weighted_intensity"][0] / result["weight"][0] == 15.0
    assert result["weighted_variance"][0] / result["weight"][0] ** 2 == 7.5


def test_arbitrary_rectangular_tiles_match_full_detector():
    """Shared corner rays make arbitrary tile boundaries numerically exact."""
    rows, columns = 3, 5
    yy, xx = np.meshgrid(
        np.arange(rows + 1),
        np.arange(columns + 1),
        indexing="ij",
    )
    rays = np.stack(
        (
            (xx - columns / 2) * 0.02,
            np.ones_like(xx),
            (yy - rows / 2) * 0.02,
        ),
        axis=-1,
    )
    rays /= np.linalg.norm(rays, axis=-1, keepdims=True)
    intensity = np.arange(1, rows * columns + 1, dtype=np.float64).reshape(
        rows, columns
    )
    variance = intensity + 0.5
    mask = np.zeros((rows, columns), dtype=bool)
    kernel = _kernel(max_depth=2)
    full = kernel.accumulate(
        intensity,
        variance,
        mask,
        rays,
        np.zeros(4),
        np.zeros(4),
    )
    tiles = (
        (0, 1, 0, 2),
        (0, 1, 2, 5),
        (1, 3, 0, 2),
        (1, 3, 2, 5),
    )
    tiled = _reduce_batches(
        kernel.accumulate(
            intensity[row_start:row_stop, column_start:column_stop],
            variance[row_start:row_stop, column_start:column_stop],
            mask[row_start:row_stop, column_start:column_stop],
            rays[
                row_start : row_stop + 1,
                column_start : column_stop + 1,
            ],
            np.zeros(4),
            np.zeros(4),
        )
        for row_start, row_stop, column_start, column_stop in tiles
    )

    for name in ("chunk_id", "local_voxel_id", "contributors"):
        np.testing.assert_array_equal(tiled[name], full[name])
    for name in (
        "weighted_intensity",
        "weighted_variance",
        "weight",
    ):
        np.testing.assert_allclose(tiled[name], full[name], rtol=2e-15)
    assert tiled["weight"].sum() == pytest.approx(rows * columns)


def test_native_xxh3_128_matches_known_vector():
    assert (
        native.xxh3_128(np.empty(0, dtype=np.uint8))
        == "99aa06d3014798d86001c324468d497f"
    )


def test_shared_correction_pipeline_propagates_factor_uncertainty():
    scan = _FakeScan([10.0])
    scan.exposure_time = np.array([2.0])
    scan.exposure_time_variance = np.array([0.25])
    config = type(
        "Config",
        (),
        {
            "corrections": CorrectionState(
                use_background=True, normalize_exposure=True
            ),
            "detector": object(),
        },
    )()
    provenance = {}
    correct = _correction_pipeline(
        config,
        scan,
        {
            "background": np.array([[2.0]]),
            "background_variance": np.array([[1.0]]),
        },
        provenance,
    )

    intensity, variance, mask = correct(
        h5_Image(np.array([[10.0]])),
        np.array([[10.0]]),
        0,
        (0, 1, 0, 1),
    )

    assert intensity[0, 0] == 4.0
    assert variance[0, 0] == 3.75
    assert not mask[0, 0]
    assert provenance["factor_uncertainty"]["exposure"] == "propagated"


class _CorrectionScan:
    """A scan carrying an exposure time and one monitor, with variances."""

    def __init__(self, image):
        self.image = image
        self.exposure_time = np.array([0.37])
        self.exposure_time_variance = np.array([0.011])
        self.ic2 = np.array([1.7e5])
        self.ic2_variance = np.array([913.0])

    def __len__(self):
        return 1

    def get_raw_img(self, index):
        return h5_Image(self.image)


class _SolidAngleDetector:
    """Just enough detector for the static per-pixel correction factor."""

    _polFactor = 0.93
    _polAxis = 0.11
    _deltaChi = 1.5707963267948966

    def __init__(self, shape):
        self.shape = shape

    def solidAngleArray(self):
        rows, columns = self.shape
        grid = np.arange(rows * columns, dtype=np.float64).reshape(rows, columns)
        return 1.0 + grid / (rows * columns)

    def polarization(self, factor, axis_offset):
        rows, columns = self.shape
        grid = np.arange(rows * columns, dtype=np.float64).reshape(rows, columns)
        return 0.5 + factor * 0.25 + axis_offset + grid / (2 * rows * columns)

    def polarizationArray(self, shape=None):
        """Mirror the conversion Detector2D_SXRD applies before pyFAI."""
        return self.polarization(
            factor=2.0 * self._polFactor - 1.0,
            axis_offset=self._deltaChi - self._polAxis,
        )


@pytest.mark.parametrize("propagate", [False, True])
def test_native_correction_is_bit_for_bit_with_the_numpy_form(monkeypatch, propagate):
    """The fused native pass must not move a single bit.

    True under ``-ffp-contract=off``, which the ``strict_fp_contract`` build
    option sets and the test workflow enables. A compiler that fuses
    ``spread += value * value * factor_variance`` into an ``fma`` (Apple clang
    on arm64) rounds once where NumPy rounds twice and moves the variance by
    one ULP; release wheels are built that way deliberately.

    Correction applies a per-pixel factor and then each scalar factor in
    turn, and NumPy does it in eight or nine full-detector passes. Doing
    the same arithmetic one pixel at a time in the extension is a pure
    memory-traffic change: every operation is element-wise, so there is
    no reduction to reassociate and no reason for a result to move.

    Both variance branches are exercised, because they differ in *order*
    -- a propagated factor scales the variance, then uses the intensity
    from before the factor was applied, then scales the intensity -- and
    an implementation that got that order wrong would still look
    plausible on the deterministic branch.
    """
    import orgui.reconstruction_job as job_module

    rows, columns = 9, 11
    rng = np.random.default_rng(20260811)
    image = rng.uniform(0.0, 5000.0, size=(rows, columns))
    # A pixel that must come back masked whatever else happens.
    image[4, 5] = np.inf
    scan = _CorrectionScan(image)
    if not propagate:
        del scan.exposure_time_variance
        del scan.ic2_variance
    config = type(
        "Config",
        (),
        {
            "corrections": CorrectionState(
                use_mask=True,
                use_solid_angle=True,
                use_polarization=True,
                normalize_exposure=True,
                monitor_corrections=("ic2",),
            ),
            "detector": _SolidAngleDetector((rows, columns)),
        },
    )()
    static_mask = np.zeros((rows, columns), dtype=bool)
    static_mask[0, 0] = True
    assets = {"mask": static_mask}

    def corrected():
        pipeline = _correction_pipeline(config, scan, assets, {})
        return pipeline.correct_frame(h5_Image(image), image, 0)

    native_intensity, native_variance, native_mask = corrected()
    monkeypatch.setattr(job_module, "_correction_extension", lambda: None)
    numpy_intensity, numpy_variance, numpy_mask = corrected()

    assert job_module._correction_extension() is None
    np.testing.assert_array_equal(native_intensity, numpy_intensity)
    np.testing.assert_array_equal(native_variance, numpy_variance)
    np.testing.assert_array_equal(native_mask, numpy_mask)
    # ... and the comparison has to be worth something: the factors must
    # actually have changed the values, and the non-finite pixel must be
    # masked while its neighbours are not.
    assert not np.array_equal(native_intensity, image)
    assert native_mask[4, 5]
    assert native_mask[0, 0]
    assert native_mask.sum() == 2


def test_shared_correction_repairs_across_detector_tile_boundaries():
    scan = _FakeScan([0.0])
    config = type(
        "Config",
        (),
        {
            "corrections": CorrectionState(
                use_mask=True,
                repair_masked_pixels=True,
                normalize_exposure=False,
            ),
            "detector": object(),
        },
    )()
    image = np.arange(49, dtype=np.float64).reshape(7, 7)
    mask = np.zeros((7, 7), dtype=bool)
    mask[3, 3] = True
    payload = h5_Image(image)
    correct = _correction_pipeline(
        config,
        scan,
        {
            "mask": mask,
            "repair": {
                "max_component_pixels": 4,
                "max_span": 3,
                "radius": 2,
                "min_valid_neighbors": 6,
                "row_gaps": np.empty((0, 2), dtype=np.int32),
                "column_gaps": np.empty((0, 2), dtype=np.int32),
            },
        },
        {},
    )

    full = correct.correct_frame(payload, image, 0)
    left = correct(payload, image[:, :4], 0, (0, 7, 0, 4))
    right = correct(payload, image[:, 4:], 0, (0, 7, 4, 7))

    for index in range(3):
        combined = np.concatenate((left[index], right[index]), axis=1)
        np.testing.assert_array_equal(combined, full[index])
    assert not full[2][3, 3]


def test_pixel_repair_propagates_interpolation_weights():
    intensity = np.arange(9, dtype=np.float64).reshape(3, 3)
    variance = np.ones((3, 3), dtype=np.float64)
    mask = np.zeros((3, 3), dtype=bool)
    mask[1, 1] = True

    _, repaired_variance, remaining, repaired = repair_intensity_variance(
        intensity,
        variance,
        mask,
        max_component_pixels=2,
        max_span=2,
        radius=2,
        min_valid_neighbors=4,
    )

    assert repaired[1, 1]
    assert not remaining[1, 1]
    assert repaired_variance[1, 1] == 0.5


def test_native_coordinate_frames_match_vlieg():
    lattice = HKLVlieg.Lattice([3.9, 4.1, 5.2], [90.0, 95.0, 110.0])
    calculator = HKLVlieg.UBCalculator(lattice, 70.0)
    u = HKLVlieg.Rotation.from_euler("xyz", [0.17, -0.08, 0.11]).as_matrix()
    calculator.setU(u)
    position = np.deg2rad([0.3, 7.0, 3.0, -4.0, 2.0, 5.0])
    primary = HKLVlieg.primBeamAngles(position)
    gamma_p = primary[2]
    delta_p = primary[1]
    ray = np.array(
        [
            np.sin(delta_p) * np.cos(gamma_p),
            np.cos(delta_p) * np.cos(gamma_p),
            np.sin(gamma_p),
        ]
    )
    rays = np.broadcast_to(ray, (2, 2, 3)).copy()
    sample_angles = position[[0, 3, 4, 5]]
    q_alpha = np.asarray(
        HKLVlieg.VliegAngles(calculator).QAlpha(*position[:3])
    ).reshape(3)
    matrices = HKLVlieg.createVliegMatrices(position)
    q_omega = matrices[3].T @ q_alpha
    q_chi = matrices[4].T @ q_omega
    q_phi = matrices[5].T @ q_chi
    q_lab = matrices[0] @ q_alpha
    expected = {
        "lab": q_lab,
        "alpha": q_alpha,
        "omega": q_omega,
        "chi": q_chi,
        "phi": q_phi,
        "crystal": np.linalg.inv(calculator.getU()) @ q_phi,
        "hkl": np.linalg.inv(calculator.getUB()) @ q_phi,
    }

    for frame, coordinates in expected.items():
        kernel = native.ReconstructionKernel(
            np.full(3, -100.0),
            np.ones(3),
            np.full(3, 200, dtype=np.int64),
            np.full(3, 16, dtype=np.int64),
            frame,
            calculator.getK(),
            np.linalg.inv(calculator.getUB()),
            np.linalg.inv(calculator.getU()),
        )
        actual = kernel.coordinate(rays, sample_angles, sample_angles, 0, 0)
        assert np.allclose(actual, coordinates, rtol=1e-13, atol=1e-13), frame


def test_native_coordinate_interpolates_detector_corner_rays():
    rays = np.array(
        [
            [[-0.2, 0.97, -0.1], [0.3, 0.93, -0.15]],
            [[-0.25, 0.94, 0.2], [0.35, 0.9, 0.25]],
        ],
        dtype=np.float64,
    )
    rays /= np.linalg.norm(rays, axis=2, keepdims=True)
    u = 0.37
    v = 0.61
    interpolated = (
        (1.0 - u) * (1.0 - v) * rays[0, 0]
        + u * (1.0 - v) * rays[1, 0]
        + (1.0 - u) * v * rays[0, 1]
        + u * v * rays[1, 1]
    )
    interpolated /= np.linalg.norm(interpolated)
    expected = interpolated - np.array([0.0, 1.0, 0.0])

    actual = _kernel(frame="lab").coordinate(
        rays,
        np.zeros(4),
        np.array([0.1, -0.2, 0.3, -0.4]),
        0,
        0,
        u,
        v,
        0.73,
    )

    assert np.allclose(actual, expected, rtol=1e-14, atol=1e-14)


def test_native_results_do_not_depend_on_thread_count():
    rng = np.random.default_rng(12)
    intensity = rng.random((8, 8))
    variance = intensity.copy()
    mask = np.zeros_like(intensity, dtype=bool)
    rays = _constant_rays(8, 8)
    serial = _kernel(threads=1).accumulate(
        intensity, variance, mask, rays, np.zeros(4), np.zeros(4)
    )
    threaded = _kernel(threads=4).accumulate(
        intensity, variance, mask, rays, np.zeros(4), np.zeros(4)
    )

    for name in serial:
        assert np.array_equal(serial[name], threaded[name]), name


def test_center_only_profiles_one_coordinate_per_valid_pixel():
    intensity = np.ones((2, 3), dtype=np.float64)
    variance = np.ones_like(intensity)
    mask = np.zeros_like(intensity, dtype=bool)
    mask[0, 1] = True

    result = _kernel(max_depth=0).accumulate(
        intensity,
        variance,
        mask,
        _constant_rays(2, 3),
        np.zeros(4),
        np.zeros(4),
        profile=True,
    )

    profile = result.pop("_profile")
    assert profile["valid_pixels"] == 5
    assert profile["coordinate_evaluations"] == 5
    assert profile["maximum_weights_per_pixel"] == 1
    # Depth zero skips the per-pixel weight vector entirely, so these two
    # counters are maintained by hand on that path and cannot be derived
    # from a weights.size() the fast path never computes.
    assert profile["voxel_weights"] == 5
    assert profile["unreduced_block_records"] == 5


def test_calibration_probe_returns_positive_byte_estimate():
    """The probe samples real pixels and reports positive per-pixel rates."""
    rows, columns = 40, 40
    mask = np.zeros((rows, columns), dtype=bool)
    mask[0, 0] = True  # one excluded pixel, well away from every tile drawn

    result = _calibration_probe(
        _kernel(max_depth=0, threads=1),
        mask,
        _constant_rays(rows, columns),
        np.zeros(4),
        np.zeros(4),
        budget_seconds=0.05,
        max_sample_pixels=2000,
        rng=np.random.default_rng(0),
    )

    assert result["kernel_threads"] == 1
    assert result["sampled_pixels"] > 0
    assert result["seconds_per_pixel"] >= 0.0
    assert result["records_per_pixel"] >= 0.0


def test_calibration_probe_respects_the_kernels_own_tile_limit():
    """The kernel refuses a tile whose worst-case footprint exceeds its
    memory budget, and that limit shrinks by the subdivision factor at
    every adaptive depth. The probe sizes its own samples, so it has to
    honour the same bound -- otherwise a deep-subdivision job raises
    inside what is a live, GUI-driven estimate rather than returning one.
    """
    rows, columns = 64, 64
    mask = np.zeros((rows, columns), dtype=bool)
    # This helper's 1 MiB budget makes the admissible tile far smaller
    # than the sample the probe would otherwise ask for, the same way a
    # high max_depth does against a real job's budget.
    kernel = _kernel(max_depth=2, threads=1)
    admissible = _admissible_sample_pixels(kernel, np.zeros(4), np.zeros(4))
    assert admissible < rows * columns

    result = _calibration_probe(
        kernel,
        mask,
        _constant_rays(rows, columns),
        np.zeros(4),
        np.zeros(4),
        budget_seconds=0.05,
        rng=np.random.default_rng(0),
    )

    assert result["sampled_pixels"] > 0


def test_sobol_tile_origins_cover_every_quadrant():
    """The whole point of sampling low-discrepancy positions rather than
    independent uniform ones is that a handful of tiles still reaches the
    whole detector. Eight origins must touch all four quadrants; eight
    independent draws frequently do not.
    """
    origins = _sobol_tile_origins(1000, 1000, 8, rng=np.random.default_rng(0))

    assert len(origins) == 8
    assert all(0.0 <= row <= 1.0 and 0.0 <= column <= 1.0 for row, column in origins)
    quadrants = {(row >= 0.5, column >= 0.5) for row, column in origins}
    assert len(quadrants) == 4


def test_sobol_tile_origins_are_reproducible_from_a_seed():
    """A pinned generator has to pin the sequence, so a job's plan can be
    reproduced."""
    first = _sobol_tile_origins(500, 500, 4, rng=np.random.default_rng(7))
    second = _sobol_tile_origins(500, 500, 4, rng=np.random.default_rng(7))

    assert first == second


def test_relative_standard_error_reports_convergence():
    """Undefined below two samples, zero when every sample agrees, and
    shrinking as agreeing samples accumulate."""
    assert _relative_standard_error([1.0]) == float("inf")
    assert _relative_standard_error([2.0, 2.0, 2.0]) == 0.0
    assert _relative_standard_error([1.0, 2.0, 3.0]) > _relative_standard_error(
        [1.9, 2.0, 2.1]
    )


def test_calibration_probe_reports_its_own_uncertainty():
    """The probe returns how well converged its density estimate is, so a
    caller can size a safety margin on measured scatter instead of a
    fixed guess."""
    rows, columns = 48, 48
    mask = np.zeros((rows, columns), dtype=bool)

    result = _calibration_probe(
        _kernel(max_depth=0, threads=1),
        mask,
        _constant_rays(rows, columns),
        np.zeros(4),
        np.zeros(4),
        budget_seconds=0.05,
        max_sample_pixels=4000,
        rng=np.random.default_rng(0),
    )

    assert result["sampled_tiles"] >= 2
    # Constant rays put every pixel in one voxel, so every tile agrees and
    # the estimate is exactly converged.
    assert result["records_per_pixel_relative_error"] < 1.0


def test_calibration_probe_bootstrap_avoids_masked_regions():
    """The bootstrap tile sets the rate the whole sized pass is derived
    from, so it must not land somewhere unrepresentative. A masked centre
    has to push it elsewhere rather than measuring excluded pixels.
    """
    rows, columns = 64, 64
    mask = np.zeros((rows, columns), dtype=bool)
    mask[16:48, 16:48] = True

    row, column = _representative_tile_origin(mask, 8)

    assert not mask[row : row + 8, column : column + 8].any()


def test_calibration_probe_scales_sample_size_with_budget():
    """A larger wall-time budget samples more pixels, not fewer."""
    rows, columns = 60, 60
    mask = np.zeros((rows, columns), dtype=bool)
    rays = _constant_rays(rows, columns)

    small_budget = _calibration_probe(
        _kernel(max_depth=0, threads=1),
        mask,
        rays,
        np.zeros(4),
        np.zeros(4),
        budget_seconds=0.01,
        max_sample_pixels=2000,
        rng=np.random.default_rng(0),
    )
    large_budget = _calibration_probe(
        _kernel(max_depth=0, threads=1),
        mask,
        rays,
        np.zeros(4),
        np.zeros(4),
        budget_seconds=0.05,
        max_sample_pixels=2000,
        rng=np.random.default_rng(0),
    )

    assert large_budget["sampled_pixels"] >= small_budget["sampled_pixels"]


def test_kernel_threads_sweep_scales_its_sample_to_a_whole_frame():
    """The sweep measures one sample tile, but the scheduler uses the
    result to decide how many frames can be in flight -- which needs the
    time to map a whole frame. Measured against the sample alone it came
    out low by the ratio of a frame to that tile, six times over on a
    real detector, which made every thread count look affordable.

    Timing is real everywhere except the GitHub macOS runner, which cannot
    hold a 20 ms sleep to anything like 20 ms (see ``_FAKE_SWEEP_CLOCK``).
    The property under test -- that the reported time scales by
    ``frame_pixels / sampled_pixels`` -- is arithmetic, so a faked clock
    tests it exactly rather than approximately.
    """
    import orgui.datautils.xrayutils.reconstruction as reconstruction_module

    fake_now = [0.0]

    class _FakeSweepKernel:
        def accumulate(self, *args, **kwargs):
            if _FAKE_SWEEP_CLOCK:
                fake_now[0] += 0.02
            else:
                time.sleep(0.02)
            return {}

    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0, threads=1)
    rows = columns = 8
    mask = np.zeros((rows, columns), dtype=bool)
    rays = _constant_rays(rows, columns)
    sampled_pixels = rows * columns

    def sweep(frame_pixels):
        return _kernel_threads_sweep(
            spec,
            grid,
            None,
            mask,
            rays,
            np.zeros(4),
            np.zeros(4),
            candidates=[1],
            frame_pixels=frame_pixels,
        )

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(
            reconstruction_module,
            "_kernel_for_grid",
            lambda *args, **kwargs: _FakeSweepKernel(),
        )
        if _FAKE_SWEEP_CLOCK:
            patch.setattr(
                reconstruction_module.time,
                "perf_counter",
                lambda: fake_now[0],
            )
        unscaled = sweep(None)
        per_frame = sweep(sampled_pixels * 100)

    # A frame a hundred times the sample must report a hundredfold time.
    # Real sleeps only hold that to within a quarter; a faked clock is exact.
    tolerance = 1e-9 if _FAKE_SWEEP_CLOCK else 0.25
    assert per_frame[1] == pytest.approx(100 * unscaled[1], rel=tolerance)


def test_kernel_threads_sweep_stops_early_on_plateau(monkeypatch):
    """Real native-kernel timing on a tiny test grid is noise-level (design
    doc Sec7), so this drives the sweep's own control flow (candidate
    loop, plateau early-stop, inherited values) with a fake, deterministic
    per-candidate kernel/timing instead of relying on real thread-scaling
    behavior.

    The clock is faked along with the kernel. Sleeping for the scripted
    durations and measuring them with the real ``perf_counter`` puts CI
    scheduling jitter inside the measured interval: the plateau compares a
    candidate against 0.9x the previous one, so an overshoot on the short
    candidate alone can invert the comparison. That was observed on a macOS
    runner, where candidate 8 was measured even though 4 had already
    plateaued.
    """
    import orgui.datautils.xrayutils.reconstruction as reconstruction_module

    call_log = []
    durations = {1: 0.15, 2: 0.03, 4: 0.05, 8: 0.001}
    fake_now = [0.0]
    monkeypatch.setattr(reconstruction_module.time, "perf_counter", lambda: fake_now[0])

    class _FakeSweepKernel:
        def __init__(self, threads):
            self._threads = threads

        def accumulate(self, *args, **kwargs):
            call_log.append(self._threads)
            fake_now[0] += durations[self._threads]
            return {}

    def fake_kernel_for_grid(
        spec, grid, ub_calculator, *, threads=None, memory_budget_bytes=None
    ):
        return _FakeSweepKernel(threads)

    monkeypatch.setattr(reconstruction_module, "_kernel_for_grid", fake_kernel_for_grid)

    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)
    mask = np.zeros((4, 4), dtype=bool)
    rays = _constant_rays(4, 4)

    results = _kernel_threads_sweep(
        spec,
        grid,
        _FakeUB(),
        mask,
        rays,
        np.zeros(4),
        np.zeros(4),
        candidates=[1, 2, 4, 8],
        tile_pixels=16,
        plateau_ratio=0.9,
    )

    assert set(results) == {1, 2, 4, 8}
    assert all(value > 0 for value in results.values())
    # A real 5x speedup 1 -> 2, then diminishing returns at 4 (slower than
    # 2) trips the plateau -- candidate 8 must not be measured at all, and
    # inherits candidate 4's measured time exactly.
    assert call_log == [1, 2, 4]
    assert results[1] > results[2]
    assert results[8] == results[4]


def test_footprint_split_conserves_weight_and_pixel_variance():
    rays = np.empty((2, 2, 3), dtype=np.float64)
    for column, x_value in enumerate((-0.2, 0.2)):
        rays[:, column, 0] = x_value
        rays[:, column, 1] = np.sqrt(1.0 - x_value**2)
        rays[:, column, 2] = 0.0
    result = _kernel(max_depth=2).accumulate(
        np.array([[10.0]]),
        np.array([[10.0]]),
        np.array([[False]]),
        rays,
        np.zeros(4),
        np.zeros(4),
    )

    assert result["weight"].size > 1
    assert np.sum(result["weight"]) == 1.0
    assert np.sum(result["weighted_intensity"]) == 10.0
    assert np.allclose(result["weighted_variance"], result["weight"] ** 2 * 10.0)
    assert np.all(result["contributors"] == 1)


class _FakeDetector:
    detector = type("Detector", (), {"shape": (1, 1)})()

    def __init__(self):
        self.ray_calculations = 0

    def primBeamPoints(self, rows, columns):
        self.ray_calculations += 1
        return np.zeros_like(rows), np.zeros_like(columns)


class _FakeUB:
    def getUB(self):
        return np.eye(3)

    def getU(self):
        return np.eye(3)

    def getK(self):
        return 1.0


class _FakeScan:
    shape = (1, 1)

    def __init__(self, values):
        self.values = values
        self.image_loads = 0

    def __len__(self):
        return len(self.values)

    def get_raw_img(self, index):
        self.image_loads += 1
        return h5_Image(np.array([[self.values[index]]], dtype=np.float64))

    def exposure_angle_bounds(self, config, fallback="stationary"):
        return np.zeros((len(self), 2, 4), dtype=np.float64)


def test_derive_grid_uses_native_corner_ray_interface():
    config = type(
        "Config",
        (),
        {
            "detector": _FakeDetector(),
            "ub_calculator": _FakeUB(),
        },
    )()
    scan = _FakeScan([1.0, 2.0])

    for frame in ("lab", "alpha", "omega", "chi", "phi", "crystal", "hkl"):
        grid = derive_grid(config, scan, frame=frame)
        assert grid.frame == frame
        assert np.all(np.isfinite(grid.minimum))
        assert np.all(np.isfinite(grid.maximum))
        assert np.all(np.asarray(grid.maximum) > np.asarray(grid.minimum))


def _router(
    boundaries, *, tmp_path, spec_digest="test-digest", budget=1024**3, resumed=()
):
    return _CheckpointRouter(
        boundaries,
        spec_digest=spec_digest,
        checkpoint_dir=tmp_path / "checkpoints",
        active_budget_bytes=budget,
        resumed=resumed,
    )


def test_map_frame_group_routes_batches_and_skips_resumed_checkpoints(tmp_path):
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)
    scan = _FakeScan([10.0, 20.0])

    def correction(payload, raw, frame, tile):
        return (
            raw.astype(np.float64),
            np.maximum(raw, 0).astype(np.float64),
            np.zeros(raw.shape, dtype=bool),
        )

    router = _router(
        {grid.grid_name: [(0, 1), (1, 2)]},
        tmp_path=tmp_path,
        resumed={(grid.grid_name, 1)},
    )
    detector_tiles = ((0, 1, 0, 1),)
    ray_arrays = _tile_ray_arrays(
        _FakeDetector(),
        detector_tiles,
        {(0, 1, 0, 1): _constant_rays(1, 1)},
    )
    kernels = _build_kernels(spec, _FakeUB())
    bounds = np.zeros((2, 2, 4), dtype=np.float64)

    for frame_index in range(2):
        payload = scan.get_raw_img(frame_index)
        _map_frame_group(
            spec,
            kernels,
            ray_arrays,
            detector_tiles,
            correction,
            [payload],
            [frame_index],
            bounds[frame_index : frame_index + 1, 0],
            bounds[frame_index : frame_index + 1, 1],
            router,
        )

    # Both frames were loaded (mapping is not resume-aware per frame; the
    # caller decides which frames are worth submitting), but only
    # checkpoint 0's file was written -- checkpoint 1 was pre-resumed.
    assert scan.image_loads == 2
    written_indices = {
        int(path.name.split("_")[0][len("ckpt") :]) for path in router.written
    }
    assert written_indices == {0}


def test_map_frame_group_corrects_once_before_tiling(tmp_path):
    class Detector:
        detector = type("PyfaiDetector", (), {"shape": (1, 2)})()

        def primBeamPoints(self, rows, columns):
            return np.zeros_like(rows), np.zeros_like(columns)

    calls = {"frame": 0, "tile": 0}

    def correction(payload, raw, frame, tile):
        calls["tile"] += 1
        raise AssertionError("tile correction must not be used")

    def correct_frame(payload, raw, frame):
        calls["frame"] += 1
        return (
            raw.astype(np.float64),
            np.maximum(raw, 0.0).astype(np.float64),
            np.zeros(raw.shape, dtype=bool),
        )

    correction.correct_frame = correct_frame
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)
    tiles = ((0, 1, 0, 1), (0, 1, 1, 2))
    router = _router({grid.grid_name: [(0, 1)]}, tmp_path=tmp_path)
    ray_arrays = _tile_ray_arrays(
        Detector(),
        tiles,
        {tiles[0]: _constant_rays(1, 1), tiles[1]: _constant_rays(1, 1)},
    )
    kernels = _build_kernels(spec, _FakeUB())
    payload = h5_Image(np.array([[10.0, 20.0]]))

    _map_frame_group(
        spec,
        kernels,
        ray_arrays,
        tiles,
        correction,
        [payload],
        [0],
        np.zeros((1, 4), dtype=np.float64),
        np.zeros((1, 4), dtype=np.float64),
        router,
    )

    assert calls == {"frame": 1, "tile": 0}
    assert len(router.written) == 1


class _RecordingRouter:
    """Collects routed batches instead of writing checkpoints."""

    def __init__(self):
        self.batches = []
        self.frames = []

    def route(self, grid_name, frame_index, batch, *, frames=1):
        self.batches.append(
            {name: np.asarray(values) for name, values in batch.items()}
        )
        self.frames.append((grid_name, int(frame_index), int(frames)))


def _totals(batches):
    """Per-voxel contributor counts and weight, summed over batches.

    Independent of how the frames were split into calls, which is exactly
    the property a frame group must not disturb.
    """
    contributors: dict[tuple[int, int], int] = {}
    weight: dict[tuple[int, int], float] = {}
    for batch in batches:
        for chunk, voxel, count, mass in zip(
            batch["chunk_id"],
            batch["local_voxel_id"],
            batch["contributors"],
            batch["weight"],
        ):
            key = (int(chunk), int(voxel))
            contributors[key] = contributors.get(key, 0) + int(count)
            weight[key] = weight.get(key, 0.0) + float(mass)
    return contributors, weight


def test_map_frame_group_reaches_the_same_voxels_as_single_frame_calls():
    """Grouping frames must not change the answer.

    It is deliberately not bit-for-bit -- several frames merge in the
    block map rather than in the tree accumulator, so sums associate
    differently -- so what is pinned here is what must not move: which
    voxels were reached, and how many samples reached each of them.
    """
    grid = _GridSpec(
        minimum=(-4.0, -4.0, -4.0),
        maximum=(4.0, 4.0, 4.0),
        step=(0.25, 0.25, 0.25),
        frame="lab",
        chunk_shape=(4, 4, 4),
    )
    frames = 4
    scan = _FakeScan([10.0, 20.0, 30.0, 40.0])
    detector_tiles = ((0, 1, 0, 1),)
    ray_arrays = _tile_ray_arrays(
        _FakeDetector(), detector_tiles, {(0, 1, 0, 1): _constant_rays(1, 1)}
    )
    # A rotation, so the frames land at different places in the volume and
    # a group has something to merge rather than nothing.
    bounds = np.zeros((frames, 2, 4), dtype=np.float64)
    bounds[:, :, 1] = np.deg2rad(0.1) * np.arange(frames)[:, None]

    def correction(payload, raw, frame, tile):
        return (
            raw.astype(np.float64),
            np.maximum(raw, 0.0).astype(np.float64),
            np.zeros(raw.shape, dtype=bool),
        )

    results = {}
    for frames_per_group in (1, 2, 4):
        spec = _ReconstructionSpec(
            grids=(grid,), max_depth=0, frames_per_group=frames_per_group
        )
        kernels = _build_kernels(spec, _FakeUB())
        router = _RecordingRouter()
        for origin in range(0, frames, frames_per_group):
            group = list(range(origin, origin + frames_per_group))
            _map_frame_group(
                spec,
                kernels,
                ray_arrays,
                detector_tiles,
                correction,
                [scan.get_raw_img(index) for index in group],
                group,
                np.ascontiguousarray(bounds[group, 0]),
                np.ascontiguousarray(bounds[group, 1]),
                router,
            )
        # One route() call per group, each declaring its own frame count,
        # so the checkpoint countdown still sees every frame exactly once.
        assert sum(entry[2] for entry in router.frames) == frames
        assert len(router.frames) == frames // frames_per_group
        results[frames_per_group] = (
            _totals(router.batches),
            sum(int(batch["chunk_id"].size) for batch in router.batches),
        )

    (reference_contributors, reference_weight), reference_rows = results[1]
    assert reference_contributors, "the fixture must reach at least one voxel"
    for frames_per_group in (2, 4):
        (contributors, weight), rows = results[frames_per_group]
        assert contributors == reference_contributors
        assert set(weight) == set(reference_weight)
        for key, value in weight.items():
            assert value == pytest.approx(reference_weight[key], rel=1e-12)
        # The point of grouping: the same samples leave the kernel as
        # fewer records, because frames sharing a voxel now merge inside
        # one call. Without this the test would still pass if grouping
        # silently degraded to mapping each frame on its own.
        assert rows < reference_rows

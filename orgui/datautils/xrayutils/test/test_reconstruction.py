"""Regression tests for native reciprocal-space reconstruction."""

from __future__ import annotations

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
    _calibration_probe,
    _files_per_job,
    _kernel_threads_sweep,
    _map_one_frame,
    _merge_sorted_batches,
    _reduce_batches,
    _tile_ray_arrays,
)


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
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
    np.testing.assert_array_equal(
        merged["weighted_intensity"], [10.0, 20.0, 34.0]
    )
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
        actual = kernel.coordinate(
            rays, sample_angles, sample_angles, 0, 0
        )
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


def test_kernel_threads_sweep_stops_early_on_plateau(monkeypatch):
    """Real native-kernel timing on a tiny test grid is noise-level (design
    doc Sec7), so this drives the sweep's own control flow (candidate
    loop, plateau early-stop, inherited values) with a fake, deterministic
    per-candidate kernel/timing instead of relying on real thread-scaling
    behavior."""
    import orgui.datautils.xrayutils.reconstruction as reconstruction_module

    call_log = []
    durations = {1: 0.15, 2: 0.03, 4: 0.05, 8: 0.001}

    class _FakeSweepKernel:
        def __init__(self, threads):
            self._threads = threads

        def accumulate(self, *args, **kwargs):
            call_log.append(self._threads)
            time.sleep(durations[self._threads])
            return {}

    def fake_kernel_for_grid(
        spec, grid, ub_calculator, *, threads=None, memory_budget_bytes=None
    ):
        return _FakeSweepKernel(threads)

    monkeypatch.setattr(
        reconstruction_module, "_kernel_for_grid", fake_kernel_for_grid
    )

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
    assert np.allclose(
        result["weighted_variance"], result["weight"] ** 2 * 10.0
    )
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


def test_map_one_frame_routes_batches_and_skips_resumed_checkpoints(tmp_path):
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
        _map_one_frame(
            spec,
            kernels,
            ray_arrays,
            detector_tiles,
            correction,
            payload,
            frame_index,
            bounds[frame_index, 0],
            bounds[frame_index, 1],
            router,
        )

    # Both frames were loaded (mapping is not resume-aware per frame; the
    # caller decides which frames are worth submitting), but only
    # checkpoint 0's file was written -- checkpoint 1 was pre-resumed.
    assert scan.image_loads == 2
    written_indices = {
        int(path.name.split("_")[0][len("ckpt"):]) for path in router.written
    }
    assert written_indices == {0}


def test_map_one_frame_corrects_once_before_tiling(tmp_path):
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

    _map_one_frame(
        spec,
        kernels,
        ray_arrays,
        tiles,
        correction,
        payload,
        0,
        np.zeros(4, dtype=np.float64),
        np.zeros(4, dtype=np.float64),
        router,
    )

    assert calls == {"frame": 1, "tile": 0}
    assert len(router.written) == 1

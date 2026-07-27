"""Regression tests for native reciprocal-space reconstruction."""

from __future__ import annotations

import json

import h5py
import numpy as np
import pytest

from orgui.backend.scans import h5_Image
from orgui.app.config_data import CorrectionState
from orgui.app.mask_config import repair_intensity_variance
from orgui.datautils.xrayutils import HKLVlieg
import orgui.datautils.xrayutils.reconstruction as reconstruction
from orgui.reconstruction_job import _correction_pipeline, derive_grid
from orgui.datautils.xrayutils.reconstruction import (
    _GridSpec,
    _ReconstructionSpec,
    _finalize_reconstruction,
    _map_frame_range,
    _map_partition,
    _merge_sorted_batches,
    _read_parquet,
    _reduce_batches,
    _reduce_partition,
    _write_manifest,
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


def test_map_partition_accepts_local_image_and_ray_caches(tmp_path):
    pytest.importorskip("pyarrow")
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)
    scan = _FakeScan([10.0])
    detector = _FakeDetector()
    payload = h5_Image(np.array([[10.0]], dtype=np.float64))

    _map_partition(
        spec,
        scan,
        detector,
        _FakeUB(),
        (0, 1),
        (0, 1, 0, 1),
        np.zeros((1, 2, 4), dtype=np.float64),
        tmp_path / "map",
        correction_pipeline=lambda payload, raw, frame, tile: (
            raw.astype(np.float64),
            np.maximum(raw, 0).astype(np.float64),
            np.zeros(raw.shape, dtype=bool),
        ),
        job_digest="cache-test",
        image_payloads={0: payload},
        corner_rays=_constant_rays(1, 1),
    )

    assert scan.image_loads == 0
    assert detector.ray_calculations == 0


def test_map_frame_range_corrects_once_before_tiling(tmp_path):
    pytest.importorskip("pyarrow")

    class Detector:
        detector = type("PyfaiDetector", (), {"shape": (1, 2)})()

        def primBeamPoints(self, rows, columns):
            return np.zeros_like(rows), np.zeros_like(columns)

    class Scan:
        def __len__(self):
            return 1

        def get_raw_img(self, index):
            return h5_Image(np.array([[10.0, 20.0]]))

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

    _map_frame_range(
        spec,
        Scan(),
        Detector(),
        _FakeUB(),
        (0, 1),
        tiles,
        np.zeros((1, 2, 4), dtype=np.float64),
        tmp_path / "map",
        correction_pipeline=correction,
        job_digest="full-frame-correction",
        corner_rays={
            tiles[0]: _constant_rays(1, 1),
            tiles[1]: _constant_rays(1, 1),
        },
    )

    assert calls == {"frame": 1, "tile": 0}


def test_streaming_accumulator_flushes_by_retained_bytes(tmp_path):
    """A larger byte ceiling combines images before Parquet serialization."""
    pytest.importorskip("pyarrow")
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)
    bounds = np.zeros((3, 2, 4), dtype=np.float64)
    tile = (0, 1, 0, 1)

    def correction(payload, raw, frame, detector_tile):
        return (
            raw.astype(np.float64),
            np.maximum(raw, 0).astype(np.float64),
            np.zeros(raw.shape, dtype=bool),
        )

    small_scan = _FakeScan([10.0, 20.0, 30.0])
    image_updates = []
    small = _map_frame_range(
        spec,
        small_scan,
        _FakeDetector(),
        _FakeUB(),
        (0, 3),
        (tile,),
        bounds,
        tmp_path / "small",
        correction_pipeline=correction,
        job_digest="small-accumulator",
        corner_rays={tile: _constant_rays(1, 1)},
        accumulation_budget_bytes=1,
        image_progress=lambda frame, retained, segments: image_updates.append(
            (frame, retained, segments)
        ),
    )
    large_scan = _FakeScan([10.0, 20.0, 30.0])
    large = _map_frame_range(
        spec,
        large_scan,
        _FakeDetector(),
        _FakeUB(),
        (0, 3),
        (tile,),
        bounds,
        tmp_path / "large",
        correction_pipeline=correction,
        job_digest="large-accumulator",
        corner_rays={tile: _constant_rays(1, 1)},
        accumulation_budget_bytes=1024**2,
    )

    assert small_scan.image_loads == 3
    assert large_scan.image_loads == 3
    assert [update[0] for update in image_updates] == [0, 1, 2]
    assert small.metadata["accumulation_segments"] == 3
    assert large.metadata["accumulation_segments"] == 1
    assert len(small.partitions) == 3
    assert len(large.partitions) == 1
    small_batch = _reduce_batches(
        _read_parquet(partition.uri) for partition in small.partitions
    )
    large_batch = _reduce_batches(
        _read_parquet(partition.uri) for partition in large.partitions
    )
    for name in small_batch:
        np.testing.assert_array_equal(small_batch[name], large_batch[name])


def test_parquet_reduce_and_hdf5_finalize(tmp_path, monkeypatch):
    pytest.importorskip("pyarrow")
    checksum_calls = []
    original_checksum = reconstruction._uri_checksum_and_size

    def counted_checksum(uri):
        checksum_calls.append(str(uri))
        return original_checksum(uri)

    monkeypatch.setattr(
        reconstruction, "_uri_checksum_and_size", counted_checksum
    )
    verification_cache = set()
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    spec = _ReconstructionSpec(
        grids=(grid,),
        max_depth=1,
        threads=2,
        compression="gzip",
    )
    scan = _FakeScan([10.0, 20.0])
    bounds = np.zeros((1, 2, 4), dtype=np.float64)
    map_manifests = []
    for frame in range(2):
        manifest = _map_partition(
            spec,
            scan,
            _FakeDetector(),
            _FakeUB(),
            (frame, frame + 1),
            (0, 1, 0, 1),
            bounds,
            tmp_path / "map",
            correction_pipeline=lambda payload, raw, frame, tile: (
                raw.astype(np.float64),
                np.maximum(raw, 0).astype(np.float64),
                np.zeros(raw.shape, dtype=bool),
            ),
            job_digest="test-job",
            verification_cache=verification_cache,
        )
        path = tmp_path / f"map-{frame}.json"
        _write_manifest(manifest, path)
        map_manifests.append(path)

    reduced = _reduce_partition(
        map_manifests,
        tmp_path / "reduced",
        verification_cache=verification_cache,
    )
    reduced_path = tmp_path / "reduced.json"
    _write_manifest(reduced, reduced_path)
    output = tmp_path / "result.h5"
    result = _finalize_reconstruction(
        [reduced_path],
        output,
        verification_cache=verification_cache,
    )

    assert result["chunks_written"] == 1
    assert len(checksum_calls) == 3
    assert all(
        partition.size_bytes
        for path in map_manifests
        for partition in reconstruction._read_manifest(path).partitions
    )
    assert all(chunk.size_bytes for chunk in reduced.chunks)
    with h5py.File(output, "r") as h5file:
        group = h5file["entry/reconstruction/results/q_lab"]
        assert group.attrs["NX_class"] == "NXdata"
        assert group["intensity"][1, 1, 1] == 15.0
        assert group["variance"][1, 1, 1] == 7.5
        assert group["weight"][1, 1, 1] == 2.0
        assert group["contributors"][1, 1, 1] == 2
        assert np.isnan(group["intensity"][0, 0, 0])
        saved = json.loads(h5file["entry/reconstruction/configuration_json"][()])
        assert saved["grids"][0]["frame"] == "lab"


def test_reduce_shards_large_chunk_and_reuses_checkpoint(tmp_path, monkeypatch):
    """Large chunks reduce by bounded local ranges and resume by shard."""
    pytest.importorskip("pyarrow")
    grid = _GridSpec(
        minimum=(0.0, 0.0, 0.0),
        maximum=(1.0, 1.0, 300000.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(1, 1, 300000),
    )
    spec = _ReconstructionSpec(
        grids=(grid,),
        max_depth=0,
        threads=1,
        compression="gzip",
    )
    manifest_paths = []
    for task_index, intensities in enumerate(
        ([2.0, 4.0, 6.0], [3.0, 5.0, 7.0])
    ):
        batch = {
            "chunk_id": np.zeros(3, dtype=np.uint64),
            "local_voxel_id": np.array(
                [1, 260000, 299999], dtype=np.uint64
            ),
            "weighted_intensity": np.asarray(
                intensities, dtype=np.float64
            ),
            "weighted_variance": np.ones(3, dtype=np.float64),
            "weight": np.ones(3, dtype=np.float64),
            "contributors": np.ones(3, dtype=np.uint64),
        }
        partition_path = tmp_path / f"map-{task_index}.parquet"
        reconstruction._write_parquet(batch, partition_path, {})
        checksum, size_bytes = reconstruction._uri_checksum_and_size(
            partition_path
        )
        manifest = reconstruction._TaskManifest(
            kind="map",
            task_id=f"map-{task_index}",
            spec_hash=spec.digest,
            status="complete",
            spec=spec.to_dict(),
            partitions=[
                reconstruction._PartitionFile(
                    grid_name=grid.grid_name,
                    bucket=0,
                    uri=str(partition_path),
                    rows=3,
                    checksum=checksum,
                    size_bytes=size_bytes,
                )
            ],
            metadata={"base_context_sha256": "same-context"},
        )
        manifest_path = tmp_path / f"map-{task_index}.json"
        _write_manifest(manifest, manifest_path)
        manifest_paths.append(manifest_path)

    progress = []
    verification_cache = set()
    reduced = _reduce_partition(
        manifest_paths,
        tmp_path / "reduced",
        verification_cache=verification_cache,
        memory_budget_bytes=1,
        progress=lambda current, total, message: progress.append(
            (current, total, message)
        ),
        checkpoint_root=tmp_path / "manifests",
    )

    assert reduced.status == "complete"
    assert len(reduced.chunks) == 2
    assert {chunk.shard_start for chunk in reduced.chunks} == {0, 262144}
    assert progress[-1][:2] == (2, 2)

    def unexpected_read(*args, **kwargs):
        raise AssertionError("completed shards must not read source Parquet")

    monkeypatch.setattr(
        reconstruction._ParquetRangeReader, "read", unexpected_read
    )
    resumed = _reduce_partition(
        manifest_paths,
        tmp_path / "reduced",
        verification_cache=verification_cache,
        memory_budget_bytes=1,
        checkpoint_root=tmp_path / "manifests",
    )
    assert resumed.status == "complete"
    assert len(resumed.chunks) == 2

    reduced_path = tmp_path / "reduced.json"
    _write_manifest(resumed, reduced_path)
    output = tmp_path / "sharded-result.h5"
    _finalize_reconstruction(
        [reduced_path],
        output,
        verification_cache=verification_cache,
    )
    with h5py.File(output, "r") as h5file:
        group = h5file["entry/reconstruction/results/q_lab"]
        np.testing.assert_array_equal(
            group["intensity"][0, 0, [1, 260000, 299999]],
            [2.5, 4.5, 6.5],
        )
        np.testing.assert_array_equal(
            group["contributors"][0, 0, [1, 260000, 299999]],
            [2, 2, 2],
        )

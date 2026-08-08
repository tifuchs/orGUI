"""Tests for the HDF5 checkpoint format and tree-merge accumulator (design
doc Sec8/Sec9): :func:`_write_checkpoint`/:func:`_read_checkpoint`/
:func:`_verify_checkpoint` and :class:`_CheckpointAccumulator`.
"""

from __future__ import annotations

import h5py
import numpy as np
import pytest

from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointAccumulator,
    _merge_sorted_batches,
    _read_checkpoint,
    _reduce_batches,
    _tree_finalize,
    _tree_insert,
    _verify_checkpoint,
    _write_checkpoint,
)


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


def _batch(chunk_id, local, intensity, contributors=None):
    size = len(local)
    return {
        "chunk_id": np.full(size, chunk_id, dtype=np.uint32),
        "local_voxel_id": np.asarray(local, dtype=np.uint32),
        "weighted_intensity": np.asarray(intensity, dtype=np.float64),
        "weighted_variance": np.asarray(intensity, dtype=np.float64),
        "weight": np.ones(size, dtype=np.float64),
        "contributors": (
            np.ones(size, dtype=np.uint32)
            if contributors is None
            else np.asarray(contributors, dtype=np.uint32)
        ),
    }


def _kernel(*, threads=1, max_depth=0, ub=None, memory_budget_bytes=256 * 1024 * 1024):
    if ub is None:
        ub = np.eye(3)
    return native.ReconstructionKernel(
        np.array([-20.0, -20.0, -20.0]),
        np.array([0.01, 0.01, 0.01]),
        np.array([4000, 4000, 4000], dtype=np.int64),
        np.array([64, 64, 64], dtype=np.int64),
        "lab",
        1.0,
        np.linalg.inv(ub),
        np.eye(3),
        max_depth,
        threads,
        16,
        memory_budget_bytes,
    )


def _constant_rays(rows, columns):
    rays = np.zeros((rows + 1, columns + 1, 3), dtype=np.float64)
    rays[..., 1] = 1.0
    return rays


# ---------------------------------------------------------------------------
# Tree-insert/finalize (extracted from _reduce_batches -- Sec9 "Rejected:
# naive incremental merge" / "Adopted" split).
# ---------------------------------------------------------------------------


def test_tree_insert_and_finalize_match_reduce_batches_directly():
    """The extracted helpers, driven by hand, match _reduce_batches on the
    same input -- the refactor this test guards is behavior-preserving by
    construction, not just by re-running the existing suite."""
    batches = [
        _batch(0, [1, 3], [10.0, 30.0], [1, 2]),
        _batch(0, [2, 3], [20.0, 4.0], [3, 4]),
        _batch(1, [0], [5.0]),
    ]

    expected = _reduce_batches(batches)

    levels: list = []
    for batch in batches:
        _tree_insert(levels, batch)
    actual = _tree_finalize(levels)

    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)


# ---------------------------------------------------------------------------
# _merge_sorted_batches (native, single-pass): output must be exactly
# trimmed to its true size, not left with unused worst-case-sized capacity.
# ---------------------------------------------------------------------------


def test_merge_sorted_batches_output_exact_size_duplicate_free():
    left = _batch(0, [1, 3, 5], [1.0, 2.0, 3.0])
    right = _batch(0, [2, 4, 6], [10.0, 20.0, 30.0])
    result = _merge_sorted_batches(left, right)
    for name in result:
        assert result[name].shape == (6,), name
    np.testing.assert_array_equal(result["local_voxel_id"], [1, 2, 3, 4, 5, 6])


def test_merge_sorted_batches_output_exact_size_all_duplicates():
    left = _batch(0, [1, 2, 3], [1.0, 2.0, 3.0], contributors=[1, 1, 1])
    right = _batch(0, [1, 2, 3], [10.0, 20.0, 30.0], contributors=[2, 2, 2])
    result = _merge_sorted_batches(left, right)
    for name in result:
        assert result[name].shape == (3,), name
    np.testing.assert_array_equal(result["local_voxel_id"], [1, 2, 3])
    np.testing.assert_array_equal(
        result["weighted_intensity"], [11.0, 22.0, 33.0]
    )
    np.testing.assert_array_equal(result["contributors"], [3, 3, 3])


def test_merge_sorted_batches_output_exact_size_partial_overlap():
    left = _batch(0, [1, 2, 4], [1.0, 2.0, 4.0])
    right = _batch(0, [2, 3, 4], [20.0, 3.0, 40.0])
    result = _merge_sorted_batches(left, right)
    for name in result:
        assert result[name].shape == (4,), name
    np.testing.assert_array_equal(result["local_voxel_id"], [1, 2, 3, 4])
    np.testing.assert_array_equal(
        result["weighted_intensity"], [1.0, 22.0, 3.0, 44.0]
    )


# ---------------------------------------------------------------------------
# _CheckpointAccumulator
# ---------------------------------------------------------------------------


def test_checkpoint_accumulator_matches_static_tree_reduce():
    """Feeding the same batches through the accumulator one at a time,
    versus _reduce_batches over the whole pre-collected list, produces
    identical output -- the accumulator is a thin incremental wrapper
    around the same primitive, not a different algorithm."""
    rng = np.random.default_rng(0)
    batches = [
        _batch(
            int(rng.integers(0, 4)),
            rng.integers(0, 50, size=20),
            rng.uniform(0, 10, size=20),
        )
        for _ in range(12)
    ]

    expected = _reduce_batches(batches)

    accumulator = _CheckpointAccumulator()
    for batch in batches:
        accumulator.insert(batch)
    actual = accumulator.finalize()

    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)


def test_checkpoint_accumulator_should_flush_tracks_byte_total():
    accumulator = _CheckpointAccumulator()
    batch = _batch(0, list(range(100)), np.ones(100))
    batch_bytes = sum(np.asarray(values).nbytes for values in batch.values())

    assert not accumulator.should_flush(batch_bytes * 2)
    accumulator.insert(batch)
    assert not accumulator.should_flush(batch_bytes * 2)
    accumulator.insert(batch)
    assert accumulator.should_flush(batch_bytes * 2)


def test_checkpoint_accumulator_finalize_resets_state():
    accumulator = _CheckpointAccumulator()
    accumulator.insert(_batch(0, [1, 2], [1.0, 2.0]))
    first = accumulator.finalize()
    assert first["chunk_id"].size == 2

    second = accumulator.finalize()
    assert second["chunk_id"].size == 0
    assert not accumulator.should_flush(1)


# ---------------------------------------------------------------------------
# HDF5 checkpoint I/O
# ---------------------------------------------------------------------------


def test_checkpoint_write_read_roundtrip_columnar_chunked(tmp_path):
    batch = _batch(0, list(range(200)), np.linspace(0, 1, 200))
    path = tmp_path / "checkpoint_000.h5"

    digest = _write_checkpoint(
        batch, path, spec_digest="deadbeef", metadata={"checkpoint_index": 0}
    )

    with h5py.File(path, "r") as h5file:
        assert h5file.attrs["spec_sha256"] == "deadbeef"
        assert h5file.attrs["xxh3_128"] == digest
        assert h5file.attrs["rows"] == 200
        assert h5file.attrs["checkpoint_index"] == 0
        for name in batch:
            dataset = h5file[name]
            assert dataset.chunks == (200,)
            # h5py's .compression convenience property does not recognize
            # dynamically-loaded third-party filters (hdf5plugin's
            # Bitshuffle-LZ4 included) and reads None even when genuinely
            # applied -- check the low-level filter list instead.
            assert dataset.id.get_create_plist().get_nfilters() >= 1

    roundtrip = _read_checkpoint(path)
    for name in batch:
        np.testing.assert_array_equal(roundtrip[name], batch[name], err_msg=name)
        assert roundtrip[name].dtype == batch[name].dtype, name
    # The narrow key/count columns must round-trip as uint32 (40-byte row
    # format), not the pre-narrowing uint64 -- a silent widening here would
    # double the on-disk/in-memory footprint _CHECKPOINT_BYTES_PER_ROW
    # assumes.
    assert roundtrip["chunk_id"].dtype == np.uint32
    assert roundtrip["local_voxel_id"].dtype == np.uint32
    assert roundtrip["contributors"].dtype == np.uint32


def test_checkpoint_write_chunks_are_capped_at_65536_rows(tmp_path):
    rows = 200_000
    batch = _batch(0, np.arange(rows) % 4096, np.ones(rows))
    path = tmp_path / "checkpoint_big.h5"

    _write_checkpoint(batch, path, spec_digest="x")

    with h5py.File(path, "r") as h5file:
        assert h5file["chunk_id"].chunks == (65536,)


def test_checkpoint_write_is_atomic_tmp_then_replace(tmp_path):
    batch = _batch(0, [1, 2, 3], [1.0, 2.0, 3.0])
    path = tmp_path / "checkpoint_atomic.h5"
    temporary = path.with_name(path.name + ".tmp")

    assert not path.exists()
    _write_checkpoint(batch, path, spec_digest="x")

    assert path.exists()
    assert not temporary.exists()


def test_checkpoint_digest_attrs_detect_spec_change(tmp_path):
    batch = _batch(0, [1, 2, 3], [1.0, 2.0, 3.0])
    path = tmp_path / "checkpoint_digest.h5"
    _write_checkpoint(batch, path, spec_digest="spec-a")

    assert _verify_checkpoint(path, spec_digest="spec-a") is True
    assert _verify_checkpoint(path, spec_digest="spec-b") is False
    assert _verify_checkpoint(tmp_path / "missing.h5", spec_digest="spec-a") is False


def test_checkpoint_full_verify_detects_corrupted_data(tmp_path):
    batch = _batch(0, [1, 2, 3], [1.0, 2.0, 3.0])
    path = tmp_path / "checkpoint_corrupt.h5"
    _write_checkpoint(batch, path, spec_digest="spec-a")

    assert _verify_checkpoint(path, spec_digest="spec-a", full=True) is True

    with h5py.File(path, "r+") as h5file:
        h5file["weighted_intensity"][0] = 999.0

    assert _verify_checkpoint(path, spec_digest="spec-a", full=False) is True
    assert _verify_checkpoint(path, spec_digest="spec-a", full=True) is False


# ---------------------------------------------------------------------------
# Validation strategy (design doc Sec9): equivalence against the current
# reduce path on synthetic data, since neither path has an external ground
# truth. This is the only correctness evidence either path has -- see the
# design doc's "Validation strategy" subsection.
# ---------------------------------------------------------------------------


def test_checkpoint_accumulator_matches_current_reduce_on_synthetic_scan():
    """Route a small deterministic synthetic 'scan' (a handful of frames
    with a fixed RNG seed, no real detector/crystal data -- mapping
    correctness depends only on geometry and the mask, not pixel values,
    per Sec6) through both the current per-flush _reduce_batches path and
    the new per-checkpoint tree-merge accumulator, and assert equivalent
    output. This is the "phys" claim of Phase 2: it establishes the new
    path preserves the old path's behavior, not that either is correct
    against an external reference, since none exists for either."""
    rows, columns = 12, 10
    rays = _constant_rays(rows, columns)
    kernel = _kernel(max_depth=2, threads=1)

    rng = np.random.default_rng(42)
    frame_batches = []
    for frame_index in range(9):
        intensity = rng.uniform(0.0, 10.0, size=(rows, columns))
        variance = np.maximum(intensity, 1.0)
        mask = rng.uniform(size=(rows, columns)) < 0.05
        # A small angular sweep per frame, like adjacent rocking-scan steps,
        # so different frames' footprints genuinely overlap and collide in
        # the same voxels (exercising the merge/dedup path, not just
        # disjoint inserts).
        omega = np.deg2rad(frame_index * 0.1)
        angles_start = np.array([0.0, omega, 0.0, 0.0])
        angles_end = np.array([0.0, omega + np.deg2rad(0.1), 0.0, 0.0])
        result = kernel.accumulate(
            intensity, variance, mask, rays, angles_start, angles_end
        )
        frame_batches.append(result)

    expected = _reduce_batches(frame_batches)

    accumulator = _CheckpointAccumulator()
    for batch in frame_batches:
        accumulator.insert(batch)
    actual = accumulator.finalize()

    assert expected["chunk_id"].size > 0, "synthetic scan produced no records at all"
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)

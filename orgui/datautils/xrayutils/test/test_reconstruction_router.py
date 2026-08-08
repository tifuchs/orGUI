"""Tests for the checkpoint router, resume discovery, and checkpoint-based
finalize (design doc Sec9/Sec10/Sec11): :class:`_CheckpointRouter`,
:func:`_discover_checkpoint_state`, :class:`_CheckpointRangeReader`, and
:func:`_finalize_reconstruction`.
"""

from __future__ import annotations

import threading
import time
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np
import pytest

from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointAccumulator,
    _CheckpointRangeReader,
    _CheckpointRouter,
    _GridSpec,
    _ReconstructionSpec,
    _discover_checkpoint_state,
    _finalize_reconstruction,
    _read_checkpoint,
    _reduce_batches,
    _write_checkpoint,
)


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


def _batch(chunk_id, local, intensity, contributors=None):
    size = len(local)
    return {
        "chunk_id": np.full(size, chunk_id, dtype=np.uint64),
        "local_voxel_id": np.asarray(local, dtype=np.uint64),
        "weighted_intensity": np.asarray(intensity, dtype=np.float64),
        "weighted_variance": np.asarray(intensity, dtype=np.float64),
        "weight": np.ones(size, dtype=np.float64),
        "contributors": (
            np.ones(size, dtype=np.uint64)
            if contributors is None
            else np.asarray(contributors, dtype=np.uint64)
        ),
    }


# ---------------------------------------------------------------------------
# _CheckpointRouter: frame -> checkpoint index lookup
# ---------------------------------------------------------------------------


def test_checkpoint_index_for_frame_matches_boundaries(tmp_path):
    router = _CheckpointRouter(
        {"hkl": [(0, 3), (3, 10), (10, 12)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=10**9,
    )
    for frame_index, expected in ((0, 0), (2, 0), (3, 1), (9, 1), (10, 2), (11, 2)):
        assert router.checkpoint_index_for_frame("hkl", frame_index) == expected


def test_checkpoint_index_for_frame_rejects_out_of_range(tmp_path):
    router = _CheckpointRouter(
        {"hkl": [(2, 5)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=10**9,
    )
    with pytest.raises(ValueError, match="precedes"):
        router.checkpoint_index_for_frame("hkl", 0)
    with pytest.raises(ValueError, match="outside"):
        router.checkpoint_index_for_frame("hkl", 5)


# ---------------------------------------------------------------------------
# _CheckpointRouter: routing, resume skip, safety-valve split
# ---------------------------------------------------------------------------


def test_router_writes_one_file_and_matches_reduce_batches(tmp_path):
    """Routing a set of per-frame batches produces one checkpoint file whose
    content matches directly reducing the same batches -- the router is a
    thin incremental wrapper, not a different merge algorithm."""
    rng = np.random.default_rng(1)
    batches = [
        _batch(0, rng.integers(0, 50, size=10), rng.uniform(0, 10, size=10))
        for _ in range(6)
    ]
    router = _CheckpointRouter(
        {"hkl": [(0, 6)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=10**9,
    )
    for frame_index, batch in enumerate(batches):
        router.route("hkl", frame_index, batch)

    assert len(router.written) == 1
    expected = _reduce_batches(batches)
    actual = _read_checkpoint(router.written[0])
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)


def test_router_skips_already_resumed_checkpoints(tmp_path):
    router = _CheckpointRouter(
        {"hkl": [(0, 2), (2, 4)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=10**9,
        resumed={("hkl", 1)},
    )
    for frame_index in range(4):
        router.route("hkl", frame_index, _batch(0, [frame_index], [1.0]))

    # Only checkpoint 0 (not pre-resumed) was ever written.
    assert len(router.written) == 1
    with h5py.File(router.written[0], "r") as h5file:
        assert h5file.attrs["checkpoint_index"] == 0


def test_router_safety_valve_splits_into_multiple_parts_summing_frames(tmp_path):
    """A tiny active_budget_bytes forces a split mid-checkpoint; the parts'
    frames_covered must still sum to the full planned frame count."""
    router = _CheckpointRouter(
        {"hkl": [(0, 5)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=1,  # any non-empty batch exceeds this
    )
    for frame_index in range(5):
        router.route(
            "hkl", frame_index, _batch(0, [frame_index], [float(frame_index)])
        )

    assert len(router.written) > 1
    total_covered = 0
    for path in router.written:
        with h5py.File(path, "r") as h5file:
            assert h5file.attrs["checkpoint_index"] == 0
            total_covered += int(h5file.attrs["frames_covered"])
    assert total_covered == 5

    boundaries = {"hkl": [(0, 5)]}
    resumed, files = _discover_checkpoint_state(tmp_path, boundaries, "d")
    assert resumed == {("hkl", 0)}
    assert sorted(files["hkl"]) == sorted(router.written)


# ---------------------------------------------------------------------------
# _CheckpointRouter: concurrent route() calls (the lock-scoping fix that
# moves accumulator.insert() outside the router's lock)
# ---------------------------------------------------------------------------


def test_router_concurrent_routes_to_same_key_lose_no_frames(tmp_path, monkeypatch):
    """route() no longer holds the router lock across insert(), so several
    frames destined for the same checkpoint can have their insert() calls
    genuinely overlap. Delay the *first* frame's insert() so it is still
    in flight when the *last* frame's route() call decrements _remaining
    to zero -- the exact window the in-flight counter exists to guard --
    and confirm the flush still waits for it, losing nothing."""
    rng = np.random.default_rng(2)
    n_frames = 8
    batches = [
        _batch(0, rng.integers(0, 50, size=10), rng.uniform(0, 10, size=10))
        for _ in range(n_frames)
    ]
    expected = _reduce_batches(batches)

    real_insert = _CheckpointAccumulator.insert
    delay_next = threading.Event()
    delay_next.set()  # delay exactly the first call that reaches insert()
    lock = threading.Lock()

    def delayed_insert(self, batch):
        should_delay = False
        with lock:
            if delay_next.is_set():
                delay_next.clear()
                should_delay = True
        real_insert(self, batch)
        if should_delay:
            # Released the accumulator's own lock already (real_insert
            # returned); hold up *returning to route()* so this frame's
            # in-flight count stays nonzero while other frames' route()
            # calls race ahead and reach their own flush check first.
            time.sleep(0.3)

    monkeypatch.setattr(_CheckpointAccumulator, "insert", delayed_insert)

    router = _CheckpointRouter(
        {"hkl": [(0, n_frames)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=10**9,
    )
    with ThreadPoolExecutor(max_workers=n_frames) as pool:
        list(
            pool.map(
                lambda item: router.route("hkl", item[0], item[1]),
                enumerate(batches),
            )
        )

    assert len(router.written) == 1
    actual = _read_checkpoint(router.written[0])
    for name in expected:
        np.testing.assert_array_equal(actual[name], expected[name], err_msg=name)
    with h5py.File(router.written[0], "r") as h5file:
        assert h5file.attrs["frames_covered"] == n_frames


def test_router_concurrent_routes_to_different_keys_do_not_block(tmp_path, monkeypatch):
    """A slow insert() cascade for one checkpoint must not stall route()
    calls for a different, unrelated checkpoint -- the reason the merge
    moved outside the router's single job-wide lock in the first place."""
    slow_key_seen = threading.Event()

    real_insert = _CheckpointAccumulator.insert

    def slow_insert(self, batch):
        real_insert(self, batch)
        if not slow_key_seen.is_set():
            slow_key_seen.set()
            time.sleep(0.5)

    monkeypatch.setattr(_CheckpointAccumulator, "insert", slow_insert)

    router = _CheckpointRouter(
        {"hkl": [(0, 1), (1, 2)]},
        spec_digest="d",
        checkpoint_dir=tmp_path,
        active_budget_bytes=10**9,
    )

    def route_slow():
        router.route("hkl", 0, _batch(0, [1], [1.0]))

    def route_fast():
        slow_key_seen.wait(timeout=5)  # ensure the slow route() has started first
        started = time.perf_counter()
        router.route("hkl", 1, _batch(1, [2], [2.0]))
        return time.perf_counter() - started

    with ThreadPoolExecutor(max_workers=2) as pool:
        slow_future = pool.submit(route_slow)
        fast_future = pool.submit(route_fast)
        fast_elapsed = fast_future.result(timeout=5)
        slow_future.result(timeout=5)

    assert fast_elapsed < 0.4  # well under the slow route's 0.5s hold
    assert len(router.written) == 2


# ---------------------------------------------------------------------------
# _discover_checkpoint_state
# ---------------------------------------------------------------------------


def test_discover_checkpoint_state_detects_partial_and_stale(tmp_path):
    boundaries = {"hkl": [(0, 3), (3, 6)]}
    complete_path = tmp_path / "hkl" / "ckpt0000_p0000.h5"
    _write_checkpoint(
        _batch(0, [1], [1.0]),
        complete_path,
        spec_digest="d",
        metadata={"checkpoint_index": 0, "frames_covered": 3},
    )
    partial_path = tmp_path / "hkl" / "ckpt0001_p0000.h5"
    _write_checkpoint(
        _batch(0, [2], [1.0]),
        partial_path,
        spec_digest="d",
        metadata={"checkpoint_index": 1, "frames_covered": 2},
    )

    resumed, files = _discover_checkpoint_state(tmp_path, boundaries, "d")

    assert resumed == {("hkl", 0)}
    assert files["hkl"] == [complete_path]
    assert partial_path.exists()  # not cleaned up without cleanup_stale=True


def test_discover_checkpoint_state_cleans_up_stale_parts(tmp_path):
    boundaries = {"hkl": [(0, 3)]}
    stale_path = tmp_path / "hkl" / "ckpt0000_p0000.h5"
    _write_checkpoint(
        _batch(0, [1], [1.0]),
        stale_path,
        spec_digest="d",
        metadata={"checkpoint_index": 0, "frames_covered": 1},
    )

    resumed, _files = _discover_checkpoint_state(
        tmp_path, boundaries, "d", cleanup_stale=True
    )

    assert resumed == set()
    assert not stale_path.exists()


def test_discover_checkpoint_state_rejects_digest_mismatch(tmp_path):
    boundaries = {"hkl": [(0, 1)]}
    path = tmp_path / "hkl" / "ckpt0000_p0000.h5"
    _write_checkpoint(
        _batch(0, [1], [1.0]),
        path,
        spec_digest="old-digest",
        metadata={"checkpoint_index": 0, "frames_covered": 1},
    )

    resumed, files = _discover_checkpoint_state(tmp_path, boundaries, "new-digest")

    assert resumed == set()
    assert files["hkl"] == []


# ---------------------------------------------------------------------------
# _CheckpointRangeReader
# ---------------------------------------------------------------------------


def test_checkpoint_range_reader_matches_direct_slicing(tmp_path):
    batch = {
        "chunk_id": np.array([0, 0, 0, 1, 1, 2], dtype=np.uint64),
        "local_voxel_id": np.array([1, 5, 9, 2, 7, 0], dtype=np.uint64),
        "weighted_intensity": np.arange(6, dtype=np.float64),
        "weighted_variance": np.arange(6, dtype=np.float64),
        "weight": np.ones(6, dtype=np.float64),
        "contributors": np.ones(6, dtype=np.uint64),
    }
    path = tmp_path / "reader.h5"
    _write_checkpoint(batch, path, spec_digest="d")

    reader = _CheckpointRangeReader(path, batch_rows=2)
    try:
        chunk0 = reader.read(0, 0, 10)
        np.testing.assert_array_equal(chunk0["local_voxel_id"], [1, 5, 9])
        chunk1 = reader.read(1, 0, 10)
        np.testing.assert_array_equal(chunk1["local_voxel_id"], [2, 7])
        chunk2 = reader.read(2, 0, 10)
        np.testing.assert_array_equal(chunk2["local_voxel_id"], [0])
    finally:
        reader.close()


def test_checkpoint_range_reader_rejects_out_of_order_reads(tmp_path):
    batch = _batch(1, [5], [1.0])
    path = tmp_path / "reader.h5"
    _write_checkpoint(batch, path, spec_digest="d")

    reader = _CheckpointRangeReader(path)
    try:
        reader.read(1, 0, 10)
        with pytest.raises(ValueError, match="sorted order"):
            reader.read(0, 0, 10)
    finally:
        reader.close()


# ---------------------------------------------------------------------------
# _finalize_reconstruction: reads directly from checkpoint files
# ---------------------------------------------------------------------------


def _grid():
    return _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )


def test_finalize_reads_single_checkpoint_per_grid(tmp_path):
    grid = _grid()
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0, compression="gzip")
    # The whole (2,2,2) grid is one chunk (chunk_shape == grid shape);
    # global voxel (1,1,1) is local_voxel_id 1*2*2 + 1*2 + 1 == 7.
    batch = _batch(0, [7], [10.0], contributors=[1])
    batch["weight"] = np.array([2.0])
    batch["weighted_variance"] = np.array([4.0])
    path = tmp_path / grid.grid_name / "ckpt0000_p0000.h5"
    _write_checkpoint(batch, path, spec_digest=spec.digest)

    result = _finalize_reconstruction(
        spec, {grid.grid_name: [path]}, tmp_path / "out.h5"
    )

    assert result["chunks_written"] == 1
    with h5py.File(result["path"], "r") as h5file:
        group = h5file[f"entry/reconstruction/results/{grid.grid_name}"]
        assert group["intensity"][1, 1, 1] == 5.0
        assert group["weight"][1, 1, 1] == 2.0
        assert group["contributors"][1, 1, 1] == 1
        assert np.isnan(group["intensity"][0, 0, 0])


def test_finalize_merges_overlapping_records_across_checkpoints(tmp_path):
    """Two checkpoints (different frame ranges) both touching the same
    voxel must have their contributions summed, not overwritten -- unlike
    the retired reduce-phase shards, checkpoints are not deduplicated
    against each other before finalize."""
    grid = _grid()
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)
    # Global voxel (1,1,1) is local_voxel_id 7 (see the single-checkpoint
    # test above for the derivation).
    first = _batch(0, [7], [3.0], contributors=[1])
    second = _batch(0, [7], [7.0], contributors=[1])
    path_a = tmp_path / grid.grid_name / "ckpt0000_p0000.h5"
    path_b = tmp_path / grid.grid_name / "ckpt0001_p0000.h5"
    _write_checkpoint(first, path_a, spec_digest=spec.digest)
    _write_checkpoint(second, path_b, spec_digest=spec.digest)

    result = _finalize_reconstruction(
        spec, {grid.grid_name: [path_a, path_b]}, tmp_path / "merged.h5"
    )

    with h5py.File(result["path"], "r") as h5file:
        group = h5file[f"entry/reconstruction/results/{grid.grid_name}"]
        assert group["weight"][1, 1, 1] == 2.0
        assert group["intensity"][1, 1, 1] == 5.0  # (3+7)/2
        assert group["contributors"][1, 1, 1] == 2


def test_finalize_rejects_unknown_grid(tmp_path):
    grid = _grid()
    spec = _ReconstructionSpec(grids=(grid,), max_depth=0)

    with pytest.raises(ValueError, match="Unknown grid"):
        _finalize_reconstruction(
            spec, {"not-a-real-grid": []}, tmp_path / "out.h5"
        )

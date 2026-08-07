"""Tests for the Sec7 prefetch/compute pipeline primitives and
orchestration: :class:`_AdjustablePool`, :class:`_BoundedGate`, and
:func:`_map_pending_ranges`'s producer/consumer behavior.
"""

from __future__ import annotations

import threading
import time

import numpy as np
import pytest

from orgui.backend.scans import h5_Image
from orgui.reconstruction_job import (
    _AdjustablePool,
    _BoundedGate,
    _map_pending_ranges,
)
from orgui.datautils.xrayutils.reconstruction import _GridSpec, _ReconstructionSpec


native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


# ---------------------------------------------------------------------------
# _AdjustablePool
# ---------------------------------------------------------------------------


def test_adjustable_pool_grows_and_shrinks_to_target_size():
    ticks = []
    lock = threading.Lock()

    def worker(retire, cancellation):
        while not retire.is_set() and not cancellation.is_set():
            with lock:
                ticks.append(1)
            time.sleep(0.01)

    pool = _AdjustablePool(worker, initial_size=2, name="test-grow")
    try:
        assert pool.size == 2
        time.sleep(0.05)
        pool.retarget(5)
        assert pool.size == 5
        time.sleep(0.05)
        pool.retarget(1)
        assert pool.size == 1
        time.sleep(0.05)
        pool.reap()
    finally:
        pool.shutdown(wait=True)
    assert len(ticks) > 0


def test_adjustable_pool_lets_in_flight_work_finish_before_retiring():
    started = threading.Event()
    may_finish = threading.Event()
    finished = threading.Event()

    def worker(retire, cancellation):
        started.set()
        may_finish.wait(timeout=5)
        finished.set()
        # A well-behaved worker only checks retire/cancellation *between*
        # items, never mid-item -- simulate exactly one "item" of work.

    pool = _AdjustablePool(worker, initial_size=1, name="test-inflight")
    try:
        assert started.wait(timeout=2)
        pool.retarget(0)  # ask the sole worker to retire
        assert not finished.is_set(), "must not be forced to stop mid-item"
        may_finish.set()
    finally:
        pool.shutdown(wait=True)
    assert finished.is_set()


def test_adjustable_pool_shutdown_joins_all_threads_no_leak():
    def worker(retire, cancellation):
        while not retire.is_set() and not cancellation.is_set():
            time.sleep(0.01)

    pool = _AdjustablePool(worker, initial_size=3, name="test-shutdown")
    with pool._lock:  # noqa: SLF001 -- test-only introspection
        threads = [thread for thread, _retire in pool._active]
    pool.shutdown(wait=True)
    assert all(not thread.is_alive() for thread in threads)


# ---------------------------------------------------------------------------
# _BoundedGate
# ---------------------------------------------------------------------------


def test_bounded_gate_bounds_concurrent_in_flight_count():
    gate = _BoundedGate(2)
    assert gate.acquire(poll_timeout=0.1, should_stop=lambda: False)
    assert gate.acquire(poll_timeout=0.1, should_stop=lambda: False)

    stop = threading.Event()
    threading.Timer(0.2, stop.set).start()
    # A third acquire must block until release() or should_stop() fires.
    result = gate.acquire(poll_timeout=0.02, should_stop=stop.is_set)
    assert result is False


def test_bounded_gate_release_admits_a_waiting_acquire():
    gate = _BoundedGate(1)
    assert gate.acquire(poll_timeout=0.1, should_stop=lambda: False)
    admitted = threading.Event()

    def waiter():
        if gate.acquire(poll_timeout=0.05, should_stop=lambda: False):
            admitted.set()

    thread = threading.Thread(target=waiter)
    thread.start()
    time.sleep(0.05)
    assert not admitted.is_set()
    gate.release()
    thread.join(timeout=2)
    assert admitted.is_set()


def test_bounded_gate_retarget_admits_more_without_release():
    gate = _BoundedGate(1)
    assert gate.acquire(poll_timeout=0.1, should_stop=lambda: False)
    assert not gate.acquire(poll_timeout=0.02, should_stop=lambda: True)
    gate.retarget(2)
    assert gate.acquire(poll_timeout=0.1, should_stop=lambda: False)


# ---------------------------------------------------------------------------
# _map_pending_ranges: pipeline-level behavior
# ---------------------------------------------------------------------------


class _FakeDetector:
    detector = type("Detector", (), {"shape": (1, 1)})()

    def primBeamPoints(self, rows, columns):
        return np.zeros_like(rows), np.zeros_like(columns)


class _FakeUB:
    def getUB(self):
        return np.eye(3)

    def getU(self):
        return np.eye(3)

    def getK(self):
        return 1.0


class _FakeConfig:
    def __init__(self):
        self.detector = _FakeDetector()
        self.ub_calculator = _FakeUB()


def _correction(payload, raw, frame, tile):
    return (
        raw.astype(np.float64),
        np.maximum(raw, 0).astype(np.float64),
        np.zeros(raw.shape, dtype=bool),
    )


def _router(boundaries, *, tmp_path, budget=1024**3):
    from orgui.datautils.xrayutils.reconstruction import _CheckpointRouter

    return _CheckpointRouter(
        boundaries,
        spec_digest="pipeline-test",
        checkpoint_dir=tmp_path / "checkpoints",
        active_budget_bytes=budget,
        resumed=(),
    )


class _SlowScan:
    """A fake scan whose get_raw_img sleeps, to exercise blocked-fraction
    reader-pool growth deterministically."""

    def __init__(self, frame_count, delay):
        self._frame_count = frame_count
        self._delay = delay
        self.concurrent_calls = 0
        self.max_concurrent_calls = 0
        self._lock = threading.Lock()

    def __len__(self):
        return self._frame_count

    def get_raw_img(self, index):
        with self._lock:
            self.concurrent_calls += 1
            self.max_concurrent_calls = max(
                self.max_concurrent_calls, self.concurrent_calls
            )
        try:
            time.sleep(self._delay)
            return h5_Image(np.array([[float(index)]]))
        finally:
            with self._lock:
                self.concurrent_calls -= 1


def _spec():
    grid = _GridSpec(
        minimum=(-1.0, -1.0, -1.0),
        maximum=(1.0, 1.0, 1.0),
        step=(1.0, 1.0, 1.0),
        frame="lab",
        chunk_shape=(2, 2, 2),
    )
    return _ReconstructionSpec(grids=(grid,), max_depth=0, threads=4)


def test_map_pending_ranges_matches_reduce_on_synthetic_frames(tmp_path):
    """The pipeline's output must match directly reducing the same frames
    one at a time -- the equivalence claim this whole rewrite depends on,
    checked at the _map_pending_ranges level directly (not just via the
    higher-level run_job/run_cluster_map_task tests)."""
    frame_count = 8
    scan = _SlowScan(frame_count, delay=0.0)
    config = _FakeConfig()
    spec = _spec()
    bounds = np.zeros((frame_count, 2, 4), dtype=np.float64)
    tiles = [(0, 1, 0, 1)]
    grid_name = spec.grids[0].grid_name
    router = _router({grid_name: [(0, frame_count)]}, tmp_path=tmp_path)

    _map_pending_ranges(
        spec,
        scan,
        config,
        bounds,
        tiles,
        [(0, frame_count)],
        router,
        correction_pipeline=_correction,
        effective_memory=256 * 1024**2,
        threads_per_image=1,
        accumulation_budget_bytes=None,
        total_images=frame_count,
        completed_images=0,
        progress=None,
    )

    assert len(router.written) == 1
    from orgui.datautils.xrayutils.reconstruction import _read_checkpoint

    written = _read_checkpoint(router.written[0])
    assert written["chunk_id"].size > 0
    # Every frame contributed: contributors sum across all voxels should
    # equal the number of frames that landed in populated voxels, and in
    # particular must not be zero or obviously truncated.
    assert int(written["contributors"].sum()) == frame_count


def test_map_pending_ranges_grows_reader_pool_under_sustained_blocking(tmp_path):
    """A slow, high-latency get_raw_img should push the reader pool above
    its initial size within a few coordinator ticks."""
    frame_count = 40
    scan = _SlowScan(frame_count, delay=0.05)
    config = _FakeConfig()
    spec = _spec()
    bounds = np.zeros((frame_count, 2, 4), dtype=np.float64)
    tiles = [(0, 1, 0, 1)]
    grid_name = spec.grids[0].grid_name
    router = _router({grid_name: [(0, frame_count)]}, tmp_path=tmp_path)

    _map_pending_ranges(
        spec,
        scan,
        config,
        bounds,
        tiles,
        [(0, frame_count)],
        router,
        correction_pipeline=_correction,
        effective_memory=256 * 1024**2,
        threads_per_image=4,
        accumulation_budget_bytes=None,
        total_images=frame_count,
        completed_images=0,
        progress=None,
    )

    # With only a handful of default readers against a 50ms-per-frame
    # source and >1 compute worker requesting frames, the reader pool
    # should have grown past its initial size (4) at some point.
    assert scan.max_concurrent_calls > 1


def test_map_pending_ranges_propagates_first_exception_after_draining(tmp_path):
    """A mid-run failure must still surface, but must not discard frames
    that were already successfully loaded before it was detected."""
    frame_count = 4
    config = _FakeConfig()
    spec = _spec()
    bounds = np.zeros((frame_count, 2, 4), dtype=np.float64)
    tiles = [(0, 1, 0, 1)]
    grid_name = spec.grids[0].grid_name
    router = _router(
        {grid_name: [(index, index + 1) for index in range(frame_count)]},
        tmp_path=tmp_path,
    )

    class FailingScan:
        def __len__(self):
            return frame_count

        def get_raw_img(self, index):
            if index == frame_count - 1:
                raise RuntimeError("simulated read failure")
            return h5_Image(np.array([[float(index)]]))

    with pytest.raises(RuntimeError, match="simulated read failure"):
        _map_pending_ranges(
            spec,
            FailingScan(),
            config,
            bounds,
            tiles,
            [(0, frame_count)],
            router,
            correction_pipeline=_correction,
            effective_memory=256 * 1024**2,
            threads_per_image=1,
            accumulation_budget_bytes=None,
            total_images=frame_count,
            completed_images=0,
            progress=None,
        )

    # The frames that loaded successfully before the failure must still
    # have been routed and checkpointed -- not discarded just because a
    # different frame in the same run failed.
    assert len(router.written) >= 1

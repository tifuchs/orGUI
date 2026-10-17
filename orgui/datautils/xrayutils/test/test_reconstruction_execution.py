"""The opt-in runtime preserves checkpoint coverage and bounded failure behavior."""

import os
import threading
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.datautils.xrayutils import reconstruction as core
from orgui.datautils.xrayutils.test import test_reconstruction_pipeline as pipeline


@pytest.mark.parametrize("frames_per_group", [1, 2])
@pytest.mark.parametrize("whole_frame_correction", [False, True])
def test_pooled_mapping_matches_legacy_checkpoints(
    tmp_path,
    monkeypatch,
    frames_per_group,
    whole_frame_correction,
):
    """Compare legacy versus pooled/filtered paths through actual checkpoint I/O."""
    original = pipeline._map_pending_ranges

    def run(*args, **kwargs):
        if not os.environ.get("ORGUI_NO_FRAME_SKIP"):
            kwargs["native_execution"] = {
                "fine_task_samples": 64,
                "max_inflight_calls": 2,
                "scratch_bytes": 1024,
            }
        return original(*args, **kwargs)

    monkeypatch.setattr(pipeline, "_map_pending_ranges", run)
    pipeline.test_tile_reach_skips_empty_calls_without_changing_checkpoints(
        tmp_path,
        monkeypatch,
        frames_per_group,
        whole_frame_correction,
    )


def test_overlapping_calls_are_bounded_and_failure_routes_nothing():
    """A failing call must drain its peers before releasing the group's inputs."""
    execution = core._NativeExecution(2, max_inflight_calls=2, scratch_bytes=0)
    first_started = threading.Event()
    second_finished = threading.Event()
    retired = threading.Event()
    calls = []
    routed = []

    class Kernel:
        def __init__(self, index):
            self.index = index

        def accumulate_group_tile(self, *args):
            """Hold one call open until a concurrent call reports failure."""
            calls.append(self.index)
            if self.index == 0:
                first_started.set()
                assert second_finished.wait(5)
                retired.set()
                return core._empty_batch()
            assert first_started.wait(5)
            second_finished.set()
            raise RuntimeError("injected mapping failure")

    grids = [SimpleNamespace(grid_name=str(i)) for i in range(3)]
    values = np.ones((1, 1))
    tile = (0, 1, 0, 1)
    try:
        with pytest.raises(RuntimeError, match="injected mapping failure"):
            core._map_frame_group(
                SimpleNamespace(grids=grids),
                {str(i): Kernel(i) for i in range(3)},
                {tile: np.ones((2, 2, 3))},
                [tile],
                None,
                None,
                [0],
                np.zeros((1, 4)),
                np.zeros((1, 4)),
                SimpleNamespace(route=lambda *a, **k: routed.append(a)),
                corrected_frames=[(values, values, np.zeros_like(values, dtype=bool))],
                native_execution=execution,
            )
        assert retired.is_set()
        assert not routed
        # The third call can be admitted only after the first has retired.
        assert calls[:2] == [0, 1]
    finally:
        execution.close()


@pytest.mark.parametrize(
    "options",
    [
        {"fine_task_samples": -1},
        {"scratch_bytes": -1},
        {"max_inflight_calls": 0},
        {"fine_task_samples": 1.5},
    ],
)
def test_invalid_execution_limits(options):
    """Invalid budgets fail before any mapping or worker allocation."""
    with pytest.raises(ValueError):
        core._NativeExecution(2, **options)

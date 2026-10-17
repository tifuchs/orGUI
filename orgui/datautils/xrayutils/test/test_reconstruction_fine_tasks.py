"""Coarse rejection and fine native tasks preserve pixel-level integration."""

import numpy as np
import pytest

from orgui.datautils.xrayutils.test.test_reconstruction_group import (
    COLUMNS,
    ROWS,
    UB,
    WAVEVECTOR,
    _angles,
    _grid,
    _payload,
    _rays,
    native,
)


def _make_kernel(minimum, step, shape, depth, threads, fine, weighting):
    return native.ReconstructionKernel(
        minimum,
        step,
        shape,
        np.array([16, 16, 16], dtype=np.int64),
        "hkl",
        WAVEVECTOR,
        np.ascontiguousarray(np.linalg.inv(UB)),
        np.eye(3),
        depth,
        threads,
        256,
        1 << 30,
        weighting,
        fine_task_samples=fine,
    )


def _compare(actual, expected, exact=False):
    for key in ("chunk_id", "local_voxel_id", "contributors"):
        np.testing.assert_array_equal(actual[key], expected[key])
    for key in ("weight", "weighted_intensity", "weighted_variance"):
        if exact:
            np.testing.assert_array_equal(actual[key], expected[key])
        else:
            np.testing.assert_allclose(
                actual[key], expected[key], rtol=1e-10, atol=1e-12
            )


@pytest.mark.parametrize("fine", [1, 64, 128, 256])
@pytest.mark.parametrize("exposure", ["stationary", "swept", "mixed"])
@pytest.mark.parametrize("depth", [0, 2])
def test_fine_tasks_match_legacy_and_are_deterministic(fine, exposure, depth):
    """Work completion order must not affect the reduction order."""
    rays = _rays()
    start, end = _angles(3, exposure != "stationary")
    if exposure == "mixed":
        end[1] = start[1]
    grid = _grid(rays, start, end)
    intensity, variance, mask = _payload(3)
    args = (intensity, variance, mask, rays, start, end)
    reference = _make_kernel(*grid, depth, 1, 0, "parameter_average").accumulate_group(
        *args
    )
    outputs = []
    for threads in (1, 4):
        kernel = _make_kernel(*grid, depth, threads, fine, "parameter_average")
        result = kernel.accumulate_group(*args, True)
        _compare(result, reference)
        assert result["_profile"]["pixels_seen"] == 3 * ROWS * COLUMNS
        assert kernel.configuration["fine_task_samples"] == fine
        outputs.append(result)
    _compare(*outputs, exact=True)


@pytest.mark.parametrize(
    "weighting", ["parameter_average", "reciprocal_volume_average"]
)
@pytest.mark.parametrize("kind", ["edge", "outside", "masked"])
def test_fine_tasks_rejection_and_tile_entry_point(weighting, kind):
    """Check narrow-grid boundaries, empty results, masks and tile offsets."""
    rays = _rays()
    start, end = _angles(2, True)
    minimum, step, shape = _grid(rays, start, end)
    if kind == "outside":
        minimum = minimum + 100
    elif kind == "edge":
        shape = np.maximum(shape // 2, 1)
    intensity, variance, mask = _payload(2)
    if kind == "masked":
        mask[:] = True
    # Non-zero source offsets exercise the non-copying group tile path.
    args = (
        list(intensity),
        list(variance),
        list(mask),
        np.ascontiguousarray(rays[1:-1, 2:-2]),
        start,
        end,
        1,
        ROWS - 1,
        2,
        COLUMNS - 2,
        True,
    )
    legacy = _make_kernel(minimum, step, shape, 2, 1, 0, weighting)
    candidate = _make_kernel(minimum, step, shape, 2, 4, 64, weighting)
    expected = legacy.accumulate_group_tile(*args)
    actual = candidate.accumulate_group_tile(*args)
    _compare(actual, expected)
    if kind == "edge":
        assert actual["weight"].size > 0
    else:
        assert actual["weight"].size == 0
    if kind == "outside":
        profile = actual["_profile"]
        assert profile["bricks"] == 0
        assert profile["rejected_coarse_bricks"] == profile["coarse_bricks"]
        assert profile["pixels_seen"] == 2 * (ROWS - 2) * (COLUMNS - 4)

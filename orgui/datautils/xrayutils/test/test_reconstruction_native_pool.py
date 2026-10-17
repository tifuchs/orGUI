"""Shared native workers preserve ordering and have bounded, explicit lifetimes."""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

from orgui.datautils.xrayutils.test.test_reconstruction_fine_tasks import _compare
from orgui.datautils.xrayutils.test.test_reconstruction_group import (
    UB,
    WAVEVECTOR,
    _angles,
    _grid,
    _payload,
    _rays,
    native,
)


def _case(context=None, fine=64):
    rays = _rays()
    start, end = _angles(3, True)
    minimum, step, shape = _grid(rays, start, end)
    kernel = native.ReconstructionKernel(
        minimum,
        step,
        shape,
        np.array([16, 16, 16], dtype=np.int64),
        "hkl",
        WAVEVECTOR,
        np.ascontiguousarray(np.linalg.inv(UB)),
        np.eye(3),
        2,
        4,
        256,
        1 << 30,
        fine_task_samples=fine,
        execution_context=context,
    )
    return kernel, (*_payload(3), rays, start, end)


@pytest.mark.parametrize("threads", [1, 4])
@pytest.mark.parametrize("scratch", [0, 1024, 1 << 20])
def test_shared_pool_concurrent_calls_match_serial(threads, scratch):
    """More calling threads than workers must not deadlock or change sums."""
    reference, args = _case()
    expected = reference.accumulate_group(*args)
    context = native.NativeExecutionContext(threads, 2, scratch)
    try:
        kernel, args = _case(context)
        assert context.threads == threads
        assert context.scratch_per_worker * threads <= scratch
        with ThreadPoolExecutor(max_workers=6) as callers:
            futures = [
                callers.submit(kernel.accumulate_group, *args) for _ in range(12)
            ]
            for future in futures:
                _compare(future.result(timeout=20), expected, exact=True)
    finally:
        context.close()
    context.close()  # Idempotent shutdown.
    with pytest.raises(RuntimeError, match="closed"):
        kernel.accumulate_group(*args)


@pytest.mark.parametrize("operation", ["close", "cancel"])
def test_shutdown_during_concurrent_submissions(operation):
    """Admission, drain and cancellation races must return to every caller."""
    context = native.NativeExecutionContext(2, 1, 1024)
    kernel, args = _case(context, fine=1)
    expected = kernel.accumulate_group(*args)
    with ThreadPoolExecutor(max_workers=6) as callers:
        futures = [callers.submit(kernel.accumulate_group, *args) for _ in range(20)]
        getattr(context, operation)()
        for future in futures:
            try:
                result = future.result(timeout=20)
            except RuntimeError as error:
                assert "closed" in str(error) or "cancelled" in str(error)
            else:
                _compare(result, expected, exact=True)
    context.close()


def test_context_is_kept_alive_by_kernel():
    """The kernel retains the pool while native pointers are in use."""
    context = native.NativeExecutionContext(2)
    kernel, args = _case(context)
    del context
    assert kernel.accumulate_group(*args)["weight"].size > 0


@pytest.mark.parametrize("threads, calls", [(0, 1), (1, 0), (1025, 1)])
def test_invalid_pool_configuration(threads, calls):
    """Reject invalid limits before creating workers."""
    with pytest.raises(ValueError):
        native.NativeExecutionContext(threads, calls)

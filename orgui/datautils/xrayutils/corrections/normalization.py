# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
r"""Counting time and incident-flux normalization.

Both of Vlieg's integrated-intensity expressions carry :math:`\Phi_0 T`, so a
frame means nothing until it is divided by its counting time and by whatever
counter the incident flux is tracked with. Two measurements are only
comparable -- to each other, and to an absolute scale -- once both have been.

This module is the policy: what a usable counter is, and what the divisor is.
It takes **values**, not a scan object. Finding the counters in a particular
beamline's scan is a data-format question and belongs to
:mod:`orgui.backend.scans` and the application layer, which is why no scan
attribute name appears here.
"""

import numpy as np

__all__ = [
    "broadcast_counter",
    "normalization_divisor",
]


def broadcast_counter(value, size, name):
    """Broadcast a scalar or per-frame counter to ``size`` values.

    :param value: A scalar, or one value per frame.
    :param int size: Number of frames.
    :param str name: Counter name, used in the error message.
    :returns: An array of shape ``(size,)``.
    :rtype: numpy.ndarray
    :raises ValueError: If the counter has neither one value nor one value
        per frame.
    """
    array = np.atleast_1d(np.asarray(value, dtype=np.float64)).ravel()
    if array.size == 1:
        return np.full(int(size), array[0], dtype=np.float64)
    if array.size != int(size):
        raise ValueError(
            f"Counter {name!r} has {array.size:d} values for {int(size):d} images"
        )
    return array


def normalization_divisor(size, exposure_time=None, monitors=None):
    """Per-frame divisor from the counting time and the monitor counters.

    :param int size: Number of frames.
    :param exposure_time: Counting time of every frame, in seconds; a scalar
        or one value per frame. ``None`` leaves the counting time out, which
        is what a scan that does not report one gets: the reciprocal-space
        reconstruction records the normalization as unavailable rather than
        failing the job, and this follows it.
    :param monitors: Mapping of counter name to its values, a scalar or one
        value per frame. Applied in iteration order.
    :returns: ``(divisor, applied)`` -- an array of shape ``(size,)`` and the
        names of the normalizations that contributed, as ``"exposure"`` and
        ``"monitor:<name>"``.
    :rtype: tuple
    :raises ValueError: If a counter has a non-positive, zero or non-finite
        value that would make the normalization undefined, or the wrong
        number of values.
    """
    divisor = np.ones(int(size), dtype=np.float64)
    applied = []

    if exposure_time is not None:
        values = broadcast_counter(exposure_time, size, "exposure_time")
        if np.any(values <= 0) or not np.all(np.isfinite(values)):
            raise ValueError("Exposure time must be finite and positive")
        divisor *= values
        applied.append("exposure")

    for name, value in (monitors or {}).items():
        values = broadcast_counter(value, size, name)
        if np.any(values == 0) or not np.all(np.isfinite(values)):
            raise ValueError(f"Monitor {name} must be finite and nonzero")
        divisor *= values
        applied.append(f"monitor:{name}")

    return divisor, applied

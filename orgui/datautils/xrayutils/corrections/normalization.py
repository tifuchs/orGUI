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
    "MONITOR_INTEGRATED",
    "MONITOR_RATE",
    "broadcast_counter",
    "frame_fluence",
    "normalization_divisor",
    "relative_frame_fluence",
]

MONITOR_RATE = "rate"
MONITOR_INTEGRATED = "integrated"


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


def _positive(value, name):
    """Return ``value`` as float array after physical-value validation."""
    array = np.asarray(value, dtype=np.float64)
    if np.any(array <= 0) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite and positive")
    return array


def _monitor_inputs(monitor, monitor_kind):
    """Validate the monitor-kind pair shared by fluence helpers."""
    if monitor_kind is None:
        if monitor is not None:
            raise ValueError("monitor_kind is required when a monitor is given")
        return None
    if monitor_kind not in (MONITOR_RATE, MONITOR_INTEGRATED):
        raise ValueError(
            f"monitor_kind must be {MONITOR_RATE!r} or "
            f"{MONITOR_INTEGRATED!r}, got {monitor_kind!r}"
        )
    if monitor is None:
        raise ValueError(f"a {monitor_kind} monitor value is required")
    return _positive(monitor, "monitor")


def frame_fluence(
    total_flux,
    exposure_time=None,
    monitor=None,
    monitor_kind=None,
    reference_monitor=None,
    reference_exposure=None,
):
    r"""Calibrated incident photons :math:`Q_f` in every frame.

    ``total_flux`` is the calibrated full beam at the sample position before
    sample clipping, in photons/s. With no monitor it is assumed constant and
    multiplied by the frame exposure. A rate-like monitor scales that product
    by ``monitor / reference_monitor``. An integrated monitor already contains
    the frame exposure, so it instead uses the calibration exposure exactly
    once:

    .. math::

        Q_f = \Phi_\mathrm{ref} T_f M_f/M_\mathrm{ref}

    for a rate, and

    .. math::

        Q_f = \Phi_\mathrm{ref} T_\mathrm{ref} U_f/U_\mathrm{ref}

    for an integrated reading. Passing ``exposure_time`` in the integrated
    case is allowed for a common calling interface but does not multiply the
    result.

    :param total_flux: Reference total incident flux, photons/s.
    :param exposure_time: Per-frame exposure, seconds. Required without a
        monitor and for a rate-like monitor.
    :param monitor: Per-frame primary-monitor reading.
    :param str monitor_kind: ``"rate"``, ``"integrated"`` or ``None``.
    :param reference_monitor: Monitor reading associated with ``total_flux``.
    :param reference_exposure: Exposure of an integrated reference reading,
        seconds. Required only for an integrated monitor.
    :returns: Incident photons per frame, broadcast over all inputs.
    :rtype: numpy.ndarray
    :raises ValueError: If a required input is absent, non-finite or not
        positive, or the monitor kind is invalid.
    """
    flux = _positive(total_flux, "total_flux")
    monitor_values = _monitor_inputs(monitor, monitor_kind)

    if monitor_kind is None:
        if reference_monitor is not None or reference_exposure is not None:
            raise ValueError("monitor references require a monitor")
        if exposure_time is None:
            raise ValueError("exposure_time is required without a monitor")
        return flux * _positive(exposure_time, "exposure_time")

    if reference_monitor is None:
        raise ValueError("reference_monitor is required for a calibrated monitor")
    reference = _positive(reference_monitor, "reference_monitor")

    if monitor_kind == MONITOR_RATE:
        if exposure_time is None:
            raise ValueError("exposure_time is required for a rate monitor")
        if reference_exposure is not None:
            raise ValueError(
                "reference_exposure applies only to an integrated monitor"
            )
        return (
            flux
            * _positive(exposure_time, "exposure_time")
            * monitor_values
            / reference
        )

    if reference_exposure is None:
        raise ValueError(
            "reference_exposure is required for an integrated monitor"
        )
    return (
        flux
        * _positive(reference_exposure, "reference_exposure")
        * monitor_values
        / reference
    )


def relative_frame_fluence(
    exposure_time=None,
    monitor=None,
    monitor_kind=None,
):
    r"""Uncalibrated relative incident fluence for every frame.

    This preserves variations between frames without claiming photon units.
    With no monitor the relative divisor is the exposure. For a rate monitor
    it is ``exposure_time * monitor``; for an integrated monitor it is the
    monitor reading alone, because that reading already integrated over the
    frame exposure.

    :param exposure_time: Per-frame exposure, seconds.
    :param monitor: Per-frame primary-monitor reading.
    :param str monitor_kind: ``"rate"``, ``"integrated"`` or ``None``.
    :returns: Relative per-frame fluence divisor.
    :rtype: numpy.ndarray
    :raises ValueError: If required inputs are absent or nonphysical.
    """
    monitor_values = _monitor_inputs(monitor, monitor_kind)
    if monitor_kind is None:
        if exposure_time is None:
            raise ValueError("exposure_time is required without a monitor")
        return _positive(exposure_time, "exposure_time")
    if monitor_kind == MONITOR_RATE:
        if exposure_time is None:
            raise ValueError("exposure_time is required for a rate monitor")
        return _positive(exposure_time, "exposure_time") * monitor_values
    return monitor_values

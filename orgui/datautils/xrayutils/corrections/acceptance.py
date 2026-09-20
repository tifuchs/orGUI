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
r"""How much of a diffraction rod a region of interest accepts.

A rocking scan does not measure a whole rod. It intercepts a slice of one,
and the slice is as long as the detector aperture is wide in the out-of-plane
direction: :math:`\Delta l = C_\mathrm{rod}\,(V_u / \lambda A_u)\,\Delta\gamma`
(E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532, equation 20). Its integrated
intensity is therefore **proportional to** :math:`\Delta\gamma`, which is why
:func:`~.measurement.angular_factor` refuses to reduce a rocking scan without
one. A stationary measurement integrates across the whole rod cross-section
and has no such factor.

With a point detector :math:`\Delta\gamma` was a slit setting: fixed for a
whole experiment, and absorbed into the overall scale factor nobody had to
know. On an area detector it is a property of the *region of interest*, and
orGUI sizes regions of interest per detector position
(:func:`orgui.app.ROIutils.calc_corrections` scales them with the projected
sample size and a parallax correction). It is therefore not constant along a
scan, and leaving it out changes the shape of a rod, not only its scale.

What this module computes is the angular height of the region of interest as
seen from the sample, in the surface frame, from the calibrated detector
geometry. It is a geometric estimate: it assumes the region is centered on
the rod and that the rod runs along the region's vertical direction. Both
hold for the regions orGUI places from a calculated reflection position, and
:func:`gamma_range` measures how well by reporting the span over the whole
region.

The in-plane counterpart, Vlieg's :math:`C_\mathrm{det}` (his section 2.4),
is the other half of the detector acceptance and is not implemented yet; a
region wide enough to contain the whole in-plane peak profile has
:math:`C_\mathrm{det} = 1`, which is what the reduction assumes.

.. note::

    Evaluate the acceptance at the actual detector-arm position whenever it
    is known. Arm invariance is not general: a rolled or oblique detector
    changes the projection of the pixel-edge rays onto :math:`\gamma`, and
    realistic counterexamples differ by several percent. Omitting the arm is
    retained only as an explicit compatibility fallback for measurements that
    predate stored arm metadata.
"""

import numpy as np

__all__ = [
    "gamma_range",
    "out_of_plane_acceptance",
    "pixel_acceptance",
]


def _surface_gamma(detector, row, column, alpha, gamma_arm=None, delta_arm=None):
    """Surface-frame exit angle at the given detector points, in radian.

    Isolates the one call into
    :meth:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD.surfaceAnglesPoint`,
    whose first argument is pyFAI dimension 1 -- the detector **row**, not a
    horizontal coordinate -- so the ordering is stated once rather than at
    every call site.

    The two coordinates are broadcast against each other first. pyFAI asserts
    that they have the same size, so a scan of per-frame region heights at
    one fixed column -- which is what a rocking integration passes -- would
    otherwise fail inside the extension with ``pos2.size == size``.
    """
    if (gamma_arm is None) != (delta_arm is None):
        raise ValueError("both gamma_arm and delta_arm must be given together")
    row, column = np.broadcast_arrays(
        np.asarray(row, dtype=np.float64),
        np.asarray(column, dtype=np.float64),
    )
    if gamma_arm is None or (
        np.ndim(gamma_arm) == 0 and np.ndim(delta_arm) == 0
    ):
        gamma, _delta = detector.surfaceAnglesPoint(
            np.ascontiguousarray(row),
            np.ascontiguousarray(column),
            alpha,
            gamma_arm,
            delta_arm,
        )
        return np.asarray(gamma, dtype=np.float64)

    # ``surfaceAnglesPoint`` accepts one arm geometry at a time. A rocking
    # reduction instead supplies one ROI and one arm position per curve, so
    # evaluate those pairs rather than forming every arm-by-ROI combination.
    row, column, alpha, gamma_arm, delta_arm = np.broadcast_arrays(
        row,
        column,
        np.asarray(alpha, dtype=np.float64),
        np.asarray(gamma_arm, dtype=np.float64),
        np.asarray(delta_arm, dtype=np.float64),
    )
    result = np.empty(row.shape, dtype=np.float64)
    for index in np.ndindex(row.shape):
        gamma, _delta = detector.surfaceAnglesPoint(
            np.array([row[index]]),
            np.array([column[index]]),
            float(alpha[index]),
            float(gamma_arm[index]),
            float(delta_arm[index]),
        )
        result[index] = np.asarray(gamma, dtype=np.float64).reshape(-1)[0]
    return result


def out_of_plane_acceptance(
    detector,
    row,
    column,
    row_size,
    alpha,
    gamma_arm=None,
    delta_arm=None,
):
    r"""Out-of-plane angular acceptance :math:`\Delta\gamma`, in radian.

    The exit-angle span between the top and bottom **edges** of a region of
    interest, evaluated at its center column. Edges, not pixel centers: a
    region of ``row_size`` rows accepts photons from ``row - row_size/2`` to
    ``row + row_size/2``, and using the centers would understate the
    acceptance by one pixel.

    The center column is the right place to measure it because the region is
    centered on the rod, and it is the range of :math:`\gamma` over which
    *the rod* crosses the aperture that sets the intercepted rod length.
    Where the detector is rotated enough that the rod does not run along the
    region's columns, this underestimates the span; :func:`gamma_range`
    reports the span over the whole region so the two can be compared.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`,
        or anything with the same ``surfaceAnglesPoint`` contract.
    :param row: Center row of the region, in pixels along pyFAI dimension 1.
        Scalar or one value per frame.
    :param column: Center column of the region, in pixels along pyFAI
        dimension 2. Scalar or one value per frame.
    :param row_size: Height of the region, in pixels. Scalar or one value per
        frame -- orGUI resizes regions along a scan, so this is normally an
        array.
    :param alpha: Incidence angle of the frame, in radian.
    :param gamma_arm: Detector arm position, in radian; ``None`` is the
        calibrated position. Give both arm angles or neither.
    :param delta_arm: Detector arm position, in radian; ``None`` is the
        calibrated position.
    :returns: :math:`\Delta\gamma` in radian, broadcast over the inputs and
        always positive.
    :rtype: numpy.ndarray
    :raises ValueError: If a region height is not positive, or exactly one
        arm angle is given.
    """
    row = np.asarray(row, dtype=np.float64)
    column = np.asarray(column, dtype=np.float64)
    row_size = np.asarray(row_size, dtype=np.float64)
    if np.any(row_size <= 0) or not np.all(np.isfinite(row_size)):
        raise ValueError(
            "the region of interest must have a positive height in pixels; "
            f"got {row_size!r}"
        )
    half = row_size / 2.0
    lower = _surface_gamma(
        detector, row - half, column, alpha, gamma_arm, delta_arm
    )
    upper = _surface_gamma(
        detector, row + half, column, alpha, gamma_arm, delta_arm
    )
    return np.abs(upper - lower)


def gamma_range(
    detector,
    row,
    column,
    row_size,
    column_size,
    alpha,
    gamma_arm=None,
    delta_arm=None,
):
    r"""Exit-angle span over a whole rectangular region, in radian.

    The corner-to-corner :math:`\gamma` extent, as opposed to the
    center-column extent :func:`out_of_plane_acceptance` returns. The two
    agree when lines of constant :math:`\gamma` run along the detector rows;
    they separate when the detector is rotated about the beam, and the ratio
    is a cheap diagnostic for whether the acceptance estimate can be trusted.

    Evaluated at the four corners: :math:`\gamma` is monotonic in both pixel
    directions over a region small enough to hold one reflection, so the
    extremes are corners.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`.
    :param row: Center row of the region, in pixels.
    :param column: Center column of the region, in pixels.
    :param row_size: Height of the region, in pixels.
    :param column_size: Width of the region, in pixels.
    :param alpha: Incidence angle of the frame, in radian.
    :param gamma_arm: Detector arm position, in radian.
    :param delta_arm: Detector arm position, in radian.
    :returns: The full span, broadcast over the inputs and always positive.
    :rtype: numpy.ndarray
    :raises ValueError: If a region size is not positive.
    """
    row = np.asarray(row, dtype=np.float64)
    column = np.asarray(column, dtype=np.float64)
    row_size = np.asarray(row_size, dtype=np.float64)
    column_size = np.asarray(column_size, dtype=np.float64)
    for name, value in (("height", row_size), ("width", column_size)):
        if np.any(value <= 0) or not np.all(np.isfinite(value)):
            raise ValueError(
                f"the region of interest must have a positive {name} in "
                f"pixels; got {value!r}"
            )
    half_row = row_size / 2.0
    half_column = column_size / 2.0
    corners = [
        _surface_gamma(
            detector,
            row + dr * half_row,
            column + dc * half_column,
            alpha,
            gamma_arm,
            delta_arm,
        )
        for dr in (-1.0, 1.0)
        for dc in (-1.0, 1.0)
    ]
    stacked = np.stack(np.broadcast_arrays(*corners))
    return np.max(stacked, axis=0) - np.min(stacked, axis=0)


def pixel_acceptance(detector, row, column, alpha, gamma_arm=None, delta_arm=None):
    r"""Exit-angle height of a single pixel, in radian.

    The :math:`\gamma` subtended by one pixel row at the given position: the
    resolution limit of :func:`out_of_plane_acceptance`, and the natural unit
    to express a region's acceptance in. A region of ``n`` rows on a flat
    detector viewed near its normal accepts about ``n`` times this; the ratio
    departs from ``n`` exactly where the detector is oblique, which is the
    effect that makes the acceptance vary along a scan.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`.
    :param row: Row of the pixel, in pixels.
    :param column: Column of the pixel, in pixels.
    :param alpha: Incidence angle of the frame, in radian.
    :param gamma_arm: Detector arm position, in radian.
    :param delta_arm: Detector arm position, in radian.
    :returns: The angular height of that pixel, broadcast over the inputs.
    :rtype: numpy.ndarray
    """
    return out_of_plane_acceptance(
        detector, row, column, 1.0, alpha, gamma_arm, delta_arm
    )

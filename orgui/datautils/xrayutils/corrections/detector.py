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
r"""Per-pixel correction factors of a detector image.

The factors that depend on where a photon landed on the detector rather than
on the scan: the solid angle each pixel subtends and the polarization factor
at its scattering angle. Both are properties of the calibrated geometry, so
they are evaluated once and reused for every frame that shares it.

They are returned as a **multiplicative** factor, the reciprocal of the
physical quantities, because that is how a corrected intensity is formed:

.. math::

    I_\mathrm{corr} = I \cdot \frac{1}{\Omega\,P}

The rocking integration, the stationary integration and the reciprocal-space
reconstruction each built this array themselves; this is the one definition
they share.

.. warning::

    :meth:`~.DetectorCalibration.Detector2D_SXRD.polarizationArray` evaluates
    the polarization at the **calibrated** detector position. That is correct
    for a detector whose arm does not move -- every pixel already carries its
    own scattering angle -- but not for a scan that drives the arm, where it
    understates the correction badly: 3 % at a scattering angle of 10
    degrees, 10 % at 18 and 33 % at 30.
    :meth:`~.DetectorCalibration.Detector2D_SXRD.polarizationAtPoints`
    follows the arm and is what such a scan needs. This function reproduces
    the historical, arm-blind behavior; see
    ``doc/design/ctr_structure_factor_scale.md`` finding F5.
"""

import numpy as np

__all__ = [
    "pixel_factors",
    "polarization_arm_correction",
    "roi_mean_inverse_polarization",
    "roi_mean_inverse_solid_angle",
]


def pixel_factors(detector, solid_angle=False, polarization=False, shape=None):
    r"""Per-pixel multiplicative correction factor of a detector image.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`
        with its calibrated geometry and polarization set.
    :param bool solid_angle: Divide by the solid angle each pixel subtends.
    :param bool polarization: Divide by the polarization factor.
    :param shape: Detector shape; taken from the detector when omitted.
    :returns: The factor for every pixel, or ``None`` when neither
        correction is enabled -- which lets a caller skip the multiplication
        entirely rather than multiply by an array of ones.
    :rtype: numpy.ndarray or None
    """
    factor = None
    if solid_angle:
        factor = 1.0 / np.asarray(
            detector.solidAngleArray(shape) if shape is not None
            else detector.solidAngleArray(),
            dtype=np.float64,
        )
    if polarization:
        values = np.asarray(
            detector.polarizationArray(shape) if shape is not None
            else detector.polarizationArray(),
            dtype=np.float64,
        )
        factor = 1.0 / values if factor is None else factor / values
    return factor


def roi_mean_inverse_solid_angle(
    detector, row, column, row_size, column_size, shape=None
):
    r"""Mean of :math:`1/\widetilde{\Omega}` over a rectangular region.

    The solid-angle content of the per-pixel correction of
    :func:`pixel_factors`, reduced onto one region of interest. A
    region-summed intensity that has been divided by the solid angle needs
    this factor divided back out before it becomes a structure factor: the
    sum over a region is already the complete angular integral, each pixel
    weighted by the solid angle it subtends, so the correction double-counts
    the detector obliquity there. See
    ``doc/design/ctr_structure_factor_scale.md`` finding F6.

    Evaluated over the **nominal** region rectangle, from the calibrated
    geometry alone, so it needs no image and can be computed per frame outside
    an integration loop. That is deliberately not identical to the mean over
    the *valid* pixels that the integration accumulates: the two differ only
    where masked pixels correlate with the detector obliquity, and
    :math:`\widetilde{\Omega}` varies by well under a percent across one
    region.

    .. note::

        Because :func:`pixel_factors` returns the solid angle and the
        polarization as one fused array, the applied region mean is
        :math:`\langle 1/(\widetilde{\Omega}P)\rangle` rather than
        :math:`\langle 1/\widetilde{\Omega}\rangle\,\langle 1/P\rangle`.
        Dividing this factor out therefore leaves the covariance of the two
        over the region, which is second order in their variation across it --
        of order :math:`10^{-6}` for a region of a hundred pixels at a metre.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`.
    :param row: Region centre row, in pixels along pyFAI dimension 1. Scalar
        or one value per frame.
    :param column: Region centre column, in pixels along pyFAI dimension 2.
    :param row_size: Region height in pixels.
    :param column_size: Region width in pixels.
    :param shape: Detector shape; taken from the detector when omitted.
    :returns: The mean of the reciprocal normalized solid angle, broadcast
        over the inputs. A region entirely off the detector contains no pixels
        and no counts, and yields ``1.0`` so that it neither rescales nor
        invalidates the zero intensity there. So does a region whose position
        is not finite: a real scan has frames where the rod never reaches the
        detector, and those arrive here as ``inf`` or ``nan`` centres while
        carrying no counts either.
    :rtype: numpy.ndarray
    :raises ValueError: If a finite region size is not positive.
    """
    solid_angle = np.asarray(
        detector.solidAngleArray(shape) if shape is not None
        else detector.solidAngleArray(),
        dtype=np.float64,
    )
    row, column, row_size, column_size = np.broadcast_arrays(
        np.asarray(row, dtype=np.float64),
        np.asarray(column, dtype=np.float64),
        np.asarray(row_size, dtype=np.float64),
        np.asarray(column_size, dtype=np.float64),
    )
    # A non-finite position or size is a frame on which the rod never reached
    # the detector. Those carry no counts, so they are skipped rather than
    # rejected; only a *finite* size that is not positive is a real error.
    defined = (
        np.isfinite(row)
        & np.isfinite(column)
        & np.isfinite(row_size)
        & np.isfinite(column_size)
    )
    for name, value in (("height", row_size), ("width", column_size)):
        if np.any(value[defined] <= 0):
            raise ValueError(
                f"the region of interest must have a positive {name} in "
                f"pixels; got {value!r}"
            )

    n_rows, n_columns = solid_angle.shape[0], solid_angle.shape[1]
    out = np.ones(row.shape, dtype=np.float64)
    for index in np.ndindex(*row.shape):
        if not defined[index]:
            continue
        r0 = int(np.floor(row[index] - row_size[index] / 2.0))
        r1 = int(np.ceil(row[index] + row_size[index] / 2.0))
        c0 = int(np.floor(column[index] - column_size[index] / 2.0))
        c1 = int(np.ceil(column[index] + column_size[index] / 2.0))
        block = solid_angle[
            max(r0, 0):min(r1, n_rows), max(c0, 0):min(c1, n_columns)
        ]
        if block.size:
            out[index] = np.mean(1.0 / block)
    return out


def _roi_sample_grid(row, column, row_size, column_size, samples):
    """Coordinates spanning one region, at most ``samples`` per direction.

    The polarization varies smoothly across a region, so a coarse sample
    gives its mean to far better accuracy than the correction itself is
    known. Capping the count keeps the cost per frame independent of how
    large the region is.

    Samples sit at the centres of equal sub-intervals, the midpoint rule, so
    a single sample lands on the region centre rather than on an edge and no
    sample count is biased towards one side.

    :returns: ``(rows, columns)`` as 1D arrays of pixel coordinates.
    :rtype: tuple
    """

    def _centres(centre, size, count):
        count = int(min(max(int(round(size)), 1), count))
        return centre + ((np.arange(count) + 0.5) / count - 0.5) * size

    return (
        _centres(row, row_size, samples),
        _centres(column, column_size, samples),
    )


def roi_mean_inverse_polarization(
    detector,
    row,
    column,
    row_size,
    column_size,
    alpha,
    gamma_arm=None,
    delta_arm=None,
    samples=17,
):
    r"""Mean of :math:`1/P` over a region, at one detector arm position.

    Uses
    :meth:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD.polarizationAtPoints`,
    which evaluates the z-axis polarization expression in the surface-frame
    angles of each pixel **at the arm position given**. That is the difference
    from :func:`pixel_factors`, whose array comes from pyFAI's detector-frame
    expression at the *calibrated* position; see
    :func:`polarization_arm_correction`.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`.
    :param row: Region centre row, in pixels along pyFAI dimension 1.
    :param column: Region centre column, in pixels along pyFAI dimension 2.
    :param row_size: Region height in pixels.
    :param column_size: Region width in pixels.
    :param alpha: Incidence angle of the frame, in radian.
    :param gamma_arm: Detector arm position, in radian; ``None`` is the
        calibrated position. Give both arm angles or neither.
    :param delta_arm: Detector arm position, in radian.
    :param int samples: Upper bound on the sample count per direction.
    :returns: The mean of :math:`1/P` over the region.
    :rtype: float
    """
    rows, columns = _roi_sample_grid(
        float(row), float(column), float(row_size), float(column_size), samples
    )
    grid_rows, grid_columns = np.meshgrid(rows, columns, indexing="ij")
    polarization = detector.polarizationAtPoints(
        np.ascontiguousarray(grid_rows.ravel()),
        np.ascontiguousarray(grid_columns.ravel()),
        float(alpha),
        gamma_arm,
        delta_arm,
    )
    return float(np.mean(1.0 / np.asarray(polarization, dtype=np.float64)))


def polarization_arm_correction(
    detector,
    row,
    column,
    row_size,
    column_size,
    alpha,
    gamma_arm,
    delta_arm,
    samples=17,
):
    r"""Factor moving a polarization correction onto the frame's arm position.

    :func:`pixel_factors` divides by the polarization of the **calibrated**
    geometry, evaluated once outside the frame loop. For a detector whose arm
    does not move that is correct -- every pixel already carries its own
    scattering angle. For a scan that drives the arm it is not: the same pixel
    looks in a different direction on every frame, and the correction comes
    out far too small. Measured at the centre of a detector at one metre, with
    the arm following a specular scan, the calibrated-position polarization is
    high by 0.8 % at a scattering angle of 5 degrees, 3.2 % at 10, 10.7 % at
    18 and 33.6 % at 30.

    This returns
    :math:`\langle 1/P_\mathrm{arm}\rangle / \langle 1/P_\mathrm{home}\rangle`
    over the region, the factor an intensity already corrected with the
    calibrated-position polarization must be multiplied by. It is **exactly
    one** when the arm sits at its calibrated position, so a fixed-arm scan is
    untouched.

    The ratio is taken between two region means rather than at the region
    centre because the polarization is not flat across a region at a large
    scattering angle: for a 100-pixel region at one metre near
    :math:`2\theta = 30` degrees it varies by about a percent from edge to
    edge. What is left out is the covariance with the solid angle over the
    region, which :func:`pixel_factors` fuses into the same array -- second
    order in the variation of both across one region.

    See ``doc/design/ctr_structure_factor_scale.md`` finding F5.

    :param detector: A
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`.
    :param row: Region centre row, in pixels along pyFAI dimension 1.
    :param column: Region centre column, in pixels along pyFAI dimension 2.
    :param row_size: Region height in pixels.
    :param column_size: Region width in pixels.
    :param alpha: Incidence angle of the frame, in radian.
    :param gamma_arm: Detector arm position of the frame, in radian.
    :param delta_arm: Detector arm position of the frame, in radian.
    :param int samples: Upper bound on the sample count per direction.
    :returns: The multiplicative correction, ``1.0`` at the calibrated
        position and ``1.0`` where the region position is not finite -- a
        frame on which the rod never reached the detector, which carries no
        counts to correct.
    :rtype: float
    """
    if not all(
        np.isfinite(float(value))
        for value in (row, column, row_size, column_size, alpha)
    ):
        return 1.0
    at_home = roi_mean_inverse_polarization(
        detector, row, column, row_size, column_size, alpha, samples=samples
    )
    at_arm = roi_mean_inverse_polarization(
        detector,
        row,
        column,
        row_size,
        column_size,
        alpha,
        gamma_arm,
        delta_arm,
        samples=samples,
    )
    return at_arm / at_home

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
"""Reducing per-pixel corrections onto a summed region of interest.

A region of interest is summed pixel by pixel, but its correction factors are
reported and applied per region. This module holds the bookkeeping that
bridges the two, and the container the applied factors are carried in.
"""

import numpy as np

__all__ = [
    "CorrectionFactors",
    "roi_mean_correction",
]


class CorrectionFactors(dict):
    """Per-image correction divisors, plus the names of what was applied.

    A plain dict of ``name -> array`` with two conveniences: :attr:`applied`
    lists the corrections that actually contributed, and :meth:`divisor`
    multiplies a chosen subset together.
    """

    def __init__(self, factors, applied):
        super().__init__(factors)
        #: Names of the corrections that contributed, in application order.
        self.applied = tuple(applied)

    def divisor(self, *names):
        """Product of the named factors, or 1.0 when none are present.

        :param names: Factor names to multiply. Missing names are skipped,
            so a caller can ask for a correction that was not enabled.
        :rtype: numpy.ndarray or float
        """
        product = 1.0
        for name in names:
            if name in self:
                product = product * self[name]
        return product


def roi_mean_correction(correction_sum, pixel_count):
    """Mean per-pixel correction over one region of interest.

    The ROI summing accumulates the correction array over the same pixels it
    sums the counts over, giving ``correction_sum`` and the number of valid
    pixels ``pixel_count``. An ROI-summed intensity is corrected by the
    *mean* of the correction across those pixels.

    The nominal ROI area must not enter here: the integrated intensity is
    already rescaled from the valid pixels to the nominal ROI area when the
    background is subtracted. Multiplying by the area a second time scaled
    every corrected intensity by the size of its ROI, and because the
    projected ROI size varies over the detector, two measurements of one rod
    taken on different parts of it were scaled apart.

    :param correction_sum: Summed correction array over the ROI.
    :param pixel_count: Number of valid pixels contributing to that sum.
    :returns: The mean correction, and 0 where no pixel was valid.
    :rtype: numpy.ndarray
    """
    correction_sum = np.asarray(correction_sum, dtype=np.float64)
    pixel_count = np.asarray(pixel_count, dtype=np.float64)
    return np.divide(
        correction_sum,
        pixel_count,
        out=np.zeros_like(correction_sum),
        where=pixel_count > 0,
    )

"""Selected axes of segmented scans stay aligned with their images."""

import numpy as np
import pytest

from orgui.backend.interlacedScanLoader import InterlacedScan


class _Segment:
    """Minimal scan segment with independently shaped motor values."""

    def __init__(self, points, *, th, mu):
        self.nopoints = points
        self.th = th
        self.mu = mu
        self.auxillary_counters = []

    def __len__(self):
        return self.nopoints


@pytest.mark.parametrize("axis", ["th", "mu"])
def test_a_scalar_selected_axis_is_expanded_per_image(axis):
    first = _Segment(2, th=2.0, mu=np.array([2.0]))
    second = _Segment(3, th=np.array(1.0), mu=1.0)

    scan = InterlacedScan([first, second], True, axis)

    np.testing.assert_allclose(scan.axis, [1.0, 1.0, 1.0, 2.0, 2.0])
    assert scan.axis.shape == (len(scan),)


def test_an_axis_with_the_wrong_number_of_values_is_rejected():
    segment = _Segment(3, th=[1.0, 2.0], mu=0.0)

    with pytest.raises(ValueError, match="2 'th' values for 3 images"):
        InterlacedScan([segment], False, "th")

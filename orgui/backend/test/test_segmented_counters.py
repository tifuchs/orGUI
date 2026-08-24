"""A segmented scan makes the counters its segments share available.

The segmented ("interlaced") scan used to report no counters at all, so a
scan combined from segments lost the monitors, the exposure time and every
other counter its backend provides -- silently, and only for scans built in
the segmented scan loader.
"""

import logging

import numpy as np
import pytest

from orgui.backend.interlacedScanLoader import InterlacedScan, shared_counters


class _Segment:
    """A scan segment with just enough of the backend interface."""

    def __init__(self, th, counters=(), declared=None, mu=0.5, marker=0):
        self.th = np.asarray(th, dtype=np.float64)
        self.mu = mu
        self.omega = -1 * self.th
        self.axisname = "th"
        self.axis = self.th
        self.nopoints = self.th.size
        self.scanname = f"scan_{marker}"
        self._marker = marker
        self._declared = list(counters) if declared is None else list(declared)
        for name, value in dict(counters).items():
            setattr(self, name, value)

    @property
    def auxillary_counters(self):
        return list(self._declared)

    def get_raw_img(self, i):
        # identifies both the segment and the position inside it
        return self._marker + i / 100

    def __len__(self):
        return self.nopoints


def test_counters_of_all_segments_are_concatenated_in_image_order():
    first = _Segment([0.0, 1.0], {"ic1": [10.0, 11.0], "diode": [1.0, 2.0]})
    second = _Segment([2.0, 3.0], {"ic1": [12.0, 13.0], "diode": [3.0, 4.0]})

    scan = InterlacedScan([first, second], False, "th")

    assert sorted(scan.auxillary_counters) == ["diode", "ic1"]
    np.testing.assert_allclose(scan.ic1, [10.0, 11.0, 12.0, 13.0])
    np.testing.assert_allclose(scan.diode, [1.0, 2.0, 3.0, 4.0])


def test_a_counter_missing_from_one_segment_is_dropped(caplog):
    """Keeping it would misalign every value after the gap."""
    first = _Segment([0.0, 1.0], {"ic1": [10.0, 11.0], "ic2": [1.0, 2.0]})
    second = _Segment([2.0, 3.0], {"ic1": [12.0, 13.0]}, declared=["ic1", "ic2"])

    with caplog.at_level(logging.INFO, logger="orgui.backend.interlacedScanLoader"):
        scan = InterlacedScan([first, second], False, "th")

    assert scan.auxillary_counters == ["ic1"]
    assert not hasattr(scan, "ic2")
    assert "'ic2'" in caplog.text


def test_a_counter_held_once_per_segment_is_expanded_per_image():
    """Segments are separate scans and may use different exposure times."""
    first = _Segment([0.0, 1.0], {"exposure_time": 0.25})
    second = _Segment([2.0, 3.0, 4.0], {"exposure_time": np.array([1.0])})

    scan = InterlacedScan([first, second], False, "th")

    np.testing.assert_allclose(
        scan.exposure_time, [0.25, 0.25, 1.0, 1.0, 1.0]
    )


def test_exposure_time_is_merged_without_being_declared():
    """The analysis reads it by name, whether a backend declares it or not."""
    first = _Segment([0.0, 1.0], {"exposure_time": [0.25, 0.25]}, declared=[])
    second = _Segment([2.0, 3.0], {"exposure_time": [1.0, 1.0]}, declared=[])

    scan = InterlacedScan([first, second], False, "th")

    np.testing.assert_allclose(scan.exposure_time, [0.25, 0.25, 1.0, 1.0])


def test_sorted_segments_keep_their_counters_aligned_with_the_images():
    first = _Segment([0.0, 2.0], {"ic1": [10.0, 12.0]}, marker=0)
    second = _Segment([1.0, 3.0], {"ic1": [21.0, 23.0]}, marker=1)

    scan = InterlacedScan([first, second], True, "th")

    np.testing.assert_allclose(scan.th, [0.0, 1.0, 2.0, 3.0])
    # the counter of image i must come from the same segment and position
    # as the image itself
    images = [scan.get_raw_img(i) for i in range(len(scan))]
    assert images == [0.0, 1.0, 0.01, 1.01]
    np.testing.assert_allclose(scan.ic1, [10.0, 21.0, 12.0, 23.0])


def test_a_counter_that_cannot_be_aligned_with_the_images_is_dropped(caplog):
    first = _Segment([0.0, 1.0], {"ic1": [10.0, 11.0]})
    second = _Segment([2.0, 3.0], {"ic1": [12.0, 13.0, 14.0]})

    with caplog.at_level(logging.INFO, logger="orgui.backend.interlacedScanLoader"):
        scan = InterlacedScan([first, second], False, "th")

    assert scan.auxillary_counters == []
    assert "'ic1'" in caplog.text


def test_a_counter_never_overwrites_the_scan_geometry():
    """A backend counter named like a motor must not redefine the axis.

    The segments here declare their motors as counters, which a backend that
    copies motor positions into its counter list does.
    """
    declared = ["ic1", "th", "mu", "axis", "title"]
    first = _Segment([0.0, 1.0], {"ic1": [10.0, 11.0]}, declared=declared)
    second = _Segment([2.0, 3.0], {"ic1": [12.0, 13.0]}, declared=declared)

    assert list(shared_counters([first, second])) == ["ic1"]

    scan = InterlacedScan([first, second], False, "th")

    np.testing.assert_allclose(scan.th, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(scan.axis, [0.0, 1.0, 2.0, 3.0])
    assert scan.mu == 0.5
    assert scan.title == "interlaced scan"


def test_multi_column_counters_are_concatenated_along_the_image_axis():
    first = _Segment([0.0, 1.0], {"vbpm": [[1.0, 2.0], [3.0, 4.0]]})
    second = _Segment([2.0], {"vbpm": [[5.0, 6.0]]})

    counters = shared_counters([first, second])

    np.testing.assert_allclose(
        counters["vbpm"], [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    )


def test_shared_counters_of_no_segments_is_empty():
    assert shared_counters([]) == {}


@pytest.mark.parametrize("sort", [False, True])
def test_counter_values_have_one_entry_per_image(sort):
    first = _Segment([0.0, 2.0, 4.0], {"ic1": [1.0, 2.0, 3.0]})
    second = _Segment([1.0, 3.0], {"ic1": [4.0, 5.0]})

    scan = InterlacedScan([first, second], sort, "th")

    for name in scan.auxillary_counters:
        assert np.asarray(getattr(scan, name)).shape[0] == len(scan)

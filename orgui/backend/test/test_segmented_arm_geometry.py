"""A segmented scan describes the detector arm of its segments.

Without this the per-pixel conversions fall back to the fixed arm position in
the configuration, so every image of a segmented scan that moves the detector
-- an XRR scan measured in several ``a2scan th tth`` segments, say -- is
converted at the wrong place. The merge is deliberately conservative: a
mismatch between the segments is reported and leaves the arm to the
configuration, rather than being resolved by guessing.
"""

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.backend.interlacedScanLoader import InterlacedScan, shared_arm_geometry
from orgui.backend.scans import scan_arm_angles, scan_exposure_angle_bounds

LOGGER = "orgui.backend.interlacedScanLoader"


class _Segment:
    """A scan segment with just enough of the backend interface."""

    def __init__(self, th, **attributes):
        self.th = np.asarray(th, dtype=np.float64)
        self.mu = 0.5
        self.omega = -1 * self.th
        self.axisname = "th"
        self.axis = self.th
        self.nopoints = self.th.size
        self.auxillary_counters = []
        for name, value in attributes.items():
            setattr(self, name, value)

    def __len__(self):
        return self.nopoints


@pytest.fixture
def caplog_loader(caplog):
    caplog.set_level(logging.INFO, logger=LOGGER)
    return caplog


def test_a_moving_arm_is_concatenated_over_the_segments():
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], delta_arm=0.0)
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0], delta_arm=0.0)

    scan = InterlacedScan([first, second], False, "th")

    np.testing.assert_allclose(scan.gamma_arm, [1.2, 2.0, 4.0, 8.0])
    assert scan.delta_arm == 0.0


def test_sorted_segments_keep_the_arm_aligned_with_the_images():
    first = _Segment([0.0, 2.0], gamma_arm=[10.0, 12.0])
    second = _Segment([1.0, 3.0], gamma_arm=[21.0, 23.0])

    scan = InterlacedScan([first, second], True, "th")

    np.testing.assert_allclose(scan.th, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(scan.gamma_arm, [10.0, 21.0, 12.0, 23.0])


def test_one_fixed_value_shared_by_every_segment_stays_a_fixed_value():
    first = _Segment([0.0, 1.0], gamma_arm=7.5)
    second = _Segment([2.0, 3.0], gamma_arm=np.float64(7.5))

    geometry = shared_arm_geometry([first, second])

    assert geometry["gamma_arm"] == 7.5
    assert np.ndim(geometry["gamma_arm"]) == 0


def test_fixed_values_that_differ_are_expanded_per_image(caplog_loader):
    """The segments were measured at different arm positions."""
    first = _Segment([0.0, 1.0], gamma_arm=5.0)
    second = _Segment([2.0, 3.0, 4.0], gamma_arm=10.0)

    scan = InterlacedScan([first, second], False, "th")

    np.testing.assert_allclose(scan.gamma_arm, [5.0, 5.0, 10.0, 10.0, 10.0])
    assert "gamma_arm" in caplog_loader.text


def test_a_fixed_value_is_not_combined_with_per_image_values(caplog_loader):
    """A fixed value often means a backend did not find the motor."""
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0])
    second = _Segment([2.0, 3.0], gamma_arm=0.0)

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, "gamma_arm")
    assert "gamma_arm" in caplog_loader.text
    assert any(record.levelno == logging.WARNING for record in caplog_loader.records)


def test_an_arm_missing_from_one_segment_is_not_merged(caplog_loader):
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0])
    second = _Segment([2.0, 3.0])

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, "gamma_arm")
    assert "gamma_arm" in caplog_loader.text


def test_an_arm_value_that_fits_neither_form_is_not_merged(caplog_loader):
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0])
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0, 9.0])

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, "gamma_arm")
    assert "gamma_arm" in caplog_loader.text


@pytest.mark.parametrize(
    "name, value",
    [("arm_angle_frame", "surface"), ("continuous_exposure", True)],
)
def test_arm_properties_shared_by_every_segment_are_taken_over(name, value):
    first = _Segment([0.0, 1.0], **{name: value})
    second = _Segment([2.0, 3.0], **{name: value})

    scan = InterlacedScan([first, second], False, "th")

    assert getattr(scan, name) == value


@pytest.mark.parametrize(
    "name, values",
    [
        ("arm_angle_frame", ("prim", "surface")),
        ("continuous_exposure", (True, False)),
    ],
)
def test_arm_properties_that_disagree_are_reported_and_left_alone(
    name, values, caplog_loader
):
    first = _Segment([0.0, 1.0], **{name: values[0]})
    second = _Segment([2.0, 3.0], **{name: values[1]})

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, name)
    assert name in caplog_loader.text
    assert any(record.levelno == logging.WARNING for record in caplog_loader.records)


def test_an_arm_property_stated_by_only_some_segments_is_reported(caplog_loader):
    first = _Segment([0.0, 1.0], arm_angle_frame="surface")
    second = _Segment([2.0, 3.0])

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, "arm_angle_frame")
    assert "arm_angle_frame" in caplog_loader.text


def test_disagreement_about_the_angle_frame_blocks_the_arm_angles(caplog_loader):
    """The motor values of the segments would be in different conventions."""
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], arm_angle_frame="prim")
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0], arm_angle_frame="surface")

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, "gamma_arm")
    assert not hasattr(scan, "arm_angle_frame")


def test_disagreement_about_a_continuous_exposure_keeps_the_arm_angles():
    """It says how to read an exposure, not where the arm was."""
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], continuous_exposure=True)
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0], continuous_exposure=False)

    scan = InterlacedScan([first, second], False, "th")

    np.testing.assert_allclose(scan.gamma_arm, [1.2, 2.0, 4.0, 8.0])
    assert not hasattr(scan, "continuous_exposure")


def test_the_arm_is_not_reported_as_an_auxiliary_counter():
    """It is geometry, not a counter, wherever a backend lists it."""
    declared = ["gamma_arm", "delta_arm", "ic1"]
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], delta_arm=0.0)
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0], delta_arm=0.0)
    for segment, values in ((first, [10.0, 11.0]), (second, [12.0, 13.0])):
        segment.ic1 = values
        segment.auxillary_counters = list(declared)

    scan = InterlacedScan([first, second], False, "th")

    assert scan.auxillary_counters == ["ic1"]
    np.testing.assert_allclose(scan.gamma_arm, [1.2, 2.0, 4.0, 8.0])


def test_the_merged_arm_reaches_the_angle_conversion():
    """``scan_arm_angles`` is what the pixel conversions ask for the arm."""
    config = SimpleNamespace(
        mu=0.0, chi=0.0, phi=0.0, gamma_arm=0.0, delta_arm=0.0
    )
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], delta_arm=0.0)
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0], delta_arm=0.0)

    scan = InterlacedScan([first, second], False, "th")
    gamma_arm, delta_arm = scan_arm_angles(scan, config)

    # scan motor values are degrees, the conversions take radians
    np.testing.assert_allclose(gamma_arm, np.deg2rad([1.2, 2.0, 4.0, 8.0]))
    np.testing.assert_allclose(delta_arm, np.zeros(4))


def test_exposure_bounds_still_come_from_each_segment_itself():
    """They already descend into the segments; the merge must not disturb it."""
    config = SimpleNamespace(
        mu=0.0, chi=0.0, phi=0.0, gamma_arm=0.0, delta_arm=0.0
    )
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], delta_arm=0.0)
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0], delta_arm=0.0)

    scan = InterlacedScan([first, second], False, "th")
    bounds = scan_exposure_angle_bounds(scan, config)

    np.testing.assert_allclose(
        bounds[:, :, 4], np.deg2rad([[1.2, 1.2], [2.0, 2.0], [4.0, 4.0], [8.0, 8.0]])
    )


def test_shared_arm_geometry_of_no_segments_is_empty():
    assert shared_arm_geometry([]) == {}


def test_a_silent_segment_agrees_with_a_stationary_arm(caplog_loader):
    """Not declaring ``continuous_exposure`` means the arm does not move."""
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0], continuous_exposure=False)
    second = _Segment([2.0, 3.0], gamma_arm=[4.0, 8.0])

    scan = InterlacedScan([first, second], False, "th")

    assert scan.continuous_exposure is False
    np.testing.assert_allclose(scan.gamma_arm, [1.2, 2.0, 4.0, 8.0])
    assert not [
        record for record in caplog_loader.records if record.levelno >= logging.WARNING
    ]


def test_an_arm_value_that_is_not_a_number_does_not_fail_the_scan(caplog_loader):
    first = _Segment([0.0, 1.0], gamma_arm=[1.2, 2.0])
    second = _Segment([2.0, 3.0], gamma_arm="unknown")

    scan = InterlacedScan([first, second], False, "th")

    assert not hasattr(scan, "gamma_arm")
    assert "gamma_arm" in caplog_loader.text

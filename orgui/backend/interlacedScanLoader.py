# /*##########################################################################
#
# Copyright (c) 2024-2025 Finn Schroeter
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
__author__ = "Finn Schroeter"
__credits__ = ["Finn Schroeter"]
__copyright__ = "Copyright 2024-2025 Finn Schroeter"
__license__ = "MIT License"
__version__ = "1.3.0"

import logging

import numpy as np

logger = logging.getLogger(__name__)

# Counters the reconstruction and the integration read from a scan by name,
# whether or not a backend declares them in ``auxillary_counters``. They are
# merged like any declared counter, so that a segmented scan normalizes the
# same way a single scan does.
_IMPLICIT_COUNTERS = ("exposure_time", "exposure_time_variance")

# Attributes the segmented scan sets for itself. A counter of one of these
# names would overwrite the scan geometry, so it is never merged. Names
# defined on the class -- its methods and properties -- are excluded as well,
# see :func:`_counter_candidates`.
_RESERVED_NAMES = frozenset(
    {
        "arm_angle_frame",
        "axis",
        "axisname",
        "continuous_exposure",
        "delta_arm",
        "gamma_arm",
        "indices",
        "mu",
        "nopoints",
        "omega",
        "sort",
        "subscans",
        "th",
        "title",
    }
)

#: Detector arm motor positions merged from the segments. Degrees, like the
#: motor values a backend reports and like ``th``/``mu`` -- only the fixed
#: values in the configuration are radians.
_ARM_ANGLE_NAMES = ("gamma_arm", "delta_arm")

#: Arm properties that say how to read the arm motors, rather than where the
#: arm was. They must agree exactly between the segments.
_ARM_FLAG_NAMES = ("arm_angle_frame", "continuous_exposure")


def _segment_length(scan):
    """Number of images in one scan segment."""
    points = getattr(scan, "nopoints", None)
    if points is None:
        points = len(scan)
    return int(points)


def _axis_column(scan, name):
    """Return one selected-axis value per image of ``scan``.

    Backends store a motor that did not move as either a scalar or a
    one-element array. If that motor is selected as the segmented-scan axis,
    its fixed position applies to every image in the segment and must be
    expanded before concatenation.

    :param scan: One scan segment.
    :param str name: Motor name, currently ``"th"`` or ``"mu"``.
    :returns: One-dimensional floating-point array aligned with the images.
    :rtype: numpy.ndarray
    :raises ValueError: If the motor values cannot be aligned with the images.
    """
    points = _segment_length(scan)
    try:
        values = np.asarray(getattr(scan, name), dtype=np.float64)
    except (AttributeError, TypeError, ValueError) as error:
        raise ValueError(
            f"Scan segment does not provide a numeric '{name}' axis."
        ) from error

    if values.size == 1:
        return np.full(points, float(values.reshape(-1)[0]))
    if values.ndim == 1 and values.size == points:
        return values
    raise ValueError(
        f"Scan segment has {values.size} '{name}' values for {points} images; "
        "the selected axis must be fixed or provide one value per image."
    )


def _counter_column(scan, name):
    """One counter value per image of ``scan``, or ``None`` if unusable.

    A counter held once for a whole segment (a scalar, e.g. an exposure time
    set for the segment as a whole) is expanded to one value per image, since
    the segments of a segmented scan may well have been measured with
    different settings.

    :param scan: One scan segment.
    :param str name: Counter name.
    :returns: Array whose first axis has one entry per image, or ``None`` when
        the segment has no such counter or holds a number of values that
        cannot be aligned with its images.
    :rtype: numpy.ndarray or None
    """
    value = getattr(scan, name, None)
    if value is None:
        return None
    array = np.asarray(value)
    if array.dtype == np.dtype(object):
        return None
    points = _segment_length(scan)
    if array.ndim == 0:
        return np.full(points, array)
    if array.shape[0] == points:
        return array
    if array.shape[0] == 1:
        return np.repeat(array, points, axis=0)
    return None


def _counter_candidates(scansegments):
    """Counter names any segment declares, in the order they are declared.

    :param scansegments: The scan segments.
    :returns: Candidate counter names, without duplicates.
    :rtype: list[str]
    """
    names = []
    for scan in scansegments:
        declared = getattr(scan, "auxillary_counters", ()) or ()
        for name in list(declared) + [f"{n}_variance" for n in declared]:
            if name in names or name in _RESERVED_NAMES:
                continue
            if hasattr(InterlacedScan, name):
                continue
            names.append(name)
    for name in _IMPLICIT_COUNTERS:
        if name not in names:
            names.append(name)
    return names


def shared_counters(scansegments, indices=None):
    """Counters that every segment of a segmented scan provides.

    A segmented ("interlaced") scan concatenates the images of several scans,
    so a counter is only meaningful for the combined scan when *every* segment
    provides it: taking a counter that only some segments have would leave
    fewer values than there are images, and misalign every value after the
    gap. Counters found in only part of the segments are therefore dropped,
    with a log message naming them.

    Backends declare their counters in :attr:`Scan.auxillary_counters` but
    typically set the attribute only where the file actually contains the
    counter, so both the declaration and the value are checked. The counters
    the analysis reads by name -- ``exposure_time`` and its variance -- are
    considered whether or not a backend declares them.

    This merges values, not units or conventions: the segments are assumed to
    have been measured with the same counters meaning the same thing, which is
    what makes them segments of one scan.

    :param scansegments: The scan segments, in the order in which their images
        are concatenated.
    :param indices: Optional index array reordering the concatenated images,
        as used when the segments are sorted by their axis value. Counter
        values are reordered with it so that they stay aligned with the
        images.
    :returns: Counter name to one value per image of the combined scan, in the
        order the segments declare them.
    :rtype: dict
    """
    scansegments = list(scansegments)
    counters = {}
    if not scansegments:
        return counters

    incomplete = []
    for name in _counter_candidates(scansegments):
        columns = [_counter_column(scan, name) for scan in scansegments]
        if any(column is None for column in columns):
            if any(column is not None for column in columns):
                incomplete.append(name)
            continue
        try:
            merged = np.concatenate(columns)
        except ValueError:
            logger.warning(
                "Counter '%s' has a different shape in different scan "
                "segments and is not available in the segmented scan.",
                name,
            )
            continue
        if indices is not None:
            merged = merged[indices]
        counters[name] = merged

    if incomplete:
        logger.info(
            "Counters %s are not available in all scan segments and are "
            "therefore not available in the segmented scan.",
            ", ".join(f"'{name}'" for name in incomplete),
        )
    return counters


def _arm_flag(scan, name):
    """One arm property of a segment, normalized for comparison.

    :param scan: One scan segment.
    :param str name: One of :data:`_ARM_FLAG_NAMES`.
    :returns: The value, or ``None`` when the segment does not state it.
    """
    value = getattr(scan, name, None)
    if value is None:
        return None
    return bool(value) if name == "continuous_exposure" else str(value)


def _merge_arm_flag(scansegments, name):
    """One arm property shared by every segment.

    :param scansegments: The scan segments.
    :param str name: One of :data:`_ARM_FLAG_NAMES`.
    :returns: ``(value, agree)``. ``value`` is ``None`` when the property is
        not merged, either because no segment states it or because the
        segments disagree; ``agree`` is False only in the latter case.
    :rtype: tuple
    """
    values = [_arm_flag(scan, name) for scan in scansegments]
    if all(value is None for value in values):
        return None, True
    if name == "continuous_exposure":
        # A scan that does not declare it has an arm that is stationary within
        # each frame, which is the documented default, so silence is not a
        # disagreement here. ``arm_angle_frame`` has no such default: a scan
        # that does not state it leaves the frame to the configuration.
        values = [False if value is None else value for value in values]
    if any(value is None for value in values):
        logger.warning(
            "'%s' is stated by only %d of the %d scan segments, so the "
            "segments do not agree on it and it is not used for the "
            "segmented scan.",
            name,
            sum(value is not None for value in values),
            len(values),
        )
        return None, False
    if any(value != values[0] for value in values):
        logger.warning(
            "The scan segments disagree about '%s' (%s), so it is not used "
            "for the segmented scan.",
            name,
            ", ".join(sorted({str(value) for value in values})),
        )
        return None, False
    return values[0], True


def _segment_name(scan, position):
    """How one scan segment is named in a message.

    :param scan: One scan segment.
    :param int position: Its index among the segments.
    :returns: The scan name the backend gave it, or its position.
    :rtype: str
    """
    name = getattr(scan, "scanname", None)
    return str(name) if name else f"segment {position + 1}"


def _arm_column(scan, name):
    """One segment's arm position, as ``(kind, values)``.

    :param scan: One scan segment.
    :param str name: One of :data:`_ARM_ANGLE_NAMES`.
    :returns: ``("fixed", float)`` for a single value covering the whole
        segment, ``("frames", array)`` for one value per image, ``("bad",
        None)`` for a value that is neither, and ``(None, None)`` when the
        segment does not provide the motor at all.
    :rtype: tuple
    """
    value = getattr(scan, name, None)
    if value is None:
        return None, None
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError):
        # A motor position that is not a number at all is reported like any
        # other unusable value, rather than failing the whole scan.
        return "bad", None
    if array.size == 1:
        return "fixed", float(array.reshape(-1)[0])
    if array.ndim == 1 and array.size == _segment_length(scan):
        return "frames", array
    return "bad", None


def _named_segments(scansegments, columns, kind):
    """The segments whose arm value is of one kind, named for a message.

    :param scansegments: The scan segments.
    :param columns: Their :func:`_arm_column` results.
    :param str kind: The kind to report.
    :returns: The names, comma separated.
    :rtype: str
    """
    return ", ".join(
        _segment_name(scan, position)
        for position, (scan, (column_kind, _values)) in enumerate(
            zip(scansegments, columns)
        )
        if column_kind == kind
    )


def shared_arm_geometry(scansegments, indices=None):
    """Detector arm position and properties shared by every segment.

    A segmented ("interlaced") scan is one scan for everything downstream, so
    the arm has to be described on the combined scan: without this the
    per-pixel conversions fall back to the fixed arm position from the
    configuration, and every image of a scan that moves the detector is
    converted at the wrong place. Reciprocal-space reconstruction is not
    affected either way -- :func:`.scans.scan_exposure_angle_bounds` already
    descends into the segments and reads each one's own arm.

    The arm is merged conservatively, because getting it wrong yields
    plausible-looking but wrong angles rather than an error:

    * ``gamma_arm`` and ``delta_arm`` must be provided by *every* segment, and
      all in the same form -- either a single fixed value per segment, or one
      value per image. A mix of the two is refused, since a fixed value often
      means that a backend did not find the motor rather than that the arm was
      parked. Fixed values that differ between segments are expanded to one
      value per image, which is where the arm actually was.
    * ``arm_angle_frame`` and ``continuous_exposure`` must match exactly.
      Segments that disagree, or that state one of them only in part, are
      reported and the property is left to the configuration. Disagreement
      about ``arm_angle_frame`` also blocks the arm angles themselves: it
      means the motor values of different segments are in different
      conventions, so concatenating them would mix the two.

    Whatever is not merged is simply not set, which leaves the previous
    behaviour -- the fixed value from the configuration -- in place.

    :param scansegments: The scan segments, in the order in which their images
        are concatenated.
    :param indices: Optional index array reordering the concatenated images,
        as used when the segments are sorted by their axis value.
    :returns: Attribute name to value, with the arm angles in degrees.
    :rtype: dict
    """
    scansegments = list(scansegments)
    geometry = {}
    if not scansegments:
        return geometry

    frames_agree = True
    for name in _ARM_FLAG_NAMES:
        value, agree = _merge_arm_flag(scansegments, name)
        if value is not None:
            geometry[name] = value
        if name == "arm_angle_frame" and not agree:
            frames_agree = False

    for name in _ARM_ANGLE_NAMES:
        columns = [_arm_column(scan, name) for scan in scansegments]
        kinds = {kind for kind, _values in columns}
        if kinds == {None}:
            continue
        if None in kinds:
            logger.warning(
                "'%s' is provided by only %d of the %d scan segments and is "
                "therefore not used for the segmented scan.",
                name,
                sum(kind is not None for kind, _values in columns),
                len(columns),
            )
            continue
        if "bad" in kinds:
            logger.warning(
                "'%s' of %s is neither a single value nor one value per image "
                "of that segment, so it is not used for the segmented scan.",
                name,
                _named_segments(scansegments, columns, "bad"),
            )
            continue
        if kinds == {"fixed", "frames"}:
            logger.warning(
                "'%s' is a single fixed value in %s and one value per image "
                "in the other segments. Such a combination is not merged, so "
                "it is not used for the segmented scan.",
                name,
                _named_segments(scansegments, columns, "fixed"),
            )
            continue
        if not frames_agree:
            logger.warning(
                "'%s' is not used for the segmented scan, because the "
                "segments disagree about 'arm_angle_frame' and their motor "
                "values may therefore not mean the same thing.",
                name,
            )
            continue

        if kinds == {"fixed"}:
            values = [value for _kind, value in columns]
            if all(value == values[0] for value in values):
                geometry[name] = values[0]
                continue
            logger.info(
                "'%s' is a fixed value that differs between the scan "
                "segments (%s), so the segmented scan carries one value per "
                "image.",
                name,
                ", ".join(f"{value:g}" for value in values),
            )
            columns = [
                ("frames", np.full(_segment_length(scan), value))
                for scan, (_kind, value) in zip(scansegments, columns)
            ]

        merged = np.concatenate([values for _kind, values in columns])
        if indices is not None:
            merged = merged[indices]
        geometry[name] = merged

    return geometry


class InterlacedScan:
    def __init__(self, scansegments, sort, ax):
        self.subscans = scansegments

        scan0 = self.subscans[0]
        lensum = 0
        ax_area = []

        self.sort = sort

        # self.axisname = scan0.axisname # todo: check if all scans have the same axis as selected in GUI window  # noqa: E501
        self.axisname = ax

        if ax == "th":
            for i in scansegments:
                lensum += i.nopoints
                ax_area.append(_axis_column(i, "th"))

            self.th = np.concatenate(ax_area)
            self.omega = -1 * self.th

            self.mu = scan0.mu  # todo: check if all scans share the same mu
            self.axis = self.th

        if ax == "mu":
            for i in scansegments:
                lensum += i.nopoints
                ax_area.append(_axis_column(i, "mu"))

            self.mu = np.concatenate(ax_area)

            self.th = scan0.th  # todo: check if all scans share the same th
            self.omega = -1 * self.th
            self.axis = self.mu

        if self.sort:
            if ax == "mu":
                self.indices = np.argsort(self.mu)
                self.mu = self.mu[self.indices]
                self.axis = self.mu
            elif ax == "th":
                self.indices = np.argsort(self.th)
                self.th = self.th[self.indices]
                self.axis = self.th
        else:
            self.indices = None

        self.title = "interlaced scan"

        self.nopoints = lensum

        # Counters are merged here rather than read on demand, so that the
        # combined scan carries them as plain attributes, exactly as a scan
        # from a backend does.
        self._shared_counters = shared_counters(self.subscans, self.indices)
        for name, values in self._shared_counters.items():
            setattr(self, name, values)

        # The detector arm is geometry rather than a counter, so it is merged
        # separately and is not reported as an auxiliary counter.
        arm = shared_arm_geometry(self.subscans, self.indices)
        for name, value in arm.items():
            setattr(self, name, value)

    def get_raw_img(self, i):
        len_previous = 0
        if self.sort:
            for k in self.subscans:
                if self.indices[i] < k.nopoints + len_previous:
                    return k.get_raw_img(self.indices[i] - len_previous)
                else:
                    len_previous += k.nopoints
        else:
            for k in self.subscans:
                if i < k.nopoints + len_previous:
                    return k.get_raw_img(i - len_previous)
                else:
                    len_previous += k.nopoints

    def get_origin_scan_nr(self, i):
        len_previous = 0
        if self.sort:
            for k in self.subscans:
                if self.indices[i] < k.nopoints + len_previous:
                    try:
                        return k.scanname
                    except Exception:
                        return
                else:
                    len_previous += k.nopoints
        else:
            for k in self.subscans:
                if i < k.nopoints + len_previous:
                    try:
                        return k.scanname
                    except Exception:
                        return
                else:
                    len_previous += k.nopoints

    @property
    def auxillary_counters(self):
        """Counters every scan segment provides, merged over the segments.

        A counter that only some of the segments have is not listed, see
        :func:`shared_counters`. Each listed counter is also set as an
        attribute holding one value per image of the combined scan.
        """
        return list(self._shared_counters)

    def __len__(self):
        return self.nopoints

    @classmethod
    def parse_h5_node(cls, obj):
        ddict = dict()
        return ddict

"""The segmented scan loader asks the backend, it does not guess.

Two routes lead to the same answer: a backend implementing ``listScans`` is
asked directly, and any other backend has its own ``parse_h5_node`` applied to
every entry of the file root. Neither route may disagree with how a single scan
is opened, because the identifiers produced here are handed to ``openScan``.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from orgui.backend import scans

orGUI = pytest.importorskip("orgui.app.orGUI", reason="Qt bindings are required")


# --- fast path -----------------------------------------------------------


class _ListingBackend:
    """A backend that implements the fast path."""

    returns = []

    @classmethod
    def listScans(cls, h5file):
        del h5file
        return cls.returns


def _fastPath(returns):
    backend = type("_Backend", (_ListingBackend,), {"returns": returns})
    return orGUI.orGUI._scanEntriesFromBackend(backend, object())


def test_bare_numbers_are_accepted():
    """The simplest thing a backend can return."""
    entries = _fastPath([1, 2, 10])

    assert [ddict["scanno"] for ddict, _ in entries] == [1, 2, 10]
    assert [label for _, label in entries] == [None, None, None]


def test_bare_names_are_accepted():
    """A backend may be addressed by a string rather than by a number."""
    entries = _fastPath(["ascan_12", "dscan_3"])

    assert [ddict["scanno"] for ddict, _ in entries] == ["ascan_12", "dscan_3"]


def test_an_identifier_is_reported_under_both_keys():
    """openScan reads 'name' for some beamtimes and 'scanno' for others."""
    ((ddict, _label),) = _fastPath(["ascan_12"])

    assert ddict["scanno"] == "ascan_12"
    assert ddict["name"] == "ascan_12"


def test_pairs_carry_a_label_for_the_dialog():
    entries = _fastPath([(1, "ascan th 0 90 90 1"), (2, "dscan mu 0 1 20 1")])

    assert [ddict["scanno"] for ddict, _ in entries] == [1, 2]
    assert [label for _, label in entries] == [
        "ascan th 0 90 90 1",
        "dscan mu 0 1 20 1",
    ]


def test_the_legacy_id31_listing_fills_the_key_openscan_reads():
    """openScan opens ch5523 with ddict['name'], not with a scan number.

    The segmented loader never set that key, so this beamtime raised KeyError
    before it could open anything. Its listing reports group names, which the
    normalisation has to put where openScan looks.
    """
    from orgui.backend.beamline import id31_tools

    class _LegacyFile:
        names = {"ascan_12": b"ascan th 0 90 90 1"}

        def __iter__(self):
            return iter(self.names)

        def __getitem__(self, path):
            name, _, field = path.partition("/")
            if name not in self.names or field != "title":
                raise KeyError(path)
            return type("_D", (), {"__getitem__": lambda s, k: self.names[name]})()

    listed = id31_tools.BlissScan.listScans(_LegacyFile())
    backend = type("_Backend", (_ListingBackend,), {"returns": listed})

    ((ddict, label),) = orGUI.orGUI._scanEntriesFromBackend(backend, object())

    assert ddict["name"] == "ascan_12"
    assert label == "ascan th 0 90 90 1"


def test_a_backend_without_the_fast_path_asks_for_the_fallback():
    """None is the signal to fall back, and must not be confused with []."""

    class _NoFastPath(scans.Scan):
        def __init__(self):
            pass

        @classmethod
        def parse_h5_node(cls, node):
            return {}

        def __len__(self):
            return 0

        def get_raw_img(self, i):
            return None

    assert orGUI.orGUI._scanEntriesFromBackend(_NoFastPath, object()) is None
    assert _fastPath([]) == []


# --- fallback ------------------------------------------------------------


class _StubItem:
    def __init__(self, name):
        self.name = name


class _StubModel:
    """Enough of the silx HDF5 tree model to walk the root's children."""

    def __init__(self, names):
        self._names = list(names)

    def rowCount(self, parent):
        del parent
        return len(self._names)

    def index(self, row, column, parent):
        del column, parent
        return row

    def data(self, index, role):
        del role
        return _StubItem(self._names[index])


class _StubNode:
    """What silx H5Node reports for a direct child of the file root."""

    def __init__(self, item):
        self.basename = item.name
        self.local_name = "/" + item.name


@pytest.fixture
def h5node(monkeypatch):
    """Wrap tree items in the stub node instead of the real silx one."""
    monkeypatch.setattr(orGUI.silx.gui.hdf5, "H5Node", _StubNode)


def _fallback(backend, names):
    return orGUI.orGUI._scanEntriesFromNodes(backend, _StubModel(names), None)


class _BlissLikeBackend:
    """Numbered by the integer in front of the dot, as the ID31 backends are."""

    @classmethod
    def parse_h5_node(cls, node):
        scanno, _subscan = node.basename.split(".")
        return {"scanno": int(scanno), "name": node.local_name}


class _NamedBackend:
    """Addressed by name, as the legacy ch5523 backend is."""

    @classmethod
    def parse_h5_node(cls, node):
        return {
            "scanno": int(node.local_name.split("_")[-1]),
            "name": node.local_name.strip("/"),
        }


def test_fallback_numbers_scans_the_way_the_backend_does(h5node):
    entries = _fallback(_BlissLikeBackend, ["1.1", "2.1", "10.1"])

    assert [ddict["scanno"] for ddict, _ in entries] == [1, 2, 10]


def test_fallback_collapses_the_subscans_of_one_measurement(h5node):
    """1.1 and 1.2 are the same scan, and must appear once."""
    entries = _fallback(_BlissLikeBackend, ["1.1", "1.2", "2.1", "2.2"])

    assert [ddict["scanno"] for ddict, _ in entries] == [1, 2]
    # the first occurrence wins, which is the fast counter subscan
    assert entries[0][0]["name"] == "/1.1"


@pytest.mark.parametrize("name", ["instrument", "title", "start_time", "entry"])
def test_fallback_skips_entries_the_backend_cannot_parse(h5node, name):
    """Raising from parse_h5_node is how a backend says "not a scan"."""
    entries = _fallback(_BlissLikeBackend, [name, "3.1"])

    assert [ddict["scanno"] for ddict, _ in entries] == [3]


def test_fallback_skips_a_backend_that_reports_no_scan_number(h5node):
    """InterlacedScan and ImportImagesScan return an empty dict."""

    class _NoScanNo:
        @classmethod
        def parse_h5_node(cls, node):
            return {}

    assert _fallback(_NoScanNo, ["1.1", "2.1"]) == []


def test_fallback_keeps_the_name_a_backend_is_opened_with(h5node):
    """ch5523 is opened by name; that name has to survive to openScan."""
    entries = _fallback(_NamedBackend, ["ascan_12", "dscan_3"])

    assert [ddict["name"] for ddict, _ in entries] == ["ascan_12", "dscan_3"]
    assert [ddict["scanno"] for ddict, _ in entries] == [12, 3]


def test_fallback_reports_no_label(h5node):
    """parse_h5_node has no title to offer, so the dialog shows the name."""
    entries = _fallback(_BlissLikeBackend, ["1.1"])

    assert entries[0][1] is None


def test_both_routes_agree_on_the_same_file(h5node):
    """The whole point: the fast path may not disagree with parse_h5_node."""
    names = ["1.1", "1.2", "2.1", "2.2", "10.1"]

    fallback = _fallback(_BlissLikeBackend, names)
    fast = _fastPath([1, 2, 10])

    assert [ddict["scanno"] for ddict, _ in fallback] == [
        ddict["scanno"] for ddict, _ in fast
    ]


# --- selected axis -------------------------------------------------------


def test_all_scalar_selected_axis_shows_a_warning(monkeypatch):
    segments = [
        SimpleNamespace(th=1.0, mu=[0.0, 1.0]),
        SimpleNamespace(th=np.array(2.0), mu=[2.0, 3.0]),
    ]
    warnings = []
    monkeypatch.setattr(
        orGUI.qutils,
        "warning_detailed_message",
        lambda *args: warnings.append(args),
    )
    parent = SimpleNamespace(
        _scanAxisHasOnlyScalars=orGUI.orGUI._scanAxisHasOnlyScalars
    )

    result = orGUI.orGUI._warnForScalarScanAxis(parent, segments, "th")

    assert result is None
    assert len(warnings) == 1
    assert warnings[0][0] is parent
    assert "'th'" in warnings[0][2]
    assert "Select 'mu'" in warnings[0][2]


def test_a_per_image_selected_axis_does_not_show_a_warning(monkeypatch):
    segments = [
        SimpleNamespace(th=1.0, mu=[0.0, 1.0]),
        SimpleNamespace(th=[2.0, 3.0], mu=[2.0, 3.0]),
    ]
    warnings = []
    monkeypatch.setattr(
        orGUI.qutils,
        "warning_detailed_message",
        lambda *args: warnings.append(args),
    )

    parent = SimpleNamespace(
        _scanAxisHasOnlyScalars=orGUI.orGUI._scanAxisHasOnlyScalars
    )
    result = orGUI.orGUI._warnForScalarScanAxis(parent, segments, "th")

    assert result is None
    assert warnings == []

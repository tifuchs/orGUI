"""Backends describe their own file layout to the segmented scan loader.

The loader used to guess the layout from the backend's class name, which
silently mishandled every backend not called ``BlissScan*`` -- including the
example backend shipped with orGUI. A backend now either implements
``Scan.listScans``, or lets the loader fall back to its own ``parse_h5_node``.
"""

import importlib.util
import pathlib

import pytest

from orgui.backend import scans
from orgui.backend.beamline import id31_tools


class _StubDataset:
    def __init__(self, value):
        self._value = value

    def __getitem__(self, key):
        assert key == ()
        return self._value


class _StubH5File:
    """Just enough of an h5py file: iteration over names and item access."""

    def __init__(self, entries):
        self._entries = dict(entries)

    def __iter__(self):
        return iter(self._entries)

    def __getitem__(self, path):
        name, _, field = path.partition("/")
        if name not in self._entries:
            raise KeyError(path)
        title = self._entries[name]
        if field != "title" or title is None:
            raise KeyError(path)
        return _StubDataset(title)


BLISS_FILE = _StubH5File(
    [
        ("2.1", b"ascan th 0 90 90 1"),
        ("2.2", b"subscan of 2"),
        ("1.1", b"ascan th 0 45 45 1"),
        ("1.2", b"subscan of 1"),
        ("10.1", b"ascan th 0 10 10 1"),
    ]
)

BLISS_BACKENDS = (id31_tools.BlissScan_EBS, id31_tools.BlissScan_EBS_p4)


def test_the_base_class_declares_the_fast_path_as_unavailable():
    """The template declares listScans; not implementing it is legal."""
    with pytest.raises(NotImplementedError):
        scans.Scan.listScans(BLISS_FILE)


def test_the_declaration_is_not_abstract():
    """An @abstractmethod would break every existing backend on construction."""
    assert "listScans" not in scans.Scan.__abstractmethods__


@pytest.mark.parametrize("backend", BLISS_BACKENDS)
def test_only_the_fast_counter_subscan_is_listed(backend):
    entries = backend.listScans(BLISS_FILE)

    assert [scanno for scanno, _ in entries] == [1, 2, 10]


@pytest.mark.parametrize("backend", BLISS_BACKENDS)
def test_scans_are_listed_in_scan_order(backend):
    """The dialog lists them by scan number, not in file order."""
    entries = backend.listScans(BLISS_FILE)

    assert entries == sorted(entries)


@pytest.mark.parametrize("backend", BLISS_BACKENDS)
def test_the_label_is_the_scan_title(backend):
    entries = dict(backend.listScans(BLISS_FILE))

    assert entries[1] == "ascan th 0 45 45 1"
    assert entries[2] == "ascan th 0 90 90 1"


@pytest.mark.parametrize("backend", BLISS_BACKENDS)
def test_a_double_digit_subscan_is_not_the_fast_counter_subscan(backend):
    """``1.10`` must not be mistaken for ``1.1``."""
    entries = backend.listScans(
        _StubH5File([("1.1", b"real"), ("1.10", b"other subscan")])
    )

    assert [scanno for scanno, _ in entries] == [1]


@pytest.mark.parametrize("backend", BLISS_BACKENDS)
@pytest.mark.parametrize("name", ["notascan", "instrument", "1.x", "x.1"])
def test_entries_that_are_not_scans_are_skipped(backend, name):
    entries = backend.listScans(_StubH5File([(name, b"title"), ("3.1", b"ok")]))

    assert [scanno for scanno, _ in entries] == [3]


@pytest.mark.parametrize("backend", BLISS_BACKENDS)
def test_a_scan_without_a_title_falls_back_to_its_name(backend):
    entries = backend.listScans(_StubH5File([("7.1", None)]))

    assert entries == [(7, "7.1")]


# --- the legacy ch5523 format, named "<command>_<scanno>" ----------------

LEGACY_FILE = _StubH5File(
    [
        ("dscan_3", b"dscan mu 0 1 20 1"),
        ("ascan_12", b"ascan th 0 90 90 1"),
        ("ascan_7", b"ascan th 0 45 45 1"),
        ("instrument", b"not a scan"),
    ]
)


def test_the_legacy_backend_lists_its_own_naming_convention():
    """It has no subscans, so Fastscan's implementation does not apply."""
    entries = id31_tools.BlissScan.listScans(LEGACY_FILE)

    assert [name for name, _ in entries] == ["dscan_3", "ascan_7", "ascan_12"]


def test_the_legacy_backend_reports_the_name_it_is_opened_with():
    """openScan passes ddict['name'] to BlissScan, not the scan number.

    Reporting the number would make orGUI look up a group that does not exist.
    """
    entries = dict(id31_tools.BlissScan.listScans(LEGACY_FILE))

    assert "ascan_12" in entries
    assert entries["ascan_12"] == "ascan th 0 90 90 1"


def test_the_legacy_backend_ignores_entries_that_are_not_scans():
    entries = id31_tools.BlissScan.listScans(LEGACY_FILE)

    assert all(name != "instrument" for name, _ in entries)


def test_the_legacy_backend_does_not_list_subscan_style_names():
    """A "1.1" name is not a scan of this format and must not be offered."""
    assert id31_tools.BlissScan.listScans(BLISS_FILE) == []


def _example_backend(filename, classname):
    """Load a shipped example backend, or skip when it is not in the tree."""
    path = (
        pathlib.Path(__file__).resolve().parents[3] / "examples" / "backend" / filename
    )
    if not path.is_file():
        pytest.skip("examples are not part of an installed orGUI")
    spec = importlib.util.spec_from_file_location(path.stem + "_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return getattr(module, classname)


def test_the_example_backend_implements_the_fast_path():
    """The reviewed case: a backend with no "BlissScan" in its name."""
    backend = _example_backend("ID31_EBS_p4_backend.py", "ID31_EBS_p4_2026")

    assert "BlissScan" not in backend.__name__
    assert backend.listScans(BLISS_FILE) == id31_tools.BlissScan_EBS_p4.listScans(
        BLISS_FILE
    )


# --- the CHESS QM2 example, SPEC style names with one scan per group ------


def _qm2():
    return _example_backend("CHESS_QM2.py", "QM2_backend_2026")


def _qm2_resolves(h5file, scanno):
    """The group QM2_backend_2026.__init__ would open for this identifier."""
    for name in h5file:
        head, _, tail = str(name).partition(".")
        if int(head) == scanno:
            return ".".join((head, tail))
    return None


QM2_FILE = _StubH5File(
    [
        ("2.1", b"ascan th 0 90 90 1"),
        ("1.1", b"ascan th 0 45 45 1"),
        ("10.1", b"dscan mu 0 1 20 1"),
    ]
)


def test_qm2_lists_the_groups_of_the_file_root():
    entries = _qm2().listScans(QM2_FILE)

    assert [scanno for scanno, _ in entries] == [1, 2, 10]
    assert dict(entries)[2] == "ascan th 0 90 90 1"


def test_qm2_reports_the_scan_number_its_constructor_takes():
    """__init__ matches int(name before the dot), so the number is the id."""
    backend = _qm2()

    for scanno, _title in backend.listScans(QM2_FILE):
        assert _qm2_resolves(QM2_FILE, scanno) is not None


def test_qm2_offers_a_scan_number_once_per_file():
    """__init__ opens the first group with a matching number.

    Offering the later subscans would let one be selected and another opened.
    """
    file = _StubH5File([("1.1", b"first"), ("1.2", b"second"), ("2.1", b"other")])

    entries = _qm2().listScans(file)

    assert [scanno for scanno, _ in entries] == [1, 2]
    assert dict(entries)[1] == "first"


def test_qm2_skips_groups_that_are_not_scans():
    """__init__ would raise on these, so they must not be offered."""
    file = _StubH5File(
        [("1.1", b"real scan"), ("instrument", b"meta"), ("title", None)]
    )

    assert [scanno for scanno, _ in _qm2().listScans(file)] == [1]

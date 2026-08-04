"""Tests for the resilience of the orGUI database against I/O errors.

The database file can become unavailable at any time, e.g. if the drive
holding it is disconnected (issue #23). orGUI must report this and keep
operating instead of terminating.
"""

import os

import h5py
import pytest
from silx.gui import qt

from orgui.app.database import DataBase, DBUnavailableError
from orgui.app.qutils import RobustHdf5TreeModel


_QAPP = None


def _getQApplication():
    """Return the application instance, keeping a reference to it."""
    global _QAPP
    _QAPP = qt.QApplication.instance() or qt.QApplication([])
    return _QAPP


def _writeError(*args, **kwargs):
    raise RuntimeError(
        "Can't decrement id ref count (file write failed: errno = 22, "
        "error message = 'Invalid argument')"
    )


@pytest.fixture
def database(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _getQApplication()
    db = DataBase(None)
    yield db
    db.closeSafe()
    # the widget is not deleted on purpose: silx parents its global wait icon
    # and thread pool to the first widget using them, deleting the widget
    # breaks the tree model of all following tests.


def test_temporary_database_is_open_after_creation(database):
    assert database.isOpen()
    assert os.path.isfile(database._filepath)


def test_close_does_not_keep_a_stale_file_handle(database):
    database.close()

    assert not database.isOpen()
    with pytest.raises(DBUnavailableError):
        database.add_nxdict({"@NX_class": "NXcollection", "counter": 1})


def test_failed_close_detaches_and_allows_a_new_database(database, tmp_path):
    nxfile = database.nxfile
    nxfile.close = _writeError  # simulate a disconnected drive

    with pytest.raises(Exception):
        database.close()

    # orGUI must not hold on to the broken file, otherwise no new database
    # can be created afterwards.
    assert not database.isOpen()
    assert database.hdf5model.rowCount() == 0

    database.createNewDBFile(str(tmp_path / "new_database.h5"))
    assert database.isOpen()
    database.add_nxdict({"@NX_class": "NXcollection", "counter": 1}, h5path="/entry")
    assert "entry/counter" in database.nxfile


def test_close_safe_reports_failure_without_raising(database):
    database.nxfile.close = _writeError

    assert database.closeSafe() is False
    assert not database.isOpen()
    assert database.closeSafe() is True  # nothing left to close


def test_fallback_database_is_created_if_the_new_database_fails(database, monkeypatch):
    monkeypatch.setattr(
        qt.QFileDialog,
        "getSaveFileName",
        staticmethod(
            lambda *args, **kwargs: (
                os.path.join(os.getcwd(), "missing_drive", "database.h5"),
                "",
            )
        ),
    )

    database.onNewDatabase()  # must not raise, the directory does not exist

    # a temporary database keeps orGUI operational
    assert database.isOpen()
    database.add_nxdict({"@NX_class": "NXcollection", "counter": 1}, h5path="/entry")


def test_robust_tree_model_removes_file_which_cannot_be_closed(tmp_path):
    _getQApplication()
    filename = str(tmp_path / "data.h5")
    with h5py.File(filename, "w") as h5:
        h5["counter"] = 1

    model = RobustHdf5TreeModel(ownFiles=True)
    model.insertFile(filename)
    assert model.rowCount() == 1

    h5file = model.data(
        model.index(0, 0), role=RobustHdf5TreeModel.H5PY_OBJECT_ROLE
    ).file
    h5file.close = _writeError

    model.removeH5pyObject(h5file)  # must not raise into the Qt event loop
    assert model.rowCount() == 0

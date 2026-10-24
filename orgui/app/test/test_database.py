"""Tests for the resilience of the orGUI database against I/O errors.

The database file can become unavailable at any time, e.g. if the drive
holding it is disconnected (issue #23). orGUI must report this and keep
operating instead of terminating.
"""

import os
import sys
import threading
from types import SimpleNamespace
from unittest.mock import Mock

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


def test_failed_integration_flush_detaches_unusable_database(database, tmp_path):
    nxfile = database.nxfile
    nxfile.flush = _writeError  # buffered write fails after dicttonx returns
    nxfile.close = _writeError  # disconnected drive also fails on close

    with pytest.raises(RuntimeError, match="file write failed"):
        database.add_nxdict({"entry": {"counter": 1}})

    assert not database.isOpen()
    assert database.hdf5model.rowCount() == 0

    database.createNewDBFile(str(tmp_path / "recovered.h5"))
    database.add_nxdict({"entry": {"counter": 2}})
    assert database.nxfile["entry/counter"][()] == 2


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


@pytest.fixture
def exception_hook(database, monkeypatch):
    """Install the production hook with dialogs captured instead of opened."""
    from orgui.app.orGUI import UncaughtHook
    import importlib

    module = importlib.import_module("orgui.app.orGUI")
    # Record the original hook before the constructor changes it.
    monkeypatch.setattr(sys, "excepthook", sys.excepthook)
    hook = UncaughtHook()
    hook.set_orgui(SimpleNamespace(database=database))
    monkeypatch.setattr(module.logger, "warning", Mock())
    monkeypatch.setattr(sys, "exit", Mock())
    yield hook, module.logger.warning
    qt.QApplication.processEvents()


def _hdf5_error(filename, error_type=OSError):
    # Reproduce Cython's traceback, which enters sys.excepthook directly
    # instead of raising through the application's try/except boundary.
    error = error_type(
        f"Write through page buffer failed (file write failed: "
        f"filename = '{filename}', errno = 22, error message = 'Invalid argument')"
    )
    try:
        exec(compile("raise error", "h5py/_objects.pyx", "exec"), {"error": error})
    except error_type:
        return sys.exc_info()


@pytest.mark.parametrize("error_type", [OSError, RuntimeError])
def test_hdf5_finalizer_failure_recovers_without_exit(
    database, exception_hook, tmp_path, error_type
):
    """HDF5 finalizer errors disable writes now and close on the GUI thread."""
    hook, warning = exception_hook
    old_file = database.nxfile
    error = _hdf5_error(database._filepath, error_type)
    hook.exception_hook(*error)
    hook.exception_hook(*error)  # one failure can produce several callbacks

    assert not database.isOpen()
    assert database.nxfile is old_file  # no HDF5 re-entry from the hook
    assert old_file.id.valid
    with pytest.raises(DBUnavailableError, match="Create or open"):
        database.add_nxdict({"entry": {"counter": 1}})
    assert not warning.called

    qt.QApplication.processEvents()
    assert database.nxfile is None
    assert database.hdf5model.rowCount() == 0
    assert not old_file.id.valid
    warning.assert_called_once()
    assert "create or open" in warning.call_args.kwargs["extra"]["description"]
    sys.exit.assert_not_called()

    database.createNewDBFile(str(tmp_path / "replacement.h5"))
    database.add_nxdict({"entry": {"counter": 2}})
    assert database.nxfile["entry/counter"][()] == 2


def test_recovery_close_error_does_not_reenter_fatal_hook(
    database, exception_hook, monkeypatch, tmp_path
):
    """A second finalizer error during failed close still permits reopening."""
    hook, warning = exception_hook
    old_file = database.nxfile
    original_close = old_file.close
    error = _hdf5_error(database._filepath)

    def failed_close():
        hook.exception_hook(*error)
        raise error[1]

    monkeypatch.setattr(old_file, "close", failed_close)
    try:
        hook.exception_hook(*error)
        qt.QApplication.processEvents()
        assert database.nxfile is None
        assert database.hdf5model.rowCount() == 0
        warning.assert_called_once()
        sys.exit.assert_not_called()
        database.createNewDBFile(str(tmp_path / "reconnected.h5"))
        database.add_nxdict({"entry": {"counter": 3}})
        assert database.nxfile["entry/counter"][()] == 3
    finally:
        original_close()


def test_queued_recovery_does_not_close_replacement(
    database, exception_hook, tmp_path
):
    """Delayed recovery and finalizers must not invalidate a new database."""
    hook, _ = exception_hook
    error = _hdf5_error(database._filepath)
    hook.exception_hook(*error)
    database.createNewDBFile(str(tmp_path / "replacement.h5"))
    replacement = database.nxfile
    qt.QApplication.processEvents()
    hook.exception_hook(*error)  # late finalizer from the previous database
    qt.QApplication.processEvents()
    assert database.nxfile is replacement
    assert database.isOpen()
    database.add_nxdict({"entry": {"counter": 4}})
    sys.exit.assert_not_called()


def test_missing_database_is_nonfatal(database, exception_hook):
    """A GUI action with no open database requests a replacement."""
    hook, warning = exception_hook
    database.close()
    try:
        database._requireOpenFile()
    except DBUnavailableError:
        hook.exception_hook(*sys.exc_info())
    qt.QApplication.processEvents()
    warning.assert_called_once()
    sys.exit.assert_not_called()


def test_database_recovery_does_not_swallow_unrelated_errors(
    database, exception_hook
):
    """Only identified HDF5 storage failures qualify for database recovery."""
    hook, _ = exception_hook
    error = _hdf5_error("another_file.h5")
    assert not hook._is_database_io_error(error[1], error[2])
    assert not hook._is_database_io_error(ValueError("bad shape"), error[2])
    assert not hook._is_database_io_error(
        OSError(f"file write failed: {database._filepath}"), None
    )
    assert database.isOpen()


def test_finalizer_error_during_flush_does_not_report_success(
    database, exception_hook, monkeypatch
):
    """An unraisable HDF5 error must also abort the enclosing integration save."""
    hook, _ = exception_hook
    error = _hdf5_error(database._filepath)
    monkeypatch.setattr(database.nxfile, "flush", lambda: hook.exception_hook(*error))
    with pytest.raises(DBUnavailableError):
        database.add_nxdict({"entry": {"counter": 1}})
    qt.QApplication.processEvents()
    assert database.nxfile is None
    sys.exit.assert_not_called()


def test_recovery_from_worker_closes_on_gui_thread(
    database, exception_hook, monkeypatch
):
    """A native callback on a worker must not close HDF5 or touch widgets."""
    hook, _ = exception_hook
    error = _hdf5_error(database._filepath)
    close = database.closeSafe
    close_threads = []

    def record_close():
        close_threads.append(threading.get_ident())
        return close()

    monkeypatch.setattr(database, "closeSafe", record_close)
    worker = threading.Thread(target=hook.exception_hook, args=error)
    worker.start()
    worker.join(timeout=5)
    assert not worker.is_alive()
    assert not close_threads
    qt.QApplication.processEvents()
    assert close_threads == [threading.get_ident()]
    sys.exit.assert_not_called()


def test_recovery_during_close_does_not_close_twice(
    database, exception_hook, monkeypatch
):
    """Processing Qt events during shutdown must not re-enter file closure."""
    hook, _ = exception_hook
    error = _hdf5_error(database._filepath)
    close = database.nxfile.close
    tracked_close = Mock(wraps=close)
    monkeypatch.setattr(database.nxfile, "close", tracked_close)
    monkeypatch.setattr(
        database, "_waitForPendingOperations", qt.QApplication.processEvents
    )
    hook.exception_hook(*error)
    database.close()
    tracked_close.assert_called_once()
    assert not database._closing
    sys.exit.assert_not_called()


def test_failed_integration_restores_controls(database, exception_hook):
    """An aborted job releases its busy flag and re-enables the GUI controls."""
    from orgui.app.orGUI import orGUI

    hook, _ = exception_hook
    controls = [qt.QWidget(), qt.QWidget()]
    owner = SimpleNamespace(
        centralWidget=lambda: controls[0],
        menuBar=lambda: controls[1],
        _integration_context=lambda: None,
    )
    error = _hdf5_error(database._filepath)

    def operation(context):
        assert not any(widget.isEnabled() for widget in controls)
        hook.exception_hook(*error)
        database._requireOpenFile()

    with pytest.raises(DBUnavailableError):
        orGUI._run_scan_integration(owner, operation)
    assert all(widget.isEnabled() for widget in controls)
    assert not owner._scan_integration_active
    qt.QApplication.processEvents()
    sys.exit.assert_not_called()


def test_open_existing_database_after_recovery(
    database, exception_hook, tmp_path
):
    """Recovery supports opening an existing file as well as creating one."""
    filename = tmp_path / "existing.h5"
    with h5py.File(filename, "w") as file:
        file["counter"] = 5
    hook, _ = exception_hook
    hook.exception_hook(*_hdf5_error(database._filepath))
    qt.QApplication.processEvents()
    database.openDBFile(str(filename))
    assert database.isOpen()
    assert database.nxfile["counter"][()] == 5


def test_programming_error_keeps_existing_fatal_policy(
    database, exception_hook, monkeypatch
):
    """The recovery branch does not silently swallow unrelated exceptions."""
    from orgui.app import qutils

    hook, _ = exception_hook
    monkeypatch.setattr(
        qutils, "critical_detailed_message", lambda *a: qt.QMessageBox.Discard
    )
    hook.exception_hook(ValueError, ValueError("unexpected bug"), None)
    sys.exit.assert_called_once_with(1)

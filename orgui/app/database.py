# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
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
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import silx.gui.hdf5
from silx.gui import icons
from silx.gui.data import DataViewerFrame
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
from silx.gui import qt
import tempfile
import silx.io.h5py_utils
from silx.io import nxdata
import silx.io.utils
import h5py

import datetime
import json
import os
import traceback
import time

from .. import resources
from .config_data import ConfigData, ConfigHandler

import logging
from silx.io.dictdump import dicttonx, nxtodict

logger = logging.getLogger(__name__)


def config_data_to_json(config):
    """Convert :class:`ConfigData` to the central JSON representation."""
    if not isinstance(config, ConfigData):
        raise TypeError("config must be a ConfigData instance")
    return config.to_json_dict()


def config_data_from_json(values):
    """Create :class:`ConfigData` from the central JSON representation."""
    return ConfigData.from_json_dict(values)


class DBCloseError(IOError):
    pass


class DBUnavailableError(IOError):
    """Raised if an operation requires a database file, but none is open."""

    pass


DEFAULT_FILTERS = {  # Filters available with h5py/libhdf5
    "Raw": None,
    "GZip": "gzip",
    "LZF": "lzf",
}

FILTERS = {**DEFAULT_FILTERS}

try:
    import hdf5plugin

    LOSSLESS_FILTERS = {
        "BZip2": hdf5plugin.BZip2(),
        "LZ4": hdf5plugin.LZ4(),
        "ZStd": hdf5plugin.Zstd(),
    }
    FILTERS.update(**LOSSLESS_FILTERS)

    BITSHUFFLE_FILTERS = {
        "Bitshuffle-lz4": hdf5plugin.Bitshuffle(cname="lz4"),
        "Bitshuffle-zstd": hdf5plugin.Bitshuffle(cname="zstd"),
    }
    FILTERS.update(**BITSHUFFLE_FILTERS)

    BLOSC_FILTERS = {}
    for cname in ("lz4", "blosclz", "lz4", "lz4hc", "snappy", "zlib", "zstd"):
        for shuffle_name, shuffle in [
            ("NoShuffle", hdf5plugin.Blosc.NOSHUFFLE),
            ("Shuffle", hdf5plugin.Blosc.SHUFFLE),
            ("BitShuffle", hdf5plugin.Blosc.BITSHUFFLE),
        ]:
            for clevel in [5]:  # (1, 3, 5, 9):
                BLOSC_FILTERS[f"Blosc-{cname}-{shuffle_name}-{clevel}"] = (
                    hdf5plugin.Blosc(cname=cname, clevel=clevel, shuffle=shuffle)
                )
    FILTERS.update(**BLOSC_FILTERS)

    BLOSC2_FILTERS = {}
    for cname in ("lz4", "blosclz", "lz4", "lz4hc", "zlib", "zstd"):
        for filters_name, filters in [
            ("NoFilter", hdf5plugin.Blosc2.NOFILTER),
            ("Shuffle", hdf5plugin.Blosc2.SHUFFLE),
            ("BitShuffle", hdf5plugin.Blosc2.BITSHUFFLE),
        ]:
            for clevel in [5]:  # (1, 3, 5, 9):
                BLOSC2_FILTERS[f"Blosc2-{cname}-{filters_name}-{clevel}"] = (
                    hdf5plugin.Blosc2(cname=cname, clevel=clevel, filters=filters)
                )
    FILTERS.update(**BLOSC2_FILTERS)


except Exception:
    logger.warning(
        "Unable to load hdf5plugin compression filters",
        exc_info=True,
        extra={"show_dialog": False},
    )


class DataBase(qt.QMainWindow):
    compression = FILTERS["Raw"]

    sigChangeRockingScan = qt.Signal(object)

    def __init__(self, plot, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.nxfile = None
        self.filedialogdir = os.getcwd()
        self.plot = plot
        self.config_target = parent
        self.config_handler = ConfigHandler(
            parent, create_dataset_args={"compression": self.compression}
        )

        self.view = silx.gui.hdf5.Hdf5TreeView()
        self.view.setSortingEnabled(True)

        # self.view.addContextMenuCallback(self.nexus_treeview_callback)
        # self.hdf5model = self.view.findHdf5TreeModel()

        customModel = Hdf5TreeModel
        customModel.NAME_COLUMN = 0
        customModel.DESCRIPTION_COLUMN = 1
        customModel.VALUE_COLUMN = 2
        customModel.SHAPE_COLUMN = 3
        customModel.TYPE_COLUMN = 4
        customModel.NODE_COLUMN = 5
        customModel.LINK_COLUMN = 6

        self.hdf5model = customModel(self.view, ownFiles=False)
        self.view.setModel(self.hdf5model)

        self.dataviewer = DataViewerFrame.DataViewerFrame()
        self.dataviewerDialog = qt.QDialog(self)
        dvlayout = qt.QVBoxLayout()
        dvlayout.addWidget(self.dataviewer)
        self.dataviewerDialog.setLayout(dvlayout)
        self.dataviewerDialog.setModal(False)

        self.view.setExpandsOnDoubleClick(False)
        self.hdf5model.setFileMoveEnabled(True)
        # self.__treeModelSorted = silx.gui.hdf5.NexusSortFilterProxyModel(self.view)
        # self.__treeModelSorted.setSourceModel(self.hdf5model)
        # self.__treeModelSorted.sort(0, qt.Qt.AscendingOrder)
        # self.__treeModelSorted.setSortCaseSensitivity(qt.Qt.CaseInsensitive)

        # self.hdfTreeView.setModel(self.__treeModelSorted)

        self.view.doubleClicked.connect(self._onNEXUSDoubleClicked)
        self.view.addContextMenuCallback(self.nexus_treeview_callback)

        try:
            self.openTemporaryDatabase()
        except Exception:
            # orGUI must start up even without a writable location for the
            # scratch database, e.g. if the drive with the working directory
            # is disconnected or read only.
            logger.error(
                "Cannot create the initial temporary database.",
                exc_info=True,
                extra={
                    "title": "Cannot create database",
                    "description": "Cannot create the initial temporary database. "
                    "Create a new database before saving any data.",
                    "show_dialog": True,
                    "parent": self,
                },
            )
        # self.add_nxdict(gauss)
        # self.add_nxdict(gauss)
        self.setCentralWidget(self.view)

        toolbar = qt.QToolBar("Database toolbar", self)

        newDatabaseAct = toolbar.addAction(
            resources.getQicon("document-nx-new"), "Create new orgui database"
        )
        newDatabaseAct.triggered.connect(self.onNewDatabase)

        loadDatabaseAct = toolbar.addAction(
            icons.getQIcon("document-open"), "Open orgui database"
        )
        loadDatabaseAct.triggered.connect(self.onOpenDatabase)

        savenewact = toolbar.addAction(
            icons.getQIcon("layer-nx"), "Select orgui database location"
        )
        savenewact.triggered.connect(self.onSaveNewDBFile)

        saveact = toolbar.addAction(
            icons.getQIcon("document-save"), "Save orgui database"
        )
        saveact.triggered.connect(self.onSaveDBFile)

        self.addToolBar(toolbar)

    def _onNEXUSDoubleClicked(self, index):
        try:
            nodes = list(self.view.selectedH5Nodes())
            if len(nodes) > 0:
                obj = nodes[0]
                if obj.ntype is h5py.Dataset:
                    roi_node = self.get_roinode(obj.h5py_object)
                    if roi_node is not None:
                        self.plot_signal_callback(roi_node, obj)
                        return
                self.plot_default(obj.h5py_object)
        except Exception as e:
            # reading the database can fail at any time, e.g. if the drive
            # with the database was disconnected.
            logger.error(
                "Cannot read the selected database entry.",
                exc_info=True,
                extra={
                    "title": "Cannot read database",
                    "description": f"Cannot read the selected database entry:\n{e}",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )

    def plot_default(self, h5py_object):
        try:
            nxdat = nxdata.get_default(h5py_object)
            # print(nxdat.axes, nxdat.signal, nxdat.title)
            # print(nxdat.signal_name, nxdat.errors, nxdat.axes_names)
            if nxdat is not None and len(nxdat.axes) == 1:
                self.plot.addCurve(
                    nxdat.axes[0],
                    nxdat.signal,
                    legend=nxdat.title,
                    xlabel=nxdat.axes_names[0],
                    ylabel=nxdat.signal_name,
                    yerror=nxdat.errors,
                )
        except Exception as e:
            traceback.print_exc()
            print(f"Cannot plot data: {e}")

    def get_roinode(self, obj):
        if silx.io.utils.get_h5_class(obj) is None:
            return None

        while obj.name != "/":
            meta = obj.attrs.get("orgui_meta", False)
            if meta and meta == "roi":
                return obj
            obj = obj.parent

    def get_scannode(self, obj):
        if silx.io.utils.get_h5_class(obj) is None:
            return None

        while obj.name != "/":
            meta = obj.attrs.get("orgui_meta", False)
            if meta and meta == "scan":
                return obj
            obj = obj.parent

    def nexus_treeview_callback(self, event):
        """Fill the context menu of the database tree view.

        Building the menu reads attributes from the database file, which can
        fail at any time, e.g. if the drive with the database was
        disconnected. Such errors must not escape into the Qt event loop.
        """
        try:
            self._fillNexusTreeviewMenu(event)
        except Exception as e:
            logger.error(
                "Cannot create the context menu of the database entry.",
                exc_info=True,
                extra={
                    "title": "Cannot read database",
                    "description": f"Cannot read the selected database entry:\n{e}",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )

    def _fillNexusTreeviewMenu(self, event):
        objects = list(event.source().selectedH5Nodes())
        if not objects:
            return
        obj = objects[0]  # for single selection
        menu = event.menu()
        action = qt.QAction("Refresh", menu)
        action.triggered.connect(lambda: self.onRefreshNode(obj))
        menu.addAction(action)
        action = qt.QAction("details", menu)
        action.triggered.connect(lambda: self.view_data_callback(obj))
        menu.addAction(action)

        if obj.ntype is h5py.Dataset:
            roi_node = self.get_roinode(obj.h5py_object)
            if roi_node is not None:
                action = qt.QAction(f"plot {obj.name}", menu)
                action.triggered.connect(
                    lambda: self.plot_signal_callback(roi_node, obj)
                )
                menu.addAction(action)

        if obj.ntype is h5py.Group:
            meta = obj.h5py_object.attrs.get("orgui_meta", False)
            if ConfigHandler.is_config_group(obj.h5py_object):
                action = qt.QAction("Load configuration", menu)
                action.triggered.connect(lambda: self.onLoadConfig(obj.h5py_object))
                menu.addAction(action)

            elif meta and "roi" in meta:
                if "rocking" in meta:
                    action = qt.QAction("Show in rocking integration", menu)
                    action.triggered.connect(lambda: self.onShowRoIntegrate(obj))
                    menu.addAction(action)
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda: self.onDeleteNode(obj.h5py_object))
                menu.addAction(action)

            elif meta and "rocking" in meta:
                action = qt.QAction("Show in rocking integration", menu)
                action.triggered.connect(lambda: self.onShowRoIntegrate(obj))
                menu.addAction(action)
                action = qt.QAction("rename", menu)
                action.triggered.connect(lambda: self.onRenameNode(obj.h5py_object))
                menu.addAction(action)
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda: self.onDeleteNode(obj.h5py_object))
                menu.addAction(action)

            if meta and "scan" in meta:
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda: self.onDeleteScan(obj.h5py_object))
                menu.addAction(action)

    def onLoadConfig(self, obj):
        """Apply a selected stored configuration to the configured GUI target."""
        try:
            self.config_handler.apply_config_group(obj, self.config_target)
        except Exception:
            logger.exception(
                "Cannot load configuration from database.",
                extra={
                    "title": "Cannot load configuration",
                    "description": f"Cannot load configuration group {obj.name}",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )

    def plot_signal_callback(self, roi_node, dataset):
        try:
            nxdat = nxdata.get_default(roi_node)
            data = dataset.h5py_object[()]
            # print(nxdat.axes, nxdat.signal, nxdat.title)
            # print(nxdat.signal_name, nxdat.errors, nxdat.axes_names)
            if nxdat is not None and len(nxdat.axes) == 1:
                self.plot.addCurve(
                    nxdat.axes[0],
                    data,
                    legend=nxdat.title + "_" + dataset.name,
                    xlabel=nxdat.axes_names[0],
                    ylabel=dataset.name,
                )
        except Exception as e:
            traceback.print_exc()
            print(f"Cannot plot data: {e}")

    def view_data_callback(self, obj):
        self.dataviewer.setData(obj)
        self.dataviewerDialog.open()

    def onSaveNewDBFile(self):
        fileTypeDict = {"NEXUS Files (*.h5)": ".h5", "All files (*)": ""}
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getSaveFileName(
            self, "Select database location", self.filedialogdir, fileTypeFilter[:-2]
        )
        if filename == "":
            return

        self.filedialogdir = os.path.splitext(filename)[0]
        # filename += fileTypeDict[filetype]

        try:
            self.saveNewDBFile(filename)
        except DBCloseError:
            logger.exception(
                "Cannot close old db file. The database might be corrupted!",
                extra={
                    "title": "The database might be corrupted",
                    "description": "Will proceed with new data base.",
                    "show_dialog": True,
                    "parent": self,
                },
            )
            try:
                self.saveNewDBFile(filename)
            except Exception:
                logger.exception(
                    "Cannot create db file.",
                    extra={
                        "title": "Cannot create new db file",
                        "show_dialog": True,
                        "parent": self,
                    },
                )
        except Exception:
            logger.exception(
                "Cannot create db file.",
                extra={
                    "title": "Cannot create new db file",
                    "show_dialog": True,
                    "parent": self,
                },
            )
        finally:
            self._ensureDatabaseAvailable()

    def onNewDatabase(self):
        fileTypeDict = {"NEXUS Files (*.h5)": ".h5", "All files (*)": ""}
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getSaveFileName(
            self, "Select database location", self.filedialogdir, fileTypeFilter[:-2]
        )
        if filename == "":
            return
        self.filedialogdir = os.path.splitext(filename)[0]
        try:
            self.createNewDBFile(filename)
        except DBCloseError:
            logger.exception(
                "Cannot close old db file. The database might be corrupted!",
                extra={
                    "title": "The database might be corrupted",
                    "description": "Will proceed with new data base.",
                    "show_dialog": True,
                    "parent": self,
                },
            )
            try:
                self.createNewDBFile(filename)
            except Exception:
                logger.exception(
                    "Cannot create db file.",
                    extra={
                        "title": "Cannot create new db file",
                        "description": f"Filename: {filename}",
                        "show_dialog": True,
                        "parent": self,
                    },
                )
        except Exception:
            logger.exception(
                "Cannot create db file.",
                extra={
                    "title": "Cannot create new db file",
                    "description": f"Filename: {filename}",
                    "show_dialog": True,
                    "parent": self,
                },
            )
        finally:
            self._ensureDatabaseAvailable()

    def onSaveDBFile(self):
        fileTypeDict = {"NEXUS Files (*.h5)": ".h5", "All files (*)": ""}
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getSaveFileName(
            self, "Save database", self.filedialogdir, fileTypeFilter[:-2]
        )
        if filename == "":
            return

        self.filedialogdir = os.path.splitext(filename)[0]
        # filename += fileTypeDict[filetype]
        try:
            self.saveDBFile(filename)

        except Exception:
            logger.exception(
                "Cannot save db file.",
                extra={
                    "title": "Cannot save data base file",
                    "show_dialog": True,
                    "description": f"Filename: {filename}",
                    "parent": self,
                },
            )

    def onOpenDatabase(self):
        fileTypeDict = {"NEXUS Files (*.h5)": ".h5", "All files (*)": ""}
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getOpenFileName(
            self, "Save database", self.filedialogdir, fileTypeFilter[:-2]
        )
        if filename == "":
            return

        self.filedialogdir = os.path.splitext(filename)[0]
        # filename += fileTypeDict[filetype]
        try:
            self.openDBFile(filename)
        except Exception:
            logger.exception(
                "Cannot open data base file.",
                extra={
                    "title": "Cannot open data base file",
                    "show_dialog": True,
                    "description": f"Filename: {filename}",
                    "parent": self,
                },
            )
        finally:
            self._ensureDatabaseAvailable()

    def saveNewDBFile(self, filename):
        alldata = nxtodict(self._requireOpenFile())
        self.createNewDBFile(filename, alldata)

    def saveDBFile(self, filename):
        alldata = nxtodict(self._requireOpenFile())
        dicttonx(
            alldata, filename, create_dataset_args={"compression": self.compression}
        )

    def createNewDBFile(self, filename, datadict=None):
        self._closePreviousDBFile()

        fileattrs = {
            "@NX_class": "NXroot",
            "@creator": f"orGUI version {__version__}",
            "@file_name": str(os.path.basename(filename)),
            "@file_time": datetime.datetime.utcnow().isoformat(),
        }
        if datadict is None:
            datadict = fileattrs
        else:
            datadict.update(fileattrs)

        dicttonx(
            datadict, filename, create_dataset_args={"compression": self.compression}
        )
        self.openDBFile(filename)

    def openDBFile(self, filename):
        self._closePreviousDBFile()
        self.nxfile = silx.io.h5py_utils.File(filename, "a")
        self._filepath = filename
        while self.hdf5model.hasPendingOperations():
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.insertH5pyObject(self.nxfile)
        self.view.expandToDepth(0)

    def _closePreviousDBFile(self):
        """Close the currently open database file, if there is any.

        :raises DBCloseError:
            If the database file cannot be closed properly. The database is
            detached from orGUI in any case, so that a new database file can
            be created or opened afterwards.
        """
        if self.nxfile is None:
            return
        try:
            self.close()
        except Exception as e:
            raise DBCloseError(
                "Cannot close previous database file.\nThe database might be corrupted."  # noqa: E501
            ) from e  # convert to common IOError since can also be RuntimeError

    def isOpen(self):
        """Whether a usable database file is currently open.

        Returns False if the database was never opened, or if the file handle
        was closed, e.g. after a failed write to a disconnected drive.

        :rtype: bool
        """
        return bool(self.nxfile)

    def _requireOpenFile(self):
        """Return the open database file.

        :raises DBUnavailableError: If no database file is open.
        """
        if not self.isOpen():
            raise DBUnavailableError(
                "No database file open. Create or open a database first."
            )
        return self.nxfile

    def openTemporaryDatabase(self):
        """Create and open a scratch database in a temporary directory.

        Used at startup and as a fallback if the database on disk became
        unavailable, so that orGUI stays operational.
        """
        # close first: this discards the previous temporary directory
        self._closePreviousDBFile()
        self.temp_directory = self._createTempDirectory()
        self.createNewDBFile(
            os.path.join(self.temp_directory.name, "orgui_database.h5")
        )

    @staticmethod
    def _createTempDirectory():
        """Create the temporary directory of the scratch database.

        The working directory is preferred, since it is usually located on the
        same drive as the measured data. If it is not writable, e.g. because
        the drive was disconnected, the system temporary directory is used.
        """
        try:
            return tempfile.TemporaryDirectory(
                dir=os.getcwd(), ignore_cleanup_errors=True
            )
        except OSError:
            logger.warning(
                "Cannot create a temporary database directory in %s. "
                "Using the system temporary directory instead.",
                os.getcwd(),
                exc_info=True,
            )
            return tempfile.TemporaryDirectory(ignore_cleanup_errors=True)

    def _ensureDatabaseAvailable(self):
        """Fall back to a temporary database if no database file is open.

        Called after a failed attempt to create or open a database file, so
        that data can still be saved somewhere and orGUI keeps operating.
        """
        if self.isOpen():
            return
        try:
            self.openTemporaryDatabase()
        except Exception:
            logger.error(
                "Cannot create a temporary database.",
                exc_info=True,
                extra={
                    "title": "No database available",
                    "description": "Cannot create a temporary database. "
                    "No data can be saved until a new database is created.",
                    "show_dialog": True,
                    "parent": self,
                },
            )
            return
        logger.warning(
            "No database file available, created the temporary database %s.",
            self._filepath,
            extra={
                "title": "Temporary database created",
                "description": "No database file is available. Created the "
                f"temporary database\n{self._filepath}\n"
                "Save the database to a permanent location to keep the data.",
                "show_dialog": True,
                "dialog_level": logging.WARNING,
                "parent": self,
            },
        )

    def add_nxdict(self, nxentry, update_mode="add", h5path="/"):
        nxfile = self._requireOpenFile()
        dicttonx(
            nxentry,
            nxfile,
            h5path=h5path,
            update_mode="add",
            create_dataset_args={"compression": self.compression},
        )
        while self.hdf5model.hasPendingOperations():
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.synchronizeH5pyObject(nxfile)
        self.view.expandToDepth(0)

    def register_external_result(
        self, path, checksum, grids, status, job_digest
    ):
        """Register a standalone reconstruction without copying its datasets."""
        self.add_nxdict(
            {
                "external_reconstructions": {
                    job_digest[:16]: {
                        "@NX_class": "NXnote",
                        "@orgui_meta": "external_reconstruction",
                        "path": os.path.abspath(path),
                        "sha256": checksum,
                        "status": status,
                        "grids_json": json.dumps(grids, sort_keys=True),
                        "job_sha256": job_digest,
                    }
                }
            }
        )

    def write_scan_config(self, scan_name, config=None):
        """Write the default config group for a scan."""
        if config is None:
            config = ConfigData.from_gui(self.config_target)
        group = self._requireOpenFile()[scan_name]
        self.config_handler.create_dataset_args = {"compression": self.compression}
        return self.config_handler.write_scan_config(group, config)

    def write_integration_config(self, integration_path, config=None):
        """Write the config group used for an integration result."""
        if config is None:
            config = ConfigData.from_gui(self.config_target)
        group = self._requireOpenFile()[integration_path]
        self.config_handler.create_dataset_args = {"compression": self.compression}
        return self.config_handler.write_integration_config(group, config)

    def onDeleteScan(self, obj):
        btn = qt.QMessageBox.question(
            self,
            "Delete scan?",
            "Are you sure that you want to delete {} from the orgui database?".format(
                obj.name.split("/")[-1]
            ),  # noqa: E501
        )
        if btn == qt.QMessageBox.Yes:
            self.onDeleteNode(obj)

    def onDeleteNode(self, obj):
        """Delete a node from the database and report errors to the user."""
        try:
            self.delete_node(obj)
        except Exception as e:
            logger.error(
                "Cannot delete entry from the database.",
                exc_info=True,
                extra={
                    "title": "Cannot delete database entry",
                    "description": f"Cannot delete entry from the database:\n{e}",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )

    def onRefreshNode(self, obj):
        """Reload a database entry in the tree view, reporting errors."""
        try:
            self.hdf5model.synchronizeH5pyObject(obj)
        except Exception as e:
            logger.error(
                "Cannot refresh the database entry.",
                exc_info=True,
                extra={
                    "title": "Cannot refresh database entry",
                    "description": f"Cannot refresh the database entry:\n{e}",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )

    def onShowRoIntegrate(self, obj):
        try:
            h5_obj = self._findRockingScanGroup(obj)
        except Exception as e:
            logger.error(
                "Cannot read the rocking scan from the database.",
                exc_info=True,
                extra={
                    "title": "Cannot read database",
                    "description": f"Cannot read the rocking scan:\n{e}",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return
        if h5_obj is None:
            return  # invalid dataset
        self.sigChangeRockingScan.emit(h5_obj.name)

    def _findRockingScanGroup(self, obj):
        """Return the rocking scan group of a database entry, or None."""
        meta = obj.h5py_object.attrs.get("orgui_meta", False)
        h5_obj = obj.h5py_object
        while h5_obj.name != "/":  # search for rocking scan group
            meta = h5_obj.attrs.get("orgui_meta", False)
            if meta and meta == "rocking":
                break
            h5_obj = h5_obj.parent
        else:
            dir_name = h5_obj.name
            msgbox = qt.QMessageBox(
                qt.QMessageBox.Critical,
                "Invalid rocking scan",
                f"Not a rocking scan: {dir_name}.",
                qt.QMessageBox.Ok,
                self,
            )
            msgbox.exec()
            return None  # invalid dataset

        return h5_obj

    def delete_node(self, obj):
        nxfile = self._requireOpenFile()
        basename = obj.name.split("/")[-1]
        objpar = obj.parent
        del objpar[basename]
        while self.hdf5model.hasPendingOperations():
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.synchronizeH5pyObject(nxfile)
        self.view.expandToDepth(0)

    def onRenameNode(self, obj):
        basename = obj.name.split("/")[-1]
        newname, success = qt.QInputDialog.getText(
            self,
            "Rename NEXUS node",
            "New name:",
            qt.QLineEdit.EchoMode.Normal,
            basename,
        )
        if success and newname != "":
            try:
                self.rename_node(obj, newname)
            except Exception as e:
                logger.error(
                    "Cannot rename entry in the database.",
                    exc_info=True,
                    extra={
                        "title": "Cannot rename database entry",
                        "description": f"Cannot rename entry {basename}:\n{e}",
                        "show_dialog": True,
                        "dialog_level": logging.WARNING,
                        "parent": self,
                    },
                )

    def rename_node(self, obj, newname):
        nxfile = self._requireOpenFile()
        basename = obj.name.split("/")[-1]
        objpar = obj.parent
        objpar.move(basename, newname)
        while self.hdf5model.hasPendingOperations():
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.synchronizeH5pyObject(nxfile)
        self.view.expandToDepth(0)

    def close(self):
        """Close the database file and detach it from the tree view.

        orGUI always detaches from the database file, even if it cannot be
        closed properly, e.g. if the drive with the database was disconnected.
        A new database can therefore be created or opened afterwards, also
        after this method raised an error.

        :raises DBCloseError:
            If pending database operations do not finish in time. The database
            file is not closed in this case, since closing it while the tree
            model reads it can crash the application.
        :raises Exception:
            The error raised by h5py if the file cannot be closed properly,
            usually a RuntimeError or an OSError. The database file might be
            corrupted in this case.
        """
        if self.nxfile is None:
            return
        # keep the file open on timeout, the tree model might still read it.
        self._waitForPendingOperations()

        nxfile = self.nxfile
        self.nxfile = None  # the file handle must not be used any longer
        try:
            self._detachFromTreeView(nxfile)
            self._waitForPendingOperations()
            try:
                nxfile.close()
            except Exception:
                logger.error(
                    "Closing of database file failed. The database file might be corrupted!",  # noqa: E501
                    exc_info=True,
                )
                raise
        finally:
            self._discardTempDirectory()

    def closeSafe(self, show_dialog=False):
        """Close the database file, without raising on error.

        Used on shutdown paths, where a failed close must not escape into the
        Qt event loop or the exception hook.

        :param bool show_dialog:
            Notify the user with a message box if the database file could not
            be closed properly.
        :returns: True if the database file was closed properly.
        :rtype: bool
        """
        try:
            self.close()
        except Exception:
            # logged as a warning on purpose: in cli context, error records
            # re-raise the exception, see orgui.logger_utils.
            logger.warning(
                "Cannot close the database file. The database might be corrupted!",
                exc_info=True,
                extra={
                    "title": "The database might be corrupted",
                    "description": "Cannot close the database file properly. "
                    "The data written since the last save might be lost.",
                    "show_dialog": show_dialog,
                    "dialog_level": logging.ERROR,
                    "parent": self,
                },
            )
            return False
        return True

    def _waitForPendingOperations(self, timeout=60.0):
        """Wait until the tree model finished its asynchronous operations.

        :param float timeout: Maximum waiting time in seconds.
        :raises DBCloseError: If the operations did not finish in time.
        """
        deadline = time.monotonic() + timeout
        while self.hdf5model.hasPendingOperations():
            if time.monotonic() > deadline:
                raise DBCloseError(
                    "Timeout on hdf5 model operation. This is probably a bug, or a "
                    "very long writing operation occurs. Please report this if no "
                    "long writing operation is running."
                )
            # the pending operations of the tree model only finish if Qt can
            # process events, see also add_nxdict.
            qt.QApplication.processEvents()
            time.sleep(0.01)

    def _detachFromTreeView(self, nxfile):
        """Remove a database file from the tree view. Does not raise.

        Removing the item reads from the file, which can fail if the file
        became unavailable. orGUI has to detach from it in any case.
        """
        # logged as warnings on purpose: in cli context, error records re-raise
        # the exception, see orgui.logger_utils.
        try:
            self.hdf5model.removeH5pyObject(nxfile)
        except Exception:
            logger.warning(
                "Cannot remove the database from the tree view. "
                "Clearing the tree view instead.",
                exc_info=True,
            )
            try:
                self.hdf5model.clear()
            except Exception:
                logger.warning("Cannot clear the database tree view.", exc_info=True)

    def _discardTempDirectory(self):
        """Delete the temporary directory of the scratch database, if any."""
        temp_directory = getattr(self, "temp_directory", None)
        if temp_directory is None:
            return
        del self.temp_directory
        try:
            temp_directory.cleanup()
        except Exception:
            logger.warning(
                "Cannot remove the temporary database directory %s.",
                temp_directory.name,
                exc_info=True,
            )

# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
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
__credits__ = ["Finn Schroeter"]
__copyright__ = "Copyright 2020-2026 Timo Fuchs"
__license__ = "MIT License"
from .. import __version__

__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import logging
from .. import logger_utils

import copy
import sys
import os
import re
from silx.gui import qt

from io import StringIO
import concurrent.futures
import threading

# from IPython import embed
import silx.gui.plot
from silx.gui.plot import items
from silx.gui.colors import Colormap
import fabio

# from silx import sx

import silx
from silx.utils.weakref import WeakMethodProxy
from silx.gui.plot.Profile import ProfileToolBar
from silx.gui.plot.tools.roi import RegionOfInterestManager
from silx.gui.plot.actions import control as control_actions

import traceback

from . import qutils, ROIutils, autoBraggSearch
from .QScanSelector import QScanSelector
from .QReflectionSelector import QReflectionSelector, QReflectionAnglesDialog
from .QUBCalculator import QUBCalculator
from .peak1Dintegr import RockingPeakIntegrator
from .ArrayTableDialog import ArrayTableDialog
from .bgroi import RectangleBgROI
from .database import DataBase, FILTERS
from ..backend.scans import SimulationScan
from ..backend import backends
from ..backend import universalScanLoader
from ..backend import interlacedScanLoader
from .. import resources

import numpy as np
from ..datautils.xrayutils import HKLVlieg, CTRcalc
from ..datautils.xrayutils import ReciprocalNavigation as rn

# legacy import:
from ..backend.beamline.id31_tools import Fastscan
from ..backend.beamline.P212_tools import FioFastsweep

logger = logging.getLogger(__name__)  # noqa

try:
    from silx.gui import console
except Exception:
    console = False

try:
    from . import _roi_sum_accel

    HAS_ACCEL = True
except Exception as exc:
    logger.debug("ROI sum accelerator unavailable: %s", exc)
    HAS_ACCEL = False

QTVERSION = qt.qVersion()
DEBUG = 0

MAX_ROIS_DISPLAY = 100

MAX_MEMORY = 6000  # MB,
try:
    import psutil

    memory_info = psutil.virtual_memory()
    avail_memory = memory_info.total / (1024**2)
    MAX_MEMORY = avail_memory * 0.7  # cap memory at 70%
except Exception:
    print("Cannot retrieve available memory size. Cap usage to 6 GB.")

silx.config.DEFAULT_PLOT_SYMBOL = "."


def _display_roi_geometry(center, left, right, top, bottom):
    """Convert detector-array ROI slices to plot rectangle geometry.

    Detector row indices increase toward the displayed bottom, while plot
    coordinates increase toward the displayed top. Therefore the array
    ``bottom`` extent is the plot ``top`` extent and vice versa.
    """
    origin = (center[0].start, center[1].start)
    size = (
        center[0].stop - center[0].start,
        center[1].stop - center[1].start,
    )
    left_size = left[0].stop - left[0].start
    right_size = right[0].stop - right[0].start
    top_size = bottom[1].stop - bottom[1].start
    bottom_size = top[1].stop - top[1].start
    return origin, size, left_size, right_size, top_size, bottom_size


class orGUI(qt.QMainWindow):
    def __init__(self, configfile, parent=None):
        """Initialize the main orGUI window.

        .. note::
           CLI-capable. CLI mode still creates this Qt main-window object, but
           it must not show blocking dialogs during construction.
        """
        qt.QMainWindow.__init__(self, parent)
        # self.setWindowIcon(resources.getQicon("orguiicon"))
        self.h5database = None  # must be a h5py file-like, by default not opened to avoid reading issues at beamtimes!  # noqa: E501
        self.images_loaded = False
        self.resetZoom = True
        self.background_image = None
        # icon = resources.getQicon("sum_image.svg")
        self.fscan = None
        self.activescanname = "scan"
        self.numberthreads = (
            int(min(os.cpu_count(), 16)) if os.cpu_count() is not None else 1
        )
        self.maxMemory = MAX_MEMORY
        self.maxROIs = MAX_ROIS_DISPLAY

        self.filedialogdir = os.getcwd()

        self.excludedImagesDialog = ArrayTableDialog(True, 1)
        self.excludedImagesDialog.setArrayData(
            np.array([-1]), editable=True, header=["image no"]
        )

        self.centralPlot = Plot2DHKL(self.newXyHKLConverter(), parent=self)
        self.centralPlot.setDefaultColormap(
            Colormap(name="inferno", normalization="log")
        )
        self.centralPlot.setCallback(self._graphCallback)
        toolbar = qt.QToolBar()
        toolbar.addAction(
            control_actions.OpenGLAction(parent=toolbar, plot=self.centralPlot)
        )
        self.centralPlot.addToolBar(toolbar)

        self.currentImageLabel = None
        self.currentAddImageLabel = None

        selectorDock = qt.QDockWidget("Scan data")
        selectorDock.setAllowedAreas(
            qt.Qt.LeftDockWidgetArea | qt.Qt.RightDockWidgetArea
        )
        self.scanSelector = QScanSelector(self)
        selectorDock.setWidget(self.scanSelector)
        self.addDockWidget(qt.Qt.LeftDockWidgetArea, selectorDock)

        self.imagepath = ""
        self.imageno = 0

        # ubWidget = qt.QSplitter(qt.Qt.Vertical)
        # ubWidget.setChildrenCollapsible(False)
        ubLayout = qt.QVBoxLayout()
        ubWidget = qt.QWidget()
        self.ubcalc = QUBCalculator(None, self)
        self.ubcalc.sigNewReflection.connect(self._onNewReflection)
        self.ubcalc.sigViewReflection.connect(self.onViewCalculatedReflection)
        self.ubcalc.viewReflectionS1Act.toggled.connect(
            lambda checked: checked
            and self.scanSelector.roiTrackingS1Act.setChecked(True)
        )
        self.ubcalc.viewReflectionS2Act.toggled.connect(
            lambda checked: checked
            and self.scanSelector.roiTrackingS2Act.setChecked(True)
        )
        self.scanSelector.roiTrackingS1Act.toggled.connect(
            lambda checked: checked and self.ubcalc.viewReflectionS1Act.setChecked(True)
        )
        self.scanSelector.roiTrackingS2Act.toggled.connect(
            lambda checked: checked and self.ubcalc.viewReflectionS2Act.setChecked(True)
        )

        self.resetIntegrPlotCurves = qt.QAction("Reset plot", self)
        self.resetIntegrPlotCurves.triggered.connect(self._removeAllIntegrPlotCurves)

        self.resetIntegrPlotCurveSet = qt.QAction("Remove curves", self)
        self.resetIntegrPlotCurveSet.triggered.connect(self._removeIntegrPlotCurveSet)

        maincentralwidget = qt.QTabWidget()

        self.integrdataPlot = silx.gui.plot.Plot1D(self)
        legendwidget = self.integrdataPlot.getLegendsDockWidget()

        plotRemovalBar = qt.QToolBar()
        plotRemovalBar.addAction(self.resetIntegrPlotCurves)
        plotRemovalBar.addAction(self.resetIntegrPlotCurveSet)
        self.integrdataPlot.addToolBar(plotRemovalBar)

        toolbar = qt.QToolBar()
        toolbar.addAction(
            control_actions.OpenGLAction(parent=toolbar, plot=self.integrdataPlot)
        )
        self.integrdataPlot.addToolBar(toolbar)

        self.integrdataPlot.addDockWidget(qt.Qt.RightDockWidgetArea, legendwidget)
        legendwidget.show()
        self.database = DataBase(self.integrdataPlot)
        dbdockwidget = qt.QDockWidget("Integrated data")
        dbdockwidget.setWidget(self.database)

        self.integrdataPlot.addDockWidget(qt.Qt.RightDockWidgetArea, dbdockwidget)

        self.roPkIntegrTab = RockingPeakIntegrator(self.database)

        self.scanSelector.sigImageNoChanged.connect(self._onSliderValueChanged)

        self.scanSelector.sigImagePathChanged.connect(self._onImagePathChanged)
        self.scanSelector.sigScanChanged.connect(self._onScanChanged)

        self.scanSelector.showMaxAct.toggled.connect(self._onMaxToggled)
        self.scanSelector.showSumAct.toggled.connect(self._onSumToggled)
        self.scanSelector.roiTrackingAct.toggled.connect(self._onROITrackingChanged)
        self.scanSelector.roiTrackingS1Act.triggered.connect(self._onROITrackingChanged)
        self.scanSelector.roiTrackingS2Act.triggered.connect(self._onROITrackingChanged)

        self.scanSelector.sigROIChanged.connect(self.updateROI)
        self.scanSelector.sigROIChanged.connect(self._onROITrackingChanged)
        self.scanSelector.sigROIintegrate.connect(self.integrateROI)
        self.scanSelector.sigSearchHKL.connect(self.onSearchHKLforStaticROI)

        self.scanSelector.excludeImageAct.toggled.connect(self._onToggleExcludeImage)

        toolbar = self.scanSelector.getScanToolbar()

        self.centralPlot.addToolBar(qt.Qt.BottomToolBarArea, toolbar)

        maincentralwidget.addTab(self.centralPlot, "Scan Image browser")
        maincentralwidget.addTab(self.integrdataPlot, "ROI integrated data")
        maincentralwidget.addTab(self.roPkIntegrTab, "Rocking scan integrate")

        self.setCentralWidget(maincentralwidget)

        # Create the object controlling the ROIs and set it up
        self.roiManager = RegionOfInterestManager(self.centralPlot)
        self.roiManager.setColor("pink")  # Set the color of ROI

        # self.roiTable = RegionOfInterestTableWidget()
        # self.roiTable.setRegionOfInterestManager(self.roiManager)

        self.rocking_rois = []

        self.roiS1 = RectangleBgROI()
        self.roiS1.setLineWidth(2)
        self.roiS1.setLineStyle("-")
        self.roiS1.setColor("red")
        self.roiS1.setBgStyle("pink", "-", 2.0)
        self.roiS1.setVisible(True)
        self.roiS1.setGeometry(origin=(0, 0), size=(0, 0))
        # S1 is also fixed roi, editing enabled when static ROI:
        self.roiS1.sigEditingFinished.connect(self._onStaticROIedited)

        self.roiS2 = RectangleBgROI()
        self.roiS2.setLineWidth(2)
        self.roiS2.setLineStyle("-")
        self.roiS2.setColor("red")
        self.roiS2.setBgStyle("pink", "-", 2.0)
        self.roiS2.setVisible(True)
        self.roiS2.setGeometry(origin=(0, 0), size=(0, 0))
        self.roiManager.addRoi(self.roiS1, useManagerColor=False)
        self.roiManager.addRoi(self.roiS2, useManagerColor=False)

        # self.reflTable.view._model.dataChanged.connect(printmodel)
        # self.reflTable.setArrayData(np.array([0,0,0,0,10,10],dtype=np.float64))
        ubDock = qt.QDockWidget("Reciprocal space navigation")
        ubDock.setAllowedAreas(
            qt.Qt.LeftDockWidgetArea
            | qt.Qt.RightDockWidgetArea
            | qt.Qt.BottomDockWidgetArea
        )

        self.reflectionSel = QReflectionSelector(self.centralPlot, self.ubcalc, self)
        self.reflectionSel.sigQueryImageChange.connect(self._onChangeImage)
        self.reflectionSel.sigQueryCenterPlot.connect(self._onCenterGraph)

        self.ubcalc.sigReflectionMismatchChanged.connect(
            self.reflectionSel.setReflectionMismatch
        )
        self.ubcalc.setReflectionHandler(self.getReflections)

        self.ubcalc.sigPlottableMachineParamsChanged.connect(self._onPlotMachineParams)
        self.ubcalc.sigReplotRequest.connect(self.updatePlotItems)
        self.allimgsum = None
        self.allimgmax = None
        self.reflectionSel.setSizePolicy(
            qt.QSizePolicy.Preferred, qt.QSizePolicy.Expanding
        )
        ubLayout.addWidget(self.reflectionSel)
        self.ubcalc.setSizePolicy(qt.QSizePolicy.Preferred, qt.QSizePolicy.Maximum)
        ubLayout.addWidget(self.ubcalc)

        ubWidget.setLayout(ubLayout)
        ubDock.setWidget(ubWidget)
        self.centralPlot.addDockWidget(qt.Qt.RightDockWidgetArea, ubDock)

        menu_bar = qt.QMenuBar()
        file = menu_bar.addMenu("&File")
        file.addAction(self.scanSelector.openFileAction)
        file.addAction(self.scanSelector.refreshFileAction)
        file.addAction(self.scanSelector.closeFileAction)
        file.addSeparator()

        self.folderToScan = file.addAction("Generate scan from images")
        self.folderToScan.triggered.connect(self._onLoadScanFromImages)

        self.interlacedScanMenu = file.addAction("Load interlaced/segmented scan")
        self.interlacedScanMenu.triggered.connect(self._onLoadInterlacedScan)
        file.addSeparator()

        self.loadImagesAct = file.addAction("reload images")
        self.loadImagesAct.triggered.connect(self._onLoadAll)

        config_menu = menu_bar.addMenu("&Config")
        loadConfigAct = qt.QAction(
            "Load config", self
        )  # connected with UBCalculator below
        # loadXtalAct = qt.QAction("Load Crystal file",self)
        machineParamsAct = qt.QAction("Machine parameters", self)
        machineParamsAct.setCheckable(True)
        xtalParamsAct = qt.QAction("Crystal parameters", self)
        xtalParamsAct.setCheckable(True)
        cpucountAct = qt.QAction("Set CPU count", self)
        roicountAct = qt.QAction("Set ROI count", self)

        loadConfigAct.triggered.connect(self.ubcalc._onLoadConfig)
        machineParamsAct.toggled.connect(
            lambda checked: self.ubcalc.machineDialog.setVisible(checked)
        )
        self.ubcalc.machineDialog.sigHide.connect(
            lambda: machineParamsAct.setChecked(False)
        )

        xtalParamsAct.toggled.connect(
            lambda checked: self.ubcalc.xtalDialog.setVisible(checked)
        )
        self.ubcalc.xtalDialog.sigHide.connect(lambda: xtalParamsAct.setChecked(False))

        cpucountAct.triggered.connect(self._onSelectCPUcount)
        roicountAct.triggered.connect(self._onSelectROIcount)

        self.autoLoadAct = qt.QAction("Auto load scans", self)
        self.autoLoadAct.setCheckable(True)
        self.autoLoadAct.setChecked(True)

        self.showExcludedImagesAct = qt.QAction("Excluded images", self)
        self.showExcludedImagesAct.setCheckable(True)
        self.showExcludedImagesAct.toggled.connect(
            lambda visible: self.excludedImagesDialog.setVisible(visible)
        )

        self.backgroundImageAct = qt.QAction("Subtract/select background image", self)
        self.backgroundImageAct.toggled.connect(self._onSetBackgroundImage)
        self.backgroundImageAct.setCheckable(True)
        self.backgroundImageAct.setChecked(False)

        self.dbCompressionAct = qt.QAction("Database compression", self)
        self.dbCompressionAct.triggered.connect(self._onChangeDBCompression)

        config_menu.addAction(loadConfigAct)
        # config_menu.addAction(loadXtalAct)
        config_menu.addSeparator()
        config_menu.addAction(machineParamsAct)
        config_menu.addAction(xtalParamsAct)
        config_menu.addSeparator()
        config_menu.addAction(cpucountAct)
        config_menu.addAction(self.dbCompressionAct)
        config_menu.addAction(self.autoLoadAct)
        config_menu.addSeparator()
        config_menu.addAction(self.showExcludedImagesAct)
        config_menu.addAction(self.backgroundImageAct)

        view_menu = menu_bar.addMenu("&View")
        showRefReflectionsAct = view_menu.addAction("reference reflections")
        showRefReflectionsAct.setCheckable(True)
        showRefReflectionsAct.setChecked(True)
        showRefReflectionsAct.toggled.connect(
            lambda checked: self.reflectionSel.setReferenceReflectionsVisible(checked)
        )

        showBraggAct = view_menu.addAction("allowed Bragg reflections")
        showBraggAct.setCheckable(True)
        showBraggAct.setChecked(False)
        showBraggAct.toggled.connect(self.onShowBragg)

        self.showROIAct = view_menu.addAction("show ROI")
        self.showROIAct.setCheckable(True)
        self.showROIAct.setChecked(False)
        self.showROIAct.toggled.connect(self.onShowROI)
        self.roivisible = False
        # self.scanSelector.showROICheckBox.addAction(showROIAct)

        self.showInterpolatedBgAct = view_menu.addAction("show interpolated bg")
        self.showInterpolatedBgAct.setCheckable(True)
        self.showInterpolatedBgAct.setChecked(False)
        self.showInterpolatedBgAct.setEnabled(False)
        self.showInterpolatedBgAct.toggled.connect(self.onShowInterpolatedBg)

        self.showCTRreflAct = view_menu.addAction("CTR reflections")
        self.showCTRreflAct.setCheckable(True)
        self.showCTRreflAct.setChecked(False)
        self.showCTRreflAct.toggled.connect(self.onShowCTRreflections)
        self.reflectionsVisible = False

        self.showMachineParamsAct = view_menu.addAction("machine parameters")
        self.showMachineParamsAct.setCheckable(True)
        self.showMachineParamsAct.setChecked(False)
        self.showMachineParamsAct.toggled.connect(self._onPlotMachineParams)

        view_menu.addSeparator()

        view_menu.addAction(roicountAct)

        view_menu.addSeparator()

        saveBraggAct = view_menu.addAction("save allowed Bragg reflections")
        saveBraggAct.setCheckable(False)
        saveBraggAct.triggered.connect(self.saveBraggRefl)
        if console:
            view_menu.addSeparator()

            custom_banner = f"""orGUI v. {__version__} console
Available variables:
orgui : top level gui
ub : gui for UB matrix and angle calculations
"""

            self.console_dockwidget = console.IPythonDockWidget(
                self,
                {"orgui": self, "ub": self.ubcalc, "ROIutils": ROIutils},
                custom_banner,
                "orGUI console",
            )

            self.console_dockwidget.setAllowedAreas(
                qt.Qt.LeftDockWidgetArea
                | qt.Qt.RightDockWidgetArea
                | qt.Qt.BottomDockWidgetArea
            )
            self.tabifyDockWidget(selectorDock, self.console_dockwidget)
            # self.addDockWidget(qt.Qt.LeftDockWidgetArea,self.console_dockwidget)
            self.console_dockwidget.setVisible(False)
            consoleViewAct = self.console_dockwidget.toggleViewAction()
            view_menu.addAction(consoleViewAct)

        ##############################

        editUAct = qt.QAction("Edit orientation matrix", self)
        editUAct.setCheckable(True)
        editUAct.toggled.connect(
            lambda checked: self.ubcalc.ueditDialog.setVisible(checked)
        )
        self.ubcalc.ueditDialog.sigHide.connect(lambda: editUAct.setChecked(False))

        calcCTRsAvailableAct = qt.QAction("Calculate available CTRs", self)
        calcCTRsAvailableAct.triggered.connect(self._onCalcAvailableCTR)
        rs = menu_bar.addMenu("&Reciprocal space")
        rs.addAction(calcCTRsAvailableAct)
        rs.addAction(editUAct)

        simul = menu_bar.addMenu("&Simulation")

        createScanAct = simul.addAction("Create dummy scan")
        createScanAct.triggered.connect(self._onCreateScan)

        helpmenu = menu_bar.addMenu("&Help")

        diffractAct = helpmenu.addAction("Diffraction geometry")
        diffractAct.triggered.connect(self._onShowDiffractionGeometry)

        helpmenu.addSeparator()

        aboutAct = helpmenu.addAction("About")
        aboutAct.triggered.connect(self._onShowAbout)

        aboutQtAct = helpmenu.addAction("About Qt")
        aboutQtAct.triggered.connect(lambda: qt.QMessageBox.aboutQt(self))

        self.setMenuBar(menu_bar)

        if configfile is not None:
            self.ubcalc.readConfig(configfile)

    def _removeAllIntegrPlotCurves(self):
        """CLI-capable: clear all integration curves from the plot widget."""
        # remove plotted curves
        # curveList = self.integrdataPlot.getAllCurves(just_legend=True)
        # for i in curveList:
        #    self.integrdataPlot.removeCurve(i)

        # remove plotted curves and curves that have been hidden
        itemList = self.integrdataPlot.getAllCurves(withhidden=True)
        for it in itemList:
            self.integrdataPlot.removeCurve(it)

    def _removeIntegrPlotCurveSet(self):
        """GUI-only: open a dialog to choose integration curves to hide/delete."""
        curveList = self.integrdataPlot.getAllCurves(just_legend=True, withhidden=True)

        # find out if curves are already hidden to later check boxes in QT GUI window
        hidden = np.zeros_like(curveList, dtype="?")
        for nr, i in enumerate(curveList):
            if self.integrdataPlot.getCurve(i).getVisibleBounds() is None:
                hidden[nr] = True

        d = QPlotDeleteWindow(curveList, hidden)
        if d.exec() == qt.QDialog.Accepted:
            if d.action == "delete":
                try:
                    for i, j in enumerate(d.boxes):
                        if j.isChecked():
                            self.integrdataPlot.removeCurve(curveList[i])
                            self.integrdataPlot.getLegendsDockWidget().updateLegends()
                except MemoryError:
                    logger.exception(
                        "Can not delete selected plots.",
                        extra={
                            "title": "Can not delete selected plots.",
                            "description": "Can not delete selected plots.",
                            "show_dialog": True,
                            "dialog_level": logging.WARNING,
                            "parent": self,
                        },
                    )
            elif d.action == "hide":
                try:
                    for i, j in enumerate(d.boxes):
                        if j.isChecked():
                            self.integrdataPlot.getCurve(curveList[i]).setVisible(False)
                            self.integrdataPlot.getLegendsDockWidget().updateLegends()
                except MemoryError:
                    logger.exception(
                        "Can not delete selected plots.",
                        extra={
                            "title": "Can not delete selected plots.",
                            "description": "Can not delete selected plots.",
                            "show_dialog": True,
                            "dialog_level": logging.WARNING,
                            "parent": self,
                        },
                    )

    def get_rocking_coordinates(
        self, H_0=None, H_1=None, maxValue=None, step_width=None, **kwargs
    ):
        """Calculate detector coordinates along a reciprocal-space line.

        The trajectory is the line :math:`\\vec{H}_0 + s\\vec{H}_1` in
        reciprocal space, where :math:`\\vec{H}_0` is a reciprocal-space origin
        vector, :math:`\\vec{H}_1` is a direction vector, and :math:`s` is a
        scalar trajectory coordinate.

        :param numpy.ndarray H_0:
            Starting reciprocal-space vector in r.l.u. Shape ``(3,)``.
        :param numpy.ndarray H_1:
            Reciprocal-space direction vector in r.l.u. Shape ``(3,)``.
        :param float maxValue:
            Maximum scalar :math:`s` value sampled along
            :math:`\\vec{H}_0 + s\\vec{H}_1`.
        :param float step_width:
            Step width in scalar ``s`` units.
        :returns:
            Reflection dictionary with detector pixel coordinates, masks, and
            trajectory values.
        :rtype: dict

        .. note::
           CLI-capable. Missing arguments are read from GUI widgets.
        """
        # going back to the more universal integration along H_0 + s*H_1 coordinates
        if H_0 is None:
            H_0 = self.scanSelector.ro_H_0_dialog.get_hkl()
        if H_1 is None:
            H_1 = (
                self.scanSelector.ro_H_1_dialog.get_hkl()
            )  # default to CTR integration
        if step_width is None:
            step_width = self.scanSelector.roscanDeltaS.value()
        if maxValue is None:
            maxValue = self.scanSelector.roscanMaxS.value()

        dc = self.ubcalc.detectorCal
        xoffset, yoffset = self.scanSelector.roioptions.get_offsets()
        xoffset = kwargs.get("xoffset", xoffset)
        yoffset = kwargs.get("yoffset", yoffset)

        step_nr = round(maxValue / step_width) + 1
        s_points = np.linspace(0, maxValue, step_nr)

        hkl_desired = (
            np.outer(H_1, s_points).T + H_0
        )  # F contiguous is faster in anglesToHKL

        refldict = self.ubcalc.calcReflection(hkl_desired)  # F contiguous is faster

        ymask1 = np.logical_and(
            refldict["xy_1"][..., 1] >= 0,
            refldict["xy_1"][..., 1] < dc.detector.shape[0],
        )
        xmask1 = np.logical_and(
            refldict["xy_1"][..., 0] >= 0,
            refldict["xy_1"][..., 0] < dc.detector.shape[1],
        )
        yxmask1 = np.logical_and(xmask1, ymask1)

        ymask2 = np.logical_and(
            refldict["xy_2"][..., 1] >= 0,
            refldict["xy_2"][..., 1] < dc.detector.shape[0],
        )
        xmask2 = np.logical_and(
            refldict["xy_2"][..., 0] >= 0,
            refldict["xy_2"][..., 0] < dc.detector.shape[1],
        )
        yxmask2 = np.logical_and(xmask2, ymask2)

        refldict["mask_1"] = yxmask1
        refldict["mask_2"] = yxmask2
        refldict["s"] = s_points

        if xoffset != 0.0 or yoffset != 0.0:
            # warnings.warn("Nonzero pixel offset selected. Experimental feature! Angles and hkl are incorrect!!!")  # noqa: E501
            refldict["xy_1"][..., 0] += xoffset
            refldict["xy_2"][..., 0] += xoffset
            refldict["xy_1"][..., 1] += yoffset
            refldict["xy_2"][..., 1] += yoffset
        refldict["H_0"] = H_0
        refldict["H_1"] = H_1
        return refldict

    def get_Bragg_rocking_coordinates(self, strainVec=None, **kwargs):
        """Calculate detector coordinates for Bragg rocking reflections.

        :param numpy.ndarray strainVec:
            Optional fractional strain applied component-wise to the direct
            lattice constants ``a`` before calculating Bragg coordinates. This
            scales each lattice constant as ``a_strained = a * (1 + strainVec)``
            and does not represent a full strain tensor. If omitted, GUI spin
            boxes are read and interpreted as percent strain.
        :returns:
            Reflection dictionary with hkl values in r.l.u., detector pixel
            coordinates, angles in rad, and masks.
        :rtype: dict

        .. note::
           CLI-capable for ``th`` scans when state is already configured.
        """
        if strainVec is None:
            strainVec = (
                np.array([h.value() for h in self.scanSelector.strain_Bragg]) / 100.0
            )

        xoffset, yoffset = self.scanSelector.roioptions.get_offsets()
        xoffset = kwargs.get("xoffset", xoffset)
        yoffset = kwargs.get("yoffset", yoffset)

        if self.fscan is not None:
            if self.fscan.axisname != "th":
                raise NotImplementedError(
                    f"Calculation of available Bragg reflections is not implemented for {self.fscan.axisname} - scans"  # noqa: E501
                )

            xtal = self.ubcalc.crystal
            ommin = np.deg2rad(np.amin(self.fscan.omega))
            ommax = np.deg2rad(np.amax(self.fscan.omega))
            dc = self.ubcalc.detectorCal
            mu = self.ubcalc.mu
            ub = self.ubcalc.ubCal
            chi = self.ubcalc.chi
            phi = self.ubcalc.phi
            xtal.setEnergy(ub.getEnergy() * 1e3)

            # apply strain:
            ub_strained = copy.deepcopy(ub)

            xtal_cp = copy.deepcopy(xtal)
            xtal_cp.a = xtal_cp.a * (1.0 + strainVec)
            ub_strained.setLattice(xtal_cp)

            if self.scanSelector.bragg_multiple_enable.isChecked():
                hkl_factor = np.array(
                    [h.value() for h in self.scanSelector.bragg_multiple]
                )
                a = xtal_cp.a
                a /= hkl_factor
                xtal_singleatom = CTRcalc.UnitCell(a, np.rad2deg(xtal_cp.alpha))
                xtal_singleatom.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
                xtal_singleatom.setEnergy(ub.getEnergy() * 1e3)
                ub_strained.setLattice(xtal_singleatom)

                try:
                    hkls, yx, angles = rn.thscanBragg(
                        xtal_singleatom,
                        ub_strained,
                        mu,
                        dc,
                        (ommin, ommax),
                        chi=chi,
                        phi=phi,
                    )
                    hkls = hkls.astype(np.float64)
                    hkls *= hkl_factor.T
                except Exception as e:
                    raise Exception("Cannot calculate Bragg reflections") from e

            else:
                try:
                    hkls, yx, angles = rn.thscanBragg(
                        xtal_cp, ub_strained, mu, dc, (ommin, ommax), chi=chi, phi=phi
                    )
                    hkls = hkls.astype(np.float64)
                except Exception as e:
                    raise Exception("Cannot calculate Bragg reflections") from e

        if hkls.size > 0:
            s_points = np.arange(hkls.shape[0])
            mask = np.ones(hkls.shape[0], dtype=bool)
        else:
            s_points = np.array([])
            mask = np.array([])

        refldict = {
            "hkl": hkls,
            "xy_1": yx[:, ::-1],
            "angles_1": angles,
            "s": s_points,
            "mask_1": mask,
        }
        if xoffset != 0.0 or yoffset != 0.0:
            # warnings.warn("Nonzero pixel offset selected. Experimental feature! Angles and hkl are incorrect!!!")  # noqa: E501
            refldict["xy_1"][..., 0] += xoffset
            refldict["xy_1"][..., 1] += yoffset

        return refldict

    def intkeys_rocking(self, refldict, **kwargs):
        """Build center ROI slice coordinates for rocking-scan integration.

        This converts detector coordinates from ``refldict`` into rectangular
        center-ROI slice tuples. The slice tuple order is detector
        ``(x, y)``/``(horizontal, vertical)``. NumPy images are indexed as
        ``image[y_slice, x_slice]``, so integration code applies these slices
        as ``key[::-1]``. Coordinates and ROI extents are clipped to the
        detector bounds.

        :param dict refldict:
            Reflection dictionary produced by a rocking coordinate helper. For
            CTR-style rocking scans, it represents points sampled along the
            reciprocal-space line :math:`\\vec{H}_0 + s\\vec{H}_1`. The selected
            intersection must provide ``xy_<intersect>`` detector coordinates
            in pixels, and normally ``mask_<intersect>`` detector-valid flags.
        :param int vsize:
            Optional vertical base ROI size in pixels. Defaults to the GUI
            value and, by default, may be adjusted by advanced ROI options such
            as detector inclination and projected sample size.
        :param int hsize:
            Optional horizontal base ROI size in pixels. Defaults to the GUI
            value and, by default, may be adjusted by advanced ROI options such
            as detector inclination and projected sample size.
        :param numpy.ndarray size_exact:
            Optional per-point ROI sizes in pixels with columns
            ``(hsize, vsize)``. If provided, advanced ROI size calculation is
            skipped and these exact per-point sizes are used.
        :param bool mask:
            Whether to apply ``mask_<intersect>`` before creating slices.
            Defaults to ``True``.
        :param bool autovsize:
            Whether to derive vertical ROI size from the median vertical
            spacing between neighboring detector coordinates. Defaults to the
            GUI setting.
        :param bool autohsize:
            Whether to derive horizontal ROI size from the median horizontal
            spacing between neighboring detector coordinates. Defaults to the
            GUI setting.
        :param int intersect:
            Detector-intersection index to use, usually ``1`` or ``2``.
            Defaults to the GUI-selected intersection.
        :returns:
            Dictionary containing ``center`` as a list of ``(x_slice, y_slice)``
            ROI coordinates, effective ``vsize`` and ``hsize`` in pixels, and
            ``size_exact`` when corrected per-point sizes were used.
        :rtype: dict

        .. note::
           CLI-capable. Defaults are read from GUI ROI controls.
        """
        vsize = kwargs.get("vsize", int(self.scanSelector.vsize.value()))
        hsize = kwargs.get("hsize", int(self.scanSelector.hsize.value()))
        size_exact = kwargs.get("size_exact", None)
        apply_mask = kwargs.get("mask", True)
        autoROIVsize = kwargs.get(
            "autovsize", self.scanSelector.autoROIVsize.isChecked()
        )
        autoROIHsize = kwargs.get(
            "autohsize", self.scanSelector.autoROIHsize.isChecked()
        )
        if self.scanSelector.intersS1Act.isChecked():
            intersect = 1
        elif self.scanSelector.intersS2Act.isChecked():
            intersect = 2
        else:
            intersect = 1  # default
        intersect = kwargs.get("intersect", intersect)

        xy = refldict[f"xy_{int(intersect)}"]

        if size_exact is None:
            roioptions = self.scanSelector.roioptions.get_parameters()
            if roioptions["DetectorInclination"] or roioptions["ProjectSampleSize"]:
                if roioptions["ProjectSampleSize"]:
                    size_exact = ROIutils.calc_corrections(
                        xy,
                        self.ubcalc.detectorCal,
                        np.array([hsize, vsize]),
                        roioptions,
                        roioptions["DetectorInclination"],
                        roioptions["factor"],
                    )
                else:
                    size_exact = ROIutils.calc_corrections(
                        xy,
                        self.ubcalc.detectorCal,
                        np.array([hsize, vsize]),
                        None,
                        roioptions["DetectorInclination"],
                        roioptions["factor"],
                    )

        if apply_mask:
            xy = xy[refldict[f"mask_{int(intersect)}"]]
            if size_exact is not None:
                size_exact = size_exact[refldict[f"mask_{int(intersect)}"]]

        step_nr = xy.shape[0]
        if step_nr == 0:
            return {"center": [], "vsize": vsize, "hsize": hsize}
        if step_nr > 1:
            if autoROIVsize:
                # dist_in_pixels = np.abs(xy[0][1] - xy[-1][1])
                dist_in_pixels = np.median(np.abs(np.diff(xy[:, 1])))
                # roi_vlength = np.ceil(dist_in_pixels/step_nr)
                vsize = int(np.ceil(dist_in_pixels))
                if size_exact is not None:
                    size_exact[:, 1] = vsize

            if autoROIHsize:
                # dist_in_pixels = np.abs(xy[0][1] - xy[-1][1])
                dist_in_pixels = np.median(np.abs(np.diff(xy[:, 0])))
                # roi_vlength = np.ceil(dist_in_pixels/step_nr)
                hsize = int(np.ceil(dist_in_pixels))
                if size_exact is not None:
                    size_exact[:, 0] = hsize

        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape

        coord_restr = np.clip(xy, [0, 0], [dethsize, detvsize])
        if size_exact is not None:
            vhalfsize = size_exact[:, 1] // 2
            hhalfsize = size_exact[:, 0] // 2
            fromcoords = np.round(coord_restr - np.vstack([hhalfsize, vhalfsize]).T)
            tocoords = np.round(coord_restr + np.vstack([hhalfsize, vhalfsize]).T)

            mask_hsize = (size_exact[:, 0] % 2).astype(bool)
            remainder_mask = coord_restr[:, 0] % 1 < 0.5
            tocoords[mask_hsize & remainder_mask, 0] += 1
            fromcoords[mask_hsize & ~remainder_mask, 0] -= 1

            mask_vsize = (size_exact[:, 1] % 2).astype(bool)
            remainder_mask = coord_restr[:, 1] % 1 < 0.5
            tocoords[mask_vsize & remainder_mask, 1] += 1
            fromcoords[mask_vsize & ~remainder_mask, 1] -= 1

        else:
            vhalfsize = vsize // 2
            hhalfsize = hsize // 2
            fromcoords = np.round(coord_restr - np.array([hhalfsize, vhalfsize]))
            tocoords = np.round(coord_restr + np.array([hhalfsize, vhalfsize]))
            if hsize % 2:
                remainder_mask = coord_restr[:, 0] % 1 < 0.5
                tocoords[remainder_mask, 0] += 1
                fromcoords[~remainder_mask, 0] -= 1
            if vsize % 2:
                remainder_mask = coord_restr[:, 1] % 1 < 0.5
                tocoords[remainder_mask, 1] += 1
                fromcoords[~remainder_mask, 1] -= 1

        fromcoords = np.clip(fromcoords, [0, 0], [dethsize, detvsize])
        tocoords = np.clip(tocoords, [0, 0], [dethsize, detvsize])

        locations = []
        for roifrom, toroi in zip(
            fromcoords, tocoords
        ):  # any way to do this with ndarray operations?
            locations.append(
                tuple(
                    slice(int(fromcoord), int(tocoord))
                    for fromcoord, tocoord in zip(roifrom, toroi)
                )
            )
        roi_dict = {"center": locations, "vsize": vsize, "hsize": hsize}

        if size_exact is not None:
            roi_dict["size_exact"] = size_exact
        return roi_dict

    def intbkgkeys_rocking(self, refldict, **kwargs):
        """Build center and background ROI slice coordinates.

        This extends :meth:`intkeys_rocking` by adding left, right, top, and
        bottom background ROIs around each center ROI. Slice tuple order stays
        in detector ``(x, y)``/``(horizontal, vertical)`` coordinates. NumPy
        images are indexed as ``image[y_slice, x_slice]``, so integration code
        applies these slices as ``key[::-1]``.

        :param dict refldict:
            Reflection dictionary produced by a rocking coordinate helper. For
            CTR-style rocking scans, it represents points sampled along the
            reciprocal-space line :math:`\\vec{H}_0 + s\\vec{H}_1`. The selected
            intersection must provide ``xy_<intersect>`` detector coordinates
            in pixels, and normally ``mask_<intersect>`` detector-valid flags.
        :param int left:
            Optional left-background width in detector pixels. Defaults to the
            GUI value.
        :param int right:
            Optional right-background width in detector pixels. Defaults to the
            GUI value.
        :param int top:
            Optional top-background height in detector pixels. Defaults to the
            GUI value.
        :param int bottom:
            Optional bottom-background height in detector pixels. Defaults to
            the GUI value.
        :param kwargs:
            Center-ROI options forwarded to :meth:`intkeys_rocking`, including
            ``vsize``, ``hsize``, ``size_exact``, ``mask``, ``autovsize``,
            ``autohsize``, and ``intersect``.
        :returns:
            Dictionary containing ``center``, ``left``, ``right``, ``top``, and
            ``bottom`` lists of ``(x_slice, y_slice)`` ROI coordinates,
            effective center-ROI ``vsize`` and ``hsize`` in pixels, and
            ``size_exact`` when corrected per-point sizes were used.
        :rtype: dict

        .. note::
           CLI-capable. Defaults are read from GUI ROI controls.
        """
        left = kwargs.get("left", int(self.scanSelector.left.value()))
        right = kwargs.get("right", int(self.scanSelector.right.value()))
        top = kwargs.get("top", int(self.scanSelector.top.value()))
        bottom = kwargs.get("bottom", int(self.scanSelector.bottom.value()))

        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape

        roi_dict = self.intkeys_rocking(refldict, **kwargs)
        crois = roi_dict["center"]
        leftrois = []
        rightrois = []
        toprois = []
        bottomrois = []

        for croi in crois:
            leftrois.append(
                (
                    slice(
                        int(np.clip(croi[0].start - left, 0, dethsize)), croi[0].start
                    ),
                    croi[1],
                )
            )
            rightrois.append(
                (
                    slice(
                        croi[0].stop, int(np.clip(croi[0].stop + right, 0, dethsize))
                    ),
                    croi[1],
                )
            )
            toprois.append(
                (
                    croi[0],
                    slice(
                        int(np.clip(croi[1].start - top, 0, detvsize)), croi[1].start
                    ),
                )
            )
            bottomrois.append(
                (
                    croi[0],
                    slice(
                        croi[1].stop, int(np.clip(croi[1].stop + bottom, 0, detvsize))
                    ),
                )
            )
        roi_dict["left"] = leftrois
        roi_dict["right"] = rightrois
        roi_dict["top"] = toprois
        roi_dict["bottom"] = bottomrois
        return roi_dict

    def rocking_extraction(self):
        """Start CTR-style rocking extraction along a reciprocal-space line.

        :math:`\\vec{H}_0` and :math:`\\vec{H}_1` are reciprocal-space vectors
        in r.l.u. read from the rocking-scan controls. They define the
        reciprocal-space line :math:`\\vec{H}_0 + s\\vec{H}_1` whose detector
        intersections are integrated. This function uses configuration from
        the UI elements, including the line vectors, intersection selection,
        ROI sizes, background sizes, masks, and correction options.

        :returns:
            Status dictionary describing success, cancellation, or error.
        :rtype: dict

        .. note::
           CLI-capable when scan, database, and ROI state are preconfigured.
        """
        logger.info("Start hklscan rocking integration")
        if self.fscan is None:  # or isinstance(self.fscan, SimulationScan):
            # In GUI mode this shows a dialog and then returns the status dict
            # below. In CLI mode logger.error intentionally raises through the
            # CLI logging handler so scripts fail loudly on missing scan state.
            logger.error(
                "No scan loaded.",
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Cannot integrate scan: No scan loaded.",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {"status": "error", "message": "no scan loaded"}

        try:
            refldict = self.get_rocking_coordinates()
        except Exception:
            logger.exception(
                "Rocking scan extraction not implemented for scan axis "
                + str(self.fscan.axisname),
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Rocking scan extraction not implemented for scan axis "  # noqa: E501
                    + str(self.fscan.axisname),
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {
                "status": "error",
                "message": "Rocking scan extraction not implemented for scan axis "
                + str(self.fscan.axisname),
                "traceback": traceback.format_exc(),
            }

        if self.scanSelector.intersS1Act.isChecked():
            intersect = 1
        elif self.scanSelector.intersS2Act.isChecked():
            intersect = 2
        else:
            intersect = 1  # default
        mask = refldict[f"mask_{intersect}"]
        xy = refldict[f"xy_{intersect}"][mask]

        refldict["angles"] = refldict[f"angles_{intersect}"][mask]
        refldict["s_masked"] = refldict["s"][mask]
        refldict["hkl_masked"] = refldict["hkl"][mask]

        roi_keys = self.intbkgkeys_rocking(refldict)
        hkl_del_gam = self.getStaticROIparams(xy)

        ro_name = "rocking_[{:.2f} {:.2f} {:.2f}]_H0_[{:.2f} {:.2f} {:.2f}]_H1".format(
            *refldict["H_0"], *refldict["H_1"]
        )

        logger.info(f"Start rocking integration of scan {ro_name}")
        return self.rocking_integrate(xy, roi_keys, hkl_del_gam, refldict, ro_name)

    def rocking_Bragg_extraction(self):
        """Integrate rocking-scan ROIs centered on Bragg peak coordinates.

        Bragg peak detector coordinates are calculated from the current
        crystal, UB, detector, strain, and ``th``-scan state. The valid
        Bragg coordinates are converted to center and
        background ROIs, then integrated across the rocking scan.

        :returns:
            Status dictionary describing success, cancellation, or error.
        :rtype: dict

        .. note::
           CLI-capable when scan, database, and ROI state are preconfigured.
        """
        logger.info("Start Bragg rocking integration")
        if self.fscan is None:  # or isinstance(self.fscan, SimulationScan):
            logger.error(
                "No scan loaded.",
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Cannot integrate scan: No scan loaded.",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {"status": "error", "message": "no scan loaded"}

        try:
            refldict = self.get_Bragg_rocking_coordinates()
        except Exception:
            logger.exception(
                "Rocking scan extraction not implemented for scan axis "
                + str(self.fscan.axisname),
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Rocking scan extraction not implemented for scan axis "  # noqa: E501
                    + str(self.fscan.axisname),
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {
                "status": "error",
                "message": "Rocking scan extraction not implemented for scan axis "
                + str(self.fscan.axisname),
                "traceback": traceback.format_exc(),
            }

        mask = refldict["mask_1"]
        xy = refldict["xy_1"][mask]

        refldict["angles"] = refldict["angles_1"][mask]
        refldict["s_masked"] = refldict["s"][mask]
        refldict["hkl_masked"] = refldict["hkl"][mask]

        roi_keys = self.intbkgkeys_rocking(
            refldict, autovsize=False, autohsize=False, intersect=1
        )
        hkl_del_gam = self.getStaticROIparams(xy)

        ro_name = "rocking_Bragg"
        return self.rocking_integrate(xy, roi_keys, hkl_del_gam, refldict, ro_name)

    def rocking_static_extraction(self, xy, hsize, vsize):
        """Integrate fixed detector-pixel ROIs through a rocking scan.

        This path is for ROIs defined directly in detector pixel coordinates
        rather than from reciprocal-space line or Bragg-coordinate searches.
        Each input coordinate is treated as a fixed ROI center for every image
        in the scan. The provided ROI sizes are passed as exact per-ROI sizes,
        so automatic ROI sizing and advanced ROI size corrections are disabled.

        :param numpy.ndarray xy:
            ROI center coordinates in detector pixels. Shape ``(N, 2)`` for
            multiple ROIs, with columns ``(x, y)``.
        :param numpy.ndarray hsize:
            Horizontal ROI sizes in pixels, one value per ROI.
        :param numpy.ndarray vsize:
            Vertical ROI sizes in pixels, one value per ROI.
        :returns:
            Status dictionary from the rocking integration.
        :rtype: dict

        .. note::
           CLI-capable when scan and database state are preconfigured.
        """
        logger.info("Start static fixed integration")
        if self.fscan is None:  # or isinstance(self.fscan, SimulationScan):
            logger.error(
                "No scan loaded.",
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Cannot integrate scan: No scan loaded.",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {"status": "error", "message": "no scan loaded"}

        hkl_del_gam = self.getStaticROIparams(xy)
        refldict = {"xy_1": xy}
        size_exact = np.vstack([hsize, vsize]).T

        roi_keys = self.intbkgkeys_rocking(
            refldict,
            autovsize=False,
            autohsize=False,
            intersect=1,
            mask=False,
            size_exact=size_exact,
        )

        ro_name = "rocking_static"
        return self.rocking_integrate(xy, roi_keys, hkl_del_gam, refldict, ro_name)

    def rocking_integrate(self, xylist, rois, hkl_del_gam, refldict, name):
        """Integrate rocking-scan images over prepared ROIs.

        For CTR-style rocking scans, ``refldict`` may contain
        :math:`\\vec{H}_0` and :math:`\\vec{H}_1` reciprocal-space vectors in
        r.l.u.; together they define the line
        :math:`\\vec{H}_0 + s\\vec{H}_1` that is recorded in the output
        trajectory.

        :param numpy.ndarray xylist:
            ROI center coordinates in detector pixels.
        :param dict rois:
            Center and background ROI slice definitions.
        :param numpy.ndarray hkl_del_gam:
            Per-point hkl, detector-angle, and pixel metadata.
        :param dict refldict:
            Reflection metadata used to describe the saved trajectory.
        :param str name:
            Dataset name stem for saved integrated data.
        :returns:
            Status dictionary describing success, cancellation, or error.
        :rtype: dict

        .. note::
           CLI-capable. Progress reporting is routed through ``logger_utils``.
        """
        logger.info(f"Start rocking integration of scan {name}")
        try:
            image = self.fscan.get_raw_img(0)
        except Exception:
            logger.exception(
                "No image found in current scan.",
                extra={
                    "title": "No image found in current scan",
                    "description": "Cannot integrate scan: No image found in current scan.",  # noqa: E501
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {
                "status": "error",
                "message": "No image found in current scan",
                "traceback": traceback.format_exc(),
            }
        if self.database.nxfile is None:
            logger.exception(
                "Cannot integrate scan: No database available.",
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Cannot integrate scan: No database available.",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return {"status": "error", "message": "No database available"}
        dc = self.ubcalc.detectorCal

        imgmask = None

        if self.scanSelector.useMaskBox.isChecked():
            if self.centralPlot.getMaskToolsDockWidget().getSelectionMask() is None:
                if logger_utils.get_logging_context() == "gui":
                    btn = qt.QMessageBox.question(
                        self,
                        "No mask available",
                        """No mask was selected with the masking tool.
        Do you want to continue without mask?""",
                    )
                    if btn != qt.QMessageBox.Yes:
                        return {
                            "status": "cancelled",
                            "message": "Reason: no mask selected",
                        }
                logger.warn("No mask was selected with the masking tool.")
            else:
                imgmask = (
                    self.centralPlot.getMaskToolsDockWidget().getSelectionMask() > 0.0
                )

        corr = (
            self.scanSelector.useSolidAngleBox.isChecked()
            or self.scanSelector.usePolarizationBox.isChecked()
        )

        C_arr = np.ones(dc.detector.shape, dtype=np.float64)
        if self.scanSelector.useSolidAngleBox.isChecked():
            C_arr /= dc.solidAngleArray()
        if self.scanSelector.usePolarizationBox.isChecked():
            C_arr /= dc.polarization(factor=dc._polFactor, axis_offset=dc._polAxis)

        def fill_counters(image, pixelavail, key, bkgkey):
            """CLI-safe: sum one center ROI and its background ROIs."""

            cimg = image[key[::-1]]

            # !!!!!!!!!! add mask here  !!!!!!!!!
            croi = np.nansum(cimg)
            cpixel = np.nansum(pixelavail[key[::-1]])
            bgroi = 0.0
            bgpixel = 0.0
            for bg in bkgkey:
                bgimg = image[bg[::-1]]
                bgroi += np.nansum(bgimg)
                bgpixel += np.nansum(pixelavail[bg[::-1]])

            return (croi, cpixel, bgroi, bgpixel)

        hkl_del_gam_1 = hkl_del_gam[0]  # needed to initialize integration

        # initialize 1d np arrays for storing roi integration counters for all images
        croi1_a = np.zeros_like(hkl_del_gam_1.shape[0], dtype=np.float64)
        cpixel1_a = np.zeros_like(hkl_del_gam_1.shape[0], dtype=np.float64)
        bgroi1_a = np.zeros_like(hkl_del_gam_1.shape[0], dtype=np.float64)
        bgpixel1_a = np.zeros_like(hkl_del_gam_1.shape[0], dtype=np.float64)

        # initialize 2d np array to store roi integration counters together for all images/ROIs  # noqa: E501
        croi1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        cpixel1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        bgroi1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        bgpixel1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        Corr_croi1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        Corr_cpixel1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        Corr_bgroi1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        Corr_bgpixel1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )

        bgimg_croi1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        bgimg_cpixel1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        bgimg_bgroi1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )
        bgimg_bgpixel1_all = np.zeros(
            (hkl_del_gam_1.shape[0],) + (xylist.shape[0],), dtype=np.float64
        )

        progress = logger_utils.create_progress_logger(
            self, len(self.fscan), "Integrating images"
        )

        background_image = self.background_image
        has_bg_img = False
        roioptions = self.scanSelector.roioptions.get_parameters()
        use_fitted_background = bool(roioptions.get("FittedBackground", False))
        fitted_background_order = int(roioptions.get("FittedBackgroundOrder", 1))
        if use_fitted_background and not HAS_ACCEL:
            logger.warning(
                "Fitted local background requires the compiled ROI accelerator; "
                "using summed background ROIs instead."
            )
        if HAS_ACCEL:
            if imgmask is not None:
                mask = np.ascontiguousarray(imgmask, dtype=bool)
            else:
                mask = np.zeros(image.img.shape, dtype=bool)
            if corr:
                C_arr = np.ascontiguousarray(C_arr, dtype=np.float64)
            else:
                C_arr = np.ones(image.img.shape, dtype=np.float64)
            C_arr[mask] = np.nan

            roi_lists_accel = []
            for roiname in ["center", "left", "right", "top", "bottom"]:
                roi_list = []
                for r in rois[roiname]:
                    roi_list.append(
                        np.array([[r[0].start, r[0].stop], [r[1].start, r[1].stop]])
                    )
                roi_list = np.ascontiguousarray(np.stack(roi_list), dtype=np.int64)
                roi_lists_accel.append(roi_list)
            if (
                background_image is not None
                and background_image.shape == image.img.shape
            ):
                if use_fitted_background:
                    logger.warning(
                        "Fitted local background is ignored when a background "
                        "image is selected."
                    )
                bg = background_image.astype(np.float64, order="C", copy=True)
                bg[mask] = np.nan
                has_bg_img = True

                def sumImage(i):
                    """CLI-safe worker: read and integrate one image with background."""
                    image = self.fscan.get_raw_img(i).img.astype(
                        np.float64, order="C", copy=True
                    )  # unlocks gil during file read

                    all_counters = np.zeros(
                        (roi_lists_accel[0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    Carr_counters = np.zeros(
                        (roi_lists_accel[0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    BgImg_counters = np.zeros(
                        (roi_lists_accel[0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    _roi_sum_accel.processImage_bg_Carr(
                        image,
                        bg,
                        mask,
                        C_arr,
                        *roi_lists_accel,
                        all_counters,
                        Carr_counters,
                        BgImg_counters,
                    )  # compiled accelerator releases the GIL
                    return all_counters, Carr_counters, BgImg_counters
            else:

                def sumImage(i):
                    """CLI-safe worker: read and integrate one image."""
                    image = self.fscan.get_raw_img(i).img.astype(
                        np.float64, order="C", copy=True
                    )  # unlocks gil during file read

                    Carr_counters = np.zeros(
                        (roi_lists_accel[0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    all_counters = np.zeros(
                        (roi_lists_accel[0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    if use_fitted_background:
                        _roi_sum_accel.processImage_polybg_Carr(
                            image,
                            mask,
                            C_arr,
                            *roi_lists_accel,
                            all_counters,
                            Carr_counters,
                            fitted_background_order,
                        )  # compiled accelerator releases the GIL
                    else:
                        _roi_sum_accel.processImage_Carr(
                            image,
                            mask,
                            C_arr,
                            *roi_lists_accel,
                            all_counters,
                            Carr_counters,
                        )  # compiled accelerator releases the GIL
                    return all_counters, Carr_counters

        else:

            def sumImage(i):
                """CLI-safe worker: read and integrate one image without acceleration."""  # noqa: E501
                image = self.fscan.get_raw_img(i).img.astype(
                    np.float64, order="C", copy=True
                )
                if (
                    background_image is not None
                    and background_image.shape == image.shape
                ):
                    np.subtract(image, background_image, out=image)
                if imgmask is not None:
                    image[imgmask] = np.nan
                    pixelavail = (~imgmask).astype(np.float64)
                else:
                    pixelavail = np.ones_like(image)
                if corr:
                    image *= C_arr

                all_counters1 = np.zeros((xylist.shape[0],) + (4,))

                for crnr in range(xylist.shape[0]):
                    # set ROI (moved to rocking-function)

                    # get roi
                    key = rois["center"][crnr]
                    bgkey = [
                        rois["left"][crnr],
                        rois["right"][crnr],
                        rois["top"][crnr],
                        rois["bottom"][crnr],
                    ]
                    # fill counters
                    counters1 = fill_counters(image, pixelavail, key, bgkey)

                    all_counters1[crnr] = counters1

                return all_counters1

        cancelled = False
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.numberthreads
        ) as executor:
            futures = {}
            for i in range(len(self.fscan)):
                futures[executor.submit(sumImage, i)] = i

            status = "no error"
            for f in concurrent.futures.as_completed(futures):  # iteration over jobs
                try:
                    i = futures[f]
                    if has_bg_img:
                        img_counters, Carr_counters, BgImg_counters = f.result()
                        bgimg_croi1_all[i] = BgImg_counters.T[0]
                        bgimg_cpixel1_all[i] = BgImg_counters.T[1]
                        bgimg_bgroi1_all[i] = BgImg_counters.T[2]
                        bgimg_bgpixel1_all[i] = BgImg_counters.T[3]
                    else:
                        img_counters, Carr_counters = f.result()

                    croi1_all[i] = img_counters.T[0]
                    cpixel1_all[i] = img_counters.T[1]
                    bgroi1_all[i] = img_counters.T[2]
                    bgpixel1_all[i] = img_counters.T[3]

                    Corr_croi1_all[i] = Carr_counters.T[0]
                    Corr_cpixel1_all[i] = Carr_counters.T[1]
                    Corr_bgroi1_all[i] = Carr_counters.T[2]
                    Corr_bgpixel1_all[i] = Carr_counters.T[3]
                    # for j in range(len(img_counters)): # iteration over ROIs
                    #    (croi1, cpixel1, bgroi1, bgpixel1) = img_counters[j]
                    #
                    #    croi1_all[i][j] = croi1
                    #    cpixel1_all[i][j] = cpixel1
                    #    bgroi1_all[i][j] = bgroi1
                    #    bgpixel1_all[i][j] = bgpixel1
                    progress.update(futures[f])
                    del f
                except concurrent.futures.CancelledError:
                    pass
                except Exception:
                    logger.warn(
                        f"Cannot read image, cancel to avoid memory leak:\n{traceback.format_exc()}"  # noqa: E501
                    )
                    [f.cancel() for f in futures]
                    cancelled = True
                    status = "error"
                    exc_info = sys.exc_info()
                if progress.wasCanceled():
                    cancelled = True
                    [f.cancel() for f in futures]
                    break

        progress.finish()
        if cancelled:
            if status == "error":
                trace = "".join(traceback.format_exception(*exc_info))
                logger.error(
                    "Error during integration",
                    exc_info=exc_info,
                    extra={
                        "title": "Error during integration",
                        "description": "Error during integration. Integration was aborted.",  # noqa: E501
                        "show_dialog": True,
                        "dialog_level": logging.ERROR,
                        "parent": self,
                    },
                )
                return {
                    "status": "error",
                    "message": "Error during integration",
                    "traceback": trace,
                }
            else:
                return {
                    "status": "cancelled",
                    "message": "Reason: Cancelled during integration",
                }

        currentPlotCount = len(self.integrdataPlot.getAllCurves())
        numberOfNewPlots = xylist.shape[0]
        maxAmountOfPlots = 30
        plotOnlyNth = (
            numberOfNewPlots // max((maxAmountOfPlots - currentPlotCount), 1)
        ) + 1

        # print('Number of integration curves: ' + str(numberOfPlots))
        # print('We can plot every ' + str(plotOnlyNth) + '-th curve.' )

        suffix = ""
        i = 0
        while (
            self.activescanname + "/measurement/" + name + suffix
            in self.database.nxfile
        ):
            suffix = f"_{i}"
            i += 1
        name = name + suffix

        auxcounters = {"@NX_class": "NXcollection"}
        for auxname in self.fscan.auxillary_counters:
            if hasattr(self.fscan, auxname):
                cntr = getattr(self.fscan, auxname)
                if cntr is not None:
                    auxcounters[auxname] = cntr

        if hasattr(self.fscan, "title"):
            title = str(self.fscan.title)
        else:
            title = f"{self.fscan.axisname}-scan"

        mu, om = self.getMuOm()
        if len(np.asarray(om).shape) == 0:
            om = np.full_like(mu, om)
        if len(np.asarray(mu).shape) == 0:
            mu = np.full_like(om, mu)

        data = {
            self.activescanname: {  # legacy, to be removed!
                "instrument": {
                    "@NX_class": "NXinstrument",
                    "positioners": {
                        "@NX_class": "NXcollection",
                        self.fscan.axisname: self.fscan.axis,
                    },
                },
                "auxillary": auxcounters,
                "measurement": {
                    "@NX_class": "NXentry",
                    "@default": name,
                    name: {
                        "@NX_class": "NXentry",
                        "@default": "rois",
                        "@orgui_meta": "rocking",
                        "rois": {
                            "@NX_class": "NXcollection",
                            "@default": None,
                            "@orgui_meta": "roi rocking",
                        },
                    },
                },
                "title": f"{title}",
                "@NX_class": "NXentry",
                "@default": f"measurement/{name}",
                "@orgui_meta": "scan",
            }
        }

        croibg1_bgimg_a = None
        croibg1_bgimg_err_a = None

        # plot and save data in database
        for d in range(croi1_all.shape[1]):
            roi_d = rois["center"][d]
            roi_size = (roi_d[0].stop - roi_d[0].start) * (
                roi_d[1].stop - roi_d[1].start
            )

            hkl_del_gam_1 = hkl_del_gam[d]

            croi1_a = croi1_all[..., d]
            cpixel1_a = cpixel1_all[..., d]
            bgroi1_a = bgroi1_all[..., d]
            bgpixel1_a = bgpixel1_all[..., d]

            Corr_croi1_a = Corr_croi1_all[..., d]
            Corr_cpixel1_a = Corr_cpixel1_all[..., d]
            Corr_bgroi1_a = Corr_bgroi1_all[..., d]
            Corr_bgpixel1_all[..., d]

            bgimg_croi1_a = bgimg_croi1_all[..., d]
            bgimg_cpixel1_a = bgimg_cpixel1_all[..., d]
            bgimg_bgroi1_a = bgimg_bgroi1_all[..., d]
            bgimg_bgpixel1_a = bgimg_bgpixel1_all[..., d]

            Corr1 = Corr_croi1_a * (
                roi_size / Corr_cpixel1_a
            )  # normalize to number of pixels of center roi (croi)

            if np.any(
                bgimg_cpixel1_a
            ):  # assume the background image has no errors (would need a separate error image for that)  # noqa: E501
                bgimg_croi1_norm = bgimg_croi1_a * (cpixel1_a / bgimg_cpixel1_a)
                if np.any(bgpixel1_a):
                    bgimg_bgroi1_norm = bgimg_bgroi1_a * (bgpixel1_a / bgimg_bgpixel1_a)

                    # method 1: simply subtract bg image from data and then subtract the remaining background  # noqa: E501
                    croibg1_a = (
                        (croi1_a - bgimg_croi1_norm)
                        - (cpixel1_a / bgpixel1_a) * (bgroi1_a - bgimg_bgroi1_norm)
                    ) * (roi_size / cpixel1_a)
                    croibg1_err_a = np.sqrt(
                        croi1_a + ((cpixel1_a / bgpixel1_a) ** 2) * bgroi1_a
                    ) * (roi_size / cpixel1_a)

                    # method 2: scale bg image croi and subtract scaled bg image croi. Use ratio of bgroi of image and bg image as scale factor.  # noqa: E501
                    factor = bgroi1_a / bgimg_bgroi1_norm
                    croibg1_bgimg_a = (croi1_a - factor * bgimg_croi1_norm) * (
                        roi_size / cpixel1_a
                    )
                    croibg1_bgimg_err_a = np.sqrt(
                        croi1_a + ((cpixel1_a / bgpixel1_a) ** 2) * bgroi1_a
                    ) * (roi_size / cpixel1_a)

                else:  # not possible if no bgroi is set.
                    croibg1_a = (croi1_a - bgimg_croi1_norm) * (roi_size / cpixel1_a)
                    croibg1_err_a = np.sqrt(croi1_a) * (roi_size / cpixel1_a)

            else:  # no background image
                if np.any(bgpixel1_a):
                    croibg1_a = (croi1_a - (cpixel1_a / bgpixel1_a) * bgroi1_a) * (
                        roi_size / cpixel1_a
                    )
                    croibg1_err_a = np.sqrt(
                        croi1_a + ((cpixel1_a / bgpixel1_a) ** 2) * bgroi1_a
                    ) * (roi_size / cpixel1_a)
                else:
                    croibg1_a = croi1_a * (roi_size / cpixel1_a)
                    croibg1_err_a = np.sqrt(croi1_a) * (roi_size / cpixel1_a)

            if corr:
                croibg1_a *= Corr1
                croibg1_err_a *= Corr1
                if croibg1_bgimg_a is not None:
                    croibg1_bgimg_a *= Corr1
                    croibg1_bgimg_err_a *= Corr1

            rod_mask1 = np.isfinite(croibg1_a)

            axis_masked = hkl_del_gam_1[:, 5][rod_mask1]

            croibg1_a_masked = croibg1_a[rod_mask1]

            croibg1_err_a_masked = croibg1_err_a[rod_mask1]

            # save data

            x, y = xylist[d]
            name1 = f"rocking_{d}"
            if "angles" in refldict:
                alpha1, delta1, gamma1, omega1, chi1, phi1 = refldict["angles"][d]
                sixc_angles_hkl = {
                    "@NX_class": "NXpositioner",
                    "alpha": np.rad2deg(alpha1),
                    "omega": np.rad2deg(omega1),
                    "theta": np.rad2deg(-1 * omega1),
                    "delta": np.rad2deg(delta1),
                    "gamma": np.rad2deg(gamma1),
                    "chi": np.rad2deg(chi1),
                    "phi": np.rad2deg(phi1),
                    "@unit": "deg",
                }
                traj1 = {
                    # "@direction" : u"Rocking scan at fixed pixel location along H_1*s + H_0 in reciprocal space",  # noqa: E501
                    "@NX_class": "NXcollection",
                    "axis": hkl_del_gam_1[:, 5],
                    "HKL_sixc_angles": sixc_angles_hkl,
                }
                # determine the type of rocking scan:
                if "H_1" in refldict:  # H_1 * s H_0 -like rocking scan (CTR scan)
                    traj1["s"] = refldict["s_masked"][d]
                    traj1["H_1"] = refldict["H_1"]
                    traj1["H_0"] = refldict["H_0"]
                    # equal refldict['hkl_masked']?
                    traj1["HKL_pk"] = (
                        refldict["H_1"] * refldict["s_masked"][d] + refldict["H_0"]
                    )
                elif "s_masked" in refldict:
                    traj1["s"] = refldict["s_masked"][d]
                    traj1["HKL_pk"] = refldict["hkl_masked"][d]
            else:
                traj1 = {
                    # "@direction" : u"Rocking scan at fixed pixel location along H_1*s + H_0 in reciprocal space",  # noqa: E501
                    "@NX_class": "NXcollection",
                    "axis": hkl_del_gam_1[:, 5],
                }

            suffix = ""
            i = 0

            while (
                self.activescanname + "/measurement/" + name + "/" + name1 + suffix
                in self.database.nxfile
            ):
                suffix = f"_{i}"
                i += 1

            availname1 = name1 + suffix

            x, y = xylist[d]  #
            # x_coord1_a = xylist[:,0]
            # y_coord1_a = xylist[:,1]

            datas1 = {
                "@NX_class": "NXdata",
                "sixc_angles": {
                    "@NX_class": "NXpositioner",
                    "alpha": np.rad2deg(mu),
                    "omega": np.rad2deg(om),
                    "theta": np.rad2deg(-1 * om),
                    "delta": np.rad2deg(hkl_del_gam_1[:, 3]),
                    "gamma": np.rad2deg(hkl_del_gam_1[:, 4]),
                    "chi": np.rad2deg(self.ubcalc.chi),
                    "phi": np.rad2deg(self.ubcalc.phi),
                    "@unit": "deg",
                },
                "hkl": {
                    "@NX_class": "NXcollection",
                    "h": hkl_del_gam_1[:, 0],
                    "k": hkl_del_gam_1[:, 1],
                    "l": hkl_del_gam_1[:, 2],
                },
                "counters": {
                    "@NX_class": "NXdetector",
                    "croibg": croibg1_a,
                    "croibg_errors": croibg1_err_a,
                    "croibg_bgimg": croibg1_bgimg_a,  # when None, will not create data set  # noqa: E501
                    "croibg_bgimg_errors": croibg1_bgimg_err_a,  # when None, will not create data set  # noqa: E501
                    "croi": croi1_a,
                    "bgroi": bgroi1_a,
                    "croi_pix": cpixel1_a,
                    "bgroi_pix": bgpixel1_a,
                    "Cfactors_croi": Corr_croi1_a,
                    "Cfactors_bgroi": Corr_bgroi1_a,
                    "bgimg_croi": bgimg_croi1_a,
                    "bgimg_bgroi": bgimg_bgroi1_a,
                },
                "pixelcoord": {
                    "@NX_class": "NXdetector",
                    "x": x,
                    "y": y,
                    "vsize": (roi_d[1].stop - roi_d[1].start),
                    "hsize": (roi_d[0].stop - roi_d[0].start),
                },
                "trajectory": traj1,
                "@signal": "counters/croibg",
                "@axes": "trajectory/axis",
                "@title": self.activescanname + "_" + availname1,
                "@orgui_meta": "roi rocking",
            }

            data[self.activescanname]["measurement"][name]["rois"]["@default"] = (
                availname1
            )
            if np.any(cpixel1_a > 0.0):
                data[self.activescanname]["measurement"][name]["rois"][availname1] = (
                    datas1
                )
                if d % plotOnlyNth == 0 and not min(croibg1_a_masked) == max(
                    croibg1_a_masked
                ):
                    self.integrdataPlot.addCurve(
                        axis_masked,
                        croibg1_a_masked,
                        legend=self.activescanname + "_" + availname1,
                        xlabel=f"trajectory/{self.fscan.axisname}",
                        ylabel="counters/croibg",
                        yerror=croibg1_err_a_masked,
                    )

        # lets keep legacy data structure for now

        data_2d_structured = {
            self.activescanname: {
                "instrument": {
                    "@NX_class": "NXinstrument",
                    "positioners": {
                        "@NX_class": "NXcollection",
                        self.fscan.axisname: self.fscan.axis,
                    },
                },
                "auxillary": auxcounters,
                "measurement": {
                    "@NX_class": "NXentry",
                    "@default": name,
                    name: {
                        "@NX_class": "NXentry",
                        "@default": "rois",
                        "@orgui_meta": "rocking",
                    },
                },
                "title": f"{title}",
                "@NX_class": "NXentry",
                "@default": f"measurement/{name}",
                "@orgui_meta": "scan",
            }
        }
        alpha = []
        theta = []
        delta = []
        gamma = []
        chi = []
        phi = []
        omega = []
        alpha_pk = []
        theta_pk = []
        delta_pk = []
        gamma_pk = []
        chi_pk = []
        phi_pk = []
        omega_pk = []
        x = []
        y = []
        h = []
        k = []
        l = []  # noqa: E741
        croibg = []
        croibg_errors = []
        croi = []
        bgroi = []
        croi_pix = []
        bgroi_pix = []
        croibg_bgimg = []
        croibg_bgimg_errors = []
        Cfactors_croi = []
        Cfactors_bgroi = []
        bgimg_croi = []
        bgimg_bgroi = []
        axis = []
        s = []
        H_0 = []
        H_1 = []
        HKL_pk = []
        vsize = []
        hsize = []

        # from IPython import embed; embed()

        optional_labels = {"s": s, "H_1": H_1, "H_0": H_0, "HKL_pk": HKL_pk}

        for sc in data[self.activescanname]["measurement"][name]["rois"]:
            if sc.startswith("@"):
                continue
            try:
                dsc = data[self.activescanname]["measurement"][name]["rois"][sc]

                # 2D arrays
                alpha.append(dsc["sixc_angles"]["alpha"])
                theta.append(dsc["sixc_angles"]["theta"])
                delta.append(dsc["sixc_angles"]["delta"])
                gamma.append(dsc["sixc_angles"]["gamma"])
                chi.append(dsc["sixc_angles"]["chi"])
                phi.append(dsc["sixc_angles"]["phi"])
                omega.append(dsc["sixc_angles"]["omega"])

                # 2D arrays
                h.append(dsc["hkl"]["h"])
                k.append(dsc["hkl"]["k"])
                l.append(dsc["hkl"]["l"])

                # 2D arrays
                croibg.append(dsc["counters"]["croibg"])
                croibg_errors.append(dsc["counters"]["croibg_errors"])
                croi.append(dsc["counters"]["croi"])
                bgroi.append(dsc["counters"]["bgroi"])
                croi_pix.append(dsc["counters"]["croi_pix"])
                bgroi_pix.append(dsc["counters"]["bgroi_pix"])
                if dsc["counters"]["croibg_bgimg"] is not None:
                    croibg_bgimg.append(dsc["counters"]["croibg_bgimg"])
                    croibg_bgimg_errors.append(dsc["counters"]["croibg_bgimg_errors"])
                Cfactors_croi.append(dsc["counters"]["Cfactors_croi"])
                Cfactors_bgroi.append(dsc["counters"]["Cfactors_bgroi"])
                bgimg_croi.append(dsc["counters"]["bgimg_croi"])
                bgimg_bgroi.append(dsc["counters"]["bgimg_bgroi"])

                # 1D arrays
                x.append(dsc["pixelcoord"]["x"])
                y.append(dsc["pixelcoord"]["y"])

                # 1D arrays
                vsize.append(dsc["pixelcoord"]["vsize"])
                hsize.append(dsc["pixelcoord"]["hsize"])

                axis.append(dsc["trajectory"]["axis"])

                for lbl in optional_labels:
                    if lbl in dsc["trajectory"]:
                        optional_labels[lbl].append(dsc["trajectory"][lbl])

                # 1d Array
                if "HKL_sixc_angles" in dsc["trajectory"]:
                    alpha_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["alpha"])
                    theta_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["theta"])
                    delta_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["delta"])
                    gamma_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["gamma"])
                    chi_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["chi"])
                    phi_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["phi"])
                    omega_pk.append(dsc["trajectory"]["HKL_sixc_angles"]["omega"])
            except Exception:
                logger.exception(
                    "Unexpected exception while creating data sets to save"
                )
                # from IPython import embed; embed()
                # sys.exit(0)

        rois = {
            "@NX_class": "NXcollection",
            "@default": "croibg",
            "@orgui_meta": "roi rocking",
            "alpha": np.vstack(alpha),
            "theta": np.vstack(theta),
            "delta": np.vstack(delta),
            "gamma": np.vstack(gamma),
            "chi": np.vstack(chi),
            "phi": np.vstack(phi),
            "omega": np.vstack(omega),
            "h": np.vstack(h),
            "k": np.vstack(k),
            "l": np.vstack(l),
            "croibg": np.vstack(croibg),
            "croibg_errors": np.vstack(croibg_errors),
            "croi": np.vstack(croi),
            "bgroi": np.vstack(bgroi),
            "croi_pix": np.vstack(croi_pix),
            "bgroi_pix": np.vstack(bgroi_pix),
            "Cfactors_croi": np.vstack(Cfactors_croi),
            "Cfactors_bgroi": np.vstack(Cfactors_bgroi),
            "bgimg_croi": np.vstack(bgimg_croi),
            "bgimg_bgroi": np.vstack(bgimg_bgroi),
            "x": np.array(x),
            "y": np.array(y),
            "vsize": np.array(vsize),
            "hsize": np.array(hsize),
            "axis": np.vstack(axis),
        }
        if alpha_pk:
            rois["alpha_pk"] = np.array(alpha_pk)
            rois["theta_pk"] = np.array(theta_pk)
            rois["delta_pk"] = np.array(delta_pk)
            rois["gamma_pk"] = np.array(gamma_pk)
            rois["chi_pk"] = np.array(chi_pk)
            rois["phi_pk"] = np.array(phi_pk)
            rois["omega_pk"] = np.array(omega_pk)

        if croibg_bgimg:
            rois["croibg_bgimg"] = np.vstack(croibg_bgimg)
            rois["croibg_bgimg_errors"] = np.vstack(croibg_bgimg_errors)

        for lbl in optional_labels:
            if optional_labels[lbl]:
                rois[lbl] = np.squeeze(np.vstack(optional_labels[lbl]))

        scsize = rois["axis"].shape[0]
        for t in rois:
            if t.startswith("@"):
                continue
            if rois[t].shape[0] != scsize:
                logger.error(
                    "Error during ro integration: roi %s does not match scan size %s. "
                    "This is likely a coding error",
                    t,
                    scsize,
                )
                return {
                    "status": "error",
                    "message": "Error during ro integration: size mismatch",
                    "traceback": "",
                }

        data_2d_structured[self.activescanname]["measurement"][name]["rois"] = rois

        self.database.add_nxdict(data_2d_structured)
        logger.info(f"Rocking integration succeeded and data saved with name {name}")
        return {"status": "success"}

    def updatePlotItems(self, recalculate=True):
        """Refresh displayed ROI and reflection overlays.

        :param bool recalculate:
            Recalculate CTR coordinates before updating reflection markers.

        .. note::
           CLI-capable, but it mutates Qt plot widgets.
        """
        if self.roivisible:
            try:
                self.updateROI()
            except Exception:
                pass

        if self.reflectionsVisible:
            if recalculate:
                try:
                    hkm = self.calculateAvailableCTR()
                    hk = np.unique(hkm[:, :2], axis=0)
                    H_0 = np.hstack((hk, np.zeros((hk.shape[0], 1))))
                    H_1 = np.array([0, 0, 1])
                    self.reflectionsToDisplay = H_0, H_1
                except Exception:
                    self.showCTRreflAct.setChecked(False)
                    self.reflectionsVisible = False
            self.updateReflections()

        if self.reflectionSel.showBraggReflections:
            try:
                self.calcBraggRefl()
            except Exception:
                pass

    def _onSetBackgroundImage(self, checked):
        """GUI-only: select or clear a background image via a file dialog."""
        if not checked:
            self.background_image = None
            self.plotImage(self.imageno)  # will not raise Exception,
            return
        extensions = {}
        for description, ext in silx.io.supported_extensions().items():
            extensions[description] = " ".join(sorted(list(ext)))

        extensions["NeXus layout from EDF files"] = "*.edf"
        extensions["NeXus layout from TIFF image files"] = "*.tif *.tiff"
        extensions["NeXus layout from CBF files"] = "*.cbf"
        extensions["NeXus layout from MarCCD image files"] = "*.mccd"

        all_supported_extensions = set()
        for name, exts in extensions.items():
            exts = exts.split(" ")
            all_supported_extensions.update(exts)
        all_supported_extensions = sorted(list(all_supported_extensions))

        filters = []
        filters.append(
            "All supported files ({})".format(" ".join(all_supported_extensions))
        )
        for name, extension in extensions.items():
            filters.append(f"{name} ({extension})")
        filters.append("All files (*)")

        fileTypeFilter = ""
        for f in filters:
            fileTypeFilter += f + ";;"

        # call dialog
        filename, _ = qt.QFileDialog.getOpenFileName(
            self, "Open background image", "", fileTypeFilter[:-2]
        )
        if filename == "":
            self.backgroundImageAct.setChecked(False)
            self.plotImage(self.imageno)  # will not raise Exception,
            return
        try:
            with fabio.open(filename) as fabf:
                self.background_image = fabf.data.astype(
                    np.float64, order="C", copy=True
                )
            self.plotImage(self.imageno)  # will not raise Exception,

        except Exception:
            traceback.print_exc()
            self.backgroundImageAct.setChecked(False)
            self.plotImage(self.imageno)  # will not raise Exception,
            return

        # dialog = ImageFileDialog.ImageFileDialog(self)
        # result = dialog.exec()
        # if result:
        #    self.background_image = dialog.selectedImage().astype(np.float64)

    def onShowBragg(self, visible):
        """GUI/CLI hint: toggle Bragg reflection overlays on the plot."""
        try:
            self.reflectionSel.setBraggReflectionsVisible(visible)
            if visible:
                self.calcBraggRefl()
        except Exception:
            logger.exception(
                "Cannot show Bragg reflections",
                extra={
                    "title": "Cannot show Bragg reflections",
                    "description": "Cannot show Bragg reflections",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            # qutils.warning_detailed_message(self, "Cannot show show Bragg reflections", "Cannot show Bragg reflections", traceback.format_exc())  # noqa: E501
            # qt.QMessageBox.critical(self,"Cannot show show Bragg reflections", "Cannot Cannot show Bragg reflections:\n%s" % traceback.format_exc())  # noqa: E501

    def onShowROI(self, visible):
        """GUI/CLI hint: toggle ROI display and refresh ROI graphics."""
        self.roivisible = visible
        self.showInterpolatedBgAct.setEnabled(visible and HAS_ACCEL)
        try:
            if self.showInterpolatedBgAct.isChecked():
                self.plotImage(self.imageno)
                self.updateROI()
                return
            self.updateROI()
        except Exception:
            logger.exception(
                "Cannot show ROI",
                extra={
                    "title": "Cannot show ROI",
                    "description": "Cannot show ROI",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            # qutils.warning_detailed_message(self, "Cannot show ROI", "Cannot show ROI", traceback.format_exc())  # noqa: E501
            # qt.QMessageBox.critical(self,"Cannot show ROI", "Cannot Cannot show ROI:\n%s" % traceback.format_exc())  # noqa: E501

    def onShowInterpolatedBg(self, visible):
        """GUI/CLI hint: toggle fitted-background preview in center ROIs."""
        try:
            self.plotImage(self.imageno)
            self.updateROI()
        except Exception:
            logger.exception(
                "Cannot show interpolated background",
                extra={
                    "title": "Cannot show interpolated background",
                    "description": "Cannot show interpolated background",
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )

    def onShowCTRreflections(self, visible):
        """GUI/CLI hint: toggle CTR reflection overlays on the plot."""
        self.reflectionsVisible = visible
        if self.reflectionsVisible:
            try:
                hkm = self.calculateAvailableCTR()
            except Exception:
                self.showCTRreflAct.setChecked(False)
                self.reflectionsVisible = False
                logger.exception(
                    "Cannot calculate CTR locations",
                    extra={
                        "title": "Cannot calculate CTR locations",
                        "description": "Cannot calculate CTR locations",
                        "show_dialog": True,
                        "dialog_level": logging.WARNING,
                        "parent": self,
                    },
                )
                # qutils.warning_detailed_message(self, "Cannot calculate CTR locations", "Cannot calculate CTR locatons", traceback.format_exc())  # noqa: E501
                # qt.QMessageBox.critical(self,"Cannot calculate CTR locatons", "Cannot calculate CTR locatons:\n%s" % traceback.format_exc())  # noqa: E501
                return
            hk = np.unique(hkm[:, :2], axis=0)
            H_0 = np.hstack((hk, np.zeros((hk.shape[0], 1))))
            H_1 = np.array([0, 0, 1])

            self.reflectionsToDisplay = H_0, H_1
        self.updateReflections()

    def _onShowAbout(self):
        """GUI-only: show the modal About dialog."""
        dial = AboutDialog(self, __version__)
        dial.exec()

    def _onShowDiffractionGeometry(self):
        """GUI-only: show or create the diffraction-geometry dialog."""
        if hasattr(self, "diffractometerdialog"):
            self.diffractometerdialog.show()
        else:
            self.diffractometerdialog = QDiffractometerImageDialog()
            self.diffractometerdialog.show()

    def _onToggleExcludeImage(self, exclude):
        """CLI-capable: update the excluded-image list for the active image."""
        currentimgno = self.scanSelector.slider.value()
        data = self.excludedImagesDialog.getData()
        imgno_in_excudearray = currentimgno in data

        if imgno_in_excudearray and exclude:
            return
        if not imgno_in_excudearray and exclude:
            data = np.hstack([data, currentimgno])
            self.excludedImagesDialog.updateArrayData(data)
        else:
            data = data[data != currentimgno]
            self.excludedImagesDialog.updateArrayData(data)

    def _onSelectCPUcount(self):
        """GUI-only: ask the user for the worker-thread count."""
        maxavail = os.cpu_count() if os.cpu_count() is not None else 1
        if "SLURM_CPUS_ON_NODE" in os.environ:
            maxavail = int(os.environ["SLURM_CPUS_ON_NODE"])

        cpus, success = qt.QInputDialog.getInt(
            self,
            "CPU count",
            f"CPU count (detected: {maxavail})",
            self.numberthreads,
            1,
        )
        if success:
            self.numberthreads = cpus

    def _onSelectROIcount(self):
        """GUI-only: ask the user for the maximum displayed ROI count."""

        rois, success = qt.QInputDialog.getInt(
            self,
            "ROI count",
            f"Max ROI count to display (current: {self.maxROIs})",
            self.maxROIs,
            1,
        )
        if success:
            self.maxROIs = rois

    def _onChangeDBCompression(self):
        """GUI-only: ask the user for the database compression filter."""
        filter_names = list(FILTERS.keys())
        currentCompression = self.database.compression
        for fn in filter_names:
            if currentCompression == FILTERS[fn]:
                break
        idx = filter_names.index(fn)
        selection, success = qt.QInputDialog.getItem(
            self,
            "Data compression settings",
            "Available data base compression methods:\n(See discussion under https://github.com/tifuchs/orGUI/issues/16)\nRecommended: Blosc-lz4-Shuffle-5",  # noqa: E501
            filter_names,
            idx,
            False,
        )
        if success:
            self.database.compression = FILTERS[selection]

    def calcBraggRefl(self):
        """Calculate and display available Bragg reflections for the scan.

        .. note::
           CLI-capable for ``th`` scans, but it updates GUI reflection state.
        """
        if self.fscan is not None:
            if self.fscan.axisname != "th":
                raise NotImplementedError(
                    f"Calculation of available Bragg reflections is not implemented for {self.fscan.axisname} - scans"  # noqa: E501
                )
            try:
                xtal = self.ubcalc.crystal
                ommin = np.deg2rad(np.amin(self.fscan.omega))
                ommax = np.deg2rad(np.amax(self.fscan.omega))
                dc = self.ubcalc.detectorCal
                mu = self.ubcalc.mu
                ub = self.ubcalc.ubCal
                chi = self.ubcalc.chi
                phi = self.ubcalc.phi
                xtal.setEnergy(ub.getEnergy() * 1e3)
                hkls, yx, angles = rn.thscanBragg(
                    xtal, ub, mu, dc, (ommin, ommax), chi=chi, phi=phi
                )
                self.reflectionSel.setBraggReflections(hkls, yx, angles)
            except Exception:
                logger.exception(
                    "Cannot calculate Bragg reflections",
                    extra={
                        "title": "Cannot calculate Bragg reflections",
                        "description": "Cannot calculate Bragg reflections",
                        "show_dialog": True,
                        "dialog_level": logging.WARNING,
                        "parent": self,
                    },
                )
                # qutils.warning_detailed_message(self, "Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections", traceback.format_exc())  # noqa: E501
                # qt.QMessageBox.critical(self,"Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections:\n%s" % traceback.format_exc())  # noqa: E501

    def autoFindBraggReference(self, **kwargs):
        """Search the active scan for a Bragg peak and seed the UB matrix.

        The search streams the max-pixel trace image by image, emits the first
        sharp post-burn-in peak candidate after a short lookahead confirmation,
        and validates that candidate immediately. If refinement or Bragg
        indexing fails, the stream continues from the next unread image until a
        valid seed is found or the scan ends.

        :param int burn_in:
            Number of initial images ignored for sharp-increase detection.
            Defaults to ``5``.
        :param float min_sharpness:
            Minimum robust level z-score for automatic candidates. Defaults to
            ``6``.
        :param float min_derivative_sharpness:
            Minimum robust first-difference z-score for automatic candidates.
            Defaults to ``6``.
        :param int lookahead:
            Number of following images read to confirm a local peak-like max
            trace candidate. Defaults to ``2``.
        :param int max_candidates:
            Maximum number of streamed candidates to validate before giving up.
            Defaults to no explicit limit.
        :param float qnorm_tolerance:
            Allowed Q-norm mismatch in Angstrom^-1 when matching Bragg
            reflections. Defaults to ``0.05``.
        :returns:
            Best UB seed, or ``None`` if no reliable seed was found.
        :rtype: orgui.app.autoBraggSearch.UBSeed or None

        .. note::
           CLI-capable when scan, mask, detector, and crystal state are already
           configured. The method updates the reference-reflection table and U.
        """
        if self.fscan is None:
            raise ValueError("No scan loaded")
        status_callback = kwargs.get("status_callback", None)

        def report(event, **fields):
            if status_callback is not None:
                status_callback(event, **fields)

        burn_in = kwargs.get("burn_in", 5)
        history = kwargs.get("history", 20)
        min_history = kwargs.get("min_history", 5)
        min_sharpness = kwargs.get("min_sharpness", 6.0)
        min_derivative_sharpness = kwargs.get("min_derivative_sharpness", 6.0)
        min_prominence_sharpness = kwargs.get("min_prominence_sharpness", 4.0)
        lookahead = kwargs.get("lookahead", 2)
        refractory = kwargs.get("refractory", 3)
        max_candidates = kwargs.get("max_candidates", None)
        mask_distance = kwargs.get("mask_distance", 3)
        qnorm_tolerance = kwargs.get("qnorm_tolerance", 0.05)
        adaptive_after_candidates = kwargs.get("adaptive_after_candidates", 5)
        adaptive_qnorm_tolerance = kwargs.get(
            "adaptive_qnorm_tolerance", max(0.25, 5.0 * qnorm_tolerance)
        )
        adaptive_assignment_pixel_tolerance = kwargs.get(
            "adaptive_assignment_pixel_tolerance",
            max(100.0, 3.0 * kwargs.get("assignment_pixel_tolerance", 30.0)),
        )
        adaptive_confirmation_pixel_tolerance = kwargs.get(
            "adaptive_confirmation_pixel_tolerance",
            max(
                30.0,
                3.0 * kwargs.get("confirmation_pixel_tolerance", 8.0),
            ),
        )
        adaptive_confirmation_image_tolerance = kwargs.get(
            "adaptive_confirmation_image_tolerance",
            max(8, 2 * kwargs.get("confirmation_image_tolerance", 3)),
        )
        adaptive_scale_detector_filter = kwargs.get(
            "adaptive_scale_detector_filter", True
        )
        adaptive_scale_detector_fraction = kwargs.get(
            "adaptive_scale_detector_fraction", 0.25
        )
        adaptive_scale_outlier_q_tolerance = kwargs.get(
            "adaptive_scale_outlier_q_tolerance",
            adaptive_qnorm_tolerance,
        )
        adaptive_scale_outlier_angle_fraction = kwargs.get(
            "adaptive_scale_outlier_angle_fraction",
            adaptive_scale_detector_fraction,
        )
        hkl_candidate_mode = kwargs.get("hkl_candidate_mode", "detector_position")
        if hkl_candidate_mode == "qnorm":
            logger.warning(
                "Automatic Bragg HKL assignment uses Q-norm search. This "
                "checks more hypotheses and may take longer than detector "
                "position matching."
            )
        elif hkl_candidate_mode != "detector_position":
            raise ValueError(
                "hkl_candidate_mode must be 'detector_position' or 'qnorm'"
            )
        axis_half_width = kwargs.get("axis_half_width", 1.0)
        fine_axis_half_width = kwargs.get("fine_axis_half_width", 0.4)
        roi_size = kwargs.get("roi_size", (80, 80))
        fine_roi_size = kwargs.get("fine_roi_size", (40, 40))

        mask = None
        mask_widget = self.centralPlot.getMaskToolsDockWidget()
        if mask_widget.getSelectionMask() is not None:
            mask = mask_widget.getSelectionMask() > 0.0
        else:
            if logger_utils.get_logging_context() == "gui":
                # GUI-only: user-triggered automatic Bragg search confirmation.
                btn = qt.QMessageBox.question(
                    self,
                    "Automatic Bragg search without mask",
                    "No detector mask is set. Automatic Bragg search will "
                    "continue without masked-pixel rejection or mask-aware "
                    "peak refinement.\n\nDo you want to continue?",
                    qt.QMessageBox.Yes | qt.QMessageBox.No,
                    qt.QMessageBox.No,
                )
                if btn != qt.QMessageBox.Yes:
                    report(
                        "failed",
                        message=(
                            "Automatic Bragg search cancelled because no "
                            "detector mask is set."
                        ),
                    )
                    return None
            logger.warning(
                "Automatic Bragg search is running without a detector mask. "
                "Masked-pixel rejection and mask-aware peak refinement will "
                "be disabled."
            )
            report(
                "warning",
                message=(
                    "No detector mask is set; masked-pixel rejection and "
                    "mask-aware peak refinement are disabled."
                ),
            )

        excluded_images = self.excludedImagesDialog.getData()
        self.ubcalc.crystal.setEnergy(self.ubcalc.ubCal.getEnergy() * 1e3)
        max_q = kwargs.get(
            "max_q",
            getattr(self.ubcalc.detectorCal, "Qmax", None),
        )

        logger.info("Start online automatic Bragg max-pixel search")
        report(
            "start",
            message=(
                f"Started automatic Bragg search: mode={hkl_candidate_mode}, images={len(self.fscan)}"  # noqa: E501
            ),
            images_read=0,
        )
        candidate_count = 0
        refined_peaks = []
        unmatched_peaks = []
        qnorm_scale = 1.0
        adaptive_matching = False
        last_scale_fit_count = 0
        scale_revalidation_done = False

        def scale_fit_detector_filter(peak, hkls, norms):
            """Filter scale-fit HKLs by broad predicted detector position."""
            if (
                not adaptive_scale_detector_filter
                or hkl_candidate_mode != "detector_position"
            ):
                return np.ones(len(hkls), dtype=bool)
            if len(hkls) == 0:
                return np.zeros(0, dtype=bool)
            try:
                refldict = self.ubcalc.calcReflection(np.asarray(hkls, dtype=float))
            except Exception:
                logger.debug(
                    "Cannot calculate vectorized detector-position filter "
                    "for adaptive Q-scale fit.",
                    exc_info=True,
                )
                return np.zeros(len(hkls), dtype=bool)
            detvsize, dethsize = self.ubcalc.detectorCal.detector.shape
            half_width = adaptive_scale_detector_fraction * np.array(
                [dethsize, detvsize], dtype=float
            )
            observed_xy = np.asarray(peak.xy, dtype=float)
            keep = np.zeros(len(hkls), dtype=bool)
            for key in ("xy_1", "xy_2"):
                xy = np.atleast_2d(np.asarray(refldict[key], dtype=float))
                finite = np.all(np.isfinite(xy), axis=1)
                in_region = np.all(
                    np.abs(xy - observed_xy) <= half_width,
                    axis=1,
                )
                keep |= finite & in_region
            return keep

        def filter_scale_fit_matches(qscale_fit):
            """Return provisional scale-fit matches without obvious outliers."""
            if qscale_fit is None:
                return []
            matches = list(qscale_fit.get("matches", ()))
            if not matches:
                return []
            scale = qscale_fit["qnorm_scale"]
            q_residual = np.array(
                [
                    abs(match["configured_norm"] * scale - match["qnorm"])
                    for match in matches
                ]
            )
            q_keep = q_residual <= adaptive_scale_outlier_q_tolerance
            if hkl_candidate_mode != "detector_position":
                for match, residual in zip(matches, q_residual):
                    match["scaled_q_residual"] = float(residual)
                return [match for match, keep in zip(matches, q_keep) if keep]

            hkls = np.asarray([match["hkl"] for match in matches], dtype=float)
            try:
                refldict = self.ubcalc.calcReflection(hkls)
                predicted_1 = np.atleast_2d(
                    np.asarray(refldict["angles_1"], dtype=float)
                )[:, [1, 2]]
                predicted_2 = np.atleast_2d(
                    np.asarray(refldict["angles_2"], dtype=float)
                )[:, [1, 2]]
            except Exception:
                logger.debug(
                    "Cannot calculate vectorized delta/gamma outlier check "
                    "for adaptive Q-scale fit.",
                    exc_info=True,
                )
                predicted_1 = np.full((len(matches), 2), np.inf)
                predicted_2 = np.full((len(matches), 2), np.inf)

            peaks = [match["peak"] for match in matches]
            mus = np.array([self.getMuOm(peak.imageno)[0] for peak in peaks])
            measured_gamma, measured_delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([peak.xy[1] for peak in peaks]),
                np.array([peak.xy[0] for peak in peaks]),
                mus,
            )
            measured = np.column_stack((measured_delta, measured_gamma))
            try:
                gamrange, delrange = self.ubcalc.detectorCal.rangegamdel(
                    float(np.nanmedian(mus))
                )
                angle_half_width = adaptive_scale_outlier_angle_fraction * (
                    np.array(
                        [
                            abs(delrange[1] - delrange[0]),
                            abs(gamrange[1] - gamrange[0]),
                        ],
                        dtype=float,
                    )
                )
            except Exception:
                logger.debug(
                    "Cannot calculate detector delta/gamma range for "
                    "adaptive Q-scale outlier check.",
                    exc_info=True,
                )
                angle_half_width = np.full(2, np.inf, dtype=float)

            residual_1 = np.abs(predicted_1 - measured)
            residual_2 = np.abs(predicted_2 - measured)
            angle_keep_1 = np.all(residual_1 <= angle_half_width, axis=1)
            angle_keep_2 = np.all(residual_2 <= angle_half_width, axis=1)
            delgam_error = np.minimum(
                np.linalg.norm(residual_1, axis=1),
                np.linalg.norm(residual_2, axis=1),
            )
            keep = q_keep & (angle_keep_1 | angle_keep_2)
            for match, qres, dgerr in zip(matches, q_residual, delgam_error):
                match["scaled_q_residual"] = float(qres)
                match["delgam_error"] = float(dgerr)
            return [match for match, keep_match in zip(matches, keep) if keep_match]

        def validate_automatic_peak(
            peak,
            current_qnorm_tolerance,
            current_assignment_pixel_tolerance,
            current_confirmation_pixel_tolerance,
            current_confirmation_image_tolerance,
        ):
            """Try HKL assignment and second-peak confirmation for ``peak``."""
            if hkl_candidate_mode == "qnorm":
                seeds = autoBraggSearch.rank_ub_seeds(
                    [peak],
                    self.ubcalc.detectorCal,
                    self.getMuOm,
                    self.ubcalc.ubCal,
                    self.ubcalc.crystal,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                    qnorm_tolerance=current_qnorm_tolerance,
                    max_q=max_q,
                    qnorm_scale=qnorm_scale,
                )
            else:
                seeds = self._rankAutomaticSeedsByDetectorPosition(
                    peak,
                    q_tolerance=current_qnorm_tolerance,
                    max_q=max_q,
                    pixel_tolerance=current_assignment_pixel_tolerance,
                    max_reflections=kwargs.get("assignment_reflections", 40),
                    qnorm_scale=qnorm_scale,
                )
            if not seeds:
                return None, None, "no HKL hypotheses in local Q region"
            report(
                "hkl_hypotheses",
                count=len(seeds),
                message=(
                    f"Testing {len(seeds)} HKL hypotheses for candidate image {peak.imageno}."  # noqa: E501
                ),
            )

            if hkl_candidate_mode == "detector_position":
                seed_hypotheses = seeds
            else:
                seed_hypotheses = seeds[: kwargs.get("max_seed_hypotheses", 12)]
            for candidate_seed in seed_hypotheses:
                confirmation = self._confirmAutomaticBraggSeed(
                    candidate_seed,
                    mask=mask,
                    excluded_images=excluded_images,
                    max_q=max_q,
                    axis_half_width=axis_half_width,
                    fine_axis_half_width=fine_axis_half_width,
                    roi_size=roi_size,
                    fine_roi_size=fine_roi_size,
                    pixel_tolerance=current_confirmation_pixel_tolerance,
                    image_tolerance=current_confirmation_image_tolerance,
                    max_reflections=kwargs.get("confirmation_reflections", 12),
                    intensity_ratio_check=kwargs.get("intensity_ratio_check", True),
                    prominence_threshold=kwargs.get(
                        "confirmation_prominence_threshold", 6.0
                    ),
                    status_callback=status_callback,
                )
                if confirmation is not None:
                    if np.isfinite(confirmation["predicted_intensity_ratio"]):
                        intensity_message = (
                            "intensity ratio obs/calc={:.3g}/{:.3g}".format(
                                confirmation["observed_intensity_ratio"],
                                confirmation["predicted_intensity_ratio"],
                            )
                        )
                    else:
                        intensity_message = (
                            "prominence z seed/confirmation={:.3g}/{:.3g}".format(
                                confirmation["seed_prominence_z"],
                                confirmation["confirmation_prominence_z"],
                            )
                        )
                    report(
                        "confirmation",
                        message=(
                            "Confirmed hkl={} with second hkl={}; "
                            "pixel error={:.3g}, image error={}, {}".format(
                                candidate_seed.hkl,
                                confirmation["hkl"],
                                confirmation["pixel_error"],
                                confirmation["image_error"],
                                intensity_message,
                            )
                        ),
                    )
                    return candidate_seed, confirmation, None
            return None, None, "no second Bragg peak confirmed the HKL assignment"

        def accept_automatic_seed(seed, confirmation):
            """Apply a confirmed automatic UB seed and add reference peaks."""
            self.ubcalc.ubCal.setU(seed.U)
            self.ubcalc.uedit.setU(seed.U)
            eventdict = {"x": seed.peak.xy[0], "y": seed.peak.xy[1]}
            refl = self.reflectionSel.addReflection(
                eventdict, seed.peak.imageno, seed.hkl
            )
            if confirmation is not None:
                confirm_eventdict = {
                    "x": confirmation["peak"].xy[0],
                    "y": confirmation["peak"].xy[1],
                }
                self.reflectionSel.addReflection(
                    confirm_eventdict,
                    confirmation["peak"].imageno,
                    confirmation["hkl"],
                )
            extra_q_tolerance = (
                adaptive_qnorm_tolerance if adaptive_matching else qnorm_tolerance
            )
            extra_count = self._addAutomaticObservedReflections(
                refined_peaks,
                seed,
                confirmation,
                q_tolerance=extra_q_tolerance,
                score_tolerance_deg=2.0,
                intensity_ratio_check=kwargs.get("intensity_ratio_check", True),
                status_callback=status_callback,
            )
            if extra_count:
                report(
                    "accepted",
                    message=(
                        f"Added {extra_count} additional observed Bragg reference "
                        "reflection(s) after U confirmation."
                    ),
                )
            self.reflectionSel.setReflectionActive(refl.identifier)
            self._onChangeImage(seed.peak.imageno)
            self._onCenterGraph(seed.peak.xy)
            self.ubcalc.updateReflectionMismatch()
            self.ubcalc.sigReplotRequest.emit(False)
            logger.info(
                "Automatic Bragg seed selected hkl=%s at image %s, xy=%s, "
                "score=%.4g, norm mismatch=%.4g Angstrom^-1.",
                seed.hkl,
                seed.peak.imageno,
                seed.peak.xy,
                seed.score,
                seed.norm_mismatch,
            )
            report(
                "success",
                message=(
                    f"Selected hkl={seed.hkl} at image={seed.peak.imageno}, xy=({seed.peak.xy[0]:.2f}, {seed.peak.xy[1]:.2f}), "  # noqa: E501
                    f"norm mismatch={seed.norm_mismatch:.4g}"
                ),
            )
            return seed

        candidates = autoBraggSearch.iter_sharp_peak_candidates(
            self.fscan,
            mask=mask,
            excluded_images=excluded_images,
            burn_in=burn_in,
            history=history,
            min_history=min_history,
            level_z=min_sharpness,
            derivative_z=min_derivative_sharpness,
            lookahead=lookahead,
            min_prominence_z=min_prominence_sharpness,
            refractory=refractory,
            mask_distance=mask_distance,
        )
        for maximum in candidates:
            report(
                "candidate",
                message=(
                    f"Candidate {candidate_count + 1}: image={maximum.imageno}, xy=({maximum.xy[0]:.2f}, {maximum.xy[1]:.2f}), "  # noqa: E501
                    f"max={maximum.value:.4g}, z={maximum.sharpness:.3g}, dz={maximum.derivative_sharpness:.3g}"  # noqa: E501
                ),
                images_read=min(
                    len(self.fscan),
                    maximum.imageno + lookahead + 1,
                ),
            )
            candidate_count += 1
            if max_candidates is not None and candidate_count > max_candidates:
                logger.warning(
                    "Automatic Bragg search stopped after %s rejected candidates.",
                    max_candidates,
                )
                report(
                    "failed",
                    message=(f"Stopped after {max_candidates} rejected candidates."),
                )
                return None
            try:
                peak = autoBraggSearch.refine_peak_3d(
                    self.fscan,
                    maximum,
                    axis_half_width=axis_half_width,
                    roi_size=roi_size,
                    fine_axis_half_width=fine_axis_half_width,
                    fine_roi_size=fine_roi_size,
                    mask=mask,
                    excluded_images=excluded_images,
                    max_workers=self.numberthreads,
                )
            except Exception:
                logger.warning(
                    "Skipping automatic Bragg candidate in image %s after "
                    "3D peak refinement failed.",
                    maximum.imageno,
                    exc_info=True,
                )
                report(
                    "rejected",
                    message=(
                        f"Rejected candidate image {maximum.imageno}: 3D refinement failed."  # noqa: E501
                    ),
                )
                continue
            if not autoBraggSearch.far_from_mask(mask, peak.xy, mask_distance):
                logger.info(
                    "Skipping automatic Bragg candidate in image %s because "
                    "the refined peak is too close to masked pixels.",
                    peak.imageno,
                )
                report(
                    "rejected",
                    message=(
                        f"Rejected candidate image {peak.imageno}: refined peak is too "
                        "close to masked pixels."
                    ),
                )
                continue
            report(
                "refined",
                message=(
                    f"Refined candidate: image={peak.imageno}, axis={peak.axis_value:.5g}, "  # noqa: E501
                    f"xy=({peak.xy[0]:.2f}, {peak.xy[1]:.2f})"
                ),
            )
            refined_peaks.append(peak)
            if (
                not scale_revalidation_done
                and len(unmatched_peaks) >= adaptive_after_candidates
                and len(unmatched_peaks) > last_scale_fit_count
            ):
                broadening_started = not adaptive_matching
                qscale_fit = autoBraggSearch.estimate_qnorm_scale(
                    unmatched_peaks,
                    self.ubcalc.detectorCal,
                    self.getMuOm,
                    self.ubcalc.ubCal,
                    self.ubcalc.crystal,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                    adaptive_qnorm_tolerance,
                    max_q=max_q,
                    hkl_filter=scale_fit_detector_filter,
                )
                last_scale_fit_count = len(unmatched_peaks)
                adaptive_matching = True
                if qscale_fit is not None:
                    qnorm_scale = qscale_fit["qnorm_scale"]
                    if broadening_started:
                        logger.warning(
                            "Automatic Bragg search broadened the Q shell "
                            "from %.4g to %.4g Angstrom^-1 after %s "
                            "unmatched candidates. Tentative reciprocal Q "
                            "scale is %.5g (direct lattice scale %.5g, "
                            "scatter %.3g).",
                            qnorm_tolerance,
                            adaptive_qnorm_tolerance,
                            len(unmatched_peaks),
                            qnorm_scale,
                            qscale_fit["direct_lattice_scale"],
                            qscale_fit["scatter"],
                        )
                        report(
                            "adaptive",
                            message=(
                                "Q shell broadened from {:.4g} to {:.4g} "
                                "Angstrom^-1 after {} unmatched candidates. "
                                "Tentative reciprocal Q scale={:.5g} "
                                "(direct lattice scale={:.5g}, scatter={:.3g}).".format(
                                    qnorm_tolerance,
                                    adaptive_qnorm_tolerance,
                                    len(unmatched_peaks),
                                    qnorm_scale,
                                    qscale_fit["direct_lattice_scale"],
                                    qscale_fit["scatter"],
                                )
                            ),
                        )
                    else:
                        report(
                            "adaptive",
                            message=(
                                "Recomputed adaptive Q scale from {} "
                                "unmatched candidates: reciprocal scale={:.5g} "
                                "(direct lattice scale={:.5g}, scatter={:.3g}).".format(
                                    len(unmatched_peaks),
                                    qnorm_scale,
                                    qscale_fit["direct_lattice_scale"],
                                    qscale_fit["scatter"],
                                )
                            ),
                        )
                    inlier_matches = filter_scale_fit_matches(qscale_fit)
                    if len(inlier_matches) >= adaptive_after_candidates:
                        qnorm_scale = float(
                            np.median(
                                [match["qnorm_scale"] for match in inlier_matches]
                            )
                        )
                        report(
                            "adaptive",
                            message=(
                                "Scale-fit outlier check kept {}/{} "
                                "provisional reflections; rechecking their "
                                "automatic UB validation with inlier Q "
                                "scale {:.5g}.".format(
                                    len(inlier_matches),
                                    len(qscale_fit["matches"]),
                                    qnorm_scale,
                                )
                            ),
                        )
                        for match in inlier_matches:
                            seed, confirmation, _reason = validate_automatic_peak(
                                match["peak"],
                                adaptive_qnorm_tolerance,
                                adaptive_assignment_pixel_tolerance,
                                adaptive_confirmation_pixel_tolerance,
                                adaptive_confirmation_image_tolerance,
                            )
                            if seed is not None:
                                return accept_automatic_seed(seed, confirmation)
                        report(
                            "adaptive",
                            message=(
                                f"Rechecked {len(inlier_matches)} scale-fit provisional "  # noqa: E501
                                "reflections; no second-peak confirmation "
                                "was found."
                            ),
                        )
                        scale_revalidation_done = True
                    else:
                        report(
                            "adaptive",
                            message=(
                                "Scale-fit outlier check kept {}/{} "
                                "provisional reflections; continuing the "
                                "stream until at least {} remain.".format(
                                    len(inlier_matches),
                                    len(qscale_fit["matches"]),
                                    adaptive_after_candidates,
                                )
                            ),
                        )
                else:
                    if broadening_started:
                        logger.warning(
                            "Automatic Bragg search broadened the Q shell "
                            "from %.4g to %.4g Angstrom^-1 after %s "
                            "unmatched candidates. No reliable Q-scale fit "
                            "was available.",
                            qnorm_tolerance,
                            adaptive_qnorm_tolerance,
                            len(unmatched_peaks),
                        )
                        report(
                            "adaptive",
                            message=(
                                f"Q shell broadened from {qnorm_tolerance:.4g} to {adaptive_qnorm_tolerance:.4g} "  # noqa: E501
                                f"Angstrom^-1 after {len(unmatched_peaks)} unmatched candidates. "  # noqa: E501
                                "No reliable Q-scale fit was available."
                            ),
                        )
                    else:
                        report(
                            "adaptive",
                            message=(
                                f"Rechecked adaptive Q-scale fit from {len(unmatched_peaks)} "  # noqa: E501
                                "unmatched candidates; no reliable Q-scale "
                                "fit was available."
                            ),
                        )

            current_qnorm_tolerance = (
                adaptive_qnorm_tolerance if adaptive_matching else qnorm_tolerance
            )
            current_assignment_pixel_tolerance = (
                adaptive_assignment_pixel_tolerance
                if adaptive_matching
                else kwargs.get("assignment_pixel_tolerance", 30.0)
            )
            current_confirmation_pixel_tolerance = (
                adaptive_confirmation_pixel_tolerance
                if adaptive_matching
                else kwargs.get("confirmation_pixel_tolerance", 8.0)
            )
            current_confirmation_image_tolerance = (
                adaptive_confirmation_image_tolerance
                if adaptive_matching
                else kwargs.get("confirmation_image_tolerance", 3)
            )

            seed, confirmation, rejection_reason = validate_automatic_peak(
                peak,
                current_qnorm_tolerance,
                current_assignment_pixel_tolerance,
                current_confirmation_pixel_tolerance,
                current_confirmation_image_tolerance,
            )
            if seed is None:
                logger.warning(
                    "Skipping automatic Bragg candidate in image %s because %s.",
                    peak.imageno,
                    rejection_reason,
                )
                report(
                    "rejected",
                    message=(
                        f"Rejected candidate image {peak.imageno}: {rejection_reason}."
                    ),
                )
                unmatched_peaks.append(peak)
                continue
            return accept_automatic_seed(seed, confirmation)

        logger.warning(
            "Automatic Bragg search reached the end of the scan without a "
            "valid Bragg seed."
        )
        report(
            "failed",
            message="Reached the end of the scan without a valid Bragg seed.",
            images_read=len(self.fscan),
        )
        return None

    def _rankAutomaticSeedsByDetectorPosition(
        self,
        peak,
        q_tolerance=0.05,
        max_q=None,
        pixel_tolerance=30.0,
        max_reflections=40,
        qnorm_scale=1.0,
    ):
        """Rank HKL hypotheses by predicted detector position.

        :param orgui.app.autoBraggSearch.RefinedPeak peak:
            Refined first peak.
        :param float q_tolerance:
            Q-shell half-width around the measured peak norm in Angstrom^-1.
        :param float max_q:
            Optional maximum allowed reflection norm in Angstrom^-1.
        :param float pixel_tolerance:
            Maximum detector-position mismatch in pixels for candidate seeds.
        :param int max_reflections:
            Maximum intensity-sorted HKLs from the Q shell to test.
        :param float qnorm_scale:
            Multiplicative scale applied to predicted reciprocal-vector norms
            before Q-shell matching.
        :returns:
            Candidate seeds sorted by detector mismatch, then intensity rank.
        :rtype: list[orgui.app.autoBraggSearch.UBSeed]
        """
        q_phi = autoBraggSearch.q_phi_from_peak(
            peak,
            self.ubcalc.detectorCal,
            self.getMuOm,
            self.ubcalc.chi,
            self.ubcalc.phi,
            self.ubcalc.ubCal.getK(),
        )
        qnorm = float(np.linalg.norm(q_phi))
        hkls, intensity, norms = autoBraggSearch.allowed_bragg_in_q_region(
            self.ubcalc.crystal,
            qnorm,
            q_tolerance,
            max_q=max_q,
            qnorm_scale=qnorm_scale,
        )
        if len(hkls) > 0:
            expanded = []
            seen = set()
            for hkl, inten, norm in zip(
                hkls[:max_reflections],
                intensity[:max_reflections],
                norms[:max_reflections],
            ):
                same_q_tolerance = max(1e-6, 1e-5 * max(1.0, abs(norm)))
                group_hkls, group_intensity, group_norms = (
                    autoBraggSearch.allowed_bragg_same_q(
                        self.ubcalc.crystal,
                        float(norm),
                        tolerance=same_q_tolerance,
                        max_q=max_q,
                    )
                )
                if len(group_hkls) == 0:
                    group_hkls = np.asarray([hkl], dtype=float)
                    group_intensity = np.asarray([inten], dtype=float)
                    group_norms = np.asarray([norm], dtype=float)
                for group_hkl, group_inten, group_norm in zip(
                    group_hkls, group_intensity, group_norms
                ):
                    key = tuple(np.round(group_hkl, 8))
                    if key in seen:
                        continue
                    seen.add(key)
                    expanded.append((group_hkl, group_inten, group_norm))
            if expanded:
                hkls = np.asarray([item[0] for item in expanded], dtype=float)
                intensity = np.asarray([item[1] for item in expanded], dtype=float)
                norms = np.asarray([item[2] for item in expanded], dtype=float)
                order = np.argsort(intensity)[::-1]
                hkls = hkls[order]
                intensity = intensity[order]
                norms = norms[order]
        previous_u = np.asarray(self.ubcalc.ubCal.getU()).copy()
        seeds = []
        try:
            for hkl, _intensity, norm in zip(
                hkls,
                intensity,
                norms,
            ):
                try:
                    refldict = self.searchPixelCoordHKL(hkl)
                except Exception:
                    logger.debug(
                        "Cannot predict first-peak HKL hypothesis %s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                solution_mismatches = []
                for solution in (1, 2):
                    if not refldict.get(f"selectable_{solution}", False):
                        continue
                    if f"imageno_{solution}" not in refldict:
                        continue
                    predicted_xy = np.asarray(refldict[f"xy_{solution}"], dtype=float)
                    predicted_imageno = refldict[f"imageno_{solution}"]
                    pixel_error = np.linalg.norm(peak.xy - predicted_xy)
                    image_error = abs(peak.imageno - int(predicted_imageno))
                    solution_mismatches.append((pixel_error + image_error, pixel_error))
                if not solution_mismatches:
                    continue
                score, pixel_error = min(solution_mismatches)
                if pixel_error > pixel_tolerance:
                    continue
                mu, omega = self.getMuOm(peak.imageno)
                gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                    np.array([peak.xy[1]]),
                    np.array([peak.xy[0]]),
                    mu,
                )
                pos = np.array(
                    [
                        mu,
                        float(delta[0]),
                        float(gamma[0]),
                        omega,
                        self.ubcalc.chi,
                        self.ubcalc.phi,
                    ]
                )
                try:
                    U = autoBraggSearch.seed_u_from_single_reflection(
                        self.ubcalc.ubCal, pos, hkl
                    )
                except Exception:
                    logger.debug(
                        "Cannot calculate one-reflection UB seed for hkl=%s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                seeds.append(
                    autoBraggSearch.UBSeed(
                        U=U,
                        hkl=np.asarray(hkl, dtype=float),
                        peak=peak,
                        score=float(score),
                        norm_mismatch=abs(float(norm) * qnorm_scale - qnorm),
                        position_mismatch=float(pixel_error),
                    )
                )
        finally:
            self.ubcalc.ubCal.setU(previous_u)
        return sorted(
            seeds,
            key=lambda seed: (
                seed.position_mismatch,
                seed.norm_mismatch,
                seed.score,
            ),
        )

    def _automaticReflectionIntensity(self, hkl):
        """Return calculated ``abs(F)**2`` for one Bragg reflection."""
        hkl = np.asarray(hkl, dtype=float)
        h = np.asarray([hkl[0]], dtype=float)
        k = np.asarray([hkl[1]], dtype=float)
        l = np.asarray([hkl[2]], dtype=float)  # noqa: E741
        if hasattr(self.ubcalc.crystal, "uc_bulk"):
            structure_factor = self.ubcalc.crystal.uc_bulk.F_uc(h, k, l)[0]
        else:
            structure_factor = self.ubcalc.crystal.F_uc(h, k, l)[0]
        return float(np.abs(structure_factor) ** 2)

    def _addAutomaticObservedReflections(
        self,
        peaks,
        seed,
        confirmation,
        q_tolerance,
        score_tolerance_deg=2.0,
        intensity_ratio_check=True,
        intensity_ratio_tolerance=0.5,
        status_callback=None,
    ):
        """Add already-observed peaks that match the confirmed UB matrix.

        :param list[orgui.app.autoBraggSearch.RefinedPeak] peaks:
            Refined peaks already observed during the automatic search.
        :param orgui.app.autoBraggSearch.UBSeed seed:
            Accepted seed reflection.
        :param dict confirmation:
            Accepted confirmation reflection dictionary.
        :param float q_tolerance:
            Maximum Q-norm mismatch in Angstrom^-1.
        :param float score_tolerance_deg:
            Maximum angular UB mismatch in deg.
        :param bool intensity_ratio_check:
            If true, require each added reflection intensity to match the mean
            observed/calculated intensity scale within tolerance.
        :param float intensity_ratio_tolerance:
            Maximum fractional deviation from the mean intensity scale.
        :param callable status_callback:
            Optional automatic-search status callback for GUI/CLI feedback.
        :returns:
            Number of additional reflections added.
        :rtype: int
        """
        if not peaks:
            return 0
        accepted = []
        for accepted_peak in (seed.peak, confirmation.get("peak")):
            if accepted_peak is not None:
                accepted.append(accepted_peak)
        intensity_scale = []
        if intensity_ratio_check:
            intensity_refs = [
                (seed.peak, seed.hkl),
                (confirmation.get("peak"), confirmation.get("hkl")),
            ]
            for ref_peak, ref_hkl in intensity_refs:
                if ref_peak is None or ref_hkl is None:
                    continue
                ref_intensity = autoBraggSearch.estimate_rocking_intensity(ref_peak)
                ref_f2 = self._automaticReflectionIntensity(ref_hkl)
                if (
                    ref_intensity is not None
                    and ref_intensity["intensity"] > 0.0
                    and ref_f2 > 0.0
                ):
                    intensity_scale.append(ref_intensity["intensity"] / ref_f2)

        def is_same_peak(peak_a, peak_b):
            return (
                peak_a.imageno == peak_b.imageno
                and np.linalg.norm(peak_a.xy - peak_b.xy) < 2.0
            )

        added = 0
        tested_peaks = []
        for peak in peaks:
            if any(is_same_peak(peak, existing) for existing in accepted):
                continue
            if any(is_same_peak(peak, existing) for existing in tested_peaks):
                continue
            tested_peaks.append(peak)
            mu, omega = self.getMuOm(peak.imageno)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([peak.xy[1]]),
                np.array([peak.xy[0]]),
                mu,
            )
            angles = np.array(
                [
                    mu,
                    float(delta[0]),
                    float(gamma[0]),
                    omega,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                ]
            )
            try:
                h, k, l = self.ubcalc.angles.anglesToHkl(*angles)  # noqa: E741
            except Exception:
                logger.debug(
                    "Cannot convert observed automatic peak at image %s, xy=%s to HKL.",
                    peak.imageno,
                    peak.xy,
                    exc_info=True,
                )
                continue
            measured_hkl = np.asarray([h, k, l], dtype=float).reshape(3)
            rounded_hkl = np.rint(measured_hkl).astype(float)
            if not np.all(np.isfinite(measured_hkl)) or np.allclose(rounded_hkl, 0.0):
                continue

            if status_callback is not None:
                status_callback(
                    "observed_hkl",
                    message=(
                        f"Observed peak image={peak.imageno}, xy=({peak.xy[0]:.2f}, {peak.xy[1]:.2f}): "  # noqa: E501
                        f"HKL from angles=({measured_hkl[0]:.3f}, {measured_hkl[1]:.3f}, {measured_hkl[2]:.3f}), rounded={rounded_hkl}."  # noqa: E501
                    ),
                )
            logger.info(
                "Automatic observed peak image=%s, xy=%s gives HKL=%s; rounded to %s.",
                peak.imageno,
                peak.xy,
                measured_hkl,
                rounded_hkl,
            )
            try:
                mismatch = self.ubcalc.getReflectionMismatch(
                    np.atleast_2d(rounded_hkl),
                    np.atleast_2d(angles),
                )
            except Exception:
                logger.debug(
                    "Cannot score observed automatic peak hkl=%s",
                    rounded_hkl,
                    exc_info=True,
                )
                continue
            angular_mismatch_deg = float(
                np.rad2deg(np.asarray(mismatch["angle_mismatch"])[0])
            )
            q_mismatch = float(np.asarray(mismatch["norm_mismatch"])[0])
            message = (
                f"Observed peak image={peak.imageno} rounded hkl={rounded_hkl}: angular "  # noqa: E501
                f"mismatch={angular_mismatch_deg:.4g} deg, Q mismatch={q_mismatch:.4g} Angstrom^-1 "  # noqa: E501
                f"(limits {score_tolerance_deg:.4g} deg, {q_tolerance:.4g} Angstrom^-1)."  # noqa: E501
            )
            logger.info(message)
            if status_callback is not None:
                status_callback("observed_hkl_test", message=message)
            if angular_mismatch_deg > score_tolerance_deg or q_mismatch > q_tolerance:
                continue
            intensity_message = ""
            peak_intensity = None
            if intensity_ratio_check:
                peak_intensity = autoBraggSearch.estimate_rocking_intensity(peak)
                peak_f2 = self._automaticReflectionIntensity(rounded_hkl)
                if (
                    peak_intensity is None
                    or peak_intensity["intensity"] <= 0.0
                    or peak_f2 <= 0.0
                    or not intensity_scale
                ):
                    logger.info(
                        "Rejected observed automatic Bragg reflection hkl=%s "
                        "at image %s because intensity scale could not be "
                        "tested.",
                        rounded_hkl,
                        peak.imageno,
                    )
                    continue
                observed_scale = peak_intensity["intensity"] / peak_f2
                mean_scale = float(np.mean(intensity_scale))
                scale_error = abs(observed_scale - mean_scale) / mean_scale
                intensity_message = f", intensity scale={observed_scale:.4g}, mean scale={mean_scale:.4g}, relerr={scale_error:.3g}"  # noqa: E501
                message = (
                    f"Observed peak image={peak.imageno} rounded hkl={rounded_hkl} intensity check: "  # noqa: E501
                    f"scale={observed_scale:.4g}, mean scale={mean_scale:.4g}, relerr={scale_error:.3g} "  # noqa: E501
                    f"(limit {intensity_ratio_tolerance:.3g})."
                )
                logger.info(message)
                if status_callback is not None:
                    status_callback("observed_intensity_test", message=message)
                if scale_error > intensity_ratio_tolerance:
                    continue
            self.reflectionSel.addReflection(
                {"x": peak.xy[0], "y": peak.xy[1]},
                peak.imageno,
                rounded_hkl,
            )
            if (
                intensity_ratio_check
                and peak_intensity is not None
                and peak_intensity["intensity"] > 0.0
            ):
                intensity_scale.append(observed_scale)
            added += 1
            logger.info(
                "Added observed automatic Bragg reflection hkl=%s at image "
                "%s, xy=%s, angular mismatch=%.4g deg, Q mismatch=%.4g "
                "Angstrom^-1%s.",
                rounded_hkl,
                peak.imageno,
                peak.xy,
                angular_mismatch_deg,
                q_mismatch,
                intensity_message,
            )
        return added

    def _automaticPeakAngles(self, peak):
        """Return Vlieg angles in rad for one refined automatic peak."""
        mu, omega = self.getMuOm(peak.imageno)
        gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
            np.array([peak.xy[1]]),
            np.array([peak.xy[0]]),
            mu,
        )
        return np.array(
            [
                mu,
                float(delta[0]),
                float(gamma[0]),
                omega,
                self.ubcalc.chi,
                self.ubcalc.phi,
            ]
        )

    def _recalculateAutomaticUFromReflections(self):
        """Recalculate U from current reference reflections without dialogs."""
        hkls, angles = self.getReflections()
        if len(hkls) < 1:
            raise ValueError("At least one reflection is required")
        if len(hkls) == 1:
            self.ubcalc.ubCal.zmodeUSingleRefl(angles[0], hkls[0])
        elif len(hkls) == 2:
            self.ubcalc.ubCal.setPrimaryReflection(angles[0], hkls[0])
            self.ubcalc.ubCal.setSecondayReflection(angles[1], hkls[1])
            self.ubcalc.ubCal.calculateU()
        else:
            self.ubcalc.ubCal.calculateUFromReflections(hkls, angles)
        self.ubcalc.uedit.setU(self.ubcalc.ubCal.getU())
        self.ubcalc.updateReflectionMismatch()
        self.ubcalc.sigReplotRequest.emit(False)

    def _automaticReflectionIntensityScale(
        self,
        refl,
        mask=None,
        axis_half_width=1.0,
        fine_axis_half_width=0.4,
        roi_size=(80, 80),
        fine_roi_size=(40, 40),
    ):
        """Return observed/calculated intensity scale for one reflection."""
        candidate = autoBraggSearch.ImageMaximum(
            int(refl.imageno),
            np.asarray(refl.xy, dtype=float),
            np.nan,
        )
        peak = autoBraggSearch.refine_peak_3d(
            self.fscan,
            candidate,
            axis_half_width=axis_half_width,
            roi_size=roi_size,
            fine_axis_half_width=fine_axis_half_width,
            fine_roi_size=fine_roi_size,
            mask=mask,
            excluded_images=self.excludedImagesDialog.getData(),
            max_workers=self.numberthreads,
        )
        intensity = autoBraggSearch.estimate_rocking_intensity(peak)
        f2 = self._automaticReflectionIntensity(refl.hkl)
        if intensity is None or intensity["intensity"] <= 0.0 or f2 <= 0.0:
            return None, peak
        return intensity["intensity"] / f2, peak

    def autoAddCalculatedBraggReflections(self, count, **kwargs):
        """Validate and add calculated Bragg reflections from the current UB.

        :param int count:
            Number of additional reflections to add.
        :param kwargs:
            Automatic Bragg options from
            :class:`orgui.app.QReflectionSelector.AutoBraggOptionsDialog`.
        :returns:
            Number of accepted reflections.
        :rtype: int

        .. note::
           GUI-only. This method mutates the reference-reflection table,
           performs local image peak searches, and emits GUI update signals.
        """
        if self.fscan is None:
            raise ValueError("No scan loaded")
        count = int(count)
        if count <= 0:
            return 0
        status_callback = kwargs.get("status_callback", None)

        def report(event, **fields):
            if status_callback is not None:
                status_callback(event, **fields)

        q_tolerance = kwargs.get(
            "additional_q_tolerance",
            kwargs.get(
                "adaptive_qnorm_tolerance",
                kwargs.get("qnorm_tolerance", 0.05),
            ),
        )
        score_tolerance_deg = kwargs.get("additional_score_tolerance_deg", 2.0)
        intensity_ratio_check = kwargs.get("intensity_ratio_check", True)
        intensity_ratio_tolerance = kwargs.get("intensity_ratio_tolerance", 0.5)
        axis_half_width = kwargs.get("axis_half_width", 1.0)
        fine_axis_half_width = kwargs.get("fine_axis_half_width", 0.4)
        roi_size = kwargs.get("roi_size", (80, 80))
        fine_roi_size = kwargs.get("fine_roi_size", (40, 40))
        pixel_tolerance = kwargs.get(
            "additional_pixel_tolerance",
            kwargs.get("confirmation_pixel_tolerance", 8.0),
        )
        image_tolerance = kwargs.get(
            "additional_image_tolerance",
            kwargs.get("confirmation_image_tolerance", 3),
        )

        mask = None
        mask_widget = self.centralPlot.getMaskToolsDockWidget()
        if mask_widget.getSelectionMask() is not None:
            mask = mask_widget.getSelectionMask() > 0.0
        else:
            logger.warning(
                "Adding calculated Bragg reflections without a detector mask."
            )
            report(
                "warning",
                message=(
                    "No detector mask is set; additional calculated "
                    "reflections will be refined without mask rejection."
                ),
            )

        intensity_scale = []
        if intensity_ratio_check:
            for refl in list(self.reflectionSel.reflections):
                try:
                    scale, _peak = self._automaticReflectionIntensityScale(
                        refl,
                        mask=mask,
                        axis_half_width=axis_half_width,
                        fine_axis_half_width=fine_axis_half_width,
                        roi_size=roi_size,
                        fine_roi_size=fine_roi_size,
                    )
                except Exception:
                    logger.debug(
                        "Cannot estimate intensity scale for hkl=%s",
                        refl.hkl,
                        exc_info=True,
                    )
                    continue
                if scale is not None and np.isfinite(scale):
                    intensity_scale.append(float(scale))
            report(
                "additional_intensity_scale",
                message=(
                    f"Built intensity scale list from {len(intensity_scale)} current reference "  # noqa: E501
                    "reflection(s)."
                ),
            )

        added = 0
        rejected = 0
        attempted_hkls = set()
        while added < count:
            candidates = self.reflectionSel.getBraggCandidates(recalculate=True)
            if not candidates:
                break
            accepted_this_round = False
            for candidate_refl in candidates:
                hkl = np.asarray(candidate_refl.hkl, dtype=float)
                key = tuple(np.round(hkl, 8))
                if key in attempted_hkls:
                    continue
                attempted_hkls.add(key)
                report(
                    "additional_candidate",
                    message=(
                        f"Testing calculated Bragg hkl={hkl} at image={candidate_refl.imageno}, "  # noqa: E501
                        f"xy=({candidate_refl.xy[0]:.2f}, {candidate_refl.xy[1]:.2f})."
                    ),
                )
                maximum = autoBraggSearch.ImageMaximum(
                    int(candidate_refl.imageno),
                    np.asarray(candidate_refl.xy, dtype=float),
                    np.nan,
                )
                try:
                    peak = autoBraggSearch.refine_peak_3d(
                        self.fscan,
                        maximum,
                        axis_half_width=axis_half_width,
                        roi_size=roi_size,
                        fine_axis_half_width=fine_axis_half_width,
                        fine_roi_size=fine_roi_size,
                        mask=mask,
                        excluded_images=self.excludedImagesDialog.getData(),
                        max_workers=self.numberthreads,
                    )
                except Exception:
                    rejected += 1
                    logger.debug(
                        "Cannot refine calculated Bragg hkl=%s",
                        hkl,
                        exc_info=True,
                    )
                    report(
                        "rejected",
                        message=(
                            f"Rejected calculated hkl={hkl}: 3D peak search failed."
                        ),
                    )
                    continue
                pixel_error = float(
                    np.linalg.norm(peak.xy - np.asarray(candidate_refl.xy, dtype=float))
                )
                image_error = abs(int(peak.imageno) - int(candidate_refl.imageno))
                if pixel_error > pixel_tolerance or image_error > image_tolerance:
                    rejected += 1
                    report(
                        "rejected",
                        message=(
                            f"Rejected calculated hkl={hkl}: refined peak moved "
                            f"{pixel_error:.4g} px and {image_error} image(s) (limits {pixel_tolerance:.4g} px, {image_tolerance})."  # noqa: E501
                        ),
                    )
                    continue

                angles = self._automaticPeakAngles(peak)
                try:
                    mismatch = self.ubcalc.getReflectionMismatch(
                        np.atleast_2d(hkl),
                        np.atleast_2d(angles),
                    )
                except Exception:
                    rejected += 1
                    logger.debug(
                        "Cannot score calculated Bragg hkl=%s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                angular_mismatch_deg = float(
                    np.rad2deg(np.asarray(mismatch["angle_mismatch"])[0])
                )
                q_mismatch = float(np.asarray(mismatch["norm_mismatch"])[0])
                report(
                    "additional_mismatch",
                    message=(
                        f"Calculated hkl={hkl}: angular mismatch={angular_mismatch_deg:.4g} deg, "  # noqa: E501
                        f"Q mismatch={q_mismatch:.4g} Angstrom^-1 (limits {score_tolerance_deg:.4g} deg, "  # noqa: E501
                        f"{q_tolerance:.4g} Angstrom^-1)."
                    ),
                )
                if (
                    angular_mismatch_deg > score_tolerance_deg
                    or q_mismatch > q_tolerance
                ):
                    rejected += 1
                    continue

                scale_message = ""
                observed_scale = None
                if intensity_ratio_check:
                    intensity = autoBraggSearch.estimate_rocking_intensity(peak)
                    f2 = self._automaticReflectionIntensity(hkl)
                    if (
                        intensity is None
                        or intensity["intensity"] <= 0.0
                        or f2 <= 0.0
                        or not intensity_scale
                    ):
                        rejected += 1
                        report(
                            "rejected",
                            message=(
                                f"Rejected calculated hkl={hkl}: intensity scale "
                                "could not be tested."
                            ),
                        )
                        continue
                    observed_scale = intensity["intensity"] / f2
                    mean_scale = float(np.mean(intensity_scale))
                    scale_error = abs(observed_scale - mean_scale) / mean_scale
                    report(
                        "additional_intensity_test",
                        message=(
                            f"Calculated hkl={hkl} intensity check: scale={observed_scale:.4g}, "  # noqa: E501
                            f"mean scale={mean_scale:.4g}, relerr={scale_error:.3g} (limit {intensity_ratio_tolerance:.3g})."  # noqa: E501
                        ),
                    )
                    if scale_error > intensity_ratio_tolerance:
                        rejected += 1
                        continue
                    scale_message = f", intensity scale={observed_scale:.4g}, mean scale={mean_scale:.4g}, relerr={scale_error:.3g}"  # noqa: E501

                self.reflectionSel.addReflection(
                    {"x": peak.xy[0], "y": peak.xy[1]},
                    peak.imageno,
                    hkl,
                )
                if intensity_ratio_check and observed_scale is not None:
                    intensity_scale.append(float(observed_scale))
                self._recalculateAutomaticUFromReflections()
                added += 1
                accepted_this_round = True
                logger.info(
                    "Added calculated Bragg reflection hkl=%s at image %s, "
                    "xy=%s, angular mismatch=%.4g deg, Q mismatch=%.4g "
                    "Angstrom^-1%s.",
                    hkl,
                    peak.imageno,
                    peak.xy,
                    angular_mismatch_deg,
                    q_mismatch,
                    scale_message,
                )
                report(
                    "accepted",
                    message=(
                        f"Added calculated hkl={hkl} at image={peak.imageno}, xy=({peak.xy[0]:.2f}, {peak.xy[1]:.2f})"  # noqa: E501
                        f"{scale_message}. Recalculated U."
                    ),
                )
                break
            if not accepted_this_round:
                break

        report(
            "additional_done",
            message=(
                f"Additional calculated Bragg pass complete: added={added}, "
                f"rejected={rejected}."
            ),
        )
        return added

    def _confirmAutomaticBraggSeed(
        self,
        seed,
        mask=None,
        excluded_images=(),
        max_q=None,
        axis_half_width=1.0,
        fine_axis_half_width=0.4,
        roi_size=(80, 80),
        fine_roi_size=(40, 40),
        pixel_tolerance=8.0,
        image_tolerance=3,
        max_reflections=12,
        status_callback=None,
        intensity_ratio_tolerance=0.5,
        intensity_ratio_check=True,
        prominence_threshold=6.0,
    ):
        """Confirm an automatic seed by finding another predicted Bragg peak.

        :param orgui.app.autoBraggSearch.UBSeed seed:
            Candidate seed from the first detected Bragg peak.
        :param numpy.ndarray mask:
            Optional boolean detector mask where ``True`` pixels are invalid.
        :param float max_q:
            Maximum reciprocal-vector norm in Angstrom^-1.
        :returns:
            Confirmation dictionary with ``hkl``, predicted reflection, and
            refined peak, or ``None``.
        :rtype: dict or None

        .. note::
           CLI-capable. Temporarily applies ``seed.U`` only while calculating
           predicted positions and always restores the previous U.
        """
        if max_q is None:
            max_q = getattr(self.ubcalc.detectorCal, "Qmax", None)
        if max_q is None:
            return None
        hkls, _intensity, _norms = autoBraggSearch.allowed_bragg_with_intensity(
            self.ubcalc.crystal, max_q
        )
        seed_intensity = autoBraggSearch.estimate_rocking_intensity(seed.peak)
        seed_f2 = self._automaticReflectionIntensity(seed.hkl)
        previous_u = np.asarray(self.ubcalc.ubCal.getU()).copy()
        try:
            self.ubcalc.ubCal.setU(seed.U)
            tested = 0
            for hkl in hkls:
                if np.allclose(hkl, seed.hkl):
                    continue
                try:
                    refldict = self.searchPixelCoordHKL(hkl)
                except Exception:
                    logger.debug(
                        "Cannot predict confirmation reflection %s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                for solution in (1, 2):
                    if not refldict.get(f"selectable_{solution}", False):
                        continue
                    if f"imageno_{solution}" not in refldict:
                        continue
                    predicted_xy = np.asarray(refldict[f"xy_{solution}"], dtype=float)
                    predicted_imageno = refldict[f"imageno_{solution}"]
                    candidate = autoBraggSearch.ImageMaximum(
                        int(predicted_imageno),
                        predicted_xy,
                        np.nan,
                    )
                    try:
                        peak = autoBraggSearch.refine_peak_3d(
                            self.fscan,
                            candidate,
                            axis_half_width=axis_half_width,
                            roi_size=roi_size,
                            fine_axis_half_width=fine_axis_half_width,
                            fine_roi_size=fine_roi_size,
                            mask=mask,
                            excluded_images=excluded_images,
                            max_workers=self.numberthreads,
                        )
                    except Exception:
                        logger.debug(
                            "Cannot refine confirmation reflection %s",
                            hkl,
                            exc_info=True,
                        )
                        continue
                    pixel_error = np.linalg.norm(peak.xy - predicted_xy)
                    image_error = abs(peak.imageno - int(predicted_imageno))
                    if (
                        pixel_error <= pixel_tolerance
                        and image_error <= image_tolerance
                    ):
                        confirm_intensity = autoBraggSearch.estimate_rocking_intensity(
                            peak
                        )
                        confirm_f2 = self._automaticReflectionIntensity(hkl)
                        intensity_ok = True
                        observed_ratio = np.nan
                        predicted_ratio = np.nan
                        ratio_error = np.nan
                        seed_prominence = np.nan
                        confirm_prominence = np.nan
                        if seed_intensity is not None:
                            seed_prominence = seed_intensity["prominence_z"]
                        if confirm_intensity is not None:
                            confirm_prominence = confirm_intensity["prominence_z"]
                        if intensity_ratio_check:
                            if (
                                seed_intensity is not None
                                and confirm_intensity is not None
                                and seed_intensity["intensity"] > 0.0
                                and confirm_intensity["intensity"] > 0.0
                                and seed_f2 > 0.0
                                and confirm_f2 > 0.0
                            ):
                                observed_ratio = (
                                    confirm_intensity["intensity"]
                                    / seed_intensity["intensity"]
                                )
                                predicted_ratio = confirm_f2 / seed_f2
                                ratio_error = (
                                    abs(observed_ratio - predicted_ratio)
                                    / predicted_ratio
                                )
                                intensity_ok = ratio_error <= intensity_ratio_tolerance
                            else:
                                intensity_ok = False
                        else:
                            intensity_ok = (
                                np.isfinite(seed_prominence)
                                and np.isfinite(confirm_prominence)
                                and seed_prominence >= prominence_threshold
                                and confirm_prominence >= prominence_threshold
                            )
                        if not intensity_ok:
                            logger.debug(
                                "Reject confirmation hkl=%s for seed hkl=%s "
                                "because intensity check failed. observed=%s "
                                "predicted=%s relerr=%s seed_prominence=%s "
                                "confirm_prominence=%s",
                                hkl,
                                seed.hkl,
                                observed_ratio,
                                predicted_ratio,
                                ratio_error,
                                seed_prominence,
                                confirm_prominence,
                            )
                            continue
                        return {
                            "hkl": np.asarray(hkl, dtype=float),
                            "solution": solution,
                            "predicted_xy": predicted_xy,
                            "predicted_imageno": int(predicted_imageno),
                            "pixel_error": pixel_error,
                            "image_error": image_error,
                            "peak": peak,
                            "seed_intensity": seed_intensity,
                            "confirmation_intensity": confirm_intensity,
                            "observed_intensity_ratio": observed_ratio,
                            "predicted_intensity_ratio": predicted_ratio,
                            "intensity_ratio_error": ratio_error,
                            "seed_prominence_z": seed_prominence,
                            "confirmation_prominence_z": confirm_prominence,
                        }
                tested += 1
                if status_callback is not None:
                    status_callback(
                        "confirmation_search",
                        message=(
                            f"Tried {tested} confirmation reflections for seed "
                            f"hkl={seed.hkl}."
                        ),
                    )
                if tested >= max_reflections:
                    break
        finally:
            self.ubcalc.ubCal.setU(previous_u)
        return None

    def saveBraggRefl(self):
        """Save calculated Bragg reflection coordinates selected in the GUI.

        .. note::
           GUI-only. This path opens confirmation and file-save dialogs.
        """
        try:
            hkls, yx, angles = self.reflectionSel.getBraggReflections()
        except ValueError:
            if self.fscan is not None:
                try:
                    xtal = self.ubcalc.crystal
                    ommin = np.deg2rad(np.amin(self.fscan.omega))
                    ommax = np.deg2rad(np.amax(self.fscan.omega))
                    dc = self.ubcalc.detectorCal
                    mu = self.ubcalc.mu
                    ub = self.ubcalc.ubCal
                    chi = self.ubcalc.chi
                    phi = self.ubcalc.phi
                    xtal.setEnergy(ub.getEnergy() * 1e3)
                    hkls, yx, angles = rn.thscanBragg(
                        xtal, ub, mu, dc, (ommin, ommax), chi=chi, phi=phi
                    )

                    # self.reflectionSel.setBraggReflections(hkls, yx, angles)
                except Exception:
                    logger.exception(
                        "Cannot calculate Bragg reflections",
                        extra={
                            "title": "Cannot calculate Bragg reflections",
                            "description": "Cannot calculate Bragg reflections",
                            "show_dialog": True,
                            "dialog_level": logging.WARNING,
                            "parent": self,
                        },
                    )
                    # qutils.warning_detailed_message(self, "Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections", traceback.format_exc())  # noqa: E501
                    # qt.QMessageBox.critical(self,"Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections:\n%s" % traceback.format_exc())  # noqa: E501
                    return
            else:
                logger.exception(
                    "Cannot calculate Bragg reflections: No scan loaded.",
                    extra={
                        "title": "Cannot calculate Bragg reflections",
                        "description": "Cannot calculate Bragg reflections: No scan loaded.",  # noqa: E501
                        "show_dialog": True,
                        "dialog_level": logging.WARNING,
                        "parent": self,
                    },
                )
                # qt.QMessageBox.critical(self,"Cannot calculate Bragg reflections", "Cannot calculate Bragg reflections:\nNo scan loaded.")  # noqa: E501
                return

        hkm = np.concatenate((hkls, yx[:, ::-1], np.rad2deg(angles)), axis=1)

        sio = StringIO()
        np.savetxt(
            sio,
            hkm,
            fmt="%.3f",
            delimiter="\t",
            header="H K L x y alpha delta gamma omega chi phi",
        )

        # Question dialog for saving the possible CTR locations
        msgbox = qt.QMessageBox(
            qt.QMessageBox.Question,
            "Saving Bragg reflection ...",
            "Found possible Bragg reflections. Do you want to save the following positions?",  # noqa: E501
            qt.QMessageBox.Yes | qt.QMessageBox.No,
            self,
        )

        msgbox.setDetailedText(sio.getvalue())

        clickedbutton = msgbox.exec()
        # Question dialog for saving the possible CTR locations
        # clickedbutton=qt.QMessageBox.question(self, 'Saving CTR locations...', 'Do you want to save the following positions: \n' + hkstring +"?");  # noqa: E501

        if clickedbutton == qt.QMessageBox.Yes:
            # File saving
            fileTypeDict = {
                "dat Files (*.dat)": ".dat",
                "txt Files (*.txt)": ".txt",
                "All files (*)": "",
            }
            fileTypeFilter = ""
            for f in fileTypeDict:
                fileTypeFilter += f + ";;"

            filename, filetype = qt.QFileDialog.getSaveFileName(
                self, "Save reflections", self.filedialogdir, fileTypeFilter[:-2]
            )
            if filename == "":
                return

            self.filedialogdir = os.path.splitext(filename)[0]
            filename += fileTypeDict[filetype]
            np.savetxt(
                filename,
                hkm,
                fmt="%.3f",
                header="H K L x y alpha delta gamma omega chi phi",
            )

    def calculateAvailableCTR(self):
        """Calculate CTR coordinates available to the current ``th`` scan.

        :returns:
            Array columns containing H, K, and detector-side information.
        :rtype: numpy.ndarray

        .. note::
           CLI-safe when a compatible scan and UB state are loaded.
        """
        if self.fscan is None:
            raise Exception("No scan selected!")
        if self.fscan.axisname != "th":
            raise NotImplementedError(
                f"Calculation of available CTRs is not implemented for {self.fscan.axisname} - scans"  # noqa: E501
            )
        xtal = self.ubcalc.crystal
        ommin = np.deg2rad(np.amin(self.fscan.omega))
        ommax = np.deg2rad(np.amax(self.fscan.omega))
        dc = self.ubcalc.detectorCal
        mu = self.ubcalc.mu
        chi = self.ubcalc.chi
        phi = self.ubcalc.phi
        ub = self.ubcalc.ubCal
        xtal.setEnergy(ub.getEnergy() * 1e3)
        hk, xmirror = rn.thscanCTRs(xtal, ub, mu, dc, (ommin, ommax), chi=chi, phi=phi)
        xmirror = np.array(xmirror).astype(np.float64)
        # making the hk list of arrays into a reasonable string
        hkm = np.concatenate(
            (np.array(hk), xmirror.reshape((1, xmirror.size)).T), axis=1
        )
        return hkm

    def _onCalcAvailableCTR(self):
        """GUI-only: calculate CTRs and optionally save them via dialogs."""
        try:
            hkm = self.calculateAvailableCTR()
        except Exception:
            qutils.warning_detailed_message(
                self,
                "Cannot calculate CTR locatons",
                "Cannot calculate CTR locatons",
                traceback.format_exc(),
            )
            # qt.QMessageBox.critical(self,"Cannot calculate CTR locatons", "Cannot calculate CTR locatons:\n%s" % traceback.format_exc())  # noqa: E501
            return
        sio = StringIO()
        np.savetxt(sio, hkm, fmt="%.3f", delimiter="\t", header="H K detectorRight")

        # Question dialog for saving the possible CTR locations
        msgbox = qt.QMessageBox(
            qt.QMessageBox.Question,
            "Saving CTR locations...",
            "Found CTRs. Do you want to save the following positions?",
            qt.QMessageBox.Yes | qt.QMessageBox.No,
            self,
        )
        msgbox.setDetailedText(sio.getvalue())

        clickedbutton = msgbox.exec()
        # Question dialog for saving the possible CTR locations
        # clickedbutton=qt.QMessageBox.question(self, 'Saving CTR locations...', 'Do you want to save the following positions: \n' + hkstring +"?");  # noqa: E501

        if clickedbutton == qt.QMessageBox.Yes:
            # File saving
            fileTypeDict = {
                "dat Files (*.dat)": ".dat",
                "txt Files (*.txt)": ".txt",
                "All files (*)": "",
            }
            fileTypeFilter = ""
            for f in fileTypeDict:
                fileTypeFilter += f + ";;"

            filename, filetype = qt.QFileDialog.getSaveFileName(
                self, "Save reflections", self.filedialogdir, fileTypeFilter[:-2]
            )
            if filename == "":
                return

            self.filedialogdir = os.path.splitext(filename)[0]
            filename += fileTypeDict[filetype]
            np.savetxt(filename, hkm, fmt="%.3f", header="H K mirror")

    def getReflections(self):
        """Return selected reference reflections for UB calculation.

        The reference reflections are taken from the reflection selector. For
        each selected reflection, the detector pixel coordinate is converted to
        surface detector angles ``delta`` and ``gamma`` in rad, and the image
        number is converted to ``omega`` in rad. The returned angle rows use
        the Vlieg six-circle order ``[alpha, delta, gamma, omega, chi, phi]``.
        In this GUI workflow, ``alpha`` is the current fixed incidence angle
        stored as ``self.ubcalc.mu``.

        :returns:
            Tuple ``(hkls, angles)``. ``hkls`` has columns ``[h, k, l]`` in
            r.l.u. ``angles`` has columns
            ``[alpha, delta, gamma, omega, chi, phi]`` in rad.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]

        .. note::
           CLI-safe when reference-selection state is already populated.
        """
        hkls = []
        angles = []
        for refl in self.reflectionSel.reflections:
            if not self.isValidImageNo(refl.imageno):
                logger.warning(
                    "Skipping stale reflection %s with image number %s "
                    "outside the active scan range.",
                    refl.hkl,
                    refl.imageno,
                )
                continue
            # print(refl.xy)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([refl.xy[1]]), np.array([refl.xy[0]]), self.ubcalc.mu
            )
            delta = float(delta[0])
            gamma = float(gamma[0])
            try:
                pos = np.array(
                    [
                        self.ubcalc.mu,
                        delta,
                        gamma,
                        self.imageNoToOmega(refl.imageno),
                        self.ubcalc.chi,
                        self.ubcalc.phi,
                    ]
                )
            except Exception:
                # from IPython import embed; embed()
                raise
            # print(pos)
            hkls.append(refl.hkl)
            angles.append(pos)
        return np.array(hkls), np.array(angles)

    def getReflectionFitData(self):
        """Return reference reflections including detector pixel positions.

        :returns:
            ``(hkls, angles, xy)`` with HKL in r.l.u., Vlieg angles in rad,
            and detector coordinates ``(x, y)`` in pixels.
        :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]
        """
        hkls = []
        angles = []
        xy = []
        for refl in self.reflectionSel.reflections:
            if not self.isValidImageNo(refl.imageno):
                continue
            mu, omega = self.getMuOm(refl.imageno)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([refl.xy[1]]),
                np.array([refl.xy[0]]),
                mu,
            )
            hkls.append(refl.hkl)
            angles.append(
                [
                    mu,
                    float(delta[0]),
                    float(gamma[0]),
                    omega,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                ]
            )
            xy.append(refl.xy)
        return np.asarray(hkls), np.asarray(angles), np.asarray(xy)

    def _onPlotMachineParams(self, enable=None):
        """GUI/CLI hint: toggle detector center and azimuth plot markers."""
        # [cp,azimxy,polax] = paramslist
        if enable is None:
            enable = self.showMachineParamsAct.isChecked()
        if enable:
            fit2DCal = self.ubcalc.detectorCal.getFit2D()
            cp = fit2DCal["centerX"], fit2DCal["centerY"]
            gam_p, _ = self.ubcalc.detectorCal.rangegamdel_p
            azimy, azimx = self.ubcalc.detectorCal.pixelsPrimeBeam(gam_p[1] / 5, 0)[0]
            self.centralPlot.addMarker(
                cp[0],
                cp[1],
                legend="CentralPixel",
                text="CP",
                color="yellow",
                symbol="+",
            )
            self.centralPlot.addMarker(
                azimx, azimy, legend="azimuth", text="Azim", color="yellow", symbol="+"
            )
        else:
            self.centralPlot.removeMarker("CentralPixel")
            self.centralPlot.removeMarker("azimuth")

    def searchPixelCoordHKL(self, hkl):
        """Find detector coordinates and image numbers for a reflection.

        This calculates the two possible diffractometer solutions for ``hkl``
        with :meth:`QUBCalculator.calcReflection`, which uses
        :meth:`HKLVlieg.VliegAngles.anglesZmode`. The angle calculation follows
        the six-circle geometry described by Lohmeier & Vlieg, J. Appl. Cryst.
        26, 706-716 (1993), https://doi.org/10.1107/S0021889893004868.

        Scan-axis mapping:

        * ``mu`` scans: the image number is selected from the solution
          ``alpha`` angle.
        * ``th`` scans: the image number is selected from ``-omega`` because
          the scan axis stores theta with the opposite sign convention used
          internally for omega.

        :param numpy.ndarray hkl:
            Reciprocal-space coordinate in r.l.u. Shape ``(3,)``.
        :returns:
            Reflection dictionary from ``calcReflection`` with added
            ``imageno_1``/``imageno_2`` entries when the candidate angle falls
            within the active scan range, and ``selectable_1``/``selectable_2``
            flags indicating whether ``xy_1``/``xy_2`` lie on the detector.
            The ``angles_1`` and ``angles_2`` rows use
            ``[alpha, delta, gamma, omega, chi, phi]`` in rad; ``xy_1`` and
            ``xy_2`` are detector pixel coordinates as ``(x, y)``.
        :rtype: dict

        .. note::
           CLI-capable for active ``mu`` and ``th`` scans. Unsupported scan
           axes are reported through logging and return ``None``.
        """
        refldict = self.ubcalc.calcReflection(hkl)
        dc = self.ubcalc.detectorCal

        if self.fscan.axisname == "mu":
            angle_idx = 0
            sign = 1.0
        elif self.fscan.axisname == "th":
            angle_idx = 3
            sign = -1.0
        else:
            logger.error(
                f"Cannot calculate reflection. {self.fscan.axisname} is no supported scan axis.",  # noqa: E501
                extra={
                    "title": "Cannot calculate reflection.",
                    "description": f"Cannot calculate reflection. {self.fscan.axisname} is no supported scan axis.",  # noqa: E501
                    "show_dialog": True,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            return
        try:
            imageno1 = self.axisToImageNo(
                np.rad2deg(refldict["angles_1"][angle_idx]) * sign
            )
            refldict["imageno_1"] = imageno1
        except Exception:
            imageno1 = None
        xy = refldict["xy_1"]
        onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and (
            xy[1] >= 0 and xy[1] < dc.detector.shape[0]
        )
        if onDetector:
            refldict["selectable_1"] = True
        else:
            refldict["selectable_1"] = False

        try:
            imageno2 = self.axisToImageNo(
                np.rad2deg(refldict["angles_2"][angle_idx]) * sign
            )
            refldict["imageno_2"] = imageno2
        except Exception:
            imageno2 = None
        xy = refldict["xy_2"]
        onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and (
            xy[1] >= 0 and xy[1] < dc.detector.shape[0]
        )
        if onDetector:
            refldict["selectable_2"] = True
        else:
            refldict["selectable_2"] = False
        return refldict

    def onSearchHKLforStaticROI(self, hkl):
        """GUI/CLI hint: search an hkl and prompt for a static ROI location."""
        try:
            refldict = self.searchPixelCoordHKL(hkl)
        except Exception as e:
            qutils.warning_detailed_message(
                self,
                "Cannot calculate location of reflection",
                f"Cannot calculate position of reflection:\n{e}",
                traceback.format_exc(),
            )
            return
        refl_dialog = QReflectionAnglesDialog(
            refldict, "Select reflection location", self
        )
        if qt.QDialog.Accepted == refl_dialog.exec():
            for i, cb in enumerate(refl_dialog.checkboxes, 1):
                if cb.isChecked():
                    xy = refldict[f"xy_{i}"]
                    self.scanSelector.set_xy_static_loc(xy[0], xy[1])
                    return

    def _onStaticROIedited(self):
        """GUI/CLI hint: synchronize ROI widget edits back to scan controls."""
        xy = self.roiS1.getCenter()
        hsize, vsize = np.round(self.roiS1.getSize())
        self.scanSelector.hsize.blockSignals(True)
        self.scanSelector.vsize.blockSignals(True)
        self.scanSelector.hsize.setValue(hsize)
        self.scanSelector.vsize.setValue(vsize)
        self.scanSelector.hsize.blockSignals(False)
        self.scanSelector.vsize.blockSignals(False)
        self.scanSelector.set_xy_static_loc(xy[0], xy[1])

    def _setCalculatedReflectionImageInfo(self, refldict):
        axisname = self.fscan.axisname
        dc = self.ubcalc.detectorCal

        if axisname == "mu":
            angle_idx = 0
            sign = 1.0
        elif axisname == "th":
            angle_idx = 3
            sign = -1.0
        else:
            qt.QMessageBox.warning(
                self,
                "Cannot calculate reflection",
                f"Cannot calculate reflection.\n{self.fscan.axisname} is no supported scan axis.",  # noqa: E501
            )
            return False

        for intersect in (1, 2):
            try:
                imageno = self.axisToImageNo(
                    np.rad2deg(refldict[f"angles_{intersect}"][angle_idx]) * sign
                )
                refldict[f"imageno_{intersect}"] = imageno
                xy = refldict[f"xy_{intersect}"]
                onDetector = (xy[0] >= 0 and xy[0] < dc.detector.shape[1]) and (
                    xy[1] >= 0 and xy[1] < dc.detector.shape[0]
                )
                refldict[f"selectable_{intersect}"] = onDetector
            except Exception:
                refldict[f"selectable_{intersect}"] = False
        return True

    def onViewCalculatedReflection(self, refldict, intersect):
        """GUI-only: switch to a calculated reflection and center the plot."""
        if not self._setCalculatedReflectionImageInfo(refldict):
            return
        if not refldict.get(f"selectable_{intersect}", False):
            qutils.warning_detailed_message(
                self,
                "Reflection is not on detector",
                "The selected reflection intersect is not on the detector.",
                "",
            )
            return
        imageno = refldict.get(f"imageno_{intersect}")
        if imageno is None:
            qutils.warning_detailed_message(
                self,
                "Reflection image is not in scan",
                "The selected reflection intersect is not in the scan range.",
                "",
            )
            return
        self._onChangeImage(imageno)
        self._onCenterGraph(refldict[f"xy_{intersect}"])

    def _onNewReflection(self, refldict):
        """GUI-only: prompt the user to add calculated reflection candidates."""
        if not self._setCalculatedReflectionImageInfo(refldict):
            return

        refl_dialog = QReflectionAnglesDialog(
            refldict,
            "Select reflections to add into list of reference reflections",
            self,
        )
        if qt.QDialog.Accepted == refl_dialog.exec():
            for i, cb in enumerate(refl_dialog.checkboxes, 1):
                if cb.isChecked():
                    xy = refldict[f"xy_{i}"]
                    eventdict = {"x": xy[0], "y": xy[1]}
                    self.reflectionSel.addReflection(
                        eventdict, refldict[f"imageno_{i}"], refldict["hkl"]
                    )

    def newXyHKLConverter(self):
        """Create a pixel-to-hkl converter bound to the active image state.

        :returns:
            Callable accepting ``x`` and ``y`` detector pixel coordinates.
        :rtype: callable

        .. note::
           CLI-safe. The returned callable reads current scan and UB state.
        """

        def xyToHKL(x, y):
            """CLI-safe: convert detector pixels to hkl and detector angles."""
            # print("xytoHKL:")
            # print("x,y = %s, %s" % (x,y))
            if self.fscan is None:
                return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])
            mu, om = self.getMuOm(self.imageno)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([y]), np.array([x]), mu
            )
            # print(self.ubcalc.detectorCal)
            # print(x,y)
            # print(self.ubcalc.detectorCal.tth(np.array([y]),np.array([x])))
            pos = [mu, delta[0], gamma[0], om, self.ubcalc.chi, self.ubcalc.phi]
            pos = HKLVlieg.crystalAngles(pos, self.ubcalc.n)
            hkl = np.concatenate(
                (
                    np.array(self.ubcalc.angles.anglesToHkl(*pos)),
                    np.rad2deg([delta[0], gamma[0]]),
                )
            )
            return hkl

        return xyToHKL

    def getMuOm(self, imageno=None):
        """Return mu and omega for an image or the whole active scan.

        :param imageno:
            Optional image index. If omitted, arrays for the full scan are
            returned where available.
        :returns:
            ``mu`` and ``omega`` angles in rad.
        :rtype: tuple

        .. note::
           CLI-safe when a scan is loaded.
        """
        if imageno is not None:
            if not self.isValidImageNo(imageno):
                raise IndexError(
                    f"Image number {imageno} is outside the active scan range."
                )
            if self.fscan.axisname == "th":
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.imageNoToAxis(imageno))
            elif self.fscan.axisname == "mu":
                mu = np.deg2rad(self.imageNoToAxis(imageno))
                om = -1 * np.deg2rad(self.fscan.th)
                if len(np.asarray(om).shape) > 0:
                    om = om[0]
            else:
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.fscan.th)
                if len(np.asarray(om).shape) > 0:
                    om = om[0]
            return mu, om
        else:
            if self.fscan.axisname == "th":
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.fscan.axis)
            elif self.fscan.axisname == "mu":
                mu = np.deg2rad(self.fscan.axis)
                om = -1 * np.deg2rad(self.fscan.th)
            else:
                mu = self.ubcalc.mu
                om = -1 * np.deg2rad(self.fscan.th)
            return mu, om

    def omegaToImageNo(self, omega):
        """Map an omega angle to the nearest image index.

        :param float omega:
            Omega angle in rad.
        :returns:
            Nearest image index.
        :rtype: int

        .. note::
           CLI-safe when a scan is loaded.
        """
        if self.fscan is not None:
            omrad = np.deg2rad(self.fscan.omega)
            ommax = np.amax(omrad)
            ommin = np.amin(omrad)
            # print(ommin,omega,ommax)
            if omega < ommin or omega > ommax:
                omdeg = np.rad2deg([ommin, omega, ommax])
                raise Exception(
                    "omega not in range: {} < {} < {}".format(*tuple(omdeg))
                )
            return np.argmin(np.abs(omrad - omega))
        else:
            raise Exception("No Scan selected")

    def isValidImageNo(self, imageno):
        """Return whether ``imageno`` addresses an image in the active scan.

        :param int imageno:
            Image index to validate.
        :returns:
            ``True`` when an active scan exists and ``imageno`` is an integer
            index within ``[0, len(scan))``.
        :rtype: bool

        .. note::
           CLI-safe.
        """
        if self.fscan is None:
            return False
        try:
            image_index = int(imageno)
        except (TypeError, ValueError):
            return False
        return image_index == imageno and 0 <= image_index < len(self.fscan)

    def imageNoToOmega(self, imageno):
        """Return omega for an image index.

        :param int imageno:
            Image index in the active scan.
        :returns:
            Omega angle in rad, ``0.0`` when no scan is loaded.
        :raises IndexError:
            If ``imageno`` is outside the active scan range.

        .. note::
           CLI-safe.
        """
        if self.fscan is not None:
            if not self.isValidImageNo(imageno):
                raise IndexError(
                    f"Image number {imageno} is outside the active scan range."
                )
            return np.deg2rad(self.fscan.omega[int(imageno)])
        else:
            return 0.0

    def imageNoToAxis(self, imageno):
        """Return the scan-axis value for an image index.

        :param int imageno:
            Image index in the active scan.
        :returns:
            Scan-axis value in the scan's stored units, usually deg.
        :raises IndexError:
            If ``imageno`` is outside the active scan range.

        .. note::
           CLI-safe.
        """
        if self.fscan is not None:
            if not self.isValidImageNo(imageno):
                raise IndexError(
                    f"Image number {imageno} is outside the active scan range."
                )
            return self.fscan.axis[int(imageno)]
        else:
            return 0.0

    def axisToImageNo(self, axisval):
        """Map a scan-axis value to the nearest image index.

        :param float axisval:
            Scan-axis value in the scan's stored units, usually deg.
        :returns:
            Nearest image index.
        :rtype: int

        .. note::
           CLI-safe when a scan is loaded.
        """
        if self.fscan is not None:
            # axis = np.deg2rad(self.fscan.axis)
            axismax = np.amax(self.fscan.axis)
            axismin = np.amin(self.fscan.axis)
            # print(ommin,omega,ommax)
            if axisval < axismin or axisval > axismax:
                axisrange = [axismin, axisval, axismax]
                raise Exception(
                    'Value of scan axis "{}" not in range: {} < {} < {}'.format(
                        *tuple([self.fscan.axisname] + axisrange)
                    )
                )
            return np.argmin(np.abs(self.fscan.axis - axisval))
        else:
            raise Exception("No Scan loaded")

    def _onCreateScan(self):
        """GUI-only: create a simulation scan from dialog-entered angles."""
        try:
            mu, om = self.getMuOm(self.imageno)
        except Exception:
            mu = self.ubcalc.mu
            om = 0.0
        th = om * -1.0
        muTh = np.rad2deg([mu, th])  # defaults if fixed
        diag = QScanCreator(muTh)
        if diag.exec() == qt.QDialog.Accepted:
            shape = self.ubcalc.detectorCal.detector.shape
            try:
                axis = diag.scanaxis.currentText()
                if axis == "theta":
                    axis = "th"
                elif axis == "mu":
                    pass
                fscan = SimulationScan(
                    shape,
                    diag.omstart.value(),
                    diag.omend.value(),
                    diag.no.value(),
                    axis,
                    diag.fixedAngle.value(),
                )
                self._onScanChanged(fscan)
            except MemoryError:
                qutils.warning_detailed_message(
                    self,
                    "Can not create simulation scan",
                    "Can not create simualtion scan. Memory is insufficient for the scan size. See details for further information.",  # noqa: E501
                    traceback.format_exc(),
                )

    def _onLoadInterlacedScan(self):
        """GUI-only: build an interlaced scan from selected HDF5 scans."""
        # GUI function to concatenate multiple scans into one
        # a user can select which scans to combine in a GUI dialog
        # uses new class interlacedScan

        # grab the selected h5 file / node from the tree in the scan selector tab
        model = self.scanSelector.hdfTreeView.model()
        selection = self.scanSelector.hdfTreeView.selectionModel()
        indexes = selection.selectedIndexes()
        if indexes == []:
            qutils.warning_detailed_message(
                self,
                "Can not create interlaced scan",
                "Can not create interlaced scan. Select a node in the tree view first!",
                traceback.format_exc(),
            )
            return
        rootI = indexes.pop(0)

        # address the root node to correctly get the scan names
        if rootI.parent().isValid():
            h5file = model.data(
                model.parent(rootI), role=silx.gui.hdf5.Hdf5TreeModel.H5PY_OBJECT_ROLE
            )
        else:
            h5file = model.data(
                rootI, role=silx.gui.hdf5.Hdf5TreeModel.H5PY_OBJECT_ROLE
            )

        isID31 = self.scanSelector.btid.currentText() in [
            "ch5523",
            "ch5700",
            "ch5918",
            "ch6392",
            "ch7131",
            "ch7149",
            "ch7856",
            "ch8153",
            "id31_default",
        ]
        kl_full = list(h5file.keys())
        kl = np.empty(0, dtype=int)
        for i in kl_full:
            if isID31:
                pattern = r"\.\d"
                result = re.findall(pattern, i)[0][1:]
                if result == "1":
                    # select only scan names which are ending on suffix '.1' (fast counters of id31 hdf5 format)  # noqa: E501
                    kl = np.append(kl, i)
            else:
                kl = np.append(kl, i)

        # separate scan nr and delete duplicates suffixes
        if isID31:
            # try to get the scan nr and '/title' from the hdf5 file
            nr = np.empty(0, dtype=int)
            name = np.empty(0, dtype=str)

            for i in kl:
                pattern = r"\d+\."
                result = re.findall(pattern, i)
                name = np.append(name, h5file[i + "/title"])
                nr = np.append(nr, int(result[0][:-1]))

            lsort = np.argsort(nr)[::1]
            nr = nr[lsort]
            name = name[lsort]
        else:
            nr = np.empty(0, dtype=int)
            name = np.empty(0, dtype=str)
            for nth, i in enumerate(kl):
                name = np.append(name, i)
                nr = np.append(
                    nr, nth + 1
                )  # create scan nr list with ascending integers, starting with 1
                # This will later be used to address the subscans, so check if your scans are handled like this!!!  # noqa: E501

        # open GUI dialog to select which scans to combine
        interlacedSelectDialog = qt.QDialog()

        llayout = qt.QGridLayout()
        llayout.addWidget(qt.QLabel("Available scans:"), 0, 0)

        a = qt.QScrollArea()
        b = qt.QFormLayout()
        box = qt.QGroupBox()
        scanBoxes = []
        # labels = []
        for i, item in enumerate(nr):
            ithScanBox = qt.QCheckBox()
            scanBoxes.append(ithScanBox)
            # labels.append(qt.QLabel('Scan '+str(item)+':'+str(name[i][2:-1])))
            if isID31:
                b.addRow(
                    qt.QLabel("Scan " + str(item) + ": " + name[i].decode()), ithScanBox
                )
            else:
                b.addRow(qt.QLabel("Scan " + str(item) + ": " + name[i]), ithScanBox)

        box.setLayout(b)
        a.setWidget(box)
        a.setWidgetResizable(True)
        llayout.addWidget(a, 1, 0, 1, -1)

        llayout.addWidget(qt.QLabel("sort the scans by axis values?"), 2, 0)

        noScans = qt.QCheckBox()
        llayout.addWidget(noScans, 2, 1)

        llayout.addWidget(qt.QLabel("Backend:"), 3, 0)

        IS_btid = qt.QComboBox()
        [IS_btid.addItem(bt) for bt in backends.fscans]
        IS_btid.setCurrentText(self.scanSelector.btid.currentText())
        llayout.addWidget(IS_btid, 3, 1)

        llayout.addWidget(qt.QLabel("scan axis:"), 4, 0)

        axisbox = qt.QComboBox()
        [axisbox.addItem(a) for a in ["th", "mu"]]
        axisbox.setCurrentText("th")
        llayout.addWidget(axisbox, 4, 1)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(
            interlacedSelectDialog.accept
        )
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(
            interlacedSelectDialog.reject
        )

        llayout.addWidget(buttons, 5, 0, -1, -1)
        interlacedSelectDialog.setLayout(llayout)
        interlacedSelectDialog.setWindowTitle("Segmented scan loader")

        if not interlacedSelectDialog.exec() == qt.QDialog.Accepted:
            return

        # generate scan objects for selected scans
        selectedScans = []
        for i, j in enumerate(scanBoxes):
            if j.isChecked():
                selectedScans.append(nr[i])
                # selectedScans.append(name[i])

        nodes = list(self.scanSelector.hdfTreeView.selectedH5Nodes())
        obj = nodes[0]

        scansegments = []
        for i in selectedScans:
            ddict = dict()
            ddict["scanno"] = int(i)
            ddict["file"] = obj.local_filename
            # ddict['node'] = kl[i]
            ddict["beamtime"] = IS_btid.currentText()
            try:
                scansegments.append(backends.openScan(IS_btid.currentText(), ddict))
            except Exception as e:
                msg = qt.QMessageBox()
                msg.setIcon(qt.QMessageBox.Warning)
                msg.setWindowTitle("Cannot open scan")
                msg.setText(f"Cannot open scan:\n{e}\nDo you want to continue?")
                msg.setDetailedText(traceback.format_exc())
                msg.setStandardButtons(qt.QMessageBox.Yes | qt.QMessageBox.No)
                msg.setDefaultButton(qt.QMessageBox.Yes)
                result = msg.exec()
                if result != qt.QMessageBox.Yes:
                    return
        if not scansegments:  # no scans loaded - abort.
            qt.QMessageBox.critical(self, "No scans loaded", "No scans were loaded.")
            return

        # create interlaced scan object
        self.scanno = 1
        self.fscan = interlacedScanLoader.InterlacedScan(
            scansegments, noScans.isChecked(), axisbox.currentText()
        )
        self.imageno = 0
        self.plotImage()
        self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
        self.activescanname = "{}-segmentedScan {} {}-{}".format(
            self.fscan.axisname,
            ",".join(str(itemNr) for itemNr in selectedScans),
            np.amin(self.fscan.axis),
            np.amax(self.fscan.axis),
        )

        # generate sum and max image
        self.images_loaded = False
        if self.fscan is not None and self.autoLoadAct.isChecked():
            self.loadAll()
            self.scanSelector.showMaxAct.setChecked(False)
            self.scanSelector.showMaxAct.setChecked(True)

    def _onLoadScanFromImages(self):
        """GUI-only: import detector images through file and setup dialogs."""
        # generates a scan from a selected folder containing raw detector images

        # generate file source selection GUI

        # create filter of scan image formats (following code is copied from silx view)
        extensions = {}
        for description, ext in silx.io.supported_extensions().items():
            extensions[description] = " ".join(sorted(list(ext)))

        extensions["NeXus layout from EDF files"] = "*.edf"
        extensions["NeXus layout from TIFF image files"] = "*.tif *.tiff"
        extensions["NeXus layout from CBF files"] = "*.cbf"
        extensions["NeXus layout from MarCCD image files"] = "*.mccd"

        all_supported_extensions = set()
        for name, exts in extensions.items():
            exts = exts.split(" ")
            all_supported_extensions.update(exts)
        all_supported_extensions = sorted(list(all_supported_extensions))

        filters = []
        filters.append(
            "All supported files ({})".format(" ".join(all_supported_extensions))
        )
        for name, extension in extensions.items():
            filters.append(f"{name} ({extension})")
        filters.append("All files (*)")

        fileTypeFilter = ""
        for f in filters:
            fileTypeFilter += f + ";;"

        # call dialog
        filename, _ = qt.QFileDialog.getOpenFileName(
            self, "Open data source", "", fileTypeFilter[:-2]
        )

        # Qt dialog returns '' if cancelled
        if filename == "":
            qt.QMessageBox.warning(
                self, "Error - Open data source", "No data source selected"
            )
            return

        # search files using ImportImagesScan backend
        importedscan = universalScanLoader.ImportImagesScan(filename)

        if importedscan.inpath is None:
            qt.QMessageBox.critical(
                self,
                "Images could not be imported",
                "The selected data source is not suitable\n"
                "It is necessary to select a file containing raw detector image(s)!",
            )
            return

        if importedscan.shape != self.ubcalc.detectorCal.detector.shape:
            qt.QMessageBox.critical(
                self,
                "Detector data mismatch",
                "The selected image data shape does not match to the detector data shape:\n"  # noqa: E501
                "Detector Size {}x{}\n"
                "Data size {}x{}\n"
                "Please first adjust the detector configuration to load this data".format(  # noqa: E501
                    *self.ubcalc.detectorCal.detector.shape, *importedscan.shape
                ),
            )
            return

        [imagePrefix, found_scanfiles] = importedscan.inpath

        # generate dialog with list of files and frames
        nrofFilesfound = len(found_scanfiles)
        messageStr = "Found " + str(nrofFilesfound) + " files in selected directory"

        if importedscan.FramesPerFile > 1:
            if nrofFilesfound == 0:
                messageStr = "No images found!!!"
                fullStr = messageStr
            elif 0 < nrofFilesfound < 4:
                messageStr += ":\n"
                for i in range(nrofFilesfound - 1):
                    messageStr += (
                        imagePrefix
                        + found_scanfiles[i]
                        + ": "
                        + str(importedscan.FramesPerFile)
                        + " frames"
                    )  ### mark expected nr when file not actually loaded
                    if i > 0:
                        messageStr += " (expected)"
                    messageStr += "\n"
                messageStr += (
                    imagePrefix
                    + found_scanfiles[nrofFilesfound - 1]
                    + ": "
                    + str(importedscan.FramesLastFile)
                    + " frames\n"
                    + str(importedscan.nopoints)
                    + " frames in total."
                )
                fullStr = messageStr
            else:
                messageStr += ":\n"
                for i in range(0, 3):
                    messageStr += (
                        imagePrefix
                        + found_scanfiles[i]
                        + ": "
                        + str(importedscan.FramesPerFile)
                        + " frames"
                    )
                    if i > 0:
                        messageStr += " (expected)"
                    messageStr += "\n"
                fullStr = messageStr
                for i in range(3, nrofFilesfound - 1):
                    fullStr += (
                        imagePrefix
                        + found_scanfiles[i]
                        + ": "
                        + str(importedscan.FramesPerFile)
                        + " frames (expected) \n"
                    )

                messageStr += (
                    "..."
                    + "\n"
                    + imagePrefix
                    + found_scanfiles[nrofFilesfound - 1]
                    + ": "
                    + str(importedscan.FramesLastFile)
                    + " frames\n"
                    + str(importedscan.nopoints)
                    + " frames in total."
                )
                fullStr += (
                    imagePrefix
                    + found_scanfiles[nrofFilesfound - 1]
                    + ": "
                    + str(importedscan.FramesLastFile)
                    + " frames\n"
                    + str(importedscan.nopoints)
                    + " frames in total."
                )

        else:
            if nrofFilesfound == 0:
                messageStr = "No images found!!!"
                fullStr = messageStr
            elif 0 < nrofFilesfound < 4:
                messageStr += ":\n"
                for i in range(nrofFilesfound - 1):
                    messageStr += imagePrefix + found_scanfiles[i] + "\n"
                messageStr += imagePrefix + found_scanfiles[nrofFilesfound - 1]
                fullStr = messageStr
            else:
                messageStr += ":\n"
                for i in range(0, 3):
                    messageStr += imagePrefix + found_scanfiles[i] + "\n"
                fullStr = messageStr
                for i in range(3, nrofFilesfound):
                    fullStr += imagePrefix + found_scanfiles[i] + "\n"
                messageStr += (
                    "..." + "\n" + imagePrefix + found_scanfiles[nrofFilesfound - 1]
                )
                fullStr += "\n" + str(importedscan.nopoints) + " frames in total."

        msg0 = qt.QMessageBox(self)
        msg0.setWindowTitle("Manual scan import")
        msg0.setText(messageStr)
        msg0.setDetailedText(fullStr)
        msg0.exec()

        # angle conversions
        try:
            mu, om = self.getMuOm(self.imageno)
        except Exception:
            mu = self.ubcalc.mu
            om = 0.0
        th = om * -1.0
        muTh = np.rad2deg([mu, th])  # defaults if fixed

        # open scan creator GUI to let the user insert missing scan angles
        diag = QImportScanCreator(muTh)
        # detector pixel nr and frame nr is adapted from opened image file
        diag.no.setValue(importedscan.nopoints)

        if diag.exec() == qt.QDialog.Accepted:
            try:
                axis = diag.scanaxis.currentText()
                if axis == "theta":
                    axis = "th"
                elif axis == "mu":
                    pass
                # pass inserted angles to scan object
                importedscan.set_axis(
                    diag.omstart.value(),
                    diag.omend.value(),
                    axis,
                    diag.fixedAngle.value(),
                )
                self._onScanChanged(importedscan)
            except MemoryError:
                qutils.warning_detailed_message(
                    self,
                    "Can not create scan",
                    "Can not create scan. Memory is insufficient for the scan size. See details for further information.",  # noqa: E501
                    traceback.format_exc(),
                )

    def _onScanChanged(self, sel_list):
        """CLI-capable: load or activate the selected scan object/list."""
        self.resetZoom = True
        # print(sel_list)
        self.activescanname = "scan"
        if isinstance(sel_list, list):
            self.sel_list = sel_list
            if len(sel_list):
                self.specfile = sel_list[0]["SourceName"]
                try:
                    self.scanno = int(float(sel_list[0]["Key"])) - 1
                    self.fscan = Fastscan(self.specfile, self.scanno)
                    self.imageno = 0
                except Exception:
                    self.scanno = 0
                    self.fscan = FioFastsweep(self.specfile)
                    self.imageno = 0
                self.reflectionSel.setImage(self.imageno)
                if self.imagepath != "":
                    self.fscan.set_image_folder(self.imagepath)
                    self.plotImage()
                    self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)

            else:
                self.scanSelector.setRange(0, 0)
                self.imageno = 0
                self.reflectionSel.setImage(self.imageno)
                # print(self.centralPlot._callback)

        elif isinstance(sel_list, universalScanLoader.ImportImagesScan):
            self.scanno = 1
            self.fscan = sel_list
            self.imageno = 0
            self.plotImage()
            self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
            self.activescanname = f"{self.fscan.axisname}-rawImport {np.amin(self.fscan.axis)}-{np.amax(self.fscan.axis)}"  # noqa: E501

            self.images_loaded = False
            if self.fscan is not None and self.autoLoadAct.isChecked():
                self.loadAll()
                self.scanSelector.showMaxAct.setChecked(False)
                self.scanSelector.showMaxAct.setChecked(True)

        elif isinstance(sel_list, SimulationScan):
            self.scanno = 1
            self.fscan = sel_list
            self.imageno = 0
            self.plotImage()
            self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
            self.activescanname = f"{self.fscan.axisname}-sim {np.amin(self.fscan.axis)}-{np.amax(self.fscan.axis)}"  # noqa: E501
        else:
            if "name" in sel_list:
                self.activescanname = sel_list["name"]
            else:
                self.activescanname = "scan"
            self.hdffile = sel_list["file"]
            # self.scanname = sel_list['name'].strip("/")
            try:
                logger.info("Loading scan...")
                if logger_utils.get_logging_context() == "gui":
                    msg = qt.QMessageBox(self)
                    msg.setWindowTitle("Loading Scan")
                    msg.setText("Loading Scan. This might take a while...")
                    msg.setStandardButtons(qt.QMessageBox.Cancel)
                    msg.setModal(True)
                    msg.show()

                if "beamtime" in sel_list:
                    self.fscan = backends.openScan(sel_list["beamtime"], sel_list)
                else:
                    self.fscan = backends.openScan(
                        self.scanSelector.btid.currentText(), sel_list
                    )

                self.plotImage()
                self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
                if logger_utils.get_logging_context() == "gui":
                    msg.hide()
                self.images_loaded = False
                if self.fscan is not None and self.autoLoadAct.isChecked():
                    self.loadAll()
                    self.scanSelector.showMaxAct.setChecked(False)
                    self.scanSelector.showMaxAct.setChecked(True)
            except Exception as exc:
                if logger_utils.get_logging_context() == "gui":
                    msg.hide()
                exc_msg = str(exc) or exc.__class__.__name__
                logger.exception(
                    "Cannot open scan",
                    extra={
                        "title": "Cannot open scan",
                        "description": f"Cannot open scan:\n{exc_msg}",
                        "show_dialog": True,
                        "parent": self,
                    },
                )
                # qutils.warning_detailed_message(self, "Cannot open scan", "Cannot open scan" , traceback.format_exc())  # noqa: E501
                # qt.QMessageBox.critical(self,"Cannot open scan", "Cannot open scan:\n%s" % traceback.format_exc())  # noqa: E501
        if hasattr(self.fscan, "name"):
            self.activescanname = self.fscan.name

    def _onImagePathChanged(self, path):
        """CLI-capable: update the image folder for the active scan."""
        # print("newpath %s" % path)
        self.imagepath = path
        if self.fscan is not None:
            self.fscan.set_image_folder(self.imagepath)
            self.plotImage()
            self.scanSelector.setAxis(self.fscan.axis, self.fscan.axisname)
            # self.scanSelector.slider.setMinimum(0)
            # self.scanSelector.slider.setMaximum(self.fscan.nopoints-1)
        else:
            self.scanSelector.setRange(0, 0)
            self.imageno = 0
            self.reflectionSel.setImage(self.imageno)
            # print(self.centralPlot._callback)

    def _onChangeImage(self, imageno):
        """GUI/CLI hint: switch the active image and replot it."""
        if self.fscan is not None:
            if not self.isValidImageNo(imageno):
                logger.warning(
                    "Skipping request to display stale image number %s outside "
                    "the active scan range.",
                    imageno,
                )
                return
            self.scanSelector.slider.setValue(imageno)
            self.plotImage(self.scanSelector.slider.value())

    def _onSliderValueChanged(self, value):
        """GUI/CLI hint: replot the image selected by the scan slider."""
        if self.fscan is not None:
            self.plotImage(value)
        # print(self.centralPlot._callback)

    def _onLoadAll(self):
        """GUI/CLI hint: reload all images and refresh max-image display."""
        self.images_loaded = False
        if self.fscan is not None:
            self.loadAll()
            self.scanSelector.showMaxAct.setChecked(False)
            self.scanSelector.showMaxAct.setChecked(True)

    def loadAll(self):
        """Load all scan images and compute summed and maximum images.

        :returns:
            ``None``. Results are stored on ``self.allimgsum`` and
            ``self.allimgmax``.

        .. note::
           CLI-capable. Progress reporting is routed through ``logger_utils``.
        """
        try:
            image = self.fscan.get_raw_img(0)
        except Exception:
            logger.exception(
                "No image found in scan. (image 0 is missing)",
                extra={
                    "title": "No image found in scan.",
                    "description": "No image found in scan. (image 0 is missing)",
                    "show_dialog": False,
                    "dialog_level": logging.ERROR,
                    "parent": self,
                },
            )
            return
        self.allimgsum = np.zeros_like(image.img, dtype=np.float64)
        self.allimgmax = np.zeros_like(image.img, dtype=np.float64)

        progress = logger_utils.create_progress_logger(
            self, len(self.fscan), "Reading images"
        )

        img_size = self.allimgsum.nbytes / (1024**2)

        min(np.floor((self.maxMemory - 2000) / (img_size * self.numberthreads)), 5)
        np.arange(len(self.fscan))
        # self.excludedImagesDialog.getData()

        lock = threading.Lock()

        self.images_loaded = True
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.numberthreads
        ) as executor:  # speedup only for the file reads
            futures = {}
            excl = self.excludedImagesDialog.getData()
            bg = self.background_image

            def readfile_max(imgno):
                """CLI-safe worker: read one image into sum and max buffers."""
                if imgno in excl:  # skip if excluded
                    return imgno
                image = self.fscan.get_raw_img(imgno).img.astype(
                    np.float64, order="C", copy=True
                )
                if bg is not None and bg.shape == image.shape:
                    if HAS_ACCEL:
                        _roi_sum_accel.calcBgSub(image, bg)
                    else:
                        np.subtract(image, bg, out=image)
                with lock:
                    self.allimgsum += image
                    np.maximum(self.allimgmax, image, out=self.allimgmax)
                return imgno

            for i in range(len(self.fscan)):
                futures[executor.submit(readfile_max, i)] = i

            for f in concurrent.futures.as_completed(futures):
                try:
                    imgno = f.result()
                    progress.update(imgno)
                except concurrent.futures.CancelledError:
                    pass

                except Exception:
                    logger.warn(f"Cannot read image:\n{traceback.format_exc()}")
                    # print("Cannot read image:\n%s" % traceback.format_exc())

                if progress.wasCanceled():
                    [f.cancel() for f in futures]
                    self.images_loaded = False
                    logger.warn(
                        "Loading of images cancelled. Max and Sum images are incomplete."  # noqa: E501
                    )
                    break
        progress.finish()

    def _onMaxToggled(self, value):
        """GUI/CLI hint: toggle display of the precomputed maximum image."""
        if self.scanSelector.showSumAct.isChecked():
            self.scanSelector.showSumAct.setChecked(False)
        if value:
            if not self.images_loaded and self.fscan is not None:
                if logger_utils.get_logging_context() == "gui":
                    btn = qt.QMessageBox.question(
                        self,
                        "Incomplete sum / max image",
                        "Sum/Max image was not loaded completely. Displayed maximum image will be incomplete! Do you want to load all images?",  # noqa: E501
                        qt.QMessageBox.Yes | qt.QMessageBox.No | qt.QMessageBox.Cancel,
                    )
                    if btn == qt.QMessageBox.Yes:
                        self.loadAll()
                    elif btn == qt.QMessageBox.Cancel:
                        self.scanSelector.showMaxAct.setChecked(False)
                        return
                else:
                    logger.warning(
                        "Maximum image requested before all images were loaded; continuing with incomplete data."  # noqa: E501
                    )
            if self.allimgmax is not None:
                self.currentAddImageLabel = self.centralPlot.addImage(
                    self.allimgmax,
                    legend="special",
                    replace=False,
                    resetzoom=False,
                    copy=True,
                    z=1,
                )
                self.centralPlot.setActiveImage(self.currentAddImageLabel)
                self.scanSelector.alphaslider.setLegend(self.currentAddImageLabel)
            else:
                self.scanSelector.showMaxAct.setChecked(False)
        else:
            if self.currentAddImageLabel is not None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
                self.centralPlot.removeImage(self.currentAddImageLabel)
                self.currentAddImageLabel = None

    def _onSumToggled(self, value):
        """GUI/CLI hint: toggle display of the precomputed summed image."""
        if self.scanSelector.showMaxAct.isChecked():
            self.scanSelector.showMaxAct.setChecked(False)
        if value:
            if not self.images_loaded and self.fscan is not None:
                if logger_utils.get_logging_context() == "gui":
                    btn = qt.QMessageBox.question(
                        self,
                        "Incomplete sum / max image",
                        "Sum/Max image was not loaded completely. Displayed sum image will be incomplete! Do you want to load all images?",  # noqa: E501
                        qt.QMessageBox.Yes | qt.QMessageBox.No | qt.QMessageBox.Cancel,
                    )
                    if btn == qt.QMessageBox.Yes:
                        self.loadAll()
                    elif btn == qt.QMessageBox.Cancel:
                        self.scanSelector.showSumAct.setChecked(False)
                        return
                else:
                    logger.warning(
                        "Summed image requested before all images were loaded; continuing with incomplete data."  # noqa: E501
                    )
            if self.allimgsum is not None:
                self.currentAddImageLabel = self.centralPlot.addImage(
                    self.allimgsum,
                    legend="special",
                    replace=False,
                    resetzoom=False,
                    copy=True,
                    z=1,
                )
                self.centralPlot.setActiveImage(self.currentAddImageLabel)
                self.scanSelector.alphaslider.setLegend(self.currentAddImageLabel)
            else:
                self.scanSelector.showSumAct.setChecked(False)
        else:
            if self.currentAddImageLabel is not None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
                self.centralPlot.removeImage(self.currentAddImageLabel)
                self.currentAddImageLabel = None

    def _roi_preview_enabled(self):
        """Return whether fitted-background ROI preview should alter the image."""
        return (
            HAS_ACCEL
            and self.roivisible
            and getattr(self, "showInterpolatedBgAct", None) is not None
            and self.showInterpolatedBgAct.isChecked()
        )

    def _set_center_roi_preview_color(self, roi):
        """Color a center ROI according to fitted-background preview state."""
        roi.setColor("blue" if self._roi_preview_enabled() else "red")
        roi.setBgStyle("pink", "-", roi.getLineWidth())

    def _roi_array_from_key(self, key):
        """Convert one ROI slice key to the compiled ``(1, 2, 2)`` format."""
        roi = np.array(
            [[[key[0].start, key[0].stop], [key[1].start, key[1].stop]]],
            dtype=np.int64,
        )
        return np.ascontiguousarray(roi)

    def _current_preview_mask(self, image_shape):
        """Return the detector mask used for fitted-background preview."""
        if not self.scanSelector.useMaskBox.isChecked():
            return np.zeros(image_shape, dtype=bool)
        mask_widget = self.centralPlot.getMaskToolsDockWidget()
        mask = mask_widget.getSelectionMask()
        if mask is None or mask.shape != image_shape:
            return np.zeros(image_shape, dtype=bool)
        return np.ascontiguousarray(mask > 0.0, dtype=bool)

    def _apply_interpolated_bg_patch(self, image, mask, ckey, bgkeys):
        """Overwrite one center ROI in ``image`` with fitted background."""
        roioptions = self.scanSelector.roioptions.get_parameters()
        fit_order = int(roioptions.get("FittedBackgroundOrder", 1))
        patch, stats = _roi_sum_accel.interpolate_polybg_croi(
            image,
            mask,
            self._roi_array_from_key(ckey),
            self._roi_array_from_key(bgkeys[0]),
            self._roi_array_from_key(bgkeys[1]),
            self._roi_array_from_key(bgkeys[2]),
            self._roi_array_from_key(bgkeys[3]),
            fit_order,
        )
        if stats["success"]:
            image[ckey[::-1]] = patch
        return stats

    def _apply_interpolated_bg_preview(self, image, image_no):
        """Overwrite visible center ROIs with local polynomial backgrounds."""
        if not self._roi_preview_enabled():
            return

        mask = self._current_preview_mask(image.shape)
        current_mode = self.scanSelector.scanstab.currentIndex()
        if current_mode in [0, 1]:
            hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(image_no)
            for hkl_del_gam in [hkl_del_gam_1, hkl_del_gam_2]:
                if hkl_del_gam[0, -1]:
                    loc = hkl_del_gam[0, 6:8]
                    self._apply_interpolated_bg_patch(
                        image,
                        mask,
                        self.intkey(loc),
                        self.bkgkeys(loc),
                    )
        elif current_mode == 2:
            refldict = self.get_rocking_coordinates()
            self._apply_interpolated_bg_preview_for_rocking(image, mask, refldict)
        elif current_mode == 3:
            refldict = self.get_Bragg_rocking_coordinates()
            if len(refldict["xy_1"]) > 0:
                roi_keys = self.intbkgkeys_rocking(
                    refldict,
                    autovsize=False,
                    autohsize=False,
                    intersect=1,
                )
                self._apply_interpolated_bg_preview_for_keys(image, mask, roi_keys)

    def _apply_interpolated_bg_preview_for_rocking(self, image, mask, refldict):
        """Overwrite displayed rocking center ROIs with fitted backgrounds."""
        roi_keys = self.intbkgkeys_rocking(refldict)
        self._apply_interpolated_bg_preview_for_keys(image, mask, roi_keys)

    def _apply_interpolated_bg_preview_for_keys(self, image, mask, roi_keys):
        """Overwrite the same subset of center ROIs shown by ROI graphics."""
        number_rois = len(roi_keys["center"])
        divider = 1
        if number_rois > self.maxROIs:
            divider = np.ceil(number_rois / self.maxROIs)
        no_rois_to_display = int(np.floor(number_rois / divider))
        for i in np.arange(no_rois_to_display) * divider:
            index = int(i)
            self._apply_interpolated_bg_patch(
                image,
                mask,
                roi_keys["center"][index],
                (
                    roi_keys["left"][index],
                    roi_keys["right"][index],
                    roi_keys["top"][index],
                    roi_keys["bottom"][index],
                ),
            )

    def plotImage(self, key=0):
        """Plot one raw image from the active scan.

        :param int key:
            Image index in the active scan.

        .. note::
           CLI-capable, but it mutates Qt plot and reflection widgets.
        """
        try:
            image = self.fscan.get_raw_img(key).img.astype(
                np.float64, order="C", copy=True
            )
            bg = self.background_image
            if bg is not None and bg.shape == image.shape:
                if HAS_ACCEL:
                    _roi_sum_accel.calcBgSub(image, bg)
                else:
                    np.subtract(image, bg, out=image)
            try:
                self._apply_interpolated_bg_preview(image, key)
            except Exception:
                logger.warning(
                    "Cannot apply interpolated background preview",
                    exc_info=True,
                )
            # if self.currentImageLabel is not None:
            #    self.centralPlot.removeImage(self.currentImageLabel)

            self.currentImageLabel = self.centralPlot.addImage(
                image,
                legend="scan_image",
                replace=False,
                resetzoom=self.resetZoom,
                copy=True,
            )
            if self.currentAddImageLabel is None:
                self.centralPlot.setActiveImage(self.currentImageLabel)
            self.resetZoom = False
            self.imageno = key
            self.reflectionSel.setImage(self.imageno)
            self.updateROI(image_changed=True)
            self._centerTrackedROI(key)
            self.updateReflections()

            mu, om = self.getMuOm(self.imageno)
            self.ubcalc.uedit.setAngles(mu, self.ubcalc.chi, self.ubcalc.phi, om)

            self.scanSelector.excludeImageAct.blockSignals(True)
            self.scanSelector.excludeImageAct.setChecked(
                key in self.excludedImagesDialog.getData()
            )
            self.scanSelector.excludeImageAct.blockSignals(False)

        except Exception:
            logger.exception(
                "Cannot plot image",
                extra={
                    "title": "Cannot plot image",
                    "description": "Cannot plot image",
                    "show_dialog": False,
                    "dialog_level": logging.ERROR,
                    "parent": self,
                },
            )

    def updateReflections(self):
        """Update CTR reflection markers for the current image.

        .. note::
           CLI-capable, but it mutates Qt plot widgets.
        """
        if not self.reflectionsVisible:
            self.centralPlot.removeCurve("all_image_reflections")
            return
        H_0, H_1 = self.reflectionsToDisplay
        # H_0 = np.array([[1,0,0], [1,1,0]])
        # H_1 = np.array([[0,0,1], [0,0,1]])

        hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(
            self.imageno, H_0, H_1, intersect=True
        )

        mask1 = hkl_del_gam_1[:, -1].nonzero()
        mask2 = hkl_del_gam_2[:, -1].nonzero()

        masked_hkl_del_gam = np.vstack((hkl_del_gam_1[mask1], hkl_del_gam_2[mask2]))

        self.centralPlot.addCurve(
            masked_hkl_del_gam[:, -3],
            masked_hkl_del_gam[:, -2],
            legend="all_image_reflections",
            linestyle=" ",
            symbol=".",
            color="y",
            resetzoom=False,
        )

    def updateROI(self, **kwargs):
        """Update visible ROI graphics for the current scan mode.

        .. note::
           CLI-capable, but it mutates Qt ROI widgets and reads GUI controls.
        """
        if not self.roivisible:
            # for roi in self.rois:
            self.roiS1.setVisible(False)
            self.roiS2.setVisible(False)
            if self.rocking_rois:
                for roi in self.rocking_rois:
                    roi.setVisible(False)
                    roi.setEditable(False)
            self.roiManager._roisUpdated()
            return
            # self.centralPlot.removeMarker('main_croi_loc')

        current_mode = self.scanSelector.scanstab.currentIndex()
        if current_mode == 0 or current_mode == 1:
            if self.rocking_rois:
                for roi in self.rocking_rois:
                    roi.setVisible(False)
                    roi.setEditable(False)
            try:
                hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(self.imageno)
            except Exception:
                # print(traceback.format_exc())
                # for roi in self.rois:
                #    roi.setVisible(False)
                self.roiS1.setVisible(False)
                self.roiS2.setVisible(False)
                # self.centralPlot.removeMarker('main_croi_loc')
                return

            if hkl_del_gam_1[0, -1]:
                if self.roivisible:
                    self.plotROI(hkl_del_gam_1[0, 6:8], self.roiS1)
                for i, spinbox in enumerate(self.scanSelector.hkl_static):
                    spinbox.setValue(hkl_del_gam_1[0, i])
            else:
                self.roiS1.setVisible(False)

            if hkl_del_gam_2[0, -1]:
                if self.roivisible:
                    self.plotROI(hkl_del_gam_2[0, 6:8], self.roiS2)
                for i, spinbox in enumerate(self.scanSelector.hkl_static):
                    spinbox.setValue(hkl_del_gam_1[0, i])
            else:
                self.roiS2.setVisible(False)

            if current_mode == 1:
                self.roiS1.setEditable(True)
            else:
                self.roiS1.setEditable(False)

        elif current_mode == 2 and not kwargs.get("image_changed", False):
            self.roiS1.setVisible(False)
            self.roiS2.setVisible(False)
            try:
                refldict = self.get_rocking_coordinates()
            except Exception:
                # print(traceback.format_exc())
                if self.rocking_rois:
                    for roi in self.rocking_rois:
                        roi.setVisible(False)
                        roi.setEditable(False)
                return
            roi_keys = self.intbkgkeys_rocking(refldict)
            self.scanSelector.autoSize_label.setText(
                "{} x {}".format(roi_keys["hsize"], roi_keys["vsize"])
            )

            number_rois = len(roi_keys["center"])
            divider = 1
            if number_rois > self.maxROIs:
                divider = np.ceil(number_rois / self.maxROIs)
            no_rois_to_display = int(np.floor(number_rois / divider))

            # lazy create ROIs
            if len(self.rocking_rois) < no_rois_to_display:
                for i in range(no_rois_to_display - len(self.rocking_rois)):
                    roi = RectangleBgROI()
                    roi.setLineWidth(1)
                    roi.setLineStyle("-")
                    roi.setColor("red")
                    roi.setBgStyle("pink", "-", 1.0)
                    roi.setVisible(False)
                    roi.setGeometry(origin=(0, 0), size=(0, 0))
                    self.rocking_rois.append(roi)
                    self.roiManager.addRoi(roi, useManagerColor=False)

            for roino, i in enumerate(np.arange(no_rois_to_display) * divider):
                self._set_center_roi_preview_color(self.rocking_rois[roino])
                ckey = roi_keys["center"][int(i)]
                leftkey = roi_keys["left"][int(i)]
                rightkey = roi_keys["right"][int(i)]
                topkey = roi_keys["top"][int(i)]
                bottomkey = roi_keys["bottom"][int(i)]

                origin, size, left, right, top, bottom = _display_roi_geometry(
                    ckey, leftkey, rightkey, topkey, bottomkey
                )
                self.rocking_rois[roino].setGeometry(
                    origin=origin,
                    size=size,
                    left=left,
                    right=right,
                    top=top,
                    bottom=bottom,
                )
                self.rocking_rois[roino].setVisible(True)
                self.rocking_rois[roino].setEditable(False)
            for roi in self.rocking_rois[no_rois_to_display:]:
                roi.setVisible(False)
                roi.setEditable(False)

        elif current_mode == 3 and not kwargs.get("image_changed", False):
            self.roiS1.setVisible(False)
            self.roiS2.setVisible(False)
            try:
                refldict = self.get_Bragg_rocking_coordinates()
                if len(refldict["xy_1"]) <= 0:
                    raise Exception("No reflections found.")
            except Exception:
                # print(traceback.format_exc())
                if self.rocking_rois:
                    for roi in self.rocking_rois:
                        roi.setVisible(False)
                        roi.setEditable(False)
                logger.exception(
                    "Cannot calculate Bragg rocking coordinates",
                    extra={
                        "title": "Cannot calculate Bragg rocking coordinates",
                        "description": "Cannot calculate Bragg rocking coordinates",
                        "show_dialog": False,
                        "dialog_level": logging.ERROR,
                        "parent": self,
                    },
                )
                return
            roi_keys = self.intbkgkeys_rocking(
                refldict, autovsize=False, autohsize=False, intersect=1
            )

            number_rois = len(roi_keys["center"])
            divider = 1
            if number_rois > self.maxROIs:
                divider = np.ceil(number_rois / self.maxROIs)
            no_rois_to_display = int(np.floor(number_rois / divider))

            # lazy create ROIs
            if len(self.rocking_rois) < no_rois_to_display:
                for i in range(no_rois_to_display - len(self.rocking_rois)):
                    roi = RectangleBgROI()
                    roi.setLineWidth(1)
                    roi.setLineStyle("-")
                    roi.setColor("red")
                    roi.setBgStyle("pink", "-", 1.0)
                    roi.setVisible(False)
                    roi.setGeometry(origin=(0, 0), size=(0, 0))
                    self.rocking_rois.append(roi)
                    self.roiManager.addRoi(roi, useManagerColor=False)

            for roino, i in enumerate(np.arange(no_rois_to_display) * divider):
                self._set_center_roi_preview_color(self.rocking_rois[roino])
                ckey = roi_keys["center"][int(i)]
                leftkey = roi_keys["left"][int(i)]
                rightkey = roi_keys["right"][int(i)]
                topkey = roi_keys["top"][int(i)]
                bottomkey = roi_keys["bottom"][int(i)]

                origin, size, left, right, top, bottom = _display_roi_geometry(
                    ckey, leftkey, rightkey, topkey, bottomkey
                )
                self.rocking_rois[roino].setGeometry(
                    origin=origin,
                    size=size,
                    left=left,
                    right=right,
                    top=top,
                    bottom=bottom,
                )
                self.rocking_rois[roino].setVisible(True)
                self.rocking_rois[roino].setEditable(True)
            for roi in self.rocking_rois[no_rois_to_display:]:
                roi.setVisible(False)
                roi.setEditable(False)

        self.roiManager._roisUpdated()
        # self.centralPlot.removeMarker('main_croi_loc')

    def getStaticROIparams(self, xy, **kwargs):
        """Calculate hkl, detector angles, and pixel metadata for fixed ROIs.

        :param xy:
            ROI center coordinates in detector pixels.
        :returns:
            Array containing hkl in r.l.u., detector angles in rad, trajectory
            coordinate, pixel coordinates, and detector mask flag.
        :rtype: numpy.ndarray

        .. note::
           CLI-safe when scan and UB state are loaded.
        """
        if self.fscan is None:
            raise Exception("No scan loaded!")
        mu, om = self.getMuOm()
        # mu_cryst = HKLVlieg.crystalAngles_singleArray(mu, self.ubcalc.n)

        if "mask" in kwargs:
            mask = kwargs["mask"]
            xy = xy[mask]

        if len(np.asarray(om).shape) == 0:
            om = np.full(len(self.fscan), om)

        hkl_del_gam = np.empty((xy.shape[0], len(self.fscan), 6), dtype=np.float64)
        for i, xy_i in enumerate(xy):
            x = np.full(len(self.fscan), xy_i[0])
            y = np.full(len(self.fscan), xy_i[1])
            gamma, delta, alpha = self.ubcalc.detectorCal.crystalAnglesPoint(
                y, x, mu, self.ubcalc.n
            )

            if len(np.asarray(alpha).shape) == 0:
                alpha = np.full(len(self.fscan), alpha)

            hkl = self.ubcalc.angles.anglesToHkl(
                alpha, delta, gamma, om, self.ubcalc.chi, self.ubcalc.phi
            )
            # for i in range(len(self.fscan)):

            hkl_del_gam[i, :, :3] = np.array(hkl).T
            hkl_del_gam[i, :, 3] = delta
            hkl_del_gam[i, :, 4] = gamma
            hkl_del_gam[i, :, 5] = self.fscan.axis
        return hkl_del_gam

    def getROIloc(self, imageno=None, H_0=None, H_1=None, **kwargs):
        """Calculate detector ROI locations for a reciprocal-space line.

        When :math:`\\vec{H}_0` and :math:`\\vec{H}_1` are provided, they
        define the line :math:`\\vec{H}_0 + s\\vec{H}_1` in reciprocal space.
        The returned arrays describe the two possible detector intersections
        of that line with the Ewald sphere.

        :param int imageno:
            Optional image index. If omitted, locations are calculated over the
            scan where applicable.
        :param numpy.ndarray H_0:
            Starting reciprocal-space vector in r.l.u. Shape ``(3,)`` or
            ``(N, 3)``.
        :param numpy.ndarray H_1:
            Reciprocal-space direction vector in r.l.u. Shape ``(3,)`` or
            ``(N, 3)``.
        :returns:
            Two arrays for the two possible detector intersections.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]

        .. note::
           CLI-capable. Missing coordinates are read from GUI controls.
        """
        if self.fscan is None:
            raise Exception("No scan loaded!")

        mu, om = self.getMuOm(imageno)
        mu_cryst = HKLVlieg.crystalAngles_singleArray(mu, self.ubcalc.n)
        dc = self.ubcalc.detectorCal
        # mu = self.ubcalc.mu
        angles = self.ubcalc.angles

        if self.scanSelector.scanstab.currentIndex() == 1 and not kwargs.get(
            "intersect", False
        ):
            if imageno is None:
                hkl_del_gam_1 = np.ones((len(self.fscan), 6), dtype=np.float64)
                x = np.full(len(self.fscan), self.scanSelector.xy_static[0].value())
                y = np.full(len(self.fscan), self.scanSelector.xy_static[1].value())

                gamma, delta, alpha = self.ubcalc.detectorCal.crystalAnglesPoint(
                    y, x, mu, self.ubcalc.n
                )
                np.arange(len(self.fscan))

                if len(np.asarray(om).shape) == 0:
                    om = np.full(len(self.fscan), om)

                if len(np.asarray(alpha).shape) == 0:
                    alpha = np.full(len(self.fscan), alpha)

                yx1 = np.vstack((y, x)).T
                yx2 = np.full_like(yx1, np.inf)

                for i in range(len(self.fscan)):
                    pos = [
                        alpha[i],
                        delta[i],
                        gamma[i],
                        om[i],
                        self.ubcalc.chi,
                        self.ubcalc.phi,
                    ]
                    # pos = np.vstack(pos).T
                    hkl = np.array(self.ubcalc.angles.anglesToHkl(*pos))
                    hkl_del_gam_1[i, :3] = hkl
                hkl_del_gam_1[:, 3] = delta
                hkl_del_gam_1[:, 4] = gamma
                hkl_del_gam_1[:, 5] = self.fscan.axis
                hkl_del_gam_2 = np.full_like(hkl_del_gam_1, -1)
            else:
                hkl_del_gam_1 = np.ones(6, dtype=np.float64)
                x = self.scanSelector.xy_static[0].value()
                y = self.scanSelector.xy_static[1].value()
                yx1 = np.zeros((1, 2))
                yx1[0][0] = y
                yx1[0][1] = x
                yx2 = np.full_like(yx1, np.inf)

                if len(np.asarray(om).shape) > 0:
                    om = om[imageno]
                if len(np.asarray(mu).shape) > 0:
                    mu = mu[imageno]
                gamma, delta, alpha = self.ubcalc.detectorCal.crystalAnglesPoint(
                    np.array([y]), np.array([x]), mu, self.ubcalc.n
                )
                gamma, delta, alpha = (
                    gamma[0],
                    delta[0],
                    alpha,
                )  # crystalAnglesPoint retains shape, even for 0d array
                pos = [alpha, delta, gamma, om, self.ubcalc.chi, self.ubcalc.phi]
                hkl_del_gam_1[:3] = np.array(self.ubcalc.angles.anglesToHkl(*pos))
                hkl_del_gam_1[3] = delta
                hkl_del_gam_1[4] = gamma
                hkl_del_gam_1[5] = self.fscan.axis[imageno]
                hkl_del_gam_2 = np.full_like(hkl_del_gam_1, -1)

        else:
            if H_0 is None or H_1 is None:
                H_1 = np.array([h.value() for h in self.scanSelector.H_1])
                H_0 = np.array([h.value() for h in self.scanSelector.H_0])

            hkl_del_gam_1, hkl_del_gam_2, Qa_1, Qa_2 = angles.anglesIntersectLineEwald(
                H_0, H_1, mu_cryst, om, self.ubcalc.phi, self.ubcalc.chi, Qalpha=True
            )
            # H, K, L ,delta_1, gamma_1, HKL_Q1[-1]=s

            delta1 = hkl_del_gam_1[..., 3]
            delta2 = hkl_del_gam_2[..., 3]
            gam1 = hkl_del_gam_1[..., 4]
            gam2 = hkl_del_gam_2[..., 4]

            Qmin, Qmax = dc.Qrange
            Qa_1_n = np.linalg.norm(Qa_1, axis=-1)
            Qa_2_n = np.linalg.norm(Qa_2, axis=-1)

            mask1 = np.logical_and(Qmin <= Qa_1_n, Qmax >= Qa_1_n)
            mask2 = np.logical_and(Qmin <= Qa_2_n, Qmax >= Qa_2_n)

            yx1 = dc.pixelsCrystalAngles(gam1, delta1, mu, self.ubcalc.n)
            yx2 = dc.pixelsCrystalAngles(gam2, delta2, mu, self.ubcalc.n)
            yx1[~mask1] = np.inf
            yx2[~mask2] = np.inf

        ymask1 = np.logical_and(yx1[..., 0] >= 0, yx1[..., 0] < dc.detector.shape[0])
        xmask1 = np.logical_and(yx1[..., 1] >= 0, yx1[..., 1] < dc.detector.shape[1])
        yxmask1 = np.logical_and(xmask1, ymask1)

        ymask2 = np.logical_and(yx2[..., 0] >= 0, yx2[..., 0] < dc.detector.shape[0])
        xmask2 = np.logical_and(yx2[..., 1] >= 0, yx2[..., 1] < dc.detector.shape[1])
        yxmask2 = np.logical_and(xmask2, ymask2)

        xy1 = yx1[..., ::-1]
        xy2 = yx2[..., ::-1]

        if (
            not kwargs.get("intersect", False)
            and self.scanSelector.scanstab.currentIndex() != 1
        ):
            xoffset, yoffset = self.scanSelector.roioptions.get_offsets()

            if xoffset != 0.0 or yoffset != 0.0:
                logger.warn(
                    "Nonzero pixel offset selected. Experimental feature! Angles and hkl are incorrect!!!"  # noqa: E501
                )
                xy1[..., 0] += xoffset
                xy2[..., 0] += xoffset
                xy1[..., 1] += yoffset
                xy2[..., 1] += yoffset

        return np.concatenate(
            (np.atleast_2d(hkl_del_gam_1), xy1, yxmask1[..., np.newaxis]), axis=-1
        ), np.concatenate(
            (np.atleast_2d(hkl_del_gam_2), xy2, yxmask2[..., np.newaxis]), axis=-1
        )

    def plotROI(self, loc, roi):
        """Set a visible ROI object around a detector-pixel location.

        :param loc:
            ROI center in detector pixels.
        :param roi:
            ROI widget to update.

        .. note::
           CLI-capable, but it mutates Qt ROI widgets.
        """

        key = self.intkey(loc)
        leftkey, rightkey, topkey, bottomkey = self.bkgkeys(loc)

        # print([(roi, roi.isEditable()) for roi in self.rois])

        origin, size, left, right, top, bottom = _display_roi_geometry(
            key, leftkey, rightkey, topkey, bottomkey
        )
        self._set_center_roi_preview_color(roi)
        roi.setGeometry(
            origin=origin, size=size, left=left, right=right, top=top, bottom=bottom
        )
        roi.setVisible(True)
        # self.roiManager._roisUpdated()

    def integrateROI(self):
        """Integrate the active ROI workflow and save data to the database.

        The active tab in ``self.scanSelector.scanstab`` selects the
        integration workflow:

        * ``hklscan`` (tab id ``0``): integrate stationary-scan ROIs whose
          detector
          coordinates are calculated from the reciprocal-space line
          :math:`\\vec{H}_0 + s\\vec{H}_1`. :math:`\\vec{H}_0` and
          :math:`\\vec{H}_1` are numpy vector values in r.l.u. read from the
          ROI controls, and the two Ewald-sphere intersections are integrated
          as separate S1/S2 trajectories.
        * ``fixed`` (tab id ``1``): integrate a stationary detector-pixel ROI
          from the ``xy_static`` controls. The same pixel coordinates are used
          through the scan, while the corresponding reciprocal-space
          coordinates and diffractometer angles are recorded for each image.
        * ``rocking hklscan`` (tab id ``2``): delegate to
          :meth:`rocking_extraction`, which integrates multiple rocking-scan
          ROIs whose coordinates are sampled along
          :math:`\\vec{H}_0 + s\\vec{H}_1`.
        * ``rocking Bragg`` (tab id ``3``): delegate to
          :meth:`rocking_Bragg_extraction`, which calculates Bragg peak
          coordinates from the current crystal, detector, UB, strain, and scan
          state before integrating valid rocking-scan ROIs.

        All modes use the current UI/database state for scan selection, ROI
        sizes, masks, background settings, and correction factors. The
        resulting intensities and metadata are written to the active Nexus
        database file.

        :returns:
            Status dictionary describing success, cancellation, or error.
        :rtype: dict

        .. note::
           CLI-capable when scan, database, and ROI state are preconfigured.
        """

        if self.scanSelector.scanstab.currentIndex() == 2:
            return self.rocking_extraction()
        elif self.scanSelector.scanstab.currentIndex() == 3:
            return self.rocking_Bragg_extraction()

        try:
            image = self.fscan.get_raw_img(0)
        except Exception:
            logger.exception(
                "Cannot perform stationary scan integration: no images found.",
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Cannot perform stationary scan integration: no images found.",  # noqa: E501
                    "show_dialog": False,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            # print("no images found! %s" % e)
            return {
                "status": "error",
                "message": "No image found in current scan",
                "traceback": traceback.format_exc(),
            }
        if self.database.nxfile is None:
            logger.exception(
                "Cannot perform stationary scan integration: no database available.",
                extra={
                    "title": "Cannot integrate scan",
                    "description": "Cannot perform stationary scan integration: no database available.",  # noqa: E501
                    "show_dialog": False,
                    "dialog_level": logging.WARNING,
                    "parent": self,
                },
            )
            # print("No database available")
            return {"status": "error", "message": "No database available"}

        logger.info("Start integration of stationary scan")
        dc = self.ubcalc.detectorCal
        # mu = self.ubcalc.mu

        H_1 = np.array([h.value() for h in self.scanSelector.H_1])
        H_0 = np.array([h.value() for h in self.scanSelector.H_0])

        vsize = int(self.scanSelector.vsize.value())
        hsize = int(self.scanSelector.hsize.value())
        vsize * hsize  # as set in GUI, no corrections

        imgmask = None

        if self.scanSelector.useMaskBox.isChecked():
            if self.centralPlot.getMaskToolsDockWidget().getSelectionMask() is None:
                if logger_utils.get_logging_context() == "gui":
                    btn = qt.QMessageBox.question(
                        self,
                        "No mask available",
                        """No mask was selected with the masking tool.
        Do you want to continue without mask?""",
                    )
                    if btn != qt.QMessageBox.Yes:
                        return {
                            "status": "cancelled",
                            "message": "Reason: no mask selected",
                        }
                logger.warn(
                    "No mask was selected with the masking tool. Continue without mask."
                )
            else:
                imgmask = (
                    self.centralPlot.getMaskToolsDockWidget().getSelectionMask() > 0.0
                )

        corr = (
            self.scanSelector.useSolidAngleBox.isChecked()
            or self.scanSelector.usePolarizationBox.isChecked()
        )

        C_arr = np.ones(dc.detector.shape, dtype=np.float64)
        if self.scanSelector.useSolidAngleBox.isChecked():
            C_arr /= dc.solidAngleArray()
        if self.scanSelector.usePolarizationBox.isChecked():
            C_arr /= dc.polarization(factor=dc._polFactor, axis_offset=dc._polAxis)

        hkl_del_gam_s1, hkl_del_gam_s2 = self.getROIloc()

        nodatapoints = len(self.fscan)
        # print(hkl_del_gam_1s.shape)

        if hkl_del_gam_s1.shape[0] == 1:
            hkl_del_gam_1 = np.zeros(
                (nodatapoints, hkl_del_gam_s1.shape[1]), dtype=np.float64
            )
            hkl_del_gam_2 = np.zeros(
                (nodatapoints, hkl_del_gam_s1.shape[1]), dtype=np.float64
            )
            hkl_del_gam_1[:] = hkl_del_gam_s1[0]
            hkl_del_gam_2[:] = hkl_del_gam_s2[0]
        else:
            hkl_del_gam_1, hkl_del_gam_2 = hkl_del_gam_s1, hkl_del_gam_s2

        dataavail = np.logical_or(hkl_del_gam_1[:, -1], hkl_del_gam_2[:, -1])

        croi1_a = np.zeros_like(dataavail, dtype=np.float64)
        cpixel1_a = np.zeros_like(dataavail, dtype=np.float64)
        bgroi1_a = np.zeros_like(dataavail, dtype=np.float64)
        bgpixel1_a = np.zeros_like(dataavail, dtype=np.float64)
        x_coord1_a = hkl_del_gam_1[:, 6]
        y_coord1_a = hkl_del_gam_1[:, 7]
        roi_hsize1_a = np.full_like(dataavail, hsize, dtype=int)
        roi_vsize1_a = np.full_like(dataavail, vsize, dtype=int)

        croi2_a = np.zeros_like(dataavail, dtype=np.float64)
        cpixel2_a = np.zeros_like(dataavail, dtype=np.float64)
        bgroi2_a = np.zeros_like(dataavail, dtype=np.float64)
        bgpixel2_a = np.zeros_like(dataavail, dtype=np.float64)

        bgimg_croi1_a = np.zeros_like(dataavail, dtype=np.float64)
        bgimg_cpixel1_a = np.zeros_like(dataavail, dtype=np.float64)
        bgimg_bgroi1_a = np.zeros_like(dataavail, dtype=np.float64)
        bgimg_bgpixel1_a = np.zeros_like(dataavail, dtype=np.float64)

        Corr_croi1_a = np.zeros_like(dataavail, dtype=np.float64)
        Corr_cpixel1_a = np.zeros_like(dataavail, dtype=np.float64)
        Corr_bgroi1_a = np.zeros_like(dataavail, dtype=np.float64)
        Corr_bgpixel1_a = np.zeros_like(dataavail, dtype=np.float64)

        bgimg_croi2_a = np.zeros_like(dataavail, dtype=np.float64)
        bgimg_cpixel2_a = np.zeros_like(dataavail, dtype=np.float64)
        bgimg_bgroi2_a = np.zeros_like(dataavail, dtype=np.float64)
        bgimg_bgpixel2_a = np.zeros_like(dataavail, dtype=np.float64)

        Corr_croi2_a = np.zeros_like(dataavail, dtype=np.float64)
        Corr_cpixel2_a = np.zeros_like(dataavail, dtype=np.float64)
        Corr_bgroi2_a = np.zeros_like(dataavail, dtype=np.float64)
        Corr_bgpixel2_a = np.zeros_like(dataavail, dtype=np.float64)

        x_coord2_a = hkl_del_gam_2[:, 6]
        y_coord2_a = hkl_del_gam_2[:, 7]
        roi_hsize2_a = np.full_like(dataavail, hsize, dtype=int)
        roi_vsize2_a = np.full_like(dataavail, vsize, dtype=int)

        progress = logger_utils.create_progress_logger(
            self, len(self.fscan), "Integrating stationary scan"
        )

        has_bg_img = False
        roioptions = self.scanSelector.roioptions.get_parameters()
        use_fitted_background = bool(roioptions.get("FittedBackground", False))
        fitted_background_order = int(roioptions.get("FittedBackgroundOrder", 1))
        if use_fitted_background and not HAS_ACCEL:
            logger.warning(
                "Fitted local background requires the compiled ROI accelerator; "
                "using summed background ROIs instead."
            )

        if imgmask is not None:
            mask = np.ascontiguousarray(imgmask, dtype=bool)
        else:
            mask = np.zeros(image.img.shape, dtype=bool)
        if corr:
            C_arr = np.ascontiguousarray(C_arr, dtype=np.float64)
            C_arr[mask] = 0.0
        else:
            C_arr = np.ones(image.img.shape, dtype=np.float64)
            C_arr[mask] = 0.0

        for i in range(len(self.fscan)):
            key = self.intkey(hkl_del_gam_1[i, 6:8])
            croi_key = np.array(
                [[key[0].start, key[0].stop], [key[1].start, key[1].stop]]
            )
            roi_hsize1_a[i] = int(np.abs(np.diff(croi_key[0])[0]))
            roi_vsize1_a[i] = int(np.abs(np.diff(croi_key[1])[0]))
            key = self.intkey(hkl_del_gam_2[i, 6:8])
            croi_key = np.array(
                [[key[0].start, key[0].stop], [key[1].start, key[1].stop]]
            )
            roi_hsize2_a[i] = int(np.abs(np.diff(croi_key[0])[0]))
            roi_vsize2_a[i] = int(np.abs(np.diff(croi_key[1])[0]))

        if HAS_ACCEL:
            roi_lists_accel = []
            for i in range(len(self.fscan)):
                roi_lists = [[], [], [], [], []]
                if hkl_del_gam_1[i, -1]:
                    key = self.intkey(hkl_del_gam_1[i, 6:8])
                    croi_key = np.array(
                        [[key[0].start, key[0].stop], [key[1].start, key[1].stop]]
                    )
                    roi_lists[0].append(croi_key)  # center
                    bkgkey = self.bkgkeys(hkl_del_gam_1[i, 6:8])
                    for r, l in zip(bkgkey, roi_lists[1:]):  # noqa: E741
                        l.append(
                            np.array([[r[0].start, r[0].stop], [r[1].start, r[1].stop]])
                        )
                else:
                    [
                        l.append(np.array([[0, 0], [0, 0]])) for l in roi_lists[1:]  # noqa: E741
                    ]  # will result in zeros, convert to np.nan later
                    roi_lists[0].append(np.array([[0, 0], [0, 0]]))
                if hkl_del_gam_2[i, -1]:
                    key = self.intkey(hkl_del_gam_2[i, 6:8])
                    croi_key = np.array(
                        [[key[0].start, key[0].stop], [key[1].start, key[1].stop]]
                    )
                    roi_lists[0].append(croi_key)  # center
                    bkgkey = self.bkgkeys(hkl_del_gam_2[i, 6:8])
                    for r, l in zip(bkgkey, roi_lists[1:]):  # noqa: E741
                        l.append(
                            np.array([[r[0].start, r[0].stop], [r[1].start, r[1].stop]])
                        )
                else:
                    [
                        l.append(np.array([[0, 0], [0, 0]])) for l in roi_lists[1:]  # noqa: E741
                    ]  # will result in zeros, convert to np.nan later
                    roi_lists[0].append(np.array([[0, 0], [0, 0]]))
                roi_lists = [
                    np.ascontiguousarray(np.stack(l), dtype=np.int64) for l in roi_lists  # noqa: E741
                ]
                roi_lists_accel.append(roi_lists)

            if (
                self.background_image is not None
                and self.background_image.shape == image.img.shape
            ):
                if use_fitted_background:
                    logger.warning(
                        "Fitted local background is ignored when a background "
                        "image is selected."
                    )
                has_bg_img = True
                background_image = self.background_image.astype(
                    np.float64, order="C", copy=True
                )
                background_image[mask] = 0.0

                def sumImage(i):
                    """CLI-safe worker: integrate one stationary image with background."""  # noqa: E501
                    all_counters = np.zeros(
                        (roi_lists_accel[i][0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    Carr_counters = np.zeros(
                        (roi_lists_accel[i][0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    BgImg_counters = np.zeros(
                        (roi_lists_accel[i][0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    if not dataavail[i]:
                        return all_counters, Carr_counters, BgImg_counters
                    image = self.fscan.get_raw_img(i).img.astype(
                        np.float64, order="C", copy=True
                    )  # unlocks gil during file read
                    _roi_sum_accel.processImage_bg_Carr(
                        image,
                        background_image,
                        mask,
                        C_arr,
                        *roi_lists_accel[i],
                        all_counters,
                        Carr_counters,
                        BgImg_counters,
                    )  # compiled accelerator releases the GIL
                    return all_counters, Carr_counters, BgImg_counters
            else:

                def sumImage(i):
                    """CLI-safe worker: integrate one stationary image."""
                    all_counters = np.zeros(
                        (roi_lists_accel[i][0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    Carr_counters = np.zeros(
                        (roi_lists_accel[i][0].shape[0],) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    if not dataavail[i]:
                        return all_counters, Carr_counters
                    image = self.fscan.get_raw_img(i).img.astype(
                        np.float64, order="C", copy=True
                    )  # unlocks gil during file read
                    if use_fitted_background:
                        _roi_sum_accel.processImage_polybg_Carr(
                            image,
                            mask,
                            C_arr,
                            *roi_lists_accel[i],
                            all_counters,
                            Carr_counters,
                            fitted_background_order,
                        )  # compiled accelerator releases the GIL
                    else:
                        _roi_sum_accel.processImage_Carr(
                            image,
                            mask,
                            C_arr,
                            *roi_lists_accel[i],
                            all_counters,
                            Carr_counters,
                        )  # compiled accelerator releases the GIL
                    return all_counters, Carr_counters

        else:  # not HAS_ACCEL
            if (
                self.background_image is not None
                and self.background_image.shape == image.img.shape
            ):
                has_bg_img = True
                background_image = self.background_image.astype(
                    np.float64, order="C", copy=True
                )
                background_image[mask] = 0.0

                def sumImage(i):
                    """CLI-safe worker: integrate one image with background."""
                    all_counters = np.zeros(
                        (2,) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    Carr_counters = np.zeros(
                        (2,) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    BgImg_counters = np.zeros(
                        (2,) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    if not dataavail[i]:
                        return all_counters, Carr_counters, BgImg_counters
                    else:
                        image = self.fscan.get_raw_img(i).img.astype(
                            np.float64, order="C", copy=True
                        )
                        if imgmask is not None:
                            image[imgmask] = np.nan
                            pixelavail = (~imgmask).astype(np.float64)
                        else:
                            pixelavail = np.ones_like(image)

                        for intersect, hkl_del_gam_current in zip(
                            range(2), [hkl_del_gam_1, hkl_del_gam_2]
                        ):
                            if hkl_del_gam_current[i, -1]:
                                key = self.intkey(hkl_del_gam_current[i, 6:8])
                                bkgkey = self.bkgkeys(hkl_del_gam_current[i, 6:8])

                                all_counters[intersect, 0] = np.nansum(image[key[::-1]])
                                Carr_counters[intersect, 0] = np.nansum(
                                    C_arr[key[::-1]]
                                )
                                BgImg_counters[intersect, 0] = np.nansum(
                                    background_image[key[::-1]]
                                )

                                cpixel1 = np.nansum(pixelavail[key[::-1]])
                                all_counters[intersect, 1] = cpixel1
                                Carr_counters[intersect, 1] = cpixel1
                                BgImg_counters[intersect, 1] = cpixel1

                                bgpixel1 = 0.0
                                for bg in bkgkey:
                                    image[bg[::-1]]
                                    all_counters[intersect, 2] += np.nansum(
                                        image[bg[::-1]]
                                    )
                                    Carr_counters[intersect, 2] += np.nansum(
                                        C_arr[bg[::-1]]
                                    )
                                    BgImg_counters[intersect, 2] += np.nansum(
                                        background_image[bg[::-1]]
                                    )
                                    bgpixel1 += np.nansum(pixelavail[bg[::-1]])

                                all_counters[intersect, 3] = bgpixel1
                                Carr_counters[intersect, 3] = bgpixel1
                                BgImg_counters[intersect, 3] = bgpixel1
                        return all_counters, Carr_counters, BgImg_counters
            else:

                def sumImage(i):
                    """CLI-safe worker: integrate one image without acceleration."""
                    all_counters = np.zeros(
                        (2,) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    Carr_counters = np.zeros(
                        (2,) + (4,), dtype=np.float64
                    )  # need gil for python object creation
                    if not dataavail[i]:
                        return all_counters, Carr_counters
                    else:
                        image = self.fscan.get_raw_img(i).img.astype(
                            np.float64, order="C", copy=True
                        )
                        if imgmask is not None:
                            image[imgmask] = np.nan
                            pixelavail = (~imgmask).astype(np.float64)
                        else:
                            pixelavail = np.ones_like(image)

                        for intersect, hkl_del_gam_current in zip(
                            range(2), [hkl_del_gam_1, hkl_del_gam_2]
                        ):
                            if hkl_del_gam_current[i, -1]:
                                key = self.intkey(hkl_del_gam_current[i, 6:8])
                                bkgkey = self.bkgkeys(hkl_del_gam_current[i, 6:8])

                                all_counters[intersect, 0] = np.nansum(image[key[::-1]])
                                Carr_counters[intersect, 0] = np.nansum(
                                    C_arr[key[::-1]]
                                )

                                cpixel1 = np.nansum(pixelavail[key[::-1]])
                                all_counters[intersect, 1] = cpixel1
                                Carr_counters[intersect, 1] = cpixel1

                                bgpixel1 = 0.0
                                for bg in bkgkey:
                                    image[bg[::-1]]
                                    all_counters[intersect, 2] += np.nansum(
                                        image[bg[::-1]]
                                    )
                                    Carr_counters[intersect, 2] += np.nansum(
                                        C_arr[bg[::-1]]
                                    )
                                    bgpixel1 += np.nansum(pixelavail[bg[::-1]])

                                all_counters[intersect, 3] = bgpixel1
                                Carr_counters[intersect, 3] = bgpixel1
                        return all_counters, Carr_counters

        cancelled = False
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.numberthreads
        ) as executor:  # speedup only for the file reads
            futures = {}
            for i in range(len(self.fscan)):
                futures[executor.submit(sumImage, i)] = i

            for f in concurrent.futures.as_completed(futures):
                try:
                    i = futures[f]
                    # (croi1, cpixel1, bgroi1, bgpixel1), (croi2, cpixel2, bgroi2, bgpixel2) = f.result()  # noqa: E501
                    if has_bg_img:
                        all_counters, Carr_counters, BgImg_counters = f.result()
                        bgimg_croi1_a[i] = BgImg_counters[0, 0]
                        bgimg_cpixel1_a[i] = BgImg_counters[0, 1]
                        bgimg_bgroi1_a[i] = BgImg_counters[0, 2]
                        bgimg_bgpixel1_a[i] = BgImg_counters[0, 3]
                        bgimg_croi2_a[i] = BgImg_counters[1, 0]
                        bgimg_cpixel2_a[i] = BgImg_counters[1, 1]
                        bgimg_bgroi2_a[i] = BgImg_counters[1, 2]
                        bgimg_bgpixel2_a[i] = BgImg_counters[1, 3]

                    else:
                        all_counters, Carr_counters = f.result()
                        bgimg_croi1_a[i] = 0.0
                        bgimg_cpixel1_a[i] = 0.0
                        bgimg_bgroi1_a[i] = 0.0
                        bgimg_bgpixel1_a[i] = 0.0
                        bgimg_croi2_a[i] = 0.0
                        bgimg_cpixel2_a[i] = 0.0
                        bgimg_bgroi2_a[i] = 0.0
                        bgimg_bgpixel2_a[i] = 0.0

                    croi1_a[i] = all_counters[0, 0]
                    cpixel1_a[i] = all_counters[0, 1]
                    bgroi1_a[i] = all_counters[0, 2]
                    bgpixel1_a[i] = all_counters[0, 3]
                    croi2_a[i] = all_counters[1, 0]
                    cpixel2_a[i] = all_counters[1, 1]
                    bgroi2_a[i] = all_counters[1, 2]
                    bgpixel2_a[i] = all_counters[1, 3]

                    Corr_croi1_a[i] = Carr_counters[0, 0]
                    Corr_cpixel1_a[i] = Carr_counters[0, 1]
                    Corr_bgroi1_a[i] = Carr_counters[0, 2]
                    Corr_bgpixel1_a[i] = Carr_counters[0, 3]
                    Corr_croi2_a[i] = Carr_counters[1, 0]
                    Corr_cpixel2_a[i] = Carr_counters[1, 1]
                    Corr_bgroi2_a[i] = Carr_counters[1, 2]
                    Corr_bgpixel2_a[i] = Carr_counters[1, 3]

                    progress.update(futures[f])
                except concurrent.futures.CancelledError:
                    pass
                except Exception:
                    logger.warn(f"Cannot read image:\n{traceback.format_exc()}")

                if progress.wasCanceled():
                    [f.cancel() for f in futures]
                    cancelled = True
                    break

        progress.finish()

        if cancelled:
            return {
                "status": "cancelled",
                "message": "Reason: Cancelled during integration",
            }

        roi_size1 = roi_hsize1_a * roi_vsize1_a
        roi_size2 = roi_hsize2_a * roi_vsize2_a

        Corr1 = Corr_croi1_a * (
            roi_size1 / Corr_cpixel1_a
        )  # normalize to number of pixels of center roi (croi)
        Corr2 = Corr_croi2_a * (roi_size2 / Corr_cpixel2_a)
        croibg1_bgimg_a = None
        croibg1_bgimg_err_a = None

        if np.any(
            bgimg_cpixel1_a
        ):  # assume the background image has no errors (would need a separate error image for that)  # noqa: E501
            bgimg_croi1_norm = bgimg_croi1_a * (cpixel1_a / bgimg_cpixel1_a)
            if np.any(bgpixel1_a):
                bgimg_bgroi1_norm = bgimg_bgroi1_a * (bgpixel1_a / bgimg_bgpixel1_a)

                # method 1: simply subtract bg image from data and then subtract the remaining background  # noqa: E501
                croibg1_a = (
                    (croi1_a - bgimg_croi1_norm)
                    - (cpixel1_a / bgpixel1_a) * (bgroi1_a - bgimg_bgroi1_norm)
                ) * (roi_size1 / cpixel1_a)
                croibg1_err_a = np.sqrt(
                    croi1_a + ((cpixel1_a / bgpixel1_a) ** 2) * bgroi1_a
                ) * (roi_size1 / cpixel1_a)

                # method 2: scale bg image croi and subtract scaled bg image croi. Use ratio of bgroi of image and bg image as scale factor.  # noqa: E501
                factor = bgroi1_a / bgimg_bgroi1_norm
                croibg1_bgimg_a = (croi1_a - factor * bgimg_croi1_norm) * (
                    roi_size1 / cpixel1_a
                )
                croibg1_bgimg_err_a = np.sqrt(
                    croi1_a + ((cpixel1_a / bgpixel1_a) ** 2) * bgroi1_a
                ) * (roi_size1 / cpixel1_a)

            else:  # not possible if no bgroi is set.
                croibg1_a = (croi1_a - bgimg_croi1_norm) * (roi_size1 / cpixel1_a)
                croibg1_err_a = np.sqrt(croi1_a) * (roi_size1 / cpixel1_a)

        else:  # no background image
            if np.any(bgpixel1_a):
                croibg1_a = (croi1_a - (cpixel1_a / bgpixel1_a) * bgroi1_a) * (
                    roi_size1 / cpixel1_a
                )
                croibg1_err_a = np.sqrt(
                    croi1_a + ((cpixel1_a / bgpixel1_a) ** 2) * bgroi1_a
                ) * (roi_size1 / cpixel1_a)
            else:
                croibg1_a = croi1_a * (roi_size1 / cpixel1_a)
                croibg1_err_a = np.sqrt(croi1_a) * (roi_size1 / cpixel1_a)

        croibg2_bgimg_a = None
        croibg2_bgimg_err_a = None
        if np.any(
            bgimg_cpixel2_a
        ):  # assume the background image has no errors (would need a separate error image for that)  # noqa: E501
            bgimg_croi2_norm = bgimg_croi2_a * (cpixel2_a / bgimg_cpixel2_a)
            if np.any(bgpixel2_a):
                bgimg_bgroi2_norm = bgimg_bgroi2_a * (bgpixel2_a / bgimg_bgpixel2_a)

                # method 1: simply subtract bg image from data and then subtract the remaining background  # noqa: E501
                croibg2_a = (
                    (croi2_a - bgimg_croi2_norm)
                    - (cpixel2_a / bgpixel2_a) * (bgroi2_a - bgimg_bgroi2_norm)
                ) * (roi_size2 / cpixel2_a)
                croibg2_err_a = np.sqrt(
                    croi2_a + ((cpixel2_a / bgpixel2_a) ** 2) * bgroi2_a
                ) * (roi_size2 / cpixel2_a)

                # method 2: scale bg image croi and subtract scaled bg image croi. Use ratio of bgroi of image and bg image as scale factor.  # noqa: E501
                factor = bgroi2_a / bgimg_bgroi2_norm
                croibg2_bgimg_a = (croi2_a - factor * bgimg_croi2_norm) * (
                    roi_size2 / cpixel2_a
                )
                croibg2_bgimg_err_a = np.sqrt(
                    croi2_a + ((cpixel2_a / bgpixel2_a) ** 2) * bgroi2_a
                ) * (roi_size2 / cpixel2_a)

            else:  # not possible if no bgroi is set.
                croibg2_a = (croi2_a - bgimg_croi2_norm) * (roi_size2 / cpixel2_a)
                croibg2_err_a = np.sqrt(croi2_a) * (roi_size2 / cpixel2_a)

        else:  # no background image
            if np.any(bgpixel2_a):
                croibg2_a = (croi2_a - (cpixel2_a / bgpixel2_a) * bgroi2_a) * (
                    roi_size2 / cpixel2_a
                )
                croibg2_err_a = np.sqrt(
                    croi2_a + ((cpixel2_a / bgpixel2_a) ** 2) * bgroi2_a
                ) * (roi_size2 / cpixel2_a)
            else:
                croibg2_a = croi2_a * (roi_size2 / cpixel2_a)
                croibg2_err_a = np.sqrt(croi2_a) * (roi_size2 / cpixel2_a)

        if corr:
            croibg1_a *= Corr1
            croibg1_err_a *= Corr1
            croibg2_a *= Corr2
            croibg2_err_a *= Corr2
            if croibg1_bgimg_a is not None:
                croibg1_bgimg_a *= Corr1
                croibg1_bgimg_err_a *= Corr1
            if croibg2_bgimg_a is not None:
                croibg2_bgimg_a *= Corr2
                croibg2_bgimg_err_a *= Corr2

        rod_mask1 = np.isfinite(croibg1_a)
        rod_mask2 = np.isfinite(croibg2_a)

        s1_masked = hkl_del_gam_1[:, 5][rod_mask1]
        s2_masked = hkl_del_gam_2[:, 5][rod_mask2]

        croibg1_a_masked = croibg1_a[rod_mask1]
        croibg2_a_masked = croibg2_a[rod_mask2]

        croibg1_err_a_masked = croibg1_err_a[rod_mask1]
        croibg2_err_a_masked = croibg2_err_a[rod_mask2]

        # name = str(H_1) + "*s+" + str(H_0)
        if self.scanSelector.scanstab.currentIndex() == 1:
            x = self.scanSelector.xy_static[0].value()
            y = self.scanSelector.xy_static[1].value()
            name1 = f"pixloc[{x:.2f} {y:.2f}]"
            name2 = (
                f"pixloc[{x:.2f} {y:.2f}]_2"  # does not exist, Just for compatibility
            )
            traj1 = {
                "@NX_class": "NXcollection",
                "@direction": "Fixed pixel coordinates",
                "s": hkl_del_gam_1[:, 5],
            }
            traj2 = {
                "@NX_class": "NXcollection",
                "@direction": "Fixed pixel coordinates",
                "s": hkl_del_gam_2[:, 5],
            }
        else:
            name1 = str(H_1) + "*s1+" + str(H_0)
            name2 = str(H_1) + "*s2+" + str(H_0)
            traj1 = {
                "@NX_class": "NXcollection",
                "@direction": "Intergrated along H_1*s + H_0 in reciprocal space",
                "H_1": H_1,
                "H_0": H_0,
                "s": hkl_del_gam_1[:, 5],
            }
            traj2 = {
                "@NX_class": "NXcollection",
                "@direction": "Intergrated along H_1*s + H_0 in reciprocal space",
                "H_1": H_1,
                "H_0": H_0,
                "s": hkl_del_gam_2[:, 5],
            }

        defaultS1 = croibg1_a_masked.size > croibg2_a_masked.size

        if hasattr(self.fscan, "title"):
            title = str(self.fscan.title)
        else:
            title = f"{self.fscan.axisname}-scan"

        mu, om = self.getMuOm()
        if len(np.asarray(om).shape) == 0:
            om = np.full_like(mu, om)
        if len(np.asarray(mu).shape) == 0:
            mu = np.full_like(om, mu)

        suffix = ""
        i = 0

        while (
            self.activescanname + "/measurement/" + name1 + suffix
            in self.database.nxfile
        ):
            suffix = f"_{i}"
            i += 1
        availname1 = name1 + suffix

        suffix = ""
        i = 0
        while (
            self.activescanname + "/measurement/" + name2 + suffix
            in self.database.nxfile
        ):
            suffix = f"_{i}"
            i += 1

        availname2 = name2 + suffix

        auxcounters = {"@NX_class": "NXcollection"}
        for auxname in self.fscan.auxillary_counters:
            if hasattr(self.fscan, auxname):
                cntr = getattr(self.fscan, auxname)
                if cntr is not None:
                    auxcounters[auxname] = cntr

        datas1 = {
            "@NX_class": "NXdata",
            "sixc_angles": {
                "@NX_class": "NXpositioner",
                "alpha": np.rad2deg(mu),
                "omega": np.rad2deg(om),
                "theta": np.rad2deg(-1 * om),
                "delta": np.rad2deg(hkl_del_gam_1[:, 3]),
                "gamma": np.rad2deg(hkl_del_gam_1[:, 4]),
                "chi": np.rad2deg(self.ubcalc.chi),
                "phi": np.rad2deg(self.ubcalc.phi),
                "@unit": "deg",
            },
            "hkl": {
                "@NX_class": "NXcollection",
                "h": hkl_del_gam_1[:, 0],
                "k": hkl_del_gam_1[:, 1],
                "l": hkl_del_gam_1[:, 2],
            },
            "counters": {
                "@NX_class": "NXdetector",
                "croibg": croibg1_a,
                "croibg_errors": croibg1_err_a,
                "croibg_bgimg": croibg1_bgimg_a,  # when None, will not create data set
                "croibg_bgimg_errors": croibg1_bgimg_err_a,  # when None, will not create data set  # noqa: E501
                "croi": croi1_a,
                "bgroi": bgroi1_a,
                "croi_pix": cpixel1_a,
                "bgroi_pix": bgpixel1_a,
                "Cfactors_croi": Corr_croi1_a,
                "Cfactors_bgroi": Corr_bgroi1_a,
                "bgimg_croi": bgimg_croi1_a,
                "bgimg_bgroi": bgimg_bgroi1_a,
            },
            "pixelcoord": {
                "@NX_class": "NXdetector",
                "x": x_coord1_a,
                "y": y_coord1_a,
                "vsize": vsize,
                "hsize": hsize,
                "vsize_corr": roi_vsize1_a,
                "hsize_corr": roi_hsize1_a,
            },
            "trajectory": traj1,
            "@signal": "counters/croibg",
            "@axes": "trajectory/s",
            "@title": self.activescanname + "_" + availname1,
            "@orgui_meta": "roi",
        }

        datas2 = {
            "@NX_class": "NXdata",
            "sixc_angles": {
                "@NX_class": "NXpositioner",
                "alpha": np.rad2deg(mu),
                "omega": np.rad2deg(om),
                "theta": np.rad2deg(-1 * om),
                "delta": np.rad2deg(hkl_del_gam_2[:, 3]),
                "gamma": np.rad2deg(hkl_del_gam_2[:, 4]),
                "chi": np.rad2deg(self.ubcalc.chi),
                "phi": np.rad2deg(self.ubcalc.phi),
                "@unit": "deg",
            },
            "hkl": {
                "@NX_class": "NXcollection",
                "h": hkl_del_gam_2[:, 0],
                "k": hkl_del_gam_2[:, 1],
                "l": hkl_del_gam_2[:, 2],
            },
            "counters": {
                "@NX_class": "NXdetector",
                "croibg": croibg2_a,
                "croibg_errors": croibg2_err_a,
                "croibg_bgimg": croibg2_bgimg_a,
                "croibg_bgimg_errors": croibg2_bgimg_err_a,
                "croi": croi2_a,
                "bgroi": bgroi2_a,
                "croi_pix": cpixel2_a,
                "bgroi_pix": bgpixel2_a,
                "Cfactors_croi": Corr_croi2_a,
                "Cfactors_bgroi": Corr_bgroi2_a,
                "bgimg_croi": bgimg_croi2_a,
                "bgimg_bgroi": bgimg_bgroi2_a,
            },
            "pixelcoord": {
                "@NX_class": "NXdetector",
                "x": x_coord2_a,
                "y": y_coord2_a,
                "vsize": vsize,
                "hsize": hsize,
                "vsize_corr": roi_vsize2_a,
                "hsize_corr": roi_hsize2_a,
            },
            "trajectory": traj2,
            "@signal": "counters/croibg",
            "@axes": "trajectory/s",
            "@title": self.activescanname + "_" + availname2,
            "@orgui_meta": "roi",
        }

        data = {
            self.activescanname: {
                "instrument": {
                    "@NX_class": "NXinstrument",
                    "positioners": {
                        "@NX_class": "NXcollection",
                        self.fscan.axisname: self.fscan.axis,
                    },
                },
                "auxillary": auxcounters,
                "measurement": {
                    "@NX_class": "NXentry",
                    "@default": availname1 if defaultS1 else availname2,
                },
                "title": f"{title}",
                "@NX_class": "NXentry",
                "@default": "measurement/%s"
                % (availname1 if defaultS1 else availname2),
                "@orgui_meta": "scan",
            }
        }

        names_to_log = ""

        if np.any(cpixel1_a > 0.0):
            self.integrdataPlot.addCurve(
                s1_masked,
                croibg1_a_masked,
                legend=self.activescanname + "_" + availname1,
                xlabel="trajectory/s",
                ylabel="counters/croibg",
                yerror=croibg1_err_a_masked,
            )

            data[self.activescanname]["measurement"][availname1] = datas1
            names_to_log += availname1
        if np.any(cpixel2_a > 0.0):
            self.integrdataPlot.addCurve(
                s2_masked,
                croibg2_a_masked,
                legend=self.activescanname + "_" + availname2,
                xlabel="trajectory/s",
                ylabel="counters/croibg",
                yerror=croibg2_err_a_masked,
            )

            data[self.activescanname]["measurement"][availname2] = datas2

            names_to_log += availname2

        self.database.add_nxdict(data)
        logger.info(f"stationary scan integrated and saved with name(s) {names_to_log}")
        return {"status": "success"}

    def _graphCallback(self, eventdict):
        """GUI-only: handle plot mouse and marker events."""
        # print(eventdict)
        if eventdict["event"] == "mouseDoubleClicked":
            # newReflection = np.array([1,1,1,self.imageno,eventdict['x'],eventdict['y']])  # noqa: E501
            if self.scanSelector.select_roi_action.isChecked():
                self.scanSelector.set_xy_static_loc(eventdict["x"], eventdict["y"])
                self.scanSelector.select_roi_action.setChecked(False)
            else:
                hkl = self.centralPlot.xyHKLConverter(eventdict["x"], eventdict["y"])[
                    :3
                ]
                self.reflectionSel.addReflection(eventdict, self.imageno, hkl)

        if eventdict["event"] == "markerMoved":
            if eventdict["label"].startswith("__"):
                return
            self.reflectionSel.moveReflection(
                eventdict["label"], [eventdict["x"], eventdict["y"]]
            )
            self.reflectionSel.setReflectionActive(eventdict["label"])
        if eventdict["event"] == "markerClicked":
            if eventdict["label"].startswith("__"):
                return
            self.reflectionSel.setReflectionActive(eventdict["label"])

    def intkey(self, coords):
        """Create a center ROI key from detector pixel coordinates.

        The returned key represents the clipped horizontal and vertical bounds
        of the center ROI around ``coords``.

        ROI sizes are read from the current ``hsize`` and ``vsize`` controls.
        In non-fixed scan modes, advanced ROI options may replace those nominal
        sizes with detector-inclination or projected-sample-size corrected
        bounds before clipping to the detector.

        :param numpy.ndarray coords:
            ROI center coordinates ``(x, y)`` in detector pixels.
        :returns:
            ROI key ``(x_bounds, y_bounds)`` clipped to the detector extent.
        :rtype: tuple[slice, slice]

        .. note::
           CLI-capable. ROI dimensions are read from GUI controls.
        """

        vsize = int(self.scanSelector.vsize.value())
        hsize = int(self.scanSelector.hsize.value())

        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape

        coord_restr = np.clip(np.asarray(coords), [0, 0], [dethsize, detvsize])

        roioptions = self.scanSelector.roioptions.get_parameters()
        current_mode = self.scanSelector.scanstab.currentIndex()
        if (
            roioptions["DetectorInclination"] or roioptions["ProjectSampleSize"]
        ) and current_mode != 1:
            if roioptions["ProjectSampleSize"]:
                size_exact = ROIutils.calc_corrections(
                    coord_restr,
                    self.ubcalc.detectorCal,
                    np.array([hsize, vsize]),
                    roioptions,
                    roioptions["DetectorInclination"],
                    roioptions["factor"],
                )
            else:
                size_exact = ROIutils.calc_corrections(
                    coord_restr,
                    self.ubcalc.detectorCal,
                    np.array([hsize, vsize]),
                    None,
                    roioptions["DetectorInclination"],
                    roioptions["factor"],
                )
            hsize = size_exact[0][0]
            vsize = size_exact[0][1]

        vhalfsize = vsize // 2
        hhalfsize = hsize // 2
        fromcoords = np.round(
            np.asarray(coord_restr) - np.array([hhalfsize, vhalfsize])
        )
        tocoords = np.round(np.asarray(coord_restr) + np.array([hhalfsize, vhalfsize]))

        if hsize % 2:
            if coord_restr[0] % 1 < 0.5:
                tocoords[0] += 1
            else:
                fromcoords[0] -= 1
        if vsize % 2:
            if coord_restr[1] % 1 < 0.5:
                tocoords[1] += 1
            else:
                fromcoords[1] -= 1

        fromcoords = np.clip(np.asarray(fromcoords), [0, 0], [dethsize, detvsize])
        tocoords = np.clip(np.asarray(tocoords), [0, 0], [dethsize, detvsize])

        loc = tuple(
            slice(int(fromcoord), int(tocoord))
            for fromcoord, tocoord in zip(fromcoords, tocoords)
        )

        # from IPython import embed; embed()

        return loc

    def bkgkeys(self, coords):
        """Create background ROI keys around a center ROI.

        Background keys are derived from the center ROI key returned by
        :meth:`intkey`. The left and right background ROIs extend horizontally
        beside the center ROI and keep the same vertical bounds. The top and
        bottom background ROIs extend vertically above and below the center ROI
        and keep the same horizontal bounds.

        Background widths and heights are read from the current ``left``,
        ``right``, ``top``, and ``bottom`` controls. All bounds are clipped to
        the detector extent.

        :param numpy.ndarray coords:
            Center ROI coordinates ``(x, y)`` in detector pixels.
        :returns:
            Left, right, top, and bottom background ROI keys.
        :rtype: tuple[tuple[slice, slice], tuple[slice, slice], tuple[slice, slice], tuple[slice, slice]]

        .. note::
           CLI-capable. Background sizes are read from GUI controls.
        """  # noqa: E501

        left = int(self.scanSelector.left.value())
        right = int(self.scanSelector.right.value())
        top = int(self.scanSelector.top.value())
        bottom = int(self.scanSelector.bottom.value())

        detvsize, dethsize = self.ubcalc.detectorCal.detector.shape

        croi = self.intkey(coords)
        croi[0]
        croi[1]

        leftkey = (
            slice(int(np.clip(croi[0].start - left, 0, dethsize)), croi[0].start),
            croi[1],
        )
        rightkey = (
            slice(croi[0].stop, int(np.clip(croi[0].stop + right, 0, dethsize))),
            croi[1],
        )

        topkey = (
            croi[0],
            slice(int(np.clip(croi[1].start - top, 0, detvsize)), croi[1].start),
        )
        bottomkey = (
            croi[0],
            slice(croi[1].stop, int(np.clip(croi[1].stop + bottom, 0, detvsize))),
        )
        return leftkey, rightkey, topkey, bottomkey

    def _onCenterGraph(self, xy):
        """GUI/CLI hint: recenter the plot axes on detector pixel coordinates."""
        # img = self.centralPlot.getImage()
        # shape = img.shape if img is not None else (100,100)

        x1, x2 = self.centralPlot.getXAxis().getLimits()
        x_center = (x1 + x2) / 2
        x1_new = x1 - x_center + xy[0]
        x2_new = x2 - x_center + xy[0]
        self.centralPlot.getXAxis().setLimits(x1_new, x2_new)

        y1, y2 = self.centralPlot.getYAxis().getLimits()
        y_center = (y1 + y2) / 2
        y1_new = y1 - y_center + xy[1]
        y2_new = y2 - y_center + xy[1]
        self.centralPlot.getYAxis().setLimits(y1_new, y2_new)

    def _onROITrackingChanged(self, *args):
        """GUI-only: recenter the detector view after ROI tracking changes."""
        self._centerTrackedROI(self.imageno)

    def _centerTrackedROI(self, image_no):
        """Center the detector view on the selected stationary ROI intersect."""
        if not self.scanSelector.roiTrackingAct.isChecked():
            return
        if self.scanSelector.scanstab.currentIndex() not in [0, 1]:
            return
        try:
            hkl_del_gam_1, hkl_del_gam_2 = self.getROIloc(image_no)
            if self.scanSelector.roiTrackingS2Act.isChecked():
                hkl_del_gam = hkl_del_gam_2
            else:
                hkl_del_gam = hkl_del_gam_1
            if hkl_del_gam[0, -1]:
                self._onCenterGraph(hkl_del_gam[0, 6:8])
        except Exception:
            logger.warning(
                "Cannot center tracked ROI",
                exc_info=True,
            )

    def closeEvent(self, event):
        """GUI/CLI hint: close the database before the main window closes."""
        self.database.close()
        super().closeEvent(event)


class Plot2DHKL(silx.gui.plot.PlotWindow):
    sigKeyPressDelete = qt.pyqtSignal()

    def __init__(self, xyHKLConverter, parent=None, backend=None):
        """GUI-only: initialize the image plot with hkl position readout."""
        self.xyHKLConverter = xyHKLConverter

        posInfo = [
            ("X", lambda x, y: x),
            ("Y", lambda x, y: y),
            ("H", lambda x, y: self.xyHKLConverter(x, y)[0]),
            ("K", lambda x, y: self.xyHKLConverter(x, y)[1]),
            ("L", lambda x, y: self.xyHKLConverter(x, y)[2]),
            ("del", lambda x, y: self.xyHKLConverter(x, y)[3]),
            ("gam", lambda x, y: self.xyHKLConverter(x, y)[4]),
            ("Data", WeakMethodProxy(self._getImageValue)),
        ]

        super().__init__(
            parent=parent,
            backend=backend,
            resetzoom=True,
            autoScale=False,
            logScale=False,
            grid=False,
            curveStyle=False,
            colormap=True,
            aspectRatio=True,
            yInverted=True,
            copy=True,
            save=True,
            print_=True,
            control=True,
            position=posInfo,
            roi=False,
            mask=True,
        )

        if parent is None:
            self.setWindowTitle("Plot2D")
        self.getXAxis().setLabel("Columns")
        self.getYAxis().setLabel("Rows")

        # if silx.config.DEFAULT_PLOT_IMAGE_Y_AXIS_ORIENTATION == 'downward':
        self.getYAxis().setInverted(True)

        self.profile = ProfileToolBar(plot=self)
        self.addToolBar(self.profile)

        self.colorbarAction.setVisible(True)
        self.getColorBarWidget().setVisible(True)

        # Put colorbar action after colormap action
        actions = self.toolBar().actions()
        for action in actions:
            if action is self.getColormapAction():
                break

        self.sigActiveImageChanged.connect(self.__activeImageChanged)

    def keyPressEvent(self, event):
        """GUI-only: emit delete-key signal for plot interactions."""
        key = event.key()
        if key == qt.Qt.Key_Delete and not event.isAutoRepeat():
            self.sigKeyPressDelete.emit()
        super().keyPressEvent(event)

    def setXyHKLconverter(self, xyHKLConverter):
        """Set the pixel-to-hkl converter used by the plot readout.

        .. note::
           CLI-capable, but normally used by GUI plot setup.
        """
        self.xyHKLConverter = xyHKLConverter

    def __activeImageChanged(self, previous, legend):
        """Handle change of active image

        :param Union[str,None] previous: Legend of previous active image
        :param Union[str,None] legend: Legend of current active image
        """
        if previous is not None:
            item = self.getImage(previous)
            if item is not None:
                item.sigItemChanged.disconnect(self.__imageChanged)

        if legend is not None:
            item = self.getImage(legend)
            item.sigItemChanged.connect(self.__imageChanged)

        positionInfo = self.getPositionInfoWidget()
        if positionInfo is not None:
            positionInfo.updateInfo()

    def __imageChanged(self, event):
        """Handle update of active image item

        :param event: Type of changed event
        """
        if event == items.ItemChangedType.DATA:
            positionInfo = self.getPositionInfoWidget()
            if positionInfo is not None:
                positionInfo.updateInfo()

    def _getImageValue(self, x, y):
        """Get status bar value of top most image at position (x, y)

        :param float x: X position in plot coordinates
        :param float y: Y position in plot coordinates
        :return: The value at that point or '-'
        """
        value = "-"
        valueZ = -float("inf")
        mask = 0
        maskZ = -float("inf")

        for image in self.getAllImages():
            data = image.getData(copy=False)
            isMask = isinstance(image, items.MaskImageData)
            if isMask:
                zIndex = maskZ
            else:
                zIndex = valueZ
            if image.getZValue() >= zIndex:
                # This image is over the previous one
                ox, oy = image.getOrigin()
                sx, sy = image.getScale()
                row, col = (y - oy) / sy, (x - ox) / sx
                if row >= 0 and col >= 0:
                    # Test positive before cast otherwise issue with int(-0.5) = 0
                    row, col = int(row), int(col)
                    if row < data.shape[0] and col < data.shape[1]:
                        v, z = data[row, col], image.getZValue()
                        if not isMask:
                            value = v
                            valueZ = z
                        else:
                            mask = v
                            maskZ = z
        if maskZ > valueZ and mask > 0:
            return value, "Masked"
        return value

    def _getImageDims(self, *args):
        """GUI/CLI hint: return active-image dimensions for plot status text."""
        activeImage = self.getActiveImage()
        if activeImage is not None and activeImage.getData(copy=False) is not None:
            dims = activeImage.getData(copy=False).shape[1::-1]
            return "x".join(str(dim) for dim in dims)
        else:
            return "-"

    def getProfileToolbar(self):
        """Profile tools attached to this plot

        See :class:`silx.gui.plot.Profile.ProfileToolBar`
        """
        return self.profile

    # @deprecated(replacement="getProfilePlot", since_version="0.5.0")
    def getProfileWindow(self):
        """Return the profile plot window.

        .. note::
           CLI-capable, but normally used by GUI profile tools.
        """
        return self.getProfilePlot()

    def getProfilePlot(self):
        """Return plot window used to display profile curve.

        :return: :class:`Plot1D`
        """
        return self.profile.getProfilePlot()


class QImportScanCreator(qt.QDialog):
    def __init__(self, defaultMuTh, parent=None):
        """GUI-only: initialize the raw-image import scan setup dialog."""
        qt.QDialog.__init__(self, parent)
        self.defaultMuTh = defaultMuTh

        layout = qt.QGridLayout()

        layout.addWidget(qt.QLabel("scan axis:"), 0, 0)
        layout.addWidget(qt.QLabel("axis start:"), 1, 0)
        layout.addWidget(qt.QLabel("axis end:"), 2, 0)
        self.fixed_label = qt.QLabel("mu (fixed):")
        layout.addWidget(self.fixed_label, 3, 0)
        layout.addWidget(qt.QLabel("no frames:"), 5, 0)

        self.scanaxis = qt.QComboBox()
        self.scanaxis.addItem("theta")
        self.scanaxis.addItem("mu")
        self.scanaxis.setCurrentIndex(0)
        self.scanaxis.currentIndexChanged.connect(self.onScanAxisChanged)

        self.omstart = qt.QDoubleSpinBox()
        self.omstart.setRange(-180, 180)
        self.omstart.setDecimals(4)
        self.omstart.setSuffix(" °")
        self.omstart.setValue(-90)

        self.omend = qt.QDoubleSpinBox()
        self.omend.setRange(-180, 180)
        self.omend.setDecimals(4)
        self.omend.setSuffix(" °")
        self.omend.setValue(90)

        self.no = qt.QSpinBox()
        self.no.setReadOnly(True)
        self.no.setRange(1, 1000000000)
        self.no.setValue(180)

        self.fixedAngle = qt.QDoubleSpinBox()
        self.fixedAngle.setRange(-180, 180)
        self.fixedAngle.setValue(self.defaultMuTh[0])

        layout.addWidget(self.scanaxis, 0, 1)
        layout.addWidget(self.omstart, 1, 1)
        layout.addWidget(self.omend, 2, 1)
        layout.addWidget(self.no, 5, 1)
        layout.addWidget(self.fixedAngle, 3, 1)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        layout.addWidget(buttons, 6, 0, -1, -1)

        test = qt.QLabel("Parameters determined from loaded scan:")
        layout.addWidget(test, 4, 0, 1, 2)

        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)

        self.setLayout(layout)

    def onScanAxisChanged(self, index):
        """GUI-only: update fixed-angle fields for the selected scan axis."""
        if index == 0:
            self.omstart.setValue(-90.0)
            self.omend.setValue(90.0)
            self.fixed_label.setText("mu (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[0])

        elif index == 1:
            self.omstart.setValue(0.0)
            self.omend.setValue(15.0)
            self.fixed_label.setText("theta (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[1])


class QPlotDeleteWindow(qt.QDialog):
    def __init__(self, curveList, hidden, parent=None):
        """GUI-only: initialize a dialog for hiding or deleting curves."""
        qt.QDialog.__init__(self, parent)
        self.curves = curveList
        self.action = None

        layout = qt.QGridLayout()

        # create 'select all' button
        self.selectAllPlotsCheckbox = qt.QCheckBox()
        layout.addWidget(qt.QLabel("select all"), 0, 0, 1, 1)
        layout.addWidget(self.selectAllPlotsCheckbox, 0, 1, 1, 1)
        self.selectAllPlotsCheckbox.stateChanged.connect(self.checkOrUncheckAll)

        # create an entry (checkbox + name) for each plot curve
        self.boxes = []
        for i, j in enumerate(self.curves):
            newbox = qt.QCheckBox()
            if hidden[i]:
                newbox.setChecked(True)
            self.boxes.append(newbox)
            layout.addWidget(qt.QLabel(j), i + 1, 0, 1, 1)
            layout.addWidget(self.boxes[i], i + 1, 1, 1, 1)

        self.buttons = qt.QDialogButtonBox()
        self.buttons.addButton("Delete curves", self.buttons.ActionRole)
        self.buttons.addButton("Hide curves", self.buttons.ActionRole)
        self.buttons.addButton(qt.QDialogButtonBox.Cancel)

        if curveList == []:
            layout.addWidget(self.buttons, 0, 0, -1, -1)
        else:
            layout.addWidget(self.buttons, i + 2, 0, -1, -1)

        self.buttons.buttons()[1].clicked.connect(self.deleteClicked)
        self.buttons.buttons()[2].clicked.connect(self.hideClicked)
        self.buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)

        self.setLayout(layout)

    def deleteClicked(self):
        """GUI-only: accept the dialog with a delete-curves action."""
        self.action = "delete"
        self.accept()

    def hideClicked(self):
        """GUI-only: accept the dialog with a hide-curves action."""
        self.action = "hide"
        self.accept()

    def checkOrUncheckAll(self):
        """GUI-only: mirror the select-all checkbox to all curve checkboxes."""
        if self.selectAllPlotsCheckbox.isChecked():
            for i in self.boxes:
                i.blockSignals(True)
                i.setChecked(True)
                i.blockSignals(False)
        elif not self.selectAllPlotsCheckbox.isChecked():
            for i in self.boxes:
                i.blockSignals(True)
                i.setChecked(False)
                i.blockSignals(False)


class QScanCreator(qt.QDialog):
    def __init__(self, defaultMuTh, parent=None):
        """GUI-only: initialize the simulation scan setup dialog."""
        qt.QDialog.__init__(self, parent)
        self.defaultMuTh = defaultMuTh

        layout = qt.QGridLayout()

        layout.addWidget(qt.QLabel("scan axis:"), 0, 0)
        layout.addWidget(qt.QLabel("axis start:"), 1, 0)
        layout.addWidget(qt.QLabel("axis end:"), 2, 0)
        layout.addWidget(qt.QLabel("no points:"), 3, 0)
        self.fixed_label = qt.QLabel("mu (fixed):")
        layout.addWidget(self.fixed_label, 4, 0)

        self.scanaxis = qt.QComboBox()
        self.scanaxis.addItem("theta")
        self.scanaxis.addItem("mu")
        self.scanaxis.setCurrentIndex(0)
        self.scanaxis.currentIndexChanged.connect(self.onScanAxisChanged)

        self.omstart = qt.QDoubleSpinBox()
        self.omstart.setRange(-180, 180)
        self.omstart.setDecimals(4)
        self.omstart.setSuffix(" °")
        self.omstart.setValue(-90)

        self.omend = qt.QDoubleSpinBox()
        self.omend.setRange(-180, 180)
        self.omend.setDecimals(4)
        self.omend.setSuffix(" °")
        self.omend.setValue(90)

        self.no = qt.QSpinBox()
        self.no.setRange(1, 1000000000)
        self.no.setValue(180)

        self.fixedAngle = qt.QDoubleSpinBox()
        self.fixedAngle.setRange(-180, 180)
        self.fixedAngle.setValue(self.defaultMuTh[0])

        layout.addWidget(self.scanaxis, 0, 1)
        layout.addWidget(self.omstart, 1, 1)
        layout.addWidget(self.omend, 2, 1)
        layout.addWidget(self.no, 3, 1)
        layout.addWidget(self.fixedAngle, 4, 1)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        layout.addWidget(buttons, 5, 0, -1, -1)

        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)

        self.setLayout(layout)

    def onScanAxisChanged(self, index):
        """GUI-only: update fixed-angle fields for the selected scan axis."""
        if index == 0:
            self.omstart.setValue(-90.0)
            self.omend.setValue(90.0)
            self.fixed_label.setText("mu (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[0])

        elif index == 1:
            self.omstart.setValue(0.0)
            self.omend.setValue(15.0)
            self.fixed_label.setText("theta (fixed):")
            self.fixedAngle.setValue(self.defaultMuTh[1])


class QDiffractometerImageDialog(qt.QDialog):
    def __init__(self, parent=None):
        """GUI-only: initialize the diffraction-geometry image dialog."""
        qt.QDialog.__init__(self, parent)
        verticalLayout = qt.QVBoxLayout(self)
        verticalLayout.setContentsMargins(0, 0, 0, 0)
        img = qutils.AspectRatioPixmapLabel(self)
        pixmp = qt.QPixmap(resources.getDiffractometerPath())
        # img.setScaledContents(False)
        img.setPixmap(pixmp)

        verticalLayout.addWidget(img)

        # reader = qt.QImageReader()
        # img = reader.read()
        # view = qt.QGraphicsView(self)
        # svg = qt.QSvgWidget(self)
        # svg.load(resources.getDiffractometerPath())

        # verticalLayout.addWidget(svg)
        self.setLayout(verticalLayout)


class AboutDialog(qt.QDialog):
    def __init__(self, version, msg="", parent=None):
        """GUI-only: initialize the About dialog content."""
        qt.QDialog.__init__(self, parent)
        layout = qt.QVBoxLayout()
        self.setWindowTitle("About orGUI")

        pixmap = resources.getSplashScreen(str(version))
        self.logo = qt.QLabel()
        app = qt.QApplication.instance()
        screenGeometry = app.primaryScreen().availableGeometry()
        splashpm = pixmap.scaledToHeight(
            int(screenGeometry.height() / 5), qt.Qt.SmoothTransformation
        )
        self.logo.setPixmap(splashpm)

        messageStr = f"orGUI version {version}"
        messageStr += msg
        messageStr += """<br> <br>
Copyright (c) 2020-2026 Timo Fuchs, published under MIT License
<br> <br>
orGUI: Orientation and Integration with 2D detectors.<br>
Zenodo. <a href=\"https://doi.org/10.5281/zenodo.12592485\">https://doi.org/10.5281/zenodo.12592485</a> <br> <br>
New software updates will be published under <a href=\"https://doi.org/10.5281/zenodo.12592485\">Zenodo</a>.
<br> <br>
Help requests can be send via Email to Timo Fuchs.
<br> <br>
"orGUI" was developed during the PhD work of Timo Fuchs,<br>
within the group of Olaf Magnussen.
"""  # noqa: E501
        self.label = qt.QLabel()
        self.label.setText(messageStr)
        self.label.setTextInteractionFlags(qt.Qt.TextBrowserInteraction)
        self.label.setTextFormat(qt.Qt.RichText)

        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok)
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.accept)
        layout.addWidget(self.logo)
        layout.addWidget(self.label)
        layout.addWidget(buttons)
        self.setLayout(layout)


class UncaughtHook(qt.QObject):
    # _exception_caught = qt.Signal(object)

    def __init__(self, *args, **kwargs):
        """GUI-only: install the Qt-aware uncaught exception hook."""
        super().__init__(*args, **kwargs)

        # this registers the exception_hook() function as hook with the Python interpreter  # noqa: E501
        sys.excepthook = self.exception_hook

        self.orgui = None

        # connect signal to execute the message box function always on main thread
        # self._exception_caught.connect(show_exception_box)

    def set_orgui(self, orgui):
        """Attach the main window used for fatal-error recovery.

        .. note::
           GUI-only. The recovery path can show modal dialogs.
        """
        self.orgui = orgui

    def exception_hook(self, exc_type, exc_value, exc_traceback):
        """Handle uncaught exceptions through the GUI fatal-error path.

        .. note::
           GUI-only. This path may show modal dialogs and terminate the process.
        """
        if issubclass(exc_type, KeyboardInterrupt):
            # ignore keyboard interrupt to support console applications
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
        else:
            log_msg = "\n".join(
                [
                    "".join(traceback.format_tb(exc_traceback)),
                    f"{exc_type.__name__}: {exc_value}",
                ]
            )
            print(f"Uncaught exception:\n {log_msg}")  # exc_info=exc_info)

            # trigger message box show
            # self._exception_caught.emit(log_msg)

            if qt.QApplication.instance() is not None:
                if self.orgui is None:
                    errorbox = qt.QMessageBox(
                        qt.QMessageBox.Critical,
                        "Uncaught Exception",
                        f"An unexpected error occured. The program will terminate now:\n{log_msg}",  # noqa: E501
                        qt.QMessageBox.Ok,
                    )
                    errorbox.exec()
                    sys.exit(1)
                else:
                    resBtn = qutils.critical_detailed_message(
                        self.orgui,
                        "Uncaught Exception",
                        "An unexpected error has occured.\norGUI will terminate now.\nDo you want to try to save the database before terminating?",  # noqa: E501
                        log_msg,
                        qt.QMessageBox.Save | qt.QMessageBox.Discard,
                    )
                    # errorbox = qt.QMessageBox(qt.QMessageBox.Critical,
                    #                          "Uncaught Exception",
                    #                          "An unexpected error occured:\n{0}\nDo you want to try to save the data before terminating?".format(log_msg),  # noqa: E501
                    #                           qt.QMessageBox.Save | qt.QMessageBox.Discard)  # noqa: E501

                    # resBtn = errorbox.exec()

                    if resBtn == qt.QMessageBox.Save:
                        try:
                            self.orgui.database.onSaveDBFile()

                        except Exception:
                            print(
                                f"Fatal error: Cannot save database:\n{traceback.format_exc()}"  # noqa: E501
                            )
                            qutils.critical_detailed_message(
                                self.orgui,
                                "Fatal error",
                                "Cannot save database.",
                                traceback.format_exc(),
                            )
                    self.orgui.database.close()
            else:
                print("No QApplication instance available.")
            sys.exit(1)


def main(configfile):
    """Start a standalone GUI application for a config file.

    :param configfile:
        Path to the orGUI configuration file.

    .. note::
       GUI-only. CLI startup is handled by ``orgui.main``.
    """

    a = qt.QApplication(["orGUI"])

    qt_exception_hook = UncaughtHook()

    mainWindow = orGUI(configfile)
    qt_exception_hook.set_orgui(mainWindow)
    mainWindow.show()
    # a.lastWindowClosed.connect(a.quit)
    return a.exec_()


if __name__ == "__main__":
    main("./config")

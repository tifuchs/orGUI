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

import traceback
import logging
from silx.gui import qt
from silx.gui import icons
from .ArrayTableDialog import ArrayEditWidget
from dataclasses import dataclass
import copy

import numpy as np
from ..datautils.xrayutils import HKLVlieg

from .. import resources
from . import qutils
from . import imagePeakFinder

logger = logging.getLogger(__name__)


def _reflection_mismatch_score(mismatch, angular_tolerance, norm_tolerance):
    """Return absolute per-reflection mismatch scores.

    A score of ``1`` corresponds to one detector-pixel-equivalent mismatch.

    :param dict mismatch:
        Reflection mismatch dictionary from
        :meth:`HKLVlieg.UBCalculator.getReflectionMismatch`.
    :param float or numpy.ndarray angular_tolerance:
        Angular tolerance in rad.
    :param float or numpy.ndarray norm_tolerance:
        Q-norm tolerance in Angstrom^-1.
    :returns:
        Maximum of angular and Q-norm mismatch in tolerance units.
    :rtype: numpy.ndarray
    """
    angle = np.asarray(mismatch["angle_mismatch"], dtype=float)
    q_norm = np.asarray(mismatch["norm_mismatch"], dtype=float)
    angular_tolerance = np.maximum(
        np.asarray(angular_tolerance, dtype=float), np.finfo(float).eps
    )
    norm_tolerance = np.maximum(
        np.asarray(norm_tolerance, dtype=float), np.finfo(float).eps
    )
    return np.maximum(angle / angular_tolerance, q_norm / norm_tolerance)


def _local_detector_pixel_angular_tolerance(detector, xy, alpha_i):
    """Calculate local detector angular resolution at reflection pixels.

    :param DetectorCalibration.Detector2D_SXRD detector:
        Detector geometry used to convert pixels to surface angles.
    :param numpy.ndarray xy:
        Detector coordinates ``(x, y)`` in pixels, shaped ``(N, 2)``.
    :param numpy.ndarray alpha_i:
        Incidence angles in rad, shaped ``(N,)``.
    :returns:
        Per-reflection angular size of the finer local detector pixel axis in
        rad. The mismatch score is scalar, so rectangular pixels must use the
        stricter axis to avoid accepting errors that exceed one pixel along the
        higher-resolution detector direction.
    :rtype: numpy.ndarray
    """
    xy = np.atleast_2d(np.asarray(xy, dtype=float))
    alpha_i = np.asarray(alpha_i, dtype=float)
    x = xy[:, 0]
    y = xy[:, 1]
    gamma, delta = detector.surfaceAnglesPoint(y, x, alpha_i)
    gamma_x, delta_x = detector.surfaceAnglesPoint(y, x + 1.0, alpha_i)
    gamma_y, delta_y = detector.surfaceAnglesPoint(y + 1.0, x, alpha_i)
    horizontal = np.hypot(delta_x - delta, gamma_x - gamma)
    vertical = np.hypot(delta_y - delta, gamma_y - gamma)
    return np.minimum(horizontal, vertical)


def _relative_mismatch_score(angle, relative_norm):
    """Return the within-table relative mismatch score.

    :param numpy.ndarray angle:
        Angular reflection mismatch in rad.
    :param numpy.ndarray relative_norm:
        Relative Q-norm mismatch.
    :returns:
        Score normalized from best to worst current reflection.
    :rtype: numpy.ndarray
    """
    def normalized(values):
        finite = np.isfinite(values)
        result = np.ones(values.shape, dtype=float)
        if not np.any(finite):
            return result
        low = np.min(values[finite])
        high = np.max(values[finite])
        if np.isclose(low, high):
            result[finite] = 0.0
        else:
            result[finite] = (values[finite] - low) / (high - low)
        return result

    return normalized(0.5 * (normalized(angle) + normalized(relative_norm)))


@dataclass
class HKLReflection:
    xy: np.ndarray
    hkl: np.ndarray
    imageno: int
    identifier : str

    @classmethod
    def fromStr(cls,string):
        h,k,l,x,y,imageno = string.split()
        refl = cls([x,y],[h,k,l],int(imageno))
        return refl

    @staticmethod
    def fromArray(arr, identifier):
        reflist = []
        if arr.ndim == 1:
            h,k,l,x,y,imageno = arr
            refl = HKLReflection([x,y],[h,k,l],int(imageno), identifier)
            return refl
        if len(identifier) != len(arr):
            raise ValueError("Mismatch between number of reflections in array and identifiers.")
        for row, ident in zip(arr, identifier):
            h,k,l,x,y,imageno = row
            reflist.append(HKLReflection([x,y],[h,k,l],int(imageno), ident))
        return reflist

    def __array__(self):
        return np.concatenate([self.hkl,self.xy,[self.imageno]])

class QReflectionSelector(qt.QWidget):
    sigQueryImageChange = qt.pyqtSignal(int)
    sigQuerySaveReflections = qt.pyqtSignal()
    sigQueryLoadReflections = qt.pyqtSignal()
    sigQueryCenterPlot = qt.pyqtSignal(object)
    #sigImagePathChanged = qt.pyqtSignal(object)
    #sigImageNoChanged = qt.pyqtSignal(object)
    def __init__(self,plot, ubcalc, parent=None):
        qt.QWidget.__init__(self, parent=None)
        layout = qt.QVBoxLayout()
        #self.setOrientation(qt.Qt.Vertical)
        self.plot = plot
        self.ubcalc = ubcalc
        self.orparent = parent
        self.pkImgDiag = PeakImgRangeDialog([0, 360], 'NoAxis')
        self.peakReflDialog = PeakReflDialog(self)
        self.autoBraggStatusDialog = AutoBraggStatusDialog(self)
        self.autoBraggOptionsDialog = AutoBraggOptionsDialog(self)

        self.reflections = []
        self.reflBragg = []

        self.nextNo = 0
        self.imageno = 0

        self._showReferenceReflections = True
        self.showBraggReflections = False

        self.activeReflection = None

        self.plot.sigKeyPressDelete.connect(self._onDelete)

        editorTabWidget = qt.QTabWidget()

        self.refleditor = ArrayEditWidget(True, 1)
        header = ['H', 'K', 'L', 'x', 'y', 'imageno']
        self.refleditor.setArrayData(np.array([]).reshape((-1,6)), editable=True, header=header)
        #self.updateEditor()
        self.refleditor.model.dataChanged.connect(self._onDataChanged)
        self.refleditor.sigRowAdded.connect(self._onRowAdded)
        self.refleditor.sigDataLoaded.connect(self.reflectionsFromEditor)
        self.refleditor.sigRowsDeleted.connect(self._onRowsDeleted)
        self.refleditor.view.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.refleditor.view.selectionModel().currentRowChanged.connect(self._onSelectedRowChanged)

        editorTabWidget.addTab(self.refleditor, "intrinsic coord")

        self.refleditor_angles = ArrayEditWidget(True, -1, False)
        header = ['H', 'K', 'L', 'alpha', 'delta', 'gamma', 'omega', 'chi', 'phi']
        self.refleditor_angles.setArrayData(np.array([]).reshape((-1,9)), editable=False, header=header)
        self.refleditor_angles.view.setSelectionBehavior(qt.QAbstractItemView.SelectRows)
        self.refleditor_angles.view.selectionModel().currentRowChanged.connect(self._onSelectedRowChanged)
        editorTabWidget.addTab(self.refleditor_angles, "SIXC angles")

        layout.addWidget(editorTabWidget)
        self.mismatchLabel = qt.QLabel(
            "Mismatch: unavailable"
        )
        self.mismatchLabel.setToolTip(
            "Mean absolute mismatch across the reference reflections. "
            "The first value is the angular difference between Q calculated "
            "from the measured angles and Q calculated from UB, in degrees. "
            "The second value is the difference in Q magnitude, in Å⁻¹."
        )
        layout.addWidget(self.mismatchLabel)

        pixelMismatchLayout = qt.QHBoxLayout()
        pixelMismatchLayout.addWidget(qt.QLabel("Resolution limit:"))
        self.pixelMismatchLimit = qt.QDoubleSpinBox()
        self.pixelMismatchLimit.setDecimals(2)
        self.pixelMismatchLimit.setRange(0.1, 20.0)
        self.pixelMismatchLimit.setSingleStep(0.1)
        self.pixelMismatchLimit.setValue(1.0)
        self.pixelMismatchLimit.setSuffix(" px")
        self.pixelMismatchLimit.setToolTip(
            "Rows at or below this detector-pixel-equivalent mismatch are "
            "colored blue; other rows use relative green-to-red coloring."
        )
        self.pixelMismatchLimit.valueChanged.connect(
            lambda _value: self.ubcalc.updateReflectionMismatch()
        )
        pixelMismatchLayout.addWidget(self.pixelMismatchLimit)
        pixelMismatchLayout.addStretch(1)
        layout.addLayout(pixelMismatchLayout)
        #self.addWidget(editorTabWidget)

        controlsLayout = qt.QHBoxLayout()

        self.editControlsGroup = qt.QGroupBox("Reflection edit")
        editControlsLayout = qt.QHBoxLayout()
        editControlsLayout.setContentsMargins(6, 3, 6, 3)
        self.editToolbar = qt.QToolBar()
        self.editToolbar.setMovable(False)
        self.editToolbar.setOrientation(qt.Qt.Horizontal)
        editControlsLayout.addWidget(self.editToolbar)
        self.editControlsGroup.setLayout(editControlsLayout)

        self.autoControlsGroup = qt.QGroupBox("Auto UB/Reflections")
        autoControlsLayout = qt.QHBoxLayout()
        autoControlsLayout.setContentsMargins(6, 3, 6, 3)
        self.autoToolbar = qt.QToolBar()
        self.autoToolbar.setMovable(False)
        self.autoToolbar.setOrientation(qt.Qt.Horizontal)
        autoControlsLayout.addWidget(self.autoToolbar)
        self.autoControlsGroup.setLayout(autoControlsLayout)

        self.setImageAct = qt.QAction(resources.getQicon('select-image'), "select Image", self)
        self.setImageAct.setToolTip("set this image for the chosen reflection")
        self.setImageAct.triggered.connect(self._onSetImage)

        self.gotoImageAct = qt.QAction(resources.getQicon('search-image'), "search image",self)
        self.gotoImageAct.setToolTip("diplay image with the chosen reflection")
        self.gotoImageAct.triggered.connect(self._onGotoImage)

        self.peakSearchInImage = qt.QAction(resources.getQicon('search-peak'), "2D peak search",self)
        self.peakSearchInImage.setToolTip("perform peak search in the image")
        self.peakSearchInImage.triggered.connect(self._onSinglePeakSearch)

        self.addBraggReflAct = qt.QAction(
            resources.getQicon("add-bragg-reflection"),
            "add Bragg reflection",
            self,
        )
        self.addBraggReflAct.setToolTip("add the next Bragg reflection, which improves the U matrix the most")
        self.addBraggReflAct.triggered.connect(self._onAddBraggReflection)

        self.autoBraggSeedAct = qt.QAction(
            resources.getQicon("auto-bragg-ub"),
            "auto Bragg/UB",
            self,
        )
        self.autoBraggSeedAct.setToolTip(
            "open automatic UB/reflection controls"
        )
        self.autoBraggSeedAct.triggered.connect(self._onAutoBraggSeed)

        self.autoBraggOptionsAct = qt.QAction(
            resources.getQicon("auto-bragg-options"),
            "auto Bragg options",
            self,
        )
        self.autoBraggOptionsAct.setToolTip(
            "configure automatic Bragg search tolerances and HKL assignment"
        )
        self.autoBraggOptionsAct.triggered.connect(
            self.autoBraggOptionsDialog.show
        )

        self.editToolbar.addAction(self.gotoImageAct)
        self.editToolbar.addSeparator()
        self.editToolbar.addAction(self.setImageAct)
        self.editToolbar.addAction(self.peakSearchInImage)

        self.autoToolbar.addAction(self.autoBraggOptionsAct)
        self.autoToolbar.addAction(self.autoBraggSeedAct)
        self.autoToolbar.addAction(self.addBraggReflAct)

        controlsLayout.addWidget(self.editControlsGroup)
        controlsLayout.addStretch(1)
        controlsLayout.addWidget(self.autoControlsGroup)
        layout.addLayout(controlsLayout)

        self.setLayout(layout)

    # GUI-only: user-triggered automatic Bragg search and UB seed.
    def _onAutoBraggSeed(self):
        """Open automatic UB/reflection controls."""
        self.autoBraggStatusDialog.show()
        self.autoBraggStatusDialog.raise_()

    def _onSinglePeakSearch(self):

        if self.activeReflection is not None:
            try:
                refl = self.reflections[self.indexReflection(self.activeReflection)]
                if self.orparent.fscan is not None:
                    self.pkImgDiag.set_axis(self.orparent.fscan.axis, self.orparent.fscan.axisname)
                    if self.pkImgDiag.exec() == qt.QDialog.Accepted:
                        roidict = self.pkImgDiag.get_roi()
                        #mu, om = self.orparent.getMuOm(refl.imageno)
                        axis = self.orparent.fscan.axis
                        ax_min = roidict['from'] + axis[refl.imageno]
                        ax_max = roidict['to'] + axis[refl.imageno]
                        xy_start = refl.xy
                        com_d = self.searchPeak(xy_start,
                                                roidict['vsize'],
                                                roidict['hsize'],
                                                [ax_min, ax_max])
                        ax_min = roidict['from_2'] + com_d['axis_com']
                        ax_max = roidict['to_2'] + com_d['axis_com']
                        xy_start = com_d['xy']
                        com_d = self.searchPeak(xy_start,
                                                roidict['vsize_2'],
                                                roidict['hsize_2'],
                                                [ax_min, ax_max])

                        refl.xy = com_d['xy']
                        refl.imageno = np.argmin(np.abs(axis - com_d['axis_com']))
                        self.updateEditor()
                        self.sigQueryImageChange.emit(refl.imageno)

            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot calculte peak position:\n%s" % e, traceback.format_exc())

    def searchPeak(self, xy_start, vsize, hsize, ax_range):
        kwargs = {
            'excluded_images' : self.orparent.excludedImagesDialog.getData(),
            'max_workers' : self.orparent.numberthreads
        }
        if self.orparent.centralPlot.getMaskToolsDockWidget().getSelectionMask() is not None:
            kwargs['mask'] = self.orparent.centralPlot.getMaskToolsDockWidget().getSelectionMask() > 0.0
        com_d = imagePeakFinder.find_COM_Image(xy_start,
                                                vsize,
                                                hsize,
                                                self.orparent.fscan,
                                                ax_range,
                                                **kwargs)
        return com_d



    def getBraggCandidates(self, recalculate=False):
        if recalculate:
            self.orparent.calcBraggRefl()
        if not self.reflBragg:
            raise ValueError("No Bragg reflections calculated or available")

        refl_new = []

        if self.reflections == []:
            qutils.warning_detailed_message(self, "Error","No Bragg reflections set yet, place at least one reference reflection", '')
            return

        refl_hkl_current = np.vstack([r.hkl for r in self.reflections])

        for refl in self.reflBragg:
            for cref in self.reflections:
                if np.allclose(refl.hkl, cref.hkl):
                    break
            else:
                refl_new.append(refl)

        if refl_new == []:
            qutils.warning_detailed_message(self, "Error","No new Bragg reflections found",'')
            return

        refl_new_hkl_Bragg = np.vstack([r.hkl for r in refl_new])

        Q_current = self.orparent.ubcalc.ubCal.lattice.getB() @ refl_hkl_current.T
        Q_candidates = self.orparent.ubcalc.ubCal.lattice.getB() @ refl_new_hkl_Bragg.T

        refl_new_norm = []

        for Qc in Q_candidates.T:
            norm = np.linalg.norm(Q_current.T - Qc, axis=1)
            refl_new_norm.append(np.sum(norm))

        refl_new_norm = np.array(refl_new_norm)
        order = np.argsort(refl_new_norm)[::-1]

        ordered_refl = []
        for o in order:
            ordered_refl.append(refl_new[o])
        return ordered_refl

    def _onAddBraggReflection(self):
        if self.peakReflDialog.isVisible():
            self.peakReflDialog.activateWindow()
            return
        if not self.reflBragg:
            qutils.warning_detailed_message(self, "Error","No Bragg reflection positions calculated.\nUse view->show Bragg reflections to calculate them", '')
            return
        self.peakReflDialog.set_candidates(self.getBraggCandidates())
        self.peakReflDialog.show()


    def _onDataChanged(self, *dat):
        data = self.refleditor.getData()
        idx_no = dat[0].row()
        identifier = self.reflections[idx_no].identifier
        if data.ndim > 1:
            refl = HKLReflection.fromArray(data[idx_no], identifier)
        else:
            refl = HKLReflection.fromArray(data, identifier)
        self.reflections[idx_no] = refl
        if self.activeReflection != identifier:
            self.setReflectionActive(identifier)
        else:
            self.redrawActiveReflection()
        self.ubcalc.updateReflectionMismatch()

    def anglesReflection(self, refl):
        if self.orparent.fscan is None:
            return np.array([np.nan,np.nan,np.nan, np.nan, np.nan, np.nan])
        if not self.orparent.isValidImageNo(refl.imageno):
            logger.warning(
                "Skipping angle calculation for stale reflection %s with "
                "image number %s outside the active scan range.",
                refl.hkl,
                refl.imageno,
            )
            return np.array([np.nan,np.nan,np.nan, np.nan, np.nan, np.nan])
        mu, om = self.orparent.getMuOm(refl.imageno)

        gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(np.array([refl.xy[1]]),np.array([refl.xy[0]]), mu)
        pos = [mu,delta[0],gamma[0],om,self.ubcalc.chi,self.ubcalc.phi]
        pos = HKLVlieg.crystalAngles(pos,self.ubcalc.n)
        return np.array(pos)


    def _onRowAdded(self, row):
        identifier = 'ref_'+str(self.nextNo)
        if self.refleditor.getData().ndim > 1:
            refl = HKLReflection.fromArray(self.refleditor.getData()[row], identifier)
        else:
            refl = HKLReflection.fromArray(self.refleditor.getData(), identifier)
        self.nextNo += 1
        angles = np.concatenate([refl.hkl,np.rad2deg(self.anglesReflection(refl))])
        angledata = self.refleditor_angles.getData()
        if angledata.ndim > 1:
            angledata = np.insert(angledata, row, angles, axis=0)
        elif angledata.ndim == 1:
            if row == 0:
                angledata = np.vstack([angles, angledata])
            elif row == 1:
                angledata = np.vstack([angledata, angles])
        else:
            angledata = angles[np.newaxis,:]
        self.refleditor_angles.updateArrayData(angledata)

        if self._showReferenceReflections:
            self.plot.addMarker(*refl.xy,legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        self.reflections.insert(row, refl)
        self.setReflectionActive(refl.identifier)
        self.ubcalc.updateReflectionMismatch()

    def _onRowsDeleted(self, rows):
        idents = []
        for r in rows:
            idents.append(self.reflections[r].identifier)
        for ident in idents:
            self.deleteReflection(ident)
        """
        print(self.refleditor_angles.getData())
        angledata = np.atleast_2d(self.refleditor_angles.getData())
        print(angledata)
        mask = np.ones(len(angledata),dtype=bool)
        print(mask)
        mask[np.array(rows)] = False
        angledata = angledata[mask]
        print(angledata)
        print(angledata.ndim )
        if angledata.ndim == 1:
            angledata = angledata[np.newaxis,:]
        elif angledata.ndim == 0:
            angledata = np.array([]).reshape((-1,9))
        self.refleditor_angles.updateArrayData(angledata)
        """

    def _onSelectedRowChanged(self, selected, deselected):
        row = selected.row()
        refl = self.reflections[row]
        if self._showReferenceReflections:
            if self.activeReflection is not None:
                self.setReflectionInactive(self.activeReflection)
            self.plot.removeMarker(refl.identifier)
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
            self.activeReflection = refl.identifier
        else:
            self.activeReflection = refl.identifier

    def setBraggReflectionsVisible(self, visible):
        if visible:
            self.showBraggReflections = True
            self.redrawBraggReflections()
        else:
            self.showBraggReflections = False
            self.clearPlotBraggReflections()

    def setBraggReflections(self, hkls, yx, angles):
        self.__cached_bragg_refl = hkls, yx, angles
        if self.showBraggReflections:
            self.clearPlotBraggReflections()
        self.reflBragg.clear()
        for i, (hkl, yxr, pos) in enumerate(zip(hkls, yx, angles)):
            identifier = 'bragg_'+str(i)
            if self.orparent.fscan.axisname == 'mu':
                angle_idx = 0
                sign = 1.
            elif self.orparent.fscan.axisname == 'th':
                angle_idx = 3
                sign = -1.
            else:
                qt.QMessageBox.warning(self,"Cannot calculate reflection","Cannot calculate reflection.\n%s is no supported scan axis." % self.orparent.fscan.axisname)
                return
            try:
                imageno1 = self.orparent.axisToImageNo(np.rad2deg(pos[angle_idx]) * sign)
            except Exception as e:
                logger.warning(
                    "Skipping Bragg reflection %s because its scan-axis "
                    "position is outside the active scan range: %s",
                    hkl,
                    e,
                )
                continue

            self.reflBragg.append(HKLReflection(yxr[::-1], hkl, imageno1, identifier))
        if self.showBraggReflections:
            self.redrawBraggReflections()
        else:
            self.clearPlotBraggReflections()

    def getBraggReflections(self):
        if hasattr(self, "__cached_bragg_refl") and self.showBraggReflections:
            return self.__cached_bragg_refl
        else:
            raise ValueError("Bragg reflections are not calculated or are not shown!")

    def setReferenceReflectionsVisible(self, visible):
        if visible:
            self._showReferenceReflections = True
            self.redrawRefReflections()
        else:
            self._showReferenceReflections = False
            self.clearPlotRefReflections()

    def setImage(self,imageno):
        self.imageno = imageno

    def _onSaveReflections(self):
        self.sigQuerySaveReflections.emit()

    def _onLoadReflections(self):
        self.sigQueryLoadReflections.emit()


    def updateEditor(self):
        if self.reflections:
            array = np.vstack(self.reflections)
            self.refleditor.updateArrayData(array)
            anglearray = []
            for refl in self.reflections:
                anglearray.append(np.concatenate([refl.hkl,np.rad2deg(self.anglesReflection(refl))]))
            anglearray = np.vstack(anglearray)
            self.refleditor_angles.updateArrayData(anglearray)
        else:
            self.refleditor.updateArrayData(np.array([]).reshape((-1,6)))
            self.refleditor_angles.updateArrayData(np.array([]).reshape((-1,9)))
        try:
            idx = self.indexReflection(self.activeReflection)
            self.refleditor.view.selectRow(idx)
            self.refleditor_angles.view.selectRow(idx)
        except:
            if self.reflections:
                self.setReflectionActive(self.reflections[0].identifier)
            else:
                self.activeReflection = None
        self.ubcalc.updateReflectionMismatch()

    def setReflectionMismatch(self, mismatch):
        """Color reflection rows by UB disagreement.

        Rows that reach the detector-pixel-equivalent resolution threshold are
        colored blue. Remaining rows use relative green-to-red ranking across
        the current reference reflections.

        :param dict mismatch:
            Result from
            :meth:`HKLVlieg.UBCalculator.getReflectionMismatch`, or ``None``
            to clear the colors.
        """
        if mismatch is None:
            self.mismatchLabel.setText(
                "Mismatch: unavailable"
            )
            self.refleditor.model.setArrayColors()
            self.refleditor_angles.model.setArrayColors()
            self.refleditor.view.viewport().update()
            self.refleditor_angles.view.viewport().update()
            return

        angle = np.asarray(mismatch["angle_mismatch"], dtype=float)
        absolute_norm = np.asarray(
            mismatch["norm_mismatch"], dtype=float
        )
        relative_norm = np.asarray(
            mismatch["relative_norm_mismatch"], dtype=float
        )
        if (
            len(angle) != len(self.reflections)
            or len(absolute_norm) != len(self.reflections)
            or len(relative_norm) != len(self.reflections)
        ):
            self.mismatchLabel.setText(
                "Mismatch: unavailable"
            )
            self.refleditor.model.setArrayColors()
            self.refleditor_angles.model.setArrayColors()
            self.refleditor.view.viewport().update()
            self.refleditor_angles.view.viewport().update()
            return

        finite = np.isfinite(angle) & np.isfinite(absolute_norm)
        if np.any(finite):
            mean_angle = np.rad2deg(np.mean(angle[finite]))
            mean_norm = np.mean(absolute_norm[finite])
            self.mismatchLabel.setText(
                f"Mismatch: {mean_angle:.4g}°; "
                f"{mean_norm:.4g} Å⁻¹"
            )
        else:
            self.mismatchLabel.setText(
                "Mismatch: unavailable"
            )

        relative_score = np.clip(
            _relative_mismatch_score(angle, relative_norm), 0.0, 1.0
        )
        relative_score[~np.isfinite(relative_score)] = 1.0
        green = np.array([205.0, 245.0, 205.0])
        red = np.array([255.0, 190.0, 190.0])
        row_colors = (
            green + relative_score[:, np.newaxis] * (red - green)
        ).astype(np.uint8)
        # Convert the angular and Q-norm errors to detector-pixel-equivalent
        # units using the local detector-angle gradient at each reflection.
        # One pixel means the UB disagreement is at the instrumental angular
        # resolution for that detector position. This absolute "good enough"
        # test only overrides the relative green-to-red ranking for rows below
        # the GUI threshold.
        xy = np.asarray([refl.xy for refl in self.reflections], dtype=float)
        alpha_i = []
        orparent = getattr(self, "orparent", None)
        for refl in self.reflections:
            if (
                orparent is not None
                and orparent.isValidImageNo(refl.imageno)
            ):
                mu, _omega = orparent.getMuOm(refl.imageno)
            else:
                mu = self.ubcalc.mu
            alpha_i.append(mu)
        angular_tolerance = _local_detector_pixel_angular_tolerance(
            self.ubcalc.detectorCal, xy, np.asarray(alpha_i)
        )
        norm_tolerance = self.ubcalc.ubCal.getK() * angular_tolerance
        pixel_score = _reflection_mismatch_score(
            mismatch, angular_tolerance, norm_tolerance
        )
        limit_widget = getattr(self, "pixelMismatchLimit", None)
        pixel_limit = 1.0 if limit_widget is None else limit_widget.value()
        resolved = np.isfinite(pixel_score) & (pixel_score <= pixel_limit)
        row_colors[resolved] = np.array([185, 220, 255], dtype=np.uint8)

        def apply_colors(editor):
            table_shape = editor.model.getData().shape
            if len(table_shape) != 2 or table_shape[0] != len(row_colors):
                editor.model.setArrayColors()
            else:
                colors = np.broadcast_to(
                    row_colors[:, np.newaxis, :], table_shape + (3,)
                ).copy()
                editor.model.setArrayColors(bgcolors=colors)
            editor.view.viewport().update()

        apply_colors(self.refleditor)
        apply_colors(self.refleditor_angles)


    def reflectionsFromEditor(self):
        data = self.refleditor.getData()
        idents = ['ref_'+str(i) for i in range(len(data))]
        self.nextNo = len(data)
        try:
            if len(data) > 1:
                reflist = HKLReflection.fromArray(data, idents)
                self.setReflections(reflist)
            else:
                self.setReflections([])
        except Exception:
            qt.QMessageBox.critical(self,"Cannot read reflections", "Cannot read reflections %s" % traceback.format_exc())
            return

    def clearPlotRefReflections(self):
        for refl in self.reflections:
            self.plot.removeMarker(refl.identifier)

    def clearPlotBraggReflections(self):
        for refl in self.reflBragg:
            self.plot.removeMarker(refl.identifier)

    def redrawBraggReflections(self):
        self.clearPlotBraggReflections()
        if self.showBraggReflections:
            for refl in self.reflBragg:
                #refl = self.reflBragg[identifier]
                #text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl)
                self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,color='g',selectable=False,draggable=False,symbol='x', text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl))


    def redrawRefReflections(self):
        self.setReflections(copy.deepcopy(self.reflections))

    def setReflections(self,refls):
        self.clearPlotRefReflections()
        self.activeReflection = None
        self.nextNo = 0
        self.reflections.clear()
        for refl in refls:
            refl.identifier = 'ref_'+str(self.nextNo)
            if self._showReferenceReflections:
                self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')

            self.reflections.append(refl)
            self.nextNo += 1
        self.updateEditor()

    def addReflection(self,eventdict,imageno,hkl=np.array([np.nan,np.nan,np.nan])):
        identifier = 'ref_'+str(self.nextNo)
        refl = HKLReflection([eventdict['x'],eventdict['y']],hkl,imageno,identifier)
        if self._showReferenceReflections:
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=refl.identifier,text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')
        self.reflections.append(refl)
        self.nextNo += 1
        self.updateEditor()
        self.setReflectionActive(identifier)
        return refl


    def indexReflection(self, identifier):
        return [ref.identifier for ref in self.reflections].index(identifier)

    def setReflectionActive(self,identifier):
        idx = self.indexReflection(identifier)
        refl = self.reflections[idx]
        if self._showReferenceReflections:
            if self.activeReflection is not None:
                self.setReflectionInactive(self.activeReflection)
            self.plot.removeMarker(identifier)
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(identifier),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')
            self.activeReflection = identifier
        else:
            self.activeReflection = identifier
        selectModel = self.refleditor.view.selectionModel()
        if idx != selectModel.currentIndex().row():
            self.refleditor.view.selectRow(idx)
            self.refleditor_angles.view.selectRow(idx)

    def setReflectionInactive(self,identifier):
        if self._showReferenceReflections:
            self.plot.removeMarker(identifier)
            refl = self.reflections[self.indexReflection(identifier)]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(identifier),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='blue',selectable=True,draggable=True,symbol='.')

    def moveReflection(self,identifier,xy):
        self.reflections[self.indexReflection(identifier)].xy = xy
        self.updateEditor()

    def redrawActiveReflection(self):
        if self._showReferenceReflections:
            self.plot.removeMarker(self.activeReflection)
            refl = self.reflections[self.indexReflection(self.activeReflection)]
            self.plot.addMarker(refl.xy[0],refl.xy[1],legend=str(self.activeReflection),text="(%0.1f,%0.1f,%0.1f)" % tuple(refl.hkl),color='red',selectable=True,draggable=True,symbol='.')

    def _onGotoImage(self):
        if self.activeReflection is not None:
            self.sigQueryImageChange.emit(self.reflections[self.indexReflection(self.activeReflection)].imageno)
            self.sigQueryCenterPlot.emit(self.reflections[self.indexReflection(self.activeReflection)].xy)

    def _onSetImage(self):
        if self.activeReflection is not None:
            self.reflections[self.indexReflection(self.activeReflection)].imageno = self.imageno
            self.updateEditor()

    def _onDelete(self):
        try:
            idx = self.indexReflection(self.activeReflection)
        except:
            return
        refl = self.reflections[idx]
        btn = qt.QMessageBox.question(self, "Delete reflection?",
                                "Do you want to delete the reflection?\n%s" % refl.hkl,
                                qt.QMessageBox.Yes | qt.QMessageBox.No, qt.QMessageBox.Yes)
        if btn == qt.QMessageBox.Yes:
            self.deleteReflection(self.activeReflection)

    def deleteReflection(self,identifier):
        try:
            del self.reflections[self.indexReflection(identifier)]
        except:
            return

        if self._showReferenceReflections:
            self.plot.removeMarker(identifier)
        if self.activeReflection == identifier:
            self.activeReflection = None
            if self.reflections:
                newactive = self.reflections[0]
                self.setReflectionActive(newactive.identifier)
        self.updateEditor()





class QReflectionAnglesDialog(qt.QDialog):
    """Displays the angles for a reflection. Allows for multiple solutions.
    
    """

    class QSelectableLabel(qt.QLabel):
        def __init__(self, *args, **kwags):
            qt.QLabel.__init__(self,*args, **kwags)
            self.setTextInteractionFlags(qt.Qt.TextSelectableByKeyboard | qt.Qt.TextSelectableByMouse)

    def __init__(self,reflection, message, parent=None):
        qt.QDialog.__init__(self, parent=None)

        layout = qt.QVBoxLayout()

        layout.addWidget(qt.QLabel(message))

        layout.addWidget(QReflectionAnglesDialog.QSelectableLabel("Reflection: %s" % reflection['hkl']))

        reflectionLayout = qt.QGridLayout()

        available_data = [k.split('_')[0] for k in reflection.keys()]
        self.header = ""
        self.idx_dict = {}
        i = 1
        if 'xy' in available_data:
            for lbl in ['x', 'y']:
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel(lbl), 0, i)
                self.idx_dict[lbl] = i
                i += 1
                self.header += lbl + '\t'
        if 'angles' in available_data:
            for lbl in ['mu', 'delta', 'gamma', 'theta', 'chi', 'phi']:
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel(lbl), 0, i)
                self.idx_dict[lbl] = i
                i += 1
                self.header += lbl + '\t'
        if 'imageno' in available_data:
            for lbl in ['imageno']:
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel(lbl), 0, i)
                self.idx_dict[lbl] = i
                i += 1
                self.header += lbl + '\t'

        self.checkboxes = []
        self.datastr = ""
        number_refl = max([int(k.split('_')[-1]) for k in reflection.keys() if '_' in k])
        for i in range(1,number_refl+1):
            box = qt.QCheckBox()
            reflectionLayout.addWidget(box, i, 0)

            if not reflection.get('selectable_%s' % i, False):
                box.setEnabled(False)

            self.checkboxes.append(box)
            if 'xy_%s' % i in reflection.keys():
                xy = reflection['xy_%s' % i]
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % xy[0]), i, self.idx_dict['x'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % xy[1]), i, self.idx_dict['y'])
                self.datastr += "%.5f\t" % xy[0]
                self.datastr += "%.5f\t" % xy[1]
            elif 'xy' in available_data:
                self.datastr += "\t\t"

            if 'angles_%s' % i in reflection.keys():
                angles = np.rad2deg(reflection['angles_%s' % i])
                angles[3] *= -1. # to theta
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[0]), i, self.idx_dict['mu'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[1]), i, self.idx_dict['delta'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[2]), i, self.idx_dict['gamma'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[3]), i, self.idx_dict['theta'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[4]), i, self.idx_dict['chi'])
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%.5f" % angles[5]), i, self.idx_dict['phi'])
                for ang in angles:
                    self.datastr += "%.5f\t" % ang
            elif 'angles' in available_data:
                self.datastr += "\t\t\t\t\t\t"

            if 'imageno_%s' % i in reflection.keys():
                imageno = reflection['imageno_%s' % i]
                reflectionLayout.addWidget(QReflectionAnglesDialog.QSelectableLabel("%s" % imageno), i, self.idx_dict['imageno'])
                self.datastr += "%s" % imageno
            elif 'imageno' in available_data:
                self.datastr += "\t"
            self.datastr += "\n"
        layout.addLayout(reflectionLayout)


        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel,
                                      qt.Qt.Horizontal)
        copybutton = buttons.addButton("Copy to clipboad",qt.QDialogButtonBox.HelpRole)
        copybutton.clicked.connect(self.copy_to_clipboard)

        layout.addWidget(buttons)

        okbtn = buttons.button(qt.QDialogButtonBox.Ok)
        cancelbtn = buttons.button(qt.QDialogButtonBox.Cancel)

        okbtn.clicked.connect(self.accept)
        cancelbtn.clicked.connect(self.reject)

        self.setLayout(layout)

    def copy_to_clipboard(self):
        app = qt.QApplication.instance()
        app.clipboard().setText(self.header + "\n" + self.datastr)


class AutoBraggStatusDialog(qt.QDialog):
    """Show live automatic Bragg-search status and summary statistics."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.refl_selector = parent
        self.setWindowTitle("Automatic Bragg/UB search")
        self.setModal(False)
        self.stats = {}

        layout = qt.QVBoxLayout()
        self.summaryLabel = qt.QLabel("Automatic Bragg search not started.")
        layout.addWidget(self.summaryLabel)

        self.logText = qt.QPlainTextEdit()
        self.logText.setReadOnly(True)
        self.logText.setMinimumSize(520, 300)
        layout.addWidget(self.logText)

        seed_layout = qt.QHBoxLayout()
        self.seedUBBtn = qt.QPushButton("seed UB")
        self.seedUBBtn.setToolTip(
            "Search scan images for a Bragg peak, seed the UB matrix, and "
            "confirm it with a second Bragg peak."
        )
        self.seedUBBtn.clicked.connect(self._onSeedUB)
        seed_layout.addWidget(self.seedUBBtn)
        self.autoAddAfterSeed = qt.QCheckBox("add calculated after seed")
        self.autoAddAfterSeed.setChecked(True)
        self.autoAddAfterSeed.setToolTip(
            "After UB seeding succeeds, automatically run the calculated "
            "Bragg-reflection add pass below."
        )
        seed_layout.addWidget(self.autoAddAfterSeed)
        seed_layout.addStretch(1)
        layout.addLayout(seed_layout)

        add_layout = qt.QHBoxLayout()
        add_layout.addWidget(qt.QLabel("Additional reflections:"))
        self.additionalCount = qt.QSpinBox()
        self.additionalCount.setRange(1, 100000)
        self.additionalCount.setValue(50)
        self.additionalCount.setToolTip(
            "Number of additional calculated Bragg reflections to validate "
            "and add from the current UB matrix."
        )
        add_layout.addWidget(self.additionalCount)
        self.addAllBragg = qt.QCheckBox("add all Bragg")
        self.addAllBragg.setToolTip(
            "Ignore the number field and keep adding calculated Bragg "
            "reflections until no more candidates pass validation."
        )
        self.addAllBragg.toggled.connect(self.additionalCount.setDisabled)
        add_layout.addWidget(self.addAllBragg)
        self.addCalculatedBtn = qt.QPushButton("Add calculated")
        self.addCalculatedBtn.setToolTip(
            "Use the current UB matrix to predict Bragg reflections, refine "
            "their real peak positions, validate mismatch and intensity, and "
            "add accepted reflections."
        )
        self.addCalculatedBtn.clicked.connect(self._onAddCalculated)
        add_layout.addWidget(self.addCalculatedBtn)
        add_layout.addStretch(1)
        layout.addLayout(add_layout)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Reset | qt.QDialogButtonBox.Close
        )
        buttons.button(qt.QDialogButtonBox.Reset).clicked.connect(self.reset)
        buttons.rejected.connect(self.hide)
        buttons.button(qt.QDialogButtonBox.Close).clicked.connect(self.hide)
        layout.addWidget(buttons)
        self.setLayout(layout)

    def reset(self):
        """Clear previous search output."""
        self.stats = {
            "images": 0,
            "candidates": 0,
            "refined": 0,
            "hkl_hypotheses": 0,
            "confirmations": 0,
            "rejected": 0,
        }
        self.logText.clear()
        self.summaryLabel.setText("Automatic Bragg search running...")

    def add_status(self, event, **fields):
        """Append one automatic-search status event.

        :param str event:
            Short event name emitted by ``autoFindBraggReference``.
        :param kwargs fields:
            Event details to display in the dialog.
        """
        if "images_read" in fields:
            self.stats["images"] = fields["images_read"]
        if event == "candidate":
            self.stats["candidates"] += 1
        elif event == "refined":
            self.stats["refined"] += 1
        elif event == "hkl_hypotheses":
            self.stats["hkl_hypotheses"] += fields.get("count", 0)
        elif event == "confirmation":
            self.stats["confirmations"] += 1
        elif event in {"rejected", "failed", "error"}:
            self.stats["rejected"] += 1

        message = fields.pop("message", None)
        if message is None:
            details = ", ".join(
                f"{key}={value}" for key, value in fields.items()
            )
            message = f"{event}: {details}" if details else event
        self.logText.appendPlainText(message)
        self.summaryLabel.setText(
            "Images read: {images}; candidates: {candidates}; refined: "
            "{refined}; HKL hypotheses: {hkl_hypotheses}; confirmations: "
            "{confirmations}; rejected: {rejected}".format(**self.stats)
        )
        app = qt.QApplication.instance()
        if app is not None:
            app.processEvents()

    def _calculated_count(self):
        """Return requested number of calculated Bragg reflections."""
        if self.addAllBragg.isChecked():
            return 100000
        return self.additionalCount.value()

    # GUI-only: user-triggered automatic Bragg search and UB seed.
    def _onSeedUB(self):
        """Search the active scan for a Bragg peak and seed U."""
        if self.refl_selector is None:
            return
        options_dialog = getattr(
            self.refl_selector, "autoBraggOptionsDialog", None
        )
        options = {} if options_dialog is None else options_dialog.get_options()
        self.reset()
        self.seedUBBtn.setEnabled(False)
        self.addCalculatedBtn.setEnabled(False)
        try:
            seed = self.refl_selector.orparent.autoFindBraggReference(
                **options,
                status_callback=self.add_status,
            )
        except Exception as e:
            self.add_status("error", message=f"Error: {e}")
            qutils.warning_detailed_message(
                self,
                "Cannot auto-find Bragg reflection",
                f"Cannot auto-find Bragg reflection:\n{e}",
                traceback.format_exc(),
            )
            return
        finally:
            self.seedUBBtn.setEnabled(True)
            self.addCalculatedBtn.setEnabled(True)
        if seed is None:
            self.add_status(
                "failed", message="No reliable Bragg reflection found."
            )
            qutils.warning_detailed_message(
                self,
                "No Bragg reflection found",
                "No reliable Bragg reflection candidate was found.",
                "",
            )
            return
        self.add_status(
            "accepted",
            message=(
                "Accepted seed hkl=%s at image %s, xy=(%.2f, %.2f)"
                % (
                    seed.hkl,
                    seed.peak.imageno,
                    seed.peak.xy[0],
                    seed.peak.xy[1],
                )
            ),
        )
        if self.autoAddAfterSeed.isChecked():
            self.seedUBBtn.setEnabled(False)
            self.addCalculatedBtn.setEnabled(False)
            try:
                self._addCalculatedFromDialog(options)
            except Exception as e:
                self.add_status("error", message=f"Error: {e}")
                qutils.warning_detailed_message(
                    self,
                    "Cannot add calculated Bragg reflections",
                    f"Cannot add calculated Bragg reflections:\n{e}",
                    traceback.format_exc(),
                )
            finally:
                self.seedUBBtn.setEnabled(True)
                self.addCalculatedBtn.setEnabled(True)

    def _addCalculatedFromDialog(self, options):
        """Run calculated Bragg-reflection addition from dialog controls."""
        added = self.refl_selector.orparent.autoAddCalculatedBraggReflections(
            self._calculated_count(),
            status_callback=self.add_status,
            **options,
        )
        self.add_status(
            "accepted",
            message="Added %s calculated Bragg reflection(s)." % added,
        )
        return added

    # GUI-only: user-triggered automatic validation of calculated reflections.
    def _onAddCalculated(self):
        """Add calculated Bragg reflections using the current UB matrix."""
        if self.refl_selector is None:
            return
        options_dialog = getattr(
            self.refl_selector, "autoBraggOptionsDialog", None
        )
        options = {} if options_dialog is None else options_dialog.get_options()
        self.addCalculatedBtn.setEnabled(False)
        self.seedUBBtn.setEnabled(False)
        try:
            self._addCalculatedFromDialog(options)
        except Exception as e:
            self.add_status("error", message=f"Error: {e}")
            qutils.warning_detailed_message(
                self,
                "Cannot add calculated Bragg reflections",
                f"Cannot add calculated Bragg reflections:\n{e}",
                traceback.format_exc(),
            )
            return
        finally:
            self.addCalculatedBtn.setEnabled(True)
            self.seedUBBtn.setEnabled(True)


class AutoBraggOptionsDialog(qt.QDialog):
    """Configure automatic Bragg-search tolerances."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Automatic Bragg options")
        layout = qt.QGridLayout()

        row = 0
        layout.addWidget(qt.QLabel("HKL assignment:"), row, 0)
        self.hklMode = qt.QComboBox()
        self.hklMode.addItem("Detector position (fast)", "detector_position")
        self.hklMode.addItem("Q norm shell (slower)", "qnorm")
        self.hklMode.setToolTip(
            "Detector position is the default. Q norm shell tests more HKL "
            "hypotheses and can take longer."
        )
        layout.addWidget(self.hklMode, row, 1, 1, 3)

        row += 1
        self.qnormTolerance = self._add_double(
            layout,
            row,
            "Q shell half-width:",
            0.05,
            0.0001,
            10.0,
            " Angstrom^-1",
            "Allowed local Q-norm shell around the measured peak.",
        )
        row += 1
        self.assignmentPixelTolerance = self._add_double(
            layout,
            row,
            "Assignment pixel tolerance:",
            30.0,
            0.1,
            10000.0,
            " px",
            "Maximum first-peak detector-position mismatch.",
        )
        row += 1
        self.confirmationPixelTolerance = self._add_double(
            layout,
            row,
            "Confirmation pixel tolerance:",
            8.0,
            0.1,
            10000.0,
            " px",
            "Maximum second-peak detector-position mismatch.",
        )
        row += 1
        self.confirmationImageTolerance = self._add_int(
            layout,
            row,
            "Confirmation image tolerance:",
            3,
            0,
            100000,
            " images",
            "Maximum image-index mismatch for second-peak confirmation.",
        )
        row += 1
        self.adaptiveAfterCandidates = self._add_int(
            layout,
            row,
            "Broaden after unmatched:",
            5,
            1,
            100000,
            " candidates",
            "Switch to broader tolerances after this many refined candidates "
            "fail HKL assignment or confirmation.",
        )
        row += 1
        self.adaptiveQnormTolerance = self._add_double(
            layout,
            row,
            "Broader Q shell half-width:",
            0.25,
            0.0001,
            10.0,
            " Angstrom^-1",
            "Q-norm shell used after the unmatched-candidate threshold.",
        )
        row += 1
        self.adaptiveAssignmentPixelTolerance = self._add_double(
            layout,
            row,
            "Broader assignment tolerance:",
            100.0,
            0.1,
            10000.0,
            " px",
            "First-peak delta/gamma detector-position tolerance used after "
            "the unmatched-candidate threshold.",
        )
        row += 1
        self.adaptiveConfirmationPixelTolerance = self._add_double(
            layout,
            row,
            "Broader confirmation tolerance:",
            30.0,
            0.1,
            10000.0,
            " px",
            "Second-peak delta/gamma detector-position tolerance used after "
            "the unmatched-candidate threshold.",
        )
        row += 1
        self.adaptiveConfirmationImageTolerance = self._add_int(
            layout,
            row,
            "Broader confirmation image tolerance:",
            8,
            0,
            100000,
            " images",
            "Second-peak image-index tolerance used after the "
            "unmatched-candidate threshold.",
        )
        row += 1
        layout.addWidget(qt.QLabel("Scale-fit detector filter:"), row, 0)
        self.adaptiveScaleDetectorFilter = qt.QCheckBox("enabled")
        self.adaptiveScaleDetectorFilter.setChecked(True)
        self.adaptiveScaleDetectorFilter.setToolTip(
            "For detector-position HKL assignment, estimate the Q scale only "
            "from reflections whose current-U detector prediction lands near "
            "the observed candidate peak."
        )
        layout.addWidget(self.adaptiveScaleDetectorFilter, row, 1)
        self.adaptiveScaleDetectorFraction = self._add_double(
            layout,
            row,
            "region half-size:",
            0.25,
            0.001,
            10.0,
            " detector",
            "Half-width and half-height of the accepted detector region as a "
            "fraction of detector size.",
            column=2,
        )
        row += 1
        self.adaptiveScaleOutlierQTolerance = self._add_double(
            layout,
            row,
            "Scale-fit Q outlier tolerance:",
            0.25,
            0.0001,
            10.0,
            " Angstrom^-1",
            "Maximum scaled Q-norm residual for provisional reflections used "
            "in the adaptive scale-fit retry.",
        )
        self.adaptiveScaleOutlierAngleFraction = self._add_double(
            layout,
            row,
            "delta/gamma fraction:",
            0.25,
            0.001,
            10.0,
            " detector",
            "Maximum delta/gamma residual as a fraction of the detector "
            "angular range for provisional scale-fit reflections.",
            column=2,
        )
        row += 1
        self.axisHalfWidth = self._add_double(
            layout,
            row,
            "Peak search axis half-width:",
            1.0,
            0.001,
            100000.0,
            " deg",
            "First-pass local 3D peak-search half-width on scan axis.",
        )
        row += 1
        self.fineAxisHalfWidth = self._add_double(
            layout,
            row,
            "Fine axis half-width:",
            0.4,
            0.001,
            100000.0,
            " deg",
            "Second-pass local 3D peak-search half-width on scan axis.",
        )
        row += 1
        self.roiVSize = self._add_int(
            layout,
            row,
            "Peak search ROI vertical:",
            80,
            1,
            100000,
            " px",
            "First-pass local peak-search vertical ROI size.",
        )
        self.roiHSize = self._add_int(
            layout,
            row,
            "horizontal:",
            80,
            1,
            100000,
            " px",
            "First-pass local peak-search horizontal ROI size.",
            column=2,
        )
        row += 1
        self.fineRoiVSize = self._add_int(
            layout,
            row,
            "Fine ROI vertical:",
            40,
            1,
            100000,
            " px",
            "Second-pass local peak-search vertical ROI size.",
        )
        self.fineRoiHSize = self._add_int(
            layout,
            row,
            "horizontal:",
            40,
            1,
            100000,
            " px",
            "Second-pass local peak-search horizontal ROI size.",
            column=2,
        )
        row += 1
        self.assignmentReflections = self._add_int(
            layout,
            row,
            "Assignment reflections:",
            40,
            1,
            100000,
            "",
            "Maximum first-peak HKL hypotheses from the local Q shell.",
        )
        row += 1
        self.maxSeedHypotheses = self._add_int(
            layout,
            row,
            "Seed hypotheses:",
            12,
            1,
            100000,
            "",
            "Maximum seed HKL hypotheses sent to second-peak confirmation.",
        )
        row += 1
        self.confirmationReflections = self._add_int(
            layout,
            row,
            "Confirmation reflections:",
            12,
            1,
            100000,
            "",
            "Maximum predicted strong reflections tested for confirmation.",
        )
        row += 1
        layout.addWidget(qt.QLabel("Intensity validation:"), row, 0)
        self.intensityRatioCheck = qt.QCheckBox("match calculated |F|^2 ratio")
        self.intensityRatioCheck.setChecked(True)
        self.intensityRatioCheck.setToolTip(
            "When enabled, require seed/confirmation integrated rocking "
            "intensity ratio to match the calculated structure-factor "
            "intensity ratio within 50%."
        )
        layout.addWidget(self.intensityRatioCheck, row, 1)
        self.confirmationProminenceThreshold = self._add_double(
            layout,
            row,
            "prominence fallback:",
            6.0,
            0.0,
            100000.0,
            " z",
            "If calculated intensity-ratio matching is disabled, require both "
            "rocking curves to exceed this robust side-band prominence.",
            column=2,
        )

        row += 1
        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Close)
        buttons.button(qt.QDialogButtonBox.Close).clicked.connect(self.hide)
        layout.addWidget(buttons, row, 0, 1, 4)
        self.setLayout(layout)

    def _add_double(
        self,
        layout,
        row,
        label,
        value,
        minimum,
        maximum,
        suffix,
        tooltip,
        column=0,
    ):
        layout.addWidget(qt.QLabel(label), row, column)
        box = qt.QDoubleSpinBox()
        box.setDecimals(4)
        box.setRange(minimum, maximum)
        box.setValue(value)
        box.setSuffix(suffix)
        box.setToolTip(tooltip)
        layout.addWidget(box, row, column + 1)
        return box

    def _add_int(
        self,
        layout,
        row,
        label,
        value,
        minimum,
        maximum,
        suffix,
        tooltip,
        column=0,
    ):
        layout.addWidget(qt.QLabel(label), row, column)
        box = qt.QSpinBox()
        box.setRange(minimum, maximum)
        box.setValue(value)
        box.setSuffix(suffix)
        box.setToolTip(tooltip)
        layout.addWidget(box, row, column + 1)
        return box

    def get_options(self):
        """Return options for ``orGUI.autoFindBraggReference``."""
        return {
            "hkl_candidate_mode": self.hklMode.currentData(),
            "qnorm_tolerance": self.qnormTolerance.value(),
            "assignment_pixel_tolerance": (
                self.assignmentPixelTolerance.value()
            ),
            "confirmation_pixel_tolerance": (
                self.confirmationPixelTolerance.value()
            ),
            "confirmation_image_tolerance": (
                self.confirmationImageTolerance.value()
            ),
            "adaptive_after_candidates": (
                self.adaptiveAfterCandidates.value()
            ),
            "adaptive_qnorm_tolerance": (
                self.adaptiveQnormTolerance.value()
            ),
            "adaptive_assignment_pixel_tolerance": (
                self.adaptiveAssignmentPixelTolerance.value()
            ),
            "adaptive_confirmation_pixel_tolerance": (
                self.adaptiveConfirmationPixelTolerance.value()
            ),
            "adaptive_confirmation_image_tolerance": (
                self.adaptiveConfirmationImageTolerance.value()
            ),
            "adaptive_scale_detector_filter": (
                self.adaptiveScaleDetectorFilter.isChecked()
            ),
            "adaptive_scale_detector_fraction": (
                self.adaptiveScaleDetectorFraction.value()
            ),
            "adaptive_scale_outlier_q_tolerance": (
                self.adaptiveScaleOutlierQTolerance.value()
            ),
            "adaptive_scale_outlier_angle_fraction": (
                self.adaptiveScaleOutlierAngleFraction.value()
            ),
            "axis_half_width": self.axisHalfWidth.value(),
            "fine_axis_half_width": self.fineAxisHalfWidth.value(),
            "roi_size": (self.roiVSize.value(), self.roiHSize.value()),
            "fine_roi_size": (
                self.fineRoiVSize.value(),
                self.fineRoiHSize.value(),
            ),
            "assignment_reflections": self.assignmentReflections.value(),
            "max_seed_hypotheses": self.maxSeedHypotheses.value(),
            "confirmation_reflections": self.confirmationReflections.value(),
            "intensity_ratio_check": self.intensityRatioCheck.isChecked(),
            "confirmation_prominence_threshold": (
                self.confirmationProminenceThreshold.value()
            ),
        }



class PeakReflDialog(qt.QDialog):


    def __init__(self, refl_selector):
        qt.QDialog.__init__(self, refl_selector)
        self.refl_selector = refl_selector
        self.candidate_reflections = []
        self.currentRefl = None
        self.idx = 0

        layout = qt.QVBoxLayout()

        toplayout = qt.QHBoxLayout()
        hkl = np.array([np.nan, np.nan, np.nan])
        self._hkl_label = qt.QLabel("New reflection: %s (%s/%s)" % (hkl, 0,0))
        toplayout.addWidget(self._hkl_label)


        self.refreshBtn = qt.QPushButton()
        self.refreshBtn.setIcon(icons.getQIcon('view-refresh'))
        self.refreshBtn.setToolTip("Refresh and reorder Bragg reflection list")
        self.refreshBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.refreshBtn.clicked.connect(self._onRefreshBragg)
        toplayout.addWidget(self.refreshBtn)

        self.gotoImageBtn = qt.QPushButton("search image")
        self.gotoImageBtn.setIcon(resources.getQicon('search-image'))
        self.gotoImageBtn.setToolTip("diplay image with the chosen reflection")
        self.gotoImageBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.gotoImageBtn.clicked.connect(self._onGotoImage)
        toplayout.addWidget(self.gotoImageBtn)

        self.next_btn = qt.QPushButton('next')
        self.next_btn.setIcon(icons.getQIcon('next'))
        self.next_btn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.next_btn.setToolTip("Goto next Bragg reflection")
        self.next_btn.clicked.connect(self._onNextBragg)
        toplayout.addWidget(self.next_btn)

        self.prev_btn = qt.QPushButton('previous')
        self.prev_btn.setIcon(icons.getQIcon('previous'))
        self.prev_btn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.prev_btn.setToolTip("Goto previous Bragg reflection")
        self.prev_btn.clicked.connect(self._onPrevBragg)
        toplayout.addWidget(self.prev_btn)

        layout.addLayout(toplayout)

        editlayout = qt.QHBoxLayout()

        editlayout.addWidget(qt.QLabel("Edit reflection: "))

        self.setImageBtn = qt.QPushButton("select Image")
        self.setImageBtn.setIcon(resources.getQicon('select-image'))
        self.setImageBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.setImageBtn.setToolTip("set this image for the chosen reflection")
        self.setImageBtn.clicked.connect(self._onSelectImg)
        editlayout.addWidget(self.setImageBtn)

        self.peakSearchBtn = qt.QPushButton("2D peak search")
        self.peakSearchBtn.setIcon(resources.getQicon('search-peak'))
        self.peakSearchBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.peakSearchBtn.setToolTip("perform peak search in the image")
        self.peakSearchBtn.clicked.connect(self._onPeakSearch)
        editlayout.addWidget(self.peakSearchBtn)

        editlayout.addWidget(qt.QLabel("Accept reflection:"))


        self.acceptReflBtn = qt.QPushButton('Add reflection')
        self.acceptReflBtn.setIcon(icons.getQIcon('selected'))
        self.acceptReflBtn.setSizePolicy(qt.QSizePolicy(qt.QSizePolicy.Fixed, qt.QSizePolicy.Minimum))
        self.acceptReflBtn.setToolTip("add reflection to the list of reference reflections")
        self.acceptReflBtn.clicked.connect(self._onAcceptRefl)
        editlayout.addWidget(self.acceptReflBtn)

        layout.addLayout(editlayout)

        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)

        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.onCancel)

        layout.addWidget(buttons)

        self.setLayout(layout)

    def _onGotoImage(self):
        if self.currentRefl is not None:
            try:
                self.refl_selector.setReflectionActive(self.currentRefl.identifier)
                self.refl_selector.sigQueryImageChange.emit(self.currentRefl.imageno)
                self.refl_selector.sigQueryCenterPlot.emit(self.currentRefl.xy)
            except Exception as e:
                print("Cannot goto reflection %s: %s" % (self.currentRefl, e))


    def _onNextBragg(self):
        if self.currentRefl is not None and self.candidate_reflections:
            if self.idx < len(self.candidate_reflections)-1:
                self.set_current_candidate(self.idx + 1)
            else:
                qt.QMessageBox.warning(self, "Max Bragg reflections reached", "The maximum of available Bragg reflections has been reached.")

    def _onPrevBragg(self):
        if self.currentRefl is not None and self.candidate_reflections:
            if self.idx > 0:
                self.set_current_candidate(self.idx - 1)
            else:
                qt.QMessageBox.warning(self, "Min Bragg reflections reached", "The lowest of available Bragg reflections has been reached.")

    def _onSelectImg(self):
        if self.currentRefl is not None and self.candidate_reflections:
            try:
                self.refl_selector.setReflectionActive(self.currentRefl.identifier)
                self.refl_selector._onSetImage()
            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot select image number:\n%s" % e, traceback.format_exc())

    def _onRefreshBragg(self):
        try:
            refl_list = self.refl_selector.getBraggCandidates(True)
            self.set_candidates(refl_list)
        except Exception as e:
            qutils.warning_detailed_message(self, "Error","Cannot refresh Bragg reflections:\n%s" % e, traceback.format_exc())

    def _onPeakSearch(self):
        if self.currentRefl is not None and self.candidate_reflections:
            try:
                self.refl_selector.setReflectionActive(self.currentRefl.identifier)
                self.refl_selector._onSinglePeakSearch()
            except Exception as e:
                qutils.warning_detailed_message(self, "Error","Cannot perform peak search:\n%s" % e, traceback.format_exc())

    def _onAcceptRefl(self):
        if self.currentRefl is not None and self.candidate_reflections:
            del self.candidate_reflections[self.idx]
            self.currentRefl = None
            if self.idx < len(self.candidate_reflections) -1:
                self.set_current_candidate(self.idx)
            elif self.candidate_reflections:
                self.set_current_candidate(self.idx-1)
            else:
                self.idx = 0
                hkl = np.array([np.nan, np.nan, np.nan])
                self._hkl_label.setText("New reflection: %s (%s/%s)" % (hkl, 0,0))


    def set_current_candidate(self, idx):
        refl = self.candidate_reflections[idx]
        if self.currentRefl is not None:
            try:
                self.refl_selector.deleteReflection(self.currentRefl.identifier)
            except Exception as e:
                print("Cannot remove reflection %s: %s" % (self.currentRefl, e))
        dd = {'x' : refl.xy[0], 'y' : refl.xy[1]}
        refl = self.refl_selector.addReflection(dd,refl.imageno,hkl=refl.hkl)
        self.candidate_reflections[idx] = refl
        self.currentRefl = refl
        self.idx = idx
        self._hkl_label.setText("New reflection: %s (%s/%s)" % (refl.hkl, idx+1, len(self.candidate_reflections)))
        self.refl_selector.sigQueryImageChange.emit(refl.imageno)
        self.refl_selector.sigQueryCenterPlot.emit(refl.xy)


    def set_candidates(self,refl_list):
        if refl_list:
            self.candidate_reflections = copy.deepcopy(refl_list)
            self.set_current_candidate(0)
        else:
            if self.currentRefl is not None:
                try:
                    self.refl_selector.deleteReflection(self.currentRefl.identifier)
                except Exception as e:
                    print("Cannot remove reflection %s: %s" % (self.currentRefl, e))
            self.idx = 0
            self.currentRefl = None
            self.candidate_reflections = []
            hkl = np.array([np.nan, np.nan, np.nan])
            self._hkl_label.setText("New reflection: %s (%s/%s)" % (hkl, 0,0))


    def onOk(self):
        self.set_candidates([])
        self.accept()

    def onCancel(self):
        self.set_candidates([])
        self.reject()

class PeakImgRangeDialog(qt.QDialog):


    def __init__(self,axis, axisname, parent=None):
        qt.QDialog.__init__(self, parent)

        layout = qt.QGridLayout()

        self._scan_label = qt.QLabel("%s-scan: from %s to %s, max image will be calculated relative to set peak position" % (axisname, axis[0], axis[-1]))
        layout.addWidget(self._scan_label,0,0, 1, -1)

        layout.addWidget(qt.QLabel("1st pass (rough alignment):"),1,0, 1, -1)

        layout.addWidget(qt.QLabel("max img range:"),2, 0)
        layout.addWidget(qt.QLabel("from (rel):"),2, 1)
        layout.addWidget(qt.QLabel("to (rel):"),2, 3)

        self.roi_from = qt.QDoubleSpinBox()
        self.roi_from.setRange(-1000000, 1000000)
        self.roi_from.setDecimals(2)
        self.roi_from.setSuffix(" °")
        self.roi_from.setValue(-2.0)
        layout.addWidget(self.roi_from, 2,2)

        self.roi_to = qt.QDoubleSpinBox()
        self.roi_to.setRange(-1000000, 1000000)
        self.roi_to.setDecimals(2)
        self.roi_to.setSuffix(" °")
        self.roi_to.setValue(2.0)
        layout.addWidget(self.roi_to, 2,4)


        layout.addWidget(qt.QLabel("ROI size:"),3, 0)
        layout.addWidget(qt.QLabel("size v:"),3, 1)
        layout.addWidget(qt.QLabel("size h:"),3, 3)

        self.vsize = qt.QSpinBox()
        self.vsize.setRange(1, 1000000)
        self.vsize.setSuffix(" pix")
        self.vsize.setValue(100)
        layout.addWidget(self.vsize, 3,2)

        self.hsize = qt.QSpinBox()
        self.hsize.setRange(1, 1000000)
        self.hsize.setSuffix(" pix")
        self.hsize.setValue(100)
        layout.addWidget(self.hsize, 3,4)



        layout.addWidget(qt.QLabel("2nd pass (fine alignment):"),4,0, 1, -1)

        layout.addWidget(qt.QLabel("max img range:"),5, 0)
        layout.addWidget(qt.QLabel("from (rel):"),5, 1)
        layout.addWidget(qt.QLabel("to (rel):"),5, 3)

        self.roi_from_2 = qt.QDoubleSpinBox()
        self.roi_from_2.setRange(-1000000, 1000000)
        self.roi_from_2.setDecimals(2)
        self.roi_from_2.setSuffix(" °")
        self.roi_from_2.setValue(-1.0)
        layout.addWidget(self.roi_from_2, 5,2)

        self.roi_to_2 = qt.QDoubleSpinBox()
        self.roi_to_2.setRange(-1000000, 1000000)
        self.roi_to_2.setDecimals(2)
        self.roi_to_2.setSuffix(" °")
        self.roi_to_2.setValue(1.0)
        layout.addWidget(self.roi_to_2, 5,4)


        layout.addWidget(qt.QLabel("ROI size:"),6, 0)
        layout.addWidget(qt.QLabel("size v:"),6, 1)
        layout.addWidget(qt.QLabel("size h:"),6, 3)

        self.vsize_2 = qt.QSpinBox()
        self.vsize_2.setRange(1, 1000000)
        self.vsize_2.setSuffix(" pix")
        self.vsize_2.setValue(60)
        layout.addWidget(self.vsize_2, 6,2)

        self.hsize_2 = qt.QSpinBox()
        self.hsize_2.setRange(1, 1000000)
        self.hsize_2.setSuffix(" pix")
        self.hsize_2.setValue(60)
        layout.addWidget(self.hsize_2, 6,4)

        buttons = qt.QDialogButtonBox(qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel)
        layout.addWidget(buttons,7,0,1,-1)

        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)

        self.setLayout(layout)

    def set_axis(self,axis, axisname):
        self._scan_label.setText("%s-scan: from %s to %s, max image will be calculated relative to set peak position" % (axisname, axis[0], axis[-1]))

    def _verify_ranges(self):
        if self.roi_from.value() > self.roi_to.value():
            raise ValueError("Invalid input of max img range: from %s > to %s" % (self.roi_from.value(),self.roi_to.value()) )
        if self.roi_from_2.value() > self.roi_to_2.value():
            raise ValueError("Invalid input of max img range: from %s > to %s" % (self.roi_from.value(),self.roi_to.value()) )
        return True

    def get_roi(self):
        self._verify_ranges()
        ddict = {
            'from' : self.roi_from.value(),
            'to' : self.roi_to.value(),
            'vsize' : self.vsize.value(),
            'hsize' : self.hsize.value(),
            'from_2' : self.roi_from_2.value(),
            'to_2' : self.roi_to_2.value(),
            'vsize_2' : self.vsize_2.value(),
            'hsize_2' : self.hsize_2.value()
        }
        return ddict


    def onOk(self):
        try:
            self._verify_ranges()
        except Exception as e:
            qt.QMessageBox.warning(self, "Invalid input", str(e))
            return
        self.accept()

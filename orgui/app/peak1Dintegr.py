# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
#
# class CurvesROIWidget:
# Copyright (c) 2004-2024 European Synchrotron Radiation Facility
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
__credits__ = []
__copyright__ = "Copyright 2020-2026 Timo Fuchs"
__license__ = "MIT License"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import dataclasses
import json
import logging
import sys
import os
from silx.gui import qt
import warnings
from silx.gui import icons

import weakref

# from IPython import embed
import silx.gui.plot


import silx
from silx.gui.plot.actions import control as control_actions

# from silx.gui.widgets.TableWidget import TableWidget
from silx.gui.plot.CurvesROIWidget import ROITable, ROI, _FloatItem, _RoiMarkerHandler
from silx.utils.weakref import WeakMethodProxy

import traceback

from . import qutils
from .config_data import (
    CURVE_CORRECTIONS_GROUP,
    CURVE_CORRECTIONS_SCHEMA_VERSION,
    ConfigData,
    CorrectionState,
    corrections_from_nxdict,
    curve_correction_record_from_nxdict,
    detector_from_nxdict,
)
from .integration_corrections import (
    FOOTPRINT_APPLY,
    FOOTPRINT_KEEP,
    FOOTPRINT_REMOVE,
    corrected_curve_from_record,
    framewise_illumination_divisor,
)
from .. import resources
from .. import logger_utils
from ..datautils.xrayutils.corrections import beamprofile
from ..datautils.xrayutils.corrections import (
    acceptance as acceptance_corrections,
    activearea as activearea_corrections,
    detector as detector_corrections,
    measurement as measurement_corrections,
    normalization as normalization_corrections,
)

import numpy as np
from scipy import interpolate as interp

from contextlib import contextmanager
from functools import partial
from typing import NamedTuple

from silx.io import dictdump
from silx.io.dictdump import h5todict

logger = logging.getLogger(__name__)


@contextmanager
def blockSignals(qobjects):
    try:
        for obj in qobjects:
            obj.blockSignals(True)
        yield
        for obj in qobjects:
            obj.blockSignals(False)
    except TypeError:
        qobjects.blockSignals(True)
        yield
        qobjects.blockSignals(False)


QTVERSION = qt.qVersion()
DEBUG = 0

MAX_ROIS_DISPLAY = 100

if hasattr(np, "trapezoid"):  # ToDo remove for orGUI release >1.5
    _trapz_impl = np.trapezoid  # numpy >= 2.0
else:
    _trapz_impl = np.trapz  # noqa: NPY201  # numpy < 2.0


def _compute_rocking_integration(
    s_array,
    axis,
    croibg_curves,
    croibg_errors_curves,
    roi_info,
    aux,
    use_lorentz,
    use_footprint,
    C_Lor=1.0,
    C_rod=1.0,
    C_flux_on_sample=1.0,
    C_illum_area=1.0,
    C_norm=1.0,
    detector_acceptance=None,
    solid_angle_mean=None,
    ctr_croibg_curves=None,
    ctr_croibg_errors_curves=None,
    angle_unit="deg",
    progress_callback=None,
    should_cancel=None,
):
    """Aggregate rocking-scan ROI counts into signal, background, and F2_hkl.

    Pure-numpy core of :meth:`RockingPeakIntegrator.integrate`, factored out
    so it can be unit tested without a running Qt application or an on-disk
    database. This function does not read or write HDF5 and does not touch
    GUI state; ``RockingPeakIntegrator.integrate`` is responsible for
    resolving inputs from the database and writing the result back.

    :param numpy.ndarray s_array:
        Rocking-scan parameter values, shape ``(n_s,)``.
    :param numpy.ndarray axis:
        Rocking axis values (deg), shape ``(n_pts,)``, shared by every
        ``s`` point.
    :param numpy.ndarray croibg_curves:
        Background-subtracted per-image ROI curve for every ``s`` point,
        shape ``(n_s, n_pts)``.
    :param numpy.ndarray croibg_errors_curves:
        1-sigma errors of ``croibg_curves``, same shape.
    :param dict roi_info:
        Mapping of ROI name (``sig_*``/``bg_*``) to a dict with
        ``from``/``to`` arrays of shape ``(n_s,)``, in the same units as
        ``axis``.
    :param dict aux:
        Mapping of auxiliary counter name to an array of shape ``(n_pts,)``.
    :param bool use_lorentz:
        Apply the Lorentz/rod-intersection correction and compute
        ``F2_hkl``/``F2_hkl_errors``.
    :param bool use_footprint:
        Apply the numerical active-area correction to the corrected
        intensities.
    :param C_Lor:
        Scalar ``1.0`` or array of shape ``(n_s, n_pts)``, Lorentz factor.
    :param C_rod:
        Scalar ``1.0`` or array of shape ``(n_s, n_pts)``, rod-intersection
        factor.
    :param C_flux_on_sample:
        Scalar ``1.0`` or array of shape ``(n_s, n_pts)``, intercepted-flux
        diagnostic used to construct ``C_illum_area``. It is stored but is
        not a second divisor.
    :param C_illum_area:
        Scalar ``1.0`` or array of shape ``(n_s, n_pts)``, illuminated-area
        factor.
    :param C_norm:
        Scalar ``1.0`` or array of shape ``(n_s, n_pts)``, the per-frame
        exposure-time and monitor divisor. It is applied **inside** the
        rocking integral, not to the result: with a varying counting time or a
        drifting monitor the quantity Vlieg's expression integrates is
        :math:`\\int N(\\omega)/(T M)\\,d\\omega`, and dividing the finished
        integral by a mean would only be equivalent for constant counters.
    :param detector_acceptance:
        Out-of-plane acceptance :math:`\\Delta\\gamma` of the region of
        interest per ``s`` point, in **radian**, shape ``(n_s,)``. A rocking
        scan intercepts a slice of rod proportional to it, so it divides
        ``F2_hkl``. ``None`` leaves it out, which reproduces the historical,
        acceptance-blind scale.
    :param solid_angle_mean:
        Region mean of :math:`1/\\widetilde{\\Omega}` per ``s`` point, shape
        ``(n_s,)``, when the solid-angle correction is already inside the
        curves. It divides ``F2_hkl``, removing it again: a region sum is the
        complete angular integral already, so the correction double-counts the
        detector obliquity in a structure factor even though it is what a
        broad, non-rod feature wants on its intensity. ``None`` when the
        correction was not applied.
    :param ctr_croibg_curves:
        Optional polarization-only per-image ROI curves for new extractions,
        with the same shape as ``croibg_curves``. When supplied, ``F2_hkl``
        is derived from this branch while the saved intensity outputs retain
        ``croibg_curves``. Legacy curves omit it and keep the historical
        solid-angle compensation path.
    :param ctr_croibg_errors_curves:
        One-sigma errors paired with ``ctr_croibg_curves``. Both CTR arrays
        must be supplied together.
    :param str angle_unit:
        Unit of ``axis``, ``'deg'`` or ``'rad'``. The published expressions
        integrate the rocking angle in radian; ``'deg'`` converts the integral
        when ``F2_hkl`` is formed, leaving the stored intensities and interval
        widths in the unit they were measured in.
    :param progress_callback:
        Optional callable invoked with the current ``s`` index after each
        point is processed.
    :param should_cancel:
        Optional zero-argument callable; if it returns truthy after a given
        ``s`` point, the remaining points are left unintegrated (matching
        the original in-GUI cancel behavior).
    :returns:
        Mapping with keys ``int_data``, ``croi``, ``croi_errors``,
        ``raw_croi``, ``raw_croi_errors``, ``bgroi``, ``bgroi_errors``,
        ``raw_bgroi``, ``raw_bgroi_errors``, ``croibg``, ``croibg_errors``,
        ``raw_croibg``, ``raw_croibg_errors``, ``auxil``, and, only when
        ``use_lorentz`` is ``True``, ``F2_hkl`` and ``F2_hkl_errors``.
    :rtype: dict
    """
    if (ctr_croibg_curves is None) != (ctr_croibg_errors_curves is None):
        raise ValueError(
            "ctr_croibg_curves and ctr_croibg_errors_curves must be given together"
        )
    int_data = {}
    for roikey in roi_info:
        if roikey.startswith("sig") or roikey.startswith("bg"):
            int_data[roikey] = {
                "cnts": [],
                "cnts_errors": [],
                "raw_cnts": [],
                "raw_cnts_errors": [],
                "int_interval": [],
                "C_norm": [],
                "C_Lor": [],
                "C_rod": [],
                "C_Lorentz_rod": [],
                "C_flux_on_sample": [],
                "C_illum_area": [],
                "auxillary": dict((a, []) for a in aux),
                "auxillary_int": dict((a, []) for a in aux),
                "auxillary_num": dict((a, []) for a in aux),
            }

    for i, s in enumerate(s_array):
        croibg = croibg_curves[i]
        croibg_errors = croibg_errors_curves[i]

        for roikey in int_data:
            roi = roi_info[roikey]

            idx_from = np.argmin(np.abs(axis - roi["from"][i]))
            idx_to = np.argmin(np.abs(axis - roi["to"][i]))
            if idx_from > idx_to:
                idx_from, idx_to = idx_to, idx_from
            int_interval = np.abs(axis[idx_to] - axis[idx_from])  # can be negative!
            sign_interval = np.sign(
                axis[idx_to] - axis[idx_from]
            )  # we force integrals positive
            int_data[roikey]["int_interval"].append(int_interval)

            # ROI boundaries are included so the integrated samples, error
            # weights, and reported interval all describe the same domain.
            roi_slice = slice(idx_from, idx_to + 1)
            roi_axis = axis[roi_slice]
            cnts = croibg[roi_slice]
            cnts_errors = croibg_errors[roi_slice]

            C_corr = np.ones(cnts.size, dtype=float)
            if not np.isscalar(C_norm) or C_norm != 1.0:
                # Per-frame, so it belongs under the integral sign; see the
                # C_norm parameter documentation.
                C_norm_roi = np.broadcast_to(
                    np.asarray(C_norm, dtype=float), croibg_curves.shape
                )[i][roi_slice]
                int_data[roikey]["C_norm"].append(np.mean(C_norm_roi))
                C_corr = C_corr * C_norm_roi
            else:
                int_data[roikey]["C_norm"].append(1.0)

            if use_lorentz:
                lorentz_roi = np.asarray(C_Lor[i][roi_slice], dtype=float)
                rod_roi = np.asarray(C_rod[i][roi_slice], dtype=float)

                def interval_mean(values):
                    if values.size < 2 or int_interval == 0.0:
                        return float(np.mean(values))
                    return float(
                        _trapz_impl(values, roi_axis)
                        * sign_interval
                        / int_interval
                    )

                int_data[roikey]["C_Lor"].append(interval_mean(lorentz_roi))
                int_data[roikey]["C_rod"].append(interval_mean(rod_roi))
                int_data[roikey]["C_Lorentz_rod"].append(
                    interval_mean(lorentz_roi * rod_roi)
                )

            if use_footprint:
                int_data[roikey]["C_flux_on_sample"].append(
                    np.mean(C_flux_on_sample[i][roi_slice])
                )
                int_data[roikey]["C_illum_area"].append(
                    np.mean(C_illum_area[i][roi_slice])
                )
                # C_illum_area is the numerical Vlieg active area and already
                # contains the beam/sample overlap integral reported as
                # C_flux_on_sample. Multiplying both would count that overlap
                # twice and overcorrect a beam that clips the sample edge.
                C_corr *= C_illum_area[i][roi_slice]

            I_raw = (
                _trapz_impl(cnts, roi_axis) * sign_interval
            )  # we force integrals positive
            I_corr = (
                _trapz_impl(cnts / C_corr, roi_axis) * sign_interval
            )  # we force integrals positive

            # The trapezoidal weights must be calculated over the local ROI,
            # since each ROI has its own endpoints.
            roi_dx = np.diff(roi_axis)
            deltaaxis = np.zeros_like(roi_axis)
            if roi_dx.size:
                deltaaxis[0] = roi_dx[0] / 2
                deltaaxis[-1] = roi_dx[-1] / 2
                deltaaxis[1:-1] = (roi_dx[:-1] + roi_dx[1:]) / 2
            I_raw_error = np.sqrt(np.sum((cnts_errors * deltaaxis) ** 2))
            # Propagate directly instead of dividing by I_raw, which can be
            # zero for an empty or cancelling ROI.
            I_corr_error = np.sqrt(np.sum(((cnts_errors / C_corr) * deltaaxis) ** 2))

            int_data[roikey]["raw_cnts"].append(I_raw)
            int_data[roikey]["raw_cnts_errors"].append(I_raw_error)

            int_data[roikey]["cnts"].append(I_corr)
            int_data[roikey]["cnts_errors"].append(I_corr_error)

            for a in aux:
                int_data[roikey]["auxillary_int"][a].append(
                    _trapz_impl(aux[a][roi_slice], roi_axis) * sign_interval
                )  # we force integrals positive
                int_data[roikey]["auxillary"][a].append(np.sum(aux[a][roi_slice]))
                int_data[roikey]["auxillary_num"][a].append(
                    float(aux[a][roi_slice].size)
                )

        if progress_callback is not None:
            progress_callback(i)
        if should_cancel is not None and should_cancel():
            break

    for roikey in int_data:
        for d in list(int_data[roikey].keys()):
            if not d.startswith("auxillary"):
                int_data[roikey][d] = np.array(int_data[roikey][d])
            else:
                for dd in list(int_data[roikey][d].keys()):
                    int_data[roikey][d][dd] = np.array(int_data[roikey][d][dd])

    # signals:

    croi = np.zeros(s_array.size, dtype=float)
    croi_errors = np.zeros(s_array.size, dtype=float)
    raw_croi = np.zeros(s_array.size, dtype=float)
    raw_croi_errors = np.zeros(s_array.size, dtype=float)
    sig_interval = np.zeros(s_array.size, dtype=float)
    C_Lorentz = np.zeros(s_array.size, dtype=float)
    C_rod_intersect = np.zeros(s_array.size, dtype=float)
    C_Lorentz_rod = np.zeros(s_array.size, dtype=float)
    aux_cnts_integral = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    aux_cnts_integral_mean = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    aux_cnts_sum = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    aux_cnts_mean = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    aux_cnts_num = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)

    for roikey in int_data:
        if roikey.startswith("sig"):
            sig_interval += int_data[roikey]["int_interval"]
            for a in aux_cnts_num:
                aux_cnts_num[a] += int_data[roikey]["auxillary_num"][a]

    for roikey in int_data:
        if roikey.startswith("sig"):
            croi += int_data[roikey]["cnts"]
            croi_errors += int_data[roikey]["cnts_errors"] ** 2
            raw_croi += int_data[roikey]["raw_cnts"]
            raw_croi_errors += int_data[roikey]["raw_cnts_errors"] ** 2
            for a in aux_cnts_sum:
                aux_cnts_sum[a] += int_data[roikey]["auxillary"][a]
                aux_cnts_mean[a] += int_data[roikey]["auxillary"][a] / aux_cnts_num[a]
                aux_cnts_integral[a] += int_data[roikey]["auxillary_int"][a]
                aux_cnts_integral_mean[a] += (
                    int_data[roikey]["auxillary_int"][a] / sig_interval
                )
            if use_lorentz:
                C_Lorentz += int_data[roikey]["C_Lor"] * (
                    int_data[roikey]["int_interval"] / sig_interval
                )
                C_rod_intersect += int_data[roikey]["C_rod"] * (
                    int_data[roikey]["int_interval"] / sig_interval
                )
                C_Lorentz_rod += int_data[roikey]["C_Lorentz_rod"] * (
                    int_data[roikey]["int_interval"] / sig_interval
                )

    raw_croi_errors = np.sqrt(raw_croi_errors)
    croi_errors = np.sqrt(croi_errors)

    bgroi = np.zeros(s_array.size, dtype=float)
    bgroi_errors = np.zeros(s_array.size, dtype=float)
    raw_bgroi = np.zeros(s_array.size, dtype=float)
    raw_bgroi_errors = np.zeros(s_array.size, dtype=float)
    bg_interval = np.zeros(s_array.size, dtype=float)
    bgaux_cnts_integral = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    bgaux_cnts_integral_mean = dict(
        (a, np.zeros(s_array.size, dtype=float)) for a in aux
    )
    bgaux_cnts_sum = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    bgaux_cnts_mean = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)
    bgaux_cnts_num = dict((a, np.zeros(s_array.size, dtype=float)) for a in aux)

    for roikey in int_data:
        if roikey.startswith("bg"):
            bg_interval += int_data[roikey]["int_interval"]
            for a in bgaux_cnts_num:
                bgaux_cnts_num[a] += int_data[roikey]["auxillary_num"][a]

    for roikey in int_data:
        if roikey.startswith("bg"):
            ratio = sig_interval / bg_interval
            bgroi += int_data[roikey]["cnts"] * ratio
            bgroi_errors += (int_data[roikey]["cnts_errors"] * ratio) ** 2
            raw_bgroi += int_data[roikey]["raw_cnts"] * ratio
            raw_bgroi_errors += (int_data[roikey]["raw_cnts_errors"] * ratio) ** 2
            for a in bgaux_cnts_sum:
                bgaux_cnts_sum[a] += int_data[roikey]["auxillary"][a]
                bgaux_cnts_mean[a] += (
                    int_data[roikey]["auxillary"][a] / bgaux_cnts_num[a]
                )
                bgaux_cnts_integral[a] += int_data[roikey]["auxillary_int"][a]
                bgaux_cnts_integral_mean[a] += (
                    int_data[roikey]["auxillary_int"][a] / bg_interval
                )

    raw_bgroi_errors = np.sqrt(raw_bgroi_errors)
    bgroi_errors = np.sqrt(bgroi_errors)

    # not divided by sig_interval - see issue #25
    croibg = croi - bgroi  # / sig_interval # already normalized bg
    croibg_errors = np.sqrt(croi_errors**2 + bgroi_errors**2)  # / sig_interval

    raw_croibg = raw_croi - raw_bgroi  # already normalized bg
    raw_croibg_errors = np.sqrt(raw_croi_errors**2 + raw_bgroi_errors**2)

    result = {
        "int_data": int_data,
        "croi": croi,
        "croi_errors": croi_errors,
        "raw_croi": raw_croi,
        "raw_croi_errors": raw_croi_errors,
        "bgroi": bgroi,
        "bgroi_errors": bgroi_errors,
        "raw_bgroi": raw_bgroi,
        "raw_bgroi_errors": raw_bgroi_errors,
        "croibg": croibg,
        "croibg_errors": croibg_errors,
        "raw_croibg": raw_croibg,
        "raw_croibg_errors": raw_croibg_errors,
    }

    auxil = {"@NX_class": "NXcollection"}
    for a in aux_cnts_sum:
        auxil[a] = {
            "@NX_class": "NXcollection",
            "csum": aux_cnts_sum[a],
            "cmean": aux_cnts_mean[a],
            "cintegral": aux_cnts_integral[a],
            "cintegral_mean": aux_cnts_integral_mean[a],
            "bgsum": bgaux_cnts_sum[a],
            "bgmean": bgaux_cnts_mean[a],
            "bgintegral": bgaux_cnts_integral[a],
            "bgintegral_mean": bgaux_cnts_integral_mean[a],
        }
    result["auxil"] = auxil

    if use_lorentz:
        # The rocking angle is the integration variable and must be in radian
        # (Vlieg eq. 42, Drnec eq. 2). The exposure and monitor divisor is
        # already inside croibg, applied per frame above, so only the angle
        # conversion is left for normalized_intensity to do here.
        intensity = measurement_corrections.normalized_intensity(
            croibg, angle_unit=angle_unit
        )
        intensity_errors = measurement_corrections.normalized_intensity(
            croibg_errors, angle_unit=angle_unit
        )
        # The joint factor is integrated with the same trapezoidal angular
        # quadrature as the counts. Multiplying independently averaged
        # Lorentz and rod terms leaves a covariance residual whenever both
        # vary through the rocking window.
        denominator = C_Lorentz_rod
        if detector_acceptance is not None:
            denominator = denominator * np.asarray(detector_acceptance, dtype=float)
        if solid_angle_mean is not None:
            denominator = denominator * np.asarray(solid_angle_mean, dtype=float)
        result["F2_hkl"] = intensity / denominator
        result["F2_hkl_errors"] = intensity_errors / denominator

        if ctr_croibg_curves is not None:
            # New extractions carry a polarization-only curve made from the
            # same base signal. Re-run the pure aggregation on that branch so
            # signal/background windows, nonuniform or reversed angular
            # quadrature and error propagation stay identical. No separately
            # estimated solid-angle mean enters this path.
            ctr_result = _compute_rocking_integration(
                s_array,
                axis,
                ctr_croibg_curves,
                ctr_croibg_errors_curves,
                roi_info,
                aux,
                use_lorentz,
                use_footprint,
                C_Lor=C_Lor,
                C_rod=C_rod,
                C_flux_on_sample=C_flux_on_sample,
                C_illum_area=C_illum_area,
                C_norm=C_norm,
                detector_acceptance=detector_acceptance,
                solid_angle_mean=None,
                angle_unit=angle_unit,
            )
            result["F2_hkl"] = ctr_result["F2_hkl"]
            result["F2_hkl_errors"] = ctr_result["F2_hkl_errors"]

    return result


class RockingPeakIntegrator(qt.QMainWindow):
    def __init__(self, database, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.database = database
        fontMetric = self.fontMetrics()
        iconSize = qt.QSize(fontMetric.height(), fontMetric.height())
        self.filedialogdir = "."
        self._currentRoInfo = {}
        self._idx = 0
        self.footprint_action = FOOTPRINT_KEEP
        self.replacement_illumination = None
        self.replacement_illumination_convention = None
        self.allow_illumination_convention_change = False

        dbdockwidget = qt.QDockWidget("Integrated data")

        self.dbside = qt.QMainWindow()

        self._addRoiDialog = ROICreatorDialog(np.array([0, 0]), "unknown")
        self._addAllRoiDialog = IntegrationEstimator(np.array([0, 0]), "unknown")

        self.integrationCorrection = IntegrationCorrectionsDialog()

        self.dbview = silx.gui.hdf5.Hdf5TreeView()
        self.dbview.setSortingEnabled(True)
        self.dbview.setModel(self.database.hdf5model)
        self.dbview.setExpandsOnDoubleClick(False)
        self.dbview.addContextMenuCallback(self.database.nexus_treeview_callback)

        self.dbside.setCentralWidget(self.dbview)

        dbdockwidget.setWidget(self.dbside)
        self.addDockWidget(qt.Qt.RightDockWidgetArea, dbdockwidget)

        self.plotROIselect = silx.gui.plot.PlotWindow(
            parent=self,
            backend=None,
            resetzoom=True,
            autoScale=True,
            logScale=True,
            grid=True,
            curveStyle=True,
            colormap=False,
            aspectRatio=False,
            yInverted=False,
            copy=True,
            save=True,
            print_=True,
            control=False,
            position=True,
            roi=False,
            mask=False,
            fit=True,
        )

        toolbar = qt.QToolBar()
        toolbar.addAction(
            control_actions.OpenGLAction(parent=toolbar, plot=self.plotROIselect)
        )
        self.plotROIselect.addToolBar(toolbar)

        self.curveSlider = qutils.DataRangeSlider()

        self.roiwidget = CurvesROIWidget(self, "ROIs", self.plotROIselect)
        self.roiwidget.roiTable.clear()

        # self.roiwidget.sigROISignal

        # self.plotROIselect.removeToolBar(self.plotROIselect.toolBar())
        # self.plotROIselect.addToolBar(qt.Qt.LeftToolBarArea ,self.plotROIselect.toolBar())  # noqa: E501

        # self.plotROintegrated = silx.gui.plot.Plot1D(self)

        self.layout = qt.QVBoxLayout()
        self.mainwidget = qt.QWidget()

        self.zoomslider = qt.QSlider()
        self.zoomslider.setOrientation(qt.Qt.Horizontal)
        self.zoomslider.setEnabled(True)
        self.zoomslider.setMinimum(0)
        self.zoomslider.setMaximum(100)
        self.zoomslider.setValue(10)

        self.zoom_menu = qt.QMenu()

        self.zoomwidget = qt.QWidgetAction(self.zoom_menu)
        self.zoomwidget.setDefaultWidget(self.zoomslider)

        self.zoom_menu.addAction(self.zoomwidget)

        # self.alpha_btn = qt.QToolButton(resources.getQicon("sum_image.png"),"slider")
        self.zoom_btn = qt.QToolButton()
        self.zoom_btn.setIcon(resources.getQicon("search"))
        self.zoom_btn.setToolTip("set automatic zoom factor")
        self.zoom_btn.setPopupMode(qt.QToolButton.InstantPopup)
        self.zoom_btn.setMenu(self.zoom_menu)

        plotselecttools_bar = qt.QHBoxLayout()
        plotselecttools_bar.addWidget(self.curveSlider)
        plotselecttools_bar.addWidget(self.zoom_btn)
        self.autozoom_checkbox = qt.QCheckBox("auto zoom")
        self.autozoom_checkbox.setChecked(True)
        plotselecttools_bar.addWidget(self.autozoom_checkbox)

        _plotlayout = qt.QVBoxLayout()
        _plotlayout.addWidget(self.plotROIselect)
        _plotlayout.addLayout(plotselecttools_bar)

        self.layout.addLayout(_plotlayout, 2)
        # self.layout.addWidget(self.curveSlider, 1)

        # ROIS

        bottom_layout = qt.QHBoxLayout()
        bottom_layout.addWidget(self.roiwidget)

        roi_edit_layout = qt.QVBoxLayout()

        anchorROIsGroup = qt.QGroupBox("Anchor ROIs")
        anchorROIsGroupLayout = qt.QHBoxLayout()

        # self.anchorROIButton = qt.QPushButton()
        # self.anchorROIButton.setIcon(resources.getQicon("anchor-ROI"))
        # self.anchorROIButton.setIconSize(iconSize)
        # self.anchorROIButton.setCheckable(True)
        #
        # self.anchorROIButton.toggled.connect(self.onAnchorBtnToggled)

        self.anchorLoadRoiBtn = qt.QPushButton()
        self.anchorLoadRoiBtn.setIcon(icons.getQIcon("document-open"))
        self.anchorLoadRoiBtn.setIconSize(iconSize)
        self.anchorLoadRoiBtn.clicked.connect(self.onAnchorLoadRoi)

        self.anchorSaveRoiBtn = qt.QPushButton()
        self.anchorSaveRoiBtn.setIcon(icons.getQIcon("document-save"))
        self.anchorSaveRoiBtn.setIconSize(iconSize)
        self.anchorSaveRoiBtn.clicked.connect(self.onAnchorSaveRoi)

        self.previousROIButton = qt.QPushButton()
        self.previousROIButton.setIcon(icons.getQIcon("previous"))
        self.previousROIButton.setIconSize(iconSize)
        self.previousROIButton.clicked.connect(self.onToPreviousAnchor)
        self.nextROIButton = qt.QPushButton()
        self.nextROIButton.setIcon(icons.getQIcon("next"))
        self.nextROIButton.setIconSize(iconSize)
        self.nextROIButton.clicked.connect(self.onToNextAnchor)

        self.fitAnchorsButton = qt.QPushButton("Fit between anchors")
        self.fitAnchorsButton.clicked.connect(self.fit_anchors_along_rod)

        # anchorROIsGroupLayout.addWidget(self.anchorROIButton)
        anchorROIsGroupLayout.addWidget(self.fitAnchorsButton)
        anchorROIsGroupLayout.addStretch()
        anchorROIsGroupLayout.addWidget(self.anchorLoadRoiBtn)
        anchorROIsGroupLayout.addWidget(self.anchorSaveRoiBtn)
        anchorROIsGroupLayout.addWidget(self.previousROIButton)
        anchorROIsGroupLayout.addWidget(self.nextROIButton)

        anchorROIsGroup.setLayout(anchorROIsGroupLayout)

        roi_edit_layout.addWidget(anchorROIsGroup)

        modifyROIsGroup = qt.QGroupBox("Modify ROIs")
        modifyROIsGroupLayout = qt.QHBoxLayout()

        self.addAllROIButton = qt.QPushButton("Add All")
        self.addAllROIButton.clicked.connect(self.onAddAllROI)

        # self.estimateROIButton.setIcon(icons.getQIcon("previous"))
        # self.estimateROIButton.setIconSize(iconSize)

        self.addROIButton = qt.QPushButton("Add")
        self.addROIButton.clicked.connect(self.onAddROI)
        # self.addROIButton.setIcon(icons.getQIcon("previous"))
        # self.addROIButton.setIconSize(iconSize)

        self.deleteROIButton = qt.QPushButton("Delete")
        self.deleteROIButton.clicked.connect(self.onDeleteROI)
        # self.deleteROIButton.setIcon(icons.getQIcon("previous"))
        # self.deleteROIButton.setIconSize(iconSize)
        self.deleteAllROIButton = qt.QPushButton("Delete All")
        self.deleteAllROIButton.clicked.connect(self.onDeleteAllROI)

        modifyROIsGroupLayout.addWidget(self.addAllROIButton)
        modifyROIsGroupLayout.addWidget(self.addROIButton)
        modifyROIsGroupLayout.addWidget(self.deleteROIButton)
        modifyROIsGroupLayout.addWidget(self.deleteAllROIButton)

        modifyROIsGroup.setLayout(modifyROIsGroupLayout)

        roi_edit_layout.addWidget(modifyROIsGroup)

        integrateOptionsGroup = qt.QGroupBox("CTR reduction")
        integrateOptionsGroupLayout = qt.QGridLayout()

        self.lorentzButton = qt.QCheckBox("Calculate CTR structure factor")
        self.lorentzButton.setToolTip(
            "Apply the rocking-mode Lorentz, rod-intersection and angular-"
            "acceptance factors. Disable this to save diagnostic intensity "
            "when Q or H is deliberately absent."
        )
        self.normalizationStatus = qt.QLabel(
            "Normalization: Unknown (no rocking scan selected)"
        )
        self.normalizationStatus.setWordWrap(True)
        self.footprintAction = qt.QComboBox()
        self.footprintAction.addItem("Keep stored correction", FOOTPRINT_KEEP)
        self.footprintAction.addItem("Apply current beam settings", FOOTPRINT_APPLY)
        self.footprintAction.addItem(
            "Remove stored correction", FOOTPRINT_REMOVE
        )
        self.footprintAction.setToolTip(
            "Every action starts from the stored base curve. Unsafe actions "
            "are disabled when correction provenance is unknown."
        )
        self.footprintStatus = qt.QLabel(
            "Saved result: Unknown (no rocking scan selected)"
        )
        self.footprintStatus.setWordWrap(True)
        self.reductionPreview = qt.QLabel("Next reduction: no scan selected")
        self.reductionPreview.setWordWrap(True)
        self.reductionDetails = qt.QLabel("")
        self.reductionDetails.setWordWrap(True)
        self.footprintOptionsButton = qt.QPushButton("Current beam settings ...")
        self.footprintOptionsButton.clicked.connect(self._showFootprintOptions)
        self.footprintAction.currentIndexChanged.connect(
            self._onFootprintActionChanged
        )
        self.lorentzButton.toggled.connect(self._onFootprintActionChanged)

        integrateOptionsGroupLayout.addWidget(self.lorentzButton, 0, 0, 1, 2)
        integrateOptionsGroupLayout.addWidget(self.normalizationStatus, 1, 0, 1, 2)
        integrateOptionsGroupLayout.addWidget(qt.QLabel("Footprint action:"), 2, 0)
        integrateOptionsGroupLayout.addWidget(self.footprintAction, 2, 1)
        integrateOptionsGroupLayout.addWidget(self.footprintStatus, 3, 0, 1, 2)
        integrateOptionsGroupLayout.addWidget(self.reductionPreview, 4, 0, 1, 2)
        integrateOptionsGroupLayout.addWidget(self.reductionDetails, 5, 0, 1, 2)
        integrateOptionsGroupLayout.addWidget(
            self.footprintOptionsButton, 6, 0, 1, 2
        )

        integrateOptionsGroup.setLayout(integrateOptionsGroupLayout)

        roi_edit_layout.addWidget(integrateOptionsGroup)

        self.integrateButton = qt.QPushButton("integrate")
        self.integrateButton.clicked.connect(self.onIntegrate)

        roi_edit_layout.addWidget(self.integrateButton)

        bottom_layout.addLayout(roi_edit_layout)

        self.layout.addLayout(bottom_layout, 1)

        # self.zoom_btn_act = qt.QWidgetAction(self)
        # self.zoom_btn_act.setDefaultWidget(self.zoom_btn)

        # self.toolbar.addAction(self.showMaxAct)
        # self.toolbar.addAction(self.showSumAct)
        # self.toolbar.addAction(self.alpha_btn_act)

        # self.layout.addWidget(self.plotROintegrated, 1)

        self.mainwidget.setLayout(self.layout)
        self.setCentralWidget(self.mainwidget)

        toolbar = qt.QToolBar("Database toolbar", self.dbside)

        loadDatabaseAct = toolbar.addAction(
            icons.getQIcon("document-open"), "Open orgui database"
        )
        loadDatabaseAct.triggered.connect(self.database.onOpenDatabase)

        savenewact = toolbar.addAction(
            icons.getQIcon("layer-nx"), "Select orgui database location"
        )
        savenewact.triggered.connect(self.database.onSaveNewDBFile)

        saveact = toolbar.addAction(
            icons.getQIcon("document-save"), "Save orgui database"
        )
        saveact.triggered.connect(self.database.onSaveDBFile)

        self.database.sigChangeRockingScan.connect(self.onChangeRockingScan)

        self.zoomslider.valueChanged.connect(self.resetXZoomScaled)

        self.dbside.addToolBar(toolbar)

        # self.curveSlider.setAxis(np.arange(100) / 10., "deg")

        self.curveSlider.sigValueChanged.connect(self.onSliderValueChanged)
        # self.roiwidget.sigROISignal.connect(lambda d: print(d))

        self.integrationCorrection.settingsChanged.connect(
            self._onFootprintActionChanged
        )
        self._refreshReductionCorrectionStatus()

    def _check_ro_present(self):
        if not self.database.nxfile:
            self._currentRoInfo = {}
            raise OSError("No database available.")
        if self._currentRoInfo and "name" in self._currentRoInfo:
            if self._currentRoInfo["name"] in self.database.nxfile:
                return True
            else:
                n = self._currentRoInfo["name"]
                self._currentRoInfo = {}
                raise ValueError(f"Scan {n} is not in database.")
        else:
            raise ValueError("No rocking scan loaded.")

    # GUI-only: user-triggered non-modal dialog.
    def _showFootprintOptions(self):
        """Show the beam-profile settings without blocking this window.

        .. note::
           GUI-only. The dialog is non-modal so the beam profile can be
           adjusted while the integrated curves stay visible.
        """
        self.integrationCorrection.show()
        self.integrationCorrection.raise_()
        self.integrationCorrection.activateWindow()

    def useSharedFootprintOptions(self, dialog):
        """Use an externally owned beam-profile dialog instead of the own one.

        The incident beam is a property of the experiment, not of one
        integration mode, so the rocking-scan and stationary-scan
        integrations are pointed at a single dialog rather than each keeping
        their own beam profile.

        :param dialog: The
            :class:`IntegrationCorrectionsDialog` to adopt.
        """
        if dialog is None or dialog is self.integrationCorrection:
            return
        previous = self.integrationCorrection
        self.integrationCorrection = dialog
        self.integrationCorrection.settingsChanged.connect(
            self._onFootprintActionChanged
        )
        if previous is not None:
            previous.deleteLater()
        self._onFootprintActionChanged()

    @staticmethod
    def _correctionUiState(record, record_present=True):
        """Describe safe rocking-reduction actions for a stored record.

        :param CurveCorrectionRecord or None record: Parsed stored record.
        :param bool record_present: Whether a correction group exists at all.
        :rtype: dict
        """
        if record is None:
            reason = (
                "an incomplete correction record"
                if record_present
                else "legacy data without a correction record"
            )
            return {
                "actions": (FOOTPRINT_KEEP,),
                "normalization": "Unknown (legacy data)",
                "footprint": "Unknown",
                "details": (
                    f"Correction provenance is unknown because this scan has {reason}. "
                    "Re-extract images to change this correction safely."
                ),
                "known_total_flux": False,
            }

        if record.algorithm.startswith("legacy_rocking_roi_"):
            return {
                "actions": (FOOTPRINT_KEEP, FOOTPRINT_APPLY),
                "normalization": "Applied later by the legacy reducer",
                "footprint": "Not applied to the saved curve",
                "details": (
                    f"Algorithm {record.algorithm}; scale {record.scale_convention}. "
                    "Current beam settings may be applied using the preserved "
                    "legacy density/area convention."
                ),
                "known_total_flux": False,
            }

        if record.algorithm.startswith("framewise_ctr_total_flux_"):
            reversible = (
                record.base_croibg is not None
                and record.base_croibg_variance is not None
            )
            known_illumination = record.illumination_status in (
                "applied",
                "not_applied",
            )
            actions = [FOOTPRINT_KEEP]
            if reversible and known_illumination:
                actions.append(FOOTPRINT_APPLY)
            if (
                reversible
                and record.illumination_status == "applied"
                and record.illumination_divisor is not None
            ):
                actions.append(FOOTPRINT_REMOVE)
            normalization = {
                "applied": "Applied during extraction",
                "not_applied": "Not applied",
                "unavailable": "Unavailable",
            }.get(record.normalization_status, "Unknown")
            footprint = {
                "applied": "Applied during extraction",
                "not_applied": "Not applied",
                "unavailable": "Unavailable",
            }.get(record.illumination_status, "Unknown")
            extra = ""
            if not reversible or not known_illumination:
                extra = " Re-extract images to change this correction safely."
            return {
                "actions": tuple(actions),
                "normalization": normalization,
                "footprint": footprint,
                "details": (
                    f"Algorithm {record.algorithm}; scale {record.scale_convention}; "
                    f"illumination convention "
                    f"{record.illumination_convention or 'none'}.{extra}"
                ),
                "known_total_flux": True,
            }

        return {
            "actions": (FOOTPRINT_KEEP,),
            "normalization": "Unknown",
            "footprint": "Unknown",
            "details": (
                f"Unsupported correction algorithm {record.algorithm!r}. "
                "Re-extract images to change this correction safely."
            ),
            "known_total_flux": False,
        }

    @staticmethod
    def _storedCurveCorrectionRecord(h5_obj):
        """Return a parsed stored record, or ``None`` when it is incomplete.

        Parsed from the live group: 2-D per-curve arrays stay unread
        ``h5py.Dataset`` objects, valid while the database file is open.
        """
        if CURVE_CORRECTIONS_GROUP not in h5_obj:
            return None
        try:
            return curve_correction_record_from_nxdict(
                h5_obj[CURVE_CORRECTIONS_GROUP]
            )
        except (KeyError, TypeError, ValueError):
            return None

    def _setFootprintActionEnabled(self, action, enabled):
        """Enable one footprint-action choice in the combo-box model."""
        index = self.footprintAction.findData(action)
        if index >= 0:
            self.footprintAction.model().item(index).setEnabled(bool(enabled))

    def _refreshReductionCorrectionStatus(self):
        """Refresh saved-result provenance and safe next actions."""
        record = None
        present = False
        if self._currentRoInfo and self.database.nxfile:
            name = self._currentRoInfo.get("name")
            if name in self.database.nxfile:
                h5_obj = self.database.nxfile[name]
                present = CURVE_CORRECTIONS_GROUP in h5_obj
                record = self._storedCurveCorrectionRecord(h5_obj)
        state = self._correctionUiState(record, present)
        self._activeCorrectionRecord = record
        self._activeCorrectionUiState = state
        allowed = state["actions"]
        for action in (FOOTPRINT_KEEP, FOOTPRINT_APPLY, FOOTPRINT_REMOVE):
            self._setFootprintActionEnabled(action, action in allowed)
        if self.footprintAction.currentData() not in allowed:
            self.footprintAction.setCurrentIndex(
                self.footprintAction.findData(FOOTPRINT_KEEP)
            )
        self.normalizationStatus.setText(
            f"Saved normalization: {state['normalization']}"
        )
        self.footprintStatus.setText(
            f"Saved footprint: {state['footprint']}"
        )
        self.reductionDetails.setText(state["details"])
        self._onFootprintActionChanged()

    def _onFootprintActionChanged(self, *args):
        """Prepare and describe the selected next-reduction footprint action."""
        action = self.footprintAction.currentData() or FOOTPRINT_KEEP
        self.footprint_action = action
        self.replacement_illumination = None
        self.replacement_illumination_convention = None
        self.allow_illumination_convention_change = False

        action_text = {
            FOOTPRINT_KEEP: "keep the stored footprint state",
            FOOTPRINT_APPLY: "apply current beam settings from the base curve",
            FOOTPRINT_REMOVE: "remove the stored footprint from the base curve",
        }[action]
        record = getattr(self, "_activeCorrectionRecord", None)
        known_total_flux = bool(
            getattr(self, "_activeCorrectionUiState", {}).get(
                "known_total_flux", False
            )
        )
        if action == FOOTPRINT_APPLY:
            action_text = (
                "Replaced here using current beam settings"
                if record is not None and record.illumination_status == "applied"
                else "Applied here using current beam settings"
            )
        f2_ready = action != FOOTPRINT_REMOVE
        if known_total_flux:
            f2_ready = record.normalization_status == "applied"
            if action == FOOTPRINT_KEEP:
                f2_ready = f2_ready and record.illumination_status == "applied"
            elif action == FOOTPRINT_APPLY:
                f2_ready = (
                    f2_ready
                    and self.integrationCorrection.horizontalInterceptedFraction()
                    is not None
                )
            else:
                f2_ready = False
        if not f2_ready and self.lorentzButton.isChecked():
            with qt.QSignalBlocker(self.lorentzButton):
                self.lorentzButton.setChecked(False)
        self.lorentzButton.setEnabled(f2_ready)
        quantity = (
            "CTR structure factor"
            if self.lorentzButton.isChecked() and f2_ready
            else "diagnostic intensity (not F²)"
        )
        warning = ""
        if action == FOOTPRINT_REMOVE and self.lorentzButton.isChecked():
            warning = " Remove leaves H absent, so F² is unavailable."
        self.reductionPreview.setText(
            f"Next reduction: {action_text}; output {quantity}.{warning}"
        )

    def _prepareFootprintAction(self, h5_obj):
        """Resolve any live illumination divisor before reading the curve."""
        action = self.footprintAction.currentData() or FOOTPRINT_KEEP
        self.footprint_action = action
        self.replacement_illumination = None
        self.replacement_illumination_convention = None
        record = self._storedCurveCorrectionRecord(h5_obj)
        if (
            action == FOOTPRINT_APPLY
            and record is not None
            and record.algorithm.startswith("framewise_ctr_total_flux_")
        ):
            horizontal = self.integrationCorrection.horizontalInterceptedFraction()
            if horizontal is None:
                raise ValueError(
                    "Applying current total-flux beam settings requires an "
                    "explicit horizontal interception choice."
                )
            if record.alpha is None:
                raise ValueError(
                    "The stored curve has no incidence angles; re-extract "
                    "images before changing its footprint."
                )
            self.replacement_illumination = framewise_illumination_divisor(
                np.asarray(record.alpha),
                self.integrationCorrection.sampleLength(),
                self.integrationCorrection.beamProfile(),
                horizontal_fraction=horizontal,
            )[0]
            self.replacement_illumination_convention = "total_flux_H"
        elif action != FOOTPRINT_KEEP and record is None:
            raise ValueError(
                "Correction provenance is unknown. Re-extract images before "
                "changing the footprint."
            )
        return record

    def onIntegrate(self):
        try:
            self._check_ro_present()
        except Exception as e:
            logger.exception(
                "Cannot integrate scan",
                extra={
                    "title": "Cannot integrate scan",
                    "show_dialog": True,
                    "description": str(e),
                    "parent": self,
                    "dialog_level": logging.WARNING,
                },
            )
            return
        try:
            self.integrate()
        except ValueError as e:
            if e.args[0] == "No rocking scan selected.":
                logger.error("No rocking scan selected.", exc_info=True)
            else:
                logger.exception(
                    "Error during scan integration",
                    extra={
                        "title": "Error during scan integration",
                        "show_dialog": True,
                        "description": str(e),
                        "parent": self,
                    },
                )
        except Exception as e:
            logger.exception(
                "Error during scan integration",
                extra={
                    "title": "Error during scan integration",
                    "show_dialog": True,
                    "description": str(e),
                    "parent": self,
                },
            )
            return

    def set_roscan(self, name):
        ro_info = self.get_rocking_scan_info(name)
        # if ro_info['name'] + '/integration' not in self.database.nxfile:
        #    roi1D_info = self.estimate_roi1D_info(ro_info)
        #    self.database.add_nxdict(roi1D_info, update_mode='modify', h5path=ro_info['name'] + '/integration')  # noqa: E501

        self._currentRoInfo = ro_info
        self._idx = 0

        self.curveSlider.setAxis(self._currentRoInfo["s"], "s")
        self.curveSlider.setIndex(0)
        self.plotRoCurve(0)
        self.plotROIselect.resetZoom()
        if self.autozoom_checkbox.isChecked():
            self.resetXZoomScaled(self.zoomslider.value())
        self._refreshReductionCorrectionStatus()

    def onAnchorSaveRoi(self):
        # GUI-only: user-triggered save dialog path.
        try:
            self._check_ro_present()
            self.saveAnchorROI(None)
        except (ValueError, OSError):
            logger.error(
                "Cannot save roi locations",
                extra={
                    "title": "Cannot save roi locations",
                    "show_dialog": True,
                    "description": "Cannot save ROI locations:\nNo ROI information availabe",  # noqa: E501
                    "parent": self,
                    "dialog_level": logging.WARNING,
                },
            )
            return

        fileTypeDictSave1D = {
            "Plain ascii file (*.dat)": "dat",
            "CSV file (*.csv)": "csv",
            "NumPy format (*.npy)": "ndarray",
        }

        fileTypeDict = {**fileTypeDictSave1D}

        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"

        filename, filetype = qt.QFileDialog.getSaveFileName(
            self, "Save ROI to file", self.filedialogdir, fileTypeFilter[:-2]
        )
        if filename == "":
            return
        self.filedialogdir = os.path.splitext(filename)[0]

        if filetype in fileTypeDictSave1D:
            fileext = fileTypeDictSave1D[filetype]
            try:
                self.saveAnchorROI(filename, fileext)
            except Exception as e:
                logger.exception(
                    "Cannot save ROIs",
                    extra={
                        "title": "Cannot save ROIs",
                        "show_dialog": True,
                        "description": f"Cannot save ROI locations: {e}",
                        "detailed_text": traceback.format_exc(),
                        "parent": self,
                    },
                )

    def saveAnchorROI(self, filename=None, fileext=None):
        """Save anchor ROI information without invoking GUI code.

        :param str filename:
            Destination path. If ``None``, only validate that anchor ROI data is
            available for saving.
        :param str fileext:
            Output format identifier such as ``"dat"``, ``"csv"``, or
            ``"ndarray"``. If omitted, it is inferred from ``filename`` when a
            filename is provided.
        :raises ValueError:
            If no rocking scan is loaded or no anchor ROI information is
            available.
        :raises Exception:
            If the requested file format is unsupported.
        """
        self._check_ro_present()
        integration_path = self._currentRoInfo["name"] + "/integration/"
        if integration_path not in self.database.nxfile:
            raise ValueError("No ROI information availabe")
        if integration_path + "peakpos" not in self.database.nxfile:
            raise ValueError("No ROI information availabe")

        if filename is None:
            return

        if fileext is None:
            _, inferred_ext = os.path.splitext(filename)
            if inferred_ext == ".npy":
                fileext = "ndarray"
            else:
                fileext = inferred_ext.lstrip(".").lower()

        self.saveROI(filename, fileext)

    def saveROI(self, filename, fileext="dat"):
        roi_info = h5todict(
            self.database.nxfile, self._currentRoInfo["name"] + "/integration/"
        )

        fmt = "%.7g"
        csvdelim = ";"

        data = []
        header = ["peakpos"]
        data.append(roi_info["peakpos"])

        for k in list(roi_info.keys()):
            if k.startswith("sig") or k.startswith("bg"):
                header.append(k + "_from")
                data.append(roi_info[k]["from"])
                header.append(k + "_to")
                data.append(roi_info[k]["to"])
                header.append(k + "_anchor")
                data.append(roi_info[k]["anchor"].astype(float))

        data = np.vstack(data).T

        if fileext == "dat":
            header_line = " ".join(header)
            np.savetxt(filename, data, header=header_line, fmt=fmt)
        elif fileext == "csv":
            header_line = csvdelim.join(header)
            np.savetxt(filename, data, header=header_line, fmt=fmt, delimiter=csvdelim)
        elif fileext == "ndarray":
            np.save(filename, {"header": header, "data": data})
        else:
            raise Exception(f"No supported file type {fileext}")

    # def onAnchorLoadRoi(self):
    # datasetDialog = DataFileDialog.DataFileDialog(self)
    # datasetDialog.setFilterMode(DataFileDialog.DataFileDialog.FilterMode.ExistingDataset)  # noqa: E501
    #
    # def customFilter(obj):
    #    print(obj.basename)
    #    return True
    #
    #    if "NX_class" in obj.attrs:
    #        if 'orgui_meta' in obj.attrs and obj.attrs['orgui_meta'] == 'rocking':
    #            if 'integration' in obj:
    #                return True
    #    return False
    #
    # if datasetDialog.exec():
    #    print(datasetDialog.selectedUrl())

    def onAnchorLoadRoi(self):
        # GUI-only: user-triggered open dialog path.
        try:
            self._check_ro_present()
        except Exception as e:
            logger.exception(
                "Cannot load roi locations",
                extra={
                    "title": "Cannot load roi locations",
                    "show_dialog": True,
                    "description": f"Cannot load ROI locations:\n{e}",
                    "parent": self,
                    "dialog_level": logging.WARNING,
                },
            )
            return
        if self._currentRoInfo:
            fileTypeDictSave1D = {
                "Plain ascii file (*.dat)": "dat",
                "CSV file (*.csv)": "csv",
                "NumPy format (*.npy)": "ndarray",
            }
            fileTypeDict = {**fileTypeDictSave1D}

            fileTypeFilter = ""
            for f in fileTypeDict:
                fileTypeFilter += f + ";;"

            filename, filetype = qt.QFileDialog.getOpenFileName(
                self, "Load ROI from file", self.filedialogdir, fileTypeFilter[:-2]
            )
            if filename == "":
                return
            self.filedialogdir = os.path.splitext(filename)[0]

            try:
                self.loadAnchorROI(filename)
            except Exception as e:
                logger.exception(
                    "Cannot load ROIs",
                    extra={
                        "title": "Cannot load ROIs",
                        "show_dialog": True,
                        "description": f"Cannot apply ROI locations to rocking scan: {e}",  # noqa: E501
                        "detailed_text": traceback.format_exc(),
                        "parent": self,
                    },
                )
                return

    def loadAnchorROI(self, filename, interpolate=True, pktol=0.5):
        """Load anchor ROI information without invoking GUI code.

        :param str filename:
            Source path containing ROI locations.
        :param bool interpolate:
            If ``True``, interpolate ROI limits onto the current rocking scan
            peak positions when needed.
        :param float pktol:
            Maximum tolerated mismatch in peak position before interpolation is
            rejected. Unit follows the rocking scan axis and is typically deg.
        """
        roi_info = self.loadROIinfo(filename)
        self.setROIinfo(roi_info, interpolate=interpolate, pktol=pktol)

    def setROIinfo(self, roi_info, interpolate=True, pktol=0.5):
        if self._currentRoInfo:
            ro_info = self.get_rocking_scan_info(self._currentRoInfo["name"])
            sc_h5 = self.database.nxfile[ro_info["name"]]["rois"]
            if ro_info["axisname"] == "mu":
                pk_pos_exact = sc_h5["alpha_pk"][()]
            elif ro_info["axisname"] == "th":
                pk_pos_exact = sc_h5["theta_pk"][()]
            elif ro_info["axisname"] == "chi":
                pk_pos_exact = sc_h5["chi_pk"][()]
            elif ro_info["axisname"] == "phi":
                pk_pos_exact = sc_h5["phi_pk"][()]
            else:
                raise ValueError(
                    "Cannot estimate peak position: unknown scan axis {}".format(
                        ro_info["axisname"]
                    )  # noqa: E501
                )

            def _set_roi_info(roi_info):
                if (
                    self._currentRoInfo["name"] + "/integration/"
                    in self.database.nxfile
                ):
                    del self.database.nxfile[
                        self._currentRoInfo["name"] + "/integration/"
                    ]

                for k in roi_info:
                    if k.startswith(("sig", "bg")):
                        roi_info[k]["@NX_class"] = "NXcollection"
                        roi_info[k]["from"] = roi_info[k]["from"].astype(np.float64)
                        roi_info[k]["to"] = roi_info[k]["to"].astype(np.float64)
                        roi_info[k]["anchor"] = roi_info[k]["anchor"].astype(bool)

                roi_info["peakpos"] = roi_info["peakpos"].astype(np.float64)
                roi_info["@NX_class"] = "NXcollection"
                roi_info["@info"] = "ROI information for the rocking scan integration"

                self.database.add_nxdict(
                    roi_info,
                    update_mode="modify",
                    h5path=ro_info["name"] + "/integration",
                )
                self.plotRoCurve(self._idx)

            pk_pos = roi_info["peakpos"]

            if pk_pos.size != pk_pos_exact.size:
                if not interpolate:
                    raise ValueError(
                        f"Scan length mismatch: requires interpolation is: {pk_pos_exact.size}, supplied roi_info: {pk_pos.size}"  # noqa: E501
                    )
            elif not np.all(np.isclose(pk_pos, pk_pos_exact)):
                if not interpolate:
                    raise ValueError("Peak position mismatch: requires interpolation")
            else:
                _set_roi_info(roi_info)
                return

            if np.abs(np.amax(pk_pos) - np.amax(pk_pos_exact)) >= pktol:
                raise ValueError(
                    "peak position mismatch tolerance exceeded at max by %s"
                    % (np.amax(pk_pos) - np.amax(pk_pos_exact))
                )

            if np.abs(np.amin(pk_pos) - np.amin(pk_pos_exact)) >= pktol:
                raise ValueError(
                    "peak position mismatch tolerance exceeded at min by %s"
                    % (np.amin(pk_pos) - np.amin(pk_pos_exact))
                )

            roi_info_interp = dict()

            roi_info_interp["peakpos"] = pk_pos_exact
            for k in roi_info:
                if k.startswith(("sig", "bg")):
                    roi_info_interp[k] = dict()
                    interfrom = interp.interp1d(
                        roi_info["peakpos"],
                        roi_info[k]["from"],
                        fill_value="extrapolate",
                    )
                    roi_info_interp[k]["from"] = interfrom(pk_pos_exact)

                    interto = interp.interp1d(
                        roi_info["peakpos"], roi_info[k]["to"], fill_value="extrapolate"
                    )
                    roi_info_interp[k]["to"] = interto(pk_pos_exact)

                    interanchor = interp.interp1d(
                        roi_info["peakpos"],
                        roi_info[k]["anchor"],
                        kind="nearest",
                        fill_value="extrapolate",
                    )
                    roi_info_interp[k]["anchor"] = interanchor(pk_pos_exact).astype(
                        bool
                    )

            _set_roi_info(roi_info_interp)

    def loadROIinfo(self, filename):
        """Load ROI information from a file.

        :param str filename:
            ROI file in ``.dat``, ``.csv``, or ``.npy`` format.
        :returns:
            Mapping with ``peakpos`` and per-ROI ``from``, ``to``, and
            ``anchor`` arrays.
        :rtype: dict
        :raises Exception:
            If the file extension is unsupported or the file contents do not
            provide the required ROI columns.
        """
        fname, fileext = os.path.splitext(filename)
        if fileext == ".dat":
            data = np.genfromtxt(filename, names=True)
            roi_info = self._load_roi_info_from_structured_array(data)
        elif fileext == ".csv":
            data = np.genfromtxt(filename, names=True, delimiter=";")
            roi_info = self._load_roi_info_from_structured_array(data)
        elif fileext == ".npy":
            roi_info = self._load_roi_info_from_npy(filename)
        else:
            raise Exception(f"Not supported file type {fileext}")
        return roi_info

    def _load_roi_info_from_npy(self, filename):
        """Load ROI information from the self-describing NumPy format.

        :param str filename:
            NumPy ``.npy`` file created by :meth:`saveROI` with
            ``fileext="ndarray"``.
        :returns:
            ROI information mapping used by :meth:`setROIinfo`.
        :rtype: dict
        :raises ValueError:
            If the file does not contain the expected header and data payload.
        """
        payload = np.load(filename, allow_pickle=True)
        if not isinstance(payload, np.ndarray) or payload.shape != ():
            raise ValueError(f"NumPy ROI file {filename} has unsupported payload shape")

        payload = payload.item()
        if not isinstance(payload, dict):
            raise ValueError(
                f"NumPy ROI file {filename} does not contain a ROI payload dictionary"
            )
        if "header" not in payload or "data" not in payload:
            raise ValueError(f"NumPy ROI file {filename} is missing header or data")

        header = tuple(payload["header"])
        data = np.asarray(payload["data"])
        if data.ndim != 2:
            raise ValueError(f"NumPy ROI file {filename} data must be a 2D array")
        if data.shape[1] != len(header):
            raise ValueError(
                f"NumPy ROI file {filename} column mismatch: data has {data.shape[1]} columns, header has {len(header)}"  # noqa: E501
            )

        structured = np.core.records.fromarrays(data.T, names=header)
        return self._load_roi_info_from_structured_array(structured)

    def _load_roi_info_from_structured_array(self, data):
        """Build ROI info from a structured array.

        :param numpy.ndarray data:
            Structured array with a ``peakpos`` column and ROI columns named
            like ``sig_1_from`` and ``bg_2_anchor``.
        :returns:
            ROI information mapping used by :meth:`setROIinfo`.
        :rtype: dict
        :raises ValueError:
            If required columns are missing.
        """
        if data.dtype.names is None or "peakpos" not in data.dtype.names:
            raise ValueError("ROI file does not contain required column peakpos")

        roi_names = []
        for k in data.dtype.names:
            if k.startswith(("sig", "bg")):
                k_sp = k.split("_")
                if len(k_sp) < 3:
                    continue
                int(k_sp[1])
                rname = k_sp[0] + "_" + k_sp[1]
                if rname not in roi_names:
                    roi_names.append(rname)

        roi_info = dict()
        roi_info["peakpos"] = data["peakpos"]
        for k in roi_names:
            for suffix in ("_from", "_to", "_anchor"):
                if k + suffix not in data.dtype.names:
                    raise ValueError(
                        "ROI file is missing required column %s" % (k + suffix)
                    )
            roi_info[k] = dict()
            roi_info[k]["from"] = data[k + "_from"]
            roi_info[k]["to"] = data[k + "_to"]
            roi_info[k]["anchor"] = data[k + "_anchor"].astype(bool)
        return roi_info

    # ddict = {
    #        "event": "indexChanged",
    #        "oldtxt": self._lineTxt,
    #        "newtxt": txt,
    #        "idx" : idx,
    #        "value" : self._data[idx],
    #        "id": id(self),
    #    }

    def onToNextAnchor(self):
        try:
            self._check_ro_present()
        except Exception:  # non essential: silent
            return
        if (
            self._currentRoInfo
            and self._currentRoInfo["name"] + "/integration/" in self.database.nxfile
        ):
            roi_info = h5todict(
                self.database.nxfile, self._currentRoInfo["name"] + "/integration/"
            )
            anchors_all = []
            for roikey in roi_info:
                if roikey.startswith("sig") or roikey.startswith("bg"):
                    rdict = roi_info[roikey]
                    anchors = np.nonzero(rdict["anchor"])[0]
                    anchors_all.append(anchors)
            anchors_all = np.sort(np.unique(np.concatenate(anchors_all)))

            anchors_next = anchors_all[anchors_all > self._idx]
            if anchors_next.size > 0:
                self.plotRoCurve(anchors_next[0])

    def onToPreviousAnchor(self):
        try:
            self._check_ro_present()
        except Exception:  # non essential: silent
            return
        if (
            self._currentRoInfo
            and self._currentRoInfo["name"] + "/integration/" in self.database.nxfile
        ):
            roi_info = h5todict(
                self.database.nxfile, self._currentRoInfo["name"] + "/integration/"
            )
            anchors_all = []
            for roikey in roi_info:
                if roikey.startswith("sig") or roikey.startswith("bg"):
                    rdict = roi_info[roikey]
                    anchors = np.nonzero(rdict["anchor"])[0]
                    anchors_all.append(anchors)
            anchors_all = np.sort(np.unique(np.concatenate(anchors_all)))

            anchors_prev = anchors_all[anchors_all < self._idx]
            if anchors_prev.size > 0:
                self.plotRoCurve(anchors_prev[-1])

    # def get_anchors(self):

    def fit_anchors_along_rod(self):
        if (
            self._currentRoInfo
            and self._currentRoInfo["name"] + "/integration/" in self.database.nxfile
        ):
            roih5grp = self.database.nxfile[
                self._currentRoInfo["name"] + "/integration/"
            ]
            roi_info = h5todict(
                self.database.nxfile, self._currentRoInfo["name"] + "/integration/"
            )
            pk_pos_exact = roi_info["peakpos"]
            self._currentRoInfo["axis"]
            s_array = self._currentRoInfo["s"]
            for roikey in roi_info:
                if roikey.startswith("sig") or roikey.startswith("bg"):
                    rdict = roi_info[roikey]
                    anchors = np.nonzero(rdict["anchor"])[0]
                    if anchors.size == 0:
                        continue
                    from_ar_anchors = rdict["from"][anchors] - pk_pos_exact[anchors]
                    to_ar_anchors = rdict["to"][anchors] - pk_pos_exact[anchors]

                    from_ar = np.zeros(rdict["from"].size)
                    to_ar = np.zeros(rdict["to"].size)

                    # constant to first anchor
                    idx_prev = 0
                    from_ar[: anchors[idx_prev]] = (
                        pk_pos_exact[: anchors[idx_prev]] + from_ar_anchors[idx_prev]
                    )
                    to_ar[: anchors[idx_prev]] = (
                        pk_pos_exact[: anchors[idx_prev]] + to_ar_anchors[idx_prev]
                    )
                    if anchors.size > 0:
                        for idx in range(1, anchors.size):
                            slope = (
                                from_ar_anchors[idx] - from_ar_anchors[idx_prev]
                            ) / (s_array[anchors[idx]] - s_array[anchors[idx_prev]])
                            from_tmp = slope * (
                                s_array[anchors[idx_prev] : anchors[idx]]
                                - s_array[anchors[idx_prev]]
                            )
                            from_ar[anchors[idx_prev] : anchors[idx]] = (
                                from_tmp + from_ar_anchors[idx_prev]
                            ) + pk_pos_exact[anchors[idx_prev] : anchors[idx]]

                            slope = (to_ar_anchors[idx] - to_ar_anchors[idx_prev]) / (
                                s_array[anchors[idx]] - s_array[anchors[idx_prev]]
                            )
                            to_tmp = slope * (
                                s_array[anchors[idx_prev] : anchors[idx]]
                                - s_array[anchors[idx_prev]]
                            )
                            to_ar[anchors[idx_prev] : anchors[idx]] = (
                                to_tmp + to_ar_anchors[idx_prev]
                            ) + pk_pos_exact[anchors[idx_prev] : anchors[idx]]

                            idx_prev = idx

                    # constant from last anchor
                    from_ar[anchors[idx_prev] :] = (
                        pk_pos_exact[anchors[idx_prev] :] + from_ar_anchors[idx_prev]
                    )
                    to_ar[anchors[idx_prev] :] = (
                        pk_pos_exact[anchors[idx_prev] :] + to_ar_anchors[idx_prev]
                    )

                    roih5grp[roikey]["from"][:] = from_ar
                    roih5grp[roikey]["to"][:] = to_ar
            self.plotRoCurve(self._idx)

    def _rocking_normalization(self, aux, size):
        """Per-frame exposure and monitor divisor of the rocking scan.

        Unlike the stationary path this cannot ask a live scan object: a
        rocking integration runs off the database, so the counters are read
        from the ``auxillary`` group that
        :meth:`orgui.app.orGUI.orGUI.rocking_integrate` copied there. Which
        counters count as a monitor is the same setting the stationary
        integration and the reconstruction use, so all three normalize
        identically.

        A counter that is simply not there is not an error -- a backend that
        does not declare ``exposure_time`` in ``auxillary_counters`` never
        stored one. The names of the factors that did apply are returned so
        they can be saved beside ``F2_hkl``; without them a rod cannot be put
        on a common scale after the fact.

        :param dict aux: Auxiliary counters, each of shape ``(n_pts,)``.
        :param int size: Number of frames of the rocking scan.
        :returns: ``(divisor, applied)`` with the divisor of shape
            ``(size,)``.
        :rtype: tuple
        """
        config_target = self.database.config_target
        correction_state = getattr(config_target, "ctr_correction_state", None)
        monitor_names = tuple(
            getattr(
                correction_state,
                "monitor_corrections",
                getattr(config_target, "reconstruction_monitor_corrections", ()),
            ) or ()
        )

        normalize_exposure = getattr(
            correction_state,
            "normalize_exposure",
            getattr(config_target, "reconstruction_normalize_exposure", True),
        )
        exposure = aux.get("exposure_time") if normalize_exposure else None
        if normalize_exposure and exposure is None:
            logger.warning(
                "The rocking scan stores no exposure_time counter, so the "
                "integrated intensities are not normalized to counting time. "
                "They are then only comparable to other scans of the same "
                "duration. The scan backend decides this by declaring "
                "'exposure_time' in auxillary_counters."
            )

        monitors = {}
        for name in monitor_names:
            if name in aux:
                monitors[name] = aux[name]
            else:
                logger.warning(
                    "Monitor counter %r is configured but was not stored with "
                    "this rocking scan; skipping it.",
                    name,
                )

        return normalization_corrections.normalization_divisor(
            size, exposure_time=exposure, monitors=monitors
        )

    def _stored_detector(self, scangroup):
        """The detector geometry the rocking curves were measured with.

        Read from the configuration stored beside the scan, never from the
        current application state. The acceptance and the solid-angle factor
        are properties of the geometry the data was *taken* with, and a
        reduction run later -- from a batch script, or after another
        calibration has been loaded -- would otherwise silently use whatever
        detector happens to be loaded. Measured on a LaNiO3 rocking scan, that
        mistake scaled every ``Delta_gamma`` by 2.3 and left no trace in the
        output.

        :param scangroup: The scan group holding ``configuration``.
        :returns: A
            :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`,
            or ``None`` when the scan stores no detector configuration.
        :rtype: object or None
        """
        path = scangroup.name + "/configuration/instrument/detector_SXRD"
        if path not in self.database.nxfile:
            logger.warning(
                "This scan stores no detector configuration, so the "
                "out-of-plane acceptance cannot be calculated from the "
                "geometry the data was measured with. F2_hkl is left on the "
                "acceptance-blind scale."
            )
            return None
        try:
            return detector_from_nxdict(h5todict(self.database.nxfile, path))
        except Exception:
            logger.exception(
                "Cannot rebuild the detector geometry stored with this scan, "
                "so the out-of-plane acceptance cannot be calculated. F2_hkl "
                "is left on the acceptance-blind scale.",
                extra={"title": "Cannot read the stored detector geometry"},
            )
            return None

    def _rocking_solid_angle_mean(self, detector, scangroup, cnters, x, y):
        """Region mean of the solid-angle correction that was applied, or None.

        The solid-angle correction is applied to the *intensity* when the
        rocking curves are extracted, which is useful there: for a broad,
        non-rod feature a differential cross-section is what is wanted. It
        must not reach a structure factor, though, because a region sum is
        already the complete angular integral with every pixel weighted by the
        solid angle it subtends. So it is measured over the same regions here
        and divided back out of ``F2_hkl``. See
        ``doc/design/ctr_structure_factor_scale.md`` finding F6.

        Whether it was applied is a property of the *extraction*, not of the
        switches in this dialog, so it is read from the configuration snapshot
        stored with the scan rather than from the current GUI state -- as is
        the detector geometry it is measured over.

        :param detector: The geometry the scan was measured with, from
            :meth:`_stored_detector`.
        :param scangroup: The scan group holding the ``configuration`` written
            when the rocking curves were extracted.
        :param cnters: The ``rois`` group of the rocking scan.
        :param x: Region centre column per ``s`` point, in pixels.
        :param y: Region centre row per ``s`` point, in pixels.
        :returns: ``(mean, applied)`` -- the per-``s`` mean of
            :math:`1/\\widetilde{\\Omega}` and whether it will be divided out,
            or ``(None, False)`` when the correction was not applied or cannot
            be established.
        :rtype: tuple
        """
        path = scangroup.name + "/configuration/orgui/integration_corrections"
        try:
            group = h5todict(self.database.nxfile, path)
            if "json" in group:
                # Configurations written before the typed layout.
                raw = group["json"]
                if isinstance(raw, bytes):
                    raw = raw.decode()
                state = CorrectionState.from_dict(json.loads(str(raw)))
            else:
                state = corrections_from_nxdict(group)
            was_applied = bool(state.use_solid_angle)
        except Exception:
            logger.warning(
                "Cannot tell from this scan's stored configuration whether the "
                "solid angle correction was applied when the rocking curves "
                "were extracted, so it is not divided out of F2_hkl. If it was "
                "applied, F2_hkl carries the detector obliquity and will not "
                "agree with a stationary integration of the same rod."
            )
            return None, False

        if not was_applied:
            return None, False

        if detector is None:
            logger.warning(
                "The solid angle correction was applied to these rocking "
                "curves, but the detector geometry stored with the scan is "
                "not available to measure it over the regions of interest, so "
                "it is not divided out of F2_hkl."
            )
            return None, False

        hsize = np.asarray(cnters["hsize"][()], dtype=float)
        vsize = np.asarray(cnters["vsize"][()], dtype=float)
        if hsize.ndim > 1:
            hsize = hsize[:, 0]
        if vsize.ndim > 1:
            vsize = vsize[:, 0]

        mean = detector_corrections.roi_mean_inverse_solid_angle(
            detector,
            np.asarray(y, dtype=float),
            np.asarray(x, dtype=float),
            np.maximum(vsize, 1.0),
            np.maximum(hsize, 1.0),
        )
        return np.asarray(mean, dtype=float), True

    def _rocking_acceptance(self, detector, cnters, x, y):
        """Out-of-plane acceptance of every region of interest, in radian.

        A rocking scan intercepts a slice of rod proportional to
        :math:`\\Delta\\gamma` (Vlieg equation 20), so ``F2_hkl`` is only on
        the same scale as a stationary measurement once it is divided by it.

        The region centre and its vertical size are stored per ``s`` point.
        Note the coordinate order: ``surfaceAnglesPoint`` takes pyFAI
        dimension 1 first, which is the detector *row*, and orGUI's ``y`` is
        the row while ``x`` is the column -- every call site in the
        application passes them in that swapped order.

        New extractions store the true detector-arm scattering angles for
        every source frame. The frame nearest the calculated peak in
        ``(alpha, theta)`` supplies the arm position used here. Old databases
        have no such arrays; they retain the historical calibration-position
        result, with a warning that makes that fallback visible.

        :param detector: The geometry the scan was measured with, from
            :meth:`_stored_detector`; ``None`` leaves the acceptance out.
        :param cnters: The ``rois`` group of the rocking scan.
        :param x: Region centre column per ``s`` point, in pixels.
        :param y: Region centre row per ``s`` point, in pixels.
        :returns: ``(acceptance, applied)`` -- the acceptance in radian of
            shape ``(n_s,)``, or ``(None, False)`` when the stored geometry
            is unavailable.
        :rtype: tuple
        """
        if detector is None:
            return None, False

        vsize = cnters["vsize"][()]
        vsize = np.asarray(vsize, dtype=float)
        if vsize.ndim > 1:
            vsize = vsize[:, 0]
        alpha_pk = np.deg2rad(np.asarray(cnters["alpha_pk"][()], dtype=float))
        gamma_arm, delta_arm = RockingPeakIntegrator._rocking_peak_arm_angles(
            cnters
        )

        acceptance = acceptance_corrections.out_of_plane_acceptance(
            detector,
            np.asarray(y, dtype=float),
            np.asarray(x, dtype=float),
            vsize,
            alpha_pk,
            gamma_arm,
            delta_arm,
        )
        return np.asarray(acceptance, dtype=float), True

    @staticmethod
    def _rocking_peak_arm_angles(cnters):
        """Stored detector arm at each calculated rocking-curve peak.

        The saved arm arrays are in the primary-beam frame and carry an
        explicit ``rad`` group attribute. A rocking scan can run in either
        ``mu`` or ``th``; selecting the stored frame nearest the calculated
        peak in both ``(alpha, theta)`` coordinates handles either mode and a
        reversed scan without relying on a separate axis-name string.

        :param cnters: The stored rocking ``rois`` group.
        :returns: ``(gamma_arm, delta_arm)`` in radian, or ``(None, None)``
            for an old database that has no arm snapshot.
        :rtype: tuple
        """
        if "gamma_arm" not in cnters or "delta_arm" not in cnters:
            logger.warning(
                "This rocking extraction stores no detector-arm positions; "
                "evaluating its out-of-plane acceptance at the detector "
                "calibration position (legacy fallback)."
            )
            return None, None

        try:
            gamma_arm = np.asarray(cnters["gamma_arm"][()], dtype=float)
            delta_arm = np.asarray(cnters["delta_arm"][()], dtype=float)
            alpha_pk = np.asarray(cnters["alpha_pk"][()], dtype=float).reshape(-1)
            theta_pk = np.asarray(cnters["theta_pk"][()], dtype=float).reshape(-1)

            def at_peak(values):
                values = np.asarray(values, dtype=float)
                if values.ndim == 0:
                    return np.full(alpha_pk.shape, float(values))
                if values.shape == alpha_pk.shape:
                    return values

                alpha = np.asarray(cnters["alpha"][()], dtype=float)
                theta = np.asarray(cnters["theta"][()], dtype=float)
                expected = (alpha_pk.size, alpha.shape[-1])
                alpha = np.broadcast_to(alpha, expected)
                theta = np.broadcast_to(theta, expected)
                values = np.broadcast_to(values, expected)

                def angular_difference(actual, target):
                    return (actual - target[:, None] + 180.0) % 360.0 - 180.0

                distance = angular_difference(alpha, alpha_pk) ** 2
                distance += angular_difference(theta, theta_pk) ** 2
                if np.any(~np.isfinite(distance).any(axis=1)):
                    raise ValueError("no finite frame lies near a calculated peak")
                indices = np.nanargmin(distance, axis=1)
                return values[np.arange(alpha_pk.size), indices]

            gamma_arm = at_peak(gamma_arm)
            delta_arm = at_peak(delta_arm)
            attrs = getattr(cnters, "attrs", {})
            unit = attrs.get("detector_arm_unit", "rad")
            if isinstance(unit, bytes):
                unit = unit.decode()
            if unit == "deg":
                gamma_arm = np.deg2rad(gamma_arm)
                delta_arm = np.deg2rad(delta_arm)
            elif unit != "rad":
                raise ValueError(f"unsupported detector-arm unit {unit!r}")
            if not np.all(np.isfinite(gamma_arm)) or not np.all(
                np.isfinite(delta_arm)
            ):
                raise ValueError("detector-arm values are not finite")
            return gamma_arm, delta_arm
        except Exception:
            logger.warning(
                "Cannot resolve the detector-arm position stored with this "
                "rocking extraction; evaluating its out-of-plane acceptance "
                "at the detector calibration position.",
                exc_info=True,
            )
            return None, None

    def integrate(self):
        """Integrate rocking-scan ROIs.

        This shared path must stay safe in both GUI and CLI startup modes.
        Progress reporting is routed through :mod:`orgui.logger_utils` so CLI
        mode logs progress instead of opening modal dialogs.
        """
        if not self._currentRoInfo:
            raise ValueError("No rocking scan selected.")
        name = self._currentRoInfo["name"]
        h5_obj = self.database.nxfile[name]
        self._prepareFootprintAction(h5_obj)
        curves = self.get_all_ro_curves()
        correction_record = curves.get("correction_record")
        versioned_total_flux = correction_record is not None
        cnters = h5_obj["rois"]
        footprint_action = self.footprint_action
        apply_legacy_footprint = (
            not versioned_total_flux and footprint_action == FOOTPRINT_APPLY
        )

        s_array = cnters["s"][()]
        alpha = np.deg2rad(cnters["alpha"][()])
        # alpha_pk = cnters['alpha_pk'][()]
        delta = np.deg2rad(cnters["delta"][()])
        gamma = np.deg2rad(cnters["gamma"][()])

        axis = curves["axis"]
        scangroup = h5_obj.parent.parent
        aux = {}
        if "auxillary" in scangroup:
            for aux_n in scangroup["auxillary"]:
                if scangroup["auxillary"][aux_n].shape == axis.shape:
                    aux[aux_n] = scangroup["auxillary"][aux_n][()]

        if (
            self.lorentzButton.isChecked()
        ):  # force Lorentz positive, sign of integrated intensities is forced positive.
            # A mu scan rocks the incidence angle, which is how a
            # reflectivity curve is measured; a th scan rocks the sample.
            # The two take different Lorentz factors from the z-axis table,
            # and neither is the stationary-scan factor. Which factors each
            # mode applies is decided in one place, by mode_components.
            if curves["axisname"] == "mu":
                mode = measurement_corrections.REFLECTIVITY_ROCKING
            elif curves["axisname"] == "th":
                mode = measurement_corrections.ROCKING
            else:
                raise NotImplementedError()
            components = measurement_corrections.mode_components(
                mode, alpha=alpha, delta=delta, gamma=gamma
            )
            C_Lor = components["C_Lorentz"]
            C_rod = components["C_rod"]
        else:
            C_Lor = 1.0
            C_rod = 1.0

        if versioned_total_flux:
            # Q and H were resolved framewise at extraction and reconstructed
            # from their stored divisors by get_all_ro_curves. Applying the
            # live dialog here would scale the same photons a second time.
            C_flux_on_sample = 1.0
            C_illum_area = 1.0
        elif apply_legacy_footprint:
            L = self.integrationCorrection.sampleLength()  # sample size, m
            # Gaussian or measured beam profile, depending on the dialog.
            # C_flux_on_sample is saved as the diagnostic numerator of the
            # single applied active-area factor.
            profile = self.integrationCorrection.beamProfile()
            C_flux_on_sample, C_illum_area = profile.corrections(alpha, L)

        else:
            C_flux_on_sample = 1.0
            C_illum_area = 1.0

        if cnters["x"].ndim > 1:
            warnings.warn(
                "You are using an old data base orGUI v1.3.0-alpha"
                "X and Y pixel coordinates of rocking scans will be incorrect"
            )
            x = cnters["x"][:, 0][
                ()
            ]  # may provide fix of database here, if anyone asks
            y = cnters["y"][:, 0][
                ()
            ]  # may provide fix of database here, if anyone asks
        else:
            x = cnters["x"][()]
            y = cnters["y"][()]

        if versioned_total_flux:
            if (
                correction_record.normalization_status != "applied"
                or curves["illumination_action"]
                not in ("applied", "replaced", "applied_here")
                or curves["illumination_convention"] != "total_flux_H"
            ):
                raise ValueError(
                    "This versioned total-flux curve lacks an applied Q or H "
                    "divisor and cannot be reduced to a CTR structure factor."
                )
            C_norm = np.ones_like(curves["croibg"], dtype=np.float64)
            normalization_applied = tuple(
                correction_record.normalization_components
            )
        else:
            C_norm, normalization_applied = self._rocking_normalization(
                aux, axis.size
            )
            C_norm = np.broadcast_to(C_norm, np.shape(curves["croibg"])).copy()

        # Only F2_hkl is divided by these, so asking for them with the
        # Lorentz switch off would warn about a detector nothing needs.
        detector_acceptance, acceptance_applied = None, False
        solid_angle_mean, solid_angle_compensated = None, False
        if self.lorentzButton.isChecked():
            detector = self._stored_detector(scangroup)
            detector_acceptance, acceptance_applied = self._rocking_acceptance(
                detector, cnters, x, y
            )
            # New extractions carry their polarization-only CTR branch
            # directly. The separately estimated solid-angle compensation is
            # retained only for legacy database files that lack that branch.
            if "ctr_croibg" not in curves and not versioned_total_flux:
                solid_angle_mean, solid_angle_compensated = (
                    self._rocking_solid_angle_mean(
                        detector, scangroup, cnters, x, y
                    )
                )

        self.database.nxfile[self._currentRoInfo["name"] + "/integration/"]
        roi_info = h5todict(
            self.database.nxfile, self._currentRoInfo["name"] + "/integration/"
        )

        progress = logger_utils.create_progress_logger(
            self, s_array.size, "Integrating rocking scans"
        )

        result = _compute_rocking_integration(
            s_array,
            axis,
            curves["croibg"],
            curves["croibg_errors"],
            roi_info,
            aux,
            self.lorentzButton.isChecked(),
            apply_legacy_footprint,
            C_Lor=C_Lor,
            C_rod=C_rod,
            C_flux_on_sample=C_flux_on_sample,
            C_illum_area=C_illum_area,
            C_norm=C_norm,
            detector_acceptance=detector_acceptance,
            solid_angle_mean=solid_angle_mean,
            ctr_croibg_curves=curves.get("ctr_croibg"),
            ctr_croibg_errors_curves=curves.get("ctr_croibg_errors"),
            angle_unit="deg",
            progress_callback=progress.update,
            should_cancel=progress.wasCanceled,
        )
        progress.finish()

        int_data = result["int_data"]
        croi = result["croi"]
        croi_errors = result["croi_errors"]
        raw_croi = result["raw_croi"]
        raw_croi_errors = result["raw_croi_errors"]
        bgroi = result["bgroi"]
        bgroi_errors = result["bgroi_errors"]
        raw_bgroi = result["raw_bgroi"]
        raw_bgroi_errors = result["raw_bgroi_errors"]
        croibg = result["croibg"]
        croibg_errors = result["croibg_errors"]
        raw_croibg = result["raw_croibg"]
        raw_croibg_errors = result["raw_croibg_errors"]
        auxil = result["auxil"]
        if self.lorentzButton.isChecked():
            F2_hkl = result["F2_hkl"]
            F2_hkl_errors = result["F2_hkl_errors"]
            if versioned_total_flux and correction_record.scale_convention.startswith(
                "total_flux_calibrated"
            ):
                provenance = correction_record.profile_provenance
                try:
                    prefactor = measurement_corrections.total_flux_prefactor(
                        float(provenance["wavelength_angstrom"]),
                        float(provenance["unitcell_area_angstrom2"]),
                    )
                    efficiency = float(
                        provenance.get("detector_efficiency_assumed", 1.0)
                    )
                    transmission = float(
                        provenance.get("external_transmission_assumed", 1.0)
                    )
                    absolute_divisor = prefactor * efficiency * transmission
                    if not np.isfinite(absolute_divisor) or absolute_divisor <= 0:
                        raise ValueError("absolute divisor is not positive")
                except (KeyError, TypeError, ValueError) as error:
                    raise ValueError(
                        "The calibrated total-flux curve lacks valid wavelength, "
                        "surface unit-cell area, or detector-response provenance."
                    ) from error
                F2_hkl = F2_hkl / absolute_divisor
                F2_hkl_errors = F2_hkl_errors / absolute_divisor

        int_data["@NX_class"] = "NXdetector"

        if "H_1" in cnters:
            H_1 = cnters["H_1"][0][()]
            H_0 = cnters["H_0"][0][()]

            traj1 = {
                "@NX_class": "NXcollection",
                "@direction": " Integrated rocking scan along H_1*s + H_0 in reciprocal space",  # noqa: E501
                "H_1": H_1,
                "H_0": H_0,
                "s": s_array,
            }

            suffix = ""
            i = 0
            name1 = str(H_1) + "*s1+" + str(H_0)
        else:
            traj1 = {
                "@NX_class": "NXcollection",
                "@direction": " Integrated rocking intensity at fixed position",
                "s": s_array,
            }

            suffix = ""
            i = 0
            name1 = "static_rocking_scan"

        while (
            self._currentRoInfo["name"] + "/measurement/" + name1 + suffix
            in self.database.nxfile
        ):
            suffix = f"_{i}"
            i += 1
        availname1 = name1 + suffix

        datas1 = {
            "@NX_class": "NXdata",
            "sixc_angles": {
                "@NX_class": "NXpositioner",
                "alpha": cnters["alpha_pk"][()],
                "omega": cnters["omega_pk"][()],
                "theta": cnters["theta_pk"][()],
                "delta": cnters["delta_pk"][()],
                "gamma": cnters["gamma_pk"][()],
                "chi": cnters["chi_pk"][()],
                "phi": cnters["phi_pk"][()],
                "@unit": "deg",
            },
            "hkl": {
                "@NX_class": "NXcollection",
                "h": cnters["HKL_pk"][:, 0][()],
                "k": cnters["HKL_pk"][:, 1][()],
                "l": cnters["HKL_pk"][:, 2][()],
            },
            "counters": {
                "@NX_class": "NXdetector",
                "croibg": croibg,
                "croibg_errors": croibg_errors,
                "croi": croi,
                "croi_errors": croi_errors,
                "bgroi": bgroi,
                "bgroi_errors": bgroi_errors,
                "raw_croibg": raw_croibg,
                "raw_croibg_errors": raw_croibg_errors,
                "raw_croi": raw_croi,
                "raw_croi_errors": raw_croi_errors,
                "raw_bgroi": raw_bgroi,
                "raw_bgroi_errors": raw_bgroi_errors,
                "integrated": int_data,
            },
            "pixelcoord": {"@NX_class": "NXdetector", "x": x, "y": y},
            "auxillary": auxil,
            "trajectory": traj1,
            "@signal": "counters/croibg",
            "@axes": "trajectory/s",
            "@title": self._currentRoInfo["name"] + "_" + availname1,
            "@orgui_meta": "roi",
        }

        config_target = self.database.config_target
        config_snapshot = ConfigData.from_gui(config_target)
        measurement = {
            "@NX_class": "NXentry",
            "@default": availname1,
            availname1: datas1,
        }
        measurement[availname1]["configuration"] = config_snapshot.to_nxdict(
            role="integration", source="integration_save"
        )

        if self.lorentzButton.isChecked():
            measurement[availname1]["counters"]["F2_hkl"] = F2_hkl
            measurement[availname1]["counters"]["F2_hkl_errors"] = F2_hkl_errors
            measurement[availname1]["@signal"] = "counters/F2_hkl"
            # What the reduction actually divided out. A saved rod cannot be
            # put on a common scale with another scan after the fact without
            # this, and which factors were available depends on the scan.
            reduction = {
                "@NX_class": "NXcollection",
                "@mode": mode,
                "@angle_unit": "rad",
                "@normalization_applied": ",".join(normalization_applied) or "none",
                "@acceptance_applied": bool(acceptance_applied),
                "@solid_angle_compensated": bool(solid_angle_compensated),
                "@photon_curve_used": "ctr_croibg" in curves,
                "@active_area_applied": bool(
                    curves["illumination_action"]
                    in ("applied", "replaced", "applied_here")
                    if versioned_total_flux
                    else apply_legacy_footprint
                ),
            }
            if versioned_total_flux:
                reduction.update({
                    "@curve_algorithm": correction_record.algorithm,
                    "@scale_convention": correction_record.scale_convention,
                    "@illumination_action": curves["illumination_action"],
                    "@illumination_convention": (
                        curves["illumination_convention"] or "none"
                    ),
                    "@normalization_source": "stored_frame_divisor",
                    "@illumination_source": "stored_frame_divisor",
                })
            if detector_acceptance is not None:
                reduction["detector_acceptance"] = detector_acceptance
                reduction["@detector_acceptance_unit"] = "rad"
            if apply_legacy_footprint:
                reduction["sample_size"] = L
                reduction["@sample_size_unit"] = "m"
            measurement[availname1]["reduction"] = reduction

        self.database.add_nxdict(
            measurement,
            update_mode="modify",
            h5path=self._currentRoInfo["name"] + "/measurement",
        )
        self._refreshReductionCorrectionStatus()

    def onSliderValueChanged(self, ddict):
        try:
            self.plotRoCurve(ddict["idx"])
        except Exception:  # non essential: silent
            try:
                self._check_ro_present()
            except Exception:
                return

    def onRoiChanged(self, roi):
        try:
            self._check_ro_present()
        except Exception as e:
            logger.exception(
                "Cannot change roi limits",
                extra={
                    "title": "Cannot change roi limits",
                    "show_dialog": True,
                    "description": str(e),
                    "parent": self,
                    "dialog_level": logging.WARNING,
                },
            )
            return
        # roi_dict = self.get_roi1D_info(self._idx)
        roi_t = self.roiwidget.roiTable.roidict[roi]
        new_roi_dict = roi_t.toDict()
        if self._currentRoInfo:
            # from IPython import embed; embed()
            if self._currentRoInfo["name"] + "/integration/" in self.database.nxfile:
                h5grp = self.database.nxfile[
                    self._currentRoInfo["name"] + "/integration/"
                ]
                if roi in h5grp:
                    if (
                        h5grp[roi]["from"][self._idx] != new_roi_dict["from"]
                        or h5grp[roi]["to"][self._idx] != new_roi_dict["to"]
                    ):
                        h5grp[roi]["from"][self._idx] = new_roi_dict["from"]
                        h5grp[roi]["to"][self._idx] = new_roi_dict["to"]
                        h5grp[roi]["anchor"][self._idx] = True
                        self.roiwidget.roiTable.isModelSetting = True
                        roi_t.setAnchor(True)
                        self.roiwidget.roiTable._updateRoiInfo(roi_t.getID())
                        self.roiwidget.roiTable.isModelSetting = False
                        # print('set true: %s -> %s' % (h5grp[roi]['anchor'][self._idx], roi_t.isAnchor()))  # noqa: E501
                    else:
                        if (
                            h5grp[roi]["anchor"][self._idx] != roi_t.isAnchor()
                        ):  # toggled anchor
                            # print('toggled: %s -> %s' % (h5grp[roi]['anchor'][self._idx], roi_t.isAnchor()))  # noqa: E501
                            h5grp[roi]["anchor"][self._idx] = roi_t.isAnchor()
                    # else:
                    #    print('New Anchor set: True')
                    #    h5grp[roi]['anchor'][self._idx] = True
                    #    roi_t.setAnchor(True)
                    # roi_t.anchor_updating = False

        # with qt.QSignalBlocker(self.anchorROIButton):
        #    self.anchorROIButton.setChecked(True)

    @staticmethod
    def _legacy_rocking_curve_group(h5_obj):
        """Return legacy ``rois`` only when the sibling record permits it.

        Stage-2 records describe the unchanged legacy extraction and therefore
        explicitly allow this compatibility path. A future normalized record
        uses the same distinct sibling name but a different algorithm; this
        reducer must refuse to reinterpret that curve as unnormalized input.
        An incomplete record beside an existing legacy group remains usable as
        legacy data, with its provenance reported as unknown.
        """
        if CURVE_CORRECTIONS_GROUP in h5_obj:
            record = h5_obj[CURVE_CORRECTIONS_GROUP]
            version = int(record.attrs.get("orgui_schema_version", 0))
            contract = record.attrs.get("orgui_curve_contract", "")
            if isinstance(contract, bytes):
                contract = contract.decode()
            try:
                algorithm = record["identity"]["algorithm"][()]
                if isinstance(algorithm, bytes):
                    algorithm = algorithm.decode()
                else:
                    algorithm = str(algorithm)
            except (KeyError, TypeError, ValueError):
                algorithm = None
            if (
                version == CURVE_CORRECTIONS_SCHEMA_VERSION
                and contract == "frame_corrections"
                and algorithm is not None
                and not algorithm.startswith("legacy_rocking_roi_")
            ):
                raise ValueError(
                    "This rocking curve uses a versioned normalized contract "
                    "that the legacy reducer must not reinterpret. Use a "
                    "reducer that supports its correction record."
                )
            if algorithm is None:
                logger.warning(
                    "The rocking curve has an incomplete correction record; "
                    "using its preserved legacy ROI curve with unknown "
                    "provenance."
                )
        return h5_obj["rois"]

    def get_ro_curve(self, idx):
        """Return one rocking curve under its stored correction contract."""
        return self.get_all_ro_curves(idx)

    def get_all_ro_curves(self, idx=None):
        """Return rocking curves without silently changing their scale.

        New total-flux records are reconstructed from the immutable base
        curve and the exact stored Q/H divisors. Legacy records retain the
        historical ``rois`` path and are normalized later by
        :meth:`_rocking_normalization`.

        :param int or None idx: Read only this curve. It is shown in its
            stored footprint state; a pending footprint action is applied
            only by :meth:`integrate`.
        """
        name = self._currentRoInfo["name"]
        h5_obj = self.database.nxfile[name]
        if CURVE_CORRECTIONS_GROUP in h5_obj:
            record = RockingPeakIntegrator._storedCurveCorrectionRecord(h5_obj)
            if record is not None and record.algorithm.startswith(
                "framewise_ctr_total_flux_"
            ):
                if idx is None:
                    footprint = {
                        "footprint_action": getattr(
                            self, "footprint_action", FOOTPRINT_KEEP
                        ),
                        "replacement_illumination": getattr(
                            self, "replacement_illumination", None
                        ),
                        "replacement_convention": getattr(
                            self, "replacement_illumination_convention", None
                        ),
                        "allow_convention_change": bool(
                            getattr(
                                self, "allow_illumination_convention_change", False
                            )
                        ),
                    }
                else:
                    # (curve, frame) arrays: read one row; 1-D divisors are
                    # shared per frame and broadcast unchanged.
                    record = dataclasses.replace(record, **{
                        field: getattr(record, field)[idx]
                        for field in (
                            "base_croibg",
                            "base_croibg_variance",
                            "normalization_divisor",
                            "illumination_divisor",
                        )
                        if np.ndim(getattr(record, field)) == 2
                    })
                    footprint = {}
                curve, errors, action, convention = corrected_curve_from_record(
                    record, **footprint
                )
                return {
                    "axisname": self._currentRoInfo["axisname"],
                    "axis": self._currentRoInfo["axis"],
                    "croibg": curve,
                    "croibg_errors": errors,
                    "correction_record": record,
                    "illumination_action": action,
                    "illumination_convention": convention,
                }

        cnters = self._legacy_rocking_curve_group(h5_obj)
        rows = () if idx is None else idx
        curve = {
            "axisname": self._currentRoInfo["axisname"],
            "axis": self._currentRoInfo["axis"],
            "croibg": cnters["croibg"][rows],
            "croibg_errors": cnters["croibg_errors"][rows],
        }
        if "ctr_croibg" in cnters and "ctr_croibg_errors" in cnters:
            curve["ctr_croibg"] = cnters["ctr_croibg"][rows]
            curve["ctr_croibg_errors"] = cnters["ctr_croibg_errors"][rows]
        return curve

    # def onAnchorBtnToggled(self, state):
    #    if self._currentRoInfo:
    #        if self._currentRoInfo['name'] + '/integration/' in self.database.nxfile:
    #            h5grp = self.database.nxfile[self._currentRoInfo['name'] + '/integration/']  # noqa: E501
    #            for k in h5grp:
    #                if k.startswith('sig') or k.startswith('bg'):
    #                    h5grp[k]['anchor'][self._idx] = state
    #
    #    print(state)

    def onAddAllROI(self):
        try:
            self._check_ro_present()
        except Exception as e:  # non essential: silent
            logger.exception(
                "Cannot add ROIs",
                extra={
                    "title": "Cannot add ROIs",
                    "show_dialog": True,
                    "description": f"Cannot add ROIs:\n{e}",
                    "parent": self,
                    "dialog_level": logging.WARNING,
                },
            )
            return
        if self._currentRoInfo:
            try:
                self._addAllRoiDialog.set_axis(
                    self._currentRoInfo["axis"], self._currentRoInfo["axisname"]
                )
                if self._addAllRoiDialog.exec() == qt.QDialog.Accepted:
                    roidict = self._addAllRoiDialog.get_ranges()

                    if (
                        self._currentRoInfo["name"] + "/integration/"
                        in self.database.nxfile
                    ):
                        del self.database.nxfile[
                            self._currentRoInfo["name"] + "/integration/"
                        ]

                    ro_info = self.get_rocking_scan_info(self._currentRoInfo["name"])
                    sc_h5 = self.database.nxfile[ro_info["name"]]["rois"]
                    if ro_info["axisname"] == "mu":
                        pk_pos_exact = sc_h5["alpha_pk"][()]
                    elif ro_info["axisname"] == "th":
                        pk_pos_exact = sc_h5["theta_pk"][()]
                    elif ro_info["axisname"] == "chi":
                        pk_pos_exact = sc_h5["chi_pk"][()]
                    elif ro_info["axisname"] == "phi":
                        pk_pos_exact = sc_h5["phi_pk"][()]
                    else:
                        raise ValueError(
                            "Cannot estimate peak position: unknown scan axis {}".format(  # noqa: E501
                                ro_info["axisname"]
                            )  # noqa: E501
                        )

                    roi1Ddict = {
                        "@NX_class": "NXcollection",
                        "@info": "ROI information for the rocking scan integration",
                        "peakpos": pk_pos_exact,
                        "sig_1": {
                            "@NX_class": "NXcollection",
                            "from": pk_pos_exact + roidict["sig_1"]["from"],
                            "to": pk_pos_exact + roidict["sig_1"]["to"],
                            "anchor": np.zeros_like(pk_pos_exact, dtype=bool),
                        },
                        "bg_1": {
                            "@NX_class": "NXcollection",
                            "from": pk_pos_exact + roidict["bg_1"]["from"],
                            "to": pk_pos_exact + roidict["bg_1"]["to"],
                            "anchor": np.zeros_like(pk_pos_exact, dtype=bool),
                        },
                        "bg_2": {
                            "@NX_class": "NXcollection",
                            "from": pk_pos_exact + roidict["bg_2"]["from"],
                            "to": pk_pos_exact + roidict["bg_2"]["to"],
                            "anchor": np.zeros_like(pk_pos_exact, dtype=bool),
                        },
                    }

                    self.database.add_nxdict(
                        roi1Ddict,
                        update_mode="modify",
                        h5path=ro_info["name"] + "/integration",
                    )
                    self.plotRoCurve(self._idx)
            except Exception as e:
                logger.exception(
                    "Cannot add ROI",
                    extra={
                        "title": "Cannot add ROI",
                        "show_dialog": True,
                        "description": f"Cannot add ROI: {e}",
                        "detailed_text": traceback.format_exc(),
                        "parent": self,
                        "dialog_level": logging.WARNING,
                    },
                )

    def onAddROI(self):
        try:
            self._check_ro_present()
        except Exception as e:  # non essential: silent
            logger.exception(
                "Cannot add ROIs",
                extra={
                    "title": "Cannot add ROIs",
                    "show_dialog": True,
                    "description": f"Cannot add ROIs:\n{e}",
                    "parent": self,
                    "dialog_level": logging.WARNING,
                },
            )
            return
        if self._currentRoInfo:
            try:
                self._addRoiDialog.set_axis(
                    self._currentRoInfo["axis"], self._currentRoInfo["axisname"]
                )
                if self._addRoiDialog.exec() == qt.QDialog.Accepted:
                    roidict = self._addRoiDialog.get_roi()
                    ro_info = self.get_rocking_scan_info(self._currentRoInfo["name"])
                    sc_h5 = self.database.nxfile[ro_info["name"]]["rois"]
                    if ro_info["axisname"] == "mu":
                        pk_pos_exact = sc_h5["alpha_pk"][()]
                    elif ro_info["axisname"] == "th":
                        pk_pos_exact = sc_h5["theta_pk"][()]
                    elif ro_info["axisname"] == "chi":
                        pk_pos_exact = sc_h5["chi_pk"][()]
                    elif ro_info["axisname"] == "phi":
                        pk_pos_exact = sc_h5["phi_pk"][()]
                    else:
                        raise ValueError(
                            "Cannot estimate peak position: unknown scan axis {}".format(  # noqa: E501
                                ro_info["axisname"]
                            )  # noqa: E501
                        )

                    roi_dict = self.get_roi1D_info(self._idx)  # check existing ROIS

                    prefix = "sig" if roidict["signal"] else "bg"
                    no = 0
                    for k in roi_dict:
                        if k.startswith(prefix):
                            no = max(no, int(k.split("_")[1]))

                    no += 1

                    roi1Ddict = {
                        "@NX_class": "NXcollection",
                        "@info": "ROI information for the rocking scan integration",
                        "peakpos": pk_pos_exact,
                        f"{prefix}_{no}": {
                            "@NX_class": "NXcollection",
                            "from": pk_pos_exact + roidict["from"],
                            "to": pk_pos_exact + roidict["to"],
                            "anchor": np.zeros_like(pk_pos_exact, dtype=bool),
                        },
                    }
                    self.database.add_nxdict(
                        roi1Ddict,
                        update_mode="modify",
                        h5path=ro_info["name"] + "/integration",
                    )
                    self.plotRoCurve(self._idx)
            except Exception as e:
                logger.exception(
                    "Cannot add ROI",
                    extra={
                        "title": "Cannot add ROI",
                        "show_dialog": True,
                        "description": f"Cannot add ROI: {e}",
                        "detailed_text": traceback.format_exc(),
                        "parent": self,
                        "dialog_level": logging.WARNING,
                    },
                )

    def onDeleteROI(self):
        try:
            self._check_ro_present()
        except Exception as e:  # non essential: silent
            logger.warning(
                "Invalid rocking scan info",
                extra={
                    "title": "Cannot delete ROIs",
                    "show_dialog": True,
                    "parent": self,
                    "description": f"Cannot delete ROIs:\n{e}",
                },
            )
            return
        if self._currentRoInfo:
            roi_name = self.roiwidget.currentRoi.getName()
            if qt.QMessageBox.Yes == qt.QMessageBox.question(
                self, "Delete ROI", f"Are you sure you want to delete ROI {roi_name}?"
            ):
                if (
                    self._currentRoInfo["name"] + "/integration/" + roi_name
                    in self.database.nxfile
                ):
                    del self.database.nxfile[
                        self._currentRoInfo["name"] + "/integration/" + roi_name
                    ]
                    self.plotRoCurve(self._idx)
                else:
                    logger.error(
                        "Cannot delete ROI",
                        extra={
                            "title": "Cannot delete ROI",
                            "show_dialog": True,
                            "description": "Cannot delete ROI: no such ROI in database",
                            "parent": self,
                            "dialog_level": logging.WARNING,
                        },
                    )

    def onDeleteAllROI(self):
        try:
            self._check_ro_present()
        except Exception as e:  # non essential: silent
            logger.warning(
                "Invalid rocking scan info",
                extra={
                    "title": "Cannot delete ROIs",
                    "show_dialog": True,
                    "parent": self,
                    "description": f"Cannot delete ROIs:\n{e}",
                },
            )
            # qt.QMessageBox.warning(self,'Cannot delete ROIs', 'Cannot delete ROIs:\n%s' % e)  # noqa: E501
            return
        if self._currentRoInfo:
            if qt.QMessageBox.Yes == qt.QMessageBox.question(
                self, "Delete All ROIs", "Are you sure you want to delete all ROIs?"
            ):
                if (
                    self._currentRoInfo["name"] + "/integration/"
                    in self.database.nxfile
                ):
                    del self.database.nxfile[
                        self._currentRoInfo["name"] + "/integration/"
                    ]
                    self.plotRoCurve(self._idx)
                else:
                    logger.error(
                        "Cannot delete ROIs",
                        extra={
                            "title": "Cannot delete ROIs",
                            "show_dialog": True,
                            "description": "Cannot delete ROI: no ROIs in database",
                            "parent": self,
                            "dialog_level": logging.WARNING,
                        },
                    )

    def plotRoCurve(self, idx):
        ro_curve = self.get_ro_curve(idx)
        s = self._currentRoInfo["s"][idx]
        if self._currentRoInfo["type"] == "hklscan":
            hkl = self._currentRoInfo["H_1"] * s + self._currentRoInfo["H_0"]
            title = "Rocking scan at s = {}, HKL = [{:.2f} {:.2f} {:.2f}]".format(
                s, *hkl
            )  # noqa: E501
        else:
            hkl = self._currentRoInfo["HKL_pk"][idx]
            title = "Rocking scan at s = {}, HKL = [{:.2f} {:.2f} {:.2f}]".format(
                s, *hkl
            )  # noqa: E501
        self.plotROIselect.setGraphTitle(title)
        # print('before table clear')
        self.roiwidget.roiTable.clear()
        # print('after table clear')

        # print('before plot clear')
        self.plotROIselect.clear()
        # print('after plot clear')

        lbl = self.plotROIselect.addCurve(
            ro_curve["axis"],
            ro_curve["croibg"],
            legend=title,
            xlabel="{} / deg".format(self._currentRoInfo["axisname"]),
            ylabel="center ROI background subtracted / cnts",
            yerror=ro_curve["croibg_errors"],
            resetzoom=False,
        )

        self.plotROIselect.setActiveCurve(lbl)
        self._idx = idx

        with qt.QSignalBlocker(self.curveSlider):
            self.curveSlider.setIndex(idx)

        roi_dict = self.get_roi1D_info(idx)
        # with qt.QSignalBlocker(self.anchorROIButton):
        #    for k in roi_dict:
        #        dk = roi_dict[k]
        #        if dk['anchor']:
        #            self.anchorROIButton.setChecked(True)
        #            break
        #    else:
        #        self.anchorROIButton.setChecked(False)

        minfrom = np.inf
        maxto = -np.inf
        for key in roi_dict:
            if key.startswith("sig"):
                roi_d = roi_dict[key]
                roi = AnchorROI(
                    key,
                    fromdata=roi_d["from"],
                    todata=roi_d["to"],
                    type_=str(roi_d["type"]),
                    anchor=roi_d["anchor"],
                )
                minfrom = min(minfrom, roi_d["from"])
                maxto = max(maxto, roi_d["to"])
                self.roiwidget.roiTable.isModelSetting = True  # workaround for a signal loop, which resets anchor status to False  # noqa: E501
                # took me 2 days to concede, do not try to fix this bug. Workaround is enough  # noqa: E501
                self.roiwidget.roiTable.addRoi(roi)
                self.roiwidget.roiTable._updateRoiInfo(roi.getID())
                roi.sigChanged.connect(partial(self.onRoiChanged, key))
                self.roiwidget.roiTable.isModelSetting = False
            elif key.startswith("bg"):  # color?
                roi_d = roi_dict[key]
                roi = AnchorROI(
                    key,
                    fromdata=roi_d["from"],
                    todata=roi_d["to"],
                    type_=str(roi_d["type"]),
                    anchor=roi_d["anchor"],
                )
                minfrom = min(minfrom, roi_d["from"])
                maxto = max(maxto, roi_d["to"])
                self.roiwidget.roiTable.isModelSetting = True  # workaround for a signal loop, which resets anchor status to False  # noqa: E501
                # took me 2 days to concede, do not try to fix this bug. Workaround is enough  # noqa: E501
                self.roiwidget.roiTable.addRoi(roi)
                self.roiwidget.roiTable._updateRoiInfo(roi.getID())
                roi.sigChanged.connect(partial(self.onRoiChanged, key))
                self.roiwidget.roiTable.isModelSetting = False

        # toDo: set from GUI
        if self.autozoom_checkbox.isChecked():
            if np.all(np.isfinite([minfrom, maxto])):
                self.resetXZoomScaled(self.zoomslider.value())

    def resetXZoomScaled(self, value):
        try:
            roi_dict = self.get_roi1D_info(self._idx)
            minfrom = np.inf
            maxto = -np.inf
            for key in roi_dict:
                if key.startswith("sig"):
                    roi_d = roi_dict[key]
                    minfrom = min(minfrom, roi_d["from"])
                    maxto = max(maxto, roi_d["to"])
                elif key.startswith("bg"):  # color?
                    roi_d = roi_dict[key]
                    minfrom = min(minfrom, roi_d["from"])
                    maxto = max(maxto, roi_d["to"])
            if np.all(np.isfinite([minfrom, maxto])):
                add_rnge = (maxto - minfrom) * ((value**1.5) / 100)
                self.plotROIselect.setGraphXLimits(minfrom - add_rnge, maxto + add_rnge)
        except Exception:
            pass
            # print("Cannot zoom:" , e)
            # can fail, add better error handling later!

    def get_roi1D_info(self, idx):
        name = self._currentRoInfo["name"]
        if "integration" in self.database.nxfile[name]:
            h5_obj = self.database.nxfile[name]["integration"]

            roi1Ddict = {}
            for k in h5_obj:
                if k.startswith("sig") or k.startswith("bg"):
                    roi1Ddict[k] = {
                        "name": k,
                        "from": h5_obj[k]["from"][idx][()],
                        "to": h5_obj[k]["to"][idx][()],
                        "type": "signal",
                        "anchor": h5_obj[k]["anchor"][idx][()],
                    }
        else:
            roi1Ddict = {}
        return roi1Ddict

    # def estimate_roi1D_info(self, ro_info):
    #    # ToDo make settings available from GUI
    #
    #    # estimate peak pos
    #    axis_pk_pos = []
    #    axis_pk_pos_idx = []
    #
    #    sc_h5 = self.database.nxfile[ro_info['name']]['rois']
    #    if ro_info['axisname'] == 'mu':
    #        pk_pos_exact = sc_h5['alpha_pk'][()]
    #    elif ro_info['axisname'] == 'th':
    #        pk_pos_exact = sc_h5['theta_pk'][()]
    #    elif ro_info['axisname'] == 'chi':
    #        pk_pos_exact = sc_h5['chi_pk'][()]
    #    elif ro_info['axisname'] == 'phi':
    #        pk_pos_exact = sc_h5['phi_pk'][()]
    #    else:
    #        raise ValueError("Cannot estimate peak position: unknown scan axis %s" % ro_info['axisname'])  # noqa: E501
    #    #idx_pk = np.argmin(np.abs(sc_h5['axis'][()] - pk_pos_exact), axis=1)
    #
    #    # ToDo make settings available from GUI
    #    sig_width = 0.1 # deg
    #    bg_1 = -0.2 , -0.1
    #    bg_2 = 0.1 , 0.2
    #
    #    roi1Ddict = {
    #        "@NX_class": u"NXcollection",
    #        "@info" : u"ROI information for the rocking scan integration",
    #        "peakpos" : pk_pos_exact,
    #        'sig_1' : {
    #                "@NX_class": u"NXcollection",
    #                'from' :  pk_pos_exact - 0.5*sig_width,
    #                'to' : pk_pos_exact + 0.5*sig_width,
    #                'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
    #            },
    #        'bg_1' : {
    #                "@NX_class": u"NXcollection",
    #                'from' :  pk_pos_exact + bg_1[0],
    #                'to' : pk_pos_exact + bg_1[1],
    #                'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
    #            },
    #        'bg_2' : {
    #                "@NX_class": u"NXcollection",
    #                'from' :  pk_pos_exact + bg_2[0],
    #                'to' : pk_pos_exact + bg_2[1],
    #                'fixed' : np.zeros_like(pk_pos_exact, dtype=bool)
    #            }
    #    }
    #
    #    return roi1Ddict

    # def set_

    def onChangeRockingScan(self, name):
        try:
            self.set_roscan(name)
        except Exception as e:
            logger.exception(
                "Invalid rocking scan info",
                extra={
                    "title": "Invalid rocking scan info",
                    "show_dialog": True,
                    "parent": self,
                    "description": f"Invalid or missig rocking scan info: {name}.\n{e}",
                },
            )
            return

    def get_rocking_scan_info(self, name):
        if name not in self.database.nxfile:
            raise ValueError(f"scan {name} is not in the database")

        h5_obj = self.database.nxfile[name]
        if "orgui_meta" not in h5_obj.attrs or h5_obj.attrs["orgui_meta"] != "rocking":
            raise ValueError(f"scan {name} is not a valid rocking scan")

        scangroup = h5_obj.parent.parent
        positioners = scangroup["instrument/positioners"]
        if len(positioners) > 1:
            raise NotImplementedError(
                "Multiple positioner changes are not yet implemented"
            )

        axisname = list(positioners.keys())[0]
        axis = positioners[axisname][()]

        ddict = {"name": name, "axisname": axisname, "axis": axis}

        if "H_1" in h5_obj["rois"]:
            H_1 = h5_obj["rois"]["H_1"][()]
            H_0 = h5_obj["rois"]["H_0"][()]

            if not (
                np.allclose(H_1.T[0], H_1.T[0][0])
                and np.allclose(H_1.T[1], H_1.T[1][0])
                and np.allclose(H_1.T[2], H_1.T[2][0])
            ):
                raise ValueError("Rocking scan H_1 mismatch: Are these multiple scans?")

            if not (
                np.allclose(H_0.T[0], H_0.T[0][0])
                and np.allclose(H_0.T[1], H_0.T[1][0])
                and np.allclose(H_0.T[2], H_0.T[2][0])
            ):
                raise ValueError("Rocking scan H_0 mismatch: Are these multiple scans?")

            ddict["H_0"] = H_0[0]
            ddict["H_1"] = H_1[0]
            ddict["type"] = "hklscan"
            ddict["HKL_pk"] = h5_obj["rois"]["HKL_pk"][()]
        elif "s" in h5_obj["rois"]:
            ddict["type"] = "static"
            ddict["HKL_pk"] = h5_obj["rois"]["HKL_pk"][()]
        else:
            raise ValueError("Invalid scan: scan has no parameter s")

        s_array = h5_obj["rois"]["s"][()]
        ddict["s"] = s_array
        return ddict


class _ShapeParameter(NamedTuple):
    """One numeric control of an analytical beam shape.

    :param str label: Text shown next to the spin box.
    :param str suffix: Unit shown in the spin box, empty when dimensionless.
    :param float default: Value the control starts at.
    :param int decimals: Digits shown after the decimal point.
    :param float minimum: Smallest accepted value.
    :param float maximum: Largest accepted value.
    :param float scale: Factor converting the shown value to the SI value
        the beam-profile factories expect, ``1e-6`` for micrometers.
    """

    label: str
    suffix: str
    default: float
    decimals: int
    minimum: float
    maximum: float
    scale: float


class _BeamShape(NamedTuple):
    """An analytical beam shape offered by the corrections dialog."""

    name: str
    factory: object
    parameters: tuple


_MICRONS = (" microns", 4, 1e-6, 1e6, 1e-6)


def _width(label, default):
    """Build a width parameter shown in micrometers."""
    suffix, decimals, minimum, maximum, scale = _MICRONS
    return _ShapeParameter(label, suffix, default, decimals, minimum, maximum, scale)


#: Analytical beam shapes, in the order they appear in the dialog. The first
#: is the default and reproduces the Gaussian correction orGUI has always
#: applied; see :mod:`orgui.datautils.xrayutils.corrections.beamprofile`.
BEAM_SHAPES = (
    _BeamShape(
        "Gaussian",
        beamprofile.gaussian_profile,
        (_width("beam size (FWHM):", 20.0),),
    ),
    _BeamShape(
        "Top hat",
        beamprofile.top_hat_profile,
        (_width("beam width:", 20.0),),
    ),
    _BeamShape(
        "Trapezoid",
        beamprofile.trapezoid_profile,
        (_width("base width:", 30.0), _width("flat width:", 10.0)),
    ),
    _BeamShape(
        "Smoothed top hat",
        beamprofile.smoothed_top_hat_profile,
        (_width("beam width:", 20.0), _width("edge sigma:", 2.0)),
    ),
    _BeamShape(
        "Generalized normal",
        beamprofile.generalized_normal_profile,
        (
            _width("beam size (FWHM):", 20.0),
            _ShapeParameter("flatness:", "", 2.0, 3, 0.1, 100.0, 1.0),
        ),
    ),
    _BeamShape(
        "Skew normal",
        beamprofile.skew_normal_profile,
        (
            _width("beam size (FWHM):", 20.0),
            _ShapeParameter("skew:", "", 0.0, 3, -50.0, 50.0, 1.0),
        ),
    ),
)

#: Largest number of numeric controls any shape in :data:`BEAM_SHAPES` needs.
_MAX_SHAPE_PARAMETERS = max(len(shape.parameters) for shape in BEAM_SHAPES)


class IntegrationCorrectionsDialog(qt.QDialog):
    """Settings of the numerical active-area correction.

    The incident beam is described either by an analytical shape from
    :data:`BEAM_SHAPES` or by a beam profile measured at the beamline. Both
    are evaluated by :mod:`orgui.datautils.xrayutils.corrections.beamprofile`, which
    evaluates the illuminated surface integral over the projected sample
    footprint. The intercepted-flux fraction is available as a diagnostic
    numerator but is not applied as a second correction. Only a measured
    profile can represent a beam that is asymmetric or has several maxima.

    User-facing units are millimeters for the sample size and micrometers
    for beam widths and offsets; :meth:`beamProfile` converts to the meters
    the correction module works in.
    """

    #: File-column meaning of the loaded beam-profile file.
    CONTENT_PROFILE = "beam profile"
    CONTENT_HEIGHT_SCAN = "height scan (-dI/dz)"
    settingsChanged = qt.Signal()

    def __init__(self, parent=None):
        qt.QDialog.__init__(self, parent)
        verticalLayout = qt.QVBoxLayout(self)
        verticalLayout.setContentsMargins(0, 0, 0, 0)
        img = qutils.AspectRatioPixmapLabel(self)
        pixmp = qt.QPixmap(resources.getPath("incident_corrections.png"))
        img.setPixmap(pixmp)
        # A schematic, not the main content: capped so a wide-aspect image
        # cannot by itself push the dialog past a normal screen's height.
        img.setMaximumHeight(110)

        verticalLayout.addWidget(img)

        self._settings_save = None
        # Tabulated measured profile, in meters and arbitrary intensity units.
        self._profile_z = None
        self._profile_intensity = None
        # Per-shape memory, so switching shapes and back keeps the values.
        self._shape_values = {
            shape.name: [p.default for p in shape.parameters] for shape in BEAM_SHAPES
        }

        sizesLayout = qt.QGridLayout()
        sizesLayout.addWidget(qt.QLabel("Sample size L:"), 0, 0)
        self.L = qt.QDoubleSpinBox()
        self.L.setRange(0.00001, 1000000)
        self.L.setDecimals(4)
        self.L.setSuffix(" mm")
        self.L.setValue(5)
        self.L.setToolTip(
            "Sample size along the beam. With the beam profile, this sets the "
            "illuminated fraction of the projected footprint, which is the "
            "active-area correction for open post-sample slits."
        )
        sizesLayout.addWidget(self.L, 0, 1)

        sizesLayout.addWidget(qt.QLabel("Sample size W:"), 1, 0)
        self.W = qt.QDoubleSpinBox()
        self.W.setRange(0.00001, 1000000)
        self.W.setDecimals(4)
        self.W.setSuffix(" mm")
        self.W.setValue(5)
        self.W.setToolTip(
            "Sample size perpendicular to the beam, in the surface plane.\n"
            "The horizontal extent of the active area is taken to be this "
            "value, i.e. the beam is assumed at least as wide as the sample, "
            "so the sample bounds the illuminated width.\n"
            "Only the absolute active area in square meter uses it; a "
            "structure factor on a relative scale is unaffected."
        )
        self._legacyWidthToolTip = self.W.toolTip()
        sizesLayout.addWidget(self.W, 1, 1)

        self.legacyFluxLabel = qt.QLabel("Legacy flux density:")
        sizesLayout.addWidget(self.legacyFluxLabel, 2, 0)
        self.beamFlux = qt.QDoubleSpinBox()
        self.beamFlux.setRange(0.0, 1e30)
        self.beamFlux.setDecimals(3)
        self.beamFlux.setSuffix(" ph/(s·mm²)")
        self.beamFlux.setValue(0.0)
        self.beamFlux.setToolTip(
            "Legacy incident photon flux density at the sample position. "
            "This remains photons/(s mm²) and is never reinterpreted as "
            "total photons/s. New total-flux calibration is configured in "
            "the parent corrections dialog."
        )
        self._legacyFluxToolTip = self.beamFlux.toolTip()
        sizesLayout.addWidget(self.beamFlux, 2, 1)

        sizesLayout.addWidget(qt.QLabel("Horizontal interception:"), 3, 0)
        self.horizontalInterception = qt.QComboBox()
        self.horizontalInterception.addItem("Not specified", "")
        self.horizontalInterception.addItem("Full beam intercepted", "full")
        self.horizontalInterception.addItem(
            "Known intercepted fraction", "fraction"
        )
        self.horizontalInterception.setToolTip(
            "A vertical beam profile cannot establish how much of the beam "
            "is intercepted horizontally. New total-flux corrections require "
            "this choice to be explicit."
        )
        sizesLayout.addWidget(self.horizontalInterception, 3, 1)

        self.horizontalFractionLabel = qt.QLabel("Horizontal fraction:")
        sizesLayout.addWidget(self.horizontalFractionLabel, 4, 0)
        self.horizontalFraction = qt.QDoubleSpinBox()
        self.horizontalFraction.setRange(0.000001, 1.0)
        self.horizontalFraction.setDecimals(6)
        self.horizontalFraction.setSingleStep(0.01)
        self.horizontalFraction.setValue(1.0)
        self.horizontalFraction.setToolTip(
            "Known fraction of the full incident beam intercepted in the "
            "horizontal direction. Must be greater than zero and at most one."
        )
        sizesLayout.addWidget(self.horizontalFraction, 4, 1)

        modeLayout = qt.QHBoxLayout()
        self.analyticalButton = qt.QRadioButton("analytical beam shape")
        self.analyticalButton.setChecked(True)
        self.analyticalButton.setToolTip(
            "Describe the beam by an analytical distribution."
        )
        self.measuredButton = qt.QRadioButton("measured beam profile")
        self.measuredButton.setToolTip(
            "Use a beam profile measured at the beamline. Required for a "
            "beam that is asymmetric or has more than one maximum."
        )
        modeLayout.addWidget(self.analyticalButton)
        modeLayout.addWidget(self.measuredButton)

        # Left column: the beam and sample. Right column: where the sample
        # sits in the beam, and the resulting preview. Side by side rather
        # than one long stack, so the dialog fits a normal screen instead of
        # running off the bottom of it.
        leftColumn = qt.QVBoxLayout()
        leftColumn.addLayout(sizesLayout)
        leftColumn.addLayout(modeLayout)
        # Only the active beam model's settings take up space; the other is
        # hidden rather than merely disabled, which used to reserve room for
        # both at once.
        self.beamModelStack = qt.QStackedWidget()
        self.beamModelStack.addWidget(self._createShapeGroup())
        self.beamModelStack.addWidget(self._createProfileGroup())
        leftColumn.addWidget(self.beamModelStack)
        leftColumn.addStretch(1)

        rightColumn = qt.QVBoxLayout()
        rightColumn.addWidget(self._createCenteringGroup())
        rightColumn.addWidget(self._createPreviewGroup())

        columns = qt.QHBoxLayout()
        columns.addLayout(leftColumn, 1)
        columns.addLayout(rightColumn, 1)
        verticalLayout.addLayout(columns)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.onCancel)
        verticalLayout.addWidget(buttons)

        self.setLayout(verticalLayout)

        self.analyticalButton.toggled.connect(self._onModeChanged)
        self.horizontalInterception.currentIndexChanged.connect(
            self._onHorizontalInterceptionChanged
        )
        self.horizontalFraction.valueChanged.connect(self._settingsChanged)
        self.L.valueChanged.connect(self._settingsChanged)
        self.W.valueChanged.connect(self._settingsChanged)
        self.beamFlux.valueChanged.connect(self._settingsChanged)
        self.profileOffset.valueChanged.connect(self._settingsChanged)
        self.profileCenter.currentIndexChanged.connect(self._settingsChanged)
        self.analyticalButton.toggled.connect(self._settingsChanged)
        self._onShapeChanged()
        self._onModeChanged()
        self._onHorizontalInterceptionChanged()
        self._settings_save = self.settings()

    def _settingsChanged(self, *args):
        """Notify owners that the effective beam settings changed."""
        self.settingsChanged.emit()

    def _onHorizontalInterceptionChanged(self, *args):
        """Enable the fraction editor only for the explicit fraction mode."""
        fraction_mode = self.horizontalInterception.currentData() == "fraction"
        self.horizontalFractionLabel.setEnabled(fraction_mode)
        self.horizontalFraction.setEnabled(fraction_mode)
        self._updateLegacyControlState()
        self._updatePreview()
        self._settingsChanged()

    def horizontalInterceptionMode(self):
        """Return ``'full'``, ``'fraction'`` or ``None`` from the UI."""
        return self.horizontalInterception.currentData() or None

    def horizontalInterceptedFraction(self):
        """Return the stated horizontal fraction, or ``None`` if unresolved."""
        mode = self.horizontalInterceptionMode()
        if mode == "full":
            return 1.0
        if mode == "fraction":
            return self.horizontalFraction.value()
        return None

    def setHorizontalInterception(self, mode, fraction=None):
        """Restore the horizontal interception convention.

        :param str or None mode: ``'full'``, ``'fraction'`` or ``None``.
        :param float or None fraction: Fraction used by ``'fraction'`` mode.
        """
        index = self.horizontalInterception.findData(mode or "")
        with blockSignals([self.horizontalInterception, self.horizontalFraction]):
            self.horizontalInterception.setCurrentIndex(max(index, 0))
            if fraction is not None:
                self.horizontalFraction.setValue(float(fraction))
        self._onHorizontalInterceptionChanged()

    def setTotalFluxMode(self, enabled):
        """Show which legacy density controls are inactive in total-flux mode."""
        self._totalFluxMode = bool(enabled)
        self._updateLegacyControlState()

    def _updateLegacyControlState(self):
        """Apply total-flux/legacy enablement without changing stored values."""
        enabled = bool(getattr(self, "_totalFluxMode", False))
        self.beamFlux.setEnabled(not enabled)
        self.legacyFluxLabel.setEnabled(not enabled)
        full_horizontal = (
            enabled and self.horizontalInterceptionMode() == "full"
        )
        self.W.setEnabled(not full_horizontal)
        if full_horizontal:
            reason = (
                "Not a total-flux scale input. Horizontal interception is "
                "explicitly the full beam."
            )
            self.W.setToolTip(reason)
        elif enabled:
            self.W.setToolTip(
                "Retained for geometry and legacy density calculations; it "
                "does not set the total-flux scale."
            )
        else:
            self.W.setToolTip(self._legacyWidthToolTip)
        if enabled:
            self.beamFlux.setToolTip(
                "Legacy density value preserved for compatibility; the "
                "active total-flux convention uses photons/s in the parent "
                "corrections dialog."
            )
        else:
            self.beamFlux.setToolTip(self._legacyFluxToolTip)

    def _createShapeGroup(self):
        """Build the analytical-shape group box."""
        self.shapeGroup = qt.QGroupBox("Analytical shape parameters")
        grid = qt.QGridLayout()

        grid.addWidget(qt.QLabel("shape:"), 0, 0)
        self.shapeSelector = qt.QComboBox()
        self.shapeSelector.addItems([shape.name for shape in BEAM_SHAPES])
        grid.addWidget(self.shapeSelector, 0, 1)

        self.shapeParameterLabels = []
        self.shapeParameters = []
        for row in range(_MAX_SHAPE_PARAMETERS):
            label = qt.QLabel("")
            spin = qt.QDoubleSpinBox()
            grid.addWidget(label, row + 1, 0)
            grid.addWidget(spin, row + 1, 1)
            self.shapeParameterLabels.append(label)
            self.shapeParameters.append(spin)
            spin.valueChanged.connect(self._onShapeValueChanged)

        self.shapeGroup.setLayout(grid)
        self.shapeSelector.currentIndexChanged.connect(self._onShapeChanged)
        return self.shapeGroup

    def _createProfileGroup(self):
        """Build the measured-profile group box."""
        self.profileGroup = qt.QGroupBox("Measured profile file")
        grid = qt.QGridLayout()

        grid.addWidget(qt.QLabel("file:"), 0, 0)
        self.profileFileEdit = qt.QLineEdit()
        self.profileFileEdit.setToolTip(
            "Text file with a position column and an intensity column. "
            "Comment lines start with '#'."
        )
        grid.addWidget(self.profileFileEdit, 0, 1)
        self.profileBrowseButton = qt.QPushButton("browse")
        grid.addWidget(self.profileBrowseButton, 0, 2)

        grid.addWidget(qt.QLabel("file contains:"), 1, 0)
        self.profileContent = qt.QComboBox()
        self.profileContent.addItems([self.CONTENT_PROFILE, self.CONTENT_HEIGHT_SCAN])
        self.profileContent.setToolTip(
            "'height scan' differentiates the transmitted intensity of a "
            "sample height scan into the beam profile, and keeps only the "
            "range over which the sample cuts into the beam."
        )
        grid.addWidget(self.profileContent, 1, 1, 1, 2)

        grid.addWidget(qt.QLabel("position column unit:"), 2, 0)
        self.profileUnit = qt.QComboBox()
        self.profileUnit.addItems(["mm", "microns"])
        grid.addWidget(self.profileUnit, 2, 1, 1, 2)

        self.profileGroup.setLayout(grid)

        # GUI-only: user-triggered file dialog path.
        self.profileBrowseButton.clicked.connect(self._onBrowseProfile)
        self.profileFileEdit.editingFinished.connect(self._onProfileFileChosen)
        self.profileContent.currentIndexChanged.connect(self.loadProfile)
        self.profileUnit.currentIndexChanged.connect(self.loadProfile)
        return self.profileGroup

    def _createCenteringGroup(self):
        """Build the group box placing the sample within the beam."""
        group = qt.QGroupBox("Sample position in the beam")
        grid = qt.QGridLayout()

        grid.addWidget(qt.QLabel("sample centered on:"), 0, 0)
        self.profileCenter = qt.QComboBox()
        self.profileCenter.addItems(["centroid", "peak", "median"])
        self.profileCenter.setToolTip(
            "Point of the beam profile the center of the sample is aligned "
            "to. 'median' is the half-cut position an edge-scan alignment "
            "converges to. All three coincide for a symmetric beam."
        )
        grid.addWidget(self.profileCenter, 0, 1)

        grid.addWidget(qt.QLabel("sample offset:"), 1, 0)
        self.profileOffset = qt.QDoubleSpinBox()
        self.profileOffset.setRange(-1000000, 1000000)
        self.profileOffset.setDecimals(4)
        self.profileOffset.setSuffix(" microns")
        self.profileOffset.setValue(0)
        self.profileOffset.setToolTip(
            "Displacement of the sample center from the reference point "
            "above, positive toward larger positions."
        )
        grid.addWidget(self.profileOffset, 1, 1)

        group.setLayout(grid)
        self.profileCenter.currentIndexChanged.connect(self._updatePreview)
        self.profileOffset.valueChanged.connect(self._updatePreview)
        return group

    def _createPreviewGroup(self):
        """Build the preview of the beam profile currently described."""
        group = qt.QGroupBox("Beam profile")
        box = qt.QVBoxLayout()
        # The plot itself is built on first display only: this dialog is
        # constructed for every session, including headless CLI runs that
        # never show it.
        self.profilePlot = None
        self._previewLayout = box

        self.profileInfo = qt.QLabel("no beam profile")
        self.profileInfo.setWordWrap(True)
        box.addWidget(self.profileInfo)

        group.setLayout(box)
        return group

    def _onModeChanged(self):
        """Show only the settings of the selected beam model.

        The other group is hidden by switching the stacked page rather than
        merely disabled, so it stops reserving layout space it is not using.
        """
        analytical = self.analyticalButton.isChecked()
        self.beamModelStack.setCurrentWidget(
            self.shapeGroup if analytical else self.profileGroup
        )
        self._updatePreview()
        self._settingsChanged()

    def _onShapeChanged(self):
        """Relabel the numeric controls for the newly selected shape."""
        shape = self.currentShape()
        stored = self._shape_values[shape.name]
        widgets = self.shapeParameters + self.shapeParameterLabels
        with blockSignals(widgets):
            for index, spin in enumerate(self.shapeParameters):
                label = self.shapeParameterLabels[index]
                if index < len(shape.parameters):
                    parameter = shape.parameters[index]
                    label.setText(parameter.label)
                    spin.setDecimals(parameter.decimals)
                    spin.setRange(parameter.minimum, parameter.maximum)
                    spin.setSuffix(parameter.suffix)
                    spin.setValue(stored[index])
                    label.setVisible(True)
                    spin.setVisible(True)
                else:
                    label.setVisible(False)
                    spin.setVisible(False)
        self._updatePreview()
        self._settingsChanged()

    def _onShapeValueChanged(self):
        """Remember the edited values for the current shape and redraw."""
        shape = self.currentShape()
        self._shape_values[shape.name] = [
            self.shapeParameters[i].value() for i in range(len(shape.parameters))
        ]
        self._updatePreview()
        self._settingsChanged()

    def currentShape(self):
        """Return the selected analytical beam shape.

        :rtype: _BeamShape
        """
        return BEAM_SHAPES[max(self.shapeSelector.currentIndex(), 0)]

    # GUI-only: creates a plot widget, on the path that displays the dialog.
    def showEvent(self, event):
        """Create the profile preview plot the first time the dialog is shown.

        .. note::
           GUI-only. Everything the corrections need is available without
           the preview, so headless use never builds a plot widget.
        """
        if self.profilePlot is None:
            self.profilePlot = silx.gui.plot.Plot1D(self)
            self.profilePlot.setGraphXLabel("position rel. to sample center / microns")
            self.profilePlot.setGraphYLabel("normalized profile / mm$^{-1}$")
            # Small enough that the two-column dialog still fits a normal
            # screen; still tall enough to read the profile shape.
            self.profilePlot.setMinimumHeight(160)
            self._previewLayout.insertWidget(0, self.profilePlot)
            self._updatePreview()
        qt.QDialog.showEvent(self, event)

    # GUI-only: user-triggered dialog path.
    def _onBrowseProfile(self):
        """Pick a beam-profile file.

        .. note::
           GUI-only. Opens a blocking file dialog and must not be called
           from CLI, batch, or other shared non-interactive code.
        """
        filename, _ = qt.QFileDialog.getOpenFileName(
            self,
            "Open beam profile",
            os.path.dirname(self.profileFileEdit.text()),
            "Text files (*.dat *.txt *.csv);;All files (*)",
        )
        if filename:
            self.profileFileEdit.setText(filename)
            self._onProfileFileChosen()

    def _onProfileFileChosen(self):
        """Load the named profile file and make it the beam model.

        Choosing a profile is an explicit request to use it. Without this,
        the analytical shape stays selected and silently keeps correcting
        with its own beam, which is a correction that can be wrong by orders
        of magnitude without anything looking amiss.
        """
        if self.loadProfile():
            self.measuredButton.setChecked(True)

    def loadProfile(self):
        """Read the beam-profile file named in the dialog.

        Failures are reported as warnings and leave the dialog without a
        profile, so that an unreadable file cannot silently be integrated
        with stale data.

        :returns: ``True`` if a profile was loaded.
        :rtype: bool
        """
        path = self.profileFileEdit.text().strip()
        self._profile_z = None
        self._profile_intensity = None
        if not path:
            self._updatePreview()
            return False
        height_scan = self.profileContent.currentText() == self.CONTENT_HEIGHT_SCAN
        z_scale = 1e-3 if self.profileUnit.currentText() == "mm" else 1e-6
        try:
            z, intensity = beamprofile.read_profile_file(
                path, z_scale=z_scale, height_scan=height_scan
            )
            # Construct once here so a malformed profile is reported while
            # the dialog is open rather than in the middle of an integration.
            beamprofile.MeasuredBeamProfile(z, intensity)
        except Exception:
            logger.warning("Cannot read beam profile %s", path, exc_info=True)
            self._updatePreview()
            return False
        self._profile_z = z
        self._profile_intensity = intensity
        self._updatePreview()
        return True

    def _updatePreview(self):
        """Redraw the profile preview and its summary line.

        The preview plot only exists once the dialog has been shown; the
        summary line is kept up to date either way.
        """
        if self.profilePlot is not None:
            self.profilePlot.remove(kind="curve")
            self.profilePlot.remove(kind="marker")
        try:
            profile = self.beamProfile()
        except (ValueError, TypeError) as error:
            self.profileInfo.setText(str(error))
            return

        z, density = profile.profile_curve()
        centroid = profile.centroid_position
        if self.profilePlot is not None:
            self.profilePlot.addCurve(z * 1e6, density * 1e-3, legend="beam profile")
            # No text: the axis label already names the origin, and a
            # second label collides with the centroid when they are close.
            self.profilePlot.addXMarker(0.0, legend="sample center", color="black")
            if np.isfinite(centroid):
                self.profilePlot.addXMarker(
                    centroid * 1e6, legend="centroid", text="centroid", color="green"
                )

        summary = (
            f"FWHM {profile.fwhm * 1e6:.1f} microns, "
            f"extent {z[0] * 1e6:.1f} to {z[-1] * 1e6:.1f} microns "
            f"relative to the sample center"
        )
        rms = getattr(profile, "rms_width", None)
        if rms is not None and np.isfinite(rms):
            summary += f", rms width {rms * 1e6:.1f} microns"
        if np.isfinite(centroid):
            summary += f", centroid at {centroid * 1e6:+.1f} microns"
        else:
            summary += ", centroid undefined (the profile has no center of mass)"
        mode = self.horizontalInterceptionMode()
        if mode == "full":
            summary += "; horizontal fraction 1 (full beam intercepted)"
        elif mode == "fraction":
            summary += f"; horizontal fraction {self.horizontalFraction.value():.6g}"
        else:
            summary += "; horizontal interception not specified"
        summary += "; vertical fraction is evaluated per frame from L and incidence"
        self.profileInfo.setText(summary)

    def measuredProfile(self):
        """Return the loaded measured beam profile.

        :returns: The tabulated profile, referenced to the sample center
            chosen in the dialog.
        :rtype: orgui.datautils.xrayutils.corrections.beamprofile.MeasuredBeamProfile
        :raises ValueError: If no profile file has been loaded.
        """
        if self._profile_z is None:
            raise ValueError(
                "No beam profile loaded. Select a valid beam profile file in "
                "the footprint correction options, or switch back to an "
                "analytical beam shape."
            )
        return beamprofile.MeasuredBeamProfile(
            self._profile_z,
            self._profile_intensity,
            center=self.profileCenter.currentText(),
            offset=self.profileOffset.value() * 1e-6,  # microns -> m
        )

    def analyticalProfile(self):
        """Return the analytical beam profile described by the dialog.

        :rtype: orgui.datautils.xrayutils.corrections.beamprofile.BeamProfile
        :raises ValueError: If the shape rejects the entered parameters.
        """
        shape = self.currentShape()
        values = [
            self.shapeParameters[index].value() * parameter.scale
            for index, parameter in enumerate(shape.parameters)
        ]
        return shape.factory(
            *values,
            center=self.profileCenter.currentText(),
            offset=self.profileOffset.value() * 1e-6,  # microns -> m
        )

    def beamProfile(self):
        """Return the incident-beam profile selected in the dialog.

        :returns: An analytical or measured beam profile, in meters.
        :rtype: orgui.datautils.xrayutils.corrections.beamprofile.BeamProfile
        :raises ValueError: If the measured profile is selected but no
            profile file has been loaded, or if the analytical parameters
            do not describe a usable profile.
        """
        if self.analyticalButton.isChecked():
            return self.analyticalProfile()
        return self.measuredProfile()

    def sampleLength(self):
        """Sample size along the beam, converted from millimeter to **meter**.

        :rtype: float
        """
        return self.L.value() * 1e-3

    def setSampleLength(self, value):
        """Set the sample size along the beam, from **meter**.

        :param float value: Sample size in meter.
        """
        self.L.setValue(value * 1e3)

    def sampleWidth(self):
        """Sample size perpendicular to the beam, in **meter**.

        The dialog shows millimeter. This is the horizontal extent the active
        area is taken to have: the beam is assumed at least as wide as the
        sample, so the sample bounds the illuminated width rather than the
        beam. Only :meth:`activeArea` uses it.

        :rtype: float
        """
        return self.W.value() * 1e-3

    def setSampleWidth(self, value):
        """Set the sample size perpendicular to the beam, from **meter**.

        :param float value: Sample size in meter.
        """
        self.W.setValue(value * 1e3)

    def activeArea(self, alpha):
        """Illuminated sample area at incidence angle ``alpha``, in m^2.

        The dimensionless active-area correction applied to an integrated
        intensity is
        :meth:`~.beamprofile.BeamProfile.illuminated_area_fraction`; this is
        the same quantity carrying its area, which is what an *absolute*
        structure factor needs (issue #15). It assumes open post-sample
        slits, and the horizontal extent of :meth:`sampleWidth`.

        :param alpha: Incidence angle(s) in radian, any array shape.
        :returns: The active area in square meter, broadcast over ``alpha``.
        :rtype: numpy.ndarray
        """
        return activearea_corrections.beam_limited_area(
            alpha, self.sampleWidth(), self.sampleLength(), self.beamProfile()
        )

    def beamFluxDensity(self):
        """Incident flux density, converted to **photons/(s m^2)**.

        The dialog shows photons/(s mm^2). ``0`` means unmeasured; pass it as
        ``flux_density`` to :func:`~.measurement.scale_factor` only once it
        is nonzero, since ``0`` there would zero the absolute scale rather
        than fall back to the relative one.

        :rtype: float
        """
        return self.beamFlux.value() * 1e6

    def setBeamFluxDensity(self, value):
        """Set the incident flux density from **photons/(s m^2)**.

        :param float value: Flux density in photons/(s m^2).
        """
        self.beamFlux.setValue(value * 1e-6)

    def settings(self):
        """Return the dialog state as a plain dict.

        :rtype: dict
        """
        shape = self.currentShape()
        return {
            "L": self.L.value(),
            "W": self.W.value(),
            "beam_flux": self.beamFlux.value(),
            "horizontal_interception": self.horizontalInterceptionMode(),
            "horizontal_intercepted_fraction": (
                self.horizontalInterceptedFraction()
            ),
            "analytical": self.analyticalButton.isChecked(),
            "shape": shape.name,
            "shape_values": list(self._shape_values[shape.name]),
            "profile_file": self.profileFileEdit.text(),
            "profile_content": self.profileContent.currentText(),
            "profile_unit": self.profileUnit.currentText(),
            "profile_center": self.profileCenter.currentText(),
            "profile_offset": self.profileOffset.value(),
        }

    def setSettings(self, settings):
        """Restore the dialog state from :meth:`settings`.

        :param dict settings: State to apply. Missing keys are left alone.
        """
        widgets = [
            self.profileFileEdit,
            self.profileContent,
            self.profileUnit,
            self.profileCenter,
            self.profileOffset,
            self.horizontalInterception,
            self.horizontalFraction,
            self.analyticalButton,
            self.shapeSelector,
        ] + self.shapeParameters
        with blockSignals(widgets):
            if "L" in settings:
                self.L.setValue(settings["L"])
            if "W" in settings:
                self.W.setValue(settings["W"])
            if "beam_flux" in settings:
                self.beamFlux.setValue(settings["beam_flux"])
            if "horizontal_interception" in settings:
                mode = settings["horizontal_interception"]
                index = self.horizontalInterception.findData(mode or "")
                self.horizontalInterception.setCurrentIndex(max(index, 0))
            if (
                "horizontal_intercepted_fraction" in settings
                and settings["horizontal_intercepted_fraction"] is not None
            ):
                self.horizontalFraction.setValue(
                    settings["horizontal_intercepted_fraction"]
                )
            if "analytical" in settings:
                self.analyticalButton.setChecked(bool(settings["analytical"]))
                self.measuredButton.setChecked(not settings["analytical"])
            if "profile_file" in settings:
                self.profileFileEdit.setText(settings["profile_file"])
            for key, widget in (
                ("shape", self.shapeSelector),
                ("profile_content", self.profileContent),
                ("profile_unit", self.profileUnit),
                ("profile_center", self.profileCenter),
            ):
                if key in settings:
                    index = widget.findText(settings[key])
                    if index >= 0:
                        widget.setCurrentIndex(index)
            if "shape_values" in settings:
                self._shape_values[self.currentShape().name] = list(
                    settings["shape_values"]
                )
            if "profile_offset" in settings:
                self.profileOffset.setValue(settings["profile_offset"])
        self._onShapeChanged()
        self._onModeChanged()
        self._onHorizontalInterceptionChanged()
        self.loadProfile()

    def onOk(self):
        self._settings_save = self.settings()
        self.accept()

    def onCancel(self):
        if self._settings_save is not None:
            self.setSettings(self._settings_save)
        self.reject()


class IntegrationEstimator(qt.QDialog):
    def __init__(self, axis, axisname, parent=None):
        qt.QDialog.__init__(self, parent)

        layout = qt.QGridLayout()

        self._scan_label = qt.QLabel(
            f"{axisname}-scan: from {axis[0]} to {axis[-1]}\nroi ranges will be set relative to calculated peak position"  # noqa: E501
        )
        layout.addWidget(self._scan_label, 0, 0, 1, -1)
        layout.addWidget(qt.QLabel("center roi:"), 1, 0)
        layout.addWidget(qt.QLabel("from (rel):"), 1, 1)
        layout.addWidget(qt.QLabel("to (rel):"), 1, 3)

        layout.addWidget(qt.QLabel("bg roi 1:"), 2, 0)
        layout.addWidget(qt.QLabel("from (rel):"), 2, 1)
        layout.addWidget(qt.QLabel("to (rel):"), 2, 3)

        layout.addWidget(qt.QLabel("bg roi 2:"), 3, 0)
        layout.addWidget(qt.QLabel("from (rel):"), 3, 1)
        layout.addWidget(qt.QLabel("to (rel):"), 3, 3)

        self.croi_from = qt.QDoubleSpinBox()
        self.croi_from.setRange(-1000000, 1000000)
        self.croi_from.setDecimals(4)
        self.croi_from.setSuffix(" °")
        self.croi_from.setValue(-0.1)
        layout.addWidget(self.croi_from, 1, 2)

        self.croi_to = qt.QDoubleSpinBox()
        self.croi_to.setRange(-1000000, 1000000)
        self.croi_to.setDecimals(4)
        self.croi_to.setSuffix(" °")
        self.croi_to.setValue(0.1)
        layout.addWidget(self.croi_to, 1, 4)

        self.bgroi1_from = qt.QDoubleSpinBox()
        self.bgroi1_from.setRange(-1000000, 1000000)
        self.bgroi1_from.setDecimals(4)
        self.bgroi1_from.setSuffix(" °")
        self.bgroi1_from.setValue(-0.3)
        layout.addWidget(self.bgroi1_from, 2, 2)

        self.bgroi1_to = qt.QDoubleSpinBox()
        self.bgroi1_to.setRange(-1000000, 1000000)
        self.bgroi1_to.setDecimals(4)
        self.bgroi1_to.setSuffix(" °")
        self.bgroi1_to.setValue(-0.12)
        layout.addWidget(self.bgroi1_to, 2, 4)

        self.bgroi2_from = qt.QDoubleSpinBox()
        self.bgroi2_from.setRange(-1000000, 1000000)
        self.bgroi2_from.setDecimals(4)
        self.bgroi2_from.setSuffix(" °")
        self.bgroi2_from.setValue(0.12)
        layout.addWidget(self.bgroi2_from, 3, 2)

        self.bgroi2_to = qt.QDoubleSpinBox()
        self.bgroi2_to.setRange(-1000000, 1000000)
        self.bgroi2_to.setDecimals(4)
        self.bgroi2_to.setSuffix(" °")
        self.bgroi2_to.setValue(0.3)
        layout.addWidget(self.bgroi2_to, 3, 4)

        layout.addWidget(
            qt.QLabel("Attention: This will reset and override all existing ROIs!"),
            4,
            0,
            1,
            -1,
        )

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        layout.addWidget(buttons, 5, 0, 1, -1)

        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)

        self.setLayout(layout)

    def set_axis(self, axis, axisname):
        self._scan_label.setText(
            f"{axisname}-scan: from {axis[0]} to {axis[-1]}, roi ranges will be set relative to calculated peak position"  # noqa: E501
        )

    def _verify_ranges(self):
        if self.croi_from.value() > self.croi_to.value():
            raise ValueError(
                f"Invalid input of center roi: from {self.croi_from.value()} > to {self.croi_to.value()}"  # noqa: E501
            )
        if self.bgroi1_from.value() > self.bgroi1_to.value():
            raise ValueError(
                f"Invalid input of bg roi 1: from {self.bgroi1_from.value()} > to {self.bgroi1_to.value()}"  # noqa: E501
            )
        if self.bgroi2_from.value() > self.bgroi2_to.value():
            raise ValueError(
                f"Invalid input of bg roi 2: from {self.bgroi2_from.value()} > to {self.bgroi2_to.value()}"  # noqa: E501
            )
        return True

    def get_ranges(self):
        self._verify_ranges()
        ddict = {
            "sig_1": {
                "name": "sig_1",
                "from": self.croi_from.value(),
                "to": self.croi_to.value(),
            },
            "bg_1": {
                "name": "bg_1",
                "from": self.bgroi1_from.value(),
                "to": self.bgroi1_to.value(),
            },
            "bg_2": {
                "name": "bg_2",
                "from": self.bgroi2_from.value(),
                "to": self.bgroi2_to.value(),
            },
        }
        return ddict

    def onOk(self):
        try:
            self._verify_ranges()
        except Exception as e:
            logger.exception(
                "Invalid input",
                extra={
                    "title": "Invalid input",
                    "show_dialog": True,
                    "parent": self,
                    "description": str(e),
                },
            )
            return
        self.accept()


class ROICreatorDialog(qt.QDialog):
    def __init__(self, axis, axisname, parent=None):
        qt.QDialog.__init__(self, parent)

        layout = qt.QGridLayout()

        self._scan_label = qt.QLabel(
            f"{axisname}-scan: from {axis[0]} to {axis[-1]}, roi ranges will be set relative to calculated peak position"  # noqa: E501
        )
        layout.addWidget(self._scan_label, 0, 0, 1, -1)

        layout.addWidget(qt.QLabel("roi:"), 1, 0)
        layout.addWidget(qt.QLabel("from (rel):"), 1, 1)
        layout.addWidget(qt.QLabel("to (rel):"), 1, 3)

        self.roi_from = qt.QDoubleSpinBox()
        self.roi_from.setRange(-1000000, 1000000)
        self.roi_from.setDecimals(4)
        self.roi_from.setSuffix(" °")
        self.roi_from.setValue(-0.1)
        layout.addWidget(self.roi_from, 1, 2)

        self.roi_to = qt.QDoubleSpinBox()
        self.roi_to.setRange(-1000000, 1000000)
        self.roi_to.setDecimals(4)
        self.roi_to.setSuffix(" °")
        self.roi_to.setValue(0.1)
        layout.addWidget(self.roi_to, 1, 4)

        checkboxlayout = qt.QHBoxLayout()

        self.signalCheckbox = qt.QCheckBox("Signal")
        self.backgroundCheckbox = qt.QCheckBox("Background")

        self.btngroup = qt.QButtonGroup()
        self.btngroup.addButton(self.signalCheckbox)
        self.btngroup.addButton(self.backgroundCheckbox)
        self.btngroup.setExclusive(True)

        self.signalCheckbox.setChecked(True)

        checkboxlayout.addWidget(qt.QLabel("ROI type:"))
        checkboxlayout.addWidget(self.signalCheckbox)
        checkboxlayout.addWidget(self.backgroundCheckbox)

        layout.addLayout(checkboxlayout, 2, 0, 1, -1)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        layout.addWidget(buttons, 3, 0, 1, -1)

        buttons.button(qt.QDialogButtonBox.Ok).clicked.connect(self.onOk)
        buttons.button(qt.QDialogButtonBox.Cancel).clicked.connect(self.reject)

        self.setLayout(layout)

    def set_axis(self, axis, axisname):
        self._scan_label.setText(
            f"{axisname}-scan: from {axis[0]} to {axis[-1]}, roi ranges will be set relative to calculated peak position"  # noqa: E501
        )

    def _verify_ranges(self):
        if self.roi_from.value() > self.roi_to.value():
            raise ValueError(
                f"Invalid input of roi: from {self.roi_from.value()} > to {self.roi_to.value()}"  # noqa: E501
            )
        return True

    def get_roi(self):
        self._verify_ranges()
        signal = self.signalCheckbox.isChecked()

        ddict = {
            "from": self.roi_from.value(),
            "to": self.roi_to.value(),
            "signal": signal,
        }
        return ddict

    def onOk(self):
        try:
            self._verify_ranges()
        except Exception as e:
            logger.exception(
                "Invalid input",
                extra={
                    "title": "Invalid input",
                    "show_dialog": True,
                    "parent": self,
                    "description": str(e),
                },
            )
            return
        self.accept()


class AnchorROI(ROI):
    # sigAnchorChanged = qt.Signal()

    def __init__(self, name, fromdata=None, todata=None, type_=None, anchor=False):
        ROI.__init__(self, name, fromdata, todata, type_)
        self._anchor = bool(anchor)
        self.anchor_updating = False

    def isAnchor(self):
        return self._anchor

    def setAnchor(self, anchor):
        # print(self.getName(), self._anchor, anchor)
        if self._anchor != bool(anchor):
            self._anchor = bool(anchor)
            self.sigChanged.emit()

    def toDict(self):
        """

        :return: dict containing the roi parameters
        """
        ddict = super().toDict()
        ddict["anchor"] = self.isAnchor()
        return ddict

    @staticmethod
    def _fromDict(dic):
        assert "name" in dic
        roi = AnchorROI(name=dic["name"])
        roi._extraInfo = {}
        for key in dic:
            if key == "from":
                roi.setFrom(dic["from"])
            elif key == "to":
                roi.setTo(dic["to"])
            elif key == "type":
                roi.setType(dic["type"])
            elif key == "anchor":
                roi.setAnchor(dic["anchor"])
            else:
                roi._extraInfo[key] = dic[key]
        return roi


class SelectableROITable(ROITable):
    COLUMNS_INDEX = dict(
        [
            ("Anchor", 0),
            ("ID", 1),
            ("ROI", 2),
            ("Type", 3),
            ("From", 4),
            ("To", 5),
            ("Raw Counts", 6),
            ("Net Counts", 7),
            ("Raw Area", 8),
            ("Net Area", 9),
        ]
    )

    COLUMNS = list(COLUMNS_INDEX.keys())

    def __init__(self, parent=None, plot=None, rois=None):
        super().__init__(parent, plot, rois)
        self.isModelSetting = False

    def _updateRoiInfo(self, roiID):
        if self._userIsEditingRoi is True:
            return
        if roiID not in self._roiDict:
            return
        super()._updateRoiInfo(roiID)

        roi = self._roiDict[roiID]
        itemID = self._getItem(name="ID", roi=roi, row=None)

        self._getItem(name="Anchor", row=itemID.row(), roi=roi)
        # print('update roi info: ', roi.getName(), roi.isAnchor())
        self.setAnchorState(roiID, roi.isAnchor())

    def setAnchorState(self, roiID, anchor):
        if roiID not in self._roiDict:
            return
        roi = self._roiDict[roiID]
        itemID = self._getItem(name="ID", roi=roi, row=None)

        itemAnchor = self._getItem(name="Anchor", row=itemID.row(), roi=roi)
        # with qt.QSignalBlocker(self):
        if anchor:
            itemAnchor.setCheckState(qt.Qt.Checked)
        else:
            itemAnchor.setCheckState(qt.Qt.Unchecked)

    def addRoi(self, roi):
        """

        :param :class:`ROI` roi: roi to add to the table
        """
        assert isinstance(roi, AnchorROI)
        self._getItem(name="ID", row=None, roi=roi)
        self._roiDict[roi.getID()] = roi
        self._markersHandler.add(roi, _ColorRoiMarkerHandler(roi, self.plot))
        self._updateRoiInfo(roi.getID())
        callback = partial(WeakMethodProxy(self._updateRoiInfo), roi.getID())
        roi.sigChanged.connect(callback)
        # set it as the active one
        self.setActiveRoi(roi)

    def setRois(self, rois, order=None):
        """Set the ROIs by providing a dictionary of ROI information.

        The dictionary keys are the ROI names.
        Each value is a sub-dictionary of ROI info with the following fields:

        - ``"from"``: x coordinate of the left limit, as a float
        - ``"to"``: x coordinate of the right limit, as a float
        - ``"type"``: type of ROI, as a string (e.g "channels", "energy")


        :param roidict: Dictionary of ROIs
        :param str order: Field used for ordering the ROIs.
             One of "from", "to", "type".
             None (default) for no ordering, or same order as specified
             in parameter ``rois`` if provided as a dict.
        """
        assert order in [None, "from", "to", "type"]
        self.clear()

        # backward compatibility since 0.10.0
        if isinstance(rois, dict):
            for roiName, roi in rois.items():
                if isinstance(roi, AnchorROI):
                    _roi = roi
                else:
                    roi["name"] = roiName
                    _roi = AnchorROI._fromDict(roi)
                self.addRoi(_roi)
        else:
            for roi in rois:
                assert isinstance(roi, AnchorROI)
                self.addRoi(roi)
        self._updateMarkers()

    def load(self, filename):
        """
        Load ROI widget information from a file storing a dict of ROI.

        :param str filename: The file from which to load ROI
        """
        roisDict = dictdump.load(filename)
        rois = []

        # Remove rawcounts and netcounts from ROIs
        for roiDict in roisDict["ROI"]["roidict"].values():
            roiDict.pop("rawcounts", None)
            roiDict.pop("netcounts", None)
            rois.append(AnchorROI._fromDict(roiDict))

        self.setRois(rois)

    def _getItem(self, name, row, roi):
        if row:
            item = self.item(row, self.COLUMNS_INDEX[name])
        else:
            item = None
        if item:
            return item
        else:
            if name == "ID":
                assert roi
                if roi.getID() in self._roiToItems:
                    return self._roiToItems[roi.getID()]
                else:
                    # create a new row
                    row = self.rowCount()
                    self.setRowCount(self.rowCount() + 1)
                    item = qt.QTableWidgetItem(
                        str(roi.getID()), type=qt.QTableWidgetItem.Type
                    )
                    self._roiToItems[roi.getID()] = item
            elif name == "ROI":
                item = qt.QTableWidgetItem(
                    roi.getName() if roi else "", type=qt.QTableWidgetItem.Type
                )
                if roi.getName().upper() in ("ICR", "DEFAULT"):
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                else:
                    item.setFlags(
                        qt.Qt.ItemIsSelectable
                        | qt.Qt.ItemIsEnabled
                        | qt.Qt.ItemIsEditable
                    )
            elif name == "Type":
                item = qt.QTableWidgetItem(type=qt.QTableWidgetItem.Type)
                item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
            elif name in ("To", "From"):
                item = _FloatItem()
                if roi.getName().upper() in ("ICR", "DEFAULT"):
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                else:
                    item.setFlags(
                        qt.Qt.ItemIsSelectable
                        | qt.Qt.ItemIsEnabled
                        | qt.Qt.ItemIsEditable
                    )
            elif name in ("Raw Counts", "Net Counts", "Raw Area", "Net Area"):
                item = _FloatItem()
                item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
            elif name == "Anchor":
                item = qt.QTableWidgetItem(type=qt.QTableWidgetItem.Type)
                if roi.getName().upper() in ("ICR", "DEFAULT"):
                    item.setFlags(qt.Qt.ItemIsSelectable | qt.Qt.ItemIsEnabled)
                else:
                    item.setFlags(
                        qt.Qt.ItemIsSelectable
                        | qt.Qt.ItemIsEnabled
                        | qt.Qt.ItemIsUserCheckable
                    )
                    item.setIcon(resources.getQicon("anchor-ROI"))
                    item.setCheckState(qt.Qt.Unchecked)
            else:
                raise ValueError("item type not recognized")

            self.setItem(row, self.COLUMNS_INDEX[name], item)
            return item

    def _itemChanged(self, item):
        def getRoi():
            IDItem = self.item(item.row(), self.COLUMNS_INDEX["ID"])
            assert IDItem
            id = int(IDItem.text())
            assert id in self._roiDict
            roi = self._roiDict[id]
            return roi

        def signalChanged(roi):
            if self.activeRoi and roi.getID() == self.activeRoi.getID():
                self.activeROIChanged.emit()

        super()._itemChanged(item)

        self._userIsEditingRoi = True
        if not self.isModelSetting:
            if item.column() == self.COLUMNS_INDEX["Anchor"]:
                roi = getRoi()
                anchor = item.checkState() == qt.Qt.Checked
                if anchor != roi.isAnchor():
                    # print('anchor changed:', getRoi().getName() ,getRoi().getID())
                    roi.setAnchor(anchor)

        self._userIsEditingRoi = False


class CurvesROIWidget(qt.QWidget):
    """
    Widget displaying a table of ROI information.

    Implements also the following behavior:

    * if the roiTable has no ROI when showing create the default ICR one

    :param parent: See :class:`QWidget`
    :param str name: The title of this widget
    """

    sigROIWidgetSignal = qt.Signal(object)
    """Signal of ROIs modifications.

    Modification information if given as a dict with an 'event' key
    providing the type of events.

    Type of events:

    - AddROI, DelROI, LoadROI and ResetROI with keys: 'roilist', 'roidict'
    - selectionChanged with keys: 'row', 'col' 'roi', 'key', 'colheader',
      'rowheader'
    """

    sigROISignal = qt.Signal(object)

    def __init__(self, parent=None, name=None, plot=None):
        super().__init__(parent)
        if name is not None:
            self.setWindowTitle(name)
        self.__lastSigROISignal = None
        """Store the last value emitted for the sigRoiSignal. In the case the
        active curve change we need to add this extra step in order to make
        sure we won't send twice the sigROISignal.
        This come from the fact sigROISignal is connected to the
        activeROIChanged signal which is emitted when raw and net counts
        values are changing but are not embed in the sigROISignal.
        """
        assert plot is not None
        self._plotRef = weakref.ref(plot)
        self._showAllMarkers = False
        self.currentROI = None

        layout = qt.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.roiTable = SelectableROITable(self, plot=plot)
        rheight = self.roiTable.horizontalHeader().sizeHint().height()
        self.roiTable.setMinimumHeight(4 * rheight)
        layout.addWidget(self.roiTable)
        self._roiFileDir = qt.QDir.home().absolutePath()
        # self._showAllCheckBox.toggled.connect(self.roiTable.showAllMarkers)
        self.roiTable.showAllMarkers(True)

        self._isConnected = False  # True if connected to plot signals
        self._isInit = False

    def getROIListAndDict(self):
        return self.roiTable.getROIListAndDict()

    def getPlotWidget(self):
        """Returns the associated PlotWidget or None

        :rtype: Union[~silx.gui.plot.PlotWidget,None]
        """
        return None if self._plotRef is None else self._plotRef()

    @property
    def roiFileDir(self):
        """The directory from which to load/save ROI from/to files."""
        if not os.path.isdir(self._roiFileDir):
            self._roiFileDir = qt.QDir.home().absolutePath()
        return self._roiFileDir

    @roiFileDir.setter
    def roiFileDir(self, roiFileDir):
        self._roiFileDir = str(roiFileDir)

    def setRois(self, rois, order=None):
        return self.roiTable.setRois(rois, order)

    def getRois(self, order=None):
        return self.roiTable.getRois(order)

    def setMiddleROIMarkerFlag(self, flag=True):
        return self.roiTable.setMiddleROIMarkerFlag(flag)

    def _add(self):
        """Add button clicked handler"""

        def getNextRoiName():
            rois = self.roiTable.getRois(order=None)
            roisNames = []
            [roisNames.append(roiName) for roiName in rois]
            nrois = len(rois)
            if nrois == 0:
                return "ICR"
            else:
                i = 1
                newroi = f"newroi {i:d}"
                while newroi in roisNames:
                    i += 1
                    newroi = f"newroi {i:d}"
                return newroi

        roi = AnchorROI(name=getNextRoiName())

        if roi.getName() == "ICR":
            roi.setType("Default")
        else:
            roi.setType(self.getPlotWidget().getXAxis().getLabel())

        xmin, xmax = self.getPlotWidget().getXAxis().getLimits()
        fromdata = xmin + 0.25 * (xmax - xmin)
        todata = xmin + 0.75 * (xmax - xmin)
        if roi.isICR():
            fromdata, dummy0, todata, dummy1 = self._getAllLimits()
        roi.setFrom(fromdata)
        roi.setTo(todata)
        self.roiTable.addRoi(roi)

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "AddROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def _del(self):
        """Delete button clicked handler"""
        self.roiTable.deleteActiveRoi()

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "DelROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def _reset(self):
        """Reset button clicked handler"""
        self.roiTable.clear()
        old = self.blockSignals(True)  # avoid several sigROISignal emission
        self._add()
        self.blockSignals(old)

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "ResetROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def _load(self):
        """Load button clicked handler"""
        dialog = qt.QFileDialog(self)
        dialog.setNameFilters(["INI File  *.ini", "JSON File *.json", "All *.*"])
        dialog.setFileMode(qt.QFileDialog.ExistingFile)
        dialog.setDirectory(self.roiFileDir)
        if not dialog.exec():
            dialog.close()
            return

        # pyflakes bug http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=666494
        outputFile = dialog.selectedFiles()[0]
        dialog.close()

        self.roiFileDir = os.path.dirname(outputFile)
        self.roiTable.load(outputFile)

        # back compatibility pymca roi signals
        ddict = {}
        ddict["event"] = "LoadROI"
        ddict["roilist"] = self.roiTable.roidict.values()
        ddict["roidict"] = self.roiTable.roidict
        self.sigROIWidgetSignal.emit(ddict)
        # end back compatibility pymca roi signals

    def load(self, filename):
        """Load ROI widget information from a file storing a dict of ROI.

        :param str filename: The file from which to load ROI
        """
        self.roiTable.load(filename)

    def _save(self):
        # GUI-only: save button handler for the ROI widget.
        """Save button clicked handler"""
        dialog = qt.QFileDialog(self)
        dialog.setNameFilters(["INI File  *.ini", "JSON File *.json"])
        dialog.setFileMode(qt.QFileDialog.AnyFile)
        dialog.setAcceptMode(qt.QFileDialog.AcceptSave)
        dialog.setDirectory(self.roiFileDir)
        if not dialog.exec():
            dialog.close()
            return

        outputFile = dialog.selectedFiles()[0]
        extension = "." + dialog.selectedNameFilter().split(".")[-1]
        dialog.close()

        if not outputFile.endswith(extension):
            outputFile += extension

        if os.path.exists(outputFile):
            try:
                os.remove(outputFile)
            except OSError:
                msg = qt.QMessageBox(self)
                msg.setIcon(qt.QMessageBox.Critical)
                msg.setText(f"Input Output Error: {sys.exc_info()[1]}")
                msg.exec()
                return
        self.roiFileDir = os.path.dirname(outputFile)
        self.save(outputFile)

    def save(self, filename):
        """Save current ROIs of the widget as a dict of ROI to a file.

        :param str filename: The file to which to save the ROIs
        """
        self.roiTable.save(filename)

    def calculateRois(self, roiList=None, roiDict=None):
        """Compute ROI information"""
        return self.roiTable.calculateRois()

    def showAllMarkers(self, _show=True):
        self.roiTable.showAllMarkers(_show)

    def _getAllLimits(self):
        """Retrieve the limits based on the curves."""
        plot = self.getPlotWidget()
        curves = () if plot is None else plot.getAllCurves()
        if not curves:
            return 1.0, 1.0, 100.0, 100.0

        xmin, ymin = None, None
        xmax, ymax = None, None

        for curve in curves:
            x = curve.getXData(copy=False)
            y = curve.getYData(copy=False)
            if xmin is None:
                xmin = x.min()
            else:
                xmin = min(xmin, x.min())
            if xmax is None:
                xmax = x.max()
            else:
                xmax = max(xmax, x.max())
            if ymin is None:
                ymin = y.min()
            else:
                ymin = min(ymin, y.min())
            if ymax is None:
                ymax = y.max()
            else:
                ymax = max(ymax, y.max())

        return xmin, ymin, xmax, ymax

    def showEvent(self, event):
        self._visibilityChangedHandler(visible=True)
        qt.QWidget.showEvent(self, event)

    def hideEvent(self, event):
        self._visibilityChangedHandler(visible=False)
        qt.QWidget.hideEvent(self, event)

    def _visibilityChangedHandler(self, visible):
        """Handle widget's visibility updates.

        It is connected to plot signals only when visible.
        """
        if visible:
            # if no ROI existing yet, add the default one
            if self.roiTable.rowCount() == 0:
                old = self.blockSignals(True)  # avoid several sigROISignal emission
                self._add()
                self.blockSignals(old)
                self.calculateRois()

    def fillFromROIDict(self, *args, **kwargs):
        self.roiTable.fillFromROIDict(*args, **kwargs)

    def _emitCurrentROISignal(self):
        ddict = {}
        ddict["event"] = "currentROISignal"
        if self.roiTable.activeRoi is not None:
            ddict["ROI"] = self.roiTable.activeRoi.toDict()
            ddict["current"] = self.roiTable.activeRoi.getName()
        else:
            ddict["current"] = None

        if self.__lastSigROISignal != ddict:
            self.__lastSigROISignal = ddict
            self.sigROISignal.emit(ddict)

    @property
    def currentRoi(self):
        return self.roiTable.activeRoi


class _ColorRoiMarkerHandler(_RoiMarkerHandler):
    def __init__(self, roi, plot):
        super().__init__(roi, plot)
        if roi.getName().startswith("sig"):
            self._color = "red"
        elif roi.getName().startswith("bg"):
            self._color = "blue"
        else:
            self._color = "black"

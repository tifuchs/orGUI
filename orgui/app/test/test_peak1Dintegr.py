"""Regression tests for rocking-scan ROI integration error propagation.

These tests exercise :func:`orgui.app.peak1Dintegr._compute_rocking_integration`
directly, without a Qt application or an on-disk database: it is the pure-numpy
aggregation core factored out of
:meth:`orgui.app.peak1Dintegr.RockingPeakIntegrator.integrate` for exactly this
purpose. See `issue #34
<https://github.com/tifuchs/orGUI/issues/34>`_ and
``doc/design/rocking_integration_error_propagation.md`` for the full analysis
and fix plan these tests pin down.

The tests pin the corrections described in the design record and protect the
saved rocking intensities and uncertainties from regression.
"""

import os
import tempfile
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from silx.io.dictdump import dicttonx

from orgui.app.config_data import (
    CURVE_CORRECTIONS_GROUP,
    CurveCorrectionRecord,
    curve_correction_record_to_nxdict,
)
from orgui.app.peak1Dintegr import (
    RockingPeakIntegrator,
    _compute_rocking_integration,
    _trapz_impl,
)
from orgui.app.integration_corrections import (
    FOOTPRINT_APPLY,
    FOOTPRINT_KEEP,
    FOOTPRINT_REMOVE,
)


def _piecewise_curve(axis, regions, background=0.0):
    """Build a curve on ``axis`` that is constant on each ``(lo, hi, value)``
    closed interval in ``regions`` and ``background`` elsewhere.

    Piecewise-constant regions make the correct trapezoidal integral over any
    sub-interval fully contained in one region hand-computable exactly
    (``value * width``), independent of grid spacing.
    """
    curve = np.full_like(axis, background, dtype=float)
    for lo, hi, value in regions:
        curve[(axis >= lo) & (axis <= hi)] = value
    return curve


def _roi_info(peakpos, rois):
    """Build a ``roi_info`` mapping like the one read from the database.

    :param peakpos: scalar or per-``s`` array added to every ``from``/``to``
        offset, matching how ``peakpos``-relative ROIs are stored.
    :param dict rois: mapping of ROI name to ``(from_offset, to_offset)``.
    """
    peakpos = np.atleast_1d(np.asarray(peakpos, dtype=float))
    info = {}
    for name, (frm, to) in rois.items():
        info[name] = {"from": peakpos + frm, "to": peakpos + to}
    return info


def _reference_integration(
    s_array, axis, croibg_curves, croibg_errors_curves, roi_info
):
    """Independently-written, correct rocking-scan aggregation.

    This is the ground truth described in
    ``doc/design/rocking_integration_error_propagation.md``: closed-interval
    trapezoidal integration (endpoint included), locally-correct trapezoid
    error weights, background subtracted as a density (``sig_interval /
    bg_interval``, not squared into each ROI's own fractional width), and
    full sum-of-squares error accumulation. It intentionally does not support
    the Lorentz or footprint corrections (`aux` counters either) - the tests
    that need those compute their own expectation inline.
    """
    n_s = s_array.size
    per_roi = {}
    for roikey, roi in roi_info.items():
        cnts = np.zeros(n_s)
        cnts_err = np.zeros(n_s)
        interval = np.zeros(n_s)
        for i in range(n_s):
            idx_from = np.argmin(np.abs(axis - roi["from"][i]))
            idx_to = np.argmin(np.abs(axis - roi["to"][i]))
            if idx_from > idx_to:
                idx_from, idx_to = idx_to, idx_from
            sl = slice(idx_from, idx_to + 1)  # closed interval
            sub_axis = axis[sl]
            sub_cnts = croibg_curves[i][sl]
            sub_err = croibg_errors_curves[i][sl]
            interval[i] = sub_axis[-1] - sub_axis[0]
            cnts[i] = _trapz_impl(sub_cnts, sub_axis)

            dx = np.diff(sub_axis)
            w = np.empty_like(sub_axis)
            if w.size == 1:
                w[:] = 0.0
            else:
                w[0] = dx[0] / 2
                w[-1] = dx[-1] / 2
                w[1:-1] = (dx[:-1] + dx[1:]) / 2
            cnts_err[i] = np.sqrt(np.sum((sub_err * w) ** 2))
        per_roi[roikey] = {"cnts": cnts, "cnts_err": cnts_err, "interval": interval}

    sig_interval = sum(per_roi[k]["interval"] for k in per_roi if k.startswith("sig"))
    bg_interval = sum(per_roi[k]["interval"] for k in per_roi if k.startswith("bg"))

    croi = sum(per_roi[k]["cnts"] for k in per_roi if k.startswith("sig"))
    croi_err2 = sum(per_roi[k]["cnts_err"] ** 2 for k in per_roi if k.startswith("sig"))

    ratio = sig_interval / bg_interval
    bgroi = sum(per_roi[k]["cnts"] * ratio for k in per_roi if k.startswith("bg"))
    bgroi_err2 = sum(
        (per_roi[k]["cnts_err"] * ratio) ** 2 for k in per_roi if k.startswith("bg")
    )

    croi_errors = np.sqrt(croi_err2)
    bgroi_errors = np.sqrt(bgroi_err2)

    croibg = croi - bgroi
    croibg_errors = np.sqrt(croi_errors**2 + bgroi_errors**2)

    return {
        "croi": croi,
        "croi_errors": croi_errors,
        "bgroi": bgroi,
        "bgroi_errors": bgroi_errors,
        "croibg": croibg,
        "croibg_errors": croibg_errors,
    }


def test_background_aggregation_matches_reference_with_two_unequal_bg_rois():
    """doc/design/rocking_integration_error_propagation.md findings #1 and #5.

    With two background ROIs of unequal width, the correct background
    estimate is the density-weighted mean across both ROIs, rescaled to the
    signal window width. The code instead divides by ``bg_interval`` twice
    (once per-ROI, once for the total), so the subtracted background - and
    therefore ``croibg`` itself, not just its error - comes out wrong. The
    integration window must also include both ROI endpoints so it matches the
    reported interval (finding #5).
    """
    axis = np.linspace(-1.0, 1.0, 401)  # dx = 0.005; ROI edges below are exact points
    regions = [
        (-0.1, 0.1, 1000.0),  # sig_1
        (-0.4, -0.2, 100.0),  # bg_1, width 0.2
        (0.2, 0.5, 300.0),  # bg_2, width 0.3 (deliberately unequal to bg_1)
    ]
    curve = _piecewise_curve(axis, regions)
    errors = np.sqrt(np.abs(curve))

    s_array = np.array([0.0])
    croibg_curves = curve[None, :]
    croibg_errors_curves = errors[None, :]

    roi_info = _roi_info(
        0.0, {"sig_1": (-0.1, 0.1), "bg_1": (-0.4, -0.2), "bg_2": (0.2, 0.5)}
    )

    got = _compute_rocking_integration(
        s_array,
        axis,
        croibg_curves,
        croibg_errors_curves,
        roi_info,
        aux={},
        use_lorentz=False,
        use_footprint=False,
    )
    expected = _reference_integration(
        s_array, axis, croibg_curves, croibg_errors_curves, roi_info
    )

    np.testing.assert_allclose(got["bgroi"], expected["bgroi"])
    np.testing.assert_allclose(got["croibg"], expected["croibg"])


def test_background_error_scales_with_exposure_time():
    """doc/design/rocking_integration_error_propagation.md finding #2.

    ``bgroi_errors`` is a Poisson-derived quantity and must scale as
    ``sqrt(T)`` with an overall count/exposure-time factor ``T``, same as
    every other error in this pipeline. The code discards the accumulated
    background variance and substitutes ``sqrt(croi_errors)`` (the *signal*
    error, square-rooted an extra time), which does not scale that way.
    """
    axis = np.linspace(-1.0, 1.0, 401)
    roi_info = _roi_info(
        0.0, {"sig_1": (-0.1, 0.1), "bg_1": (-0.4, -0.2), "bg_2": (0.2, 0.4)}
    )
    s_array = np.array([0.0])

    def _curves(exposure_factor):
        regions = [
            (-0.1, 0.1, 1000.0 * exposure_factor),
            (-0.4, -0.2, 100.0 * exposure_factor),
            (0.2, 0.4, 300.0 * exposure_factor),
        ]
        curve = _piecewise_curve(axis, regions)
        errors = np.sqrt(np.abs(curve))  # Poisson: sigma ~ sqrt(rate * T)
        return curve[None, :], errors[None, :]

    c1, e1 = _curves(1.0)
    c100, e100 = _curves(100.0)

    r1 = _compute_rocking_integration(
        s_array, axis, c1, e1, roi_info, aux={}, use_lorentz=False, use_footprint=False
    )
    r100 = _compute_rocking_integration(
        s_array, axis, c100, e100, roi_info,
        aux={}, use_lorentz=False, use_footprint=False,
    )

    ratio = r100["bgroi_errors"] / r1["bgroi_errors"]
    np.testing.assert_allclose(ratio, 10.0, rtol=0.05)  # sqrt(100) == 10


def test_F2_hkl_errors_uses_corrected_croibg_errors_under_footprint_correction():
    """doc/design/rocking_integration_error_propagation.md finding #3.

    ``F2_hkl`` is built from the footprint-corrected ``croibg``; the paired
    error must be built from the correspondingly-corrected ``croibg_errors``,
    not the uncorrected ``raw_croibg_errors``. The two only agree when the
    footprint correction is identically 1 (i.e. the footprint option is
    off), which is presumably why this was never noticed.
    """
    axis = np.linspace(-1.0, 1.0, 401)
    regions = [(-0.1, 0.1, 1000.0), (-0.4, -0.2, 100.0), (0.2, 0.4, 300.0)]
    curve = _piecewise_curve(axis, regions)
    errors = np.sqrt(np.abs(curve))

    s_array = np.array([0.0])
    croibg_curves = curve[None, :]
    croibg_errors_curves = errors[None, :]

    roi_info = _roi_info(
        0.0, {"sig_1": (-0.1, 0.1), "bg_1": (-0.4, -0.2), "bg_2": (0.2, 0.4)}
    )

    # A footprint correction that is not identically 1 anywhere in the ROIs.
    C_flux_on_sample = np.full((1, axis.size), 0.5)
    C_illum_area = np.full((1, axis.size), 0.5)
    C_Lor = np.full((1, axis.size), 3.0)
    C_rod = np.full((1, axis.size), 1.0)

    result = _compute_rocking_integration(
        s_array,
        axis,
        croibg_curves,
        croibg_errors_curves,
        roi_info,
        aux={},
        use_lorentz=True,
        use_footprint=True,
        C_Lor=C_Lor,
        C_rod=C_rod,
        C_flux_on_sample=C_flux_on_sample,
        C_illum_area=C_illum_area,
    )

    # C_flux_on_sample is the numerator already contained in C_illum_area;
    # only the latter is an intensity divisor. With both arrays at 0.5 the
    # correction is therefore 2, not 4.
    np.testing.assert_allclose(result["croibg"], result["raw_croibg"] / 0.5)

    # F2_hkl = croibg / denom for some denom (C_Lorentz * C_rod_intersect);
    # recover it from the already-computed F2_hkl and croibg rather than
    # re-deriving C_Lorentz/C_rod_intersect's exact ROI-weighted average.
    denom = result["croibg"] / result["F2_hkl"]
    expected_F2_hkl_errors = result["croibg_errors"] / denom

    np.testing.assert_allclose(result["F2_hkl_errors"], expected_F2_hkl_errors)


def test_integrated_error_is_finite_when_raw_signal_integral_is_zero():
    """doc/design/rocking_integration_error_propagation.md finding #4.

    ``I_corr_error`` is computed as ``(I_raw_error / I_raw) * I_corr``, which
    is ``0/0 -> nan`` whenever an ROI's raw trapezoidal integral is exactly
    zero (e.g. an all-zero-count ROI, or one straddling a sign change). A
    propagated error must stay finite.
    """
    axis = np.linspace(-1.0, 1.0, 401)
    # sig_1's window (-0.1, 0.1) is left out of `regions`, so it is all-zero.
    regions = [(-0.4, -0.2, 100.0), (0.2, 0.4, 300.0)]
    curve = _piecewise_curve(axis, regions)
    errors = np.sqrt(np.abs(curve))

    s_array = np.array([0.0])
    roi_info = _roi_info(
        0.0, {"sig_1": (-0.1, 0.1), "bg_1": (-0.4, -0.2), "bg_2": (0.2, 0.4)}
    )

    result = _compute_rocking_integration(
        s_array,
        axis,
        curve[None, :],
        errors[None, :],
        roi_info,
        aux={},
        use_lorentz=False,
        use_footprint=False,
    )

    assert np.all(np.isfinite(result["croibg_errors"])), result["croibg_errors"]


def test_normalization_is_applied_inside_the_rocking_integral():
    """A varying counting time must be divided out per frame.

    Vlieg's rocking expression integrates ``N(omega)/(T M)`` over the rocking
    angle, so the divisor belongs under the integral sign. Dividing the
    finished integral by a mean counting time is only the same thing when the
    counting time is constant, and a scan whose exposure drifts is exactly the
    case the normalization exists for.

    The curve is flat at ``10`` and the exposure ramps as ``T = 2 + omega``,
    which makes both readings closed forms rather than a re-implementation of
    the code:

    * per frame, as required:
      ``Int 10/(2+omega) domega = 10 ln(2.4/1.6)`` over ``[-0.4, 0.4]``
    * divided afterwards by the mean exposure, which must not be what comes
      out: ``10 * 0.8 / 2 = 4`` exactly, since the mean of ``T`` over a
      symmetric interval is 2.

    They differ by 1.4 %, so the wrong one cannot pass.
    """
    axis = np.linspace(-1.0, 1.0, 401)
    curve = _piecewise_curve(axis, [(-0.5, 0.5, 10.0)])
    exposure = 2.0 + axis

    result = _compute_rocking_integration(
        np.array([0.0]),
        axis,
        curve[None, :],
        np.sqrt(curve)[None, :],
        _roi_info(0.0, {"sig_1": (-0.4, 0.4)}),
        {},
        False,
        False,
        C_norm=exposure[None, :],
        angle_unit="rad",
    )

    per_frame = 10.0 * np.log(2.4 / 1.6)
    np.testing.assert_allclose(result["croibg"], per_frame, rtol=1e-5)
    assert not np.isclose(per_frame, 4.0, rtol=1e-3)


def test_the_rocking_integral_is_converted_to_radian_for_f2():
    """``F2_hkl`` carries the radian integral; the stored intensity does not.

    The published expressions integrate the rocking angle in radian, but the
    motor axis is in degrees and ``croibg`` is kept in the unit it was
    measured in, so the ``180/pi`` shows up as the ratio between them.
    """
    axis = np.linspace(-1.0, 1.0, 401)
    curve = _piecewise_curve(axis, [(-0.5, 0.5, 10.0)])
    ones = np.ones((1, axis.size))

    result = _compute_rocking_integration(
        np.array([0.0]),
        axis,
        curve[None, :],
        np.sqrt(curve)[None, :],
        _roi_info(0.0, {"sig_1": (-0.4, 0.4)}),
        {},
        True,
        False,
        C_Lor=ones,
        C_rod=ones,
        angle_unit="deg",
    )

    np.testing.assert_allclose(
        result["F2_hkl"] / result["croibg"], np.deg2rad(1.0), rtol=1e-12
    )


def test_the_acceptance_divides_f2_and_only_f2():
    """``Delta_gamma`` scales the structure factor, not the intensity.

    The intensity column stays the measured integral; the acceptance is part
    of forming ``F2_hkl``, because it is the rod slice the region intercepted
    rather than anything about the counts.
    """
    axis = np.linspace(-1.0, 1.0, 401)
    curve = _piecewise_curve(axis, [(-0.5, 0.5, 10.0)])
    ones = np.ones((1, axis.size))
    acceptance = np.array([np.deg2rad(0.4)])

    common = dict(
        C_Lor=ones,
        C_rod=ones,
        angle_unit="rad",
    )
    args = (
        np.array([0.0]),
        axis,
        curve[None, :],
        np.sqrt(curve)[None, :],
        _roi_info(0.0, {"sig_1": (-0.4, 0.4)}),
        {},
        True,
        False,
    )
    without = _compute_rocking_integration(*args, **common)
    with_it = _compute_rocking_integration(
        *args, detector_acceptance=acceptance, **common
    )

    np.testing.assert_allclose(with_it["croibg"], without["croibg"], rtol=1e-12)
    np.testing.assert_allclose(
        with_it["F2_hkl"] * acceptance, without["F2_hkl"], rtol=1e-12
    )
    np.testing.assert_allclose(
        with_it["F2_hkl_errors"] * acceptance, without["F2_hkl_errors"], rtol=1e-12
    )


def test_joint_lorentz_rod_factor_uses_the_count_quadrature():
    """Correlated factors are averaged as a product on a nonuniform axis."""
    axis = np.array([-0.5, -0.2, 0.0, 0.1, 0.5])
    curve = np.full(axis.shape, 10.0)
    lorentz = np.array([[1.0, 1.4, 2.0, 2.5, 4.0]])
    rod = np.array([[4.0, 2.5, 2.0, 1.4, 1.0]])
    result = _compute_rocking_integration(
        np.array([0.0]),
        axis,
        curve[None, :],
        np.sqrt(curve)[None, :],
        _roi_info(0.0, {"sig_1": (-0.5, 0.5)}),
        {},
        True,
        False,
        C_Lor=lorentz,
        C_rod=rod,
        angle_unit="rad",
    )

    joint = _trapz_impl((lorentz * rod)[0], axis) / (axis[-1] - axis[0])
    separate = (
        _trapz_impl(lorentz[0], axis)
        * _trapz_impl(rod[0], axis)
        / (axis[-1] - axis[0]) ** 2
    )
    assert not np.isclose(joint, separate)
    np.testing.assert_allclose(result["F2_hkl"], result["croibg"] / joint)

    reversed_result = _compute_rocking_integration(
        np.array([0.0]),
        axis[::-1],
        curve[None, ::-1],
        np.sqrt(curve)[None, ::-1],
        _roi_info(0.0, {"sig_1": (-0.5, 0.5)}),
        {},
        True,
        False,
        C_Lor=lorentz[:, ::-1],
        C_rod=rod[:, ::-1],
        angle_unit="rad",
    )
    np.testing.assert_allclose(reversed_result["F2_hkl"], result["F2_hkl"])
    np.testing.assert_allclose(
        reversed_result["F2_hkl_errors"], result["F2_hkl_errors"]
    )


def test_polarization_only_rocking_branch_controls_f2():
    """Solid-angle changes remain in intensity but cannot enter new F2."""
    axis = np.linspace(-0.5, 0.5, 101)
    photon_curve = np.full((1, axis.size), 10.0)
    intensity_curve = photon_curve * 2.5
    errors = np.sqrt(photon_curve)
    common = dict(
        s_array=np.array([0.0]),
        axis=axis,
        croibg_errors_curves=errors,
        roi_info=_roi_info(0.0, {"sig_1": (-0.4, 0.4)}),
        aux={},
        use_lorentz=True,
        use_footprint=False,
        C_Lor=np.ones_like(photon_curve),
        C_rod=np.ones_like(photon_curve),
        ctr_croibg_curves=photon_curve,
        ctr_croibg_errors_curves=errors,
        angle_unit="rad",
    )

    with_solid_angle = _compute_rocking_integration(
        croibg_curves=intensity_curve, **common
    )
    without = _compute_rocking_integration(
        croibg_curves=photon_curve, **common
    )

    np.testing.assert_allclose(
        with_solid_angle["croibg"], 2.5 * without["croibg"]
    )
    np.testing.assert_allclose(with_solid_angle["F2_hkl"], without["F2_hkl"])
    np.testing.assert_allclose(
        with_solid_angle["F2_hkl_errors"], without["F2_hkl_errors"]
    )


class _RowOnlyDetector:
    """A detector whose exit angle depends only on pyFAI dimension 1.

    ``surfaceAnglesPoint`` takes the *row* first, so a gamma built from its
    first argument alone is a direct probe of whether the caller got the
    coordinate order right.
    """

    PER_ROW = 1e-4

    def surfaceAnglesPoint(self, x, y, alpha_i, gamma_arm=None, delta_arm=None):
        gamma = self.PER_ROW * np.asarray(x, dtype=float)
        return gamma, np.zeros_like(gamma)


def _stub(detector=None, monitors=()):
    """Minimal stand-in for the parts of the integrator the helpers touch."""
    ubcalc = None if detector is None else SimpleNamespace(detectorCal=detector)
    config_target = SimpleNamespace(
        ubcalc=ubcalc, reconstruction_monitor_corrections=tuple(monitors)
    )
    return SimpleNamespace(database=SimpleNamespace(config_target=config_target))


def test_the_acceptance_reads_the_row_from_y_and_the_column_from_x():
    """orGUI's ``y`` is the detector row, and pyFAI takes the row first.

    ``detvsize, dethsize = detector.shape`` at ``orGUI.py:1096`` fixes ``y``
    as the row and ``vsize`` as its extent, and every ``surfaceAnglesPoint``
    call in the application passes the two coordinates swapped. Getting this
    backwards yields a plausible but wrong acceptance, so it is pinned here:
    the result must follow ``vsize`` and ``y`` and ignore ``x``.
    """
    detector = _RowOnlyDetector()
    vsize = np.array([10.0, 40.0, 80.0])
    cnters = {
        "vsize": vsize,
        "alpha_pk": np.zeros(3),
    }

    acceptance, applied = RockingPeakIntegrator._rocking_acceptance(
        _stub(), detector, cnters,
        x=np.array([5.0, 300.0, 470.0]), y=np.full(3, 250.0)
    )

    assert applied is True
    np.testing.assert_allclose(acceptance, vsize * detector.PER_ROW, rtol=1e-12)


def test_a_missing_detector_leaves_the_acceptance_out_rather_than_failing():
    """A scan storing no detector geometry must not lose the integration."""
    cnters = {"vsize": np.array([40.0]), "alpha_pk": np.zeros(1)}

    acceptance, applied = RockingPeakIntegrator._rocking_acceptance(
        _stub(), None, cnters, x=np.array([100.0]), y=np.array([200.0])
    )

    assert acceptance is None
    assert applied is False


def test_the_rocking_normalization_uses_the_stored_counters():
    """Exposure and the configured monitors multiply into one divisor."""
    aux = {
        "exposure_time": np.array([2.0, 4.0]),
        "mondio": np.array([10.0, 5.0]),
        "unused": np.array([7.0, 7.0]),
    }

    divisor, applied = RockingPeakIntegrator._rocking_normalization(
        _stub(monitors=("mondio",)), aux, 2
    )

    np.testing.assert_allclose(divisor, [20.0, 20.0], rtol=1e-12)
    assert applied == ["exposure", "monitor:mondio"]


def test_a_missing_exposure_counter_is_skipped_and_recorded():
    """A backend that declares no exposure_time still integrates.

    ``P212_tools`` and the base ``Scan`` return an empty
    ``auxillary_counters``, so nothing was stored to normalize by. That is a
    scale the user has to know about, not a reason to fail the job, and
    ``applied`` is what records it.
    """
    divisor, applied = RockingPeakIntegrator._rocking_normalization(
        _stub(monitors=("mondio",)), {}, 3
    )

    np.testing.assert_allclose(divisor, np.ones(3), rtol=1e-12)
    assert applied == []


def test_the_detector_comes_from_the_scan_not_from_the_application():
    """The acceptance must use the geometry the data was measured with.

    ``Delta_gamma`` is a property of the detector the rocking curves were
    recorded on. Reading it from whatever calibration the application happens
    to hold makes a reduction run later -- from a batch script, or after a
    different calibration was loaded -- silently wrong: on a real LaNiO3 scan
    that scaled every acceptance by 2.3 with nothing in the output to show
    for it.
    """
    pyFAI = pytest.importorskip("pyFAI")
    from silx.io.dictdump import dicttonx

    from orgui.app.config_data import detector_to_nxdict
    from orgui.datautils.xrayutils import DetectorCalibration

    stored = DetectorCalibration.Detector2D_SXRD()
    stored.detector = pyFAI.detectors.Detector(
        pixel1=172e-6, pixel2=172e-6, max_shape=(619, 487)
    )
    stored.dist = 0.5
    stored.poni1, stored.poni2 = 0.05, 0.04
    stored.rot1 = stored.rot2 = stored.rot3 = 0.0
    stored.set_energy(15.0)
    stored.setAzimuthalReference(np.deg2rad(90.0))
    stored.setPolarization(0.0, 1.0)

    with tempfile.TemporaryDirectory() as folder:
        path = os.path.join(folder, "scan.h5")
        dicttonx(
            {"configuration": {"instrument": {
                "detector_SXRD": detector_to_nxdict(stored)}}},
            path,
            h5path="/61.1",
            update_mode="add",
        )
        with h5py.File(path, "r") as handle:
            integrator = SimpleNamespace(
                database=SimpleNamespace(
                    nxfile=handle,
                    # A *different* geometry in the application, which must
                    # not be the one that gets used.
                    config_target=SimpleNamespace(
                        ubcalc=SimpleNamespace(detectorCal="wrong detector")
                    ),
                )
            )
            got = RockingPeakIntegrator._stored_detector(
                integrator, handle["/61.1"]
            )

    assert got is not None
    assert got != "wrong detector"
    assert got.dist == pytest.approx(0.5)
    assert got.poni1 == pytest.approx(0.05)
    assert got.poni2 == pytest.approx(0.04)


def test_a_scan_without_a_stored_detector_gives_none():
    """An older database has no geometry to read; that is not a crash."""
    with tempfile.TemporaryDirectory() as folder:
        path = os.path.join(folder, "scan.h5")
        with h5py.File(path, "w") as handle:
            handle.create_group("/61.1")
        with h5py.File(path, "r") as handle:
            integrator = SimpleNamespace(
                database=SimpleNamespace(nxfile=handle, config_target=None)
            )
            assert RockingPeakIntegrator._stored_detector(
                integrator, handle["/61.1"]
            ) is None


def _scan_with_corrections(folder, corrections_group):
    """Write one scan whose stored configuration holds ``corrections_group``."""
    from silx.io.dictdump import dicttonx

    path = os.path.join(folder, "scan.h5")
    dicttonx(
        {"configuration": {"orgui": {
            "integration_corrections": corrections_group}}},
        path,
        h5path="/61.1",
        update_mode="add",
    )
    return path


def _solid_angle_detector():
    """A calibrated geometry the solid-angle correction can be measured on."""
    pyFAI = pytest.importorskip("pyFAI")
    from orgui.datautils.xrayutils import DetectorCalibration

    detector = DetectorCalibration.Detector2D_SXRD()
    detector.detector = pyFAI.detectors.Detector(
        pixel1=172e-6, pixel2=172e-6, max_shape=(619, 487)
    )
    detector.dist = 0.5
    detector.poni1, detector.poni2 = 0.05, 0.04
    detector.rot1 = detector.rot2 = detector.rot3 = 0.0
    detector.set_energy(15.0)
    return detector


def test_rocking_acceptance_uses_the_stored_arm_position():
    """Pin the rolled-detector counterexample from the physics review.

    Frame diffraction angles are stored in degrees, while the explicitly
    unit-tagged arm snapshot and the acceptance API use radians. At this
    deliberately oblique geometry, using the home arm understates the accepted
    gamma span by about 3.4 percent.
    """
    from orgui.datautils.xrayutils.corrections import acceptance

    detector = _solid_angle_detector()
    pixel = detector.detector.pixel1
    detector.dist = 0.15
    detector.poni1 = detector.detector.shape[0] * pixel / 2.0
    detector.poni2 = detector.detector.shape[1] * pixel / 2.0
    detector.rot3 = np.deg2rad(30.0)
    detector.reset()
    detector._cached_array = {}

    row = np.array([500.0])
    column = np.array([440.0])
    height = np.array([60.0])
    alpha_deg = np.array([10.0])
    gamma_arm = np.deg2rad([[0.0, 30.0, 60.0]])
    delta_arm = np.deg2rad([[0.0, 40.0, 50.0]])
    counters = {
        "vsize": height,
        "alpha_pk": alpha_deg,
        "theta_pk": np.array([5.1]),
        "alpha": np.array([[9.0, 10.0, 11.0]]),
        "theta": np.array([[4.0, 5.0, 6.0]]),
        "gamma_arm": gamma_arm,
        "delta_arm": delta_arm,
    }

    expected = acceptance.out_of_plane_acceptance(
        detector,
        row,
        column,
        height,
        np.deg2rad(alpha_deg),
        gamma_arm[:, 1],
        delta_arm[:, 1],
    )
    actual, applied = RockingPeakIntegrator._rocking_acceptance(
        _stub(), detector, counters, x=column, y=row
    )

    assert applied is True
    np.testing.assert_allclose(actual, expected, rtol=1e-10)


def test_legacy_rocking_acceptance_reports_calibration_fallback(caplog):
    """An old database remains usable but does not hide missing arm history."""
    detector = _RowOnlyDetector()
    counters = {"vsize": np.array([40.0]), "alpha_pk": np.zeros(1)}

    with caplog.at_level("WARNING"):
        actual, applied = RockingPeakIntegrator._rocking_acceptance(
            _stub(), detector, counters, x=np.array([240.0]), y=np.array([300.0])
        )

    assert applied is True
    np.testing.assert_allclose(actual, 40.0 * detector.PER_ROW, rtol=1e-12)
    assert "calibration position (legacy fallback)" in caplog.text


def test_stage2_record_explicitly_allows_the_legacy_curve(tmp_path):
    """The new provenance sibling does not alter current numerical defaults."""
    path = tmp_path / "legacy_curve_record.h5"
    with h5py.File(path, "w") as database:
        scan = database.create_group("scan")
        rois = scan.create_group("rois")
        record = scan.create_group("ctr_curve_v3")
        record.attrs["orgui_schema_version"] = 3
        record.attrs["orgui_curve_contract"] = "frame_corrections"
        identity = record.create_group("identity")
        identity.create_dataset("algorithm", data="legacy_rocking_roi_v1")

        selected = RockingPeakIntegrator._legacy_rocking_curve_group(scan)
        assert selected.name == rois.name


def test_versioned_normalized_curve_is_not_consumed_as_legacy(tmp_path):
    """Future Q/H-normalized curves require their explicit dispatch path."""
    path = tmp_path / "normalized_curve_record.h5"
    with h5py.File(path, "w") as database:
        scan = database.create_group("scan")
        scan.create_group("rois")
        record = scan.create_group("ctr_curve_v3")
        record.attrs["orgui_schema_version"] = 3
        record.attrs["orgui_curve_contract"] = "frame_corrections"
        identity = record.create_group("identity")
        identity.create_dataset("algorithm", data="total_flux_framewise_v1")

        with pytest.raises(ValueError, match="must not reinterpret"):
            RockingPeakIntegrator._legacy_rocking_curve_group(scan)


def test_versioned_total_flux_curve_uses_stored_divisors_once(tmp_path):
    """The supported reducer dispatch ignores the legacy sibling values."""
    base = np.array([[80.0, 96.0, 120.0], [40.0, 48.0, 60.0]])
    variance = np.array([[16.0, 36.0, 64.0], [4.0, 9.0, 16.0]])
    q = np.array([2.0, 4.0, 5.0])
    h = np.array([0.5, 0.25, 1.0])
    record = CurveCorrectionRecord(
        algorithm="framewise_ctr_total_flux_v1",
        output_quantity="rocking_ctr_photon_curve",
        scale_convention="total_flux_calibrated",
        normalization_status="applied",
        normalization_divisor=q,
        normalization_unit="photons",
        illumination_status="applied",
        illumination_divisor=h,
        illumination_convention="total_flux_H",
        base_croibg=base,
        base_croibg_variance=variance,
    )
    path = tmp_path / "total_flux_curve.h5"
    dicttonx(
        {
            "scan": {
                "@NX_class": "NXcollection",
                "rois": {
                    "@NX_class": "NXcollection",
                    "croibg": np.full_like(base, -999.0),
                    "croibg_errors": np.ones_like(base),
                },
                CURVE_CORRECTIONS_GROUP: curve_correction_record_to_nxdict(
                    record
                ),
            }
        },
        path,
    )

    with h5py.File(path, "r") as handle:
        integrator = SimpleNamespace(
            _currentRoInfo={
                "name": "/scan",
                "axisname": "mu",
                "axis": np.arange(3.0),
            },
            database=SimpleNamespace(nxfile=handle),
        )
        curves = RockingPeakIntegrator.get_all_ro_curves(integrator)

    np.testing.assert_allclose(curves["croibg"], base / q / h)
    np.testing.assert_allclose(curves["croibg_errors"], np.sqrt(variance) / q / h)
    assert curves["correction_record"].algorithm == record.algorithm
    assert curves["illumination_action"] == "applied"


def test_one_total_flux_curve_is_read_lazily_in_its_stored_state(tmp_path):
    """Selecting a curve reads its row only and ignores a pending action."""
    base = np.array([[80.0, 96.0, 120.0], [40.0, 48.0, 60.0]])
    variance = np.array([[16.0, 36.0, 64.0], [4.0, 9.0, 16.0]])
    q = np.array([2.0, 4.0, 5.0])
    h = np.array([[0.5, 0.25, 1.0], [0.8, 0.4, 0.2]])
    record = CurveCorrectionRecord(
        algorithm="framewise_ctr_total_flux_v1",
        output_quantity="rocking_ctr_photon_curve",
        scale_convention="total_flux_calibrated",
        normalization_status="applied",
        normalization_divisor=q,
        normalization_unit="photons",
        illumination_status="applied",
        illumination_divisor=h,
        illumination_convention="total_flux_H",
        base_croibg=base,
        base_croibg_variance=variance,
    )
    path = tmp_path / "lazy_curve.h5"
    dicttonx(
        {"scan": {CURVE_CORRECTIONS_GROUP: curve_correction_record_to_nxdict(record)}},
        path,
    )

    with h5py.File(path, "r") as handle:
        integrator = SimpleNamespace(
            _currentRoInfo={"name": "/scan", "axisname": "mu", "axis": np.arange(3.0)},
            database=SimpleNamespace(nxfile=handle),
            footprint_action=FOOTPRINT_APPLY,
            replacement_illumination=None,
        )
        stored = RockingPeakIntegrator._storedCurveCorrectionRecord(handle["scan"])
        assert isinstance(stored.base_croibg, h5py.Dataset)
        assert isinstance(stored.illumination_divisor, h5py.Dataset)
        curve = RockingPeakIntegrator.get_all_ro_curves(integrator, 1)

    np.testing.assert_allclose(curve["croibg"], base[1] / q / h[1])
    np.testing.assert_allclose(curve["croibg_errors"], np.sqrt(variance[1]) / q / h[1])
    assert curve["illumination_action"] == "applied"


def test_unknown_legacy_curve_only_allows_keep_and_requests_reextraction():
    """Unknown provenance cannot expose apply/remove as safe operations."""
    state = RockingPeakIntegrator._correctionUiState(
        None, record_present=False
    )

    assert state["actions"] == (FOOTPRINT_KEEP,)
    assert state["normalization"] == "Unknown (legacy data)"
    assert "Re-extract images" in state["details"]


def test_known_legacy_curve_allows_current_legacy_footprint_only():
    """A typed legacy record retains the old, explicitly labeled workflow."""
    record = CurveCorrectionRecord(
        algorithm="legacy_rocking_roi_v1",
        output_quantity="legacy_roi_curve",
        scale_convention="legacy_relative",
    )
    state = RockingPeakIntegrator._correctionUiState(record)

    assert state["actions"] == (FOOTPRINT_KEEP, FOOTPRINT_APPLY)
    assert state["normalization"] == "Applied later by the legacy reducer"
    assert FOOTPRINT_REMOVE not in state["actions"]


def test_reversible_total_flux_curve_exposes_keep_apply_and_remove():
    """All three actions are safe only with a known reversible base curve."""
    record = CurveCorrectionRecord(
        algorithm="framewise_ctr_total_flux_v1",
        output_quantity="rocking_ctr_photon_curve",
        scale_convention="total_flux_relative",
        normalization_status="applied",
        illumination_status="applied",
        illumination_divisor=np.ones(3),
        illumination_convention="total_flux_H",
        base_croibg=np.ones((2, 3)),
        base_croibg_variance=np.ones((2, 3)),
    )
    state = RockingPeakIntegrator._correctionUiState(record)

    assert state["actions"] == (
        FOOTPRINT_KEEP,
        FOOTPRINT_APPLY,
        FOOTPRINT_REMOVE,
    )
    assert state["normalization"] == "Applied during extraction"
    assert state["footprint"] == "Applied during extraction"
    assert "total_flux_H" in state["details"]


@pytest.mark.parametrize("applied_at_extraction", [True, False])
def test_the_typed_corrections_group_is_read_back(applied_at_extraction):
    """The F6 compensation must survive the layout it is stored in.

    Whether the solid-angle correction was applied to the intensity is read
    from the configuration written at extraction time. When that group
    changed from a single JSON string to typed datasets, this reader kept
    parsing the JSON and silently stopped compensating -- the unit tests
    passed because the layout and the reduction were only ever tested apart.
    So this writes the group with the *current* writer and reads it back
    through the reducer.
    """
    from orgui.app.config_data import CorrectionState, corrections_to_nxdict

    detector = _solid_angle_detector()
    state = CorrectionState(use_solid_angle=applied_at_extraction)
    cnters = {"hsize": np.full(3, 20.0), "vsize": np.full(3, 5.0)}

    with tempfile.TemporaryDirectory() as folder:
        path = _scan_with_corrections(folder, corrections_to_nxdict(state))
        with h5py.File(path, "r") as handle:
            integrator = SimpleNamespace(
                database=SimpleNamespace(nxfile=handle, config_target=None)
            )
            mean, compensated = RockingPeakIntegrator._rocking_solid_angle_mean(
                integrator, detector, handle["/61.1"], cnters,
                x=np.full(3, 240.0), y=np.full(3, 300.0),
            )

    assert compensated is applied_at_extraction
    if applied_at_extraction:
        assert mean is not None
        assert np.all(np.isfinite(mean)) and np.all(mean > 0.0)
    else:
        assert mean is None


def test_a_legacy_json_corrections_group_is_still_read():
    """Databases written before the typed layout must keep reducing."""
    import json

    detector = _solid_angle_detector()
    cnters = {"hsize": np.full(2, 20.0), "vsize": np.full(2, 5.0)}
    legacy = {"json": json.dumps({"use_solid_angle": True})}

    with tempfile.TemporaryDirectory() as folder:
        path = _scan_with_corrections(folder, legacy)
        with h5py.File(path, "r") as handle:
            integrator = SimpleNamespace(
                database=SimpleNamespace(nxfile=handle, config_target=None)
            )
            mean, compensated = RockingPeakIntegrator._rocking_solid_angle_mean(
                integrator, detector, handle["/61.1"], cnters,
                x=np.full(2, 240.0), y=np.full(2, 300.0),
            )

    assert compensated is True
    assert mean is not None and np.all(mean > 0.0)

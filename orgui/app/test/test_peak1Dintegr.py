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

from types import SimpleNamespace

import numpy as np

from orgui.app.peak1Dintegr import (
    RockingPeakIntegrator,
    _compute_rocking_integration,
    _trapz_impl,
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
        _stub(detector), cnters, x=np.array([5.0, 300.0, 470.0]), y=np.full(3, 250.0)
    )

    assert applied is True
    np.testing.assert_allclose(acceptance, vsize * detector.PER_ROW, rtol=1e-12)


def test_a_missing_detector_leaves_the_acceptance_out_rather_than_failing():
    """CLI use without a calibration must not lose the integration."""
    cnters = {"vsize": np.array([40.0]), "alpha_pk": np.zeros(1)}

    acceptance, applied = RockingPeakIntegrator._rocking_acceptance(
        _stub(None), cnters, x=np.array([100.0]), y=np.array([200.0])
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

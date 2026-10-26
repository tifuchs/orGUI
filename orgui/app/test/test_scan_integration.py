"""Routing, batching and persistence regressions for scan extraction."""

from collections import Counter
from dataclasses import replace
from types import SimpleNamespace

import h5py
import numpy as np
import pytest
from silx.io.dictdump import dicttonx

from orgui import logger_utils
from orgui.app import scan_integration as integration
from orgui.app.config_data import (
    CorrectionState,
    ROIState,
    roi_from_nxdict,
    roi_to_nxdict,
)
from orgui.app.peak1Dintegr import RockingPeakIntegrator


class ReductionDriver:
    """Run the shared reduction adapter without creating a main window."""

    integrate = RockingPeakIntegrator.integrate
    _integrate_curve_tile = RockingPeakIntegrator._integrate_curve_tile
    get_all_ro_curves = RockingPeakIntegrator.get_all_ro_curves
    _prepareFootprintAction = RockingPeakIntegrator._prepareFootprintAction
    _rocking_normalization = RockingPeakIntegrator._rocking_normalization
    _rocking_acceptance = RockingPeakIntegrator._rocking_acceptance
    _rocking_solid_angle_mean = RockingPeakIntegrator._rocking_solid_angle_mean
    _legacy_rocking_curve_group = staticmethod(
        RockingPeakIntegrator._legacy_rocking_curve_group
    )
    _storedCurveCorrectionRecord = staticmethod(
        RockingPeakIntegrator._storedCurveCorrectionRecord
    )

    def _refreshReductionCorrectionStatus(self):
        pass

    def _stored_detector(self, group):
        return None


class Scan:
    axisname = "th"
    auxillary_counters = ()
    supports_concurrent_read = True

    def __init__(self, frames=7):
        self.axis = np.arange(frames, dtype=float)
        self.reads = Counter()

    def __len__(self):
        return len(self.axis)

    def get_raw_img(self, i):
        """Produce two distinguishable peaks with a constant background."""
        self.reads[i] += 1
        image = np.full((12, 12), 5.0)
        image[2:4, 2:4] = 100 + i
        image[2:4, 7:9] = 200 + i
        return SimpleNamespace(img=image)


class Database:
    compression = None

    def __init__(self, file):
        self.nxfile = file

    def isOpen(self):
        """Return whether the test database is usable."""
        return bool(self.nxfile.id.valid)

    def _requireOpenFile(self):
        return self.nxfile


@pytest.fixture
def context(tmp_path, monkeypatch):
    """Provide a captured test job with an isolated writable HDF5 file."""
    monkeypatch.setattr(logger_utils, "_LOGGING_CONTEXT", "cli")
    scan = Scan()
    file = h5py.File(tmp_path / "batch.h5", "w")
    database = Database(file)
    options = {
        key: False
        for key in (
            "mask",
            "solid_angle",
            "polarization",
            "lorentz",
            "footprint",
            "normalization",
        )
    }
    options.update(
        region={"hsize": 2, "vsize": 2, "left": 1, "right": 1, "top": 0, "bottom": 0},
        advanced={
            "fitted_background": False,
            "detector_inclination": False,
            "project_sample_size": False,
        },
    )
    selector = integration.prepared_selector(
        options, 0, [0, 0, 0], [0, 0, 1], [3, 3], None
    )
    snapshot = SimpleNamespace(
        corrections=CorrectionState(),
        unit_cell=SimpleNamespace(uc_area=1.0),
        ub_calculator=SimpleNamespace(getLambda=lambda: 1.0),
        to_nxdict=lambda **kwargs: {"@NX_class": "NXcollection"},
    )

    def save(data, description):
        """Save test results to the isolated output file."""
        dicttonx(data, file, update_mode="add")
        file.flush()

    captured = integration._CapturedScan(scan, scan.get_raw_img(0))
    ctx = integration.IntegrationContext(
        fscan=captured,
        scanSelector=selector,
        ubcalc=SimpleNamespace(
            detectorCal=SimpleNamespace(detector=SimpleNamespace(shape=(12, 12))),
            chi=0.0,
            phi=0.0,
        ),
        database=database,
        config_snapshot=snapshot,
        activescanname="scan",
        background_image=None,
        numberthreads=2,
        owner=None,
        integrdataPlot=SimpleNamespace(
            getAllCurves=lambda: [], addCurve=lambda *a, **kw: None
        ),
        getMuOm=lambda: (np.zeros(len(scan)), -np.deg2rad(scan.axis)),
        getArmAngles=lambda: (np.zeros(len(scan)), np.zeros(len(scan))),
        getROIloc=None,
        get_detector_mask=lambda shape: None,
        _repair_config_for_image=lambda shape: (False, None, [], []),
        _saveIntegrationResult=save,
    )
    yield ctx, scan, file
    file.close()


def _prepared_line(ctx, x=3, identifier="first"):
    centers = [(slice(x - 1, x + 1), slice(2, 4))]
    rois = {
        "center": centers,
        "left": [(slice(x - 2, x - 1), slice(2, 4))],
        "right": [(slice(x + 1, x + 2), slice(2, 4))],
        "top": [(slice(x - 1, x + 1), slice(2, 2))],
        "bottom": [(slice(x - 1, x + 1), slice(4, 4))],
    }
    geometry = np.zeros((1, len(ctx.fscan), 6))
    geometry[0, :, 0] = x
    geometry[0, :, 3:5] = 0.1
    geometry[0, :, 5] = ctx.fscan.axis
    return {
        "id": identifier,
        "xy": np.array([[x, 3.0]]),
        "rois": rois,
        "name": identifier,
        "geometry": geometry,
        "reflections": {
            "H_0": np.array([x, 0, 0]),
            "H_1": np.array([0, 0, 1]),
            "angles": np.array([[0.0, 0.1, 0.1, 0.0, 0.0, 0.0]]),
            "s_masked": np.array([1.0]),
            "hkl_masked": np.array([[x, 0, 1]]),
        },
    }


@pytest.mark.parametrize("total_flux", [False, True])
@pytest.mark.parametrize("angle_mode", ["counter", "auto", "fixed_incidence"])
def test_shape_extraction_and_bounded_curve_replacement(
    context, monkeypatch, total_flux, angle_mode,
):
    """Both conventions save actual angles and re-read/rebuild shape divisors."""
    from orgui.app.config_data import (
        CURVE_CORRECTIONS_GROUP, curve_correction_record_from_nxdict,
    )
    from orgui.app.sample_interception_config import shape_frame_factors
    from orgui.datautils.xrayutils.corrections.beamprofile import gaussian_profile

    ctx, scan, file = context
    settings = {
        "version": 1, "enabled": True,
        "shape": {"kind": "rectangle", "dimensions_m": [.01, .01]},
        "reference_incidence_deg": .36, "normal_rotation_confirmed": True,
        "azimuth_source": "phi", "horizontal": {
            "analytical": True, "shape": "Top hat", "shape_values": [20000],
        },
    }
    scan.phi = np.linspace(0, 90, len(scan))
    if angle_mode == "auto":
        settings["azimuth_source"] = "auto"
    elif angle_mode == "fixed_incidence":
        settings.update(incidence_source="fixed", fixed_incidence_deg=0.6)
    scan.exposure_time = np.ones(len(scan))
    profile = gaussian_profile(160e-6)
    footprint = SimpleNamespace(beamProfile=lambda: profile, sampleLength=lambda: .01)
    options = ctx.scanSelector.get_integration_options()
    options["footprint"] = options["normalization"] = True
    ctx.config_snapshot.corrections = CorrectionState(
        sample_interception=settings,
        total_flux_calibrated=False if total_flux else None,
    )
    ctx = replace(
        ctx,
        getMuOm=lambda: (np.full(len(scan), np.deg2rad(0.36)), -np.deg2rad(scan.axis)),
        scanSelector=integration.prepared_selector(
            options,
            0,
            [0, 0, 0],
            [0, 0, 1],
            [3, 3],
            footprint,
        ),
    )
    from orgui.app import sample_interception_config

    evaluated_shapes = []

    def evaluate(*args, **kwargs):
        """Count the frame vectors sent to the real overlap calculation."""
        evaluated_shapes.append(np.shape(args[2]))
        return shape_frame_factors(*args, **kwargs)

    monkeypatch.setattr(sample_interception_config, "shape_frame_factors", evaluate)
    first = _prepared_line(ctx)
    for key in first["rois"]:
        first["rois"][key] *= 3
    first["xy"] = np.repeat(first["xy"], 3, axis=0)
    first["geometry"] = np.repeat(first["geometry"], 3, axis=0)
    for key in ("s_masked", "angles", "hkl_masked"):
        first["reflections"][key] = np.repeat(first["reflections"][key], 3, axis=0)
    result = integration.integrate_rocking_scan(
        ctx, lines=[first, _prepared_line(ctx, 8, "second")],
        batch_options=integration.BatchOptions(curves=2),
    )
    assert result["status"] == "success", result
    assert evaluated_shapes == [(len(scan),)]
    for name, rows in (("first", 3), ("second", 1)):
        saved = curve_correction_record_from_nxdict(
            file[f"scan/measurement/{name}"][CURVE_CORRECTIONS_GROUP]
        )
        expected = integration.integration_corrections.frame_correction_policy(
            ctx.fscan, ctx.config_snapshot.corrections, len(scan),
            use_normalization=True, use_illumination=True,
            alpha=saved.alpha, omega=ctx.getMuOm()[1],
            beam_profile=profile, sample_length=.01,
        )
        for field in ("illumination_divisor", "intercepted_fraction",
                      "normalization_divisor"):
            np.testing.assert_array_equal(
                getattr(saved, field), getattr(expected, field)
            )
        for field in ("sample_incidence_rad", "sample_azimuth_rad",
                      "interception_valid", "interception_error_H"):
            np.testing.assert_array_equal(saved.profile_provenance[field],
                                          expected.interception_provenance[field])
        assert saved.base_croibg.shape == (rows, len(scan))
    group = file["scan/measurement/first"]
    record = curve_correction_record_from_nxdict(group[CURVE_CORRECTIONS_GROUP])
    assert record.algorithm == (
        "shape_interception_total_flux_v1"
        if total_flux
        else "shape_interception_legacy_v1"
    )
    np.testing.assert_allclose(
        record.profile_provenance["sample_azimuth_rad"],
        -np.deg2rad(scan.axis) if angle_mode == "auto" else np.deg2rad(scan.phi),
    )
    np.testing.assert_allclose(
        record.profile_provenance["sample_incidence_rad"],
        np.deg2rad(0.6 if angle_mode == "fixed_incidence" else 0.36),
    )
    driver = ReductionDriver()
    driver._currentRoInfo = {"name": group.name, "axisname": "th", "axis": scan.axis}
    driver.database = SimpleNamespace(nxfile=file)
    driver.footprintAction = SimpleNamespace(currentData=lambda: "apply")
    replacement_settings = dict(settings, orientation_deg=45)
    driver.integrationCorrection = SimpleNamespace(
        sampleInterceptionSettings=lambda: replacement_settings,
        beamProfile=lambda: profile,
    )
    driver._prepareFootprintAction(group, lazy=True)
    curve = driver.get_all_ro_curves(rows=slice(0, 1))
    factors = shape_frame_factors(
        replacement_settings, profile,
        np.deg2rad(0.6) if angle_mode == "fixed_incidence" else record.alpha[0:1],
        -np.deg2rad(scan.axis) if angle_mode == "auto" else np.deg2rad(scan.phi),
    )
    divisor = factors[0 if total_flux else 2]
    np.testing.assert_allclose(curve["croibg"], record.base_croibg[0:1]/divisor)
    np.testing.assert_allclose(
        curve["croibg_errors"], np.sqrt(record.base_croibg_variance[0:1]) / divisor
    )


@pytest.mark.parametrize("total_flux", [False, True])
def test_stationary_shape_divisor_is_applied_once(context, monkeypatch, total_flux):
    """Saved stationary signal/errors and record agree under either convention."""
    from orgui.app.config_data import (
        CURVE_CORRECTIONS_GROUP, curve_correction_record_from_nxdict,
    )
    from orgui.app.integration_corrections import corrected_curve_from_record
    from orgui.datautils.xrayutils.corrections.beamprofile import gaussian_profile
    from orgui.app import sample_interception_config

    evaluate = sample_interception_config.shape_frame_factors
    evaluations = []

    def counted(*args, **kwargs):
        """Record overlap evaluations across lines and successive jobs."""
        evaluations.append(np.shape(args[2]))
        return evaluate(*args, **kwargs)

    monkeypatch.setattr(sample_interception_config, "shape_frame_factors", counted)

    ctx, scan, file = context
    scan.phi = np.linspace(0, 90, len(scan))
    scan.exposure_time = np.arange(1., len(scan)+1)
    ctx.config_snapshot.corrections = CorrectionState(
        sample_interception={
            "version": 1, "enabled": True,
            "shape": {"kind": "rectangle", "dimensions_m": [.01, .01]},
            "reference_incidence_deg": .36, "normal_rotation_confirmed": True,
            "azimuth_source": "phi", "horizontal": {
                "analytical": True, "shape": "Top hat", "shape_values": [20000],
            },
        },
        total_flux_calibrated=False if total_flux else None,
    )
    profile = gaussian_profile(160e-6)
    options = ctx.scanSelector.get_integration_options()
    options["footprint"] = options["normalization"] = True
    ctx = replace(
        ctx,
        getMuOm=lambda: (np.full(len(scan), np.deg2rad(0.36)), -np.deg2rad(scan.axis)),
        scanSelector=integration.prepared_selector(
            options,
            0,
            [0, 0, 0],
            [0, 0, 1],
            [3, 3],
            SimpleNamespace(beamProfile=lambda: profile, sampleLength=lambda: 0.01),
        ),
    )

    def geometry(line):
        array = np.zeros((len(scan), 9))
        array[:, 0] = line["H_0"][0]
        array[:, 2] = array[:, 5] = scan.axis
        array[:, 3:5] = .1
        array[:, 6:8] = [line["H_0"][0], 3]
        array[:, -1] = 1
        other = array.copy()
        other[:, -1] = 0
        return array, other

    result = integration.integrate_stationary_scan(
        ctx, lines=[
            {"id": "first", "H_0": [3, 0, 0], "H_1": [0, 0, 1]},
            {"id": "second", "H_0": [8, 0, 0], "H_1": [0, 0, 1]},
        ],
        geometry_provider=geometry,
    )
    assert result["status"] == "success", result
    assert evaluations == [(len(scan),)]
    group = next(
        group for group in file["scan/measurement"].values()
        if CURVE_CORRECTIONS_GROUP in group
    )
    record = curve_correction_record_from_nxdict(group[CURVE_CORRECTIONS_GROUP])
    signal, errors, _, _ = corrected_curve_from_record(record)
    np.testing.assert_allclose(group["counters/ctr_croibg"][()], signal)
    np.testing.assert_allclose(group["counters/ctr_croibg_errors"][()], errors)
    name = "C_illumination" if total_flux else "C_illum_area"
    np.testing.assert_allclose(group["counters"][name][()], record.illumination_divisor)

    # A subsequent job must evaluate fresh settings even with the same context.
    ctx.config_snapshot.corrections.sample_interception["orientation_deg"] = 45
    result = integration.integrate_stationary_scan(
        ctx, lines=[{"id": "third", "H_0": [3, 0, 0], "H_1": [0, 0, 1]}],
        geometry_provider=geometry,
    )
    assert result["status"] == "success", result
    assert evaluations == [(len(scan),), (len(scan),)]


@pytest.mark.parametrize("accelerator", [False, True])
@pytest.mark.parametrize("frames", [1, 3, 64])
def test_rocking_batch_routing_and_singleton_shapes(
    context, monkeypatch, accelerator, frames
):
    """Pin line routing, source reads and singleton shapes across batch sizes."""
    ctx, scan, file = context
    if accelerator and not integration.HAS_ACCEL:
        pytest.skip("Compiled ROI accelerator unavailable")
    monkeypatch.setattr(integration, "HAS_ACCEL", accelerator)
    result = integration.integrate_rocking_scan(
        ctx,
        lines=[_prepared_line(ctx), _prepared_line(ctx, 8, "second")],
        batch_options=integration.BatchOptions(frames=frames, curves=1),
    )
    assert result["status"] == "success", result
    for name, counts in (("first", 95), ("second", 195)):
        group = file[f"scan/measurement/{name}"]
        assert group.attrs["orgui_line_id"] == name
        assert group["rois/s"].shape == (1,)
        assert group["rois/H_0"].shape == (1, 3)
        np.testing.assert_allclose(group["rois/croibg"][0], 4 * (counts + scan.axis))
    assert scan.reads == {i: 1 for i in range(len(scan))}
    assert "_orgui_integration_work" not in file


@pytest.mark.parametrize("mode", ["stationary", "single", "rocking"])
@pytest.mark.parametrize("total_flux", [False, True])
def test_interception_only_for_frames_with_roi_data(
    context, monkeypatch, mode, total_flux,
):
    """Union all lines/intersections, preserving alignment and usable results."""
    import json
    from orgui.app import sample_interception_config
    from orgui.app.config_data import (
        CURVE_CORRECTIONS_GROUP, curve_correction_record_from_nxdict,
    )
    from orgui.app.integration_corrections import corrected_curve_from_record
    from orgui.datautils.xrayutils.corrections.beamprofile import gaussian_profile

    ctx, scan, file = context
    scan.exposure_time = np.arange(1., len(scan) + 1)
    alpha = np.deg2rad(np.linspace(.2, .8, len(scan)))
    omega = -np.deg2rad(scan.axis)
    settings = {
        "version": 1, "enabled": True,
        "shape": {"kind": "rectangle", "dimensions_m": [.01, .01]},
        "offset_m": [.001, .002], "azimuth_source": "auto",
        "horizontal": {
            "analytical": True, "shape": "Top hat", "shape_values": [20000],
        },
    }
    ctx.config_snapshot.corrections = CorrectionState(
        sample_interception=settings,
        total_flux_calibrated=False if total_flux else None,
    )
    profile = gaussian_profile(160e-6)
    options = ctx.scanSelector.get_integration_options()
    options["footprint"] = options["normalization"] = True
    usable = np.zeros((len(scan), 2), dtype=bool)
    usable[[1, 2], 0] = True
    usable[[2, 4, 5], 1] = True

    def geometry(line):
        """Include an S2-only frame and frames needed only by the later line."""
        column = int(line["id"] == "second")
        array = np.zeros((len(scan), 9))
        array[:, 0] = line["H_0"][0]
        array[:, 2] = array[:, 5] = scan.axis
        array[:, 3:5] = .1
        array[:, 6:8] = [line["H_0"][0], 3]
        array[:, -1] = usable[:, column]
        other = array.copy()
        other[:, -1] = 0
        if column:
            array[5, -1] = 0
            other[5, -1] = 1
        return array, other

    lines = [
        {"id": "first", "H_0": [3, 0, 0], "H_1": [0, 0, 1]},
        {"id": "second", "H_0": [8, 0, 0], "H_1": [0, 0, 1]},
    ]
    ctx = replace(
        ctx, getMuOm=lambda: (alpha, omega), getROIloc=lambda: geometry(lines[0]),
        scanSelector=integration.prepared_selector(
            options, 0, [3, 0, 0], [0, 0, 1], [3, 3],
            SimpleNamespace(beamProfile=lambda: profile, sampleLength=lambda: .01),
        ),
    )
    overlap = sample_interception_config.overlap
    evaluated = []

    def counted(*args, **kwargs):
        """Record the actual acquisition angles entering numerical overlap."""
        evaluated.append((args[3].copy(), args[4].copy()))
        return overlap(*args, **kwargs)

    monkeypatch.setattr(sample_interception_config, "overlap", counted)
    if mode == "rocking":
        prepare = integration._prepare_rocking_worker

        def prepare_sparse(*args, **kwargs):
            """Emulate per-frame masked detector ROIs at the worker boundary."""
            worker, *factors = prepare(*args, **kwargs)

            def read(index):
                """Leave valid zero-signal frames in the extracted curves."""
                groups = worker(index)
                for group in groups:
                    group[~usable[index]] = 0
                # Equal peak and background density: zero net signal is data.
                if index == 1:
                    groups[0][0, 0] = groups[0][0, 2] * (
                        groups[0][0, 1] / groups[0][0, 3]
                    )
                return groups

            return read, *factors

        monkeypatch.setattr(integration, "_prepare_rocking_worker", prepare_sparse)
        result = integration.integrate_rocking_scan(
            ctx, lines=[_prepared_line(ctx), _prepared_line(ctx, 8, "second")],
            batch_options=integration.BatchOptions(frames=2, curves=1),
        )
    else:
        result = integration.integrate_stationary_scan(
            ctx, lines=None if mode == "single" else lines,
            geometry_provider=geometry,
        )
    assert result["status"] == "success", result
    required = usable[:, 0] if mode == "single" else usable.any(axis=1)
    assert len(evaluated) == 1
    np.testing.assert_array_equal(evaluated[0][0], alpha[required])
    np.testing.assert_array_equal(evaluated[0][1], omega[required])
    monkeypatch.setattr(sample_interception_config, "overlap", overlap)
    for group in file["scan/measurement"].values():
        if CURVE_CORRECTIONS_GROUP not in group:
            continue
        record = curve_correction_record_from_nxdict(group[CURVE_CORRECTIONS_GROUP])
        reference = integration.integration_corrections.frame_correction_policy(
            ctx.fscan, ctx.config_snapshot.corrections, len(scan),
            use_normalization=True, use_illumination=True, alpha=record.alpha,
            omega=omega, beam_profile=profile, sample_length=.01,
        )
        divisor = np.asarray(record.illumination_divisor)
        mask = np.broadcast_to(required, divisor.shape)
        np.testing.assert_array_equal(
            record.profile_provenance["interception_required"], mask
        )
        assert np.isnan(divisor[~mask]).all()
        np.testing.assert_array_equal(
            divisor[mask], reference.illumination_divisor[mask]
        )
        np.testing.assert_array_equal(
            record.normalization_divisor, reference.normalization_divisor
        )
        alignment = json.loads(record.profile_provenance["sample_interception_json"])
        assert alignment["reference_incidence_deg"] == pytest.approx(.2)
        # Unused factors change, but no saved usable intensity or error changes.
        full_record = replace(
            record, illumination_divisor=reference.illumination_divisor
        )
        for actual, expected in zip(corrected_curve_from_record(record)[:2],
                                    corrected_curve_from_record(full_record)[:2]):
            np.testing.assert_array_equal(actual, expected)


def test_line_configuration_round_trip():
    """Ordered line identities survive typed config serialization."""
    lines = [
        {"id": "rod_a", "H_0": [1, 0, 0], "H_1": [0, 0, 1]},
        {"id": "rod_b", "H_0": [0, 1, 0], "H_1": [0, 0, 1]},
    ]
    restored = roi_from_nxdict(roi_to_nxdict(ROIState(lines={"rocking": lines})))
    assert restored.lines == {"rocking": lines}
    assert roi_from_nxdict({}).lines == {}


def test_shape_policy_cache_cancellation_and_validation(context, monkeypatch):
    """Cancelled overlap is not cached; reuse still validates captured state."""
    ctx, scan, _ = context
    cache = {}
    ctx = replace(ctx, _frame_policy_cache=cache)
    state = CorrectionState(sample_interception={"enabled": True})
    cancelled = [True]
    finished = []
    calls = []
    progress = SimpleNamespace(
        update=lambda value: None,
        wasCanceled=lambda: cancelled[0],
        finish=lambda: finished.append(True),
    )
    monkeypatch.setattr(
        logger_utils, "create_progress_logger", lambda *a, **kw: progress
    )

    def evaluate(*args, **kwargs):
        """Exercise cancellation before returning a reusable policy."""
        calls.append(True)
        kwargs["progress"](1, len(scan))
        return object()

    monkeypatch.setattr(integration.integration_corrections,
                        "frame_correction_policy", evaluate)
    kwargs = dict(use_illumination=True, use_normalization=False,
                  alpha=np.full(len(scan), 0.01))
    with pytest.raises(integration.IntegrationCancelled):
        integration._frame_policy_with_progress(
            ctx, ctx.fscan, state, len(scan), **kwargs
        )
    assert not cache
    cancelled[0] = False
    policy = integration._frame_policy_with_progress(
        ctx, ctx.fscan, state, len(scan), **kwargs
    )
    assert integration._frame_policy_with_progress(
        ctx, ctx.fscan, state, len(scan), **kwargs
    ) is policy
    assert len(calls) == len(finished) == 2

    def invalidate():
        """Reject a stale captured context on a cache hit."""
        raise integration.IntegrationCancelled()

    with pytest.raises(integration.IntegrationCancelled):
        integration._frame_policy_with_progress(
            replace(ctx, validate=invalidate), ctx.fscan, state, len(scan), **kwargs
        )
    assert len(calls) == 2


def test_rocking_writer_keeps_normalization_names_as_shared_metadata(context):
    """Encoded component names must not be appended as ROI curve rows."""
    from orgui.app.database import RockingBatchWriter
    from orgui.app.config_data import (
        CURVE_CORRECTIONS_GROUP,
        CurveCorrectionRecord,
        curve_correction_record_to_nxdict,
        curve_correction_record_from_nxdict,
    )

    ctx, scan, file = context
    components = ("exposure", "primary_monitor:monitor:rate")
    writer = RockingBatchWriter(ctx.database, len(scan), 2)
    for value in (1.0, 2.0):
        base = np.full((1, len(scan)), value)
        record = CurveCorrectionRecord(
            algorithm="framewise_ctr_total_flux_v1",
            output_quantity="rocking_ctr_photon_curve",
            scale_convention="total_flux_calibrated",
            normalization_status="applied",
            normalization_divisor=np.ones(len(scan)),
            normalization_components=components,
            base_croibg=base,
            base_croibg_variance=base,
        )
        writer.append_tile(
            "ctr",
            {"rois": {"croibg": base},
             CURVE_CORRECTIONS_GROUP: curve_correction_record_to_nxdict(record)},
        )
    stored = writer.group["results/ctr"][CURVE_CORRECTIONS_GROUP]
    restored = curve_correction_record_from_nxdict(stored)
    assert restored.normalization_components == components
    np.testing.assert_equal(restored.base_croibg[:, 0], [1.0, 2.0])
    writer.close(state="aborted")


@pytest.mark.parametrize("accelerator", [False, True])
def test_stationary_lines_read_once_and_route(context, monkeypatch, accelerator):
    """Pin stationary multi-line routing and one read per source frame."""
    ctx, scan, file = context
    if accelerator and not integration._roi_sum_accel.HAS_ACCEL_BACKEND:
        pytest.skip("Compiled ROI accelerator unavailable")
    monkeypatch.setattr(integration, "HAS_ACCEL", accelerator)
    lines = [
        {"id": "first", "H_0": [3, 0, 0], "H_1": [0, 0, 1]},
        {"id": "second", "H_0": [8, 0, 0], "H_1": [0, 0, 1]},
    ]

    def geometry(line):
        """Return distinguishable trajectories with an invalid second intersection."""
        array = np.zeros((len(scan), 9))
        array[:, 0] = line["H_0"][0]
        array[:, 2] = array[:, 5] = scan.axis
        array[:, 3:5] = 0.1
        array[:, 6:8] = [line["H_0"][0], 3]
        array[:, -1] = 1
        other = array.copy()
        other[:, -1] = 0
        return array, other

    result = integration.integrate_stationary_scan(
        ctx, lines=lines, geometry_provider=geometry
    )
    assert result["status"] == "success", result
    groups = [
        g for g in file["scan/measurement"].values() if g.attrs.get("orgui_line_id")
    ]
    for identifier, counts in (("first", 95), ("second", 195)):
        group = next(g for g in groups if g.attrs["orgui_line_id"] == identifier)
        np.testing.assert_allclose(
            group["counters/croibg"][()], 4 * (counts + scan.axis)
        )
    assert scan.reads == {i: 1 for i in range(len(scan))}


@pytest.mark.parametrize("total_flux", [False, True])
@pytest.mark.parametrize("lorentz", [False, True])
def test_reduction_tiles_match_full_nested_results(
    context, monkeypatch, total_flux, lorentz
):
    """Bounds, errors and auxiliary outputs use global curve order across tiles."""
    from orgui.app.config_data import (
        ConfigData,
        CurveCorrectionRecord,
        CURVE_CORRECTIONS_GROUP,
        curve_correction_record_to_nxdict,
    )

    ctx, scan, file = context
    entry = _prepared_line(ctx)
    for key in entry["rois"]:
        entry["rois"][key] *= 5
    entry["xy"] = np.repeat(entry["xy"], 5, axis=0)
    entry["geometry"] = np.repeat(entry["geometry"], 5, axis=0)
    for key in ("s_masked", "angles", "hkl_masked"):
        entry["reflections"][key] = np.repeat(entry["reflections"][key], 5, axis=0)
    result = integration.integrate_rocking_scan(
        ctx, lines=[entry], batch_options=integration.BatchOptions(frames=3, curves=2)
    )
    assert result["status"] == "success", result
    group = file["scan/measurement/first"]
    group["rois/alpha"][:] = 1.2
    if total_flux:
        del group[CURVE_CORRECTIONS_GROUP]
        record = CurveCorrectionRecord(
            algorithm="framewise_ctr_total_flux_v1",
            output_quantity="rocking_ctr_photon_curve",
            scale_convention="total_flux_calibrated",
            normalization_status="applied",
            normalization_divisor=np.arange(1.0, 8.0),
            illumination_status="applied",
            illumination_convention="total_flux_H",
            illumination_divisor=np.linspace(0.2, 1.0, 35).reshape(5, 7),
            base_croibg=group["rois/ctr_croibg"][()],
            base_croibg_variance=group["rois/ctr_croibg_errors"][()] ** 2,
            alpha=np.deg2rad(group["rois/alpha"][()]),
            profile_provenance={
                "wavelength_angstrom": 1.0,
                "unitcell_area_angstrom2": 1.0,
            },
        )
        dicttonx(
            {CURVE_CORRECTIONS_GROUP: curve_correction_record_to_nxdict(record)}, group
        )
    roi_info = {
        "sig_1": {"from": np.array([1.0, 2.0, 1.0, 2.0, 3.0]), "to": np.full(5, 4.0)},
        "bg_1": {"from": np.zeros(5), "to": np.ones(5)},
        "bg_2": {"from": np.full(5, 5.0), "to": np.full(5, 6.0)},
    }
    dicttonx({"integration": roi_info}, group, update_mode="add")
    dicttonx(
        {"exposure_time": np.arange(1.0, 8.0)},
        file["scan/auxillary"],
        update_mode="add",
    )
    driver = ReductionDriver()
    driver._currentRoInfo = {"name": group.name, "axisname": "th", "axis": scan.axis}
    driver.footprintAction = SimpleNamespace(currentData=lambda: "keep")
    driver.lorentzButton = SimpleNamespace(isChecked=lambda: lorentz)
    saved = []
    driver.database = SimpleNamespace(
        nxfile=file,
        config_target=None,
        add_nxdict=lambda data, **kwargs: saved.append(data),
    )
    monkeypatch.setattr(ConfigData, "from_gui", lambda gui: ctx.config_snapshot)
    detector_reads = []
    record_reads = []
    read_record = RockingPeakIntegrator._storedCurveCorrectionRecord

    def stored_detector(scangroup):
        """Count geometry rebuilds across complete reduction runs."""
        detector_reads.append(scangroup.name)
        return None

    def stored_record(curvegroup):
        """Count correction parsing across complete reduction runs."""
        record_reads.append(curvegroup.name)
        return read_record(curvegroup)

    monkeypatch.setattr(driver, "_stored_detector", stored_detector)
    monkeypatch.setattr(driver, "_storedCurveCorrectionRecord", stored_record)
    monkeypatch.setattr(
        RockingPeakIntegrator, "_storedCurveCorrectionRecord",
        staticmethod(stored_record),
    )
    for size in (5, 1, 2):
        driver.rocking_batch_options = integration.BatchOptions(curves=size)
        driver.integrate()

    assert len(detector_reads) == (3 if lorentz else 0)
    assert len(record_reads) == 3

    def compare(a, b):
        """Compare every nested reduced array and metadata field."""
        assert a.keys() == b.keys()
        for key in a:
            if isinstance(a[key], dict):
                compare(a[key], b[key])
            elif isinstance(a[key], np.ndarray):
                np.testing.assert_allclose(a[key], b[key])
            else:
                assert a[key] == b[key]

    compare(saved[0], saved[1])
    compare(saved[0], saved[2])


def test_failed_source_does_not_publish_or_leave_staging(context):
    """Abort extraction cleanly when a source image cannot be read."""
    ctx, scan, file = context
    original = scan.get_raw_img

    def broken(index):
        """Inject an I/O failure into an otherwise usable test job."""
        if index == 4:
            raise OSError("source disconnected")
        return original(index)

    scan.get_raw_img = broken
    result = integration.integrate_rocking_scan(
        ctx,
        lines=[_prepared_line(ctx)],
        batch_options=integration.BatchOptions(frames=2, curves=1),
    )
    assert result["status"] == "error"
    assert not result["paths"]
    assert "scan" not in file
    assert "_orgui_integration_work" not in file


def test_cancelled_assembly_publishes_only_completed_lines(context, monkeypatch):
    """Remove a cancelled line before it can be published."""
    ctx, scan, file = context
    progress = SimpleNamespace(
        value=0,
        update=lambda value: setattr(progress, "value", value),
        wasCanceled=lambda: progress.value >= len(scan) + 1,
        finish=lambda: None,
    )
    monkeypatch.setattr(
        logger_utils, "create_progress_logger", lambda *a, **kw: progress
    )
    result = integration.integrate_rocking_scan(
        ctx, lines=[_prepared_line(ctx), _prepared_line(ctx, 8, "second")]
    )
    assert result["status"] == "cancelled"
    assert not result["paths"]
    assert "scan" not in file
    assert "_orgui_integration_work" not in file


def _run_preflight_mode(ctx, mode):
    if mode == "rocking":
        return integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    if mode == "single":
        return integration.integrate_stationary_scan(ctx)
    return integration.integrate_stationary_scan(
        ctx, lines=[{"id": "first", "H_0": [3, 0, 0], "H_1": [0, 0, 1]}],
        geometry_provider=lambda line: pytest.fail("geometry before preflight"),
    )


@pytest.mark.parametrize("mode", ["rocking", "stationary", "single"])
@pytest.mark.parametrize("problem, expected", [
    ("missing_exposure", "exposure_time is required"),
    ("invalid_exposure", "exposure_time must"),
    ("missing_monitor", "no primary monitor"),
    ("invalid_monitor", "monitor must"),
    ("counter_shape", "broadcast"),
    ("missing_flux", "requires total_incident_flux"),
    ("missing_reference", "reference_monitor is required"),
    ("missing_reference_exposure", "reference_exposure is required"),
])
def test_normalization_preflight_precedes_extraction(context, mode, problem, expected):
    """Validate every frame's counters before workers or geometry are started."""
    ctx, scan, file = context
    scan.ic2 = np.ones(len(scan))
    scan.exposure_time = np.ones(len(scan))
    state = CorrectionState(
        primary_monitor="ic2", primary_monitor_kind="rate",
        total_flux_calibrated=False,
    )
    if problem == "missing_exposure":
        del scan.exposure_time
    elif problem == "invalid_exposure":
        scan.exposure_time[-1] = np.nan
    elif problem == "missing_monitor":
        del scan.ic2
    elif problem == "invalid_monitor":
        scan.ic2[-1] = 0.0
    elif problem == "counter_shape":
        scan.ic2 = np.ones(len(scan) - 1)
    else:
        state.total_flux_calibrated = True
        if problem != "missing_flux":
            state.total_incident_flux = 1000.0
        if problem == "missing_reference_exposure":
            state.primary_monitor_kind = "integrated"
            state.monitor_reference_reading = 10.0
    ctx.config_snapshot.corrections = state
    options = ctx.scanSelector.get_integration_options()
    options["normalization"] = True
    ctx.scanSelector.get_integration_options = lambda: options
    reads = scan.reads.copy()
    result = _run_preflight_mode(ctx, mode)
    assert result["status"] == "error"
    assert result["message"].startswith("Correction preflight failed:")
    assert expected in result["message"]
    assert scan.reads == reads
    assert not list(file)


@pytest.mark.parametrize("mode", ["rocking", "stationary", "single"])
@pytest.mark.parametrize("problem, expected", [
    ("horizontal", "horizontal_interception"),
    ("shape", "unsupported sample-interception settings version"),
    ("azimuth", "missing_motor"),
    ("profile", "finite positive peak density"),
])
def test_footprint_preflight_precedes_extraction(context, mode, problem, expected):
    """Footprint settings fail early without performing overlap quadrature."""
    from orgui.datautils.xrayutils.corrections.beamprofile import top_hat_profile

    ctx, scan, file = context
    profile = top_hat_profile(300e-6)
    settings = {}
    if problem in {"shape", "azimuth"}:
        settings = {
            "enabled": True, "version": 1 if problem == "azimuth" else 99,
            "shape": {"kind": "rectangle", "dimensions_m": [.01, .01]},
            "azimuth_source": "missing_motor" if problem == "azimuth" else "auto",
            "normal_rotation_confirmed": True,
            "horizontal": {
                "analytical": True, "shape": "Top hat", "shape_values": [20000],
            },
        }
    elif problem == "profile":
        profile = SimpleNamespace(peak_density=np.nan)
    ctx.config_snapshot.corrections = CorrectionState(
        total_flux_calibrated=False, sample_interception=settings,
    )
    options = ctx.scanSelector.get_integration_options()
    options["footprint"] = True
    ctx = replace(
        ctx, getMuOm=lambda: (np.full(len(scan), .01), np.zeros(len(scan))),
        scanSelector=integration.prepared_selector(
            options, 0, [3, 0, 0], [0, 0, 1], [3, 3],
            SimpleNamespace(beamProfile=lambda: profile, sampleLength=lambda: .01),
        ),
    )
    reads = scan.reads.copy()
    result = _run_preflight_mode(ctx, mode)
    assert result["status"] == "error"
    assert result["message"].startswith("Correction preflight failed:")
    assert expected in result["message"]
    assert scan.reads == reads
    assert not list(file)


@pytest.mark.parametrize("normalize, kind", [(False, "rate"), (True, "integrated")])
def test_preflight_preserves_optional_exposure(context, normalize, kind):
    """Disabled normalization and integrated monitors need no frame exposure."""
    ctx, scan, _ = context
    scan.ic2 = np.ones(len(scan))
    ctx.config_snapshot.corrections = CorrectionState(
        primary_monitor="ic2", primary_monitor_kind=kind,
        total_flux_calibrated=False,
    )
    options = ctx.scanSelector.get_integration_options()
    options["normalization"] = normalize
    ctx.scanSelector.get_integration_options = lambda: options
    result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "success", result


@pytest.mark.parametrize("existing", [False, True])
def test_missing_exposure_preflight_keeps_database_usable(
    context, monkeypatch, existing
):
    """A correction failure does no extraction and leaves the database usable."""
    ctx, scan, file = context
    if existing:
        file["scan/title"] = "previous scan"
        file["scan/measurement/previous"] = [123.0]
    scan.ic2 = np.ones(len(scan))
    ctx.config_snapshot.corrections = CorrectionState(
        primary_monitor="ic2", primary_monitor_kind="rate",
        total_flux_calibrated=False,
    )
    options = ctx.scanSelector.get_integration_options()
    options["normalization"] = True
    ctx.scanSelector.get_integration_options = lambda: options
    refreshed = []
    ctx.database.sync_integration_results = lambda: refreshed.append(True)

    def finish():
        assert "_orgui_integration_work" not in file
        assert refreshed == [True]

    monkeypatch.setattr(
        logger_utils, "create_progress_logger",
        lambda *a, **kw: SimpleNamespace(
            update=lambda value: None, wasCanceled=lambda: False, finish=finish
        ),
    )
    reads_before = scan.reads.copy()
    result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "error"
    assert "exposure_time is required for a rate monitor" in result["traceback"]
    assert not result["paths"]
    assert scan.reads == reads_before
    assert refreshed == []
    if existing:
        assert file["scan/title"].asstr()[()] == "previous scan"
        np.testing.assert_array_equal(file["scan/measurement/previous"], [123.0])
        assert list(file["scan/measurement"]) == ["previous"]
    else:
        assert "scan" not in file

    # The same database supports a retry once the missing metadata is supplied.
    scan.exposure_time = np.ones(len(scan))
    refreshed.clear()
    result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "success"
    assert result["paths"] == ["/scan/measurement/first"]


def test_cancel_after_publication_keeps_completed_line(context, monkeypatch):
    """Cancellation during a later line leaves a previously published line usable."""
    from orgui.app.database import RockingBatchWriter

    ctx, scan, file = context
    cancelled = [False]
    progress = SimpleNamespace(
        update=lambda value: None, wasCanceled=lambda: cancelled[0], finish=lambda: None
    )
    monkeypatch.setattr(
        logger_utils, "create_progress_logger", lambda *a, **kw: progress
    )
    publish = RockingBatchWriter.publish

    def completed(writer, *args, **kwargs):
        """Cancel immediately after a complete line is published."""
        result = publish(writer, *args, **kwargs)
        cancelled[0] = True
        return result

    monkeypatch.setattr(RockingBatchWriter, "publish", completed)
    result = integration.integrate_rocking_scan(
        ctx, lines=[_prepared_line(ctx), _prepared_line(ctx, 8, "second")]
    )
    assert result["status"] == "cancelled"
    assert result["paths"] == ["/scan/measurement/first"]
    np.testing.assert_allclose(
        file["scan/measurement/first/rois/croibg"][0], 4 * (95 + scan.axis)
    )
    assert "second" not in file["scan/measurement"]
    assert "_orgui_integration_work" not in file


def test_failed_publication_removes_new_scan_metadata(context, monkeypatch):
    """A failed result move must roll back the newly published scan shell."""
    ctx, _, file = context
    move = h5py.Group.move

    def fail_result_move(group, source, target):
        if "/results/" in source:
            raise OSError("result move failed")
        return move(group, source, target)

    monkeypatch.setattr(h5py.Group, "move", fail_result_move)
    result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "error"
    assert "result move failed" in result["traceback"]
    assert not result["paths"]
    assert "scan" not in file
    assert "_orgui_integration_work" not in file


@pytest.mark.parametrize("accelerator", [False, True])
def test_mask_and_background_image_across_rocking_tiles(
    context, monkeypatch, accelerator
):
    """Masked signal scaling and background-image subtraction retain their formula."""
    ctx, scan, file = context
    if accelerator and not integration._roi_sum_accel.HAS_ACCEL_BACKEND:
        pytest.skip("Compiled ROI accelerator unavailable")
    monkeypatch.setattr(integration, "HAS_ACCEL", accelerator)
    options = ctx.scanSelector.get_integration_options()
    options["mask"] = True
    mask = np.zeros((12, 12), dtype=bool)
    mask[2, 2] = True
    background = np.full((12, 12), 2.0)
    background[2:4, 2:4] = 10.0
    ctx = replace(
        ctx,
        background_image=background,
        scanSelector=integration.prepared_selector(
            options, 0, [0, 0, 0], [0, 0, 1], [3, 3], None
        ),
        get_detector_mask=lambda shape: mask,
    )
    result = integration.integrate_rocking_scan(
        ctx,
        lines=[_prepared_line(ctx), _prepared_line(ctx, 8, "second")],
        batch_options=integration.BatchOptions(frames=3, curves=1),
    )
    assert result["status"] == "success", result
    np.testing.assert_allclose(
        file["scan/measurement/first/rois/croibg"][0], 4 * (87 + scan.axis)
    )
    np.testing.assert_allclose(
        file["scan/measurement/second/rois/croibg"][0], 4 * (195 + scan.axis)
    )
    assert scan.reads == {i: 1 for i in range(len(scan))}


def test_masked_line_does_not_shift_other_line(context):
    """Dropping masked rows preserves other line IDs and sampling indices."""
    ctx, scan, file = context
    options = ctx.scanSelector.get_integration_options()
    options["mask"] = True
    mask = np.zeros((12, 12), dtype=bool)
    mask[2:4, 2:4] = True
    ctx = replace(
        ctx,
        scanSelector=integration.prepared_selector(
            options, 0, [0, 0, 0], [0, 0, 1], [3, 3], None
        ),
        get_detector_mask=lambda shape: mask,
    )
    second = _prepared_line(ctx, 8, "second")
    second["point_indices"] = np.array([19])
    result = integration.integrate_rocking_scan(
        ctx, lines=[_prepared_line(ctx), second]
    )
    assert result["status"] == "success", result
    assert result["paths"] == ["/scan/measurement/second"]
    assert "first" not in file["scan/measurement"]
    assert file["scan/measurement/second/rois/point_index"][0] == 19


@pytest.mark.parametrize("compression", [None, "gzip", "lzf"])
def test_writer_reuses_compression_and_persists_batch_options(context, compression):
    """Final groups preserve the database's existing compression choice."""
    ctx, scan, file = context
    ctx.database.compression = compression
    result = integration.integrate_rocking_scan(
        ctx,
        lines=[_prepared_line(ctx)],
        batch_options=integration.BatchOptions(frames=3, curves=1, memory_mib=8),
    )
    assert result["status"] == "success", result
    group = file["scan/measurement/first"]
    assert group["rois/croibg"].compression == compression
    assert group.attrs["orgui_frame_batch"] == 3
    assert group.attrs["orgui_memory_mib"] == 8


def test_reduction_cancellation_prevents_any_save(context, monkeypatch):
    """A cancellation cannot reach the old partial-array aggregation branch."""
    ctx, scan, file = context
    result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "success"
    group = file["scan/measurement/first"]
    driver = ReductionDriver()
    driver.database = SimpleNamespace(
        nxfile=file, add_nxdict=lambda *a, **kw: pytest.fail("partial save")
    )
    driver._currentRoInfo = {"name": group.name, "axisname": "th", "axis": scan.axis}
    driver.footprintAction = SimpleNamespace(currentData=lambda: "keep")
    progress = SimpleNamespace(wasCanceled=lambda: True, finish=lambda: None)
    monkeypatch.setattr(
        logger_utils, "create_progress_logger", lambda *a, **kw: progress
    )
    with pytest.raises(integration.IntegrationCancelled):
        driver.integrate()


def test_live_footprint_replacement_reads_only_requested_rows(context, monkeypatch):
    """Live beam settings must not force a full incidence or base-curve read."""
    from orgui.app.config_data import (
        CurveCorrectionRecord,
        CURVE_CORRECTIONS_GROUP,
        curve_correction_record_to_nxdict,
    )
    from orgui.app.integration_corrections import framewise_illumination_divisor
    from orgui.datautils.xrayutils.corrections.beamprofile import top_hat_profile

    ctx, scan, file = context
    base = np.arange(1.0, 36.0).reshape(5, 7)
    alpha = np.linspace(0.01, 0.03, 35).reshape(5, 7)
    record = CurveCorrectionRecord(
        algorithm="framewise_ctr_total_flux_v1",
        output_quantity="rocking_ctr_photon_curve",
        scale_convention="total_flux_relative",
        normalization_status="applied",
        normalization_divisor=np.arange(1.0, 8.0),
        illumination_status="applied",
        illumination_convention="total_flux_H",
        illumination_divisor=np.ones((5, 7)),
        alpha=alpha,
        base_croibg=base,
        base_croibg_variance=base,
    )
    dicttonx(
        {"curve": {CURVE_CORRECTIONS_GROUP: curve_correction_record_to_nxdict(record)}},
        file,
    )
    driver = ReductionDriver()
    driver.database = SimpleNamespace(nxfile=file)
    driver._currentRoInfo = {"name": "/curve", "axisname": "th", "axis": scan.axis}
    driver.footprintAction = SimpleNamespace(currentData=lambda: "apply")
    profile = top_hat_profile(3e-4)
    driver.integrationCorrection = SimpleNamespace(
        beamProfile=lambda: profile,
        sampleLength=lambda: 0.01,
        horizontalInterceptedFraction=lambda: 1.0,
    )
    original = h5py.Dataset.__getitem__
    reads = []

    def read(dataset, key):
        """Reject unbounded reads of per-curve test data."""
        if dataset.ndim == 2:
            assert isinstance(key, slice), (dataset.name, key)
            assert key.start == 1 and key.stop == 3
            reads.append(dataset.name)
        return original(dataset, key)

    monkeypatch.setattr(h5py.Dataset, "__getitem__", read)
    driver._prepareFootprintAction(file["curve"], lazy=True)
    curve = driver.get_all_ro_curves(rows=slice(1, 3))
    illumination = framewise_illumination_divisor(
        alpha[1:3], 0.01, profile, horizontal_fraction=1.0
    )[0]
    np.testing.assert_allclose(
        curve["croibg"], base[1:3] / np.arange(1.0, 8.0) / illumination
    )
    assert reads


def test_plugin_filter_and_failed_writer_setup(context, monkeypatch):
    """Plugin filter mappings work, and setup failures leave no temporary job."""
    from orgui.app.database import RockingBatchWriter

    plugin = pytest.importorskip("hdf5plugin")
    ctx, scan, file = context
    ctx.database.compression = plugin.LZ4()
    result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "success", result
    dataset = file["scan/measurement/first/rois/croibg"]
    assert dataset.id.get_create_plist().get_filter(0)[0] == plugin.LZ4.filter_id
    original = h5py.Group.create_dataset

    def broken(group, name, *a, **kw):
        """Inject an I/O failure into an otherwise usable test job."""
        if name == "counters/2":
            raise OSError("cannot allocate staging dataset")
        return original(group, name, *a, **kw)

    monkeypatch.setattr(h5py.Group, "create_dataset", broken)
    with pytest.raises(OSError, match="cannot allocate"):
        RockingBatchWriter(ctx.database, 7, 1)
    assert "_orgui_integration_work" not in file


def test_main_window_guard_restores_controls_and_blocks_reentry(context):
    """Progress event processing cannot start a second job against captured state."""
    from orgui.app.orGUI import orGUI

    ctx, _, _ = context

    class Control:
        enabled = True

        def isEnabled(self):
            """Return the captured control state."""
            return self.enabled

        def setEnabled(self, value):
            """Update the captured control state."""
            self.enabled = value

    central, menu = Control(), Control()
    window = SimpleNamespace(
        centralWidget=lambda: central,
        menuBar=lambda: menu,
        _integration_context=lambda: ctx,
    )

    def operation(captured):
        """Verify that the outer job rejects reentry during progress events."""
        assert captured is ctx
        assert not central.enabled and not menu.enabled
        second = orGUI._run_scan_integration(window, lambda ctx: pytest.fail("reentry"))
        assert second["status"] == "error"
        return {"status": "success"}

    assert orGUI._run_scan_integration(window, operation)["status"] == "success"
    assert central.enabled and menu.enabled
    assert not window._scan_integration_active
    with pytest.raises(OSError):
        orGUI._run_scan_integration(
            window, lambda ctx: (_ for _ in ()).throw(OSError("failed"))
        )
    assert central.enabled and menu.enabled
    assert not window._scan_integration_active

"""Stationary and rocking ROI integration outside the main window.

Geometry uses radians and reciprocal-lattice units; image/ROI coordinates
use pixels. This module assembles existing correction physics without
changing its definitions. Only the collecting thread writes output data.
"""

from __future__ import annotations

import concurrent.futures
import copy
import logging
import traceback
from numbers import Integral
from dataclasses import dataclass, replace
from types import SimpleNamespace

import numpy as np
from silx.gui import qt

from .. import logger_utils
from . import ROIutils, integration_corrections
from .config_data import (
    CURVE_CORRECTIONS_GROUP,
    CurveCorrectionRecord,
    curve_correction_record_to_nxdict,
)
from ..datautils.xrayutils.corrections import detector as detector_corrections

logger = logging.getLogger(__name__)
try:
    from . import _roi_sum_accel

    HAS_ACCEL = _roi_sum_accel.HAS_ACCEL_BACKEND
except ImportError:
    _roi_sum_accel = None
    HAS_ACCEL = False


def _bounded_results(worker, indices, workers):
    """Yield indexed results with at most ``workers`` live futures."""
    indices = iter(indices)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        pending = {}
        try:
            for _ in range(workers):
                i = next(indices, None)
                if i is None:
                    break
                pending[pool.submit(worker, i)] = i
            while pending:
                done, _ = concurrent.futures.wait(
                    pending, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for future in done:
                    i = pending.pop(future)
                    result = future.result()
                    del future
                    yield i, result
                    del result
                    j = next(indices, None)
                    if j is not None:
                        pending[pool.submit(worker, j)] = j
                done.clear()
        finally:
            for future in pending:
                future.cancel()


@dataclass(frozen=True)
class IntegrationContext:
    """Prepared app inputs and callbacks, with no worker widget access.

    ``scanSelector`` is a value-only adapter, not a Qt widget. ``ubcalc``
    holds calibrated numerical state. Callbacks execute on the collector
    except the captured scan's image reader.
    """

    fscan: object
    scanSelector: object
    ubcalc: object
    database: object
    config_snapshot: object
    activescanname: str
    background_image: object
    numberthreads: int
    owner: object
    integrdataPlot: object
    getMuOm: object
    getArmAngles: object
    getROIloc: object
    get_detector_mask: object
    _repair_config_for_image: object
    _saveIntegrationResult: object
    armFrameGroups: object = None
    validate: object = None

    def intkey(self, coordinates):
        """Return clipped center ROI bounds in pixel (x, y) order."""
        return _intkey(self, coordinates)

    def bkgkeys(self, coordinates):
        """Return clipped background ROI bounds in pixel (x, y) order."""
        return _bkgkeys(self, coordinates)

    def _polarizationArmFactor(self, *args):
        return _polarization_arm_factor(self, *args)


class _Value:
    def __init__(self, value):
        self._value = value

    def value(self):
        """Return the captured numeric value."""
        return self._value

    def isChecked(self):
        """Return the captured switch value."""
        return bool(self._value)


@dataclass(frozen=True)
class BatchOptions:
    """Execution limits; sizes do not affect scientific results.

    :param int frames: Maximum source frames buffered at once.
    :param int curves: Maximum complete rocking curves assembled at once.
    :param float memory_mib: Working-array budget in MiB, excluding detector
        images, shared maps, application state and HDF5 caches.
    """

    frames: int = 64
    curves: int = 32
    memory_mib: float = 256.0

    def sizes(self, frame_count, roi_count):
        """Return batch sizes within the working-array budget."""
        if (
            not isinstance(self.frames, Integral)
            or not isinstance(self.curves, Integral)
            or frame_count < 1
            or roi_count < 1
            or self.frames < 1
            or self.curves < 1
            or not np.isfinite(self.memory_mib)
            or self.memory_mib <= 0
        ):
            raise ValueError("Batch sizes and memory budget must be positive")
        budget = self.memory_mib * 2**20
        # Include in-flight worker results alongside the collector's buffer,
        # and the legacy per-ROI dictionaries alongside their stacked tile.
        frames = min(self.frames, frame_count, int(budget // (roi_count * 192)))
        curves = min(self.curves, roi_count, int(budget // (frame_count * 1280)))
        if frames < 1 or curves < 1:
            raise ValueError("Memory budget cannot hold one frame and curve")
        return frames, curves


class IntegrationCancelled(Exception):
    """Cooperative cancellation that unwinds before saving partial results."""


class _CapturedScan:
    def __init__(self, source, image):
        self._source = source
        self._image = image

    def __len__(self):
        return len(self._source)

    def __getattr__(self, name):
        return getattr(self._source, name)

    def get_raw_img(self, index):
        """Reuse the setup frame; read other source frames on demand."""
        return self._image if index == 0 else self._source.get_raw_img(index)


class _PlotBudget:
    def __init__(self, plot, maximum=30):
        self._plot = plot
        self._remaining = max(0, maximum - len(plot.getAllCurves()))

    def getAllCurves(self):
        """Read plot count on the collecting thread."""
        return self._plot.getAllCurves()

    def addCurve(self, *args, **kwargs):
        """Keep plot payload bounded across all curve tiles."""
        if self._remaining:
            self._plot.addCurve(*args, **kwargs)
            self._remaining -= 1


def normalize_lines(lines):
    """Validate ordered line definitions with vectors in r.l.u.

    :returns: Independent definitions containing id, H_0 and H_1.
    :raises ValueError: For duplicate identifiers or invalid vectors.
    """
    output = []
    identifiers = set()
    for i, line in enumerate(lines):
        identifier = str(line.get("id", f"line_{i + 1}"))
        if not identifier or identifier in identifiers:
            raise ValueError("Line identifiers must be nonempty and unique")
        h0 = np.asarray(line["H_0"], dtype=float)
        h1 = np.asarray(line["H_1"], dtype=float)
        if (
            h0.shape != (3,)
            or h1.shape != (3,)
            or not np.all(np.isfinite([h0, h1]))
            or not np.any(h1)
        ):
            raise ValueError("Lines require finite 3-vectors and nonzero H_1")
        identifiers.add(identifier)
        output.append({"id": identifier, "H_0": h0.tolist(), "H_1": h1.tolist()})
    return output


def _available_name(database, scan_name, stem):
    name, index = stem, 0
    while f"/{scan_name}/measurement/{name}" in database.nxfile:
        name = f"{stem}_{index}"
        index += 1
    return name


def _slice_reflections(reflections, start, stop):
    result = dict(reflections)
    for key in ("angles", "s_masked", "hkl_masked"):
        if key in result:
            result[key] = np.asarray(result[key])[start:stop]
    return result


def _rectangle_key(rectangle):
    return tuple(slice(int(bounds[0]), int(bounds[1])) for bounds in rectangle)


def _validate_context(context):
    if context.validate is not None:
        context.validate()


def _frame_policy_with_progress(context, scan, state, size, **kwargs):
    """Resolve shape overlap with collector-side progress and cancellation."""
    if not (
        kwargs.get("use_illumination")
        and state.sample_interception.get("enabled", False)
    ):
        return integration_corrections.frame_correction_policy(
            scan, state, size, **kwargs
        )
    total = np.size(kwargs["alpha"])
    progress = logger_utils.create_progress_logger(
        context.owner, total, "Calculating sample interception"
    )

    def update(completed, count):
        """Validate captured state whenever GUI progress processes events."""
        progress.update(completed)
        _validate_context(context)
        if progress.wasCanceled():
            raise IntegrationCancelled()
        return True

    try:
        return integration_corrections.frame_correction_policy(
            scan, state, size, progress=update, **kwargs
        )
    finally:
        progress.finish()


def _resolve_mask(context):
    """Resolve a missing requested mask once before any image workers start."""
    options = context.scanSelector.get_integration_options()
    if not options["mask"]:
        return context
    shape = context.ubcalc.detectorCal.detector.shape
    if context.get_detector_mask(shape) is not None:
        return context
    # GUI-only: preserve the existing interactive missing-mask choice.
    if logger_utils.get_logging_context() == "gui":
        answer = qt.QMessageBox.question(
            context.owner,
            "No mask available",
            "No mask was selected with the masking tool.\n"
            "Do you want to continue without mask?",
        )
        if answer != qt.QMessageBox.Yes:
            raise IntegrationCancelled("No mask selected")
    logger.warning("No mask was selected; continuing without a mask")
    options["mask"] = False
    selector = copy.copy(context.scanSelector)
    selector.useMaskBox = _Value(False)
    selector.get_integration_options = lambda: copy.deepcopy(options)
    snapshot = copy.copy(context.config_snapshot)
    snapshot.corrections = replace(snapshot.corrections, use_mask=False)
    return replace(context, scanSelector=selector, config_snapshot=snapshot)


def integrate_stationary_scan(context, *, lines=None, geometry_provider=None):
    """Integrate multiple stationary lines, retaining counters in memory.

    :param lines: Ordered id/H_0/H_1 definitions in r.l.u.; absent definitions
        retain the current single-line or fixed-pixel workflow.
    :param geometry_provider: Callback returning the two (frame, 9)
        trajectories for one line, including pixel coordinates and validity.
    :returns: Existing status dictionary, including partial saved-line status
        on cancellation or errors.
    """
    try:
        return _integrate_stationary_job(
            _resolve_mask(context), lines, geometry_provider
        )
    except IntegrationCancelled:
        return {"status": "cancelled", "message": "No mask selected"}
    except Exception:
        details = traceback.format_exc()
        logger.warning("Stationary integration aborted:\n%s", details)
        return {
            "status": "error",
            "message": "Stationary integration aborted",
            "traceback": details,
        }


def _integrate_stationary_job(context, lines, geometry_provider):
    lines = normalize_lines(lines or [])
    if not lines or context.scanSelector.scanstab.currentIndex() == 1:
        return _integrate_stationary_line(context)
    frame_count = len(context.fscan)
    if not frame_count:
        return {"status": "error", "message": "No source frames"}
    trajectories = []
    for line in lines:
        pair = geometry_provider(line)
        trajectories.append(
            tuple(
                np.broadcast_to(array, (frame_count, array.shape[-1])) for array in pair
            )
        )
    slots = 2 * len(lines)
    rectangles = []
    for frame in range(frame_count):
        regions = [[] for _ in range(5)]
        for pair in trajectories:
            for trajectory in pair:
                valid = bool(trajectory[frame, -1])
                xy = trajectory[frame, 6:8] if valid else np.zeros(2)
                keys = (context.intkey(xy), *context.bkgkeys(xy))
                for index, key in enumerate(keys):
                    bounds = [[part.start, part.stop] for part in key]
                    regions[index].append(bounds if valid else [[0, 0], [0, 0]])
        rectangles.append(
            [np.ascontiguousarray(array, dtype=np.int64) for array in regions]
        )
    first = rectangles[0]
    regions = {
        name: [_rectangle_key(rectangle) for rectangle in first[i]]
        for i, name in enumerate(("center", "left", "right", "top", "bottom"))
    }
    prepared = _prepare_rocking_worker(
        context,
        np.zeros((slots, 2)),
        regions,
        frame_rois=rectangles,
    )
    if isinstance(prepared, dict):
        return prepared
    worker = prepared[0]
    probe = worker(0)
    arrays = [np.zeros((frame_count, slots, 4), dtype=float) for _ in probe]
    for array, counters in zip(arrays, probe):
        array[0] = counters
    del probe
    workers = (
        context.numberthreads
        if getattr(context.fscan, "supports_concurrent_read", True)
        else 1
    )
    progress = logger_utils.create_progress_logger(
        context.owner, frame_count, "Integrating stationary lines"
    )
    results = _bounded_results(worker, range(1, frame_count), max(1, workers))
    processed = 1
    try:
        for frame, counters in results:
            for array, values in zip(arrays, counters):
                array[frame] = values
            processed += 1
            progress.update(processed)
            _validate_context(context)
            if progress.wasCanceled():
                return {
                    "status": "cancelled",
                    "message": "Cancelled during integration",
                }
    finally:
        results.close()
        progress.finish()
    completed = []
    for index, line in enumerate(lines):
        _validate_context(context)

        def counter_rows(index=index):
            """Yield one line's S1/S2 counters in source-frame order."""
            for frame in range(frame_count):
                yield (
                    frame,
                    tuple(array[frame, 2 * index : 2 * index + 2] for array in arrays),
                )

        def save(data, description, line=line):
            """Attach the stable line ID before saving its ordinary results."""
            for key, value in data[context.activescanname]["measurement"].items():
                if not key.startswith("@"):
                    value["@orgui_line_id"] = line["id"]
            return context._saveIntegrationResult(data, description)

        result = _integrate_stationary_line(
            replace(context, _saveIntegrationResult=save),
            counters=counter_rows(),
            geometry=trajectories[index],
            line=line,
        )
        if result["status"] != "success":
            result["completed_lines"] = completed
            return result
        completed.append(line["id"])
    return {"status": "success", "completed_lines": completed}


def integrate_rocking_scan(
    context,
    xylist=None,
    rois=None,
    refldict=None,
    name=None,
    *,
    geometry=None,
    geometry_provider=None,
    lines=None,
    batch_options=None,
):
    """Extract all prepared rocking ROIs once and save bounded curve tiles.

    :param context: Captured :class:`IntegrationContext`.
    :param xylist: Pixel centers for the compatible single-line entry point.
    :param rois: Center/background slice lists, in pixel (x, y) order.
    :param geometry: Optional caller-supplied (ROI, frame, 6) geometry with
        hkl in r.l.u., delta/gamma in radians and scan axis in its stored units.
    :param geometry_provider: Callback taking a bounded array of pixel centers.
    :param lines: Prepared per-line dictionaries; vectors are optional for
        Bragg/fixed-pixel workflows. Each contains xy, rois, reflections, name.
    :param batch_options: Optional :class:`BatchOptions` execution limits.
    :returns: Existing status dictionary plus published result paths.
    """
    from .database import RockingBatchWriter

    if lines is None:
        lines = [
            {
                "xy": xylist,
                "rois": rois,
                "reflections": refldict,
                "name": name,
                "geometry": geometry,
            }
        ]
        for key in ("point_indices", "effective_delta_s"):
            if key in refldict:
                lines[0][key] = refldict[key]
    skipped = [
        entry.get("id", entry["name"]) for entry in lines if not len(entry["xy"])
    ]
    for identifier in skipped:
        logger.warning("Skipping line %s: no detector-valid ROIs", identifier)
    lines = [entry for entry in lines if len(entry["xy"])]
    if not lines:
        return {"status": "error", "message": "No detector-valid rocking ROIs"}
    xy = np.concatenate([entry["xy"] for entry in lines])
    regions = {
        key: [region for entry in lines for region in entry["rois"][key]]
        for key in ("center", "left", "right", "top", "bottom")
    }
    count = len(context.fscan)
    if count < 1:
        return {"status": "error", "message": "No source frames"}
    options = batch_options or BatchOptions()
    batch_frames, batch_curves = options.sizes(count, len(xy))
    workers = (
        context.numberthreads
        if getattr(context.fscan, "supports_concurrent_read", True)
        else 1
    )
    progress = logger_utils.create_progress_logger(
        context.owner, count + len(xy), "Extracting rocking scan batches"
    )
    writer = None
    state = "failed"
    try:
        context = _resolve_mask(context)
        prepared = _prepare_rocking_worker(context, xy, regions)
        if isinstance(prepared, dict):
            return prepared
        worker, polarization, background_polarization = prepared
        probe = worker(0)
        groups = len(probe)
        writer = RockingBatchWriter(
            context.database,
            count,
            len(xy),
            background=groups == 3,
        )
        writer.group.attrs["frame_batch"] = batch_frames
        writer.group.attrs["curve_batch"] = batch_curves
        writer.group.attrs["memory_mib"] = options.memory_mib
        writer.group.attrs["workers"] = max(1, workers)
        logger.info(
            "Rocking batches: %d frames, %d curves, %d ROIs, %d workers",
            batch_frames,
            batch_curves,
            len(xy),
            workers,
        )
        processed = 1
        progress.update(processed)
        for start in range(0, count, batch_frames):
            stop = min(start + batch_frames, count)
            buffer = np.empty((groups, stop - start, len(xy), 4), dtype=float)
            indices = range(max(start, 1), stop) if start == 0 else range(start, stop)
            if start == 0:
                buffer[:, 0] = probe
                del probe
            results = _bounded_results(worker, indices, max(1, workers))
            try:
                for i, result in results:
                    buffer[:, i - start] = result
                    processed += 1
                    progress.update(processed)
                    _validate_context(context)
                    if progress.wasCanceled():
                        raise IntegrationCancelled()
            finally:
                results.close()
            if progress.wasCanceled():
                raise IntegrationCancelled()
            _validate_context(context)
            writer.write_frames(start, buffer)
            del buffer
        del worker, prepared
        writer.group.attrs["state"] = "assembling"
        # Preserve the reader's two-parent scan lookup and scan metadata.
        scan_metadata = {
            context.activescanname: {
                "@NX_class": "NXentry",
                "@orgui_meta": "scan",
                "title": str(getattr(context.fscan, "title", "Rocking scan")),
                "instrument": {
                    "@NX_class": "NXinstrument",
                    "positioners": {
                        "@NX_class": "NXcollection",
                        context.fscan.axisname: context.fscan.axis,
                    },
                },
                "auxillary": {
                    "@NX_class": "NXcollection",
                    **{
                        key: getattr(context.fscan, key)
                        for key in context.fscan.auxillary_counters
                        if getattr(context.fscan, key, None) is not None
                    },
                },
                "configuration": context.config_snapshot.to_nxdict(
                    role="scan", source="scan_import"
                ),
                "measurement": {"@NX_class": "NXentry"},
            }
        }
        error = context._saveIntegrationResult(scan_metadata, "rocking scan metadata")
        if error is not None:
            return error
        offset = 0
        for entry in lines:
            size = len(entry["xy"])
            result_name = _available_name(
                context.database, context.activescanname, entry["name"]
            )
            for start in range(0, size, batch_curves):
                stop = min(start + batch_curves, size)
                if progress.wasCanceled():
                    raise IntegrationCancelled()
                tile_xy = np.asarray(entry["xy"])[start:stop]
                supplied = entry.get("geometry")
                tile_geometry = (
                    np.asarray(supplied[start:stop])
                    if supplied is not None
                    else geometry_provider(tile_xy)
                )
                tile_rois = {
                    key: value[start:stop]
                    for key, value in entry["rois"].items()
                    if key in regions
                }
                counters = writer.read_curves(offset + start, offset + stop)
                selected = np.flatnonzero(np.any(counters[1] > 0, axis=0))
                if not selected.size:
                    progress.update(count + offset + stop)
                    continue
                reflections = _slice_reflections(entry["reflections"], start, stop)
                for key in ("angles", "s_masked", "hkl_masked"):
                    if key in reflections:
                        reflections[key] = reflections[key][selected]
                tile_rois = {
                    key: [value[i] for i in selected]
                    for key, value in tile_rois.items()
                }
                tile = _assemble_rocking_tile(
                    context,
                    tile_xy[selected],
                    tile_rois,
                    tile_geometry[selected],
                    reflections,
                    result_name,
                    [values[:, selected] for values in counters],
                    polarization[offset + start : offset + stop][selected],
                    background_polarization[offset + start : offset + stop][selected],
                    row_offset=start,
                )
                if "status" in tile:
                    return tile
                original = np.asarray(entry.get("point_indices", np.arange(size)))
                tile["rois"]["point_index"] = original[start:stop][selected]
                if "effective_delta_s" in entry:
                    tile["@orgui_effective_delta_s"] = entry["effective_delta_s"]
                if "intersection" in entry:
                    tile["@orgui_intersection"] = entry["intersection"]
                writer.append_tile(result_name, tile)
                progress.update(count + offset + stop)
                _validate_context(context)
                del tile, tile_geometry, counters
            _validate_context(context)
            if progress.wasCanceled():
                raise IntegrationCancelled()
            published = writer.publish(
                context.activescanname, result_name, line_id=entry.get("id")
            )
            if published is None:
                identifier = entry.get("id", entry["name"])
                skipped.append(identifier)
                logger.warning("Skipping line %s: all ROIs are masked out", identifier)
            offset += size
        state = "complete"
        if not writer.published:
            return {
                "status": "error",
                "message": "All rocking ROIs are masked out",
                "paths": [],
                "skipped_lines": skipped,
            }
        return {
            "status": "success",
            "paths": list(writer.published),
            "skipped_lines": skipped,
        }
    except IntegrationCancelled:
        state = "cancelled"
        return {
            "status": "cancelled",
            "message": "Cancelled during rocking integration",
            "paths": [] if writer is None else list(writer.published),
        }
    except Exception:
        message = traceback.format_exc()
        logger.warning("Rocking integration aborted:\n%s", message)
        return {
            "status": "error",
            "message": "Rocking integration aborted",
            "traceback": message,
            "paths": [] if writer is None else list(writer.published),
        }
    finally:
        progress.finish()
        if writer is not None:
            try:
                writer.close(state)
            except Exception:
                logger.warning("Cannot clean integration staging data", exc_info=True)


def prepared_selector(options, mode, h0, h1, xy, footprint):
    """Build a value-only adapter for existing integration calculations.

    :param options: Captured integration options, including pixel ROI sizes.
    :param int mode: Existing stationary/fixed/rocking tab identifier.
    :param h0: Reciprocal-space origin in r.l.u.
    :param h1: Reciprocal-space direction in r.l.u.
    :param xy: Fixed ROI coordinates in detector pixels.
    """
    options = copy.deepcopy(options)
    selector = SimpleNamespace(
        get_integration_options=lambda: copy.deepcopy(options),
        roioptions=SimpleNamespace(get_parameters=lambda: options["advanced"]),
        scanstab=SimpleNamespace(currentIndex=lambda: mode),
        H_0=[_Value(v) for v in h0],
        H_1=[_Value(v) for v in h1],
        xy_static=[_Value(v) for v in xy],
        correctionsDialog=SimpleNamespace(footprintOptions_shared=lambda: footprint),
    )
    for name, value in options["region"].items():
        setattr(selector, name, _Value(value))
    for attr, key in (
        ("useMaskBox", "mask"),
        ("useSolidAngleBox", "solid_angle"),
        ("usePolarizationBox", "polarization"),
    ):
        setattr(selector, attr, _Value(options[key]))
    return selector


def _rocking_arm_snapshot(gamma_arm, delta_arm, curve_shape):
    """Build the unit-tagged per-frame arm fields saved with rocking ROIs."""
    curve_shape = tuple(curve_shape)
    if len(curve_shape) != 2:
        raise ValueError("rocking arm snapshot shape must be (curves, frames)")
    frame_count = curve_shape[1]
    gamma_arm = np.broadcast_to(np.asarray(gamma_arm, dtype=np.float64), (frame_count,))
    delta_arm = np.broadcast_to(np.asarray(delta_arm, dtype=np.float64), (frame_count,))
    return {
        # scan_arm_angles already converted these to true primary-beam
        # scattering angles. Radian storage feeds the geometry API directly.
        "@detector_arm_unit": "rad",
        "@detector_arm_angle_frame": "prim",
        "gamma_arm": np.broadcast_to(gamma_arm, curve_shape).copy(),
        "delta_arm": np.broadcast_to(delta_arm, curve_shape).copy(),
    }


def _correction_region_counters(correction, mask, center, backgrounds):
    """Sum one correction array over the exact valid ROI pixels.

    ``center`` and every entry of ``backgrounds`` use orGUI's ``(x, y)``
    slice order. The returned four counters match the accelerated ROI-sum
    contract: center sum/count and combined-background sum/count.
    """
    correction = np.asarray(correction, dtype=np.float64)
    valid = ~np.asarray(mask, dtype=bool)

    def region_values(region):
        """Return correction sum and valid pixel count for one rectangle."""
        values = correction[region[::-1]]
        region_valid = valid[region[::-1]] & np.isfinite(values)
        return float(np.sum(values[region_valid])), float(np.sum(region_valid))

    center_sum, center_pixels = region_values(center)
    background_sum = 0.0
    background_pixels = 0.0
    for region in backgrounds:
        summed, pixels = region_values(region)
        background_sum += summed
        background_pixels += pixels
    return np.array(
        [center_sum, center_pixels, background_sum, background_pixels],
        dtype=np.float64,
    )


def _warn_masked_peak_scaling(valid_pixels, nominal_pixels, context):
    """Warn once that nominal-area scaling cannot reconstruct masked peaks."""
    valid_pixels = np.asarray(valid_pixels, dtype=np.float64)
    nominal_pixels = np.asarray(nominal_pixels, dtype=np.float64)
    incomplete = (nominal_pixels > 0) & (valid_pixels < nominal_pixels)
    if np.any(incomplete):
        logger.warning(
            "%s contains masked or missing center-ROI pixels. Scaling by "
            "nominal ROI area divided by valid-pixel count preserves a flat "
            "density, but it is not a physical recovery of peak intensity "
            "hidden by detector gaps or masks.",
            context,
        )


def _curve_profile_provenance(state):
    """Flat, self-contained beam-profile provenance for a curve record."""
    result = {
        "analytical": state.beam_shape_analytical,
        "shape": state.beam_shape_name,
        "shape_values": np.asarray(state.beam_shape_values, dtype=np.float64)
        if state.beam_shape_values
        else None,
        "profile_file": state.beam_profile_file,
        "profile_content": state.beam_profile_content,
        "profile_unit": state.beam_profile_unit,
        "profile_center": state.beam_profile_center,
        "profile_offset_um": state.beam_profile_offset_um,
        "profile_positions_m": np.asarray(
            state.beam_profile_positions_m, dtype=np.float64
        )
        if state.beam_profile_positions_m
        else None,
        "profile_density_per_m": np.asarray(
            state.beam_profile_density_per_m, dtype=np.float64
        )
        if state.beam_profile_density_per_m
        else None,
    }
    return {name: value for name, value in result.items() if value is not None}


def _intkey(self, coords):
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
        roioptions["detector_inclination"] or roioptions["project_sample_size"]
    ) and current_mode != 1:
        if roioptions["project_sample_size"]:
            size_exact = ROIutils.calc_corrections(
                coord_restr,
                self.ubcalc.detectorCal,
                np.array([hsize, vsize]),
                roioptions,
                roioptions["detector_inclination"],
                roioptions["factor"],
            )
        else:
            size_exact = ROIutils.calc_corrections(
                coord_restr,
                self.ubcalc.detectorCal,
                np.array([hsize, vsize]),
                None,
                roioptions["detector_inclination"],
                roioptions["factor"],
            )
        hsize = size_exact[0][0]
        vsize = size_exact[0][1]

    vhalfsize = vsize // 2
    hhalfsize = hsize // 2
    fromcoords = np.round(np.asarray(coord_restr) - np.array([hhalfsize, vhalfsize]))
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


def _bkgkeys(self, coords):
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


def _polarization_arm_factor(self, dc, row, column, row_size, column_size, alpha):
    """Per-frame factor moving the polarization onto the real arm position.

    The per-pixel polarization array is built once, from the calibrated
    geometry. That is correct for a detector whose arm does not move, but
    on a scan that drives the arm the same pixel looks in a different
    direction on every frame and the correction comes out far too small --
    10 % at a scattering angle of 18 degrees, 33 % at 30. This returns
    what the already-corrected intensity has to be multiplied by; see
    ``doc/design/ctr_structure_factor_scale.md`` finding F5.

    The factor is exactly ``1.0`` wherever the arm sits at its calibrated
    reference, because both evaluations then use the same geometry, so a
    fixed-arm scan is untouched without needing to be special-cased. A
    constant arm is evaluated once and broadcast, which keeps the cost off
    the common path.

    A rocking or reflectivity scan tracks one region across many frames
    with the arm (and, for a mu scan, alpha itself) different on every
    one -- so nothing is constant, only the region is. That is the case
    :func:`~.corrections.detector.polarization_arm_correction_frames`
    batches into one call instead of one per frame; evaluating the
    per-frame geometry separately for every frame of a mu scan with
    thousands of s points used to make the reduction take minutes. Only
    when the region also changes frame to frame -- a stationary
    integration tracking a rod across the detector -- is there no shared
    region left to batch on, and this falls back to the per-frame loop.

    The arm position comes from :meth:`getArmAngles`, so this shares its
    convention with every other arm consumer in the application: a scan
    that knows nothing about an arm reports zero, which is the calibrated
    reference for the default calibration and therefore leaves the
    correction at one.

    :param dc: The calibrated
        :class:`~orgui.datautils.xrayutils.DetectorCalibration.Detector2D_SXRD`.
    :param row: Region centre row per frame, in pixels (orGUI's ``y``).
    :param column: Region centre column per frame, in pixels (``x``).
    :param row_size: Region height per frame, in pixels.
    :param column_size: Region width per frame, in pixels.
    :param alpha: Incidence angle per frame, in radian.
    :returns: The factor per frame.
    :rtype: numpy.ndarray
    """
    gamma_arm, delta_arm = self.getArmAngles()
    row, column, row_size, column_size, alpha, gamma_arm, delta_arm = (
        np.broadcast_arrays(
            np.asarray(row, dtype=np.float64),
            np.asarray(column, dtype=np.float64),
            np.maximum(np.asarray(row_size, dtype=np.float64), 1.0),
            np.maximum(np.asarray(column_size, dtype=np.float64), 1.0),
            np.asarray(alpha, dtype=np.float64),
            np.asarray(gamma_arm, dtype=np.float64),
            np.asarray(delta_arm, dtype=np.float64),
        )
    )
    if not row.size:
        return np.ones(row.shape, dtype=np.float64)

    def _is_constant(values):
        return np.all(values == values.flat[0])

    if _is_constant(alpha) and _is_constant(gamma_arm) and _is_constant(delta_arm):
        constant = (
            _is_constant(row)
            and _is_constant(column)
            and (_is_constant(row_size) and _is_constant(column_size))
        )
        if constant:
            index = np.unravel_index(0, row.shape)
            factor = detector_corrections.polarization_arm_correction(
                dc,
                row[index],
                column[index],
                row_size[index],
                column_size[index],
                alpha[index],
                float(gamma_arm[index]),
                float(delta_arm[index]),
            )
            return np.full(row.shape, factor)

    region_constant = (
        _is_constant(row)
        and _is_constant(column)
        and _is_constant(row_size)
        and _is_constant(column_size)
    )
    if region_constant:
        index = np.unravel_index(0, row.shape)
        return detector_corrections.polarization_arm_correction_frames(
            dc,
            row[index],
            column[index],
            row_size[index],
            column_size[index],
            alpha.ravel(),
            gamma_arm.ravel(),
            delta_arm.ravel(),
        ).reshape(row.shape)

    factor = np.ones(row.shape, dtype=np.float64)
    for index in np.ndindex(*row.shape):
        factor[index] = detector_corrections.polarization_arm_correction(
            dc,
            row[index],
            column[index],
            row_size[index],
            column_size[index],
            alpha[index],
            float(gamma_arm[index]),
            float(delta_arm[index]),
        )
    return factor


def _integrate_stationary_line(self, *, counters=None, geometry=None, line=None):
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
                "parent": self.owner,
            },
        )
        # print("no images found! %s" % e)
        return {
            "status": "error",
            "message": "No image found in current scan",
            "traceback": traceback.format_exc(),
        }
    if not self.database.isOpen():
        logger.error(
            "Cannot perform stationary scan integration: no database available.",
            extra={
                "title": "Cannot integrate scan",
                "description": "Cannot perform stationary scan integration: no database available.",  # noqa: E501
                "show_dialog": False,
                "dialog_level": logging.WARNING,
                "parent": self.owner,
            },
        )
        # print("No database available")
        return {"status": "error", "message": "No database available"}

    logger.info("Start integration of stationary scan")
    dc = self.ubcalc.detectorCal
    # mu = self.ubcalc.mu

    H_1 = np.array([h.value() for h in self.scanSelector.H_1])
    H_0 = np.array([h.value() for h in self.scanSelector.H_0])
    if line is not None:
        H_0, H_1 = np.asarray(line["H_0"]), np.asarray(line["H_1"])

    vsize = int(self.scanSelector.vsize.value())
    hsize = int(self.scanSelector.hsize.value())
    vsize * hsize  # as set in GUI, no corrections

    imgmask = None

    if self.scanSelector.useMaskBox.isChecked():
        imgmask = self.get_detector_mask(image.img.shape)
        if imgmask is None:
            # GUI-only: compatibility path for direct private helper calls.
            if logger_utils.get_logging_context() == "gui":
                btn = qt.QMessageBox.question(
                    self.owner,
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

    use_solid_angle = self.scanSelector.useSolidAngleBox.isChecked()
    use_polarization = self.scanSelector.usePolarizationBox.isChecked()
    corr = use_solid_angle or use_polarization

    # One definition of the per-pixel factors, shared with the rocking
    # integration and the reciprocal-space reconstruction. It returns None
    # when neither correction is enabled, so that the reconstruction can
    # skip its multiplication; here the array of ones is required, because
    # the branch below that rebuilds it runs only under HAS_ACCEL and the
    # NumPy-only path would otherwise be handed None.
    C_arr = detector_corrections.pixel_factors(
        dc,
        solid_angle=use_solid_angle,
        polarization=use_polarization,
    )
    if C_arr is None:
        C_arr = np.ones(dc.detector.shape, dtype=np.float64)
    P_arr = detector_corrections.pixel_factors(dc, polarization=use_polarization)
    if P_arr is None:
        P_arr = np.ones(dc.detector.shape, dtype=np.float64)

    hkl_del_gam_s1, hkl_del_gam_s2 = (
        geometry if geometry is not None else self.getROIloc()
    )

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
    roi_x_start1_a = np.zeros_like(dataavail, dtype=int)
    roi_x_stop1_a = np.zeros_like(dataavail, dtype=int)
    roi_y_start1_a = np.zeros_like(dataavail, dtype=int)
    roi_y_stop1_a = np.zeros_like(dataavail, dtype=int)

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
    roi_x_start2_a = np.zeros_like(dataavail, dtype=int)
    roi_x_stop2_a = np.zeros_like(dataavail, dtype=int)
    roi_y_start2_a = np.zeros_like(dataavail, dtype=int)
    roi_y_stop2_a = np.zeros_like(dataavail, dtype=int)

    progress = logger_utils.create_progress_logger(
        self.owner, len(self.fscan), "Integrating stationary scan"
    )

    has_bg_img = (
        self.background_image is not None
        and self.background_image.shape == image.img.shape
    )
    roioptions = self.scanSelector.roioptions.get_parameters()
    use_fitted_background = bool(roioptions.get("fitted_background", False))
    fitted_background_order = int(roioptions.get("fitted_background_order", 1))
    if use_fitted_background and not HAS_ACCEL:
        logger.warning(
            "Fitted local background requires the compiled ROI accelerator; "
            "using summed background ROIs instead."
        )
    if use_fitted_background and fitted_background_order >= 1:
        logger.warning(
            "Fitted local background order %d underestimates the "
            "background error: the saved uncertainty still assumes an "
            "unweighted flat-background sample and does not propagate "
            "the polynomial fit's covariance.",
            fitted_background_order,
        )
    repair_enabled, repair, row_gaps, col_gaps = self._repair_config_for_image(
        image.img.shape
    )
    if repair_enabled and use_fitted_background:
        logger.warning(
            "Pixel repair is disabled for fitted local background; "
            "using the original mask for this integration."
        )
        repair_enabled = False

    if imgmask is not None:
        mask = np.ascontiguousarray(imgmask, dtype=bool)
    else:
        mask = np.zeros(image.img.shape, dtype=bool)
    if corr:
        C_arr = np.ascontiguousarray(C_arr, dtype=np.float64)
    else:
        C_arr = np.ones(image.img.shape, dtype=np.float64)
    if not repair_enabled:
        C_arr[mask] = 0.0

    for i in range(len(self.fscan)):
        key = self.intkey(hkl_del_gam_1[i, 6:8])
        croi_key = np.array([[key[0].start, key[0].stop], [key[1].start, key[1].stop]])
        roi_hsize1_a[i] = int(np.abs(np.diff(croi_key[0])[0]))
        roi_vsize1_a[i] = int(np.abs(np.diff(croi_key[1])[0]))
        roi_x_start1_a[i], roi_x_stop1_a[i] = croi_key[0]
        roi_y_start1_a[i], roi_y_stop1_a[i] = croi_key[1]
        key = self.intkey(hkl_del_gam_2[i, 6:8])
        croi_key = np.array([[key[0].start, key[0].stop], [key[1].start, key[1].stop]])
        roi_hsize2_a[i] = int(np.abs(np.diff(croi_key[0])[0]))
        roi_vsize2_a[i] = int(np.abs(np.diff(croi_key[1])[0]))
        roi_x_start2_a[i], roi_x_stop2_a[i] = croi_key[0]
        roi_y_start2_a[i], roi_y_stop2_a[i] = croi_key[1]

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
                    l.append(np.array([[0, 0], [0, 0]]))
                    for l in roi_lists[1:]  # noqa: E741
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
                    l.append(np.array([[0, 0], [0, 0]]))
                    for l in roi_lists[1:]  # noqa: E741
                ]  # will result in zeros, convert to np.nan later
                roi_lists[0].append(np.array([[0, 0], [0, 0]]))
            roi_lists = [
                np.ascontiguousarray(np.stack(l), dtype=np.int64)
                for l in roi_lists  # noqa: E741
            ]
            roi_lists_accel.append(roi_lists)

    if counters is None:
        if HAS_ACCEL:
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
                    if repair_enabled:
                        _roi_sum_accel.processImage_repair_bg_Carr(
                            image,
                            background_image,
                            mask,
                            C_arr,
                            *roi_lists_accel[i],
                            row_gaps,
                            col_gaps,
                            all_counters,
                            Carr_counters,
                            BgImg_counters,
                            repair.max_component_pixels,
                            repair.max_span,
                            repair.radius,
                            repair.min_valid_neighbors,
                        )  # compiled accelerator releases the GIL
                    else:
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
                    elif repair_enabled:
                        _roi_sum_accel.processImage_repair_Carr(
                            image,
                            mask,
                            C_arr,
                            *roi_lists_accel[i],
                            row_gaps,
                            col_gaps,
                            all_counters,
                            Carr_counters,
                            repair.max_component_pixels,
                            repair.max_span,
                            repair.radius,
                            repair.min_valid_neighbors,
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
    try:
        results = (
            counters
            if counters is not None
            else _bounded_results(
                sumImage,
                range(len(self.fscan)),
                max(1, self.numberthreads)
                if getattr(self.fscan, "supports_concurrent_read", True)
                else 1,
            )
        )
        for i, result in results:
            try:
                if has_bg_img:
                    all_counters, Carr_counters, BgImg_counters = result
                    bgimg_croi1_a[i] = BgImg_counters[0, 0]
                    bgimg_cpixel1_a[i] = BgImg_counters[0, 1]
                    bgimg_bgroi1_a[i] = BgImg_counters[0, 2]
                    bgimg_bgpixel1_a[i] = BgImg_counters[0, 3]
                    bgimg_croi2_a[i] = BgImg_counters[1, 0]
                    bgimg_cpixel2_a[i] = BgImg_counters[1, 1]
                    bgimg_bgroi2_a[i] = BgImg_counters[1, 2]
                    bgimg_bgpixel2_a[i] = BgImg_counters[1, 3]

                else:
                    all_counters, Carr_counters = result
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
                progress.update(i + 1)
            except Exception:
                logger.warning("Cannot read image:\n%s", traceback.format_exc())
            if progress.wasCanceled():
                cancelled = True
                break
    finally:
        if hasattr(results, "close"):
            results.close()
        progress.finish()
    if cancelled:
        return {"status": "cancelled", "message": "Cancelled during integration"}

    roi_size1 = roi_hsize1_a * roi_vsize1_a
    roi_size2 = roi_hsize2_a * roi_vsize2_a

    # Polarization-only factors must be accumulated over the same valid
    # pixels as the combined correction. This creates the direct CTR
    # photon branch and avoids estimating it later as mean(S*P)/mean(S).
    polarization_counters = np.zeros((nodatapoints, 2, 4), dtype=np.float64)
    if use_polarization and HAS_ACCEL and repair_enabled:
        for i in range(nodatapoints):
            if not dataavail[i]:
                continue
            dummy_counters = np.zeros((2, 4), dtype=np.float64)
            _roi_sum_accel.processImage_repair_Carr(
                np.ones(image.img.shape, dtype=np.float64),
                mask,
                np.ascontiguousarray(P_arr, dtype=np.float64),
                *roi_lists_accel[i],
                row_gaps,
                col_gaps,
                dummy_counters,
                polarization_counters[i],
                repair.max_component_pixels,
                repair.max_span,
                repair.radius,
                repair.min_valid_neighbors,
            )
    else:
        for i in range(nodatapoints):
            for intersect, hkl_del_gam_current in enumerate(
                (hkl_del_gam_1, hkl_del_gam_2)
            ):
                if not hkl_del_gam_current[i, -1]:
                    continue
                coordinates = hkl_del_gam_current[i, 6:8]
                polarization_counters[i, intersect] = _correction_region_counters(
                    P_arr,
                    mask,
                    self.intkey(coordinates),
                    self.bkgkeys(coordinates),
                )

    P_croi1 = integration_corrections.roi_mean_correction(
        polarization_counters[:, 0, 0], polarization_counters[:, 0, 1]
    )
    P_bgroi1 = integration_corrections.roi_mean_correction(
        polarization_counters[:, 0, 2], polarization_counters[:, 0, 3]
    )
    P_croi2 = integration_corrections.roi_mean_correction(
        polarization_counters[:, 1, 0], polarization_counters[:, 1, 1]
    )
    P_bgroi2 = integration_corrections.roi_mean_correction(
        polarization_counters[:, 1, 2], polarization_counters[:, 1, 3]
    )
    _warn_masked_peak_scaling(
        cpixel1_a,
        np.where(hkl_del_gam_1[:, -1], roi_size1, 0),
        "Stationary S1 extraction",
    )
    _warn_masked_peak_scaling(
        cpixel2_a,
        np.where(hkl_del_gam_2[:, -1], roi_size2, 0),
        "Stationary S2 extraction",
    )

    # Mean correction over the valid pixels of the center ROI. The ROI sum
    # of the correction array must not be rescaled to the nominal ROI area
    # here: croibg already carries that (roi_size / cpixel) factor, so
    # including it again multiplied every corrected intensity by the ROI
    # area. Because the projected ROI size varies over the detector, that
    # scaled two measurements of one rod differently.
    Corr1 = integration_corrections.roi_mean_correction(Corr_croi1_a, Corr_cpixel1_a)
    Corr2 = integration_corrections.roi_mean_correction(Corr_croi2_a, Corr_cpixel2_a)
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
            # NOTE: this error term reuses the unscaled method-1 formula and
            # does not propagate `factor`. It is only exact when the
            # background image is spatially flat across both the center and
            # background ROI footprints; for a structured background image it
            # underestimates or overestimates the true error.
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

    base_signal1 = np.asarray(croibg1_a, dtype=np.float64).copy()
    base_error1 = np.asarray(croibg1_err_a, dtype=np.float64).copy()
    base_signal2 = np.asarray(croibg2_a, dtype=np.float64).copy()
    base_error2 = np.asarray(croibg2_err_a, dtype=np.float64).copy()

    # Geometrical, footprint and normalization corrections. The numerical
    # active area and normalization are intensity divisors. The
    # intercepted-flux fraction is stored only as the numerator from which
    # the active area was constructed, not divided out a second time.
    # Stationary F2 divides by 1/sin(gamma) and has no rod-interception
    # factor (Vlieg 1997, equation 54).
    mu_all, om_all = self.getMuOm()
    alpha_all = np.broadcast_to(
        np.atleast_1d(np.asarray(mu_all, dtype=np.float64)), (nodatapoints,)
    )
    options = self.scanSelector.get_integration_options()
    config_snapshot = self.config_snapshot
    beam_profile = None
    sample_size = None
    if options["footprint"]:
        footprint_dialog = self.scanSelector.correctionsDialog.footprintOptions_shared()  # noqa: E501
        beam_profile = footprint_dialog.beamProfile()
        sample_size = footprint_dialog.sampleLength()  # m
    frame_policy = _frame_policy_with_progress(
        self,
        self.fscan,
        config_snapshot.corrections,
        nodatapoints,
        use_normalization=options["normalization"],
        use_illumination=options["footprint"],
        alpha=alpha_all,
        beam_profile=beam_profile,
        sample_length=sample_size,
    )
    versioned_policy = frame_policy.new_contract or bool(
        frame_policy.interception_provenance
    )
    normalization = (
        frame_policy.normalization_divisor
        if frame_policy.normalization_status == "applied"
        else None
    )
    normalization_applied = list(frame_policy.normalization_components)

    pol_arm1 = np.ones(nodatapoints, dtype=np.float64)
    pol_arm2 = np.ones(nodatapoints, dtype=np.float64)
    if use_polarization:
        # Corr1/Corr2 carry the polarization of the calibrated geometry;
        # move it onto the arm position of each frame (finding F5). This is
        # exactly 1 for a detector whose arm does not move, and it is the
        # reflectivity case -- where the arm follows 2*alpha -- that needs
        # it most.
        pol_arm1 = self._polarizationArmFactor(
            dc, y_coord1_a, x_coord1_a, roi_vsize1_a, roi_hsize1_a, alpha_all
        )
        pol_arm2 = self._polarizationArmFactor(
            dc, y_coord2_a, x_coord2_a, roi_vsize2_a, roi_hsize2_a, alpha_all
        )

    (
        croibg1_a,
        croibg1_err_a,
        ctr_croibg1_a,
        ctr_croibg1_err_a,
    ) = integration_corrections.pixel_correction_branches(
        base_signal1, base_error1, Corr1, P_croi1, pol_arm1
    )
    (
        croibg2_a,
        croibg2_err_a,
        ctr_croibg2_a,
        ctr_croibg2_err_a,
    ) = integration_corrections.pixel_correction_branches(
        base_signal2, base_error2, Corr2, P_croi2, pol_arm2
    )
    combined_croi_factor1 = Corr1 * pol_arm1
    combined_croi_factor2 = Corr2 * pol_arm2
    combined_bgroi_factor1 = (
        integration_corrections.roi_mean_correction(Corr_bgroi1_a, Corr_bgpixel1_a)
        * pol_arm1
    )
    combined_bgroi_factor2 = (
        integration_corrections.roi_mean_correction(Corr_bgroi2_a, Corr_bgpixel2_a)
        * pol_arm2
    )
    polarization_croi_factor1 = P_croi1 * pol_arm1
    polarization_croi_factor2 = P_croi2 * pol_arm2
    polarization_bgroi_factor1 = P_bgroi1 * pol_arm1
    polarization_bgroi_factor2 = P_bgroi2 * pol_arm2
    if croibg1_bgimg_a is not None:
        croibg1_bgimg_a *= combined_croi_factor1
        croibg1_bgimg_err_a *= combined_croi_factor1
    if croibg2_bgimg_a is not None:
        croibg2_bgimg_a *= combined_croi_factor2
        croibg2_bgimg_err_a *= combined_croi_factor2

    correction_factors = []
    for hkl_del_gam in (hkl_del_gam_1, hkl_del_gam_2):
        correction_factors.append(
            integration_corrections.stationary_correction_factors(
                alpha_all,
                hkl_del_gam[:, 3],
                hkl_del_gam[:, 4],
                use_lorentz=options["lorentz"],
                use_footprint=(options["footprint"] and not versioned_policy),
                beam_profile=beam_profile,
                sample_size=sample_size,
                normalization=normalization,
                illumination_divisor=(
                    frame_policy.illumination_divisor
                    if versioned_policy
                    and frame_policy.illumination_status == "applied"
                    else None
                ),
                illumination_convention=frame_policy.illumination_convention
                or "total_flux_H",
            )
        )
    factors1, factors2 = correction_factors
    if factors1.applied:
        logger.info(
            "Stationary scan corrections applied: %s%s",
            ", ".join(factors1.applied),
            (
                f" (normalization: {', '.join(normalization_applied)})"
                if normalization_applied
                else ""
            ),
        )

    # Reversible photon-counting base, after polarization but before the
    # framewise normalization and illumination divisors. Detector solid
    # angle remains confined to the diagnostic intensity branch.
    base_croibg1 = np.asarray(ctr_croibg1_a, dtype=np.float64).copy()
    base_croibg1_variance = np.square(np.asarray(ctr_croibg1_err_a, dtype=np.float64))
    base_croibg2 = np.asarray(ctr_croibg2_a, dtype=np.float64).copy()
    base_croibg2_variance = np.square(np.asarray(ctr_croibg2_err_a, dtype=np.float64))

    croibg1_a, croibg1_err_a = integration_corrections.apply_stationary_corrections(
        croibg1_a, croibg1_err_a, factors1
    )
    croibg2_a, croibg2_err_a = integration_corrections.apply_stationary_corrections(
        croibg2_a, croibg2_err_a, factors2
    )
    ctr_croibg1_a, ctr_croibg1_err_a = (
        integration_corrections.apply_stationary_corrections(
            ctr_croibg1_a, ctr_croibg1_err_a, factors1
        )
    )
    ctr_croibg2_a, ctr_croibg2_err_a = (
        integration_corrections.apply_stationary_corrections(
            ctr_croibg2_a, ctr_croibg2_err_a, factors2
        )
    )
    if croibg1_bgimg_a is not None:
        croibg1_bgimg_a, croibg1_bgimg_err_a = (
            integration_corrections.apply_stationary_corrections(
                croibg1_bgimg_a, croibg1_bgimg_err_a, factors1
            )
        )
    if croibg2_bgimg_a is not None:
        croibg2_bgimg_a, croibg2_bgimg_err_a = (
            integration_corrections.apply_stationary_corrections(
                croibg2_bgimg_a, croibg2_bgimg_err_a, factors2
            )
        )

    F2_hkl1 = F2_hkl1_err = F2_hkl2 = F2_hkl2_err = None
    if options["lorentz"]:
        try:
            common_scale = {
                "wavelength": config_snapshot.ub_calculator.getLambda(),
                "unitcell_area": config_snapshot.unit_cell.uc_area,
            }
            F2_hkl1, F2_hkl1_err = integration_corrections.structure_factor_from_policy(
                ctr_croibg1_a,
                ctr_croibg1_err_a,
                factors1,
                frame_policy,
                **common_scale,
            )
            F2_hkl2, F2_hkl2_err = integration_corrections.structure_factor_from_policy(
                ctr_croibg2_a,
                ctr_croibg2_err_a,
                factors2,
                frame_policy,
                **common_scale,
            )
        except ValueError:
            if not versioned_policy:
                raise
            logger.warning(
                "The explicit CTR normalization is incomplete; saving "
                "diagnostic intensity without labeling it F2_hkl.",
                exc_info=True,
            )

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
        name2 = f"pixloc[{x:.2f} {y:.2f}]_2"  # does not exist, Just for compatibility
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
    gamma_arm_all, delta_arm_all = self.getArmAngles()
    gamma_arm_all = np.broadcast_to(
        np.asarray(gamma_arm_all, dtype=np.float64), (nodatapoints,)
    ).copy()
    delta_arm_all = np.broadcast_to(
        np.asarray(delta_arm_all, dtype=np.float64), (nodatapoints,)
    ).copy()

    suffix = ""
    i = 0

    while (
        self.activescanname + "/measurement/" + name1 + suffix in self.database.nxfile
    ):
        suffix = f"_{i}"
        i += 1
    availname1 = name1 + suffix

    suffix = ""
    i = 0
    while (
        self.activescanname + "/measurement/" + name2 + suffix in self.database.nxfile
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
            "ctr_croibg": ctr_croibg1_a,
            "ctr_croibg_errors": ctr_croibg1_err_a,
            "croibg_bgimg": croibg1_bgimg_a,  # when None, will not create data set
            "croibg_bgimg_errors": croibg1_bgimg_err_a,  # when None, will not create data set  # noqa: E501
            "croi": croi1_a,
            "bgroi": bgroi1_a,
            "croi_pix": cpixel1_a,
            "bgroi_pix": bgpixel1_a,
            "Cfactors_croi": Corr_croi1_a,
            "Cfactors_bgroi": Corr_bgroi1_a,
            "Cfactor_croi": combined_croi_factor1,
            "Cfactor_bgroi": combined_bgroi_factor1,
            "Pfactor_croi": polarization_croi_factor1,
            "Pfactor_bgroi": polarization_bgroi_factor1,
            "bgimg_croi": bgimg_croi1_a,
            "bgimg_bgroi": bgimg_bgroi1_a,
            # None entries do not create a data set, so only the
            # corrections that were enabled are stored.
            "F2_hkl": F2_hkl1,
            "F2_hkl_errors": F2_hkl1_err,
            "C_Lorentz": factors1.get("C_Lorentz"),
            "C_flux_on_sample": factors1.get("C_flux_on_sample"),
            "C_illum_area": factors1.get("C_illum_area"),
            "C_illumination": factors1.get("C_illumination"),
            "C_norm": factors1.get("C_norm"),
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
        "@signal": ("counters/F2_hkl" if F2_hkl1 is not None else "counters/croibg"),
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
            "ctr_croibg": ctr_croibg2_a,
            "ctr_croibg_errors": ctr_croibg2_err_a,
            "croibg_bgimg": croibg2_bgimg_a,
            "croibg_bgimg_errors": croibg2_bgimg_err_a,
            "croi": croi2_a,
            "bgroi": bgroi2_a,
            "croi_pix": cpixel2_a,
            "bgroi_pix": bgpixel2_a,
            "Cfactors_croi": Corr_croi2_a,
            "Cfactors_bgroi": Corr_bgroi2_a,
            "Cfactor_croi": combined_croi_factor2,
            "Cfactor_bgroi": combined_bgroi_factor2,
            "Pfactor_croi": polarization_croi_factor2,
            "Pfactor_bgroi": polarization_bgroi_factor2,
            "bgimg_croi": bgimg_croi2_a,
            "bgimg_bgroi": bgimg_bgroi2_a,
            # None entries do not create a data set, so only the
            # corrections that were enabled are stored.
            "F2_hkl": F2_hkl2,
            "F2_hkl_errors": F2_hkl2_err,
            "C_Lorentz": factors2.get("C_Lorentz"),
            "C_flux_on_sample": factors2.get("C_flux_on_sample"),
            "C_illum_area": factors2.get("C_illum_area"),
            "C_illumination": factors2.get("C_illumination"),
            "C_norm": factors2.get("C_norm"),
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
        "@signal": ("counters/F2_hkl" if F2_hkl2 is not None else "counters/croibg"),
        "@axes": "trajectory/s",
        "@title": self.activescanname + "_" + availname2,
        "@orgui_meta": "roi",
    }

    versioned_policy = frame_policy.new_contract or bool(
        frame_policy.interception_provenance
    )
    profile_provenance = _curve_profile_provenance(config_snapshot.corrections)
    profile_provenance.update(frame_policy.interception_provenance)
    profile_provenance.update(
        {
            "wavelength_angstrom": config_snapshot.ub_calculator.getLambda(),
            "unitcell_area_angstrom2": config_snapshot.unit_cell.uc_area,
            "detector_efficiency_assumed": 1.0,
            "external_transmission_assumed": 1.0,
        }
    )

    def versioned_stationary_curve(
        factors,
        base_croibg,
        base_croibg_variance,
        croi,
        bgroi,
        combined_croi,
        combined_bgroi,
        polarization_croi,
        polarization_bgroi,
        x,
        y,
        width,
        height,
        x_start,
        x_stop,
        y_start,
        y_stop,
    ):
        """Build the non-legacy sibling branch for one trajectory."""
        record = CurveCorrectionRecord(
            algorithm=(
                (
                    "shape_interception_total_flux_v1"
                    if frame_policy.interception_provenance
                    else "framewise_ctr_total_flux_v1"
                )
                if frame_policy.new_contract
                else (
                    "shape_interception_legacy_v1"
                    if versioned_policy
                    else "legacy_stationary_roi_v2"
                )
            ),
            output_quantity=(
                "stationary_ctr_photon_curve"
                if frame_policy.new_contract
                else (
                    "stationary_shape_density_curve"
                    if versioned_policy
                    else "stationary_roi_intensity"
                )
            ),
            scale_convention=(
                frame_policy.scale_convention
                if versioned_policy
                else (
                    "legacy_density_area" if options["footprint"] else "legacy_relative"
                )
            ),
            normalization_status=frame_policy.normalization_status,
            illumination_status=frame_policy.illumination_status,
            pixel_correction_status=("applied" if corr else "not_applied"),
            normalization_divisor=frame_policy.normalization_divisor,
            normalization_unit=frame_policy.normalization_unit,
            normalization_components=frame_policy.normalization_components,
            illumination_divisor=frame_policy.illumination_divisor,
            illumination_convention=frame_policy.illumination_convention,
            vertical_intercepted_fraction=(frame_policy.vertical_intercepted_fraction),
            horizontal_intercepted_fraction=(
                frame_policy.horizontal_intercepted_fraction
            ),
            intercepted_fraction=frame_policy.intercepted_fraction,
            alpha=alpha_all,
            base_croi=croi,
            base_croi_variance=croi,
            base_bgroi=bgroi,
            base_bgroi_variance=bgroi,
            base_croibg=base_croibg,
            base_croibg_variance=base_croibg_variance,
            combined_croi_factor=combined_croi,
            combined_bgroi_factor=combined_bgroi,
            polarization_croi_factor=polarization_croi,
            polarization_bgroi_factor=polarization_bgroi,
            gamma_arm=gamma_arm_all,
            delta_arm=delta_arm_all,
            roi_x=x,
            roi_y=y,
            roi_width=width,
            roi_height=height,
            roi_x_start=x_start,
            roi_x_stop=x_stop,
            roi_y_start=y_start,
            roi_y_stop=y_stop,
            lorentz_mode=("stationary" if options["lorentz"] else None),
            profile_provenance=profile_provenance,
        )
        return curve_correction_record_to_nxdict(record)

    datas1[CURVE_CORRECTIONS_GROUP] = versioned_stationary_curve(
        factors1,
        base_croibg1,
        base_croibg1_variance,
        croi1_a,
        bgroi1_a,
        combined_croi_factor1,
        combined_bgroi_factor1,
        polarization_croi_factor1,
        polarization_bgroi_factor1,
        x_coord1_a,
        y_coord1_a,
        roi_hsize1_a,
        roi_vsize1_a,
        roi_x_start1_a,
        roi_x_stop1_a,
        roi_y_start1_a,
        roi_y_stop1_a,
    )
    datas2[CURVE_CORRECTIONS_GROUP] = versioned_stationary_curve(
        factors2,
        base_croibg2,
        base_croibg2_variance,
        croi2_a,
        bgroi2_a,
        combined_croi_factor2,
        combined_bgroi_factor2,
        polarization_croi_factor2,
        polarization_bgroi_factor2,
        x_coord2_a,
        y_coord2_a,
        roi_hsize2_a,
        roi_vsize2_a,
        roi_x_start2_a,
        roi_x_stop2_a,
        roi_y_start2_a,
        roi_y_stop2_a,
    )
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
            "configuration": config_snapshot.to_nxdict(
                role="scan", source="scan_import"
            ),
            "@NX_class": "NXentry",
            "@default": "measurement/%s" % (availname1 if defaultS1 else availname2),
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
        data[self.activescanname]["measurement"][availname1]["configuration"] = (
            config_snapshot.to_nxdict(role="integration", source="integration_save")
        )
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
        data[self.activescanname]["measurement"][availname2]["configuration"] = (
            config_snapshot.to_nxdict(role="integration", source="integration_save")
        )

        names_to_log += availname2

    error = self._saveIntegrationResult(data, f"the scan integration {names_to_log}")
    if error is not None:
        return error
    logger.info(f"stationary scan integrated and saved with name(s) {names_to_log}")
    return {"status": "success"}


def _prepare_rocking_worker(self, xylist, rois, *, frame_rois=None):
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
                "parent": self.owner,
            },
        )
        return {
            "status": "error",
            "message": "No image found in current scan",
            "traceback": traceback.format_exc(),
        }
    if not self.database.isOpen():
        logger.error(
            "Cannot integrate scan: No database available.",
            extra={
                "title": "Cannot integrate scan",
                "description": "Cannot integrate scan: No database available.",
                "show_dialog": True,
                "dialog_level": logging.WARNING,
                "parent": self.owner,
            },
        )
        raise ValueError("No database available")
    dc = self.ubcalc.detectorCal

    imgmask = None

    if self.scanSelector.useMaskBox.isChecked():
        imgmask = self.get_detector_mask(image.img.shape)
        if imgmask is None:
            # GUI-only: compatibility path for direct private helper calls.
            if logger_utils.get_logging_context() == "gui":
                btn = qt.QMessageBox.question(
                    self.owner,
                    "No mask available",
                    """No mask was selected with the masking tool.
    Do you want to continue without mask?""",
                )
                if btn != qt.QMessageBox.Yes:
                    return {
                        "status": "cancelled",
                        "message": "Reason: no mask selected",
                    }
            logger.warning("No mask was selected with the masking tool.")

    use_solid_angle = self.scanSelector.useSolidAngleBox.isChecked()
    use_polarization = self.scanSelector.usePolarizationBox.isChecked()
    corr = use_solid_angle or use_polarization
    mask = (
        np.ascontiguousarray(imgmask, dtype=bool)
        if imgmask is not None
        else np.zeros(image.img.shape, dtype=bool)
    )

    # One definition of the per-pixel factors, shared with the rocking
    # integration and the reciprocal-space reconstruction. It returns None
    # when neither correction is enabled, so that the reconstruction can
    # skip its multiplication; here the array of ones is required, because
    # the branch below that rebuilds it runs only under HAS_ACCEL and the
    # NumPy-only path would otherwise be handed None.
    C_arr = detector_corrections.pixel_factors(
        dc,
        solid_angle=use_solid_angle,
        polarization=use_polarization,
    )
    if C_arr is None:
        C_arr = np.ones(dc.detector.shape, dtype=np.float64)
    P_arr = detector_corrections.pixel_factors(dc, polarization=use_polarization)
    if P_arr is None:
        P_arr = np.ones(dc.detector.shape, dtype=np.float64)

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

    background_image = self.background_image
    has_bg_img = False
    roioptions = self.scanSelector.roioptions.get_parameters()
    use_fitted_background = bool(roioptions.get("fitted_background", False))
    fitted_background_order = int(roioptions.get("fitted_background_order", 1))
    if use_fitted_background and not HAS_ACCEL:
        logger.warning(
            "Fitted local background requires the compiled ROI accelerator; "
            "using summed background ROIs instead."
        )
    if use_fitted_background and fitted_background_order >= 1:
        logger.warning(
            "Fitted local background order %d underestimates the "
            "background error: the saved uncertainty still assumes an "
            "unweighted flat-background sample and does not propagate "
            "the polynomial fit's covariance.",
            fitted_background_order,
        )
    repair_enabled = False
    roi_lists_accel = None
    if HAS_ACCEL:
        repair_enabled, repair, row_gaps, col_gaps = self._repair_config_for_image(
            image.img.shape
        )
        if repair_enabled and use_fitted_background:
            logger.warning(
                "Pixel repair is disabled for fitted local background; "
                "using the original mask for this integration."
            )
            repair_enabled = False
        if corr:
            C_arr = np.ascontiguousarray(C_arr, dtype=np.float64)
        else:
            C_arr = np.ones(image.img.shape, dtype=np.float64)
        if not repair_enabled:
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
        if background_image is not None and background_image.shape == image.img.shape:
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
                if repair_enabled:
                    _roi_sum_accel.processImage_repair_bg_Carr(
                        image,
                        bg,
                        mask,
                        C_arr,
                        *(roi_lists_accel if frame_rois is None else frame_rois[i]),
                        row_gaps,
                        col_gaps,
                        all_counters,
                        Carr_counters,
                        BgImg_counters,
                        repair.max_component_pixels,
                        repair.max_span,
                        repair.radius,
                        repair.min_valid_neighbors,
                    )  # compiled accelerator releases the GIL
                else:
                    _roi_sum_accel.processImage_bg_Carr(
                        image,
                        bg,
                        mask,
                        C_arr,
                        *(roi_lists_accel if frame_rois is None else frame_rois[i]),
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
                        *(roi_lists_accel if frame_rois is None else frame_rois[i]),
                        all_counters,
                        Carr_counters,
                        fitted_background_order,
                    )  # compiled accelerator releases the GIL
                elif repair_enabled:
                    _roi_sum_accel.processImage_repair_Carr(
                        image,
                        mask,
                        C_arr,
                        *(roi_lists_accel if frame_rois is None else frame_rois[i]),
                        row_gaps,
                        col_gaps,
                        all_counters,
                        Carr_counters,
                        repair.max_component_pixels,
                        repair.max_span,
                        repair.radius,
                        repair.min_valid_neighbors,
                    )  # compiled accelerator releases the GIL
                else:
                    _roi_sum_accel.processImage_Carr(
                        image,
                        mask,
                        C_arr,
                        *(roi_lists_accel if frame_rois is None else frame_rois[i]),
                        all_counters,
                        Carr_counters,
                    )  # compiled accelerator releases the GIL
                return all_counters, Carr_counters

    else:
        has_bg_img = (
            background_image is not None and background_image.shape == image.img.shape
        )

        def sumImage(i):
            """CLI-safe worker: read and integrate one image without acceleration."""  # noqa: E501
            image = self.fscan.get_raw_img(i).img.astype(
                np.float64, order="C", copy=True
            )
            if imgmask is not None:
                image[imgmask] = np.nan
                pixelavail = (~imgmask).astype(np.float64)
            else:
                pixelavail = np.ones_like(image)

            # Keep raw counts and correction factors separate, matching
            # the accelerated path.  Applying C_arr to image here would
            # corrupt the Poisson counts used for error propagation.
            all_counters = np.zeros((xylist.shape[0], 4), dtype=np.float64)
            Carr_counters = np.zeros_like(all_counters)
            correction_image = C_arr.copy()
            if imgmask is not None:
                correction_image[imgmask] = np.nan

            if has_bg_img:
                background = background_image.astype(np.float64, order="C", copy=True)
                if imgmask is not None:
                    background[imgmask] = np.nan
                BgImg_counters = np.zeros_like(all_counters)

            for crnr in range(xylist.shape[0]):
                # set ROI (moved to rocking-function)

                # get roi
                key = (
                    rois["center"][crnr]
                    if frame_rois is None
                    else _rectangle_key(frame_rois[i][0][crnr])
                )
                bgkey = [
                    rois["left"][crnr],
                    rois["right"][crnr],
                    rois["top"][crnr],
                    rois["bottom"][crnr],
                ]
                if frame_rois is not None:
                    bgkey = [_rectangle_key(array[crnr]) for array in frame_rois[i][1:]]
                # fill counters
                all_counters[crnr] = fill_counters(image, pixelavail, key, bgkey)
                Carr_counters[crnr] = fill_counters(
                    correction_image, pixelavail, key, bgkey
                )
                if has_bg_img:
                    BgImg_counters[crnr] = fill_counters(
                        background, pixelavail, key, bgkey
                    )

            if has_bg_img:
                return all_counters, Carr_counters, BgImg_counters
            return all_counters, Carr_counters

    if frame_rois is not None:
        return sumImage, None, None
    polarization_counters = np.zeros((xylist.shape[0], 4), dtype=np.float64)
    if use_polarization and HAS_ACCEL and repair_enabled:
        dummy_counters = np.zeros_like(polarization_counters)
        _roi_sum_accel.processImage_repair_Carr(
            np.ones(image.img.shape, dtype=np.float64),
            mask,
            np.ascontiguousarray(P_arr, dtype=np.float64),
            *roi_lists_accel,
            row_gaps,
            col_gaps,
            dummy_counters,
            polarization_counters,
            repair.max_component_pixels,
            repair.max_span,
            repair.radius,
            repair.min_valid_neighbors,
        )
    else:
        for crnr in range(xylist.shape[0]):
            polarization_counters[crnr] = _correction_region_counters(
                P_arr,
                mask,
                rois["center"][crnr],
                [
                    rois["left"][crnr],
                    rois["right"][crnr],
                    rois["top"][crnr],
                    rois["bottom"][crnr],
                ],
            )
    P_croi = integration_corrections.roi_mean_correction(
        polarization_counters[:, 0], polarization_counters[:, 1]
    )
    P_bgroi = integration_corrections.roi_mean_correction(
        polarization_counters[:, 2], polarization_counters[:, 3]
    )

    return sumImage, P_croi, P_bgroi


def _assemble_rocking_tile(
    self,
    xylist,
    rois,
    hkl_del_gam,
    refldict,
    name,
    counters,
    P_croi,
    P_bgroi,
    *,
    row_offset=0,
):
    """Construct existing output quantities for a bounded group of curves."""
    dc = self.ubcalc.detectorCal
    options = self.scanSelector.get_integration_options()
    use_polarization = options["polarization"]
    corr = options["solid_angle"] or use_polarization
    rocking_roi_edges = {
        edge: np.asarray(
            [getattr(region[axis], bound) for region in rois["center"]], dtype=np.int64
        )
        for edge, axis, bound in (
            ("x_start", 0, "start"),
            ("x_stop", 0, "stop"),
            ("y_start", 1, "start"),
            ("y_stop", 1, "stop"),
        )
    }
    croi1_all = counters[0]
    cpixel1_all = counters[1]
    bgroi1_all = counters[2]
    bgpixel1_all = counters[3]
    Corr_croi1_all = counters[4]
    Corr_cpixel1_all = counters[5]
    Corr_bgroi1_all = counters[6]
    Corr_bgpixel1_all = counters[7]
    bgimg_croi1_all = counters[8]
    bgimg_cpixel1_all = counters[9]
    bgimg_bgroi1_all = counters[10]
    bgimg_bgpixel1_all = counters[11]
    nominal_rocking_pixels = np.asarray(
        [
            (region[0].stop - region[0].start) * (region[1].stop - region[1].start)
            for region in rois["center"]
        ],
        dtype=np.float64,
    )
    _warn_masked_peak_scaling(
        cpixel1_all,
        np.broadcast_to(nominal_rocking_pixels, cpixel1_all.shape),
        "Rocking extraction",
    )

    currentPlotCount = len(self.integrdataPlot.getAllCurves())
    numberOfNewPlots = xylist.shape[0]
    maxAmountOfPlots = 30
    plotOnlyNth = (
        numberOfNewPlots // max((maxAmountOfPlots - currentPlotCount), 1)
    ) + 1

    # print('Number of integration curves: ' + str(numberOfPlots))
    # print('We can plot every ' + str(plotOnlyNth) + '-th curve.' )

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
    gamma_arm, delta_arm = self.getArmAngles()
    gamma_arm = np.broadcast_to(
        np.asarray(gamma_arm, dtype=np.float64), (len(self.fscan),)
    ).copy()
    delta_arm = np.broadcast_to(
        np.asarray(delta_arm, dtype=np.float64), (len(self.fscan),)
    ).copy()

    config_snapshot = self.config_snapshot
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
            "configuration": config_snapshot.to_nxdict(
                role="scan", source="scan_import"
            ),
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
        roi_size = (roi_d[0].stop - roi_d[0].start) * (roi_d[1].stop - roi_d[1].start)

        hkl_del_gam_1 = hkl_del_gam[d]

        croi1_a = croi1_all[..., d]
        cpixel1_a = cpixel1_all[..., d]
        bgroi1_a = bgroi1_all[..., d]
        bgpixel1_a = bgpixel1_all[..., d]

        Corr_croi1_a = Corr_croi1_all[..., d]
        Corr_cpixel1_a = Corr_cpixel1_all[..., d]
        # NOTE: Corr_bgroi1_a is stored as Cfactors_bgroi but is not applied.
        # The background is corrected with the center ROI's mean factor; see
        # the audit note in the changelog.
        Corr_bgroi1_a = Corr_bgroi1_all[..., d]

        bgimg_croi1_a = bgimg_croi1_all[..., d]
        bgimg_cpixel1_a = bgimg_cpixel1_all[..., d]
        bgimg_bgroi1_a = bgimg_bgroi1_all[..., d]
        bgimg_bgpixel1_a = bgimg_bgpixel1_all[..., d]

        # Mean correction over the valid pixels of the center ROI; see
        # the stationary path for why the ROI area must not appear here.
        Corr1 = integration_corrections.roi_mean_correction(
            Corr_croi1_a, Corr_cpixel1_a
        )

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
                # NOTE: this error term reuses the unscaled method-1 formula and
                # does not propagate `factor`. It is only exact when the
                # background image is spatially flat across both the center and
                # background ROI footprints; for a structured background image it
                # underestimates or overestimates the true error.
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

        base_croibg1 = np.asarray(croibg1_a, dtype=np.float64).copy()
        base_croibg1_err = np.asarray(croibg1_err_a, dtype=np.float64).copy()
        pol_arm1 = np.ones_like(base_croibg1)
        if use_polarization:
            # Corr1 carries the polarization of the calibrated geometry;
            # move it onto the arm position of each frame (finding F5).
            # Exactly 1 for a detector whose arm does not move.
            pol_arm1 = self._polarizationArmFactor(
                dc,
                xylist[d][1],
                xylist[d][0],
                roi_d[1].stop - roi_d[1].start,
                roi_d[0].stop - roi_d[0].start,
                mu,
            )

        (
            croibg1_a,
            croibg1_err_a,
            ctr_croibg1_a,
            ctr_croibg1_err_a,
        ) = integration_corrections.pixel_correction_branches(
            base_croibg1,
            base_croibg1_err,
            Corr1,
            P_croi[d],
            pol_arm1,
        )
        combined_croi_factor = Corr1 * pol_arm1
        combined_bgroi_factor = (
            integration_corrections.roi_mean_correction(
                Corr_bgroi1_a, Corr_bgpixel1_all[..., d]
            )
            * pol_arm1
        )
        polarization_croi_factor = P_croi[d] * pol_arm1
        polarization_bgroi_factor = P_bgroi[d] * pol_arm1
        if croibg1_bgimg_a is not None:
            croibg1_bgimg_a = croibg1_bgimg_a * combined_croi_factor
            croibg1_bgimg_err_a = croibg1_bgimg_err_a * combined_croi_factor

        rod_mask1 = np.isfinite(croibg1_a)

        axis_masked = hkl_del_gam_1[:, 5][rod_mask1]

        croibg1_a_masked = croibg1_a[rod_mask1]

        croibg1_err_a_masked = croibg1_err_a[rod_mask1]

        # save data

        x, y = xylist[d]
        name1 = f"rocking_{d + row_offset}"
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
                "ctr_croibg": ctr_croibg1_a,
                "ctr_croibg_errors": ctr_croibg1_err_a,
                "croibg_bgimg": croibg1_bgimg_a,  # when None, will not create data set  # noqa: E501
                "croibg_bgimg_errors": croibg1_bgimg_err_a,  # when None, will not create data set  # noqa: E501
                "croi": croi1_a,
                "bgroi": bgroi1_a,
                "croi_pix": cpixel1_a,
                "bgroi_pix": bgpixel1_a,
                "Cfactors_croi": Corr_croi1_a,
                "Cfactors_bgroi": Corr_bgroi1_a,
                "Cfactor_croi": combined_croi_factor,
                "Cfactor_bgroi": combined_bgroi_factor,
                "Pfactor_croi": polarization_croi_factor,
                "Pfactor_bgroi": polarization_bgroi_factor,
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
            "configuration": config_snapshot.to_nxdict(
                role="integration", source="integration_save"
            ),
        }

        data[self.activescanname]["measurement"][name]["rois"]["@default"] = availname1
        if np.any(cpixel1_a > 0.0):
            data[self.activescanname]["measurement"][name]["rois"][availname1] = datas1
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
                    "configuration": config_snapshot.to_nxdict(
                        role="integration", source="integration_save"
                    ),
                },
            },
            "title": f"{title}",
            "configuration": config_snapshot.to_nxdict(
                role="scan", source="scan_import"
            ),
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
    ctr_croibg = []
    ctr_croibg_errors = []
    croi = []
    bgroi = []
    croi_pix = []
    bgroi_pix = []
    croibg_bgimg = []
    croibg_bgimg_errors = []
    Cfactors_croi = []
    Cfactors_bgroi = []
    Cfactor_croi = []
    Cfactor_bgroi = []
    Pfactor_croi = []
    Pfactor_bgroi = []
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
            ctr_croibg.append(dsc["counters"]["ctr_croibg"])
            ctr_croibg_errors.append(dsc["counters"]["ctr_croibg_errors"])
            croi.append(dsc["counters"]["croi"])
            bgroi.append(dsc["counters"]["bgroi"])
            croi_pix.append(dsc["counters"]["croi_pix"])
            bgroi_pix.append(dsc["counters"]["bgroi_pix"])
            if dsc["counters"]["croibg_bgimg"] is not None:
                croibg_bgimg.append(dsc["counters"]["croibg_bgimg"])
                croibg_bgimg_errors.append(dsc["counters"]["croibg_bgimg_errors"])
            Cfactors_croi.append(dsc["counters"]["Cfactors_croi"])
            Cfactors_bgroi.append(dsc["counters"]["Cfactors_bgroi"])
            Cfactor_croi.append(dsc["counters"]["Cfactor_croi"])
            Cfactor_bgroi.append(dsc["counters"]["Cfactor_bgroi"])
            Pfactor_croi.append(dsc["counters"]["Pfactor_croi"])
            Pfactor_bgroi.append(dsc["counters"]["Pfactor_bgroi"])
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
            logger.exception("Unexpected exception while creating data sets to save")
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
        "ctr_croibg": np.vstack(ctr_croibg),
        "ctr_croibg_errors": np.vstack(ctr_croibg_errors),
        "croi": np.vstack(croi),
        "bgroi": np.vstack(bgroi),
        "croi_pix": np.vstack(croi_pix),
        "bgroi_pix": np.vstack(bgroi_pix),
        "Cfactors_croi": np.vstack(Cfactors_croi),
        "Cfactors_bgroi": np.vstack(Cfactors_bgroi),
        "Cfactor_croi": np.vstack(Cfactor_croi),
        "Cfactor_bgroi": np.vstack(Cfactor_bgroi),
        "Pfactor_croi": np.vstack(Pfactor_croi),
        "Pfactor_bgroi": np.vstack(Pfactor_bgroi),
        "bgimg_croi": np.vstack(bgimg_croi),
        "bgimg_bgroi": np.vstack(bgimg_bgroi),
        "x": np.array(x),
        "y": np.array(y),
        "vsize": np.array(vsize),
        "hsize": np.array(hsize),
        "axis": np.vstack(axis),
        # True scattering angles are repeated for each extracted curve
        # because its peak can occur at a different source frame.
        **_rocking_arm_snapshot(gamma_arm, delta_arm, np.vstack(alpha).shape),
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
            rois[lbl] = (
                np.vstack(optional_labels[lbl]).reshape(-1)
                if lbl == "s"
                else np.vstack(optional_labels[lbl])
            )

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
    options = self.scanSelector.get_integration_options()
    beam_profile = sample_length = None
    if options["footprint"]:
        footprint_dialog = self.scanSelector.correctionsDialog.footprintOptions_shared()
        beam_profile = footprint_dialog.beamProfile()
        sample_length = footprint_dialog.sampleLength()
    frame_policy = _frame_policy_with_progress(
        self,
        self.fscan,
        config_snapshot.corrections,
        rois["axis"].shape[1],
        use_normalization=options["normalization"],
        use_illumination=options["footprint"],
        alpha=np.deg2rad(rois["alpha"]),
        beam_profile=beam_profile,
        sample_length=sample_length,
    )
    versioned_policy = frame_policy.new_contract or bool(
        frame_policy.interception_provenance
    )
    profile_provenance = _curve_profile_provenance(config_snapshot.corrections)
    profile_provenance.update(frame_policy.interception_provenance)
    profile_provenance.update(
        {
            "wavelength_angstrom": config_snapshot.ub_calculator.getLambda(),
            "unitcell_area_angstrom2": config_snapshot.unit_cell.uc_area,
            "detector_efficiency_assumed": 1.0,
            "external_transmission_assumed": 1.0,
        }
    )
    # The preserved ``rois`` group remains legacy-compatible. An explicit
    # primary-monitor/total-flux setup activates the sibling framewise
    # contract; its Q/H arrays are applied by the rocking reducer before
    # angular aggregation.
    curve_record = CurveCorrectionRecord(
        algorithm=(
            ("shape_interception_total_flux_v1" if frame_policy.interception_provenance
             else "framewise_ctr_total_flux_v1")
            if frame_policy.new_contract
            else ("shape_interception_legacy_v1" if versioned_policy
                  else "legacy_rocking_roi_v2")
        ),
        output_quantity=(
            "rocking_ctr_photon_curve"
            if frame_policy.new_contract
            else ("rocking_shape_density_curve" if versioned_policy
                  else "rocking_roi_curve")
        ),
        scale_convention=(
            frame_policy.scale_convention
            if versioned_policy
            else "legacy_unnormalized"
        ),
        normalization_status=(
            frame_policy.normalization_status
            if versioned_policy
            else "not_applied"
        ),
        illumination_status=(
            frame_policy.illumination_status
            if versioned_policy
            else "not_applied"
        ),
        pixel_correction_status=("applied" if corr else "not_applied"),
        normalization_divisor=(
            frame_policy.normalization_divisor if versioned_policy else None
        ),
        normalization_unit=(
            frame_policy.normalization_unit if versioned_policy else None
        ),
        normalization_components=(
            frame_policy.normalization_components if versioned_policy else ()
        ),
        illumination_divisor=(
            frame_policy.illumination_divisor if versioned_policy else None
        ),
        illumination_convention=(
            frame_policy.illumination_convention if versioned_policy else None
        ),
        vertical_intercepted_fraction=(
            frame_policy.vertical_intercepted_fraction
            if versioned_policy
            else None
        ),
        horizontal_intercepted_fraction=(
            frame_policy.horizontal_intercepted_fraction
            if versioned_policy
            else None
        ),
        intercepted_fraction=(
            frame_policy.intercepted_fraction if versioned_policy else None
        ),
        alpha=np.deg2rad(rois["alpha"]),
        base_croi=rois["croi"],
        base_croi_variance=rois["croi"],
        base_bgroi=rois["bgroi"],
        base_bgroi_variance=rois["bgroi"],
        base_croibg=rois["ctr_croibg"],
        base_croibg_variance=np.square(rois["ctr_croibg_errors"]),
        combined_croi_factor=rois["Cfactor_croi"],
        combined_bgroi_factor=rois["Cfactor_bgroi"],
        polarization_croi_factor=rois["Pfactor_croi"],
        polarization_bgroi_factor=rois["Pfactor_bgroi"],
        gamma_arm=rois["gamma_arm"],
        delta_arm=rois["delta_arm"],
        roi_x=rois["x"],
        roi_y=rois["y"],
        roi_width=rois["hsize"],
        roi_height=rois["vsize"],
        roi_x_start=rocking_roi_edges["x_start"],
        roi_x_stop=rocking_roi_edges["x_stop"],
        roi_y_start=rocking_roi_edges["y_start"],
        roi_y_stop=rocking_roi_edges["y_stop"],
        profile_provenance=profile_provenance,
    )
    data_2d_structured[self.activescanname]["measurement"][name][
        CURVE_CORRECTIONS_GROUP
    ] = curve_correction_record_to_nxdict(curve_record)

    return data_2d_structured[self.activescanname]["measurement"][name]


def fixed_roi_geometry(self, xy, **kwargs):
    """Calculate hkl, detector angles, and pixel metadata for fixed ROIs.

    :param xy:
        ROI center coordinates in detector pixels.
    :returns:
        Array of shape (ROI, frame, 6): h, k, l in r.l.u., delta and
        gamma in radians, and the source scan axis.
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

    count = len(self.fscan)
    mu_frames = np.broadcast_to(np.asarray(mu, dtype=np.float64), (count,))
    # A fixed pixel looks in a different direction on every frame once the
    # detector arm moves, so the conversion follows the arm frame by frame.
    arm_groups = self.armFrameGroups(count)

    hkl_del_gam = np.empty((xy.shape[0], count, 6), dtype=np.float64)
    for i, xy_i in enumerate(xy):
        x = np.full(count, xy_i[0])
        y = np.full(count, xy_i[1])
        gamma = np.empty(count, dtype=np.float64)
        delta = np.empty(count, dtype=np.float64)
        alpha = np.empty(count, dtype=np.float64)
        for selection, gamma_arm, delta_arm in arm_groups:
            gamma_s, delta_s, alpha_s = self.ubcalc.detectorCal.crystalAnglesPoint(
                np.atleast_1d(y[selection]),
                np.atleast_1d(x[selection]),
                mu_frames[selection],
                self.ubcalc.n,
                gamma_arm,
                delta_arm,
            )
            gamma[selection] = gamma_s
            delta[selection] = delta_s
            alpha[selection] = alpha_s

        hkl = self.ubcalc.angles.anglesToHkl(
            alpha, delta, gamma, om, self.ubcalc.chi, self.ubcalc.phi
        )
        # for i in range(len(self.fscan)):

        hkl_del_gam[i, :, :3] = np.array(hkl).T
        hkl_del_gam[i, :, 3] = delta
        hkl_del_gam[i, :, 4] = gamma
        hkl_del_gam[i, :, 5] = self.fscan.axis
    return hkl_del_gam

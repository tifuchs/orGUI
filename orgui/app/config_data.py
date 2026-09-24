# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#
# ###########################################################################*/
"""Structured application configuration snapshots for orGUI databases."""

from __future__ import annotations

import configparser
import datetime
import json
from dataclasses import dataclass, field
from typing import Any

import h5py
import numpy as np
from silx.io.dictdump import dicttonx, nxtodict

from ..datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg

SCHEMA_VERSION = 1

#: Layout version of the integration_corrections group. Version 2 introduced
#: the typed layout; version 3 adds explicit total-flux and primary-monitor
#: settings without reinterpreting the legacy flux-density fields.
CORRECTIONS_SCHEMA_VERSION = 3

#: Layout version and group name of an applied per-curve correction record.
#: This deliberately is not named ``rois``: an older reducer must not mistake
#: a future frame-normalized curve for the legacy unnormalized input.
CURVE_CORRECTIONS_SCHEMA_VERSION = 3
CURVE_CORRECTIONS_GROUP = "ctr_curve_v3"

CORRECTION_STATUSES = frozenset(
    {"applied", "not_applied", "unavailable", "unknown"}
)

#: Layout version of the roi_integration group.
ROI_SCHEMA_VERSION = 1


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        return value.decode()
    return value


@dataclass
class CorrectionState:
    """Persist the correction selections shared by orGUI workflows.

    Large mask, background, and variance arrays are stored in a job asset
    bundle. The corresponding fields here identify datasets in that bundle.
    """

    use_mask: bool = False
    use_background: bool = False
    use_solid_angle: bool = False
    use_polarization: bool = False
    # Tri-state on purpose: ``None`` is "this configuration predates the
    # switch being stored", which must leave the GUI as it is rather than
    # silently turning the correction off.
    use_lorentz: bool | None = None
    use_footprint: bool | None = None
    use_normalization: bool | None = None
    repair_masked_pixels: bool = False
    repair_max_component_pixels: int | None = None
    repair_max_span: int | None = None
    repair_radius: int | None = None
    repair_min_valid_neighbors: int | None = None
    repair_use_pyfai_gaps: bool = True
    repair_gap_size_px: int = 1
    normalize_exposure: bool = True
    monitor_corrections: tuple[str, ...] = ()
    # Older prepared reconstruction jobs ignored use_normalization and the
    # primary-monitor contract. Missing means their original legacy behavior.
    shared_frame_normalization: bool = False
    excluded_frames: tuple[int, ...] = ()
    mask_asset: str | None = None
    background_asset: str | None = None
    background_variance_asset: str | None = None
    uncertainty_provenance: dict[str, Any] = field(default_factory=dict)
    # Sample size, in meter, along and across the beam: inputs to the beam
    # footprint dialog, not the footprint itself. ``None`` means "this
    # configuration predates them, or the dialog was never opened", which
    # must leave the dialog at its own defaults rather than zeroing it out.
    sample_length_m: float | None = None
    sample_width_m: float | None = None
    # Incident flux density Phi_0, in photons per second per square meter
    # (see corrections.measurement.scale_factor). Optional: only needed for
    # an absolutely scaled structure factor (issue #15); a relative one does
    # not use it.
    beam_flux_density: float | None = None
    # New total-flux convention. These fields describe requested settings;
    # the exact values actually used for a curve live in CurveCorrectionRecord.
    # ``None`` means unknown/not configured, never physical zero.
    total_incident_flux: float | None = None
    total_flux_calibrated: bool | None = None
    primary_monitor: str | None = None
    primary_monitor_kind: str | None = None
    primary_monitor_unit: str | None = None
    monitor_reference_reading: float | None = None
    monitor_reference_exposure_s: float | None = None
    horizontal_interception: str | None = None
    horizontal_intercepted_fraction: float | None = None
    # The beam shape the footprint dialog describes the incident beam with --
    # analytical or measured, and every input either needs. Kept in the same
    # units the dialog itself shows them in (IntegrationCorrectionsDialog.
    # settings()), not converted to SI: an analytical shape's numeric values
    # are only meaningful together with its name (a width in micrometer for
    # one shape, a dimensionless flatness or skew for another), so there is
    # no single physical unit to convert them to. ``None``/empty means "this
    # configuration predates them, or the dialog was never opened".
    beam_shape_analytical: bool | None = None
    beam_shape_name: str | None = None
    beam_shape_values: tuple[float, ...] = ()
    beam_profile_file: str | None = None
    beam_profile_content: str | None = None
    beam_profile_unit: str | None = None
    beam_profile_center: str | None = None
    beam_profile_offset_um: float | None = None
    # Embedded, normalized measured profile. A path alone is not sufficient
    # provenance once the source file moves or changes.
    beam_profile_positions_m: tuple[float, ...] = ()
    beam_profile_density_per_m: tuple[float, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-compatible correction-state dictionary."""
        result = _json_value(self.__dict__)
        result["monitor_corrections"] = list(self.monitor_corrections)
        result["excluded_frames"] = list(self.excluded_frames)
        return result

    @classmethod
    def from_dict(cls, values):
        """Build correction state from a JSON-compatible dictionary."""
        values = dict(values or {})
        values["monitor_corrections"] = tuple(
            values.get("monitor_corrections", ())
        )
        values["excluded_frames"] = tuple(
            int(value) for value in values.get("excluded_frames", ())
        )
        values["beam_shape_values"] = tuple(
            float(value) for value in values.get("beam_shape_values", ())
        )
        values["beam_profile_positions_m"] = tuple(
            float(value) for value in values.get("beam_profile_positions_m", ())
        )
        values["beam_profile_density_per_m"] = tuple(
            float(value)
            for value in values.get("beam_profile_density_per_m", ())
        )
        return cls(**values)


_TOTAL_FLUX_STATE_FIELDS = (
    "total_incident_flux",
    "total_flux_calibrated",
    "primary_monitor",
    "primary_monitor_kind",
    "primary_monitor_unit",
    "monitor_reference_reading",
    "monitor_reference_exposure_s",
    "horizontal_interception",
    "horizontal_intercepted_fraction",
)


@dataclass
class CurveCorrectionRecord:
    """Applied correction provenance for one saved curve branch.

    Arrays are stored exactly as used, with radians at geometry boundaries.
    The record is separate from :class:`CorrectionState`: that class captures
    requested settings, while this one records what happened to the saved
    numbers. Optional fields remain absent rather than being guessed as unity.

    ``base_croibg`` and its variance are the reversible scalar curve before
    normalization and illumination divisors. They are not raw detector images
    and cannot undo arbitrary pixel reweighting.
    """

    algorithm: str
    output_quantity: str
    scale_convention: str
    normalization_status: str = "unknown"
    illumination_status: str = "unknown"
    pixel_correction_status: str = "unknown"
    normalization_divisor: Any = None
    normalization_unit: str | None = None
    normalization_components: tuple[str, ...] = ()
    illumination_divisor: Any = None
    illumination_convention: str | None = None
    vertical_intercepted_fraction: Any = None
    horizontal_intercepted_fraction: Any = None
    intercepted_fraction: Any = None
    alpha: Any = None
    base_croi: Any = None
    base_croi_variance: Any = None
    base_bgroi: Any = None
    base_bgroi_variance: Any = None
    base_croibg: Any = None
    base_croibg_variance: Any = None
    combined_croi_factor: Any = None
    combined_bgroi_factor: Any = None
    polarization_croi_factor: Any = None
    polarization_bgroi_factor: Any = None
    gamma_arm: Any = None
    delta_arm: Any = None
    roi_x: Any = None
    roi_y: Any = None
    roi_width: Any = None
    roi_height: Any = None
    roi_x_start: Any = None
    roi_x_stop: Any = None
    roi_y_start: Any = None
    roi_y_stop: Any = None
    detector_acceptance: Any = None
    lorentz_mode: str | None = None
    profile_provenance: dict[str, Any] = field(default_factory=dict)


def _as_text(value):
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def _none_to_nan(value):
    return np.nan if value is None else value


def _string_array(strings):
    encoded = [str(value).encode("utf-8") for value in strings]
    if not encoded:
        return np.empty((0, 0), dtype=np.uint8)
    width = max(len(value) for value in encoded)
    data = np.zeros((len(encoded), width), dtype=np.uint8)
    for index, value in enumerate(encoded):
        data[index, : len(value)] = np.frombuffer(value, dtype=np.uint8)
    return data


def _read_string_array(data):
    array = np.asarray(data)
    if array.dtype.kind in "iu" and array.ndim == 2:
        array = array.astype(np.uint8, copy=False)
        return [bytes(row[row != 0]).decode("utf-8") for row in array]
    return [_as_text(value) for value in data]


def detector_to_nxdict(detector):
    """Serialize a :class:`Detector2D_SXRD` using its central NX serializer."""
    return detector.toNXdict()["detector_SXRD"]


def detector_from_nxdict(nxdict):
    """Load a :class:`Detector2D_SXRD` through the central NX reader."""
    return DetectorCalibration.loadNXdict(nxdict)


def unit_cell_to_nxdict(unit_cell):
    """Serialize a :class:`UnitCell` without fit parameters.

    Lattice lengths are Angstrom and lattice angles are degrees, matching
    :class:`UnitCell` construction.
    """
    return {
        "@NX_class": "NXsample",
        "unit_cell_abc": np.asarray(unit_cell.a, dtype=np.float64),
        "unit_cell_alphabetagamma": np.rad2deg(
            np.asarray(unit_cell.alpha, dtype=np.float64)
        ),
        "unit_cell_volume": float(unit_cell.volume),
        "orgui_unit_cell": {
            "@NX_class": "NXcollection",
            "class_name": "UnitCell",
            "name": unit_cell.name,
            "atom_names": _string_array(unit_cell.names),
            "basis": np.asarray(unit_cell.basis, dtype=np.float64),
            "fractional_coordinates": np.asarray(
                unit_cell.basis[:, 1:4], dtype=np.float64
            )
            if unit_cell.basis.size
            else np.empty((0, 3), dtype=np.float64),
            "debye_waller": np.asarray(unit_cell.basis[:, 4:6], dtype=np.float64)
            if unit_cell.basis.size
            else np.empty((0, 2), dtype=np.float64),
            "occupancies": np.asarray(unit_cell.basis[:, 6], dtype=np.float64)
            if unit_cell.basis.size
            else np.empty((0,), dtype=np.float64),
            "layers": np.asarray(unit_cell.basis[:, 7], dtype=np.float64)
            if unit_cell.basis.size and unit_cell.basis.shape[1] > 7
            else np.empty((0,), dtype=np.float64),
            "coherent_domain_matrices": np.asarray(
                unit_cell.coherentDomainMatrix, dtype=np.float64
            ),
            "coherent_domain_occupancies": np.asarray(
                unit_cell.coherentDomainOccupancy, dtype=np.float64
            ),
        },
    }


def unit_cell_from_nxdict(nxdict):
    """Load a :class:`UnitCell` from a config sample dictionary."""
    sample = nxdict
    if "sample" in sample:
        sample = sample["sample"]
    ucdata = sample.get("orgui_unit_cell", {})
    a = np.asarray(sample["unit_cell_abc"], dtype=np.float64)
    alpha = np.asarray(sample["unit_cell_alphabetagamma"], dtype=np.float64)
    name = _as_text(ucdata.get("name", "unnamed"))
    unit_cell = CTRcalc.UnitCell(a, alpha, name=name)
    names = _read_string_array(ucdata.get("atom_names", []))
    basis = np.asarray(ucdata.get("basis", np.empty((0, 8))), dtype=np.float64)
    for atom_name, row in zip(names, basis):
        unit_cell.addAtom(
            atom_name,
            row[1:4],
            _none_to_nan(row[4]),
            _none_to_nan(row[5]),
            _none_to_nan(row[6]),
            row[7] if row.size > 7 else 0,
        )
    if "coherent_domain_matrices" in ucdata:
        unit_cell.coherentDomainMatrix = [
            np.asarray(matrix, dtype=np.float64)
            for matrix in ucdata["coherent_domain_matrices"]
        ]
    if "coherent_domain_occupancies" in ucdata:
        unit_cell.coherentDomainOccupancy = [
            float(occupancy) for occupancy in ucdata["coherent_domain_occupancies"]
        ]
    return unit_cell


def reflections_to_nxdict(reflections, ub_calculator=None):
    """Serialize reference reflections.

    Reflection HKL values are in r.l.u., detector coordinates are pixels, and
    six-circle angle snapshots are stored in degrees as audit metadata.
    """
    hkls = []
    xy = []
    imageno = []
    identifiers = []
    angles = []
    for idx, reflection in enumerate(reflections or []):
        hkls.append(np.asarray(reflection.hkl, dtype=np.float64))
        xy.append(np.asarray(reflection.xy, dtype=np.float64))
        imageno.append(int(reflection.imageno))
        identifiers.append(getattr(reflection, "identifier", f"ref_{idx}"))
        angle = np.full((6,), np.nan)
        angles.append(np.rad2deg(np.asarray(angle, dtype=np.float64)))
    return {
        "@NX_class": "NXdata",
        "hkl": np.asarray(hkls, dtype=np.float64).reshape((-1, 3)),
        "xy": np.asarray(xy, dtype=np.float64).reshape((-1, 2)),
        "image_number": np.asarray(imageno, dtype=np.int64),
        "identifier": _string_array(identifiers),
        "sixc_angles": {
            "@NX_class": "NXpositioner",
            "@unit": "deg",
            "angles": np.asarray(angles, dtype=np.float64).reshape((-1, 6)),
            "axis_names": _string_array(
                ["alpha", "delta", "gamma", "omega", "chi", "phi"]
            ),
        },
    }


def reflections_from_nxdict(nxdict):
    """Load reference reflections from a config dictionary."""
    from .QReflectionSelector import HKLReflection

    refl = nxdict.get("reference_reflections", nxdict)
    hkls = np.asarray(refl.get("hkl", np.empty((0, 3))), dtype=np.float64).reshape(
        (-1, 3)
    )
    xy = np.asarray(refl.get("xy", np.empty((0, 2))), dtype=np.float64).reshape((-1, 2))
    imageno = np.asarray(refl.get("image_number", []), dtype=np.int64)
    identifiers = _read_string_array(refl.get("identifier", []))
    return [
        HKLReflection(point, hkl, int(img), ident)
        for point, hkl, img, ident in zip(xy, hkls, imageno, identifiers)
    ]


#: Region-of-interest option group -> the unit attributes written beside its
#: values. Anything absent here is dimensionless.
_ROI_UNITS = {
    "region": {"@unit": "px"},
    "advanced": {"@sample_size_unit": "m", "@offset_unit": "px"},
    "rocking_scan": {"@unit": "rlu"},
}


def _nx_group(mapping, units=None):
    """A NeXus subgroup from a flat mapping, skipping ``None`` values."""
    group = {"@NX_class": "NXcollection"}
    group.update(units or {})
    for key, value in mapping.items():
        if value is None:
            continue
        group[key] = value
    return group


def _plain(value):
    """Strip the numpy and bytes wrappers that HDF5 hands back.

    Live datasets with two or more dimensions (per-curve, per-frame arrays)
    are returned unread so callers can slice one curve lazily.
    """
    if isinstance(value, h5py.Dataset):
        return value if value.ndim >= 2 else _plain(value[()])
    if isinstance(value, bytes):
        return value.decode()
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return _plain(value[()])
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read_group(nxdict, name):
    """Values of one subgroup, without its NeXus bookkeeping keys."""
    group = nxdict.get(name) or {}
    return {
        key: _plain(value)
        for key, value in group.items()
        if not key.startswith("@")
    }


def corrections_to_nxdict(state):
    """Serialize a :class:`CorrectionState` as a typed NeXus group.

    Replaces the single opaque JSON string this used to be written as. Every
    value is its own dataset, so a stored configuration can be read in an
    HDF5 browser and units can sit beside the numbers that have them.

    ``None`` and empty sequences are written by *omission*: HDF5 has no null,
    and absence already means "not recorded" to the reader.

    :param CorrectionState state: The state to serialize.
    :rtype: dict
    :raises ValueError: If ``uncertainty_provenance`` is not flat.
    """
    if any(
        isinstance(value, (dict, list, tuple))
        for value in state.uncertainty_provenance.values()
    ):
        raise ValueError(
            "uncertainty_provenance must be a flat mapping of scalars; a "
            "nested value cannot be written as a NeXus group"
        )
    if len(state.beam_profile_positions_m) != len(
        state.beam_profile_density_per_m
    ):
        raise ValueError(
            "embedded beam profile positions and density must have equal length"
        )
    nxdict = {
        "@NX_class": "NXcollection",
        "@orgui_schema_version": CORRECTIONS_SCHEMA_VERSION,
        "switches": _nx_group({
            "use_mask": state.use_mask,
            "use_background": state.use_background,
            "use_solid_angle": state.use_solid_angle,
            "use_polarization": state.use_polarization,
            "use_lorentz": state.use_lorentz,
            "use_footprint": state.use_footprint,
            "use_normalization": state.use_normalization,
        }),
        "normalization": _nx_group(
            {
                "normalize_exposure": state.normalize_exposure,
                "shared_frame_normalization": state.shared_frame_normalization,
                "primary_monitor": state.primary_monitor,
                "primary_monitor_kind": state.primary_monitor_kind,
                "primary_monitor_unit": state.primary_monitor_unit,
                "monitor_reference_reading": state.monitor_reference_reading,
                "monitor_reference_exposure": (
                    state.monitor_reference_exposure_s
                ),
            },
            units={"@monitor_reference_exposure_unit": "s"},
        ),
        "pixel_repair": _nx_group({
            "enabled": state.repair_masked_pixels,
            "max_component_pixels": state.repair_max_component_pixels,
            "max_span": state.repair_max_span,
            "radius": state.repair_radius,
            "min_valid_neighbors": state.repair_min_valid_neighbors,
            "use_pyfai_gaps": state.repair_use_pyfai_gaps,
            "gap_size_px": state.repair_gap_size_px,
        }),
        "assets": _nx_group({
            "mask": state.mask_asset,
            "background": state.background_asset,
            "background_variance": state.background_variance_asset,
        }),
        "footprint": _nx_group(
            {
                "sample_length": state.sample_length_m,
                "sample_width": state.sample_width_m,
                "beam_flux_density": state.beam_flux_density,
                "total_incident_flux": state.total_incident_flux,
                "total_flux_calibrated": state.total_flux_calibrated,
                "horizontal_interception": state.horizontal_interception,
                "horizontal_intercepted_fraction": (
                    state.horizontal_intercepted_fraction
                ),
            },
            units={
                "@sample_size_unit": "m",
                "@beam_flux_density_unit": "1/(s m^2)",
                "@total_incident_flux_unit": "photons/s",
            },
        ),
        # The dialog's own display units, not SI: see CorrectionState.
        "beam_shape": _nx_group(
            {
                "analytical": state.beam_shape_analytical,
                "shape": state.beam_shape_name,
                "profile_file": state.beam_profile_file,
                "profile_content": state.beam_profile_content,
                "profile_unit": state.beam_profile_unit,
                "profile_center": state.beam_profile_center,
                "profile_offset": state.beam_profile_offset_um,
            },
            units={
                "@profile_offset_unit": "micron",
                "@profile_position_unit": "m",
                "@profile_density_unit": "1/m",
            },
        ),
    }
    if state.monitor_corrections:
        nxdict["normalization"]["monitor_corrections"] = _string_array(
            list(state.monitor_corrections)
        )
    if state.excluded_frames:
        nxdict["excluded_frames"] = np.asarray(
            state.excluded_frames, dtype=np.int64
        )
    if state.beam_shape_values:
        nxdict["beam_shape"]["shape_values"] = np.asarray(
            state.beam_shape_values, dtype=np.float64
        )
    if state.beam_profile_positions_m:
        nxdict["beam_shape"]["profile_positions"] = np.asarray(
            state.beam_profile_positions_m, dtype=np.float64
        )
    if state.beam_profile_density_per_m:
        nxdict["beam_shape"]["profile_density"] = np.asarray(
            state.beam_profile_density_per_m, dtype=np.float64
        )
    if state.uncertainty_provenance:
        nxdict["uncertainty_provenance"] = _nx_group(
            dict(state.uncertainty_provenance)
        )
    return nxdict


def corrections_from_nxdict(nxdict):
    """Rebuild a :class:`CorrectionState` from its NeXus group.

    Datasets this version does not know are ignored, so a configuration
    written by a newer orGUI still loads.

    :param dict nxdict: The ``integration_corrections`` group.
    :rtype: CorrectionState
    """
    switches = _read_group(nxdict, "switches")
    normalization = _read_group(nxdict, "normalization")
    repair = _read_group(nxdict, "pixel_repair")
    assets = _read_group(nxdict, "assets")
    excluded = _plain(nxdict.get("excluded_frames"))
    values = {
        "use_mask": bool(switches.get("use_mask", False)),
        "use_background": bool(switches.get("use_background", False)),
        "use_solid_angle": bool(switches.get("use_solid_angle", False)),
        "use_polarization": bool(switches.get("use_polarization", False)),
        "repair_masked_pixels": bool(repair.get("enabled", False)),
        "repair_use_pyfai_gaps": bool(repair.get("use_pyfai_gaps", True)),
        "repair_gap_size_px": int(repair.get("gap_size_px", 1)),
        "normalize_exposure": bool(
            normalization.get("normalize_exposure", True)
        ),
        "shared_frame_normalization": bool(
            normalization.get("shared_frame_normalization", False)
        ),
        # Written with _string_array, i.e. a uint8 matrix, so it needs the
        # matching reader rather than a plain tuple().
        "monitor_corrections": tuple(
            _read_string_array(normalization.get("monitor_corrections", []))
        ),
        "excluded_frames": tuple(int(v) for v in (excluded or ())),
        "uncertainty_provenance": _read_group(nxdict, "uncertainty_provenance"),
    }
    for name in (
        "primary_monitor",
        "primary_monitor_kind",
        "primary_monitor_unit",
    ):
        if name in normalization:
            values[name] = str(normalization[name])
    if "monitor_reference_reading" in normalization:
        values["monitor_reference_reading"] = float(
            normalization["monitor_reference_reading"]
        )
    if "monitor_reference_exposure" in normalization:
        values["monitor_reference_exposure_s"] = float(
            normalization["monitor_reference_exposure"]
        )
    for name in ("use_lorentz", "use_footprint", "use_normalization"):
        if name in switches:
            values[name] = bool(switches[name])
    for name, key in (
        ("repair_max_component_pixels", "max_component_pixels"),
        ("repair_max_span", "max_span"),
        ("repair_radius", "radius"),
        ("repair_min_valid_neighbors", "min_valid_neighbors"),
    ):
        if key in repair:
            values[name] = int(repair[key])
    for name, key in (
        ("mask_asset", "mask"),
        ("background_asset", "background"),
        ("background_variance_asset", "background_variance"),
    ):
        if key in assets:
            values[name] = str(assets[key])
    footprint = _read_group(nxdict, "footprint")
    for name, key in (
        ("sample_length_m", "sample_length"),
        ("sample_width_m", "sample_width"),
        ("beam_flux_density", "beam_flux_density"),
        ("total_incident_flux", "total_incident_flux"),
        ("horizontal_intercepted_fraction", "horizontal_intercepted_fraction"),
    ):
        if key in footprint:
            values[name] = float(footprint[key])
    if "total_flux_calibrated" in footprint:
        values["total_flux_calibrated"] = bool(
            footprint["total_flux_calibrated"]
        )
    if "horizontal_interception" in footprint:
        values["horizontal_interception"] = str(
            footprint["horizontal_interception"]
        )
    beam_shape = _read_group(nxdict, "beam_shape")
    if "analytical" in beam_shape:
        values["beam_shape_analytical"] = bool(beam_shape["analytical"])
    for name, key in (
        ("beam_shape_name", "shape"),
        ("beam_profile_file", "profile_file"),
        ("beam_profile_content", "profile_content"),
        ("beam_profile_unit", "profile_unit"),
        ("beam_profile_center", "profile_center"),
    ):
        if key in beam_shape:
            values[name] = str(beam_shape[key])
    if "profile_offset" in beam_shape:
        values["beam_profile_offset_um"] = float(beam_shape["profile_offset"])
    if "shape_values" in beam_shape:
        values["beam_shape_values"] = tuple(
            float(v) for v in np.atleast_1d(beam_shape["shape_values"])
        )
    if "profile_positions" in beam_shape:
        values["beam_profile_positions_m"] = tuple(
            float(v) for v in np.atleast_1d(beam_shape["profile_positions"])
        )
    if "profile_density" in beam_shape:
        values["beam_profile_density_per_m"] = tuple(
            float(v) for v in np.atleast_1d(beam_shape["profile_density"])
        )
    return CorrectionState(**values)


def _group_metadata(group, name, default=None):
    """Read one ``@`` metadata value from a dict or live HDF5 group."""
    if hasattr(group, "attrs"):
        return _plain(group.attrs.get(name.removeprefix("@"), default))
    return _plain(group.get(name, default))


def _validate_correction_status(name, value):
    if value not in CORRECTION_STATUSES:
        raise ValueError(
            f"{name} must be one of {sorted(CORRECTION_STATUSES)}, got {value!r}"
        )


def curve_correction_record_to_nxdict(record):
    """Serialize an applied :class:`CurveCorrectionRecord` to NeXus.

    The returned group is intended to live under
    :data:`CURVE_CORRECTIONS_GROUP`, beside rather than inside legacy
    ``rois``/``counters`` data. This makes its different normalization
    contract visible to both humans and dispatch code.

    :param CurveCorrectionRecord record: Applied curve provenance.
    :rtype: dict
    :raises ValueError: If a status or profile value is not representable.
    """
    for name in (
        "normalization_status",
        "illumination_status",
        "pixel_correction_status",
    ):
        _validate_correction_status(name, getattr(record, name))
    if any(
        isinstance(value, (dict, list, tuple))
        for value in record.profile_provenance.values()
    ):
        raise ValueError("profile_provenance must be a flat mapping of scalars")

    result = {
        "@NX_class": "NXcollection",
        "@orgui_schema_version": CURVE_CORRECTIONS_SCHEMA_VERSION,
        "@orgui_curve_contract": "frame_corrections",
        "identity": _nx_group({
            "algorithm": record.algorithm,
            "output_quantity": record.output_quantity,
            "scale_convention": record.scale_convention,
        }),
        "normalization": _nx_group({
            "status": record.normalization_status,
            "divisor": record.normalization_divisor,
            "divisor_unit": record.normalization_unit,
        }),
        "illumination": _nx_group(
            {
                "status": record.illumination_status,
                "divisor": record.illumination_divisor,
                "convention": record.illumination_convention,
                "vertical_intercepted_fraction": (
                    record.vertical_intercepted_fraction
                ),
                "horizontal_intercepted_fraction": (
                    record.horizontal_intercepted_fraction
                ),
                "intercepted_fraction": record.intercepted_fraction,
                "alpha": record.alpha,
            },
            units={"@alpha_unit": "rad"},
        ),
        "base": _nx_group({
            "croi": record.base_croi,
            "croi_variance": record.base_croi_variance,
            "bgroi": record.base_bgroi,
            "bgroi_variance": record.base_bgroi_variance,
            "croibg": record.base_croibg,
            "croibg_variance": record.base_croibg_variance,
        }),
        "pixel_corrections": _nx_group({
            "status": record.pixel_correction_status,
            "combined_croi": record.combined_croi_factor,
            "combined_bgroi": record.combined_bgroi_factor,
            "polarization_croi": record.polarization_croi_factor,
            "polarization_bgroi": record.polarization_bgroi_factor,
        }),
        "geometry": _nx_group(
            {
                "gamma_arm": record.gamma_arm,
                "delta_arm": record.delta_arm,
                "roi_x": record.roi_x,
                "roi_y": record.roi_y,
                "roi_width": record.roi_width,
                "roi_height": record.roi_height,
                "roi_x_start": record.roi_x_start,
                "roi_x_stop": record.roi_x_stop,
                "roi_y_start": record.roi_y_start,
                "roi_y_stop": record.roi_y_stop,
                "detector_acceptance": record.detector_acceptance,
                "lorentz_mode": record.lorentz_mode,
            },
            units={
                "@detector_arm_unit": "rad",
                "@roi_position_unit": "px",
                "@detector_acceptance_unit": "rad",
            },
        ),
        "profile": _nx_group(
            dict(record.profile_provenance),
            units={
                "@profile_position_unit": "m",
                "@profile_density_unit": "1/m",
                "@profile_offset_unit": "micron",
            },
        ),
    }
    if record.normalization_components:
        result["normalization"]["components"] = _string_array(
            record.normalization_components
        )
    return result


def curve_correction_record_from_nxdict(nxdict):
    """Load a :class:`CurveCorrectionRecord` from its versioned group.

    :param dict nxdict: Contents of :data:`CURVE_CORRECTIONS_GROUP`.
    :rtype: CurveCorrectionRecord
    :raises ValueError: If this is not a supported versioned curve record.
    """
    version = int(_group_metadata(nxdict, "@orgui_schema_version", 0))
    contract = _as_text(_group_metadata(nxdict, "@orgui_curve_contract", ""))
    if version != CURVE_CORRECTIONS_SCHEMA_VERSION or contract != "frame_corrections":
        raise ValueError(
            f"unsupported curve correction record version/contract: "
            f"{version}/{contract!r}"
        )
    identity = _read_group(nxdict, "identity")
    normalization = _read_group(nxdict, "normalization")
    illumination = _read_group(nxdict, "illumination")
    base = _read_group(nxdict, "base")
    pixel = _read_group(nxdict, "pixel_corrections")
    geometry = _read_group(nxdict, "geometry")
    components = normalization.get("components", ())
    values = {
        "algorithm": str(identity["algorithm"]),
        "output_quantity": str(identity["output_quantity"]),
        "scale_convention": str(identity["scale_convention"]),
        "normalization_status": str(normalization.get("status", "unknown")),
        "illumination_status": str(illumination.get("status", "unknown")),
        "pixel_correction_status": str(pixel.get("status", "unknown")),
        "normalization_components": tuple(_read_string_array(components)),
        "profile_provenance": _read_group(nxdict, "profile"),
    }
    if "divisor_unit" in normalization:
        values["normalization_unit"] = str(normalization["divisor_unit"])
    if "convention" in illumination:
        values["illumination_convention"] = str(illumination["convention"])
    groups = (
        (normalization, {"normalization_divisor": "divisor"}),
        (
            illumination,
            {
                "illumination_divisor": "divisor",
                "vertical_intercepted_fraction": "vertical_intercepted_fraction",
                "horizontal_intercepted_fraction": "horizontal_intercepted_fraction",
                "intercepted_fraction": "intercepted_fraction",
                "alpha": "alpha",
            },
        ),
        (
            base,
            {
                "base_croi": "croi",
                "base_croi_variance": "croi_variance",
                "base_bgroi": "bgroi",
                "base_bgroi_variance": "bgroi_variance",
                "base_croibg": "croibg",
                "base_croibg_variance": "croibg_variance",
            },
        ),
        (
            pixel,
            {
                "combined_croi_factor": "combined_croi",
                "combined_bgroi_factor": "combined_bgroi",
                "polarization_croi_factor": "polarization_croi",
                "polarization_bgroi_factor": "polarization_bgroi",
            },
        ),
        (
            geometry,
            {
                "gamma_arm": "gamma_arm",
                "delta_arm": "delta_arm",
                "roi_x": "roi_x",
                "roi_y": "roi_y",
                "roi_width": "roi_width",
                "roi_height": "roi_height",
                "roi_x_start": "roi_x_start",
                "roi_x_stop": "roi_x_stop",
                "roi_y_start": "roi_y_start",
                "roi_y_stop": "roi_y_stop",
                "detector_acceptance": "detector_acceptance",
                "lorentz_mode": "lorentz_mode",
            },
        ),
    )
    for group, mapping in groups:
        for field_name, dataset_name in mapping.items():
            if dataset_name in group:
                values[field_name] = group[dataset_name]
    record = CurveCorrectionRecord(**values)
    for name in (
        "normalization_status",
        "illumination_status",
        "pixel_correction_status",
    ):
        _validate_correction_status(name, getattr(record, name))
    return record


def curve_record_kind(container):
    """Classify a loaded curve container without guessing missing provenance.

    :param container: Mapping or HDF5 group containing a saved curve.
    :returns: ``"versioned"``, ``"legacy"``, or ``"unknown"``.
    :rtype: str
    """
    if CURVE_CORRECTIONS_GROUP in container:
        group = container[CURVE_CORRECTIONS_GROUP]
        version = int(_group_metadata(group, "@orgui_schema_version", 0))
        contract = _as_text(
            _group_metadata(group, "@orgui_curve_contract", "")
        )
        if (
            version == CURVE_CORRECTIONS_SCHEMA_VERSION
            and contract == "frame_corrections"
        ):
            return "versioned"
        return "unknown"
    if "rois" in container or "croibg" in container:
        return "legacy"
    counters = container.get("counters") if hasattr(container, "get") else None
    if counters is not None and "croibg" in counters:
        return "legacy"
    return "unknown"


def roi_to_nxdict(state):
    """Serialize a :class:`ROIState` as a typed NeXus group.

    :param ROIState state: The settings to serialize.
    :rtype: dict
    """
    nxdict = {
        "@NX_class": "NXcollection",
        "@orgui_schema_version": ROI_SCHEMA_VERSION,
    }
    for name in ("region", "advanced", "rocking_scan"):
        values = getattr(state, name)
        if values:
            nxdict[name] = _nx_group(dict(values), _ROI_UNITS[name])
    return nxdict


def roi_from_nxdict(nxdict):
    """Rebuild a :class:`ROIState` from its NeXus group.

    :param dict nxdict: The ``roi_integration`` group, or ``None`` for a
        configuration written before these settings were stored.
    :rtype: ROIState
    """
    nxdict = nxdict or {}
    return ROIState(
        region=_read_group(nxdict, "region"),
        advanced=_read_group(nxdict, "advanced"),
        rocking_scan=_read_group(nxdict, "rocking_scan"),
    )


@dataclass
class ROIState:
    """Region-of-interest settings that decide what a scan integrates.

    Held as the same dictionaries :class:`~orgui.app.QScanSelector.QScanSelector`
    already speaks, so a new control in the options dialog reaches the file
    without a change here. The NeXus layout below is typed and carries units;
    this is only the carrier.

    An empty dictionary means "not recorded", which is what every
    configuration written before these settings were stored looks like.
    """

    #: Nominal sizes and the automatic-sizing switches, pixels.
    region: dict = field(default_factory=dict)
    #: The advanced options dialog: sample size in meter, offsets in pixels.
    advanced: dict = field(default_factory=dict)
    #: Rocking-scan ``s`` sampling, r.l.u. Saved with a *rocking* scan,
    #: ``delta_s`` is the effective value after the resolution clipping of
    #: ``onRoSChanged`` rather than what was typed. Saved with a stationary
    #: one it is whatever the control held, since nothing there clips it and
    #: the value is unused -- so do not read it as an effective sampling
    #: without checking which mode wrote it.
    rocking_scan: dict = field(default_factory=dict)

    def is_empty(self):
        """True when nothing was recorded, so nothing should be restored."""
        return not (self.region or self.advanced or self.rocking_scan)


@dataclass
class ConfigData:
    """Physical application state persisted with scans and integrations."""

    detector: DetectorCalibration.Detector2D_SXRD
    unit_cell: CTRcalc.UnitCell
    ub_calculator: HKLVlieg.UBCalculator
    mu: float = 0.0
    chi: float = 0.0
    phi: float = 0.0
    # Fixed detector arm position, radians, for backends that do not read the
    # arm motors. A scan supplying its own per-frame values overrides these.
    gamma_arm: float = 0.0
    delta_arm: float = 0.0
    # Which angles the arm motors report: 'prim' for the true scattering
    # angles gamma_p/delta_p, 'surface' for six-circle gamma/delta.
    arm_angle_frame: str = "prim"

    refraction_index: float = 1.0
    reference_reflections: list = field(default_factory=list)
    corrections: CorrectionState = field(default_factory=CorrectionState)
    roi: ROIState = field(default_factory=ROIState)
    orgui: dict = field(default_factory=dict)

    @classmethod
    def from_ini(cls, filename):
        """Read legacy INI config into a structured config object.

        Config ``SDD`` and ``pixelsize`` are read in meters, detector center is
        read in pixels, energy is keV, and lattice constants are Angstrom.
        """
        config = configparser.ConfigParser()
        config.read(filename)
        machine = config["Machine"]
        lattice = config["Lattice"]
        diffrac = config["Diffractometer"]

        energy = machine.getfloat("E", 78.0)
        cell = CTRcalc.UnitCell(
            [
                lattice.getfloat("a1", 1.0),
                lattice.getfloat("a2", 1.0),
                lattice.getfloat("a3", 1.0),
            ],
            [
                lattice.getfloat("alpha1", 90.0),
                lattice.getfloat("alpha2", 90.0),
                lattice.getfloat("alpha3", 90.0),
            ],
        )
        cell.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
        cell.setEnergy(energy * 1e3)
        ub_calculator = HKLVlieg.UBCalculator(cell, energy)
        ub_calculator.defaultU()

        detector = DetectorCalibration.Detector2D_SXRD()
        detector.setFit2D(
            machine.getfloat("SDD", 0.729) * 1e3,
            machine.getfloat("cpx", 731),
            machine.getfloat("cpy", 1587),
            pixelX=machine.getfloat("pixelsize", 172e-6) * 1e6,
            pixelY=machine.getfloat("pixelsize", 172e-6) * 1e6,
        )
        detector.wavelength = ub_calculator.getLambda() * 1e-10
        detector.detector.shape = (
            int(machine.getfloat("sizey", 3000)),
            int(machine.getfloat("sizex", 3000)),
        )
        detector.detector.max_shape = detector.detector.shape
        detector.setAzimuthalReference(
            np.deg2rad(diffrac.getfloat("azimuthal_reference", 0))
        )
        detector.setPolarization(
            np.deg2rad(diffrac.getfloat("polarization_axis", 0)),
            diffrac.getfloat("polarization_factor", 0),
        )
        # Optional moveable detector arm. Absent keys mean a calibration taken
        # with the arm at zero and no arm motion, i.e. exactly the static
        # geometry every existing config describes.
        detarm = config["DetectorArm"] if config.has_section("DetectorArm") else {}
        detector.setArmReference(
            gamma_arm=np.deg2rad(float(detarm.get("gamma_arm_0", 0.0))),
            delta_arm=np.deg2rad(float(detarm.get("delta_arm_0", 0.0))),
        )
        return cls(
            detector=detector,
            unit_cell=cell,
            ub_calculator=ub_calculator,
            mu=np.deg2rad(diffrac.getfloat("mu", 0.05)),
            chi=np.deg2rad(diffrac.getfloat("chi", 0.0)),
            phi=np.deg2rad(diffrac.getfloat("phi", 0.0)),
            gamma_arm=np.deg2rad(float(detarm.get("gamma_arm", 0.0))),
            delta_arm=np.deg2rad(float(detarm.get("delta_arm", 0.0))),
            arm_angle_frame=str(detarm.get("angle_frame", "prim")),
            refraction_index=1.0 - lattice.getfloat("refractionindex", 0.0),
        )

    readConfig = from_ini

    @classmethod
    def from_gui(cls, gui):
        """Capture config state from an orGUI main window or QUBCalculator."""
        ub_widget = getattr(gui, "ubcalc", gui)
        reflections = []
        if hasattr(gui, "reflectionSel"):
            reflections = list(gui.reflectionSel.reflections)
        cell = unit_cell_from_nxdict(unit_cell_to_nxdict(ub_widget.crystal))
        ub_calculator = HKLVlieg.UBCalculator(cell, ub_widget.ubCal.getEnergy())
        ub_calculator.setU(ub_widget.ubCal.getU())
        corrections = CorrectionState()
        roi = ROIState()
        if hasattr(gui, "scanSelector"):
            options = gui.scanSelector.get_integration_options()
            repair = getattr(getattr(gui, "maskManager", None), "settings", None)
            repair = getattr(repair, "pixel_repair", None)
            repair_enabled = bool(getattr(repair, "enabled", False))
            excluded = getattr(gui, "excludedImagesDialog", None)
            excluded = () if excluded is None else excluded.getData()
            # Read only if the dialog was already created (i.e. opened at
            # least once); a config capture must not build GUI it does not
            # otherwise need. An unopened dialog has nothing to record, which
            # is the same "not recorded" state as a config predating it.
            footprint_dialog = getattr(
                getattr(gui.scanSelector, "correctionsDialog", None),
                "footprintOptions",
                None,
            )
            sample_length_m = None
            sample_width_m = None
            beam_flux_density = None
            beam_shape = {}
            beam_profile_positions_m = ()
            beam_profile_density_per_m = ()
            if footprint_dialog is not None:
                sample_length_m = footprint_dialog.sampleLength()
                sample_width_m = footprint_dialog.sampleWidth()
                beam_flux_density = footprint_dialog.beamFluxDensity()
                # The dialog's own settings() dict, in its own display units;
                # see CorrectionState.beam_shape_values.
                beam_shape = footprint_dialog.settings()
                if not beam_shape.get("analytical", True):
                    try:
                        positions, density = (
                            footprint_dialog.measuredProfile().profile_curve()
                        )
                    except ValueError:
                        # An unresolved imported path stays unresolved. Config
                        # capture must not make an otherwise-unused correction
                        # fatal or invent an empty physical profile.
                        pass
                    else:
                        beam_profile_positions_m = tuple(
                            float(value) for value in positions
                        )
                        beam_profile_density_per_m = tuple(
                            float(value) for value in density
                        )
            corrections = CorrectionState(
                use_mask=bool(options.get("mask", False)) or repair_enabled,
                use_background=getattr(gui, "background_image", None)
                is not None,
                use_solid_angle=bool(options.get("solid_angle", False)),
                use_polarization=bool(options.get("polarization", False)),
                use_lorentz=bool(options.get("lorentz", False)),
                use_footprint=bool(options.get("footprint", False)),
                use_normalization=bool(options.get("normalization", False)),
                repair_masked_pixels=repair_enabled,
                repair_max_component_pixels=getattr(
                    repair, "max_component_pixels", None
                ),
                repair_max_span=getattr(repair, "max_span", None),
                repair_radius=getattr(repair, "radius", None),
                repair_min_valid_neighbors=getattr(
                    repair, "min_valid_neighbors", None
                ),
                repair_use_pyfai_gaps=bool(
                    getattr(repair, "use_pyfai_gaps", True)
                ),
                repair_gap_size_px=int(
                    getattr(repair, "gap_size_px", 1)
                ),
                normalize_exposure=bool(
                    getattr(
                        getattr(gui, "ctr_correction_state", None),
                        "normalize_exposure",
                        getattr(gui, "reconstruction_normalize_exposure", True),
                    )
                ),
                monitor_corrections=tuple(
                    getattr(
                        getattr(gui, "ctr_correction_state", None),
                        "monitor_corrections",
                        getattr(gui, "reconstruction_monitor_corrections", ()),
                    )
                ),
                shared_frame_normalization=True,
                excluded_frames=tuple(
                    sorted(
                        int(value)
                        for value in np.asarray(excluded).ravel()
                        if int(value) >= 0
                    )
                ),
                sample_length_m=sample_length_m,
                sample_width_m=sample_width_m,
                beam_flux_density=beam_flux_density,
                beam_shape_analytical=beam_shape.get("analytical"),
                beam_shape_name=beam_shape.get("shape"),
                beam_shape_values=tuple(beam_shape.get("shape_values", ())),
                beam_profile_file=beam_shape.get("profile_file"),
                beam_profile_content=beam_shape.get("profile_content"),
                beam_profile_unit=beam_shape.get("profile_unit"),
                beam_profile_center=beam_shape.get("profile_center"),
                beam_profile_offset_um=beam_shape.get("profile_offset"),
                beam_profile_positions_m=beam_profile_positions_m,
                beam_profile_density_per_m=beam_profile_density_per_m,
            )
            # Stage 5 has a numerical total-flux path before Stage 6 adds its
            # widgets. Preserve explicitly loaded/programmatic version-3
            # settings on the GUI object instead of dropping them when a new
            # extraction snapshot is captured.
            total_flux_state = getattr(gui, "ctr_correction_state", None)
            if total_flux_state is not None:
                for name in _TOTAL_FLUX_STATE_FIELDS:
                    setattr(corrections, name, getattr(total_flux_state, name))
            roi = ROIState(
                region=dict(options.get("region", {})),
                advanced=dict(options.get("advanced", {})),
                rocking_scan=dict(options.get("rocking_scan", {})),
            )
        return cls(
            detector=ub_widget.detectorCal,
            unit_cell=cell,
            ub_calculator=ub_calculator,
            mu=getattr(ub_widget, "mu", 0.0),
            chi=getattr(ub_widget, "chi", 0.0),
            phi=getattr(ub_widget, "phi", 0.0),
            gamma_arm=getattr(ub_widget, "gamma_arm", 0.0),
            delta_arm=getattr(ub_widget, "delta_arm", 0.0),
            arm_angle_frame=getattr(ub_widget, "arm_angle_frame", "prim"),
            refraction_index=getattr(ub_widget, "n", 1.0),
            reference_reflections=reflections,
            corrections=corrections,
            roi=roi,
        )

    def apply_to_gui(self, gui):
        """Apply this config to an orGUI main window or QUBCalculator."""
        ub_widget = getattr(gui, "ubcalc", gui)
        ub_widget.detectorCal = detector_from_nxdict(detector_to_nxdict(self.detector))
        ub_widget.crystal = unit_cell_from_nxdict(unit_cell_to_nxdict(self.unit_cell))
        # Keep the not-yet-widget-backed total-flux settings available to the
        # extraction policy. Copy through the public JSON representation so a
        # later GUI edit cannot mutate the loaded ConfigData instance.
        gui.ctr_correction_state = CorrectionState.from_dict(
            self.corrections.to_dict()
        )
        ub_widget.ubCal = HKLVlieg.UBCalculator(
            ub_widget.crystal, self.ub_calculator.getEnergy()
        )
        ub_widget.ubCal.setU(self.ub_calculator.getU())
        # calcReflection reads the orientation through ``angles``; it must
        # not keep the replaced calculator (peak positions used the old U).
        ub_widget.angles = HKLVlieg.VliegAngles(ub_widget.ubCal)
        ub_widget.mu = self.mu
        ub_widget.chi = self.chi
        ub_widget.phi = self.phi
        ub_widget.gamma_arm = self.gamma_arm
        ub_widget.delta_arm = self.delta_arm
        ub_widget.arm_angle_frame = self.arm_angle_frame
        ub_widget.n = self.refraction_index
        if hasattr(ub_widget, "uedit"):
            ub_widget.uedit.setU(ub_widget.ubCal.getU())
        if hasattr(ub_widget, "crystalparams"):
            ub_widget.crystalparams.setValues(ub_widget.crystal, ub_widget.n)
        if hasattr(ub_widget, "machineParams"):
            ub_widget.machineParams.setValues(
                {
                    "diffractometer": {
                        "mu": self.mu,
                        "phi": self.phi,
                        "chi": self.chi,
                    },
                    "source": {"E": ub_widget.detectorCal.energy},
                    "SXRD_geometry": ub_widget.detectorCal,
                }
            )
        if hasattr(gui, "reflectionSel"):
            gui.reflectionSel.setReflections(self.reference_reflections)
        if hasattr(gui, "scanSelector"):
            options = {
                "mask": self.corrections.use_mask,
                "solid_angle": self.corrections.use_solid_angle,
                "polarization": self.corrections.use_polarization,
            }
            # Only switches this configuration actually recorded. A file
            # written before they were stored leaves them as the user has
            # them, rather than silently turning a correction off.
            for name, value in (
                ("lorentz", self.corrections.use_lorentz),
                ("footprint", self.corrections.use_footprint),
                ("normalization", self.corrections.use_normalization),
            ):
                if value is not None:
                    options[name] = value
            for name in ("region", "advanced", "rocking_scan"):
                values = getattr(self.roi, name)
                if values:
                    options[name] = dict(values)
            gui.scanSelector.set_integration_options(options)
            # Only if this configuration actually recorded them: a file
            # written before these existed must leave the dialog at whatever
            # the user already has, not reset it to its own defaults.
            if (
                self.corrections.sample_length_m is not None
                or self.corrections.sample_width_m is not None
                or self.corrections.beam_flux_density is not None
                or self.corrections.beam_shape_name is not None
                or self.corrections.beam_profile_file is not None
                or self.corrections.horizontal_interception is not None
            ):
                corrections_dialog = getattr(
                    gui.scanSelector, "correctionsDialog", None
                )
                if corrections_dialog is not None:
                    footprint_dialog = corrections_dialog.footprintOptions_shared()
                    if self.corrections.sample_length_m is not None:
                        footprint_dialog.setSampleLength(
                            self.corrections.sample_length_m
                        )
                    if self.corrections.sample_width_m is not None:
                        footprint_dialog.setSampleWidth(
                            self.corrections.sample_width_m
                        )
                    if self.corrections.beam_flux_density is not None:
                        footprint_dialog.setBeamFluxDensity(
                            self.corrections.beam_flux_density
                        )
                    if self.corrections.horizontal_interception is not None:
                        footprint_dialog.setHorizontalInterception(
                            self.corrections.horizontal_interception,
                            self.corrections.horizontal_intercepted_fraction,
                        )
                    # The dialog's own setSettings(), in its own display
                    # units; only the keys this configuration actually
                    # recorded are passed, so the rest of the dialog is left
                    # alone (setSettings()'s own contract).
                    beam_shape = {}
                    for key, value in (
                        ("analytical", self.corrections.beam_shape_analytical),
                        ("shape", self.corrections.beam_shape_name),
                        ("profile_file", self.corrections.beam_profile_file),
                        ("profile_content", self.corrections.beam_profile_content),
                        ("profile_unit", self.corrections.beam_profile_unit),
                        ("profile_center", self.corrections.beam_profile_center),
                        ("profile_offset", self.corrections.beam_profile_offset_um),
                    ):
                        if value is not None:
                            beam_shape[key] = value
                    if self.corrections.beam_shape_values:
                        beam_shape["shape_values"] = list(
                            self.corrections.beam_shape_values
                        )
                    if beam_shape:
                        footprint_dialog.setSettings(beam_shape)
        gui.reconstruction_normalize_exposure = self.corrections.normalize_exposure
        gui.reconstruction_monitor_corrections = self.corrections.monitor_corrections
        selector = getattr(gui, "scanSelector", None)
        corrections_dialog = getattr(selector, "correctionsDialog", None)
        refresh = getattr(corrections_dialog, "refresh", None)
        if refresh is not None:
            refresh()
        if (
            hasattr(gui, "maskManager")
            and self.corrections.repair_max_component_pixels is not None
        ):
            from .mask_config import PixelRepairSettings

            gui.maskManager.set_pixel_repair_settings(
                PixelRepairSettings(
                    enabled=self.corrections.repair_masked_pixels,
                    max_component_pixels=(
                        self.corrections.repair_max_component_pixels
                    ),
                    max_span=self.corrections.repair_max_span,
                    radius=self.corrections.repair_radius,
                    min_valid_neighbors=(
                        self.corrections.repair_min_valid_neighbors
                    ),
                    use_pyfai_gaps=(
                        self.corrections.repair_use_pyfai_gaps
                    ),
                    gap_size_px=self.corrections.repair_gap_size_px,
                )
            )
        if hasattr(gui, "excludedImagesDialog"):
            excluded = np.asarray(
                self.corrections.excluded_frames or (-1,), dtype=np.int64
            )
            gui.excludedImagesDialog.updateArrayData(excluded)
        if hasattr(ub_widget, "updateReflectionMismatch"):
            ub_widget.updateReflectionMismatch()
        # Every interactive path that changes mu/chi/phi, the UB matrix or
        # the detector geometry (_onMachineParamsChanged, _onCrystalParamsChanged,
        # _onAlignU) emits these right after updateReflectionMismatch(), which
        # documents exactly that contract; skipping it here left the ROI and
        # reflection overlays, and the Q-plot, showing the geometry from
        # before the config was loaded; the angle readout in the caller
        # (`orgui.mu`, etc.) was already correct -- only the cached HKL and
        # its dependent plot elements were stale.
        if hasattr(ub_widget, "sigPlottableMachineParamsChanged"):
            ub_widget.sigPlottableMachineParamsChanged.emit()
        if hasattr(ub_widget, "sigReplotRequest"):
            ub_widget.sigReplotRequest.emit(True)

    def to_nxdict(self, role="scan", source=None):
        """Convert this config to a NeXus-compatible nested dictionary."""
        nxdict = {
            "@NX_class": "NXcollection",
            "@orgui_meta": "config",
            "@orgui_config_role": role,
            "@orgui_schema_version": SCHEMA_VERSION,
            # datetime.timezone.utc, not the datetime.UTC alias, which only
            # exists from Python 3.11 while orGUI supports 3.10.
            "@orgui_config_created": datetime.datetime.now(
                datetime.timezone.utc
            ).isoformat(),
            "instrument": {
                "@NX_class": "NXinstrument",
                "detector_SXRD": detector_to_nxdict(self.detector),
            },
            "sample": unit_cell_to_nxdict(self.unit_cell),
            "reference_reflections": reflections_to_nxdict(
                self.reference_reflections, self.ub_calculator
            ),
            "orgui": {
                "@NX_class": "NXcollection",
                "diffractometer": {
                    "@NX_class": "NXcollection",
                    "mu": self.mu,
                    "chi": self.chi,
                    "phi": self.phi,
                    "gamma_arm": self.gamma_arm,
                    "delta_arm": self.delta_arm,
                    "arm_angle_frame": self.arm_angle_frame,
                    "@unit": "rad",
                },
                "source": {
                    "@NX_class": "NXcollection",
                    "energy": self.ub_calculator.getEnergy(),
                    "@energy_unit": "keV",
                    "wavelength": self.ub_calculator.getLambda(),
                    "@wavelength_unit": "Angstrom",
                },
                "refraction_index": self.refraction_index,
                "integration_corrections": corrections_to_nxdict(
                    self.corrections
                ),
                "roi_integration": roi_to_nxdict(self.roi),
                **self.orgui,
            },
        }
        nxdict["sample"]["orientation_matrix"] = self.ub_calculator.getU()
        nxdict["sample"]["ub_matrix"] = self.ub_calculator.getUB()
        if source is not None:
            nxdict["@orgui_config_source"] = source
        return nxdict

    @classmethod
    def from_nxdict(cls, nxdict):
        """Load config state from a NeXus-compatible dictionary."""
        detector = detector_from_nxdict(nxdict["instrument"]["detector_SXRD"])
        unit_cell = unit_cell_from_nxdict(nxdict["sample"])
        source = nxdict.get("orgui", {}).get("source", {})
        energy = float(source.get("energy", detector.energy))
        ub_calculator = HKLVlieg.UBCalculator(unit_cell, energy)
        ub_calculator.setU(np.asarray(nxdict["sample"]["orientation_matrix"]))
        diffrac = nxdict.get("orgui", {}).get("diffractometer", {})
        corrections_group = nxdict.get("orgui", {}).get(
            "integration_corrections", {}
        ) or {}
        if "json" in corrections_group:
            # Configurations written before the typed layout. Still read, so
            # that existing databases keep loading; never written any more.
            corrections = CorrectionState.from_dict(
                json.loads(_as_text(_plain(corrections_group["json"])) or "{}")
            )
        else:
            corrections = corrections_from_nxdict(corrections_group)
        roi = roi_from_nxdict(nxdict.get("orgui", {}).get("roi_integration"))
        return cls(
            detector=detector,
            unit_cell=unit_cell,
            ub_calculator=ub_calculator,
            mu=float(diffrac.get("mu", 0.0)),
            chi=float(diffrac.get("chi", 0.0)),
            phi=float(diffrac.get("phi", 0.0)),
            gamma_arm=float(diffrac.get("gamma_arm", 0.0)),
            delta_arm=float(diffrac.get("delta_arm", 0.0)),
            arm_angle_frame=_as_text(diffrac.get("arm_angle_frame", "prim")),
            refraction_index=float(
                nxdict.get("orgui", {}).get("refraction_index", 1.0)
            ),
            reference_reflections=reflections_from_nxdict(nxdict),
            corrections=corrections,
            roi=roi,
        )

    def to_json_dict(self) -> dict[str, Any]:
        """Serialize this configuration through the central NeXus schema."""
        return _json_value(self.to_nxdict(role="reconstruction"))

    @classmethod
    def from_json_dict(cls, values):
        """Deserialize a central configuration JSON dictionary."""
        return cls.from_nxdict(dict(values))


class ConfigHandler:
    """Read and write orGUI config groups in an HDF5 database."""

    def __init__(self, gui=None, create_dataset_args=None):
        self.gui = gui
        self.create_dataset_args = create_dataset_args or {}

    @staticmethod
    def is_config_group(group):
        """Return ``True`` for loadable orGUI config HDF5 groups."""
        return (
            isinstance(group, h5py.Group)
            and group.attrs.get("orgui_meta") == "config"
        )

    def write_scan_config(self, scan_group, config, source="scan_import"):
        """Write ``configuration`` under a scan group."""
        return self._write_config(scan_group, config, "scan", source)

    def write_integration_config(
        self, integration_group, config, source="integration_save"
    ):
        """Write ``configuration`` under an integration result group."""
        return self._write_config(integration_group, config, "integration", source)

    def _write_config(self, group, config, role, source):
        if "configuration" in group:
            del group["configuration"]
        dicttonx(
            {"configuration": config.to_nxdict(role=role, source=source)},
            group,
            update_mode="add",
            create_dataset_args=self.create_dataset_args,
        )
        return group["configuration"]

    def load_config_group(self, group):
        """Read a marked config group and return :class:`ConfigData`."""
        if not self.is_config_group(group):
            raise ValueError(f"Not an orGUI config group: {group.name}")
        return ConfigData.from_nxdict(nxtodict(group))

    def apply_config_group(self, group, gui=None):
        """Load a marked config group and apply it to a GUI object."""
        config = self.load_config_group(group)
        target = gui or self.gui
        if target is None:
            raise ValueError("No GUI target is configured for loading config groups.")
        config.apply_to_gui(target)
        return config

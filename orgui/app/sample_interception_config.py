"""Resolve saved shape/profile settings at the application unit boundary.

Settings explicitly declare a motor as rotation about the surface normal; no
SPEC name or six-circle solver angle is inferred to mean that rotation.
Profile offsets align the shape origin at the configured reference incidence.
The resulting rotation-axis alignment is fixed in laboratory beam coordinates.
"""

import copy
import logging

import numpy as np

from ..backend.scans import sample_azimuth
from ..datautils.xrayutils.corrections import beamprofile as bp
from ..datautils.xrayutils.corrections.sample_interception import SampleShape, overlap

logger = logging.getLogger(__name__)


def profile_from_settings(settings):
    """Build a normalized profile from dialog-unit settings and embedded data.

    Widths/offsets are micrometres; dimensionless shape values stay dimensionless.
    Embedded coordinates are metres relative to the configured sample centre.
    """
    center = settings.get("profile_center", "centroid")
    offset = float(settings.get("profile_offset", 0) or 0) * 1e-6
    if settings.get("analytical", True):
        factories = {
            "Gaussian": (bp.gaussian_profile, (1e-6,)),
            "Top hat": (bp.top_hat_profile, (1e-6,)),
            "Trapezoid": (bp.trapezoid_profile, (1e-6, 1e-6)),
            "Smoothed top hat": (bp.smoothed_top_hat_profile, (1e-6, 1e-6)),
            "Generalized normal": (bp.generalized_normal_profile, (1e-6, 1)),
            "Skew normal": (bp.skew_normal_profile, (1e-6, 1)),
        }
        name = settings.get("shape")
        if name not in factories:
            raise ValueError(f"unknown analytical beam shape {name!r}")
        factory, scales = factories[name]
        values = settings.get("shape_values")
        if values is None or len(values) != len(scales):
            raise ValueError("beam shape parameters do not match its shape")
        return factory(
            *[v * s for v, s in zip(values, scales)], center=center, offset=offset
        )
    if settings.get("positions_m") is not None:
        positions = np.asarray(settings["positions_m"], dtype=float)
        density = np.asarray(settings.get("density_per_m", ()), dtype=float)
        # Embedded positions already include reference/offset. Reconstruct
        # that exact origin instead of centering the same data a second time.
        profile = bp.MeasuredBeamProfile(positions, density)
        return bp.MeasuredBeamProfile(positions, density, offset=-profile.sample_center)
    path = settings.get("profile_file")
    if not path:
        raise ValueError("a measured beam requires embedded data or a profile file")
    unit = settings.get("profile_unit", "mm")
    if unit not in {"mm", "microns", "micron", "um"}:
        raise ValueError("profile coordinates must be mm or micrometres")
    positions, density = bp.read_profile_file(
        path,
        z_scale=1e-3 if unit == "mm" else 1e-6,
        height_scan="scan" in settings.get("profile_content", "").lower(),
    )
    return bp.MeasuredBeamProfile(positions, density, center=center, offset=offset)


def vertical_settings(state):
    """Return existing vertical settings, preserving embedded centred coordinates."""
    settings = {
        "analytical": state.beam_shape_analytical,
        "shape": state.beam_shape_name,
        "shape_values": list(state.beam_shape_values),
        "profile_file": state.beam_profile_file,
        "profile_content": state.beam_profile_content,
        "profile_unit": state.beam_profile_unit,
        "profile_center": state.beam_profile_center,
        "profile_offset": state.beam_profile_offset_um,
    }
    if state.beam_profile_positions_m:
        settings["positions_m"] = list(state.beam_profile_positions_m)
        settings["density_per_m"] = list(state.beam_profile_density_per_m)
    return {k: v for k, v in settings.items() if v is not None}


def embed_profile(settings, profile):
    """Copy profile settings with measured normalized data relative to its centre."""
    settings = copy.deepcopy(settings)
    if isinstance(profile, bp.MeasuredBeamProfile):
        positions, density = profile.profile_curve()
        settings["positions_m"] = positions.tolist()
        settings["density_per_m"] = density.tolist()
    return settings


def shape_from_settings(settings):
    """Build a surface shape from saved SI dimensions/vertices."""
    if settings.get("version") != 1:
        raise ValueError("unsupported sample-interception settings version")
    shape = settings.get("shape", {})
    return SampleShape(
        shape.get("kind"),
        tuple(shape.get("dimensions_m", ())),
        tuple(map(tuple, shape.get("vertices_m", ()))),
    )


def shape_frame_factors(settings, vertical, alpha, azimuth, *, progress=None):
    """Resolve shape H, f_hit and peak-referenced area fraction per frame.

    ``alpha`` and raw ``azimuth`` are radians. Profile offsets align the shape
    origin at ``reference_incidence_deg``; the laboratory axis stays fixed.
    Invalid angles/zero overlap are recorded as NaN divisors, never as unity.
    """
    shape = shape_from_settings(settings)
    if not settings.get("normal_rotation_confirmed", False):
        raise ValueError(
            "declare the azimuth source to be rotation about the surface normal"
        )
    horizontal_settings = settings.get("horizontal")
    if not horizontal_settings:
        raise ValueError("exact 2D interception requires a horizontal beam profile")
    horizontal = profile_from_settings(horizontal_settings)
    alpha, azimuth = np.broadcast_arrays(alpha, azimuth)
    sign = settings.get("azimuth_sign", 1)
    if sign not in {-1, 1}:
        raise ValueError("azimuth sign must be +1 or -1")
    psi = sign * (azimuth - np.deg2rad(settings.get("reference_azimuth_deg", 0)))
    offset = tuple(settings.get("offset_m", (0, 0)))
    if len(offset) != 2:
        raise ValueError("shape offset must have two coordinates, in metres")
    reference_alpha = np.deg2rad(settings.get("reference_incidence_deg", 0))
    if not np.isfinite(reference_alpha) or not 0 < reference_alpha <= np.pi / 2:
        raise ValueError("reference incidence must be in (0, 90] degrees")
    valid = np.isfinite(alpha) & (alpha > 0) & (alpha <= np.pi / 2) & np.isfinite(psi)
    if not valid.any():
        raise ValueError(
            "shape interception needs at least one frame with valid angles"
        )
    illumination = np.full(alpha.shape, np.nan)
    fraction = np.full(alpha.shape, np.nan)
    area = np.full(alpha.shape, np.nan)
    error = np.full(alpha.shape, np.nan)
    result = overlap(
        shape,
        vertical,
        horizontal,
        alpha[valid],
        psi[valid],
        offset=offset,
        orientation=np.deg2rad(settings.get("orientation_deg", 0)),
        axis_vertical=-offset[0] * np.sin(reference_alpha),
        axis_horizontal=-offset[1],
        rtol=float(settings.get("rtol", 1e-9)),
        progress=progress,
    )
    illumination[valid], fraction[valid] = result.illumination, result.fraction
    area[valid], error[valid] = result.effective_area, result.error
    valid &= fraction > 0
    illumination[~valid], area[~valid] = np.nan, np.nan
    if not valid.any():
        raise ValueError("no frame has positive sample/beam overlap")
    threshold = float(settings.get("overspill_threshold", 0.01))
    if not 0 <= threshold <= 1:
        raise ValueError("overspill threshold must be between zero and one")
    if np.any(~valid):
        logger.warning(
            "Shape interception excludes %d frame values (NaN divisors)",
            np.count_nonzero(~valid),
        )
    if np.any(1 - fraction[valid] > threshold):
        logger.warning(
            "Sample/beam overspill exceeds %.3g%% in %d frame values",
            100 * threshold,
            np.count_nonzero(1 - fraction[valid] > threshold),
        )
    return illumination, fraction, area / shape.area, valid, error


def shape_policy_inputs(scan, state, alpha, *, count=None):
    """Return measured frame azimuth [rad] and self-contained shape settings."""
    settings = copy.deepcopy(state.sample_interception)
    if count is None:
        count = np.asarray(alpha).shape[-1] if np.ndim(alpha) else len(scan)
    azimuth = sample_azimuth(scan, settings, count)
    return azimuth, settings


def replacement_shape_divisor(settings, vertical, record):
    """Rebuild a shape divisor using saved acquisition angles, never solver angles.

    Changing a moving readback source requires re-extraction. Sign/reference
    and geometry may change from the immutable base using the same source.
    """
    import json

    if record.alpha is None:
        raise ValueError("saved incidence is required for shape replacement")
    alpha = np.asarray(record.alpha)
    source = settings.get("azimuth_source")
    if source == "fixed":
        azimuth = np.deg2rad(settings["fixed_azimuth_deg"])
    else:
        old = json.loads(
            record.profile_provenance.get("sample_interception_json", "{}")
        )
        if (
            source != old.get("azimuth_source")
            or settings.get("azimuth_unit", "deg") != old.get("azimuth_unit", "deg")
            or "sample_azimuth_rad" not in record.profile_provenance
        ):
            raise ValueError("changing a moving azimuth source requires re-extraction")
        azimuth = np.asarray(record.profile_provenance["sample_azimuth_rad"])
    factors = shape_frame_factors(settings, vertical, alpha, azimuth)
    return (
        factors[0] if record.scale_convention.startswith("total_flux_") else factors[2]
    )

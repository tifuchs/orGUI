"""Resolve saved shape/profile settings at the application unit boundary.

Automatic angles use acquisition mu/omega, with explicit overrides available.
The flat-surface model treats azimuth as rotation about the surface normal.
Profile offsets align the shape origin at the configured reference incidence.
The resulting rotation-axis alignment is fixed in laboratory beam coordinates.
"""

import copy
import logging
from hashlib import sha256
from io import BytesIO
from pathlib import Path

import numpy as np

from ..backend.scans import sample_azimuth, sample_incidence
from ..datautils.xrayutils.corrections import beamprofile as bp
from ..datautils.xrayutils.corrections.sample_interception import SampleShape, overlap

logger = logging.getLogger(__name__)


def _alignment_axis_offset(settings):
    shape = settings.get("shape", {})
    dimensions = shape.get("dimensions_m", ())
    if (shape.get("kind") == "rectangle" and len(dimensions) == 2
            and dimensions[1] > dimensions[0]):
        return 90
    return 0


def parallel_omega_reading(settings):
    """Return omega [deg] when the long edge/polygon x axis is parallel.

    For ``edge_angle = sign * (omega - reference) + orientation``, the
    parallel reading is ``reference - sign * orientation`` (sign is +/-1).
    Do not replace the placement reference with this value for off-centre
    samples: their origin displacement rotates about the original reference.
    No periodic wrapping is applied, preserving polygon and motor conventions.
    If rectangle width exceeds length, its long edge is the local y direction;
    add 90 degrees to the old x-axis orientation before converting.
    """
    sign = settings.get("azimuth_sign", 1)
    if sign not in {-1, 1}:
        raise ValueError("azimuth sign must be +1 or -1")
    return (settings.get("reference_azimuth_deg", 0)
            - sign * (settings.get("orientation_deg", 0)
                      + _alignment_axis_offset(settings)))


def orientation_for_parallel_reading(settings, reading):
    """Return mounting orientation [deg] for a parallel omega reading [deg].

    Keep the existing placement reference and offset coordinates unchanged.
    This is the inverse of :func:`parallel_omega_reading` for either sign.
    """
    sign = settings.get("azimuth_sign", 1)
    if sign not in {-1, 1}:
        raise ValueError("azimuth sign must be +1 or -1")
    return (sign * (settings.get("reference_azimuth_deg", 0) - reading)
            - _alignment_axis_offset(settings))


def _profile_file_data(settings, base_path=None):
    path = settings.get("profile_file")
    if not path:
        raise ValueError("a measured beam requires embedded data or a profile file")
    path = Path(path)
    base = settings.get("profile_base")
    if base is not None and not path.is_absolute():
        base = Path(base)
        if not base.is_absolute():
            base = Path(base_path or Path.cwd()) / base
        path = base / path
    try:
        content = path.read_bytes()
    except OSError as error:
        raise ValueError(
            f"Cannot read measured beam profile {path}: {error}"
        ) from error
    expected = settings.get("profile_sha256")
    if expected is not None and sha256(content).hexdigest() != expected:
        raise ValueError(f"Measured beam profile content changed: {path}")
    return content


def profile_file_identity(settings, *, base_path=None):
    """Return SHA-256 of a measured profile, verifying any saved identity.

    :param dict settings: File settings with optional ``profile_base``.
    :param base_path: Directory relative bases resolve against (job JSON parent).
        Without an explicit base, historical working-directory resolution stays.
    :raises ValueError: If the source is missing or its saved identity differs.
    """
    return sha256(_profile_file_data(settings, base_path)).hexdigest()


def profile_from_settings(settings, *, base_path=None):
    """Build a normalized profile from dialog-unit settings and embedded data.

    Widths/offsets are micrometres; dimensionless shape values stay dimensionless.
    Embedded coordinates are metres relative to the configured sample centre.
    ``profile_base`` resolves relative to ``base_path`` (the job JSON parent),
    or the working directory outside a job. ``profile_sha256`` checks the exact
    bytes consumed. Embedded data take precedence over display-only filenames.
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
    unit = settings.get("profile_unit", "mm")
    if unit not in {"mm", "microns", "micron", "um"}:
        raise ValueError("profile coordinates must be mm or micrometres")
    positions, density = bp.read_profile_file(
        BytesIO(_profile_file_data(settings, base_path)),
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
        "profile_base": getattr(state, "beam_profile_base", None),
        "profile_sha256": getattr(state, "beam_profile_sha256", None),
        "profile_storage": getattr(state, "beam_profile_storage", None),
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


def validate_beam_profile(profile):
    """Validate normalized density [1/m] without evaluating sample overlap.

    :raises ValueError: For a nonpositive peak or negative/nonfinite density.
    """
    if not np.isfinite(profile.peak_density) or profile.peak_density <= 0:
        raise ValueError("beam profile must have finite positive peak density")
    density = profile.density_at(profile.integration_points)
    if np.any(density < 0) or not np.all(np.isfinite(density)):
        raise ValueError("beam density must be finite and non-negative")


def _shape_frame_inputs(settings, alpha, azimuth):
    """Validate surface geometry [m/rad] without computing frame factors."""
    shape = shape_from_settings(settings)
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
    if "reference_incidence_deg" not in settings:
        candidates = alpha[np.isfinite(alpha) & (alpha > 0) & (alpha <= np.pi / 2)]
        if candidates.size:
            reference_alpha = candidates.flat[0]
    if not np.isfinite(reference_alpha) or not 0 < reference_alpha <= np.pi / 2:
        raise ValueError("reference incidence must be in (0, 90] degrees")
    valid = np.isfinite(alpha) & (alpha > 0) & (alpha <= np.pi / 2) & np.isfinite(psi)
    if not valid.any():
        raise ValueError(
            "shape interception needs at least one frame with valid angles"
        )
    geometry = {
        "offset": offset,
        "orientation": np.deg2rad(settings.get("orientation_deg", 0)),
        "axis_vertical": -offset[0] * np.sin(reference_alpha),
        "axis_horizontal": -offset[1],
        "rtol": float(settings.get("rtol", 1e-9)),
    }
    if not np.all(np.isfinite([*offset, geometry["orientation"],
                               settings.get("reference_azimuth_deg", 0)])):
        raise ValueError("finite angles and finite placement are required")
    if not 1e-12 <= geometry["rtol"] <= 1e-3:
        raise ValueError("integration rtol must be between 1e-12 and 1e-3")
    if not 0 <= float(settings.get("overspill_threshold", 0.01)) <= 1:
        raise ValueError("overspill threshold must be between zero and one")
    return shape, horizontal, alpha, psi, valid, geometry


def resolve_shape_inputs(scan, settings, alpha, *, count=None, omega=None,
                         base_path=None):
    """Validate and freeze shape inputs without footprint quadrature.

    Acquisition ``alpha`` and optional ``omega`` are radians. Returned settings
    retain SI shape dimensions, the first physical reference incidence in
    degrees, and embedded horizontal measured data with its aligned origin.

    :param base_path: Job directory for explicitly relative profile references.
    :returns: Resolved incidence/azimuth arrays [rad] and copied shape settings.
    :raises ValueError: For invalid profiles, geometry or unavailable angles.
    """
    alpha, azimuth, settings = sample_angle_inputs(
        scan, settings, alpha, count=count, omega=omega
    )
    if not settings.get("horizontal"):
        raise ValueError("exact 2D interception requires a horizontal beam profile")
    horizontal = profile_from_settings(settings["horizontal"], base_path=base_path)
    validate_beam_profile(horizontal)
    settings["horizontal"] = embed_profile(settings["horizontal"], horizontal)
    _shape_frame_inputs(settings, alpha, azimuth)
    return alpha, azimuth, settings


def shape_frame_factors(settings, vertical, alpha, azimuth, *, progress=None):
    """Resolve shape H, f_hit and peak-referenced area fraction per frame.

    ``alpha`` and raw ``azimuth`` are radians. Profile offsets align the shape
    origin at ``reference_incidence_deg``; the laboratory axis stays fixed.
    The shape angle is ``orientation_deg + sign * (azimuth - reference)``;
    ``orientation_deg`` is the UI rotation offset. Only the readback difference
    rotates the sample-fixed origin displacement, preserving its orbit.
    Invalid angles/zero overlap are recorded as NaN divisors, never as unity.
    """
    shape, horizontal, alpha, psi, valid, geometry = _shape_frame_inputs(
        settings, alpha, azimuth
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
        progress=progress,
        **geometry,
    )
    illumination[valid], fraction[valid] = result.illumination, result.fraction
    area[valid], error[valid] = result.effective_area, result.error
    valid &= fraction > 0
    illumination[~valid], area[~valid] = np.nan, np.nan
    if not valid.any():
        raise ValueError("no frame has positive sample/beam overlap")
    threshold = float(settings.get("overspill_threshold", 0.01))
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


def shape_policy_inputs(scan, state, alpha, *, count=None, omega=None):
    """Return measured frame azimuth [rad] and self-contained shape settings."""
    _, azimuth, settings = sample_angle_inputs(
        scan, state.sample_interception, alpha, count=count, omega=omega
    )
    return azimuth, settings


def sample_angle_inputs(scan, settings, alpha, *, count=None, omega=None):
    """Return footprint incidence/azimuth [rad] and resolved shape settings.

    Automatic alignment uses the first valid frame incidence, saved explicitly
    for later replacement. Existing numeric reference incidences are retained.
    """
    settings = copy.deepcopy(settings)
    if count is None:
        count = np.asarray(alpha).shape[-1] if np.ndim(alpha) else len(scan)
    alpha = sample_incidence(scan, settings, count, alpha)
    azimuth = sample_azimuth(scan, settings, count, omega=omega)
    if "reference_incidence_deg" not in settings:
        valid = alpha[np.isfinite(alpha) & (alpha > 0) & (alpha <= np.pi / 2)]
        if valid.size:
            settings["reference_incidence_deg"] = float(np.rad2deg(valid.flat[0]))
    return alpha, azimuth, settings


def replacement_shape_divisor(settings, vertical, record):
    """Rebuild a shape divisor using saved acquisition angles, never solver angles.

    Changing a moving readback source requires re-extraction. Sign/reference
    and geometry may change from the immutable base using the same source.
    """
    import json

    if record.alpha is None:
        raise ValueError("saved incidence is required for shape replacement")
    old = json.loads(
        record.profile_provenance.get("sample_interception_json", "{}")
    )
    alpha = np.asarray(record.profile_provenance.get(
        "sample_incidence_rad", record.alpha
    ))
    incidence_source = settings.get("incidence_source", "auto")
    if incidence_source == "fixed":
        alpha = np.deg2rad(settings["fixed_incidence_deg"])
    elif (incidence_source != old.get("incidence_source", "auto")
          or (incidence_source != "auto"
              and settings.get("incidence_unit", "deg")
              != old.get("incidence_unit", "deg"))):
        raise ValueError("changing a moving incidence source requires re-extraction")
    settings = copy.deepcopy(settings)
    if "reference_incidence_deg" not in settings:
        if "reference_incidence_deg" in old:
            settings["reference_incidence_deg"] = old["reference_incidence_deg"]
    source = settings.get("azimuth_source", "auto")
    if source == "fixed":
        azimuth = np.deg2rad(settings["fixed_azimuth_deg"])
    else:
        if (
            source != old.get("azimuth_source", "auto")
            or (source != "auto" and settings.get("azimuth_unit", "deg")
                != old.get("azimuth_unit", "deg"))
            or "sample_azimuth_rad" not in record.profile_provenance
        ):
            raise ValueError("changing a moving azimuth source requires re-extraction")
        azimuth = np.asarray(record.profile_provenance["sample_azimuth_rad"])
    factors = shape_frame_factors(settings, vertical, alpha, azimuth)
    return (
        factors[0] if record.scale_convention.startswith("total_flux_") else factors[2]
    )

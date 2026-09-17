# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
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
r"""Adapter between the integration workflow and the correction physics.

Every correction factor itself lives in
:mod:`orgui.datautils.xrayutils.corrections`, which takes numbers and returns
numbers and knows nothing about scans, configurations or widgets. What is
left here is the translation in the other direction: turning *which
corrections were switched on* and *what a loaded scan calls its counters*
into arguments for those functions, and assembling the result into the bundle
that is stored beside the integrated intensity.

The legacy normalization and numerical active-area factor are intensity
divisors:

.. math::

    I_\mathrm{corr} = \frac{I}{C_\mathrm{norm}\,C_\mathrm{area}}
    \qquad
    F^2_{hkl} = \frac{I_\mathrm{corr}}{L_\mathrm{stationary}}

``C_flux_on_sample`` is retained as a diagnostic component of
``C_illum_area``. It must not be divided out separately: the numerical active
area already contains the same beam/sample overlap integral. Which angular
factors a stationary measurement applies -- and that it has no
rod-interception factor, unlike a rocking scan -- is decided by
:func:`~orgui.datautils.xrayutils.corrections.measurement.mode_components`,
the one place that mapping exists.

The exposure and monitor normalization mirrors the reciprocal-space
reconstruction (:mod:`orgui.reconstruction_job`), so a stationary integration
and a reconstruction of the same scan are normalized identically.

An explicit total-flux configuration instead stores incident photons per frame
``Q`` and the dimensionless illumination divisor ``H``. That versioned path
does not reinterpret the legacy density/area fields, and reversible changes
always restart from the stored polarization-only scalar curve.

This module holds no Qt state and reads only public scan attributes, so it is
safe in CLI and batch use.
"""

from dataclasses import dataclass

import numpy as np

from ..datautils.xrayutils.corrections import activearea, measurement
from ..datautils.xrayutils.corrections.normalization import (
    frame_fluence,
    normalization_divisor as _divisor_from_counters,
    relative_frame_fluence,
)
from ..datautils.xrayutils.corrections.roi import (  # noqa: F401
    CorrectionFactors,
    roi_mean_correction,
)

__all__ = [
    "CorrectionFactors",
    "FrameCorrectionPolicy",
    "FOOTPRINT_APPLY",
    "FOOTPRINT_KEEP",
    "FOOTPRINT_REMOVE",
    "apply_stationary_corrections",
    "corrected_curve_from_record",
    "frame_correction_policy",
    "monitor_counter_candidates",
    "normalization_divisor",
    "pixel_correction_branches",
    "roi_mean_correction",
    "stationary_correction_factors",
    "structure_factor",
    "structure_factor_from_policy",
]

FOOTPRINT_KEEP = "keep"
FOOTPRINT_APPLY = "apply"
FOOTPRINT_REMOVE = "remove"


@dataclass(frozen=True)
class FrameCorrectionPolicy:
    """Exact framewise normalization and illumination decision.

    ``new_contract`` distinguishes an explicitly configured primary-monitor
    or total-flux convention from the legacy product of reconstruction
    monitor counters. Optional arrays remain ``None`` when deliberately not
    applied; statuses record that decision without inventing unity data.
    """

    new_contract: bool
    scale_convention: str
    normalization_status: str
    normalization_divisor: object = None
    normalization_unit: str | None = None
    normalization_components: tuple[str, ...] = ()
    illumination_status: str = "not_applied"
    illumination_divisor: object = None
    illumination_convention: str | None = None
    vertical_intercepted_fraction: object = None
    horizontal_intercepted_fraction: object = None
    intercepted_fraction: object = None
    calibrated: bool = False

    @property
    def ctr_scale_ready(self):
        """Whether both required :math:`Q` and :math:`H` were applied."""
        return (
            self.normalization_status == "applied"
            and self.illumination_status == "applied"
        )


def _explicit_total_flux_contract(state):
    """Return whether ``state`` explicitly selects the new convention."""
    return any(
        getattr(state, name, None) is not None
        for name in (
            "total_incident_flux",
            "total_flux_calibrated",
            "primary_monitor",
            "primary_monitor_kind",
            "horizontal_interception",
            "horizontal_intercepted_fraction",
        )
    )


def _horizontal_fraction(state):
    """Resolve an explicit horizontal-interception setting."""
    mode = getattr(state, "horizontal_interception", None)
    if mode == "full":
        return 1.0
    if mode == "fraction":
        value = getattr(state, "horizontal_intercepted_fraction", None)
        if value is None:
            raise ValueError(
                "horizontal_intercepted_fraction is required for fraction mode"
            )
        return float(value)
    raise ValueError(
        "total-flux illumination requires horizontal_interception to be "
        "'full' or 'fraction'"
    )


def frame_correction_policy(
    scan,
    state,
    size,
    *,
    use_normalization,
    use_illumination,
    alpha=None,
    beam_profile=None,
    sample_length=None,
):
    r"""Resolve the exact framewise :math:`Q` and :math:`H` policy.

    Explicit primary-monitor/total-flux fields select the new convention.
    Otherwise the established exposure-times-multiple-monitors and legacy
    numerical active-area convention is retained. A calibrated flux is used
    only when marked calibrated and every required reference value exists.

    :param scan: Loaded scan providing exposure and counter attributes.
    :param state: Correction-state-like object with the version-3 fields.
    :param int size: Number of frames.
    :param bool use_normalization: Apply the requested frame normalization.
    :param bool use_illumination: Apply the requested footprint correction.
    :param alpha: Incidence angle(s), radian; required for illumination.
    :param beam_profile: Vertical beam profile; required for illumination.
    :param sample_length: Sample length along the beam, meter.
    :returns: Fully resolved policy with the arrays actually applied.
    :rtype: FrameCorrectionPolicy
    """
    new_contract = _explicit_total_flux_contract(state)
    if not new_contract:
        divisor = None
        components = ()
        norm_status = "not_applied"
        if use_normalization:
            divisor, applied = normalization_divisor(
                scan,
                bool(getattr(state, "normalize_exposure", True)),
                tuple(getattr(state, "monitor_corrections", ()) or ()),
                size,
            )
            components = tuple(applied)
            norm_status = "applied" if applied else "unavailable"

        illumination = vertical = None
        illum_status = "not_applied"
        if use_illumination:
            if alpha is None or beam_profile is None or sample_length is None:
                raise ValueError(
                    "legacy footprint correction needs alpha, beam profile "
                    "and sample length"
                )
            vertical, illumination = beam_profile.corrections(
                alpha, sample_length
            )
            illum_status = "applied"
        return FrameCorrectionPolicy(
            new_contract=False,
            scale_convention="legacy_relative",
            normalization_status=norm_status,
            normalization_divisor=divisor,
            normalization_unit=(
                "legacy_counter_product" if divisor is not None else None
            ),
            normalization_components=components,
            illumination_status=illum_status,
            illumination_divisor=illumination,
            illumination_convention=(
                "legacy_C_illum_area" if illumination is not None else None
            ),
            vertical_intercepted_fraction=vertical,
        )

    calibrated = bool(getattr(state, "total_flux_calibrated", False))
    normalization = None
    components = ()
    norm_status = "not_applied"
    norm_unit = None
    if use_normalization:
        exposure = getattr(scan, "exposure_time", None)
        monitor_name = getattr(state, "primary_monitor", None)
        monitor_kind = getattr(state, "primary_monitor_kind", None)
        monitor = None
        if monitor_name is not None:
            if not hasattr(scan, monitor_name):
                raise ValueError(
                    f"Active scan has no primary monitor named {monitor_name!r}"
                )
            monitor = getattr(scan, monitor_name)
        if calibrated:
            total_flux = getattr(state, "total_incident_flux", None)
            if total_flux is None:
                raise ValueError(
                    "a calibrated total-flux correction requires "
                    "total_incident_flux"
                )
            normalization = frame_fluence(
                total_flux,
                exposure_time=exposure,
                monitor=monitor,
                monitor_kind=monitor_kind,
                reference_monitor=getattr(
                    state, "monitor_reference_reading", None
                ),
                reference_exposure=getattr(
                    state, "monitor_reference_exposure_s", None
                ),
            )
            norm_unit = "photons"
        else:
            normalization = relative_frame_fluence(
                exposure_time=exposure,
                monitor=monitor,
                monitor_kind=monitor_kind,
            )
            norm_unit = "relative_fluence"
        normalization = np.broadcast_to(
            np.asarray(normalization, dtype=np.float64), (int(size),)
        ).copy()
        norm_status = "applied"
        if monitor_name is None:
            components = ("exposure",)
        elif monitor_kind == "rate":
            components = (
                "exposure",
                f"primary_monitor:{monitor_name}:{monitor_kind}",
            )
        else:
            components = (f"primary_monitor:{monitor_name}:{monitor_kind}",)

    illumination = vertical = intercepted = None
    horizontal = None
    illum_status = "not_applied"
    if use_illumination:
        if alpha is None or beam_profile is None or sample_length is None:
            raise ValueError(
                "total-flux illumination needs alpha, beam profile and "
                "sample length"
            )
        horizontal = _horizontal_fraction(state)
        vertical = np.asarray(
            beam_profile.flux_on_sample(alpha, sample_length), dtype=np.float64
        )
        illumination = activearea.illumination_divisor(
            alpha,
            sample_length,
            beam_profile,
            horizontal_fraction=horizontal,
        )
        intercepted = activearea.intercepted_fraction(
            alpha,
            sample_length,
            beam_profile,
            horizontal_fraction=horizontal,
        )
        illum_status = "applied"

    if calibrated:
        scale = "total_flux_calibrated"
    else:
        scale = "total_flux_relative"
    if norm_status != "applied" or illum_status != "applied":
        scale += "_incomplete"
    return FrameCorrectionPolicy(
        new_contract=True,
        scale_convention=scale,
        normalization_status=norm_status,
        normalization_divisor=normalization,
        normalization_unit=norm_unit,
        normalization_components=components,
        illumination_status=illum_status,
        illumination_divisor=illumination,
        illumination_convention=(
            "total_flux_H" if illumination is not None else None
        ),
        vertical_intercepted_fraction=vertical,
        horizontal_intercepted_fraction=horizontal,
        intercepted_fraction=intercepted,
        calibrated=calibrated,
    )


def corrected_curve_from_record(
    record,
    *,
    footprint_action=FOOTPRINT_KEEP,
    replacement_illumination=None,
    replacement_convention=None,
    allow_convention_change=False,
):
    """Rebuild a corrected scalar curve from its immutable base values.

    Every action starts from ``base_croibg`` and its variance, so applying or
    replacing a footprint repeatedly cannot compound scaling. Unknown or
    unavailable provenance is rejected rather than guessed.

    :returns: ``(curve, errors, illumination_status, convention)``.
    :rtype: tuple
    """
    if footprint_action not in (
        FOOTPRINT_KEEP,
        FOOTPRINT_APPLY,
        FOOTPRINT_REMOVE,
    ):
        raise ValueError(f"unknown footprint action {footprint_action!r}")
    if record.base_croibg is None or record.base_croibg_variance is None:
        raise ValueError("curve record has no reversible base curve")

    curve = np.asarray(record.base_croibg, dtype=np.float64).copy()
    variance = np.asarray(record.base_croibg_variance, dtype=np.float64).copy()
    if record.normalization_status == "applied":
        if record.normalization_divisor is None:
            raise ValueError("applied normalization has no stored divisor")
        divisor = np.asarray(record.normalization_divisor, dtype=np.float64)
        if np.any(divisor == 0) or not np.all(np.isfinite(divisor)):
            raise ValueError("normalization divisor must be finite and nonzero")
        curve = curve / divisor
        variance = variance / np.square(divisor)
    elif record.normalization_status != "not_applied":
        raise ValueError(
            "normalization provenance is insufficient for safe re-reduction"
        )

    stored_status = record.illumination_status
    stored_convention = record.illumination_convention
    if footprint_action == FOOTPRINT_KEEP:
        if stored_status == "applied":
            if record.illumination_divisor is None:
                raise ValueError("applied illumination has no stored divisor")
            illumination = np.asarray(
                record.illumination_divisor, dtype=np.float64
            )
            status = "applied"
            convention = stored_convention
        elif stored_status == "not_applied":
            illumination = None
            status = "not_applied"
            convention = None
        else:
            raise ValueError(
                "illumination provenance is insufficient for safe re-reduction"
            )
    elif footprint_action == FOOTPRINT_REMOVE:
        if stored_status != "applied" or record.illumination_divisor is None:
            raise ValueError("there is no known stored illumination to remove")
        illumination = None
        status = "removed"
        convention = stored_convention
    else:
        if stored_status not in ("applied", "not_applied"):
            raise ValueError(
                "unknown stored illumination cannot be safely replaced"
            )
        if replacement_illumination is None or replacement_convention is None:
            raise ValueError(
                "applying illumination requires its divisor and convention"
            )
        if (
            stored_status == "applied"
            and stored_convention != replacement_convention
            and not allow_convention_change
        ):
            raise ValueError(
                "changing illumination convention requires an explicit "
                "reconstruction from the stored base curve"
            )
        illumination = np.asarray(replacement_illumination, dtype=np.float64)
        status = "replaced" if stored_status == "applied" else "applied_here"
        convention = replacement_convention

    if illumination is not None:
        if np.any(illumination <= 0) or not np.all(np.isfinite(illumination)):
            raise ValueError("illumination divisor must be finite and positive")
        curve = curve / illumination
        variance = variance / np.square(illumination)
    return curve, np.sqrt(variance), status, convention


def pixel_correction_branches(
    base_intensity,
    base_errors,
    combined_factor,
    polarization_factor,
    polarization_arm_factor=1.0,
):
    r"""Build diagnostic and CTR branches from one scalar ROI signal.

    The diagnostic intensity retains the established combined per-pixel
    correction (solid angle times polarization). The CTR photon branch uses
    the independently accumulated polarization-only mean, so detector solid
    angle never has to be divided out through a separately estimated scalar.
    Both branches start from the same background-subtracted signal and apply
    the actual-arm polarization ratio exactly once.

    This remains the existing ROI-mean approximation: it does not claim that
    a geometric mean correction equals a signal-weighted per-pixel sum.

    :param base_intensity: Background-subtracted ROI signal before pixel
        corrections.
    :param base_errors: One-sigma errors paired with ``base_intensity``.
    :param combined_factor: Mean combined pixel correction over the valid
        center-ROI pixels.
    :param polarization_factor: Mean inverse-polarization correction over the
        same valid center-ROI pixels; one when polarization is disabled.
    :param polarization_arm_factor: Ratio moving the home-geometry
        polarization mean to the actual arm position; one when polarization
        is disabled.
    :returns: ``(intensity, errors, ctr_intensity, ctr_errors)``.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]
    """
    base_intensity = np.asarray(base_intensity, dtype=np.float64)
    base_errors = np.asarray(base_errors, dtype=np.float64)
    arm = np.asarray(polarization_arm_factor, dtype=np.float64)
    combined_scale = np.asarray(combined_factor, dtype=np.float64) * arm
    polarization_scale = np.asarray(polarization_factor, dtype=np.float64) * arm
    return (
        base_intensity * combined_scale,
        base_errors * combined_scale,
        base_intensity * polarization_scale,
        base_errors * polarization_scale,
    )


def monitor_counter_candidates(scan):
    """Counter names of ``scan`` usable as a monitor normalization.

    A counter qualifies when it is a one-dimensional array with one value per
    image, which is what a divisive normalization needs.

    :param scan: A loaded scan object.
    :returns: Sorted counter names.
    :rtype: list[str]
    """
    names = []
    candidates = set(getattr(scan, "auxillary_counters", ()) or ())
    try:
        length = len(scan)
    except TypeError:
        length = None
    for name in candidates:
        value = getattr(scan, name, None)
        if value is None:
            continue
        array = np.atleast_1d(np.asarray(value))
        if array.ndim != 1:
            continue
        if length is not None and array.size not in (1, length):
            continue
        names.append(name)
    return sorted(names)


def normalization_divisor(scan, normalize_exposure, monitor_corrections, size):
    """Exposure and monitor divisor for every image of a scan.

    Pulls the counters off the scan object and hands them to
    :func:`~orgui.datautils.xrayutils.corrections.normalization.normalization_divisor`,
    which owns the arithmetic and the validation. Mirrors the reciprocal-space
    reconstruction, which divides each frame by its exposure time and by every
    configured monitor counter.

    :param scan: Scan object providing ``exposure_time`` and the monitor
        counters by attribute.
    :param bool normalize_exposure: Divide by ``scan.exposure_time`` when the
        scan provides it.
    :param monitor_corrections: Iterable of monitor counter names.
    :param int size: Number of images, used to broadcast scalar counters.
    :returns: ``(divisor, applied)`` -- an array of shape ``(size,)`` and the
        names of the normalizations that contributed.
    :rtype: tuple
    :raises ValueError: If a counter is missing, or has a non-positive or
        non-finite value that would make the normalization undefined.
    """
    exposure = None
    if normalize_exposure:
        # A scan backend that reports no exposure time is not an error; the
        # reconstruction records the normalization as unavailable rather than
        # failing the job, and this follows it.
        exposure = getattr(scan, "exposure_time", None)

    monitors = {}
    for name in monitor_corrections:
        if not hasattr(scan, name):
            raise ValueError(f"Active scan has no monitor counter named {name!r}")
        monitors[name] = getattr(scan, name)

    return _divisor_from_counters(
        size, exposure_time=exposure, monitors=monitors
    )


def stationary_correction_factors(
    alpha,
    delta,
    gamma,
    use_lorentz=False,
    use_footprint=False,
    beam_profile=None,
    sample_size=None,
    normalization=None,
    solid_angle_mean=None,
    illumination_divisor=None,
):
    r"""Correction divisors for one stationary-scan trajectory.

    :param alpha: Incidence angle per image, in radian. This is the ``mu``
        circle, and the angle the footprint corrections depend on.
    :param delta: In-plane detector angle per image, in radian.
    :param gamma: Out-of-plane detector angle per image, in radian.
    :param bool use_lorentz: Add ``C_Lorentz`` -- the *stationary-mode*
        factor :math:`1/\sin\gamma`, not the rocking-scan one, as
        :func:`~orgui.datautils.xrayutils.corrections.measurement.mode_components`
        decides. Stationary integration has no rod-interception factor.
    :param bool use_footprint: Add the ``C_illum_area`` divisor and its
        diagnostic numerator ``C_flux_on_sample``.
    :param beam_profile: A
        :class:`~orgui.datautils.xrayutils.corrections.beamprofile.BeamProfile`,
        required when ``use_footprint`` is set.
    :param float sample_size: Sample size along the beam in meters, required
        when ``use_footprint`` is set.
    :param normalization: Optional per-image exposure and monitor divisor
        from :func:`normalization_divisor`, stored as ``C_norm``.
    :param solid_angle_mean: Optional per-image region mean of
        :math:`1/\widetilde{\Omega}` from
        :func:`~orgui.datautils.xrayutils.corrections.detector.roi_mean_inverse_solid_angle`,
        stored as ``C_solid_angle``. Pass it when the solid-angle correction
        was applied to the intensity, so that :func:`structure_factor` can
        divide it back out; a region sum is already a complete angular
        integral and must not carry it (finding F6).
    :param illumination_divisor: Optional explicit total-flux illumination
        :math:`H`. Mutually exclusive with the legacy ``use_footprint`` path.
    :returns: The factors, each broadcast to the shape of ``alpha``.
    :rtype: CorrectionFactors
    :raises ValueError: If the footprint correction is requested without a
        beam profile or a sample size.
    """
    alpha = np.asarray(alpha, dtype=np.float64)
    factors = {}
    applied = []

    if normalization is not None:
        factors["C_norm"] = np.broadcast_to(
            np.asarray(normalization, dtype=np.float64), alpha.shape
        ).copy()
        applied.append("normalization")

    if solid_angle_mean is not None:
        factors["C_solid_angle"] = np.broadcast_to(
            np.asarray(solid_angle_mean, dtype=np.float64), alpha.shape
        ).copy()
        applied.append("solid_angle")

    if illumination_divisor is not None:
        if use_footprint:
            raise ValueError(
                "explicit total-flux illumination and legacy footprint are "
                "mutually exclusive"
            )
        factors["C_illumination"] = np.broadcast_to(
            np.asarray(illumination_divisor, dtype=np.float64), alpha.shape
        ).copy()
        applied.append("total_flux_illumination")

    if use_footprint:
        if beam_profile is None:
            raise ValueError(
                "the footprint correction needs a beam profile; configure one "
                "in the integration corrections dialog"
            )
        if sample_size is None or not sample_size > 0:
            raise ValueError(
                f"the footprint correction needs a positive sample size, got "
                f"{sample_size!r}"
            )
        flux, area = beam_profile.corrections(alpha, sample_size)
        factors["C_flux_on_sample"] = np.broadcast_to(flux, alpha.shape).copy()
        factors["C_illum_area"] = np.broadcast_to(area, alpha.shape).copy()
        applied.append("footprint")

    if use_lorentz:
        components = measurement.mode_components(
            measurement.STATIONARY, alpha=alpha, delta=delta, gamma=gamma
        )
        for name, value in components.items():
            factors[name] = np.broadcast_to(value, alpha.shape).copy()
        applied.append("lorentz")

    return CorrectionFactors(factors, applied)


def apply_stationary_corrections(intensity, errors, factors):
    """Divide an intensity and its errors by the intensity-level factors.

    ``C_flux_on_sample`` is not an additional divisor. The numerical active
    area ``C_illum_area`` already contains that overlap integral. The Lorentz
    factor is applied separately by :func:`structure_factor`.

    :param intensity: Integrated intensity per image.
    :param errors: 1-sigma errors of ``intensity``.
    :param CorrectionFactors factors: Divisors to apply.
    :returns: ``(intensity, errors)`` corrected.
    :rtype: tuple of numpy.ndarray
    """
    divisor = factors.divisor(
        "C_norm", "C_illum_area", "C_illumination"
    )
    return np.asarray(intensity) / divisor, np.asarray(errors) / divisor


def structure_factor(intensity, errors, factors):
    r"""Form :math:`F^2_{hkl}` from an already corrected intensity.

    :math:`F^2 = I_\mathrm{corr} /
    (L_\mathrm{stationary}\,C_\mathrm{solid\,angle})`. Unlike a rocking
    scan, stationary area-detector integration has no rod-interception
    factor.

    ``C_solid_angle`` is present only when the solid-angle correction was
    applied to the intensity, and dividing by it removes that correction
    again. A region-summed intensity is already the complete angular
    integral, with every pixel weighted by the solid angle it subtends, so
    the correction double-counts the detector obliquity in a structure
    factor -- while remaining useful on the intensity itself for broad,
    non-rod features, where a differential cross-section is the goal. See
    ``doc/design/ctr_structure_factor_scale.md`` finding F6.

    :param intensity: Corrected intensity per image.
    :param errors: 1-sigma errors of ``intensity``.
    :param CorrectionFactors factors: Must contain ``C_Lorentz``; divides by
        ``C_solid_angle`` as well when it is present.
    :returns: ``(F2_hkl, F2_hkl_errors)``.
    :rtype: tuple of numpy.ndarray
    :raises KeyError: If the Lorentz factors are absent.
    """
    divisor = factors["C_Lorentz"] * factors.divisor("C_solid_angle")
    return np.asarray(intensity) / divisor, np.asarray(errors) / divisor


def structure_factor_from_policy(
    intensity,
    errors,
    factors,
    policy,
    *,
    wavelength=None,
    unitcell_area=None,
):
    """Form relative or calibrated stationary CTR :math:`|F|^2`.

    The calibrated branch is reachable only when the explicit total-flux
    policy applied both Q and H and the crystallographic scale inputs are
    supplied. Detector efficiency and external transmission are unity in this
    first wired implementation and must be recorded by the caller.
    """
    if policy.new_contract and not policy.ctr_scale_ready:
        raise ValueError(
            "new CTR normalization must apply both frame fluence Q and "
            "illumination H before a structure factor can be formed"
        )
    result, result_errors = structure_factor(intensity, errors, factors)
    if not policy.calibrated:
        return result, result_errors
    if wavelength is None or unitcell_area is None:
        raise ValueError(
            "calibrated total-flux F2 needs wavelength and surface "
            "unit-cell area"
        )
    prefactor = measurement.total_flux_prefactor(wavelength, unitcell_area)
    return result / prefactor, result_errors / prefactor

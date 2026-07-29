# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
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
"""Automatic UB seeding and reflection validation workflow."""

import logging

import numpy as np
from silx.gui import qt

from .. import logger_utils
from . import autoBraggSearch


logger = logging.getLogger(__name__)


class AutoBraggWorkflow:
    """Run automatic Bragg/UB workflows against an existing ``orGUI`` window.

    The workflow keeps orchestration separate from the main window while
    delegating access to scan, detector, UB, and reflection state back to the
    supplied ``orGUI`` instance.
    """

    def __init__(self, gui):
        """Create the workflow adapter.

        :param gui: Existing ``orGUI`` main-window instance.
        """
        self.gui = gui

    def __getattr__(self, name):
        """Delegate unknown attributes to the wrapped ``orGUI`` instance."""
        return getattr(self.gui, name)

    def autoFindBraggReference(self, **kwargs):
        """Search the active scan for a Bragg peak and seed the UB matrix.

        The search streams the max-pixel trace image by image, emits the first
        sharp post-burn-in peak candidate after a short lookahead confirmation,
        and validates that candidate immediately. If refinement or Bragg
        indexing fails, the stream continues from the next unread image until a
        valid seed is found or the scan ends.

        :param int burn_in:
            Number of initial images ignored for sharp-increase detection.
            Defaults to ``5``.
        :param float min_sharpness:
            Minimum robust level z-score for automatic candidates. Defaults to
            ``6``.
        :param float min_derivative_sharpness:
            Minimum robust first-difference z-score for automatic candidates.
            Defaults to ``6``.
        :param int lookahead:
            Number of following images read to confirm a local peak-like max
            trace candidate. Defaults to ``2``.
        :param int max_candidates:
            Maximum number of streamed candidates to validate before giving up.
            Defaults to no explicit limit.
        :param float qnorm_tolerance:
            Allowed Q-norm mismatch in Angstrom^-1 when matching Bragg
            reflections. Defaults to ``0.05``.
        :returns:
            Best UB seed, or ``None`` if no reliable seed was found.
        :rtype: orgui.app.autoBraggSearch.UBSeed or None

        .. note::
           CLI-capable when scan, mask, detector, and crystal state are already
           configured. The method updates the reference-reflection table and U.
        """
        if self.fscan is None:
            raise ValueError("No scan loaded")
        status_callback = kwargs.get("status_callback", None)

        def report(event, **fields):
            if status_callback is not None:
                status_callback(event, **fields)

        burn_in = kwargs.get("burn_in", 5)
        history = kwargs.get("history", 20)
        min_history = kwargs.get("min_history", 5)
        min_sharpness = kwargs.get("min_sharpness", 6.0)
        min_derivative_sharpness = kwargs.get("min_derivative_sharpness", 6.0)
        min_prominence_sharpness = kwargs.get("min_prominence_sharpness", 4.0)
        lookahead = kwargs.get("lookahead", 2)
        refractory = kwargs.get("refractory", 3)
        max_candidates = kwargs.get("max_candidates", None)
        mask_distance = kwargs.get("mask_distance", 3)
        qnorm_tolerance = kwargs.get("qnorm_tolerance", 0.05)
        adaptive_after_candidates = kwargs.get("adaptive_after_candidates", 5)
        adaptive_qnorm_tolerance = kwargs.get(
            "adaptive_qnorm_tolerance", max(0.25, 5.0 * qnorm_tolerance)
        )
        adaptive_assignment_pixel_tolerance = kwargs.get(
            "adaptive_assignment_pixel_tolerance",
            max(100.0, 3.0 * kwargs.get("assignment_pixel_tolerance", 30.0)),
        )
        adaptive_confirmation_pixel_tolerance = kwargs.get(
            "adaptive_confirmation_pixel_tolerance",
            max(
                30.0,
                3.0 * kwargs.get("confirmation_pixel_tolerance", 8.0),
            ),
        )
        adaptive_confirmation_image_tolerance = kwargs.get(
            "adaptive_confirmation_image_tolerance",
            max(8, 2 * kwargs.get("confirmation_image_tolerance", 3)),
        )
        adaptive_scale_detector_filter = kwargs.get(
            "adaptive_scale_detector_filter", True
        )
        adaptive_scale_detector_fraction = kwargs.get(
            "adaptive_scale_detector_fraction", 0.25
        )
        adaptive_scale_outlier_q_tolerance = kwargs.get(
            "adaptive_scale_outlier_q_tolerance",
            adaptive_qnorm_tolerance,
        )
        adaptive_scale_outlier_angle_fraction = kwargs.get(
            "adaptive_scale_outlier_angle_fraction",
            adaptive_scale_detector_fraction,
        )
        hkl_candidate_mode = kwargs.get("hkl_candidate_mode", "detector_position")
        if hkl_candidate_mode == "qnorm":
            logger.warning(
                "Automatic Bragg HKL assignment uses Q-norm search. This "
                "checks more hypotheses and may take longer than detector "
                "position matching."
            )
        elif hkl_candidate_mode != "detector_position":
            raise ValueError(
                "hkl_candidate_mode must be 'detector_position' or 'qnorm'"
            )
        axis_half_width = kwargs.get("axis_half_width", 1.0)
        fine_axis_half_width = kwargs.get("fine_axis_half_width", 0.4)
        roi_size = kwargs.get("roi_size", (80, 80))
        fine_roi_size = kwargs.get("fine_roi_size", (40, 40))

        mask = self.get_detector_mask()
        if mask is None:
            if logger_utils.get_logging_context() == "gui":
                # GUI-only: user-triggered automatic Bragg search confirmation.
                btn = qt.QMessageBox.question(
                    self.gui,
                    "Automatic Bragg search without mask",
                    "No detector mask is set. Automatic Bragg search will "
                    "continue without masked-pixel rejection or mask-aware "
                    "peak refinement.\n\nDo you want to continue?",
                    qt.QMessageBox.Yes | qt.QMessageBox.No,
                    qt.QMessageBox.No,
                )
                if btn != qt.QMessageBox.Yes:
                    report(
                        "failed",
                        message=(
                            "Automatic Bragg search cancelled because no "
                            "detector mask is set."
                        ),
                    )
                    return None
            logger.warning(
                "Automatic Bragg search is running without a detector mask. "
                "Masked-pixel rejection and mask-aware peak refinement will "
                "be disabled."
            )
            report(
                "warning",
                message=(
                    "No detector mask is set; masked-pixel rejection and "
                    "mask-aware peak refinement are disabled."
                ),
            )

        excluded_images = self.excludedImagesDialog.getData()
        self.ubcalc.crystal.setEnergy(self.ubcalc.ubCal.getEnergy() * 1e3)
        max_q = kwargs.get(
            "max_q",
            getattr(self.ubcalc.detectorCal, "Qmax", None),
        )

        logger.info("Start online automatic Bragg max-pixel search")
        report(
            "start",
            message=(
                f"Started automatic Bragg search: mode={hkl_candidate_mode}, images={len(self.fscan)}"  # noqa: E501
            ),
            images_read=0,
        )
        candidate_count = 0
        refined_peaks = []
        unmatched_peaks = []
        qnorm_scale = 1.0
        adaptive_matching = False
        last_scale_fit_count = 0
        scale_revalidation_done = False

        def scale_fit_detector_filter(peak, hkls, norms):
            """Filter scale-fit HKLs by broad predicted detector position."""
            if (
                not adaptive_scale_detector_filter
                or hkl_candidate_mode != "detector_position"
            ):
                return np.ones(len(hkls), dtype=bool)
            if len(hkls) == 0:
                return np.zeros(0, dtype=bool)
            try:
                refldict = self.ubcalc.calcReflection(np.asarray(hkls, dtype=float))
            except Exception:
                logger.debug(
                    "Cannot calculate vectorized detector-position filter "
                    "for adaptive Q-scale fit.",
                    exc_info=True,
                )
                return np.zeros(len(hkls), dtype=bool)
            detvsize, dethsize = self.ubcalc.detectorCal.detector.shape
            half_width = adaptive_scale_detector_fraction * np.array(
                [dethsize, detvsize], dtype=float
            )
            observed_xy = np.asarray(peak.xy, dtype=float)
            keep = np.zeros(len(hkls), dtype=bool)
            for key in ("xy_1", "xy_2"):
                xy = np.atleast_2d(np.asarray(refldict[key], dtype=float))
                finite = np.all(np.isfinite(xy), axis=1)
                in_region = np.all(
                    np.abs(xy - observed_xy) <= half_width,
                    axis=1,
                )
                keep |= finite & in_region
            return keep

        def filter_scale_fit_matches(qscale_fit):
            """Return provisional scale-fit matches without obvious outliers."""
            if qscale_fit is None:
                return []
            matches = list(qscale_fit.get("matches", ()))
            if not matches:
                return []
            scale = qscale_fit["qnorm_scale"]
            q_residual = np.array(
                [
                    abs(match["configured_norm"] * scale - match["qnorm"])
                    for match in matches
                ]
            )
            q_keep = q_residual <= adaptive_scale_outlier_q_tolerance
            if hkl_candidate_mode != "detector_position":
                for match, residual in zip(matches, q_residual):
                    match["scaled_q_residual"] = float(residual)
                return [match for match, keep in zip(matches, q_keep) if keep]

            hkls = np.asarray([match["hkl"] for match in matches], dtype=float)
            try:
                refldict = self.ubcalc.calcReflection(hkls)
                predicted_1 = np.atleast_2d(
                    np.asarray(refldict["angles_1"], dtype=float)
                )[:, [1, 2]]
                predicted_2 = np.atleast_2d(
                    np.asarray(refldict["angles_2"], dtype=float)
                )[:, [1, 2]]
            except Exception:
                logger.debug(
                    "Cannot calculate vectorized delta/gamma outlier check "
                    "for adaptive Q-scale fit.",
                    exc_info=True,
                )
                predicted_1 = np.full((len(matches), 2), np.inf)
                predicted_2 = np.full((len(matches), 2), np.inf)

            peaks = [match["peak"] for match in matches]
            mus = np.array([self.getMuOm(peak.imageno)[0] for peak in peaks])
            measured_gamma, measured_delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([peak.xy[1] for peak in peaks]),
                np.array([peak.xy[0] for peak in peaks]),
                mus,
            )
            measured = np.column_stack((measured_delta, measured_gamma))
            try:
                gamrange, delrange = self.ubcalc.detectorCal.rangegamdel(
                    float(np.nanmedian(mus))
                )
                angle_half_width = adaptive_scale_outlier_angle_fraction * (
                    np.array(
                        [
                            abs(delrange[1] - delrange[0]),
                            abs(gamrange[1] - gamrange[0]),
                        ],
                        dtype=float,
                    )
                )
            except Exception:
                logger.debug(
                    "Cannot calculate detector delta/gamma range for "
                    "adaptive Q-scale outlier check.",
                    exc_info=True,
                )
                angle_half_width = np.full(2, np.inf, dtype=float)

            residual_1 = np.abs(predicted_1 - measured)
            residual_2 = np.abs(predicted_2 - measured)
            angle_keep_1 = np.all(residual_1 <= angle_half_width, axis=1)
            angle_keep_2 = np.all(residual_2 <= angle_half_width, axis=1)
            delgam_error = np.minimum(
                np.linalg.norm(residual_1, axis=1),
                np.linalg.norm(residual_2, axis=1),
            )
            keep = q_keep & (angle_keep_1 | angle_keep_2)
            for match, qres, dgerr in zip(matches, q_residual, delgam_error):
                match["scaled_q_residual"] = float(qres)
                match["delgam_error"] = float(dgerr)
            return [match for match, keep_match in zip(matches, keep) if keep_match]

        def validate_automatic_peak(
            peak,
            current_qnorm_tolerance,
            current_assignment_pixel_tolerance,
            current_confirmation_pixel_tolerance,
            current_confirmation_image_tolerance,
        ):
            """Try HKL assignment and second-peak confirmation for ``peak``."""
            if hkl_candidate_mode == "qnorm":
                seeds = autoBraggSearch.rank_ub_seeds(
                    [peak],
                    self.ubcalc.detectorCal,
                    self.getMuOm,
                    self.ubcalc.ubCal,
                    self.ubcalc.crystal,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                    qnorm_tolerance=current_qnorm_tolerance,
                    max_q=max_q,
                    qnorm_scale=qnorm_scale,
                )
            else:
                seeds = self._rankAutomaticSeedsByDetectorPosition(
                    peak,
                    q_tolerance=current_qnorm_tolerance,
                    max_q=max_q,
                    pixel_tolerance=current_assignment_pixel_tolerance,
                    max_reflections=kwargs.get("assignment_reflections", 40),
                    qnorm_scale=qnorm_scale,
                )
            if not seeds:
                return None, None, "no HKL hypotheses in local Q region"
            report(
                "hkl_hypotheses",
                count=len(seeds),
                message=(
                    f"Testing {len(seeds)} HKL hypotheses for candidate image {peak.imageno}."  # noqa: E501
                ),
            )

            if hkl_candidate_mode == "detector_position":
                seed_hypotheses = seeds
            else:
                seed_hypotheses = seeds[: kwargs.get("max_seed_hypotheses", 12)]
            for candidate_seed in seed_hypotheses:
                confirmation = self._confirmAutomaticBraggSeed(
                    candidate_seed,
                    mask=mask,
                    excluded_images=excluded_images,
                    max_q=max_q,
                    axis_half_width=axis_half_width,
                    fine_axis_half_width=fine_axis_half_width,
                    roi_size=roi_size,
                    fine_roi_size=fine_roi_size,
                    pixel_tolerance=current_confirmation_pixel_tolerance,
                    image_tolerance=current_confirmation_image_tolerance,
                    max_reflections=kwargs.get("confirmation_reflections", 12),
                    intensity_ratio_check=kwargs.get("intensity_ratio_check", True),
                    prominence_threshold=kwargs.get(
                        "confirmation_prominence_threshold", 6.0
                    ),
                    status_callback=status_callback,
                )
                if confirmation is not None:
                    if np.isfinite(confirmation["predicted_intensity_ratio"]):
                        intensity_message = (
                            "intensity ratio obs/calc={:.3g}/{:.3g}".format(
                                confirmation["observed_intensity_ratio"],
                                confirmation["predicted_intensity_ratio"],
                            )
                        )
                    else:
                        intensity_message = (
                            "prominence z seed/confirmation={:.3g}/{:.3g}".format(
                                confirmation["seed_prominence_z"],
                                confirmation["confirmation_prominence_z"],
                            )
                        )
                    report(
                        "confirmation",
                        message=(
                            "Confirmed hkl={} with second hkl={}; "
                            "pixel error={:.3g}, image error={}, {}".format(
                                candidate_seed.hkl,
                                confirmation["hkl"],
                                confirmation["pixel_error"],
                                confirmation["image_error"],
                                intensity_message,
                            )
                        ),
                    )
                    return candidate_seed, confirmation, None
            return None, None, "no second Bragg peak confirmed the HKL assignment"

        def accept_automatic_seed(seed, confirmation):
            """Apply a confirmed automatic UB seed and add reference peaks."""
            self.ubcalc.ubCal.setU(seed.U)
            self.ubcalc.uedit.setU(seed.U)
            eventdict = {"x": seed.peak.xy[0], "y": seed.peak.xy[1]}
            refl = self.reflectionSel.addReflection(
                eventdict, seed.peak.imageno, seed.hkl
            )
            if confirmation is not None:
                confirm_eventdict = {
                    "x": confirmation["peak"].xy[0],
                    "y": confirmation["peak"].xy[1],
                }
                self.reflectionSel.addReflection(
                    confirm_eventdict,
                    confirmation["peak"].imageno,
                    confirmation["hkl"],
                )
            extra_q_tolerance = (
                adaptive_qnorm_tolerance if adaptive_matching else qnorm_tolerance
            )
            extra_count = self._addAutomaticObservedReflections(
                refined_peaks,
                seed,
                confirmation,
                q_tolerance=extra_q_tolerance,
                score_tolerance_deg=2.0,
                intensity_ratio_check=kwargs.get("intensity_ratio_check", True),
                status_callback=status_callback,
            )
            if extra_count:
                report(
                    "accepted",
                    message=(
                        f"Added {extra_count} additional observed Bragg reference "
                        "reflection(s) after U confirmation."
                    ),
                )
            self.reflectionSel.setReflectionActive(refl.identifier)
            self._onChangeImage(seed.peak.imageno)
            self._onCenterGraph(seed.peak.xy)
            self.ubcalc.updateReflectionMismatch()
            self.ubcalc.sigReplotRequest.emit(False)
            logger.info(
                "Automatic Bragg seed selected hkl=%s at image %s, xy=%s, "
                "score=%.4g, norm mismatch=%.4g Angstrom^-1.",
                seed.hkl,
                seed.peak.imageno,
                seed.peak.xy,
                seed.score,
                seed.norm_mismatch,
            )
            report(
                "success",
                message=(
                    f"Selected hkl={seed.hkl} at image={seed.peak.imageno}, xy=({seed.peak.xy[0]:.2f}, {seed.peak.xy[1]:.2f}), "  # noqa: E501
                    f"norm mismatch={seed.norm_mismatch:.4g}"
                ),
            )
            return seed

        candidates = autoBraggSearch.iter_sharp_peak_candidates(
            self.fscan,
            mask=mask,
            background_image=self.background_image,
            excluded_images=excluded_images,
            burn_in=burn_in,
            history=history,
            min_history=min_history,
            level_z=min_sharpness,
            derivative_z=min_derivative_sharpness,
            lookahead=lookahead,
            min_prominence_z=min_prominence_sharpness,
            refractory=refractory,
            mask_distance=mask_distance,
            max_workers=self.numberthreads,
        )
        for maximum in candidates:
            report(
                "candidate",
                message=(
                    f"Candidate {candidate_count + 1}: image={maximum.imageno}, xy=({maximum.xy[0]:.2f}, {maximum.xy[1]:.2f}), "  # noqa: E501
                    f"max={maximum.value:.4g}, z={maximum.sharpness:.3g}, dz={maximum.derivative_sharpness:.3g}"  # noqa: E501
                ),
                images_read=min(
                    len(self.fscan),
                    maximum.imageno + lookahead + 1,
                ),
            )
            candidate_count += 1
            if max_candidates is not None and candidate_count > max_candidates:
                logger.warning(
                    "Automatic Bragg search stopped after %s rejected candidates.",
                    max_candidates,
                )
                report(
                    "failed",
                    message=(f"Stopped after {max_candidates} rejected candidates."),
                )
                return None
            try:
                peak = autoBraggSearch.refine_peak_3d(
                    self.fscan,
                    maximum,
                    axis_half_width=axis_half_width,
                    roi_size=roi_size,
                    fine_axis_half_width=fine_axis_half_width,
                    fine_roi_size=fine_roi_size,
                    mask=mask,
                    background_image=self.background_image,
                    excluded_images=excluded_images,
                    max_workers=self.numberthreads,
                )
            except Exception:
                logger.warning(
                    "Skipping automatic Bragg candidate in image %s after "
                    "3D peak refinement failed.",
                    maximum.imageno,
                    exc_info=True,
                )
                report(
                    "rejected",
                    message=(
                        f"Rejected candidate image {maximum.imageno}: 3D refinement failed."  # noqa: E501
                    ),
                )
                continue
            if not autoBraggSearch.far_from_mask(mask, peak.xy, mask_distance):
                logger.info(
                    "Skipping automatic Bragg candidate in image %s because "
                    "the refined peak is too close to masked pixels.",
                    peak.imageno,
                )
                report(
                    "rejected",
                    message=(
                        f"Rejected candidate image {peak.imageno}: refined peak is too "
                        "close to masked pixels."
                    ),
                )
                continue
            report(
                "refined",
                message=(
                    f"Refined candidate: image={peak.imageno}, axis={peak.axis_value:.5g}, "  # noqa: E501
                    f"xy=({peak.xy[0]:.2f}, {peak.xy[1]:.2f})"
                ),
            )
            refined_peaks.append(peak)
            if (
                not scale_revalidation_done
                and len(unmatched_peaks) >= adaptive_after_candidates
                and len(unmatched_peaks) > last_scale_fit_count
            ):
                broadening_started = not adaptive_matching
                qscale_fit = autoBraggSearch.estimate_qnorm_scale(
                    unmatched_peaks,
                    self.ubcalc.detectorCal,
                    self.getMuOm,
                    self.ubcalc.ubCal,
                    self.ubcalc.crystal,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                    adaptive_qnorm_tolerance,
                    max_q=max_q,
                    hkl_filter=scale_fit_detector_filter,
                )
                last_scale_fit_count = len(unmatched_peaks)
                adaptive_matching = True
                if qscale_fit is not None:
                    qnorm_scale = qscale_fit["qnorm_scale"]
                    if broadening_started:
                        logger.warning(
                            "Automatic Bragg search broadened the Q shell "
                            "from %.4g to %.4g Angstrom^-1 after %s "
                            "unmatched candidates. Tentative reciprocal Q "
                            "scale is %.5g (direct lattice scale %.5g, "
                            "scatter %.3g).",
                            qnorm_tolerance,
                            adaptive_qnorm_tolerance,
                            len(unmatched_peaks),
                            qnorm_scale,
                            qscale_fit["direct_lattice_scale"],
                            qscale_fit["scatter"],
                        )
                        report(
                            "adaptive",
                            message=(
                                "Q shell broadened from {:.4g} to {:.4g} "
                                "Angstrom^-1 after {} unmatched candidates. "
                                "Tentative reciprocal Q scale={:.5g} "
                                "(direct lattice scale={:.5g}, scatter={:.3g}).".format(
                                    qnorm_tolerance,
                                    adaptive_qnorm_tolerance,
                                    len(unmatched_peaks),
                                    qnorm_scale,
                                    qscale_fit["direct_lattice_scale"],
                                    qscale_fit["scatter"],
                                )
                            ),
                        )
                    else:
                        report(
                            "adaptive",
                            message=(
                                "Recomputed adaptive Q scale from {} "
                                "unmatched candidates: reciprocal scale={:.5g} "
                                "(direct lattice scale={:.5g}, scatter={:.3g}).".format(
                                    len(unmatched_peaks),
                                    qnorm_scale,
                                    qscale_fit["direct_lattice_scale"],
                                    qscale_fit["scatter"],
                                )
                            ),
                        )
                    inlier_matches = filter_scale_fit_matches(qscale_fit)
                    if len(inlier_matches) >= adaptive_after_candidates:
                        qnorm_scale = float(
                            np.median(
                                [match["qnorm_scale"] for match in inlier_matches]
                            )
                        )
                        report(
                            "adaptive",
                            message=(
                                "Scale-fit outlier check kept {}/{} "
                                "provisional reflections; rechecking their "
                                "automatic UB validation with inlier Q "
                                "scale {:.5g}.".format(
                                    len(inlier_matches),
                                    len(qscale_fit["matches"]),
                                    qnorm_scale,
                                )
                            ),
                        )
                        for match in inlier_matches:
                            seed, confirmation, _reason = validate_automatic_peak(
                                match["peak"],
                                adaptive_qnorm_tolerance,
                                adaptive_assignment_pixel_tolerance,
                                adaptive_confirmation_pixel_tolerance,
                                adaptive_confirmation_image_tolerance,
                            )
                            if seed is not None:
                                return accept_automatic_seed(seed, confirmation)
                        report(
                            "adaptive",
                            message=(
                                f"Rechecked {len(inlier_matches)} scale-fit provisional "  # noqa: E501
                                "reflections; no second-peak confirmation "
                                "was found."
                            ),
                        )
                        scale_revalidation_done = True
                    else:
                        report(
                            "adaptive",
                            message=(
                                "Scale-fit outlier check kept {}/{} "
                                "provisional reflections; continuing the "
                                "stream until at least {} remain.".format(
                                    len(inlier_matches),
                                    len(qscale_fit["matches"]),
                                    adaptive_after_candidates,
                                )
                            ),
                        )
                else:
                    if broadening_started:
                        logger.warning(
                            "Automatic Bragg search broadened the Q shell "
                            "from %.4g to %.4g Angstrom^-1 after %s "
                            "unmatched candidates. No reliable Q-scale fit "
                            "was available.",
                            qnorm_tolerance,
                            adaptive_qnorm_tolerance,
                            len(unmatched_peaks),
                        )
                        report(
                            "adaptive",
                            message=(
                                f"Q shell broadened from {qnorm_tolerance:.4g} to {adaptive_qnorm_tolerance:.4g} "  # noqa: E501
                                f"Angstrom^-1 after {len(unmatched_peaks)} unmatched candidates. "  # noqa: E501
                                "No reliable Q-scale fit was available."
                            ),
                        )
                    else:
                        report(
                            "adaptive",
                            message=(
                                f"Rechecked adaptive Q-scale fit from {len(unmatched_peaks)} "  # noqa: E501
                                "unmatched candidates; no reliable Q-scale "
                                "fit was available."
                            ),
                        )

            current_qnorm_tolerance = (
                adaptive_qnorm_tolerance if adaptive_matching else qnorm_tolerance
            )
            current_assignment_pixel_tolerance = (
                adaptive_assignment_pixel_tolerance
                if adaptive_matching
                else kwargs.get("assignment_pixel_tolerance", 30.0)
            )
            current_confirmation_pixel_tolerance = (
                adaptive_confirmation_pixel_tolerance
                if adaptive_matching
                else kwargs.get("confirmation_pixel_tolerance", 8.0)
            )
            current_confirmation_image_tolerance = (
                adaptive_confirmation_image_tolerance
                if adaptive_matching
                else kwargs.get("confirmation_image_tolerance", 3)
            )

            seed, confirmation, rejection_reason = validate_automatic_peak(
                peak,
                current_qnorm_tolerance,
                current_assignment_pixel_tolerance,
                current_confirmation_pixel_tolerance,
                current_confirmation_image_tolerance,
            )
            if seed is None:
                logger.warning(
                    "Skipping automatic Bragg candidate in image %s because %s.",
                    peak.imageno,
                    rejection_reason,
                )
                report(
                    "rejected",
                    message=(
                        f"Rejected candidate image {peak.imageno}: {rejection_reason}."
                    ),
                )
                unmatched_peaks.append(peak)
                continue
            return accept_automatic_seed(seed, confirmation)

        logger.warning(
            "Automatic Bragg search reached the end of the scan without a "
            "valid Bragg seed."
        )
        report(
            "failed",
            message="Reached the end of the scan without a valid Bragg seed.",
            images_read=len(self.fscan),
        )
        return None

    def _rankAutomaticSeedsByDetectorPosition(
        self,
        peak,
        q_tolerance=0.05,
        max_q=None,
        pixel_tolerance=30.0,
        max_reflections=40,
        qnorm_scale=1.0,
    ):
        """Rank HKL hypotheses by predicted detector position.

        :param orgui.app.autoBraggSearch.RefinedPeak peak:
            Refined first peak.
        :param float q_tolerance:
            Q-shell half-width around the measured peak norm in Angstrom^-1.
        :param float max_q:
            Optional maximum allowed reflection norm in Angstrom^-1.
        :param float pixel_tolerance:
            Maximum detector-position mismatch in pixels for candidate seeds.
        :param int max_reflections:
            Maximum intensity-sorted HKLs from the Q shell to test.
        :param float qnorm_scale:
            Multiplicative scale applied to predicted reciprocal-vector norms
            before Q-shell matching.
        :returns:
            Candidate seeds sorted by detector mismatch, then intensity rank.
        :rtype: list[orgui.app.autoBraggSearch.UBSeed]
        """
        q_phi = autoBraggSearch.q_phi_from_peak(
            peak,
            self.ubcalc.detectorCal,
            self.getMuOm,
            self.ubcalc.chi,
            self.ubcalc.phi,
            self.ubcalc.ubCal.getK(),
        )
        qnorm = float(np.linalg.norm(q_phi))
        hkls, intensity, norms = autoBraggSearch.allowed_bragg_in_q_region(
            self.ubcalc.crystal,
            qnorm,
            q_tolerance,
            max_q=max_q,
            qnorm_scale=qnorm_scale,
        )
        if len(hkls) > 0:
            expanded = []
            seen = set()
            for hkl, inten, norm in zip(
                hkls[:max_reflections],
                intensity[:max_reflections],
                norms[:max_reflections],
            ):
                same_q_tolerance = max(1e-6, 1e-5 * max(1.0, abs(norm)))
                group_hkls, group_intensity, group_norms = (
                    autoBraggSearch.allowed_bragg_same_q(
                        self.ubcalc.crystal,
                        float(norm),
                        tolerance=same_q_tolerance,
                        max_q=max_q,
                    )
                )
                if len(group_hkls) == 0:
                    group_hkls = np.asarray([hkl], dtype=float)
                    group_intensity = np.asarray([inten], dtype=float)
                    group_norms = np.asarray([norm], dtype=float)
                for group_hkl, group_inten, group_norm in zip(
                    group_hkls, group_intensity, group_norms
                ):
                    key = tuple(np.round(group_hkl, 8))
                    if key in seen:
                        continue
                    seen.add(key)
                    expanded.append((group_hkl, group_inten, group_norm))
            if expanded:
                hkls = np.asarray([item[0] for item in expanded], dtype=float)
                intensity = np.asarray([item[1] for item in expanded], dtype=float)
                norms = np.asarray([item[2] for item in expanded], dtype=float)
                order = np.argsort(intensity)[::-1]
                hkls = hkls[order]
                intensity = intensity[order]
                norms = norms[order]
        previous_u = np.asarray(self.ubcalc.ubCal.getU()).copy()
        seeds = []
        try:
            for hkl, _intensity, norm in zip(
                hkls,
                intensity,
                norms,
            ):
                try:
                    refldict = self.searchPixelCoordHKL(hkl)
                except Exception:
                    logger.debug(
                        "Cannot predict first-peak HKL hypothesis %s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                solution_mismatches = []
                for solution in (1, 2):
                    if not refldict.get(f"selectable_{solution}", False):
                        continue
                    if f"imageno_{solution}" not in refldict:
                        continue
                    predicted_xy = np.asarray(refldict[f"xy_{solution}"], dtype=float)
                    predicted_imageno = refldict[f"imageno_{solution}"]
                    pixel_error = np.linalg.norm(peak.xy - predicted_xy)
                    image_error = abs(peak.imageno - int(predicted_imageno))
                    solution_mismatches.append((pixel_error + image_error, pixel_error))
                if not solution_mismatches:
                    continue
                score, pixel_error = min(solution_mismatches)
                if pixel_error > pixel_tolerance:
                    continue
                mu, omega = self.getMuOm(peak.imageno)
                gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                    np.array([peak.xy[1]]),
                    np.array([peak.xy[0]]),
                    mu,
                )
                pos = np.array(
                    [
                        mu,
                        float(delta[0]),
                        float(gamma[0]),
                        omega,
                        self.ubcalc.chi,
                        self.ubcalc.phi,
                    ]
                )
                try:
                    U = autoBraggSearch.seed_u_from_single_reflection(
                        self.ubcalc.ubCal, pos, hkl
                    )
                except Exception:
                    logger.debug(
                        "Cannot calculate one-reflection UB seed for hkl=%s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                seeds.append(
                    autoBraggSearch.UBSeed(
                        U=U,
                        hkl=np.asarray(hkl, dtype=float),
                        peak=peak,
                        score=float(score),
                        norm_mismatch=abs(float(norm) * qnorm_scale - qnorm),
                        position_mismatch=float(pixel_error),
                    )
                )
        finally:
            self.ubcalc.ubCal.setU(previous_u)
        return sorted(
            seeds,
            key=lambda seed: (
                seed.position_mismatch,
                seed.norm_mismatch,
                seed.score,
            ),
        )

    def _automaticReflectionIntensity(self, hkl):
        """Return calculated ``abs(F)**2`` for one Bragg reflection."""
        hkl = np.asarray(hkl, dtype=float)
        h = np.asarray([hkl[0]], dtype=float)
        k = np.asarray([hkl[1]], dtype=float)
        l = np.asarray([hkl[2]], dtype=float)  # noqa: E741
        if hasattr(self.ubcalc.crystal, "uc_bulk"):
            structure_factor = self.ubcalc.crystal.uc_bulk.F_uc(h, k, l)[0]
        else:
            structure_factor = self.ubcalc.crystal.F_uc(h, k, l)[0]
        return float(np.abs(structure_factor) ** 2)

    def _addAutomaticObservedReflections(
        self,
        peaks,
        seed,
        confirmation,
        q_tolerance,
        score_tolerance_deg=2.0,
        intensity_ratio_check=True,
        intensity_ratio_tolerance=0.5,
        status_callback=None,
    ):
        """Add already-observed peaks that match the confirmed UB matrix.

        :param list[orgui.app.autoBraggSearch.RefinedPeak] peaks:
            Refined peaks already observed during the automatic search.
        :param orgui.app.autoBraggSearch.UBSeed seed:
            Accepted seed reflection.
        :param dict confirmation:
            Accepted confirmation reflection dictionary.
        :param float q_tolerance:
            Maximum Q-norm mismatch in Angstrom^-1.
        :param float score_tolerance_deg:
            Maximum angular UB mismatch in deg.
        :param bool intensity_ratio_check:
            If true, require each added reflection intensity to match the mean
            observed/calculated intensity scale within tolerance.
        :param float intensity_ratio_tolerance:
            Maximum fractional deviation from the mean intensity scale.
        :param callable status_callback:
            Optional automatic-search status callback for GUI/CLI feedback.
        :returns:
            Number of additional reflections added.
        :rtype: int
        """
        if not peaks:
            return 0
        accepted = []
        for accepted_peak in (seed.peak, confirmation.get("peak")):
            if accepted_peak is not None:
                accepted.append(accepted_peak)
        intensity_scale = []
        if intensity_ratio_check:
            intensity_refs = [
                (seed.peak, seed.hkl),
                (confirmation.get("peak"), confirmation.get("hkl")),
            ]
            for ref_peak, ref_hkl in intensity_refs:
                if ref_peak is None or ref_hkl is None:
                    continue
                ref_intensity = autoBraggSearch.estimate_rocking_intensity(ref_peak)
                ref_f2 = self._automaticReflectionIntensity(ref_hkl)
                if (
                    ref_intensity is not None
                    and ref_intensity["intensity"] > 0.0
                    and ref_f2 > 0.0
                ):
                    intensity_scale.append(ref_intensity["intensity"] / ref_f2)

        def is_same_peak(peak_a, peak_b):
            return (
                peak_a.imageno == peak_b.imageno
                and np.linalg.norm(peak_a.xy - peak_b.xy) < 2.0
            )

        added = 0
        tested_peaks = []
        for peak in peaks:
            if any(is_same_peak(peak, existing) for existing in accepted):
                continue
            if any(is_same_peak(peak, existing) for existing in tested_peaks):
                continue
            tested_peaks.append(peak)
            mu, omega = self.getMuOm(peak.imageno)
            gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
                np.array([peak.xy[1]]),
                np.array([peak.xy[0]]),
                mu,
            )
            angles = np.array(
                [
                    mu,
                    float(delta[0]),
                    float(gamma[0]),
                    omega,
                    self.ubcalc.chi,
                    self.ubcalc.phi,
                ]
            )
            try:
                h, k, l = self.ubcalc.angles.anglesToHkl(*angles)  # noqa: E741
            except Exception:
                logger.debug(
                    "Cannot convert observed automatic peak at image %s, xy=%s to HKL.",
                    peak.imageno,
                    peak.xy,
                    exc_info=True,
                )
                continue
            measured_hkl = np.asarray([h, k, l], dtype=float).reshape(3)
            rounded_hkl = np.rint(measured_hkl).astype(float)
            if not np.all(np.isfinite(measured_hkl)) or np.allclose(rounded_hkl, 0.0):
                continue

            if status_callback is not None:
                status_callback(
                    "observed_hkl",
                    message=(
                        f"Observed peak image={peak.imageno}, xy=({peak.xy[0]:.2f}, {peak.xy[1]:.2f}): "  # noqa: E501
                        f"HKL from angles=({measured_hkl[0]:.3f}, {measured_hkl[1]:.3f}, {measured_hkl[2]:.3f}), rounded={rounded_hkl}."  # noqa: E501
                    ),
                )
            logger.info(
                "Automatic observed peak image=%s, xy=%s gives HKL=%s; rounded to %s.",
                peak.imageno,
                peak.xy,
                measured_hkl,
                rounded_hkl,
            )
            try:
                mismatch = self.ubcalc.getReflectionMismatch(
                    np.atleast_2d(rounded_hkl),
                    np.atleast_2d(angles),
                )
            except Exception:
                logger.debug(
                    "Cannot score observed automatic peak hkl=%s",
                    rounded_hkl,
                    exc_info=True,
                )
                continue
            angular_mismatch_deg = float(
                np.rad2deg(np.asarray(mismatch["angle_mismatch"])[0])
            )
            q_mismatch = float(np.asarray(mismatch["norm_mismatch"])[0])
            message = (
                f"Observed peak image={peak.imageno} rounded hkl={rounded_hkl}: angular "  # noqa: E501
                f"mismatch={angular_mismatch_deg:.4g} deg, Q mismatch={q_mismatch:.4g} Angstrom^-1 "  # noqa: E501
                f"(limits {score_tolerance_deg:.4g} deg, {q_tolerance:.4g} Angstrom^-1)."  # noqa: E501
            )
            logger.info(message)
            if status_callback is not None:
                status_callback("observed_hkl_test", message=message)
            if angular_mismatch_deg > score_tolerance_deg or q_mismatch > q_tolerance:
                continue
            intensity_message = ""
            peak_intensity = None
            if intensity_ratio_check:
                peak_intensity = autoBraggSearch.estimate_rocking_intensity(peak)
                peak_f2 = self._automaticReflectionIntensity(rounded_hkl)
                if (
                    peak_intensity is None
                    or peak_intensity["intensity"] <= 0.0
                    or peak_f2 <= 0.0
                    or not intensity_scale
                ):
                    logger.info(
                        "Rejected observed automatic Bragg reflection hkl=%s "
                        "at image %s because intensity scale could not be "
                        "tested.",
                        rounded_hkl,
                        peak.imageno,
                    )
                    continue
                observed_scale = peak_intensity["intensity"] / peak_f2
                mean_scale = float(np.mean(intensity_scale))
                scale_error = abs(observed_scale - mean_scale) / mean_scale
                intensity_message = f", intensity scale={observed_scale:.4g}, mean scale={mean_scale:.4g}, relerr={scale_error:.3g}"  # noqa: E501
                message = (
                    f"Observed peak image={peak.imageno} rounded hkl={rounded_hkl} intensity check: "  # noqa: E501
                    f"scale={observed_scale:.4g}, mean scale={mean_scale:.4g}, relerr={scale_error:.3g} "  # noqa: E501
                    f"(limit {intensity_ratio_tolerance:.3g})."
                )
                logger.info(message)
                if status_callback is not None:
                    status_callback("observed_intensity_test", message=message)
                if scale_error > intensity_ratio_tolerance:
                    continue
            self.reflectionSel.addReflection(
                {"x": peak.xy[0], "y": peak.xy[1]},
                peak.imageno,
                rounded_hkl,
            )
            if (
                intensity_ratio_check
                and peak_intensity is not None
                and peak_intensity["intensity"] > 0.0
            ):
                intensity_scale.append(observed_scale)
            added += 1
            logger.info(
                "Added observed automatic Bragg reflection hkl=%s at image "
                "%s, xy=%s, angular mismatch=%.4g deg, Q mismatch=%.4g "
                "Angstrom^-1%s.",
                rounded_hkl,
                peak.imageno,
                peak.xy,
                angular_mismatch_deg,
                q_mismatch,
                intensity_message,
            )
        return added

    def _automaticPeakAngles(self, peak):
        """Return Vlieg angles in rad for one refined automatic peak."""
        mu, omega = self.getMuOm(peak.imageno)
        gamma, delta = self.ubcalc.detectorCal.surfaceAnglesPoint(
            np.array([peak.xy[1]]),
            np.array([peak.xy[0]]),
            mu,
        )
        return np.array(
            [
                mu,
                float(delta[0]),
                float(gamma[0]),
                omega,
                self.ubcalc.chi,
                self.ubcalc.phi,
            ]
        )

    def _recalculateAutomaticUFromReflections(self):
        """Recalculate U from current reference reflections without dialogs."""
        hkls, angles = self.getReflections()
        if len(hkls) < 1:
            raise ValueError("At least one reflection is required")
        if len(hkls) == 1:
            self.ubcalc.ubCal.zmodeUSingleRefl(angles[0], hkls[0])
        elif len(hkls) == 2:
            self.ubcalc.ubCal.setPrimaryReflection(angles[0], hkls[0])
            self.ubcalc.ubCal.setSecondayReflection(angles[1], hkls[1])
            self.ubcalc.ubCal.calculateU()
        else:
            self.ubcalc.ubCal.calculateUFromReflections(hkls, angles)
        self.ubcalc.uedit.setU(self.ubcalc.ubCal.getU())
        self.ubcalc.updateReflectionMismatch()
        self.ubcalc.sigReplotRequest.emit(False)

    def _automaticReflectionIntensityScale(
        self,
        refl,
        mask=None,
        axis_half_width=1.0,
        fine_axis_half_width=0.4,
        roi_size=(80, 80),
        fine_roi_size=(40, 40),
    ):
        """Return observed/calculated intensity scale for one reflection."""
        candidate = autoBraggSearch.ImageMaximum(
            int(refl.imageno),
            np.asarray(refl.xy, dtype=float),
            np.nan,
        )
        peak = autoBraggSearch.refine_peak_3d(
            self.fscan,
            candidate,
            axis_half_width=axis_half_width,
            roi_size=roi_size,
            fine_axis_half_width=fine_axis_half_width,
            fine_roi_size=fine_roi_size,
            mask=mask,
            background_image=self.background_image,
            excluded_images=self.excludedImagesDialog.getData(),
            max_workers=self.numberthreads,
        )
        intensity = autoBraggSearch.estimate_rocking_intensity(peak)
        f2 = self._automaticReflectionIntensity(refl.hkl)
        if intensity is None or intensity["intensity"] <= 0.0 or f2 <= 0.0:
            return None, peak
        return intensity["intensity"] / f2, peak

    def autoAddCalculatedBraggReflections(self, count, **kwargs):
        """Validate and add calculated Bragg reflections from the current UB.

        :param int count:
            Number of additional reflections to add.
        :param kwargs:
            Automatic Bragg options from
            :class:`orgui.app.QReflectionSelector.AutoBraggOptionsDialog`.
        :returns:
            Number of accepted reflections.
        :rtype: int

        .. note::
           GUI-only. This method mutates the reference-reflection table,
           performs local image peak searches, and emits GUI update signals.
        """
        if self.fscan is None:
            raise ValueError("No scan loaded")
        count = int(count)
        if count <= 0:
            return 0
        status_callback = kwargs.get("status_callback", None)

        def report(event, **fields):
            if status_callback is not None:
                status_callback(event, **fields)

        q_tolerance = kwargs.get(
            "additional_q_tolerance",
            kwargs.get(
                "adaptive_qnorm_tolerance",
                kwargs.get("qnorm_tolerance", 0.05),
            ),
        )
        score_tolerance_deg = kwargs.get("additional_score_tolerance_deg", 2.0)
        intensity_ratio_check = kwargs.get("intensity_ratio_check", True)
        intensity_ratio_tolerance = kwargs.get("intensity_ratio_tolerance", 0.5)
        axis_half_width = kwargs.get("axis_half_width", 1.0)
        fine_axis_half_width = kwargs.get("fine_axis_half_width", 0.4)
        roi_size = kwargs.get("roi_size", (80, 80))
        fine_roi_size = kwargs.get("fine_roi_size", (40, 40))
        pixel_tolerance = kwargs.get(
            "additional_pixel_tolerance",
            kwargs.get("confirmation_pixel_tolerance", 8.0),
        )
        image_tolerance = kwargs.get(
            "additional_image_tolerance",
            kwargs.get("confirmation_image_tolerance", 3),
        )

        mask = self.get_detector_mask()
        if mask is None:
            logger.warning(
                "Adding calculated Bragg reflections without a detector mask."
            )
            report(
                "warning",
                message=(
                    "No detector mask is set; additional calculated "
                    "reflections will be refined without mask rejection."
                ),
            )

        intensity_scale = []
        if intensity_ratio_check:
            for refl in list(self.reflectionSel.reflections):
                try:
                    scale, _peak = self._automaticReflectionIntensityScale(
                        refl,
                        mask=mask,
                        axis_half_width=axis_half_width,
                        fine_axis_half_width=fine_axis_half_width,
                        roi_size=roi_size,
                        fine_roi_size=fine_roi_size,
                    )
                except Exception:
                    logger.debug(
                        "Cannot estimate intensity scale for hkl=%s",
                        refl.hkl,
                        exc_info=True,
                    )
                    continue
                if scale is not None and np.isfinite(scale):
                    intensity_scale.append(float(scale))
            report(
                "additional_intensity_scale",
                message=(
                    f"Built intensity scale list from {len(intensity_scale)} current reference "  # noqa: E501
                    "reflection(s)."
                ),
            )

        added = 0
        rejected = 0
        attempted_hkls = set()
        while added < count:
            candidates = self.reflectionSel.getBraggCandidates(recalculate=True)
            if not candidates:
                break
            accepted_this_round = False
            for candidate_refl in candidates:
                hkl = np.asarray(candidate_refl.hkl, dtype=float)
                key = tuple(np.round(hkl, 8))
                if key in attempted_hkls:
                    continue
                attempted_hkls.add(key)
                report(
                    "additional_candidate",
                    message=(
                        f"Testing calculated Bragg hkl={hkl} at image={candidate_refl.imageno}, "  # noqa: E501
                        f"xy=({candidate_refl.xy[0]:.2f}, {candidate_refl.xy[1]:.2f})."
                    ),
                )
                maximum = autoBraggSearch.ImageMaximum(
                    int(candidate_refl.imageno),
                    np.asarray(candidate_refl.xy, dtype=float),
                    np.nan,
                )
                try:
                    peak = autoBraggSearch.refine_peak_3d(
                        self.fscan,
                        maximum,
                        axis_half_width=axis_half_width,
                        roi_size=roi_size,
                        fine_axis_half_width=fine_axis_half_width,
                        fine_roi_size=fine_roi_size,
                        mask=mask,
                        background_image=self.background_image,
                        excluded_images=self.excludedImagesDialog.getData(),
                        max_workers=self.numberthreads,
                    )
                except Exception:
                    rejected += 1
                    logger.debug(
                        "Cannot refine calculated Bragg hkl=%s",
                        hkl,
                        exc_info=True,
                    )
                    report(
                        "rejected",
                        message=(
                            f"Rejected calculated hkl={hkl}: 3D peak search failed."
                        ),
                    )
                    continue
                pixel_error = float(
                    np.linalg.norm(peak.xy - np.asarray(candidate_refl.xy, dtype=float))
                )
                image_error = abs(int(peak.imageno) - int(candidate_refl.imageno))
                if pixel_error > pixel_tolerance or image_error > image_tolerance:
                    rejected += 1
                    report(
                        "rejected",
                        message=(
                            f"Rejected calculated hkl={hkl}: refined peak moved "
                            f"{pixel_error:.4g} px and {image_error} image(s) (limits {pixel_tolerance:.4g} px, {image_tolerance})."  # noqa: E501
                        ),
                    )
                    continue

                angles = self._automaticPeakAngles(peak)
                try:
                    mismatch = self.ubcalc.getReflectionMismatch(
                        np.atleast_2d(hkl),
                        np.atleast_2d(angles),
                    )
                except Exception:
                    rejected += 1
                    logger.debug(
                        "Cannot score calculated Bragg hkl=%s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                angular_mismatch_deg = float(
                    np.rad2deg(np.asarray(mismatch["angle_mismatch"])[0])
                )
                q_mismatch = float(np.asarray(mismatch["norm_mismatch"])[0])
                report(
                    "additional_mismatch",
                    message=(
                        f"Calculated hkl={hkl}: angular mismatch={angular_mismatch_deg:.4g} deg, "  # noqa: E501
                        f"Q mismatch={q_mismatch:.4g} Angstrom^-1 (limits {score_tolerance_deg:.4g} deg, "  # noqa: E501
                        f"{q_tolerance:.4g} Angstrom^-1)."
                    ),
                )
                if (
                    angular_mismatch_deg > score_tolerance_deg
                    or q_mismatch > q_tolerance
                ):
                    rejected += 1
                    continue

                scale_message = ""
                observed_scale = None
                if intensity_ratio_check:
                    intensity = autoBraggSearch.estimate_rocking_intensity(peak)
                    f2 = self._automaticReflectionIntensity(hkl)
                    if (
                        intensity is None
                        or intensity["intensity"] <= 0.0
                        or f2 <= 0.0
                        or not intensity_scale
                    ):
                        rejected += 1
                        report(
                            "rejected",
                            message=(
                                f"Rejected calculated hkl={hkl}: intensity scale "
                                "could not be tested."
                            ),
                        )
                        continue
                    observed_scale = intensity["intensity"] / f2
                    mean_scale = float(np.mean(intensity_scale))
                    scale_error = abs(observed_scale - mean_scale) / mean_scale
                    report(
                        "additional_intensity_test",
                        message=(
                            f"Calculated hkl={hkl} intensity check: scale={observed_scale:.4g}, "  # noqa: E501
                            f"mean scale={mean_scale:.4g}, relerr={scale_error:.3g} (limit {intensity_ratio_tolerance:.3g})."  # noqa: E501
                        ),
                    )
                    if scale_error > intensity_ratio_tolerance:
                        rejected += 1
                        continue
                    scale_message = f", intensity scale={observed_scale:.4g}, mean scale={mean_scale:.4g}, relerr={scale_error:.3g}"  # noqa: E501

                self.reflectionSel.addReflection(
                    {"x": peak.xy[0], "y": peak.xy[1]},
                    peak.imageno,
                    hkl,
                )
                if intensity_ratio_check and observed_scale is not None:
                    intensity_scale.append(float(observed_scale))
                self._recalculateAutomaticUFromReflections()
                added += 1
                accepted_this_round = True
                logger.info(
                    "Added calculated Bragg reflection hkl=%s at image %s, "
                    "xy=%s, angular mismatch=%.4g deg, Q mismatch=%.4g "
                    "Angstrom^-1%s.",
                    hkl,
                    peak.imageno,
                    peak.xy,
                    angular_mismatch_deg,
                    q_mismatch,
                    scale_message,
                )
                report(
                    "accepted",
                    message=(
                        f"Added calculated hkl={hkl} at image={peak.imageno}, xy=({peak.xy[0]:.2f}, {peak.xy[1]:.2f})"  # noqa: E501
                        f"{scale_message}. Recalculated U."
                    ),
                )
                break
            if not accepted_this_round:
                break

        report(
            "additional_done",
            message=(
                f"Additional calculated Bragg pass complete: added={added}, "
                f"rejected={rejected}."
            ),
        )
        return added

    def _confirmAutomaticBraggSeed(
        self,
        seed,
        mask=None,
        excluded_images=(),
        max_q=None,
        axis_half_width=1.0,
        fine_axis_half_width=0.4,
        roi_size=(80, 80),
        fine_roi_size=(40, 40),
        pixel_tolerance=8.0,
        image_tolerance=3,
        max_reflections=12,
        status_callback=None,
        intensity_ratio_tolerance=0.5,
        intensity_ratio_check=True,
        prominence_threshold=6.0,
    ):
        """Confirm an automatic seed by finding another predicted Bragg peak.

        :param orgui.app.autoBraggSearch.UBSeed seed:
            Candidate seed from the first detected Bragg peak.
        :param numpy.ndarray mask:
            Optional boolean detector mask where ``True`` pixels are invalid.
        :param float max_q:
            Maximum reciprocal-vector norm in Angstrom^-1.
        :returns:
            Confirmation dictionary with ``hkl``, predicted reflection, and
            refined peak, or ``None``.
        :rtype: dict or None

        .. note::
           CLI-capable. Temporarily applies ``seed.U`` only while calculating
           predicted positions and always restores the previous U.
        """
        if max_q is None:
            max_q = getattr(self.ubcalc.detectorCal, "Qmax", None)
        if max_q is None:
            return None
        hkls, _intensity, _norms = autoBraggSearch.allowed_bragg_with_intensity(
            self.ubcalc.crystal, max_q
        )
        seed_intensity = autoBraggSearch.estimate_rocking_intensity(seed.peak)
        seed_f2 = self._automaticReflectionIntensity(seed.hkl)
        previous_u = np.asarray(self.ubcalc.ubCal.getU()).copy()
        try:
            self.ubcalc.ubCal.setU(seed.U)
            tested = 0
            for hkl in hkls:
                if np.allclose(hkl, seed.hkl):
                    continue
                try:
                    refldict = self.searchPixelCoordHKL(hkl)
                except Exception:
                    logger.debug(
                        "Cannot predict confirmation reflection %s",
                        hkl,
                        exc_info=True,
                    )
                    continue
                for solution in (1, 2):
                    if not refldict.get(f"selectable_{solution}", False):
                        continue
                    if f"imageno_{solution}" not in refldict:
                        continue
                    predicted_xy = np.asarray(refldict[f"xy_{solution}"], dtype=float)
                    predicted_imageno = refldict[f"imageno_{solution}"]
                    candidate = autoBraggSearch.ImageMaximum(
                        int(predicted_imageno),
                        predicted_xy,
                        np.nan,
                    )
                    try:
                        peak = autoBraggSearch.refine_peak_3d(
                            self.fscan,
                            candidate,
                            axis_half_width=axis_half_width,
                            roi_size=roi_size,
                            fine_axis_half_width=fine_axis_half_width,
                            fine_roi_size=fine_roi_size,
                            mask=mask,
                            background_image=self.background_image,
                            excluded_images=excluded_images,
                            max_workers=self.numberthreads,
                        )
                    except Exception:
                        logger.debug(
                            "Cannot refine confirmation reflection %s",
                            hkl,
                            exc_info=True,
                        )
                        continue
                    pixel_error = np.linalg.norm(peak.xy - predicted_xy)
                    image_error = abs(peak.imageno - int(predicted_imageno))
                    if (
                        pixel_error <= pixel_tolerance
                        and image_error <= image_tolerance
                    ):
                        confirm_intensity = autoBraggSearch.estimate_rocking_intensity(
                            peak
                        )
                        confirm_f2 = self._automaticReflectionIntensity(hkl)
                        intensity_ok = True
                        observed_ratio = np.nan
                        predicted_ratio = np.nan
                        ratio_error = np.nan
                        seed_prominence = np.nan
                        confirm_prominence = np.nan
                        if seed_intensity is not None:
                            seed_prominence = seed_intensity["prominence_z"]
                        if confirm_intensity is not None:
                            confirm_prominence = confirm_intensity["prominence_z"]
                        if intensity_ratio_check:
                            if (
                                seed_intensity is not None
                                and confirm_intensity is not None
                                and seed_intensity["intensity"] > 0.0
                                and confirm_intensity["intensity"] > 0.0
                                and seed_f2 > 0.0
                                and confirm_f2 > 0.0
                            ):
                                observed_ratio = (
                                    confirm_intensity["intensity"]
                                    / seed_intensity["intensity"]
                                )
                                predicted_ratio = confirm_f2 / seed_f2
                                ratio_error = (
                                    abs(observed_ratio - predicted_ratio)
                                    / predicted_ratio
                                )
                                intensity_ok = ratio_error <= intensity_ratio_tolerance
                            else:
                                intensity_ok = False
                        else:
                            intensity_ok = (
                                np.isfinite(seed_prominence)
                                and np.isfinite(confirm_prominence)
                                and seed_prominence >= prominence_threshold
                                and confirm_prominence >= prominence_threshold
                            )
                        if not intensity_ok:
                            logger.debug(
                                "Reject confirmation hkl=%s for seed hkl=%s "
                                "because intensity check failed. observed=%s "
                                "predicted=%s relerr=%s seed_prominence=%s "
                                "confirm_prominence=%s",
                                hkl,
                                seed.hkl,
                                observed_ratio,
                                predicted_ratio,
                                ratio_error,
                                seed_prominence,
                                confirm_prominence,
                            )
                            continue
                        return {
                            "hkl": np.asarray(hkl, dtype=float),
                            "solution": solution,
                            "predicted_xy": predicted_xy,
                            "predicted_imageno": int(predicted_imageno),
                            "pixel_error": pixel_error,
                            "image_error": image_error,
                            "peak": peak,
                            "seed_intensity": seed_intensity,
                            "confirmation_intensity": confirm_intensity,
                            "observed_intensity_ratio": observed_ratio,
                            "predicted_intensity_ratio": predicted_ratio,
                            "intensity_ratio_error": ratio_error,
                            "seed_prominence_z": seed_prominence,
                            "confirmation_prominence_z": confirm_prominence,
                        }
                tested += 1
                if status_callback is not None:
                    status_callback(
                        "confirmation_search",
                        message=(
                            f"Tried {tested} confirmation reflections for seed "
                            f"hkl={seed.hkl}."
                        ),
                    )
                if tested >= max_reflections:
                    break
        finally:
            self.ubcalc.ubCal.setU(previous_u)
        return None


def auto_find_bragg_reference(gui, **kwargs):
    """Search the active scan for a Bragg reference and seed the UB matrix."""
    return AutoBraggWorkflow(gui).autoFindBraggReference(**kwargs)


def auto_add_calculated_bragg_reflections(gui, count, **kwargs):
    """Validate and add calculated Bragg reflections from the current UB."""
    return AutoBraggWorkflow(gui).autoAddCalculatedBraggReflections(count, **kwargs)

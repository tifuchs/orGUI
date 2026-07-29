Automatic UB Setting
====================

The automatic UB workflow searches a scan for a Bragg peak, assigns a plausible
``hkl``, constructs an initial ``U`` matrix, and validates that orientation by
finding another Bragg peak predicted by the trial matrix. The goal is to obtain
a reliable starting orientation without first manually identifying two
reference reflections.

The controls are in the ``Auto UB/Reflections`` group of the
``Reciprocal space navigation`` panel.

Prerequisites
-------------

Before starting the automatic search, configure the same scientific state that
is required for manual UB work:

1. Load the scan that should be searched.
2. Load or enter the detector calibration.
3. Load or enter the crystal lattice and structure model. Calculated
   reflection norms and relative intensities come from this model.
4. Set the X-ray energy or wavelength.
5. Set the diffractometer geometry in the advanced UB settings. For surface
   diffraction this usually means the grazing-incidence / z-mode geometry.
6. Set a detector mask if one is available. The automatic search can run
   without a mask, but masked-pixel rejection and mask-aware peak refinement are
   disabled. In GUI mode orGUI asks before continuing without a mask.

The default ``hkl`` assignment path uses the current ``U`` matrix to predict
detector positions. This does not require a correct final orientation, but the
current matrix should be a reasonable geometry convention for the mounted
sample. The slower Q-norm path can be selected when the current detector
prediction is not useful.

Overview of the Automatic Workflow
----------------------------------

The ``auto Bragg/UB`` play button opens the automatic UB/reflection status
dialog. The dialog can be used in two ways:

* Click ``seed UB`` to run the full automatic UB seeding workflow.
* Click ``Add calculated`` to add calculated Bragg reflections from an already
  seeded ``UB`` matrix.

The ``seed UB`` button performs the following sequence:

1. Stream through the scan images and watch the one-dimensional max-pixel trace.
   Images are read only until a sharp candidate is detected and validated.
2. Refine the candidate by a local three-dimensional peak search over detector
   position and scan-axis position.
3. Convert the refined peak to a measured momentum-transfer norm and generate
   candidate ``hkl`` assignments.
4. For each candidate ``hkl``, calculate a one-reflection ``U`` matrix using
   the selected default geometry.
5. Predict additional allowed Bragg reflections from the trial ``U`` matrix.
6. Refine the predicted confirmation reflection by another local
   three-dimensional peak search.
7. Accept the seed only if the confirmation peak is close enough in detector
   position and scan image.
8. Validate the seed and confirmation intensities against the calculated
   structure-factor intensities, when intensity validation is enabled.
9. Add the seed and confirmation reflections to the reference-reflection table,
   set the accepted ``U`` matrix, and update the mismatch display.
10. Revisit other refined candidate peaks that were already observed during the
    search. For each, calculate ``hkl`` from the final angles, round to the
    nearest integer, and add it only if mismatch and intensity checks pass.

The status dialog stays open during this process and prints the candidate
images, refined positions, adaptive broadening events, confirmation attempts,
intensity checks, and accepted reflections.

Automatic Bragg Options
-----------------------

The gear button opens ``Automatic Bragg options``. The options are grouped
below by the stage of the algorithm where they apply.

HKL Assignment
~~~~~~~~~~~~~~

``HKL assignment``
   Selects how the first refined peak is assigned to candidate Miller indices.

   ``Detector position (fast)``
      The default path. orGUI calculates allowed reflections near the measured
      Q shell and uses the current ``U`` matrix to predict their detector
      positions. Candidates whose predicted detector positions are near the
      observed peak are tried first. If one ``hkl`` has multiple
      ``delta``/``gamma`` solutions, each selectable solution is considered.

   ``Q norm shell (slower)``
      Uses only the measured Q norm to collect candidate reflections. This can
      test many more hypotheses because it does not filter by predicted
      detector position. It is useful when the current ``U`` matrix is too poor
      for detector-position prediction.

``Q shell half-width``
   The initial half-width, in :math:`\mathrm{\AA}^{-1}`, around the measured
   peak Q norm used to select possible Bragg reflections.

``Assignment pixel tolerance``
   Used in the detector-position assignment path. A first-peak candidate must
   have a predicted detector position within this distance, in pixels, of the
   observed refined peak.

``Assignment reflections``
   Maximum number of first-peak reflection hypotheses taken from the local Q
   shell before ranking and testing.

``Seed hypotheses``
   Maximum number of first-peak ``hkl`` hypotheses sent to the confirmation
   stage. In detector-position mode, same-Q ambiguity is considered before this
   limit is applied.

Confirmation Peak Search
~~~~~~~~~~~~~~~~~~~~~~~~

After a one-reflection ``U`` matrix is calculated, orGUI predicts other allowed
Bragg reflections and tries to verify one of them in the data.

``Confirmation pixel tolerance``
   Maximum allowed detector-position difference, in pixels, between the
   predicted confirmation position and the refined confirmation peak.

``Confirmation image tolerance``
   Maximum allowed image-index difference between the predicted confirmation
   image and the refined confirmation peak.

``Confirmation reflections``
   Maximum number of predicted strong reflections to test for confirmation of
   one seed hypothesis. The predicted reflections are ordered by calculated
   relative intensity.

Adaptive Broadening and Q-Scale Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If early candidates fail because the configured lattice constants or orientation
are not close enough, the search broadens the matching tolerances instead of
stopping immediately.

``Broaden after unmatched``
   Number of refined candidates that may fail before the broader tolerances are
   activated.

``Broader Q shell half-width``
   Q-shell half-width, in :math:`\mathrm{\AA}^{-1}`, used after broadening.
   This broader range is also used when checking already observed peaks after
   the seed has been confirmed.

``Broader assignment tolerance``
   Detector-position tolerance, in pixels, used for first-peak assignment after
   broadening.

``Broader confirmation tolerance``
   Detector-position tolerance, in pixels, used for the second Bragg peak after
   broadening.

``Broader confirmation image tolerance``
   Image-index tolerance used for the second Bragg peak after broadening.

``Scale-fit detector filter``
   Enables a broad detector-position filter when estimating whether the
   configured reciprocal-lattice scale is systematically wrong. In detector
   position mode, only provisional reflections whose current-``U`` prediction
   lands in a broad region around the observed peak are used for the scale
   estimate.

``region half-size``
   Half-width and half-height of that broad detector region, expressed as a
   fraction of the detector size. For example, ``0.25 detector`` means a region
   extending one quarter of the detector width and height around the observed
   peak.

``Scale-fit Q outlier tolerance``
   Maximum scaled Q-norm residual, in :math:`\mathrm{\AA}^{-1}`, for
   provisional reflections kept in the adaptive scale-fit retry.

``delta/gamma fraction``
   Maximum detector-angle residual for provisional scale-fit reflections,
   expressed as a fraction of the detector angular range. This removes obvious
   detector-position outliers before retrying validation with the tentative
   scale.

Peak Refinement Windows
~~~~~~~~~~~~~~~~~~~~~~~

Every candidate peak is refined by a local three-dimensional search. The same
settings are used for the seed candidate, confirmation peak, already observed
candidate peaks, and automatically added calculated reflections.

``Peak search axis half-width``
   Coarse half-width of the local scan-axis search window, usually in degrees
   of the scan motor.

``Fine axis half-width``
   Fine half-width of the second local scan-axis search window around the
   coarse result.

``Peak search ROI vertical`` and ``horizontal``
   Coarse detector ROI size in pixels.

``Fine ROI vertical`` and ``horizontal``
   Fine detector ROI size in pixels for the second pass.

Larger windows are more tolerant of a poor prediction, but they require more
image access and increase the chance that the local search locks onto the wrong
feature.

Intensity Validation
~~~~~~~~~~~~~~~~~~~~

``Intensity validation``
   When enabled, the seed and confirmation reflections must have integrated
   rocking-curve intensities whose ratio is consistent with the calculated
   structure-factor intensity ratio ``abs(F)**2``. The default tolerance is a
   50 percent relative error.

   For additional observed or calculated reflections, orGUI builds a list of
   observed intensity scales,

   .. math::

      S_i = I_{\mathrm{obs},i} / |F_i|^2,

   from the already accepted reflections. A new reflection is accepted only if
   its scale is within 50 percent of the current mean scale. Accepted
   reflections are appended to the scale list before testing the next
   reflection.

``prominence fallback``
   Used only when calculated intensity-ratio matching is disabled. In that
   mode, both the seed and confirmation rocking curves must have at least this
   robust side-band prominence, in z-score units.

Status Dialog and Additional Reflections
----------------------------------------

The ``auto Bragg/UB`` play button opens a non-modal status dialog. It reports
running statistics and appends one line for each important event: candidate
detection, peak refinement, HKL hypotheses, confirmation attempts, adaptive
broadening, intensity tests, accepted reflections, and failures.

``seed UB``
   Starts the automatic max-trace peak search, one-reflection ``U`` seeding,
   confirmation-reflection search, and optional post-seed reflection addition.

``add calculated after seed``
   When checked, orGUI automatically runs the calculated Bragg-reflection add
   pass after ``seed UB`` succeeds. This is enabled by default. Disable it when
   you want to inspect the seed and confirmation reflections before adding more
   references.

``Additional reflections``
   Number of additional calculated Bragg reflections to add from the current
   ``UB`` matrix. The default is ``50``.

``add all Bragg``
   Ignores the number field and keeps testing calculated Bragg reflections
   until no further candidates pass peak refinement, mismatch, and optional
   intensity validation.

``Add calculated``
   Calculates allowed Bragg reflections from the current ``UB`` matrix using
   the same candidate ordering as the manual ``add Bragg reflection`` button.
   For each candidate, orGUI refines the real peak position, checks detector
   movement, checks angular and Q mismatch, applies the same optional intensity
   scale validation, adds accepted reflections, and recalculates ``U`` after
   every accepted reflection.

``Reset``
   Clears the status log and statistics in the dialog. It does not remove
   reference reflections and does not change the current ``UB`` matrix.

What Is Added to the Reflection Table
-------------------------------------

When the automatic seed succeeds, orGUI adds at least two reflections:

1. The seed reflection used to construct the first trial ``U`` matrix.
2. The independently found confirmation reflection that verified the trial
   matrix.

orGUI may then add additional observed peaks that were already refined during
the max-trace search. These are not assigned by a Q-shell search. Instead,
orGUI calculates continuous ``hkl`` from the final measured angles, rounds to
the nearest integer ``hkl``, and accepts the reflection only if the angular
mismatch, Q mismatch, and optional intensity-scale check pass.

The ``Add calculated`` button can add more reflections after the initial seed.
Those reflections start from calculated Bragg positions rather than from
previously observed max-trace candidates. The measured peak position is still
found by local 3D peak refinement before the reflection is accepted.

Failure Modes and Practical Checks
----------------------------------

If no candidate is accepted, inspect the status dialog for the first repeated
failure mode:

``no HKL hypotheses in local Q region``
   The measured peak Q norm is not close to any allowed reflection under the
   current lattice constants and Q-shell tolerance. Broadening or correcting
   the lattice constants may be required.

``no second Bragg peak confirmed``
   A seed candidate could be constructed, but its predicted confirmation peak
   was not found within the detector/image tolerances. Check the geometry mode,
   current ``U`` convention, and confirmation tolerances.

``intensity check failed``
   The geometry looked plausible, but the observed rocking-curve intensity was
   inconsistent with the calculated relative intensity. This can indicate a
   wrong ``hkl`` assignment, a poor background estimate, a saturated or masked
   peak, or a structure model that does not describe the measured sample.

``running without a detector mask``
   The search can proceed, but bright bad pixels or detector artifacts may be
   selected as candidates. Set a mask when possible.

After automatic setup, validate the result exactly as for manual UB setting:
inspect the reference-reflection mismatch table, step through the scan with
calculated Bragg or CTR overlays enabled, and confirm that predicted positions
track real peaks across the detector and scan axis.

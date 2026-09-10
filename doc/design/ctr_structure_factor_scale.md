# One structure-factor scale for rocking scans, stationary scans and reflectivity

> **Status: analysis complete, reduction landed and wired into both
> integration paths.** This is the review document for
> [issue #82](https://github.com/tifuchs/orGUI/issues/82) ("Regression tests
> and validation of equivalence of rocking and stationary scan integration")
> and for the physics half of
> [issue #15](https://github.com/tifuchs/orGUI/issues/15) ("Calculate
> quantitatively exact structure factors").
>
> It records why a rocking scan and a stationary scan of the same rod did
> **not** produce the same `F2_hkl` in orGUI, quantifies each reason on
> simulated data, and states what is still needed for an absolute scale and
> for absolute reflectivity. Everything below was verified numerically against
> the code, not by inspection alone.
>
> F1, F2, F3, F5 and F6 are now **applied**, which changed saved numbers in
> both modes; the two paths agree to `1e-6` on simulated data, limited by the
> trapezoidal sampling of the rocking profile. F4 and F7 remain open and are
> described below as they stand. Findings are written in the present tense
> of the analysis; section 5 says what each one's status is now, and
> [`ctr_structure_factor_handover.md`](ctr_structure_factor_handover.md)
> section 5 says how each was wired.

## 1. What the two papers require

E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532 gives the integrated intensity
for each measurement mode of a six-circle/z-axis diffractometer. J. Drnec
*et al.*, *J. Appl. Cryst.* **47** (2014) 365 rewrites the same expressions
for a two-dimensional detector and states the equivalence claim of issue #82
directly (their section 3.4 and Fig. 8, right: rocking scans at low *l* and a
stationary scan at high *l* overlap once both are reduced properly).

**Rocking scan** (Vlieg eq. 42/43, Drnec eq. 3):

```
I_omega = Phi_0 T_omega (r_e^2 A lambda^2 / A_u^2) |F|^2 P L_phi C_rod C_det C_beam Delta_gamma
```

with `L_phi = 1/(sin(delta) cos(beta_in) cos(gamma))` (eq. 16), `C_rod` the rod
interception, which reduces to `cos(gamma)` in the z-axis mode (eq. 23), and
`Delta_gamma` the **out-of-plane angular acceptance of the detector aperture**.
`I_omega` is the integral of the per-frame counts over the rocking angle **in
radian**, divided by the per-frame counting time.

**Stationary measurement** (Vlieg eq. 54, Drnec eq. 5):

```
I_s = Phi_0 T_s (r_e^2 A lambda^2 / A_u^2) |F|^2 P L_s C_det,s C_beam
```

with `L_s = 1/sin(beta_out)`, in the z-axis mode `1/sin(gamma)`. There is **no**
rod-interception factor and **no** `Delta_gamma`: a stationary measurement
intercepts the whole rod cross-section rather than a slice of it.

**Specular reflectivity** (Vlieg section 3.2, eq. 62/63) is the stationary case
at `gamma = alpha` with a beam-limited active area, `A_r = A_0/sin(alpha)`:

```
R = I_s / (Phi_0 T A_0 C_beam) = r_e^2 lambda^2 P_r |F|^2 / (A_u^2 sin^2(alpha))
```

One power of `1/sin(alpha)` is the footprint, the other is `L_s`. Nothing
beyond `|F|^2` is needed, so **an absolutely scaled `|F|^2` yields an absolute
reflectivity for free** - see section 6.

Their ratio (Vlieg eq. 64/65) is the identity issue #82 is really about:

```
I_s / I_omega = T omega_0 sin(delta) cos(beta_in) / (Delta_gamma sin(gamma))
```

## 2. What orGUI did before this work

| | rocking (`peak1Dintegr`) | stationary (`orGUI.integrateROI` + `integration_corrections`) | reconstruction (`reconstruction_job`) |
|---|---|---|---|
| solid angle | per-pixel array, applied as ROI mean | per-pixel array, applied as ROI mean | per pixel |
| polarization | per-pixel array, applied as ROI mean | per-pixel array, applied as ROI mean | per pixel |
| exposure / monitor | **not applied** | `normalization_divisor` | applied |
| footprint `C_beam` | `beamprofile` `C_illum_area` | `beamprofile` `C_illum_area` | not applied |
| `C_area` (`1/sin(delta)`) | not applied (documented as deliberate) | not applied (documented as deliberate) | not applied |
| Lorentz | `1/(sin(delta) cos(alpha) cos(gamma))` | `1/sin(gamma)` | not applied (correct: the voxel binning absorbs it, Drnec section 4) |
| rod interception | `cos(gamma)` | none (correct) | none |
| `Delta_gamma` | **not applied** | not applicable (correct) | `Delta_l` not applied |
| rocking angle unit | **degrees** | n/a | n/a |
| `C_det` | not modelled | not modelled | not modelled |

The two Lorentz factors, the rod interception and the absence of a stationary
`C_rod` are all correct and match the ANA/ROD z-axis table. The disagreement
comes from the three cells in bold.

## 3. Findings

Ranked by how badly they distort the saved numbers.

### F1 - rocking scans carry no exposure-time or monitor normalization

`orgui/app/integration_corrections.py:normalization_divisor` is wired into the
stationary path only (`orGUI.py:6361`). `peak1Dintegr.py` contains no
occurrence of `exposure`, `monitor` or `normaliz` at all. Both published
expressions carry `Phi_0 T`, so a rocking integral and a stationary frame are
not comparable until both are divided by counting time and monitor.

Worse than the mode mismatch: within a *single* rocking data set, rocking scans
taken at different *l* with different counting times or a drifting ring current
are put on different scales. Nothing in the saved data records this.

### F2 - the rocking integral is taken in degrees

`_compute_rocking_integration` trapezoid-integrates `croibg` over `axis`, which
is the rocking motor position in degrees (`peak1Dintegr.py:235`-`243`). Vlieg
eq. 42 and Drnec eq. 2 integrate in radian. Factor `180/pi = 57.3`, constant,
so it is invisible inside one data set and only shows up when comparing to
another mode or to an absolute scale.

### F3 - no out-of-plane detector acceptance `Delta_gamma` for rocking scans

This is the physically interesting one. A rocking scan intercepts a slice of
rod of length `Delta_l = C_rod (V_u / lambda A_u) Delta_gamma` (Vlieg eq. 20),
so its integrated intensity is *proportional to the vertical angular size of
the region of interest*. A stationary measurement has no such factor.

Before this branch, `Delta_gamma` appeared nowhere in orGUI: `grep -rn
"acceptance"` over `orgui/` returned nothing. It is now estimated by
`corrections/acceptance.py` (section 4.4), though the integration paths do
not call it yet. In the point-detector world it was a fixed slit setting and
disappeared into the overall scale factor. On an area detector with ROIs sized
per detector position by `orgui/app/ROIutils.py:calc_corrections` (projected
sample size plus a parallax correction, so it changes with `delta` and
`gamma`), it is **not constant along a scan**, and the omission changes the
*shape* of the rod, not just its scale.

`test_a_resized_region_of_interest_distorts_the_rocking_rod` simulates a rod
measured with the acceptance sweeping from 0.7x to 1.6x of nominal and shows
the distortion carried straight into `F2_hkl`: a factor 2.3 across the rod.

### F1+F2+F3 combined - the measured gap

`orgui/app/test/test_scan_mode_equivalence.py::test_rocking_and_stationary_paths_differ_by_the_missing_normalizations`
simulates a Pt(111) `(1, 0, l)` rod at 17.7 keV from a known `|F|^2(l)`, using
Vlieg eq. 42 and 54 as the forward model, and pushes the result through both of
orGUI's current paths. The ratio of the two `F2_hkl` is, to seven digits,

```
F2_rocking / F2_stationary = T_omega * monitor_omega * Delta_gamma_in_degrees
```

(measured 0.5249999866 against a predicted 0.525). The `180/pi` of F2 and the
`Delta_gamma` in radian of F3 combine into the acceptance expressed in degrees;
the counting time and monitor of F1 sit in front. For a realistic setup - a
100-pixel-tall ROI at 1 m with 172 um pixels is `Delta_gamma ~ 1 degree`, with
a per-frame time of a few tenths of a second - the two modes land within a
factor of a few of each other, which is exactly the regime in which the
discrepancy looks like a plausible scale factor rather than a bug.

### F4 - the active area is a dimensionless fraction, and `C_area` is skipped

`geometrycorrections.area_correction` (`1/sin(delta)`) exists but is
deliberately not applied; the numerical `beamprofile` factors are used instead.
`BeamProfile.illuminated_area_fraction` returns
`C_flux_on_sample / (p_max * L sin(alpha))` - the beam profile integrated over
the projected sample, normalized to a uniform beam. It is a *fraction*, not an
area.

This is **not** a mode-equivalence problem: the active area enters both
expressions identically and cancels from their ratio
(`test_footprint_area_cancels_between_the_modes`). It is a problem for
issue #15, and it hides an undocumented assumption:

* With **open post-sample slits and an area detector** - orGUI's usual
  configuration - the illuminated footprint defines the active area, it has no
  `delta` dependence, and skipping `C_area` is right.
  `test_the_stationary_path_recovers_the_rod_up_to_one_constant` confirms the
  stationary path then recovers the rod shape to machine precision.
* With **slits narrower than the footprint**,
  `C_area = 1/(sin(delta) cos(alpha - beta_in))` applies and varies along the
  rod. On the simulated rod above, `delta` moves from 16.77 to 16.95 degrees
  between `l = 0.4` and `l = 3.0`, a 1.1 % shape error - small here, much
  larger between rods at different in-plane momentum transfer.

Neither the assumption nor which case a given data set is in is recorded
anywhere in the output.

### F5 - the polarization correction is evaluated at the calibrated arm position

`orGUI.py:1590` and `orGUI.py:5752` build the per-pixel correction array once,
outside the frame loop, from `dc.polarizationArray()`. That is pyFAI's
detector-frame expression evaluated at the *calibrated* geometry.
`DetectorCalibration.polarizationAtPoints` exists precisely because that is
wrong for a moving detector arm, and says so in its own docstring; per-frame
arm angles are already resolved by `backend/scans.py:scan_arm_angles` and used
by `getROIloc` and the cursor readout.

Measured at the detector centre for a specular scan whose arm follows
`gamma = 2 alpha`, horizontally polarized beam:

| 2theta | P (calibrated position) | P (arm-following) | error |
|---|---|---|---|
| 1 deg | 1.00000 | 0.99970 | +0.03 % |
| 4 deg | 1.00000 | 0.99513 | +0.49 % |
| 10 deg | 1.00000 | 0.96985 | +3.11 % |
| 18 deg | 1.00000 | 0.90451 | +10.56 % |
| 30 deg | 1.00000 | 0.75000 | +33.33 % |

For a **fixed** arm the static array is correct - each pixel already carries
its own scattering angle, and the apparent `alpha` dependence of the ANA
z-axis expression is only a change of frame. The bug is scoped to scans that
move the arm, which is exactly the reflectivity case of section 6.

**Fixed**, as a per-frame multiplicative factor rather than by rebuilding the
array. `detector.polarization_arm_correction` returns
`<1/P_arm> / <1/P_home>` over the region of interest, which is what an
intensity already corrected at the calibrated position has to be multiplied
by, and both integration paths apply it. Three properties made this the
proportionate fix rather than a per-frame rebuild of the full array:

* It is **exactly 1** when the arm sits at its calibrated reference, because
  both evaluations then use the same geometry. A fixed-arm scan is therefore
  bit-identical, and no switch or special case is needed to keep it that way.
* It costs a region-sized evaluation per frame instead of a detector-sized
  one, capped at a 17x17 sample of the region because the polarization is
  smooth across it. A constant arm is evaluated once and broadcast.
* It leaves the fused per-pixel array, the native ROI accumulation and the
  stored `Cfactors_croi` untouched.

Re-measured through the shipped factor, along a specular scan whose arm
follows `gamma = 2 alpha`, it reproduces the table above: `+0.03 %` at
`2theta = 1`, `+0.49 %` at 4, `+3.11 %` at 10, `+10.57 %` at 18 and
`+33.35 %` at 30 degrees.

Two things it does **not** cover. The polarization and the solid angle are
fused into one array by `pixel_factors`, so the ratio of region means leaves
their covariance over the region - second order in the variation of both
across it. And the **reciprocal-space reconstruction** applies the
polarization per pixel rather than as a region mean, so this factor does not
apply to it; that path would need the array itself rebuilt per frame, which is
a detector-sized evaluation per frame and belongs with the rest of the
reconstruction work in section 7.

### F6 - solid-angle correction applied to an already-summed ROI

For a ROI-summed intensity, the raw sum over pixels *is* the angular integral
`Phi_0 T Int Int (dsigma/dOmega) dgamma dpsi`; no solid-angle correction is
needed. orGUI multiplies the sum by the ROI-mean of `1/solidAngleArray`
(`integration_corrections.roi_mean_correction`), which introduces a factor
`Omega_ref / <Omega_ROI>` that varies over the detector - a few tenths of a
percent for a detector at 1 m, several percent for a close-in detector.

The solid-angle array and an angle-derived `Delta_gamma` describe the same
geometry. Applying one without the other double-counts, so the pair has to be
decided together with F3.

**Resolved: the solid-angle correction stays available as a switch, scales the
intensity, and is divided back out when `F2_hkl` is formed.** It is kept
because it is the right correction for a *broad or diffuse* feature, where a
differential cross-section rather than an integrated rod intensity is wanted;
it is removed from the structure factor because a rod is integrated by summing
a region, which is already the complete angular integral. The per-pixel
reconstruction keeps applying it without compensation, under its own switch
(`reconstruction_job.py:1254`).

An earlier revision of this document had the correction dropped from ROI
integration altogether. That was reversed: it discarded a capability that has
a real use, and the compensation costs one geometric factor. The physics below
is unchanged - what changed is where the factor is removed, not whether.

Three steps to that conclusion:

1. **A stationary ROI sum is a complete integral.** With the reflection inside
   the region - open post-sample slits and a region sized from the projected
   sample, orGUI's usual configuration - every photon is counted once, already
   weighted by the solid angle of the pixel that caught it. There is nothing
   to correct.
2. **Closing the slits would not make it the right correction either.** A
   truncated integral is Vlieg's `C_det` (F7), an aperture-against-profile
   model. A per-pixel obliquity weight does not repair a cut-off tail, so the
   conclusion does not rest on the open-slit assumption holding.
3. **A rocking sum is a complete integral too**, over a `Delta_gamma` slice of
   the rod. The slice is what `Delta_gamma` accounts for; there is still no
   per-pixel weighting to undo.

The compensation is
`detector.roi_mean_inverse_solid_angle`, the mean of `1/Omega~` over the
nominal region rectangle, evaluated from the calibrated geometry alone so it
needs no image and can be computed per frame outside the integration loop.
Two approximations come with that, both far below the effect being removed:
it is the mean over the *nominal* region rather than over the valid pixels the
integration accumulated, and because `pixel_factors` fuses the solid angle
with the polarization into one array, the applied mean is
`<1/(Omega~ P)>` rather than `<1/Omega~><1/P>`, leaving their covariance over
the region - second order in the variation across it, of order `1e-6`.

What the correction cannot do is stay in *uncompensated*, because the
two modes do not carry it symmetrically: for a rocking scan the
omega-integrated counts in pixel row *j* go as that row's gamma height
`dgamma_j`, so a corrected sum pairs with `Sum_j dgamma_j / Omega~_j` and the
factor cancels against the acceptance, while a stationary sum has no
`Delta_gamma` to cancel against and the factor survives into the saved number.
The two modes would still disagree by `<1/Omega~>` - 0.7 % at 1 m, 7 % for a
detector at 0.3 m - after F1-F3 were fixed.

Measured on a calibrated 172 um detector with a 60-row region, where
`Omega~ = solidAngleArray` is `cos^3` of the incidence angle on the detector
face to `1e-9`, and `nominal = n * pixel / dist`:

| dist | theta | `Delta_gamma_edge`/nominal | exact/nominal | exact/`Delta_gamma_edge` | `1/<Omega~>` |
|---|---|---|---|---|---|
| 1.0 m | 2.7 deg | 0.99774 | 1.00112 | 1.003383 | 1.003383 |
| 0.3 m | 5.6 deg | 0.99039 | 1.00476 | 1.014512 | 1.014509 |
| 0.3 m | 9.0 deg | 0.97549 | 1.01236 | 1.037795 | 1.037788 |

The last two columns agreeing to `5e-6` is the cancellation. In closed form
`Delta_gamma_edge = nominal cos^2(theta)` and `Omega~ = cos^3(theta)`, so the
partner of a corrected sum is `nominal / cos(theta)`.

Two consequences that survive the decision, because they are what would have
gone wrong under the other one:

* **`Delta_gamma` does not go away.** It is a *size* factor - how many rows of
  rod the region accepted - while the solid angle is an *obliquity* factor.
  They overlap only in the obliquity. `acceptance.out_of_plane_acceptance` is
  still what the rocking path divides by.
* **A constant nominal `n * pixel / dist` is not the "pixel-derived"
  acceptance.** It is wrong by `1/cos(theta)`: 0.1 % at 1 m, 1.2 % at 9
  degrees. `pixel_acceptance` times the row count is the same angle-derived
  quantity as `out_of_plane_acceptance`, not a cheaper alternative to it.

On the user-facing side the switch and its config key are untouched, so
nothing about existing configurations changes. Note that `useSolidAngleBox`
is *also* the store the reciprocal-space reconstruction persists its own
solid-angle setting through: `config_data.py:383`/`:467` map the `solidAngle`
integration option onto `CorrectionState.use_solid_angle`, and
`ReconstructionDialog` mirrors its checkbox via
`scanSelector.get/set_integration_options`. Removing the key would therefore
have silently disabled the correction in the one path that must keep it - a
trap worth recording, because it is invisible from either file alone.

What the switch means is now stated in its tooltip and in the `SOLA` badge,
since "applied to the intensity but not to `F2_hkl`" is not something a user
can infer from a checkbox.

### F7 - `C_det` is assumed to be 1 in both modes

Neither path models the in-plane detector acceptance (Vlieg sections 2.4/2.5,
eq. 27/28). It is close to 1 for a wide ROI and a sharp in-plane profile, but
it is **not the same** for the two modes when the rod is broad: Drnec's Fig. 18
shows direct stationary integration underestimating `|F|` at low *l* by up to a
factor of 2 relative to reciprocal-space integration, for exactly this reason.
This is the one remaining mechanism that can make the two modes disagree
*after* F1-F3 are fixed, and it is a modelling problem, not a normalization
one.

## 4. What landed

### 4.1 One home for the corrections

Before this branch the answer to "which factors apply to this scan" existed in
four places, "what is the active area" in three, "what is the exposure and
monitor divisor" in three, and the five-line per-pixel solid-angle and
polarization block in three verbatim copies. Adding the reduction as one more
module would have made that five, four, four and three.

Everything now lives in `orgui/datautils/xrayutils/corrections/`, split by
*what a factor depends on*, which is also what makes each piece testable on
its own:

| module | depends on |
|---|---|
| `geometry.py` | diffractometer angles only (the ANA/ROD z-axis table) |
| `beamprofile.py` | the vertical profile of the incident beam |
| `activearea.py` | beam, sample and slits |
| `detector.py` | the calibrated detector geometry, per pixel |
| `normalization.py` | counting time and monitor values |
| `roi.py` | reducing per-pixel factors onto a summed region |
| `measurement.py` | which of the above a scan mode applies |

The package is physics: numbers in, numbers out. Nothing in it reads a scan
object, a configuration or a widget. `geometrycorrections.py` and
`beamprofile.py` at the old level were released under those names and remain
as re-export aliases; `test_corrections_package.py` checks that the aliases
hand out the *same objects*, not merely importable ones.

`orgui/app/integration_corrections.py` is now the adapter in the other
direction -- switch states and scan counter names in, arguments out -- and
that boundary is written into `orgui/app/AGENTS.md`. The callers were rewired:
`orGUI.py` builds both of its per-pixel arrays with
`corrections.detector.pixel_factors`, `reconstruction_job.py` uses the same
function while keeping its own native-fused application, and
`peak1Dintegr.py` asks `measurement.mode_components` for its Lorentz and
rod-interception factors instead of selecting them itself.

### 4.2 The reduction

`corrections/measurement.py` implements

```
I = Phi_0 (r_e^2 A lambda^2 / A_u^2) |F|^2 eta C_det
```

* `normalized_intensity(counts, exposure_time, monitor, angle_unit)` - F1 and
  F2 in one place; `angle_unit="deg"` converts a rocking integral to radian.
* `mode_components(mode, ...)` - the named factors a mode applies, and the
  single place that mapping exists. Both integration paths call it.
* `angular_factor(mode, ...)` - `eta`, the product of those and, for the
  rocking modes, the acceptance. They require `detector_acceptance` and raise
  without it (F3); `STATIONARY` *rejects* it.
* `scale_factor(wavelength, unitcell_area, active_area, flux_density)` - the
  prefactor, in 1/s. With the defaults it is the relative scale that makes the
  modes agree; with the measured flux and area it is the absolute one.
* `structure_factor_squared` / `integrated_intensity` - the reduction and its
  exact inverse. The inverse is public because it is the forward model the
  regression tests simulate with.
* `reflectivity_from_structure_factor` / `structure_factor_from_reflectivity` -
  Vlieg eq. 63, generalized to a non-specular exit angle.

Units are mixed on purpose and stated at every boundary: wavelength and
unit-cell area in Angstrom (repository convention, and what `Lattice.uc_area`
returns), areas and flux density in meter (what the config carries). Getting
that wrong is a factor of `1e20`, so
`test_scale_factor_carries_its_documented_units` pins it.

### 4.3 The active area, and a function that should not exist

An earlier draft of this work added
`active_area_footprint(alpha, w, h, L) = w * min(L, h/sin(alpha))`. It is
redundant. `BeamProfile.illuminated_area_fraction` is already the mean of
`p(z)/p_max` over the projected sample, so the absolute area is

```
A(alpha) = beam_width * sample_length * illuminated_area_fraction(alpha, L)
```

and the `min()` form is exactly this for a `top_hat_profile`, and only for
that profile -- verified to `3.3e-16`. Extending the existing one-dimensional
area machinery to a real area is therefore one multiplication by the beam
width, which is what `activearea.beam_limited_area` does. The closed form was
deleted rather than moved, and `test_corrections_activearea.py` pins both the
identity and the fact that a Gaussian of the same width departs from it by up
to 19 %, so it does not come back as a "simplification".

`activearea.slit_limited_area` covers the other limit, built on the same
`geometry.area_correction` the z-axis table lists, so `1/sin(delta)` also has
one definition.

**No existing number changed.** The move is behavior-preserving: the full
suite reports the same 80 failures (all from the unbuilt native extension)
before and after, with 782 passing against 770.

### Test coverage

`test_corrections_measurement.py` writes every expectation out from the papers
rather than calling the module back - the z-axis Lorentz factors, the linear
dependence on `Delta_gamma`, the mode ratio of eq. 65, the unit conversions,
and the round trip through the forward model.

The strongest of them is `test_reflectivity_reproduces_the_fresnel_asymptote`:
for a semi-infinite substrate far from a bulk Bragg peak the truncation-rod
amplitude is `|F| = rho_e A_u / q_z`, and feeding that into eq. 63 must
reproduce the textbook `(q_c / 2 q_z)^4` with `q_c^2 = 16 pi r_e rho_e`. The
two are algebraically identical, so the test validates the absolute scale, the
Angstrom/meter conversions and both powers of `sin` end to end. It agrees to
`rtol = 1e-12`.

`test_scan_mode_equivalence.py` simulates the real thing: a Pt(111) `(1, 0, l)`
rod, angles from `HKLVlieg.VliegAngles.anglesZmode`, measured as rocking scans
(different counting time, different monitor, acceptance varying 0.7x to 1.6x
along the rod, flat background, real `_compute_rocking_integration`
aggregation) and as a stationary *l* scan (different counting time and monitor,
real `integration_corrections` chain). Reduced through the new module, both
recover the input `|F|^2` - the stationary path exactly, the rocking one to
`1e-6`, limited only by the trapezoidal sampling of the rocking profile.

`test_rocking_and_stationary_paths_differ_by_the_missing_normalizations` is a
**characterization** test: it asserts today's gap is exactly
`T * monitor * Delta_gamma_in_degrees`. It must be replaced by a plain equality
when the GUI paths adopt the reduction. Its docstring says so.

### 4.4 The Delta_gamma estimator

`corrections/acceptance.py` supplies the factor F3 says is missing. It reads
the exit angle at the top and bottom **edges** of a region of interest with
`Detector2D_SXRD.surfaceAnglesPoint` and returns the span:

* `out_of_plane_acceptance(detector, row, column, row_size, alpha, ...)` -
  `Delta_gamma` in radian, evaluated at the region's centre column, which is
  where the rod crosses the aperture. Vectorized over frames, since orGUI
  resizes regions along a scan.
* `gamma_range(...)` - the span over the whole rectangle. It equals the
  centre-column value when iso-gamma lines run along the detector rows, and
  exceeds it on a rolled detector, so their ratio is a cheap check on whether
  the estimate can be trusted.
* `pixel_acceptance(...)` - the one-row case; the resolution limit, and the
  natural unit to quote an acceptance in.

Edges, not pixel centres: a region of `n` rows accepts `n` pixels' worth of
angle, not `n - 1`. Using centres would understate every rocking acceptance
by one pixel - 1 % on a 100-row region, 10 % on a 10-row one.

Two things fell out of writing the tests that were not obvious beforehand:

* **pyFAI requires both pixel-coordinate arrays to have the same size.** A
  scalar column with a per-frame `row_size` - exactly what a rocking
  integration passes - failed inside the extension with
  `assert pos2.size == size`. The estimator broadcasts the two coordinates
  before the call.
* **`Delta_gamma` is nearly invariant under a detector-arm rotation.** Driving
  the gamma arm shifts the exit angle at the region centre by the full arm
  angle but leaves the span unchanged to nine digits, because the rotation is
  about the axis gamma is measured around; a 40-degree delta-arm rotation
  moves it by one part in `1e5`. This is the opposite of the polarization
  factor of F5, where ignoring the arm is a 10 % error at 18 degrees. An
  acceptance computed without arm bookkeeping is still usable.

Measured against a calibrated 172 um detector at 1 m, a 60-row region accepts
0.5913 degrees, matching `n * pixel / dist` to 2e-3, and shrinks by a few per
mille towards the detector edge - the obliquity that makes `Delta_gamma` vary
along a scan in the first place.

## 5. Status of each finding

**F1, F2, F3 and F6 are applied.** Steps 1-5 below are done; the wiring of
each is in [`ctr_structure_factor_handover.md`](ctr_structure_factor_handover.md)
section 5, and the equivalence is asserted by
`test_scan_mode_equivalence.py::test_rocking_and_stationary_paths_agree`.
Step 6, the real-data check, is the one that remains, together with F4 and
F7.

Recorded in the order they had to be done, because the ordering was itself a
result -- F6 had to be settled before F3 could be wired:

1. **Normalize rocking scans.** Give `RockingPeakIntegrator.integrate` the same
   `normalization_divisor` the stationary path uses. The per-frame counting
   time, not the sum over the scan: the rocking angle is the integration
   variable, the time is not. (F1)
2. **Integrate in radian**, or divide by `180/pi` at the end. (F2)
3. **Divide by `Delta_gamma`.** (F3) The estimator now exists — section 4.4,
   `corrections/acceptance.py` — so this is a wiring step rather than a design
   one:
   * `out_of_plane_acceptance` takes the span at the ROI centre column, which
     is the resolution of the rotated-ROI question raised when this section was
     written: the corner-to-corner span overestimates the acceptance at the
     `delta` of the rod, so it is reported separately by `gamma_range` as a
     diagnostic rather than used.
   * Pass the per-frame arm angles from `scan_arm_angles` where they are known.
     They turn out to matter far less than expected (section 4.4), but they
     cost nothing.
   * Settle the solid-angle compensation first (step 4): `Delta_gamma` and
     the solid-angle correction each carry the same obliquity, and leaving
     both in double-counts it.
4. **Divide the solid-angle correction back out of `F2_hkl`, in both modes.**
   (F6) A ROI sum is already a complete angular integral, so the correction
   must not reach a structure factor - but it stays on the intensity, where a
   broad or diffuse feature needs it. This also touches `orGUI.integrateROI`,
   and it changes stationary `F2_hkl` by `<1/Omega~>` whenever the switch was
   on - 0.7 % at 1 m, 7 % at 0.3 m - so it belongs in the same breaking commit
   as the rocking wiring, not a later one. The reconstruction path keeps
   applying it uncompensated and is unaffected.
5. **Record the mode and its inputs in the saved data**, next to `F2_hkl`: the
   acceptance, the normalization that was applied, and the active-area
   assumption of F4. Without them a saved rod cannot be put on a common scale
   after the fact.
6. **Then verify on real data. Not done — this is the remaining step.** The
   overlap region of a rocking scan and a stationary scan on the same rod
   (Drnec Fig. 8, right) is the acceptance test. Simulation cannot catch F7,
   an incorrect `Delta_gamma` definition, or a beamline that reports counting
   time in the wrong place. The simulated equivalence test does cover the F6
   compensation
   (`test_the_solid_angle_correction_does_not_reach_the_structure_factor`),
   but it supplies the region mean itself rather than reading a detector, so
   the geometry evaluation is pinned separately in
   `test_corrections_detector.py`.

## 6. What is needed to close issue #15, and reflectivity

Issue #15 asks for
`I_sc = (Phi_0 T r_e^2 A_0 lambda^2 / A_u^2) |F|^2 P L_s C_area C_det C_beam`
with the constants put back. After section 5, what is still missing is only the
constants:

* `lambda`, `A_u`, `r_e` - already available (`UBCalculator.getLambda`,
  `Lattice.uc_area`, `measurement.CLASSICAL_ELECTRON_RADIUS`).
* `Phi_0` and `T` - user input. `T` is in the scan; `Phi_0` needs a flux
  measurement and a field to put it in.
* `A` - `activearea.beam_limited_area(alpha, beam_width, sample_length,
  profile)` returns it in square meter, as `w * L * C_illum_area(alpha, L)`
  (section 4.3; `C_illum_area` alone is a fraction, not an area - F4). Of its
  inputs only the **horizontal beam width** is missing: the footprint dialog
  already asks for the sample size and the beam profile, and the profile
  carries the vertical direction. `activearea.slit_limited_area` covers the
  narrow-slit case instead.
* `C_det` - F7, unmodelled.

**Reflectivity comes for free.**
`R = r_e^2 lambda^2 P_r |F|^2 / (A_u^2 sin(alpha) sin(beta_out))` needs nothing
beyond an absolutely scaled `|F|^2`, so specular reflectivity and the crystal
truncation rods can be refined together against one model.

The data container is already in place for this. `CTRplotutil` gives every
`CTR` a frozen `MeasurementReduction` whose `quantity` is either
`"structure_factor"` or `"reflectivity"` (the dimensionless R), and a
`CTRScanGeometry` recording the z-mode scan rule (`fixed="in"/"out"/"eq"` plus
the fixed angle) - which is exactly `alpha` and `beta_out`. What is missing is
only the conversion, and that is what
`measurement.structure_factor_from_reflectivity` now supplies: given
`Lattice.uc_area` and the wavelength, a reflectivity rod can be turned into an
`|F|^2` rod and joined to the truncation rods. Today every fitting, scaling,
averaging and export path in `CTRplotutil` rejects reflectivity explicitly
(`CTRplotutil.py:1108`, `:1497`), which is the right default while no
conversion existed; the natural next step is a `CTR` method that performs it
rather than relaxing those guards.

Three caveats, all now in the module docstrings:

1. It is the **kinematic** result. It fails near a bulk Bragg peak and below
   the critical angle, where refraction and multiple scattering take over -
   which is where the interesting part of a reflectivity curve usually is.
   Compare against the DWBA machinery already in `CTRdwba.py` there.
2. **A reflectivity scan moves the detector arm**, which is where F5 bit
   hardest: 10 % at `2theta = 18` degrees, 33 % at 30. That is now corrected
   for the two direct-space integration paths, but *not* for the
   reciprocal-space reconstruction, which applies the polarization per pixel
   (see F5).
3. **Off-specular**, `R` is the fraction of the incident flux scattered into
   that rod. That is well defined for a truncation rod integrated across its
   cross-section, but not for diffuse scattering, where only a differential
   cross-section is meaningful. `reflectivity_from_structure_factor` takes a
   `beta_out` for the first case and should not be used for the second.

Also note the ANA "reflectivity rocking scan" Lorentz factor `1/sin(2 alpha)`,
which `peak1Dintegr` already selects for a `mu` scan: that is for rocking
*through* the specular ridge and is a different measurement from a stationary
specular scan. Both are legitimate; they are not interchangeable, and only the
stationary one is what `reflectivity_from_structure_factor` reduces.

## 7. Open questions

Nothing here now blocks the wiring. The two that used to head this list are
settled: `Delta_gamma` for a rotated ROI in section 4.4 (the centre-column
span, with the corner-to-corner span kept as a diagnostic), and the F6
solid-angle double count in F6 itself (ROI-summed integration stops applying
it, in both modes).

* **The stored intensity column changes too, not only `F2_hkl`.** Dropping
  the solid-angle correction moves it by `<1/Omega~>`, so external
  post-processing scripts reading the intensity column - not just those
  reading `F2_hkl` - see the break. It is the same factor either way; there
  is no variant of F6 that leaves the intensity column alone and still puts
  the modes on one scale.
* **The reciprocal-space reconstruction as a third route.** Drnec section 4
  shows that voxel binning absorbs the Lorentz factor, so a reconstructed map
  needs `C_area`, `C_beam`, `P` and `Delta_l` but no Lorentz factor - which is
  what `reconstruction_job.py` does, minus `Delta_l` and `C_beam`. Landing it
  on the same absolute scale as the two direct-space modes is a third, larger
  piece of work and was not analysed here.
* **`C_det` (F7)** is the only remaining mechanism that can break the mode
  equivalence after F1-F3, and it is the one that cannot be validated on
  simulated data. It needs the real-data comparison of section 5, step 6.
* **Error propagation through the new factors.** `measurement` reduces
  intensities; the errors follow the same divisors, but a `Delta_gamma`
  estimated from the geometry has an uncertainty of its own that nothing
  currently tracks.

## 8. Reproducing the numbers

```powershell
pytest orgui/datautils/xrayutils/test/test_corrections_measurement.py
pytest orgui/datautils/xrayutils/test/test_corrections_acceptance.py
pytest orgui/app/test/test_scan_mode_equivalence.py
```

The F6 cancellation and the table in it are
`test_corrections_acceptance.py::test_the_solid_angle_correction_and_the_acceptance_carry_one_obliquity`.

None of these needs the native extension. Much of the rest of the suite does,
and a checkout without it reports around 80 failures that have nothing to do
with this work — `ctr_structure_factor_handover.md` section 3 has the
out-of-tree build recipe that turns those into a clean
1120 passed / 3 skipped run.

## 9. Where the implementation stands

This document is the physics. For the state of the branch — what is wired and
what is not, the next commit in detail, the build recipe, and the predictions
made here that turned out wrong — see
[`ctr_structure_factor_handover.md`](ctr_structure_factor_handover.md).

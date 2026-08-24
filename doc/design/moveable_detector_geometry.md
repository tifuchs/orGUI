# Moveable detector geometry: pyFAI goniometer vs. the Vlieg 6-circle

Status: assessment; Sec. 7 steps A (detector arm), B (Q path) and the
half-pixel fix are implemented. Reconstruction (D), scan/config plumbing for
per-frame arm angles (part of C), the GUI (E) and the pyFAI
`GoniometerRefinement` wrapper (F) are not.
Scope: rotating detector arm (`gamma_arm`, `delta_arm`) now, translation later.

## 1. Verdict up front

Adopting the pyFAI *geometry* parametrisation for a moving detector does **not**
break orGUI's 6-circle equations. The reason is a stronger statement than
"approximately compatible", and it was verified numerically rather than assumed:

> A rigid rotation of the detector about the sample position is represented
> **exactly** by pyFAI's `(rot1, rot2, rot3)`, with `dist`, `poni1`, `poni2`
> left untouched.

Because pyFAI places the sample at the origin and computes a pixel position as
`R(rot1,rot2,rot3) @ [p1 - poni1, p2 - poni2, dist]`, rotating the detector
about the sample is a left-multiplication `R_arm @ R`, and `R_arm @ R` is again
a rotation that the `rot3*rot2*rot1` triple can express. The sample-to-PONI
distance is preserved, so `dist` does not change, and `poni1/poni2` are
detector-internal so they do not change either. Verified to machine precision
(~2e-16 in the outgoing unit ray) for combinations up to `gamma_arm = 40 deg`,
`delta_arm = 75 deg`, on top of a home geometry with non-zero tilts.

What is *not* safe, and is the main thing this document argues for, is the
naive shortcut `gamma = gamma_arm + gamma_pixel`. Rotations do not add.
Measured error of the additive assumption, Pilatus 1M at 0.8 m:

| arm position | max error in gamma | max error in delta |
| --- | --- | --- |
| `gamma_arm=10 deg, delta_arm=20 deg` | 0.060 deg | 0.215 deg |
| `gamma_arm=30 deg, delta_arm=60 deg` | 0.204 deg | 1.414 deg |

That is one to two orders of magnitude above the per-pixel angular resolution.
The `gamma_arm` / `delta_arm` split is therefore not cosmetic bookkeeping — it
is required for correctness, and the conversion between arm and measured angles
must go through matrices, never through angle addition.

`pyFAI.goniometer` itself is a *fitting and bookkeeping* layer, not a geometry
model (Sec. 3). Its value to orGUI is `GoniometerRefinement` — calibrating one
shared geometry against several images taken at different arm positions. Its
`GeometryTransformation` formula mechanism is usable but awkward for a
two-circle arm; `BaseTransformation` (a plain Python callable) is the better
hook.

## 2. What orGUI assumes today

The static-detector assumption is not localised in one place; it is a contract
running through four layers.

**Layer 1 — pixel to lab angles.** `Detector2D_SXRD`
(`orgui/datautils/xrayutils/DetectorCalibration.py:77`) subclasses
`pyFAI.geometry.Geometry` and adds the SXRD angle conventions on top:

- `primBeamAngles` (`:251`) / `primBeamPoints` (`:273`) turn pyFAI's
  `(tth, chi)` into lab-frame `(gamma_p, delta_p)`.
- `surfaceAngles` (`:283`) / `surfaceAnglesPoint` (`:328`) tilt those into the
  alpha (surface) frame, giving the `(gamma, delta)` that the Vlieg equations
  consume.
- `pixelsTthChi` (`:386`), `pixelsPrimeBeam` (`:548`) and
  `pixelsSurfaceAngles` (`:563`) are the analytic inverses.

The whole class is one geometry with one cached array set. `_cached_array`,
`_alpha_i` and `_n_refr` memoise per-pixel angle arrays under the assumption
that geometry only changes when the user edits calibration.

**Layer 2 — the frame/Q conversion.** `orgui/app/qconversion.py` is the newer,
faster path and never forms `(gamma, delta)` at all: `qComponents` (`:438`)
calls `detectorCal.calc_pos_zyx(...)` (`:748`) and contracts the resulting
position with a precomputed `(G, c)` pair from `conversionCoefficients`
(`:225`). This layer is the *easiest* to make moveable, because a detector
rotation is exactly a rotation of the position vectors it already works with.
`detectorAngles` (`:596`) is the inverse (Q back to `delta, gamma`).

**Layer 3 — reciprocal-space reconstruction.** This is where the static
assumption is hardest-coded. `_detector_corner_rays`
(`orgui/datautils/xrayutils/reconstruction.py:271`) computes lab-frame unit
rays at pixel corners **once per mapping run**, `_tile_ray_arrays` (`:1517`)
caches them per detector tile, and `_map_frame_group` (`:1565`) shares that
read-only array across every compute worker for every frame. The native kernel
receives per-frame `angles_start` / `angles_end` of shape `(frames, 4)` —
`(alpha, omega, chi, phi)`, sample circles only. In
`orgui/datautils/xrayutils/cpp/reciprocal_reconstruction_cpp.cpp`,
`coordinate_at` (`:1773`) forms `Q = K * (ray - y_hat)` and then
`apply_frame_rotation` (`:1730`) applies the sample rotations. The detector's
contribution is entirely inside the frozen `ray`.

**Layer 4 — scan and config.** `scan_exposure_angle_bounds`
(`orgui/backend/scans.py:298`) builds the `(frames, 2, 4)` angle bounds from
`scan.mu` / `scan.omega` plus fixed `chi`/`phi` from config. The `Scan` base
class (`:83`) documents only `th`/`om`/`mu` as required. Detector motors have
no place to live.

## 3. What `pyFAI.goniometer` actually is

Worth being precise, because the name oversells it.

`Goniometer` (`pyFAI/goniometer.py:395`) holds a detector, a wavelength, a
parameter vector `param`, and a `trans_function`. The entire "geometry" is:

```
trans_function(param, motor_positions) -> (dist, poni1, poni2, rot1, rot2, rot3)
```

`get_ai(position)` (`:475`) calls it and builds a fresh `AzimuthalIntegrator`.
That is all. There is **no** kinematic model, no axis directions, no
rotation-centre concept, no notion of which motor is which. pyFAI does not know
your arm has two circles; it knows a formula that produces PONI parameters.

Three flavours of `trans_function`:

- `GeometryTransformation` (`:132`) — six `numexpr` strings, one per PONI
  parameter, over named parameters and named motor positions. Serialisable to
  JSON. This is what the pyFAI goniometer tutorials use, typically with trivial
  formulas like `rot2 = pos * pi / 180`.
- `ExtendedTransformation` (`:275`) — same, plus wavelength as a fitted
  quantity.
- `BaseTransformation` (`:72`) — an arbitrary Python callable. Not
  serialisable (`to_dict` raises), but unrestricted.

So the answer to "which rotations and translations are allowed" is: **all of
them, and none of them specifically.** The allowed set is whatever the
underlying PONI parametrisation can express, which is:

- **Rotations:** the full `SO(3)` of detector orientation, via
  `R = rot3 * rot2 * rot1` (`pyFAI/geometry/core.py:2656`), where `rot1` is a
  left-handed rotation about lab axis 1, `rot2` left-handed about axis 2, and
  `rot3` right-handed about axis 3 (the beam). Equivalently
  `R = Rz(rot3) * Ry(-rot2) * Rx(-rot1)`, an intrinsic Z-Y-X Euler sequence,
  with the usual gimbal degeneracy at `rot2 = +/- 90 deg`.
- **Translations:** only along the detector normal, through `dist`, plus the
  two in-plane offsets `poni1`, `poni2`. There is **no free translation of the
  detector in the lab frame.** This is the real limitation and it is what makes
  future translation support a genuine (if modest) extension rather than a
  parameter change — see Sec. 8.

`GoniometerRefinement` (`:777`) is the part that motivates the whole exercise:
`SingleGeometry` (`:602`) pairs an image with a motor position, `extract_cp`
(`:675`) picks calibrant rings, and `refine2` / `refine3` (`:910`, `:962`) fit
the shared `param` vector across all positions at once. That is exactly the
"calibration tools already exist" benefit — and note it is script/API-level;
`pyFAI-calib2` does not drive it.

## 4. The bridge between the two conventions

Both conventions had to be pinned down before compatibility could be checked.

**orGUI / Vlieg lab frame** (from `calculate_q_phi`,
`orgui/datautils/xrayutils/HKLVlieg.py:40`, and `QAlpha`, `:1362`): beam along
`+y`, surface normal along `+z`, in-plane horizontal `+x`. The outgoing wave
vector is

```
k_f / K = (sin(delta) cos(gamma), cos(delta) cos(gamma), sin(gamma))
        = DELTA(delta) @ GAMMA(gamma) @ y_hat
```

with `DELTA = Rz(-delta)` and `GAMMA = Rx(gamma)` (`:70`, `:76`). Note the
composition order: `GAMMA` acts first, so **`delta` is the base circle about
the vertical `z`, and the `gamma` circle rides on the delta arm.**

**pyFAI lab frame:** `calc_pos_zyx` returns `(t3, t1, t2)` = (along beam, slow,
fast).

**The relabelling**, with `dchi` the orGUI azimuthal reference
(`setAzimuthalReference`, `DetectorCalibration.py:236`; `azimuthal_reference`
in the config files):

```
x_Vlieg =  cos(dchi) t1 + sin(dchi) t2
y_Vlieg =  t3
z_Vlieg = -sin(dchi) t1 + cos(dchi) t2
```

Verified against `primBeamPoints` to 1.1e-16 for `dchi` = 0, 90 deg and an
arbitrary 21.2 deg, with non-zero `rot1/rot2/rot3`. For the standard SXRD
setting `azimuthal_reference = 90` (as in `examples/config_minimal`) this is
simply `x = t2`, `y = t3`, `z = -t1`.

## 5. Compatibility check (the part that was asked for)

With `M` the relabelling matrix above, `R_0` the pyFAI rotation of the home
calibration, and `(gamma_0, delta_0)` the arm reading **at which that
calibration was taken** (Sec. 6.1), a detector arm at
`(gamma_arm, delta_arm)` gives

```
R_arm(g, d) = DELTA(d) * GAMMA(g)
R_new = M^T * R_arm(gamma_arm, delta_arm) * R_arm(gamma_0, delta_0)^T * M * R_0
(rot1, rot2, rot3) = euler_ZYX_decompose(R_new)
dist, poni1, poni2 unchanged
```

At `(gamma_arm, delta_arm) = (gamma_0, delta_0)` this collapses to `R_0`, so a
detector parked at its calibration position behaves byte-identically to today.

Feeding that back through `Detector2D_SXRD.primBeamPoints` reproduces the
directly rotated rays to **2.5e-16**. The 6-circle equations are untouched:
they still receive `(gamma, delta)` per pixel with the same meaning as today.

For the special case `dchi = 90 deg` and an untilted home geometry the
decomposition has a closed form (verified to 3.6e-16 over `gamma_arm` in
+/-60 deg, `delta_arm` in +/-85 deg):

```
rot1 = atan2(sin(delta_arm), cos(delta_arm) cos(gamma_arm))
rot2 = asin(sin(gamma_arm) cos(delta_arm))
rot3 = atan2(sin(delta_arm) sin(gamma_arm), cos(gamma_arm))
```

These are `numexpr`-expressible, so they drop straight into a pyFAI
`GeometryTransformation` if that route is preferred. Three observations:

- `rot3` is **not** zero. A `delta`-base / `gamma`-riding arm induces a detector
  roll about the beam when expressed in pyFAI's fixed axis triple — 26.6 deg at
  `gamma_arm=30, delta_arm=60`. Anyone tempted to write `rot1 = delta;
  rot2 = gamma; rot3 = 0` into a goniometer file will be wrong by that much.
- That naive assignment *is* exact for the **opposite** mounting order
  (`GAMMA @ DELTA`, i.e. gamma as the base circle) — verified exactly. So the
  correct formula depends on the physical stacking of the two circles on the
  real machine. This is a per-beamline fact that must be captured in
  configuration, not hard-coded.
- Gimbal lock sits at `rot2 = +/- 90 deg`, which corresponds to
  `gamma_arm = +/- 90 deg` — out of reach for any realistic arm, but the
  decomposition should still use a matrix routine that fails loudly rather than
  silently.

**Because the closed form only holds for the untilted, `dchi = 90` case, the
implementation should compose and decompose matrices, not evaluate formulas.**
A real home geometry has non-zero `rot1/rot2/rot3` from calibration, and those
do not factor out.

## 6. `gamma_arm` / `delta_arm`: where the split has to live

Making the distinction precise:

- `gamma_arm`, `delta_arm` — **mechanical**, one pair per frame, read from the
  motors. They orient the detector body.
- `gamma`, `delta` — **measured**, one pair per *pixel* per frame. They are
  what the Vlieg equations consume, and today they come out of
  `surfaceAnglesPoint`.

Today the two coincide only in the degenerate sense that the arm sits at a
fixed position folded into the calibration.

### 6.1 What the arm angles are referenced against

**Not the beam, and not the central pixel. The reference is the arm position at
which the PONI calibration was taken.**

`gamma_arm` and `delta_arm` are motor readouts. Their zero is whatever the
beamline homed them to, and nothing ties that zero to the optical axis. Even if
the arm genuinely were at a true mechanical zero, the direct beam would still
not sit at the PONI unless `rot1 = rot2 = 0`. The PONI is the foot of the
perpendicular from the sample to the detector plane; the beam centre is where
the undeflected beam lands. Those are four different points that get casually
conflated (Pilatus 1M, 0.8 m, `poni1 = poni2 = 0.0866 m`):

| home tilts | PONI pixel | beam-centre pixel | separation | `(gamma, delta)` at the PONI |
| --- | --- | --- | --- | --- |
| `rot1=rot2=rot3=0` | (503.5, 503.5) | (503.5, 503.5) | 0 px | (-0.006, +0.006) deg |
| `rot1=0.02, rot2=-0.015, rot3=0.05` | (503.5, 503.5) | (433.7, 410.5) | 116 px | (-0.807, +1.194) deg |
| `rot1=0.15, rot2=-0.10, rot3=0.30` | (503.5, 503.5) | (31.5, -199.5) | 847 px | (-2.878, +9.921) deg |

The detector's own geometric centre — (521.5, 490.5) here — is unrelated to
either. So with entirely ordinary calibration tilts the PONI pixel already
reads about a degree of `delta`, and the beam centre has walked 116 px away.
"Zero at the central pixel" is not a property the system has, at any arm
position.

The workable convention is therefore relative: store `(gamma_0, delta_0)`, the
arm reading at calibration time, next to the PONI, and use the Sec. 5 formula.
Consequences:

- The arm's absolute zero never has to mean anything, and the beam never has to
  be on the detector.
- Any existing calibration keeps working: an absent `(gamma_0, delta_0)`
  defaults to the current arm reading, giving identity.
- The relative rotation is `R_arm(g,d) * R_arm(g_0,d_0)^T`, **not**
  `R_arm(g - g_0, d - d_0)`. Rotations do not subtract any more than they add:
  the discrepancy is 0.26 deg for a `(0.5, 0.5) deg` home arm position and
  2.76 deg for `(2, 5) deg`.

#### The degeneracy that decides what is worth fitting

A motor-zero offset and a detector *mounting* misalignment are not independent
parameters, and which of them is recoverable follows from the composition
order. Write the detector body orientation as `R_arm * R_mount`, with
`R_mount` fixed (detector relative to the arm). For orGUI's `DELTA @ GAMMA`
ordering:

- A **riding-circle** (`gamma`) zero offset is **exactly** absorbable into
  `R_mount`, since `D(d) G(g + g_0) R_mount = D(d) G(g) [G(g_0) R_mount]`.
  Verified to 2.2e-16. It is structurally unrecoverable from ring or reflection
  calibration at *any* number of arm positions. Do not add it as a fit
  parameter — it is a gauge choice, already absorbed in the home PONI.
- A **base-circle** (`delta`) zero offset sits between the two circles,
  `D(d) D(d_0) G(g) R_mount`, and cannot be pushed past `G(g)`. It is therefore
  recoverable — but only from data in which `gamma` actually moves:

| `delta_0` error | gamma span | residual orientation error |
| --- | --- | --- |
| 0.1 deg | 0 deg | 0.000 deg |
| 0.1 deg | 30 deg | 0.052 deg |
| 0.5 deg | 10 deg | 0.087 deg |
| 0.5 deg | 30 deg | 0.259 deg |

With `gamma` fixed the residual is identically zero: at a single `gamma`, the
base-circle zero is just as degenerate as the riding-circle one.

If the machine stacks the circles the other way round (`GAMMA @ DELTA`, gamma
as the base), the two roles swap. This is the same per-beamline fact as
Sec. 9 question 1, and it decides the fit parameter list, so it is worth
settling before writing the refinement wrapper.

Practical upshot for calibration: one calibrant image at one arm position
cannot separate arm zero from mount tilt — and does not need to, because the
combination is all that enters the prediction. A multi-position
`GoniometerRefinement` (Sec. 7F) buys exactly one extra number, the base-circle
zero, and only if the ring images span several riding-circle positions.

### 6.2 Which frame the arm angles live in

This is the subtlety behind the "if the motors are independent of alpha"
caveat, and it matters more than it looks.

orGUI's `(gamma, delta)` are defined in the **alpha frame** — `calculate_q_phi`
computes `k_f - k_i` with `k_i = ALPHA^T y_hat`, so the surface tilt is already
taken out. The detector arm, by contrast, is usually bolted to the floor and
therefore rotates in the **lab frame**, independent of alpha. The existing
helpers `primBeamAngles` (`HKLVlieg.py:112`) and `vliegDiffracAngles` (`:121`)
already convert single directions between these two frames.

They are not sufficient for arm angles. Converting a *frame* is not converting
a *direction*: if the arm is mounted on the alpha circle, its lab-frame
rotation is `ALPHA * DELTA * GAMMA * ALPHA^T`, and forcing that back into a
`DELTA(delta_p) * GAMMA(gamma_p)` pair leaves a residual roll about the
outgoing beam that a direction-only conversion silently discards:

| alpha | arm `(gamma, delta)` | residual roll |
| --- | --- | --- |
| 0.1 deg | (10, 30) deg | 0.051 deg |
| 1.0 deg | (30, 60) deg | 0.995 deg |
| 5.0 deg | (30, 60) deg | 4.878 deg |
| 15.0 deg | (30, 60) deg | 13.99 deg |

At true grazing incidence (alpha ~ 0.1-1 deg) this is 0.05-1 deg of roll —
small, but at a 5 deg detector half-diagonal it still displaces edge pixels by
~0.09 deg, comfortably above pixel resolution. **Compose matrices; do not
convert arm angles through `primBeamAngles`.**

Recommendation: define the arm rotation in the **lab frame** and apply it
between `primBeamPoints` and `surfaceAngles`. That keeps the alpha conversion,
the refraction correction and every downstream Vlieg call byte-identical to
today, and it matches how floor-mounted arms actually behave. Machines whose
detector arm rides on the alpha circle get a config flag that pre-multiplies by
`ALPHA * R_arm * ALPHA^T`.

### 6.3 API shape

**Implemented** (`Detector2D_SXRD`, camelCase to match `setAzimuthalReference`
and `getPolarization`):

```python
Detector2D_SXRD.armRotation(gamma_arm, delta_arm)   # -> 3x3, DELTA @ GAMMA
det.setArmReference(gamma_arm=..., delta_arm=...)   # once, with the PONI
det.setArm(gamma_arm=..., delta_arm=...)            # per measurement, radians
det.setArm(R_arm)                                   # or a raw rotation matrix
det.armParam                                        # arm-adjusted pyFAI geometry
```

Both `_R_arm` and `_R_arm_reference` default to the **identity**, i.e.
`(gamma_0, delta_0) = (0, 0)`, so an untouched detector is the static geometry.
Every existing method keeps its signature and its meaning, and takes its
original code path verbatim while `armIsIdentity()` holds; a detector at its
calibrated arm position produces bit-identical arrays. This is the property
that keeps the change from being a rewrite, and it is asserted by
`test_identity_arm_leaves_every_routed_result_untouched`.

## 7. Code change inventory

Ordered by increasing pain. Sizes are for the rotation-only case.

**A. `Detector2D_SXRD` — the arm itself.** *(small, self-contained)*
Add arm state plus `R_arm`, and an `effective_param()` returning
`(dist, poni1, poni2, rot1', rot2', rot3')`. Route `primBeamAngles`,
`primBeamPoints`, `surfaceAnglesPointParam`, `pixelsTthChi` and the
`rangegamdel*` properties through it. `surfaceAnglesPointParam`
(`DetectorCalibration.py:339`) already takes an explicit trial `param`, so it
is the template for how the rest should look.

Two traps here:

- pyFAI's `Geometry.reset()` clears `_cached_array` **and calls `gc.collect()`**
  (`pyFAI/geometry/core.py`, `reset`). Setting `rot1` per frame would trigger a
  full garbage collection per frame. The arm must **not** be implemented by
  assigning to the pyFAI properties. Keep `R_arm` separate and pass an explicit
  `param` to `calc_pos_zyx`.
- `_alpha_i` and `_n_refr` are orGUI attributes that `reset()` does not clear.
  They are currently safe only because the guards also test `_cached_array`.
  Any new arm state must be added to those guards deliberately.

**B. `qconversion.py` — the frame/Q path.** *(small)*
`qComponents` (`:438`) already passes an explicit `param` array to
`calc_pos_zyx` (`:748`); it just needs the arm-corrected one. `detectorAngles`
(`:596`) needs the inverse arm rotation applied to the direction before
extracting `(delta, gamma)`. `_IntegratorCache` (`:351`) keys on the geometry
tuple, so it stays correct automatically — but it caches exactly one
integrator, so a moving detector will thrash it. Make it a small LRU if
`integrateImage` is ever used per frame.

**C. Scan and config plumbing.** *(medium, mostly mechanical, wide)*
`Scan` gains optional `gamma_arm` / `delta_arm` arrays and
`gamma_arm_bounds_rad` / `delta_arm_bounds_rad`, defaulting to the config's
fixed values when a backend does not supply them — mirroring exactly how
`chi`/`phi` are handled now. `scan_exposure_angle_bounds`
(`orgui/backend/scans.py:298`) widens from 4 to 6 angles. Config gains a
`[Detector]` section for home arm position and the circle stacking order.
`ScanReference` (`:359`) and the NX serialisation in
`DetectorCalibration.toNXdict` (`:83`) need the new fields, with defaults so
old configs and stored calibrations keep loading — this is the backwards
compatibility requirement from the root `AGENTS.md`, and it is easy to satisfy
since an absent arm means identity.

**D. Reconstruction.** *(the real work)*
The frozen `ray_arrays` are the blocker. Three options:

1. *Rebuild rays per frame.* Simple, correct, and unusably slow: a full
   pixel-corner pass through pyFAI per frame, and it destroys the "build once,
   share read-only across workers" property that `_map_frame_group` documents.
   Rejected.
2. *Rotate rays in Python per frame.* One `(rows+1, columns+1, 3)` matrix
   multiply per frame per tile. Cheap relative to the kernel call, and it needs
   no C++ change at all. Correct for pure rotation. **This is the right first
   implementation** — it is also the natural fallback when the native module is
   unavailable.
3. *Push the arm into the kernel.* Widen `angles_start`/`angles_end` from 4 to
   6 and apply `R_det(t)` to the ray inside `coordinate_at`
   (`reciprocal_reconstruction_cpp.cpp:1773`), before `Q = K(ray - y_hat)`. The
   temporal machinery already exists — `frame_rotation` (`:1687`) interpolates
   linearly in `t`, and the kernel already subdivides 8-way rather than 4-way
   for a moving exposure (`_admissible_sample_pixels`, `reconstruction.py:376`).
   So a detector that moves *during* an exposure costs nothing new
   architecturally. This is the correct end state, but it is a native-code
   change with a rebuild and a matching Python fallback to keep numerically
   consistent (see `orgui/datautils/xrayutils/AGENTS.md`, "Hotspots").

Option 2 first, option 3 when continuous-arm-motion scans are actually needed.
Option 2 changes only *which* rays go in, so checkpoint format, routing, and
resume logic are untouched.

Also in this layer: the automatic reciprocal-space volume selection added in
`2368fdd`/`c43eec2` estimates grid extents from the detector's angular range.
`rangegamdel` (`DetectorCalibration.py:524`) becomes per-frame, so the estimate
must take the union over the arm trajectory, not the home position.

**E. GUI.** *(medium)*
`QUBCalculator` geometry fitting (`:1025-1140`) fits the six PONI parameters
against reflections; with an arm it must fit the *home* geometry while holding
each reflection's arm position fixed. `surfaceAnglesPointParam` already accepts
a trial `param`, so the change is to compose the per-reflection arm rotation
into the trial before calling it — genuinely small. Reflections stored in the
database must record their arm position, or old reflections become
uninterpretable. The `GeometryTabs` widget (`:67`) and the PONI import/export
(`:1179`, `:2263`) need to make clear whether they show home or effective
geometry; exporting *effective* geometry to a PONI file would produce a file
that is silently wrong at any other arm position.

**F. Calibration via pyFAI.** *(optional, and the actual payoff)*
Wrap `GoniometerRefinement` with a `BaseTransformation` whose callable does the
matrix composition of Sec. 5. `BaseTransformation.to_dict` raises, so orGUI
serialises the home geometry plus the stacking order itself rather than using
pyFAI's `.json` goniometer file. If a portable pyFAI goniometer file is wanted,
the closed forms in Sec. 5 give a `GeometryTransformation` equivalent, valid
under its stated restrictions.

## 8. Translation (future)

pyFAI's PONI parametrisation cannot express a lab-frame detector translation:
`dist` moves the detector along its own normal only. Supporting a translating
detector needs a translation vector carried alongside the PONI parameters, with
pixel position `p_lab = R_arm (R_0 p_det) + t_arm`.

The encouraging part is that layers B and D need less than expected:

- `qconversion.py` normalises positions itself; adding `t_arm` before
  normalisation is a two-line change.
- The kernel's `ray_at` (`reciprocal_reconstruction_cpp.cpp:1667`)
  **already re-normalises** after bilinear interpolation. So the corner arrays
  can hold *positions* instead of unit rays with no kernel-math change beyond
  the interpolation being of positions. This is strictly *more* accurate than
  today: bilinear interpolation of corner positions is exact for a flat
  detector, whereas bilinear interpolation of unit rays is an approximation
  that the renormalisation only partly repairs.

Switching `_detector_corner_rays` (`reconstruction.py:271`) to emit positions
rather than unit rays is therefore worth doing **as part of the rotation work**,
not deferred: it costs nothing, improves accuracy slightly, and removes the
main obstacle to translation later. It does change numerical output at the
1e-6-ish level, so it needs its own regression test and a `fix(phys)` commit
separate from the arm feature.

What translation *does* break irreparably is the assumption that the sample
sits at the rotation centre. `(gamma, delta)` per pixel remains well defined (it
is still a direction), but `GeometryCorrection.lorentzFactor` and friends
(`HKLVlieg.py:1773`) and the solid-angle term in `correctionArray`
(`DetectorCalibration.py:610`) need review, because pixel solid angle then
depends on the translation.

## 9. Open questions for the beamline

1. **Circle stacking order.** Is `delta` the base circle with `gamma` riding on
   it (matching orGUI's `DELTA @ GAMMA`), or the reverse? Sec. 5 shows the two
   differ by up to 26 deg in `rot3`. This must be answered per beamline before
   any code is written.
2. **Is the arm on the alpha circle?** Determines lab-frame vs alpha-frame arm
   composition (Sec. 6.2).
3. **Do the arm axes actually intersect at the sample?** The whole
   "`dist` unchanged" result assumes they do. A sphere-of-confusion offset
   shows up as an apparent translation and belongs in the Sec. 8 work.
4. **Are the motor readouts already `gamma`/`delta`, or raw motors?** Where the
   conversion lives (backend vs. core) should follow whatever `orgui/backend/`
   already does for sample angles.

## 10. Suggested sequencing

1. Positions-instead-of-rays in `_detector_corner_rays` (+ regression test).
   Independent, small, unlocks Sec. 8.
2. `R_arm` in `Detector2D_SXRD` with identity default, plus round-trip tests
   against `pixelsSurfaceAngles` at non-zero arm positions. Nothing else
   changes yet — the whole app still runs at the default arm.
3. `qconversion.py` arm-aware. Extends `test_q_conversion.py` naturally.
4. Scan/config plumbing, 4 -> 6 angles, with old-config defaults.
5. Reconstruction option 2 (Python-side ray rotation).
6. GUI: reflection storage, fitting, PONI export semantics.
7. pyFAI `GoniometerRefinement` wrapper.
8. Kernel option 3, only if continuous arm motion during exposure is needed.

Steps 1-3 are individually testable against the existing regression suites
(`test_DetectorCalibration.py`, `test_q_conversion.py`) and none of them change
behaviour at the default arm position — which is the property that makes this
safe to land incrementally.

## Appendix: reproducing the numbers

Every figure quoted above came from short scripts against the installed pyFAI
(2026.5.0) and the current `Detector2D_SXRD`. The essential checks, worth
turning into tests when the work starts:

1. **Axis relabelling** — build a detector with non-zero tilts, compare
   `primBeamPoints` against `calc_pos_zyx` mapped through `M(dchi)`. Expect
   ~1e-16.
2. **Exact representability** — rotate home-position rays by
   `DELTA @ GAMMA`, and separately build a detector from the recomposed
   `(rot1, rot2, rot3)`; compare rays. Expect ~2e-16.
3. **Closed forms** — sweep `gamma_arm`, `delta_arm` and compare the Sec. 5
   expressions against a matrix Euler decomposition. Expect ~4e-16.
4. **Additivity error** — the Sec. 1 table; the check that justifies the whole
   `_arm` split.
5. **Residual roll** — the Sec. 6.2 table; the check that justifies composing
   matrices instead of reusing `primBeamAngles`.
6. **Arm-offset degeneracy** — the Sec. 6.1 tables; the check that decides
   which motor zeros are worth carrying as fit parameters.

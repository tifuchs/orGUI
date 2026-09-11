# CTR structure-factor scale: implementation status and handover

> **Status as of 2026-09-11.** Branch `claude/ctr-structure-factor-9633bc`,
> nothing pushed. `git log --oneline master..HEAD` is the commit list; a
> count written here goes stale on the commit that writes it.
>
> The physics analysis is complete and quantified, and the reduction is now
> **wired in**: a rocking scan and a stationary scan of the same rod come out
> with the same `F2_hkl`, asserted by
> `test_scan_mode_equivalence.py::test_rocking_and_stationary_paths_agree`.
> This **changed saved numbers in both modes**. The **real-data check has now
> been done** — LaNiO3 scan 61, both modes run from the raw images, agreeing
> to a median 1.033 along the CTR; see
> [`ctr_structure_factor_scale.md`](ctr_structure_factor_scale.md) section
> 5.1. [#82](https://github.com/tifuchs/orGUI/issues/82) is closed up to
> `C_det` (F7), which that check bounds but does not model.
> [#15](https://github.com/tifuchs/orGUI/issues/15) still needs the two
> absolute-scale inputs of section 6.
>
> This document is the handover: what exists, how to run it, what to do next,
> and which of my predictions turned out wrong. The physics itself is in
> [`ctr_structure_factor_scale.md`](ctr_structure_factor_scale.md) — read that
> first for *why*; this one is *where things stand*.

## 1. The one-paragraph summary

A rocking scan and a stationary scan of the same rod used to differ by exactly
`T_omega * monitor_omega * Delta_gamma_in_degrees`, measured to seven digits on
simulated data. Three things were missing from the rocking path:
exposure/monitor normalization, integration in radian rather than degrees, and
division by the out-of-plane detector acceptance. A fourth, F6, affected
**both** modes: a ROI sum is already a complete angular integral, so the
solid-angle correction had to be divided back out of `F2_hkl`. It stays
applied to the *intensity*, where a broad or diffuse feature needs it. All
four are now handled, and the two modes agree to `1e-6` on simulated data —
the residual is the trapezoidal sampling of the rocking profile, not a
correction factor.

F5, the arm-blind polarization, is fixed too, in its own commit. It is
independent of mode equivalence — it cancels between the modes at the same
reflection — but it is up to a 33 % error on a scan that drives the detector
arm, which is exactly the reflectivity case.

## 2. What is on the branch

| commit | what it did |
|---|---|
| `777d887` | `refactor: collect correction factors into one corrections package` |
| `4193c80` | `feat: reduce integrated intensities to \|F_hkl\|^2 on one scale` |
| `e25b8df` | `docs: record the rocking/stationary structure-factor scale analysis` |
| `6959191` | `feat: estimate the out-of-plane detector acceptance` |
| `ecb3bb5` | `docs: settle F6 and add the structure-factor physics reference` |
| `0af9cb7` | `feat(phys)!: put rocking and stationary integration on one structure-factor scale` |
| *(this one)* | `fix(phys)!: follow the detector arm in the polarization correction` |

The first three were split out of one working tree at the end, which required
*staged versions* of four files: `peak1Dintegr.py`, `integration_corrections.py`,
`corrections/__init__.py` and `test_corrections_package.py` all reference
`measurement`, which only exists from `4193c80`. In `777d887` they keep the old
inline Lorentz selection with only the import path updated. Each commit was
verified green on its own, not just the tip.

### 2.1 The package

```
orgui/datautils/xrayutils/corrections/
  geometry.py       219   z-axis Lorentz x3, rod interception, area factor
  beamprofile.py   1049   beam profile shapes and their integrals
  activearea.py     190   active area in m^2, slit- and beam-limited
  acceptance.py     255   Delta_gamma, gamma_range, pixel_acceptance
  detector.py       319   per-pixel solid angle and polarization, their
                          region means, and the polarization arm correction
  normalization.py  103   counting time and monitor, from values
  roi.py             93   CorrectionFactors, roi_mean_correction
  measurement.py    605   mode dispatch, master equation, reflectivity
```

`orgui/datautils/xrayutils/geometrycorrections.py` and `beamprofile.py` remain
as re-export aliases because both were released under those names.
`test_corrections_package.py` asserts the aliases hand out the *same objects*,
not merely that they import.

### 2.2 What is wired

| caller | uses the package for | still does its own thing |
|---|---|---|
| `orGUI.integrateROI` (stationary) | `pixel_factors`, `mode_components`, `normalization_divisor`, `C_illum_area`, `roi_mean_inverse_solid_angle`, `polarization_arm_correction` | — |
| `orGUI.rocking_integrate` | `pixel_factors`, `polarization_arm_correction` | — |
| `peak1Dintegr.integrate` (rocking) | `mode_components`, `normalization_divisor`, `normalized_intensity`, `out_of_plane_acceptance`, `roi_mean_inverse_solid_angle` | — |
| `reconstruction_job` | `pixel_factors` (solid angle applied and **not** compensated; polarization **not** arm-corrected) | own native-fused application, own normalization loop |

Still uncalled outside the tests: `measurement.structure_factor_squared`,
`measurement.angular_factor` and `activearea.*`. That is deliberate rather
than a gap — the two integration paths form `F2_hkl` on a *relative* scale by
dividing out the mode-dependent factors they already hold as
interval-weighted means, and `structure_factor_squared` additionally divides
by the absolute prefactor, which needs the issue #15 inputs. `angular_factor`
rebuilds `eta` from point angles, which is the wrong thing for a path that
has ROI-weighted means of each component.

The rocking path could not use `normalization_divisor(scan, ...)` as the
stationary path does: it runs off the database and has no scan object. It
builds the divisor from the stored `auxillary` counters instead, via
`_rocking_normalization`, using the same monitor-name setting
(`reconstruction_monitor_corrections`) the other two paths use. Exposure time
is only there if the backend declares `exposure_time` in
`auxillary_counters` — ID31 does, `P212_tools` and the base `Scan` do not, and
a missing counter is skipped and recorded rather than failing the job.

## 3. Getting a green test run

**A fresh checkout reports ~80 failures.** They are all the unbuilt native
extension, not the branch. With the extensions built the suite is
**1120 passed, 3 skipped, 0 failed** (the 3 skips are a missing `arviz`).
Building also *unlocks* roughly 230 tests that are not collected at all
without it, so the unbuilt number is not simply "the green ones".

The environment on this machine is Python 3.14.7 (miniforge) with a working
orGUI already installed in `miniforge3\Lib\site-packages`. **Do not
`pip install -e .`** — it rebinds `import orgui` for that interpreter and the
repository owner wants the installed copy left alone. Build out of tree and
stage instead:

```powershell
# build dir must be SHORT: a path under the session scratchpad overruns
# MAX_PATH and meson dies in check_clock_skew with a FileNotFoundError
$vc = 'C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat'
$mf = 'C:\Users\timof\miniforge3'
$pre = "$mf;$mf\Scripts;$mf\Library\bin"
cmd /v:on /c "`"$vc`" >nul 2>&1 && set PATH=$pre;!PATH! && meson setup C:\Users\timof\obuild --buildtype=release"
cmd /v:on /c "`"$vc`" >nul 2>&1 && set PATH=$pre;!PATH! && meson compile -C C:\Users\timof\obuild"
cmd /v:on /c "`"$vc`" >nul 2>&1 && set PATH=$pre;!PATH! && meson install -C C:\Users\timof\obuild --destdir C:\Users\timof\ostage"

cd C:\Users\timof\ostage\Lib\site-packages
python -m pytest orgui/app/test orgui/datautils/xrayutils/test -q
```

Four traps, each of which cost an attempt:

1. `cl` is not on PATH; `vcvars64.bat` is required. It prints a harmless
   `'vswhere.exe' is not recognized` line and works anyway.
2. **`%PATH%` expands at parse time.** `vcvars && set PATH=<conda>;%PATH%`
   silently discards everything vcvars added, and meson then reports
   `Unknown compiler(s)`. Hence `cmd /v:on` and `!PATH!`.
3. meson and python are not on cmd's PATH even though they are on the bash
   shell's. The interpreter is in the miniforge root, the entry points in
   `Scripts`.
4. MAX_PATH, as above.

Do **not** copy the built `.pyd` into `orgui/` to make the worktree importable.
`install_subdir('orgui', ...)` ships whatever is in the source tree, so a stray
extension is then packaged over every future build — a trap this repository has
already been bitten by once.

`meson install --destdir` cannot write outside the staging tree, so the
installed copy is provably untouched; verify by checking that
`corrections/` is *absent* from `site-packages\orgui\datautils\xrayutils`,
not by comparing mtimes.

## 4. Rules this work established

Two of these came from the repository owner during the session and are now
also in the `AGENTS.md` files. They are the constraints a follow-up must keep.

* **`datautils` holds self-consistent physics modules. No UI or UI state may
  leak in.** No widget, no configuration object, and no *scan object* in
  `orgui/datautils/`. This is why `corrections/normalization.py` takes
  exposure and monitor *values* while `app/integration_corrections.py` keeps
  the `scan`-attribute lookup, and why the `use_lorentz` / `use_footprint`
  switches stayed on the app side.
* **One definition per factor.** Before this branch the mode→factor mapping
  lived in four places, the active area in three, the exposure/monitor divisor
  in three, and the per-pixel array in three verbatim copies. When adding a
  correction, put the physics in `corrections/` and call it; do not compute a
  factor inline in `orGUI.py` or `peak1Dintegr.py`.
* **Behaviour-preserving moves are verified against a failure-count baseline,
  not by inspection.** The refactor was checked at 80 failed / 782 passed
  before and after.

## 5. The wiring commit, as landed

`feat(phys)!` — it changed saved numbers in **both** modes. What it did:

1. **Normalization.** `peak1Dintegr.RockingPeakIntegrator._rocking_normalization`
   builds the per-frame exposure/monitor divisor from the stored `auxillary`
   counters (see section 2.2 for why not from a scan object) and it is applied
   **inside** the rocking integral, not to the finished integral: with a
   varying counting time or a drifting monitor the quantity Vlieg integrates is
   `Int N(omega)/(T M) d omega`, and dividing the result by a mean is only
   equivalent for constant counters. It rides the existing `C_corr` machinery,
   which already carried `C_illum_area` per `(s, omega)` point, so the error
   propagation followed for free.
2. **Radian.** `measurement.normalized_intensity(..., angle_unit="deg")`
   converts when `F2_hkl` is formed. The stored `croibg` and `int_interval`
   stay in the unit they were measured in — changing those would change the
   meaning of two saved columns for no gain.
3. **Acceptance.** `_rocking_acceptance` calls `out_of_plane_acceptance` with
   the region centre and vertical size stored per `s` point. Two traps here:
   the coordinate order (`surfaceAnglesPoint` takes pyFAI dimension 1 first,
   which is the *row*, and orGUI's `y` is the row while `x` is the column —
   every call site in the application passes them swapped), and `vsize` being
   the row extent, which follows from `detvsize, dethsize = detector.shape`
   at `orGUI.py:1096`. It is evaluated at the calibrated arm position, which
   costs nothing measurable because the span is arm-invariant to nine digits.
   With no calibrated detector reachable it warns and leaves `F2_hkl` on the
   acceptance-blind scale rather than failing the integration.
4. **Solid angle compensated, not removed.** The switch, its config key and
   its badge are untouched, and it still scales the intensity — a broad or
   diffuse feature wants a differential cross-section, and that capability was
   worth keeping. What changed is that `F2_hkl` divides it back out, via the
   new `detector.roi_mean_inverse_solid_angle`, measured over the same regions
   from the calibrated geometry alone (no image, so it runs outside the
   integration loop). Two things to know: the compensation is approximate at
   the `1e-6` level, because `pixel_factors` fuses the solid angle with the
   polarization so the applied mean is `<1/(Omega~ P)>` rather than a product
   of means; and for a rocking scan the correction was applied at *extraction*
   time, so whether to compensate is read from the configuration stored with
   the scan (`configuration/orgui/integration_corrections/json`) rather than
   from the current switch. An older database where that cannot be read warns
   and is left uncompensated rather than guessed at.

   A related trap, found by breaking it first: `useSolidAngleBox` is *also*
   the store the reconstruction persists its own solid-angle setting through
   (`config_data.py:383`/`:467` map the `solidAngle` integration option onto
   `CorrectionState.use_solid_angle`, and `ReconstructionDialog` mirrors its
   checkbox via `scanSelector.get/set_integration_options`). Removing the key
   would have silently disabled the correction in the one path that must keep
   it, and that is invisible from either file alone.
5. **Provenance.** A `reduction` group beside `F2_hkl` records the mode, the
   angle unit, which normalizations applied, whether the acceptance was
   applied, whether the solid-angle correction was compensated, the
   active-area assumption, and the acceptance array itself.

`test_rocking_and_stationary_paths_differ_by_the_missing_normalizations`
became `test_rocking_and_stationary_paths_agree`, and
`test_a_resized_region_of_interest_distorts_the_rocking_rod` became
`..._no_longer_distorts_...`. Both keep the *old* behaviour as a contrast
assertion via a `reduce=False` switch on the helper, so the factor that used to
be left behind cannot come back unnoticed.

## 6. After that

* **Absolute scale (#15)** needs two user inputs that have no config field yet:
  the incident flux density `Phi_0`, and the horizontal beam width for
  `activearea.beam_limited_area`. The footprint dialog already asks for the
  sample size and the beam profile. Everything else — `lambda`, `A_u`, `r_e` —
  is available.
* **Reflectivity comes for free** once `|F|^2` is absolute;
  `measurement.reflectivity_from_structure_factor` is the conversion. **F5 is
  now fixed**, which mattered most here: a reflectivity scan drives the arm,
  and the arm-blind polarization was a 10 % error at `2theta = 18` degrees and
  33 % at 30. The correction applies to the two direct-space integration
  paths; the reciprocal-space reconstruction still evaluates the polarization
  per pixel at the calibrated position.
* **`C_det` (F7)** is the only mechanism that can still break mode equivalence
  after the wiring, and it cannot be validated on simulated data. The
  real-data overlap comparison has now been run and **bounds** it: not
  detectable against 3 % scatter along a CTR whose regions cover the peak, a
  factor of 5 at the Bragg peak on the same rod. Modelling it is still open,
  and the bound is an upper limit for one rod on one sample, not a general
  result.

## 7. Predictions that turned out wrong

Recorded because each cost time and none was obvious in advance.

* **`active_area_footprint` should live next to the existing area code.** It
  should not exist at all. `w * min(L, h/sin(alpha))` is exactly
  `w * L * illuminated_area_fraction(alpha, L)` for a top-hat profile and only
  for that profile — verified to `3.3e-16`. A Gaussian of the same width
  departs by up to 19 %. The function was deleted, not moved.
* **"The refactor will remove ~100 lines from `orGUI.py`."** It is 27 added,
  10 removed — net *longer*, because the explanatory comment is worth more
  than the ten lines saved. Real length reduction there means extracting the
  `integrateROI` / `rocking_integrate` driver bodies, which is separate work.
* **"A moving detector arm changes `Delta_gamma`."** It does not. Driving the
  gamma arm shifts the exit angle at the ROI centre by the full arm angle but
  leaves the *span* identical to nine digits, because the rotation is about the
  axis gamma is measured around. A 40-degree delta-arm rotation moves it by one
  part in `1e5`. This is the opposite of the polarization factor, and it means
  an acceptance computed without arm bookkeeping is still usable.
* **"The run-to-run test-count drift is flakiness."** It was the interpreter
  changing under the session when Python 3.14 was installed mid-work.
* **`if C_arr is None: C_arr = np.ones(...)` in `orGUI.py` is dead code.** It
  is load-bearing: the branch that rebuilds `C_arr` sits inside `if HAS_ACCEL:`,
  so the NumPy-only path would be handed `None`. It is commented as such now.
* **pyFAI accepts mismatched coordinate array shapes.** It asserts
  `pos2.size == size`, so a scalar column with a per-frame region height — what
  the rocking path will pass — dies inside the extension.
  `acceptance._surface_gamma` broadcasts before the call.

## 8. Deliberately not done

* **The reconstruction's correction pass was not unified.** It is fused into
  `apply_correction_factors` with a bit-for-bit native/NumPy contract and
  variance propagation. It shares the *definition* (`pixel_factors`) but keeps
  its own streaming application; that was the right boundary.
* **F5 was not fixed while moving the code.** `pixel_factors` reproduces the
  historical arm-blind behaviour exactly. Changing it was a numerical fix that
  belonged in its own `phys` commit, not smuggled into a refactor -- and it
  landed as one. Note what it did *not* need: rebuilding the per-pixel array
  per frame. Because the integration paths reduce the polarization to a region
  mean anyway, the fix is a per-frame ratio of two region means, which is
  exactly 1 at the calibrated arm position and so leaves every fixed-arm scan
  bit-identical. `pixel_factors` itself is unchanged.
* **`meson.build` was not changed.** Its `exclude_directories: ['__pycache__']`
  only excludes the top-level directory, so 72 stale `cpython-312.pyc` files
  are sitting in the installed copy. Inert under 3.14, and the repository owner
  had previously declined a `meson.build` change for the related stale-`.pyd`
  problem, so it was left alone and raised separately.

## 9. Reproducing the key numbers

```powershell
pytest orgui/datautils/xrayutils/test/test_corrections_measurement.py
pytest orgui/datautils/xrayutils/test/test_corrections_acceptance.py
pytest orgui/datautils/xrayutils/test/test_corrections_activearea.py
pytest orgui/datautils/xrayutils/test/test_corrections_package.py
pytest orgui/app/test/test_scan_mode_equivalence.py
```

These five do not need the native extension. Everything else in the suite may,
so use section 3 before concluding anything from a failure.

The real-data check of `ctr_structure_factor_scale.md` section 5.1 needs data
that is not in the repository (LaNiO3 scan 61). To repeat it on another pair
of scans, run orGUI with `--nogui -i <script>` and, in that script: read the
configuration, load the scan, set the mask, load the reference reflections,
`fitExperiment`, then integrate the same rod twice — once stationary, once
rocking — before reducing the rocking curves. Four settings decide whether the
comparison means anything, and all four are now stored with the integration,
so on a database written after this branch they can be read back out of
`configuration/orgui/roi_integration` rather than guessed:

* the background margins, which dominate everything else (section 5.1);
* `auto_hsize`/`auto_vsize` — the reference used the *projected* horizontal
  size for the rocking pass, not the automatic one-pixel rule;
* the projected sample size, in meter;
* the reduction windows, which must be the same for both modes; a signal
  window several times wider than the peak shifts the ratio by tens of per
  cent on its own.

Compare on `F2_hkl`, not on the raw curves: match each rocking point to the
stationary points within a few thousandths of `s`, require signal-to-noise of
at least 5 on both sides, and drop regions whose `croi_pix` falls below the
nominal `vsize * hsize`, which is how a detector gap shows up.

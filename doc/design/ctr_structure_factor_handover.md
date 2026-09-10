# CTR structure-factor scale: implementation status and handover

> **Status as of 2026-09-09.** Branch `claude/ctr-structure-factor-9633bc`,
> four commits ahead of `master`, nothing pushed.
>
> The physics analysis is complete and quantified, the reduction and the
> acceptance estimator exist and are tested, and the correction factors have
> one home. **No saved number has changed yet**: the GUI integration paths do
> not call any of the new reduction. Issues
> [#82](https://github.com/tifuchs/orGUI/issues/82) and
> [#15](https://github.com/tifuchs/orGUI/issues/15) are therefore *analysed and
> equipped* but not closed.
>
> This document is the handover: what exists, how to run it, what to do next,
> and which of my predictions turned out wrong. The physics itself is in
> [`ctr_structure_factor_scale.md`](ctr_structure_factor_scale.md) — read that
> first for *why*; this one is *where things stand*.

## 1. The one-paragraph summary

A rocking scan and a stationary scan of the same rod do not currently produce
the same `F2_hkl` in orGUI. The ratio is exactly
`T_omega * monitor_omega * Delta_gamma_in_degrees`, measured to seven digits on
simulated data. Three things are missing from the rocking path: exposure/monitor
normalization, integration in radian rather than degrees, and division by the
out-of-plane detector acceptance. All three now have working implementations in
`orgui/datautils/xrayutils/corrections/`; none of them is wired in. A fourth,
F6, affects **both** modes: a ROI sum is already a complete angular integral,
so the solid-angle correction has to stop being applied to it, or the two
modes stay apart by `<1/Omega~>` even after the rocking path is fixed.

## 2. What is on the branch

| commit | what it did |
|---|---|
| `777d887` | `refactor: collect correction factors into one corrections package` |
| `4193c80` | `feat: reduce integrated intensities to \|F_hkl\|^2 on one scale` |
| `e25b8df` | `docs: record the rocking/stationary structure-factor scale analysis` |
| `6959191` | `feat: estimate the out-of-plane detector acceptance` |

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
  detector.py        89   per-pixel solid angle and polarization
  normalization.py  103   counting time and monitor, from values
  roi.py             93   CorrectionFactors, roi_mean_correction
  measurement.py    605   mode dispatch, master equation, reflectivity
```

`orgui/datautils/xrayutils/geometrycorrections.py` and `beamprofile.py` remain
as re-export aliases because both were released under those names.
`test_corrections_package.py` asserts the aliases hand out the *same objects*,
not merely that they import.

### 2.2 What is wired, and what is not

| caller | uses the package for | still does its own thing |
|---|---|---|
| `orGUI.integrateROI` (stationary) | `pixel_factors`, `mode_components`, `normalization_divisor`, `C_illum_area` | — |
| `orGUI.rocking_integrate` | `pixel_factors` | — |
| `peak1Dintegr.integrate` (rocking) | `mode_components` | **no normalization, degrees, no `Delta_gamma`** |
| `reconstruction_job` | `pixel_factors` | own native-fused application, own normalization loop |

Nothing calls `measurement.structure_factor_squared`,
`measurement.normalized_intensity`, `measurement.angular_factor`,
`activearea.*` or `acceptance.*` outside the tests. That is the gap to close.

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

## 5. The next commit, in detail

Wiring the rocking path, plus the F6 solid-angle removal, which also touches
the stationary path. This **changes saved numbers in both modes**, so it is
`feat(phys)!` with a `BREAKING CHANGE:` footer per the repository's commit
convention.

Steps 1-3 are `peak1Dintegr.RockingPeakIntegrator.integrate`; step 4 is that
plus `orGUI.integrateROI`:

1. Build the exposure/monitor divisor the way the stationary path does, via
   `integration_corrections.normalization_divisor(scan, ...)`. Use the
   **per-frame** counting time, not the sum over the scan: the rocking angle is
   the integration variable, the time is not.
2. Convert the trapezoid integral from degrees to radian —
   `measurement.normalized_intensity(..., angle_unit="deg")` does both this and
   step 1.
3. Divide by `acceptance.out_of_plane_acceptance(detector, row, column,
   row_size, alpha, gamma_arm, delta_arm)`. The ROI centre and size per `s`
   point are already in the database group `integration/`; the arm angles come
   from `orgui.backend.scans.scan_arm_angles`.
4. Stop applying the solid-angle correction to ROI sums — **in both modes**,
   so this touches `orGUI.integrateROI` as well. F6 is resolved this way (see
   the physics doc): a ROI sum is already a complete angular integral, and
   leaving the correction in means the modes still differ by `<1/Omega~>` —
   0.7 % at 1 m, 7 % at 0.3 m — because a rocking scan's solid-angle factor
   cancels against its acceptance and a stationary scan's has nothing to
   cancel against. In practice: drop the `solid_angle=` argument at
   `orGUI.py:1595` and `:5765`, then the `useSolidAngleBox` widget, the
   `solidAngle` key in `get_integration_options` and the `SOLA` badge.
   `set_integration_options` ignores keys it has no branch for, so old
   configuration files still load. The reconstruction path has its own switch
   (`reconstruction_job.py:1254`) and must keep it.
5. Store the acceptance, the applied normalization and the active-area
   assumption next to `F2_hkl`. Without them a saved rod cannot be put on a
   common scale after the fact.

Then `orgui/app/test/test_scan_mode_equivalence.py::test_rocking_and_stationary_paths_differ_by_the_missing_normalizations`
**must flip** from a characterization test to a plain equality. Its docstring
says so. It currently asserts the gap *is* `T * monitor * Delta_gamma_deg`; if
it starts failing after a wiring change, that is the change working.

**F6 is settled** (physics doc, F6): ROI-summed integration stops applying the
solid-angle correction, in both modes. Two things fell out of measuring it
that matter for step 3: `Delta_gamma` does **not** become unnecessary — it is a
size factor, the solid angle an obliquity factor, and they overlap only in the
obliquity — and a constant nominal `n * pixel / dist` is **not** a valid
"pixel-derived" acceptance, being wrong by `1/cos(theta)` (1.2 % at 9 degrees).
`out_of_plane_acceptance` is still what the rocking path divides by.

## 6. After that

* **Absolute scale (#15)** needs two user inputs that have no config field yet:
  the incident flux density `Phi_0`, and the horizontal beam width for
  `activearea.beam_limited_area`. The footprint dialog already asks for the
  sample size and the beam profile. Everything else — `lambda`, `A_u`, `r_e` —
  is available.
* **Reflectivity comes for free** once `|F|^2` is absolute;
  `measurement.reflectivity_from_structure_factor` is the conversion. But fix
  **F5** first: a reflectivity scan drives the detector arm, and the
  polarization is still evaluated at the calibrated position, which is a 10 %
  error at `2theta = 18` degrees and 33 % at 30.
* **`C_det` (F7)** is the only mechanism that can still break mode equivalence
  after the wiring, and it cannot be validated on simulated data. It needs the
  real-data overlap comparison.

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
  historical arm-blind behaviour exactly. Changing it is a numerical fix that
  belongs in its own `phys` commit, not smuggled into a refactor.
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

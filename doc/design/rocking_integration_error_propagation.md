# Rocking-scan and stationary-scan integration: error-propagation and background bugs, and a plan to fix them

> **Status: Steps 0-6 complete; quantitative fixes for findings #1-#6 applied, #7-#8 documented as known approximations per Step 6.** This document
> is the plan referenced from [issue #34](https://github.com/tifuchs/orGUI/issues/34)
> ("Untested error propagation in rocking integration"). It records every bug
> found in `orgui/app/peak1Dintegr.py` and the image-level ROI integration in
> `orgui/app/orGUI.py`, verified numerically (not just by inspection), and lays
> out the order to fix them in and the tests that should pin each fix down.
> `orgui/app/peak1Dintegr.py` is not touched by any work in flight on other
> branches as of this writing (2026-08-12) — the analysis below is against
> `master` at that date.
>
> The pure-numpy aggregation core of `RockingPeakIntegrator.integrate()`
> (`peak1Dintegr.py:1017`-`1245` as originally numbered) has been extracted,
> behavior-preserving, into a module-level function
> `_compute_rocking_integration` so it is callable without Qt or an on-disk
> database. `orgui/app/test/test_peak1Dintegr.py` (new) exercises it directly
> and pins findings #1-#4 with tests that **fail against current `master`**,
> confirmed by running them: `pytest orgui/app/test/test_peak1Dintegr.py -q`
> reports 4 failed. The full `orgui/app/test/` suite (261 tests) still passes,
> confirming the extraction changed nothing about `integrate()`'s behavior.
> Finding #5 is not separately pinned — it is entangled with #1 in the
> background-aggregation test (see that test's docstring) since both bugs
> touch the same trapezoidal-slicing code path; a dedicated #5-only test can
> be added once #1 is fixed and the interaction can be isolated.

Neither `orgui/app/peak1Dintegr.py` nor the ROI-sum integration path in
`orgui/app/orGUI.py` has any test coverage today. `orgui/app/test/` covers the
accelerators (`test_roi_sum_accel.py`), peak search, config, and masks, but
nothing exercises `RockingPeakIntegrator.integrate()` or the `sumImage`
closures in `orGUI.py`. That absence is exactly why issue #34 exists, and it
is why every fix below should land with a regression test rather than trusting
review.

## Severity ranking

Ordered by how badly wrong the *saved numbers* are, not by how easy the fix
is:

1. **§1 — background under-subtracted with >1 bg ROI.** Wrong intensities
   (`croibg`, `F2_hkl`), not just wrong errors. Silent, and hits the default
   ROI configuration (`Add All` creates `bg_1` and `bg_2`).
2. **§2 — bg error accumulator discarded and replaced with a value that
   doesn't scale with counting time.** Wrong errors, silently, and worse for
   the weak/background-dominated reflections that need error bars most.
3. **§6 — non-accelerated rocking `sumImage` return-arity mismatch.** Crashes
   for any ROI count other than 2, and silently reads garbage for exactly 2.
   Numpy-only installs (`HAS_ACCEL == False`) hit this on every rocking
   integration.
4. **§4 — `I_corr_error` computed as a ratio through a value that can be
   zero.** Produces `inf`/`nan` in saved data rather than a wrong-but-finite
   number.
5. **§3 — `F2_hkl_errors` uses the uncorrected `raw_croibg_errors`.** Wrong
   only when the footprint correction is enabled; silently correct when it
   isn't, which is presumably why it was never noticed.
6. **§5 — trapezoid error weights wrong at ROI boundaries, and integration
   window off by one grid point relative to `int_interval`.** Systematic
   over-estimate of the error (~4-5% for a 20-point ROI, worse for smaller
   ROIs) plus a background-scaling bias when signal and background ROIs have
   different point counts.
7. **§7 — bg-image "method 2" error reuses the method-1 formula.** An
   approximation that is only exact for a spatially flat background image;
   not documented as such.
8. **§8 — fitted polynomial background variance.** Exact at `fit_order = 0`
   only; underestimated at order 1 and 2, since the legacy
   `(cpixel/bgpixel)² · bgroi` formula assumes no extrapolation bias and
   Poisson-weighted samples, and the polynomial fit does neither.

§1-§6 are argued to be outright bugs (the code does not implement what its own
comments/nearby code say it should). §7-§8 are known approximations that
should be documented as such, not necessarily rewritten — see "Out of scope
for the initial fix" below.

## Finding detail

> **Line numbers below are pre-extraction (Step 0 moved this logic into
> `_compute_rocking_integration`, see the plan).** Every `peak1Dintegr.py:NNNN`
> reference in this section describes where the code lived in
> `RockingPeakIntegrator.integrate()` at the time of the initial analysis
> (2026-08-12, `master`). After Step 0, the same lines live inside the
> module-level `_compute_rocking_integration` function
> (`peak1Dintegr.py:84` onward) with the same relative order and the same
> bugs, unchanged — only their line numbers and indentation shifted. Search
> by the quoted code rather than trusting the absolute numbers; they will
> drift further as each fix step lands.

### §1 — background double-thinned when more than one bg ROI is anchored

`orgui/app/peak1Dintegr.py:1200` (context: `integrate()`, the `bg`-summation
loop starting at `peak1Dintegr.py:1198`):

```python
ratio = (sig_interval / bg_interval) * (
    int_data[roikey]["int_interval"] / bg_interval
)
bgroi += int_data[roikey]["cnts"] * ratio
```

Intended quantity: mean background density across *all* bg ROIs at this
`s`-point, rescaled to the width of the signal window:

```
B_hat = (sum_i I_i / sum_i b_i) * s = sum_i I_i * (s / b)
```

The code multiplies by an extra `int_data[roikey]["int_interval"] /
bg_interval` — the fractional width of *this one* bg ROI relative to the
total bg width. With `N` bg ROIs of equal width that extra factor is `1/N`
for each of them, so the sum comes out `1/N` of the correct background.
`onAddAllROI` (`peak1Dintegr.py:1459` in `orGUI.py` — actually
`RockingPeakIntegrator.onAddAllROI`) creates `bg_1` and `bg_2` by default, so
**the default configuration halves the subtracted background**, inflating
every `croibg`/`F2_hkl` value produced through the "Add All" + "integrate"
path.

Verified numerically (one bg ROI of width 0.1 and count 50, one more of width
0.1 and count 50, both mapped to a 0.2-wide signal window):

```
bg as coded : 50.0   correct: 100.0   ratio: 0.5
```

**Fix:** drop the intra-loop `int_data[roikey]["int_interval"] / bg_interval`
factor — `ratio = sig_interval / bg_interval` is already summing over all
`roikey` starting with `"bg"`, so each ROI's raw count already carries its own
weight through `int_data[roikey]["cnts"]`. Apply the same fix to
`bgroi_errors`/`raw_bgroi_errors` accumulation on the neighboring lines
(§2 below covers the separate `np.sqrt` bug there, but the `ratio` fix must
land first or the error fix will encode the same 1/N bias).

### §2 — bg error accumulator discarded

`orgui/app/peak1Dintegr.py:1219-1220`:

```python
raw_bgroi_errors = np.sqrt(raw_croi_errors)
bgroi_errors = np.sqrt(croi_errors)
```

The loop immediately above (`peak1Dintegr.py:1198-1217`) accumulates
`bgroi_errors += (int_data[roikey]["cnts_errors"] * ratio) ** 2` (and the
`raw_` equivalent) across all bg ROIs — that is the right pattern, sum of
variances then one final `sqrt`. But the two lines above throw that away and
replace it with the square root of the unrelated, already-square-rooted
*signal* error. The result does not scale with exposure time the way a
Poisson-derived quantity must:

```
T (relative exposure)   true sigma_bg   coded sigma_bg
1                        3.50            3.16
10                      11.07            5.62
100                     35.00           10.00
1000                   110.68           17.78
```

and it feeds directly into `croibg_errors = sqrt(croi_errors**2 +
bgroi_errors**2)` (`peak1Dintegr.py:1224`), so the underestimate for weak,
background-dominated reflections runs from ~1% at short exposure to over 5%
here and grows without bound as exposure increases, since `bgroi_errors`
stops scaling with `sqrt(T)` at all.

**Fix:** `raw_bgroi_errors = np.sqrt(raw_bgroi_errors)` and `bgroi_errors =
np.sqrt(bgroi_errors)` — i.e. finish the sum-of-squares that was already
built, instead of recomputing from the wrong variable. Must be fixed together
with §1 (the `ratio` bug feeds the same accumulator).

### §3 — `F2_hkl_errors` uses the wrong (uncorrected) numerator

`orgui/app/peak1Dintegr.py:1230-1231`:

```python
F2_hkl = croibg / (C_Lorentz * C_rod_intersect)
F2_hkl_errors = raw_croibg_errors / (C_Lorentz * C_rod_intersect)
```

In the affected implementation, `croibg` has both then-applied footprint
factors (`C_flux_on_sample`, `C_illum_area`) folded in when the footprint
button is checked; `raw_croibg_errors` does not
— it is the error of the *uncorrected* integral. The two only agree when the
footprint correction is off (`C_corr == 1` everywhere), which is presumably
why this was never flagged.

**Fix:** `F2_hkl_errors = croibg_errors / (C_Lorentz * C_rod_intersect)`.

### §4 — `I_corr_error` as a ratio through a possible zero

`orgui/app/peak1Dintegr.py:1092-1095`:

```python
I_raw_error = np.sqrt(
    np.sum((cnts_errors * deltaaxis[idx_from:idx_to]) ** 2)
)  # to be checked!
I_corr_error = (I_raw_error / I_raw) * I_corr
```

`I_raw` can legitimately be zero or near-zero — an empty ROI
(`idx_from == idx_to`, see §5) or a background window whose counts fluctuate
about zero after any subtraction upstream. `I_corr_error` then becomes `inf`
or `nan` and propagates into the saved HDF5 dataset. It is also an unnecessary
detour: `C_corr` is already known per-point, so the corrected error can be
propagated directly without going through the raw-to-corrected ratio.

**Fix:** propagate directly, mirroring the `I_raw_error` computation but
including `C_corr`:

```python
I_corr_error = np.sqrt(
    np.sum(((cnts_errors / C_corr) * deltaaxis[idx_from:idx_to]) ** 2)
) * sign_interval  # sign doesn't matter for an error but keep consistent typing
```

(`sign_interval` multiplication on an error term is a no-op since errors are
non-negative by construction here; keep it out or make explicit that errors
are unsigned — decide during implementation and add a comment, don't leave it
implicit.)

### §5 — trapezoid error weights wrong at ROI boundaries; `int_interval`/slice mismatch

`orgui/app/peak1Dintegr.py:1036-1040` builds `deltaaxis` once for the *whole*
scan axis:

```python
dx = np.diff(axis)
deltaaxis = np.empty_like(axis)
deltaaxis[0] = dx[0] / 2
deltaaxis[-1] = dx[-1] / 2
deltaaxis[1:-1] = (dx[:-1] + dx[1:]) / 2
```

then slices it per-ROI (`deltaaxis[idx_from:idx_to]` at
`peak1Dintegr.py:1093`). For a sub-interval, the correct trapezoidal weight is
`dx/2` at *both ends of that sub-interval*, but the sliced array only gets the
true `dx/2` weight if the ROI happens to start/end at the axis boundary — for
any interior ROI (the normal case) both endpoints get the "interior" weight
`(dx_{k-1}+dx_{k+1})/2`, roughly twice too large. Verified for a 20-point ROI
on a uniform 0.01-wide grid:

```
sigma_I coded: 0.04472   true: 0.04301   overestimate factor: 1.0398 (~4%)
```

Smaller ROIs (fewer points) are proportionally worse, since the two
mis-weighted endpoints are a larger fraction of the total.

Separately, in the same block: `cnts = croibg[idx_from:idx_to]`
(`peak1Dintegr.py:1065`) is a half-open slice that drops the point at
`idx_to`, while `int_interval = np.abs(axis[idx_to] - axis[idx_from])`
(`peak1Dintegr.py:1059`) spans the *closed* interval including it. The
integral is therefore computed over one fewer grid step than the width it is
later normalized against (`sig_interval`, `bg_interval` in §1, and
`C_Lor`/`C_rod`/footprint means at `peak1Dintegr.py:1070-1079`). This biases
the §1 background scaling whenever signal and background ROIs span different
numbers of points, and it biases every `np.mean(...[idx_from:idx_to])`
correction-factor average by the same one-point truncation.

**Fix, two parts:**
- Compute a local `deltaaxis` from the sliced sub-axis
  (`axis[idx_from:idx_to+1]`) rather than slicing the globally-built array —
  or equivalently, correct the two endpoint weights after slicing. Prefer the
  local recomputation: it's less error-prone than patching endpoints in place
  and keeps this function's error propagation self-contained.
- Change every `idx_from:idx_to` ROI slice in `integrate()` (counts, weights,
  correction factors, auxiliary counters — `peak1Dintegr.py:1065-1113`) to
  `idx_from:idx_to + 1`, and add a regression assertion that
  `axis[idx_from:idx_to+1].size` matches the number of trapezoid points used
  for `int_interval`. This is a wider-reaching change than §1-§4 — audit
  every consumer of `idx_from`/`idx_to` in the function before changing the
  slice bound, since `_trapz_impl(cnts, axis[idx_from:idx_to])` calls
  (`peak1Dintegr.py:1086`, `1089`, `1105`) must move together with `cnts`, and
  the `C_corr` array size (`peak1Dintegr.py:1068`, sized from `cnts.size`)
  must track it too.

### §6 — non-accelerated rocking `sumImage` returns the wrong arity

`orgui/app/orGUI.py:1783-1819` (`sumImage` closure used when `HAS_ACCEL` is
`False`, in the rocking-scan integration path) returns a single array
`all_counters1` of shape `(n_rois, 4)`. The consumer expects two:

```python
img_counters, Carr_counters = f.result()  # orGUI.py:1840
```

For any `n_rois != 2` this raises `ValueError: too many values to unpack` (or
too few). For exactly `n_rois == 2` it silently succeeds by unpacking two rows
of the counters array as if they were `(img_counters, Carr_counters)` —
producing counters that are neither. This path is only reached when the
compiled ROI accelerator (`_roi_sum_accel`) is unavailable
(`HAS_ACCEL == False`, e.g. a pure-Python / non-compiled install), so it is
plausible this was never exercised in the maintainer's own dev environment,
which matches "less tested" in the issue.

Contrast with the stationary-scan fallback at `orGUI.py:5765` and neighbors,
which correctly returns `all_counters, Carr_counters, BgImg_counters` (three
arrays) and is consumed with matching unpacking (`orGUI.py:5900`-ish region,
already checked in the earlier pass over this file). The rocking closure looks
like an incomplete port of that pattern — it neither computes correction-array
sums (`Carr_counters`) nor separates background-image counters, and also does
background subtraction and correction-array multiplication in-place on the
raw image (`orGUI.py:1792`, `1799`) before summing, which destroys the Poisson
assumption `sqrt(croi)` needs downstream and can make `cnts` negative.

**Fix:** rewrite the `HAS_ACCEL == False` rocking `sumImage` to match the
structure of the stationary-scan fallback: sum raw counts and correction-array
sums *separately* (never multiply the image by `C_arr` before computing
`croi`), return `(all_counters, Carr_counters)` — and `BgImg_counters` too if
this path is meant to support the background-image workflow at all; if not,
that should be an explicit `NotImplementedError`/warning rather than a silent
gap. Confirm against `fill_counters` (`orGUI.py`, `def fill_counters` near
where `sumImage` is defined) whether it already does the right per-ROI, no-
premultiply summation and can simply be reused instead of hand-rewriting.

### §7 — bg-image "method 2" error is a flat-background approximation

`orgui/app/orGUI.py` (rocking-image path around what was `orGUI.py:1860` and
its stationary-scan mirror around `orGUI.py:5288` prior to line drift; find via
`croibg1_bgimg_err_a = np.sqrt(` — currently `orGUI.py:2026` and
`orGUI.py:5918`/`5961`): "method 2" scales the background image by
`factor = bgroi / bgimg_bgroi_norm` before subtracting it, but its error term
reuses the method-1 expression `sqrt(croi + (cpixel/bgpixel)^2 * bgroi)`,
ignoring `factor`. Propagating through `factor` properly gives something
closer to `sqrt(croi + (bgimg_croi_norm/bgimg_bgroi_norm)^2 * bgroi)` — the
coded version is only exact when the background image is spatially flat
across both the center and background ROI footprints (so that `croi`'s
"as-if-Poisson" term absorbs the discrepancy). Not necessarily wrong for a
real flat-field background, but currently undocumented as an approximation.

**Fix (documentation-only, not correctness-critical):** add a docstring/inline
comment at both call sites stating the flat-background-image assumption
explicitly, so a future reader — or a user with a genuinely structured
background image — knows the caveat. A follow-on ticket can implement the
`factor`-aware propagation if someone needs it; don't block the §1-§6 fixes on
this.

### §8 — fitted polynomial background variance underestimated above order 0

`orgui/app/cpp/roi_sum_cpp.cpp:1789` (`processImage_polybg_Carr`,
`fit_local_polynomial_background`/`stored_background` machinery starting
around `roi_sum_cpp.cpp:1256`):

```cpp
stored_background = fitted_background * fit.sample_count / center.pixels;
```

This is constructed specifically so the legacy
`(cpixel/bgpixel)^2 * bgroi` error formula downstream in `orGUI.py`
reproduces the correct variance — but that identity only holds when the
"background" really is a simple unweighted count over `bgpixel` pixels, i.e.
`fit_order == 0` (flat/mean background, no per-pixel extrapolation). At
`fit_order` 1 or 2 the polynomial extrapolates from the four side ROIs into
the center ROI footprint; the true variance of that extrapolated estimate
depends on the fit's covariance (via the local design matrix) and is larger
than the sampled-mean variance the stored formula assumes, and the fit itself
is unweighted least squares rather than Poisson-weighted. Net effect: errors
are **underestimated** whenever `FittedBackgroundOrder >= 1` is selected in
the ROI options.

**Fix (larger, defer):** this needs either (a) propagating the fit's actual
covariance through to a per-ROI background variance in the C++ kernel, or (b)
a documented caveat plus a GUI-level warning next to the "Fitted local
background" order selector, similar to the existing `HAS_ACCEL`/repair
warnings already in `orGUI.py` (e.g. `orGUI.py:1473`-style pattern: "Fitted
local background requires..."). Recommend (b) for the initial pass — it's a
one-line warning next to existing warnings of the same shape — and open a
follow-up issue for (a), since a correct covariance propagation touches the
C++ kernel and needs its own accuracy validation against synthetic Poisson
data.

## Plan

Work in the order below; each step is small enough to review independently
and each should land with the test that pins it, in the same commit or PR
step, per the repository's "add a focused regression test or precise
docstring" rule.

### Step 0 — regression test scaffold (done)

Rather than fighting `RockingPeakIntegrator.integrate()`'s Qt/HDF5
dependencies (`ConfigData.from_gui`, `self.database.nxfile`,
`self.lorentzButton.isChecked()`, ...) in a test fixture, the pure-numpy
aggregation core was extracted into a module-level function,
`_compute_rocking_integration(s_array, axis, croibg_curves,
croibg_errors_curves, roi_info, aux, use_lorentz, use_footprint, C_Lor=1.0,
C_rod=1.0, C_flux_on_sample=1.0, C_illum_area=1.0, progress_callback=None,
should_cancel=None)`, defined just above `class RockingPeakIntegrator` in
`peak1Dintegr.py`. `integrate()` now calls it and unpacks the returned dict
(`int_data`, `croi(_errors)`, `bgroi(_errors)`, `croibg(_errors)`,
`raw_*` variants, `auxil`, and `F2_hkl(_errors)` when Lorentz is on) instead
of doing the aggregation inline. This is a pure move — every line of the
aggregation logic, **bugs included**, was carried over unchanged; nothing
about `integrate()`'s behavior changed. Confirmed by running the full
`orgui/app/test/` suite before and after: 261 tests, all passing both times.

`orgui/app/test/test_peak1Dintegr.py` (new) calls
`_compute_rocking_integration` directly — no `QApplication`, no HDF5 file —
with hand-constructed, piecewise-constant synthetic curves whose correct
trapezoidal integral is exactly `value * width` regardless of grid
alignment, and independently re-derives the correct aggregation in a small
reference implementation (`_reference_integration`, also in the test file)
matching this document's description of what the fixed code should do. Four
tests currently pin four of the five findings above:

- `test_background_aggregation_matches_reference_with_two_unequal_bg_rois` —
  §1 (and, entangled, §5: see the test's docstring) — asserts `bgroi` and
  `croibg` against the reference. **Currently fails**: coded `bgroi` is
  `24.18` against a reference of `43.6` for the two-unequal-bg-ROI case
  worked out in "Finding detail" §1 above.
- `test_background_error_scales_with_exposure_time` — §2 — asserts
  `bgroi_errors` scales as `sqrt(T)` between two exposure levels 100x apart.
  **Currently fails**: actual ratio is `~3.16` (`100**0.25`) against the
  expected `10.0` (`100**0.5`), i.e. `bgroi_errors` scales as `T**0.25`, not
  `T**0.5`, exactly the mismatch predicted in "Finding detail" §2.
- `test_F2_hkl_errors_uses_corrected_croibg_errors_under_footprint_correction`
  — §3 — asserts `F2_hkl_errors` against `croibg_errors / denom` (`denom`
  recovered from the already-computed `F2_hkl`/`croibg`). **Currently
  fails**: off by the footprint correction factor, as expected.
- `test_integrated_error_is_finite_when_raw_signal_integral_is_zero` — §4 —
  asserts `np.all(np.isfinite(croibg_errors))` for an all-zero-count signal
  ROI. **Currently fails**: `croibg_errors` is `nan`, with the exact
  `0/0`-division `RuntimeWarning` predicted in "Finding detail" §4 surfacing
  in the test output.

§5 (endpoint trapezoid weights, `idx_from:idx_to` vs. the closed interval
`int_interval` is measured against) is not pinned by its own isolated test:
every synthetic curve here is piecewise-*constant*, which makes the
half-open-slice truncation and the ratio bug (§1) both visible in the same
`bgroi`/`croibg` comparison, and they cannot be cleanly separated without
first fixing one of them. Once §1 is fixed (Step 1 below),
`test_background_aggregation_matches_reference_with_two_unequal_bg_rois`
will *still* fail (from §5 alone) until §5 is also fixed (Step 4) — which is
correct and expected, not a test bug. A dedicated #5-only regression (e.g.
comparing `int_interval` against the actual number of points integrated)
can be added once #1 lands and the residual gap is entirely attributable to
#5.

**Verified**: `pytest orgui/app/test/test_peak1Dintegr.py -q` →
`4 failed, 2 warnings`. `pytest orgui/app/test/ -q` → `4 failed, 261 passed`
(the 4 are the new, intentionally-failing tests above; nothing pre-existing
regressed). `ruff check orgui/app/peak1Dintegr.py
orgui/app/test/test_peak1Dintegr.py` → clean.

### Step 1 — fix §1 and §2 together (background ratio + error accumulator)

These two share the same accumulator loop
(`peak1Dintegr.py:1192`-`1220`) and should land as one change:
correct the `ratio` expression, and stop discarding the accumulated
`bgroi_errors`/`raw_bgroi_errors`. Commit as `fix(app,phys): correct rocking
background subtraction and error accumulation` per the `phys` scope rule in
`AGENTS.md` (this changes measurement results). Re-run Step 0's tests.

### Step 2 — fix §3 (`F2_hkl_errors`)

One-line change once §1/§2 land (the corrected `croibg_errors` needs to exist
first for this to be meaningful). Same `phys` scope.

### Step 3 — fix §4 (`I_corr_error`)

Direct propagation instead of the raw/corrected ratio. Verify against the
zero-crossing ROI test case from Step 0. `phys` scope.

### Step 4 — fix §5 (trapezoid weights and slice/interval consistency)

The riskiest single change — touches every consumer of `idx_from`/`idx_to` in
`integrate()`. Do this as its own commit, separate from §1-§4, so a revert is
cheap if the slice-width change has a consumer this plan missed. Grep
`idx_from` and `idx_to` in `peak1Dintegr.py` immediately before writing the
patch to get the current, authoritative list of consumers — don't rely on the
line numbers in this document, which will drift as §1-§4 land first. `phys`
scope.

### Step 5 — fix §6 (non-accelerated rocking `sumImage`)

Independent of §1-§5 (different file, different code path). Needs its own
test: a rocking-integration test run with `HAS_ACCEL` patched/mocked to
`False` (or run in an environment without the compiled extension, if CI has
one) to actually exercise this closure — Step 0's default test run will only
hit this path if the test environment lacks the compiled accelerator, so
don't assume Step 0 alone covers it. `fix(app)` scope (not `phys` — this is a
crash/data-integrity fix in a path that, per the arity mismatch, cannot
currently produce trustworthy numbers at all for `n_rois != 2`).

### Step 6 — document §7 and §8 as known approximations

Docstring/comment additions at the four sites named in §7, plus a GUI warning
next to the fitted-background-order control per §8's option (b). No `phys`
commit needed for §7 if it's comment-only; a one-line GUI warning for §8 is
`fix(app)` or `docs(app)`, whichever the actual diff ends up looking like —
decide at implementation time. Open a separate tracking issue for §8's proper
covariance-propagation fix (option (a)) rather than doing it inline here.

### Step 7 — changelog and topic docs

Per `doc/AGENTS.md`, this changes "calculation results... or anything that
changes previously produced numbers" — add an `## [Unreleased]` bullet to
`CHANGELOG.md` once §1-§5 land, under whatever theme heading matches existing
style (likely "Scientific and analysis additions:" or a bug-fix-focused
heading if one exists by then — check current `CHANGELOG.md` structure at
implementation time). Consider whether `doc/source/image_integration.rst`
(`image_integration.rst:43`-`62`, which already describes the rocking
hklscan/Bragg tabs) needs a short note about background-ROI count handling —
likely yes for §1, since users configuring more than 2 bg ROIs manually should
know the old and new behavior differ.

## Out of scope for the initial fix

- §7 and §8's *quantitative* fixes (as opposed to documenting them as known
  approximations) are separate, larger pieces of work — flagged above but not
  planned in detail here.
- Any change to the `IntegrationCorrectionsDialog` footprint-correction
  physics itself (`peak1Dintegr.py:991`-`1010`) — not touched by issue #34 and
  not found to be wrong during this pass; only its *interaction* with
  `F2_hkl_errors` (§3) was.
- Re-deriving the Lorentz-factor formulas (`C_Lor`, `C_rod` at
  `peak1Dintegr.py:976`-`989`) — out of scope for an error-propagation issue;
  flag separately if their correctness is ever questioned.
- Any change to the compiled accelerator's Poisson formula for the
  already-correct paths (`processImage_bg_Carr`, `processImage_Carr`, etc.) —
  verified correct during this pass (`(cpixel/bgpixel)^2` scaling, correction
  array applied after the fact) and should not be touched.

## 2026-08 footprint-convention follow-up

The statement above that footprint physics was not found wrong is superseded
by a later comparison with Vlieg (1997), equations (40)--(41) and (54).
`C_illum_area = C_flux_on_sample / (p_max L sin(alpha))` is already the
numerically illuminated surface integral, up to an angle-independent scale.
The integration multiplied it by `C_flux_on_sample` again, so an offset beam
entered quadratically. Rocking and stationary integration now divide only by
`C_illum_area`; the flux fraction remains saved as diagnostic provenance.
The same audit found that stationary equation (54) has no rod-interception
factor, so stationary `F2_hkl` now divides by its Lorentz factor alone.

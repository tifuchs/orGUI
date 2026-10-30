# Water DWBA profiling and optimization

> Superseded architecture (2026-10-29): the fitted sharp/tanh paths now use
> onset-only optical interfaces and analytic water perturbations, including
> the infinite Hosemann sum. The numerical-profile fitting and quadrature-cache
> optimizations recorded below are historical. See the current water section
> in `doc/source/dwba.rst` and `_CTRwater_analytic.py`. Explicit erf validation
> retains its continuous numerical solver.


## Current analytic implementation

Reproduction: `python -m examples.CTR.benchmark_water_DWBA --output tmp/water_analytic_compact --profile`.
The measurements below use the same 180-point 00L example and one thread
for BLAS/OpenMP in the benchmark process. Total time includes preparation;
the water column is the median of seven **uncached** perturbation calls.
Machine load varied, so small differences in total times are not speedup claims.

| Surface | Model | Optical interfaces | First total (ms) | Water kernel (ms) |
|---|---|---:|---:|---:|
| rough | dry | 16 | 131.7 | — |
| rough | step_tanh | 30 | 248.5 | 5.29 |
| rough | step_1ML_tanh | 30 | 277.2 | 9.54 |
| rough | layered | 30 | 169.9 | 6.20 |
| flat | dry | 10 | 110.8 | — |
| flat | step_tanh | 11 | 235.7 | 1.07 |
| flat | step_1ML_tanh | 11 | 109.8 | 1.26 |
| flat | layered | 11 | 98.9 | 1.20 |

Flat water adds exactly one interface to the dry solid reference. Rough
water adds the 14 physical terrace onsets. Adsorption does not increase
that count, and reducing paracrystal damping does not enumerate more layers.
The fast water path has no quadrature or continuous-field sampling.
The dominant work is now native solid scattering and structural preparation.

The refreshed RuO2/TiO2 water notebook compares the same physical step
against explicit numerical erf confirmation: maximum 00L intensity
difference is 0.0494% above its stated intensity floor; median 0.00103%.
01L/11L maxima are 0.124%/0.0881%. These are example-specific validation
results, not universal error bounds for the constant-medium approximation.

## Historical numerical-profile measurements

Measured 2026-10-29 using the RuO2/TiO2 four-model comparison: 20 keV,
180 points on 00L through L=3, equal local incidence/exit angles, s
polarization, W=-6, and alpha=0.15 (rough) or alpha=0 (flat). Wet cases use
3 Angstrom onset, 0.60 Angstrom initial width, molecular O+2H factors,
and onset-centred adsorption. Step A0 is zero or one molecule per surface
cell. Layered water has 2.8 Angstrom spacing, 0.40 Angstrom cumulative
broadening and relative first-layer weight one.

## Method and reproducibility

`examples/CTR/benchmark_water_DWBA.py` profiles dry, step, step+1ML and
layered water on both surfaces. `--include-erf` adds numerical erf optical
references for the two step models. It times preparation and evaluation
separately, records warm calls and fit-like changes, and exports complex
amplitudes and independent cProfile reports. The 01L/11L measurements follow
00L and may reuse fields. Model construction/import time is excluded.

```powershell
python -m examples.CTR.benchmark_water_DWBA --output tmp/water_profile --include-erf --profile
```

Wall times varied substantially during the work, including the unchanged
dry control. Therefore the table below uses **interleaved original/optimized
kernels in the same process**, alternating execution order between models.
BLAS/OpenMP thread environment variables were one for that benchmark process;
no library-wide or system-wide thread setting was changed. First-call times
are single observations; warm times are medians of two prepared calls.
The machine was not otherwise isolated, so small changes are noise rather
than evidence of speedup. These absolute times should not be compared with
older notebook timings made under different load/thread conditions.

Raw evidence is in `tmp/water_perf_before/` (all 36 rod/model combinations
and 12 cProfile reports), `tmp/water_perf_after_initial/` (first trial), and
`tmp/water_perf_paired/` (paired timings, complex amplitudes and copies of
the evaluated kernel definitions). The temporary pairing driver reconstructs
the original changed functions and verifies them against saved original
amplitudes before comparing the optimized functions.

## Paired 00L results

All times are seconds. "Repeated" means unchanged water and prepared optics;
atoms are still evaluated on every call.

| Surface | Model / optical reference | First before | First after | Repeated before | Repeated after |
|---|---|---:|---:|---:|---:|
| Rough | Dry | 0.35 | 0.28 | 0.206 | 0.348 |
| Rough | Step / tanh | 18.35 | 14.25 | 0.759 | 0.876 |
| Rough | Step + 1 ML / tanh | 24.37 | 25.05 | 0.705 | 0.604 |
| Rough | Layered / sharp | 12.25 | 10.52 | 12.945 | 0.201 |
| Rough | Step / erf | 5.98 | 5.35 | 0.656 | 0.936 |
| Rough | Step + 1 ML / erf | 24.55 | 9.84 | 0.928 | 0.613 |
| Flat | Dry | 0.15 | 0.15 | 0.099 | 0.103 |
| Flat | Step / tanh | 2.99 | 2.81 | 0.563 | 0.646 |
| Flat | Step + 1 ML / tanh | 3.13 | 2.94 | 0.430 | 0.463 |
| Flat | Layered / sharp | 0.87 | 0.80 | 1.235 | 0.105 |
| Flat | Step / erf | 2.66 | 2.72 | 0.635 | 0.563 |
| Flat | Step + 1 ML / erf | 7.29 | 2.44 | 0.468 | 0.398 |

The robust improvements are roughly **64x/12x for repeated rough/flat
layered water**, and **2.5x/3x for first rough/flat adsorbed erf calculations**.
Plain rough tanh improves modestly. No first-call improvement is established
for the adsorbed tanh case; dry kernels were not changed.

Changing the water profile invalidates its contrast integral. For example,
the optimized rough layered calculation still took about 11 seconds for a
changed initial width against frozen optics in this run. Unchanged-call
speedups must not be presented as speedups for fitting water parameters.
The reusable benchmark applies two distinct width updates when measuring
frozen versus automatic optics, preventing the second measurement from
reusing the first contrast. The original paired run's auxiliary automatic
timings followed the same update and are excluded here; in particular the
layered equimolar reference does not depend on the initial layer width.

## Bottlenecks and retained changes

1. **Layered-water contrast.** The original rough profile performed about
   29,000 Gaussian density evaluations for one specular scan. It rebuilt
   molecular factors and layer arrays for individual quadrature intervals.
   Adjacent intervals in one optical medium now share bounded batches.
   All original Gauss nodes, weights, boundaries and convergence tolerances
   remain in use. A four-entry cache reuses converged contrast when fields
   and physical water are unchanged; atomic coordinates/crystallinity remain
   live. Water settings, factors, widths, coverage, placements, weights,
   optical profile, angles, polarization and reference area participate in
   cache invalidation. An ambient mismatch still raises before cache lookup.

2. **Tanh internal fields.** Plain rough tanh originally spent most of its
   first evaluation sampling Chebyshev expansions, not solving the boundary
   matching problem. The polynomial basis is now built once per coordinate
   and contracted across field/angle coefficients. A non-BLAS contraction
   avoids making many small segments depend on worker-thread scheduling.
   A BLAS prototype had unstable end-to-end timing and was not retained.
   The analytic solution, initial expansion order 16, tail/error checks and
   exact analytic fallback are unchanged. An order-8 prototype was also
   discarded because it added work when the stricter tail check failed.

3. **Erf optical fields.** ODE right-hand sides repeatedly request scalar
   densities. Positive-width scalar samples now vectorize molecular factor
   terms and terrace weights; vector-grid and zero-width paths retain their
   previous formulas. This especially helps the adsorbed numerical reference.
   Solver tolerances and integration method are unchanged.

4. **Dry and non-specular rods.** Dry time is largely the native atomic
   kernel and model validation. Water's direct non-specular amplitude is
   zero. Those paths are not the main water bottleneck and were not broadly
   refactored.

## Accuracy and remaining work

Across the 12 paired specular comparisons, maximum absolute differences were
4.86e-12 in total complex reflection amplitude, 1.39e-10 electrons in F_h,
and 4.89e-12 in reference reflection. The initial batching/sampling trial
also checked all 36 combinations of models and 00L/01L/11L against the saved
original complex arrays. Comparison is performed on amplitudes, not only
intensities, so phase errors cannot hide.

Focused tests cover cache reuse with live atomic changes, water-width
invalidation, bounded eviction, ambient rejection, both polarizations,
scalar/vector density equivalence including legacy forms, independent
sampled-field quadrature and analytic/numerical field comparisons.

Validation: 112 focused water/continuous tests passed during development.
On the final source, 101 DWBA/continuous/migration tests and a further 28
targeted water tests passed. Ruff passed across `orgui` and the benchmark;
Sphinx built successfully with the two existing documentation warnings.

The main remaining target is **rough adsorbed tanh**, which constructs
hundreds of finite analytic segments. A representative profile had 745
nonconstant segments; fundamental-matrix/hypergeometric evaluation dominates.
Batching segment basis construction or reusing angle-dependent constants
is a better next target than further polynomial sampling tweaks. Reducing
the number of segments changes the optical approximation and needs its own
accuracy study. This work does not establish a universal total-speed
advantage of tanh over erf.

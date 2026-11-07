# Finite-film dissolution performance assessment

Date: 2026-11-07. Scientific implementation: `898bca1`, including final DWBA
master `5fc381b`. This assessment changes no runtime scientific code or fits.

## Measurement protocol

The supplied C1 loader verifies the SHA256 of both `0056_C1.xtal` and its
companion JSON. It restores the reference cell, energy, native incoherence,
atomic occupancies and Debye-Waller factors. The benchmark uses one identity
orientation domain, incidence angle 0.34 degrees, 512 points along
`(h, k, l) = (1, -1, 2.4 ... 2.8)`, s incident/outgoing polarization and no
extra bulk damping. Numerical-library thread counts are set to one. Loading
is excluded from evaluation timings; standalone evaluation includes its own
batch entry and cleanup. Three fresh-model timings are reported as medians
with ranges. Batched repeat/channel/rod timings are single observations.

These are per-call costs, not a complete 5,596-observation, four-domain,
resolution-convolved fitting objective. The two ensemble models have different
physical support policies, so their timing ratio is a workload comparison, not
proof of physical equivalence. C1 has 55 nonzero dissolution states spanning
indices -49 through +5; the legacy retained array has 32 entries, 31 nonzero.
No small probabilities were removed to improve throughput.

Fresh repetitions and the unchanged in-batch repeat must produce bitwise
identical reflectivities within each model. Work-count instrumentation runs
separately from wall-clock timings. Profiling times include profiler overhead;
cumulative times overlap and must not be added. Memory cases run in separate
fresh processes and distinguish resident memory from object lifetime.

Reproduce the timing/profiling assessment from the checkout:

```powershell
python -m benchmarks.benchmark_dissolution_performance `
  --request-dir F:/fuchs-4355-a/processed/LaNiO3/CTR/results/0056/dwba_ensemble/installed_api_verification `
  --output output/dissolution/performance_c1.json
```

Use a new output filename for every run. Individual sections can be repeated
with `--section fix`, `models`, `profile` or `scaling`. For separate memory
processes use `--section memory --kind legacy|dissolution --policy
shared|state_specific`. The raw JSON contains samples, work counts, profiles,
native binary hash and host information. It is kept under ignored `output/`;
the source benchmark and this record are versioned.

## Optical-position correction

The benchmark extracts the actual pre-fix combining functions from `898bca1^`
and times old/new implementations on the same captured C1 arrays, alternating
their order over five repetitions. Each helper sample covers ten repetitions
of the complete captured workload and is divided by ten. The public-API sample
uses the same arrays as a standalone helper benchmark; it is not an additional
cost incurred by the tagged DWBA path.

Tagged combining increases from 1.795 to 3.647 ms per captured workload; public
combining increases from 1.215 to 3.351 ms. Native end-to-end 512-ray medians are
0.545 s before and 0.558 s after, with ranges 0.361–0.662 s and 0.360–0.578 s.
Host variation exceeds that median difference, so this does not establish a
precise end-to-end percentage penalty. The helper overhead is measurable and
small in absolute terms. The scientific correction intentionally changes C1
native reflectivity here by up to `1.57e-4` relative; old/new outputs are not
required to agree. Repeated outputs within each corrected model are checked.

## C1 evaluation costs

Times are seconds on this Windows/Python 3.14.9 host (AMD family 26, model 36).

| DWBA model | Fresh median (range) | Batch first call | Identical batch repeat | New exit polarization | New rod |
| --- | --- | --- | --- | --- | --- |
| Legacy, shared | 0.596 (0.572–0.683) | 0.640 | 0.348 | 0.406 | 0.419 |
| Dissolution, shared | 14.316 (13.896–14.590) | 14.834 | 10.744 | 14.134 | 11.759 |
| Legacy, state-specific | 8.669 (8.558–9.521) | 8.398 | 5.974 | 6.692 | 5.120 |
| Dissolution, state-specific | 11.824 (11.797–12.205) | 12.483 | 9.419 | 11.246 | 10.159 |

Batch entry costs another 0.005–0.008 s for legacy or 0.047–0.062 s for
dissolution. Compared within each policy, fresh dissolution is about 24.0x
slower for shared fields and 1.36x for state-specific fields. The large shared
gap reflects the legacy model's optimized surface-height correction basis,
which dissolution does not use. Changing the field policy to obtain a faster
calculation changes the optical approximation and is not a performance fix.

For 512-point kinematic `F2`, fresh medians are 0.114 s for legacy (0.087–0.126)
and 0.729 s for dissolution (0.676–0.773): about 6.4x slower. Dissolution
evaluates complete native crystals for each state, including common material.

Measured build counts in a separate two-call, 32-point immutable batch:

| Model | First snapshots / packings | Repeated snapshots / packings |
| --- | --- | --- |
| Legacy, shared | 1 / 2 | 0 / 0 |
| Dissolution, shared | 111 / 56 | 110 / 55 |
| Legacy, state-specific | 32 / 32 | 0 / 0 |
| Dissolution, state-specific | 56 / 56 | 55 / 55 |

The dissolution wrapper retains state objects but does not enter each private
state's native DWBA batch context. Snapshot and packing reuse therefore never
activate for those states. Shared evaluation also snapshots every state for
the union grid, then again for scattering. It rebuilds that union/regridding
for each channel or rod and prepares the common optical fields separately on
each private state. Master's general height-domain path already demonstrates
batch-managed state snapshots, shared preparation and explicit cache cleanup.

The instrumented shared dissolution call takes 18.5 s: native record scattering
accounts for 7.61 s, optical snapshots for 5.42 s cumulative, and field retrieval
and construction for 3.13 s cumulative. Those are measured hot paths, not
promised optimization speedups. Python profiling overhead makes these numbers
inappropriate substitutes for the uninstrumented timings above.

The state-specific profiled call totals 14.3 s, including 6.70 s in native
record scattering, 2.58 s in optical snapshots, 2.90 s in field work, and
1.53 s in atom packing. Both paths therefore have costs beyond Python setup;
eliminating redundant snapshots alone will not make them as fast as the
legacy shared correction basis.

Single-call dissolution ray-count probes (same 55 states):

| Rays | Shared, seconds | State-specific, seconds |
| --- | --- | --- |
| 32 | 4.777 | 4.653 |
| 512 | 14.316 | 11.824 |
| 2,048 | 45.113 | 42.525 |

The 512-point row is the three-call median above; the other rows are single
fresh calls. There is a substantial setup floor and a substantial ray-dependent
cost. These measurements do not justify linear extrapolation to the full fit.
The sequential assessment process peaked at 3,213 MiB across all workloads;
that is not a per-call memory measurement.

## Memory and private-state lifetime

Separate fresh processes evaluated one 512-ray batch. RSS values include the
Python runtime and imported libraries; each baseline was about 132 MiB.

| Model | Inside batch, MiB | After batch exit, MiB | After dropping result/model and collecting, MiB | Process peak, MiB |
| --- | --- | --- | --- | --- |
| Legacy, shared | 159 | 159 | 142 | 174 |
| Dissolution, shared | 578 | 577 | 172 | 606 |
| Legacy, state-specific | 576 | 192 | 177 | 592 |
| Dissolution, state-specific | 513 | 512 | 178 | 541 |

Weak references show that all 55 private dissolution states remain alive after
batch exit and disappear after explicit cyclic collection. Legacy
state-specific retains 31 private state objects until collection, but its
explicit cleanup releases most of their large caches at batch exit. This is
delayed reclamation, not evidence of an uncollectable leak. The dissolution
path should adopt the same cache-release lifecycle rather than relying on
frequent global `gc.collect()` calls in fitting code. Allocators may retain
some pages even after objects are collected, so RSS is not an object-size sum.

## Optimization priorities

1. Put each private dissolution state's DWBA object into the native immutable
   batch lifecycle, as the master height-domain implementation does. Reuse its
   optical snapshot, packing and compatible native caches across channels and
   rods. Explicitly release private field tables at batch exit, including errors.
2. Compute the shared union grid and regridded state atomic models once per
   immutable batch. Keep geometry- and polarization-dependent field preparation
   separate from that invariant material assembly, preparing each common
   geometry/channel field once and reusing it for all compatible states.
3. After measuring these changes, assess common-material or incremental-state
   acceleration. Kinematic common substrate amplitudes can be reused; DWBA
   reuse must respect each state's optical field and reference subtraction.
   The native surface-only correction basis cannot be assumed valid across
   mixed-interface dissolution or exhausted states.
4. Only then consider reducing Python coalescing overhead. The current exact
   nonzero rule must continue to retain signed shares and empty media.

Validation for any optimization must compare uncached and cached native
amplitudes, coherent/incoherent endpoints and mixed intensities for both field
policies; include s/p channels, different rods, changed geometry and model
parameters outside a batch, failed candidates, complete exhaustion, and the
optical-support regression. Add deterministic work-count and cache-lifetime
checks. Preserve existing numerical budgets and material inventory; skipping
rare states requires a separately justified amplitude/moment error bound.

## Kinematic root-cause follow-up

The kinematic regression is an evaluation-strategy issue in the dissolution
wrapper. `FiniteFilmDissolution._iter_states` calls `path.state(index)` and
`state.crystal.F(...)` for every nonzero state. Each call reconstructs and
evaluates the entire retained crystal. It does not reuse common material as
the legacy `KinematicIncoherentContext`/surface correction path does.

On C1, one call performs 55 complete crystal evaluations, 55 identical bulk
evaluations, and 832 native unit-cell evaluations. The legacy profile performs
one complete crystal evaluation, one bulk evaluation and 123 unit-cell
evaluations. Construction invokes `state()` 56 times: the 55 ensemble states
plus a nominal state used while establishing the path's initial inventory.
The lazy coherent context is not causing an additional unused full evaluation.

The follow-up profiled dissolution call totals 0.839 s:

| Work | Time, seconds | Count |
| --- | --- | --- |
| Native unit-cell scattering | 0.525 self | 832 |
| Native bulk scattering | 0.071 self | 55 |
| Complete state construction | 0.151 cumulative | 56 |
| Path/support snapshot | 0.041 cumulative | 1 |

Cumulative construction/snapshot times overlap slightly and should not be
summed. No optical-profile construction or electric-field calculation occurs:
the timed diagnostic replaces those entry points with exceptions and still
completes successfully. DWBA object/geometry bookkeeping exists in generic
state construction, but it is not the expensive field-solving path.

A temporary C1-only memoization control hashes exact native scattering inputs
within one evaluation and reuses the original native function's result. It
does not replace any scattering formula, skip states, change weights, or keep
anything across fit proposals. It finds only 116 distinct unit-cell evaluations
among the 832 calls (716 redundant calls, 86%), and one distinct bulk evaluation
among 55. All five paired trials produce bitwise-identical `F2` arrays.

| Evaluation | Median, seconds | Range, seconds |
| --- | --- | --- |
| Current kinematic dissolution | 0.713 | 0.471–0.775 |
| Temporary exact-input memoization | 0.197 | 0.177–0.208 |

This is a measured 3.63x reduction before eliminating repeated state
construction. An initial exploratory cached call took 0.091 s, but the paired
five-trial result is the representative measurement. The cache control still
constructs all states, hashes arrays, runs Python dispatch and sums complete
crystal amplitudes. It is not a proposed production-wide UnitCell cache; the
diagnostic explicitly requires native C1 cells without custom form factors.

For kinematic optimization, prioritize one bulk/common-material evaluation,
native per-layer/domain amplitudes reused across retained-height states, and
cumulative sums over removable material intervals plus each state's exact
termination correction. Keep complete crystal construction for independent
oracles, inventories, and optical evaluation. Coherent and incoherent mixing
must continue to use complete amplitudes, with reference-area scaling, outer
domains, attenuation, signed records, and exposure policy preserved. This
follow-up diagnoses and measures the issue; runtime implementation is unchanged.

Reproduce with the benchmark's `--section kinematic` option and a new output
filename. The run stores paired timings, exact cache misses and the profile;
the original measurements are in `output/dissolution/kinematic_redundancy.json`.

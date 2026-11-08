# Finite-film dissolution performance assessment

Date: 2026-11-07. Original assessment: `898bca1`, including final DWBA
master `5fc381b`. Historical measurements below predate the production
optimizations recorded in the final section. No fits or input data were changed.

The staged follow-up is documented in the
[kinematic and DWBA performance implementation plan](finite_film_dissolution_performance_plan.md).

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

## DWBA root-cause follow-up

`benchmarks/audit_dissolution_dwba_reuse.py` adds temporary native reuse
controls, without editing the runtime implementation or C++ kernel:

- `state_batches` enters every private state's existing native batch context,
  shares the immutable cell-profile cache, and clears private persistent caches
  at exit. The production wrapper enters only the averaged crystal's batch.
- `shared_prepared` additionally reuses the common prepared fields, replacing
  only state identity and regridded atomic/reference shares. It asserts exact
  equality of union-grid boundaries and values before reuse. It retains the
  original all-state union grid and all probabilities.

First-call work counts on C1 (32-ray diagnostic, same 55 nonzero states):

| Shared-field control | Snapshot builds | Packings | Preparations | Field solves | Distinct field-solver inputs | Native scattering calls |
| --- | --- | --- | --- | --- | --- | --- |
| Current implementation | 111 | 56 | 57 | 114 | 4 | 56 |
| Native state batches | 56 | 56 | 57 | 114 | 4 | 56 |
| Batches + common preparation | 56 | 56 | 2 | 4 | 4 | 56 |

On an identical second call, current code rebuilds 110 snapshots and 55
packings; both controls rebuild neither. Common preparation reduces the second
call to one explicit preparation. The remaining first-call pair prepares an
initial averaged reference and then the refined union reference; a future
implementation can order assembly before field preparation to avoid that
initial unused reference preparation.

The state-specific case is different: its 112 solver calls have 112 distinct
inputs. Those fields cannot be shared merely because atomic material overlaps.
Its redundant repeated snapshots/packings still disappear with native batches.

Paired 512-ray timings, three alternating-order trials for fresh medians:

| Policy/control | Fresh median, seconds | Identical in-batch repeat, seconds |
| --- | --- | --- |
| Shared, current | 7.438 | 5.984 |
| Shared, native state batches | 6.535 | 4.020 |
| Shared, batches + common preparation | 4.904 | 3.924 |
| State-specific, current | 6.298 | 4.959 |
| State-specific, native state batches | 5.968 | 3.762 |

Repeat timings are single observations. Batch entry is about 0.023–0.028 s for
current code and 0.047–0.051 s for the controls, which eagerly create the private
states. These paired-run timings differ from the earlier host timings; compare
within a run rather than interpreting that change as a code speedup. C1's mixed,
coherent and incoherent outputs are bitwise unchanged in these controls. A
further 24 comparisons cover s/s, s/p, p/s, p/p, mixed specular/non-specular
coordinates, damping, interface removal and near-exhaustion states; the maximum
relative change is `1.78e-15` (acceptance `rtol=1e-10`, `atol=1e-28`).

### Why native batching does not remove the main 512-ray cost

The final master dispatcher in `cpp/CTRdwba_cpp.hpp` selects its bounded native
finite-scattering cache only for calls of at most **128 rays**, with a 1 MiB
limit per cache entry. This is unrelated to any distribution state-count limit.
At 512 rays the uncached kernel is selected even with a native batch active.
Measured shared-preparation control counts confirm it:

| Rays | First / repeated finite-cache hits | First / repeated normal-phase exponentials |
| --- | --- | --- |
| 32 | 0 / 55 | 2,350,848 / 1,067,520 |
| 512 | 0 / 0 | 37,613,568 / 37,613,568 |

The 512-ray legacy shared-height implementation performs two native scattering
calls and 2,134,528 normal-phase exponentials, versus 56 calls and 37,613,568 for
dissolution even after these controls. Thus Python caching alone cannot close
the gap. A material-aware common/per-layer amplitude basis on the shared fields
is the larger algorithmic opportunity. The existing surface-only basis cannot
be applied blindly across mixed-interface removal and exhaustion. Neither the
128-ray threshold nor memory bound was changed in this assessment.

### Unused averaged-crystal scattering in state-specific mode

The wrapper always calls `average.evaluate_prepared`. In state-specific mode
none of that result's scattered amplitudes or intensities are consumed: the
coherent intensity is computed from the state amplitude sum, and only the
averaged result's geometry-derived `amplitude_prefactor` is used at return.
The `state_batches_no_average` diagnostic substitutes an object exposing only
the native prefactor property, so any other result access would fail. C1 is
bitwise unchanged and eight boundary/polarization checks agree within
`1.45e-15`; packings and scattering calls drop from 56 to 55. This control
retains the averaged preparation and phase reference, and does not assume that
the distinct state fields can be merged. Its one-pair timing is not used to
claim a stable percentage speedup.

### Reproduction and implementation boundary

At the historical `0a38763` revision, run
`python -m benchmarks.audit_dissolution_dwba_reuse` with `--request-dir`
as above, `--policy shared` or `state_specific`, and a new `--output` JSON.
Optional `--points` and `--repeats` control the timing workload. The helper
`_counts(policy, mode, points)` provides the 32/512-ray native-work probes;
`_control("state_batches_no_average")` isolates the unused averaged scattering.
Raw outputs are `output/dissolution/dwba_reuse_shared.json`,
`dwba_reuse_state_specific.json`, `dwba_reuse_native_limits.json` and
`dwba_reuse_legacy_native_work.json`. Native state/field reuse and removal of
unused work are diagnostic controls only, not yet production changes. Any
implementation must also validate changed rods/channels within a batch,
parameter updates between batches, exception cleanup and cache lifetime.

## First implementation stages

Stages 1 and 2 of the performance plan are implemented on the dissolution
worktree. These are numerical work-reuse changes. Material selection, physical
exhaustion, numerical support, probabilities and all optical positions remain
unchanged. No data or scales were fitted. The native extension is unchanged:
SHA256 `4373ef0907b7178b4659173040f926cfe635bd0d853d7fd2e550cbfa1a96062f`.

### Implemented changes

- Kinematic evaluation validates material selections without building a full
  crystal for every height. It evaluates the protected bulk/material once,
  accumulates removable-interval amplitudes and growth prefixes, and adds each
  exact termination/reference correction. All scattering and area/domain
  conversion still use native UnitCell/FrozenMaterial/crystal evaluators.
  Unknown cell subclasses or custom evaluators keep the complete-state path.
- DWBA enters the native batch lifecycle for every private state. One shared
  cell-profile cache serves the immutable ensemble. Shared-reference mode
  constructs the all-state union grid and rebins each decomposition once;
  all states use the same prepared wavefields. The last shared ray preparation
  is reused, with bounded native preparation caches for state-specific mode.
- State-specific mode no longer evaluates the unused averaged atomic
  scattering. A private geometry-only helper supplies the existing amplitude
  prefactor. Distinct state wavefields, the averaged phase frame and native
  reference subtraction are retained.
- All private average/state caches are cleared at outer batch exit, including
  after errors. No cache survives across standalone fit proposals.

### Paired C1 measurements

`benchmarks/benchmark_dissolution_optimized.py` compares the optimized path with
frozen native complete-state oracles. It alternates timing order, records source,
native and C1 hashes, and counts native work separately from timing. Loading is
excluded; fresh calls include snapshot construction, preparation and cleanup.
Inputs are the same 512-point C1 rod and single-thread policy described above.

| Mode | Complete-state median | Optimized median | Ratio |
| --- | ---: | ---: | ---: |
| Kinematic fresh, five pairs | 0.36686 s | 0.04803 s | 7.64x |
| Shared DWBA fresh, three pairs | 7.64163 s | 5.12783 s | 1.49x |
| Shared DWBA repeat | 6.44930 s | 4.07152 s | 1.58x |
| State-specific DWBA fresh, three pairs | 6.43319 s | 6.05091 s | 1.06x |
| State-specific DWBA repeat | 5.45972 s | 3.78691 s | 1.44x |

Kinematic ranges were 0.36326–0.43288 s complete and 0.04724–0.05216 s optimized.
A separate same-host five-trial comparison measured legacy at 0.05024 s median,
versus optimized dissolution at 0.04809 s. Its complete-state median was
0.36207 s; compare within each run rather than treating the earlier 0.7 s
observation as a fixed machine speed.

DWBA timings are preliminary Stage 2 measurements, not final acceptance.
Shared fresh ranged 7.58417–7.86855 s complete and 5.08671–9.79295 s optimized;
the first optimized call was the high outlier. Repeat ranges were
6.17872–6.58582 s and 3.99948–4.18506 s. State-specific fresh ranged
6.38010–11.94844 s complete and 5.89238–6.09333 s optimized; repeat ranges were
5.30776–10.53486 s and 3.71965–3.89172 s. Retain these outliers in the report;
final acceptance needs the planned five-pair trials after the later stages.

C1 kinematic intensity differences are at most `1.44e-14` relative. Both DWBA
policies produce bitwise-identical reflectivity, coherent/incoherent intensities
and squared prefactors in these runs. Repeat results are bitwise equal to fresh
results. Tests also compare the optical grid values and boundaries exactly.

| Native work, C1 | Before | After |
| --- | ---: | ---: |
| Kinematic bulk / unit-cell calls | 55 / 832 | 1 / 116 |
| Shared DWBA first snapshot / packing builds | 111 / 56 | 56 / 56 |
| Shared DWBA first field solves | 114 | 4 |
| Shared DWBA repeat snapshot / packing / field solves | repeated | 0 / 0 / 0 |
| State-specific first packing / scattering calls | 56 / 56 | 55 / 55 |
| State-specific distinct first field solves | 112 | 112 |
| State-specific repeat snapshot / packing / field solves | repeated | 0 / 0 / 0 |

The shared path still performs 56 native scattering calls and 37,613,568 normal
exponentials at 512 rays. State-specific still performs 55 scattering calls and
36,331,008 normal exponentials. Native finite-cache hits remain zero at this
ray count. This is why Stage 2 does not meet the final shared-DWBA target.

The focused tests cover complete-state amplitude comparisons, unequal lateral
areas/reference frames, mixed interfaces, explicit finite substrate, growth and
exhaustion, common outer domains, all polarization pairs, mixing endpoints,
scalar/broadcast HKL, custom-cell fallback, native work counts and nested/error
batch cleanup. Existing live-update and optical-position regressions remain
part of the wider validation.

Reproduce with `python -m benchmarks.benchmark_dissolution_optimized`, supplying
the hash-checking `--request-dir` used above, `--mode kinematic`, `shared` or
`state_specific`, `--repeats` and a new `--output`. Raw files are
`output/dissolution/optimized-kinematic.json`, `optimized-shared.json` and
`optimized-state-specific.json`. Source hashes in those files describe the
measured implementation, before its commit was created.

### Full kinematic replay and regression results

The fixed-input C1 replay covers 5,596 non-specular observations, 26,042 original
resolution nodes and all four recorded domain rotations. It averages domain
intensities, applies the original resolution matrix, takes square roots and
retains the original rod scales. No refitting or support changes were made.

| Full replay | Evaluation time | chi-squared / N |
| --- | ---: | ---: |
| Complete-state dissolution oracle | 66.20480 s | 8.394242357554576 |
| Optimized dissolution | 5.00442 s | 8.394242357554575 |
| Legacy model, different physical boundary | 9.04367 s | 8.394242571229753 |

This single replay is about 13.2x faster than the complete-state calculation;
it is an objective-level check, not a five-pair performance acceptance study.
Maximum node-intensity relative difference is `3.22e-14`; maximum residual
change is `2.00e-14` sigma. Tightening the tail tolerances produced no observable
change in this replay. The legacy saved amplitude also agrees within `1.14e-13`
absolute, confirming the original grid/domain/scale reconstruction.

Run `python -m benchmarks.validate_dissolution_c1 --root <LaNiO3-directory>
--compare-complete-states --output <new-report.json>` to reproduce. Raw report:
`output/dissolution/optimized-full-c1.json`. Loading and material-statistics
reporting are excluded from the evaluation times; the four native evaluations
and resolution/domain averaging are included.

Validation completed:

- Full `test_CTR*.py` suite: 1,351 passed, 3 skipped, 193 passing subtests.
- Final performance-regression module: 50 passed, including four additional
  cases for changing channels/ray blocks and explicit `cache=False` within a
  batch. The explicit bypass still disables native preparation caching.
- Ruff passes for the scientific-core tree and both changed benchmark drivers;
  `git diff --check` passes. Release notes were regenerated from the changelog.
- No C++ changes were made, so the previously verified native extension remains
  the evaluated binary. Isolated-wheel verification belongs to the final stage.

### Remaining stages

The shared-field material amplitude basis (Stage 3), further field-independent
state-specific reuse (Stage 4), large-ray native cache/block policy (Stage 5)
and final isolated-wheel/full-DWBA-objective validation (Stage 6) are still open.
Do not label the entire performance plan complete or infer full-objective DWBA
throughput from these 512-ray measurements.

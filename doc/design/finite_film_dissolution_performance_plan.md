# Kinematic and DWBA dissolution performance implementation plan

Status: planned; diagnostic controls exist, production optimization has not
started. Date: 2026-11-07. Baseline: `0a38763`, including authoritative DWBA
master `5fc381b` and optical-position correction `898bca1`.

Remove repeated full-state work while preserving the implemented finite-film
model and native scattering. Implement kinematic amplitude reuse first, then
DWBA lifecycle/preparation reuse, then the larger shared-field amplitude basis.
Treat state-specific field costs separately. Finish with measured large-scan
and complete-objective validation before declaring the fitting path ready.

The evidence and reproducible diagnostics are in
[the performance assessment](finite_film_dissolution_performance.md).
The physical contract remains in
[the dissolution implementation plan](finite_film_dissolution_implementation_plan.md)
and [its validation record](finite_film_dissolution_validation.md).

## Scope and invariants

This is a performance change, not a new dissolution or optical model.

- Preserve the selected finite owners, mixed-interface populations, signed
  reference records, exposure policy, absolute coordinates and lattice/domain
  conventions. Dissolvable substrate still requires an explicit finite Film.
- Preserve exhaustion probability, all retained nonzero probabilities, tail
  certificates, normalization and realized moments. Do not shorten support,
  relax tolerances or enlarge the Film to obtain speed.
- Preserve the corrected optical positions, all-state shared union grid,
  intentionally empty media, native reference subtraction and common phase
  convention. Never choose a different optical policy automatically.
- Keep native unit-cell/bulk scattering and the native DWBA kernel as the
  scientific implementation. Reorganize their inputs and reuse their outputs;
  do not reproduce their scattering formulas in a separate Python backend.
- Kinematic outputs remain squared amplitudes per reference lateral cell.
  DWBA outputs remain the existing dimensionless observables. Preserve area
  normalization, coherent outer domains, polarization, damping and broadcasting.
- Preserve both mixture endpoints and their interpolation. In particular,
  dissolution's coherent endpoint is the bounded ensemble mean amplitude;
  the wrapped legacy crystal's amplitude is not a valid shortcut for it.
- Preserve rejection of unsupported topology, inconsistent geometry, invalid
  inventory and resource-limit violations. Acceleration must not bypass checks.
- No data refit, fit replacement, main-checkout edit or master merge is included.

All caches are evaluation-local or owned by an explicit immutable native batch.
Fit proposals outside a batch rebuild the necessary state. Initially add no
persistent cross-proposal cache and no general global UnitCell memoization.

## Implementation sequence

| Stage | Deliverable | Main files | Exit condition |
| --- | --- | --- | --- |
| 0 | Frozen native oracles and work counters | dissolution tests and benchmarks | Correctness and performance baselines reproducible |
| 1 | Material selection plan and kinematic amplitude stream | `CTRdissolution.py`, `CTRincoherent.py` | One bulk evaluation; no complete crystal per height in the fast path |
| 2 | Native DWBA state lifecycle and one shared preparation | `_CTRdissolution_dwba.py`, narrowly `CTRdwba.py` | Repeated snapshots/packings eliminated; duplicate shared field solves removed |
| 3 | Shared-field material amplitude basis | `_CTRdissolution_dwba.py`, native record adapters | Native scattering work follows distinct material terms, not repeated complete states |
| 4 | State-specific field-independent reuse | DWBA preparation/packing and native helpers as needed | Distinct fields preserved; remaining reusable work removed |
| 5 | Bounded large-ray execution and memory cleanup | native dispatcher/cache only if justified | Stable memory and throughput across ray/block sizes |
| 6 | Full validation and performance acceptance | tests, C1 drivers, documentation | Scientific gates and measured performance goals satisfied |

Keep these as small Conventional Commits with `phys` scope for implementation
and scientific regression changes. Stage 1 is useful independently of DWBA.
Stage 2 can land independently of the shared-field basis. Each stage must retain
a correct native fallback while later stages are developed.

## Stage 0 — establish independent oracles

1. Retain the current complete-state kinematic evaluation as a private oracle:
   `DissolutionPath.state(index).crystal.F(...)`, followed by the current weighted
   coherent/incoherent sums. Keep complete states available for inventories,
   inspection and DWBA even when kinematic fitting no longer constructs them.
2. Retain current complete-state DWBA evaluation as the reference for optimized
   preparations/bases. Add an internal test-only selection mechanism, not a
   serialized model setting or a second user-facing physics mode.
3. Add work-count assertions for bulk calls, distinct native contributions,
   state construction, snapshot/packing builds, field solves and native record
   work. Count outside wall-clock benchmark regions.
4. Save exact source/native hashes and C1 input hashes with timings. Preserve the
   final master changes; check for new uncommitted work before implementation.
5. Preserve independent plain-Film, exposed-interface and explicit substrate-Film
   oracles. Comparing two paths that share a newly edited selector is not enough
   to validate that selector.

Baseline indicators: C1 kinematic evaluation makes 55 bulk and 832 unit-cell
calls; only one bulk and 116 complete unit-cell input sets are distinct. Shared
DWBA makes 114 field-solver calls for four distinct inputs and 56 native
scattering calls. These indicate redundant work, not fixed target counts for a
new per-layer decomposition.

## Stage 1 — kinematic material amplitudes

### Separate selection from complete-crystal construction

Refactor the internals of `DissolutionPath` to expose an immutable, private
material selection plan built once per evaluation. It must describe protected
records, ordered removable intervals, continued-growth intervals and each
height's exact termination additions/subtractions. Use the same selection
rules as `state()`; do not substitute a nominal-height-only model.

Compute exposure eligibility and inventory validity once while planning.
Use cumulative species inventories where appropriate, but still reject every
state that the existing full-state constructor rejects. Repeated construction
of SXRDCrystal/DWBA objects must not be needed to validate or select material.

Suggested private interfaces are an amplitude-plan object and
`_iter_kinematic_amplitudes(...)`; names may follow surrounding conventions.
Leave public `state()`, statistics, persistence and fitted parameters compatible.

### Evaluate native contributions once

1. Broadcast HKL through the existing native crystal convention.
2. Evaluate the protected semi-infinite bulk once, retaining its reference-area
   factor, lattice repeat, coherent domains and attenuation.
3. Evaluate the other protected material once. For each outer domain, transform
   HKL before calling native scattering; use the canonical component-domain
   scaling helpers rather than copying their scale-factor arithmetic.
4. Evaluate each required removable/growth interval contribution once with the
   native cell/domain evaluator. Reuse exact immutable cell/geometry terms
   within this request. Include occupancy, signed weights, crystallinity,
   reference transforms, area and domain identity in their meaning.
5. Walk retained heights in order, updating a running material amplitude by
   adding/removing interval contributions. Add that state's exact native
   termination correction, including the unreconstructed fallback where needed.
6. Stream the complete state amplitude into the existing probability-validated
   accumulator. Preserve `sum(p*A)`, `sum(p*abs(A)**2)` and the mixing fraction.

Do not repeatedly build a full state and hash all its arrays as the permanent
solution. The diagnostic memoization still incurs that overhead. Prefer stable
term identifiers inside the immutable request plan; do not retain object-id
keys across requests. Unknown cell subclasses or custom form-factor behavior
must use a documented correct fallback until their reuse contract is validated.

Bound temporary memory: avoid a dense states × terms × rays tensor. Stream
height prefixes and termination corrections; if needed, process ray blocks
without changing normalization or the phase/reference convention. Scattering
may exploit exact zeros, but optical bookkeeping records remain intact and no
small nonzero probability or material tail is treated as zero.

Stage gate: one bulk call per request; no complete-state construction in the
eligible amplitude loop; native contribution evaluations scale with distinct
required terms rather than the sum of terms across all states. All current
material/inventory/shape/mixture tests must still pass.

## Stage 2 — DWBA batches, fields and unused work

1. Use an exception-safe `ExitStack` for the averaged crystal and private state
   native batch contexts. Share compatible immutable cell-profile caches using
   the authoritative master's existing approach. Nested batches share one
   scope; private caches are released on normal exit and exceptions.
2. Retain each state's `(records, atomic_model, optical_profile)` once per
   batch. Retain compatible packing tables. Do not invoke an uncached snapshot
   once for the shared grid and again for scattering.
3. In shared mode, form the all-state union grid before field preparation.
   Preserve every representable split and the existing near-coincident-boundary
   rule. Regrid each atomic model once onto this grid. Cache the grid and those
   models for the immutable batch.
4. Prepare common fields once per distinct geometry/polarization request.
   Reuse the resulting preparation for each state by changing only the native
   state identities and atomic/reference model, as in the master's general
   shared-state implementation. Do not solve an intermediate unrefined
   averaged field merely to obtain geometry.
5. In state-specific mode, retain distinct optical profiles/fields. Reuse each
   state's compatible snapshot, packing and cached fields across calls;
   equality of material alone is not a field-cache key.
6. Remove unused averaged-crystal scattering in state-specific mode. Factor the
   existing native amplitude-prefactor property into a small private helper of
   the prepared geometry/reference, and have the existing property delegate to
   it. Preserve the averaged preparation, returned result metadata and phase
   origin. Do not ship the diagnostic's placeholder result object.
7. Explicitly clear private persistent field/preparation caches at batch exit,
   in addition to closing native batch caches. Do not add routine global
   garbage collection as the memory-management solution.

Stage gates for identical repeated requests: zero snapshot and packing builds,
zero new field solves, and no private field-cache growth after scope exit.
For the C1 shared first request, field solves should follow distinct common
incoming/outgoing fields, not 55 states; the target is two solves after avoiding
the intermediate averaged preparation. Respect native de-duplication when the
two requests are identical. State-specific mode must perform no unused
averaged-crystal native scattering.

## Stage 3 — shared-field material amplitude basis

Stage 2 alone leaves 56 full native scattering evaluations. Replace this repeated
work for eligible shared-field calculations using the Stage 1 material plan.

Pack protected material, unique removable/growth terms and termination terms
into a material-aware basis of native records. Evaluate common contributions
once on the unchanged shared fields and obtain per-term native amplitudes.
Use interval prefixes and signed termination corrections to reconstruct each
complete state, then apply the current shared-reference subtraction, common
phase and probability/coherence sums. Keep the averaged coherent endpoint
consistent with the current implementation and independent oracle.

Do not directly reuse the legacy surface-height basis across mixed interfaces:
it assumes a common Film/interface that the dissolution model can remove.
Owner/provenance, exact geometry and termination validity must define terms.
If native record output needs an internal extension, preserve the existing
kernel as the uncached reference; do not create another scattering formula.

Implement eligibility in two steps:

- First support fully ordered records and non-specular cases where the native
  planar-reference term vanishes. Fixed layer crystallinity must still weight
  the correct atomic terms. Use native specularity classification. The C1
  non-specular workload must use this fast path.
- For partially ordered specular records, retain each state's effective-density
  reference shares. These can depend on its optical partition even on common
  fields; they cannot be replaced by the averaged reference. Extend reference
  integration with native per-state shares only after independent validation.
  Until then, use the optimized complete-state path for ineligible requests.

An initial mixed specular/non-specular request may take the correct fallback;
later point partitioning must preserve order, shapes and coherent sums. Record
fast-path/fallback reasons in developer diagnostics. Do not hide implementation
errors behind a broad exception-to-fallback handler.

Stage gate for eligible C1 shared work: bounded native calls per ray block
(aim for one basis call plus at most one averaged/common call), rather than a
call per height. Phase-exponential work should follow unique material terms.
Verify reconstructed complex amplitudes and both intensity endpoints before
accepting a faster mixed intensity.

## Stage 4 — state-specific reuse without merging fields

Profile after Stage 2. Keep genuinely different field solves and propagation
phases separate. Their cost is part of the selected optical approximation.

Reuse static atom/geometry descriptions, domain transforms and compatible
packing layouts across states. Where profiling justifies it, extend native
precomputation of field-independent atomic factors and lateral contributions.
Trace every proposed key against native inputs: coordinates, energy, atomic
factors, displacement/occupancy, crystallinity, transforms and selected records.
Keep field-dependent normal phases, boundary integration, bulk attenuation and
state reference subtraction outside an incompatible cache entry.

Retain the original native path as an oracle and fallback. Do not share a bulk
DWBA amplitude across states with different fields. Do not claim that the
state-specific path should match shared-field throughput. Profile remaining
native work before adding C++ complexity.

## Stage 5 — dense scans and bounded memory

The current native finite cache has a 128-ray dispatch threshold and a 1 MiB
entry bound; these are unrelated to support-state limits. Do not simply raise
either value. First measure Stages 1–4 at 32, 127, 128, 129, 512 and 2,048 rays.

Compare an appropriate uncached dense path with bounded ray blocks scheduled
to reuse compatible fields/native sums across channels. Preserve the public
prepared result, global phase origin, point ordering, scalar/broadcast shape,
and intensity-before-resolution processing. A ray block must not cause a
different material support or optical grid. If changing native dispatch or
cache sizing is beneficial, document measured memory and fresh/repeated-call
tradeoffs and retain a bounded uncached fallback.

Measure first call, repeated calls, changed exit channel, changed rod and changed
parameters outside a batch. Test nested/failed batches. Use weak references and
cache-size counters to verify lifetime; RSS alone is allocator-dependent. No
unbounded retained states × channels × rods cache is acceptable.

## Stage 6 — correctness and performance acceptance

### Scientific checks

- Compare complex native state amplitudes where applicable, coherent and
  incoherent endpoints, mixed intensities and prefactors. Pure cache reuse
  should remain bitwise equal where the current diagnostics establish that.
  For reordered sums, retain existing fixture tolerances and inspect errors
  near zeros separately; do not replace them with a looser global tolerance.
- Keep the optical-support observable gate at `1e-9`. Keep the diagnostic
  DWBA reuse comparison at `rtol=1e-10`, `atol=1e-28` on its existing fixtures.
  Add scale-aware amplitude checks near zeros with an explicitly justified
  absolute budget, rather than dividing by vanishing amplitudes.
- Exercise full exhaustion, partial mixed-interface removal, explicit finite
  substrate, multiple Films, incomplete motifs, growth, unequal lateral areas,
  nonuniform spacing, signed references, tiny nonzero tails and empty media.
- Exercise all polarization pairs, common phase, specular/non-specular mixtures,
  near-Bragg points, scalar and broadcast inputs, rotations, coherent outer
  domains, layer order, energy, attenuation and invalid-candidate recovery.
- Compare after live changes to rates, offset, incoherent fraction, atom positions,
  occupancy, displacement factors, geometry and domain rotations. Include
  configuration/deepcopy/pickle/native-save round trips with no retained caches.
- Run focused dissolution/optical/native-reuse tests, then the full CTR suite.
  If shared native helpers change, broaden to their other scientific users.
  Fix local lint only; avoid unrelated refactoring.

Rebuild the extension after every C++ change before trusting tests. For final
package validation, build without in-place source binaries shadowing the fresh
extension, assert a single native archive entry and matching build hash, then
test an isolated wheel installation. Do not replace the normal installation.

### Performance goals

These are proposed engineering acceptance goals, not measured speedups. Use
paired same-host runs with identical inputs/native builds and at least five
alternating trials. Report medians, ranges, native work counts and memory.
Keep timing thresholds out of noisy unit-test CI; enforce deterministic work
counts there. Measure snapshot/preparation work as part of each fresh call.
Treat a paired median slowdown above 5% that exceeds observed repeat variability
as a regression. If variability masks the comparison, repeat under more stable
conditions rather than adjusting the scientific acceptance tolerances.

| Workload | Goal |
| --- | --- |
| C1 kinematic, 512 rays | At least 5x faster than current dissolution and within 1.5x of paired legacy throughput; aim for about 0.1 s on the assessment host. The 0.2 s memoization control is an intermediate result, not completion. |
| C1 shared DWBA, 512 rays | At least 4x faster fresh than current dissolution and within 3x of paired legacy shared throughput after the material basis; Stage 2's approximately 1.5x diagnostic gain is not completion. |
| State-specific DWBA | Eliminate redundant builds and unused average scattering; target at least 1.25x faster identical batch repeats, with no material fresh-call regression. Set further native targets from the Stage 4 profile rather than assuming fields can be shared. |
| Memory/lifetime | Release private field/native caches at batch exit; no monotonic retained-cache growth across repeated fit proposals or exceptions. |
| Full C1 objective | Demonstrate the speedup on all recorded domains, rods and resolution nodes, preserving original intensity/resolution ordering and fixed scales. No extrapolation from 512 rays alone. |

Validate fixed-parameter C1 predictions against the current dissolution oracle,
not against legacy physics as if they were identical models. Preserve canonical
data hashes and output files. Include broader retained supports up to supported
resource limits. If the performance goals are missed, identify the remaining
native work and leave that stage incomplete; never meet a target by changing
the physical model or numerical error budgets.

## Delivery

Update this plan's stage status and the performance assessment with measured
results as each stage lands. Update relevant API/topic documentation when
cache/batch behavior changes; add a concise changelog entry and regenerate
release notes for user-visible throughput changes. Retain the slow independent
oracles, work-count regressions and benchmark drivers. This plan itself does
not implement the optimizations or authorize a data refit.

# Complete finite film dissolution implementation plan

Status: implemented in `codex/finite-film-dissolution`; see
`finite_film_dissolution_validation.md` for executed checks and limits.
Date: 2026-11-07.

Implement an opt-in native model that permits complete dissolution of designated
finite Films, including their material represented inside EpitaxyInterface
components. Preserve the semi-infinite substrate. If substrate dissolution is
required, represent the removable substrate explicitly as a finite Film using
the substrate unit cell. Use native kinematic and DWBA scattering throughout.

This plan separates the physical extension from numerical support selection.
Neither removing the existing SkellamSurface guard nor increasing Film width
implements complete dissolution correctly. No existing fit is to be overwritten
or refitted as part of implementation validation.

## Scientific contract

The user established two requirements: complete film dissolution must be
representable, and dissolvable substrate must be an explicit finite Film. The
following boundary and interface policies are proposed implementation choices;
they must be explicit in configuration and validated before production fitting.

| Concern | Contract |
| --- | --- |
| Material ownership | A finite Film owns its generated interior and its assigned material contributions in adjacent interfaces. Ownership is by component identity, not element name or shared unit-cell object. |
| Protected substrate | The semi-infinite bulk never becomes removable through a surface parameter or an interface association. |
| Several Films | A front can traverse an explicitly selected, connected sequence of Films. It cannot jump across a protected Film to dissolve material below it. |
| Exhaustion | Proposed physical policy: accumulate all probability beyond the last removable state into the exhausted state. Do not discard or conditionally renormalize this probability. |
| Remaining substrate | Preserve existing lower-material occupancy and geometry. Removing upper material does not replenish missing lower material or restore a pristine substrate. |
| Coordinates | Preserve absolute positions, lattice conventions, lateral registry, layer cycles, strain, and component area factors when assembling states. |
| Exposed structure | Select a termination valid for the remaining material and thickness. The exhausted state cannot retain a translated LaNiO3 surface motif. |
| Conservation | For every finite owner and chemical species, initial inventory plus explicit additions equals remaining inventory plus removed inventory. An exhausted owner has zero remaining inventory within certified numerical support. |
| Compatibility | Existing SkellamSurface, PoissonSurface, SurfaceHeightDomains, files, and defaults retain their current behavior. The new physical policy is opt-in. |

Use Angstrom for absolute heights and displacements, dimensionless probabilities,
seconds and inverse seconds for the existing rate/time parameters, and integer
structural-state indices for front progression. Convert chemical inventories to
atoms per common lateral reference area before comparing components. Anomalous
scattering factors and Debye-Waller factors are not chemical atom counts.

## Baseline and motivation

The 0056/C1 sidecar and native reload give:

| Quantity | Value |
| --- | --- |
| Nominal saved Film width | 38 structural planes |
| Generated Film above interface support | 25 planes |
| Interface upper stacking boundary | 22.9615451034 Angstrom |
| Nominal Film surface boundary | 70.0521379260 Angstrom |
| Mean surface height | 64.1986672656 Angstrom |
| Parent mean displacement | -3.10755838349 planes |
| Raw Skellam variance | 7.32090470521 planes squared |
| Variance including fractional recentering | 7.48873437393 planes squared |
| RMS height | 5.15464027028 Angstrom |
| Incoherent fraction | 0.205425229351 |
| Asymmetry | -1 |
| Retained height interval | -25 to +6 |

For C1, probability below -25 is 2.3556437895e-10 and expected attempted removal
beyond that boundary is 3.0650543186e-10 plane equivalents per lateral column.
The immediate support ceiling is numerical. At fixed mean and asymmetry,
interval tightening with existing native bounds permits approximately 6.2559
Angstrom RMS; sharper tails with all existing fit audit criteria permit
approximately 6.4506 Angstrom. These are distribution-level estimates, not
validated scattering allowances or unconstrained fitted roughness values.

Reproducibility sources:

- `output/skellam_support_assessment/audit.py` and `audit.json` in this checkout
  contain calculations and input hashes. They are local evidence, not portable
  test dependencies.
- Canonical analysis root: `F:/fuchs-4355-a/processed/LaNiO3/CTR`.
- Analysis inputs: `settings/1V2/0056_C1.xtal`, `settings/1V2/0056_C1.json`,
  `results/0056/SUMMARY.html`, `fits/master/model_0053/skellam_nonspec_model.py`,
  and `fits/master/model_0055/wide_model.py` under that root.
- Related designs: [recentered Skellam](recentered_skellam_surface_implementation_plan.md),
  [multiple Film epitaxy](multifilm_epitaxy_api_implementation_plan.md), and
  [incoherent models](incoherent_ctr_models.md). Read current code alongside
  historical plans; the latter include superseded implementation stages.

The prior focused baseline passed 118 tests across distribution, Skellam, and
terrace-DWBA suites, with persistence tests rerun using a writable temporary
directory. This is baseline evidence, not evidence for the proposed model.

## State definition and interface treatment

### Material ownership and signed corrections

EpitaxyInterface currently adds upper material, adds lower material at its
interface geometry, and subtracts the actual underlying reference material.
Some contributions carry negative algebraic weights. They are scattering and
density bookkeeping, not negative atom populations available for dissolution.

Add an immutable inventory view with stable owner IDs, interface side, source
structural interval, element occupancies, absolute transforms, lateral area,
and reference-replacement provenance. Physical populations must be distinct
from algebraic reference corrections. Keep the atomic and optical provenance
consistent. Reconstruct each state from remaining physical populations and the
reference corrections needed to represent those populations exactly once.

A Film can contribute to interfaces on both sides. Identify all of those
contributions explicitly. Do not infer ownership from names such as `LaNiO3`,
from coincident layer labels, or from equal chemical formulas. A finite STO
Film and the semi-infinite STO substrate have distinct ownership.

### Ordered dissolution front

Introduce a front grid mapping each integer state index to native structural
intervals, remaining material, exposed termination, and absolute front height.
Index zero matches the existing nominal Film boundary. Above that boundary,
the C1-compatible policy continues the top Film's existing layer cycle and
termination templates, retaining the growth states of the recentered family.
Those states contain explicit added inventory; they are not pure dissolution.

Below zero, traverse actual material intervals from top to bottom. In a uniform
LaNiO3 Film, one index step reproduces the present c/2 spacing. In interfaces or
another Film, derive height from the native geometry and explicit layer
transitions. Never multiply every state index by the top Film's c/2.

The first implementation supports compatible ordered layer cycles and
separated interface regions. Unequal or ambiguous interface-cycle mappings,
overlapping interfaces, or a disconnected removable set must fail with a
specific validation error. Expanding these cases requires a separate physical
definition, not automatic interpolation of layer labels.

### Proposed interface closure

Use a frozen, averaged interface for the first implementation: front progression
removes designated upper/lower owner populations in their structural intervals,
while protected populations retain their original occupancies and transforms.
This is a material-selective model, not a claim that the mixed interface is
a collection of uniquely known microscopic terraces. It may expose a porous or
compositionally depleted region rather than pristine STO.

The present averaged interface does not specify a joint distribution of buried
interface height and dissolution-front height. Do not infer lateral correlations
from its one-point occupancy. A model with correlated microscopic interface
terraces would require joint states and is deferred. Confirm the adequacy of
the frozen averaged closure before interpreting interface-sensitive fitted data.

For the first version, remaining atoms do not mechanically relax or restack
after dissolution. Default exposure leaves surviving atoms unchanged; an
optional termination bank can supply an explicit exposed reconstruction.
Existing C1 termination motifs must be recovered where sufficient Film remains.
If a motif extends into a removed or protected region, select a defined thin-film
or interface termination or reject that state. Never silently trim a motif,
subtract absent atoms, or manufacture a complete termination from insufficient
inventory. Multi-plane motif behavior near exhaustion is a required design
fixture before the state builder is finalized.

### Physical boundary and realized statistics

Let J be the existing recentered Skellam parent index and j_min the exhausted
state. The proposed physical mapping is K = max(J, j_min), so

    P(K = j_min) = P(J <= j_min)
    P(K = j) = P(J = j), j > j_min.

Compute the exhausted-state probability directly from a stable CDF. It is a
retained physical state, not omitted numerical mass. This is a static capped
distribution; it is not an absorbing kinetic process, a reflected random walk,
or a conditioned Skellam distribution.

For each K, obtain physical height z_K from the front grid. Report E[z_K],
Var(z_K), RMS height, exhausted probability, exposed-material probabilities,
and per-species remaining/removed/added inventories. Parent mean and raw variance
remain distinct parameter names. Fractional recentering and the boundary change
realized moments, and nonuniform structural spacing changes their conversion.
If a single geometric surface height is ambiguous in a porous interface state,
label z_K as the front coordinate and also report exposed material information;
do not advertise it as a uniquely measured atomic surface height.

## Native architecture and API

Proposed public names below are new API, not existing methods.

    Film and interface material ownership
                    |
                    v
        DissolutionPath and front grid <--- parent height distribution
                    |                                  |
                    +-----------> bounded state weights + tail certificate
                                          |
                                          v
                         complete native crystal states
                              /                    \
                     native F(h,k,l)          native DWBA fields
                              \                    /
                          ensemble endpoint accumulation
                                          |
                          native optimizer and resolution

Add `CTRdissolution.py` for the path, physical state assembly, and inventory
contract. Keep distribution mathematics in `CTRdistributions.py`. Register an
opt-in `FiniteFilmDissolution` ensemble in `CTRincoherent.py` with a stable
configuration type such as `finite_film_dissolution`.

| Proposed object | Responsibility |
| --- | --- |
| `DissolutionPath` | Resolve selected Film IDs, interface-side ownership, protected boundary, front grid, exposure policies, and allowed upper extension. |
| `DissolutionState` | Immutable state identity, absolute geometry, physical inventory, exposed materials, and native crystal construction. |
| `DissolutionStatistics` | Parent parameters, realized physical moments, species inventories, exhausted probability, and support diagnostics. |
| `SupportCertificate` | Raw retained probability, numerical omitted tails, moment bounds, normalization effects, and reason for any resource/admissibility failure. |
| `FiniteFilmDissolution` | Fit parameter plumbing, state weights, native kinematic/DWBA evaluation, persistence, and evaluation-local caching. |

The configuration must explicitly store selected Film IDs, ownership bindings,
front/upper-extension policy, `boundary_policy="accumulate_at_exhaustion"`,
interface closure, termination policies, profile parameters, support settings,
incoherent fraction, DWBA optical-reference policy, and schema version. Reject
duplicate/missing IDs, ambiguous ownership, and unsupported policy versions.
Do not silently infer substrate-removal permission on load.

### Kinematic integration

For complete native state amplitudes F_j, calculate

    F_coh = sum_j p_j F_j
    I_inc = sum_j p_j abs(F_j)**2
    I_mix = (1-kappa)*abs(F_coh)**2 + kappa*I_inc.

Keep substrate/Film/interface interference inside every F_j. A mixed model
exposes intensity rather than a unique complex amplitude; diagnostic coherent
amplitude is explicitly labeled as such.

`CoherentStateEnsembleModel._evaluate_F2` currently uses the wrapped crystal as
its coherent endpoint rather than the reconstructed state sum. The new model
must supply a consistent new endpoint through an explicit override or narrowly
scoped endpoint hook. Do not inherit that behavior blindly or change the
legacy wrappers' coherent endpoints. Recover existing C1 results within the
certified truncation/normalization error when interface-reaching mass is tiny.

### DWBA integration

Reuse native state preparation, optical profiles, polarization handling,
reference subtraction, and kernels. Generalize state provision separately
from the current single-surface `_flat_crystal` helper and Film-only prefix
optimization. The initial correct implementation can construct full native
state crystals; optimize only after equivalence is established.

For `state_specific`, solve each state's native optical fields, align all
amplitudes to the same ambient phase gauge, and mix full reflection amplitudes.
For `shared`, explicitly build the averaged physical optical density and apply
the existing shared-field approximation consistently to each state. The shared
coherent endpoint must agree with native scattering of the averaged atomic
structure on those same fields, within numerical truncation error. Do not
require shared and state-specific optical policies to agree in general.

Do not retain pristine or unchanged C1 fields after an interface inventory
changes. Both atomic records and optical density must use the same state.
Integrate the new model with native DWBA optimizer routing and the existing
intensity-to-observable/resolution sequence; preserve the original C1
kinematic fit convention when replaying that fit.

### Caching and persistence

Start with evaluation-local or explicit batch caches. Cache identity must
include ownership/path topology, state index, parent parameters, termination
policy, atom occupancies/DWF, geometry/strain/registry, crystallinity, optical
policy, energy, domains, and incidence/exit geometry as applicable. Exceptions
must leave source crystals unchanged and invalidate incomplete state caches.

Round-trip new settings through the native incoherent configuration and fit
serialization paths. A crystal file alone is not a full ensemble definition;
require its model configuration companion. Preserve legacy files and defaults,
deepcopy/pickle behavior, optimizer parameter registration, and worker loading.

## Numerical support as a separate workstream

First improve interval selection using the existing native probability and
first-moment bounds and their safety margins. Search for a smaller certified
interval instead of only expanding the initial 8-sigma radius. Separately add
stable one-sided Poisson tails and sharper two-sided Skellam evaluation, with
certified remainders and cancellation checks. Do not relax existing tolerances.

For the bounded physical model, all lower tail mass belongs to the exhausted
state. Only unrepresented upper states and unresolved interface/material
support contribute to numerical omission. Audit probability, physical-height
first/second moments, species inventory, and scattering separately. A first
moment bound does not certify variance or incoherent intensity; a small tail
probability alone does not certify DWBA accuracy.

Any normalization of retained numerical state weights must be explicit,
reported once, and bounded in its effects on moments and observables. It must
never absorb the physical exhausted mass into other states. Resource limits
must fail diagnostically rather than silently clip states. Distinguish an
optional state-count resource ceiling from legacy `exact_layer_count` semantics.

The interface itself uses finite numerical support. Report its residual material
uncertainty independently so the depth of complete dissolution is not defined
by an unaudited 8-sigma window or arbitrary per-bin occupancy cutoff. The
exhausted endpoint removes all modeled finite-owner material; changes under
tighter interface support must converge within declared inventory budgets.

## Implementation sequence and acceptance gates

| Stage | Work and principal files | Completion gate |
| --- | --- | --- |
| 0 | Freeze source/input hashes; add portable synthetic fixtures and an optional read-only C1 replay harness. | Baseline loads reproduce crystal, metadata, rotations, support, and existing native predictions. No test requires the F drive. |
| 1 | Add material provenance/inventory views to `CTRfilm.py` and minimal adapters to `CTRdwba.py`/`CTRoptics.py` as needed. | Reassembling unchanged populations reproduces atom inventory, absolute geometry, complex F, optical profiles, and DWBA. Signed reference corrections never appear as removable atoms. |
| 2 | Implement `DissolutionPath` and deterministic native states in `CTRdissolution.py`. Resolve thin-film terminations and the frozen interface closure using explicit fixtures. | Intact, partially removed, interface-entering, exhausted, and explicit-substrate-Film states pass species and geometry checks. Protected material is unchanged. |
| 3 | Add bounded profile mapping, realized statistics, and support certificates in `CTRdistributions.py`/`CTRdissolution.py`. | Physical exhausted probability is retained exactly to floating-point accuracy; moments and numerical tails agree with independent sums. |
| 4 | Register `FiniteFilmDissolution`; implement native kinematic endpoints, configuration, and fitting parameters in `CTRincoherent.py`, with narrow `CTRopt.py` changes if needed. | Complete-state oracle and endpoint/mixing tests pass; invalid proposals recover cleanly. |
| 5 | Generalize native DWBA state provision and support both optical-reference policies. | Native independent-state oracles, phase, density, polarization, specular/off-specular, and weak-optics tests pass. |
| 6 | Complete serialization, deepcopy/worker tests, cache mutation tests, and separate support-selection improvements. | Round trips and live mutations agree with fresh evaluations; legacy behavior remains unchanged. |
| 7 | Add examples, user docs, changelog, installed-build check, and read-only C1 validation. | Error budgets below pass; no physical-model discrepancy is hidden by refitting scales or parameters. |
| 8 | Prepare a separately numbered analysis fit, only after model validation and an explicit refit request. | New fit records physical policies, realized moments, exposed-substrate probability, inventory, and convergence evidence. |

Numerical interval tightening can be developed independently of stages 1-5;
do not couple its release to a change in physical boundary semantics. No stage
requires C++ changes initially. If profiling later justifies native kernel
changes, rebuild and compare Python/native results before trusting benchmarks.

## Regression coverage and quantitative validation

Add focused modules for dissolution inventory/states and dissolution ensembles,
reusing independent fixture builders rather than reproducing implementation
branches in the test oracle. Extend existing distribution, Skellam, multilayer,
incoherent, optical-inventory, and DWBA suites where contracts are shared.

| Test family | Required cases |
| --- | --- |
| Ownership | Upper and lower interface contributions; a Film between two interfaces; identical-composition owners; explicit STO Film over protected STO; no double removal of reference corrections. |
| Geometry | Zero removal, each layer transition, cyclically renamed layers, absolute registry, nonuniform spacing, interface strain/offset, removal into an interface, and a physically equivalent finite-substrate split. |
| Chemistry | Alternating LaO/NiO2 inventories, partial occupancies, oxygen linkage, multi-plane surface motifs, no ghost termination at exhaustion, and explicit added inventory for positive states. |
| Boundary | All probability at exhaustion; probability spanning exhaustion; no spurious conditioning; exactly zero finite-owner inventory in the exhausted state; protected owner unchanged in every state. |
| Distribution | Both one-sided limits, two-sided rates, zero variance, fractional shifts, asymmetry endpoints, broad tails, moment certification, and continuity at support/state transitions. |
| Kinematic | Complex native F against independently constructed states, full-state interference, kappa 0/0.205425229351/1, four C1 rotations, and lateral-area conversion. |
| DWBA | Both optical policies, shared-profile consistency, common phase origin, polarization channels, critical-angle regions, deep minima, specular/off-specular points, and weak-optics agreement with kinematic results. |
| State lifecycle | Live changes to each dependency, batch versus fresh evaluation, rejected-candidate recovery, source immutability, deepcopy/pickle, worker execution, and configuration/text-fit round trips. |
| Failure behavior | Missing owner, removable bulk, protected gap, incompatible cycles, overlapping interfaces, unspecified exhausted termination, unsupported overlayer/water treatment, and resource exhaustion. |

Proposed numerical acceptance targets, to be calibrated against independent
native reference calculations in stage 0:

- Deterministic inventory and geometry equivalence: relative error at most
  1e-12, with explicit absolute floors for zero populations and coordinates.
- Deterministic complex-amplitude and endpoint equivalence: relative error at
  most 1e-9 plus a scale-aware absolute bound. Use a fixed reference amplitude
  scale near extinctions, not relative error divided by a vanishing amplitude.
- Probability normalization error at most 1e-12 for test fixtures; separately
  report all certified omitted mass and normalization-induced moment changes.
- Preserve the C1 helper's existing numerical audit: omitted mass <1e-7,
  omitted absolute first moment <1e-5 planes, retained parent-mean error <1e-6
  planes wherever that unchanged parent-distribution audit applies. Do not
  apply a parent-mean equality test to a physically capped distribution.
- Require a separately declared physical-height second-moment error budget.
  Tightening that budget and upper/interface support must stabilize realized
  RMS and species inventories; first-moment success cannot substitute for it.
- For C1 numerical convergence, proposed maximum change in predicted fitted
  observable is 1e-3 of each supplied standard error, checked at fixed physical
  parameters and fixed nuisance scales after original resolution processing.
  Record the resulting objective change. Compare sharper-support calculations
  to still tighter references, not only to the current implementation.

Physical differences from new interface or exhaustion states must be reported
separately from numerical errors. These targets are implementation acceptance
criteria, not changes to the fit's error bars or fit weights. If a reference
cannot meet them, document the conditioning and establish a justified budget
before accepting results; do not loosen tolerances to get a test to pass.

Run narrow new tests first, then affected existing tests and `ruff check orgui`.
Use a fresh workspace-local pytest base directory where necessary. Verify both
source and installed-package paths for the C1 replay. The canonical JSON is
required to restore energy, reference cell, attenuation, incoherence, and all
four rotations; the xtal file alone is insufficient. Preserve the original
resolution kernel and observable conversion for each forward model.

## Documentation and release

Add a synthetic example demonstrating complete upper-Film dissolution and a
second example with a finite substrate Film. Include inventories and exhausted
probability, not just roughness curves. Document parent versus realized moments,
the frozen interface closure, termination rules, and both DWBA optical policies
in the existing CTR and DWBA topic pages. Add Sphinx-friendly docstrings to every
new public function, including units and failure conditions.

Add user-facing release notes through `CHANGELOG.md` and regenerate
`doc/source/release_notes.rst`; do not edit the generated file directly. Use
Conventional Commits with `phys` scope for changes affecting generated data or
measurement results. No changelog entry is needed for this plan alone.

## Deliberate limits and future work

The first version is a static ensemble, not time-dependent dissolution kinetics
or an off-rod diffuse-scattering model. It does not infer microscopic lateral
interface correlations, dissolve the semi-infinite bulk, restore depleted STO,
or relax newly exposed atoms. Attached water/overlayers require an explicit
state-placement policy and tests before support is claimed; otherwise reject
those configurations in the new model without changing legacy behavior.

The critical scientific choices to revisit are interface correlations, exposed
termination chemistry, and whether substrate relaxation is needed. The main
engineering tradeoff is complete-state native evaluation versus faster additive
prefix evaluation. Start with complete states as the correctness reference;
retain the same reference when optimizing and when extending to additional
materials or more complex interface geometry.

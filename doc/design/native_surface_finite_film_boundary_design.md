# ADR: Finite-film boundaries in the native surface model

**Status:** Implemented in the dissolution worktree; validation recorded below.
**Date:** 2026-11-07.
**Decider:** Project maintainer.

## Context and decision

### Complete substrate preservation (2026-11-08)

The protected reference is the complete original substrate/interface profile,
not pristine sharp bulk alone. The original frozen closure already retains
that profile at every dissolution depth. Signed reference terms describe its
initial departure from the sharp substrate; they are not substrate dissolution.
The maintainer explicitly confirmed preservation of the original profile on
both sides of z=0. An attempted bare-bulk-only correction removed the real
upper substrate extension and was reverted before commit, together with its
experimental `substrate_policy` API and erroneous refit warning. Density
differences from a sharp reference do not establish invalidity of earlier fits.

Regressions now isolate complete STO cells, including oxygen and signed
references, and compare density and complex amplitude before, during and after
removal. Exhaustion must recover the same original rough-substrate density and
native DWBA scattering for both optical-reference and coherence policies.
Removable substrate material still requires an explicit finite substrate Film.

Protection is independence from parameters that do not belong to the substrate
or its interface, not a permanently cached substrate. Across evaluations,
surface width, asymmetry, mean, numerical support and unrelated Film/overlayer
changes must retain the same protected coordinates, occupancies and intrinsic
optical profile. Changes to the substrate/interface itself still regenerate
those records. DWBA illumination can change with the overlying material even
when the substrate material is unchanged.

The averaged crystal now inserts protected records once with their original
weights and applies probabilities only to removable populations. This removes
distribution-dependent duplication and summation of signed substrate records.
No interface coordinates, probability tolerances or scattering kernels change.
Regressions check exact protected atom records, density, amplitude and optical
profiles under unrelated parameter changes, as well as sensitivity to actual
interface changes. An explicit finite substrate partition is checked against
the original nominal density and, with bulk damping disabled, its complex
amplitude after accounting for the common origin translation. The latter
qualification matters: native numerical bulk damping is not applied to finite
Film terms, so repartitioning at nonzero damping is not an exact kinematic
amplitude identity. Density comparisons retain the saved damping setting,
which does not enter their evaluation.

Validation of this follow-up: 258 focused tests passed, covering the unchanged
interface anchors and registration, native dissolution, boundaries, overlayers,
substrate independence and finite-substrate partitioning. Local Ruff and diff
whitespace checks passed. No geometry-test expectations were changed.

Keep `Film`, `EpitaxyInterface`, `SkellamSurface` and `PoissonSurface` as the
structural model. Make complete finite-material dissolution an explicit boundary
policy of a native surface within that stack. Use `SurfaceHeightDomains` for
coherent/incoherent intensity mixing. New models must not need a
`FiniteFilmDissolution` wrapper, a second set of surface parameters, or a second
crystal definition.

`FiniteFilmDissolution` has never been pushed or published. Remove its API,
registration and separate execution path outright when this refactor lands.
There is no deprecation period, compatibility adapter or requirement to load its
experimental configuration/pickle format. Compatibility requirements below
apply to the established native Film/surface/domain APIs.

This changes ownership and evaluation plumbing, not the physical boundary law.
The proposed finite boundary still accumulates exhausted probability and removes
finite material, including its mixed-interface share. New surfaces select all
finite Films by default; an explicit boundary can restrict that selection.
The semi-infinite bulk remains protected. Dissolvable substrate material must
be represented by an explicit finite `Film` of that material.

The current implementation already uses native cells for scattering; the
duplication is in model ownership and orchestration:

- `FiniteFilmDissolution.__init__` copies the surface distribution into a local
  `[t, k_up, k_down, offset, incoherent_fraction]` basis.
- `DissolutionPath.__init__` rejects fitted surface parameters and substitutes a
  deterministic profile in a copied crystal before laying out material.
- Its support helper accepts only `RecenteredSkellamProfile`, although the native
  geometry can come from a Poisson surface.
- `FiniteFilmDissolution` and `DissolutionDWBA` supply separate evaluation,
  statistics, persistence and batch entry points alongside the native surface
  and `SurfaceHeightDomains` paths.

Retain the material ownership, signed-reference accounting, native amplitude
reuse and regression oracles developed there. Move those capabilities behind a
common native stack/state interface rather than discarding them.

## Constraints

1. Native surface parameters are the sole source of height probabilities in a
   new model. Fitting a surface changes its coherent, incoherent, density and
   DWBA predictions together.
2. Old saved models without boundary metadata, and surfaces constructed with
   explicit `boundary=LegacySurfaceBoundary()`, retain their current behavior, support rules and
   numerical tolerances. New constructors default to all finite Films. In particular, retain the
   existing Skellam finite-support rejection in legacy mode. Do not introduce a
   new Poisson rejection or new tail selection merely by loading an old model.
3. Preserve the chosen exhaustion, interface and exposure laws. No silent
   clipping, conditioning on survival, tolerance relaxation, interface movement,
   or increase in nominal Film width.
4. Preserve native reference-area scaling, signed interface subtraction,
   crystallinity, outer-domain transforms, absolute coordinates and phase.
5. Keep native kinematic/DWBA kernels. A shared evaluation plan is temporary
   derived data, with no fit parameters or independently editable geometry.
6. No refit, private dataset dependency or migration of stored files as part of
   this design. Repository tests and benchmark acceptance use synthetic models.

## Ownership and public API

| Owner | Authoritative data |
| --- | --- |
| `Film` | Material, nominal structural width, lattice/cycle and finite inventory |
| `EpitaxyInterface` | Native mixing geometry and signed replacement records |
| `SkellamSurface` | `t`, `k_up`, `k_down`, `offset`, termination bank and surface order |
| `PoissonSurface` | `W`, `alpha`, `offset`, termination bank and surface order |
| `RoughSurface` settings | Boundary policy and numerical height-support settings |
| `SXRDCrystal` | Component topology, reference cell/area, bulk and bound material evaluation |
| `SurfaceHeightDomains` | Incoherent fraction and DWBA optical-reference selection |

`PoissonHeightDomains` remains the existing compatibility name with its Poisson
target check. Neither domain wrapper gains surface-height fit parameters.

### Proposed constructors

The following names and keywords are proposed API, not executable examples of
the current implementation. Policy classes are immutable settings, not models.

```python
from orgui.datautils.xrayutils.CTRfilm import (
    Film, SkellamSurface, PoissonSurface, FiniteFilmBoundary,
)
from orgui.datautils.xrayutils.CTRdistributions import (
    RecenteredSkellamProfile, PoissonProfile, HeightSupport,
)
from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRincoherent import SurfaceHeightDomains

film = Film(material_cell, width_layers=40, name="film")
surface = SkellamSurface(
    termination_cell,
    name="surface",
    profile=RecenteredSkellamProfile(t=1., k_up=0., k_down=5., offset=-3.),
    boundary=FiniteFilmBoundary(
        films=("film",),
        boundary_policy="accumulate_at_exhaustion",
        interface_closure="frozen_material_selective",
        exposure_policy="native_or_unreconstructed",
        upper_policy="continue_top_film",
    ),
    height_support=HeightSupport(
        tail_probability=1e-10,
        moment_tolerance=1e-10,
        second_moment_tolerance=1e-8,
        max_states=4096,
    ),
)
crystal = SXRDCrystal(substrate, film, surface, stacking=[1, 2])
crystal.setEnergy(10000.)  # eV

# Native coherent evaluation includes the same finite-material boundary.
coherent_intensity = crystal.F2(h, k, ell)

# Optional intensity mixing; surface parameters remain on crystal["surface"].
domains = SurfaceHeightDomains(
    crystal, surface="surface", incoherent_fraction=0.2,
    optical_reference="shared",
)
mixed_intensity = domains.F2(h, k, ell)
```

For Poisson growth/etching, replace only the surface construction:

```python
surface = PoissonSurface(
    termination_cell,
    name="surface",
    profile=PoissonProfile(mean_change=-5., alpha=0.7, offset=-2.),
    boundary=FiniteFilmBoundary(films=("film",)),
)
```

The omitted finite-boundary policy fields have exactly the defaults shown in
the Skellam example. Persist policy values explicitly. Omitting `boundary`
selects `FiniteFilmBoundary()`: all finite Films are removable. Explicit
`boundary=LegacySurfaceBoundary()` selects legacy behavior, as do old saved files without boundary
metadata. The numerical defaults apply to the finite-boundary path only.

`films=None` (the default) infers all finite Films from the live stack. It stays
automatic across text and HDF5 round trips rather than storing a frozen list.
An explicit nonempty tuple identifies existing components by stable name and
must form a connected upper sequence. Empty tuples are invalid. Material ownership in
interfaces is derived from native adjacency, not a second editable material
mapping. Persisted owner information is an audit/validation record. An explicit
finite substrate Film can be included by name; `bulk` cannot.

Do not automatically create a substrate Film or redistribute thickness when a
boundary is enabled. A separate repartition operation would need to demonstrate
unchanged species inventory, coordinates, phase and scattering before any
dissolution is applied.

### Fitting and numerical settings

Surface fitting keeps the existing native parameter names and machinery. Film
width, interface parameters and atom parameters also remain on their native
components. The domain wrapper retains only its incoherent-fraction basis.
Changing `boundary` changes a physical policy; it is not a continuous fit
parameter.

Extend the native surface's numerical-settings storage rather than adding a
second authoritative profile to a wrapper. `HeightSupport` is a validated view
of those settings. The profile returned by `_make_profile()` is an immutable
evaluation snapshot derived from native basis and settings.

- With `height_support=None`, retain tolerances supplied by the native profile;
  supply named defaults only for settings that profile did not previously have.
- If both an explicit profile and explicit support settings specify the same
  tolerance differently, reject the ambiguity. Do not silently prefer or relax
  one. A deliberate settings update is validated and invalidates evaluation
  caches.
- Keep the existing Skellam `profile_settings` reader as a compatibility reader;
  serialize one resolved numerical-settings record in the new format.
- `max_states` remains a resource-error threshold. For finite boundaries,
  `exact_layer_count` cannot activate the legacy Poisson approximation that
  aggregates occupied heights. Keep every retained nonzero bounded state, or
  raise the resource error. Legacy surface-domain behavior remains unchanged.

### Inspection

Add stack-bound inspection methods rather than putting another crystal on the
surface:

```python
statistics = crystal.surface_height_statistics(surface="surface")
state = crystal.surface_height_state(surface="surface", index=-12)
```

The state index is dimensionless, relative to the nominal top Film boundary.
It counts the existing ordered structural intervals; it is not an Angstrom
depth or an assumed constant plane spacing. State inspection may explicitly
materialize a complete native crystal. Routine evaluation must not do so for
every height.

Statistics retain parent and realized index/physical-height moments, exhausted
probability, exposed owners, species inventories per reference lateral cell,
normalization factors and tail/interface audit bounds. Physical heights are in
Angstrom. Wrapper statistics delegate to this same native API.

## One structural-state provider

Introduce a private stack-bound `_SurfaceStatePlan` shared by coherent evaluation,
domain averaging and DWBA. Its finite-material implementation absorbs the
useful parts of `DissolutionPath`; its legacy implementation adapts the existing
flat-height corrections. This is an evaluation abstraction, not a new saved
structure or a generic framework for unrelated disorder models.

It supplies: component provenance, protected and removable records, ordered
intervals, state probabilities, deterministic termination/reference corrections,
averaged material records, optical-grid requirements and inventory diagnostics.
The complete-state materializer remains an inspection/reference operation.

### Binding and stacking

The native surface must be bound to its Film and full lower stack before its
finite support can be evaluated. Refactor the current stacking sequence into
geometry binding and statistical material assembly:

1. Validate component ordering, names, reference transforms, Film/interface
   bindings and material ownership using nominal native geometry.
2. Build structural intervals and the exhaustion index from that geometry.
3. Construct the surface's actual live profile and bounded state distribution.
4. Assemble the effective material contributions and optical inventory.

Do not implement step 1 by replacing the live profile with a flat profile,
zeroing its fit parameters, or catching and suppressing the old finite-support
exception. Internal geometry-only binding must have an explicit contract and
must not recursively invoke statistical scattering.

Retain one active finite-boundary surface initially. Overlayers above it are
common native components at the realized mean front, outside the removable
inventory; they retain their own weights and domains. Reject unsupported
attached water, overlapping interfaces, incompatible cycles,
non-unit component weights, differing outer domains and depth-dependent order
where the current finite-material engine does not define them. The design is
general across native Poisson/Skellam families and selected finite Films; it
does not invent joint distributions for multiple rough boundaries.

### Distribution-independent exhaustion

Generalize `capped_skellam_states` behind a discrete-height distribution contract.
Required operations are integer PMF, stable CDF/survival probability, expected
height, variance, and certified upper-tail mass/absolute first/second moments.
Implement adapters for the two native profile types, with independent tests.
Do not depend on the Skellam-specific `mu_up`, `mu_down` or `_shift` fields in
the common state selector.

For parent integer height `H` and geometric exhaustion index `m`, the physical
state is `J = max(H, m)`. Its mass at `m` is `P(H <= m)`; masses above `m` retain
their parent values before the separately reported numerical-tail treatment.
This is accumulation at exhaustion, not a distribution conditioned on survival
and not a bounded kinetic process.

A fixed parent mean does not fix the realized bounded mean. Zero mean index
change also does not imply conserved mass when different intervals contain
different species or populations. Compute material balance from native state
inventories, not from a height variance or plane count alone.

Poisson must retain its signed step-plus-Poisson law, including fractional
deterministic shifts and `alpha=0`. Its expected height is `W + offset`;
`PoissonProfile.mean` is currently the legacy alias for `W` and must not be used
as the general actual mean. Use `expected_height_change`. Its variance is
`alpha * abs(W) + f * (1-f)`, where `f` is the fractional part of the magnitude
of the deterministic structural cursor. Skellam retains its existing
recentering and interpolation variance. Do not serialize either family as the
other or equate fitted means/variances across them.

Preserve the existing explicit non-exhausted upper-tail normalization rule and
certificate when replacing the current finite-dissolution implementation. Do not
modify its error budgets in this refactor. Report the retained mass and normalization factor;
never rescale exhausted probability. Separate any future numerical-support
improvement from this structural migration.

## Native observables and coherence

### Kinematic scattering and density

`SXRDCrystal.F` evaluates the bounded mean amplitude from the effective native
material records. `F2` is its squared magnitude. `SurfaceHeightDomains` uses the
same state plan:

`I = (1-kappa) * |sum(p_j * A_j)|^2 + kappa * sum(p_j * |A_j|^2)`.

Retain the existing optimized protected-material evaluation and material-prefix
sums. Preserve custom-cell/native complete-state fallbacks. No optical profile
or field solver is permitted on this kinematic path.

`evaluate_kinematic` and `F_surf` must reflect the same effective material.
Preserve component provenance in their decomposition: selected Film and
interface contributions change where their material is removed. Do not hide
all cross-interface changes in an unexplained surface-only amplitude.

`zDensity_G`, physical optical inventory and `optical_profile` must consume the
same effective material definition. Low crystallinity changes scattering order,
not physical material inventory. Source components remain the native editable
definitions; generated evaluation records are temporary views, never serialized
as extra Films or substituted into `uc_surface_list`.

### DWBA

The native crystal's coherent DWBA evaluates the averaged bounded material on
its common optical reference. Shared-reference domain averaging uses the same
all-state union grid, including exhausted-state boundaries and empty media.
Thus its `kappa=0` result must agree with native crystal DWBA. Do not let native
coherent evaluation silently use a different slab partition.

For `optical_reference="state_specific"`, preserve distinct state fields and
their common phase convention. Its coherent endpoint remains the weighted sum
of those state amplitudes and need not equal averaged-reference native DWBA.
Do not collapse these two physical optical policies during unification.

`SurfaceHeightDWBA` consumes the common state provider; finite boundaries no
longer require a separate public `DissolutionDWBA` engine. Private optimized
strategies may differ by record topology and optical policy. Share the native
packing, field preparation, reference subtraction, amplitude prefactor and
batch lifecycle. Do not replace the native kernel or infer a new reference
profile merely from the parent height histogram.

Keep the current snapshot/packing/shared-field reuse and exception cleanup.
Further shared-material DWBA acceleration is a later optimization of this
provider, not another public model. Distinct state-specific fields are not
eligible for amplitude reuse solely because atomic coordinates match.

## Persistence and replacement of the unpublished implementation

### New configuration

Store a versioned `height_boundary` record and resolved `height_support` settings
on the native surface. The boundary record contains type, selected Film names
and the explicit physical policy values; it contains no fitted surface basis.
The surface's existing parameter record remains authoritative. Domain-wrapper
configuration stores its surface name, incoherent fraction and optical policy.

Extend native text parsing/writing and HDF5 fit settings together. Loading a
crystal without its incoherent wrapper must still preserve the finite-boundary
coherent structure. Reject unknown boundary versions, missing referenced Films,
ambiguous names and mismatched ownership audits.

An optional comment alone is insufficient forward-version protection: old
readers may ignore it and evaluate a different physical model. Require an
explicit version/capability guard that makes unsupported boundary-bearing text
files fail in supported older readers, or introduce a guarded format version.
If that cannot be guaranteed for historical readers, disable legacy-text export
of boundary-bearing models and use the versioned structured format until it can.
Never offer a silent legacy export that drops the boundary policy.

### Remove `FiniteFilmDissolution`

Delete the class, its public export, the `finite_film_dissolution` registry entry
and wrapper-specific configuration handling. Merge the material accounting into
the shared native state provider, then remove `_CTRdissolution_dwba.py` and the
separate DWBA class. Move or rename retained internals so there is no public
parallel dissolution model. Do not leave an alias or constructor that secretly
builds another wrapper.

Rewrite examples and tests to construct native Film/surface stacks and standard
domain wrappers. Update benchmarks to use those APIs. Preserve scientific oracle
fixtures in tests; preserving their expected physics does not require retaining
the old production class. Remove outdated unpublished API documentation and
rewrite its Unreleased changelog entry around the native implementation.

Any local experimental parameter sets needed for a before/after comparison can
be converted explicitly in test fixtures or a one-off script. Do not introduce a
public migration API or automatic loader for them. The parameter correspondence
for those comparisons is:

| Unpublished wrapper source | Native destination |
| --- | --- |
| Wrapper `t`, `k_up`, `k_down`, `offset` | Target `SkellamSurface` basis |
| Wrapper `basis_0`, fit bounds, priors and errors | Corresponding native surface parameter records |
| Wrapper incoherent fraction | `SurfaceHeightDomains.incoherent_fraction` |
| Selected Films and physical policies | Target surface boundary settings |
| Tail/moment/resource settings | Native surface height-support settings |
| Optical-reference selection | Standard domain wrapper |
| Interface owner map | Validated stack-derived ownership audit |

Use the wrapper's actual values for comparison fixtures; the wrapped surface may
contain older values. If local fitted models are converted, preserve their
current/baseline values, bounds, priors, errors and parameter/covariance mappings
explicitly. Reject unresolved constraints rather than silently dropping them or
refitting. This is a validation concern, not a supported legacy API contract.

Do not automatically convert an experimental wrapper on a `PoissonSurface`
into Poisson parameters: the wrapper may have an independent two-sided Skellam
law. Such a mismatch requires an explicit surface-family replacement preserving
the termination geometry, or is rejected. Conversely, native Poisson models
remain native Poisson models when opting into a boundary.

There is no requirement to preserve the old five-entry wrapper basis, constructor,
import path, model-type key, pickles or mutable-array access. New standard domain
wrappers have their usual fraction basis; scripts fit height parameters on the
native surface. Unknown experimental model-type keys should fail clearly instead
of silently selecting a different physical model.

Continue to test established native text/HDF5/configuration/deepcopy/pickle paths,
including new round trips with the explicit boundary setting. Their compatibility
is independent of removing the unpublished wrapper.

## Options and trade-offs

| Option | Complexity | Maintenance cost | Evaluation scaling | Repository fit |
| --- | --- | --- | --- | --- |
| Keep the separate wrapper, share a few helpers | Low initially | Two parameter/persistence paths remain | Current optimizations retained | Existing code, but does not meet the requested ownership model |
| Only relax surface support and extend the top Film | Low superficially | Physically inconsistent at interfaces | Cheap but wrong inventory | Rejected: moves/adds material or leaves interface material behind |
| Native boundary settings plus one stack/state provider | Moderate; binding and persistence work required | One authoritative native structure | Retains prefix and native DWBA reuse | Recommended; builds on existing Film/surface/domain APIs |
| Replace all structures with a universal ensemble graph | High | Broad framework and migration burden | Unproven | Unnecessary scope for this change |

The recommended option makes coherent-only models self-contained and keeps fit
ownership understandable. It requires more than renaming the current wrapper:
stack binding, all observables, decomposition and save/load must be integrated.
The private provider is necessary derived bookkeeping; deleting that bookkeeping
would lose material ownership and signed-interface correctness.

## Implementation sequence and acceptance

1. **Freeze native oracles.** Keep current finite-dissolution predictions and
   independent plain-Film/explicit-interface fixtures. Cover both native profile
   families without private datasets. No expected-number updates to hide a
   migration discrepancy.
2. **Generalize distribution support.** Add the two native adapters and bounded
   selector. Verify exact discrete limits, fractional offsets, upper-tail bounds,
   moment certificates, exhaustion mass and resource failures. Leave old support
   methods unchanged.
3. **Separate geometry binding from statistical assembly.** Introduce the common
   state provider; preserve interval geometry, signed references and exposure
   rules without fake profiles or duplicate fit parameters.
4. **Integrate coherent native observables.** Implement boundary settings and
   `F`, `F2`, decomposition, density, optics and inspection through the same plan.
   Gate release on consistency across these entry points.
5. **Unify domain and DWBA evaluation.** Route standard wrappers through the
   provider; retain both optical policies, all mixture endpoints and native
   batch/field reuse. Keep optimized kinematics independent of DWBA work.
6. **Implement persistence and remove the unpublished API.** Test guarded text,
   HDF5, config/deepcopy/pickle and native optimizer registration. Delete the
   wrapper, registry entry and separate DWBA path; port examples and regression
   fixtures without a compatibility layer.
7. **Benchmark and document.** Update synthetic examples, user docs and release
   notes once implemented. Measure fresh and batched calls, changed profiles,
   channels and ray blocks, memory bounds and exception cleanup. Verify there
   is no remaining production dependency on the removed wrapper.

Required gates include:

- New constructor defaults remove all finite Films; explicit selection protects
  lower Films. Legacy models retain their serialized interpretation.
- Skellam and signed Poisson cover zero variance, growth, etching, fractional
  shifts and complete exhaustion; actual Poisson means are not confused with `W`.
- One authoritative parameter update changes every relevant observable. No stale
  profile after fit updates, geometry changes, serialization or batch exit.
- State inventories, means/variances and parent/realized distinctions match the
  present physical model. Semi-infinite bulk is never removed implicitly.
- Include multiple selected Films, mixed interfaces, protected material,
  unequal lateral areas/spacings, nonidentity common domains, incomplete motifs,
  signed references, constant crystallinity, empty media and tiny retained mass.
- Compare complex amplitudes and coherent/incoherent intensities, not only a
  fitted objective. Keep existing regression tolerances; shared optical-support
  observable differences remain within `1e-9`, and exact DWBA reuse comparisons
  retain `rtol=1e-10, atol=1e-28` on the established synthetic fixtures.
- Coherent kinematic equality with the standard wrapper at zero incoherence;
  shared-DWBA equality on one common grid; intentionally distinct state-specific
  optics. Exercise all polarization pairs and scalar/broadcast inputs.
- Preserve one kinematic bulk evaluation, avoid full crystal creation per height
  on the native fast path, and require no optical solves for kinematics. Repeated
  immutable batches add no redundant snapshots/packing/field solves.
- Compare runtime with the current optimized implementation using at least five
  alternating trials on repository fixtures. Investigate regressions exceeding
  5% beyond run variability. No performance acceptance through reduced support.
- Rebuild and test the native extension if C++ changes; validate an isolated
  package without a stale in-place extension shadowing it.

## Review points

The maintainer accepted native surface ownership and standard domain wrappers,
then requested implementation in the worktree. The subsequent default-policy
decision makes all finite Films removable when `boundary` is omitted. Explicit
`FiniteFilmBoundary(films=...)` limits that range; the bulk remains protected.
Removal of the unpublished wrapper without API preservation was
explicitly agreed. No data were refitted or saved experimental models modified.

Related records: [original physical implementation](finite_film_dissolution_implementation_plan.md),
[performance plan](finite_film_dissolution_performance_plan.md), and
[measured implementation](finite_film_dissolution_performance.md).

## Implementation and validation

### Common-overlayer regression fix

The original dissolution engine rejected any component above its target. Making
that engine the native default exposed a regression: the established native
surface/domain model already allowed a common overlayer at the mean surface
height. The planner now partitions the stack at the selected surface, computes
the realized mean from the bounded distribution, and binds the suffix with
native stacking at that height. Suffix components remain common in all states;
they are not included in removable-film selection or material inventory audits.
The original native components retain their scattering, density, weights,
outer domains and optical references. Kinematic evaluation still computes the
bulk and common suffix once. Separately stacked analytic water retains its
native reference; shared analytic-water ensemble fields remain unsupported.

Validation: 26 overlayer regressions cover both height laws, native Film and
UnitCell caps, zero-crystallinity caps, separate layered/step water, fixed mean
placement through exhaustion, native complete-state DWBA references, mixing,
weights/domains, persistence and one-bulk/no-optics kinematic evaluation. Across
the focused scattering, water, optical and native API runs, 619 unique tests
and 118 subtests passed. Ruff and diff whitespace checks passed.

### Default-policy follow-up

New native surface constructors now default to `FiniteFilmBoundary()`. Its
omitted argument and explicit `boundary=None` are equivalent. Legacy behavior
requires `LegacySurfaceBoundary()`, which is also the resolved policy when
loading old text, parameter dictionaries or pickles without boundary metadata.
New legacy saves record their policy explicitly. Dispatch uses the policy type.

The finite policy's
`films=None` selects every finite Film in the current stack. A named nonempty
upper sequence restricts removal; semi-infinite bulk is always protected.
Text/HDF5 round trips preserve automatic selection. Explicit `boundary=LegacySurfaceBoundary()`
and old persisted models without boundary metadata retain legacy behavior.

The default change does not relax topology validation. In particular, some
independently displaced interface geometries have nonmonotonic removal-front
heights and remain unsupported, as do water and depth-dependent order. Existing
tests of those native legacy behaviors now request `boundary=LegacySurfaceBoundary()` explicitly.

Validation after this follow-up: 34 native-boundary tests passed, including
default finite-substrate removal, explicit limits, broad signed Poisson/Skellam
laws, mixed interfaces, persistence and rejection of unordered fronts. Across
the broad CTR/scattering run and affected-file reruns, 1,400 unique tests and
193 subtests passed, with 3 skipped. Ruff and the synthetic example passed;
Sphinx reported its four existing warnings. The isolated-wheel and performance
results below predate this default-policy follow-up; no native kernels changed.

The implemented API is `FiniteFilmBoundary` on either native surface, with
immutable `HeightSupport` budgets. Surface fitting parameters remain native;
`SurfaceHeightDomains` retains only its existing incoherent-fraction parameter.
`FiniteFilmDissolution`, its registration, `CTRdissolution.py`, and the separate
DWBA class/module have been removed without aliases or config migration.

`_CTRsurface_states._SurfaceStatePlan` derives finite material states from the
editable native stack. It has no independent fitting basis or policy settings.
Legacy surfaces retain their established correction stream. Boundary-enabled
surfaces bind nominal geometry without replacing their profiles or zeroing their
parameters. The discrete-law adapter handles signed Poisson exactly, including
its offset and fractional-step variance, before the common bounded selector.

Native coherent amplitude, component decomposition, density, and DWBA now use
the bounded material. The standard domain facade retains the optimized prefix
amplitudes and both DWBA optical policies. Custom cell evaluators use complete
states where linear averaging cannot be assumed. Native coherent and shared
domain DWBA use identical union grids. Public stratified-profile and wavefield
calls preserve their explicit interfaces; sample-center midpoints cannot recover
those interfaces in general.

Native text stores boundary/support JSON behind the mandatory
`FiniteFilmBoundaryV1` / `RequiresFiniteFilmBoundary` guard. Established surface
readers reject the second token instead of ignoring the physical policy.
HDF5 settings normalize scalar storage types. Bounded numerical budgets have a
single stored authority, `height_support`; ordinary Skellam models retain their
existing `profile_settings`. Ownership audits are derived from the live topology
and exposed in statistics, rather than stored as a second editable ownership
definition. Native load/save validates the selected Film names and policies.

Validation includes independent plain-Film endpoints, protected mixed-interface
populations, explicitly removable substrate Films, unequal reference areas,
domain transforms, signed references, complete-state amplitude oracles, both
native distributions, optimizer updates, guarded text, HDF5, config, deepcopy,
pickle, all four polarization channels, batch reuse, and exception cleanup.
The native kinematic fast paths evaluate bulk once and perform no optical
preparation. The final wheel is checked for duplicate members, removed modules,
exact Python-source identity and agreement with the freshly built native binary.

Executed checks: **1,389 tests and 193 subtests passed** in the full selected
CTR/scattering regression run (three skips); **323 tests and 37 subtests passed**
from the isolated wheel. The focused native-boundary tests also verify realized
index moments separately from parent moments. `ruff check orgui` and
`git diff --check` pass. The synthetic example runs with both protected and
explicitly finite substrate inventories. Sphinx builds successfully with four
pre-existing documentation warnings outside the changed sections.

### Synthetic performance checks

The repository-contained command is:

```powershell
python -m benchmarks.benchmark_surface_boundary --points 512 --repeats 5
```

For a 25-plane synthetic Film and 67 retained states, an initial three-trial run
measured median kinematic times of **41.2 ms** for native coherent evaluation and
**33.2 ms** for mixed domains, versus **64.4 ms** for the complete-state oracle.
Shared DWBA measured **6.48 s** versus **12.10 s** for that oracle. State-specific
DWBA measured **8.56 s** versus **8.23 s**; no speedup is established for that
case. These are synthetic end-to-end timings, not a new fitted-data benchmark.

A separate five-trial alternating comparison loaded the exact optimized
pre-refactor implementation from `3268ae421` and used the same native kernels,
fixtures and tolerances. All three policies produced bitwise-identical
intensities. Median seconds were:

| Path | Before | Native API refactor |
| --- | ---: | ---: |
| Kinematic domains | 0.0837 | 0.0754 |
| Shared DWBA domains | 5.765 | 6.411 |
| State-specific DWBA domains | 7.208 | 7.527 |

Other Python jobs were running on the host. Shared-DWBA trial ranges overlapped
substantially (before 5.60–9.08 s; after 6.03–6.94 s), so the median shift does
not establish an intrinsic slowdown. An instrumented before/after comparison
at 32 ray points found **identical work counts**, including native exponentials:

| First evaluation | Shared | State-specific |
| --- | ---: | ---: |
| Atomic snapshots | 68 | 68 |
| Atom packings | 68 | 67 |
| Native evaluations | 68 | 67 |
| Field solves | 4 | 136 |

Repeated immutable batches added **zero snapshots, packings or field solves**
and retained the same native-cache hits as before. Timing evidence supports
preserved computational work; it does not establish a new state-specific DWBA
speed improvement. No support reduction or tolerance relaxation was used.

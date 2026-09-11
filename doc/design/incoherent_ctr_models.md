# Incoherent CTR models: design and implementation plan

Status: draft, updated 2026-09-11. Revised the same day after an
implementation-readiness review against the code: the height-state index
convention, the split between coherent and incoherent retention policies, and
the treatment of `CTROptimizer.startp` were all corrected before any code was
written. Each correction is stated where it applies rather than collected in a
separate errata section.

This record defines an extensible, opt-in API for incoherent CTR calculations.
The first implementation is a surface whose local height follows the existing
`PoissonProfile`, but whose laterally large height domains do not all scatter
coherently. It also defines the partial-coherence extension between the
existing coherent limit and a fully incoherent domain average.

The design deliberately does **not** implement Sinha diffuse scattering, a
height-height power spectral density, or an off-rod diffuse line shape. It is a
domain-mixture model for kinematical CTR `F2`. Domain boundaries and scattering
that they place outside the measured CTR acceptance are outside its scope.

## Decisions

1. Existing `SXRDCrystal.F` and `PoissonSurface.F_uc` behavior remains fully
   coherent and is the default. Merely constructing a `PoissonSurface` must not
   opt a calculation into incoherent averaging.
2. An incoherent kinematical `IncoherentF2Model` consumes complete coherent
   state amplitudes and returns `F2 = abs(F) ** 2`, in the squared units of the
   existing kinematical structure factor. It has no `F` method because a mixed
   state has no unique complex amplitude. `F2` is not detector counts and is
   not reflectivity.
3. In a kinematical calculation, incoherent averaging acts on `F2` for the
   **complete coherent crystal state**, `abs(F_total_n) ** 2`. It does not act
   on the isolated surface correction and it does not use optical reflectivity.
4. DWBA is not supported in the first implementation. A future DWBA extension
   would mix `abs(r_total_n) ** 2` for independently prepared height states;
   it must not mix `abs(F_h_n) ** 2` or reuse one height-averaged optical field.
5. The public separation between coherent and incoherent models must not force
   the implementation to deep-copy or reevaluate the common bulk/Film stack
   for every height. `SXRDCrystal` will expose an immutable decomposition of a
   kinematical evaluation, and `PoissonSurface` will expose flat-height
   corrections. `PoissonHeightDomains` combines those results through the
   general context API.
6. The initial partial-coherence parameter is an effective
   `incoherent_fraction`, denoted by `kappa` below. It is dimensionless and lies
   in `[0, 1]`. Using the coherent crystal directly is the default;
   constructing an incoherent wrapper is explicit opt-in.
7. A beam coherence length is not stored on `PoissonSurface`. Surface-domain
   correlation is a sample property, while beam mutual coherence, footprint,
   and detector acceptance are measurement properties. A later physical
   coherence kernel may combine them to calculate a point-dependent `kappa`.
8. The new module is `CTRincoherent`, not `CTRintensity`. It owns a general
   `IncoherentModel` wrapper/parameter contract, the kinematical
   `IncoherentF2Model` quantity contract, a registry for stable model type
   names, and the Poisson height-domain implementation.
9. An `IncoherentModel` wraps one primary coherent `SXRDCrystal` and presents
   the fit-parameter methods already consumed by `CTROptimizer`.
   `IncoherentF2Model` adds `F2(h, k, l)`. The kinematical optimizer accepts
   either the coherent crystal or an `IncoherentF2Model` as its first argument;
   it must not contain Poisson-specific branches.
10. Arbitrary lists of incoherent models are not accepted because sequential
   incoherent averaging has no generally valid physical meaning. A joint model
   must make the correlations between its state variables explicit.
11. Height-state enumeration has a fixed, non-fit configuration setting
    `exact_layer_count`, initially defaulting to 10. If the candidate
    distribution contains no more than this number of layers, retain every
    state. For a wider distribution, choose a contiguous retained interval
    from the calculated Poisson probability masses and the profile's
    cumulative tail-probability target. Do not use measured CTR values or
    uncertainties to choose the forward-model support. Renormalize the
    retained probabilities to sum to one and report the excluded probability
    mass. `exact_layer_count` is a policy switch, not a hard cap: more than ten
    states may be retained when the probability-mass criterion requires them.
    Here "exact" means that every state in the profile's finite candidate
    support is evaluated; it does not claim evaluation of an infinite Poisson
    tail.

12. Height states are indexed by the structural layer `n` of the top filled
    layer, and the mass of that state is `PoissonProfile.probability(n + 1)`.
    The `+ 1` follows from `occupancy(n) = P(H > n)`: layer `n` is the top
    filled layer exactly when the signed height change equals `n + 1`.
13. `CTROptimizer.startp` keeps its existing meaning as the model-block
    preparation snapshot. This feature does not redefine it into a full
    optimizer vector, because that would silently change its length for every
    existing fit which uses callbacks or fitted resolution.

## Existing coherent behavior

`PoissonProfile.occupancy` describes cumulative material occupancy at each
structural-layer offset, `occupancy(n) = P(H > n)` for the signed height
change `H`. `surface_occupancy` takes the difference between successive
cumulative occupancies and therefore supplies the exposed fraction of each
height state.

Two consequences fix the index convention used throughout this record. Layer
`n` is the top filled layer exactly when `H == n + 1`, so the exposed mass of
state `n` is `probability(n + 1)`, not `probability(n)`. And
`surface_occupancy` is not merely a stylistic alternative to that expression:
the two agree to floating-point roundoff on every interior bin, but
`surface_occupancy` sets its terminal bin to `occupancy[-1]` and therefore
folds the complete upper tail into the highest represented height.

`PoissonSurface.createLayers` currently assigns these fractions to coherent
domain occupancies. It also assigns the cumulative rough-Film correction and
the termination-specific replacement of exposed Film material. `F_uc` sums all
of those contributions as one complex amplitude. `SXRDCrystal.F` then adds the
bulk and every top-level component before the caller takes an absolute value.

For flat-height amplitudes `A_n(Q)` and the coherent path's exposed fractions
`q_n`, the current calculation is equivalent, up to the configured Poisson
tail truncation, to

```text
A_coherent(Q)  = sum_n q_n A_n(Q)
F2_coherent(Q) = abs(A_coherent(Q)) ** 2.
```

`q_n` is written separately from the ensemble probability `p_n` used from the
next section onwards, and the two are not interchangeable. `q_n` is what
`createLayers` assigns today: `surface_occupancy` masses, upper tail folded
into the terminal bin, left unnormalized after thresholding. `p_n` is the
incoherent ensemble's `probability(n + 1)` mass over the retained interval,
renormalized to one. They agree on every interior retained bin and differ at
the boundary. That difference is exactly why the live coherent evaluation, and
not the finite state-sum reconstruction, is the authoritative `kappa = 0`
endpoint.

This is the characteristic-function form used by the classical coherent CTR
roughness treatments. Harada writes the damping factor as the squared modulus
of the Fourier sum over relative areas at each step height. Dale *et al.* use
the same squared expectation for discrete height distributions. The latter
paper calls the measured scattering incoherent in its introductory wording,
but its height-distribution factor is still the coherent amplitude average in
the terminology of this design.

## Height-domain ensemble

Let `p_n` be the represented probability of height state `n`, and let `A_n` be
the complete coherent amplitude of the crystal when that height state covers
the local coherent patch. Each `A_n` includes:

- the semi-infinite bulk amplitude;
- the underlying Film and all other common components below the target
  surface;
- the flat Film-height correction for state `n`;
- the termination-specific surface slab for state `n`, including its
  interference with the bulk and Film;
- crystal area scaling, component weight, and supported coherent-domain
  transforms.

The limiting squared structure factors are

```text
F2_coherent   = abs(sum_n p_n A_n) ** 2
F2_incoherent = sum_n p_n abs(A_n) ** 2.
```

The second expression is the large-domain limit: each coherence patch sees one
flat height and their scattering strengths add without cross terms. It is not

```text
abs(F_bulk_and_film) ** 2 + sum_n p_n abs(F_surface_n) ** 2,
```

because that incorrect expression removes bulk--surface and Film--surface
interference inside each domain.

### Partial coherence

Use

```text
F2_kappa = (1 - kappa) F2_coherent_live
           + kappa F2_incoherent_retained,
```

where `kappa = 0` is fully coherent and `kappa = 1` is fully
incoherent. `incoherent_fraction` is preferred over an unqualified
`coherence` parameter so that the endpoint convention is visible at every
call site.

For an ideal complete normalized height distribution,
`F2_coherent_live = abs(sum_n p_n A_n) ** 2` and the same interpolation can be
written

```text
F2_kappa = F2_coherent
           + kappa * (
               sum_n p_n abs(A_n) ** 2
               - abs(sum_n p_n A_n) ** 2
             ).
```

Equivalently in that complete-support limit, the height-state coherence
matrix is

```text
C = (1 - kappa) outer(p, p) + kappa diag(p),
F2_kappa = A.conj() @ C @ A.
```

This matrix is positive semidefinite for valid probabilities and
`0 <= kappa <= 1`, so the model cannot create negative `F2`.

This equality is exact for a complete normalized height distribution. The
finite numerical ensemble renormalizes its retained Poisson probability
masses, so the state weights used for the incoherent endpoint also form a
normalized distribution. In production, the live coherent result is the
authoritative `kappa=0` endpoint; the state-sum reconstruction is a diagnostic
which must agree within a Q-dependent truncation tolerance justified from the
excluded mass and a state-amplitude envelope, or from an extended-support
convergence check. Excluded probability alone is not an amplitude-error
bound. Because the live coherent endpoint is not replaced by the finite
state-sum reconstruction, the coherence-matrix expression above is the ideal
complete-support interpretation; the implemented partial model is the convex
interpolation of the exact live coherent endpoint and the normalized finite
incoherent endpoint.

The interpolation follows from averaging finite coherence patches. If `f_n`
is the fraction of height `n` inside one patch, `E[f_n] = p_n`, and

```text
Cov(f_n, f_m) = kappa * (p_n delta_nm - p_n p_m),
```

then `E[abs(sum_n f_n A_n) ** 2]` is exactly `F2_kappa`. For a simple patch
containing `N_eff` independent, equal-area domains, `kappa = 1 / N_eff`.
This interpretation gives the required limits without claiming that one
universal formula maps a quoted coherence length to `kappa`.

### Coherence-length extension

The Poisson profile is a one-point height distribution. A physical
partial-coherence calculation additionally requires the lateral two-point
correlation of the domain field. On the CTR, a reduced scalar can be written
schematically as

```text
             integral d2r K_Q(r) rho_height(r)
kappa(Q) =  -----------------------------------,
                     integral d2r K_Q(r)
```

where `rho_height` is the normalized surface-domain correlation and `K_Q`
contains the incident mutual coherence, illuminated-footprint autocorrelation,
and reciprocal-space acceptance. With this convention:

- domains much larger than the effective coherence area give `kappa -> 1`;
- many independent domains within a coherence area give `kappa -> 0`.

The projected transverse coherence is generally anisotropic and can vary with
incidence/exit geometry. Detector acceptance and bandwidth can also make the
effective mixture reflection dependent. Therefore:

- the first implementation accepts a scalar `incoherent_fraction`;
- a later `SurfaceCoherenceKernel` may return `kappa` for each requested point;
- domain correlation lengths belong to a sample-side domain model;
- beam coherence and acceptance belong to CTR measurement/resolution metadata;
- a CTR-only scalar fit cannot separately identify domain size and beam
  coherence length. AFM, transverse scans, rocking widths, or reciprocal-space
  maps are needed to constrain one of them.

Longitudinal coherence is omitted initially. For atomic step-height
differences it will normally be much longer than the corresponding optical
path difference. If it becomes relevant, it must damp individual off-diagonal
height-state terms as a function of both `(n, m)` and scattering geometry;
that cannot in general be represented by one scalar `kappa`.

## Why the model returns `F2`, not intensity or reflectivity

The mathematical incoherent operation squares the complete coherent amplitude
of each state and removes the appropriate cross terms. The public quantity
must still identify which forward amplitude was squared.

For the kinematical CTR model, `SXRDCrystal.F` is the reference-lateral-cell-
normalized structure factor. Its squared quantity is therefore

```text
F2 = abs(F) ** 2,
```

in squared structure-factor units (nominally electrons squared where `F` is in
electrons). It is not observable detector counts: incident flux, footprint,
polarization/Lorentz factors, detector response, acquisition time, background,
and the fitted experimental scale are not part of `F2`.

It is also not reflectivity. Reflectivity is a dimensionless intensity ratio
formed from the optical reflection amplitude. The quantities are:

| Forward model | Coherent state amplitude | Domain-averaged output |
|---|---|---|
| Kinematical CTR | `F_total_n` | `F2 = sum_n p_n * abs(F_total_n) ** 2` |
| DWBA reflectivity | `r_total_n = r_0_n + r_h_n` | `R = sum_n p_n * abs(r_total_n) ** 2` |
| Isolated DWBA matrix-element diagnostic | `F_h_n` | not a reflectivity observable |

The first implementation is specifically an incoherent `F2` model. It must
not call `SXRDCrystal.specular_reflectivity`, which is an optical profile
calculation unrelated to the kinematical `SXRDCrystal.F` path.

`F2` must never change meaning according to an optimizer mode. A future
incoherent DWBA implementation should expose a separately named `R` boundary,
or a typed forward-result object whose quantity is explicitly `"R"`; it must
not return reflectivity from a method named `F2`.

DWBA cannot use the proposed common-prefix optimization without further
physics. Changing the flat surface height changes the one-dimensional optical
reference profile, its Fresnel amplitude `r_0`, its internal fields, and the
atomic/reference decomposition. A physically correct DWBA domain ensemble
would prepare and evaluate each height state separately and then mix the
resulting `abs(r_total_n) ** 2`. Until that exists, combining an incoherent
surface `F2` model with DWBA must raise a clear `NotImplementedError` or
`ValueError`; silently reverting to coherent behavior is forbidden.

## API design

### Module and responsibility boundary

Add `orgui.datautils.xrayutils.CTRincoherent`. Keep `SXRDCrystal` as the
coherent amplitude model and make each incoherent model a squared-structure-
factor wrapper around one coherent crystal. The module contains:

- the abstract `IncoherentModel` wrapper/parameter contract;
- the `IncoherentF2Model` kinematical quantity contract;
- shared parameter metadata and validation;
- the model-type registry used by configuration/persistence;
- a lazy kinematical evaluation context which caches coherent components;
- reusable state-ensemble mixing helpers; and
- `PoissonHeightDomains`, the first concrete model.

The kinematical optimizer receives either an `SXRDCrystal` or an
`IncoherentF2Model` as the existing first positional model argument. The
registry is not a service locator in the numerical hot path. It exists so a
saved type name can be converted to a class and so external packages can add
models without editing an `if/elif` chain.

Callers that do not supply an incoherent model continue to use
`abs(crystal.F(...)) ** 2`. Thus old scripts, saved crystal files, optimizer
construction, and direct `F` calls remain coherent without a compatibility
switch.

Add `SXRDCrystal.F2(h, k, l)` as a thin coherent convenience returning
`abs(SXRDCrystal.F(h, k, l)) ** 2`. The coherent crystal and incoherent wrapper
then share an `F2` evaluation boundary, while only the coherent crystal exposes
the phase-bearing `F` method. Existing callers may continue taking `abs(F) ** 2`.

### Contract class

Use an abstract base class rather than only a `Protocol`. The contract should
reuse the fit API already implemented by `SXRDCrystal` and
`CTRutil.LinearFitFunctions`, rather than inventing an `expose_parameter`
dialect. The quantity-neutral base class owns one primary coherent crystal,
model-local fit parameters, and the concatenation of both parameter blocks.
Quantity-specific subclasses define numerical evaluation methods.

The planned public surface is:

```python
class IncoherentModel(LinearFitFunctions, ABC):
    model_type: ClassVar[str]
    supported_forward_models: ClassVar[frozenset[str]]

    def __init__(self, crystal: SXRDCrystal, *, name: str = "incoherent"): ...

    @property
    def coherent_model(self) -> SXRDCrystal: ...

    def addFitParameter(
        self,
        parameter,
        limits=(-np.inf, np.inf),
        **keyargs,
    ) -> Parameter: ...

    @property
    def fitparnames(self) -> list[str]: ...

    @property
    def priors(self) -> list[object]: ...

    def getInitialParameters(self, force_recalculate=False) -> np.ndarray: ...
    def getStartParamAndLimits(
        self,
        force_recalculate=False,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]: ...
    def setParameters(self, values: ArrayLike) -> None: ...
    def setFitParameters(self, values):  # raises; use setParameters
        ...
    def getFitErrors(self) -> np.ndarray: ...
    def setFitErrors(self, errors: ArrayLike | None) -> None: ...
    def parameter_list(self) -> list[Parameter]: ...
    def parametersToDict(self) -> dict: ...
    def parametersFromDict(self, data, override_values=True) -> None: ...
    def clearParameters(self) -> None: ...
    def validate(self, *, forward_model: str) -> None: ...
    def to_config(self) -> dict: ...


class IncoherentF2Model(IncoherentModel, ABC):
    output_quantity: ClassVar[str] = "F2"

    def F2(self, h, k, l) -> np.ndarray: ...
```

The model-local parameter API follows existing usage exactly:

```python
model.addFitParameter(
    "incoherent_fraction",
    limits=(0.0, 1.0),
    name="rough_surface incoherent_fraction",
)
```

`addFitParameter` itself is the inherited `LinearFitFunctions` implementation,
which already resolves a string through `parameterLookup` and rejects a name
already present in `self.fitparnames`. Because `fitparnames` is a composite
here, that inherited check spans the wrapped crystal's names too, so the
invariant that local and coherent names stay unique after concatenation is
enforced at the point of registration rather than at `prepareFit`.

`PoissonHeightDomains` defines `parameterLookup` and `parameterLookup_inv` for
its local numerical basis, just as the existing classes do. Its initial basis
contains `incoherent_fraction`; calling `addFitParameter` promotes that basis
entry into the fitted vector. A setting therefore remains fixed unless the
caller explicitly adds it as a fit parameter. Structural/numerical policy
settings such as `surface` and `exact_layer_count` are constructor/configuration
fields, not basis entries, and never appear in the optimizer parameter vector.

The wrapper's public parameter order is always:

```text
[incoherent-model-local parameters] [wrapped coherent-crystal parameters]
```

`getStartParamAndLimits`, `getInitialParameters`, `fitparnames`, `priors`,
`parameter_list`, `setParameters`, `setFitErrors`, and `getFitErrors` all use
that same order. The base class splits and delegates the coherent tail to the
wrapped `SXRDCrystal`. This makes the wrapper directly compatible with the
optimizer's existing model contract.

Following the existing composite-model convention, `parametersToDict` returns
both a model-local subtree and the wrapped crystal's parameter subtree;
`parametersFromDict` restores both, and `clearParameters` clears both. The
initial dictionary configuration stores only the local subtree because the
factory receives the coherent crystal separately, but that is a configuration-
layer choice rather than a different public parameter API.

The contract has these strict invariants:

- values, lower bounds, upper bounds, names, priors, parameter records, and
  errors all have the same total length and stable order;
- local parameter semantics, default names, limits, and priors follow
  `LinearFitFunctions.addFitParameter`;
- local and wrapped-crystal parameter names are unique after concatenation;
- starts lie within limits and updates validate their complete input shape
  before changing either local or coherent state;
- `setFitErrors(None)` clears errors in both parameter blocks;
- after errors are cleared, `getFitErrors()` follows the existing convention
  and raises `ValueError` rather than returning `None`;
- `validate` resolves referenced components by stable name on
  `self.coherent_model` and checks calculation-mode capabilities; and
- `IncoherentF2Model.F2` returns a real, finite, nonnegative squared structure factor
  broadcast-compatible with `h`, `k`, and `l`. It never returns `sqrt(F2)`,
  detector counts, reflectivity, or a complex amplitude.

The base class should implement the parameter methods as template methods and
delegate validation to `_validate_model(forward_model)`.
`IncoherentF2Model` implements the public `F2` validation boundary and delegates
numerical work to `_evaluate_F2(context)`. This gives third-party models one
conformance path without making them reproduce vector validation and coherent-
tail delegation.

### Evaluation context and reusable state ensembles

The generic model contract must not be tied to `PoissonSurface`. Its numerical
input is a lazy `KinematicIncoherentContext` containing the wrapped coherent
crystal, the requested HKL arrays, and an evaluation-local cache:

```python
context = KinematicIncoherentContext(self.coherent_model, h, k, l)
context.coherent  # lazy KinematicAmplitudeResult
context.component("rough_surface")
complete_states = context.iter_component_states(
    "rough_surface",
    state_evaluator,
)
```

`IncoherentF2Model.F2(h, k, l)` constructs this context and calls the protected
model hook. `CTROptimizer` uses the same entry point, including inside
resolution quadrature. The context centralizes reference-area scaling,
component weights, outer coherent-domain transforms, shapes, and caching; an
incoherent model must not duplicate them.

`state_evaluator` is called by the context at every required transformed HKL;
it is not an amplitude precomputed only at the original HKL. This is necessary
because an outer coherent-domain transform changes the coordinates at which a
component state must be evaluated. The public production iterator is
`iter_component_states`: each yield contains one layer number, one normalized
probability, and the complete amplitude array for that height. One-height
streaming is the initial memory contract; batching is deferred until profiling
shows that it is needed.

Provide an optional `CoherentStateEnsembleModel(IncoherentF2Model)` base for
the common case in which a model generates probabilities `p_n` and complete
amplitudes `A_n`. It implements validated coherent, incoherent, and coherence-
matrix sums. `PoissonHeightDomains` subclasses it. A model with different
kinematical physics may subclass `IncoherentF2Model`; a future quantity such as
reflectivity starts from the quantity-neutral `IncoherentModel` instead.

This admits, without optimizer changes, future registered models for discrete
chemical terminations, lateral phase domains, orientational/mosaic domains,
or a physically calculated mutual-coherence matrix. Those models may expose
zero, one, or many fit parameters through the same block contract.

The optimizer accepts one root model, not a list. For example, two rough
interfaces cannot generally be evaluated by applying two independent mixers
in sequence. Their joint states need `p_nm`, or a documented independence
factorization, and their partial coherence needs the corresponding joint
coherence matrix. Such behavior belongs in a registered joint/composite model.

### Model registration

Registration maps the class's single stable `model_type` persistence key to
one contract-compliant class:

```python
@register_incoherent_model
class PoissonHeightDomains(CoherentStateEnsembleModel):
    model_type = "poisson_height_domains"

    def __init__(
        self,
        crystal,
        *,
        surface,
        incoherent_fraction,
        exact_layer_count=10,
        name="incoherent",
    ): ...
    ...

model = create_incoherent_model(
    "poisson_height_domains",
    crystal,
    surface="rough_surface",
    incoherent_fraction=0.35,
)
```

Direct construction remains the ordinary Python API. Registration must reject
duplicate keys, abstract classes, classes not derived from `IncoherentModel`,
concrete classes without a stable `output_quantity`, and built-in key
replacement. Public discovery uses immutable `IncoherentModelInfo` records
containing the key, class, short description, supported forward models, and
`output_quantity`. `available_incoherent_models()` returns a read-only mapping
from registry key to that record, so key membership and metadata lookup are
unambiguous.
Deserialization may instantiate only registered keys; it must not import an
arbitrary class path from a file. The dictionary round-trip entry point is
`create_incoherent_model_from_config(crystal, config)`; it validates the
registered type, reconstructs the non-basis settings, and restores only the
wrapper-local parameter subtree. Full GUI/session wiring is deferred.

Third-party registration is process-local. It need not use packaging entry
points in the first implementation. Entry-point discovery can be added later
without changing the optimizer contract.

### `CTROptimizer` contract

Use the existing optimizer constructor shape. The incoherent wrapper is the
forward model:

```python
fit = CTROptimizer(
    domains,
    ctrs,
)
```

`CTROptimizer` deep-copies the complete forward-model argument once. Store that
owned object as `self.model`. For a wrapper, set `self.xtal` to
`self.model.coherent_model`; for an ordinary coherent fit, `self.model` and
`self.xtal` are the same object. This preserves the established meaning of
`optimizer.xtal`, and existing callbacks and crystal constraints continue to
receive an `SXRDCrystal`. Component selection inside the copied crystal is by
stable name, not a cached component identity.

No `incoherent_model=` constructor keyword or setter is needed. To change
between coherent and incoherent fitting, construct the optimizer with either
the original `SXRDCrystal` or an `IncoherentF2Model` wrapping it. This keeps
one authoritative forward model and prevents mismatched crystal/model pairs.

The complete fit-vector layout is:

```text
[resolution] [callbacks] [subclass blocks]
[incoherent-model-local parameters] [wrapped crystal parameters]
```

For the base `CTROptimizer`, `subclass blocks` is empty. For
`CTROptAngleCorrection`, its existing phase/amplitude block remains in that
position. Relative ordering of every pre-existing parameter is unchanged when
the first model argument is an ordinary `SXRDCrystal`, and crystal parameters
remain the tail. The incoherent wrapper supplies its combined local/crystal
block through the methods the optimizer already calls; the optimizer does not
need Poisson-specific parameter slicing.

At `prepareFit`, the optimizer checks each block and the concatenated vector:

- wrapper values, bounds, names, parameter records, and priors agree in length;
- full optimizer values, bounds, and names agree in length;
- the full parameter names are unique;
- all starts satisfy bounds;
- wrapper validation succeeds for its coherent crystal; and
- the selected calculation consumes the wrapper's declared `output_quantity`;
  an `F2` wrapper is rejected when DWBA is enabled.

The optimizer's `get_parameters`, `set_parameters`, `get_bounds`, `set_errors`,
`fitparnames`, and `priors` plumbing uses `self.model`, so it sees the wrapper
as one fitted model. Physics callbacks and crystal constraints use
`self.xtal`. The wrapper performs the local/coherent parameter split
internally. For a Poisson model with one fitted fraction and two crystal
parameters, and no other optional blocks, the observable contract is:

```text
n_parameters                 = 3
fitparnames                  = ["rough_surface incoherent_fraction",
                                <crystal name 0>, <crystal name 1>]
get_parameters()             = [0.35, <crystal value 0>, <crystal value 1>]
get_bounds()[0]              = [0.0,  <crystal lower 0>, <crystal lower 1>]
get_bounds()[1]              = [1.0,  <crystal upper 0>, <crystal upper 1>]
```

`CTROptimizer.n_parameters` is the length of this frozen vector and is valid
after `prepareFit`. `get_parameters()` is the authoritative live parameter
vector and reflects changes made to the model after preparation.

`startp` is deliberately left alone. It stays the preparation-time snapshot of
the **model block only**, exactly as `prepareFit` sets it today, and
`lower_bounds`/`higher_bounds` likewise stay the unprefixed model bounds while
`self.bounds` carries the full prefixed pair. Redefining `startp` as the full
vector would silently change its length for every existing fit which registers
a callback or fits resolution, which is a breaking change this feature has no
reason to make. For a wrapper, `startp` therefore has the wrapper's
local-plus-coherent length. Its docstring gains a note that
`get_parameters()` is the authoritative full vector; deprecating `startp`
outright is a separate API change. Before preparation, callers can use
`len(model.getInitialParameters())` for the wrapper's total, but should not
infer the optimizer total because callbacks, subclass blocks, and resolution
settings can add prefixes.

`prepareFit` records the model parameter-name/count/bounds signature. If
callers add, remove, rename, or re-bound fit parameters through
`optimizer.model` afterward, evaluation and parameter setters must raise that
`prepareFit()` is required; they must not slice a stale layout silently.

`CTROptimizer.priors` retains its existing model-block scope: for an
`IncoherentF2Model`, it is the wrapper's local-plus-coherent prior list. It does
not cover resolution, callback, or angle-correction prefixes because those
APIs currently expose no priors. Extending priors to the full optimizer vector
is a separate API change and is not implied by `n_parameters`.

Direct calculation and resolution sampling use one canonical `F2` callback.
Conversion to the stored structure-factor representation is:

```text
predicted |F| = sqrt(resolution_operator(F2)).
```

Resolution is applied to `F2` before the square root. Existing helpers named
`CTRresolution.sample_intensity` and `fast_convolve_intensity` provide the
correct numerical boundary when their input is explicitly documented as
`F2`; their generic names do not change the quantity.
`sample_structure_factor` keeps its public name and `|F|` output, but should
prefer a model's `F2` method. It retains `abs(F) ** 2` as a compatibility
fallback for existing crystal-like objects which implement only `F`.

### End-to-end usage mockup

The following is a design mockup, not executable code yet. Assume `bulk_uc`,
`film_uc`, `termination_cells`, `h`, `k`, `l`, and `ctrs` have been constructed
by the existing APIs.

```python
import numpy as np

from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRdistributions import PoissonProfile
from orgui.datautils.xrayutils.CTRfilm import Film, PoissonSurface
from orgui.datautils.xrayutils.CTRincoherent import (  # planned module
    PoissonHeightDomains,
    available_incoherent_models,
    create_incoherent_model_from_config,
)
from orgui.datautils.xrayutils.CTRopt import CTROptimizer


# This remains an ordinary coherent structural model. Constructing a
# PoissonSurface does not opt into incoherent averaging.
film = Film(film_uc, name="film")
film.basis[0] = 12.0
surface = PoissonSurface(
    termination_cells,
    profile=PoissonProfile(mean_change=1.5, alpha=0.6, offset=0.0),
    name="rough_surface",
)
# Stacking levels are ordering keys and must be a numpy array; a plain list
# fails when the constructor reorders them. The values match the convention
# used by the existing tests for a Film plus surface pair.
crystal = SXRDCrystal(bulk_uc, film, surface, stacking=np.array([1, 2]))

F_coherent = crystal.F(h, k, l)
F2_coherent = crystal.F2(h, k, l)
np.testing.assert_allclose(F2_coherent, np.abs(F_coherent) ** 2)


# Explicit fully incoherent large-domain limit. Its local setting is fixed
# unless addFitParameter is called for that setting.
large_domains = PoissonHeightDomains(
    crystal,
    surface="rough_surface",
    incoherent_fraction=1.0,
    exact_layer_count=10,  # default; shown here to make the policy explicit
)
assert len(large_domains.getInitialParameters()) == len(crystal.fitparnames)
F2_incoherent = large_domains.F2(h, k, l)


# A partially coherent model with one optimizer-visible parameter.
mixed_domains = PoissonHeightDomains(
    crystal,
    surface="rough_surface",
    incoherent_fraction=0.35,
)
mixed_domains.addFitParameter(
    "incoherent_fraction",
    limits=(0.0, 1.0),
    name="rough_surface incoherent_fraction",
)
assert mixed_domains.fitparnames[0] == "rough_surface incoherent_fraction"

F2_mixed = mixed_domains.F2(h, k, l)
F_magnitude_mixed = np.sqrt(F2_mixed)  # only at a |F| API boundary


# The wrapper is the optimizer's model. The existing model-copy operation
# copies both the wrapper and its coherent crystal.
fit = CTROptimizer(
    mixed_domains,
    ctrs,
)
fit.prepareFit()
assert fit.xtal is fit.model.coherent_model
assert fit.n_parameters == len(fit.get_parameters())
prepared_start = fit.startp.copy()  # legacy model-block snapshot; the full
                                    # vector is get_parameters()
assert fit.fitparnames[-(len(crystal.fitparnames) + 1)] == (
    "rough_surface incoherent_fraction"
)
prediction = fit.flat_prediction()


# Registry construction is useful for a UI or dictionary configuration, but is not
# required for direct Python use.
assert "poisson_height_domains" in available_incoherent_models()
config = mixed_domains.to_config()
restored = create_incoherent_model_from_config(crystal, config)
```

Advanced model implementations can inspect the planned coherent decomposition
without reconstructing the crystal:

```python
evaluation = crystal.evaluate_kinematic(h, k, l)
surface_part = next(
    part for part in evaluation.components
    if part.name == "rough_surface"
)
common_amplitude = evaluation.total - surface_part.amplitude

height_states = surface.flat_domain_corrections(h, k, l)
height_states.layer_numbers
height_states.probabilities
height_states.raw_retained_probability
height_states.excluded_probability
for state in height_states.iter_states():
    state.layer_number
    state.probability
    state.amplitude

# Explicit diagnostic/test materialization, not the production hot path.
all_corrections = height_states.as_array()
```

`PoissonHeightDomains` intentionally has no `F` method. Requesting one is a
category error because the relative phase between incoherent patches is not
defined. The initial model advertises only kinematical support, so this fails
during preparation rather than silently applying the wrong averaging:

```python
fit.set_dwba(True)
fit.prepareFit()
# ValueError: poisson_height_domains does not support the DWBA forward model
```

### Persistence

Do not put incoherent-model settings or beam coherence in `.xtal` files.
`.xtal` describes the coherent sample structure, while the effective mixture
depends on the measurement coherence and resolution. The first implementation
provides a plain dictionary contract only:

```python
config = model.to_config()
# {
#     "type": "poisson_height_domains",
#     "settings": {
#         "surface": "rough_surface",
#         "exact_layer_count": 10,
#     },
#     "parameters": <local parametersToDict subtree>,
# }

restored = create_incoherent_model_from_config(crystal, config)
```

The local parameter subtree contains `basis_0` and any registered `Parameter`
records, so it is the single persisted source for `incoherent_fraction`
whether fixed or fitted. Do not duplicate the fraction under `settings`.
`exact_layer_count` is a non-fit numerical setting and therefore belongs in
`settings`.

The factory receives an already constructed coherent crystal, resolves `type`
through the registry, wraps the crystal, and restores the stored local subtree
through the same private local helper used by the wrapper's composite
`parametersFromDict`. The public `parametersFromDict` remains a full
local-plus-coherent round trip. This increment does not choose or modify a GUI
session format; a session layer may later store this dictionary or omit it for
coherent behavior. A future physical kernel may refer separately to sample
domain-correlation parameters and measurement coherence metadata.

## Reusable kinematical decomposition

Strict separation of the public models does not require repeated evaluation of
the common crystal. Refactor `SXRDCrystal.F` around an immutable result:

```python
@dataclass(frozen=True)
class KinematicComponentAmplitude:
    index: int
    name: str
    amplitude: complex | np.ndarray


@dataclass(frozen=True)
class KinematicAmplitudeResult:
    bulk: complex | np.ndarray
    components: tuple[KinematicComponentAmplitude, ...]
    total: complex | np.ndarray
```

Proposed public method:

```python
result = crystal.evaluate_kinematic(h, k, l)
```

Its rules are:

- input coordinates and output units are identical to `SXRDCrystal.F`;
- `bulk` includes reference-area scaling;
- every component amplitude includes its reference-area scaling, crystal
  weight, and outer coherent-domain transforms;
- `total == bulk + sum(component.amplitude for component in components)`;
- returned arrays are read-only or otherwise documented as caller-owned;
- `SXRDCrystal.F` becomes a thin wrapper returning `.total`;
- the refactor must be bitwise equal where practical and numerically equal in
  every existing CTR test.

The decomposition is useful beyond this feature, but its first consumer is the
`CTRincoherent` evaluation context. For a target surface at index `s`, it
obtains the shared
amplitude once:

```text
A_common = result.total - result.components[s].amplitude.
```

It then combines `A_common` with each flat-height correction of the target
surface. A raw surface correction is not yet a crystal contribution. For outer
crystal domain `d`, the state contribution must be

```text
Delta A_n = area_scale * component_weight
            * sum_d outer_domain_occupancy_d
                    * Delta F_n(M_d @ hkl).

A_n = A_common + Delta A_n.
```

Factor the existing area/weight/domain loop into a shared evaluator used by
both `evaluate_kinematic` and the context's component-state operation. Do not
duplicate this logic in `PoissonHeightDomains`. This avoids recalculating the
bulk and ordinary Film components for every Poisson state without dropping an
established crystal-level transform.

### Flat-height surface evaluation

Add a `PoissonSurface` method returning coherent corrections rather than
complete crystal amplitudes:

```python
states = surface.flat_domain_corrections(
    h,
    k,
    l,
    exact_layer_count=10,
)

states.layer_numbers   # structural-layer offsets
states.probabilities   # normalized p_n for the retained Poisson support
states.raw_retained_probability
states.excluded_probability
states.excluded_lower_probability
states.excluded_upper_probability
states.iter_states()   # one correction Delta F_n at a time
states.as_array()      # explicit optional materialization for diagnostics
```

Each returned correction for state `n` contains the flat Film-height change and
the termination replacement for that height, both at occupancy one. It must
not contain `p_n`; probability is applied only by
`CoherentStateEnsembleModel`.

For sorted height states, calculate the Film part cumulatively:

1. Evaluate every distinct generated Film layer amplitude once.
2. Use a prefix sum starting at `min(0, lowest retained state)` -- the sharp
   Film boundary, or below it when the surface is etched -- to construct the
   flat Film correction for every state. Starting at the lowest retained state
   instead would drop the Film layers between the boundary and that state.
3. Evaluate each termination replacement once per termination-cycle state and
   translated height.
4. Add the appropriate termination replacement to each Film prefix.

This changes the expensive scaling from repeated construction of a complete
crystal for every height to one common-crystal evaluation plus one pass over
the represented layers and height states.

The first implementation should reuse the current layer-cycle,
growth/etching, offset, strain, termination-bank, and profile support logic. It
must not rederive layer numbering independently. It must, however, obtain the
height-state masses as `PoissonProfile.probability(layer_numbers + 1)` rather
than from `surface_occupancy`; otherwise the finite-support tail is folded
into the highest flat-height state before the incoherent model applies its own
cutoff policy. At `kappa = 1` that folded bin would become a single spurious
flat state carrying the whole tail mass, which is a physically wrong result
that still looks plausible. Because the interior bins of the two expressions
agree to roundoff, a regression test must assert both halves of the relation:
interior agreement, and deliberate disagreement at the terminal bin.

The retained coherent component set and the retained height-state ensemble are
different objects and must not be conflated. `createLayers` retains the union
of exposed-surface and Film-correction layers and never renormalizes; the
incoherent ensemble retains exposed-height states only, and renormalizes. A
third range is independent of both: the Film correction of flat state `n`
covers every structural layer between the sharp boundary and `n`, so the
cumulative Film prefix must span
`[min(0, lowest_retained_state), highest_retained_state]` taken from the full
candidate support rather than from either retained mask.

`PoissonHeightDomains.exact_layer_count` defaults to 10 and is a positive
integer, non-fit setting. If the candidate support contains at most this many
nonzero states, retain them all. For a wider support, retain a contiguous
interval containing the mode of the calculated masses -- `argmax` of
`probability(layer_numbers + 1)`, ties resolved to the lower index, never
`floor(rate)` -- and at least `exact_layer_count` states, then expand it using
the adjacent calculated probability masses until the profile's cumulative
`tail_probability` target is met. This is a
probability-data cutoff, not a dependency on observed CTR values. It does not
impose a ten-state cap.

Normalize the retained masses so `sum(states.probabilities) == 1` within
floating-point tolerance. Record the raw retained mass and excluded lower and
upper tail masses before normalization. Compare coherent reconstructions with
a Q-dependent tolerance derived from an amplitude envelope or an
extended-support convergence check, not directly with the dimensionless tail
probability.

### Evaluation lifetime and cache invalidation

The mandatory cache is an evaluation-local decomposition: for one requested
HKL batch, bulk and common components are calculated once and reused for every
height state. This removes the multiplicative height-state overhead without a
long-lived invalidation hazard.

Do not add a global array cache in the first increment. Optimizer parameters,
component weights, coherent-domain transforms, stacking, attenuation,
reference-cell transforms, energy, and HKL arrays can all invalidate it. If
profiling later shows that cross-iteration reuse matters, add component
fingerprints or mutation generations before retaining arrays across calls.
Object identity or only `_basis_created` is not an adequate cache key.

### Stacking boundary

The first implementation targets one `PoissonSurface` which is the highest
material component in stacking order. This is the case in which all material
below the surface is common and only the surface correction varies.

If a component is stacked above the target surface, its absolute position or
structure may depend on the selected height. Treating it as part of
`A_common` would be wrong. Initially reject that topology with an actionable
error. A later extension can split the evaluation into:

- a cached common prefix below the target surface;
- the target flat-height correction;
- a state-dependent suffix restacked and evaluated for each height.

Multiple incoherently rough interfaces also require a joint height-state
distribution or an explicit independence assumption and are out of the first
scope.

## Validation invariants

### Scientific endpoint tests

1. **No opt-in, no change:** existing `SXRDCrystal.F`, optimizer predictions,
   and resolution calculations remain numerically unchanged.
2. **Coherent reconstruction:** `sum(p_n * A_n)` agrees with the existing
   coherent `PoissonSurface` result within a Q-dependent truncation tolerance
   justified by an extended-support reference or amplitude-weighted tail.
3. **Fully incoherent endpoint:** `kappa=1` agrees with a direct weighted sum
   of `F2` values from independently constructed deterministic flat-height
   crystals.
4. **Deterministic height:** when only one `p_n` is nonzero, predictions are
   independent of `kappa`.
5. **Two-height anti-Bragg case:** equal populations of two height states whose
   amplitudes differ by a sign cancel coherently; the fully incoherent result
   remains the single-flat-state `F2`, and the partial model fills the
   minimum linearly with `kappa`. This is reachable with a real profile rather
   than a synthetic ensemble: `PoissonProfile(mean_change=0.5, alpha=0.0)` has
   `rate == 0` and a one-half deterministic step fraction, so it populates
   structural layers `-1` and `0` at exactly `0.5` each, and there
   `probability(n + 1)` and `surface_occupancy(n)` coincide exactly.
6. **Nonnegative `F2`:** random valid amplitudes, probabilities, and `kappa`
   values never produce negative `F2` beyond roundoff.
7. **Complete-state interference:** a fixture with a nonzero bulk amplitude
   distinguishes the correct `sum(p_n * abs(A_n)**2)` result from the
   incorrect sum of isolated surface squared amplitudes.

### Existing surface behavior

Exercise all of the current high-risk cases:

- positive growth and negative etching;
- fractional deterministic step and nonzero `offset`;
- multi-layer Film cycles and termination-specific slabs;
- non-unit surface/Film out-of-plane repeat ratios;
- coherent-domain transforms and occupancies;
- reference-area scaling and component weights;
- Poisson probability cutoff, normalization, and tail metadata;
- copied crystals, fitted termination parameters, and coupled parameters.

### Observable and resolution tests

- `SXRDCrystal.F2(h, k, l)` equals `abs(SXRDCrystal.F(h, k, l)) ** 2`
  and retains the reference-cell normalization of `F`.
- No resolution: stored structure-factor prediction is `sqrt(F2_kappa)`.
- Fast convolution: convolve `F2_kappa`, then take the square root.
- Quadrature sampling: evaluate `F2_kappa` at every sampled HKL, integrate, then
  take the square root.
- Attempting to use the Poisson model with DWBA raises explicitly before any
  mixed prediction is returned.
- Attempting to request `F` from the incoherent model is unsupported by API,
  not implemented as `sqrt(F2)` with a fabricated phase.

### Contract and optimizer tests

- A dummy registered wrapper with zero local fit parameters contributes only
  its wrapped crystal's vector, bounds, names, priors, and errors.
- Dummy models with one and several local fit parameters prove that every
  contract array has the declared length and stable order.
- The complete optimizer vector has
  `[resolution][callbacks][subclass][incoherent][crystal]` order, while a fit
  with an ordinary coherent crystal retains its exact existing layout.
- Optimizer `set_parameters` and `set_errors` deliver the combined model block
  to the wrapper; the wrapper sends only its coherent tail to `SXRDCrystal`.
- Bad shapes, duplicate names, inverted bounds, starts outside bounds,
  nonfinite values, invalid errors, and missing model priors for unbounded
  parameters fail during `prepareFit` with model and parameter names in the
  error.
- Optimizer construction deep-copies the model. Mutating the caller's original
  model afterward does not alter the optimizer-owned model.
- Registry lookup and configuration round-trips preserve settings, local
  parameter order, bounds, current values, and priors supported by the existing
  `Parameter.asdict` format. Duplicate or unknown type keys fail
  deterministically.
- A minimal third-party test model can be registered and fitted without
  changes to `CTRopt.py`, proving that the optimizer has no Poisson-specific
  branch.

## Implementation plan

Implement this as a sequence of reviewable increments. Each increment has a
test gate; do not start optimizer integration until the coherent decomposition
and deterministic-height oracle agree.

### Preparatory fix -- error routing in `evaluateStatistics`

Files:

- `orgui/datautils/xrayutils/CTRopt.py`
- `orgui/datautils/xrayutils/test/test_CTRopt.py`

Work:

1. `CTROptimizer.evaluateStatistics` calls `self.xtal.setFitErrors(errors[3:])`
   or `self.xtal.setFitErrors(errors)` directly. That hardcodes the resolution
   prefix width and skips registered callbacks entirely, so with a callback
   present the crystal receives the callback's error slice. Replace both
   branches with `self.set_errors(errors)`, which is the existing splitter for
   resolution, callbacks, subclass blocks, and the model.
2. This is a pre-existing defect independent of incoherent models, so it lands
   as its own commit ahead of increment 1: bisectable on its own, and the
   wrapper work does not carry an unrelated behavior change.
3. It also composes with increment 5 at no extra cost. Once `_set_model_errors`
   forwards to `self.model`, a wrapper's combined local-plus-coherent error
   block reaches the wrapper through the same splitter, with no second site to
   update.
4. The repository AGENTS.md places error propagation under the `phys` scope, so
   the commit is `fix(phys)`. The method is already `DeprecationWarning`-marked
   and stays that way.

Gate:

- A fit with a registered callback and fitted resolution routes each error
  slice to its owner, asserted per block.
- With no callbacks and no fitted resolution, `evaluateStatistics` behavior is
  numerically unchanged.

### Increment 0 -- characterize the existing coherent model

Files:

- `orgui/datautils/xrayutils/test/test_CTRcalc.py`
- new `orgui/datautils/xrayutils/test/_poisson_oracle.py` for the shared
  deterministic flat-height fixtures. The repository has no `conftest.py` and
  its CTR tests are `unittest`-style classes run under pytest, so the shared
  fixtures are a plain importable module in the test package, not pytest
  fixtures.

Work:

1. Add a small Film plus `PoissonSurface` fixture with a short retained height
   support and at least two distinct termination cells.
2. Construct independent deterministic flat-height crystals for every retained
   height from explicit Film layers and the matching termination component.
   These are the correctness oracle; they must use neither `PoissonSurface`
   occupancy assembly nor the new state-extraction code.
3. Record coherent results for positive growth, negative etching, fractional
   deterministic steps, nonzero offset, multi-layer cycles, component weight,
   reference-area scaling, and an outer coherent-domain transform.
4. Add the equal-population, two-height anti-Bragg fixture and a fixture with a
   nonzero bulk amplitude which detects accidental averaging of only the
   surface correction.

Gate:

- Production code is unchanged.
- All current `test_CTRcalc.py` tests and the new characterization fixtures
  pass.

### Increment 1 -- add the common coherent `F2` boundary

File:

- `orgui/datautils/xrayutils/CTRcalc.py`

Work:

1. Add immutable `KinematicComponentAmplitude` and
   `KinematicAmplitudeResult` records.
2. Factor the existing component area-scale, weight, and outer-domain loop into
   one private evaluator. It must evaluate transformed HKL coordinates before
   applying each domain occupancy.
3. Add `SXRDCrystal.evaluate_kinematic(h, k, l)`. Its result contains the
   scaled bulk amplitude, one scaled amplitude per top-level component in
   stable crystal order, and their total.
4. Make `SXRDCrystal.F` a thin wrapper returning the result's total.
5. Add `SXRDCrystal.F2` returning `abs(F) ** 2`, documented as squared
   structure factor in the squared units of `F`, not counts or reflectivity.
6. Preserve scalar/array broadcasting, stacking, attenuation, reference frame,
   and contiguous-array behavior already exercised by the CTR tests.

Gate:

- `evaluate_kinematic(...).total` is numerically equal to the pre-refactor
  `F(...)` for all increment-0 fixtures.
- `F2(...) == abs(F(...)) ** 2` for scalar and array HKL input.
- The full `test_CTRcalc.py` suite passes before continuing.

### Increment 2 -- factor Poisson height-state construction

Files:

- `orgui/datautils/xrayutils/CTRfilm.py`
- new `orgui/datautils/xrayutils/test/test_CTRfilm.py`, importing the
  deterministic oracle from `test/_poisson_oracle.py`

Work:

1. Extract the profile-support calculations currently embedded in
   `PoissonSurface.createLayers` into one private helper which applies **no**
   retention mask. It returns a frozen candidate record: support bounds,
   layer numbers, cumulative material occupancy, exposed-height probabilities
   (`probability(layer_numbers + 1)`), Film correction occupancy, and the tail
   probability. Selection is then two separate callers of that record:

   - the coherent selector reproduces today's union mask
     `(exposed > tail_probability) | (abs(film_correction) > tail_probability)`
     with no renormalization, so `createLayers` behavior is unchanged; and
   - the incoherent selector applies the `exact_layer_count` policy below to
     the exposed-height states alone and renormalizes.

   Do not give the helper a single "retained mask" output. The two policies
   retain different sets for different reasons, and merging them is what would
   silently make the coherent path adopt the incoherent cutoff or the reverse.
2. Add an immutable flat-height result containing layer numbers, normalized
   probabilities, raw retained mass, excluded lower and upper mass, their sum
   as `excluded_probability`, `iter_states()`, and an explicit diagnostic
   `as_array()` -- the complete attribute set sketched under "Flat-height
   surface evaluation". Returned arrays must be caller-owned or read-only; the
   production path must not require full amplitude materialization.
3. Add
   `PoissonSurface.flat_domain_corrections(h, k, l, *, exact_layer_count=10)`.
   For state `n`, the result is the occupancy-one component correction
   relative to the same sharp Film boundary used by current `F_uc`;
   probability is not folded into the amplitude.
4. Build each deterministic state from the material indicator for that height,
   the corresponding termination slab, and the matching negative Film
   termination. Reuse the current generated cells, strain, translation,
   growth/etching, and termination-cycle conventions.
5. Use cumulative Film-layer amplitudes so all state corrections are produced
   in one pass. Do not construct or deep-copy an `SXRDCrystal` per height. The
   prefix range is taken from the candidate support, not from either retained
   mask: it must cover every structural layer between the sharp Film boundary
   and the highest retained state, including layers whose own exposed
   probability is negligible.
6. Keep the calculation side-effect free with respect to persistent coherent
   domain matrices. Parameter changes must be picked up through the same
   synchronization path as `createLayers`.
7. Accept the wrapper's `exact_layer_count` policy value, default 10. For no
   more than that many candidate states, retain all states. For wider support,
   choose a contiguous, mode-containing interval from exact
   `PoissonProfile.probability` masses and expand it until the profile's
   cumulative tail target is met. The mode is `argmax` of that calculated
   probability array with ties resolved to the lower index, never
   `floor(rate)`: for `alpha < 1` the distribution is a two-point step mixture
   convolved with the Poisson and its maximum need not sit at the Poisson
   mode. Never use measured CTR data to choose the support.
8. Renormalize retained probabilities exactly once, after selection. Report
   raw retained and excluded tail masses, but do not use a dimensionless
   probability by itself as an amplitude tolerance. Do not accidentally fold
   the tail into the boundary state through `surface_occupancy` and then
   renormalize it a second time.

Gate:

- Every returned state correction agrees with the independent deterministic
  crystal oracle from increment 0 after crystal-level scaling is applied.
- Interior bins of `surface_occupancy(n)` equal `probability(n + 1)` to
  roundoff, and the terminal bin deliberately does not. Both halves are
  asserted so a later refactor cannot silently exchange the two expressions.
- The coherent and incoherent selectors retain their own sets from one shared
  candidate record, and `createLayers` output is unchanged by the refactor.
- The probability-weighted coherent state sum reconstructs the current
  `PoissonSurface` amplitude within a Q-dependent tolerance established by an
  extended-support convergence fixture or an amplitude-weighted tail bound.
- Supports below, at, and above `exact_layer_count=10` verify all-state versus
  probability-cutoff selection, normalization, and reported excluded mass.
- Positive/negative growth, offsets, layer cycles, termination fits, and copied
  crystals are covered.

### Increment 3 -- implement the general incoherent wrapper contract

Files:

- new `orgui/datautils/xrayutils/CTRincoherent.py`
- new `orgui/datautils/xrayutils/test/test_CTRincoherent.py`

Work:

1. Implement quantity-neutral `IncoherentModel(LinearFitFunctions, ABC)` with
   one primary `coherent_model` and a validation hook. Add
   `IncoherentF2Model` with the `F2` boundary and abstract `_evaluate_F2` hook.
2. Compose the existing fit API in local-then-coherent order. Call
   `super().__init__()` so the inherited empty `basis`, `basis_0` and
   `parameters` exist even for a wrapper with no local parameters. Note that
   `LinearFitFunctions` already implements `getInitialParameters`,
   `getStartParamAndLimits`, `setFitErrors`, `getFitErrors`, `fitparnames`,
   `priors` and `parameter_list` as **local-only** methods. They must
   therefore all be overridden as composites which take the local block from
   an explicit `LinearFitFunctions.<name>(self, ...)` call and append the
   coherent tail. Inheriting one of them by accident does not raise; it
   returns a silently short vector, which is the failure this contract exists
   to prevent. Only `setParameters` and `validate` are genuinely new:
   `setParameters` exists on `SXRDCrystal`, which is not a
   `LinearFitFunctions` subclass, and has no base implementation to reuse.
   Override the inherited `setFitParameters` to raise `NotImplementedError`
   naming `setParameters` as the composite entry point. Nothing calls it on a
   wrapper -- `SXRDCrystal.setParameters` only calls `setFitParameters` on its
   own component unit cells, and a wrapper is never a component -- so raising
   costs nothing and converts a silent wrong-length write into an immediate
   error.
3. Validate complete vector lengths before mutating either block. Preserve the
   wrapped crystal as the tail so existing optimizer parameter assumptions
   remain valid.
   Treat zero-length local or coherent blocks explicitly: do not call a child
   `getFitErrors()` for a block with no fit parameters, because existing
   implementations may raise when no errors are set; concatenate an empty
   array for that block instead.
4. Reuse `Parameter` and `addFitParameter` for model-local state. Implement
   `parametersToDict`, `parametersFromDict`, and `clearParameters` as composite
   operations over the local state and wrapped crystal, following Film and
   other existing owning-model conventions. Keep private local-only
   serialization helpers for the dictionary config factory; do not add a
   parallel public parameter metadata class.
5. Implement `KinematicIncoherentContext`. It lazily evaluates the coherent
   decomposition once per HKL request and applies component state evaluators
   at every required outer-domain-transformed coordinate. Its production
   iterator yields exactly one complete height-state amplitude at a time.
6. Add `CoherentStateEnsembleModel(IncoherentF2Model)` for validated
   probability/amplitude and coherence-matrix sums. Keep the ensemble base
   general; it must not import or test for `PoissonSurface`.
7. Add `register_incoherent_model`, `available_incoherent_models`, and
   `create_incoherent_model`, plus
   `create_incoherent_model_from_config(crystal, config)`. Use
   `cls.model_type` as the only key source; reject duplicate/unknown keys,
   abstract or wrong base classes, and arbitrary import paths. The config
   factory restores only wrapper-local state around the supplied coherent
   crystal.

Gate:

- Synthetic zero-, one-, and multi-local-parameter wrapper classes satisfy the
  complete fit contract. One conformance test catches every forgotten
  composite override at once: `fitparnames`, `getInitialParameters()`,
  `parameter_list()`, `priors` and `getStartParamAndLimits()[0]` all have the
  same length, and the tail of each equals the wrapped crystal's own entries
  in order.
- `setFitParameters` raises on a wrapper instead of writing the local block.
- Copying and serialization retain local parameters and the wrapped coherent
  tail without identity-based component references.
- A synthetic two-state ensemble verifies coherent, incoherent, and partial
  limits independently of the Poisson implementation.

### Increment 4 -- implement `PoissonHeightDomains`

Files:

- `orgui/datautils/xrayutils/CTRincoherent.py`
- `orgui/datautils/xrayutils/test/test_CTRincoherent.py`

Work:

1. Register `PoissonHeightDomains` under `poisson_height_domains` and give it a
   one-value local basis with `parameterLookup = {"incoherent_fraction": 0}`.
   Initialize both `basis` and `basis_0` from the constructor value and expose
   `incoherent_fraction` as a property whose getter returns
   `float(self.basis[0])` and whose setter validates a finite value in
   `[0, 1]` before writing `self.basis[0]`. A numpy scalar view is not
   available for this; the point is that there is no second stored scalar
   which can drift away from fitted state. Seeding `basis_0` as well is what
   makes `updateFromParameters` restore the configured fraction when it is not
   a fit parameter.
   Add the positive-integer non-fit setting `exact_layer_count`, default 10.
2. Accept the coherent crystal as the first positional argument and a stable
   target surface name. Validate exactly one matching `PoissonSurface`, a
   finite fraction in `[0, 1]`, and the supported stacking topology. For the
   first implementation, the target must be the final component in resolved
   stacking order, must already be bound immediately above its Film, and may
   have no same-level or later component whose position could depend on the
   selected height.
3. Use the context to evaluate the common crystal amplitude once and to turn
   raw flat-height corrections into complete state amplitudes with all area,
   weight, and outer-domain factors applied.
   Calls at outer-domain-transformed HKL must return identical layer-number and
   probability vectors; only their amplitudes may depend on the transformed
   coordinates. Reject inconsistent state metadata.
4. Stream one height at a time while accumulating the fully incoherent value
   `F2_incoherent = sum(p_n * abs(A_n) ** 2)` without requiring an
   `(N_height, N_q)` complex array.
5. Use the coherent result already cached in the context as the exact coherent
   endpoint:

   ```text
   F2_coherent = abs(context.coherent.total) ** 2
   F2 = (1 - kappa) * F2_coherent
        + kappa * F2_incoherent.
   ```

   This is numerically identical to `coherent_model.F2(h, k, l)`, but the
   implementation must use the cached total rather than call `F2` and trigger
   a second bulk and Film evaluation. The independently reconstructed coherent
   amplitude remains a diagnostic against the Q-dependent truncation
   tolerance; it is not allowed to perturb the `kappa=0` endpoint.
6. Return finite nonnegative `float64` `F2`. Reject a request for `F` by not
   defining that method.

Gate:

- `kappa=0` equals the wrapped coherent crystal's `F2`.
- `kappa=1` equals the deterministic-crystal `F2` average.
- Intermediate values are the convex interpolation, including the anti-Bragg
  minimum.
- A one-height distribution is independent of `kappa`.
- Bulk--surface interference and outer-domain transforms match the direct
  oracle.

### Increment 5 -- integrate `F2` with resolution and `CTROptimizer`

Files:

- `orgui/datautils/xrayutils/CTRresolution.py`
- `orgui/datautils/xrayutils/CTRopt.py`
- `orgui/datautils/xrayutils/test/test_CTRresolution.py`
- `orgui/datautils/xrayutils/test/test_CTRopt.py`
- `orgui/datautils/xrayutils/test/test_CTRopt_dwba.py`

Work:

1. Change `CTRresolution.sample_structure_factor` to prefer a callable
   `F2(h, k, l)` and retain `abs(F) ** 2` as the compatibility fallback. Keep
   its public return value as resolution-broadened `|F|`.
2. Add one private optimizer helper for kinematical `F2`. It prefers
   `model.F2(...)` and falls back to `abs(model.F(...)) ** 2` for existing
   crystal-like objects and test doubles which implement only `F`. Route the
   no-resolution prediction through `sqrt(helper(...))`.
3. For sampled resolution, evaluate `F2` at every quadrature point before
   integration. For fast convolution, fill the cached input collection with
   `sqrt(F2)` and keep calling `CTRresolution.fast_convolve`, which squares
   its input, convolves, and takes one square root. That is exactly
   "resolution applied to `F2` before the square root" and it introduces no
   second collection-building path, so `_require_structure_factors` and
   `preserve_measurement_metadata` keep working unchanged. The cost is one
   redundant square-root/square round trip per point, at roundoff level.
4. Continue accepting ordinary coherent crystals as the first optimizer
   argument. Accept an `IncoherentF2Model` through the same argument; do not
   add an `incoherent_model=` keyword or a Poisson-specific branch.
5. Store the owned forward model as `self.model`. Keep `self.xtal` as the
   primary coherent `SXRDCrystal` (`self.model.coherent_model` for a wrapper),
   so existing callbacks, constraints, and public `optimizer.xtal` access do
   not receive a new object type. Use `self.model` for forward evaluation and
   model-parameter methods.
6. Let the wrapper's existing fit methods supply local plus coherent model
   parameters. Verify total values, bounds, names, priors, and error slices in
   the existing optimizer order:

   ```text
   [resolution] [callbacks] [subclass]
   [incoherent local] [coherent crystal]
   ```

7. Add `CTROptimizer.n_parameters`, derived from the prepared full vector, so
   that `n_parameters == len(get_parameters()) == len(fitparnames) ==
   len(bounds[0])`. Leave `startp`, `lower_bounds` and `higher_bounds` with
   their existing model-block scope and length; only their docstrings change,
   to point at `get_parameters()` as the authoritative full vector. Record the
   prepared model name/count/bounds signature and reject later structural
   parameter mutations until `prepareFit()` is called again.
8. Keep `CTROptimizer.priors` model-scoped, matching current behavior. Require
   the wrapper's local-plus-coherent priors to match only its combined model
   block and document that optional optimizer prefixes have no prior API.
9. Validate `output_quantity == "F2"` on the kinematical path. Reject an
   incoherent `F2` wrapper before entering any DWBA code; never reinterpret its
   result as reflectivity.
10. Preserve cache invalidation after local or coherent parameters change.

Gate:

- With an ordinary `SXRDCrystal`, predictions and parameter layout remain
  numerically unchanged.
- Existing third-party/test crystal-like objects providing only `F` continue
  to work in direct and resolution-sampled fits.
- At `prepareFit`, `bounds`, `fitparnames`, `get_parameters()` and
  `n_parameters` describe the same full vector with every optional prefix
  combination. `startp` is explicitly excluded from that equality and instead
  satisfies `len(startp) == len(model.getInitialParameters())`, unchanged for
  an ordinary coherent fit with callbacks or fitted resolution.
- `priors` remains explicitly model-scoped and the wrapper's prior length
  matches its local-plus-coherent model block.
- With a wrapper, no-resolution, sampled-resolution, and fast-convolution
  predictions equal direct `F2` reference calculations.
- Resolution is demonstrably applied before `sqrt`.
- Callback, resolution, angle-correction, wrapper-local, and crystal error
  slices reach the correct owners.
- Registered fit callbacks and crystal displacement constraints receive
  `optimizer.xtal`, not the wrapper, while model parameter methods use
  `optimizer.model`.
- DWBA plus an `F2` wrapper fails during `prepareFit` with an actionable error.

### Increment 6 -- dictionary configuration and user documentation

Files:

- `orgui/datautils/xrayutils/CTRincoherent.py`
- `orgui/datautils/xrayutils/test/test_CTRincoherent.py`
- `doc/source/ctr_structure_factors.rst`
- `CHANGELOG.md`
- an example under `examples/CTR/`

Work:

1. Implement `to_config()` as a plain dictionary containing the registered
   type, non-basis settings (`surface` and `exact_layer_count`), and the local
   subtree of the wrapper's `parametersToDict` payload. The local basis is the
   only saved source for `incoherent_fraction`; do not duplicate it in
   settings.
2. Implement `create_incoherent_model_from_config(crystal, config)`. It wraps
   the supplied coherent crystal and restores the local subtree without
   replaying a second coherent parameter copy. Validate unknown keys, missing
   settings, and unregistered model types. Do not add this dictionary to a GUI
   or optimizer-session file format in the first implementation.
3. Document `F`, `F2`, `r`, `R`, and detector counts as distinct quantities.
4. Add the coherent/partial/incoherent Poisson example from this record and
   show `addFitParameter` plus optimizer usage.
5. Add a `feat(phys)` changelog entry only when the implementation and tests
   land.

Gate:

- Dictionary round-trip preserves wrapper type, target name,
  `exact_layer_count`, fixed fraction, fitted values, limits, errors, and
  priors which the existing `Parameter` serialization supports. This feature
  does not broaden prior serialization.
- Old sessions and `.xtal` files remain untouched and load unchanged because
  this increment does not integrate the new dictionary contract with them.

### Increment 7 -- normalize the commit history

Do this last, after increment 6 and before the branch is integrated or
pushed. The commits written during increments 0--6 carry multi-paragraph
bodies which restate their own diffs. They are too long.

The convention for this repository, on top of the Conventional Commits rule
in the root `AGENTS.md`:

- a subject line;
- optionally a body of **zero to three lines**, separated from the subject by
  one blank line;
- optionally footers, separated from the body by one blank line, with tokens
  using `-` in place of whitespace (`Reviewed-by:`, `Refs:`) and
  `BREAKING CHANGE` as the one permitted exception;
- everything wrapped at 72 columns.

The body says what changed and why, not how. Omit it when the change is
simple: a body which restates the subject is a sign that there should not be
one.

Work:

1. Reword the branch's own commits to that shape. They are unpushed, so this
   is a local rebase and rewrites no published history. Check that first with
   `git log origin/poisson-incoh..HEAD`; if any commit has been published,
   leave it and note the exception here instead.
2. Keep every subject's existing Conventional Commits type and scope. This is
   a length and layout change, not a reclassification.
3. Preserve the `BREAKING CHANGE:` footer on any commit which carries one.

Gate:

- No commit body exceeds three lines and no line exceeds 72 columns.
- `git log origin/poisson-incoh..HEAD` contains the same set of changes as
  before the rewording, verified by comparing the tree at the branch tip.

### Deferred increment -- physical coherence kernel

Do not include this in the first implementation. A later change may add a
sample domain-correlation model and measurement-side mutual-coherence input,
with sample-plane lengths in Angstrom and point-dependent `kappa(Q)` including
footprint and reciprocal-space acceptance. Keep scalar
`incoherent_fraction` as the explicit phenomenological model.

### Deferred increment -- GUI/session persistence

Choose a concrete session owner and schema only when a UI or saved optimizer
workflow needs the feature. That layer may persist the dictionary returned by
`to_config()` and interpret an absent dictionary as the existing coherent
behavior. It must continue storing the coherent structure through the
canonical `.xtal` path rather than embedding a second crystal copy.

### Verification sequence

Run the narrowest test after each increment, then the combined scientific
suite:

```text
pytest orgui/datautils/xrayutils/test/test_CTRcalc.py
pytest orgui/datautils/xrayutils/test/test_CTRfilm.py
pytest orgui/datautils/xrayutils/test/test_CTRincoherent.py
pytest orgui/datautils/xrayutils/test/test_CTRresolution.py
pytest orgui/datautils/xrayutils/test/test_CTRopt.py
pytest orgui/datautils/xrayutils/test/test_CTRopt_dwba.py
ruff check orgui
```

### Definition of done

The first implementation is complete only when:

- coherent behavior remains the default and existing coherent tests pass;
- `F2` has one documented meaning throughout the kinematical path;
- the fully incoherent result matches independent flat-height crystals;
- the optimizer consumes the wrapper through its existing model/parameter
  contract;
- resolution acts on `F2` before conversion to `|F|`;
- unsupported DWBA and stacking cases fail explicitly;
- bulk and common Film amplitudes are evaluated once per HKL batch, verified by
  call-count instrumentation as well as timing; and
- no long-lived cache with incomplete invalidation is introduced.

## Performance acceptance

For `N_q` requested points and `N_h` represented heights, the implementation
must not evaluate the semi-infinite bulk or the common Film `N_h` times.
Expected work is approximately

```text
one common coherent evaluation
+ one evaluation of each distinct generated layer/termination amplitude
+ O(N_h * N_q) vectorized combination and `F2` accumulation.
```

Benchmark the direct and quadrature-resolution paths with positive growth and
negative etching. Record peak memory as well as time. The first implementation
streams exactly one height-amplitude array at a time while accumulating
`sum(p_n * abs(A_n)**2)` and the diagnostic `sum(p_n * A_n)`; it does not
materialize an `(N_h, N_q)` complex array in the production path. Any future
batching change requires profiling evidence and must preserve the public
one-state iterator or add a separate opt-in bulk API.

## Literature basis

- J. Harada, [*Evaluation of the roughness of a crystal surface by X-ray
  scattering. I. Theoretical considerations*](https://doi.org/10.1107/S0108767392003246),
  Acta Cryst. A **48**, 764--771 (1992). The coherent CTR roughness factor is
  the squared modulus of the height-population Fourier sum.
- D. Dale, A. Fleet, Y. Suzuki, and J. D. Brock,
  [*X-ray scattering from real surfaces: Discrete and continuous components
  of roughness*](https://doi.org/10.1103/PhysRevB.74.085419), Phys. Rev. B
  **74**, 085419 (2006). Discrete height distributions enter through the
  characteristic function and show the coherent anti-Bragg cancellation used
  in the validation fixture.
- K. Y. C. Lee *et al.*,
  [*Synchrotron X-ray study of lung surfactant-specific protein SP-B in lipid
  monolayers*](https://doi.org/10.1016/S0006-3495(01)75724-4), Biophys. J.
  **81**, 572--585 (2001). The appendix distinguishes coherent effective-medium
  averaging from area-weighted intensity averaging when domains exceed the
  beam coherence length.
- H.-J. Lee *et al.*,
  [*Characterizing Pattern Structures Using X-Ray
  Reflectivity*](https://www.nist.gov/publications/characterizing-pattern-structures-using-x-ray-reflectivity)
  (2008). Pattern periods spanning the beam coherence scale demonstrate the
  boundary of the coherent effective-medium approximation.
- T. S. Lyford, S. P. Collins, P. F. Fewster, and P. A. Thomas,
  [*X-ray investigation of lateral hetero-structures of inversion domains in
  LiNbO3, KTiOPO4 and KTiOAsO4*](https://doi.org/10.1107/S2053273315001503),
  Acta Cryst. A **71**, 255--267 (2015). Their diffraction model combines
  coherent and incoherent contributions with an empirical coefficient and
  shows that effective coherence can depend on reflection, bandwidth, and
  detector acceptance.
- I. A. Vartanyants and I. K. Robinson,
  [*Origins of decoherence in coherent X-ray diffraction
  experiments*](https://doi.org/10.1016/S0030-4018(03)01558-X), Optics
  Communications **222**, 29--50 (2003). The beam is fundamentally described
  by a mutual-intensity function, which need not reduce to one coherence
  length.

## Settled names

The architectural names are decisions in this record:

- module: `CTRincoherent`;
- wrapper/parameter contract: `IncoherentModel`;
- kinematical quantity contract: `IncoherentF2Model`;
- first implementation: `PoissonHeightDomains`;
- squared-structure-factor method: `F2(h, k, l)` on both `SXRDCrystal` and
  `IncoherentF2Model`;
- optimizer construction: `CTROptimizer(incoherent_model, ctrs)`, using the
  existing first positional model argument;
- registry key: `poisson_height_domains`;
- coherent decomposition: `SXRDCrystal.evaluate_kinematic`;
- flat-state evaluation: `PoissonSurface.flat_domain_corrections`; and
- evaluation context: `KinematicIncoherentContext`.

The amplitude/`F2`/reflectivity boundary, coherent default, complete-state
averaging, DWBA exclusion, parameter-block contract, registry role,
common-prefix reuse, and ownership of coherence metadata are decisions.

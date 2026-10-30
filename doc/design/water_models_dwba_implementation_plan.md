# Water models, rough-surface attachment, and DWBA: implementation plan

> Superseded architecture (2026-10-29): the fitted sharp/tanh paths now use
> onset-only optical interfaces and analytic water perturbations, including
> the infinite Hosemann sum. The numerical-profile fitting and quadrature-cache
> optimizations recorded below are historical. See the current water section
> in `doc/source/dwba.rst` and `_CTRwater_analytic.py`. Explicit erf validation
> retains its continuous numerical solver.


Status: implemented in the working tree; validation and implementation notes below.
Date: 2026-10-28.

This plan combines [issue #83: water on a rough surface](https://github.com/tifuchs/orGUI/issues/83)
and [issue #81: DWBA of water models](https://github.com/tifuchs/orGUI/issues/81),
including the agreed addition of an analytic tanh optical reference and a
semi-infinite water incident medium.

## 1. Scope and decisions

- `WaterModel` owns the intrinsic density profile, form-factor convention, and
  fit parameters. Support erf step, layered, single excess layer, and new tanh
  profiles. Existing oxygen variants remain explicit material approximations.
- `PoissonSurface` may own one water model. It places that model at every
  terrace using the same height distribution as the solid. Keep standalone
  water and its current stacking behavior available.
- Water occupies the upper half-space. The optical profile approaches bulk
  water as z tends to positive infinity; do not append a vacuum cap.
- The core DWBA result is normalized locally at the immersed interface.
  External cell windows, a remote air-water interface, and macroscopic water
  path attenuation are separate experimental propagation problems.
- Use an explicit reference/contrast decomposition. A sharp equimolar water
  reference is the baseline. The tanh profile additionally supports an analytic
  reference where the *complete dielectric profile* permits it.
- Preserve dry-crystal behavior and existing saved-model loading. Diagnose and
  document verified legacy water defects separately from adding new behavior.
- This work includes coherent specular water scattering and its optical effect
  on all rods, plus correlated water placement in the existing kinematic
  incoherent-height model. It does not add diffuse scattering or incoherent DWBA.
- Do not assign crystallinity to water automatically. Solid crystallinity
  continues to scale only the solid perturbation, leaving optical density intact.

## 2. Current implementation and affected files

| File | Current behavior and required extension |
|---|---|
| `orgui/datautils/xrayutils/CTRuc.py` | Defines `WaterModel`, its density, amplitudes, parameters, and text format. Some plain water amplitudes require unavailable molecular data; oxygen variants contain density/amplitude inconsistencies to characterize. |
| `orgui/datautils/xrayutils/CTRfilm.py` | `PoissonSurface` already shares terrace geometry between coherent and flat-height calculations. Add water ownership and expose reusable physical terrace placement. |
| `orgui/datautils/xrayutils/CTRcalc.py` | Combines component amplitudes, density, optical profiles, parameters, and serialization. Include owned water once and preserve stacking semantics. |
| `orgui/datautils/xrayutils/CTRoptics.py` | Samples water on the structural grid; provides slab fields and profile boundaries. Add explicit water asymptotes and continuous-reference field support. |
| `orgui/datautils/xrayutils/CTRdwba.py` | Rejects water, assembles atomic reference shares, appends vacuum, and uses vacuum-based observable normalization. Add continuous contributions and medium-aware geometry. |
| `orgui/datautils/xrayutils/cpp/CTRcalc_cpp.cpp` | Native layered solver already uses the incident index in the conserved tangential wavevector and incident normal wavevector. Preserve its established branch and polarization conventions. |
| `orgui/datautils/xrayutils/cpp/CTRdwba_cpp.hpp` | Atomic DWBA uses slab plane-wave channels and polarization factors. Extend carefully for non-vacuum ambient and continuous-reference fields. |
| `orgui/datautils/xrayutils/CTRincoherent.py` | Flat-height states must include attached water at the same height as their solid termination. |
| `orgui/datautils/xrayutils/CTRopt.py` | Exercise owned-water parameter collection, frozen optical references, preparation, and batch evaluation. |

Read the applicable root, scientific-core, and documentation `AGENTS.md` files
before implementation. Preserve unrelated working-tree changes.

## 3. Mathematical contracts

### 3.1 Coordinates, units, and optical convention

Use z increasing from substrate into water. All internal profile lengths are
Angstrom, q is in inverse Angstrom, and internal angles are radians. Existing
water lengths stored as fractional lattice coordinates are converted once at
the evaluation boundary. Use physical Cartesian heights and reciprocal vectors;
do not assume `q_z = 2*pi*l/a[2]` for arbitrary reference transforms.

Photon energy remains eV at `setEnergy`, consistently with the existing API.
The repository uses `n = 1 - delta - i*beta` and negative-root Renaud normal
wavevectors. Retain those signs throughout fields, phases, and amplitudes.

For incident medium water and its specified glancing angle alpha_w,

```math
k_0 = 2\pi/\lambda,\qquad
k_\parallel = k_0 n_w\cos\alpha_w,\qquad
\kappa_j^2 = k_0^2 n_j^2-k_\parallel^2.
```

The existing field solver uses `kz = -kappa`, with the physical absorption
branch. Normalize incident amplitude at a finite reference plane in homogeneous
water. There is no requirement for a finite-amplitude source at infinity in an
absorbing half-space.

For externally supplied geometry, propagate the tangential wavevector through
the actual entrance interface; a real external angle is not automatically a
real angle in absorbing water. Keep that conversion explicit.

### 3.2 Terrace following

Let H_m denote the physical water attachment height of terrace m and p_m its
probability. Let z_w be the intrinsic water offset relative to that height.

```math
\overline{\rho}_w(z)=\sum_m p_m\rho_w(z-H_m-z_w),\qquad
\sum_m p_m=1.
```

For identical translated profiles, the coherent kinematic amplitude is

```math
F_w(q_z)=F_w^{\rm local}(q_z)e^{iq_z z_w}
\sum_m p_m e^{iq_z H_m}.
```

Apply this relation to the intrinsic profile without a second Gaussian
roughness factor. `sigma_0` describes intrinsic width; terrace roughness comes
from p_m. Surface probabilities are geometrical, not multiplied by termination
atomic occupancy or crystallinity.

Derive H_m from the shared terrace geometry including strain, layer order,
stacking offsets, and the structural top boundary. Do not replace it with an
atom-coordinate maximum or a mean-height approximation. The layer index for
the top filled solid layer and the height-change index differ by one in the
current probability helpers; pin this distinction in tests.

For incoherent kinematic domains, use

```math
I=\sum_m p_m\left|F_{{\rm solid},m}+F_{w,m}\right|^2.
```

Water follows the same state as the solid. It must not be left in the common
fixed-height amplitude or averaged independently before squaring.

### 3.3 Reference and perturbation

Define `rho` as the model's complex scattering density, with molecular or
atomic number-density normalization distinguished from electron density.

```math
\rho_w=\rho_{w,0}+\Delta\rho_w,\qquad
\Delta\rho_w\longrightarrow0\quad(z\to+\infty).
```

The reference tends to the *exact model bulk-water value*. Add that contribution
to the solid optical reference, keeping separate tags for subtraction. Where
a sharp water reference is used, choose its dividing surface by

```math
\int_{-\infty}^{\infty}
[\rho_w(z)-\rho_\infty\Theta(z-z_G)]\,dz=0.
```

This zero-excess condition fixes z_G; the decomposition identity alone does not.
Apply identical terrace shifts and probabilities to actual and reference water.
Subtract analytically or in a stable combined integrand before evaluating
semi-infinite tails; do not subtract separately divergent bulk-water integrals.

For the erf step, z_G = z_0. For a normalized Gaussian G and layered model,

```math
\rho_w(z)=\rho_\infty d\sum_{j\ge0}w_j
G(z-z_0-jd;s_j),\qquad
s_j^2=s_0^2+j s_b^2,\quad w_0=A_0,\quad w_{j>0}=1,
```

the transform, with `exp(+i*q*z)`, is

```math
\widetilde\rho_w(q)=\rho_\infty d\,e^{iqz_0-q^2s_0^2/2}
\left[A_0+\frac{t(q)}{1-t(q)}\right],\qquad
t(q)=e^{iqd-q^2s_b^2/2},
```

and the equimolar position for decaying layering is

```math
z_G=z_0+d(1/2-A_0)-\frac{s_b^2}{2d}.
```

The successive-broadening term is required. Before implementation, this formula
was checked against Gaussian cumulative integrals and the small-q transform at
three parameter sets. Add durable regression tests. Use `expm1` and small-q
limits to avoid cancellation. An indefinitely periodic model with `s_b=0` does
not approach homogeneous water; retain any existing kinematic support, but
reject that non-decaying case for homogeneous-ambient DWBA unless a separate
periodic-ambient model is implemented.

For a step plus an excess Gaussian of integrated density Gamma,

```math
\rho_w=\rho_{\rm step}+\Gamma G(z-z_1;s_1),\qquad
z_G=z_0-\Gamma/\rho_\infty.
```

Map each legacy `A0` to its actual meaning: a dimensionless first-layer weight
and an excess amount per lateral unit cell are different quantities. The simple
shift above assumes excess and bulk have the same complex scattering contrast;
otherwise use a real material dividing surface and retain the unmatched complex
excess in the perturbation.

### 3.4 Tanh profile and exactness

The exactly solvable scalar/s-polarized Epstein transition is defined in
dielectric function, not by independently interpolating n:

```math
\varepsilon_0(z)=\varepsilon_-+
\frac{\varepsilon_+-\varepsilon_-}{2}
\left[1+\tanh\left(\frac{z-z_0}{w}\right)\right],\qquad
\varepsilon=n^2.
```

Use `watertype="tanh"`. Preserve RMS-width semantics by
`w = sqrt(12)*s_0/pi`, where s_0 is the physical width in Angstrom. z_0 is the
zero-excess dividing surface of the symmetric transition.

For lossless propagating endpoints, the scalar intensity oracle is

```math
R=\frac{\sinh^2[\pi w(\kappa_- - \kappa_+)/2]}
        {\sinh^2[\pi w(\kappa_- + \kappa_+)/2]}.
```

Production needs complex r, t, and internal fields, including absorption and
phase. Use stable special-function representations and verify them against a
numerical boundary-value solution. Do not substitute complex wavevectors into
the lossless intensity formula and treat it as the full solution.

Keep the distinction between the exact wave potential `epsilon` and the
first-order X-ray relation between density and n explicit. Define the tanh
reference in epsilon with the established endpoint indices, then derive n and
the matching contrast consistently. Do not silently change all existing optical
density conversions, or construct a second inconsistent tanh in n or density.
Any retained weak-susceptibility approximation must be named and tested.

Exact analytic evaluation applies when the total reference is a supported
Epstein segment with valid boundary matching. A weighted sum of shifted tanh
profiles, or overlapping film/water transitions, is generally not another
Epstein profile. Never fit such a sum to one tanh as an unreported shortcut.
Use a converged continuous/numerical reference in those cases. Multiple analytic
segments may be combined only with explicit endpoint field matching; multiplying
isolated-interface reflectivities is not a multilayer solution.

Exact s-polarized support is required. p polarization must use its proper
electromagnetic equation, numerically unless a separately verified analytic
solution is available. Unpolarized results combine independently computed s and
p intensities. If physical continuum equals reference, the water perturbation
vanishes within the chosen optical/scattering convention.

### 3.5 DWBA matrix elements and observables

For a slab reference, reuse the four channel coefficients and evaluate

```math
F_w^{\rm DWBA}=A_{\rm ref}\sum_{s,\mu,\nu}C_{s,\mu\nu}
\int_{z_s^-}^{z_s^+}\Delta\overline\rho_w(z)
e^{iq_{s,\mu\nu}z}\,dz.
```

Gaussian and erf terms admit analytic finite-interval integrals. Treat complex
arguments, near-zero q, and semi-infinite tails stably. Tail selection must be
controlled by an error estimate rather than a fixed number of extra Gaussians.

Inside a continuous tanh reference, fields are not four constant-wavevector
plane waves. Evaluate the distorted-field product against contrast using the
continuous field evaluator and controlled quadrature or a verified analytic
integral. Atomic contributions inside that region also require those fields at
their physical positions; changing only the Fresnel coefficient is insufficient.
Keep the existing efficient homogeneous/slab channels outside such regions.

The laterally averaged water directly contributes only to the specular coherent
amplitude. Its optical reference affects all rods. The exact specular mask and
reference transforms must agree with the rest of DWBA.

Derive the local-ambient Green-function normalization and polarization factors
before changing `DWBAResult.amplitude_prefactor`. Replace vacuum-only incident
and exit quantities consistently, including native vector factors and geometry.
For s polarization the exit normal wavevector enters the conversion; p and
flux observables require the corresponding field/admittance convention. Do not
implement this by blindly multiplying the old result by n_w or its real part.
Keep externally calibrated intensity separate from locally normalized complex
reflection amplitudes in absorbing water.

## 4. Proposed API and lifecycle

These names are implementation targets; follow existing constructor conventions
when adding the optional keyword:

```python
water = WaterModel(a, alpha, watertype="tanh", name="water")
surface.attach_water(water)
# Equivalent optional PoissonSurface(..., water=water) constructor support.
```

- `surface.water` is absent/None by default. `attach_water(None)` detaches it.
  Require unique ownership; reject registering the same instance both as owned
  water and a top-level crystal component. Treat ambiguous multiple independent
  water half-spaces as unsupported until their mixture semantics are explicit.
- Extend a general owned-component traversal rather than putting water in an
  API that promises only `UnitCell` objects. Apply energy, fit values/errors,
  copying, and persistence exactly once. Namespace child parameters by owner.
- Expose an immutable terrace-state description containing physical attachment
  heights, state identifiers, probabilities, and excluded-tail diagnostics.
  Reuse geometry from `_layer_domain_geometry` and probability candidates;
  do not derive probabilities from rounded atomic occupancies.
- Keep current coherent and incoherent truncation contracts distinct. Water's
  bulk asymptote must remain normalized even when finite terrace support is
  truncated. Prefer analytic characteristic functions/CDFs where possible;
  otherwise expand support to a stated tolerance and account for omitted mass.
  Do not silently renormalize the existing solid coherent amplitudes.
- Add explicit geometry metadata identifying incident-medium angles versus
  external/vacuum geometry. Existing dry calls keep their numerical meaning.
  For water, provide a clear local-angle entry point; HKL convenience preparation
  must perform the corresponding medium-aware inversion or fail clearly when
  the external propagation geometry is missing.
- Preserve immutable prepared references. In automatic mode, changes to water
  shape, asymptote, energy, placement, or terrace distribution rebuild affected
  optics. In explicitly frozen-reference fitting, evaluate changed density
  against the stored reference, subject to compatibility checks. A changed bulk
  asymptote cannot be treated as a localized perturbation of the old ambient.

## 5. Implementation stages and acceptance gates

### Stage 0: characterize existing water behavior

Add `test_CTRwater.py` with small deterministic fixtures. Record working legacy
`F_uc` results, standalone stacking, energy dependence, area normalization,
parameter order, and text/HDF5 round trips. Compare each density transform with
its amplitude, subtracting the bulk asymptote correctly.

Audit oxygen-variant scale factors, Gaussian widths, unexplained density phase
factors, real/reference-space transforms, zero-width limits, and disabled plain
variants. Record verified defects and separate their repairs from compatibility
fixtures. Do not bake an inconsistent formula into the new shared primitives.

Gate: documented per-variant normalization and known discrepancies, with focused
tests sufficient to identify intentional numerical changes.

### Stage 1: shared intrinsic water profiles

Refactor only the water portion of `CTRuc.py`, extracting a private helper module
if needed. Implement shared density, Fourier transform, bulk asymptote, dividing
surface, and interval-integral primitives. Add tanh density/transform support.

For previously unavailable plain profiles, provide an explicitly documented
independent-atom H2O approximation using bundled O and H factors. Keep oxygen
approximations named; never advertise either as the unavailable measured
molecular form factor. Use the same selected convention in density, kinematic
amplitude, optical constants, and DWBA. Include analytic small-q and sharp-step
limits and parameter validation.

Gate: transform/density agreement, equimolar identities, physical units, finite
edge behavior, and legacy compatibility or documented defect corrections.

### Stage 2: surface-owned water and kinematic integration

Implement the shared terrace description and water attachment in `CTRfilm.py`.
Extend `CTRcalc.py` ownership, amplitude/density assembly, fitting, energy, and
persistence. Add water to matching flat states in `CTRincoherent.py` without
duplicating its contribution in the common amplitude.

Gate: explicit independently assembled flat-terrace ensembles agree with
coherent and incoherent calculations for growth and etching, multiple termination
types, strain, offsets, and tilted repeats. Zero roughness recovers flat water;
copied/restored models keep independent ownership and live parameter behavior.

### Stage 3: water ambient and reference assembly

Extend `CTRoptics.py` and `CTRdwba.py` reference construction with exact
semi-infinite water endpoints and tagged water reference contributions. Resolve
profile boundaries geometrically, rather than letting coarse sample midpoints
move an equimolar interface. Preserve negative-z tails near the onset; the
current sampling start at nominal onset is insufficient for full transitions.

Derive and implement the ambient-aware geometry, s/p normalization, polarization
vectors, and observable contracts. Retain existing vacuum results. Provide
local incident-medium geometry and explicit errors for unsupported external
geometry conversion.

Gate: water/solid Fresnel and field-boundary tests for s and p, absorption,
critical-region behavior, no spurious top reflection, and invariance under
extending the homogeneous-water tail. Uniform-medium contrast gives zero
reflection. Verify source/reference-plane shifts by phase and attenuation.

### Stage 4: continuous water contrast with slab references

Introduce water descriptors beside atomic descriptors, with separate packing
paths. Preserve the public total/reference/contrast identities; document how
legacy names such as `F_atomic` represent continuum contributions, or add a
clear compatible alias without duplicating stored arrays.

Implement water interval integrals and complete reference subtraction. Keep
native atomic evaluation for the existing slab path. Support mixed rod arrays,
diagnostics, unpolarized evaluation, prepared references, and optimizer batches.

Gate: direct quadrature oracles, high-angle kinematic limits, exact cancellation
when actual water equals reference, zero direct off-specular water contribution,
and finite values near critical angles. Changing water optics must affect
non-specular solid amplitudes even though direct water scattering is zero there.

### Stage 5: analytic tanh and general continuous fields

Add an immutable continuous-field evaluator alongside `LayeredElectricField`.
Implement exact complex s-polarized Epstein amplitudes and internal fields, plus
a converged numerical s/p solver for validation and general profiles. The
numerical method must retain the continuous material profile rather than replace
it with an uncontrolled single effective tanh.

Add analytic eligibility detection with inspectable solver metadata. Match
eligible segments using fields and derivatives/admittances at explicit finite
planes. Bound tail-matching error. Use the full numerical reference for
overlapping transitions and terrace-averaged tanh profiles.

Extend atomic DWBA evaluation to sampled continuous fields with correct vector
normalization; preserve homogeneous bulk analytic sums. Add contrast integration
against continuous fields. Profile evaluation and field integration tolerances
must be independent of the plotted density grid.

Gate: complex r/t phase and fields agree with an independent numerical solution;
lossless scalar reflectivity matches the sinh formula; sharp-width and
zero-contrast limits recover Fresnel and uniform propagation. Verify flux
conservation where applicable, absorption branches, p fallback, and combined
solid/water DWBA against a converged fine-slab reference. Include an overlapping
transition and a multi-terrace case that cannot use the isolated tanh shortcut.

### Stage 6: persistence, cache, fit, and performance regression

Complete model text and fit-dictionary/HDF5 persistence for owned water, profile
type, and form-factor selection. Old files imply no owned water. Round-trip
plain and oxygen variants, tanh, shared parameter names, fit errors, and legacy
surface constructors. Keep prepared fields and caches runtime-only.

Test automatic invalidation, explicit frozen-reference fitting, changed ambient
rejection, batch reuse, and unchanged solid-crystallinity cache behavior. Measure
representative specular scans through L=3 and non-specular rods; compare accuracy
at tightened tolerances before choosing defaults. Bound live temporary memory
and stream terrace/layer sums when possible.

Gate: no stale optical fields after supported mutations, reproducible saved
models, and documented convergence/performance for representative examples.

### Stage 7: documentation and executed example

Update `doc/source/ctr_structure_factors.rst` and `doc/source/dwba.rst` with the
attachment API, profile normalization, water ambient, local angle convention,
reference choices, exactness limits, and external-cell boundary. Extend the
DWBA design math where the ambient and continuous-field contracts change.

Add a short RuO2/TiO2 water notebook comparing dry, standalone water, and
surface-following water. Show intrinsic/averaged density and 00L reflectivity
through L=3, with the incident medium and angle convention visible. Include
erf/tanh comparison at equal RMS width and one layered-water example. Do not
compare dry and wet curves under an unstated change of experimental geometry.

After implementation, add an Unreleased changelog entry and regenerate
`doc/source/release_notes.rst` through `doc/generate_release_notes.py`. This plan
alone does not claim a released feature or require a changelog entry.

Gate: notebook executes from a clean kernel, contains finite verified outputs,
and the Sphinx build adds no new warnings.

## 6. Validation order

1. New profile and terrace tests, followed by existing water/optical-profile
   tests and relevant CTR calculation, Poisson surface, and incoherent tests.
2. New ambient and DWBA tests, followed by existing `test_CTRdwba.py`,
   `test_CTRopt_dwba.py`, and `test_CTRcrystallinity.py` coverage.
3. Rebuild the native extension whenever C++ changes; verify tests load that
   build. Compare supported Python/native water paths where both exist.
4. Run `ruff check orgui`, relevant broader scientific regression tests, the
   executed notebook, and the Sphinx documentation build.

Use analytic limits and independent quadrature/boundary-value solutions as
oracles, not duplicate implementations of the same formula. Near zeros compare
absolute complex-amplitude error; elsewhere include relative error. Set numerical
defaults from convergence sweeps, including weak broadening, long layering tails,
critical angles, absorption, and cancellation-dominated specular points.

All seven implementation stages plus characterization are required to complete
the agreed scope. A sampled-only tanh curve, water kinematics without DWBA, or an
isolated-interface intensity formula without internal fields is an intermediate
milestone, not completion.

## 7. References

- [Issue #83](https://github.com/tifuchs/orGUI/issues/83): water placement on a
  rough surface.
- [Issue #81](https://github.com/tifuchs/orGUI/issues/81): water DWBA references.
- [Regan et al., Gaussian layering model](https://liquids.seas.harvard.edu/peter/old.pdf/metal_14.pdf).
- [Diamant and Fernandez-Guasti, Epstein profile and sharp-interface limit](https://luz.izt.uam.mx/mfg/arti/11-15/refl-disc-jo2014.pdf).
- [X-Ray Data Booklet, complex optical-index convention](https://xdb.lbl.gov/Section4/Sec_4-1.html).


## 8. Implementation record (2026-10-28)

The shared profiles live in `_CTRwater.py`; `CTRcontinuous.py` supplies exact
Epstein s fields and the independent numerical Maxwell solver. Slab water
contrast is in `_CTRwater_dwba.py`, and full continuous-reference assembly and
atomic/continuum evaluation are in `_CTRwater_continuous.py`.

`PoissonSurface.attach_water`, constructor `water=`, and `terrace_distribution`
provide ownership and physical placement. The existing flat-state termination
correction now includes its matching water, so `CTRincoherent.py` needs no
separate water path. Model text and HDF5 fit dictionaries preserve water type,
number density, fit values and errors. Explicit parameter names retain the
existing shared-name fitting convention rather than introducing implicit
namespaces; independent water parameters should receive distinct names.

Both reference paths retain exact water endpoints and frozen-reference ambient
checks. Generated HKL/CTRScanGeometry scans use the local incident medium;
measured six-circle geometry requires `angles_in_ambient=True` to assert that
external cell optics have already been handled. Optimizer-generated geometry
sets this explicitly, while measured external records are rejected for wet
models. Native p-polarization normal factors now include the ambient index.

The original implementation used exact complex boundary amplitudes and dense
Maxwell fields during DWBA quadrature. The performance revision below replaces
that redundant ODE solve with exact vectorized hypergeometric fields. General overlapping film/terrace
profiles and p polarization use numerical Maxwell integration. Gaussian water
layers and finite atomic blocks stream through bounded temporary arrays.
Quadrature resolves narrow intrinsic features, compares successive Gauss orders
and raises on nonconvergence instead of silently accepting an unresolved peak.

Regression oracles include density-derivative Fourier transforms and equimolar
integrals for every profile, independent flat Film ensembles for growth and
etching, HDF5/text/copy ownership, live and frozen fits, unpolarized sums, exact
Epstein versus Maxwell fields, uniform/sharp limits and lossless flux, and
native fine-slab DWBA for both isolated and terrace-averaged tanh references in
s and p. The latter caught and fixed a reciprocal p contraction sign and a
bulk-tail branch mismatch when water absorbs more strongly than the substrate.

The executed `examples/CTR/RuO2_TiO2_water_DWBA.ipynb` and its Sphinx copy show
water density and local 00L reflectivity through L=3. A 216-point run at 20 keV
with BLAS threads bounded took approximately 0.2 s dry, 0.6 s fixed erf, 9 s
terrace erf, 23 s terrace tanh and 25 s terrace layered on the development
machine. These are illustrative timings, not a performance guarantee. Shape
integration and field tolerances are independent of the plotted density grid.

Deliberate limits remain those of the agreed scope: one water half-space,
identity outer water domains, positive water weight, decaying layering for a
homogeneous ambient, semi-infinite bulk mode for continuous references, and no
incoherent DWBA, diffuse water scattering or implicit external-cell propagation.

Validation: 407 focused scientific tests and 136 subtests passed; `ruff check
orgui` and `git diff --check` passed. The native extension was rebuilt before
testing the changed p-polarization kernel. Existing NumPy deprecation warnings
remain in CTRplotutil. Both notebook copies contain executed outputs.

The updated DWBA LaTeX document compiled successfully with the desktop compiler.


## 9. Tanh performance revision

Profiling identified repeated reconstruction of water factors and parameters in
the numerical RHS as the main cost: roughly 28 of 34 profiled seconds in the
216-point terrace example. The frozen reference now precomputes its constants;
its RHS evaluates only a vectorized tanh sum and a lookup of solid density.
Homogeneous intervals propagate analytically (s and p); numerical integration
remains only where the total reference varies. The monotonic water endpoint
span bounds the residual tail below 1e-15 of bulk density before treating an
interval as homogeneous.

Eligible s-polarized Epstein references no longer run a second ODE for internal
fields. Connection formulas keep hypergeometric arguments <=1/2; vectorized
series evaluate the field and its analytic derivative together. A condition
estimate routes exceptional ill-conditioned entries to arbitrary precision.
The numerical solver remains an independent oracle and the p fallback.

Continuum quadrature batches nodes, reuses coincident incident/exit fields, and
avoids threaded BLAS dispatch for small contractions. A four-entry cache stores
only continuum integrals, keyed by the frozen reference and decomposition,
angles/polarizations, attenuation and live water state. It neither freezes
atomic parameters nor accepts stale water changes, and clear_cache clears it.

The example now defaults to tanh for both fixed and terrace-following water;
erf and layered water are explicit opt-in comparisons. No approximation of a
terrace mixture by one effective tanh was introduced.


Measured after this revision (216 points, default library threading): the
executed notebook reports 3.50 s / 1.15 s for first/reused fixed-tanh scans,
4.07 s / 1.11 s for terrace-tanh scans, and 3.22 ms for 1000 exact single-interface
coefficients. A separate scan comparison with the pre-revision amplitudes had
maximum abs(error)/max(1, abs(F)) of 6.6e-9 for isolated tanh and 2.8e-8 for
terrace tanh, within the integration contract. Timings vary with machine load.

Validation: 412 tests and 136 subtests passed; Ruff and diff whitespace checks
passed. Added tests compare vectorized hypergeometric derivatives to mpmath,
forbid ODE calls on the exact path, verify homogeneous s/p propagation against
numerical integration, and exercise integral reuse, water invalidation and
live atomic changes. Both example and documentation notebooks were re-executed.


## 10. Piecewise analytic fitting reference and erf confirmation

The earlier full numerical tanh mixture is superseded by finite Epstein
segments for s polarization. Cuts include solid discontinuities and midplanes
between adjacent water onsets. Each segment retains the nearest onset's center
and intrinsic width and matches the full physical dielectric at its two ends.
Two independent hypergeometric solutions match E/E' at the finite cuts; no ODE
is used. Constant intervals use exact homogeneous propagation. All overlapping
physical water steps remain in the physical density, and the reference residual
is assigned to the water DWBA contrast. The arithmetic stores epsilon-1 so the
density reconstruction does not lose precision by subtracting near-unit indices.
The homogeneous bulk tail now explicitly evaluates the substrate side of its
matching plane, which may coincide with a solid discontinuity.

An explicit `erf` continuum uses the same forward molecule factor and intrinsic
RMS width as tanh; the legacy `step` keeps atomic form-factor smoothing and its
sharp reference. The full averaged erf reference is solved numerically for
confirmation. P polarization remains numerical for both profiles. Exact here
means an analytic solution of the declared piecewise s reference, or a converged
numerical erf reference, not exact atomistic DWBA.

Analytic internal fields use checked spectral interpolation for repeated
samples, with direct analytic fallback. A compiled series kernel was measured
and rejected because it was slower than the NumPy implementation. Total first
DWBA evaluation is not guaranteed faster for tanh: internal field construction
can outweigh the savings in boundary matching. The notebook reports separate
optical/full/reused times rather than assuming a speed ordering. Comparisons
use a common phase/normalization plane in absorbing water.


Validation of this revision: 427 focused tests and 136 subtests passed, including
independent ODE and fine-slab comparisons, lossless flux, critical-asymptote
limits, no-ODE assertions for Poisson tanh s fields, separated/overlapping
terraces, full erf profile identity and tightened-tolerance convergence. Ruff
passed. Both notebook copies were executed; the LaTeX math document compiled.

The executed 216-point RuO2/TiO2 comparison measured tanh optics/full/reused at
0.33/3.89/0.58 seconds and erf at 0.30/2.00/0.49 seconds on this run. Thus this
implementation does not yet establish a total-time advantage for tanh. Maximum
absolute reflectivity difference was 2.6e-6; maximum relative difference was
4.32% where the erf reflectivity exceeds 1e-10. Tightening the erf field tolerance
from 2e-10 to 2e-12 changed the sampled complex reference reflection amplitude
by at most 1.72e-8. These are example-specific measurements, not accuracy or
performance guarantees for arbitrary profiles.

## 11. Three-profile CTR audit, 2026-10-29

The executed example and its identical documentation copy now include dry,
terrace-following tanh, erf, and layered water. The intrinsic continuum midpoint
is 3 Angstrom above each terrace. For layered water this same parameter sets
the first Gaussian layer center, not its equimolar dividing plane. The Gaussian
spacing/initial width/cumulative broadening are 2.8/1.0/1.4 Angstrom, with
first-layer weight 0.8. Density profiles reach a bulk-water plateau.

All three rods, 00L, 01L, and 11L, use 216 points through L=3, s polarization,
equal local glancing angles, and one common amplitude-normalization plane.
No curve is individually rescaled. Maximum tanh/erf intensity differences are
6.05%, 0.0409%, and 0.0245%, respectively; median absolute differences are
0.138%, 0.000622%, and 0.000167%. Relative comparisons exclude intensities below
1e-8 of the maximum erf intensity on each rod. The largest 00L discrepancy is
near a minimum. The maximum absolute 00L reflectivity difference is 4.69e-6.
Tightening erf's Maxwell tolerance changes sampled complex reference reflection
by at most 1.79e-8. These differences include physical tanh/erf shape differences,
not just solver error. Layered water is a third physical profile.

The audit added erf/layered coverage to frozen-reference fitting, automatic
preparation, optimizer A-B-A scans, HDF5/copy persistence, and high-angle limits.
An independent Simpson integration of layered-water contrast against sampled
native fields verifies the DWBA integral for two broadening lengths, including
a long layering tail. All earlier density/Fourier, equimolar, terrace, native
fine-slab, polarization, analytic/no-ODE and numerical convergence tests remain.

Completion assessment:

- Erf solves the full current averaged reference numerically to controlled
  tolerance. It is not an exact closed-form solution or exact atomistic optics.
- Tanh uses exact finite Epstein segments for s polarization, including field
  matching and the residual density in DWBA. P polarization remains numerical.
- Layered water uses sharp equimolar reference optics plus localized DWBA
  contrast. Full continuous Maxwell optics of the oscillating density is not
  implemented for this profile.
- For films/terraces, adding water to constant solid slabs still creates the
  previously identified ramps in n0. A complementary smooth solid/water
  reference has not been implemented; shifting the onset outward reduces
  overlap without changing that construction.
- Analytic s propagation is implemented, but an overall tanh fitting-speed
  advantage remains unproven. On this notebook execution, full/reused 00L
  times were 8.71/1.66 s tanh, 5.48/1.38 s erf and 20.62/11.82 s layered.
  Machine load affects timings. Non-specular rods were all below 0.6 s after
  preparing 00L; optical caches were shared between those scans.

Thus all three profiles are usable through the stated paths, but the physical
reference construction and the requested speed objective remain open. Passing
regressions does not resolve those design/performance limitations.

Final validation: all CTR test modules passed (902 tests and 193 subtests;
3 skipped). `ruff check orgui` passed. Both notebook copies contain the same
executed outputs. The Sphinx HTML build succeeded with one existing warning
about the unrelated sensor-thickness notebook missing from a toctree.

## 12. Physical-family consolidation

The public family list is now `step` and `layered`. A newly constructed step
has an erf number-density onset plus optional adsorbed Gaussian, with A0
counting excess molecules per lateral cell. It retains the existing one-layer
offset rule for the Gaussian center; A0=0 removes adsorption. Layered A0
continues to mean a relative first-layer weight. Both default to independent
O + 2H factors, including their atomic-size convolution in density.

`optical_reference` is a separate choice: `tanh` (default for step), `erf`
(full numerical reference), or `sharp` (currently required for layered).
Changing this choice does not change physical density, kinematic amplitude,
or form-factor convention. Tanh references include adsorption feature planes
in their finite-interval matching. Numerical references include the full
Gaussian excess and cannot use the monotonic-endpoint shortcut over a peak.
Both include all physical/reference differences in DWBA.

Deprecated aliases remain accepted with DeprecationWarning. The three oxygen
aliases map to canonical families while preserving scaled O, O, or O2- factors
and their original adsorption-position conventions. The old `1layer` alias
retains sharp optics and adsorption. Old `tanh`/`erf` aliases preserve their
forward-only factors and physical onset shapes. An unversioned plain `step`
file retains sharp optics and `adsorption=False`, even if its formerly ignored
A0 is nonzero or supplied as the default 1 for a four-parameter file. Users can
explicitly enable adsorption on those imports. This also preserves legacy
fit scans of a historically inert A0. The A0 fit-name alias is accepted without
changing old generated `rho_0` labels.

Text/WAT files now write canonical families and version-2 JSON `water_settings`
comments. HDF5 stores the same settings, including compatibility conventions.
Do not strip the metadata or expect old readers that ignore it to reproduce
new step semantics. Six old type names remain import aliases, not public
family choices. Density and form-factor compatibility settings are deliberate
metadata; a change of optical solver does not reset them.

Before editing, all eight previous types were sampled for physical density,
three complex kinematic amplitudes and three complex DWBA amplitudes. Importing
their old text after consolidation preserved density/kinematics to rounding
precision and DWBA amplitudes to at most 2.2e-19 absolute in that check. Golden
kinematic values, old four-parameter files, inactive-A0 fit state, canonical
and oxygen-convention HDF5/XPR/WAT round trips, density/Fourier normalization,
and cache invalidation are covered by regression tests. New adsorbed steps
are checked against independent fine-slab DWBA in both polarizations with
flat and Poisson surfaces. Refining the independent slab grid from 1501 to
3001 nodes resolved the adsorption-peak discretization error without relaxing
the comparison tolerance. No-ODE tanh checks and tightened-tolerance erf
checks also cover the new adsorbed steps.

The notebook and generated xtal examples now compare two optical references
for the same physical step with A0=0.8, plus molecular layering. Earlier
sections' tanh-versus-erf metrics refer to their original different physical
profiles and are historical measurements. The separate n0 ramp issue and the
unproven total tanh speed advantage are not resolved by this consolidation.

Validation of the consolidation: 929 CTR tests and 193 subtests passed, with
3 skips. Four additional canonical/oxygen HDF5 and error-text storage cases
passed after that suite was collected; the full 19-test migration module also
passed on the final source. Ruff passed, the LaTeX math compiled, and both
notebook copies contain identical executed outputs. In this same-density
comparison the maximum tanh/erf intensity differences above the stated rod
cutoff are 0.000381% (00L), 0.000441% (01L), and 0.000306% (11L). The maximum
absolute 00L reflectivity difference is 5.03e-8. These are example-specific
results, not universal DWBA accuracy guarantees.

## 13. Profiling the consolidated models

The four-case rough/flat example and both step optical references were
profiled separately on 00L, 01L and 11L. See
[water_dwba_performance.md](water_dwba_performance.md) for paired timings,
retained optimizations, complex-amplitude comparisons and remaining work.
The largest gain is reuse of unchanged layered-water contrast; this does
not remove the cost of fitting the water profile itself. Adsorbed tanh's
first-call cost remains dominated by its many analytic segments.

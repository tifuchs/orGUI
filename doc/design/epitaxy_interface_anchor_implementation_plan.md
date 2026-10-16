# Epitaxy interface anchor implementation plan

Status: Stages 1 and 2 implemented; final Stage 2 validation recorded below.
Assessment date: 2026-10-16.
Assessed revision: `2f7b1b63df79a988f7b80791b503e1fb7b49551a`.
The remote was fetched during assessment; local `master` and `origin/master`
matched at that revision. Existing unrelated working-tree edits must be preserved.

Implement the two discontinuities in separate stages. Stage 1 corrects the
upper-material anchor reconstructed from cyclic layer domains. Stage 2 removes
parameter-dependent jumps in lower-bulk registration while preserving exact
subtraction of the material already supplied by the underlying crystal.
Stage 2 builds on Stage 1, but each stage needs its own reproducer, numerical
comparison, tests, and release-note entry. Passing Stage 1 does not establish
continuity at lower-bulk registration thresholds.

## Stage 1 implementation decisions (historical stage boundary)

The production correction uses a private `_ideal_top_anchor` method in
`CTRfilm.py`. Interior intervals are half-open. The exact final support endpoint
is valid and uses the last interval with fraction one; nonfinite and
out-of-support coordinates are rejected. Upper intervals must be finite and
positive within one repeat. No public parameters or persistence fields changed.

Film W retains its existing nominal average-spacing convention. Stage 1 does
not redefine it as an exact height or invert the unequal-spacing boundary map.
The lower-bulk anchor, subtraction mask, and lower domains remain unchanged.
The remaining Stage 2 threshold issue is explicitly documented in the topic
guide and changelog.

Maintained regressions live in `test_CTRinterface_anchors.py`. They cover
zero-skew layer counts and density/complex-F convergence; cyclic and unequal
interval geometry; affine endpoints and invalid support; coupling, mismatch,
and offset; relabelling and explicit transitions; construction history and
copies; unchanged signed lower domains; upper site accounting under support
extension; and text geometry plus HDF5 companion fit-setting roundtrips.
Historical source-0035 replay remains unverified; synthetic odd/even registry
stacks provide the maintained fixtures.

Validation on 2026-10-16:

- Before the correction, the new zero-skew tests failed for reversed-order
  widths 17, 18, and 32; the corresponding default-order tests passed.
- All 40 new regression cases pass. For nominal widths 17, 18, and 32, the
  corrected reversed-order counts are 11, 12, and 26 at all three skew signs
  (negative, zero, and positive). The one-repeat height loss is gone.
- The combined anchor, CTR calculation, Film, optical-profile, DWBA, and DWBA
  optimizer suites pass: 308 tests and 99 subtests. The run used a freshly built
  CPython 3.14 native CTR extension selected explicitly in the test process;
  no installed binary was replaced. The default local extension discovery
  encountered incompatible pre-existing build artifacts before this selection.
- `ruff check orgui` and `git diff --check` pass. Existing NumPy deprecation
  warnings in `CTRplotutil.py` remain outside this change.
- Release notes were regenerated. A scoped Sphinx build of the changed CTR
  guide and release notes succeeds, with six warnings from omitted surrounding
  documents/toctrees and one existing release-note reference warning. The full
  site build remains unverified because importing the unrelated GUI API crashes
  in `pyFAI.units.ne_equation` with Windows exception `0xc06d007f`.

Detailed local validation logs are in `tmp/epitaxy_stage1_native_tests.txt` and
`tmp/epitaxy_stage1_scoped_sphinx.txt`. The Stage 2 lower-registration threshold
was unresolved at that stage and was not included in its completion claim.

## Stage 2 implementation decisions and evidence

The lower cycle origin is derived from the underlying stack, independently of
skew. With absolute origin `R`, first support index `j0`, lower repeat `a_b`, and
`N` layers, the lower support starts at `R + j0*a_b/N`. Equivalently the local
anchor subtracted from generated domains is `-j0*a_b/N - (R - below_H)`.
The coupled field starts at this same support boundary. Its accumulated
displacement propagates to the Film-side nominal origin as before.

Direct bulk includes native repeat zero in `F_bulk`. Its reference is the next
occurrence of the requested lower-cycle layer above that native termination;
assuming all bulk cells already end at zero would leave a duplicated native
cell in synthetic and user-defined positive-origin cells. Above a finite Film,
the actual upper boundary and next cyclic layer define the reference. An explicit
transition to another starting layer adds its forward phase offset, preserving
the lower lattice under nonuniform intervals. Layer renaming does not change
native geometry. Standalone interfaces without crystal context use a sharp
lower lattice at negative structural indices relative to the supplied height.

`SXRDCrystal.apply_stacking` binds preceding components and their weights before
generating each interface. This derived context is rebuilt on stacking and is
not added to persistence formats. Subtraction copies actual source atom rows,
domain matrices, and occupancies inside the half-open ideal lower support
window, including both bulk and finite registry Film when they overlap it.
It no longer uses `P_top >= 0.5`. Equal-material source corrections are merged
into their lower addition record; other source corrections have separate
records included in amplitudes, density, optical profiles, and DWBA. Energy and
reference-frame changes invalidate these derived copies.

The support window covers the union of statistical support and the fixed sharp
boundary. Large positive/negative skew otherwise leaves a fully occupied gap
or an unremoved part of bulk. `fixed_ucs` is rounded outward to complete cycles
and expanded to include the sharp boundary when needed. Its stored parameter
value is not rewritten. Film W retains the Stage 1 nominal average-spacing
meaning, including the existing integer count quantization.

The finite replacement uses the existing unattenuated surface-component
convention; the infinite bulk attenuation convention is unchanged. The
independent amplitude oracle compares the finite signed correction, so the
common deep bulk cancels without comparing an arbitrary slab to an attenuated
infinite sum. Tail convergence assumes the lower template and the actual deep
lower stack describe the same material and registry. In particular, a Film
whose native first layer is at 0.1 but is placed at absolute zero must be
paired with a correspondingly registered bulk in a tail-convergence fixture.
Otherwise varying support changes the physical amount of replaced material.

Validation on 2026-10-16:

- Six initial regressions failed before editing the implementation: both
  `S = +/-2/3` thresholds with coupling 0, 0.4, and 1.
- In the three-layer width-0.5, Film-W-32.25 reproducer, the old positive
  threshold produced a 6.58069633 Angstrom height jump and about 77.63 electrons
  of complex-amplitude change. At `2/3 +/-1e-8`, the corrected Film has 16 layers
  on both sides and no physical-height jump. Maximum complex-F and density
  differences are 3.31e-7 electrons and 3.04e-8 electrons/Angstrom^3. They decrease
  proportionally to epsilon from 1e-4 through 1e-10; nominal-origin displacement
  is continuous too (9.87e-8 Angstrom across the 1e-8 pair).
- Explicit signed atom coordinates and occupancies match a separately
  enumerated finite bulk/registry replacement. Complex amplitudes agree within
  1e-9 electrons across the sampled rods. Cases include unequal intervals,
  nonzero native origins, explicit transitions, coupling 0/0.4/1, mismatch,
  offset, and two registry Film thicknesses.
- All reachable nearest-repeat ties for width 1 in two-, three-, and four-layer
  cycles are covered in every rotation, including exact ties, reverse updates,
  and copies. Additional tests cover CDF=0.5, support changes, large skew,
  differing cycle labels, energy refresh, and text/HDF5 persistence at a
  lower-reference threshold. DWBA checks both merged and separate source
  correction records.
- The final combined native run passes 374 tests and 99 subtests, including
  63 Stage 2 registration cases and the Stage 1 suite. Existing NumPy
  deprecation warnings remain unchanged. `ruff check orgui` and
  `git diff --check` pass. The scoped Sphinx build succeeds with the same
  seven scoped/existing warnings as Stage 1. Logs are in
  `tmp/epitaxy_stage2_native_tests.txt` and
  `tmp/epitaxy_stage2_scoped_sphinx.txt`.

Historical source-0035 replay remains unavailable. The full Sphinx build still
has the unrelated GUI/pyFAI import crash described under Stage 1; the changed
guide and generated release notes are validated with the scoped build.

## Evidence and limits

The proposal's original `ruo2_model` loader and source-0035 historical models
were not replayed. Assessment used synthetic two-layer Si/C cells with matched
normal repeats of 6.58069633 Angstrom, Skellam width 0.15, coupling 1, offset 0,
and a sharp Poisson surface. Reversing the bulk registry with
`translate_layered([0, 0, -1])` produces interface order `[2, 1]` through native
stacking. Only interface `S` is registered as a fit parameter; Film `basis` and
`basis_0` both hold the requested width.

| Nominal Film width in layers | Current generated layers at negative/positive skew | Candidate generated layers at negative/zero/positive skew |
| --- | --- | --- |
| 17 | 12 / 10 | 11 / 11 / 11 |
| 18 | 13 / 11 | 12 / 12 / 12 |
| 32 | 27 / 25 | 26 / 26 / 26 |

For skew +/-1e-8, the original reversed stack loses 6.58069633 Angstrom of
surface height. An in-memory application of the proposed upper-anchor formula
removes that jump. Across three sampled rods, the maximum complex-amplitude
sign-pair difference decreases from approximately 1.41e-4 electrons at skew
+/-1e-4 to 1.41e-8 electrons at +/-1e-8, then approaches floating-point error.
The sampled default-order prediction is unchanged.

The prototype passed 156 comparisons against cumulative boundary interpolation,
including rotated and relabelled two- and three-layer cycles, unequal spacings,
and nonzero native origins. Twelve coupling/lattice/offset zero-crossing
combinations and odd/even registry construction-history checks passed. With the
candidate active, `test_CTRcalc.py` and `test_CTRfilm.py` reported 151 passed,
5 skipped, and 99 passing subtests. These are feasibility results, not complete
material-accounting or persistence validation.

A distinct failure remains with a uniform three-layer cycle, the same repeat,
interface width 0.5, Film width 32.25, coupling 1, and offset 0:

| Quantity | S = 2/3 - 1e-8 | S = 2/3 + 1e-8 |
| --- | --- | --- |
| Structural support index range | -9 through 17 | -9 through 17 |
| Generated Film layers | 16 | 16 |
| Lower registration anchor in Angstrom | 19.74208899 | 26.32278532 |
| Film top in Angstrom | 74.5812250733 | 68.0005287433 |

This one-repeat jump survives Stage 1. Identical support indices and Film layer
counts isolate it from support truncation and Film rounding. The maximum
sampled complex-amplitude difference is about 77.63 electrons.

Local assessment artifacts are `tmp/epitaxy_anchor_assessment.py` and
`tmp/epitaxy_anchor_assessment_results.txt`. These are scratch evidence tied to
the assessed pre-fix revision; the script's in-memory source substitution does
not target the corrected implementation. The maintained tests now contain the
required synthetic fixtures and do not import from `tmp`.

## Coordinate contract shared by both stages

| Quantity | Units and meaning |
| --- | --- |
| `loc = mu1 - mu2` | Statistical mean in structural-layer intervals. It is not an absolute height. |
| `unitcells` | Integer structural support indices despite the variable name. |
| `x = loc - unitcells[0]` | Continuous structural interval coordinate from generated support start. |
| `uc.layerpos[layer_id]` | Native fractional normal-repeat coordinate of a layer boundary. |
| Ordered layer spacing | Fractional normal-repeat distance to the next boundary, including cycle wrap. |
| Domain matrix | Affine transform of native fractional coordinates into generated fractional coordinates. |
| `a3_top`, `a3_bottom` | Normal repeat lengths in Angstrom under the existing implementation's geometry convention. |
| `loc_absolute` | Current statistical/stack reference in Angstrom, assigned from the component below. |
| `stacking_loc_absolute` | Film-side nominal width origin in Angstrom, including accumulated displacement and offset. |
| `stacking_height_absolute` | Actual upper support boundary in Angstrom; Film starts here. |
| Film `W` | Requested nominal width in structural layers, not the remaining generated layer count. |

Preserve supported cyclic rotations and explicit `LayerTransition` mappings.
Layer IDs are labels, not numerical coordinates. Relabel coordinates, cycle
order, and transition maps together when testing invariance. A physical change
of termination may change scattering and must not be confused with relabelling.

Preserve the current normal-coordinate convention; do not broaden this work
into a lattice-frame or arbitrary-rotation redesign. Keep the statistical
boundary, Film width origin, and physical support boundary distinct.

## Stage 1 Correct the upper anchor for cyclic layer order

### 1.1 Establish maintained failing regressions

Add focused tests under `orgui/datautils/xrayutils/test/`, preferably a dedicated
`test_CTRinterface_anchors.py` for the new matrix of cases. Reuse existing
fixtures where practical without moving unrelated tests.

1. Reproduce the table above through `SXRDCrystal.setParameters` and
   `apply_stacking`, using both default and reversed orders and widths
   17, 18, and 32. Include zero and sign pairs at 1e-4, 1e-6, 1e-8, 1e-10,
   1e-12, and 1e-14. Do not use numerical collapse at 1e-16 as success.
2. Test the anchor independently of Film rounding. Build monotonically
   unwrapped boundary positions by accumulating physical layer intervals, then
   interpolate that boundary list. Compare against the generated anchor just
   below, at, and above structural-layer and cycle boundaries.
3. Cover all rotations of two- and three-layer cycles, unequal spacing,
   nonzero native origins, arbitrary numeric IDs, and explicit transitions.
4. Exercise native updates in ascending and descending order, fresh models,
   and deep copies. Verify equality at identical parameters and preservation
   of source coordinates, layer metadata, and fixed reference widths.

Record the expected failures before editing the production formula. Keep the
Stage 2 counterexample separately identified; it must not be mistaken for a
Stage 1 acceptance failure or silently converted into a permanent expected result.

### 1.2 Apply the localized geometry correction

In `EpitaxyInterface.createInterfaceCells` in `CTRfilm.py`, reuse the existing
`top_layer_space` computed from the actual unwrapped stacking order:

```python
structural_index = int(np.floor(loc_rescaled))
cycle_index, order_index = divmod(structural_index, n_layers)
fraction = loc_rescaled - structural_index
layer_id = self.layer_order[order_index]
matrix = ideal_top[order_index][cycle_index]
native_z = self.uc_top.layerpos[layer_id]
interval = top_layer_space[order_index]
ideal_top_loc = a3_top * (
    matrix[2, 3] + matrix[2, 2] * (native_z + fraction * interval)
)
```

Use the half-open interval convention: an exact boundary selects the interval
starting there, with fraction zero. Validate that the coordinate is finite,
the selected domain exists, and the ordered intervals are positive and span a
valid cycle. Negative Python indices must not silently select the wrong domain.
If the coordinate reaches the terminal support endpoint, either evaluate the
last interval's upper endpoint explicitly or raise a documented support error;
choose based on which endpoint states the support API legitimately permits.
Do not extrapolate or clamp an invalid support coordinate.

Replace the incorrect comment claiming every `top_layer_offset` is zero.
A small private helper is appropriate if it makes independent anchor testing
clearer; avoid a public API or general domain-generation refactor.

Leave lower-bulk registration and its subtraction domains unchanged in this
stage. Trace the resulting upper anchor through `ideal_upper`, accumulated
strain displacement, `stacking_loc_absolute`, and Film support subtraction.
Do not compensate by changing Film width, skew sign, or stacking order.

### 1.3 Validate ownership and compatibility

For coupling 0, an intermediate value, and 1, check matched and mismatched
normal repeats and zero/nonzero offset. Compare occupied atom positions,
occupancies, density, and complex amplitudes for interface plus remaining Film.
Assert that their adjoining boundaries meet and that changing support length
does not duplicate fully occupied upper layers.

Include odd/even registry Films and both PoissonSurface and explicit UnitCell
terminations. Preserve default-order predictions where the old anchor was
consistent; unequal-spacing and nonzero-origin cases may intentionally change.
Record those changes against the independent geometry oracle rather than
requiring agreement with an incorrect old result.

Review Film's current conversion
`W - n_layers * (below_H - below_loc) / a3` separately from the anchor fix.
It uses average layer spacing. State whether W retains that nominal convention
for unequal-spacing cycles or requires inversion of an ordered boundary map.
Do not silently redefine W in Stage 1. If an additional width correction is
necessary, specify and test it explicitly before claiming unequal-spacing
end-to-end acceptance.

### Stage 1 completion criteria

- The two-layer zero-skew loss is removed at all three requested widths;
  exact zero agrees with the limiting geometry.
- Independent anchor interpolation passes all valid cyclic boundary cases.
- Complex amplitudes and density converge for shrinking sign pairs; finite
  opposite skews are not required to produce identical predictions.
- Fixed-bulk subtraction coordinates remain unchanged by this patch.
- Maintained regressions, existing CTR tests, and applicable persistence checks
  pass; prediction changes and the remaining Stage 2 issue are documented.

## Stage 2 Keep lower-bulk registration fixed across statistical thresholds

### 2.1 Isolate all affected registration boundaries

Start from the Stage 1 implementation and add the three-layer counterexample
above as a failing test. Use Film width 32.25 to stay away from a Film rounding
tie. Include exact S = 2/3, reverse update direction, and fresh construction.

Instrument or expose test diagnostics for `loc`, support bounds, the selected
integer lower anchor, all three absolute boundaries, Film count, and signed
lower-material domains. Identify crossings where
`loc_rescaled / n_layers = integer + 1/2`, including negative means and Python's
ties-to-even behavior at exact values. Choose widths that make each crossing
reachable within the supported skew range; do not test invalid skew values.

Separate fixed-support sweeps from sweeps that change truncated support. Extend
the matrix to rotated cycles, odd/even layer counts, transitions, unequal
spacing, and nonzero native origins. Record which jumps arise from lower-anchor
selection, sharp-subtraction membership, and actual support changes.

### 2.2 Establish the lower-material reference and subtraction rule

The current formula
`round(loc_rescaled / n_layers) * a3_bottom` lets the statistical mean select a
different lower-lattice registration. Replacing it with the continuous upper
anchor would misregister subtraction. Changing `round` to `floor`, introducing
an epsilon, or freezing a cached first anchor would merely move the jump or make
results depend on construction history.

Derive registration from the actual underlying component, its fixed termination,
cyclic state, and explicit transition. Keep that physical reference independent
of skew. Treat movement of the generated support window as an indexing change:
if its first index changes, compensate the local support origin so every
retained unstrained lower-lattice reference site stays fixed. Coupled material
positions may still move continuously as the occupancy-dependent strain changes.

Before choosing replacement arithmetic, write down and test these rules:

1. Which lower atom sites are already supplied by the underlying stack, and
   which generated negative domains remove those exact sites? Verify both a
   direct bulk interface and one following a finite registry Film.
2. How does a support index map to the lower cell's native layer and integer
   repeat under rotated order and explicit transitions?
3. Where is the coupled displacement field anchored, and how does its absolute
   reference remain stable when support is extended or shortened?
4. How do independent upper and lower lattices, explicit offset, and the
   Film-side nominal origin relate to that fixed reference?

Audit `sharp_top = probability_top >= 0.5` and the associated `sharp_bottom`
subtraction mask. A statistical occupancy threshold is not sufficient proof
that an already-counted bulk atom exists at a generated site. Define the mask
from the sharp underlying material's actual occupied sites if the audit shows
the current classification changes that reference. Test mask and anchor changes
together; a continuous reported height alone does not prove correct cancellation.

Use a finite, explicitly enumerated lower stack as an independent oracle.
Enumerate signed atom positions and occupancies, replace only the relevant
interface window, and compare with the assembled component representation.
For amplitudes, either compare differences where common deep bulk cancels or
use the same finite/attenuated bulk convention on both sides. Do not compare an
arbitrary finite slab directly with an infinite undamped bulk sum.

The replacement formula is a Stage 2 derivation deliverable, not a settled
design decision in this plan. If the desired reference cannot be determined
from existing component metadata, specify the smallest necessary internal state
addition and its persistence behavior before implementation. Preserve existing
public parameter meanings wherever possible.

### 2.3 Implement and test the complete registration change

Keep changes centered on `CTRfilm.py`. Modify `CTRstacking.py`, `CTRcalc.py`, or
`CTRuc.py` only if the reference derivation demonstrates missing information.
Do not redesign the whole stacking API as part of this repair.

Validate the independent-lattice endpoint at coupling 0, the coupled endpoint
at 1, and interpolation at intermediate coupling. Check offset accommodation
at both boundaries. Recalculate physical support, accumulated displacement,
and Film-side origin from the same fixed-reference convention on every update.

Compare signed domain sites and occupancies before comparing density and F.
Require exact site registration within numerical tolerance, cancellation of
only the underlying material that is replaced, and no missing or duplicate
occupied layers at the interface/Film join. Include all Stage 1 cases so this
repair cannot reintroduce the cyclic upper-anchor defect.

### Stage 2 completion criteria

- The three-layer S = 2/3 counterexample is continuous, including the exact
  threshold, without changing the requested Film width.
- Other reachable lower-registration thresholds show no finite jump caused
  solely by changing the integer reference or subtraction classification.
- Negative domains remain registered to actual underlying occupied sites;
  the independent signed-material oracle agrees with the assembled stack.
- Support-window changes alter predictions only within a justified truncation
  budget, not by a full repeat or a fully occupied layer.
- Results are independent of update direction, construction history, and copy
  history, with persistence and Stage 1 regressions passing.

## Numerical tolerances and validation sequence

At identical parameters, begin with absolute tolerances of 1e-9 Angstrom for
coordinates and 1e-9 electrons for complex amplitudes, then justify adjustments
from lattice size, atom count, and summation error. Specify density tolerances
separately in electrons per Angstrom cubed. Check real and imaginary amplitudes
or the norm of their difference, not only intensity or amplitude magnitude.

For different finite skews, measure convergence over several epsilon values
above the floating-point floor. Do not require strict monotonic decrease once
roundoff dominates. Isolate legitimate discrete Film count changes from
reference-induced jumps by using fixtures away from count-rounding ties, then
test the intended quantization behavior separately.

Vary tail probability and support independently. Relate omitted occupancy to
an amplitude bound using omitted atom contributions on the sampled rods, and
test density on a fixed absolute-z grid covering both supports. Do not assume
tail probability itself is an absolute tolerance for height or electron count.
For a coupled field, bound or measure any additional error accumulated in the
displacement through truncated tails.

Run checks in this order for each stage:

```powershell
python -m pytest orgui/datautils/xrayutils/test/test_CTRinterface_anchors.py -q
python -m pytest orgui/datautils/xrayutils/test/test_CTRcalc.py orgui/datautils/xrayutils/test/test_CTRfilm.py -q
python -m pytest orgui/datautils/xrayutils/test/test_CTRoptical_profile.py orgui/datautils/xrayutils/test/test_CTRdwba.py orgui/datautils/xrayutils/test/test_CTRopt_dwba.py -q
ruff check orgui/datautils/xrayutils
```

Use the configured scientific Python environment; on the assessed Windows host,
Miniforge's `Library/bin` must be on PATH for NumPy linear algebra. Fix lint
findings introduced by this work without unrelated scientific cleanup. Expand
testing further only when changed dependencies or failures warrant it.

For text and existing HDF5 model/fit persistence, verify layer state, transitions,
`basis_0`, fixed widths, and regenerated predictions after roundtrip. Compare
fresh and restored objects at the same parameters. Historical predictions need
not match a corrected calculation, but saved parameter values must not be
silently rewritten. Replay source-0035 cases if distributable; otherwise retain
equivalent synthetic registry fixtures and state that historical replay remains
unverified.

## Delivery and compatibility

Keep two separately reviewable changes, with suggested commit subjects:

- `fix(phys): reconstruct epitaxy upper anchors in cyclic layer intervals`
- `fix(phys): preserve lower bulk registration across interface skew changes`

Each change should contain its failing-then-passing regression, localized
implementation, before/after numerical evidence, public docstrings, and the
relevant additions to `doc/source/ctr_structure_factors.rst` and `CHANGELOG.md`.
Regenerate `doc/source/release_notes.rst` with
`python doc/generate_release_notes.py`; never edit it manually. Build Sphinx
when the topic/reference documentation is changed.

Document supported cyclic/transition semantics, boundary origins, and Film W
semantics. State that affected geometries, calculated amplitudes, and fitted
parameters can change. Existing saved fits remain historical records; users
should explicitly recompute affected predictions and reassess fits. If final
changes break generated-data or post-processing compatibility, use `!` and a
`BREAKING CHANGE:` footer as required by the repository convention.

Stage 1 may be delivered independently once its criteria pass, with the remaining
lower-registration limitation stated. Do not label the whole proposal complete
until Stage 2 and the shared accounting/persistence checks pass. Record any
remaining unequal-spacing width decision or unavailable historical fixture
explicitly rather than weakening the acceptance tests.

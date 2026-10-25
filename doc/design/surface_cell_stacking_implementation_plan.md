# Surface cell stacking implementation plan

Date: 2026-10-24.

Status: stages 0–6 implemented and validated in the working tree on 2026-10-24.
Stage 7 remains a separate optional feature. The staged design below records the
implementation contract; the completion notes record the resulting behavior.

## Implementation and validation record

- Miller input now defaults to `cell_mode="surface"`; an explicit matrix still
  defaults to `"lattice"`. Both modes remain explicitly selectable.
- `SurfaceCellSpec.bulk_transform` retains the periodic enumeration basis.
  Projected positions and Wyckoff derivatives retain their image representatives;
  source occupancy and isotropic displacement values survive conversion.
- `UnitCell.stacking_translation` is an immutable fractional tuple and
  `repeat_cartesian` supplies its physical vector. Native text files, copies,
  layer splitting and supercells preserve the state. Nonvertical-cell atom,
  lattice and domain records use full floating-point serialization precision.
- Bulk scattering reuses each backend's accelerated motif calculation and
  evaluates the vector geometric series in the shared wrapper. DWBA already
  accepted a Cartesian repeat vector; its packing now supplies the stored
  repeat. No native ABI or new dependency was necessary. The native extension
  was rebuilt locally before validation.
- Films, interfaces, Poisson height-state corrections, layered translations,
  real-space copies and Fourier density include lateral offsets. Stack state
  carries Cartesian registry separately from normal placement. Normal-only
  strain and rigid gaps retain their distinct meanings. A component above a
  Poisson surface uses nominal mean registry, matching the existing mean-height
  layout convention, rather than a separate upper stack per realization.
- Geometry refresh and prepared-DWBA fingerprints include repeat changes.
  ROD/SUR export rejects a nonvertical repeat before opening the destination.
  Kinematic bulk coherent domains that change a nonvertical repeat are rejected;
  differing bulk repeats require separate cells. Legacy vertical behavior stays
  on its existing path.
- Existing rutile Python/notebook examples were updated and executed. The new
  `examples/CTR/surface_stacking_example.py` was executed against the user's
  published TbScO3 Pnma CIF without modifying the external research directory.
  All four cuts produced 90-degree angles and `(0, 0.4759479411538949, 1)`
  repeat translations. Complex amplitudes agreed with lattice mode at identical
  parent-frame scattering points.
- Validation: 471 tests and 102 subtests passed across the ten relevant CTR
  suites, including 32 new stacking tests. The stacking suite was rerun after
  the final domain-serialization precision change: all 32 passed. These cover
  reference-index conversion across NumPy/Numba/C++, finite Cartesian sums,
  projection, Wyckoff fitting after supercell/save, partial-film continuation,
  normal gaps/strain, interface cancellation, positive/negative Poisson states,
  Fourier density, the DWBA vacuum limit and stale-cache rejection.
- `ruff check orgui` and checks of both Python examples passed. Sphinx built
  successfully with three warnings outside the changed topic pages: geometry
  emphasis markup, an unlisted parallax notebook, and the existing release-note
  anchor. Release notes were regenerated from the changelog.
- Optional stage 7 (physical miscut terraces, width statistics and coherence)
  is not implemented. Crystal repeat shear alone does not imply those effects.

## Outcome

Build surface coordinate cells with `a` and `b` in the surface plane and `c`
perpendicular to that plane, while preserving the original crystal through an
explicit repeat translation. Use the same repeat for real-space atom placement
and reciprocal-space phases. Preserve the existing oblique lattice-cell mode,
full-matrix input, CIF occupancies, displacement import, and Wyckoff fitting.

Implement this as several independently reviewable changes. Keep the new mode
opt-in until films, interfaces, roughness, density, and DWBA support it. During
development, unsupported consumers must reject a nonvertical repeat rather
than silently calculate vertical stacking. Add no dependencies.

## Scientific contract

Use column vectors throughout this contract. Let the existing periodic basis
be `P = [A, B, C_bulk]`, expressed in a common Cartesian frame in Angstrom.

```text
n = cross(A, B) / norm(cross(A, B))
d = dot(C_bulk, n) > 0
C_surface = d * n
S = [A, B, C_surface]
t = solve(S, C_bulk) = (sx, sy, 1)
T = S @ t
```

`S` is the surface coordinate basis; `T` is the physical upward crystal
repeat. The two in-plane vectors need not be orthogonal: only alpha and beta
are constrained to 90 degrees. The normal repeat height is `d`, not the length
of `C_bulk`.

Convert motif positions without strain or changes of origin:

```text
f_surface = solve(S, P @ f_native)
          = (x_native + sx*z_native,
             y_native + sy*z_native,
             z_native)
```

Preserve explicit origins in the same Cartesian frame when converting an
actual terminated cell. A rigid rotation into orGUI's canonical lattice frame
must rotate positions, reference-frame relations, and Cartesian displacement
tensors consistently.

For Cartesian `Q` in inverse Angstrom and `T` in Angstrom, the repeat phase is
`dot(Q, T)` in radians. With reciprocal coordinates in the surface basis it is
`2*pi*(h*sx + k*sy + l)`. `refHKLTransform` converts external indices into the
cell's reciprocal coordinates first; the repeat phase is applied once after
that conversion.

For a downward bulk series including repeat zero:

```text
r_n = r_0 - n*T
F_bulk = F_motif / (1 - exp(-eta - 1j*dot(Q, T)))
```

`eta` is dimensionless attenuation per normal repeat, with the existing public
reference-cell attenuation convention converted explicitly. Retain the
existing termination and repeat-zero phase origin. A finite stack uses a
finite sum or explicit layer copies when their motifs differ.

Do not wrap new fractional `x` or `y` independently during conversion merely
to make them lie in `[0, 1)`: moving individual atoms by lateral cell vectors
changes a finite motif amplitude at fractional in-plane indices. Preserve
unwrapped representatives or track their image offsets. Likewise, moving a
whole motif across a `z` repeat boundary requires the complete repeat vector.

## Proposed API and compatibility

```python
# Final default for Miller-based construction, after all release gates pass.
options = CTRsymmetry.surface_terminations(
    cif_path, miller=(1, 0, 1), cell_mode="surface"
)
cell = options[0].build_unitcell(
    "TbScO3", displacement="isotropic_equivalent"
)
cell.stacking_translation  # (sx, sy, 1.0), fractional in this cell

# Preserve the original periodic lattice-cell representation.
options = CTRsymmetry.surface_terminations(
    cif_path, miller=(1, 0, 1), cell_mode="lattice"
)

# An omitted mode with a full matrix preserves that matrix's geometry.
options = CTRsymmetry.surface_terminations(cif_path, transform=T_bulk)

# Explicitly request a perpendicular coordinate cell from a bulk transform.
options = CTRsymmetry.surface_terminations(
    cif_path, transform=T_bulk, cell_mode="surface"
)
```

- Internally resolve omitted `cell_mode` according to the input form. During
  the opt-in stages, omission retains current lattice behavior for all forms.
- Store a validated immutable three-component `stacking_translation`; initially
  require finite `sx`, `sy` and `sz == 1` within numerical tolerance. Existing
  cells and files default to `(0, 0, 1)`. Do not add an independently fitted
  stacking parameter in this work.
- Retain the original integer bulk transform separately from the projected
  surface coordinate transform, which need not be integer. Specify their
  meanings in metadata and docstrings; never overwrite the original transform
  and then use the projected one as a periodic supercell transform.
- Enumerate terminations using the original periodic structure. Convert each
  motif afterwards. Preserve cut identity and chemical layer order across
  modes; do not rerun periodic symmetry detection on the projected box as
  though it had ordinary vertical periodicity.
- A projected box has physical lattice translations `A`, `B`, and `T`.
  Algorithms assuming translations `A`, `B`, and `C_surface` need adaptation
  or an explicit unsupported-operation error.
- New files must retain enough precision for repeat phases and round trips.
  Use explicit format/capability metadata. Old files remain readable; do not
  promise that old program versions can read new nonvertical-repeat files.
  Exporters unable to represent stacking must convert exactly or reject the
  export. They must not silently omit the translation.

## Stage 0 Establish reference cases and invariants

Deliverable: focused scientific fixtures and independent reference sums.

1. Preserve baseline calculations for current lattice mode, including existing
   rutile and STO examples, without refreshing their expected outputs merely
   to accommodate a changed result.
2. Add a synthetic oblique cell with analytically known shear, fractional
   in-plane scattering points, and explicit finite atom sums.
3. Add a portable TbScO3 geometry fixture using the published lattice constants;
   keep tests independent of a user's Zotero or research directory. Include a
   small source-property fixture for partial occupancy and displacement import.
4. Define comparisons at identical Cartesian `Q`, with aligned origins and
   lateral-area normalization. Raw amplitudes from different index conventions
   or different motif representatives are not interchangeable references.
5. Cover a nonrectangular in-plane cell and a case with zero lateral shift.

Exit criteria: fixtures demonstrate the expected phase and coordinate
identities independently of production stacking helpers. Record current test
failures separately from feature regressions.

## Stage 1 Add repeat state and persistence

Primary areas: `CTRuc.py`, `CTRstacking.py` where cycle state needs it,
`CTRsymmetry.py`, serialization and DWBA geometry fingerprints.

1. Introduce validated repeat state and one canonical conversion from
   fractional repeat to physical Cartesian repeat. Keep existing default
   arithmetic where practical to minimize numerical drift.
2. Preserve state through layer splitting, cloning, deep copy, optimizer
   copies, and multiprocessing serialization.
3. Add file round trips and loss-prevention checks for legacy exporters.
4. Include repeat state in all relevant geometry caches and stale-preparation
   checks. Avoid mutable array aliases that bypass invalidation.
5. Add temporary capability guards to consumers not yet migrated.

Exit criteria: legacy files and vertical-repeat cells retain their behavior;
nonvertical state survives supported copies and saves; unsupported operations
fail explicitly. No new default is enabled.

## Stage 2 Generate surface cells and update bulk scattering

Primary areas: `_CTRsurface.py`, `CTRsymmetry.py`, `CTRuc.py`,
`_CTRcalc_accel.py`, `cpp/CTRcalc_cpp.cpp`.

1. Add explicit `cell_mode="surface"` and construct `S`, `T`, transformed atom
   representatives, and source metadata from the original periodic motif.
2. Convert Wyckoff affine expressions and their coordinate derivatives along
   with positions. Retain source occupancy and isotropic-equivalent B values.
   Preserve Cartesian ADP tensors in metadata under frame rotations.
3. Use the complete repeat phase in Python, Numba, and C++ bulk sums, preserving
   attenuation units and the existing repeat-zero convention.
4. Audit coherent-domain behavior: currently the bulk lattice sum follows its
   documented repeat convention outside the motif domain sum. Preserve legacy
   behavior; reject configurations requiring differing bulk repeat vectors
   until they have an explicitly supported formulation.
5. Treat anisotropic displacement semantics explicitly. Legacy `iDW`/`oDW`
   values must not silently acquire a different physical axis. In surface mode
   document axes relative to the surface normal. Use isotropic cases for exact
   old/new equivalence tests; test directional cases against a Cartesian
   reference and reject unrepresentable tensor conversions.

Exit criteria: Cartesian positions and bulk complex amplitudes agree with
independent references, including off-specular and fractional in-plane points.
Both reference transforms and all available backends agree. Rebuild native
extensions before native validation. Source-property and Wyckoff tests pass.

## Stage 3 Support finite films and cell operations

Primary areas: `CTRuc.py`, `CTRfilm.py`, symmetry metadata helpers.

1. Generate whole-repeat offsets using the full translation, reusing existing
   affine domain translations. Select partial layers from the physical cycle
   without assuming that all interlayer spacings are equal.
2. Distinguish continuation through a crystal repeat from a rigid normal
   displacement of a component. Changing a vacuum/interface gap does not
   implicitly advance its lateral registry.
3. Update termination changes and layered translations. Preserve the absolute
   origin and phase when changing the representative repeat.
4. Update supercell positions and repeat metadata together. For diagonal
   repeat counts `(Na, Nb, Nc)` in the surface coordinate basis, the new
   fractional repeat is `(Nc*sx/Na, Nc*sy/Nb, 1)` and atom copies include the
   intervening lateral offsets. Audit Wyckoff coupling and image offsets.
5. Define strain behavior through the actual deformation: apply it to both
   motif positions and the physical repeat. A fit that only stretches the
   normal direction keeps lateral registry fixed; regenerating from a changed
   parent lattice recomputes the derived projection. Do not mix those operations.

Exit criteria: explicit Cartesian sums match generated finite films for zero,
one, several, and partial repeat counts; termination shifts preserve intended
registry; supercells preserve appropriately normalized amplitudes; fit updates
are reflected after cloning and caching.

## Stage 4 Support interfaces and roughness constructions

Primary areas: `EpitaxyInterface`, `PoissonSurface`, stacking state propagation
in `CTRcalc.py` and `CTRstacking.py`.

1. Carry enough registry information through component boundaries to continue
   a cycle correctly, alongside layer identity and physical height. Do not
   infer registry from height alone when different cycles or origins exist.
2. Update interface layer generation, lower-bulk subtraction, registration,
   strain interpolation, and rigid offsets using full position vectors.
3. Update terrace-height candidates and termination corrections in the current
   Poisson model so each selected height has the corresponding registry.
4. Retain existing probability and coherence conventions. This stage does not
   introduce Petach terrace-width distributions or step-edge jaggedness.

Exit criteria: sharp-interface additions/subtractions cancel in both positions
and complex amplitudes; zero-roughness limits recover the corresponding flat
film; equivalent descriptions of a termination agree; existing interface and
roughness regression tests remain valid in lattice mode.

## Stage 5 Complete density and DWBA support

Primary areas: `CTRuc.py`, `CTRoptics.py`, `CTRdwba.py`,
`cpp/CTRcalc_cpp.cpp`, `cpp/CTRdwba_cpp.hpp` as required.

1. Repeat density profiles by physical normal height and include lateral
   Fourier phases for nonzero in-plane indices. Lateral shifts leave the
   planar average unchanged when normal geometry is unchanged.
2. Use the stored physical repeat in DWBA preparation. Its native kernel
   already accepts a three-component repeat and evaluates `Q dot T`; reuse
   that mechanism rather than multiplying a separate correction into it.
3. Verify optical layer depths, volume/area normalization, normal strain,
   attenuation, orientation constraints, and interface boundary placement.
4. Include repeat changes in prepared-state validity and cache keys. Validate
   complex normal wavevectors and depth-dependent wavefields using the actual
   translated positions, not a kinematic intensity multiplier.
5. Remove capability guards only for paths covered by the completed tests.

Exit criteria: density integrals preserve areal electron content; shifted
stacks have the expected in-plane Fourier phase; DWBA agrees with explicit
finite/reference calculations and its existing kinematic-limit conventions;
changing repeat geometry cannot reuse stale prepared results.

## Stage 6 Enable the default and migrate examples

Primary areas: `_CTRsurface.py`, `examples/CTR`, existing STO research examples
when included in the implementation task, topic documentation and changelog.

1. Make Miller/in-plane construction default to surface mode. Keep full-matrix
   input in lattice mode unless surface mode is explicitly requested.
2. Update examples to state their intended mode. Preserve deliberate legacy
   fits with explicit lattice mode and add a perpendicular TbScO3 example
   showing both geometry and stacking translation.
3. Document that the coordinate `c` need not be a bulk lattice translation,
   that gamma need not be 90 degrees, and that reciprocal indices change with
   the coordinate basis. Show how to compare the same Cartesian scattering
   points across modes.
4. Update `doc/source/wyckoff_fitting.rst` and the relevant structure-factor
   documentation; add a compatibility/migration changelog entry and regenerate
   `doc/source/release_notes.rst` using its generator.
5. Treat the default change as a behavior change even though the physical
   structure is preserved: users supplying unchanged numeric HKL arrays may
   otherwise evaluate different physical points. Follow the repository's
   breaking-change convention if released compatibility is affected.

Exit criteria: the complete supported workflow runs in surface mode; examples
execute; legacy mode reproduces baseline results; docs describe actual
behavior; the default switch introduces no silent unsupported paths.

## Stage 7 Optional miscut model as a separate feature

This stage is not required to release perpendicular surface cells.

[Petach et al., Physical Review B 95, 184104 (2017)](https://doi.org/10.1103/PhysRevB.95.184104)
provides the motivating vector phase: Eq. (6) is
`theta = qx*M*a - qz*a = Q dot S`, with terrace translation
`S = (M*a, 0, -a)` in that paper's cubic coordinate frame. Our generalized
crystal-repeat formulation uses the same translation identity, but `S` and `T`
represent different physical operations.

1. Define terrace geometry, termination sequence, coherence, width statistics,
   and step-edge statistics independently of the bulk repeat. A CIF-derived
   lateral repeat component does not imply physical miscut or roughness.
2. Begin with deterministic identical terraces and coherent amplitude sums;
   compare with the paper's ideal limit. For alternating terminations, retain
   separate motif amplitudes instead of assuming one common factor.
3. Add disorder averages and the sharp, broad, and diffuse intensity components
   only after reproducing the paper's cubic reference cases.
4. Expose the measured quantity explicitly: individual subrod, integrated sharp
   peaks, or sharp plus broad scattering. Include coherence and detector
   acceptance/resolution where required.
5. Derive the noncubic generalization before using it. Do not transplant the
   integrated roughness factors as universal multipliers: their assumptions
   include wide terraces, appropriate coherence, and slowly varying structure
   factors over the integrated peaks. Avoid double-counting morphology already
   represented by an existing roughness model.

Exit criteria: reproduce appropriate published limiting cases and independent
real-space sums; document integration and coherence assumptions. This feature
has its own release scope and cannot delay stages 0 through 6.

## Validation and release discipline

Run narrow tests at each stage, expanding only when new paths are involved.
The relevant existing suites include `test_CTRsurface.py`,
`test_CTRsymmetry.py`, `test_CTRcalc.py`, `test_CTRfilm.py`,
`test_CTRinterface_registration.py`, `test_CTRinterface_anchors.py`,
`test_CTRoptical_profile.py`, `test_CTRdwba.py`, and `test_CTRopt_dwba.py`.
Add focused regressions to the corresponding areas rather than tests that
merely restate a helper's implementation.

- Use mixed absolute/relative tolerances for complex amplitudes near zeros.
  Compare away from singular unattenuated Bragg points; test limits separately.
- Include nonzero shifts in both lateral directions, fractional HKL, multiple
  origins, nonuniform layer spacing, and an unchanged zero-shift control.
- Check Python, Numba, and rebuilt C++ paths independently where available;
  report unavailable backends rather than claiming their validation.
- Run local Ruff checks, then `ruff check orgui` at the integration gate.
  Execute affected examples and build Sphinx after user-facing documentation
  changes. Do not run unrelated detector/reconstruction suites solely because
  this feature changes CTR geometry unless shared code is also modified.
- Use stage-sized Conventional Commits with the `phys` scope for scientific
  behavior changes. Keep unrelated working-tree changes out of the work.

The release gate is scientific equivalence at fixed physical coordinates plus
correct real-space stacking, not a visual check that the cell angles are 90
degrees. A failed gate keeps surface mode opt-in; it never triggers an implicit
strain, a silent fallback to vertical stacking, or deletion of repeat metadata.

# Reciprocal-space volume weighting

Status: first production and synthetic-validation stages implemented. The
root-cell centred-secant Jacobian mode is available explicitly as
``reciprocal_volume_average``; ``parameter_average`` remains the compatibility
default. Independent per-leaf and Gaussian-quadrature reference paths now
bound its error on small cases. Representative real-job validation and any
adaptive hybrid remain pending.

## Decision summary

The reciprocal-space reconstruction should continue to produce a local
corrected intensity field rather than raw counts per output voxel. For a
continuous scan, however, the average in each voxel should be weighted by the
reciprocal-space volume represented by each detector-pixel/exposure cell.

For detector-cell coordinates `(u, v)` and normalized exposure position `t`,
let `y(u, v, t)` be either Cartesian `Q` or `(H, K, L)`. The source-cell
volume is

```text
delta_V = integral |det(dy / d(u, v, t))| du dv dt.
```

For the part of source cell `i` mapped into output voxel `v`, call that volume
`delta_V_iv`. The proposed estimator is

```text
I_v   = sum_i(delta_V_iv I_i) / sum_i(delta_V_iv)
var_v = sum_i(delta_V_iv^2 var_i) / sum_i(delta_V_iv)^2,
```

where `I_i` has already had background, mask, detector, exposure, and monitor
corrections applied. Leaves from one source pixel must still be combined before
global accumulation so that subdivision does not make correlated pieces look
independent.

The first production implementation reuses the existing `weight` accumulator
for reciprocal volume rather than adding another floating-point field to every
record. A cell-level approximation is the low-cost first step. Full per-leaf
volumes should be an accuracy mode or an automatic fallback after the cell
approximation has been measured against them.

A genuinely stationary exposure has zero three-dimensional geometrical volume.
It must remain on the existing two-dimensional/parametric path until a physical
resolution thickness is defined. Adjacent-frame spacing is a quadrature weight,
not an instrumental resolution width.

## Why this is needed

The current reconstruction correctly maps detector-pixel footprints and
exposure sweeps into regular reciprocal-space voxels. Its adaptive subdivision
conserves unit source-pixel weight, and finalization calculates

```text
I_v = sum(w_iv I_i) / sum(w_iv).
```

The weights are fractions of the parameter cell `(u, v, t)`: a stationary
split produces four children of weight `1/4`, and a moving split produces eight
children of weight `1/8`. They answer which voxels a source pixel overlaps, but
not how much reciprocal-space volume each overlap occupies. See
`orgui/datautils/xrayutils/cpp/reciprocal_reconstruction_cpp.cpp`,
`subdivide_stationary` and `subdivide_moving`, and the current accumulator
contract in `doc/source/reciprocal_space_reconstruction.rst`.

This is a common and useful BINoculars-style average-pixel map. It is not,
strictly, a volume average in reciprocal space when the transformation
Jacobian varies across detector position or scan angle.

Mapping and intensity normalization remain separate operations:

- the coordinate transform says where a measurement contributes;
- the detector solid-angle correction converts integrated pixel counts to a
  local differential intensity;
- reciprocal-volume weight says how much of the output volume that corrected
  local estimate represents.

Using solid angle and reciprocal volume together is not double counting. They
appear on opposite sides of the estimator: solid angle corrects the measured
dependent variable, while `delta_V` supplies the integration measure used to
average that variable in `Q` or `HKL`.

## Coordinate and unit definitions

The native kernel uses

```text
Q = k (outgoing_unit_ray - incident_unit_ray),
k = 2 pi / wavelength,
```

and then rotates `Q` into the selected sample-fixed frame. The `hkl` frame
applies `(UB)^-1`. These are the existing conventions in
`reciprocal_reconstruction_cpp.cpp::coordinate_at` and `HKLVlieg.py`.

Let

```text
y(u, v, t) = Q(u, v, t)                  for a Cartesian Q frame,
y(u, v, t) = (UB)^-1 Q(u, v, t)          for HKL.
```

The local Jacobian is

```text
J = [dy/du, dy/dv, dy/dt].
```

For one normalized pixel/exposure cell,

```text
delta_V = integral_[0,1]^3 |det(J)| du dv dt.
```

Its units are:

- `Angstrom^-3` for Cartesian Q frames;
- `r.l.u.^3` for `hkl`.

For sufficiently small cells, the first-order approximation is

```text
delta_V ~= abs((delta_y_u cross delta_y_v) dot delta_y_t).
```

The existing geometry-matched step estimator already constructs the three
finite-difference vectors and places them into a local Jacobian in
`orgui/reconstruction_job.py::estimate_geometry_steps`. It currently uses
`J J^T / 12` to estimate projected axis widths. Its determinant is the missing
source-cell volume diagnostic.

### Output voxel volume

For an orthogonal Cartesian Q grid,

```text
V_voxel_Q = delta_q1 delta_q2 delta_q3       [Angstrom^-3].
```

For an HKL grid, the numerical coordinate volume is

```text
V_voxel_HKL = delta_H delta_K delta_L        [r.l.u.^3].
```

The same HKL voxel has physical Q-space volume

```text
V_voxel_Q = |det(UB)| delta_H delta_K delta_L
          = (2 pi)^3 / V_direct_cell * delta_H delta_K delta_L.
```

`U` is a rotation, so it does not change the determinant. `B` and `UB` include
the repository's established `2 pi` convention.

The product of the three complete grid ranges is only the bounding-box volume.
It is not the volume measured by the experiment. The measured support is the
union of the source cells. The sum of their contributions includes repeated
coverage and can exceed the unique volume of a voxel; that multiplicity is
useful for averaging but must not be labelled unique coverage.

## Estimator

For source pixel/exposure `i`, first form corrected differential intensity

```text
R_i = (counts_i - background_i)
      / (incident_flux_i exposure_i efficiency_i polarization_i solid_angle_i).
```

The current implementation may use relative rather than absolute factors, in
which case `R_i` remains on a relative scale. That changes the global scale,
not the volume-weighting algorithm.

Let `delta_V_iv` be the reciprocal-space volume of the intersection between
source cell `i` and output voxel `v`. Then accumulate

```text
S_I(v) = sum_i delta_V_iv R_i
S_V(v) = sum_i delta_V_iv
S_e(v) = sum_i delta_V_iv^2 var(R_i).
```

Finalization gives

```text
I(v)   = S_I(v) / S_V(v)
var(v) = S_e(v) / S_V(v)^2.
```

This has three useful invariants:

1. A constant corrected scattering field reconstructs to the same constant,
   independently of detector distance, scan step, voxel step, or local
   Jacobian.
2. Subdividing a source cell without changing its mapped support does not
   change its total reciprocal-volume contribution.
3. A near-singular Jacobian contributes little volume rather than causing a
   division by a small `delta_V`.

Do not construct `R_i / delta_V_i` and then average it. That is numerically
unstable near a Lorentz/Jacobian degeneracy and represents a different output
quantity.

### Repeated observations and statistical weighting

Volume weighting defines the geometrical average. It is not necessarily the
statistically optimal combination of repeated measurements with very different
precision. The first implementation should keep the scientific change narrow
and use volume weights alone, matching the variance propagation above.

An eventual precision-weighted repeated-observation mode would need to define
its estimator explicitly. Inverse variance estimated from the same Poisson
count can bias an average, so it must not be added casually as another factor
in `weight`.

## Stationary and stepped scans

For a stationary exposure,

```text
dy/dt = 0
```

and therefore `delta_V = 0`. One detector image samples a two-dimensional
patch of the Ewald sphere. Assigning it to a layer of three-dimensional output
voxels does not give that layer a physical thickness.

Three cases must remain distinct:

1. **Continuous exposure with exact angle bounds.** The exposure sweeps a
   physical three-dimensional cell. Volume weighting is well defined.
2. **Step scan with stationary exposures at different motor positions.** The
   frame centres provide quadrature samples of a three-dimensional scan, but
   the spacing between centres is not the resolution of one exposure. A
   Voronoi interval may be used for numerical integration if it is labelled as
   such.
3. **Single stationary frame.** A three-dimensional thickness requires an
   explicit resolution model, including the relevant beam divergence, energy
   bandwidth, detector point-spread, sample mosaicity, and calibration
   covariance.

Until cases 2 and 3 have explicit contracts, `volume_weighted` reconstruction
should require a nonzero exposure sweep and give a clear error otherwise. The
legacy parametric-average mode remains available for visualization and
backwards compatibility.

The current geometry-step estimator uses adjacent-frame displacement when an
exposure is stationary. That remains a useful grid-step heuristic, but it must
not be reused silently as physical source-cell thickness.

## Kernel implementation choices

### Option A: one volume per source cell

Compute one volume from the eight pixel/exposure root corners and distribute it
through the existing leaf fractions:

```text
delta_V_iv ~= w_iv delta_V_i.
```

The root corners are already evaluated for a moving cell. A centred finite
difference can average the four opposing edges in each parametric direction
before taking the determinant. This avoids choosing one privileged corner.

That reuse applies at adaptive depths 1 and above. The depth-0 centre-only fast
path currently evaluates only one coordinate per pixel. A volume-weighted
depth-0 path would need extra geometry evaluations or an analytic Jacobian, so
the low percentage below does not apply to centre-only mapping.

Advantages:

- no additional coordinate transformations;
- one determinant per source pixel;
- existing `VoxelWeight` and record layouts can be retained;
- expected end-to-end cost is approximately 1-3%.

Limit:

- assumes the Jacobian is constant across one pixel/exposure cell;
- does not allocate warped-cell volume exactly between intersected voxels.

For typical detector pixels and short continuous exposures this approximation
is expected to be much more accurate than the measurement resolution, but that
is a hypothesis to measure, not a reason to omit validation.

### Option B: one affine volume per adaptive leaf

Carry the transformed coordinate with each `LatticeVoxel`. When a leaf is
accepted, estimate its volume from its four or eight corner coordinates and
store that volume in `VoxelWeight.weight`.

Advantages:

- follows local Jacobian variation;
- allocates volume at the same resolution used for voxel-boundary decisions;
- needs no coordinate evaluations beyond those already performed by the
  subdivision.

Costs:

- `LatticeVoxel` grows by one three-vector;
- recursion-stack and lattice working sets grow;
- one determinant or equivalent triple product is evaluated per leaf;
- expected end-to-end cost is approximately 10-25% at continuous Balanced
  depth.

This is the recommended quantitative reference implementation.

### Option C: polyhedral or Gaussian integration

An eight-corner hexahedron can be decomposed into tetrahedra, or `|det(J)|` can
be evaluated with Gaussian quadrature. These approaches are appropriate for a
small validation subset, not as the first production path.

The corner polyhedron is still an approximation to the smooth normalized-ray
and trigonometric scan mapping. High-order quadrature that adds coordinate
evaluations can dominate the kernel rather than improve a limiting error.

Expected cost ranges from 20-50% end-to-end for several tetrahedral triple
products per leaf to roughly 1.5-2.5x total runtime for quadrature that adds
many coordinate evaluations.

### Recommended staged implementation

1. Add a mode field with explicit `parameter_average` and
   `reciprocal_volume_average` values.
2. Implement Option A and expose cell-volume diagnostics.
3. [Implemented] Implement Option B in the subdivision benchmark/reference
   path.
4. [Synthetic complete; real jobs pending] Compare A and B on synthetic fields
   and representative real jobs.
5. Keep A as the production path if its difference from B remains below the
   statistical validation threshold. Otherwise use B, or select B only for
   cells whose Jacobian-variation diagnostic exceeds a measured tolerance.
6. [Implemented] Keep Option C as a validation oracle on bounded subsets.

An automatic hybrid is preferable to making every ordinary cell pay the full
per-leaf cost. Its switching metric and tolerance must come from the A/B
comparison; they are deliberately not fixed in this document before those
measurements exist.

## Record and file-format design

The current native record contains three doubles:

```text
weighted_intensity, weighted_variance, weight.
```

In `reciprocal_volume_average` mode, reinterpret them as

```text
sum(delta_V I), sum(delta_V^2 var), sum(delta_V).
```

This keeps the 40-byte record and all existing tree-merge machinery. The final
`weight` dataset becomes reciprocal-volume coverage with units determined by
the output frame. Add attributes that state:

- reconstruction averaging mode;
- weight quantity and units;
- whether exposure bounds were exact or inferred;
- source-volume approximation (`cell`, `leaf`, or validation quadrature);
- whether the solid angle was relative or absolute;
- schema/scientific-contract version.

Do not add a second double to every intermediate record merely to preserve the
old parametric weight. Restoring a 48-byte record is expected to reintroduce a
measurable memory-bandwidth penalty. If the old diagnostic is needed, make it
an explicitly requested output and measure its cost.

`contributors` retains its current meaning: the number of independent source
pixels that reached a voxel. It is not a volume and should not be used as one.

The new mode changes numerical results and weight units, so prepared jobs,
checkpoint parts, and final files need a schema/scientific-contract version
bump. A job must not resume across modes or volume algorithms.

## Performance estimate and first measurement

The root-cell implementation was benchmarked after rebuilding the native
extension on 2026-09-20. A 256 x 256 continuous-exposure synthetic detector,
HKL output, Balanced depth 2, one native thread, and 4096-pixel work blocks
gave these medians over six interleaved timed runs after warm-up:

| Mode | Median kernel time |
|---|---:|
| ``parameter_average`` | 286.027 ms |
| ``reciprocal_volume_average`` | 284.984 ms |

The observed difference was -0.36%, which is timing noise: for this case the
root-cell determinant has no measurable performance cost. This agrees with the
design expectation that all eight coordinates are already evaluated by the
continuous adaptive splitter. It is a kernel microbenchmark, not yet a
full-pipeline or multi-platform result, so retain a 1-3% planning allowance.

The remaining figures are predictions for algorithms not yet implemented,
based on existing reconstruction measurements.

The current continuous-exposure Balanced path performs about 149 coordinate
evaluations and costs about 4.29 microseconds per valid pixel in the subdivision
benchmark. Mapping has also been measured as the dominant stage in a
kernel-bound job. Conversely, record reduction is memory-bandwidth-bound.

| Volume algorithm | Estimated kernel cost | Estimated end-to-end cost |
|---|---:|---:|
| Root-cell determinant, existing leaf fractions | +1-4% | **+1-3%** |
| Adaptive hybrid | +4-12% | **+3-8%** |
| One affine determinant per leaf | +15-35% | **+10-25%** |
| Several tetrahedra per leaf | +30-70% | **+20-50%** |
| Added high-order coordinate quadrature | 2-4x | **1.5-2.5x total** |

The whole-job conversion assumes mapping remains 60-80% of wall time. Coarse
grids or slow storage reduce the visible percentage; fine continuous maps at
Balanced depth are the regime most likely to realize the upper end.

The first two rows assume a continuous scan at adaptive depth 1 or higher,
where the root corners are already part of subdivision. At depth 0, deriving a
cell volume from eight corners replaces the one-coordinate centre fast path and
could make the kernel several times slower. Centre-only mapping is already a
lossy visualization mode and need not support quantitative volume weighting in
the first implementation.

Adding another double to each intermediate record would grow it from 40 to
approximately 48 bytes. The inverse change, 48 to 40 bytes, previously improved
realistic mapping throughput by 12-24%. Reusing `weight` avoids that likely
penalty and avoids approximately 20% more scratch-record traffic.

Illustrative scaling for a job whose current runtime is `T`:

```text
root-cell:       1.01T to 1.03T
hybrid:          1.03T to 1.08T
per-leaf:        1.10T to 1.25T
high quadrature: 1.5T  to 2.5T
```

The later-mode estimates must be replaced by results from
`benchmarks/benchmark_reconstruction_subdivision.py` and the full-pipeline
benchmark once an implementation exists. Record coordinate evaluations,
leaves per pixel, records per pixel, nanoseconds per pixel, frames per second,
peak resident memory, scratch bytes, and final-file bytes.

## First synthetic accuracy measurement

An independent Python oracle now validates the native root-cell mode without
sharing its determinant or subdivision implementation. It provides:

- uniform dyadic leaves with a centred-secant determinant and exact mapped
  centre for each leaf (Option B reference);
- tensor-product Gauss-Legendre integration of ``abs(det(J))``, with numerical
  centred derivatives, on a bounded pixel subset (Option C reference);
- conditioning, normalized-determinant, leaf-variation, and orientation-sign
  diagnostics.

The analytic regression suite covers affine determinants, curved-map
convergence, a folded mapping hidden by the root corners, output rotations,
linear Q-to-HKL volume scaling, scan-direction reversal, and stationary zero
volume.

The first native comparison used an 8 x 8 detector tile at pixel origin
``(1000, 1000)``, 0.172 mm pixels, 500 mm sample-detector distance, HKL output,
and subdivision depth 2. Quadrature order 3 was evaluated on eight pixels.
The scan angle below is the continuous exposure width.

| Sweep | Max root/quadrature volume error | Max leaf/quadrature volume error | Root/leaf total volume | Map p99 relative intensity | Map p99 pull |
|---:|---:|---:|---:|---:|---:|
| 0.01 deg | 2.12e-8 | 1.70e-8 | 0.999999990 | 1.14e-4 | 0.00281 sigma |
| 0.1 deg | 5.30e-7 | 4.81e-8 | 0.999999519 | 1.15e-4 | 0.00255 sigma |
| 1.0 deg | 5.08e-5 | 3.19e-6 | 0.999952399 | 5.57e-6 | 0.000521 sigma |

At 0.1 deg and depth 3, the maximum leaf/quadrature volume error fell to
2.41e-8; the map p99 pull remained small at 0.00311 sigma. No orientation-sign
change occurred in any of these cells. At the 0.1 deg production-like case,
the root Jacobian condition number was 1.71--1.73 and the normalized
determinant was about 0.96, so this was a well-conditioned geometry rather than
a degeneracy test.

These results support retaining the root-cell implementation for short
continuous exposures on this synthetic geometry. They do not justify making
the mode the default: representative calibrated jobs, ill-conditioned
geometries, detector-distance comparisons, and full-pipeline measurements are
still required. The Python reference is intentionally slow (about 151 ms for
64 pixels at depth 2 versus about 1 ms in the native mapper) and is an accuracy
oracle, not a production Option B performance measurement.

Reproduce or extend the measurement with
``benchmarks/benchmark_reconstruction_volume_weighting.py``. Pass
``--synthetic-sweep-degrees`` for a synthetic sweep or ``--job`` and ``--frame``
for a prepared real job. The JSON report contains the exact native-module path
so a stale extension is visible.

## Numerical and scientific risks

### Jacobian degeneracy

The determinant approaches zero where the detector and scan directions fail to
span three independent reciprocal-space directions. Use `abs(det(J))` for
volume, but also record a conditioning diagnostic. Do not repair a degenerate
cell by imposing a minimum volume: that creates artificial reciprocal-space
support.

### Warped and folded cells

The mapping is smooth for the supported flat fixed detector and ordinary short
exposures. Nevertheless, the correct mathematical volume is the integral of
`abs(det(J))`, not the absolute value of an integral of a signed determinant.
If determinant signs change within one source cell, force refinement or flag
the cell rather than allowing cancellation.

### Detector geometry

Accurate reciprocal coordinates depend on beam centre, sample-detector
distance, pixel pitch, detector tilt/rotation, wavelength, sample angles, and
`UB`. A volume-weighted accumulator cannot repair calibration error.

The native reconstruction currently refuses a moving detector arm because the
fixed detector ray array cannot represent it. Volume weighting must preserve
that refusal until per-frame/per-exposure arm geometry is implemented. Detector
translations would additionally change the per-pixel solid angle.

### Resolution versus sampling

Source-cell volume is a sampling measure, not a complete instrumental
resolution model. The existing geometry-matched step estimator deliberately
omits beam divergence, energy bandwidth, detector point-spread, mosaicity, and
calibration covariance. Those terms should eventually enter a resolution
covariance or convolution kernel, not be hidden inside `delta_V`.

### Absolute scale

The current detector solid-angle array is normally relative. It corrects shape
across the detector but does not alone produce intensity per absolute steradian.
Comparisons across detector distances or configurations require an absolute
solid angle, calibrated incident flux, detector efficiency, external
transmission, and any experiment-specific acceptance factors.

## Validation plan

### Mathematical unit tests

1. An affine synthetic mapping must reproduce its analytic determinant exactly.
2. Applying a rotation to Cartesian Q must leave `delta_V` unchanged.
3. Transforming Q volumes into HKL must multiply by `1 / |det(UB)|`.
4. Reversing the scan direction must change determinant sign but not volume.
5. A stationary exposure must return zero three-dimensional volume.
6. Splitting a cell must conserve the sum of child volumes to the tolerance of
   the selected volume algorithm.
7. Constant corrected intensity must remain constant after arbitrary
   subdivision and overlapping coverage.

### Synthetic forward tests

1. Reconstruct a constant isotropic scatterer across the full detector. The
   corrected map must be flat versus detector position and distance.
2. Reconstruct analytic linear and Gaussian fields in Q and compare voxel
   averages with numerical reference integration.
3. Repeat with reversed scan direction, changed scan step, and changed output
   voxel step.
4. Exercise a locally ill-conditioned geometry and verify that it is flagged
   without infinities or artificial minimum volume.
5. Verify variance propagation with repeated Poisson simulations.

### Real-data convergence tests

1. Compare root-cell, per-leaf, and bounded high-order reference volumes.
2. Run Balanced and High subdivision on representative short windows.
3. Express differences in units of each voxel's propagated standard error,
   matching the existing subdivision validation.
4. Require the production approximation to stay below `0.05 sigma` median and
   below a separately reported tail threshold; also report voxels reached by
   only one method.
5. Check invariance under detector-distance/configuration changes where the
   same physical reciprocal-space region was measured.
6. Compare a reciprocal-space integrated reflection against the established
   direct rocking/ROI result after each route's documented normalization.

### Performance tests

Measure at least:

- stationary and continuous exposure models;
- depths 0 through 3, with Balanced as the production focus;
- a fine kernel-bound grid and a coarse correction/I/O-bound grid;
- one and many native threads;
- kernel-only and complete mapping-pipeline throughput;
- peak memory and checkpoint size.

No performance prediction in this document should become a release claim
until those measurements replace it.

## Rollout and compatibility

1. [Implemented] Land the volume calculation behind an explicit mode.
2. [Implemented] Keep existing jobs on `parameter_average` unless migrated
   deliberately.
3. [Implemented] Version the job descriptor, checkpoint identity, final-file
   provenance, and `weight` units together.
4. [Implemented] Expose the mode in the GUI with the continuous-exposure and
   depth requirements stated explicitly. It remains opt-in while broader A/B
   validation is pending.
5. If the quantitative mode becomes the default, document the result change in
   the user-facing reconstruction page and changelog. Do not silently reinterpret
   existing saved maps.

## Open decisions

1. Whether the production path is root-cell, per-leaf, or an automatic hybrid.
2. The hybrid's Jacobian-variation metric and measured threshold.
3. Whether step scans receive a labelled Voronoi quadrature volume or remain on
   the legacy path until a resolution kernel exists.
4. Whether absolute solid angle is required for the first volume-weighted mode
   or remains a separate calibrated-scale option.
5. Whether unique coverage fraction is worth its additional geometric and file
   cost; summed reciprocal-volume coverage is sufficient for the estimator.
6. How an eventual statistical repeated-observation weight composes with the
   geometrical volume weight without introducing count-dependent bias.

## References

- S. Gorfman *et al.*, "Measurements and scaling of X-ray total scattering
  from single crystals", *J. Appl. Cryst.* **59** (2026) 1224-1237. The reciprocal
  volume of one detector-pixel/rotation cell is defined by the determinant of
  `d(Qx, Qy, Qz) / d(x, y, phi)`, and voxel intensity is normalized by the sum
  of those volumes. DOI: <https://doi.org/10.1107/S1600576726006187>.
- S. Roobol *et al.*, "BINoculars: data reduction and analysis software for
  two-dimensional detectors in surface X-ray diffraction", *J. Appl. Cryst.*
  **48** (2015) 1324-1329. Describes the established average-pixel reciprocal
  binning model: <https://doi.org/10.1107/S1600576715009607>.
- D. Kriegner, E. Wintersberger and J. Stangl, "xrayutilities: a versatile tool
  for reciprocal space conversion of scattering data recorded with linear and
  area detectors", *J. Appl. Cryst.* **46** (2013) 1162-1170. Documents exact
  flat-detector geometry and the importance of detector misalignment:
  <https://doi.org/10.1107/S0021889813017214>.

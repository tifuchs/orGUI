# Sample interception: numerical and persistence contract

The implementation follows the refined MISC-011 request and plan under
`artifacts/sample-interception/`. Units at the core boundary are metres,
radians and normalized inverse-metre densities. UI shape dimensions are mm;
profile widths/offsets are micrometres. Internal geometry is
`R(psi) [c + R(theta0) u]`, with sample-fixed displacement c. Axis coordinates
in the centred profile frames are `-c_x sin(alpha_ref)` and `-c_y`; the alignment
stays fixed in the laboratory. Profile offsets then retain their existing
positive sample-displacement sign. The application requires an explicit,
confirmed normal-rotation motor name/unit/sign/reference or an explicit fixed
angle. Six-circle solver angles are never guessed as that motor.

H is evaluated directly as a surface integral. A polygon is sectioned in x,
all transverse interior intervals are integrated with the horizontal CDF,
and vertical density weights the remaining one-dimensional quadrature. Circle
sections are exact. Aligned rectangles use a product of interval probabilities
except for tiny projected intervals, which retain direct surface integration
to avoid cancellation. Vertex x values, vertical knots/quantiles and horizontal
transitions crossing edges split quadrature. Characteristic quantiles guide
integration; they are not truncation limits. Measured tables alone have finite,
zero-exterior support. Numerical error applies to H and excludes uncertainty
in measured beam profiles and geometry. Exact duplicate angle/alignment tuples
are cached within one call, without rounding. A callback permits cancellation
between frames; no UI enters the numerical module.

Reconstruction preparation now resolves and validates inputs without evaluating
H. It freezes reference incidence, horizontal measured data and selected custom
azimuth snapshots as before. Preview validates inputs and reports pending
evaluation; job status uses settings or recorded provenance without reopening
the scan to resolve factors. Full quadrature, zero-overlap checks and application
remain at execution start, with progress/cancellation callbacks. No prepared
job schema or correction convention changes. This removes the per-frame
quadrature cost from configuration iterations on large scans; asset snapshots
and bounded checkpoint estimation still run during preparation.

Polygon quadrature retains physical x breakpoints and maps each piece onto a
separate unit interval. The crossing-edge pairs are fixed within a piece;
their transverse heights are interpolated directly in its local parameter.
This avoids repeatedly subtracting almost equal global x coordinates and
dividing by a tiny edge x span near beam alignment. Merely rescaling QUADPACK's
integration interval would leave that ill-conditioned intersection calculation
in place. No acquisition angle is snapped to an aligned value, and the surface
Jacobian, probability definitions and geometry conventions are unchanged.

Total-flux uses Q and H framewise. Legacy shape density uses peak reference
`Aeff = H/(pz_peak ph_peak)`, divided by geometric sample area for the relative
correction. It is distinct from the total-flux fraction. Existing 1D functions
and their endpoint/interpolation behavior are preserved. The new public
measured-density API rejects negative data and uses zero exterior support.
Old custom subclasses still instantiate, and report unsupported 2D capabilities
when called. Horizon/invalid angles are NaN, zero overlap has f_hit=0 and a NaN
application divisor, and no usable frames fails. Overspill warnings aggregate
instead of producing per-frame dialogs.

CorrectionState gains optional `sample_interception` JSON. All previous typed
fields retain their units and dataset names; an additional NeXus settings group
carries the version-1 model and requested-settings layout becomes 4 only when
the model exists. Applied curve records store model JSON, embedded normalized
horizontal data, raw readback radians, validity, H error and method; existing
vertical provenance and immutable signal/variance remain authoritative.

**Compatibility finding that changed the plan:** advancing the curve container
schema to 4 was unsafe. The previous rocking reducer guards only schema 3 and
would otherwise fall through to its old ROI branch. Keep the curve schema/group
at 3 and use distinct `shape_interception_*_v1` algorithms: the old guard now
refuses reinterpretation. A regression executes the actual prior guard from
the baseline revision against a real HDF5 shape record. Requested-settings
version 4 is independent of that applied curve container.

Keep/remove/replace reconstruct from immutable base data. Replacement uses
recorded source incidence/readback, slices per-curve data before computation,
and preserves Q and the selected convention. A changed moving source or unit
requires extraction. The reduction output records action/convention/source and
current embedded settings when replaced; the source record remains unchanged.

Independent tensor Gauss-Legendre surface integration, analytic interval
probabilities and uniform-area references verify rotating offsets, symmetry,
grazing limits, winding/concavity and probability normalization. The reviewed
10 mm square / 0.36 degree / 160 micrometre Gaussian references are retained
as explicit regressions, including the 20 mm top-hat distinction. End-to-end
stationary and rocking tests check both density and total-flux conventions,
error propagation, HDF5 persistence and bounded replacement reads.

Initial support uses representative frame angles; it makes no exposure-motion
quadrature claim. Flat surfaces and separable independent beam densities are
assumptions. Reconstruction remains without illumination. Native ROI acceleration
is optional: the wrapper's `HAS_ACCEL_BACKEND` must be checked, because importing
the wrapper alone succeeds even when its compiled functions are unavailable.

## Validation and timing (2026-10-01)

The final feature/configuration gate passed 201 tests (5 optional accelerator
skips). Reconstruction, HKL, detector, CTR and backend regression suites passed
367 tests (13 optional-backend skips, 131 subtests). A final scalar-incidence
and frame-policy check passed 18 tests. The independent calibrated pixel-ray
forward test with a rotating, displaced shape passed and recovered F² within
3e-10 relative tolerance; its tensor surface integral converged between orders
80 and 120 before use as the forward reference. Ruff and diff whitespace checks
passed. Sphinx built successfully; warnings refer to existing topic pages,
notebooks and release-note references, with none in the added topic section.

The broad app/profile run initially passed 533 tests and identified 13 missing
accelerator checks plus one new smoothed-top-hat issue. Its open-quantile and
infinite-interval behavior was fixed and all analytical-family checks passed.
With explicit `ORGUI_ACCEL_BACKEND=numba` and a writable local
`NUMBA_CACHE_DIR`, the ROI/extraction/correction run passed 58 tests, including
all scan extraction tests. Six unrelated repair/polynomial-background tests
still require the absent C++ backend. No production environment setting changed.

On this Windows/Miniforge environment, a 37-angle square/Gaussian/top-hat sweep
took roughly 1.8–2.0 s; a 101-point measured vertical profile took 2.6–3.3 s.
One thousand exactly repeated geometries took 0.05–0.07 s. Maximum reported
quadrature error on H in those runs was about 3.2e-13. These are observed local
timings, not a cross-machine guarantee. Large measured tables or polygons cost
more; collector-side progress/cancellation is available for extraction and
preview. The UI preview was rendered and inspected, and a focused widget test
checks source-frame selection, radian/degree boundaries and plot API behavior.

## Sequential performance follow-up (2026-10-01)

Gaussian density now uses its normalized exponential directly, and Gaussian
interval probabilities reuse the frozen mean and standard deviation. Uniform
profiles use the exact clipped interval length divided by beam width. Other
analytical families and custom profiles retain their general distribution
path. Frozen distribution quantiles are cached in raw profile coordinates;
each access returns a fresh array shifted by the current sample centre.
Horizontal breakpoints are fetched once per frame rather than once per edge.

Circular quadrature no longer uses the vertices of the display outline. It
uses `x = cx + R sin(theta)`, with `theta = pi (t - 1/2)` and exact half-chord
`R cos(theta)`. Its normalized Jacobian is `(pi/2) cos(theta)`; multiplying by
the final `2 R pz_peak` scale preserves the original surface integral and H
units. This removes the square-root endpoint singularity and cancellation
for narrow horizontal beams. Physical profile transitions still split the
quadrature, and integration tolerance/convergence checks remain unchanged.
No multiprocessing or threaded frame execution was added.

The paired Windows/Ryzen AI 9 HX 370 benchmark uses 37 synthetic frames,
three warmed repetitions, Python 3.14.7, NumPy 2.5.3 and SciPy 1.18.0, with
`rtol=1e-9`. Median application-policy elapsed time divided by frame count
excludes image I/O, normalization and ROI integration:

| Case | Before (ms/frame) | After (ms/frame) | Speedup |
| --- | ---: | ---: | ---: |
| Aligned square | 0.2465 | 0.0786 | 3.1x |
| Rotating square | 28.7072 | 2.3769 | 12.1x |
| Rotating displaced rectangle | 32.6564 | 8.0198 | 4.1x |
| Circle | 175.4348 | 4.2098 | 41.7x |
| Rotating concave polygon | 67.3011 | 17.7160 | 3.8x |
| Square, 101-point measured vertical profile | 42.5613 | 4.5230 | 9.4x |
| Square, 501-point measured vertical profile | 121.3759 | 13.3423 | 9.1x |
| Identical rotated geometry repeated | 0.8770 | 0.1228 | 7.1x |

The reproducible driver is `benchmarks/benchmark_sample_interception.py`;
machine-specific JSON captures live in ignored `benchmarks/baselines/`.
The correction/integration validation passed 319 tests with the optional
Numba ROI backend enabled only for that test process, and Ruff passed.
New independent references cover Gaussian disk probability, narrow uniform
strips through a disk, displaced circles integrated in polar coordinates,
profile offsets, support edges, broadcasting and tail/tiny probabilities.
All benchmark geometries were also compared with the original HEAD source;
the maximum relative change in H was `1e-15`.
The package-export assertion was also brought up to date with the existing
public `sample_interception` module.

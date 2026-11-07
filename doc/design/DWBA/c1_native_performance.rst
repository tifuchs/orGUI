C1 native DWBA performance audit
===============================

Scope and numerical contract
----------------------------

This audit uses the unchanged LaNiO3 0056/C1 crystal and JSON sidecar supplied
with ``installed_api_verification/PERFORMANCE_REQUEST.md`` on 2026-11-07.
The standalone loader checks both source hashes and restores the STO reference,
20 keV energy, 0.34 degree incidence (converted to radians), empirical kinematic
attenuation, native distribution and selected optical policy. No fitting or
data rewriting is performed. The 512-point rod uses one identity domain;
the saved 12-point sample uses all four rotations, both exit polarizations
and both bulk attenuation choices, 0 and 0.005.

The optimization is confined to the native layered record-amplitude kernel.
It does not change terrace selection, field solving, optical profiles,
probability normalization, ambient phase alignment, partial coherence,
polarization/domain intensity combination, or resolution convolution.
Positions are Cartesian Angstrom and momentum components are inverse Angstrom.

For each output record, identical scattering states at the same exact depth
and in the same medium share a lateral sum:

.. math::

   \sum_j w_j f(Q_\parallel,Q_z)
      e^{i(Q_x x_j+Q_y y_j+Q_z z)}
   = f(Q_\parallel,Q_z)e^{iQ_z z}
      \sum_j w_j e^{i(Q_x x_j+Q_y y_j)}.

Scattering-state keys retain the complete form-factor row, both anisotropic
DWF parameters and occupancy. Position keys compare floating-point bit
patterns; there is no closeness tolerance. Signed weights stay in the lateral
sum. The terminal bulk numerator has a separate accumulator and uses the
original semi-infinite denominator. Explicit bulk attenuation is included in
each repeat's weight, independently of optical absorption.

All new tables are local to a single native invocation. Lateral sums are
reused only while the in-plane momentum is bit-identical, including within
a rod. Their scratch storage does not grow with the scan-point count; it
does grow with the packed finite atoms and explicit substrate repeats.
No fitted-parameter cache persists across calls or batches. There is no new
public batching API. Existing batches still reuse terrace copies, atomic
packing and optical fields across channels and rods.

Private terrace states also register ``clear_cache`` with the batch's
``ExitStack``. Their crystal/state reference cycles can outlive the batch;
the callback releases the large optical field tables immediately, including
on exceptions. The wrapped user's crystal keeps its persistent field cache.

The regrouping changes floating-point summation order. Bitwise identity with
the old atom loop is therefore not promised; comparisons retain the requested
``rtol=1e-9, atol=1e-28`` for optical quantities. No tolerance is relaxed.

Reproduction
------------

``benchmarks/benchmark_dwba_c1.py`` accepts ``--request-dir`` and a new
``--output`` JSON path. It writes all endpoint arrays to an adjacent NPZ.
``--native`` and ``--ensemble`` optionally select the saved native extension
and ensemble Python module for A/B execution without replacing the user's
installation. ``--compare`` points to
the baseline JSON/NPZ and checks every optical array, including parameter
proposals. Run the before and after commands sequentially.

Each policy has three fresh 512-point calls and two calls within one batch.
Five sequences each make three calls outside batches: unchanged parameters,
incidence increments of 0.00001 degrees, DWF increments of 0.0001 Angstrom
squared, a common surface-coverage factor of 0.9999, and a surface-variance
factor of 0.9999. The coverage proposal preserves oxygen/cation ratios; the
variance proposal moves inward from the saved upper support limit. These
are legal representative perturbations, not a reconstruction of the historical
27-coordinate optimizer or a refit. The saved model has no active fit-parameter
list, so the benchmark updates its native parameter bases and refreshes them
through ``setFitParameters([])`` outside the evaluation batch.

The process benchmark runs eight complete 512-point state-specific DWF
proposals at each of 1, 2, 4 and 8 workers. All BLAS/OpenMP environment limits
are set to one before NumPy imports. Reports separate model loading from
evaluation and include pool wall time and per-process peak working sets.
The reported sum of worker peaks is an upper bound, not a synchronized
measurement of total resident memory.

``benchmarks/profile_dwba_c1.py`` writes cold/warm cProfile files, readable
function tables and the native kernel's private per-thread work counters.
The latter count input instances, grouped motifs, phase exponentials and
evaluated motif/channels; they avoid noisy timing assertions in regression
tests.

Measured results and validation
-------------------------------

The host has 12 physical / 24 logical CPU cores and 33.36 GB RAM, running
Windows 11 and Python 3.14.9. The native build uses MSVC 19.51, release
optimization, without additional BLAS/OpenMP workers. The original version
label remains ``1.5.1.dev382+ga5a6afe16.d20261107``. Actual imported paths and
SHA256 hashes are recorded in ``output/dwba_performance/before.json``,
``after.json`` and ``installed_audit.json``. The installed wheel's native
binary hash matches the benchmarked checkout binary exactly.

Host load changed materially during the long suites: available RAM dropped
to about 1 GB and the baseline's later calls slowed down. Consequently,
``benchmark_dwba_c1_paired.py`` alternates fresh subprocesses, using A/B,
B/A, A/B order, and compares all returned optical arrays for every pair.
These matched measurements are the primary speed comparison; they are not
ratios against the historical 10.648-second measurement in the request.

.. list-table:: Alternating subprocess medians, three samples, seconds
   :header-rows: 1

   * - Workload
     - Before
     - After
     - Speedup
   * - Fresh 512-point state-specific call
     - 16.916
     - 7.867
     - 2.15x
   * - Legal DWF proposal, 512 points
     - 16.121
     - 7.422
     - 2.17x
   * - 12 points, four domains, two channels, two damping controls
     - 16.983
     - 13.047
     - 1.30x

Fresh-process model loading/proposal setup takes approximately 0.065 seconds,
reported separately; evaluation timings include the native batch lifecycle.
The small multi-channel workload does not reach the 2x target. Its remaining
limit is mostly repeated Python optical-profile/state construction, not
crystal copying alone. Sharing unchanged component profiles across terraces
would be a separate, carefully validated follow-up; no such persistent model
cache is introduced here.

For transparency, the earlier complete-suite medians are below. They use
identical workloads before/after but are separated in time and are affected
by the changing host load. The apparent small-sample slowdown in that suite
is why the alternating subprocess measurements were added.

.. list-table:: Complete suites, seconds
   :header-rows: 1

   * - Workload
     - Before
     - After
   * - State-specific fresh 512
     - 10.271
     - 8.732
   * - State-specific repeat inside batch
     - 8.309
     - 4.090
   * - Unchanged model, outside batches
     - 10.515
     - 7.324
   * - Geometry proposal
     - 10.525
     - 7.408
   * - DWF proposal
     - 10.540
     - 7.464
   * - Occupancy proposal
     - 10.478
     - 7.656
   * - Roughness/profile proposal
     - 10.777
     - 7.903
   * - Shared fresh 512
     - 0.667
     - 0.364
   * - Shared repeat inside batch
     - 0.547
     - 0.254
   * - 12-point multi-domain/channel sample
     - 10.954
     - 12.275

.. list-table:: Eight complete 512-point jobs, including pool startup
   :header-rows: 1

   * - Workers
     - Before wall seconds
     - After wall seconds
     - Before worker peak bound, GB
     - After worker peak bound, GB
   * - 1
     - 123.01
     - 71.72
     - 1.667
     - 0.594
   * - 2
     - 108.98
     - 38.30
     - 2.354
     - 1.065
   * - 4
     - 68.86
     - 25.72
     - 2.746
     - 1.889
   * - 8
     - 55.81
     - 22.80
     - 3.531
     - 3.532

On this host, four workers provide a useful memory/time compromise: eight
are only about 13 percent faster and nearly double the worker peak bound.
The callback reduces retained memory over successive jobs, not the first
job's required live field storage. Eight workers each run just one job here,
so their peak bound is essentially unchanged. The complete-suite parent
peak was approximately 1.773 GB before and 0.687 GB after. No cyclic-GC
settings, kernel weights, ray flags or support limits were changed.

.. list-table:: Cold 512-point cProfile breakdown, seconds
   :header-rows: 1

   * - Work
     - Before
     - After
   * - Crystal copies
     - 0.392
     - 0.385
   * - Atomic/optical-profile snapshots
     - 2.720
     - 2.835
   * - Reference object construction
     - 0.003
     - 0.004
   * - Field solving
     - 0.602
     - 0.548
   * - Atomic packing and factor-table extraction
     - 0.680
     - 0.779
   * - Native amplitude kernel
     - 12.415
     - 4.053
   * - Remaining geometry, field assembly and Python orchestration
     - approximately 1.0
     - approximately 1.0

These instrumented times are diagnostic, not the headline benchmark. The
native kernel still includes complex form-factor lookup/calculation. The
profiler's ``--kernel-only`` controls separate factor-table work from the
finite and substrate atomic sums without using reduced models for the
end-to-end benchmark. Both cold and repeated candidate calls process
113,576 atomic instances across the common model and all 31 active terraces,
grouped into 31,372 lateral motifs. Thus the improvement reduces work within
each evaluation; it does not depend on memoizing identical predictions.

On the mean-state preparation, warmed empty-instance controls cost 0.025
seconds before and 0.036 seconds after; additionally removing the atom/factor
tables reduces these to 0.022 and 0.020 seconds. This bounds factor/scattering-
state setup and cache lookup separately from the expensive atomic sums.
These small timings are noisy on this host and are not additive estimates
or speedup claims. Full raw controls and cache hit counts are retained in
``kernel_before/kernel.json`` and ``kernel_after/kernel.json``.

The full suite and paired runs compare mixed intensity, both endpoint
intensities and squared amplitude prefactors. Maximum observed relative
difference is ``7.70e-13``. The strict saved 12-point optical comparison
passes with the original tolerances. The change intentionally gives up
bit-for-bit identity with the old summation, while repeated calls to the
new implementation remain exactly identical.

Validation: 151 focused tests pass from an isolated wheel installation in
``output/dwba_performance/installed/site-packages``; another 358 core,
optimizer, optical-inventory, water and form-factor-cache tests pass from the
checkout (plus 81 subtests). ``ruff check orgui`` and the benchmark-script
lint checks pass. New tests independently sum complex atom phases, retain
signed weights and near-equal scattering states, check changing in-plane
momentum, count native work, and check cache release after normal and
exceptional batch exits. Existing tests cover coherent/incoherent endpoints,
Born and one-height limits, polarization, domain rotation and coincident
boundaries.

The supplied ``verify.py`` also passes against the asserted installed build:
the near-coincident-boundary reproducer, the saved four-domain sample, and
endpoint/batch/parameter checks on 45 points across nine rods under both
policies. Its saved-sample maximum optical relative error is ``5.60e-13``.
Its fresh state-specific 512-point median is 7.941 seconds and repeat is
4.741 seconds; shared median is 0.506 seconds and repeat is 0.299 seconds.
Its reported historical speedup is not used as a controlled A/B measurement.

Artifacts and commands
----------------------

All local measurement artifacts are under ``output/dwba_performance``:
``before.json/.npz``, ``after.json/.npz``, ``paired/paired.json``,
``profile_before_full``, ``profile_after``, ``summary.json``,
``installed_audit.json``, and the installed/broad pytest logs.
Saved baseline native and ensemble modules in that directory permit exact
A/B reproduction. They are local audit artifacts, not a second implementation
shipped in the package.

The final verification and copies of the reports are exported to the canonical
CTR directory
``data_out/1V2/0056/diagnostics/performance_trial_20261107_motifs_verified``.
An earlier elevated run resolved the old installation and is explicitly
excluded in ``performance_trial_20261107_motifs/NOT_CANDIDATE.txt``. The
``verify_dwba_c1_installed.py`` launcher asserts the install path and ensemble
hash before invoking the supplied, unchanged driver and checks its reported
hash afterward. It changes import routing only. The successful run is also
retained locally under ``output/dwba_performance/verification``.

Use the supplied request directory as ``REQUEST`` below, with the same
Miniforge Python executable as the request::

    python benchmarks/benchmark_dwba_c1.py --request-dir REQUEST --output NEW/before.json --native output/dwba_performance/baseline_CTRcalc_cpp.cp314-win_amd64.pyd --ensemble output/dwba_performance/baseline_CTRincoherent_dwba.py
    python benchmarks/benchmark_dwba_c1.py --request-dir REQUEST --output NEW/after.json --compare NEW/before.json
    python benchmarks/benchmark_dwba_c1_paired.py --request-dir REQUEST --baseline-native output/dwba_performance/baseline_CTRcalc_cpp.cp314-win_amd64.pyd --baseline-ensemble output/dwba_performance/baseline_CTRincoherent_dwba.py --output NEW/paired
    python benchmarks/verify_dwba_c1_installed.py --site-packages output/dwba_performance/installed/site-packages --request-dir REQUEST --output NEW/verification --ensemble-sha256 e9fce724e2f02c0f6150e9b6e6835050c8ae214b6db5b74ac77654b498a9a371

The wheel is installed into an isolated directory for validation; the user's
main Miniforge installation is not replaced. Source model hashes are checked
by the supplied loader on every load. Historical predictions and the original
``calculate.py`` checkpoints are not overwritten.

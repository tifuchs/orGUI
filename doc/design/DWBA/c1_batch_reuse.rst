C1 DWBA batch-reuse follow-up
============================

Scope
-----

This follows ``c1_native_performance.rst`` and compares against that report's
optimized implementation, not the older installed version in the original
request. The unchanged C1 loader restores the saved crystal, JSON sidecar,
20 keV energy, STO reference, 0.34 degree incidence and full height support.
No fitting, source-data changes or historical-output replacement is involved.
Lengths remain Angstrom, wavevectors inverse Angstrom and internal angles
radians. Optical policies and intensity-combination order remain unchanged.

Implementation
--------------

* Exact unit-cell optical profiles are shared across terrace snapshots within
  one batch. Keys cover basis/factors, energy, lattice/volume/area, coherent
  domains, physical-plane metadata, layer order and crystallinity labels.
  Unknown subclasses and instance-level optical overrides bypass reuse.
  Masked layer shares and bulk profiles have separate keys. Outer transforms,
  signed component weights, slab assembly and optical discontinuities remain
  state-specific. No tolerance-based merging is introduced.
* Small native scan blocks retain finite atomic sums by output record,
  optical medium and each of the four wave channels, before multiplication
  by optical coefficients. Exact inputs include positions, signed weights,
  atom indices, DWF/occupancy, form factors, interfaces, in-plane momentum,
  indexed normal wavevectors and which channels were actually evaluated.
  Changing polarization recalculates coefficients. Changing attenuation
  recalculates all substrate terms and reference subtraction. Neither control
  is folded into a multiplicative approximation. If separating optical and
  atomic factors would overflow or underflow, evaluation falls back to the
  original scaled multiplication order and does not retain that finite sum.
* Native storage is limited to one entry of at most 1 MiB per state and scans
  of at most 128 points. Larger requests use the complete uncached native
  instantiation. These limits affect retention only; no atoms, terraces,
  observations or resolution coordinates are removed.
* Explicit batches retain up to four assembled field requests per state.
  Standalone single-channel calls avoid these field/scattering tables while
  still sharing cell optics across heights. Nested scopes retain their
  caller's policy. All new data is released at outer scope exit, including
  exceptions. No new cross-proposal cache is introduced.
* The shared-field correction basis keeps its existing complete-amplitude
  reuse. It skips extra native/field scratch and cross-terrace profile keys
  when no complete terrace states are needed. Copying a crystal drops opaque
  native scratch rather than attempting to serialize it.

``cpp/CTRdwba_records.hpp`` contains one source of scientific arithmetic,
compiled with and without cache bookkeeping. The uncached expansion's
arithmetic and control-flow tokens match the frozen original kernel exactly;
``uncached_expansion.diff`` records the comments-only difference. This keeps
cache objects and their cleanup machinery out of the uncached compilation.

The existing ``model.dwba.batch()`` API enables reuse automatically. Keep
domain changes outside batches, and put outgoing-polarization and damping
requests for a domain inside the same batch. The benchmark uses exactly this
public pattern. No new public batching API or analysis-side physics is needed.

Measured result
---------------

Three alternating fresh-process before/after pairs, with BLAS/OpenMP threads
fixed to one, produced these median wall times. The speedup column is the
median of paired ratios, rather than the ratio of independent time medians.

.. list-table:: Paired timings
   :header-rows: 1

   * - Workload
     - Before (s)
     - After (s)
     - Speedup
   * - Fresh state-specific 512-point s-to-s
     - 8.105
     - 6.934
     - 1.09x
   * - Legal DWF proposal, fresh 512 points
     - 5.649
     - 5.182
     - 1.09x
   * - 12 points, four domains, two exits, two damping choices
     - 8.506
     - 6.533
     - 1.30x
   * - Shared-field 512 points (seven fresh models per process)
     - 0.318
     - 0.319
     - Approximately 1x

The multi-channel workload takes about 23% less elapsed time. The DWF
proposal comparison shows a more modest gain of about 9%. Fresh 512-point
results were sensitive to host load across successive series; a large,
reliable single-channel speedup is not established. The additional 2x
end-to-end target is not achieved. Shared-field throughput is preserved in
the longer comparison, which uses three alternating pairs with seven fresh
models per child process to reduce sensitivity to individual short calls.

The retained finite sums target small scan blocks. A 26,042-coordinate scan
submitted as one call still uses the complete uncached native kernel. A
streaming evaluator that processes controls together would be needed to
extend this finite-sum benefit to large grids without retaining large tables.

Paired 512-point optical arrays are bit-identical to the baseline. The largest
relative difference in the multi-channel arrays against the preceding
optimized baseline is 8.484e-14. The final installed verification differs
from the original saved optical predictions by at most 5.210e-13. Both pass
the unchanged ``rtol=1e-9, atol=1e-28`` checks. Repeated calls within a batch
remain bit-identical.

Fresh-process peak RSS is approximately 442 MB for both implementations on
512 points. The multi-channel sample grows from approximately 215 MB to
238 MB. Intermediate unbounded-cache experiments were rejected because they
increased fresh-call time and memory; ``paired`` and ``paired_bounded`` are
exploratory artifacts, not the delivered result. ``paired_delivered`` records
the final native kernel/state-specific comparisons; ``shared_confirmed``
records the final shared-only dispatch adjustment. Earlier ``paired_final``,
``paired_verified`` and ``paired_shared_final`` directories are intermediate
measurements and are not the delivered comparison.

The host has 12 physical/24 logical processors and about 33.36 GB RAM.
Background load varied substantially: separated runs must not be treated as
controlled speedups. Full parameter-change and worker sweeps are retained as
additional correctness/scaling evidence; paired measurements above are the
primary performance result.

The final full validation sweep used one fresh evaluation per workload;
the frozen baseline sweep used three. All unchanged, geometry, DWF,
occupancy and profile proposals passed optical comparisons. Their final
single-worker times were respectively 7.854, 8.956, 8.742, 8.196 and 8.932 s.
State-specific fresh/first-batch/repeated-batch times were
8.582/8.509/5.386 s; shared-field times were 0.340/0.308/0.197 s.
These separated observations are not controlled incremental speedups.

.. list-table:: Eight-job process scaling, raw full-sweep observations
   :header-rows: 1

   * - Workers
     - Before wall (s)
     - Final wall (s)
     - Final mean setup/job (s)
     - Final mean evaluation/job (s)
     - Final peak RSS sum bound (GB)
   * - 1
     - 89.109
     - 68.565
     - 0.069
     - 8.192
     - 0.594
   * - 2
     - 48.838
     - 37.290
     - 0.072
     - 8.515
     - 1.067
   * - 4
     - 34.210
     - 23.417
     - 0.084
     - 9.928
     - 1.894
   * - 8
     - 24.552
     - 16.409
     - 0.174
     - 12.456
     - 3.540

Memory is a sum of process peaks, an upper bound rather than a simultaneous
measurement. Baseline bounds were 0.593/1.066/1.891/3.534 GB. Eight workers
gave the shortest total time on this host, while individual evaluations
became slower as workers competed for resources. ``delivered_summary.json``
and the raw ``before.json`` / ``delivered_after.json`` retain these results.

Validation and artifacts
------------------------

Artifacts live under ``output/dwba_reuse``. The ``baseline`` directory freezes
the prior optimized native binary, DWBAState, ensemble facade and native
sources. The benchmark loader restores all three implementation modules,
avoiding an accidental mixture of baseline native code and candidate Python.

New regression tests cover exact cross-copy optical reuse, physical key
dependencies, all four polarizations and both damping choices, signed finite
weights, bulk/reference terms, one-bit invalidation, absent optical channels,
memory/point-limit fallbacks, nested scope policy and normal/exception cleanup.
Existing endpoint, Born, domain-rotation, partial-crystallinity and
near-coincident-boundary tests remain required.

The final isolated wheel passed 537 tests plus 81 subtests, including the
28 new reuse regressions. The wheel is under ``wheels_release`` and installed
under ``installed_release/site-packages``. The main Miniforge installation
was not changed. Actual imported paths and hashes, rather than the unchanged
version label, identify the tested implementation.

The original verification driver passed against that installation in a fresh
``verification_release`` directory, including its original failure
reproducer, saved 12-ray/four-domain predictions and endpoint, repeat and
post-batch mutation checks on 45 points across nine rods for both policies.
Its historical 26.104 s / current 10.424 s ratio is not a controlled
incremental speedup and is not used in the paired table above.
``source_hashes_release.json`` records matching checkout/installed hashes,
the wheel hash and unchanged crystal/sidecar hashes;
``verification_release/import_audit.json`` records the actual imports.
The supplied driver's field named ``native_module_sha256`` identifies the
Python ensemble facade; the import audit separately identifies the compiled
extension.

Verified artifacts are exported to the fresh diagnostics directory::

    F:/fuchs-4355-a/processed/LaNiO3/CTR/data_out/1V2/0056/diagnostics/performance_trial_20261107_batch_reuse

Reproduce the paired benchmark from the repository root::

    & C:/Users/timof/miniforge3/python.exe benchmarks/benchmark_dwba_c1_paired.py `
      --request-dir F:/fuchs-4355-a/processed/LaNiO3/CTR/results/0056/dwba_ensemble/installed_api_verification `
      --baseline-native output/dwba_reuse/baseline/_CTRcalc_cpp.cp314-win_amd64.pyd `
      --baseline-ensemble output/dwba_reuse/baseline/_CTRincoherent_dwba.py `
      --output output/dwba_reuse/paired_new

Choose a fresh output directory. ``benchmark_dwba_c1.py`` additionally covers
unchanged, geometry, DWF, occupancy and roughness proposals outside batches,
both optical policies, and eight-job runs with 1/2/4/8 workers. Its JSON records
imported source paths/hashes, setup versus evaluation time and peak RSS.
Use ``--cases shared --inner-repeats 7`` with the paired driver to reproduce
the longer shared-field comparison.

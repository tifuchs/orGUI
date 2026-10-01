# Rocking batch integration implementation plan

Date: 2026-09-30.

Status: implemented and regression-validated in the working tree.
This document retains the reviewed design and records the implementation below.

The parent [multiple-line integration plan](multiple_line_scan_integration_plan.md)
contains the approved stationary in-memory approach and the rocking memory
assessment. Stationary results remain in RAM. This plan applies bounded
persistence to rocking extraction and bounded reads to its later reduction.

## Outcome and constraints

Integrate all selected reciprocal-space lines in one pass through source
images, without allocations proportional to all frames times all sampled
points. Retain complete rocking curves on disk for interactive peak/background
selection and re-reduction. Keep the existing ROI kernels, geometry formulas,
correction branches, float64 precision, error propagation, and per-line
NeXus result format.

No detector image stack is retained. Workers perform image reads/integration;
a single collecting thread owns output writes. Source backends declaring
`supports_concurrent_read = False` must be read serially.

## Integration module boundary

Create `orgui/app/scan_integration.py` for both stationary and rocking
integration orchestration. This is an explicit user requirement: reduce
`orGUI.py` while adding multiple-line support and rocking batches.

- Move the large stationary integration body, rocking worker/collector logic,
  curve assembly and shared integration helpers into this module. Proposed
  public entry points are `integrate_stationary_scan()` and
  `integrate_rocking_scan()`; finalize signatures during extraction.
- Keep `orGUI.integrateROI()`, `rocking_extraction()`,
  `rocking_Bragg_extraction()`, `rocking_static_extraction(xy, hsize, vsize)`
  and `rocking_integrate(xylist, rois, hkl_del_gam, refldict, name)` as thin compatible
  workflow wrappers. They read widget state, prepare immutable options/config
  snapshots and line definitions, dispatch, and apply plotting/UI updates.
  Preserve scripting entry points, arguments and returned status dictionaries.
- Give the new module explicit prepared input/context objects and narrowly
  scoped callbacks for image access, geometry, progress, persistence and
  sampled plot output. Do not pass the entire main window as a replacement
  for a defined interface or scatter direct widget access throughout it.
  App-level context may contain configuration snapshots and scan access;
  numerical `datautils` APIs must continue receiving numbers/arrays only.
- Numerical geometry remains in the existing geometry modules; widget-driven
  preparation stays in the wrappers. Move integration-specific numerical
  trajectory assembly out of `orGUI.py` where it can accept explicit values
  and call the established geometry APIs, including bounded curve geometry.
- Keep correction definitions in `datautils/xrayutils/corrections`, app
  policy adaptation in `integration_corrections.py`, HDF5 lifecycle helpers
  in `database.py`, and interactive angular reduction in `peak1Dintegr.py`.
  The new file coordinates those components rather than duplicating them.
- Shared integration functions must be safe for CLI/headless use. Workers
  must not access Qt widgets or write output HDF5. The collecting path sends
  progress/results through the established main-thread workflow and handles
  logging/cancellation consistently.
- Add Sphinx-friendly docstrings for every public function, including units
  and array/index conventions at relevant boundaries. Avoid cyclic imports:
  the new module must not import the `orGUI` main-window class.

First extract the existing single-line behavior behind these boundaries and
verify equivalence, then generalize routing and batching. Do not combine this
move with unrelated formatting, variable renaming or physics changes. Keep
focused tests of the new module independent of a live main window where
possible, plus wrapper tests for existing dispatch and scripting behavior.

### Required compatibility adapter for geometry

The current public `rocking_integrate()` requires an already materialized
`hkl_del_gam` argument. Preserve that signature and honor caller-supplied
geometry by taking views/slices; do not recompute or discard it. Calls using
this legacy API retain the RAM cost the caller has already incurred.

The built-in CTR, Bragg and fixed-pixel extraction wrappers must instead
dispatch directly to the new prepared-request API with a lazy geometry
provider. They must stop calling `getStaticROIparams(xy)` for all ROIs before
dispatch, or batching would still allocate the complete `(R, F, 6)` array.
Keep `getStaticROIparams()` itself compatible for scripts. Its actual six
columns are `(h, k, l, delta, gamma, scan_axis)`; its current docstring's
mention of pixel coordinates/mask flags does not describe the returned
array. Pixel centers and detector-valid flags are separate prepared inputs.

## API review findings and resolved details

Reviewed against the working-tree APIs on 2026-09-30. These are implementation
requirements beyond changing the worker loop:

| Current API or behavior | Required plan detail |
| --- | --- |
| `get_rocking_scan_info()` validates one constant H_0/H_1 and assumes `h5_obj.parent.parent` is the scan | Publish one ordinary `/scan/measurement/<line_name>` rocking group per line. Do not put mixed lines into a single `rois` group or insert an extra published job-container level. |
| `rocking_static_extraction()` and Bragg extraction also call `rocking_integrate()` | All three workflows use the new engine. Line vectors remain optional for Bragg/fixed-pixel batches; preserve their existing synthetic `s`, peak metadata and ROI-size behavior. |
| `intkeys_rocking()` derives automatic sizes from neighboring coordinates and reads GUI defaults | Resolve size/projection/offset settings for each line before flattening. Explicitly pass its intersection and retain the same validity ordering. Never infer neighboring spacing across different lines. |
| `ROIState` and option round-trip currently cover only region/advanced/rocking_scan dictionaries | Line collections and batch options need an explicit additive serialization/restore path, not merely new keys in `get_integration_options()`. Absence must restore existing single-line behavior. |
| `_prepareFootprintAction()` loads `record.alpha` to build replacement illumination | Capture the live beam/profile action once; evaluate replacement illumination over each curve tile rather than materializing the full 2-D record. |
| Acceptance and peak-arm helpers also load full arrays | Slice their inputs or add private row-aware adapters; changing only the curve reader leaves an unbounded reduction path. |
| The core reducer indexes ROI bounds and correction arrays with local row numbers | Slice every row-dependent input, preserve full-line anchor interpolation, offset progress, and merge all nested result fields. |
| Rocking read exceptions abort extraction before saving | Abort the batch job on the first read/integration failure. Do not introduce stationary-style skipped-image processing or publish fill values as a completed rocking result. |
| `add_nxdict()` flushes and synchronizes the GUI tree; no batch writer exists yet | Add lifecycle/slice-write helpers and use tree synchronization only at publication. These are proposed new APIs, not existing capabilities. |
| Deleting HDF5 datasets does not guarantee file shrinkage | Same-file scratch can leave unused allocated disk space after cleanup; this is an explicit storage tradeoff, described below. |

## Batch definitions and routing

Use two batch dimensions for different purposes:

1. **Frame batch:** up to `B` source frames, integrating all valid sampled ROIs
   from all selected lines. Worker output is raw/correction/background
   counters, indexed by frame and a stable global ROI slot.
2. **Curve batch:** up to `K` ROI rows, each containing all `F` source frames.
   Used for geometry, curve construction, correction/provenance, output
   writing, and later angular reduction. Never split an angular integral at a
   frame-batch boundary.

Prepare an ordered line table with `line_id`, `H_0`, `H_1`, sampling settings,
selected intersection, and the existing config snapshot. Keep existing
single-line calls as an adapter. Current rocking intersection selection stays
unchanged; handling both intersections is not required to implement batches.

For each line, prepare coordinates/ROI rectangles using
`get_rocking_coordinates()` and `intbkgkeys_rocking()`. Concatenate valid
sampled points into the kernel's existing `(R, 2, 2)` ROI inputs. Store an
explicit mapping from each global slot to line id, original sampled-point
index, intersection and `s`. Retain per-line start/stop offsets and the
original detector-valid mask. Empty lines are reported and skipped without
shifting other identities. Completion order must never determine a slot.

Keep the first implementation's ROI/background/advanced settings and selected
intersection shared across requested CTR lines, matching existing controls.
Resolve the effective sampling for each line using the existing detector
resolution rule without repeatedly mutating GUI controls. Persist each line's
effective `delta_s`, not just the initial requested value. Bragg/fixed-pixel
preparation uses its current override rules and does not require H_0/H_1.

Add optional ordered line-definition fields to `ROIState` and its dedicated
NeXus serializers, config capture/restore, and selector APIs; version that
additive ROI configuration change. Existing `region`, `advanced`,
`rocking_scan` and `LegacyKeyDict` aliases keep their meanings. Batch sizes
and memory budget are execution options saved as job provenance, not new
scientific correction fields or changes to `ctr_curve_v3` algorithm labels.

## Storage choice and lifecycle

Use temporary chunked datasets inside the active database, under a reserved
top-level work group such as `/_orgui_integration_work/<job_id>`. This avoids
introducing an external scratch-file lifecycle and reuses the selected HDF5
compression filters. Dataset names and lifecycle attributes in this namespace
are internal, not a replacement for the public curve schema.

Store counters as `(ROI, frame)` datasets, matching the existing curve reader
orientation. Raw center/background sums, pixel counts, correction sums/counts
are required; background-image counters are created only when that image is
used. Final legacy fields that currently contain zeros without a background
image must still have their existing meanings and presence. Do not apply
normalization or illumination prematurely to temporary counters.

Proposed initial chunk shape is `(min(R, 32), min(F, 64))`, tuned by measured
write/read behavior and the current compression filter requirements. Select
chunk shape per dataset where dimensions differ. Cap HDF5 chunk-cache use;
do not allocate a large independent cache for every dataset/line. Write
contiguous bounded tiles rather than all ROI rows in a single transposed
copy. File compression affects I/O and disk space, not batch sizing in RAM.

The job records config/source identity, line mapping, total frame count,
batch sizes, committed frame ranges, failed-frame information, and state:
`extracting`, `assembling`, `complete`, `cancelled` or `failed`. Mark a frame
batch committed only after its counters and frame-status data are written and
successfully flushed. Completed ranges are diagnostic in version one;
automatic resume is deferred until input/config identity validation exists.

Build final per-line rocking groups within staging. Publish each fully built
line at its ordinary measurement path using an HDF5 group move after a
successful flush; retain the current unique-name/suffix behavior. Update scan
defaults and the database tree only for published complete groups. A move is
a metadata operation, not a claim of a crash-safe multi-line transaction.
If publication fails after some lines complete, report those complete lines
and the remaining failure explicitly; never overwrite or delete earlier
measurements. Internal work groups must be excluded from normal result
discovery and reduction.

Delete temporary counters after successful publication. On handled
cancellation/failure, stop work and clean only this job's unpublished work
when the file is usable, leaving prior results intact. If cleanup or flush
fails, report the remaining job path and preserve incomplete status where
possible. On reopening, identify abandoned work groups and offer explicit
cleanup through a GUI-only action; do not present them as valid measurements.

Same-file temporary storage is retained as the initial implementation choice,
but deleting it means removing dataset links, not guaranteeing reclaimed
filesystem space. HDF5 may retain free space for reuse and the physical file
may stay larger than the final result. Log this distinction; do not promise
that cleanup shrinks the database or perform automatic repacking. If bounded
final physical file size becomes a requirement, change the counters to a
separate temporary HDF5 file with explicit cleanup; that is a storage-policy
decision independent of scientific batching. A local h5py allocation/deletion
check during this review confirmed unchanged file size after deletion,
flush and close when later allocations prevented tail truncation.

Use existing `@NX_class`, `@orgui_meta`, `@unit` and related dictionary metadata
through the established NeXus serializer for small metadata groups; raw h5py
attributes use the names without `@`. Normalize `database.compression` using
the same supported filter semantics as the current `dicttonx` path, including
None, built-in strings and plugin filters, and test the direct dataset path.
Use a reserved work ancestor to suppress all unfinished results in context
menus/discovery, even if staged descendants have their eventual rocking
metadata. The tree does not currently hide underscore-prefixed groups by
contract. Reuse `_requireOpenFile()`, flush failures and `closeSafe()` behavior.

Before publishing a line, ensure the real scan parent contains its existing
`instrument/positioners`, auxiliary counters, title and scan config; ensure
the line group contains integration config. Freeze the source scan/geometry
objects and database identity for the job, and prevent GUI close/replacement
or scientific state changes from invalidating them during event processing.
Final paths keep `orgui_meta = rocking`, the existing `rois`/curve-record
groups, `H_0/H_1` shape `(P, 3)` when present and `s` shape `(P,)`, even for
one-point lines. Preserve separate peak-vector and per-frame field shapes.
Write job/line identity as additive metadata, retaining existing result names
for single-line calls and suffixing collisions for multiple definitions.

## Stage 1: extract bounded frame batches

1. Freeze line/ROI/config state before launching workers. Share the existing
   detector mask, correction maps, background image, and repair settings.
2. Preallocate temporary counter datasets; do not call
   `getStaticROIparams()` for every ROI/frame pair at this stage. Raw image
   integration requires prepared static ROI rectangles, not the complete
   reciprocal-space trajectory array.
3. Submit at most the bounded in-flight limit, initially the worker count.
   Never queue all scan frames or retain completed futures scan-wide.
4. Read each source frame once and call the existing accelerator with all
   prepared ROI slots. Generalize the NumPy fallback to the same mapping.
5. Copy returned counters to that frame's location in the current batch,
   record success/error status, and release the future and its result.
   A read/integration exception cancels outstanding work and aborts the job,
   matching current rocking behavior. Persist diagnostics when possible but
   do not publish any extraction output from a failed frame pass.
6. Write and flush the completed batch, then release/reuse the buffer before
   advancing. Report progress by completed frames, not submission indices.

Batch storage is `O(B * R)`, worker image memory is `O(W * detector_pixels)`,
and prepared ROI metadata is `O(R)`. File reads still scale with `F`, not
`F * lines`. No output HDF5 API may be called from an image worker.

## Stage 2: assemble bounded groups of complete curves

1. Visit each line and read at most `K` temporary ROI rows across all frames.
2. Compute the corresponding geometry as `(K, F, 6)` using a localized
   extension/private helper of `getStaticROIparams()`. Preserve scanned-arm,
   refraction, angle signs and units. Shared per-frame axes may remain in RAM.
3. Run existing background subtraction, masked-area scaling, polarization
   and uncertainty calculations over these complete curves. Preserve
   whole-curve decisions such as `np.any(bgpixel)` and alternative
   background-image branches; a different frame batch must not select a
   different algorithm for part of a curve.
4. Build the existing diagnostic and CTR photon branches plus versioned base
   variance and pixel factors. Keep all current normalization/illumination
   application points and algorithm/scale labels. Continue calling physics
   in `datautils`; do not duplicate its formulas in the writer.
5. Write slices directly into the final per-line `rois` and `ctr_curve_v3`
   datasets in staging. Preserve units, shapes, optional fields, trajectory
   and config snapshots. Framewise policy computation must use correctly
   sliced per-curve alpha and the full source frame metadata; do not pass a
   smaller frame count to helpers that index full scan counters.
6. Release all tile data, including legacy dictionary assembly temporaries.
   Remove the scan-wide per-ROI dictionaries and `np.vstack` pass from the
   new path. Keep plotting bounded using the existing curve sampling policy.
7. Verify dataset shapes, completion masks and record availability before
   publishing the line. Release its temporary counters when safe.

Do not remove compatibility branches or replace repeated saved arrays with
new link/broadcast conventions in this change. Such disk-format optimization
needs a separate audit of existing readers and scripts.

## Stage 3: keep later angular reduction bounded

`RockingPeakIntegrator.integrate()` currently calls `get_all_ro_curves()` and
can materialize all source curves. Bounded extraction alone is insufficient
if reduction restores the same memory scaling.

Add a bounded curve-row reader/iterator that handles both legacy `rois` and
versioned records. Slice every two-dimensional correction field with the
same ROI range; one-dimensional per-frame divisors stay shared. Support old
files with no line table or job-state attributes. Keep `get_ro_curve(idx)`
and explicit `get_all_ro_curves()` scripting behavior compatible, while the
normal reduction path uses bounded reads.

Run `_compute_rocking_integration()` on complete curves in each tile, with
the corresponding peak/background bounds, anchor/interpolation selections,
monitor factors, exact stored Q/H policy and selected footprint action.
Resolve curve-dependent selections using the full line's scalar metadata
before tiling so interpolation does not reset at tile boundaries. Apply
normalization/illumination at the same stage and exactly once. Retain the
current degree-to-radian conversion and Lorentz/acceptance conventions.

Combine only the reduced per-ROI scalar results (or write them in bounded
tiles if needed), then save using existing result/error handling. Check all
remaining GUI paths for accidental full materialization of curves or
two-dimensional alpha/illumination arrays.

Specifically, the tile adapter must provide `s_array[K]`, curve/error arrays
`(K, F)`, `roi_info[roi][from/to][K]`, tiled Lorentz/rod/illumination factors,
and acceptance/solid-angle factors `[K]`. Keep source `axis[F]` and auxiliary
counters `[F]` shared. Preserve scalar-factor support and legacy x/y/size
fallback shapes. Adapt `_rocking_peak_arm_angles()`, `_rocking_acceptance()`
and `_rocking_solid_angle_mean()` to tile data while retaining their unit
attributes, calibrated-arm fallback and warning behavior. This applies to
the reversed-axis and mu/th peak-selection rules as well.

Read `_compute_rocking_integration()`'s actual nested return schema when
merging: aggregate scalar result arrays and every `int_data` ROI field,
including nested auxiliary counters. Do not merge only `croibg`/`F2_hkl`.
Progress callbacks need a tile row offset. Its current cancellation callback
breaks the row loop before aggregation and can leave incompatible lengths;
the new driver must use an explicit cancellation signal that unwinds before
aggregation/publication, or a focused core fix with a regression test. Keep
the public single-call helper behavior compatible wherever possible. On
reduction cancellation, do not save partially merged scalar results.

Preserve total-flux validation, separate diagnostic/CTR branches, stored
algorithm dispatch, and the calibrated absolute prefactor applied after the
core in `integrate()`. Keep the footprint action validation identical. The
new pipeline may capture scalar policy globally and compute numerical
factors per tile, but it must not weaken unavailable-provenance checks.

## Memory budget and initial settings

Start with `B = 64` frames and `K = 32` curves, reducing either size to meet a
dedicated working-memory budget. These are performance defaults, not fixed
scientific settings; changing them must not change results. Do not repurpose
`maxMemory` without auditing its existing MB convention and other consumers.

For 10,000 frames and 100 rods with 61 valid points each (`R = 6100`):

- A 64-frame buffer of twelve float64 counters is 35.7 MiB.
- A 32-curve complete-frame tile costs 195.3 MiB at a conservative 80-field
  allocation estimate; begin with at least 256 MiB reserved for this stage
  and measure real correction/geometry allocations.
- Six full geometry fields within that tile occupy 14.6 MiB; include this in
  the tile estimate rather than counting it again.
- Detector images, conversion/loading copies, correction maps, background,
  masks, repair buffers, GUI data and bounded HDF5 caches remain separate.
  Reduce worker count as needed; respect backend concurrency constraints.

Log selected `B`, `K`, workers, ROI count and estimated working allocations.
Reject settings that cannot accommodate one worker/frame/curve after shared
buffers are accounted for. Stage 1 and stage 2 run sequentially, so do not
reserve both full working peaks simultaneously. Temporary raw counters alone
require approximately 3.64 GiB without a background image or 5.45 GiB with
one, uncompressed; disk planning must include final outputs and temporary
overlap. Do not estimate free disk requirements from the reported 1 GB
compressed file size.

## Implementation order and affected files

1. **Pin current behavior:** focused baselines for single-line extraction,
   failure fills, background branches, geometry and reduction. Review
   existing local changes in `peak1Dintegr.py` and its tests; preserve them.
2. **Extract the integration module:** move existing single-line orchestration
   to `orgui/app/scan_integration.py`, define prepared inputs/callbacks and
   keep compatible thin `orGUI.py` wrappers. Verify stationary and rocking
   equivalence before batching.
3. **Define routing and batch options:** line/ROI mapping and batch planning
   in `scan_integration.py`, workflow preparation in `orGUI.py`, minimal
   config serialization in `config_data.py`, and multiple-line controls in
   `QScanSelector.py`.
4. **Add staging writer:** focused helpers in `database.py` for chunked
   allocation, bounded slice writes, flush/error handling, result publication
   and cleanup. Reuse `_requireOpenFile()` and existing compression settings;
   avoid GUI tree refresh or `add_nxdict()` after every numeric slice.
5. **Batch extraction:** replace scan-wide counters/futures with frame batches
   in `scan_integration.py`, using existing kernels and fallback semantics.
6. **Batch assembly:** compute geometry and corrected curves per ROI tile,
   write the existing output schema directly, and remove the duplicate
   scan-wide assembly representation.
7. **Batch reduction:** update normal `peak1Dintegr.py` reduction to use the
   bounded reader for both legacy and versioned inputs.
8. **Finish workflow:** progress/cancellation and incomplete-job discovery,
   compatible line metadata/names, bounded plotting, and shared CLI errors
   through `logger_utils`. Do not introduce modal dialogs in shared paths.
9. **Validate and document:** targeted regressions, representative peak-RAM
   measurements, `doc/source/image_integration.rst`, and `CHANGELOG.md`.
   Regenerate release notes through `doc/generate_release_notes.py`.

## Acceptance checks

- Single-line batched output and angular reduction agree with the baseline
  for raw/corrected signals, uncertainties, trajectories and stored factors.
  Exact sums/counts should match where the same kernel is used; establish
  justified numerical tolerances for geometry/correction/reduction results.
- Multi-line extraction equals separate per-line extraction with identical
  settings, including different valid-point counts, empty lines, and stable
  original point identities after filtering.
- Results are invariant under `B`/`K` choices, final short batches and
  deliberately shuffled worker completion. Include a batch size of one.
- Exercise accelerator/fallback, masks, eligible repair, supported fitted
  background behavior, background images, changing detector arms,
  normalization and legacy/explicit illumination policies.
- Later reduction agrees for peak/background windows across tile boundaries,
  interpolation anchors, stored/replaced footprint actions and re-reduction.
- Source-image read count is one per frame during extraction (apart from
  existing setup/preview reads); outstanding futures and buffered frames
  never exceed their limits.
- Cancellation and injected read/write/flush failures return explicit status,
  leave existing measurements usable, and never publish partial line groups.
  Reopening an abandoned staging group must not expose it to the reducer.
- Old databases and single-line configuration files still load/reduce.
  New per-line result groups are readable through existing public layouts.
- Cover existing direct `rocking_integrate()` calls with deliberately supplied
  geometry, and the CTR/Bragg/fixed-pixel wrappers with lazy geometry. Test
  singleton/empty lines, duplicate names, S2 selection, per-line automatic
  sizing and line config capture/restore. Verify correct scan-parent paths.
- Exercise replacement-footprint, acceptance and peak-arm paths with a dataset
  spy that rejects unbounded 2-D reads. Cover all nested reduction outputs,
  progress offsets and cancellation inside a tile, not just between tiles.
- Check built-in and available plugin compression paths, staging exclusions
  in context menus, and physical/logical disk-space reporting after cleanup.
- Instrument maximum allocated buffer shapes and measure process peak RAM
  at fixed `B`/`K`/workers while increasing `F`/`R`. RAM must follow bounded
  batch terms plus line/frame metadata, rather than a full `F * R` payload.
  Use generated small fixtures for routine tests and a separate scaling run
  for the 10,000-frame/100-rod case.

Run the narrowest relevant app tests first, then `pytest orgui/app/test` and
`ruff check orgui` as warranted. If geometry helpers change, also run
`test_HKLcalc.py` and `test_DetectorCalibration.py`. Fix local findings without
scientific or lint-only refactoring.

## Implementation and validation record

Implemented on 2026-09-30 in the working tree:

- `scan_integration.py` owns stationary extraction/assembly, rocking workers,
  frame batches and curve tiles. `orGUI.py` retains compatible entry points
  and captures numerical state before workers/progress events. Its net size
  is approximately 2,480 lines smaller.
- Ordered stationary/rocking line tables are editable in the selector and
  persisted by additive ROI schema version 2. Older config restores select
  the single-line controls. Fixed and Bragg dispatch remain unchanged.
- `RockingBatchWriter` uses same-file chunked scratch, payload flush before
  committed-frame markers, per-line publication, compression-filter reuse,
  job-specific cleanup and reader/context-menu exclusions. Published groups
  record batch sizes/budget, IDs, original point indices and sampling step.
- `rocking_tiles.py` supplies bounded HDF5 row views and nested scalar-result
  assembly. Later reduction tiles ROI bounds, correction records, live
  illumination replacement and acceptance/arm inputs. Cancellation raises
  before partial aggregation/save. The existing numerical core is unchanged
  relative to the user's pre-existing physics edits.
- Frame budget estimation includes both the collector buffer and in-flight
  results (192 bytes/frame/ROI). Curve estimation reserves 1,280 bytes per
  curve/frame for geometry, counters, temporary legacy dictionaries, stacking
  and corrections. This is more conservative than the provisional 80-field
  estimate above: at F=10,000 and 256 MiB the default K=32 is capped to 20.
  Detector images/maps, file caches and application state remain additional.

Synthetic scaling experiment (`tmp/measure_rocking_batches.py`, 64 MiB budget,
B=64, K=32, two workers, gzip, 1,000 frames, streamed 12x12 images):

| Curves | Traced Python/NumPy peak | Process RSS growth | Result |
| ---: | ---: | ---: | --- |
| 500 | 36.9 MiB | 96.2 MiB | complete |
| 1,000 | 19.6 MiB | 87.7 MiB | complete |

The first run includes import/allocator/cache startup; these numbers are not
a detector-sized production RAM guarantee. Doubling the stored curve payload
does not double working RAM at fixed tiles. This validates bounded assembly
on generated input; a real 10,000-frame/100-line acquisition remains an
explicit production-scale validation task, particularly for detector image
copies, repair and illumination options. No representative raw scan was
provided. Same-file cleanup still does not promise physical file shrinkage;
automatic resume/repack remain deferred.

Validation: the app suite plus HKL/detector geometry checks passed (581 tests,
one skip, 239 subtests). Final focused runs cover additional row-lazy footprint
replacement, plugin compression/setup failure, GUI reentry/control restoration
and HDF5 line-table round trips. `ruff check orgui` passes. The subprocess
Numba check requires PYTHONPATH pointing to this checkout; otherwise it loads
an older installed copy outside the writable workspace.

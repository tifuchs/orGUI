# Multiple-line scan integration plan

Date: 2026-09-30.

Status: stationary in-memory architecture and rocking batch processing implemented
in the working tree. The concrete rocking implementation sequence is in
[`rocking_batch_integration_implementation_plan.md`](rocking_batch_integration_implementation_plan.md).
The implementation retains the existing correction formulas and per-line
results. Validation and measured synthetic memory are recorded in the rocking plan.

## Scope and scientific invariants

Extend stationary integration from one reciprocal-space line to an ordered
collection of lines. Each line is defined by `H_0 + s * H_1`, with `H_0` and
`H_1` in r.l.u. Retain the existing geometry, S1/S2 Ewald intersections,
detector coordinate conventions, ROI/background rules, corrections,
uncertainty propagation, and configuration snapshots.

Keep the current architecture and compiled image integration functions. This
is a routing and workflow extension, not a physics or storage-format rewrite.
Move stationary and rocking integration orchestration into the new app module
`orgui/app/scan_integration.py` as part of the implementation. Keep thin
compatible workflow wrappers in `orGUI.py` and retain scientific definitions
in their existing geometry/correction modules.

Multiple reciprocal-space line selection applies to the stationary `hklscan`
and rocking `hklscan` workflows (current tab ids 0 and 2). Preserve stationary
`fixed` ROI behavior (tab 1) and rocking Bragg behavior (tab 3); their adapters
may use the shared engine without requiring line vectors. Capture the mode
explicitly before dispatch: `getROIloc()` currently branches on the live tab
even when H_0/H_1 are supplied, so explicit vectors alone do not select the
stationary line workflow.

## Agreed stationary decisions

- Stream detector images; keep extracted counters and trajectory metadata in
  RAM until extraction completes. Do not write each frame to the database.
- Read each source frame once for integration and process every requested line
  and intersection in that call. Share detector-sized correction maps, masks,
  and any background image across lines.
- Assign each line a stable identity and preserve its definition in the saved
  metadata. Route results by explicit frame, line, and intersection indices,
  independently of worker completion order.
- Use array storage conceptually indexed as
  `(frame, line, intersection, counter)`. The kernel can continue receiving
  flattened ROI arrays; use a deterministic slot mapping such as
  `slot = 2 * line_index + intersection_index` for S1/S2.
- Keep invalid intersections in stable slots with explicit validity flags.
  Do not compact them without retaining an explicit reverse mapping.
- Bound work submission and release completed futures after copying their
  counters. Avoid retaining duplicate results for the entire scan.
- After extraction, perform correction/output assembly and database saving
  one line at a time to limit temporary allocations. Preserve the existing
  per-trajectory datasets and correction records where possible.
- Existing single-line controls/configurations remain usable and retain their
  previous result naming and numerical behavior. The multiple-line workflow
  should default to the existing single line.
- Preserve `integrateROI()`'s current mode dispatch and the fixed-pixel
  stationary workflow. Multiple-line definitions must not silently override
  those modes.

### Implementation sequence

1. Extract integration orchestration into `orgui/app/scan_integration.py`
   behind thin existing `orGUI.py` entry points, and introduce an ordered
   collection of line definitions at the app/workflow boundary, with a
   single-line adapter for existing callers. Keep GUI and config objects out
   of `datautils`.
2. Prepare each line's S1/S2 trajectories using the existing geometry helpers;
   assemble validity and ROI arrays using the explicit slot mapping.
3. Generalize stationary worker output allocation and result distribution to
   all slots, including the NumPy fallback. Retain the current compiled
   functions, which already accept an arbitrary number of ROIs.
4. Generalize polarization, background, correction, and uncertainty stages
   over lines/intersections without changing their formulas.
5. Add line selection/editing and per-line output/provenance using the
   existing workflow. Decide the minimal configuration representation during
   implementation; preserve old configuration/calibration compatibility.
6. Verify single-line equivalence, distinct multi-line signals, S1/S2
   routing, invalid intersections, out-of-order completion, masking,
   background images, and both accelerator/fallback paths. Retain existing
   cancellation and database-error handling.
7. Measure peak RAM on a representative scan, then update user documentation
   and the changelog for the implemented behavior.

### Stationary memory estimate

Conservative case: 10,000 frames, 100 lines, both intersections retained,
64-bit numeric arrays. This is 2,000,000 frame/line/intersection records.
The following extends current allocations rather than assuming a minimal
intensity-only output.

| Allocation | Fields per record | RAM |
| --- | ---: | ---: |
| Geometry, pixel coordinates, validity | 9 | 137.3 MiB |
| Center and four background rectangles, four int64 bounds each | 20 | 305.2 MiB |
| Raw, correction-map, background-image counters | 12 | 183.1 MiB |
| ROI sizes and integration boundaries | 6 | 91.6 MiB |
| Polarization counters | 4 | 61.0 MiB |
| Subtotal | 51 | 778.2 MiB |

Allow roughly 1-2 GiB for scan-wide results/geometry including derived
quantities and temporary arrays. Background counters are included even when
unused. This is an allocation estimate, not a measured implementation peak.

Image memory scales with active workers and detector pixels, not line count:
`workers * pixels * 8` bytes for one float64 image per worker. At 16 workers
this is 0.30 GiB for the example 1475 x 1679 detector, 0.50 GiB for 2048 x
2048, or 0.95 GiB for 8 megapixels. Loading/conversion copies, detector caches,
maps, masks, background and optional repair buffers need additional headroom.
About 2-4 GiB incremental RAM is a planning allowance for a streamed 2-8
megapixel case, beyond the existing app footprint. Available RAM and options
must be checked when validating it.

Do not retain the full detector stack: 10,000 float64 images would require
184.5 GiB at the example detector size or 596.0 GiB at 8 megapixels.

## Rocking scan assessment

### Why line count alone is insufficient

Stationary integration produces up to two ROI observations per line per
frame. Rocking extraction samples many fixed detector ROIs along each line
and retains an entire framewise rocking curve for every sampled point.
The final angular reduction in `peak1Dintegr.py` happens later; the extraction
must preserve the curves for peak/background selection and reprocessing.

Let `F` be frames, `L` lines, and `P_j` the number of valid sampled detector
points on line `j` for the selected intersection. The relevant size is
`F * sum(P_j)`, or `F * L * P` for equally sampled lines. Current rocking
extraction selects one intersection. Processing both would add their valid
point counts and can approximately double memory.

`get_rocking_coordinates()` requests
`round(maxValue / step_width) + 1` points before detector-valid filtering.
The controls initially set max s = 6 and delta s = 0.1, giving 61 requested
points. Actual settings, the automatic step adjustment, and detector coverage
determine the valid count; 61 is an illustration, not a universal default
after scan/geometry changes.

### Current allocation behavior

In `orGUI.rocking_integrate()`:

- `getStaticROIparams()` retains six float64 geometry fields for every
  curve/frame pair.
- Twelve float64 `(frame, ROI)` counter arrays are allocated, including
  background-image arrays even when no background image is used.
- All submitted futures remain in a dictionary, retaining another eight
  counter fields per pair, or twelve with a background image.
- Per-ROI legacy dictionaries retain derived curves/factors and angle arrays
  while the final structured result is built. They reference some existing
  arrays and own others; they are not all independent copies.
- The `np.vstack` conversion creates full arrays for the output while the
  original arrays/dictionaries remain live. For a mu/omega scan this is
  approximately 27 full curve/frame fields, or 29 with background-image
  alternative curves. Chi/phi and peak metadata can be smaller per-ROI
  fields rather than full frame arrays.
- The versioned correction record shares several structured arrays in RAM
  but adds base variance and radian incidence arrays. Enabled framewise
  illumination may add more full arrays. Shared arrays must not be counted
  twice just because they occur at multiple dictionary paths.
- Only `data_2d_structured` is saved. The transient per-ROI legacy dictionary
  is used to assemble it, so it increases RAM without itself becoming a
  second complete saved representation.

Extraction counters plus geometry alone need `18 * 8 = 144` bytes per
curve/frame pair. A planning estimate for the current assembly peak is
roughly 65-80 float64-sized fields per pair (520-640 bytes), including retained
futures, per-ROI derived arrays, structured copies, and correction temporaries.
This range is inferred from allocations, not a profiled peak or a guaranteed
upper bound. Detector buffers, existing GUI data, allocator overhead, and
plotting are additional. Worker images normally expire before output
assembly, so their full peak need not coincide with the assembly peak.

For F = 10,000 and L = 100, one selected intersection:

| Valid points per line | Curve/frame pairs | Counters + geometry | Estimated current assembly peak |
| ---: | ---: | ---: | ---: |
| 1 | 1 million | 0.13 GiB | 0.48-0.60 GiB |
| 50 | 50 million | 6.71 GiB | 24.2-29.8 GiB |
| 61 | 61 million | 8.18 GiB | 29.5-36.4 GiB |
| 100 | 100 million | 13.41 GiB | 48.4-59.6 GiB |
| 500 | 500 million | 67.06 GiB | 242-298 GiB |

Thus simultaneous in-memory extension of the current rocking architecture
to 100 rods is not a safe general target for 16-32 GiB machines. Even 64 GiB
can be marginal around 100 valid points per line once app memory and enabled
options are included. Fewer frames or fewer valid points reduce these
estimates linearly.

### Compressed database size

A reported database size of about 1 GB does not specify decompressed array
size, the number of scans it contains, or how much is loaded simultaneously.
Repeated frame axes, metadata and counters can compress strongly. RAM also
includes working arrays and copies absent from disk. Conversely, loading a
database does not necessarily materialize all its arrays: the current
versioned rocking reader preserves 2-D `h5py.Dataset` objects and supports
row selection for an individual curve.

For a representative database, sum `dataset.size * dataset.dtype.itemsize`
for the numeric arrays in the relevant scan/rocking groups and compare it
with `dataset.id.get_storage_size()`. Account separately for aliases/hard
links and multiple saved correction branches, then profile the extraction
and output-assembly phases. No representative 1 GB rocking database was
provided or benchmarked for this assessment. On-disk size alone must not
drive a RAM policy.

### Agreed rocking direction

Keep the existing per-image integration kernel and preserved scientific curve
format. Use bounded frame batches for extraction into temporary chunked
datasets, followed by bounded groups of complete curves for correction and
final output. HDF5 writes belong to the collecting thread, not image workers.
Read each detector image once during extraction, with explicit
`(line_id, point_index, intersection)` routing across all selected lines.

Keep full source curves for later peak/background selection and reprocessing;
batch the later angular reduction as well. Preserve completion/cancellation
semantics and old single-line files. The linked implementation plan specifies
temporary storage, publication, memory limits, compatibility and validation.

## Code references at assessment time

- `orgui/app/orGUI.py`: `integrateROI`, `getROIloc`,
  `get_rocking_coordinates`, `getStaticROIparams`, `rocking_integrate`,
  `_rocking_arm_snapshot`.
- `orgui/app/cpp/roi_sum_cpp.cpp`: `processImage_Carr`,
  `processImage_bg_Carr` and related repair/fitted-background variants;
  ROI input is `(N, 2, 2)`, counter output is `(N, 4)`.
- `orgui/app/peak1Dintegr.py`: `_storedCurveCorrectionRecord`,
  `get_all_ro_curves` and later angular reduction.
- `orgui/app/config_data.py`: `curve_correction_record_to_nxdict` and
  versioned curve records.

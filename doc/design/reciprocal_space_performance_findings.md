# Reciprocal-space mapping: performance findings and open items

Companion to `reciprocal_space_scratch_architecture.md`, which describes what
the pipeline *is*. This one records what was measured, what those measurements
settled, and what is still open. Two work streams have grown out of it into
plans of their own: `reciprocal_space_mapping_locality_and_subdivision.md`
(complete) and `reciprocal_space_mapping_serial_fraction.md` (next).

> **Read `reciprocal_space_mapping_locality_and_subdivision.md` alongside
> this.** Since these numbers were taken, the subdivision was rewritten
> (corner-passing, 1.65-2.98x at depth >= 1), the memory prechecks were
> re-bounded, and mapping now maps several adjacent frames in one native
> call by default. The timings below are the **per-frame** pipeline and
> still describe it; where a grouped run differs, it is noted inline. That
> document also carries the measurement discipline this one's "Measurement
> traps" section started, and it is the more complete of the two on that
> subject.

## How to read the numbers here

Everything below was measured on one machine against one dataset:

| | |
|---|---|
| CPU | AMD Ryzen AI 9 HX 370, 12 cores / 24 threads |
| Cache | L1d 48 KiB per core, L2 1 MiB per core, L3 24 MiB total |
| RAM | 31 GiB |
| Storage | NVMe, 967 MB/s cold on this dataset (19.4 frames/s) |
| Dataset | `39_1-rsmap`, La3Ni2O7, CHESS QM2, 2527 x 2463 = 6,224,001 px/frame, 49.8 MB/frame |
| Job | `center` accuracy (depth 0), 2000^3 hkl grid, 64^3 chunks, 24 threads, 10 GiB budget |

The *reasoning* generalises. The *constants* are calibrated to this machine —
notably the 16384-pixel work block (an L2-sized working set), the 25%
accumulator share, and one record per pixel as the per-frame record estimate.
On a machine with a different L2, or a detector with very different record
density, they should be re-measured rather than trusted.

## Where the time goes

Per frame, one kernel thread, before this round of work (native grid):

| Stage | ms/frame | Share |
|---|---|---|
| load | 48.3 | 4.4% |
| correct_frame | 42.3 | 3.9% |
| slice/copy per tile | 24.0 | 2.2% |
| **kernel** | **853.4** | **78.3%** |
| reduce (merge tile batches) | 121.2 | 11.1% |
| **total** | **1089.3** | 0.92 frames/s |

The kernel dominates a single frame, but it releases the GIL; the other ~114 ms
(10.5%) does not. That is the ceiling on image-level parallelism — see open
item 1.

**This table is the per-frame path.** A grouped run redistributes it rather
than shrinking it: `load` and `correct_frame` move out of the compute worker
into a prepare pool, where they overlap the next group's kernel instead of
queueing behind it, and `reduce` shrinks because one merged batch replaces
`F` per-frame batches. The kernel's own per-sample work is unchanged — which
is why the 34% record saving converts to only ~12% wall time.

## Subdivision hierarchy, with footprints

One image passes through seven levels. Sizes are for the job above.

| # | Scope | Structure | Footprint |
|---|---|---|---|
| 1 | per pixel | `std::vector<VoxelWeight>` + recursion stack | tiny; **absent at depth 0** (fast path); 1.2-3.5 KiB of stack at the deepest setting since corner-passing replaced the dense cache |
| 2 | per block | `std::pmr::map<RecordKey, RecordAccum>` (arena-backed) | ~0.5 MiB |
| 3 | per block | `std::vector<Record>`, sorted snapshot of (2) — **all blocks live at once** | records x 40 B |
| 4 | per call | loser-tree merge of block runs -> six NumPy arrays | ~17 MiB per tile |
| 5 | per frame | `tile_batches` + `_reduce_batches` | ~100 MB |
| 6 | **per checkpoint** | `_CheckpointAccumulator._levels`, binary-counter tree, spans ~85 frames | **the dominant term** |
| 7 | finalization | `_CheckpointRangeReader` + one output chunk at a time | ~1 GiB |

Note there is **no thread-level structure**: a worker thread processes blocks in
sequence, reusing one arena. The loser tree merges *per-block* runs at the
*call* level.

With frame grouping the *scope* of levels 2-5 changes but their structure
does not. A "block" becomes a brick in (row, column, frame) holding the same
sample count, so level 2's working set is unchanged and level 3 falls with
the record density. Level 5 becomes per *group* — one merged batch for `F`
frames, measured 0.575x the sum of the `F` it replaces — and level 6 gains
a new resident term of its own: the group buffer, `F` full-detector
corrected frames, which is the largest single allocation in the pipeline
and what bounds how many group calls can run at once.

Other fixed costs: detector ray cache 143 MiB for the job; per worker thread a
1.13 MiB arena; grid chunks 32,768 of 7 MiB dense each.

## Settled — please do not re-litigate without new evidence

**The accumulator container stays `std::pmr::map`.** Measured alternatives at
realistic redundancy (block 4096, ~7:1 in-block):

| Container | ns/pixel | Bytes/key |
|---|---|---|
| `pmr::map` + last-key cache (current) | **30.1** | 72 |
| open-addressing hash + sort | 27.6 | ~80 effective |
| tlx B+ tree, key -> index + side array | 32.8 | 55 |
| tlx B+ tree holding values | 41.4 | 80 |
| sorted array, dedup on insert | 126.5 | 40 |
| vector + sort | 40.4 | 40 |

The B+ tree's memory win (~1.3x) costs speed and needs a side array, because a
B-tree relocates values on split and the last-key cache holds a pointer. The
sorted-run idea fails because the key stream is *not* monotone at any block
size: `backscans/pixel` stays at 0.97-0.996 from 1024 to 65536-pixel blocks. A
`(chunk, local)` key blocks reciprocal space in three dimensions, so one voxel
of movement in x jumps the key by a whole chunk row.

**Tiling does not reduce memory pressure.** 1 / 6 / 12 tiles measured 9.61 /
10.12 / 9.81 GiB peak — flat. Tiling subdivides *one frame's* native call; the
memory that matters accumulates *across* frames in the router, which no tiling
arrangement touches. Fewer tiles is better for throughput — *within a
limit*, see the next entry, which was measured later and bounds this one.

**Band height is flat, and "fewer tiles is better" has a limit.** Swept in
2026-08 at four frames per group over 240 / 421 / 506 / 600 rows,
throughput is flat across the region that matters and the height the memory
budget already produces sits in the middle of it:

| bands | tiles | concurrent calls x threads | median ms/frame |
|---|---|---|---|
| 240 | 11 | 4 x 6 | 236.7 |
| **421 (the default)** | 7 | **3 x 8** | **205.4** |
| 506 | 5 | 3 x 8 | 207.1 |
| 600 | 4 | 2 x 12 | 237.1 |

Note the last row against the claim above it: **four tiles is worse than
seven**. Fewer tiles is better only while the tile stays small enough for
the memory budget to afford the same concurrency — past that, a larger tile
costs a whole concurrent call and the tile-count saving does not come close
to paying for it. Between 421 and 506 rows, where concurrency is equal, the
two are indistinguishable, which is the form the original finding should
have taken.

**Frame grouping pays only where the kernel dominates and records
dominate the kernel.** Measured in 2026-08 across two grids and three
depths: the record saving transfers everywhere and even improves with a
coarser grid (0.450x at 1000^3 depth 0), but the *time* saving appears
only at depth 0 on the fine 2000^3 grid (0.872x, every pair). Depth 1 and
2 are nulls because each (pixel, frame) still walks its own subdivision
tree — ~31 evaluations per pixel at continuous depth 1, ~149 at depth 2,
against ~1 at depth 0 — so record density stops controlling kernel time.
Coarse-grid depth 0 is a null for the opposite reason: it is not
kernel-bound at all. Detail, and what it retires, in the locality
document's phase 4.

**Nothing above `balanced` changes a reconstructed intensity.** Depth 3
differs from depth 5 by 0.013-0.024 of a voxel's own error bar, needing
~10^5 contributors per voxel to reach 1 sigma against the ~1,000-2,500 a
full job delivers. Depth 0 -> 1 is the only clearly real step. See the
locality document, "The prior question, answered".

Band height and group size also substitute for each other almost exactly —
at depth 2, widening the bands makes the group-size gate halve the group
and the two cancel — which is why neither repays hard optimisation. Detail
in the locality document's phase 3 step 7.

**The transition-rate correlation is container cost, not branch
misprediction.** `voxel_id`'s cost is flat across transition rates (17.3 ->
17.4 ns/px from 0% to 100%); what scales is the `pmr::map` insert (0.5 -> 90
ns/px). An earlier session's 0.9959 correlation was measuring unique-key
working set.

**Reciprocal arithmetic was implemented, measured, and rejected.** Replacing the
ray-normalisation and grid-step divisions with reciprocal multiplies is worth
~35% in a cache-resident microbenchmark but only ~3% on real data, where those
divisions hide under memory stalls. It shifts a coordinate by at most 1.28e-12
voxel widths — zero assignment changes across 9.2M real records and 600M
uniform samples, but 4.7% for coordinates sitting exactly on a voxel edge,
which lattice-aligned synthetic data produces routinely.

## Open items

### 1. GIL-serialised per-frame Python work

The largest remaining structural limit. Load + correct + per-tile copy is
~114 ms of a 1089 ms frame and holds the GIL. Amdahl on a 10.5% serial fraction
predicts a 6.1x ceiling at 16 workers; measured 4.2x (26% efficiency; 14% at 24
workers on the coarser grid). Roughly two thirds of the loss is explained by
that serial fraction, the rest most likely memory bandwidth.

Candidate directions: move the per-tile slice/copy into the native call, or
release the GIL across correction. Note the memory work reduced worker counts,
so re-measure the balance before assuming the old figures still hold.

**Partly addressed, structurally rather than by tuning.** The grouped
scheduler hoists `load` and `correct_frame` into a prepare pool, so one
group's GIL-held Python overlaps the previous group's native call instead of
competing with it. That is the direction this item asked for, and it is what
turned frame grouping from a 1.4x regression into a 0.88x win. It applies
only above one frame per group; the per-frame path still has the serial
fraction described above, and the per-tile slice/copy is still in Python on
both paths.

**Now the binding constraint, and measured rather than inferred.** Phase 4
turned this up by accident. On the current (coarse) reference job at depth
0 the whole pipeline runs at ~110 ms/frame while using **7-8 of 24 cores** —
sampled per process during live runs, not derived from a model. Mapping is
not the limit there at all; loading and correcting are. The Amdahl estimate
above put the ceiling at 6.1x and measured 4.2x, which understated how
lopsided this gets once the kernel is cheap: a coarse grid, a fast disk and
depth 0 leave two thirds of the machine idle. **This is the largest
remaining win in the mapping phase**, and it is now better evidenced than
anything in the locality document's remaining ideas.

Worth **up to 2.1x, not more**: cold reads cap the pipeline at 51.5
ms/frame (49.8 MB at 967 MB/s) however free the CPU work becomes, and an
earlier draft that reasoned from idle cores alone and claimed 3x was
wrong. This item is now a plan of its own —
`reciprocal_space_mapping_serial_fraction.md` — and further work on it
belongs there rather than here.

**Measured, and the diagnosis in this item is wrong.** The serial-fraction
document's step 1 put a probe on the GIL itself: a Python thread wanting
to run waits **23%** of its wall time during a depth-0 run (1% idle), on
both grids. The pipeline is not GIL-bound, the prepare pool that was
supposed to be blocked is idle 87% of the time, and "loading and
correcting hold the GIL, so adding workers does not help" is not what
holds this pipeline at 6.4 of 24 cores. What does is memory: the same
correction costs 2.3x the cycles inside the pipeline that it costs alone
on one thread. The useful consequence is that **work removed anywhere
helps, including off the critical path** — which is not what a serial
fraction predicts. Two changes came out of it (0.84x and 0.92x); see that
document.

### 2. `bounded_block_size` still clamps memory with a depth-blind constant

The Python side now caps the work block against the arena the kernel will
reserve (`work_block_memory_cap`). The kernel's own `bounded_block_size()`
still divides its budget by a magic `512` that ignores `worst_leaves`, so the
two disagree about what a block costs.

**The conclusion this item drew has been superseded.** It argued the
authoritative cap belongs in the kernel, because the kernel knows whether an
exposure is stationary while the Python side must assume the worst. Both
memory prechecks are now bounded by the record ceiling instead of the leaf
count, and that bound saturates above depth 0 — so neither side depends on
the exposure model any more and the asymmetry disappeared rather than being
relocated. See the locality document's phase 2.

Still open, and narrower than it was: `bounded_block_size`'s magic `512` is
untouched, and a third site — the detector band height in tile planning —
was found in 2026-08 still using `8**depth` and has been corrected. Any
remaining depth-blind constant in this area is worth checking against the
same record ceiling.

### 3. Depth 3+ block sizes are extrapolated

The `16384 >> depth` rule was measured at depths 0, 1 and 2 (optima 16384, 8192,
4096, consistent at 12 and 24 threads). Depths 3-8 follow by extrapolation only.

**Mostly retired by the depth-convergence result.** Depths above 2 change
no intensity a measurement can resolve, so tuning their block size
optimises settings nobody should be running for production. Worth doing
only if someone is using them deliberately to check convergence, and even
then the runtime is dominated by the subdivision itself rather than by
block-size choice.

### 4. `image_workers = min(len(pending_ranges), ...)`

Concurrency is capped by the *scheduling* granularity. Production sizes
`frame_batch` for roughly 4x as many ranges as workers, so it is fine mid-job,
but the tail of every job runs with fewer workers by construction. It is also
the single most effective way to invalidate a benchmark — see below.

**Gone above one frame per group.** The grouped scheduler does not use
`image_workers` at all: work is a queue of frame groups drained by a fixed
number of concurrent native calls, so scheduling granularity no longer bounds
concurrency. The per-frame path is unchanged, and the benchmark trap still
applies to it.

### 5. Finalization scaling at full job length

Shrinking the accumulator share raises the checkpoint *part* count, and
finalization k-way merges across all parts per output chunk. Measured flat over
2 checkpoints (finalize 38.7-40.6 s across a 4x accumulator range), but a real
job has 43 checkpoints, where the run count would grow roughly 86 -> 300, i.e.
log2 6.4 -> 8.2. Still on a non-dominant term against HDF5 compression, but
unverified at scale.

### 6. The probe's reported error is conservative by construction

`records_per_pixel_relative_error` treats the Sobol tiles as independent while
the sequence is stratified, so it overstates true scatter by 2-4x (reported
5.5-7.5% against observed 1.6-3.4%). Safe for a margin, not calibrated. A
calibrated figure needs randomised QMC — several independent scrambles, variance
across their means — which is too noisy at the 10-15 tiles a 0.1 s budget buys.

### 7. Two fixes are unexercised on this hardware

The rebalance unit fix (sample tile -> whole frame) and the adaptive re-check
cadence only change a decision when the delivery rate is high enough that a
candidate's requirement exceeds its worker ceiling. Here mapping is
compute-bound — ~5 frames/s against 19 available — so the requirement stays at
its floor either way. They need a machine where I/O outpaces compute.

## Measurement traps

Every one of these produced a confident, wrong conclusion during this work.

**`image_workers = min(len(pending_ranges), ...)`.** Passing one pending range
pins the compute pool to a single worker regardless of thread settings. This
invalidated two separate rounds of conclusions, including a recommended
accumulator constant that had to be retracted: with one worker, shrinking the
accumulator looked like it *halved* peak memory; with a real pool it *raised*
it, because faster mapping puts more frames in flight.

**Harness costs inside the timed region.** An early tiling comparison timed its
own `np.ascontiguousarray` copies and reported 2.4x where the real figure was
1.13x.

**Unequal work between arms.** Synthetic per-tile geometry that gave each
"tile" the full angular range produced 8x the records in the split arm.

**Overriding depth without overriding block size.** `max_depth` alone leaves
`work_block_pixels` at its depth-0 value, inflating the kernel's arena from
4.8 GB to 77 GB and making depth 8 look broken. Always derive both.

**This machine is bimodal.** Identical binaries measure 42 or 71 ns/pixel
depending on thermal/clock state. Single samples at high thread counts are
worthless; alternate arms and take best-of-N.

**End-to-end pipeline timings are worse than bimodal.** Single runs vary by
up to 30% — the same configuration measured 183, 202 and 342 ms/frame on
different occasions — and absolute values drift between *sessions*, so a
number written down yesterday cannot be compared with one taken today. Only
paired, interleaved ratios carry, nothing else may run on the machine while
measuring (a `pytest` run in another shell is enough), and a clean
separation between arms is worth more than a difference of means. Two rounds
of conclusions in the locality document were withdrawn for this. Its
"Measuring any of this" section is the full version.

**The reference job's own configuration drifts.** `39_1-rsmap.json` on
disk is no longer the job the table at the top of this document
describes: it now carries `high` accuracy and a 1000^3 grid, every step
2x the documented one. That is a different regime — at depth 0 it maps at
~110 ms/frame using 7 of 24 cores, bound by loading and correcting rather
than by the kernel, so phase 3's whole result is invisible on it. A 2026-08
sweep was run and had to be reinterpreted. Print the grid before trusting
any comparison against numbers written here.

**At depth 0 the page cache is the noise.** A 234-frame window is 11.7 GB
against 31 GB of RAM and the run itself holds several more, so how much
of the window survives between repeats varies. One fine-grid depth-0 arm
spread 176.5 to 340.6 ms/frame — a factor of 1.9 — where depth >= 1 arms
hold ±7%. Six pairs, not three, at depth 0.

**"Nothing else may run" is worth verifying, and worth measuring
correctly.** This machine carries ~1.8-2.9 cores of persistent background
load. Sampling it per run turns that from an invisible confounder into a
number you can gate on — but define it as *the sum of other user
processes*, not as machine-busy minus the benchmark's own. Kernel,
interrupt and storage-driver time is attributed to System rather than to
the process that caused it, and mapping is I/O-heavy, so the
total-minus-child definition charges the benchmark's own I/O to
"foreign": it reported 7.2 cores of interference where per-process
accounting showed 2.9.

**A run that maps nothing looks like a fast run.** Roughly 2 of ~40
pipeline runs in 2026-08 routed zero records, wrote a checkpoint with 0
rows and `frames_covered` set to the full count, and exited 0 — at 50.5
ms/frame against ~250, i.e. a 5x "speed-up". Always compare a run's
checkpoint fingerprint against its partner's and discard mismatches; a
timing-only harness will silently take these as wins.

*It happened again during the serial-fraction work, once in about 50
runs, and `benchmark_reconstruction_ab.py`'s fingerprint check caught it
automatically — which is what that check is for.* Two details worth
adding to the description above. It was on the **unchanged** arm of the
comparison, so whatever it is has nothing to do with the changes being
measured. And it was **46.4 ms/frame against 71-100 for its neighbours**,
under half rather than a fifth: the "5x faster" signature is not
reliable, because how fast a real run is depends on the page cache. The
fingerprint is the reliable signal; the clock is not.

`_map_pending_ranges` now warns, on both schedulers and through the
progress channel as well as `warnings`, when a whole run routes nothing.
That makes the next occurrence visible; it does not explain it.

**The cause is still unknown, and the search was stopped deliberately
rather than exhausted.** Everything below was measured, not reasoned
about. What it adds up to is a reframing worth having before anyone
starts again: *every static input the run is built from is
deterministic*, so this is not a bad value read at setup — it is a race
or a load-dependent condition during execution.

Bit-identical across repeated fresh job loads, so none of them can be
the poison:

| input | repeats |
|---|---|
| solid angle and polarization (`static_factor`) | 12 |
| exposure angle bounds, grid minimum/step/shape | 20 |
| UB, U, wavevector K (`det(U)` = 1.000000000) | 8 |
| detector corner rays (all unit norm) | 6 |
| static mask (0.0857 masked) | 12 |

Also eliminated:

- **Not a swallowed exception.** Both schedulers re-raise the first
  recorded exception, and the failing runs exited 0 with no traceback, so
  correction did not raise.
- **Not the native memory precheck.** It throws; it has no path that
  returns empty.
- **Not the pixel repair plan**, despite being the best structural
  candidate: a shared native object called concurrently by every prepare
  thread, whose return value *replaces* the mask outright. Its
  `apply_inplace` is `const`, holds no `mutable` members, and keeps all
  scratch local to the call. Thread-safe.
- **Not reproducible on demand.** 108 attempts across four
  configurations — subprocess and in-process, coarse grid and fine,
  grouped scheduler and per-frame, depths 0 and 1 — returned the
  identical record count and the identical 0.0857 masked fraction every
  time. Against an observed rate of ~2 in 40, that is not merely bad
  luck: whatever triggers it is not uniformly random per run.

**A second, milder failure of the same family was caught in 2026-08, and
it is the sharpest evidence yet.** One run in twelve of an interleaved
sweep reported a total weighted intensity of **2,240,706.86 against
1,839,128.67** — 22% high — and a total weighted variance of **236.23
against 136.96**, 72% high, while the other eleven agreed with each other
to one or two units in the last place. In the same run:

| quantity | anomalous run | the other eleven |
|---|---|---|
| voxel fingerprint | identical | identical |
| contributors | 1,329,409,714 | 1,329,409,714 |
| total weight | 1,329,409,714.0 | 1,329,409,714.0 |
| rows written | 47,226,553 | 47,226,553 (three others) |
| **weighted intensity** | **2,240,706.86** | 1,839,128.665149 |
| **weighted variance** | **236.23** | 136.960292 |

Read what that excludes. Identical voxels and identical contributor
counts mean the geometry, the angle bounds, the grid *and the mask* were
all right: the same detector samples reached the same voxels. Only the
*values* they carried were wrong. So this is not a geometry failure at
all — it is the pixel data or the factors applied to it, which is the
same half of the pipeline the zero-record failure points at, reached from
the opposite direction.

And the two totals moved by **different** factors, 1.22 against 1.72.
A wrong image would move both alike, since the variance is the clipped
image itself. A wrong *scalar factor* moves the intensity by `k` and the
variance by `k**2`, which is the shape actually observed. The per-frame
exposure and monitor factors are the only per-frame scalars in the
pipeline, and they are read from arrays captured once at pipeline
construction — so on the face of it they cannot vary, which is exactly
what makes this worth writing down rather than explaining away.

Two honesty notes. The run was on the **unchanged** arm of the
comparison, so it is not caused by whatever was being measured. But it
was on a build carrying the fused native correction pass
(`apply_correction_factors`), which is verified bit-for-bit against the
NumPy form by test and over a 78-frame real window — verified, not
exonerated. Anyone hunting this should not assume that pass is innocent
merely because its unit test passes.

**The instrument to use next time is one bit.** Only two things can
produce an empty result: every pixel masked, or every coordinate outside
the grid. So record the masked fraction when it happens. Normal is
0.0857; ~1.0 means the raw data or the correction went non-finite, and
0.0857-with-no-records means the pixels were fine and the geometry
applied during mapping is at fault. That single number halves the search
space, and the guard above is what will tell you a run is worth looking
at.

*Add a second: the per-frame factor actually applied.* The 22% run above
says the failure mode is not always all-or-nothing, and that a run can
look completely healthy on every count anyone currently checks —
records, voxels, contributors, masked fraction — while carrying wrong
values. Logging the exposure and monitor factor per frame, or simply the
per-frame corrected sum, would separate "the image was wrong" from "the
factor was wrong" the first time it happens again. Total weighted
intensity per checkpoint is worth comparing between runs for the same
reason; `benchmark_reconstruction_ab.py` now does exactly that, which is
how this one was found.

Both observed failures came late in long, memory-heavy sequences with
desktop applications resident, which the reproduction loops could not
recreate. A separate one-off hang of the same pipeline under `pytest`,
also not reproduced, may or may not be related.

**The hang happened twice more in 2026-08, and this time it was caught
alive.** Both were the pipeline benchmark on unchanged code. The state
of the hung process is specific enough to be worth writing down, because
it rules out most of where one would look first:

| | |
|---|---|
| checkpoints | **all three written**, full size, on disk |
| threads | 58, still alive |
| resident memory | 883 MB |
| CPU consumed before the hang | 188 s |
| CPU while hung | **0.00 cores**, indefinitely |

So mapping had finished and every record had been flushed; the process
then stopped with its pools still up and nothing running. That points at
the coordinator's exit condition or a pool join, not at mapping,
correction or the kernel — and it is a *different* failure from the
zero-record one, which completes and exits 0. `benchmark_reconstruction_ab.py`
now kills and retries a run that exceeds a per-run timeout, so a hang
costs a run rather than a sweep; it still deserves a stack dump the next
time one is caught alive (`py-spy dump --pid`), which is the one thing
not done here.

## Calibration constants and what they rest on

| Constant | Value | Basis |
|---|---|---|
| `WORK_BLOCK_PRESETS["medium"]` | 16384 px | L2-sized working set; measured optimum at depths 0-2, both thread counts |
| block scaling | `>> depth` | measured optima 16384 / 8192 / 4096 at depths 0 / 1 / 2 |
| `CHECKPOINT_MEMORY_SHARE` | 0.25 | peak tracks 2x the accumulator budget; 4x shrink cost +6% scratch |
| `_FRAME_RECORD_BYTES_PER_PIXEL` | 2 x 40 B | bounds measured density 0.46-0.87 rec/px, doubled for the merge |
| `_RESERVED_RECORDS_PER_PIXEL` | 4 | 6-8x headroom over measured density; heap fallback covers the rest |
| `_REBALANCE_INITIAL_SECONDS` | 30 | first check well inside a short job; backs off x2 to the old 600 |
| probe `budget_seconds` | 0.1 | holds 0.04-0.10 s at every depth 0-8; GUI limit is 0.5 s |
| `_GROUP_ADVANCE_VOXEL_LIMIT` | 1.0 voxel | striding sweep: the grouping gain is already gone by one voxel of travel per frame |
| `_GROUP_MINIMUM_CONCURRENT_CALLS` | 3 | **purely empirical**: 4 frames/3 calls measured 0.87x, 8/2 was 1.01x, 8/1 was 1.24x |
| `_GROUP_SIZE_CANDIDATES` | 1, 2, 4, 8, 16 | powers of two spanning the useful range; 4 is what the reference job picks |
| `_GROUP_PIPELINE_*_READAHEAD` | 1 to 3 groups | one is the minimum for double-buffering; more costs a full group buffer each and buys nothing measured |

The group constants are the least well-founded in this table. The voxel limit
comes from a striding sweep on one scan; the concurrency floor of 3 is fitted
to a single job's throughput curve and reproduces its optimum but has not
been tested anywhere else. Treat both as calibrated to this machine and this
scan in the same way the work-block size is.

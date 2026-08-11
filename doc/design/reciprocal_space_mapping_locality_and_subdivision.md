# Reciprocal-space mapping: locality, subdivision cost, and the scheduling that would carry them

> **Status: phases 1-3 implemented, phase 4 not.** Frame grouping is on by
> default and chooses its own group size per job. Each step below says which
> it is and what it measured. Everything here was measured against real
> beamtime data, and where a prediction was made before implementing, the
> step records whether it held — including where it did not.

Third companion to `reciprocal_space_scratch_architecture.md` (what the pipeline
*is*) and `reciprocal_space_performance_findings.md` (what has been measured
about it). This one holds two findings that neither of those settles, the
places where they meet, and what the scheduler and the memory budget would have
to do to carry them.

- **A — cross-frame locality.** The native work block deduplicates voxels within
  one image, mostly by accident. Extending it across adjacent images cuts
  emitted records by ~42% and *shrinks* the cache working set. **The kernel and
  the pipeline wiring, the scheduler and the per-job group size are all
  implemented.** Records through the router fall to 0.657x at the group size
  chosen for the reference job, and end-to-end mapping to **0.88x**. See
  phase 3.
- **B — subdivision cost.** Above depth 3 the shared-corner cache switched off
  and the kernel did several times more ray work than it needed to. A
  corner-passing recursion removes that, bit-for-bit. **Phase 1 is implemented**
  — 1.65-2.98x stationary and 1.85-2.85x continuous at depths 1-4 on real data,
  every output array bit-identical. See the phase 1 steps for what landed and
  the one case that did not.

B landed first: it is smaller, it changes one file, its output is bit-for-bit
identical, and it is what makes A viable at depth > 0 at all — see "Where the
two meet". The measurements in "Finding B" below are the *pre-change* state,
kept as the record of why the work was done.

### Where this stands

Everything is on branch `perf/subdivision-corner-passing`, unpushed.

| phase | state |
|---|---|
| 1 — corner-passing subdivision | done, bit-for-bit, 1.65–2.98x |
| 2 — the `worst_leaves` memory bound | done |
| 3 — frame-group bricks at depth 0 | done, steps 1–7; **0.88x end-to-end** |
| 4 — frame-group bricks above depth 0 | **not started** |

What a reader should know before touching it:

- **Grouping is on by default** and sizes itself per job. On the reference
  scan it picks four frames per call and maps at 0.88x. Jobs outside the
  regime (coarse scans, interlaced order, tight memory) get one frame per
  call and are unaffected.
- **Output moved.** Contributions merge inside the kernel rather than in
  the checkpoint accumulator, so sums associate differently — in the last
  bits, and better conditioned. Voxels reached and contributor counts are
  unchanged, and tests pin that. This is a `phys` change and it applies to
  every job, not only grouped ones.
- **Read "Measuring any of this" before running a benchmark.** Two rounds
  of conclusions in this document had to be withdrawn for measurement
  error, not for reasoning error.
- **The next real question is phase 4's first measurement**, not phase 4's
  design: `F = 1` against automatic `F` at depth 1 and 2, interleaved. The
  incidental depth-2 numbers suggest grouping does less there.

## How this was measured

All numbers are from the `39_1-rsmap` job used throughout the findings document
(La3Ni2O7, CHESS QM2, Pilatus 6M, 2527 x 2463, 3651 frames, omega rotation
0.1 deg/frame over 365 deg, `center` accuracy = depth 0, 2000^3 hkl grid with
steps [0.00457, 0.00208, 0.00265]).

Which voxel a pixel reaches depends only on geometry and the mask, never on the
recorded intensity, so the depth-0 pixel-to-voxel map was reproduced in NumPy
with the affine ray-to-coordinate transform extracted from `kernel.coordinate()`
itself — every sign and frame convention comes from the kernel rather than from
a re-derivation. Validated against a real `accumulate()` call on frame 1800:

| | native kernel | model |
|---|---|---|
| block records, 16384-pixel runs | 2,690,636 | 2,698,000 (+0.3%) |
| records after the cross-block merge | 2,577,418 | 2,577,295 (-0.005%) |

Depth scaling (finding B) was measured directly through the kernel's own
`profile=True` counters, not modelled.

## Test configuration: continuous exposures

Use `angle_fallback="midpoint"` as the standard configuration for this work.
QM2 rotates continuously, so it is the physically right model, and it is free at
depth 0: `_midpoint_bounds` centres each sweep on the stationary angle to within
1.3e-7 rad, and depth 0 samples only `transforms[size/2]`. Re-running the whole
locality measurement under continuous bounds reproduced the stationary result to
one part in 10^6 (43,087,060 records against 43,087,014 over a 16-frame window;
32x32x16 relative 0.5749 either way).

It stops being free from depth 1 up, where `stationary=false` switches the
kernel from a 4-corner quadtree to an 8-corner octree. That is exactly the path
with no hand-written fast case and the most to gain from finding B, so testing
under continuous bounds exercises the harder half.

## Finding A: adjacent frames are as close as adjacent pixels

Median displacement of a pixel's reciprocal-space coordinate, in voxels:

| next column | next row | **next frame** |
|---|---|---|
| 0.49 | 0.80 | **0.71** (p10 0.22, p90 1.48) |

So (row, column, frame) is a roughly isotropic sampling lattice, and the natural
work block is a **brick in all three**. The current block is a flat run of 16384
pixels — on a 2463-column detector, a 6.65 x 2463 x 1 strip, aspect ratio near
370:1 through a volume that wants 1:1:1.

Exhaustive — every block of a real edge-inclusive partition of a 16-frame window
over full frames, with the job's real mask. No sampling.

| block shape | records / valid sample | vs. current | live map |
|---|---|---|---|
| 7 x 2463 x 1 (**current flat run**) | 0.4741 | 1.000 | 522 KiB |
| 128 x 128 x 1 | 0.4580 | 0.966 | 467 KiB |
| 64 x 64 x 4 | 0.3056 | 0.645 | 307 KiB |
| 32 x 64 x 8 | 0.2809 | 0.593 | 294 KiB |
| **32 x 32 x 16** | **0.2725** | **0.575** | **282 KiB** |
| union over the whole window (floor) | 0.2602 | 0.549 | — |

Stable to three decimals at scan positions 0, 900 and 2700.

- **Fixing the block shape within one image is worth 3.4%.** A flat run does
  waste its perimeter, but a square detector tile is not where the volume is.
- **A 32x32x16 brick captures 94% of what a 16-frame window offers.** The
  redundancy is entirely local in frame index: after the first frame each
  additional one contributes ~1.404 M new voxels against its own 2.577 M, and
  that ratio is flat from F = 2. It needs a cache-sized brick, not a buffer.
- **The working set falls.** Same 16384 samples, fewer distinct voxels, about
  half the map nodes. The brick does not compete with the L2-sized work-block
  calibration; it relaxes it.

### The condition the whole idea rests on

Striding the frame window stands in for a coarser scan:

| deg/frame | shift (voxels) | 32x32x16 vs. flat run |
|---|---|---|
| 0.10 (this job) | 0.69 | 0.65 |
| 0.20 | 1.38 | 1.00 |
| 0.40 | 2.75 | 1.02 |
| 0.80 | 5.51 | 1.02 |

**Past roughly one voxel of travel per frame the gain is gone**, and a brick
becomes marginally worse than the current layout because it trades in-image
block area for frame extent that no longer overlaps. This is a statement about
grid step versus angular step: refine the grid and you move into the regime;
coarsen the scan and you leave it. `F = 1` — bit-for-bit today's behaviour —
must stay a reachable outcome, not a degenerate case of new machinery.

## Finding B: the subdivision is adaptive, but has a cliff at depth 3

A cell is split only if its corners disagree on the voxel, so work tracks the
length of the voxel-boundary curves crossing the pixel footprint: O(2^depth),
not O(4^depth). Frame 1800, 512x512 tile, one thread,
`work_block_pixels = 16384 >> depth`:

| depth | evals/px | x prev | leaves/px | rec/px | ns/px | x prev |
|---|---|---|---|---|---|---|
| 0 | 1.00 | — | 1.0 | 0.489 | 83.8 | — |
| 1 | 11.43 | 11.43 | 3.9 | 0.570 | 543.2 | 6.48 |
| 2 | 26.30 | 2.30 | 11.8 | 0.611 | 615.7 | 1.13 |
| 3 | 58.55 | 2.23 | 30.2 | 0.632 | 3090.4 | 5.02 |
| **4** | **398.21** | **6.80** | 69.8 | 0.642 | 9215.8 | 2.98 |
| 5 | 863.88 | 2.17 | 151.8 | 0.647 | 19341.4 | 2.10 |
| 6 | 1810.67 | 2.10 | 318.7 | 0.650 | 39801.6 | 2.06 |

The steady trend is x2.1-2.3 per level. **Except at 3 to 4, where it is x6.80** —
`begin_pixel`'s `if (max_depth > 3) { side = 0; return; }`. Above depth 3 the
shared-corner cache switches off and every corner of every cell is recomputed.
Redundancy factor 6.80 / 2.15 ~ **3.2x**, which is what you expect when each
dyadic lattice point is shared by up to four cells. `ACCURACY_DEPTHS` maps
`very_high` to 4 and `maximum` to 5: both run in the uncached regime.

### Continuous exposures are far worse

The table above is the stationary quadtree. Under `midpoint` bounds the kernel
runs the 8-corner octree instead, and both terms get worse: the tree grows x4
per level rather than x2 (a 3D domain cut by 2D surfaces), and up to eight cells
share a lattice point rather than four. Same real frame, 256x256 tile, one
thread, both arms measured by `benchmarks/benchmark_reconstruction_subdivision.py`:

| depth | stationary evals/px | ns/px | continuous evals/px | ns/px |
|---|---|---|---|---|
| 0 | 1.00 | 78.9 | 1.00 | 71.5 |
| 1 | 11.37 | 553.2 | 33.84 | 1467.2 |
| 2 | 26.05 | 600.4 | 150.81 | 8466.3 |
| 3 | 57.87 | 3086.8 | 672.92 | 41343.1 |
| 4 | 393.06 | 8952.2 | **13777.53** | **281743.5** |

The cliff at 3 to 4 is **x20.5** against a x4.5 trend — redundancy factor ~4.6,
consistent with eight-way sharing. In absolute terms depth 4 with continuous
exposures costs 282 us per pixel: 29 minutes of single-thread CPU for one frame,
74 hours of CPU for this job. **Continuous exposures above depth 3 are not
usable today**, which matters because continuous is the physically right model
for a rotation scan.

Two side observations from the same sweep. A swept exposure genuinely reaches
more voxels — 0.886 records/pixel at depth 4 against 0.605 stationary, and
0.684 against 0.538 already at depth 1 — so switching a job to `midpoint`
bounds raises record volume by up to ~46%, which the finding-A memory
accounting would have to carry. And the depth-4 moving maximum is 2325 leaves
for a single pixel against the `8^4 = 4096` the precheck assumes, so the bound
is nearer the truth in 3D than in 2D but the mean (1432) is still far below it.

Three things are wrong with the cache, in increasing order of importance.

**It is 3D even for stationary exposures.** `side^3` entries, but with `t` pinned
only the plane `it = (side-1)/2` is touched — 289 of 4913 at depth 3. The index
`(iu*side + iv)*side + it` strides 408 B, so every access lands on its own cache
line: 289 points cost 289 lines where a contiguous 2D array would need ~110. Per
worker thread the array is 3.0 / 17.1 / 115.2 KiB at depths 1 / 2 / 3 — past L1d
(48 KiB) at depth 3, and 842 KiB / 6.3 MiB at depths 4 / 5, which is why it is
capped at all.

**It is dense, but the walk is sparse.** ~20% of the lattice is touched at depth
3, ~11% at depth 6. Density is the wrong bet and gets worse with depth.

**The generic path pays for its machinery regardless.** Depth 2 (hand-written
`split_pixel_stationary_depth2`) does **2.3x the evaluations of depth 1 in 1.13x
the time**. Backing a raw evaluation out at ~23 ns from the uncached depths 4-6,
the generic path carries roughly 250 ns/pixel of stack traffic, cache indexing
and generation bookkeeping that the depth-2 special case simply does not have.

The cache is still a net win where enabled — about 28% at depth 3 — but it buys
far less than the structure it displaces could.

### Sub-pixel geometry is not the problem

For the record, since it comes up: `ray_at` bilinearly blends the four corner
unit rays rather than calling the detector model at the sub-pixel position.
Measured against exact `primBeamPoints` over 21,605 pixels spanning the frame,
the error is **2.1e-5 voxels median, 2.7e-5 max**, against a one-pixel step of
0.79 voxels — one part in 30,000 of a pixel. Replacing it would cost ~61 ns per
point (`verify_geometry_cost.py`) against ~7 ns, i.e. ~3.6x slower at depth 4,
to fix nothing measurable. The motor angles, meanwhile, are already exact:
`frame_rotation` interpolates the *angle* and takes its sine and cosine, so
every sub-sample gets a true rotation, and the dyadic `t` lattice is hit exactly
by `llround`. Keep the bilinear rays.

## Where the two meet: the per-pixel state budget

`ray_at(u, v)` is frame-invariant — bilinear interpolation of detector corner
rays plus a normalise. Only `M_f . ray + b_f` and `voxel_id` depend on the frame.
So in principle a frame-group brick could evaluate each lattice ray once and
reuse it across all F frames.

In practice, only at depth 0:

- *Which* lattice points are needed is frame-dependent — the tree follows the
  voxel boundaries, and those rotate. Only the shallow levels are reliably
  shared.
- A per-pixel dyadic ray cache costs `(2^(depth+1)+1)^2 x 24 B` **per pixel**:
  6.9 KiB at depth 3, 26 KiB at depth 4, 101 KiB at depth 5. Per pixel, so it can
  only be alive for one pixel at a time. That forces **pixel-outer, frame-inner**
  iteration at depth > 0 — the opposite of the frames-outer/pixels-inner order
  that makes depth 0 stream well. The brick loop needs two modes.
- The payoff is small anyway: the unit ray is ~30% of an evaluation, so
  cross-frame caching saves at most ~28%, against ~2.5x from finding B.

This is the dependency between the two findings. Corner-passing keeps the
per-pixel subdivision state at 1.2-3.5 KiB of stack instead of a 6.9-101 KiB
dense array, which is what makes pixel-outer/frame-inner iteration affordable
inside a brick at all. **At depth 0, cache centre rays for the brick (24 KiB for
32x32) and iterate frames-outer. Above depth 0, iterate pixels-outer and do not
build a cross-frame ray cache.**

## Memory

### What finding B costs and frees

Frees, per worker thread: the whole `PixelCoordinateCache` — 3.0 / 17.1 /
115.2 KiB at depths 1 / 2 / 3, and the 842 KiB / 6.3 MiB it would need at depths
4 / 5 if it were not disabled. Costs: recursion frames of 144 B (stationary) or
432 B (moving) per level, bounded by `max_depth <= 8`, so 1.2 / 3.5 KiB.

Nothing else moves. `work_block_memory_cap`, `_RESERVED_RECORDS_PER_PIXEL` and
the arena are untouched.

### The `worst_leaves` bound is what actually caps depth

Separate from both findings, and now well-evidenced. `native_bytes_per_pixel =
128 + 2 * children**depth * 40` appears twice, deliberately mirrored: in the
kernel's `accumulate` precheck and in `_frame_parallelism`. Measured maximum
leaves for any *single* pixel:

| depth | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| measured max | 55 | 142 | 319 | 679 |
| measured mean | 30.2 | 69.8 | 151.8 | 318.7 |
| `4^depth` assumed | 64 | 256 | 1024 | 4096 |

The real bound grows ~2.1x per level, because a pixel's leaf count is set by
boundary length, not area. Two compounding conservatisms: `4^depth` (or `8^depth`
on the Python side, which cannot know whether a range is stationary), and
assuming *every* pixel simultaneously hits its worst case. At depth 3 that is
5248 B/px, so a full frame "needs" 32.7 GB and the layout is forced into tiny
tiles; a 512x512 tile at depth 4 could not run at all until
`memory_budget_bytes` was raised to 2^50 for the sweep above.

Because `largest_tile_pixels` feeds `worker_memory`, which divides into
`image_workers`, this bound throttles concurrency across the whole depth range.
This is findings open item 2, it is a **memory** fix rather than a compute one,
and it should be its own commit — not folded into either finding.

### What finding A costs and frees

Per-pixel constants, for the 6,224,001-pixel detector:

| term | B/px | per frame |
|---|---|---|
| raw decoded image (`_RAW_IMAGE_BYTES_PER_PIXEL`) | 8 | 49.8 MB |
| corrected intensity + variance + mask | 17 | 105.8 MB |
| `_FRAME_RECORD_BYTES_PER_PIXEL` (batch + merge copy) | 80 | 497.9 MB |
| **`python_memory` per in-flight frame** | **105** | **653.5 MB** |

A frame group changes the shape of this:

- **Resident group buffer: F x 105.8 MB.** Only the corrected arrays need to
  stay; the raw image is released after correction. 846 MB at F = 8, 1.69 GB at
  F = 16. This is per *group*, not per worker — a group has ~6000 bricks and can
  occupy every thread on its own.
- **Correction transient: (correction workers) x 155.6 MB**, unchanged in
  character from today.
- **`_FRAME_RECORD_BYTES_PER_PIXEL` does not multiply by F.** One merged batch
  per group replaces F per-frame batches, and it is 0.573x the size of their
  sum.

The term that *does* grow is level 3 of the findings document's subdivision
table — every block's sorted `std::vector<Record>` live at once. Today ~108 MB
per concurrent frame; for a 16-frame group over full frames it would be 991 MB.

**Bound it with the detector tiling that already exists.** Bricks covering one
row band for all F frames, emitted and routed before the next band. A 256-row
band at F = 16 is 630,528 x 16 x 0.2725 x 40 B = **110 MB** — the same figure as
today. Frames stay resident across bands (you cannot load a band without the
frame), so the group buffer is unaffected; only the record footprint is banded.
Tiling therefore acquires a real purpose it does not have today, where the
findings measured it flat at 9.6-10.1 GiB across 1/6/12 tiles.

Net for the job: records through the checkpoint layer fall from 9.41e9 to
5.39e9, 376 GB to 216 GB. The checkpoint accumulator is the dominant memory
term and what caps concurrent frames at nine, so either the same
`CHECKPOINT_MEMORY_SHARE` holds ~1.75x as many frames per checkpoint — fewer
parts, cheaper finalization, which speaks to findings open item 5 — or the share
can shrink and hand the budget to the pipeline.

## Scheduling

Frame grouping does not fit the current Sec7 pipeline; it replaces the question
that pipeline exists to answer.

**Today.** `image_workers` compute workers each run a kernel with
`kernel_threads` native threads, `image_workers x kernel_threads ~ spec.threads`,
and `_map_pending_ranges` rebalances the pair live against the reader pool's
measured delivery rate, with `_kernel_threads_sweep` and asymmetric hysteresis.
The pair exists because one frame's kernel call cannot use all the threads
profitably and the memory budget caps how many frames can be in flight.

**With groups.** One `accumulate_group()` call over a 16-frame group has ~6000
bricks — enough to saturate 24 threads by itself. So `kernel_threads =
spec.threads` and one group at a time becomes the natural layout, and the joint
balancing problem largely dissolves. What replaces it is a simpler question:
**how many prefetch/correct workers keep one all-threads group call fed?**

The measured single-thread frame cost gives the answer directly. Python work is
load 68.9 ms + correct 30.1 ms = 99.0 ms per frame, GIL-held. Native work is
~490 ms per frame of block wall, GIL-released, and finding A should take the
record-dependent part of it down. At 24 threads a group's kernel occupies
`490 / 24 ~ 20 ms` of wall time per frame, so keeping it fed needs
`ceil(99.0 / 20.4) = 5` concurrent load-and-correct workers. Their CPU share is
only the correction part — 5 x 30/99 ~ 1.5 cores of the 24.

Consequences worth stating plainly:

- **The GIL-serialised fraction stops competing with the parallel region.**
  Findings open item 1 predicts a 6.1x ceiling at 16 workers from a 10.5% serial
  fraction and measured 4.2x. Double-buffering groups against a correct pool
  attacks that structurally rather than by tuning the split.
- **`image_workers = min(len(pending_ranges), ...)` dissolves.** Findings open
  item 4 — the tail of every job running with fewer workers by construction —
  goes away, because parallelism now lives inside one group rather than across
  scheduling ranges.
- **`_BoundedGate` retargets** from `image_workers + _PREFETCH_QUEUE_SLACK` to
  `F + slack`: it bounds frames in the group buffer.
- **The live rebalance survives, with a different variable.** `kernel_threads`
  is pinned at the budget, so what remains adaptive is `F` and the reader pool
  size, against memory pressure and delivery rate. The existing hysteresis and
  backoff machinery (`_REBALANCE_INITIAL_SECONDS`, `_REBALANCE_BACKOFF`) applies
  unchanged; only the quantity it steers changes.
- **A group is a synchronisation point.** Group N+1's kernel cannot start until
  all F of its frames are loaded and corrected. Double-buffering hides this only
  if the correct pool keeps up; if it does not, the pipeline stalls once per
  group rather than degrading smoothly. This is the main scheduling risk, and it
  argues for the smaller F (8) until measured.

## Measuring any of this

Read this before running a benchmark. Most of it was learned by getting a
measurement wrong first, and two rounds of conclusions had to be withdrawn
because of it.

### The instruments

| script | measures | typical cost |
|---|---|---|
| `benchmark_reconstruction_subdivision.py` | kernel cost and evaluations per pixel by depth and exposure, with XXH3 fingerprints and a `--compare` mode | seconds |
| `benchmark_reconstruction_group.py` | one `accumulate_group` call against per-image `accumulate`, by group size | seconds |
| `benchmark_reconstruction_pipeline.py` | the whole mapping phase — reader pool, correction, kernel, router, checkpoint writes | ~1 min per 234-frame arm at depth 0 |

`benchmark_reconstruction_mapping.py` is dead: it imports `_map_frame_range`
and the Parquet writers, both retired. Do not resurrect it; the pipeline
benchmark replaces it against the pipeline that exists.

### Two traps that silently invalidate a run

**`PYTHONPATH` must point at the checkout.** `python benchmarks/<script>.py`
puts `benchmarks/` on `sys.path`, *not* the repository root, so a
pip-installed `orgui` shadows the working tree and the run measures a
binary you did not build. Both reconstruction benchmarks print the resolved
extension path — check it. This has cost a whole measurement round.

**The reference job is `schema_version: 5` and `JOB_SCHEMA_VERSION` has
moved on**, so `read_job` rejects it. Copy the JSON and bump the field;
nothing else in the dict changed.

### Run-to-run noise is the hard part

Single runs of the pipeline benchmark vary by up to **30%**. The same
`F = 4` configuration measured 183.4, 202.0 and 342.4 ms/frame on different
occasions. Consequences, all of them learned the hard way:

- **Interleave the arms** (A, B, A, B, ...) with two or three repeats.
  Never one run of each in sequence — ordering effects and drift are
  confounded with the arm.
- **Nothing else may run on the machine.** A `pytest` run in another shell
  is enough to poison a sweep.
- **Only paired ratios are stable.** Absolute values drift between
  sessions: `F = 1` measured 233 ms/frame one day and 278 the next, in both
  arms alike, while the ratio held at 0.87. Never compare a number taken
  today against one written down yesterday.
- **Prefer a clean separation over a mean.** The `F = 4` result was
  trustworthy because every `F = 4` run beat every `F = 1` run, not
  because the means differed.

Two conclusions drawn from single sequential runs — a "1.47x regression"
at `F = 8` and a "0.85x win" from finer bands — both dissolved when
measured properly. Anything in this document from before that was
understood is flagged where it appears.

Baselines must come from the same checkout, captured immediately before
the change: `git stash push -- <changed files>`, run, `git stash pop`.
Never compare across builds.

### Measure the right quantity

**Records entering the router, not rows written to checkpoints.** The
checkpoint accumulator already merges across frames, so it collapses the
same redundancy later and more expensively and the saving never appears on
disk — rows moved 0.97x at `F = 8` where routed records moved 0.598x.

**Do not concatenate checkpoint parts to summarise them.** A 234-frame
window is ~350 M rows; reading it all at once swapped the machine on top
of the mapping run's own on-budget ~9.4 GB peak. Stream in row chunks.

**Where bit-for-bit is impossible, pin the contributor-weighted key sum.**
Frame grouping changes floating-point association by design, so a digest
moves for legitimate reasons. `sum over rows of h(key) * contributors`,
accumulated in wrapping `uint64`, is a sum and therefore independent of row
order and of how the accumulator split its parts. It equals the same sum
over distinct voxels with their total contributor counts, so it pins both
the voxel set and the contributor apportionment. Two baseline runs with
visibly different flush partitions (352.3 M against 350.1 M rows) agreed on
it exactly, which is what makes it usable as an acceptance criterion.

### Sizing a run

234 frames is six scheduling ranges of the reference job and about a minute
per arm at depth 0 — enough for three checkpoints to complete and flush.
Depth 2 is roughly seven times the kernel cost, so use 78 frames there
(~2 min per arm). A three-repeat, two-arm interleaved comparison is
therefore 6–15 minutes; budget for it rather than trying to shortcut it.

## Plan

### Phase 1 — Corner-passing subdivision (finding B)

Independent of everything else, bit-for-bit, one function.

**Mechanism.** A cell carries its own corner values; splitting evaluates only
the points the split creates. Stationary (quadtree over u, v): 4 corners in,
3x3 lattice out, **5 new evaluations** shared by 4 children. Moving (octree over
u, v, t): 8 corners in, 3x3x3 lattice, **19 new** shared by 8 children. Recursive
descent rather than the explicit stack, so the 9 or 27 lattice values live in the
caller's frame and children read them directly instead of each pushed cell
carrying a copy. This is `split_pixel_stationary_depth2` generalised — the
fastest per-evaluation path in the sweep, and an existing proof the structure
works.

What is given up is sharing across *parents*: a point on the boundary between
two cells with different parents is evaluated twice. Bounded, and the depth-2
path already accepts it.

**Lazy corners.** Neither `split_pixel` nor `same_lattice_voxel` exits early on
the first disagreeing corner. Under corner-passing that costs nothing for a
cell's inherited corners, which are already in hand, but it does mean the 5 (or
19) new points of a split are all materialised even when the first two settle a
child. Evaluate them on demand, as `split_pixel_stationary_depth2` already does
with its `evaluated` bitmask.

**Predicted, as acceptance criteria.** With `S` splits, the stationary quadtree
has `4S+1` cells: uncached probing costs `4 x cells ~ 16S`, corner-passing `5S`,
perfect sharing `~4S`. The octree has `8S+1` cells: `64S`, `19S`, `~8S`.

| exposure | depth | evals/px now | predicted | note |
|---|---|---|---|---|
| stationary | 3 | 57.87 (misses only) | ~72 | count *rises*, time falls |
| stationary | 4 | 393.06 (all probes) | ~155 | 2.5x |
| continuous | 3 | 672.92 (misses only) | ~1600 | count rises 2.4x |
| continuous | 4 | 13777.53 (all probes) | ~4100 | 3.4x |

Note the asymmetry: corner-passing recovers less of the ideal in 3D (19S against
~8S) than in 2D (5S against ~4S), because 19 new points serve 8 children rather
than 5 serving 4. It should still win below depth 4 in both, because the dense
cache's per-probe cost is paid on *all* `16S`/`64S` probes while corner-passing
pays only for evaluations — but **the depth <= 3 continuous case is the one that
could plausibly regress**, and it is the case to measure first rather than last.

Time: **2-2.5x at stationary depth >= 4**, **~3x at continuous depth >= 4**,
and 1.5-2x at depth 3 where the probe overhead disappears. Outside those bands
the model is wrong and that is worth understanding before proceeding.

**Bit-for-bit is the acceptance criterion.** Same lattice points, same all-same
test, same centre fallback, same weights; `voxel_id` is deterministic, so
re-evaluating a shared point yields the identical value. This is a pure
performance change with no `phys` impact. One thing to preserve deliberately:
leaf push order. The current stack pushes children 3->0 and pops from the back,
visiting 0,1,2,3 depth-first; `weights` is later sorted by voxel with a
non-stable `std::sort` and summed in vector order, so a different visit order
would change FP association. Recursing children in index order reproduces it.

**Steps.**

1. *Scaffolding, before touching the kernel.* **Done.**

   Golden digests turned out to be the wrong thing to check in: floating-point
   contraction differs between compilers, so a digest captured on one build
   fails on another for reasons that have nothing to do with the subdivision.
   Split in two instead.

   - `orgui/datautils/xrayutils/test/test_reconstruction_subdivision.py` —
     an independent Python implementation of the adaptive rule, which obtains
     every voxel through the kernel's own `coordinate()` so it cannot drift on
     frame conventions or rounding, and reproduces only the *algorithm*: which
     lattice points are visited, which cells resolve, what weight each leaf
     carries. It matches the kernel exactly at depths 0-4 stationary and 0-3
     moving, covering the centre-only fast path, `split_pixel_stationary_depth2`
     and generic `split_pixel` in both its cached and uncached regimes.
     Alongside it: weight/intensity conservation at every depth, thread- and
     block-size independence on full output arrays, sorted-unique output keys,
     and depth convergence. Portable, and it catches an algorithmic error rather
     than merely a change. 32 tests, 1.4 s.
   - `benchmarks/benchmark_reconstruction_subdivision.py` — the depth/exposure
     cost sweep plus XXH3 fingerprints of all six output arrays, with a
     `--compare` mode that reports bitwise agreement and per-case speed-up
     against a captured baseline. Runs on the synthetic geometry or a real job.
     Baselines live in `benchmarks/baselines/` (gitignored — they are per build
     and per machine).

   One thing the reference must *not* assert, and now documents: the emitted
   voxel count is not monotone in depth. A cell still unresolved at `max_depth`
   contributes at its centre, so raising the depth replaces one centre sample
   with four elsewhere and a voxel the coarser sample happened to hit can go
   unvisited. Real frames average that away; a small tile does not.
2. *Stationary recursion.* **Done.** `subdivide_stationary` /
   `split_pixel_stationary`: the 2D special case of the same shape — 9-node
   lattice, 5 new per split, 4 children — with `t` pinned so the rotation is
   chosen once by the caller rather than recovered from `t` per node.

   | exposure | depth | evals/px was | now | ns/px was | now | speed-up |
   |---|---|---|---|---|---|---|
   | stationary | 1 | 11.32 | 11.32 | 707.7 | 430.0 | **1.65x** |
   | stationary | 2 | 25.89 | 25.89 | 652.0 | 601.5 | 1.08x |
   | stationary | 3 | 57.44 | 64.61 | 3061.8 | 1502.7 | **2.04x** |
   | stationary | 4 | 389.67 | 142.87 | 9378.9 | 3147.0 | **2.98x** |

   **`split_pixel_stationary_depth2` is kept, not deleted.** Routing depth 2
   through the general recursion measured a consistent 10% *regression*
   (0.90x, best-of-9, on the default `balanced` preset), because that path
   holds one flat 5x5 lattice for the whole pixel and so shares nodes between
   *sibling subtrees* — which corner-passing gives up. 22.02 evaluations per
   pixel against 20.82. Depth 2 therefore still routes to it and the row above
   is unchanged behaviour.

   That result is worth more than the code it saves: it says the flat lazy
   lattice is the better structure *where it fits*. The reason to abandon it
   was never the idea, it was the shape — the removed cache was 3D even for
   stationary exposures, so at depth 3 it spent 115 KiB to use 289 scattered
   entries. A flat **2D** lattice is 4.6 KiB contiguous at depth 3, 17 KiB at
   depth 4, 68 KiB at depth 5. Generalising `split_pixel_stationary_depth2` to
   arbitrary depth might beat corner-passing across the stationary range the
   way it does at depth 2. Worth measuring; not attempted here.
3. *Moving recursion.* **Done, and taken first** — it is where the prize is,
   and it was the case that might have regressed. `subdivide_moving` /
   `split_pixel_moving`, 27-node lattice indexed `(iu*3 + iv)*3 + it`, 19 new
   per split, children visited in index order so leaf emission and therefore
   floating-point association are unchanged. `accumulate_block` routes
   `!stationary` to it; `split_pixel` and the dense cache remain, now used only
   by the stationary path.

   Measured on a real frame, 128x128 tile, one thread, against a baseline
   captured from the same checkout immediately before the change:

   | exposure | depth | evals/px was | now | ns/px was | now | speed-up |
   |---|---|---|---|---|---|---|
   | continuous | 1 | 33.91 | 33.91 | 1510.2 | 815.0 | **1.85x** |
   | continuous | 2 | 152.12 | 194.43 | 8913.9 | 4291.3 | **2.08x** |
   | continuous | 3 | 682.52 | 1009.03 | 42748.9 | 22161.5 | **1.93x** |
   | continuous | 4 | 14031.56 | 4662.31 | 283696.2 | 99668.8 | **2.85x** |

   All fingerprints identical. Stationary unchanged (0.97-1.10x, noise).

   Three things this settled. The feared depth <= 3 regression did not happen:
   evaluations *rise* by 24-48% from cross-parent duplication, yet time still
   nearly halves, because the dense cache pays its per-probe cost on all `64S`
   probes while corner-passing pays only per evaluation. Depth 1 is the clean
   case — one split, so corner-passing is exactly optimal and the evaluation
   count is unchanged; the entire 1.85x there is cache overhead removed. And
   the cross-parent loss is far smaller in practice than the worst-case `19S`
   against `~8S` model suggests, because not every cell splits.

   Deliberately deferred: lazy corner materialisation. All 19 new points are
   computed eagerly; only leaf children that disagree early would benefit, and
   the structural change was worth measuring on its own first.
4. *Delete.* **Done.** `PixelCoordinateCache`, `CachedVoxel`,
   `cached_voxel_key`, `begin_pixel`, `struct Cell`, generic `split_pixel`, the
   `max_depth > 3` cutoff and the per-block `stack`/cache scratch are all gone.
   `split_pixel_stationary_depth2` stays, for the reason in step 2. Net +133
   lines (329 added, 196 removed).

   A block's only remaining scratch is the leaf vector: subdivision state now
   lives on the recursion stack, 9 `LatticeVoxel` per level stationary and 27
   moving, so ~1.2 KiB / ~3.5 KiB at the deepest setting against a per-worker
   array that was 115 KiB at depth 3 and would have been 6.3 MiB at depth 5.

5. *Measure.* **Done.** The cliff is gone and `coordinate_evaluations` now
   means one thing at every depth — distinct evaluations — where before it
   counted misses below depth 4 and all probes above, which is what made the
   3-to-4 step read as x6.80. Scaling is now uniformly geometric:

   | depth | stationary evals/px | x prev | continuous evals/px | x prev |
   |---|---|---|---|---|
   | 1 | 9.94 | — | 31.05 | — |
   | 2 | 20.83 | 2.10 | 148.66 | 4.79 |
   | 3 | 48.16 | 2.31 | 695.15 | 4.68 |
   | 4 | 102.40 | 2.13 | 3038.42 | 4.37 |
   | 5 | 212.91 | 2.08 | 12739.36 | 4.19 |
   | 6 | 435.95 | 2.05 | 52205.21 | 4.10 |

   x2 per level for a 1D boundary set sampled in 2D, x4 for a 2D boundary set
   sampled in 3D — which is what the adaptive rule should cost and, before
   this, only did below depth 4.

### Phase 2 — The `worst_leaves` bound (memory, standalone). Done.

Both prechecks are now bounded by `reserved_records_per_pixel`, the constant
the arena already used, hoisted to file scope in the kernel so the two users
share one definition. Measured on the real detector against the job's own
10 GiB budget:

| depth | tiles | before | after |
|---|---|---|---|
| 0 | 1 | 5 workers, 1.93 GiB | unchanged |
| 2 | 1 | 1 worker, 9.63 GiB | 3 workers, 3.21 GiB |
| 3 | 1 | 1 worker, 31.03 GiB | 3 workers, 3.21 GiB |
| 5 | 1 | 1 worker, 476.2 GiB | 3 workers, 3.21 GiB |
| 3 | 4 bands | 1 worker, 9.63 GiB | 7 workers, 1.38 GiB |
| 5 | 4 bands | 1 worker, 119.6 GiB | 7 workers, 1.38 GiB |

It also closed the drift risk findings open item 2 raised, but not the way
that item anticipated. The concern was that a spec cannot know whether a
frame range's exposure rotates, so the Python side had to assume `8**depth`
where the kernel knew it was `4**depth` — and the conclusion drawn was that
the authoritative cap therefore belonged in the kernel. With the bound
saturating above depth 0, neither side depends on the exposure model at all,
so the asymmetry disappears rather than being relocated. `_frame_parallelism`
no longer takes `stationary`, and the two call sites that computed it only to
pass it no longer do.

What was *not* done: a per-thread term for the leaf buffer. It is the one
genuinely exposure- and depth-dependent allocation left, and it is the reason
the original bound looked plausible — but it is bounded by thread count
rather than tile size, and at 24 threads it reaches ~14 MB at moving depth 6
and ~230 MB at depth 8. Too small to model against a multi-gigabyte budget.

### Phase 2 — original reasoning, retained

Replace `children**depth` in both mirrored sites with a bound derived from
boundary growth (~2^depth) plus headroom, keeping the two in step. Findings open
item 2 argues the authoritative cap belongs in the kernel, which knows
`stationary` at `accumulate()` time while the Python side must assume the worst.
Unblocks realistic tile sizes — and therefore `image_workers` — from depth 3 up.
Not required by either finding; do not fold it in.

### Phase 3 — Frame-group bricks at depth 0 (finding A)

1. *Kernel.* **Done.** `accumulate_group()` takes `(frames, rows, columns)`
   arrays and `(frames, 4)` angle bounds and returns the same six columnar
   arrays merged across the group. The block queue is a brick index;
   `block_results`, the per-worker arena and the loser-tree merge are
   unchanged. The caller sets the group size purely by how many frames it
   passes, and one brick holds the same sample count a single-image block did,
   so the frame extent is spent out of the detector-plane extent rather than
   added to it.

   Measured on the real job, 512 x 2463 tile, depth 0, against the same frames
   mapped one at a time:

   | group | brick | records | time |
   |---|---|---|---|
   | 1 | 128x128x1 | x0.966 | x0.903 |
   | 2 | 91x90x2 | x0.752 | x0.752 |
   | 4 | 64x64x4 | x0.645 | x0.715 |
   | 8 | 45x45x8 | x0.595 | x0.606 |
   | 16 | 32x32x16 | **x0.575** | **x0.621** |

   Record ratios at one thread, time ratios at twelve. The record figure
   matches the exhaustive geometric prediction of 0.5749 made before any of it
   was implemented, which retires the main open question of finding A.

   It also settles the one the plan flagged as unpredictable: **the kernel
   itself is 1.6x faster**, not merely the stages downstream. About 42% fewer
   distinct keys are inserted into the block map, and the pixel centre ray --
   pure detector geometry, unchanged by the rotation -- is computed once per
   brick instead of once per (pixel, frame). A single-frame group is already
   1.1x faster from the block shape alone.

   One correction to the plan: **`F = 1` is not bit-for-bit with today.** A
   brick partitions the detector into rectangles where a block partitions it
   into runs of the flattened image, so pixels group into accumulators
   differently and their sums associate differently. Same voxels, same
   contributor counts, totals equal to rounding. The fallback is
   scientifically identical rather than byte-identical — and slightly better,
   at 3.4% fewer records.
2. *Brick loop, depth 0.* Precompute the brick's unit centre rays once into a
   24 KiB stack buffer, iterate frames-outer / pixels-inner. Note separately
   that at depth 0 the centre ray is frame-invariant for the whole job — pure
   detector geometry, currently recomputed 3651 times per pixel. A
   whole-detector centre-ray cache (149 MB, alongside the existing 143 MiB
   corner cache) removes it entirely, with or without grouping.
3. *Banding.* **Done, and it needed no new code.** The detector tiling that
   already exists is the banding: `_execution_layout` gives this job six
   full-width bands of ~421 x 2463 (1.04 Mpx), and `_map_frame_group` maps
   one band across all F frames before moving to the next. At F = 8 a band
   is 8.3 M samples, which the kernel's own precheck passes comfortably
   (208 B/sample against 10 GiB / 8.3 M = 1290 B available), and its record
   footprint lands near 93 MB — the same place as the ~110 MB figure above.
   Tiling has acquired the real purpose this document predicted it would.
4. *Wiring it into the mapping path.* **Done.** `_map_one_frame` is gone;
   `_map_frame_group` replaces it and is the only path, with F = 1 the
   degenerate case rather than a separate branch. `_ReconstructionSpec`
   carries `frames_per_group`. The pieces that had to move with it:

   - `_CheckpointRouter.route()` takes `frames=`, because the checkpoint
     countdown and the `frames_covered` attribute that makes a part
     resumable both count frames, not `route()` calls. A group declaring
     one frame would leave its checkpoint permanently short and
     unresumable. It also rejects a group straddling a checkpoint
     boundary, which is otherwise a silent corruption of that bookkeeping.
   - `_frame_groups` cuts groups at the group size and at every grid's
     checkpoint boundaries. Scheduling ranges are already contiguous and
     exclusion-free, so `excluded_frames` needs no separate handling —
     which is a weaker statement than this document's risk item assumed,
     and worth having checked rather than assumed.
   - `_angles_advance_monotonically` drops the job to F = 1 when the
     angles do not advance monotonically, which is the interlaced-scan
     answer this document called for.
   - The `_BoundedGate` bounds in-flight *groups*, not frames. A reader
     loads a whole group before queueing it, so a per-frame gate would let
     several readers each hold a partial group with none able to finish.
   - `_frame_parallelism` charges a worker F x the native tile working set
     (mirroring the kernel's own precheck, which is written against
     `frames * tile_pixels`) and F x the corrected-frame buffer, but only
     1 x the record term — a group emits one merged batch, not F of them.

   **`F = 1` is not bit-for-bit with the old per-image path**, as step 1
   already established. Same voxels, same contributor counts, totals equal
   to rounding. That is a `phys`-scope change to every job, not only to
   jobs that opt into grouping.
5. *Scheduler.* **Done**, as `_map_frame_groups_streamed`, and it makes
   grouping pay: **0.866x at F = 4** on the real job, measured against the
   per-frame pipeline. One correction moved, one premise of this document
   fell over, and the second is the more useful result.

   **What worked: hoisting correction out of the compute worker.** Loading
   and correcting a frame is GIL-held Python; mapping is not. Correction
   now runs in a prepare pool, so group N+1's Python work overlaps group
   N's native call instead of queueing behind it. Prepare workers
   cooperate on *one* group at a time rather than each taking a group of
   their own: a group buffer is F full-detector corrected frames, the
   largest resident term in the pipeline, so per-worker groups would
   multiply it by the pool size. Frame granularity keeps live buffers down
   to the pipeline depth while still spreading a group's F loads and
   corrections across the whole pool.

   **What was wrong: "one group call saturates the thread budget".** This
   document argued that a group call has ~6000 bricks and so needs no
   help, which makes `kernel_threads = spec.threads` with one call in
   flight the natural layout. Implemented exactly that way it measured
   *1.38x worse* than the old scheduler at F = 4 — 342.4 ms/frame against
   248.1. Brick count is not the constraint; whatever the call does around
   its parallel region is. Three concurrent calls of eight threads beat
   one call of twenty-four, and `_group_pipeline_layout` now buys
   concurrency from memory much as `_frame_parallelism` did, the
   difference being that a group buffer is shared by the prepare pool
   rather than owned by one frame worker.

   Interleaved repeats, 234 frames, the job's own budget and tiling:

   | arm | layout | ms/frame | vs F = 1 |
   |---|---|---|---|
   | F = 1 (per-frame pipeline) | 8 workers x 1 | 232.8 | 1.000 |
   | **F = 4** | **3 calls x 8 threads** | **201.5** | **0.866** |
   | F = 4, 256-row bands | 4 calls x 6 threads | 205.8 | 0.884 |
   | F = 8, 256-row bands | 2 calls x 12 threads | 236.1 | 1.014 |
   | F = 8, the job's own bands | 1 call x 24 threads | 288.5 | 1.24 |

   Every F = 4 run was faster than every F = 1 run, so the separation
   survives the noise even though each arm spreads about +-7%.

   **F = 4 is the optimum here because the two effects cross.** Records
   fall monotonically with F (0.776 / 0.657 / 0.598 at F = 2 / 4 / 8) but
   concurrency falls with it too, since both the group buffer and the
   call's native working set scale with F. At F = 8 the job's own tiling
   affords one concurrent call, and that costs more than the extra 6% of
   records is worth. Banding more finely buys the concurrency back — 256-row
   bands restore two calls and take F = 8 from 1.24x to 1.014x — but not
   enough to beat F = 4. **Band height and group size trade against each
   other for the same memory**, which is a knob this document had not
   identified, and picking the pair rather than picking F alone is what
   step 6 should actually do.
6. *Choosing F per job.* **Done.** `spec.frames_per_group` now defaults to
   `None`, meaning "measure this job and choose"; an explicit integer still
   overrides. Two gates, and the second is the one that is easy to get
   wrong.

   **The regime gate.** `_frame_advance_voxels` asks the kernel where a
   sample of pixels lands on one frame and on the next, and takes the
   difference in units of the grid step. Geometry only — no image data and
   no I/O — so it costs a few hundred `coordinate()` calls at job start.
   Past `_GROUP_ADVANCE_VOXEL_LIMIT` (one voxel) consecutive frames tile
   rather than overlap and the answer is `F = 1`. Sampled at three places
   in the scan, not one, so a non-uniform angular step is caught. The
   monotonicity check for interlaced scans gates both this and any
   explicit request.

   **The concurrency gate, which is what density alone would get wrong.**
   Records fall monotonically with F, so choosing on density picks the
   largest group that fits memory — which on this job is the configuration
   that maps *slowest*. The rule is instead the largest F that still leaves
   `_GROUP_MINIMUM_CONCURRENT_CALLS` (3) native calls affordable.

   On the reference job this picks **F = 4 with 3 concurrent calls of 8
   threads**, on every run, which is the configuration the sweep in step 5
   found by hand. Interleaved against an explicit `F = 1`, three pairs:

   | pair | automatic | F = 1 | ratio |
   |---|---|---|---|
   | 1 | 262.3 | 311.2 | 0.843 |
   | 2 | 223.4 | 255.6 | 0.874 |
   | 3 | 250.3 | 268.5 | 0.932 |

   **0.88x**, with records at 0.657x. Note both arms are slower in absolute
   terms than the step 5 table taken a day earlier (278 ms/frame for `F = 1`
   against 233); the machine drifts between sessions, and only the paired
   ratio is stable. That is the whole reason the arms are interleaved.

   One case worth knowing about, pinned by a test: a **lab-frame grid**
   measures zero travel, because Q in the lab frame does not depend on the
   sample angles at all. Every frame reaches the same voxels, so grouping
   helps most there and the gate correctly does not stand in the way.

7. *Choosing the band height too.* **Measured, and there is nothing to
   choose.** Step 5 showed band height and group size trade against each
   other for the same memory, so re-tiling to afford a larger group looked
   like the obvious next lever. It is not: swept at `F = 4` on the
   reference job, band height is flat across the region that matters and
   the existing default sits in the middle of it.

   | bands | tiles | layout | median ms/frame |
   |---|---|---|---|
   | 240 | 11 | 4 x 6 | 236.7 |
   | **421 (the default)** | 7 | **3 x 8** | **205.4** |
   | 506 | 5 | 3 x 8 | 207.1 |
   | 600 | 4 | 2 x 12 | 237.1 |

   Both directions lose, and for different reasons. Finer bands buy a
   fourth concurrent call but take two threads off each, and 4 x 6 is
   worse than 3 x 8 — more calls, each with less work to spread. Coarser
   bands cost a call outright. **`F = 8` cannot be rescued this way at
   all**: its group buffer does not shrink with the band height, so it
   holds at two concurrent calls from 200 rows to 2527 and never reaches
   the floor of three.

   So the pair does not need choosing jointly; the band height that
   maximises concurrency at the chosen `F` is already what the memory
   budget produces.

   **One real defect surfaced while sweeping it, and is fixed.** Tile
   planning still sized bands by `8**depth` bytes per pixel — the
   worst-case leaf count phase 2 replaced in both memory prechecks and
   left behind at this third site. It is wildly conservative: at depth 2
   it claims 5248 bytes for a pixel that costs about 106, banding this
   detector into 13 strips of 194 rows, and at depth 5 it asks for 2.6 MB
   per pixel and collapses a band to **one row — 2527 native calls per
   frame instead of six**. Bounding it by the same record ceiling the
   arena and both prechecks use makes the band height depth-independent.

   Measured at depth 2, interleaved, the fix is a **wash** (619 against
   629 ms/frame, 0.98x) — the 5-20% predicted from the depth-0 sweep did
   not appear. The reason is visible in the run: wider bands make the
   group-size gate drop `F` from 4 to 2, and the two effects cancel almost
   exactly. That is worth stating as a result in its own right — **band
   height and group size substitute for each other**, which is why neither
   is worth optimising hard, and why the joint search this step set out to
   build would have found a flat surface.

   The fix is kept for the depth-5 pathology and for consistency with
   phase 2, not for throughput. A test pins the band count against depth.

### Measurement noise on this machine

See "Measuring any of this" above: single runs vary by up to 30%, arms must
be interleaved, and only paired ratios mean anything. The step 5 table was
taken that way. Two earlier conclusions were not, and were withdrawn.

### What the wired pipeline measured

`benchmarks/benchmark_reconstruction_pipeline.py` runs the real
`_map_pending_ranges` — reader pool, correction, kernel, router, checkpoint
writes — over a bounded window and reports both throughput and the records
*entering* the checkpoint layer. 234 frames of the `39_1-rsmap` job (six
scheduling ranges, three checkpoints, depth 0, stationary bounds, 24
threads, the job's own 10 GB budget), against a baseline captured from the
same checkout immediately before the change:

| F | records/frame | vs F=1 | ms/frame | vs baseline | workers |
|---|---|---|---|---|---|
| baseline (`accumulate`) | 2,584,531 | 1.000 | 236.8 | 1.00 | 8 |
| 1 | 2,584,531 | 1.000 | 263.1 | 1.11 | 8 |
| 2 | 2,005,606 | 0.776 | 245.8 | 1.04 | 6 |
| 4 | 1,699,144 | 0.657 | 248.1 | 1.05 | 3 |
| 8 | 1,545,463 | **0.598** | 347.4 | **1.47** | **1** |

**The record prediction held.** 0.598 at F = 8 against the 0.575 the kernel
benchmark measured on a single tile and the 0.5749 the exhaustive geometric
model predicted before any of it existed. Finding A survives contact with
the whole pipeline.

**The throughput prediction did not**, on the per-frame scheduler. The
band stated before implementing was 0.7-1.0x; F = 8 landed at 1.47x. The
cause is visible in the last column: `_frame_parallelism` charges a worker
for F resident frames, so `image_workers` collapses from 8 to 1, and at one
worker a group's GIL-held correction has nothing to overlap with. This is
the synchronisation point this document flagged as the main scheduling
risk, arriving where it was predicted to — and it is what step 5 below
fixes, after which F = 4 reaches 0.866x. Read the two together: on this
scheduler the wiring alone was not worth turning on, and
`frames_per_group` stayed at 1 until step 5 landed.

The individual figures in the table above are single runs, taken before
the run-to-run spread on this machine was understood, and should be read
as the shape of the effect rather than as values — see "Measurement noise"
under step 5. The record column is unaffected: it is deterministic.

Two smaller things the measurement settled. Records *written* to checkpoint
parts barely move with F (0.97x at F = 8) because the checkpoint
accumulator already merges across frames — it collapses the same
redundancy later and more expensively, so the saving never appears there.
Records entering the router is the quantity that matters, and it is what
the table above reports. And F = 1 through `accumulate_group` routes
*exactly* the same record count as `accumulate` (604,780,339, to the
record): the 3.4% block-shape gain is real in block records before the
cross-block merge, but one call's merged output is the set of voxels that
frame reached, which does not depend on how the frame was partitioned.

### Phase 4 — Frame-group bricks above depth 0. Not started.

Everything in phase 3 was designed and measured at depth 0. Grouping is not
*disabled* above it — the gates run at any depth and the streamed scheduler
works — but nothing about it has been tuned there, and the little evidence
that exists is not encouraging.

**What is already known, from the depth-2 runs taken while fixing band
height.** At depth 2 on the reference job the group-size gate picks
`F = 2`, not the `F = 4` it picks at depth 0, and mapping runs at
~620 ms/frame against ~205 at depth 0. The reason `F` falls is structural:
a group call's native working set is `F x tile pixels x bytes-per-pixel`,
and the per-pixel term roughly quadruples from depth 0 to depth 2, so the
same memory affords half the group at the same concurrency. Grouping and
adaptive depth compete for the same budget.

That was measured incidentally, on two arms that differed in band height as
well, so treat it as an observation rather than a result. **The first real
phase 4 measurement is the obvious one nobody has taken: `F = 1` against
automatic `F` at depth 1 and depth 2, interleaved, three repeats.** Until
that exists there is no evidence that grouping helps or hurts above depth 0,
and the honest reading of the numbers above is that it does less.

Only after that is the design work worth doing:

- **Pixel-outer / frame-inner iteration.** Above depth 0 the per-pixel
  subdivision state cannot be shared across frames (see "Where the two
  meet"), so the brick loop needs its second mode. Corner-passing is what
  makes this affordable at all — 1.2–3.5 KiB of stack instead of a
  6.9–101 KiB dense array per pixel.
- **No cross-frame ray cache.** The payoff is capped near 28% and the
  cache is per pixel, so it can only live for one pixel at a time.
- **Brick dimensions re-derived** against the phase 2 record bound rather
  than the leaf count.
- **Whether the record saving even survives.** At higher depth a pixel
  resolves into a footprint rather than a point, so cross-frame overlap
  should be at least as good — but "should" is doing all the work in that
  sentence, and the records-converge-at-depth-3 observation below suggests
  the marginal voxel is increasingly a partial-volume weight rather than a
  new voxel. Measure the routed-record ratio at depth 1–3 before assuming
  finding A transfers.

There is also a prior question worth settling first, because it may make
the whole of phase 4 moot: **records converge at depth 3** (0.632 → 0.650
across depths 3–6 while cost rises 3x / 6x / 13x). If `very_high` and
`maximum` do not change a reconstructed intensity measurably, the right
answer is to discourage them rather than to optimise them.

## Risks and open questions

**Interlaced scans.** `orgui/backend/interlacedScanLoader.py` exists, so frames
adjacent in file order are not always adjacent in angle. Groups must be formed by
angle adjacency; where that conflicts with sequential I/O or the resume
frame-range structure, the honest answer is `F = 1`. Check monotonicity of
`exposure_angle_bounds` per job rather than assuming it.

**Non-uniform angular steps.** The same problem in milder form; the per-job probe
handles it if it samples more than one place in the scan.

**Bit-for-bit.** Phase 1 must be bit-identical and should be held to that — a
discrepancy means the lattice or the visit order diverged, and both are bugs.
Phase 3 is *not* bit-identical: contributions from several frames merge in the
block map instead of in the tree accumulator. Identical arithmetic, fewer
intermediate roundings, so slightly better conditioned, but any test pinning an
XXH3 digest of record values will move. That is a `phys`-scope change and needs
saying so.

**The group synchronisation point**, as described under Scheduling. *Real, and
resolved by moving correction into the prepare pool (phase 3 step 5).* The
preference for F = 8 was wrong on this job: F = 4 measured better, because
concurrency is bought from the same memory the group buffer spends.

**The Scheduling section above is superseded in one respect.** Its claim that a
group call saturates the thread budget on its own, and that `kernel_threads =
spec.threads` with one group in flight is therefore the natural layout, was
implemented and measured 1.38x worse than keeping several concurrent calls. It
is retained as the reasoning that led to the experiment; phase 3 step 5 records
what replaced it. The prediction that the pool would need about five workers was
sound in spirit — the prepare pool grows to its ceiling at F = 4 — but
`image_workers` does not dissolve, it becomes concurrent group calls.

**MSVC inlining** of the recursion is an assumption; check the measurement, not
the intent. The 19-point octree bookkeeping is the fiddly part.

**Depth > 0 for finding A** is unmeasured, as is the interaction between swept
exposures and cross-frame overlap above depth 0 — each frame then covers a
0.71-voxel slab rather than a plane, and adjacent slabs tile rather than overlap.

**`excluded_frames`** must be skipped when forming groups, not merely masked
inside them. *Resolved, and more cheaply than expected:* `_included_ranges`
already splits scheduling ranges at every excluded frame, so a range is
contiguous and exclusion-free by construction and `_frame_groups` only has
to cut at checkpoint boundaries.

**Records converge at depth 3.** 0.632 -> 0.642 -> 0.647 -> 0.650 across depths
3-6, while cost rises 3x / 6x / 13x. Everything above depth 3 buys finer
partial-volume weights, not new voxels. Whether `very_high` and `maximum` change
a reconstructed intensity measurably is a separate question, but it is worth
asking before optimising them — it bears directly on how much phases 1 and 2 are
worth.

## Also

The benchmark mechanics that used to live here — which script to run, the
`PYTHONPATH` and schema-version traps, the noise discipline, why routed
records rather than rows on disk, and the contributor-weighted key sum used
where bit-for-bit is impossible — are now collected under "Measuring any of
this" above, before the plan rather than after it.

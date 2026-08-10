# Reciprocal-space mapping: locality, subdivision cost, and the scheduling that would carry them

> **Status: idea for later. Nothing here is implemented, and no decision has
> been taken.** The measurements are real and were taken against real beamtime
> data; the designs sketched from them are proposals. Treat this as a record of
> an investigation, not a plan of record.

Third companion to `reciprocal_space_scratch_architecture.md` (what the pipeline
*is*) and `reciprocal_space_performance_findings.md` (what has been measured
about it). This one holds two findings that neither of those settles, the
places where they meet, and what the scheduler and the memory budget would have
to do to carry them.

- **A — cross-frame locality.** The native work block deduplicates voxels within
  one image, mostly by accident. Extending it across adjacent images cuts
  emitted records by ~42% and *shrinks* the cache working set. **Not started.**
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

### Phase 2 — The `worst_leaves` bound (memory, standalone)

Replace `children**depth` in both mirrored sites with a bound derived from
boundary growth (~2^depth) plus headroom, keeping the two in step. Findings open
item 2 argues the authoritative cap belongs in the kernel, which knows
`stationary` at `accumulate()` time while the Python side must assume the worst.
Unblocks realistic tile sizes — and therefore `image_workers` — from depth 3 up.
Not required by either finding; do not fold it in.

### Phase 3 — Frame-group bricks at depth 0 (finding A)

1. *Kernel.* `accumulate_group()` taking `(frames, rows, columns)` arrays and
   `(frames, 2, 4)` angle bounds, returning the same six columnar arrays merged
   across the group. The block queue becomes a 3D brick index; `block_results`,
   the per-worker arena and the loser-tree merge are unchanged.
2. *Brick loop, depth 0.* Precompute the brick's unit centre rays once into a
   24 KiB stack buffer, iterate frames-outer / pixels-inner. Note separately
   that at depth 0 the centre ray is frame-invariant for the whole job — pure
   detector geometry, currently recomputed 3651 times per pixel. A
   whole-detector centre-ray cache (149 MB, alongside the existing 143 MiB
   corner cache) removes it entirely, with or without grouping.
3. *Banding.* Bricks iterate row band by row band so level-3 records stay near
   today's ~110 MB, as above.
4. *Choosing F per job.* Extend the calibration probe: run `accumulate_group` on
   one sample brick at F in {1, 2, 4, 8, 16}, take the smallest F past which the
   marginal density gain falls below a few per cent, capped by memory. Five small
   native calls, and it measures the quantity that actually matters rather than
   deriving it from the angular step. `F = 1` must reproduce today exactly.
5. *Scheduler.* As above: pin `kernel_threads`, size the correct pool from the
   measured Python-to-native ratio, retarget `_BoundedGate`, repoint the live
   rebalance at `F`.

### Phase 4 — Frame-group bricks above depth 0

Only after phases 1 and 3. Pixel-outer / frame-inner, no cross-frame ray cache,
brick dimensions re-derived against the corrected leaf bound from phase 2. All
of finding A is measured at depth 0 only; at higher depth a pixel resolves into
a footprint rather than a point, so overlap should be at least as good, but that
is an expectation, not a result.

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

**The group synchronisation point**, as described under Scheduling. Prefer F = 8
until measured.

**MSVC inlining** of the recursion is an assumption; check the measurement, not
the intent. The 19-point octree bookkeeping is the fiddly part.

**Depth > 0 for finding A** is unmeasured, as is the interaction between swept
exposures and cross-frame overlap above depth 0 — each frame then covers a
0.71-voxel slab rather than a plane, and adjacent slabs tile rather than overlap.

**`excluded_frames`** must be skipped when forming groups, not merely masked
inside them.

**Records converge at depth 3.** 0.632 -> 0.642 -> 0.647 -> 0.650 across depths
3-6, while cost rises 3x / 6x / 13x. Everything above depth 3 buys finer
partial-volume weights, not new voxels. Whether `very_high` and `maximum` change
a reconstructed intensity measurably is a separate question, but it is worth
asking before optimising them — it bears directly on how much phases 1 and 2 are
worth.

## Also

`benchmarks/benchmark_reconstruction_mapping.py` no longer imports: it wants
`_map_frame_range`, which was replaced by `_map_one_frame`.

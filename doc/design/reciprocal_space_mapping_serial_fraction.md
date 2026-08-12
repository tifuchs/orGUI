# Reciprocal-space mapping: the per-frame serial fraction

> **Status: planned, nothing implemented.** This is the next feature in
> the mapping phase. It exists because phase 4 of
> `reciprocal_space_mapping_locality_and_subdivision.md` measured
> something it was not looking for: at low adaptive depth the pipeline
> leaves two thirds of the machine idle, and the kernel is not what is
> holding it back.

Fourth document in the set, after `reciprocal_space_scratch_architecture.md`
(what the pipeline is), `reciprocal_space_performance_findings.md` (what has
been measured), and `reciprocal_space_mapping_locality_and_subdivision.md`
(the locality and subdivision work, complete). This one takes over findings
open item 1, which has been the largest structural limit in that document
since it was written and is now the largest *measured* one.

## The finding

Sampled per process during live runs on the reference job at depth 0:

| | |
|---|---|
| wall time | ~110 ms/frame |
| cores in use | **7-8 of 24** |
| cold I/O floor | 51.5 ms/frame (49.8 MB at 967 MB/s) |

Mapping is not the constraint in this regime. Loading and correcting are,
and they hold the GIL, so adding workers does not help. The old Amdahl
estimate — 10.5% serial fraction, 6.1x ceiling at 16 workers, 4.2x
measured — understated how lopsided this becomes once the kernel is cheap.

**The headroom is about 2.1x, not more.** An earlier draft of this claim
said 3x, reasoning from idle cores alone; that is wrong, because cold
reads cap the pipeline at 51.5 ms/frame however free the CPU work
becomes. 110 to 51.5 is 2.14x and that is the whole prize. It shrinks
further whenever the page cache holds the window, which on this machine
it partly does — see the depth-0 noise trap in the findings document.

## Where this pays, and where it does not

Grouping (phase 3 of the locality document) pays where the kernel
dominates: depth 0 on a grid fine enough to make record density the cost.
This work pays where the kernel is cheap: depth 0 on a coarse grid, which
is `center` accuracy on a routine job — the common case, and the one that
should be fastest. Above depth 1 neither matters: the subdivision walk
costs 600-2300 ms/frame and a 60 ms serial fraction disappears beside it.

So the two are complementary rather than competing, and both are depth-0
work. That is worth knowing before optimising either further: **the whole
remaining opportunity in mapping is at depth 0**, because depths above 2
change no intensity a measurement can resolve.

## What has to be measured first

The stage table in the findings document (load 48.3 / correct 42.3 /
slice-copy 24.0 / kernel 853.4 / reduce 121.2 ms) predates grouping, the
prepare pool and the current grid, and was taken in a kernel-dominated
regime. **It cannot be used to apportion 110 ms/frame today.** Re-measure
before designing anything:

1. Per-stage wall and CPU time at depth 0 on both grids, inside the real
   scheduler rather than a harness, distinguishing
   - HDF5 read and decompression (GIL released inside h5py),
   - correction arithmetic (several full-detector NumPy passes),
   - the per-tile `np.stack` / `ascontiguousarray` copy,
   - time blocked on the GIL itself.
2. How much of `load` is genuinely I/O-bound. If it is already at the
   cold-read floor, correction is the only thing left worth attacking.

The GIL-blocked fraction is the number that decides this. If it is small,
the pipeline is I/O-bound, the ceiling is nearer 1.2x than 2.1x, and this
whole document should be closed the way phase 4 was.

## Candidate mechanisms

In increasing order of invasiveness. Each carries its own predicted band,
stated before implementing so it can be wrong on the record.

### A. Drop the per-tile stack and copy

`_map_frame_group` builds `np.ascontiguousarray(np.stack(...))` for
intensity, variance and mask, per tile, per group. Tiles partition the
detector, so this copies every corrected frame exactly once — 105 MB per
frame — in Python, under the GIL, purely to hand the kernel a contiguous
`(frames, rows, columns)` buffer.

`accumulate_group` could instead take the whole corrected frame plus the
tile rectangle and read the sub-rectangle itself; it already walks bricks
within a tile, so a row stride costs it nothing.

**Predicted: 0.80-0.90x** at depth 0, from removing ~24 ms of the ~110.
Findings open item 1 named this direction; it is the cheapest of the
three and does not touch any arithmetic, so it should be bit-for-bit.

### B. Fuse and release the correction

Correction is several full-detector NumPy passes — background subtract,
static factor, exposure and monitor scaling, the finiteness check, mask
combination. Each is memory-bound over ~50 MB, and NumPy holds the GIL
across the small ones. Fusing them into a single pass in the existing C++
extension would both cut memory traffic and release the GIL, letting the
prepare pool actually run in parallel.

**Predicted: 0.65-0.80x** on top of A. **This is a `phys` change** — it
moves arithmetic that produces the corrected intensities and variances
into another language and another association order. It needs the same
treatment phase 3 got: voxels, contributor counts and totals pinned, and
the last-bit movement stated. Do not fold it in with A.

### C. Extend the prepare pool to the per-frame path

The grouped scheduler already hoists load and correct into a prepare
pool; the per-frame scheduler does not, so a job outside the grouping
regime — a coarse scan, an interlaced one — still queues its GIL-held
Python behind its own native call. With B done, this is mostly
bookkeeping.

**Predicted: little on the reference job** (which groups), and up to the
same 0.87x grouping bought, on jobs that do not.

## Ordering

A, then measure; B, then measure; C last. A is bit-for-bit and cheap. B
is where the prize is and where the risk is. C only pays for jobs the
first two do not already cover.

Do not do all three and measure once. Phase 4 is the standing argument:
a saving of known size need not convert, and finding that out one change
at a time is what makes it interpretable.

## Risks

**The ceiling is low and the noise is high.** 2.1x is the best case, and
depth 0 is the regime the findings document flags as dominated by page
cache state — one arm spread 176.5 to 340.6 ms/frame. Six interleaved
pairs minimum, and the foreign-load guard on every run. It is entirely
possible for a real 15% win to be invisible here.

**B changes reconstructed values.** Not in which voxels are reached, but
in the last bits of the corrected intensities that feed them. That is
acceptable on the same grounds the phase 3 change was — the feature has
never been released — but it must be said, and it must not be discovered
by a user.

**I/O may already be the wall.** If step 1 of the measurement shows load
sitting at the cold-read floor with correction overlapped, there is
nothing here worth building, and saying so is the correct outcome.

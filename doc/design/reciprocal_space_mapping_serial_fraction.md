# Reciprocal-space mapping: the per-frame serial fraction

> **Status: measured and complete. A and B are implemented (0.84x and
> 0.92x, compounding to about 0.77x at depth 0); C was built, measured as
> a regression, and reverted.** The premise the document was written on —
> that the pipeline is held back by GIL-serialised per-frame Python — was
> measured false at step 1, and what replaced it is that this machine is
> short of memory bandwidth rather than of cores. Both changes that
> survived are traffic reductions, and the one that was purely a
> scheduling change is the one that failed. B's bit-for-bit property is
> qualified: it holds under `-ffp-contract=off` (the `strict_fp_contract`
> build option, which the test workflow sets and release wheels do not),
> and to within one ULP of the variance on compilers that fuse
> multiply-add. See B.
>
> It exists because phase 4 of
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

## What had to be measured first

The stage table in the findings document (load 48.3 / correct 42.3 /
slice-copy 24.0 / kernel 853.4 / reduce 121.2 ms) predates grouping, the
prepare pool and the current grid, and was taken in a kernel-dominated
regime. **It cannot be used to apportion 110 ms/frame today.** So the
first step was to re-measure, inside the real scheduler rather than a
harness, distinguishing HDF5 read and decompression, correction
arithmetic, the per-tile copy, the native call, the cross-tile merge and
the checkpoint route — and, for each, wall time against time actually
spent executing.

**The GIL-blocked fraction was named as the number that decides this: if
it were small, the pipeline is bound by something this document cannot
move and the whole thing should be closed the way phase 4 was.** It was
measured. It is small. What the same measurement turned up instead is
that the largest per-frame Python cost is not where this document assumed
it was, so the ordering below survives while the reasoning under it does
not. Both are recorded next.

## Step 1 — where the time actually goes

Instrument: `benchmarks/benchmark_reconstruction_stages.py`, which runs
the real `_map_pending_ranges` and wraps each stage's existing call site
from outside — no instrumentation compiled into the pipeline, no
arithmetic touched. On-core time is `QueryThreadCycleTime`, not
`time.thread_time`: the latter is `GetThreadTimes` on Windows and
quantises to the ~15.6 ms scheduler tick, which is larger than most of
these stages.

Reference job as it stands on disk (1000^3 grid), `--depth 0`, 234
frames, the scheduler's own choice of four frames per call and three
concurrent calls of eight threads. **Milliseconds per frame, summed over
every thread** — so a column may exceed the 111.5 ms/frame the run took.

| stage | wall | on-core | blocked |
|---|---|---|---|
| load (`get_raw_img`) | 100.2 | 63.4 | 36.8 |
| correct (`correct_frame`) | 132.8 | 89.3 | 43.5 |
| **`_map_frame_group`** | **257.9** | **114.3** | **143.5** |
| — per-tile assembly | 88.9 | 43.0 | 45.9 |
| — of which `np.stack` | 70.3 | 27.7 | 42.6 |
| — the native call | 103.0 | 22.4 | 80.6 |
| — cross-tile merge | 21.6 | 15.2 | 6.4 |
| — checkpoint route | 44.4 | 33.8 | 10.6 |

Whole run: **111.5 ms/frame, 6.40 of 24 cores, 713 ms of CPU per frame**,
of which 267 ms is on Python threads and the remaining ~446 ms is the
kernel's own native worker pool. (The native call's own on-core column is
the *calling* thread's only; its 80.6 ms of "blocked" is a thread waiting
on its native workers, not on the GIL.)

The same measurement on the documented 2000^3 grid, where the kernel
dominates: 280.0 ms/frame, 6.19 cores, 1734 ms of CPU per frame, and the
route becomes the largest single term (262.1 wall / 232.5 on-core) as
records rise 5.2x. Per-tile assembly is essentially unchanged in absolute
terms (95.4 / 53.8 / 41.7), as it should be — it is a function of the
detector, not of the grid.

### The GIL is not the constraint

A probe thread runs bursts of pure Python — no I/O, nothing that could
release the GIL voluntarily — several times the interpreter's switch
interval, and compares each burst's wall time against the cycles it
executed. That is the blocked fraction a Python thread wanting to run
actually experiences:

| | |
|---|---|
| idle machine | **1%** |
| during the run, coarse grid | **23%** |
| during the run, fine grid | **22%** |

So a Python thread gets the GIL roughly three times out of four, and the
total GIL-held work is on the order of 25-30 ms per frame against 713 ms
of CPU per frame. **The pipeline is not GIL-bound**, and the premise this
document was written on — "loading and correcting hold the GIL, so adding
workers does not help" — is wrong as stated.

### Correction is already overlapped

The prepare pool the grouped scheduler introduced does what phase 3
claimed for it, and more completely than anyone checked. Load plus
correct is 233 ms/frame of thread-wall spread over the 14-16 prepare
workers the pool grows to: **about 13% occupancy**. The compute worker,
by contrast, is busy 257.9 / 3 workers / 111.5 ms = **77% of the time**.
The critical path is the compute worker, and correction is not on it.

### What *is* expensive is doing the same arithmetic under load

| stage | uncontended, one thread | in the pipeline | inflation |
|---|---|---|---|
| load | 32.1 ms/frame on-core | 63.4 | 2.0x |
| correct | 38.5 | 89.3 | 2.3x |

Same code, same frames, twice the cycles. Cycles are counted while the
thread is on a core, so this is not GIL waiting — it is memory: 24
threads sharing one L3 and one controller. That matters twice over. It
says the machine's spare cores are worth less than their count suggests,
and it says **removing memory traffic is worth more than removing
instructions**.

### What this retires, and what survives

- **Closed: the GIL argument.** 23% is not nothing, but it is not what
  holds the pipeline at 6.4 cores, and no amount of releasing it converts
  into the 2.1x this document opened with.
- **Closed: "the prepare pool cannot run in parallel".** It can, and does;
  it is idle 87% of the time.
- **Survives, for a different reason: A.** The per-tile assembly is 88.9
  ms/frame of wall on the one thread that is 77% busy — 34% of that
  thread's occupied time, and the largest thing on the critical path
  after the native call itself. It is not a GIL problem; it is 105 MB per
  frame of pointless copying on the busiest thread, and 210 MB/frame of
  memory traffic on a machine that is measurably short of bandwidth.
- **In doubt: B.** Its stated mechanism (release the GIL so the prepare
  pool can run) is measured false. Its unstated one (cut memory traffic)
  survives, but the pool it would speed up is 13% occupied. See B below,
  which is reconsidered against this rather than against the plan.
- **The 2.1x headroom is not reachable by this route.** Reaching the
  51.5 ms/frame cold-read floor needs the CPU work to become nearly free,
  and 713 ms of CPU per frame at 2x inflation is not nearly free.

## Candidate mechanisms

In increasing order of invasiveness. Each carries its own predicted band,
stated before implementing so it can be wrong on the record.

### A. Drop the per-tile stack and copy — **done, 0.840x**

`_map_frame_group` built `np.ascontiguousarray(np.stack(...))` for
intensity, variance and mask, per tile, per group. Tiles partition the
detector, so this copied every corrected frame exactly once — 105 MB per
frame — in Python, on the compute worker, purely to hand the kernel a
contiguous `(frames, rows, columns)` buffer.

`accumulate_group_tile` now takes the whole corrected frames plus the
tile rectangle and reads the sub-rectangle in place. The brick loop
already walked rows within a tile, so the change is one stride: a
`GroupPixels` view carries one pointer per frame at the tile's first
pixel and the distance between rows in the caller's own array, and
`accumulate_group` builds the same view over its contiguous buffer. One
loop, two ways in, and on the reference job's full-width bands the row
stride is *identical* either way — the tile is the frame's whole width,
so only the base offset differs.

**Predicted 0.80-0.90x** by this document, and **0.72-0.90x** from the
step-1 stage table, which put the assembly at 34% of the busy time of the
one thread that is 77% busy. Measured, six interleaved pairs, 234 frames,
depth 0, coarse grid:

| pair | baseline | tile view | ratio |
|---|---|---|---|
| 1 | 106.5 | 92.2 | 0.866 |
| 2 | 105.8 | 91.4 | 0.864 |
| 3 | 108.4 | 91.0 | 0.840 |
| 4 | 103.1 | 85.7 | 0.831 |
| 5 | 103.5 | 87.0 | 0.841 |
| 6 | 102.2 | 85.0 | 0.832 |

**0.840x, and every tile-view run beat every baseline run** (85.0-92.2
against 102.2-108.4). Foreign load 2.16-2.71 cores throughout, inside
the machine's usual band.

**Bit-for-bit, as required.** Identical routed record count, identical
checkpoint rows, contributors and voxel fingerprint, and identical
totals to the last bit (`weight`, `weighted_intensity` and
`weighted_variance` compared by `repr`) over a 78-frame window; the
kernel-level test pins the same thing directly, against the copied form,
at depths 0-2, stationary and moving, for a tile in the middle of the
frame as well as one at its edge.

Re-measured with the stage instrument afterwards, the pipeline moved as
the model said it would: 111.5 -> 91.4 ms/frame, cores in use 6.40 ->
8.40, compute-worker occupancy 77% -> 71%, GIL contention 23% -> 19%.
Python-thread CPU fell 267 -> 246 ms/frame. Note that *process* CPU per
frame did not fall — it rose slightly, to 768 ms — which is the memory
story again: at higher concurrency every remaining stage costs more
cycles for the same work.

One thing this did not need: a memory-model change. `_frame_parallelism`
never charged for the per-tile buffer (it charges the whole detector's
corrected arrays per in-flight frame, which is 6x larger), so removing it
makes the existing estimate slightly more conservative rather than
wrong.

### B. Fuse the correction

Correction is several full-detector NumPy passes — background subtract,
static factor, exposure and monitor scaling, the finiteness check, mask
combination. Each is memory-bound over ~50 MB, and NumPy takes and drops
the GIL between them.

**The stated mechanism was wrong.** This document predicted 0.65-0.80x
from "release the GIL, letting the prepare pool actually run in
parallel". Step 1 measured the GIL at 23% and the prepare pool at 13%
occupancy: the pool already runs in parallel, and there is no serialised
correction to release.

**The unstated mechanism survives, and was measured before building
anything.** A ceiling probe — `--correction minimal` in the stage
benchmark, which keeps the mask and the finiteness check and drops the
scaling arithmetic entirely, so it routes the identical records and only
carries wrong intensities — bounds from above what any fusion could buy:

| pair | full | no scaling | ratio |
|---|---|---|---|
| 1 | 98.2 | 83.7 | 0.852 |
| 2 | 91.0 | 83.6 | 0.918 |
| 3 | 90.9 | 85.1 | 0.937 |
| 4 | 85.9 | 83.9 | 0.976 |

**0.928x, clean separation.** The prediction stated before running that
probe was 0.95-1.00x — *that is also wrong*, in the other direction. Not
doing ~700 MB/frame of memory traffic is worth about 7%, on a pipeline
whose prepare pool is idle 87% of the time. Which is the step-1 result
again from a different angle: this machine is short of memory bandwidth,
not of cores, so **work removed anywhere helps, including off the
critical path.**

That is what justified building B after its own premise had failed, and
it is why what was built is narrower than what was proposed: one native
pass for the arithmetic that costs, and not a native rewrite of a
correction step that is mostly allocation.

**What was built.** `apply_correction_factors` applies the per-pixel
static factor, then each scalar factor in turn, then folds the finiteness
check into the mask — one pass, one pixel at a time, with the GIL
released for its whole duration. What stayed in NumPy is the part that is
allocation rather than arithmetic (the raw copy, the variance
initialisation, the background subtraction) and the repair plan, which is
already native and has to run between the two.

**It is bit-for-bit, which the plan assumed it could not be.** Correction
is entirely element-wise: there is no reduction, so there is nothing to
reassociate, and doing the same operations on the same values in the same
order gives the same bits. Two details make that true rather than nearly
true. `factor ** 2` is `pow` in Python and a multiplication in C++, and
those are not guaranteed to agree in the last bit, so the caller passes
the squared factor it already computed rather than letting the kernel
recompute it. And the propagated-uncertainty branch scales the variance,
*then* uses the intensity from before the factor was applied, *then*
scales the intensity — an order an independent implementation would very
plausibly get wrong, so the test exercises both branches. The end-to-end
run agrees on voxels, contributors and all three totals by `repr`.

**It holds only where the compiler does not contract, and release wheels
are built where it does.** *Found 2026-08, when the app suite first ran
on macOS in CI.* The propagated branch computes `spread += value * value
* factor_variance`, and clang defaults to `-ffp-contract=fast`, so on
arm64 that becomes a single `fma`: one rounding where NumPy does two. The
result differs by **one unit in the last place** — 19 of 99 elements, max
relative difference 2.36e-16 against a double epsilon of 2.22e-16, max
absolute 5.4e-20 on values of ~1e-4. Only variance moves; intensity is
pure multiplication and stayed identical. The FMA result is the *more*
accurate of the two, so this is a code-generation choice, not an error.

The `strict_fp_contract` Meson option (`-ffp-contract=off`, default
**false**) exists so the test suite can hold the invariant, and
`tests.yml` sets it. Release wheels are deliberately built without it, to
keep the compiler's own code generation. So the claim above should be
read as: **bit-for-bit under `-ffp-contract=off`, and to within one ULP
of the variance otherwise.** Anyone comparing a shipped macOS arm64 build
against the NumPy form directly should expect the last bit to move; x86
and MSVC do not contract by default and are unaffected either way.

**Measured, six interleaved pairs twice over**, because the first sweep
ran into rising background load and the second into both that and a
zero-record run:

| sweep | ratios |
|---|---|
| 1 | 0.923, 0.924, 0.975, 1.152, 0.940, 0.780 |
| 2 | 0.885, 0.939, 0.915, 0.912 (+2 discarded) |

**Median 0.923 over the ten valid pairs, nine of them below 1.0** — but
**no clean separation**, and by this project's own standard that is worth
less than the number suggests. What makes it believable anyway is that an
independent measurement predicted it before it existed: the ceiling probe
separated cleanly at 0.928, and a fused pass that recovers essentially
all of a 0.928 ceiling is what 0.923 says. Two measurements of different
things agreeing is the evidence here, not either one alone.

The prediction stated before implementing was 0.94-0.97x, reasoning that
the fused pass still moves intensity and variance once and so recovers
only part of the ceiling. It recovered all of it, and slightly more —
most likely because the probe arm still made NumPy calls that take and
drop the GIL, where the fused pass releases it once.

The second sweep also cut correction's own cost, visible in the stage
instrument: 89.3 -> 50.8 ms/frame of on-core CPU, a 43% cut in the stage
itself, which is the mechanism doing what it was built to do even where
the wall time is noisy.

### C. Extend the prepare pool to the per-frame path

The grouped scheduler already hoists load and correct into a prepare
pool; the per-frame scheduler does not, so a job outside the grouping
regime — a coarse scan, an interlaced one — still queues its GIL-held
Python behind its own native call. With B done, this is mostly
bookkeeping.

**Predicted: little on the reference job** (which groups), and up to the
same 0.87x grouping bought, on jobs that do not. Sharpened from the stage
table before implementing: **0.85-0.95x on the per-frame arm**, since
that path spends 68.8 ms/frame of compute-worker wall on correction, 17%
of the occupied time of a thread that is 81% busy.

**Built, measured, and reverted.** `--group 1`, six interleaved pairs,
234 frames, depth 0:

| pair | correction in worker | correction in pool | ratio |
|---|---|---|---|
| 1 | 76.6 | 72.4 | 0.946 |
| 2 | 69.5 | 73.7 | 1.060 |
| 3 | 69.8 | 72.3 | 1.036 |
| 4 | 69.5 | 71.9 | 1.035 |
| 5 | 69.3 | 72.5 | 1.046 |
| 6 | 69.7 | 71.1 | 1.020 |

**1.035x — a regression**, in five pairs of six; the one pair that
favoured it is the first, against that arm's slowest run. Outside the
predicted band and on the wrong side of 1.

**Why, and it is this document's GIL mistake in miniature.** The
prediction assumed correction was *serialised* on the per-frame path
because it sits on the compute worker. It was not: that path runs six
compute workers, so correction was already spread across six threads.
Hoisting moved it into a reader pool that had shrunk to **one** worker —
the pool sizes itself on a blocked-fraction signal, compute was never
starved, so it had every reason to shrink — and it made each queued frame
hold its corrected arrays (105 MB) instead of its raw image (50 MB),
buying fewer frames in flight for the same budget.

The grouped scheduler's prepare pool pays because a group is *one* native
call, so its correction really was queued behind it. The per-frame path
had no such queue to remove. **The code is not kept**; what is kept is
the measurement, and the rule it illustrates: on this pipeline, "move the
Python work off the compute thread" is only worth anything where the
compute thread is a bottleneck of one.

## Ordering

A, then measure; B, then measure; C last. A is bit-for-bit and cheap. B
is where the prize is and where the risk is. C only pays for jobs the
first two do not already cover.

Do not do all three and measure once. Phase 4 is the standing argument:
a saving of known size need not convert, and finding that out one change
at a time is what makes it interpretable.

*This held, and it earned its keep.* Measuring one change at a time is
what made B's failed premise visible without taking A's win down with it,
and it is the only reason C's regression is attributable to C rather than
to a bundle.

## Where this leaves the pipeline

| | before | after A | after A and B |
|---|---|---|---|
| ms/frame, depth 0, coarse grid | 111.5 | 91.4 | ~84 |
| cores in use of 24 | 6.40 | 8.40 | 8.2 |
| GIL contention | 23% | 19% | 17% |
| correction, on-core ms/frame | 89.3 | 93.7 | 50.8 |

Absolute values across those columns come from different runs and drift
with the page cache — the ratios are the measured quantities (0.840 and
0.923, each over six interleaved pairs) and the table is the shape, not
the evidence.

**What is left, and what is not.** The document opened claiming 2.1x of
headroom bounded by a 51.5 ms/frame cold-read floor. A and B together
take about 23% of the wall, and the remaining gap to that floor is not
reachable by this route: at 8 cores of 24 with the GIL at 17%, what
limits the pipeline is that 700-odd ms of CPU per frame costs twice its
uncontended cycles. The next real lever would be *less work*, not more
overlap — the checkpoint route is now the largest Python-side term on the
compute worker (36 ms/frame of on-core time, and 232 on the fine grid),
and nothing in this document touched it.

**The one general result worth carrying forward:** on this machine,
removing memory traffic pays wherever it happens, including on threads
that are idle 87% of the time, while moving work between threads pays
only where it was genuinely serialised. Every measurement in this
document is an instance of that, including the two that came out wrong.

## Risks

**The ceiling is low and the noise is high.** 2.1x is the best case, and
depth 0 is the regime the findings document flags as dominated by page
cache state — one arm spread 176.5 to 340.6 ms/frame. Six interleaved
pairs minimum, and the foreign-load guard on every run. It is entirely
possible for a real 15% win to be invisible here.

*Both of these happened.* The noise was the binding difficulty
throughout: A separated cleanly, B did not and rests on agreeing with an
independently measured ceiling, and one sweep had to be repeated for
background load and another partly discarded for a zero-record run.

**B changes reconstructed values.** Not in which voxels are reached, but
in the last bits of the corrected intensities that feed them. That is
acceptable on the same grounds the phase 3 change was — the feature has
never been released — but it must be said, and it must not be discovered
by a user.

*It does not.* Scoping B to the element-wise arithmetic made it
bit-for-bit, verified by test on both variance branches and end to end on
a real 78-frame window. No reconstructed value moves, so nothing needs
saying to a user beyond that mapping got faster. The `phys` scope stays
on the commit anyway: the arithmetic that produces intensities and
variances moved into another language, and that is what the scope is for.

**I/O may already be the wall.** If step 1 of the measurement shows load
sitting at the cold-read floor with correction overlapped, there is
nothing here worth building, and saying so is the correct outcome.

*It was not the wall, and neither was the GIL.* Step 1 found correction
fully overlapped, as this risk anticipated, but load nowhere near the
cold-read floor and the machine short of memory bandwidth instead. The
document was right that the answer might not be what it assumed, and
wrong about which way.

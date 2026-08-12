# Reciprocal-space reconstruction: scratch-file and threading architecture

Status: single-node and cluster implementation complete (Phases 1-3 of the
implementation plan: calibration probe/file-count formula, HDF5 checkpoint
format/tree-merge accumulator, run_job integration, and Sec13 cluster
coordination via `run_cluster_map_task`/`run_cluster_finalize`/
`generate_cluster_scripts`). Cluster array size is never read from a
scheduler environment variable (most schedulers besides Slurm do not
expose one to a running task); it is baked into the generated commands as
an explicit `--total-tasks` argument instead.

Sec7 (thread allocation policy) is fully implemented, in two parts. Phase
4a: the previously fully-serial per-frame mapping loop is now a
producer/consumer pipeline (`_AdjustablePool`, `_BoundedGate`,
`_map_pending_ranges` in `orgui/reconstruction_job.py`) -- a small
prefetch-reader pool that adapts continuously via a blocked-fraction
signal (asymmetric hysteresis: grow eagerly, shrink cautiously) feeding a
depth-bounded queue that a compute-worker pool drains. Phase 4b: when
`ReconstructionJob.threads_per_image` is `None` ("automatic", the new
default -- a concrete int remains a supported pinned advanced override),
`_map_pending_ranges` seeds `kernel_threads` from how many frames the
memory budget can hold in flight (`ceil(thread_budget /
affordable_workers)`, from the same `_frame_parallelism` the pool is
sized with) and then rebalances `kernel_threads`/`image_workers` live
against the reader pool's own measured frame-delivery rate, via a
wall-clock `_kernel_threads_sweep` (`orgui/datautils/xrayutils/reconstruction.py`)
and the design doc's own joint-balancing rule (largest feasible
`kernel_threads` from a thread-budget-derived candidate set, hysteresis
against noise). A `kernel_threads` change swaps the whole compute-worker
generation (each worker's native kernel is built once, at worker-start);
an `image_workers`-only change just resizes the existing pool. This is a
schema-affecting change: `JOB_SCHEMA_VERSION` 4 -> 5, since
`threads_per_image`'s type is now `int | None` and its default changed
from `4` to `None`.

> **Revised after measurement.** The seed and the rebalance cadence
> described above both changed once they were measured against real
> beamtime data; see
> `reciprocal_space_performance_findings.md`.
>
> The seed was an I/O-optimistic `kernel_threads=1`. That is only
> sensible if memory affords enough concurrent frames to occupy the
> thread budget, and it does not: the budget caps concurrent frames at
> nine at `center` accuracy and five at `balanced`, so one native thread
> each ran a 24-thread job on five threads — measured 1.85x slower than
> using them all, and the rebalance could not correct it until its first
> interval had elapsed.
>
> The cadence was a fixed `_REBALANCE_INTERVAL_SECONDS`, illustratively
> 600s. One fixed value serves neither end of the range: a job shorter
> than the interval never re-evaluated at all, and every longer job spent
> its first ten minutes on the seed. It is now
> `_REBALANCE_INITIAL_SECONDS` (30s), doubling by `_REBALANCE_BACKOFF` up
> to `_REBALANCE_MAXIMUM_SECONDS` (the same 600s ceiling) while nothing
> changes, and resetting whenever the pair does change.
>
> The joint-balancing rule itself also compared a frames-per-second
> delivery rate against the time to map the sweep's *sample tile* rather
> than a whole frame, understating the required worker count by the ratio
> of the two (6x on a 6.2-megapixel detector). `_kernel_threads_sweep`
> now takes `frame_pixels` and scales its result accordingly.

**Sec7 now has a second scheduler beside the one described above.** Where a
job's geometry makes consecutive images overlap in reciprocal space,
mapping maps several of them in one native call
(`_map_frame_groups_streamed`), and that is a different scheduling problem:
one call covers `F` images, so there is no `image_workers`/`kernel_threads`
pair left to balance. It is replaced by a fixed number of concurrent group
calls, fed by a pool that loads *and corrects* — moving the GIL-held Python
work off the compute worker so it overlaps the previous call instead of
queueing behind it. The group size is measured per job. The pipeline
described above is what runs when a job maps one image per call, which is
what happens outside that regime. See
[`reciprocal_space_mapping_locality_and_subdivision.md`](reciprocal_space_mapping_locality_and_subdivision.md),
phase 3, which also records the two premises of Sec7 that this measured
false.

Sec14 (UI reporting) is implemented: `ReconstructionDialog`'s "Output
grids" tab live-reports estimated checkpoint file and cluster-job counts
(`ReconstructionDialog._refresh_file_count_summary`, reusing
`estimate_checkpoint_plan` from Phase 1 unchanged). This completes every
phase of the implementation plan; the redesign described by this document
is done.

**Reading this document after the fact.** It records the design as
reasoned *before* the pipeline was profiled against real beamtime data.
Several of its quantitative choices did not survive that measurement, and
the passages affected now carry an inline "Revised after measurement"
note pointing at
[`reciprocal_space_performance_findings.md`](reciprocal_space_performance_findings.md),
which holds the measurements, the settled questions, and what remains
open. The architecture — checkpoint scratch format, tree-merge
accumulator, prefetch/compute pipeline, cluster coordination — is
unchanged and still described accurately here. What changed are constants
and two rules: the `kernel_threads` seed, the rebalance cadence and its
unit error, how the memory budget is divided, the native work-block size
and detector tiling, and the record row size.

One structural addition postdates that: the unit of native work is no
longer necessarily one image, so Sec7's thread policy is one of two
schedulers rather than the only one, and Sec9's accumulator receives one
batch per *group* rather than per frame. Both are described in
[`reciprocal_space_mapping_locality_and_subdivision.md`](reciprocal_space_mapping_locality_and_subdivision.md);
neither changes the scratch format, the resume model, or the cluster
coordination this document specifies.
Relates to: issue #71, `orgui/datautils/xrayutils/reconstruction.py`,
`orgui/reconstruction_job.py`, `orgui/reconstruction_cluster.py`,
`orgui/app/ReconstructionDialog.py`

**On the numbers in this document.** It contains many specific measured
values — thread counts, RAM budgets, I/O throughput figures, per-call CPU
timings, compression ratios — all measured on one development machine
during this investigation. They are evidence for the architectural
conclusions (which mechanism to build, which formula to use, which
format/algorithm to adopt), **not defaults, thresholds, or configuration
values to hardcode**. None of the following may be fixed in the
implementation:

- **Thread budget** (`total_threads`) and **RAM budget**
  (`ram_budget_bytes`) — always the user's own resource request for that
  job/node. `24` appears often below because it's this document's own test
  environment, not a recommended value.
- **I/O throughput** / achievable frame delivery rate — discovered live
  (§7), since it depends on the real storage, network, and concurrent
  cluster load at run time, none of which any benchmark run ahead of time
  can represent.
- **Per-call compute cost** (`per_call_time(kernel_threads)`,
  mapping/reduction timings) — depends on the actual CPU, the actual
  crystal orientation, and the actual grid/depth settings of each job;
  measured by the §6 probe for every job, not read from this document.
- **Tail-flush latency** (§10) — likewise measured locally per job, not
  carried over from this investigation's numbers.

What *is* meant to carry over as a fixed decision is the qualitative,
structural conclusions — HDF5 over Parquet (§8), columnar over compound
layout (§8), XXH3 over SHA-256 (§8), tree-merge over naive incremental
merge (§9) — because those are comparisons of one approach against another
expected to hold on ordinary hardware, not absolute throughput figures to
target or replicate.

## 1. Problem statement

The current out-of-core reconstruction pipeline (introduced in `d2421e1`,
"feat: add centralized reciprocal-space reconstruction") ties scratch-file
count to **grid structure** (`partition_chunk_span` for map-phase Parquet
partitions, and one reduce-phase Parquet file per `(chunk_id, shard)` pair)
rather than to **data volume versus available memory**, which was the
original intent. On a real production job
(`39_1-rsmap`, CHESS QM2, La3Ni2O7, Pilatus 6M, 1000x1000x1500 `hkl` grid,
3651 frames) this produced 4096 reduce-phase shard files from 58 map-phase
tasks, and a live re-run at `max_depth=0` crashed with
`OSError: [Errno 24] Too many open files` partway through the reduce phase
(measured: mapping completed in 882.7 s, reduce crashed 17 s later at
shard 63/4096).

Root cause: `_ParquetRangeReader` objects are cached per worker per
partition URI in `reduce_plan_batch()` (`reconstruction.py:1294`), and are
only released when the `(grid_name, bucket)` group changes
(`readers.clear()`). With `partition_chunk_span=4096` the entire job's grid
fit in one bucket, so the group never changed and the cache never cleared:
16 reduce workers each accumulated up to 58 open partition readers,
~928 concurrently open files, past the platform's file-handle limit
(most acute on Windows, whose default per-process limit is far below
Linux's).

Separately, `threads_per_image` (a per-job setting, inherited unchanged
across accuracy presets) does not scale with `max_depth`. At `max_depth=0`
the kernel does ~1.8 CPU-seconds of work per frame, but `threads_per_image=8`
(tuned for the job's original `max_depth=3` run) capped concurrent frame
processing to `spec.threads // 8` workers, wasting most of the available
thread budget on 8-way parallelism for work too small to benefit from it.
Measured: mapping phase took 882.7 s versus a ~70-140 s estimate from
isolated per-frame kernel-cost measurements at proper parallelism.

## 2. Goals

- Scratch/checkpoint file count should be driven by
  **data volume vs. memory budget**, not by output-grid chunking.
- File count should be small: on the order of a user-configured
  checkpoint count (default 10), not thousands.
- Checkpointing should still support resuming an interrupted job.
- Fix the file-handle exhaustion bug as a byproduct of the redesign, not as
  a separate patch.
- Address the `threads_per_image` mismatch: image-level and kernel-level
  parallelism should be chosen from a measured, settings-aware cost
  estimate, not a static inherited value.
- No resource budget or environment-dependent rate — thread count, RAM
  budget, I/O throughput, per-call compute cost — may be hardcoded from
  this document's own measurements. Every one is either the user's actual
  resource request or a quantity the implementation measures live, on
  whatever machine and job are actually running (see the note at the top
  of this document).

## 3. Non-goals

- Changing the scientific algorithm (adaptive footprint splitting itself,
  voxel binning, variance propagation) is out of scope. This document is
  about the I/O/scheduling layer around the existing kernel.
- Redesigning the GUI beyond the reporting fields specified in §14.
- A general-purpose file-format migration outside reconstruction scratch
  data (Parquet remains available/optional elsewhere if needed).

## 4. User-facing parameters

Two parameters, both per job:

- **`ram_budget_bytes`** — combined memory for the whole job: all threads,
  one cluster node. Already exists conceptually as
  `runtime_memory_bytes`/`memory_override_bytes`; this redesign makes it
  actually drive file count, not just gate the native kernel's memory
  precheck.
- **`checkpoint_count`** — minimum number of resumable checkpoint files.
  **Default 10.** A floor, not a target: actual file count can exceed it
  if data volume requires more (see §5).

## 5. File-count formula

```
files_per_job = max(checkpoint_count, ceil(job_data_bytes / ram_budget_bytes))
```

`job_data_bytes` is the **reduced** (post in-frame-dedup) record volume for
the job's frame-range slice — not raw per-pixel/per-corner records, and
explicitly decoupled from output-grid chunk count.

Checkpoint frame ranges are **contiguous**, computed once and persisted
(see §11) — not recomputed on every resume, because the estimator in §6 has
sampling variance and recomputation could shift boundaries between runs
with identical settings.

### Cluster case

Each cluster job (SGE/Slurm array task, one node) is assigned a disjoint
contiguous frame-range slice up front, and applies the formula
independently to its own slice's `job_data_bytes`. Total files across the
cluster = `sum(files_per_job over all nodes)`. No cross-node coordination
is needed to compute this (see §13).

## 6. Estimating `job_data_bytes` and per-frame cost: the calibration probe

A single mechanism serves both the file-count formula (§5) and live UI
responsiveness (per-image time estimates), because a single native
`accumulate(..., profile=True)` call already returns both
`block_mapping_cpu_seconds` (time) and `reduced_block_records`
(post-dedup record count -> bytes) from the same sampled pixels.

### Why it doesn't need a real loaded image

Coordinate-evaluation cost depends only on geometry (corner rays, angle
bounds, grid, UB) and the mask, not on pixel intensity values. The probe
uses the **real mask** and **real angle bounds** for a representative
frame, but **dummy zero intensity/variance** — decoupling the probe from
image I/O and disk-speed variability entirely.

### Two-pass, budget-bounded (~100 ms) design

1. **Bootstrap pass** (a few ms): run the kernel on one tiny tile
   (e.g. 16x16 = 256 pixels) to get a preliminary per-pixel rate
   `r0 = t0 / 256`.
2. **Sized pass** (remaining budget): compute a target sample size
   `N = remaining_budget / r0` (clamped to a sane range), split across
   several small tiles **scattered across the detector** (stratified, not
   one contiguous block — per-pixel split cost varies spatially, e.g. near
   voxel boundaries or the specular rod). Run one `accumulate()` call per
   tile, sum `block_mapping_cpu_seconds` and pixel/record counts across all
   of them.

This self-scales usefully with `max_depth`: at high depth (measured
~64 us/valid-pixel at `max_depth=3` on the real job) the 100 ms budget only
affords a small, noisier sample — acceptable, because at that cost regime
the qualitative answer ("expensive") is unambiguous even with noise. At low
depth (~0.3 us/valid-pixel at `max_depth=0`) the same budget affords a much
larger, more precise sample, which is where fine distinctions (depth 0 vs
1) are more likely to matter.

### Known limitation

A small tile under-counts cross-tile-boundary duplicate voxel hits relative
to a full frame, so the record-count half of the estimate is slightly
pessimistic (more records than reality). Should be surfaced in the UI as a
conservative estimate, not treated as a correctness bug.

### Multi-grid extension

The probe naturally extends to jobs with multiple simultaneous output
grids (`spec.grids` is already a tuple): each grid has its own kernel
(`_kernel_for_grid`), and running each grid's kernel on the same sampled
tile costs little extra since no new image data or geometry sampling is
needed. Per-grid `job_data_bytes` and per-grid `files_per_job` fall out
directly (§12).

## 7. Thread allocation policy: `image_workers`, `kernel_threads`, and I/O prefetch

Replaces the static, depth-unaware `threads_per_image` setting responsible
for the mapping-phase slowdown in §1 (measured 882.7 s against an estimated
~70-140 s, traced to `threads_per_image=8` capping frame-level concurrency
to `spec.threads // 8 = 3` workers at `spec.threads=24`).

> **Revised after measurement.** The `work_block_pixels=4096` used
> throughout this section's sweeps is no longer a fixed value: the block
> is now chosen by cache scale and adaptive depth
> (`WORK_BLOCK_PRESETS`, default `medium` = 16384 pixels at `center`,
> halving per depth — so 4096 at `balanced`, which is what this section
> measured). Detector tiles are also full-width row bands of equal height
> rather than a square grid. Both change the absolute timings tabulated
> below, though not the shape of the thread-scaling argument they were
> gathered to make.

> **This section now describes one of two schedulers.** Mapping maps
> several adjacent frames in one native call where the job's geometry
> supports it, and a frame group is a different scheduling problem: one
> call covers `F` frames, so the joint `image_workers`/`kernel_threads`
> balance below is replaced by a fixed number of concurrent group calls
> fed by a prepare pool that also does the correction. The rule that
> survives is the *metric* — aggregate throughput per committed thread,
> not per-call latency — and the blocked-fraction signal, which now sizes
> the prepare pool. Everything below still describes the path taken when
> a job maps one frame per call, which is what happens outside the regime
> (coarse scans, interlaced frame order, tight memory).
>
> One claim in this section was tested directly and failed. A group call
> has thousands of bricks, which suggested it would saturate the thread
> budget on its own and make `kernel_threads = spec.threads` with one call
> in flight the natural layout. Measured, that is **1.38x worse** than
> three concurrent calls of eight threads. Available work is not what
> limits thread scaling here. See
> `reciprocal_space_mapping_locality_and_subdivision.md`, phase 3 step 5.

### The metric: aggregate throughput per committed thread, not per-call latency

A policy that picks `kernel_threads` to minimize a single call's wall time
picks wrong, because it ignores what those threads aren't doing (mapping
other frames). Measured directly at `max_depth=0` (corrected
`work_block_pixels=4096` — an initial pass with `work_block_pixels=1<<30`,
intended as "unlimited," accidentally disabled all internal blocking
(`blocks = ceil(pixels/block_size)` collapsed to 1 for every tile tested)
and gave misleadingly flat results; corrected below):

| pixels     | t=1     | t=4     | t=8     | t=16    | t=24    |
|-----------:|--------:|--------:|--------:|--------:|--------:|
| 65,536     | 6.08ms  | 2.41ms  | 2.76ms  | 4.14ms  | 3.79ms  |
| 262,144    | 18.45ms | 7.82ms  | 5.86ms  | 6.12ms  | 7.28ms  |
| 1,048,576  | 48.05ms | 17.44ms | 14.99ms | 13.48ms | 13.41ms |

Real speedup exists up to ~4-8 threads, then flattens hard. Converting the
1,048,576-pixel row to aggregate throughput for a 24-thread budget — the
scale this design will actually be implemented and tested against
(`image_workers = 24 / kernel_threads`, throughput = `image_workers /
call_time`):

| `kernel_threads` | ms/call | `image_workers` @24 threads | aggregate frames/s-equiv |
|------------------:|--------:|------------------------------:|---------------------------:|
| 1                  | 48.05   | 24                             | ~500                       |
| 4                  | 17.44   | 6                              | ~344                       |
| 8                  | 14.99   | 3                              | ~200                       |
| 24                 | 13.41   | 1                              | ~75                        |

(These ratios — 2.5-7x between `kernel_threads=1` and the larger values —
are independent of `total_threads`: `image_workers` scales linearly with
the budget for any fixed `kernel_threads`, so the throughput *ratio*
between two `kernel_threads` choices cancels the budget out. Only the
absolute frames/s-equiv numbers change with thread count, not the shape
of the comparison.)

Naively, `kernel_threads=1` looks like the winner by 2.5-7x. **That
conclusion is wrong** — it silently assumes all 24 threads can be kept
fed as 24 independent single-threaded workers, which only holds if
something can actually deliver that many concurrent frames. The next
subsection shows that this assumption is genuinely uncertain at this
thread scale, not clearly true or false. The threads that can't be used
as additional `image_workers` aren't protected by setting
`kernel_threads=1` — they just sit idle, when the measured data above
shows they could instead be speeding up the calls already in flight for
free. See "Adopted" below for the corrected, unified rule.

### At low depth, I/O is the real ceiling — the reverse holds at higher depth

The ~2000 frames/s-equivalent figure above is compute-only and, **at
`max_depth=0` specifically**, does not survive contact with real reads.
Measured against the real scan data (concurrent `fabio.open()` reads via
`ThreadPoolExecutor`, same 64-file subset reused across all trials):

| readers | frames/s | MB/s   |
|--------:|---------:|-------:|
| 1       | 26.6     | 662.5  |
| 2       | 92.4     | 2301.3 |
| 4       | 155.0    | 3859.0 |
| **8**   | **176.6**| **4395.9** |
| 16      | 170.3    | 4238.7 |
| 32      | 168.8    | 4202.3 |

Read throughput saturates at 8 concurrent readers, ~176 frames/s. At a
24-thread budget — the scale this design will actually run on, not a
large cluster node — compute at `kernel_threads=1` (all 24 threads spent
as `image_workers`) only reaches ~83 frames/s (24 / ~290 ms-per-frame,
extrapolated to full-frame scale), which doesn't even reach that I/O
figure. At this thread scale compute and I/O sit much closer together
than a large-budget illustration would suggest — plausibly on either
side, depending on the real (unmeasurable-in-advance, and itself
unreliable as measured here) I/O rate. This is weaker and more honest
than claiming I/O is clearly the ceiling at `max_depth=0`: at 24 threads
it may not be, which is exactly why the rule below discovers the regime
live rather than assuming a direction.

This uncertainty does not extend to `max_depth>0` in the same way. Per-call
compute cost grows steeply with depth (the §17-confirmed sweep: ~4.9x per
depth step at fixed tile size), while the achievable I/O rate does not
move — and at only 24 threads to begin with, there is far less compute
headroom to spare before the crossover. At `max_depth=3`, 24
single-threaded workers manage only ~0.066 frames/s in aggregate (worked
out below) — nowhere close to any plausible I/O rate, however low — so
compute, not I/O, is unambiguously the ceiling from `max_depth>=1`
onward at this thread scale. The rule adopted below is built to discover
which regime applies live, precisely because the answer is genuinely
close at `max_depth=0` and not close at all past it.

**This number must not be hardcoded into the policy — the measurement
itself is unreliable, for two independent reasons:**

1. The benchmark reused the same 64 files across every concurrency trial.
   After the `readers=1` pass, those files were very likely sitting in the
   OS page cache for every subsequent trial, so everything from
   `readers=2` onward may be measuring cached-read plus CBF-decode
   throughput, not storage I/O at all.
2. Even with clean methodology (distinct, never-touched files per trial),
   this was run against a local disk copy (`D:\sarker-4910-f`), not the
   real NFS mount (`/nfs/chess/...`) the production job reads from. NFS
   has its own client/server caching layers, network latency, and
   contention from other users' concurrent jobs on shared beamline
   storage — none of which a single-user local-machine benchmark can
   reproduce. A "corrected" version of this same benchmark would still
   only describe local SSD behavior.

### Adopted: match compute throughput to observed I/O rate, spend surplus threads on `kernel_threads`

Given the I/O ceiling cannot be reliably pre-measured — not as an
implementation gap, but because NFS caching and multi-user cluster load
make advance measurement structurally unreliable — the policy discovers it
live rather than pre-computing it, and uses it as the target that
`kernel_threads`/`image_workers` are jointly balanced against. Three
pieces:

1. **A small fixed-default prefetch pool** (e.g. 4 workers), decoupled
   entirely from `image_workers`/`kernel_threads` sizing, replacing
   today's fully serial per-task frame loop (`_map_frame_range` does not
   overlap image N+1's read with image N's processing at all — the only
   overlap that exists today is incidental, from different concurrent
   tasks happening to be at different stages). Mirrors EpiQ-Map's
   one-frame-ahead prefetch pattern (`ThreadPoolExecutor(max_workers=1)`
   in `auto_rsm.py`'s `iter_frames()`), generalized to a small pool. The
   default size is chosen because diminishing returns past a handful of
   concurrent reads is a common pattern across storage backends, not
   because 4 is known to be optimal for this job's actual NFS mount — it
   isn't claimed to be.
2. **Live I/O-rate adaptation.** Track, during actual job execution, how
   often compute workers sit blocked waiting on a prefetched frame (I/O
   can't keep up -> grow the prefetch pool, up to a cap) versus how often
   the prefetch queue is comfortably full ahead of compute (I/O keeps up
   easily). This produces a running estimate of the achievable frame
   delivery rate against the real production storage, not a synthetic
   benchmark that cannot represent it.
3. **Joint `kernel_threads`/`image_workers` balancing against that
   estimate.** The §6 probe is extended to measure `per_call_time` at a
   small candidate set of `kernel_threads` values (1, 2, 4, 8, ... up to
   where the sweep in this section shows returns flatten), not just one.
   Given the current I/O-rate estimate `r` (frames/s) from (2), choose the
   *largest* `kernel_threads` from that candidate set for which
   `ceil(r * per_call_time(kernel_threads)) * kernel_threads <=
   total_threads` still holds — i.e. the most threads concentrated per
   call that still leaves enough `image_workers` to keep up with `r`. Any
   thread capacity beyond what's needed to match `r` goes toward
   `kernel_threads`, not toward more `image_workers` sitting idle waiting
   for frames that aren't arriving any faster.

   > **Note on units.** `r` is frames per second, so the quantity it
   > multiplies must be the time to map a *whole frame*. The
   > implementation measured a single sample tile (`_KERNEL_SWEEP_TILE_PIXELS`,
   > ~1 megapixel) and used that directly, understating the requirement
   > by the ratio of a frame to that tile — 6x on a 6.2-megapixel
   > detector, enough that every candidate looked feasible and the choice
   > fell through to the "largest" tie-break rather than the measurement.
   > `_kernel_threads_sweep` now scales its result to a frame.

This is one continuous rule, not a per-depth special case, and it produces
both ends of the range from the same computation:

- **At `max_depth=0`**, with a 24-thread budget, take a purely
  illustrative I/O-rate estimate of 20 frames/s (a made-up number for
  exposition — not a measurement, unlike the unreliable-but-real ~176
  figure above, chosen here specifically because it's low enough to be
  achievable within 24 threads at every candidate `kernel_threads`,
  unlike 176). Matching 20 frames/s needs `ceil(20 * 0.290) = 6`
  `image_workers` at `kernel_threads=1` (6 threads used), `ceil(20 *
  0.104) = 3` at `kernel_threads=4` (12 threads used), `ceil(20 * 0.089) =
  2` at `kernel_threads=8` (16 threads used) — and `kernel_threads=16`
  would need 2 workers x 16 threads = 32, exceeding the 24-thread budget,
  so it's infeasible. The rule picks `kernel_threads=8` here: the largest
  candidate that still fits, using 16 of 24 threads with headroom to
  spare. The exact numbers depend entirely on the live-measured rate, not
  this illustration — a higher live estimate would push the choice back
  toward lower `kernel_threads`/more `image_workers`, a lower one further
  toward concentrating threads per call.
- **At `max_depth=3`**, per-call time balloons to ~366 CPU-seconds at
  `kernel_threads=1`. 24 such single-threaded workers could only manage
  ~0.066 frames/s in aggregate — far below any plausible I/O rate, however
  low — so no candidate `image_workers` count at `kernel_threads=1` can
  satisfy the balancing condition in (3). The rule is then forced toward
  the largest `kernel_threads` the candidate set offers, i.e.
  `threads_per_image` "wins decisively," exactly as expected, derived from
  the mechanism rather than asserted as a separate case.

### Prefetch queue depth: bounded by compute's consumption rate, not by I/O's delivery rate

The prefetch pool from (1) has two distinct parameters that must not be
conflated: how many *reader threads* it runs, and how many *decoded
images* it is allowed to hold buffered ahead of compute. The
`blocked_fraction` signal in (2) correctly governs the first. It says
nothing about the second, and leaving the second unbounded is a real
failure mode in the compute-bound regime: at `max_depth=3`, I/O can
deliver a frame roughly every ~43-60 ms while a single-threaded compute
call takes ~3 s (measured, 65,536 px tile) — over a 3 s window I/O could
hand over ~50 decoded images. Prefetching that far ahead would be wrong:
each decoded Pilatus 6M frame is tens of MB, and buffering dozens of them
achieves nothing (compute can only consume one at a time per active
worker) while consuming real memory that competes with §10's checkpoint
accumulators for the same `ram_budget_bytes`.

The queue only needs to stay far enough ahead that a compute worker never
finds it empty when it asks — which is a function of how many workers are
consuming from it, not of how fast it could theoretically be filled:

```
N ~= image_workers + small_constant   (e.g. +2 to +4)
```

Implemented as an ordinary bounded queue: reader threads block on `put()`
once it holds `N` images, resuming only as compute workers dequeue. This
is passive backpressure, not another instance of the live-adaptation loop
in (2) — it needs no monitoring of its own, and it automatically produces
the right behavior in both regimes without knowing which one it's in:

- **Compute-bound** (`max_depth=3`, illustratively `image_workers=1` at a
  24-thread budget once `kernel_threads` claims most or all of it per the
  rule above — the confirmed sweep shows depth-3's plateau is still
  climbing at 16 threads, so at a 24-thread total budget `kernel_threads`
  could plausibly claim the entire allocation): `N` is small (~3-5), the
  queue fills almost immediately, and the (nominally larger, e.g.
  4-worker) prefetch pool sits mostly blocked on `put()`, consuming no CPU
  regardless of how much faster the storage could go. No 50-image
  buildup.
- **I/O-bound** (`max_depth=0`, `image_workers` potentially up to the full
  24-thread budget): the queue rarely reaches `N` at all, since compute
  drains it as fast or faster than the (I/O-limited) pool can fill it —
  the cap is simply inactive, and the `blocked_fraction` mechanism in (2)
  is what's actually governing behavior here, as designed.

This adds a third constraint to the balancing rule in (3), alongside the
thread-budget and I/O-rate-matching conditions: `image_workers *
bytes_per_decoded_image <= image_buffer_budget_bytes`, a slice of
`ram_budget_bytes` reserved separately from §10's checkpoint-accumulator
allocation. In practice this rarely binds at a 24-thread scale
(`image_workers` is bounded by `total_threads=24` regardless, so total
buffer need tops out around `24 * ~50 MB ~= 1.2 GB` even in the most
I/O-bound case, against typical multi-GB-to-tens-of-GB budgets), but it is
an explicit, checked bound rather than an assumption, since it is exactly
what prevents the compute-bound case from buffering images nobody is
ready to process.

## 8. Checkpoint file format: HDF5, not Parquet

Benchmarked against the real record schema (`chunk_id`, `local_voxel_id`,
`weighted_intensity`, `weighted_variance`, `weight`, `contributors`,
48 bytes/row, 800k rows / 38.4 MB raw, matching one frame's depth-0 reduced
output):

> **Now 40 bytes/row.** `chunk_id`, `local_voxel_id` and `contributors`
> were narrowed from 64-bit to 32-bit (`_CHECKPOINT_BYTES_PER_ROW`), which
> changed both the checkpoint format and the final output's
> `contributors` dtype. The relative comparison below is unaffected; only
> the absolute row size is. The same correction applies to §12's
> "fixed 48-bytes-per-row schema".

| backend                | write MB/s | read MB/s | file size |
|-------------------------|-----------:|----------:|----------:|
| parquet:zstd (current)  |        334 |       775 |   22.2 MB |
| parquet:lz4             |        401 |      1598 |   24.1 MB |
| parquet:snappy          |        428 |      1753 |   24.4 MB |
| hdf5:gzip4              |         85 |       398 |   18.9 MB |
| hdf5:lzf                |        475 |       494 |   21.0 MB |
| **hdf5:bitshuffle-lz4** |    **905** |   **840** | **19.7 MB** |
| hdf5:none               |       2251 |      1726 |   44.1 MB |

`bitshuffle-lz4` — already used by `_finalize_reconstruction` for the final
output, via the already-present `hdf5plugin` dependency — writes ~2.3x
faster than the current `parquet:zstd`, reads comparably, and produces
smaller files. Switching removes the `pyarrow` dependency entirely
(currently optional-but-required: `_require_pyarrow()` raises if it's
missing) since HDF5 is already a hard dependency everywhere else in orGUI
(config, job assets, final output). h5py's resizable datasets
(`maxshape=(None,)` + `.resize()`) also fit the append-as-you-go checkpoint
model better than Parquet's row-group/file-per-flush model.

### Layout: separate 1D datasets per column, chunk=65536 rows

Benchmarked columnar (six separate datasets, one per record field) against
a single compound/table dataset, at four chunk sizes (16384-262144 rows),
same schema and row count as above:

| layout                     | write MB/s | read MB/s | size     |
|-----------------------------|-----------:|----------:|---------:|
| columnar, chunk=16384       |      851.9 |     865.1 | 19.35 MB |
| compound, chunk=16384       |      886.0 |     792.7 | 19.44 MB |
| **columnar, chunk=65536**   |  **973.5** | **935.5** | 19.36 MB |
| compound, chunk=65536       |      908.6 |     904.0 | 19.45 MB |
| columnar, chunk=131072      |      915.1 |     778.7 | 19.38 MB |
| compound, chunk=131072      |      867.4 |     876.5 | 19.47 MB |
| columnar, chunk=262144      |      887.9 |     762.2 | 19.42 MB |
| compound, chunk=262144      |      792.4 |     784.1 | 19.50 MB |

File size is essentially identical across layouts (~19.4 MB either way) —
mixed-type compound records were expected to compress worse than pure
same-type columns under bitshuffle; the measured difference is under 1%,
not a real factor. Performance differences between layouts are a few
percent either direction depending on chunk size; columnar at
`chunk=65536` is the best-measured configuration.

Adopted: **columnar, chunk=65536 rows** — the measured winner, but decided
on a second ground independent of the numbers: the in-memory
representation used everywhere else in this design (the native kernel's
output, `_merge_sorted_batches`, the §9 tree-merge `levels` structure) is
already a dict of separate arrays. A compound dataset would require
assembling a structured array before every checkpoint write and
disassembling it after every resume/finalize read — real conversion cost
on a layout that doesn't even win on raw throughput. Columnar writes each
already-existing array straight through. `chunk=65536` also gives a
sensible chunk count per checkpoint (~11-15 chunks for the ~740k-985k row
checkpoints in the §15 worked example) — fine enough to avoid forcing
full-file reads for partial access during resume verification, coarse
enough not to fragment.

### Checksums: SHA-256 -> XXH3-128, computed inline

Scratch-file verification currently uses `hashlib.sha256()`
(`_uri_checksum_and_size`), not the fast native XXH3 binding already used
for scientific-input fingerprinting. Benchmarked:

| buffer  | SHA-256          | XXH3-128 (native)   | speedup |
|---------|------------------:|---------------------:|--------:|
| 256 MB  | 0.108 s (2.4 GB/s) | 0.0066 s (39.0 GB/s) |   16.5x |
| 1024 MB | 0.436 s (2.3 GB/s) | 0.0272 s (37.7 GB/s) |   16.0x |
| 4096 MB | 1.812 s (2.3 GB/s) | 0.1202 s (34.1 GB/s) |   15.1x |

At the real job's scale (~66.5 GB compressed scratch volume total), SHA-256
costs ~29 s CPU-bound, XXH3 ~1.9 s — both trivial next to reconstruction
wall time, so the real win isn't the hash speed, it's that going from 4096
files to ~10-30 also removes ~400x the per-file overhead (opens,
Python-level calls, thread submission).

Checksums should be computed **inline during write**, not by re-reading
the file afterward (`_uri_checksum_and_size`'s current pattern) — this
removes an entire extra disk-I/O pass, independent of which hash algorithm
is used, and matters more on network scratch storage (this job's actual
NFS mount) than on local SSD.

## 9. Checkpoint write path: per-checkpoint tree-merge accumulator

### Rejected: naive incremental merge

A single lock-protected accumulator, with each finishing frame merged in
one at a time via `_merge_sorted_batches`, is **quadratic** in
frames-per-checkpoint: merge cost scales with current accumulator size, so
total cost over N frames is O(R*N^2/2). Measured directly (batch size
739,791 rows, matching real depth-0 output):

| frames | naive merge | tree merge (`_reduce_batches`) | speedup |
|-------:|------------:|--------------------------------:|--------:|
|      8 |     0.28 s  |                          0.20 s |    1.4x |
|     16 |     1.08 s  |                          0.52 s |    2.1x |
|     32 |     3.90 s  |                          1.29 s |    3.0x |
|     64 |    17.00 s  |                          3.20 s |    5.3x |
|    128 |    67.21 s  |                          8.81 s |    7.6x |

Extrapolated to a real checkpoint (~365 frames at `files_per_job=10` on a
3651-frame job): naive ~546 s/checkpoint (~1.9x the checkpoint's own
mapping cost at depth 0 -- it would dominate the job), tree ~31 s/checkpoint
(~10% addition on top of mapping) — a per-frame cost (~85 ms amortized)
comparable in magnitude to a single frame's image load and correction time,
which is why §7 sizes `image_workers` from total per-frame pipeline cost
rather than kernel cost alone.

### Adopted: tree/binary-counter accumulator, per (checkpoint, grid)

State per `(checkpoint, grid)` pair:

- one lock
- `levels: list[batch | None]` — the same binary-counter structure
  `_reduce_batches` already builds internally, but persisted across time
  and fed incrementally instead of built once from a pre-collected list
- a running byte-total, for the budget check below

On a worker finishing a frame's (already kernel-side-reduced) batch:

1. Acquire the checkpoint's lock.
2. Run `_reduce_batches`'s existing insert loop: descend levels, merging
   and cascading exactly as it does today, placing the result at the
   first empty level (or appending a new one).
3. Update the byte-total; if it crosses this checkpoint's budget share
   (§10), trigger a flush before releasing.
4. Release the lock.

Insert cost is *amortized* O(log N) per call (a cascade to depth k happens
once per 2^k insertions), so total lock-held time across N insertions is
O(N log N), matching the measured tree numbers. Mapping — the expensive,
fully parallel part — happens entirely outside the lock, so contention
stays low relative to how often workers need it.

On checkpoint completion (frame range exhausted) or a mid-range budget
split: fold remaining levels together
(`_reduce_batches`'s existing finalization step), write to the
checkpoint's HDF5 file via write-`.tmp`-then-`.replace()` (§11), reset
`levels`.

### Validation strategy: equivalence against the current reduce, not an external ground truth

Neither the current per-flush reduce logic nor this replacement has a
real validation dataset — the current implementation isn't independently
verified either, though it isn't known to be wrong. That rules out
"regression testing against known-correct output" as an option for either
side. What's achievable instead: generate a small synthetic dataset from
deterministic random noise (fixed RNG seed, a handful of frames, a small
grid — no need for real detector/crystal data, per §6's finding that
mapping cost and correctness depend only on geometry and mask, not pixel
values), run it through both the current `_reduce_batches`/per-flush path
and the new per-checkpoint tree-merge accumulator, and assert the two
produce identical (or, if floating-point summation order differs,
equivalent within floating-point tolerance) `(chunk_id, local_voxel_id) ->
(weighted_intensity, weighted_variance, weight, contributors)` output.
This establishes that the new path preserves the old path's behavior,
which is the honest and available claim — not that either is correct
against some external reference, since no such reference exists to test
against.

## 10. Concurrency and memory division

Two independent axes, deliberately not conflated:

- **`image_workers`**: how many frames are processed concurrently. Sized
  per §7 (calibration probe + I/O-prefetch adaptation) — a
  scheduling/throughput decision.
- **`max_concurrent_active_checkpoints`**: how many checkpoints have
  non-empty in-memory accumulator state at once. **Small, fixed, default
  2** (current + next), not user-facing. Checkpoints are processed in
  frame-range order; workers are allowed to start claiming the next
  checkpoint's frames slightly before the current one fully drains, to
  avoid idle-worker gaps at boundaries, without letting the active set
  grow unbounded.

```
active_budget_per_accumulator =
    ram_budget_bytes / (max_concurrent_active_checkpoints * number_of_grids)
```

Dividing by 2 (not `image_workers`, and not `files_per_job`) is not just
simpler, it's more memory-efficient: checkpoints that haven't started yet
hold zero memory, so there's no reason to reserve a slice of the budget for
all of them upfront. A smaller divisor gives each active accumulator more
headroom, meaning fewer safety-valve splits in practice.

> **Revised after measurement.** This formula hands the accumulators the
> *whole* budget: two active checkpoints, each allowed half of it. The
> frame pipeline was sized against the whole budget independently, so
> nothing constrained their sum, and measured peak resident memory
> reached 12.4-14.5 GiB against a 10 GiB budget.
>
> `split_memory_budget()` now divides the budget once —
> `CHECKPOINT_MEMORY_SHARE` (0.25) to the accumulators, the rest to the
> pipeline — so the numerator above is a quarter of `ram_budget_bytes`,
> not all of it. The reasoning for the divisor of 2 is unchanged; only
> what is being divided.
>
> The "more headroom means fewer splits" argument holds but is worth
> less than it sounds: quartering the accumulator share grew a checkpoint
> set by 6%, because voxel overlap is concentrated among *neighbouring*
> frames — the crystal rotates, so frames far apart touch largely
> disjoint reciprocal space, and a smaller window already captures nearly
> all the deduplication. Mapping was also slightly *faster* with the
> smaller share. See `reciprocal_space_performance_findings.md`.

Worked example (real job, one arbitrary choice of `ram_budget_bytes=24 GB`
picked to illustrate the formula, not a recommended value — see the note
at the top of this document — with `files_per_job=10` from the checkpoint
floor dominating): each active checkpoint gets `24/2 = 12 GB`. The
calibrated per-checkpoint estimate was ~12.97 GB raw — right at the edge,
meaning this specific combination would occasionally trip the safety valve
near the end of a checkpoint, producing 11 files for that run instead of
10. That's the safety valve working as intended at the margin, not a bug:
the calibration is close but not exact, and the system self-corrects
instead of overrunning the real memory budget.

### When `max_concurrent_active_checkpoints=2` isn't enough

The "current + next" window only avoids idle-worker gaps if checkpoint N
finishes committing (its final tree-merge fold and write, §9) *before*
checkpoint N+1's frames are fully claimed and drained. If N+1 drains
faster than N's tail-flush completes, workers run out of claimable work
and stall — the window can't advance to N+2 until N is actually done, not
just fully assigned. This is a real risk exactly in the regime this
document keeps returning to: at `max_depth=0`, §7's rule pushes
`image_workers` high and per-frame cost is tens of milliseconds, so a
365-frame checkpoint can drain in a couple hundred milliseconds — very
plausibly faster than N's tail-flush, which is small relative to the
whole checkpoint's merge cost (most of the tree-merge already happened
incrementally during accumulation) but not necessarily small in absolute
terms.

Unlike the I/O rate in §7, this isn't something that needs live
adaptation. I/O rate depends on external, unpredictable factors (network
conditions, other users' concurrent cluster load) that can only be
discovered at runtime. This doesn't: `image_workers` (§7),
frames-per-checkpoint (from `files_per_job` and total frame count), and
per-frame cost (§6 probe) are all already known before the job runs, so
whether checkpoint drain time (`frame_batch_size * per_frame_time /
image_workers`) will outrun the estimated tail-flush latency is a
comparison that can be made once, at prepare time, from numbers already
in hand — no second runtime control loop needed, *provided* tail-flush
latency itself has an estimate. It didn't have one; here's one that
doesn't hardcode this session's specific measured numbers.

**Tail-flush latency formula.** The tail fold combines whatever tree
levels remain when the checkpoint's last frame lands — in the worst case
(a single large final fold), that touches close to the checkpoint's full
row count. A conservative (safe-to-overestimate) bound:

```
tail_flush_latency ~= checkpoint_final_rows / merge_rows_per_second
```

where:

- `checkpoint_final_rows = (job_data_bytes / files_per_job) / bytes_per_row`
  — already derivable from quantities §5/§6 establish (the per-checkpoint
  share of `job_data_bytes`) and §8's fixed 48-bytes-per-row schema. No
  new measurement needed for this half.
- `merge_rows_per_second` — the throughput of `_merge_sorted_batches`
  itself, which is a property of the machine running the job, not of this
  session's benchmark runs. It should be measured with a small local
  calibration call (merge two ~100k-row synthetic batches, time it, divide)
  at prepare time — the same "measure it on this machine" pattern already
  used by the §6 probe and the §7 `per_call_time(kernel_threads)` sweep,
  applied once more to a third operation. This is a cheap, sub-second
  calibration (the §9 benchmark's own numbers, e.g. ~739,791-row batches
  merging in a fraction of a second, show `_merge_sorted_batches`
  throughput is high relative to typical checkpoint row counts — but that
  specific rate is this test machine's, illustrative of the *order of
  magnitude* only, not a value to hardcode into the formula).

This keeps both halves of the comparison — checkpoint drain time and
tail-flush latency — built from quantities either already known or cheap
to measure locally, portable to whatever machine actually runs the job.

**Resolution**: a prepare-time check, not a runtime one, and not a
primary UI setting:

- If the computed checkpoint drain time comes out below the estimated
  tail-flush latency, auto-bump `max_concurrent_active_checkpoints` (e.g.
  2 -> 3 or 4, capped low) as part of preparing the job — automatic, not
  something the user needs to notice or request.
- The consequence (smaller `active_budget_per_accumulator`, more frequent
  safety-valve splits, likely a few more files than the calibrated
  `files_per_job` estimate) should be visible in §14's UI reporting, not
  silent.
- A manual advanced override stays available for cases the heuristic gets
  wrong, following the same precedent as `advanced_depth` overriding
  `accuracy` elsewhere in the job config — but it is not a primary
  setting, since the default path is expected to get this right without
  the user needing to know the parameter exists.

## 11. Resume and staleness

### Invalidation

`checkpoint_count` becomes a new field on `_ReconstructionSpec`, alongside
`ram_budget_bytes` (already present as `memory_budget_bytes`). Both then
participate in the existing `spec.digest` mechanism unchanged — the same
mechanism that already invalidates manifests on any scientific-setting
change today. No new invalidation philosophy, just two new fields feeding
an existing one.

Checkpoint frame-range boundaries are computed once at prepare time (from
the calibration estimate) and **persisted**, not recomputed on resume —
recomputing from a fresh (sampling-variant) calibration run could shift
boundaries between two runs with identical user settings, spuriously
invalidating good checkpoints.

**No backward compatibility with existing prepared jobs.** This code
hasn't shipped to users yet, so job files prepared under the current
schema (`JOB_SCHEMA_VERSION=3`, Parquet-based, `partition_chunk_span`,
static `threads_per_image`) are not a migration concern. A schema version
bump simply makes them unreadable by the new code — re-preparing from
scratch is an acceptable answer, not something this design needs to
smooth over.

### Atomicity

Checkpoint HDF5 files use the same write-`.tmp`-then-`.replace()` pattern
already used twice in the codebase (`_write_manifest`,
`_finalize_reconstruction`) — a crash mid-write leaves an orphaned `.tmp`,
never a corrupt file at the real name.

### Simplification: no separate manifest registry

A checkpoint file's presence at its final name already implies
completeness (atomic rename guarantees this). Each checkpoint stores its
own spec digest as an HDF5 attribute (same pattern as
`process.attrs["spec_sha256"]` on the final output today). Resume becomes:
for each persisted checkpoint boundary, does the file exist at its final
name with a matching digest attribute? Given file counts are now small
(10-30s, not thousands) and checksums are now XXH3-fast, resume can also
cheaply re-verify each existing checkpoint's data checksum every time as a
belt-and-suspenders check against silent corruption — something too slow
to do routinely under the old design. This removes an entire layer of
JSON-manifest-list bookkeeping (`job.map_manifests`) that only existed
because the old design had thousands of files to track.

## 12. Multi-grid jobs

`files_per_job` is computed **per grid**, not shared — grids can have very
different data volumes (a coarse overview grid vs. a fine ROI grid), and
forcing shared boundaries would either over-fragment the smaller grid or
under-protect the larger one's memory footprint. Total file count for the
job is the sum across grids (e.g. 3 grids at the checkpoint floor = 30
files — still comfortably "tens, not thousands").

Different grids' checkpoint boundaries may differ freely. A worker
processing frame F routes that frame's per-grid results to whichever of
*that grid's own* checkpoint boundaries currently contains F (a per-grid
lookup) — this doesn't complicate frame-level scheduling. The §9
tree-merge accumulator generalizes from "per checkpoint" to "per
(checkpoint, grid)" with no new mechanism. The §10 memory-division formula
gets a `number_of_grids` factor since all requested grids' accumulators
for the active checkpoint(s) are alive in memory simultaneously (one
corrected image feeds a kernel call per grid before moving to the next
frame).

## 13. Cluster (multi-node) coordination

No new synchronization primitive needed. Frame ranges are disjoint across
nodes (assigned up front), and each frame's contribution to the final grid
is independent of every other frame's — so nodes need zero communication
during mapping. Each node runs the full single-node design (§5-11)
independently against its own slice and its own `ram_budget_bytes`.

`reconstruction_cluster.py`'s `generate_cluster_scripts` already submits a
finalize job with a scheduler-level dependency on the whole map array
completing (`qsub -hold_jid` / `sbatch --dependency=afterok:...`) — this is
the only synchronization point, and it already exists. The redesign
changes what finalize processes (a small number of HDF5 checkpoints
instead of thousands of Parquet shards, merged via the same §9 tree
structure) but not the barrier mechanism itself.

Remaining concerns are existing patterns, not new risks:

- **Filename collisions on shared scratch**: avoided by construction —
  each node's checkpoint files are named from its own frame-range/job
  identity, so disjoint slices can't collide.
- **Partial writes visible to other processes**: covered by the same
  atomic rename from §11, which the current code already trusts for
  manifests on the same (NFS) shared storage.
- **A node crashing mid-slice**: handled by the node's own resume logic
  (§11) plus the scheduler simply not firing the finalize dependency until
  every array task reports success; a failed node just gets resubmitted.

## 14. UI reporting

Implemented (`ReconstructionDialog._refresh_file_count_summary`,
`orgui/app/ReconstructionDialog.py`). At prepare time, before the job is
committed, the calibration probe (§6, via `estimate_checkpoint_plan`)
runs and the UI reports, live-updating as the user adjusts the grid
table, `ram_budget_bytes` / `checkpoint_count`, accuracy, angle fallback,
or cluster node count:

- `number_of_files` (total, summed across grids -- and, when the Cluster
  tab's node count is more than 1, a second reported total for that
  many nodes, since the same prepared job can still be run either
  locally or as a cluster job)
- `number_of_jobs` (cluster array size, shown alongside the cluster
  total, not unconditionally -- a job's cluster settings default to a
  node count regardless of whether the user ends up running it locally)
- `files_per_job` (per node; per-grid breakdown shown only when more
  than one grid is requested)

This goes in `ReconstructionDialog`'s **"Output grids"** tab
(`_grid_tab()`), alongside the existing per-axis grid bounds — the natural
place for it today, since file/job counts are a direct consequence of the
grid the user is defining right there, not because it's the ideal final
home. `ReconstructionDialog`'s tab layout more broadly (Experiment /
Output grids / Performance / Cluster / Job and output / Preview and
status) is expected to need a rework once this whole redesign lands, given
how much of what "Performance" and "Cluster" mean is changing — that
rework is out of scope here and left for later.

## 15. Worked example (real job)

The `ram_budget_bytes` values in the table below (4/8/13/24/128 GB) are
arbitrary points chosen to show the formula's shape across a range, not
recommended budgets — see the note at the top of this document.
`39_1-rsmap`, Pilatus 6M, La3Ni2O7, 3651 frames, `max_depth=0`. Measured
`final_records = 739,791` for a representative frame -> 35.5 MB/frame
(post-dedup) -> ~129.7 GB total job data (pre-final-merge; the conservative
planning number, since cross-frame overlap reduction can't be assumed
ahead of time). Notably this volume is roughly depth-independent
(`final_records` at `max_depth=3` was 985,539 for the same frame, same
order of magnitude) — confirming the formula correctly decouples file
count from both grid-chunk count and accuracy preset, the two things it
was conflated with before.

| `ram_budget_bytes` | `files_for_memory` | `files_per_job` (floor=10) |
|--------------------:|--------------------:|-----------------------------:|
|                4 GB |                  33 |                           33 |
|                8 GB |                  17 |                           17 |
|               13 GB |                  10 |                           10 |
|               24 GB |                   6 |                       **10** |
|              128 GB |                   2 |                       **10** |

Split 4 ways across a cluster (~913 frames/node): ~32.4 GB/node; at 8
GB/node budget, `ceil(32.4/8)=5` -> checkpoint floor wins -> 10 files/node
x 4 nodes = **40 files total**, versus today's 4096.

## 16. Migration summary

| | current | proposed |
|---|---|---|
| map-phase file boundary | grid-chunk buckets (`partition_chunk_span`) | contiguous frame-range checkpoints |
| reduce-phase files | one per `(chunk_id, shard)`, up to thousands | none — folded into the checkpoint write itself |
| format | Parquet (optional dependency) | HDF5 (`bitshuffle-lz4`, already a hard dependency) |
| checksum | SHA-256, re-read after write | XXH3-128, computed inline |
| cross-frame merge | ad hoc, per flush cycle | tree/binary-counter, per (checkpoint, grid) |
| resume tracking | JSON manifest list per partition/shard | checkpoint file presence + embedded digest attribute |
| `threads_per_image` | static, inherited across accuracy presets | `kernel_threads`/`image_workers` jointly balanced live against an adaptively-measured I/O rate, plus a small adaptive prefetch pool (§7) |
| images per native call | one | measured per job; several adjacent images share one call where the scan's geometry makes them overlap in reciprocal space |

## 17. Open items

- ~~Run the §7 `per_call_time(kernel_threads)` sweep at `max_depth=2/3`~~ —
  done. Confirmed: at a 256x256 (65,536 px) tile, depth=2 scales cleanly to
  16 threads (0.599s->0.105s) and depth=3 likewise (2.917s->0.459s), both
  flattening between 16 and 24 threads; the same tile size at `max_depth=0`
  showed no real scaling at all (noise-level, ~6ms flat). Confirms the
  plateau genuinely shifts higher with depth rather than being asserted by
  extrapolation — the balancing rule's candidate `kernel_threads` set needs
  to be depth-specific (e.g. roughly {1,2,4,8} at depth 0 vs {1,4,8,16,24+}
  at depth 3), which §7's design already accommodates by treating
  `per_call_time(kernel_threads)` as a per-depth input. Not yet measured at
  full-frame (6.2M px) scale, where the depth-2/3 plateau likely sits
  higher still than 16.
- ~~The exact form of the §7 live I/O-rate and `kernel_threads`/
  `image_workers` adaptation~~ — resolved: `blocked_fraction` signal,
  time-based re-evaluation cadence, asymmetric hysteresis bands (grow the
  prefetch pool eagerly, shrink it cautiously), and a separate
  compute-bounded prefetch *queue depth* (`N ~= image_workers +
  small_constant`) distinct from prefetch *pool size*. The specific
  threshold values (e.g. 20%/2% blocked-fraction bands, exact cadence in
  seconds) are illustrative, not tuned against real production load — that
  residual calibration is implementation work, not an open design
  question.
- ~~Whether `max_concurrent_active_checkpoints` should ever be exposed as
  an advanced override~~ — resolved: a prepare-time check (checkpoint
  drain time vs. tail-flush latency) auto-adjusts it when needed, surfaced
  in §14's UI reporting; a manual advanced override remains available
  (matching the `advanced_depth`/`accuracy` precedent) for cases the
  heuristic gets wrong, but it is not a primary setting. Tail-flush
  latency now has an explicit formula (§10) built from §5/§6/§8 quantities
  plus one additional cheap local calibration (`_merge_sorted_batches`
  throughput on this machine) — not hardcoded to this session's measured
  numbers.
- **Frame grouping changes what §7 balances, and §7 has not been rewritten
  around it.** Above one image per call the joint
  `image_workers`/`kernel_threads` problem is replaced by concurrent group
  calls plus a prepare pool; §7 carries a note saying so but still reads
  as though the joint balance is the only policy. Rewriting it properly
  needs the phase 4 measurements first, since the balance above depth 0 is
  unmeasured. Meanwhile the two paths coexist and the per-frame one is
  unchanged.
- **The group-size and concurrency constants rest on one job.** The
  one-voxel travel limit comes from a striding sweep on a single scan; the
  floor of three concurrent calls is fitted to one throughput curve. Both
  reproduce that job's optimum and neither has been tested on a second
  dataset or a second machine.
- ~~Exact HDF5 layout for checkpoint files~~ — resolved: columnar (separate
  1D datasets per record field), `chunk=65536` rows. Measured against a
  compound/table dataset across four chunk sizes; columnar at 65536 was
  both the best-performing configuration and the one requiring no
  array-assembly/disassembly conversion against the in-memory
  representation already used everywhere else in the design (§8).

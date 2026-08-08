Reciprocal-Space Reconstruction
===============================

The reciprocal-space reconstruction workflow maps corrected detector images
into one or more regular three-dimensional HKL or momentum-transfer grids. It
is designed for data sets that are much larger than memory. Pixel geometry,
angle transforms, adaptive footprint splitting, and local accumulation run in
a C++17 extension. Python coordinates scan loading, corrections, resumable
checkpointed HDF5 scratch records, and final NeXus/HDF5 output.

The experiment configuration is not entered a second time. Reconstruction uses
the active orGUI scan, detector calibration, crystal, UB matrix, mask,
background, correction selections, exclusions, CPU limit, memory limit, and
HDF5 filter registry. Preparing a job freezes those settings into a
checksummed, resumable snapshot.

Installation
------------

The native reconstruction extension is built with orGUI. ``hdf5plugin`` makes
the optional HDF5 filters registered by orGUI -- including the
``bitshuffle-lz4`` filter used by the out-of-core checkpoint scratch format --
available:

.. code-block:: bash

   pip install "orGUI[reconstruction]"

The reconstruction feature does not use Numba, OpenMP, or TBB. The C++ kernel
uses a bounded ``std::thread`` worker pool.

Opening the Workflow
--------------------

Load the scan and experiment configuration normally, then open
``Reciprocal space -> Reconstruct reciprocal space``. The dialog contains:

* ``Experiment``: active scientific state, corrections, exclusions, and
  exposure-bound policy;
* ``Output grids``: coordinate systems, bounds, steps, interval counts, HDF5
  settings, and geometry-matched step estimation;
* ``Performance``: footprint accuracy, CPU and memory allocation, and advanced
  task layout;
* ``Job and output``: job descriptor, scratch directory, and final HDF5 path;
* ``Preview and status``: estimates, prepared JSON, status, results, and
  user-facing errors.

Preparing Versus Running
~~~~~~~~~~~~~~~~~~~~~~~~

``Preview`` reads the current live orGUI state and estimates grid sizes and
execution layout without freezing anything.

``Prepare Job`` writes the JSON descriptor and immutable HDF5 asset bundle.
Large masks, backgrounds, background variances, and repair geometry are stored
in the asset bundle rather than duplicated in JSON. The JSON and asset bundle
are checksummed.

``Run Locally`` prepares the current settings and executes them. ``Open Job``
opens an existing descriptor, while ``Resume`` verifies its sources and
checksums and continues incomplete work. A completed result is registered in
the active orGUI database as an external result; the large arrays remain in the
standalone file.

Experiment Parameters
---------------------

Active Experiment
~~~~~~~~~~~~~~~~~

The read-only overview shows the active scan name, frame count, detector shape,
X-ray energy in keV, UB matrix, correction state, and global CPU and memory
limits. Active masks and backgrounds are shown as resolved inputs. Asset paths
such as ``/mask`` refer to datasets that will be created in the prepared asset
bundle.

``Refresh from active orGUI state``
   Re-reads the central state. If no grid exists, it also derives the initial
   HKL grid and chooses writable default paths under the current working
   directory.

Corrections and Exclusions
~~~~~~~~~~~~~~~~~~~~~~~~~~

``Mask and repair``
   Opens the shared orGUI mask and pixel-repair settings. Reconstruction uses
   the same active detector mask as image integration. A matching mask is
   required when mask correction or repair is enabled.

``Background``
   Opens the shared background-image control. The active background and
   optional background variance are frozen into the job asset bundle.

``Excluded frames``
   Opens the shared frame-exclusion editor. Excluded frames are omitted from
   mapping, task layout, geometry-step estimation, and coverage accounting.

``Missing exact angle bounds``
   Selects the policy used when a scan backend has no exact encoder start/end
   positions. ``Stationary exposure`` uses identical start and end positions.
   ``Midpoint inference`` explicitly estimates boundaries halfway between
   adjacent nominal positions. Exact backend bounds always take precedence.

``User note``
   Optional free text stored in the job descriptor and final provenance.

``Normalize by exposure time``
   Divides each frame by its exposure time when the active scan backend
   provides one. Enabled by default. This setting is specific to
   reconstruction and does not affect ROI/CTR image integration.

``Monitor corrections``
   Comma-separated scan counter names applied as divisive monitor
   normalizations, each with uncertainty propagated when the backend
   exposes a matching ``<name>_variance`` counter. This setting is specific
   to reconstruction and does not affect ROI/CTR image integration.

Output Grid Parameters
----------------------

Each row describes one independent output coordinate system.

``Name``
   Optional HDF5 group name. Empty names are generated from the frame name and
   made unique.

``Frame``
   Coordinate system. ``hkl`` is in reciprocal lattice units. Q frames are in
   :math:`\mathrm{\mathring{A}}^{-1}`:

   ``lab``
      Laboratory axes.

   ``alpha``
      Frame after the alpha rotation.

   ``omega``
      Frame after alpha and omega rotations.

   ``chi``
      Frame after alpha, omega, and chi rotations.

   ``phi``
      Frame after the full sample rotation.

   ``crystal``
      Crystal-oriented Cartesian Q coordinates using the orientation matrix.

``Min 1``, ``Min 2``, ``Min 3``
   Lower voxel edges. Units are r.l.u. for HKL and
   :math:`\mathrm{\mathring{A}}^{-1}` for Q.

``Max 1``, ``Max 2``, ``Max 3``
   Exclusive upper grid bounds in the same units as the minima.

``Step 1``, ``Step 2``, ``Step 3``
   Voxel widths. Editing a step recalculates its interval count.

``Intervals 1``, ``Intervals 2``, ``Intervals 3``
   Editable numbers of voxels along each axis. Editing an interval count
   recalculates the corresponding step from the current bounds.

``Est. size``
   Dense uncompressed payload estimate:

   .. math::

      N_1 N_2 N_3 \times 32\ \mathrm{bytes}.

   The 32 bytes comprise float64 intensity, variance, and weight plus uint64
   contributors. Sparse HDF5 chunk allocation and compression usually reduce
   payload size; axes, metadata, and HDF5 structures add a smaller overhead.

Below the grid table, an estimated checkpoint scratch-file and cluster-job
count is shown, from a short live calibration probe against the active scan
and detector geometry. It reports the single-node file count (with a
per-grid breakdown when more than one grid is defined) and, if the Cluster
tab's node count is greater than one, a second estimate for that many
cluster nodes -- the same prepared job can still be run either locally or
as a cluster job. This updates automatically as grids, accuracy, angle
fallback, thread/memory budgets, ``Checkpoint count``, or the cluster node
count change; it requires an active scan and at least one output grid, and
is an estimate, not a guarantee -- the checkpoint-count floor and available
memory may still produce more files than shown.

Adding and Removing Grids
~~~~~~~~~~~~~~~~~~~~~~~~~

``Add derived HKL grid``
   Adds a grid covering the active detector and scan in HKL.

``Add derived Q grid``
   Prompts for a Q reference frame and adds its derived coverage.

``Remove selected grid``
   Removes every row containing a selected table cell.

Automatic Coverage and Initial Steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coverage derivation transforms all four detector corners at the start and end
of every exposure. Per-axis minima and maxima bound those transformed points.
The initial sampling count is:

.. code-block:: python

   max(64, min(512, max(detector_rows, detector_columns, frame_count)))

The initial step is the axis extent divided by that scalar count. One step of
padding is added on both sides, so the resulting grid has approximately two
more intervals. The rule is a coverage heuristic, not an instrumental
resolution model.

Geometry-Matched Step Estimator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Estimate geometry-matched steps`` applies to selected grid rows, or all rows
when nothing is selected. It samples a 5 by 5 set of detector pixels and up to
32 included scan frames. At each point it evaluates finite-difference
generators for one detector-row pixel, one detector-column pixel, and the
exposure sweep. A stationary exposure uses the adjacent-frame center
displacement for the third direction.

With those generators as columns of a local Jacobian :math:`J`, uniform pixel
and scan-cell widths give:

.. math::

   \Sigma_y = \frac{J J^\mathsf{T}}{12}.

The local axis estimates are
:math:`\sqrt{\operatorname{diag}(\Sigma_y)}`. The estimator dialog exposes the
percentile applied independently to the sampled axis estimates. The default is
10 percent. Lower percentiles select finer steps and protect more locally
high-resolution regions; higher percentiles produce smaller grids.

This is a geometry-matched sampling estimate. It deliberately excludes beam
divergence, energy bandwidth, detector point-spread, sample mosaicity, and
calibration covariance.

HDF5 Settings
-------------

``HDF5 settings`` is available from the grid tab and the main Configuration
menu. These values are shared by all reconstruction grids.

``Chunk shape (voxels)``
   Three spatial HDF5 chunk lengths. A chunk is the unit written by the
   finalizer. Every populated chunk is written once; untouched chunks remain
   unallocated. One float64 dataset chunk must remain below the HDF5 4 GiB
   limit.

``Active database filter``
   Compression used by normal orGUI database writes and, by default, by the
   reconstruction.

``Override for reconstruction output``
   When enabled, ``Reconstruction filter`` replaces the active database filter
   only for prepared reconstruction output. The selected registered filter is
   frozen into the job.

Accuracy Parameters
-------------------

Each detector pixel and exposure sweep is treated as a parametric cell.
Adaptive splitting tests transformed cell corners against voxel boundaries.
A leaf wholly inside one voxel contributes its full weight. Otherwise it is
subdivided until it becomes single-voxel or reaches the selected depth, where
its centroid is assigned. Total source-pixel weight is conserved.

``Center only (depth 0)``
   Assign one centroid per pixel. Fastest, without finite-footprint splitting.

``Low (depth 1)``
   One adaptive subdivision level.

``Balanced (depth 2)``
   Default compromise between boundary accuracy and runtime.

``High (depth 3)``
   Finer boundary subdivision at substantially higher cost.

``Very high (depth 4)``
   Intended for validation or demanding coarse-voxel boundaries.

``Maximum (depth 5)``
   Highest UI setting. Benchmark memory and runtime before large production
   jobs.

Stationary cells split into four spatial children per level. Swept cells split
in detector row, detector column, and exposure position and can create eight
children per level. Consequently, worst-case work grows rapidly with depth.

Performance Parameters
----------------------

Detected values appear in ``Detected execution layout``. An unchecked
``Override`` editor displays the detected value but does not freeze it into the
job.

``Total thread budget``
   Optional per-job replacement for ``orGUI.numberthreads``. It is divided
   between concurrent image workers and native threads within each image.

``Native threads per image``
   C++ threads assigned to one image. Remaining threads can run other images
   concurrently. By default (``Override`` unchecked) this is chosen
   automatically: mapping starts image-parallel (one native thread per image)
   and periodically rebalances native threads per image against concurrent
   images while running, using the real, currently measured image delivery
   rate for this job's storage -- not a value that can be predicted or
   configured ahead of time. Check ``Override`` to pin a fixed value for the
   whole job instead (the previous, always-static behavior). Either way the
   actual value is bounded by the total thread and memory budgets.

Mapping overlaps image loading with image processing: a small prefetch-reader
pool loads images ahead of the compute workers that process them, growing and
shrinking live based on how often compute sits waiting for an image. Progress
messages report the current prefetch reader count alongside concurrent image
workers and native threads per image.

``Total memory budget``
   Optional per-job replacement for ``orGUI.maxMemory``, in MiB. It bounds
   image workers, native working memory, retained records, sorting, and Arrow
   conversion.

``Accumulation per worker``
   Optional conservative per-image-worker memory margin used only to size
   ``concurrent image workers`` (larger values leave more headroom and admit
   fewer concurrent workers). It no longer sizes a Parquet flush buffer -- the
   checkpoint scratch format has no per-worker segment file of its own to
   tune.

``Frames per task``
   Consecutive included frames in one deterministic resumable map task.
   Automatic sizing targets several tasks per concurrent image worker and caps
   a task at 64 frames.

``Tile rows`` and ``Tile columns``
   Detector dimensions passed to one native-kernel call. Tiles bound image,
   variance, mask, and corner-ray arrays and may end at arbitrary pixel
   boundaries. All tiles of one loaded image are processed without reloading
   or reapplying full-image corrections.

``Native block pixels``
   Flattened pixels in one C++ scheduling and local-reduction block inside a
   detector tile. It is not the tile area. For example, a 1024 by 1024 tile
   and a 65,536-pixel block create 16 blocks. Multiple blocks are necessary
   for native threads to share a tile; very small blocks add scheduling,
   sorting, and merge overhead.

``Checkpoint count``
   Minimum number of resumable HDF5 checkpoint files per output grid,
   default 10. A floor, not a target: the actual count exceeds it whenever
   the estimated data volume would not otherwise fit the memory budget
   across that many files. Checkpoint frame-range boundaries are computed
   once from a calibration estimate when the job is prepared and never
   recomputed on resume.

The detected summary also reports frame tasks, detector tiles, total map tasks,
concurrent image workers, memory per image, and effective accumulation memory.

Correction and Variance Model
-----------------------------

Corrections are evaluated once for each loaded full image, before detector
tiling:

1. Use variance supplied by the image backend. If unavailable, use the flagged
   approximation :math:`v=\max(I,0)` from the current raw image.
2. Subtract the active background. Add background variance when supplied;
   otherwise record that the background was treated as deterministic.
3. Apply the active detector mask and optional pixel repair.
4. Apply inverse solid-angle and inverse polarization factors.
5. Normalize by exposure time when enabled and available.
6. Apply explicitly selected divisive monitor counters.
7. Mask non-finite corrected intensity or variance.

For a deterministic multiplicative factor :math:`c`:

.. math::

   I' = cI,\qquad v' = c^2v.

When factor variance :math:`v_c` is available:

.. math::

   v' = c^2v + I^2v_c.

The constant detector mask is analyzed once. Repairable components and their
neighbors are stored in a native repair plan reused for all images. Isolated
pixels use the local valid-neighbor median; accepted small components use
inverse-distance weights. Interpolation variance uses the squared weights.
Repair introduces spatial covariance that is not stored in the final file.

Voxel Accumulation
------------------

Leaves from the same source pixel entering the same voxel are combined before
global accumulation, preserving their correlation. For independent source
pixels the kernel accumulates:

.. math::

   S_I = \sum wI,\qquad
   S_v = \sum w^2v,\qquad
   S_w = \sum w.

Final voxel datasets are:

.. math::

   I_\mathrm{voxel} = \frac{S_I}{S_w},\qquad
   v_\mathrm{voxel} = \frac{S_v}{S_w^2}.

Empty intensity and variance voxels are NaN. ``weight`` stores :math:`S_w`;
``contributors`` stores the number of independent source pixels. Adaptive
footprint splitting creates cross-voxel covariance, so ``variance`` contains
marginal variances only.

Out-of-Core Execution
---------------------

The execution stages are:

1. Verify the job JSON, scan reference, source fingerprint, and asset checksum.
2. Precompute detector corner-ray lattices for bounded tiles.
3. Divide included frames into a small number of contiguous per-grid
   checkpoint frame ranges, sized from a live calibration estimate of data
   volume versus the memory budget (a floor of ``Checkpoint count`` files,
   exceeded only when data volume requires it). Boundaries are computed once
   when the job is prepared and persisted in the job JSON.
4. Load and correct images in parallel.
5. Map tiles in the native kernel, releasing the Python GIL, and route each
   frame's reduced records directly to the checkpoint accumulator for the
   grid and frame range it belongs to.
6. Merge each frame's records into its checkpoint's in-memory tree as it
   arrives (an amortized :math:`O(\log N)` insert per frame, not a full
   re-sort). When a checkpoint's frame range finishes -- or its own share of
   the memory budget is exceeded first, in which case it splits into an
   extra file -- write it as an immutable, checksummed HDF5 checkpoint file.
7. Once every planned checkpoint is present and its checksum verified,
   stream each output grid's populated spatial chunks directly from its
   checkpoint files (merging any records checkpoints have in common),
   writing each chunk to the standalone HDF5 file exactly once.
8. Validate the output checksum, register it as an external result, then
   remove the checkpoint files and the asset bundle.

Interrupted jobs retain checkpoint data. Cleanup occurs only after successful
HDF5 close, validation, and checksum calculation. Cleanup errors are recorded
but do not invalidate a completed scientific result.

Resuming a job re-verifies each planned checkpoint's presence and checksum
directly -- there is no separate manifest registry to consult. A planned
checkpoint already fully covered by valid, checksummed files on disk is
skipped; anything else is remapped from scratch.

Paths
-----

``Job JSON``
   Resumable descriptor containing frozen scientific settings, paths, build
   metadata, task state, and checksums.

``Scratch directory``
   Writable per-job directory for the asset bundle and checkpoint files. It
   should normally be on fast local storage.

``Standalone HDF5``
   Final NeXus-style output file. Defaults are created under the current
   working directory because raw-data locations are often read-only.

Command-Line Execution
----------------------

The CLI consumes the exact job prepared by the UI; it does not provide a
parallel experiment-configuration interface:

.. code-block:: bash

   orGUI rsmap run JOB.json
   orGUI rsmap resume JOB.json
   orGUI rsmap status JOB.json

The direct alias is equivalent:

.. code-block:: bash

   orGUI-rsmap run JOB.json
   orGUI-rsmap resume JOB.json
   orGUI-rsmap status JOB.json

``run`` and ``resume`` both verify and continue deterministic task state.
``status`` is read-only and prints the current descriptor status as JSON.

Cluster Batch Execution
-----------------------

The **Cluster** tab generates an SGE or Slurm batch bundle from the same
prepared job JSON used for local execution. It does not introduce another
experiment configuration. Each array element runs the full single-node
mapping pipeline against its own disjoint, equal share of the scan's
frames, computed purely from its position in the array (``task_index``)
and the array's total size (``total_tasks``) -- not from the scan's frame
count or any grid structure. No node communicates with any other node
while mapping; a single dependent finalizer job merges every node's
checkpoint files once all of them are complete. The bundle contains:

``*-map.sge`` or ``*-map.slurm``
   A job array with one element per requested node. Each element reads the
   immutable job and asset bundle, computes its own frame slice, and
   writes only its own checkpoint files under a node-specific scratch
   subdirectory. Array elements never update the shared job JSON, so they
   may execute concurrently.

``*-finalize.sge`` or ``*-finalize.slurm``
   A single dependent job. It verifies that every node's planned
   checkpoints are present and checksummed, merges their checkpoint files
   directly into the final HDF5 file, and records completion in the job
   JSON.

``*-submit.sh``
   A convenience submission wrapper. For SGE it captures ``qsub -terse``
   and submits the finalizer with ``-hold_jid``. For Slurm it captures
   ``sbatch --parsable`` and uses ``afterok`` on the complete array. The
   finalizer always verifies every node, including on SGE installations
   where a job hold records completion rather than successful exit.

The scripts require the job JSON, scratch directory, immutable assets, scan
source, and final output directory to be reachable at the same absolute paths
from every compute node. Scratch should normally reside on a high-throughput
shared filesystem or node-local storage explicitly staged by site-specific
setup commands.

**Array size is never read from the scheduler at run time.** Most
schedulers besides Slurm do not reliably expose a running task's total
array size as an environment variable (SGE, PBS Pro, and LSF expose only
each task's own index). To stay portable across schedulers, orGUI never
relies on one: the array size chosen in the **Number of nodes** field is
baked directly into the generated ``cluster-map``/``cluster-finalize``
commands as an explicit ``--total-tasks`` argument, the same way
``--cpus``/``--memory-gib`` already are. Only each array element's own
index comes from the scheduler (``SGE_TASK_ID``, ``SLURM_ARRAY_TASK_ID``),
normalized to a 0-based ``task_index`` by the generated script.

The generated worker commands are also available for inspection and manual
scheduler integration:

.. code-block:: bash

   orGUI rsmap cluster-map JOB.json --task-index 0 --total-tasks 8 --cpus 4 --memory-gib 16
   orGUI rsmap cluster-finalize JOB.json --total-tasks 8 --cpus 24 --memory-gib 64
   orGUI rsmap cluster-scripts JOB.json --output-directory batch

Scheduler Parameters
~~~~~~~~~~~~~~~~~~~~

``Scheduler``
   ``SGE`` (default) or ``Slurm``. SGE arrays are one-based and converted to
   orGUI's zero-based task index. Slurm arrays are generated zero-based.

``Job name``
   Scheduler-safe base name for the map array and finalizer. Letters, numbers,
   periods, underscores, and hyphens are accepted.

``Queue / partition``
   Optional SGE ``-q`` queue or Slurm ``--partition``.

``Project / account``
   Optional SGE ``-P`` project or Slurm ``--account``.

``Script directory``
   Destination for the three generated scripts. It defaults beneath the
   current writable working directory.

``Working directory``
   Shared directory selected with ``cd`` before environment setup and Python
   execution.

``Python executable``
   Python command used to invoke ``orgui.reconstruction_cli`` after environment
   setup. ``python`` is the portable default; an absolute cluster-environment
   path can be supplied.

``Setup commands``
   Verbatim shell lines placed after ``set -euo pipefail`` and ``cd``. Use
   these for module loading and conda or virtual-environment activation.

``Number of nodes``
   How many array elements (nodes) to request. Each node maps an equal,
   disjoint share of the scan's frames. See the array-size note above --
   this value is never inferred from the scheduler.

``CPUs / slots per mapping task``
   Native C++ threads used by each array element. This allocation is
   independent of the number of simultaneously running array elements.

``Memory per mapping task``
   Total RAM budget passed to one mapping process. Slurm receives it through
   ``--mem``. SGE memory complexes are normally per-slot, so orGUI divides the
   total by the slot count and rounds upward.

``Mapping wall time``
   Per-element SGE ``h_rt`` or Slurm ``--time`` limit.

``Maximum concurrent tasks``
   Optional array throttle: SGE ``-tc`` or the Slurm ``%N`` array suffix.
   Zero leaves concurrency to site policy.

``Reduction CPUs / slots``
   Accepted for symmetry with the mapping allocation; the finalizer's own
   chunk-merge loop is currently single-threaded, matching the single-node
   finalize path, so this is presently informational rather than sized
   against an internal worker pool.

``Reduction memory``
   Informational RAM budget for the finalizer process. Finalize already
   bounds memory per checkpoint file via its own streaming range reader
   regardless of this value, matching the single-node path.

``Reduction wall time``
   SGE ``h_rt`` or Slurm ``--time`` for the dependent job.

``SGE parallel environment``
   Name requested through ``-pe``; ``smp`` is the default. It must match the
   target site's configured shared-memory parallel environment.

``SGE memory resource``
   Per-slot consumable complex used in ``-l`` requests. ``h_vmem`` is the
   default, but sites may require ``mem_free`` or another name.

``Extra map/finalizer directives``
   Optional scheduler-specific header lines. Each non-empty line must start
   with ``#$`` for SGE or ``#SBATCH`` for Slurm.

See the `Grid Engine qsub reference
<https://gridengine.eu/mangridengine/htmlman1/qsub.html>`_ and the `Slurm job
array reference <https://slurm.schedmd.com/job_array.html>`_ for scheduler
semantics.

Job Descriptor Reference
------------------------

Users normally create this file through the UI. The fields are documented for
inspection, scheduling, and provenance:

``schema_version``
   Exact reconstruction job schema. Unsupported schemas are rejected.

``config``
   Authoritative central ``ConfigData`` snapshot containing detector
   calibration, crystal, UB, diffractometer values, and correction state.

``scan_reference``
   Serializable reference for the active scan, including supported slices,
   interlaced scans, manual scans, simulations, and external backend files.

``grids``
   List of output grid dictionaries: minimum, maximum, step, frame, name, and
   chunk shape.

``scratch_path``, ``output_path``
   Absolute local paths described above.

``compression``
   Frozen name from the central HDF5 filter registry.

``assets_path``, ``assets_sha256``
   Immutable job asset bundle and SHA-256 checksum.

``source_fingerprint_sha256``
   Digest of the serialized scan reference used to detect changed sources.

``build_metadata``
   orGUI, NumPy, h5py, compiler, and native build information.

``runtime_threads``, ``runtime_memory_bytes``
   Global orGUI defaults captured when the job was prepared.

``thread_override``, ``memory_override_bytes``
   Optional per-job replacements for the captured defaults.

``threads_per_image``
   Requested native threads per concurrent image, or ``null`` for automatic
   (the default): chosen at run time and periodically rebalanced live against
   the measured image delivery rate, instead of a single fixed value.

``accumulation_budget_bytes``
   Optional retained-record bytes per image worker.

``angle_fallback``
   ``stationary`` or explicitly selected ``midpoint`` inference.

``accuracy``
   Named footprint-depth preset.

``advanced_depth``
   Internal legacy field in the current descriptor schema. UI-prepared jobs
   leave it null and use ``accuracy``.

``frame_batch``
   Optional frames-per-task override.

``tile_shape``
   Optional detector tile rows and columns.

``work_block_pixels``
   Optional native scheduling-block override.

``checkpoint_count``
   Minimum number of resumable HDF5 checkpoint files per output grid. A
   floor, not a target: the actual count exceeds it when the estimated
   data volume would not otherwise fit the memory budget across that many
   files.

``checkpoint_plan``
   Per-grid list of contiguous frame-range checkpoint boundaries, computed
   once from a calibration estimate when the job is prepared and never
   recomputed on resume.

``user_note``
   Optional free-text provenance.

``status``
   Current job state, such as ``prepared``, ``mapping``, ``finalizing``,
   or ``complete``.

``output_sha256``
   Final standalone-file checksum after verified completion.

``correction_provenance``
   Recorded variance source, factor-uncertainty assumptions, repair-plan
   configuration, and image-processing provenance.

``cleanup_errors``
   Nonfatal failures encountered while deleting verified scratch outputs.

``cluster_settings``
   Frozen scheduler, environment, map-array resource, reduction-resource, and
   optional directive settings used to regenerate the batch bundle.

Output File Layout
------------------

The final file contains an ``NXentry`` at ``/entry``. The reconstruction
``NXprocess`` stores the frozen compute configuration, scientific context,
provenance, marginal-variance warning, and central orGUI scan configuration.

Each selected coordinate system is an ``NXdata`` group below
``/entry/reconstruction/results``. It contains:

* voxel-center axes ``h``, ``k``, ``l`` in r.l.u. or ``qx``, ``qy``, ``qz`` in
  :math:`\mathrm{\mathring{A}}^{-1}`;
* float64 ``intensity``;
* float64 marginal ``variance``;
* float64 ``weight``;
* uint64 ``contributors``;
* coordinate-frame, signal, axes, and units attributes.

All scientific arrays default to float64 except the contributor count. The
file is conventional HDF5 and can be read with h5py, silx, and NeXus-aware
tools.

Practical Guidance
------------------

* Start with ``Balanced`` depth and the 10-percent geometry estimator.
* Check interval counts and the dense size estimate before preparing.
* Put scratch data on fast local SSD when possible.
* Keep enough native blocks per tile to occupy all native threads.
* Prefer automatic tiling and task sizing until representative benchmarks
  justify overrides.
* Use ``status`` before resuming a job copied from another machine.
* Treat the variance arrays as marginal uncertainties, not a complete
  covariance representation.

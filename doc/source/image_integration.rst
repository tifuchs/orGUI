Image Integration
=================

orGUI integrates detector images by calculating detector coordinates for the
requested reciprocal-space feature, placing a region of interest (ROI) at that
position, and summing the detector intensity with optional background and
correction handling. Integration results are written to the active Nexus
database file.

Stationary Reciprocal-Space Integration
---------------------------------------

The ``hklscan`` tab integrates a stationary scan along a reciprocal-space line:

.. math::

   \vec{H}(s) = \vec{H}_0 + s\vec{H}_1

The vectors :math:`\vec{H}_0` and :math:`\vec{H}_1` are entered in reciprocal
lattice units. For each image, orGUI calculates the intersections of this line
with the Ewald sphere and converts the valid intersections to detector
coordinates. The two possible intersections are handled as separate S1 and S2
trajectories.

This mode is the usual choice for crystal truncation rods in stationary
rotation scans. To preview the calculated ROI positions, enable
``View -> Show ROI``.

Fixed-Pixel Integration
-----------------------

The ``fixed`` tab integrates a detector ROI at a fixed pixel position through
the scan. The same detector coordinates are used for each image, while orGUI
records the corresponding reciprocal-space coordinates and diffractometer angles
for the active scan state.

This mode is useful when the desired detector region is selected directly from
the image rather than from a reciprocal-space trajectory.

Rocking Reciprocal-Space Integration
------------------------------------

The ``rocking hklscan`` tab integrates multiple ROIs through a rocking scan.
The ROI centers are sampled along the reciprocal-space line
:math:`\vec{H}_0 + s\vec{H}_1`; each sampled coordinate is converted to a
detector position and integrated across the rocking scan.

Rocking integration is often a better choice for quantitative CTR extraction
because the rocking dimension helps separate the CTR signal from broad or
environmental background. It is more computationally expensive than stationary
integration.

ROI summing uses the C++ acceleration backend by default. See
:doc:`acceleration_backends` for the process-global backend selector and the
optional Numba ROI backend.

Rocking Bragg Integration
-------------------------

The ``rocking Bragg`` tab calculates detector coordinates for Bragg peak
positions from the current crystal, UB matrix, detector calibration, strain, and
scan state. Valid Bragg positions are then integrated as rocking-scan ROIs.

Masks, Backgrounds, and Corrections
-----------------------------------

The image-view mask tool can be used to mask detector pixels before
integration. Center ROIs and background ROIs are defined from the ROI controls
on the scan selector. Background subtraction is performed from neighboring
background regions or from a compatible background image when configured.

Solid-angle and polarization correction factors can be applied during
integration. The current implementation records the applied ROI sizes,
reciprocal-space coordinates, detector coordinates, and relevant scan metadata
with the integrated intensities.

Detector Masks
~~~~~~~~~~~~~~

Detector masks mark invalid pixels. In orGUI, mask values are interpreted as:

* zero or ``False``: valid detector pixel
* nonzero or ``True``: invalid detector pixel

Masks can be created or edited interactively with the image-view mask tool.
They can also be loaded from the config file with ``[Mask]``. Config masks are
loaded with FabIO, so EDF, NumPy ``.npy``, and other FabIO-readable 2D image
formats are supported.

When ``Use pixel mask`` is enabled for integration, masked pixels are excluded
from center ROI sums, background ROI sums, background-image counters, fitted
background samples, and correction counters unless pixel repair explicitly
repairs a tiny signal-ROI defect as described below.

.. _pixel-repair-algorithm:

Pixel Repair Algorithm
~~~~~~~~~~~~~~~~~~~~~~

Pixel repair is an optional, conservative correction for tiny masked defects
inside signal ROIs. It is configured by ``[Mask.PixelRepair]`` and is disabled
unless that section exists and ``enabled = True``. It requires the C++ ROI
backend. If the C++ backend is unavailable, orGUI logs a warning and continues
without repair.

Repair is intended for isolated dead pixels or very small clusters. It is not
image-wide inpainting, spline fitting, or detector-gap interpolation.

Candidate Selection
^^^^^^^^^^^^^^^^^^^

The repair code first treats the boolean detector mask as connected components
of masked pixels. A component is rejected when any of the following is true:

* the component intersects or touches an analytical detector gap interval;
* the component contains more than ``max_component_pixels`` pixels;
* the component row span or column span exceeds ``max_span`` pixels;
* fewer than ``min_valid_neighbors`` original-valid neighboring pixels are
  found within ``radius`` pixels;
* the valid neighbors are not distributed on enough sides of the defect.

Only original-valid neighboring pixels are used for interpolation. Other masked
pixels and ``NaN`` image values are ignored.

Detector Gap Exclusion
^^^^^^^^^^^^^^^^^^^^^^

When ``use_pyfai_gaps = True``, orGUI excludes analytical detector gaps from
repair. The pyFAI detector mask is reduced to full-width row gaps and
full-height column gaps, stored as ``int32`` ``[start, stop)`` intervals.
Repair candidates that intersect or touch one of these intervals are rejected.

If pyFAI cannot provide a useful analytical gap mask, orGUI falls back to
detector module-size metadata when available. In that metadata fallback path,
``gap_size_px`` supplies the detector gap width only when the detector metadata
do not provide one. It is not used to grow or dilate arbitrary masked pixels.

Interpolation
^^^^^^^^^^^^^

Accepted defects are filled from local original-valid neighbors:

* isolated one-pixel defects use the local valid-neighbor median;
* tiny multi-pixel clusters use an inverse-distance weighted mean.

No polynomial, spline, cross-gap, or image-wide interpolation is performed.

Counter Accounting
^^^^^^^^^^^^^^^^^^

Repaired pixels are counted only where their use is scientifically narrow and
explicit:

* repaired pixels may contribute to the signal ROI raw intensity sum;
* repaired pixels may contribute to the finite correction sum for the signal
  ROI;
* repaired pixels do not contribute to background ROI sums;
* repaired pixels do not contribute to background-image counters;
* repaired pixels do not contribute to fitted-background samples;
* repaired pixels do not alter fitted-background preview images.

For fitted local background integration, repair is currently disabled and orGUI
logs a warning. The original mask is used for the fitted-background path.

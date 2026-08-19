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

.. warning::

   A ***critical bug*** was fixed that affects rocking-scan ROI integration.
   Versions up to and including v1.5.0
   under-subtracted the background when more than one one-dimensional
   background ROI was used (including the default ``Add All`` configuration).
   They also used inconsistent ROI endpoint weights and, in some cases,
   incorrect propagated uncertainties. Consequently, previously saved rocking
   ``croibg`` and ``F2_hkl`` intensities and their errors are incorrect. The effect is most
   severe for background-dominated reflections and integrations with multiple
   background ROIs (underestimates the background).

   To recover an affected result, repeat the final rocking integration with the
   corrected version. The Nexus database retains the source rocking ROI curves
   (``rois/croibg`` and ``rois/croibg_errors``) and the one-dimensional ROI
   definitions, so the detector images do not need to be integrated again when
   those source datasets are still available. Retain or clearly label the old
   measurement rather than treating it as comparable with the regenerated one.

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

The Corrections Dialog
~~~~~~~~~~~~~~~~~~~~~~

All corrections are configured in one non-modal dialog, opened with
``Corrections and normalization ...`` on the ``ROI integration`` tab. It stays
open while you work, because an integration is often repeated with a
correction switched on or off.

Below the button, a status line abbreviates the enabled corrections in colour:
``MASK`` (pixel mask), ``SOLA`` (solid angle), ``POL`` (polarization),
``LOR`` (Lorentz factor and rod interception), ``FOOT`` (beam footprint) and
``NORM`` (exposure and monitor normalization). With nothing enabled it reads
*no corrections applied*. Hovering the line explains each abbreviation.

Detector masks are still drawn with the mask tool of the image view, and the
reciprocal-space reconstruction keeps its own settings dialog; both are
reachable from buttons inside the corrections dialog, and the mask, solid
angle, polarization and normalization settings are shared between the two
workflows rather than duplicated.

Geometrical Corrections
~~~~~~~~~~~~~~~~~~~~~~~

``Lorentz and rod interception`` divides out the Lorentz factor :math:`L` and
the rod interception :math:`C_\mathrm{rod}`, and stores the structure factor

.. math:: F^2_{hkl} = \frac{I_\mathrm{corr}}{L \cdot C_\mathrm{rod}}

The factors are those tabulated for the **z-axis geometry** in Appendix A of
the ANA/ROD manual (E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532), with
:math:`\alpha` the incidence angle (the ``mu`` circle), :math:`\delta` the
in-plane and :math:`\gamma` the out-of-plane detector angle:

.. math::

   L_\mathrm{rocking} &= \frac{1}{\sin\delta\,\cos\alpha\,\cos\gamma} \\
   L_\mathrm{reflectivity} &= \frac{1}{\sin 2\alpha} \\
   L_\mathrm{stationary} &= \frac{1}{\sin\gamma} \\
   C_\mathrm{rod} &= \cos\gamma

The three Lorentz factors are alternatives, not factors to be combined, and
orGUI selects between them by how the scan was measured:

* a **rocking scan** about the sample rotation (``th``) uses
  :math:`L_\mathrm{rocking}`,
* a **rocking scan** about the incidence angle (``mu``), which is how a
  reflectivity curve is measured, uses :math:`L_\mathrm{reflectivity}`,
* a **stationary scan** (``hklscan`` and ``fixed``), integrated across the rod
  on the area detector, uses :math:`L_\mathrm{stationary}`.

Using the rocking-scan factor for a stationary scan, or the reverse, is a
silent scaling error, which is why the choice follows the integration mode
rather than a user setting.

The manual also lists an area correction :math:`1/\sin\delta` "ignoring
footprint and sample size". orGUI does not apply it: the footprint corrections
below evaluate the beam profile and the finite sample size numerically
instead, which is the row the manual marks as calculated numerically.

Exposure and Monitor Normalization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Normalize integrated intensities`` divides every image by its exposure time
and by each configured monitor counter, the same normalization the
reciprocal-space reconstruction applies, so an integration and a
reconstruction of one scan are normalized identically. The settings are
shared: changing them here changes them for the reconstruction and the other
way round.

A scan backend that provides no exposure time is not an error; the exposure
part of the normalization is then skipped. A monitor counter that is missing,
zero or non-finite does stop the integration, since it would otherwise scale
intensities by infinity. The dialog lists the counters of the loaded scan that
can be used.

Footprint Corrections
~~~~~~~~~~~~~~~~~~~~~

``Beam footprint`` applies two corrections that depend on the incidence angle
:math:`\alpha` and on the
vertical profile :math:`p(z)` of the incident beam, normalized to
:math:`\int p(z)\,\mathrm{d}z = 1`. With a sample of length :math:`L` along the
beam, its projected size across the beam is :math:`h(\alpha) = L\sin\alpha`,
and with the sample centered at :math:`z_0`:

.. math::

   C_\mathrm{flux} = \int_{z_0-h/2}^{z_0+h/2} p(z)\,\mathrm{d}z
   \qquad
   C_\mathrm{area} = \frac{C_\mathrm{flux}}{p_\mathrm{max}\,h}

``C_flux_on_sample`` is the fraction of the incident flux that strikes the
sample, correcting the beam overspill at grazing incidence. ``C_illum_area``
is the mean of :math:`p(z)/p_\mathrm{max}` over the projected sample
footprint, the active surface area seen by the measurement; it tends to 1 when
the sample is small compared with the beam and falls off as
:math:`1/\sin\alpha` once the beam is fully on the sample. Integrated
intensities are divided by the product of the two, and both factors are stored
alongside the result.

Only the cumulative integral of :math:`p(z)` and its peak value enter these
expressions, so any distribution can describe the beam. The ``options``
button next to the checkbox chooses between an analytical shape and a
measured profile.

Analytical Beam Shapes
~~~~~~~~~~~~~~~~~~~~~~

All widths are entered in micrometers.

``Gaussian``
   The default, and numerically identical to the correction orGUI applied up
   to and including v1.5.0.

``Top hat``
   A beam of uniform intensity, as defined by a slit. ``C_flux_on_sample``
   is then ``L sin(alpha) / width`` until the sample intercepts the whole
   beam, and ``C_illum_area`` stays at 1 over the same range. This is the
   linear geometrical correction of Gibaud, Vignaud & Sinha (1993),
   *Acta Cryst.* A49, 642, equation (12).

``Trapezoid``
   Two slits of different aperture, entered as the base and flat-top widths.
   A flat top of zero gives the triangle of two matched slits.

``Smoothed top hat``
   A slit-defined beam whose edges are blurred by the finite source size:
   flat in the middle with error-function flanks of the given sigma.

``Generalized normal``
   Density proportional to ``exp(-|z / scale| ** flatness)``. A flatness of
   2 is exactly Gaussian, 1 gives a cusp, and large values approach a top
   hat. This is the one-parameter family for a focused beam that is neither
   Gaussian nor flat.

``Skew normal``
   An asymmetric profile with a Gaussian core, for a beam whose two flanks
   differ. A skew of 0 is Gaussian.

Each shape is scaled to the requested full width at half maximum, so
switching between them changes the shape rather than the width.

Measured Beam Profiles
~~~~~~~~~~~~~~~~~~~~~~

``measured beam profile`` reads a profile measured at the beamline from a
two-column text file and integrates it numerically. Use it when the beam is
asymmetric or has more than one maximum, which no analytical shape with one
or two parameters can represent.

Measuring the Beam Profile
~~~~~~~~~~~~~~~~~~~~~~~~~~

A sample height scan measures the beam profile directly. As the sample edge
cuts into the beam, the transmitted intensity is the beam profile integrated
over the part of the beam that still passes the sample, so the profile is the
negated derivative of the measured curve:

.. math::

   I(z) = I_0 \int_z^\infty p(z')\,\mathrm{d}z'
   \qquad\Longrightarrow\qquad
   p(z) \propto -\frac{\mathrm{d}I}{\mathrm{d}z}

Differentiating assumes no beam shape at all, unlike fitting a width to the
transmitted curve, and the constant offset of a partially transparent or
background-affected scan drops out. Because the sample itself is the cutting
edge and the diode integrates the whole transmitted beam, the result needs no
deconvolution of a detector aperture.

Set ``file contains`` to ``height scan (-dI/dz)`` to have orGUI differentiate
a height scan on loading; it then keeps only the range over which the sample
cuts into the beam, discarding the mirrored second edge a scan that moves the
sample fully through the beam records afterwards. Set it to ``beam profile``
for a file that already holds a profile.

``examples/beam_profiles`` contains a measured profile in the expected
format.

Placing the Sample in the Beam
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A profile carries no absolute position, so ``sample centered on`` selects the
point of the beam that the center of the sample sits on: the intensity
centroid, the maximum, or the median (the half-cut position an edge-scan
alignment converges to). All three coincide for a symmetric beam, but for an
asymmetric one they give measurably different corrections. ``sample offset``
displaces the sample further, in micrometers, and applies to analytical
shapes as well, which is how an alignment error is modeled.

The preview at the bottom of the dialog plots whichever profile is currently
described, normalized and against position relative to the sample center. A
black line marks that center and a green line the centroid of the beam, so
the effect of the centering choice is visible directly: the two coincide for
a symmetric beam and separate for an asymmetric one. The line below the plot
reports the full width at half maximum, the rms width, the extent of the
profile and the centroid offset.

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

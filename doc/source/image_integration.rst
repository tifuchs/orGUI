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

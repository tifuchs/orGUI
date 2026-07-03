Setting the UB Matrix
=====================

The UB matrix describes how reciprocal lattice coordinates are oriented on the
diffractometer. The ``B`` matrix is determined by the lattice constants and
converts ``(h, k, l)`` from reciprocal lattice units to Cartesian
reciprocal-space coordinates. The ``U`` matrix describes how the crystal is
mounted relative to the diffractometer sample circles.

orGUI provides two complementary workflows:

``Manual reflection setting``
   The user selects known reference reflections in detector images, assigns
   ``hkl`` values, and calculates ``U`` from one or more measured reflections.
   This is the most transparent workflow and is the best way to diagnose
   detector, scan-axis, or indexing problems.

``Automatic UB setting``
   orGUI streams through the scan images, detects a sharp Bragg-like peak,
   assigns candidate ``hkl`` values, builds a one-reflection trial ``U`` matrix,
   confirms that matrix with another predicted Bragg peak, and optionally adds
   more validated reflections.

.. toctree::
   :maxdepth: 2

   manual_reflection_setting
   automatic_ub_setting

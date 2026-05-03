Setting the UB Matrix
=====================

The UB matrix describes how reciprocal lattice coordinates are oriented on the
diffractometer. The ``B`` matrix is determined by the lattice constants and
converts ``(h, k, l)`` from reciprocal lattice units to Cartesian
reciprocal-space coordinates. The ``U`` matrix describes how the crystal is
mounted relative to the diffractometer sample circles.

Reference-Reflection Method
---------------------------

The usual workflow is to determine ``U`` from one or more observed reference
reflections:

1. Open a scan and display an image that contains a known reflection.
2. Double-click the reflection position in the image. This creates a reference
   marker in the ``Reciprocal space navigation`` window.
3. Move the marker if needed. The marker position defines the detector
   coordinates, and the selected image number defines the sample rotation angle
   for that observation.
4. Enter or correct the reflection ``h``, ``k``, and ``l`` values in the
   reference-reflection table.
5. Repeat for a second linearly independent reflection when possible.
6. Click ``Calculate U``.

Two linearly independent reflections fully determine the orientation matrix.
For the common surface-diffraction case, orGUI can also estimate an orientation
from a single reference reflection by assuming that the ``L`` direction points
toward the ``z``/azimuth reference direction. This is useful for quick setup
when a clear reflection near ``L = 0`` is visible, but a two-reflection
orientation is the more general procedure.

Expert Matrix Editing
---------------------

The menu entry ``Reciprocal Space -> Edit orientation matrix`` opens the expert
path for editing the matrix directly. This is intended for cases where the
orientation can be inferred from known symmetry or from an already established
experiment geometry. It is less forgiving than the reference-reflection method,
especially in grazing-incidence geometries.

Validation
----------

After calculating the matrix, enable ``View -> CTR reflections`` and step
through the active scan. The calculated CTR or Bragg positions should follow the
features in the detector images. If they do not, check the assigned ``hkl``
values, the selected image number of each reference reflection, the detector
calibration, and the beamline scan-axis convention.

The underlying calculation follows the Busing and Levy orientation-matrix
method adapted to the orGUI diffractometer geometry. The observed detector
coordinates and scan angles define momentum-transfer vectors in the inner
sample-circle frame; ``U`` is then obtained by matching those vectors to the
corresponding reciprocal lattice vectors.

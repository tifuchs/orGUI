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

Reference-Reflection Mismatch Coloring
--------------------------------------

The reference-reflection table reports the mean mismatch between the measured
momentum-transfer vectors and the vectors calculated from the current ``UB``
matrix. The angular mismatch :math:`\Delta\theta_i` is the angle between the
two vectors for reflection :math:`i`. The Q-norm mismatch
:math:`\Delta Q_i` is the absolute difference between their magnitudes, in
:math:`\mathrm{\AA}^{-1}`.

Rows are colored in two steps. First, orGUI computes a relative mismatch
ranking across the currently selected reference reflections. The angular
mismatch :math:`\Delta\theta_i` and the relative Q-norm mismatch
:math:`\Delta Q_i / \lVert Q_{\mathrm{UB},i} \rVert` are normalized
independently over the current table:

.. math::

   r_{\theta,i}
      =
      \frac{
         \Delta\theta_i - \min_j(\Delta\theta_j)
      }{
         \max_j(\Delta\theta_j) - \min_j(\Delta\theta_j)
      },

.. math::

   r_{Q,i}
      =
      \frac{
         \Delta Q_i / \lVert Q_{\mathrm{UB},i} \rVert
         -
         \min_j\left(
            \Delta Q_j / \lVert Q_{\mathrm{UB},j} \rVert
         \right)
      }{
         \max_j\left(
            \Delta Q_j / \lVert Q_{\mathrm{UB},j} \rVert
         \right)
         -
         \min_j\left(
            \Delta Q_j / \lVert Q_{\mathrm{UB},j} \rVert
         \right)
      }.

If all finite values in one channel are equal, that normalized channel is set
to zero for those finite rows. The relative color score is the mean of the two
normalized channels:

.. math::

   r_i = \frac{1}{2}\left(r_{\theta,i} + r_{Q,i}\right).

This score maps the best current agreement to green and the worst current
agreement to red. It is a comparative diagnostic within the current reference
set, so it remains useful when the fit is not yet limited by instrumental
resolution.

Second, orGUI checks whether an individual reflection is already within the
local detector-pixel resolution. This is an absolute test, not a relative
ranking. At the reflection pixel position :math:`(x_i, y_i)`, orGUI evaluates
the detector angle transform at the center pixel and at the two neighboring
pixels:

.. math::

   \boldsymbol{\theta}_i
      &= \boldsymbol{\theta}(x_i, y_i), \\
   \boldsymbol{\theta}_{x,i}
      &= \boldsymbol{\theta}(x_i + 1, y_i), \\
   \boldsymbol{\theta}_{y,i}
      &= \boldsymbol{\theta}(x_i, y_i + 1),

where :math:`\boldsymbol{\theta} = (\delta, \gamma)` in radians. The local
angular size of one detector pixel is then

.. math::

   \theta_{\mathrm{pix},i}
      =
      \max\left(
         \left\lVert
            \boldsymbol{\theta}_{x,i} - \boldsymbol{\theta}_i
         \right\rVert,
         \left\lVert
            \boldsymbol{\theta}_{y,i} - \boldsymbol{\theta}_i
         \right\rVert
      \right).

Using the incident wavevector magnitude :math:`K`, this local angular
tolerance is converted to a Q-norm tolerance:

.. math::

   Q_{\mathrm{pix},i} = K\,\theta_{\mathrm{pix},i}.

The pixel-equivalent mismatch score is

.. math::

   s_i =
   \max\left(
      \frac{\Delta\theta_i}{\theta_{\mathrm{pix},i}},
      \frac{\Delta Q_i}{Q_{\mathrm{pix},i}}
   \right).

Thus :math:`s_i \le 1` means that the reflection mismatch is within roughly one
local detector pixel of the current instrumental angular resolution. The
``Resolution limit`` control in the reference-reflection panel sets the maximum
allowed value of :math:`s_i` for this absolute test. Reflections below that
limit are marked blue to indicate that the agreement is resolution-limited.
All other reflections keep the relative green-to-red coloring from
:math:`r_i`.

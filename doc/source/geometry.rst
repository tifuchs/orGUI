orGUI Geometry
==============

orGUI was developed for high-energy surface X-ray diffraction (HESXRD) and
transmission surface diffraction workflows where a large, stationary 2D detector
records many reciprocal-space features during a sample rotation. The geometry
combines a pyFAI-style detector description with a surface-diffraction
diffractometer convention based on Busing and Levy and on Lohmeier and Vlieg.

.. image:: ../../orgui/resources/icons/diffractometer_v3.png
   :height: 360px
   :alt: orGUI diffractometer geometry

Reference Frames
----------------

The diffractometer origin is the sample position. In the laboratory frame,
``y`` points along the incident beam, ``x`` points perpendicular to the beam in
the radial direction, and ``z`` completes the right-handed frame. The azimuth
parameter in the machine settings rotates the diffractometer frame around the
incident beam direction; for an azimuth of 90 degrees, ``z`` points upward.

The sample is positioned by the ``alpha``, ``omega``, ``chi``, and ``phi``
circles. orGUI follows the Lohmeier and Vlieg angle convention for the shared
axes, but uses a ``phi`` circle around the ``x`` axis. The theta scan axis used
at some beamlines is mapped to the internal omega angle as
:math:`\theta = -\omega`.

Detector Geometry
-----------------

The detector model follows the pyFAI geometry convention. A detector pixel is
described relative to the point of normal incidence, the sample-detector
distance, and the detector rotations. The detector calibration therefore
provides the mapping between detector pixel coordinates and the scattering
angles observed by the area detector.

For surface diffraction, orGUI converts the detector scattering direction into
the in-plane and out-of-plane detector angles ``delta`` and ``gamma``. The
conversion accounts for the azimuthal reference direction of the surface normal
and for the current ``alpha`` setting. These angles are then used together with
the sample-circle angles in the reciprocal-space calculation.

Diffraction Equation
--------------------

The reciprocal lattice vector is written in reciprocal lattice coordinates as
:math:`\vec{H} = (h, k, l)^T`. The crystal ``B`` matrix converts this vector
from reciprocal lattice units to a Cartesian reciprocal-space vector, and the
orientation matrix ``U`` attaches the crystal to the inner sample-circle frame.

In orGUI notation, the sample rotations transform the crystal vector into the
``alpha`` frame. The momentum transfer vector is calculated from the incoming
and outgoing wavevectors in the same frame. The diffraction condition is met
when the transformed reciprocal lattice vector equals the measured momentum
transfer vector.

For crystal truncation rod and arbitrary line integration, the reciprocal-space
trajectory is described as

.. math::

   \vec{H}(s) = \vec{H}_0 + s\vec{H}_1

where :math:`\vec{H}_0` and :math:`\vec{H}_1` are vectors in reciprocal lattice
units. orGUI calculates the intersections of this line with the Ewald sphere and
converts the corresponding detector angles back to detector pixel coordinates.

Momentum Transfer and the Q-Plot
--------------------------------

The central plot of orGUI normally shows the detector image in **pixel
coordinates**; the reciprocal-space position of a pixel is reported in
reciprocal lattice units as the cursor is moved. The ``Q-plot`` toolbar action
switches the display to **momentum transfer**, the representation commonly used
for grazing-incidence work. This is a change of the displayed coordinates only,
it does not alter the loaded data.

Momentum transfer is split into a component parallel to the surface and one
along the surface normal,

* :math:`q_\parallel` -- in-plane momentum transfer, labelled ``qip`` and
  plotted as ``q_par``,
* :math:`q_\perp` -- out-of-plane momentum transfer, labelled ``qoop`` and
  plotted as ``q_perp``.

orGUI can calculate these two quantities in two different ways. Both describe
the same reciprocal space, but they differ in what they return and in when they
are used.

Per pixel with ``QAlpha``
~~~~~~~~~~~~~~~~~~~~~~~~~

This is the route used everywhere else in orGUI, for example for reflection
positions and for the rod integration. The detector calibration converts pixel
coordinates into the surface angles ``gamma`` and ``delta``, and the Vlieg
diffraction equation converts those angles into a Cartesian momentum transfer
vector in the ``alpha`` frame,

.. math::

   \begin{aligned}
   Q_x &= K \sin(\delta)\cos(\gamma) \\
   Q_y &= K \left[\cos(\delta)\cos(\gamma) - \cos(\alpha)\right] \\
   Q_z &= K \left[\sin(\gamma) + \sin(\alpha)\right]
   \end{aligned}

with :math:`K = 2\pi/\lambda`. The out-of-plane component is :math:`Q_z`, and
the in-plane component is the radial component in the surface plane:

.. code-block:: python

   import numpy as np

   detectorCal = ubcalc.detectorCal   # DetectorCalibration.Detector2D_SXRD
   angles = ubcalc.angles             # HKLVlieg.VliegAngles
   alpha_i = ubcalc.mu                # incidence angle in rad

   gamma, delta = detectorCal.surfaceAngles(alpha_i, shape=image.shape)
   Qxyz = angles.QAlpha(alpha_i, delta, gamma)   # inverse Angstrom

   q_perp = Qxyz[..., 2]                          # out-of-plane
   q_par = np.hypot(Qxyz[..., 0], Qxyz[..., 1])   # in-plane, unsigned

The result has the same shape as the detector image: every pixel keeps its own
momentum transfer. The sampling is therefore irregular and curved in
:math:`(q_\parallel, q_\perp)`, which is exactly what is needed to look up
where a reflection falls, but it cannot be handed to an image widget directly.

On a regular grid with ``FiberIntegrator``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   The ``Q-plot`` and the :mod:`orgui.app.qconversion` module are
   **experimental**. They live in the application layer, not in
   ``orgui.datautils.xrayutils``, precisely so that they are not mistaken for
   the production reciprocal-space code. The conventions described here may
   still change.

This is the route used by the ``Q-plot`` action. pyFAI's ``FiberIntegrator``
computes the same per-pixel momentum transfer and additionally **rebins** the
intensities onto a regular grid, which is what makes the result displayable as
an image. orGUI drives it through :mod:`orgui.app.qconversion`, which supplies
pyFAI with orGUI's own angle conventions:

.. code-block:: python

   from orgui.app import qconversion

   res = qconversion.integrateImage(
       detectorCal,
       image,
       alpha_i,                 # incidence angle in rad
       frame="Q_alpha",         # see "Reciprocal-space frames" below
   )

   res.intensity     # rebinned image, shape (npt_oop, npt_ip)
   res.inplane       # q_par axis
   res.outofplane    # q_perp axis

``FiberIntegrator`` was introduced in pyFAI 2025.1. With an older pyFAI the
``Q-plot`` action is unavailable and reports an error instead.

Which one to use
~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * -
     - ``QAlpha``
     - ``FiberIntegrator``
   * - Returns
     - one :math:`(q_\parallel, q_\perp)` per pixel
     - intensities rebinned onto a regular grid
   * - Sampling
     - irregular, curved
     - regular, uniform
   * - In-plane sign
     - unsigned, :math:`q_\parallel \ge 0`
     - signed
   * - Requires
     - nothing beyond orGUI
     - pyFAI >= 2025.1
   * - Use for
     - reflection positions, integration, any calculation
     - displaying a whole image in reciprocal space

Use ``QAlpha`` whenever the momentum transfer of a specific pixel or reflection
is needed. Use ``FiberIntegrator`` when a complete image has to be shown or
exported in reciprocal-space coordinates, and accept that the rebinning
quantises positions to the width of one grid cell.

Reciprocal-space frames
~~~~~~~~~~~~~~~~~~~~~~~

``QAlpha`` returns the momentum transfer in the **alpha frame**, the frame of
the sample surface. The drop-down next to the ``Q-plot`` action selects which
frame the image is displayed in; the remaining frames are reached by undoing
the sample rotations, following
:math:`\vec{H} = UB^{-1}\,\Phi^{-1}\,X^{-1}\,\Omega^{-1}\,\vec{Q}_\alpha`.

.. list-table::
   :header-rows: 1
   :widths: 20 46 34

   * - Frame
     - Reached from the alpha frame by
     - Defined for
   * - ``Q_alpha``
     - nothing, this is the default
     - any image
   * - ``Q_lab``
     - the ``alpha`` rotation
     - any image
   * - ``Q_omega``
     - undoing ``omega``
     - a single image only
   * - ``Q_chi``
     - undoing ``omega`` and ``chi``
     - a single image only
   * - ``Q_phi``
     - undoing ``omega``, ``chi`` and ``phi``
     - a single image only
   * - ``Q_cryst``
     - additionally undoing the orientation matrix ``U``
     - a single image only

For every frame the out-of-plane component is the ``z`` axis of that frame and
the in-plane component is the radial component in its ``xy`` plane. Maximum and
sum images combine many ``omega`` angles, so the frames that undo ``omega`` are
refused for them. ``Q_cryst`` additionally needs the orientation matrix.

``Q_cryst`` is the Cartesian reciprocal-space frame of the crystal, so it holds
:math:`B\vec{H}`. Multiplying ``Q_phi`` by :math:`UB^{-1}`, or equivalently
``Q_cryst`` by :math:`B^{-1}`, gives the reciprocal lattice coordinates
:math:`\vec{H} = (h, k, l)^T`, the same result as
:meth:`~orgui.datautils.xrayutils.HKLVlieg.VliegAngles.anglesToHkl`.

Why orGUI supplies its own pyFAI units
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pyFAI's built-in fiber units cannot express orGUI's geometry, for two
independent reasons.

**The azimuth is continuous, ``sample_orientation`` is not.** pyFAI's
``sample_orientation`` is the EXIF dihedral group acting on the image array, so
it only spans quarter turns. It describes how the image array is laid out, not
where the sample sits. orGUI's azimuthal reference rotates the ``alpha``
rotation axis about the incident beam and takes arbitrary values, so it is
passed through ``tilt_angle`` instead, while ``sample_orientation`` stays at the
identity:

.. math::

   \mathrm{tilt\_angle} = -\left(\mathrm{azimuth} + \frac{\pi}{2}\right)

The quarter turn is the offset between orGUI's azimuthal reference and pyFAI's
fiber convention. Encoding the azimuth in ``sample_orientation`` instead would
only be correct for azimuths that are exact multiples of 90 degrees.

**The rotations are composed in the opposite order.** pyFAI combines the
incidence and tilt rotations extrinsically, about fixed axes, as
:math:`R_x(-\mathrm{tilt})\,R_y(\mathrm{incidence})`. In orGUI the azimuth
rotates the ``alpha`` axis itself, so the incidence rotation has to be applied
about the *already rotated* axis,
:math:`R_y(\mathrm{incidence})\,R_x(-\mathrm{tilt})`. The two agree only when
one of the two angles vanishes.

:mod:`~orgui.app.qconversion` therefore builds its own
``pyFAI.units.UnitFiber`` objects that implement orGUI's convention, and hands
them to ``FiberIntegrator``. pyFAI still performs the rebinning; only the
coordinate definition is replaced.

With that in place the two routes agree to numerical precision. This is
verified in ``orgui/app/test/test_q_conversion.py``, which compares them per
pixel for flat and tilted detectors, at several incidence angles and at
azimuths that are deliberately not multiples of 90 degrees. The same file
checks that ``Q_phi`` reproduces ``anglesToHkl``, that ``Q_cryst`` reproduces
:math:`B\vec{H}`, that every frame preserves :math:`|\vec{Q}|`, and that a
single bright pixel is placed by the rebinning within one grid cell of the
position predicted by ``QAlpha``.

Performance
~~~~~~~~~~~

Detector images are large, so the conversion is collapsed into a single affine
relation before any pixel data is touched. With
:math:`n = |(x, y, z)|` every component reduces to

.. math::

   q_j = k \left( \frac{G_{j0} x + G_{j1} y + G_{j2} z}{n} - c_j \right)

where the matrix :math:`G` and the offset :math:`\vec{c}` absorb the sample
orientation map, the beam and incidence rotations, the axis relabelling and the
frame rotation. A compiled kernel,
``orgui/app/cpp/qconversion_cpp.cpp``, evaluates both displayed quantities in a
single pass; a plain numpy implementation is used when the extension has not
been built. The result is cached, because pyFAI evaluates the in-plane and the
out-of-plane unit separately but passes the same pixel positions to both, so an
image is converted once rather than twice.

On a 6.2 megapixel Pilatus 6M the compiled kernel converts a full image in
roughly 50 ms, and the second unit evaluation is served from the cache.

.. note::

   ``numexpr`` is deliberately not used here. Version 2.11.0 returns wrong
   results for a small fraction of multi-threaded evaluations of expressions of
   this kind, which is not acceptable for a coordinate transform.

Further Reading
---------------

The implemented geometry is summarized in chapter 5, especially section 5.3, of
Timo Fuchs' dissertation. The angle convention builds on:

* W. R. Busing and H. A. Levy, Acta Crystallographica 22, 457-464 (1967).
* M. Lohmeier and E. Vlieg, Journal of Applied Crystallography 26, 706-716
  (1993), https://doi.org/10.1107/S0021889893004868.

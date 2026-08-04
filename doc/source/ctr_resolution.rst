CTR Resolution Modeling
=======================

``CTRresolution`` supplies optional instrumental-resolution broadening for
calculated crystal truncation rods. It does not change ``SXRDCrystal.F``.
``CTRopt`` can opt into the same sampled broadening for calculated fit
amplitudes.

Intensity and output convention
-------------------------------

Resolution acts on calculated intensity rather than on the complex amplitude:

.. math::

   I_{\mathrm{res}}(L_i) =
   \frac{\int R_i(L-L_i)\lvert F(H_i,K_i,L)\rvert^2\,\mathrm{d}L}
        {\int R_i(L-L_i)\,\mathrm{d}L}.

The returned ``CTRCollection`` follows the existing ``sfI`` convention and
stores the effective amplitude

.. math::

   F_{\mathrm{eff}}(L_i) = \sqrt{I_{\mathrm{res}}(L_i)}.

Phase and experimental-error information cannot be retained by this operation
and is therefore absent from the result. Input collections are never modified.

Resolution functions
--------------------

The built-in functions are ``BoxResolution`` and ``GaussianResolution``. Their
effective width is in reciprocal lattice units:

.. math::

   \Delta L = \Delta L_0 + \Delta L_1\lvert\sin(\gamma)\rvert
   + \Delta L_2\lvert\sin(\delta)\rvert.

Here ``gamma`` and ``delta`` are Vlieg detector angles in radians. For a box,
``DeltaL`` is the complete support width. For a Gaussian it is the full width
at half maximum (FWHM). Setting both angle-dependent contributions to zero
creates a constant function; setting all contributions to zero disables
broadening.

Angle-dependent models require angle records on their input collection. They
can be calculated for all rods at once:

.. code-block:: python

   ub = HKLVlieg.UBCalculator(crystal.reference_uc, energy=70.0)
   ub.defaultU_GID()
   angle_calculator = HKLVlieg.VliegAngles(ub)
   ctrs.calcAnglesZmode(
       angle_calculator,
       fixedangle=np.deg2rad(0.1),
       fixed="in",
   )

Fast convolution on existing points
-----------------------------------

``fast_convolve`` uses only the L coordinates and amplitudes already present
in each rod. It does not interpolate and does not require an equidistant grid.
Local composite-trapezoidal weights account for unequal point spacing, so a
dense part of a scan is not overrepresented merely because it has more points.

At the first and last points the kernel is truncated to the available data and
renormalized. This preserves constant intensity without assuming either zero
intensity or an extrapolation outside the measured range. L coordinates must
be finite and unique within a rod, but they may be in any order.

.. code-block:: python

   from orgui.datautils.xrayutils import CTRresolution

   model = CTRresolution.GaussianResolution(
       delta_l_0=0.015,
       delta_l_1=0.08,
   )
   broadened = CTRresolution.fast_convolve(calculated_ctrs, model)

The fast method approximates the convolution represented by the existing
sampling. If the kernel is narrower than the local point spacing, it may only
contain the central point and consequently cannot reproduce unresolved detail.

Sampled structure factors
-------------------------

``sample_structure_factor`` reevaluates ``crystal.F`` around every requested
point and is slower but does not depend on the density of the input L samples.
It uses Gauss-Legendre quadrature for a box and Gauss-Hermite quadrature for an
unbounded Gaussian. The quadrature order must be a positive odd integer and
defaults to 25.

.. code-block:: python

   sampled = CTRresolution.sample_structure_factor(
       data_positions,
       crystal,
       model,
       quadrature_order=25,
   )

The central point's gamma determines the width of all quadrature samples for
that point; the same applies to delta when ``delta_l_2`` is nonzero. H and K
remain fixed. This implementation therefore models only out-of-plane L
resolution.

Fitting resolution widths
-------------------------

Use :meth:`CTRopt.CTROptimizer.fit_resolution` to include the three resolution
widths directly in an optimizer's parameter vector. They are always the first
three parameters, in ``delta_l_0``, ``delta_l_1``, ``delta_l_2`` order, and use
r.l.u. bounds. Angle-dependent widths require cached Z-mode angle records;
calculate them once because they depend on the fixed CTR geometry and
experimental setup, not on the fitted widths.

.. code-block:: python

   from orgui.datautils.xrayutils import CTRopt, CTRresolution

   optimizer = CTRopt.CTROptAngleCorrection(crystal, measured_ctrs)
   optimizer.calc_resolution_angles_zmode(angle_calculator, fixedangle=0.1)
   optimizer.fit_resolution(
       CTRresolution.GaussianResolution(0.015, 0.08, 0.02),
       lower_bounds=(0.0, 0.0, 0.0),
       higher_bounds=(0.1, 0.2, 0.1),
   )
   optimizer.prepareFit()
   parameters = optimizer.get_parameters()

For ``CTROptAngleCorrection``, fit traces label these leading values as
``resolution_delta_l_0``, ``resolution_delta_l_1``, and
``resolution_delta_l_2``. Their estimated standard errors are stored in
``optimizer.resolution_errors`` after statistics are evaluated.

``fit_resolution(..., calculation="sample")`` is the default and evaluates
the crystal at quadrature points. Use ``calculation="convolve"`` to convolve
calculated values on the existing L grid instead; this is typically faster for
fits but inherits the sampling limitations of :func:`fast_convolve`. In both
modes, the optimizer retains a separate ``calculated_CTRs`` collection whose
``sfI`` arrays are updated once per parameter assignment; experimental CTRs
are not modified.

A complete runnable comparison using the bundled CTR reference data is in
``examples/CTR/ctr_resolution_example.ipynb``.

API reference
-------------

.. autoclass:: orgui.datautils.xrayutils.CTRresolution.ResolutionFunction
   :members: width, weights, quadrature
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRresolution.BoxResolution
   :members:
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRresolution.GaussianResolution
   :members:
   :member-order: bysource

.. autofunction:: orgui.datautils.xrayutils.CTRresolution.fast_convolve

.. autofunction:: orgui.datautils.xrayutils.CTRresolution.sample_structure_factor

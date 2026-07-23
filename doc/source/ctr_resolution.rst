CTR Resolution Modeling
=======================

``CTRresolution`` supplies optional instrumental-resolution broadening for
calculated crystal truncation rods. It does not change ``SXRDCrystal.F`` or
enable broadening implicitly in fitting code.

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

   \Delta L = \Delta L_0 + \Delta L_1\lvert\sin(\gamma)\rvert.

Here ``gamma`` is the Vlieg out-of-plane detector angle in radians. For a box,
``DeltaL`` is the complete support width. For a Gaussian it is the full width
at half maximum (FWHM). Setting ``delta_l_1=0`` creates a constant function;
setting both contributions to zero disables broadening.

Gamma-dependent models require angle records on their input collection. They
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
that point. H and K remain fixed. This first implementation therefore models
only out-of-plane L resolution; a future in-plane model can use the H, K, L,
and complete angle context accepted by ``ResolutionFunction.width``.

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

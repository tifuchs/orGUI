Distorted-Wave Born Approximation (DWBA)
========================================

This chapter documents the atomistic DWBA calculation for CTRs. Structure
factors and matrix elements use electrons per reference lateral cell, following
the normalization conventions in :doc:`ctr_structure_factors`.

Unit-cell distorted-wave Born amplitude
----------------------------------------

The DWBA interface is a lazy runtime state owned by each crystal.  It is not
written to ``.xtal`` files.  Configure a non-specular scan once, then evaluate
any compatible rod in reference-cell r.l.u.:

.. code-block:: python

   crystal.dwba.set_ctr_geometry(alpha_i=np.deg2rad(0.2))
   result_20L = crystal.dwba.evaluate(2.0, 0.0, L)
   result_00L = crystal.dwba.evaluate(0.0, 0.0, L_specular)

Exactly one of ``alpha_i``, ``alpha_f``, or ``equal_angles=True`` configures a
rule.  Angles are radians.  ``rods=None`` establishes the default
non-specular rule; ``rods=[(h, k), ...]`` establishes overrides.  The
``(0, 0, L)`` rod always uses equal angles.  ``set_orientation(U)`` replaces
the default identity orientation used by the internal ``UBCalculator`` and
``VliegAngles`` objects.

``prepare`` returns an immutable ``PreparedCTR`` that freezes geometry,
polarization, the optical reference :math:`n_0`, generated-record topology,
and incident and reciprocal exit fields.  ``evaluate`` caches preparations
automatically, while ``evaluate_prepared`` deliberately continues using the
retained preparation's frozen :math:`n_0` as atomic coordinates, species,
occupancies, displacement factors, component weights, and domain occupancies
change.  Energy, lattice, reference transform, orientation, reference area,
record set, domain matrix, or stacking-geometry changes invalidate retained
geometry.  ``cache_info`` reports the split geometry, optical-reference,
unique-angle field, and 32-entry preparation caches; ``clear_cache`` clears
them.  Repeated glancing angles are solved and stored once, with native index
maps selecting their field columns for each scan point.

Measured geometries are available without changing a persistent rod rule.
``prepare_from_glancing(h, k, alpha_i, alpha_f)`` derives ``l`` from the Ewald
condition, and ``prepare_from_vlieg(alpha, delta, gamma, omega, chi, phi)``
derives hkl using the configured orientation.  Their ``evaluate_...``
counterparts prepare and evaluate in one call.  Inputs broadcast; scalar input
produces scalar amplitude properties.

The measured Vlieg entry points currently require ``chi=phi=0`` (within
``1e-12 rad`` of zero). These inner sample circles rotate the surface normal,
so nonzero values cannot be interpreted by simply using ``alpha`` and
``gamma`` as glancing angles in the configured DWBA surface frame. They raise
``ValueError``; general six-circle frame conversion is not implemented.
Kinematic ``CTR.calcAnglesZmode`` calculations retain support for nonzero
sample circles.

When generating angle records with a calculator configured for the bulk cell,
pass the reference-to-bulk mapping explicitly::

   angles = ctr.calcAnglesZmode(
       bulk_vlieg_angles,
       fixedangle=np.deg2rad(0.2),
       hkl_transform=crystal.uc_bulk.refHKLTransform,
   )
   prepared = crystal.dwba.prepare_from_vlieg(
       *(angles[name] for name in ("alpha", "delta", "gamma", "omega", "chi", "phi"))
   )

The calculator and DWBA state must use the same bulk lattice, energy, and
orientation. ``hkl_transform`` maps reference HKL columns into calculator HKL
columns, both in r.l.u.; omitting it preserves the existing convention that
the calculator already uses the CTR reference lattice. ``CTRCollection``
forwards the same transform to each rod. The DWBA measured-angle path applies
the inverse reference transform when recovering HKL.

The optical solver uses the Renaud-signed normal wavevector

.. math::

   k_{z,j}=-k_0\sqrt{n_j^2-n_0^2\cos^2\alpha},
   \qquad \Re k_{z,j}\leq0,\quad\Im k_{z,j}\geq0,

and stores downward and upward field amplitudes as ``A_plus`` and
``A_minus``.  For channel signs
:math:`\sigma_i,\sigma_f\in\{+1,-1\}` the native kernel uses

.. math::

   Q_{z,j}^{\sigma_i\sigma_f}
   =-\left(\sigma_i k_{zj,i}+\sigma_f k_{zj,f}\right).

The reciprocal exit field is solved through the same complex profile without
conjugating its amplitudes.  Explicit ``"s"`` and ``"p"`` fields use the
bilinear reciprocal-field contraction.  Each channel receives its own
generally complex :math:`Q^2`, analytic Waasmaier--Kirfel factor, anisotropic
Debye--Waller factor, and reference-position phase.  Every transformed atom is
assigned to the optical medium containing its physical surface-normal
coordinate; an atom exactly on an interface belongs to the upper medium.
The ``++``, ``+-``, ``-+``, and ``--`` amplitudes are summed coherently.  Bulk
repeats crossing finite graded media are evaluated explicitly.  Once a full
repeat lies in the terminal substrate, its closed-form lattice tail contains
only ``++`` because :math:`A_i^-=A_f^-=0` there.

Atomic ``UnitCell``, ``Film``, ``EpitaxyInterface``, and ``PoissonSurface``
components participate in this calculation.  Their stacking-generated cells,
signed domain occupancies, crystal weights and domains, and
:math:`A_{\mathrm{ref}}/A_{\mathrm{uc}}` lateral-area conversion are retained.
Continuous ``WaterModel`` components, Python special form-factor callbacks,
surface-normal tilts, changed projected areas, and in-plane strain transforms
are rejected explicitly.

At a specular point, the kernel subtracts the Fourier amplitude of each
record's share of the same piecewise-constant :math:`\delta+i\beta` reference
used by the electric-field solver.  Finite layers use a stable finite-slab
integral and the substrate uses the terminal half-space integral.  This planar
reference amplitude is exactly zero for non-specular points; no near-specular
interpolation is applied.

The optional empirical ``bulk_attenuation`` exponent defaults to zero because
the complex DWBA wavevectors already contain physical absorption.  A nonzero
value remains available for kinematical comparisons and is applied
continuously within the top cell and between repeats.  At zero attenuation an
exact Bragg pole is rejected.  ``CTRutil.set_atten_from_dwba`` estimates the
corresponding scalar kinematical attenuation at a fixed incidence angle;
``CTRutil.attenuation_from_dwba`` supplies its broadcast diagnostic form.

Returned quantities and the prefactor
-------------------------------------

``DWBAResult`` returns one matrix element; every other amplitude is derived
from it by a single prefactor.  Writing :math:`\mathbf{h}` for the rod, the matrix
element :math:`F_{\mathbf{h}}` is converted into a dimensionless reflection
coefficient by

.. math::

   r_{\mathbf{h}}=
   \frac{2\pi i r_e}{A_{\mathrm{ref}}\,\kappa_f}\,F_{\mathbf{h}},
   \qquad
   \kappa_f=k_0\sin\alpha_f,

with :math:`r_e` the classical electron radius and :math:`\kappa_f` the
**exit** normal wavevector in vacuum, which off specular differs from the
incident one.  :math:`\kappa_f` is positive, whereas the design note
``doc/design/DWBA/dwba_bulk_specular_math.tex`` writes the same relation as its
Eq. (1.33) with a minus sign because it uses the negative-root convention
:math:`k_{z,0,f}=-\kappa_f`.  The two are identical.

The classical electron radius therefore enters **once**, inside this
prefactor.  The separate convention, in which a bare matrix element becomes a
cross-section kernel through :math:`r_e^2|F_{\mathbf{h}}|^2`, applies to
``F_h`` and not to any reflection coefficient; the two must never be mixed.

Only ``F_h``, ``unperturbed_amplitude``, and the per-record arrays inside
``contributions`` are stored.  The native kernel computes exactly two arrays
per generated record — the atomic and the planar-reference amplitude — and
those, together with the Fresnel :math:`r_0` the field solver already
produced, are the only independent quantities.  Everything else below is
derived arithmetic evaluated on access, which costs microseconds and avoids
retaining redundant arrays during a fit.  The chain ``F_h`` →
``scattered_amplitude`` → ``total_amplitude`` → ``F_effective`` means these
members are not independent.

``F_h``
   *Stored.*  The coherent DWBA contrast matrix element
   :math:`F_{\mathbf{h}}`.  At the specular rod it already contains the
   constant :math:`\delta_m+i\beta_m` slab terms; on every genuine CTR it
   cannot, because the laterally uniform optical reference has no Fourier
   component there.
``F_atomic`` and ``F_reference``
   *Derived.*  The coherent sums of the actual atomic amplitude and the planar
   optical-reference amplitude.  ``F_reference`` is **positive** and
   enters by subtraction, ``F_h == F_atomic - F_reference``; the design note
   writes the same physics as one Fourier coefficient whose contrast
   :math:`\Delta\rho_f` already carries :math:`-\overline{\rho}_{f,m}`.  It is
   identically zero on every non-specular point, where all records share one
   zero array rather than each holding its own.
``contributions``
   *Stored.*  An ordered tuple of immutable ``DWBAContribution`` records: bulk
   first, followed by components in crystal order and their generated cells in
   model order.  Each record identifies its component, generated-record index
   and role, optional structural layer, and exposes read-only ``F_atomic`` and
   ``F_reference`` amplitudes plus the derived ``F_h``.  Roles are ``bulk``,
   ``unit_cell_layer``, ``film_layer``, ``interface_top``,
   ``interface_bottom``, ``surface_termination``, ``covered_film``, and
   ``sharp_film_correction``.
``unperturbed_amplitude``
   *Stored.*  The Fresnel reference amplitude :math:`r_0` for specular points and zero
   for non-specular points.
``scattered_amplitude``
   *Derived.*  The dimensionless scattered amplitude :math:`r_{\mathbf{h}}`, i.e. the prefactor
   above applied to ``F_h``.
``total_amplitude``
   *Derived.*  :math:`r=r_0+r_{\mathbf{h}}`, that is
   ``unperturbed_amplitude + scattered_amplitude``.
``F_effective``
   *Derived.*  ``total_amplitude`` with the prefactor inverted,

   .. math::

      F_{\mathrm{eff}}=
      \frac{A_{\mathrm{ref}}\,\kappa_f}{2\pi i r_e}\,\left(r_0+r_{\mathbf{h}}\right),

   i.e. the structure factor a kinematic analysis would infer from this
   reflectivity.  Use it to compare a DWBA model against kinematically reduced
   data.  Note it is a *DWBA* output, not the kinematical model's own
   structure factor; that is ``SXRDCrystal.F``, which the diagnostics notebook
   writes :math:`F_{\mathrm{kin}}`.  The two are compared with each other, not
   interchangeable.

   It is **not** linear in the model electron density, because
   :math:`r_0` is zeroth order and exact in the reference to all orders, so it
   is not the DWBA structure factor of the sample either; that is ``F_h``.  It
   vanishes as :math:`\alpha_f\to0` and is meaningful only well above the
   critical angle, and off specular it degenerates to ``F_h``.  Its
   justification is the large-:math:`Q_z` limit, where
   :math:`r_0\to4\pi r_e\overline{\rho}_f/Q_z^2` and it tends to
   :math:`F_{\mathbf{h}}+A_{\mathrm{ref}}\overline{\rho}_f/(iQ_z)`,
   the full kinematic structure factor of the actual density.
``reflectivity``
   *Derived.*  The coherent specular reflectivity :math:`|r_0+r_{\mathbf{h}}|^2`.  This is the
   primary intensity observable.  On the specular rod at small angles it is
   *exact* up to the first-order truncation in the contrast: the problem is
   genuinely one-dimensional there, so the scalar reduction involves no
   approximation and the ``s`` and ``p`` channels become degenerate, while
   :math:`r_0` carries the exact Fresnel and refraction response of the
   prepared reference to all orders.  The residual error is therefore second
   order in the density contrast and shrinks as the angle decreases; below the
   critical angle the exact reflectivity is unity and this expression
   reproduces it.  A reference whose density is wrong by ten percent gives a
   relative error near :math:`10^{-4}` at one tenth of the critical angle,
   rising to a few percent near :math:`\alpha_c`, and falling quadratically as
   the contrast error is reduced.  Squaring the amplitude rather than expanding
   it also keeps the result non-negative.  ``reflectivity`` requires an
   entirely specular, same-polarization, semi-infinite result.

Quantities that are one expression away are deliberately not provided, so that
there is exactly one way to spell each of them:

.. code-block:: python

   abs(result.F_h) ** 2                      # squared structure factor
   abs(result.scattered_amplitude) ** 2      # squared scattered amplitude
   r_e**2 * abs(result.F_h) ** 2             # cross-section kernel per cell,
                                             # not integrated over acceptance,
                                             # footprint, coherence, resolution,
                                             # flux, exposure, or efficiency

   # strictly linearised reflectivity, a regime diagnostic: a large difference
   # from result.reflectivity means |r_h| is not small against |r_0|, so the
   # first-order treatment is being pushed.  It can go negative for that reason.
   (abs(result.unperturbed_amplitude) ** 2
    + 2.0 * np.real(np.conj(result.unperturbed_amplitude)
                    * result.scattered_amplitude))

``crystal.dwba.reflectivity(..., polarization="unpolarized")`` evaluates the
``s`` and ``p`` channels independently and averages their reflectivities
incoherently.  The default ``bulk_mode="semi_infinite"`` performs the physical
atom-minus-reference calculation and includes the unperturbed specular
reflection.  Diagnostic ``bulk_mode="unit_cell"`` instead evaluates exactly
one bulk repeat plus every finite generated record, sets every
``F_reference`` to zero, and omits the unperturbed reflection.  Experimental
acceptance integration remains a separate operation.

.. _dwba-fitting:

Fitting measured CTRs with DWBA
-------------------------------

``CTROptimizer`` can evaluate the live DWBA model instead of the default
kinematical structure factor. The switch is explicit and runtime-only:

.. code-block:: python

   from orgui.datautils.xrayutils import CTRopt

   fit = CTRopt.CTROptimizer(
       crystal,
       measured_ctrs,
       scale_policy={"F": "scaled", "R": "fixed"},
   )
   fit.set_dwba()
   fit.prepareFit()

   prediction = fit.flat_prediction()
   normalized_residual = fit.weighted_residues()

Configure the crystal's energy, lattice, reference transform, and DWBA
orientation before constructing the optimizer: the optimizer owns deep copies
of both the crystal and the measurements. Each measured CTR must carry a
:class:`~orgui.datautils.xrayutils.CTRplotutil.MeasurementReduction` with a
non-``None``
:class:`~orgui.datautils.xrayutils.CTRplotutil.PolarizationReduction`.
The incident ``s_fraction`` and outgoing analyser setting select only the
required one, two, or four polarization channels; their field intensities are
summed incoherently.

The prediction remains in each dataset's stored quantity. If
:math:`R_{\mathrm{res}}` is the direct or resolution-broadened field intensity,
:math:`c` is ``DWBAResult.amplitude_prefactor`` at the central geometry, and
:math:`P` is the conventional polarization factor stored for F data, then

.. math::

   \widehat R = R_{\mathrm{res}}, \qquad
   \widehat F = \frac{\sqrt{R_{\mathrm{res}}/P}}{|c|}.

``P=None`` means :math:`P=1`. Corrected reflectivity data cannot carry P:
stored R is compared directly with :math:`\widehat R`. The optimizer never
converts measured values or uncertainties between F and R. Signed and zero
observations are therefore retained, while uncertainties must be finite and
strictly positive. Analytical ``fixed``, ``scaled``, and ``global`` policies
are applied to the prediction only; see :doc:`ctr_structure_factors` for the
common fitting-result and grouping contract.

Geometry is measurement-owned. Without resolution, and for convolution on the
measured L grid, point-aligned six-circle records in radians are sufficient.
Alternatively attach a z-mode scan rule before constructing the optimizer:

.. code-block:: python

   from orgui.datautils.xrayutils.CTRplotutil import CTRScanGeometry

   for rod in measured_ctrs:
       rod.scan_geometry = CTRScanGeometry(
           fixed="in", angle=np.deg2rad(0.3), mirrorx=False
       )

Fixed-width quadrature sampling requires this rule because displaced L points
need new physical angles. If both records and a rule are present, the rule must
reproduce the measured central physical geometry and scattering branch.
Every record is checked against its H, K, and L coordinates; periodic rotations
such as an omega shift by :math:`2\pi` remain equivalent. The existing
``chi=phi=0`` limitation of the measured-angle DWBA path still applies.

Resolution always averages field intensity before the F/R representation is
formed. Measured-grid ``calculation="convolve"`` and fixed-width
``calculation="sample"`` are supported. Fitting widths with regenerated
sampling grids is not; use fitted-width convolution instead. Details and code
are in :doc:`ctr_resolution`.

Fitting uses physical ``bulk_mode="semi_infinite"``. A fixed nonnegative
empirical exponent can be supplied with
``fit.set_dwba(bulk_attenuation=value)`` for controlled kinematical
comparisons. ``get_dwba()`` returns both current settings. One outer
``DWBAState.batch()`` is used for every objective evaluation, sharing the live
atomic snapshot and atom packing across rods and requested polarization
channels. Persistent preparation and field reuse remains owned by
``DWBAState``.

``CTROptAngleCorrection`` deliberately rejects ``set_dwba`` because the
meaning of its empirical angle correction for mixed F/R DWBA predictions is
not defined. Continuous-density ``WaterModel`` components, nonzero chi/phi,
and the other core limitations listed above also remain unsupported.

Comparing a DWBA rod with a kinematical rod
--------------------------------------------

The :doc:`RuO2/TiO2 DWBA diagnostics tutorial
<dwba/dwba_diagnostics>` works through the refractive-index split,
specular Fresnel-plus-contrast amplitude, critical edge, extended specular
rod, and fixed-incidence non-specular rods for flat and rough film models.
Its executed plots compare each DWBA observable with the corresponding
kinematical curve.

``DWBAResult.F_h`` and ``SXRDCrystal.F`` are both amplitudes in
electrons per reference lateral cell, and off specular they are directly
comparable.  (Use ``F_effective`` instead when comparing on the specular rod,
where ``F_h`` alone omits the unperturbed Fresnel amplitude.)  They are not identical even in the weak-scattering limit, and
four separate terms account for the difference.  Three are physics and one is
a difference of convention, so a comparison that does not separate them can be
read the wrong way round.

Polarization
   The DWBA matrix element contains the bilinear polarization contraction; for
   ``"s"`` this is :math:`P_{ss}=\cos\varphi`, available as
   ``PreparedCTR.cos_azimuth``.  ``SXRDCrystal.F`` is a bare structure factor
   with no polarization factor, so the DWBA amplitude is smaller by
   :math:`\cos\varphi`.  For an in-plane scattering angle of
   :math:`12^\circ` this is a two percent effect.

Optical field amplitude
   Every record is weighted by :math:`A_i^{\sigma_i}A_f^{\sigma_f}` of the
   medium that contains it, not of the ambient.  A film that is denser than
   its substrate has its own, larger critical angle, and an incidence angle
   chosen as a safe multiple of the *substrate* critical angle can still sit on
   the flank of the *film* resonance.  Inspect ``field_i.n`` for the
   per-medium :math:`\delta` and ``field_i.A_plus`` and ``field_i.A_minus``
   for the branch amplitudes actually used.

Bulk truncation
   The kinematical bulk series is damped by the empirical scalar
   ``SXRDCrystal.atten``; the DWBA series is damped by
   :math:`\mathrm{Im}\,k_z`.  Off Bragg this changes the bulk term by a small
   percentage, which is invisible while the bulk term dominates the rod.  It
   is not invisible where the bulk term and a surface correction nearly
   cancel: a rough or partly dissolved surface has such near-cancellations
   between its Bragg poles, and there a one percent change of the bulk term
   can be the size of the whole remaining amplitude.  Match the two
   prescriptions with ``CTRutil.attenuation_from_dwba`` before comparing.

Refraction
   The internal normal wavevectors are smaller than the vacuum ones, so a DWBA
   feature sits at slightly larger nominal :math:`l`.  The shift is a fixed
   small number in r.l.u. and is negligible on a smooth part of a rod.  It
   dominates at any sharp feature, because :math:`|\mathrm{d}F/\mathrm{d}l|`
   is then large compared with :math:`|F|`.  The shift is medium dependent:
   a bulk Bragg pole follows the substrate :math:`\delta`, while the film-
   dominated rod between the poles follows the film :math:`\delta`.

``examples/CTR/RuO2_TiO2_nonspecular_DWBA.ipynb`` works this decomposition
through for the RuO2/TiO2 Poisson dissolution family on the ``(0, 1, L)``,
``(1, 1, L)``, and ``(2, 0, L)`` rods, and shows that once the four terms are
accounted for, the corrected kinematical rod converges onto the DWBA rod in the
Born limit.

DWBA diagnostics tutorial
-------------------------

The executed tutorial below works through the refractive-index decomposition,
specular reflectivity, and non-specular CTRs for flat and rough RuO2/TiO2 film
systems. Its rendered plots include the corresponding kinematical curves.

.. toctree::
   :maxdepth: 1

   dwba/dwba_diagnostics

Kinematical attenuation and valid L range
-----------------------------------------

The attenuation tutorial uses a RuO2/TiO2 model with both surface and
interfacial roughness set to zero. It answers how to set the empirical
kinematical `atten` exponent from the DWBA wavevectors and how to derive
a defensible minimum exit angle and minimum :math:`L` for a
fixed-incidence rod.

.. toctree::
   :maxdepth: 1

   dwba/kinematic_attenuation

API reference
-------------

.. autoclass:: orgui.datautils.xrayutils.CTRdwba.DWBAContribution
   :members:
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRdwba.DWBAResult
   :members:
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRdwba.DWBAState
   :members: set_orientation, set_ctr_geometry, prepare, prepare_from_glancing, prepare_from_vlieg, evaluate, evaluate_prepared, evaluate_from_glancing, evaluate_from_vlieg, reflectivity, cache_info, clear_cache
   :member-order: bysource

.. automethod:: orgui.datautils.xrayutils.CTRopt.CTROptimizer.set_dwba

.. automethod:: orgui.datautils.xrayutils.CTRopt.CTROptimizer.get_dwba

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

``DWBAResult`` separates the matrix element and observable field amplitudes:

``F_contrast``
   The coherent DWBA contrast structure factor :math:`F_{\Delta,\mathrm{DWBA}}`.
``F_atomic`` and ``F_reference``
   The coherent sums of the actual atomic amplitude and the subtracted planar
   optical-reference amplitude.  Thus
   ``F_contrast == F_atomic - F_reference``.
``contributions``
   An ordered tuple of immutable ``DWBAContribution`` records: bulk first,
   followed by components in crystal order and their generated cells in model
   order.  Each record identifies its component, generated-record index and
   role, optional structural layer, and exposes read-only ``F_atomic``,
   ``F_reference``, and ``F_contrast`` amplitudes.  Roles are ``bulk``,
   ``unit_cell_layer``, ``film_layer``, ``interface_top``,
   ``interface_bottom``, ``surface_termination``, ``covered_film``, and
   ``sharp_film_correction``.
``unperturbed_amplitude``
   The Fresnel reference amplitude :math:`r_0` for specular points and zero
   for non-specular points.
``scattered_amplitude``
   The dimensionless first-Born scattered far-field amplitude

   .. math::

      a_{\mathrm{scattered}}=
      \frac{2\pi i r_e}{k_0\sin\alpha_f\,A_{\mathrm{ref}}}
      F_{\Delta,\mathrm{DWBA}}.
``total_amplitude``
   ``unperturbed_amplitude + scattered_amplitude``.
``structure_factor_squared`` and ``scattered_amplitude_squared``
   The squared moduli of the two corresponding amplitudes.
``differential_cross_section_kernel``
   :math:`r_e^2|F_{\Delta,\mathrm{DWBA}}|^2` per reference lateral cell.  It is
   not an intensity integrated over detector acceptance, footprint,
   coherence, resolution, flux, exposure, or efficiency.
``reflectivity`` and ``first_order_reflectivity``
   The coherent :math:`|r_0+a_{\mathrm{scattered}}|^2` and its formal
   first-order truncation.  These properties require an entirely specular,
   same-polarization, semi-infinite result.

``crystal.dwba.reflectivity(..., polarization="unpolarized")`` evaluates the
``s`` and ``p`` channels independently and averages their reflectivities
incoherently.  The default ``bulk_mode="semi_infinite"`` performs the physical
atom-minus-reference calculation and includes the unperturbed specular
reflection.  Diagnostic ``bulk_mode="unit_cell"`` instead evaluates exactly
one bulk repeat plus every finite generated record, sets every
``F_reference`` to zero, and omits the unperturbed reflection.  Experimental
acceptance integration remains a separate operation.

Comparing a DWBA rod with a kinematical rod
--------------------------------------------

The :doc:`RuO2/TiO2 DWBA diagnostics tutorial
<dwba/dwba_diagnostics>` works through the refractive-index split,
specular Fresnel-plus-contrast amplitude, critical edge, extended specular
rod, and fixed-incidence non-specular rods for flat and rough film models.
Its executed plots compare each DWBA observable with the corresponding
kinematical curve.

``DWBAResult.F_contrast`` and ``SXRDCrystal.F`` are both amplitudes in
electrons per reference lateral cell, and off specular they are directly
comparable.  They are not identical even in the weak-scattering limit, and
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

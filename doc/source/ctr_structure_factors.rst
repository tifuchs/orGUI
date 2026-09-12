CTR Structure-Factor Normalization
==================================

The CTR calculation API uses canonical structure-factor units throughout:
every public ``F`` or ``F_uc`` method returns a complex scattering amplitude
in electrons.

.. warning::

   This convention changes amplitudes produced by older releases that divided
   unit-cell factors by area or volume. Legacy reference amplitudes and fitted
   experimental scale factors must be migrated by the corresponding constant
   normalization factor. Intensities change by the square of that factor.

CTR NeXus angle and index persistence
----------------------------------------

``CTR.toNXdict`` writes the independent H, K, and L arrays in reference-cell
r.l.u. and the six-circle angle arrays in radians. New payloads carry
``@orgui_ctr_schema=2`` and ``sixc_angles/@unit="rad"``. The loader preserves
each angle column for any number of scan points.

``CTR.fromNXdict(payload, angle_units=None)`` reads the declared angle unit
from versioned payloads. For unversioned payloads, the default preserves
numeric angles as radians, because historical orGUI writers incorrectly
labeled their radian data as degrees. To load genuine external degree-valued
data, use ``angle_units="deg"`` explicitly; ``angle_units="rad"`` also provides
an explicit override. ``CTRCollection.fromNXdict`` forwards this option to
every rod. Unit overrides take precedence over stored labels, and conversions
apply only to the six-circle angles, not the structure-factor phase.

Historical writers also copied H into the K array. The loader preserves the
stored K values: there is no generally reliable way to reconstruct the
original K without independent index information. Correct those files using
trustworthy source metadata before scientific analysis.

Measurement quantity, polarization, and scan metadata
------------------------------------------------------

The :doc:`ctr_data_workflow` tutorial demonstrates the preferred loading and
plotting API, NeXus persistence, and explicit migration of legacy ANAROD data.

For new text datasets, use the one-line loader
``CTRCollection.fromCTRFile(path)``. The ``.ctr`` format keeps the useful
ANAROD-style whitespace-separated H, K, L, value, and uncertainty columns but
omits the unused mode column. Leading ``# key: value`` lines declare schema 2,
the stored quantity, polarization state, and column names. An optional trailing
``polarization_factor`` column is retained point by point and split with each
rod. Legacy ANAROD files remain the responsibility of
``CTRCollection.fromANAROD``; the new loader rejects them rather than guessing
missing reduction metadata.

Each ``CTR`` owns a frozen ``MeasurementReduction``. Its ``quantity`` is the
sole definition of what ``CTR.sfI`` and ``CTR.err`` contain:

* ``"structure_factor"`` stores F and its uncertainty in the existing
  arbitrary structure-factor units;
* ``"reflectivity"`` stores the corrected dimensionless field-intensity ratio
  R and its absolute uncertainty.

Legacy constructors and files without this metadata default to structure
factor with unknown polarization provenance. New measurements can record the
incident mixture and outgoing analysis explicitly:

.. code-block:: python

   import numpy as np
   from orgui.datautils.xrayutils.CTRplotutil import (
       CTR,
       CTRScanGeometry,
       MeasurementReduction,
       PolarizationReduction,
   )

   reduction = MeasurementReduction(
       quantity="structure_factor",
       polarization=PolarizationReduction(
           s_fraction=0.7,
           outgoing="unanalysed",
           polarization_factor=np.asarray([0.98, 0.97, 0.96]),
       ),
   )
   rod = CTR(
       (1, 0),
       l,
       F,
       sigma_F,
       reduction=reduction,
       scan_geometry=CTRScanGeometry(
           fixed="in", angle=np.deg2rad(0.3), mirrorx=False
       ),
   )

``s_fraction`` is the incoherent incident s fraction in the local Renaud
basis: one is pure s and zero is pure p. ``outgoing`` is ``"s"``, ``"p"``,
or ``"unanalysed"``. For structure-factor data,
``polarization_factor`` stores the conventional pointwise intensity factor P,
not its reciprocal. ``None`` records P=1. A supplied factor must be finite,
strictly positive, and aligned with the CTR points. Reflectivity cannot carry
this factor because stored R must already contain the selected experimental
correction.

``CTRScanGeometry`` records only the z-mode scan rule. Fixed incidence or exit
uses ``fixed="in"`` or ``"out"`` and an angle in radians in
``(0, pi/2]``. Equal-angle scans use ``fixed="eq", angle=None``. ``mirrorx``
selects the negative-delta scattering branch. Stored six-circle angle records
remain the measured central geometry; the scan rule does not replace them.

Schema-2 NeXus output preserves both metadata records and stores P as a
pointwise dataset. Existing schema-2 payloads that predate these fields still
load with the legacy defaults. ``CTR.fromArray`` and both ANAROD import paths
accept ``reduction=`` and ``scan_geometry=`` overrides.

Copying and point selection preserve the metadata, applying cuts to P and the
six-circle records together with the values and errors. Geometry- or
reduction-aware binning is not defined yet and raises instead of independently
averaging quantities that do not preserve the measurement equation. Existing
kinematical fitting and crystal scaling, difference calculations, collection
scaling, structure-factor resolution convolution, phase/complex-F operations,
ANAROD F export, and symmetry averaging reject reflectivity explicitly. Plot
panels select labels from the stored quantity and F/R datasets cannot share one
axis.

CTR fit predictions and statistics
----------------------------------

Constructing a ``CTROptimizer`` does not evaluate the model. Call
``prepareFit()`` after defining fitted parameters, callbacks, resolution-width
parameters, or displacement constraints. Prediction, residual, likelihood,
and statistics methods reject an unprepared layout with a message directing
the caller back to ``prepareFit()``. Fixed scale-policy and resolution-model
changes only invalidate the current result; the next evaluation refreshes it.

``flat_prediction(x=None, specular=True)`` is the quantity-neutral flattened
prediction API. It returns the final model values after analytical scaling, in
the same rod and point order as the stored data. ``flat_Fcalc`` remains the
F-only spelling and rejects a collection containing reflectivity. The public
``calculated_CTRs`` collection carries the same final prediction arrays and no
copied observation uncertainties. ``Rfactor`` selects only F data and
``Rfactor_R`` selects only R data; a quantity absent from the collection
returns ``None`` without warning.

The built-in kinematical optimizer still accepts structure-factor data only.
These quantity-aware outputs also define the common result contract used by
models that can predict reflectivity.

Call ``optimizer.set_dwba()`` to opt into live DWBA predictions. In that mode,
each rod's polarization reduction and measured six-circle records or z-mode
scan rule determine the field-intensity calculation. Structure-factor rods
receive an effective F prediction using their conventional pointwise P;
corrected-reflectivity rods receive direct R. Stored observations and errors
are never converted between these representations. The complete geometry,
polarization, bulk, and limitation contract is documented in
:ref:`the DWBA fitting section <dwba-fitting>`.

Analytical scale policies are configured explicitly:

.. code-block:: python

   optimizer.set_scale_policies({
       "F": "global",       # one scale shared by the default F group
       (0, 0): "fixed",      # use full ctr_id when hk is ambiguous
   })
   optimizer.prepareFit()
   prediction = optimizer.flat_prediction()
   diagnostics = optimizer.statistics()

``"fixed"`` uses a multiplier of one, ``"scaled"`` fits one analytical scale
per selected rod, and ``"global"`` fits one scale for the nonempty global
group. Analytical scales count as fitted parameters in the reported degrees
of freedom. ``statistics()["covariance"]`` is scaled by the reduced
chi-square, so its diagonal square root equals the reported parameter errors.
When the degrees of freedom are nonpositive or the local covariance is not
estimable, the unavailable statistics and parameter errors are ``None`` and a
warning is emitted. Raw chi-square and the available F/R diagnostics remain
reported. The p-value retains the existing chi-square interpretation and is
therefore heuristic when empirical per-rod weights differ from one.

Signed intensity-to-amplitude conversion
----------------------------------------

``CTR.convertToF()`` converts an intensity in arbitrary F-squared units to a
signed amplitude in place:

.. math::

   g(I) = \operatorname{sign}(I)\sqrt{|I|}.

Negative and zero measurements are retained, and the original intensity is
recoverable as ``g(I) * abs(g(I))``. This signed representation is not a phase
or an unbiased physical amplitude estimate. The conversion adds no
normalization or reflectivity-to-structure-factor physics, and therefore
rejects reflectivity-tagged datasets. Make a deep copy first if the original
intensity arrays are also needed.

For a finite input uncertainty :math:`\sigma_I > 0`, the existing ``err``
array is replaced by the effective symmetric interval half-width

.. math::

   \sigma_{F,\mathrm{eff}}
   = \frac{g(I + \sigma_I) - g(I - \sigma_I)}{2}.

This convention is finite at zero,
:math:`\sigma_{F,\mathrm{eff}}(0)=\sqrt{\sigma_I}`, and approaches ordinary
first-order propagation at high signal-to-noise. It summarizes the transformed
interval around the unchanged central value ``g(I)``; it does not make the
near-zero distribution exactly Normal. The implementation avoids subtracting
nearly equal roots for strong measurements. Invalid uncertainties raise before
the CTR is modified.

By default, nonfinite central intensities are removed while all finite signed
points remain aligned with L, polarization factors, angles, and other
pointwise fields. Pass ``excludeInvalid=False`` to preserve the historical
choice of retaining nonfinite points. Auxiliary ``bgI`` and ``ctrI`` intensity
counters receive the same signed-square-root mapping. The collection method
forwards the option and converts each member in place.

Unit-cell amplitudes
--------------------

``UnitCell.F_uc`` returns

.. math::

   F_{\mathrm{uc}}(\mathbf{Q})
   = \sum_i o_i f_i(\mathbf{Q})
     \exp\left(i\mathbf{Q}\cdot\mathbf{r}_i\right),

including coherent-domain occupancies and displacement factors. No unit-cell
area or volume normalization is applied.

``Film.F_uc`` and ``PoissonSurface.F_uc`` sum their generated layer amplitudes
and return electrons for one lateral unit cell of their source
``UnitCell``.

Semi-infinite bulk amplitude
-----------------------------

``UnitCell.F_bulk`` sums ``UnitCell.F_uc_bulk_direct`` over an infinite stack
of bulk repeats along the out-of-plane direction using a closed-form
geometric series with attenuation:

.. math::

   F_{\mathrm{bulk}}(\mathbf{Q})
   = \frac{F_{\mathrm{uc,bulk}}(\mathbf{Q})}
     {1 - \exp\left(-2\pi i\, l_{\mathrm{bulk}} - \mathrm{atten}\right)},

Here ``atten`` is the empirical, dimensionless amplitude-decay exponent per
reference-cell out-of-plane repeat. The
:doc:`RuO2/TiO2 kinematical-attenuation tutorial
<dwba/kinematic_attenuation>` shows how to derive it from the DWBA
wavevectors at fixed incidence and how to determine the corresponding
minimum reliable exit angle and :math:`L` range.

where :math:`l_{\mathrm{bulk}}` is the out-of-plane reciprocal index *after*
conversion from the reference unit cell via ``refHKLTransform``, i.e. the
third component of ``refHKLTransform @ (h, k, l)``. This is the same index
used to phase every atom in ``F_uc_bulk_direct``, so the phase advance per
bulk repeat in the denominator is consistent with the periodicity actually
being summed.

.. warning::

   Versions up to and including v1.5.0 used the raw, untransformed ``l`` in
   this denominator instead of the reference-transformed index. This was only
   correct when ``refHKLTransform``'s third row equals ``(0, 0, 1)``, i.e.
   only when the bulk cell's own out-of-plane reciprocal axis exactly
   coincides with the reference cell's, in both direction and length. Since
   ``refHKLTransform = B_mat_inv @ rotMatrix @ uc.B_mat`` (see
   ``UnitCell.setReferenceUnitCell``), this held only for the default case of
   a component using its own bulk cell as the reference (no ``reference_uc``
   set). Any explicit ``reference_uc`` whose out-of-plane reciprocal axis
   differs from the bulk's was affected — including a plain scale difference
   between the reference and bulk out-of-plane axis length, not only a
   rotated or reindexed reference.

Distorted-wave Born approximation
----------------------------------

The DWBA calculation, observable-amplitude conventions, comparison with the
kinematical CTR, and worked RuO2/TiO2 tutorial are documented in
:doc:`dwba`.

Surface structure on a rough Film
---------------------------------

A ``PoissonSurface`` is stacked immediately above a ``Film`` and replaces the
Film structure only on the fraction that is truly exposed.  Let
:math:`\theta_i` be the cumulative material occupancy of structural layer
:math:`i`, ordered from bottom to top.  Its exposed surface occupancy is

.. math::

   s_i = \theta_i - \theta_{i+1},

where the occupancy above the highest represented layer is set to zero.  The
part of layer :math:`i` that is covered by another layer remains Film material
with occupancy

.. math::

   c_i = \theta_i - s_i.

The sharp Film already supplies occupancy :math:`\chi_{i<0}` below its nominal
boundary.  Before applying a reconstructed surface structure,
``PoissonSurface`` adds the rough-Film correction

.. math::

   \Delta F_{\mathrm{rough}}
   = \sum_i \left(\theta_i-\chi_{i<0}\right)F_{\mathrm{Film},i}.

The exposed fraction is subsequently replaced by a termination-specific
surface slab as described below.  Film and surface structures can consequently
occur at the same terrace height on complementary lateral fractions.

.. figure:: _static/poisson_surface_occupancies.svg
   :alt: Stacked bars of covered Film and exposed surface occupancy versus structural-layer offset.
   :align: center
   :width: 92%

   Occupancy decomposition for ``PoissonProfile(mean_change=2, alpha=0.5)``.
   Each bar has total height :math:`\theta_i`; its orange segment is the true
   surface fraction :math:`s_i`, and its blue segment is covered Film
   :math:`c_i`.

Termination-specific relaxed surface slabs
-------------------------------------------

A reconstructed surface generally does not repeat the bulk layer cycle.  The
``layer`` column therefore has a narrower meaning for a surface slab: it
selects the Film termination to which the *complete slab* belongs.  Internal
planes of that slab are distinguished by their z coordinates, not by cycling
layer identifiers.  ``UnitCell.as_surface_termination`` makes this explicit by
assigning one termination identifier to every atom and setting
``layer_behavior="select"``.  Selection masks an inactive cell during the
calculation; it does not overwrite atomic occupancy fit parameters.

A ``PoissonSurface`` owns exactly one complete unit cell for every member of
the underlying Film's primitive stacking cycle.  For a two-layer RuO2 Film,
this means two termination cells regardless of how many bulk cells deep each
surface slab is.  The surface-normal length of every termination cell may be
an integer multiple of the Film c axis, while its lateral lattice constants and
lattice angles must match the Film.

Surface-supercell generation example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CTRutil.generate_surface_termination_cells`` automates the affine cycling,
top-plane selection, naming, and whole-slab relabeling.  For example, a
two-cell-deep RuO2 slab and both of its surface terminations can be built as

.. code-block:: python

   from orgui.datautils.xrayutils.CTRdistributions import PoissonProfile
   from orgui.datautils.xrayutils.CTRfilm import PoissonSurface
   from orgui.datautils.xrayutils.CTRutil import (
       generate_surface_termination_cells,
   )

   slab = ruo2_surface.supercell((1, 1, 2), symmetry="independent")
   terminations = generate_surface_termination_cells(
       slab,
       ruo2_film,
       name_template="RuO2_termination_{layer:g}",
   )

   surface = PoissonSurface(
       terminations,
       profile=PoissonProfile(mean_change=-0.55, alpha=0.5),
   )

Here ``ruo2_film`` may be the underlying ``Film``, its ``UnitCell``, its
``LayerCycle``, or simply the ordered tuple ``(1, 2)``.  The number of internal
layers in ``slab`` must be an integer multiple of this primitive Film-cycle
length.  ``symmetry="independent"`` gives repeated copies independent Wyckoff
sites; use ``symmetry="preserve"`` when all repeats should retain shared
symmetry parameters.

The helper chooses which original plane is at the top of each generated slab
and then assigns the complete slab to that Film-cycle state.  Atoms at the
selected top plane can consequently be relaxed independently of their former
internal layer numbers.  The returned termination cells are separate objects
and can carry independent coordinate, displacement-factor, occupancy, and
Wyckoff fit parameters.  The same mapping contract can be consumed by other
surface-roughness models; it is not specific to ``PoissonSurface``.

At exposed terrace :math:`i`, let :math:`t_i` be the corresponding Film-cycle
termination.  CTRcalc generates a bulk-reference slab with the same depth and
alignment as the selected surface slab and adds

.. math::

   \Delta F_{\mathrm{termination}}
   = \sum_i s_i
     \left(
       F_{\mathrm{surface\ slab},t_i}
       - F_{\mathrm{Film\ slab},t_i}
     \right).

Thus a multi-layer relaxed slab replaces, rather than duplicates, the same
depth of Film material.  If every termination slab is identical to its Film
reference, this term cancels numerically in ``F_uc``, ``zDensity_G``, and the
optical profile.  A zero-width distribution selects and replaces only the
single sharp termination.

The immediately underlying ``Film`` and its current stacking phase are inferred
when ``SXRDCrystal`` applies stacking; no Film reference is serialized.
Termination banks are stored in ``.xtal`` and ``.xpr`` files as
``TerminationUnitCell <layer> <name>`` sections.  A legacy single ``UnitCell``
surface remains readable and is expanded automatically into one termination
variant per Film-cycle state.

Epitaxy-interface strain coupling and offset
---------------------------------------------

An ``EpitaxyInterface`` can join materials with different lateral areas.  Its
canonical lateral cell is the lower unit cell.  Internally it combines the
upper and lower amplitudes as

.. math::

   F_{\mathrm{interface}}
   = A_{\mathrm{lower}}
     \left(
       \frac{F_{\mathrm{upper}}}{A_{\mathrm{upper}}}
       + \frac{F_{\mathrm{lower}}}{A_{\mathrm{lower}}}
     \right).

The result is therefore in electrons per lower interface cell.

The lower material is a correction to the semi-infinite bulk, whose atoms
remain at their unstrained positions below the nominal interface boundary.
For every represented lower-material layer, the interface therefore adds its
distributed occupancy at the strain-field position and subtracts the
corresponding sharp-bulk occupancy at the original bulk-lattice position.
This remains true when the statistical support extends below zero.  The
subtraction is deliberately *not* strained: it removes the density already
present in the sharp, semi-infinite bulk before the distributed replacement
is added.

Coordinate model
^^^^^^^^^^^^^^^^

The dimensionless ``strain_coupling`` parameter, written :math:`\kappa`,
controls the out-of-plane lattice transition, with
:math:`0\leq\kappa\leq1`.  It does not control whether film and bulk
amplitudes scatter coherently; those amplitudes are always summed as complex
structure factors.  Let :math:`z_{0,m}` be an atom's position on the
independent lattice of material :math:`m`, and let :math:`z_1` be its position
in the fully strain-coupled, zero-offset field.  Strain coupling alone gives
the linear interpolation

.. math::

   z_m(\kappa,0)
   = z_{0,m} + \kappa\left(z_1-z_{0,m}\right).

For the fully strain-coupled discrete profile, every structural layer supplies
its own physical interval.  If :math:`d_{\mathrm{top},i}` and
:math:`d_{\mathrm{bottom},i}` are the upper- and lower-lattice intervals for
layer :math:`i`, the coupled interval is

.. math::

   d_{1,i}
   = P_i d_{\mathrm{top},i}
     + (1-P_i)d_{\mathrm{bottom},i},

and the layer origins are the cumulative sum of these intervals.  The
calculation must therefore accumulate all structural layers independently;
one layer's strain cannot represent a complete multi-layer unit-cell cycle.

The ``offset`` parameter :math:`o` is expressed in fractional lower-bulk
:math:`c` coordinates.  Its physical displacement is

.. math::

   \Delta = o\,c_{\mathrm{bulk}}.

Positive :math:`o` points toward increasing :math:`z`.  Let :math:`P_i` be
the upper-material occupancy at represented interface layer :math:`i`.  The
profile used for the offset field is normalized over the complete represented
statistical support,

.. math::

   \widehat P_i
   = \frac{P_i-P_{\mathrm{deep\ bulk}}}
          {P_{\mathrm{film}}-P_{\mathrm{deep\ bulk}}},
   \qquad
   \widehat P_{\mathrm{deep\ bulk}}=0,\quad
   \widehat P_{\mathrm{film}}=1.

Using this occupancy coordinate, the generated upper- and lower-material
positions are

.. math::

   z_{\mathrm{top}}(\kappa,o)
   &= z_{0,\mathrm{top}}
      + \kappa(z_1-z_{0,\mathrm{top}})
      + \Delta\left[(1-\kappa)+\kappa\widehat{P}(z)\right], \\
   z_{\mathrm{bottom}}(\kappa,o)
   &= z_{0,\mathrm{bottom}}
      + \kappa(z_1-z_{0,\mathrm{bottom}})
      + \Delta\kappa\widehat{P}(z).

The offset terms have a useful decomposition:

.. math::

   u_{\mathrm{top}}(z)
   &= \Delta(1-\kappa) + \kappa\Delta\widehat P(z),\\
   u_{\mathrm{bottom}}(z)
   &= \kappa\Delta\widehat P(z).

The first term is a rigid registry offset between two independent lattices.
The second is a shared displacement field that expands or contracts the
interface as a whole.  Their difference is independent of :math:`z`,

.. math::

   u_{\mathrm{top}}-u_{\mathrm{bottom}}=(1-\kappa)\Delta,

while the shared displacement accumulated from deep bulk to film is

.. math::

   u_\kappa(\mathrm{film})
   -u_\kappa(\mathrm{deep\ bulk})
   =\kappa\Delta.

Equivalently, wherever the occupancy profile is differentiable, the additional
local strain due to the offset is

.. math::

   \epsilon_o(z)
   = \frac{\mathrm{d}u_C}{\mathrm{d}z}
   = \kappa\Delta\frac{\mathrm{d}\widehat P}{\mathrm{d}z}.

This separates three named regimes:

.. list-table::
   :header-rows: 1
   :widths: 18 25 32 25 22

   * - ``strain_coupling``
     - Regime
     - Atomic positions
     - Residual registry offset
     - Shared expansion
   * - ``0``
     - Independent-lattice limit
     - Both materials retain their independent ideal lattice spacings.
     - :math:`\Delta`
     - :math:`0`
   * - ``0.5``
     - Partially strain-coupled interface
     - Half of the strain and offset field is applied.
     - :math:`\Delta/2`
     - :math:`\Delta/2`
   * - ``1``
     - Fully strain-coupled interface
     - Both distributed materials follow the same strain field.
     - :math:`0`
     - :math:`\Delta`

Thus :math:`\kappa=0` is the **independent-lattice limit**,
:math:`0<\kappa<1` is a **partially strain-coupled interface**, and
:math:`\kappa=1` is the **fully strain-coupled interface**.  In the
independent-lattice limit the lower material remains on its bulk lattice and
the upper material remains on its own lattice, translated rigidly by
:math:`\Delta`; there is no strain field.  In the fully strain-coupled limit,
coincident lattices receive identical occupancy-mediated displacements and
therefore remain coincident.  Intermediate strain coupling continuously
partitions :math:`\Delta` between rigid registry and shared interface strain.

The strain-coupled field is anchored to the unstrained deep bulk at the lower
edge of the represented statistical support.  Its accumulated displacement is
propagated upward rather than removed by re-anchoring at the nominal boundary.
The complete Film and all subsequently stacked surface, water, and other
components receive the full :math:`\Delta` translation without a thickness
change.  The fixed sharp-bulk subtraction, semi-infinite bulk, and nominal
statistical boundary remain unchanged.

The conventional epitaxial *relaxation degree* uses the opposite endpoint
direction: :math:`R=0` denotes a pseudomorphic, fully strained layer and
:math:`R=1` a layer at its relaxed lattice constant.  Consequently,
``strain_coupling`` is qualitatively analogous to :math:`1-R`, but it is not
the same model: :math:`\kappa` controls an occupancy-mediated, depth-dependent
out-of-plane coordinate field shared by both interface materials, rather than
a uniform in-plane film lattice parameter.  See `Zhylik et al., Journal of
Applied Crystallography 46, 919--925 (2013)
<https://doi.org/10.1107/S0021889813011907>`_ for the conventional relaxation
definition and its diffraction treatment.

Newly serialized interfaces use the explicit four-column header
``Width/cells Skew/cells StrainCoupling Offset/bulk_frac``.  ``.xtal`` and
``.xpr`` models written by orGUI v1.5.0 used the two-column header
``Width/cells Skew/cells``; these records remain readable and are interpreted
as :math:`\kappa=1` and :math:`o=0`, preserving the fully strain-coupled,
zero-offset model.  Saving the reconstructed model writes the current
four-column form.

The optional companion ``.h5`` file contains fit definitions rather than a
separate crystal model.  A v1.5.0 two-element interface baseline stored there
is likewise expanded to ``[Width, Skew, 1, 0]`` when loaded.  Existing Width
and Skew fit parameters retain their original indices; ``strain_coupling`` and
``offset`` become fit parameters only when explicitly selected.

RuO2 on TiO2 example
^^^^^^^^^^^^^^^^^^^^

The following density profiles use a 15 nm RuO2 film on TiO2 with a 1 nm
Skellam interface width.  The columns use offsets :math:`o=0`, :math:`0.02`,
and :math:`0.17`.  For the TiO2 lower-bulk lattice
:math:`c_{\mathrm{bulk}}=6.5807` Angstrom, these correspond to
:math:`\Delta=0`, :math:`0.13`, and :math:`1.12` Angstrom.  Each column
vertically separates the :math:`\kappa=0`, :math:`\kappa=0.5`, and
:math:`\kappa=1` densities.
Gray is the summed density, blue is the TiO2 contribution including the
semi-infinite bulk and lower-interface correction, and orange is the RuO2
interface-plus-film contribution.

.. figure:: _static/epitaxy_strain_coupling_offset_density.png
   :alt: RuO2 on TiO2 interface densities and shared displacement profiles for three offsets and three strain-coupling values.
   :align: center
   :width: 100%

   Strain coupling partitions the physical offset into rigid registry and
   shared interface expansion.  Arrows above the density panels indicate the
   full bulk-to-film displacement :math:`\Delta`.  The lower panels show the
   shared displacement :math:`u_\kappa(z)` obtained from the generated
   lower-material domain transforms; its film-side plateau is
   :math:`\kappa\Delta`.

The :math:`\kappa=0` density peaks show two independently periodic lattices:
the TiO2 peaks retain the lower-bulk spacing, while the RuO2 peaks retain the
RuO2 film spacing and are shifted rigidly by :math:`\Delta`.  No density
feature is bent into a strain field.  At :math:`\kappa=1` the offset is
accumulated through the occupancy profile, so both distributed materials
follow the same displacement and the lower-panel curve rises from zero to
:math:`\Delta`.  The :math:`\kappa=0.5` row is the exact intermediate case:
half of the offset remains as registry separation and half appears as a
measurable expansion of the interface.

The lower curves are plotted at the represented lower-material layer
positions rather than inferred from peak maxima.  They therefore provide a
direct diagnostic of the coordinate transforms used by ``F_uc``,
``zDensity_G``, and the optical profile.

Crystal composition
-------------------

``SXRDCrystal`` automatically uses its bulk unit cell as the reciprocal-space
and lateral-area reference unless ``reference_uc`` is supplied explicitly.
The constructor propagates that reference to every source and generated layer
unit cell.

For each crystal component :math:`j`, ``SXRDCrystal.F`` evaluates

.. math::

   F_{\mathrm{crystal}}
   = \frac{A_{\mathrm{ref}}}{A_{\mathrm{bulk}}}F_{\mathrm{bulk}}
     + \sum_j
       w_j d_j
       \frac{A_{\mathrm{ref}}}{A_j}F_j ,

where:

* :math:`A_{\mathrm{ref}}` is ``reference_uc.uc_area``;
* :math:`A_j` is the component ``uc_area``;
* :math:`w_j` is the dimensionless crystal-component weight;
* :math:`d_j` is a dimensionless coherent-domain occupancy.

Thus ``SXRDCrystal.F`` returns electrons per reference lateral cell and is
invariant when a component is replaced by an equivalent in-plane supercell.
No illuminated footprint, detector response, or experimental scale factor is
included. Calculated intensity is proportional to
:math:`|F_{\mathrm{crystal}}|^2`.

Amplitudes, squared structure factors, and reflectivity
--------------------------------------------------------

Four quantities appear along the CTR path and are deliberately not
interchangeable.

``F``
   The coherent complex structure factor returned by
   :meth:`~orgui.datautils.xrayutils.CTRcalc.SXRDCrystal.F`, in electrons per
   lateral cell of the reference unit cell.  It carries a phase.

``F2``
   The squared structure factor
   :meth:`~orgui.datautils.xrayutils.CTRcalc.SXRDCrystal.F2`, equal to
   ``abs(F) ** 2`` for a coherent crystal, in the squared units of ``F``.  It
   is real and nonnegative and carries no phase.

``r``
   The optical reflection amplitude of the DWBA path, a dimensionless complex
   ratio formed from the electron-density profile.

``R``
   Reflectivity, the dimensionless intensity ratio ``abs(r) ** 2``.

None of these is detector counts.  Incident flux, illuminated footprint,
polarization and Lorentz factors, detector response, acquisition time,
background, and the fitted experimental scale all sit outside them.

An incoherent model has a ``F2`` but no ``F``: a mixed state has no unique
complex amplitude, so requesting one is a category error rather than a
missing feature.  ``F2`` never changes meaning with the optimizer mode; a
reflectivity is returned from a separately named boundary, never from a
method called ``F2``.

Incoherent height domains
---------------------------------------------

A :class:`~orgui.datautils.xrayutils.CTRfilm.PoissonSurface` describes a
distribution of surface heights.  By default those heights are added as
*amplitudes*, which is the coherent limit: every height lies inside one
coherence patch, and the height distribution enters through the squared
modulus of its characteristic function.

When the lateral height domains are large compared with the projected
coherence area, each patch instead sees a single flat height and the patches
add in intensity.  The two limits are

.. math::

   F^2_{\mathrm{coherent}} = \left| \sum_n p_n A_n \right|^2,
   \qquad
   F^2_{\mathrm{incoherent}} = \sum_n p_n \left| A_n \right|^2,

where :math:`p_n` is the probability of height state :math:`n` and
:math:`A_n` is the **complete** crystal amplitude when that height covers the
patch: bulk, Film, the flat-height Film correction, and the exposed
termination together.  Averaging only the surface correction would drop the
bulk-surface and Film-surface interference inside each domain and is a
different quantity.

:class:`~orgui.datautils.xrayutils.CTRincoherent.PoissonHeightDomains`
interpolates between them with a dimensionless ``incoherent_fraction``
:math:`\kappa` in ``[0, 1]``:

.. math::

   F^2_\kappa = (1 - \kappa) F^2_{\mathrm{coherent}}
                + \kappa F^2_{\mathrm{incoherent}}.

For a patch containing :math:`N_{\mathrm{eff}}` independent equal-area
domains, :math:`\kappa \approx 1 / N_{\mathrm{eff}}`.  A CTR-only fit cannot
separate domain size from beam coherence length; AFM, transverse scans,
rocking widths, or reciprocal-space maps are needed to constrain one of them.

Constructing a ``PoissonSurface`` does **not** opt a calculation into
incoherent averaging.  Wrapping the crystal is the explicit opt-in, and
existing scripts, saved crystals, and optimizer setups stay coherent::

    from orgui.datautils.xrayutils.CTRincoherent import PoissonHeightDomains
    from orgui.datautils.xrayutils.CTRopt import CTROptimizer

    # The crystal on its own is unchanged and fully coherent.
    coherent = crystal.F2(h, k, l)

    domains = PoissonHeightDomains(
        crystal,
        surface="rough_surface",
        incoherent_fraction=0.35,
    )
    mixed = domains.F2(h, k, l)

The fraction is a fixed setting until it is explicitly added as a fit
parameter, after which it occupies the first entry of the model block::

    domains.addFitParameter(
        "incoherent_fraction",
        limits=(0.0, 1.0),
        name="rough_surface incoherent_fraction",
    )

    fit = CTROptimizer(domains, ctrs)
    fit.prepareFit()
    assert fit.xtal is fit.model.coherent_model

The wrapper is passed through the optimizer's existing model argument.
``optimizer.xtal`` remains the coherent crystal, so registered fit callbacks
and displacement constraints keep receiving an ``SXRDCrystal``, while
``optimizer.model`` is the fitted forward model.  Resolution is applied to
``F2`` before the conversion back to a stored ``|F|``.

The first implementation is kinematical only.  Combining an incoherent
``F2`` model with DWBA raises during ``prepareFit``: changing the flat
surface height changes the optical reference profile and its internal fields,
so a correct DWBA ensemble would have to prepare and evaluate each height
separately and mix ``abs(r) ** 2``.

The target surface does not have to be the topmost component: a water layer
or a cap may be stacked above it.  Anything above the surface is placed once
at the surface's *mean* height and is common to every domain, which is how
the coherent model already treats it, so the :math:`\kappa = 0` endpoint is
unchanged.  Note the approximation this carries: in a strict large-domain
limit an overlayer would follow each domain's own height rather than the
mean.

Height states are indexed by the structural layer :math:`n` of the top filled
layer, and the mass of that state is ``probability(n + 1)``: layer :math:`n`
is the top filled layer exactly when the signed height change equals
:math:`n + 1`.  The retained interval is chosen from the calculated
probability masses and the profile's tail target, never from measured CTR
values, and is renormalized once.  ``exact_layer_count`` (default 10) retains
every state for a narrow distribution; it is a policy switch, not a cap.

API reference
-------------

.. autoclass:: orgui.datautils.xrayutils.CTRplotutil.PolarizationReduction

.. autoclass:: orgui.datautils.xrayutils.CTRplotutil.MeasurementReduction

.. autoclass:: orgui.datautils.xrayutils.CTRplotutil.CTRScanGeometry

.. autoclass:: orgui.datautils.xrayutils.CTRdistributions.SurfaceProfile
   :members: support, occupancy, correction, surface_occupancy
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRdistributions.PoissonProfile
   :members:
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRcalc.SXRDCrystal
   :members: F, F2, F_surf, evaluate_kinematic, setGlobalReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRincoherent.IncoherentModel
   :members: coherent_model, setParameters, validate, to_config
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRincoherent.IncoherentF2Model
   :members: F2
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRincoherent.CoherentStateEnsembleModel
   :members: incoherent_fraction, accumulate_states
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRincoherent.PoissonHeightDomains
   :members: incoherent_fraction, target_surface
   :member-order: bysource

.. autofunction:: orgui.datautils.xrayutils.CTRincoherent.available_incoherent_models

.. autofunction:: orgui.datautils.xrayutils.CTRincoherent.create_incoherent_model

.. autofunction:: orgui.datautils.xrayutils.CTRincoherent.create_incoherent_model_from_config

.. autoclass:: orgui.datautils.xrayutils.CTRuc.UnitCell
   :members: F_uc, F_bulk, setReferenceUnitCell, supercell, affine_layer_transform, as_surface_termination
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.Film
   :members: F_uc, uc_area, setReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.EpitaxyInterface
   :members: F_uc, uc_area, setReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.PoissonSurface
   :members: F_uc, uc_area, setReferenceUnitCell, flat_domain_corrections
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.FlatHeightCorrections
   :members: excluded_probability, iter_states, as_array
   :member-order: bysource

.. autofunction:: orgui.datautils.xrayutils.CTRutil.generate_surface_termination_cells

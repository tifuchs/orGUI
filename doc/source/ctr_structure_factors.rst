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

An ``EpitaxyInterface`` can join materials with different lateral areas. Its
canonical lateral cell is the lower unit cell. Internally it combines the
upper and lower amplitudes as

.. math::

   F_{\mathrm{interface}}
   = A_{\mathrm{lower}}
     \left(
       \frac{F_{\mathrm{upper}}}{A_{\mathrm{upper}}}
       + \frac{F_{\mathrm{lower}}}{A_{\mathrm{lower}}}
     \right).

The result is therefore in electrons per lower interface cell.

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

API reference
-------------

.. autoclass:: orgui.datautils.xrayutils.CTRcalc.SXRDCrystal
   :members: F, F_surf, setGlobalReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRuc.UnitCell
   :members: F_uc, F_bulk, setReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.Film
   :members: F_uc, uc_area, setReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.EpitaxyInterface
   :members: F_uc, uc_area, setReferenceUnitCell
   :member-order: bysource

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.PoissonSurface
   :members: F_uc, uc_area, setReferenceUnitCell
   :member-order: bysource

Wyckoff Symmetry And Fitting
============================

The CTR symmetry helpers make Wyckoff information available as metadata on
ordinary ``UnitCell`` atoms. Structure-factor calculations still use normal
``UnitCell.basis`` coordinates; symmetry is used only for queries, file
persistence, and optional fit-parameter construction.

The symmetry construction path uses optional PyXtal support. Importing
``CTRcalc``, ``CTRuc``, or ``CTRsymmetry`` does not import PyXtal. PyXtal is
only required when assigning symmetry from a structure seed or querying
Wyckoff tables. Files that already contain saved symmetry metadata can be
loaded and fitted without PyXtal.

Building A Surface Symmetry Model
---------------------------------

Use ``CTRsymmetry.model_from_seed`` to assign parent-crystal Wyckoff sites and
then build a commensurate CTR surface unit cell:

.. code-block:: python

   import numpy as np
   from orgui.datautils.xrayutils import CTRsymmetry

   surface_spec = CTRsymmetry.SurfaceCellSpec(
       parent_a=(4.653255, 4.653255, 2.9692),
       parent_alpha=(90.0, 90.0, 90.0),
       transform=np.asarray(
           [
               [1.0, 0.0, 1.0],
               [-1.0, 0.0, 1.0],
               [0.0, 1.0, 0.0],
           ]
       ),
       layer_origins=(0.0, 0.5),
       translation_range=1,
   )

   model = CTRsymmetry.model_from_seed(parent_structure, surface_spec, tol=1e-3)
   unitcell = model.build_unitcell("RuO2")

The ``transform`` columns are surface-cell vectors expressed in parent
conventional fractional coordinates. Parent coordinates are mapped to the
surface cell as:

.. math::

   p_{\mathrm{surface}} = S^{-1}(p_{\mathrm{parent}} - p_{\mathrm{origin}}).

If PyXtal standardizes or reorients the input structure, the surface transform
must be expressed in that PyXtal parent basis.

Querying Wyckoff Metadata
-------------------------

``UnitCell`` exposes resolved metadata:

.. code-block:: python

   unitcell.wyckoff_sites()
   unitcell.wyckoff_couplings("O_4f")
   unitcell.wyckoff_site_couplings("O_4f")
   unitcell.atom_wyckoff_metadata(0)

``wyckoff_sites()`` returns site dictionaries with the site id, element,
Wyckoff label, free variables, generated atom indices, space group, and status.
The status is one of:

``metadata_only``
   Symmetry metadata is present, but no Wyckoff-related fit parameter has been
   added.

``symmetry_preserving``
   A true Wyckoff variable fit was added with ``addWyckoffParameter``.

``site_displaced``
   A representative-site shift was added with ``addWyckoffShift``. This moves
   all atoms in the assigned site together, but may lower the crystallographic
   site symmetry.

``partially_overridden``
   A manual atom-coordinate parameter touches atoms in the site.

Two Fit Modes
-------------

``addWyckoffParameter`` fits a true Wyckoff variable such as ``u``, ``v``, or
``w``. It preserves the assigned Wyckoff-site symmetry because all equivalent
atom coordinates are updated through the stored affine Wyckoff couplings:

.. code-block:: python

   unitcell.addWyckoffParameter(
       "O_4f",
       "u",
       absolute_limits=(0.28, 0.34),
   )

The fitted value is a delta from the stored variable value. For rutile oxygen,
coordinates such as ``u``, ``1-u``, ``0.5+u``, and ``0.5-u`` therefore move
with the correct signs.

``addWyckoffShift`` fits a parent conventional fractional displacement of the
representative site coordinate:

.. code-block:: python

   unitcell.addWyckoffShift(
       "O_4f",
       "x",
       limits=(-0.02, 0.02),
   )

Here ``"x"``, ``"y"``, and ``"z"`` are parent conventional fractional axes,
not generated surface-cell axes. The shift is propagated through the stored
space-group operation for each generated atom and then through the surface-cell
transform. This intentionally keeps the original atom assignment and atom count
while allowing the site to move away from its exact Wyckoff constraint.

Plural helpers are available for fitting several variables or shifts:

.. code-block:: python

   unitcell.addWyckoffParameters("Mg_8d")          # all free variables
   unitcell.addWyckoffShifts("Mg_8d", axes=("x", "z"))

Container Classes
-----------------

``Film`` and ``PoissonSurface`` forward Wyckoff helper calls to their single
template unit cell:

.. code-block:: python

   film.addWyckoffParameter("O_4f", "u", limits=(-0.05, 0.05))
   poisson.addWyckoffShift("O_4f", "x", limits=(-0.02, 0.02))

``EpitaxyInterface`` follows the same selector convention as its existing
``addFitParameter`` API. Select the internal unit cell with ``unitcell``:

.. code-block:: python

   interface.addWyckoffParameter(
       "O_4f",
       "u",
       limits=(-0.05, 0.05),
       unitcell="top",
   )

   interface.addWyckoffShift(
       "O_4f",
       "x",
       limits=(-0.02, 0.02),
       unitcell=["top", "bottom"],
   )

When a list is supplied, one parameter is added to each selected internal unit
cell. They remain separate parameters unless a higher-level ``SXRDCrystal``
coupled parameter links them.

Linking Across Unit Cells
-------------------------

``SXRDCrystal`` provides the same two Wyckoff helper modes at the crystal
level. These helpers create one coupled crystal parameter and ordinary
relative fit parameters inside each selected unit cell.

For a symmetry-preserving Wyckoff variable:

.. code-block:: python

   crystal.addWyckoffParameter(
       {
           "RuO2": ("O_4f", "u"),
           "TiO2": ("O_4f", "u"),
       },
       name="shared_rutile_oxygen_u",
       limits=(-0.05, 0.05),
   )

For a representative-site shift:

.. code-block:: python

   crystal.addWyckoffShift(
       {
           "RuO2": ("O_4f", "x"),
           "TiO2": ("O_4f", "x"),
       },
       name="shared_oxygen_x_shift",
       limits=(-0.02, 0.02),
   )

For ``EpitaxyInterface`` components, select the internal unit cell in the
value dictionary, or use a ``(component, unitcell)`` key:

.. code-block:: python

   crystal.addWyckoffParameter(
       {
           "interface": {
               "site_id": "O_4f",
               "variable": "u",
               "unitcell": ("top", "bottom"),
           }
       },
       name="shared_interface_u",
       limits=(-0.05, 0.05),
   )

   crystal.addWyckoffShift(
       {("interface", "top"): ("O_4f", "x")},
       name="top_interface_oxygen_shift",
       limits=(-0.02, 0.02),
   )

Persistence
-----------

Symmetry-bearing ``.xtal`` and ``.xpr`` files store resolved plain-text tables
after the layer metadata. The persisted information includes:

* space group number and symbol;
* surface transform and origin;
* Wyckoff site ids, labels, variables, and representative parent coordinates;
* generated atom metadata;
* Wyckoff variable couplings;
* representative-site shift couplings.

Reloaded files can be queried and fitted with ``addWyckoffParameter`` and
``addWyckoffShift`` without PyXtal because the resolved couplings are already
present in the file.

API Reference
-------------

.. automodule:: orgui.datautils.xrayutils.CTRsymmetry
   :members: SurfaceCellSpec, SurfaceSymmetryModel, possible_wyckoff_positions,
             sites_from_seed, model_from_seed, surface_unitcell_from_seed
   :member-order: bysource
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRuc.UnitCell
   :members: wyckoff_sites, wyckoff_couplings, wyckoff_site_couplings,
             atom_wyckoff_metadata, addWyckoffParameter, addWyckoffParameters,
             addWyckoffShift, addWyckoffShifts
   :member-order: bysource
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.Film
   :members: addWyckoffParameter, addWyckoffParameters, addWyckoffShift,
             addWyckoffShifts
   :member-order: bysource
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.EpitaxyInterface
   :members: addWyckoffParameter, addWyckoffParameters, addWyckoffShift,
             addWyckoffShifts
   :member-order: bysource
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRfilm.PoissonSurface
   :members: addWyckoffParameter, addWyckoffParameters, addWyckoffShift,
             addWyckoffShifts
   :member-order: bysource
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRcalc.SXRDCrystal
   :members: addWyckoffParameter, addWyckoffShift
   :member-order: bysource
   :no-index:

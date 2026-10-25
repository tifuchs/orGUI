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

Enumerating Bulk Surface Terminations
-------------------------------------

``CTRsymmetry.surface_terminations`` returns inspectable bulk-cut candidates.
It accepts either Miller indices in the **input CIF setting**, with an optional
preferred in-plane direction, or an explicit full surface matrix:

.. code-block:: python

   from orgui.datautils.xrayutils import CTRsymmetry

   options = CTRsymmetry.surface_terminations(
       "TbScO3_published_Pnma.cif",
       miller=(1, 0, 1),       # (101) Pnma corresponds to (110) Pbnm
       in_plane=(1, 0, -1),
       layer_tolerance=0.25,  # Angstrom along the true surface normal
   )

   # Alternative: keep this exact basis and repeat; no cell reduction.
   matrix_options = CTRsymmetry.surface_terminations(
       "TbScO3_published_Pnma.cif",
       transform=[[1, 0, 1], [0, 1, 0], [-1, 0, 1]],
       layer_tolerance=0.25,
   )

   for option in options:
       print(option.id, option.cut_interval, option.equivalent_cuts)
       print([layer["composition"] for layer in option.layer_sequence])

   # Select a candidate after inspecting its composition and stacking.
   substrate = options[0].build_unitcell(
       "TbScO3_surface",
       displacement="isotropic_equivalent",
   )
   print(substrate.stacking_translation)  # (sx, sy, 1), fractional in this cell
   print(substrate.repeat_cartesian)       # physical crystal repeat in Angstrom

Miller-index construction defaults to ``cell_mode="surface"``: a and b lie in
the surface plane and c is perpendicular to both. Alpha and beta are 90 degrees;
gamma is determined by the in-plane lattice and need not be 90 degrees. This
changes the coordinate basis without straining the crystal. The physical upward
repeat is ``a*sx + b*sy + c``, stored as ``stacking_translation=(sx, sy, 1)``.
The coordinate c alone need not be a crystal lattice translation.

Use ``cell_mode="lattice"`` with Miller indices to retain the original oblique
periodic cell. Explicit ``transform`` matrices default to lattice mode;
``cell_mode="surface"`` explicitly projects their third vector onto the normal.

``transform`` columns are direct surface vectors in parent fractional
coordinates. For this enumeration API it must be an integer 3-by-3 matrix with
positive determinant. Do not combine it with ``miller`` or ``in_plane``.
Automatic orientation can choose a different basis or repeat from a previously
hand-written matrix, so use the matrix form to preserve existing rod indexing.

Reciprocal indices change with the coordinate basis. To compare the same
physical scattering points in the two modes, convert parent-cell indices using
``hkl_cell = option.surface_spec.transform.T @ hkl_parent`` (three-row arrays).
Do not reuse unchanged numeric HKL arrays when migrating an existing fit.
``refHKLTransform`` still maps reference indices into each cell's coordinates;
the repeat phase is then evaluated once in those coordinates.

Atoms are grouped by adjacent height separations along the true surface normal,
including across the periodic boundary. ``layer_tolerance`` is the maximum
separation within a group, not a maximum total group thickness. A cut is placed
at the midpoint of each remaining gap. Distorted or buckled chemical layers
can split into several groups; inspect the result and choose an appropriate
tolerance rather than assuming equally spaced or chemically defined planes.

``unique=True`` merges candidates only under bulk symmetries preserving the
oriented surface normal, species, occupancies and known displacement data.
Anisotropic tensors are transformed when testing equivalence. The complete
stacking participates in this check, not only the outermost composition.
``unique=False`` returns all candidate gaps before merging. This enumerates
bulk-truncated cuts for the requested cell and grouping, not reconstructed
surfaces, adsorption, or lateral film/substrate registries.

Each candidate exposes:

* ``surface_spec``: resolved parent lattice, coordinate transform, origin and
  layer boundaries. In surface mode ``bulk_transform`` retains the original
  integer periodic basis for symmetry enumeration. Parent translation bounds
  are calculated automatically.
* ``cut_interval``: the open gap in unshifted oriented fractional z. Its upper
  endpoint can exceed 1 for a wraparound gap. ``equivalent_cuts`` includes all
  gaps merged into the candidate, including the representative.
* ``layer_sequence``: layers from the surface inward, with occupancy-weighted
  ``composition``, mean ``height_angstrom``, ``spacing_to_next_angstrom``, and
  individual ``atoms``. The final spacing continues into the next bulk repeat.
* ``outermost_atoms``: the first layer's atom dictionaries, with source labels,
  occupancies, fractional coordinates and displacement properties.

The candidate spec describes a wrapped 0 <= z < 1 repeat. ``build_unitcell``
returns an independent cell translated to -1 <= z < 0, with the selected cut at
z=0 and crystal on its negative side. Layer heights use this same interface
origin and perpendicular distance, even when C is oblique. The built cell
retains Wyckoff metadata. In surface mode this translation follows the full
repeat, including its lateral component. Fractional x and y may lie outside
[0, 1); preserving these representatives preserves phase at fractional in-plane
indices. Candidate IDs are deterministic for unchanged inputs
and settings; they are not persistent identifiers across changes of setting,
cell, or grouping tolerance. No preferred termination or surface energy is
inferred.

Preserving CIF Source Properties
--------------------------------

``CTRsymmetry.structure_from_cif(path, block=None)`` uses the existing pymatgen
dependency and retains partial occupancies, source labels, published Ueq/B,
and anisotropic displacement tensors. Bibliography-only blocks are ignored;
multiple structural blocks require ``block="name"``. The returned pymatgen
Structure carries ``B_iso`` and ``U_cart`` site properties in Angstrom squared,
plus ``cif_label``, ``cif_block`` and ``displacement_anisotropic`` provenance.
``U_cart`` is Cartesian and rotated for every symmetry image. Published Ueq
is retained independently of the tensor-derived equivalent value.

Use this loader before modifying the lattice, and retain the site properties
when constructing or copying structures. Ordinary pymatgen Structures can
also provide ``B_iso``, ``U_iso`` or ``U_cart`` properties explicitly. Missing
values remain missing. Unknown, nonpositive or overfull CIF occupancies and
mixed-species sites are rejected rather than normalized or reduced to a single
element. Values that differ within a geometrically assigned Wyckoff orbit
raise an error instead of using just one representative.

Source occupancy and displacement values are now inherited by
``sites_from_seed``, ``model_from_seed`` and ``surface_unitcell_from_seed``.
``iDW`` and ``oDW`` are B parameters in Angstrom squared; U is converted using
:math:`B = 8\pi^2 U`. Explicit ``occ``, ``iDW`` and ``oDW`` still override the
source. The historical missing-data fallback B=0.5 is retained by these older
helpers, with a warning; ``fallback_B=None`` disables it. New termination
builders require an explicit ``fallback_B`` if source displacement data are
absent. Both explicit B overrides can also be supplied.

Anisotropic data require ``displacement="isotropic_equivalent"`` or explicit
in-plane and out-of-plane B overrides. The equivalent policy uses published
Ueq/B when present, otherwise one third of the Cartesian tensor trace, and
assigns the resulting B to both directions. General anisotropic scattering is
not implemented by this option. Full tensors remain available on the source
and candidates for inspection and equivalence checks, in the input structure's
Cartesian frame. They are provenance, not surface-frame anisotropic scattering
parameters. Explicit ``iDW``/``oDW`` overrides in surface mode are referenced to
the surface plane/normal. Isotropic-equivalent B is independent of this rotation.

The example ``examples/CTR/surface_stacking_example.py`` accepts a CIF path and
compares both modes at identical physical scattering points. It defaults to
(101), suitable for the published TbScO3 Pnma setting.

This corrects previously silent changes to partial occupancies and displacement
parameters and can change computed intensities. To reproduce the old uniform
assignment deliberately, pass ``occ=1.0, iDW=0.5, oDW=0.5``.

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

The transform is expressed in the input parent basis. When PyXtal standardizes
or reorients the structure, orGUI maps the assigned sites back to that basis.
Use ``translation_range=None`` to calculate complete per-atom translation
bounds automatically instead of specifying a search range manually.

Querying Wyckoff Metadata
-------------------------

``UnitCell`` exposes resolved metadata:

.. code-block:: python

   unitcell.wyckoff_sites()
   unitcell.wyckoff("O_4f")
   unitcell.wyckoff_couplings("O_4f")
   unitcell.wyckoff_site_couplings("O_4f")
   unitcell.atom_wyckoff_metadata(0)

``wyckoff_sites()`` returns site dictionaries with the site id, element,
Wyckoff label, free variables, generated atom indices, site occupancy and
Debye-Waller parameters, space group, and status. ``wyckoff(site_id)`` returns
the same dictionary for one site.

Direct convenience setters use symmetric signatures for atom and site data:

.. code-block:: python

   # Surface-cell atom coordinates or per-atom physical parameters.
   unitcell.set_wyckoff_atom_parameter("O_4f", "occ", 0.5)
   unitcell.set_wyckoff_atom_parameter("O_4f", "z", oxygen_z_values)

   # Representative coordinate in the parent conventional cell.
   unitcell.set_wyckoff_site_parameter("O_4f", "x", 0.31)
   unitcell.set_wyckoff_site_parameter("O_4f", "iDW", 0.6)

The atom setter accepts ``x``, ``y``, ``z``, ``iDW``, ``oDW``, and ``occ``.
Coordinate values are in surface-cell fractional units and may be supplied per
atom; ``iDW``, ``oDW``, and ``occ`` are scalar site-wide values. The site setter
accepts ``x``, ``y``, and ``z`` in parent conventional fractional units and
propagates the displacement through the stored symmetry couplings. It also
accepts scalar site-wide ``iDW``, ``oDW``, and ``occ`` values.
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
       limits=(0.28, 0.34),
   )

The fitted value and limits are absolute Wyckoff-variable coordinates. The
implementation converts the absolute value to an internal delta from the
site's reference value and applies that delta through the affine couplings.
For rutile oxygen, coordinates such as ``u``, ``1-u``, ``0.5+u``, and
``0.5-u`` therefore move with the correct signs.

``addWyckoffShift`` fits a relative parent conventional fractional displacement
from the representative symmetry-site coordinate:

.. code-block:: python

   unitcell.addWyckoffShift(
       "O_4f",
       "x",
       limits=(-0.02, 0.02),
   )

Here ``"x"``, ``"y"``, and ``"z"`` are parent conventional fractional axes,
not generated surface-cell axes, and the parameter is a relative shift rather
than an absolute coordinate. The shift is propagated through the stored
space-group operation for each generated atom and then through the surface-cell
transform. This intentionally keeps the original atom assignment and atom
count while allowing the site to move away from its exact Wyckoff constraint.

Plural helpers are available for fitting several variables or shifts:

.. code-block:: python

   unitcell.addWyckoffParameters("Mg_8d")          # all free variables
   unitcell.addWyckoffShifts("Mg_8d", axes=("x", "z"))

Container Classes
-----------------

``Film`` and ``PoissonSurface`` forward Wyckoff helper calls to their single
template unit cell:

.. code-block:: python

   film.addWyckoffParameter("O_4f", "u", limits=(0.25, 0.35))
   poisson.addWyckoffShift("O_4f", "x", limits=(-0.02, 0.02))

``EpitaxyInterface`` follows the same selector convention as its existing
``addFitParameter`` API. Select the internal unit cell with ``unitcell``:

.. code-block:: python

   interface.addWyckoffParameter(
       "O_4f",
       "u",
       limits=(0.25, 0.35),
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
level. ``addWyckoffParameter`` distributes one absolute value to the selected
unit cells, while ``addWyckoffShift`` distributes one relative displacement.
If selected sites begin with different values, the coupled starting value is
their mean; setting the coupled parameter assigns the same absolute Wyckoff
variable to every selected site.

For a symmetry-preserving Wyckoff variable:

.. code-block:: python

   crystal.addWyckoffParameter(
       {
           "RuO2": ("O_4f", "u"),
           "TiO2": ("O_4f", "u"),
       },
       name="shared_rutile_oxygen_u",
       limits=(0.25, 0.35),
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
       limits=(0.25, 0.35),
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
* deduplicated matrices for Wyckoff-variable and representative-site
  couplings, referenced by the generated atoms.

The matrix representation avoids repeating one text row for every nonzero
atom/coordinate coupling. The reader also accepts the older expanded
``wyckoff_couplings`` and ``wyckoff_site_couplings`` tables.

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

.. autofunction:: orgui.datautils.xrayutils.CTRsymmetry.surface_terminations
   :no-index:

.. autofunction:: orgui.datautils.xrayutils.CTRsymmetry.structure_from_cif
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRsymmetry.SurfaceTermination
   :members: build_unitcell, outermost_atoms
   :no-index:

.. autoclass:: orgui.datautils.xrayutils.CTRuc.UnitCell
   :members: wyckoff_sites, wyckoff, set_wyckoff_atom_parameter,
             set_wyckoff_site_parameter, wyckoff_couplings, wyckoff_site_couplings,
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

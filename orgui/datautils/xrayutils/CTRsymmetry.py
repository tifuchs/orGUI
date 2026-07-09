# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
"""Symmetry metadata helpers for CTR surface unit cells.

The construction helpers in this module can use PyXtal to assign parent-cell
atoms to Wyckoff positions. The generated surface atoms are ordinary
``CTRuc.UnitCell`` atoms; symmetry information is kept as metadata so the
structure-factor calculation remains unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import itertools
import math

import numpy as np

from .HKLVlieg import Lattice


COORDINATE_NAMES = ("x", "y", "z")
DEFAULT_VARIABLE_NAMES = ("u", "v", "w")
SYMMETRY_SECTION_HEADERS = {
    "surface_transform:",
    "wyckoff_sites:",
    "wyckoff_atoms:",
    "wyckoff_couplings:",
    "wyckoff_site_couplings:",
}


@dataclass(frozen=True)
class AffineExpression:
    """Represent one fractional coordinate as an affine expression.

    :param float constant:
        Constant coordinate offset in fractional coordinates.
    :param dict coefficients:
        Mapping from variable name to fractional-coordinate coefficient.
    """

    constant: float
    coefficients: dict[str, float] = field(default_factory=dict)

    def evaluate(self, variables):
        """Evaluate the expression for a variable dictionary.

        :param dict variables:
            Mapping from variable name to current fractional-coordinate value.
        :returns:
            Coordinate value in fractional units.
        :rtype:
            float
        """
        value = self.constant
        for variable, factor in self.coefficients.items():
            value += factor * variables[variable]
        return value

    def shifted(self, shift):
        """Return a copy with an added constant offset.

        :param float shift:
            Offset in fractional coordinates.
        :returns:
            Shifted expression.
        :rtype:
            AffineExpression
        """
        return AffineExpression(self.constant + shift, dict(self.coefficients))


@dataclass(frozen=True)
class SurfaceCellSpec:
    """Describe an oriented surface cell derived from a parent cell.

    :param parent_a:
        Parent conventional-cell lengths in Angstrom.
    :param parent_alpha:
        Parent conventional-cell angles in degrees.
    :param transform:
        3-by-3 matrix whose columns are the surface-cell vectors expressed in
        parent conventional fractional coordinates.
    :param origin:
        Surface-cell origin in parent conventional fractional coordinates.
    :param z_bounds:
        Inclusive lower and exclusive upper bounds in surface fractional z.
    :param layer_origins:
        Optional fractional z positions used as layer starts. If provided,
        atoms are assigned to intervals between successive origins.
    :param translation_range:
        Parent-cell integer translations tested in each direction before the
        surface cell is wrapped and deduplicated.
    """

    parent_a: tuple[float, float, float]
    parent_alpha: tuple[float, float, float]
    transform: np.ndarray
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    z_bounds: tuple[float, float] = (0.0, 1.0)
    layer_origins: tuple[float, ...] | None = None
    translation_range: int | tuple[int, int] = 1

    def __post_init__(self):
        transform = np.asarray(self.transform, dtype=np.float64)
        if transform.shape != (3, 3):
            raise ValueError("Surface-cell transform must be a 3-by-3 matrix.")
        if abs(np.linalg.det(transform)) < 1e-12:
            raise ValueError("Surface-cell transform must be invertible.")
        object.__setattr__(self, "transform", transform)
        object.__setattr__(
            self,
            "parent_a",
            tuple(float(v) for v in self.parent_a),
        )
        object.__setattr__(
            self,
            "parent_alpha",
            tuple(float(v) for v in self.parent_alpha),
        )
        object.__setattr__(self, "origin", tuple(float(v) for v in self.origin))
        object.__setattr__(self, "z_bounds", tuple(float(v) for v in self.z_bounds))
        if self.layer_origins is not None:
            object.__setattr__(
                self,
                "layer_origins",
                tuple(float(v) % 1.0 for v in self.layer_origins),
            )

    @property
    def inverse_transform(self):
        """Return the parent-to-surface fractional-coordinate transform.

        :returns:
            Inverse of ``transform``.
        :rtype:
            numpy.ndarray
        """
        return np.linalg.inv(self.transform)

    def parent_translations(self):
        """Generate parent conventional-cell translations to search.

        :returns:
            Iterator over integer translation vectors in parent fractional
            coordinates.
        :rtype:
            iterator
        """
        if isinstance(self.translation_range, int):
            lo = -self.translation_range
            hi = self.translation_range
        else:
            lo, hi = self.translation_range
        for translation in itertools.product(range(lo, hi + 1), repeat=3):
            yield np.asarray(translation, dtype=np.float64)

    def surface_lattice_parameters(self):
        """Return surface-cell lengths and angles.

        :returns:
            Tuple ``(a, alpha)`` where ``a`` contains lengths in Angstrom and
            ``alpha`` contains angles in degrees.
        :rtype:
            tuple
        """
        parent_lattice = Lattice(self.parent_a, self.parent_alpha)
        surface_matrix = parent_lattice.R_mat @ self.transform
        lengths = np.linalg.norm(surface_matrix, axis=0)
        alpha = _angle_between(surface_matrix[:, 1], surface_matrix[:, 2])
        beta = _angle_between(surface_matrix[:, 2], surface_matrix[:, 0])
        gamma = _angle_between(surface_matrix[:, 0], surface_matrix[:, 1])
        return lengths, np.asarray([alpha, beta, gamma], dtype=np.float64)

    def parent_to_surface(self, parent_fractional):
        """Transform parent conventional fractional coordinates to surface coordinates.

        :param parent_fractional:
            Fractional coordinate in the parent conventional cell.
        :returns:
            Fractional coordinate in the surface cell.
        :rtype:
            numpy.ndarray
        """
        parent_fractional = np.asarray(parent_fractional, dtype=np.float64)
        return self.inverse_transform @ (parent_fractional - np.asarray(self.origin))


@dataclass(frozen=True)
class WyckoffSiteSpec:
    """Describe one affine Wyckoff site in the parent conventional cell.

    :param str site_id:
        Stable site identifier used by query and parameter APIs.
    :param str element:
        Element symbol for generated atoms.
    :param str wyckoff_label:
        Crystallographic Wyckoff label, for example ``"4f"``.
    :param tuple coordinates:
        Parent-cell coordinate expressions for all symmetry-generated atoms.
    :param dict variables:
        Current values of the independent Wyckoff variables.
    :param float occ:
        Site occupancy.
    :param float iDW:
        In-plane Debye-Waller parameter.
    :param float oDW:
        Out-of-plane Debye-Waller parameter.
    """

    site_id: str
    element: str
    wyckoff_label: str
    coordinates: tuple[tuple[AffineExpression, AffineExpression, AffineExpression], ...]
    representative_parent_fractional: tuple[float, float, float] | None = None
    operation_matrices: tuple[np.ndarray, ...] = ()
    variables: dict[str, float] = field(default_factory=dict)
    occ: float = 1.0
    iDW: float = 0.5
    oDW: float = 0.5

    def parent_positions(self):
        """Evaluate the parent conventional-cell positions.

        :returns:
            List of fractional coordinates in the parent conventional cell.
        :rtype:
            list
        """
        positions = []
        for coordinate in self.coordinates:
            positions.append(
                np.asarray(
                    [expression.evaluate(self.variables) for expression in coordinate],
                    dtype=np.float64,
                )
            )
        return positions


@dataclass(frozen=True)
class WyckoffCoupling:
    """Describe one atom-coordinate dependency on a Wyckoff variable.

    :param int atom_index:
        Index of the generated atom in ``UnitCell.basis``.
    :param str coordinate:
        Coordinate name, one of ``"x"``, ``"y"``, or ``"z"``.
    :param str variable:
        Wyckoff variable name.
    :param float constant:
        Constant part in surface fractional coordinates.
    :param float factor:
        Variable coefficient in surface fractional coordinates.
    :param str site_id:
        Site identifier owning the coupling.
    """

    atom_index: int
    coordinate: str
    variable: str
    constant: float
    factor: float
    site_id: str

    def value(self, variable_value):
        """Evaluate this coordinate for a variable value.

        :param float variable_value:
            Variable value in fractional coordinates.
        :returns:
            Surface fractional coordinate.
        :rtype:
            float
        """
        return self.constant + self.factor * variable_value


@dataclass(frozen=True)
class WyckoffSiteCoupling:
    """Describe one atom-coordinate dependency on site representative motion.

    :param int atom_index:
        Index of the generated atom in ``UnitCell.basis``.
    :param str coordinate:
        Surface coordinate name, one of ``"x"``, ``"y"``, or ``"z"``.
    :param str axis:
        Parent conventional representative coordinate, one of ``"x"``,
        ``"y"``, or ``"z"``.
    :param float factor:
        Surface-coordinate coefficient for a parent representative-coordinate
        displacement.
    :param str site_id:
        Site identifier owning the coupling.
    """

    atom_index: int
    coordinate: str
    axis: str
    factor: float
    site_id: str


@dataclass(frozen=True)
class GeneratedWyckoffAtom:
    """Metadata for one generated surface-cell atom.

    :param int atom_index:
        Index in ``UnitCell.basis``.
    :param str element:
        Element symbol.
    :param str site_id:
        Source Wyckoff site identifier.
    :param str wyckoff_label:
        Source Wyckoff label.
    :param numpy.ndarray parent_fractional:
        Parent conventional fractional coordinate including the selected
        parent-cell translation.
    :param numpy.ndarray surface_fractional:
        Wrapped surface-cell fractional coordinate.
    :param int layer:
        Assigned layer index stored in ``UnitCell.basis[:, 7]``.
    :param tuple couplings:
        Coordinate couplings for this atom.
    """

    atom_index: int
    element: str
    site_id: str
    wyckoff_label: str
    parent_fractional: np.ndarray
    surface_fractional: np.ndarray
    layer: int
    couplings: tuple[WyckoffCoupling, ...] = ()
    site_couplings: tuple[WyckoffSiteCoupling, ...] = ()


@dataclass
class SurfaceSymmetryModel:
    """Hold generated surface-cell atoms and their symmetry metadata."""

    surface_spec: SurfaceCellSpec
    sites: tuple[WyckoffSiteSpec, ...]
    atoms: list[GeneratedWyckoffAtom] = field(default_factory=list)
    spacegroup_number: int | None = None
    spacegroup_symbol: str | None = None

    def wyckoff_sites(self, parameters=None):
        """Return summary dictionaries for generated Wyckoff sites.

        :param dict parameters:
            Optional ``UnitCell.parameters`` dictionary used to report whether
            individual atom-coordinate parameters break a site's symmetry.
        :returns:
            List of site summary dictionaries.
        :rtype:
            list
        """
        output = []
        for site in self.sites:
            atoms = [
                atom.atom_index
                for atom in self.atoms
                if atom.site_id == site.site_id
            ]
            variables = dict(site.variables)
            output.append(
                {
                    "site_id": site.site_id,
                    "element": site.element,
                    "wyckoff_label": site.wyckoff_label,
                    "variables": variables,
                    "representative_parent_fractional": (
                        None
                        if site.representative_parent_fractional is None
                        else tuple(site.representative_parent_fractional)
                    ),
                    "atom_indices": atoms,
                    "spacegroup_number": self.spacegroup_number,
                    "spacegroup_symbol": self.spacegroup_symbol,
                    "status": self.symmetry_status(site.site_id, parameters),
                }
            )
        return output

    def wyckoff_couplings(self, site_id=None):
        """Return Wyckoff variable couplings for generated atom coordinates.

        :param str site_id:
            Optional site identifier. If omitted, all couplings are returned.
        :returns:
            List of :class:`WyckoffCoupling` objects.
        :rtype:
            list
        """
        couplings = []
        for atom in self.atoms:
            if site_id is None or atom.site_id == site_id:
                couplings.extend(atom.couplings)
        return couplings

    def wyckoff_site_couplings(self, site_id=None):
        """Return site-displacement couplings for generated atom coordinates.

        :param str site_id:
            Optional site identifier. If omitted, all couplings are returned.
        :returns:
            List of :class:`WyckoffSiteCoupling` objects.
        :rtype:
            list
        """
        couplings = []
        for atom in self.atoms:
            if site_id is None or atom.site_id == site_id:
                couplings.extend(atom.site_couplings)
        return couplings

    def atom_wyckoff_metadata(self, atom_index):
        """Return metadata for one generated atom.

        :param int atom_index:
            Index in ``UnitCell.basis``.
        :returns:
            Atom metadata or ``None`` when the atom has no symmetry metadata.
        :rtype:
            GeneratedWyckoffAtom or None
        """
        for atom in self.atoms:
            if atom.atom_index == atom_index:
                return atom
        return None

    def symmetry_status(self, site_id, parameters=None):
        """Return whether a site is still fully symmetry-preserving.

        :param str site_id:
            Site identifier.
        :param dict parameters:
            Optional ``UnitCell.parameters`` dictionary.
        :returns:
            ``"metadata_only"``, ``"symmetry_preserving"``, or
            ``"partially_overridden"``.
        :rtype:
            str
        """
        if not parameters:
            return "metadata_only"

        site_pairs = {
            (atom.atom_index, parameter_index)
            for atom in self.atoms
            if atom.site_id == site_id
            for parameter_index in (1, 2, 3)
        }
        if not site_pairs:
            return "metadata_only"

        preserving = False
        displaced = False
        for kind in ("absolute", "relative"):
            for parameter in parameters.get(kind, []):
                pairs = set(zip(parameter.indices[0], parameter.indices[1]))
                wyckoff_settings = parameter.settings.get("wyckoff", {})
                if (
                    kind == "relative"
                    and wyckoff_settings.get("site_id") == site_id
                    and pairs <= site_pairs
                ):
                    if wyckoff_settings.get("kind") == "coordinate":
                        preserving = True
                    elif wyckoff_settings.get("kind") == "site_displacement":
                        displaced = True
                    else:
                        return "partially_overridden"
                elif pairs & site_pairs:
                    return "partially_overridden"
        if displaced:
            return "site_displaced"
        if preserving:
            return "symmetry_preserving"
        return "metadata_only"

    def build_unitcell(self, name, layer_behavior="ignore"):
        """Create a ``CTRuc.UnitCell`` and attach this model as metadata.

        :param str name:
            Unit-cell name.
        :param str layer_behavior:
            Layer behavior passed to ``CTRuc.UnitCell``.
        :returns:
            Generated unit cell containing ordinary CTR atoms.
        :rtype:
            CTRuc.UnitCell
        """
        from .CTRuc import UnitCell

        a, alpha = self.surface_spec.surface_lattice_parameters()
        unitcell = UnitCell(a, alpha, name=name, layer_behavior=layer_behavior)
        generated = generate_surface_atoms(self.surface_spec, self.sites)
        for atom in generated:
            site = self._site(atom.site_id)
            unitcell.addAtom(
                atom.element,
                atom.surface_fractional,
                site.iDW,
                site.oDW,
                site.occ,
                layer=atom.layer,
            )
        if self.surface_spec.layer_origins is not None:
            unitcell.layerpos = {
                float(index + 1): origin
                for index, origin in enumerate(self.surface_spec.layer_origins)
            }
        self.atoms = [
            GeneratedWyckoffAtom(
                atom.atom_index,
                atom.element,
                atom.site_id,
                atom.wyckoff_label,
                atom.parent_fractional,
                atom.surface_fractional,
                atom.layer,
                atom.couplings,
                atom.site_couplings,
            )
            for atom in generated
        ]
        unitcell.symmetry_metadata = self
        return unitcell

    def _site(self, site_id):
        for site in self.sites:
            if site.site_id == site_id:
                return site
        raise KeyError(site_id)


def possible_wyckoff_positions(spacegroup, style="pyxtal"):
    """Return possible Wyckoff positions for a space group using PyXtal.

    PyXtal is imported lazily by this function.

    :param int spacegroup:
        International space-group number.
    :param str style:
        PyXtal Wyckoff setting style.
    :returns:
        List of dictionaries with ``label``, ``multiplicity``, and ``dof``.
    :rtype:
        list
    :raises ImportError:
        If PyXtal is not installed.
    """
    Group = _import_pyxtal_group()
    group = Group(spacegroup, style=style)
    positions = []
    for wp in group:
        label = wp.get_label()
        positions.append(
            {
                "label": label,
                "multiplicity": int(label[:-1]),
                "letter": label[-1],
                "dof": int(wp.get_dof()),
            }
        )
    return positions


def sites_from_seed(
    seed,
    tol=1e-4,
    a_tol=5.0,
    backend="pymatgen",
    style="pyxtal",
    iDW=0.5,
    oDW=0.5,
    occ=1.0,
    variable_names=DEFAULT_VARIABLE_NAMES,
):
    """Assign parent atoms to Wyckoff sites using PyXtal.

    PyXtal is imported lazily by this function. ``seed`` may be any input
    accepted by ``pyxtal.pyxtal().from_seed()``, including a CIF path,
    pymatgen ``Structure``, or ASE ``Atoms``.

    :param seed:
        Structure seed accepted by PyXtal.
    :param float tol:
        Coordinate tolerance passed to PyXtal.
    :param float a_tol:
        Angle tolerance in degrees passed to PyXtal.
    :param str backend:
        PyXtal seed backend.
    :param str style:
        PyXtal Wyckoff setting style.
    :param float iDW:
        In-plane Debye-Waller parameter assigned to generated sites.
    :param float oDW:
        Out-of-plane Debye-Waller parameter assigned to generated sites.
    :param float occ:
        Occupancy assigned to generated sites.
    :param tuple variable_names:
        Names for Wyckoff free variables in order of PyXtal free coordinates.
    :returns:
        Tuple ``(sites, spacegroup_number, spacegroup_symbol)``.
    :rtype:
        tuple
    :raises ImportError:
        If PyXtal is not installed.
    """
    pyxtal = _import_pyxtal_class()
    crystal = pyxtal()
    crystal.from_seed(seed, tol=tol, a_tol=a_tol, backend=backend, style=style)
    sites = tuple(
        _site_from_pyxtal_site(
            atom_site,
            variable_names=variable_names,
            iDW=iDW,
            oDW=oDW,
            occ=occ,
        )
        for atom_site in crystal.atom_sites
    )
    return sites, int(crystal.group.number), crystal.group.symbol


def model_from_seed(
    seed,
    surface_spec,
    tol=1e-4,
    a_tol=5.0,
    backend="pymatgen",
    style="pyxtal",
    iDW=0.5,
    oDW=0.5,
    occ=1.0,
    variable_names=DEFAULT_VARIABLE_NAMES,
):
    """Build a surface symmetry model from a PyXtal-compatible seed.

    :param seed:
        Structure seed accepted by PyXtal.
    :param SurfaceCellSpec surface_spec:
        Surface-cell transform and parent lattice.
    :returns:
        Surface symmetry model ready to build a ``UnitCell``.
    :rtype:
        SurfaceSymmetryModel
    """
    sites, spacegroup_number, spacegroup_symbol = sites_from_seed(
        seed,
        tol=tol,
        a_tol=a_tol,
        backend=backend,
        style=style,
        iDW=iDW,
        oDW=oDW,
        occ=occ,
        variable_names=variable_names,
    )
    return SurfaceSymmetryModel(
        surface_spec,
        sites,
        spacegroup_number=spacegroup_number,
        spacegroup_symbol=spacegroup_symbol,
    )


def surface_unitcell_from_seed(
    seed,
    surface_spec,
    name,
    layer_behavior="ignore",
    **kwargs,
):
    """Build a generated CTR surface ``UnitCell`` from a PyXtal-compatible seed.

    :param seed:
        Structure seed accepted by PyXtal.
    :param SurfaceCellSpec surface_spec:
        Surface-cell transform and parent lattice.
    :param str name:
        Unit-cell name.
    :param str layer_behavior:
        Layer behavior passed to ``CTRuc.UnitCell``.
    :param kwargs:
        Additional keyword arguments passed to :func:`model_from_seed`.
    :returns:
        Generated CTR unit cell with attached symmetry metadata.
    :rtype:
        CTRuc.UnitCell
    """
    model = model_from_seed(seed, surface_spec, **kwargs)
    return model.build_unitcell(name, layer_behavior=layer_behavior)


def rutile_110_surface_spec(parent_a, parent_alpha=(90.0, 90.0, 90.0)):
    """Return the standard rutile ``(110)`` commensurate surface-cell spec.

    The surface basis columns are ``[1 -1 0]``, ``[0 0 1]``, and ``[1 1 0]``
    in parent conventional fractional coordinates.

    :param parent_a:
        Parent conventional-cell lengths in Angstrom.
    :param parent_alpha:
        Parent conventional-cell angles in degrees.
    :returns:
        Surface-cell specification.
    :rtype:
        SurfaceCellSpec
    """
    transform = np.asarray(
        [
            [1.0, 0.0, 1.0],
            [-1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    return SurfaceCellSpec(
        tuple(parent_a),
        tuple(parent_alpha),
        transform,
        layer_origins=(0.0, 0.5),
        translation_range=1,
    )


def generate_surface_atoms(surface_spec, sites, tolerance=1e-8):
    """Generate surface-cell atom metadata from parent Wyckoff sites.

    :param SurfaceCellSpec surface_spec:
        Surface-cell transform and bounds.
    :param tuple sites:
        Parent-cell Wyckoff site specifications.
    :param float tolerance:
        Fractional-coordinate tolerance for bounds and deduplication.
    :returns:
        Generated atom metadata. ``atom_index`` values match the returned list
        order and therefore the order used by ``SurfaceSymmetryModel``.
    :rtype:
        list
    """
    generated = []
    seen = set()
    inverse = surface_spec.inverse_transform
    origin = np.asarray(surface_spec.origin, dtype=np.float64)
    z_min, z_max = surface_spec.z_bounds

    for site in sites:
        for coordinate_index, coordinate in enumerate(site.coordinates):
            for translation in surface_spec.parent_translations():
                parent_exprs = tuple(
                    expression.shifted(translation[index])
                    for index, expression in enumerate(coordinate)
                )
                transform_exprs = tuple(
                    expression.shifted(translation[index] - origin[index])
                    for index, expression in enumerate(coordinate)
                )
                surface_exprs = _transform_expressions(inverse, transform_exprs)
                surface_values = np.asarray(
                    [
                        expression.evaluate(site.variables)
                        for expression in surface_exprs
                    ],
                    dtype=np.float64,
                )
                if surface_values[2] < z_min - tolerance:
                    continue
                if surface_values[2] >= z_max - tolerance:
                    continue

                wrapped_values, wrapped_exprs = _wrap_surface_expressions(
                    surface_values,
                    surface_exprs,
                )
                key = (
                    site.site_id,
                    tuple(np.round(wrapped_values / tolerance).astype(np.int64)),
                )
                if key in seen:
                    continue
                seen.add(key)

                layer = _assign_layer(wrapped_values[2], surface_spec.layer_origins)
                atom_index = len(generated)
                couplings = []
                for axis, expression in enumerate(wrapped_exprs):
                    for variable, factor in expression.coefficients.items():
                        if not math.isclose(factor, 0.0, abs_tol=tolerance):
                            couplings.append(
                                WyckoffCoupling(
                                    atom_index=atom_index,
                                    coordinate=COORDINATE_NAMES[axis],
                                    variable=variable,
                                    constant=expression.constant,
                                    factor=factor,
                                    site_id=site.site_id,
                                )
                            )
                site_couplings = []
                if coordinate_index < len(site.operation_matrices):
                    operation = np.asarray(
                        site.operation_matrices[coordinate_index],
                        dtype=np.float64,
                    )
                    site_displacements = inverse @ operation[:3, :3]
                    for parent_axis, axis_name in enumerate(COORDINATE_NAMES):
                        for surface_axis, coordinate_name in enumerate(
                            COORDINATE_NAMES
                        ):
                            factor = site_displacements[surface_axis, parent_axis]
                            if not math.isclose(factor, 0.0, abs_tol=tolerance):
                                site_couplings.append(
                                    WyckoffSiteCoupling(
                                        atom_index=atom_index,
                                        coordinate=coordinate_name,
                                        axis=axis_name,
                                        factor=float(factor),
                                        site_id=site.site_id,
                                    )
                                )
                generated.append(
                    GeneratedWyckoffAtom(
                        atom_index=atom_index,
                        element=site.element,
                        site_id=site.site_id,
                        wyckoff_label=site.wyckoff_label,
                        parent_fractional=np.asarray(
                            [
                                expression.evaluate(site.variables)
                                for expression in parent_exprs
                            ],
                            dtype=np.float64,
                        ),
                        surface_fractional=wrapped_values,
                        layer=layer,
                        couplings=tuple(couplings),
                        site_couplings=tuple(site_couplings),
                    )
                )

    generated.sort(
        key=lambda atom: (
            -atom.layer,
            atom.site_id,
            round(atom.surface_fractional[2], 10),
            round(atom.surface_fractional[1], 10),
            round(atom.surface_fractional[0], 10),
        )
    )
    return [
        GeneratedWyckoffAtom(
            atom_index=index,
            element=atom.element,
            site_id=atom.site_id,
            wyckoff_label=atom.wyckoff_label,
            parent_fractional=atom.parent_fractional,
            surface_fractional=atom.surface_fractional,
            layer=atom.layer,
            couplings=tuple(
                WyckoffCoupling(
                    atom_index=index,
                    coordinate=coupling.coordinate,
                    variable=coupling.variable,
                    constant=coupling.constant,
                    factor=coupling.factor,
                    site_id=coupling.site_id,
                )
                for coupling in atom.couplings
            ),
            site_couplings=tuple(
                WyckoffSiteCoupling(
                    atom_index=index,
                    coordinate=coupling.coordinate,
                    axis=coupling.axis,
                    factor=coupling.factor,
                    site_id=coupling.site_id,
                )
                for coupling in atom.site_couplings
            ),
        )
        for index, atom in enumerate(generated)
    ]


def symmetry_metadata_to_lines(model):
    """Serialize symmetry metadata as editable plain-text table lines.

    :param SurfaceSymmetryModel model:
        Symmetry model to serialize.
    :returns:
        Lines without trailing newline characters.
    :rtype:
        list
    """
    lines = []
    if model.spacegroup_number is not None:
        symbol = model.spacegroup_symbol or ""
        lines.append(f"spacegroup: {model.spacegroup_number} {symbol}".rstrip())
    lines.append("surface_transform:")
    for row in model.surface_spec.transform:
        lines.append("  " + _format_values(row))
    lines.append("surface_origin: " + _format_values(model.surface_spec.origin))
    lines.append("wyckoff_sites:")
    lines.append(
        "  site_id element wyckoff_label variables "
        "representative_x representative_y representative_z occ iDW oDW"
    )
    for site in model.sites:
        variables = _format_variables(site.variables)
        representative = site.representative_parent_fractional
        representative_values = (
            ["-", "-", "-"]
            if representative is None
            else _format_values(representative).split()
        )
        lines.append(
            "  "
            + " ".join(
                [
                    site.site_id,
                    site.element,
                    site.wyckoff_label,
                    variables,
                    *representative_values,
                    _format_float(site.occ),
                    _format_float(site.iDW),
                    _format_float(site.oDW),
                ]
            )
        )
    lines.append("wyckoff_atoms:")
    lines.append(
        "  atom_index element site_id wyckoff_label "
        "parent_x parent_y parent_z surface_x surface_y surface_z layer"
    )
    for atom in model.atoms:
        values = [
            str(atom.atom_index),
            atom.element,
            atom.site_id,
            atom.wyckoff_label,
            *_format_values(atom.parent_fractional).split(),
            *_format_values(atom.surface_fractional).split(),
            _format_float(atom.layer),
        ]
        lines.append("  " + " ".join(values))
    lines.append("wyckoff_couplings:")
    lines.append("  atom_index coordinate site_id variable constant factor")
    for coupling in model.wyckoff_couplings():
        lines.append(
            "  "
            + " ".join(
                [
                    str(coupling.atom_index),
                    coupling.coordinate,
                    coupling.site_id,
                    coupling.variable,
                    _format_float(coupling.constant),
                    _format_float(coupling.factor),
                ]
            )
        )
    lines.append("wyckoff_site_couplings:")
    lines.append("  atom_index coordinate site_id axis factor")
    for coupling in model.wyckoff_site_couplings():
        lines.append(
            "  "
            + " ".join(
                [
                    str(coupling.atom_index),
                    coupling.coordinate,
                    coupling.site_id,
                    coupling.axis,
                    _format_float(coupling.factor),
                ]
            )
        )
    return lines


def symmetry_metadata_from_lines(lines, unitcell=None):
    """Deserialize plain-text table lines into a symmetry metadata model.

    This parser has no PyXtal dependency and is used when loading `.xtal` and
    `.xpr` files containing resolved symmetry metadata.

    :param list lines:
        Symmetry metadata lines.
    :param CTRuc.UnitCell unitcell:
        Optional unit cell providing fallback lattice and atom values.
    :returns:
        Parsed symmetry model, or ``None`` when no symmetry data is present.
    :rtype:
        SurfaceSymmetryModel or None
    """
    cleaned = [line.strip() for line in lines if line.strip()]
    if not cleaned:
        return None

    spacegroup_number = None
    spacegroup_symbol = None
    transform = np.identity(3, dtype=np.float64)
    origin = (0.0, 0.0, 0.0)
    site_specs = []
    atoms = []
    couplings_by_atom = {}
    site_couplings_by_atom = {}
    sections = _collect_sections(cleaned)

    if "spacegroup" in sections:
        parts = sections["spacegroup"][0].split()
        if parts:
            spacegroup_number = int(parts[0])
            if len(parts) > 1:
                spacegroup_symbol = " ".join(parts[1:])

    if "surface_transform" in sections:
        transform = np.asarray(
            [
                [float(value) for value in row.split()]
                for row in sections["surface_transform"]
            ],
            dtype=np.float64,
        )
    if "surface_origin" in sections:
        origin = tuple(float(value) for value in sections["surface_origin"][0].split())

    for row in _data_rows(sections.get("wyckoff_sites", [])):
        parts = row.split()
        variables = _parse_variables(parts[3])
        if len(parts) >= 10:
            representative = (
                None
                if parts[4] == "-"
                else tuple(float(value) for value in parts[4:7])
            )
            occ, iDW, oDW = (float(value) for value in parts[7:10])
        else:
            representative = None
            occ, iDW, oDW = (float(value) for value in parts[4:7])
        site_specs.append(
            WyckoffSiteSpec(
                site_id=parts[0],
                element=parts[1],
                wyckoff_label=parts[2],
                coordinates=(),
                representative_parent_fractional=representative,
                variables=variables,
                occ=occ,
                iDW=iDW,
                oDW=oDW,
            )
        )

    for row in _data_rows(sections.get("wyckoff_couplings", [])):
        parts = row.split()
        coupling = WyckoffCoupling(
            atom_index=int(parts[0]),
            coordinate=parts[1],
            site_id=parts[2],
            variable=parts[3],
            constant=float(parts[4]),
            factor=float(parts[5]),
        )
        couplings_by_atom.setdefault(coupling.atom_index, []).append(coupling)

    for row in _data_rows(sections.get("wyckoff_site_couplings", [])):
        parts = row.split()
        coupling = WyckoffSiteCoupling(
            atom_index=int(parts[0]),
            coordinate=parts[1],
            site_id=parts[2],
            axis=parts[3],
            factor=float(parts[4]),
        )
        site_couplings_by_atom.setdefault(coupling.atom_index, []).append(coupling)

    for row in _data_rows(sections.get("wyckoff_atoms", [])):
        parts = row.split()
        atom_index = int(parts[0])
        atoms.append(
            GeneratedWyckoffAtom(
                atom_index=atom_index,
                element=parts[1],
                site_id=parts[2],
                wyckoff_label=parts[3],
                parent_fractional=np.asarray(parts[4:7], dtype=np.float64),
                surface_fractional=np.asarray(parts[7:10], dtype=np.float64),
                layer=int(float(parts[10])),
                couplings=tuple(couplings_by_atom.get(atom_index, ())),
                site_couplings=tuple(site_couplings_by_atom.get(atom_index, ())),
            )
        )

    parent_a = (1.0, 1.0, 1.0)
    parent_alpha = (90.0, 90.0, 90.0)
    if unitcell is not None:
        parent_a = tuple(float(value) for value in unitcell.a)
        parent_alpha = tuple(float(value) for value in np.rad2deg(unitcell.alpha))
    surface_spec = SurfaceCellSpec(
        parent_a,
        parent_alpha,
        transform,
        origin=origin,
        layer_origins=_layer_origins_from_unitcell(unitcell),
    )
    return SurfaceSymmetryModel(
        surface_spec,
        tuple(site_specs),
        atoms=atoms,
        spacegroup_number=spacegroup_number,
        spacegroup_symbol=spacegroup_symbol,
    )


def _import_pyxtal_class():
    try:
        from pyxtal import pyxtal
    except ImportError as exc:
        raise ImportError(
            "PyXtal is required for symmetry construction. Install orGUI with "
            "the optional 'symmetry' dependencies to use this feature."
        ) from exc
    return pyxtal


def _import_pyxtal_group():
    try:
        from pyxtal.symmetry import Group
    except ImportError as exc:
        raise ImportError(
            "PyXtal is required for Wyckoff lookup. Install orGUI with the "
            "optional 'symmetry' dependencies to use this feature."
        ) from exc
    return Group


def _site_from_pyxtal_site(atom_site, variable_names, iDW, oDW, occ):
    wp = atom_site.wp
    representative = _wrap_fractional(atom_site.position)
    free_values = list(wp.get_free_xyzs(atom_site.position))
    if len(free_values) > len(variable_names):
        raise ValueError(
            f"Not enough Wyckoff variable names for {wp.get_label()}: "
            f"{len(free_values)} required."
        )
    variables = {
        variable_names[index]: float(value)
        for index, value in enumerate(free_values)
    }
    base_position, coefficients = _primary_position_affine(
        wp,
        free_values,
        tuple(variables),
    )
    primary = tuple(
        AffineExpression(base_position[axis], dict(coefficients[axis]))
        for axis in range(3)
    )
    coordinates = tuple(
        _apply_operation_to_expressions(operation.affine_matrix, primary)
        for operation in wp.ops
    )
    operation_matrices = _site_operation_matrices(
        wp,
        representative,
        coordinates,
        variables,
    )
    element = str(atom_site.specie)
    return WyckoffSiteSpec(
        site_id=f"{element}_{wp.get_label()}",
        element=element,
        wyckoff_label=wp.get_label(),
        coordinates=coordinates,
        representative_parent_fractional=tuple(
            float(value) for value in representative
        ),
        operation_matrices=operation_matrices,
        variables=variables,
        occ=occ,
        iDW=iDW,
        oDW=oDW,
    )


def _site_operation_matrices(wp, representative, coordinates, variables):
    group = _import_pyxtal_group()(wp.number)
    general_operations = group[0].ops
    operation_matrices = []
    for coordinate in coordinates:
        target = _wrap_fractional(
            [expression.evaluate(variables) for expression in coordinate]
        )
        for operation in general_operations:
            matrix = np.asarray(operation.affine_matrix, dtype=np.float64)
            generated = _wrap_fractional(
                matrix[:3, :3] @ representative + matrix[:3, 3]
            )
            if _fractional_positions_close(generated, target):
                operation_matrices.append(matrix)
                break
        else:
            raise ValueError(
                f"Could not match a full symmetry operation for Wyckoff "
                f"position {wp.get_label()} at {target}."
            )
    return tuple(operation_matrices)


def _primary_position_affine(wp, free_values, variable_names):
    if not free_values:
        return np.asarray(wp.get_position_from_free_xyzs([]), dtype=np.float64), (
            {},
            {},
            {},
        )
    zeros = np.zeros(len(free_values), dtype=np.float64)
    base_position = np.asarray(wp.get_position_from_free_xyzs(zeros), dtype=np.float64)
    coefficients = []
    for axis in range(3):
        coefficients.append({})
    step = 1e-5
    for index, variable in enumerate(variable_names):
        probe = np.zeros(len(free_values), dtype=np.float64)
        probe[index] = step
        position = np.asarray(wp.get_position_from_free_xyzs(probe), dtype=np.float64)
        delta = (position - base_position) / step
        for axis, factor in enumerate(delta):
            if not math.isclose(factor, 0.0, abs_tol=1e-12):
                coefficients[axis][variable] = float(factor)
    return base_position, tuple(coefficients)


def _apply_operation_to_expressions(matrix, expressions):
    matrix = np.asarray(matrix, dtype=np.float64)
    transformed = []
    for row in matrix[:3]:
        constant = row[3]
        coefficients = {}
        for factor, expression in zip(row[:3], expressions):
            constant += factor * expression.constant
            for variable, coefficient in expression.coefficients.items():
                coefficients[variable] = coefficients.get(variable, 0.0) + (
                    factor * coefficient
                )
        transformed.append(AffineExpression(float(constant), coefficients))
    return tuple(transformed)


def _angle_between(v1, v2):
    cosine = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _transform_expressions(matrix, expressions):
    transformed = []
    for row in matrix:
        constant = 0.0
        coefficients = {}
        for factor, expression in zip(row, expressions):
            constant += factor * expression.constant
            for variable, coefficient in expression.coefficients.items():
                coefficients[variable] = coefficients.get(variable, 0.0) + (
                    factor * coefficient
                )
        transformed.append(AffineExpression(constant, coefficients))
    return tuple(transformed)


def _wrap_surface_expressions(values, expressions):
    wrapped_values = np.mod(values, 1.0)
    wrapped_values[np.isclose(wrapped_values, 1.0, atol=1e-12)] = 0.0
    wrapped_expressions = []
    for value, wrapped_value, expression in zip(values, wrapped_values, expressions):
        shift = wrapped_value - value
        wrapped_expressions.append(expression.shifted(round(shift)))
    return wrapped_values, tuple(wrapped_expressions)


def _wrap_fractional(values):
    wrapped = np.mod(np.asarray(values, dtype=np.float64), 1.0)
    wrapped[np.isclose(wrapped, 1.0, atol=1e-12)] = 0.0
    return wrapped


def _fractional_positions_close(first, second, tolerance=1e-6):
    diff = np.mod(np.asarray(first) - np.asarray(second) + 0.5, 1.0) - 0.5
    return np.allclose(diff, 0.0, atol=tolerance)


def _assign_layer(z_value, layer_origins):
    if not layer_origins:
        return 0
    origins = np.asarray(sorted(layer_origins), dtype=np.float64)
    z_value = z_value % 1.0
    matches = np.flatnonzero(z_value + 1e-12 >= origins)
    if matches.size == 0:
        return len(origins)
    return int(matches[-1] + 1)


def _format_float(value):
    return f"{float(value):.12g}"


def _format_values(values):
    return " ".join(_format_float(value) for value in values)


def _format_variables(variables):
    if not variables:
        return "-"
    return ",".join(
        f"{name}={_format_float(value)}" for name, value in variables.items()
    )


def _parse_variables(text):
    if text == "-":
        return {}
    variables = {}
    for item in text.split(","):
        name, value = item.split("=", 1)
        variables[name] = float(value)
    return variables


def _collect_sections(lines):
    sections = {}
    current = None
    for line in lines:
        if line.startswith("spacegroup:"):
            sections["spacegroup"] = [line.split(":", 1)[1].strip()]
            current = None
        elif line.startswith("surface_origin:"):
            sections["surface_origin"] = [line.split(":", 1)[1].strip()]
            current = None
        elif line in SYMMETRY_SECTION_HEADERS:
            current = line[:-1]
            sections[current] = []
        elif current is not None:
            sections[current].append(line)
    return sections


def _data_rows(rows):
    for row in rows:
        stripped = row.strip()
        if not stripped:
            continue
        if stripped.split()[0] in {
            "site_id",
            "atom_index",
        }:
            continue
        yield stripped


def _layer_origins_from_unitcell(unitcell):
    if unitcell is None or not getattr(unitcell, "layerpos", None):
        return None
    return tuple(unitcell.layerpos[layer] for layer in sorted(unitcell.layerpos))

"""Material-resolved finite-film dissolution using native crystal states.

The interface closure is frozen and material selective. A structural front
removes finite-owner populations; it never replenishes the protected substrate.
All amplitudes and optical profiles are evaluated by ordinary native crystals.
"""

import copy
from dataclasses import dataclass

import numpy as np

from .CTRcalc import SXRDCrystal
from .CTRdistributions import (
    DEFAULT_TAIL_PROBABILITY, PoissonProfile, RecenteredSkellamProfile,
    SkellamProfile,
)
from .CTRfilm import (
    EpitaxyInterface, Film, RoughSurface, SkellamSurface, _FrozenMaterial,
)


@dataclass(frozen=True)
class MaterialRecord:
    """One native domain with owner and algebraic-reference provenance.

    ``lower_A`` and ``upper_A`` are structural boundaries, not atom extrema.
    Protected reference records may have negative weights; removable records
    are combined with their exact replacements before state selection.
    """

    owner: str | None
    role: str
    cell: object
    matrix: np.ndarray
    weight: float
    lower_A: float
    upper_A: float
    interval: tuple


@dataclass(frozen=True)
class DissolutionState:
    """Complete native state and species inventories per reference area."""

    index: int
    front_height_A: float
    exposed_owners: tuple
    inventory: dict
    crystal: SXRDCrystal


def _geometry_key(record):
    cell = record.cell
    return (record.owner, tuple(cell.names), cell.basis.tobytes(),
            cell.a.tobytes(), cell.alpha.tobytes(),
            record.matrix.tobytes(), cell.refHKLTransform.tobytes(),
            repr(getattr(cell, "_optical_plane_geometry", None)),
            tuple(cell.getLayerCrystallinity(x)
                  for x in cell.crystallinity_layers))


def _inventory(records, reference_area):
    result = {}
    for record in records:
        if record.owner is None:
            continue
        owner = result.setdefault(record.owner, {})
        scale = record.weight * reference_area/record.cell.uc_area
        for name, occupancy in zip(record.cell.names, record.cell.basis[:, 6]):
            owner[name] = owner.get(name, 0.) + float(scale*occupancy)
    return result


def _copy_geometry(source, target):
    target.dwba.set_orientation(source.dwba.orientation)
    target.dwba._default_geometry = copy.deepcopy(source.dwba._default_geometry)
    target.dwba._rod_geometry = copy.deepcopy(source.dwba._rod_geometry)
    target.dwba.set_field_cache_capacity(source.dwba._field_cache_capacity)


class DissolutionPath:
    """Snapshot a connected finite-material path above protected bulk.

    :param SXRDCrystal crystal: Native source; never mutated.
    :param films: Stable names of a connected upper sequence of finite Films.
    :param str surface: Name of the top RoughSurface carrying native motifs.
    :param str exposure_policy:
        ``preserve`` retains unmodified survivors. ``native_or_unreconstructed``
        uses the native top-Film motif only when its complete depth is present;
        otherwise it explicitly exposes the unmodified survivors.
    :param dict interface_owners:
        Optional persisted mapping from interface names to [lower, upper]
        owner names. The protected bulk is represented by None.
    :raises ValueError: For ambiguous topology, geometry, or ownership.
    """

    def __init__(self, crystal, *, films, surface,
                 exposure_policy="native_or_unreconstructed",
                 interface_owners=None):
        if exposure_policy not in {"preserve", "native_or_unreconstructed"}:
            raise ValueError("Unknown exposure_policy")
        self.films = tuple(films)
        if not self.films or len(set(self.films)) != len(self.films):
            raise ValueError("films must contain distinct finite Film names")
        self.exposure_policy = exposure_policy
        self.source = copy.deepcopy(crystal, {id(crystal.dwba): None})
        _copy_geometry(crystal, self.source)
        components = self.source.uc_surface_list
        names = [c.name for c in components]
        if len(set(names)) != len(names):
            raise ValueError("Component names must be unique")
        matches = [c for c in components if c.name == surface]
        if len(matches) != 1 or not isinstance(matches[0], RoughSurface):
            raise ValueError("surface must identify one top RoughSurface")
        self.surface = matches[0]
        if any(self.surface.parameters.values()):
            raise ValueError("Fit parent height parameters on the dissolution wrapper")
        if self.surface.water is not None:
            raise ValueError("Dissolution water placement is not defined")
        if self.surface.crystallinity_profile is not None:
            raise ValueError("Dissolution depth-dependent surface order is not defined")
        # The distribution belongs to the ensemble. Lay out a deterministic
        # nominal source without asking the old surface to etch below its Film.
        self.surface.profile = (RecenteredSkellamProfile()
                                if isinstance(self.surface, SkellamSurface)
                                else PoissonProfile(0))
        self.surface.basis = np.array([
            getattr(self.surface.profile, key) for key in self.surface._basis_fields])
        self.surface.basis_0 = self.surface.basis.copy()
        order = np.argsort(self.source.uc_stacking)
        if len(set(self.source.uc_stacking)) != len(components):
            raise ValueError("Dissolution requires an ordered single stack")
        ordered = [components[i] for i in order]
        if ordered[-1] is not self.surface:
            raise ValueError("Overlayers above the dissolution surface are unsupported")
        if any(not isinstance(c, (Film, EpitaxyInterface, RoughSurface))
               for c in ordered):
            raise ValueError("Dissolution requires Film/interface components")
        if any(getattr(c, "crystallinity_profiles", {}) for c in ordered):
            raise ValueError(
                "Dissolution depth-dependent interface order is not defined")
        available = [c.name for c in ordered if isinstance(c, Film)]
        if any(name not in available for name in self.films):
            raise ValueError(
                "Only explicit finite Films may dissolve; bulk is protected")
        selected = [name for name in available if name in self.films]
        if selected != available[-len(selected):]:
            raise ValueError("Removable Films must form a connected upper sequence")
        self.films = tuple(selected)
        if not np.all(self.source.weights == 1):
            raise ValueError("Dissolution requires unit component weights")
        self.domains = copy.deepcopy(self.source.domains[order[-1]])
        for domains in self.source.domains:
            if (len(domains) != len(self.domains)
                    or any(w != v or not np.array_equal(a, b)
                           for (a, w), (b, v) in zip(domains, self.domains))):
                raise ValueError("Dissolution requires common outer domains")
        self.source.apply_stacking()
        self.top_film = self.surface.underlying_film
        if self.top_film is None or self.top_film.name != self.films[-1]:
            raise ValueError("Surface must be bound to the highest removable Film")
        self.interface_owners = {}
        self.interface_support_audit = {}
        last_film = None
        for i, component in enumerate(ordered[:-1]):
            if isinstance(component, Film):
                last_film = component.name
            else:
                if i+1 >= len(ordered) or not isinstance(ordered[i+1], Film):
                    raise ValueError(
                        "Each interface must be followed by its upper Film")
                self.interface_owners[component.name] = [last_film, ordered[i+1].name]
                count = len(component.layer_order)
                interface_profile = SkellamProfile(
                    component.basis[0], component.basis[1],
                    getattr(component.profile, "tail_probability",
                            DEFAULT_TAIL_PROBABILITY))
                up, down = interface_profile.parameters(count)
                parent = RecenteredSkellamProfile(1, up, down, up-down)
                offsets = np.concatenate([
                    cell._crystallinity_domain_offsets
                    for cell in component.top_layers])
                low, high = int(offsets.min()), int(offsets.max())
                probability, moment = parent.omitted_bounds(
                    min(low, 0), max(high, 0))
                # Sum of missing cumulative occupancies is bounded by
                # E[|X|;outside] + max(|edge|)*P(outside).
                self.interface_support_audit[component.name] = dict(
                    offsets=[low, high], omitted_probability_bound=probability,
                    omitted_absolute_first_moment_bound=moment,
                    missing_plane_inventory_bound=(
                        moment + max(abs(low), abs(high)+1)*probability))
        if interface_owners is not None and interface_owners != self.interface_owners:
            raise ValueError("Persisted interface ownership disagrees with the stack")
        interfaces = [c for c in ordered if isinstance(c, EpitaxyInterface)]
        for lower, upper in zip(interfaces, interfaces[1:]):
            if lower._support_bounds[1] > upper._support_bounds[0] + 1e-9:
                raise ValueError("Overlapping interface support is unsupported")
        records = []
        for component in ordered[:-1]:
            if isinstance(component, Film):
                for cell in component.layer_ucs:
                    records.extend(self._domains(
                        cell, component.unitcell, component.name, "film",
                        component.name))
            else:
                lower_owner, upper_owner = self.interface_owners[component.name]
                for side, cells, parent, owner in (
                    ("top", component.top_layers, component.uc_top, upper_owner),
                    ("bottom", component.bottom_layers, component.uc_bottom,
                     lower_owner),
                ):
                    for cell in cells:
                        records.extend(self._domains(
                            cell, parent, owner, side, component.name,
                            component.layer_order))
                for cell in component._reference_layers:
                    records.extend(self._domains(
                        cell, cell, None, "reference", component.name))
        # Exact native replacement records cancel their original populations.
        # No negative atom population may be removed as if it were material.
        combined = {}
        protected = []
        for record in records:
            if record.owner not in self.films:
                protected.append(record)
                continue
            key = _geometry_key(record)
            if key in combined:
                old = combined[key]
                combined[key] = MaterialRecord(
                    old.owner, old.role, old.cell, old.matrix,
                    old.weight+record.weight, old.lower_A, old.upper_A,
                    old.interval)
            else:
                combined[key] = record
        removable = list(combined.values())
        if any(r.weight < -1e-12 for r in removable):
            raise ValueError("Unmatched signed reference in removable inventory")
        removable = [r for r in removable if r.weight > 0]
        if not removable:
            raise ValueError("Selected Films have no represented material")
        intervals = {}
        for r in removable:
            low, high = intervals.get(r.interval, (np.inf, -np.inf))
            intervals[r.interval] = min(low, r.lower_A), max(high, r.upper_A)
        self._intervals = sorted(intervals, key=lambda key: intervals[key][0])
        self._rank = {key: i for i, key in enumerate(self._intervals)}
        self.minimum = -len(self._intervals)
        self.nominal_height_A = float(self.surface.below_H)
        self._heights = [intervals[self._intervals[0]][0]]
        self._heights.extend(intervals[key][1] for key in self._intervals)
        if np.any(np.diff(self._heights) < -1e-9):
            raise ValueError("Interface front heights are not ordered")
        self._heights[-1] = self.nominal_height_A
        self.protected_owner = (available[-len(self.films)-1]
                                if len(available) > len(self.films) else "__bulk__")
        self.records = tuple(protected + removable)
        self._protected = tuple(protected)
        self._removable = tuple(removable)
        self._interval_records = {
            key: tuple(r for r in removable if r.interval == key)
            for key in self._intervals}
        self.initial_inventory = self.state(0).inventory

    def _domains(self, cell, parent, owner, role, component_name, order=None):
        layer = float(cell.layers[0])
        sequence = list(parent.layer_cycle.layers)
        positions = parent.layerpos
        i = sequence.index(layer)
        spacing = (positions[sequence[(i+1) % len(sequence)]] - positions[layer]) % 1
        if spacing == 0:
            spacing = 1.
        provenance = getattr(cell, "_dissolution_reference_owners", [])
        for j, (matrix, weight) in enumerate(zip(
                cell.coherentDomainMatrix, cell.coherentDomainOccupancy)):
            actual_owner = owner
            reference = provenance[j] if j < len(provenance) else None
            if reference is not None:
                kind, name = reference
                if kind == "film" and name not in {
                        c.name for c in self.source.uc_surface_list
                        if isinstance(c, Film)}:
                    raise ValueError("Nested interface replacement is unsupported")
                actual_owner = name if kind == "film" else None
            low = (matrix[2, 2]*positions[layer]+matrix[2, 3])*cell.a[2]
            high = low + matrix[2, 2]*spacing*cell.a[2]
            if high <= low:
                raise ValueError("Structural intervals must have positive height")
            # The paired interface sides share explicit cycle/index identity,
            # even when their frozen geometries have different displacements.
            interval = ((component_name, j, list(order).index(layer))
                        if order is not None and reference is None else
                        (component_name, float(low)))
            mat = np.array(matrix, copy=True)
            mat.flags.writeable = False
            yield MaterialRecord(actual_owner, role, cell, mat, float(weight),
                                 float(low), float(high), interval)

    def _complete_motif(self, height, depth):
        """Check material, not the arbitrary Film/interface component boundary.

        Only exactly full populations qualify; no occupancy cutoff is used.
        Native stacking roundoff in dimensionless matrices is allowed at 1e-12.
        The actual records and probabilities are never changed by this check.
        """
        for index in range(height-depth, height):
            if index >= 0:
                continue  # Explicit native top-Film growth.
            rank = len(self._intervals) + index
            if rank < 0:
                return False
            records = self._interval_records[self._intervals[rank]]
            if len(records) != 1:
                return False
            record = records[0]
            geometry = self.surface._layer_domain_geometry(index)
            expected = geometry.film_uc
            matrix = self.surface._with_below_translation(
                geometry.film_matrix, expected)
            if (record.owner != self.top_film.name or record.weight != 1.
                    or not np.array_equal(record.cell.basis, expected.basis)
                    or not np.array_equal(record.cell.names, expected.names)
                    or not np.array_equal(record.cell.a, expected.a)
                    or not np.array_equal(record.cell.alpha, expected.alpha)
                    or not np.allclose(record.matrix, matrix, rtol=0, atol=1e-12)
                    or any(record.cell.getLayerCrystallinity(layer)
                           != expected.getLayerCrystallinity(layer)
                           for layer in expected.crystallinity_layers)):
                return False
        return True

    def _extra_records(self, height):
        surface = self.surface
        records = []
        for index in range(max(0, height)):
            geometry = surface._layer_domain_geometry(index)
            matrix = surface._with_below_translation(
                geometry.film_matrix, geometry.film_uc)
            records.append(MaterialRecord(
                self.top_film.name, "growth", geometry.film_uc, matrix, 1.,
                0., 0., ("growth", index)))
        depth = max(len(c.crystallinity_layers)
                    for c in surface.termination_cells.values())
        if (self.exposure_policy == "native_or_unreconstructed"
                and self._complete_motif(height, depth)):
            geometry = surface._layer_domain_geometry(height-1)
            for role, cell, matrix, weight in (
                ("termination", geometry.surface_uc, geometry.surface_matrix,
                 surface._termination_domain_occupancy[geometry.layer_id]),
                ("termination_reference", geometry.reference_uc,
                 geometry.reference_matrix, -1.),
            ):
                records.append(MaterialRecord(
                    self.top_film.name, role, cell,
                    surface._with_below_translation(matrix, cell), weight,
                    0., 0., (role, height)))
        return records

    def _crystal(self, records):
        components = {}
        # Keep domains grouped to avoid one expensive optical component per
        # atomic plane. The native UnitCell still owns all scattering math.
        groups = {}
        for record in records:
            # Zero-weight protected domains still anchor native optical slab
            # boundaries. Removing them changes the native profile grouping.
            group = (self.surface.name if record.role in {
                "growth", "termination", "termination_reference"}
                else record.interval[0])
            key = (group, id(record.cell))
            if key not in groups:
                cell = copy.copy(record.cell)
                cell.basis = cell.basis.copy()
                cell.basis_0 = cell.basis_0.copy()
                cell.name = f"dissolution_{len(groups)}"
                cell.coherentDomainMatrix = []
                cell.coherentDomainOccupancy = []
                groups[key] = cell
                components.setdefault(group, []).append(cell)
            cell = groups[key]
            cell.coherentDomainMatrix.append(record.matrix)
            cell.coherentDomainOccupancy.append(record.weight)
        source = self.source
        cells = [_FrozenMaterial(parts, f"dissolution_{i}", source.reference_area)
                 for i, parts in enumerate(components.values())]
        crystal = SXRDCrystal(
            copy.deepcopy(source.uc_bulk), *cells,
            reference_uc=source.reference_uc,
            reference_rotation=source.reference_rotation, atten=source.atten)
        crystal.weights = np.ones(len(cells))
        crystal.weights_0 = crystal.weights.copy()
        crystal.domains = [copy.deepcopy(self.domains) for _ in cells]
        _copy_geometry(source, crystal)
        return crystal

    def state(self, index):
        """Build one complete state; index zero is the nominal Film boundary."""
        if int(index) != index or index < self.minimum:
            raise ValueError("State index is below exhaustion or is not integer")
        index = int(index)
        keep = len(self._intervals) + min(index, 0)
        retained = [r for r in self._removable if self._rank[r.interval] < keep]
        records = [*self._protected, *retained, *self._extra_records(index)]
        inventory = _inventory(records, self.source.reference_area)
        for owner in self.films:
            inventory.setdefault(owner, {})
        if any(v < -1e-10 for values in inventory.values() for v in values.values()):
            raise ValueError("Exposed motif creates negative material inventory")
        if index <= 0:
            height = self._heights[keep]
        else:
            height = (self.surface.below_H
                      + (index//len(self.surface.layer_order)
                         + self.surface.layerpos[index % len(self.surface.layer_order)]
                         - self.surface.layerpos[0])*self.top_film.unitcell.a[2])
        exposed = (tuple(sorted({r.owner for r in retained
                                 if self._rank[r.interval] == keep-1}))
                   if keep else (self.protected_owner,))
        if index > 0:
            exposed = (self.top_film.name,)
        elif keep:
            interval = self._intervals[keep-1]
            protected = {r.owner or "__bulk__" for r in self._protected
                         if r.interval == interval and r.weight > 0}
            exposed = tuple(sorted(set(exposed) | protected))
        return DissolutionState(index, float(height), exposed, inventory,
                                self._crystal(records))

    def averaged_crystal(self, indices, probabilities):
        """Build the native averaged atomic/optical structure on shared fields."""
        records = []
        for index, probability in zip(indices, probabilities):
            if probability == 0:
                continue
            keep = len(self._intervals) + min(int(index), 0)
            selected = [*self._protected,
                        *(r for r in self._removable if self._rank[r.interval] < keep),
                        *self._extra_records(int(index))]
            records.extend(MaterialRecord(
                r.owner, r.role, r.cell, r.matrix, r.weight*probability,
                r.lower_A, r.upper_A, r.interval) for r in selected)
        return self._crystal(records)

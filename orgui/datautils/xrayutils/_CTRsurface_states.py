"""Material-resolved finite-film dissolution using native crystal states.

The interface closure is frozen and material selective. A structural front
removes finite-owner populations. The original substrate profile, including
its interface extension and signed reference terms, is fixed in every state.
Amplitudes and optical profiles use native unit-cell and crystal evaluators.
"""

import copy
from dataclasses import dataclass

import numpy as np

from .CTRcalc import (
    SXRDCrystal, KinematicAmplitudeResult, KinematicComponentAmplitude,
    _broadcast_kinematic_hkl,
)
from .CTRdistributions import (
    DEFAULT_TAIL_PROBABILITY, RecenteredSkellamProfile, bounded_height_states,
    SkellamProfile,
)
from .CTRfilm import (
    EpitaxyInterface, Film, RoughSurface, _FrozenMaterial,
)
from .CTRuc import UnitCell


@dataclass(frozen=True)
class _MaterialRecord:
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
    optical_bookkeeping: bool = False


@dataclass(frozen=True)
class _SurfaceHeightState:
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


def _native_linear_cells(cells):
    return all(type(cell) is UnitCell and not cell._special_formfactors_present
               and "F_uc" not in cell.__dict__ and "F_bulk" not in cell.__dict__
               for cell in cells)


class _SurfaceStatePlan:
    """Snapshot a connected finite-material path above protected bulk.

    :param SXRDCrystal crystal: Native source; never mutated.
    :param str surface: Name of the RoughSurface carrying native motifs.
    :raises ValueError: For ambiguous topology, geometry, or ownership.
    """

    def __init__(self, crystal, *, surface):
        self.source = copy.deepcopy(crystal, {id(crystal.dwba): None})
        _copy_geometry(crystal, self.source)
        components = self.source.uc_surface_list
        names = [c.name for c in components]
        if len(set(names)) != len(names):
            raise ValueError("Component names must be unique")
        matches = [c for c in components if c.name == surface]
        if len(matches) != 1 or not isinstance(matches[0], RoughSurface):
            raise ValueError("surface must identify one RoughSurface")
        self.surface = matches[0]
        if self.surface.water is not None:
            raise ValueError("Dissolution water placement is not defined")
        if self.surface.crystallinity_profile is not None:
            raise ValueError("Dissolution depth-dependent surface order is not defined")
        if not self.surface._has_finite_boundary():
            raise ValueError("A finite Film boundary is required")
        self.films = self.surface.boundary.films
        self.exposure_policy = self.surface.boundary.exposure_policy
        self.profile = self.surface._make_profile()
        order = np.argsort(self.source.uc_stacking)
        if len(set(self.source.uc_stacking)) != len(components):
            raise ValueError("Dissolution requires an ordered single stack")
        ordered = [components[i] for i in order]
        if sum(isinstance(c, RoughSurface) for c in ordered) != 1:
            raise ValueError("A finite boundary requires exactly one surface")
        surface_slot = ordered.index(self.surface)
        overlayer_indices = order[surface_slot+1:]
        ordered = ordered[:surface_slot+1]
        if any(not isinstance(c, (Film, EpitaxyInterface, RoughSurface))
               for c in ordered):
            raise ValueError("Dissolution requires Film/interface components")
        if any(getattr(c, "crystallinity_profiles", {}) for c in ordered):
            raise ValueError(
                "Dissolution depth-dependent interface order is not defined")
        available = [c.name for c in ordered if isinstance(c, Film)]
        if not available:
            raise ValueError("Dissolution requires at least one finite Film")
        if self.films is None:
            self.films = tuple(available)
        if any(name not in available for name in self.films):
            raise ValueError(
                "Only finite Films below the target surface may dissolve; "
                "bulk is protected")
        selected = [name for name in available if name in self.films]
        if selected != available[-len(selected):]:
            raise ValueError("Removable Films must form a connected upper sequence")
        self.films = tuple(selected)
        if not np.all(self.source.weights[order[:surface_slot+1]] == 1):
            raise ValueError("Dissolution requires unit component weights")
        self.domains = copy.deepcopy(self.source.domains[order[surface_slot]])
        for index in order[:surface_slot+1]:
            domains = self.source.domains[index]
            if (len(domains) != len(self.domains)
                    or any(w != v or not np.array_equal(a, b)
                           for (a, w), (b, v) in zip(domains, self.domains))):
                raise ValueError("Dissolution requires common outer domains")
        # First bind only nominal geometry. The realized mean requires the
        # lower material inventory, so resolve the overlayer anchor afterwards.
        self.source.apply_stacking(_surface_heights={})
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
                component._ensure_optical_bookkeeping()
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
                # Keep the complete initial substrate/interface profile at
                # all heights, including positive-z populations and signed
                # reference terms. Its difference from pristine sharp bulk
                # is initial interface structure, not dissolved substrate.
                protected.append(record)
                continue
            key = _geometry_key(record)
            if key in combined:
                old = combined[key]
                combined[key] = _MaterialRecord(
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
        self.initial_inventory = self._selected_material(0)[3]
        self._overlayers = []
        if len(overlayer_indices):
            _, indices, probabilities, _ = self._snapshot()
            self.overlayer_anchor_A = float(probabilities @ np.array([
                self._front_height(int(index)) for index in indices]))
            self.source.apply_stacking(
                _surface_heights={self.surface.name: self.overlayer_anchor_A})
            self._overlayers = [
                (components[i], self.source.weights[i], self.source.domains[i])
                for i in overlayer_indices]

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
            marked = getattr(cell, "_optical_bookkeeping_domains", ())
            bookkeeping = bool(marked[j]) if j < len(marked) else False
            yield _MaterialRecord(actual_owner, role, cell, mat, float(weight),
                                 float(low), float(high), interval, bookkeeping)

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
            records.append(_MaterialRecord(
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
                records.append(_MaterialRecord(
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
            # Preserve protected zero-weight domains for their geometry;
            # only explicitly marked subtraction slots are optically inert.
            group = (self.surface.name if record.role in {
                "growth", "termination", "termination_reference"}
                else record.interval[0])
            key = (group, id(record.cell))
            if key not in groups:
                cell = copy.copy(record.cell)
                cell.basis = cell.basis.copy()
                cell.basis_0 = cell.basis_0.copy()
                cell.coherentDomainMatrix = []
                cell.coherentDomainOccupancy = []
                cell._optical_bookkeeping_domains = []
                groups[key] = cell
                components.setdefault(group, []).append(cell)
            cell = groups[key]
            cell.coherentDomainMatrix.append(record.matrix)
            cell.coherentDomainOccupancy.append(record.weight)
            cell._optical_bookkeeping_domains.append(record.optical_bookkeeping)
        source = self.source
        cells = [_FrozenMaterial(parts, name, source.reference_area)
                 for name, parts in components.items()]
        count = len(cells)
        # Overlayers are common native components at the realized mean front,
        # not removable material or height-dependent terrace decorations.
        cells.extend(copy.deepcopy(component) for component, _, _ in self._overlayers)
        crystal = SXRDCrystal(
            copy.deepcopy(source.uc_bulk), *cells,
            reference_uc=source.reference_uc,
            reference_rotation=source.reference_rotation, atten=source.atten)
        crystal.enable_uc_stacking = False
        crystal.weights = np.r_[np.ones(count),
                                [weight for _, weight, _ in self._overlayers]]
        crystal.weights_0 = crystal.weights.copy()
        crystal.domains = ([copy.deepcopy(self.domains) for _ in range(count)]
                           + [copy.deepcopy(domains)
                              for _, _, domains in self._overlayers])
        _copy_geometry(source, crystal)
        return crystal

    def _selected_material(self, index):
        """Select and validate material without constructing scattering objects."""
        if int(index) != index or index < self.minimum:
            raise ValueError("State index is below exhaustion or is not integer")
        index = int(index)
        keep = len(self._intervals) + min(index, 0)
        retained = [r for r in self._removable if self._rank[r.interval] < keep]
        extra = self._extra_records(index)
        records = [*self._protected, *retained, *extra]
        inventory = _inventory(records, self.source.reference_area)
        for owner in self.films:
            inventory.setdefault(owner, {})
        if any(v < -1e-10 for values in inventory.values() for v in values.values()):
            raise ValueError("Exposed motif creates negative material inventory")
        return keep, retained, extra, inventory

    def _iter_kinematic_amplitudes(self, indices, h, k, l):  # noqa: E741
        """Stream complete native amplitudes using material prefixes.

        HKL is in the configured reference frame (r.l.u.); amplitudes are in
        electrons per reference lateral cell. Native component-domain scaling
        and FrozenMaterial's cell-area conversion remain the only area factors.
        The plan and all cell copies live only for this request.
        """
        indices = list(indices)
        if any(b <= a for a, b in zip(indices, indices[1:])):
            raise ValueError("Kinematic state indices must be strictly increasing")
        selections = []
        cells = {id(r.cell): r.cell for r in self.records}
        for index in indices:
            keep, _, extra, _ = self._selected_material(index)
            selections.append((int(index), keep, extra))
            cells.update((id(r.cell), r.cell) for r in extra)
        # Preserve custom-cell/callback behavior through the complete-state
        # oracle rather than assuming that an unknown evaluator is additive.
        bulk = self.source.uc_bulk
        native_cells = [bulk, *cells.values()]
        if not _native_linear_cells(native_cells):
            for index in indices:
                yield int(index), self.state(index).crystal.F(h, k, l)
            return
        if not indices:
            return

        coordinates, shape = _broadcast_kinematic_hkl(h, k, l)
        hkl = np.vstack(coordinates)
        common = self._crystal(self._protected)
        protected = common.F(*coordinates)
        # Removable components share outer domains. Their native carrier must
        # not borrow the independent domains of a common overlayer.
        carrier = common if self._protected else self._crystal(
            self._removable[:1])
        templates = {}
        for key, original in cells.items():
            cell = copy.copy(original)
            cell.basis = cell.basis.copy()
            cell.basis_0 = cell.basis_0.copy()
            cell.setReferenceUnitCell(self.source.reference_uc,
                                      self.source.reference_rotation)
            templates[key] = cell

        def _amplitude(records):
            grouped = {}
            for record in records:
                key = id(record.cell)
                if key not in grouped:
                    cell = copy.copy(templates[key])
                    cell.coherentDomainMatrix = []
                    cell.coherentDomainOccupancy = []
                    grouped[key] = cell
                grouped[key].coherentDomainMatrix.append(record.matrix)
                grouped[key].coherentDomainOccupancy.append(record.weight)
            if not grouped:
                return 0.
            component = _FrozenMaterial(
                grouped.values(), "kinematic_terms", self.source.reference_area)
            return sum(carrier._iter_scaled_contributions(
                0, hkl, lambda _unused, h, k, ell: component.F_uc(h, k, ell)))

        retained_amplitude = np.zeros_like(protected)
        growth_amplitude = np.zeros_like(protected)
        previous_keep = previous_growth = 0
        for index, keep, extra in selections:
            for rank in range(previous_keep, keep):
                retained_amplitude += _amplitude(
                    self._interval_records[self._intervals[rank]])
            growth_amplitude += _amplitude(
                r for r in extra if r.role == "growth"
                and r.interval[1] >= previous_growth)
            termination = _amplitude(r for r in extra if r.role != "growth")
            amplitude = (protected + retained_amplitude
                         + growth_amplitude + termination)
            yield index, amplitude.reshape(shape)
            previous_keep, previous_growth = keep, max(0, index)

    def _front_height(self, index):
        """Return the physical structural-front coordinate in Angstrom."""
        if index <= 0:
            return self._heights[len(self._intervals)+index]
        return (self.surface.below_H
                + (index//len(self.surface.layer_order)
                   + self.surface.layerpos[index % len(self.surface.layer_order)]
                   - self.surface.layerpos[0])*self.top_film.unitcell.a[2])

    def state(self, index):
        """Build one complete state; index zero is the nominal Film boundary."""
        keep, retained, extra, inventory = self._selected_material(index)
        index = int(index)
        records = [*self._protected, *retained, *extra]
        height = self._front_height(index)
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
        return _SurfaceHeightState(index, float(height), exposed, inventory,
                                self._crystal(records))

    def averaged_crystal(self, indices, probabilities):
        """Build the native averaged atomic/optical structure on shared fields."""
        # Protected material is independent of the surface distribution. Keep
        # its original coordinates and weights once, including signed interface
        # references; only removable populations carry state probabilities.
        records = list(self._protected)
        for index, probability in zip(indices, probabilities):
            if probability == 0:
                continue
            keep = len(self._intervals) + min(int(index), 0)
            selected = [*(r for r in self._removable if self._rank[r.interval] < keep),
                        *self._extra_records(int(index))]
            records.extend(_MaterialRecord(
                r.owner, r.role, r.cell, r.matrix, r.weight*probability,
                r.lower_A, r.upper_A, r.interval, r.optical_bookkeeping)
                           for r in selected)
        return self._crystal(records)

    def _snapshot(self):
        indices, probabilities, certificate = bounded_height_states(
            self.profile, self.minimum, support=self.surface.height_support)
        return self, indices, probabilities, certificate

    def _average(self):
        _, indices, probabilities, _ = self._snapshot()
        return self.averaged_crystal(indices, probabilities)

    def _evaluate_kinematic(self, h, k, ell):
        average = self._average()
        # Only averaged material must commute with probability weighting.
        # Native overlayers are common components evaluated without averaging.
        cells = [average.uc_bulk, *(cell for c in average.uc_surface_list
                                   if isinstance(c, _FrozenMaterial)
                                   for cell in c.cells)]
        if _native_linear_cells(cells):
            return average.evaluate_kinematic(h, k, ell)
        # Unknown evaluators need not commute with probability averaging.
        # Keep their complete-state behavior and component decomposition.
        _, indices, probabilities, _ = self._snapshot()
        bulk = total = 0.
        components = {}
        for index, probability in zip(indices, probabilities):
            if not probability:
                continue
            result = self.state(int(index)).crystal.evaluate_kinematic(h, k, ell)
            bulk = bulk + probability * result.bulk
            total = total + probability * result.total
            for component in result.components:
                components[component.name] = (components.get(component.name, 0.)
                                              + probability * component.amplitude)
        return KinematicAmplitudeResult(bulk, tuple(
            KinematicComponentAmplitude(index, name, value)
            for index, (name, value) in enumerate(components.items())), total)

    def _optical_snapshot(self):
        """Retain analytic references; otherwise use the all-state union grid."""
        from ._CTRincoherent_dwba import _shared_grid, _regrid_model
        _, indices, probabilities, _ = self._snapshot()
        average = self.averaged_crystal(indices, probabilities)
        records, model, profile = average.dwba._atomic_model_snapshot()
        if profile.continuous is not None:
            # Native analytic water requires its onset/reference metadata.
            # A sampled union grid cannot replace that physical reference.
            # State-specific ensembles use each native reference separately;
            # shared continuous references retain their explicit rejection.
            return records, model, profile
        boundaries = [self.state(int(index)).crystal.dwba
                      ._atomic_model_snapshot()[2].boundaries
                      for index, probability in zip(indices, probabilities)
                      if probability]
        grid = _shared_grid(profile, np.concatenate(boundaries))
        return records, _regrid_model(model, profile, grid.values[:, 0]), grid
    def statistics(self):
        """Return realized front moments and species inventories, with tail audit.

        Heights are Angstrom; atom counts are per reference lateral area.
        A porous interface's front coordinate is not a unique atomic surface.
        """
        from dataclasses import asdict
        path, indices, probabilities, certificate = self._snapshot()
        states = [path.state(int(j)) for j in indices]
        heights = np.array([s.front_height_A for s in states])
        mean_index = float(probabilities @ indices)
        variance_index = float(probabilities @ (indices-mean_index)**2)
        mean = float(probabilities @ heights)
        variance = float(probabilities @ (heights-mean)**2)
        inventory = {}
        exposed = {}
        for state, p in zip(states, probabilities):
            for owner, atoms in state.inventory.items():
                values = inventory.setdefault(owner, {})
                for species, number in atoms.items():
                    values[species] = values.get(species, 0.) + float(p*number)
            for owner in state.exposed_owners:
                exposed[owner] = exposed.get(owner, 0.) + float(p)
        # Report additions/removals state by state: changes of termination
        # chemistry must not be hidden by cancellation in the ensemble mean.
        added, removed = {}, {}
        for state, p in zip(states, probabilities):
            for owner in path.films:
                initial = path.initial_inventory.get(owner, {})
                current = state.inventory.get(owner, {})
                for species in initial.keys() | current.keys():
                    delta = current.get(species, 0.) - initial.get(species, 0.)
                    for target, value in ((added, max(delta, 0.)),
                                          (removed, max(-delta, 0.))):
                        atoms = target.setdefault(owner, {})
                        atoms[species] = atoms.get(species, 0.) + float(p*value)
        # For omitted positive indices, |z| <= |z0| + a_c*(J+1).
        # This conservative bound covers unequal top-cycle spacings. Include
        # the explicitly reported redistribution of upper-tail probability.
        scale = float(path.top_film.unitcell.a[2])
        intercept = abs(path.nominal_height_A) + scale
        first_bound = (intercept*certificate.omitted_probability
                       + scale*certificate.omitted_first_moment)
        second_bound = (intercept**2*certificate.omitted_probability
                        + 2*intercept*scale*certificate.omitted_first_moment
                        + scale**2*certificate.omitted_second_moment)
        first_bound += certificate.omitted_probability*float(np.max(abs(heights)))
        second_bound += certificate.omitted_probability*float(np.max(heights**2))
        return dict(parent_mean_index=self.profile.expected_height_change,
                    parent_variance_index2=self._parent_variance(),
                    realized_mean_index=mean_index,
                    realized_variance_index2=variance_index,
                    mean_front_height_A=mean, variance_front_height_A2=variance,
                    rms_front_height_A=float(np.sqrt(variance)),
                    exhausted_probability=float(probabilities[0]),
                    exhausted_index=path.minimum, exposed_owner_probabilities=exposed,
                    initial_inventory=path.initial_inventory,
                    remaining_inventory=inventory, added_inventory=added,
                    removed_inventory=removed, support=asdict(certificate),
                    height_first_moment_error_bound_A=first_bound,
                    height_second_moment_error_bound_A2=second_bound,
                    interface_support_audit=path.interface_support_audit,
                    interface_owners=path.interface_owners)

    def _parent_variance(self):
        from .CTRdistributions import PoissonProfile
        if isinstance(self.profile, PoissonProfile):
            _, _, fraction = self.profile._step_parameters()
            return self.profile.rate + fraction * (1 - fraction)
        return self.profile.variance

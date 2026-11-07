"""Native complete-state DWBA for statistical terrace ensembles.

Lengths are Angstrom, angles radians, hkl reference-cell r.l.u. All mixing
uses dimensionless reflection amplitudes in one ambient phase gauge. Terrace
caches last only for one explicit batch or standalone request.
"""

import copy
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, replace
from types import MethodType

import numpy as np

from .CTRdistributions import PoissonProfile
from .CTRdwba import DWBAState, _PreparedGeometry, _hash_arrays, _readonly_result
from .CTRoptics import StratifiedProfile
from .CTRuc import WaterModel


def _flat_profile(surface):
    return PoissonProfile(0, offset=surface._dwba_height + 1)


def _flat_crystal(model, height):
    """Freeze common components at their mean positions and replace the target.

    The existing surface assembler builds deterministic signed Film removal
    and termination records, including local-depth crystallinity. Copying
    avoids temporary mutations of source cells or of their cached wavefields.
    """
    original = model.coherent_model
    crystal = copy.deepcopy(original, {id(original.dwba): None})
    crystal.enable_uc_stacking = False
    crystal._stack_policy = None
    target = next(c for c in crystal.uc_surface_list if c.name == model.surface)
    target._dwba_height = int(height)
    target._make_profile = MethodType(_flat_profile, target)
    target.createLayers()
    crystal.dwba.set_orientation(original.dwba.orientation)
    crystal.dwba.set_field_cache_capacity(original.dwba._field_cache_capacity)
    return crystal


def _geometry(prepared):
    return _PreparedGeometry(
        prepared.shape, prepared.scalar, prepared.hkl.reshape(3, -1),
        prepared.alpha_i, prepared.alpha_f, prepared.vlieg_angles,
    )


def _phase_to(prepared, common):
    """Translate both ambient wave origins before coherently adding states."""
    def _origin(field):
        return np.asarray(field.z_reference).reshape(-1)[0]

    dz_i = _origin(prepared.field_i) - _origin(common.field_i)
    dz_f = _origin(prepared.field_f) - _origin(common.field_f)
    phase = np.exp(1j * common.k0 * common.reference.n[0] * (
        np.sin(common.alpha_i) * dz_i + np.sin(common.alpha_f) * dz_f))
    return phase.reshape(common.shape)


def _shared_grid(shared, state_boundaries):
    """Split shared optics at representable state interfaces (Angstrom).

    Keep every shared optical boundary at its original position. A state-only
    split is numerically coincident when no float lies strictly between it
    and a neighboring boundary. Prefer the shared boundary in that case;
    never merge optical media using a distance or optical-contrast tolerance.
    Resolvable state splits retain the piecewise-constant reference shares.
    """
    boundaries = np.union1d(shared.boundaries, state_boundaries)
    optical = np.isin(boundaries, shared.boundaries)
    kept = []
    for index, boundary in enumerate(boundaries):
        if kept and np.nextafter(boundaries[kept[-1]], np.inf) >= boundary:
            if not optical[index]:
                continue
            if not optical[kept[-1]]:
                kept[-1] = index
                continue
            # Two invalid optical interfaces must still reach strict validation.
        kept.append(index)
    boundaries = boundaries[kept]
    centers = np.r_[boundaries[0] - 1,
                    boundaries[:-1] + np.diff(boundaries) / 2,
                    boundaries[-1] + 1]
    values = np.column_stack((centers, shared.values[
        np.searchsorted(shared.boundaries, centers), 1:]))
    return StratifiedProfile(values, boundaries)


def _shared_preparation(state, common, atomic_model, profile):
    """Integrate a state's density on fields of the averaged optical medium.

    Split the averaged slabs at representable state boundaries, preserving n(z).
    Rebin state reference shares by exact piecewise-constant lookup, not
    interpolation. This lets the existing native crystallinity machinery
    evaluate c*A_state + (1-c)*R_state on the shared fields, including layer
    overrides, while the common reference is subtracted exactly once later.
    """
    shared = common.reference.profile
    if common.reference.continuous is not None or profile.continuous is not None:
        raise NotImplementedError(
            "Shared terrace fields require a layered optical reference; "
            "use state_specific for continuous water references."
        )
    grid = _shared_grid(shared, profile.boundaries)
    centers = grid.values[:, 0]
    state_indices = np.searchsorted(profile.boundaries, centers)
    reference = state._reference(grid)
    # Stored shares run from ambient to substrate, opposite to profile z.
    indices = (len(profile.values) - 1 - state_indices)[::-1]
    shares = atomic_model.reference_shares[:, indices, :]
    layer_shares = atomic_model.layer_reference_shares[:, indices, :]
    model = replace(
        atomic_model, reference_shares=shares,
        layer_reference_shares=layer_shares,
        decomposition_fingerprint=_hash_arrays(shares, layer_shares),
    )
    return state._prepare_explicit(
        _geometry(common), reference, model,
        common.polarization_i, common.polarization_f, 1e-7, 1e-8,
    )


def _regrid_model(model, profile, centers):
    """Rebin planar density shares exactly, with ambient-first array ordering."""
    indices = (len(profile.values) - 1
               - np.searchsorted(profile.boundaries, centers))[::-1]
    shares = np.ascontiguousarray(model.reference_shares[:, indices, :])
    layers = np.ascontiguousarray(model.layer_reference_shares[:, indices, :])
    return replace(model, reference_shares=shares, layer_reference_shares=layers,
                   decomposition_fingerprint=_hash_arrays(shares, layers))


class _SharedTerraceBasis:
    """Linear flat-height corrections on one shared field, without state copies.

    For fully ordered dry records, or off specular with fixed layer order,
    the planar state density never enters F_atomic. Evaluate each Film layer
    and termination once, then form complete states with prefix sums. Common
    substrate/Film interference remains in the native amplitude.
    """

    def __init__(self, model, prepared, heights, records):
        target = model.target_surface()
        template = next(r for r in records if r.component_index >= 0
                        and r.component_name == model.surface)
        self.component_index = template.component_index
        self.state = DWBAState(model.coherent_model)
        basis_records = [records[0]]
        self.film_rows = {}
        self.termination_rows = {}
        occupancy = target.underlying_film.unitcell.coherentDomainOccupancy[0]

        def _append(cell, matrix, weight):
            cell = copy.copy(cell)
            # Only scattering occupancy is weighted. Optical fields remain
            # those of the unmodified physical density. Off specular the
            # planar (1-c)*R term vanishes exactly.
            cell.basis = cell.basis.copy()
            cell.basis[:, 6] *= [cell.getLayerCrystallinity(layer)
                                 for layer in cell._crystallinity_labels()]
            cell.coherentDomainMatrix = [target._with_below_translation(matrix, cell)]
            cell.coherentDomainOccupancy = [weight]
            index = len(basis_records)
            basis_records.append(replace(template, cell=cell, record_index=index))
            return index

        for height in range(min(0, int(heights[0]) + 1),
                            max(-1, int(heights[-1])) + 1):
            geometry = target._layer_domain_geometry(height)
            self.film_rows[height] = _append(
                geometry.film_uc, geometry.film_matrix, occupancy)
        for height in heights:
            geometry = target._layer_domain_geometry(int(height))
            surface = _append(geometry.surface_uc, geometry.surface_matrix,
                              target._termination_domain_occupancy[geometry.layer_id])
            reference = _append(geometry.reference_uc, geometry.reference_matrix,
                                -occupancy)
            self.termination_rows[int(height)] = surface, reference
        self.records = tuple(basis_records)
        self.packed = self.state._pack_atomic_records(prepared, self.records)
        # The coherent result already contains the complete substrate. The
        # native finite-record kernel accepts an empty bulk atom list.
        for name in ("bulk_positions", "bulk_atom_index", "bulk_weight",
                     "bulk_repeat_coordinate"):
            self.packed[name] = self.packed[name][:0]

    def _intensity(self, coherent, heights, probabilities, attenuation):
        prepared = coherent.prepared
        shares = np.zeros((len(self.records), len(prepared.reference.n), 2))
        atomic_model = replace(prepared._atomic_model, reference_shares=shares)
        raw_prepared = replace(prepared, _atomic_model=atomic_model)
        atomic, _ = self.state._evaluate_raw(
            raw_prepared, self.records, self.packed, "semi_infinite", attenuation)
        old_surface = sum(part.F_atomic for part in coherent.contributions
                          if part.component_index == self.component_index)
        common = coherent.total_amplitude - coherent.amplitude_prefactor * old_surface
        # Prefixes are anchored at height -1 (the sharp Film boundary).
        prefixes = {-1: np.zeros(prepared.alpha_i.size, dtype=complex)}
        running = prefixes[-1]
        for height in range(0, max(-1, int(heights[-1])) + 1):
            running = running + atomic[self.film_rows[height]]
            prefixes[height] = running
        running = prefixes[-1]
        for height in range(-2, min(-1, int(heights[0])) - 1, -1):
            running = running - atomic[self.film_rows[height + 1]]
            prefixes[height] = running
        intensity = np.zeros(prepared.shape)
        for height, probability in zip(heights, probabilities):
            surface, reference = self.termination_rows[int(height)]
            correction = prefixes[int(height)] + atomic[surface] + atomic[reference]
            amplitude = common + coherent.amplitude_prefactor * correction.reshape(
                prepared.shape)
            intensity += probability * np.abs(amplitude)**2
        return intensity


@dataclass(frozen=True)
class DWBAEnsembleResult:
    """One polarization channel's mixed field intensity, without an amplitude.

    ``reflectivity`` is dimensionless ``|r|^2`` (no normal-flux correction), as
    in DWBAResult. ``amplitude_prefactor`` converts field intensity to an
    effective squared structure factor by dividing by its squared modulus.
    Convolve ``reflectivity`` before that conversion or taking a square root.
    ``coherent_intensity`` and ``incoherent_intensity`` are the endpoints;
    the latter is None when the shared coherent fast path skips the states.
    ``prepared`` describes the common geometry and averaged optical reference;
    state-specific optical fields are evaluation-local.
    """

    prepared: object
    reflectivity: object
    amplitude_prefactor: object
    coherent_intensity: object
    incoherent_intensity: object
    optical_reference: str


class SurfaceHeightDWBA:
    """DWBA evaluation facade returned by SurfaceHeightDomains.dwba.

    Geometry, orientation, and field-cache controls belong to the wrapped
    crystal's ``dwba`` state. Terrace snapshots live only inside ``batch``;
    individual evaluations automatically open and close their own batch.
    """

    def __init__(self, model):
        self._model = model
        self._scope = None

    def __getstate__(self):
        # Optimizers deepcopy fitted models. Never copy a live ExitStack or
        # carry its frozen terrace states into a new owner.
        return {"_model": self._model, "_scope": None}

    @contextmanager
    def batch(self):
        """Reuse terraces, packing and fields across channels and rods.

        The model must not be mutated inside this scope. Nested scopes share
        the outer snapshot. All terrace data is released at exit, including
        after exceptions; the next batch reads current fitted parameters.
        """
        if self._scope is not None:
            yield self
            return
        self._model.validate(forward_model="dwba")
        with ExitStack() as stack:
            stack.enter_context(self._model.coherent_model.dwba.batch())
            self._scope = {"stack": stack}
            try:
                yield self
            finally:
                self._scope = None

    def _selection(self):
        if "selection" not in self._scope:
            model = self._model
            target = model.target_surface()
            candidates = target._profile_candidates()
            low, high, probabilities, *_ = target._height_state_selection(
                candidates, model.exact_layer_count)
            if abs(float(np.sum(probabilities)) - 1) > model.probability_tolerance:
                raise ValueError("Terrace probabilities must sum to one")
            heights = candidates.layer_numbers[low:high + 1]
            active = probabilities > 0
            self._scope["selection"] = heights[active], probabilities[active]
        return self._scope["selection"]

    def _shared_basis(self, prepared, heights):
        if "basis_eligibility" not in self._scope:
            native = self._model.coherent_model.dwba
            records, _, _ = native._atomic_model_snapshot()
            target = self._model.target_surface()
            dry = (target.crystallinity_profile is None and all(
                not isinstance(record.cell, WaterModel) for record in records))
            ordered = dry and all(
                all(record.cell.getLayerCrystallinity(layer) == 1
                    for layer in record.cell.crystallinity_layers)
                for record in records)
            self._scope["basis_eligibility"] = dry, ordered, records
        dry, ordered, records = self._scope["basis_eligibility"]
        if not dry or (np.any(prepared.is_specular) and not ordered):
            return None
        if "basis" not in self._scope:
            self._scope["basis"] = _SharedTerraceBasis(
                self._model, prepared, heights, records)
        return self._scope["basis"]

    def _states(self, heights):
        if "states" not in self._scope:
            states = []
            for height in heights:
                state = _flat_crystal(self._model, height).dwba
                self._scope["stack"].enter_context(state.batch())
                # Private states form reference cycles with their crystals.
                # Release field arrays before waiting for cyclic collection.
                self._scope["stack"].callback(state.clear_cache)
                _, atomic, profile = state._atomic_model_snapshot()
                states.append((state, atomic, profile))
            self._scope["states"] = states
        return self._scope["states"]

    def _shared_states(self, common, states):
        native = self._model.coherent_model.dwba
        if "shared" not in self._scope:
            profile = common.reference.profile
            if common.reference.continuous is not None or any(
                p.continuous is not None for _, _, p in states
            ):
                raise NotImplementedError(
                    "Shared terrace fields require a layered optical reference; "
                    "use state_specific for continuous water references.")
            grid = _shared_grid(profile, np.concatenate(
                [p.boundaries for _, _, p in states]))
            centers = grid.values[:, 0]
            reference = native._reference(grid)
            atomic = _regrid_model(common._atomic_model, profile, centers)
            models = [_regrid_model(m, p, centers) for _, m, p in states]
            self._scope["shared"] = reference, atomic, models
        reference, atomic, models = self._scope["shared"]
        shared = native._prepare_explicit(
            _geometry(common), reference, atomic,
            common.polarization_i, common.polarization_f, 1e-7, 1e-8)
        for (state, _, _), atomic in zip(states, models):
            yield state, replace(
                shared, crystal_identity=id(state._crystal),
                bulk_identity=id(state._crystal.uc_bulk), _atomic_model=atomic)

    def _state_intensities(self, coherent, attenuation):
        model = self._model
        prepared = coherent.prepared
        heights, probabilities = self._selection()
        if model.optical_reference == "shared":
            basis = self._shared_basis(prepared, heights)
            if basis is not None:
                return coherent.reflectivity, basis._intensity(
                    coherent, heights, probabilities, attenuation)
        states = self._states(heights)
        if model.optical_reference == "shared":
            preparations = self._shared_states(prepared, states)
        else:
            preparations = ((state, state._prepare_geometry(
                _geometry(prepared), polarization_i=prepared.polarization_i,
                polarization_f=prepared.polarization_f, geometry_rtol=1e-7,
                geometry_atol=1e-8, cache=False)) for state, _, _ in states)
        amplitude_sum = np.zeros(prepared.shape, dtype=complex)
        incoherent = np.zeros(prepared.shape, dtype=float)
        prefactor = coherent.amplitude_prefactor
        reference_amplitude = coherent.F_reference
        for probability, (state, state_prepared) in zip(probabilities, preparations):
            result = state.evaluate_prepared(
                state_prepared, bulk_attenuation=attenuation)
            phase = _phase_to(state_prepared, prepared)
            if model.optical_reference == "shared":
                amplitude = coherent.unperturbed_amplitude + prefactor * (
                    result.F_atomic * phase - reference_amplitude)
            else:
                amplitude = result.total_amplitude * phase
            amplitude_sum += probability * amplitude
            incoherent += probability * np.abs(amplitude)**2
        coherent_intensity = (coherent.reflectivity
                              if model.optical_reference == "shared"
                              else np.abs(amplitude_sum)**2)
        return coherent_intensity, incoherent

    def _evaluate(self, method, args, kwargs):
        with self.batch():
            return self._evaluate_in_batch(method, args, kwargs)

    def _evaluate_in_batch(self, method, args, kwargs):
        model = self._model
        fraction = model._validated_fraction()
        kwargs = dict(kwargs)
        bulk_mode = kwargs.pop("bulk_mode", "semi_infinite")
        if bulk_mode != "semi_infinite":
            raise ValueError("Terrace intensities require semi_infinite bulk mode")
        attenuation = kwargs.pop("bulk_attenuation", 0.0)
        native = model.coherent_model.dwba
        prepared = getattr(native, method)(*args, **kwargs)
        coherent = native.evaluate_prepared(prepared, bulk_attenuation=attenuation)
        coherent_intensity = coherent.reflectivity
        incoherent = None
        if fraction != 0 or model.optical_reference == "state_specific":
            coherent_intensity, incoherent = self._state_intensities(
                coherent, attenuation)
        intensity = (coherent_intensity if incoherent is None else
                     (1 - fraction) * coherent_intensity + fraction * incoherent)
        if not np.all(np.isfinite(intensity)) or np.any(intensity < 0):
            raise ValueError("Terrace DWBA intensity must be finite and nonnegative")

        def _freeze(value):
            return (None if value is None else _readonly_result(
                value, prepared.shape, prepared.scalar))

        return DWBAEnsembleResult(
            prepared, _freeze(intensity), coherent.amplitude_prefactor,
            _freeze(coherent_intensity), _freeze(incoherent),
            model.optical_reference)

    def evaluate(self, h, k, l, **kwargs):  # noqa: E741
        """Evaluate configured hkl (r.l.u.) with native DWBA channel options."""
        return self._evaluate("prepare", (h, k, l), kwargs)

    def evaluate_from_glancing(self, h, k, alpha_i, alpha_f, **kwargs):
        """Evaluate h, k in r.l.u. and glancing angles in radians."""
        return self._evaluate(
            "prepare_from_glancing", (h, k, alpha_i, alpha_f), kwargs)

    def evaluate_from_vlieg(self, alpha, delta, gamma, omega, chi, phi, **kwargs):
        """Evaluate six measured Vlieg angles in radians."""
        return self._evaluate(
            "prepare_from_vlieg", (alpha, delta, gamma, omega, chi, phi), kwargs)

    def _polarized(self, evaluator, args, s_fraction, outgoing, kwargs):
        native = self._model.coherent_model.dwba
        with self.batch():
            return native._polarization_reflectivity(
                evaluator, args, s_fraction=s_fraction, outgoing=outgoing,
                kwargs=kwargs)

    def polarization_reflectivity(
        self, h, k, l, *, s_fraction=1.0, outgoing="unanalysed", **kwargs  # noqa: E741
    ):
        """Sum incident/exit channel intensities at hkl in r.l.u."""
        return self._polarized(
            self.evaluate, (h, k, l), s_fraction, outgoing, kwargs)

    def polarization_reflectivity_from_glancing(
        self, h, k, alpha_i, alpha_f, *, s_fraction=1.0,
        outgoing="unanalysed", **kwargs
    ):
        """Sum channels using glancing angles in radians and h, k in r.l.u."""
        return self._polarized(self.evaluate_from_glancing,
                               (h, k, alpha_i, alpha_f),
                               s_fraction, outgoing, kwargs)

    def polarization_reflectivity_from_vlieg(
        self, alpha, delta, gamma, omega, chi, phi, *, s_fraction=1.0,
        outgoing="unanalysed", **kwargs
    ):
        """Sum channels using six Vlieg angles in radians."""
        return self._polarized(self.evaluate_from_vlieg,
                               (alpha, delta, gamma, omega, chi, phi),
                               s_fraction, outgoing, kwargs)

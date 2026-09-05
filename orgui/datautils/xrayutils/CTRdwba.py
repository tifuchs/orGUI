"""Distorted-wave Born approximation state and observable amplitudes."""

from collections import OrderedDict
from dataclasses import dataclass
import hashlib

import numpy as np

from ._CTRnative import HAS_CPP_ACCEL, _CTRcalc_cpp
from .CTRfilm import EpitaxyInterface, Film, PoissonSurface
from .CTRoptics import (
    HC_KEV_ANGSTROM,
    LayeredElectricField,
    StratifiedProfile,
    homogeneous_bulk_profile,
    profile_boundaries,
    solve_electric_field,
    top_layer_spacing,
)
from .CTRuc import UnitCell, WaterModel, _coherent_domain_arrays
from .HKLVlieg import UBCalculator, VliegAngles


__all__ = [
    "DWBAContribution",
    "DWBAResult",
    "DWBAState",
    "PreparedCTR",
    "PreparedOpticalReference",
]


_CLASSICAL_ELECTRON_RADIUS_ANGSTROM = 2.8179403262e-5
_PREPARATION_CACHE_SIZE = 32
_FIELD_CACHE_SIZE = 4096
_OPTICAL_DELTA_TOLERANCE = 1e-9
_OPTICAL_BETA_TOLERANCE = 1e-9


def _readonly_array(value, dtype=None):
    array = np.array(value, dtype=dtype, order="C", copy=True)
    array.flags.writeable = False
    return array


def _readonly_result(value, shape, scalar, dtype=None):
    array = np.asarray(value, dtype=dtype).reshape(shape)
    if scalar:
        return array.item()
    return _readonly_array(array)


def _hash_arrays(*values):
    digest = hashlib.blake2b(digest_size=20)
    for value in values:
        array = np.ascontiguousarray(value)
        digest.update(array.dtype.str.encode("ascii"))
        digest.update(np.asarray(array.shape, dtype=np.int64).tobytes())
        digest.update(array.tobytes())
    return digest.digest()


def _format_hkl(hkl, index):
    return "(" + ", ".join(f"{value:g}" for value in hkl[:, index]) + ")"


@dataclass(frozen=True)
class _LiveAtomicRecord:
    component_index: int
    component_name: str
    component_type: str
    record_index: int
    record_name: str
    role: str
    layer: float | None
    cell: UnitCell
    component_weight: float
    outer_domains: tuple
    identity_token: tuple


@dataclass(frozen=True)
class _PreparedAtomicModel:
    descriptors: tuple
    reference_shares: np.ndarray
    geometry_fingerprint: bytes
    decomposition_fingerprint: bytes


def _record_layer(cell):
    layers = np.unique(np.asarray(cell.basis[:, 7], dtype=np.float64))
    if layers.size == 1:
        return float(layers[0])
    return None


def _copy_parent_domains(parent, layer):
    layer.coherentDomainMatrix = parent.coherentDomainMatrix
    layer.coherentDomainOccupancy = parent.coherentDomainOccupancy
    return layer


def _component_record_cells(component):
    if isinstance(component, UnitCell):
        layers = component.split_in_layers(ordered=True)
        if component.layer_behavior == "select":
            selected = component.build_selected_basis()[0]
            selected_layers = set(np.asarray(selected[:, 7], dtype=np.float64))
            layers = OrderedDict(
                (layer_id, layer)
                for layer_id, layer in layers.items()
                if float(layer_id) in selected_layers
            )
        return [
            ("unit_cell_layer", _copy_parent_domains(component, layer))
            for layer in layers.values()
        ]
    if isinstance(component, Film):
        if np.any(component._basis_created != component.basis):
            component.createLayers()
        return [("film_layer", cell) for cell in component.layer_ucs]
    if isinstance(component, EpitaxyInterface):
        if np.any(component._basis_created != component.basis):
            component.createInterfaceCells()
        return [
            *(("interface_top", cell) for cell in component.top_layers),
            *(("interface_bottom", cell) for cell in component.bottom_layers),
        ]
    if isinstance(component, PoissonSurface):
        if np.any(component._basis_created != component.basis):
            component.createLayers()
        records = []
        for surface_cell, film_cell in zip(
            component.layer_ucs, component.film_layer_ucs
        ):
            records.append(("surface_termination", surface_cell))
            records.append(("covered_film", film_cell))
        records.extend(
            ("sharp_film_correction", cell)
            for cell in component._film_termination_ucs.values()
        )
        return records
    if isinstance(component, WaterModel):
        raise NotImplementedError(
            "DWBA does not support continuous-density WaterModel components."
        )
    raise NotImplementedError(
        f"DWBA does not support {type(component).__name__} components."
    )


def _record_descriptor(record):
    return (
        record.component_index,
        record.component_name,
        record.component_type,
        record.record_index,
        record.record_name,
        record.role,
        record.layer,
    )


def _hash_record_geometry(records):
    digest = hashlib.blake2b(digest_size=20)
    for record in records:
        digest.update(repr(record.identity_token).encode("utf-8"))
        cell = record.cell
        digest.update(
            _hash_arrays(
                np.asarray(cell.basis.shape, dtype=np.int64),
                np.asarray(cell.basis[:, 7], dtype=np.float64),
                cell.B_mat,
                cell.R_mat,
                cell.R_mat_inv,
                cell.refHKLTransform,
                cell.refRealTransform,
                np.asarray([getattr(cell, "_E", np.nan)], dtype=np.float64),
                np.asarray(cell.coherentDomainMatrix, dtype=np.float64).reshape(
                    (-1, 3, 4)
                ),
                np.asarray(
                    [domain[0] for domain in record.outer_domains],
                    dtype=np.float64,
                ).reshape((-1, 3, 3)),
            )
        )
    return digest.digest()


def _combine_tagged(values, shares, z_tolerance=0.6):
    values = np.asarray(values, dtype=np.float64)
    shares = np.asarray(shares, dtype=np.float64)
    if len(values) == 0:
        return values.reshape((0, 3)), shares[:, :0]
    order = np.argsort(values[:, 0], kind="stable")
    values = values[order]
    shares = shares[:, order]
    group_starts = [0]
    for index in range(1, len(values)):
        if values[index, 0] - values[group_starts[-1], 0] > z_tolerance:
            group_starts.append(index)
    group_stops = [*group_starts[1:], len(values)]
    combined = np.empty((len(group_starts), 3), dtype=np.float64)
    combined_shares = np.zeros(
        (shares.shape[0], len(group_starts), 2), dtype=np.float64
    )
    for group_index, (start, stop) in enumerate(zip(group_starts, group_stops)):
        combined_shares[:, group_index] = shares[:, start:stop].sum(axis=1)
        combined[group_index, 0] = values[start, 0]
        combined[group_index, 1:] = combined_shares[:, group_index].sum(axis=0)
    return np.ascontiguousarray(combined), np.ascontiguousarray(combined_shares)


def _stratify_tagged(profile, shares, delta_tolerance, beta_tolerance):
    profile = np.asarray(profile, dtype=np.float64)
    shares = np.asarray(shares, dtype=np.float64)
    boundaries = profile_boundaries(profile)
    if delta_tolerance is None or len(profile) <= 3:
        return (
            StratifiedProfile(
                np.ascontiguousarray(profile.copy()),
                np.ascontiguousarray(boundaries),
            ),
            np.ascontiguousarray(shares.copy()),
        )
    if beta_tolerance is None:
        beta_tolerance = delta_tolerance
    if delta_tolerance < 0.0 or beta_tolerance < 0.0:
        raise ValueError("delta and beta tolerances must be non-negative.")

    values = [profile[0].copy()]
    tagged = [shares[:, 0].copy()]
    inherited_boundaries = [boundaries[0]]
    start = 1
    last_finite = len(profile) - 2
    while start <= last_finite:
        stop = start
        delta_min = delta_max = profile[start, 1]
        beta_min = beta_max = profile[start, 2]
        while stop < last_finite:
            candidate = stop + 1
            next_delta_min = min(delta_min, profile[candidate, 1])
            next_delta_max = max(delta_max, profile[candidate, 1])
            next_beta_min = min(beta_min, profile[candidate, 2])
            next_beta_max = max(beta_max, profile[candidate, 2])
            if (
                next_delta_max - next_delta_min > delta_tolerance
                or next_beta_max - next_beta_min > beta_tolerance
            ):
                break
            stop = candidate
            delta_min, delta_max = next_delta_min, next_delta_max
            beta_min, beta_max = next_beta_min, next_beta_max
        lower = boundaries[start - 1]
        upper = boundaries[stop]
        thicknesses = np.diff(boundaries[start - 1 : stop + 1])
        merged_shares = np.average(
            shares[:, start : stop + 1], axis=1, weights=thicknesses
        )
        merged = profile[start].copy()
        merged[0] = 0.5 * (lower + upper)
        merged[1:] = merged_shares.sum(axis=0)
        values.append(merged)
        tagged.append(merged_shares)
        inherited_boundaries.append(upper)
        start = stop + 1
    values.append(profile[-1].copy())
    tagged.append(shares[:, -1].copy())
    return (
        StratifiedProfile(
            np.ascontiguousarray(np.asarray(values, dtype=np.float64)),
            np.ascontiguousarray(np.asarray(inherited_boundaries, dtype=np.float64)),
        ),
        np.ascontiguousarray(np.stack(tagged, axis=1)),
    )


@dataclass(frozen=True)
class PreparedOpticalReference:
    """Immutable optical reference shared by compatible CTR preparations.

    :param StratifiedProfile profile:
        Frozen layer-center optical constants and physical boundaries.
    :param numpy.ndarray n:
        Complex refractive indices ordered from ambient to substrate.
    :param float energy_eV:
        Photon energy in eV.
    :param float k0:
        Vacuum wavevector magnitude in inverse Angstrom.
    :param bytes fingerprint:
        Exact content fingerprint used by the automatic cache.
    """

    profile: StratifiedProfile
    n: np.ndarray
    energy_eV: float
    k0: float
    fingerprint: bytes


@dataclass(frozen=True)
class PreparedCTR:
    """Immutable DWBA geometry and frozen optical reference for one scan."""

    shape: tuple
    scalar: bool
    hkl: np.ndarray
    vlieg_angles: np.ndarray | None
    Q_parallel: np.ndarray
    is_specular: np.ndarray
    alpha_i: np.ndarray
    alpha_f: np.ndarray
    cos_azimuth: np.ndarray
    sin_azimuth: np.ndarray
    polarization_i: str
    polarization_f: str
    reference: PreparedOpticalReference
    field_i: LayeredElectricField
    field_f: LayeredElectricField
    field_i_index: np.ndarray
    field_f_index: np.ndarray
    crystal_identity: int
    bulk_identity: int
    energy_eV: float
    B_mat: np.ndarray
    R_mat: np.ndarray
    R_mat_inv: np.ndarray
    ref_hkl_transform: np.ndarray
    orientation: np.ndarray
    reference_area: float
    _atomic_model: _PreparedAtomicModel

    @property
    def k0(self):
        """Return the frozen vacuum wavevector in inverse Angstrom."""
        return self.reference.k0


@dataclass(frozen=True)
class DWBAContribution:
    """One generated atomic record's coherent DWBA amplitudes.

    ``component_index`` is ``-1`` for the bulk and otherwise indexes
    ``SXRDCrystal.uc_surface_list``. Amplitudes use electrons per configured
    reference lateral cell and have the same scalar or array shape as the
    parent :class:`DWBAResult`.
    """

    component_index: int
    component_name: str
    component_type: str
    record_index: int
    record_name: str
    role: str
    layer: float | None
    F_atomic: complex | np.ndarray
    F_reference: complex | np.ndarray
    F_contrast: complex | np.ndarray


@dataclass(frozen=True)
class DWBAResult:
    """DWBA matrix element and precisely normalized derived observables."""

    prepared: PreparedCTR
    F_contrast: complex | np.ndarray
    unperturbed_amplitude: complex | np.ndarray
    scattered_amplitude: complex | np.ndarray
    total_amplitude: complex | np.ndarray
    bulk_mode: str
    contributions: tuple[DWBAContribution, ...]

    @property
    def F_atomic(self):
        """Return the coherent sum of actual-density record amplitudes."""
        return sum(
            (contribution.F_atomic for contribution in self.contributions),
            start=0j,
        )

    @property
    def F_reference(self):
        """Return the coherent planar-reference amplitude that was subtracted."""
        return sum(
            (contribution.F_reference for contribution in self.contributions),
            start=0j,
        )

    @property
    def structure_factor_squared(self):
        """Return ``abs(F_contrast)**2`` in squared-electron units."""
        return np.abs(self.F_contrast) ** 2

    @property
    def scattered_amplitude_squared(self):
        """Return the squared modulus of the dimensionless scattered field."""
        return np.abs(self.scattered_amplitude) ** 2

    @property
    def differential_cross_section_kernel(self):
        """Return ``r_e**2 * abs(F_contrast)**2`` per reference cell.

        This is not integrated over footprint, coherence, detector acceptance,
        or instrumental resolution.
        """
        return (
            _CLASSICAL_ELECTRON_RADIUS_ANGSTROM**2
            * self.structure_factor_squared
        )

    def _require_specular_reflection(self):
        if self.bulk_mode != "semi_infinite":
            raise ValueError("Reflectivity requires a semi-infinite bulk result.")
        if not np.all(self.prepared.is_specular):
            raise ValueError("Reflectivity is defined only for specular results.")
        if self.prepared.polarization_i != self.prepared.polarization_f:
            raise ValueError(
                "Reflectivity requires identical incident and analyzed "
                "polarizations."
            )

    @property
    def reflectivity(self):
        """Return coherent specular reflectivity ``abs(total_amplitude)**2``."""
        self._require_specular_reflection()
        return np.abs(self.total_amplitude) ** 2

    @property
    def first_order_reflectivity(self):
        """Return reflectivity with ``abs(scattered_amplitude)**2`` omitted."""
        self._require_specular_reflection()
        return np.abs(self.unperturbed_amplitude) ** 2 + 2.0 * np.real(
            np.conj(self.unperturbed_amplitude) * self.scattered_amplitude
        )


@dataclass(frozen=True)
class _GeometryRule:
    fixed: str
    angle: float
    mirrorx: bool


@dataclass(frozen=True)
class _PreparedGeometry:
    shape: tuple
    scalar: bool
    hkl_reference: np.ndarray
    alpha_i: np.ndarray
    alpha_f: np.ndarray
    vlieg_angles: np.ndarray | None


class DWBAState:
    """Resolve, prepare, cache, and evaluate DWBA scans for one crystal."""

    def __init__(self, crystal):
        self._crystal = crystal
        self._orientation = _readonly_array(np.identity(3), np.float64)
        self._default_geometry = None
        self._rod_geometry = {}
        self._geometry_revision = 0
        self._geometry_cache = OrderedDict()
        self._reference_cache = OrderedDict()
        self._field_cache = OrderedDict()
        self._prepared_cache = OrderedDict()
        self._cache_counters = {
            "geometry_hits": 0,
            "geometry_misses": 0,
            "reference_hits": 0,
            "reference_misses": 0,
            "field_hits": 0,
            "field_misses": 0,
            "prepared_hits": 0,
            "prepared_misses": 0,
        }

    @property
    def orientation(self):
        """Return the configured 3-by-3 orientation matrix."""
        return self._orientation

    def set_orientation(self, U):
        """Set the orientation matrix used by Vlieg angle calculations.

        :param numpy.ndarray U:
            Finite, nonsingular 3-by-3 orientation matrix.
        """
        U = np.asarray(U, dtype=np.float64)
        if U.shape != (3, 3) or not np.all(np.isfinite(U)):
            raise ValueError("U must be a finite 3-by-3 matrix.")
        if abs(np.linalg.det(U)) <= np.finfo(np.float64).eps:
            raise ValueError("U must be nonsingular.")
        self._orientation = _readonly_array(U)
        self._geometry_revision += 1
        self._geometry_cache.clear()
        self._prepared_cache.clear()

    def set_ctr_geometry(
        self,
        *,
        alpha_i=None,
        alpha_f=None,
        equal_angles=False,
        rods=None,
        mirrorx=False,
    ):
        """Set a persistent Vlieg z-mode geometry rule.

        Exactly one of ``alpha_i``, ``alpha_f``, or ``equal_angles=True`` is
        required. Angles are in radians. ``rods=None`` sets the default rule
        for non-specular points; otherwise ``rods`` contains ``(h, k)`` pairs.
        """
        selected = int(alpha_i is not None) + int(alpha_f is not None) + int(
            bool(equal_angles)
        )
        if selected != 1:
            raise ValueError(
                "Set exactly one of alpha_i, alpha_f, or equal_angles=True."
            )
        if equal_angles:
            rule = _GeometryRule("eq", 0.0, bool(mirrorx))
        else:
            angle = float(alpha_i if alpha_i is not None else alpha_f)
            if not np.isfinite(angle) or angle <= 0.0 or angle > 0.5 * np.pi:
                raise ValueError("The fixed glancing angle must be in (0, pi/2].")
            rule = _GeometryRule(
                "in" if alpha_i is not None else "out",
                angle,
                bool(mirrorx),
            )
        if rods is None:
            self._default_geometry = rule
        else:
            normalized = []
            for rod in rods:
                if len(rod) != 2:
                    raise ValueError("Each rod must be an (h, k) pair.")
                key = (float(rod[0]), float(rod[1]))
                if not np.all(np.isfinite(key)):
                    raise ValueError("Rod coordinates must be finite.")
                normalized.append(key)
            if not normalized:
                raise ValueError("rods must not be empty.")
            for key in normalized:
                self._rod_geometry[key] = rule
        self._geometry_revision += 1
        self._geometry_cache.clear()
        self._prepared_cache.clear()

    def clear_cache(self):
        """Clear all automatic DWBA geometry, reference, field, and scan caches."""
        self._geometry_cache.clear()
        self._reference_cache.clear()
        self._field_cache.clear()
        self._prepared_cache.clear()

    def _energy_changed(self):
        self.clear_cache()

    def cache_info(self):
        """Return DWBA cache sizes and cumulative hit/miss counters."""
        return {
            **self._cache_counters,
            "geometry_size": len(self._geometry_cache),
            "reference_size": len(self._reference_cache),
            "field_size": len(self._field_cache),
            "prepared_size": len(self._prepared_cache),
            "prepared_capacity": _PREPARATION_CACHE_SIZE,
            "field_capacity": _FIELD_CACHE_SIZE,
        }

    def _require_available(self):
        if not HAS_CPP_ACCEL or not hasattr(
            _CTRcalc_cpp, "unitcell_F_DWBA_records"
        ):
            raise RuntimeError(
                "DWBA requires the native CTR extension with DWBA support."
            )
        if not hasattr(self._crystal.uc_bulk, "_E"):
            raise ValueError("Set the bulk unit-cell energy before using DWBA.")
        for component in self._crystal.uc_surface_list:
            if isinstance(component, WaterModel):
                raise NotImplementedError(
                    "DWBA does not support continuous-density WaterModel "
                    "components."
                )
            if not isinstance(
                component, UnitCell | Film | EpitaxyInterface | PoissonSurface
            ):
                raise NotImplementedError(
                    "DWBA supports only UnitCell, Film, EpitaxyInterface, "
                    f"and PoissonSurface components, not {type(component).__name__}."
                )

    def _live_records(self):
        self._require_available()
        if self._crystal.enable_uc_stacking:
            self._crystal.apply_stacking()
        bulk = self._crystal.uc_bulk
        records = [
            _LiveAtomicRecord(
                component_index=-1,
                component_name=bulk.name,
                component_type=type(bulk).__name__,
                record_index=0,
                record_name=bulk.name,
                role="bulk",
                layer=None,
                cell=bulk,
                component_weight=1.0,
                outer_domains=((np.identity(3, dtype=np.float64), 1.0),),
                identity_token=("bulk", id(bulk)),
            )
        ]
        for component_index, (component, component_weight, outer_domains) in enumerate(
            zip(
                self._crystal.uc_surface_list,
                self._crystal.weights,
                self._crystal.domains,
            )
        ):
            if not np.isfinite(component_weight):
                raise ValueError("DWBA component weights must be finite.")
            outer_domains = tuple(outer_domains)
            for record_index, (role, cell) in enumerate(
                _component_record_cells(component)
            ):
                records.append(
                    _LiveAtomicRecord(
                        component_index=component_index,
                        component_name=component.name,
                        component_type=type(component).__name__,
                        record_index=record_index,
                        record_name=cell.name,
                        role=role,
                        layer=_record_layer(cell),
                        cell=cell,
                        component_weight=float(component_weight),
                        outer_domains=outer_domains,
                        identity_token=(
                            "surface",
                            component_index,
                            id(component),
                            role,
                            record_index,
                            _record_layer(cell),
                        ),
                    )
                )
        return tuple(records)

    def _tagged_component_profile(self, component, records, record_count):
        values = []
        shares = []
        for global_index, record in records:
            profile = record.cell.optical_profile()
            if not len(profile):
                continue
            for outer_matrix, outer_occupancy in record.outer_domains:
                if not np.isfinite(outer_occupancy):
                    raise ValueError(
                        "DWBA crystal-domain occupancies must be finite."
                    )
                position_transform = self._live_position_transform(
                    record.cell, outer_matrix
                )
                normal_stretch = float(position_transform[2, 2])
                scale = (
                    record.component_weight
                    * float(outer_occupancy)
                    / normal_stretch
                )
                transformed = np.array(profile, copy=True)
                transformed[:, 0] *= normal_stretch
                transformed[:, 1:] *= scale
                tagged = np.zeros(
                    (record_count, len(profile), 2), dtype=np.float64
                )
                tagged[global_index] = transformed[:, 1:]
                values.append(transformed)
                shares.append(tagged)
        if not values:
            return (
                np.empty((0, 3), dtype=np.float64),
                np.empty((record_count, 0, 2), dtype=np.float64),
            )
        values = np.concatenate(values, axis=0)
        shares = np.concatenate(shares, axis=1)
        if not isinstance(component, UnitCell):
            values, shares = _combine_tagged(values, shares)
        return values, shares

    def _atomic_model_snapshot(self):
        records = self._live_records()
        record_count = len(records)
        if not self._crystal.uc_surface_list:
            profile = homogeneous_bulk_profile(self._crystal.uc_bulk)
            shares = np.zeros((1, 2, 2), dtype=np.float64)
            shares[0, 0] = profile.values[0, 1:]
            stratified = profile
        else:
            bulk_profile = self._crystal.uc_bulk.optical_profile_asbulk()
            bulk_shares = np.zeros(
                (record_count, len(bulk_profile), 2), dtype=np.float64
            )
            bulk_shares[0] = bulk_profile[:, 1:]
            component_values = [bulk_profile]
            component_shares = [bulk_shares]
            for component_index, component in enumerate(
                self._crystal.uc_surface_list
            ):
                selected = [
                    (global_index, record)
                    for global_index, record in enumerate(records)
                    if record.component_index == component_index
                ]
                values, tagged = self._tagged_component_profile(
                    component, selected, record_count
                )
                if len(values):
                    component_values.append(values)
                    component_shares.append(tagged)
            structural, shares = _combine_tagged(
                np.concatenate(component_values, axis=0),
                np.concatenate(component_shares, axis=1),
            )
            dz = top_layer_spacing(structural)
            vacuum = np.array(
                [[structural[-1, 0] + dz, 0.0, 0.0]], dtype=np.float64
            )
            structural = np.concatenate((structural, vacuum), axis=0)
            shares = np.concatenate(
                (shares, np.zeros((record_count, 1, 2), dtype=np.float64)),
                axis=1,
            )
            stratified, shares = _stratify_tagged(
                structural,
                shares,
                _OPTICAL_DELTA_TOLERANCE,
                _OPTICAL_BETA_TOLERANCE,
            )
        if not np.allclose(
            shares.sum(axis=0), stratified.values[:, 1:], rtol=2e-14, atol=1e-18
        ):
            raise RuntimeError(
                "DWBA record reference shares do not reconstruct the optical profile."
            )
        geometry_fingerprint = _hash_record_geometry(records)
        digest = hashlib.blake2b(digest_size=20)
        digest.update(geometry_fingerprint)
        digest.update(_hash_arrays(shares))
        for record in records:
            digest.update(repr(_record_descriptor(record)).encode("utf-8"))
        model = _PreparedAtomicModel(
            descriptors=tuple(_record_descriptor(record) for record in records),
            reference_shares=_readonly_array(shares[:, ::-1], np.float64),
            geometry_fingerprint=geometry_fingerprint,
            decomposition_fingerprint=digest.digest(),
        )
        return records, model, stratified

    def _vlieg_angles(self):
        energy_eV = float(self._crystal.uc_bulk._E)
        ub = UBCalculator(self._crystal.uc_bulk, energy_eV * 1e-3)
        ub.setLambda(HC_KEV_ANGSTROM / (energy_eV * 1e-3))
        ub.setU(self._orientation)
        return VliegAngles(ub)

    def _reference(self, profile):
        values = _readonly_array(profile.values, np.float64)
        boundaries = _readonly_array(profile.boundaries, np.float64)
        energy_eV = float(self._crystal.uc_bulk._E)
        energy_bytes = np.asarray([energy_eV], dtype=np.float64)
        simplification = np.asarray(
            [_OPTICAL_DELTA_TOLERANCE, _OPTICAL_BETA_TOLERANCE],
            dtype=np.float64,
        )
        fingerprint = _hash_arrays(
            values, boundaries, energy_bytes, simplification
        )
        cached = self._reference_cache.get(fingerprint)
        if cached is not None:
            self._cache_counters["reference_hits"] += 1
            self._reference_cache.move_to_end(fingerprint)
            return cached
        self._cache_counters["reference_misses"] += 1
        wavelength = HC_KEV_ANGSTROM / (energy_eV * 1e-3)
        reference = PreparedOpticalReference(
            profile=StratifiedProfile(values, boundaries),
            n=_readonly_array(
                (1.0 - values[:, 1] - 1j * values[:, 2])[::-1],
                np.complex128,
            ),
            energy_eV=energy_eV,
            k0=2.0 * np.pi / wavelength,
            fingerprint=fingerprint,
        )
        self._reference_cache[fingerprint] = reference
        while len(self._reference_cache) > _PREPARATION_CACHE_SIZE:
            self._reference_cache.popitem(last=False)
        return reference

    def _field_cache_key(self, reference, polarization, angle):
        return (
            reference.fingerprint,
            polarization,
            np.asarray(float(angle), dtype=np.float64).tobytes(),
        )

    @staticmethod
    def _slice_field(field, index):
        return LayeredElectricField(
            z_interfaces=field.z_interfaces,
            z_reference=field.z_reference,
            n=field.n,
            kz=_readonly_array(np.asarray(field.kz)[:, index]),
            A_plus=_readonly_array(np.asarray(field.A_plus)[:, index]),
            A_minus=_readonly_array(np.asarray(field.A_minus)[:, index]),
            r_S=complex(np.asarray(field.r_S)[index]),
            t_S=complex(np.asarray(field.t_S)[index]),
            polarization=field.polarization,
            angle_shape=(),
            _A_plus_over_kz=_readonly_array(
                np.asarray(field._A_plus_over_kz)[:, index]
            ),
            _A_minus_over_kz=_readonly_array(
                np.asarray(field._A_minus_over_kz)[:, index]
            ),
        )

    @staticmethod
    def _assemble_fields(fields, shape):
        first = fields[0]
        def stack(name):
            return _readonly_array(
                np.stack([np.asarray(getattr(field, name)) for field in fields], axis=1)
                .reshape((len(first.n),) + shape)
            )

        return LayeredElectricField(
            z_interfaces=first.z_interfaces,
            z_reference=first.z_reference,
            n=first.n,
            kz=stack("kz"),
            A_plus=stack("A_plus"),
            A_minus=stack("A_minus"),
            r_S=_readonly_array(
                np.asarray([field.r_S for field in fields], dtype=np.complex128)
                .reshape(shape)
            ),
            t_S=_readonly_array(
                np.asarray([field.t_S for field in fields], dtype=np.complex128)
                .reshape(shape)
            ),
            polarization=first.polarization,
            angle_shape=shape,
            _A_plus_over_kz=stack("_A_plus_over_kz"),
            _A_minus_over_kz=stack("_A_minus_over_kz"),
        )

    def _fields(self, reference, angles, polarization):
        flat = np.asarray(angles, dtype=np.float64).reshape(-1)
        _, first_index, inverse_sorted = np.unique(
            flat.view(np.uint64), return_index=True, return_inverse=True
        )
        order = np.argsort(first_index)
        unique_values = flat[first_index[order]]
        sorted_to_order = np.empty(order.size, dtype=np.int64)
        sorted_to_order[order] = np.arange(order.size, dtype=np.int64)
        point_index = sorted_to_order[inverse_sorted]
        missing_values = []
        missing_keys = []
        for angle in unique_values:
            key = self._field_cache_key(reference, polarization, angle)
            if key in self._field_cache:
                self._cache_counters["field_hits"] += 1
                self._field_cache.move_to_end(key)
            else:
                missing_keys.append(key)
                missing_values.append(float(angle))
        if missing_values:
            self._cache_counters["field_misses"] += len(missing_values)
            solved = solve_electric_field(
                reference.profile.values,
                reference.energy_eV,
                np.asarray(missing_values),
                polarization,
                boundaries=reference.profile.boundaries,
            )
            if not np.array_equal(solved.n, reference.n):
                raise RuntimeError("The solved field does not use its frozen n0.")
            for index, key in enumerate(missing_keys):
                self._field_cache[key] = self._slice_field(solved, index)
            while len(self._field_cache) > _FIELD_CACHE_SIZE:
                self._field_cache.popitem(last=False)
        fields = [
            self._field_cache[
                self._field_cache_key(reference, polarization, angle)
            ]
            for angle in unique_values
        ]
        return (
            self._assemble_fields(fields, (len(fields),)),
            _readonly_array(point_index, np.int64),
        )

    def _geometry_key(self, h, k, l, geometry_rtol, geometry_atol):  # noqa: E741
        return (
            _hash_arrays(h, k, l, self._orientation),
            float(self._crystal.uc_bulk._E),
            self._geometry_revision,
            float(geometry_rtol),
            float(geometry_atol),
            _hash_arrays(
                self._crystal.uc_bulk.B_mat,
                self._crystal.uc_bulk.refHKLTransform,
            ),
        )

    def _rule_for_constant_rod(self, h, k):
        for rod, rule in self._rod_geometry.items():
            if h == rod[0] and k == rod[1]:
                return rule
        return self._default_geometry

    def _geometry_from_hkl(  # noqa: E741
        self, h, k, l, geometry_rtol, geometry_atol  # noqa: E741
    ):
        arrays = np.broadcast_arrays(
            np.asarray(h, dtype=np.float64),
            np.asarray(k, dtype=np.float64),
            np.asarray(l, dtype=np.float64),
        )
        if arrays[0].size == 0 or not all(
            np.all(np.isfinite(array)) for array in arrays
        ):
            raise ValueError("DWBA coordinates must be finite and nonempty.")
        key = self._geometry_key(*arrays, geometry_rtol, geometry_atol)
        cached = self._geometry_cache.get(key)
        if cached is not None:
            self._cache_counters["geometry_hits"] += 1
            self._geometry_cache.move_to_end(key)
            return cached
        self._cache_counters["geometry_misses"] += 1
        shape = arrays[0].shape
        scalar = arrays[0].ndim == 0
        hkl_reference = np.vstack([array.ravel() for array in arrays])
        hkl_bulk = self._crystal.uc_bulk.refHKLTransform @ hkl_reference
        specular = np.logical_and(
            np.isclose(hkl_reference[0], 0.0, rtol=0.0, atol=1e-14),
            np.isclose(hkl_reference[1], 0.0, rtol=0.0, atol=1e-14),
        )
        groups = []
        if np.any(specular):
            groups.append((specular, _GeometryRule("eq", 0.0, False)))
        nonspecular = ~specular
        if np.any(nonspecular):
            unique_rods = np.unique(
                hkl_reference[:2, nonspecular].T, axis=0
            )
            resolved_rules = {
                self._rule_for_constant_rod(float(h_value), float(k_value))
                for h_value, k_value in unique_rods
            }
            if None in resolved_rules:
                raise ValueError(
                    "No CTR geometry is configured for one or more "
                    "non-specular rods."
                )
            if len(resolved_rules) != 1:
                raise ValueError(
                    "Mixed hkl arrays requiring conflicting rod-specific "
                    "geometry rules must be split by rod."
                )
            groups.append((nonspecular, resolved_rules.pop()))
        angles = np.empty((hkl_reference.shape[1], 6), dtype=np.float64)
        calculator = self._vlieg_angles()
        for mask, rule in groups:
            with np.errstate(invalid="ignore"):
                solved = calculator.anglesZmode(
                    hkl_bulk[:, mask],
                    rule.angle,
                    fixed=rule.fixed,
                    chi=0.0,
                    phi=0.0,
                    mirrorx=rule.mirrorx,
                )
            angles[mask] = np.asarray(solved, dtype=np.float64).reshape(-1, 6)
        invalid = ~np.all(np.isfinite(angles), axis=1)
        if np.any(invalid):
            index = int(np.flatnonzero(invalid)[0])
            raise ValueError(
                "The configured CTR geometry has no solution for hkl="
                f"{_format_hkl(hkl_reference, index)}."
            )
        prepared = _PreparedGeometry(
            shape=shape,
            scalar=scalar,
            hkl_reference=_readonly_array(hkl_reference),
            alpha_i=_readonly_array(angles[:, 0]),
            alpha_f=_readonly_array(angles[:, 2]),
            vlieg_angles=_readonly_array(angles.reshape(shape + (6,))),
        )
        self._geometry_cache[key] = prepared
        while len(self._geometry_cache) > _PREPARATION_CACHE_SIZE:
            self._geometry_cache.popitem(last=False)
        return prepared

    def _prepare_explicit(
        self,
        geometry,
        reference,
        atomic_model,
        polarization_i,
        polarization_f,
        geometry_rtol,
        geometry_atol,
    ):
        if polarization_i not in {"s", "p"}:
            raise ValueError("polarization_i must be 's' or 'p'.")
        if polarization_f not in {"s", "p"}:
            raise ValueError("polarization_f must be 's' or 'p'.")
        hkl_reference = geometry.hkl_reference
        alpha_i = np.asarray(geometry.alpha_i).reshape(-1)
        alpha_f = np.asarray(geometry.alpha_f).reshape(-1)
        if np.any(alpha_i <= 0.0) or np.any(alpha_i > 0.5 * np.pi):
            raise ValueError("alpha_i must lie in (0, pi/2] radians.")
        if np.any(alpha_f <= 0.0) or np.any(alpha_f > 0.5 * np.pi):
            raise ValueError("alpha_f must lie in (0, pi/2] radians.")
        hkl_bulk = self._crystal.uc_bulk.refHKLTransform @ hkl_reference
        Q_cartesian = self._orientation @ self._crystal.uc_bulk.B_mat @ hkl_bulk
        Q_parallel = np.ascontiguousarray(Q_cartesian[:2].T)
        k0 = reference.k0
        expected_Qz = k0 * (np.sin(alpha_i) + np.sin(alpha_f))
        if not np.allclose(
            Q_cartesian[2], expected_Qz, rtol=geometry_rtol, atol=geometry_atol
        ):
            maximum_error = float(np.max(np.abs(Q_cartesian[2] - expected_Qz)))
            raise ValueError(
                "The transformed hkl and glancing angles violate the Ewald "
                f"normal geometry (maximum |delta Qz|={maximum_error:.6g} "
                "Angstrom^-1)."
            )
        incident_parallel = k0 * np.cos(alpha_i)
        exit_parallel = k0 * np.cos(alpha_f)
        Q_parallel_squared = np.sum(Q_parallel**2, axis=1)
        denominator = 2.0 * incident_parallel * exit_parallel
        cos_azimuth = np.ones_like(denominator)
        regular = denominator > np.finfo(np.float64).eps * k0**2
        cos_azimuth[regular] = (
            incident_parallel[regular] ** 2
            + exit_parallel[regular] ** 2
            - Q_parallel_squared[regular]
        ) / denominator[regular]
        geometry_scale = np.maximum(
            np.maximum(incident_parallel, exit_parallel),
            np.sqrt(Q_parallel_squared),
        )
        tolerance = geometry_rtol + np.divide(
            geometry_atol,
            geometry_scale,
            out=np.full_like(geometry_scale, geometry_atol),
            where=geometry_scale > 0.0,
        )
        if np.any(np.abs(cos_azimuth[regular]) > 1.0 + tolerance[regular]):
            raise ValueError(
                "The in-plane scattering vector has no physical azimuth for "
                "the supplied angles."
            )
        singular = ~regular
        singular_error = np.abs(
            Q_parallel_squared - (incident_parallel - exit_parallel) ** 2
        )
        singular_tolerance = (geometry_atol + geometry_rtol * geometry_scale) ** 2
        if np.any(singular & (singular_error > singular_tolerance)):
            raise ValueError("The normal-incidence in-plane geometry is invalid.")
        cos_azimuth = np.clip(cos_azimuth, -1.0, 1.0)
        sin_azimuth = np.sqrt(np.maximum(0.0, 1.0 - cos_azimuth**2))
        field_i, field_i_index = self._fields(
            reference, alpha_i, polarization_i
        )
        field_f, field_f_index = self._fields(
            reference, alpha_f, polarization_f
        )
        if field_i.n is not reference.n and not np.array_equal(
            field_i.n, reference.n
        ):
            raise RuntimeError("The incident field does not share the frozen n0.")
        if field_f.n is not reference.n and not np.array_equal(
            field_f.n, reference.n
        ):
            raise RuntimeError("The exit field does not share the frozen n0.")
        if np.any(np.asarray(field_i.A_minus)[-1] != 0.0) or np.any(
            np.asarray(field_f.A_minus)[-1] != 0.0
        ):
            raise RuntimeError(
                "The terminal substrate must satisfy A_i_minus=A_f_minus=0."
            )
        specular_tolerance = 64.0 * np.finfo(np.float64).eps * max(k0, 1.0)
        is_specular = np.linalg.norm(Q_parallel, axis=1) <= specular_tolerance
        return PreparedCTR(
            shape=geometry.shape,
            scalar=geometry.scalar,
            hkl=_readonly_array(hkl_reference.reshape((3,) + geometry.shape)),
            vlieg_angles=geometry.vlieg_angles,
            Q_parallel=_readonly_array(Q_parallel),
            is_specular=_readonly_array(is_specular, np.int64),
            alpha_i=_readonly_array(alpha_i),
            alpha_f=_readonly_array(alpha_f),
            cos_azimuth=_readonly_array(cos_azimuth),
            sin_azimuth=_readonly_array(sin_azimuth),
            polarization_i=polarization_i,
            polarization_f=polarization_f,
            reference=reference,
            field_i=field_i,
            field_f=field_f,
            field_i_index=field_i_index,
            field_f_index=field_f_index,
            crystal_identity=id(self._crystal),
            bulk_identity=id(self._crystal.uc_bulk),
            energy_eV=reference.energy_eV,
            B_mat=_readonly_array(self._crystal.uc_bulk.B_mat, np.float64),
            R_mat=_readonly_array(self._crystal.uc_bulk.R_mat, np.float64),
            R_mat_inv=_readonly_array(self._crystal.uc_bulk.R_mat_inv, np.float64),
            ref_hkl_transform=_readonly_array(
                self._crystal.uc_bulk.refHKLTransform, np.float64
            ),
            orientation=self._orientation,
            reference_area=float(self._crystal.reference_area),
            _atomic_model=atomic_model,
        )

    def _prepared_key(
        self,
        geometry,
        reference,
        atomic_model,
        polarization_i,
        polarization_f,
        geometry_rtol,
        geometry_atol,
    ):
        return (
            _hash_arrays(
                geometry.hkl_reference,
                geometry.alpha_i,
                geometry.alpha_f,
            ),
            _hash_arrays(
                self._crystal.uc_bulk.B_mat,
                self._crystal.uc_bulk.R_mat,
                self._crystal.uc_bulk.R_mat_inv,
                self._crystal.uc_bulk.refHKLTransform,
                self._orientation,
                np.asarray([self._crystal.reference_area]),
            ),
            reference.fingerprint,
            atomic_model.decomposition_fingerprint,
            polarization_i,
            polarization_f,
            float(geometry_rtol),
            float(geometry_atol),
            self._geometry_revision,
        )

    def _prepare_geometry(
        self,
        geometry,
        *,
        polarization_i,
        polarization_f,
        geometry_rtol,
        geometry_atol,
        cache,
    ):
        _, atomic_model, profile = self._atomic_model_snapshot()
        reference = self._reference(profile)
        key = self._prepared_key(
            geometry,
            reference,
            atomic_model,
            polarization_i,
            polarization_f,
            geometry_rtol,
            geometry_atol,
        )
        if cache and key in self._prepared_cache:
            self._cache_counters["prepared_hits"] += 1
            self._prepared_cache.move_to_end(key)
            return self._prepared_cache[key]
        self._cache_counters["prepared_misses"] += 1
        prepared = self._prepare_explicit(
            geometry,
            reference,
            atomic_model,
            polarization_i,
            polarization_f,
            geometry_rtol,
            geometry_atol,
        )
        if cache:
            self._prepared_cache[key] = prepared
            while len(self._prepared_cache) > _PREPARATION_CACHE_SIZE:
                self._prepared_cache.popitem(last=False)
        return prepared

    def prepare(
        self,
        h,
        k,
        l,  # noqa: E741
        *,
        polarization_i="s",
        polarization_f="s",
        geometry_rtol=1e-7,
        geometry_atol=1e-8,
        cache=True,
    ):
        """Prepare and optionally cache DWBA state resolved from hkl.

        The configured CTR geometry supplies the glancing angles. Inputs
        broadcast, and ``h``, ``k``, and ``l`` are reference-cell r.l.u.

        :param h: Reference in-plane index in r.l.u.
        :param k: Reference in-plane index in r.l.u.
        :param l: Reference out-of-plane index in r.l.u.
        :param str polarization_i: Incident ``"s"`` or ``"p"`` polarization.
        :param str polarization_f: Reciprocal-exit polarization.
        :param bool cache: Reuse and retain the exact preparation when true.
        :returns: Immutable geometry, optical reference, and field snapshot.
        :rtype: PreparedCTR
        """
        self._require_available()
        geometry = self._geometry_from_hkl(
            h, k, l, geometry_rtol, geometry_atol
        )
        return self._prepare_geometry(
            geometry,
            polarization_i=polarization_i,
            polarization_f=polarization_f,
            geometry_rtol=geometry_rtol,
            geometry_atol=geometry_atol,
            cache=bool(cache),
        )

    def _l_from_glancing(self, h, k, alpha_i, alpha_f):
        arrays = np.broadcast_arrays(
            np.asarray(h, dtype=np.float64),
            np.asarray(k, dtype=np.float64),
            np.asarray(alpha_i, dtype=np.float64),
            np.asarray(alpha_f, dtype=np.float64),
        )
        if arrays[0].size == 0 or not all(
            np.all(np.isfinite(array)) for array in arrays
        ):
            raise ValueError("Coordinates and glancing angles must be finite.")
        transform = (
            self._orientation
            @ self._crystal.uc_bulk.B_mat
            @ self._crystal.uc_bulk.refHKLTransform
        )
        normal = np.asarray(transform[2], dtype=np.float64)
        if abs(normal[2]) <= np.finfo(np.float64).eps * np.linalg.norm(normal):
            raise ValueError(
                "The reference out-of-plane axis has no surface-normal component."
            )
        energy_eV = float(self._crystal.uc_bulk._E)
        k0 = 2.0 * np.pi / (HC_KEV_ANGSTROM / (energy_eV * 1e-3))
        h_array, k_array, ai, af = arrays
        Qz = k0 * (np.sin(ai) + np.sin(af))
        return (Qz - normal[0] * h_array - normal[1] * k_array) / normal[2]

    def prepare_from_glancing(
        self,
        h,
        k,
        alpha_i,
        alpha_f,
        *,
        polarization_i="s",
        polarization_f="s",
        geometry_rtol=1e-7,
        geometry_atol=1e-8,
        cache=True,
    ):
        """Prepare from h, k and measured incident and exit glancing angles.

        The out-of-plane index is derived from the Ewald normal condition.
        All angles are radians; coordinates and angles broadcast together.

        :param h: Reference in-plane index in r.l.u.
        :param k: Reference in-plane index in r.l.u.
        :param alpha_i: Incident glancing angle in radians.
        :param alpha_f: Exit glancing angle in radians.
        :returns: Immutable preparation whose ``hkl[2]`` contains derived l.
        :rtype: PreparedCTR
        """
        l_value = self._l_from_glancing(h, k, alpha_i, alpha_f)
        arrays = np.broadcast_arrays(
            np.asarray(h, dtype=np.float64),
            np.asarray(k, dtype=np.float64),
            np.asarray(l_value, dtype=np.float64),
            np.asarray(alpha_i, dtype=np.float64),
            np.asarray(alpha_f, dtype=np.float64),
        )
        geometry = _PreparedGeometry(
            shape=arrays[0].shape,
            scalar=arrays[0].ndim == 0,
            hkl_reference=_readonly_array(
                np.vstack([array.ravel() for array in arrays[:3]])
            ),
            alpha_i=_readonly_array(arrays[3].ravel()),
            alpha_f=_readonly_array(arrays[4].ravel()),
            vlieg_angles=None,
        )
        return self._prepare_geometry(
            geometry,
            polarization_i=polarization_i,
            polarization_f=polarization_f,
            geometry_rtol=geometry_rtol,
            geometry_atol=geometry_atol,
            cache=bool(cache),
        )

    def prepare_from_vlieg(
        self,
        alpha,
        delta,
        gamma,
        omega,
        chi,
        phi,
        *,
        polarization_i="s",
        polarization_f="s",
        geometry_rtol=1e-7,
        geometry_atol=1e-8,
        cache=True,
    ):
        """Prepare from measured six-circle Vlieg angles.

        All six angles are radians and broadcast together. The configured
        orientation converts them to reference-cell hkl.

        :param alpha: Vlieg incidence circle in radians.
        :param delta: Vlieg detector delta circle in radians.
        :param gamma: Vlieg exit circle in radians.
        :param omega: Vlieg sample omega circle in radians.
        :param chi: Vlieg sample chi circle in radians.
        :param phi: Vlieg sample phi circle in radians.
        :returns: Immutable measured-geometry preparation.
        :rtype: PreparedCTR
        """
        arrays = np.broadcast_arrays(
            *[
                np.asarray(value, dtype=np.float64)
                for value in (alpha, delta, gamma, omega, chi, phi)
            ]
        )
        if arrays[0].size == 0 or not all(
            np.all(np.isfinite(array)) for array in arrays
        ):
            raise ValueError("Vlieg angles must be finite and nonempty.")
        shape = arrays[0].shape
        bulk_hkl = np.vstack(
            self._vlieg_angles().anglesToHkl(
                *[array.ravel() for array in arrays]
            )
        )
        hkl_reference = np.linalg.solve(
            self._crystal.uc_bulk.refHKLTransform, bulk_hkl
        )
        vlieg = np.stack(arrays, axis=-1)
        geometry = _PreparedGeometry(
            shape=shape,
            scalar=arrays[0].ndim == 0,
            hkl_reference=_readonly_array(hkl_reference),
            alpha_i=_readonly_array(arrays[0].ravel()),
            alpha_f=_readonly_array(arrays[2].ravel()),
            vlieg_angles=_readonly_array(vlieg),
        )
        return self._prepare_geometry(
            geometry,
            polarization_i=polarization_i,
            polarization_f=polarization_f,
            geometry_rtol=geometry_rtol,
            geometry_atol=geometry_atol,
            cache=bool(cache),
        )

    def _position_transform(self, prepared, cell, outer_matrix):
        reference_q = (
            prepared.orientation
            @ prepared.B_mat
            @ prepared.ref_hkl_transform
        )
        return self._position_transform_from_reference(
            reference_q, cell, outer_matrix
        )

    def _live_position_transform(self, cell, outer_matrix):
        bulk = self._crystal.uc_bulk
        reference_q = (
            self._orientation @ bulk.B_mat @ bulk.refHKLTransform
        )
        return self._position_transform_from_reference(
            reference_q, cell, outer_matrix
        )

    def _position_transform_from_reference(
        self, reference_q, cell, outer_matrix
    ):
        outer_matrix = np.asarray(outer_matrix, dtype=np.float64)
        if outer_matrix.shape != (3, 3) or not np.all(np.isfinite(outer_matrix)):
            raise ValueError("Crystal domain matrices must be finite 3-by-3 arrays.")
        q_transform = (
            cell.B_mat
            @ cell.refHKLTransform
            @ outer_matrix
            @ np.linalg.inv(reference_q)
        )
        position_transform = q_transform.T
        if not np.all(np.isfinite(position_transform)):
            raise ValueError("DWBA record coordinate transforms must be finite.")
        self._validate_planar_transform(position_transform)
        return position_transform

    @staticmethod
    def _validate_planar_transform(position_transform):
        """Validate one physical Cartesian transform for layered optics."""
        position_transform = np.asarray(position_transform, dtype=np.float64)
        tolerance = 2e-10 * max(1.0, np.linalg.norm(position_transform))
        if (
            np.any(np.abs(position_transform[2, :2]) > tolerance)
            or position_transform[2, 2] <= tolerance
        ):
            raise NotImplementedError(
                "DWBA crystal domains must preserve the planar surface normal."
            )
        lateral = position_transform[:2, :2]
        projected_jacobian = abs(np.linalg.det(lateral))
        if not np.isclose(projected_jacobian, 1.0, rtol=2e-10, atol=2e-12):
            raise NotImplementedError(
                "DWBA does not yet support transformed lateral-area Jacobians."
            )
        if not np.allclose(
            lateral.T @ lateral,
            np.identity(2),
            rtol=2e-10,
            atol=2e-12,
        ):
            raise NotImplementedError(
                "DWBA does not yet support in-plane strain transforms."
            )

    def _pack_atomic_records(self, prepared, records):
        atom_basis = []
        form_factors = []
        finite_positions = []
        finite_atom_index = []
        finite_record_index = []
        finite_weight = []
        bulk_positions = []
        bulk_atom_index = []
        bulk_weight = []
        bulk_repeat_coordinate = []
        bulk_repeat = None

        for global_record_index, record in enumerate(records):
            cell = record.cell
            basis, factors, _ = cell.build_selected_basis()
            if cell._special_formfactors_present:
                raise NotImplementedError(
                    "DWBA does not support Python special form-factor callbacks."
                )
            basis = np.asarray(basis, dtype=np.float64)
            factors = np.asarray(factors, dtype=np.float64)
            atom_offset = sum(len(value) for value in atom_basis)
            atom_basis.append(basis)
            form_factors.append(factors)
            atom_indices = np.arange(len(basis), dtype=np.int64) + atom_offset
            domain_matrices, domain_occupancies = _coherent_domain_arrays(
                cell.coherentDomainMatrix,
                cell.coherentDomainOccupancy,
            )
            if len(domain_matrices) != len(domain_occupancies):
                raise ValueError(
                    f"UnitCell {cell.name!r} has inconsistent coherent domains."
                )
            if not np.all(np.isfinite(domain_matrices)) or not np.all(
                np.isfinite(domain_occupancies)
            ):
                raise ValueError("DWBA coherent-domain values must be finite.")
            for outer_matrix, outer_occupancy in record.outer_domains:
                if not np.isfinite(outer_occupancy):
                    raise ValueError("DWBA crystal-domain occupancies must be finite.")
                position_transform = self._position_transform(
                    prepared, cell, outer_matrix
                )
                if record.component_index == -1:
                    repeat = position_transform @ np.asarray(cell.R_mat[:, 2])
                    if bulk_repeat is None:
                        bulk_repeat = repeat
                    elif not np.allclose(bulk_repeat, repeat, rtol=2e-12, atol=2e-13):
                        raise ValueError(
                            "All bulk domains must use one physical repeat vector."
                        )
                for matrix, occupancy in zip(
                    domain_matrices, domain_occupancies
                ):
                    effective = cell.R_mat_inv @ matrix[:, :3] @ cell.R_mat
                    self._validate_planar_transform(
                        position_transform @ matrix[:, :3]
                    )
                    fractional = basis[:, 1:4] @ effective.T + matrix[:, 3]
                    cell_cartesian = fractional @ cell.R_mat.T
                    positions = cell_cartesian @ position_transform.T
                    scale = (
                        prepared.reference_area
                        / cell.uc_area
                        * record.component_weight
                        * float(outer_occupancy)
                        * float(occupancy)
                    )
                    if record.component_index == -1:
                        bulk_positions.append(positions)
                        bulk_atom_index.append(atom_indices)
                        bulk_weight.append(np.full(len(basis), scale))
                        bulk_repeat_coordinate.append(fractional[:, 2])
                    else:
                        finite_positions.append(positions)
                        finite_atom_index.append(atom_indices)
                        finite_record_index.append(
                            np.full(len(basis), global_record_index, dtype=np.int64)
                        )
                        finite_weight.append(np.full(len(basis), scale))

        def concatenate(values, shape, dtype):
            if not values:
                return np.empty(shape, dtype=dtype)
            return np.ascontiguousarray(np.concatenate(values, axis=0), dtype=dtype)

        if bulk_repeat is None:
            raise ValueError("The DWBA bulk must contain at least one domain.")
        return {
            "atom_basis": concatenate(atom_basis, (0, 8), np.float64),
            "form_factors": concatenate(form_factors, (0, 13), np.float64),
            "finite_positions": concatenate(
                finite_positions, (0, 3), np.float64
            ),
            "finite_atom_index": concatenate(
                finite_atom_index, (0,), np.int64
            ),
            "finite_record_index": concatenate(
                finite_record_index, (0,), np.int64
            ),
            "finite_weight": concatenate(finite_weight, (0,), np.float64),
            "bulk_positions": concatenate(bulk_positions, (0, 3), np.float64),
            "bulk_atom_index": concatenate(bulk_atom_index, (0,), np.int64),
            "bulk_weight": concatenate(bulk_weight, (0,), np.float64),
            "bulk_repeat_coordinate": concatenate(
                bulk_repeat_coordinate, (0,), np.float64
            ),
            "bulk_repeat": np.ascontiguousarray(bulk_repeat, dtype=np.float64),
        }

    def _validate_prepared(self, prepared):
        if not isinstance(prepared, PreparedCTR):
            raise TypeError("prepared must be a PreparedCTR instance.")
        if prepared.crystal_identity != id(self._crystal):
            raise ValueError("The prepared CTR belongs to another crystal.")
        if prepared.bulk_identity != id(self._crystal.uc_bulk):
            raise ValueError("The prepared CTR belongs to another bulk cell.")
        frozen_values = (
            ("energy", prepared.energy_eV, getattr(self._crystal.uc_bulk, "_E", None)),
            ("B_mat", prepared.B_mat, self._crystal.uc_bulk.B_mat),
            ("R_mat", prepared.R_mat, self._crystal.uc_bulk.R_mat),
            ("R_mat_inv", prepared.R_mat_inv, self._crystal.uc_bulk.R_mat_inv),
            (
                "refHKLTransform",
                prepared.ref_hkl_transform,
                self._crystal.uc_bulk.refHKLTransform,
            ),
            ("orientation", prepared.orientation, self._orientation),
            ("reference_area", prepared.reference_area, self._crystal.reference_area),
        )
        for name, frozen, current in frozen_values:
            if current is None or not np.array_equal(frozen, current):
                raise ValueError(
                    f"The crystal {name} changed; prepare the CTR again."
                )
        records = self._live_records()
        if tuple(_record_descriptor(record) for record in records) != (
            prepared._atomic_model.descriptors
        ) or _hash_record_geometry(records) != (
            prepared._atomic_model.geometry_fingerprint
        ):
            raise ValueError(
                "The crystal's generated DWBA record geometry changed; "
                "prepare the CTR again."
            )
        return records

    def evaluate_prepared(
        self,
        prepared,
        *,
        bulk_mode="semi_infinite",
        bulk_attenuation=0.0,
    ):
        """Evaluate current atomic parameters on a frozen DWBA preparation.

        No wavefield calculation occurs. ``bulk_attenuation`` is an optional
        empirical exponent per bulk repeat and defaults to zero because
        physical absorption is already present in the complex wavevectors.

        :param PreparedCTR prepared: Retained preparation from this crystal.
        :param str bulk_mode: ``"unit_cell"`` or ``"semi_infinite"``.
        :param float bulk_attenuation: Non-negative exponent per bulk repeat.
        :returns: Coherent contrast and normalized observable amplitudes.
        :rtype: DWBAResult
        """
        self._require_available()
        if bulk_mode not in {"unit_cell", "semi_infinite"}:
            raise ValueError("bulk_mode must be 'unit_cell' or 'semi_infinite'.")
        if not np.isfinite(bulk_attenuation) or bulk_attenuation < 0.0:
            raise ValueError("bulk_attenuation must be finite and non-negative.")
        records = self._validate_prepared(prepared)
        packed = self._pack_atomic_records(prepared, records)
        media = len(prepared.reference.n)

        def field_array(value):
            return np.ascontiguousarray(np.asarray(value).reshape(media, -1))

        atomic, reference = _CTRcalc_cpp.unitcell_F_DWBA_records(
            prepared.Q_parallel,
            prepared.is_specular,
            prepared.reference.n,
            prepared.field_i_index,
            prepared.field_f_index,
            field_array(prepared.field_i.kz),
            field_array(prepared.field_i.A_plus),
            field_array(prepared.field_i.A_minus),
            field_array(prepared.field_i._A_plus_over_kz),
            field_array(prepared.field_i._A_minus_over_kz),
            field_array(prepared.field_f.kz),
            field_array(prepared.field_f.A_plus),
            field_array(prepared.field_f.A_minus),
            field_array(prepared.field_f._A_plus_over_kz),
            field_array(prepared.field_f._A_minus_over_kz),
            prepared.field_i.z_reference,
            prepared.field_f.z_reference,
            prepared.field_i.z_interfaces,
            prepared.alpha_i,
            prepared.alpha_f,
            prepared.cos_azimuth,
            prepared.sin_azimuth,
            prepared.k0,
            prepared.polarization_i,
            prepared.polarization_f,
            prepared._atomic_model.reference_shares,
            prepared.reference_area,
            bulk_mode == "semi_infinite",
            float(bulk_attenuation),
            packed["atom_basis"],
            packed["form_factors"],
            packed["finite_positions"],
            packed["finite_atom_index"],
            packed["finite_record_index"],
            packed["finite_weight"],
            packed["bulk_positions"],
            packed["bulk_atom_index"],
            packed["bulk_weight"],
            packed["bulk_repeat_coordinate"],
            packed["bulk_repeat"],
        )
        atomic = np.asarray(atomic, dtype=np.complex128)
        reference = np.asarray(reference, dtype=np.complex128)
        contributions = []
        for descriptor, atomic_row, reference_row in zip(
            prepared._atomic_model.descriptors, atomic, reference
        ):
            contrast_row = atomic_row - reference_row
            contributions.append(
                DWBAContribution(
                    component_index=descriptor[0],
                    component_name=descriptor[1],
                    component_type=descriptor[2],
                    record_index=descriptor[3],
                    record_name=descriptor[4],
                    role=descriptor[5],
                    layer=descriptor[6],
                    F_atomic=_readonly_result(
                        atomic_row, prepared.shape, prepared.scalar, np.complex128
                    ),
                    F_reference=_readonly_result(
                        reference_row,
                        prepared.shape,
                        prepared.scalar,
                        np.complex128,
                    ),
                    F_contrast=_readonly_result(
                        contrast_row,
                        prepared.shape,
                        prepared.scalar,
                        np.complex128,
                    ),
                )
            )
        coherent_contrast = sum(
            (contribution.F_contrast for contribution in contributions),
            start=0j,
        )
        F_contrast = _readonly_result(
            coherent_contrast, prepared.shape, prepared.scalar, np.complex128
        )
        F_flat = np.asarray(F_contrast, dtype=np.complex128).reshape(-1)
        kappa_f = prepared.k0 * np.sin(prepared.alpha_f)
        scattered_flat = (
            2j
            * np.pi
            * _CLASSICAL_ELECTRON_RADIUS_ANGSTROM
            * F_flat
            / (kappa_f * prepared.reference_area)
        )
        reflection = np.asarray(
            prepared.field_i.r_S, dtype=np.complex128
        ).reshape(-1)[prepared.field_i_index]
        points = prepared.alpha_i.size
        unperturbed_flat = np.zeros(points, dtype=np.complex128)
        if (
            bulk_mode == "semi_infinite"
            and prepared.polarization_i == prepared.polarization_f
        ):
            mask = prepared.is_specular.astype(bool)
            unperturbed_flat[mask] = reflection[mask]
        total_flat = unperturbed_flat + scattered_flat
        return DWBAResult(
            prepared=prepared,
            F_contrast=F_contrast,
            unperturbed_amplitude=_readonly_result(
                unperturbed_flat, prepared.shape, prepared.scalar
            ),
            scattered_amplitude=_readonly_result(
                scattered_flat, prepared.shape, prepared.scalar
            ),
            total_amplitude=_readonly_result(
                total_flat, prepared.shape, prepared.scalar
            ),
            bulk_mode=bulk_mode,
            contributions=tuple(contributions),
        )

    def evaluate(self, h, k, l, **kwargs):  # noqa: E741
        """Prepare configured hkl if needed and return a DWBA result.

        Coordinates are reference-cell r.l.u. and broadcast as in
        :meth:`prepare`.

        :returns: Current atomic evaluation on the cached preparation.
        :rtype: DWBAResult
        """
        evaluation_keys = {"bulk_mode", "bulk_attenuation"}
        evaluation = {
            key: kwargs.pop(key) for key in tuple(kwargs) if key in evaluation_keys
        }
        prepared = self.prepare(h, k, l, **kwargs)
        return self.evaluate_prepared(prepared, **evaluation)

    def evaluate_from_glancing(self, h, k, alpha_i, alpha_f, **kwargs):
        """Prepare measured glancing geometry and return a DWBA result.

        Incident and exit glancing angles are radians.

        :rtype: DWBAResult
        """
        evaluation_keys = {"bulk_mode", "bulk_attenuation"}
        evaluation = {
            key: kwargs.pop(key) for key in tuple(kwargs) if key in evaluation_keys
        }
        prepared = self.prepare_from_glancing(
            h, k, alpha_i, alpha_f, **kwargs
        )
        return self.evaluate_prepared(prepared, **evaluation)

    def evaluate_from_vlieg(
        self, alpha, delta, gamma, omega, chi, phi, **kwargs
    ):
        """Prepare measured Vlieg geometry and return a DWBA result.

        All six Vlieg angles are radians.

        :rtype: DWBAResult
        """
        evaluation_keys = {"bulk_mode", "bulk_attenuation"}
        evaluation = {
            key: kwargs.pop(key) for key in tuple(kwargs) if key in evaluation_keys
        }
        prepared = self.prepare_from_vlieg(
            alpha, delta, gamma, omega, chi, phi, **kwargs
        )
        return self.evaluate_prepared(prepared, **evaluation)

    def reflectivity(
        self,
        h,
        k,
        l,  # noqa: E741
        *,
        polarization="s",
        first_order=False,
        **kwargs,
    ):
        """Return specular coherent or formal first-order reflectivity.

        hkl coordinates are reference-cell r.l.u. ``"unpolarized"`` averages
        independently evaluated s and p reflectivities incoherently.

        :param str polarization: ``"s"``, ``"p"``, or ``"unpolarized"``.
        :param bool first_order: Omit the squared scattered amplitude if true.
        :returns: Scalar or broadcast reflectivity.
        """
        if polarization not in {"s", "p", "unpolarized"}:
            raise ValueError("polarization must be 's', 'p', or 'unpolarized'.")

        def one(pol):
            result = self.evaluate(
                h,
                k,
                l,
                polarization_i=pol,
                polarization_f=pol,
                **kwargs,
            )
            return (
                result.first_order_reflectivity
                if first_order
                else result.reflectivity
            )

        if polarization == "unpolarized":
            value = 0.5 * (one("s") + one("p"))
        else:
            value = one(polarization)
        if np.ndim(value) == 0:
            return float(value)
        return _readonly_array(value, np.float64)

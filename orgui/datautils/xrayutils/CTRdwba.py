"""Distorted-wave Born approximation state and observable amplitudes."""

from collections import OrderedDict
from dataclasses import dataclass
import hashlib

import numpy as np

from ._CTRnative import HAS_CPP_ACCEL, _CTRcalc_cpp
from .CTRoptics import (
    HC_KEV_ANGSTROM,
    LayeredElectricField,
    StratifiedProfile,
    solve_electric_field,
)
from .CTRuc import _coherent_domain_arrays
from .HKLVlieg import UBCalculator, VliegAngles


__all__ = [
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
    array = np.ascontiguousarray(value, dtype=dtype)
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

    @property
    def k0(self):
        """Return the frozen vacuum wavevector in inverse Angstrom."""
        return self.reference.k0


@dataclass(frozen=True)
class DWBAResult:
    """DWBA matrix element and precisely normalized derived observables."""

    prepared: PreparedCTR
    F_contrast: complex | np.ndarray
    unperturbed_amplitude: complex | np.ndarray
    scattered_amplitude: complex | np.ndarray
    total_amplitude: complex | np.ndarray
    bulk_mode: str

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
            _CTRcalc_cpp, "unitcell_F_DWBA_bulk"
        ):
            raise RuntimeError(
                "DWBA requires the native CTR extension with DWBA support."
            )
        if self._crystal.uc_surface_list:
            raise NotImplementedError(
                "DWBA currently supports only the bulk unit cell."
            )
        if not hasattr(self._crystal.uc_bulk, "_E"):
            raise ValueError("Set the bulk unit-cell energy before using DWBA.")

    def _vlieg_angles(self):
        energy_eV = float(self._crystal.uc_bulk._E)
        ub = UBCalculator(self._crystal.uc_bulk, energy_eV * 1e-3)
        ub.setLambda(HC_KEV_ANGSTROM / (energy_eV * 1e-3))
        ub.setU(self._orientation)
        return VliegAngles(ub)

    def _reference(self):
        self._require_available()
        profile = self._crystal._optical_reference_profile(
            _OPTICAL_DELTA_TOLERANCE, _OPTICAL_BETA_TOLERANCE
        )
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
        )

    def _prepared_key(
        self,
        geometry,
        reference,
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
        reference = self._reference()
        key = self._prepared_key(
            geometry,
            reference,
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
        reference = self._reference()
        h_array, k_array, ai, af = arrays
        Qz = reference.k0 * (np.sin(ai) + np.sin(af))
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
        self._validate_prepared(prepared)
        if bulk_mode not in {"unit_cell", "semi_infinite"}:
            raise ValueError("bulk_mode must be 'unit_cell' or 'semi_infinite'.")
        if not np.isfinite(bulk_attenuation) or bulk_attenuation < 0.0:
            raise ValueError("bulk_attenuation must be finite and non-negative.")
        bulk = self._crystal.uc_bulk
        basis, form_factors, _ = bulk.build_selected_basis()
        if bulk._special_formfactors_present:
            raise NotImplementedError(
                "DWBA does not support Python special form-factor callbacks."
            )
        domain_matrix, domain_occupancy = _coherent_domain_arrays(
            bulk.coherentDomainMatrix,
            bulk.coherentDomainOccupancy,
        )
        media = len(prepared.reference.n)

        def field_array(value):
            return np.ascontiguousarray(np.asarray(value).reshape(media, -1))

        medium_index = np.full(len(basis), media - 1, dtype=np.int64)
        amplitude = _CTRcalc_cpp.unitcell_F_DWBA_bulk(
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
            prepared.alpha_i,
            prepared.alpha_f,
            prepared.cos_azimuth,
            prepared.sin_azimuth,
            prepared.k0,
            prepared.polarization_i,
            prepared.polarization_f,
            medium_index,
            bulk_mode == "semi_infinite",
            float(bulk_attenuation),
            np.ascontiguousarray(basis, dtype=np.float64),
            np.ascontiguousarray(form_factors, dtype=np.float64),
            np.ascontiguousarray(bulk.R_mat, dtype=np.float64),
            np.ascontiguousarray(bulk.R_mat_inv, dtype=np.float64),
            domain_matrix,
            domain_occupancy,
        )
        amplitude *= prepared.reference_area / bulk.uc_area
        F_contrast = _readonly_result(
            amplitude, prepared.shape, prepared.scalar, np.complex128
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

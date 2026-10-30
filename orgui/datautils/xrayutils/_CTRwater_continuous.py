"""Explicit numerical erf confirmation and continuous DWBA matrix elements.

The production sharp/tanh paths are in _CTRwater_analytic. No sampled or
finite-segment tanh reference is constructed here.
"""

import copy
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .CTRcontinuous import solve_continuous_field
from .CTRoptics import StratifiedProfile
from ._CTRwater_dwba import _scale
from ._CTRwater import _integration_breaks, _step_density_function


@dataclass(frozen=True)
class _WaterReference:
    water: object
    weight: float
    water_index: int
    boundaries: np.ndarray
    shares: np.ndarray
    values: np.ndarray
    fingerprint: bytes

    def __post_init__(self):
        density = (self.shares[..., 0] + 1j*self.shares[..., 1]) / _scale(
            self.water._E
        )
        solid = np.delete(density, self.water_index, axis=0).sum(axis=0)
        solid.flags.writeable = False
        object.__setattr__(self, "_solid_density", solid)
        object.__setattr__(self, "_bulk_density", self.weight*self.water.bulk_density)
        object.__setattr__(self, "_scale", _scale(self.water._E))
        object.__setattr__(self, "_physical_density", _step_density_function(
            self.water, averaged=True
        ))

    def _water_density(self, z):
        return self.weight * self._physical_density(z)

    def _constant_permittivity(self, left, right):
        if self.water.adsorption and self.water.basis[4]:
            # An adsorption peak is not monotonic: equal endpoints do not
            # bound its interior. Never skip its Maxwell propagation.
            return None
        # Water is monotonic. Its endpoint difference bounds its variation
        # throughout this interval, which contains no solid discontinuity.
        # 1e-15 of water density is below the existing tail and ODE tolerances.
        water = self._water_density(np.array([left, right]))
        if abs(water[1] - water[0]) <= 1e-15 * abs(self._bulk_density):
            return complex(self.epsilon((left + right) / 2))
        return None

    def atomic_density(self, z):
        """Return tagged solid reference densities at physical z in Angstrom."""
        indices = np.searchsorted(self.boundaries, z)
        result = (
            self.shares[:, indices, 0] + 1j * self.shares[:, indices, 1]
        ) / _scale(self.water._E)
        result = np.array(result, copy=True)
        result[self.water_index] = 0.0
        return result

    def epsilon(self, z):
        """Return the full physical step/adsorption reference permittivity."""
        z = np.asarray(z)
        density = self._solid_density[np.searchsorted(self.boundaries, z)]
        density = density + self._water_density(z)
        return (1 - self._scale * density) ** 2

    def reference_density(self, z):
        """Return tagged densities in electrons per cubic Angstrom."""
        density = self.atomic_density(z)
        density[self.water_index] = self._water_density(z)
        return density

    def bounds(self):
        """Return matching planes enclosing all transitions in Angstrom."""
        low, high = self.water.transition_extent()
        return min(low, self.boundaries[0]), max(high, self.boundaries[-1])

    def solve(self, alpha, polarization, *, rtol=2e-10):
        """Solve the immutable reference for local angles in radians."""
        bottom = 1 - complex(*self.values[0, 1:])
        top = 1 - complex(*self.values[-1, 1:])
        return solve_continuous_field(
            self.epsilon,
            self.bounds(),
            bottom,
            top,
            self.water._E,
            alpha,
            polarization,
            z_reference=self.bounds()[1],
            interfaces=np.unique(np.r_[
                self.boundaries, _integration_breaks(self.water, *self.bounds())
            ]) if self.water.adsorption and self.water.basis[4] else self.boundaries,
            constant_permittivity=self._constant_permittivity,
            rtol=rtol,
        )


def _make_continuous(profile, shares, records):
    water_index = next(i for i, r in enumerate(records) if r.role == "water")
    record = records[water_index]
    record.cell._validate_settings()
    if record.cell.optical_reference != "erf" or record.cell.basis[2] == 0:
        from ._CTRwater_analytic import _make_analytic
        return _make_analytic(profile, shares, records)
    water = copy.copy(record.cell)
    water.__dict__ = {
        k: copy.deepcopy(v) for k, v in water.__dict__.items() if k != "_surface_owner"
    }
    for value in water.__dict__.values():
        if isinstance(value, np.ndarray):
            value.flags.writeable = False
    placements = getattr(water, "_placements", None)
    fingerprint = repr(
        (
            water.type,
            water._settings(),
            water.basis.tolist(),
            water.pw,
            water.pos_absolute,
            None if placements is None else [p.tolist() for p in placements],
            record.component_weight,
        )
    ).encode()
    spec = _WaterReference(
        water,
        record.component_weight,
        water_index,
        profile.boundaries.copy(),
        shares.copy(),
        profile.values.copy(),
        fingerprint,
    )
    for value in (spec.boundaries, spec.shares, spec.values):
        value.flags.writeable = False
    values = profile.values.copy()
    n = np.sqrt(spec.epsilon(values[:, 0]))
    # Endpoint rows are semi-infinite, not finite samples of a smooth tail.
    values[1:-1, 1] = 1 - n[1:-1].real
    values[1:-1, 2] = -n[1:-1].imag
    updated = shares.copy()
    updated[water_index] = values[:, 1:] - np.delete(updated, water_index, axis=0).sum(
        axis=0
    )
    return StratifiedProfile(values, profile.boundaries, spec), updated


def _channel_atoms(
    prepared, positions, atom_indices, weights, packed, ii, jf, field_i, field_f
):
    """Evaluate local distorted waves at physical atom positions."""
    if not len(positions):
        return np.empty((0, len(ii)), complex)
    zi, pi, mi = field_i.plane_waves(positions[:, 2])
    zf, pf, mf = (
        (zi, pi, mi) if field_i is field_f else field_f.plane_waves(positions[:, 2])
    )
    zi, pi, mi = zi[:, ii], pi[:, ii], mi[:, ii]
    zf, pf, mf = zf[:, jf], pf[:, jf], mf[:, jf]
    qpar = prepared.Q_parallel
    qpar2 = np.sum(qpar * qpar, axis=1)
    phase = np.exp(
        1j * (positions[:, 0, None] * qpar[:, 0] + positions[:, 1, None] * qpar[:, 1])
    )
    basis = packed["atom_basis"][atom_indices]
    factors = packed["form_factors"][atom_indices]
    result = np.zeros(phase.shape, complex)
    si, sf = np.sin(prepared.alpha_i), np.sin(prepared.alpha_f)
    for sign_i, ai in ((1, pi), (-1, mi)):
        for sign_f, af in ((1, pf), (-1, mf)):
            qz = -(sign_i * zi + sign_f * zf)
            q2 = qpar2 + qz * qz
            form = (factors[:, 10] + factors[:, 11] + 1j * factors[:, 12])[
                :, None
            ] + np.zeros_like(q2)
            for j in range(5):
                form += factors[:, j, None] * np.exp(-factors[:, j + 5, None] * q2)
            form *= (
                np.exp(
                    -(basis[:, 4, None] * qpar2 + basis[:, 5, None] * qz * qz)
                    / (16 * np.pi * np.pi)
                )
                * basis[:, 6, None]
            )
            coeff = ai * af
            pol_i, pol_f = prepared.polarization_i, prepared.polarization_f
            if pol_i == pol_f == "s":
                coeff *= prepared.cos_azimuth
            elif pol_i == "s":
                coeff *= -sign_f * sf * prepared.sin_azimuth
            elif pol_f == "s":
                coeff *= -sign_i * si * prepared.sin_azimuth
            else:
                normal = (
                    prepared.k0**2
                    * prepared.reference.n[0] ** 2
                    * np.cos(prepared.alpha_i)
                    * np.cos(prepared.alpha_f)
                )
                coeff *= (
                    si
                    * sf
                    * (-sign_i * sign_f * prepared.cos_azimuth + normal / (zi * zf))
                )
            result += coeff * form * phase * weights[:, None]
    return result


def _evaluate_continuous(
    prepared, records, packed, bulk_attenuation, *, integral_cache=None
):
    spec = prepared.reference.continuous
    fi, ff = prepared.field_i, prepared.field_f
    ii, jf = prepared.field_i_index, prepared.field_f_index
    count, points = len(records), len(ii)
    atomic, reference = (
        np.zeros((count, points), complex),
        np.zeros((count, points), complex),
    )
    chunk = max(1, min(256, 262144 // max(points, 1)))
    for start in range(0, len(packed["finite_positions"]), chunk):
        selection = slice(start, start + chunk)
        finite = _channel_atoms(
            prepared,
            packed["finite_positions"][selection],
            packed["finite_atom_index"][selection],
            packed["finite_weight"][selection],
            packed,
            ii,
            jf,
            fi,
            ff,
        )
        np.add.at(atomic, packed["finite_record_index"][selection], finite)
    low, high = spec.bounds()
    repeat = packed["bulk_repeat"]
    base = packed["bulk_positions"]
    repeats = max(0, int(np.ceil((np.max(base[:, 2]) - low) / repeat[2])))
    for j in range(repeats + 1):
        positions = base - j * repeat
        weights = packed["bulk_weight"] * np.exp(
            bulk_attenuation * (packed["bulk_repeat_coordinate"] - j)
        )
        amplitude = _channel_atoms(
            prepared,
            positions,
            packed["bulk_atom_index"],
            weights,
            packed,
            ii,
            jf,
            fi,
            ff,
        )
        if j == repeats:
            # Reuse the field solver's radiation branch, including a complex
            # tangential wavevector in absorbing water. A bare principal sqrt
            # can select the opposite branch for a less absorbing substrate.
            kz_i = fi.plane_waves(np.asarray(low - 1.0))[0][ii]
            kz_f = ff.plane_waves(np.asarray(low - 1.0))[0][jf]
            qz = -(kz_i + kz_f)
            phase = prepared.Q_parallel @ repeat[:2] + qz * repeat[2]
            denominator = -np.expm1(-1j * phase - bulk_attenuation)
            if np.any(np.abs(denominator) < 1e-14):
                raise ValueError("Semi-infinite bulk has an exact Bragg pole.")
            amplitude /= denominator
        atomic[0] += amplitude.sum(axis=0)
    selected = np.flatnonzero(prepared.is_specular)
    water = records[spec.water_index].cell
    if not np.isclose(
        water.bulk_density * records[spec.water_index].component_weight,
        spec.water.bulk_density * spec.weight,
        rtol=1e-12,
        atol=1e-14,
    ):
        raise ValueError("Water ambient changed; prepare the CTR again.")
    if not len(selected):
        return atomic, reference
    placements = getattr(water, "_placements", None)
    integral_key = (
        prepared.reference.fingerprint,
        prepared._atomic_model.decomposition_fingerprint,
        prepared.alpha_i.tobytes(),
        prepared.alpha_f.tobytes(),
        prepared.is_specular.tobytes(),
        prepared.reference_area,
        prepared.polarization_i,
        prepared.polarization_f,
        bulk_attenuation,
        repeat[2],
        water.type,
        tuple(water._settings().items()),
        water.basis.tobytes(),
        water.f.tobytes(),
        water._E,
        water.pw,
        water.pos_absolute,
        water.volume,
        water.uc_area,
        records[spec.water_index].component_weight,
        None if placements is None else tuple(v.tobytes() for v in placements),
    )
    if integral_cache is not None and integral_key in integral_cache:
        cached_reference, cached_water = integral_cache[integral_key]
        integral_cache.move_to_end(integral_key)
        atomic[spec.water_index] += cached_water
        return atomic, cached_reference.copy()
    live_low, live_high = water.transition_extent()
    live_weight = records[spec.water_index].component_weight
    if water.type == "step":
        physical_density = _step_density_function(water, averaged=True)

        def live_density(z):
            """Evaluate the same continuum used by the immutable optical reference."""
            return live_weight * physical_density(z)
    else:

        def live_density(z):
            """Evaluate a changed physical profile against the frozen reference."""
            return live_weight * water.zDensity_G(z, 0.0, 0.0)

    low, high = min(low, live_low), max(high, live_high)
    edges = np.unique(
        np.r_[
            low,
            spec.boundaries[(spec.boundaries > low) & (spec.boundaries < high)],
            getattr(spec, "_edges", []),
            _integration_breaks(water, low, high),
            _integration_breaks(spec.water, low, high),
            high,
        ]
    )
    qmax = max(1.0, 2 * prepared.k0 * np.max(np.sin(prepared.alpha_i)))

    def integrate(order):
        """Return reference and localized water amplitudes at one Gauss order."""
        ref_integral = np.zeros_like(reference)
        water_integral = np.zeros(points, complex)
        nodes, weights = leggauss(order)
        grids, quadrature_weights = [], []
        for left, right in zip(edges[:-1], edges[1:]):
            panels = max(1, int(np.ceil((right - left) * qmax / 3)))
            steps = np.linspace(left, right, panels + 1)
            half = np.diff(steps) / 2
            grids.append(
                ((steps[:-1] + steps[1:])[:, None] / 2 + half[:, None] * nodes).ravel()
            )
            quadrature_weights.append((half[:, None] * weights).ravel())
        grid = np.concatenate(grids)
        all_weights = np.concatenate(quadrature_weights)
        chunk = max(1, min(512, 131072 // max(points, 1)))
        for start in range(0, len(grid), chunk):
            z = grid[start : start + chunk]
            w = all_weights[start : start + chunk]
            ei, ni = fi.sample(z)
            ef, nf = (ei, ni) if fi is ff else ff.sample(z)
            tangent_sign = -1 if prepared.polarization_i == "p" else 1
            product = (
                tangent_sign * ei[:, ii[selected]] * ef[:, jf[selected]]
                + ni[:, ii[selected]] * nf[:, jf[selected]]
            )
            if prepared.polarization_i != prepared.polarization_f:
                product *= 0.0
            density = spec.reference_density(z)
            density[0] *= np.exp(bulk_attenuation * z / repeat[2])
            integral = (
                np.einsum("rz,zp->rp", density * w, product) * prepared.reference_area
            )
            ref_integral[:, selected] += integral
            live = live_density(z)
            # Compute water contrast directly, retaining the finite reference as a
            # diagnostic baseline rather than separate divergent upper integrals.
            water_integral[selected] += (
                np.einsum("z,zp->p", (live - density[spec.water_index]) * w, product)
                * prepared.reference_area
            )
        return ref_integral, water_integral

    previous = integrate(20)
    for order in (40, 80):
        current = integrate(order)
        if all(
            np.all(np.abs(a - b) <= 1e-9 + 1e-7 * np.abs(a))
            for a, b in zip(current, previous)
        ):
            break
        previous = current
    else:
        raise RuntimeError("Continuous DWBA quadrature did not converge.")
    reference += current[0]
    atomic[spec.water_index] += current[1]
    # Homogeneous lower half-space, with only transmitted branches.
    # Use the substrate side of the matching plane, which can coincide with
    # a solid discontinuity. The field is continuous but its local kappa is not.
    below = np.nextafter(low, -np.inf)
    zi, pi, _ = fi.plane_waves(np.array(below))
    zf, pf, _ = ff.plane_waves(np.array(below))
    q = -(zi[ii[selected]] + zf[jf[selected]]) - 1j * bulk_attenuation / repeat[2]
    ei, ni = fi.sample(np.array(low))
    ef, nf = ff.sample(np.array(low))
    tangent_sign = -1 if prepared.polarization_i == "p" else 1
    prod = (
        tangent_sign * ei[ii[selected]] * ef[jf[selected]]
        + ni[ii[selected]] * nf[jf[selected]]
    )
    if prepared.polarization_i != prepared.polarization_f:
        prod *= 0.0
    density = spec.reference_density(below)
    density[0] *= np.exp(bulk_attenuation * low / repeat[2])
    density[spec.water_index] = 0.0
    reference[:, selected] += (
        density[:, None] * prod / (1j * q) * prepared.reference_area
    )
    atomic[spec.water_index] += reference[spec.water_index]
    if integral_cache is not None:
        # These depend on the frozen optics and live water, not on atomic
        # coordinates/crystallinity. Keep bounded copies, never mutable outputs.
        integral_cache[integral_key] = (
            reference.copy(),
            atomic[spec.water_index].copy(),
        )
        while len(integral_cache) > 4:
            integral_cache.popitem(last=False)
    return atomic, reference

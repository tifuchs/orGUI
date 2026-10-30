"""Constant-medium DWBA water amplitudes and analytic Epstein interfaces.

Lengths are Angstrom, q is inverse Angstrom, and n=1-delta-i*beta.
Molecular layers are perturbations; they never create optical slabs. The
four plane-wave channels use the water-side medium at each terrace onset,
as for an atomic form factor assigned to its host medium. This is a local
constant-field approximation, not integration through a continuous profile.
"""

from dataclasses import dataclass, replace

import numpy as np
from scipy.special import expit, loggamma

from .CTRoptics import StratifiedProfile, solve_electric_field
from ._CTRwater import _excess, _factors, _parameters
from ._CTRwater_dwba import _scale


def _tanh_ratio(q, sigma):
    x = np.sqrt(3.0) * sigma * np.asarray(q, complex)
    small = abs(x) < 1e-3
    y = np.where(small, 1.0, np.where(x.real < 0, -x, x))
    return np.where(
        small, 1 - x*x/6 + 7*x**4/360 - 31*x**6/15120,
        2*y*np.exp(-y)/(-np.expm1(-2*y)),
    )


def _tanh_minus_one(q, sigma):
    x = np.sqrt(3.0) * sigma * np.asarray(q, complex)
    return np.where(
        abs(x) < 1e-3, -x*x/6 + 7*x**4/360 - 31*x**6/15120,
        _tanh_ratio(q, sigma) - 1,
    )


def _contrast_transform(water, q, reference_z, reference_sigma=0.0):
    """Localized density transform per area, with its finite q=0 limit.

Subtract the reference half-space analytically BEFORE taking q=0. The
Hosemann series uses exp(i*q*d-b**2*q**2/2) and is evaluated once, without
enumerating layers. Its localized analytic continuation is used for complex
DWBA q; the separate half-space transforms need not converge there.
"""
    q = np.asarray(q, complex)
    if np.all(q == 0):
        # Specular reciprocal fields give two identically zero-q channels.
        # Their value is exactly the integrated excess, even in a frozen fit.
        return np.broadcast_to(
            water.bulk_density*(reference_z-water.reference_position), q.shape)
    safe = np.where(q == 0, 1.0, q)
    z, d, s, b, a = _parameters(water)
    # Work relative to the reference plane, avoiding large absolute phases.
    dz = z - reference_z
    terms = _factors(water)
    f0 = sum(c for c, _ in terms)
    df = sum(c * np.expm1(-v*q*q) for c, v in terms)
    f = f0 + df
    gaussian = np.exp(1j*q*dz - 0.5*s*s*q*q)
    if water.type == "layered":
        x = 1j*q*d - 0.5*b*b*q*q
        small = abs(x) < 0.05
        # -d/expm1(x)-i/q, expanded without subtracting two poles.
        xx = np.where(small, 1.0, x)
        regular = np.where(
            small,
            (b*b/(2*d))/(1 + 1j*b*b*q/(2*d))
            + d*(0.5 - x/12 + x**3/720 - x**5/30240 + x**7/1209600),
            -d/np.expm1(xx) - 1j/safe,
        )
        shape = gaussian * (regular + d*(a-1))
        shape += 1j*np.expm1(1j*q*dz - 0.5*s*s*q*q)/safe
        shape = np.where(q == 0, d*(a-0.5)+b*b/(2*d)-dz, shape)
    else:
        gamma, offset = _excess(water, s, a)
        if water.density_profile == "tanh":
            shape = 1j*(np.expm1(1j*q*dz)
                        + np.exp(1j*q*dz)*_tanh_minus_one(q, s))/safe
        else:
            shape = 1j*np.expm1(1j*q*dz - 0.5*s*s*q*q)/safe
        shape = np.where(q == 0, -dz, shape)
        shape += gamma/water.pw * gaussian*np.exp(1j*q*offset)
    # Convert the hard reference subtraction to a tanh reference if requested.
    result = water.pw*(f*shape + 1j*df/safe
                      - 1j*f0*_tanh_minus_one(q, reference_sigma)/safe)
    return result


def _epstein_coefficients(upper_kz, lower_kz, sigma):
    """Return downward-incidence s r,t at the transition center.

Gamma recurrences remove the poles at a critical wavevector and at zero
contrast. sigma is the RMS width of the tanh derivative, in Angstrom.
"""
    width = np.sqrt(12.0)*sigma/np.pi
    a, b = -0.5j*width*lower_kz, -0.5j*width*upper_kz
    same = upper_kz == lower_kz
    denominator = np.where(same, 1.0, a+b)
    r = (b-a)/denominator * np.exp(
        loggamma(1-2*b) - loggamma(1+2*b)
        + 2*loggamma(1+a+b) - 2*loggamma(1+a-b)
    )
    t = 2*b/denominator * np.exp(
        2*loggamma(1+a+b) - loggamma(1+2*a) - loggamma(1+2*b)
    )
    return np.where(same, 0.0, r), np.where(same, 1.0, t)


def _interface_fields(values, boundaries, widths, energy, angles, polarization):
    """Compose physical interface S matrices; never subdivide the profile."""
    field = solve_electric_field(values, energy, angles, polarization,
                                 boundaries=boundaries)
    if not np.any(widths):
        return field
    if polarization != "s":
        raise ValueError(
            "Analytic tanh interfaces support s polarization only; use sharp "
            "for analytic p DWBA or explicitly select erf for numerical validation."
        )
    k = field.kz
    count = len(k)-1
    r = np.empty((count, k.shape[1]), complex)
    t = np.empty_like(r)
    rr = np.empty_like(r)
    tt = np.empty_like(r)
    for j, sigma in enumerate(widths[::-1]):
        if sigma:
            r[j], t[j] = _epstein_coefficients(k[j], k[j+1], sigma)
            rr[j], tt[j] = _epstein_coefficients(k[j+1], k[j], sigma)
        else:
            same = k[j] == k[j+1]
            den = np.where(same, 1.0, k[j]+k[j+1])
            r[j] = np.where(same, 0.0, (k[j]-k[j+1])/den)
            rr[j] = -r[j]
            t[j], tt[j] = 1+r[j], 1-r[j]
    # R[j] is the load in medium j at its upper reference plane.
    load = np.zeros_like(k)
    denominator = np.empty_like(r)
    for j in range(count-1, -1, -1):
        denominator[j] = 1-rr[j]*load[j+1]
        reflected = r[j] + tt[j]*load[j+1]*t[j]/denominator[j]
        distance = field.z_reference[j] - field.z_interfaces[j]
        load[j] = reflected*np.exp(2j*k[j]*distance)
    plus = np.ones_like(k)
    for j in range(count):
        distance = field.z_reference[j] - field.z_interfaces[j]
        plus[j+1] = plus[j]*np.exp(1j*k[j]*distance)*t[j]/denominator[j]
    minus = plus*load
    po = np.divide(plus, k, out=np.zeros_like(k), where=k != 0)
    mo = np.divide(minus, k, out=np.zeros_like(k), where=k != 0)
    if np.any(k[1:-1] == 0):
        # Internal critical plane-wave bases are singular. No ODE fallback:
        # callers can use a neighboring angle or the sharp transfer solver.
        raise ValueError("A tanh stack has an exactly critical internal medium.")
    arrays = (plus, minus, po, mo)
    for value in arrays:
        value.flags.writeable = False
    return replace(field, A_plus=plus, A_minus=minus,
                   _A_plus_over_kz=po, _A_minus_over_kz=mo,
                   r_S=minus[0], t_S=plus[-1])


@dataclass(frozen=True)
class _AnalyticWaterReference:
    """Frozen onset-only optics, separate from live molecular parameters."""

    water_index: int
    centers: np.ndarray
    probabilities: np.ndarray
    heights: np.ndarray
    sigma: float
    energy: float
    values: np.ndarray
    boundaries: np.ndarray
    widths: np.ndarray
    shares: np.ndarray
    fingerprint: bytes

    def solve(self, angles, polarization):
        """Return constant-medium fields with analytic interface amplitudes."""
        return _interface_fields(self.values, self.boundaries, self.widths,
                                 self.energy, angles, polarization)

    def reference_density(self, z):
        """Return the declared onset-only reference, in electrons/Angstrom^3."""
        z = np.asarray(z)
        indices = np.searchsorted(self.boundaries, z)
        result = ((self.shares[:, indices, 0]+1j*self.shares[:, indices, 1])
                  / _scale(self.energy)).copy()
        fractions = ((z[..., None] > self.centers) if not self.sigma else
                     expit(np.pi*(z[..., None]-self.centers)
                           / (np.sqrt(3)*self.sigma)))
        bulk = complex(*self.shares[self.water_index, -1])/_scale(self.energy)
        result[self.water_index] = bulk*np.sum(fractions*self.probabilities, axis=-1)
        return result

    def atomic_density(self, z):
        """Return tagged solid densities for plotting, excluding water."""
        result = self.reference_density(z)
        result[self.water_index] = 0.0
        return result

    def bounds(self):
        """Return a plotting window; these bounds never enter the solver."""
        radius = max(1.0, 18*np.sqrt(12)*self.sigma/np.pi)
        return (min(self.boundaries[0], self.centers.min()-radius),
                max(self.boundaries[-1], self.centers.max()+radius))

    def epsilon(self, z):
        """Return the reference index squared for density diagnostics."""
        return (1-_scale(self.energy)*self.reference_density(z).sum(axis=0))**2


def _make_analytic(profile, shares, records):
    index = next(i for i, r in enumerate(records) if r.role == "water")
    water = records[index].cell
    heights, probabilities = getattr(water, "_placements", None) or ([0.0], [1.0])
    heights, probabilities = np.array(heights), np.array(probabilities)
    centers = water.reference_position + heights
    sigma = _parameters(water)[2] if water.optical_reference == "tanh" else 0.0
    widths = np.zeros(len(profile.boundaries))
    widths[np.searchsorted(profile.boundaries, centers)] = sigma
    arrays = [centers, probabilities, heights, profile.values.copy(),
              profile.boundaries.copy(), widths, shares.copy()]
    for value in arrays:
        value.flags.writeable = False
    fingerprint = b"analytic-water-v1" + b"".join(x.tobytes() for x in arrays)
    spec = _AnalyticWaterReference(index, *arrays[:3], sigma, water._E,
                                   *arrays[3:], fingerprint)
    return StratifiedProfile(profile.values, profile.boundaries, spec), shares


def _water_contrast(prepared, record, record_index, **unused):
    """Evaluate one closed water transform per terrace and DWBA channel."""
    water = record.cell
    shares = prepared._atomic_model.reference_shares[record_index]
    if not np.isclose(water.bulk_density*record.component_weight,
                      complex(*shares[0])/_scale(water._E), rtol=1e-12, atol=1e-14):
        raise ValueError("Water ambient changed; prepare the CTR again.")
    result = np.zeros(len(prepared.alpha_i), complex)
    selected = np.flatnonzero(prepared.is_specular)
    if not len(selected) or prepared.polarization_i != prepared.polarization_f:
        return result
    spec = prepared.reference.continuous
    if not isinstance(spec, _AnalyticWaterReference):
        raise ValueError(
            "Analytic water DWBA requires a prepared water onset reference."
        )
    live = getattr(water, "_placements", None) or ([0.0], [1.0])
    if not (np.array_equal(live[0], spec.heights)
            and np.array_equal(live[1], spec.probabilities)):
        raise ValueError("Water terraces changed; prepare the CTR again.")
    fi, ff = prepared.field_i, prepared.field_f
    ii, jf = prepared.field_i_index[selected], prepared.field_f_index[selected]
    media = len(prepared.reference.n)
    # Like an atom's Gaussian form factor, the entire localized perturbation
    # uses the constant wavevector on the water side of its reference onset.
    host = media-1-np.searchsorted(spec.boundaries, spec.centers, side="right")
    ki, kf = fi.kz[host[:, None], ii], ff.kz[host[:, None], jf]
    for sign_i, ai in ((1, fi.A_plus), (-1, fi.A_minus)):
        for sign_f, af in ((1, ff.A_plus), (-1, ff.A_minus)):
            coeff = ai[host[:, None], ii]*af[host[:, None], jf]
            if not np.any(coeff):
                continue
            if prepared.polarization_i == "p":
                normal = (prepared.k0**2*prepared.reference.n[0]**2
                          * np.cos(prepared.alpha_i[selected])
                          * np.cos(prepared.alpha_f[selected]))
                coeff *= (np.sin(prepared.alpha_i[selected])
                          * np.sin(prepared.alpha_f[selected])
                          * (-sign_i*sign_f + normal/(ki*kf)))
            q = -(sign_i*ki + sign_f*kf)
            phase = (sign_i*ki*(fi.z_reference[host]-spec.centers)[:, None]
                     + sign_f*kf*(ff.z_reference[host]-spec.centers)[:, None])
            transform = _contrast_transform(
                water, q, (spec.centers-spec.heights)[:, None], spec.sigma)
            result[selected] += np.sum(
                spec.probabilities[:, None]*coeff*np.exp(1j*phase)*transform, axis=0)
    return result*prepared.reference_area*record.component_weight

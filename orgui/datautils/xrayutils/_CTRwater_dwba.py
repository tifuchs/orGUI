"""Water reference assembly and numerical integrals for independent oracles.

Production sharp/tanh water matrix elements live in _CTRwater_analytic.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss

from .CTRoptics import HC_KEV_ANGSTROM, StratifiedProfile, profile_boundaries
from ._CTRwater import _integration_breaks


def _scale(energy):
    wavelength = HC_KEV_ANGSTROM / (energy * 1e-3)
    return 2.8179403262e-5 * wavelength * wavelength / (2 * np.pi)


def _with_water_reference(values, shares, records):
    old_boundaries = (values.boundaries if isinstance(values, StratifiedProfile)
                      else profile_boundaries(values))
    water_records = [(i, r) for i, r in enumerate(records) if r.role == "water"]
    boundaries = list(old_boundaries)
    for _, record in water_records:
        water = record.cell
        center = water.reference_position  # also validates a homogeneous ambient
        heights, _ = getattr(water, "_placements", None) or ([0.0], [1.0])
        boundaries.extend(center + np.asarray(heights))
    boundaries = np.unique(boundaries)
    centers = np.r_[
        boundaries[0] - 1, 0.5 * (boundaries[:-1] + boundaries[1:]), boundaries[-1] + 1
    ]
    indices = np.searchsorted(old_boundaries, centers)
    tagged = shares[:, indices].copy()
    for i, record in water_records:
        water = record.cell
        heights, probabilities = getattr(water, "_placements", None) or ([0.0], [1.0])
        fraction = sum(
            p * (centers > water.reference_position + h)
            for h, p in zip(heights, probabilities)
        )
        density = water.bulk_density * record.component_weight * _scale(water._E)
        tagged[i, :, 0] = fraction * density.real
        tagged[i, :, 1] = fraction * density.imag
    return StratifiedProfile(
        np.column_stack((centers, tagged.sum(axis=0))), boundaries
    ), tagged


def _water_contrast(
    prepared, record, record_index, *, rtol=1e-7, atol=1e-9, cache=None
):
    """Converge the localized integral independently of the plotting grid."""
    water = record.cell
    shares = prepared._atomic_model.reference_shares[record_index]
    if not np.isclose(
        water.bulk_density * record.component_weight,
        complex(*shares[0]) / _scale(water._E), rtol=1e-12, atol=1e-14,
    ):
        raise ValueError("Water ambient changed; prepare the CTR again.")
    if not np.any(prepared.is_specular):
        return np.zeros(len(prepared.alpha_i), complex)
    if cache is not None:
        placements = getattr(water, "_placements", None)
        key = (
            prepared.reference.fingerprint, shares.tobytes(),
            prepared.alpha_i.tobytes(), prepared.alpha_f.tobytes(),
            prepared.is_specular.tobytes(), prepared.reference_area,
            prepared.polarization_i, prepared.polarization_f,
            water.type, tuple(water._settings().items()), water.basis.tobytes(),
            water.f.tobytes(), water._E, water.pw, water.pos_absolute,
            water.volume, water.uc_area, record.component_weight, rtol, atol,
            None if placements is None else tuple(
                np.asarray(v).tobytes() for v in placements
            ),
        )
        if key in cache:
            cache.move_to_end(key)
            return cache[key].copy()
    previous = _water_contrast_order(prepared, record, record_index, 16)
    for order in (32, 64):
        current = _water_contrast_order(prepared, record, record_index, order)
        if np.all(np.abs(current - previous) <= atol + rtol * np.abs(current)):
            if cache is not None:
                cache[key] = current.copy()
                while len(cache) > 4:
                    cache.popitem(last=False)
            return current
        previous = current
    raise RuntimeError("Water DWBA quadrature did not converge; check profile widths.")


def _water_contrast_order(prepared, record, record_index, order):
    """Integrate density minus frozen reference before any half-space limits."""
    water = record.cell
    points = len(prepared.alpha_i)
    result = np.zeros(points, complex)
    selected = np.flatnonzero(prepared.is_specular)
    rho_bulk = water.bulk_density * record.component_weight
    scale = _scale(water._E)
    refshares = prepared._atomic_model.reference_shares[record_index]
    frozen_bulk = complex(*refshares[0]) / scale
    if not np.isclose(rho_bulk, frozen_bulk, rtol=1e-12, atol=1e-14):
        raise ValueError("Water ambient changed; prepare the CTR again.")
    if not len(selected):
        return result
    low, high = water.transition_extent()
    boundaries = prepared.reference.profile.boundaries
    low = min(low, boundaries[0])
    high = max(high, boundaries[-1])
    breaks = np.unique(
        np.r_[
            low,
            boundaries[(boundaries > low) & (boundaries < high)],
            _integration_breaks(water, low, high),
            high,
        ]
    )
    fi, ff = prepared.field_i, prepared.field_f
    ii, jf = prepared.field_i_index[selected], prepared.field_f_index[selected]
    media = len(prepared.reference.n)
    ai, af = prepared.alpha_i[selected], prepared.alpha_f[selected]
    si, sf = np.sin(ai), np.sin(af)
    ki = np.asarray(fi.kz).reshape(media, -1)
    kf = np.asarray(ff.kz).reshape(media, -1)
    ip = np.asarray(fi.A_plus).reshape(media, -1)
    im = np.asarray(fi.A_minus).reshape(media, -1)
    fp = np.asarray(ff.A_plus).reshape(media, -1)
    fm = np.asarray(ff.A_minus).reshape(media, -1)
    xg, wg = leggauss(order)
    grids = {}
    for left, right in zip(breaks[:-1], breaks[1:]):
        medium = media - 1 - np.searchsorted(boundaries, 0.5 * (left + right))
        # Each panel resolves at least four radians of the fastest channel.
        qmax = np.max(np.abs(ki[medium, ii])) + np.max(np.abs(kf[medium, jf]))
        panels = max(1, int(np.ceil((right - left) * max(qmax, 1.0) / 4)))
        edges = np.linspace(left, right, panels + 1)
        half = 0.5 * np.diff(edges)
        z = (0.5 * (edges[:-1] + edges[1:])[:, None] + half[:, None] * xg).ravel()
        weights = (half[:, None] * wg).ravel()
        grids.setdefault(medium, []).append((z, weights))
    for medium, pieces in grids.items():
        # Preserve every quadrature panel/node, but sample adjacent panels in
        # bounded batches instead of rebuilding water factors for each panel.
        z = np.concatenate([piece[0] for piece in pieces])
        weights = np.concatenate([piece[1] for piece in pieces])
        density = np.empty(len(z), complex)
        for start in range(0, len(z), 512):
            sl = slice(start, start + 512)
            density[sl] = record.component_weight * water.zDensity_G(z[sl], 0.0, 0.0)
        density -= complex(*refshares[medium]) / scale
        weighted = weights * density
        if not np.any(weighted):
            continue
        for sign_i, amp_i in ((1, ip), (-1, im)):
            for sign_f, amp_f in ((1, fp), (-1, fm)):
                coeff = amp_i[medium, ii] * amp_f[medium, jf]
                if prepared.polarization_i != prepared.polarization_f:
                    continue
                if prepared.polarization_i == "p":
                    normal = (
                        prepared.k0**2
                        * prepared.reference.n[0] ** 2
                        * np.cos(ai)
                        * np.cos(af)
                    )
                    coeff *= (
                        si
                        * sf
                        * (
                            -sign_i * sign_f
                            + normal / (ki[medium, ii] * kf[medium, jf])
                        )
                    )
                q = -(sign_i * ki[medium, ii] + sign_f * kf[medium, jf])
                phase = (
                    sign_i * ki[medium, ii] * fi.z_reference[medium]
                    + sign_f * kf[medium, jf] * ff.z_reference[medium]
                )
                for start in range(0, len(selected), 128):
                    sl = slice(start, start + 128)
                    integrals = np.zeros(len(q[sl]), complex)
                    for first in range(0, len(z), 512):
                        block = slice(first, first + 512)
                        integrals += np.einsum(
                            "z,zp->p", weighted[block],
                            np.exp(1j * (z[block, None] * q[sl] + phase[sl])),
                            optimize=False,
                        )
                    result[selected[sl]] += (
                        coeff[sl] * integrals * prepared.reference_area
                    )
    return result

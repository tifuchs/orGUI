"""Continuous water profiles in Angstrom and electrons/Angstrom cubed.

The plain profiles use independent O + 2H atoms, not measured molecular
scattering factors. Legacy oxygen conventions are retained on import.
"""

import numpy as np
from scipy.special import erf, expit


def _parameters(water):
    water._validate_settings()
    values = np.asarray(water.basis, dtype=float)
    if values.shape != (5,) or not np.all(np.isfinite(values)):
        raise ValueError("Water parameters must be five finite scalars.")
    z, d, s, b, a = values
    if d <= 0 or s < 0 or b < 0 or a < 0 or not np.isfinite(water.pw) or water.pw <= 0:
        raise ValueError(
            "Water spacing/density must be positive; widths/A0 nonnegative."
        )
    height = water.volume / water.uc_area
    return water.pos_absolute + z * height, d * height, s * height, b * height, a


def _factors(water):
    if not hasattr(water, "f"):
        raise ValueError("Set the WaterModel energy before evaluation.")
    if water.form_factor == "forward":
        # A continuum susceptibility profile; no additional atomic convolution.
        return [(sum(c for c, _ in _molecular_factors(water)), 0.0)]
    if water.form_factor in {"oxygen", "scaled_oxygen", "oxide"}:
        index = 2 if water.form_factor == "oxide" else 0
        scale = 10 / 8 if water.form_factor == "scaled_oxygen" else 1.0
        row = water.f[index]
        return [(scale * (row[10] + row[11] + 1j * row[12]), 0.0)] + [
            (scale * row[j], row[j + 5]) for j in range(5)
        ]
    return _molecular_factors(water)


def _molecular_factors(water):
    terms = []
    for row, scale in zip(water.f[:2], (1.0, 2.0)):
        # Forward normalization removes the small Waasmaier fit residual at
        # q=0, keeping molecular density consistent with 10 electrons/H2O.
        electrons = 8.0 if scale == 1.0 else 1.0
        constant = electrons - np.sum(row[:5]) + row[11] + 1j * row[12]
        terms.append((scale * constant, 0.0))
        terms.extend((scale * row[j], row[j + 5]) for j in range(5))
    return terms


def _bulk_density(water):
    return water.pw * sum(c for c, _ in _factors(water))


def _excess(water, s, a):
    if water.type != "step" or not water.adsorption or a == 0:
        return 0.0, 0.0
    # Legacy A0 is the number of excess scattering units per lateral cell.
    gamma = a / water.uc_area
    if water.adsorption_position == "onset" or s == 0:
        return gamma, 0.0
    relative = gamma / (np.sqrt(2 * np.pi) * s * water.pw)
    p0, p1, p2 = 0.85886939, 1.03841354, -1.35438352
    offset = s * (
        np.sqrt(2 * np.log(2))
        + 0.5 * p0 * (1 - erf((np.log10(relative) - p2) / (np.sqrt(2) * p1)))
    )
    return gamma, offset


def _reference_position(water):
    z, d, s, b, a = _parameters(water)
    if water.type.startswith("layered"):
        if b == 0:
            raise ValueError(
                "Layered DWBA water requires sigma_bar > 0 for a homogeneous ambient."
            )
        return z + d * (0.5 - a) - b * b / (2 * d)
    gamma, _ = _excess(water, s, a)
    return z - gamma / water.pw


def _fourier(water, q):
    """Amplitude per lateral cell; the half-space transform excludes q=0."""
    q = np.asarray(q, dtype=complex)
    if np.any(q == 0):
        raise ValueError("The semi-infinite water amplitude is undefined at q_z=0.")
    z, d, s, b, a = _parameters(water)
    f = sum(c * np.exp(-v * q * q) for c, v in _factors(water))
    if water.type.startswith("layered"):
        x = 1j * q * d - 0.5 * b * b * q * q
        shape = d * np.exp(-0.5 * s * s * q * q) * (a - 1 - 1 / np.expm1(x))
    elif water.density_profile == "tanh":
        x = np.sqrt(3.0) * s * q
        # x/sinh(x), stable in the sharp limit and at large positive real x.
        small = np.abs(x) < 1e-4
        safe = np.where(small, 1.0, x)
        positive = np.where(safe.real < 0, -safe, safe)
        ratio = 2 * positive * np.exp(-positive) / (-np.expm1(-2 * positive))
        ratio = np.where(small, 1 - x * x / 6 + 7 * x**4 / 360, ratio)
        shape = 1j * ratio / q
        gamma, offset = _excess(water, s, a)
        shape += gamma / water.pw * np.exp(-0.5*s*s*q*q + 1j*q*offset)
    else:
        gamma, offset = _excess(water, s, a)
        shape = np.exp(-0.5 * s * s * q * q) * (
            1j / q + gamma / water.pw * np.exp(1j * q * offset)
        )
    return water.uc_area * water.pw * f * shape * np.exp(1j * q * z)


def _gaussian(x, variance):
    if variance <= 0:
        raise ValueError(
            "A zero-width delta layer cannot be sampled; use a positive sigma_0."
        )
    return np.exp(-0.5 * x * x / variance) / np.sqrt(2 * np.pi * variance)


def _step(x, variance):
    if variance == 0:
        return np.where(x < 0, 0.0, np.where(x > 0, 1.0, 0.5))
    return 0.5 * (1 + erf(x / np.sqrt(2 * variance)))


def _step_density_function(water, *, averaged=False):
    """Freeze factors and lengths for repeated physical step-density sampling."""
    z0, _, sigma, _, a = _parameters(water)
    terms = _factors(water)
    gamma, offset = _excess(water, sigma, a)
    placements = getattr(water, "_placements", None) if averaged else None
    heights, probabilities = placements or ([0.0], [1.0])
    centers = z0 + np.asarray(heights)
    probabilities = np.asarray(probabilities)
    pw = water.pw
    tanh = water.density_profile == "tanh"
    width = np.sqrt(12) * sigma / np.pi
    coefficients = np.asarray([c for c, _ in terms])
    variances = sigma*sigma + 2*np.asarray([v for _, v in terms])
    positive_widths = np.all(variances > 0)
    if positive_widths:
        inverse_widths = 1/np.sqrt(2*variances)
        gaussian_weights = gamma*coefficients/np.sqrt(2*np.pi*variances)

    def density(z):
        """Evaluate physical density, optionally averaged over frozen terraces."""
        x = np.asarray(z)[..., None] - centers
        if np.ndim(z) == 0 and positive_widths:
            # ODE references query one z at a time. Evaluate all molecular
            # factor terms/terraces together instead of a Python loop per term.
            xx = x[:, None]
            shape = (
                expit(2*xx/width) if tanh and width
                else 0.5*(1+erf(xx*inverse_widths))
            )
            values = pw*coefficients*shape
            if gamma:
                values += gaussian_weights*np.exp(-0.5*(xx-offset)**2/variances)
            return np.sum(values*probabilities[:, None])
        result = np.zeros(x.shape, complex)
        for coefficient, atomic_width in terms:
            variance = sigma*sigma + 2*atomic_width
            shape = (
                expit(2*x/width) if tanh and width else _step(x, variance)
            )
            result += pw*coefficient*shape
            if gamma:
                result += gamma*coefficient*_gaussian(x-offset, variance)
        return np.sum(result*probabilities, axis=-1)

    return density


def _density(water, z):
    z = np.asarray(z, dtype=float)
    if not np.all(np.isfinite(z)):
        raise ValueError("Water density coordinates must be finite.")
    z0, d, s, b, a = _parameters(water)
    rho = np.zeros(z.shape, complex)
    if z.size == 0:
        return rho
    if water.type == "step":
        return _step_density_function(water)(z)
    terms = _factors(water)
    if water.type.startswith("layered"):
        # Include all Gaussian centers whose eight-sigma tail reaches the grid.
        # Solving nd - 8*sqrt(s^2+nb^2+2v) > zmax avoids a fixed extra-layer count.
        xmax = max(0.0, float(np.max(z)) - z0)
        v0 = s * s + 2 * max(v for _, v in terms)
        root = (8 * b + np.sqrt(64 * b * b + 4 * d * (xmax + 8 * np.sqrt(v0)))) / (
            2 * d
        )
        count = int(np.ceil(root * root)) + 2
        if count > 200000:
            raise ValueError(
                "Water layering requires too many terms; "
                "increase spacing or reduce broadening."
            )
        indices = np.arange(0 if a else 1, count)
        coefficients = np.array([c for c, _ in terms], complex)
        variances = (
            s * s + indices[:, None] * b * b + 2 * np.array([v for _, v in terms])
        )
        if np.any(variances <= 0):
            raise ValueError(
                "A zero-width delta layer cannot be sampled; use a positive sigma_0."
            )
        centers = np.broadcast_to(z0 + indices[:, None] * d, variances.shape).ravel()
        weights = (
            np.where(indices == 0, a, 1.0)[:, None]
            * water.pw
            * d
            * coefficients
            / np.sqrt(2 * np.pi * variances)
        ).ravel()
        variances = variances.ravel()
        flat, output = z.ravel(), rho.ravel()
        # Bound temporary Gaussian matrices and stream long layering tails.
        for start in range(0, len(centers), 256):
            block = slice(start, start + 256)
            for first in range(0, flat.size, 2048):
                selected = slice(first, first + 2048)
                exponent = -0.5 * (flat[selected][None, :] - centers[block, None]) ** 2
                exponent /= variances[block, None]
                output[selected] += np.sum(
                    weights[block, None] * np.exp(exponent), axis=0
                )
    return rho


def _extent(water, tolerance=1e-10):
    """Finite transition interval, with homogeneous tails to stated tolerance."""
    if not np.isfinite(tolerance) or not 0 < tolerance < 1:
        raise ValueError("Water tail tolerance must lie strictly between zero and one.")
    z, d, s, b, a = _parameters(water)
    width = max(s, np.sqrt(2 * max(v for _, v in _factors(water))))
    low = z - 12 * max(width, 0.01)
    high = z + 12 * max(width, 0.01)
    if water.density_profile == "tanh" or water.optical_reference == "tanh":
        radius = max(
            12 * max(width, 0.01), -np.log(tolerance) * np.sqrt(12) * s / (2 * np.pi)
        )
        low, high = z - radius, z + radius
    if water.type.startswith("layered"):
        _reference_position(water)  # reject non-decaying layering
        # Dominant reciprocal-lattice harmonic decays as exp(-2*pi^2*b^2*z/d^3).
        high += max(12 * b * b / d, -np.log(tolerance) * d**3 / (2 * np.pi**2 * b * b))
    gamma, offset = _excess(water, s, a)
    if gamma:
        high += abs(offset)
    zg = _reference_position(water)
    return min(low, zg - 0.01), max(high, zg + 0.01)


def _integration_breaks(water, low, high):
    """Resolve intrinsic features even when their width is below a wave period."""
    z, d, s, b, a = _parameters(water)
    heights, _ = getattr(water, "_placements", None) or ([0.0], [1.0])
    variances = np.unique(s * s + 2 * np.array([v for _, v in _factors(water)]))
    points = []
    offsets = np.array([-8, -4, -1, 0, 1, 4, 8])
    for height in heights:
        center = z + height
        count = 1
        if water.type.startswith("layered"):
            # Broad later layers are resolved by the <=4 A wave panels.
            # Seed narrow peaks so two orders cannot both miss a peak.
            narrow_count = int(
                np.ceil(max(0.0, min(0.5, d / 4) ** 2 - variances[0]) / (b * b))
            )
            count = max(1, min(narrow_count, int(np.ceil((high - center) / d)) + 1))
        for j in range(count):
            widths = np.sqrt(variances + j * b * b)
            points.extend((center + j * d + widths[:, None] * offsets).ravel())
        gamma, offset = _excess(water, s, a)
        if gamma:
            points.extend(
                (center + offset + np.sqrt(variances)[:, None] * offsets).ravel()
            )
    points = np.unique(points)
    return points[(points > low) & (points < high)]

"""Continuous planar optical fields, including the exact Epstein transition.

Indices follow n=1-delta-i*beta; z increases from substrate to ambient.
Incident electric amplitude is one at the specified finite reference plane.
"""

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import expit, loggamma

from .CTRoptics import HC_KEV_ANGSTROM


def _kappa(epsilon, parallel, k0):
    value = k0 * np.sqrt(np.asarray(epsilon, complex) - parallel**2)
    return np.where(
        (value.imag > 0) | ((value.imag == 0) & (value.real < 0)), -value, value
    )


@dataclass(frozen=True)
class ContinuousElectricField:
    """Immutable continuous field with a local-ambient normalization.

    ``sample(z)`` returns tangential and normal electric components, with
    z in Angstrom and the angle axis last. ``r_S`` uses the Renaud p sign.
    """

    r_S: np.ndarray
    t_S: np.ndarray
    alpha: np.ndarray
    n_ambient: complex
    n_substrate: complex
    k0: float
    z_reference: float
    polarization: str
    solver: str
    _sample: object = field(repr=False)
    _waves: object = field(repr=False, default=None)

    def __post_init__(self):
        for name in ("r_S", "t_S", "alpha"):
            value = np.array(getattr(self, name), copy=True)
            value.flags.writeable = False
            object.__setattr__(self, name, value)

    def sample(self, z):
        """Return (tangential, normal) electric fields at z in Angstrom."""
        z = np.asarray(z, float)
        if not np.all(np.isfinite(z)):
            raise ValueError("Field sample coordinates must be finite.")
        return self._sample(z)

    def plane_waves(self, z):
        """Return local (kz, A_plus, A_minus) at z in Angstrom.

        This is a local decomposition of the field and its derivative, not a
        claim that a continuous profile has constant wavevectors between z's.
        """
        return self._waves(np.asarray(z, float))


def _hyp2f1_value_derivative(a, b, c, x):
    """Return 2F1 and x*d(2F1)/dx for 0 <= x <= 1/2.

    Connection formulas keep the argument in this rapidly convergent disk.
    The derivative is accumulated from the same series, without differencing.
    Ill-conditioned entries use arbitrary precision, never a silent truncation.
    """
    shape = np.broadcast_shapes(np.shape(a), np.shape(b), np.shape(c), np.shape(x))
    term = np.ones(shape, complex)
    value = term.copy()
    derivative = np.zeros(shape, complex)
    absolute_sum = np.ones(shape)
    absolute_derivative = np.zeros(shape)
    for n in range(1, 513):
        term *= ((a + n - 1) * (b + n - 1) / ((c + n - 1) * n)) * x
        value += term
        derivative += n * term
        size = np.abs(term)
        absolute_sum += size
        absolute_derivative += n * size
        converged = (size <= 2e-15 * np.maximum(np.abs(value), 1e-300)) & (
            n * size <= 2e-15 * np.maximum(np.abs(derivative), 1e-300)
        )
        if np.all(converged):
            break
    bad = (
        ~converged
        | ~np.isfinite(value)
        | ~np.isfinite(derivative)
        | (absolute_sum * 2e-16 > 1e-11 * np.maximum(np.abs(value), 1e-300))
        | (absolute_derivative * 2e-16 > 1e-11 * np.maximum(np.abs(derivative), 1e-300))
    )
    if np.any(bad):
        import mpmath as mp

        arrays = [np.broadcast_to(v, shape) for v in (a, b, c, x)]
        with mp.workdps(50):
            for index in zip(*np.nonzero(bad)):
                aa, bb, cc, xx = (complex(v[index]) for v in arrays)
                value[index] = complex(mp.hyp2f1(aa, bb, cc, xx))
                derivative[index] = complex(
                    xx * aa * bb / cc * mp.hyp2f1(aa + 1, bb + 1, cc + 1, xx)
                )
    return value, derivative


def solve_continuous_field(
    epsilon,
    bounds,
    n_substrate,
    n_ambient,
    energy_eV,
    alpha,
    polarization="s",
    *,
    z_reference=0.0,
    rtol=2e-10,
    atol=2e-12,
    interfaces=(),
    constant_permittivity=None,
    _segment_propagator=None,
):
    """Solve Maxwell's planar s or p field for a supplied dielectric profile.

    :param callable epsilon: Relative permittivity at physical z (Angstrom).
    :param tuple bounds: Finite (lower, upper) matching planes in Angstrom.
    :param complex n_substrate: Homogeneous lower-half-space index.
    :param complex n_ambient: Homogeneous upper-half-space index.
    :param float energy_eV: Photon energy in eV.
    :param alpha: Glancing angles in radians in the incident medium.
    :param str polarization: ``s`` or ``p``; p solves the magnetic-field equation.
    :param float z_reference: Incident plane-wave phase origin in Angstrom.
    :param interfaces: Known discontinuities at which to split integration.
    :param callable constant_permittivity: Optional callback ``(left, right)``
        returning a constant epsilon for a known homogeneous interval, or None.
        Homogeneous intervals propagate analytically, including p admittance.
    :returns: ContinuousElectricField with dense field evaluation.
    """
    angles = np.atleast_1d(np.asarray(alpha, float))
    lower, upper = map(float, bounds)
    if (
        angles.ndim != 1
        or not np.all(np.isfinite(angles))
        or np.any((angles <= 0) | (angles > np.pi / 2))
    ):
        raise ValueError("Angles must be a finite 1D array in (0, pi/2].")
    if (
        not np.isfinite(energy_eV)
        or energy_eV <= 0
        or not lower < upper
        or not np.all(np.isfinite([lower, upper, z_reference]))
    ):
        raise ValueError(
            "Positive energy and increasing finite matching planes required."
        )
    if polarization not in {"s", "p"}:
        raise ValueError("polarization must be s or p.")
    if (
        not np.all(np.isfinite([n_substrate, n_ambient, rtol, atol]))
        or rtol <= 0
        or atol <= 0
        or n_substrate == 0
        or n_ambient == 0
    ):
        raise ValueError("Finite nonzero indices and positive tolerances required.")
    k0 = 2 * np.pi / (HC_KEV_ANGSTROM / (energy_eV * 1e-3))
    parallel = complex(n_ambient) * np.cos(angles)
    kb = _kappa(n_substrate**2, parallel, k0)
    kt = _kappa(n_ambient**2, parallel, k0)
    count = len(angles)
    p = polarization == "p"
    initial = np.r_[np.ones(count, complex), 1j * kb / (n_substrate**2 if p else 1)]

    def rhs(z, y):
        """Return the first-order Maxwell system at physical z."""
        eps = complex(epsilon(z))
        first, second = y[:count], y[count:]
        if p:
            return np.r_[
                eps * second, k0 * k0 * (parallel * parallel / eps - 1) * first
            ]
        return np.r_[second, k0 * k0 * (parallel * parallel - eps) * first]

    edges = np.unique(
        np.r_[
            lower,
            np.asarray(interfaces)[
                (np.asarray(interfaces) > lower) & (np.asarray(interfaces) < upper)
            ],
            upper,
        ]
    )
    solutions = []
    state = initial
    for left, right in zip(edges[:-1], edges[1:]):
        constant = (
            constant_permittivity(left, right)
            if constant_permittivity is not None
            else None
        )
        if constant is not None:

            def propagate(z, origin=left, y=state.copy(), eps=complex(constant)):
                """Propagate one homogeneous Maxwell segment exactly."""
                z = np.asarray(z)
                distance = (z - origin).ravel()[None, :]
                kap2 = k0 * k0 * (eps - parallel * parallel)
                phase = np.sqrt(kap2)[:, None] * distance
                cosine = np.cos(phase)
                sine_over_k = distance * np.sinc(phase / np.pi)
                admittance = eps if p else 1.0
                first = (
                    y[:count, None] * cosine
                    + admittance * y[count:, None] * sine_over_k
                )
                second = (
                    -kap2[:, None] * y[:count, None] * sine_over_k / admittance
                    + y[count:, None] * cosine
                )
                return np.concatenate((first, second)).reshape((2 * count,) + z.shape)

            solutions.append(propagate)
            state = propagate(right)
            continue
        if _segment_propagator is not None:
            propagate = _segment_propagator(left, right, state.copy(), parallel, k0)
            solutions.append(propagate)
            state = propagate(right)
            continue
        solution = solve_ivp(
            rhs,
            (left, right),
            state,
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=True,
        )
        if not solution.success:
            raise RuntimeError(
                f"Continuous optical integration failed: {solution.message}"
            )
        solutions.append(solution.sol)
        state = solution.y[:, -1]
    incident = 0.5 * (
        state[:count] + state[count:] * (n_ambient**2 if p else 1) / (1j * kt)
    )
    reflected = 0.5 * (
        state[:count] - state[count:] * (n_ambient**2 if p else 1) / (1j * kt)
    )
    normalization = np.exp(1j * kt * (upper - z_reference)) / incident
    r = reflected / incident * np.exp(2j * kt * (upper - z_reference))
    # In the Renaud p basis magnetic amplitudes have the same reflection sign.
    t = (
        normalization
        * np.exp(-1j * kb * (lower - z_reference))
        * (n_ambient / n_substrate if p else 1)
    )

    def states(z):
        """Evaluate matched field states with the angle axis first."""
        flat = z.ravel()
        y = np.empty((2 * count, flat.size), complex)
        below, above = flat < lower, flat > upper
        y[:, below] = initial[:, None] * np.tile(
            np.exp(1j * kb[:, None] * (flat[below] - lower)), (2, 1)
        )
        inc = incident[:, None] * np.exp(1j * kt[:, None] * (flat[above] - upper))
        ref = reflected[:, None] * np.exp(-1j * kt[:, None] * (flat[above] - upper))
        y[:, above] = np.vstack(
            (inc + ref, 1j * kt[:, None] * (inc - ref) / (n_ambient**2 if p else 1))
        )
        inside = np.flatnonzero(~(below | above))
        indices = np.minimum(
            np.searchsorted(edges, flat[inside], side="right") - 1, len(solutions) - 1
        )
        for index in np.unique(indices):
            selected = inside[indices == index]
            y[:, selected] = solutions[index](flat[selected])
        y *= np.tile(normalization, 2)[:, None]
        return y

    def permittivity(z):
        """Sample the profile with homogeneous asymptotes outside bounds."""
        flat = z.ravel()
        try:
            eps = np.broadcast_to(epsilon(flat), flat.shape).astype(complex).copy()
        except (TypeError, ValueError):
            eps = np.array([epsilon(v) for v in flat], complex)
        eps[flat < lower] = n_substrate**2
        eps[flat > upper] = n_ambient**2
        return eps

    def sample(z):
        """Return electric components at physical z in Angstrom."""
        y = states(z)
        if p:
            eps = permittivity(z)
            tangent = (n_ambient * y[count:] / (1j * k0)).T
            normal = (n_ambient * parallel[:, None] * y[:count] / eps).T
        else:
            tangent = y[:count].T
            normal = np.zeros_like(tangent)
        return tangent.reshape(z.shape + (count,)), normal.reshape(z.shape + (count,))

    def waves(z):
        """Decompose the local field and derivative into reciprocal branches."""
        y = states(z)
        eps = permittivity(z)
        kap = _kappa(eps[:, None], parallel, k0)
        if p:
            tangential = n_ambient * y[count:].T / (1j * k0 * np.sin(angles))
            difference = (
                kap * n_ambient * y[:count].T / (k0 * eps[:, None] * np.sin(angles))
            )
        else:
            tangential = y[:count].T
            difference = y[count:].T / (1j * kap)
        shape = z.shape + (count,)
        return (
            -kap.reshape(shape),
            (0.5 * (tangential + difference)).reshape(shape),
            (0.5 * (tangential - difference) * (-1 if p else 1)).reshape(shape),
        )

    return ContinuousElectricField(
        r,
        t,
        angles,
        n_ambient,
        n_substrate,
        k0,
        z_reference,
        polarization,
        "numerical",
        sample,
        waves,
    )


def solve_tanh_field(
    n_substrate,
    n_ambient,
    z0,
    sigma,
    energy_eV,
    alpha,
    polarization="s",
    *,
    method="auto",
    rtol=2e-10,
):
    """Solve an Epstein transition in epsilon=n**2, with RMS width in Angstrom.

    s polarization uses exact complex gamma/hypergeometric amplitudes and
    fields. p polarization and ``method='numerical'`` use Maxwell integration.
    Angles are radians in the ambient; incident phase is referenced to z0.
    """
    angles = np.atleast_1d(np.asarray(alpha, float))
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("Tanh sigma must be positive and finite.")
    if (
        not np.all(np.isfinite([z0, n_substrate, n_ambient, rtol]))
        or rtol <= 0
        or n_substrate == 0
        or n_ambient == 0
    ):
        raise ValueError("Finite position/indices and positive tolerance required.")
    if method not in {"auto", "analytic", "numerical"}:
        raise ValueError("method must be auto, analytic, or numerical.")
    w = np.sqrt(12.0) * sigma / np.pi
    eps_b, eps_t = complex(n_substrate) ** 2, complex(n_ambient) ** 2

    def epsilon(z):
        """Return the Epstein dielectric profile at z in Angstrom."""
        return eps_b + (eps_t - eps_b) * expit(2 * (z - z0) / w)

    if method == "numerical" or polarization == "p":
        if method == "analytic" and polarization == "p":
            raise ValueError("Exact tanh fields are available for s polarization only.")
        return solve_continuous_field(
            epsilon,
            (z0 - 18 * w, z0 + 18 * w),
            n_substrate,
            n_ambient,
            energy_eV,
            angles,
            polarization,
            z_reference=z0,
            rtol=rtol,
        )
    if (
        polarization != "s"
        or angles.ndim != 1
        or np.any((angles <= 0) | (angles > np.pi / 2))
        or not np.all(np.isfinite(angles))
    ):
        raise ValueError("Finite glancing angles and s/p polarization required.")
    if energy_eV <= 0 or not np.isfinite(energy_eV):
        raise ValueError("Photon energy must be positive and finite.")
    k0 = 2 * np.pi / (HC_KEV_ANGSTROM / (energy_eV * 1e-3))
    parallel = complex(n_ambient) * np.cos(angles)
    kb, kt = _kappa(eps_b, parallel, k0), _kappa(eps_t, parallel, k0)
    a, b = 0.5j * w * kb, 0.5j * w * kt
    log_inc = (
        loggamma(1 + 2 * a) + loggamma(2 * b) - loggamma(a + b) - loggamma(1 + a + b)
    )
    t = np.exp(-log_inc)
    if eps_b == eps_t:
        r = np.zeros_like(t)
    else:
        log_ref = (
            loggamma(1 + 2 * a)
            + loggamma(-2 * b)
            - loggamma(1 + a - b)
            - loggamma(a - b)
        )
        r = np.exp(log_ref - log_inc)

    def field_and_derivative(z):
        """Evaluate exact fields and analytic derivatives in vectorized blocks."""
        flat = z.ravel()
        values = np.empty((flat.size, len(angles)), complex)
        derivatives = np.empty_like(values)
        # Bound intermediate series arrays for long scans and dense grids.
        chunk = max(1, 65536 // len(angles))
        for first in range(0, flat.size, chunk):
            sl = slice(first, first + chunk)
            positions = flat[sl]
            v = (positions - z0) / w
            ev, dv = values[sl], derivatives[sl]
            above = (v > 18) | (eps_b == eps_t)
            below = (v < -18) & ~above
            inc = np.exp(1j * (positions[above, None] - z0) * kt)
            ref = r * np.exp(-1j * (positions[above, None] - z0) * kt)
            ev[above], dv[above] = inc + ref, 1j * kt * (inc - ref)
            ev[below] = t * np.exp(1j * (positions[below, None] - z0) * kb)
            dv[below] = 1j * kb * ev[below]
            lower = (v <= 0) & ~below & ~above
            upper = (v > 0) & ~above
            x = expit(2 * v[lower, None])
            f, df = _hyp2f1_value_derivative(a + b, 1 + a + b, 1 + 2 * a, x)
            factor = t * np.exp(a * np.log(x) + b * np.log1p(-x))
            ev[lower] = factor * f
            dv[lower] = 2 / w * factor * ((a * (1 - x) - b * x) * f + (1 - x) * df)
            u = expit(-2 * v[upper, None])
            f_inc, d_inc = _hyp2f1_value_derivative(1 + a - b, a - b, 1 - 2 * b, u)
            f_ref, d_ref = _hyp2f1_value_derivative(a + b, 1 + a + b, 1 + 2 * b, u)
            inc = np.exp(-b * np.log(u))
            ref = r * np.exp(b * np.log(u))
            prefactor = np.exp(a * np.log1p(-u))
            combined = inc * f_inc + ref * f_ref
            ev[upper] = prefactor * combined
            dv[upper] = (
                2
                / w
                * prefactor
                * (
                    a * u * combined
                    - (1 - u) * (inc * (-b * f_inc + d_inc) + ref * (b * f_ref + d_ref))
                )
            )
        shape = z.shape + (len(angles),)
        return values.reshape(shape), derivatives.reshape(shape)

    def sample(z):
        """Return electric components at physical z in Angstrom."""
        result, _ = field_and_derivative(z)
        return result, np.zeros_like(result)

    def waves(z):
        """Decompose the exact field using its analytic spatial derivative."""
        values, derivatives = field_and_derivative(z)
        kap = _kappa(np.asarray(epsilon(z))[..., None], parallel, k0)
        difference = derivatives / (1j * kap)
        return -kap, 0.5 * (values + difference), 0.5 * (values - difference)

    return ContinuousElectricField(
        r,
        t,
        angles,
        n_ambient,
        n_substrate,
        k0,
        z0,
        polarization,
        "analytic_epstein",
        sample,
        waves,
    )

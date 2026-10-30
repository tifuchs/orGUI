"""Finite Epstein segments; epsilon=A+B*expit(2*(z-center)/width).

Two analytic solutions form the field/derivative fundamental matrix. Matching
finite segments solves the declared piecewise profile, not a sum of overlapping
infinite Epstein interfaces. Coordinates and widths are in Angstrom.
"""

from dataclasses import replace

import numpy as np
from scipy.fft import dct
from scipy.special import expit, loggamma

from .CTRcontinuous import _hyp2f1_value_derivative, _kappa, solve_continuous_field


def _fundamental(z, bottom, top, center, width, parallel, k0):
    """Return [E,E'] columns for two independent analytic solutions."""
    z = np.atleast_1d(z)
    kb, kt = _kappa(bottom, parallel, k0), _kappa(top, parallel, k0)
    # Choose the larger asymptotic wavevector as the independent Frobenius
    # basis. The opposite endpoint can be exactly critical; handle that rare
    # connection-formula limit at high precision below.
    flip = abs(kb) < abs(kt)
    direction = np.where(flip, -1.0, 1.0)
    a = 0.5j * width * np.where(flip, kt, kb)
    b = 0.5j * width * np.where(flip, kb, kt)
    a = a[:, None] * np.array([1.0, -1.0])
    b = b[:, None]
    v = (z[:, None, None] - center) * direction[None, :, None] / width
    value = np.empty((len(z), len(parallel), 2), complex)
    derivative = np.empty_like(value)
    cancellation = np.zeros(value.shape, bool)
    for mirrored in (False, True):
        angles = np.flatnonzero(flip == mirrored)
        if not len(angles):
            continue
        sign = -1 if mirrored else 1
        zz = sign * (z - center) / width
        aa, bb = a[angles], b[angles]
        lower, upper = np.flatnonzero(zz <= 0), np.flatnonzero(zz > 0)
        if len(lower):
            x = expit(2 * zz[lower, None, None])
            f, df = _hyp2f1_value_derivative(aa + bb, 1 + aa + bb, 1 + 2 * aa, x)
            pref = np.exp(aa * np.log(x) + bb * np.log1p(-x))
            selected = np.ix_(lower, angles, np.arange(2))
            value[selected] = pref * f
            derivative[selected] = (
                sign * 2 / width * pref * ((aa * (1 - x) - bb * x) * f + (1 - x) * df)
            )
        if len(upper):
            u = expit(-2 * zz[upper, None, None])
            # Connect to u=1-x; every series argument remains <= 1/2.
            with np.errstate(invalid="ignore", divide="ignore"):
                ci = np.exp(
                    loggamma(1 + 2 * aa)
                    + loggamma(2 * bb)
                    - loggamma(aa + bb)
                    - loggamma(1 + aa + bb)
                )
                cr = np.exp(
                    loggamma(1 + 2 * aa)
                    + loggamma(-2 * bb)
                    - loggamma(1 + aa - bb)
                    - loggamma(aa - bb)
                )
            fi, di = _hyp2f1_value_derivative(1 + aa - bb, aa - bb, 1 - 2 * bb, u)
            fr, dr = _hyp2f1_value_derivative(aa + bb, 1 + aa + bb, 1 + 2 * bb, u)
            inc, ref = ci * np.exp(-bb * np.log(u)), cr * np.exp(bb * np.log(u))
            pref = np.exp(aa * np.log1p(-u))
            selected = np.ix_(upper, angles, np.arange(2))
            value[selected] = pref * (inc * fi + ref * fr)
            derivative[selected] = (
                sign
                * 2
                / width
                * pref
                * (
                    aa * u * (inc * fi + ref * fr)
                    - (1 - u) * (inc * (-bb * fi + di) + ref * (bb * fr + dr))
                )
            )
            cancellation[selected] = 2e-16 * abs(pref) * (
                abs(inc * fi) + abs(ref * fr)
            ) > 1e-11 * np.maximum(abs(value[selected]), 1e-300)
    bad = cancellation | ~np.isfinite(value) | ~np.isfinite(derivative)
    # At a critical asymptote the gamma connection terms separately diverge,
    # while their sum is finite. Direct 2F1 evaluates the removable limit.
    if np.any(bad):
        import mpmath as mp

        with mp.workdps(50):
            for iz, ia, branch in zip(*np.nonzero(bad)):
                aa, bb = complex(a[ia, branch]), complex(b[ia, 0])
                vv = mp.mpf(float(v[iz, ia, 0]))
                xxx = 1 / (1 + mp.exp(-2 * vv))
                fff = mp.hyp2f1(aa + bb, 1 + aa + bb, 1 + 2 * aa, xxx)
                dff = (
                    (aa + bb)
                    * (1 + aa + bb)
                    / (1 + 2 * aa)
                    * mp.hyp2f1(1 + aa + bb, 2 + aa + bb, 2 + 2 * aa, xxx)
                )
                pref = xxx**aa * (1 - xxx) ** bb
                value[iz, ia, branch] = complex(pref * fff)
                derivative[iz, ia, branch] = (
                    complex(
                        2
                        / width
                        * pref
                        * ((aa * (1 - xxx) - bb * xxx) * fff + xxx * (1 - xxx) * dff)
                    )
                    * direction[ia]
                )
    return np.stack((value, derivative), axis=-2)


def _solve_piecewise(
    epsilon,
    edges,
    segments,
    n_bottom,
    n_top,
    energy,
    alpha,
    polarization="s",
    *,
    rtol=2e-10,
):
    """Match finite analytic segments, or solve their p equation numerically.

    Each segment is (A, B, center, width) for epsilon=1+A+B*logistic.
    The conserved parallel wavevector always uses the true upper ambient.
    """

    def constant(left, right):
        """Return the exact epsilon of a homogeneous interval, or None."""
        index = np.searchsorted(edges, (left + right) / 2) - 1
        a, b, _, _ = segments[index]
        return 1 + a if b == 0 else None

    def segment(left, right, state, parallel, k0):
        """Construct the analytic state propagator on one finite segment."""
        index = np.searchsorted(edges, (left + right) / 2) - 1
        a, b, center, width = segments[index]
        matrix = _fundamental([left], 1 + a, 1 + a + b, center, width, parallel, k0)[0]
        count = len(parallel)
        coefficients = np.linalg.solve(matrix, state.reshape(2, count).T[..., None])

        def exact(z):
            """Evaluate matched field and derivative at Angstrom coordinates."""
            z = np.asarray(z)
            matrix = _fundamental(
                z.ravel(), 1 + a, 1 + a + b, center, width, parallel, k0
            )
            values = (matrix @ coefficients)[..., 0]
            return values.transpose(2, 1, 0).reshape((2 * count,) + z.shape)

        endpoint = exact(right)
        expansion = None

        def propagate(z):
            """Sample the analytic segment with checked interpolation."""
            # Matching uses exact endpoint values. Repeated quadrature/atomic
            # samples use a checked spectral expansion of those analytic fields,
            # not an ODE or a sliced approximation to the optical potential.
            nonlocal expansion
            z = np.asarray(z)
            if z.ndim == 0 and z == right:
                return endpoint
            if expansion is None:
                midpoint, half = (left + right) / 2, (right - left) / 2
                check_x = np.linspace(-1, 1, 11)
                check = exact(midpoint + half * check_x)
                scale = np.maximum(np.max(abs(check), axis=1), 1e-300)
                for order in (16, 32, 64, 128, 256):
                    nodes = np.cos(np.pi * (np.arange(order) + 0.5) / order)
                    values = exact(midpoint + half * nodes)
                    trial = (dct(values, type=2, axis=1) / order).T.copy()
                    trial[0] *= 0.5
                    tail = np.sum(abs(trial[-4:]), axis=0)
                    error = np.max(
                        abs(np.polynomial.chebyshev.chebval(check_x, trial) - check),
                        axis=1,
                    )
                    if np.all(tail < 2e-13 * scale) and np.all(error < 2e-12 * scale):
                        expansion = trial
                        break
                if expansion is None:
                    # Pathological widths/angles retain direct analytic fields.
                    # Never hide a failed interpolation behind looser accuracy.
                    return exact(z)
            x = (2 * z.ravel() - left - right) / (right - left)
            # Build each polynomial only once per coordinate. Clenshaw's
            # broadcast recurrence otherwise repeats it for every angle/field
            # component and allocates a full sample-by-component array per term.
            basis = np.polynomial.chebyshev.chebvander(x, len(expansion) - 1)
            # A non-BLAS contraction avoids worker-thread overhead for the many
            # small segments encountered with adsorbed layers and terraces.
            return np.einsum("dc,zd->cz", expansion, basis, optimize=False).reshape(
                (2 * count,) + z.shape
            )

        return propagate

    result = solve_continuous_field(
        epsilon,
        (edges[0], edges[-1]),
        n_bottom,
        n_top,
        energy,
        alpha,
        polarization,
        z_reference=edges[-1],
        interfaces=edges,
        rtol=rtol,
        constant_permittivity=constant,
        _segment_propagator=segment if polarization == "s" else None,
    )
    return replace(
        result,
        solver=(
            "analytic_piecewise_epstein"
            if polarization == "s"
            else "numerical_piecewise_p"
        ),
    )

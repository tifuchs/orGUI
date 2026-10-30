"""Exact Epstein amplitudes checked against independent Maxwell integration."""

import numpy as np
import pytest

from ..CTRcontinuous import solve_tanh_field
from ..CTRoptics import solve_electric_field


@pytest.mark.parametrize("absorption", [0.0, 2e-7])
def test_epstein_complex_fields_against_ode(absorption):
    """Gamma-function phases and hypergeometric fields match a boundary solve."""
    nb, nt = 1 - 8e-6 - 1j * absorption, 1 - 2e-6 - 1j * absorption / 3
    angles = np.array([0.002, 0.006, 0.03])
    args = (nb, nt, 1.3, 0.8, 10000.0, angles)
    analytic = solve_tanh_field(*args)
    numerical = solve_tanh_field(*args, method="numerical")
    np.testing.assert_allclose(analytic.r_S, numerical.r_S, rtol=1e-7, atol=2e-9)
    np.testing.assert_allclose(analytic.t_S, numerical.t_S, rtol=1e-7, atol=2e-9)
    z = np.linspace(-5.0, 7.0, 21)
    np.testing.assert_allclose(
        analytic.sample(z)[0], numerical.sample(z)[0], rtol=1e-7, atol=2e-9
    )


def test_epstein_lossless_intensity_and_flux():
    """Independent sinh reflectivity and conserved flux check the amplitudes."""
    nb, nt = 1 - 8e-6, 1 - 2e-6
    angles = np.array([0.01, 0.03, 0.1])
    sigma = 1.2
    f = solve_tanh_field(nb, nt, 0.0, sigma, 10000.0, angles)
    kt = f.k0 * nt * np.sin(angles)
    kb = f.k0 * np.sqrt(nb * nb - (nt * np.cos(angles)) ** 2)
    width = np.sqrt(12) * sigma / np.pi
    expected = (
        np.sinh(np.pi * width * (kt - kb) / 2) / np.sinh(np.pi * width * (kt + kb) / 2)
    ) ** 2
    np.testing.assert_allclose(abs(f.r_S) ** 2, expected, rtol=1e-9)
    np.testing.assert_allclose(
        abs(f.r_S) ** 2 + kb / kt * abs(f.t_S) ** 2, 1.0, atol=1e-12
    )


@pytest.mark.parametrize("pol", ["s", "p"])
def test_tanh_sharp_limit_matches_native(pol):
    """Both polarization equations recover the existing water-ambient solver."""
    nb, nt = 1 - 8e-6 - 2e-7j, 1 - 2e-6 - 1e-8j
    angles = np.array([0.002, 0.01, 0.1])
    f = solve_tanh_field(nb, nt, 0.0, 1e-4, 10000.0, angles, pol)
    profile = np.array([[-1.0, 1 - nb.real, -nb.imag], [1.0, 1 - nt.real, -nt.imag]])
    native = solve_electric_field(profile, 10000.0, angles, pol)
    np.testing.assert_allclose(f.r_S, native.r_S, rtol=1e-6, atol=2e-8)


@pytest.mark.parametrize("pol", ["s", "p"])
def test_uniform_medium_and_reference_translation(pol):
    """No-contrast water has no reflection; translating an interface fixes phase."""
    n = 1 - 2e-6 - 1e-8j
    angle = np.array([0.01, 0.07])
    uniform = solve_tanh_field(n, n, 2.0, 0.6, 10000.0, angle, pol)
    np.testing.assert_allclose(uniform.r_S, 0.0, atol=1e-10)
    np.testing.assert_allclose(uniform.t_S, 1.0, atol=1e-9)
    original = solve_tanh_field(n - 5e-6, n, 0.0, 0.6, 10000.0, angle, pol)
    shifted = solve_tanh_field(n - 5e-6, n, 3.0, 0.6, 10000.0, angle, pol)
    np.testing.assert_allclose(shifted.r_S, original.r_S, rtol=1e-8, atol=1e-10)
    # Each field is normalized at its own translated reference plane.
    np.testing.assert_allclose(
        shifted.sample(np.array([1.0, 3.0, 5.0]))[0],
        original.sample(np.array([-2.0, 0.0, 2.0]))[0],
        rtol=1e-8,
        atol=1e-9,
    )


def test_exact_tanh_fields_do_not_run_an_ode(monkeypatch):
    """Exact s coefficients, internal fields and derivatives need no integration."""
    from .. import CTRcontinuous

    def forbidden(*args, **kwargs):
        """Fail if the exact path requests numerical propagation."""
        raise AssertionError("The analytic tanh path must not integrate an ODE")

    monkeypatch.setattr(CTRcontinuous, "solve_ivp", forbidden)
    field = solve_tanh_field(
        1 - 8e-6 - 2e-7j, 1 - 2e-6 - 1e-8j, 0.0, 1.0, 10000.0, [0.002, 0.07]
    )
    z = np.linspace(-25.0, 25.0, 101)
    values = field.sample(z)[0]
    _, forward, backward = field.plane_waves(z)
    np.testing.assert_allclose(forward + backward, values, atol=1e-13)


def test_vectorized_hypergeometric_derivative_against_high_precision():
    """Series and its analytic derivative match independent arbitrary precision."""
    import mpmath as mp
    from ..CTRcontinuous import _hyp2f1_value_derivative

    a = np.array([0.2 + 0.7j, 2.0 + 5.0j, 0.1 + 0.02j])
    b, c = a + 1, 1 + 0.7 * a
    x = np.array([1e-8, 0.1, 0.49, 0.5])[:, None]
    value, derivative = _hyp2f1_value_derivative(a, b, c, x)
    with mp.workdps(50):
        for i, xx in enumerate(x[:, 0]):
            for j in range(len(a)):
                aa, bb, cc = map(complex, (a[j], b[j], c[j]))
                expected = complex(mp.hyp2f1(aa, bb, cc, float(xx)))
                expected_d = complex(
                    xx * mp.diff(lambda v: mp.hyp2f1(aa, bb, cc, v), mp.mpf(float(xx)))
                )
                np.testing.assert_allclose(
                    value[i, j], expected, rtol=2e-12, atol=1e-14
                )
                np.testing.assert_allclose(
                    derivative[i, j], expected_d, rtol=2e-12, atol=1e-14
                )


@pytest.mark.parametrize("pol", ["s", "p"])
def test_homogeneous_segment_propagation_matches_ode(pol):
    """Exact homogeneous transfer preserves continuous fields and p admittance."""
    from ..CTRcontinuous import solve_continuous_field

    eps = (1 - 3e-6 - 2e-8j) ** 2
    args = (
        lambda z: eps,
        (-6.0, 4.0),
        1 - 8e-6 - 1e-7j,
        1 - 2e-6 - 1e-8j,
        10000.0,
        [0.002, 0.05],
        pol,
    )
    ordinary = solve_continuous_field(*args)
    fast = solve_continuous_field(*args, constant_permittivity=lambda a, b: eps)
    np.testing.assert_allclose(fast.r_S, ordinary.r_S, rtol=1e-7, atol=1e-9)
    np.testing.assert_allclose(
        fast.sample(np.linspace(-10, 8, 20))[0],
        ordinary.sample(np.linspace(-10, 8, 20))[0],
        atol=1e-8,
    )


@pytest.mark.parametrize("absorption", [0.0, 2e-7])
def test_piecewise_epstein_fields_and_matching(absorption, monkeypatch):
    """Finite, overlapping-width segments match an independent Maxwell solve."""
    from scipy.special import expit
    from .. import CTRcontinuous
    from .._CTRepstein import _solve_piecewise

    bottom, top = (1 - 8e-6 - 1j * absorption) ** 2, (1 - 2e-6) ** 2
    edges = np.array([-9.0, -0.4, 1.2, 9.0])
    centers = np.array([-1.0, 0.5, 2.0])
    weights = np.array([0.2, 0.5, 0.3])
    width = 1.3

    def target(z):
        """Evaluate a physical overlapping mixture for endpoint matching."""
        return bottom + (top - bottom) * np.sum(
            weights * expit(2 * (np.asarray(z)[..., None] - centers) / width), axis=-1
        )

    segments = []
    for left, right, center in zip(edges[:-1], edges[1:], centers):
        x = expit(2 * (np.array([left, right]) - center) / width)
        b = (target(right) - target(left)) / (x[1] - x[0])
        segments.append((target(left) - 1 - b * x[0], b, center, width))

    def epsilon(z):
        """Evaluate the declared piecewise potential independently."""
        z = np.asarray(z)
        i = np.clip(np.searchsorted(edges, z, side="right") - 1, 0, 2)
        a, b, center, w = np.asarray(segments)[i].T
        return 1 + a + b * expit(2 * (z - center.real) / w.real)

    angles = np.array([0.001, 0.0035, 0.01, 0.15])
    args = (epsilon, edges, segments, np.sqrt(bottom), np.sqrt(top), 10000.0, angles)
    numerical = CTRcontinuous.solve_continuous_field(
        epsilon,
        (edges[0], edges[-1]),
        np.sqrt(bottom),
        np.sqrt(top),
        10000.0,
        angles,
        interfaces=edges,
        z_reference=edges[-1],
    )

    def forbidden(*args, **kwargs):
        """Reject any ODE on the piecewise analytic branch."""
        raise AssertionError("Piecewise s fields must not run an ODE")

    monkeypatch.setattr(CTRcontinuous, "solve_ivp", forbidden)
    field = _solve_piecewise(*args)
    np.testing.assert_allclose(field.r_S, numerical.r_S, rtol=2e-6, atol=3e-9)
    z = np.r_[np.linspace(-12, 12, 117), edges[1:-1] - 1e-8, edges[1:-1] + 1e-8]
    np.testing.assert_allclose(
        field.sample(z)[0], numerical.sample(z)[0], rtol=2e-8, atol=3e-9
    )
    # Repeated samples exercise the validated analytic-field interpolation.
    np.testing.assert_array_equal(field.sample(z)[0], field.sample(z)[0])
    if not absorption:
        kt = field.k0 * np.sqrt(top) * np.sin(angles)
        kb = CTRcontinuous._kappa(bottom, np.sqrt(top) * np.cos(angles), field.k0)
        np.testing.assert_allclose(
            abs(field.r_S) ** 2 + kb.real / kt.real * abs(field.t_S) ** 2,
            1.0,
            atol=2e-10,
        )


def test_epstein_fundamental_at_critical_asymptote():
    """The removable gamma singularity must give finite independent fields."""
    from .._CTRepstein import _fundamental

    matrix = _fundamental(
        np.array([-1.0, 0.0, 1.0]), 0.99, 1.0, 0.0, 1.0, np.array([1.0 + 0j]), 2.0
    )
    assert np.all(np.isfinite(matrix))
    determinant = np.linalg.det(matrix)
    np.testing.assert_allclose(determinant, determinant[0, 0], atol=1e-12)
    assert abs(determinant[0, 0]) > 0.1

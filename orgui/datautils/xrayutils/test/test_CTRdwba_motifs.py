"""Independent atom-sum and work-count checks for native lateral motifs."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from .. import CTRuc
from .test_CTRdwba import _synthetic_record_call
from .test_CTRincoherent_dwba import _hkl, _model


@pytest.mark.parametrize("signed", [False, True])
def test_lateral_motif_matches_direct_complex_phase_sum(signed):
    """A motif retains every weighted lateral site and all four wave channels."""
    count = 64
    rng = np.random.default_rng(753)
    positions = np.column_stack((rng.normal(size=(count, 2)),
                                 np.full(count, -.43)))
    weights = rng.uniform(.2, .8, count)
    if signed:
        weights[::2] *= -1
    q = np.array([.37, -.19])
    ki, kf = -1.07 + .09j, -1.18 + .07j
    ai, af = [.73 - .04j, .21 + .03j], [.69 + .02j, -.12 + .09j]
    basis = np.array([[6., 0., 0., 0., .23, .31, .67]])
    factors = np.zeros((1, 13))
    factors[0, [0, 5, 10, 11, 12]] = [1.7, .42, .8, .13, .21]
    actual = _synthetic_record_call(
        q_parallel=q, kz_i=[-1.21 + .04j, ki], kz_f=[-1.46 + .03j, kf],
        A_i_plus=[.81 + .06j, ai[0]], A_i_minus=[-.17 + .08j, ai[1]],
        A_f_plus=[.77 - .05j, af[0]], A_f_minus=[.14 + .07j, af[1]],
        basis=basis, factors=factors, finite_positions=positions,
        finite_records=np.ones(count, dtype=np.int64), finite_weights=weights,
    )[0][1, 0]
    expected = 0j
    for i, si in enumerate((1, -1)):
        for j, sf in enumerate((1, -1)):
            qz = -(si * ki + sf * kf)
            factor = 1.7 * np.exp(-.42 * (q @ q + qz * qz)) + .93 + .21j
            factor *= np.exp(-(.23 * (q @ q) + .31 * qz * qz)
                             / (16 * np.pi**2)) * .67
            for position, weight in zip(positions, weights):
                expected += (weight * ai[i] * af[j] * np.cos(.41) * factor
                             * np.exp(1j * (q @ position[:2] + qz * position[2])))
    assert_allclose(actual, expected, rtol=3e-14, atol=1e-28)
    work = CTRuc._CTRcalc_cpp._dwba_motif_work()
    assert work["instances"] == count
    assert work["motifs"] == 1
    assert work["motif_channels"] == 4
    assert work["normal_exponentials"] == 4
    assert work["lateral_exponentials"] == count


def test_motif_keys_keep_records_and_near_equal_scattering_states():
    """No tolerance may merge distinct depths, occupancies or DWF parameters."""
    basis = np.tile([6., 0., 0., 0., .2, .3, .7], (4, 1))
    basis[1, 4] = np.nextafter(.2, 1.)
    basis[2, 6] = np.nextafter(.7, 1.)
    positions = np.tile([.13, .27, -.5], (5, 1))
    positions[3, 2] = np.nextafter(-.5, 0.)
    factors = np.zeros((4, 13))
    factors[:, 10] = 1.
    _synthetic_record_call(
        q_parallel=[.17, .08], kz_i=[-1 + .02j, -1.3 + .05j],
        kz_f=[-1.1 + .03j, -1.4 + .06j],
        A_i_plus=[1., .8], A_i_minus=[.1, .2],
        A_f_plus=[1., .7], A_f_minus=[.2, .1],
        reference_shares=np.zeros((3, 2, 2)), basis=basis, factors=factors,
        finite_positions=positions, finite_atoms=np.array([0, 1, 2, 3, 0]),
        finite_records=np.array([1, 1, 1, 1, 2]),
    )
    assert CTRuc._CTRcalc_cpp._dwba_motif_work()["motifs"] == 5


def test_repeated_coordinates_reuse_lateral_work():
    """Identical in-plane momenta require just one set of lateral exponentials."""
    crystal = _model().coherent_model
    with crystal.dwba.batch():
        crystal.dwba.evaluate(1., 0., .45)
        single = CTRuc._CTRcalc_cpp._dwba_motif_work()
        crystal.dwba.evaluate(1., 0., np.full(32, .45))
        repeated = CTRuc._CTRcalc_cpp._dwba_motif_work()
    assert repeated["lateral_exponentials"] == single["lateral_exponentials"]
    assert repeated["normal_exponentials"] == 32 * single["normal_exponentials"]
    assert repeated["motif_channels"] == 32 * single["motif_channels"]


@pytest.mark.parametrize("channels", [("s", "s"), ("s", "p"),
                                      ("p", "s"), ("p", "p")])
def test_changing_in_plane_vector_matches_separate_calls(channels):
    """Reuse on a rod must refresh at every change of its in-plane vector."""
    model = _model(policy="state_specific")
    points = np.array(np.broadcast_arrays(*_hkl())).T
    # Repeat and revisit a rod, including specular and near-Bragg points.
    points = points[[0, 0, 2, 1, 2, 6]]
    options = dict(polarization_i=channels[0], polarization_f=channels[1],
                   bulk_attenuation=.005)
    with model.dwba.batch():
        actual = model.dwba.evaluate(*points.T, **options)
        separate = [model.dwba.evaluate(*point, **options) for point in points]
    for name in ("reflectivity", "coherent_intensity", "incoherent_intensity"):
        assert_allclose(getattr(actual, name), [getattr(r, name) for r in separate],
                        rtol=1e-12, atol=1e-28)

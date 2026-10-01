"""Count-statistical propagation through a split reciprocal-space integral."""

import numpy as np
import pytest

native = pytest.importorskip(
    "orgui.datautils.xrayutils._reciprocal_reconstruction_cpp"
)


def _dense_indices(batch, shape, chunk):
    chunks = np.ceil(shape / chunk).astype(int)
    chunk_index = np.array(np.unravel_index(batch["chunk_id"], chunks)).T
    local = np.array(np.unravel_index(batch["local_voxel_id"], chunk)).T
    return tuple((chunk_index * chunk + local).T)


@pytest.mark.parametrize("depth", [0, 2, 3])
@pytest.mark.parametrize("moving", [False, True])
@pytest.mark.parametrize(
    "weighting", ["parameter_average", "reciprocal_volume_average"]
)
def test_backprojection_matches_independent_pixel_impulses(depth, moving, weighting):
    """The adjoint preserves each pixel's full contribution before variance squaring."""
    if weighting == "reciprocal_volume_average" and (not moving or depth == 0):
        pytest.skip("Volume weighting needs a moving exposure and cell corners")
    shape = np.array([17, 19, 13])
    chunk = np.array([4, 8, 4])
    kernel = native.ReconstructionKernel(
        np.full(3, -.2), np.full(3, .025), shape, chunk, "hkl", 1.,
        np.eye(3), np.eye(3), depth, 1, 16, 4 * 1024**2, weighting,
    )
    rays = np.zeros((3, 4, 3))
    rays[..., 0] = np.linspace(-.075, .075, 4)[None, :]
    rays[..., 1] = 1.
    rays[..., 2] = np.linspace(-.06, .06, 3)[:, None]
    rays /= np.linalg.norm(rays, axis=2, keepdims=True)
    start = np.array([.08, -.025, .02, -.02])
    end = start.copy()
    if moving:
        end[1] += .05
    mask = np.zeros((2, 3), dtype=bool)
    mask[0, 2] = True
    coefficients = np.random.default_rng(82).uniform(.5, 2., size=shape)
    back = kernel.backproject_coefficients(coefficients, mask, rays, start, end)
    expected = np.zeros_like(back)
    for row, column in np.ndindex(mask.shape):
        impulse = np.zeros(mask.shape)
        impulse[row, column] = 1.
        batch = kernel.accumulate(
            impulse, np.zeros_like(impulse), mask, rays, start, end
        )
        expected[row, column] = np.sum(
            coefficients[_dense_indices(batch, shape, chunk)]
            * batch["weighted_intensity"]
        )
    np.testing.assert_allclose(back, expected, rtol=2e-12, atol=1e-18)
    assert back[0, 2] == 0
    assert np.count_nonzero(back) > 0
    counts = np.array([[10., 20., 30.], [8., 12., 4.]])
    batch = kernel.accumulate(counts, counts, mask, rays, start, end)
    functional = coefficients[_dense_indices(batch, shape, chunk)]
    np.testing.assert_allclose(np.sum(back * counts),
                               np.sum(functional * batch["weighted_intensity"]),
                               rtol=2e-12)
    variance = np.sum(back**2 * counts)
    marginal_sum = np.sum(functional**2 * batch["weighted_variance"])
    assert variance >= marginal_sum * (1 - 2e-12)
    if depth > 0:
        assert variance > marginal_sum * 1.05  # actual cross-voxel covariance


def test_backprojection_rejects_nonfinite_coefficients():
    """An invalid functional cannot silently produce invalid count uncertainties."""
    kernel = native.ReconstructionKernel(
        np.zeros(3), np.ones(3), np.ones(3, dtype=int),
        np.ones(3, dtype=int), "lab", 1., np.eye(3), np.eye(3),
    )
    rays = np.zeros((2, 2, 3))
    rays[..., 1] = 1.
    with pytest.raises(ValueError, match="finite"):
        kernel.backproject_coefficients(
            np.full((1, 1, 1), np.nan), np.zeros((1, 1), dtype=bool),
            rays, np.zeros(4), np.zeros(4),
        )

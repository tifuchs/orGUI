"""Tests for the independent reciprocal-volume validation oracle."""

import numpy as np
import pytest

from benchmarks.reconstruction_volume_oracle import (
    _dyadic_leaf_oracle,
    _gauss_legendre_volume,
)


def test_affine_oracles_reproduce_the_analytic_determinant():
    matrix = np.array(
        [
            [1.2, -0.3, 0.2],
            [0.1, 0.9, -0.4],
            [-0.2, 0.3, 1.4],
        ]
    )
    offset = np.array([2.0, -1.0, 0.5])

    def coordinate(u, v, t):
        return offset + matrix @ np.array([u, v, t])

    expected = abs(float(np.linalg.det(matrix)))
    leaves = _dyadic_leaf_oracle(coordinate, 3)
    quadrature = _gauss_legendre_volume(coordinate, order=3)

    assert leaves["root"]["volume"] == pytest.approx(expected, rel=1e-14)
    assert leaves["leaf_total_volume"] == pytest.approx(expected, rel=1e-14)
    assert quadrature["volume"] == pytest.approx(expected, rel=1e-10)
    assert not leaves["orientation_sign_change"]
    assert not quadrature["orientation_sign_change"]


def test_leaf_oracle_converges_for_a_curved_mapping():
    curvature = 0.8

    def coordinate(u, v, t):
        return np.array([u, v, t * (1.0 + curvature * u**2)])

    exact = 1.0 + curvature / 3.0
    depth_one = _dyadic_leaf_oracle(coordinate, 1)
    depth_three = _dyadic_leaf_oracle(coordinate, 3)
    quadrature = _gauss_legendre_volume(coordinate, order=3)
    root_error = abs(depth_one["root"]["volume"] - exact)
    depth_one_error = abs(depth_one["leaf_total_volume"] - exact)
    depth_three_error = abs(depth_three["leaf_total_volume"] - exact)

    assert depth_one_error < root_error
    assert depth_three_error < depth_one_error
    assert quadrature["volume"] == pytest.approx(exact, rel=1e-9)


def test_oracle_detects_a_fold_hidden_by_the_root_corners():
    def coordinate(u, v, t):
        return np.array([(u - 0.5) ** 2, v, t])

    leaves = _dyadic_leaf_oracle(coordinate, 3)
    quadrature = _gauss_legendre_volume(coordinate, order=4)

    assert leaves["root"]["volume"] == pytest.approx(0.0, abs=1e-15)
    assert leaves["leaf_total_volume"] == pytest.approx(0.5, rel=1e-14)
    assert leaves["orientation_sign_change"]
    # A finite Gauss rule cannot integrate the absolute-value cusp exactly,
    # but it must recover the folded volume rather than the root's zero.
    assert quadrature["volume"] == pytest.approx(0.5, rel=0.05)
    assert quadrature["orientation_sign_change"]


def test_volume_is_invariant_under_output_rotation():
    matrix = np.diag([0.7, 1.3, 2.1])
    angle = 0.71
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    def original(u, v, t):
        return matrix @ np.array([u, v, t])

    def rotated(u, v, t):
        return rotation @ original(u, v, t)

    original_volume = _dyadic_leaf_oracle(original, 2)["leaf_total_volume"]
    rotated_volume = _dyadic_leaf_oracle(rotated, 2)["leaf_total_volume"]
    assert rotated_volume == pytest.approx(original_volume, rel=1e-14)


def test_reversing_scan_direction_changes_sign_but_not_volume():
    matrix = np.array(
        [
            [1.0, 0.2, -0.1],
            [0.1, 1.3, 0.4],
            [-0.2, 0.1, 0.9],
        ]
    )

    def forward(u, v, t):
        return matrix @ np.array([u, v, t])

    def reverse(u, v, t):
        return matrix @ np.array([u, v, 1.0 - t])

    forward_leaves = _dyadic_leaf_oracle(forward, 2)
    reverse_leaves = _dyadic_leaf_oracle(reverse, 2)
    forward_quadrature = _gauss_legendre_volume(forward, order=3)
    reverse_quadrature = _gauss_legendre_volume(reverse, order=3)

    assert reverse_leaves["leaf_total_volume"] == pytest.approx(
        forward_leaves["leaf_total_volume"], rel=1e-14
    )
    assert np.sum(reverse_leaves["leaf_signed_volumes"]) == pytest.approx(
        -np.sum(forward_leaves["leaf_signed_volumes"]), rel=1e-14
    )
    assert reverse_quadrature["volume"] == pytest.approx(
        forward_quadrature["volume"], rel=1e-10
    )
    assert reverse_quadrature["signed_volume"] == pytest.approx(
        -forward_quadrature["signed_volume"], rel=1e-10
    )


def test_linear_coordinate_transform_scales_volume_by_its_determinant():
    q_matrix = np.array(
        [
            [1.1, 0.2, -0.1],
            [0.0, 0.8, 0.3],
            [0.2, -0.2, 1.4],
        ]
    )
    q_to_hkl = np.diag([0.5, 0.25, 2.0])

    def q_coordinate(u, v, t):
        return q_matrix @ np.array([u, v, t])

    def hkl_coordinate(u, v, t):
        return q_to_hkl @ q_coordinate(u, v, t)

    q_volume = _dyadic_leaf_oracle(q_coordinate, 2)["leaf_total_volume"]
    hkl_volume = _dyadic_leaf_oracle(hkl_coordinate, 2)[
        "leaf_total_volume"
    ]
    assert hkl_volume == pytest.approx(
        abs(np.linalg.det(q_to_hkl)) * q_volume, rel=1e-14
    )


def test_stationary_mapping_has_zero_three_dimensional_volume():
    def coordinate(u, v, _t):
        return np.array([u, v, u + v])

    leaves = _dyadic_leaf_oracle(coordinate, 2)
    quadrature = _gauss_legendre_volume(coordinate, order=2)
    assert leaves["root"]["volume"] == pytest.approx(0.0, abs=1e-15)
    assert leaves["leaf_total_volume"] == pytest.approx(0.0, abs=1e-15)
    assert quadrature["volume"] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("depth", [-1, 6])
def test_leaf_oracle_rejects_unbounded_depth(depth):
    with pytest.raises(ValueError, match="depth"):
        _dyadic_leaf_oracle(lambda u, v, t: (u, v, t), depth)

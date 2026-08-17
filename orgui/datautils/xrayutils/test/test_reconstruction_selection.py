"""Tests for automatic reciprocal-space volume selection."""

from types import SimpleNamespace

import numpy as np
import pytest

import orgui.datautils.xrayutils.reconstruction as reconstruction_module
from orgui.datautils.xrayutils import CTRcalc
from orgui.datautils.xrayutils.reconstruction import _ReconstructionSpec
from orgui.reconstruction_selection import (
    derive_bragg_grids,
    derive_ctr_grids,
    hkl_to_frame_matrix,
    sample_hkl_coverage,
)


def _unit_cell():
    """Primitive cubic P1 cell: every integer reflection is allowed."""
    cell = CTRcalc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
    cell.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    cell.setEnergy(70e3)
    return cell


class _UBCalculator:
    """Identity orientation with a diagonal, anisotropic B matrix."""

    @staticmethod
    def getK():
        return 1.0

    @staticmethod
    def getUB():
        return np.diag([2.0, 2.0, 0.5])

    @staticmethod
    def getU():
        return np.eye(3)


def _config():
    return SimpleNamespace(
        detector=_Detector(),
        unit_cell=_unit_cell(),
        ub_calculator=_UBCalculator(),
        corrections=SimpleNamespace(excluded_frames=()),
    )


class _Detector:
    detector = SimpleNamespace(shape=(3, 4))

    @staticmethod
    def primBeamPoints(rows, columns):
        return np.zeros_like(rows), np.zeros_like(columns)


class _Scan:
    def __len__(self):
        return 3

    @staticmethod
    def exposure_angle_bounds(config, fallback):
        bounds = np.zeros((3, 2, 4), dtype=np.float64)
        bounds[:, 1, 0] = 1.0
        return bounds


class _Kernel:
    """Reports the exposure parameter and the angle sweep as coordinates."""

    def __init__(self, *args):
        self.frame = args[4]

    @staticmethod
    def coordinate(rays, start, end, row, column, u, v, t):
        angle = start[0] + t * (end[0] - start[0])
        return np.asarray((rays[0, 0, 0], rays[0, 0, 1], 4.0 * angle))


def _grid_cloud(h_range, k_range, l_range, count=41):
    """Dense hkl cloud filling the given closed ranges.

    Dense enough that every feature box below contains interior samples, so
    that the reachability test never turns on a sample landing exactly on a
    box edge.
    """
    axes = [
        np.linspace(lower, upper, count)
        for lower, upper in (h_range, k_range, l_range)
    ]
    return np.stack(
        [values.ravel() for values in np.meshgrid(*axes, indexing="ij")],
        axis=-1,
    )


def test_hkl_frame_matrix_is_identity():
    """The hkl grid frame is r.l.u. already, so no transform applies."""
    np.testing.assert_allclose(
        hkl_to_frame_matrix(_UBCalculator(), "hkl"), np.eye(3)
    )


def test_crystal_frame_matrix_matches_kernel_convention():
    """q_crystal = U^-1 UB h, the inverse of the kernel's own transform."""
    matrix = hkl_to_frame_matrix(_UBCalculator(), "crystal")
    np.testing.assert_allclose(matrix, np.diag([2.0, 2.0, 0.5]))
    # Round trip through the kernel's hkl = UB^-1 q and crystal = U^-1 q.
    hkl = np.asarray([1.0, -2.0, 3.0])
    q_lab = _UBCalculator.getUB() @ hkl
    np.testing.assert_allclose(
        matrix @ hkl, np.linalg.inv(_UBCalculator.getU()) @ q_lab
    )


def test_q_crystal_is_accepted_as_a_synonym():
    """Grid specs spell the crystal frame either way."""
    np.testing.assert_allclose(
        hkl_to_frame_matrix(_UBCalculator(), "q_crystal"),
        hkl_to_frame_matrix(_UBCalculator(), "crystal"),
    )


@pytest.mark.parametrize("frame", ["lab", "omega", "phi", "chi", "alpha"])
def test_sample_rotating_frames_are_rejected(frame):
    """A fixed hkl feature has no fixed box in a sample-rotating frame."""
    with pytest.raises(ValueError, match="does not rotate with the sample"):
        hkl_to_frame_matrix(_UBCalculator(), frame)


def test_coverage_samples_every_included_frame_and_exposure_endpoint(
    monkeypatch,
):
    """Both exposure endpoints of every included frame reach the cloud."""
    monkeypatch.setattr(
        reconstruction_module,
        "_native_module",
        lambda: SimpleNamespace(ReconstructionKernel=_Kernel),
    )
    config = _config()
    coverage = sample_hkl_coverage(
        config, _Scan(), detector_samples=3, frame_samples=8
    )
    # 3 rows x 3 sampled columns of pixels, 3 frames, 2 exposure endpoints.
    assert coverage.shape == (3 * 3 * 3 * 2, 3)
    # The mock kernel maps t onto the third axis as 4 * angle, angle in [0, 1].
    np.testing.assert_allclose(np.unique(coverage[:, 2]), [0.0, 4.0])


def test_excluded_frames_are_left_out_of_the_coverage(monkeypatch):
    """Frames the user excluded cannot select a feature."""
    monkeypatch.setattr(
        reconstruction_module,
        "_native_module",
        lambda: SimpleNamespace(ReconstructionKernel=_Kernel),
    )
    config = _config()
    config.corrections = SimpleNamespace(excluded_frames=(0, 2))
    coverage = sample_hkl_coverage(
        config, _Scan(), detector_samples=3, frame_samples=8
    )
    assert coverage.shape == (3 * 3 * 1 * 2, 3)


def test_ctr_grids_are_columns_over_the_measured_l_range():
    """Each rod spans the measured L and its own H/K half-widths."""
    coverage = _grid_cloud((-1.2, 1.2), (-1.2, 1.2), (0.2, 3.4))
    grids = derive_ctr_grids(
        _config(),
        _Scan(),
        step=(0.01, 0.01, 0.005),
        half_width=(0.1, 0.15, 0.0),
        coverage=coverage,
    )
    assert [grid.name for grid in grids] == [
        "ctr_m1_m1",
        "ctr_m1_0",
        "ctr_m1_1",
        "ctr_0_m1",
        "ctr_0_0",
        "ctr_0_1",
        "ctr_1_m1",
        "ctr_1_0",
        "ctr_1_1",
    ]
    rod = next(grid for grid in grids if grid.name == "ctr_1_m1")
    np.testing.assert_allclose(rod.minimum, (0.9, -1.15, 0.2))
    np.testing.assert_allclose(rod.maximum, (1.1, -0.85, 3.4))
    assert rod.frame == "hkl"
    assert all(grid.step == (0.01, 0.01, 0.005) for grid in grids)


def test_ctr_half_width_extends_the_measured_l_range():
    """A non-zero third half-width pads the rod beyond the measured L."""
    coverage = _grid_cloud((-0.2, 0.2), (-0.2, 0.2), (1.0, 2.0))
    (rod,) = derive_ctr_grids(
        _config(),
        _Scan(),
        step=(0.01, 0.01, 0.01),
        half_width=(0.05, 0.05, 0.25),
        coverage=coverage,
    )
    assert rod.name == "ctr_0_0"
    np.testing.assert_allclose(rod.minimum, (-0.05, -0.05, 0.75))
    np.testing.assert_allclose(rod.maximum, (0.05, 0.05, 2.25))


def test_unreached_rods_are_dropped():
    """H and K limits may ask for rods the scan never touches."""
    coverage = _grid_cloud((-0.2, 0.2), (-0.2, 0.2), (0.5, 2.0))
    grids = derive_ctr_grids(
        _config(),
        _Scan(),
        step=(0.01, 0.01, 0.01),
        half_width=(0.05, 0.05, 0.0),
        h_limits=3,
        k_limits=3,
        coverage=coverage,
    )
    assert [grid.name for grid in grids] == ["ctr_0_0"]


def test_symmetric_limits_expand_to_a_signed_range():
    """A single H limit means -H..H, the default symmetric selection."""
    coverage = _grid_cloud((-2.5, 2.5), (-0.2, 0.2), (0.5, 2.0))
    grids = derive_ctr_grids(
        _config(),
        _Scan(),
        step=(0.01, 0.01, 0.01),
        half_width=(0.3, 0.3, 0.0),
        h_limits=1,
        k_limits=(0, 0),
        coverage=coverage,
    )
    assert [grid.name for grid in grids] == ["ctr_m1_0", "ctr_0_0", "ctr_1_0"]


def test_ctr_grids_in_the_crystal_frame_carry_inverse_angstrom_bounds():
    """The rod box is the transformed box, in Angstrom^-1."""
    coverage = _grid_cloud((0.8, 1.2), (-0.2, 0.2), (0.0, 2.0))
    (rod,) = derive_ctr_grids(
        _config(),
        _Scan(),
        step=(0.01, 0.01, 0.01),
        half_width=(0.05, 0.05, 0.0),
        frame="crystal",
        coverage=coverage,
    )
    assert rod.name == "ctr_1_0"
    assert rod.frame == "crystal"
    # B = diag(2, 2, 0.5): the rod center moves to q = 2 Angstrom^-1, its
    # H/K half-widths stay the given Angstrom^-1 values, and the measured
    # L range 0..2 r.l.u. becomes 0..1 Angstrom^-1 along the third axis.
    np.testing.assert_allclose(rod.minimum, (1.95, -0.05, 0.0))
    np.testing.assert_allclose(rod.maximum, (2.05, 0.05, 1.0))


def test_bragg_grids_are_boxes_centered_on_each_reflection():
    """Every reachable integer reflection becomes its own symmetric box."""
    coverage = _grid_cloud((-0.2, 1.2), (-0.2, 0.2), (-0.2, 1.2))
    grids = derive_bragg_grids(
        _config(),
        _Scan(),
        step=(0.02, 0.02, 0.02),
        half_width=0.15,
        coverage=coverage,
    )
    assert [grid.name for grid in grids] == [
        "bragg_0_0_0",
        "bragg_0_0_1",
        "bragg_1_0_0",
        "bragg_1_0_1",
    ]
    peak = next(grid for grid in grids if grid.name == "bragg_1_0_1")
    np.testing.assert_allclose(peak.minimum, (0.85, -0.15, 0.85))
    np.testing.assert_allclose(peak.maximum, (1.15, 0.15, 1.15))


def test_bragg_selection_honours_explicit_index_limits():
    """Index limits narrow the enumeration before the coverage test."""
    coverage = _grid_cloud((-1.2, 1.2), (-1.2, 1.2), (-1.2, 1.2))
    grids = derive_bragg_grids(
        _config(),
        _Scan(),
        step=(0.02, 0.02, 0.02),
        half_width=0.1,
        h_limits=(1, 1),
        k_limits=0,
        l_limits=(-1, 0),
        coverage=coverage,
    )
    assert [grid.name for grid in grids] == ["bragg_1_0_m1", "bragg_1_0_0"]


def test_selection_beyond_the_grid_limit_is_refused():
    """A runaway selection fails loudly instead of queuing thousands of grids."""
    coverage = _grid_cloud((-2.5, 2.5), (-2.5, 2.5), (0.5, 2.0))
    with pytest.raises(ValueError, match="above the limit"):
        derive_ctr_grids(
            _config(),
            _Scan(),
            step=(0.01, 0.01, 0.01),
            half_width=(0.3, 0.3, 0.0),
            coverage=coverage,
            max_grids=4,
        )


def test_selected_grids_form_one_parallel_reconstruction_spec():
    """One selection is one job: unique names, one frame, one set of units."""
    coverage = _grid_cloud((-1.2, 1.2), (-1.2, 1.2), (0.2, 3.4))
    grids = derive_ctr_grids(
        _config(),
        _Scan(),
        step=(0.02, 0.02, 0.01),
        half_width=(0.1, 0.1, 0.0),
        coverage=coverage,
    )
    spec = _ReconstructionSpec(grids=tuple(grid.to_spec() for grid in grids))
    assert len(spec.grids) == len(grids)
    assert {grid.frame for grid in spec.grids} == {"hkl"}
    assert len({grid.grid_name for grid in spec.grids}) == len(grids)


def test_zero_half_width_on_every_axis_is_rejected():
    """A box needs extent; a bare reflection position is not a grid."""
    coverage = _grid_cloud((-0.2, 0.2), (-0.2, 0.2), (-0.2, 0.2))
    with pytest.raises(ValueError, match="no extent"):
        derive_bragg_grids(
            _config(),
            _Scan(),
            step=(0.01, 0.01, 0.01),
            half_width=0.0,
            coverage=coverage,
        )

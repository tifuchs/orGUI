"""Tests for reconstruction exposure-angle preparation."""

import numpy as np
import pytest

from types import SimpleNamespace

from orgui.backend import backends
from orgui.backend.interlacedScanLoader import InterlacedScan
from orgui.backend.scans import (
    ScanReference,
    SimulationScan,
    load_scan_backend_file,
    scan_exposure_angle_bounds,
)


def test_stationary_exposure_bounds_broadcast_fixed_angles():
    """Stationary exposures have identical start and end positions."""
    scan = SimulationScan((1, 1), 0.0, 10.0, 2, axis="mu", fixed=3.0)
    bounds = scan_exposure_angle_bounds(
        scan,
        SimpleNamespace(mu=0.0, chi=0.4, phi=0.5),
    )

    assert bounds.shape == (2, 2, 6)
    np.testing.assert_allclose(bounds[:, 0], bounds[:, 1])
    np.testing.assert_allclose(
        bounds[:, 0],
        [
            [0.0, np.deg2rad(-3.0), 0.4, 0.5, 0.0, 0.0],
            [np.deg2rad(10.0), np.deg2rad(-3.0), 0.4, 0.5, 0.0, 0.0],
        ],
    )


def test_swept_exposure_bounds_use_adjacent_midpoints():
    """Continuous sweeps use midpoint edges and extrapolated end edges."""
    scan = SimulationScan((1, 1), 0.0, -20.0, 3, axis="th", fixed=0.1)
    bounds = scan_exposure_angle_bounds(
        scan,
        SimpleNamespace(mu=0.1, chi=0.3, phi=0.4),
        fallback="midpoint",
    )

    np.testing.assert_allclose(
        bounds[:, :, 1],
        np.deg2rad([[-5.0, 5.0], [5.0, 15.0], [15.0, 25.0]]),
    )
    np.testing.assert_allclose(bounds[:, :, 0], np.deg2rad(0.1))


def _arm_scan(gamma_deg, delta_deg, points=3):
    scan = SimulationScan((1, 1), 0.0, 2.0, points, axis="th", fixed=0.0)
    scan.gamma_arm = np.asarray(gamma_deg, dtype=float)
    scan.delta_arm = np.asarray(delta_deg, dtype=float)
    return scan


def test_moving_arm_is_stationary_within_each_exposure_by_default():
    """Continuity is declared, never inferred from the motor moving."""
    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = _arm_scan([1.0, 2.0, 3.0], [10.0, 20.0, 30.0])

    bounds = scan_exposure_angle_bounds(scan, config)

    assert bounds.shape == (3, 2, 6)
    np.testing.assert_allclose(bounds[:, 0, 4:], bounds[:, 1, 4:])
    np.testing.assert_allclose(
        bounds[:, 0, 4], np.deg2rad([1.0, 2.0, 3.0])
    )
    np.testing.assert_allclose(
        bounds[:, 0, 5], np.deg2rad([10.0, 20.0, 30.0])
    )


def test_continuous_exposure_sweeps_the_arm_across_each_frame():
    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = _arm_scan([1.0, 2.0, 3.0], [10.0, 20.0, 30.0])
    scan.continuous_exposure = True

    bounds = scan_exposure_angle_bounds(scan, config)

    np.testing.assert_allclose(
        bounds[:, :, 4], np.deg2rad([[0.5, 1.5], [1.5, 2.5], [2.5, 3.5]])
    )
    # the sample circles keep the caller's own fallback, here stationary
    np.testing.assert_allclose(bounds[:, 0, :4], bounds[:, 1, :4])


def test_surface_frame_arm_readouts_are_converted_with_the_scan_alpha():
    """A beamline reporting six-circle arm angles declares it once.

    The stored values are the true scattering angles either way, so everything
    downstream keeps working in one convention.
    """
    from orgui.datautils.xrayutils import DetectorCalibration

    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    alpha_deg = 5.0
    scan = SimulationScan((1, 1), alpha_deg, alpha_deg, 2, axis="mu", fixed=0.0)
    scan.gamma_arm = np.array([12.0, 12.0])
    scan.delta_arm = np.array([30.0, 30.0])
    scan.arm_angle_frame = "surface"

    bounds = scan_exposure_angle_bounds(scan, config)

    expected = DetectorCalibration.armAnglesFromSurface(
        np.deg2rad(12.0), np.deg2rad(30.0), np.deg2rad(alpha_deg)
    )
    np.testing.assert_allclose(bounds[:, :, 4], expected[0])
    np.testing.assert_allclose(bounds[:, :, 5], expected[1])
    # and it really is a conversion, not a pass-through
    assert not np.isclose(bounds[0, 0, 4], np.deg2rad(12.0))


def test_scan_arm_angles_agrees_with_the_exposure_bounds():
    """One resolver, so the cursor readout and the mapper cannot disagree.

    ``scan_arm_angles`` is what the GUI asks for the arm of a single frame;
    it has to give the same answer as the bounds array built for mapping.
    """
    from orgui.backend.scans import scan_arm_angles

    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = _arm_scan([1.0, 2.0, 3.0], [10.0, 20.0, 30.0])

    gamma_arm, delta_arm = scan_arm_angles(scan, config)
    bounds = scan_exposure_angle_bounds(scan, config)

    np.testing.assert_allclose(gamma_arm, bounds[:, 0, 4])
    np.testing.assert_allclose(delta_arm, bounds[:, 0, 5])
    np.testing.assert_allclose(gamma_arm, np.deg2rad([1.0, 2.0, 3.0]))


def test_scan_arm_angles_falls_back_to_zero_without_an_arm():
    """A scan and config that know nothing about an arm give a parked one."""
    from orgui.backend.scans import scan_arm_angles

    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = SimulationScan((1, 1), 0.0, 2.0, 3, axis="th", fixed=0.0)

    gamma_arm, delta_arm = scan_arm_angles(scan, config)

    np.testing.assert_array_equal(gamma_arm, np.zeros(3))
    np.testing.assert_array_equal(delta_arm, np.zeros(3))


def test_unknown_arm_angle_frame_is_rejected():
    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = SimulationScan((1, 1), 0.0, 1.0, 2, axis="th", fixed=0.0)
    scan.arm_angle_frame = "detector"

    with pytest.raises(ValueError, match="arm_angle_frame"):
        scan_exposure_angle_bounds(scan, config)


def test_explicit_arm_bounds_win_over_every_fallback():
    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = _arm_scan([1.0, 2.0, 3.0], [10.0, 20.0, 30.0], points=2)
    scan.gamma_arm = np.array([1.0, 2.0])
    scan.delta_arm = np.array([10.0, 20.0])
    scan.continuous_exposure = True
    scan.gamma_arm_bounds_rad = np.array([[0.1, 0.2], [0.2, 0.3]])

    bounds = scan_exposure_angle_bounds(scan, config, fallback="midpoint")

    np.testing.assert_allclose(bounds[:, :, 4], scan.gamma_arm_bounds_rad)


def test_config_supplies_the_arm_when_the_backend_does_not():
    config = SimpleNamespace(
        mu=0.0, chi=0.0, phi=0.0, gamma_arm=0.25, delta_arm=-0.5
    )
    scan = SimulationScan((1, 1), 0.0, 2.0, 2, axis="th", fixed=0.0)

    bounds = scan_exposure_angle_bounds(scan, config)

    np.testing.assert_allclose(bounds[:, :, 4], 0.25)
    np.testing.assert_allclose(bounds[:, :, 5], -0.5)


def test_four_column_exposure_bounds_from_an_older_backend_still_work():
    """A backend written before the arm keeps working unchanged."""
    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    scan = SimulationScan((1, 1), 0.0, 2.0, 2, axis="th", fixed=0.0)
    scan.exposure_bounds_rad = np.arange(2 * 2 * 4, dtype=float).reshape(2, 2, 4)

    bounds = scan_exposure_angle_bounds(scan, config)

    assert bounds.shape == (2, 2, 6)
    np.testing.assert_allclose(bounds[:, :, :4], scan.exposure_bounds_rad)
    np.testing.assert_allclose(bounds[:, :, 4:], 0.0)


def test_reconstruction_refuses_a_moving_arm_rather_than_ignoring_it():
    """Dropping the arm columns would map a wrong volume that looks right."""
    from orgui.datautils.xrayutils.reconstruction import _sample_angle_bounds

    parked = np.zeros((3, 2, 6))
    parked[..., 4] = 0.3
    np.testing.assert_allclose(
        _sample_angle_bounds(parked), np.zeros((3, 2, 4))
    )

    moving = parked.copy()
    moving[2, :, 4] = 0.4
    with pytest.raises(NotImplementedError, match="moving"):
        _sample_angle_bounds(moving)


def test_exact_bounds_override_fallback_and_survive_interlacing():
    config = SimpleNamespace(mu=0.0, chi=0.0, phi=0.0)
    first = SimulationScan((1, 1), 0.0, 1.0, 2)
    second = SimulationScan((1, 1), 2.0, 3.0, 2)
    first.omega_bounds_rad = np.array([[0.0, 0.01], [0.01, 0.02]])
    second.omega_bounds_rad = np.array([[0.02, 0.03], [0.03, 0.04]])
    scan = InterlacedScan([first, second], False, "th")

    bounds = scan_exposure_angle_bounds(scan, config, fallback="midpoint")

    np.testing.assert_allclose(
        bounds[:, :, 1],
        np.concatenate(
            [first.omega_bounds_rad, second.omega_bounds_rad]
        ),
    )


def test_scan_reference_reopens_simulated_and_interlaced_scans():
    first = SimulationScan((2, 3), 0.0, 1.0, 2)
    second = SimulationScan((2, 3), 2.0, 3.0, 2)
    scan = InterlacedScan([first, second], False, "th")

    reopened = ScanReference.from_dict(
        ScanReference.from_scan(scan).to_dict()
    ).open()

    assert isinstance(reopened, InterlacedScan)
    assert len(reopened) == 4
    np.testing.assert_allclose(reopened.axis, scan.axis)


def test_scan_reference_reopens_run_path_backend(tmp_path):
    backend_file = tmp_path / "custom_backend.py"
    source_file = tmp_path / "data.h5"
    source_file.touch()
    backend_file.write_text(
        "\n".join(
            [
                "import numpy as np",
                "from orgui.backend.scans import Scan, h5_Image",
                "",
                "class CustomScan(Scan):",
                "    def __init__(self, source, number=None):",
                "        self.hdffilepath_orNode = source",
                "        self.scanno = number",
                "        self.axisname = 'th'",
                "        self.axis = np.array([0.0])",
                "        self.th = self.axis",
                "        self.omega = -self.axis",
                "        self.mu = 0.0",
                "    def __len__(self):",
                "        return 1",
                "    def get_raw_img(self, index):",
                "        return h5_Image(np.ones((1, 1)))",
                "    @classmethod",
                "    def parse_h5_node(cls, node):",
                "        return {}",
            ]
        ),
        encoding="utf-8",
    )
    _, scan_class = load_scan_backend_file(backend_file)
    scan = scan_class(str(source_file), 7)

    reference = ScanReference.from_scan(scan)
    reopened = ScanReference.from_dict(reference.to_dict()).open()

    assert reference.kind == "backend_file"
    assert reference.parameters["backend_file"] == str(backend_file.absolute())
    assert type(reopened).__qualname__ == "CustomScan"
    assert reopened.scanno == 7


def test_old_run_path_reference_uses_loaded_backend(tmp_path, monkeypatch):
    backend_file = tmp_path / "legacy_backend.py"
    source_file = tmp_path / "data.h5"
    source_file.touch()
    backend_file.write_text(
        "\n".join(
            [
                "import numpy as np",
                "from orgui.backend.scans import Scan, h5_Image",
                "",
                "class LegacyScan(Scan):",
                "    def __init__(self, source):",
                "        self.hdffilepath_orNode = source",
                "        self.axisname = 'th'",
                "        self.axis = np.array([0.0])",
                "        self.th = self.axis",
                "        self.omega = -self.axis",
                "        self.mu = 0.0",
                "    def __len__(self):",
                "        return 1",
                "    def get_raw_img(self, index):",
                "        return h5_Image(np.ones((1, 1)))",
                "    @classmethod",
                "    def parse_h5_node(cls, node):",
                "        return {}",
            ]
        ),
        encoding="utf-8",
    )
    _, scan_class = load_scan_backend_file(backend_file)
    monkeypatch.setitem(backends.fscans, "legacy-test", scan_class)
    reference = ScanReference(
        kind="backend",
        module="<run_path>",
        class_name="LegacyScan",
        parameters={"source": str(source_file), "scan_number": None},
        source_fingerprints=(ScanReference._fingerprint(source_file),),
    )

    assert type(reference.open()).__qualname__ == "LegacyScan"

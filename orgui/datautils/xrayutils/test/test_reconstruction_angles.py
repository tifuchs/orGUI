"""Tests for reconstruction exposure-angle preparation."""

import numpy as np

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

    assert bounds.shape == (2, 2, 4)
    np.testing.assert_allclose(bounds[:, 0], bounds[:, 1])
    np.testing.assert_allclose(
        bounds[:, 0],
        [
            [0.0, np.deg2rad(-3.0), 0.4, 0.5],
            [np.deg2rad(10.0), np.deg2rad(-3.0), 0.4, 0.5],
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

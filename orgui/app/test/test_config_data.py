import h5py
import numpy as np
from silx.io.dictdump import dicttonx, nxtodict
import pytest

from orgui.app.QReflectionSelector import HKLReflection
from orgui.app.config_data import ConfigData, ConfigHandler
from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg


def _make_config():
    unit_cell = CTRcalc.UnitCell([3.0, 4.0, 5.0], [90.0, 91.0, 120.0], name="bulk")
    unit_cell.addAtom("Pt", [0.0, 0.5, 0.25], 0.1, 0.2, 0.8, 1)
    unit_cell.addAtom("O", [0.25, 0.25, 0.5], 0.3, 0.4, 0.9, 2)
    unit_cell.coherentDomainMatrix = [
        np.vstack((np.eye(3).T, np.array([0.0, 0.0, 0.0]))).T,
        np.array(
            [
                [1.0, 0.0, 0.0, 0.5],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.25],
            ]
        ),
    ]
    unit_cell.coherentDomainOccupancy = [0.7, 0.3]

    ub_calculator = HKLVlieg.UBCalculator(unit_cell, 70.0)
    ub_calculator.defaultU()

    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(729.0, 731.0, 1587.0, pixelX=172.0, pixelY=172.0)
    detector.set_wavelength(ub_calculator.getLambda() * 1e-10)
    detector.detector.shape = (120, 240)
    detector.detector.max_shape = (120, 240)
    detector.detector.binning = (2, 3)
    detector.set_roi([10, 90], [20, 140])
    detector.setAzimuthalReference(0.2)
    detector.setPolarization(0.3, 0.9)

    reflections = [
        HKLReflection([11.0, 12.0], [1.0, 0.0, 0.0], 5, "ref_a"),
        HKLReflection([21.0, 22.0], [0.0, 1.0, 0.0], 6, "ref_b"),
    ]
    return ConfigData(
        detector,
        unit_cell,
        ub_calculator,
        mu=0.1,
        chi=0.2,
        phi=0.3,
        refraction_index=0.999,
        reference_reflections=reflections,
    )


def test_config_data_round_trips_through_nexus_dict(tmp_path):
    config = _make_config()
    filename = tmp_path / "config.h5"

    dicttonx({"configuration": config.to_nxdict(role="scan")}, filename)
    loaded_dict = nxtodict(filename)["configuration"]
    loaded = ConfigData.from_nxdict(loaded_dict)

    assert loaded_dict["@orgui_meta"] == "config"
    assert loaded_dict["@orgui_config_role"] == "scan"
    assert loaded_dict["@orgui_schema_version"] == 1
    assert np.allclose(
        loaded.detector.get_config()["dist"], config.detector.get_config()["dist"]
    )
    assert np.allclose(loaded.detector.wavelength, config.detector.wavelength)
    assert np.allclose(
        loaded.detector.getPolarization(), config.detector.getPolarization()
    )
    assert np.allclose(
        loaded.detector.getAzimuthalReference(),
        config.detector.getAzimuthalReference(),
    )
    assert tuple(loaded.detector.detector._binning) == tuple(
        config.detector.detector._binning
    )
    assert tuple(loaded.detector.detector.shape) == tuple(
        config.detector.detector.shape
    )
    assert np.allclose(loaded.detector._roi, config.detector._roi)

    assert loaded.unit_cell.name == config.unit_cell.name
    assert np.allclose(loaded.unit_cell.a, config.unit_cell.a)
    assert np.allclose(loaded.unit_cell.alpha, config.unit_cell.alpha)
    assert loaded.unit_cell.names == config.unit_cell.names
    assert np.allclose(loaded.unit_cell.basis, config.unit_cell.basis)
    assert np.allclose(
        loaded.unit_cell.coherentDomainMatrix,
        config.unit_cell.coherentDomainMatrix,
    )
    assert np.allclose(
        loaded.unit_cell.coherentDomainOccupancy,
        config.unit_cell.coherentDomainOccupancy,
    )

    assert np.allclose(loaded.ub_calculator.getU(), config.ub_calculator.getU())
    assert np.allclose(loaded.ub_calculator.getUB(), config.ub_calculator.getUB())
    assert np.isclose(
        loaded.ub_calculator.getEnergy(), config.ub_calculator.getEnergy()
    )
    assert np.isclose(
        loaded.ub_calculator.getLambda(), config.ub_calculator.getLambda()
    )
    assert [refl.identifier for refl in loaded.reference_reflections] == [
        "ref_a",
        "ref_b",
    ]
    assert np.allclose(loaded.reference_reflections[0].hkl, [1.0, 0.0, 0.0])
    assert np.allclose(loaded.reference_reflections[1].xy, [21.0, 22.0])


def test_config_handler_writes_scan_and_integration_config(tmp_path):
    config = _make_config()
    filename = tmp_path / "database.h5"

    with h5py.File(filename, "w") as h5file:
        scan = h5file.create_group("scan_1")
        integration = h5file.create_group("scan_1/measurement/result")
        handler = ConfigHandler()

        scan_config = handler.write_scan_config(scan, config)
        integration_config = handler.write_integration_config(integration, config)

        assert ConfigHandler.is_config_group(scan_config)
        assert ConfigHandler.is_config_group(integration_config)
        assert scan_config.attrs["orgui_config_role"] == "scan"
        assert integration_config.attrs["orgui_config_role"] == "integration"
        assert not ConfigHandler.is_config_group(scan)
        loaded = handler.load_config_group(scan_config)
        assert np.allclose(loaded.ub_calculator.getU(), config.ub_calculator.getU())


def test_config_data_writes_with_blosc_compression(tmp_path):
    hdf5plugin = pytest.importorskip("hdf5plugin")
    config = _make_config()
    filename = tmp_path / "compressed_config.h5"
    compression = hdf5plugin.Blosc(
        cname="lz4", shuffle=hdf5plugin.Blosc.SHUFFLE, clevel=5
    )

    dicttonx(
        {"configuration": config.to_nxdict(role="scan")},
        filename,
        create_dataset_args=compression,
    )
    loaded = ConfigData.from_nxdict(nxtodict(filename)["configuration"])

    assert loaded.unit_cell.names == config.unit_cell.names
    assert [refl.identifier for refl in loaded.reference_reflections] == [
        "ref_a",
        "ref_b",
    ]

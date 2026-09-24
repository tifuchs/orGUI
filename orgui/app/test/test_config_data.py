import h5py
import json

import numpy as np
from silx.gui import qt
from silx.io.dictdump import dicttonx, nxtodict
import pytest
from types import SimpleNamespace

from orgui.app.QReflectionSelector import HKLReflection
from orgui.app.config_data import CorrectionState, ConfigData, ConfigHandler
from orgui.app.database import config_data_from_json, config_data_to_json
from orgui.app.peak1Dintegr import IntegrationCorrectionsDialog
from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg
from orgui.reconstruction_job import _snapshot_assets


@pytest.fixture(scope="session")
def qapp():
    """The Qt application, kept referenced for the whole test session.

    A ``QApplication`` that is not held on to is garbage-collected, and
    creating a widget afterwards aborts the interpreter.
    """
    application = qt.QApplication.instance()
    if application is None:
        application = qt.QApplication([])
    return application


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
        corrections=CorrectionState(
            use_mask=True,
            use_background=True,
            use_solid_angle=True,
            normalize_exposure=False,
            monitor_corrections=("mondio",),
            excluded_frames=(2, 7),
            uncertainty_provenance={"background": "measured"},
        ),
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
    assert loaded.corrections == config.corrections


def test_config_data_round_trips_through_database_json():
    config = _make_config()
    loaded = config_data_from_json(config_data_to_json(config))

    assert loaded.corrections == config.corrections
    assert np.allclose(loaded.ub_calculator.getUB(), config.ub_calculator.getUB())


@pytest.fixture
def legacy_corrections_database(tmp_path):
    """Database snapshot using the pre-typed correction-settings group."""
    config = _make_config()
    values = {
        "use_mask": True,
        "use_background": True,
        "use_solid_angle": True,
        "use_polarization": True,
        "use_lorentz": True,
        "use_footprint": False,
        "use_normalization": True,
        "normalize_exposure": True,
        "monitor_corrections": ["mondio"],
        "sample_length_m": 0.004,
        "sample_width_m": 0.008,
        "beam_flux_density": 2.5e12,
    }
    nxdict = config.to_nxdict(role="scan")
    nxdict["orgui"]["integration_corrections"] = {
        "@NX_class": "NXcollection",
        "json": json.dumps(values),
    }
    filename = tmp_path / "legacy_corrections.h5"
    dicttonx({"configuration": nxdict}, filename)
    return filename, values


def test_pre_typed_database_correction_settings_still_load(
    legacy_corrections_database,
):
    """Pin dispatch of an old database containing one opaque JSON dataset."""
    filename, values = legacy_corrections_database

    stored = nxtodict(filename)["configuration"]
    loaded = ConfigData.from_nxdict(stored)

    assert "json" in stored["orgui"]["integration_corrections"]
    serialized = loaded.corrections.to_dict()
    assert {name: serialized[name] for name in values} == values
    assert loaded.corrections.monitor_corrections == ("mondio",)


def test_enabled_pixel_repair_implies_mask_correction():
    config = _make_config()
    repair = SimpleNamespace(
        enabled=True,
        max_component_pixels=4,
        max_span=3,
        radius=2,
        min_valid_neighbors=6,
        use_pyfai_gaps=True,
        gap_size_px=1,
    )
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=config.detector,
            crystal=config.unit_cell,
            ubCal=config.ub_calculator,
            mu=config.mu,
            chi=config.chi,
            phi=config.phi,
            n=config.refraction_index,
        ),
        scanSelector=SimpleNamespace(
            get_integration_options=lambda: {
                "mask": False,
                "solid_angle": False,
                "polarization": False,
            }
        ),
        maskManager=SimpleNamespace(
            settings=SimpleNamespace(pixel_repair=repair)
        ),
        excludedImagesDialog=SimpleNamespace(
            getData=lambda: np.empty(0, dtype=np.int64)
        ),
        reconstruction_normalize_exposure=False,
        reconstruction_monitor_corrections=("mondio",),
    )

    captured = ConfigData.from_gui(gui)

    assert captured.corrections.repair_masked_pixels is True
    assert captured.corrections.use_mask is True
    assert captured.corrections.normalize_exposure is False
    assert captured.corrections.monitor_corrections == ("mondio",)
    assert captured.corrections.shared_frame_normalization is True

    gui.ctr_correction_state = CorrectionState(
        normalize_exposure=True, monitor_corrections=("current",)
    )
    captured = ConfigData.from_gui(gui)
    assert captured.corrections.normalize_exposure is True
    assert captured.corrections.monitor_corrections == ("current",)


def test_from_gui_captures_the_footprint_dialogs_inputs(qapp):
    """L, W and the beam flux are captured only if the dialog was opened.

    An unopened footprint dialog has nothing to record, which must be the
    same "not recorded" state as a config written before these fields
    existed -- not a silent zero.
    """
    config = _make_config()
    unopened_gui = SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=config.detector,
            crystal=config.unit_cell,
            ubCal=config.ub_calculator,
            mu=config.mu,
            chi=config.chi,
            phi=config.phi,
            n=config.refraction_index,
        ),
        scanSelector=SimpleNamespace(
            get_integration_options=lambda: {
                "mask": False,
                "solid_angle": False,
                "polarization": False,
            },
            correctionsDialog=SimpleNamespace(footprintOptions=None),
        ),
        excludedImagesDialog=SimpleNamespace(
            getData=lambda: np.empty(0, dtype=np.int64)
        ),
    )
    unopened = ConfigData.from_gui(unopened_gui)
    assert unopened.corrections.sample_length_m is None
    assert unopened.corrections.sample_width_m is None
    assert unopened.corrections.beam_flux_density is None
    assert unopened.corrections.beam_shape_name is None
    assert unopened.corrections.beam_shape_values == ()

    footprint_dialog = IntegrationCorrectionsDialog()
    try:
        footprint_dialog.L.setValue(2.5)
        footprint_dialog.W.setValue(7.5)
        footprint_dialog.beamFlux.setValue(4.2)
        footprint_dialog.shapeSelector.setCurrentIndex(
            footprint_dialog.shapeSelector.findText("Trapezoid")
        )
        base, flat = footprint_dialog.shapeParameters[:2]
        base.setValue(90.0)
        flat.setValue(30.0)
        footprint_dialog.profileCenter.setCurrentIndex(
            footprint_dialog.profileCenter.findText("median")
        )
        footprint_dialog.profileOffset.setValue(-12.5)
        gui = SimpleNamespace(**unopened_gui.__dict__)
        gui.scanSelector = SimpleNamespace(
            get_integration_options=unopened_gui.scanSelector.get_integration_options,
            correctionsDialog=SimpleNamespace(footprintOptions=footprint_dialog),
        )

        captured = ConfigData.from_gui(gui)

        assert captured.corrections.sample_length_m == pytest.approx(2.5e-3)
        assert captured.corrections.sample_width_m == pytest.approx(7.5e-3)
        assert captured.corrections.beam_flux_density == pytest.approx(4.2e6)
        assert captured.corrections.beam_shape_analytical is True
        assert captured.corrections.beam_shape_name == "Trapezoid"
        assert captured.corrections.beam_shape_values == (90.0, 30.0)
        assert captured.corrections.beam_profile_center == "median"
        assert captured.corrections.beam_profile_offset_um == pytest.approx(-12.5)
    finally:
        footprint_dialog.deleteLater()


def test_from_gui_embeds_a_loaded_measured_beam_profile():
    """A future reader must not depend on the original profile file path."""
    config = _make_config()
    positions = np.array([-2e-4, 0.0, 2e-4])
    density = np.array([500.0, 4000.0, 500.0])
    profile = SimpleNamespace(profile_curve=lambda: (positions, density))
    footprint = SimpleNamespace(
        sampleLength=lambda: 3e-3,
        sampleWidth=lambda: 8e-3,
        beamFluxDensity=lambda: 2e12,
        settings=lambda: {
            "analytical": False,
            "profile_file": "moved/or/replaced.dat",
            "profile_content": "intensity",
            "profile_unit": "mm",
        },
        measuredProfile=lambda: profile,
    )
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=config.detector,
            crystal=config.unit_cell,
            ubCal=config.ub_calculator,
            mu=config.mu,
            chi=config.chi,
            phi=config.phi,
            n=config.refraction_index,
        ),
        scanSelector=SimpleNamespace(
            get_integration_options=lambda: {},
            correctionsDialog=SimpleNamespace(footprintOptions=footprint),
        ),
    )

    captured = ConfigData.from_gui(gui).corrections

    assert captured.beam_profile_positions_m == pytest.approx(positions)
    assert captured.beam_profile_density_per_m == pytest.approx(density)


def test_apply_to_gui_restores_the_footprint_dialogs_inputs(qapp):
    """Restoring a config sets L, W and the beam flux on the shared dialog."""
    config = _make_config()
    config.corrections = CorrectionState(
        sample_length_m=4e-3,
        sample_width_m=9e-3,
        beam_flux_density=3e12,
        beam_shape_analytical=True,
        beam_shape_name="Trapezoid",
        beam_shape_values=(90.0, 30.0),
        beam_profile_center="median",
        beam_profile_offset_um=-12.5,
    )
    footprint_dialog = IntegrationCorrectionsDialog()
    try:
        corrections_dialog = SimpleNamespace(
            footprintOptions_shared=lambda: footprint_dialog
        )
        gui = SimpleNamespace(
            ubcalc=SimpleNamespace(
                detectorCal=config.detector,
                crystal=config.unit_cell,
                ubCal=config.ub_calculator,
                mu=config.mu,
                chi=config.chi,
                phi=config.phi,
                n=config.refraction_index,
            ),
            scanSelector=SimpleNamespace(
                get_integration_options=lambda: {},
                set_integration_options=lambda options: None,
                correctionsDialog=corrections_dialog,
            ),
        )

        config.apply_to_gui(gui)

        assert footprint_dialog.sampleLength() == pytest.approx(4e-3)
        assert footprint_dialog.sampleWidth() == pytest.approx(9e-3)
        assert footprint_dialog.beamFluxDensity() == pytest.approx(3e12)
        assert footprint_dialog.analyticalButton.isChecked() is True
        assert footprint_dialog.currentShape().name == "Trapezoid"
        assert footprint_dialog.shapeParameters[0].value() == pytest.approx(90.0)
        assert footprint_dialog.shapeParameters[1].value() == pytest.approx(30.0)
        assert footprint_dialog.profileCenter.currentText() == "median"
        assert footprint_dialog.profileOffset.value() == pytest.approx(-12.5)
    finally:
        footprint_dialog.deleteLater()


def test_apply_to_gui_sets_reconstruction_normalization_attributes():
    config = _make_config()
    config.corrections = CorrectionState(
        normalize_exposure=False, monitor_corrections=("mondio",)
    )
    refreshed = []
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=config.detector,
            crystal=config.unit_cell,
            ubCal=config.ub_calculator,
            mu=config.mu,
            chi=config.chi,
            phi=config.phi,
            n=config.refraction_index,
        ),
        scanSelector=SimpleNamespace(
            set_integration_options=lambda options: None,
            correctionsDialog=SimpleNamespace(
                refresh=lambda: refreshed.append(True)
            ),
        ),
    )

    config.apply_to_gui(gui)

    assert gui.reconstruction_normalize_exposure is False
    assert gui.reconstruction_monitor_corrections == ("mondio",)
    assert gui.ctr_correction_state.normalize_exposure is False
    assert gui.ctr_correction_state.monitor_corrections == ("mondio",)
    assert refreshed == [True]


def test_apply_to_gui_peak_angles_use_the_loaded_orientation():
    """Loaded U must reach ``angles``, which calcReflection reads."""
    config = _make_config()
    rotation = np.deg2rad(68.0)
    config.ub_calculator.setU(
        np.array(
            [
                [np.cos(rotation), -np.sin(rotation), 0.0],
                [np.sin(rotation), np.cos(rotation), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
    )
    stale = HKLVlieg.UBCalculator(config.unit_cell, 70.0)
    stale.defaultU()
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=config.detector,
            crystal=config.unit_cell,
            ubCal=stale,
            angles=HKLVlieg.VliegAngles(stale),
        ),
    )

    config.apply_to_gui(gui)

    hkl = np.array([[1.0], [0.0], [1.0]])
    expected = HKLVlieg.VliegAngles(config.ub_calculator).anglesZmode(
        hkl, config.mu, "in", config.chi, config.phi
    )
    np.testing.assert_allclose(
        gui.ubcalc.angles.anglesZmode(hkl, config.mu, "in", config.chi, config.phi),
        expected,
    )


def test_total_flux_state_survives_before_its_widgets_exist():
    """Stage-5 programmatic settings round-trip through the GUI snapshot."""
    config = _make_config()
    config.corrections = CorrectionState(
        total_incident_flux=2.5e10,
        total_flux_calibrated=True,
        primary_monitor="ic2",
        primary_monitor_kind="rate",
        primary_monitor_unit="count/s",
        monitor_reference_reading=4200.0,
        horizontal_interception="fraction",
        horizontal_intercepted_fraction=0.85,
    )
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(
            detectorCal=config.detector,
            crystal=config.unit_cell,
            ubCal=config.ub_calculator,
            mu=config.mu,
            chi=config.chi,
            phi=config.phi,
            n=config.refraction_index,
        ),
        scanSelector=SimpleNamespace(
            get_integration_options=lambda: {},
            set_integration_options=lambda options: None,
        ),
    )

    config.apply_to_gui(gui)
    captured = ConfigData.from_gui(gui).corrections

    assert captured.total_incident_flux == pytest.approx(2.5e10)
    assert captured.total_flux_calibrated is True
    assert captured.primary_monitor == "ic2"
    assert captured.primary_monitor_kind == "rate"
    assert captured.monitor_reference_reading == pytest.approx(4200.0)
    assert captured.horizontal_interception == "fraction"
    assert captured.horizontal_intercepted_fraction == pytest.approx(0.85)


def test_apply_to_gui_requests_a_replot(qapp):
    """Loading a config must refresh HKL-dependent plot elements too.

    Every interactive path that changes mu/chi/phi, the UB matrix or the
    detector geometry (``_onMachineParamsChanged``, ``_onCrystalParamsChanged``,
    ``_onAlignU``) emits ``sigReplotRequest`` right after
    ``updateReflectionMismatch()`` -- that method's own docstring documents
    the contract. Without it here, ``mu``/``chi``/``phi`` and the UB matrix
    were already correct immediately after loading (they are read fresh, not
    cached), but the ROI and reflection overlays and the Q-plot kept showing
    the geometry from before the config was loaded -- what looked like
    "loading a config does not update the angles to hkl conversion."
    """
    config = _make_config()

    class _UBStub(qt.QObject):
        sigPlottableMachineParamsChanged = qt.pyqtSignal()
        sigReplotRequest = qt.pyqtSignal(bool)

        def __init__(self):
            super().__init__()
            self.detectorCal = config.detector
            self.crystal = config.unit_cell
            self.ubCal = config.ub_calculator
            self.mu = 0.0
            self.chi = 0.0
            self.phi = 0.0
            self.n = 1.0

        def updateReflectionMismatch(self):
            pass

    ub_widget = _UBStub()
    gui = SimpleNamespace(ubcalc=ub_widget)

    replot_calls = []
    ub_widget.sigReplotRequest.connect(replot_calls.append)
    plottable_calls = []
    ub_widget.sigPlottableMachineParamsChanged.connect(
        lambda: plottable_calls.append(True)
    )

    config.apply_to_gui(gui)

    assert replot_calls == [True]
    assert plottable_calls == [True]
    assert ub_widget.mu == config.mu, "the angle itself must also be updated"


def test_snapshot_assets_serializes_active_mask(tmp_path):
    config = _make_config()
    config.corrections = CorrectionState(use_mask=True)
    mask = np.zeros(config.detector.detector.shape, dtype=bool)
    mask[3, 4] = True
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(detectorCal=config.detector),
        get_detector_mask=lambda shape: mask,
    )
    assets = tmp_path / "job-assets.nxs"

    _snapshot_assets(gui, config, assets)

    assert config.corrections.mask_asset == "/mask"
    with h5py.File(assets, "r") as h5file:
        np.testing.assert_array_equal(h5file["mask"][()], mask)


def test_snapshot_assets_rejects_missing_enabled_mask(tmp_path):
    config = _make_config()
    config.corrections = CorrectionState(use_mask=True)
    gui = SimpleNamespace(
        ubcalc=SimpleNamespace(detectorCal=config.detector),
        get_detector_mask=lambda shape: None,
    )

    with pytest.raises(ValueError, match="no active mask matches"):
        _snapshot_assets(gui, config, tmp_path / "job-assets.nxs")


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

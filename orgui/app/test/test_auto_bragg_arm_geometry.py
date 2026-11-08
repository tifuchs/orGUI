"""Automatic orientation must use the measured frame's detector geometry."""

import copy
from types import MethodType, SimpleNamespace

import numpy as np
import pytest

from orgui.app import autoBraggSearch
from orgui.app.autoBraggWorkflow import AutoBraggWorkflow
from orgui.app.orGUI import orGUI
from orgui.app.QUBCalculator import QUBCalculator
from orgui.app.QReflectionSelector import QReflectionSelector
from orgui.backend.interlacedScanLoader import InterlacedScan
from orgui.backend.scans import SimulationScan
from orgui.datautils.xrayutils import DetectorCalibration, HKLVlieg
from orgui.datautils.xrayutils import ReciprocalNavigation as rn


def _detector(home=(0.0, 0.0)):
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(400, 300, 500, tilt=3, tiltPlanRotation=20,
                      pixelX=200, pixelY=200)
    detector.detector.shape = (1100, 800)
    detector.setArmReference(gamma_arm=home[0], delta_arm=home[1])
    return detector


def _gui(scan, detector, ub=None, mu=0.01, chi=0.02, phi=-0.03):
    if ub is None:
        ub = HKLVlieg.UBCalculator(
            HKLVlieg.Lattice([4, 5, 6], [90, 90, 90]), 70
        )
        ub.setU(np.eye(3))
    calculator = SimpleNamespace(
        detectorCal=detector, ubCal=ub, crystal=ub.lattice,
        angles=HKLVlieg.VliegAngles(ub), n=1.0, mu=mu, chi=chi, phi=phi,
        gamma_arm=0.0, delta_arm=0.0, arm_angle_frame="prim",
        uedit=SimpleNamespace(setU=lambda _u: None),
        updateReflectionMismatch=lambda: None,
        sigReplotRequest=SimpleNamespace(emit=lambda _value: None),
    )
    detector.wavelength = ub.getLambda() * 1e-10
    gui = SimpleNamespace(fscan=scan, ubcalc=calculator)
    calculator.mainGui = gui
    calculator.calcReflection = MethodType(QUBCalculator.calcReflection, calculator)
    calculator.getReflectionMismatch = ub.getReflectionMismatch
    for name in ("getMuOm", "getArmAngles", "isValidImageNo", "imageNoToAxis",
                 "axisToImageNo", "getReflections", "getReflectionFitData",
                 "searchPixelCoordHKL", "armFrameGroups", "calcBraggRefl",
                 "_qRangePerFrame"):
        setattr(gui, name, MethodType(getattr(orGUI, name), gui))
    reflections = []

    def add_reflection(event, imageno, hkl):
        """Record a reference reflection in the test selector."""
        reflection = SimpleNamespace(
            xy=np.array([event["x"], event["y"]]), imageno=imageno,
            hkl=np.asarray(hkl),
        )
        reflections.append(reflection)
        return reflection

    gui.reflectionSel = SimpleNamespace(
        reflections=reflections, addReflection=add_reflection
    )
    return gui


def _peak(imageno, xy):
    candidate = autoBraggSearch.ImageMaximum(imageno, np.asarray(xy), 100.0)
    return autoBraggSearch.RefinedPeak(candidate, candidate.xy, imageno, 0.0, 100.0)


@pytest.mark.parametrize("metadata", ["missing", "missing_zero", "scalar",
                                      "array", "surface", "segmented"])
def test_peak_q_uses_scan_arm_and_full_sample_rotation(metadata):
    """Compare with rotating a calibrated outgoing ray in the laboratory."""
    home = np.deg2rad([2.0, -4.0])
    detector = _detector(home)
    scan = SimulationScan((1, 1), 0.3, 0.9, 3, axis="mu", fixed=12.5)
    if metadata == "scalar":
        scan.gamma_arm, scan.delta_arm = 3.0, 7.0
    elif metadata in {"array", "surface"}:
        scan.gamma_arm = np.array([1.0, 3.0, 6.0])
        scan.delta_arm = np.array([2.0, 7.0, 11.0])
        if metadata == "surface":
            scan.arm_angle_frame = "surface"
    elif metadata == "segmented":
        first = SimulationScan((1, 1), 0.3, 0.9, 2, axis="mu", fixed=12.5)
        second = SimulationScan((1, 1), 0.5, 0.7, 2, axis="mu", fixed=12.5)
        first.gamma_arm, first.delta_arm = [1.0, 6.0], [2.0, 11.0]
        second.gamma_arm, second.delta_arm = [3.0, 4.0], [7.0, 8.0]
        scan = InterlacedScan([first, second], True, "mu")
        np.testing.assert_array_equal(scan.gamma_arm, [1., 3., 4., 6.])
    gui = _gui(scan, detector)
    # Test the config fallback separately from omission at the calibration.
    if metadata != "missing_zero":
        gui.ubcalc.gamma_arm, gui.ubcalc.delta_arm = np.deg2rad([1.5, 5.0])
    xy = np.array([347.3, 219.7])
    gamma0, delta0 = detector.primBeamPoints(xy[1:2], xy[0:1])
    home_ray = np.array([np.sin(delta0[0]) * np.cos(gamma0[0]),
                         np.cos(delta0[0]) * np.cos(gamma0[0]),
                         np.sin(gamma0[0])])
    for index in range(len(scan)):
        peak = _peak(index, xy)
        mu, omega = gui.getMuOm(index)
        if metadata in {"missing", "missing_zero"}:
            gamma_arm, delta_arm = (
                np.deg2rad([1.5, 5.0]) if metadata == "missing" else (0., 0.)
            )
        else:
            gamma_arm = np.deg2rad(
                np.broadcast_to(scan.gamma_arm, (len(scan),))[index]
            )
            delta_arm = np.deg2rad(
                np.broadcast_to(scan.delta_arm, (len(scan),))[index]
            )
            if metadata == "surface":
                # Surface outgoing ray rotated by incidence into the beam frame.
                ray = HKLVlieg.calcALPHA(mu) @ (
                    HKLVlieg.calcDELTA(delta_arm)
                    @ HKLVlieg.calcGAMMA(gamma_arm) @ np.array([0., 1., 0.])
                )
                gamma_arm = np.arcsin(ray[2])
                delta_arm = np.arctan2(ray[0], ray[1])
        outgoing = (HKLVlieg.calcDELTA(delta_arm) @ HKLVlieg.calcGAMMA(gamma_arm)
                    @ (HKLVlieg.calcDELTA(home[1])
                       @ HKLVlieg.calcGAMMA(home[0])).T @ home_ray)
        q_lab = gui.ubcalc.ubCal.getK() * (outgoing - [0., 1., 0.])
        # Laboratory -> surface -> phi; independent of pixel-to-surface angles.
        expected = (HKLVlieg.calcPHI(gui.ubcalc.phi).T
                    @ HKLVlieg.calcCHI(gui.ubcalc.chi).T
                    @ HKLVlieg.calcOMEGA(omega).T
                    @ HKLVlieg.calcALPHA(mu).T @ q_lab)
        measured = autoBraggSearch.q_phi_from_peak(
            peak, detector, gui.getMuOm, gui.ubcalc.chi, gui.ubcalc.phi,
            gui.ubcalc.ubCal.getK(), arm_resolver=gui.getArmAngles,
        )
        np.testing.assert_allclose(measured, expected, atol=2e-12)
        calibrated_q = autoBraggSearch.q_phi_from_peak(
            peak, detector, gui.getMuOm, gui.ubcalc.chi, gui.ubcalc.phi,
            gui.ubcalc.ubCal.getK(),
        )
        assert np.linalg.norm(calibrated_q - expected) > 0.1
        np.testing.assert_allclose(
            AutoBraggWorkflow(gui)._automaticPeakAngles(peak),
            autoBraggSearch.angles_from_peak(
                peak, detector, gui.getMuOm, gui.ubcalc.chi, gui.ubcalc.phi,
                gui.getArmAngles,
            ), atol=1e-14,
        )


def test_omitted_resolver_and_fixed_home_preserve_calibration_geometry():
    """Preserve calibrated geometry for legacy callers and fixed-home scans."""
    home = np.deg2rad([2.0, -4.0])
    detector = _detector(home)
    peak = _peak(0, [340.0, 210.0])
    def resolve(_index):
        """Resolve fixed calibration-home arm angles for a test frame."""
        return 0.01, -0.2
    calibrated = autoBraggSearch.q_phi_from_peak(
        peak, detector, resolve, 0.02, -0.03, 20.0
    )
    explicit = autoBraggSearch.q_phi_from_peak(
        peak, detector, resolve, 0.02, -0.03, 20.0,
        arm_resolver=lambda _index: home,
    )
    np.testing.assert_allclose(explicit, calibrated, atol=1e-13)
    zero = autoBraggSearch.q_phi_from_peak(
        peak, detector, resolve, 0.02, -0.03, 20.0,
        arm_resolver=lambda _index: (0.0, 0.0),
    )
    assert np.linalg.norm(zero - calibrated) > 0.1


def _indexed_peaks():
    detector = _detector(np.deg2rad([2.0, -4.0]))
    scan = SimulationScan((1, 1), -90.0, 0.0, 3, fixed=0.5)
    gui = _gui(scan, detector, mu=np.deg2rad(0.5))
    hkls = np.array([[1., 1., 1.], [2., 1., 1.], [1., 2., 2.]])
    positions = gui.ubcalc.angles.anglesZmode(
        hkls.T, gui.ubcalc.mu, "in", gui.ubcalc.chi, gui.ubcalc.phi
    )
    scan.th = scan.axis = -np.rad2deg(positions[:, 3])
    scan.omega = -scan.th
    scan.gamma_arm = [1., 4., 7.]
    scan.delta_arm = [2., 7., 11.]
    peaks = []
    for index, pos in enumerate(positions):
        yx = detector.pixelsSurfaceAngles(
            pos[2], pos[1], pos[0], *gui.getArmAngles(index)
        )[0]
        peaks.append(_peak(index, yx[::-1]))
    return gui, hkls, positions, peaks


def test_seed_scoring_and_scale_estimation_follow_frame_arm(monkeypatch):
    """Index measured peaks and estimate scale using their frame arm positions."""
    gui, hkls, _positions, peaks = _indexed_peaks()
    ub = gui.ubcalc.ubCal
    norms = np.linalg.norm((ub.lattice.B_mat @ hkls.T).T, axis=1)

    def shell(_xtal, qnorm, tolerance, **_kwargs):
        """Supply the known reflection shell for seed matching."""
        keep = np.abs(norms - qnorm) <= tolerance
        return hkls[keep], norms[keep]

    # Only reflection enumeration is substituted; geometry, U solver and
    # indexing scores remain real numerical calculations without native F_uc.
    monkeypatch.setattr(autoBraggSearch, "allowed_bragg_by_qnorm", shell)
    def intensity_shell(*args, **kwargs):
        """Supply known reflection intensities for scale estimation."""
        indices, lengths = shell(*args, **kwargs)
        return indices, np.ones(len(indices)), lengths

    monkeypatch.setattr(autoBraggSearch, "allowed_bragg_in_q_region", intensity_shell)
    monkeypatch.setattr(autoBraggSearch, "allowed_bragg_same_q", intensity_shell)
    seeds = autoBraggSearch.rank_ub_seeds(
        peaks, gui.ubcalc.detectorCal, gui.getMuOm, ub, ub.lattice,
        gui.ubcalc.chi, gui.ubcalc.phi, qnorm_tolerance=1e-7,
        arm_resolver=gui.getArmAngles,
    )
    assert len(seeds) == len(peaks)
    # The score divides inverse-projection roundoff by the tight Q tolerance.
    assert seeds[0].score < 1e-3
    np.testing.assert_allclose(seeds[0].U, ub.getU(), atol=1e-10)
    # The old calibrated-arm path misses these known Q shells completely.
    assert not autoBraggSearch.rank_ub_seeds(
        peaks, gui.ubcalc.detectorCal, gui.getMuOm, ub, ub.lattice,
        gui.ubcalc.chi, gui.ubcalc.phi, qnorm_tolerance=1e-7,
    )
    position_seeds = AutoBraggWorkflow(gui)._rankAutomaticSeedsByDetectorPosition(
        peaks[1], q_tolerance=1e-7, pixel_tolerance=1e-5,
    )
    assert len(position_seeds) == 1
    np.testing.assert_array_equal(position_seeds[0].hkl, hkls[1])
    np.testing.assert_allclose(position_seeds[0].U, ub.getU(), atol=1e-10)
    scale = autoBraggSearch.estimate_qnorm_scale(
        peaks, gui.ubcalc.detectorCal, gui.getMuOm, ub, ub.lattice,
        gui.ubcalc.chi, gui.ubcalc.phi, 1e-7, arm_resolver=gui.getArmAngles,
    )
    assert scale is not None
    assert scale["qnorm_scale"] == pytest.approx(1.0, abs=1e-10)


def test_prediction_observed_insertion_and_reference_recalculation_roundtrip():
    """Preserve frame geometry through prediction, insertion and U recalculation."""
    gui, hkls, positions, peaks = _indexed_peaks()
    workflow = AutoBraggWorkflow(gui)
    # Search must reproject at the predicted frame, independently of the
    # currently displayed image or the calibration's arm reference.
    gui.imageno = 2
    reflection = gui.searchPixelCoordHKL(hkls[1])
    assert reflection["imageno_1"] == 1
    assert reflection["selectable_1"]
    np.testing.assert_allclose(reflection["xy_1"], peaks[1].xy, atol=1e-8)
    seed = autoBraggSearch.UBSeed(np.eye(3), hkls[0], peaks[0], 0.0, 0.0)
    for peak, hkl in zip(peaks[:2], hkls[:2]):
        gui.reflectionSel.addReflection(dict(zip(("x", "y"), peak.xy)),
                                        peak.imageno, hkl)
    added = workflow._addAutomaticObservedReflections(
        peaks, seed, {"peak": peaks[1], "hkl": hkls[1]}, 1e-7,
        intensity_ratio_check=False,
    )
    assert added == 1
    measured_hkls, measured_positions = gui.getReflections()
    np.testing.assert_array_equal(measured_hkls, hkls)
    np.testing.assert_allclose(measured_positions, positions, atol=1e-12)
    np.testing.assert_allclose(gui.getReflectionFitData()[1], positions, atol=1e-12)
    selector = SimpleNamespace(orparent=gui, ubcalc=gui.ubcalc)
    for reflection, position in zip(gui.reflectionSel.reflections, positions):
        np.testing.assert_allclose(
            QReflectionSelector.anglesReflection(selector, reflection), position,
            atol=1e-12,
        )
    gui.ubcalc.ubCal.setU(HKLVlieg.calcOMEGA(0.1))
    workflow._recalculateAutomaticUFromReflections()
    np.testing.assert_allclose(gui.ubcalc.ubCal.getU(), np.eye(3), atol=1e-10)


def test_confirmation_refinement_starts_at_frame_arm_prediction(monkeypatch):
    """Start confirmation refinement at the frame-specific detector prediction."""
    gui, hkls, _positions, peaks = _indexed_peaks()
    gui.background_image, gui.numberthreads = None, 1
    workflow = AutoBraggWorkflow(gui)
    workflow._automaticReflectionIntensity = lambda _hkl: 1.0
    q_limits = []

    def allowed(_xtal, max_q):
        """Return selected reflection candidates for the test geometry."""
        q_limits.append(max_q)
        return hkls, np.ones(3), np.ones(3)

    monkeypatch.setattr(autoBraggSearch, "allowed_bragg_with_intensity", allowed)
    monkeypatch.setattr(
        autoBraggSearch, "estimate_rocking_intensity",
        lambda _peak: {"intensity": 100., "prominence_z": 10.},
    )
    attempted = []

    def refine(_scan, candidate, **_kwargs):
        """Capture the predicted starting point for peak refinement."""
        attempted.append(candidate)
        return peaks[candidate.imageno]

    monkeypatch.setattr(autoBraggSearch, "refine_peak_3d", refine)
    seed = autoBraggSearch.UBSeed(np.eye(3), hkls[0], peaks[0], 0., 0.)
    confirmation = workflow._confirmAutomaticBraggSeed(
        seed, pixel_tolerance=1e-5, image_tolerance=0,
    )
    assert confirmation is not None
    np.testing.assert_array_equal(confirmation["hkl"], hkls[1])
    np.testing.assert_allclose(attempted[0].xy, peaks[1].xy, atol=1e-8)
    assert confirmation["pixel_error"] < 1e-8
    assert q_limits == [max(
        gui.ubcalc.detectorCal.qRangeAtArm(*gui.getArmAngles(index))[1]
        for index in range(len(gui.fscan))
    )]


def test_calculated_candidates_include_reflections_beyond_calibration(monkeypatch):
    """Do not drop peaks before the automatic additional-reflection pass."""
    gui, hkls, positions, _peaks = _indexed_peaks()
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(1000, 50, 50, pixelX=100, pixelY=100)
    detector.detector.shape = (100, 100)
    detector.wavelength = gui.ubcalc.ubCal.getLambda() * 1e-10
    gui.ubcalc.detectorCal = detector
    # Park each frame's narrow detector directly on its known outgoing ray.
    primary_positions = np.array([HKLVlieg.primBeamAngles(p) for p in positions])
    gui.fscan.gamma_arm = np.rad2deg(primary_positions[:, 2])
    gui.fscan.delta_arm = np.rad2deg(primary_positions[:, 1])
    norms = np.linalg.norm((gui.ubcalc.crystal.B_mat @ hkls.T).T, axis=1)
    assert np.min(norms) > detector.Qmax
    gui.ubcalc.crystal.setEnergy = lambda _energy: None

    def allowed(_xtal, maxQ, **_kwargs):
        """Return selected reflection candidates for the test geometry."""
        return hkls[norms < maxQ]

    monkeypatch.setattr(rn, "allowedReflections_G", allowed)
    captured = []
    gui.reflectionSel.setBraggReflections = lambda *args: captured.append(args)
    gui.calcBraggRefl()
    assert captured
    candidate_hkls, yx, candidate_angles = captured[0]
    # Only the middle omega is strictly inside this scan's existing bounds.
    np.testing.assert_array_equal(candidate_hkls, hkls[:1])
    np.testing.assert_allclose(candidate_angles, positions[:1], atol=1e-12)
    expected_yx = detector.pixelsSurfaceAngles(
        positions[0, 2], positions[0, 1], positions[0, 0], *gui.getArmAngles(0)
    )
    np.testing.assert_allclose(yx, expected_yx, atol=1e-8)
    assert np.all((yx >= 0) & (yx < 100))
    # The same candidate vanishes through the old calibration-only path.
    legacy_hkls, _yx, _angles = rn.thscanBragg(
        gui.ubcalc.crystal, gui.ubcalc.ubCal, gui.ubcalc.mu, detector,
        (positions[:, 3].min(), positions[:, 3].max()),
        chi=gui.ubcalc.chi, phi=gui.ubcalc.phi,
    )
    assert len(legacy_hkls) == 0


def test_fixed_arm_calculated_candidates_match_calibration_path(monkeypatch):
    """Retain the legacy candidate set when the arm stays at calibration."""
    gui, hkls, positions, _peaks = _indexed_peaks()
    gui.fscan.gamma_arm, gui.fscan.delta_arm = 2.0, -4.0
    gui.ubcalc.crystal.setEnergy = lambda _energy: None
    norms = np.linalg.norm((gui.ubcalc.crystal.B_mat @ hkls.T).T, axis=1)
    monkeypatch.setattr(
        rn, "allowedReflections_G",
        lambda _xtal, maxQ, **_kwargs: hkls[norms < maxQ],
    )
    captured = []
    gui.reflectionSel.setBraggReflections = lambda *args: captured.append(args)
    gui.calcBraggRefl()
    assert captured
    legacy = rn.thscanBragg(
        gui.ubcalc.crystal, gui.ubcalc.ubCal, gui.ubcalc.mu,
        gui.ubcalc.detectorCal,
        (positions[:, 3].min(), positions[:, 3].max()),
        chi=gui.ubcalc.chi, phi=gui.ubcalc.phi,
    )
    assert len(legacy[0]) > 0
    for current, calibrated in zip(captured[0], legacy):
        np.testing.assert_allclose(current, calibrated, atol=1e-10)


def test_mu_scan_reference_rows_retain_each_incidence_angle():
    """Keep each mu-scan reference aligned with its own incidence angle."""
    scan = SimulationScan((1, 1), 0.3, 0.9, 3, axis="mu", fixed=12.5)
    scan.gamma_arm, scan.delta_arm = [1., 4., 7.], [2., 7., 11.]
    gui = _gui(scan, _detector())
    expected = []
    for index in range(3):
        mu, omega = gui.getMuOm(index)
        pos = [mu, 0.05, 0.03, omega, gui.ubcalc.chi, gui.ubcalc.phi]
        yx = gui.ubcalc.detectorCal.pixelsSurfaceAngles(
            pos[2], pos[1], mu, *gui.getArmAngles(index)
        )[0]
        gui.reflectionSel.addReflection({"x": yx[1], "y": yx[0]}, index, [1, 1, 1])
        expected.append(pos)
    np.testing.assert_allclose(gui.getReflections()[1], expected, atol=1e-12)


@pytest.mark.parametrize("arm_metadata", ["moving", "surface", "fixed_home"])
def test_experiment_fit_uses_fixed_measured_arms(arm_metadata, monkeypatch):
    """Recover distance from lab-ray data while preserving known arm motors."""
    home = np.deg2rad([2., -4.])
    true_detector = _detector(home)
    starting_detector = copy.deepcopy(true_detector)
    # pyFAI's custom deepcopy does not preserve orGUI's arm reference.
    starting_detector.setArmReference(true_detector.getArmReference())
    starting_detector.dist *= 0.85
    starting_distance = starting_detector.dist
    scan = SimulationScan((1, 1), -20., 20., 5, fixed=0.4)
    gui = _gui(scan, starting_detector, mu=np.deg2rad(0.4))
    if arm_metadata == "fixed_home":
        # The configured fixed arm is exactly where calibration was measured.
        gui.ubcalc.gamma_arm, gui.ubcalc.delta_arm = home
        arm_readbacks = np.broadcast_to(home, (5, 2)).copy()
    else:
        scan.gamma_arm = np.array([1., 3., 6., 8., -2.])
        scan.delta_arm = np.array([2., 7., 11., -5., 4.])
        if arm_metadata == "surface":
            scan.arm_angle_frame = "surface"
        arm_readbacks = np.deg2rad(np.column_stack((scan.gamma_arm, scan.delta_arm)))
    xy = np.array([[180., 210.], [570., 300.], [220., 700.],
                   [620., 820.], [420., 450.]])
    gamma0, delta0 = true_detector.primBeamPoints(
        xy[:, 1].copy(), xy[:, 0].copy(), *home
    )
    home_rays = np.column_stack((np.sin(delta0) * np.cos(gamma0),
                                 np.cos(delta0) * np.cos(gamma0), np.sin(gamma0)))
    expected_u = HKLVlieg.Rotation.from_euler("xyz", [.1, -.2, .3]).as_matrix()
    expected_arms = []
    q_phi = []
    for index, (gamma_arm, delta_arm) in enumerate(arm_readbacks):
        mu, omega = gui.getMuOm(index)
        if arm_metadata == "surface":
            ray = (HKLVlieg.calcALPHA(mu) @ HKLVlieg.calcDELTA(delta_arm)
                   @ HKLVlieg.calcGAMMA(gamma_arm) @ np.array([0., 1., 0.]))
            gamma_arm, delta_arm = np.arcsin(ray[2]), np.arctan2(ray[0], ray[1])
        expected_arms.append((gamma_arm, delta_arm))
        outgoing = (HKLVlieg.calcDELTA(delta_arm) @ HKLVlieg.calcGAMMA(gamma_arm)
                    @ (HKLVlieg.calcDELTA(home[1])
                       @ HKLVlieg.calcGAMMA(home[0])).T @ home_rays[index])
        q_lab = gui.ubcalc.ubCal.getK() * (outgoing - [0., 1., 0.])
        q_phi.append(HKLVlieg.calcPHI(gui.ubcalc.phi).T
                     @ HKLVlieg.calcCHI(gui.ubcalc.chi).T
                     @ HKLVlieg.calcOMEGA(omega).T
                     @ HKLVlieg.calcALPHA(mu).T @ q_lab)
    expected_arms = np.asarray(expected_arms)
    hkls = np.linalg.solve(expected_u @ gui.ubcalc.crystal.B_mat,
                           np.asarray(q_phi).T).T
    # Scramble reference order and interleave a stale reference: the snapshot
    # must keep arm positions aligned with the returned pixels and HKLs.
    order = [3, 0, 4, 1, 2]
    for index in order:
        gui.reflectionSel.addReflection(
            {"x": xy[index, 0], "y": xy[index, 1]}, index, hkls[index]
        )
    gui.reflectionSel.reflections.insert(
        2, SimpleNamespace(xy=np.zeros(2), imageno=99, hkl=np.ones(3))
    )
    assert len(gui.getReflectionFitData()) == 3
    snapshot = gui.getReflectionFitData(include_arm_angles=True)
    np.testing.assert_array_equal(snapshot[0], hkls[order])
    np.testing.assert_allclose(snapshot[3], expected_arms[order], atol=1e-14)
    gui.reflectionSel.updateEditor = lambda: None
    gui.ubcalc.crystalparams = SimpleNamespace(setValues=lambda *_args: None)
    gui.ubcalc.machineParams = SimpleNamespace(setValues=lambda *_args: None)
    fixed_config = gui.ubcalc.gamma_arm, gui.ubcalc.delta_arm
    fixed_scan_arms = copy.deepcopy(arm_readbacks)
    fixed_reference = starting_detector.getArmReference()
    starting_geometry = np.asarray(starting_detector.param).copy()
    trial_calls = []
    convert = starting_detector.surfaceAnglesPointParam

    def trial_angles(rows, columns, alpha, param, gamma_arm=None, delta_arm=None):
        # Every objective call sees the original shared geometry and the same
        # arm reference; only the explicit trial parameter array can change.
        """Check shared geometry stays fixed during trial-angle conversion."""
        np.testing.assert_array_equal(starting_detector.param, starting_geometry)
        np.testing.assert_array_equal(
            starting_detector.getArmReference(), fixed_reference
        )
        trial_calls.append((np.asarray(param).copy(), gamma_arm, delta_arm))
        return convert(rows, columns, alpha, param, gamma_arm, delta_arm)

    monkeypatch.setattr(starting_detector, "surfaceAnglesPointParam", trial_angles)
    from orgui.app import QUBCalculator as qub_module

    least_squares = qub_module.optimize.least_squares
    optimizer_calls = []

    def checked_optimizer(fun, x0, **kwargs):
        """Verify only the requested detector distance is optimized."""
        optimizer_calls.append(np.asarray(x0).copy())
        # Distance is the sole requested fit variable; arm gamma/delta add no
        # entries to x0 or the optimizer's bounds.
        np.testing.assert_array_equal(x0, [starting_distance])
        np.testing.assert_array_equal(kwargs["bounds"][0], [np.finfo(float).eps])
        np.testing.assert_array_equal(kwargs["bounds"][1], [35.])
        np.testing.assert_allclose(fun([true_detector.dist]), 0., atol=1e-12)
        return least_squares(fun, x0, **kwargs)

    monkeypatch.setattr(qub_module.optimize, "least_squares", checked_optimizer)
    rms = QUBCalculator.fitExperiment(gui.ubcalc, {"dist"})
    assert len(optimizer_calls) == 1
    assert rms < 1e-10
    assert starting_detector.dist == pytest.approx(true_detector.dist, abs=1e-9)
    np.testing.assert_allclose(gui.ubcalc.ubCal.getU(), expected_u, atol=1e-9)
    np.testing.assert_array_equal(starting_detector.getArmReference(), fixed_reference)
    np.testing.assert_array_equal(
        (gui.ubcalc.gamma_arm, gui.ubcalc.delta_arm), fixed_config
    )
    if arm_metadata != "fixed_home":
        np.testing.assert_array_equal(
            np.deg2rad(np.column_stack((scan.gamma_arm, scan.delta_arm))),
            fixed_scan_arms,
        )
    actual_arms = np.unique(np.array([call[1:] for call in trial_calls]), axis=0)
    np.testing.assert_allclose(
        actual_arms, np.unique(expected_arms, axis=0), atol=1e-14
    )
    assert len(np.unique([call[0][0] for call in trial_calls])) > 1
    for parameter, _gamma, _delta in trial_calls:
        np.testing.assert_array_equal(parameter[1:], starting_geometry[1:6])

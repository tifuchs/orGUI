import numpy as np

from orgui.app import autoBraggSearch
from orgui.backend.scans import SimulationScan
from orgui.datautils.xrayutils import CTRcalc, HKLVlieg


class CountingScan(SimulationScan):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.read_images = []

    def get_raw_img(self, i):
        self.read_images.append(i)
        return super().get_raw_img(i)


def test_scan_image_maxima_ranks_sharp_post_burn_in_peak():
    scan = SimulationScan((8, 8), 0.0, 9.0, 10)
    for i in range(len(scan)):
        scan.images[i] = np.full((8, 8), 10.0 + i)
    scan.images[7, 3, 4] = 500.0

    maxima = autoBraggSearch.scan_image_maxima(scan, burn_in=3, history=4)

    assert maxima[0].imageno == 7
    np.testing.assert_allclose(maxima[0].xy, [4.0, 3.0])
    assert maxima[0].sharpness > 100.0


def test_scan_image_maxima_rejects_pixels_near_mask():
    scan = SimulationScan((8, 8), 0.0, 9.0, 10)
    scan.images[7, 3, 4] = 500.0
    scan.images[8, 6, 6] = 400.0
    mask = np.zeros((8, 8), dtype=bool)
    mask[3, 5] = True

    maxima = autoBraggSearch.scan_image_maxima(
        scan, mask=mask, burn_in=1, history=4, mask_distance=1
    )

    assert all(maximum.imageno != 7 for maximum in maxima)


def test_iter_sharp_peak_candidates_yields_before_scan_end():
    scan = CountingScan((8, 8), 0.0, 19.0, 20)
    scan.images[:] = 10.0
    scan.images[7, 3, 4] = 500.0

    candidates = autoBraggSearch.iter_sharp_peak_candidates(
        scan,
        burn_in=3,
        history=4,
        min_history=4,
        level_z=6.0,
        derivative_z=6.0,
        lookahead=1,
    )

    candidate = next(candidates)

    assert candidate.imageno == 7
    np.testing.assert_allclose(candidate.xy, [4.0, 3.0])
    assert max(scan.read_images) == 8
    assert 9 not in scan.read_images


def test_iter_sharp_peak_candidates_can_continue_after_rejected_candidate():
    scan = CountingScan((8, 8), 0.0, 19.0, 20)
    scan.images[:] = 10.0
    scan.images[7, 3, 4] = 500.0
    scan.images[12, 5, 6] = 700.0

    candidates = autoBraggSearch.iter_sharp_peak_candidates(
        scan,
        burn_in=3,
        history=4,
        min_history=4,
        level_z=6.0,
        derivative_z=6.0,
        lookahead=1,
        refractory=1,
    )

    first = next(candidates)
    second = next(candidates)

    assert first.imageno == 7
    assert second.imageno == 12
    assert max(scan.read_images) == 13


def test_seed_u_from_single_reflection_uses_existing_zmode_solver():
    lattice = HKLVlieg.Lattice([3.0, 4.0, 5.0], [90.0, 90.0, 90.0])
    ub = HKLVlieg.UBCalculator(lattice, 70.0)
    ub.defaultU()
    hkl = np.array([1.0, 1.0, 0.0])
    pos = np.array(
        [
            np.deg2rad(0.2),
            np.deg2rad(8.0),
            np.deg2rad(4.0),
            np.deg2rad(-12.0),
            0.0,
            0.0,
        ]
    )

    expected = HKLVlieg.UBCalculator(lattice, 70.0)
    expected.defaultU()
    expected.zmodeUSingleRefl(pos, hkl)

    seeded_u = autoBraggSearch.seed_u_from_single_reflection(ub, pos, hkl)

    np.testing.assert_allclose(seeded_u, expected.getU(), atol=1e-12)


def test_rank_ub_seeds_prefers_reflection_with_matching_q_norm():
    xtal = CTRcalc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
    xtal.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    xtal.setEnergy(70000.0)
    ub = HKLVlieg.UBCalculator(xtal, 70.0)
    ub.defaultU()
    hkl = np.array([1.0, 0.0, 0.0])
    q_norm = np.linalg.norm(xtal.B_mat @ hkl)
    delta = 2.0 * np.arcsin(q_norm / (2.0 * ub.getK()))

    class FakeDetector:
        def surfaceAnglesPoint(self, y, x, mu):
            return np.array([0.0]), np.array([delta])

    peak = autoBraggSearch.RefinedPeak(
        candidate=autoBraggSearch.ImageMaximum(
            imageno=0, xy=np.array([2.0, 3.0]), value=100.0
        ),
        xy=np.array([2.0, 3.0]),
        imageno=0,
        axis_value=0.0,
        value=100.0,
    )

    seeds = autoBraggSearch.rank_ub_seeds(
        [peak],
        FakeDetector(),
        lambda _imageno: (0.0, 0.0),
        ub,
        xtal,
        0.0,
        0.0,
        qnorm_tolerance=1e-6,
        max_q=3.0,
    )

    assert seeds
    np.testing.assert_allclose(
        np.linalg.norm(xtal.B_mat @ seeds[0].hkl),
        q_norm,
    )


def test_allowed_bragg_in_q_region_sorts_by_intensity_and_deduplicates():
    xtal = CTRcalc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
    xtal.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    xtal.setEnergy(70000.0)
    qnorm = np.linalg.norm(xtal.B_mat @ np.array([1.0, 0.0, 0.0]))

    hkls, intensity, norms = autoBraggSearch.allowed_bragg_in_q_region(
        xtal, qnorm, 1e-6, max_q=3.0
    )

    assert len(hkls) == 3
    assert intensity[0] > 0.0
    np.testing.assert_allclose(norms, qnorm)


def test_allowed_bragg_same_q_returns_full_same_norm_group():
    xtal = CTRcalc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
    xtal.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    xtal.setEnergy(70000.0)
    qnorm = np.linalg.norm(xtal.B_mat @ np.array([1.0, 0.0, 0.0]))

    hkls, intensity, norms = autoBraggSearch.allowed_bragg_same_q(
        xtal, qnorm, tolerance=1e-6, max_q=3.0
    )

    assert len(hkls) == 3
    assert intensity[0] > 0.0
    np.testing.assert_allclose(norms, qnorm)
    assert {tuple(hkl) for hkl in hkls} == {
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, 1),
    }


def test_deduplicate_only_merges_sign_inverted_equal_f_and_q():
    hkls = np.array(
        [
            [-1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 1.0],
            [1.0, 2.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ]
    )
    structure_factor = np.array(
        [
            2.0 + 0.0j,
            2.0 + 0.0j,
            3.0 + 0.0j,
            2.0 + 0.0j,
            2.0 + 0.0j,
            2.0 + 0.0j,
        ]
    )
    norms = np.array([1.5, 1.5, 1.5, 2.0, 2.0, 1.6])

    deduped_hkls, _structure_factor, deduped_norms = (
        autoBraggSearch._deduplicate_structure_equivalent(hkls, structure_factor, norms)
    )

    assert len(deduped_hkls) == 5
    assert any(np.allclose(hkl, [-1.0, 1.0, 1.0]) for hkl in deduped_hkls)
    assert any(np.allclose(hkl, [1.0, -1.0, 1.0]) for hkl in deduped_hkls)
    assert any(np.allclose(hkl, [1.0, 2.0, 0.0]) for hkl in deduped_hkls)
    assert any(np.allclose(hkl, [2.0, 1.0, 0.0]) for hkl in deduped_hkls)
    kept_shifted_q = [
        norm
        for hkl, norm in zip(deduped_hkls, deduped_norms)
        if np.allclose(hkl, [1.0, 1.0, 1.0])
    ]
    np.testing.assert_allclose(kept_shifted_q, [1.6])


def test_allowed_bragg_in_q_region_accepts_reciprocal_scale():
    xtal = CTRcalc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
    xtal.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    xtal.setEnergy(70000.0)
    configured_qnorm = np.linalg.norm(xtal.B_mat @ np.array([1.0, 0.0, 0.0]))
    measured_qnorm = 1.08 * configured_qnorm

    hkls_without_scale, _intensity, _norms = autoBraggSearch.allowed_bragg_in_q_region(
        xtal, measured_qnorm, 1e-6, max_q=3.0
    )
    hkls_with_scale, _intensity, norms = autoBraggSearch.allowed_bragg_in_q_region(
        xtal,
        measured_qnorm,
        1e-6,
        max_q=3.0,
        qnorm_scale=1.08,
    )

    assert len(hkls_without_scale) == 0
    assert len(hkls_with_scale) == 3
    np.testing.assert_allclose(norms, configured_qnorm)


def test_estimate_qnorm_scale_reports_direct_lattice_scale():
    xtal = CTRcalc.UnitCell([4.0, 4.0, 4.0], [90.0, 90.0, 90.0])
    xtal.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
    xtal.setEnergy(70000.0)
    ub = HKLVlieg.UBCalculator(xtal, 70.0)
    ub.defaultU()
    configured_qnorm = np.linalg.norm(xtal.B_mat @ np.array([1.0, 0.0, 0.0]))
    reciprocal_scale = 1.08
    delta = 2.0 * np.arcsin(reciprocal_scale * configured_qnorm / (2.0 * ub.getK()))

    class FakeDetector:
        def surfaceAnglesPoint(self, y, x, mu):
            return np.array([0.0]), np.array([delta])

    peaks = [
        autoBraggSearch.RefinedPeak(
            candidate=autoBraggSearch.ImageMaximum(
                imageno=i, xy=np.array([2.0, 3.0]), value=100.0
            ),
            xy=np.array([2.0, 3.0]),
            imageno=i,
            axis_value=0.0,
            value=100.0,
        )
        for i in range(3)
    ]

    estimate = autoBraggSearch.estimate_qnorm_scale(
        peaks,
        FakeDetector(),
        lambda _imageno: (0.0, 0.0),
        ub,
        xtal,
        0.0,
        0.0,
        tolerance=0.2,
        max_q=3.0,
    )

    assert estimate is not None
    np.testing.assert_allclose(estimate["qnorm_scale"], reciprocal_scale, rtol=1e-12)
    np.testing.assert_allclose(estimate["direct_lattice_scale"], 1.0 / reciprocal_scale)


def test_estimate_rocking_intensity_reports_prominence_z():
    peak = autoBraggSearch.RefinedPeak(
        candidate=autoBraggSearch.ImageMaximum(
            imageno=0, xy=np.array([2.0, 3.0]), value=100.0
        ),
        xy=np.array([2.0, 3.0]),
        imageno=0,
        axis_value=0.0,
        value=100.0,
        rocking_curve=np.array([10.0, 11.0, 10.0, 60.0, 11.0, 10.0, 11.0]),
        rocking_axis=np.arange(7.0),
    )

    estimate = autoBraggSearch.estimate_rocking_intensity(peak)

    assert estimate is not None
    assert estimate["intensity"] > 0.0
    assert estimate["prominence_z"] > 6.0

"""Polarization extraction, approximation and configuration regressions."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from orgui.app import scan_integration as integration
from orgui.app.config_data import (
    CorrectionState, corrections_from_nxdict, corrections_to_nxdict,
)
from orgui.app.test.test_scan_integration import context, _prepared_line  # noqa: F401

# Pytest injects the imported fixture into the matching test parameter.
# ruff: noqa: F811


@pytest.mark.parametrize("repair_enabled", [False, True])
def test_polarization_counter_matches_previous_pass(repair_enabled):
    """Fixed acceptance reproduces native repaired and NumPy masked sums."""
    if repair_enabled and integration._roi_sum_accel.ROI_ACCEL_BACKEND != "cpp":
        pytest.skip("Pixel repair requires the C++ backend")
    mask = np.zeros((20, 24), dtype=bool)
    mask[8, 8] = True
    mask[7, 5] = True
    mask[2:5, 17:20] = True
    correction = 1 + np.arange(mask.size).reshape(mask.shape) / 1000
    repair = SimpleNamespace(max_component_pixels=4, max_span=3,
                             radius=2, min_valid_neighbors=6)
    gaps = np.array([[17, 19]], dtype=np.int32)
    center = (slice(6, 12), slice(6, 12))
    backgrounds = [(slice(4, 7), slice(6, 12))] * 4
    count = integration._prepare_polarization_counter(
        correction, mask, repair if repair_enabled else None, gaps, gaps,
    )
    if repair_enabled:
        arrays = [np.array([[[s.start, s.stop] for s in region]], dtype=np.int64)
                  for region in [center, *backgrounds]]
        expected = np.zeros((1, 4))
        integration._roi_sum_accel.processImage_repair_Carr(
            np.ones(mask.shape), mask, correction, *arrays, gaps, gaps,
            np.zeros((1, 4)), expected, repair.max_component_pixels,
            repair.max_span, repair.radius, repair.min_valid_neighbors,
        )
        expected = expected[0]
    else:
        expected = integration._correction_region_counters(
            correction, mask, center, backgrounds,
        )
    np.testing.assert_allclose(count(center, backgrounds), expected, rtol=1e-14)


@pytest.mark.parametrize("workflow", ["single", "lines", "rocking"])
@pytest.mark.parametrize("mode", ["pixel", "center"])
@pytest.mark.parametrize("backend", ["numpy", "cpp", "repair"])
def test_polarization_extraction(context, monkeypatch, workflow, mode, backend):
    """Factors survive all routes without a post-extraction repair pass."""
    if backend != "numpy" and integration._roi_sum_accel.ROI_ACCEL_BACKEND != "cpp":
        pytest.skip("C++ backend unavailable")
    monkeypatch.setattr(integration, "HAS_ACCEL", backend != "numpy")
    repair_enabled = backend == "repair"
    ctx, scan, file = context
    mask = np.zeros((12, 12), dtype=bool)
    mask[2, 2] = True
    gaps = np.empty((0, 2), dtype=np.int32)
    repair = SimpleNamespace(max_component_pixels=4, max_span=3,
                             radius=1, min_valid_neighbors=2)
    options = ctx.scanSelector.get_integration_options()
    options.update(polarization=True, polarization_mode=mode,
                   solid_angle=True, mask=True)
    polarized_requests = []

    def factors(detector, solid_angle=False, polarization=False):
        """Distinct constant factors detect omission and double application."""
        polarized_requests.append(polarization)
        if not solid_angle and not polarization:
            return None
        return np.full((12, 12), (3. if solid_angle else 1.)
                       * (2. if polarization else 1.))

    monkeypatch.setattr(integration.detector_corrections, "pixel_factors", factors)
    monkeypatch.setattr(integration.detector_corrections,
                        "polarization_arm_correction", lambda *a, **k: 1.)
    dc = ctx.ubcalc.detectorCal
    dc.polarizationAtPointsFrames = lambda row, col, alpha, *a: np.full(
        (len(alpha), len(row)), .5)
    first = np.zeros((len(scan), 9))
    first[:, 2] = first[:, 5] = scan.axis
    first[:, 3:5] = .1
    first[:, 6:8] = [3., 3.]
    first[:, -1] = 1.
    second = first.copy()
    second[:, -1] = 0.
    ctx = replace(
        ctx, getROIloc=lambda: (first, second),
        get_detector_mask=lambda shape: mask,
        _repair_config_for_image=lambda shape: (repair_enabled, repair, gaps, gaps),
        scanSelector=integration.prepared_selector(
            options, 0, [0, 0, 0], [0, 0, 1], [3, 3], None),
    )
    repairs = []
    original = integration._roi_sum_accel.repair_masked_pixels_inplace

    def counted(*args):
        """Record acceptance construction, which must not scale with frames."""
        repairs.append(True)
        return original(*args)

    monkeypatch.setattr(integration._roi_sum_accel, "repair_masked_pixels_inplace",
                        counted)
    if workflow == "single":
        result = integration.integrate_stationary_scan(ctx)
    elif workflow == "lines":
        result = integration.integrate_stationary_scan(
            ctx, lines=[{"id": "a", "H_0": [0, 0, 0], "H_1": [0, 0, 1]}],
            geometry_provider=lambda line: (first, second),
        )
    else:
        result = integration.integrate_rocking_scan(ctx, lines=[_prepared_line(ctx)])
    assert result["status"] == "success", result
    assert len(repairs) == int(repair_enabled and mode == "pixel")
    assert scan.reads == {i: 1 for i in range(len(scan))}
    if mode == "center":
        assert not any(polarized_requests)
    group = next(iter(file["scan/measurement"].values()))
    pixel = group["ctr_curve_v3/pixel_corrections"]
    np.testing.assert_allclose(pixel["polarization_croi"][()], 2.)
    np.testing.assert_allclose(pixel["combined_croi"][()], 6.)


def test_polarization_mode_persistence_and_legacy_default():
    """Old configs retain pixel averaging and new configs retain the mode."""
    assert corrections_from_nxdict({}).polarization_mode == "pixel"
    state = CorrectionState(use_polarization=True, polarization_mode="center")
    restored = corrections_from_nxdict(corrections_to_nxdict(state))
    assert restored.polarization_mode == "center"
    assert restored.use_polarization

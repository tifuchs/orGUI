"""One simulated detector data set reduced in angular and HKL coordinates.

The forward model is a normalized Gaussian rod cross-section in HKL, with
literal SI scattering scale, pixel solid angles and exposure fluence. It does
not call the production intensity or structure-factor forward functions.
Actual detector corrections, native mapping, checkpoints and HDF5 voxel
means are used before the reciprocal-space integration is tested.
"""

from dataclasses import dataclass
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from orgui.app import integration_corrections as ic
from orgui.app.config_data import CorrectionState
from orgui.app.peak1Dintegr import _compute_rocking_integration
from orgui.backend.scans import h5_Image
from orgui.datautils.xrayutils import DetectorCalibration, HKLVlieg
from orgui.datautils.xrayutils.corrections import acceptance, measurement
from orgui.datautils.xrayutils.reconstruction import (
    _CheckpointRouter, _GridSpec, _ReconstructionSpec, _build_kernels,
    _finalize_reconstruction, _map_frame_group, _tile_ray_arrays,
)
from orgui.reconstruction_job import _correction_pipeline


pytest.importorskip("orgui.datautils.xrayutils._reciprocal_reconstruction_cpp")


@dataclass
class _Scan:
    exposure_time: np.ndarray
    monitor: np.ndarray
    auxillary_counters = ("monitor",)

    def __len__(self):
        return self.exposure_time.size


def _simulate(ell, f2, monitor_kind, *, reverse=False):
    """Return physical detector frames sampling the same Gaussian CTR."""
    lattice = HKLVlieg.Lattice([2.7748, 2.7748, 6.7964], [90, 90, 120])
    ub = HKLVlieg.UBCalculator(lattice, 17.7)
    ub.defaultU_GID()
    angles = HKLVlieg.VliegAngles(ub)
    position = angles.anglesZmode(
        np.array([[1.0], [0.0], [ell]]), np.deg2rad(0.6), fixed="in"
    )[0]
    alpha, delta, gamma, omega, chi, phi = position
    _, delta_p, gamma_p, *_ = HKLVlieg.primBeamAngles(position)

    shape, pixel_m, distance_m = (241, 161), 172e-6, 4.0
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.detector = pytest.importorskip("pyFAI").detectors.Detector(
        pixel1=pixel_m, pixel2=pixel_m, max_shape=shape
    )
    detector.poni1, detector.poni2 = np.asarray(shape) * pixel_m / 2
    detector.dist = distance_m
    detector.rot1 = detector.rot2 = detector.rot3 = 0.0
    detector.setAzimuthalReference(np.pi / 2)
    detector.setPolarization(0.0, 1.0)
    detector.rot1, detector.rot2, detector.rot3 = detector.paramAtArm(
        float(gamma_p), float(delta_p)
    )[3:6]
    detector.setArmReference(gamma_arm=float(gamma_p), delta_arm=float(delta_p))
    detector.reset()
    detector._cached_array = {}
    rows, columns = np.indices(shape, dtype=float)
    gg, dd = detector.surfaceAnglesPoint(rows, columns, alpha)
    # A rigidly rotated flat detector has the same absolute solid angles.
    yy = (rows + 0.5) * pixel_m - detector.poni1
    xx = (columns + 0.5) * pixel_m - detector.poni2
    solid_angle = pixel_m**2 * distance_m / (
        distance_m**2 + xx**2 + yy**2
    )**1.5
    polarization = 1.0 - (
        np.sin(alpha) * np.cos(dd) * np.cos(gg)
        + np.cos(alpha) * np.sin(gg)
    )**2
    np.testing.assert_allclose(detector.polarizationArray(), polarization, rtol=2e-7)

    axis = np.linspace(-0.006, 0.006, 161)
    if reverse:
        axis = axis[::-1]
    exposure = 0.3 + 0.7 * np.arange(axis.size) / (axis.size - 1)
    relative_flux = 1.0 + 0.19 * np.sin(axis / 0.001)
    flux, reference, reference_exposure = 3.7e10, 1.8e6, 0.5
    monitor = reference * relative_flux
    if monitor_kind == "integrated":
        monitor = monitor * exposure / reference_exposure
    scan = _Scan(exposure, monitor)
    # Both monitor conventions represent the same actual photons per frame.
    fluence = flux * exposure * relative_flux
    state = CorrectionState(
        use_normalization=True, shared_frame_normalization=True,
        use_solid_angle=True, use_polarization=True,
        total_incident_flux=flux, total_flux_calibrated=True,
        primary_monitor="monitor", primary_monitor_kind=monitor_kind,
        primary_monitor_unit="counts/s" if monitor_kind == "rate" else "counts",
        monitor_reference_reading=reference,
        monitor_reference_exposure_s=(
            reference_exposure if monitor_kind == "integrated" else None
        ),
    )
    area_m2, electron_radius_m = lattice.uc_area * 1e-20, 2.8179403262e-15
    illumination, sigma_hk = 1.3, 0.0006
    raw, bounds = [], []
    half_step = (axis[1] - axis[0]) / 2
    # Counts are exposure averages, not point samples. Nine-point midpoint
    # quadrature resolves each small continuous rocking exposure.
    offsets = ((np.arange(9) + 0.5) / 9 - 0.5) * 2 * half_step
    for index, offset in enumerate(axis):
        density = np.zeros(shape)
        for suboffset in offsets:
            h, k, _ = angles.anglesToHkl(
                alpha, dd, gg, omega + offset + suboffset, chi, phi
            )
            density += np.exp(-((h - 1)**2 + k**2) / (2 * sigma_hk**2))
        density /= 9 * 2 * np.pi * sigma_hk**2
        raw.append(fluence[index] * illumination * electron_radius_m**2
                   / area_m2 * f2 * density * polarization * solid_angle)
        bounds.append([
            [alpha, omega + offset - half_step, chi, phi],
            [alpha, omega + offset + half_step, chi, phi],
        ])
    angular_scale = electron_radius_m**2 * (ub.getLambda() * 1e-10)**2 / area_m2**2
    return SimpleNamespace(
        detector=detector, ub=ub, scan=scan, state=state, raw=np.array(raw),
        bounds=np.array(bounds), axis=axis, position=position, fluence=fluence,
        polarization=polarization, illumination=illumination, scale=angular_scale,
        area=lattice.uc_area, reference_solid_angle=pixel_m**2 / distance_m**2,
    )


def _angular_reductions(data):
    """Stationary ROI and rocking aggregation of the same detector frames."""
    alpha, delta, gamma, *_ = data.position
    corrected_pixels = data.raw / data.polarization
    middle = len(data.scan) // 2
    factors = ic.stationary_correction_factors(
        alpha, delta, gamma, use_lorentz=True,
        normalization=data.fluence[middle],
        illumination_divisor=data.illumination,
    )
    stationary, errors = ic.apply_stationary_corrections(
        np.sum(corrected_pixels[middle]), 0.0, factors
    )
    stationary = ic.structure_factor(stationary, errors, factors)[0] / data.scale

    # A shorter vertical aperture keeps the transverse peak fully captured
    # throughout rocking, while the stationary ROI spans the whole profile.
    curve = corrected_pixels[:, 100:141].sum(axis=(1, 2))[None, :]
    dg = acceptance.out_of_plane_acceptance(
        data.detector, 120.0, 80.0, 41.0, alpha
    )
    size = curve.shape
    result = _compute_rocking_integration(
        np.array([0.0]), data.axis, curve, np.zeros(size),
        {"sig_1": {"from": np.array([data.axis.min()]),
                   "to": np.array([data.axis.max()])}}, {}, True, False,
        C_Lor=np.full(size, 1 / (np.sin(delta) * np.cos(alpha) * np.cos(gamma))),
        C_rod=np.full(size, np.cos(gamma)),
        C_norm=data.fluence[None, :] * data.illumination,
        detector_acceptance=np.asarray(dg).reshape(1), angle_unit="rad",
    )
    return float(stationary), float(result["F2_hkl"][0] / data.scale)


def _mapped_reduction(data, ell, tmp_path, *, step, repeats, weighting_mode):
    """Run production corrections, native mapping and HDF5 finalization."""
    pytest.importorskip("orgui.datautils.xrayutils._reciprocal_reconstruction_cpp")
    extent = 0.003
    grid = _GridSpec(
        (1 - extent, -extent, ell - 0.001),
        (1 + extent, extent, ell + 0.001), (step, step, 0.001), "hkl",
    )
    spec = _ReconstructionSpec(
        (grid,), max_depth=2, compression="gzip", weighting_mode=weighting_mode,
    )
    scan = _Scan(np.tile(data.scan.exposure_time, repeats),
                 np.tile(data.scan.monitor, repeats))
    config = SimpleNamespace(detector=data.detector, corrections=data.state)
    provenance = {}
    pipeline = _correction_pipeline(config, scan, {}, provenance)
    assert provenance["normalization_unit"] == "photons"
    tiles = [(0, data.raw.shape[1], 0, data.raw.shape[2])]
    rays = _tile_ray_arrays(data.detector, tiles)
    kernels = _build_kernels(spec, data.ub)
    router = _CheckpointRouter(
        {grid.grid_name: [(0, len(scan))]}, spec_digest=spec.digest,
        checkpoint_dir=tmp_path / "checkpoints", active_budget_bytes=32 * 1024**2,
    )
    for frame in range(len(scan)):
        index = frame % len(data.scan)
        bounds = data.bounds[index]
        _map_frame_group(
            spec, kernels, rays, tiles, pipeline, [h5_Image(data.raw[index])],
            [frame], bounds[None, 0], bounds[None, 1], router,
        )
    output = tmp_path / "map.h5"
    _finalize_reconstruction(spec, {grid.grid_name: router.written}, output)
    with h5py.File(output, "r") as file:
        group = file[f"entry/reconstruction/results/{grid.grid_name}"]
        assert group.attrs["coordinate_frame"] == "hkl"
        assert np.all(group["weight"][()] > 0)  # complete selected volume
        return measurement.reciprocal_map_structure_factor_squared(
            group["intensity"][()], grid.step, unitcell_area=data.area,
            reference_solid_angle=data.reference_solid_angle,
            illumination_divisor=data.illumination,
        )


@pytest.mark.parametrize("ell, f2", [(2.0, 420.0), (3.0, 130.0), (4.0, 900.0)])
@pytest.mark.parametrize("monitor_kind, reverse, step, repeats, weighting_mode", [
    ("rate", False, 0.0002, 1, "reciprocal_volume_average"),
    ("integrated", True, 0.0003, 2, "parameter_average"),
])
def test_all_three_paths_recover_the_same_simulated_rod(
    ell, f2, monitor_kind, reverse, step, repeats, weighting_mode, tmp_path,
):
    """Absolute scale survives grid changes, reversed scans and repeat frames."""
    data = _simulate(ell, f2, monitor_kind, reverse=reverse)
    stationary, rocking = _angular_reductions(data)
    mapped = _mapped_reduction(
        data, ell, tmp_path, step=step, repeats=repeats, weighting_mode=weighting_mode,
    )
    # Finite pixel/sweep resolution, HKL voxelization and the angular ROI's
    # center-ray approximation are distinct quadratures of the same density.
    # The tolerance is 0.1%, not a fitted scale or a per-mode renormalization.
    np.testing.assert_allclose([stationary, rocking, mapped], f2, rtol=1e-3)
    np.testing.assert_allclose([rocking, mapped], stationary, rtol=1e-3)


def test_reciprocal_grid_refinement_is_inside_the_equivalence_tolerance(tmp_path):
    """Voxel refinement resolves the integral before asserting mode agreement."""
    data = _simulate(3.0, 130.0, "rate")
    values = []
    for index, step in enumerate((0.0006, 0.0003, 0.00015)):
        folder = tmp_path / str(index)
        folder.mkdir()
        values.append(_mapped_reduction(
            data, 3.0, folder, step=step, repeats=1,
            weighting_mode="reciprocal_volume_average",
        ))
    np.testing.assert_allclose(values, 130.0, rtol=3e-4)
    assert abs(values[-1] - values[-2]) / 130.0 < 1e-4

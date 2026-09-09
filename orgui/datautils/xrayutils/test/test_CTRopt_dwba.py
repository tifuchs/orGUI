"""Focused integration tests for live DWBA predictions in ``CTROptimizer``."""

from unittest import mock

import numpy as np
import pytest

from .. import CTRcalc, CTRopt, CTRplotutil, CTRresolution, CTRuc


def _crystal(*, fit_x=False):
    cell = CTRuc.UnitCell([3.0, 3.0, 4.0], [90.0, 90.0, 90.0])
    cell.addAtom("C", [0.17, 0.23, 0.31], 0.18, 0.37, 1.0)
    cell.addAtom("O", [0.61, 0.47, 0.72], 0.12, 0.29, 1.0)
    cell.setEnergy(10000.0)
    cell.f[0, 11] = 0.27
    cell.f[0, 12] = 0.19
    if fit_x:
        cell.addFitParameter(([0], "x"), limits=(0.05, 0.45), name="x_position")
    return CTRcalc.SXRDCrystal(cell)


def _generated_angles(
    crystal, ctr, geometry=None, h=None, k=None, l_values=None
):
    if geometry is None:
        geometry = ctr.scan_geometry
    if h is None:
        h = ctr.harr
    if k is None:
        k = ctr.karr
    if l_values is None:
        l_values = ctr.l
    fixed_angle = 0.0 if geometry.fixed == "eq" else geometry.angle
    return CTRplotutil._calculate_angles_zmode(
        h,
        k,
        l_values,
        crystal.dwba._vlieg_angles(),
        fixed_angle,
        fixed=geometry.fixed,
        hkl_transform=crystal.uc_bulk.refHKLTransform,
        mirrorx=geometry.mirrorx,
    )


def _ctr(
    crystal,
    *,
    hk=(0.0, 0.0),
    quantity="reflectivity",
    s_fraction=1.0,
    outgoing="s",
    polarization_factor=None,
    geometry=None,
    records=True,
    keep_rule=True,
    values=None,
    l_values=None,
):
    if geometry is None:
        geometry = CTRplotutil.CTRScanGeometry("eq")
    if l_values is None:
        l_values = np.array([0.12, 0.24, 0.39])
    else:
        l_values = np.asarray(l_values, dtype=np.float64)
    if values is None:
        values = np.array([0.8, 1.1, 1.4])
    reduction = CTRplotutil.MeasurementReduction(
        quantity,
        CTRplotutil.PolarizationReduction(
            s_fraction, outgoing, polarization_factor
        ),
    )
    ctr = CTRplotutil.CTR(
        hk,
        l_values,
        np.asarray(values, dtype=np.float64),
        np.array([0.1, 0.2, 0.3]),
        reduction=reduction,
        scan_geometry=geometry,
    )
    if records:
        ctr.angles = _generated_angles(crystal, ctr)
    if not keep_rule:
        ctr.scan_geometry = None
    return ctr


def _direct_intensity_and_prefactor(state, ctr, angles):
    polarization = ctr.reduction.polarization
    intensity = None
    prefactor = None
    with state.batch():
        for polarization_i, polarization_f, weight in state._polarization_pairs(
            polarization.s_fraction, polarization.outgoing
        ):
            result = state.evaluate_from_vlieg(
                *(angles[name] for name in angles.dtype.names),
                polarization_i=polarization_i,
                polarization_f=polarization_f,
            )
            contribution = weight * result.reflectivity
            intensity = contribution if intensity is None else intensity + contribution
            if prefactor is None:
                prefactor = result.amplitude_prefactor
    return np.asarray(intensity), np.asarray(prefactor)


def test_dwba_configuration_round_trips_and_angle_optimizer_rejects_it():
    """The small public switch round-trips and excludes angle correction."""
    crystal = _crystal()
    ctrs = CTRplotutil.CTRCollection([_ctr(crystal)])
    optimizer = CTRopt.CTROptimizer(crystal, ctrs)

    assert optimizer.get_dwba() == {"enabled": False, "bulk_attenuation": 0.0}
    optimizer.set_dwba(bulk_attenuation=0.25)
    assert optimizer.get_dwba() == {"enabled": True, "bulk_attenuation": 0.25}
    with pytest.raises(TypeError, match="enabled must be boolean"):
        optimizer.set_dwba(enabled="yes")
    with pytest.raises(ValueError, match="finite nonnegative"):
        optimizer.set_dwba(bulk_attenuation=-0.1)

    angle_optimizer = CTRopt.CTROptAngleCorrection(crystal, ctrs)
    with pytest.raises(ValueError, match="not supported.*CTROptAngleCorrection"):
        angle_optimizer.set_dwba()

    with mock.patch.object(
        optimizer.xtal.dwba,
        "evaluate_from_vlieg",
        wraps=optimizer.xtal.dwba.evaluate_from_vlieg,
    ) as evaluate:
        optimizer.prepareFit()
    assert evaluate.call_args.kwargs["bulk_mode"] == "semi_infinite"
    assert evaluate.call_args.kwargs["bulk_attenuation"] == 0.25
    optimizer.set_dwba(False)
    with pytest.raises(RuntimeError, match="requires prepareFit"):
        optimizer.flat_prediction()
    angle_optimizer.useAnglecorr = True
    with pytest.raises(ValueError, match="not supported.*CTROptAngleCorrection"):
        angle_optimizer.set_dwba()


def test_mixed_f_and_r_predictions_keep_observations_in_stored_quantities():
    """F and R are predicted independently without transforming observations."""
    crystal = _crystal()
    factor = np.array([0.73, 0.81, 0.92])
    rods = CTRplotutil.CTRCollection(
        [
            _ctr(
                crystal,
                quantity="structure_factor",
                polarization_factor=factor,
                values=[-1.0, 0.0, 2.0],
            ),
            _ctr(crystal, quantity="reflectivity", values=[0.0, 0.4, 1.2]),
        ]
    )
    optimizer = CTRopt.CTROptimizer(
        crystal, rods, scale_policy={"F": "fixed", "R": "fixed"}
    )
    optimizer.set_dwba()
    stored = [(ctr.sfI.copy(), ctr.err.copy()) for ctr in optimizer.CTRs]
    optimizer.prepareFit()

    expected = []
    for ctr in optimizer.CTRs:
        intensity, prefactor = _direct_intensity_and_prefactor(
            optimizer.xtal.dwba, ctr, ctr.angles
        )
        if ctr.reduction.quantity == "reflectivity":
            expected.append(intensity)
        else:
            expected.append(
                np.sqrt(intensity / ctr.reduction.polarization.polarization_factor)
                / np.abs(prefactor)
            )
    np.testing.assert_allclose(optimizer.flat_prediction(), np.concatenate(expected))
    for ctr, (values, errors) in zip(optimizer.CTRs, stored):
        np.testing.assert_array_equal(ctr.sfI, values)
        np.testing.assert_array_equal(ctr.err, errors)
    np.testing.assert_allclose(
        optimizer.residues(),
        np.concatenate(
            [ctr.sfI - prediction for ctr, prediction in zip(optimizer.CTRs, expected)]
        ),
    )
    assert all(ctr.err is None for ctr in optimizer.calculated_CTRs)


def test_global_reflectivity_scale_multiplies_live_predictions():
    """A shared analytical R scale is fitted after live DWBA prediction."""
    crystal = _crystal()
    rods = [_ctr(crystal), _ctr(crystal)]
    for rod in rods:
        intensity, _ = _direct_intensity_and_prefactor(
            crystal.dwba, rod, rod.angles
        )
        rod.sfI = 2.5 * intensity
    optimizer = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection(rods),
        scale_policy={"R": "global"},
    )
    optimizer.set_dwba()
    optimizer.prepareFit()

    np.testing.assert_allclose(optimizer.flat_prediction(), optimizer.flat_data()[0])
    np.testing.assert_allclose(optimizer.residues(), 0.0, atol=1e-20)
    assert optimizer._fitted_scale_count() == 1


@pytest.mark.parametrize(
    ("s_fraction", "outgoing", "expected_calls"),
    [
        (1.0, "s", 1),
        (1.0, "unanalysed", 2),
        (0.35, "s", 2),
        (0.35, "unanalysed", 4),
    ],
)
def test_dwba_evaluates_only_requested_polarization_pairs(
    s_fraction, outgoing, expected_calls
):
    """Zero-weight incident and unrequested outgoing channels are skipped."""
    crystal = _crystal()
    rod = _ctr(
        crystal,
        s_fraction=s_fraction,
        outgoing=outgoing,
        records=False,
    )
    optimizer = CTRopt.CTROptimizer(
        crystal, CTRplotutil.CTRCollection([rod]), scale_policy={"R": "fixed"}
    )
    optimizer.set_dwba()
    state = optimizer.xtal.dwba
    with mock.patch.object(
        state, "evaluate_from_vlieg", wraps=state.evaluate_from_vlieg
    ) as evaluate:
        optimizer.prepareFit()
    assert evaluate.call_count == expected_calls


def test_one_outer_batch_reuses_the_atomic_snapshot_across_rods():
    """One objective evaluation shares snapshot and packing across rods."""
    crystal = _crystal()
    rods = CTRplotutil.CTRCollection(
        [_ctr(crystal, records=False), _ctr(crystal, records=False)]
    )
    optimizer = CTRopt.CTROptimizer(
        crystal, rods, scale_policy={"R": "fixed"}
    )
    optimizer.set_dwba()
    before = optimizer.xtal.dwba.cache_info()
    optimizer.prepareFit()
    after = optimizer.xtal.dwba.cache_info()
    assert after["snapshot_builds"] - before["snapshot_builds"] == 1
    assert after["packing_builds"] - before["packing_builds"] == 1


def test_measured_grid_resolution_convolves_field_intensity_before_reduction():
    """Measured-grid convolution operates directly on field intensity."""
    crystal = _crystal()
    rod = _ctr(crystal, records=True, keep_rule=False)
    optimizer = CTRopt.CTROptimizer(
        crystal, CTRplotutil.CTRCollection([rod]), scale_policy={"R": "fixed"}
    )
    resolution = CTRresolution.GaussianResolution(0.18)
    optimizer.set_dwba()
    optimizer.set_resolution(resolution, calculation="convolve")
    optimizer.prepareFit()

    fitted_rod = optimizer.CTRs[0]
    intensity, _ = _direct_intensity_and_prefactor(
        optimizer.xtal.dwba, fitted_rod, fitted_rod.angles
    )
    expected = CTRresolution.fast_convolve_intensity(
        fitted_rod.harr,
        fitted_rod.karr,
        fitted_rod.l,
        intensity,
        resolution,
        fitted_rod.angles,
    )
    np.testing.assert_allclose(optimizer.flat_prediction(), expected)


def test_sampling_generates_displaced_geometry_and_broadens_intensity():
    """Sampling evaluates displaced-L geometry through the CTR scan rule."""
    crystal = _crystal()
    geometry = CTRplotutil.CTRScanGeometry("in", angle=0.04)
    rod = _ctr(
        crystal,
        hk=(0.2, 0.0),
        geometry=geometry,
        records=False,
        l_values=[0.6, 0.8, 1.0],
    )
    optimizer = CTRopt.CTROptimizer(
        crystal, CTRplotutil.CTRCollection([rod]), scale_policy={"R": "fixed"}
    )
    resolution = CTRresolution.BoxResolution(0.03)
    optimizer.set_dwba()
    optimizer.set_resolution(resolution, calculation="sample")
    optimizer.prepareFit()

    fitted_rod = optimizer.CTRs[0]

    def intensity(h, k, l):  # noqa: E741
        angles = _generated_angles(
            optimizer.xtal, fitted_rod, h=h, k=k, l_values=l
        )
        values, _ = _direct_intensity_and_prefactor(
            optimizer.xtal.dwba, fitted_rod, angles
        )
        return values

    central = _generated_angles(optimizer.xtal, fitted_rod)
    expected = CTRresolution.sample_intensity(
        fitted_rod.harr,
        fitted_rod.karr,
        fitted_rod.l,
        intensity,
        resolution,
        central,
    )
    np.testing.assert_allclose(optimizer.flat_prediction(), expected, rtol=2e-13)


def test_sampling_geometry_accepts_periodic_records_and_rejects_other_branch():
    """Central validation tolerates 2-pi rotations but rejects mirror branches."""
    crystal = _crystal()
    geometry = CTRplotutil.CTRScanGeometry("in", angle=0.04)
    periodic = _ctr(
        crystal,
        hk=(0.2, 0.0),
        geometry=geometry,
        l_values=[0.6, 0.8, 1.0],
    )
    periodic.angles["omega"] += 2.0 * np.pi
    accepted = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([periodic]),
        scale_policy={"R": "fixed"},
    )
    accepted.set_dwba()
    accepted.set_resolution(
        CTRresolution.GaussianResolution(0.01), calculation="sample"
    )
    accepted.prepareFit()

    wrong_branch = _ctr(
        crystal,
        hk=(0.2, 0.0),
        geometry=geometry,
        l_values=[0.6, 0.8, 1.0],
    )
    wrong_branch.scan_geometry = CTRplotutil.CTRScanGeometry(
        "in", angle=0.04, mirrorx=True
    )
    rejected = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([wrong_branch]),
        scale_policy={"R": "fixed"},
    )
    rejected.set_dwba()
    rejected.set_resolution(
        CTRresolution.GaussianResolution(0.01), calculation="sample"
    )
    with pytest.raises(ValueError, match="central measured geometry|branch"):
        rejected.prepareFit()


def test_central_angle_records_must_reconstruct_the_measured_hkl():
    """Valid Vlieg records for another coordinate cannot enter the fit."""
    crystal = _crystal()
    rod = _ctr(
        crystal,
        hk=(0.2, 0.0),
        geometry=CTRplotutil.CTRScanGeometry("in", angle=0.04),
        l_values=[0.6, 0.8, 1.0],
    )
    rod.angles["omega"] += 0.02
    optimizer = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([rod]),
        scale_policy={"R": "fixed"},
    )
    optimizer.set_dwba()

    with pytest.raises(ValueError, match="do not reproduce every measured H, K"):
        optimizer.prepareFit()


def test_dwba_preparation_rejects_missing_metadata_and_unsupported_sampling():
    """Preparation reports missing reduction/geometry and deferred sampling."""
    crystal = _crystal()
    legacy = CTRplotutil.CTR(
        (0.0, 0.0), [0.2], [1.0], [0.1], scan_geometry=CTRplotutil.CTRScanGeometry("eq")
    )
    missing_polarization = CTRopt.CTROptimizer(
        crystal, CTRplotutil.CTRCollection([legacy])
    )
    missing_polarization.set_dwba()
    with pytest.raises(ValueError, match="requires polarization-reduction"):
        missing_polarization.prepareFit()

    no_geometry = _ctr(crystal, records=False)
    no_geometry.scan_geometry = None
    missing_geometry = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([no_geometry]),
        scale_policy={"R": "fixed"},
    )
    missing_geometry.set_dwba()
    with pytest.raises(ValueError, match="central angle records or CTRScanGeometry"):
        missing_geometry.prepareFit()

    records_only = _ctr(crystal, records=True, keep_rule=False)
    sampling = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([records_only]),
        scale_policy={"R": "fixed"},
    )
    sampling.set_dwba()
    sampling.set_resolution(
        CTRresolution.GaussianResolution(0.01), calculation="sample"
    )
    with pytest.raises(ValueError, match="use calculation='convolve'"):
        sampling.prepareFit()

    fitted_sampling = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([_ctr(crystal, records=False)]),
        scale_policy={"R": "fixed"},
    )
    fitted_sampling.set_dwba()
    fitted_sampling.fit_resolution(
        CTRresolution.GaussianResolution(0.01),
        [0.0, 0.0, 0.0],
        [0.1, 0.1, 0.1],
        calculation="sample",
    )
    with pytest.raises(ValueError, match="fitted-width resolution sampling"):
        fitted_sampling.prepareFit()


@pytest.mark.parametrize(
    ("attribute", "values", "message"),
    [
        ("sfI", [1.0, np.nan, -1.0], "observations must be finite real"),
        (
            "err",
            [0.1, 0.0, 0.3],
            "uncertainties must be finite and strictly positive",
        ),
        ("err", [0.1, 0.2], "must be one-dimensional and point-aligned"),
    ],
)
def test_dwba_preparation_rejects_malformed_measurement_arrays(
    attribute, values, message
):
    """Malformed arrays fail explicitly without filtering signed values."""
    crystal = _crystal()
    rod = _ctr(crystal, values=[-1.0, 0.0, 2.0])
    setattr(rod, attribute, np.asarray(values, dtype=np.float64))
    optimizer = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([rod]),
        scale_policy={"R": "fixed"},
    )
    optimizer.set_dwba()

    with pytest.raises(ValueError, match=message):
        optimizer.prepareFit()


def test_live_dwba_evaluation_is_deterministic_for_a_b_a_parameters():
    """Live parameter evaluation returns exactly to its first prediction."""
    crystal = _crystal(fit_x=True)
    rod = _ctr(
        crystal,
        hk=(0.2, 0.0),
        geometry=CTRplotutil.CTRScanGeometry("in", angle=0.04),
        records=False,
        l_values=[0.6, 0.8, 1.0],
    )
    optimizer = CTRopt.CTROptimizer(
        crystal, CTRplotutil.CTRCollection([rod]), scale_policy={"R": "fixed"}
    )

    optimizer.set_dwba()
    optimizer.prepareFit()
    first_a = optimizer.flat_prediction([0.17]).copy()
    prediction_b = optimizer.flat_prediction([0.34]).copy()
    second_a = optimizer.flat_prediction([0.17]).copy()

    assert not np.array_equal(first_a, prediction_b)
    np.testing.assert_array_equal(first_a, second_a)

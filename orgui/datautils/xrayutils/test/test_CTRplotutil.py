"""Regression coverage for CTR persistence and measurement metadata."""

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest
from silx.io.dictdump import dicttonx, nxtodict

from ... import util as data_util
from .. import CTRplotutil
from ..CTRplotutil import (
    CTR,
    CTRCollection,
    CTRScanGeometry,
    MeasurementReduction,
    PolarizationReduction,
)


_ANGLE_NAMES = ("alpha", "delta", "gamma", "omega", "chi", "phi")


def _rod(npoints):
    rod = CTR((1.0, -2.0), np.linspace(0.1, 0.9, npoints),
              np.arange(npoints) + 1.0, err=np.full(npoints, 0.1))
    rod.angles = np.rec.fromarrays(
        [np.linspace(0.01, 0.02, npoints) + i * 0.1 for i in range(6)],
        names=_ANGLE_NAMES,
    )
    rod.bgI = np.arange(npoints) * 0.2
    return rod


@pytest.mark.parametrize("npoints", [1, 2, 6, 9])
def test_nexus_round_trip_preserves_all_angle_columns_and_k(npoints):
    """Point counts must not be confused with the six angle columns."""
    original = _rod(npoints)
    payload = original.toNXdict()
    assert payload["@orgui_ctr_schema"] == 2
    assert payload["sixc_angles"]["@unit"] == "rad"
    restored = CTR.fromNXdict(payload)
    assert restored.hk == original.hk
    for attr in ("harr", "karr", "l", "sfI", "err", "bgI"):
        np.testing.assert_array_equal(getattr(restored, attr), getattr(original, attr))
    for name in _ANGLE_NAMES:
        np.testing.assert_array_equal(restored.angles[name], original.angles[name])


def test_nexus_file_round_trip_preserves_schema_and_angles(tmp_path):
    """The schema marker and radian unit survive the actual NeXus writer."""
    original = _rod(2)
    path = str(tmp_path / "ctr_angles.nxs")
    dicttonx({"entry": CTRCollection([original], name="test").toNXdict()}, path)
    payload = nxtodict(path)["entry"]
    restored = CTRCollection.fromNXdict(payload)[0]
    np.testing.assert_array_equal(restored.karr, original.karr)
    for name in _ANGLE_NAMES:
        np.testing.assert_array_equal(restored.angles[name], original.angles[name])


@pytest.mark.parametrize("label", ["deg", "rad", None])
def test_unversioned_nexus_preserves_legacy_radian_values(label):
    """Old orGUI's degree label does not trigger an implicit conversion."""
    original = _rod(2)
    payload = original.toNXdict()
    del payload["@orgui_ctr_schema"]
    if label is None:
        del payload["sixc_angles"]["@unit"]
    else:
        payload["sixc_angles"]["@unit"] = label
    restored = CTR.fromNXdict(payload)
    for name in _ANGLE_NAMES:
        np.testing.assert_array_equal(restored.angles[name], original.angles[name])


@pytest.mark.parametrize("versioned", [False, True])
def test_external_degree_conversion_and_collection_forwarding(versioned):
    """Explicit input units override metadata and reach every collection rod."""
    original = _rod(2)
    payload = CTRCollection([original], name="test").toNXdict()
    rod_payload = payload[repr(original)]
    if not versioned:
        del rod_payload["@orgui_ctr_schema"]
    for name in _ANGLE_NAMES:
        rod_payload["sixc_angles"][name] = np.rad2deg(original.angles[name])
    # The explicit override takes priority even over a new writer's rad label.
    restored = CTRCollection.fromNXdict(payload, angle_units="deg")[0]
    for name in _ANGLE_NAMES:
        np.testing.assert_allclose(restored.angles[name], original.angles[name])
    if versioned:
        rod_payload["sixc_angles"]["@unit"] = "deg"
        automatic = CTR.fromNXdict(rod_payload)
        for name in _ANGLE_NAMES:
            np.testing.assert_allclose(automatic.angles[name], original.angles[name])


def test_legacy_k_values_are_not_guessed():
    """A historical corrupt K cannot be distinguished from a real H=K rod."""
    payload = _rod(2).toNXdict()
    del payload["@orgui_ctr_schema"]
    payload["hkl"]["k"] = payload["hkl"]["h"].copy()
    restored = CTR.fromNXdict(payload)
    np.testing.assert_array_equal(restored.karr, payload["hkl"]["k"])


def test_nexus_rejects_bad_units_and_mismatched_angle_lengths():
    """Invalid angle metadata fails before creating misleading records."""
    payload = _rod(2).toNXdict()
    with pytest.raises(ValueError, match="angle_units"):
        CTR.fromNXdict(payload, angle_units="turns")
    payload["sixc_angles"]["alpha"] = [0.1]
    with pytest.raises(ValueError, match="same shape as L"):
        CTR.fromNXdict(payload)


def _polarized_reduction(factor=(0.8, 0.9, 1.0)):
    """Return structure-factor metadata with pointwise conventional P."""
    return MeasurementReduction(
        "structure_factor",
        PolarizationReduction(0.7, "unanalysed", factor),
    )


def test_convert_to_f_retains_signed_values_and_transforms_interval_errors():
    """Signed conversion is reversible and finite through an interval crossing."""
    intensity = np.array([-9.0, 0.0, 1.0, 9.0])
    rod = CTR(
        (1.0, 0.0),
        [0.1, 0.2, 0.3, 0.4],
        intensity,
        [1.0, 4.0, 2.0, 1.0],
        reduction=_polarized_reduction([0.7, 0.8, 0.9, 1.0]),
        scan_geometry=CTRScanGeometry("in", 0.04),
    )
    rod.bgI = np.array([-4.0, 0.0, 9.0, 16.0])
    rod.ctrI = np.array([-1.0, 1.0, 4.0, 25.0])

    rod.convertToF()

    np.testing.assert_array_equal(rod.sfI, [-3.0, 0.0, 1.0, 3.0])
    np.testing.assert_allclose(rod.sfI * np.abs(rod.sfI), intensity)
    np.testing.assert_allclose(
        rod.err,
        [
            (np.sqrt(10.0) - np.sqrt(8.0)) / 2.0,
            2.0,
            (np.sqrt(3.0) + 1.0) / 2.0,
            (np.sqrt(10.0) - np.sqrt(8.0)) / 2.0,
        ],
    )
    np.testing.assert_array_equal(rod.bgI, [-2.0, 0.0, 3.0, 4.0])
    np.testing.assert_array_equal(rod.ctrI, [-1.0, 1.0, 2.0, 5.0])
    assert rod.reduction == _polarized_reduction([0.7, 0.8, 0.9, 1.0])
    assert rod.scan_geometry == CTRScanGeometry("in", 0.04)


def test_convert_to_f_uses_stable_high_signal_uncertainty():
    """Strong-signal propagation must not subtract nearly equal roots."""
    maximum = np.finfo(np.float64).max
    intensity = np.array([-maximum, maximum])
    sigma = np.array([1e300, 1e300])
    rod = CTR((1.0, 0.0), [0.1, 0.2], intensity, sigma)

    rod.convertToF()

    expected_linear = sigma / (2.0 * np.sqrt(np.abs(intensity)))
    np.testing.assert_allclose(rod.err, expected_linear, rtol=1e-15)
    assert rod.err[0] == rod.err[1]
    assert np.all(rod.err > 0.0)


def test_convert_to_f_excludes_only_nonfinite_points_and_keeps_alignment():
    """Default filtering retains negative/zero data and selects all metadata."""
    reduction = _polarized_reduction([0.7, 0.8, 0.9, 1.0])
    rod = CTR(
        (1.0, 0.0),
        [0.1, 0.2, 0.3, 0.4],
        [-4.0, np.nan, 0.0, np.inf],
        reduction=reduction,
    )

    rod.convertToF()

    np.testing.assert_array_equal(rod.l, [0.1, 0.3])
    np.testing.assert_array_equal(rod.sfI, [-2.0, 0.0])
    np.testing.assert_array_equal(
        rod.reduction.polarization.polarization_factor, [0.7, 0.9]
    )


@pytest.mark.parametrize("error", ([1.0, 0.0], [1.0, np.nan], [1.0]))
def test_convert_to_f_rejects_invalid_uncertainty_without_mutating(error):
    """Invalid interval inputs fail atomically instead of inventing errors."""
    rod = CTR((1.0, 0.0), [0.1, 0.2], [-4.0, 9.0], error)
    original_values = rod.sfI.copy()
    original_errors = rod.err.copy()

    with pytest.raises(ValueError, match="uncert"):
        rod.convertToF()

    np.testing.assert_array_equal(rod.sfI, original_values)
    np.testing.assert_array_equal(rod.err, original_errors)


def test_collection_convert_to_f_forwards_filtering_and_rejects_reflectivity():
    """The collection forwards compatibility options to every member."""
    retained_nan = CTR((1.0, 0.0), [0.1, 0.2], [-1.0, np.nan])
    zero = CTR((2.0, 0.0), [0.1], [0.0])
    collection = CTRCollection([retained_nan, zero])

    collection.convertToF(excludeInvalid=False)

    np.testing.assert_array_equal(retained_nan.sfI[:1], [-1.0])
    assert np.isnan(retained_nan.sfI[1])
    np.testing.assert_array_equal(zero.sfI, [0.0])

    reflectivity = CTR(
        (3.0, 0.0),
        [0.1],
        [0.5],
        reduction=MeasurementReduction("reflectivity"),
    )
    with pytest.raises(ValueError, match="structure-factor data only"):
        reflectivity.convertToF()


def test_measurement_metadata_is_validated_immutable_and_value_comparable():
    """Array metadata is copied, read-only, comparable, and unhashable."""
    source = np.array([0.8, 0.9, 1.0])
    polarization = PolarizationReduction(0.7, "unanalysed", source)
    reduction = MeasurementReduction("structure_factor", polarization)
    equivalent = _polarized_reduction()

    source[0] = 99.0
    assert reduction == equivalent
    assert not reduction.polarization.polarization_factor.flags.writeable
    assert reduction != MeasurementReduction(
        "structure_factor", PolarizationReduction(0.6, "unanalysed", [0.8, 0.9, 1.0])
    )
    with pytest.raises(TypeError):
        hash(polarization)
    with pytest.raises(TypeError):
        hash(reduction)

    geometry = CTRScanGeometry("in", 0.05, mirrorx=True)
    assert geometry == CTRScanGeometry("in", 0.05, mirrorx=True)
    assert isinstance(hash(geometry), int)


@pytest.mark.parametrize(
    "factory, message",
    [
        (lambda: PolarizationReduction(-0.1, "s"), "s_fraction"),
        (lambda: PolarizationReduction(True, "s"), "s_fraction"),
        (lambda: PolarizationReduction(0.5, "circular"), "outgoing"),
        (lambda: PolarizationReduction(0.5, "p", [1.0, 0.0]), "strictly positive"),
        (lambda: MeasurementReduction("intensity"), "quantity"),
        (lambda: MeasurementReduction("structure_factor", object()), "polarization"),
        (lambda: CTRScanGeometry("in"), "angle"),
        (lambda: CTRScanGeometry("out", 0.0), "angle"),
        (lambda: CTRScanGeometry("out", True), "angle"),
        (lambda: CTRScanGeometry("eq", 0.1), "None"),
        (lambda: CTRScanGeometry("bad", None), "fixed"),
    ],
)
def test_metadata_records_reject_invalid_values(factory, message):
    """Invalid reduction and scan-rule records fail at construction."""
    with pytest.raises((TypeError, ValueError), match=message):
        factory()


def test_ctr_metadata_properties_validate_shape_quantity_and_reset():
    """CTR assignment applies shape and quantity checks and names the rod."""
    rod = _rod(3)
    assert rod.reduction == MeasurementReduction()
    assert rod.scan_geometry is None
    with pytest.raises(TypeError, match="<CTR.*reduction"):
        rod.reduction = None
    with pytest.raises(TypeError, match="<CTR.*scan_geometry"):
        rod.scan_geometry = "fixed-in"
    with pytest.raises(ValueError, match="<CTR.*same shape"):
        rod.reduction = _polarized_reduction([0.8, 0.9])
    with pytest.raises(ValueError, match="<CTR.*already be corrected"):
        rod.reduction = MeasurementReduction(
            "reflectivity", PolarizationReduction(1.0, "s", [1.0, 1.0, 1.0])
        )

    rod.reduction = _polarized_reduction()
    rod.scan_geometry = CTRScanGeometry("out", 0.04)
    rod.scan_geometry = None
    rod.reduction = MeasurementReduction()
    assert rod.scan_geometry is None
    assert rod.reduction == MeasurementReduction()


def test_nexus_round_trip_preserves_reduction_factor_and_scan_rule():
    """Schema 2 stores conventional P and the scalar scan rule in radians."""
    original = CTR(
        (1.0, -2.0),
        [0.1, 0.2, 0.3],
        [1.0, 2.0, 3.0],
        [0.1, 0.1, 0.1],
        reduction=_polarized_reduction(),
        scan_geometry=CTRScanGeometry("in", 0.03, mirrorx=True),
    )
    payload = original.toNXdict()
    stored_factor = payload["measurement_reduction"]["polarization"][
        "polarization_factor"
    ]
    np.testing.assert_array_equal(stored_factor, [0.8, 0.9, 1.0])

    restored = CTR.fromNXdict(payload)
    assert restored.reduction == original.reduction
    assert restored.scan_geometry == original.scan_geometry


def test_nexus_file_round_trip_preserves_measurement_metadata(tmp_path):
    """Metadata survives the actual NeXus writer rather than only dictionaries."""
    original = CTR(
        (1.0, 0.0),
        [0.1, 0.2, 0.3],
        [1.0, 2.0, 3.0],
        reduction=_polarized_reduction(),
        scan_geometry=CTRScanGeometry("eq"),
    )
    path = str(tmp_path / "ctr_metadata.nxs")
    dicttonx({"entry": CTRCollection([original], name="test").toNXdict()}, path)
    restored = CTRCollection.fromNXdict(nxtodict(path)["entry"])[0]
    assert restored.reduction == original.reduction
    assert restored.scan_geometry == original.scan_geometry


def test_existing_schema_2_without_fitting_metadata_uses_legacy_defaults():
    """Early schema-2 payloads remain valid without the new metadata groups."""
    payload = _rod(2).toNXdict()
    del payload["measurement_reduction"]
    restored = CTR.fromNXdict(payload)
    assert restored.reduction == MeasurementReduction()
    assert restored.scan_geometry is None


def test_cut_and_deepcopy_preserve_aligned_metadata_without_aliasing():
    """Selections slice factors and angle records with the measured points."""
    original = _rod(3)
    original.reduction = _polarized_reduction()
    original.scan_geometry = CTRScanGeometry("in", 0.05)
    copied = copy.deepcopy(original)
    copied.cut(1, 3)

    np.testing.assert_array_equal(copied.l, original.l[1:3])
    np.testing.assert_array_equal(
        copied.angles["gamma"], original.angles["gamma"][1:3]
    )
    np.testing.assert_array_equal(
        copied.reduction.polarization.polarization_factor, [0.9, 1.0]
    )
    assert copied.scan_geometry == original.scan_geometry
    assert copied.reduction.polarization.polarization_factor is not (
        original.reduction.polarization.polarization_factor
    )


def test_plain_array_and_anarod_imports_accept_measurement_metadata():
    """Plain import overrides are propagated and full-table P is split by rod."""
    geometry = CTRScanGeometry("out", 0.04)
    one_rod = np.array(
        [[1.0, 0.0, 0.1, 2.0, 0.2], [1.0, 0.0, 0.2, 3.0, 0.3]]
    )
    reduction = MeasurementReduction(
        "structure_factor", PolarizationReduction(1.0, "s", [0.8, 0.9])
    )
    restored = CTR.fromArray(
        one_rod, reduction=reduction, scan_geometry=geometry
    )
    assert restored.reduction == reduction
    assert restored.scan_geometry == geometry

    two_rods = np.vstack((one_rod, one_rod + [1.0, 1.0, 0.0, 0.0, 0.0]))
    table_reduction = MeasurementReduction(
        "structure_factor",
        PolarizationReduction(1.0, "s", [0.8, 0.9, 1.0, 1.1]),
    )
    collection = CTRCollection.fromANAROD(
        two_rods, reduction=table_reduction, scan_geometry=geometry
    )
    assert all(rod.scan_geometry == geometry for rod in collection)
    np.testing.assert_array_equal(
        collection[0].reduction.polarization.polarization_factor, [0.8, 0.9]
    )
    np.testing.assert_array_equal(
        collection[1].reduction.polarization.polarization_factor, [1.0, 1.1]
    )


def test_ctr_file_loader_reads_metadata_and_splits_pointwise_factor(tmp_path):
    """The compact CTR format restores reduction metadata in one call."""
    path = tmp_path / "measurement.ctr"
    path.write_text(
        "# orgui_ctr_schema: 2\n"
        "# quantity: structure_factor\n"
        "# s_fraction: 0.75\n"
        "# outgoing: unanalysed\n"
        "# columns: H K L F_HKL errorF polarization_factor\n"
        "1 0 0.1 10 1 0.8\n"
        "1 0 0.2 11 1.1 0.9\n"
        "2 0 0.1 20 2 1.0\n",
        encoding="utf-8",
    )

    collection = CTRCollection.fromCTRFile(path)

    assert collection.name == path.name
    assert collection.getHKList() == [(1.0, 0.0), (2.0, 0.0)]
    np.testing.assert_array_equal(collection[0].sfI, [10.0, 11.0])
    np.testing.assert_array_equal(collection[0].err, [1.0, 1.1])
    assert collection[0].reduction.quantity == "structure_factor"
    assert collection[0].reduction.polarization.s_fraction == 0.75
    assert collection[0].reduction.polarization.outgoing == "unanalysed"
    np.testing.assert_array_equal(
        collection[0].reduction.polarization.polarization_factor, [0.8, 0.9]
    )
    np.testing.assert_array_equal(
        collection[1].reduction.polarization.polarization_factor, [1.0]
    )


def test_ctr_file_loader_rejects_legacy_anarod_file(tmp_path):
    """The new loader must not silently guess metadata for legacy files."""
    path = tmp_path / "legacy.dat"
    path.write_text(
        "# H K L F_HKL errorF mode\n1 0 0.1 10 1 -3\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="fromANAROD"):
        CTRCollection.fromCTRFile(path)


def test_ctr_file_loader_requires_trailing_polarization_factor(tmp_path):
    """A pointwise extension cannot displace value or uncertainty columns."""
    path = tmp_path / "reordered.ctr"
    path.write_text(
        "# orgui_ctr_schema: 2\n"
        "# quantity: structure_factor\n"
        "# columns: H K L F_HKL polarization_factor errorF\n"
        "1 0 0.1 10 0.8 1\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="optional trailing"):
        CTRCollection.fromCTRFile(path)


def test_metadata_dependent_binning_and_reflectivity_consumers_reject(tmp_path):
    """Legacy operations fail instead of inventing reduction semantics."""
    legacy = CTR((1.0, 0.0), [0.1, 0.2], [1.0, 2.0], [0.1, 0.1])
    assert legacy.generateAverage(nbins=1).reduction == MeasurementReduction()

    metadata_rod = copy.deepcopy(legacy)
    metadata_rod.reduction = MeasurementReduction(
        "structure_factor", PolarizationReduction(1.0, "s")
    )
    with pytest.raises(NotImplementedError, match="measurement-reduction"):
        metadata_rod.generateAverage(nbins=1)

    reflectivity = CTR(
        (1.0, 0.0),
        [0.1, 0.2],
        [1.0, 2.0],
        [0.1, 0.1],
        reduction=MeasurementReduction(
            "reflectivity", PolarizationReduction(1.0, "s")
        ),
    )
    with pytest.raises(ValueError, match="structure-factor data only"):
        reflectivity.get_scale(None)
    collection = CTRCollection([reflectivity])
    with pytest.raises(ValueError, match="collection scaling"):
        collection *= 2.0
    with pytest.raises(ValueError, match="CTR differences"):
        reflectivity.generateDifference(legacy)
    with pytest.raises(ValueError, match="symmetry averaging"):
        data_util.averageCTRs([[reflectivity]])
    with pytest.raises(ValueError, match="ANAROD structure-factor export"):
        collection.toANAROD(tmp_path / "reflectivity.dat")


def test_reflectivity_rejects_phase_and_complex_structure_factor_operations():
    """A reflectivity tag cannot imply a complex structure factor."""
    reduction = MeasurementReduction(
        "reflectivity", PolarizationReduction(1.0, "s")
    )
    with pytest.raises(ValueError, match="phase information"):
        CTR((1.0, 0.0), [0.1], [0.5], phi=[0.0], reduction=reduction)

    reflectivity = CTR((1.0, 0.0), [0.1], [0.5], reduction=reduction)
    with pytest.raises(ValueError, match="phase assignment"):
        reflectivity.setPhase([0.0])
    with pytest.raises(ValueError, match="complex structure factors"):
        reflectivity.getComplexSF()


def test_angle_generation_rejects_a_non_round_tripping_solver():
    """The shared helper validates every reconstructed reference coordinate."""

    class NonRoundTrippingAngles:
        @staticmethod
        def anglesZmode(hkl, fixedangle, **kwargs):
            return np.zeros((hkl.shape[1], 6))

        @staticmethod
        def anglesToHkl(alpha, delta, gamma, omega, chi, phi):
            zeros = np.zeros_like(alpha)
            return zeros, zeros, zeros

    rod = CTR((1.0, 0.0), [0.1, 0.2], [1.0, 1.0])
    with pytest.raises(ValueError, match="do not reproduce every requested"):
        rod.calcAnglesZmode(NonRoundTrippingAngles(), fixedangle=0.05)


def test_quantity_aware_plot_labels_and_axis_separation():
    """Mixed F/R panels use distinct labels and never share their y axis."""
    structure_factor = CTR((1.0, 0.0), [0.1], [1.0])
    reflectivity = CTR(
        (2.0, 0.0),
        [0.1],
        [1.0],
        reduction=MeasurementReduction(
            "reflectivity", PolarizationReduction(1.0, "s")
        ),
    )
    figure = CTRplotutil.ctrfigure()
    figure.addCollection(CTRCollection([structure_factor, reflectivity]))
    figure.generateCTRplot(cols=2)
    assert figure.axes[0].get_ylabel() == "Structure factor / arb. units"
    assert figure.axes[1].get_ylabel() == "Reflectivity / dimensionless"
    assert not figure.axes[0].get_shared_y_axes().joined(
        figure.axes[0], figure.axes[1]
    )
    plt.close(figure)

    structure_factor.setToDefaultID()
    reflectivity.hk = structure_factor.hk
    reflectivity.setToDefaultID()
    figure = CTRplotutil.ctrfigure()
    figure.addCTR(structure_factor)
    with pytest.raises(ValueError, match="different quantities"):
        figure.addCTR(reflectivity)
    plt.close(figure)

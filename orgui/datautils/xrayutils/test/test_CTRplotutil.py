"""Regression coverage for CTR NeXus persistence and angle metadata."""

import numpy as np
import pytest
from silx.io.dictdump import dicttonx, nxtodict

from ..CTRplotutil import CTR, CTRCollection


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

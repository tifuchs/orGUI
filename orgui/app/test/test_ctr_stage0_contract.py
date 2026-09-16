"""Frozen CTR extraction behavior before the targeted correction rework.

The compact JSON fixture records the numerical contract at commit
``d52ab8f5234e``.  These are characterization tests: they protect legacy
dispatch and make later convention changes explicit, but they are not an
independent validation of the CTR physics.
"""

import json
from pathlib import Path

import numpy as np

from orgui.app import integration_corrections as stationary
from orgui.app.config_data import CorrectionState
from orgui.app.peak1Dintegr import _compute_rocking_integration


_FIXTURE = Path(__file__).with_name("fixtures") / "ctr_stage0_contract.json"


def _contract():
    """Return the immutable Stage-0 characterization data."""
    return json.loads(_FIXTURE.read_text(encoding="utf-8"))


def test_stationary_contract_at_stage0():
    """Pin stationary normalization, Lorentz reduction, errors, and units."""
    case = _contract()["stationary"]
    alpha, delta, gamma = (
        np.asarray(case[name], dtype=float)
        for name in ("alpha", "delta", "gamma")
    )
    factors = stationary.stationary_correction_factors(
        alpha,
        delta,
        gamma,
        use_lorentz=True,
        normalization=np.asarray(case["normalization"], dtype=float),
    )
    intensity, errors = stationary.apply_stationary_corrections(
        case["counts"], case["errors"], factors
    )
    f2, f2_errors = stationary.structure_factor(intensity, errors, factors)

    expected = case["expected"]
    np.testing.assert_allclose(intensity, expected["intensity"], rtol=1e-14)
    np.testing.assert_allclose(
        errors, expected["intensity_errors"], rtol=1e-14
    )
    np.testing.assert_allclose(f2, expected["F2_hkl"], rtol=1e-14)
    np.testing.assert_allclose(
        f2_errors, expected["F2_hkl_errors"], rtol=1e-14
    )


def test_rocking_contract_at_stage0():
    """Pin two-background aggregation and current rocking F2 scaling."""
    case = _contract()["rocking"]
    curves = np.asarray(case["curves"], dtype=float)
    shape = curves.shape
    roi = {
        name: {
            "from": np.full(shape[0], limits[0], dtype=float),
            "to": np.full(shape[0], limits[1], dtype=float),
        }
        for name, limits in case["roi"].items()
    }
    result = _compute_rocking_integration(
        np.asarray(case["s"], dtype=float),
        np.asarray(case["axis"], dtype=float),
        curves,
        np.sqrt(curves),
        roi,
        aux={},
        use_lorentz=True,
        use_footprint=False,
        C_Lor=np.broadcast_to(
            np.asarray(case["C_Lorentz"], dtype=float)[:, None], shape
        ),
        C_rod=np.broadcast_to(
            np.asarray(case["C_rod"], dtype=float)[:, None], shape
        ),
        C_norm=np.broadcast_to(
            np.asarray(case["normalization"], dtype=float)[:, None], shape
        ),
        detector_acceptance=np.asarray(case["acceptance"], dtype=float),
        angle_unit=case["units"]["axis"],
    )

    for name, expected in case["expected"].items():
        np.testing.assert_allclose(result[name], expected, rtol=1e-14)


def test_legacy_correction_settings_keep_their_old_meanings():
    """Pin the opaque-JSON settings used by pre-typed database groups."""
    values = _contract()["legacy_corrections"]

    state = CorrectionState.from_dict(values)

    serialized = state.to_dict()
    assert {name: serialized[name] for name in values} == values
    assert state.use_background is True
    assert state.monitor_corrections == ("mondio",)
    # These fields retain the old density/area convention until a versioned
    # replacement exists: meter, meter, and photons / (s m^2), respectively.
    assert state.sample_length_m == 0.004
    assert state.sample_width_m == 0.008
    assert state.beam_flux_density == 2.5e12

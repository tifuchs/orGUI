"""Audit native optical-profile invariance under interface-support expansion.

This diagnostic verifies one support-expansion fixture without changing any
tolerance. Run with python -m benchmarks.audit_dissolution_dwba_partition.
The synthetic fixtures require no measured data.
"""

import argparse
import copy
import json
from unittest.mock import patch

import numpy as np

from orgui.datautils.xrayutils import CTRdwba
from orgui.datautils.xrayutils.CTRdistributions import (
    RecenteredSkellamProfile, SkellamProfile,
)
from orgui.datautils.xrayutils.CTRincoherent import FiniteFilmDissolution
from orgui.datautils.xrayutils.test.test_CTRdissolution import _crystal, _hkl


def _errors(actual, reference):
    return dict(max_relative=float(np.max(abs(actual/reference-1))),
                max_scale_normalized=float(np.max(abs(actual-reference))
                                           / np.max(abs(reference))))


def _crystals(offset=.03):
    source = _crystal(width=40, interface=True)
    source["interface"].basis[3] = offset
    source["interface"].basis_0[3] = offset
    source["interface"].profile = SkellamProfile(.1, tail_probability=1e-14)
    source["interface"].createInterfaceCells()
    surface = source["surface"]
    surface.unitcell.basis[:, 6] *= .5
    surface.unitcell.basis_0[:] = surface.unitcell.basis
    surface.set_ucs(surface.unitcell)
    source.apply_stacking()
    expanded = copy.deepcopy(source)
    expanded["interface"].fixed_ucs = 20
    expanded["interface"].createInterfaceCells()
    expanded.apply_stacking()
    return source, expanded


def _report(source, expanded):
    report = dict(
        interpretation="Synthetic support audit; not a universal DWBA error bound.",
        native_kinematic=_errors(expanded.F(*_hkl()), source.F(*_hkl())),
        native_dwba=_errors(expanded.dwba.evaluate(*_hkl()).reflectivity,
                            source.dwba.evaluate(*_hkl()).reflectivity))
    for policy in ("shared", "state_specific"):
        models = [FiniteFilmDissolution(
            c, films=["film"], surface="surface", optical_reference=policy,
            profile=RecenteredSkellamProfile(offset=-18.5),
            incoherent_fraction=.2) for c in (source, expanded)]
        report[policy] = _errors(models[1].dwba.evaluate(*_hkl()).reflectivity,
                                  models[0].dwba.evaluate(*_hkl()).reflectivity)
    return report


def _legacy_anchor_probe(original):
    """Restore the old lowest-input coordinate to reproduce the fixed defect.

    This is a causal experiment, not a production profile algorithm. It keeps
    all density/reference shares, all empty groups, and the 0.6-Angstrom grouping
    tolerance. Production now uses the lowest nonzero contributor instead.
    """
    def combine(values, shares, z_tolerance=.6, density_count=None):
        result, tagged = original(values, shares, z_tolerance, density_count)
        if not len(values):
            return result, tagged
        order = np.argsort(values[:, 0], kind="stable")
        values = values[order]
        starts = [0]
        for index in range(1, len(values)):
            if values[index, 0]-values[starts[-1], 0] > z_tolerance:
                starts.append(index)
        for group, start in enumerate(starts):
            result[group, 0] = values[start, 0]
        return result, tagged
    return combine


def main():
    """Audit native states, optionally isolate zero-weight position selection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnose-zero-anchors", action="store_true")
    args = parser.parse_args()
    source, expanded = _crystals()
    report = _report(source, expanded)
    for name in ("native_dwba", "shared", "state_specific"):
        assert report[name]["max_relative"] < 1e-9
    if args.diagnose_zero_anchors:
        with patch.object(CTRdwba, "_combine_tagged", _legacy_anchor_probe(
                CTRdwba._combine_tagged)):
            legacy = _report(*_crystals())
            narrow = source.dwba._atomic_model_snapshot(raw_structure=True)
            wide = expanded.dwba._atomic_model_snapshot(raw_structure=True)
            with patch.object(CTRdwba, "_OPTICAL_DELTA_TOLERANCE", 0.), \
                    patch.object(CTRdwba, "_OPTICAL_BETA_TOLERANCE", 0.):
                no_simplification = _report(*_crystals())
        np.testing.assert_equal(narrow.shape, wide.shape)
        np.testing.assert_allclose(narrow[:, 1:], wide[:, 1:],
                                   rtol=1e-12, atol=1e-18)
        changed = np.flatnonzero(abs(narrow[:, 0]-wide[:, 0]) > 1e-9)
        for name in ("native_dwba", "shared", "state_specific"):
            np.testing.assert_allclose(no_simplification[name]["max_relative"],
                                       legacy[name]["max_relative"], rtol=1e-10)
        report = dict(
            legacy_anchor=legacy, legacy_no_simplification=no_simplification,
            corrected=report,
            zero_interface_displacement=_report(*_crystals(offset=0.)),
            moved_centers_A=[dict(narrow=float(narrow[i, 0]),
                                   expanded=float(wide[i, 0])) for i in changed],
            max_raw_delta_beta_change=float(np.max(abs(narrow[:, 1:]-wide[:, 1:]))),
            probe_scope="Temporary Python patch only; native kernel and physical "
                        "state weights, coordinates and tolerances are unchanged.")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

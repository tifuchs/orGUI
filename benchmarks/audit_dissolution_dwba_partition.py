"""Report inherited native optical-profile sensitivity to interface support.

This diagnostic does not declare convergence or change any tolerance. Run from
the checkout with python -m benchmarks.audit_dissolution_dwba_partition.
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
        interpretation="Numerical diagnostic, not a passed DWBA support gate.",
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


def _active_anchor_probe(original):
    """Change only the representative position of a nonempty merged group.

    This is a causal experiment, not a production profile algorithm. It keeps
    all density/reference shares, all empty groups, and the 0.6-Angstrom grouping
    tolerance. Signed nonzero shares qualify even if their total density cancels.
    """
    def combine(values, shares, z_tolerance=.6, density_count=None):
        result, tagged = original(values, shares, z_tolerance, density_count)
        if not len(values):
            return result, tagged
        order = np.argsort(values[:, 0], kind="stable")
        values, shares = values[order], shares[:, order]
        starts = [0]
        for index in range(1, len(values)):
            if values[index, 0]-values[starts[-1], 0] > z_tolerance:
                starts.append(index)
        for group, (start, stop) in enumerate(zip(starts, [*starts[1:], len(values)])):
            active = np.flatnonzero(np.any(shares[:, start:stop] != 0, axis=(0, 2)))
            if len(active):
                result[group, 0] = values[start+active[0], 0]
        return result, tagged
    return combine


def main():
    """Audit native states, optionally isolate zero-weight position selection."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnose-zero-anchors", action="store_true")
    args = parser.parse_args()
    source, expanded = _crystals()
    report = _report(source, expanded)
    if args.diagnose_zero_anchors:
        narrow = source.dwba._atomic_model_snapshot(raw_structure=True)
        wide = expanded.dwba._atomic_model_snapshot(raw_structure=True)
        np.testing.assert_equal(narrow.shape, wide.shape)
        np.testing.assert_allclose(narrow[:, 1:], wide[:, 1:],
                                   rtol=1e-12, atol=1e-18)
        changed = np.flatnonzero(abs(narrow[:, 0]-wide[:, 0]) > 1e-9)
        with patch.object(CTRdwba, "_OPTICAL_DELTA_TOLERANCE", 0.), \
                patch.object(CTRdwba, "_OPTICAL_BETA_TOLERANCE", 0.):
            no_simplification = _report(*_crystals())
        with patch.object(CTRdwba, "_combine_tagged", _active_anchor_probe(
                CTRdwba._combine_tagged)):
            active_anchor = _report(*_crystals())
        for name in ("native_dwba", "shared", "state_specific"):
            assert active_anchor[name]["max_relative"] < 1e-9
            np.testing.assert_allclose(no_simplification[name]["max_relative"],
                                       report[name]["max_relative"], rtol=1e-10)
        report = dict(
            original=report, no_simplification=no_simplification,
            active_anchor_only=active_anchor,
            zero_interface_displacement=_report(*_crystals(offset=0.)),
            moved_centers_A=[dict(narrow=float(narrow[i, 0]),
                                   expanded=float(wide[i, 0])) for i in changed],
            max_raw_delta_beta_change=float(np.max(abs(narrow[:, 1:]-wide[:, 1:]))),
            probe_scope="Temporary Python patch only; native kernel and physical "
                        "state weights, coordinates and tolerances are unchanged.")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

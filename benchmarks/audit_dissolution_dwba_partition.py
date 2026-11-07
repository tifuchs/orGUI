"""Report inherited native optical-profile sensitivity to interface support.

This diagnostic does not declare convergence or change any tolerance. Run from
the checkout with python -m benchmarks.audit_dissolution_dwba_partition.
The synthetic fixtures require no measured data.
"""

import copy
import json

import numpy as np

from orgui.datautils.xrayutils.CTRdistributions import (
    RecenteredSkellamProfile, SkellamProfile,
)
from orgui.datautils.xrayutils.CTRincoherent import FiniteFilmDissolution
from orgui.datautils.xrayutils.test.test_CTRdissolution import _crystal, _hkl


def _errors(actual, reference):
    return dict(max_relative=float(np.max(abs(actual/reference-1))),
                max_scale_normalized=float(np.max(abs(actual-reference))
                                           / np.max(abs(reference))))


def main():
    """Compare unwrapped native crystals and both dissolution field policies."""
    source = _crystal(width=40, interface=True)
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
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Synthetic complete-film dissolution; no measured data or fit is modified."""

import copy
import json

import numpy as np

from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRdistributions import RecenteredSkellamProfile
from orgui.datautils.xrayutils.CTRfilm import Film, SkellamSurface
from orgui.datautils.xrayutils.CTRincoherent import FiniteFilmDissolution
from orgui.datautils.xrayutils.CTRuc import UnitCell


def example(removable_substrate=False):
    """Build a synthetic stack with a fixed, explicitly finite inventory.

    All lattice constants are Angstrom and energy is eV. The optional carbon
    substrate Film is specified from the start; it is not an automatic
    repartition of an existing bulk structure.
    """
    bulk = UnitCell([3., 3., 6.], [90., 90., 90.], name="protected_bulk")
    bulk.addAtom("C", [0., 0., 0.], .1, .1, 1., layer=0)
    material = UnitCell([3., 3., 6.], [90., 90., 90.], name="material")
    material.addAtom("O", [0., 0., 0.], .1, .1, .7, layer=0)
    components, owners = [], []
    if removable_substrate:
        components.append(Film(copy.deepcopy(bulk), width_layers=3,
                               name="finite_substrate"))
        owners.append("finite_substrate")
    components.extend([Film(material, width_layers=5, name="film"),
                       SkellamSurface(copy.deepcopy(material), name="surface")])
    owners.append("film")
    crystal = SXRDCrystal(bulk, *components,
                          stacking=np.arange(1, len(components)+1))
    crystal.setEnergy(10000.)
    crystal.dwba.set_ctr_geometry(alpha_i=.03)
    model = FiniteFilmDissolution(
        crystal, films=owners, surface="surface",
        profile=RecenteredSkellamProfile(1., 0., 8., -5.),
        incoherent_fraction=.2)
    print(json.dumps(model.statistics(), indent=2))
    print("Kinematic intensity:", model.F2(1., 0., .7))
    print("DWBA reflectivity:", model.dwba.evaluate(1., 0., .7).reflectivity)


if __name__ == "__main__":
    example()
    example(removable_substrate=True)

"""Compare perpendicular and periodic surface cells at the same physical Q.

For the published TbScO3 Pnma CIF, run::

    python examples/CTR/surface_stacking_example.py path/to/TbScO3_published_Pnma.cif

The default (101) indices refer to the input CIF setting. Supply another plane
with ``--miller h k l``. No preferred termination or surface energy is inferred.
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from orgui.datautils.xrayutils import CTRsymmetry  # noqa: E402


def main():
    """Build both representations and verify their complex bulk amplitudes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif", type=Path)
    parser.add_argument("--miller", nargs=3, type=int, default=(1, 0, 1))
    args = parser.parse_args()

    # New default: a and b in-plane, c perpendicular; gamma may differ from 90.
    surfaces = CTRsymmetry.surface_terminations(args.cif, miller=args.miller)
    # Compatibility option: retain the original periodic lattice cell.
    lattices = CTRsymmetry.surface_terminations(
        args.cif, miller=args.miller, cell_mode="lattice",
    )
    parent_hkl = np.array([[0.31, 1.2], [0.41, -0.2], [0.37, 1.43]])
    for option, legacy in zip(surfaces, lattices, strict=True):
        cell = option.build_unitcell(
            f"surface_{option.id}", displacement="isotropic_equivalent",
        )
        old = legacy.build_unitcell(
            f"lattice_{legacy.id}", displacement="isotropic_equivalent",
        )
        print(option.id, [layer["composition"] for layer in option.layer_sequence])
        print("  lengths / Angstrom:", cell.a)
        print("  angles / degrees:", np.rad2deg(cell.alpha))
        print("  stacking_translation:", cell.stacking_translation)
        print("  physical repeat / Angstrom:", cell.repeat_cartesian)

        # Transform the SAME parent reciprocal points into each coordinate
        # basis. Reusing identical numeric h,k,l across modes changes physical Q.
        hkl_surface = option.surface_spec.transform.T @ parent_hkl
        hkl_lattice = legacy.surface_spec.transform.T @ parent_hkl
        cell.setEnergy(15000)  # eV for the existing CTR API
        old.setEnergy(15000)
        # Equal attenuation per repeat; avoid singular unattenuated Bragg points.
        np.testing.assert_allclose(
            cell.F_bulk(*hkl_surface, atten=0.1),
            old.F_bulk(*hkl_lattice, atten=0.1),
            rtol=1e-9, atol=1e-8,
        )
    print("Both representations agree at identical physical scattering points.")


if __name__ == "__main__":
    main()

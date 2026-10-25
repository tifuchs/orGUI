"""Generate rutile (110) CTR surface cells from PyXtal Wyckoff metadata.

This example shows the intended workflow for symmetry-aware CTR models:

1. start from a parent crystal structure,
2. enumerate distinct bulk terminations using a matrix or Miller indices,
3. build a candidate, or preserve an existing fixed-origin surface cell,
4. inspect couplings and optionally add a symmetry-preserving fit parameter.

The generated atoms are normal CTR atoms. Symmetry is metadata until a user
explicitly adds a Wyckoff parameter.
"""

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from pymatgen.core import Lattice, Structure  # noqa: E402

from orgui.datautils.xrayutils import CTRsymmetry  # noqa: E402


def rutile_parent_structure(metal, parent_a, oxygen_u, dw=None):
    """Create a conventional rutile parent structure.

    :param str metal:
        Metal element symbol, for example ``"Ru"`` or ``"Ti"``.
    :param tuple parent_a:
        Parent conventional rutile lattice constants in Angstrom.
    :param float oxygen_u:
        Oxygen Wyckoff ``4f`` internal coordinate.
    :param dw:
        Optional isotropic B parameter in Angstrom squared, stored on each
        source site. None leaves displacement data unspecified.
    :returns:
        Pymatgen structure accepted by PyXtal.
    :rtype:
        pymatgen.core.Structure
    """
    u = oxygen_u
    return Structure(
        Lattice.tetragonal(parent_a[0], parent_a[2]),
        [metal, metal, "O", "O", "O", "O"],
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [u, u, 0.0],
            [1.0 - u, 1.0 - u, 0.0],
            [0.5 - u, 0.5 + u, 0.5],
            [0.5 + u, 0.5 - u, 0.5],
        ],
        site_properties={} if dw is None else {"B_iso": [dw] * 6},
    )


def make_rutile_110_surface(name, metal, parent_a, oxygen_u, dw):
    """Create the original fixed-origin rutile (110) ``UnitCell``.

    Keep the established z=0 cut and two-layer grouping for existing fits.
    ``surface_terminations`` below demonstrates midpoint cuts separately.

    :param str name:
        Unit-cell name.
    :param str metal:
        Metal element symbol, for example ``"Ru"`` or ``"Ti"``.
    :param tuple parent_a:
        Parent conventional rutile lattice constants in Angstrom.
    :param float oxygen_u:
        Oxygen Wyckoff ``4f`` internal coordinate.
    :param float dw:
        In-plane and out-of-plane B parameter in Angstrom squared.
    :returns:
        Generated CTR unit cell with attached symmetry metadata.
    :rtype:
        orgui.datautils.xrayutils.CTRuc.UnitCell
    """
    parent = rutile_parent_structure(metal, parent_a, oxygen_u, dw=dw)
    surface_spec = replace(
        CTRsymmetry.rutile_110_surface_spec(parent_a), translation_range=None,
    )
    return CTRsymmetry.surface_unitcell_from_seed(
        parent, surface_spec, name, tol=1e-3, fallback_B=None,
    )


def print_wyckoff_summary(unitcell):
    """Print generated sites and affine couplings.

    :param orgui.datautils.xrayutils.CTRuc.UnitCell unitcell:
        Generated unit cell with attached symmetry metadata.
    """
    print(f"\n{unitcell.name}")
    print("surface cell a/Angstrom:", [round(value, 5) for value in unitcell.a])
    print("atoms:", len(unitcell.names))

    for site in unitcell.wyckoff_sites():
        print(
            f"site {site['site_id']} ({site['wyckoff_label']}): "
            f"atoms={site['atom_indices']} status={site['status']}"
        )

    print("oxygen u couplings:")
    for coupling in unitcell.wyckoff_couplings("O_4f"):
        expression = f"{coupling.constant:+.5f} {coupling.factor:+.1f}*u"
        print(
            f"  atom {coupling.atom_index:02d} "
            f"{unitcell.names[coupling.atom_index]}.{coupling.coordinate} = "
            f"{expression}"
        )


def main():
    """Run the rutile (110) Wyckoff metadata example."""
    parent = rutile_parent_structure(
        "Ru", (4.653255, 4.653255, 2.9692), 0.30569, dw=0.5033,
    )
    options = CTRsymmetry.surface_terminations(
        parent,
        # Right-handed candidate basis: A=[-1,1,0], B=[0,0,1], C=[1,1,0].
        transform=[[-1, 0, 1], [1, 0, 1], [0, 1, 0]],
        cell_mode="lattice",  # also the default for a full matrix
        layer_tolerance=0.25,  # Angstrom along the true surface normal
    )
    for option in options:
        print(option.id, "equivalent cuts:", option.equivalent_cuts)
        print("layers from surface inward:", [
            layer["composition"] for layer in option.layer_sequence
        ])
    # First candidate is an illustration, not a predicted stable termination.
    candidate = options[0].build_unitcell("RuO2_midpoint_cut")
    print("Candidate z range (crystal below interface at zero):",
          candidate.basis[:, 3].min(), candidate.basis[:, 3].max())

    # Miller indices default to c along the surface normal. The crystal repeat
    # can have a lateral component, stored separately from this coordinate c.
    automatic = CTRsymmetry.surface_terminations(
        parent, miller=(1, 1, 0), in_plane=(1, -1, 0),
    )
    print("Automatic surface basis (columns):", automatic[0].surface_spec.transform)
    surface = automatic[0].build_unitcell("RuO2_surface")
    print("Surface angles (degrees):", np.rad2deg(surface.alpha))
    print("Crystal repeat (surface fractional):", surface.stacking_translation)
    # Use cell_mode="lattice" with Miller indices to recover the oblique basis.

    # Preserve the existing fixed-origin examples and their fit conventions.
    ruo2 = make_rutile_110_surface(
        name="RuO2",
        metal="Ru",
        parent_a=(4.653255, 4.653255, 2.9692),
        oxygen_u=0.30569,
        dw=0.5033,
    )
    tio2 = make_rutile_110_surface(
        name="TiO2",
        metal="Ti",
        parent_a=(4.653255, 4.653255, 2.9692),
        oxygen_u=0.30458,
        dw=0.3719,
    )

    print_wyckoff_summary(ruo2)
    print_wyckoff_summary(tio2)

    ruo2.addWyckoffParameter(
        "O_4f",
        "u",
        absolute_limits=(0.28, 0.34),
    )
    print("\nAfter adding symmetry-preserving RuO2 oxygen u parameter:")
    print(ruo2.wyckoff_sites()[1])
    print("initial fit value is absolute u:", ruo2.getInitialParameters())

    ruo2.setFitParameters([0.31569])
    print("first few oxygen coordinates after delta_u = +0.01:")
    oxygen_indices = sorted(
        {coupling.atom_index for coupling in ruo2.wyckoff_couplings("O_4f")}
    )
    for index in oxygen_indices[:4]:
        xyz = [round(value, 5) for value in ruo2.basis[index, 1:4]]
        print(index, ruo2.names[index], xyz)

    tio2.addFitParameter(
        ([0], ["z"]),
        limits=(0.0, 1.0),
        name="manual_top_oxygen_z",
    )
    print("\nManual single-atom TiO2 fit intentionally lowers symmetry:")
    print(tio2.wyckoff_sites()[1])


if __name__ == "__main__":
    main()

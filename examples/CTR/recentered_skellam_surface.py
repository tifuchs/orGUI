"""Compare recentered height profiles and coherent CTRs at fixed mean coverage.

Run from a checkout with ``python examples/CTR/recentered_skellam_surface.py``.
Use ``--output path.png`` to save a figure without opening a plotting window.
Time is in seconds, rates in inverse seconds, and height changes in structural
layers. The illustrative cell has one equivalent carbon layer per repeat, so
coverage conservation is also material conservation.
"""

import argparse
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from orgui.datautils.xrayutils.CTRcalc import (  # noqa: E402
    Film, RecenteredSkellamProfile, SkellamSurface, SXRDCrystal, UnitCell,
)


def make_crystal(profile):
    """Build an 80-layer carbon Film with the supplied height profile."""
    bulk = UnitCell([3., 3., 3.], [90., 90., 90.], name="bulk")
    bulk.addAtom("C", [0., 0., 0.], .1, .1, 1., layer=0)
    film_cell = UnitCell.fromStr(bulk.toStr())
    film_cell.name = "film_cell"
    surface_cell = UnitCell.fromStr(bulk.toStr())
    surface_cell.name = "surface_cell"
    film = Film(film_cell, width_layers=80, name="film")
    surface = SkellamSurface(surface_cell, profile=profile, name="surface")
    # Weak bulk convergence damping avoids dominating the identical-material
    # Film/bulk cancellation away from Bragg poles in this roughness example.
    crystal = SXRDCrystal(
        bulk, film, surface, stacking=np.array([1, 2]), atten=1e-8)
    crystal.setEnergy(20000.)  # UnitCell energy API uses eV.
    return crystal


def main():
    """Plot probabilities, occupancies, and CTRs; report conservation errors."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Save the figure to this path")
    args = parser.parse_args()
    if args.output:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    profiles = {
        "Flat": RecenteredSkellamProfile(),
        "Symmetric": RecenteredSkellamProfile(1, .4, .4),
        "Asymmetric": RecenteredSkellamProfile(1, .1, .7),
        "Recentered Poisson": RecenteredSkellamProfile(1, 0, .8),
    }
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    layers = np.arange(-8, 9)
    ell = np.linspace(.15, 2.85, 651)
    for label, profile in profiles.items():
        low, high = profile.support()
        coverage = profile.correction(np.arange(low, high + 1)).sum()
        print(f"{label}: coverage residual={coverage:.3e} layers, "
              f"variance={profile.variance:.6g} layers^2")
        axes[0].plot(layers, profile.probability(layers), ".-", label=label)
        axes[1].step(layers, profile.occupancy(layers), where="mid", label=label)
        amplitude = make_crystal(profile).F(np.zeros_like(ell), np.zeros_like(ell), ell)
        axes[2].semilogy(ell, abs(amplitude)**2, label=label)
    axes[0].set(xlabel="Height change (layers)", ylabel="Terrace probability")
    axes[1].set(xlabel="Layer offset", ylabel="Material occupancy")
    axes[2].set(xlabel="L (r.l.u.)", ylabel="Coherent |F|² (electrons²)",
                title="(0, 0, L) rod")
    for axis in axes:
        axis.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.output, dpi=160)
        print(f"Saved {args.output.resolve()}")
    else:
        plt.show()


if __name__ == "__main__":
    main()

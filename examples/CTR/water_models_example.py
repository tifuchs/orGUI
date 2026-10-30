"""Create terrace-following step and layered water model files.

Run from the repository root with
``python -m examples.CTR.water_models_example``. Output files are written
beside this script. Lengths below are in Angstrom; WaterModel stores them
as fractions of its surface-normal repeat height, volume / uc_area.
"""

from pathlib import Path

import numpy as np

from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRfilm import PoissonSurface
from orgui.datautils.xrayutils.CTRuc import WaterModel


def create_example(kind, optical_reference=None, A0=0.8):
    """Return RuO2/TiO2 with the requested water profile at 20 keV.

    :param str kind: A type listed in ``WaterModel.supportedTypes``.
    :param str optical_reference: ``tanh``, ``erf``, or ``sharp``; optional.
    :param float A0: Adsorbed molecules per cell (step) or first-layer weight
        (layered). Set to zero for a step without adsorption.
    """
    folder = Path(__file__).resolve().parent
    crystal = SXRDCrystal.fromFile(folder / "RuO2_TiO2_Poisson_etching.xtal")
    surface = next(c for c in crystal.uc_surface_list if isinstance(c, PoissonSurface))
    surface.basis[1] = 0.15  # Smooth terraces, as in the DWBA notebook.
    cell = crystal.uc_bulk
    height = cell.volume / cell.uc_area
    water = WaterModel(
        cell.a.copy(), np.rad2deg(cell.alpha), watertype=kind,
        optical_reference=optical_reference, name=f"water_{kind}"
    )
    water.setWaterParameters(
        zpos=3.0 / height,
        layer_spacing=2.8 / height,
        sigma_0=1.0 / height,
        sigma_bar=1.4 / height,
        A0=A0,
    )
    water.pw = 0.0334  # Water molecules / Angstrom^3.
    surface.attach_water(water)
    crystal.setEnergy(20000.0)  # Core API takes eV; .xtal displays keV.
    return crystal


def main():
    """Write complete files and verify their serialized physical settings."""
    folder = Path(__file__).resolve().parent
    for suffix, kind, reference in (
        ("tanh", "step", "tanh"),
        ("erf", "step", "erf"),
        ("layered", "layered", "sharp"),
    ):
        crystal = create_example(kind, reference)
        path = folder / f"RuO2_TiO2_water_{suffix}.xtal"
        crystal.toFile(path)
        restored = SXRDCrystal.fromFile(path)
        surface = next(
            c for c in restored.uc_surface_list if isinstance(c, PoissonSurface)
        )
        assert surface.water.type == kind
        assert surface.water.pw == 0.0334
        original = next(
            c for c in crystal.uc_surface_list if isinstance(c, PoissonSurface)
        ).water
        assert surface.water._settings() == original._settings()
        # Text parameters are rounded to five decimal places.
        np.testing.assert_allclose(surface.water.basis, original.basis, atol=5.1e-6)
        print(path.name)
        print(f"AttachedWater {original.name}\n{original.toStr(False)}")


if __name__ == "__main__":
    main()

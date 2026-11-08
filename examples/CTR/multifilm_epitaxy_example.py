"""Construct and fit a validated two-film stack using synthetic layered cells.

Run from a checkout with orGUI installed, or with the checkout on PYTHONPATH::

    python examples/CTR/multifilm_epitaxy_example.py

The C/Si/Ge cells below illustrate the API; they are not physical crystal
structures. Distances are Angstrom and Film widths count structural layers.
"""

import copy

import numpy as np

from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal, StackValidationError
from orgui.datautils.xrayutils.CTRdistributions import SkellamProfile
from orgui.datautils.xrayutils.CTRfilm import EpitaxyInterface, Film
from orgui.datautils.xrayutils.CTRuc import UnitCell


def _cell(name, element, layer_ids=(1, 2)):
    cell = UnitCell([3., 3., 6.], [90.] * 3, name=name)
    for layer, z in zip(layer_ids, (0., .5)):
        cell.addAtom(element, [0., 0., z], .1, .1, 1., layer=layer)
        cell.layerpos[float(layer)] = z
    return cell


def _interface(top, bottom, name, layer_transition=None):
    interface = EpitaxyInterface(
        copy.deepcopy(top), copy.deepcopy(bottom), name=name,
        profile=SkellamProfile(.5, tail_probability=1e-10),
        layer_transition=layer_transition,
    )
    # Existing fitting conventions use basis_0 for fixed profile parameters.
    interface.basis_0[:] = interface.basis
    return interface


def main():
    """Evaluate a two-film model, couple disorder, and reject a thin Film."""
    substrate = _cell("substrate", "C", layer_ids=(42, 7))
    substrate_template = _cell("substrate_template", "C")
    material_A = _cell("A", "Si")
    material_B = _cell("B", "Ge")
    crystal = SXRDCrystal.from_stack(substrate, [
        # The same planes have different IDs in the substrate and template.
        # A lower termination at z=0 (42) continues at z=0.5 (2), and
        # termination at z=0.5 (7) continues at next-repeat z=0 (1).
        _interface(material_A, substrate_template, "substrate_A",
                   layer_transition={42: 2, 7: 1}),
        Film(copy.deepcopy(material_A), name="film_A", width_layers=40),
        _interface(material_B, material_A, "A_B"),
        Film(copy.deepcopy(material_B), name="film_B", width_layers=40),
    ])
    crystal.setEnergy(18000.)  # This existing API uses eV.
    for component in crystal.uc_surface_list:
        print(f"{component.name}: generated support "
              f"{component.pos_absolute:.3f} to "
              f"{component.height_absolute:.3f} Angstrom")

    crystal.addFitParameter({
        "substrate_A": {"atoms": [0], "par": ["oDW"],
                        "settings": {"unitcell": "top"}},
        "film_A": {"atoms": [0], "par": ["oDW"]},
        "A_B": {"atoms": [0], "par": ["oDW"],
                "settings": {"unitcell": "bottom"}},
    }, name="A_disorder", limits=(0., 2.))
    crystal.getStartParamAndLimits()
    crystal.setParameters([.2])
    ell = np.linspace(.13, 2.87, 30)
    amplitude = crystal.F(1., 0., ell)
    print("First three complex amplitudes:", amplitude[:3])

    # Native text retains strict validation. Fit parameter bindings are stored
    # separately through parametersToDict()/parametersFromDict() (e.g. HDF5).
    restored = SXRDCrystal.fromStr(crystal.toStr(showErrors=False))
    restored.parametersFromDict(crystal.parametersToDict())
    np.testing.assert_allclose(restored.F(1., 0., ell), amplitude)

    thin = copy.deepcopy(crystal)
    thin["film_A"].basis[0] = 12
    try:
        thin.F(1., 0., ell)
    except StackValidationError as exc:
        print("Expected thin-film rejection:", exc)
    else:
        raise AssertionError("Overlapping interfaces were accepted")


if __name__ == "__main__":
    main()

"""Optical references must conserve each internal plane's areal electron count."""

import copy

from ..CTRfilm import LegacySurfaceBoundary

import numpy as np
import pytest

from ..CTRcalc import Film, PoissonProfile, PoissonSurface, SXRDCrystal, UnitCell
from ..CTRcrystallinity import CrystallinityProfile
from ..CTRoptics import HC_KEV_ANGSTROM
from ..CTRutil import generate_surface_termination_cells


def _crystal(strain=1.0, repeats=1, roughness=0.0):
    cell = UnitCell([3.0, 3.0, 8.0], [90.0, 90.0, 90.0], name="four_planes")
    for index, element in enumerate(["La", "Ni", "La", "Ni"]):
        cell.addAtom(element, [0.0, 0.0, index / 4], 0.3, 0.3, 1.0, layer=index + 1)
        cell.layerpos[float(index + 1)] = index / 4
    bulk = copy.deepcopy(cell)
    cell.coherentDomainMatrix[0][2, 2] = strain
    film = Film(cell, width_layers=16, name="film")
    bank = generate_surface_termination_cells(
        cell.supercell((1, 1, repeats)), cell.layer_cycle
    )
    # This oracle audits legacy termination records and depth-dependent order.
    surface = PoissonSurface(
        bank, boundary=LegacySurfaceBoundary(),
        profile=PoissonProfile(roughness), name="surface")
    crystal = SXRDCrystal(bulk, film, surface, stacking=np.array([1, 2]))
    crystal.setEnergy(20000.0)
    return crystal, surface


def _assert_inventory(crystal):
    records, model, profile = crystal.dwba._atomic_model_snapshot()
    wavelength = HC_KEV_ANGSTROM / (crystal.uc_bulk._E * 1e-3)
    factor = 2.8179403262e-5 * wavelength**2 / (2 * np.pi)
    widths = np.diff(profile.boundaries)
    checked = 0
    for slot, (index, layer) in enumerate(model.layer_reference_layout):
        record = records[index]
        if record.role != "surface_termination":
            continue
        cell = record.cell
        mask = cell._crystallinity_labels() == layer
        forward = (
            cell.f[:, :5].sum(axis=1)
            + cell.f[:, 10]
            + cell.f[:, 11]
            + 1j * cell.f[:, 12]
        )
        atomic = np.sum(cell.basis[mask, 6] * forward[mask]) / cell.uc_area
        atomic *= sum(cell.coherentDomainOccupancy) * record.component_weight
        optical = model.layer_reference_shares[slot, ::-1]
        reference = np.dot(optical[1:-1, 0] + 1j * optical[1:-1, 1], widths) / factor
        np.testing.assert_allclose(
            reference,
            atomic,
            atol=2e-12,
            rtol=2e-12,
            err_msg=f"{record.record_name}, group {layer}",
        )
        checked += abs(atomic) > 1e-12
    assert checked >= 4


@pytest.mark.parametrize("strain", [1.0, 0.93, 1.08])
@pytest.mark.parametrize("repeats", [1, 2])
@pytest.mark.parametrize("roughness", [0.0, 1.3, -1.3])
def test_internal_plane_reference_inventory(strain, repeats, roughness):
    """Strained multi-plane slabs conserve both dispersive and absorptive inventory."""
    crystal, surface = _crystal(strain, repeats, roughness)
    _assert_inventory(crystal)
    surface.setCrystallinityProfile(CrystallinityProfile(PoissonProfile(2.0)))
    _assert_inventory(crystal)


@pytest.mark.parametrize("extension", ["xtal", "xpr"])
def test_optical_planes_survive_native_reload(tmp_path, extension):
    """Native persistence retains internal origins independently of the selector."""
    crystal, surface = _crystal(1.07, 2)
    # Atomic displacements do not redefine the structural optical-plane origins.
    for cell in surface.termination_cells.values():
        cell.basis[:, 3] += 0.037
        cell.basis_0[:] = cell.basis
    surface.setCrystallinityProfile(PoissonProfile(2.0))
    path = str(tmp_path / f"planes.{extension}")
    crystal.toFile(path)
    restored = SXRDCrystal.fromFile(path)
    restored.setEnergy(20000.0)
    _assert_inventory(restored)
    np.testing.assert_allclose(restored.optical_profile(), crystal.optical_profile())


def test_legacy_order_labels_recover_internal_planes():
    """Older terminations with surviving order labels recover their plane geometry."""
    crystal, _ = _crystal()
    text = "\n".join(
        line
        for line in crystal.toStr().splitlines()
        if not line.startswith("optical_plane:")
    )
    restored = SXRDCrystal.fromStr(text)
    restored.setEnergy(20000.0)
    _assert_inventory(restored)


def test_collapsed_supercell_keeps_each_physical_plane():
    """Repeating a collapsed termination repeats optical geometry and order IDs."""
    crystal, surface = _crystal()
    cell = surface.termination_cells[4.0]
    repeated = cell.supercell((1, 1, 2)).as_surface_termination(4.0)
    repeated.setEnergy(20000.0)
    original = cell.optical_profile()
    expected = np.concatenate((original, original + [cell.a[2], 0.0, 0.0]))
    np.testing.assert_allclose(repeated.optical_profile(), expected)

"""Read-only C1 replay on its original resolution grid and fixed rod scales.

Run with --root pointing to the LaNiO3 analysis directory. No fit is performed.
Outputs are written only to --output. The private dataset is not distributed.
"""

# Thread limits must precede NumPy/native-library imports.
# ruff: noqa: E402

import argparse
import copy
from contextlib import nullcontext
import hashlib
import json
import os
from pathlib import Path
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
from scipy.sparse import load_npz

from ._surface_boundary import _bounded_model
from orgui.datautils.xrayutils import CTRincoherent
from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRincoherent import (
    create_incoherent_model_from_config,
)
from orgui.datautils.xrayutils.CTRuc import UnitCell


def main():
    """Replay all four recorded domain rotations, without reoptimizing scales."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interface-fixed-ucs", type=int)
    parser.add_argument("--compare-complete-states", action="store_true",
                        help="Compare optimized kinematics to full native states")
    args = parser.parse_args()
    ctr = args.root / "CTR"
    xtal = ctr / "settings/1V2/0056_C1.xtal"
    metadata = xtal.with_suffix(".json")
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in (xtal, metadata)}
    assert hashes[xtal.name] == (
        "0e8c5f28b843d0a7877afa89e50322ed9ca40fe9ea63bd1cb3fc8956bc639b17")
    assert hashes[metadata.name] == (
        "63e5cb0951fa17401e6d7e183411183179ade2bdc963adc07e2cc9ddb8713adc")
    meta = json.loads(metadata.read_text())
    result = np.load(ctr / "data_out/1V2/0056/diagnostics/C1/result.npz")
    resolution = args.root / "FIT/FIT-046_1V2-bulk-dwf-reference"
    grid = np.load(resolution / "resolution_grid.npz")
    keep = np.any(grid["data"][:, :2] != 0, axis=1)
    data = grid["data"][keep]
    np.testing.assert_array_equal(data, result["data"])
    kernel = load_npz(resolution / "resolution_kernel.npz")[keep]
    used = np.asarray(kernel.getnnz(axis=0) > 0).ravel()
    kernel, hkl = kernel[:, used], grid["hkl"][used]
    crystal = SXRDCrystal.fromFile(str(xtal))
    crystal.setGlobalReferenceUnitCell(UnitCell(
        meta["reference_a_A"], meta["reference_angles_deg"], name="reference"))
    crystal.setEnergy(meta["energy_eV"])
    crystal.atten = meta["kinematic_attenuation"]
    crystal.domains = [[(np.eye(3), 1.)] for _ in crystal.uc_surface_list]
    config = meta["native_incoherent_config"]
    config["parameters"]["basis_0"] = np.asarray(config["parameters"]["basis_0"])
    old = create_incoherent_model_from_config(crystal, config)
    new = _bounded_model(
        crystal, films=["LaNiO3"], surface="LaNiO3_surface",
        incoherent_fraction=meta["physical"]["incoherent_fraction"])
    tight = _bounded_model(
        crystal, films=new.target_surface().boundary.films, surface=new.surface,
        profile=new.target_surface()._make_profile(),
        incoherent_fraction=new.incoherent_fraction, tail_probability=1e-12,
        moment_tolerance=1e-12, second_moment_tolerance=1e-12)
    models = [("legacy", old), ("dissolution", new), ("tight", tight)]
    if args.compare_complete_states:
        models.append(("complete", new))
    if args.interface_fixed_ucs is not None:
        expanded = copy.deepcopy(crystal)
        interface = expanded["STO_to_LaNiO3"]
        interface.fixed_ucs = args.interface_fixed_ucs
        interface.createInterfaceCells()
        models.append(("wider_interface", _bounded_model(
            expanded, films=new.target_surface().boundary.films, surface=new.surface,
        profile=new.target_surface()._make_profile(),
            incoherent_fraction=new.incoherent_fraction)))
    predictions, timings, raw_intensities = {}, {}, {}
    for name, model in models:
        start = time.perf_counter()
        intensities = []
        for number, rotation in enumerate(meta["domain_rotations"]):
            print(f"{name}: domain {number+1}/4", flush=True)
            model.coherent_model.domains = [
                [(np.asarray(rotation), 1.)]
                for _ in model.coherent_model.uc_surface_list]
            if name == "complete":
                from .benchmark_dissolution_optimized import _oracle
                context = _oracle("kinematic")
            else:
                context = nullcontext()
            with context:
                intensities.append(model.F2(*hkl.T))
        predictions[name] = np.sqrt(kernel @ np.mean(intensities, axis=0))
        timings[name] = time.perf_counter()-start
        if args.compare_complete_states and name in {"dissolution", "complete"}:
            raw_intensities[name] = np.asarray(intensities)
    _, groups = np.unique(data[:, :2], axis=0, return_inverse=True)
    scale = result["scales"][groups]
    residuals = {k: (scale*v-data[:, 3])/data[:, 4]
                 for k, v in predictions.items()}
    delta = residuals["dissolution"]-residuals["legacy"]
    numerical = residuals["tight"]-residuals["dissolution"]
    report = dict(
        input_hashes=hashes, import_path=CTRincoherent.__file__,
        observations=len(data), resolution_points=len(hkl), rotations=4,
        fixed_scales=result["scales"].tolist(), statistics=new.statistics(),
        evaluation_seconds=timings,
        legacy_saved_amplitude_max_abs=float(np.max(abs(
            predictions["legacy"]-result["amplitude"]))),
        chi2_N={k: float(v@v/len(v)) for k, v in residuals.items()},
        physical_change_max_sigma=float(np.max(abs(delta))),
        tighter_support_change_max_sigma=float(np.max(abs(numerical))))
    if args.compare_complete_states:
        for values in (raw_intensities, predictions):
            np.testing.assert_allclose(values["dissolution"], values["complete"],
                                       rtol=2e-12, atol=1e-9)
        reference = raw_intensities["complete"]
        nonzero = reference != 0
        report["optimized_node_intensity_max_relative_difference"] = float(
            np.max(abs(raw_intensities["dissolution"][nonzero]
                       / reference[nonzero]-1)))
        report["optimized_change_max_sigma"] = float(np.max(abs(
            residuals["dissolution"]-residuals["complete"])))
    if "wider_interface" in residuals:
        difference = residuals["wider_interface"]-residuals["dissolution"]
        report["interface_fixed_ucs"] = args.interface_fixed_ucs
        report["wider_interface_change_max_sigma"] = float(np.max(abs(difference)))
        report["wider_interface_statistics"] = models[-1][1].statistics()
    np.testing.assert_allclose(predictions["legacy"], result["amplitude"],
                               rtol=1e-9, atol=1e-9)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    assert report["tighter_support_change_max_sigma"] < 1e-3
    if "wider_interface" in residuals:
        assert report["wider_interface_change_max_sigma"] < 1e-3


if __name__ == "__main__":
    main()

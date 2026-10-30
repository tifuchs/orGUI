"""Benchmark water DWBA on the four-case notebook's 180-point CTR grid.

Run ``python -m examples.CTR.benchmark_water_DWBA --output tmp/water_before
--include-erf --profile`` from the repository root. Timings exclude model
construction and imports. JSON timings, complex-amplitude NPZ snapshots and
optional independent cold cProfile runs are written beneath the output folder.
"""

import argparse
import cProfile
import json
import os
from pathlib import Path
import platform
import pstats
import time

import numpy as np

from examples.CTR.water_models_example import create_example
from orgui.datautils.xrayutils.CTRfilm import PoissonSurface


def make_model(case, roughness):
    """Build a fresh 20 keV model; roughness is the Poisson fraction alpha."""
    family = "layered" if case == "layered" else "step"
    reference = (
        "sharp" if family == "layered" else ("erf" if case.endswith("erf") else "tanh")
    )
    a0 = 1.0 if "1ML" in case or family == "layered" else 0.0
    model = create_example(family, reference, A0=a0)
    surface = next(c for c in model.uc_surface_list if isinstance(c, PoissonSurface))
    surface.basis[1] = roughness
    if case == "dry":
        surface.attach_water(None)
    else:
        water = surface.water
        height = water.volume / water.uc_area
        water.setWaterParameters(
            3.0 / height, 2.8 / height, 0.60 / height, 0.40 / height, a0
        )
        water.adsorption_position = "onset"
    model.apply_stacking()
    return model


def benchmark(output, include_erf=False, profile=False, repeats=3):
    """Write timings in seconds and amplitude snapshots for all selected cases."""
    output.mkdir(parents=True, exist_ok=True)
    probe = make_model("dry", 0.0)
    height = probe.uc_bulk.volume / probe.uc_bulk.uc_area
    feature = height / 2.8
    ell = np.unique(
        np.r_[
            np.geomspace(0.006, 0.12, 18),
            np.linspace(0.13, 3.0, 121),
            np.linspace(feature - 0.18, feature + 0.18, 41),
        ]
    )
    cases = ["dry", "step_tanh", "step_1ML_tanh", "layered"]
    if include_erf:
        cases += ["step_erf", "step_1ML_erf"]
    rows, amplitudes = [], {"L": ell}
    for surface_name, alpha in (("rough", 0.15), ("flat", 0.0)):
        for case in cases:
            model = make_model(case, alpha)
            for h, k in ((0, 0), (0, 1), (1, 1)):
                model.dwba.set_ctr_geometry(equal_angles=True, rods=[(h, k)])
                hh, kk = np.full_like(ell, h), np.full_like(ell, k)
                start = time.perf_counter()
                prepared = model.dwba.prepare(hh, kk, ell)
                preparation = time.perf_counter() - start
                start = time.perf_counter()
                result = model.dwba.evaluate_prepared(prepared)
                evaluation = time.perf_counter() - start
                warm_prepared, warm_auto = [], []
                for _ in range(repeats):
                    start = time.perf_counter()
                    again = model.dwba.evaluate_prepared(prepared)
                    warm_prepared.append(time.perf_counter() - start)
                    np.testing.assert_array_equal(
                        again.total_amplitude, result.total_amplitude
                    )
                    start = time.perf_counter()
                    again = model.dwba.evaluate(hh, kk, ell)
                    warm_auto.append(time.perf_counter() - start)
                    np.testing.assert_array_equal(
                        again.total_amplitude, result.total_amplitude
                    )
                key = f"{surface_name}_{case}_{h}{k}L"
                amplitudes[key] = np.stack(
                    (result.total_amplitude, result.F_h, result.unperturbed_amplitude)
                )
                row = dict(
                    surface=surface_name,
                    model=case,
                    rod=f"{h}{k}L",
                    points=len(ell),
                    prepare_s=preparation,
                    evaluate_s=evaluation,
                    first_total_s=preparation + evaluation,
                    warm_prepared_s=float(np.median(warm_prepared)),
                    warm_auto_s=float(np.median(warm_auto)),
                    optical_interfaces=len(prepared.reference.n)-1,
                )
                if (h, k) == (0, 0):
                    from orgui.datautils.xrayutils.CTRoptics import solve_electric_field
                    angles = np.unique(prepared.alpha_i)
                    reference = prepared.reference
                    start = time.perf_counter()
                    if reference.continuous is None:
                        solve_electric_field(
                            reference.profile.values, reference.energy_eV,
                            angles, "s", boundaries=reference.profile.boundaries,
                        )
                    else:
                        reference.continuous.solve(angles, "s")
                    row["uncached_optical_solve_s"] = time.perf_counter()-start
                    if case != "dry" and not case.endswith("erf"):
                        from orgui.datautils.xrayutils._CTRwater_analytic import (
                            _water_contrast,
                        )
                        records = model.dwba._live_records()
                        index = next(i for i, r in enumerate(records)
                                     if r.role == "water")
                        kernel_times = []
                        for _ in range(7):
                            start = time.perf_counter()
                            _water_contrast(prepared, records[index], index)
                            kernel_times.append(time.perf_counter()-start)
                        row["uncached_water_kernel_s"] = float(np.median(kernel_times))
                    # Fit-like changes: distinguish frozen optics from an
                    # automatic preparation after changing the physical water.
                    old_atomic = model.uc_bulk.basis[0, 3]
                    try:
                        model.uc_bulk.basis[0, 3] += 0.01
                        start = time.perf_counter()
                        model.dwba.evaluate_prepared(prepared)
                        row["atomic_fit_s"] = time.perf_counter() - start
                    finally:
                        model.uc_bulk.basis[0, 3] = old_atomic
                    surface = next(
                        c
                        for c in model.uc_surface_list
                        if isinstance(c, PoissonSurface)
                    )
                    if surface.water is not None:
                        water = surface.water
                        old_width = water.basis[2]
                        try:
                            water.basis[2] *= 1.01
                            start = time.perf_counter()
                            model.dwba.evaluate_prepared(prepared)
                            row["water_fit_frozen_s"] = time.perf_counter() - start
                            # A second distinct update prevents the automatic
                            # path from hitting the just-computed frozen entry
                            # when the optical reference itself is unchanged.
                            water.basis[2] *= 1.01
                            start = time.perf_counter()
                            model.dwba.evaluate(hh, kk, ell)
                            row["water_fit_refreshed_s"] = time.perf_counter() - start
                        finally:
                            water.basis[2] = old_width
                rows.append(row)
                print(json.dumps(row), flush=True)
            if profile:
                fresh = make_model(case, alpha)
                fresh.dwba.set_ctr_geometry(equal_angles=True, rods=[(0, 0)])
                profiler = cProfile.Profile()
                profiler.runcall(
                    fresh.dwba.evaluate, np.zeros_like(ell), np.zeros_like(ell), ell
                )
                stem = output / f"{surface_name}_{case}"
                profiler.dump_stats(str(stem) + ".prof")
                with Path(str(stem) + ".txt").open("w", encoding="utf-8") as stream:
                    stats = pstats.Stats(profiler, stream=stream)
                    stats.sort_stats("cumulative").print_stats(35)
                    stats.sort_stats("tottime").print_stats(25)
                print(f"Profile saved: {stem.name}", flush=True)
            # Retain completed rows even if a later optional case fails.
            payload = dict(
                python=platform.python_version(),
                numpy=np.__version__,
                thread_environment={
                    name: os.environ.get(name) for name in (
                        "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"
                    )
                },
                repeats=repeats,
                rows=rows,
            )
            (output / "timings.json").write_text(
                json.dumps(payload, indent=2), encoding="utf-8"
            )
            np.savez_compressed(output / "amplitudes.npz", **amplitudes)


def main():
    """Parse benchmark options and run independent rough/flat water models."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--include-erf", action="store_true")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    benchmark(args.output, args.include_erf, args.profile, args.repeats)


if __name__ == "__main__":
    main()

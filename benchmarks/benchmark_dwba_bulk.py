"""Non-gating serial native unit-cell DWBA throughput benchmark.

The scattering geometry is generated with the persistent DWBA state, so the
benchmark measures the same call
sequence a grazing-incidence CTR scan uses: one preparation for the rod,
followed by repeated evaluations of the current atomic model. The timed
comparison reports the cost of the corrected optional empirical attenuation
against the default path, which omits its intra-cell weights entirely.
"""

import argparse
import copy
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRfilm import EpitaxyInterface, Film, SkellamProfile
from orgui.datautils.xrayutils.CTRuc import UnitCell

#: Rod and fixed incidence angle of the benchmarked grazing-incidence scan.
BENCHMARK_ROD = (1.0, 0.0)
BENCHMARK_ALPHA_I_DEG = 0.2
#: L range that stays inside the Ewald sphere for this cell and energy.
BENCHMARK_L_RANGE = (0.1, 3.0)


def _timed_evaluations(crystal, prepared, repetitions, attenuation):
    crystal.dwba.evaluate_prepared(
        prepared, bulk_attenuation=attenuation
    )
    start = time.perf_counter()
    for _ in range(repetitions):
        crystal.dwba.evaluate_prepared(
            prepared, bulk_attenuation=attenuation
        )
    return time.perf_counter() - start


def _benchmark(points, atoms, repetitions, trials, attenuation):
    rng = np.random.default_rng(20260831)
    cell = UnitCell([3.8, 3.8, 5.6], [90.0, 90.0, 90.0])
    for position in rng.random((atoms, 3)):
        cell.addAtom("Si", position, 0.15, 0.22, 1.0)
    cell.setEnergy(12000.0)
    h, k = BENCHMARK_ROD
    L = np.linspace(BENCHMARK_L_RANGE[0], BENCHMARK_L_RANGE[1], points)
    film = Film(copy.deepcopy(cell), name="film")
    film.basis[0] = 8.0
    interface = EpitaxyInterface(
        copy.deepcopy(cell),
        copy.deepcopy(cell),
        profile=SkellamProfile(0.75),
        name="interface",
    )
    models = (
        ("bulk", SXRDCrystal(cell)),
        ("film", SXRDCrystal(cell, film, stacking=np.array([1]))),
        (
            "interface",
            SXRDCrystal(cell, interface, stacking=np.array([1])),
        ),
    )
    print(
        f"({h:g} {k:g} L) rod, alpha_i={BENCHMARK_ALPHA_I_DEG:g} deg, "
        f"L in [{BENCHMARK_L_RANGE[0]:g}, {BENCHMARK_L_RANGE[1]:g}]"
    )
    for label, crystal in models:
        start = time.perf_counter()
        crystal.dwba.set_ctr_geometry(
            alpha_i=np.deg2rad(BENCHMARK_ALPHA_I_DEG)
        )
        prepared = crystal.dwba.prepare(h, k, L)
        preparation = time.perf_counter() - start
        record_count = len(crystal.dwba.evaluate_prepared(prepared).contributions)

        timings = {0.0: [], attenuation: []}
        for trial in range(trials):
            order = (0.0, attenuation) if trial % 2 == 0 else (attenuation, 0.0)
            for value in order:
                timings[value].append(
                    _timed_evaluations(
                        crystal, prepared, repetitions, value
                    )
                )
        zero_elapsed = float(np.median(timings[0.0]))
        attenuated_elapsed = float(np.median(timings[attenuation]))
        evaluations = points * repetitions
        print(f"[{label}] {record_count} records")
        print(
            f"  dwba.prepare: {preparation:.6f} s, "
            f"{points / preparation:,.0f} points/s"
        )
        print(
            "  dwba.evaluate_prepared, no empirical attenuation: "
            f"{zero_elapsed:.6f} s, {evaluations / zero_elapsed:,.0f} points/s"
        )
        print(
            f"  dwba.evaluate_prepared, attenuation={attenuation:g}: "
            f"{attenuated_elapsed:.6f} s, "
            f"{evaluations / attenuated_elapsed:,.0f} points/s"
        )
        print(
            "  active empirical-attenuation cost: "
            f"{attenuated_elapsed / zero_elapsed:.4f}x "
            f"({(attenuated_elapsed / zero_elapsed - 1.0) * 100.0:+.2f}%)"
        )


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=4096)
    parser.add_argument("--atoms", type=int, default=128)
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--trials", type=int, default=7)
    parser.add_argument("--attenuation", type=float, default=0.01)
    arguments = parser.parse_args()
    if arguments.trials < 1:
        parser.error("--trials must be at least 1")
    if arguments.attenuation <= 0.0:
        parser.error("--attenuation must be positive")
    _benchmark(
        arguments.points,
        arguments.atoms,
        arguments.repetitions,
        arguments.trials,
        arguments.attenuation,
    )


if __name__ == "__main__":
    _main()

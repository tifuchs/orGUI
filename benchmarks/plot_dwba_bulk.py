"""Non-gating serial native bulk-DWBA throughput benchmark."""

import argparse
import time

import numpy as np

from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRoptics import HC_KEV_ANGSTROM
from orgui.datautils.xrayutils.CTRuc import UnitCell


def _plot(points, atoms, repetitions):
    rng = np.random.default_rng(20260831)
    cell = UnitCell([3.8, 3.8, 5.6], [90.0, 90.0, 90.0])
    for position in rng.random((atoms, 3)):
        cell.addAtom("Si", position, 0.15, 0.22, 1.0)
    cell.setEnergy(12000.0)

    crystal = SXRDCrystal(cell)

    alpha_i = np.full(points, 0.18)
    alpha_f = np.linspace(0.08, 0.48, points)
    azimuth = np.linspace(0.02, 0.7, points)
    wavelength = HC_KEV_ANGSTROM / (cell._E * 1e-3)
    k0 = 2.0 * np.pi / wavelength
    incident_parallel = k0 * np.cos(alpha_i)
    exit_parallel = k0 * np.cos(alpha_f)
    Q_parallel = np.sqrt(
        incident_parallel**2
        + exit_parallel**2
        - 2.0 * incident_parallel * exit_parallel * np.cos(azimuth)
    )
    h = Q_parallel * cell.a[0] / (2.0 * np.pi)
    l_value = (
        k0
        * (np.sin(alpha_i) + np.sin(alpha_f))
        * cell.a[2]
        / (2.0 * np.pi)
    )
    prepared = crystal.dwba.prepare_from_glancing(
        h, 0.0, alpha_i, alpha_f
    )
    np.testing.assert_allclose(prepared.hkl[2], l_value)
    crystal.dwba.evaluate_prepared(prepared)

    start = time.perf_counter()
    for _ in range(repetitions):
        crystal.dwba.evaluate_prepared(prepared)
    elapsed = time.perf_counter() - start
    evaluations = points * repetitions
    print(
        f"{points} points x {atoms} atoms x {repetitions} repetitions: "
        f"{elapsed:.6f} s, {evaluations / elapsed:,.0f} points/s"
    )


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=4096)
    parser.add_argument("--atoms", type=int, default=128)
    arguments = parser.parse_args()
    _plot(arguments.points, arguments.atoms, arguments.repetitions)


if __name__ == "__main__":
    _main()

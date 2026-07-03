"""Compare CTR NumPy and Numba timings on the current host.

Run from the repository root:

    python benchmarks/verify_ctr_numba_host.py

For host-to-host comparisons, run the same command on each machine. To avoid
reusing stale compiled kernels across machines, set a host-local cache first:

    export NUMBA_CACHE_DIR=/tmp/numba-cache-$USER-$(hostname)
    rm -rf "$NUMBA_CACHE_DIR"
    python benchmarks/verify_ctr_numba_host.py

To check whether host-specific CPU code generation is responsible for timing
differences, also test with:

    NUMBA_CPU_NAME=generic python benchmarks/verify_ctr_numba_host.py
"""

from __future__ import annotations

import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

import numpy as np


def find_repo_root() -> Path:
    """Return the repository root containing ``orgui`` and ``pyproject.toml``."""
    path = Path(__file__).resolve()
    for candidate in (path.parent, *path.parents):
        if (candidate / "orgui").is_dir() and (
            candidate / "pyproject.toml"
        ).exists():
            return candidate
    raise RuntimeError("Could not find the orGUI repository root")


REPO_ROOT = find_repo_root()
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
os.environ.setdefault(
    "MPLCONFIGDIR",
    f"/tmp/matplotlib-{os.environ.get('USER', 'orgui')}",
)

from orgui.datautils import util  # noqa: E402
from orgui.datautils.xrayutils import CTRcalc, CTRuc  # noqa: E402


def print_command_output(command: list[str]) -> None:
    """Print command output if the command is available."""
    try:
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
        )
    except FileNotFoundError:
        print(f"{command[0]}: not available")
        return
    if completed.stdout:
        print(completed.stdout.strip())
    if completed.stderr:
        print(completed.stderr.strip())


def median_time(func, repeats: int) -> float:
    """Return median runtime in seconds after one warmup call."""
    func()
    timings = []
    for _ in range(repeats):
        start = time.perf_counter()
        func()
        timings.append(time.perf_counter() - start)
    return statistics.median(timings)


def load_crystal():
    """Load the bundled CTR regression crystal."""
    xpr_path = (
        REPO_ROOT
        / "orgui"
        / "datautils"
        / "xrayutils"
        / "test"
        / "testdata"
        / "0V12_calculated.xpr"
    )
    xtal = CTRcalc.SXRDCrystal.fromFile(xpr_path)
    pt100 = CTRcalc.UnitCell(
        [3.9242, 3.9242, 3.9242],
        [90.0, 90.0, 90.0],
    )
    xtal.setGlobalReferenceUnitCell(
        pt100,
        util.z_rotation(np.deg2rad(45.0)),
    )
    return xtal


def set_backend(backend: str) -> None:
    """Select the CTR acceleration backend."""
    CTRuc.set_accel_backend(backend)


def print_environment() -> None:
    """Print host, CPU, package, and relevant environment information."""
    print("=== Host ===")
    print("host", platform.node())
    print("platform", platform.platform())
    print("python", platform.python_version())
    print("repo", REPO_ROOT)

    print("\n=== CPU ===")
    print_command_output(
        [
            "bash",
            "-lc",
            "lscpu | egrep 'Model name|CPU\\(s\\)|Thread|Core|Socket|Flags|NUMA'",
        ]
    )

    print("\n=== Packages ===")
    try:
        import numba

        print("numba", numba.__version__)
        try:
            print("numba_num_threads", numba.get_num_threads())
        except Exception as exc:
            print("numba_num_threads", f"unavailable: {exc}")
    except Exception as exc:
        print("numba", f"unavailable: {exc}")
    print("numpy", np.__version__)

    print("\n=== Environment ===")
    for name in (
        "NUMBA_CPU_NAME",
        "NUMBA_CACHE_DIR",
        "NUMBA_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
    ):
        print(name, os.environ.get(name))


def run_benchmark(repeats: int = 5) -> None:
    """Run CTR backend timing comparisons."""
    xtal = load_crystal()
    print("\n=== CTR Benchmark ===")
    print("default_backend", CTRuc.CTR_ACCEL_BACKEND)
    print("cpp_accel_module", CTRuc.HAS_CPP_ACCEL)
    print("numba_accel_module", CTRuc.HAS_NUMBA_ACCEL)
    backends = ["numpy"]
    if CTRuc.HAS_CPP_ACCEL:
        backends.append("cpp")
    if CTRuc.HAS_NUMBA_ACCEL:
        backends.append("numba")

    for n_points in (2_000, 20_000, 200_000):
        h = np.zeros(n_points, dtype=np.float64)
        k = np.zeros(n_points, dtype=np.float64)
        l_values = np.ascontiguousarray(
            np.linspace(0.05, 7.0, n_points),
            dtype=np.float64,
        )

        print(f"\nN {n_points}")
        for name, func in (
            ("bulk.F_uc", lambda: xtal.uc_bulk.F_uc(h, k, l_values)),
            (
                "bulk.F_bulk",
                lambda: xtal.uc_bulk.F_bulk(h, k, l_values, xtal.atten),
            ),
            ("xtal.F", lambda: xtal.F(h, k, l_values)),
        ):
            timings = {}
            for backend in backends:
                set_backend(backend)
                timings[backend] = median_time(func, repeats)

            text = f"{name:12s}"
            for backend in backends:
                text += f" {backend}={timings[backend]:.6g}s"
            for backend in backends:
                if backend == "numpy":
                    continue
                text += (
                    f" {backend}_speedup="
                    f"{timings['numpy'] / timings[backend]:.2f}x"
                )
            print(text)


def main() -> None:
    """Print environment information and run the benchmark."""
    print_environment()
    run_benchmark()


if __name__ == "__main__":
    main()

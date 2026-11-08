"""Paired C1 timings against frozen complete-state native oracles.

Run as a module. Input hashes are checked by the supplied C1 loader. Timings
include batch setup and cleanup, exclude data loading, and never refit data.
"""

import argparse
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
from unittest.mock import patch

# This module sets numerical-library thread limits before importing NumPy.
from .benchmark_dissolution_performance import _model, _rays, _setup, _summary
from .audit_dissolution_dwba_reuse import _counts
import numpy as np

from orgui.datautils.xrayutils import CTRdwba
from orgui.datautils.xrayutils._CTRsurface_states import _SurfaceStatePlan
from orgui.datautils.xrayutils._CTRincoherent_dwba import SurfaceHeightDWBA
from orgui.datautils.xrayutils.test._dissolution_dwba_oracle import (
    _complete_dwba_batch, _complete_dwba_oracle,
)


def _kinematic_oracle(self, indices, h, k, ell):
    for index in indices:
        yield int(index), self.state(index).crystal.F(h, k, ell)


def _oracle(mode):
    stack = ExitStack()
    if mode == "kinematic":
        stack.enter_context(patch.object(_SurfaceStatePlan,
                                        "_iter_kinematic_amplitudes",
                                        _kinematic_oracle))
    else:
        stack.enter_context(patch.object(SurfaceHeightDWBA, "batch",
                                        _complete_dwba_batch))
        stack.enter_context(patch.object(SurfaceHeightDWBA, "_evaluate_in_batch",
                                        _complete_dwba_oracle))
    return stack


def _values(model, mode, hkl):
    if mode == "kinematic":
        return model.F2(*hkl)
    result = model.dwba.evaluate(*hkl)
    return np.array([result.reflectivity, result.coherent_intensity,
                     result.incoherent_intensity, abs(result.amplitude_prefactor)**2])


def main():
    """Compare identical dissolution physics before and after reuse."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=["kinematic", "shared", "state_specific"],
                        required=True)
    parser.add_argument("--points", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.output.exists() or args.repeats < 1 or args.points < 1:
        parser.error("Use a new output file, positive repeats and positive points")
    _setup(args.request_dir, None)
    root = Path(__file__).resolve().parents[1]
    files = [root / "orgui/datautils/xrayutils" / name for name in (
        "_CTRsurface_states.py", "CTRincoherent.py", "CTRdwba.py",
        "_CTRincoherent_dwba.py", "test/_dissolution_dwba_oracle.py")]
    files.append(Path(CTRdwba._CTRcalc_cpp.__file__))
    data = args.request_dir.resolve().parents[3] / "settings/1V2"
    files.extend(data / name for name in ("0056_C1.xtal", "0056_C1.json"))
    report = dict(mode=args.mode, points=args.points, threads=1,
                  platform=platform.platform(), python=platform.python_version(),
                  commit=subprocess.check_output(
                      ["git", "rev-parse", "HEAD"], text=True).strip(),
                  sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in files}, timings={}, max_relative_difference=0.)
    policy = "shared" if args.mode == "kinematic" else args.mode
    times = {mode: {"fresh": [], "repeat": []} for mode in ("optimized", "complete")}
    rays = _rays(args.points)
    for iteration in range(args.repeats):
        results = {}
        # Alternate order to reduce warm-host bias.
        order = ("optimized", "complete") if iteration % 2 == 0 else (
            "complete", "optimized")
        for mode in order:
            model = _model("dissolution", policy)
            with (_oracle(args.mode) if mode == "complete" else ExitStack()):
                start = time.perf_counter()
                results[mode] = _values(model, args.mode, rays)
                times[mode]["fresh"].append(time.perf_counter()-start)
                if args.mode != "kinematic":
                    with model.dwba.batch():
                        _values(model, args.mode, rays)
                        start = time.perf_counter()
                        repeated = _values(model, args.mode, rays)
                        times[mode]["repeat"].append(time.perf_counter()-start)
                    np.testing.assert_array_equal(repeated, results[mode])
        actual, expected = results["optimized"], results["complete"]
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-28)
        nonzero = expected != 0
        difference = float(np.max(abs(actual[nonzero]/expected[nonzero]-1)))
        report["max_relative_difference"] = max(
            difference, report["max_relative_difference"])
        print(f"Pair {iteration+1}/{args.repeats} passed", flush=True)
    report["timings"] = {mode: {key: _summary(values) for key, values in rows.items()
                               if values} for mode, rows in times.items()}
    if args.mode != "kinematic":
        report["work_counts"] = _counts(policy, "baseline", points=args.points)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2))
    print(json.dumps({key: report[key] for key in (
        "timings", "max_relative_difference")}, indent=2))


if __name__ == "__main__":
    main()

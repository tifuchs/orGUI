"""Alternate baseline/candidate fresh-process C1 timings on a changing host."""

import argparse
import json
from pathlib import Path
import statistics
import subprocess
import sys
import time

from benchmark_dwba_c1 import _load, _optical, _proposal, _rays, _sample, _setup


def main():
    """Measure three paired fresh, DWF-proposal and full-sample evaluations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--baseline-native", type=Path, required=True)
    parser.add_argument("--baseline-ensemble", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", choices=["fresh", "dwf", "sample", "shared"])
    parser.add_argument("--backend", choices=["before", "after"])
    parser.add_argument("--inner-repeats", type=int, default=1)
    parser.add_argument("--cases", nargs="+",
                        choices=["fresh", "dwf", "sample", "shared"],
                        default=["fresh", "dwf", "sample", "shared"])
    args = parser.parse_args()
    if args.inner_repeats < 1:
        parser.error("--inner-repeats must be positive")
    if args.case:
        import numpy as np
        import psutil
        before = args.backend == "before"
        _setup(args.request_dir, args.baseline_native if before else None,
               args.baseline_ensemble if before else None)
        samples, setups = [], []
        for _ in range(args.inner_repeats):
            if args.case == "sample":
                elapsed, values = _sample()
            else:
                start = time.perf_counter()
                model = _load("shared" if args.case == "shared" else "state_specific")
                if args.case == "dwf":
                    _proposal(model, "dwf", 0)
                setups.append(time.perf_counter() - start)
                start = time.perf_counter()
                values = _optical(model.dwba.evaluate(*_rays()))
                elapsed = time.perf_counter() - start
            samples.append(elapsed)
        np.save(args.output.with_suffix(".npy"), values)
        memory = psutil.Process().memory_info()
        args.output.write_text(json.dumps(dict(
            setup_seconds=statistics.median(setups) if setups else None,
            evaluation_seconds=statistics.median(samples), evaluation_samples=samples,
            peak_rss_bytes=getattr(memory, "peak_wset", memory.rss))))
        return
    if args.output.exists():
        parser.error("Choose a new output directory")
    args.output.mkdir(parents=True)
    import numpy as np
    report = {}
    for case in args.cases:
        report[case] = []
        for trial in range(3):
            pair = {}
            order = ("before", "after") if trial % 2 == 0 else ("after", "before")
            for backend in order:
                output = args.output / f"{case}_{trial}_{backend}.json"
                subprocess.run([
                    sys.executable, str(Path(__file__).resolve()),
                    "--request-dir", str(args.request_dir),
                    "--baseline-native", str(args.baseline_native.resolve()),
                    "--baseline-ensemble", str(args.baseline_ensemble.resolve()),
                    "--output", str(output), "--case", case, "--backend", backend,
                    "--inner-repeats", str(args.inner_repeats),
                ], check=True)
                pair[backend] = json.loads(output.read_text())
            before, after = [np.load(args.output / f"{case}_{trial}_{b}.npy")
                             for b in ("before", "after")]
            np.testing.assert_allclose(after, before, rtol=1e-9, atol=1e-28)
            pair["max_relative_error"] = float(np.max(
                abs(after - before) / np.maximum(abs(before), 1e-300)))
            report[case].append(pair)
            (args.output / "paired.json").write_text(json.dumps(report, indent=2))
            print(case, trial, pair, flush=True)


if __name__ == "__main__":
    main()

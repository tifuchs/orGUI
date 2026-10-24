"""Compare Python, serial C++, and threaded C++ overlap on identical batches.

Run ``python -m benchmarks.benchmark_sample_interception_native --frames 256``.
Lengths are metres and angles radians; no detector I/O or GUI work is timed.
"""

import argparse
from dataclasses import replace
import json
import os
from pathlib import Path
import platform
from time import perf_counter

import numpy as np

from .benchmark_sample_interception import _cases
from orgui.datautils.xrayutils.corrections import beamprofile as bp
from orgui.datautils.xrayutils.corrections import sample_interception as si


def _all_cases(frames, measured_points):
    """Add non-Gaussian profiles on both axes to the existing geometry cases."""
    circle = None
    for case in _cases(frames, measured_points):
        yield case
        if case.name == "circle":
            circle = case
    profiles = {
        "trapezoid": lambda width: bp.trapezoid_profile(width, width * .35),
        "triangle": bp.triangular_profile,
        "smoothed": lambda width: bp.smoothed_top_hat_profile(width, width * .08),
        "gennorm": lambda width: bp.generalized_normal_profile(width, 3.5),
        "skewnorm": lambda width: bp.skew_normal_profile(width, -2.5),
    }
    for name, make in profiles.items():
        yield replace(circle, name=f"circle_{name}", vertical=make(160e-6),
                      horizontal=make(.0028))


def main():
    """Validate values, time warmed batches, and count adaptive fallbacks."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=256)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 4, 0])
    parser.add_argument("--measured-points", type=int, nargs="+", default=[101, 501])
    parser.add_argument("--cases", nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.frames < 2 or args.repeats < 1 or min(args.workers) < 0
            or min(args.measured_points) < 3):
        parser.error("positive sizes, >= 3 knots and nonnegative workers required")
    native = si._sample_interception_cpp
    if native is None:
        parser.error("build the native interception extension first")
    report = dict(python=platform.python_version(), platform=platform.platform(),
                  numpy=np.__version__, logical_cpus=os.cpu_count(),
                  frames=args.frames, repeats=args.repeats, cases=[])
    scalar = si._frame_overlap
    fallback = []

    def counted(*values):
        """Count scalar fallbacks without changing their numerical behavior."""
        fallback.append(1)
        return scalar(*values)

    try:
        for case in _all_cases(args.frames, args.measured_points):
            if case.baseline or (args.cases and case.name not in args.cases):
                continue
            timings = {}
            reference = None
            for workers in [None, *args.workers]:
                si._sample_interception_cpp = None if workers is None else native
                si._frame_overlap = counted

                def run():
                    """Evaluate the same shape, profiles and acquisition frames."""
                    return si.overlap(case.shape, case.vertical, case.horizontal,
                                      case.alpha, case.azimuth, offset=case.offset,
                                      workers=0 if workers is None else workers)

                fallback.clear()
                actual = run()
                fallback_count = len(fallback)
                if reference is None:
                    reference = actual
                else:
                    np.testing.assert_allclose(actual.illumination,
                                               reference.illumination, rtol=1e-10)
                elapsed = []
                for _ in range(args.repeats):
                    start = perf_counter()
                    run()
                    elapsed.append(perf_counter() - start)
                name = "python" if workers is None else f"native_{workers}"
                timings[name] = dict(median_s=float(np.median(elapsed)),
                                     scalar_frames=fallback_count)
            report["cases"].append(dict(name=case.name, timings=timings))
            print(json.dumps(report["cases"][-1]), flush=True)
    finally:
        si._sample_interception_cpp, si._frame_overlap = native, scalar
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

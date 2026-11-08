"""Synthetic native surface-boundary timings against complete-state oracles.

Run as a module. No external data, fit or modified scientific tolerances are
required. Timings include fresh material snapshots and batch cleanup.
"""

import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np  # noqa: E402

from ._surface_boundary import _bounded_model  # noqa: E402
from orgui.datautils.xrayutils.CTRdistributions import RecenteredSkellamProfile  # noqa: E402
from orgui.datautils.xrayutils._CTRsurface_states import _SurfaceStatePlan  # noqa: E402
from orgui.datautils.xrayutils._CTRincoherent_dwba import SurfaceHeightDWBA  # noqa: E402
from orgui.datautils.xrayutils.test.test_CTRdissolution import _crystal  # noqa: E402
from orgui.datautils.xrayutils.test._dissolution_dwba_oracle import (  # noqa: E402
    _complete_dwba_batch,
    _complete_dwba_oracle,
)


def _complete_amplitudes(self, indices, h, k, ell):
    for index in indices:
        yield int(index), self.state(index).crystal.F(h, k, ell)


def main():
    """Report paired native boundary timings on a synthetic 25-plane Film."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument(
        "--mode",
        choices=["all", "kinematic", "shared", "state_specific"],
        default="all",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.points < 1 or args.repeats < 1:
        parser.error("points and repeats must be positive")
    if args.output is not None and args.output.exists():
        parser.error("Choose a new output file")
    hkl = (
        np.ones(args.points),
        np.zeros(args.points),
        np.linspace(0.4, 2.9, args.points),
    )
    modes = (
        ["kinematic", "shared", "state_specific"] if args.mode == "all" else [args.mode]
    )
    report = dict(points=args.points, repeats=args.repeats, threads=1, modes={})
    for mode in modes:
        model = _bounded_model(
            _crystal(width=25),
            films=("film",),
            surface="surface",
            profile=RecenteredSkellamProfile(1.0, 2.0, 5.0, -4.0),
            incoherent_fraction=0.2,
            optical_reference="shared" if mode == "kinematic" else mode,
        )
        native = model.coherent_model
        times = {key: [] for key in ("native_coherent", "domains", "complete_states")}
        error = 0.0
        for trial in range(args.repeats):
            values = {}
            order = list(times) if trial % 2 == 0 else list(reversed(times))
            for kind in order:
                with ExitStack() as stack:
                    if kind == "complete_states":
                        if mode == "kinematic":
                            stack.enter_context(
                                patch.object(
                                    _SurfaceStatePlan,
                                    "_iter_kinematic_amplitudes",
                                    _complete_amplitudes,
                                )
                            )
                        else:
                            stack.enter_context(
                                patch.object(
                                    SurfaceHeightDWBA, "batch", _complete_dwba_batch
                                )
                            )
                            stack.enter_context(
                                patch.object(
                                    SurfaceHeightDWBA,
                                    "_evaluate_in_batch",
                                    _complete_dwba_oracle,
                                )
                            )
                    target = native if kind == "native_coherent" else model
                    start = time.perf_counter()
                    value = (
                        target.F2(*hkl)
                        if mode == "kinematic"
                        else target.dwba.evaluate(*hkl).reflectivity
                    )
                    times[kind].append(time.perf_counter() - start)
                    values[kind] = value
            np.testing.assert_allclose(
                values["domains"], values["complete_states"], rtol=1e-10, atol=1e-28
            )
            error = max(
                error,
                float(
                    np.max(abs(values["domains"] - values["complete_states"]))
                    / np.max(abs(values["complete_states"]))
                ),
            )
        report["modes"][mode] = dict(
            states=len(model._snapshot()[1]),
            max_scale_normalized_difference=error,
            seconds=times,
            median_seconds={
                key: float(np.median(value)) for key, value in times.items()
            },
        )
        print(mode, report["modes"][mode]["median_seconds"], flush=True)
    if args.output is not None:
        args.output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    else:
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

"""Time footprint factors per frame, without detector I/O or ROI integration.

Run from the repository root with::

    python -m benchmarks.benchmark_sample_interception --output \
        benchmarks/baselines/sample_interception.json

Inputs use metres and radians at the scientific boundary. Each batch starts
with an empty overlap cache. Repeated geometry is cached only within a batch;
the scalar measurement starts a separate call for each frame. Batch statistics
are elapsed time divided by frame count, not a distribution of frame latencies.
"""

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import subprocess
from time import perf_counter
from types import SimpleNamespace

import numpy as np
import scipy

from orgui.app.integration_corrections import frame_correction_policy
from orgui.datautils.xrayutils.corrections import beamprofile as bp
from orgui.datautils.xrayutils.corrections.sample_interception import (
    SampleShape,
    overlap,
)


@dataclass
class _Case:
    name: str
    shape: SampleShape
    vertical: object
    horizontal: object
    horizontal_settings: dict
    alpha: np.ndarray
    azimuth: np.ndarray
    offset: tuple = (0.0, 0.0)
    baseline: bool = False


def _cases(frames, measured_points):
    vertical = bp.gaussian_profile(160e-6)
    wide = bp.top_hat_profile(0.02)
    narrow = bp.gaussian_profile(2.8e-3)
    wide_settings = {
        "analytical": True, "shape": "Top hat", "shape_values": [20000],
    }
    narrow_settings = {
        "analytical": True, "shape": "Gaussian", "shape_values": [2800],
    }
    square = SampleShape("rectangle", (0.01, 0.01))
    alpha = np.full(frames, np.deg2rad(0.36))
    azimuth = np.deg2rad(np.linspace(0, 90, frames))
    sweep_alpha = np.deg2rad(np.linspace(0.05, 2.0, frames))
    yield _Case(
        "legacy_1d_gaussian", square, vertical, wide, wide_settings,
        sweep_alpha, np.zeros(frames), baseline=True,
    )
    yield _Case(
        "aligned_square", square, vertical, wide, wide_settings,
        sweep_alpha, np.zeros(frames),
    )
    yield _Case(
        "rotating_square", square, vertical, wide, wide_settings, alpha, azimuth,
    )
    yield _Case(
        "rotating_offset_rectangle", SampleShape("rectangle", (0.01, 0.006)),
        vertical, narrow, narrow_settings, alpha, azimuth, offset=(0.002, 0.001),
    )
    yield _Case(
        "circle", SampleShape("circle", (0.01,)), vertical, narrow,
        narrow_settings, sweep_alpha, np.zeros(frames),
    )
    # Centred U shape: some surface sections have two disjoint intervals.
    vertices = (
        (-0.005, -0.005), (0.005, -0.005), (0.005, -0.002),
        (-0.002, -0.002), (-0.002, 0.002), (0.005, 0.002),
        (0.005, 0.005), (-0.005, 0.005),
    )
    yield _Case(
        "rotating_concave_polygon", SampleShape("polygon", vertices=vertices),
        vertical, narrow, narrow_settings, alpha, azimuth,
    )
    for points in measured_points:
        z = np.linspace(-480e-6, 480e-6, points)
        sigma = 160e-6 / np.sqrt(8 * np.log(2))
        measured = bp.MeasuredBeamProfile(z, np.exp(-0.5 * (z / sigma) ** 2))
        yield _Case(
            f"rotating_square_measured_{points}", square, measured, wide,
            wide_settings, alpha, azimuth,
        )
    yield _Case(
        "repeated_rotated_square", square, vertical, wide, wide_settings,
        alpha, np.full(frames, np.deg2rad(37.0)),
    )


def _settings(case, rtol):
    return {
        "version": 1,
        "enabled": not case.baseline,
        "shape": {
            "kind": case.shape.kind,
            "dimensions_m": list(case.shape.dimensions),
            "vertices_m": list(case.shape.vertices),
        },
        "offset_m": list(case.offset),
        "reference_incidence_deg": 0.36,
        "normal_rotation_confirmed": True,
        "azimuth_source": "phi",
        "azimuth_unit": "rad",
        "azimuth_sign": 1,
        "horizontal": case.horizontal_settings,
        "rtol": rtol,
    }


def _calls(case, rtol, scalar=False):
    middle = len(case.alpha) // 2
    alpha = case.alpha[middle:middle + 1] if scalar else case.alpha
    azimuth = case.azimuth[middle:middle + 1] if scalar else case.azimuth
    state = SimpleNamespace(sample_interception=_settings(case, rtol))
    scan = SimpleNamespace(phi=azimuth)

    def _core():
        if case.baseline:
            return case.vertical.corrections(alpha, 0.01)
        return overlap(
            case.shape, case.vertical, case.horizontal, alpha, azimuth,
            offset=case.offset,
            axis_vertical=-case.offset[0] * np.sin(np.deg2rad(0.36)),
            axis_horizontal=-case.offset[1], rtol=rtol,
        )

    def _policy():
        return frame_correction_policy(
            scan, state, alpha.size, use_normalization=False,
            use_illumination=True, alpha=alpha, beam_profile=case.vertical,
            sample_length=0.01,
        )

    return _core, _policy


def _time(call, repeats, frames):
    call()  # Untimed warm-up, including lazy profile setup.
    elapsed = []
    for _ in range(repeats):
        start = perf_counter()
        call()
        elapsed.append(perf_counter() - start)
    per_frame = np.asarray(elapsed) * 1000 / frames
    return {
        "median_ms_per_frame": float(np.median(per_frame)),
        "min_ms_per_frame": float(np.min(per_frame)),
        "max_ms_per_frame": float(np.max(per_frame)),
        "elapsed_seconds": elapsed,
    }


def main():
    """Benchmark warmed numerical factors and the application policy adapter."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frames", type=int, default=37)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--measured-points", type=int, nargs="+", default=[101, 501])
    parser.add_argument("--cases", nargs="+", help="Run only these case names")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.frames < 2 or args.repeats < 1 or min(args.measured_points) < 3:
        parser.error("frames >= 2, repeats >= 1 and measured points >= 3 required")
    # Retain warning checks in the adapter; suppress only console emission.
    import logging
    logging.getLogger("orgui.app.sample_interception_config").addHandler(
        logging.NullHandler()
    )
    logging.getLogger("orgui.app.sample_interception_config").propagate = False
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True,
    ).stdout.strip()
    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "revision": revision,
        "environment": {
            "platform": platform.platform(), "processor": platform.processor(),
            "logical_cpus": os.cpu_count(), "python": platform.python_version(),
            "numpy": np.__version__, "scipy": scipy.__version__,
        },
        "frames": args.frames, "repeats": args.repeats, "rtol": args.rtol,
        "selected_cases": args.cases,
        "scope": (
            "Single process; prebuilt vertical profile; core also prebuilds "
            "shape/horizontal profile. Policy includes horizontal construction, "
            "validation and provenance. No I/O, GUI, ROI or normalization. "
            "Scalar uses middle frame; batch cache resets on every call."
        ),
        "cases": [],
    }
    print("case | core batch ms/frame | policy batch ms/frame | scalar policy ms",
          flush=True)
    for case in _cases(args.frames, args.measured_points):
        if args.cases and case.name not in args.cases:
            continue
        core, policy = _calls(case, args.rtol)
        reference = core()
        applied = policy()
        if case.baseline:
            np.testing.assert_allclose(applied.illumination_divisor, reference[1])
        else:
            np.testing.assert_allclose(
                applied.intercepted_fraction, reference.fraction, rtol=1e-12,
            )
            np.testing.assert_allclose(
                applied.illumination_divisor,
                reference.effective_area / case.shape.area, rtol=1e-12,
            )
        _, scalar_policy = _calls(case, args.rtol, scalar=True)
        result = {
            "name": case.name, "settings": _settings(case, args.rtol),
            "vertical_profile_points": len(case.vertical.integration_points),
            "alpha_rad": case.alpha.tolist(), "azimuth_rad": case.azimuth.tolist(),
            "unique_geometries": len(set(zip(case.alpha, case.azimuth))),
            "core_batch": _time(core, args.repeats, args.frames),
            "policy_batch": _time(policy, args.repeats, args.frames),
            "policy_single_frame": _time(scalar_policy, args.repeats, 1),
        }
        report["cases"].append(result)
        medians = [result[key]["median_ms_per_frame"] for key in (
            "core_batch", "policy_batch", "policy_single_frame",
        )]
        print(f"{case.name} | " + " | ".join(f"{v:.4f}" for v in medians),
              flush=True)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()

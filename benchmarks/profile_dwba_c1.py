"""Profile the hash-checked C1 model with the current or a saved native kernel."""

import argparse
import cProfile
import json
from pathlib import Path
import pstats
import time

from benchmark_dwba_c1 import _load, _rays, _setup


def _kernel_parts(model, native):
    """Isolate work on the complete model's mean-state 512-point preparation.

    Empty-instance variants are profiling controls, not alternative predictions.
    Their difference from empty atom tables isolates factor setup/lookups; the
    remaining empty-table cost is native field/geometry setup and orchestration.
    """
    state = model.coherent_model.dwba
    result = {}
    with state.batch():
        prepared = state.prepare(*_rays())
        records, _, _ = state._atomic_model_snapshot()
        packed = state._pack_atomic_records(prepared, records)
        for mode in ("complete", "no_finite", "no_bulk", "no_instances", "no_atoms"):
            arrays = dict(packed)
            for name in arrays:
                if ((mode in ("no_finite", "no_instances", "no_atoms")
                     and name.startswith("finite_"))
                    or (mode in ("no_bulk", "no_instances", "no_atoms")
                        and name.startswith("bulk_") and name != "bulk_repeat")
                    or (mode == "no_atoms" and name in ("atom_basis", "form_factors"))):
                    arrays[name] = arrays[name][:0]
            state._evaluate_raw(prepared, records, arrays, "semi_infinite", 0.)
            native.reset_form_factor_cache_stats()
            start = time.perf_counter()
            for _ in range(5):
                state._evaluate_raw(prepared, records, arrays, "semi_infinite", 0.)
            result[mode] = dict(mean_seconds=(time.perf_counter() - start) / 5,
                                form_factor_cache=native.form_factor_cache_stats())
    return result


def main():
    """Save cold/warm Python profiles and native phase-work diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--ensemble", type=Path)
    parser.add_argument("--kernel-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    _setup(args.request_dir, args.native, args.ensemble)
    from orgui.datautils.xrayutils import CTRdwba
    native = CTRdwba._CTRcalc_cpp
    if args.kernel_only:
        report = _kernel_parts(_load(), native)
        (args.output / "kernel.json").write_text(json.dumps(report, indent=2))
        return
    original = CTRdwba.DWBAState._evaluate_raw
    work = []

    def _record(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        if hasattr(native, "_dwba_motif_work"):
            work.append(native._dwba_motif_work())
        return result

    CTRdwba.DWBAState._evaluate_raw = _record
    model = _load()
    report = {}
    try:
        with model.dwba.batch():
            for label in ("cold", "warm"):
                work.clear()
                profiler = cProfile.Profile()
                profiler.runcall(model.dwba.evaluate, *_rays())
                profiler.dump_stats(str(args.output / f"{label}.prof"))
                stats = pstats.Stats(profiler)
                with (args.output / f"{label}.txt").open("w") as stream:
                    stats.stream = stream
                    stats.strip_dirs().sort_stats("cumulative").print_stats(50)
                report[label] = dict(total_seconds=stats.total_tt,
                                     native_work=list(work))
    finally:
        CTRdwba.DWBAState._evaluate_raw = original
    (args.output / "profile.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

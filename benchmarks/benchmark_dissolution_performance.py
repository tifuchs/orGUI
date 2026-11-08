"""Read-only C1 dissolution and optical-position performance assessment.

Run as a module with --request-dir pointing to the supplied hash-checking C1
loader. Timings exclude loading, use one numerical-library thread, and never
fit or modify input data. Use --section to repeat an individual measurement.
"""

import argparse
import ast
from collections import Counter
from contextlib import ExitStack
import cProfile
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import pstats
import subprocess
import time
from unittest.mock import patch
import weakref

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np  # noqa: E402
import psutil  # noqa: E402

from .benchmark_dwba_c1 import _load, _setup  # noqa: E402
from orgui.datautils.xrayutils import CTRdwba, CTRoptics  # noqa: E402
from ._surface_boundary import _bounded_model  # noqa: E402
from orgui.datautils.xrayutils.CTRuc import UnitCell, _coherent_domain_arrays  # noqa: E402


def _rays(points):
    h = np.ones(points)
    return h, -h, np.linspace(2.4, 2.8, points)


def _model(kind, policy):
    legacy = _load(policy)
    if kind == "legacy":
        return legacy
    return _bounded_model(
        legacy.coherent_model, films=["LaNiO3"], surface="LaNiO3_surface",
        incoherent_fraction=legacy.incoherent_fraction, optical_reference=policy)


def _timed(function):
    start = time.perf_counter()
    result = function()
    return time.perf_counter() - start, result


def _summary(samples):
    return dict(seconds=samples, median_seconds=float(np.median(samples)),
                min_seconds=min(samples), max_seconds=max(samples))


def _old_function(filename, name):
    # Extract just the actual pre-fix function, not an approximation of its cost.
    source = subprocess.check_output(
        ["git", "show", "898bca1^:" + filename], text=True)
    node = next(n for n in ast.parse(source).body
                if isinstance(n, ast.FunctionDef) and n.name == name)
    scope = dict(np=np)
    exec(compile(ast.Module(body=[node], type_ignores=[]), filename, "exec"), scope)
    return scope[name]


def _fix_cost(repeats):
    captured = []
    original = CTRdwba._combine_tagged

    def capture(values, shares, *args, **kwargs):
        captured.append((values.copy(), shares.copy(), args, kwargs))
        return original(values, shares, *args, **kwargs)

    crystal = _load("shared").coherent_model
    with patch.object(CTRdwba, "_combine_tagged", capture):
        crystal.dwba._atomic_model_snapshot(raw_structure=True)
    old_tagged = _old_function(
        "orgui/datautils/xrayutils/CTRdwba.py", "_combine_tagged")
    old_public = _old_function(
        "orgui/datautils/xrayutils/CTRoptics.py", "combine_profiles")
    report = dict(captured_calls=len(captured),
                  rows=[len(v) for v, _, _, _ in captured],
                  share_shapes=[s.shape for _, s, _, _ in captured])
    for label, old, new in (("tagged", old_tagged, original),
                            ("public", old_public, CTRoptics.combine_profiles)):
        times = dict(before=[], after=[])
        for trial in range(repeats):
            for key in (("before", "after") if trial % 2 == 0
                        else ("after", "before")):
                function = old if key == "before" else new

                def run():
                    for _ in range(10):
                        for values, shares, args, kwargs in captured:
                            if label == "tagged":
                                function(values, shares, *args, **kwargs)
                            else:
                                function(values)

                elapsed, _ = _timed(run)
                times[key].append(elapsed / 10)
        report[label] = {key: _summary(value) for key, value in times.items()}
    # Paired end-to-end native evaluation; old and new intentionally differ in
    # optics, so record that difference rather than assert identical physics.
    times, results = dict(before=[], after=[]), {}
    for trial in range(repeats):
        for key in (("before", "after") if trial % 2 == 0 else ("after", "before")):
            crystal = _load("shared").coherent_model
            with patch.object(CTRdwba, "_combine_tagged",
                              old_tagged if key == "before" else original):
                elapsed, value = _timed(lambda: crystal.dwba.evaluate(*_rays(512)))
            times[key].append(elapsed)
            results[key] = value.reflectivity
    report["native_512"] = {key: _summary(value) for key, value in times.items()}
    report["native_512"]["relative_observable_change"] = float(np.max(
        abs(results["after"] / results["before"] - 1)))
    return report


def _work_counts(model):
    counts = Counter()
    with ExitStack() as stack:
        for name, key in (("_atomic_model_snapshot", "snapshot_builds"),
                          ("_pack_atomic_records", "packing_builds")):
            original = getattr(CTRdwba.DWBAState, name)

            def count(self, *args, _original=original, _key=key, **kwargs):
                before = self._cache_counters[_key]
                value = _original(self, *args, **kwargs)
                counts[_key] += self._cache_counters[_key] - before
                return value

            stack.enter_context(patch.object(CTRdwba.DWBAState, name, count))
        with model.dwba.batch():
            first = model.dwba.evaluate(*_rays(32))
            cold = dict(counts)
            counts.clear()
            second = model.dwba.evaluate(*_rays(32))
            warm = dict(counts)
        np.testing.assert_array_equal(first.reflectivity, second.reflectivity)
    return dict(batch_first=cold, batch_repeat=warm)


def _evaluate(kind, policy, repeats):
    samples, checks = [], []
    for _ in range(repeats):
        model = _model(kind, policy)
        elapsed, value = _timed(lambda: model.dwba.evaluate(*_rays(512)))
        samples.append(elapsed)
        checks.append(value.reflectivity)
    for value in checks[1:]:
        np.testing.assert_array_equal(value, checks[0])
    model = _model(kind, policy)
    start = time.perf_counter()
    with model.dwba.batch():
        entry = time.perf_counter() - start
        stages = {}
        for label, rays, options in (
                ("first", _rays(512), {}), ("repeat", _rays(512), {}),
                ("exit_p", _rays(512), dict(polarization_f="p")),
                ("new_rod", (np.ones(512), np.zeros(512),
                             np.linspace(2.4, 2.8, 512)), {})):
            elapsed, value = _timed(lambda: model.dwba.evaluate(*rays, **options))
            stages[label] = elapsed
            if label in {"first", "repeat"}:
                np.testing.assert_array_equal(value.reflectivity, checks[0])
    result = dict(fresh=_summary(samples), batch_entry_seconds=entry,
                  batch_calls_seconds=stages,
                  work_counts=_work_counts(_model(kind, policy)))
    if kind == "dissolution":
        _, indices, probabilities, _ = model._snapshot()
        result["states"] = len(indices)
        result["nonzero_states"] = int(np.count_nonzero(probabilities))
        result["index_range"] = indices[[0, -1]].tolist()
    else:
        candidates = model.target_surface()._profile_candidates()
        low, high, masses, *_ = model.target_surface()._height_state_selection(
            candidates, model.exact_layer_count)
        result["states"] = high - low + 1
        result["nonzero_states"] = int(np.count_nonzero(masses))
    return result


def _profile(policy):
    model = _model("dissolution", policy)
    profiler = cProfile.Profile()
    profiler.runcall(model.dwba.evaluate, *_rays(512))
    stats = pstats.Stats(profiler)
    rows = [dict(file=Path(file).name, line=line, function=name, calls=value[1],
                 self_seconds=value[2], cumulative_seconds=value[3])
            for (file, line, name), value in stats.stats.items()]
    return dict(instrumented_total_seconds=stats.total_tt,
                cumulative=sorted(rows, key=lambda r: -r["cumulative_seconds"])[:25],
                self_time=sorted(rows, key=lambda r: -r["self_seconds"])[:20])


def _memory(kind, policy):
    # Invoke each memory case in a fresh process; peak RSS is process-wide.
    process = psutil.Process()
    gc.collect()
    baseline = process.memory_info().rss
    model = _model(kind, policy)
    with model.dwba.batch():
        result = model.dwba.evaluate(*_rays(512))
        states = model.dwba._scope.get("states")
        references = ([weakref.ref(s) for s in states.values()]
                      if isinstance(states, dict) else
                      [weakref.ref(s[0]) for s in states] if states else [])
        del states
        inside = process.memory_info().rss
    after = process.memory_info().rss
    alive = sum(reference() is not None for reference in references)
    del result, model
    gc.collect()
    return dict(baseline_rss_bytes=baseline, in_batch_rss_bytes=inside,
                after_batch_rss_bytes=after,
                after_gc_rss_bytes=process.memory_info().rss,
                private_states_after_batch=alive,
                private_states_after_gc=sum(r() is not None for r in references))


def _kinematic_analysis(repeats):
    """Isolate repeated native evaluations with a temporary exact C1 cache.

    This diagnostic is intentionally not a production/general UnitCell cache.
    It assumes the validated C1 native cells without custom form-factor hooks.
    All patches and stored results are confined to one F2 call.
    """
    originals = dict(unit=UnitCell.F_uc, bulk=UnitCell.F_bulk)

    def evaluate(model, memoize):
        counts, caches = Counter(), dict(unit={}, bulk={})

        def wrap(kind):
            def call(cell, h, k, ell, *args):
                counts[kind + "_calls"] += 1
                if not memoize:
                    return originals[kind](cell, h, k, ell, *args)
                assert type(cell) is UnitCell
                assert not cell._special_formfactors_present
                basis, factors, _, order = cell._coherent_basis()
                domains, weights = _coherent_domain_arrays(
                    cell.coherentDomainMatrix, cell.coherentDomainOccupancy)
                arrays = (basis, factors, cell.refHKLTransform, cell.B_mat,
                          cell.R_mat, cell.R_mat_inv, domains, weights,
                          cell.stacking_translation, h, k, ell)
                key = (order, cell.uc_area, args, *(
                    (np.asarray(a).shape, np.asarray(a).dtype.str,
                     np.asarray(a).tobytes()) for a in arrays))
                if key not in caches[kind]:
                    counts[kind + "_misses"] += 1
                    caches[kind][key] = originals[kind](cell, h, k, ell, *args)
                return caches[kind][key]
            return call

        with patch.object(UnitCell, "F_uc", wrap("unit")), \
                patch.object(UnitCell, "F_bulk", wrap("bulk")), \
                patch.object(CTRdwba, "solve_electric_field",
                             side_effect=AssertionError("Unexpected DWBA work")), \
                patch.object(UnitCell, "optical_profile",
                             side_effect=AssertionError("Unexpected optical work")):
            elapsed, result = _timed(lambda: model.F2(*_rays(512)))
        return elapsed, result, dict(counts)

    times = dict(baseline=[], memoized=[])
    counts = {}
    reference = _model("dissolution", "shared").F2(*_rays(512))
    for trial in range(repeats):
        order = ("baseline", "memoized") if trial % 2 == 0 else ("memoized", "baseline")
        for kind in order:
            model = _model("dissolution", "shared")
            elapsed, result, counts[kind] = evaluate(model, kind == "memoized")
            np.testing.assert_array_equal(result, reference)
            times[kind].append(elapsed)
    model = _model("dissolution", "shared")
    profiler = cProfile.Profile()
    profiler.runcall(model.F2, *_rays(512))
    stats = pstats.Stats(profiler)
    rows = [dict(file=Path(file).name, line=line, function=name, calls=value[1],
                 self_seconds=value[2], cumulative_seconds=value[3])
            for (file, line, name), value in stats.stats.items()]
    return dict(timings={key: _summary(values) for key, values in times.items()},
                counts=counts, bitwise_equal=True, optical_work_forbidden=True,
                instrumented_total_seconds=stats.total_tt,
                profile=sorted(rows, key=lambda r: -r["cumulative_seconds"])[:35])


def main():
    """Measure separately the numerical fix, complete models and hot paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--section", choices=["all", "fix", "models", "profile",
                                             "scaling", "memory", "kinematic"],
                        default="all")
    parser.add_argument("--kind", choices=["legacy", "dissolution"],
                        default="dissolution")
    parser.add_argument("--policy", choices=["shared", "state_specific"],
                        default="shared")
    args = parser.parse_args()
    if args.output.exists() or args.repeats < 1:
        parser.error("Use a new output file and positive repeat count")
    _setup(args.request_dir, None)
    report = dict(platform=platform.platform(), python=platform.python_version(),
                  cpu=platform.processor(), threads=1, repeats=args.repeats,
                  commit=subprocess.check_output(
                      ["git", "rev-parse", "HEAD"], text=True).strip(),
                  native_sha256=hashlib.sha256(
                      Path(CTRdwba._CTRcalc_cpp.__file__).read_bytes()).hexdigest(),
                  request_dir=str(args.request_dir), measurements={})

    def save(key, value):
        report["measurements"][key] = value
        report["process_peak_rss_bytes"] = psutil.Process().memory_info().peak_wset
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2))
        print(key, "completed", flush=True)

    if args.section in {"all", "fix"}:
        save("optical_fix", _fix_cost(max(5, args.repeats)))
    if args.section in {"all", "models"}:
        for policy in ("shared", "state_specific"):
            for kind in ("legacy", "dissolution"):
                save(f"{kind}_{policy}", _evaluate(kind, policy, args.repeats))
        for kind in ("legacy", "dissolution"):
            times = []
            for _ in range(args.repeats):
                model = _model(kind, "shared")
                elapsed, values = _timed(lambda: model.F2(*_rays(512)))
                assert np.all(np.isfinite(values))
                times.append(elapsed)
            save(f"kinematic_{kind}", _summary(times))
    if args.section in {"all", "profile"}:
        for policy in ("shared", "state_specific"):
            save(f"profile_{policy}", _profile(policy))
    if args.section in {"all", "scaling"}:
        for policy in ("shared", "state_specific"):
            for points in (32, 2048):
                model = _model("dissolution", policy)
                elapsed, value = _timed(lambda: model.dwba.evaluate(*_rays(points)))
                assert np.all(np.isfinite(value.reflectivity))
                save(f"scaling_{policy}_{points}", dict(seconds=elapsed))
    if args.section == "memory":
        save(f"memory_{args.kind}_{args.policy}", _memory(args.kind, args.policy))
    if args.section == "kinematic":
        save("kinematic_redundancy", _kinematic_analysis(max(5, args.repeats)))


if __name__ == "__main__":
    main()

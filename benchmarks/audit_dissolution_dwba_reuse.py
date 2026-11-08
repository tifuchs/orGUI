"""Diagnostic reuse controls for dissolution; runtime implementation unchanged.

Run as a module with the hash-checking C1 --request-dir and a new --output.
Native batches and shared prepared fields are temporary controls, not a new
solver. Timings and instrumented work counts are collected separately.
"""

import argparse
from collections import Counter
from contextlib import contextmanager, ExitStack
from dataclasses import replace
import json
from pathlib import Path
import time
from unittest.mock import patch

from .benchmark_dissolution_performance import (
    _model, _rays, _setup, _summary, _timed,
)
import numpy as np

from orgui.datautils.xrayutils import CTRdwba, _CTRincoherent_dwba as dissolution
from orgui.datautils.xrayutils._CTRincoherent_dwba import _regrid_model, _shared_grid


@contextmanager
def _control(mode):
    if mode == "baseline":
        yield
        return
    original_batch = dissolution.SurfaceHeightDWBA.batch

    @contextmanager
    def batch(self):
        if self._scope is not None:
            yield self
            return
        with original_batch(self):
            with ExitStack() as stack:
                average = self._scope["average"].dwba
                for index, probability in zip(self._scope["indices"],
                                              self._scope["probabilities"]):
                    if not probability:
                        continue
                    state = self._scope["path"].state(int(index)).crystal.dwba
                    self._scope["states"][int(index)] = state
                    stack.enter_context(state.batch())
                    state._batch_profiles = average._batch_profiles
                    stack.callback(state.clear_cache)
                if mode == "state_batches_no_average":
                    assert self._model.optical_reference == "state_specific"

                    class PrefactorOnly:
                        # Fail if any other result is consumed. Reuse the native
                        # geometry-only prefactor property without its formula.
                        amplitude_prefactor = CTRdwba.DWBAResult.amplitude_prefactor

                        def __init__(self, prepared):
                            self.prepared = prepared

                    def prefactor_only(prepared, **kwargs):
                        return PrefactorOnly(prepared)

                    stack.enter_context(patch.object(
                        average, "evaluate_prepared", prefactor_only))
                yield self

    rebinned = {}

    def common_preparation(state, common, atomic_model, profile):
        # The production caller has already formed the all-state union grid.
        # Assert that reusing it introduces no optical-interface change.
        grid = _shared_grid(common.reference.profile, profile.boundaries)
        np.testing.assert_array_equal(grid.boundaries,
                                      common.reference.profile.boundaries)
        np.testing.assert_array_equal(grid.values, common.reference.profile.values)
        key = (state, common.reference.fingerprint)
        if key not in rebinned:
            rebinned[key] = _regrid_model(atomic_model, profile, grid.values[:, 0])
        return replace(common, crystal_identity=id(state._crystal),
                       bulk_identity=id(state._crystal.uc_bulk),
                       _atomic_model=rebinned[key])

    with ExitStack() as patches:
        patches.enter_context(patch.object(
            dissolution.SurfaceHeightDWBA, "batch", batch))
        if mode == "shared_prepared":
            patches.enter_context(patch.object(
                dissolution, "_shared_preparation", common_preparation))
        yield


def _values(result):
    return np.array([result.reflectivity, result.coherent_intensity,
                     result.incoherent_intensity])


def _check(actual, expected):
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-28)
    nonzero = expected != 0
    return float(np.max(abs(actual[nonzero] / expected[nonzero] - 1)))


def _counts(policy, mode, points=32):
    counts, unique_fields = Counter(), set()
    with ExitStack() as stack:
        for name, key in (("_atomic_model_snapshot", "snapshot_builds"),
                          ("_pack_atomic_records", "packing_builds")):
            original = getattr(CTRdwba.DWBAState, name)

            def count(self, *args, _original=original, _key=key, **kwargs):
                before = self._cache_counters[_key]
                result = _original(self, *args, **kwargs)
                counts[_key] += self._cache_counters[_key] - before
                return result

            stack.enter_context(patch.object(CTRdwba.DWBAState, name, count))
        for name in ("_prepare_explicit", "_evaluate_raw"):
            original = getattr(CTRdwba.DWBAState, name)

            def count(*args, _original=original, _key=name, **kwargs):
                counts[_key] += 1
                result = _original(*args, **kwargs)
                if _key == "_evaluate_raw":
                    work = CTRdwba._CTRcalc_cpp._dwba_motif_work()
                    for item in ("finite_cache_hits", "lateral_exponentials",
                                 "normal_exponentials"):
                        counts[item] += work[item]
                return result

            stack.enter_context(patch.object(CTRdwba.DWBAState, name, count))
        original_solve = CTRdwba.solve_electric_field

        def solve(profile, energy, angles, polarization, **kwargs):
            counts["field_solver_calls"] += 1
            unique_fields.add((float(energy), polarization,
                               np.asarray(profile).tobytes(),
                               np.asarray(angles).tobytes(),
                               np.asarray(kwargs["boundaries"]).tobytes()))
            return original_solve(profile, energy, angles, polarization, **kwargs)

        stack.enter_context(patch.object(CTRdwba, "solve_electric_field", solve))
        model = _model("dissolution", policy)
        with _control(mode), model.dwba.batch():
            first = _values(model.dwba.evaluate(*_rays(points)))
            cold = dict(counts, unique_field_solves=len(unique_fields))
            counts.clear()
            unique_fields.clear()
            repeated = _values(model.dwba.evaluate(*_rays(points)))
            warm = dict(counts, unique_field_solves=len(unique_fields))
        _check(repeated, first)
    return dict(first=cold, repeat=warm)


def _edge_checks(policy, modes):
    from ._surface_boundary import _bounded_model
    from orgui.datautils.xrayutils.CTRdistributions import RecenteredSkellamProfile
    from orgui.datautils.xrayutils.test.test_CTRdissolution import _crystal, _hkl
    checked, maximum = 0, 0.
    for near_boundary in (False, True):
        model = _bounded_model(
            _crystal(interface=True), films=["film"], surface="surface",
            profile=RecenteredSkellamProfile(offset=-5.5), optical_reference=policy,
            incoherent_fraction=.3)
        if near_boundary:
            model.target_surface().basis[3] = (
                model.coherent_model._surface_state_plan().minimum + .5)
        for incoming, outgoing in (("s", "s"), ("s", "p"),
                                   ("p", "s"), ("p", "p")):
            options = dict(polarization_i=incoming, polarization_f=outgoing,
                           bulk_attenuation=.025)
            expected = _values(model.dwba.evaluate(*_hkl(), **options))
            for mode in modes[1:]:
                with _control(mode):
                    actual = _values(model.dwba.evaluate(*_hkl(), **options))
                maximum = max(maximum, _check(actual, expected))
                checked += 1
    return dict(comparisons=checked, max_relative_difference=maximum)


def main():
    """Compare native reuse controls while retaining every physical state."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--policy", choices=["shared", "state_specific"], required=True)
    parser.add_argument("--points", type=int, default=512)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.output.exists() or args.repeats < 1 or args.points < 2:
        parser.error("Use a new output, positive repeats and at least two points")
    _setup(args.request_dir, None)
    modes = ["baseline", "state_batches"]
    if args.policy == "shared":
        modes.append("shared_prepared")
    report = dict(policy=args.policy, points=args.points, repeats=args.repeats,
                  note="Temporary controls only; runtime code unchanged.", modes={})

    def save():
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2))

    samples = {mode: [] for mode in modes}
    expected, maximum = None, {mode: 0. for mode in modes}
    for trial in range(args.repeats):
        for mode in (modes if trial % 2 == 0 else modes[::-1]):
            model = _model("dissolution", args.policy)
            with _control(mode):
                elapsed, result = _timed(
                    lambda: model.dwba.evaluate(*_rays(args.points)))
            actual = _values(result)
            if expected is None:
                expected = actual
            maximum[mode] = max(maximum[mode], _check(actual, expected))
            samples[mode].append(elapsed)
            print("fresh", trial, mode, elapsed, flush=True)
    for mode in modes:
        model = _model("dissolution", args.policy)
        stages = {}
        with _control(mode):
            start = time.perf_counter()
            with model.dwba.batch():
                entry = time.perf_counter() - start
                for label in ("first", "repeat"):
                    elapsed, result = _timed(
                        lambda: model.dwba.evaluate(*_rays(args.points)))
                    maximum[mode] = max(
                        maximum[mode], _check(_values(result), expected))
                    stages[label] = elapsed
        report["modes"][mode] = dict(fresh=_summary(samples[mode]),
                                    batch_entry_seconds=entry, batch_seconds=stages,
                                    max_relative_difference=maximum[mode],
                                    counts=_counts(args.policy, mode))
        save()
        print("batch and counts", mode, "completed", flush=True)
    report["edge_checks"] = _edge_checks(args.policy, modes)
    save()
    print("edge checks", report["edge_checks"], flush=True)


if __name__ == "__main__":
    main()

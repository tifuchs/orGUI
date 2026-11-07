"""Controlled C1 throughput audit using the supplied, hash-checked native loader.

Run from the checkout with --request-dir pointing to installed_api_verification
and --output pointing to a new JSON file. Optional --native loads a saved native
extension for an A/B run without replacing the installation. All model edits
are small in-memory proposals outside batches; source data is read-only.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np  # noqa: E402 -- thread limits must precede NumPy import
import psutil  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _import(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _setup(request, native, ensemble=None):
    global _loader, _root
    _loader = _import(Path(request) / "benchmark_current_model.py", "c1_loader")
    _root = Path(request).resolve().parents[3]
    # Later performance baselines also change DWBAState. Restore that complete
    # module before importing the matching ensemble, rather than timing a mix.
    if ensemble and (Path(ensemble).parent / "CTRdwba.py").exists():
        from orgui.datautils import xrayutils
        saved_dwba = _import(
            Path(ensemble).parent / "CTRdwba.py", "orgui.datautils.xrayutils.CTRdwba")
        xrayutils.CTRdwba = saved_dwba
    if ensemble:
        from orgui.datautils import xrayutils
        from orgui.datautils.xrayutils import CTRincoherent
        saved = _import(
            ensemble, "orgui.datautils.xrayutils._CTRincoherent_dwba")
        xrayutils._CTRincoherent_dwba = saved
        CTRincoherent.SurfaceHeightDWBA = saved.SurfaceHeightDWBA
        CTRincoherent.DWBAEnsembleResult = saved.DWBAEnsembleResult
    if native:
        from orgui.datautils.xrayutils import CTRdwba, CTRoptics, CTRuc
        binary = _import(native, "_CTRcalc_cpp")
        for module in (CTRdwba, CTRoptics, CTRuc):
            module._CTRcalc_cpp = binary


def _load(policy="state_specific"):
    return _loader.load_model(_root, policy)


def _rays():
    return np.ones(512), -np.ones(512), np.linspace(2.4, 2.8, 512)


def _optical(result):
    return np.array([result.reflectivity, result.coherent_intensity,
                     result.incoherent_intensity, abs(result.amplitude_prefactor)**2])


def _proposal(model, kind, index):
    surface = model.target_surface()
    if kind == "geometry":
        model.coherent_model.dwba.set_ctr_geometry(
            alpha_i=np.deg2rad(.34 + (index + 1) * 1e-5))
    elif kind == "dwf":
        surface.unitcell.basis_0[:, 4:6] += 1e-4
        surface.setFitParameters([])
    elif kind == "occupancy":
        # A common coverage scale preserves the saved oxygen/cation ratios.
        surface.unitcell.basis_0[:, 6] *= .9999
        surface.setFitParameters([])
    elif kind == "profile":
        # Decrease variance from the upper finite-support bound, never widen it.
        surface.basis_0[2] *= .9999
        surface.setFitParameters([])
    model.validate(forward_model="dwba")


def _job(index):
    start = time.perf_counter()
    model = _load()
    setup = time.perf_counter() - start
    _proposal(model, "dwf", index)
    start = time.perf_counter()
    result = model.dwba.evaluate(*_rays())
    elapsed = time.perf_counter() - start
    memory = psutil.Process().memory_info()
    return dict(setup_seconds=setup, evaluation_seconds=elapsed,
                peak_rss_bytes=getattr(memory, "peak_wset", memory.rss),
                checksum=float(np.sum(result.reflectivity)))


def _sample():
    old = np.load(_root / "data_out/1V2/0056/diagnostics/"
                  "installed_api_verification/batch_999.npz")
    metadata = json.loads((_root / "settings/1V2/0056_C1.json").read_text())
    model = _load()
    values = np.zeros((len(old["coordinates"]), 3))
    controls = []
    start = time.perf_counter()
    for rotation in metadata["domain_rotations"]:
        crystal = model.coherent_model
        crystal.domains = [[(np.array(rotation), 1.)]
                           for _ in crystal.uc_surface_list]
        crystal.dwba.clear_cache()
        domain = []
        with model.dwba.batch():
            for outgoing in ("s", "p"):
                for attenuation in (0., metadata["kinematic_attenuation"]):
                    result = model.dwba.evaluate(
                        *old["coordinates"].T, polarization_f=outgoing,
                        bulk_attenuation=attenuation)
                    domain.append(_optical(result))
        controls.append(domain)
    elapsed = time.perf_counter() - start
    controls = np.asarray(controls)
    # Same intensity order: equal four-domain average, then sum exit channels.
    for extra, columns in ((0, [1, 3, 4]), (1, [2, 7, 8])):
        values[:] = controls[:, extra::2, :3, :].sum(axis=(0, 1)).T / 4
        np.testing.assert_allclose(values, old["values"][:, columns],
                                   rtol=1e-9, atol=1e-28)
    return elapsed, controls


def main():
    """Measure complete native calls, legal proposals and process throughput."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native", type=Path)
    parser.add_argument("--ensemble", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--workers", type=int, nargs="*", choices=range(1, 9),
                        default=[1, 2, 4, 8])
    parser.add_argument("--compare", type=Path)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if args.output.exists():
        parser.error("Choose a new output file")
    native = str(args.native.resolve()) if args.native else None
    ensemble = str(args.ensemble.resolve()) if args.ensemble else None
    _setup(args.request_dir, native, ensemble)
    from orgui.datautils.xrayutils import CTRdwba, _CTRincoherent_dwba
    sources = [CTRdwba.__file__, _CTRincoherent_dwba.__file__,
               CTRdwba._CTRcalc_cpp.__file__]
    report = dict(sources={str(p): hashlib.sha256(Path(p).read_bytes()).hexdigest()
                           for p in sources}, points=512, repeats=args.repeats,
                  threads=1, policies={}, proposals={}, scaling={})
    arrays = {}

    def _save():
        args.output.parent.mkdir(parents=True, exist_ok=True)
        memory = psutil.Process().memory_info()
        report["parent_peak_rss_bytes"] = getattr(memory, "peak_wset", memory.rss)
        args.output.write_text(json.dumps(report, indent=2))
        np.savez(args.output.with_suffix(".npz"), **arrays)

    for policy in ("shared", "state_specific"):
        times, setups = [], []
        for trial in range(args.repeats):
            start = time.perf_counter()
            model = _load(policy)
            setups.append(time.perf_counter() - start)
            start = time.perf_counter()
            result = model.dwba.evaluate(*_rays())
            times.append(time.perf_counter() - start)
            arrays[f"{policy}_fresh_{trial}"] = _optical(result)
        model = _load(policy)
        batch = []
        with model.dwba.batch():
            for trial in range(2):
                start = time.perf_counter()
                result = model.dwba.evaluate(*_rays())
                batch.append(time.perf_counter() - start)
                arrays[f"{policy}_batch_{trial}"] = _optical(result)
        np.testing.assert_array_equal(arrays[f"{policy}_batch_0"],
                                      arrays[f"{policy}_batch_1"])
        report["policies"][policy] = dict(setup_seconds=setups, fresh_seconds=times,
                                          batch_seconds=batch)
        print(policy, report["policies"][policy], flush=True)
        _save()
    for kind in ("unchanged", "geometry", "dwf", "occupancy", "profile"):
        model = _load()
        times = []
        for trial in range(args.repeats):
            _proposal(model, kind, trial)
            start = time.perf_counter()
            result = model.dwba.evaluate(*_rays())
            times.append(time.perf_counter() - start)
            arrays[f"{kind}_{trial}"] = _optical(result)
        report["proposals"][kind] = times
        print(kind, times, flush=True)
        _save()
    samples = []
    for trial in range(args.repeats):
        elapsed, arrays[f"sample_{trial}"] = _sample()
        samples.append(elapsed)
    report["sample_seconds"] = samples
    print("sample", samples, flush=True)
    _save()
    for workers in args.workers:
        available = psutil.virtual_memory().available
        start = time.perf_counter()
        with ProcessPoolExecutor(max_workers=workers, initializer=_setup,
                                 initargs=(args.request_dir, native, ensemble)) as pool:
            jobs = list(pool.map(_job, range(8)))
        report["scaling"][workers] = dict(
            wall_seconds=time.perf_counter() - start, jobs=jobs,
            available_ram_at_start_bytes=available,
            peak_worker_rss_sum_bound=sum(sorted(
                (job["peak_rss_bytes"] for job in jobs), reverse=True)[:workers]))
        print("workers", workers, report["scaling"][workers], flush=True)
        _save()
    if args.compare:
        old = np.load(args.compare.with_suffix(".npz"))
        errors = {}
        for key, array in arrays.items():
            np.testing.assert_allclose(array, old[key], rtol=1e-9, atol=1e-28)
            errors[key] = float(np.max(abs(array - old[key]) /
                                       np.maximum(abs(old[key]), 1e-300)))
        report["comparison_max_relative_error"] = max(errors.values())
    report["status"] = "passed"
    _save()


if __name__ == "__main__":
    main()

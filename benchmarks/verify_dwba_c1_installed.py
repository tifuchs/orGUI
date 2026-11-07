"""Run the supplied C1 verification driver against an asserted wheel install.

This launcher changes import routing only. The supplied driver's model loading,
scientific calculation, tolerances and benchmark workloads run unchanged.
"""

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import runpy
import sys


def main():
    """Assert imported paths/hashes before and after the installed verification."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-packages", type=Path, required=True)
    parser.add_argument("--request-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ensemble-sha256", required=True)
    args = parser.parse_args()
    site = args.site_packages.resolve()
    request = args.request_dir.resolve()
    output = args.output.resolve()
    if output.exists():
        parser.error("Choose a new output directory; checkpoints must not be reused")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    sys.path.insert(0, str(site))
    modules = [importlib.import_module("orgui.datautils.xrayutils." + name)
               for name in ("CTRdwba", "_CTRincoherent_dwba", "_CTRcalc_cpp")]
    sources = {}
    for module in modules:
        path = Path(module.__file__).resolve()
        assert path.is_relative_to(site), path
        sources[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    ensemble_hash = sources[str(Path(modules[1].__file__).resolve())]
    assert ensemble_hash == args.ensemble_sha256, ensemble_hash
    audit = dict(version=importlib.metadata.version("orgui"), sources=sources)
    output.mkdir(parents=True)
    (output / "import_audit.json").write_text(json.dumps(audit, indent=2))
    print("ASSERTED_INSTALL", json.dumps(audit), flush=True)
    sys.argv = [str(request / "verify.py"), str(request.parents[3]), str(output)]
    runpy.run_path(sys.argv[0], run_name="__main__")
    report = json.loads((output / "verification.json").read_text())
    assert report["native_module_sha256"] == ensemble_hash
    assert report["status"] == "passed"


if __name__ == "__main__":
    main()

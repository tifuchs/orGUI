"""Check imports needed by the Read the Docs autodoc build."""

import importlib
import os
import sys
import traceback

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, REPO_ROOT)

AUTODOC_MODULES = [
    "orgui.main",
    "orgui.app.orGUI",
]


def _main():
    if sys.platform.startswith("linux"):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

    failed = False
    for module_name in AUTODOC_MODULES:
        print(f"Checking import: {module_name}", flush=True)
        try:
            importlib.import_module(module_name)
        except Exception:
            failed = True
            print(f"FAILED import: {module_name}", file=sys.stderr)
            traceback.print_exc()
        else:
            print(f"OK import: {module_name}", flush=True)

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(_main())

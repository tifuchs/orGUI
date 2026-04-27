# AGENTS.md

## Purpose

This repository contains `orGUI`, a Python GUI and analysis toolkit for surface X-ray diffraction data from large 2D detectors. It combines:

- reciprocal-space geometry and diffractometer calculations
- detector calibration and detector-to-angle / detector-to-Q transforms
- beamline and file-format loading
- GUI workflows for orientation finding and ROI-based integration

For this project, scientific correctness is more important than stylistic cleanup. Prefer small, well-verified changes over broad refactors. Units matter, but they are one part of the broader correctness story alongside geometry conventions, detector behavior, file-format compatibility, and regression coverage.

## Repository map

- `orgui/datautils/xrayutils/`
  Core scientific code. This is the highest-risk area.
  - `HKLVlieg.py`: lattice, UB, diffractometer geometry, angle transforms
  - `DetectorCalibration.py`: detector geometry, pixel/angle/Q conversions, pyFAI integration
  - `ReciprocalNavigation.py`: allowed reflections / CTR navigation helpers
  - `CTRcalc.py`: main CTR crystal entry point; defines `SXRDCrystal`, imports and connects the lower-level CTR stack, and is the usual entry point into `CTRuc.py` and `CTRfilm.py`
  - `CTRuc.py`: unit-cell and water-model definitions; atomic positions, structure factors, density helpers, fit parameters
  - `CTRfilm.py`: higher-level film and epitaxy objects built on top of unit cells, including stacked/interface models such as `Film` and `EpitaxyInterface`
  - `CTRutil.py`: shared CTR parameter utilities and parsing helpers used by the CTR modules
  - `CTRplotutil.py`: CTR data containers, plotting helpers, and ANAROD-style import/export paths used in tests and analysis workflows
  - `_CTRcalc_accel.py`: optional numba-accelerated kernels for selected CTR calculations
  - `element_data.py`, `unitcells/`: scattering-factor data and bundled structure/unit-cell reference files
- `orgui/app/`
  GUI components and user workflows.
  - `QUBCalculator.py`: crystal + machine parameter UI, source/detector conversions, config integration
  - `database.py`: config loading, detector/crystal object setup, Nexus/HDF5 handling
  - `QScanSelector.py`, `ROIutils.py`, `orGUI.py`: integration and image workflows
- `orgui/backend/`
  Scan loading abstractions and beamline-specific adapters.
- `examples/`
  Best source of current config conventions and user-facing units.
- `orgui/datautils/xrayutils/test/`
  Regression tests for lattice math, detector conversions, and CTR calculations.

Within the CTR stack, the typical dependency direction is:

- `CTRcalc.py` -> `CTRuc.py` and `CTRfilm.py`
- `CTRfilm.py` -> `CTRuc.py`
- `CTRplotutil.py` consumes CTR objects for analysis, plotting, and test/reference comparisons

## Working principles

- Do not change scientific conventions silently.
- Treat explicit scale factors and coordinate transforms with care until you confirm the source and target meaning.
- Keep user-facing scientific behavior stable unless the task explicitly asks for a behavior change.
- When changing geometry or conversion code, preserve backwards compatibility for existing config files and stored calibration data if possible.
- Prefer adding comments and docstrings that state units at function boundaries over broad rewrites.
- Follow PEP 8 for new and modified code where practical, but do not force broad style churn in stable scientific code just to satisfy formatting preferences.

## Unit conventions

Unit documentation is currently incomplete in the codebase. When adding or modifying code, make units explicit where they help avoid ambiguity, especially at config, UI, file-format, and library boundaries.

### Default unit list

Unless a function or file format explicitly requires something else, use these defaults in new documentation and new user-facing labels:

- momentum transfer `Q`: `Angstrom^-1`
- reciprocal lattice vector lengths: `Angstrom^-1`
- reciprocal lattice coordinates `(h, k, l)`: `r.l.u.`
- wavelength: `Angstrom`
- X-ray energy: `keV`
- direct lattice constants `a, b, c`: `Angstrom`
- lattice angles `alpha, beta, gamma`: `deg`
- diffractometer angles (`mu`, `omega`, `delta`, `gamma`, `chi`, `phi`, azimuth, polarization axis): `deg` in config/UI, `rad` in internal geometry/math unless documented otherwise
- sample-detector distance `SDD`: `m` in config, often converted internally for Fit2D / pyFAI-facing APIs
- detector pixel size: `m` in config, often converted to `um` for Fit2D-style APIs and display
- detector center coordinates `cpx`, `cpy`: pixels
- image / ROI coordinates: pixels
- electron density: `electrons / Angstrom^3`
- structure factor / intensity: arbitrary units unless a calibrated scale is explicitly introduced

These defaults are guidance, not a mandate to rename established variables or rewrite stable scientific code just to restate units.

### Known conversion hotspots

Be especially careful in these files:

- `orgui/app/QUBCalculator.py`
  - energy `E` is user-facing in `keV`
  - wavelength is user-facing in `Angstrom`
  - the pyFAI geometry model wavelength is written in meters via `wavelength * 1e-10`
- `orgui/app/database.py`
  - config `SDD` is read in meters, then passed to Fit2D-style APIs in millimeters
  - config `pixelsize` is read in meters, then passed to Fit2D-style APIs in micrometers
  - detector wavelength is stored in meters, while UB calculator wavelength is in Angstrom
- `orgui/datautils/xrayutils/HKLVlieg.py`
  - lattice lengths are in Angstrom
  - reciprocal vectors include the `2*pi` factor
  - many angular APIs work in radians internally
- `orgui/datautils/xrayutils/DetectorCalibration.py`
  - pyFAI geometry fields use pyFAI conventions
  - some helpers export to mm or nm for external tools

If you touch one side of a conversion, inspect the matching inverse conversion in the same workflow when one exists.

## Preferred documentation pattern

Documentation is planned to use Sphinx with Markedly Structured Text. New docstrings and docstring edits should therefore use explicit Sphinx-friendly Markedly Structured Text formatting instead of NumPy-style sections.

Every new public function should now have a docstring. In this repository, every function whose name does not start with `_` is considered public.

For any scientific function that crosses a unit boundary, prefer docstrings like:

```python
def set_Xray_source(source_config):
    """Set the X-ray source.

    :param dict source_config:
        Either ``{"E": ...}`` in keV or ``{"wavelength": ...}`` in Angstrom.
    :raises ValueError:
        If neither energy nor wavelength is provided.
    """
```

For GUI-only public functions, keep the normal summary line for Sphinx output and add the GUI-only marker as a nearby code comment. If useful for generated documentation, also add a Sphinx note:

```python
# GUI-only: user-triggered dialog path. Must not run in CLI or batch mode.
def open_detector_dialog(self):
    """Open the detector selection dialog.

    .. note::
       GUI-only. This path may create blocking dialogs and must not be
       called from CLI, batch, or other shared non-interactive code.
    """
```

Explanation:

- use the `# GUI-only: ...` comment as the canonical code-review marker
- keep the first docstring line as a normal summary so Sphinx renders useful API docs
- use `.. note::` only when the GUI-only restriction should also appear in generated documentation

Good additions:

- mention whether angles are `deg` or `rad`
- mention whether `Q` is `Angstrom^-1` or `r.l.u.`
- mention whether a value is "user-facing/config/UI" or "internal/pyFAI/Fit2D"
- use explicit fields such as `:param:`, `:type:`, `:returns:`, `:rtype:`, and `:raises:` where helpful

## Change guidance by area

### Geometry, lattice, reciprocal-space math

- Treat `HKLVlieg.py`, `ReciprocalNavigation.py`, and `CTR*` modules as scientific core code.
- Preserve numerical conventions unless the change is deliberate and tested.
- Avoid changing sign conventions, axis ordering, or inclusion of `2*pi` without tracing downstream use.
- Review geometry assumptions, not only units.
- Add or update regression tests for any change in these files.

### Detector calibration and image geometry

- Confirm pyFAI / Fit2D unit expectations before editing conversion code.
- Preserve round-trip behavior for pixel-to-angle and angle-to-pixel functions.
- If adding a new exported quantity, document its unit at the point of creation.

### GUI code

- UI labels and tooltips should display units whenever a numerical field would be ambiguous.
- If a GUI action edits scientific parameters, trace through to the backend object to confirm unit consistency.
- Avoid changing GUI behavior and numerical behavior in the same patch unless necessary.
- Mark functions, methods, or local branches that are intentionally interactive as `GUI-only`.
- Unmarked code must be treated as shared code and therefore remain safe for CLI/headless execution.
- `QMessageBox`, `QProgressDialog`, and other modal blocking UI are allowed only in code explicitly marked `GUI-only`.
- Prefer a nearby comment such as `# GUI-only: user-triggered dialog path.`
- If the constraint should also appear in generated documentation, use a Sphinx-friendly docstring note such as `.. note:: GUI-only. Must not be called from CLI or batch mode.`

## Startup and logging model

The repository has two startup modes selected in `orgui/main.py`:

- GUI mode
  Starts a normal Qt application and calls `logger_utils.set_logging_context("gui")` before launching the main window.
- CLI / headless mode
  Calls `logger_utils.set_logging_context("cli")`, sets `QT_QPA_PLATFORM=minimal`, and still creates a Qt application plus an `orGUI` main-window object for scripting access.

This means CLI mode is not "Qt-free". Shared code can still run inside a Qt application even when no GUI should block the session.

### Logging behavior

The special handler logic lives in `orgui/logger_utils.py`.

- In GUI context:
  - `MessageBoxHandler` is attached to the `orgui` logger hierarchy.
  - A dialog is shown only when the log record has `extra={"show_dialog": True}`.
  - `dialog_level` controls whether the dialog is shown as information, warning, or error.
- In CLI context:
  - `CLIExceptionHandler` is attached to the same logger hierarchy.
  - Any `ERROR` or higher record becomes an exception.
  - `show_dialog` is ignored in CLI mode.
  - `WARNING` and `INFO` records stay non-fatal log output.

### Safe logging rules

- For shared code that may run in either GUI or CLI mode, do not create direct `QMessageBox` instances or other modal blocking dialogs. In such paths, blocking UI in CLI mode is a bug. Use logging and the context-aware handlers in `logger_utils` instead.
- Use `extra={"show_dialog": True, "parent": ..., "title": ..., "description": ...}` only when a GUI dialog is explicitly desired.
- If a condition should remain non-fatal in CLI mode, log it at `WARNING` or lower.
- If a condition should abort CLI execution, log it at `ERROR` or use `logger.exception(...)`.
- Do not rely on `dialog_level=logging.WARNING` to make an `ERROR` record non-fatal in CLI mode; CLI behavior is driven by the record level, not the dialog level.
- Be careful with `logger.exception(...)` in shared code: it logs at `ERROR` and therefore raises in CLI mode.

### Dialog safety

Blocking dialogs in headless or automation paths are a fatal risk because CLI mode still has a Qt event environment but must not block on interactive UI.

- Avoid direct `qt.QMessageBox` calls in code paths that may run from CLI scripts, batch mode, worker threads, or non-interactive flows.
- Direct message boxes are acceptable in clearly GUI-only interaction paths where user confirmation is the whole point.
- When in doubt, route user-visible failure reporting through the logger handler instead of constructing dialogs directly.
- Never create modal dialogs from worker threads.

### Progress reporting

- Prefer `logger_utils.create_progress_logger(...)` for new long-running work.
- In GUI mode it creates a `QProgressDialog`.
- In CLI mode it creates a logging-based progress reporter that emits `PROGRESS:...` messages instead of opening UI.
- Avoid raw `QProgressDialog` in shared or batch-capable paths; it is only appropriate for clearly GUI-only workflows.
- Update and finish progress objects from the controlling thread, not from worker threads.

### File loading / backends

- New beamline loaders should normalize loaded metadata into the project's existing conventions rather than introducing backend-specific units deep in the app.
- Document any assumptions about angle direction, detector orientation, or energy units near the loader.

## Validation expectations

For changes in scientific core code, run the narrowest relevant tests first, then broaden if needed:

- `ruff check orgui`
- `pytest orgui/datautils/xrayutils/test/test_HKLcalc.py`
- `pytest orgui/datautils/xrayutils/test/test_DetectorCalibration.py`
- `pytest orgui/datautils/xrayutils/test/test_CTRcalc.py`

When changing config-loading or unit-conversion code, also inspect:

- `examples/config_minimal`
- other `examples/config_*` files

If a change cannot be covered by an existing test, add at least one regression test or a precise docstring/comment explaining the invariant being preserved.

Lint findings should be fixed when they are local to the change, but lint cleanup must not drive risky scientific refactors.

## Safe change strategy

When working on a bug or feature:

1. Identify whether the values involved are user-facing, config-facing, internal scientific values, or external library values.
2. Write down the relevant conventions at each boundary before changing code: units, coordinate frame, angle convention, and any library-specific expectations.
3. Verify both forward and inverse transforms if they exist.
4. Prefer localized fixes over repository-wide "cleanup" in numerical code.
5. Update docs, labels, or comments when a convention was previously implicit.

## Things to avoid

- Do not replace explicit conversion factors with "cleaner" code unless the replacement is unit-audited.
- Do not rename variables in scientific code if that makes established literature conventions harder to follow.
- Do not remove old config behavior casually; many users likely depend on existing beamline and config conventions.
- Do not prioritize lint-only edits over preserving readable scientific intent.

## Good first documentation targets

If the task is primarily documentation-oriented, these are high-value places to improve:

- config file field units in `examples/config_*`
- docstrings for geometry conversion methods in `HKLVlieg.py`
- docstrings for detector conversion methods in `DetectorCalibration.py`
- user-facing tooltips and labels in `QUBCalculator.py`

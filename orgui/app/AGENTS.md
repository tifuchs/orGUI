# AGENTS.md

## Scope

This directory contains GUI components and user workflows:

- `QUBCalculator.py`: crystal and machine parameter UI, source/detector
  conversions, config integration.
- `database.py`: config loading, detector/crystal setup, Nexus/HDF5 handling.
- `QScanSelector.py`, `ROIutils.py`, `orGUI.py`: integration and image
  workflows.
- Dialogs and tests under this directory cover user-facing behavior and GUI
  regressions.

Use the repository root instructions together with this file.

## GUI And Shared-Code Boundary

CLI/headless mode still creates a Qt application and an `orGUI` main-window
object for scripting access. Code in this directory is not automatically
GUI-only just because it imports Qt.

- Mark intentionally interactive functions, methods, or branches with a nearby
  comment such as `# GUI-only: user-triggered dialog path.`
- Unmarked code must be safe for CLI/headless execution.
- Direct `QMessageBox`, `QProgressDialog`, and other modal blocking UI are
  acceptable only in clearly marked GUI-only paths.
- Do not create modal dialogs from worker threads.
- For shared paths, report user-visible failures through `logger_utils` instead
  of constructing dialogs directly.

For GUI-only public functions, keep a normal Sphinx summary and add a note only
when the restriction should appear in generated docs:

```python
# GUI-only: user-triggered dialog path. Must not run in CLI or batch mode.
def open_detector_dialog(self):
    """Open the detector selection dialog.

    .. note::
       GUI-only. This path may create blocking dialogs and must not be called
       from CLI, batch, or other shared non-interactive code.
    """
```

## Logging And Progress

The startup and logging model is defined in `orgui/main.py` and
`orgui/logger_utils.py`.

- GUI mode attaches `MessageBoxHandler`; a dialog is shown only when a log
  record has `extra={"show_dialog": True}`.
- CLI mode attaches `CLIExceptionHandler`; any `ERROR` or higher record becomes
  an exception, and `show_dialog` is ignored.
- If a shared-code condition should remain non-fatal in CLI mode, log it at
  `WARNING` or lower.
- If a condition should abort CLI execution, log it at `ERROR` or use
  `logger.exception(...)`.
- Prefer `logger_utils.create_progress_logger(...)` for new long-running work.
  It creates a `QProgressDialog` in GUI mode and logging-based `PROGRESS:...`
  messages in CLI mode.

## Unit And Config Boundaries

UI labels and tooltips should show units whenever a numeric field would be
ambiguous. If a GUI action edits scientific parameters, trace through to the
backend object and confirm unit consistency.

Known hotspots:

- `QUBCalculator.py`: energy `E` is user-facing in keV; wavelength is
  user-facing in Angstrom; pyFAI wavelength is written in meters via
  `wavelength * 1e-10`.
- `database.py`: config `SDD` is read in meters, then passed to Fit2D-style APIs
  in millimeters.
- `database.py`: config `pixelsize` is read in meters, then passed to
  Fit2D-style APIs in micrometers.
- `database.py`: detector wavelength is stored in meters, while UB calculator
  wavelength is in Angstrom.
- Image and ROI coordinates are pixels.

If you touch one side of a conversion, inspect the matching inverse conversion
in the same workflow when one exists. Avoid changing GUI behavior and numerical
behavior in the same patch unless necessary.

## Backend And Example Interactions

New or changed beamline/file loading should normalize metadata into the
project's existing conventions rather than introducing backend-specific units
deep in the app.

Document assumptions about angle direction, detector orientation, or energy
units near the loader or config boundary. For config changes, inspect
`examples/config_minimal` and the other `examples/config_*` files.

## Tests

Choose checks based on the touched behavior:

- `pytest orgui/app/test`
- `ruff check orgui/app`
- For config or unit conversions, inspect `examples/config_*` and consider the
  detector/geometry tests under `orgui/datautils/xrayutils/test`.

Add a focused regression test when a GUI workflow or config conversion bug can
be reproduced without broad UI automation.

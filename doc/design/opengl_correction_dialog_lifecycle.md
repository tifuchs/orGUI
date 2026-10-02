# OpenGL correction-dialog lifecycle investigation

## Reopening the options dialog

The reported failure followed loading an HDF5 configuration with OpenGL
plotting enabled. The second opening of the correction-options window was
blank, with Qt reporting a non-OpenGL surface and failure to create its RHI.

A native Windows reproduction created an OpenGL main plot, showed the raster
options dialog, and then created the OpenGL beam previews in its child dialog.
Reopening the options dialog reproduced those Qt messages. Giving the
top-level beam dialog `WA_NativeWindow` isolates its preview hierarchy from
the options window's surface selection. Embedded horizontal-profile controls
stay ordinary child widgets. The previews retain the selected silx backend.

The regression in `test_integration_options_dialog.py` now covers an HDF5
roundtrip, an existing main plot, all four previews, colormap editing,
reopening, preservation of settings and controls, and deferred destruction.
The native OpenGL case verifies valid rendering contexts; it is skipped on
headless Qt platforms.

The separate pyFAI warning came from nested NeXus attributes in detector
configuration dictionaries. Recursive attribute removal before detector
construction preserves calibration values and avoids passing those attributes
to the detector factory.

## Capture-triggered native fault

Investigative processes also caused Windows access-violation popups. These
are a separate failure from the raster/OpenGL surface mismatch.

In the installed environment (Python 3.14.7, PySide6 6.11.2, silx 3.2.0a0),
opening the diagnostic footprint colormap editor, hiding it, and calling
`beam.grab().save(...)` reproducibly corrupts a reference. The eventual crash
occurs during Python garbage collection. Without the beam capture, the same
three-cycle workflow using the reported HDF5 configuration, colormap editor,
OpenGL previews, deferred window deletion, and explicit garbage collection
exits successfully. Capturing only the options window also succeeds.

The native debugger identified the stale entry as
`silx.gui.plot.ColorBar.ColorScaleBar.resizeEvent`. A hardware watchpoint
caught its reference count changing from one to zero in
`shiboken6!Sbk_GetPyOverride`, while `QWidget::grab` was recursively sending
pending resize events. The function remained referenced by its class
dictionary; the later `dict_traverse` access violation is the consequence.

PySide's [6.11.2 override-dispatch source](https://github.com/pyside/pyside-setup/blob/v6.11.2/sources/shiboken6/libshiboken/basewrapper.cpp)
decrements the override reference on its pending-error path.
[Override lookup](https://github.com/pyside/pyside-setup/blob/v6.11.2/sources/shiboken6/libshiboken/bindingmanager.cpp)
returns the function from the bound method. The watchpoint and these sources
indicate a binding reference-ownership defect. The preceding error that takes
this dispatch path has not yet been identified. Basic standalone silx plots
and standalone colorbars did not reproduce the fault, so those simpler cases
do not establish an independent upstream reproducer.

No library-memory patch, reference-retention workaround, backend substitution,
or environment upgrade has been applied. Native investigative processes use
process-local Windows error-mode and Windows Error Reporting flags to suppress
crash dialogs. Capture of the beam widget remains excluded from the successful
workflow checks; its failure is unresolved.

## Validation

- Real stored configuration: three native OpenGL reopen cycles, including
  colormap use and explicit cleanup, passed without capture or Qt RHI errors.
- Expanded native HDF5/backend regression: two cases passed.
- Relevant GUI/config and detector checks: 166 passed, one headless OpenGL
  skip, and 46 subtests passed. Existing numerical and dependency warnings
  remain.
- Ruff on the modified Python files and `git diff --check` passed.

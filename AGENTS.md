# AGENTS.md

## Purpose

This repository contains `orGUI`, a Python GUI and analysis toolkit for surface X-ray diffraction data from large 2D detectors. It combines:

- reciprocal-space geometry and diffractometer calculations
- detector calibration and detector-to-angle / detector-to-Q transforms
- beamline and file-format loading
- GUI workflows for orientation finding and ROI-based integration

For this project, scientific correctness is more important than stylistic cleanup. Prefer small, well-verified changes over broad refactors.

## Repository map

- `orgui/datautils/xrayutils/`
  Core scientific code. This is the highest-risk area.
  - `HKLVlieg.py`: lattice, UB, diffractometer geometry, angle transforms
  - `DetectorCalibration.py`: detector geometry, pixel/angle/Q conversions, pyFAI integration
  - `ReciprocalNavigation.py`: allowed reflections / CTR navigation helpers
  - `CTRcalc.py`, `CTRuc.py`, `CTRfilm.py`, `CTRplotutil.py`: CTR and crystal structure calculations
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

## Working principles

- Do not change units silently.
- Do not "simplify" code that contains explicit scale factors until you confirm the source and target units.
- Keep user-facing scientific behavior stable unless the task explicitly asks for a behavior change.
- When changing geometry or conversion code, preserve backwards compatibility for existing config files and stored calibration data if possible.
- Prefer adding comments and docstrings that state units at function boundaries over broad rewrites.

## Unit conventions

Unit documentation is currently incomplete in the codebase. When adding or modifying code, make units explicit in docstrings, comments, variable names, and UI labels.

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

If you touch one side of a conversion, inspect the matching inverse conversion in the same workflow.

## Preferred documentation pattern

Documentation is planned to use Sphinx with Markedly Structured Text. New docstrings and docstring edits should therefore use explicit Sphinx-friendly Markedly Structured Text formatting instead of NumPy-style sections.

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
- Add or update regression tests for any change in these files.

### Detector calibration and image geometry

- Confirm pyFAI / Fit2D unit expectations before editing conversion code.
- Preserve round-trip behavior for pixel-to-angle and angle-to-pixel functions.
- If adding a new exported quantity, document its unit at the point of creation.

### GUI code

- UI labels and tooltips should display units whenever a numerical field would be ambiguous.
- If a GUI action edits scientific parameters, trace through to the backend object to confirm unit consistency.
- Avoid changing GUI behavior and numerical behavior in the same patch unless necessary.

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
2. Write down the units at each boundary before changing code.
3. Verify both forward and inverse transforms if they exist.
4. Prefer localized fixes over repository-wide "cleanup" in numerical code.
5. Update docs, labels, or comments when a unit was previously implicit.

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

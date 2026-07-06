# AGENTS.md

## Scope

This directory contains the highest-risk scientific code:

- `HKLVlieg.py`: lattice, UB, diffractometer geometry, angle transforms.
- `DetectorCalibration.py`: detector geometry, pixel/angle/Q conversions,
  pyFAI integration.
- `ReciprocalNavigation.py`: allowed reflections and CTR navigation helpers.
- `CTRcalc.py`: main CTR crystal entry point; defines `SXRDCrystal`.
- `CTRuc.py`: unit-cell and water-model definitions, atomic positions,
  structure factors, density helpers, fit parameters.
- `CTRfilm.py`: film and epitaxy objects built on top of unit cells.
- `CTRutil.py`, `CTRplotutil.py`: CTR parsing, data containers, plotting, and
  ANAROD-style import/export paths.
- `_CTRcalc_accel.py`: optional numba-accelerated kernels.
- `element_data.py`, `unitcells/`: scattering-factor data and bundled reference
  structure files.
- `test/`: regression tests for lattice math, detector conversions, and CTR
  calculations.

Use the repository root instructions together with this file.

## Scientific Change Workflow

Before changing geometry, detector, reciprocal-space, or CTR behavior:

1. Identify whether each value is user-facing, config-facing, internal
   scientific state, or external-library state.
2. Write down the relevant units, coordinate frame, angle convention, and
   library-specific expectation at each boundary.
3. Verify forward and inverse transforms when both exist.
4. Keep the fix localized. Do not combine numerical behavior changes with broad
   cleanup.
5. Add or update a regression test when behavior changes.

Preserve numerical conventions unless the change is deliberate and tested.
Avoid changing sign conventions, axis ordering, or inclusion of `2*pi` without
tracing downstream use.

## Unit Conventions

Unless a function or file format explicitly requires something else, use these
defaults in new documentation and user-facing labels:

- momentum transfer `Q`: `Angstrom^-1`
- reciprocal lattice vector lengths: `Angstrom^-1`
- reciprocal lattice coordinates `(h, k, l)`: `r.l.u.`
- wavelength: `Angstrom`
- X-ray energy: `keV`
- direct lattice constants `a, b, c`: `Angstrom`
- lattice angles `alpha`, `beta`, `gamma`: `deg`
- diffractometer angles: `deg` in config/UI, `rad` in internal geometry/math
  unless documented otherwise
- sample-detector distance `SDD`: `m` in config, often converted internally for
  Fit2D/pyFAI-facing APIs
- detector pixel size: `m` in config, often converted to `um` for Fit2D-style
  APIs and display
- detector center coordinates `cpx`, `cpy`: pixels
- image and ROI coordinates: pixels
- electron density: `electrons / Angstrom^3`
- structure factor and intensity: arbitrary units unless a calibrated scale is
  explicitly introduced

These defaults are guidance, not a mandate to rename established variables or
rewrite stable code just to restate units.

## Hotspots

- `HKLVlieg.py`: lattice lengths are in Angstrom; reciprocal vectors include
  the `2*pi` factor; many angular APIs work in radians internally.
- `DetectorCalibration.py`: pyFAI geometry fields use pyFAI conventions; some
  helpers export to mm or nm for external tools.
- Detector calibration should preserve round-trip behavior for pixel-to-angle
  and angle-to-pixel functions.
- If adding an exported quantity, document its unit at the point of creation.

Within the CTR stack, the typical dependency direction is:

- `CTRcalc.py` -> `CTRuc.py` and `CTRfilm.py`
- `CTRfilm.py` -> `CTRuc.py`
- `CTRplotutil.py` consumes CTR objects for analysis, plotting, and
  test/reference comparisons

## Documentation

For public scientific functions, document units and coordinate conventions at
the boundary. Mention whether angles are degrees or radians, whether `Q` is
`Angstrom^-1` or `r.l.u.`, and whether a value is user/config-facing or internal
to pyFAI/Fit2D.

Use Sphinx-friendly fields such as `:param:`, `:type:`, `:returns:`, `:rtype:`,
and `:raises:` where they clarify the contract.

## Tests

Run the narrowest relevant regression test first:

- `pytest orgui/datautils/xrayutils/test/test_HKLcalc.py`
- `pytest orgui/datautils/xrayutils/test/test_DetectorCalibration.py`
- `pytest orgui/datautils/xrayutils/test/test_CTRcalc.py`

Use `ruff check orgui/datautils/xrayutils` for local lint checks. If a change
touches config-facing unit conversions, also inspect `examples/config_minimal`
and the other `examples/config_*` files.

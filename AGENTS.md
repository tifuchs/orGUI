# AGENTS.md

## Purpose

`orGUI` is a Python GUI and analysis toolkit for surface X-ray
diffraction data from large 2D detectors. It combines reciprocal-space
geometry, detector calibration, beamline/file loading, and GUI workflows for
orientation finding and ROI-based integration.

Scientific correctness is more important than stylistic cleanup. Prefer small,
well-verified changes over broad refactors. Units matter, but they are only one
part of correctness alongside geometry conventions, detector behavior,
file-format compatibility, and regression coverage.

## Instruction Tree

Read the most local `AGENTS.md` for the files you touch:

- `AGENTS.md`: repository-wide rules.
- `orgui/app/AGENTS.md`: GUI workflows, config loading, dialogs, progress, and
  user-facing unit boundaries.
- `orgui/datautils/xrayutils/AGENTS.md`: scientific core, detector geometry,
  reciprocal-space math, CTR calculations, and related tests.

Do not add another nested `AGENTS.md` unless it supplies critical local guidance
that would otherwise distract unrelated work.

## Working Principles

- Do not change scientific conventions silently.
- Treat explicit scale factors and coordinate transforms as scientific behavior.
  Before editing them, identify the source and target units/conventions.
- Keep user-facing scientific behavior stable unless the task explicitly asks
  for a behavior change.
- Preserve backwards compatibility for existing config files and stored
  calibration data when possible.
- Prefer comments and docstrings that state units at function boundaries over
  broad rewrites.
- Follow PEP 8 for new and modified code where practical, but do not force
  broad style churn in stable scientific code.

## Documentation

Documentation is planned to use Sphinx with Markedly Structured Text. New
docstrings and docstring edits should use Sphinx-friendly formatting rather than
NumPy-style sections.

Every new public function should have a docstring. In this repository, every
function whose name does not start with `_` is considered public.

For scientific functions that cross a unit boundary, document the boundary:

```python
def set_Xray_source(source_config):
    """Set the X-ray source.

    :param dict source_config:
        Either ``{"E": ...}`` in keV or ``{"wavelength": ...}`` in Angstrom.
    :raises ValueError:
        If neither energy nor wavelength is provided.
    """
```

Good documentation targets:

- config field units in `examples/config_*`
- geometry conversion methods in `HKLVlieg.py`
- detector conversion methods in `DetectorCalibration.py`
- user-facing labels and tooltips in `QUBCalculator.py`

## Validation

Run the narrowest relevant checks first, then broaden when the change warrants
it:

- `ruff check orgui`
- `pytest orgui/datautils/xrayutils/test/test_HKLcalc.py`
- `pytest orgui/datautils/xrayutils/test/test_DetectorCalibration.py`
- `pytest orgui/datautils/xrayutils/test/test_CTRcalc.py`

When changing config-loading or unit-conversion behavior, inspect
`examples/config_minimal` and the other `examples/config_*` files.

If a scientific change cannot be covered by an existing test, add a focused
regression test or a precise docstring/comment explaining the invariant being
preserved. Fix lint findings local to your change, but do not let lint cleanup
drive risky scientific refactors.

## Commits

Use Conventional Commits for commit and merge messages, following the repository
convention in issue #38.

- Allowed types: `build`, `ci`, `docs`, `feat`, `fix`, `perf`, `refactor`,
  `style`, and `test`.
- Use the `phys` scope for changes that can affect measurement results,
  including counters, error propagation, physics math, or generated data.
- Mark breaking changes with `!` and include a `BREAKING CHANGE:` footer when
  old post-processing scripts or generated data compatibility may be affected.

## Avoid

- Do not replace explicit conversion factors with "cleaner" code unless the
  replacement is unit-audited.
- Do not rename variables in scientific code if that makes established
  literature conventions harder to follow.
- Do not remove old config behavior casually; many users likely depend on
  existing beamline and config conventions.
- Do not prioritize lint-only edits over preserving readable scientific intent.

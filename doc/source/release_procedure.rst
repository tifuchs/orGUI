Release Procedure
=================

Use this checklist for releases from ``master``. Versions are taken from Git
tags through ``setuptools_scm``. Tags must use the form ``vX.Y.Z``.

1. Confirm the release candidate
--------------------------------

Check that the intended release commit is clean and tested:

.. code-block:: powershell

   git status --short
   git log --oneline --decorate -5
   ruff check orgui
   pytest orgui/datautils/xrayutils/test/test_HKLcalc.py
   pytest orgui/datautils/xrayutils/test/test_DetectorCalibration.py
   pytest orgui/datautils/xrayutils/test/test_CTRcalc.py
   python -m build
   python -m twine check dist/*

Manual GUI/science smoke test:

* install from a clean environment or from the built wheel
* start the GUI and confirm clean startup
* load at least one ``examples/`` configuration
* check source/detector units: energy in keV, wavelength in Angstrom, SDD in m,
  pixel size in m, and detector center in pixels
* load a representative image or scan
* run one known UB/orientation calculation
* run one representative integration workflow
* compare CTR or structure-factor output with a previous release or reference
  when scientific code changed
* start CLI/headless mode when shared loading or logging code changed

2. Update ``CHANGELOG.md`` and write release notes
--------------------------------------------------

2.1 Automatic ``CHANGELOG`` using cliff
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Choose the bump type manually: ``patch`` for bug fixes, documentation,
packaging, and compatibility fixes; ``minor`` for backwards-compatible
features. Replace ``patch`` with ``minor`` when appropriate. A major version
change implies substantial API or other user-facing breaking changes and
requires a dedicated release workflow.

Preview the next version:

.. code-block:: powershell

   git cliff --bumped-version --bump patch

Prepend only unreleased changes to ``CHANGELOG.md``:

.. code-block:: powershell

   git cliff --unreleased --bump patch --prepend CHANGELOG.md

Then inspect the changelog. This command prepends the new release section, but
manual cleanup is still required:

* remove the duplicated ``# Changelog`` title and introductory text if they
  appear below the new release section
* add the new compare link near the bottom of ``CHANGELOG.md`` before the
  existing release links, for example:

  .. code-block:: markdown

     [1.4.2]: https://github.com/tifuchs/orGUI/compare/v1.4.1..v1.4.2

2.2 Add release notes
~~~~~~~~~~~~~~~~~~~~~

Edit the new release section by hand. Place the user-facing release notes
directly below the new ``## [X.Y.Z]`` heading and commit range, and above the
generated commit groups such as ``### Added`` or ``### Fixed``.

Include only high-signal notes:

* one short sentence classifying the release, such as bugfix, maintenance, or
  feature release
* important user-visible changes
* upgrade or dependency warnings
* known issues
* scientific correctness fixes that affect calculations or previous results

Example:

.. code-block:: markdown

   ## [1.4.1] (2026-03-27)

   [260ab5d](https://github.com/tifuchs/orGUI/commit/260ab5d)...[3f5ce1a](https://github.com/tifuchs/orGUI/commit/3f5ce1a)

   This is a bugfix and maintenance release.

   **Attention: pyFAI versions 2025.12.0, 2025.12.1 and 2026.2.0 cause a
   software crash when selecting a detector in the machine parameters.**

   Two critical bugs were fixed that affect physics calculations:

   - Allowed Bragg reflections close to the search bounds were sometimes not
     found.
   - The complex phase of the Fourier components rho_G of the electron density
     along z was incorrectly calculated.

   ### Added

   - Generated commit list starts here.

Regenerate the documentation release-notes page:

.. code-block:: powershell

   python doc\generate_release_notes.py

Test the Sphinx documentation build:

.. code-block:: powershell

   python -m sphinx -b html doc\source doc\build\html

3. Commit, tag, and push the release
------------------------------------

.. code-block:: powershell

   git diff -- CHANGELOG.md
   git add CHANGELOG.md
   git commit -m "docs: update changelog for vX.Y.Z"
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   git show --stat vX.Y.Z
   git push origin master
   git push origin vX.Y.Z

Use the version from the new ``CHANGELOG.md`` release section for ``vX.Y.Z``.
The tag push starts ``.github/workflows/release.yml``.

4. Verify the workflow
----------------------

Check that GitHub Actions:

* builds the wheel and source distribution
* runs ``twine check``
* extracts the release notes from ``CHANGELOG.md``
* creates the GitHub release
* publishes to PyPI

If publishing to PyPI has already succeeded, do not reuse the same version
number for a fix.

5. Verify the published package
-------------------------------

.. code-block:: powershell

   python -m venv .venv-release-check
   .\.venv-release-check\Scripts\python -m pip install --upgrade pip
   .\.venv-release-check\Scripts\python -m pip install "orGUI[full]==X.Y.Z"
   .\.venv-release-check\Scripts\python -c "import orgui; print(orgui.__version__)"

Confirm the GitHub release, PyPI page, attached artifacts, and Read the Docs
build.

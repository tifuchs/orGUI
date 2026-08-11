# AGENTS.md

## Scope

This directory contains the Sphinx documentation source and the tooling that
keeps it in sync with releases:

- `source/*.rst`: hand-maintained topic and reference pages (`geometry.rst`,
  `image_integration.rst`, `reciprocal_space_reconstruction.rst`,
  `ctr_structure_factors.rst`, `beamline_backends.rst`, ...), wired together by
  the toctree in `source/index.rst`.
- `source/api.rst`: autodoc/autosummary entry point into the `orgui` package.
- `source/release_notes.rst`: **generated**. Do not hand-edit; see below.
- `source/release_procedure.rst`: the authoritative step-by-step release
  checklist (changelog, tagging, publishing, verification).
- `source/developer.rst`: developer-facing toctree root.
- `source/conf.py`: Sphinx configuration (`furo` theme, `nbsphinx`, autodoc).
- `generate_release_notes.py`: extracts hand-written release notes from
  `../CHANGELOG.md` into `source/release_notes.rst`.
- `design/*.md`: engineering design records, **not** part of the Sphinx
  toctree and not user-facing. They record why a subsystem is built the way
  it is, what was measured, and what remains open — including decisions that
  were tried and reversed. Written for whoever picks the work up next, not
  for users; a user-facing consequence belongs in `source/*.rst` and
  `../CHANGELOG.md` as well.
- `build/`, `_build/`: build output, not source.

Use the repository root `AGENTS.md` together with this file.

## Two Documentation Tracks

1. **Topic/reference docs** (`source/*.rst`): describe features, workflows,
   geometry conventions, and APIs. Hand-maintained, evolve with the code they
   document.
2. **Changelog / release notes** (`../CHANGELOG.md` ->
   `generate_release_notes.py` -> `source/release_notes.rst`): a chronological
   record of user-facing changes, curated at release time per
   `source/release_procedure.rst`.

3. **Design records** (`design/*.md`): why a subsystem is the way it is, what
   was measured, what is still open. Not published, not in the toctree.

A given code change may need updates on either track, both, or neither.
Internal refactors typically need neither — but a refactor that invalidates a
measurement or a decision recorded under `design/` should say so there, and a
prediction that turned out wrong is worth recording rather than deleting.

## Updating Topic Pages

- Write plain reST matching the existing files. (The repository root
  `AGENTS.md` describes Markedly Structured Text as a future direction for
  *docstrings*; the `source/*.rst` files today are conventional Sphinx reST —
  match what is already there rather than introducing a new dialect.)
- Before creating a new page, look for an existing one that already covers the
  topic (e.g. reconstruction changes usually belong in
  `reciprocal_space_reconstruction.rst`, CTR model changes in
  `ctr_structure_factors.rst`). Prefer extending an existing page.
- If you do add a new page, add it to the toctree in `source/index.rst` —
  pages not listed there are unreachable from the built docs.
- State units and coordinate conventions the same way the code does. Do not
  invent a different convention for the doc page than the one in the
  docstring or `orgui/datautils/xrayutils/AGENTS.md` / `orgui/app/AGENTS.md`.
- Keep edits localized to the feature that changed. Don't reflow or rewrite
  unrelated prose on the same page.

## Updating the Changelog

`CHANGELOG.md` has an `## [Unreleased]` section at the top. Add short,
high-signal prose bullets there for user-facing changes, grouped under a
theme heading consistent with the existing style (e.g. "Scientific and
analysis additions:", "CLI and scripting highlights:", "Other user-visible
changes:") — see the existing `[1.5.0]` section for the pattern.

This is **not** the same format as the auto-generated `### Added` / `###
Fixed` / `### Removed` groups that `git cliff` inserts at release time (see
`cliff.toml` and `source/release_procedure.rst` section 2.1) — those are
generated from conventional-commit messages only when a release is cut. Do
not hand-write `###`-style group headings into the Unreleased section.

Before adding a bullet, check whether the change is already covered by an
existing Unreleased bullet — do not duplicate.

Never hand-edit `source/release_notes.rst`. After changing
`CHANGELOG.md`, regenerate it:

```powershell
python doc\generate_release_notes.py
```

## What Counts As User-Facing

Document it (changelog and/or topic page):

- GUI behavior, dialogs, CLI flags/behavior
- config file fields or file-format changes (see `examples/config_*`)
- calculation results, scientific conventions, or anything that changes
  previously produced numbers
- public API functions/classes and their documented contracts
- installation/packaging changes

Usually skip:

- pure internal refactors with no observable behavior change
- test-only changes
- CI/build tooling changes
- comment- or docstring-only clarifications that don't change a documented
  contract

`cliff.toml`'s commit parsers are a coarse heuristic for release-time grouping
(e.g. it groups anything matching `test` under "Fixed") — use judgment about
user impact when deciding whether something needs a changelog/doc entry, not
the cliff group label.

## Regenerating and Validating

```powershell
python doc\generate_release_notes.py
python -m sphinx -b html doc\source doc\build\html
```

Treat new Sphinx build warnings (broken references, toctree issues) as
signal, not noise.

## Avoid

- Do not invent user-facing prose for a change you can't verify from its
  commit message and diff. If the intent is unclear, say so rather than
  guessing.
- Do not touch scientific unit/convention statements in doc pages without
  checking the underlying code and its local `AGENTS.md`.
- Do not edit or reflow already-released `CHANGELOG.md` sections.
- Do not remove or rewrite existing hand-written notes when adding new ones.

## Automated Doc-Sync

Two automations walk commits since the last release tag and draft
`CHANGELOG.md` / topic-page updates for gaps they find. Both follow the rules
in this file and never touch `master` or hand-edit `source/release_notes.rst`
directly; this file is the durable reference they read for conventions.

- **`.github/workflows/docs-sync.yml`** (durable): runs daily on a schedule
  (`push` is not a supported trigger event for `anthropics/claude-code-action`)
  against the current tip of `master`, plus manual `workflow_dispatch`. A
  `check-open-pr` gate job runs first and skips the sync entirely if a PR
  from `docs/auto-sync` is already open — this avoids force-pushing new
  commits underneath a PR while it's still awaiting review (e.g. if nobody
  gets to it for a week). It force-updates the `docs/auto-sync` branch with
  any missing changelog/doc entries and never opens a pull request itself —
  review and merge that branch (and open the PR) manually. Requires a
  `CLAUDE_CODE_OAUTH_TOKEN` repository secret (generated via `claude
  setup-token` from a Claude Pro/Max subscription, not a metered API key).
- A local, session-scoped `CronCreate` task (set up ad hoc in a Claude Code
  session, not checked into the repo) may also run the same procedure against
  a local sibling worktree for the same `docs/auto-sync` branch, without
  pushing. It expires after 7 days or when the session ends — treat the
  GitHub Actions workflow as the source of truth if the two ever disagree.

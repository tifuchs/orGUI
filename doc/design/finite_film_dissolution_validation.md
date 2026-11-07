# Finite-film dissolution implementation and validation

Date: 2026-11-07. Branch: `codex/finite-film-dissolution`.

## Scope and provenance

The two dissolution implementation/diagnostic commits have been rebased onto
authoritative local master `5fc381b0815b3d92095acbe36d02c7c4fdb7b460`.
Its finalized prerequisites are `af86a67` (crystallinity and optical inventory),
`7980c2e` (DWBA terrace ensembles), and `5fc381b` (native reuse and batching).
The initial prerequisite snapshot `eeeb30b` is no longer in branch ancestry.
The old branch is preserved at
`codex/backup-finite-film-pre-rebase-10dc5e2`.

Master's `_CTRincoherent_dwba.py`, `cpp/CTRdwba_cpp.hpp` dispatcher and
`test/test_CTRdwba_reuse.py` are unchanged by the replay. The dissolution guard
now lives in master's `cpp/CTRdwba_records.hpp`: the homogeneous-bulk shortcut
requires both one atomic record and exactly two optical media. The small
`CTRdwba.py` additions only admit and enumerate frozen material components.
The original checkout and canonical fit files were not edited.

## Validation after rebasing onto final DWBA master

- Rebuilt the native extension against the finalized master headers, including
  the relocated homogeneous-bulk guard. All 254 dissolution, DWBA, reuse,
  motif, incoherent-DWBA, optical-profile and optical-inventory tests pass from
  source and independently from an isolated wheel installation. These include
  all 55 dissolution regressions. Both runs report 28 existing NumPy deprecation
  warnings.
- `ruff check orgui --no-cache` and `git diff --check master` pass. Regenerating
  release notes produces no additional changes.
- The zero-anchor causal diagnostic passes unchanged: final master's native
  reuse changes do not remove the optical-center discrepancy described below.
  The coordinate-only probe still reduces all three relative discrepancies
  below `1e-9`; this is not a production repair.
- The original checkout's HEAD, tracked status and binary tracked-diff SHA256
  match their pre-rebase snapshot. Master remains unchanged. No refit or new C1
  replay was performed during this rebase; the data results below are historical.

The first wheel attempt exposed a local packaging hazard: Meson's source-tree
installation also copied an old in-place `.pyd`, producing a duplicate archive
entry that shadowed the newly compiled extension. That wheel failed two native
extreme-range reuse tests and is not the validated artifact. Temporarily moving
the worktree's in-place binaries outside `orgui` during wheel construction,
then restoring them, produced one native archive entry with the same SHA256 as
the fresh build. No packaging source or master files were changed.

Validated wheel: `orgui-1.5.1.dev387+gd8a81c26c-cp314-cp314-win_amd64.whl`.
SHA256: `aeb065ac2390a02e26f623799222190fbc0ae8b7f16e4c374d6c7d6a18cf796c`.
Native extension SHA256:
`4373ef0907b7178b4659173040f926cfe635bd0d853d7fd2e550cbfa1a96062f`.
The wheel precedes the documentation-only amendment recording these results.

## Implemented contract

- `FiniteFilmDissolution` wraps the native crystal and exposes fit parameters
  `t`, `k_up`, `k_down`, `offset` and `incoherent_fraction`.
- Named finite owners include their populations in both adjacent interfaces.
  Algebraic reference subtractions retain their provenance. The semi-infinite
  bulk is protected; removable substrate requires an explicit finite Film.
- All probability below exhaustion accumulates at exhaustion. This is an
  explicit physical cap, not conditioning or numerical tail clipping.
- The frozen, material-selective closure preserves surviving populations and
  vacancies. It does not replace removed LaNiO3 with pristine STO.
- Native termination motifs require complete, exactly full top-owner material
  with matching Film chemistry and geometry. Eligibility crosses component
  labels when those physical records agree. No occupancy threshold makes a
  partially occupied population full. Incomplete motifs use the explicit
  unreconstructed-survivor policy. The `preserve` policy is also available.
- Native cells calculate complex kinematic amplitudes and DWBA. Both optical
  field policies and coherent/incoherent mixing use complete states. Existing
  legacy surfaces retain their support restriction.
- Upper-tail mass and first/second index-moment certificates, normalization
  factor, physical-height moment bounds, separate interface-support audit,
  exhausted probability and species inventories are reported. Certificates
  concern truncation/redistribution, not floating-point rounding or a universal
  DWBA error bound. `max_states` fails explicitly rather than clipping.
- `tight_support()` is a separate opt-in numerical diagnostic using unchanged
  native probability/first-moment budgets. No legacy tolerance was relaxed.

The frozen native-cell adapter preserves optical plane geometry and component
grouping. Zero-weight protected records are retained; the optical-position
correction below prevents them from anchoring nonempty merged groups, while
preserving wholly empty groups.

## Executed tests before the rebase

- 55 focused dissolution cases pass from source and the rebuilt isolated
  installed wheel. They cover independent PMFs/moments; plain-Film state
  amplitudes; interface populations; explicit finite substrate; a Film between
  two interfaces; unequal lateral areas; nonuniform spacing; inventory and
  source immutability; support convergence; full coherent/incoherent mixing;
  DWBA polarization, common phase, shared/state-specific fields; native
  optimizer updates; JSON, native text plus companion metadata, deepcopy and
  pickle; live changes and rejected-candidate recovery.
- The complete scientific-core run passed 2,053 tests and 247 subtests, with
  seven skips. Its sole failure was a missing native ROI extension in the new
  source worktree. After copying the already-built extension locally, that
  reconstruction test passed. Later dissolution-only changes were checked by
  the focused suite; legacy implementation was unchanged by those changes.
- `ruff check orgui` and lint on the example/validation scripts pass.
- The synthetic example runs for both ordinary Film removal and an explicit
  substrate Film. Changelog-generated release notes were refreshed. Sphinx
  HTML builds succeed; unrelated pre-existing warnings remain.

The pre-rebase wheel was built locally without changing the normal installed orGUI:
`orgui-1.5.1.dev383+geeeb30b38.d20261107-cp314-cp314-win_amd64.whl`.
SHA256: `8d5914ff7c4fc740d4e941931bdc32fb2d64ca37ba378d82c8c9cf7c84ffcdbd`.

## C1 fixed-parameter replay

`benchmarks/validate_dissolution_c1.py` requires the companion metadata and
checks these input hashes:

- `0056_C1.xtal`: `0e8c5f28b843d0a7877afa89e50322ed9ca40fe9ea63bd1cb3fc8956bc639b17`
- `0056_C1.json`: `63e5cb0951fa17401e6d7e183411183179ade2bdc963adc07e2cc9ddb8713adc`

The replay uses all 5,596 non-specular observations, 26,042 original resolution
nodes, all four recorded domain rotations, original intensity convolution and
the nine saved rod scales. No structural parameter or nuisance scale is fitted.

The legacy replay reproduces saved amplitudes to maximum absolute error
`1.1368683772161603e-13`, with chi-squared/N `8.394242571229753`.
The new model gives `8.394242357554576`; its maximum prediction change is
`1.5610775954133516e-6` of a supplied standard error. Tightening parent numerical
budgets to `1e-12` changes no predicted point. Source and isolated installed
results agree. Current realized mean/RMS are 64.1986672656/5.15464027028 Angstrom.
Full exhaustion has probability `9.6845727817e-28` for the original interface
representation (exhausted index -49).

Increasing interface support from offsets [-12, 11] to [-24, 23], with
`fixed_ucs=24`, changes the full resolution-processed C1 prediction by at most
`4.9586388329e-9` supplied standard errors. Chi-squared/N becomes
`8.394242357436742`. The mean changes by approximately `4.1e-13` Angstrom,
RMS by `2e-15` Angstrom, and remaining species inventories by at most
`1.1e-11` atoms per reference lateral cell. The more deeply represented
exhausted index is -61, with probability `1.1495901666e-38`; both exhaustion
probabilities are negligible for C1. The separate interface missing-inventory
bound decreases from approximately `5.05e-5` to `2.67e-12` plane equivalents.
This passes the declared `1e-3`-standard-error kinematic convergence criterion.
It does not establish convergence for a future broader fit or for DWBA.

The old -25 boundary and complete exhaustion are different events. Representative
unbounded parents at fixed C1 mean/asymmetry give:

| Parent RMS (Angstrom) | P(H < -25) | E[(-25-H)+] (planes) | Full-exhaustion probability |
| --- | --- | --- | --- |
| 5.15464 | 2.35564e-10 | 3.06505e-10 | 9.68457e-28 |
| 5.98506 | 9.22799e-9 | 1.29976e-8 | 1.15841e-23 |
| 7.55706 | 9.63678e-7 | 1.58113e-6 | 2.83838e-18 |
| 9.43618 | 3.11942e-5 | 6.16222e-5 | 5.36490e-14 |
| 13.3320 | 1.41502e-3 | 4.01372e-3 | 5.51888e-9 |
| 18.8453 | 1.47114e-2 | 6.39211e-2 | 1.01331e-5 |

Plane-equivalent excess is geometric, not a chemical atom count. The wrapper's
species inventories provide the chemical accounting. These broader values
are audits, not fitted results or certificates for an unconstrained roughness.

## Remaining physical and numerical limits

The interface retains finite native numerical support. Report its uncertainty
separately; broad distributions require renewed checks of interface support,
inventories, realized front moments and scattering. In a porous interface the
structural front is not a unique atomic surface height. Full exhaustion means
all represented selected-owner material has been removed.

The synthetic DWBA interface-support fixture now passes the proposed `1e-9`
relative-observable gate after the optical-position correction below. Before
that correction its native DWBA error was `2.0573e-4`, inherited also by the
shared and state-specific dissolution ensembles. The correction does not
establish invariance under arbitrary repartition: existing grouping membership
still depends on input origins, and genuine nonzero tails affect optics.
Optical-profile convergence remains a prerequisite for a DWBA refit or a claim
of DWBA-invariant repartition. No numerical tolerance was relaxed or native
scattering kernel replaced.

Water, other overlayers, overlapping interface supports, non-unit component
weights, different outer domains across components and depth-dependent
surface/interface crystallinity profiles are explicitly unsupported. Constant
native crystallinity is supported. No automatic migration enlarges the Film,
moves the substrate boundary, restores STO vacancies or changes saved fits.

Stage 8 (a new analysis fit) was not performed. No claim of an unconstrained
C1 roughness or production DWBA fit readiness follows from this implementation.

## Root cause identified: zero-weight records select optical positions

The support discrepancy above is now isolated. Reproduce the causal experiment
with:

```powershell
python -m benchmarks.audit_dissolution_dwba_partition --diagnose-zero-anchors
```

`EpitaxyInterface.createInterfaceCells` retains algebraic lower-material
subtraction domains with weight `-0.0` when the native underlying-component
replacement is used. `_optical_profile` emits these as zero-density rows at
their unstrained positions. Extending interface support extends these empty
records too. At one representative plane in the expanded fixture:

| Record | Position (Angstrom) | Weight | Optical contribution |
| --- | --- | --- | --- |
| Upper material | 30.18 | 1 | Nonzero |
| Lower material | 30.18 | 0 | Zero |
| Unused sharp-boundary subtraction | 30.00 | -0.0 | Zero |

`CTRdwba._combine_tagged` sorts records by position, coalesces origins within
0.6 Angstrom, sums their density/reference shares, and selects the *lowest*
position without checking whether that row contributes any material. Therefore
the empty row moves the merged optical center from 30.18 to 30.00 Angstrom.
The public `CTRoptics.combine_profiles` has the same position-selection rule.
`profile_boundaries` then places interfaces at midpoints of these displaced
centers. The reference fields and DWBA amplitudes consequently change despite
unchanged atomic positions and occupancies.

In this fixture twelve centers move by 0.18 Angstrom. The raw delta/beta arrays
are exactly unchanged (`max_raw_delta_beta_change == 0`). The discrepancy is
already present before optical slab simplification and native kernel evaluation.

The diagnostic temporarily changes only the representative coordinate of each
nonempty coalesced group: it selects the lowest position with a nonzero tagged
share. It retains all shares, signed contributions, empty groups and the 0.6
Angstrom grouping tolerance. It changes no source atoms, distributions, native
kernel, fit parameters or physical boundary policy.

| Calculation | Original maximum relative support change | Coordinate-only probe |
| --- | --- | --- |
| Ordinary native DWBA | 2.05732655e-4 | 7.10543e-14 |
| Dissolution, shared fields | 1.48178726e-4 | 6.67244e-14 |
| Dissolution, state-specific fields | 1.24291133e-4 | 9.43690e-14 |

Two controls support the causal attribution: setting both optical
simplification tolerances to zero leaves the original discrepancies unchanged;
setting interface displacement to zero makes occupied and empty positions
coincide and reduces all three discrepancies below 4e-15 without the probe.

That experiment identified the cause before the runtime repair. Simply deleting
every zero-density row would not be justified: signed reference contributions
and wholly empty physical regions must survive. The production correction is
recorded below; the diagnostic can still reproduce the legacy coordinate rule
using a temporary patch.

## Optical-position correction

Both `CTRdwba._combine_tagged` and `CTRoptics.combine_profiles` now retain the
lowest contributing origin within each existing coalesced group. Tagged
assembly examines every individual signed share, including areal inventory
slots beyond `density_count`, before cancellation. Public profile assembly
examines each input delta/beta row. The test is exact nonzero, with no density
threshold. Wholly empty groups keep their original lowest coordinate. Group
membership, the 0.6 Angstrom tolerance, sums, populations, boundaries between
materials and all physical dissolution policies are otherwise unchanged;
optical midpoint boundaries are recomputed from the corrected centers.

The independent pre-fix tests reproduce coordinate and observable failures.
Regression coverage checks signed cancellation, tiny nonzero tails, separate
reference-only inventory, empty media, public/tagged profile coordinates and
midpoint boundaries, native complex amplitudes and intensities for s and p,
and shared/state-specific dissolution at incoherent fractions 0, 0.2 and 1.
Repartitioned stacking coordinates are compared to a four-machine-epsilon
roundoff budget scaled by the coordinate range, in Angstrom; observables use
the original `1e-9` gate.

The updated diagnostic verifies the production fix and reproduces the old
coordinate rule only in its temporary legacy control. Maximum relative
support errors are `7.10543e-14` for native DWBA, `6.67244e-14` for shared
dissolution fields, and `9.43690e-14` for state-specific fields. The legacy
control reproduces the original table above. This is a correction to numerical
optical assembly, independent of adopting a dissolution model. It changes no
kinematic atoms, fit data, distribution weights or probability/moment budgets.

Validation: the complete `test_CTR*.py` run passed 1,303 tests and 193 subtests,
with three skips. Its single failure was an existing self-snapshot expecting
an empty subtraction origin at 0 Angstrom. Inspection of its component rows
confirmed that the lowest contributor is at 0.007419160637154795 Angstrom; only
that coordinate expectation was updated, with all density and kinematic
expectations retained. The corrected case and 21 new regressions then passed
together. One additional same-origin signed-cancellation regression also
passed, bringing coverage to 22 new cases and 1,305 distinct passing CTR cases
across the broad run and focused reruns. `ruff check orgui` (plus the audit
script) and `git diff --check` pass. C++ sources are unchanged by this fix;
tests use the freshly rebuilt final-master extension identified above. No
data were refitted.

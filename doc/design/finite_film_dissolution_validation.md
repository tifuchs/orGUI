# Finite-film dissolution implementation and validation

Date: 2026-11-07. Branch: `codex/finite-film-dissolution`.

## Scope and provenance

The worktree starts at repository commit
`a5a6afe16b072ab4092c6a16aec54c886b9f805e`, followed by prerequisite snapshot
`eeeb30b`. That snapshot preserves the initial uncommitted scientific-core
changes, including Skellam, multilayer, crystallinity and native DWBA work.
Implementation is a separate commit after the snapshot. The original checkout
and canonical fit files were not edited by this implementation.

The original checkout subsequently changed in `CTRdwba.py`,
`_CTRincoherent_dwba.py`, `cpp/CTRdwba_cpp.hpp` and `test/test_CTRdwba_reuse.py`.
Those later native-reuse changes remain outside this frozen worktree baseline.
Integration must retain the dissolution fix that disables the homogeneous-bulk
shortcut when the reference contains more than two optical media, even if the
kernel has been moved to another header by that later work (currently
`cpp/CTRdwba_records.hpp`, whose `homogeneous_bulk_only` expression still needs
the `media == 2` condition).

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
grouping. In particular, zero-weight protected records can anchor native
optical slab boundaries and must not be dropped casually.

## Executed tests

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

The final wheel was built locally without changing the normal installed orGUI:
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

**DWBA interface-support invariance has not passed the proposed 1e-9 gate.**
The reproducible `benchmarks/audit_dissolution_dwba_partition.py` diagnostic
shows that widening interface support can change native optical profiles even
for ordinary unwrapped crystals whose kinematic amplitudes agree. In its
synthetic motif fixture the native DWBA maximum relative change is 2.0573e-4;
the maximum error divided by the maximum reference intensity is 3.1052e-7.
The new shared and state-specific ensembles also inherit optical grouping
sensitivity. This is separate from exact-state/native-DWBA agreement, which
passes. The tolerance was not relaxed, and no replacement optical solver was
introduced. Optical-profile convergence remains a prerequisite for a DWBA
refit or a claim of DWBA-invariant repartition.

Water, other overlayers, overlapping interface supports, non-unit component
weights, different outer domains across components and depth-dependent
surface/interface crystallinity profiles are explicitly unsupported. Constant
native crystallinity is supported. No automatic migration enlarges the Film,
moves the substrate boundary, restores STO vacancies or changes saved fits.

Stage 8 (a new analysis fit) was not performed. No claim of an unconstrained
C1 roughness or production DWBA fit readiness follows from this implementation.

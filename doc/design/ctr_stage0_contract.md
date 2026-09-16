# CTR targeted implementation: Stage 0 contract

Status: implementation note, 2026-09-16. Numerical behavior was captured at
commit `d52ab8f5234e`, before Stage 1 of
[`ctr_targeted_implementation_plan.md`](ctr_targeted_implementation_plan.md).

This note records what the characterization fixture means and which metadata
is not available. It does not claim that the frozen numbers independently
validate the physics.

## Frozen extraction behavior

`orgui/app/test/fixtures/ctr_stage0_contract.json` and
`test_ctr_stage0_contract.py` pin two compact cases:

- A stationary extraction with detector counts and one-sigma count errors,
  radian diffraction angles, an exposure-times-monitor divisor, and the
  current stationary Lorentz-to-`F2_hkl` conversion.
- A mechanical rocking extraction with a degree-valued motor axis, two
  background regions, per-frame normalization inside the angular integral,
  radian detector acceptance, and the current relative `F2_hkl` scale.

The fixture also pins a pre-typed correction-settings dictionary. Its sample
dimensions are meters and `beam_flux_density` remains photons per second per
square meter. Missing fields remain missing; zero is not substituted for an
unknown calibration. Existing tests in `test_config_data.py` and
`test_peak1Dintegr.py` cover typed NeXus round trips, legacy JSON database
dispatch, and both stored correction-group layouts used by the reducer.

Strict expected-failure tests record the two Stage 1 defects without changing
application code: the arm-following polarization evaluator ignores a nonzero
polarization axis, and rocking acceptance ignores stored actual-arm angles.

## Settings-loader audit

All `examples/config_*` files use the legacy INI loader. They define detector,
lattice and diffractometer values; none selects a normalization monitor or
stores integration-correction provenance. `ConfigData.from_ini` therefore
starts with an empty `CorrectionState`. Current database snapshots write the
typed `integration_corrections` layout (schema version 2), while
`ConfigData.from_nxdict` still recognizes the older group containing one
opaque JSON dataset. Neither old layout distinguishes rate-like from
frame-integrated monitors.

## Backend counter-semantics audit

The implementation assumption is that a selected incident-beam monitor is
**rate-like**. It is therefore multiplied by exposure time in the fluence
divisor. This matches the current normalization arithmetic. The present
backend contract names arrays but does not persist a counter kind or unit, so
the assumption must remain visible in provenance until a beamline-specific
calibration establishes the counter-to-total-flux conversion.

| Backend | Persisted candidates | Meaning established by code | Rate/integrated status |
|---|---|---|---|
| Base `Scan`, imported images | none | no monitor metadata | unavailable |
| Built-in P212 | implicit `exposure_time` when present | frame exposure, seconds by loader convention | not a flux monitor; P212 flux counters are unavailable |
| ID31 fast/BLISS variants | `exposure_time`, `time`/`elapsed_time`/`epoch`, `srcur`, `mondio`, plus electrochemistry values | exposure is seconds; time/epoch are coordinates; `srcur` and `mondio` are often divided by their scan mean and are therefore dimensionless relative factors | selected beam monitors are treated as rate-like relative-flux proxies |
| Example CHESS QM2 | `diode`, `ic1`, `ic2`, `emon`, `nemon`, `pemon` | raw SPEC columns; `ic2` is ion chamber 2 and is the expected primary incident-beam monitor | `ic2` is treated as rate-like; a later calibrated factor converts its reading to total photons/s |
| Example P212 | `beckvolt_1`, `beckvolt_2`, `eh3_entrance`, `oh2_diode1`, `oh2_diode2`, `petracurrent`, `timestamp` | raw scan columns; timestamp is not a flux monitor | a selected beam monitor is treated as rate-like; its calibration is unresolved |
| Segmented/interlaced scans | intersection of segment counters, plus implicit exposure time | values are concatenated/reordered; semantics are inherited | selected monitors remain rate-like; units and calibration are not merged |

`current`, `potential`, `scaled_potv2f`, motor positions and timestamps are not
incident-flux monitors merely because they appear in `auxillary_counters`.
Until Stage 3 adds an explicit monitor unit and calibration relation, the
legacy behavior is exactly exposure time times the selected rate-like monitor
product recorded by `monitor_corrections`; it must be shown as legacy rather
than silently assigned an absolute photon scale. For QM2, the intended primary
monitor is `ic2`, with a future scalar calibration relating its rate to total
incident flux in photons/s.

## Stage 1 inputs now pinned

- Polarization: at the calibrated arm position,
  `polarizationAtPoints` must match `polarizationArray` for nonzero axes and
  mixed polarization fractions. The axis-zero baseline already passes.
- Acceptance: scan arm values are degrees at the storage boundary and must be
  converted to radians before calling `out_of_plane_acceptance`. The fallback
  for old scans with no arm values must remain explicit.

The rolled-detector case uses the published review parameters: 172 micrometer
pixels, 0.15 meter distance, 30 degree detector roll, ROI center `(row,
column) = (500, 440)`, height 60 pixels, incidence 10 degrees and actual
`(gamma, delta)` arm angles `(30, 40)` degrees. It is a counterexample to arm
invariance, not a statement about a typical experimental error.

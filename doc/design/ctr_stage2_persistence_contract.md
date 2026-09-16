# CTR Stage 2 persistence contract

Status: implemented, 2026-09-16.

Stage 2 introduces persistence and dispatch boundaries only. It deliberately
does not enable the total-flux calculation or move normalization/footprint
between extraction and reduction. Existing numerical defaults and legacy
dataset meanings remain unchanged.

## Requested settings

`configuration/orgui/integration_corrections` now has schema version 3. The
legacy `beam_flux_density` field remains photons/(s m2). The distinct
`total_incident_flux` field is photons/s and is never inferred from that
density. Primary-monitor name, rate/integrated kind, unit, calibration reading
and calibration exposure have separate fields. Horizontal interception is an
explicit mode plus optional fraction. Missing values mean unknown or not
configured; they do not mean zero flux, full interception or a unity divisor.

Measured vertical profiles store their normalized positions [m] and density
[1/m] in addition to the source path and import settings. Analytical profiles
remain reproducible from their named shape and parameters.

## Applied curve record

New extractions retain all established groups and add a sibling named
`ctr_curve_v3`. It is not named `rois` or `counters`, so an older reducer cannot
mistake this branch for the legacy unnormalized input. The record separates:

- identity: algorithm, output quantity and scale convention;
- normalization: explicit status, exact divisor and contributing counters;
- illumination: explicit status, exact divisor, named legacy/new convention,
  overlap components and alpha in radians when available;
- base: center/background sums and variances plus the reversible background-
  subtracted scalar curve and variance before normalization/illumination;
- pixel corrections: the established combined center/background factors,
  with separate polarization-only fields reserved for Stage 4;
- geometry: actual detector arms in radians and ROI location/size in pixels;
- profile: embedded profile provenance rather than only an external path.

Correction statuses are exactly `applied`, `not_applied`, `unavailable` or
`unknown`. An absent dataset is never read as unity.

Stationary records capture the scalar curve immediately before the existing
normalization and active-area divisors. Rocking extraction currently applies
neither divisor, so both statuses are explicitly `not_applied`; its current
legacy reducer remains responsible for them until Stage 5.

## Compatibility and dispatch

The existing `rois`, `croibg`, `Cfactors_croi` and `Cfactors_bgroi` datasets
retain their meanings. Current Stage 2 records identify their algorithm as
`legacy_stationary_roi_v1` or `legacy_rocking_roi_v1`, allowing the existing
paths to reproduce the same numbers.

The rocking loader checks the sibling record before returning the legacy ROI
curve. A record with a non-legacy algorithm (for example the later framewise
total-flux convention) is refused by the legacy reducer instead of being
silently reinterpreted. Incomplete older records remain accessible through the
legacy diagnostic path with unknown provenance. Reduction outputs continue to
be added as new measurement groups, so re-reduction does not overwrite the
source extraction record.

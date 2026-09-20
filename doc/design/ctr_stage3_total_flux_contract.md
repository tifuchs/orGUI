# CTR Stage 3 total-flux numerical contract

Status: implemented, 2026-09-16.

Stage 3 adds pure numerical building blocks for calibrated total incident flux.
It does not change either integration path, the current GUI/config behavior or
the legacy flux-density APIs. Stage 5 will select and persist this convention
during extraction.

## Per-frame incident photons

`corrections.normalization.frame_fluence` returns `Q_f` in photons. A constant
total flux `Phi` [photons/s] gives `Q_f = Phi T_f`. A rate-like primary monitor
uses

```text
Q_f = Phi_ref T_f M_f / M_ref.
```

An integrated monitor already contains the exposure and instead uses

```text
Q_f = Phi_ref T_ref U_f / U_ref.
```

There is no second multiplication by `T_f` in the integrated case. Monitor kind
is explicit; the numerical core does not infer it from a counter name or unit.
The beamline counters expected for the first wiring pass, including QM2 `ic2`,
are rate-like after applying their configured conversion to total incident
flux. `relative_frame_fluence` applies the same exposure rule without claiming
photon units when no absolute calibration exists.

All fluxes, exposures, monitor readings and calibration references must be
finite and strictly positive. Calibration references correspond to the stated
reference flux; no scan-mean normalization is introduced.

## Beam/sample illumination

For a normalized vertical beam profile, sample length `L` [m], incidence angle
`alpha` [rad] and explicit horizontal intercepted fraction `f_x`,

```text
f_z   = profile.flux_on_sample(alpha, L)
f_hit = f_z f_x
H     = f_hit / sin(alpha).
```

`corrections.activearea.intercepted_fraction` returns `f_hit`, while
`illumination_divisor` returns the dimensionless `H`. These are different
quantities: photons physically hitting the sample are `Q_f f_hit`, whereas the
retained CTR reduction divides counts by `Q_f H`.

`BeamProfile.flux_over_sine` evaluates the ratio directly. Built-in measured,
distribution and analytical Gaussian profiles continue it to the finite
grazing-incidence limit, avoiding numerical zero divided by zero. Existing
external profile subclasses remain instantiable through a non-abstract default
implementation.

The vertical profile supplies no horizontal geometry. `f_x` is therefore an
explicit number in `(0, 1]`; one means the full horizontal beam is intercepted.
No general two-dimensional clipping model is implied.

## Structure-factor scale

`corrections.measurement.total_flux_prefactor` returns

```text
K = r_e^2 lambda^2 / A_u^2,
```

with wavelength in Angstrom and surface unit-cell area in square Angstrom.
`structure_factor_squared_from_photon_yield` reduces
`Y = N_net / (Q_f H)` with the existing scan-mode angular factor and optional
detector efficiency. `photon_yield_from_structure_factor` is its exact forward
inverse.

The total-flux identity is pinned independently for a separable 2-D Gaussian:
total flux times `H` equals peak flux density times the existing effective beam
area. Thus `Q H` replaces the legacy density-times-area fluence; it is not an
extra footprint correction. `scale_factor`, `structure_factor_squared` and
`integrated_intensity` retain their legacy density semantics unchanged.

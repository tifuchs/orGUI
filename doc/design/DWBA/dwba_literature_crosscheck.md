# Literature cross-check for the bulk DWBA specular note

Companion to `dwba_bulk_specular_math.tex`. This note records what external
literature does and does not corroborate in that derivation, so the provenance
of each step is traceable without repeating the search.

## Headline

No published DWBA treatment of the **specular** CTR carrying a mean-density
subtraction was found. That is a gap in coverage, not a contradiction:

- The DWBA-CTR literature (Vineyard, Kaganer, Stepanov, Holy) computes
  **non-specular** rods, where Eq. `half-space-transform` makes the
  mean-density term vanish identically. Nobody needed it.
- The specular reflectivity literature (Parratt, Nevot-Croce, Sinha, Fenter)
  works with **continuum** density profiles, where the crystalline lattice sum
  `1/(1 - exp(-i Qz c))` never appears, so the double count cannot arise.

The note sits in the gap between the two. Every individual ingredient is
independently confirmed; the combination is not. Validation therefore has to be
numerical, see "Recommended benchmarks" below.

## 1. Closest match: Kaganer, Phys. Rev. B 75, 245425 (2007)

*Crystal truncation rods in kinematical and dynamical x-ray diffraction
theories*, arXiv:cond-mat/0702679. Effectively the same construction as the
design note, benchmarked against multi-beam dynamical theory.

| Design note | Kaganer (2007) |
| --- | --- |
| Reference = homogeneous half-space at cell-average density (Eqs. `reference-index`, `mean-density-from-cell`) | "In the zeroth order, the scattering problem is solved for a uniform medium having the same polarizability chi_0 as the crystal under investigation" |
| Reciprocity far-field Green function (Eq. `green-far-field`) | "The DWBA formulation in terms of the reciprocity theorem in electrodynamics is most straightforward" (cites Dmitrienko & Kaganer 1987; Kaganer, Stepanov & Koehler, PRB 52, 16369 (1995)) |
| Four Renaud channels with complex per-channel `Qz` (Eqs. `channel-Qz`, `four-channel-atomic-amplitude`) | His Eq. (4): `E_DWBA = sum_{i,j=1,2} D_i^in D_j^out E_kin(hk L_ij)`, complex `L_ij` from in-medium `kz` replacing real `L` |
| `delta r = (2 pi i r_e / kappa_0) F / A_ref` (Eq. `reflection-correction`) | His Eq. (2): `E_kin = r_e lambda / (a^2 gamma_out) * i F_hkL / (exp(2 pi i L) - 1)` |

### The prefactor identity

With `lambda = 2 pi / k0`, `gamma_out = sin(alpha_f)`, `kappa_f = k0 gamma_out`
and `a^2 = A_ref`:

```
r_e lambda / (a^2 gamma_out)  =  2 pi r_e / (kappa_f A_ref)
```

so

```
E_kin = (2 pi i r_e / kappa_f) * (1 / A_ref) * F_hkL / (exp(i Qz c) - 1)
```

which is Eq. `reflection-correction` with the atomic lattice sum inserted. The
agreement is exact, including the factor `i` and the absence of any residual
`k0^2 / (4 pi)`. This independently confirms invariant 1 of the design note
against a result that has itself been checked against dynamical theory.

Two convention differences, both harmless but both worth respecting:

- **Origin.** `1/(exp(i Qz c) - 1)` versus `1/(1 - exp(-i Qz c))` differ by one
  factor `exp(-i Qz c)`, that is, a one-cell shift of the unit-cell origin.
  Because `delta r` is added coherently to `r0`, that phase is not free.
- **Exit angle.** Kaganer's normalizing wavevector is the *exit* normal
  component. On the specular rod it equals `kappa_0`, so the design note is
  correct as written, but the same normalization must not be reused
  off-specular.

### Why Kaganer never needs the subtraction

All his rods are non-specular (11L, 1-3L, 22L). His zeroth-order reference
carries `chi_0`, so on a specular rod his `E_kin` would double count exactly the
term the design note removes. His silence is consistent with, and is external
evidence for, invariant 3 (off-specular invariance).

### Caveat he raises that the design note originally did not

The mathematical surface plane must not bisect an atom. The kinematical and
dynamical calculations otherwise agree only down to about `1e-7` of the peak
intensity. For GaAs(001), atoms must be placed at fractional depths
`1/8, 3/8, 5/8, 7/8` rather than `0, 1/4, 1/2, 3/4`. This bears directly on the
design note, because `rho_bar_f * Theta(-z)` is a sharp step at `z = 0` while
`rho_f(r)` is a sum of smooth atomic clouds: the cancellation in
Eq. `contrast-bulk-cellwise` is exact only if the top cell's electron density
lies entirely in `z < 0`. Now recorded as invariant 8.

## 2. Half-space transform and the Fresnel limit

The `1/(i Qz)` regularization of the semi-infinite integral is standard in the
reflectivity literature, reached by integration by parts rather than by a
convergence factor. Vaknin, *X-ray diffraction techniques for liquid surfaces
and interfaces* (arXiv:cond-mat/0101142), Eq. (41):

```
integral rho(z) exp(i Qz z) dz = -(1 / (i Qz)) integral (d rho / dz) exp(i Qz z) dz
```

For a sharp step of density `rho` this gives amplitude `rho / (i Qz)` and hence
`R_BA = 16 pi^2 r_e^2 rho^2 / Qz^4 = (Qc / 2 Qz)^4`, that is,
`|r_BA| = 4 pi r_e rho / Qz^2`.

The design note's phase check, `delta r = pi r_e rho_f / kappa_0^2` with
`Qz = 2 kappa_0`, is `4 pi r_e rho_f / Qz^2`. Exact magnitude agreement. The
sign is also consistent with the exact Fresnel expansion
`r = (kz0 - kz1) / (kz0 + kz1) ~ + delta k0^2 / (2 kappa_0^2)`, so invariant 5
holds as stated.

## 3. Contrast relative to the reference is the orthodox DWBA rule

- **Vineyard**, PRB 26, 4146 (1982). Origin of the surface DWBA: the reference
  is a homogeneous dielectric slab carrying the crystal's *average*
  polarizability, and the resulting distorted wave then illuminates the actual
  crystal. The design note's split is Vineyard's split.
- **Sinha, Sirota, Garoff & Stanley**, PRB 38, 2297 (1988). Reference is the
  ideal sharp interface; perturbation is the difference from it. The specular
  channel does receive a first-order DWBA correction, yielding Nevot-Croce.
  Precedent for Eq. `total-reflection-amplitude`, `r_total = r0 + delta r`.
- **Holy & Baumbach**, PRB 49, 10668 (1994); Pietsch, Holy & Baumbach,
  *High-Resolution X-Ray Scattering* (Springer 2004). Perturbation
  `delta chi = chi(r) - <chi(z)>`: specular in the reference, diffuse in the
  perturbation.
- **Zhou & Chen**, Phys. Rep. 257, 223 (1995). One-dimensional specular DWBA
  with the perturbation taken relative to a reference profile; the strict 1-D
  analogue of Eq. `contrast-dwba-definition`.

## 4. Strongest conceptual corroboration

Standard dynamical theory expands `chi(r) = sum_H chi_H exp(i H.r)`. The `H = 0`
term produces refraction and Fresnel specular reflection; `H != 0` produces
Bragg and CTR scattering.

The design note's `A_ref * rho_bar_f / (i Qz)` *is* the `H = 0` term: it is
precisely the `Qz -> 0` divergence of the kinematic lattice sum,

```
F_uc / (1 - exp(-i Qz c))  ->  A_ref c rho_bar_f / (i Qz c)
```

So "subtract the mean before adding `r0`" is the DWBA restatement of standard
dynamical-theory bookkeeping. This framing is harder to argue with than the
half-space-transform argument alone and is worth leading with in discussion.

## 5. Recommended benchmarks

Since no paper does the specular DWBA-CTR directly, validation is numerical.

- **Holy & Fewster**, J. Appl. Cryst. 41, 18-26 (2008), *Crystal truncation rod
  X-ray scattering: exact dynamical calculation*. Slab decomposition with a
  Fresnel reflection/transmission matrix formalism, "correct both in the
  reciprocal-lattice points and between them", any polarization, non-coplanar.
  A slab-dynamical calculation of a bulk-terminated crystal along 00L contains
  Fresnel *and* CTR with no double counting, so it is the natural benchmark for
  Eq. `contrast-bulk-separated`.
- **Caticha**, PRB 47, 76 (1993) and PRB 49, 33 (1994). Darwin dynamical theory
  over the whole angular range 0 to 90 degrees, explicitly covering the far
  tails between Bragg peaks. The reference for whether the specular rod joins
  continuously onto Fresnel.
- **Colella**, PRB 43, 13827 (1991); **Takahashi & Nakatani**, Surf. Sci. (1995)
  and PRB 62, 3630 (2000). Darwin-approach dynamical CTR theory whose exact
  expressions reduce to the conventional ones for CTR scattering, specular
  reflectivity, and Bragg reflection *separately*. That reduction is the
  published statement nearest to the design note's decomposition.

## Source list

| Work | Link |
| --- | --- |
| Kaganer, PRB 75, 245425 (2007) | <https://arxiv.org/abs/cond-mat/0702679> |
| Vineyard, PRB 26, 4146 (1982) | <https://doi.org/10.1103/PhysRevB.26.4146> |
| Sinha et al., PRB 38, 2297 (1988) | <https://doi.org/10.1103/PhysRevB.38.2297> |
| Holy & Baumbach, PRB 49, 10668 (1994) | <https://doi.org/10.1103/PhysRevB.49.10668> |
| Kaganer, Stepanov & Koehler, PRB 52, 16369 (1995) | <https://doi.org/10.1103/PhysRevB.52.16369> |
| Zhou & Chen, Phys. Rep. 257, 223 (1995) | <https://doi.org/10.1016/0370-1573(94)00110-O> |
| Holy & Fewster, J. Appl. Cryst. 41, 18 (2008) | <https://doi.org/10.1107/S0021889807049886> |
| Caticha, PRB 47, 76 (1993) | <https://doi.org/10.1103/PhysRevB.47.76> |
| Colella, PRB 43, 13827 (1991) | <https://doi.org/10.1103/PhysRevB.43.13827> |
| Takahashi & Nakatani, PRB 62, 3630 (2000) | <https://doi.org/10.1103/PhysRevB.62.3630> |
| Renaud, Lazzari & Leroy, Surf. Sci. Rep. 64, 255 (2009) | <https://doi.org/10.1016/j.surfrep.2009.07.002> |
| Vaknin, liquid surfaces review | <https://arxiv.org/abs/cond-mat/0101142> |

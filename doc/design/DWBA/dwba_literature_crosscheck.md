# Literature cross-check for the bulk DWBA specular note

Companion to `dwba_bulk_specular_math.tex`. This note records what external
literature does and does not corroborate in that derivation, so the provenance
of each step is traceable without repeating the search.

Equation and invariant numbers below refer to the current build of
`dwba_bulk_specular_math.tex` (Ch. 1 derives the master formula; Ch. 2 the
four channels; Ch. 3 the specular coefficient; Ch. 4 the model mapping,
observables and invariants).

## Headline

No published DWBA treatment of the **specular** CTR carrying a mean-density
subtraction was found. That is a gap in coverage, not a contradiction:

- The DWBA-CTR literature (Vineyard, Kaganer, Stepanov, Holy) computes
  **non-specular** rods, where Eq. (1.27) makes the mean-density term vanish
  identically. Nobody needed it.
- The specular reflectivity literature (Parratt, Nevot-Croce, Sinha, Fenter)
  works with **continuum** density profiles, where the crystalline lattice sum
  `1/(1 - exp(-i Qz c))` never appears, so the double count cannot arise.

The note sits in the gap between the two. Every individual ingredient is
independently confirmed; the combination is not. Validation therefore has to be
numerical, see "Recommended benchmarks" below.

## 0. Sign conventions: read this first

Two conventions must be held fixed when comparing any formula below with the
design note, because they compensate each other and an error in one is
invisible if the other is also wrong.

1. **Negative-root `kz`.** Eq. (1.5) fixes `Re kz <= 0`. A wave travelling
   towards `+z` (upward, out of the sample) therefore carries `exp(+i kz z)`,
   the opposite pairing from a positive-root convention. On the specular rod
   `Qz = -2 kz > 0` by Eq. (2.4).
2. **Crystallographic time convention.** Eq. (1.1) uses `exp[i(wt - k.r)]`, so
   outgoing waves go as `exp(-i k0 r)/r` and the Weyl angular-spectrum
   prefactor is `-i/2pi`, i.e. the complex conjugate of the familiar
   `exp(+i k0 R)` form.

Consequently Eq. (1.33) reads `r_h = -2 pi i re F_h / (A_ref kz,0,f)` with a
**minus** sign, while sources that use a positive root quote the same physical
quantity with a **plus** sign. Both are correct in their own convention. The
design note never needs a conversion because Secs. 1.4-1.6 are carried out
natively in convention 1; the Born-Fresnel check, Eq. (1.34), pins the sign
independently of any external source.

## 1. Closest match: Kaganer, Phys. Rev. B 75, 245425 (2007)

*Crystal truncation rods in kinematical and dynamical x-ray diffraction
theories*, arXiv:cond-mat/0702679. Effectively the same construction as the
design note, benchmarked against multi-beam dynamical theory.

| Design note | Kaganer (2007) |
| --- | --- |
| Reference = homogeneous half-space at cell-average density (Eqs. (1.3), (1.8)) | "In the zeroth order, the scattering problem is solved for a uniform medium having the same polarizability chi_0 as the crystal under investigation" |
| Reciprocity via the reciprocal exit field (Eq. (1.22)) | "The DWBA formulation in terms of the reciprocity theorem in electrodynamics is most straightforward" (cites Dmitrienko & Kaganer 1987; Kaganer, Stepanov & Koehler, PRB 52, 16369 (1995)) |
| Four Renaud channels with complex per-channel `Qz` (Eqs. (2.1), (2.4), (2.11)) | His Eq. (4): `E_DWBA = sum_{i,j=1,2} D_i^in D_j^out E_kin(hk L_ij)`, complex `L_ij` from in-medium `kz` replacing real `L` |
| `r_h = -2 pi i re F_h / (A_ref kz,0,f)` (Eq. (1.33)) | His Eq. (2): `E_kin = re lambda / (a^2 gamma_out) * i F_hkL / (exp(2 pi i L) - 1)` |

### The prefactor identity

With `lambda = 2 pi / k0`, `gamma_out = sin(alpha_f)`, `kappa_f = k0 gamma_out`
and `a^2 = A_ref`:

```
re lambda / (a^2 gamma_out)  =  2 pi re / (kappa_f A_ref)
```

so

```
E_kin = (2 pi i re / kappa_f) * (1 / A_ref) * F_hkL / (exp(i Qz c) - 1)
```

which is Eq. (1.33) with the atomic lattice sum inserted. The agreement is
exact in magnitude and in the factor `i`, with no residual `k0^2 / (4 pi)`;
the apparent sign difference is entirely the `kappa_f = -kz,0,f` convention of
Sec. 0 above. This independently confirms invariant 3 (Born-Fresnel amplitude)
against a result that has itself been checked against dynamical theory.

Two convention differences, both harmless but both worth respecting:

- **Origin.** `1/(exp(i Qz c) - 1)` versus `1/(1 - exp(-i Qz c))` differ by one
  factor `exp(-i Qz c)`, that is, a one-cell shift of the unit-cell origin.
  Because `r_0` (Eq. (1.33) at `h = 0`) is added coherently to the unperturbed
  `r0` of Eq. (4.1), that phase is not free.
- **Exit angle.** Kaganer's normalizing wavevector is the *exit* normal
  component. On the specular rod it equals the incident one, so the design note
  is correct as written; but `kz,0,f` in Eq. (1.33) must be the genuine exit
  value off-specular, as Eq. (1.30) requires.

### Why Kaganer never needs the subtraction

All his rods are non-specular (11L, 1-3L, 22L). His zeroth-order reference
carries `chi_0`, so on a specular rod his `E_kin` would double count exactly the
term the design note keeps inside `F_0`. His silence is consistent with, and is
external evidence for, invariant 6 (off-specular invariance) and its
justification in Eq. (1.27).

### Caveat he raises that the design note originally did not

The mathematical surface plane must not bisect an atom. The kinematical and
dynamical calculations otherwise agree only down to about `1e-7` of the peak
intensity. For GaAs(001), atoms must be placed at fractional depths
`1/8, 3/8, 5/8, 7/8` rather than `0, 1/4, 1/2, 3/4`. This bears directly on the
design note, because the reference term of Eq. (3.5) is a sharp step at the
slab boundary while the atomic part of `F_0` is a sum of smooth atomic clouds:
the cancellation is exact only if the top cell's electron density lies entirely
inside the slab. This is the "within the atomic-truncation accuracy" clause of
invariant 5.

## 2. Half-space transform and the Fresnel limit

The `1/(i Qz)` regularization of the semi-infinite integral, Eq. (3.4), is
standard in the reflectivity literature, reached by integration by parts rather
than by a convergence factor. Vaknin, *X-ray diffraction techniques for liquid
surfaces and interfaces* (arXiv:cond-mat/0101142), Eq. (41):

```
integral rho(z) exp(i Qz z) dz = -(1 / (i Qz)) integral (d rho / dz) exp(i Qz z) dz
```

For a sharp step of density `rho` this gives amplitude `rho / (i Qz)` and hence
`R_BA = 16 pi^2 re^2 rho^2 / Qz^4 = (Qc / 2 Qz)^4`, that is,
`|r_BA| = 4 pi re rho / Qz^2`.

This is exactly Eq. (1.34), including its positive sign, which the design note
derives directly rather than importing. The sign is also consistent with the
exact Fresnel expansion `r = (kz0 - kz1) / (kz0 + kz1) ~ + delta k0^2 / (2 kz0^2)`
for `n < 1`. Invariant 3 is this check.

## 3. Contrast relative to the reference is the orthodox DWBA rule

- **Vineyard**, PRB 26, 4146 (1982). Origin of the surface DWBA: the reference
  is a homogeneous dielectric slab carrying the crystal's *average*
  polarizability, and the resulting distorted wave then illuminates the actual
  crystal. The design note's split, Eq. (1.9), is Vineyard's split.
- **Sinha, Sirota, Garoff & Stanley**, PRB 38, 2297 (1988). Reference is the
  ideal sharp interface; perturbation is the difference from it. The specular
  channel does receive a first-order DWBA correction, yielding Nevot-Croce.
  Precedent for Eq. (4.1), `r_total = r0 + r_0`.
- **Holy & Baumbach**, PRB 49, 10668 (1994); Pietsch, Holy & Baumbach,
  *High-Resolution X-Ray Scattering* (Springer 2004). Perturbation
  `delta chi = chi(r) - <chi(z)>`: specular in the reference, diffuse in the
  perturbation.
- **Zhou & Chen**, Phys. Rep. 257, 223 (1995). One-dimensional specular DWBA
  with the perturbation taken relative to a reference profile; the strict 1-D
  analogue of Eq. (1.9).

## 4. Strongest conceptual corroboration

Standard dynamical theory expands `chi(r) = sum_H chi_H exp(i H.r)`. The `H = 0`
term produces refraction and Fresnel specular reflection; `H != 0` produces
Bragg and CTR scattering.

This is precisely the structure of Eqs. (1.26)-(1.28): the design note's Fourier
expansion of `delta n^2` along the surface lattice *is* that expansion, and the
statement that the prepared reference `n0^2(z)` enters the `h = 0` coefficient
and no other is the statement that `chi_0` is an `H = 0` object. The reference
term of Eq. (3.5), `A_ref rho_bar_f / (i Qz)`, is the `Qz -> 0` divergence of
the kinematic lattice sum,

```
F_uc / (1 - exp(-i Qz c))  ->  A_ref c rho_bar_f / (i Qz c)
```

as Eq. (3.12) makes explicit. So "the delta, beta terms live inside `F_0`" is
the DWBA restatement of standard dynamical-theory bookkeeping. This framing is
harder to argue with than the half-space-transform argument alone and is worth
leading with in discussion.

## 5. Provenance of the 1D stratified Green function

Renaud's Eq. (100) quotes the one-dimensional Green function of a stratified
medium, including its factor `1/kz`, without derivation, citing Morse &
Feshbach for the general construction plus several reflectivity-specific
applications. Those sources are:

| Work | Role |
| --- | --- |
| Morse & Feshbach, *Methods of Theoretical Physics* (1953) | The general method: two homogeneous solutions satisfying the boundary condition at each end, divided by their Wronskian |
| Zhou & Chen, Phys. Rep. 257, 223 (1995) | Review-length derivation for X-ray/neutron reflectometry; the most pedagogical source |
| Andreev, Michette & Renwick, J. Mod. Opt. 35, 1667 (1988) | "Reflectivity and angular spectrum of scattered radiation" - the angular-spectrum route |
| Sears, PRB 48, 17477 (1993) | Same construction for neutron reflection |
| Feranchuk, PRB 52, 9214 (1995) | Graded-interface case worked through in detail |
| Maradudin & Mills, PRB 11, 1392 (1975) | Earliest application to EM scattering from a rough semi-infinite medium |

Section 1.5 of the design note now carries out that construction explicitly, in
orGUI's own conventions, so the derivation no longer depends on any of these.
The key intermediate result, Eq. (1.23), is that the Wronskian equals
`2 i kz,0` *independently of the prepared stack* - the `r0` terms cancel
identically. That is the origin of every factor `1/kz` in this subject, and it
is checkable in code: invariant 13.

Note also that Renaud's Eq. (100) is the `q_par = 0` specialization. The
general-`q` form, Eq. (1.15), is what a CTR at `Q_par = h != 0` requires, and it
is obtained by the same construction with no extra assumptions; see Eq. (1.30).

## 6. What the far-field route costs

Renaud's Eqs. (78)-(79) obtain the Green function by expanding
`|r - r'| ~ r - kf.r'/k0`. That single substitution does two jobs with very
different error budgets:

- **Phase.** The neglected term is quadratic, `~ k0 L^2 / R`, which for a
  grazing-incidence footprint of `L ~ 1 mm` to `1 cm` at `lambda ~ 1 A` and
  `R ~ 1 m` is `1e4` to `1e7`. The far-field form of the Green function is
  therefore never valid on the specular rod. Eq. (1.25) records this; Eq. (1.24)
  avoids it.
- **Polarization direction.** The neglected term is linear, `~ L / R ~ 1e-2`,
  and is not multiplied by `k0`. Fixing the exit polarization basis at the
  nominal `kf`, as Eq. (2.7) does, therefore remains accurate even where the
  phase treatment above is mandatory.

The two must not be conflated: only the first requires the boundary-value
construction of Sec. 1.5.

## 7. Recommended benchmarks

Since no paper does the specular DWBA-CTR directly, validation is numerical.

- **Holy & Fewster**, J. Appl. Cryst. 41, 18-26 (2008), *Crystal truncation rod
  X-ray scattering: exact dynamical calculation*. Slab decomposition with a
  Fresnel reflection/transmission matrix formalism, "correct both in the
  reciprocal-lattice points and between them", any polarization, non-coplanar.
  A slab-dynamical calculation of a bulk-terminated crystal along 00L contains
  Fresnel *and* CTR with no double counting, so it is the natural benchmark for
  Eq. (3.10).
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
| Andreev, Michette & Renwick, J. Mod. Opt. 35, 1667 (1988) | <https://doi.org/10.1080/09500348814551781> |
| Sears, PRB 48, 17477 (1993) | <https://doi.org/10.1103/PhysRevB.48.17477> |
| Feranchuk, PRB 52, 9214 (1995) | <https://doi.org/10.1103/PhysRevB.52.9214> |
| Maradudin & Mills, PRB 11, 1392 (1975) | <https://doi.org/10.1103/PhysRevB.11.1392> |
| Morse & Feshbach, Methods of Theoretical Physics (McGraw-Hill, 1953) | - |
| Holy & Fewster, J. Appl. Cryst. 41, 18 (2008) | <https://doi.org/10.1107/S0021889807049886> |
| Caticha, PRB 47, 76 (1993) | <https://doi.org/10.1103/PhysRevB.47.76> |
| Colella, PRB 43, 13827 (1991) | <https://doi.org/10.1103/PhysRevB.43.13827> |
| Takahashi & Nakatani, PRB 62, 3630 (2000) | <https://doi.org/10.1103/PhysRevB.62.3630> |
| Renaud, Lazzari & Leroy, Surf. Sci. Rep. 64, 255 (2009) | <https://doi.org/10.1016/j.surfrep.2009.07.002> |
| Vaknin, liquid surfaces review | <https://arxiv.org/abs/cond-mat/0101142> |

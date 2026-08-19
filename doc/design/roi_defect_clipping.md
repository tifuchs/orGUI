# Handling masked-pixel defects inside signal ROIs

Status: proposed, nothing implemented. Written 2026-08-19 after the ROI
normalization audit that followed the two `Corr_croi` bugs.

## The problem

An ROI-summed intensity is scaled from the valid pixels to the nominal ROI
area:

```python
croibg = (croi - (cpixel / bgpixel) * bgroi) * (roi_size / cpixel)
```

The background term is correct: `bgroi / bgpixel` is a density, and a density
extrapolates properly. The `roi_size / cpixel` factor is different. It applies
to the *background-subtracted* signal, and so assumes the net signal is spread
uniformly over the ROI. A Bragg or CTR reflection is not: it is concentrated
near the ROI centre, because the ROI is placed on the calculated reflection
position.

So the factor is a uniform-density extrapolation of something that is not
uniform, and it is applied silently. Nothing in the output says a point was
clipped, and nothing rejects one.

## What was measured

A 10x10 signal ROI with four surrounding background ROIs, a linear background
and a Gaussian peak of sigma = 2 px at the ROI centre, compared against the
true integrated intensity. `repair` is the existing
`Mask.PixelRepair` path at its default settings (4 px component, span 3,
radius 2, 6 valid neighbours).

| defect | current (rescale) | no rescale | repair + rescale |
| --- | --- | --- | --- |
| none | 0.0% | 0.0% | 0.0% |
| 2x2 on peak core | -10.9% | -14.5% | **-4.1%** |
| 2x2 off peak | **-0.9%** | -4.9% | -3.6% |
| 3x3 on core (exceeds repair limit) | -24.0% | -30.8% | -24.0% |
| 1 px column through peak | -12.6% | -21.4% | -12.6% |

A second geometry, a concentrated peak filling 16 of the 100 ROI pixels:

| defect | current (rescale) | no rescale |
| --- | --- | --- |
| 4 px off peak | +4.7% | +0.5% |
| 20 px stripe off peak | **+27.3%** | +1.9% |
| 20 px stripe through peak | -37.5% | -50.0% |

Conclusions from those numbers:

1. **Dropping the rescaling is not the fix.** For a peak that fills the ROI
   the rescaling is the better of the two in every case. It only misbehaves
   when the ROI is loose around a concentrated peak, where it inflates an
   off-peak defect by tens of percent.
2. **A defect on the peak core costs 10-30% and no post-hoc scaling recovers
   it.** Only repair does, and only while the defect is small enough.
3. **Repair silently does nothing when its limits are exceeded.** In the first
   attempt at this measurement, 4 px wide row defects were never repaired
   because their span of 4 exceeds `max_span = 3`; the run looked identical to
   the unrepaired one. There is no signal to the user that repair declined.

## Proposal

### 1. Flag, and optionally reject, clipped signal ROIs

Highest value for the least work, and the only item that removes the silent
part of the problem.

orGUI has a shortcut most integration codes do not: the signal ROI is centred
on the *calculated* reflection position, so whether a defect sits on the peak
is directly computable without fitting anything. A criterion of the form

* any masked pixel within `r` pixels of the ROI centre, or
* masked fraction of the signal ROI above `x`

separates the damaging cases (rows 2, 4, 5 of the first table) from the benign
ones (row 3) without a model.

Implementation notes:

* Store `masked_fraction` (and a `core_masked` flag) as counters next to
  `croi_pix`. The fraction is already derivable from the stored `croi_pix`,
  `hsize_corr` and `vsize_corr`, but nothing surfaces it and nothing acts on
  it.
* Rejection means writing `NaN`, which the existing
  `rod_mask = np.isfinite(croibg)` already drops from the saved rod.
* Both thresholds belong in the advanced ROI options next to the repair
  settings, and rejection should default to off so that existing behaviour is
  unchanged until a user opts in.

### 2. Inflate the uncertainty with the masked fraction

Cheap, needs no threshold, and composes with everything else.

`croibg_errors` is currently multiplied by the same `roi_size / cpixel`, which
captures the counting statistics lost with the masked pixels but not the
*systematic* error of the extrapolation. Adding a systematic term that grows
with the masked fraction -- at its simplest, a term proportional to
`masked_fraction * croibg` added in quadrature -- lets any downstream fit
down-weight clipped points automatically, with nothing to tune and no points
discarded.

This is the item to do alongside 1: together they remove the silent bias
without requiring any peak model.

### 3. Profile-weighted rescaling

Replace the uniform factor with

```
roi_size / cpixel   ->   sum(profile) / sum(profile[valid])
```

which reduces to the current factor for a flat profile and is correct for a
peaked one. It fixes both failure modes at once, and needs a profile, which
leads to item 4.

### 4. Profile recovery from neighbouring frames

The proper fix, and the only one that genuinely recovers an on-core defect.

In a rotation scan the same reflection crosses many frames and a defect clips
it on only some of them. Accumulating the 2D peak profile from the unclipped
frames gives both the profile for item 3 and a principled rejection criterion.
This is in essence what XDS and DIALS do with profile fitting. Estimated at a
day or so of work; it subsumes items 1 to 3.

### 5. Widen repair, but never for detector gaps

Repair works -- it is the only thing that recovered an on-core defect in the
table above -- and its current limits are conservative enough that roughly
only a 2x2 qualifies. Worth loosening for isolated dead pixels.

Detector gaps must stay excluded. They are long stripes where interpolation is
fabrication rather than estimation, and those points belong in the rejection
path of item 1 instead. The existing `use_pyfai_gaps` handling already refuses
them; that behaviour should not be relaxed.

Repair should also report when it declines a defect, so that the silent no-op
seen while measuring the table above cannot happen unnoticed.

## Recommendation

Do 1 and 2 together; they are roughly an afternoon and they remove the silent
bias. Treat 3 and 4 as a separate piece of work, justified only if the CTR data
needs it.

## Caveat on priority

On the FeReO4 scan 39 data the background is around 87% of the signal ROI
counts over much of the 1 0 L rod, so for those points the clipping bias is
small next to the background-subtraction uncertainty. It matters most for the
strong reflections near the Bragg peaks -- which is also where they carry the
most weight in a structure-factor fit.

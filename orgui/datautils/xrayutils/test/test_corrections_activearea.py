"""Regression tests for the illuminated active surface area.

:mod:`orgui.datautils.xrayutils.corrections.activearea` builds the area
:math:`A` that multiplies every integrated-intensity expression of E. Vlieg,
*J. Appl. Cryst.* **30** (1997) 532, in the two limits an experiment can be
in. See ``doc/design/ctr_structure_factor_scale.md`` finding F4.
"""

import numpy as np
import pytest

from orgui.datautils.xrayutils.corrections import activearea, geometry
from orgui.datautils.xrayutils.corrections.beamprofile import (
    gaussian_profile,
    top_hat_profile,
)

#: Horizontal beam size, vertical beam size and sample length, in meter.
W, H, L = 300e-6, 20e-6, 5e-3


def test_beam_limited_area_reduces_to_the_closed_form_for_a_top_hat():
    """``w L C_illum`` is ``w min(L, h/sin alpha)`` for a uniform beam.

    This identity is why the module offers no separate closed form. It also
    fixes the meaning of ``illuminated_area_fraction``: it is the mean of
    ``p(z)/p_max`` over the projected sample, so multiplying by the sample
    length and the beam width is all that turns it into an absolute area.
    """
    alpha = np.deg2rad(np.array([0.05, 0.1, 0.23, 0.5, 1.0, 2.0, 5.0, 15.0]))

    from_profile = activearea.beam_limited_area(alpha, W, L, top_hat_profile(H))
    closed_form = W * activearea.footprint_length(alpha, H, L)

    np.testing.assert_allclose(from_profile, closed_form, rtol=1e-12)


def test_a_real_beam_is_not_the_closed_form():
    """A Gaussian of the same width differs enough to matter.

    The closed form is a top-hat statement, not a general one: taking it for
    a focused beam is a percent-level error at high incidence and much worse
    where the beam spills off the sample.
    """
    alpha = np.deg2rad(np.array([0.05, 0.23, 2.0, 15.0]))

    gauss = activearea.beam_limited_area(alpha, W, L, gaussian_profile(H))
    top_hat = activearea.beam_limited_area(alpha, W, L, top_hat_profile(H))

    ratio = gauss / top_hat
    assert ratio.min() < 0.85, "grazing incidence must show the profile tails"
    assert ratio.max() > 1.06, "a flooded sample must show the peak density"


def test_footprint_length_is_clipped_by_the_sample():
    """Below the flooding angle the beam sets the length, above it the sample."""
    flooding = np.arcsin(H / L)
    alpha = np.array([flooding / 2.0, flooding * 2.0])

    length = activearea.footprint_length(alpha, H, L)

    np.testing.assert_allclose(length[0], L, rtol=1e-12)
    np.testing.assert_allclose(length[1], H / np.sin(alpha[1]), rtol=1e-12)


def test_slit_limited_area_is_the_z_axis_table_entry():
    """``s1 s2 / sin(delta)``, with the cosine that the z-axis mode drops.

    :func:`~.geometry.area_correction` is the tabulated row; this only turns
    it into an area and reinstates ``cos(alpha - beta_in)`` for a geometry
    where the incidence angle is not the angle to the surface.
    """
    delta = np.deg2rad(np.array([5.0, 17.0, 40.0]))
    s1, s2 = 200e-6, 1e-3

    z_axis = activearea.slit_limited_area(delta, s1, s2)

    np.testing.assert_allclose(
        z_axis, s1 * s2 * geometry.area_correction(delta), rtol=1e-12
    )
    # beta_in = alpha is the z-axis mode, so the cosine is one.
    alpha = np.deg2rad(0.6)
    np.testing.assert_allclose(
        activearea.slit_limited_area(delta, s1, s2, alpha=alpha, beta_in=alpha),
        z_axis,
        rtol=1e-12,
    )
    tilted = activearea.slit_limited_area(
        delta, s1, s2, alpha=alpha, beta_in=np.deg2rad(10.0)
    )
    np.testing.assert_allclose(
        tilted, z_axis / np.cos(alpha - np.deg2rad(10.0)), rtol=1e-12
    )


def test_the_two_limits_disagree_and_that_is_the_experimental_question():
    """They are alternatives, not factors, and only one can be right.

    The slit-limited area varies along a rod because it depends on delta;
    the beam-limited one does not. Which applies is set by the slit
    settings, and orGUI's open-slit configuration is the second.
    """
    delta = np.deg2rad(np.array([16.77, 16.86, 16.95]))
    alpha = np.deg2rad(np.full(delta.size, 0.6))

    slits = activearea.slit_limited_area(delta, 200e-6, 1e-3)
    beam = activearea.beam_limited_area(alpha, W, L, top_hat_profile(H))

    assert slits.max() / slits.min() > 1.005
    np.testing.assert_allclose(beam, beam[0], rtol=1e-12)


def test_sizes_must_be_physical():
    """A zero or negative size is a units mistake, not a value."""
    alpha = np.deg2rad(1.0)
    with pytest.raises(ValueError, match="beam_width must be positive"):
        activearea.beam_limited_area(alpha, 0.0, L, top_hat_profile(H))
    with pytest.raises(ValueError, match="sample_length must be positive"):
        activearea.beam_limited_area(alpha, W, -1.0, top_hat_profile(H))
    with pytest.raises(ValueError, match="needs a beam profile"):
        activearea.beam_limited_area(alpha, W, L, None)
    with pytest.raises(ValueError, match="beam_height must be positive"):
        activearea.footprint_length(alpha, 0.0, L)
    with pytest.raises(ValueError, match="slit_width must be positive"):
        activearea.slit_limited_area(np.deg2rad(17.0), 0.0, 1e-3)
    with pytest.raises(ValueError, match="beta_in needs alpha"):
        activearea.slit_limited_area(
            np.deg2rad(17.0), 200e-6, 1e-3, beta_in=np.deg2rad(1.0)
        )

"""Rocking and stationary integration must give the same structure factor.

`Issue #82 <https://github.com/tifuchs/orGUI/issues/82>`_: once every scaling
factor and experimental correction is accounted for, a rocking scan and a
stationary area-detector measurement of the same reflection must yield the
same :math:`|F_{hkl}|^2`; only the resolution differs.

These tests build simulated data for exactly that comparison. One rod with a
known :math:`|F_{hkl}|^2(l)` is measured twice -- as a series of rocking scans
and as a stationary *l* scan -- using
:func:`orgui.datautils.xrayutils.corrections.measurement.integrated_intensity` as
the forward model, i.e. Vlieg's equations 42 and 54 written out. The simulated
counts are then reduced twice: through
:mod:`orgui.datautils.xrayutils.corrections.measurement`, which must recover the
input, and through the two correction paths orGUI ships today
(:func:`orgui.app.peak1Dintegr._compute_rocking_integration` and
:mod:`orgui.app.integration_corrections`), whose disagreement is pinned to the
exact factor it is.

Everything that is not under test is left out deliberately. The polarization
factor and the footprint correction are absent from the simulation because
both depend only on the incidence and detector angles: they are identical for
the two modes at the same reflection, cancel from their ratio, and are
divided out per pixel long before the code under test here. Background
subtraction and error propagation are covered by ``test_peak1Dintegr.py``; a
flat background is included only so that the real aggregation code path is
exercised.

``doc/design/ctr_structure_factor_scale.md`` records the analysis.
"""

import numpy as np
import pytest

from orgui.app import integration_corrections as ic
from orgui.app.peak1Dintegr import _compute_rocking_integration
from orgui.datautils.xrayutils import HKLVlieg
from orgui.datautils.xrayutils.corrections import measurement as ii

#: Fixed incidence angle of the simulated z-axis scans, in radian.
ALPHA_IN = np.deg2rad(0.6)

#: Per-frame counting time and monitor of the two simulated measurements.
ROCKING_EXPOSURE, ROCKING_MONITOR = 0.4, 3.0
STATIONARY_EXPOSURE, STATIONARY_MONITOR = 2.0, 7.0

#: Flat background under every simulated rocking frame, in counts.
BACKGROUND = 50.0

#: Illuminated area (square meter) and flux density (photons / s / m^2).
ACTIVE_AREA, FLUX_DENSITY = 5.73e-7, 1e16


@pytest.fixture(scope="module")
def rod():
    """A (1, 0, l) rod of Pt(111) in the z-axis geometry at 17.7 keV.

    :returns: ``(l, alpha, delta, gamma, f2, wavelength, unitcell_area)`` with
        the angles in radian, ``f2`` the true :math:`|F_{hkl}|^2` in electron
        units squared, the wavelength in Angstrom and the unit-cell area in
        square Angstrom.
    """
    lattice = HKLVlieg.Lattice([2.7748, 2.7748, 6.7964], [90.0, 90.0, 120.0])
    ub = HKLVlieg.UBCalculator(lattice, 17.7)
    ub.defaultU_GID()
    angles = HKLVlieg.VliegAngles(ub)

    ell = np.linspace(0.4, 3.0, 12)
    hkl = np.vstack([np.ones_like(ell), np.zeros_like(ell), ell])
    computed = angles.anglesZmode(hkl, ALPHA_IN, fixed="in")
    alpha, delta, gamma = computed[:, 0], computed[:, 1], computed[:, 2]

    # A crystal-truncation-rod-like profile: sharp at the bulk Bragg
    # positions, weak in between, so a scale error that varies along the rod
    # is visible as a shape change and not only as an offset.
    f2 = 100.0 / (np.sin(np.pi * ell / 2.0) ** 2 + 0.02)
    return ell, alpha, delta, gamma, f2, ub.getLambda(), lattice.uc_area


def _simulate_stationary(rod):
    """Counts of a stationary *l* scan of the rod.

    :returns: The per-frame background-subtracted counts, in the units
        :meth:`orgui.app.orGUI.orGUI.integrateROI` hands to
        :mod:`orgui.app.integration_corrections`.
    """
    _, alpha, delta, gamma, f2, wavelength, uc_area = rod
    intensity = ii.integrated_intensity(
        f2,
        ii.STATIONARY,
        gamma=gamma,
        wavelength=wavelength,
        unitcell_area=uc_area,
        active_area=ACTIVE_AREA,
        flux_density=FLUX_DENSITY,
    )
    return intensity * STATIONARY_EXPOSURE * STATIONARY_MONITOR


def _simulate_rocking(rod, acceptance, axis=None, width=0.06):
    """Per-frame rocking curves of the rod.

    The rocking profile integrates to one over the rocking angle *in radian*,
    so the trapezoidal integral of the returned curves over the axis in
    degrees is the Vlieg integrated intensity times ``180/pi``.

    :param rod: The ``rod`` fixture.
    :param acceptance: Out-of-plane acceptance of the region of interest per
        rod point, in radian.
    :param axis: Rocking axis in degrees; a default symmetric axis is used
        when omitted.
    :param float width: Standard deviation of the rocking profile, in degrees.
    :returns: ``(axis, curves)`` with ``curves`` of shape ``(n_l, n_axis)``.
    """
    _, alpha, delta, gamma, f2, wavelength, uc_area = rod
    if axis is None:
        axis = np.linspace(-1.0, 1.0, 801)
    intensity = ii.integrated_intensity(
        f2,
        ii.ROCKING,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=acceptance,
        wavelength=wavelength,
        unitcell_area=uc_area,
        active_area=ACTIVE_AREA,
        flux_density=FLUX_DENSITY,
    )
    profile = np.exp(-0.5 * (axis / width) ** 2) / (width * np.sqrt(2.0 * np.pi))
    profile = profile * np.rad2deg(1.0)  # unit integral in radian
    counts = intensity * ROCKING_EXPOSURE * ROCKING_MONITOR
    return axis, counts[:, None] * profile[None, :] + BACKGROUND


def _roi_info(size):
    """Signal and background windows, in degrees, for every rod point.

    The background window sits ten standard deviations out, where the
    simulated profile has died away, so ``croibg`` is the signal alone.
    """
    return {
        "sig_1": {"from": np.full(size, -0.5), "to": np.full(size, 0.5)},
        "bg_1": {"from": np.full(size, -1.0), "to": np.full(size, -0.6)},
    }


def _orgui_rocking_f2(rod, acceptance, reduce=True, solid_angle_mean=None):
    """``F2_hkl`` as :mod:`orgui.app.peak1Dintegr` computes it.

    :param rod: The ``rod`` fixture.
    :param acceptance: Out-of-plane acceptance per rod point, in radian. It
        both shapes the simulated measurement and, when ``reduce`` is set, is
        divided back out.
    :param bool reduce: Pass the exposure/monitor divisor, the acceptance and
        the radian conversion, as :meth:`RockingPeakIntegrator.integrate`
        does. ``False`` reproduces the historical scale, which is what the
        distortion test below needs to compare against.
    """
    ell, alpha, delta, gamma, _, _, _ = rod
    axis, photon_curves = _simulate_rocking(rod, acceptance)
    curves = photon_curves
    if solid_angle_mean is not None:
        # As the per-pixel correction does when the switch is on: it scales
        # the diagnostic intensity, region by region. The independently
        # accumulated CTR branch retains the photon counts.
        curves = curves * np.asarray(solid_angle_mean)[:, None]
    shape = curves.shape
    lorentz = np.broadcast_to(
        (1.0 / (np.sin(delta) * np.cos(alpha) * np.cos(gamma)))[:, None], shape
    )
    rod_interception = np.broadcast_to(np.cos(gamma)[:, None], shape)
    extra = {}
    if reduce:
        extra = dict(
            C_norm=np.full(shape, ROCKING_EXPOSURE * ROCKING_MONITOR),
            detector_acceptance=acceptance,
            ctr_croibg_curves=photon_curves,
            ctr_croibg_errors_curves=np.sqrt(np.abs(photon_curves)),
            angle_unit="deg",
        )
    else:
        extra = dict(angle_unit="rad")  # no conversion: integrate as measured
    result = _compute_rocking_integration(
        ell,
        axis,
        curves,
        np.sqrt(np.abs(curves)),
        _roi_info(ell.size),
        {},
        True,
        False,
        C_Lor=lorentz,
        C_rod=rod_interception,
        **extra,
    )
    return result["F2_hkl"]


def _orgui_stationary_f2(rod, solid_angle_mean=None):
    """``F2_hkl`` as :mod:`orgui.app.integration_corrections` computes it.

    :param solid_angle_mean: When given, the region-mean solid-angle
        correction is applied only to the diagnostic intensity. The direct
        photon branch used for the structure factor remains unchanged.
    """
    ell, alpha, delta, gamma, _, _, _ = rod
    photon_counts = _simulate_stationary(rod)
    combined = 1.0 if solid_angle_mean is None else solid_angle_mean
    _, _, ctr_counts, ctr_errors = ic.pixel_correction_branches(
        photon_counts,
        np.sqrt(photon_counts),
        combined,
        1.0,
    )
    factors = ic.stationary_correction_factors(
        alpha,
        delta,
        gamma,
        use_lorentz=True,
        normalization=np.full(
            ell.size, STATIONARY_EXPOSURE * STATIONARY_MONITOR
        ),
    )
    ctr_intensity, ctr_errors = ic.apply_stationary_corrections(
        ctr_counts, ctr_errors, factors
    )
    return ic.structure_factor(ctr_intensity, ctr_errors, factors)[0]


def test_the_unified_reduction_recovers_one_structure_factor(rod):
    """Both simulated measurements reduce to the input ``|F|^2``.

    This is what issue #82 asks for, and it holds with the acceptance of the
    region of interest deliberately different at every point of the rod and
    the two measurements taken with different counting times and monitors.
    """
    ell, alpha, delta, gamma, f2, wavelength, uc_area = rod
    acceptance = np.deg2rad(0.35) * np.linspace(0.7, 1.6, ell.size)
    scale = dict(
        wavelength=wavelength,
        unitcell_area=uc_area,
        active_area=ACTIVE_AREA,
        flux_density=FLUX_DENSITY,
    )

    axis, curves = _simulate_rocking(rod, acceptance)
    result = _compute_rocking_integration(
        ell,
        axis,
        curves - BACKGROUND,
        np.sqrt(np.abs(curves)),
        _roi_info(ell.size),
        {},
        False,
        False,
    )
    f2_rocking = ii.structure_factor_squared(
        ii.normalized_intensity(
            result["croibg"],
            exposure_time=ROCKING_EXPOSURE,
            monitor=ROCKING_MONITOR,
            angle_unit="deg",
        ),
        ii.ROCKING,
        alpha=alpha,
        delta=delta,
        gamma=gamma,
        detector_acceptance=acceptance,
        **scale,
    )
    f2_stationary = ii.structure_factor_squared(
        ii.normalized_intensity(
            _simulate_stationary(rod),
            exposure_time=STATIONARY_EXPOSURE,
            monitor=STATIONARY_MONITOR,
        ),
        ii.STATIONARY,
        gamma=gamma,
        **scale,
    )

    np.testing.assert_allclose(f2_stationary, f2, rtol=1e-12)
    np.testing.assert_allclose(f2_rocking, f2, rtol=1e-6)
    np.testing.assert_allclose(f2_rocking, f2_stationary, rtol=1e-6)


def test_the_stationary_path_recovers_the_rod_up_to_one_constant(rod):
    """orGUI's stationary path already has the right shape.

    Its Lorentz factor is the published one and its normalization divides by
    the counting time and the monitor, so what it produces differs from the
    true ``|F|^2`` by a single number -- the absolute scale of issue #15 --
    and not by anything that varies along the rod.
    """
    f2 = rod[4]

    ratio = _orgui_stationary_f2(rod) / f2

    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-12)


def test_rocking_and_stationary_paths_agree(rod):
    """The two integration paths of orGUI now land on one scale.

    This is issue #82. The same rod, simulated as a rocking scan and as a
    stationary scan with *different* counting times and monitors, goes through
    the two real integration paths and comes out with the same ``F2_hkl``.

    This test used to characterize the gap between the paths, which was
    exactly ``exposure * monitor * Delta_gamma_in_degrees`` -- the rocking
    path applied neither the exposure/monitor normalization nor the
    out-of-plane acceptance, and integrated the rocking angle in degrees.
    Those three are the wiring this asserts is in place; the factor they used
    to leave behind is kept below as the thing that must not come back.
    """
    ell = rod[0]
    acceptance = np.full(ell.size, np.deg2rad(0.35))

    rocking = _orgui_rocking_f2(rod, acceptance)
    stationary = _orgui_stationary_f2(rod)

    np.testing.assert_allclose(rocking / stationary, 1.0, rtol=1e-6)

    # The historical scale, for contrast: without the reduction the rocking
    # path overshoots by the acceptance expressed in degrees times the
    # counting time and monitor it never divided out.
    historical = _orgui_rocking_f2(rod, acceptance, reduce=False)
    gap = ROCKING_EXPOSURE * ROCKING_MONITOR * np.rad2deg(acceptance)
    np.testing.assert_allclose(historical / stationary, gap, rtol=1e-6)


def test_a_resized_region_of_interest_no_longer_distorts_the_rocking_rod(rod):
    """The acceptance divisor removes a *shape* error, not just a scale one.

    :func:`orgui.app.ROIutils.calc_corrections` sizes regions of interest
    from the projected sample size and the parallax at each detector
    position, so their out-of-plane acceptance changes along a scan. Without
    the ``1/Delta_gamma`` divisor that change is carried straight into
    ``F2_hkl`` and the same rod comes out with a different shape; with it the
    rod is independent of how the regions were sized.

    Both halves are asserted, because the second is what makes the first
    worth having: a factor 2.3 of distortion across this rod, cured.
    """
    ell = rod[0]
    fixed = np.full(ell.size, np.deg2rad(0.35))
    resized = np.deg2rad(0.35) * np.linspace(0.7, 1.6, ell.size)

    # With the reduction: the resized regions give the same rod.
    np.testing.assert_allclose(
        _orgui_rocking_f2(rod, resized) / _orgui_rocking_f2(rod, fixed),
        1.0,
        rtol=1e-6,
    )

    # Without it: the region sizing leaks into the rod shape.
    carried = _orgui_rocking_f2(rod, resized, reduce=False) / _orgui_rocking_f2(
        rod, fixed, reduce=False
    )
    shape_change = carried / carried[0]
    np.testing.assert_allclose(
        shape_change, np.linspace(0.7, 1.6, ell.size) / 0.7, rtol=1e-6
    )
    assert shape_change.max() / shape_change.min() > 2.0


def test_the_stationary_path_assumes_a_slit_independent_active_area(rod):
    """What orGUI's missing area correction implies about the setup.

    orGUI divides out the numerical beam-profile factor but not Vlieg's
    :math:`C_\\mathrm{area} = 1/(\\sin\\delta\\cos(\\alpha-\\beta_
    \\mathrm{in}))`. That is right for an area detector with open post-sample
    slits, where the illuminated footprint and not the slits defines the
    active area -- the simulations above, which come out exact. It is wrong
    for a slit-limited setup, where the active area varies with
    :math:`\\delta` along the rod and the omission shows up as a
    rod-dependent shape error.
    """
    _, alpha, delta, gamma, f2, _, _ = rod
    slit_limited = 1.0 / (np.sin(delta) * np.cos(alpha - alpha))

    ratio = _orgui_stationary_f2(rod) * slit_limited / f2

    spread = ratio.max() / ratio.min() - 1.0
    assert spread > 1e-3, "delta must vary enough along the rod to see this"


def test_the_solid_angle_correction_does_not_reach_the_structure_factor(rod):
    """Switching it on changes the intensity but not ``F2_hkl``, in both modes.

    The detector solid-angle correction stays available because it is the
    right thing for a broad, non-rod feature, where a differential cross
    section is what is wanted. For a rod it is not: a region sum is already
    the complete angular integral, each pixel weighted by the solid angle it
    subtends. So the reduction divides it back out, and a structure factor is
    the same number whether or not the switch was on.

    The factor here varies along the rod, as a real one does -- the obliquity
    grows as the reflection moves up the detector -- so a leftover would show
    up as a shape error and not merely a scale.
    """
    ell = rod[0]
    acceptance = np.full(ell.size, np.deg2rad(0.35))
    solid_angle = 1.0 + 0.07 * np.linspace(0.0, 1.0, ell.size)

    np.testing.assert_allclose(
        _orgui_stationary_f2(rod, solid_angle_mean=solid_angle),
        _orgui_stationary_f2(rod),
        rtol=1e-12,
    )
    np.testing.assert_allclose(
        _orgui_rocking_f2(rod, acceptance, solid_angle_mean=solid_angle),
        _orgui_rocking_f2(rod, acceptance),
        rtol=1e-12,
    )

    # And the modes still agree with the correction enabled.
    np.testing.assert_allclose(
        _orgui_rocking_f2(rod, acceptance, solid_angle_mean=solid_angle)
        / _orgui_stationary_f2(rod, solid_angle_mean=solid_angle),
        1.0,
        rtol=1e-6,
    )

"""Regression tests for the z-axis geometrical correction factors.

These pin the entries of the z-axis column of Appendix A of the ANA/ROD
manual against literal transcriptions of the formulas, and pin the two
factors that orGUI's rocking-scan integration already applied before
:mod:`orgui.datautils.xrayutils.geometrycorrections` existed, so that
factoring them out did not change any result.
"""

import numpy as np
import pytest

from orgui.datautils.xrayutils import geometrycorrections as gc

#: Angles in radian, spanning grazing incidence to a normal scattering angle.
ALPHA = np.deg2rad(np.array([0.2, 0.5, 1.0, 3.0]))
DELTA = np.deg2rad(np.array([5.0, 20.0, 45.0, 80.0]))
GAMMA = np.deg2rad(np.array([1.0, 3.0, 10.0, 30.0]))


def test_rocking_scan_lorentz_matches_the_table():
    r"""``1 / (sin delta cos alpha cos gamma)``."""
    expected = 1.0 / (np.sin(DELTA) * np.cos(ALPHA) * np.cos(GAMMA))

    np.testing.assert_allclose(
        gc.lorentz_rocking_scan(DELTA, ALPHA, GAMMA), expected, rtol=1e-12
    )


def test_reflectivity_rocking_lorentz_matches_the_table():
    r"""``1 / sin(2 alpha)``."""
    np.testing.assert_allclose(
        gc.lorentz_reflectivity_rocking_scan(ALPHA),
        1.0 / np.sin(2.0 * ALPHA),
        rtol=1e-12,
    )


def test_stationary_lorentz_matches_the_table():
    r"""``1 / sin gamma``, the stationary-mode entry.

    This is a different factor from the rocking-scan one, and depends on the
    exit angle alone.
    """
    np.testing.assert_allclose(
        gc.lorentz_stationary(GAMMA), 1.0 / np.sin(GAMMA), rtol=1e-12
    )


def test_rod_interception_and_area_correction_match_the_table():
    r"""``cos gamma`` and ``1 / sin delta``."""
    np.testing.assert_allclose(gc.rod_interception(GAMMA), np.cos(GAMMA), rtol=1e-12)
    np.testing.assert_allclose(
        gc.area_correction(DELTA), 1.0 / np.sin(DELTA), rtol=1e-12
    )


def test_rocking_factors_reproduce_the_previous_implementation():
    """The formulas orGUI's rocking integration applied before the refactor.

    ``peak1Dintegr`` computed these inline; extracting them must not have
    changed a single value.
    """
    previous_th = np.abs(1 / (np.sin(DELTA) * np.cos(ALPHA) * np.cos(GAMMA)))
    previous_mu = np.abs(1 / np.sin(2 * ALPHA))
    previous_rod = np.cos(GAMMA)

    np.testing.assert_array_equal(
        gc.lorentz_rocking_scan(DELTA, ALPHA, GAMMA), previous_th
    )
    np.testing.assert_array_equal(
        gc.lorentz_reflectivity_rocking_scan(ALPHA), previous_mu
    )
    np.testing.assert_array_equal(gc.rod_interception(GAMMA), previous_rod)


def test_lorentz_factors_stay_positive_for_a_reversed_scan():
    """Integrated intensities are forced positive whichever way a scan runs."""
    assert np.all(gc.lorentz_rocking_scan(-DELTA, ALPHA, GAMMA) > 0)
    assert np.all(gc.lorentz_reflectivity_rocking_scan(-ALPHA) > 0)
    assert np.all(gc.lorentz_stationary(-GAMMA) > 0)


def test_lorentz_factor_dispatches_on_the_mode():
    """The three factors are alternatives selected by measurement mode."""
    np.testing.assert_allclose(
        gc.lorentz_factor(gc.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA),
        gc.lorentz_rocking_scan(DELTA, ALPHA, GAMMA),
    )
    np.testing.assert_allclose(
        gc.lorentz_factor(gc.REFLECTIVITY_ROCKING, alpha=ALPHA),
        gc.lorentz_reflectivity_rocking_scan(ALPHA),
    )
    np.testing.assert_allclose(
        gc.lorentz_factor(gc.STATIONARY, gamma=GAMMA), gc.lorentz_stationary(GAMMA)
    )
    # The stationary and rocking factors genuinely differ; a mix-up would be
    # a silent scaling error in the structure factors.
    assert not np.allclose(
        gc.lorentz_factor(gc.STATIONARY, gamma=GAMMA),
        gc.lorentz_factor(gc.ROCKING, alpha=ALPHA, delta=DELTA, gamma=GAMMA),
    )


def test_lorentz_factor_rejects_unknown_modes_and_missing_angles():
    """A missing angle must fail loudly rather than default to something."""
    with pytest.raises(ValueError, match="unknown Lorentz mode"):
        gc.lorentz_factor("rscan", alpha=ALPHA)
    with pytest.raises(ValueError, match="missing gamma"):
        gc.lorentz_factor(gc.STATIONARY)
    with pytest.raises(ValueError, match="missing delta"):
        gc.lorentz_factor(gc.ROCKING, alpha=ALPHA, gamma=GAMMA)

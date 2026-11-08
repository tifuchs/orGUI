"""ROI estimator covariance rather than independent-count substitutions."""
import numpy as np

from ..corrections.roi import roi_difference_variance


def test_shared_donors_cancel_background_noise():
    """A repaired mean minus that same donor mean has identically zero noise."""
    # S=(X+Y)/2 and B=X+Y for independent Poisson mean10 donors.
    assert roi_difference_variance(5., 20., 10., .5) == 0.
    # S=X+Y+(X+Y)/2, B=Y: coefficients are1.5X+.5Y.
    assert roi_difference_variance(45., 10., 15., 1.) == 25.


def test_variance_roundoff_does_not_hide_inconsistent_moments():
    """Only machine-roundoff negative results are clipped."""
    result = roi_difference_variance([1., 1.], [1., 1.],
                                     [1.+np.finfo(float).eps, 2.], 1.)
    assert result[0] == 0.
    assert np.isnan(result[1])

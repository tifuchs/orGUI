"""Regression tests for the CTR fitting optimizers.

Increment 1 of the CTR optimizer rework covers the callback ``set_errors``
arity repair. The fixtures here are shared with the later increments of
``doc/design/dwba_ctr_fitting_implementation_plan.md``.
"""

import unittest

import numpy as np

from .. import CTRopt, CTRplotutil, CTRresolution


def _angles(gamma, delta=None):
    """Return six-circle angle records in rad for one CTR."""
    gamma = np.asarray(gamma, dtype=np.float64)
    zeros = np.zeros_like(gamma)
    if delta is None:
        delta = zeros
    else:
        delta = np.asarray(delta, dtype=np.float64)
    return np.rec.fromarrays(
        [zeros, delta, gamma, zeros, zeros, zeros],
        names="alpha,delta,gamma,omega,chi,phi",
    )


class FitCrystal:
    """Minimal crystal interface used to drive the optimizers.

    The model is ``F(l) = offset + slope*l + sum_i p_i * l**(i+2)``, in
    arbitrary structure-factor units. ``offset`` and ``slope`` are written by
    fit callbacks, so callback parameters change the model independently of the
    crystal parameters and the weighted-residual Jacobian keeps full rank.

    :param sequence parameters:
        Initial crystal fit parameters.
    :param sequence names:
        Optional fit parameter names; defaults to ``xtal_0``, ``xtal_1``, ...
    """

    def __init__(self, parameters=(1.0,), names=None):
        self.parameters = np.asarray(parameters, dtype=np.float64)
        if names is None:
            names = [f"xtal_{i}" for i in range(self.parameters.size)]
        self.fitparnames = list(names)
        self.priors = []
        self.errors = None
        self.error_calls = []
        self.offset = 1.0
        self.slope = 0.5

    def F(self, h, k, l):  # noqa: N802,E741
        """Return the calculated structure factor at the given HKL."""
        lvalues = np.asarray(l, dtype=np.float64)
        value = self.offset + self.slope * lvalues
        for i, parameter in enumerate(self.parameters):
            value = value + parameter * lvalues ** (i + 2)
        return value

    def getStartParamAndLimits(self):  # noqa: N802
        """Return the start parameters and their lower and upper limits."""
        return (
            self.parameters.copy(),
            np.full(self.parameters.size, 0.1),
            np.full(self.parameters.size, 10.0),
        )

    def getInitialParameters(self):  # noqa: N802
        """Return the current crystal fit parameters."""
        return self.parameters.copy()

    def setParameters(self, parameters):  # noqa: N802
        """Set the crystal fit parameters."""
        self.parameters = np.asarray(parameters, dtype=np.float64)

    def setFitErrors(self, errors):  # noqa: N802
        """Store the crystal parameter errors and record the call."""
        self.errors = np.asarray(errors, dtype=np.float64)
        self.error_calls.append(self.errors)


def _fixture_ctrs():
    """Return one specular and two non-specular rods with distinct weights.

    The constructor's ``|h| + |k|`` sort places the specular rod first, which
    is the only reason the shared-scale branch of ``CTROptAngleCorrection``---
    which appends individually scaled specular entries ahead of the shared
    ones---still emits results in ``fit.CTRs`` order.
    """
    specification = (
        ((0.0, 0.0), [3.0, 5.0, 9.0], 1.0),
        ((1.0, 0.0), [2.0, 4.0, 8.0], 2.0),
        ((2.0, 0.0), [1.0, 3.0, 7.0], 0.5),
    )
    lvalues = np.array([0.5, 1.5, 2.5])
    ctrs = []
    for hk, sfI, weight in specification:
        sfI = np.asarray(sfI, dtype=np.float64)
        err = 0.1 * (1.0 + np.arange(sfI.size))
        ctr = CTRplotutil.CTR(hk, lvalues.copy(), sfI, err)
        ctr.weight = weight
        ctr.angles = _angles(0.1 * np.arange(sfI.size), 0.2 * np.arange(sfI.size))
        ctrs.append(ctr)
    return CTRplotutil.CTRCollection(ctrs)


def _register_callback(optimizer, name, n_pars):
    """Register a callback and return it with a log of its parameter calls.

    :param CTRopt.CTROptimizer optimizer:
        Optimizer to register the callback on.
    :param str name:
        Callback name.
    :param int n_pars:
        One to drive ``offset`` alone, two to drive ``offset`` and ``slope``.
    :returns:
        The registered callback and the list recording its parameter calls.
    :rtype: tuple[CTRopt.FitCallback, list]
    """
    calls = []
    init = [1.0, 0.5][:n_pars]

    def apply(xtal, x):
        calls.append(np.copy(x))
        xtal.offset = x[0]
        if n_pars > 1:
            xtal.slope = x[1]

    optimizer.register_fit_callback(
        apply, [0.1] * n_pars, [5.0] * n_pars, init, name=name
    )
    return optimizer.callbacks[0], calls


class TestCallbackErrors(unittest.TestCase):
    """Increment 1: the callback ``set_errors`` arity repair."""

    def test_statistics_sets_callback_errors(self):
        """``statistics`` must reach a registered callback without raising.

        ``CTROptAngleCorrection.set_errors`` used to call
        ``cb.set_errors(self.xtal, slice)`` against the one-argument
        ``FitCallback.set_errors``, so requesting statistics for any fit with a
        registered callback raised ``TypeError``.
        """
        optimizer = CTRopt.CTROptAngleCorrection(FitCrystal(), _fixture_ctrs())
        callback, calls = _register_callback(optimizer, "cb", 2)
        optimizer.prepareFit()
        x = optimizer.get_parameters()
        self.assertEqual(x.size, 3)

        stat = optimizer.statistics(x)

        self.assertEqual(stat["nodatapoints"], 9)
        self.assertEqual(optimizer.errors.size, 3)
        self.assertTrue(np.all(np.isfinite(optimizer.errors)))
        np.testing.assert_allclose(callback.errors, optimizer.errors[:2])
        np.testing.assert_allclose(optimizer.xtal.errors, optimizer.errors[2:])
        self.assertTrue(calls)

    def test_set_errors_slices_by_parameter_layout(self):
        """Each consumer receives its own segment of the error vector.

        The layout is resolution, callbacks in ``self.callbacks`` order, angle
        correction, then crystal. The two angle-correction entries are dropped:
        the subclass keeps no angle-correction error state.
        """
        optimizer = CTRopt.CTROptAngleCorrection(FitCrystal(), _fixture_ctrs())
        optimizer.useAnglecorr = True
        optimizer.fit_resolution(
            CTRresolution.BoxResolution(0.1, 0.0, 0.0),
            lower_bounds=[0.0, 0.0, 0.0],
            higher_bounds=[1.0, 1.0, 1.0],
        )
        first, _ = _register_callback(optimizer, "cb_first", 2)
        second, _ = _register_callback(optimizer, "cb_second", 1)
        optimizer.prepareFit()

        self.assertEqual(optimizer.callbacks, [second, first])
        self.assertEqual(optimizer.get_parameters().size, 9)

        errors = 100.0 + np.arange(9, dtype=np.float64)
        optimizer.set_errors(errors)

        np.testing.assert_allclose(optimizer.resolution_errors, errors[:3])
        np.testing.assert_allclose(second.errors, errors[3:4])
        np.testing.assert_allclose(first.errors, errors[4:6])
        np.testing.assert_allclose(optimizer.xtal.errors, errors[8:])
        self.assertEqual(len(optimizer.xtal.error_calls), 1)


if __name__ == "__main__":
    unittest.main()

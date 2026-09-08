"""Regression tests for the CTR fitting optimizers.

Increment 1 of the CTR optimizer rework covers the callback ``set_errors``
arity repair. The fixtures here are shared with the later increments of
``doc/design/dwba_ctr_fitting_implementation_plan.md``.
"""

import unittest

import numpy as np

from .. import CTRcalc, CTRopt, CTRplotutil, CTRresolution, CTRsymmetry


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


def _register_callback(optimizer, name, n_pars, inert=False):
    """Register a callback and return it with a log of its parameter calls.

    :param CTRopt.CTROptimizer optimizer:
        Optimizer to register the callback on.
    :param str name:
        Callback name.
    :param int n_pars:
        One to drive ``offset`` alone, two to drive ``offset`` and ``slope``.
    :param bool inert:
        Record calls without touching the crystal. Use this with a real
        :class:`CTRcalc.SXRDCrystal`, which has no such attributes.
    :returns:
        The registered callback and the list recording its parameter calls.
    :rtype: tuple[CTRopt.FitCallback, list]
    """
    calls = []
    init = [1.0, 0.5][:n_pars]

    def apply(xtal, x):
        calls.append(np.copy(x))
        if inert:
            return
        xtal.offset = x[0]
        if n_pars > 1:
            xtal.slope = x[1]

    optimizer.register_fit_callback(
        apply, [0.1] * n_pars, [5.0] * n_pars, init, name=name
    )
    return optimizer.callbacks[0], calls


def _wyckoff_surface_cell(name="surface"):
    """Return a two-atom surface unit cell carrying one Wyckoff site.

    The symmetry metadata couples atom 0's ``x`` coordinate to the site
    variable ``u``, which is what ``addWyckoffParameter`` fits. Lattice
    constants are in Angstrom and angles in degrees.
    """
    unitcell = CTRcalc.UnitCell([3.0, 4.0, 5.0], [90.0, 90.0, 90.0], name=name)
    unitcell.addAtom("C", [0.25, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
    unitcell.addAtom("O", [0.75, 0.5, 0.5], 0.1, 0.1, 1.0, layer=1)
    unitcell.layerpos[1.0] = 0.0
    unitcell.symmetry_metadata = CTRsymmetry.SurfaceSymmetryModel(
        CTRsymmetry.SurfaceCellSpec(
            (3.0, 4.0, 5.0),
            (90.0, 90.0, 90.0),
            np.identity(3),
            layer_origins=(0.0,),
        ),
        (
            CTRsymmetry.WyckoffSiteSpec(
                site_id="C_1",
                element="C",
                wyckoff_label="1a",
                coordinates=(),
                representative_parent_fractional=(0.25, 0.0, 0.0),
                variables={"u": 0.25},
                occ=1.0,
                iDW=0.1,
                oDW=0.1,
            ),
        ),
        atoms=[
            CTRsymmetry.GeneratedWyckoffAtom(
                atom_index=0,
                element="C",
                site_id="C_1",
                wyckoff_label="1a",
                parent_fractional=np.array([0.25, 0.0, 0.0]),
                surface_fractional=np.array([0.25, 0.0, 0.0]),
                layer=1,
                couplings=(
                    CTRsymmetry.WyckoffCoupling(
                        atom_index=0,
                        coordinate="x",
                        variable="u",
                        constant=0.0,
                        factor=1.0,
                        site_id="C_1",
                    ),
                ),
                site_couplings=(
                    CTRsymmetry.WyckoffSiteCoupling(
                        atom_index=0,
                        coordinate="x",
                        axis="x",
                        factor=1.0,
                        site_id="C_1",
                    ),
                ),
            ),
        ],
    )
    return unitcell


def _parameter_crystal(wyckoff="parameter"):
    """Return a crystal exercising every populated ``SXRDCrystal`` bucket.

    One parameter of each kind, in the order the crystal reports them:

    - ``coupled_u``: crystal-level coupled parameter (``parameters["coupled"]``)
    - ``w_surface``: unit-cell weight (``parameters["weight"]``)
    - ``abs_z``: unit-cell absolute parameter (``parameters["absolute"]``)
    - ``rel_x``: unit-cell relative parameter (``parameters["relative"]``)
    - ``wyck``: unit-cell Wyckoff site variable or site shift

    ``parameters["domain"]`` stays empty: nothing in the package appends to it,
    so a domain parameter can only reach a crystal through deserialization,
    and ``getFitErrors`` raises for one.

    :param str wyckoff:
        ``"parameter"`` fits the site variable ``u`` through
        ``addWyckoffParameter``; ``"shift"`` fits a rigid site shift along
        ``x`` through ``addWyckoffShift``.
    :returns: The crystal and its surface unit cell.
    :rtype: tuple[CTRcalc.SXRDCrystal, CTRcalc.UnitCell]
    """
    if wyckoff not in {"parameter", "shift"}:
        raise ValueError("wyckoff must be 'parameter' or 'shift'")
    bulk = CTRcalc.UnitCell([3.0, 4.0, 5.0], [90.0, 90.0, 90.0], name="bulk")
    bulk.addAtom("C", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0, layer=1)
    bulk.layerpos[1.0] = 0.0
    surface = _wyckoff_surface_cell()
    crystal = CTRcalc.SXRDCrystal(bulk, surface)
    crystal.setEnergy(20.0)

    surface.addFitParameter(([1], "z"), limits=(-1.0, 1.0), name="abs_z")
    surface.addRelParameter(([1], "x"), [1.0], (-0.5, 0.5), name="rel_x")
    if wyckoff == "parameter":
        surface.addWyckoffParameter("C_1", "u", limits=(0.1, 0.4), name="wyck")
    else:
        surface.addWyckoffShift("C_1", "x", limits=(-0.2, 0.2), name="wyck")
    crystal.addWeightParameter(
        ["surface"], [1.0], limits=(0.0, 2.0), name="w_surface"
    )
    crystal.addWyckoffParameter(
        {"surface": ("C_1", "u")}, limits=(0.1, 0.4), name="coupled_u"
    )
    return crystal, surface


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


class TestCallbackBounds(unittest.TestCase):
    """Increment 1: the ``FitCallback`` bounds ordering repair."""

    def test_callback_bounds_keep_the_lower_upper_order(self):
        """Callback bounds must reach the optimizer as ``(lower, upper)``.

        ``FitCallback.__init__`` used to store ``(bounds_high, bounds_low)``,
        and ``prepareFit`` prepends ``bounds[0]`` to the lower array, so the
        assembled problem had a lower bound above its upper bound and a solver
        rejected it.
        """
        optimizer = CTRopt.CTROptAngleCorrection(FitCrystal(), _fixture_ctrs())
        callback, _ = _register_callback(optimizer, "cb", 2)

        np.testing.assert_allclose(callback.bounds[0], [0.1, 0.1])
        np.testing.assert_allclose(callback.bounds[1], [5.0, 5.0])

        optimizer.prepareFit()
        lower, higher = optimizer.get_bounds()

        np.testing.assert_allclose(lower[:2], [0.1, 0.1])
        np.testing.assert_allclose(higher[:2], [5.0, 5.0])
        self.assertTrue(np.all(lower <= higher))

    def test_callback_bound_size_errors_name_the_failing_side(self):
        """A wrong-size bound sequence must name the side that is wrong."""
        with self.assertRaisesRegex(ValueError, "lower bounds"):
            CTRopt.FitCallback(lambda xtal, x: None, [0.1, 0.2], [5.0], [1.0])
        with self.assertRaisesRegex(ValueError, "upper bounds"):
            CTRopt.FitCallback(lambda xtal, x: None, [0.1], [5.0, 6.0], [1.0])


class TestCrystalParameterLayout(unittest.TestCase):
    """The optimizer's layout against the real ``SXRDCrystal`` parameter API.

    These use a real crystal rather than a test double because the parameter
    kinds---crystal-level coupled and weight parameters, and the unit cell's
    absolute, relative, and Wyckoff parameters---are reordered by
    ``getStartParamAndLimits``, and ``setFitErrors`` undoes that ordering. A
    double cannot exercise that.
    """

    def test_optimizer_bounds_follow_crystal_parameter_limits(self):
        """Bounds reach the optimizer in order and keep lower below upper."""
        crystal, _ = _parameter_crystal()
        optimizer = CTRopt.CTROptAngleCorrection(crystal, _fixture_ctrs())
        optimizer.prepareFit()

        start, lower, higher = optimizer.xtal.getStartParamAndLimits()
        np.testing.assert_allclose(optimizer.get_bounds()[0], lower)
        np.testing.assert_allclose(optimizer.get_bounds()[1], higher)
        np.testing.assert_allclose(optimizer.get_parameters(), start)
        self.assertTrue(np.all(lower <= higher))

    def test_fitparnames_cover_every_crystal_parameter_kind(self):
        """All five populated buckets appear once, in crystal order."""
        for wyckoff in ("parameter", "shift"):
            with self.subTest(wyckoff=wyckoff):
                crystal, _ = _parameter_crystal(wyckoff)
                optimizer = CTRopt.CTROptAngleCorrection(crystal, _fixture_ctrs())
                optimizer.prepareFit()

                self.assertEqual(
                    optimizer.fitparnames,
                    ["coupled_u", "w_surface", "abs_z", "rel_x", "wyck"],
                )
                self.assertEqual(
                    len(optimizer.fitparnames), optimizer.get_parameters().size
                )
                self.assertEqual(optimizer.xtal.parameters["domain"], [])

    def test_error_slices_reach_each_crystal_parameter_kind(self):
        """``set_errors`` routes one error to each parameter object.

        The crystal reorders the flat vector through
        ``fit_metadata_cache["par_idx_sortarray"]`` before splitting it per
        unit cell, so a per-kind assertion is the only way to see the whole
        round trip.
        """
        for wyckoff in ("parameter", "shift"):
            with self.subTest(wyckoff=wyckoff):
                crystal, _ = _parameter_crystal(wyckoff)
                optimizer = CTRopt.CTROptAngleCorrection(crystal, _fixture_ctrs())
                optimizer.prepareFit()

                optimizer.set_errors(np.array([0.11, 0.22, 0.33, 0.44, 0.55]))

                surface = optimizer.xtal.uc_surface_list[0]
                recorded = {
                    par.name: par.error
                    for par in (
                        optimizer.xtal.parameters["coupled"]
                        + optimizer.xtal.parameters["weight"]
                        + surface.parameters["absolute"]
                        + surface.parameters["relative"]
                    )
                }
                self.assertEqual(
                    recorded,
                    {
                        "coupled_u": 0.11,
                        "w_surface": 0.22,
                        "abs_z": 0.33,
                        "rel_x": 0.44,
                        "wyck": 0.55,
                    },
                )

    def test_weight_errors_follow_the_coupled_block(self):
        """A weight parameter must not read a coupled parameter's error.

        ``setFitParameters`` and ``setLimits`` start the weight block at
        ``number_coupled`` in the reordered vector; ``setFitErrors`` started it
        at zero, so with both kinds present every weight parameter reported the
        first coupled parameter's error instead of its own.
        """
        crystal, _ = _parameter_crystal()
        crystal.setFitErrors([0.11, 0.22, 0.33, 0.44, 0.55])

        coupled = crystal.parameters["coupled"][0]
        weight = crystal.parameters["weight"][0]
        self.assertEqual(coupled.name, "coupled_u")
        self.assertEqual(weight.name, "w_surface")
        self.assertAlmostEqual(coupled.error, 0.11)
        self.assertAlmostEqual(weight.error, 0.22)
        np.testing.assert_allclose(
            crystal.werrors[weight.indices], np.abs(weight.factors) * 0.22
        )

    def test_crystal_keeps_its_own_slice_under_a_full_parameter_layout(self):
        """Resolution, callbacks and angle correction must not shift the tail.

        With all three active the crystal owns only the trailing five entries;
        everything earlier belongs to the other consumers.
        """
        crystal, _ = _parameter_crystal()
        optimizer = CTRopt.CTROptAngleCorrection(crystal, _fixture_ctrs())
        optimizer.useAnglecorr = True
        optimizer.fit_resolution(
            CTRresolution.BoxResolution(0.1, 0.0, 0.0),
            lower_bounds=[0.0, 0.0, 0.0],
            higher_bounds=[1.0, 1.0, 1.0],
        )
        callback, _ = _register_callback(optimizer, "cb", 2, inert=True)
        optimizer.prepareFit()

        self.assertEqual(optimizer.get_parameters().size, 12)
        errors = 100.0 + np.arange(12, dtype=np.float64)
        optimizer.set_errors(errors)

        np.testing.assert_allclose(optimizer.resolution_errors, errors[:3])
        np.testing.assert_allclose(callback.errors, errors[3:5])
        np.testing.assert_allclose(
            optimizer.xtal.parameters["coupled"][0].error, errors[7]
        )
        np.testing.assert_allclose(
            optimizer.xtal.uc_surface_list[0].parameters["absolute"][0].error,
            errors[9],
        )


if __name__ == "__main__":
    unittest.main()

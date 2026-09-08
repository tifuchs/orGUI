"""Regression tests for the CTR fitting optimizers.

Increments 1--7 of the CTR optimizer rework cover the callback error path,
parameter-name layout, value-preserving optimizer class split, and legacy
objective characterization, calculation adapters, prediction-side analytical
scales, and explicit scale policies. The fixtures here are shared with later
increments of ``doc/design/dwba_ctr_fitting_implementation_plan.md``.
"""

import unittest
from unittest import mock

import numpy as np
import pytest

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


class _ConstantFitCrystal(FitCrystal):
    """Fit double returning a unit structure factor at every coordinate."""

    def F(self, h, k, l):  # noqa: N802,E741
        """Return unit structure factors matching the supplied L coordinates."""
        return np.ones_like(np.asarray(l, dtype=np.float64))


class _VanishingFitCrystal(FitCrystal):
    """Fit double whose prediction vanishes only when its parameter is zero."""

    def F(self, h, k, l):  # noqa: N802,E741
        """Return ``theta * (l + 1)`` at the supplied L coordinates."""
        return self.parameters[0] * (np.asarray(l, dtype=np.float64) + 1.0)


class _ConstraintFitCrystal(FitCrystal):
    """Fit double exposing the surface displacement constraint interface."""

    def __init__(self):
        super().__init__()
        self.uc_bulk = mock.Mock()
        self.uc_bulk.basis = np.array(
            [[0.0, 0.0, 0.0, 0.0, 0.20, 0.30]], dtype=np.float64
        )
        self._surface_basis = np.array(
            [
                [0.0, 0.0, 0.0, 2.0, 0.05, 0.08],
                [0.0, 0.0, 0.0, 1.0, 0.10, 0.12],
            ],
            dtype=np.float64,
        )

    def getSurfaceDWConstraintEnable(self):  # noqa: N802
        """Enable both surface rows for displacement constraints."""
        return np.array([True, True])

    def getSurfaceBasis(self):  # noqa: N802
        """Return the surface basis consumed by constraint calculation."""
        return self._surface_basis


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


def _unit_model_ctrs(specification):
    """Build unit-error, unit-weight rods for constant-model regressions."""
    ctrs = []
    for hk, values in specification:
        observations = np.asarray(values, dtype=np.float64)
        ctr = CTRplotutil.CTR(
            hk,
            np.arange(observations.size, dtype=np.float64),
            observations,
            np.ones_like(observations),
        )
        ctr.weight = 1.0
        ctrs.append(ctr)
    return CTRplotutil.CTRCollection(ctrs)


def _register_callback(optimizer, name, n_pars, inert=False, parnames=None):
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
        apply,
        [0.1] * n_pars,
        [5.0] * n_pars,
        init,
        name=name,
        parnames=parnames,
    )
    return optimizer.callbacks[0], calls


class _FakePopulation:
    """Small pygmo population substitute for trace-layout tests."""

    def __init__(self, parameters, fitness):
        self._parameters = np.asarray(parameters, dtype=np.float64)
        self._fitness = np.asarray(fitness, dtype=np.float64).reshape(-1, 1)

    @property
    def champion_x(self):
        """Return the parameter vector with the lowest fitness."""
        return self._parameters[np.argmin(self._fitness[:, 0])]

    def get_x(self):
        """Return every population parameter vector."""
        return self._parameters

    def get_f(self):
        """Return every population fitness value."""
        return self._fitness


class _FakeIsland:
    """Small pygmo island substitute for trace-layout tests."""

    def __init__(self, population):
        self._population = population

    def get_population(self):
        """Return the island's population."""
        return self._population


class _FakeArchipelago(list):
    """Small pygmo archipelago substitute for trace-layout tests."""

    def get_champions_f(self):
        """Return the best fitness from every island."""
        return [
            [np.min(island.get_population().get_f()[:, 0])] for island in self
        ]


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


_LEGACY_FLAT_FCALC = np.array(
    [
        2.020833333333333,
        5.3888888888888875,
        11.451388888888888,
        1.5572519083969467,
        4.1526717557251915,
        8.82442748091603,
        1.17375,
        3.1300000000000003,
        6.65125,
    ]
)

_BASE_INCREMENT5_REFERENCE = {
    "weighted_residues": [
        9.791666666666671,
        -1.9444444444444402,
        -8.171296296296289,
        6.261403558598435,
        -1.0795523376893883,
        -3.8863884156817914,
        -1.2285980323116268,
        -0.45961940777125576,
        0.8220116331293595,
    ],
    "weighted_residues2": [
        95.87673611111119,
        3.780864197530847,
        66.77008316186543,
        39.205174523629125,
        1.165433249810623,
        15.104014917545623,
        1.5094531250000007,
        0.21124999999999983,
        0.6757031249999966,
    ],
    "residues": [
        0.7268041237113407,
        -0.2886597938144324,
        -1.8195876288659782,
        0.42647058823529393,
        -0.14705882352941213,
        -0.7941176470588243,
        -0.22204472843450485,
        -0.16613418530351431,
        0.4456869009584654,
    ],
    "Rfactor": 0.12540821749218442,
    "fitness": 224.29871241149283,
    "chi2_red": 28.037339051436604,
    "covariance": [[0.003918139665577524]],
    "fit_errors": [0.33144261985248585],
    "nparameters": 1,
}

_ANGLE_INCREMENT5_REFERENCES = {
    (False, False): {
        "weighted_residues": [
            9.791666666666671,
            -1.9444444444444402,
            -8.171296296296289,
            12.175324675324672,
            2.900432900432897,
            0.7756132756132704,
            -1.956168831168832,
            -1.7748917748917759,
            -1.4727633477633482,
        ],
        "weighted_residues2": [
            95.8767361111112,
            3.780864197530848,
            66.77008316186544,
            148.23853094956985,
            8.412511009913588,
            0.601575953307547,
            3.826596496036434,
            3.1502408125784784,
            2.1690318785151046,
        ],
        "residues": [
            0.7268041237113407,
            -0.2886597938144324,
            -1.8195876288659782,
            0.6563593932322052,
            0.3127187864644103,
            0.12543757292882063,
            -0.4218203033838974,
            -0.7654609101516923,
            -0.9527421236872815,
        ],
        "Rfactor": 0.15817348573575216,
        "fitness": 332.8261705704285,
        "scaled_errors": [
            0.07422680412371135,
            0.1484536082474227,
            0.2226804123711341,
            0.10781796966161027,
            0.21563593932322053,
            0.32345390898483084,
            0.10781796966161027,
            0.21563593932322053,
            0.32345390898483084,
        ],
        "log_prob": -158.8930549945955,
        "chi2_red": 41.60327132130356,
        "covariance": [[0.003105474966227227]],
        "fit_errors": [0.35944111840671117],
        "nparameters": 1,
    },
    (False, True): {
        "weighted_residues": [
            3.95681001721084,
            -0.5064800670980821,
            -3.534743294629085,
            6.986951803638445,
            7.779918531930959,
            -2.172487624896587,
            -3.6447549598168343,
            -0.8575372622484914,
            -0.946378129577784,
        ],
        "weighted_residues2": [
            15.656345512300046,
            0.2565220583676778,
            12.494410158925277,
            48.81749550636652,
            60.52713236348277,
            4.719702480328815,
            13.284238717109414,
            0.7353701561446379,
            0.8956315641431449,
        ],
        "residues": [
            0.2278989259663886,
            -0.07942766948096658,
            -0.8959459668551952,
            0.33584989408253163,
            1.0160850643985304,
            -0.25006232738410006,
            -0.6236448962933355,
            -0.3741112284119019,
            -0.5582297053586887,
        ],
        "Rfactor": 0.09435610305535341,
        "fitness": 157.38684851716832,
        "scaled_errors": [
            0.05759663086554629,
            0.15682289322076137,
            0.25346846777149357,
            0.090873159671776,
            0.24829235791747495,
            0.3062741836174822,
            0.08675780903984016,
            0.23930494561288085,
            0.33695273969576667,
        ],
        "log_prob": -70.9472482493825,
        "chi2_red": 26.231141419528054,
        "covariance": [
            [0.033723271892569386, -0.012719795766604647, -0.030032290819859638],
            [-0.012719795766604645, 0.006727675930074963, 0.013448686123837829],
            [-0.030032290819859635, 0.013448686123837829, 0.031465060957814026],
        ],
        "fit_errors": [
            0.9405317188394994,
            0.42008882244895696,
            0.9084957147716717,
        ],
        "nparameters": 3,
    },
    (True, False): {
        "weighted_residues": [
            9.791666666666671,
            -1.9444444444444402,
            -8.171296296296289,
            8.854961832061067,
            -1.5267175572519125,
            -5.496183206106875,
            -0.8687500000000002,
            -0.3249999999999999,
            0.5812499999999986,
        ],
        "weighted_residues2": [
            95.8767361111112,
            3.780864197530848,
            66.77008316186544,
            78.41034904725828,
            2.330866499621247,
            30.208029835091246,
            0.7547265625000005,
            0.10562499999999994,
            0.33785156249999837,
        ],
        "residues": [
            0.7268041237113407,
            -0.2886597938144324,
            -1.8195876288659782,
            0.42647058823529393,
            -0.14705882352941213,
            -0.7941176470588243,
            -0.22204472843450485,
            -0.16613418530351431,
            0.4456869009584654,
        ],
        "Rfactor": 0.12540821749218442,
        "fitness": 278.5751319774783,
        "scaled_errors": [
            0.07422680412371135,
            0.1484536082474227,
            0.2226804123711341,
            0.0963235294117647,
            0.1926470588235294,
            0.28897058823529415,
            0.12779552715654952,
            0.25559105431309903,
            0.3833865814696486,
        ],
        "log_prob": -131.93930216195574,
        "chi2_red": 34.82189149718479,
        "covariance": [[0.0026845497859506883]],
        "fit_errors": [0.3057467928943254],
        "nparameters": 1,
    },
    (True, True): {
        "weighted_residues": [
            3.95681001721084,
            -0.5064800670980821,
            -3.534743294629085,
            4.17260504177382,
            5.033171809019492,
            -6.9043326380949885,
            -2.425322668570648,
            0.321379995701755,
            0.8328225750077206,
        ],
        "weighted_residues2": [
            15.656345512300046,
            0.2565220583676778,
            12.494410158925277,
            17.4106328346363,
            25.33281845910854,
            47.669809177463705,
            5.88219004668265,
            0.10328510163726007,
            0.6935934414424905,
        ],
        "residues": [
            0.2278989259663886,
            -0.07942766948096658,
            -0.8959459668551952,
            0.17469613880547152,
            0.5757653261004405,
            -0.974254670684676,
            -0.48994288399810015,
            0.17907619877320524,
            0.6534145985187099,
        ],
        "Rfactor": 0.1027384789514959,
        "fitness": 125.49960679056394,
        "scaled_errors": [
            0.05759663086554629,
            0.15682289322076137,
            0.25346846777149357,
            0.08373480694027359,
            0.22878826630502203,
            0.2822154498493247,
            0.10100571160018999,
            0.2786050799182137,
            0.3922891970793734,
        ],
        "log_prob": -55.21436746612869,
        "chi2_red": 20.916601131760657,
        "covariance": [
            [0.03573138020249195, -0.013249923383573824, -0.03198400403470676],
            [-0.013249923383573824, 0.006922791246815929, 0.013894263605147884],
            [-0.03198400403470676, 0.013894263605147884, 0.03285443965538678],
        ],
        "fit_errors": [
            0.8645108607883628,
            0.3805276087040722,
            0.8289771525137394,
        ],
        "nparameters": 3,
    },
}


def _increment6_reference(
    weighted_residues,
    residues,
    flat_fcalc,
    rfactor,
    fitness,
    covariance,
    fit_errors,
    *,
    nparameters=1,
    log_prob=None,
):
    """Build a compact numerical characterization after the scale correction."""
    reference = {
        "weighted_residues": weighted_residues,
        "weighted_residues2": np.square(weighted_residues),
        "residues": residues,
        "flat_Fcalc": flat_fcalc,
        "Rfactor": rfactor,
        "fitness": fitness,
        "chi2_red": fitness / (9 - nparameters),
        "covariance": covariance,
        "fit_errors": fit_errors,
        "nparameters": nparameters,
    }
    if log_prob is not None:
        reference["likelihood_errors"] = [0.1, 0.2, 0.3] * 3
        reference["log_prob"] = log_prob
    return reference


_BASE_LEGACY_REFERENCE = _increment6_reference(
    [
        11.089494163424124,
        -0.21400778210116655,
        -5.71984435797665,
        6.658359573818854,
        -0.5502776507288288,
        -3.1365826091543325,
        -1.1830969490669867,
        -0.39895129677840246,
        0.9079581237025696,
    ],
    [
        1.1089494163424125,
        -0.04280155642023331,
        -1.7159533073929953,
        0.4708171206225682,
        -0.07782101167315147,
        -0.6653696498054469,
        -0.16731517509727634,
        -0.11284046692607008,
        0.3852140077821007,
    ],
    [
        1.8910505836575875,
        5.042801556420233,
        10.715953307392995,
        1.5291828793774318,
        4.0778210116731515,
        8.665369649805447,
        1.1673151750972763,
        3.11284046692607,
        6.614785992217899,
    ],
    0.11302575504910131,
    212.5972762645914,
    [[0.004225212366360498]],
    [0.3350874215648266],
)

_ANGLE_INCREMENT6_REFERENCES = {
    (False, False): _increment6_reference(
        [
            14.191217342968315,
            3.9216231239577537,
            0.13896609227348478,
            8.382434685936628,
            -2.1567537520844926,
            -6.3887344821196965,
            -2.904391328515843,
            -3.039188438021123,
            -3.2638502871965906,
        ],
        [
            1.4191217342968314,
            0.7843246247915507,
            0.04168982768204543,
            0.4191217342968314,
            -0.21567537520844926,
            -0.9583101723179546,
            -0.5808782657031686,
            -1.2156753752084493,
            -1.9583101723179546,
        ],
        [1.5808782657031686, 4.215675375208449, 8.958310172317955] * 3,
        0.18078826861483893,
        360.84668951886846,
        [[0.002912941697409698]],
        [0.3624785111340676],
        log_prob=-173.34580412901403,
    ),
    (False, True): _increment6_reference(
        [
            9.506296825308462,
            4.928640390376935,
            3.611102342924158,
            3.048178322842107,
            3.935753592035627,
            -8.7948647625582,
            -4.67615719147934,
            -1.8546716661215514,
            -2.4512354229742552,
        ],
        [
            0.7383976365153861,
            1.042361689173422,
            1.2343706415994742,
            0.13131183968621873,
            0.46325395961686056,
            -1.2769310383261878,
            -0.7692819735730783,
            -0.8416016717833195,
            -1.5661807943811779,
        ],
        [1.5808782657031686, 4.215675375208449, 8.958310172317955] * 3,
        0.18440897989307575,
        261.14723734373257,
        [
            [0.038710415383091536, -0.01332873676059801, -0.03404602704901221],
            [-0.01332873676059801, 0.006745271067149528, 0.013600952381716463],
            [-0.034046027049012204, 0.013600952381716462, 0.03385340332920189],
        ],
        [1.2980188772203765, 0.5418346771726456, 1.213859049622607],
        nparameters=3,
        log_prob=-123.4960780414461,
    ),
    (True, False): _increment6_reference(
        [
            11.089494163424124,
            -0.21400778210116655,
            -5.71984435797665,
            9.416342412451364,
            -0.7782101167315147,
            -4.435797665369646,
            -0.8365758754863817,
            -0.2821011673151752,
            0.6420233463035011,
        ],
        _BASE_LEGACY_REFERENCE["residues"],
        _BASE_LEGACY_REFERENCE["flat_Fcalc"],
        0.11302575504910131,
        265.8803501945525,
        [[0.002860868542478485]],
        [0.30835221945341146],
        log_prob=-125.86263446685605,
    ),
    (True, True): _increment6_reference(
        [
            4.258356416588224,
            -0.21114807320227363,
            -3.146454413688498,
            4.70061591801101,
            5.54849997099951,
            -6.016572229359093,
            -2.2039595208052933,
            0.5353884626319467,
            1.1558003010370153,
        ],
        [
            0.3307660565655486,
            -0.04465585735947819,
            -1.0755416447829376,
            0.20249685499922854,
            0.6530806671182203,
            -0.8735492848857973,
            -0.36257684684545377,
            0.24294533282367903,
            0.7384815904086945,
        ],
        _BASE_LEGACY_REFERENCE["flat_Fcalc"],
        0.10346174713399309,
        123.639094296015,
        [
            [0.036781230132263185, -0.01364452181463681, -0.03293541971137809],
            [-0.01364452181463681, 0.007128572724362648, 0.01431332513038369],
            [-0.03293541971137808, 0.01431332513038369, 0.033844956383305896],
        ],
        [0.8705934738868573, 0.3832688950035545, 0.8351207251769125],
        nparameters=3,
        log_prob=-54.74200651758731,
    ),
}


def _increment7_reference(
    settings,
    weighted_residues,
    fitness,
    covariance,
    fit_errors,
    log_prob,
    *,
    residues=None,
    rfactor=None,
    pvalue=0.0,
):
    """Update the Increment 6 snapshot for linear angle-optimizer weights."""
    reference = _ANGLE_INCREMENT6_REFERENCES[settings].copy()
    reference.update(
        {
            "weighted_residues": weighted_residues,
            "weighted_residues2": np.square(weighted_residues),
            "fitness": fitness,
            "chi2_red": fitness / (9 - reference["nparameters"]),
            "covariance": covariance,
            "fit_errors": fit_errors,
            "log_prob": log_prob,
            "pvalue": pvalue,
        }
    )
    if residues is not None:
        reference["residues"] = residues
    if rfactor is not None:
        reference["Rfactor"] = rfactor
    return reference


_ANGLE_LEGACY_REFERENCES = {
    (False, False): _increment7_reference(
        (False, False),
        [
            14.191217342968315,
            3.9216231239577537,
            0.13896609227348478,
            5.927276409279118,
            -1.5250552034484748,
            -4.517517475507163,
            -4.107429607225916,
            -4.298061507656975,
            -4.61578134170874,
        ],
        331.30520041998636,
        [[0.004369412546133979]],
        [0.4253835180332156],
        -158.57505957957298,
    ),
    (False, True): _increment7_reference(
        (False, True),
        [
            9.32619684343576,
            4.752251990923697,
            3.379195012698496,
            1.9257652516342938,
            2.558891224669818,
            -6.6049793416623865,
            -6.733342134496628,
            -2.7391635521660676,
            -3.642030487150156,
        ],
        240.9683530245108,
        [
            [0.03412554948754149, -0.013677283319011595, -0.029881054116630553],
            [-0.013677283319011595, 0.008771974751680327, 0.015515018585677717],
            [-0.029881054116630553, 0.015515018585677717, 0.033308496359875685],
        ],
        [1.1706961928748005, 0.5935442006406657, 1.1565965523582658],
        -113.40663588183521,
        residues=[
            0.7244084456248574,
            1.0050571801320123,
            1.1550985598864791,
            0.11732264879568999,
            0.42594945057545086,
            -1.3562031200391829,
            -0.783271164463607,
            -0.8789061808247292,
            -1.645452876094173,
        ],
        rfactor=0.1850488189800105,
    ),
    (True, False): _increment7_reference(
        (True, False),
        _BASE_LEGACY_REFERENCE["weighted_residues"],
        212.5972762645914,
        [[0.004225212366316432]],
        [0.33508742156307925],
        -99.2210975018755,
    ),
    (True, True): _increment7_reference(
        (True, True),
        [
            4.258356416588224,
            -0.21114807320227363,
            -3.146454413688498,
            3.323837391379014,
            3.923381954907116,
            -4.254359022878479,
            -3.1168694452441534,
            0.7571536249921801,
            1.6345484611214534,
        ],
        85.57865534319242,
        [
            [0.03951006763610368, -0.01674102685774627, -0.0361966169600124],
            [-0.01674102685774627, 0.010369758737461677, 0.019408419295404156],
            [-0.0361966169600124, 0.019408419295404156, 0.04144650043841291],
        ],
        [0.7506906221180976, 0.3845841929897126, 0.7688666742030562],
        -35.71178704117603,
        pvalue=2.220446049250313e-16,
    ),
}


def _legacy_characterization_optimizer(
    optimizer_type, scaleindividual=None, use_angle_correction=None
):
    """Return a prepared optimizer with identifiable legacy angle inputs."""
    optimizer = optimizer_type(FitCrystal(), _fixture_ctrs())
    if scaleindividual is not None:
        optimizer.scaleindividual = scaleindividual
    if use_angle_correction is not None:
        optimizer.useAnglecorr = use_angle_correction
    omega_values = (
        (-2.0, 0.0, 2.0),
        (-1.0, 1.0, 3.0),
        (-2.5, 0.5, 2.5),
    )
    for ctr, omega in zip(optimizer.CTRs, omega_values):
        ctr.angles["omega"] = omega
    if isinstance(optimizer, CTRopt.CTROptAngleCorrection):
        optimizer.prepareFit(start=[0.3, 0.4])
    else:
        optimizer.prepareFit()
    return optimizer


class TestLegacyOptimizerCharacterization(unittest.TestCase):
    """Pin intentional numerical changes throughout formula centralization."""

    def _assert_legacy_outputs(self, optimizer, expected):
        """Assert the legacy objective, output, and reporting definitions."""
        parameters = optimizer.get_parameters()

        # Legacy scaling and weight powers; increments 6 and 7 update these.
        np.testing.assert_allclose(
            optimizer.weighted_residues(parameters),
            expected["weighted_residues"],
        )
        np.testing.assert_allclose(
            optimizer.weighted_residues2(parameters),
            expected["weighted_residues2"],
        )
        np.testing.assert_allclose(optimizer.fitness(parameters), [expected["fitness"]])

        # These observation-scaled outputs migrate in increments 6 and 8.
        np.testing.assert_allclose(optimizer.residues(parameters), expected["residues"])
        np.testing.assert_allclose(
            optimizer.flat_Fcalc(parameters), expected["flat_Fcalc"]
        )
        np.testing.assert_allclose(optimizer.Rfactor(parameters), expected["Rfactor"])

        if "likelihood_errors" in expected:
            residuals, errors = optimizer.weighted_residues_errors(parameters)
            np.testing.assert_allclose(residuals, expected["weighted_residues"])
            np.testing.assert_allclose(errors, expected["likelihood_errors"])
            np.testing.assert_allclose(
                optimizer.log_prob(parameters), expected["log_prob"]
            )

        statistics = optimizer.statistics(parameters)
        self.assertEqual(
            set(statistics),
            {
                "Chisqr",
                "nodatapoints",
                "Chisqr_red",
                "noparameters",
                "pvalue",
                "Rfactor",
                "covariance",
            },
        )
        np.testing.assert_allclose(statistics["Chisqr"], expected["fitness"])
        self.assertEqual(statistics["nodatapoints"], 9)
        np.testing.assert_allclose(statistics["Chisqr_red"], expected["chi2_red"])
        self.assertEqual(statistics["noparameters"], parameters.size)
        self.assertEqual(parameters.size, expected["nparameters"])
        np.testing.assert_allclose(
            statistics["pvalue"], expected.get("pvalue", 0.0)
        )
        np.testing.assert_allclose(statistics["Rfactor"], expected["Rfactor"])
        # Increment 8 replaces the unscaled covariance/scaled-error pairing.
        np.testing.assert_allclose(
            statistics["covariance"], expected["covariance"], rtol=1e-6
        )
        np.testing.assert_allclose(optimizer.errors, expected["fit_errors"], rtol=1e-6)

    def test_base_optimizer_legacy_outputs(self):
        """Pin every public numerical output of the base optimizer."""
        optimizer = _legacy_characterization_optimizer(CTRopt.CTROptimizer)
        self._assert_legacy_outputs(optimizer, _BASE_LEGACY_REFERENCE)

    def test_angle_optimizer_legacy_outputs(self):
        """Pin shared/individual scaling with angle correction on and off."""
        for settings, expected in _ANGLE_LEGACY_REFERENCES.items():
            scaleindividual, use_angle_correction = settings
            with self.subTest(
                scaleindividual=scaleindividual,
                use_angle_correction=use_angle_correction,
            ):
                optimizer = _legacy_characterization_optimizer(
                    CTRopt.CTROptAngleCorrection,
                    scaleindividual,
                    use_angle_correction,
                )
                self._assert_legacy_outputs(optimizer, expected)

    def test_rod_weight_power_is_linear_in_both_classes(self):
        """Doubling a rod weight doubles its chi-square influence."""

        def rod_chisqr(optimizer_type, weight):
            ctrs = _fixture_ctrs()
            ctrs[1].weight = weight
            optimizer = optimizer_type(FitCrystal(), ctrs)
            if isinstance(optimizer, CTRopt.CTROptAngleCorrection):
                optimizer.scaleindividual = True
            optimizer.prepareFit()
            parameters = optimizer.get_parameters()
            return np.sum(optimizer.weighted_residues2(parameters)[3:6])

        base_unit = rod_chisqr(CTRopt.CTROptimizer, 1.0)
        base_doubled = rod_chisqr(CTRopt.CTROptimizer, 2.0)
        angle_unit = rod_chisqr(CTRopt.CTROptAngleCorrection, 1.0)
        angle_doubled = rod_chisqr(CTRopt.CTROptAngleCorrection, 2.0)

        np.testing.assert_allclose(base_doubled / base_unit, 2.0)
        np.testing.assert_allclose(angle_doubled / angle_unit, 2.0)

    def test_increment6_angle_weight_response_was_fourfold(self):
        """Keep the superseded quadratic response explicit in characterization."""
        legacy_middle = np.sum(
            _ANGLE_INCREMENT6_REFERENCES[(True, False)][
                "weighted_residues2"
            ][3:6]
        )
        linear_middle = np.sum(
            _ANGLE_LEGACY_REFERENCES[(True, False)][
                "weighted_residues2"
            ][3:6]
        )
        unit_weight_value = linear_middle / 2.0

        np.testing.assert_allclose(legacy_middle / unit_weight_value, 4.0)

    def test_scale_estimator_multiplies_the_prediction(self):
        """The analytical scale produces the accepted effective prediction."""
        crystal = _ConstantFitCrystal()
        ctrs = _unit_model_ctrs((((1.0, 0.0), [1.0, 2.0]),))
        optimizer = CTRopt.CTROptimizer(crystal, ctrs)
        optimizer.prepareFit()

        np.testing.assert_allclose(
            optimizer.flat_Fcalc(optimizer.get_parameters()), [3.0 / 2.0] * 2
        )

    def test_shared_scale_rfactor_keeps_legacy_index_misalignment(self):
        """Pin D3 until the subclass R-factor override is removed."""
        crystal = _ConstantFitCrystal()
        ctrs = _unit_model_ctrs(
            (
                ((0.0, 0.0), [3.0, 3.0]),
                ((1.0, 0.0), [1.0, 2.0]),
                ((2.0, 0.0), [1.0, 4.0]),
            )
        )
        optimizer = CTRopt.CTROptAngleCorrection(crystal, ctrs)
        optimizer.scaleindividual = False
        optimizer.set_scale_policy((0.0, 0.0), "scaled")
        optimizer.prepareFit()

        # D3 repeats the specular rod and omits the last non-specular rod.
        # Increment 8 deletes this override rather than preserving the defect.
        self.assertAlmostEqual(
            optimizer.Rfactor(optimizer.get_parameters()),
            0.21428571428571427,
        )

    def test_evaluate_statistics_remains_deprecated(self):
        """The excluded legacy statistics path still emits its warning."""
        optimizer = CTRopt.CTROptimizer(FitCrystal(), _fixture_ctrs())
        optimizer.prepareFit()
        with self.assertWarnsRegex(
            DeprecationWarning, "evaluateStatistics is deprecated"
        ):
            optimizer.evaluateStatistics(optimizer.get_parameters())


class TestLegacyCalculationAdapters(unittest.TestCase):
    """Increment 5: share inputs while legacy public outputs still differ."""

    def test_shared_inputs_apply_the_angle_hook_once(self):
        """One record owns each model, observation, and uncertainty array."""
        ctrs = _unit_model_ctrs((((1.0, 0.0), [1.0, 2.0]),))
        ctrs[0].angles = _angles([0.0, 0.0])
        ctrs[0].angles["omega"] = [-1.0, 1.0]
        optimizer = CTRopt.CTROptAngleCorrection(_ConstantFitCrystal(), ctrs)
        optimizer.useAnglecorr = True
        optimizer.prepareFit(start=[0.3, 0.4])

        (calculation,) = optimizer._calculation_inputs()
        correction = optimizer.get_anglecorrection(
            optimizer.CTRs[0].angles["omega"]
        )

        np.testing.assert_allclose(calculation.prediction, [1.0, 1.0])
        np.testing.assert_allclose(calculation.angle_correction, correction)
        np.testing.assert_allclose(
            calculation.observation, np.array([1.0, 2.0]) * correction
        )
        np.testing.assert_allclose(calculation.uncertainty, correction)

    def test_adapters_preserve_incompatible_legacy_predictions(self):
        """Individual flattened and shared residual scales stay distinct."""
        ctrs = _unit_model_ctrs(
            (
                ((1.0, 0.0), [1.0, 2.0]),
                ((2.0, 0.0), [1.0, 4.0]),
            )
        )
        optimizer = CTRopt.CTROptAngleCorrection(_ConstantFitCrystal(), ctrs)
        optimizer.scaleindividual = False
        optimizer.prepareFit()
        parameters = optimizer.get_parameters()

        # The inherited adapter now follows the explicit global policy.
        np.testing.assert_allclose(
            optimizer.flat_Fcalc(parameters),
            [2.0, 2.0, 2.0, 2.0],
        )

        # The subclass adapter uses the same prediction-side global scale.
        calculations = optimizer._legacy_angle_calculations(
            corrected_scale_errors=False
        )
        np.testing.assert_allclose(
            np.concatenate(
                [
                    calculation.scaled_prediction
                    for calculation in calculations
                ]
            ),
            [2.0, 2.0, 2.0, 2.0],
        )

    def test_apply_corrections_consumes_the_shared_input_records(self):
        """The temporary adapter preserves individual correction export."""
        optimizer = _legacy_characterization_optimizer(
            CTRopt.CTROptAngleCorrection,
            scaleindividual=True,
            use_angle_correction=True,
        )
        calculations = optimizer._scaled_calculations()
        expected_values = []
        expected_errors = []
        for calculation in calculations:
            expected_values.append(
                calculation.values.observation / calculation.scale
            )
            expected_errors.append(
                calculation.values.uncertainty / calculation.scale
            )

        optimizer.applyCorrections()

        for ctr, values, errors in zip(
            optimizer.CTRs, expected_values, expected_errors
        ):
            np.testing.assert_allclose(ctr.sfI, values)
            np.testing.assert_allclose(ctr.err, errors)
        self.assertEqual(optimizer.amp, 0.0)


class TestScaleEstimation(unittest.TestCase):
    """Increment 6: analytical scales multiply calculated predictions."""

    def test_prepare_and_direct_evaluation_describe_zero_prediction_norm(self):
        """A transient zero model fails clearly and a later trial can recover."""
        ctrs = _unit_model_ctrs((((1.0, 0.0), [1.0, 2.0]),))
        optimizer = CTRopt.CTROptimizer(_VanishingFitCrystal((0.0,)), ctrs)

        with self.assertRaisesRegex(
            ValueError,
            r"rod .*current parameters.*weighted prediction norm is zero",
        ):
            optimizer.prepareFit()
        with self.assertRaisesRegex(
            ValueError,
            r"rod .*current parameters.*weighted prediction norm is zero",
        ):
            optimizer.weighted_residues([0.0])

        self.assertTrue(np.isinf(optimizer.fitness([0.0])[0]))
        self.assertTrue(np.isfinite(optimizer.fitness([1.0])[0]))

    def test_nonpositive_scale_fails_but_fixed_policy_skips_estimation(self):
        """Zero observations violate fitted scales but remain valid when fixed."""
        ctrs = _unit_model_ctrs((((1.0, 0.0), [0.0, 0.0]),))
        fitted = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)
        with self.assertRaisesRegex(
            ValueError,
            r"rod .*current parameters.*analytical scale is nonpositive",
        ):
            fitted.prepareFit()

        fixed = CTRopt.CTROptimizer(
            _ConstantFitCrystal(), ctrs, scale_policy={"F": "fixed"}
        )
        fixed.prepareFit()
        np.testing.assert_allclose(fixed.flat_Fcalc(), [1.0, 1.0])

    def test_global_scale_failure_names_the_group(self):
        """Shared scale errors identify the affected global group."""
        ctrs = _unit_model_ctrs(
            (((1.0, 0.0), [0.0, 0.0]), ((2.0, 0.0), [0.0, 0.0]))
        )
        optimizer = CTRopt.CTROptimizer(
            _ConstantFitCrystal(), ctrs, scale_policy={"F": "global"}
        )
        with self.assertRaisesRegex(
            ValueError,
            r"global group .*current parameters.*analytical scale is nonpositive",
        ):
            optimizer.prepareFit()

    def test_scalar_angle_objective_rejects_only_scale_failures(self):
        """Likelihood trials map scale failures to minus infinity."""
        ctrs = _unit_model_ctrs((((1.0, 0.0), [0.0, 0.0]),))
        optimizer = CTRopt.CTROptAngleCorrection(_ConstantFitCrystal(), ctrs)
        self.assertEqual(optimizer.log_prob([1.0]), -np.inf)

        with mock.patch.object(
            optimizer,
            "weighted_residues_errors",
            side_effect=ValueError("unrelated failure"),
        ):
            with self.assertRaisesRegex(ValueError, "unrelated failure"):
                optimizer.log_prob([1.0])

        base = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)
        with mock.patch.object(
            base, "weighted_residues2", side_effect=ValueError("input failure")
        ):
            with self.assertRaisesRegex(ValueError, "input failure"):
                base.fitness([1.0])

    def test_likelihood_uses_supplied_uncertainties_unchanged(self):
        """Fitted scale and empirical correction do not inflate log sigma."""
        ctrs = _unit_model_ctrs((((1.0, 0.0), [1.0, 2.0]),))
        ctrs[0].angles = _angles([0.0, 0.0])
        ctrs[0].angles["omega"] = [-1.0, 1.0]
        optimizer = CTRopt.CTROptAngleCorrection(_ConstantFitCrystal(), ctrs)
        optimizer.useAnglecorr = True
        optimizer.prepareFit(start=[0.3, 0.4])

        _, errors = optimizer.weighted_residues_errors(
            optimizer.get_parameters()
        )
        np.testing.assert_array_equal(errors, optimizer.CTRs[0].err)


class TestScalePolicy(unittest.TestCase):
    """Increment 6: explicit policies are the sole scale-group authority."""

    def test_defaults_overrides_partial_updates_and_restoration(self):
        """Configuration round-trips and partial updates preserve overrides."""
        ctrs = _unit_model_ctrs(
            (((1.0, 0.0), [1.0, 2.0]), ((2.0, 0.0), [1.0, 4.0]))
        )
        optimizer = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)
        first_id = optimizer.CTRs[0].ctr_id
        second_id = optimizer.CTRs[1].ctr_id

        self.assertEqual(optimizer.get_scale_policy(), {"F": "scaled", "R": "fixed"})
        optimizer.set_scale_policies(
            {"F": "global", first_id: "fixed", second_id: "scaled"}
        )
        optimizer.set_scale_policies({"R": "scaled"})
        self.assertEqual(optimizer.get_scale_policy(first_id), "fixed")
        self.assertEqual(optimizer.get_scale_policy((2.0, 0.0)), "scaled")

        saved = optimizer.get_scale_policy()
        restored = CTRopt.CTROptimizer(
            _ConstantFitCrystal(), ctrs, scale_policy=saved
        )
        self.assertEqual(restored.get_scale_policy(), saved)

    def test_ambiguous_unknown_and_invalid_selectors_are_rejected(self):
        """Rod shorthand never silently resolves the first matching dataset."""
        ctrs = _unit_model_ctrs(
            (((1.0, 0.0), [1.0]), ((1.0, 0.0), [2.0]))
        )
        optimizer = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)

        with self.assertRaisesRegex(ValueError, "Ambiguous CTR shorthand"):
            optimizer.set_scale_policy((1.0, 0.0), "fixed")
        optimizer.set_scale_policy(optimizer.CTRs[0].ctr_id, "fixed")
        with self.assertRaisesRegex(ValueError, "Unknown CTR identifier"):
            optimizer.get_scale_policy((9.0, 9.0))
        with self.assertRaisesRegex(ValueError, "scale quantity must"):
            optimizer.set_scale_policy_default("structure_factor", "fixed")
        with self.assertRaisesRegex(ValueError, "scale policy must"):
            optimizer.set_scale_policy_default("F", "fitted")

    def test_global_group_rejects_actual_mixed_quantities_atomically(self):
        """There is one global scale group, never implicit F/R subgroups."""
        ctrs = _unit_model_ctrs(
            (((1.0, 0.0), [1.0]), ((2.0, 0.0), [2.0]))
        )
        ctrs[1].reduction = CTRplotutil.MeasurementReduction("reflectivity")
        optimizer = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)
        optimizer.set_scale_policy_default("F", "global")

        with self.assertRaisesRegex(ValueError, "cannot mix"):
            optimizer.set_scale_policy_default("R", "global")
        self.assertEqual(optimizer.get_scale_policy(), {"F": "global", "R": "fixed"})
        with self.assertRaisesRegex(ValueError, "cannot mix"):
            optimizer.scaleindividual = False

    def test_scaleindividual_is_write_only_and_preserves_overrides(self):
        """Compatibility assignments update defaults without shadow state."""
        ctrs = _unit_model_ctrs(
            (((1.0, 0.0), [1.0]), ((2.0, 0.0), [2.0]))
        )
        optimizer = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)
        first_id = optimizer.CTRs[0].ctr_id
        optimizer.set_scale_policy(first_id, "fixed")

        optimizer.scaleindividual = False
        self.assertEqual(
            optimizer.get_scale_policy(),
            {"F": "global", "R": "global", first_id: "fixed"},
        )
        optimizer.scaleindividual = True
        self.assertEqual(optimizer.get_scale_policy(first_id), "fixed")
        with self.assertRaisesRegex(AttributeError, "get_scale_policy"):
            _ = optimizer.scaleindividual

    def test_policy_changes_take_effect_after_preparation(self):
        """Evaluation resolves current policies instead of cached grouping."""
        ctrs = _unit_model_ctrs((((1.0, 0.0), [1.0, 2.0]),))
        optimizer = CTRopt.CTROptimizer(_ConstantFitCrystal(), ctrs)
        optimizer.prepareFit()
        np.testing.assert_allclose(optimizer.flat_Fcalc(), [1.5, 1.5])

        optimizer.set_scale_policy_default("F", "fixed")
        np.testing.assert_allclose(optimizer.flat_Fcalc(), [1.0, 1.0])
        optimizer.set_scale_policy((1.0, 0.0), "global")
        np.testing.assert_allclose(optimizer.flat_Fcalc(), [1.5, 1.5])

    def test_angle_adapters_resolve_fixed_scaled_and_global_policies(self):
        """Retained subclass methods no longer read scaleindividual."""
        optimizer = _legacy_characterization_optimizer(
            CTRopt.CTROptAngleCorrection,
            scaleindividual=True,
            use_angle_correction=True,
        )
        optimizer.set_scale_policy(optimizer.CTRs[0].ctr_id, "fixed")
        optimizer.set_scale_policy(optimizer.CTRs[1].ctr_id, "global")
        optimizer.set_scale_policy(optimizer.CTRs[2].ctr_id, "global")
        parameters = optimizer.get_parameters()

        self.assertEqual(optimizer.residues(parameters).size, 9)
        self.assertEqual(optimizer.weighted_residues(parameters).size, 9)
        self.assertEqual(optimizer.weighted_residues_errors(parameters)[0].size, 9)
        self.assertTrue(np.isfinite(optimizer.Rfactor(parameters)))
        optimizer.applyCorrections()
        self.assertEqual(optimizer.amp, 0.0)


class TestOptimizerClassSplit(unittest.TestCase):
    """Increment 3: generic machinery lives on ``CTROptimizer``."""

    def test_base_optimizer_runs_a_callback_fit_end_to_end(self):
        """Callbacks, statistics, and trace export work without the subclass."""
        az = pytest.importorskip("arviz")
        optimizer = CTRopt.CTROptimizer(FitCrystal(), _fixture_ctrs())
        callback, calls = _register_callback(optimizer, "cb", 1)
        optimizer.prepareFit()
        parameters = optimizer.get_parameters()

        self.assertEqual(optimizer.callback_names, ["cb"])
        self.assertEqual(optimizer.fitparnames, ["cb", "xtal_0"])
        self.assertEqual(parameters.size, 2)
        self.assertEqual(len(optimizer.fitness(parameters)), 1)
        stat = optimizer.statistics(parameters)
        self.assertEqual(stat["nodatapoints"], 9)
        self.assertTrue(np.all(np.isfinite(callback.errors)))
        self.assertTrue(np.all(np.isfinite(optimizer.xtal.errors)))
        self.assertTrue(calls)

        population = _FakePopulation(
            np.vstack((parameters, parameters + 0.01)), [1.0, 2.0]
        )
        archipelago = _FakeArchipelago([_FakeIsland(population)])
        trace = object()
        with mock.patch.object(az, "from_dict", return_value=trace):
            self.assertIs(optimizer.set_archi_result(archipelago), trace)

        optimizer.unregister_fit_callback("cb")
        self.assertEqual(optimizer.callback_names, [])
        with self.assertRaisesRegex(ValueError, "not a registered callback"):
            optimizer.unregister_fit_callback("cb")

    def test_base_optimizer_exposes_constraint_methods(self):
        """Disabled constraints return an empty float vector and flat fitness."""
        optimizer = CTRopt.CTROptimizer(FitCrystal(), _fixture_ctrs())
        optimizer.prepareFit()
        parameters = optimizer.get_parameters()

        constraints = optimizer.get_inequalconstraints()
        self.assertEqual(constraints.dtype, np.dtype(np.float64))
        self.assertEqual(constraints.size, 0)
        self.assertEqual(optimizer.get_nic(), 0)
        np.testing.assert_allclose(
            optimizer.fitness(parameters),
            [np.sum(optimizer.weighted_residues2(parameters))],
        )

    def test_constraints_extend_fitness_on_the_base_class(self):
        """Enabled displacement constraints follow the scalar objective."""
        optimizer = CTRopt.CTROptimizer(
            _ConstraintFitCrystal(), _fixture_ctrs()
        )
        optimizer.dw_zconstraints = True
        optimizer.prepareFit()
        parameters = optimizer.get_parameters()
        constraints = optimizer.get_inequalconstraints()

        self.assertEqual(optimizer.get_nic(), constraints.size)
        self.assertGreater(constraints.size, 0)
        np.testing.assert_allclose(
            optimizer.fitness(parameters),
            np.concatenate(
                ([np.sum(optimizer.weighted_residues2(parameters))], constraints)
            ),
        )

    def test_subclass_bounds_and_parameters_are_unchanged_by_the_split(self):
        """Every optional block keeps its legacy position and bounds."""
        for fit_resolution in (False, True):
            for with_callback in (False, True):
                for use_angle_correction in (False, True):
                    with self.subTest(
                        fit_resolution=fit_resolution,
                        with_callback=with_callback,
                        use_angle_correction=use_angle_correction,
                    ):
                        optimizer = CTRopt.CTROptAngleCorrection(
                            FitCrystal(
                                (1.0, 2.0), names=["xtal_a", "xtal_b"]
                            ),
                            _fixture_ctrs(),
                        )
                        optimizer.useAnglecorr = use_angle_correction
                        parameter_blocks = []
                        lower_blocks = []
                        upper_blocks = []
                        if fit_resolution:
                            optimizer.fit_resolution(
                                CTRresolution.BoxResolution(0.11, 0.12, 0.13),
                                lower_bounds=[0.01, 0.02, 0.03],
                                higher_bounds=[0.5, 0.6, 0.7],
                            )
                            parameter_blocks.append([0.11, 0.12, 0.13])
                            lower_blocks.append([0.01, 0.02, 0.03])
                            upper_blocks.append([0.5, 0.6, 0.7])
                        if with_callback:
                            _register_callback(optimizer, "cb", 2)
                            parameter_blocks.append([1.0, 0.5])
                            lower_blocks.append([0.1, 0.1])
                            upper_blocks.append([5.0, 5.0])
                        if use_angle_correction:
                            parameter_blocks.append([0.25, 0.4])
                            lower_blocks.append([-1.0, 0.0])
                            upper_blocks.append([1.0, 0.8])
                        parameter_blocks.append([1.0, 2.0])
                        lower_blocks.append([0.1, 0.1])
                        upper_blocks.append([10.0, 10.0])

                        optimizer.prepareFit(
                            phaselim=[-1.0, 1.0],
                            amplim=[0.0, 0.8],
                            start=[0.25, 0.4],
                        )

                        np.testing.assert_allclose(
                            optimizer.get_parameters(),
                            np.concatenate(parameter_blocks),
                        )
                        lower, upper = optimizer.get_bounds()
                        np.testing.assert_allclose(
                            lower, np.concatenate(lower_blocks)
                        )
                        np.testing.assert_allclose(
                            upper, np.concatenate(upper_blocks)
                        )
                        for ctr in optimizer.CTRs:
                            np.testing.assert_allclose(
                                ctr.invrelerrsqrd_weight,
                                ctr.weight * ctr.err**-2,
                            )

    def test_deleted_members_are_gone(self):
        """The obsolete correction and commented plotting members stay absent."""
        optimizer = CTRopt.CTROptAngleCorrection(
            FitCrystal(), _fixture_ctrs()
        )
        self.assertFalse(hasattr(optimizer, "get_anglecorrection_"))
        self.assertFalse(hasattr(optimizer, "defaultCTRplotsettings"))
        self.assertFalse(hasattr(optimizer, "plotParametersetCTR"))


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


class TestParameterNameLayout(unittest.TestCase):
    """Increment 2: fit trace names follow the optimizer vector layout."""

    @staticmethod
    def _optimizer(callback_count, fit_resolution, use_angle_correction):
        optimizer = CTRopt.CTROptAngleCorrection(
            FitCrystal((1.0, 2.0), names=["xtal_a", "xtal_b"]),
            _fixture_ctrs(),
        )
        optimizer.useAnglecorr = use_angle_correction
        if fit_resolution:
            optimizer.fit_resolution(
                CTRresolution.BoxResolution(0.1, 0.0, 0.0),
                lower_bounds=[0.0, 0.0, 0.0],
                higher_bounds=[1.0, 1.0, 1.0],
            )
        if callback_count >= 1:
            _register_callback(optimizer, "first", 2)
        if callback_count == 2:
            _register_callback(
                optimizer, "second", 1, parnames=["second_value"]
            )
        return optimizer

    def test_fitparnames_match_parameter_vector_position_by_position(self):
        """Every optional parameter block contributes names in vector order."""
        for fit_resolution in (False, True):
            for callback_count in (0, 1, 2):
                for use_angle_correction in (False, True):
                    with self.subTest(
                        fit_resolution=fit_resolution,
                        callback_count=callback_count,
                        use_angle_correction=use_angle_correction,
                    ):
                        optimizer = self._optimizer(
                            callback_count,
                            fit_resolution,
                            use_angle_correction,
                        )
                        optimizer.prepareFit()

                        expected = []
                        if fit_resolution:
                            expected += [
                                "resolution_delta_l_0",
                                "resolution_delta_l_1",
                                "resolution_delta_l_2",
                            ]
                        if callback_count == 1:
                            expected += ["first_0", "first_1"]
                        elif callback_count == 2:
                            expected += ["second_value", "first_0", "first_1"]
                        if use_angle_correction:
                            expected += [
                                "anglecorrection_phase",
                                "anglecorrection_amplitude",
                            ]
                        expected += ["xtal_a", "xtal_b"]

                        self.assertEqual(optimizer.fitparnames, expected)
                        self.assertEqual(
                            len(optimizer.fitparnames),
                            optimizer.get_parameters().size,
                        )

    def test_callback_parameter_names_validate_length(self):
        """An explicit callback name list must cover every callback value."""
        with self.assertRaisesRegex(ValueError, "Number of parameter names"):
            CTRopt.FitCallback(
                lambda xtal, x: None,
                [0.1, 0.1],
                [5.0, 5.0],
                [1.0, 2.0],
                name="callback",
                parnames=["only_one"],
            )

    def test_duplicate_parameter_names_are_rejected(self):
        """Duplicate callback/crystal names cannot silently drop a trace."""
        optimizer = CTRopt.CTROptAngleCorrection(
            FitCrystal(names=["duplicate"]), _fixture_ctrs()
        )
        _register_callback(optimizer, "duplicate", 1)

        with self.assertRaisesRegex(ValueError, "Duplicate fit parameter names"):
            optimizer.prepareFit()

    def test_archi_trace_columns_follow_the_actual_layout(self):
        """Each exported trace column carries its named parameter values."""
        az = pytest.importorskip("arviz")
        optimizer = self._optimizer(2, True, True)
        optimizer.prepareFit()
        npars = len(optimizer.fitparnames)
        first = np.arange(3 * npars, dtype=np.float64).reshape(3, npars)
        second = first + 1000.0
        archi = _FakeArchipelago(
            [
                _FakeIsland(_FakePopulation(first, [3.0, 1.0, 2.0])),
                _FakeIsland(_FakePopulation(second, [6.0, 4.0, 5.0])),
            ]
        )
        optimizer.statistics = mock.Mock(
            return_value={"nodatapoints": 9, "covariance": np.identity(npars)}
        )
        captured = {}
        trace = object()

        def capture_from_dict(params, attrs):
            captured["params"] = params
            captured["attrs"] = attrs
            return trace

        with mock.patch.object(az, "from_dict", side_effect=capture_from_dict):
            actual = optimizer.set_archi_result(archi)

        self.assertIs(actual, trace)
        self.assertEqual(
            list(captured["params"]), optimizer.fitparnames + ["chisqr"]
        )
        for column, name in enumerate(optimizer.fitparnames):
            np.testing.assert_allclose(
                captured["params"][name],
                np.stack((first[:, column], second[:, column])),
            )
        np.testing.assert_allclose(
            captured["params"]["chisqr"],
            [[3.0, 1.0, 2.0], [6.0, 4.0, 5.0]],
        )
        self.assertEqual(captured["attrs"], {"nodatapoints": 9})

    def test_archi_result_rejects_a_name_vector_length_mismatch(self):
        """Trace export fails explicitly when its positional schema is stale."""
        optimizer = self._optimizer(1, False, False)
        optimizer.prepareFit()
        parameters = optimizer.get_parameters()[None, :]
        archi = _FakeArchipelago(
            [_FakeIsland(_FakePopulation(parameters, [1.0]))]
        )
        optimizer.fitparnames.pop()
        optimizer.statistics = mock.Mock()

        with self.assertRaisesRegex(ValueError, "name count.*vector length"):
            optimizer.set_archi_result(archi)
        optimizer.statistics.assert_not_called()


class TestStoredQuantityGuards(unittest.TestCase):
    """Kinematical optimizers must not reinterpret reflectivity as F."""

    def test_both_kinematical_optimizers_reject_reflectivity(self):
        """Preparation names the unsupported reflectivity dataset."""
        ctrs = _fixture_ctrs()
        ctrs[1].reduction = CTRplotutil.MeasurementReduction(
            "reflectivity", CTRplotutil.PolarizationReduction(1.0, "s")
        )
        for optimizer_type in (
            CTRopt.CTROptimizer,
            CTRopt.CTROptAngleCorrection,
        ):
            with self.subTest(optimizer_type=optimizer_type.__name__):
                optimizer = optimizer_type(FitCrystal(), ctrs)
                with self.assertRaisesRegex(
                    ValueError, "<CTR.*kinematical CTR fitting"
                ):
                    optimizer.prepareFit()


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

# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import numpy as np
import copy
import warnings

# from functools import partial
from .. import util

from .CTRcalc import SXRDCrystal
from . import CTRresolution

from scipy import stats

from collections.abc import Callable


class CTROptimizer:
    def __init__(self, xtal, CTRs):
        self.CTRs = copy.deepcopy(CTRs)
        self.CTRs.sort(key=lambda x: abs(x.hk[0]) + abs(x.hk[1]))
        self.xtal = copy.deepcopy(xtal)
        self.scaling = util.get_scale_chi2
        self.resolution = None
        self.resolution_calculation = "sample"
        self._fit_resolution = False
        self._resolution_bounds = None
        self.resolution_errors = None
        self.calculated_CTRs = None
        self._resolution_input_ctrs = None
        self.dw_zconstraints = False
        self.nic = 0
        self.callbacks = []

    def register_fit_callback(
        self,
        function: Callable,
        bounds_low: list,
        bounds_high: list,
        init: list,
        **kwargs,
    ):
        """Register an additional parameterized crystal callback."""
        callback = FitCallback(function, bounds_low, bounds_high, init, **kwargs)
        self.callbacks.insert(0, callback)
        return callback.name

    @property
    def callback_names(self):
        """Names of registered callbacks in parameter-vector order."""
        return [callback.name for callback in self.callbacks]

    def unregister_fit_callback(self, name: str):
        """Remove one registered callback by name."""
        try:
            idx = self.callback_names.index(name)
        except ValueError as error:
            raise ValueError("%s is not a registered callback") from error
        del self.callbacks[idx]

    def _validate_kinematical_input(self):
        """Reject stored quantities unsupported by the kinematical model."""
        for ctr in self.CTRs:
            if ctr.reduction.quantity != "structure_factor":
                raise ValueError(
                    f"{ctr!r}: kinematical CTR fitting supports "
                    "structure-factor data only."
                )

    def set_resolution(self, resolution, calculation=None):
        """Set the L-direction resolution model used for calculated CTRs.

        :param CTRresolution.ResolutionFunction resolution:
            Resolution model evaluated in r.l.u. along each CTR. Pass ``None``
            to use unbroadened structure factors.
        :param str calculation:
            Optional resolution calculation mode: ``"sample"`` evaluates
            structure factors at quadrature points, while ``"convolve"``
            convolves calculated values on the existing L grid.
        :raises TypeError:
            If ``resolution`` is neither ``None`` nor a resolution model.
        """
        if resolution is not None and not isinstance(
            resolution, CTRresolution.ResolutionFunction
        ):
            raise TypeError("resolution must be a ResolutionFunction or None")
        self.resolution = resolution
        if calculation is not None:
            self.set_resolution_calculation(calculation)
        self._invalidate_resolution_cache()

    def set_resolution_calculation(self, calculation):
        """Set how resolution-broadened calculated CTRs are evaluated.

        :param str calculation:
            ``"sample"`` for quadrature sampling of ``crystal.F`` or
            ``"convolve"`` for fast convolution on the existing L grid.
        :raises ValueError:
            If ``calculation`` is not a supported mode.
        """
        if calculation not in {"sample", "convolve"}:
            raise ValueError("calculation must be 'sample' or 'convolve'")
        self.resolution_calculation = calculation
        self._invalidate_resolution_cache()

    def calc_resolution_angles_zmode(
        self, vliegangles, fixedangle=np.deg2rad(0.1), fixed="in", **keyargs
    ):
        """Calculate and cache Z-mode angles used by an angle-dependent model.

        The cached angle records are attached to the optimizer's CTR geometry.
        Their values are in rad and remain valid while the HKL coordinates and
        experimental geometry are unchanged.

        :param HKLVlieg.VliegAngles vliegangles:
            Angle calculator configured for the CTR reference lattice.
        :param float fixedangle:
            Fixed incidence or exit angle in rad.
        :param str fixed:
            Whether ``fixedangle`` fixes the incident (``"in"``) or exit
            (``"out"``) beam angle.
        :returns:
            Cached structured angle records, one per CTR.
        :rtype: list[numpy.recarray]
        """
        angles = self.CTRs.calcAnglesZmode(
            vliegangles, fixedangle=fixedangle, fixed=fixed, **keyargs
        )
        self._resolution_input_ctrs = None
        self._invalidate_resolution_cache()
        return angles

    def fit_resolution(
        self, resolution, lower_bounds, higher_bounds, calculation="sample"
    ):
        """Include all three resolution widths in the fit parameter array.

        Resolution parameters occupy the leading three entries of the
        optimizer parameter array in this order: ``delta_l_0``, ``delta_l_1``,
        and ``delta_l_2``. The widths and their bounds are in r.l.u.

        :param CTRresolution.ResolutionFunction resolution:
            Initial box or Gaussian resolution model.
        :param sequence lower_bounds:
            Three finite lower bounds for the resolution widths in r.l.u.
        :param sequence higher_bounds:
            Three finite upper bounds for the resolution widths in r.l.u.
        :param str calculation:
            ``"sample"`` for quadrature sampling or ``"convolve"`` for
            fast convolution on the measured L grid.
        :raises TypeError:
            If ``resolution`` is not a resolution model.
        :raises ValueError:
            If either bound sequence does not contain three valid bounds.
        """
        self.set_resolution(resolution, calculation=calculation)
        lower_bounds = np.asarray(lower_bounds, dtype=np.float64)
        higher_bounds = np.asarray(higher_bounds, dtype=np.float64)
        if lower_bounds.shape != (3,) or higher_bounds.shape != (3,):
            raise ValueError("Resolution bounds must each contain three values")
        if (
            not np.all(np.isfinite(lower_bounds))
            or not np.all(np.isfinite(higher_bounds))
            or np.any(lower_bounds < 0.0)
            or np.any(lower_bounds > higher_bounds)
        ):
            raise ValueError("Resolution bounds must be finite and nonnegative")
        self._fit_resolution = True
        self._resolution_bounds = (lower_bounds, higher_bounds)

    def _resolution_parameters(self):
        """Return the resolution-width fit parameters in r.l.u."""
        return np.asarray(
            (
                self.resolution.delta_l_0,
                self.resolution.delta_l_1,
                self.resolution.delta_l_2,
            ),
            dtype=np.float64,
        )

    def _set_resolution_parameters(self, parameters):
        """Set the fitted resolution widths from a three-value parameter slice."""
        values = np.asarray(parameters, dtype=np.float64)
        if values.shape != (3,):
            raise ValueError("Resolution fit parameters must contain three values")
        self.resolution = type(self.resolution)(*values)
        self._invalidate_resolution_cache()

    def _invalidate_resolution_cache(self):
        """Discard calculated amplitudes after geometry or model changes."""
        self.calculated_CTRs = None

    def _resolution_input_collection(self):
        """Return the reusable unbroadened input for fast convolution."""
        if self._resolution_input_ctrs is None:
            self._resolution_input_ctrs = copy.deepcopy(self.CTRs)
        return self._resolution_input_ctrs

    def _update_resolution_cache(self):
        """Calculate and cache resolution-broadened amplitudes for all CTRs.

        ``sample_structure_factor`` currently returns a collection rather than
        accepting an output buffer. A flattened-array kernel would avoid this
        remaining collection allocation in future performance work.
        """
        if self.resolution is None:
            self.calculated_CTRs = None
            return
        if self.resolution_calculation == "sample":
            self.calculated_CTRs = CTRresolution.sample_structure_factor(
                self.CTRs, self.xtal, self.resolution
            )
            return

        input_ctrs = self._resolution_input_collection()
        for source, calculated in zip(self.CTRs, input_ctrs):
            calculated.sfI = np.abs(self.xtal.F(source.harr, source.karr, source.l))
        self.calculated_CTRs = CTRresolution.fast_convolve(
            input_ctrs, self.resolution
        )

    def _append_resolution_bounds(self, bounds):
        """Prefix native resolution bounds when resolution fitting is enabled."""
        if not self._fit_resolution:
            return bounds
        return (
            np.concatenate((self._resolution_bounds[0], bounds[0])),
            np.concatenate((self._resolution_bounds[1], bounds[1])),
        )

    def _split_resolution_parameters(self, parameters):
        """Split leading resolution parameters from an optimizer vector."""
        if not self._fit_resolution:
            return parameters
        self._set_resolution_parameters(parameters[:3])
        return parameters[3:]

    def _calculated_amplitude(self, ctr, index):
        """Return the calculated amplitude for one CTR, with resolution."""
        if self.resolution is None:
            return np.abs(self.xtal.F(ctr.harr, ctr.karr, ctr.l))
        if self.calculated_CTRs is None:
            self._update_resolution_cache()
        return self.calculated_CTRs[index].sfI

    def _model_parameters(self):
        """Return the crystal-owned tail of the fit parameter vector."""
        return self.xtal.getInitialParameters()

    def _set_model_parameters(self, parameters):
        """Set parameters owned by the fitted crystal model."""
        self.xtal.setParameters(parameters)

    def _prepend_model_bounds(self, bounds):
        """Add subclass-owned parameters ahead of the crystal bounds."""
        return bounds

    def _set_model_errors(self, errors):
        """Forward the model-owned error tail to the crystal."""
        self.xtal.setFitErrors(errors)

    def _model_parameter_names(self):
        """Return names for subclass and crystal model parameters."""
        return list(self.xtal.fitparnames)

    def _residual_weight(self, ctr):
        """Return the base optimizer's legacy squared-residual weight."""
        return ctr.weight * ctr.err**-2

    def _fit_parameter_names(self):
        """Build names in the same order as the fit parameter vector."""
        names = []
        if self._fit_resolution:
            names.extend(
                (
                    "resolution_delta_l_0",
                    "resolution_delta_l_1",
                    "resolution_delta_l_2",
                )
            )
        for callback in self.callbacks:
            names.extend(callback.parnames)
        names.extend(self._model_parameter_names())
        if len(names) != len(set(names)):
            duplicates = sorted(
                name for name in set(names) if names.count(name) > 1
            )
            raise ValueError(
                "Duplicate fit parameter names are not allowed: "
                + ", ".join(duplicates)
            )
        return names

    def prepareFit(self):
        """Prepare model, callback, constraint, and resolution fit state.

        When a resolution model is configured, this validates any required
        cached angle records and calculates the initial ``calculated_CTRs``
        cache before the optimizer begins evaluating trial parameters.
        """
        self._validate_kinematical_input()
        self.startp, self.lower_bounds, self.higher_bounds = (
            self.xtal.getStartParamAndLimits()
        )
        self.bounds = self._prepend_model_bounds(
            (self.lower_bounds, self.higher_bounds)
        )
        for ctr in self.CTRs:
            ctr.invrelerrsqrd_weight = self._residual_weight(ctr)
        for callback in reversed(self.callbacks):
            self.bounds = (
                np.concatenate((callback.bounds[0], self.bounds[0])),
                np.concatenate((callback.bounds[1], self.bounds[1])),
            )
        self.bounds = self._append_resolution_bounds(self.bounds)
        if self.dw_zconstraints:
            self.nic = self.get_inequalconstraints().size
        else:
            self.nic = 0
        self.fitparnames = self._fit_parameter_names()
        self.priors = self.xtal.priors
        self._update_resolution_cache()

    def get_bounds(self):
        return self.bounds

    def get_parameters(self):
        parameters = self._model_parameters()
        for callback in reversed(self.callbacks):
            parameters = np.concatenate(
                (callback.get_parameters(self.xtal), parameters)
            )
        if self._fit_resolution:
            parameters = np.concatenate((self._resolution_parameters(), parameters))
        return parameters

    def set_parameters(self, x):
        """Set resolution, callback, and model parameters in layout order."""
        x = self._split_resolution_parameters(x)
        counter = 0
        for callback in self.callbacks:
            callback.set_parameters(
                self.xtal, x[counter : counter + callback.n_pars]
            )
            counter += callback.n_pars
        self._set_model_parameters(x[counter:])
        self._update_resolution_cache()

    def set_errors(self, xerror):
        """Split fitted errors across resolution, callbacks, and model."""
        self.errors = xerror
        if self._fit_resolution:
            self.resolution_errors = xerror[:3]
            xerror = xerror[3:]
        counter = 0
        for callback in self.callbacks:
            callback.set_errors(xerror[counter : counter + callback.n_pars])
            counter += callback.n_pars
        self._set_model_errors(xerror[counter:])

    def weighted_residues2(self, x=None):
        if x is not None:
            self.set_parameters(x)
        residues = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            scale = self.scaling(F_theo, ctr.sfI, ctr.err)  # scale CTR
            residues.append(
                (ctr.invrelerrsqrd_weight / scale**2)
                * ((ctr.sfI * scale - F_theo) ** 2)
            )
        return np.concatenate(residues)

    def residues(self, x=None):
        if x is not None:
            self.set_parameters(x)
        residues = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            scale = self.scaling(F_theo, ctr.sfI, ctr.err)  # scale CTR
            residues.append(ctr.sfI * scale - F_theo)
        return np.concatenate(residues)

    def flat_data(self, specular=True):
        dat = []
        err = []
        for i, ctr in enumerate(
            filter(lambda x: specular or x.hk != (0, 0), self.CTRs)
        ):
            dat.append(ctr.sfI)
            err.append(ctr.err)
        return np.concatenate(dat), np.concatenate(err)

    @property
    def nopoints(self):
        F, err = self.flat_data()
        return F.size

    def flat_Fcalc(self, x=None):
        if x is not None:
            self.set_parameters(x)
        F = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            scale = self.scaling(F_theo, ctr.sfI, ctr.err)  # scale CTR
            F.append(F_theo / scale)
        return np.concatenate(F)

    def Rfactor(self, x=None):
        if x is not None:
            self.set_parameters(x)
        residues = []
        Fobs = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            scale = self.scaling(F_theo, ctr.sfI, ctr.err)  # scale CTR
            residues.append(np.abs(ctr.sfI * scale - F_theo))
            Fobs.append(np.abs(ctr.sfI * scale))
        residues = np.sum(np.concatenate(residues))
        Fobs = np.sum(np.concatenate(Fobs))
        return residues / Fobs

    def weighted_residues(self, x=None):
        if x is not None:
            self.set_parameters(x)
        residues = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            scale = self.scaling(F_theo, ctr.sfI, ctr.err)  # scale CTR
            residues.append(
                np.sqrt(ctr.weight) * ((ctr.sfI * scale - F_theo) / (ctr.err * scale))
            )
        return np.concatenate(residues)

    def get_inequalconstraints(self):
        """Return crystal displacement constraints for the current model."""
        if not self.dw_zconstraints:
            return np.array([], dtype=np.float64)
        dwc_enable = self.xtal.getSurfaceDWConstraintEnable()
        sur_basis = self.xtal.getSurfaceBasis()[dwc_enable]
        sur_basis = sur_basis[np.argsort(sur_basis[:, 3])]
        iDWconstr = np.diff(
            self.xtal.uc_bulk.basis[0, 4] - sur_basis[:, 4], prepend=0
        )
        oDWconstr = np.diff(
            self.xtal.uc_bulk.basis[0, 5] - sur_basis[:, 5], prepend=0
        )
        return np.concatenate((iDWconstr, oDWconstr))

    def get_nic(self):
        """Return the number of active inequality constraints."""
        return self.nic

    def fitness(self, x):
        objective = np.sum(self.weighted_residues2(x))
        if self.dw_zconstraints:
            return np.concatenate(
                ([objective], self.get_inequalconstraints())
            )
        return [objective]

    def statistics(self, x):
        self.set_parameters(x)
        residues2 = self.weighted_residues2()

        Rfactor = self.Rfactor()

        stat = dict()

        # variance = np.concatenate([ctr.err**2 for ctr in self.CTRs])
        # varmat_i = np.diag(1/variance)

        chi2_result = np.sum(residues2)
        pvalue = 1 - stats.chi2.cdf(chi2_result, residues2.size - x.size)
        chi2_red = chi2_result / (residues2.size - x.size)

        pcov = util.leastsq_covariance(self.weighted_residues, x)

        self.errors = np.sqrt(np.diag(pcov) * chi2_red)
        self.set_errors(self.errors)
        self.set_parameters(x)

        stat["Chisqr"] = chi2_result
        stat["nodatapoints"] = residues2.size
        stat["Chisqr_red"] = chi2_red
        stat["noparameters"] = x.size
        stat["pvalue"] = pvalue
        stat["Rfactor"] = Rfactor
        stat["covariance"] = pcov

        return stat
    def set_archi_result(self, archi):
        """Convert an archipelago population to an ArviZ fit trace."""
        islandid = int(np.argmin([f[0] for f in archi.get_champions_f()]))
        minisland = archi[islandid]
        pop_min = minisland.get_population()

        res = pop_min.champion_x
        if len(self.fitparnames) != np.asarray(res).size:
            raise ValueError(
                "Fit parameter name count does not match the champion parameter "
                "vector length."
            )

        stat = self.statistics(res)
        popsize = archi[0].get_population().get_f().shape[0]
        params = {
            name: np.empty((len(archi), popsize))
            for name in self.fitparnames
        }
        params["chisqr"] = np.empty((len(archi), popsize))

        for i, island in enumerate(archi):
            population = island.get_population()
            for j, parameter in enumerate(self.fitparnames):
                params[parameter][i] = population.get_x()[:, j]
            params["chisqr"][i] = population.get_f()[:, 0]

        stat.pop("covariance", np.array([]))

        import arviz as az

        return az.from_dict(params, attrs=stat)

    def evaluateStatistics(self, x):
        warnings.warn(
            "usage of evaluateStatistics is deprecated, use CTROptimizer.statistics instead!",  # noqa: E501
            DeprecationWarning,
        )

        self.set_parameters(x)
        residues2 = self.weighted_residues2()

        # variance = np.concatenate([ctr.err**2 for ctr in self.CTRs])
        # varmat_i = np.diag(1/variance)

        chi2_result = np.sum(residues2)

        pvalue = 1 - stats.chi2.cdf(chi2_result, residues2.size - x.size)

        chi2_red = chi2_result / (residues2.size - x.size)

        pcov = util.leastsq_covariance(self.residues, x)

        errors = np.sqrt(np.diag(pcov) * chi2_red)
        if self._fit_resolution:
            self.resolution_errors = errors[:3]
            self.xtal.setFitErrors(errors[3:])
        else:
            self.xtal.setFitErrors(errors)
        self.set_parameters(x)

        return chi2_result, chi2_red, pvalue, residues2.size

    def printStatistics(self, x):
        # chi2_result, chi2_red , pvalue, nodatapoints = self.evaluateStatistics(x)
        stat = self.statistics(x)
        print(
            "Chisqr = {:.4f}, Chisqr_red = {:.4f}, R-factor = {:.4f} ,p-value = {:.6f}, n_refl = {}".format(  # noqa: E501
                stat["Chisqr"],
                stat["Chisqr_red"],
                stat["Rfactor"],
                stat["pvalue"],
                stat["nodatapoints"],
            )
        )

    def get_name(self):
        return "CTR optimizer"

    def setCTRPlotSettings(self, lrange, plotsize, **settings):
        self.lrange = lrange
        self.plotsize = plotsize
        self.settings = settings


class FitCallback:
    global_counter = 1

    def __init__(
        self,
        function: Callable,
        bounds_low: list,
        bounds_high: list,
        init: list,
        **kwargs,
    ):
        """Wrap a callback that contributes parameters to a fit.

        :param Callable function:
            Function taking the :class:`SXRDCrystal` and an array of callback
            parameters.
        :param sequence bounds_low:
            Lower bounds for the callback parameters.
        :param sequence bounds_high:
            Upper bounds for the callback parameters.
        :param sequence init:
            Initial callback parameter values.
        :param str name:
            Optional callback name.
        :param sequence parnames:
            Optional parameter names. Defaults to the callback name for one
            parameter and ``<name>_<index>`` for multiple parameters.
        :raises ValueError:
            If a bound or parameter-name count does not match ``init``.
        """
        self.name = kwargs.get("name", f"default-{FitCallback.global_counter}")
        FitCallback.global_counter += 1

        self.inital = np.asarray(init)
        self.n_pars = self.inital.size
        parnames = kwargs.get("parnames")
        if parnames is None:
            if self.n_pars == 1:
                parnames = [self.name]
            else:
                parnames = [f"{self.name}_{i}" for i in range(self.n_pars)]
        self.parnames = list(parnames)
        if len(self.parnames) != self.n_pars:
            raise ValueError(
                "Number of parameter names does not match number of initial "
                "parameters."
            )

        self.current_values = np.copy(self.inital)

        low_bnds = np.asarray(bounds_low)
        if low_bnds.size != self.n_pars:
            raise ValueError(
                "Number of lower bounds does not match number of initial parameters."
            )
        high_bnds = np.asarray(bounds_high)
        if high_bnds.size != self.n_pars:
            raise ValueError(
                "Number of upper bounds does not match number of initial parameters."
            )

        # (lower, upper), matching the optimizer's own bounds convention:
        # prepareFit prepends bounds[0] to the lower and bounds[1] to the
        # upper array.
        self.bounds = (low_bnds, high_bnds)

        self.function = function

    def __call__(self, xtal: SXRDCrystal, x: list) -> None:
        self.current_values = np.copy(x)
        return self.function(xtal, x)

    def get_parameters(self, xtal: SXRDCrystal = None):
        return self.current_values

    def set_parameters(self, xtal: SXRDCrystal, x: list):
        self.current_values = np.copy(x)
        return self.function(xtal, x)

    def set_errors(self, xerror):
        self.errors = xerror

    def get_bounds(self):
        return self.bounds

    def __repr__(self):
        return f"<{type(self).__name__} : {self.name} ({self.n_pars} pars)>"


class CTROptAngleCorrection(CTROptimizer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.scaleindividual = False
        self.useAnglecorr = False
        self.phasevelocity = 1.0

    def prepareFit(self, phaselim=[0, 2 * np.pi], amplim=[0, 0.75], start=[0.0, 0.0]):
        """Prepare angle-correction and optional resolution fit state.

        Resolution setup validates cached angles and calculates the initial
        ``calculated_CTRs`` cache. Its three native parameters prefix all
        callback, angle-correction, and crystal parameters.
        """
        self.phaselim = phaselim
        self.amplim = amplim
        self.phase, self.amp = start
        super().prepareFit()

    def _model_parameters(self):
        """Prepend active angle-correction values to crystal parameters."""
        parameters = super()._model_parameters()
        if self.useAnglecorr:
            parameters = np.concatenate(([self.phase, self.amp], parameters))
        return parameters

    def _set_model_parameters(self, parameters):
        """Split active angle correction from crystal parameters."""
        if self.useAnglecorr:
            self.phase, self.amp = parameters[:2]
            parameters = parameters[2:]
        super()._set_model_parameters(parameters)

    def _prepend_model_bounds(self, bounds):
        """Prepend active angle-correction bounds to crystal bounds."""
        if not self.useAnglecorr:
            return bounds
        return (
            np.concatenate(([self.phaselim[0], self.amplim[0]], bounds[0])),
            np.concatenate(([self.phaselim[1], self.amplim[1]], bounds[1])),
        )

    def _set_model_errors(self, errors):
        """Drop angle-correction errors before forwarding the crystal tail."""
        if self.useAnglecorr:
            errors = errors[2:]
        super()._set_model_errors(errors)

    def _model_parameter_names(self):
        """Return active angle-correction and crystal parameter names."""
        names = []
        if self.useAnglecorr:
            names.extend(
                (
                    "anglecorrection_phase",
                    "anglecorrection_amplitude",
                )
            )
        names.extend(super()._model_parameter_names())
        return names

    def _residual_weight(self, ctr):
        """Return the angle optimizer's legacy residual weight."""
        return np.sqrt(ctr.weight) / ctr.err

    def applyCorrections(self):
        F_obs = []
        F_t = []
        F_err = []
        if self.useAnglecorr:
            for i, ctr in enumerate(self.CTRs):
                F_theo = self._calculated_amplitude(ctr, i)
                if hasattr(ctr, "angles"):
                    anglecorr = self.get_anglecorrection(ctr.angles["omega"])
                else:
                    anglecorr = 1.0
                ctr *= anglecorr
                if self.scaleindividual or ctr.hk == (0, 0):
                    scale = self.scaling(F_theo, ctr.sfI, ctr.err)  # scale CTR
                    ctr *= scale
                else:
                    F_obs.append(ctr.sfI)
                    F_t.append(F_theo)
                    F_err.append(ctr.err)
            if not self.scaleindividual:
                scale = self.scaling(
                    np.concatenate(F_t), np.concatenate(F_obs), np.concatenate(F_err)
                )
                for i, ctr in enumerate(filter(lambda x: x.hk != (0, 0), self.CTRs)):
                    ctr *= scale
                    ctr.invrelerrsqrd_weight = np.sqrt(ctr.weight) / ctr.err

            self.amp = 0.0
        else:
            warnings.warn("Angle correction was not enabled. Skip applyCorrections.")

    def log_prob(self, x):
        resid, err = self.weighted_residues_errors(x)
        return -0.5 * np.sum(resid**2 + np.log(2 * np.pi * err**2))

    def get_anglecorrection(
        self, omega, x=None
    ):  # improved the fit a bit, but not significantly...
        if self.useAnglecorr:
            if x is not None:
                self.phase, self.amp = x[:2]
            return np.sqrt(
                1 + self.amp * np.sin(self.phasevelocity * (omega + self.phase))
            )
        else:
            return 1.0

    def weighted_residues2(self, x):
        return self.weighted_residues(x) ** 2

    def residues(self, x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            if hasattr(ctr, "angles"):
                anglecorr = self.get_anglecorrection(ctr.angles["omega"])
            else:
                anglecorr = 1.0
            F_obs_corr = ctr.sfI * anglecorr
            F_err_corr = ctr.err * anglecorr
            if self.scaleindividual or ctr.hk == (0, 0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr)  # scale CTR
                residues.append(F_obs_corr * scale - F_theo)
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            return np.concatenate(residues)
        else:
            scale = self.scaling(
                np.concatenate(F_t), np.concatenate(F_obs), np.concatenate(F_err)
            )
            for i, ctr in enumerate(filter(lambda x: x.hk != (0, 0), self.CTRs)):
                residues.append(F_obs[i] * scale - F_t[i])
            return np.concatenate(residues)
        return np.concatenate(residues)

    def Rfactor(self, x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            if hasattr(ctr, "angles"):
                anglecorr = self.get_anglecorrection(ctr.angles["omega"])
            else:
                anglecorr = 1.0
            F_obs_corr = ctr.sfI * anglecorr
            F_err_corr = ctr.err * anglecorr
            if self.scaleindividual or ctr.hk == (0, 0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr)  # scale CTR
                residues.append(F_obs_corr * scale - F_theo)
                F_obs.append(F_obs_corr * scale)
                F_err.append(F_err_corr)
                F_t.append(F_theo)
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            residues = np.concatenate(residues)
            F_obs = np.concatenate(F_obs)
        else:
            scale = self.scaling(
                np.concatenate(F_t), np.concatenate(F_obs), np.concatenate(F_err)
            )
            for i, ctr in enumerate(filter(lambda x: x.hk != (0, 0), self.CTRs)):
                residues.append(F_obs[i] * scale - F_t[i])
                F_obs[i] *= scale
            residues = np.concatenate(residues)
            F_obs = np.concatenate(F_obs)
        residues = np.sum(np.abs(residues))
        return residues / np.sum(np.asarray(F_obs))

    def weighted_residues(self, x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            if hasattr(ctr, "angles"):
                anglecorr = self.get_anglecorrection(ctr.angles["omega"])
            else:
                anglecorr = 1.0
            F_obs_corr = ctr.sfI * anglecorr
            F_err_corr = ctr.err * anglecorr
            if self.scaleindividual or ctr.hk == (0, 0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr)  # scale CTR
                residues.append(
                    (ctr.weight / (F_err_corr * scale)) * (F_obs_corr * scale - F_theo)
                )
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            return np.concatenate(residues)
        else:
            _, err = self.flat_data(False)
            scale = self.scaling(
                np.concatenate(F_t), np.concatenate(F_obs), np.asarray(err)
            )
            for i, ctr in enumerate(filter(lambda x: x.hk != (0, 0), self.CTRs)):
                residues.append(
                    (ctr.weight / (scale * F_err[i])) * (F_obs[i] * scale - F_t[i])
                )
            return np.concatenate(residues)

    def weighted_residues_errors(self, x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        F_err = []
        scaled_errors = []
        for i, ctr in enumerate(self.CTRs):
            F_theo = self._calculated_amplitude(ctr, i)
            if hasattr(ctr, "angles"):
                anglecorr = self.get_anglecorrection(ctr.angles["omega"])
            else:
                anglecorr = 1.0
            F_obs_corr = ctr.sfI * anglecorr
            F_err_corr = ctr.err * anglecorr
            if self.scaleindividual or ctr.hk == (0, 0):
                scale = self.scaling(F_theo, F_obs_corr, F_err_corr)  # scale CTR
                residues.append(
                    (ctr.weight / (F_err_corr * scale)) * (F_obs_corr * scale - F_theo)
                )
                scaled_errors.append(F_err_corr * scale)
            else:
                F_obs.append(F_obs_corr)
                F_t.append(F_theo)
                F_err.append(F_err_corr)
        if self.scaleindividual:
            return np.concatenate(residues), np.concatenate(scaled_errors)
        else:
            _, err = self.flat_data(False)
            scale = self.scaling(
                np.concatenate(F_t), np.concatenate(F_obs), np.asarray(err)
            )
            for i, ctr in enumerate(filter(lambda x: x.hk != (0, 0), self.CTRs)):
                residues.append(
                    (ctr.weight / (scale * F_err[i])) * (F_obs[i] * scale - F_t[i])
                )
                scaled_errors.append(F_err[i] * scale)
            return np.concatenate(residues), np.concatenate(scaled_errors)

    def statistics(self, x=None):
        if x is None:
            x = self.get_parameters()

        # self.xtal.setParameters(x)
        residues2 = self.weighted_residues2(x)

        Rfactor = self.Rfactor(x)

        stat = dict()

        # variance = np.concatenate([ctr.err**2 for ctr in self.CTRs])
        # varmat_i = np.diag(1/variance)

        chi2_result = np.sum(residues2)
        pvalue = 1 - stats.chi2.cdf(chi2_result, residues2.size - x.size)
        chi2_red = chi2_result / (residues2.size - x.size)

        pcov = util.leastsq_covariance(self.weighted_residues, x)

        errors = np.sqrt(np.diag(pcov) * chi2_red)

        self.set_errors(errors)
        self.set_parameters(x)

        stat["Chisqr"] = chi2_result
        stat["nodatapoints"] = residues2.size
        stat["Chisqr_red"] = chi2_red
        stat["noparameters"] = x.size
        stat["pvalue"] = pvalue
        stat["Rfactor"] = Rfactor
        stat["covariance"] = pcov

        return stat

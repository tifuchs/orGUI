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

import copy
import warnings
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from scipy import stats

# from functools import partial
from .. import util

from .CTRcalc import SXRDCrystal
from . import CTRresolution

@dataclass(frozen=True)
class _CTRCalculation:
    """Unscaled values for one measured CTR in its legacy representation."""

    ctr: object
    prediction: np.ndarray
    observation: np.ndarray
    uncertainty: np.ndarray
    angle_correction: object


@dataclass(frozen=True)
class _LegacyScaledCalculation:
    """Temporary scaled-prediction view used by legacy public adapters."""

    values: _CTRCalculation
    scale: float

    @property
    def scaled_prediction(self):
        """Return the analytical scale multiplied by the prediction."""
        return self.values.prediction * self.scale

    @property
    def residual(self):
        """Return observation minus scaled prediction."""
        return self.values.observation - self.scaled_prediction


class _ScaleEstimationError(ValueError):
    """Raised when a fitted analytical scale is invalid."""


_SCALE_POLICIES = frozenset({"fixed", "scaled", "global"})
_SCALE_QUANTITIES = {
    "F": "structure_factor",
    "R": "reflectivity",
}
_SCALE_QUANTITY_KEYS = {value: key for key, value in _SCALE_QUANTITIES.items()}


class CTROptimizer:
    def __init__(self, xtal, CTRs, *, scale_policy=None):
        self.CTRs = copy.deepcopy(CTRs)
        self.CTRs.sort(key=lambda x: abs(x.hk[0]) + abs(x.hk[1]))
        self.xtal = copy.deepcopy(xtal)
        self._scale_policy_defaults = {
            "structure_factor": "scaled",
            "reflectivity": "fixed",
        }
        self._scale_policy_overrides = {}
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
        if scale_policy is not None:
            self.set_scale_policies(scale_policy)

    @staticmethod
    def _validate_scale_policy_value(policy):
        """Return a validated scale-policy value."""
        if not isinstance(policy, str) or policy not in _SCALE_POLICIES:
            raise ValueError(
                "scale policy must be 'fixed', 'scaled', or 'global'"
            )
        return policy

    @staticmethod
    def _scale_quantity(quantity):
        """Translate a public scale quantity selector to stored metadata."""
        try:
            return _SCALE_QUANTITIES[quantity]
        except (KeyError, TypeError) as error:
            raise ValueError("scale quantity must be 'F' or 'R'") from error

    def _resolve_scale_rod_id(self, rodid):
        """Resolve one full or unique shorthand rod identifier."""
        full_matches = [ctr for ctr in self.CTRs if ctr.ctr_id == rodid]
        if len(full_matches) == 1:
            return full_matches[0].ctr_id
        if len(full_matches) > 1:
            raise ValueError(f"Ambiguous CTR identifier {rodid!r}")

        is_hk = (
            isinstance(rodid, tuple)
            and len(rodid) == 2
            and all(np.isscalar(value) for value in rodid)
        )
        if is_hk:
            hk_matches = [ctr for ctr in self.CTRs if ctr.hk == rodid]
            if len(hk_matches) == 1:
                return hk_matches[0].ctr_id
            if len(hk_matches) > 1:
                raise ValueError(
                    f"Ambiguous CTR shorthand {rodid!r}; use the full ctr_id"
                )
        raise ValueError(f"Unknown CTR identifier {rodid!r}")

    def _resolved_scale_policy(
        self, ctr, *, defaults=None, overrides=None
    ):
        """Return one CTR's resolved policy from candidate configuration."""
        if defaults is None:
            defaults = self._scale_policy_defaults
        if overrides is None:
            overrides = self._scale_policy_overrides
        return overrides.get(ctr.ctr_id, defaults[ctr.reduction.quantity])

    def _validate_scale_configuration(self, defaults, overrides):
        """Reject a global group containing unlike stored quantities."""
        global_quantities = {
            ctr.reduction.quantity
            for ctr in self.CTRs
            if self._resolved_scale_policy(
                ctr, defaults=defaults, overrides=overrides
            )
            == "global"
        }
        if len(global_quantities) > 1:
            raise ValueError(
                "The global scale group cannot mix structure-factor (F) and "
                "reflectivity (R) CTRs"
            )

    def set_scale_policies(self, policies):
        """Partially update quantity defaults and per-rod scale policies.

        :param dict policies:
            Mapping with optional ``"F"`` and ``"R"`` defaults plus CTR
            identifiers. Policy values are ``"fixed"``, ``"scaled"``, or
            ``"global"``.
        :raises TypeError:
            If ``policies`` is not a dictionary.
        :raises ValueError:
            If a selector, rod identifier, policy, or resulting global group
            is invalid.
        """
        if not isinstance(policies, dict):
            raise TypeError("scale policies must be provided as a dictionary")
        defaults = self._scale_policy_defaults.copy()
        overrides = self._scale_policy_overrides.copy()
        for key, policy in policies.items():
            policy = self._validate_scale_policy_value(policy)
            if key in _SCALE_QUANTITIES:
                defaults[self._scale_quantity(key)] = policy
            else:
                overrides[self._resolve_scale_rod_id(key)] = policy
        self._validate_scale_configuration(defaults, overrides)
        self._scale_policy_defaults = defaults
        self._scale_policy_overrides = overrides

    def set_scale_policy_default(self, quantity, policy):
        """Set the default scale policy for F or R datasets.

        :param str quantity: ``"F"`` or ``"R"``.
        :param str policy: ``"fixed"``, ``"scaled"``, or ``"global"``.
        """
        self._scale_quantity(quantity)
        self.set_scale_policies({quantity: policy})

    def set_scale_policy(self, rodid, policy):
        """Set an explicit scale policy for one CTR.

        :param rodid:
            A unique ``(h, k)`` shorthand or full ``ctr.ctr_id``.
        :param str policy: ``"fixed"``, ``"scaled"``, or ``"global"``.
        """
        policy = self._validate_scale_policy_value(policy)
        full_id = self._resolve_scale_rod_id(rodid)
        overrides = self._scale_policy_overrides.copy()
        overrides[full_id] = policy
        self._validate_scale_configuration(
            self._scale_policy_defaults, overrides
        )
        self._scale_policy_overrides = overrides

    def get_scale_policy(self, rodid=None):
        """Return the scale configuration or one rod's resolved policy.

        :param rodid:
            Optional unique ``(h, k)`` shorthand or full ``ctr.ctr_id``.
        :returns:
            The complete restorable configuration when ``rodid`` is omitted,
            otherwise the selected rod's resolved policy.
        """
        if rodid is not None:
            full_id = self._resolve_scale_rod_id(rodid)
            ctr = next(ctr for ctr in self.CTRs if ctr.ctr_id == full_id)
            return self._resolved_scale_policy(ctr)
        policies = {
            _SCALE_QUANTITY_KEYS[quantity]: policy
            for quantity, policy in self._scale_policy_defaults.items()
        }
        policies.update(self._scale_policy_overrides)
        return policies

    @property
    def scaleindividual(self):
        """Reject reads of the write-only compatibility setting."""
        raise AttributeError(
            "scaleindividual is write-only; use get_scale_policy()"
        )

    @scaleindividual.setter
    def scaleindividual(self, individual):
        """Redirect legacy grouping assignments to both quantity defaults."""
        if not isinstance(individual, bool | np.bool_):
            raise TypeError("scaleindividual must be a boolean")
        policy = "scaled" if individual else "global"
        self.set_scale_policies({"F": policy, "R": policy})

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

    def _angle_correction(self, ctr):
        """Return the multiplicative legacy correction for one measured CTR."""
        return 1.0

    def _calculation_inputs(self, apply_angle_correction=True):
        """Return shared model, observation, and uncertainty values per CTR.

        ``apply_angle_correction=False`` is a temporary legacy adapter for
        inherited methods whose historical output ignored the subclass angle
        correction. It is private migration state and is removed when the
        common public result contract lands.
        """
        calculations = []
        for index, ctr in enumerate(self.CTRs):
            if apply_angle_correction:
                angle_correction = self._angle_correction(ctr)
            else:
                angle_correction = 1.0
            calculations.append(
                _CTRCalculation(
                    ctr=ctr,
                    prediction=self._calculated_amplitude(ctr, index),
                    observation=ctr.sfI * angle_correction,
                    uncertainty=ctr.err * angle_correction,
                    angle_correction=angle_correction,
                )
            )
        return tuple(calculations)

    def _effective_objective_weight(self, ctr):
        """Return the scale estimator's effective per-point rod weight."""
        return ctr.weight

    def _scale_group_label(self, calculations, policy):
        """Describe a fitted scale group for an actionable error message."""
        if policy == "scaled":
            return f"rod {calculations[0].ctr.ctr_id!r}"
        identifiers = [calculation.ctr.ctr_id for calculation in calculations]
        return f"global group {identifiers!r}"

    def _estimate_scale(self, calculations, policy):
        """Estimate one positive prediction multiplier for a fitted group."""
        weighted_prediction_norm = 0.0
        weighted_overlap = 0.0
        for calculation in calculations:
            weight = self._effective_objective_weight(calculation.ctr)
            inverse_variance = calculation.uncertainty**-2
            weighted_prediction_norm += np.sum(
                weight * calculation.prediction**2 * inverse_variance
            )
            weighted_overlap += np.sum(
                weight
                * calculation.prediction
                * calculation.observation
                * inverse_variance
            )

        label = self._scale_group_label(calculations, policy)
        parameters = np.asarray(self.get_parameters()).tolist()
        if not np.isfinite(weighted_prediction_norm) or (
            weighted_prediction_norm <= 0.0
        ):
            if weighted_prediction_norm == 0.0:
                reason = "weighted prediction norm is zero"
            elif not np.isfinite(weighted_prediction_norm):
                reason = "weighted prediction norm is non-finite"
            else:
                reason = "weighted prediction norm is nonpositive"
            raise _ScaleEstimationError(
                f"Scale estimation failed for {label} at the current "
                f"parameters {parameters}: {reason}. Choose different "
                "starting parameters or use a fixed scale when appropriate."
            )

        scale = weighted_overlap / weighted_prediction_norm
        if not np.isfinite(scale) or scale <= 0.0:
            raise _ScaleEstimationError(
                f"Scale estimation failed for {label} at the current "
                f"parameters {parameters}: analytical scale is nonpositive "
                f"or non-finite (a={scale!r}). Choose different starting "
                "parameters or use a fixed scale when appropriate."
            )
        return scale

    def _scaled_calculations(self, apply_angle_correction=True):
        """Return policy-grouped calculations with prediction-side scales."""
        self._validate_scale_configuration(
            self._scale_policy_defaults, self._scale_policy_overrides
        )
        calculations = self._calculation_inputs(
            apply_angle_correction=apply_angle_correction
        )
        scales = [None] * len(calculations)
        global_indices = []
        try:
            for index, calculation in enumerate(calculations):
                policy = self._resolved_scale_policy(calculation.ctr)
                if policy == "fixed":
                    scales[index] = 1.0
                elif policy == "scaled":
                    scales[index] = self._estimate_scale(
                        (calculation,), policy
                    )
                else:
                    global_indices.append(index)

            if global_indices:
                members = tuple(calculations[index] for index in global_indices)
                global_scale = self._estimate_scale(members, "global")
                for index in global_indices:
                    scales[index] = global_scale
        except _ScaleEstimationError:
            self._invalidate_resolution_cache()
            raise

        return tuple(
            _LegacyScaledCalculation(calculation, scale)
            for calculation, scale in zip(calculations, scales)
        )

    def _legacy_individual_calculations(self):
        """Return policy-scaled values without the angle-correction hook."""
        return self._scaled_calculations(apply_angle_correction=False)

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
            ctr.invrelerrsqrd_weight = ctr.weight * ctr.err**-2
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
        self._scaled_calculations()

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
        return self.weighted_residues(x) ** 2

    def residues(self, x=None):
        if x is not None:
            self.set_parameters(x)
        return np.concatenate(
            [
                calculation.residual
                for calculation in self._legacy_individual_calculations()
            ]
        )

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
        return np.concatenate(
            [
                calculation.scaled_prediction
                for calculation in self._legacy_individual_calculations()
            ]
        )

    def Rfactor(self, x=None):
        if x is not None:
            self.set_parameters(x)
        residues = []
        Fobs = []
        for calculation in self._legacy_individual_calculations():
            residues.append(np.abs(calculation.residual))
            Fobs.append(np.abs(calculation.values.observation))
        residues = np.sum(np.concatenate(residues))
        Fobs = np.sum(np.concatenate(Fobs))
        return residues / Fobs

    def weighted_residues(self, x=None):
        if x is not None:
            self.set_parameters(x)
        residues = []
        for calculation in self._legacy_individual_calculations():
            ctr = calculation.values.ctr
            residues.append(
                np.sqrt(ctr.weight)
                * calculation.residual
                / calculation.values.uncertainty
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
        try:
            objective = np.sum(self.weighted_residues2(x))
        except _ScaleEstimationError:
            objective = np.inf
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

    def _angle_correction(self, ctr):
        """Return the legacy empirical correction for one measured CTR."""
        if hasattr(ctr, "angles"):
            return self.get_anglecorrection(ctr.angles["omega"])
        return 1.0

    def _legacy_angle_calculations(self, corrected_scale_errors):
        """Return policy-scaled, angle-aware temporary calculations.

        :param bool corrected_scale_errors:
            Retained only for compatibility with the temporary adapter's
            increment-5 signature. Scale estimation now follows one path.
        """
        return self._scaled_calculations()

    def applyCorrections(self):
        if self.useAnglecorr:
            calculations = self._scaled_calculations()
            for calculation in calculations:
                ctr = calculation.values.ctr
                ctr *= calculation.values.angle_correction / calculation.scale
                ctr.invrelerrsqrd_weight = ctr.weight * ctr.err**-2

            self.amp = 0.0
        else:
            warnings.warn("Angle correction was not enabled. Skip applyCorrections.")

    def log_prob(self, x):
        try:
            resid, err = self.weighted_residues_errors(x)
        except _ScaleEstimationError:
            return -np.inf
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
        return np.concatenate(
            [
                calculation.residual
                for calculation in self._legacy_angle_calculations(
                    corrected_scale_errors=True
                )
            ]
        )

    def Rfactor(self, x):
        self.set_parameters(x)
        residues = []
        F_obs = []
        F_t = []
        global_calculations = []
        for calculation in self._legacy_angle_calculations(
            corrected_scale_errors=True
        ):
            F_obs.append(calculation.values.observation)
            F_t.append(calculation.values.prediction)
            if self._resolved_scale_policy(calculation.values.ctr) == "global":
                global_calculations.append(calculation)
            else:
                residues.append(calculation.residual)
        if not global_calculations:
            residues = np.concatenate(residues)
            F_obs = np.concatenate(F_obs)
        else:
            scale = global_calculations[0].scale
            for i, calculation in enumerate(global_calculations):
                # Preserve the temporary legacy adapter's collection-index
                # behavior until the subclass override is removed in increment 8.
                residues.append(F_obs[i] - scale * F_t[i])
            residues = np.concatenate(residues)
            F_obs = np.concatenate(F_obs)
        residues = np.sum(np.abs(residues))
        return residues / np.sum(np.asarray(F_obs))

    def weighted_residues(self, x):
        self.set_parameters(x)
        return np.concatenate(
            [
                (
                    np.sqrt(calculation.values.ctr.weight)
                    / calculation.values.uncertainty
                )
                * calculation.residual
                for calculation in self._legacy_angle_calculations(
                    corrected_scale_errors=False
                )
            ]
        )

    def weighted_residues_errors(self, x):
        self.set_parameters(x)
        calculations = self._legacy_angle_calculations(
            corrected_scale_errors=False
        )
        residues = [
            (
                np.sqrt(calculation.values.ctr.weight)
                / calculation.values.uncertainty
            )
            * calculation.residual
            for calculation in calculations
        ]
        supplied_errors = [calculation.values.ctr.err for calculation in calculations]
        return np.concatenate(residues), np.concatenate(supplied_errors)

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

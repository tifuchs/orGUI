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
from .CTRincoherent import IncoherentModel
from . import CTRplotutil, CTRresolution

@dataclass(frozen=True)
class _CTRCalculation:
    """Unscaled values for one measured CTR."""

    ctr: object
    prediction: np.ndarray
    observation: np.ndarray
    uncertainty: np.ndarray
    angle_correction: object


@dataclass(frozen=True)
class _ScaledCalculation:
    """Prediction and fitted analytical scale for one measured CTR."""

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
_DWBA_ANGLE_NAMES = ("alpha", "delta", "gamma", "omega", "chi", "phi")


class CTROptimizer:
    def __init__(self, xtal, CTRs, *, scale_policy=None):
        """Fit one forward model against a CTR collection.

        :param xtal:
            Forward model. Either a coherent ``SXRDCrystal`` or an
            ``IncoherentF2Model`` wrapping one. The model is deep-copied once,
            so later changes to the caller's object do not reach the fit.
        :param CTRplotutil.CTRCollection CTRs:
            Measured CTRs.
        :param scale_policy:
            Optional initial scale policy mapping.
        """
        self.CTRs = copy.deepcopy(CTRs)
        self.CTRs.sort(key=lambda x: abs(x.hk[0]) + abs(x.hk[1]))
        # `self.model` is the fitted forward model and `self.xtal` stays the
        # primary coherent crystal, so registered callbacks, displacement
        # constraints, and public `optimizer.xtal` access keep receiving an
        # SXRDCrystal whether or not the fit is wrapped.
        self.model = copy.deepcopy(xtal)
        if isinstance(self.model, IncoherentModel):
            self.xtal = self.model.coherent_model
        else:
            self.xtal = self.model
        self.n_parameters = None
        self._prepared_signature = None
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
        self._resolution_calculated_ctrs = None
        self._resolution_input_ctrs = None
        self._dwba_enabled = False
        self._dwba_bulk_attenuation = 0.0
        self._prepared = False
        self._dw_zconstraints = False
        self.nic = 0
        self.callbacks = []
        self.errors = None
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
        self._invalidate_calculated_results()

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
        self._invalidate_calculated_results()

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

    def set_dwba(self, enabled=True, *, bulk_attenuation=0.0):
        """Configure live semi-infinite DWBA predictions for fitting.

        The optimizer uses the copied crystal's existing DWBA orientation,
        energy, lattice, and reference transform. Experimental observations
        and uncertainties remain in each CTR's stored representation.

        :param bool enabled: Enable or disable DWBA prediction.
        :param float bulk_attenuation:
            Optional nonnegative empirical exponent per bulk repeat.
        :raises TypeError: If ``enabled`` is not boolean.
        :raises ValueError: If ``bulk_attenuation`` is invalid.
        """
        if not isinstance(enabled, bool | np.bool_):
            raise TypeError("enabled must be boolean")
        attenuation = np.asarray(bulk_attenuation)
        if attenuation.ndim != 0 or not np.isrealobj(attenuation):
            raise ValueError(
                "bulk_attenuation must be a finite nonnegative scalar"
            )
        attenuation = float(attenuation)
        if not np.isfinite(attenuation) or attenuation < 0.0:
            raise ValueError(
                "bulk_attenuation must be a finite nonnegative scalar"
            )
        enabled = bool(enabled)
        if (
            enabled != self._dwba_enabled
            or attenuation != self._dwba_bulk_attenuation
        ):
            self._dwba_enabled = enabled
            self._dwba_bulk_attenuation = attenuation
            self._require_reprepare()

    def get_dwba(self):
        """Return the current DWBA fitting configuration."""
        return {
            "enabled": self._dwba_enabled,
            "bulk_attenuation": self._dwba_bulk_attenuation,
        }

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
        self._require_reprepare()
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
        self._require_reprepare()

    @property
    def dw_zconstraints(self):
        """Whether Debye-Waller displacement constraints are fitted."""
        return self._dw_zconstraints

    @dw_zconstraints.setter
    def dw_zconstraints(self, enabled):
        enabled = bool(enabled)
        if enabled != self._dw_zconstraints:
            self._dw_zconstraints = enabled
            self._require_reprepare()

    def _validate_kinematical_input(self):
        """Reject stored quantities unsupported by the kinematical model."""
        for ctr in self.CTRs:
            if ctr.reduction.quantity != "structure_factor":
                raise ValueError(
                    f"{ctr!r}: kinematical CTR fitting supports "
                    "structure-factor data only."
                )
        if isinstance(self.model, IncoherentModel):
            if self.model.output_quantity != "F2":
                raise ValueError(
                    f"{self.model.model_type} returns "
                    f"{self.model.output_quantity!r}, which the kinematical "
                    "path cannot consume as a squared structure factor"
                )
            self.model.validate(forward_model="kinematical")

    def _kinematic_F2(self, h, k, l):  # noqa: N802,E741
        """Return the forward model's squared structure factor.

        Prefers ``F2``. Crystal-like objects and test doubles which implement
        only ``F`` keep working through the squared-modulus fallback.

        :returns:
            Real squared structure factor in the squared units of ``F``.
        """
        squared = getattr(self.model, "F2", None)
        if callable(squared):
            return np.asarray(squared(h, k, l), dtype=np.float64)
        return np.abs(self.model.F(h, k, l)) ** 2

    def _validate_measurement_input(self):
        """Require aligned finite real data and positive uncertainties."""
        for ctr in self.CTRs:
            arrays = {
                "H": np.asarray(ctr.harr),
                "K": np.asarray(ctr.karr),
                "L": np.asarray(ctr.l),
                "observations": np.asarray(ctr.sfI),
                "uncertainties": np.asarray(ctr.err),
            }
            shape = arrays["L"].shape
            if len(shape) != 1 or any(
                values.shape != shape for values in arrays.values()
            ):
                raise ValueError(
                    f"{ctr!r}: H, K, L, observations, and uncertainties "
                    "must be one-dimensional and point-aligned"
                )
            for name in ("H", "K", "L", "observations"):
                values = arrays[name]
                if not np.isrealobj(values) or not np.all(np.isfinite(values)):
                    raise ValueError(f"{ctr!r}: {name} must be finite real values")
            uncertainty = arrays["uncertainties"]
            if (
                not np.isrealobj(uncertainty)
                or not np.all(np.isfinite(uncertainty))
                or np.any(uncertainty <= 0.0)
            ):
                raise ValueError(
                    f"{ctr!r}: uncertainties must be finite and strictly positive"
                )

    def _validate_dwba_input(self):
        """Validate measurement metadata required by live DWBA prediction."""
        if isinstance(self.model, IncoherentModel):
            # Rejected before any DWBA code runs. Changing the flat surface
            # height changes the optical reference profile and its internal
            # fields, so a kinematical F2 mixture is not a reflectivity and
            # must never be reinterpreted as one.
            raise ValueError(
                f"{self.model.model_type} does not support the DWBA forward "
                f"model: it returns {self.model.output_quantity}, not a "
                "reflectivity"
            )
        if not isinstance(self.xtal, SXRDCrystal):
            raise TypeError("DWBA fitting requires an SXRDCrystal model")
        if self._fit_resolution and self.resolution_calculation == "sample":
            raise ValueError(
                "DWBA fitting does not support fitted-width resolution sampling; "
                "use calculation='convolve' or fixed sampling widths"
            )
        for ctr in self.CTRs:
            if ctr.reduction.polarization is None:
                raise ValueError(
                    f"{ctr!r}: DWBA fitting requires polarization-reduction "
                    "metadata"
                )
            angles = getattr(ctr, "angles", None)
            if angles is None and ctr.scan_geometry is None:
                raise ValueError(
                    f"{ctr!r}: DWBA fitting requires central angle records "
                    "or CTRScanGeometry"
                )
            if (
                self.resolution is not None
                and self.resolution_calculation == "sample"
                and ctr.scan_geometry is None
            ):
                raise ValueError(
                    f"{ctr!r}: DWBA resolution sampling requires "
                    "CTRScanGeometry; use calculation='convolve' with central "
                    "angle records"
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
        if calculation is not None and calculation not in {"sample", "convolve"}:
            raise ValueError("calculation must be 'sample' or 'convolve'")
        self.resolution = resolution
        if calculation is not None:
            self.resolution_calculation = calculation
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
        self._require_reprepare()

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
        self._resolution_calculated_ctrs = None
        self._invalidate_calculated_results()

    def _invalidate_calculated_results(self):
        """Discard the last complete public prediction collection."""
        self.calculated_CTRs = None

    def _require_reprepare(self):
        """Invalidate layout-dependent state after a fit-definition change."""
        if hasattr(self, "_prepared"):
            self._prepared = False
        self._invalidate_resolution_cache()

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
            self._resolution_calculated_ctrs = None
            return
        if self.resolution_calculation == "sample":
            self._resolution_calculated_ctrs = CTRresolution.sample_structure_factor(
                self.CTRs, self.model, self.resolution
            )
            return

        input_ctrs = self._resolution_input_collection()
        # `fast_convolve` squares its input, convolves, and takes one square
        # root, so feeding it sqrt(F2) applies the resolution to F2 before the
        # square root without a second collection-building path.
        for source, calculated in zip(self.CTRs, input_ctrs):
            calculated.sfI = np.sqrt(
                self._kinematic_F2(source.harr, source.karr, source.l)
            )
        self._resolution_calculated_ctrs = CTRresolution.fast_convolve(
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

    def _dwba_measured_angles(self, ctr):
        """Return validated point-aligned measured Vlieg records in rad."""
        angles = getattr(ctr, "angles", None)
        if angles is None:
            return None
        angles = np.asarray(angles)
        if angles.dtype.names != _DWBA_ANGLE_NAMES or angles.shape != ctr.l.shape:
            raise ValueError(
                f"{ctr!r}: central angle records must be point-aligned Vlieg "
                "records with alpha, delta, gamma, omega, chi, and phi"
            )
        if not all(np.all(np.isfinite(angles[name])) for name in _DWBA_ANGLE_NAMES):
            raise ValueError(f"{ctr!r}: central angle records must be finite")
        return angles

    def _dwba_generated_angles(self, ctr, h, k, l):  # noqa: E741
        """Generate point-aligned Vlieg records from one CTR scan rule."""
        geometry = ctr.scan_geometry
        if geometry is None:
            raise ValueError(f"{ctr!r}: CTRScanGeometry is required")
        fixed_angle = 0.0 if geometry.fixed == "eq" else geometry.angle
        try:
            return CTRplotutil._calculate_angles_zmode(
                h,
                k,
                l,
                self.xtal.dwba._vlieg_angles(),
                fixed_angle,
                fixed=geometry.fixed,
                chi=0.0,
                phi=0.0,
                hkl_transform=self.xtal.uc_bulk.refHKLTransform,
                mirrorx=geometry.mirrorx,
            )
        except (TypeError, ValueError) as error:
            raise type(error)(f"{ctr!r}: {error}") from error

    def _dwba_central_angles(self, ctr):
        """Return authoritative measured or generated central angle records."""
        measured = self._dwba_measured_angles(ctr)
        if measured is not None:
            return measured
        return self._dwba_generated_angles(ctr, ctr.harr, ctr.karr, ctr.l)

    def _validate_dwba_central_coordinates(self, ctr, angles):
        """Require every central Vlieg record to reconstruct the measured HKL."""
        polarization = ctr.reduction.polarization
        polarization_i, polarization_f, _ = self.xtal.dwba._polarization_pairs(
            polarization.s_fraction, polarization.outgoing
        )[0]
        prepared = self.xtal.dwba.prepare_from_vlieg(
            *(angles[name] for name in _DWBA_ANGLE_NAMES),
            polarization_i=polarization_i,
            polarization_f=polarization_f,
        )
        expected_hkl = np.vstack((ctr.harr, ctr.karr, ctr.l))
        if not np.allclose(
            prepared.hkl.reshape(3, -1), expected_hkl, rtol=1e-7, atol=1e-8
        ):
            raise ValueError(
                f"{ctr!r}: central angle records do not reproduce every "
                "measured H, K, and L"
            )

    @staticmethod
    def _dwba_beam_directions(angles):
        """Return incident and outgoing unit directions in the sample frame."""
        alpha = angles["alpha"]
        delta = angles["delta"]
        gamma = angles["gamma"]
        omega = angles["omega"]
        cos_omega = np.cos(omega)
        sin_omega = np.sin(omega)

        incident_alpha = np.stack(
            (np.zeros_like(alpha), np.cos(alpha), -np.sin(alpha)), axis=-1
        )
        outgoing_alpha = np.stack(
            (
                np.sin(delta) * np.cos(gamma),
                np.cos(delta) * np.cos(gamma),
                np.sin(gamma),
            ),
            axis=-1,
        )

        def rotate_to_sample(vectors):
            rotated = np.empty_like(vectors)
            rotated[..., 0] = (
                cos_omega * vectors[..., 0] - sin_omega * vectors[..., 1]
            )
            rotated[..., 1] = (
                sin_omega * vectors[..., 0] + cos_omega * vectors[..., 1]
            )
            rotated[..., 2] = vectors[..., 2]
            return rotated

        return rotate_to_sample(incident_alpha), rotate_to_sample(outgoing_alpha)

    def _validate_dwba_sampling_geometry(self, ctr, central_angles):
        """Require a scan rule reproducing the central measured geometry."""
        generated = self._dwba_generated_angles(
            ctr, ctr.harr, ctr.karr, ctr.l
        )
        if getattr(ctr, "angles", None) is None:
            return
        central_directions = self._dwba_beam_directions(central_angles)
        generated_directions = self._dwba_beam_directions(generated)
        equivalent = all(
            np.allclose(central, solved, rtol=1e-7, atol=1e-8)
            for central, solved in zip(central_directions, generated_directions)
        )
        if not equivalent:
            raise ValueError(
                f"{ctr!r}: CTRScanGeometry does not reproduce the central "
                "measured geometry or scattering branch"
            )

    def _dwba_intensity_and_prefactor(self, ctr, angles):
        """Evaluate requested polarization pairs and return R and c."""
        reduction = ctr.reduction.polarization
        state = self.xtal.dwba
        intensity = None
        prefactor = None
        for polarization_i, polarization_f, weight in state._polarization_pairs(
            reduction.s_fraction, reduction.outgoing
        ):
            result = state.evaluate_from_vlieg(
                *(angles[name] for name in _DWBA_ANGLE_NAMES),
                polarization_i=polarization_i,
                polarization_f=polarization_f,
                bulk_mode="semi_infinite",
                bulk_attenuation=self._dwba_bulk_attenuation,
            )
            contribution = weight * np.asarray(result.reflectivity)
            intensity = (
                contribution if intensity is None else intensity + contribution
            )
            if prefactor is None:
                prefactor = np.asarray(result.amplitude_prefactor)
        intensity = np.asarray(intensity, dtype=np.float64)
        if not np.all(np.isfinite(intensity)) or np.any(intensity < 0.0):
            raise ValueError(
                f"{ctr!r}: DWBA field intensity is not nonnegative finite"
            )
        return intensity, prefactor

    def _dwba_prediction(self, ctr):
        """Return one live DWBA prediction in the CTR's stored quantity."""
        central_angles = self._dwba_central_angles(ctr)
        self._validate_dwba_central_coordinates(ctr, central_angles)
        if self.resolution is None:
            intensity, prefactor = self._dwba_intensity_and_prefactor(
                ctr, central_angles
            )
            broadened = intensity
        elif self.resolution_calculation == "convolve":
            intensity, prefactor = self._dwba_intensity_and_prefactor(
                ctr, central_angles
            )
            broadened = CTRresolution.fast_convolve_intensity(
                ctr.harr,
                ctr.karr,
                ctr.l,
                intensity,
                self.resolution,
                central_angles,
            )
        else:
            self._validate_dwba_sampling_geometry(ctr, central_angles)
            sampled_prefactor = {}

            def intensity_at(h, k, l):  # noqa: E741
                angles = self._dwba_generated_angles(ctr, h, k, l)
                values, prefactors = self._dwba_intensity_and_prefactor(ctr, angles)
                points = ctr.l.size
                order = values.size // points
                sampled_prefactor["central"] = np.asarray(prefactors).reshape(
                    points, order
                )[:, order // 2]
                return values

            broadened = CTRresolution.sample_intensity(
                ctr.harr,
                ctr.karr,
                ctr.l,
                intensity_at,
                self.resolution,
                central_angles,
            )
            prefactor = sampled_prefactor["central"]

        if ctr.reduction.quantity == "reflectivity":
            return np.ascontiguousarray(broadened)
        polarization_factor = ctr.reduction.polarization.polarization_factor
        if polarization_factor is None:
            polarization_factor = 1.0
        prediction = np.sqrt(broadened / polarization_factor) / np.abs(prefactor)
        return np.ascontiguousarray(prediction)

    def _calculated_value(self, ctr, index):
        """Return the calculated value for one CTR in its stored quantity."""
        if self._dwba_enabled:
            return self._dwba_prediction(ctr)
        if self.resolution is None:
            return np.sqrt(self._kinematic_F2(ctr.harr, ctr.karr, ctr.l))
        if self._resolution_calculated_ctrs is None:
            self._update_resolution_cache()
        return self._resolution_calculated_ctrs[index].sfI

    def _angle_correction(self, ctr):
        """Return the multiplicative legacy correction for one measured CTR."""
        return 1.0

    def _calculation_inputs(self):
        """Return shared model, observation, and uncertainty values per CTR.
        """
        if self._dwba_enabled:
            with self.xtal.dwba.batch():
                predictions = [
                    self._calculated_value(ctr, index)
                    for index, ctr in enumerate(self.CTRs)
                ]
        else:
            predictions = [
                self._calculated_value(ctr, index)
                for index, ctr in enumerate(self.CTRs)
            ]
        calculations = []
        for ctr, prediction in zip(self.CTRs, predictions):
            angle_correction = self._angle_correction(ctr)
            calculations.append(
                _CTRCalculation(
                    ctr=ctr,
                    prediction=prediction,
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

    def _scaled_calculations(self):
        """Return policy-grouped calculations with prediction-side scales."""
        self._validate_scale_configuration(
            self._scale_policy_defaults, self._scale_policy_overrides
        )
        calculations = self._calculation_inputs()
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
            _ScaledCalculation(calculation, scale)
            for calculation, scale in zip(calculations, scales)
        )

    def _require_prepared(self):
        """Reject evaluation before the current fit layout is prepared."""
        if not self._prepared:
            raise RuntimeError(
                "CTR fit evaluation requires prepareFit() after constructing "
                "or changing the fit definition"
            )

    def _publish_calculated_ctrs(self, calculations):
        """Publish final scaled predictions without observation errors."""
        predicted = copy.deepcopy(self.CTRs)
        for ctr, calculation in zip(predicted, calculations):
            ctr.sfI = np.ascontiguousarray(calculation.scaled_prediction)
            ctr.err = None
            ctr.withErr = False
        name = getattr(self.CTRs, "name", "calculated CTRs")
        self.calculated_CTRs = CTRplotutil.CTRCollection(predicted, name=name)

    def _evaluate(self, x=None):
        """Evaluate and publish the common final-result pipeline."""
        self._require_prepared()
        self._invalidate_calculated_results()
        try:
            if x is not None:
                self.set_parameters(x)
            elif not self._dwba_enabled and self.resolution is not None and (
                self._resolution_calculated_ctrs is None
            ):
                self._update_resolution_cache()
            calculations = self._scaled_calculations()
            self._publish_calculated_ctrs(calculations)
        except Exception:
            self._invalidate_calculated_results()
            raise
        return calculations

    def _model_parameters(self):
        """Return the model-owned tail of the fit parameter vector.

        For a wrapper this is its local block followed by the wrapped
        crystal's, which the wrapper splits internally.
        """
        return self.model.getInitialParameters()

    def _set_model_parameters(self, parameters):
        """Set parameters owned by the fitted forward model."""
        self.model.setParameters(parameters)

    def _prepend_model_bounds(self, bounds):
        """Add subclass-owned parameters ahead of the model bounds."""
        return bounds

    def _set_model_errors(self, errors):
        """Forward the model-owned error tail to the forward model."""
        self.model.setFitErrors(errors)

    def _model_parameter_names(self):
        """Return names for subclass and forward-model parameters."""
        return list(self.model.fitparnames)

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

        This validates the data and parameter layout and publishes the initial
        final prediction in :attr:`calculated_CTRs`.
        """
        self._prepared = False
        self._invalidate_resolution_cache()
        if not self.CTRs:
            raise ValueError("Cannot prepare a fit without CTR data")
        for ctr in self.CTRs:
            if ctr.l.size == 0 or ctr.sfI.size == 0:
                raise ValueError(f"Cannot prepare a fit with empty CTR {ctr!r}")
            if ctr.err is None:
                raise ValueError(f"Cannot prepare a fit without errors for {ctr!r}")
        try:
            self._validate_measurement_input()
            if self._dwba_enabled:
                self._validate_dwba_input()
            else:
                self._validate_kinematical_input()
            # `startp`, `lower_bounds` and `higher_bounds` keep their existing
            # model-block scope. `self.bounds`, `fitparnames`,
            # `get_parameters()` and `n_parameters` describe the full vector.
            self.startp, self.lower_bounds, self.higher_bounds = (
                self.model.getStartParamAndLimits()
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
            # Model-scoped, matching current behavior: resolution, callback
            # and angle-correction prefixes expose no prior API.
            self.priors = self.model.priors
            self._prepared = True
            self.n_parameters = len(self.get_parameters())
            if len(self.fitparnames) != self.n_parameters or len(
                self.bounds[0]
            ) != self.n_parameters:
                raise ValueError(
                    "Prepared fit vector is inconsistent: "
                    f"{self.n_parameters} values, {len(self.fitparnames)} "
                    f"names, {len(self.bounds[0])} bounds"
                )
            self._prepared_signature = self._model_signature()
            self._evaluate()
        except Exception:
            self._prepared = False
            self._invalidate_calculated_results()
            raise

    def _model_signature(self):
        """Return the structural layout of the model parameter block.

        Names, count, and bounds. A change to any of them means the prepared
        vector no longer describes the model, so slicing it would silently
        address the wrong parameters.
        """
        start, lower, upper = self.model.getStartParamAndLimits()
        return (
            tuple(self.model.fitparnames),
            len(start),
            tuple(np.asarray(lower, dtype=np.float64).tolist()),
            tuple(np.asarray(upper, dtype=np.float64).tolist()),
        )

    def _require_current_signature(self):
        """Reject a stale layout after structural parameter changes."""
        if self._prepared_signature is None:
            return
        if self._model_signature() != self._prepared_signature:
            self._prepared = False
            raise ValueError(
                "Fit parameters were added, removed, renamed, or re-bounded "
                "on optimizer.model after prepareFit(). Call prepareFit() "
                "again before evaluating or setting parameters."
            )

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
        self._require_current_signature()
        self._invalidate_calculated_results()
        x = self._split_resolution_parameters(x)
        counter = 0
        for callback in self.callbacks:
            callback.set_parameters(
                self.xtal, x[counter : counter + callback.n_pars]
            )
            counter += callback.n_pars
        self._set_model_parameters(x[counter:])
        if not self._dwba_enabled:
            self._update_resolution_cache()

    def set_errors(self, xerror):
        """Split fitted errors across resolution, callbacks, and model."""
        self.errors = xerror
        if xerror is None:
            self.resolution_errors = None
            for callback in self.callbacks:
                callback.set_errors(None)
            self._set_model_errors(None)
            return
        xerror = np.asarray(xerror)
        if self._fit_resolution:
            self.resolution_errors = xerror[:3]
            xerror = xerror[3:]
        counter = 0
        for callback in self.callbacks:
            callback.set_errors(xerror[counter : counter + callback.n_pars])
            counter += callback.n_pars
        self._set_model_errors(xerror[counter:])

    def weighted_residues2(self, x=None):
        """Return squared weighted residuals in flattened CTR order."""
        return self.weighted_residues(x) ** 2

    def residues(self, x=None):
        """Return observation-minus-prediction values in flattened CTR order."""
        return np.concatenate(
            [calculation.residual for calculation in self._evaluate(x)]
        )

    def flat_data(self, specular=True):
        """Return stored observations and uncertainties in flattened CTR order."""
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
        """Number of stored specular and nonspecular data points."""
        F, err = self.flat_data()
        return F.size

    def flat_prediction(self, x=None, *, specular=True):
        """Return final predictions in flattened CTR order.

        :param array-like x: Optional complete fitted parameter vector.
        :param bool specular: Include the ``(0, 0)`` rod when true.
        :returns: Final F or R predictions after analytical scaling.
        :rtype: numpy.ndarray
        """
        calculations = self._evaluate(x)
        selected = [
            calculation.scaled_prediction
            for calculation in calculations
            if specular or calculation.values.ctr.hk != (0, 0)
        ]
        if not selected:
            return np.empty(0, dtype=np.float64)
        return np.concatenate(selected)

    def flat_Fcalc(self, x=None):
        """Return final F predictions, rejecting collections containing R."""
        if any(
            ctr.reduction.quantity != "structure_factor" for ctr in self.CTRs
        ):
            raise ValueError(
                "flat_Fcalc only supports structure-factor data; use "
                "flat_prediction for reflectivity or mixed collections"
            )
        return self.flat_prediction(x)

    def _rfactor(self, calculations, quantity, label):
        """Return an R factor for one stored quantity selection."""
        selected = [
            calculation
            for calculation in calculations
            if calculation.values.ctr.reduction.quantity == quantity
        ]
        if not selected:
            return None
        denominator = sum(
            np.sum(np.abs(calculation.values.observation))
            for calculation in selected
        )
        if denominator == 0.0:
            warnings.warn(
                f"Cannot calculate {label}: the observed-value denominator is zero",
                RuntimeWarning,
            )
            return None
        numerator = sum(
            np.sum(np.abs(calculation.residual)) for calculation in selected
        )
        return numerator / denominator

    def Rfactor(self, x=None):
        """Return the structure-factor R diagnostic, or ``None`` without F."""
        return self._rfactor(
            self._evaluate(x), "structure_factor", "Rfactor"
        )

    def Rfactor_R(self, x=None):
        """Return the reflectivity R diagnostic, or ``None`` without R."""
        return self._rfactor(self._evaluate(x), "reflectivity", "Rfactor_R")

    def weighted_residues(self, x=None):
        """Return normalized residuals with linear per-rod objective weights."""
        return np.concatenate(
            [
                np.sqrt(calculation.values.ctr.weight)
                * calculation.residual
                / calculation.values.uncertainty
                for calculation in self._evaluate(x)
            ]
        )

    def weighted_residues_errors(self, x=None):
        """Return weighted residuals and raw supplied uncertainties."""
        calculations = self._evaluate(x)
        residues = [
            np.sqrt(calculation.values.ctr.weight)
            * calculation.residual
            / calculation.values.uncertainty
            for calculation in calculations
        ]
        errors = [calculation.values.ctr.err for calculation in calculations]
        return np.concatenate(residues), np.concatenate(errors)

    def log_prob(self, x):
        """Return the Gaussian log likelihood for a parameter vector."""
        try:
            resid, err = self.weighted_residues_errors(x)
        except _ScaleEstimationError:
            return -np.inf
        return -0.5 * np.sum(resid**2 + np.log(2 * np.pi * err**2))

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

    def _fitted_scale_count(self):
        """Return the number of independently fitted analytical scales."""
        scaled = sum(
            self._resolved_scale_policy(ctr) == "scaled" for ctr in self.CTRs
        )
        global_group = any(
            self._resolved_scale_policy(ctr) == "global" for ctr in self.CTRs
        )
        return scaled + int(global_group)

    def statistics(self, x=None):
        """Return fit diagnostics and scaled local covariance.

        Analytical scale parameters count toward the degrees of freedom even
        though they are eliminated from the numerical optimizer vector.
        Covariance and parameter errors are reported only when estimable.
        """
        if x is None:
            x = self.get_parameters()
        x = np.asarray(x, dtype=np.float64)
        calculations = self._evaluate(x)
        weighted = np.concatenate(
            [
                np.sqrt(calculation.values.ctr.weight)
                * calculation.residual
                / calculation.values.uncertainty
                for calculation in calculations
            ]
        )
        chi2_result = np.sum(weighted**2)
        noparameters = x.size + self._fitted_scale_count()
        nu = weighted.size - noparameters
        result = {
            "Chisqr": chi2_result,
            "nodatapoints": weighted.size,
            "Chisqr_red": None,
            "noparameters": noparameters,
            "pvalue": None,
            "Rfactor": self._rfactor(
                calculations, "structure_factor", "Rfactor"
            ),
            "Rfactor_R": self._rfactor(
                calculations, "reflectivity", "Rfactor_R"
            ),
            "covariance": None,
        }
        if nu <= 0:
            warnings.warn(
                "Reduced chi-square and covariance are unavailable because "
                "the fit has no positive degrees of freedom",
                RuntimeWarning,
            )
            self.set_errors(None)
            return result

        chi2_red = chi2_result / nu
        result["Chisqr_red"] = chi2_red
        result["pvalue"] = stats.chi2.sf(chi2_result, nu)

        if x.size == 0:
            result["covariance"] = np.empty((0, 0), dtype=np.float64)
            self.set_errors(np.empty(0, dtype=np.float64))
            return result

        try:
            local_covariance = np.asarray(
                util.leastsq_covariance(self.weighted_residues, x),
                dtype=np.float64,
            )
            if (
                local_covariance.shape != (x.size, x.size)
                or not np.all(np.isfinite(local_covariance))
                or np.linalg.matrix_rank(local_covariance) < x.size
            ):
                raise np.linalg.LinAlgError("rank-deficient covariance")
            reported_covariance = chi2_red * local_covariance
            diagonal = np.diag(reported_covariance)
            if np.any(diagonal < 0.0):
                raise np.linalg.LinAlgError("negative covariance diagonal")
            result["covariance"] = reported_covariance
            self.set_errors(np.sqrt(diagonal))
        except (ValueError, np.linalg.LinAlgError):
            warnings.warn(
                "Fit covariance and parameter errors are unavailable because "
                "the local weighted-residual Jacobian is rank deficient",
                RuntimeWarning,
            )
            self.set_errors(None)
        finally:
            self._evaluate(x)
        return result

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

        attrs = {
            key: value
            for key, value in stat.items()
            if key != "covariance" and value is not None
        }

        import arviz as az

        return az.from_dict(params, attrs=attrs)

    def evaluateStatistics(self, x):
        warnings.warn(
            "evaluateStatistics is deprecated and marked for removal; use "
            "CTROptimizer.statistics instead!",
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
        # Splitting the vector here by hand hardcoded the three-entry
        # resolution prefix and skipped every registered callback, so with a
        # callback present the crystal was handed the callback's error slice.
        # `set_errors` is the same splitter `statistics` uses.
        self.set_errors(errors)
        self.set_parameters(x)

        return chi2_result, chi2_red, pvalue, residues2.size

    def printStatistics(self, x):
        """Print a compact summary of current fit diagnostics."""
        stat = self.statistics(x)
        display = {
            key: "n/a" if value is None else f"{value:.6g}"
            for key, value in stat.items()
            if key != "covariance"
        }
        print(
            f"Chisqr = {display['Chisqr']}, "
            f"Chisqr_red = {display['Chisqr_red']}, "
            f"R-factor(F) = {display['Rfactor']}, "
            f"R-factor(R) = {display['Rfactor_R']}, "
            f"p-value = {display['pvalue']}, "
            f"n_refl = {display['nodatapoints']}"
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
        self._use_anglecorr = False
        self.phasevelocity = 1.0
        self.phase_error = None
        self.amp_error = None

    def set_dwba(self, enabled=True, *, bulk_attenuation=0.0):
        """Reject the unsupported combination of DWBA and angle correction."""
        raise ValueError(
            "DWBA fitting is not supported by CTROptAngleCorrection; use "
            "CTROptimizer"
        )

    @property
    def useAnglecorr(self):
        """Whether empirical angle-correction parameters are fitted."""
        return self._use_anglecorr

    @useAnglecorr.setter
    def useAnglecorr(self, enabled):
        enabled = bool(enabled)
        if enabled != self._use_anglecorr:
            self._use_anglecorr = enabled
            self._require_reprepare()

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
        """Split angle-correction errors from the crystal-owned tail."""
        if errors is None:
            self.phase_error = None
            self.amp_error = None
            super()._set_model_errors(None)
            return
        if self.useAnglecorr:
            self.phase_error, self.amp_error = errors[:2]
            errors = errors[2:]
        else:
            self.phase_error = None
            self.amp_error = None
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

    def applyCorrections(self):
        """Apply the fitted angle correction permanently to stored data."""
        if self.useAnglecorr:
            calculations = self._evaluate()
            for calculation in calculations:
                ctr = calculation.values.ctr
                ctr *= calculation.values.angle_correction / calculation.scale
                ctr.invrelerrsqrd_weight = ctr.weight * ctr.err**-2

            self.amp = 0.0
            self._require_reprepare()
        else:
            warnings.warn("Angle correction was not enabled. Skip applyCorrections.")

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

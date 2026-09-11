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
"""Incoherent kinematical CTR models.

An incoherent model wraps one coherent :class:`CTRcalc.SXRDCrystal` and
returns a squared structure factor for a mixed state. Constructing a wrapper
is the only way to opt in: ``SXRDCrystal.F`` and every existing calculation
stay fully coherent.

Quantities, following ``doc/design/incoherent_ctr_models.md``:

- ``F`` is the coherent complex amplitude and exists only on the crystal. A
  mixed state has no unique complex amplitude, so no wrapper defines ``F``.
- ``F2`` is ``abs(F) ** 2`` in the squared units of ``F``, nominally electrons
  squared per reference lateral cell. It is not detector counts and not a
  reflectivity.

The incoherent average acts on the amplitude of the *complete* coherent
crystal state, not on the isolated surface correction, so bulk-surface and
Film-surface interference stays inside each domain.
"""

__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import inspect
from abc import ABC, abstractmethod
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from .CTRcalc import SXRDCrystal
from .CTRfilm import PoissonSurface
from .CTRutil import LinearFitFunctions

__all__ = [
    "CoherentStateEnsembleModel",
    "IncoherentF2Model",
    "IncoherentModel",
    "IncoherentModelInfo",
    "IncoherentState",
    "KinematicIncoherentContext",
    "PoissonHeightDomains",
    "available_incoherent_models",
    "create_incoherent_model",
    "create_incoherent_model_from_config",
    "register_incoherent_model",
]


@dataclass(frozen=True)
class IncoherentState:
    """One complete coherent state of the wrapped crystal.

    :param int index:
        Position of the state in its source's own ordering.
    :param float probability:
        Normalized weight of this state in the mixture.
    :param amplitude:
        Complete crystal amplitude when this state covers the coherent patch,
        in electrons per reference lateral cell.
    """

    index: int
    probability: float
    amplitude: object


class KinematicIncoherentContext:
    """Lazy, evaluation-local decomposition for one HKL request.

    The coherent decomposition is calculated at most once per context, so the
    semi-infinite bulk and every common component are evaluated once no matter
    how many states a model streams. Nothing is cached beyond the context,
    which is what keeps the cache free of invalidation hazards: optimizer
    parameters, component weights, domain transforms, stacking, attenuation,
    and the coordinates themselves all change between requests.

    :param CTRcalc.SXRDCrystal crystal:
        The wrapped coherent crystal.
    :param numpy.ndarray h:
        Reference-frame reciprocal coordinate in r.l.u.
    :param numpy.ndarray k:
        Reference-frame reciprocal coordinate in r.l.u.
    :param numpy.ndarray l:
        Reference-frame reciprocal coordinate in r.l.u.
    """

    def __init__(self, crystal, h, k, l):  # noqa: E741
        self._crystal = crystal
        self._h = h
        self._k = k
        self._l = l
        self._coherent = None
        self._hkl = None

    @property
    def crystal(self):
        """Return the wrapped coherent crystal."""
        return self._crystal

    @property
    def coordinates(self):
        """Return the requested ``(h, k, l)`` coordinates in r.l.u."""
        return self._h, self._k, self._l

    @property
    def coherent(self):
        """Return the coherent decomposition, evaluating it on first use.

        :rtype: CTRcalc.KinematicAmplitudeResult
        """
        if self._coherent is None:
            self._coherent = self._crystal.evaluate_kinematic(
                self._h, self._k, self._l
            )
        return self._coherent

    def component(self, name):
        """Return one component's amplitude record, selected by stable name.

        :param str name:
            Component name on the wrapped crystal.
        :rtype: CTRcalc.KinematicComponentAmplitude
        """
        return self.coherent.component(name)

    def common_amplitude(self, name):
        """Return the crystal amplitude with one component removed.

        This is the part shared by every state of that component: the bulk,
        the underlying Film, and all other components.

        :param str name:
            Component name on the wrapped crystal.
        :returns:
            Complex amplitude in electrons per reference lateral cell.
        """
        return self.coherent.total - self.component(name).amplitude

    def _stacked_hkl(self):
        """Return the stacked ``(3, N)`` coordinates used by domain transforms."""
        if self._hkl is None:
            self._hkl = np.vstack((self._h, self._k, self._l))
        return self._hkl

    def iter_component_states(self, name, state_evaluator):
        """Yield one complete state amplitude of a component at a time.

        ``state_evaluator`` is called once per outer coherent domain, at that
        domain's transformed coordinates, and must return an object exposing
        ``iter_states()`` over items carrying ``probability`` and
        ``amplitude``. Calling it at the transformed coordinates is required:
        an outer domain changes where a component state has to be evaluated,
        so an amplitude precomputed only at the original coordinates would be
        placed wrongly.

        Every domain must report the same states in the same order. Only the
        amplitudes may differ between domains, because only the coordinates
        do.

        :param str name:
            Component name on the wrapped crystal.
        :param callable state_evaluator:
            Called as ``state_evaluator(component, h, k, l)``.
        :returns:
            Generator of :class:`IncoherentState`.
        :raises ValueError:
            If the domains disagree on the state metadata.
        """
        part = self.component(name)
        common = self.common_amplitude(name)
        component = self._crystal.uc_surface_list[part.index]
        hkl = self._stacked_hkl()

        factors = list(self._crystal._component_domain_factors(part.index))
        sources = []
        for matrix, scale in factors:
            transformed = np.dot(matrix, hkl)
            sources.append(
                (
                    scale,
                    state_evaluator(
                        component,
                        transformed[0],
                        transformed[1],
                        transformed[2],
                    ),
                )
            )

        iterators = [source.iter_states() for _, source in sources]
        scales = [scale for scale, _ in sources]
        for index, group in enumerate(zip(*iterators)):
            reference = group[0]
            for other in group[1:]:
                same_layer = getattr(other, "layer_number", None) == getattr(
                    reference, "layer_number", None
                )
                if not same_layer or not np.isclose(
                    other.probability, reference.probability
                ):
                    raise ValueError(
                        "Coherent domains of component "
                        f"{name!r} disagree on state {index}: only the state "
                        "amplitudes may depend on the domain transform"
                    )
            delta = None
            for scale, state in zip(scales, group):
                scaled = scale * state.amplitude
                delta = scaled if delta is None else delta + scaled
            yield IncoherentState(
                index=index,
                probability=float(reference.probability),
                amplitude=common if delta is None else common + delta,
            )


class IncoherentModel(LinearFitFunctions, ABC):
    """Quantity-neutral wrapper around one coherent crystal.

    The wrapper presents the fit API the optimizer already consumes, with its
    own local parameters first and the wrapped crystal's parameters as the
    tail::

        [incoherent-model-local parameters] [wrapped coherent-crystal parameters]

    ``LinearFitFunctions`` already implements ``getInitialParameters``,
    ``getStartParamAndLimits``, ``setFitErrors``, ``getFitErrors``,
    ``fitparnames``, ``priors`` and ``parameter_list`` as **local-only**
    methods. Every one of them is overridden here as a composite. Inheriting
    one by accident does not raise; it silently returns a vector one block
    short, which is the failure this contract exists to prevent.

    :param CTRcalc.SXRDCrystal crystal:
        The coherent model to wrap.
    :param str name:
        Wrapper name, used as the default prefix for local parameter names.
    """

    model_type = None
    supported_forward_models = frozenset({"kinematical"})
    output_quantity = None

    parameterLookup = {}
    parameterLookup_inv = {}

    def __init__(self, crystal, *, name="incoherent"):
        super().__init__()
        if not isinstance(crystal, SXRDCrystal):
            raise TypeError(
                "An incoherent model wraps one SXRDCrystal, not "
                f"{type(crystal).__name__}"
            )
        self._coherent_model = crystal
        self.name = name

    @property
    def coherent_model(self):
        """Return the wrapped coherent crystal."""
        return self._coherent_model

    # -- parameter block sizes ------------------------------------------

    @property
    def n_local_parameters(self):
        """Return the number of wrapper-local fit parameters."""
        return len(LinearFitFunctions.parameter_list(self))

    @property
    def n_coherent_parameters(self):
        """Return the number of wrapped-crystal fit parameters."""
        return len(self._coherent_model.fitparnames)

    def _split(self, values, what):
        """Split a full vector into its local and coherent blocks."""
        values = np.asarray(values)
        expected = self.n_local_parameters + self.n_coherent_parameters
        if values.shape[0] != expected:
            raise ValueError(
                f"{what} has {values.shape[0]} entries but this model has "
                f"{expected}: {self.n_local_parameters} local and "
                f"{self.n_coherent_parameters} coherent"
            )
        local = self.n_local_parameters
        return values[:local], values[local:]

    # -- composite fit API ----------------------------------------------

    def getStartParamAndLimits(self, force_recalculate=False):  # noqa: N802
        """Return start values and bounds, local block first."""
        local = LinearFitFunctions.getStartParamAndLimits(
            self, force_recalculate
        )
        coherent = self._coherent_model.getStartParamAndLimits(
            force_recalculate
        )
        return tuple(
            np.concatenate(
                (
                    np.asarray(a, dtype=np.float64),
                    np.asarray(b, dtype=np.float64),
                )
            )
            for a, b in zip(local, coherent)
        )

    def getInitialParameters(self, force_recalculate=False):  # noqa: N802
        """Return the full start vector, local block first."""
        return self.getStartParamAndLimits(force_recalculate)[0]

    def setParameters(self, values):  # noqa: N802
        """Set the full parameter vector, local block first.

        The complete input is validated before either block is written, so a
        wrong length leaves the model untouched rather than half updated.

        :param values:
            Full parameter vector.
        :raises ValueError:
            If the length does not match the model.
        """
        local, coherent = self._split(values, "parameter vector")
        if local.size:
            LinearFitFunctions.setFitParameters(self, local)
        self._coherent_model.setParameters(coherent)

    def setFitParameters(self, values):  # noqa: N802
        """Refuse the local-only setter inherited from ``LinearFitFunctions``.

        :raises NotImplementedError:
            Always. ``setParameters`` is the composite entry point; the
            inherited name would write only the local block and silently
            discard the coherent tail.
        """
        raise NotImplementedError(
            f"{type(self).__name__}.setFitParameters would set only the "
            "wrapper-local block. Use setParameters, which covers the local "
            "parameters and the wrapped crystal."
        )

    def setLimits(self, limits):  # noqa: N802
        """Set bounds for both blocks, local block first."""
        local, coherent = self._split(limits, "bounds array")
        if local.size:
            LinearFitFunctions.setLimits(self, local)
        self._coherent_model.setLimits(coherent)

    def setFitErrors(self, errors):  # noqa: N802
        """Set or clear errors on both blocks."""
        if errors is None:
            LinearFitFunctions.setFitErrors(self, None)
            self._coherent_model.setFitErrors(None)
            return
        local, coherent = self._split(errors, "error vector")
        if local.size:
            LinearFitFunctions.setFitErrors(self, local)
        self._coherent_model.setFitErrors(coherent)

    def getFitErrors(self):  # noqa: N802
        """Return errors for both blocks, local block first.

        A block with no fit parameters is skipped rather than queried: the
        existing implementations raise when no errors are set, which for an
        empty block would be a meaningless failure. A non-empty block whose
        errors were never set still raises, as it does today.
        """
        blocks = []
        if self.n_local_parameters:
            blocks.append(np.asarray(LinearFitFunctions.getFitErrors(self)))
        else:
            blocks.append(np.array([], dtype=np.float64))
        if self.n_coherent_parameters:
            blocks.append(np.asarray(self._coherent_model.getFitErrors()))
        else:
            blocks.append(np.array([], dtype=np.float64))
        return np.concatenate(blocks)

    @property
    def fitparnames(self):
        """Return parameter names, local block first."""
        return list(LinearFitFunctions.fitparnames.fget(self)) + list(
            self._coherent_model.fitparnames
        )

    @property
    def priors(self):
        """Return priors, local block first."""
        return list(LinearFitFunctions.priors.fget(self)) + list(
            self._coherent_model.priors
        )

    def parameter_list(self):
        """Return parameter records, local block first."""
        return list(LinearFitFunctions.parameter_list(self)) + list(
            self._coherent_model.parameter_list()
        )

    # -- composite serialization ----------------------------------------

    def _local_parameters_to_dict(self):
        """Return only the wrapper-local parameter subtree."""
        return LinearFitFunctions.parametersToDict(self)

    def _local_parameters_from_dict(self, data, override_values=True):
        """Restore only the wrapper-local parameter subtree."""
        LinearFitFunctions.parametersFromDict(self, data, override_values)

    def parametersToDict(self):  # noqa: N802
        """Return the local subtree and the wrapped crystal's subtree."""
        return {
            "local": self._local_parameters_to_dict(),
            "coherent": self._coherent_model.parametersToDict(),
        }

    def parametersFromDict(self, data, override_values=True):  # noqa: N802
        """Restore both subtrees."""
        self._local_parameters_from_dict(data["local"], override_values)
        self._coherent_model.parametersFromDict(
            data["coherent"], override_values
        )

    def clearParameters(self):  # noqa: N802
        """Clear both blocks."""
        LinearFitFunctions.clearParameters(self)
        self._coherent_model.clearParameters()

    # -- validation and configuration -----------------------------------

    def validate(self, *, forward_model):
        """Check that this model can serve the selected forward model.

        :param str forward_model:
            Name of the calculation the optimizer will run.
        :raises ValueError:
            If the forward model is unsupported or the wrapped crystal does
            not satisfy the model's own requirements.
        """
        if forward_model not in self.supported_forward_models:
            supported = ", ".join(sorted(self.supported_forward_models))
            raise ValueError(
                f"{self.model_type or type(self).__name__} does not support "
                f"the {forward_model} forward model; supported: {supported}"
            )
        self._validate_model(forward_model)

    def _validate_model(self, forward_model):
        """Model-specific validation hook. Override as needed."""

    def _config_settings(self):
        """Return the non-basis settings to persist. Override as needed."""
        return {}

    def to_config(self):
        """Return a plain dictionary describing this model.

        The local parameter subtree is the only stored source for local basis
        values, whether they are fitted or fixed; settings never duplicate
        them. The wrapped crystal is not included: the factory receives it
        separately, so the coherent structure keeps its own canonical file.

        :rtype: dict
        """
        if not self.model_type:
            raise ValueError(
                f"{type(self).__name__} has no model_type and cannot be "
                "serialized; register it first"
            )
        return {
            "type": self.model_type,
            "settings": dict(self._config_settings()),
            "parameters": self._local_parameters_to_dict(),
        }


class IncoherentF2Model(IncoherentModel, ABC):
    """Kinematical incoherent model returning a squared structure factor."""

    output_quantity = "F2"

    def F2(self, h, k, l):  # noqa: N802,E741
        """Return the mixed squared structure factor.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Real, finite, nonnegative squared structure factor in the squared
            units of ``SXRDCrystal.F``.
        :rtype: numpy.ndarray
        :raises ValueError:
            If the model produced a negative or non-finite result.
        """
        context = KinematicIncoherentContext(self._coherent_model, h, k, l)
        value = np.asarray(self._evaluate_F2(context), dtype=np.float64)
        if not np.all(np.isfinite(value)):
            raise ValueError(
                f"{type(self).__name__} produced a non-finite F2"
            )
        if np.any(value < 0.0):
            raise ValueError(
                f"{type(self).__name__} produced a negative F2"
            )
        return value

    @abstractmethod
    def _evaluate_F2(self, context):  # noqa: N802
        """Return the mixed ``F2`` for one evaluation context."""


class CoherentStateEnsembleModel(IncoherentF2Model, ABC):
    """Mix the complete coherent amplitudes of a set of states.

    A subclass supplies normalized probabilities and complete state
    amplitudes; this class owns the validated coherent, incoherent, and
    partially coherent sums. It is deliberately generic: it neither imports
    nor tests for any particular surface model.

    The endpoints are

    .. code-block:: text

        F2_coherent   = abs(sum_n p_n A_n) ** 2
        F2_incoherent = sum_n p_n abs(A_n) ** 2

    and the model returns their convex interpolation. The ``kappa = 0``
    endpoint uses the live coherent evaluation already cached in the context,
    not the state-sum reconstruction, so a finite state support cannot perturb
    it.
    """

    #: Absolute tolerance on the sum of the streamed probabilities.
    probability_tolerance = 1e-9

    @property
    @abstractmethod
    def incoherent_fraction(self):
        """Return the mixing fraction ``kappa`` in ``[0, 1]``."""

    @abstractmethod
    def _iter_states(self, context):
        """Yield :class:`IncoherentState` for one evaluation context."""

    def _validated_fraction(self):
        """Return the mixing fraction after range and finiteness checks."""
        kappa = float(self.incoherent_fraction)
        if not np.isfinite(kappa) or not 0.0 <= kappa <= 1.0:
            raise ValueError(
                f"incoherent_fraction must be finite and within [0, 1], "
                f"got {kappa}"
            )
        return kappa

    def accumulate_states(self, context):
        """Stream the states once, returning the incoherent sum and diagnostics.

        Exactly one state amplitude is held at a time: the squared sum and the
        diagnostic coherent sum are accumulated in place.

        :param KinematicIncoherentContext context:
            Evaluation context.
        :returns:
            Incoherent ``F2``, the reconstructed coherent amplitude, and the
            total streamed probability.
        :rtype: tuple
        :raises ValueError:
            If a probability is invalid or the masses do not sum to one.
        """
        incoherent = None
        reconstructed = None
        mass = 0.0
        count = 0
        for state in self._iter_states(context):
            probability = float(state.probability)
            if not np.isfinite(probability) or probability < 0.0:
                raise ValueError(
                    f"State {count} has an invalid probability {probability}"
                )
            amplitude = np.asarray(state.amplitude)
            if not np.all(np.isfinite(amplitude)):
                raise ValueError(f"State {count} has a non-finite amplitude")
            squared = probability * np.abs(amplitude) ** 2
            weighted = probability * amplitude
            incoherent = squared if incoherent is None else incoherent + squared
            reconstructed = (
                weighted if reconstructed is None else reconstructed + weighted
            )
            mass += probability
            count += 1
        if count == 0:
            raise ValueError(
                f"{type(self).__name__} produced no states to mix"
            )
        if abs(mass - 1.0) > self.probability_tolerance:
            raise ValueError(
                f"State probabilities sum to {mass}, not one; they must be "
                "normalized over the retained states exactly once"
            )
        return incoherent, reconstructed, mass

    def _evaluate_F2(self, context):  # noqa: N802
        """Interpolate between the live coherent and incoherent endpoints."""
        kappa = self._validated_fraction()
        coherent = np.abs(context.coherent.total) ** 2
        if kappa == 0.0:
            # Identical to the wrapped crystal's F2, without streaming any
            # state or triggering a second bulk and Film evaluation.
            return coherent
        incoherent, _, _ = self.accumulate_states(context)
        return (1.0 - kappa) * coherent + kappa * incoherent

    def coherence_reconstruction_error(self, h, k, l):  # noqa: E741
        """Return the diagnostic gap between the two coherent expressions.

        The state-sum reconstruction of the coherent amplitude is a
        diagnostic against the finite-support truncation; it never replaces
        the live coherent endpoint. The excluded probability alone is not an
        amplitude bound, so compare this against a tolerance justified by an
        extended-support reference.

        :returns:
            Maximum absolute difference between the live coherent amplitude
            and the probability-weighted state sum.
        :rtype: float
        """
        context = KinematicIncoherentContext(self._coherent_model, h, k, l)
        _, reconstructed, _ = self.accumulate_states(context)
        return float(
            np.max(np.abs(context.coherent.total - reconstructed))
        )


@dataclass(frozen=True)
class IncoherentModelInfo:
    """Public description of one registered incoherent model.

    :param str key:
        Stable persistence key, taken from the class ``model_type``.
    :param type model_class:
        The registered class.
    :param str description:
        One-line summary from the class docstring.
    :param frozenset supported_forward_models:
        Forward models the class accepts.
    :param str output_quantity:
        Name of the quantity the class returns, such as ``"F2"``.
    """

    key: str
    model_class: type
    description: str
    supported_forward_models: frozenset
    output_quantity: str


_REGISTRY = {}
_BUILTIN_KEYS = set()


def register_incoherent_model(cls=None, *, builtin=False):
    """Register one incoherent model class under its ``model_type``.

    The class attribute is the only key source, so a saved configuration and
    the class cannot disagree. Usable directly or as a decorator.

    :param type cls:
        Concrete :class:`IncoherentModel` subclass.
    :param bool builtin:
        Mark the entry as shipped with the package, which then refuses to be
        replaced by a later registration.
    :returns:
        The class, so this works as a decorator.
    :raises TypeError:
        If the class is not a concrete ``IncoherentModel`` subclass.
    :raises ValueError:
        If the key is missing, already used, or the class has no output
        quantity.
    """
    if cls is None:
        def decorator(inner):
            return register_incoherent_model(inner, builtin=builtin)

        return decorator

    if not (isinstance(cls, type) and issubclass(cls, IncoherentModel)):
        raise TypeError(
            "Only IncoherentModel subclasses can be registered, not "
            f"{cls!r}"
        )
    if inspect.isabstract(cls):
        raise TypeError(
            f"{cls.__name__} is abstract and cannot be registered"
        )
    key = getattr(cls, "model_type", None)
    if not isinstance(key, str) or not key:
        raise ValueError(
            f"{cls.__name__} needs a non-empty model_type to be registered"
        )
    if not getattr(cls, "output_quantity", None):
        raise ValueError(
            f"{cls.__name__} needs a stable output_quantity to be registered"
        )
    if key in _REGISTRY:
        if key in _BUILTIN_KEYS:
            raise ValueError(
                f"{key!r} is a built-in incoherent model and cannot be "
                "replaced"
            )
        raise ValueError(f"An incoherent model is already registered as {key!r}")

    summary = (inspect.getdoc(cls) or "").strip().splitlines()
    _REGISTRY[key] = IncoherentModelInfo(
        key=key,
        model_class=cls,
        description=summary[0] if summary else "",
        supported_forward_models=frozenset(cls.supported_forward_models),
        output_quantity=cls.output_quantity,
    )
    if builtin:
        _BUILTIN_KEYS.add(key)
    return cls


def unregister_incoherent_model(key):
    """Remove one non-built-in registration.

    :param str key:
        Registry key.
    :raises KeyError:
        If the key is not registered.
    :raises ValueError:
        If the key names a built-in model.
    """
    if key not in _REGISTRY:
        raise KeyError(f"No incoherent model registered as {key!r}")
    if key in _BUILTIN_KEYS:
        raise ValueError(f"{key!r} is a built-in incoherent model")
    del _REGISTRY[key]


def available_incoherent_models():
    """Return a read-only mapping from registry key to model description.

    :rtype: types.MappingProxyType
    """
    return MappingProxyType(dict(_REGISTRY))


def _registered(key):
    """Return one registry record, or raise with the available keys."""
    try:
        return _REGISTRY[key]
    except (KeyError, TypeError):
        available = ", ".join(sorted(_REGISTRY)) or "none"
        raise ValueError(
            f"Unknown incoherent model {key!r}. Registered: {available}"
        ) from None


def create_incoherent_model(key, crystal, **keyargs):
    """Build a registered incoherent model around a coherent crystal.

    Direct construction remains the ordinary Python API; this exists so a
    stored type name can be turned into a class without importing an
    arbitrary path from a file.

    :param str key:
        Registry key.
    :param CTRcalc.SXRDCrystal crystal:
        Coherent crystal to wrap.
    :param keyargs:
        Forwarded to the model constructor.
    :rtype: IncoherentModel
    """
    return _registered(key).model_class(crystal, **keyargs)


def create_incoherent_model_from_config(crystal, config):
    """Rebuild a model from :meth:`IncoherentModel.to_config`.

    Only registered keys can be instantiated, and only wrapper-local state is
    restored: the coherent crystal is supplied by the caller and keeps the
    parameters it already has.

    :param CTRcalc.SXRDCrystal crystal:
        Coherent crystal to wrap.
    :param dict config:
        Dictionary produced by ``to_config``.
    :rtype: IncoherentModel
    :raises ValueError:
        If the type is unknown or the dictionary carries unexpected keys.
    """
    if not isinstance(config, dict):
        raise ValueError("An incoherent model configuration must be a dict")
    unexpected = set(config) - {"type", "settings", "parameters"}
    if unexpected:
        raise ValueError(
            "Unexpected keys in incoherent model configuration: "
            + ", ".join(sorted(unexpected))
        )
    if "type" not in config:
        raise ValueError("An incoherent model configuration needs a type")

    record = _registered(config["type"])
    settings = config.get("settings") or {}
    model = record.model_class(crystal, **dict(settings))
    parameters = config.get("parameters")
    if parameters:
        model._local_parameters_from_dict(parameters)
    return model


@register_incoherent_model(builtin=True)
class PoissonHeightDomains(CoherentStateEnsembleModel):
    """Laterally large height domains on a Poisson-roughened surface.

    Each coherence patch sees one flat height drawn from the surface's
    ``PoissonProfile``. The states are the complete coherent amplitudes of the
    crystal at those heights, so bulk-surface and Film-surface interference
    stays inside each domain, and ``incoherent_fraction`` interpolates between
    the coherent limit and their squared average.

    Constructing this wrapper is the opt-in: the ``PoissonSurface`` it targets
    stays an ordinary coherent structural model on its own.

    :param CTRcalc.SXRDCrystal crystal:
        Coherent crystal to wrap.
    :param str surface:
        Stable name of the target ``PoissonSurface`` component.
    :param float incoherent_fraction:
        Dimensionless mixing fraction ``kappa`` in ``[0, 1]``. Zero is fully
        coherent, one is the large-domain limit.
    :param int exact_layer_count:
        Retain every height state when the profile populates no more than
        this many. A policy switch, not a cap.
    :param str name:
        Wrapper name, the default prefix for its local parameter.
    :raises ValueError:
        If the target is missing, duplicated, not a ``PoissonSurface``, not
        the topmost component, or the fraction is out of range.
    """

    model_type = "poisson_height_domains"
    parameterLookup = {"incoherent_fraction": 0}
    parameterLookup_inv = {0: "incoherent_fraction"}

    def __init__(
        self,
        crystal,
        *,
        surface,
        incoherent_fraction=0.0,
        exact_layer_count=10,
        name="incoherent",
    ):
        super().__init__(crystal, name=name)
        self.surface = str(surface)
        if (
            int(exact_layer_count) != exact_layer_count
            or exact_layer_count < 1
        ):
            raise ValueError("exact_layer_count must be a positive integer")
        self.exact_layer_count = int(exact_layer_count)
        # One local basis entry, so a fixed fraction and a fitted one are the
        # same stored value. `basis_0` is seeded too, because
        # `updateFromParameters` rebuilds `basis` from it whenever the
        # fraction is not itself a fit parameter.
        self.basis = np.array([0.0])
        self.basis_0 = np.array([0.0])
        self.errors = None
        self.incoherent_fraction = incoherent_fraction
        self._validate_target()

    @property
    def incoherent_fraction(self):
        """Return the dimensionless mixing fraction held in the basis."""
        return float(self.basis[0])

    @incoherent_fraction.setter
    def incoherent_fraction(self, value):
        """Set the mixing fraction, which must be finite and in ``[0, 1]``."""
        value = float(value)
        if not np.isfinite(value) or not 0.0 <= value <= 1.0:
            raise ValueError(
                f"incoherent_fraction must be finite and within [0, 1], "
                f"got {value}"
            )
        self.basis[0] = value
        self.basis_0[0] = value

    def target_surface(self):
        """Return the targeted ``PoissonSurface`` component.

        Selection is by stable name: component identities do not survive the
        deep copy the optimizer performs once per fit.

        :rtype: CTRfilm.PoissonSurface
        :raises ValueError:
            If no component, or more than one, carries the name.
        """
        components = self._coherent_model.uc_surface_list
        matches = [
            component
            for component in components
            if getattr(component, "name", None) == self.surface
        ]
        if not matches:
            available = ", ".join(
                str(getattr(component, "name", "?"))
                for component in components
            )
            raise ValueError(
                f"No component named {self.surface!r} on the crystal. "
                f"Available: {available}"
            )
        if len(matches) > 1:
            raise ValueError(
                f"{len(matches)} components are named {self.surface!r}; "
                "the target must be identified by a unique name"
            )
        return matches[0]

    def _validate_target(self):
        """Check the target type and the supported stacking topology."""
        target = self.target_surface()
        if not isinstance(target, PoissonSurface):
            raise ValueError(
                f"Component {self.surface!r} is a "
                f"{type(target).__name__}, not a PoissonSurface"
            )

        crystal = self._coherent_model
        ordered = list(crystal.uc_surface_list_ordered)
        if not ordered or ordered[-1] is not target:
            raise ValueError(
                f"{self.surface!r} must be the topmost component in stacking "
                "order. A component above it could have a position which "
                "depends on the selected height, so it cannot be treated as "
                "part of the common amplitude"
            )
        levels = np.asarray(crystal.uc_stacking_ordered)
        if levels.size and int(np.sum(levels == levels[-1])) != 1:
            raise ValueError(
                f"{self.surface!r} shares its stacking level with another "
                "component, whose position could depend on the selected "
                "height"
            )

    def _validate_model(self, forward_model):
        """Check that the target is bound to the Film it corrects."""
        target = self._bound_target()
        if target.underlying_film is None:
            raise ValueError(
                f"{self.surface!r} is not stacked immediately above a Film, "
                "so it has no sharp boundary to correct"
            )

    def _bound_target(self):
        """Return the target after applying any pending stacking."""
        crystal = self._coherent_model
        if crystal.enable_uc_stacking:
            crystal.apply_stacking()
        return self.target_surface()

    def _iter_states(self, context):
        """Stream complete crystal amplitudes, one flat height at a time."""

        def state_evaluator(component, h, k, l):  # noqa: E741
            return component.flat_domain_corrections(
                h, k, l, exact_layer_count=self.exact_layer_count
            )

        return context.iter_component_states(self.surface, state_evaluator)

    def _config_settings(self):
        """Return the non-basis settings. The fraction is not among them."""
        return {
            "surface": self.surface,
            "exact_layer_count": self.exact_layer_count,
        }

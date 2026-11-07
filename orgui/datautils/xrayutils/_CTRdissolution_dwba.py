"""Complete finite-material states evaluated with the native DWBA solver."""

from contextlib import contextmanager

import numpy as np

from ._CTRincoherent_dwba import (
    DWBAEnsembleResult, SurfaceHeightDWBA, _geometry, _phase_to,
    _shared_preparation, _shared_grid, _regrid_model,
)
from .CTRdwba import _readonly_result


class DissolutionDWBA(SurfaceHeightDWBA):
    """Native polarization and geometry entry points for dissolution states.

    Batches snapshot the model and must not contain parameter mutations.
    Every standalone call takes a fresh snapshot, including after an error.
    """

    @contextmanager
    def batch(self):
        """Hold one immutable material snapshot across channels and rods."""
        if self._scope is not None:
            yield self
            return
        path, indices, probabilities, certificate = self._model._snapshot()
        self._scope = dict(path=path, indices=indices, probabilities=probabilities,
                           certificate=certificate, states={})
        try:
            self._scope["average"] = path.averaged_crystal(indices, probabilities)
            with self._scope["average"].dwba.batch():
                yield self
        finally:
            self._scope = None

    def _evaluate_in_batch(self, method, args, kwargs):
        kwargs = dict(kwargs)
        if kwargs.pop("bulk_mode", "semi_infinite") != "semi_infinite":
            raise ValueError("Dissolution protects a semi-infinite bulk")
        attenuation = kwargs.pop("bulk_attenuation", 0.)
        average = self._scope["average"].dwba
        prepared = getattr(average, method)(*args, **kwargs)
        shared = self._model.optical_reference == "shared"
        # All states and the averaged crystal must use the same slab splits.
        # In particular the bare-bulk state introduces a boundary absent from
        # a finite-film profile. A per-state union would change the numerical
        # semi-infinite layer partition between terms of a coherent sum.
        if shared:
            boundaries = []
            for index, probability in zip(self._scope["indices"],
                                          self._scope["probabilities"]):
                if not probability:
                    continue
                states = self._scope["states"]
                if int(index) not in states:
                    crystal = self._scope["path"].state(int(index)).crystal
                    states[int(index)] = crystal.dwba
                profile = states[int(index)]._atomic_model_snapshot()[2]
                boundaries.extend(profile.boundaries)
            _, atomic, profile = average._atomic_model_snapshot()
            grid = _shared_grid(prepared.reference.profile, np.asarray(boundaries))
            atomic = _regrid_model(atomic, profile, grid.values[:, 0])
            prepared = average._prepare_explicit(
                _geometry(prepared), average._reference(grid), atomic,
                prepared.polarization_i, prepared.polarization_f, 1e-7, 1e-8)
        coherent = average.evaluate_prepared(prepared, bulk_attenuation=attenuation)
        amplitude_sum = np.zeros(prepared.shape, dtype=complex)
        incoherent = np.zeros(prepared.shape, dtype=float)
        for index, probability in zip(self._scope["indices"],
                                      self._scope["probabilities"]):
            if not probability:
                continue
            states = self._scope["states"]
            if int(index) not in states:
                states[int(index)] = self._scope["path"].state(int(index)).crystal.dwba
            state = states[int(index)]
            if shared:
                _, atomic, profile = state._atomic_model_snapshot()
                state_prepared = _shared_preparation(state, prepared, atomic, profile)
            else:
                state_prepared = state._prepare_geometry(
                    _geometry(prepared), polarization_i=prepared.polarization_i,
                    polarization_f=prepared.polarization_f, geometry_rtol=1e-7,
                    geometry_atol=1e-8, cache=False)
            result = state.evaluate_prepared(state_prepared,
                                             bulk_attenuation=attenuation)
            phase = _phase_to(state_prepared, prepared)
            field = (coherent.unperturbed_amplitude + coherent.amplitude_prefactor
                     * (result.F_atomic*phase - coherent.F_reference)
                     if shared else result.total_amplitude*phase)
            amplitude_sum += probability*field
            incoherent += probability*np.abs(field)**2
        coherent_intensity = (coherent.reflectivity if shared
                              else np.abs(amplitude_sum)**2)
        fraction = self._model._validated_fraction()
        intensity = (1-fraction)*coherent_intensity + fraction*incoherent
        if not np.all(np.isfinite(intensity)) or np.any(intensity < 0):
            raise ValueError(
                "Dissolution DWBA intensity must be finite and nonnegative")

        def freeze(value):
            return _readonly_result(value, prepared.shape, prepared.scalar)

        return DWBAEnsembleResult(
            prepared, freeze(intensity), coherent.amplitude_prefactor,
            freeze(coherent_intensity), freeze(incoherent),
            self._model.optical_reference)

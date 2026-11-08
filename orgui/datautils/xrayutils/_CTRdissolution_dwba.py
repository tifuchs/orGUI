"""Complete finite-material states evaluated with the native DWBA solver."""

from contextlib import contextmanager, ExitStack
from dataclasses import replace

import numpy as np

from ._CTRincoherent_dwba import (
    DWBAEnsembleResult, SurfaceHeightDWBA, _geometry, _phase_to,
    _shared_grid, _regrid_model,
)
from .CTRdwba import _amplitude_prefactor, _readonly_result


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
            with ExitStack() as stack:
                average = self._scope["average"].dwba
                stack.callback(average.clear_cache)
                stack.enter_context(average.batch())
                for index, probability in zip(indices, probabilities):
                    if not probability:
                        continue
                    state = path.state(int(index)).crystal.dwba
                    self._scope["states"][int(index)] = state
                    stack.callback(state.clear_cache)
                    stack.enter_context(state.batch())
                    # The native cache keys contain the complete cell optics.
                    # Ownership lasts only for this immutable ensemble batch.
                    state._batch_profiles = average._batch_profiles
                yield self
        finally:
            self._scope = None

    def _shared_prepared(self, prepared):
        """Rebin all atomic models once and share the union-grid wavefields."""
        average = self._scope["average"].dwba
        union = self._scope.get("union")
        fingerprint = prepared.reference.fingerprint
        if union is None or union[0] != fingerprint:
            models = {index: state._atomic_model_snapshot()[1:]
                      for index, state in self._scope["states"].items()}
            _, atomic, profile = average._atomic_model_snapshot()
            if (prepared.reference.continuous is not None
                    or any(p.continuous is not None for _, p in models.values())):
                raise NotImplementedError(
                    "Shared dissolution fields require a layered optical reference")
            boundaries = np.concatenate([p.boundaries for _, p in models.values()])
            grid = _shared_grid(prepared.reference.profile, boundaries)
            atomic = _regrid_model(atomic, profile, grid.values[:, 0])
            rebinned = {index: _regrid_model(model, p, grid.values[:, 0])
                        for index, (model, p) in models.items()}
            union = fingerprint, average._reference(grid), atomic, rebinned
            self._scope["union"] = union
            self._scope.pop("shared_prepared", None)
        cached = self._scope.get("shared_prepared")
        if cached is None or cached[0] is not prepared:
            common = average._prepare_explicit(
                _geometry(prepared), union[1], union[2],
                prepared.polarization_i, prepared.polarization_f, 1e-7, 1e-8)
            # Keep one ray block, not an unbounded cache of full preparations.
            self._scope["shared_prepared"] = prepared, common
        return self._scope["shared_prepared"][1], union[3]

    def _evaluate_in_batch(self, method, args, kwargs):
        kwargs = dict(kwargs)
        if kwargs.pop("bulk_mode", "semi_infinite") != "semi_infinite":
            raise ValueError("Dissolution protects a semi-infinite bulk")
        attenuation = kwargs.pop("bulk_attenuation", 0.)
        cache_prepared = bool(kwargs.get("cache", True))
        average = self._scope["average"].dwba
        prepared = getattr(average, method)(*args, **kwargs)
        shared = self._model.optical_reference == "shared"
        if shared:
            prepared, models = self._shared_prepared(prepared)
            coherent = average.evaluate_prepared(prepared,
                                                 bulk_attenuation=attenuation)
        prefactor = _amplitude_prefactor(prepared)
        amplitude_sum = np.zeros(prepared.shape, dtype=complex)
        incoherent = np.zeros(prepared.shape, dtype=float)
        for index, probability in zip(self._scope["indices"],
                                      self._scope["probabilities"]):
            if not probability:
                continue
            state = self._scope["states"][int(index)]
            if shared:
                # Every state has the same bulk frame and all-state union grid.
                # Only its native record decomposition and owner tokens differ.
                state_prepared = replace(
                    prepared, crystal_identity=id(state._crystal),
                    bulk_identity=id(state._crystal.uc_bulk),
                    _atomic_model=models[int(index)])
            else:
                state_prepared = state._prepare_geometry(
                    _geometry(prepared), polarization_i=prepared.polarization_i,
                    polarization_f=prepared.polarization_f, geometry_rtol=1e-7,
                    geometry_atol=1e-8, cache=cache_prepared)
            result = state.evaluate_prepared(state_prepared,
                                             bulk_attenuation=attenuation)
            phase = _phase_to(state_prepared, prepared)
            field = (coherent.unperturbed_amplitude + prefactor
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
            prepared, freeze(intensity), prefactor,
            freeze(coherent_intensity), freeze(incoherent),
            self._model.optical_reference)

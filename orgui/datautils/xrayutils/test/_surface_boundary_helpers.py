"""Native stack builders shared by physical boundary regression tests."""

import copy

import numpy as np

from ..CTRdistributions import HeightSupport
from ..CTRfilm import FiniteFilmBoundary
from ..CTRincoherent import SurfaceHeightDomains


def _bounded_model(crystal, *, films, surface, profile=None,
                   incoherent_fraction=0., optical_reference="state_specific",
                   exposure_policy="native_or_unreconstructed", **budgets):
    target = crystal[surface]
    if profile is not None:
        target.profile = profile
        target.basis = np.array([getattr(profile, key) for key in target._basis_fields])
        target.basis_0 = target.basis.copy()
    target.boundary = FiniteFilmBoundary(tuple(films), exposure_policy=exposure_policy)
    target.height_support = HeightSupport(**budgets)
    return SurfaceHeightDomains(crystal, surface=surface,
                                incoherent_fraction=incoherent_fraction,
                                optical_reference=optical_reference)


def _path(crystal, *, films, surface, exposure_policy="native_or_unreconstructed"):
    native = copy.deepcopy(crystal)
    _bounded_model(native, films=films, surface=surface,
                   exposure_policy=exposure_policy)
    return native._surface_state_plan(surface)

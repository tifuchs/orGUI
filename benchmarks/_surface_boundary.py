"""Build independent native boundary fixtures for paired benchmarks."""

import copy

import numpy as np

from orgui.datautils.xrayutils.CTRdistributions import HeightSupport
from orgui.datautils.xrayutils.CTRfilm import FiniteFilmBoundary
from orgui.datautils.xrayutils.CTRincoherent import SurfaceHeightDomains


def _bounded_model(crystal, *, films, surface, profile=None,
                   incoherent_fraction=0., optical_reference="state_specific",
                   exposure_policy="native_or_unreconstructed", **budgets):
    crystal = copy.deepcopy(crystal)
    target = crystal[surface]
    if profile is not None:
        target.profile = profile
        target.basis = np.array([getattr(profile, key) for key in target._basis_fields])
        target.basis_0 = target.basis.copy()
    profile = target._make_profile()
    budgets.setdefault("tail_probability", profile.tail_probability)
    budgets.setdefault("moment_tolerance", getattr(profile, "moment_tolerance", 1e-10))
    target.boundary = FiniteFilmBoundary(tuple(films), exposure_policy=exposure_policy)
    target.height_support = HeightSupport(**budgets)
    return SurfaceHeightDomains(crystal, surface=surface,
                                incoherent_fraction=incoherent_fraction,
                                optical_reference=optical_reference)

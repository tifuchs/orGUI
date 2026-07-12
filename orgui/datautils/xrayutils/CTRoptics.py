"""Optical constants derived from CTR structural models."""

import numpy as np

from .CTRutil import atomic_number


def optical_profile(unit_cell, include_coherent_domains=True):
    """Return homogeneous optical constants for every layer and domain.

    The returned C-contiguous ``(N, 3)`` float64 array has columns ``z`` in
    Angstrom, ``delta``, and ``beta``.  Each row is a domain-transformed layer
    contribution; coherent-domain occupancy is already applied to ``delta``
    and ``beta`` in the convention ``n = 1-delta-i*beta``.

    :param UnitCell unit_cell:
        Unit cell with energy-dependent scattering factors populated by
        :meth:`UnitCell.setEnergy`.
    :param bool include_coherent_domains:
        Apply the unit cell's coherent-domain transforms and occupancies.
    :returns:
        One row for every structural layer and coherent domain.
    :rtype: numpy.ndarray
    :raises ValueError:
        If the unit-cell energy has not populated anomalous scattering factors.
    """
    if not hasattr(unit_cell, "_E") or not hasattr(unit_cell, "f"):
        raise ValueError(
            "Set the UnitCell energy before calculating optical constants."
        )

    layers = unit_cell.split_in_layers(ordered=True)
    layer_ids = tuple(layers)
    if not layer_ids:
        return np.empty((0, 3), dtype=np.float64)

    if hasattr(unit_cell, "_optical_layer_origin"):
        if len(layer_ids) != 1:
            raise ValueError("Layer optical metadata requires exactly one layer.")
        origins = np.asarray([unit_cell._optical_layer_origin], dtype=np.float64)
        thickness_fraction = np.asarray(
            [unit_cell._optical_layer_thickness_fraction], dtype=np.float64
        )
    else:
        origins = []
        for layer_id in layer_ids:
            origin = unit_cell.layerpos.get(float(layer_id))
            if origin is None:
                atom_mask = unit_cell.basis[:, 7] == layer_id
                origin = float(np.mean(unit_cell.basis[atom_mask, 3]))
            origins.append(float(origin))
        origins = np.asarray(origins, dtype=np.float64)
        thickness_fraction = np.empty(len(layer_ids), dtype=np.float64)
        for index in range(len(layer_ids)):
            spacing = origins[(index + 1) % len(layer_ids)] - origins[index]
            while spacing <= 0.0:
                spacing += 1.0
            thickness_fraction[index] = spacing

    wavelength = 12.398419843320026 / (unit_cell._E * 1e-3)
    scale = 2.8179403262e-5 * wavelength**2 / (2.0 * np.pi)
    if include_coherent_domains:
        domain_matrices = unit_cell.coherentDomainMatrix
        domain_occupancies = unit_cell.coherentDomainOccupancy
    else:
        domain_matrices = [np.hstack((np.identity(3), np.zeros((3, 1))))]
        domain_occupancies = [1.0]
    profile = np.empty((len(layer_ids) * len(domain_matrices), 3), dtype=np.float64)
    row = 0
    for layer_index, (layer_id, layer_cell) in enumerate(layers.items()):
        atomic_numbers = np.asarray(
            [atomic_number(name) for name in layer_cell.names], dtype=np.float64
        )
        forward_factor = np.sum(
            layer_cell.basis[:, 6]
            * (atomic_numbers + layer_cell.f[:, 11] + 1j * layer_cell.f[:, 12])
        )
        layer_volume = unit_cell.volume * thickness_fraction[layer_index]
        optical_density = scale * forward_factor / layer_volume
        origin_vector = np.array([0.0, 0.0, origins[layer_index]])
        for matrix, occupancy in zip(
            domain_matrices, domain_occupancies
        ):
            domain_matrix = unit_cell.R_mat_inv @ matrix[:, :-1] @ unit_cell.R_mat
            domain_origin = domain_matrix @ origin_vector + matrix[:, -1]
            contribution = float(occupancy) * optical_density
            profile[row] = [
                domain_origin[2] * unit_cell.a[2],
                contribution.real,
                contribution.imag,
            ]
            row += 1
    return profile


def optical_profile_asbulk(unit_cell, noUC=30):
    """Return the semi-infinite bulk approximation used by CTR densities.

    This mirrors :meth:`UnitCell.zDensity_G_asbulk`: the native bulk unit
    cell is repeated towards negative z, with ``noUC`` cells explicitly
    represented.  Coherent-domain transforms are intentionally not used for
    this bulk baseline.

    :param UnitCell unit_cell:
        Bulk unit cell with energy-dependent scattering factors populated.
    :param int noUC:
        Number of unit cells to represent below the termination.
    :returns:
        C-contiguous ``(N, 3)`` array with columns ``z`` in Angstrom,
        ``delta``, and ``beta``.
    :rtype: numpy.ndarray
    :raises ValueError:
        If ``noUC`` is not a positive integer.
    """
    if not isinstance(noUC, int | np.integer) or noUC <= 0:
        raise ValueError("noUC must be a positive integer.")

    unit_profile = optical_profile(unit_cell, include_coherent_domains=False)

    profiles = []
    for number in range(noUC):
        profile = unit_profile.copy()
        profile[:, 0] -= number * unit_cell.a[2]
        profiles.append(profile)
    return np.ascontiguousarray(np.concatenate(profiles, axis=0))


def water_optical_profile(water_model, noUC=30, z_step=0.6):
    """Sample continuous water optical constants towards positive z.

    The sampled density is calculated directly from
    :meth:`WaterModel.zDensity_G` at ``h = k = 0``.  This is the continuous
    analogue of :func:`optical_profile_asbulk`, whose repeated bulk cells
    extend towards negative z.

    :param WaterModel water_model:
        Water model with energy-dependent scattering factors populated.
    :param int noUC:
        Length of the sampled water region in water-model unit cells.
    :param float z_step:
        Uniform z sampling interval in Angstrom.
    :returns:
        C-contiguous ``(N, 3)`` array with columns ``z`` in Angstrom,
        ``delta``, and ``beta``.
    :rtype: numpy.ndarray
    :raises ValueError:
        If energy, sample length, or sample spacing is invalid.
    """
    if not hasattr(water_model, "_E") or not hasattr(water_model, "f"):
        raise ValueError(
            "Set the WaterModel energy before calculating optical constants."
        )
    if not isinstance(noUC, int | np.integer) or noUC <= 0:
        raise ValueError("noUC must be a positive integer.")
    if z_step <= 0.0:
        raise ValueError("z_step must be positive.")

    z_start = water_model.pos_absolute + water_model.basis[0] * water_model.a[2]
    z_stop = z_start + noUC * water_model.a[2]
    z = np.arange(z_start, z_stop + 0.5 * z_step, z_step, dtype=np.float64)
    rho_f = water_model.zDensity_G(z, 0.0, 0.0)
    wavelength = 12.398419843320026 / (water_model._E * 1e-3)
    scale = 2.8179403262e-5 * wavelength**2 / (2.0 * np.pi)
    profile = np.empty((z.size, 3), dtype=np.float64)
    profile[:, 0] = z
    profile[:, 1] = scale * rho_f.real
    profile[:, 2] = scale * rho_f.imag
    return np.ascontiguousarray(profile)


def add_structural_to_sampled_profile(
    structural_profile, *sampled_profiles, z_tolerance=0.6
):
    """Add discrete structural layers to a continuous sampled profile.

    A continuous profile must retain each of its sampling points.  Therefore
    it cannot use :func:`combine_profiles`, which coalesces nearby origins.
    Instead, each structural contribution within ``z_tolerance`` Angstrom is
    added to the nearest sampled point; structural contributions outside the
    sampled region keep their original position.

    :param numpy.ndarray structural_profile:
        Discrete layer contributions with columns z, delta, and beta.
    :param sampled_profiles:
        Continuous sampled profiles with the same columns.
    :param float z_tolerance:
        Maximum structural-to-sample separation in Angstrom.
    :returns:
        C-contiguous combined profile sorted by z.
    :rtype: numpy.ndarray
    """
    if z_tolerance < 0.0:
        raise ValueError("z_tolerance must be non-negative.")
    samples = [profile for profile in sampled_profiles if profile.size]
    if not samples:
        return np.ascontiguousarray(structural_profile)
    sampled = np.ascontiguousarray(np.concatenate(samples, axis=0), dtype=np.float64)
    sampled = sampled[np.argsort(sampled[:, 0], kind="stable")]
    remainder = []
    for row in structural_profile:
        right = np.searchsorted(sampled[:, 0], row[0])
        candidates = [
            index for index in (right - 1, right) if 0 <= index < len(sampled)
        ]
        nearest = min(candidates, key=lambda index: abs(sampled[index, 0] - row[0]))
        if abs(sampled[nearest, 0] - row[0]) <= z_tolerance:
            sampled[nearest, 1:] += row[1:]
        else:
            remainder.append(row)
    if remainder:
        profile = np.concatenate((np.asarray(remainder), sampled), axis=0)
    else:
        profile = sampled
    return np.ascontiguousarray(profile[np.argsort(profile[:, 0], kind="stable")])


def combine_profiles(*profiles, z_tolerance=0.6):
    """Combine and coalesce homogeneous optical-profile contributions.

    :param profiles:
        Arrays returned by :func:`optical_profile`, optionally empty.
    :param float z_tolerance:
        Maximum separation in Angstrom between layer origins that represent
        one homogeneous optical layer.  The lowest origin in each merged
        group is retained as its coordinate.
    :returns:
        C-contiguous ``(N, 3)`` array sorted by z, with delta and beta summed
        for origins separated by at most ``z_tolerance``.
    :rtype: numpy.ndarray
    """
    if z_tolerance < 0.0:
        raise ValueError("z_tolerance must be non-negative.")
    populated = [profile for profile in profiles if profile.size]
    if not populated:
        return np.empty((0, 3), dtype=np.float64)
    stacked = np.ascontiguousarray(np.concatenate(populated, axis=0), dtype=np.float64)
    order = np.argsort(stacked[:, 0], kind="stable")
    stacked = stacked[order]
    group_starts = np.r_[0, np.flatnonzero(np.diff(stacked[:, 0]) > z_tolerance) + 1]
    group_stops = np.r_[group_starts[1:], len(stacked)]
    combined = np.empty((len(group_starts), 3), dtype=np.float64)
    for index, (start, stop) in enumerate(zip(group_starts, group_stops)):
        group = stacked[start:stop]
        combined[index] = [group[0, 0], group[:, 1].sum(), group[:, 2].sum()]
    return np.ascontiguousarray(combined)

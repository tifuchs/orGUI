"""Optical constants derived from CTR structural models."""

from dataclasses import dataclass

import numpy as np

from .CTRutil import atomic_number


HC_KEV_ANGSTROM = 12.398419843320026


@dataclass(frozen=True)
class Wavefield:
    """Unperturbed one-dimensional wavefield in Renaud notation.

    Slabs are ordered from the incident medium towards the substrate.
    ``A_plus`` and ``A_minus`` are downward- and upward-propagating electric
    field amplitudes referenced to ``z_reference`` for each slab.

    :param numpy.ndarray z:
        Sample positions in Angstrom, ordered from substrate to ambient.
    :param numpy.ndarray psi:
        Complex electric field sampled at ``z``.
    :param numpy.ndarray z_interfaces:
        Interface positions in Angstrom, ordered from ambient to substrate.
    :param numpy.ndarray z_reference:
        Local amplitude-reference position for each slab in Angstrom.
    :param numpy.ndarray n:
        Complex slab refractive indices, ambient first.
    :param numpy.ndarray kz:
        Renaud normal wavevectors in inverse Angstrom.
    :param numpy.ndarray A_plus:
        Downward-propagating slab amplitudes.
    :param numpy.ndarray A_minus:
        Upward-propagating slab amplitudes.
    :param complex r_S:
        Specular reflection amplitude in the incident medium.
    :param complex t_S:
        Transmission amplitude at the substrate boundary.
    :param str polarization:
        ``"s"`` or ``"p"``.
    """

    z: np.ndarray
    psi: np.ndarray
    z_interfaces: np.ndarray
    z_reference: np.ndarray
    n: np.ndarray
    kz: np.ndarray
    A_plus: np.ndarray
    A_minus: np.ndarray
    r_S: complex
    t_S: complex
    polarization: str


@dataclass(frozen=True)
class StratifiedProfile:
    """Optical media centers and their physical interface positions.

    :param numpy.ndarray values:
        ``(N, 3)`` array containing layer-center z, delta, and beta.
    :param numpy.ndarray boundaries:
        ``N - 1`` interfaces in Angstrom, ordered from substrate to ambient.
    """

    values: np.ndarray
    boundaries: np.ndarray

    @property
    def z(self):
        """Return layer-center z coordinates in Angstrom."""
        return self.values[:, 0]

    @property
    def delta(self):
        """Return layer delta values."""
        return self.values[:, 1]

    @property
    def beta(self):
        """Return layer beta values."""
        return self.values[:, 2]


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


def water_optical_profile(water_model, noUC=30, z_step=None, z_origin=None):
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
        Uniform z sampling interval in Angstrom. Defaults to the water-model
        lattice height when no atomistic profile supplies a spacing.
    :param float z_origin:
        Optional origin of the atomistic z grid in Angstrom. When supplied,
        the first water sample is aligned to that grid at or above the water
        onset.
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
    if z_step is None:
        z_step = water_model.a[2]
    if z_step <= 0.0:
        raise ValueError("z_step must be positive.")

    water_onset = water_model.pos_absolute + water_model.basis[0] * water_model.a[2]
    z_start = water_onset
    if z_origin is not None:
        grid_index = np.ceil((water_onset - z_origin) / z_step - 1e-12)
        z_start = z_origin + grid_index * z_step
    z_stop = water_onset + noUC * water_model.a[2]
    z = np.arange(z_start, z_stop + 0.5 * z_step, z_step, dtype=np.float64)
    rho_f = water_model.zDensity_G(z, 0.0, 0.0)
    wavelength = 12.398419843320026 / (water_model._E * 1e-3)
    scale = 2.8179403262e-5 * wavelength**2 / (2.0 * np.pi)
    profile = np.empty((z.size, 3), dtype=np.float64)
    profile[:, 0] = z
    profile[:, 1] = scale * rho_f.real
    profile[:, 2] = scale * rho_f.imag
    return np.ascontiguousarray(profile)


def top_layer_spacing(profile):
    """Return the spacing between the two highest optical-layer origins.

    :param numpy.ndarray profile:
        Atomistic optical profile with z in the first column.
    :returns:
        Positive top-layer spacing in Angstrom.
    :rtype: float
    :raises ValueError:
        If fewer than two distinct layer origins are available.
    """
    z = np.unique(np.asarray(profile[:, 0], dtype=np.float64))
    if z.size < 2:
        raise ValueError("At least two atomic layer origins are required for dz.")
    dz = float(z[-1] - z[-2])
    if dz <= 0.0:
        raise ValueError("Atomic layer spacing must be positive.")
    return dz


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


def profile_boundaries(profile):
    """Return interfaces centered between adjacent optical-profile samples.

    The profile coordinates identify the centers at which ``delta`` and
    ``beta`` apply.  Each interface is therefore the midpoint of two adjacent
    center coordinates.  The first and last media remain semi-infinite.

    :param numpy.ndarray profile:
        ``(N, 3)`` profile with strictly increasing layer-center z positions
        in Angstrom, delta, and beta.
    :returns:
        ``N - 1`` interfaces in Angstrom, ordered from substrate to ambient.
    :rtype: numpy.ndarray
    :raises ValueError:
        If the profile is not a valid ordered ``(N, 3)`` array.
    """
    profile = np.asarray(profile, dtype=np.float64)
    if profile.ndim != 2 or profile.shape[1] != 3 or len(profile) < 2:
        raise ValueError("profile must be an (N, 3) array with at least two rows.")
    if not np.all(np.diff(profile[:, 0]) > 0.0):
        raise ValueError("profile z positions must be strictly increasing.")
    return np.ascontiguousarray(0.5 * (profile[:-1, 0] + profile[1:, 0]))


def stratify_profile(profile, delta_tolerance=None, beta_tolerance=None):
    """Create finite optical layers and their midpoint interfaces.

    When simplification is requested, adjacent finite media are merged without
    moving their outer interfaces.  Their optical constants are weighted by
    the physical widths between those inherited interfaces.  The substrate
    and incident medium are never merged.

    :param numpy.ndarray profile:
        ``(N, 3)`` profile with layer-center z in Angstrom, delta, and beta.
    :param float delta_tolerance:
        Maximum delta range within one merged layer. Pass ``None`` to retain
        every sampled medium.
    :param float beta_tolerance:
        Maximum beta range within one merged layer. Defaults to
        ``delta_tolerance`` when simplification is requested.
    :returns:
        Optical media and their explicit physical boundaries.
    :rtype: StratifiedProfile
    :raises ValueError:
        If the profile or tolerances are invalid.
    """
    profile = np.asarray(profile, dtype=np.float64)
    boundaries = profile_boundaries(profile)
    if delta_tolerance is None:
        return StratifiedProfile(
            np.ascontiguousarray(profile.copy()), boundaries
        )
    if beta_tolerance is None:
        beta_tolerance = delta_tolerance
    if delta_tolerance < 0.0 or beta_tolerance < 0.0:
        raise ValueError("delta and beta tolerances must be non-negative.")
    if len(profile) <= 3:
        return StratifiedProfile(
            np.ascontiguousarray(profile.copy()), boundaries
        )

    values = [profile[0].copy()]
    inherited_boundaries = [boundaries[0]]
    start = 1
    last_finite = len(profile) - 2
    while start <= last_finite:
        stop = start
        delta_min = delta_max = profile[start, 1]
        beta_min = beta_max = profile[start, 2]
        while stop < last_finite:
            candidate = stop + 1
            next_delta_min = min(delta_min, profile[candidate, 1])
            next_delta_max = max(delta_max, profile[candidate, 1])
            next_beta_min = min(beta_min, profile[candidate, 2])
            next_beta_max = max(beta_max, profile[candidate, 2])
            if (
                next_delta_max - next_delta_min > delta_tolerance
                or next_beta_max - next_beta_min > beta_tolerance
            ):
                break
            stop = candidate
            delta_min, delta_max = next_delta_min, next_delta_max
            beta_min, beta_max = next_beta_min, next_beta_max

        lower_boundary = boundaries[start - 1]
        upper_boundary = boundaries[stop]
        layer_boundaries = boundaries[start - 1 : stop + 1]
        thicknesses = np.diff(layer_boundaries)
        merged = profile[start].copy()
        merged[0] = 0.5 * (lower_boundary + upper_boundary)
        merged[1:] = np.average(
            profile[start : stop + 1, 1:], axis=0, weights=thicknesses
        )
        values.append(merged)
        inherited_boundaries.append(upper_boundary)
        start = stop + 1

    values.append(profile[-1].copy())
    return StratifiedProfile(
        np.ascontiguousarray(np.asarray(values, dtype=np.float64)),
        np.ascontiguousarray(np.asarray(inherited_boundaries, dtype=np.float64)),
    )


def simplify_profile(profile, delta_tolerance=1e-9, beta_tolerance=None):
    """Merge adjacent optically similar finite layers.

    A merged layer receives thickness-weighted delta and beta, preserving the
    integrated optical constants through the combined finite thickness. The
    semi-infinite substrate and incident-medium rows are always retained.
    Grouping is bounded by the full delta/beta range in each group, preventing
    a long sequence of small changes from drifting beyond the tolerances.

    :param numpy.ndarray profile:
        ``(N, 3)`` profile with strictly increasing z, delta, and beta.
    :param float delta_tolerance:
        Maximum delta range within one merged layer.
    :param float beta_tolerance:
        Maximum beta range within one merged layer. Defaults to
        ``delta_tolerance``.
    :returns:
        Simplified C-contiguous optical profile.
    :rtype: numpy.ndarray
    :raises ValueError:
        If the profile or tolerances are invalid.
    """
    return stratify_profile(
        profile, delta_tolerance, beta_tolerance
    ).values


def _normal_wavevector(n, k0, alpha_rad):
    """Return normal propagation constants for incidence from medium zero."""
    k_parallel_over_k0 = n[0] * np.cos(alpha_rad)
    argument = n**2 - k_parallel_over_k0**2
    q = k0 * np.sqrt(argument)
    q = np.where(q.imag > 0.0, -q, q)
    q = np.where((q.imag == 0.0) & (q.real < 0.0), -q, q)
    # Avoid cancellation in n_0**2 - (n_0*cos(alpha))**2 at tiny angles.
    q[0] = k0 * n[0] * np.sin(alpha_rad)
    return q


def _admittance(n, q, polarization):
    """Return scalar optical admittance for s or p polarization."""
    if polarization == "s":
        return q
    if polarization == "p":
        return n**2 / q
    raise ValueError("polarization must be 's' or 'p'.")


def solve_wavefield(
    profile, energy_eV, alpha, polarization="s", boundaries=None
):
    """Solve the layered unperturbed wavefield by stable Parratt recursion.

    The profile rows define optical media at layer centers. By default, their
    interfaces are placed at adjacent-center midpoints. Explicit boundaries
    preserve the original outer interfaces of a simplified profile. The first
    and last rows represent the semi-infinite substrate and incident media.
    Angles are glancing angles measured from the surface inside the incident
    medium. The conserved tangential wavevector is therefore
    ``k_x = n_0 k_0 cos(alpha)``, where ``n_0`` is the last profile row.

    :param numpy.ndarray profile:
        ``(N, 3)`` optical profile with z in Angstrom, delta, and beta.
    :param float energy_eV:
        X-ray photon energy in eV.
    :param float alpha:
        Glancing angle in degrees inside the incident medium.
    :param str polarization:
        ``"s"`` or ``"p"``.
    :param numpy.ndarray boundaries:
        Optional ``N - 1`` interfaces in Angstrom, ordered from substrate to
        ambient. Defaults to adjacent-center midpoints.
    :returns:
        Slab amplitudes and sampled electric field.
    :rtype: Wavefield
    :raises ValueError:
        If the profile, energy, angle, or polarization is invalid.
    """
    profile = np.asarray(profile, dtype=np.float64)
    if profile.ndim != 2 or profile.shape[1] != 3 or len(profile) < 2:
        raise ValueError("profile must be an (N, 3) array with at least two rows.")
    default_boundaries = profile_boundaries(profile)
    if boundaries is None:
        boundaries = default_boundaries
    else:
        boundaries = np.asarray(boundaries, dtype=np.float64)
        if boundaries.ndim != 1 or len(boundaries) != len(profile) - 1:
            raise ValueError("boundaries must be a one-dimensional (N - 1) array.")
        if not np.all(np.diff(boundaries) > 0.0):
            raise ValueError("boundaries must be strictly increasing.")
        if not np.all(profile[:-1, 0] < boundaries) or not np.all(
            boundaries < profile[1:, 0]
        ):
            raise ValueError("each boundary must lie between adjacent layer centers.")
    if energy_eV <= 0.0:
        raise ValueError("energy_eV must be positive.")
    if alpha <= 0.0 or alpha > 90.0:
        raise ValueError("alpha must be in the interval (0, 90] degrees.")
    if polarization not in {"s", "p"}:
        raise ValueError("polarization must be 's' or 'p'.")

    wavelength = HC_KEV_ANGSTROM / (energy_eV * 1e-3)
    k0 = 2.0 * np.pi / wavelength
    alpha_rad = np.deg2rad(alpha)

    z_ascending = profile[:, 0]
    n_ascending = 1.0 - profile[:, 1] - 1j * profile[:, 2]
    n = np.ascontiguousarray(n_ascending[::-1])
    z_interfaces = np.ascontiguousarray(boundaries[::-1])
    q = _normal_wavevector(n, k0, alpha_rad)
    kz = -q
    admittance = _admittance(n, q, polarization)

    interface_r = (admittance[:-1] - admittance[1:]) / (
        admittance[:-1] + admittance[1:]
    )
    reflection = np.empty(len(interface_r), dtype=np.complex128)
    reflection[-1] = interface_r[-1]
    for interface in range(len(interface_r) - 2, -1, -1):
        thickness = z_interfaces[interface] - z_interfaces[interface + 1]
        phase = np.exp(-1j * q[interface + 1] * thickness)
        reflected_below = reflection[interface + 1] * phase**2
        reflection[interface] = (
            interface_r[interface] + reflected_below
        ) / (1.0 + interface_r[interface] * reflected_below)

    A_plus = np.empty(len(n), dtype=np.complex128)
    A_minus = np.empty(len(n), dtype=np.complex128)
    z_reference = np.empty(len(n), dtype=np.float64)
    A_plus[0] = 1.0
    A_minus[0] = reflection[0]
    z_reference[0] = z_interfaces[0]

    down_at_interface = A_plus[0]
    up_at_interface = A_minus[0]
    for interface in range(len(interface_r)):
        field = down_at_interface + up_at_interface
        derivative = admittance[interface] * (
            down_at_interface - up_at_interface
        )
        A_plus[interface + 1] = 0.5 * (
            field + derivative / admittance[interface + 1]
        )
        A_minus[interface + 1] = 0.5 * (
            field - derivative / admittance[interface + 1]
        )
        z_reference[interface + 1] = z_interfaces[interface]
        if interface + 1 < len(interface_r):
            thickness = z_interfaces[interface] - z_interfaces[interface + 1]
            phase = np.exp(-1j * q[interface + 1] * thickness)
            down_at_interface = A_plus[interface + 1] * phase
            up_at_interface = A_minus[interface + 1] / phase

    psi = np.empty(len(n), dtype=np.complex128)
    for ascending_index, z_value in enumerate(z_ascending):
        medium = len(n) - 1 - ascending_index
        depth = z_reference[medium] - z_value
        psi[ascending_index] = (
            A_plus[medium] * np.exp(-1j * q[medium] * depth)
            + A_minus[medium] * np.exp(1j * q[medium] * depth)
        )
    psi = np.ascontiguousarray(psi)

    return Wavefield(
        z=np.ascontiguousarray(z_ascending),
        psi=psi,
        z_interfaces=z_interfaces,
        z_reference=np.ascontiguousarray(z_reference),
        n=n,
        kz=np.ascontiguousarray(kz),
        A_plus=np.ascontiguousarray(A_plus),
        A_minus=np.ascontiguousarray(A_minus),
        r_S=complex(A_minus[0]),
        t_S=complex(A_plus[-1]),
        polarization=polarization,
    )

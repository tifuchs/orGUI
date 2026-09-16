"""Independent deterministic flat-height crystals for Poisson surface tests.

Shared by ``test_CTRcalc.py`` and ``test_CTRfilm.py``. The package has no
``conftest.py`` and its CTR tests are ``unittest`` classes, so these are plain
importable builders rather than pytest fixtures.

The oracle represents one flat surface height as a plain :class:`CTRfilm.Film`
of the corresponding thickness. It never calls ``PoissonSurface``, so a
placement or occupancy error in the Poisson assembly cannot cancel out of a
comparison against it.

Height convention, matching ``doc/design/incoherent_ctr_models.md``: a state is
named by the structural layer ``n`` of its top filled layer. Layers ``j < 0``
lie inside the sharp Film boundary and ``j >= 0`` above it, so state ``n``
corresponds to a Film of ``w_base + n + 1`` layers, and the signed height
change is ``H = n + 1``.

Units: lattice constants and heights in Angstrom, ``h``, ``k`` and ``l`` in
r.l.u. of the reference cell, amplitudes in electrons per reference lateral
cell.
"""

import copy

import numpy as np

from .. import CTRcalc, CTRfilm

__all__ = [
    "coherent_reference",
    "flat_height_crystal",
    "height_states",
    "layered_cell",
    "poisson_crystal",
]


def layered_cell(n_layers=2, name="layered", a=3.0, c=6.0):
    """Return a unit cell whose structural layers are evenly spaced in z.

    ``UnitCell.addAtom`` defaults ``layerpos`` to ``0.0`` for every new layer,
    which collapses the whole cycle onto one stacking position: every Film
    layer is then placed a full unit cell apart instead of a fraction of one,
    and any comparison against an independently stacked Film disagrees. The
    explicit ``layerpos`` assignment below is therefore load bearing, not
    decoration.

    :param int n_layers:
        Number of structural layers in the cycle.
    :param str name:
        Unit cell name.
    :param float a:
        In-plane lattice constant in Angstrom.
    :param float c:
        Out-of-plane lattice constant in Angstrom.
    :returns:
        Cell with ``n_layers`` layers at fractional heights ``i / n_layers``.
    :rtype: CTRcalc.UnitCell
    """
    unitcell = CTRcalc.UnitCell([a, a, c], [90.0, 90.0, 90.0], name=name)
    for index in range(n_layers):
        z = index / n_layers
        element = "C" if index % 2 == 0 else "O"
        unitcell.addAtom(element, [0.0, 0.0, z], 0.1, 0.1, 1.0, layer=index)
        unitcell.layerpos[float(index)] = z
    return unitcell


def _crystal(
    bulk_cell,
    components,
    stacking,
    reference_uc=None,
    weights=None,
    domains=None,
):
    """Assemble one crystal and apply any requested weights and domains.

    ``domains`` is applied to every top-level component, not to one of them.
    A flat-height crystal has a single component where the crystal under test
    has two, so applying it selectively would make the two sides differ by the
    transform itself rather than by the surface model.
    """
    keyargs = {"stacking": np.asarray(stacking)}
    if reference_uc is not None:
        keyargs["reference_uc"] = reference_uc
    crystal = CTRcalc.SXRDCrystal(bulk_cell, *components, **keyargs)
    if weights is not None:
        crystal.weights = np.asarray(weights, dtype=np.float64)
        crystal.weights_0 = np.copy(crystal.weights)
    if domains is not None:
        for index in range(len(components)):
            crystal.setDomain(index, list(domains))
    return crystal


def poisson_crystal(
    profile,
    w_base=4.0,
    n_layers=2,
    termination_cells=None,
    reference_uc=None,
    weights=None,
    domains=None,
):
    """Return the crystal under test and its Poisson surface.

    :param CTRdistributions.PoissonProfile profile:
        Height distribution driving the surface.
    :param float w_base:
        Base Film width in structural layers.
    :param int n_layers:
        Structural layers per unit cell.
    :param dict termination_cells:
        Optional explicit termination bank. Defaults to cells generated from
        the Film's own cell, which makes the termination replacement cancel.
    :param CTRcalc.UnitCell reference_uc:
        Optional reference cell, to exercise reference-area scaling.
    :param sequence weights:
        Optional component weights for ``[film, surface]``.
    :param sequence domains:
        Optional ``(matrix, occupancy)`` list applied to every component.
    :returns:
        The crystal and its surface component.
    :rtype: tuple[CTRcalc.SXRDCrystal, CTRfilm.PoissonSurface]
    """
    film = CTRfilm.Film(layered_cell(n_layers, "film"), name="film")
    film.basis[0] = w_base
    source = (
        termination_cells
        if termination_cells is not None
        else layered_cell(n_layers, "film")
    )
    surface = CTRfilm.PoissonSurface(source, profile=profile, name="surface")
    crystal = _crystal(
        layered_cell(n_layers, "bulk"),
        (film, surface),
        [1, 2],
        reference_uc=reference_uc,
        weights=weights,
        domains=domains,
    )
    return crystal, surface


def flat_height_crystal(
    n,
    w_base=4.0,
    n_layers=2,
    reference_uc=None,
    weights=None,
    domains=None,
):
    """Return the independent crystal whose surface is flat at top layer ``n``.

    Built only from a plain Film of ``w_base + n + 1`` layers, with no
    ``PoissonSurface`` anywhere.

    :param int n:
        Structural layer of the top filled layer.
    :param float w_base:
        Base Film width in structural layers.
    :param int n_layers:
        Structural layers per unit cell.
    :param CTRcalc.UnitCell reference_uc:
        Optional reference cell, matching ``poisson_crystal``.
    :param sequence weights:
        Optional single-element component weight for the Film.
    :param sequence domains:
        Optional ``(matrix, occupancy)`` list applied to the Film.
    :returns:
        Crystal with a sharp surface at the requested height.
    :rtype: CTRcalc.SXRDCrystal
    :raises ValueError:
        If the requested height etches away more than the base Film.
    """
    width = w_base + n + 1
    if width < 0:
        raise ValueError(
            f"flat height {n} etches below the base film of {w_base} layers"
        )
    film = CTRfilm.Film(layered_cell(n_layers, "film"), name="film")
    film.basis[0] = width
    return _crystal(
        layered_cell(n_layers, "bulk"),
        (film,),
        [1],
        reference_uc=reference_uc,
        weights=weights,
        domains=domains,
    )


def height_states(profile, minimum_probability=0.0):
    """Return the represented structural layers and their exposed fractions.

    These are the coherent path's ``q_n``: ``surface_occupancy`` masses, with
    the upper tail folded into the terminal bin. They are not the normalized
    ensemble probabilities of the incoherent model.

    :param CTRdistributions.PoissonProfile profile:
        Height distribution.
    :param float minimum_probability:
        Discard states at or below this exposed fraction.
    :returns:
        Structural layer numbers and their exposed fractions.
    :rtype: tuple[numpy.ndarray, numpy.ndarray]
    """
    lower, upper = profile.support()
    layers = np.arange(lower, upper + 1)
    exposed = profile.surface_occupancy(layers)
    keep = exposed > minimum_probability
    return layers[keep], exposed[keep]


def coherent_reference(profile, h, k, l, **keyargs):  # noqa: E741
    """Return ``sum_n q_n A_n`` over independent flat-height crystals.

    :param CTRdistributions.PoissonProfile profile:
        Height distribution.
    :param numpy.ndarray h:
        Reference-frame reciprocal coordinate in r.l.u.
    :param numpy.ndarray k:
        Reference-frame reciprocal coordinate in r.l.u.
    :param numpy.ndarray l:
        Reference-frame reciprocal coordinate in r.l.u.
    :param keyargs:
        Forwarded to :func:`flat_height_crystal`.
    :returns:
        Coherent amplitude average in electrons per reference lateral cell.
    :rtype: numpy.ndarray
    """
    layers, exposed = height_states(profile)
    total = np.zeros_like(np.asarray(l, dtype=np.float64), dtype=np.complex128)
    for layer, probability in zip(layers, exposed):
        crystal = flat_height_crystal(int(layer), **keyargs)
        total = total + probability * crystal.F(h, k, l)
    return total


def flat_height_amplitudes(profile, h, k, l, **keyargs):  # noqa: E741
    """Return the layers, exposed fractions, and per-state amplitudes.

    The squared-amplitude average ``sum_n q_n abs(A_n) ** 2`` built from these
    is the fully incoherent limit which later increments must reproduce.

    :returns:
        Structural layers, exposed fractions, and one amplitude array per
        state.
    :rtype: tuple[numpy.ndarray, numpy.ndarray, list[numpy.ndarray]]
    """
    layers, exposed = height_states(profile)
    amplitudes = [
        flat_height_crystal(int(layer), **keyargs).F(h, k, l)
        for layer in layers
    ]
    return layers, exposed, amplitudes


def minimum_base_width(profile, margin=2.0):
    """Return a base Film width which survives the deepest etched state.

    :param CTRdistributions.PoissonProfile profile:
        Height distribution.
    :param float margin:
        Extra layers kept below the deepest state.
    :returns:
        Base width in structural layers.
    :rtype: float
    """
    layers, _ = height_states(profile)
    return float(max(margin, -int(layers[0]) + margin))


def termination_bank(n_layers=2, scale=1.0, displacement=0.0):
    """Return a distinct termination bank for every layer of the cycle.

    :param int n_layers:
        Structural layers per unit cell.
    :param float scale:
        Occupancy multiplier applied to every termination atom.
    :param float displacement:
        Fractional z displacement applied to every termination atom.
    :returns:
        Mapping from structural layer to termination cell.
    :rtype: dict
    """
    source = layered_cell(n_layers, "termination")
    if scale != 1.0 or displacement != 0.0:
        source.basis[:, 6] = source.basis[:, 6] * scale
        source.basis[:, 3] = source.basis[:, 3] + displacement
    film_cell = layered_cell(n_layers, "film")
    return CTRfilm.generate_surface_termination_cells(
        source, np.asarray(list(film_cell.split_in_layers()))
    )


def deep_copy_cells(cells):
    """Return an independent copy of a termination bank."""
    return {layer: copy.deepcopy(cell) for layer, cell in cells.items()}

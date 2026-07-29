"""Generate the RuO2/TiO2 strain-coupling documentation figure."""

from __future__ import annotations

import copy
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from orgui.datautils.xrayutils import CTRcalc


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
MODEL_PATH = (
    REPOSITORY_ROOT
    / "examples"
    / "CTR"
    / "RuO2_TiO2_Poisson_etching.xtal"
)
OUTPUT_PATH = (
    REPOSITORY_ROOT
    / "doc"
    / "source"
    / "_static"
    / "epitaxy_strain_coupling_offset_density.png"
)

FILM_THICKNESS_ANGSTROM = 150.0
INTERFACE_WIDTH_ANGSTROM = 10.0
TAIL_PROBABILITY = 1e-4
OFFSETS = (0.0, 0.02, 0.17)
STRAIN_COUPLINGS = (0.0, 0.5, 1.0)
Z = np.linspace(-52.0, 60.0, 6500)

COLORS = {
    "total": "#3f454c",
    "bottom": "#2878b5",
    "top": "#e1812c",
    "displacement": "#6f4c9b",
}


def _prepare_template():
    template = CTRcalc.SXRDCrystal.fromFile(MODEL_PATH)
    interface = template["TiO2toRuO2"]
    film = template["RuO2"]
    width_cells = INTERFACE_WIDTH_ANGSTROM / interface.uc_bottom.a[2]
    layer_spacing = film.unitcell.a[2] / len(film.layers)
    film_layers = int(round(FILM_THICKNESS_ANGSTROM / layer_spacing))
    return template, width_cells, film_layers


def _make_model(
    template,
    width_cells,
    film_layers,
    strain_coupling,
    offset,
):
    bulk = copy.deepcopy(template["bulk"])
    interface = copy.deepcopy(template["TiO2toRuO2"])
    film = copy.deepcopy(template["RuO2"])

    interface.profile = CTRcalc.SkellamProfile(
        width=width_cells,
        asymmetry=0.0,
        tail_probability=TAIL_PROBABILITY,
    )
    interface.basis[:] = [width_cells, 0.0, strain_coupling, offset]
    interface.basis_0[:] = interface.basis

    film.basis[0] = film_layers
    film.basis_0[:] = film.basis

    model = CTRcalc.SXRDCrystal(
        bulk,
        interface,
        film,
        stacking=np.array([1, 2]),
        atten=template.atten,
    )
    model.apply_stacking()
    return model, bulk, interface, film


def _densities(bulk, interface, film):
    bottom = np.real(bulk.zDensity_G_asbulk(Z, 0, 0))
    bottom += sum(
        np.real(layer.zDensity_G(Z, 0, 0))
        for layer in interface.bottom_layers
    )
    top = sum(
        np.real(layer.zDensity_G(Z, 0, 0))
        for layer in interface.top_layers
    )
    top += np.real(film.zDensity_G(Z, 0, 0))
    return bottom, top, bottom + top


def _shared_displacement(reference, shifted, c_bulk):
    """Return lower-material addition positions and offset displacements."""
    positions = []
    displacements = []
    occupancies = []
    for base_layer, shifted_layer in zip(
        reference.bottom_layers,
        shifted.bottom_layers,
    ):
        domain_occupancy = np.asarray(
            base_layer.coherentDomainOccupancy,
            dtype=float,
        )
        for index in np.flatnonzero(domain_occupancy > 0.0):
            base_matrix = base_layer.coherentDomainMatrix[index]
            shifted_matrix = shifted_layer.coherentDomainMatrix[index]
            positions.append(base_matrix[2, 3] * c_bulk)
            displacements.append(
                (shifted_matrix[2, 3] - base_matrix[2, 3]) * c_bulk
            )
            occupancies.append(domain_occupancy[index])

    positions = np.asarray(positions)
    displacements = np.asarray(displacements)
    occupancies = np.asarray(occupancies)
    order = np.argsort(positions)
    return positions[order], displacements[order], occupancies[order]


def _annotate_offset(ax, delta):
    start_x = -28.0
    end_x = 28.0
    arrow_y = 3.02
    ax.annotate(
        "",
        xy=(end_x, arrow_y),
        xytext=(start_x, arrow_y),
        arrowprops={
            "arrowstyle": "->",
            "color": COLORS["displacement"],
            "linewidth": 1.7,
        },
        annotation_clip=False,
    )
    ax.text(
        (start_x + end_x) / 2.0,
        arrow_y + 0.15,
        rf"deep bulk $\longrightarrow$ film: "
        rf"$\Delta={delta:.2f}\,\mathrm{{\AA}}$",
        color=COLORS["displacement"],
        ha="center",
        va="bottom",
        fontsize=9,
    )


def _plot():
    template, width_cells, film_layers = _prepare_template()
    c_bulk = template["bulk"].a[2]
    data = {}
    interfaces = {}

    for offset in OFFSETS:
        for strain_coupling in STRAIN_COUPLINGS:
            _, bulk, interface, film = _make_model(
                template,
                width_cells,
                film_layers,
                strain_coupling,
                offset,
            )
            data[offset, strain_coupling] = _densities(
                bulk,
                interface,
                film,
            )
            interfaces[offset, strain_coupling] = interface

    density_max = max(
        np.max(total)
        for bottom, top, total in data.values()
    )
    density_scale = 0.82 / density_max
    row_baselines = dict(zip(STRAIN_COUPLINGS, (0.0, 1.0, 2.0)))

    figure = plt.figure(figsize=(18.0, 8.2), layout="constrained")
    grid = figure.add_gridspec(2, 3, height_ratios=(4.4, 1.15))
    density_axes = []
    displacement_axes = []

    for column, offset in enumerate(OFFSETS):
        density_ax = figure.add_subplot(grid[0, column])
        displacement_ax = figure.add_subplot(
            grid[1, column],
            sharex=density_ax,
        )
        density_axes.append(density_ax)
        displacement_axes.append(displacement_ax)

        for strain_coupling in STRAIN_COUPLINGS:
            bottom, top, total = data[offset, strain_coupling]
            baseline = row_baselines[strain_coupling]
            density_ax.plot(
                Z,
                baseline + density_scale * total,
                color=COLORS["total"],
                linewidth=1.45,
                label=(
                    "total"
                    if column == 0 and strain_coupling == 0.0
                    else None
                ),
                zorder=3,
            )
            density_ax.plot(
                Z,
                baseline + density_scale * bottom,
                color=COLORS["bottom"],
                linewidth=1.0,
                alpha=0.92,
                label=(
                    r"TiO$_2$ contribution"
                    if column == 0 and strain_coupling == 0.0
                    else None
                ),
                zorder=2,
            )
            density_ax.plot(
                Z,
                baseline + density_scale * top,
                color=COLORS["top"],
                linewidth=1.0,
                alpha=0.92,
                label=(
                    r"RuO$_2$ contribution"
                    if column == 0 and strain_coupling == 0.0
                    else None
                ),
                zorder=2,
            )
            density_ax.text(
                -50.0,
                baseline + 0.09,
                rf"$\kappa={strain_coupling:g}$",
                ha="left",
                va="bottom",
                fontsize=10,
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.88,
                },
                zorder=5,
            )

        delta = offset * c_bulk
        _annotate_offset(density_ax, delta)
        density_ax.axvline(
            0.0,
            color="#777777",
            linestyle=":",
            linewidth=0.9,
            zorder=1,
        )
        density_ax.set_xlim(Z[0], Z[-1])
        density_ax.set_ylim(-0.08, 3.38)
        density_ax.set_yticks([0.0, 1.0, 2.0])
        density_ax.set_yticklabels([])
        density_ax.set_title(
            rf"$o={offset:g}$"
            + "\n"
            + rf"$\Delta=o\,c_{{\mathrm{{TiO_2}}}}={delta:.2f}\,$Å",
            fontsize=12,
        )
        density_ax.grid(axis="x", alpha=0.15)
        density_ax.tick_params(axis="x", labelbottom=False)

        for strain_coupling in STRAIN_COUPLINGS:
            reference = interfaces[0.0, strain_coupling]
            shifted = interfaces[offset, strain_coupling]
            positions, displacement, occupancy = _shared_displacement(
                reference,
                shifted,
                c_bulk,
            )
            visible = occupancy >= TAIL_PROBABILITY
            displacement_ax.plot(
                positions[visible],
                displacement[visible],
                marker="o",
                markersize=2.7,
                linewidth=1.2,
                color=COLORS["displacement"],
                alpha=0.35 + 0.6 * strain_coupling,
                label=(
                    rf"$\kappa={strain_coupling:g}$"
                    if column == 2
                    else None
                ),
            )
            plateau = strain_coupling * delta
            displacement_ax.text(
                58.0,
                plateau + 0.025,
                rf"${plateau:.2f}$",
                color=COLORS["displacement"],
                ha="right",
                va="bottom",
                fontsize=8.5,
            )

        displacement_ax.axhline(0.0, color="#777777", linewidth=0.7)
        displacement_ax.axvline(
            0.0,
            color="#777777",
            linestyle=":",
            linewidth=0.9,
        )
        displacement_ax.set_ylim(-0.06, 1.24)
        displacement_ax.set_xlabel(r"$z$ / Å")
        displacement_ax.grid(alpha=0.15)

    density_axes[0].set_ylabel(
        r"Vertically separated $\rho_{00}(z)$"
        "\n"
        r"(common density scale)"
    )
    displacement_axes[0].set_ylabel(
        r"Shared offset"
        "\n"
        r"$u_C(z)$ / Å"
    )
    density_axes[0].legend(
        loc="upper right",
        frameon=False,
        fontsize=9,
    )
    displacement_axes[-1].legend(
        loc="upper left",
        frameon=False,
        fontsize=8.5,
    )
    figure.suptitle(
        r"RuO$_2$/TiO$_2$: density and occupancy-mediated "
        r"interface expansion"
        "\n"
        r"15 nm RuO$_2$ film; 1 nm interface width",
        fontsize=15,
    )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PATH, dpi=180, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    _plot()

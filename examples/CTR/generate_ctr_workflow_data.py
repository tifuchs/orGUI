"""Generate the synthetic measurements used by the CTR workflow notebook.

The output is deterministic and uses the metadata-aware, ANAROD-like orGUI
CTR text format.  Run this script after installing the current orGUI checkout::

    python examples/CTR/generate_ctr_workflow_data.py
"""

from pathlib import Path

import numpy as np

from orgui.datautils.xrayutils import CTRcalc, CTRopt, CTRplotutil, CTRuc, unitcells


ENERGY_EV = 20_000.0
PT_BULK_DW = 0.4353
PT_SURFACE_EXPANSION = 0.02
PT_SURFACE_DW_FACTOR = 1.4
PT_AMPLITUDE_SCALE = 3.7
PT_MEAN_COUNTS = 200_000.0
REFLECTIVITY_INCIDENT_COUNTS = 1.0e8
RANDOM_SEED = 20260909


def build_pt100_truth():
    """Return the Pt(100) model used to generate the ``(2 0 L)`` rod.

    The Pt(100) planes in the conventional bundled cell are separated by
    ``c / 2``.  Moving the explicit outer plane by ``0.01 c`` therefore
    expands the final interplanar spacing by two percent.  Both in-plane and
    out-of-plane surface Debye--Waller values are 1.4 times the bulk value.

    :returns: Pt(100) bulk with one explicit, relaxed surface plane.
    :rtype: CTRcalc.SXRDCrystal
    """
    bulk = unitcells.unitcell("Pt100")
    surface = CTRuc.UnitCell(
        bulk.a, np.rad2deg(bulk.alpha), name="Pt_surface"
    )
    z_shift = 0.5 * PT_SURFACE_EXPANSION
    surface_dw = PT_BULK_DW * PT_SURFACE_DW_FACTOR
    surface.addAtom("Pt", [0.0, 0.0, z_shift], surface_dw, surface_dw, 1.0)
    surface.addAtom("Pt", [0.5, 0.5, z_shift], surface_dw, surface_dw, 1.0)
    crystal = CTRcalc.SXRDCrystal(bulk, surface)
    crystal.setEnergy(ENERGY_EV)
    return crystal


def generate_pt100_measurement(rng):
    """Return noisy, arbitrarily scaled Pt(100) structure-factor data.

    The conversion from model intensity to counts is chosen to give 200,000
    counts per point on average. The square-root amplitudes and propagated
    uncertainties are then both multiplied by ``PT_AMPLITUDE_SCALE`` to
    emulate an unscaled experiment.

    :param numpy.random.Generator rng: Reproducible random-number generator.
    :returns: Numerical ``H K L F errorF P`` table.
    :rtype: numpy.ndarray
    """
    ell = np.linspace(0.02, 6.0, 600)
    crystal = build_pt100_truth()
    amplitude = np.abs(
        crystal.F(np.full_like(ell, 2.0), np.zeros_like(ell), ell)
    )
    intensity = amplitude**2
    count_scale = PT_MEAN_COUNTS / np.mean(intensity)
    counts = rng.poisson(count_scale * intensity)
    measured = PT_AMPLITUDE_SCALE * np.sqrt(counts / count_scale)
    uncertainty = np.full_like(
        measured, PT_AMPLITUDE_SCALE / (2.0 * np.sqrt(count_scale))
    )
    return np.column_stack(
        (
            np.full_like(ell, 2.0),
            np.zeros_like(ell),
            ell,
            measured,
            uncertainty,
            np.ones_like(ell),
        )
    )


def generate_ruo2_tio2_measurement(rng, model_path):
    """Return noisy absolute RuO2/TiO2 low-L reflectivity data.

    The source model has a fully Poissonian etched surface with ``W = -6``
    structural layers and ``alpha = 1``.  Counts are sampled from the absolute
    DWBA reflectivity and divided by the incident count, so no arbitrary scale
    is introduced.

    :param numpy.random.Generator rng: Reproducible random-number generator.
    :param path-like model_path: RuO2/TiO2 crystal model file.
    :returns: Numerical ``H K L R errorR`` table.
    :rtype: numpy.ndarray
    """
    crystal = CTRcalc.SXRDCrystal.fromFile(model_path)
    ell = np.linspace(0.015, 0.595, 240)
    reduction = CTRplotutil.MeasurementReduction(
        "reflectivity", CTRplotutil.PolarizationReduction(1.0, "s")
    )
    rod = CTRplotutil.CTR(
        (0, 0),
        ell,
        np.ones_like(ell),
        np.ones_like(ell),
        reduction=reduction,
        scan_geometry=CTRplotutil.CTRScanGeometry("eq"),
    )
    optimizer = CTRopt.CTROptimizer(
        crystal,
        CTRplotutil.CTRCollection([rod]),
        scale_policy={"R": "fixed"},
    )
    optimizer.set_dwba()
    optimizer.prepareFit()
    absolute_reflectivity = optimizer.flat_prediction()
    counts = rng.poisson(REFLECTIVITY_INCIDENT_COUNTS * absolute_reflectivity)
    measured = counts / REFLECTIVITY_INCIDENT_COUNTS
    uncertainty = np.sqrt(np.maximum(counts, 1.0)) / REFLECTIVITY_INCIDENT_COUNTS
    return np.column_stack(
        (
            np.zeros_like(ell),
            np.zeros_like(ell),
            ell,
            measured,
            uncertainty,
        )
    )


def write_ctr_file(path, table, *, quantity, columns, outgoing="s"):
    """Write one schema-2 orGUI CTR table.

    :param path-like path: Destination filename.
    :param numpy.ndarray table: Point-aligned numerical data.
    :param str quantity: ``"structure_factor"`` or ``"reflectivity"``.
    :param str columns: Whitespace-separated column names.
    :param str outgoing: Analysed outgoing polarization channel.
    """
    header = "\n".join(
        (
            "orgui_ctr_schema: 2",
            f"quantity: {quantity}",
            "s_fraction: 1.0",
            f"outgoing: {outgoing}",
            f"columns: {columns}",
        )
    )
    np.savetxt(path, table, header=header, fmt="%.10g")


def main():
    """Regenerate both tracked CTR workflow measurement files."""
    repository_root = Path(__file__).resolve().parents[2]
    output_directory = repository_root / "doc" / "source" / "_static"
    rng = np.random.default_rng(RANDOM_SEED)
    write_ctr_file(
        output_directory / "ctr_pt100_20l.ctr",
        generate_pt100_measurement(rng),
        quantity="structure_factor",
        columns="H K L F_HKL errorF polarization_factor",
    )
    write_ctr_file(
        output_directory / "ctr_ruo2_tio2_00l.ctr",
        generate_ruo2_tio2_measurement(
            rng, Path(__file__).with_name("RuO2_TiO2_Poisson_etching.xtal")
        ),
        quantity="reflectivity",
        columns="H K L R errorR",
    )


if __name__ == "__main__":
    main()

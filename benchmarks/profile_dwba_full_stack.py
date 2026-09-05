"""Profile the atomistic DWBA kernel on a complete layered stack.

This non-gating diagnostic constructs, from bottom to top, a semi-infinite
bulk, a graded epitaxy interface, a finite Film, and a PoissonSurface.  It
reports public preparation/evaluation timings and isolates the native finite,
bulk, form-factor, and specular-reference work by timing controlled variants
of the same packed calculation.

The private calls in this file are intentional: this is a profiler for the
implementation rather than a public API example.
"""

import argparse
import cProfile
import io
from pathlib import Path
import pstats
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from orgui.datautils.xrayutils import CTRuc
from orgui.datautils.xrayutils.CTRcalc import SXRDCrystal
from orgui.datautils.xrayutils.CTRdistributions import (
    PoissonProfile,
    SkellamProfile,
)
from orgui.datautils.xrayutils.CTRfilm import (
    EpitaxyInterface,
    Film,
    PoissonSurface,
)
from orgui.datautils.xrayutils.CTRuc import UnitCell


_ENERGY_EV = 12_000.0
_LATTICE = [3.8, 3.8, 5.6]
_ANGLES = [90.0, 90.0, 90.0]


def _cell(name, atoms, seed, elements):
    rng = np.random.default_rng(seed)
    cell = UnitCell(_LATTICE, _ANGLES, name=name)
    for index, position in enumerate(rng.random((atoms, 3))):
        cell.addAtom(
            elements[index % len(elements)],
            position,
            0.15,
            0.22,
            1.0,
            layer=0,
        )
    cell.setEnergy(_ENERGY_EV)
    return cell


def _full_stack(atoms, film_layers, interface_width, surface_mean):
    bulk = _cell("bulk", atoms, 101, ("Ti", "O"))
    interface_top = _cell("interface_top", atoms, 102, ("Ru", "O"))
    interface_bottom = _cell("interface_bottom", atoms, 103, ("Ti", "O"))
    film_cell = _cell("film_cell", atoms, 104, ("Ru", "O"))
    termination = _cell("termination", atoms, 105, ("Ru", "O"))

    interface = EpitaxyInterface(
        interface_top,
        interface_bottom,
        profile=SkellamProfile(interface_width),
        name="graded_interface",
    )
    film = Film(film_cell, name="film")
    film.basis[0] = float(film_layers)
    surface = PoissonSurface(
        termination,
        profile=PoissonProfile(surface_mean, 0.75),
        name="poisson_surface",
    )
    crystal = SXRDCrystal(
        bulk,
        interface,
        film,
        surface,
        stacking=np.array([0, 1, 2]),
        name="full_dwba_profile_stack",
    )
    crystal.setEnergy(_ENERGY_EV)
    return crystal


def _field_array(prepared, value):
    media = len(prepared.reference.n)
    return np.ascontiguousarray(np.asarray(value).reshape(media, -1))


def _call_native(
    prepared,
    packed,
    *,
    native_module=None,
    specular=None,
    reference_shares=None,
    zero_reflected=False,
):
    if native_module is None:
        native_module = CTRuc._CTRcalc_cpp
    if specular is None:
        specular = prepared.is_specular
    if reference_shares is None:
        reference_shares = prepared._atomic_model.reference_shares

    incident_minus = _field_array(prepared, prepared.field_i.A_minus)
    incident_minus_over_kz = _field_array(
        prepared, prepared.field_i._A_minus_over_kz
    )
    exit_minus = _field_array(prepared, prepared.field_f.A_minus)
    exit_minus_over_kz = _field_array(
        prepared, prepared.field_f._A_minus_over_kz
    )
    if zero_reflected:
        incident_minus = np.zeros_like(incident_minus)
        incident_minus_over_kz = np.zeros_like(incident_minus_over_kz)
        exit_minus = np.zeros_like(exit_minus)
        exit_minus_over_kz = np.zeros_like(exit_minus_over_kz)

    return native_module.unitcell_F_DWBA_records(
        prepared.Q_parallel,
        specular,
        prepared.reference.n,
        prepared.field_i_index,
        prepared.field_f_index,
        _field_array(prepared, prepared.field_i.kz),
        _field_array(prepared, prepared.field_i.A_plus),
        incident_minus,
        _field_array(prepared, prepared.field_i._A_plus_over_kz),
        incident_minus_over_kz,
        _field_array(prepared, prepared.field_f.kz),
        _field_array(prepared, prepared.field_f.A_plus),
        exit_minus,
        _field_array(prepared, prepared.field_f._A_plus_over_kz),
        exit_minus_over_kz,
        prepared.field_i.z_reference,
        prepared.field_f.z_reference,
        prepared.field_i.z_interfaces,
        prepared.alpha_i,
        prepared.alpha_f,
        prepared.cos_azimuth,
        prepared.sin_azimuth,
        prepared.k0,
        prepared.polarization_i,
        prepared.polarization_f,
        reference_shares,
        prepared.reference_area,
        True,
        0.0,
        packed["atom_basis"],
        packed["form_factors"],
        packed["finite_positions"],
        packed["finite_atom_index"],
        packed["finite_record_index"],
        packed["finite_weight"],
        packed["bulk_positions"],
        packed["bulk_atom_index"],
        packed["bulk_weight"],
        packed["bulk_repeat_coordinate"],
        packed["bulk_repeat"],
    )


def _without_instances(packed, finite=False, bulk=False, atoms=False):
    result = dict(packed)
    if finite:
        result.update(
            finite_positions=np.empty((0, 3), dtype=np.float64),
            finite_atom_index=np.empty(0, dtype=np.int64),
            finite_record_index=np.empty(0, dtype=np.int64),
            finite_weight=np.empty(0, dtype=np.float64),
        )
    if bulk:
        result.update(
            bulk_positions=np.empty((0, 3), dtype=np.float64),
            bulk_atom_index=np.empty(0, dtype=np.int64),
            bulk_weight=np.empty(0, dtype=np.float64),
            bulk_repeat_coordinate=np.empty(0, dtype=np.float64),
        )
    if atoms:
        result.update(
            atom_basis=np.empty((0, 8), dtype=np.float64),
            form_factors=np.empty((0, 13), dtype=np.float64),
        )
    return result


def _interleaved_median_times(functions, repetitions, trials):
    """Time functions in shuffled trial order to limit thermal-order bias."""
    for function in functions.values():
        function()
    rng = np.random.default_rng(20260903)
    samples = {label: [] for label in functions}
    labels = list(functions)
    for _ in range(trials):
        rng.shuffle(labels)
        for label in labels:
            start = time.perf_counter()
            for _ in range(repetitions):
                functions[label]()
            samples[label].append(
                (time.perf_counter() - start) / repetitions
            )
    return {
        label: float(np.median(values)) for label, values in samples.items()
    }


def _profile_public_evaluation(crystal, prepared, repetitions):
    profiler = cProfile.Profile()
    profiler.enable()
    for _ in range(repetitions):
        crystal.dwba.evaluate_prepared(prepared)
    profiler.disable()
    stream = io.StringIO()
    stats = pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats("cumulative")
    stats.print_stats(18)
    return stream.getvalue()


def _prepare_rod(crystal, points, specular):
    L = np.linspace(0.15, 3.0, points)
    if specular:
        crystal.dwba.set_ctr_geometry(equal_angles=True)
        return crystal.dwba.prepare(0.0, 0.0, L)
    crystal.dwba.set_ctr_geometry(alpha_i=np.deg2rad(0.2), rods=[(1.0, 0.0)])
    return crystal.dwba.prepare(1.0, 0.0, L)


def _main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--points", type=int, default=1024)
    parser.add_argument("--atoms", type=int, default=32)
    parser.add_argument("--film-layers", type=int, default=16)
    parser.add_argument("--interface-width", type=float, default=1.25)
    parser.add_argument("--surface-mean", type=float, default=2.5)
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--trials", type=int, default=5)
    arguments = parser.parse_args()

    crystal = _full_stack(
        arguments.atoms,
        arguments.film_layers,
        arguments.interface_width,
        arguments.surface_mean,
    )
    start = time.perf_counter()
    prepared = _prepare_rod(crystal, arguments.points, specular=False)
    prepare_elapsed = time.perf_counter() - start
    records = crystal.dwba._validate_prepared(prepared)
    packed = crystal.dwba._pack_atomic_records(prepared, records)
    terminal_boundary = prepared.field_i.z_interfaces[-1]
    maximum_bulk_z = np.max(packed["bulk_positions"][:, 2])
    terminal_repeat = max(
        0,
        int(
            np.floor(
                (maximum_bulk_z - terminal_boundary) / packed["bulk_repeat"][2]
            )
            + 1
        ),
    )

    no_finite = _without_instances(packed, finite=True)
    no_bulk = _without_instances(packed, bulk=True)
    no_scatter = _without_instances(packed, finite=True, bulk=True)
    geometry_only = _without_instances(
        packed, finite=True, bulk=True, atoms=True
    )
    reference_zero = np.zeros_like(prepared._atomic_model.reference_shares)

    repetitions = arguments.repetitions
    trials = arguments.trials
    functions = {
        "public evaluate_prepared": lambda: crystal.dwba.evaluate_prepared(
            prepared
        ),
        "validate prepared": lambda: crystal.dwba._validate_prepared(prepared),
        "pack atomic records": lambda: crystal.dwba._pack_atomic_records(
            prepared, records
        ),
        "native full": lambda: _call_native(
            prepared, packed, reference_shares=reference_zero
        ),
        "native without finite instances": lambda: _call_native(
            prepared, no_finite, reference_shares=reference_zero
        ),
        "native without bulk instances": lambda: _call_native(
            prepared, no_bulk, reference_shares=reference_zero
        ),
        "native without scattering instances": lambda: _call_native(
            prepared, no_scatter, reference_shares=reference_zero
        ),
        "native geometry only": lambda: _call_native(
            prepared, geometry_only, reference_shares=reference_zero
        ),
        "native incident/exit ++ only": lambda: _call_native(
            prepared,
            packed,
            reference_shares=reference_zero,
            zero_reflected=True,
        ),
    }
    timings = _interleaved_median_times(functions, repetitions, trials)

    specular = _prepare_rod(crystal, arguments.points, specular=True)
    specular_records = crystal.dwba._validate_prepared(specular)
    specular_packed = crystal.dwba._pack_atomic_records(specular, specular_records)
    specular_no_scatter = _without_instances(
        specular_packed, finite=True, bulk=True
    )
    specular_timings = _interleaved_median_times(
        {
            "reference": lambda: _call_native(specular, specular_no_scatter),
            "no reference": lambda: _call_native(
                specular,
                specular_no_scatter,
                specular=np.zeros_like(specular.is_specular),
            ),
        },
        repetitions,
        trials,
    )
    specular_reference = specular_timings["reference"]
    specular_without_reference = specular_timings["no reference"]

    print("Full-stack atomistic DWBA profile")
    print(
        f"  points={arguments.points}, atoms/cell={arguments.atoms}, "
        f"records={len(records)}, media={len(prepared.reference.n)}"
    )
    print(
        f"  packed atoms={len(packed['atom_basis'])}, "
        f"finite instances={len(packed['finite_positions'])}, "
        f"bulk instances={len(packed['bulk_positions'])}"
    )
    finite_terms = arguments.points * 4 * len(packed["finite_positions"])
    bulk_terms = (
        arguments.points
        * 4
        * len(packed["bulk_positions"])
        * terminal_repeat
    )
    print(
        f"  finite channel terms={finite_terms:,}, "
        f"bulk explicit channel terms={bulk_terms:,}, "
        f"terminal repeats={terminal_repeat}"
    )
    print(f"  prepare non-specular rod: {prepare_elapsed:.6f} s")
    for label, elapsed in timings.items():
        print(f"  {label:37s} {elapsed:.6f} s")

    native = timings["native full"]
    finite_cost = native - timings["native without finite instances"]
    bulk_cost = native - timings["native without bulk instances"]
    factor_cost = (
        timings["native without scattering instances"]
        - timings["native geometry only"]
    )
    reflected_cost = native - timings["native incident/exit ++ only"]
    print("Isolated native costs (differences of matched calls)")
    print(
        f"  finite four-channel accumulation: {finite_cost:.6f} s "
        f"({100.0 * finite_cost / native:.1f}% of native)"
    )
    print(
        f"  bulk explicit/tail accumulation:  {bulk_cost:.6f} s "
        f"({100.0 * bulk_cost / native:.1f}% of native)"
    )
    print(
        f"  form-factor grouping/tables:       {factor_cost:.6f} s "
        f"({100.0 * factor_cost / native:.1f}% of native)"
    )
    print(
        f"  three reflected-channel increment: {reflected_cost:.6f} s "
        f"({100.0 * reflected_cost / native:.1f}% of native)"
    )
    print(
        "  specular reference integral:       "
        f"{specular_reference - specular_without_reference:.6f} s"
    )
    print("cProfile: repeated public non-specular evaluation")
    print(_profile_public_evaluation(crystal, prepared, repetitions))


if __name__ == "__main__":
    _main()

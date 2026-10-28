"""Physical and API regressions for separated multilayer epitaxy stacks."""

import copy
import json
import pickle

import numpy as np
import pytest

from .. import CTRcalc, CTRfilm
from ..CTRdistributions import SkellamProfile
from ..CTRstacking import StackValidationError


def _cell(name, element="C", repeat=6.0, positions=(0., .5),
          shift=(0., 0., 1.), lateral=3.):
    cell = CTRcalc.UnitCell([lateral, lateral, repeat], [90.] * 3, name=name)
    for layer, z in enumerate(positions, 1):
        cell.addAtom(element, [0., 0., z], .1, .1, 1., layer=layer)
        cell.layerpos[float(layer)] = z
    cell.stacking_translation = shift
    return cell


def _parts(widths=(40, 40), roughness=.5, tail=1e-10,
           repeats=(6., 6., 6., 6.), positions=(0., .5), shift=(0., 0., 1.),
           lateral=(3., 3., 3., 3.)):
    bulk = _cell("bulk", repeat=repeats[0], positions=positions, shift=shift,
                 lateral=lateral[0])
    below = bulk
    components = []
    for i, width in enumerate(widths):
        above = _cell(f"material_{i}", ("Si", "Ge", "Al")[i % 3],
                      repeats[(i + 1) % len(repeats)], positions, shift,
                      lateral[(i + 1) % len(lateral)])
        interface = CTRfilm.EpitaxyInterface(
            copy.deepcopy(above), copy.deepcopy(below),
            name=f"interface_{i}",
            profile=SkellamProfile(roughness, tail_probability=tail),
        )
        interface.basis_0[:] = interface.basis
        film = CTRfilm.Film(copy.deepcopy(above), name=f"film_{i}")
        film.basis[:] = film.basis_0[:] = width
        components.extend((interface, film))
        below = above
    return bulk, components


def _crystal(*, strict=False, **kwargs):
    bulk, components = _parts(**kwargs)
    if strict:
        crystal = CTRcalc.SXRDCrystal.from_stack(bulk, components)
    else:
        crystal = CTRcalc.SXRDCrystal(
            bulk, *components, stacking=np.arange(1, len(components) + 1),
        )
    crystal.setEnergy(18000.)
    return crystal


def _hkl():
    ell = np.linspace(.13, 2.87, 30)
    return np.ones_like(ell), np.zeros_like(ell), ell


@pytest.mark.parametrize("strict", [False, True])
@pytest.mark.parametrize("lateral", [(3., 3., 3., 3.), (3., 4., 5., 6.)])
@pytest.mark.parametrize("widths", [(40, 40), (39, 41), (40, 41, 42)])
def test_sharp_stack_matches_explicit_finite_sites(widths, strict, lateral):
    """Finite material sites start above native bulk repeat zero at z=6 A."""
    crystal = _crystal(widths=widths, roughness=0., strict=strict, lateral=lateral)
    h, k, ell = _hkl()
    expected = crystal.uc_bulk.F_bulk(h, k, ell, crystal.atten)
    height = 6.
    for i, count in enumerate(widths):
        side = lateral[i + 1]
        atom = CTRcalc.UnitCell([side, side, 6.], [90.] * 3, name="atom")
        atom.addAtom(("Si", "Ge", "Al")[i % 3], [0., 0., 0.], .1, .1, 1.)
        atom.setReferenceUnitCell(crystal.uc_bulk)
        atom.setEnergy(18000.)
        sites = height + np.arange(count) * 3.
        expected += (3. / side)**2 * atom.F_uc(h, k, ell) * np.sum(
            np.exp(2j * np.pi * ell[:, None] * sites / 6.), axis=1,
        )
        height += count * 3.
    np.testing.assert_allclose(crystal.F(h, k, ell), expected, atol=2e-10)


def test_separated_stack_tail_convergence():
    """Extending separated transition tails leaves the physical model stable."""
    coarse = _crystal(tail=1e-6)
    fine = _crystal(tail=1e-10)
    np.testing.assert_allclose(coarse.F(*_hkl()), fine.F(*_hkl()),
                               rtol=0, atol=1e-6)
    assert coarse["film_1"].height_absolute == pytest.approx(
        fine["film_1"].height_absolute)


@pytest.mark.parametrize("width", [0., 2.5, 40.])
def test_width_keyword_initializes_fit_basis(width):
    """The convenience keyword preserves nominal fractional layer widths."""
    film = CTRfilm.Film(_cell("cell"), width_layers=width)
    np.testing.assert_array_equal(film.basis, [width])
    np.testing.assert_array_equal(film.basis_0, [width])


@pytest.mark.parametrize("width", [-1., np.nan, np.inf])
def test_width_keyword_rejects_invalid_values(width):
    """Invalid initial widths fail before creating a partial Film."""
    with pytest.raises(ValueError, match="width_layers"):
        CTRfilm.Film(_cell("cell"), width_layers=width)


@pytest.mark.parametrize("case", ["empty", "missing", "reverse", "names",
                                  "reserved", "reuse", "surface"])
def test_invalid_topology(case):
    """Strict construction rejects ambiguous or unsupported topology."""
    bulk, components = _parts()
    if case == "empty":
        components = []
    elif case == "missing":
        components.pop()
    elif case == "reverse":
        components.reverse()
    elif case == "names":
        components[-1].name = components[0].name
    elif case == "reserved":
        components[0].name = "bulk"
    elif case == "reuse":
        components[2:] = components[:2]
    else:
        components.append(CTRfilm.PoissonSurface(_cell("surface")))
    with pytest.raises(StackValidationError):
        CTRcalc.SXRDCrystal.from_stack(bulk, components)


@pytest.mark.parametrize("kwargs", [{"overlap": "allow"},
                                    {"stacking": np.arange(4)},
                                    {"enable_stacking": False}])
def test_factory_rejects_conflicting_options(kwargs):
    """Strict semantics cannot be bypassed through constructor overrides."""
    with pytest.raises(StackValidationError):
        CTRcalc.SXRDCrystal.from_stack(*_parts(), **kwargs)


@pytest.mark.parametrize("width", [12, 8, 6])
def test_strict_rejects_overlapping_or_consumed_film(width):
    """Both overlap and insufficient lower support are named errors."""
    with pytest.raises(StackValidationError, match="interface|film"):
        _crystal(strict=True, widths=(width, 40))


@pytest.mark.parametrize("observable", ["F", "F2", "F_surf", "density", "optics"])
def test_all_observables_revalidate_and_recover(observable):
    """Invalid layouts cannot be published and do not poison later layouts."""
    crystal = _crystal(strict=True)

    def evaluate():
        if observable == "density":
            return crystal.zDensity_G(np.linspace(-40, 280, 501), 0., 0.)
        if observable == "optics":
            return crystal.optical_profile()
        return getattr(crystal, observable)(*_hkl())

    expected = evaluate()
    crystal["film_0"].basis[0] = 12
    with pytest.raises(StackValidationError, match="overlaps"):
        evaluate()
    crystal["film_0"].basis[0] = 40
    np.testing.assert_allclose(evaluate(), expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("case", ["weights", "domains", "levels", "disabled",
                                  "cell", "negative", "template_domain"])
def test_live_definition_mutations_are_checked(case):
    """Mutable legacy attributes cannot silently weaken a strict crystal."""
    crystal = _crystal(strict=True)
    if case == "weights":
        crystal.weights[0] = .5
    elif case == "domains":
        crystal.domains[0] = [(np.eye(3), .5)]
    elif case == "levels":
        crystal.uc_stacking[1] = crystal.uc_stacking[0]
    elif case == "disabled":
        crystal.enable_uc_stacking = False
    elif case == "cell":
        crystal["film_0"].unitcell.stacking_translation = (.2, 0., 1.)
    elif case == "template_domain":
        crystal["film_0"].unitcell.coherentDomainMatrix[0][2, 2] = 1.1
    else:
        crystal["film_0"].basis[0] = -1
    with pytest.raises(StackValidationError):
        crystal.F(*_hkl())


@pytest.mark.parametrize("method", ["deepcopy", "pickle", "text"])
def test_copy_and_text_preserve_geometry_and_policy(method):
    """Copied models retain strict rejection as well as valid amplitudes."""
    crystal = _crystal(strict=True)
    expected = crystal.F(*_hkl())
    if method == "deepcopy":
        restored = copy.deepcopy(crystal)
    elif method == "pickle":
        restored = pickle.loads(pickle.dumps(crystal))
    else:
        restored = CTRcalc.SXRDCrystal.fromStr(crystal.toStr(showErrors=False))
    np.testing.assert_array_equal(restored.F(*_hkl()), expected)
    restored["film_0"].basis[0] = 12
    with pytest.raises(StackValidationError):
        restored.F(*_hkl())
    np.testing.assert_array_equal(crystal.F(*_hkl()), expected)


@pytest.mark.parametrize("policy", [{"version": 2, "overlap": "error"},
                                    {"version": 1, "overlap": "allow"}, {}])
def test_unsupported_policy_is_not_silently_downgraded(policy):
    """Unknown text/HDF5 policies fail before updating model parameters."""
    crystal = _crystal()
    text = crystal.toStr(False)
    with pytest.raises(StackValidationError):
        CTRcalc.SXRDCrystal.fromStr(
            text + "\n// StackPolicy " + json.dumps(policy))
    parameters = crystal.parametersToDict()
    parameters["stack_policy"] = policy
    with pytest.raises(StackValidationError):
        crystal.parametersFromDict(parameters)


def test_hdf5_policy_roundtrip(tmp_path):
    """Companion settings upgrade a legacy-loaded model without losing checks."""
    from silx.io import dictdump

    crystal = _crystal(strict=True)
    target = tmp_path / "settings.h5"
    dictdump.dicttoh5(crystal.parametersToDict(), str(target))
    parameters = dictdump.h5todict(str(target))
    restored = _crystal()
    restored.parametersFromDict(parameters)
    np.testing.assert_array_equal(restored.F(*_hkl()), crystal.F(*_hkl()))
    restored["film_0"].basis[0] = 12
    with pytest.raises(StackValidationError):
        restored.F(*_hkl())


def test_legacy_files_keep_legacy_semantics():
    """No policy metadata means the established unvalidated constructor."""
    crystal = _crystal(widths=(12, 40))
    restored = CTRcalc.SXRDCrystal.fromStr(crystal.toStr(False))
    assert restored._stack_policy is None
    np.testing.assert_array_equal(restored.F(*_hkl()), crystal.F(*_hkl()))
    # Standalone objects restored from older pickles lack the new derived
    # bounds. Their existing position setter must still translate domains.
    interface = restored["interface_0"]
    del interface._support_bounds
    del interface._replacement_bounds
    old_position = interface.pos_absolute
    interface.pos_absolute = old_position + 1.
    assert interface.pos_absolute == pytest.approx(old_position + 1.)


@pytest.mark.parametrize("coupling", [0., .4, 1.])
@pytest.mark.parametrize("offset", [-.2, .2])
@pytest.mark.parametrize("positions,shift", [
    ((0., .5), (0., 0., 1.)),
    ((.1, .3, .8), (.5, .25, 1.)),
])
def test_strict_preserves_strained_and_cyclic_legacy_geometry(
    coupling, offset, positions, shift,
):
    """Validation leaves phase, unequal intervals and strain fields unchanged."""
    bulk, parts = _parts(widths=(80, 81, 82), repeats=(6., 5.2, 7.1, 6.8),
                         positions=positions, shift=shift)
    for interface in parts[::2]:
        interface.basis[1:] = [.1, coupling, offset]
        interface.basis_0[:] = interface.basis
    legacy = CTRcalc.SXRDCrystal(bulk, *parts, stacking=np.arange(1, 7))
    strict = CTRcalc.SXRDCrystal.from_stack(copy.deepcopy(bulk), copy.deepcopy(parts))
    for crystal in (legacy, strict):
        crystal.setEnergy(18000.)
    np.testing.assert_array_equal(strict.F(*_hkl()), legacy.F(*_hkl()))
    z = np.linspace(-50, 800, 701)
    np.testing.assert_array_equal(strict.zDensity_G(z, 0., 0.),
                                  legacy.zDensity_G(z, 0., 0.))
    np.testing.assert_array_equal(strict.optical_profile(), legacy.optical_profile())


@pytest.mark.parametrize("mutation", ["roughness", "skew", "offset", "tail"])
def test_profile_changes_revalidate_support(mutation):
    """Overlap is checked from generated geometry, not only Film width."""
    crystal = _crystal(strict=True, widths=(40, 80))
    interface = crystal["interface_1"]
    if mutation == "roughness":
        interface.basis[0] = 5.
    elif mutation == "skew":
        interface.basis[1] = -100.
    elif mutation == "offset":
        interface.basis[2:] = [0., -30.]
    else:
        # At W=16, loose supports are disjoint; tighter tails overlap.
        crystal = _crystal(strict=True, widths=(16, 40), tail=1e-3)
        for interface in crystal.uc_surface_list[::2]:
            interface.profile = SkellamProfile(.5, tail_probability=1e-10)
    with pytest.raises(StackValidationError):
        crystal.F(*_hkl())


def test_coupled_material_parameter_survives_copy_and_settings(tmp_path):
    """One atomic parameter reaches the Film and both matching interface sides."""
    from silx.io import dictdump

    crystal = _crystal(strict=True)
    crystal.addFitParameter({
        "interface_0": {"atoms": [0], "par": ["oDW"],
                        "settings": {"unitcell": "top"}},
        "film_0": {"atoms": [0], "par": ["oDW"]},
        "interface_1": {"atoms": [0], "par": ["oDW"],
                        "settings": {"unitcell": "bottom"}},
    }, name="A_disorder", limits=(0., 2.))
    crystal["interface_1"].addFitParameter("S", name="upper_skew", limits=(-.2, .2))
    crystal.getStartParamAndLimits()
    crystal.setParameters([.3, .05])
    target = tmp_path / "coupled.h5"
    dictdump.dicttoh5(crystal.parametersToDict(), str(target))
    restored = CTRcalc.SXRDCrystal.fromStr(crystal.toStr(False))
    restored.parametersFromDict(dictdump.h5todict(str(target)))
    for candidate in (crystal, copy.deepcopy(crystal),
                      pickle.loads(pickle.dumps(crystal)), restored):
        candidate.setParameters([.7, -.05])
        values = [candidate["interface_0"].uc_top.basis[0, 5],
                  candidate["film_0"].unitcell.basis[0, 5],
                  candidate["interface_1"].uc_bottom.basis[0, 5]]
        np.testing.assert_allclose(values, [.7] * 3)
        np.testing.assert_allclose(candidate.F(*_hkl()), crystal.F(*_hkl()),
                                   rtol=0, atol=1e-10)


def test_dwba_multilayer_records_and_rejected_candidate_recovery():
    """Signed optical shares and prepared DWBA survive valid-invalid-valid use."""
    from .. import CTRuc

    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA requires the native CTR extension")
    crystal = _crystal(strict=True)
    prepared = crystal.dwba.prepare_from_glancing(0., 0., .04, .04)
    shares = prepared._atomic_model.reference_shares
    np.testing.assert_allclose(shares.sum(axis=0),
                               prepared.reference.profile.values[::-1, 1:],
                               rtol=2e-14, atol=1e-18)
    expected = crystal.dwba.evaluate_prepared(prepared)
    np.testing.assert_allclose(expected.F_h,
                               expected.F_atomic - expected.F_reference)
    crystal["film_0"].basis[0] = 12
    with pytest.raises(StackValidationError):
        crystal.dwba.evaluate_prepared(prepared)
    with pytest.raises(StackValidationError):
        crystal.dwba.prepare_from_glancing(0., 0., .04, .04)
    crystal["film_0"].basis[0] = 40
    assert crystal.dwba.evaluate_prepared(prepared).F_h == expected.F_h
    assert crystal.dwba.prepare_from_glancing(0., 0., .04, .04) is prepared


def test_dwba_sharp_lower_material_cancels_for_every_interface():
    """Both buried interfaces cancel their sharp lower correction in DWBA."""
    from .. import CTRuc

    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA requires the native CTR extension")
    crystal = _crystal(strict=True, roughness=0.)
    prepared = crystal.dwba.prepare_from_glancing(
        np.array([0., .6]), 0., np.array([.04, .04]), np.array([.04, .05]))
    result = crystal.dwba.evaluate_prepared(prepared)
    seen = set()
    for record, share, part in zip(crystal.dwba._live_records(),
                                    prepared._atomic_model.reference_shares,
                                    result.contributions):
        if record.role != "interface_bottom":
            continue
        seen.add(record.component_name)
        np.testing.assert_allclose(share, 0., atol=1e-18)
        np.testing.assert_allclose(part.F_atomic, 0., atol=2e-12)
        np.testing.assert_allclose(part.F_reference, 0., atol=2e-12)
    assert seen == {"interface_0", "interface_1"}
    np.testing.assert_allclose(
        sum(part.F_h for part in result.contributions), result.F_h,
        rtol=1e-13, atol=1e-12)


def test_optimizer_reports_invalid_candidate_and_recovers():
    """The optimizer's copied model propagates geometry errors explicitly."""
    from .. import CTRopt, CTRplotutil

    crystal = _crystal(strict=True)
    crystal["film_0"].addFitParameter("W", name="middle_width", limits=(6., 60.))
    h, k, ell = _hkl()
    measured = np.abs(crystal.F(h, k, ell))
    rod = CTRplotutil.CTR((1., 0.), ell, measured, np.ones_like(ell))
    optimizer = CTRopt.CTROptimizer(
        crystal, CTRplotutil.CTRCollection([rod]), scale_policy={"F": "fixed"})
    optimizer.prepareFit()
    expected = optimizer.flat_Fcalc([40.])
    with pytest.raises(StackValidationError):
        optimizer.flat_Fcalc([12.])
    np.testing.assert_array_equal(optimizer.flat_Fcalc([40.]), expected)
    assert crystal["film_0"].basis[0] == 40.


def test_adjacent_support_boundaries_and_insufficient_gap():
    """Half-open supports may touch but must not cross."""
    # At this fixed support, both interfaces extend 6 A to each side of
    # their reference boundary. A four-layer Film gives a 12 A separation.
    bulk, parts = _parts(widths=(4, 40), roughness=0.)
    for interface in parts[::2]:
        interface.fixed_ucs = 2
    crystal = CTRcalc.SXRDCrystal.from_stack(bulk, parts)
    first, second = crystal.uc_surface_list[::2]
    assert first._support_bounds[1] == pytest.approx(second._support_bounds[0])
    crystal["film_0"].basis[0] = 3.
    with pytest.raises(StackValidationError, match="overlaps"):
        crystal.apply_stacking()


def test_reference_area_scaling_preserves_physical_scattering_vector():
    """Doubling the reference lateral cell scales every amplitude by two."""
    crystal = _crystal(strict=True, roughness=0.)
    h, k, ell = _hkl()
    expected = crystal.F(h, k, ell)
    reference = CTRcalc.UnitCell([6., 3., 6.], [90.] * 3, name="reference")
    crystal.setGlobalReferenceUnitCell(reference)
    np.testing.assert_allclose(crystal.F(2 * h, k, ell), 2 * expected,
                               rtol=1e-12, atol=1e-10)


def test_dwba_vacuum_limit_matches_multilayer_kinematic_amplitude():
    """A frozen vacuum wavefield reduces to the coherent finite/bulk sum."""
    from .. import CTRuc

    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA requires the native CTR extension")
    crystal = _crystal(strict=True, roughness=0.)
    cells = [crystal.uc_bulk]
    for interface, film in zip(crystal.uc_surface_list[::2],
                                crystal.uc_surface_list[1::2]):
        cells.extend((interface.uc_top, interface.uc_bottom, film.unitcell))
    for cell in cells:
        cell.basis[:, 6] = 0.
    prepared = crystal.dwba.prepare_from_glancing(.6, 0., .04, .05)
    np.testing.assert_array_equal(prepared.reference.n,
                                   np.ones_like(prepared.reference.n))
    for cell in cells:
        cell.basis[:, 6] = 1.
    actual = crystal.dwba.evaluate_prepared(
        prepared, bulk_attenuation=crystal.atten)
    ell = prepared.k0 * (np.sin(.04) + np.sin(.05)) * 6. / (2 * np.pi)
    # DWBA's incoming/outgoing waves have unit phase at the top optical
    # boundary; kinematic phases use the crystal origin at z=0.
    origin_phase = np.exp(-2j * np.pi * ell / 6. * prepared.field_i.z_reference[0])
    expected = crystal.F(.6, 0., ell)[0] * prepared.cos_azimuth * origin_phase
    np.testing.assert_allclose(actual.F_h, expected, rtol=1e-11, atol=1e-10)


def test_dwba_cached_batch_does_not_conceal_invalid_geometry():
    """Even a forbidden mid-batch mutation cannot bypass strict validation."""
    from .. import CTRuc

    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA requires the native CTR extension")
    crystal = _crystal(strict=True)
    with crystal.dwba.batch():
        prepared = crystal.dwba.prepare_from_glancing(0., 0., .04, .04)
        crystal["film_0"].basis[0] = 12
        with pytest.raises(StackValidationError):
            crystal.dwba.evaluate_prepared(prepared)
        assert crystal.dwba._batch_snapshot is None
    crystal["film_0"].basis[0] = 40
    assert np.isfinite(crystal.dwba.evaluate_prepared(prepared).F_h)


@pytest.mark.parametrize("name,column,change", [
    ("film_0", 0, 1.), ("film_1", 0, 1.),
    ("interface_0", 1, .1), ("interface_1", 1, .1),
])
def test_each_component_invalidates_changed_dwba_preparation(name, column, change):
    """Geometry changes in either Film or interface invalidate frozen geometry."""
    from .. import CTRuc

    if not CTRuc.HAS_CPP_ACCEL:
        pytest.skip("DWBA requires the native CTR extension")
    crystal = _crystal(strict=True)
    prepared = crystal.dwba.prepare_from_glancing(.6, 0., .04, .05)
    expected = crystal.dwba.evaluate_prepared(prepared).F_h
    old = crystal[name].basis[column]
    crystal[name].basis[column] += change
    with pytest.raises(ValueError, match="geometry changed"):
        crystal.dwba.evaluate_prepared(prepared)
    updated = crystal.dwba.prepare_from_glancing(.6, 0., .04, .05)
    assert updated is not prepared
    assert np.isfinite(crystal.dwba.evaluate_prepared(updated).F_h)
    crystal[name].basis[column] = old
    assert crystal.dwba.evaluate_prepared(prepared).F_h == expected


def test_missing_settings_metadata_does_not_disable_strict_policy():
    """Legacy companion settings cannot downgrade an explicitly strict model."""
    crystal = _crystal(strict=True)
    parameters = crystal.parametersToDict()
    del parameters["stack_policy"]
    crystal.parametersFromDict(parameters)
    crystal["film_0"].basis[0] = 12
    with pytest.raises(StackValidationError):
        crystal.F(*_hkl())

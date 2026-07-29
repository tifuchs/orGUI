"""Cancelling a parameter dialog must not leave the edited configuration active.

The machine and crystal parameter editors apply every change immediately, so
restoring the widget values is not sufficient: the saved parameters have to be
applied again when the dialog is cancelled or reset.
"""

import pytest

QUBCalculator = pytest.importorskip(
    "orgui.app.QUBCalculator", reason="Qt bindings are required"
)


class _StubSignal:
    def __init__(self):
        self.emitted = []

    def emit(self, *args):
        self.emitted.append(args)


class _StubMachineParams:
    def __init__(self):
        self.restored = []
        self.sigMachineParamsChanged = _StubSignal()

    def setValues(self, params):
        self.restored.append(params)


class _StubMachineDialog:
    def __init__(self, saved):
        self.savedParams = saved
        self.machineparams = _StubMachineParams()


class _StubCrystalParams:
    def __init__(self):
        self.restored = []
        self.sigCrystalParamsChanged = _StubSignal()

    def setValues(self, crystal, n):
        self.restored.append((crystal, n))


class _StubCrystalDialog:
    def __init__(self, saved):
        self.savedParams = saved
        self.crystalparams = _StubCrystalParams()


def _machine_params():
    return {
        "diffractometer": {"mu": 0.1, "chi": 0.2, "phi": 0.3},
        "source": {"E": 78.0},
        "SXRD_geometry": object(),
    }


def test_reset_reapplies_the_saved_machine_parameters():
    saved = _machine_params()
    dialog = _StubMachineDialog(saved)

    QUBCalculator.QMachineParametersDialog.resetParameters(dialog)

    assert dialog.machineparams.restored == [saved]
    # applied again, otherwise the edited configuration would stay active
    assert dialog.machineparams.sigMachineParamsChanged.emitted == [(saved,)]


def test_machine_parameters_are_reapplied_unrounded():
    """The saved snapshot is emitted, not the value read back from a spin box."""
    saved = _machine_params()
    dialog = _StubMachineDialog(saved)

    QUBCalculator.QMachineParametersDialog.resetParameters(dialog)

    (emitted,) = dialog.machineparams.sigMachineParamsChanged.emitted[0]
    assert emitted is saved


def test_reset_without_a_snapshot_does_nothing():
    dialog = _StubMachineDialog(None)

    QUBCalculator.QMachineParametersDialog.resetParameters(dialog)

    assert dialog.machineparams.restored == []
    assert dialog.machineparams.sigMachineParamsChanged.emitted == []


def test_reset_reapplies_the_saved_crystal():
    crystal, refraction_index = object(), 0.9999
    dialog = _StubCrystalDialog((crystal, refraction_index))

    QUBCalculator.QCrystalParameterDialog.resetParameters(dialog)

    assert dialog.crystalparams.restored == [(crystal, refraction_index)]
    assert dialog.crystalparams.sigCrystalParamsChanged.emitted == [
        (crystal, refraction_index)
    ]


def test_crystal_reset_without_a_snapshot_does_nothing():
    dialog = _StubCrystalDialog(None)

    QUBCalculator.QCrystalParameterDialog.resetParameters(dialog)

    assert dialog.crystalparams.restored == []
    assert dialog.crystalparams.sigCrystalParamsChanged.emitted == []

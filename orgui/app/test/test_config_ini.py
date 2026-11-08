import numpy as np
import pytest

from orgui.app.config_data import ConfigData
from orgui.datautils.xrayutils import CTRcalc, DetectorCalibration


def _config(tmp_path, machine='', lattice=''):
    path = tmp_path / 'config.ini'
    path.write_text(f'[Machine]\n{machine}\n[Lattice]\n{lattice}\n[Diffractometer]\n')
    return path


def test_poni_overrides_machine_and_resolves_relative_path(tmp_path, monkeypatch):
    """Resolve PONI beside the INI and use its geometry and wavelength."""
    detector = DetectorCalibration.Detector2D_SXRD()
    detector.setFit2D(540, 20, 30, pixelX=100, pixelY=100)
    detector.wavelength = 1.2e-10
    detector.save(str(tmp_path / 'det.poni'))
    path = _config(tmp_path, 'poni=det.poni\nE=77\nSDD=0.9')
    monkeypatch.chdir(tmp_path.parent)
    result = ConfigData.from_ini(path)
    assert result.detector.dist == pytest.approx(detector.dist)
    assert result.detector.wavelength == pytest.approx(1.2e-10)
    assert result.ub_calculator.getLambda() == pytest.approx(1.2)
    assert result.unit_cell._E == pytest.approx(result.detector.energy * 1000)


def test_external_xtal_and_lattice_override(tmp_path, monkeypatch):
    """Preserve crystal atoms and apply only a complete lattice override."""
    bulk = CTRcalc.UnitCell([3, 4, 5], [90, 90, 90])
    bulk.addAtom('O', [0.2, 0.3, 0.4], 0.1, 0.1, 1)
    bulk.setEnergy(70000)
    CTRcalc.SXRDCrystal(bulk).toFile(str(tmp_path / 'cell.xtal'))
    path = _config(tmp_path, lattice='crystal=cell.xtal\na1=8')
    monkeypatch.chdir(tmp_path.parent)
    result = ConfigData.from_ini(path)
    np.testing.assert_allclose(result.unit_cell.a, [3, 4, 5])
    np.testing.assert_allclose(result.unit_cell.basis[0, 1:4], [0.2, 0.3, 0.4])
    path = _config(tmp_path, lattice=('crystal=cell.xtal\na1=6\na2=7\na3=8\n'
                                      'alpha1=90\nalpha2=91\nalpha3=92'))
    result = ConfigData.from_ini(path)
    np.testing.assert_allclose(result.unit_cell.a, [6, 7, 8])
    np.testing.assert_allclose(result.unit_cell.basis[0, 1:4], [0.2, 0.3, 0.4])


def test_bundled_crystal_and_inline_defaults(tmp_path):
    """Load bundled crystals while retaining partial inline lattice defaults."""
    result = ConfigData.from_ini(_config(tmp_path, lattice='crystal=Pt100'))
    np.testing.assert_allclose(result.unit_cell.a, [3.9242] * 3)
    assert len(result.unit_cell.basis) == 4
    result = ConfigData.from_ini(_config(tmp_path, lattice='a1=2'))
    np.testing.assert_allclose(result.unit_cell.a, [2, 1, 1])


@pytest.mark.parametrize('machine,lattice', [
    ('poni=missing.poni', ''), ('', 'crystal=missing.bul'),
    ('', 'crystal=missing.xtal'),
])
def test_missing_requested_file_fails(tmp_path, machine, lattice):
    """Reject missing requested files instead of using scientific defaults."""
    with pytest.raises((OSError, ValueError)):
        ConfigData.from_ini(_config(tmp_path, machine, lattice))


def test_invalid_poni_fails(tmp_path):
    """Reject a file that does not contain a usable PONI calibration."""
    (tmp_path / 'bad.poni').write_text('not a calibration')
    with pytest.raises((OSError, ValueError, KeyError)):
        ConfigData.from_ini(_config(tmp_path, 'poni=bad.poni'))



def test_invalid_crystal_fails(tmp_path):
    """Report an invalid requested crystal with its file context."""
    (tmp_path / 'bad.xtal').write_text('not a crystal')
    with pytest.raises(ValueError, match="Cannot read crystal file"):
        ConfigData.from_ini(_config(tmp_path, lattice='crystal=bad.xtal'))

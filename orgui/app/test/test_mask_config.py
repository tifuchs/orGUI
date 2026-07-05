import configparser

import numpy as np
import pytest

from orgui.app.mask_config import MaskManager, parse_mask_settings


def test_mask_without_pixel_repair_disables_repair(tmp_path):
    mask_file = tmp_path / "mask.npy"
    np.save(mask_file, np.zeros((2, 3), dtype=bool))
    config = configparser.ConfigParser()
    config.read_dict({"Mask": {"mask": str(mask_file)}})

    settings = parse_mask_settings(config)

    assert settings.mask == str(mask_file)
    assert settings.pixel_repair is None
    manager = MaskManager()
    manager.settings = settings
    assert manager.repair_enabled(cpp_available=True) is False


def test_malformed_optional_pixel_repair_values_warn_and_default(caplog):
    config = configparser.ConfigParser()
    config.read_dict(
        {
            "Mask": {},
            "Mask.PixelRepair": {
                "enabled": "maybe",
                "max_component_pixels": "-1",
                "radius": "wide",
            }
        }
    )

    settings = parse_mask_settings(config)

    assert settings.pixel_repair.enabled is False
    assert settings.pixel_repair.max_component_pixels == 4
    assert settings.pixel_repair.radius == 2
    assert "Invalid mask configuration value" in caplog.text


def test_mask_manager_loads_numpy_mask_with_fabio(tmp_path):
    mask_file = tmp_path / "mask.npy"
    np.save(mask_file, np.array([[0, 1], [2, 0]], dtype=np.uint8))

    manager = MaskManager()
    mask = manager.load_mask(mask_file)

    np.testing.assert_array_equal(mask, [[False, True], [True, False]])
    assert mask.dtype == bool


def test_mask_manager_loads_edf_mask_with_fabio(tmp_path):
    fabio = pytest.importorskip("fabio")
    mask_file = tmp_path / "mask.edf"
    fabio.edfimage.edfimage(
        data=np.array([[0, 1], [3, 0]], dtype=np.int8)
    ).write(mask_file)

    manager = MaskManager()
    mask = manager.load_mask(mask_file)

    np.testing.assert_array_equal(mask, [[False, True], [True, False]])
    assert mask.dtype == bool


def test_repair_enabled_with_missing_cpp_warns_only(caplog):
    config = configparser.ConfigParser()
    config.read_dict({"Mask": {}, "Mask.PixelRepair": {"enabled": "True"}})
    manager = MaskManager()
    manager.settings = parse_mask_settings(config)

    assert manager.repair_enabled(cpp_available=False) is False
    assert "requires the C++ ROI backend" in caplog.text


def test_pyfai_gap_intervals_are_int32_from_full_stripes():
    class Detector:
        def calc_mask(self):
            mask = np.zeros((6, 7), dtype=np.int8)
            mask[2:4, :] = 1
            mask[:, 5:6] = 1
            mask[0, 0] = 1
            return mask

    config = configparser.ConfigParser()
    config.read_dict({"Mask.PixelRepair": {"enabled": "True"}})
    manager = MaskManager()
    manager.settings = parse_mask_settings(config)

    row_gaps, col_gaps = manager.pyfai_gap_intervals(Detector())

    assert row_gaps.dtype == np.int32
    assert col_gaps.dtype == np.int32
    np.testing.assert_array_equal(row_gaps, [[2, 4]])
    np.testing.assert_array_equal(col_gaps, [[5, 6]])


def test_pyfai_gap_intervals_fall_back_to_detector_module_metadata():
    class Detector:
        shape = (7, 8)
        module_size = (2, 3)
        MODULE_GAP = (1, 2)

    config = configparser.ConfigParser()
    config.read_dict({"Mask": {}, "Mask.PixelRepair": {"enabled": "True"}})
    manager = MaskManager()
    manager.settings = parse_mask_settings(config)

    row_gaps, col_gaps = manager.pyfai_gap_intervals(Detector())

    np.testing.assert_array_equal(row_gaps, [[2, 3], [5, 6]])
    np.testing.assert_array_equal(col_gaps, [[3, 5]])


def test_gap_size_px_is_used_when_module_gap_width_is_missing():
    class Detector:
        shape = (8, 9)
        module_size = (3, 4)

    config = configparser.ConfigParser()
    config.read_dict(
        {
            "Mask": {},
            "Mask.PixelRepair": {"enabled": "True", "gap_size_px": "2"},
        }
    )
    manager = MaskManager()
    manager.settings = parse_mask_settings(config)

    row_gaps, col_gaps = manager.pyfai_gap_intervals(Detector())

    np.testing.assert_array_equal(row_gaps, [[3, 5]])
    np.testing.assert_array_equal(col_gaps, [[4, 6]])


def test_mask_manager_shape_mismatch_returns_none(caplog):
    manager = MaskManager()
    manager.set_mask(np.zeros((2, 3), dtype=bool))

    assert manager.get_mask((3, 2)) is None
    assert "does not match image shape" in caplog.text

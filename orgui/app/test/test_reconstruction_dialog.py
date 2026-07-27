"""Regression tests for reconstruction-dialog error containment."""

import logging
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np
from silx.gui import qt

import orgui.app.ReconstructionDialog as reconstruction_dialog_module
from orgui.app.ReconstructionDialog import ReconstructionDialog
from orgui.app.config_data import CorrectionState
from orgui.app.database import FILTERS
from orgui.app.HDF5SettingsDialog import HDF5SettingsDialog


def _dialog(tmp_path):
    app = qt.QApplication.instance() or qt.QApplication([])
    parent = qt.QWidget()
    orgui = SimpleNamespace(
        fscan=None,
        filedialogdir=str(tmp_path),
        numberthreads=1,
        maxMemory=128,
        database=SimpleNamespace(compression=FILTERS["Raw"]),
        _onShowMaskConfig=lambda: None,
        backgroundImageAct=SimpleNamespace(trigger=lambda: None),
        excludedImagesDialog=SimpleNamespace(show=lambda: None),
    )
    dialog = ReconstructionDialog(orgui, parent=parent)
    dialog._test_app = app
    dialog._test_parent = parent
    return dialog


def test_no_scan_actions_are_reported_without_raising(tmp_path, caplog):
    """User actions without an active scan must remain non-fatal."""
    dialog = _dialog(tmp_path)
    caplog.set_level(logging.WARNING)

    dialog.add_derived_grid("hkl")
    dialog.preview()
    dialog.prepare()
    dialog.run_local()
    dialog.resume()
    dialog.refresh_live_state()

    assert "No active scan" in dialog.experiment_summary.toPlainText()
    records = [
        record for record in caplog.records if hasattr(record, "show_dialog")
    ]
    assert len(records) >= 5
    assert all(record.show_dialog is True for record in records)
    assert all(record.dialog_level == logging.WARNING for record in records)

    dialog.close()
    dialog._test_parent.close()


def test_auxiliary_callback_failures_are_contained(
    tmp_path, monkeypatch, caplog
):
    """File, shared-setting, and grid callback failures must not escape."""
    dialog = _dialog(tmp_path)
    caplog.set_level(logging.WARNING)

    monkeypatch.setattr(
        qt.QFileDialog,
        "getExistingDirectory",
        Mock(side_effect=RuntimeError("file dialog failed")),
    )
    assert (
        dialog._browse(
            dialog.scratch_path,
            kind="directory",
            save=False,
        )
        is False
    )

    dialog._invoke_ui_action(
        "Cannot open test settings",
        Mock(side_effect=RuntimeError("settings failed")),
    )

    dialog.grid_table.insertRow(0)
    dialog.orgui.fscan = object()
    dialog.preview()

    records = [
        record for record in caplog.records if hasattr(record, "show_dialog")
    ]
    assert len(records) == 3
    assert "empty cell" in dialog.preview_output.toPlainText()

    dialog.close()
    dialog._test_parent.close()


def test_default_paths_use_working_directory(tmp_path, monkeypatch):
    """Generated writable paths must not default beside input data."""
    dialog = _dialog(tmp_path)
    history_file = tmp_path / "La3Ni2O7"
    history_file.touch()
    dialog.orgui.filedialogdir = str(history_file)
    working_directory = tmp_path / "working"
    working_directory.mkdir()
    monkeypatch.chdir(working_directory)

    dialog._set_default_paths("39_1")

    assert dialog.job_path.text() == str(
        working_directory / "39_1-rsmap.json"
    )
    assert dialog.scratch_path.text() == str(
        working_directory / ".39_1-rsmap-scratch"
    )
    assert dialog.output_path.text() == str(
        working_directory / "39_1-rsmap.h5"
    )

    dialog.close()
    dialog._test_parent.close()


def test_detected_performance_values_are_visible_without_becoming_overrides(
    tmp_path, monkeypatch
):
    """Detected execution values populate disabled advanced editors."""
    dialog = _dialog(tmp_path)
    detected = {
        "thread_budget": 16,
        "native_threads_per_image": 4,
        "memory_budget_MiB": 2048.0,
        "accumulation_budget_MiB_per_worker": 256.0,
        "frames_per_task": 1,
        "detector_tile_shape": (1024, 768),
        "native_work_block_pixels": 65536,
        "parquet_chunk_span": 256,
        "frame_tasks": 20,
        "detector_tiles": 12,
        "map_tasks": 240,
        "parallel_layouts": [
            {
                "exposure": "stationary",
                "concurrent_image_workers": 8,
                "native_threads_per_image": 1,
                "tiles_per_image": 12,
                "memory_per_image_MiB": 256.0,
                "accumulation_MiB_per_worker": 256.0,
            }
        ],
    }
    monkeypatch.setattr(
        reconstruction_dialog_module,
        "reconstruction_execution_settings",
        lambda job, scan=None, config=None: detected,
    )

    dialog._show_execution_settings(object())

    assert dialog.thread_override[2].value() == 16
    assert dialog.threads_per_image.value() == 4
    assert dialog.memory_override[2].value() == 2048
    assert dialog.accumulation_memory[2].value() == 256
    assert dialog.tile_rows[2].value() == 1024
    assert dialog.tile_columns[2].value() == 768
    assert not dialog.thread_override[1].isChecked()
    assert not dialog.tile_rows[1].isChecked()
    assert '"concurrent_image_workers": 8' in (
        dialog.performance_summary.toPlainText()
    )

    dialog.close()
    dialog._test_parent.close()


def test_settings_are_grouped_and_have_tooltips(tmp_path):
    """Reconstruction settings use concise grouped labels and help text."""
    dialog = _dialog(tmp_path)

    group_titles = {
        group.title() for group in dialog.findChildren(qt.QGroupBox)
    }
    assert {
        "Active experiment",
        "Corrections and exclusions",
        "Exposure and job metadata",
        "Grid definitions",
        "Accuracy",
        "Parallel execution and memory",
        "Advanced settings",
        "Detected execution layout",
        "Job descriptor",
        "Storage",
    } <= group_titles

    labels = [
        label.text() for label in dialog.findChildren(qt.QLabel)
    ]
    assert not any(label.startswith("Advanced ") for label in labels)
    assert not any(label.startswith("One-job ") for label in labels)
    assert [
        dialog.accuracy.itemText(index)
        for index in range(dialog.accuracy.count())
    ] == [
        "Center only (depth 0)",
        "Low (depth 1)",
        "Balanced (depth 2)",
        "High (depth 3)",
        "Very high (depth 4)",
        "Maximum (depth 5)",
    ]
    assert [
        reconstruction_dialog_module.ACCURACY_DEPTHS[
            dialog.accuracy.itemData(index)
        ]
        for index in range(dialog.accuracy.count())
    ] == list(range(6))
    assert dialog.accuracy.currentData() == "balanced"

    setting_controls = [
        dialog.angle_fallback,
        dialog.user_note,
        dialog.grid_table,
        dialog.accuracy,
        dialog.thread_override[0],
        dialog.threads_per_image,
        dialog.memory_override[0],
        dialog.accumulation_memory[0],
        dialog.frame_batch[0],
        dialog.tile_rows[0],
        dialog.tile_columns[0],
        dialog.work_block[0],
        dialog.partition_span[0],
        dialog.performance_summary,
        dialog.job_path,
        dialog.scratch_path,
        dialog.output_path,
    ]
    assert all(control.toolTip() for control in setting_controls)
    assert all(
        dialog.grid_table.horizontalHeaderItem(column).toolTip()
        for column in range(dialog.grid_table.columnCount())
    )

    dialog.close()
    dialog._test_parent.close()


def test_json_output_has_its_own_tab_and_is_selected_when_updated(tmp_path):
    """JSON and status output must not compress the settings tabs."""
    dialog = _dialog(tmp_path)

    assert dialog.tabs.indexOf(dialog.output_tab) >= 0
    assert dialog.tabs.tabText(
        dialog.tabs.indexOf(dialog.output_tab)
    ) == "Preview and status"
    assert dialog.preview_output.parentWidget() is dialog.output_tab

    dialog._show_output('{"status": "prepared"}')

    assert dialog.tabs.currentWidget() is dialog.output_tab
    assert dialog.preview_output.toPlainText() == '{"status": "prepared"}'

    dialog.close()
    dialog._test_parent.close()


def test_grid_numbers_are_compact_without_losing_values_and_chunks_are_global(
    tmp_path,
):
    """The table display is compact while its editable values remain exact."""
    dialog = _dialog(tmp_path)
    dialog.orgui.reconstruction_chunk_shape = (32, 64, 128)
    grid = reconstruction_dialog_module.ReconstructionGrid(
        minimum=(0.1234567890123, -1.234567890123, 2.345678901234),
        maximum=(1.234567890123, 2.345678901234, 3.456789012345),
        step=(0.001234567890123, 0.002345678901234, 0.003456789012345),
        frame="hkl",
    )

    dialog._append_grid(grid)

    assert dialog.grid_table.columnCount() == 15
    assert dialog.grid_table.item(0, 2).text() == "0.123457"
    assert (
        dialog.grid_table.item(0, 2).data(qt.Qt.EditRole)
        == grid.minimum[0]
    )
    rebuilt = dialog._grids()[0]
    assert rebuilt.minimum == grid.minimum
    assert rebuilt.maximum == grid.maximum
    assert rebuilt.step == grid.step
    assert rebuilt.chunk_shape == (32, 64, 128)
    assert all(
        int(dialog.grid_table.item(0, column).data(qt.Qt.EditRole)) > 0
        for column in range(11, 14)
    )
    assert "uncompressed" in dialog.grid_table.item(0, 14).text()

    dialog.grid_table.item(0, 11).setData(qt.Qt.EditRole, 100)

    rebuilt = dialog._grids()[0]
    assert rebuilt.to_spec().shape[0] == 100
    assert rebuilt.step[0] == np.nextafter(
        (grid.maximum[0] - grid.minimum[0]) / 100, np.inf
    )

    dialog.close()
    dialog._test_parent.close()


def test_hdf5_settings_expose_global_chunks_and_compression_override(tmp_path):
    """Shared HDF5 settings cover chunking and optional output compression."""
    app = qt.QApplication.instance() or qt.QApplication([])
    dialog = HDF5SettingsDialog(
        FILTERS["Raw"],
        chunk_shape=(32, 64, 128),
        compression_override="GZip",
    )
    dialog._test_app = app

    assert dialog.chunk_shape == (32, 64, 128)
    assert dialog.database_compression_name == "Raw"
    assert dialog.compression_override == "GZip"
    dialog.override_compression.setChecked(False)
    assert dialog.compression_override is None

    dialog.close()


def test_geometry_resolution_dialog_exposes_percentile(tmp_path):
    """The local-Jacobian estimator exposes its robust percentile."""
    app = qt.QApplication.instance() or qt.QApplication([])
    dialog = reconstruction_dialog_module._GeometryResolutionDialog()
    dialog._test_app = app

    assert dialog.percentile == 10.0
    dialog.percentile_editor.setValue(25.0)
    assert dialog.percentile == 25.0

    dialog.close()


def test_geometry_step_estimate_updates_steps_intervals_and_size(
    tmp_path, monkeypatch
):
    """Applying a geometry estimate keeps all derived grid fields consistent."""
    dialog = _dialog(tmp_path)
    dialog.orgui.fscan = object()
    dialog._append_grid(
        reconstruction_dialog_module.ReconstructionGrid(
            minimum=(0.0, 0.0, 0.0),
            maximum=(1.0, 2.0, 3.0),
            step=(0.5, 0.5, 0.5),
            frame="hkl",
        )
    )
    captured = {}

    monkeypatch.setattr(
        reconstruction_dialog_module.ConfigData,
        "from_gui",
        lambda gui: object(),
    )

    def estimate(config, scan, *, frame, percentile):
        captured.update(frame=frame, percentile=percentile)
        return (0.1, 0.2, 0.3)

    monkeypatch.setattr(
        reconstruction_dialog_module, "estimate_geometry_steps", estimate
    )

    dialog._apply_geometry_steps([0], 15.0)

    assert captured == {"frame": "hkl", "percentile": 15.0}
    grid = dialog._grids()[0]
    np.testing.assert_allclose(grid.step, (0.1, 0.2, 0.3))
    assert grid.to_spec().shape == (10, 10, 10)
    assert "uncompressed" in dialog.grid_table.item(0, 14).text()

    dialog.close()
    dialog._test_parent.close()


def test_live_summary_resolves_active_mask_instead_of_showing_null_asset(
    tmp_path,
):
    """The live overview distinguishes active inputs from frozen assets."""
    dialog = _dialog(tmp_path)
    mask = np.zeros((5, 7), dtype=bool)
    mask[1, 2] = True
    mask[3, 4] = True
    dialog.orgui.get_detector_mask = lambda shape: mask
    config = SimpleNamespace(
        corrections=CorrectionState(use_mask=True)
    )

    summary = dialog._live_correction_summary(config, mask.shape)

    assert "mask_asset" not in summary
    assert "background_asset" not in summary
    assert "background_variance_asset" not in summary
    assert summary["active_inputs"]["mask"] == {
        "status": "active",
        "masked_pixels": 2,
        "job_asset_path": "/mask",
    }
    assert summary["active_inputs"]["background"] == {
        "status": "disabled"
    }

    dialog.close()
    dialog._test_parent.close()

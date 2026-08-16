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
from orgui.reconstruction_job import ReconstructionGrid


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
    dialog.create_cluster_scripts()
    dialog.resume()
    dialog.refresh_live_state()

    assert "No active scan" in dialog.experiment_summary.toPlainText()
    records = [
        record for record in caplog.records if hasattr(record, "show_dialog")
    ]
    assert len(records) >= 6
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
    assert dialog.threads_per_image[2].value() == 4
    assert dialog.memory_override[2].value() == 2048
    assert dialog.accumulation_memory[2].value() == 256
    assert dialog.tile_rows[2].value() == 1024
    assert dialog.tile_columns[2].value() == 768
    assert not dialog.thread_override[1].isChecked()
    assert not dialog.threads_per_image[1].isChecked()
    assert not dialog.tile_rows[1].isChecked()
    assert '"concurrent_image_workers": 8' in (
        dialog.performance_summary.toPlainText()
    )

    dialog.close()
    dialog._test_parent.close()


def test_file_count_summary_prompts_for_scan_and_grid(tmp_path):
    """Sec14: the estimate must guide the user rather than fail silently
    or raise before a scan/grid exists to estimate from."""
    dialog = _dialog(tmp_path)

    assert "Load a scan" in dialog.file_count_summary.text()

    class _LenScan:
        def __len__(self):
            return 4

    dialog.orgui.fscan = _LenScan()
    dialog._refresh_file_count_summary()
    assert "Add an output grid" in dialog.file_count_summary.text()

    dialog.close()
    dialog._test_parent.close()


def test_file_count_summary_reports_single_node_multi_grid_and_cluster(
    tmp_path, monkeypatch
):
    """Sec14: number_of_files, per-grid breakdown (only when there is more
    than one grid), and the cluster-node multiplication must all be
    reflected live, using the dialog's own current settings."""
    dialog = _dialog(tmp_path)

    class _LenScan:
        def __len__(self):
            return 4

    dialog.orgui.fscan = _LenScan()
    fake_config = SimpleNamespace(
        corrections=SimpleNamespace(excluded_frames=())
    )
    monkeypatch.setattr(
        reconstruction_dialog_module.ConfigData,
        "from_gui",
        classmethod(lambda cls, gui: fake_config),
    )
    monkeypatch.setattr(
        reconstruction_dialog_module,
        "_node_excluded_frames",
        lambda *args, **kwargs: set(),
    )

    calls = []

    def fake_estimate(config, scan, grids, **kwargs):
        calls.append(kwargs)
        per_grid = {
            "q_lab": {
                "job_data_bytes_estimate": 1024.0,
                "files_per_job": 3,
            }
        }
        if len(grids) > 1:
            per_grid["q_hkl"] = {
                "job_data_bytes_estimate": 2048.0,
                "files_per_job": 5,
            }
        return {
            "per_grid": per_grid,
            "files_total": sum(
                result["files_per_job"] for result in per_grid.values()
            ),
        }

    monkeypatch.setattr(
        reconstruction_dialog_module, "estimate_checkpoint_plan", fake_estimate
    )

    def grid(name, frame):
        return ReconstructionGrid(
            minimum=(-1.0, -1.0, -1.0),
            maximum=(1.0, 1.0, 1.0),
            step=(0.1, 0.1, 0.1),
            frame=frame,
            name=name,
        )

    dialog.cluster_array_task_count.setValue(1)
    dialog._append_grid(grid("q_lab", "lab"))

    assert "3" in dialog.file_count_summary.text()
    assert "cluster" not in dialog.file_count_summary.text().lower()
    assert "Per grid" not in dialog.file_count_summary.text()

    dialog.cluster_array_task_count.setValue(3)
    assert f"{3 * 3}" in dialog.file_count_summary.text()
    assert "3-node" in dialog.file_count_summary.text()

    dialog._append_grid(grid("q_hkl", "hkl"))
    assert "Per grid" in dialog.file_count_summary.text()
    assert "q_lab" in dialog.file_count_summary.text()
    assert "q_hkl" in dialog.file_count_summary.text()

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
        "Scheduler",
        "Python environment",
        "Mapping array",
        "Reduction and finalization",
        "Scheduler-specific settings",
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
    assert dialog.work_block.currentData() == "medium"
    assert "starting cache-scale" in dialog.work_block.toolTip()

    setting_controls = [
        dialog.use_pixel_mask,
        dialog.use_solid_angle,
        dialog.use_polarization,
        dialog.angle_fallback,
        dialog.user_note,
        dialog.grid_table,
        dialog.accuracy,
        dialog.thread_override[0],
        dialog.threads_per_image[0],
        dialog.memory_override[0],
        dialog.accumulation_memory[0],
        dialog.frame_batch[0],
        dialog.tile_rows[0],
        dialog.tile_columns[0],
        dialog.work_block,
        dialog.checkpoint_count,
        dialog.performance_summary,
        dialog.cluster_scheduler,
        dialog.cluster_job_name,
        dialog.cluster_script_directory,
        dialog.cluster_working_directory,
        dialog.cluster_python,
        dialog.cluster_environment,
        dialog.cluster_array_task_count,
        dialog.cluster_array_cpus,
        dialog.cluster_array_memory,
        dialog.cluster_array_walltime,
        dialog.cluster_array_concurrency,
        dialog.cluster_summary,
        dialog.cluster_reduce_cpus,
        dialog.cluster_reduce_memory,
        dialog.cluster_reduce_walltime,
        dialog.cluster_sge_pe,
        dialog.cluster_sge_memory,
        dialog.cluster_array_directives,
        dialog.cluster_reduce_directives,
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


def test_cluster_tab_uses_sge_defaults_and_separate_reduce_resources(tmp_path):
    """Cluster settings default to SGE and preserve independent resources."""
    dialog = _dialog(tmp_path)

    dialog.cluster_array_cpus.setValue(4)
    dialog.cluster_array_memory.setValue(16)
    dialog.cluster_reduce_cpus.setValue(32)
    dialog.cluster_reduce_memory.setValue(128)
    settings = dialog._cluster_settings()

    assert settings.scheduler == "sge"
    assert settings.array_cpus == 4
    assert settings.array_memory_gib == 16
    assert settings.reduce_cpus == 32
    assert settings.reduce_memory_gib == 128
    assert any(
        dialog.tabs.tabText(index) == "Cluster"
        for index in range(dialog.tabs.count())
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


def test_dialog_layout_is_bounded_and_workflow_ordered(tmp_path):
    """Dense controls must not force the dialog wider than its target size."""
    dialog = _dialog(tmp_path)
    qt_major = int(qt.qVersion().split(".", maxsplit=1)[0])
    target_width = 900 if qt_major >= 6 else 1100
    assert dialog.width() == target_width
    dialog.resize(target_width, 760)
    dialog.show()
    dialog._test_app.processEvents()

    assert dialog.minimumSizeHint().width() <= 1200
    assert dialog.width() <= 1200
    assert [
        dialog.tabs.tabText(index)
        for index in range(dialog.tabs.count())
    ] == [
        "Experiment",
        "Output grids",
        "Performance",
        "Job and output",
        "Cluster",
        "Preview and status",
    ]
    assert (
        dialog.grid_table.horizontalHeader().sectionResizeMode(0)
        == qt.QHeaderView.Interactive
    )
    assert dialog.hdf5_summary.wordWrap()
    assert len(dialog.findChildren(qt.QDialogButtonBox)) == 2

    dialog.tabs.setCurrentIndex(2)
    dialog._test_app.processEvents()
    performance_tab = dialog.tabs.currentWidget()
    assert all(
        group.geometry().bottom() <= performance_tab.height()
        for group in performance_tab.findChildren(qt.QGroupBox)
    )

    dialog.close()
    dialog._test_parent.close()


def test_open_job_restores_all_editable_job_settings(tmp_path, monkeypatch):
    """Opening a job must not leave controls at unrelated live defaults."""
    dialog = _dialog(tmp_path)
    job_path = tmp_path / "prepared.json"
    dialog.job_path.setText(str(job_path))
    dialog.memory_override[1].setChecked(True)
    dialog.threads_per_image[1].setChecked(True)
    job = SimpleNamespace(
        output_path=str(tmp_path / "result.h5"),
        scratch_path=str(tmp_path / "scratch"),
        accuracy="high",
        advanced_depth=None,
        angle_fallback="midpoint",
        user_note="restored note",
        checkpoint_count=17,
        thread_override=12,
        memory_override_bytes=None,
        frame_batch=7,
        tile_shape=(32, 48),
        work_block_pixels="small",
        threads_per_image=None,
        accumulation_budget_bytes=96 * 1024**2,
        config_data=SimpleNamespace(
            corrections=CorrectionState(
                use_mask=True,
                use_solid_angle=True,
                use_polarization=True,
                normalize_exposure=False,
                monitor_corrections=("monitor", "ring"),
            )
        ),
        grids=[
            ReconstructionGrid(
                minimum=(0.0, 0.0, 0.0),
                maximum=(1.0, 2.0, 3.0),
                step=(0.1, 0.2, 0.3),
                frame="hkl",
                name="opened",
            ).__dict__
        ],
        compression="Raw",
        cluster_settings={},
    )
    monkeypatch.setattr(dialog, "_browse", lambda *args, **kwargs: True)
    monkeypatch.setattr(
        reconstruction_dialog_module, "read_job", lambda path: job
    )
    monkeypatch.setattr(
        reconstruction_dialog_module,
        "job_status",
        lambda path: {"status": "prepared"},
    )
    shared_options = {
        "mask": False,
        "solidAngle": False,
        "polarization": False,
    }
    dialog.orgui.scanSelector = SimpleNamespace(
        get_integration_options=lambda: dict(shared_options),
        set_integration_options=lambda values: shared_options.update(values),
    )
    monkeypatch.setattr(dialog, "_show_execution_settings", Mock())

    dialog.open_job()

    assert dialog.accuracy.currentData() == "high", (
        dialog.preview_output.toPlainText()
    )
    assert dialog.angle_fallback.currentData() == "midpoint"
    assert dialog.user_note.text() == "restored note"
    assert dialog.use_pixel_mask.isChecked()
    assert dialog.use_solid_angle.isChecked()
    assert dialog.use_polarization.isChecked()
    assert shared_options == {
        "mask": True,
        "solidAngle": True,
        "polarization": True,
    }
    assert not dialog.normalize_exposure.isChecked()
    assert dialog.monitor_corrections.text() == "monitor, ring"
    assert dialog.orgui.reconstruction_monitor_corrections == (
        "monitor",
        "ring",
    )
    assert dialog.checkpoint_count.value() == 17
    assert dialog._optional_value(dialog.thread_override) == 12
    assert dialog._optional_value(dialog.memory_override) is None
    assert dialog._optional_value(dialog.frame_batch) == 7
    assert dialog._optional_value(dialog.tile_rows) == 32
    assert dialog._optional_value(dialog.tile_columns) == 48
    assert dialog.work_block.currentData() == "small"
    assert dialog._optional_value(dialog.threads_per_image) is None
    assert dialog._optional_value(dialog.accumulation_memory) == 96
    assert dialog._grids()[0].name == "opened"
    assert dialog.output_path.text() == job.output_path
    assert dialog.scratch_path.text() == job.scratch_path

    dialog._set_accuracy("low", 7)
    assert dialog.accuracy.currentText() == "Advanced (depth 7)"
    assert dialog._accuracy_settings() == ("low", 7)

    dialog.close()
    dialog._test_parent.close()


def test_reconstruction_correction_switches_sync_with_integration_options(
    tmp_path,
):
    """Scientific correction switches must reflect and update shared state."""
    dialog = _dialog(tmp_path)
    shared_options = {
        "mask": True,
        "solidAngle": False,
        "polarization": True,
        "advanced": {"unchanged": True},
    }
    dialog.orgui.scanSelector = SimpleNamespace(
        get_integration_options=lambda: dict(shared_options),
        set_integration_options=lambda values: shared_options.update(values),
    )

    dialog._sync_integration_options()

    assert dialog.use_pixel_mask.isChecked()
    assert not dialog.use_solid_angle.isChecked()
    assert dialog.use_polarization.isChecked()

    dialog.use_solid_angle.setChecked(True)
    dialog.use_polarization.setChecked(False)

    assert shared_options == {
        "mask": True,
        "solidAngle": True,
        "polarization": False,
        "advanced": {"unchanged": True},
    }

    dialog.close()
    dialog._test_parent.close()

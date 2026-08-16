"""GUI editor for centralized reciprocal-space reconstruction jobs."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np
from silx.gui import qt

from .. import logger_utils
from ..reconstruction_cluster import (
    ClusterSettings,
    generate_cluster_scripts,
)
from ..reconstruction_job import (
    ACCURACY_DEPTHS,
    WORK_BLOCK_PRESETS,
    ReconstructionGrid,
    _node_excluded_frames,
    derive_grid,
    estimate_checkpoint_plan,
    estimate_geometry_steps,
    job_status,
    prepare_job,
    read_job,
    reconstruction_execution_settings,
    run_job,
)
from .config_data import ConfigData
from .database import FILTERS
from .HDF5SettingsDialog import (
    HDF5SettingsDialog,
    compression_filter_name,
)


_GRID_COLUMNS = (
    "Name",
    "Frame",
    "Min 1",
    "Min 2",
    "Min 3",
    "Max 1",
    "Max 2",
    "Max 3",
    "Step 1",
    "Step 2",
    "Step 3",
    "Intervals 1",
    "Intervals 2",
    "Intervals 3",
    "Est. size",
)
_FRAMES = ("hkl", "lab", "alpha", "omega", "chi", "phi", "crystal")
logger = logging.getLogger(__name__)


class _GridNumberItem(qt.QTableWidgetItem):
    """Display a compact number while retaining its exact editable value."""

    def __init__(self, value):
        super().__init__()
        self.setData(qt.Qt.EditRole, float(value))

    def data(self, role):
        if role == qt.Qt.DisplayRole:
            value = super().data(qt.Qt.EditRole)
            try:
                return format(float(value), ".6g")
            except (TypeError, ValueError):
                return value
        return super().data(role)


class _GeometryResolutionDialog(qt.QDialog):
    """Choose robust sampling for the local-Jacobian step estimate."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Geometry-matched grid steps")
        layout = qt.QVBoxLayout(self)
        explanation = qt.QLabel(
            "Local detector-row, detector-column, and scan-direction "
            "Jacobians are sampled across the experiment. A lower percentile "
            "selects finer steps and protects more high-resolution regions."
        )
        explanation.setWordWrap(True)
        layout.addWidget(explanation)
        form = qt.QFormLayout()
        self.percentile_editor = qt.QDoubleSpinBox()
        self.percentile_editor.setRange(0.1, 100.0)
        self.percentile_editor.setDecimals(1)
        self.percentile_editor.setSingleStep(1.0)
        self.percentile_editor.setValue(10.0)
        self.percentile_editor.setSuffix(" %")
        tooltip = (
            "Percentile of sampled local one-sigma axis resolutions. Lower "
            "values produce finer grids; 10% is a conservative default."
        )
        self.percentile_editor.setToolTip(tooltip)
        percentile_label = qt.QLabel("Resolution percentile:")
        percentile_label.setToolTip(tooltip)
        form.addRow(percentile_label, self.percentile_editor)
        layout.addLayout(form)
        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    @property
    def percentile(self):
        """Return the selected local-resolution percentile."""
        return self.percentile_editor.value()


def _format_size(size_bytes):
    for suffix in ("B", "KiB", "MiB", "GiB", "TiB", "PiB"):
        if size_bytes < 1024 or suffix == "PiB":
            return f"{size_bytes:.3g} {suffix}"
        size_bytes /= 1024


class _ReconstructionCancelled(RuntimeError):
    """Signal cancellation between reconstruction tasks."""


class ReconstructionDialog(qt.QDialog):
    """Configure, prepare, run, and resume centralized reconstruction jobs."""

    def __init__(self, orgui, parent=None):
        super().__init__(parent or orgui)
        self.orgui = orgui
        if not hasattr(self.orgui, "reconstruction_chunk_shape"):
            self.orgui.reconstruction_chunk_shape = (64, 64, 64)
        if not hasattr(self.orgui, "reconstruction_compression_override"):
            self.orgui.reconstruction_compression_override = None
        if not hasattr(self.orgui, "reconstruction_normalize_exposure"):
            self.orgui.reconstruction_normalize_exposure = True
        if not hasattr(self.orgui, "reconstruction_monitor_corrections"):
            self.orgui.reconstruction_monitor_corrections = ()
        self.setWindowTitle("Reciprocal-space reconstruction")
        self.resize(900, 760)
        layout = qt.QVBoxLayout(self)
        self.tabs = qt.QTabWidget()
        self.tabs.setUsesScrollButtons(True)
        self.tabs.setElideMode(qt.Qt.ElideRight)
        self.tabs.addTab(self._scrollable_tab(self._data_tab()), "Experiment")
        self.tabs.addTab(self._grid_tab(), "Output grids")
        self.tabs.addTab(self._performance_tab(), "Performance")
        self.tabs.addTab(self._paths_tab(), "Job and output")
        self.tabs.addTab(self._cluster_tab(), "Cluster")
        self.output_tab = qt.QWidget()
        output_layout = qt.QVBoxLayout(self.output_tab)
        self.preview_output = qt.QPlainTextEdit()
        self.preview_output.setReadOnly(True)
        self.preview_output.setToolTip(
            "Preview estimates, prepared-job JSON, execution results, and "
            "status or error details."
        )
        output_layout.addWidget(self.preview_output)
        self.tabs.addTab(self.output_tab, "Preview and status")
        for signal in (
            self.checkpoint_count.valueChanged,
            self.thread_override[1].toggled,
            self.thread_override[2].valueChanged,
            self.memory_override[1].toggled,
            self.memory_override[2].valueChanged,
            self.cluster_array_task_count.valueChanged,
            self.accuracy.currentIndexChanged,
            self.angle_fallback.currentIndexChanged,
        ):
            signal.connect(self._refresh_file_count_summary)
        self._refresh_file_count_summary()
        layout.addWidget(self.tabs, stretch=1)
        footer = qt.QVBoxLayout()
        current_actions = qt.QDialogButtonBox()
        prepared_actions = qt.QDialogButtonBox(qt.QDialogButtonBox.Close)
        action_rows = (
            (
                current_actions,
                (
                    ("Preview", self.preview),
                    ("Prepare Job", self.prepare),
                    ("Run Locally", self.run_local),
                    ("Create Cluster Scripts", self.create_cluster_scripts),
                ),
            ),
            (
                prepared_actions,
                (("Open Job", self.open_job), ("Resume", self.resume)),
            ),
        )
        self.action_buttons = {}
        for button_box, actions in action_rows:
            for label, slot in actions:
                button = button_box.addButton(
                    label, qt.QDialogButtonBox.ActionRole
                )
                self.action_buttons[label] = button
                button.clicked.connect(slot)
                button.setToolTip(
                    {
                        "Preview": (
                            "Estimate coverage, storage, and execution layout."
                        ),
                        "Prepare Job": (
                            "Freeze the current settings into a resumable job."
                        ),
                        "Run Locally": (
                            "Prepare and execute the configured job locally."
                        ),
                        "Create Cluster Scripts": (
                            "Prepare the job and create an SGE or Slurm map "
                            "array with a dependent reduction/finalization job."
                        ),
                        "Open Job": (
                            "Select an existing reconstruction job JSON file."
                        ),
                        "Resume": (
                            "Verify and continue the selected prepared job."
                        ),
                    }[label]
                )
            footer.addWidget(button_box)
        prepared_actions.rejected.connect(self.reject)
        layout.addLayout(footer)
        self.refresh_live_state()

    @staticmethod
    def _scrollable_tab(widget):
        """Wrap a settings page so smaller screens can reach every control."""
        widget.setSizePolicy(
            qt.QSizePolicy.Preferred,
            qt.QSizePolicy.Ignored,
        )
        scroll = qt.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(qt.QFrame.NoFrame)
        scroll.setWidget(widget)
        return scroll

    def _data_tab(self):
        widget = qt.QWidget()
        layout = qt.QVBoxLayout(widget)

        experiment_group = qt.QGroupBox("Active experiment")
        experiment_layout = qt.QVBoxLayout(experiment_group)
        self.experiment_summary = qt.QPlainTextEdit()
        self.experiment_summary.setReadOnly(True)
        self.experiment_summary.setToolTip(
            "Live summary of the active scan, geometry, UB matrix, corrections, "
            "and global runtime limits."
        )
        experiment_layout.addWidget(self.experiment_summary)
        refresh = qt.QPushButton("Refresh from active orGUI state")
        refresh.setToolTip(
            "Reload the summary and default grid from the current orGUI state."
        )
        refresh.clicked.connect(self.refresh_live_state)
        experiment_layout.addWidget(refresh)
        layout.addWidget(experiment_group)

        correction_group = qt.QGroupBox("Corrections and exclusions")
        correction_layout = qt.QVBoxLayout(correction_group)
        correction_switches = qt.QGridLayout()
        shared_tooltip = (
            "This correction is shared with the ROI integration options. "
            "Changing it here also updates the main integration panel."
        )
        self.use_pixel_mask = qt.QCheckBox("Use pixel mask")
        self.use_pixel_mask.setToolTip(
            shared_tooltip
            + " Masked-pixel repair keeps the mask active when repair is enabled."
        )
        self.use_solid_angle = qt.QCheckBox("Solid angle correction")
        self.use_solid_angle.setToolTip(shared_tooltip)
        self.use_polarization = qt.QCheckBox("Polarization correction")
        self.use_polarization.setToolTip(shared_tooltip)
        for control, key, row, column in (
            (self.use_pixel_mask, "mask", 0, 0),
            (self.use_solid_angle, "solidAngle", 0, 1),
            (self.use_polarization, "polarization", 1, 0),
        ):
            control.toggled.connect(
                lambda checked, option=key: self._set_integration_option(
                    option, checked
                )
            )
            correction_switches.addWidget(control, row, column)
        correction_switches.setColumnStretch(2, 1)
        correction_layout.addLayout(correction_switches)

        shared = qt.QHBoxLayout()
        mask = qt.QPushButton("Mask and repair")
        mask.setToolTip(
            "Open the shared detector-mask and masked-pixel repair settings."
        )
        mask.clicked.connect(
            lambda: self._invoke_ui_action(
                "Cannot open mask settings",
                self.orgui._onShowMaskConfig,
            )
        )
        background = qt.QPushButton("Background")
        background.setToolTip(
            "Select or clear the shared background image and its variance."
        )
        background.clicked.connect(
            lambda: self._invoke_ui_action(
                "Cannot open background settings",
                self.orgui.backgroundImageAct.trigger,
            )
        )
        exclusions = qt.QPushButton("Excluded frames")
        exclusions.setToolTip(
            "Choose scan frames that must not contribute to reconstruction."
        )
        exclusions.clicked.connect(
            lambda: self._invoke_ui_action(
                "Cannot open excluded-frame settings",
                self.orgui.excludedImagesDialog.show,
            )
        )
        shared.addWidget(mask)
        shared.addWidget(background)
        shared.addWidget(exclusions)
        shared.addStretch(1)
        correction_layout.addLayout(shared)
        layout.addWidget(correction_group)
        self._sync_integration_options()

        metadata_group = qt.QGroupBox("Exposure and job metadata")
        form = qt.QFormLayout(metadata_group)
        self.angle_fallback = qt.QComboBox()
        self.angle_fallback.addItem("Stationary exposure", "stationary")
        self.angle_fallback.addItem(
            "Midpoint inference (explicit fallback)", "midpoint"
        )
        self.angle_fallback.setSizeAdjustPolicy(
            qt.QComboBox.AdjustToMinimumContentsLengthWithIcon
        )
        self.angle_fallback.setMinimumContentsLength(24)
        self._add_form_row(
            form,
            "Missing exact angle bounds:",
            self.angle_fallback,
            "Choose how exposure bounds are represented when the scan backend "
            "does not provide exact encoder start and end positions.",
        )
        self.user_note = qt.QLineEdit()
        self.user_note.setPlaceholderText("Optional note stored with the job")
        self._add_form_row(
            form,
            "User note:",
            self.user_note,
            "Optional free-text note embedded in the job descriptor and output.",
        )
        layout.addWidget(metadata_group)

        normalization_group = qt.QGroupBox("Exposure and monitor normalization")
        normalization_form = qt.QFormLayout(normalization_group)
        self.normalize_exposure = qt.QCheckBox("Normalize by exposure time")
        self.normalize_exposure.setChecked(
            bool(self.orgui.reconstruction_normalize_exposure)
        )
        self.normalize_exposure.toggled.connect(
            self._on_normalize_exposure_changed
        )
        self._add_form_row(
            normalization_form,
            "",
            self.normalize_exposure,
            "Divide each frame by its exposure time when the scan backend "
            "provides one. Applies only to reciprocal-space reconstruction.",
        )
        self.monitor_corrections = qt.QLineEdit()
        self.monitor_corrections.setPlaceholderText(
            "Optional scan counters, comma-separated"
        )
        self.monitor_corrections.setText(
            ", ".join(self.orgui.reconstruction_monitor_corrections)
        )
        self.monitor_corrections.editingFinished.connect(
            self._on_monitor_corrections_changed
        )
        self._add_form_row(
            normalization_form,
            "Monitor corrections:",
            self.monitor_corrections,
            "Counters applied as divisive monitor normalizations. Applies "
            "only to reciprocal-space reconstruction.",
        )
        layout.addWidget(normalization_group)
        return widget

    def _set_integration_option(self, key, checked):
        """Update one correction shared with the main integration panel."""
        selector = getattr(self.orgui, "scanSelector", None)
        if selector is not None:
            selector.set_integration_options({key: bool(checked)})

    def _sync_integration_options(self):
        """Refresh shared correction switches from the integration panel."""
        selector = getattr(self.orgui, "scanSelector", None)
        if selector is None:
            return
        options = selector.get_integration_options()
        for control, key in (
            (self.use_pixel_mask, "mask"),
            (self.use_solid_angle, "solidAngle"),
            (self.use_polarization, "polarization"),
        ):
            with qt.QSignalBlocker(control):
                control.setChecked(bool(options.get(key, False)))

    def _on_normalize_exposure_changed(self, checked):
        self.orgui.reconstruction_normalize_exposure = bool(checked)

    def _on_monitor_corrections_changed(self):
        self.orgui.reconstruction_monitor_corrections = tuple(
            value.strip()
            for value in self.monitor_corrections.text().split(",")
            if value.strip()
        )

    def _grid_tab(self):
        widget = qt.QWidget()
        layout = qt.QVBoxLayout(widget)
        grid_group = qt.QGroupBox("Grid definitions")
        grid_layout = qt.QVBoxLayout(grid_group)
        self.grid_table = qt.QTableWidget(0, len(_GRID_COLUMNS))
        self.grid_table.setHorizontalHeaderLabels(_GRID_COLUMNS)
        self.grid_table.setToolTip(
            "Edit coordinate bounds and voxel spacing for each output "
            "coordinate system."
        )
        header_tooltips = (
            "Unique output-grid name.",
            "HKL or momentum-transfer reference frame.",
            "Lower bound of coordinate axis 1.",
            "Lower bound of coordinate axis 2.",
            "Lower bound of coordinate axis 3.",
            "Upper bound of coordinate axis 1.",
            "Upper bound of coordinate axis 2.",
            "Upper bound of coordinate axis 3.",
            "Voxel spacing along coordinate axis 1.",
            "Voxel spacing along coordinate axis 2.",
            "Voxel spacing along coordinate axis 3.",
            "Editable number of voxel intervals along coordinate axis 1.",
            "Editable number of voxel intervals along coordinate axis 2.",
            "Editable number of voxel intervals along coordinate axis 3.",
            "Estimated dense, uncompressed size of intensity, variance, "
            "weight, and contributor datasets. Sparse allocation and "
            "compression normally reduce the final file size.",
        )
        for column, tooltip in enumerate(header_tooltips):
            self.grid_table.horizontalHeaderItem(column).setToolTip(tooltip)
        header = self.grid_table.horizontalHeader()
        header.setSectionResizeMode(qt.QHeaderView.Interactive)
        for column in range(len(_GRID_COLUMNS)):
            header.resizeSection(
                column,
                130 if column in {0, 14} else 82,
            )
        header.setStretchLastSection(False)
        self.grid_table.cellChanged.connect(self._on_grid_cell_changed)
        grid_layout.addWidget(self.grid_table)
        buttons = qt.QVBoxLayout()
        grid_definition_buttons = qt.QHBoxLayout()
        grid_setting_buttons = qt.QHBoxLayout()
        add_hkl = qt.QPushButton("Add derived HKL grid")
        add_hkl.setToolTip(
            "Estimate editable HKL bounds and spacing from the active scan."
        )
        add_hkl.clicked.connect(lambda: self.add_derived_grid("hkl"))
        add_q = qt.QPushButton("Add derived Q grid")
        add_q.setToolTip(
            "Add an editable momentum-transfer grid in a selected frame."
        )
        add_q.clicked.connect(self._add_q_grid)
        remove = qt.QPushButton("Remove selected grid")
        remove.setToolTip("Remove every grid row containing a selected cell.")
        remove.clicked.connect(self._remove_grid)
        estimate_steps = qt.QPushButton("Estimate geometry-matched steps")
        estimate_steps.setToolTip(
            "Estimate axis steps from local detector and scan geometry "
            "Jacobians, then apply them to selected grids or all grids."
        )
        estimate_steps.clicked.connect(self._estimate_geometry_steps)
        hdf5_settings = qt.QPushButton("HDF5 settings")
        hdf5_settings.setToolTip(
            "Set one chunk shape for all output grids and optionally override "
            "the active database compression."
        )
        hdf5_settings.clicked.connect(self._edit_hdf5_settings)
        grid_definition_buttons.addWidget(add_hkl)
        grid_definition_buttons.addWidget(add_q)
        grid_definition_buttons.addWidget(remove)
        grid_definition_buttons.addStretch(1)
        grid_setting_buttons.addWidget(estimate_steps)
        grid_setting_buttons.addWidget(hdf5_settings)
        grid_setting_buttons.addStretch(1)
        buttons.addLayout(grid_definition_buttons)
        buttons.addLayout(grid_setting_buttons)
        grid_layout.addLayout(buttons)
        self.hdf5_summary = qt.QLabel()
        self.hdf5_summary.setWordWrap(True)
        self.hdf5_summary.setToolTip(
            "Shared HDF5 chunk shape and compression selection for every grid."
        )
        grid_layout.addWidget(self.hdf5_summary)
        self._refresh_hdf5_summary()
        self.file_count_summary = qt.QLabel()
        self.file_count_summary.setWordWrap(True)
        self.file_count_summary.setToolTip(
            "Estimated checkpoint scratch-file and cluster-job counts, from "
            "a short live calibration probe against the active scan and "
            "detector geometry (design doc Sec6/Sec14). An estimate, not a "
            "guarantee -- the checkpoint-count floor and available memory "
            "may still produce more files than shown here."
        )
        grid_layout.addWidget(self.file_count_summary)
        layout.addWidget(grid_group)
        return widget

    def _refresh_hdf5_summary(self):
        chunk = self.orgui.reconstruction_chunk_shape
        override = self.orgui.reconstruction_compression_override
        if override is None:
            database = getattr(self.orgui, "database", None)
            compression = (
                f"active database filter "
                f"({compression_filter_name(database.compression)})"
                if database is not None
                else "active database filter"
            )
        else:
            compression = f"override ({override})"
        self.hdf5_summary.setText(
            f"All grids: chunks {chunk[0]} × {chunk[1]} × {chunk[2]} voxels; "
            f"compression: {compression}."
        )

    def _edit_hdf5_settings(self):
        try:
            database = self.orgui.database
            dialog = HDF5SettingsDialog(
                database.compression,
                self.orgui.reconstruction_chunk_shape,
                self.orgui.reconstruction_compression_override,
                self,
            )
            if dialog.exec() != qt.QDialog.Accepted:
                return
            database.compression = FILTERS[dialog.database_compression_name]
            self.orgui.reconstruction_chunk_shape = dialog.chunk_shape
            self.orgui.reconstruction_compression_override = (
                dialog.compression_override
            )
            self._refresh_hdf5_summary()
        except Exception as error:
            self._report_failure("Cannot update HDF5 settings", error)

    def _refresh_file_count_summary(self):
        """Live-update the checkpoint file/job count estimate (design doc
        Sec14). Never raises -- called reactively from several widgets
        across tabs, not from a button slot, so failures are reported
        in-place instead of via :meth:`_report_failure`.
        """
        try:
            if self.orgui.fscan is None:
                self.file_count_summary.setText(
                    "Load a scan to estimate checkpoint file and job counts."
                )
                return
            if self.grid_table.rowCount() == 0:
                self.file_count_summary.setText(
                    "Add an output grid to estimate checkpoint file and job "
                    "counts."
                )
                return
            grids = self._grids()
            config = ConfigData.from_gui(self.orgui)
            scan = self.orgui.fscan
            accuracy, advanced_depth = self._accuracy_settings()
            depth = (
                ACCURACY_DEPTHS[accuracy]
                if advanced_depth is None
                else advanced_depth
            )
            threads = (
                self._optional_value(self.thread_override)
                or self.orgui.numberthreads
            )
            ram_budget_bytes = (
                self._optional_value(self.memory_override)
                or self.orgui.maxMemory
            ) * 1024 * 1024
            checkpoint_count = self.checkpoint_count.value()
            angle_fallback = self.angle_fallback.currentData()

            single = estimate_checkpoint_plan(
                config,
                scan,
                grids,
                max_depth=depth,
                threads=threads,
                ram_budget_bytes=ram_budget_bytes,
                checkpoint_count=checkpoint_count,
                angle_fallback=angle_fallback,
                budget_seconds=0.05,
            )
            lines = [
                f"Estimated checkpoint files (single node): "
                f"{single['files_total']}."
            ]
            if len(single["per_grid"]) > 1:
                lines.append(
                    "  Per grid: "
                    + ", ".join(
                        f"{name}: {result['files_per_job']} files "
                        f"(~{_format_size(result['job_data_bytes_estimate'])})"
                        for name, result in single["per_grid"].items()
                    )
                )

            array_task_count = self.cluster_array_task_count.value()
            if array_task_count > 1:
                node_excluded = _node_excluded_frames(
                    len(scan),
                    config.corrections.excluded_frames,
                    total_tasks=array_task_count,
                    task_index=0,
                )
                cluster = estimate_checkpoint_plan(
                    config,
                    scan,
                    grids,
                    max_depth=depth,
                    threads=threads,
                    ram_budget_bytes=ram_budget_bytes,
                    checkpoint_count=checkpoint_count,
                    angle_fallback=angle_fallback,
                    extra_excluded_frames=node_excluded,
                    budget_seconds=0.05,
                )
                lines.append(
                    f"If run as a {array_task_count}-node cluster job: "
                    f"~{cluster['files_total'] * array_task_count} checkpoint "
                    f"files total across {array_task_count} map jobs "
                    f"(~{cluster['files_total']} per node)."
                )
            self.file_count_summary.setText("\n".join(lines))
        except Exception as error:
            self.file_count_summary.setText(
                f"Cannot estimate checkpoint file counts: {error}"
            )

    def _estimate_geometry_steps(self):
        try:
            if self.orgui.fscan is None:
                raise ValueError(
                    "Load a scan before estimating geometry-matched steps"
                )
            rows = sorted(
                {index.row() for index in self.grid_table.selectedIndexes()}
            )
            if not rows:
                rows = list(range(self.grid_table.rowCount()))
            if not rows:
                raise ValueError("Add an output grid before estimating steps")
            dialog = _GeometryResolutionDialog(self)
            if dialog.exec() != qt.QDialog.Accepted:
                return
            self._apply_geometry_steps(rows, dialog.percentile)
        except Exception as error:
            self._report_failure(
                "Cannot estimate geometry-matched grid steps", error
            )

    def _apply_geometry_steps(self, rows, percentile):
        config = ConfigData.from_gui(self.orgui)
        estimates = {}
        with qt.QSignalBlocker(self.grid_table):
            for row in rows:
                frame = self.grid_table.item(row, 1).text().strip()
                if frame not in estimates:
                    estimates[frame] = estimate_geometry_steps(
                        config,
                        self.orgui.fscan,
                        frame=frame,
                        percentile=percentile,
                    )
                for axis, step in enumerate(estimates[frame]):
                    self.grid_table.item(row, 8 + axis).setData(
                        qt.Qt.EditRole, step
                    )
                self._update_grid_row(row)

    def _performance_tab(self):
        widget = qt.QWidget()
        layout = qt.QVBoxLayout(widget)

        accuracy_group = qt.QGroupBox("Accuracy")
        accuracy_form = qt.QFormLayout(accuracy_group)
        self.accuracy = qt.QComboBox()
        self.accuracy.addItem("Center only (depth 0)", "center")
        self.accuracy.addItem("Low (depth 1)", "low")
        self.accuracy.addItem("Balanced (depth 2)", "balanced")
        self.accuracy.addItem("High (depth 3)", "high")
        self.accuracy.addItem("Very high (depth 4)", "very_high")
        self.accuracy.addItem("Maximum (depth 5)", "maximum")
        self.accuracy.setSizeAdjustPolicy(
            qt.QComboBox.AdjustToMinimumContentsLengthWithIcon
        )
        self.accuracy.setMinimumContentsLength(24)
        self.accuracy.setCurrentIndex(
            self.accuracy.findData("balanced")
        )
        self._add_form_row(
            accuracy_form,
            "Footprint preset:",
            self.accuracy,
            "Select the adaptive pixel-footprint subdivision depth. Higher "
            "depths resolve voxel boundaries more accurately but require "
            "substantially more computation.",
        )
        layout.addWidget(accuracy_group)

        execution_group = qt.QGroupBox("Parallel execution and memory")
        execution_form = qt.QFormLayout(execution_group)
        thread_tooltip = (
            "Total CPU-thread budget shared by concurrent images and the "
            "native threads working on each image."
        )
        self.thread_override = self._optional_spin(
            1, 4096, thread_tooltip
        )
        self._add_form_row(
            execution_form,
            "Total thread budget:",
            self.thread_override[0],
            thread_tooltip,
        )
        threads_per_image_tooltip = (
            "Native C++ threads assigned to one image; concurrent image "
            "workers use the remaining total thread and memory budgets. "
            "Unchecked (default): chosen automatically and rebalanced "
            "live during mapping against the real measured image "
            "delivery rate. Checked: pinned to a fixed value for the "
            "whole run."
        )
        self.threads_per_image = self._optional_spin(
            1, 4096, threads_per_image_tooltip
        )
        self._add_form_row(
            execution_form,
            "Native threads per image:",
            self.threads_per_image[0],
            threads_per_image_tooltip,
        )
        memory_tooltip = (
            "Maximum total RAM available to this reconstruction job."
        )
        self.memory_override = self._optional_spin(
            1, 1024 * 1024, memory_tooltip, suffix=" MiB"
        )
        self._add_form_row(
            execution_form,
            "Total memory budget:",
            self.memory_override[0],
            memory_tooltip,
        )
        accumulation_tooltip = (
            "Maximum reduced-record RAM retained by each image worker before "
            "it writes a larger Parquet segment. The global memory budget "
            "also reserves space for sorting and Arrow conversion."
        )
        self.accumulation_memory = self._optional_spin(
            1, 1024 * 1024, accumulation_tooltip, suffix=" MiB"
        )
        self._add_form_row(
            execution_form,
            "Accumulation per worker:",
            self.accumulation_memory[0],
            accumulation_tooltip,
        )
        layout.addWidget(execution_group)

        advanced_group = qt.QGroupBox("Advanced settings")
        advanced_form = qt.QFormLayout(advanced_group)
        frame_tooltip = (
            "Number of consecutive scan frames assigned to each map task."
        )
        self.frame_batch = self._optional_spin(1, 100000, frame_tooltip)
        self._add_form_row(
            advanced_form,
            "Frames per task:",
            self.frame_batch[0],
            frame_tooltip,
        )
        tile_tooltip = (
            "Rectangular detector tiles may end at any pixel boundary; "
            "neighboring tiles share the same detector corner rays."
        )
        self.tile_rows = self._optional_spin(1, 100000, tile_tooltip)
        self._add_form_row(
            advanced_form, "Tile rows:", self.tile_rows[0], tile_tooltip
        )
        self.tile_columns = self._optional_spin(1, 100000, tile_tooltip)
        self._add_form_row(
            advanced_form,
            "Tile columns:",
            self.tile_columns[0],
            tile_tooltip,
        )
        block_tooltip = (
            "Selects the starting cache-scale for one native work block. "
            "A block's working set is the pixel data it reads plus the "
            "accumulator it builds, and throughput is best when that fits "
            "the CPU cache: 'Minimum' is about an L1 data cache, 'Medium' "
            "about a 1 MiB L2. The final pixel count is calculated "
            "automatically from the adaptive depth and capped to fit the "
            "memory budget."
        )
        self.work_block = qt.QComboBox()
        self.work_block.setSizeAdjustPolicy(
            qt.QComboBox.AdjustToMinimumContentsLengthWithIcon
        )
        self.work_block.setMinimumContentsLength(20)
        for name, pixels in sorted(
            WORK_BLOCK_PRESETS.items(), key=lambda item: item[1]
        ):
            self.work_block.addItem(
                f"{name.capitalize()} ({pixels:,} px at center accuracy)", name
            )
        self.work_block.setCurrentIndex(self.work_block.findData("medium"))
        self._add_form_row(
            advanced_form,
            "Native block pixels:",
            self.work_block,
            block_tooltip,
        )
        checkpoint_tooltip = (
            "Minimum number of resumable HDF5 checkpoint files per output "
            "grid. A floor, not a target: the actual count exceeds it when "
            "the estimated data volume would not otherwise fit the memory "
            "budget across that many files."
        )
        self.checkpoint_count = qt.QSpinBox()
        self.checkpoint_count.setRange(1, 100000)
        self.checkpoint_count.setValue(10)
        self.checkpoint_count.setToolTip(checkpoint_tooltip)
        self._add_form_row(
            advanced_form,
            "Checkpoint count:",
            self.checkpoint_count,
            checkpoint_tooltip,
        )
        layout.addWidget(advanced_group)

        detected_group = qt.QGroupBox("Detected execution layout")
        detected_layout = qt.QVBoxLayout(detected_group)
        self.performance_summary = qt.QPlainTextEdit()
        self.performance_summary.setReadOnly(True)
        self.performance_summary.setMaximumHeight(190)
        self.performance_summary.setToolTip(
            "Resolved task, tiling, thread, and memory values for this job."
        )
        detected_layout.addWidget(self.performance_summary)
        note = qt.QLabel(
            "Unset values are derived from the active CPU count, orGUI memory "
            "budget, detector shape, and output chunk geometry. Unchecked "
            "editors display the detected values; enable Override to freeze "
            "an edited value into this job."
        )
        note.setWordWrap(True)
        note.setToolTip(
            "Enable Override beside a setting only when its detected value "
            "should be replaced for this job."
        )
        detected_layout.addWidget(note)
        layout.addWidget(detected_group)
        layout.addStretch(1)
        return widget

    def _paths_tab(self):
        widget = qt.QWidget()
        layout = qt.QVBoxLayout(widget)
        self.job_path = qt.QLineEdit()
        self.scratch_path = qt.QLineEdit()
        self.output_path = qt.QLineEdit()

        job_group = qt.QGroupBox("Job descriptor")
        job_form = qt.QFormLayout(job_group)
        job_tooltip = (
            "Resumable JSON descriptor containing the frozen experiment and "
            "execution settings."
        )
        self._add_form_row(
            job_form,
            "Job JSON:",
            self._path_row(
                self.job_path,
                "job",
                save=True,
                tooltip=job_tooltip,
            ),
            job_tooltip,
        )
        layout.addWidget(job_group)

        storage_group = qt.QGroupBox("Storage")
        storage_form = qt.QFormLayout(storage_group)
        scratch_tooltip = (
            "Writable directory for immutable Parquet partitions and the "
            "checksummed job asset bundle."
        )
        self._add_form_row(
            storage_form,
            "Scratch directory:",
            self._path_row(
                self.scratch_path,
                "directory",
                tooltip=scratch_tooltip,
            ),
            scratch_tooltip,
        )
        output_tooltip = (
            "Final standalone NeXus/HDF5 reciprocal-space reconstruction file."
        )
        self._add_form_row(
            storage_form,
            "Standalone HDF5:",
            self._path_row(
                self.output_path,
                "hdf5",
                save=True,
                tooltip=output_tooltip,
            ),
            output_tooltip,
        )
        note = qt.QLabel(
            "Scratch data are retained after interruption and removed "
            "automatically only after successful HDF5 validation."
        )
        note.setWordWrap(True)
        note.setToolTip(
            "Interrupted jobs keep their scratch data so they can be resumed."
        )
        storage_form.addRow(note)
        layout.addWidget(storage_group)
        layout.addStretch(1)
        return widget

    def _cluster_tab(self):
        widget = qt.QWidget()
        outer = qt.QVBoxLayout(widget)
        scroll = qt.QScrollArea()
        scroll.setWidgetResizable(True)
        contents = qt.QWidget()
        layout = qt.QVBoxLayout(contents)

        scheduler_group = qt.QGroupBox("Scheduler")
        scheduler_form = qt.QFormLayout(scheduler_group)
        self.cluster_scheduler = qt.QComboBox()
        self.cluster_scheduler.addItem("Sun/Grid Engine (SGE)", "sge")
        self.cluster_scheduler.addItem("Slurm", "slurm")
        self._add_form_row(
            scheduler_form,
            "Scheduler:",
            self.cluster_scheduler,
            "Generate SGE qsub scripts or Slurm sbatch scripts. SGE is the "
            "default and primary target.",
        )
        self.cluster_job_name = qt.QLineEdit("orgui-rsmap")
        self._add_form_row(
            scheduler_form,
            "Job name:",
            self.cluster_job_name,
            "Scheduler-safe base name used for the map array and finalizer.",
        )
        self.cluster_queue = qt.QLineEdit()
        self._add_form_row(
            scheduler_form,
            "Queue / partition:",
            self.cluster_queue,
            "Optional SGE queue (-q) or Slurm partition (--partition).",
        )
        self.cluster_account = qt.QLineEdit()
        self._add_form_row(
            scheduler_form,
            "Project / account:",
            self.cluster_account,
            "Optional SGE project (-P) or Slurm account (--account).",
        )
        layout.addWidget(scheduler_group)

        environment_group = qt.QGroupBox("Python environment")
        environment_form = qt.QFormLayout(environment_group)
        self.cluster_script_directory = qt.QLineEdit()
        self._add_form_row(
            environment_form,
            "Script directory:",
            self._path_row(
                self.cluster_script_directory,
                "directory",
                tooltip=(
                    "Directory receiving map, finalizer, and submission scripts."
                ),
            ),
            "Directory receiving map, finalizer, and submission scripts.",
        )
        self.cluster_working_directory = qt.QLineEdit(str(Path.cwd()))
        self._add_form_row(
            environment_form,
            "Working directory:",
            self._path_row(
                self.cluster_working_directory,
                "directory",
                tooltip="Shared working directory selected before Python starts.",
            ),
            "Shared working directory selected before Python starts.",
        )
        self.cluster_python = qt.QLineEdit("python")
        self._add_form_row(
            environment_form,
            "Python executable:",
            self.cluster_python,
            "Python command available after environment setup, normally "
            "'python' from the activated orGUI environment.",
        )
        self.cluster_environment = qt.QPlainTextEdit()
        self.cluster_environment.setPlaceholderText(
            "module load ...\nsource /path/to/venv/bin/activate"
        )
        self.cluster_environment.setMaximumHeight(90)
        self._add_form_row(
            environment_form,
            "Setup commands:",
            self.cluster_environment,
            "Shell commands run before every task, for example module loads "
            "and conda or virtual-environment activation.",
        )
        layout.addWidget(environment_group)

        map_group = qt.QGroupBox("Mapping array")
        map_form = qt.QFormLayout(map_group)
        self.cluster_array_task_count = qt.QSpinBox()
        self.cluster_array_task_count.setRange(1, 100000)
        self.cluster_array_task_count.setValue(4)
        self._add_form_row(
            map_form,
            "Number of nodes:",
            self.cluster_array_task_count,
            "How many array elements (nodes) to request. Each node maps an "
            "equal, disjoint share of the scan's frames -- the array size "
            "is never inferred from the scheduler at run time, since most "
            "schedulers besides Slurm do not expose it to a running task.",
        )
        self.cluster_array_cpus = qt.QSpinBox()
        self.cluster_array_cpus.setRange(1, 4096)
        self.cluster_array_cpus.setValue(4)
        self._add_form_row(
            map_form,
            "CPUs / slots per task:",
            self.cluster_array_cpus,
            "Native C++ threads used by each independent map-array task.",
        )
        self.cluster_array_memory = qt.QDoubleSpinBox()
        self.cluster_array_memory.setRange(0.25, 1024 * 1024)
        self.cluster_array_memory.setDecimals(2)
        self.cluster_array_memory.setValue(16.0)
        self.cluster_array_memory.setSuffix(" GiB")
        self._add_form_row(
            map_form,
            "Memory per task:",
            self.cluster_array_memory,
            "Total RAM budget for one array task. For SGE, the generated "
            "h_vmem request is divided across the requested slots.",
        )
        self.cluster_array_walltime = qt.QLineEdit("24:00:00")
        self._add_form_row(
            map_form,
            "Wall time:",
            self.cluster_array_walltime,
            "Scheduler wall-time limit for each mapping task.",
        )
        self.cluster_array_concurrency = qt.QSpinBox()
        self.cluster_array_concurrency.setRange(0, 1000000)
        self.cluster_array_concurrency.setValue(0)
        self._add_form_row(
            map_form,
            "Maximum concurrent tasks:",
            self.cluster_array_concurrency,
            "Optional array throttle (-tc or %N). Zero leaves the scheduler "
            "limit unchanged.",
        )
        self.cluster_summary = qt.QPlainTextEdit()
        self.cluster_summary.setReadOnly(True)
        self.cluster_summary.setMaximumHeight(75)
        self._add_form_row(
            map_form,
            "Detected array:",
            self.cluster_summary,
            "One array element is created for each detected deterministic "
            "map task. Prepare the job to update this count.",
        )
        layout.addWidget(map_group)

        reduce_group = qt.QGroupBox("Reduction and finalization")
        reduce_form = qt.QFormLayout(reduce_group)
        self.cluster_reduce_cpus = qt.QSpinBox()
        self.cluster_reduce_cpus.setRange(1, 4096)
        self.cluster_reduce_cpus.setValue(24)
        self._add_form_row(
            reduce_form,
            "CPUs / slots:",
            self.cluster_reduce_cpus,
            "Independent CPU count for parallel shard reduction; it need not "
            "match the map-array task size.",
        )
        self.cluster_reduce_memory = qt.QDoubleSpinBox()
        self.cluster_reduce_memory.setRange(0.25, 1024 * 1024)
        self.cluster_reduce_memory.setDecimals(2)
        self.cluster_reduce_memory.setValue(64.0)
        self.cluster_reduce_memory.setSuffix(" GiB")
        self._add_form_row(
            reduce_form,
            "Memory:",
            self.cluster_reduce_memory,
            "Total RAM budget divided among parallel reducer workers.",
        )
        self.cluster_reduce_walltime = qt.QLineEdit("24:00:00")
        self._add_form_row(
            reduce_form,
            "Wall time:",
            self.cluster_reduce_walltime,
            "Wall-time limit for reduction and final HDF5 creation.",
        )
        layout.addWidget(reduce_group)

        advanced_group = qt.QGroupBox("Scheduler-specific settings")
        advanced_form = qt.QFormLayout(advanced_group)
        self.cluster_sge_pe = qt.QLineEdit("smp")
        self._add_form_row(
            advanced_form,
            "SGE parallel environment:",
            self.cluster_sge_pe,
            "SGE parallel environment requested with -pe, commonly 'smp'.",
        )
        self.cluster_sge_memory = qt.QLineEdit("h_vmem")
        self._add_form_row(
            advanced_form,
            "SGE memory resource:",
            self.cluster_sge_memory,
            "SGE per-slot complex used for memory requests, commonly h_vmem "
            "or mem_free.",
        )
        self.cluster_array_directives = qt.QPlainTextEdit()
        self.cluster_array_directives.setMaximumHeight(70)
        self._add_form_row(
            advanced_form,
            "Extra map directives:",
            self.cluster_array_directives,
            "Optional scheduler directives, one '#$' or '#SBATCH' line each.",
        )
        self.cluster_reduce_directives = qt.QPlainTextEdit()
        self.cluster_reduce_directives.setMaximumHeight(70)
        self._add_form_row(
            advanced_form,
            "Extra finalizer directives:",
            self.cluster_reduce_directives,
            "Optional scheduler directives for reduction/finalization.",
        )
        layout.addWidget(advanced_group)
        layout.addStretch(1)
        scroll.setWidget(contents)
        outer.addWidget(scroll)
        return widget

    def _cluster_settings(self):
        return ClusterSettings(
            scheduler=self.cluster_scheduler.currentData(),
            job_name=self.cluster_job_name.text().strip(),
            script_directory=self.cluster_script_directory.text().strip(),
            working_directory=self.cluster_working_directory.text().strip(),
            python_executable=self.cluster_python.text().strip(),
            environment_setup=self.cluster_environment.toPlainText(),
            queue=self.cluster_queue.text().strip(),
            account=self.cluster_account.text().strip(),
            array_task_count=self.cluster_array_task_count.value(),
            array_cpus=self.cluster_array_cpus.value(),
            array_memory_gib=self.cluster_array_memory.value(),
            array_walltime=self.cluster_array_walltime.text().strip(),
            array_concurrency=self.cluster_array_concurrency.value(),
            reduce_cpus=self.cluster_reduce_cpus.value(),
            reduce_memory_gib=self.cluster_reduce_memory.value(),
            reduce_walltime=self.cluster_reduce_walltime.text().strip(),
            sge_parallel_environment=self.cluster_sge_pe.text().strip(),
            sge_memory_resource=self.cluster_sge_memory.text().strip(),
            extra_array_directives=(
                self.cluster_array_directives.toPlainText()
            ),
            extra_reduce_directives=(
                self.cluster_reduce_directives.toPlainText()
            ),
        )

    def _set_cluster_settings(self, values):
        settings = ClusterSettings.from_dict(values)
        self.cluster_scheduler.setCurrentIndex(
            self.cluster_scheduler.findData(settings.scheduler)
        )
        self.cluster_job_name.setText(settings.job_name)
        self.cluster_script_directory.setText(settings.script_directory)
        self.cluster_working_directory.setText(settings.working_directory)
        self.cluster_python.setText(settings.python_executable)
        self.cluster_environment.setPlainText(settings.environment_setup)
        self.cluster_queue.setText(settings.queue)
        self.cluster_account.setText(settings.account)
        self.cluster_array_task_count.setValue(settings.array_task_count)
        self.cluster_array_cpus.setValue(settings.array_cpus)
        self.cluster_array_memory.setValue(settings.array_memory_gib)
        self.cluster_array_walltime.setText(settings.array_walltime)
        self.cluster_array_concurrency.setValue(
            settings.array_concurrency
        )
        self.cluster_reduce_cpus.setValue(settings.reduce_cpus)
        self.cluster_reduce_memory.setValue(settings.reduce_memory_gib)
        self.cluster_reduce_walltime.setText(settings.reduce_walltime)
        self.cluster_sge_pe.setText(settings.sge_parallel_environment)
        self.cluster_sge_memory.setText(settings.sge_memory_resource)
        self.cluster_array_directives.setPlainText(
            settings.extra_array_directives
        )
        self.cluster_reduce_directives.setPlainText(
            settings.extra_reduce_directives
        )

    @staticmethod
    def _set_control_tooltip(control, tooltip):
        widgets = control if isinstance(control, tuple) else (control,)
        for child in widgets:
            child.setToolTip(tooltip)

    def _add_form_row(self, form, label, control, tooltip):
        label_widget = qt.QLabel(label)
        label_widget.setToolTip(tooltip)
        self._set_control_tooltip(control, tooltip)
        form.addRow(label_widget, control)

    def _optional_spin(self, minimum, maximum, tooltip, suffix=""):
        container = qt.QWidget()
        layout = qt.QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        enabled = qt.QCheckBox("Override")
        editor = qt.QSpinBox()
        editor.setRange(minimum, maximum)
        editor.setSuffix(suffix)
        editor.setEnabled(False)
        enabled.toggled.connect(editor.setEnabled)
        layout.addWidget(enabled)
        layout.addWidget(editor)
        layout.addStretch(1)
        control = (container, enabled, editor)
        self._set_control_tooltip(control, tooltip)
        return control

    def _optional_combo(self, items, tooltip, current=None):
        """An override checkbox in front of a named-choice combo.

        Same ``(container, checkbox, editor)`` shape as
        :meth:`_optional_spin`, so the same value and detected-value
        helpers handle both.
        """
        container = qt.QWidget()
        layout = qt.QHBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        enabled = qt.QCheckBox("Override")
        editor = qt.QComboBox()
        editor.setSizeAdjustPolicy(
            qt.QComboBox.AdjustToMinimumContentsLengthWithIcon
        )
        editor.setMinimumContentsLength(20)
        for label, value in items:
            editor.addItem(label, value)
        if current is not None:
            index = editor.findData(current)
            if index >= 0:
                editor.setCurrentIndex(index)
        editor.setEnabled(False)
        enabled.toggled.connect(editor.setEnabled)
        layout.addWidget(enabled)
        layout.addWidget(editor)
        layout.addStretch(1)
        control = (container, enabled, editor)
        self._set_control_tooltip(control, tooltip)
        return control

    def _path_row(self, editor, kind, save=False, tooltip=""):
        widget = qt.QWidget()
        layout = qt.QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        widget.setToolTip(tooltip)
        editor.setToolTip(tooltip)
        layout.addWidget(editor)
        button = qt.QPushButton("Browse…")
        button.setToolTip(tooltip)
        button.clicked.connect(
            lambda: self._browse(editor, kind=kind, save=save)
        )
        layout.addWidget(button)
        return widget

    def _browse(self, editor, *, kind, save):
        try:
            directory = str(self._history_directory())
            if kind == "directory":
                value = qt.QFileDialog.getExistingDirectory(
                    self, "Select scratch directory", directory
                )
            else:
                file_filter = (
                    "JSON job (*.json)"
                    if kind == "job"
                    else "HDF5/NeXus (*.h5 *.hdf5 *.nxs)"
                )
                dialog = (
                    qt.QFileDialog.getSaveFileName
                    if save
                    else qt.QFileDialog.getOpenFileName
                )
                value, _ = dialog(self, "Select file", directory, file_filter)
            if value:
                editor.setText(value)
                self.orgui.filedialogdir = str(Path(value).parent)
                return True
            return False
        except Exception as error:
            self._report_failure("Cannot select reconstruction path", error)
            return False

    def _history_directory(self):
        history = Path(str(self.orgui.filedialogdir)).expanduser()
        try:
            if history.is_dir():
                return history
            if history.parent.is_dir():
                return history.parent
        except OSError:
            pass
        return Path.cwd()

    def _set_default_paths(self, stem):
        base = Path.cwd()
        if not self.job_path.text():
            self.job_path.setText(str(base / f"{stem}-rsmap.json"))
        if not self.scratch_path.text():
            self.scratch_path.setText(str(base / f".{stem}-rsmap-scratch"))
        if not self.output_path.text():
            self.output_path.setText(str(base / f"{stem}-rsmap.h5"))
        if not self.cluster_script_directory.text():
            self.cluster_script_directory.setText(
                str(base / f"{stem}-rsmap-cluster")
            )

    def _show_output(self, description):
        self.preview_output.setPlainText(description)
        self.tabs.setCurrentWidget(self.output_tab)

    def _report_message(self, title, description):
        self._show_output(description)
        logger.warning(
            title,
            extra={
                "title": title,
                "description": description,
                "show_dialog": True,
                "dialog_level": logging.WARNING,
                "parent": self,
            },
        )

    def _report_failure(self, title, error):
        description = str(error) or type(error).__name__
        self._show_output(description)
        logger.warning(
            title,
            exc_info=True,
            extra={
                "title": title,
                "description": description,
                "show_dialog": True,
                "dialog_level": logging.WARNING,
                "parent": self,
            },
        )

    def _invoke_ui_action(self, title, action):
        try:
            action()
        except Exception as error:
            self._report_failure(title, error)

    def _live_correction_summary(self, config, detector_shape):
        correction = config.corrections
        summary = correction.to_dict()
        for asset_field in (
            "mask_asset",
            "background_asset",
            "background_variance_asset",
        ):
            summary.pop(asset_field, None)

        mask = (
            self.orgui.get_detector_mask(detector_shape)
            if correction.use_mask
            else None
        )
        if not correction.use_mask:
            mask_status = {"status": "disabled"}
        elif mask is None:
            mask_status = {"status": "enabled but unavailable"}
        else:
            mask_status = {
                "status": "active",
                "masked_pixels": int(np.count_nonzero(mask)),
                "job_asset_path": "/mask",
            }

        background = (
            getattr(self.orgui, "background_image", None)
            if correction.use_background
            else None
        )
        if not correction.use_background:
            background_status = {"status": "disabled"}
        elif background is None:
            background_status = {"status": "enabled but unavailable"}
        else:
            background_status = {
                "status": "active",
                "job_asset_path": "/background",
                "variance": (
                    "available"
                    if getattr(self.orgui, "background_variance", None)
                    is not None
                    else "deterministic"
                ),
            }
        summary["active_inputs"] = {
            "mask": mask_status,
            "background": background_status,
        }
        return summary

    def refresh_live_state(self):
        """Refresh the experiment and correction summary from orGUI."""
        try:
            self._sync_integration_options()
            if self.orgui.fscan is None:
                self.experiment_summary.setPlainText("No active scan.")
                return
            config = ConfigData.from_gui(self.orgui)
            detector_shape = self.orgui.ubcalc.detectorCal.detector.shape
            summary = {
                "scan": getattr(
                    self.orgui.fscan,
                    "name",
                    getattr(self.orgui.fscan, "title", ""),
                ),
                "frames": len(self.orgui.fscan),
                "detector_shape": detector_shape,
                "energy_keV": self.orgui.ubcalc.ubCal.getEnergy(),
                "UB": self.orgui.ubcalc.ubCal.getUB().tolist(),
                "corrections": self._live_correction_summary(
                    config, detector_shape
                ),
                "threads": self.orgui.numberthreads,
                "memory_MiB": self.orgui.maxMemory,
            }
            self.experiment_summary.setPlainText(
                json.dumps(summary, indent=2, sort_keys=True)
            )
            stem = "".join(
                character if character.isalnum() or character in "-_" else "_"
                for character in str(summary["scan"])
            ).strip("_") or "scan"
            self._set_default_paths(stem)
            if self.grid_table.rowCount() == 0:
                self.add_derived_grid("hkl")
        except Exception as error:
            self._report_failure("Cannot refresh reconstruction settings", error)

    def _add_q_grid(self):
        try:
            frame, accepted = qt.QInputDialog.getItem(
                self,
                "Q coordinate frame",
                "Frame:",
                list(_FRAMES[1:]),
                0,
                False,
            )
            if accepted:
                self.add_derived_grid(frame)
        except Exception as error:
            self._report_failure("Cannot add Q output grid", error)

    def add_derived_grid(self, frame):
        """Add a grid derived from current scan coverage."""
        if self.orgui.fscan is None:
            self._report_message(
                "No scan loaded",
                "Load a scan before deriving reciprocal-space grid coverage.",
            )
            return
        try:
            config = ConfigData.from_gui(self.orgui)
            grid = derive_grid(config, self.orgui.fscan, frame=frame)
            self._append_grid(grid)
        except Exception as error:
            self._report_failure("Cannot derive output grid", error)

    def _append_grid(self, grid):
        row = self.grid_table.rowCount()
        with qt.QSignalBlocker(self.grid_table):
            self.grid_table.insertRow(row)
            values = [
                grid.name or "",
                grid.frame,
                *grid.minimum,
                *grid.maximum,
                *grid.step,
            ]
            for column, value in enumerate(values):
                item = (
                    qt.QTableWidgetItem(str(value))
                    if column < 2
                    else _GridNumberItem(value)
                )
                self.grid_table.setItem(row, column, item)
            self._update_grid_row(row)
        self._refresh_file_count_summary()

    def _on_grid_cell_changed(self, row, column):
        try:
            self._update_grid_row(row, changed_column=column)
        except (TypeError, ValueError, OverflowError):
            size_item = self.grid_table.item(row, 14)
            if size_item is not None:
                with qt.QSignalBlocker(self.grid_table):
                    size_item.setText("invalid")
        self._refresh_file_count_summary()

    def _update_grid_row(self, row, changed_column=None):
        with qt.QSignalBlocker(self.grid_table):
            intervals = []
            for axis in range(3):
                minimum = float(
                    self.grid_table.item(row, 2 + axis).data(qt.Qt.EditRole)
                )
                maximum = float(
                    self.grid_table.item(row, 5 + axis).data(qt.Qt.EditRole)
                )
                extent = maximum - minimum
                if not np.isfinite(extent) or extent <= 0:
                    raise ValueError("Grid extent must be finite and positive")
                interval_column = 11 + axis
                if changed_column == interval_column:
                    value = float(
                        self.grid_table.item(row, interval_column).data(
                            qt.Qt.EditRole
                        )
                    )
                    count = int(value)
                    if count < 1 or value != count:
                        raise ValueError(
                            "Grid interval counts must be positive integers"
                        )
                    step = np.nextafter(extent / count, np.inf)
                    self.grid_table.item(row, 8 + axis).setData(
                        qt.Qt.EditRole, step
                    )
                else:
                    step = float(
                        self.grid_table.item(row, 8 + axis).data(
                            qt.Qt.EditRole
                        )
                    )
                    if not np.isfinite(step) or step <= 0:
                        raise ValueError(
                            "Grid steps must be finite and positive"
                        )
                    count = int(np.ceil(extent / step))
                    item = self.grid_table.item(row, interval_column)
                    if item is None:
                        item = qt.QTableWidgetItem()
                        self.grid_table.setItem(row, interval_column, item)
                    item.setData(qt.Qt.EditRole, count)
                intervals.append(count)

            estimated_bytes = int(np.prod(intervals, dtype=object)) * 32
            size_item = self.grid_table.item(row, 14)
            if size_item is None:
                size_item = qt.QTableWidgetItem()
                size_item.setFlags(
                    size_item.flags() & ~qt.Qt.ItemIsEditable
                )
                self.grid_table.setItem(row, 14, size_item)
            size_item.setText(f"{_format_size(estimated_bytes)} uncompressed")

    def _remove_grid(self):
        try:
            rows = sorted(
                {index.row() for index in self.grid_table.selectedIndexes()},
                reverse=True,
            )
            for row in rows:
                self.grid_table.removeRow(row)
            self._refresh_file_count_summary()
        except Exception as error:
            self._report_failure("Cannot remove output grid", error)

    def _grids(self):
        grids = []
        for row in range(self.grid_table.rowCount()):
            items = [
                self.grid_table.item(row, column)
                for column in range(self.grid_table.columnCount())
            ]
            if any(item is None for item in items):
                raise ValueError(f"Grid row {row + 1} contains an empty cell")
            values = [item.text().strip() for item in items[:2]]
            if values[1] not in _FRAMES:
                raise ValueError(
                    f"Grid row {row + 1} has an unknown frame: {values[1]}"
                )
            grid = ReconstructionGrid(
                minimum=tuple(
                    float(item.data(qt.Qt.EditRole)) for item in items[2:5]
                ),
                maximum=tuple(
                    float(item.data(qt.Qt.EditRole)) for item in items[5:8]
                ),
                step=tuple(
                    float(item.data(qt.Qt.EditRole)) for item in items[8:11]
                ),
                chunk_shape=tuple(
                    self.orgui.reconstruction_chunk_shape
                ),
                frame=values[1],
                name=values[0] or None,
            )
            intervals = tuple(
                int(item.data(qt.Qt.EditRole)) for item in items[11:14]
            )
            if grid.to_spec().shape != intervals:
                raise ValueError(
                    f"Grid row {row + 1} interval counts are inconsistent "
                    "with its bounds and steps"
                )
            if any(step <= 0 for step in grid.step):
                raise ValueError(f"Grid row {row + 1} steps must be positive")
            if any(
                upper <= lower
                for lower, upper in zip(grid.minimum, grid.maximum)
            ):
                raise ValueError(
                    f"Grid row {row + 1} maxima must exceed its minima"
                )
            grids.append(grid)
        if not grids:
            raise ValueError("At least one output grid is required")
        return grids

    @staticmethod
    def _optional_value(control):
        _, enabled, editor = control
        if not enabled.isChecked():
            return None
        if isinstance(editor, qt.QComboBox):
            return editor.currentData()
        return editor.value()

    @staticmethod
    def _set_detected_value(control, value):
        _, enabled, editor = control
        if enabled.isChecked() or value is None:
            return
        if isinstance(editor, qt.QComboBox):
            # A combo offers named scales; the resolved pixel count this
            # would report belongs in the performance summary instead.
            return
        editor.setValue(
            max(editor.minimum(), min(editor.maximum(), int(round(value))))
        )

    @staticmethod
    def _set_optional_value(control, value):
        """Restore an optional override without losing automatic mode."""
        _, enabled, editor = control
        is_override = value is not None
        if is_override:
            if isinstance(editor, qt.QComboBox):
                index = editor.findData(value)
                if index < 0:
                    editor.addItem(f"Custom ({value})", value)
                    index = editor.count() - 1
                editor.setCurrentIndex(index)
            else:
                editor.setValue(int(round(value)))
        enabled.setChecked(is_override)
        editor.setEnabled(is_override)

    @staticmethod
    def _set_combo_value(combo, value):
        """Select a named value, retaining custom values from saved jobs."""
        index = combo.findData(value)
        if index < 0:
            combo.addItem(f"Custom ({value})", value)
            index = combo.count() - 1
        combo.setCurrentIndex(index)

    def _set_accuracy(self, accuracy, advanced_depth):
        """Restore a preset or an advanced depth from a prepared job."""
        if advanced_depth is None:
            index = self.accuracy.findData(accuracy)
            if index < 0:
                raise ValueError(f"Unknown accuracy preset: {accuracy}")
        else:
            depth = int(advanced_depth)
            if not 0 <= depth <= 8:
                raise ValueError(
                    "Advanced split depth must be between 0 and 8"
                )
            if accuracy not in ACCURACY_DEPTHS:
                raise ValueError(f"Unknown accuracy preset: {accuracy}")
            data = f"advanced:{accuracy}:{depth}"
            index = self.accuracy.findData(data)
            if index < 0:
                self.accuracy.addItem(
                    f"Advanced (depth {depth})", data
                )
                index = self.accuracy.count() - 1
        self.accuracy.setCurrentIndex(index)

    def _accuracy_settings(self):
        """Return the preset and optional advanced depth selected in the UI."""
        accuracy = self.accuracy.currentData()
        if isinstance(accuracy, str) and accuracy.startswith("advanced:"):
            _, preset, depth = accuracy.split(":", 2)
            return preset, int(depth)
        return accuracy, None

    def _show_execution_settings(self, job, *, scan=None, config=None):
        settings = reconstruction_execution_settings(
            job, scan=scan, config=config
        )
        tile_rows, tile_columns = settings["detector_tile_shape"]
        for control, value in (
            (self.thread_override, settings["thread_budget"]),
            (self.memory_override, settings["memory_budget_MiB"]),
            (
                self.accumulation_memory,
                settings["accumulation_budget_MiB_per_worker"],
            ),
            (self.frame_batch, settings["frames_per_task"]),
            (self.tile_rows, tile_rows),
            (self.tile_columns, tile_columns),
            (self.threads_per_image, settings["native_threads_per_image"]),
        ):
            self._set_detected_value(control, value)
        self.performance_summary.setPlainText(
            json.dumps(settings, indent=2, sort_keys=True)
        )
        self.cluster_summary.setPlainText(
            json.dumps(
                {
                    "array_tasks": self.cluster_array_task_count.value(),
                    "cpus_per_task": self.cluster_array_cpus.value(),
                    "maximum_concurrent_tasks": (
                        self.cluster_array_concurrency.value() or "scheduler"
                    ),
                },
                indent=2,
            )
        )
        return settings

    def _prepare(self):
        if self.orgui.fscan is None:
            raise ValueError("Load a scan before preparing a reconstruction job")
        tile_rows = self._optional_value(self.tile_rows)
        tile_columns = self._optional_value(self.tile_columns)
        if (tile_rows is None) != (tile_columns is None):
            raise ValueError("Override both detector tile dimensions together")
        job_path = self.job_path.text().strip()
        scratch_path = self.scratch_path.text().strip()
        output_path = self.output_path.text().strip()
        if not job_path:
            raise ValueError("Select a job JSON path")
        if not scratch_path:
            raise ValueError("Select a scratch directory")
        if not output_path:
            raise ValueError("Select an output HDF5 path")
        accuracy, advanced_depth = self._accuracy_settings()
        job = prepare_job(
            self.orgui,
            job_path,
            grids=self._grids(),
            scratch_path=scratch_path,
            output_path=output_path,
            accuracy=accuracy,
            advanced_depth=advanced_depth,
            compression_override=(
                self.orgui.reconstruction_compression_override
            ),
            angle_fallback=self.angle_fallback.currentData(),
            user_note=self.user_note.text(),
            thread_override=self._optional_value(self.thread_override),
            memory_override_bytes=(
                None
                if self._optional_value(self.memory_override) is None
                else self._optional_value(self.memory_override) * 1024 * 1024
            ),
            frame_batch=self._optional_value(self.frame_batch),
            tile_shape=(
                None
                if tile_rows is None
                else (tile_rows, tile_columns)
            ),
            work_block_pixels=self.work_block.currentData(),
            checkpoint_count=self.checkpoint_count.value(),
            threads_per_image=self._optional_value(self.threads_per_image),
            accumulation_budget_bytes=(
                None
                if self._optional_value(self.accumulation_memory) is None
                else self._optional_value(self.accumulation_memory)
                * 1024
                * 1024
            ),
            cluster_settings=self._cluster_settings().to_dict(),
        )
        self._show_execution_settings(
            job,
            scan=self.orgui.fscan,
            config=job.config_data,
        )
        return job

    @qt.Slot()
    def preview(self):
        """Show live grid and resource estimates without freezing state."""
        try:
            if self.orgui.fscan is None:
                raise ValueError(
                    "Load a scan before previewing a reconstruction job"
                )
            grids = self._grids()
            accuracy, advanced_depth = self._accuracy_settings()
            depth = (
                ACCURACY_DEPTHS[accuracy]
                if advanced_depth is None
                else advanced_depth
            )
            grid_rows = []
            final_bytes = 0
            chunk_count = 0
            for grid in grids:
                shape = tuple(
                    int(np.ceil((upper - lower) / step))
                    for lower, upper, step in zip(
                        grid.minimum, grid.maximum, grid.step
                    )
                )
                voxels = int(np.prod(shape, dtype=np.int64))
                chunks = int(
                    np.prod(
                        [
                            np.ceil(size / chunk)
                            for size, chunk in zip(shape, grid.chunk_shape)
                        ]
                    )
                )
                final_bytes += voxels * 32
                chunk_count += chunks
                grid_rows.append({**grid.__dict__, "shape": shape})
            corrections = ConfigData.from_gui(self.orgui).corrections
            included_frames = len(self.orgui.fscan) - len(
                {
                    frame
                    for frame in corrections.excluded_frames
                    if 0 <= frame < len(self.orgui.fscan)
                }
            )
            detector_pixels = int(
                np.prod(self.orgui.ubcalc.detectorCal.detector.shape)
            )
            result = {
                "grids": grid_rows,
                "frames": included_frames,
                "threads": self._optional_value(self.thread_override)
                or self.orgui.numberthreads,
                "native_threads_per_image": (
                    self._optional_value(self.threads_per_image) or "automatic"
                ),
                "memory_MiB": self._optional_value(self.memory_override)
                or self.orgui.maxMemory,
                "accumulation_MiB_per_worker": self._optional_value(
                    self.accumulation_memory
                )
                or "automatic",
                "estimated_spatial_chunks": chunk_count,
                "uncompressed_final_GiB": final_bytes / 1024**3,
                "footprint_leaf_upper_bound": (
                    included_frames
                    * detector_pixels
                    * min(8**depth, 4096)
                ),
            }
            self._show_output(
                json.dumps(result, indent=2, sort_keys=True)
            )
        except Exception as error:
            self._report_failure("Cannot preview reconstruction", error)

    @qt.Slot()
    def prepare(self):
        """Freeze current state and save a prepared job."""
        try:
            job = self._prepare()
            self._show_output(
                json.dumps(job.to_dict(), indent=2, sort_keys=True)
            )
        except Exception as error:
            self._report_failure("Cannot prepare reconstruction job", error)

    def _run_path(self, path):
        path = str(path).strip()
        if not path:
            raise ValueError("Select a prepared job JSON path")
        progress = None

        def update(value, maximum, message):
            progress.total = maximum
            if hasattr(progress, "dialog"):
                progress.dialog.setMaximum(maximum)
            progress.update(value, message)
            if progress.wasCanceled():
                raise _ReconstructionCancelled(
                    "Reconstruction cancelled between tasks"
                )

        try:
            progress = logger_utils.create_progress_logger(
                self, 1, "Reciprocal-space reconstruction"
            )
            result = run_job(
                path,
                progress=update,
            )
            try:
                self._register_external_result(result)
            except Exception as error:
                self._report_failure(
                    "Reconstruction completed but result registration failed",
                    error,
                )
            self._show_output(
                json.dumps(result, indent=2, sort_keys=True)
            )
        finally:
            if progress is not None:
                try:
                    progress.finish()
                except Exception:
                    logger.warning(
                        "Cannot close reconstruction progress reporting.",
                        exc_info=True,
                    )

    @qt.Slot()
    def run_local(self):
        """Freeze live state and execute the resulting job locally."""
        try:
            self._prepare()
            self._run_path(self.job_path.text())
        except _ReconstructionCancelled as error:
            self._report_message("Reconstruction cancelled", str(error))
        except Exception as error:
            self._report_failure("Cannot run reconstruction", error)

    @qt.Slot()
    def create_cluster_scripts(self):
        """Prepare the current job and create scheduler batch scripts."""
        try:
            if self.orgui.fscan is None:
                job_path = self.job_path.text().strip()
                if not job_path:
                    raise ValueError("Select a prepared job JSON path")
                job = read_job(job_path)
            else:
                job = self._prepare()
                job_path = self.job_path.text().strip()
            result = generate_cluster_scripts(job_path, job)
            self._show_output(
                json.dumps(result, indent=2, sort_keys=True)
            )
        except Exception as error:
            self._report_failure("Cannot create cluster scripts", error)

    @qt.Slot()
    def open_job(self):
        """Open and display an existing prepared job."""
        try:
            if not self._browse(self.job_path, kind="job", save=False):
                return
            job_path = self.job_path.text().strip()
            job = read_job(job_path)
            self._show_output(
                json.dumps(job_status(job_path), indent=2)
            )
            self.output_path.setText(job.output_path)
            self.scratch_path.setText(job.scratch_path)
            self._set_accuracy(job.accuracy, job.advanced_depth)
            angle_index = self.angle_fallback.findData(job.angle_fallback)
            if angle_index < 0:
                raise ValueError(
                    "The prepared job has an unknown angle fallback: "
                    f"{job.angle_fallback}"
                )
            self.angle_fallback.setCurrentIndex(angle_index)
            self.user_note.setText(job.user_note)
            corrections = job.config_data.corrections
            self.use_pixel_mask.setChecked(corrections.use_mask)
            self.use_solid_angle.setChecked(corrections.use_solid_angle)
            self.use_polarization.setChecked(corrections.use_polarization)
            self.normalize_exposure.setChecked(
                corrections.normalize_exposure
            )
            self.monitor_corrections.setText(
                ", ".join(corrections.monitor_corrections)
            )
            self._on_monitor_corrections_changed()
            self.checkpoint_count.setValue(job.checkpoint_count)
            self._set_optional_value(
                self.thread_override, job.thread_override
            )
            self._set_optional_value(
                self.memory_override,
                None
                if job.memory_override_bytes is None
                else job.memory_override_bytes / 1024**2,
            )
            self._set_optional_value(self.frame_batch, job.frame_batch)
            tile_shape = job.tile_shape or (None, None)
            self._set_optional_value(self.tile_rows, tile_shape[0])
            self._set_optional_value(self.tile_columns, tile_shape[1])
            self._set_combo_value(
                self.work_block,
                "medium" if job.work_block_pixels is None else job.work_block_pixels,
            )
            self._set_optional_value(
                self.threads_per_image, job.threads_per_image
            )
            self._set_optional_value(
                self.accumulation_memory,
                None
                if job.accumulation_budget_bytes is None
                else job.accumulation_budget_bytes / 1024**2,
            )
            chunk_shapes = {
                tuple(values["chunk_shape"]) for values in job.grids
            }
            if len(chunk_shapes) != 1:
                raise ValueError(
                    "The prepared job contains different per-grid HDF5 chunk "
                    "shapes and cannot be represented by the global setting."
                )
            self.orgui.reconstruction_chunk_shape = chunk_shapes.pop()
            active_compression = compression_filter_name(
                self.orgui.database.compression
            )
            self.orgui.reconstruction_compression_override = (
                None
                if job.compression == active_compression
                else job.compression
            )
            self._refresh_hdf5_summary()
            self._set_cluster_settings(job.cluster_settings)
            self.grid_table.setRowCount(0)
            for values in job.grids:
                self._append_grid(ReconstructionGrid(**values))
            self._show_execution_settings(job)
        except Exception as error:
            self._report_failure("Cannot open reconstruction job", error)

    @qt.Slot()
    def resume(self):
        """Verify and resume the selected prepared job."""
        try:
            self._run_path(self.job_path.text())
        except _ReconstructionCancelled as error:
            self._report_message("Reconstruction cancelled", str(error))
        except Exception as error:
            self._report_failure("Cannot resume reconstruction", error)

    def _register_external_result(self, result):
        self.orgui.database.register_external_result(
            result["output_path"],
            result["output_sha256"],
            result["grids"],
            result["status"],
            result["job_sha256"],
        )

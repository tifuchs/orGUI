"""Shared HDF5 output settings dialog."""

from math import prod

import numpy as np
from silx.gui import qt

from .database import FILTERS


def compression_filter_name(value):
    """Return the registered name for an HDF5 compression filter."""
    for name, available in FILTERS.items():
        try:
            if value is available or value == available:
                return name
        except Exception:
            continue
    raise ValueError("The active HDF5 compression filter is not registered")


class HDF5SettingsDialog(qt.QDialog):
    """Edit shared reciprocal-space HDF5 storage settings."""

    def __init__(
        self,
        active_compression,
        chunk_shape=(64, 64, 64),
        compression_override=None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("HDF5 settings")
        layout = qt.QVBoxLayout(self)

        chunk_group = qt.QGroupBox("Spatial chunks")
        chunk_form = qt.QFormLayout(chunk_group)
        chunk_row = qt.QWidget()
        chunk_layout = qt.QHBoxLayout(chunk_row)
        chunk_layout.setContentsMargins(0, 0, 0, 0)
        chunk_tooltip = (
            "HDF5 spatial chunk dimensions in reciprocal-space voxels, shared "
            "by every output grid."
        )
        self.chunk_editors = []
        for axis, value in zip("XYZ", chunk_shape):
            label = qt.QLabel(axis)
            label.setToolTip(chunk_tooltip)
            editor = qt.QSpinBox()
            editor.setRange(1, 1_000_000)
            editor.setValue(int(value))
            editor.setToolTip(
                f"HDF5 chunk length along reciprocal-space axis {axis}."
            )
            chunk_layout.addWidget(label)
            chunk_layout.addWidget(editor)
            self.chunk_editors.append(editor)
        chunk_layout.addStretch(1)
        chunk_label = qt.QLabel("Chunk shape (voxels):")
        chunk_label.setToolTip(chunk_tooltip)
        chunk_row.setToolTip(chunk_tooltip)
        chunk_form.addRow(chunk_label, chunk_row)
        layout.addWidget(chunk_group)

        compression_group = qt.QGroupBox("Compression")
        compression_form = qt.QFormLayout(compression_group)
        self.database_compression = qt.QComboBox()
        self.database_compression.addItems(FILTERS)
        self.database_compression.setCurrentText(
            compression_filter_name(active_compression)
        )
        self.database_compression.setToolTip(
            "Compression used by normal orGUI database writes."
        )
        database_label = qt.QLabel("Active database filter:")
        database_label.setToolTip(self.database_compression.toolTip())
        compression_form.addRow(database_label, self.database_compression)
        self.override_compression = qt.QCheckBox(
            "Override for reconstruction output"
        )
        self.override_compression.setToolTip(
            "When disabled, reconstruction output uses the active orGUI "
            "database compression filter."
        )
        compression_form.addRow(self.override_compression)
        self.output_compression = qt.QComboBox()
        self.output_compression.addItems(FILTERS)
        selected = (
            compression_override
            if compression_override is not None
            else compression_filter_name(active_compression)
        )
        self.output_compression.setCurrentText(selected)
        self.output_compression.setEnabled(compression_override is not None)
        self.override_compression.setChecked(
            compression_override is not None
        )
        self.override_compression.toggled.connect(
            self.output_compression.setEnabled
        )
        self.output_compression.setToolTip(
            "Compression filter stored in the prepared job and used for its "
            "final HDF5 datasets."
        )
        output_label = qt.QLabel("Reconstruction filter:")
        output_label.setToolTip(self.output_compression.toolTip())
        compression_form.addRow(output_label, self.output_compression)
        layout.addWidget(compression_group)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    @property
    def chunk_shape(self):
        """Return the selected three-dimensional chunk shape in voxels."""
        return tuple(editor.value() for editor in self.chunk_editors)

    @property
    def compression_override(self):
        """Return the selected override name, or ``None`` for the active filter."""
        return (
            self.output_compression.currentText()
            if self.override_compression.isChecked()
            else None
        )

    @property
    def database_compression_name(self):
        """Return the selected active database compression-filter name."""
        return self.database_compression.currentText()

    def accept(self):
        """Validate settings and close the dialog."""
        chunk_bytes = prod(self.chunk_shape) * np.dtype(np.float64).itemsize
        if chunk_bytes >= 2**32:
            qt.QMessageBox.warning(
                self,
                "Invalid HDF5 chunk shape",
                "One float64 dataset chunk must be smaller than 4 GiB. "
                f"The selected shape requires {chunk_bytes / 1024**3:.2f} GiB.",
            )
            return
        super().accept()

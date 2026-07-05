"""GUI dialog for detector mask settings."""

from __future__ import annotations

import logging

import numpy as np
from silx.gui import qt

from .mask_config import PixelRepairSettings

logger = logging.getLogger(__name__)


class MaskConfigDialog(qt.QDialog):
    """Edit detector mask and pixel-repair settings.

    .. note::
       GUI-only. This dialog opens file selectors and must not be called from
       CLI or batch code.
    """

    def __init__(self, mask_manager, parent=None):
        super().__init__(parent)
        self.mask_manager = mask_manager
        self.setWindowTitle("Mask")

        layout = qt.QVBoxLayout(self)

        mask_group = qt.QGroupBox("Invalid/dead-pixel mask")
        mask_layout = qt.QGridLayout(mask_group)
        self.maskPathEdit = qt.QLineEdit()
        self.maskPathEdit.setReadOnly(True)
        load_button = qt.QPushButton("Load mask")
        clear_button = qt.QPushButton("Clear")
        mask_layout.addWidget(qt.QLabel("Mask"), 0, 0)
        mask_layout.addWidget(self.maskPathEdit, 0, 1)
        mask_layout.addWidget(load_button, 0, 2)
        mask_layout.addWidget(clear_button, 0, 3)
        layout.addWidget(mask_group)

        repair_group = qt.QGroupBox("Pixel repair")
        repair_layout = qt.QFormLayout(repair_group)
        self.repairEnabledBox = qt.QCheckBox("Enable repair")
        self.maxComponentSpin = _spin(1, 100, 4)
        self.maxSpanSpin = _spin(1, 100, 3)
        self.radiusSpin = _spin(1, 100, 2)
        self.minNeighborsSpin = _spin(1, 100, 6)
        self.usePyfaiGapsBox = qt.QCheckBox()
        self.usePyfaiGapsBox.setChecked(True)
        self.gapSizeSpin = _spin(0, 100, 6)
        repair_layout.addRow(self.repairEnabledBox)
        repair_layout.addRow("Max component pixels", self.maxComponentSpin)
        repair_layout.addRow("Max span (px)", self.maxSpanSpin)
        repair_layout.addRow("Radius (px)", self.radiusSpin)
        repair_layout.addRow("Min valid neighbors", self.minNeighborsSpin)
        repair_layout.addRow("Use pyFAI gaps", self.usePyfaiGapsBox)
        repair_layout.addRow("Fallback gap size (px)", self.gapSizeSpin)
        layout.addWidget(repair_group)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Apply
        )
        buttons.accepted.connect(self.accept)
        buttons.button(qt.QDialogButtonBox.Apply).clicked.connect(self.apply)
        layout.addWidget(buttons)

        load_button.clicked.connect(self._on_load_mask)
        clear_button.clicked.connect(self._on_clear_mask)
        self.refresh()

    def refresh(self):
        """Refresh widgets from the manager state."""

        settings = self.mask_manager.settings
        self.maskPathEdit.setText(settings.mask or "")
        repair = settings.pixel_repair or PixelRepairSettings()
        self.repairEnabledBox.setChecked(repair.enabled)
        self.maxComponentSpin.setValue(repair.max_component_pixels)
        self.maxSpanSpin.setValue(repair.max_span)
        self.radiusSpin.setValue(repair.radius)
        self.minNeighborsSpin.setValue(repair.min_valid_neighbors)
        self.usePyfaiGapsBox.setChecked(repair.use_pyfai_gaps)
        self.gapSizeSpin.setValue(repair.gap_size_px)

    def apply(self):
        """Apply current repair settings to the mask manager."""

        self.mask_manager.set_pixel_repair_settings(
            PixelRepairSettings(
                enabled=self.repairEnabledBox.isChecked(),
                max_component_pixels=self.maxComponentSpin.value(),
                max_span=self.maxSpanSpin.value(),
                radius=self.radiusSpin.value(),
                min_valid_neighbors=self.minNeighborsSpin.value(),
                use_pyfai_gaps=self.usePyfaiGapsBox.isChecked(),
                gap_size_px=self.gapSizeSpin.value(),
            )
        )

    def accept(self):
        """Apply settings and close the dialog."""

        self.apply()
        super().accept()

    # GUI-only: user-triggered mask file selector.
    def _on_load_mask(self):
        filename, _ = qt.QFileDialog.getOpenFileName(
            self,
            "Load detector mask",
            "",
            "Mask arrays (*.npy *.txt *.dat);;All files (*)",
        )
        if not filename:
            return
        try:
            mask = self.mask_manager.load_mask(filename)
        except Exception as exc:
            qt.QMessageBox.warning(self, "Cannot load mask", str(exc))
            return
        self.maskPathEdit.setText(filename)
        self._sync_plot_mask(mask)

    def _on_clear_mask(self):
        self.mask_manager.set_mask(None)
        self.mask_manager.settings.mask = None
        self.maskPathEdit.clear()
        self._sync_plot_mask(None)

    def _sync_plot_mask(self, mask):
        parent = self.parent()
        if parent is None or not hasattr(parent, "centralPlot"):
            return
        if mask is None:
            parent.centralPlot.setSelectionMask(None)
        else:
            parent.centralPlot.setSelectionMask(
                np.ascontiguousarray(mask, dtype=np.uint8)
            )


def _spin(minimum, maximum, value):
    spin = qt.QSpinBox()
    spin.setRange(minimum, maximum)
    spin.setValue(value)
    return spin

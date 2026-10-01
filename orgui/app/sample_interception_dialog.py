"""Optional sample-shape editor; scientific calculations remain in datautils."""

import copy
import logging
from pathlib import Path

import numpy as np
from silx.gui import qt
from silx.gui.plot import Plot1D, Plot2D

from .. import logger_utils
from ..backend.scans import sample_azimuth
from .sample_interception_config import (
    embed_profile,
    profile_from_settings,
    shape_frame_factors,
    shape_from_settings,
)

logger = logging.getLogger(__name__)


class SampleInterceptionDialog(qt.QDialog):
    """Edit shape geometry [mm/degrees] and an independent horizontal profile."""

    def __init__(self, footprint, settings=None):
        super().__init__(footprint)
        from .peak1Dintegr import IntegrationCorrectionsDialog

        self.footprint = footprint
        self.original = copy.deepcopy(settings or {})
        self.setWindowTitle("2D sample interception")
        self.resize(950, 750)
        layout = qt.QVBoxLayout(self)
        self.enabled = qt.QCheckBox(
            "Use exact 2D shape overlap (requires both beam profiles)"
        )
        self.enabled.setChecked(self.original.get("enabled", False))
        layout.addWidget(self.enabled)
        self.tabs = qt.QTabWidget()
        layout.addWidget(self.tabs)
        geometry = qt.QWidget()
        form = qt.QFormLayout(geometry)
        self.kind = qt.QComboBox()
        self.kind.addItems(["rectangle", "circle", "polygon"])
        self.kind.setCurrentText(
            self.original.get("shape", {}).get("kind", "rectangle")
        )
        form.addRow("Sample shape", self.kind)
        self.values = {}
        for key, label, unit, default in (
            ("length", "Rectangle length / circle diameter", " mm", 10),
            ("width", "Rectangle width", " mm", 10),
            ("offset_x", "Origin offset from axis, sample x", " mm", 0),
            ("offset_y", "Origin offset from axis, sample y", " mm", 0),
            ("orientation_deg", "Shape orientation at reference azimuth", " deg", 0),
            (
                "reference_incidence_deg",
                "Reference incidence for profile alignment",
                " deg",
                0.36,
            ),
            ("reference_azimuth_deg", "Reference readback azimuth", " deg", 0),
            ("fixed_azimuth_deg", "Fixed azimuth, if explicitly selected", " deg", 0),
            ("overspill_threshold", "Overspill warning threshold", " %", 1),
        ):
            spin = qt.QDoubleSpinBox()
            spin.setRange(-1e6, 1e6)
            spin.setDecimals(6)
            spin.setSuffix(unit)
            spin.setValue(self.original.get(key, default))
            self.values[key] = spin
            form.addRow(label, spin)
        dimensions = self.original.get("shape", {}).get("dimensions_m", [0.01, 0.01])
        if dimensions:
            self.values["length"].setValue(dimensions[0] * 1e3)
        if len(dimensions) > 1:
            self.values["width"].setValue(dimensions[1] * 1e3)
        offset = self.original.get("offset_m", [0, 0])
        for key, value in zip(("offset_x", "offset_y"), offset):
            self.values[key].setValue(value * 1e3)
        self.values["overspill_threshold"].setValue(
            self.original.get("overspill_threshold", 0.01) * 100
        )
        self.source = qt.QLineEdit(self.original.get("azimuth_source", ""))
        self.source.setPlaceholderText("Exact motor/counter name, or fixed")
        form.addRow("Actual sample azimuth source", self.source)
        self.unit = qt.QComboBox()
        self.unit.addItems(["deg", "rad"])
        self.unit.setCurrentText(self.original.get("azimuth_unit", "deg"))
        form.addRow("Readback unit", self.unit)
        self.sign = qt.QComboBox()
        self.sign.addItems(["+1", "-1"])
        self.sign.setCurrentText(f"{int(self.original.get('azimuth_sign', 1)):+d}")
        form.addRow("Readback to geometric rotation sign", self.sign)
        self.confirm = qt.QCheckBox(
            "This source represents rotation about the surface normal"
        )
        self.confirm.setChecked(self.original.get("normal_rotation_confirmed", False))
        form.addRow(self.confirm)
        self.vertices = qt.QPlainTextEdit()
        self.vertices.setPlaceholderText(
            "Polygon vertices: x y in mm, one pair per line; # comments allowed"
        )
        self.vertices.setMaximumHeight(90)
        self.vertices.setPlainText(
            "\n".join(
                f"{x * 1e3:g} {y * 1e3:g}"
                for x, y in self.original.get("shape", {}).get("vertices_m", [])
            )
        )
        form.addRow("Polygon (supplied coordinate origin is retained)", self.vertices)
        import_button = qt.QPushButton("Import polygon text / CSV…")
        import_button.clicked.connect(self._importPolygon)
        form.addRow(import_button)
        note = qt.QLabel(
            "Profile centre/offset aligns the shape origin at the reference "
            "incidence/azimuth. "
            "The offset rotates with the sample; the beam alignment stays fixed. "
            "Legacy shape area uses peak flux density and geometric sample "
            "area as its reference."
        )
        note.setWordWrap(True)
        form.addRow(note)
        self.tabs.addTab(geometry, "Shape and alignment")
        self.horizontal = IntegrationCorrectionsDialog(self, profile_only=True)
        self.horizontal.setWindowFlags(qt.Qt.Widget)
        if self.original.get("horizontal"):
            self.horizontal.setSettings(self.original["horizontal"])
        self.tabs.addTab(self.horizontal, "Horizontal beam profile")
        preview = qt.QWidget()
        preview_layout = qt.QVBoxLayout(preview)
        frame_controls = qt.QHBoxLayout()
        self.scan_preview = qt.QCheckBox("Preview loaded source frames")
        self.frame_index = qt.QSpinBox()
        self.frame_index.setPrefix("Frame ")
        main_getter = getattr(footprint.parent(), "_mainWindow", lambda: None)
        self.main_window = main_getter()
        scan = getattr(self.main_window, "fscan", None)
        self.scan_preview.setEnabled(scan is not None)
        self.frame_index.setRange(0, max(0, len(scan) - 1) if scan is not None else 0)
        frame_controls.addWidget(self.scan_preview)
        frame_controls.addWidget(self.frame_index)
        preview_layout.addLayout(frame_controls)
        preview_controls = qt.QHBoxLayout()
        self.preview_alpha = qt.QDoubleSpinBox()
        self.preview_alpha.setRange(0.000001, 90)
        self.preview_alpha.setDecimals(6)
        self.preview_alpha.setValue(0.36)
        self.preview_alpha.setSuffix(" deg incidence")
        self.preview_azimuth = qt.QDoubleSpinBox()
        self.preview_azimuth.setRange(-1e6, 1e6)
        self.preview_azimuth.setSuffix(" deg readback azimuth")
        preview_controls.addWidget(self.preview_alpha)
        preview_controls.addWidget(self.preview_azimuth)
        refresh = qt.QPushButton("Calculate preview")
        refresh.clicked.connect(self._preview)
        preview_controls.addWidget(refresh)
        preview_layout.addLayout(preview_controls)
        self.overlap_plot = Plot1D()
        self.overlap_plot.setGraphXLabel("Readback azimuth (deg)")
        self.overlap_plot.setGraphYLabel("f_hit (fraction)")
        self.overlap_plot.setGraphYLabel("H (dimensionless)", axis="right")
        preview_layout.addWidget(self.overlap_plot)
        self.footprint_plot = Plot2D()
        self.footprint_plot.setKeepDataAspectRatio(True)
        self.footprint_plot.setGraphXLabel("Along beam x (mm)")
        self.footprint_plot.setGraphYLabel("Transverse y (mm)")
        preview_layout.addWidget(self.footprint_plot)
        self.tabs.addTab(preview, "Overlap diagnostics")
        self.status = qt.QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Ok | qt.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self._acceptValidated)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def settings(self):
        """Return SI shape geometry and dialog-unit/embedded horizontal profile."""
        result = {
            k: w.value()
            for k, w in self.values.items()
            if k not in {"length", "width", "offset_x", "offset_y"}
        }
        result.update(
            {
                "version": 1,
                "enabled": self.enabled.isChecked(),
                "azimuth_source": self.source.text().strip(),
                "azimuth_unit": self.unit.currentText(),
                "azimuth_sign": int(self.sign.currentText()),
                "normal_rotation_confirmed": self.confirm.isChecked(),
                "offset_m": [
                    self.values[k].value() * 1e-3 for k in ("offset_x", "offset_y")
                ],
                "overspill_threshold": self.values["overspill_threshold"].value() / 100,
                "rtol": self.original.get("rtol", 1e-9),
            }
        )
        kind = self.kind.currentText()
        result["shape"] = {"kind": kind}
        if kind == "polygon":
            rows = [
                line.split("#")[0].replace(",", " ").split()
                for line in self.vertices.toPlainText().splitlines()
            ]
            result["shape"]["vertices_m"] = [
                [float(v) * 1e-3 for v in row] for row in rows if row
            ]
        else:
            keys = ["length", "width"] if kind == "rectangle" else ["length"]
            result["shape"]["dimensions_m"] = [
                self.values[k].value() * 1e-3 for k in keys
            ]
        horizontal = self.horizontal.settings()
        horizontal.pop("sample_interception", None)
        try:
            result["horizontal"] = embed_profile(
                horizontal, self.horizontal.beamProfile()
            )
        except ValueError:
            if result["enabled"]:
                raise
            result["horizontal"] = horizontal
        return result

    def _acceptValidated(self):
        try:
            settings = self.settings()
            if settings["enabled"]:
                if not settings["azimuth_source"]:
                    raise ValueError(
                        "provide an actual azimuth source or explicitly select fixed"
                    )
                shape_frame_factors(
                    settings,
                    self.footprint.beamProfile(),
                    np.deg2rad(self.preview_alpha.value()),
                    np.deg2rad(self.preview_azimuth.value()),
                )
        except (ValueError, NotImplementedError) as error:
            self.status.setText(str(error))
            return
        self.accept()

    def _importPolygon(self):
        # GUI-only: explicitly invoked file dialog.
        path, _ = qt.QFileDialog.getOpenFileName(
            self, "Import polygon", "", "Text/CSV (*.txt *.csv *.dat)"
        )
        if path:
            try:
                self.vertices.setPlainText(Path(path).read_text())
                self.kind.setCurrentText("polygon")
            except OSError as error:
                self.status.setText(str(error))

    # GUI-only: explicitly requested preview, with cancellable progress.
    def _preview(self):
        if getattr(self, "_preview_running", False):
            return
        self._preview_running = True
        progress = None
        try:
            settings = self.settings()
            shape = shape_from_settings(settings)
            vertical = self.footprint.beamProfile()
            angles = np.linspace(
                self.preview_azimuth.value() - 90, self.preview_azimuth.value() + 90, 37
            )
            alpha = np.deg2rad(self.preview_alpha.value())
            selected_alpha = alpha
            selected_azimuth = np.deg2rad(self.preview_azimuth.value())
            azimuth = np.deg2rad(angles)
            xlabel = "Readback azimuth (deg)"
            if self.scan_preview.isChecked():
                scan = getattr(self.main_window, "fscan", None)
                if scan is None:
                    raise ValueError("load a scan before previewing source frames")
                alpha = np.broadcast_to(self.main_window.getMuOm()[0], (len(scan),))
                azimuth = sample_azimuth(scan, settings, len(scan))
                index = min(self.frame_index.value(), len(scan) - 1)
                selected_alpha, selected_azimuth = alpha[index], azimuth[index]
                angles = np.arange(len(scan))
                xlabel = "Source frame index"
            count = np.size(azimuth)
            if logger_utils.get_logging_context() in {"gui", "cli"}:
                progress = logger_utils.create_progress_logger(
                    self, count, "Sample interception preview"
                )
            else:
                progress = logger_utils.LogProgress(
                    __name__, count, "Sample interception preview"
                )
            self.tabs.setEnabled(False)

            def update(completed, total):
                """Report bounded work and honour the preview cancellation button."""
                progress.update(completed)
                return not progress.wasCanceled()

            h, fraction, _, valid, _ = shape_frame_factors(
                settings, vertical, alpha, azimuth, progress=update
            )
            self.overlap_plot.clear()
            self.overlap_plot.setGraphXLabel(xlabel)
            self.overlap_plot.addCurve(angles, fraction, legend="f_hit")
            self.overlap_plot.addCurve(angles, h, legend="H", yaxis="right")
            angle = settings["azimuth_sign"] * np.deg2rad(
                np.rad2deg(selected_azimuth) - settings["reference_azimuth_deg"]
            )
            outline = shape.outline(
                angle, settings["offset_m"], np.deg2rad(settings["orientation_deg"])
            )
            bounds = np.column_stack((outline.min(axis=0), outline.max(axis=0)))
            span = np.maximum(bounds[:, 1] - bounds[:, 0], 1e-6)
            lower, upper = bounds[:, 0] - 0.2 * span, bounds[:, 1] + 0.2 * span
            x, y = (
                np.linspace(lower[0], upper[0], 161),
                np.linspace(lower[1], upper[1], 161),
            )
            horizontal = profile_from_settings(settings["horizontal"])
            zv = -settings["offset_m"][0] * np.sin(
                np.deg2rad(settings["reference_incidence_deg"])
            )
            yh = -settings["offset_m"][1]
            density = horizontal.density_at(y[:, None] + yh) * vertical.density_at(
                x[None, :] * np.sin(selected_alpha) + zv
            )
            self.footprint_plot.clear()
            self.footprint_plot.addImage(
                density,
                origin=tuple(lower * 1e3),
                scale=tuple((upper - lower) / 160 * 1e3),
                legend="Beam density (1/m²)",
            )
            closed = np.vstack((outline, outline[0])) * 1e3
            self.footprint_plot.addCurve(
                closed[:, 0], closed[:, 1], legend="Sample", color="red"
            )
            self.footprint_plot.addMarker(0, 0, legend="Rotation axis", symbol="+")
            self.status.setText(
                f"{np.count_nonzero(valid)} valid preview angles. "
                f"Outline at incidence {np.rad2deg(selected_alpha):g} deg, "
                f"readback {np.rad2deg(selected_azimuth):g} deg. "
                "Representative frame angles are used; measured profiles "
                "have zero density outside saved support."
            )
        except (ValueError, NotImplementedError, InterruptedError) as error:
            self.status.setText(str(error))
        finally:
            if progress is not None:
                progress.finish()
            self.tabs.setEnabled(True)
            self._preview_running = False

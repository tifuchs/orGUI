"""Sample and overlap pages embedded in the shared correction editor."""

import copy
import logging
from pathlib import Path

import numpy as np
from silx.gui import qt
from silx.gui.plot import Plot1D, Plot2D

from .. import logger_utils
from .sample_interception_config import (
    embed_profile,
    profile_from_settings,
    shape_frame_factors,
    shape_from_settings,
    sample_angle_inputs,
    parallel_omega_reading,
    orientation_for_parallel_reading,
)

logger = logging.getLogger(__name__)


class SampleInterceptionWidget(qt.QWidget):
    """Edit shape geometry [mm/degrees] using the owner's size controls."""

    def __init__(self, footprint, settings=None):
        super().__init__(footprint)
        from .peak1Dintegr import IntegrationCorrectionsDialog

        self.footprint = footprint
        self.original = copy.deepcopy(settings or {})
        self.configured = bool(self.original)
        layout = qt.QVBoxLayout(self)
        self.enabled = qt.QCheckBox(
            "Use exact 2D shape overlap (requires both beam profiles)"
        )
        self.enabled.setChecked(self.original.get("enabled", False))
        self.tabs = qt.QTabWidget()
        layout.addWidget(self.tabs)
        geometry = qt.QWidget()
        self.geometry = geometry
        form = qt.QFormLayout(geometry)
        self.form = form
        form.addRow(self.enabled)
        self.kind = qt.QComboBox()
        self.kind.addItems(["rectangle", "circle", "polygon"])
        self.kind.setCurrentText(
            self.original.get("shape", {}).get("kind", "rectangle")
        )
        form.addRow("Sample shape", self.kind)
        # These are the very same controls exposed by the legacy API.
        self.values = {"length": footprint.L, "width": footprint.W}
        self.length_label = qt.QLabel("Sample length L")
        self.width_label = qt.QLabel("Sample width W")
        form.addRow(self.length_label, footprint.L)
        form.addRow(self.width_label, footprint.W)
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
        self.import_button = import_button = qt.QPushButton(
            "Import polygon text / CSV…"
        )
        import_button.clicked.connect(self._importPolygon)
        form.addRow(import_button)
        self.angle_summary = qt.QLabel(
            "Incidence: mu (configuration / scan axis). "
            "Azimuth: omega = -th (scan)."
        )
        self.angle_summary.setWordWrap(True)
        form.addRow("Angles from scan / config", self.angle_summary)
        self.override_button = qt.QPushButton("Manual angle overrides…")
        self.override_button.setCheckable(True)
        self.override_panel = qt.QWidget()
        override_form = qt.QFormLayout(self.override_panel)
        self.angle_modes = {}
        self.angle_sources = {}
        self.angle_units = {}
        self.fixed_angles = {}
        for name in ("incidence", "azimuth"):
            mode = qt.QComboBox()
            mode.addItem(
                "mu (configuration / scan axis)" if name == "incidence" else
                "omega = -th (scan)", "auto"
            )
            mode.addItem("Counter", "counter")
            mode.addItem("Fixed", "fixed")
            source = qt.QLineEdit()
            source.setPlaceholderText("Exact loaded counter name")
            unit = qt.QComboBox()
            unit.addItems(["deg", "rad"])
            fixed = qt.QDoubleSpinBox()
            fixed.setRange(-1e6, 1e6)
            fixed.setDecimals(6)
            fixed.setSuffix(" deg")
            self.angle_modes[name] = mode
            self.angle_sources[name] = source
            self.angle_units[name] = unit
            self.fixed_angles[name] = fixed
            override_form.addRow(name.capitalize(), mode)
            override_form.addRow("Counter", source)
            override_form.addRow("Counter unit", unit)
            override_form.addRow("Fixed angle", fixed)
            mode.currentIndexChanged.connect(self._angleModeChanged)
            source.textChanged.connect(self._changed)
            unit.currentIndexChanged.connect(self._changed)
            fixed.valueChanged.connect(self._changed)
        self.sign = qt.QComboBox()
        self.sign.addItems(["+1", "-1"])
        self.sign.setToolTip(
            "edge angle = sign × (azimuth reading − parallel reading). "
            "+1 means increasing readings turn right-handed about the surface "
            "normal; -1 reverses that direction. Acquisition omega is -th. "
            "Alignment readings are always in degrees, even for radian counters."
        )
        override_form.addRow("Azimuth reading to geometric rotation sign", self.sign)
        self.source = self.angle_sources["azimuth"]
        self.unit = self.angle_units["azimuth"]
        form.addRow(self.override_button)
        form.addRow(self.override_panel)
        self.override_button.toggled.connect(self._updateGeometryControls)
        self.parallel_reading = qt.QDoubleSpinBox()
        self.parallel_reading.setRange(-1e12, 1e12)
        self.parallel_reading.setDecimals(12)
        self.parallel_reading.setSuffix(" deg")
        self.parallel_label = qt.QLabel()
        self.parallel_label.setWordWrap(True)
        form.addRow(self.parallel_label, self.parallel_reading)
        self.alignment_note = qt.QLabel()
        self.alignment_note.setWordWrap(True)
        form.addRow(self.alignment_note)
        self.off_centre_button = qt.QToolButton()
        self.off_centre_button.setText("Off-centre sample")
        self.off_centre_button.setToolButtonStyle(qt.Qt.ToolButtonTextBesideIcon)
        self.off_centre_button.setCheckable(True)
        self.off_centre_button.toggled.connect(self._updateGeometryControls)
        form.addRow(self.off_centre_button)
        self.off_centre_panel = qt.QWidget()
        self.placement_form = placement_form = qt.QFormLayout(self.off_centre_panel)
        form.addRow(self.off_centre_panel)
        self.reference_override = qt.QCheckBox("Set profile-alignment incidence")
        self.reference_override.setToolTip(
            "Otherwise use the first valid source-frame incidence. "
            "Only affects alignment for a nonzero sample x offset."
        )
        self.reference_override.toggled.connect(self._angleModeChanged)
        self.diagnostic_options = qt.QWidget()
        diagnostic_form = qt.QFormLayout(self.diagnostic_options)
        for key, label, unit, default in (
            ("offset_x", "Shape origin from axis, reference beam x", " mm", 0),
            ("offset_y", "Shape origin from axis, reference beam y", " mm", 0),
            ("reference_azimuth_deg", "Placement reference azimuth reading", " deg", 0),
            ("orientation_deg", "Shape orientation at placement reference", " deg", 0),
            ("reference_incidence_deg", "Reference incidence for profile alignment",
             " deg", 0.36),
            ("overspill_threshold", "Overspill warning threshold", " %", 1),
        ):
            spin = qt.QDoubleSpinBox()
            spin.setRange(-1e12, 1e12)
            spin.setDecimals(12)
            spin.setSuffix(unit)
            spin.setValue(self.original.get(key, default))
            self.values[key] = spin
            if key == "reference_incidence_deg":
                placement_form.addRow(self.reference_override)
            target_form = (diagnostic_form if key == "overspill_threshold"
                           else placement_form)
            target_form.addRow(label, spin)
        self.values["orientation_deg"].setToolTip(
            "Angle of the rectangle's length direction (or polygon x axis) "
            "from the projected beam at the placement reference reading. "
            "Positive is right-handed about the surface normal."
        )
        self.values["reference_azimuth_deg"].setToolTip(
            "Azimuth reading at which the offset coordinates and profile "
            "alignment are defined. Changing this changes the origin's orbit."
        )
        placement_note = qt.QLabel(
            "The displacement rotates with the sample from the placement reference; "
            "beam alignment stays fixed. Collapsing this section preserves all values."
        )
        placement_note.setWordWrap(True)
        placement_form.addRow(placement_note)
        self.sample_page = qt.QScrollArea()
        self.sample_page.setWidgetResizable(True)
        self.sample_page.setFrameShape(qt.QFrame.NoFrame)
        self.sample_page.setWidget(geometry)
        self.tabs.addTab(self.sample_page, "Sample")
        self.horizontal = IntegrationCorrectionsDialog(self, profile_only=True)
        self.horizontal.setWindowFlags(qt.Qt.Widget)
        self.tabs.addTab(self.horizontal, "Horizontal beam")
        self.preview = preview = qt.QWidget()
        self.preview_layout = preview_layout = qt.QVBoxLayout(preview)
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
        preview_layout.addWidget(self.diagnostic_options)
        self.overlap_plot = None
        self.footprint_plot = None
        self.tabs.addTab(preview, "Diagnostics")
        self.status = qt.QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        self._horizontal_defaults = self.horizontal.settings()
        self.setSettings(settings or {})
        self.enabled.toggled.connect(self._geometryChanged)
        self.kind.currentIndexChanged.connect(self._geometryChanged)
        for key, spin in self.values.items():
            if key in {"length", "width"}:
                spin.valueChanged.connect(self._dimensionsChanged)
            else:
                if key in {"reference_azimuth_deg", "orientation_deg",
                           "reference_incidence_deg", "offset_x", "offset_y"}:
                    spin.valueChanged.connect(self._placementChanged)
                else:
                    spin.valueChanged.connect(self._changed)
        self.vertices.textChanged.connect(self._changed)
        self.sign.currentIndexChanged.connect(self._signChanged)
        self.parallel_reading.valueChanged.connect(self._parallelChanged)
        self.horizontal.settingsChanged.connect(self._changed)
        self.tabs.currentChanged.connect(self._diagnosticsSelected)
        self.refreshSourceFrames()

    def _diagnosticsSelected(self, index):
        if self.tabs.widget(index) is self.preview:
            self._ensurePreviewPlots()

    def _ensurePreviewPlots(self):
        # GUI-only: plots are needed only for the diagnostics tab/preview.
        if self.overlap_plot is not None:
            return
        self.overlap_plot = Plot1D()
        self.overlap_plot.setGraphXLabel("Readback azimuth (deg)")
        self.overlap_plot.setGraphYLabel("f_hit (fraction)")
        self.overlap_plot.setGraphYLabel("H (dimensionless)", axis="right")
        self.preview_layout.addWidget(self.overlap_plot)
        self.footprint_plot = Plot2D()
        self.footprint_plot.setKeepDataAspectRatio(True)
        self.footprint_plot.setGraphXLabel("Along beam x (mm)")
        self.footprint_plot.setGraphYLabel("Transverse y (mm)")
        self.preview_layout.addWidget(self.footprint_plot)

    def setSettings(self, settings):
        """Restore shape state; enabled dimensions override legacy size copies.

        Geometry lengths are metres in saved settings and millimetres in the
        shared controls. An empty model retains the legacy rectangular sizes.
        """
        self._restoring = True
        try:
            self.original = copy.deepcopy(settings)
            # Keep full saved precision. Display conversion must not rewrite
            # the placement reference/orientation or offset orbit on load.
            self._placement_values = {
                "reference_azimuth_deg": settings.get("reference_azimuth_deg", 0),
                "orientation_deg": settings.get("orientation_deg", 0),
                "reference_incidence_deg": settings.get(
                    "reference_incidence_deg", 0.36
                ),
                "offset_m": list(settings.get("offset_m", (0, 0))),
            }
            self.configured = bool(settings)
            self.enabled.setChecked(settings.get("enabled", False))
            shape = settings.get("shape", {})
            self.kind.setCurrentText(shape.get("kind", "rectangle"))
            dimensions = shape.get("dimensions_m", [])
            if settings.get("enabled", False) and dimensions:
                self.values["length"].setValue(dimensions[0] * 1e3)
                if len(dimensions) > 1:
                    self.values["width"].setValue(dimensions[1] * 1e3)
            defaults = {
                "offset_x": 0, "offset_y": 0, "orientation_deg": 0,
                "reference_incidence_deg": 0.36, "reference_azimuth_deg": 0,
                "overspill_threshold": 1,
            }
            for key, default in defaults.items():
                self.values[key].setValue(settings.get(key, default))
            for key, value in zip(
                ("offset_x", "offset_y"), settings.get("offset_m", [0, 0])
            ):
                self.values[key].setValue(value * 1e3)
            self.values["overspill_threshold"].setValue(
                settings.get("overspill_threshold", 0.01) * 100
            )
            for name in ("incidence", "azimuth"):
                source = settings.get(f"{name}_source", "auto")
                mode = source if source in {"auto", "fixed"} else "counter"
                self.angle_modes[name].setCurrentIndex(
                    self.angle_modes[name].findData(mode)
                )
                self.angle_sources[name].setText(
                    source if mode == "counter" else ""
                )
                self.angle_units[name].setCurrentText(
                    settings.get(f"{name}_unit", "deg")
                )
                self.fixed_angles[name].setValue(settings.get(f"fixed_{name}_deg", 0))
            self.sign.setCurrentText(f"{int(settings.get('azimuth_sign', 1)):+d}")
            self.override_button.setChecked(any(
                mode.currentData() != "auto" for mode in self.angle_modes.values()
            ) or settings.get("azimuth_sign", 1) != 1)
            self.reference_override.setChecked("reference_incidence_deg" in settings)
            self.off_centre_button.setChecked(any(self._placement_values["offset_m"]))
            self.vertices.setPlainText("\n".join(
                f"{x * 1e3:.17g} {y * 1e3:.17g}"
                for x, y in shape.get("vertices_m", [])
            ))
            self.horizontal.setSettings(
                settings.get("horizontal", self._horizontal_defaults)
            )
            self._syncAlignment()
            self._updateGeometryControls()
            self._updateAngleSummary()
        finally:
            self._restoring = False

    def _changed(self, *args):
        if getattr(self, "_restoring", False):
            return
        self.configured = True
        self._updateAngleSummary()
        self.footprint._settingsChanged()

    def _geometryChanged(self, *args):
        if not getattr(self, "_restoring", False):
            self._syncAlignment()
        self._updateGeometryControls()
        self.footprint._updatePreview()
        self._changed()

    def _angleModeChanged(self, *args):
        if not hasattr(self, "off_centre_button"):
            return
        self._updateGeometryControls()
        self._changed()

    def _syncAlignment(self):
        settings = self._alignmentSettings()
        self._alignment_reading = parallel_omega_reading(settings)
        blocker = qt.QSignalBlocker(self.parallel_reading)
        self.parallel_reading.setValue(self._alignment_reading)
        del blocker

    def _parallelChanged(self, reading):
        if getattr(self, "_restoring", False):
            return
        settings = self._alignmentSettings()
        self._placement_values["orientation_deg"] = orientation_for_parallel_reading(
            settings, reading
        )
        blocker = qt.QSignalBlocker(self.values["orientation_deg"])
        self.values["orientation_deg"].setValue(self._placement_values["orientation_deg"])
        del blocker
        self._syncAlignment()
        self._changed()

    def _placementChanged(self, *args):
        if getattr(self, "_restoring", False):
            return
        for key in ("reference_azimuth_deg", "orientation_deg",
                    "reference_incidence_deg", "offset_x", "offset_y"):
            if self.sender() is self.values[key]:
                if key.startswith("offset_"):
                    index = 0 if key == "offset_x" else 1
                    self._placement_values["offset_m"][index] = (
                        self.values[key].value() * 1e-3
                    )
                else:
                    self._placement_values[key] = self.values[key].value()
                break
        self._syncAlignment()
        self._updateGeometryControls()
        self._changed()

    def _alignmentSettings(self):
        return dict(
            self._placement_values,
            azimuth_sign=int(self.sign.currentText()),
            shape={"kind": self.kind.currentText(), "dimensions_m": [
                self.values[key].value() * 1e-3 for key in ("length", "width")
            ]},
        )

    def _dimensionsChanged(self, *args):
        if getattr(self, "_restoring", False):
            return
        self._syncAlignment()

    def _signChanged(self, *args):
        if getattr(self, "_restoring", False):
            return
        if not any(self._placement_values["offset_m"]):
            # The centred control defines an alignment reading, so reversing
            # the motor direction retains that reading. Off-centre references
            # retain their separate placement/orientation parameters instead.
            self._parallelChanged(self._alignment_reading)
        else:
            self._syncAlignment()
            self._changed()

    def _updateAngleSummary(self):
        base = getattr(self, "_angle_summary_base",
                       "Incidence: mu (configuration / scan axis). "
                       "Azimuth: omega = -th (scan).")
        overrides = []
        for name, mode in self.angle_modes.items():
            if mode.currentData() == "fixed":
                overrides.append(f"{name}: {self.fixed_angles[name].value():g} deg")
            elif mode.currentData() == "counter":
                overrides.append(
                    f"{name}: counter {self.angle_sources[name].text().strip() or '?'} "
                    f"({self.angle_units[name].currentText()})"
                )
        if overrides:
            base += " Footprint overrides: " + "; ".join(overrides) + "."
        base += f" Geometric rotation sign: {self.sign.currentText()}."
        self.angle_summary.setText(base)
        self.override_button.setText(
            "Manual angle overrides (active)…"
            if overrides or self.sign.currentText() == "-1" else
            "Manual angle overrides…"
        )

    def _updateGeometryControls(self):
        exact = self.enabled.isChecked()
        kind = self.kind.currentText() if exact else "rectangle"
        self.kind.setEnabled(exact)
        off_centre = any(self._placement_values["offset_m"])
        self.off_centre_button.setArrowType(
            qt.Qt.DownArrow if self.off_centre_button.isChecked() else qt.Qt.RightArrow
        )
        for widget in (self.kind, self.alignment_note, self.angle_summary,
                       self.override_button, self.off_centre_button):
            widget.setVisible(exact)
            label = self.form.labelForField(widget)
            if label is not None:
                label.setVisible(exact)
        self.off_centre_panel.setVisible(exact and self.off_centre_button.isChecked())
        self.parallel_reading.setVisible(exact and not off_centre and kind != "circle")
        self.parallel_label.setVisible(exact and not off_centre and kind != "circle")
        self.parallel_label.setText(
            "Omega reading when the polygon's x axis is parallel to the beam (deg)"
            if kind == "polygon" else
            "Omega reading when the sample’s long edge is parallel to the beam (deg)"
        )
        self.alignment_note.setText(
            "Off-centre placement is active. Expand Off-centre sample to edit "
            "the placement reference and shape orientation."
            if off_centre else
            "A centred circle needs no angular alignment."
            if kind == "circle" else
            "Parallel means a 0° edge angle; perpendicular means 90°. "
            "Enter the selected azimuth source's reading in degrees. "
            "For polygons, the reference direction is the supplied x axis."
        )
        self.values["orientation_deg"].setVisible(kind != "circle")
        self.placement_form.labelForField(self.values["orientation_deg"]).setVisible(
            kind != "circle"
        )
        self.values["reference_incidence_deg"].setEnabled(
            self.reference_override.isChecked()
        )
        self.override_panel.setVisible(exact and self.override_button.isChecked())
        override_form = self.override_panel.layout()
        for name, mode in self.angle_modes.items():
            for widget, visible in (
                (self.angle_sources[name], mode.currentData() == "counter"),
                (self.angle_units[name], mode.currentData() == "counter"),
                (self.fixed_angles[name], mode.currentData() == "fixed"),
            ):
                widget.setVisible(visible)
                override_form.labelForField(widget).setVisible(visible)
        self.length_label.setText("Circle diameter" if kind == "circle"
                                  else "Sample length L")
        for widget in (self.length_label, self.values["length"]):
            widget.setVisible(kind != "polygon")
        for widget in (self.width_label, self.values["width"]):
            widget.setVisible(kind == "rectangle")
        for widget in (self.vertices, self.form.labelForField(self.vertices),
                       self.import_button):
            widget.setVisible(exact and kind == "polygon")
        self.footprint._updateLegacyControlState()

    def settings(self):
        """Return SI shape geometry and dialog-unit/embedded horizontal profile."""
        result = {
            k: w.value()
            for k, w in self.values.items()
            if k not in {"length", "width", "offset_x", "offset_y",
                         "reference_incidence_deg"}
        }
        result.update(
            {
                "version": 1,
                "enabled": self.enabled.isChecked(),
                "azimuth_sign": int(self.sign.currentText()),
                # Retained for readers of older version-1 settings.
                "normal_rotation_confirmed": True,
                "overspill_threshold": self.values["overspill_threshold"].value() / 100,
                "rtol": self.original.get("rtol", 1e-9),
            }
        )
        result.update(copy.deepcopy(self._placement_values))
        result.pop("reference_incidence_deg")
        if self.reference_override.isChecked():
            result["reference_incidence_deg"] = self._placement_values[
                "reference_incidence_deg"
            ]
        for name, mode in self.angle_modes.items():
            result[f"{name}_source"] = (
                self.angle_sources[name].text().strip()
                if mode.currentData() == "counter" else
                "fixed" if mode.currentData() == "fixed" else "auto"
            )
            result[f"{name}_unit"] = self.angle_units[name].currentText()
            result[f"fixed_{name}_deg"] = self.fixed_angles[name].value()
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

    def validate(self):
        """Validate the enabled shape and both profiles without opening a dialog.

        :raises ValueError: For invalid geometry, missing sources or profiles.
        """
        settings = self.settings()
        if settings["enabled"]:
            for name in ("incidence", "azimuth"):
                if not settings[f"{name}_source"]:
                    raise ValueError(f"provide a {name} override counter name")
            shape_frame_factors(
                settings, self.footprint.beamProfile(),
                np.deg2rad(self.preview_alpha.value()),
                np.deg2rad(self.preview_azimuth.value()),
            )

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

    def refreshSourceFrames(self):
        """Refresh source-frame availability from the currently loaded scan."""
        scan = getattr(self.main_window, "fscan", None)
        self.scan_preview.setEnabled(scan is not None)
        self.frame_index.setRange(0, max(0, len(scan) - 1) if scan is not None else 0)
        axis = getattr(scan, "axisname", None)
        sources = {
            "incidence": ("mu (scan axis)" if axis == "mu" else
                          "mu (configuration)"),
            "azimuth": ("omega = -th (scan axis)" if axis == "th" else
                        "omega = -th (scan readback)" if scan is not None else
                        "omega = -th (scan)"),
        }
        for name, mode in self.angle_modes.items():
            # Display text may change on scan load; the saved source stays auto.
            blocker = qt.QSignalBlocker(mode)
            mode.setItemText(mode.findData("auto"), sources[name])
            del blocker
        self._angle_summary_base = (
            f"Incidence: {sources['incidence']}. Azimuth: {sources['azimuth']}."
        )
        if scan is not None and hasattr(self.main_window, "getMuOm"):
            try:
                alpha, omega = self.main_window.getMuOm()
                self._angle_summary_base = (
                    f"Incidence: {sources['incidence']}, "
                    f"{np.rad2deg(np.ravel(alpha)[0]):g} deg; "
                    f"azimuth: {sources['azimuth']}, "
                    f"{np.rad2deg(np.ravel(omega)[0]):g} deg. "
                    f"Scan axis: {getattr(scan, 'axisname', 'unknown')}."
                )
            except (AttributeError, ValueError, IndexError) as error:
                self._angle_summary_base = f"Scan angles unavailable: {error}"
        self._updateAngleSummary()

    # GUI-only: explicitly requested preview, with cancellable progress.
    def _preview(self):
        if getattr(self, "_preview_running", False):
            return
        self._preview_running = True
        progress = None
        try:
            self._ensurePreviewPlots()
            self.refreshSourceFrames()
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
                mu, omega = self.main_window.getMuOm()
                alpha, azimuth, settings = sample_angle_inputs(
                    scan, settings, mu,
                    count=len(scan), omega=omega,
                )
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
                np.deg2rad(settings.get(
                    "reference_incidence_deg", np.rad2deg(np.ravel(alpha)[0])
                ))
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

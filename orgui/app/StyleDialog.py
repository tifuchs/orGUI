"""Runtime Qt widget-style selection dialog."""

from silx.gui import qt
from silx.gui.utils import glutils

from .ui_preferences import OPENGL_KEY, STYLE_KEY


class StyleDialog(qt.QDialog):
    """Select and apply an installed Qt widget style."""

    def __init__(self, application=None, plots=(), settings=None, parent=None):
        """Create a dialog for the given Qt application.

        :param QApplication application: Application whose widget style is
            changed. The current Qt application is used when omitted.
        :param tuple plots: Plot widgets whose current backend is reported.
        :param QSettings settings: Optional store used for persistent settings.
        :param QWidget parent: Optional parent widget.
        :raises RuntimeError: If no Qt application is available.
        """
        super().__init__(parent)
        self.application = application or qt.QApplication.instance()
        if self.application is None:
            raise RuntimeError("A QApplication is required to select a style")
        self.plots = tuple(plots)
        self.settings = settings if settings is not None else qt.QSettings()

        self.setWindowTitle("Application style")
        layout = qt.QVBoxLayout(self)
        form = qt.QFormLayout()
        self.style_selector = qt.QComboBox()
        self.style_selector.addItems(qt.QStyleFactory.keys())
        self.style_selector.setCurrentText(self.application.style().objectName())
        self.style_selector.setToolTip(
            "Qt widget styles available in this installation."
        )
        form.addRow("Style:", self.style_selector)

        self.opengl_selector = qt.QCheckBox("OpenGL rendering")
        self._plot_opengl_states = tuple(
            self._plot_uses_opengl(plot) for plot in self.plots
        )
        self.opengl_selector.setChecked(
            bool(self._plot_opengl_states) and all(self._plot_opengl_states)
        )
        self.opengl_selector.setToolTip(
            "Use OpenGL for the scan, integrated-data, and rocking plots "
            "after restarting orGUI."
        )
        form.addRow(self.opengl_selector)
        layout.addLayout(form)

        self.backend_notice = qt.QLabel(
            "Plot rendering changes take effect after restarting orGUI."
        )
        self.backend_notice.setWordWrap(True)
        layout.addWidget(self.backend_notice)

        self.remember_selector = qt.QCheckBox("remember this selection")
        self.remember_selector.setToolTip(
            "Use this style and rendering backend on future launches."
        )
        layout.addWidget(self.remember_selector)

        buttons = qt.QDialogButtonBox(
            qt.QDialogButtonBox.Apply | qt.QDialogButtonBox.Close
        )
        buttons.button(qt.QDialogButtonBox.Apply).clicked.connect(self.apply)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def apply(self):
        """Apply the style now and save a requested backend for next launch.

        Switching a visible Qt widget hierarchy between raster and OpenGL
        rendering can recreate its native top-level window.  The backend is
        therefore never changed in place; a remembered selection is applied
        during startup, before any plot widget is constructed.
        """
        style = qt.QStyleFactory.create(self.style_selector.currentText())
        if style is not None:
            self.application.setStyle(style)

        use_opengl = self.opengl_selector.isChecked()
        if use_opengl and self.plots:
            result = glutils.isOpenGLAvailable()
            if not result:
                qt.QMessageBox.critical(
                    self,
                    "OpenGL rendering is not available",
                    result.error,
                )
                return False

        backend_changed = any(
            state != use_opengl for state in self._plot_opengl_states
        )

        if self.remember_selector.isChecked():
            self.settings.setValue(STYLE_KEY, self.style_selector.currentText())
            self.settings.setValue(OPENGL_KEY, use_opengl)
            self.settings.sync()
            if backend_changed:
                self.backend_notice.setText(
                    "Plot rendering preference saved; restart orGUI to apply it."
                )
            else:
                self.backend_notice.setText(
                    "Plot rendering preference saved; the plots already use it."
                )
        elif backend_changed:
            self.backend_notice.setText(
                "Plot rendering was not changed. Select “remember this "
                "selection” and apply, then restart orGUI."
            )
        else:
            self.backend_notice.setText("Plot rendering is unchanged.")
        return True

    @staticmethod
    def _plot_uses_opengl(plot):
        """Return whether a silx plot is currently rendered with OpenGL."""
        return "opengl" in type(plot.getBackend()).__name__.lower()

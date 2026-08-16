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
        :param tuple plots: Plot widgets whose backend is changed together.
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
        self.opengl_selector.setChecked(
            bool(self.plots)
            and all(self._plot_uses_opengl(plot) for plot in self.plots)
        )
        self.opengl_selector.setToolTip(
            "Use OpenGL for the scan, integrated-data, and rocking plots."
        )
        form.addRow(self.opengl_selector)
        layout.addLayout(form)

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
        """Apply the selected style and rendering backend to this session."""
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

        for plot in self.plots:
            plot.setBackend("opengl" if use_opengl else "matplotlib")

        if self.remember_selector.isChecked():
            self.settings.setValue(STYLE_KEY, self.style_selector.currentText())
            self.settings.setValue(OPENGL_KEY, use_opengl)
            self.settings.sync()
        return True

    @staticmethod
    def _plot_uses_opengl(plot):
        """Return whether a silx plot is currently rendered with OpenGL."""
        return "opengl" in type(plot.getBackend()).__name__.lower()

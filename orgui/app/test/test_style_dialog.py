"""Tests for runtime Qt widget-style selection."""

from types import SimpleNamespace

from silx.gui import qt

from orgui.app.StyleDialog import StyleDialog
from orgui.app.ui_preferences import OPENGL_KEY, STYLE_KEY


class OpenGLBackend:
    """Minimal backend marker used to test backend selection."""


class MatplotlibBackend:
    """Minimal backend marker used to test backend selection."""


class UnavailableOpenGL:
    """Minimal unavailable OpenGL result returned by silx."""

    error = "OpenGL is unavailable"

    def __bool__(self):
        return False


class PlotStub:
    """Minimal silx plot interface used by the style dialog."""

    def __init__(self, backend=MatplotlibBackend):
        self.backend = backend()
        self.backend_names = []

    def getBackend(self):
        """Return the current test backend."""
        return self.backend

    def setBackend(self, name):
        """Record and expose the requested backend."""
        self.backend_names.append(name)
        self.backend = OpenGLBackend() if name == "opengl" else MatplotlibBackend()


def _settings(tmp_path):
    return qt.QSettings(str(tmp_path / "preferences.ini"), qt.QSettings.IniFormat)


def test_style_dialog_lists_and_applies_available_styles():
    """The dialog exposes the installed styles and applies a selection."""
    application = qt.QApplication.instance() or qt.QApplication([])
    original_style = application.style().objectName()
    dialog = StyleDialog(application=application)

    available_styles = list(qt.QStyleFactory.keys())
    assert dialog.style_selector.count() == len(available_styles)
    assert {
        dialog.style_selector.itemText(index)
        for index in range(dialog.style_selector.count())
    } == set(available_styles)

    selected_style = available_styles[-1]
    dialog.style_selector.setCurrentText(selected_style)
    dialog.apply()
    assert application.style().objectName().lower() == selected_style.lower()

    application.setStyle(qt.QStyleFactory.create(original_style))


def test_style_dialog_defers_backend_change_until_restart(monkeypatch, tmp_path):
    """A remembered OpenGL choice is saved without mutating live windows."""
    application = qt.QApplication.instance() or qt.QApplication([])
    plots = [PlotStub(OpenGLBackend) for _ in range(3)]
    settings = _settings(tmp_path)
    dialog = StyleDialog(
        application=application,
        plots=plots,
        settings=settings,
    )

    dialog.opengl_selector.setChecked(False)
    dialog.remember_selector.setChecked(True)
    assert dialog.apply() is True
    assert all(not plot.backend_names for plot in plots)
    assert settings.value(OPENGL_KEY, type=bool) is False
    assert "restart orGUI" in dialog.backend_notice.text()

    monkeypatch.setattr(
        "orgui.app.StyleDialog.glutils.isOpenGLAvailable",
        lambda: SimpleNamespace(error=""),
    )
    dialog.opengl_selector.setChecked(True)
    assert dialog.apply() is True
    assert all(not plot.backend_names for plot in plots)
    assert settings.value(OPENGL_KEY, type=bool) is True


def test_style_dialog_only_persists_when_requested(tmp_path):
    """A session-only choice does not overwrite stored UI preferences."""
    application = qt.QApplication.instance() or qt.QApplication([])
    settings = _settings(tmp_path)
    settings.setValue(STYLE_KEY, "previous-style")
    settings.setValue(OPENGL_KEY, True)
    dialog = StyleDialog(application=application, settings=settings)

    dialog.opengl_selector.setChecked(False)
    assert dialog.apply() is True
    assert settings.value(STYLE_KEY, type=str) == "previous-style"
    assert settings.value(OPENGL_KEY, type=bool) is True

    dialog.remember_selector.setChecked(True)
    assert dialog.apply() is True
    assert settings.value(STYLE_KEY, type=str) == dialog.style_selector.currentText()
    assert settings.value(OPENGL_KEY, type=bool) is False


def test_style_dialog_does_not_change_or_save_on_unavailable_opengl(
    monkeypatch, tmp_path
):
    """Unavailable OpenGL leaves all plots and preferences unchanged."""
    application = qt.QApplication.instance() or qt.QApplication([])
    plots = [PlotStub() for _ in range(3)]
    settings = _settings(tmp_path)
    dialog = StyleDialog(application=application, plots=plots, settings=settings)
    dialog.opengl_selector.setChecked(True)
    dialog.remember_selector.setChecked(True)
    monkeypatch.setattr(
        "orgui.app.StyleDialog.glutils.isOpenGLAvailable",
        UnavailableOpenGL,
    )
    monkeypatch.setattr(qt.QMessageBox, "critical", lambda *args: None)

    assert dialog.apply() is False
    assert all(not plot.backend_names for plot in plots)
    assert not settings.contains(STYLE_KEY)
    assert not settings.contains(OPENGL_KEY)

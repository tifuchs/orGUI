"""Tests for versioned portable GUI preferences."""

from silx.gui import qt

from orgui.app.ui_preferences import (
    APPLICATION_NAME,
    OPENGL_KEY,
    ORGANIZATION_NAME,
    STYLE_KEY,
    apply_saved_ui_preferences,
    configure_application_identity,
)


def _settings(tmp_path):
    return qt.QSettings(str(tmp_path / "preferences.ini"), qt.QSettings.IniFormat)


def test_configure_application_identity_uses_portable_names():
    """QSettings uses the stable application identity on every platform."""
    application = qt.QApplication.instance() or qt.QApplication([])

    configure_application_identity(application)

    assert application.organizationName() == ORGANIZATION_NAME
    assert application.applicationName() == APPLICATION_NAME


def test_saved_preferences_apply_valid_style_and_opengl(tmp_path):
    """Saved versioned values configure startup before plots are created."""
    application = qt.QApplication.instance() or qt.QApplication([])
    original_style = application.style().objectName()
    settings = _settings(tmp_path)
    style = next(iter(qt.QStyleFactory.keys()))
    settings.setValue(STYLE_KEY, style)
    settings.setValue(OPENGL_KEY, True)

    assert apply_saved_ui_preferences(application, settings=settings) is True
    assert application.style().objectName().lower() == style.lower()

    application.setStyle(qt.QStyleFactory.create(original_style))


def test_invalid_or_missing_preferences_leave_defaults_unchanged(tmp_path):
    """Unknown styles and absent backend settings safely retain defaults."""
    application = qt.QApplication.instance() or qt.QApplication([])
    original_style = application.style().objectName()
    settings = _settings(tmp_path)
    settings.setValue(STYLE_KEY, "not-an-installed-style")

    assert apply_saved_ui_preferences(application, settings=settings) is False
    assert application.style().objectName() == original_style


def test_opengl_command_line_override_wins_over_saved_choice(tmp_path):
    """The explicit command-line OpenGL request overrides a saved false value."""
    application = qt.QApplication.instance() or qt.QApplication([])
    settings = _settings(tmp_path)
    settings.setValue(OPENGL_KEY, False)

    assert apply_saved_ui_preferences(application, settings=settings) is False
    assert apply_saved_ui_preferences(
        application, force_opengl=True, settings=settings
    ) is True

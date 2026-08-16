"""Versioned, per-user preferences for the graphical interface."""

from silx.gui import qt


ORGANIZATION_NAME = "orSXRD_Soft"
APPLICATION_NAME = "orGUI"
STYLE_KEY = "ui/v1/style"
OPENGL_KEY = "ui/v1/opengl"


def configure_application_identity(application):
    """Configure the portable QSettings identity for an application."""
    application.setOrganizationName(ORGANIZATION_NAME)
    application.setApplicationName(APPLICATION_NAME)


def apply_saved_ui_preferences(application, force_opengl=False, settings=None):
    """Apply saved GUI preferences before the main window is constructed.

    :param QApplication application: Qt application receiving the widget style.
    :param bool force_opengl: Whether a command-line request overrides the
        saved OpenGL preference for this launch.
    :param QSettings settings: Optional settings store, used by tests.
    :returns: Whether newly constructed plots should use OpenGL.
    :rtype: bool
    """
    settings = settings if settings is not None else qt.QSettings()
    style_name = settings.value(STYLE_KEY, "", type=str)
    if style_name:
        style = qt.QStyleFactory.create(style_name)
        if style is not None:
            application.setStyle(style)

    return force_opengl or settings.value(OPENGL_KEY, False, type=bool)

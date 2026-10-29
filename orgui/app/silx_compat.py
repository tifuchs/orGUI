"""Compatibility fixes for silx plotting used by the application."""

from functools import wraps
from io import BytesIO

from silx.gui import qt
from silx.gui.plot import PlotWindow


def _install_png_buffer_compat():
    """Allow silx plot copy/print actions to export PNGs to Python buffers."""
    original = PlotWindow.saveGraph
    if getattr(original, "_orgui_png_buffer_compat", False) is True:
        return

    @wraps(original)
    def _save_graph(self, filename, fileFormat=None, dpi=None):
        try:
            return original(self, filename, fileFormat=fileFormat, dpi=dpi)
        except TypeError:
            if fileFormat != "png" or not isinstance(filename, BytesIO):
                raise

        # silx 3.1.3 passes BytesIO directly to QPixmap.save, although Qt
        # requires a QIODevice. Keep working silx implementations unchanged.
        # Grab the same central widget so the snapshot includes the colorbar.
        image = self.centralWidget().grab().toImage()
        if dpi is not None:
            # Qt stores integer dots per metre; one inch is 0.0254 metres.
            image.setDotsPerMeterX(round(dpi / 0.0254))
            image.setDotsPerMeterY(round(dpi / 0.0254))
        buffer = qt.QBuffer()
        if not buffer.open(qt.QIODevice.WriteOnly):
            return False
        try:
            if not image.save(buffer, "PNG"):
                return False
            filename.write(bytes(buffer.data()))
            return True
        finally:
            buffer.close()

    _save_graph._orgui_png_buffer_compat = True
    PlotWindow.saveGraph = _save_graph

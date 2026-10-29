"""Regression coverage for silx PNG export with PySide6."""

from io import BytesIO
from unittest.mock import Mock

import pytest
from silx.gui import qt
from silx.gui.plot import Plot1D, Plot2D, PlotWindow
from silx.gui.plot.actions.io import CopyAction, _plotAsPNG

from orgui.app.silx_compat import _install_png_buffer_compat


@pytest.fixture
def qapp():
    """Keep the Qt application alive throughout each test."""
    return qt.QApplication.instance() or qt.QApplication([])


@pytest.fixture
def compatible_plots(monkeypatch):
    """Install the compatibility wrapper without leaking it to other tests."""
    monkeypatch.setattr(PlotWindow, "saveGraph", PlotWindow.saveGraph)
    _install_png_buffer_compat()


@pytest.mark.parametrize("plot_class", [PlotWindow, Plot1D, Plot2D])
def test_copy_plot_png(qapp, compatible_plots, monkeypatch, plot_class):
    """The real silx copy action delivers a valid image to the clipboard."""
    plot = plot_class(backend="matplotlib")
    clipboard = Mock()
    monkeypatch.setattr(qt.QApplication, "clipboard", lambda: clipboard)
    try:
        plot.resize(640, 480)
        if plot_class is Plot2D:
            plot.addImage([[1., 2.], [3., 4.]])
        else:
            plot.addCurve([0, 1], [1, 2])
        plot.show()
        qapp.processEvents()
        # This helper is shared by the copy and print actions.
        png = _plotAsPNG(plot)
        assert png.startswith(b"\x89PNG\r\n\x1a\n")
        image = qt.QImage.fromData(png, "PNG")
        assert not image.isNull()
        assert image.size() == plot.centralWidget().grab().size()
        CopyAction(plot).copyPlot()
        clipboard.setImage.assert_called_once()
        assert not clipboard.setImage.call_args.args[0].isNull()
    finally:
        plot.close()


def test_png_buffer_dpi(qapp, compatible_plots):
    """PNG buffers preserve explicit resolution metadata."""
    plot = PlotWindow(backend="matplotlib")
    try:
        output = BytesIO()
        assert plot.saveGraph(output, fileFormat="png", dpi=150)
        image = qt.QImage.fromData(output.getvalue(), "PNG")
        assert not image.isNull()
        assert image.dotsPerMeterX() == pytest.approx(150 / 0.0254, abs=1)
        assert image.dotsPerMeterY() == pytest.approx(150 / 0.0254, abs=1)
    finally:
        plot.close()


@pytest.mark.parametrize(
    "destination, format_", [("plot.png", "png"), (BytesIO(), "svg")],
)
def test_other_exports_delegate(monkeypatch, destination, format_):
    """File paths and non-PNG output retain the original save behavior."""
    original = Mock(return_value=False)
    monkeypatch.setattr(PlotWindow, "saveGraph", original)
    _install_png_buffer_compat()
    plot = object()
    assert PlotWindow.saveGraph(plot, destination, format_, 120) is False
    original.assert_called_once_with(
        plot, destination, fileFormat=format_, dpi=120,
    )


def test_working_png_export_is_preserved(monkeypatch):
    """A silx implementation that already supports buffers is used as-is."""
    original = Mock(return_value=True)
    monkeypatch.setattr(PlotWindow, "saveGraph", original)
    _install_png_buffer_compat()
    output = BytesIO()
    plot = object()
    assert PlotWindow.saveGraph(plot, output, "png") is True
    original.assert_called_once_with(plot, output, fileFormat="png", dpi=None)
    installed = PlotWindow.saveGraph
    _install_png_buffer_compat()
    assert PlotWindow.saveGraph is installed


def test_unrelated_type_error_is_preserved(monkeypatch):
    """Errors outside the PNG buffer path are not swallowed."""
    monkeypatch.setattr(PlotWindow, "saveGraph", Mock(side_effect=TypeError))
    _install_png_buffer_compat()
    with pytest.raises(TypeError):
        PlotWindow.saveGraph(object(), "plot.png", "png")

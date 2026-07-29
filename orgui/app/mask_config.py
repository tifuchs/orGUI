"""CLI-safe detector mask configuration and loading helpers."""

from __future__ import annotations

from dataclasses import dataclass
import configparser
import logging
import os
from pathlib import Path

import fabio
import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class PixelRepairSettings:
    """Pixel-repair settings loaded from ``[Mask.PixelRepair]``.

    :param bool enabled:
        Whether masked-pixel repair is requested.
    :param int max_component_pixels:
        Maximum connected masked defect size in pixels.
    :param int max_span:
        Maximum row or column span of a repaired defect in pixels.
    :param int radius:
        Neighbor-search radius in pixels.
    :param int min_valid_neighbors:
        Minimum original-valid neighbors required for interpolation.
    :param bool use_pyfai_gaps:
        Whether analytical pyFAI detector gaps must be excluded from repair.
    :param int gap_size_px:
        Fallback detector-gap width in pixels when metadata are incomplete.
    """

    enabled: bool = False
    max_component_pixels: int = 4
    max_span: int = 3
    radius: int = 2
    min_valid_neighbors: int = 6
    use_pyfai_gaps: bool = True
    gap_size_px: int = 6


@dataclass
class MaskSettings:
    """Detector-mask settings loaded from ``[Mask]``.

    :param str mask:
        Optional mask path. Nonzero or ``True`` pixels are invalid.
    :param PixelRepairSettings pixel_repair:
        Optional repair settings. Missing ``[Mask.PixelRepair]`` disables repair.
    """

    mask: str | None = None
    pixel_repair: PixelRepairSettings | None = None


def _warn_default(option, value, default, exc):
    logger.warning(
        "Invalid mask configuration value %s=%r; using default %r. %s",
        option,
        value,
        default,
        exc,
    )


def _get_int(section, option, default, minimum=0):
    raw = section.get(option, fallback=None)
    if raw is None:
        return default
    try:
        value = int(raw)
        if value < minimum:
            raise ValueError(f"value must be >= {minimum}")
        return value
    except Exception as exc:
        _warn_default(option, raw, default, exc)
        return default


def _get_bool(section, option, default):
    raw = section.get(option, fallback=None)
    if raw is None:
        return default
    try:
        parser = configparser.ConfigParser()
        parser.add_section("value")
        parser.set("value", option, raw)
        return parser.getboolean("value", option)
    except Exception as exc:
        _warn_default(option, raw, default, exc)
        return default


def parse_mask_settings(config, base_dir=None):
    """Parse detector-mask settings from an existing config parser.

    :param configparser.ConfigParser config:
        Parsed application configuration.
    :param str base_dir:
        Directory used to resolve relative mask paths.
    :returns:
        Parsed mask settings. Missing optional sections are valid.
    :rtype: MaskSettings
    """

    if "Mask" not in config:
        return MaskSettings()

    mask_path = config["Mask"].get("mask", fallback=None)
    if mask_path is not None:
        mask_path = mask_path.strip() or None
    if mask_path and base_dir and not os.path.isabs(mask_path):
        mask_path = os.path.join(base_dir, mask_path)

    repair = None
    if "Mask.PixelRepair" in config:
        section = config["Mask.PixelRepair"]
        defaults = PixelRepairSettings()
        repair = PixelRepairSettings(
            enabled=_get_bool(section, "enabled", defaults.enabled),
            max_component_pixels=_get_int(
                section,
                "max_component_pixels",
                defaults.max_component_pixels,
                minimum=1,
            ),
            max_span=_get_int(section, "max_span", defaults.max_span, minimum=1),
            radius=_get_int(section, "radius", defaults.radius, minimum=1),
            min_valid_neighbors=_get_int(
                section,
                "min_valid_neighbors",
                defaults.min_valid_neighbors,
                minimum=1,
            ),
            use_pyfai_gaps=_get_bool(
                section,
                "use_pyfai_gaps",
                defaults.use_pyfai_gaps,
            ),
            gap_size_px=_get_int(
                section,
                "gap_size_px",
                defaults.gap_size_px,
                minimum=0,
            ),
        )
    return MaskSettings(mask=mask_path, pixel_repair=repair)


class MaskManager:
    """Manage detector mask state without depending on Qt widgets."""

    def __init__(self):
        self.settings = MaskSettings()
        self._mask = None

    def load_config(self, filename):
        """Load mask settings and optional mask data from a config file.

        :param str filename:
            Config file path.
        :returns:
            Loaded boolean mask, or ``None`` if no mask is configured.
        :rtype: numpy.ndarray or None
        """

        config = configparser.ConfigParser()
        config.read(filename)
        self.settings = parse_mask_settings(config, Path(filename).resolve().parent)
        if self.settings.mask:
            self.load_mask(self.settings.mask)
        else:
            self._mask = None
        return self._mask

    def load_mask(self, filename):
        """Load an invalid-pixel mask from disk.

        :param str filename:
            Path to a FabIO-readable 2D mask image, for example EDF or NumPy
            ``.npy``. Nonzero pixels are invalid.
        :returns:
            Boolean mask with ``True`` for invalid pixels.
        :rtype: numpy.ndarray
        """

        path = Path(filename)
        with fabio.open(str(path)) as image:
            mask = image.data
        if mask.ndim != 2:
            raise TypeError(
                f"Path {path} identifies a {mask.ndim}D array, but a 2D mask "
                "is expected"
            )
        if mask.dtype.kind == "b":
            mask = mask.astype(np.int8)
        elif mask.dtype.kind not in "fui":
            raise TypeError(
                f"Path {path} identifies a {mask.dtype.kind!r}-kind array, "
                "but a numerical mask is expected"
            )
        self.set_mask(mask)
        self.settings.mask = str(path)
        return self._mask

    def set_mask(self, mask):
        """Set the active invalid-pixel mask.

        :param numpy.ndarray mask:
            Array where nonzero or ``True`` pixels are invalid.
        """

        if mask is None:
            self._mask = None
            return
        arr = np.asarray(mask)
        if arr.ndim != 2:
            raise ValueError("Mask must be a 2D array")
        self._mask = np.ascontiguousarray(arr > 0, dtype=bool)

    def get_mask(self, shape=None):
        """Return the active invalid-pixel mask.

        :param tuple shape:
            Optional expected detector image shape.
        :returns:
            Boolean mask or ``None``.
        :rtype: numpy.ndarray or None
        """

        if self._mask is None:
            return None
        if shape is not None and tuple(self._mask.shape) != tuple(shape):
            logger.warning(
                "Configured detector mask shape %s does not match image shape %s; "
                "ignoring mask.",
                self._mask.shape,
                shape,
            )
            return None
        return self._mask

    def set_pixel_repair_settings(self, settings):
        """Set pixel-repair settings explicitly."""

        self.settings.pixel_repair = settings

    def repair_enabled(self, cpp_available=True):
        """Return whether pixel repair should run for this process."""

        repair = self.settings.pixel_repair
        if repair is None or not repair.enabled:
            return False
        if not cpp_available:
            logger.warning(
                "Pixel repair is enabled in config but requires the C++ ROI "
                "backend; continuing without repair."
            )
            return False
        return True

    def pyfai_gap_intervals(self, detector):
        """Return analytical pyFAI detector gap intervals.

        The returned arrays are ``int32`` ``[start, stop)`` intervals for rows
        and columns. pyFAI Dectris masks are reduced to full-width or
        full-height masked stripes, so ordinary bad-pixel masks do not become
        detector gaps.

        :param detector:
            pyFAI detector object.
        :returns:
            ``(row_intervals, column_intervals)``.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """

        repair = self.settings.pixel_repair or PixelRepairSettings()
        if not repair.use_pyfai_gaps or detector is None:
            return _empty_intervals(), _empty_intervals()
        mask = None
        if hasattr(detector, "calc_mask"):
            try:
                mask = detector.calc_mask()
            except Exception:
                logger.warning("Unable to calculate pyFAI detector gap mask.")
        if mask is None:
            mask = getattr(detector, "mask", None)
        if mask is None:
            return _metadata_gap_intervals(detector, repair.gap_size_px)
        mask = np.asarray(mask) > 0
        if mask.ndim != 2:
            return _metadata_gap_intervals(detector, repair.gap_size_px)
        row_intervals = _runs(np.all(mask, axis=1))
        col_intervals = _runs(np.all(mask, axis=0))
        if row_intervals.size == 0 and col_intervals.size == 0:
            return _metadata_gap_intervals(detector, repair.gap_size_px)
        return row_intervals, col_intervals


def _empty_intervals():
    return np.empty((0, 2), dtype=np.int32)


def _runs(flags):
    starts = []
    start = None
    for idx, value in enumerate(np.asarray(flags, dtype=bool)):
        if value and start is None:
            start = idx
        elif not value and start is not None:
            starts.append((start, idx))
            start = None
    if start is not None:
        starts.append((start, len(flags)))
    if not starts:
        return _empty_intervals()
    return np.asarray(starts, dtype=np.int32)


def _metadata_gap_intervals(detector, fallback_gap_size):
    module_size = getattr(detector, "module_size", None)
    if module_size is None:
        module_size = getattr(detector, "MODULE_SIZE", None)
    shape = getattr(detector, "shape", None)
    if shape is None:
        shape = getattr(detector, "max_shape", None)
    if module_size is None or shape is None:
        return _empty_intervals(), _empty_intervals()

    module_gap = getattr(detector, "MODULE_GAP", None)
    if module_gap is None:
        gap_rows = gap_cols = int(fallback_gap_size)
    else:
        gap_rows = int(module_gap[0]) if module_gap[0] else int(fallback_gap_size)
        gap_cols = int(module_gap[1]) if module_gap[1] else int(fallback_gap_size)

    return (
        _module_gap_runs(int(shape[0]), int(module_size[0]), gap_rows),
        _module_gap_runs(int(shape[1]), int(module_size[1]), gap_cols),
    )


def _module_gap_runs(axis_size, module_pixels, gap_pixels):
    if axis_size <= 0 or module_pixels <= 0 or gap_pixels <= 0:
        return _empty_intervals()
    intervals = []
    start = module_pixels
    while start < axis_size:
        stop = min(axis_size, start + gap_pixels)
        intervals.append((start, stop))
        start = stop + module_pixels
    if not intervals:
        return _empty_intervals()
    return np.asarray(intervals, dtype=np.int32)

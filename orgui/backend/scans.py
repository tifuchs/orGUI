# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"


import numpy as np
from abc import ABC, abstractmethod
from dataclasses import dataclass
import importlib
import os
from pathlib import Path
import runpy


class Scan(ABC):
    """Base class for every backend.

    All methods in this class must be populated for a working backend.
    """

    @abstractmethod
    def __init__(self, hdffilepath_orNode=None, scanno=None):
        """Constructor of the class.

        hdffilepath_orNode is an either an object selected from the file browser,
        where data scan be accessed, or a file path. Example on how to open

        .. code-block:: python
            if isinstance(hdffilepath_orNode, str): # check if this is a file path or already an opened h5node
                with silx.io.open(hdffilepath_orNode) as f: # file path, so we have to open the file manually
                    self.th = f[scanno ,'th'] # read data from the file, see below for the required entries
                    ...
            else:       # file already open, so just read the data from the file
                hdffilepath = hdffilepath_orNode.local_filename # to get the filename
                f = hdffilepath_orNode.file
                self.th = f[scanno ,'th']
                ...

        A backend must populate the following entries
        *  axisname : 'th' or 'mu', defines the scan axis
        *  axis : data value of the axis (should be a numpy array)
        *  title : a meaningful title of the scan (e.g. "ascan th 0 90 90 1")
        *  th : motor value (as array)
        *  om : = -th (just populate the value)
        *  mu : motor value
        with actual data
        """  # noqa: E501
        self.axisname = None  # either "th" or "mu", this defines the scan axis
        self.axis = None  # value of either "th" or "mu"

        self.th = 0.0  # or self.mu, depending on scanaxis
        self.omega = -self.th  # = -1*th
        self.title = "generic_scan"
        # for mu-scan you must provide a value for omega/theta

        # optional: provide a unique identifier, which is used as key in h5 database to identify  # noqa: E501
        # the dataset
        # self.name = "identifier"

    @property
    def auxillary_counters(self):
        """Optional: provide a list of counters or motor names, that should
        be copied into the orGUI data base for further processing.

        e.g. return ['exposure_time', 'elapsed_time','time', 'srcur', 'mondio', 'epoch']

        after each integration, orGUI will search for these counter names in the Scan
        object and copy the entries into the database.
        """
        return []

    @classmethod
    @abstractmethod
    def parse_h5_node(cls, node):
        pass

    @classmethod
    def listScans(cls, h5file):
        """Optional: the scans an opened HDF5 file contains.

        Used by the segmented ("interlaced") scan loader to fill its selection
        dialog, which has to know about every scan in the file before any scan
        object exists.

        Return whatever your ``__init__`` accepts as its scan identifier, which
        is the same thing :meth:`parse_h5_node` reports for a single node:

        .. code-block:: python

            return [1, 2, 10]                  # scan numbers
            return ["ascan_12", "dscan_3"]     # scan names

        Optionally pair each identifier with a label shown in the dialog:

        .. code-block:: python

            return [(1, "ascan th 0 90 90 1"), (2, "dscan mu 0 1 20 1")]

        Do not use ``float`` for BLISS style ``"<scan>.<subscan>"`` names:
        ``1.1`` and ``1.10`` are the same float. Integers, strings and numpy
        integers are all fine, the identifier is passed on untouched and is
        only converted to ``str`` for display.

        Implementing this is optional. Without it orGUI calls
        :meth:`parse_h5_node` on every entry of the file root instead, which
        needs no extra code in the backend but is slower, and which relies on
        ``parse_h5_node`` raising for entries that are not scans.

        :param h5file: an open h5py-like file object.
        :returns: identifiers, or ``(identifier, label)`` pairs.
        :rtype: list
        """
        raise NotImplementedError(f"{cls.__name__} does not implement listScans")

    @abstractmethod
    def __len__(self):
        """returns the number of entries in the scan."""
        raise NotImplementedError()

    @abstractmethod
    def get_raw_img(self, i):
        """This should return a populated image class such as h5_Image (see below)

        Only the image data is required as numpy array, accessible as h5_image.img
        """
        raise NotImplementedError()

    def exposure_angle_bounds(self, config, fallback="stationary"):
        """Return exposure bounds in Vlieg order and radians.

        :param ConfigData config:
            Active orGUI configuration providing fixed ``mu``, ``chi``, and
            ``phi`` values.
        :param str fallback:
            ``"stationary"`` or the explicitly requested ``"midpoint"``.
        :returns:
            Array shaped ``(frames, 2, 4)`` ordered
            ``alpha, omega, chi, phi``.
        """
        return scan_exposure_angle_bounds(self, config, fallback=fallback)


class h5_Image:
    def __init__(self, data, variance=None, processing_provenance=None):
        """Only the image data is required as numpy array.
        motors and counters do not need to be populated.
        """
        self.img = data
        self.variance = variance
        self.processing_provenance = dict(processing_provenance or {})
        self.motors = dict()
        self.counters = dict()


def load_scan_backend_file(filename, class_name=None):
    """Load one scan backend class from an existing Python backend file.

    :param filename:
        Backend Python file understood by orGUI's backend loader.
    :param str class_name:
        Optional qualified class name required by a serialized scan reference.
    :returns:
        ``(namespace_name, scan_class)``.
    :raises ValueError:
        If no unique matching :class:`Scan` implementation is present.
    """
    filename = str(Path(filename).absolute())
    namespace = runpy.run_path(filename)
    candidates = []
    seen = set()
    for name, value in namespace.items():
        try:
            is_scan = issubclass(value, Scan) and value is not Scan
        except TypeError:
            continue
        if not is_scan or id(value) in seen:
            continue
        if class_name is not None and value.__qualname__ != class_name:
            continue
        seen.add(id(value))
        candidates.append((name, value))
    qualifier = f" named {class_name!r}" if class_name is not None else ""
    if not candidates:
        raise ValueError(
            f"Found no Scan class{qualifier} in backend file {filename}"
        )
    if len(candidates) > 1:
        raise ValueError(
            f"Found more than one Scan class{qualifier} in backend file "
            f"{filename}. Only one is permitted"
        )
    name, scan_class = candidates[0]
    scan_class._orgui_backend_file = filename
    return name, scan_class


class SimulationScan(Scan):
    def __init__(self, detshape, axismin, axismax, points, axis="th", fixed=0.0):
        self.shape = detshape
        self.axisname = axis
        self.axis = np.linspace(axismin, axismax, points)
        if axis == "th":
            self.th = self.axis
            self.omega = -1 * self.th
            self.mu = fixed
        elif axis == "mu":
            self.mu = self.axis
            self.th = fixed
            self.omega = -1 * self.th
        else:
            raise ValueError(f"{axis} is not an implemented scan axis.")
        self.nopoints = points

        self.images = np.zeros((points, *detshape)) + 10
        self.title = f"sim ascan {self.axisname} {axismin} {axismax} {points}"

    def __len__(self):
        return self.nopoints

    def get_raw_img(self, i):
        return h5_Image(self.images[i])

    def set_raw_img(self, i, data):  # for intensity simulation in the future.
        self.images[i] = data

    def parse_h5_node(cls, node):  # unused
        pass


def _frame_values(values, count, name):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim == 0 or values.size == 1:
        return np.full(count, float(values.reshape(-1)[0]))
    if values.ndim == 1 and values.size == count:
        return values
    raise ValueError(f"{name} must be scalar or contain one value per frame")


def _midpoint_bounds(values):
    values = np.unwrap(np.asarray(values, dtype=np.float64))
    if values.size == 1:
        return np.repeat(values[:, None], 2, axis=1)
    edges = np.empty(values.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = values[0] - 0.5 * (values[1] - values[0])
    edges[-1] = values[-1] + 0.5 * (values[-1] - values[-2])
    return np.column_stack((edges[:-1], edges[1:]))


def scan_exposure_angle_bounds(scan, config, fallback="stationary"):
    """Return centralized exposure bounds for any loaded orGUI scan.

    Scan motor positions are degrees. Fixed diffractometer values from
    :class:`ConfigData` and the returned bounds are radians.
    """
    if fallback not in {"stationary", "midpoint"}:
        raise ValueError("fallback must be 'stationary' or 'midpoint'")
    if hasattr(scan, "subscans"):
        bounds = np.concatenate(
            [
                scan_exposure_angle_bounds(child, config, fallback=fallback)
                for child in scan.subscans
            ]
        )
        indices = getattr(scan, "indices", None)
        if indices is not None:
            bounds = bounds[np.asarray(indices, dtype=np.int64)]
        return np.ascontiguousarray(bounds)
    count = len(scan)
    if hasattr(scan, "mu"):
        alpha = np.deg2rad(_frame_values(scan.mu, count, "alpha"))
    else:
        alpha = _frame_values(config.mu, count, "alpha")
    omega = np.deg2rad(
        _frame_values(getattr(scan, "omega", 0.0), count, "omega")
    )
    centers = np.column_stack(
        (
            alpha,
            omega,
            np.full(count, config.chi),
            np.full(count, config.phi),
        )
    )
    bounds = np.repeat(centers[:, None, :], 2, axis=1)
    explicit = getattr(scan, "exposure_bounds_rad", None)
    if explicit is not None:
        explicit = np.asarray(explicit, dtype=np.float64)
        if explicit.shape != bounds.shape:
            raise ValueError("scan exposure_bounds_rad has an invalid shape")
        return np.ascontiguousarray(explicit)
    for axis, name in enumerate(("alpha", "omega", "chi", "phi")):
        explicit_axis = getattr(scan, f"{name}_bounds_rad", None)
        if explicit_axis is None:
            continue
        explicit_axis = np.asarray(explicit_axis, dtype=np.float64)
        if explicit_axis.shape != (count, 2):
            raise ValueError(f"scan {name}_bounds_rad has an invalid shape")
        bounds[:, :, axis] = explicit_axis
    if fallback == "midpoint":
        angle_names = ("alpha", "omega", "chi", "phi")
        for axis in range(4):
            explicit_axis = getattr(
                scan, f"{angle_names[axis]}_bounds_rad", None
            )
            if explicit_axis is None:
                bounds[:, :, axis] = _midpoint_bounds(centers[:, axis])
    return np.ascontiguousarray(bounds)


@dataclass(frozen=True)
class ScanReference:
    """Serializable reference capable of reopening an orGUI scan."""

    kind: str
    module: str
    class_name: str
    parameters: dict
    source_fingerprints: tuple[dict, ...] = ()

    @staticmethod
    def _fingerprint(path):
        path = Path(path).absolute()
        stat = path.stat()
        return {
            "path": str(path),
            "size": int(stat.st_size),
            "mtime_ns": int(stat.st_mtime_ns),
        }

    @classmethod
    def from_scan(cls, scan):
        """Create a central reference for a currently loaded scan."""
        scan_class = type(scan)
        module = scan_class.__module__
        class_name = scan_class.__qualname__
        backend_file = getattr(scan_class, "_orgui_backend_file", None)
        if backend_file is None and module in {"<run_path>", "__main__"}:
            initializer = getattr(scan_class, "__init__", None)
            backend_file = getattr(initializer, "__globals__", {}).get(
                "__file__"
            )
        if backend_file is not None:
            backend_file = str(Path(backend_file).absolute())
            if not Path(backend_file).is_file():
                backend_file = None
        if isinstance(scan, SimulationScan):
            fixed = scan.mu if scan.axisname == "th" else scan.th
            return cls(
                "simulation",
                module,
                class_name,
                {
                    "detshape": list(scan.shape),
                    "axismin": float(scan.axis[0]),
                    "axismax": float(scan.axis[-1]),
                    "points": len(scan),
                    "axis": scan.axisname,
                    "fixed": float(np.asarray(fixed).reshape(-1)[0]),
                },
            )
        if module.endswith("interlacedScanLoader"):
            return cls(
                "interlaced",
                module,
                class_name,
                {
                    "scans": [
                        reference.to_dict()
                        for reference in map(cls.from_scan, scan.subscans)
                    ],
                    "sort": bool(scan.sort),
                    "axis": scan.axisname,
                },
            )
        if module.endswith("universalScanLoader"):
            paths = [
                str(Path(directory, filename).absolute())
                for directory, filenames in [scan.inpath]
                for filename in filenames
            ]
            return cls(
                "manual_images",
                module,
                class_name,
                {
                    "filename": str(Path(scan.filename).absolute()),
                    "axis": scan.axisname,
                    "axismin": float(np.asarray(scan.axis)[0]),
                    "axismax": float(np.asarray(scan.axis)[-1]),
                    "fixed": float(
                        np.asarray(
                            scan.mu if scan.axisname == "th" else scan.th
                        ).reshape(-1)[0]
                    ),
                },
                tuple(cls._fingerprint(path) for path in paths),
            )
        source = getattr(scan, "hdffilepath_orNode", None)
        if source is not None and not isinstance(source, str | os.PathLike):
            source = getattr(source, "local_filename", None)
        if source is None:
            source = getattr(scan, "fastscan_specfile", None)
        if source is None:
            source = getattr(scan, "filename", None)
        if source is None:
            filenames = getattr(scan, "filenames", None)
            if filenames:
                source = filenames[0]
        if source is None:
            raise ValueError(
                f"{module}.{class_name} does not expose a reloadable source"
            )
        source = str(Path(source).absolute())
        scan_number = getattr(scan, "scanno", None)
        if scan_number is None:
            scan_number = getattr(scan, "scanno1", None)
            if (
                isinstance(scan_number, str)
                and scan_number.endswith(".1")
                and hasattr(scan, "scanno2")
            ):
                scan_number = scan_number[:-2]
        if scan_number is None:
            scan_name = getattr(scan, "scanname", None)
            if isinstance(scan_name, str):
                suffix = scan_name.rsplit("_", 1)[-1]
                try:
                    scan_number = int(suffix.split(".", 1)[0])
                except ValueError:
                    scan_number = None
        parameters = {"source": source, "scan_number": scan_number}
        if backend_file is not None:
            parameters["backend_file"] = backend_file
        constructor_args = getattr(scan, "_scan_reference_args", None)
        if constructor_args is not None:
            parameters["constructor_args"] = list(constructor_args)
        if hasattr(scan, "offsetindex") and (
            int(scan.offsetindex) > 0
            or len(scan) < int(getattr(scan, "scandatapoints", len(scan)))
        ):
            parameters["slice"] = [
                int(scan.offsetindex),
                int(scan.offsetindex) + len(scan),
            ]
        fingerprints = [cls._fingerprint(source)]
        for filename in getattr(scan, "filenames", ()):
            candidate = Path(filename).absolute()
            if candidate.is_file() and str(candidate) != source:
                fingerprints.append(cls._fingerprint(candidate))
        if backend_file is not None and backend_file != source:
            fingerprints.append(cls._fingerprint(backend_file))
        return cls(
            "backend_file" if backend_file is not None else "backend",
            module,
            class_name,
            parameters,
            tuple(fingerprints),
        )

    def to_dict(self):
        """Return the JSON-compatible scan reference."""
        return {
            "kind": self.kind,
            "module": self.module,
            "class_name": self.class_name,
            "parameters": self.parameters,
            "source_fingerprints": list(self.source_fingerprints),
        }

    @classmethod
    def from_dict(cls, values):
        """Build a scan reference from JSON-compatible values."""
        values = dict(values)
        values["source_fingerprints"] = tuple(
            values.get("source_fingerprints", ())
        )
        return cls(**values)

    def verify(self):
        """Verify that referenced source files have not changed."""
        for expected in self.source_fingerprints:
            actual = self._fingerprint(expected["path"])
            if actual != expected:
                raise RuntimeError(
                    f"Scan source changed after job preparation: "
                    f"{expected['path']}"
                )

    def open(self):
        """Reopen this scan through the existing central backend class."""
        self.verify()
        if self.kind == "backend_file":
            _, scan_class = load_scan_backend_file(
                self.parameters["backend_file"], self.class_name
            )
        elif self.module == "<run_path>":
            # Jobs created before backend-file provenance was recorded can
            # still run in the GUI process where that backend remains loaded.
            from . import backends

            candidates = {
                id(value): value
                for value in backends.fscans.values()
                if value.__qualname__ == self.class_name
            }
            if len(candidates) != 1:
                raise RuntimeError(
                    "This job references a transient custom backend but does "
                    "not record its Python file. Load the same backend in "
                    "orGUI and retry, or prepare a new reconstruction job."
                )
            scan_class = next(iter(candidates.values()))
        else:
            scan_class = importlib.import_module(self.module)
            for component in self.class_name.split("."):
                scan_class = getattr(scan_class, component)
        if self.kind == "simulation":
            return scan_class(**self.parameters)
        if self.kind == "interlaced":
            scans = [
                ScanReference.from_dict(value).open()
                for value in self.parameters["scans"]
            ]
            return scan_class(
                scans, self.parameters["sort"], self.parameters["axis"]
            )
        if self.kind == "manual_images":
            scan = scan_class(self.parameters["filename"])
            scan.set_axis(
                self.parameters["axismin"],
                self.parameters["axismax"],
                self.parameters["axis"],
                self.parameters["fixed"],
            )
            return scan
        constructor_args = self.parameters.get("constructor_args")
        if constructor_args is not None:
            scan = scan_class(*constructor_args)
        else:
            number = self.parameters.get("scan_number")
            scan = (
                scan_class(self.parameters["source"], number)
                if number is not None
                else scan_class(self.parameters["source"])
            )
        if "slice" in self.parameters:
            scan = scan.slice(*self.parameters["slice"])
        return scan

# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#
# ###########################################################################*/
"""Structured application configuration snapshots for orGUI databases."""

from __future__ import annotations

import configparser
import datetime
from dataclasses import dataclass, field

import h5py
import numpy as np
from silx.io.dictdump import dicttonx, nxtodict

from ..datautils.xrayutils import CTRcalc, DetectorCalibration, HKLVlieg

SCHEMA_VERSION = 1


def _as_text(value):
    if isinstance(value, bytes):
        return value.decode()
    return str(value)


def _none_to_nan(value):
    return np.nan if value is None else value


def _string_array(strings):
    encoded = [str(value).encode("utf-8") for value in strings]
    if not encoded:
        return np.empty((0, 0), dtype=np.uint8)
    width = max(len(value) for value in encoded)
    data = np.zeros((len(encoded), width), dtype=np.uint8)
    for index, value in enumerate(encoded):
        data[index, : len(value)] = np.frombuffer(value, dtype=np.uint8)
    return data


def _read_string_array(data):
    array = np.asarray(data)
    if array.dtype == np.uint8 and array.ndim == 2:
        return [bytes(row[row != 0]).decode("utf-8") for row in array]
    return [_as_text(value) for value in data]


def detector_to_nxdict(detector):
    """Serialize a :class:`Detector2D_SXRD` using its central NX serializer."""
    return detector.toNXdict()["detector_SXRD"]


def detector_from_nxdict(nxdict):
    """Load a :class:`Detector2D_SXRD` through the central NX reader."""
    return DetectorCalibration.loadNXdict(nxdict)


def unit_cell_to_nxdict(unit_cell):
    """Serialize a :class:`UnitCell` without fit parameters.

    Lattice lengths are Angstrom and lattice angles are degrees, matching
    :class:`UnitCell` construction.
    """
    return {
        "@NX_class": "NXsample",
        "unit_cell_abc": np.asarray(unit_cell.a, dtype=np.float64),
        "unit_cell_alphabetagamma": np.rad2deg(
            np.asarray(unit_cell.alpha, dtype=np.float64)
        ),
        "unit_cell_volume": float(unit_cell.volume),
        "orgui_unit_cell": {
            "@NX_class": "NXcollection",
            "class_name": "UnitCell",
            "name": unit_cell.name,
            "atom_names": _string_array(unit_cell.names),
            "basis": np.asarray(unit_cell.basis, dtype=np.float64),
            "fractional_coordinates": np.asarray(
                unit_cell.basis[:, 1:4], dtype=np.float64
            )
            if unit_cell.basis.size
            else np.empty((0, 3), dtype=np.float64),
            "debye_waller": np.asarray(unit_cell.basis[:, 4:6], dtype=np.float64)
            if unit_cell.basis.size
            else np.empty((0, 2), dtype=np.float64),
            "occupancies": np.asarray(unit_cell.basis[:, 6], dtype=np.float64)
            if unit_cell.basis.size
            else np.empty((0,), dtype=np.float64),
            "layers": np.asarray(unit_cell.basis[:, 7], dtype=np.float64)
            if unit_cell.basis.size and unit_cell.basis.shape[1] > 7
            else np.empty((0,), dtype=np.float64),
            "coherent_domain_matrices": np.asarray(
                unit_cell.coherentDomainMatrix, dtype=np.float64
            ),
            "coherent_domain_occupancies": np.asarray(
                unit_cell.coherentDomainOccupancy, dtype=np.float64
            ),
        },
    }


def unit_cell_from_nxdict(nxdict):
    """Load a :class:`UnitCell` from a config sample dictionary."""
    sample = nxdict
    if "sample" in sample:
        sample = sample["sample"]
    ucdata = sample.get("orgui_unit_cell", {})
    a = np.asarray(sample["unit_cell_abc"], dtype=np.float64)
    alpha = np.asarray(sample["unit_cell_alphabetagamma"], dtype=np.float64)
    name = _as_text(ucdata.get("name", "unnamed"))
    unit_cell = CTRcalc.UnitCell(a, alpha, name=name)
    names = _read_string_array(ucdata.get("atom_names", []))
    basis = np.asarray(ucdata.get("basis", np.empty((0, 8))), dtype=np.float64)
    for atom_name, row in zip(names, basis):
        unit_cell.addAtom(
            atom_name,
            row[1:4],
            _none_to_nan(row[4]),
            _none_to_nan(row[5]),
            _none_to_nan(row[6]),
            row[7] if row.size > 7 else 0,
        )
    if "coherent_domain_matrices" in ucdata:
        unit_cell.coherentDomainMatrix = [
            np.asarray(matrix, dtype=np.float64)
            for matrix in ucdata["coherent_domain_matrices"]
        ]
    if "coherent_domain_occupancies" in ucdata:
        unit_cell.coherentDomainOccupancy = [
            float(occupancy) for occupancy in ucdata["coherent_domain_occupancies"]
        ]
    return unit_cell


def reflections_to_nxdict(reflections, ub_calculator=None):
    """Serialize reference reflections.

    Reflection HKL values are in r.l.u., detector coordinates are pixels, and
    six-circle angle snapshots are stored in degrees as audit metadata.
    """
    hkls = []
    xy = []
    imageno = []
    identifiers = []
    angles = []
    for idx, reflection in enumerate(reflections or []):
        hkls.append(np.asarray(reflection.hkl, dtype=np.float64))
        xy.append(np.asarray(reflection.xy, dtype=np.float64))
        imageno.append(int(reflection.imageno))
        identifiers.append(getattr(reflection, "identifier", f"ref_{idx}"))
        angle = np.full((6,), np.nan)
        angles.append(np.rad2deg(np.asarray(angle, dtype=np.float64)))
    return {
        "@NX_class": "NXdata",
        "hkl": np.asarray(hkls, dtype=np.float64).reshape((-1, 3)),
        "xy": np.asarray(xy, dtype=np.float64).reshape((-1, 2)),
        "image_number": np.asarray(imageno, dtype=np.int64),
        "identifier": _string_array(identifiers),
        "sixc_angles": {
            "@NX_class": "NXpositioner",
            "@unit": "deg",
            "angles": np.asarray(angles, dtype=np.float64).reshape((-1, 6)),
            "axis_names": _string_array(
                ["alpha", "delta", "gamma", "omega", "chi", "phi"]
            ),
        },
    }


def reflections_from_nxdict(nxdict):
    """Load reference reflections from a config dictionary."""
    from .QReflectionSelector import HKLReflection

    refl = nxdict.get("reference_reflections", nxdict)
    hkls = np.asarray(refl.get("hkl", np.empty((0, 3))), dtype=np.float64).reshape(
        (-1, 3)
    )
    xy = np.asarray(refl.get("xy", np.empty((0, 2))), dtype=np.float64).reshape((-1, 2))
    imageno = np.asarray(refl.get("image_number", []), dtype=np.int64)
    identifiers = _read_string_array(refl.get("identifier", []))
    return [
        HKLReflection(point, hkl, int(img), ident)
        for point, hkl, img, ident in zip(xy, hkls, imageno, identifiers)
    ]


@dataclass
class ConfigData:
    """Physical application state persisted with scans and integrations."""

    detector: DetectorCalibration.Detector2D_SXRD
    unit_cell: CTRcalc.UnitCell
    ub_calculator: HKLVlieg.UBCalculator
    mu: float = 0.0
    chi: float = 0.0
    phi: float = 0.0
    refraction_index: float = 1.0
    reference_reflections: list = field(default_factory=list)
    orgui: dict = field(default_factory=dict)

    @classmethod
    def from_ini(cls, filename):
        """Read legacy INI config into a structured config object.

        Config ``SDD`` and ``pixelsize`` are read in meters, detector center is
        read in pixels, energy is keV, and lattice constants are Angstrom.
        """
        config = configparser.ConfigParser()
        config.read(filename)
        machine = config["Machine"]
        lattice = config["Lattice"]
        diffrac = config["Diffractometer"]

        energy = machine.getfloat("E", 78.0)
        cell = CTRcalc.UnitCell(
            [
                lattice.getfloat("a1", 1.0),
                lattice.getfloat("a2", 1.0),
                lattice.getfloat("a3", 1.0),
            ],
            [
                lattice.getfloat("alpha1", 90.0),
                lattice.getfloat("alpha2", 90.0),
                lattice.getfloat("alpha3", 90.0),
            ],
        )
        cell.addAtom("Pt", [0.0, 0.0, 0.0], 0.1, 0.1, 1.0)
        cell.setEnergy(energy * 1e3)
        ub_calculator = HKLVlieg.UBCalculator(cell, energy)
        ub_calculator.defaultU()

        detector = DetectorCalibration.Detector2D_SXRD()
        detector.setFit2D(
            machine.getfloat("SDD", 0.729) * 1e3,
            machine.getfloat("cpx", 731),
            machine.getfloat("cpy", 1587),
            pixelX=machine.getfloat("pixelsize", 172e-6) * 1e6,
            pixelY=machine.getfloat("pixelsize", 172e-6) * 1e6,
        )
        detector.wavelength = ub_calculator.getLambda() * 1e-10
        detector.detector.shape = (
            int(machine.getfloat("sizey", 3000)),
            int(machine.getfloat("sizex", 3000)),
        )
        detector.detector.max_shape = detector.detector.shape
        detector.setAzimuthalReference(
            np.deg2rad(diffrac.getfloat("azimuthal_reference", 0))
        )
        detector.setPolarization(
            np.deg2rad(diffrac.getfloat("polarization_axis", 0)),
            diffrac.getfloat("polarization_factor", 0),
        )
        return cls(
            detector=detector,
            unit_cell=cell,
            ub_calculator=ub_calculator,
            mu=np.deg2rad(diffrac.getfloat("mu", 0.05)),
            chi=np.deg2rad(diffrac.getfloat("chi", 0.0)),
            phi=np.deg2rad(diffrac.getfloat("phi", 0.0)),
            refraction_index=1.0 - lattice.getfloat("refractionindex", 0.0),
        )

    readConfig = from_ini

    @classmethod
    def from_gui(cls, gui):
        """Capture config state from an orGUI main window or QUBCalculator."""
        ub_widget = getattr(gui, "ubcalc", gui)
        reflections = []
        if hasattr(gui, "reflectionSel"):
            reflections = list(gui.reflectionSel.reflections)
        cell = unit_cell_from_nxdict(unit_cell_to_nxdict(ub_widget.crystal))
        ub_calculator = HKLVlieg.UBCalculator(cell, ub_widget.ubCal.getEnergy())
        ub_calculator.setU(ub_widget.ubCal.getU())
        return cls(
            detector=ub_widget.detectorCal,
            unit_cell=cell,
            ub_calculator=ub_calculator,
            mu=getattr(ub_widget, "mu", 0.0),
            chi=getattr(ub_widget, "chi", 0.0),
            phi=getattr(ub_widget, "phi", 0.0),
            refraction_index=getattr(ub_widget, "n", 1.0),
            reference_reflections=reflections,
        )

    def apply_to_gui(self, gui):
        """Apply this config to an orGUI main window or QUBCalculator."""
        ub_widget = getattr(gui, "ubcalc", gui)
        ub_widget.detectorCal = detector_from_nxdict(detector_to_nxdict(self.detector))
        ub_widget.crystal = unit_cell_from_nxdict(unit_cell_to_nxdict(self.unit_cell))
        ub_widget.ubCal = HKLVlieg.UBCalculator(
            ub_widget.crystal, self.ub_calculator.getEnergy()
        )
        ub_widget.ubCal.setU(self.ub_calculator.getU())
        ub_widget.mu = self.mu
        ub_widget.chi = self.chi
        ub_widget.phi = self.phi
        ub_widget.n = self.refraction_index
        if hasattr(ub_widget, "uedit"):
            ub_widget.uedit.setU(ub_widget.ubCal.getU())
        if hasattr(ub_widget, "crystalparams"):
            ub_widget.crystalparams.setValues(ub_widget.crystal, ub_widget.n)
        if hasattr(ub_widget, "machineParams"):
            ub_widget.machineParams.setValues(
                {
                    "diffractometer": {
                        "mu": self.mu,
                        "phi": self.phi,
                        "chi": self.chi,
                    },
                    "source": {"E": ub_widget.detectorCal.energy},
                    "SXRD_geometry": ub_widget.detectorCal,
                }
            )
        if hasattr(gui, "reflectionSel"):
            gui.reflectionSel.setReflections(self.reference_reflections)
        if hasattr(ub_widget, "updateReflectionMismatch"):
            ub_widget.updateReflectionMismatch()

    def to_nxdict(self, role="scan", source=None):
        """Convert this config to a NeXus-compatible nested dictionary."""
        nxdict = {
            "@NX_class": "NXcollection",
            "@orgui_meta": "config",
            "@orgui_config_role": role,
            "@orgui_schema_version": SCHEMA_VERSION,
            "@orgui_config_created": datetime.datetime.now(
                datetime.UTC
            ).isoformat(),
            "instrument": {
                "@NX_class": "NXinstrument",
                "detector_SXRD": detector_to_nxdict(self.detector),
            },
            "sample": unit_cell_to_nxdict(self.unit_cell),
            "reference_reflections": reflections_to_nxdict(
                self.reference_reflections, self.ub_calculator
            ),
            "orgui": {
                "@NX_class": "NXcollection",
                "diffractometer": {
                    "@NX_class": "NXcollection",
                    "mu": self.mu,
                    "chi": self.chi,
                    "phi": self.phi,
                    "@unit": "rad",
                },
                "source": {
                    "@NX_class": "NXcollection",
                    "energy": self.ub_calculator.getEnergy(),
                    "@energy_unit": "keV",
                    "wavelength": self.ub_calculator.getLambda(),
                    "@wavelength_unit": "Angstrom",
                },
                "refraction_index": self.refraction_index,
                **self.orgui,
            },
        }
        nxdict["sample"]["orientation_matrix"] = self.ub_calculator.getU()
        nxdict["sample"]["ub_matrix"] = self.ub_calculator.getUB()
        if source is not None:
            nxdict["@orgui_config_source"] = source
        return nxdict

    @classmethod
    def from_nxdict(cls, nxdict):
        """Load config state from a NeXus-compatible dictionary."""
        detector = detector_from_nxdict(nxdict["instrument"]["detector_SXRD"])
        unit_cell = unit_cell_from_nxdict(nxdict["sample"])
        source = nxdict.get("orgui", {}).get("source", {})
        energy = float(source.get("energy", detector.energy))
        ub_calculator = HKLVlieg.UBCalculator(unit_cell, energy)
        ub_calculator.setU(np.asarray(nxdict["sample"]["orientation_matrix"]))
        diffrac = nxdict.get("orgui", {}).get("diffractometer", {})
        return cls(
            detector=detector,
            unit_cell=unit_cell,
            ub_calculator=ub_calculator,
            mu=float(diffrac.get("mu", 0.0)),
            chi=float(diffrac.get("chi", 0.0)),
            phi=float(diffrac.get("phi", 0.0)),
            refraction_index=float(
                nxdict.get("orgui", {}).get("refraction_index", 1.0)
            ),
            reference_reflections=reflections_from_nxdict(nxdict),
        )


class ConfigHandler:
    """Read and write orGUI config groups in an HDF5 database."""

    def __init__(self, gui=None, create_dataset_args=None):
        self.gui = gui
        self.create_dataset_args = create_dataset_args or {}

    @staticmethod
    def is_config_group(group):
        """Return ``True`` for loadable orGUI config HDF5 groups."""
        return (
            isinstance(group, h5py.Group)
            and group.attrs.get("orgui_meta") == "config"
        )

    def write_scan_config(self, scan_group, config, source="scan_import"):
        """Write ``configuration`` under a scan group."""
        return self._write_config(scan_group, config, "scan", source)

    def write_integration_config(
        self, integration_group, config, source="integration_save"
    ):
        """Write ``configuration`` under an integration result group."""
        return self._write_config(integration_group, config, "integration", source)

    def _write_config(self, group, config, role, source):
        if "configuration" in group:
            del group["configuration"]
        dicttonx(
            {"configuration": config.to_nxdict(role=role, source=source)},
            group,
            update_mode="add",
            create_dataset_args=self.create_dataset_args,
        )
        return group["configuration"]

    def load_config_group(self, group):
        """Read a marked config group and return :class:`ConfigData`."""
        if not self.is_config_group(group):
            raise ValueError(f"Not an orGUI config group: {group.name}")
        return ConfigData.from_nxdict(nxtodict(group))

    def apply_config_group(self, group, gui=None):
        """Load a marked config group and apply it to a GUI object."""
        config = self.load_config_group(group)
        target = gui or self.gui
        if target is None:
            raise ValueError("No GUI target is configured for loading config groups.")
        config.apply_to_gui(target)
        return config

# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
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
"""Example ID31 EBS backend for Pilatus 4 / 4M HDF5 scans.

This file is meant to be loaded through an orGUI config ``[backend] file = ...``
entry or through the "Load backend file" GUI action. It intentionally contains
exactly one :class:`orgui.backend.scans.Scan` subclass.
"""

import copy
import os
import traceback
import warnings

import numpy as np
import scipy.interpolate
from silx.io import dictdump
import silx.io.h5py_utils

from orgui.backend.scans import Scan


class ID31DiffractLinTilt:
    """Convert the ID31 ``linai`` linear translation to incident angle ``mu``.

    :param float muoffset:
        Experiment-specific alignment offset in deg.
    """

    def __init__(self, muoffset=-0.06442416994811659):
        self.config = {
            "a": 938,
            "b": 400,
            "muoffset": muoffset,
        }

    @property
    def a(self):
        """Return linkage length ``a`` in mm."""
        return self.config.get("a", float)

    @property
    def b(self):
        """Return linkage length ``b`` in mm."""
        return self.config.get("b", float)

    @property
    def muoffset(self):
        """Return the experiment-specific ``mu`` offset in deg."""
        return self.config.get("muoffset", float)

    def c(self, a=None, b=None):
        """Return the derived linkage length in mm."""
        a = self.a if a is None else a
        b = self.b if b is None else b
        return np.sqrt(np.square(a) + np.square(b))

    def d(self, a=None, b=None):
        """Return the derived linkage angle in rad."""
        a = self.a if a is None else a
        b = self.b if b is None else b
        return np.arctan(b / a)

    def calc_from_real(self, positions_dict):
        """Convert linear translation in mm to tilt angle in deg.

        :param dict positions_dict:
            Must contain ``{"linear": ...}``.
        :returns:
            ``{"tilt": ...}`` with tilt in deg.
        :rtype: dict
        """
        a, b = self.a, self.b
        c, d = self.c(a, b), self.d(a, b)
        a2, c2 = np.square(a), np.square(c)
        linear = positions_dict["linear"]
        bc2 = np.square(b + linear)
        tilt = np.arccos((a2 + c2 - bc2) / (2 * a * c)) - d
        return {"tilt": np.rad2deg(tilt)}

    def calc_to_real(self, positions_dict):
        """Convert tilt angle in deg to linear translation in mm.

        :param dict positions_dict:
            Must contain ``{"tilt": ...}``.
        :returns:
            ``{"linear": ...}`` with translation in mm.
        :rtype: dict
        """
        a, b = self.a, self.b
        c, d = self.c(a, b), self.d(a, b)
        a2, c2 = np.square(a), np.square(c)
        tilt = np.deg2rad(positions_dict["tilt"])
        bc = np.sqrt(a2 + c2 - 2 * a * c * np.cos(tilt + d))
        return {"linear": bc - b}

    def linai_to_mu(self, linai):
        """Convert ``linai`` in mm to incident angle ``mu`` in deg."""
        tilt = self.calc_from_real({"linear": linai})
        return tilt["tilt"] - self.muoffset

    def mu_to_linai(self, mu):
        """Convert incident angle ``mu`` in deg to ``linai`` in mm."""
        tilts = {"tilt": mu + self.muoffset}
        return self.calc_to_real(tilts)["linear"]


class h5_Image:
    """Minimal image container returned by the backend scan API."""

    def __init__(self, data):
        """Store detector image data as a NumPy-like array."""
        self.img = data
        self.motors = {}
        self.counters = {}


class ID31_EBS_p4_2026(Scan):
    """Load ID31 EBS Pilatus 4 / 4M scans from BLISS HDF5 files.

    :param hdffilepath_orNode:
        Either a path to an HDF5 file or an opened silx HDF5 node.
    :param int scanno:
        Scan number to load.
    :param bool loadimg:
        If true, load all detector frames into memory. If false, images are
        read lazily from HDF5 in :meth:`get_raw_img`.
    :param bool saveh5data:
        If true, keep the parsed HDF5 dictionaries on the scan object.
    :param float muoffset:
        Experiment-specific ``linai`` to ``mu`` offset in deg. The default
        value ``-2.31`` matches the CH8153 branch configuration.
    """

    def __init__(
        self,
        hdffilepath_orNode=None,
        scanno=None,
        loadimg=True,
        saveh5data=False,
        muoffset=-2.31,
    ):
        self.cameras = []
        self.loadimg = loadimg
        self.muoffset = muoffset
        data_1 = None
        data_2 = None
        excludenames = None if self.loadimg else ["p4_lima1"]

        if hdffilepath_orNode is None:
            return
        self.hdffilepath_orNode = hdffilepath_orNode

        if isinstance(hdffilepath_orNode, str):
            filepath, _filename = os.path.split(hdffilepath_orNode)
            _unused, filename_noext = os.path.split(filepath)
            self.filename_base = filename_noext
            with silx.io.h5py_utils.File(hdffilepath_orNode) as h5file:
                data_1, data_2 = self._read_scan_groups(
                    h5file, scanno, excludenames
                )
        else:
            hdffilepath = hdffilepath_orNode.local_filename
            filepath, _filename = os.path.split(hdffilepath)
            _unused, filename_noext = os.path.split(filepath)
            self.filename_base = filename_noext
            data_1, data_2 = self._read_scan_groups(
                hdffilepath_orNode.file, scanno, excludenames
            )

        if saveh5data:
            self.data_1 = data_1
            self.data_2 = data_2

        if "title" in data_1:
            self.title = data_1["title"]

        self.positioners = data_1["instrument"]["positioners"]
        self._read_detector_data(data_1)
        self._read_scan_axis(data_1)
        self._read_metadata_counters(data_1, data_2)

        self.offsetindex = 0
        self.scandatapoints = self.nopoints
        self.imageno = np.arange(self.nopoints)
        self.omega = -1 * self.th
        self.axis = getattr(self, self.axisname)

    def _read_scan_groups(self, h5file, scanno, excludenames):
        for group_name in h5file:
            scansuffix = group_name.split("_")[-1]
            scanname_nosuffix = "_".join(group_name.split("_")[:-1])
            scanno_s, _subscanno = scansuffix.split(".")
            if int(scanno_s) == scanno:
                break
        else:
            raise OSError(f"Scan number {scanno} not found in file")

        self.scanno1 = str(scanno) + "." + "1"
        self.scanno2 = str(scanno) + "." + "2"
        self.scanname_1 = scanname_nosuffix + "_" + self.scanno1
        self.scanname_2 = scanname_nosuffix + "_" + self.scanno2
        self.scanname = self.scanname_1
        self.name = self.scanname

        if self.scanname_1 in h5file:
            data_1 = dictdump.h5todict(
                h5file, self.scanname_1, exclude_names=excludenames
            )
            if self.scanname_2 in h5file:
                data_2 = dictdump.h5todict(
                    h5file, self.scanname_2, exclude_names=excludenames
                )
            else:
                data_2 = None
        else:
            data_1 = dictdump.h5todict(
                h5file, self.scanno1, exclude_names=excludenames
            )
            if self.scanno2 in h5file:
                data_2 = dictdump.h5todict(
                    h5file, self.scanno2, exclude_names=excludenames
                )
            else:
                data_2 = None

        measurement = data_1["measurement"]
        if "p4_lima1" in measurement:
            self.cameras.append("p4_lima1")
            self.nopoints = measurement["p4_lima1"].shape[0]
        if "mpx" in measurement:
            self.cameras.append("mpx")
            self.nopoints = measurement["mpx"].shape[0]

        return data_1, data_2

    def _read_detector_data(self, data_1):
        measurement = data_1["measurement"]
        if "p4_lima1" in measurement:
            if "p4_lima1" not in self.cameras:
                self.cameras.append("p4_lima1")
            self.nopoints = measurement["p4_lima1"].shape[0]
            if self.loadimg:
                self.p4 = measurement["p4_lima1"][()]

        if "mpx" in measurement:
            self.mpx = measurement["mpx"][()]
            self.nopoints = measurement["mpx"].shape[0]
            if "mpx" not in self.cameras:
                self.cameras.append("mpx")

    def _read_scan_axis(self, data_1):
        measurement = data_1["measurement"]
        if "th" in measurement:
            self.th = measurement["th"][: self.nopoints]
            if "th_trig" in measurement and "th_delta" in measurement:
                self.th = (
                    measurement["th_trig"][: self.nopoints]
                    + measurement["th_delta"][: self.nopoints] / 2
                )
            self.axisname = "th"
            self.mu = self.positioners["mu"]
        elif "nth" in measurement:
            self.th = measurement["nth"][: self.nopoints]
            self.axisname = "th"
            self.mu = self.positioners["nai"] * -1
        elif "uth" in measurement:
            self.th = measurement["uth"][: self.nopoints]
            self.axisname = "th"
            self.mu = self.positioners["mu"] * -1
        elif "mu" in measurement:
            self.mu = measurement["mu"][: self.nopoints]
            self.axisname = "mu"
            self.th = self.positioners["th"]
        elif "nai" in measurement:
            self.mu = measurement["nai"][: self.nopoints]
            self.axisname = "mu"
            self.th = self.positioners["nth"]
        elif "linai" in measurement:
            lintomu = ID31DiffractLinTilt(muoffset=self.muoffset)
            self.linai = measurement["linai"][: self.nopoints]
            if "linai_trig" in measurement and "linai_delta" in measurement:
                self.linai = (
                    measurement["linai_trig"][: self.nopoints]
                    + measurement["linai_delta"][: self.nopoints] / 2
                )
            self.mu = lintomu.linai_to_mu(self.linai)
            self.axisname = "mu"
            self.th = self.positioners["th"]
        else:
            self.axisname = "time"
            self.th = self.positioners["th"]
            self.mu = self.positioners["mu"]

    def _read_metadata_counters(self, data_1, data_2):
        measurement = data_1["measurement"]
        if "potv" in measurement:
            self.potv = measurement["potv"][: self.nopoints]
        if "scaled_potv2f" in measurement:
            self.scaled_potv2f = measurement["scaled_potv2f"][: self.nopoints]

        if "srcur" in measurement:
            self.srcur = measurement["srcur"][: self.nopoints]
        else:
            self.srcur = np.ones(self.nopoints)

        if "timer_trig" in measurement:
            self.time = measurement["timer_trig"][: self.nopoints]
        elif "elapsed_time" in measurement:
            self.time = measurement["elapsed_time"][: self.nopoints]
        else:
            raise OSError("Cannot find time counter in scan")

        if "epoch_trig" in measurement:
            self.epoch = measurement["epoch_trig"][: self.nopoints]
        elif "epoch" in measurement:
            self.epoch = measurement["epoch"][: self.nopoints]
        else:
            raise OSError("Cannot find epoch counter in scan")

        self._read_slow_counters(data_1, data_2)

        if "mondio" in measurement:
            self.mondio = measurement["mondio"][: self.nopoints]
        if not hasattr(self, "mondio"):
            self.mondio = np.ones(self.nopoints)

        self.mondio = self.mondio / np.mean(self.mondio)
        self.srcur = self.srcur / np.mean(self.srcur)

        if "timer_delta" in measurement:
            self.exposure_time = measurement["timer_delta"][: self.nopoints]
        elif "sec" in measurement:
            self.exposure_time = measurement["sec"][: self.nopoints]
        else:
            raise OSError("Cannot find exposure time in scan")
        self.relative_exposure = self.exposure_time / np.mean(self.exposure_time)

    def _read_slow_counters(self, data_1, data_2):
        if data_2 is not None:
            for counter in ("mondio", "potential", "current"):
                if counter in data_2["measurement"]:
                    try:
                        interpolator = scipy.interpolate.interp1d(
                            data_2["measurement"]["epoch"],
                            data_2["measurement"][counter],
                        )
                        setattr(self, counter, interpolator(self.epoch))
                    except ValueError:
                        warnings.warn(
                            "Cannot interpolate "
                            f"{counter} in scan {self.scanno2}\n"
                            f"{traceback.format_exc()}"
                        )
                        setattr(self, counter, data_2["measurement"][counter])
        else:
            measurement = data_1["measurement"]
            if "potential" in measurement:
                self.potential = measurement["potential"][: self.nopoints]
            if "current" in measurement:
                self.current = measurement["current"][: self.nopoints]

    @classmethod
    def parse_h5_node(cls, obj):
        """Parse the scan number from an ID31 BLISS HDF5 node.

        :param obj:
            silx HDF5 tree object selected by the GUI.
        :returns:
            Dictionary with ``scanno`` and ``name`` keys.
        :rtype: dict
        """
        scanname = obj.local_name
        if "_" in scanname:
            scansuffix = scanname.split("_")[-1]
        elif "/" in scanname:
            scansuffix = scanname.split("/")[-1]
        else:
            scansuffix = scanname
        scanno, _subscanno = scansuffix.split(".")
        return {"scanno": int(scanno), "name": obj.local_name}

    @property
    def auxillary_counters(self):
        """Return counters copied to the orGUI integration database."""
        return [
            "current",
            "potential",
            "exposure_time",
            "elapsed_time",
            "time",
            "srcur",
            "mondio",
            "epoch",
            "scaled_potv2f",
        ]

    def get_raw_img(self, img):
        """Return detector image ``img`` as an :class:`h5_Image`.

        :param int img:
            Zero-based image index.
        """
        if not hasattr(self, "p4"):
            if isinstance(self.hdffilepath_orNode, str):
                with silx.io.h5py_utils.File(self.hdffilepath_orNode, "r") as h5file:
                    if self.scanname_1 in h5file:
                        data_1 = h5file[self.scanname_1]
                    else:
                        data_1 = h5file[self.scanno1]
                    image = data_1["measurement"]["p4_lima1"][img][()]
            else:
                h5file = self.hdffilepath_orNode.file
                if self.scanname_1 in h5file:
                    data_1 = h5file[self.scanname_1]
                else:
                    data_1 = h5file[self.scanno1]
                image = data_1["measurement"]["p4_lima1"][img][()]
        else:
            image = self.p4[img]
        return h5_Image(image)

    def slice(self, startno, endno):
        """Return a scan object containing images ``startno`` through ``endno``.

        The returned object preserves the original scan metadata and slices
        NumPy arrays along their first dimension.
        """
        if startno > endno:
            startno, endno = endno, startno

        fscan = ID31_EBS_p4_2026()
        for key, value in self.__dict__.items():
            if isinstance(value, np.ndarray) and value.ndim > 0:
                fscan.__dict__[key] = copy.deepcopy(value[startno:endno])
            else:
                fscan.__dict__[key] = copy.deepcopy(value)

        if self.loadimg is False:
            with silx.io.h5py_utils.File(self.hdffilepath_orNode, "r") as h5file:
                if self.scanname_1 in h5file:
                    data_1 = h5file[self.scanname_1]
                else:
                    data_1 = h5file[self.scanno1]
                fscan.p4 = data_1["measurement"]["p4_lima1"][startno:endno][()]

        fscan.offsetindex = copy.deepcopy(startno)
        fscan.nopoints = copy.deepcopy(endno - startno)
        return fscan

    def get_p3_img(self, img):
        """Return an image with ID31 counters attached for normalization."""
        imgdata = self.get_raw_img(img)
        imgdata.counters["Time"] = self.time[img]
        imgdata.counters["exposure"] = self.exposure_time[img]
        imgdata.counters["TrigTime"] = self.relative_exposure[img]
        imgdata.counters["imageno"] = self.imageno[img]

        if self.axisname == "th":
            imgdata.counters["th"] = self.th[img]
            imgdata.counters["om"] = self.omega[img]
        else:
            imgdata.counters["th"] = self.th
            imgdata.counters["om"] = self.omega

        if self.axisname == "mu":
            imgdata.counters["mu"] = self.mu[img]
        else:
            imgdata.counters["mu"] = self.mu

        if self.srcur is not None:
            imgdata.counters["srcur"] = self.srcur[img]
        if self.mondio is not None:
            imgdata.counters["mondio"] = self.mondio[img]
        return imgdata

    def __getitem__(self, key):
        """Return image ``key`` with ID31 counters attached."""
        return self.get_p3_img(key)

    def __len__(self):
        """Return the number of detector images in the scan."""
        return self.nopoints

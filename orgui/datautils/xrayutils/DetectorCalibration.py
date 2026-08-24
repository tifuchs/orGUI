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
import pyFAI
from pyFAI import geometry
from pyFAI.utils.mathutil import binning
from scipy.spatial.transform import Rotation

# import pyFAI.azimuthalIntegrator
import copy
import warnings

from .HKLVlieg import (
    calcDELTA,
    calcGAMMA,
    crystalAngles_singleArray,
    primBeamAngles,
    vacAngles_singleArray,
)


def load(ponifile):
    det = Detector2D_SXRD()
    det.load(ponifile)
    return det


def loadNXdict(nxdict):
    det = Detector2D_SXRD()
    det.fromNXdict(nxdict)
    return det


def _nx_to_python(value):
    if isinstance(value, dict):
        return {key: _nx_to_python(item) for key, item in value.items()}
    if isinstance(value, np.ndarray) and value.shape == ():
        return _nx_to_python(value.item())
    if isinstance(value, bytes):
        return value.decode()
    return value


"""

Attention!!!!

There is something wierd with pyFAI

x,y seem to be inverted

"""


def azimuthRelabel(deltaChi):
    """Matrix taking a pyFAI laboratory vector into the Vlieg laboratory frame.

    pyFAI orders its laboratory components ``(t1, t2, t3)`` -- slow detector
    dimension, fast detector dimension, along the beam -- while the Vlieg
    diffractometer equations in :mod:`~orgui.datautils.xrayutils.HKLVlieg` use
    ``(x, y, z)`` with ``y`` along the beam and ``z`` along the surface normal.
    The azimuthal reference fixes the roll between the two.

    :param float deltaChi:
        Azimuthal reference of the detector calibration, in rad, as returned by
        :meth:`Detector2D_SXRD.getAzimuthalReference`.
    :returns: 3x3 matrix mapping ``(t1, t2, t3)`` onto ``(x, y, z)``.
    :rtype: numpy.ndarray
    """
    cosine = np.cos(deltaChi)
    sine = np.sin(deltaChi)
    return np.array(
        [
            [cosine, sine, 0.0],
            [0.0, 0.0, 1.0],
            [-sine, cosine, 0.0],
        ]
    )


def _primBeamFromTthAzimuth(tth, azimuth):
    """Laboratory-frame ``(gamma_p, delta_p)`` from pyFAI scattering angles.

    :param tth: scattering angle in rad.
    :param azimuth: pyFAI azimuth already offset by the azimuthal reference,
        in rad.
    :returns: ``(gamma_p, delta_p)`` in rad.
    :rtype: tuple
    """
    sintth = np.sin(tth)
    costth = np.cos(tth)
    delta_p = np.arctan2(sintth * np.sin(azimuth), costth)
    gamma_p = np.arctan2(sintth * np.cos(azimuth) * np.cos(delta_p), costth)
    return gamma_p, delta_p


def pyFAIRotationMatrix(rot1, rot2, rot3):
    """Detector rotation matrix of a pyFAI geometry.

    Reproduces :meth:`pyFAI.geometry.Geometry.rotation_matrix` as a free
    function, so that a geometry can be rotated without an instance to hang it
    on. pyFAI builds ``rot3 @ rot2 @ rot1``, which is the intrinsic Z-Y-X Euler
    sequence ``Rz(rot3) Ry(-rot2) Rx(-rot1)``.

    :param float rot1, rot2, rot3: pyFAI detector rotations in rad.
    :returns: 3x3 rotation matrix in pyFAI's laboratory frame.
    :rtype: numpy.ndarray
    """
    return Rotation.from_euler("ZYX", [rot3, -rot2, -rot1]).as_matrix()


def pyFAIRotationAngles(matrix, tolerance=1e-9):
    """Decompose a pyFAI detector rotation matrix into ``(rot1, rot2, rot3)``.

    Inverse of :func:`pyFAIRotationMatrix`.

    :param numpy.ndarray matrix: 3x3 rotation matrix in pyFAI's laboratory frame.
    :param float tolerance:
        How close ``|rot2|`` may come to ``pi/2`` before the decomposition is
        rejected, in rad.
    :returns: ``(rot1, rot2, rot3)`` in rad.
    :rtype: tuple[float, float, float]
    :raises ValueError:
        At the gimbal lock of this Euler sequence, ``|rot2| = pi/2``, where
        ``rot1`` and ``rot3`` are no longer separable. Raised rather than
        silently returning the degenerate triple SciPy would pick.
    """
    with warnings.catch_warnings():
        # SciPy warns and zeroes the third angle at gimbal lock; we would
        # rather report it, so silence the warning and test for it below.
        warnings.simplefilter("ignore", UserWarning)
        rot3, minus_rot2, minus_rot1 = Rotation.from_matrix(
            np.asarray(matrix, dtype=float)
        ).as_euler("ZYX")
    rot1, rot2 = -minus_rot1, -minus_rot2
    if abs(abs(rot2) - np.pi / 2) < tolerance:
        raise ValueError(
            "Detector orientation is at the gimbal lock of the pyFAI rotation "
            f"parameterization (rot2 = {rot2:.6f} rad, i.e. +/- pi/2), where "
            "rot1 and rot3 are degenerate. For a detector arm this means an "
            "out-of-plane angle of +/- 90 degrees."
        )
    return float(rot1), float(rot2), float(rot3)


def _validatedRotation(matrix, name):
    """Check that ``matrix`` is a proper rotation, and return it as float64."""
    matrix = np.ascontiguousarray(matrix, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must be a 3x3 matrix, got shape {matrix.shape}")
    if not np.allclose(matrix @ matrix.T, np.identity(3), atol=1e-8):
        raise ValueError(f"{name} must be orthogonal")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1e-8):
        raise ValueError(
            f"{name} must be a proper rotation with determinant +1, "
            f"got {np.linalg.det(matrix):.6f}"
        )
    return matrix


def armAnglesFromSurface(gamma_arm, delta_arm, alpha_i):
    """Convert six-circle arm angles into true scattering angles.

    The arm API works in the primary-beam frame, i.e. in ``gamma_p`` and
    ``delta_p`` (see :func:`armRotation`). A diffractometer whose arm motors
    are read out in the six-circle surface frame -- ``gamma``, ``delta``,
    the angles the Vlieg equations use -- converts them here once, at the
    point where the motor values enter, rather than everywhere downstream.

    Reuses :func:`~orgui.datautils.xrayutils.HKLVlieg.primBeamAngles`, the
    same surface-to-primary-beam conversion the per-pixel code path uses. The
    arm has only two degrees of freedom, so its orientation is fully fixed by
    the direction of its central ray and nothing is lost in the conversion.

    :param float gamma_arm: out-of-plane arm angle in the surface frame, rad.
    :param float delta_arm: in-plane arm angle in the surface frame, rad.
    :param float alpha_i: incidence angle in rad.
    :returns: ``(gamma_p_arm, delta_p_arm)`` in rad.
    :rtype: tuple[float, float]

    This assumes the arm is mounted independently of the incidence circle, so
    that it rotates in the laboratory frame. For an arm carried *by* the alpha
    circle the laboratory rotation is ``ALPHA @ R_arm @ ALPHA^T`` instead,
    which differs from this by a roll about the outgoing beam of up to a
    degree at grazing incidence and much more beyond it.
    """
    _, delta_p, gamma_p = primBeamAngles(
        [alpha_i, delta_arm, gamma_arm, 0.0, 0.0, 0.0]
    )[:3]
    return gamma_p, delta_p


def armRotation(gamma_arm=0.0, delta_arm=0.0):
    """Rotation matrix of a two-circle detector arm, in the Vlieg frame.

    The angles are the **true scattering angles** of the arm's central ray,
    measured from the primary beam: ``gamma_arm`` is ``gamma_p`` and
    ``delta_arm`` is ``delta_p``. They do not depend on the incidence angle,
    which is what a floor-mounted arm reads. Six-circle surface-frame motor
    readouts have to be converted first, with :func:`armAnglesFromSurface`.

    The composition ``DELTA(delta) @ GAMMA(gamma)`` puts ``delta`` on the base
    circle about the vertical axis with the ``gamma`` circle riding on the
    delta arm. A diffractometer that stacks its circles the other way round
    should build its own matrix and pass it in directly; the two orderings
    differ by up to tens of degrees of detector roll.

    :param float gamma_arm: out-of-plane arm angle ``gamma_p`` in rad.
    :param float delta_arm: in-plane arm angle ``delta_p`` in rad.
    :returns: 3x3 rotation matrix in the Vlieg laboratory frame.
    :rtype: numpy.ndarray
    """
    return calcDELTA(delta_arm) @ calcGAMMA(gamma_arm)


def armAdjustedParam(param, arm, reference=None, deltaChi=0.0):
    """pyFAI geometry of a detector whose arm has moved.

    A detector arm rotates the detector rigidly about the sample, which leaves
    the sample-to-PONI distance and the detector-internal PONI offsets alone, so
    only the three rotations change. This is a pure function of its arguments:
    it is called per frame, concurrently, against a geometry shared read-only
    between workers, and must never mutate anything.

    :param param:
        Home geometry ``[dist, poni1, poni2, rot1, rot2, rot3]`` in m and rad,
        i.e. the calibrated geometry at the arm position given by ``reference``.
    :param numpy.ndarray arm:
        3x3 rotation of the detector body in the Vlieg laboratory frame, from
        :func:`armRotation`.
    :param numpy.ndarray reference:
        3x3 arm rotation at which ``param`` was calibrated. ``None`` means the
        identity, i.e. a calibration taken with the arm at zero.
    :param float deltaChi: azimuthal reference of the calibration, in rad.
    :returns: ``(dist, poni1, poni2, rot1, rot2, rot3)`` in m and rad.
    :rtype: tuple
    :raises ValueError: At the gimbal lock of the pyFAI rotations.
    """
    param = np.asarray(param, dtype=float)
    arm = np.asarray(arm, dtype=float)
    if reference is not None:
        arm = arm @ np.asarray(reference, dtype=float).T
    relabel = azimuthRelabel(deltaChi)
    rotation = relabel.T @ arm @ relabel @ pyFAIRotationMatrix(*param[3:6])
    return (*param[:3], *pyFAIRotationAngles(rotation))


class Detector2D_SXRD(geometry.Geometry):
    """pyFAI detector geometry carrying the SXRD angle conventions.

    The calibrated pyFAI geometry describes the detector at the *home* arm
    position, recorded by :meth:`setArmReference`. A detector arm is a scanned
    motor, so its per-frame position is **not** stored here: it is passed to the
    conversions as ``gamma_arm`` / ``delta_arm``, and turned into a geometry by
    :meth:`paramAtArm`. This object is shared read-only between reconstruction
    workers, so nothing about the moving arm may become mutable state on it.
    """

    def __init__(self, *args, **keyargs):
        super().__init__(*args, **keyargs)
        self._R_arm_reference = np.identity(3)
        self.setAzimuthalReference(0)
        self.setPolarization(0, 0)

    # ------------------------------------------------------------------
    # Detector arm
    # ------------------------------------------------------------------

    def setArmReference(self, R_arm=None, *, gamma_arm=None, delta_arm=None):
        """Record the arm position at which this geometry was calibrated.

        The calibration is only ever interpreted relative to this reference, so
        the arm's absolute motor zero never has to mean anything and the direct
        beam never has to be on the detector. Defaults to the identity, i.e. a
        calibration taken at ``(gamma_arm, delta_arm) = (0, 0)``.

        :param numpy.ndarray R_arm:
            3x3 rotation of the detector body at calibration time. Mutually
            exclusive with the angle arguments.
        :param float gamma_arm:
            Out-of-plane arm angle ``gamma_p`` at calibration, in rad.
        :param float delta_arm:
            In-plane arm angle ``delta_p`` at calibration, in rad.
        :raises ValueError: If the rotation is not a proper rotation matrix.
        """
        if R_arm is not None and (gamma_arm is not None or delta_arm is not None):
            raise ValueError(
                "Give the arm reference either as a matrix or as arm angles, not both"
            )
        if R_arm is None:
            if gamma_arm is None and delta_arm is None:
                raise ValueError(
                    "setArmReference requires a matrix or at least one arm angle"
                )
            R_arm = armRotation(gamma_arm or 0.0, delta_arm or 0.0)
        self._R_arm_reference = _validatedRotation(R_arm, "arm reference")

    def getArmReference(self):
        """Arm rotation at which this geometry was calibrated.

        :rtype: numpy.ndarray
        """
        return self._R_arm_reference.copy()

    def _armRelative(self, gamma_arm, delta_arm):
        """Rotation from the calibrated arm position to the requested one.

        :returns: 3x3 matrix, or ``None`` for "at the calibrated position".
        :rtype: numpy.ndarray or None
        :raises ValueError: If exactly one of the two angles is given.
        """
        if gamma_arm is None and delta_arm is None:
            return None
        if gamma_arm is None or delta_arm is None:
            raise ValueError(
                "Give both gamma_arm and delta_arm, or neither. One alone is "
                "ambiguous, because the calibration reference is a rotation "
                "rather than a pair of independent offsets."
            )
        return armRotation(gamma_arm, delta_arm) @ self._R_arm_reference.T

    def paramAtArm(self, gamma_arm=None, delta_arm=None):
        """pyFAI geometry of this detector with its arm at the given position.

        Pure: it neither mutates this object nor touches its caches, so it is
        safe to call per frame and from several workers at once.

        ``None`` means *at the calibrated arm position*, which returns the home
        geometry unchanged. Note that this is not the same as passing ``0.0``
        unless the calibration was taken with the arm at zero -- ``None`` is
        "wherever the calibration was taken", ``0.0`` is "the motor reads zero".

        :param float gamma_arm: out-of-plane arm angle ``gamma_p`` in rad.
        :param float delta_arm: in-plane arm angle ``delta_p`` in rad.
        :returns: ``(dist, poni1, poni2, rot1, rot2, rot3)`` in m and rad.
        :rtype: tuple
        :raises ValueError: At the gimbal lock of the pyFAI rotations.
        """
        relative = self._armRelative(gamma_arm, delta_arm)
        if relative is None:
            return tuple(np.asarray(self.param[:6], dtype=float))
        return armAdjustedParam(self.param, relative, deltaChi=self._deltaChi)

    def toNXdict(self):
        """To be used with silx.io.dictdump.dicttonx
        to save and load the data from any nexus file.

        """
        pyFAI_dict = _nx_to_python(self.get_config())
        detector_config = dict(pyFAI_dict.get("detector_config", {}))
        orientation = detector_config.get("orientation")
        if hasattr(orientation, "value"):
            detector_config["orientation"] = int(orientation.value)
        pyFAI_dict["detector_config"] = detector_config

        if hasattr(self, "_roi"):
            pyFAI_dict["poni1"] = self.poni1_max
            pyFAI_dict["poni2"] = self.poni2_max
        if self.detector.shape is not None:
            pyFAI_dict["shape"] = np.asarray(self.detector.shape, dtype=np.int64)
        if self.detector.max_shape is not None:
            pyFAI_dict["max_shape"] = np.asarray(
                self.detector.max_shape, dtype=np.int64
            )

        nxdict = {
            "detector_SXRD": {
                "config": {**pyFAI_dict},
                "title": "pyFAI detector calibration",
                "@NX_class": "NXcollection",
                "azimuth": self._deltaChi,
                "polarization": self._polFactor,
                "polarization_axis": self._polAxis,
                "binning": np.array(self.detector.binning),
                # The arm position this calibration was taken at. The scanned
                # arm itself belongs to the scan, not here.
                "arm_reference": np.asarray(self._R_arm_reference),
            },
            "@NX_class": "NXcollection",
            "@creator": f"datautils v {__version__}",
        }

        if hasattr(self, "_roi"):
            nxdict["detector_SXRD"]["roi"] = np.array(self._roi)
        return nxdict

    def fromNXdict(self, nxdict):
        """Reads detector config from nexus dict written by toNXdict().
        Use silx.io.dictdump.nxtodict to read these nexus files.

        First searches for a "detector_SXRD" entry in the passed dict, otherwise
        assumes that the config is in the current group.
        """
        if "detector_SXRD" in nxdict:
            detdict = nxdict["detector_SXRD"]
        else:
            detdict = nxdict
        config = _nx_to_python(dict(detdict["config"]))
        for entry in list(config):
            if entry.startswith("@"):
                del config[entry]
        shape = config.pop("shape", None)
        max_shape = config.pop("max_shape", None)
        if "detector" in config:
            config["detector"] = str(config["detector"])
        self.set_config(config)
        if max_shape is not None:
            self.detector.max_shape = tuple(max_shape)
        if shape is not None:
            self.detector.shape = tuple(shape)
        self.setAzimuthalReference(detdict["azimuth"])
        self.setPolarization(detdict["polarization_axis"], detdict["polarization"])
        # Calibrations written before the moveable detector arm carry no
        # reference; an absent one is the identity, i.e. calibrated at arm zero.
        if "arm_reference" in detdict:
            self.setArmReference(np.asarray(detdict["arm_reference"], dtype=float))
        if "binning" in detdict:
            self.detector.binning = tuple(detdict["binning"])
        if "roi" in detdict:
            roi = detdict["roi"]
            self.set_roi(list(roi[0]), list(roi[1]))
        return self

    def copy_bin_roi_applied(self):
        """Returns a new Detector2D_SXRD instance, where the current binning and shape
        have been set to the source image file. This is intended for further binning
        of images which already have been processed with ``pply_bin_roi_to_image` once.
        """
        new_det = Detector2D_SXRD()
        ddict = self.toNXdict()
        new_det.fromNXdict(ddict)
        binns = new_det.detector.binning
        new_det.detector.max_shape = new_det.detector.shape
        new_det.poni1_max = new_det.poni1
        new_det.poni2_max = new_det.poni2
        new_det.detector.pixel1 *= binns[0]
        new_det.detector.pixel2 *= binns[1]
        new_det.detector.binning = (1, 1)
        new_det.set_roi([0, new_det.detector.shape[0]], [0, new_det.detector.shape[1]])
        new_det.reset()
        return new_det

    def __bin_img(self, dat):
        dim1, dim2 = self.detector.binning
        if dim1 == 1 and dim2 == 1:
            return dat
        t1 = dat.shape[0] % dim1
        t2 = dat.shape[1] % dim2
        if t1 == 0:
            t1 = None
        if t2 == 0:
            t2 = None
        dat_binned = binning(dat[:-t1, :-t2], (dim1, dim2))
        return dat_binned

    def set_roi(self, range1, range2):
        self._roi = range1, range2
        offset1 = range1[0]
        offset2 = range2[0]
        self.roi_range = range1, range2
        if not hasattr(self, "poni1_max"):
            self.poni1_max = self.poni1
        if not hasattr(self, "poni2_max"):
            self.poni2_max = self.poni2
        self.poni1 = self.poni1_max - offset1 * self.detector.pixel1
        self.poni2 = self.poni2_max - offset2 * self.detector.pixel2
        self.detector.shape = (range1[1] - range1[0], range2[1] - range2[0])
        self.reset()

    def apply_bin_roi_to_image(self, img):
        mask = np.isfinite(img).astype(np.int32)
        img_binned = self.__bin_img(np.nan_to_num(img))
        img_cnts = self.__bin_img(mask).astype(np.float64)
        img_cnts[~np.isfinite(img_cnts)] = np.nan
        img_binned *= np.prod(self.detector.binning) / img_cnts
        if hasattr(self, "roi_range"):
            range1, range2 = self.roi_range
            return img_binned[range1[0] : range1[1], range2[0] : range2[1]]
        else:
            return img_binned

    def getBornAgain(self):
        """Only for surface facing upwards (along y), i.e. azimuth = 90°."""
        ba_config = {}

        ba_config["NbinsX"] = self.detector.shape[1]
        ba_config["Width"] = (self.detector.shape[1] * self.detector.pixel2) * 1e3  # mm
        ba_config["NbinsY"] = self.detector.shape[0]
        ba_config["Height"] = (
            self.detector.shape[0] * self.detector.pixel1
        ) * 1e3  # mm
        ba_config["wavelength"] = self.get_wavelength() * 1e9  # nm
        _cal = self.getFit2D()
        ba_config["Distance"] = _cal["directDist"]  # mm
        ba_config["v0(vertBeamPos)"] = (
            (self.detector.shape[0] - _cal["centerY"]) * self.detector.pixel1 * 1e3
        )  # mm
        ba_config["u0(horizontBeamPos)"] = (
            _cal["centerX"] * self.detector.pixel2 * 1e3
        )  # mm
        return ba_config

    # in rad!!!
    def setAzimuthalReference(self, deltaChi):
        self._deltaChi = copy.copy(deltaChi)
        self._cached_array = {}

    def getAzimuthalReference(self):
        return self._deltaChi

    def setPolarization(self, angle, factor):
        self._polAxis = angle
        self._polFactor = factor
        self._cached_array = {}

    def getPolarization(self):
        return self._polAxis, self._polFactor

    def polarizationAtPoints(self, x, y, alpha_i, gamma_arm=None, delta_arm=None):
        r"""Polarization correction at the given pixels, for one frame.

        Evaluates the z-axis expression of the ANA/ROD manual directly in
        the surface-frame angles,

        .. math::

            P = p_h \left[1 - (\sin\alpha \cos\delta \cos\gamma
                + \cos\alpha \sin\gamma)^2\right]
              + (1 - p_h)\left[1 - \sin^2\delta \cos^2\gamma\right]

        with :math:`p_h` the fraction of horizontally polarized light.

        Unlike :meth:`polarizationArray`, this follows the detector arm.
        pyFAI derives its scattering angles from the calibrated geometry
        alone, so with a moving arm it reports only the angles the detector
        subtends at its calibration position and the correction comes out
        far too small: on a specular scan reaching a scattering angle of 18
        degrees it gives 0.6% where the true correction is 10%.

        :param x: detector coordinates along pyFAI dimension 1, in pixels.
        :param y: detector coordinates along pyFAI dimension 2, in pixels.
        :param float alpha_i: incidence angle of this frame, in rad.
        :param float gamma_arm, delta_arm: detector arm position of this
            frame as true scattering angles, in rad; ``None`` is the
            calibrated position.
        :returns: The polarization correction at each pixel.
        :rtype: numpy.ndarray
        """
        gamma, delta = self.surfaceAnglesPoint(x, y, alpha_i, gamma_arm, delta_arm)
        fraction = self._polFactor
        p_hor = (
            1.0
            - (
                np.sin(alpha_i) * np.cos(delta) * np.cos(gamma)
                + np.cos(alpha_i) * np.sin(gamma)
            )
            ** 2
        )
        p_ver = 1.0 - (np.sin(delta) ** 2) * (np.cos(gamma) ** 2)
        return fraction * p_hor + (1.0 - fraction) * p_ver

    def polarizationArray(self, shape=None):
        r"""Polarization correction of every detector pixel.

        pyFAI evaluates the polarization in its own detector frame, so both
        of its parameters have to be translated from orGUI's diffractometer
        conventions before they are handed over:

        *Orientation.* pyFAI's ``axis_offset`` places the polarization node
        at its own azimuth ``chi = -axis_offset``, while orGUI measures
        azimuth in the diffractometer frame, ``azimuth = chi + deltaChi``
        (see :meth:`primBeamAngles`). The node therefore has to be requested
        at ``axis_offset = deltaChi - polAxis``. Leaving the azimuthal
        reference out rotates the polarization plane by that angle -- for
        the common vertical-scattering setup, ``deltaChi = 90`` degrees,
        which swaps the polarized and unpolarized directions entirely.

        *Amount.* ``polarization_factor`` is the **fraction of horizontal
        polarization** used by :meth:`.HKLVlieg.UBCalculator.polarization`
        and by the ANA/ROD manual, whereas pyFAI's ``factor`` is
        ``(Ih - Iv) / (Ih + Iv)``. The two are related by
        ``factor = 2 * fraction - 1``; they agree only for a fully
        horizontally polarized beam.

        With both translations applied this reproduces
        :meth:`.HKLVlieg.UBCalculator.polarization`, the z-axis expression
        of the ANA/ROD manual, to machine precision.

        :param shape: Detector shape; taken from the detector when omitted.
        :returns: 2D polarization correction array.
        :rtype: numpy.ndarray
        """
        return self.polarization(
            shape=shape,
            factor=2.0 * self._polFactor - 1.0,
            axis_offset=self._deltaChi - self._polAxis,
        )

    def primBeamAngles(self, shape=None):
        """gives angles in laboratory reference frame."""
        if (
            self._cached_array.get("gamma_p") is None
            or self._cached_array.get("delta_p") is None
        ):
            azimuth = (
                self.center_array(shape, unit=pyFAI.units.CHI_RAD) + self._deltaChi
            )
            tth = self.center_array(shape, unit=pyFAI.units.TTH_RAD)
            sintth = np.sin(tth)
            costth = np.cos(tth)

            delta_p = np.arctan2(sintth * np.sin(azimuth), costth)
            gamma_p = np.arctan2(sintth * np.cos(azimuth) * np.cos(delta_p), costth)

            self._cached_array["gamma_p"] = gamma_p
            self._cached_array["delta_p"] = delta_p

        return self._cached_array["gamma_p"], self._cached_array["delta_p"]

    # numpy array!!!
    def primBeamPoints(self, x, y, gamma_arm=None, delta_arm=None):
        """Laboratory-frame detector angles of the given pixels.

        :param x: detector coordinates along pyFAI dimension 1, in pixels.
        :param y: detector coordinates along pyFAI dimension 2, in pixels.
        :param float gamma_arm, delta_arm:
            Detector arm position of this frame, as true scattering angles
            ``gamma_p``/``delta_p`` in rad. ``None`` (both) means the
            calibrated arm position, which takes the original code path and
            reproduces a static detector exactly. See :meth:`paramAtArm`.
        :returns: ``(gamma_p, delta_p)`` in rad.
        :rtype: tuple
        """
        relative = self._armRelative(gamma_arm, delta_arm)
        if relative is None:
            azimuth = self.chi(x, y) + self._deltaChi
            tth = self.tth(x, y)
        else:
            tth, azimuth = self._tthAzimuthAtArm(x, y, gamma_arm, delta_arm)

        return _primBeamFromTthAzimuth(tth, azimuth)

    def _tthAzimuthAtArm(self, x, y, gamma_arm, delta_arm):
        """Scattering angles of given pixels with the arm at a given position.

        :returns: ``(tth, azimuth)`` in rad, azimuth already offset by the
            azimuthal reference.
        :rtype: tuple
        """
        t3, t1, t2 = self.calc_pos_zyx(
            d1=np.asarray(x),
            d2=np.asarray(y),
            param=np.asarray(self.paramAtArm(gamma_arm, delta_arm), dtype=float),
            do_parallax=True,
        )
        tth = np.arctan2(np.sqrt(t1 * t1 + t2 * t2), t3)
        return tth, np.arctan2(t1, t2) + self._deltaChi

    def surfaceAngles(self, alpha_i, shape=None):
        """Angles in the reference frame, where the crystal is tilted by alpha_i.

        :param float alpha_i: incidence angle in rad.
        :param shape: detector shape, defaults to the calibrated one.
        :returns: ``(gamma, delta)`` surface-frame angles in rad.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]

        ``delta`` is recovered with an arcsin and is therefore only defined
        for ``|delta| < pi/2``. A detector large enough, or close enough, to
        reach 90 deg in-plane saturates instead of wrapping, and the
        conversion stops being invertible by :meth:`pixelsSurfaceAngles`.
        """
        if hasattr(self, "_alpha_i"):
            if self._alpha_i == alpha_i:
                if (
                    self._cached_array.get("gamma") is not None
                    and self._cached_array.get("delta") is not None
                ):
                    return self._cached_array.get("gamma"), self._cached_array.get(
                        "delta"
                    )
        self._alpha_i = alpha_i
        gamma_p, delta_p = self.primBeamAngles(shape)
        gamma = np.arcsin(
            np.cos(alpha_i) * np.sin(gamma_p)
            - np.sin(alpha_i) * np.cos(delta_p) * np.cos(gamma_p)
        )
        delta = np.arcsin(
            (np.sin(delta_p) * np.cos(gamma_p)) / np.cos(gamma)
        )  # evil, revise!!!
        self._cached_array["gamma"] = gamma
        self._cached_array["delta"] = delta
        return self._cached_array.get("gamma"), self._cached_array.get("delta")

    def crystalAngles(self, alpha_i, refraction_index, shape=None):
        if hasattr(self, "_n_refr"):
            if self._n_refr == refraction_index:
                if self._cached_array.get("gamma_cry") is not None and hasattr(
                    self, "_alpha_i_ref"
                ):
                    return (
                        self._cached_array.get("gamma_cry"),
                        self._cached_array.get("delta"),
                        self._alpha_i_ref,
                    )
        self._n_refr = refraction_index
        gamma, delta = self.surfaceAngles(alpha_i, shape)
        self._cached_array["gamma_cry"] = crystalAngles_singleArray(gamma, self._n_refr)
        self._alpha_i_ref = crystalAngles_singleArray(self._alpha_i, self._n_refr)
        return (
            self._cached_array.get("gamma_cry"),
            self._cached_array.get("delta"),
            self._alpha_i_ref,
        )

    def surfaceAnglesPoint(self, x, y, alpha_i, gamma_arm=None, delta_arm=None):
        """Surface-frame detector angles of the given pixels.

        :param x: detector coordinates along pyFAI dimension 1, in pixels.
        :param y: detector coordinates along pyFAI dimension 2, in pixels.
        :param alpha_i: incidence angle in rad.
        :param float gamma_arm, delta_arm:
            Detector arm position of this frame, as true scattering angles
            ``gamma_p``/``delta_p`` in rad; ``None`` means the calibrated
            position. See :meth:`paramAtArm`.
        :returns: ``(gamma, delta)`` in rad.
        :rtype: tuple
        """
        gamma_p, delta_p = self.primBeamPoints(x, y, gamma_arm, delta_arm)
        gamma = np.arcsin(
            np.cos(alpha_i) * np.sin(gamma_p)
            - np.sin(alpha_i) * np.cos(delta_p) * np.cos(gamma_p)
        )
        delta = np.arcsin(
            (np.sin(delta_p) * np.cos(gamma_p)) / np.cos(gamma)
        )  # evil, revise!!!
        return gamma, delta

    def surfaceAnglesPointParam(
        self, x, y, alpha_i, param, gamma_arm=None, delta_arm=None
    ):
        """Calculate surface-frame detector angles for trial geometry.

        This follows the pyFAI refinement pattern by passing trial geometry
        parameters directly into the coordinate calculation. It does not
        modify this geometry object or reset its caches.

        :param numpy.ndarray x:
            Detector coordinates along pyFAI dimension 1, in pixels.
        :param numpy.ndarray y:
            Detector coordinates along pyFAI dimension 2, in pixels.
        :param numpy.ndarray alpha_i:
            Incidence angles in rad.
        :param numpy.ndarray param:
            Trial **home** pyFAI geometry
            ``[dist, poni1, poni2, rot1, rot2, rot3]`` in m and rad, i.e. the
            geometry at the calibrated arm position. The arm below is folded
            into it, so a refinement fits the home geometry while each
            reflection keeps the arm position it was measured at.
        :param float gamma_arm, delta_arm:
            Arm position of the frame each point was measured on, in rad;
            ``None`` means the calibrated position.
        :returns:
            ``(gamma, delta)`` surface-frame detector angles in rad.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """
        relative = self._armRelative(gamma_arm, delta_arm)
        if relative is not None:
            param = armAdjustedParam(param, relative, deltaChi=self._deltaChi)
        t3, t1, t2 = self.calc_pos_zyx(
            d1=np.asarray(x),
            d2=np.asarray(y),
            param=np.asarray(param, dtype=float),
            do_parallax=True,
        )
        tth = np.arctan2(np.sqrt(t1 * t1 + t2 * t2), t3)
        azimuth = np.arctan2(t1, t2) + self._deltaChi
        delta_p = np.arctan2(np.sin(tth) * np.sin(azimuth), np.cos(tth))
        gamma_p = np.arctan2(
            np.sin(tth) * np.cos(azimuth) * np.cos(delta_p),
            np.cos(tth),
        )
        alpha_i = np.asarray(alpha_i)
        gamma = np.arcsin(
            np.cos(alpha_i) * np.sin(gamma_p)
            - np.sin(alpha_i) * np.cos(delta_p) * np.cos(gamma_p)
        )
        delta = np.arcsin(np.sin(delta_p) * np.cos(gamma_p) / np.cos(gamma))
        return gamma, delta

    def crystalAnglesPoint(
        self, x, y, alpha_i, refraction_index, gamma_arm=None, delta_arm=None
    ):
        """Refraction-corrected detector angles of the given pixels.

        Used by stationary-scan integration, where a detector tracking the
        reflection keeps it at fixed ``(x, y)`` while the exit angles at that
        pixel, and hence its ``hkl``, change frame by frame.

        :param float gamma_arm, delta_arm:
            Detector arm position of this frame, as true scattering angles
            ``gamma_p``/``delta_p`` in rad; ``None`` means the calibrated
            position.
        :returns: ``(gamma_cry, delta, alpha_i_cry)`` in rad.
        :rtype: tuple
        """
        gamma, delta = self.surfaceAnglesPoint(x, y, alpha_i, gamma_arm, delta_arm)
        gamma_cry = crystalAngles_singleArray(gamma, refraction_index)
        alpha_i_cry = crystalAngles_singleArray(alpha_i, refraction_index)
        return gamma_cry, delta, alpha_i_cry

    def pixelsTthChi(self, tth, chi, gamma_arm=None, delta_arm=None):
        """Detector coordinates of the given pyFAI scattering angles.

        Inverse of the pyFAI forward direction, and therefore the exact
        inverse of :meth:`primBeamPoints` and :meth:`surfaceAnglesPoint`.

        :param numpy.ndarray tth: scattering angles in rad.
        :param numpy.ndarray chi: pyFAI azimuthal angles in rad.
        :param float gamma_arm, delta_arm:
            Detector arm position of this frame, as true scattering angles
            ``gamma_p``/``delta_p`` in rad; ``None`` means the calibrated
            position. A given scattering direction lands on a
            different pixel once the arm moves, which is what lets a rocking
            curve at fixed exit angles be extracted from a moving detector.
        :returns:
            Array shaped ``(*tth.shape, 2)`` holding pixel coordinates as
            ``(dim1, dim2)``, i.e. ``(slow, fast)``, using the same
            pixel-centre indexing that :meth:`primBeamPoints` consumes.
            Coordinates outside the detector are returned unclipped.
        :rtype: numpy.ndarray

        The inversion assumes a uniform pixel pitch. That holds for the
        regular and module detectors pyFAI models by a uniform grid,
        including ones whose gaps live only in the mask, but not for a
        detector supplying a genuinely non-uniform ``get_pixel_corners``.
        """
        tth = np.atleast_1d(np.asarray(tth))
        chi = np.atleast_1d(np.asarray(chi))
        shape = tth.shape
        assert tth.shape == chi.shape

        # from here everything flat
        # tanchi = np.tan(chi.flatten())
        # nu = np.tan(tth.flatten()) * np.sin(chi.flatten())

        # M = np.empty((np.prod(shape), 2, 2))
        # A = np.empty((np.prod(shape), 2))
        chi = chi.flatten()
        tth = tth.flatten()

        sinchi = np.sin(chi)
        coschi = np.cos(chi)

        sintth = np.sin(tth)
        costth = np.cos(tth)

        nu = sintth * sinchi

        # Detector rotation with the arm where this frame had it. ``dist`` and
        # the PONI offsets used below are unchanged by the arm, which rotates
        # the detector rigidly about the sample.
        R = pyFAIRotationMatrix(*self.paramAtArm(gamma_arm, delta_arm)[3:6])

        a = R[0, 0] * coschi - R[1, 0] * sinchi
        b = R[0, 1] * coschi - R[1, 1] * sinchi
        c = R[0, 0] * costth - R[2, 0] * nu
        d = R[0, 1] * costth - R[2, 1] * nu

        A1 = (R[1, 2] * sinchi - R[0, 2] * coschi) * self.dist
        A2 = (R[2, 2] * nu - R[0, 2] * costth) * self.dist

        """
        The rest is a fast way to solve
        M @ ptilde = A
        with

        M = np.empty((np.prod(shape), 2, 2))
        A = np.empty((np.prod(shape), 2))

        A[:, 0] = A1
        A[:, 1] = A2

        M[:, 0, 0] = a
        M[:, 0, 1] = b
        M[:, 1, 0] = c
        M[:, 1, 1] = d

        linalg.solve doesn't work since singular matrices are present

        linalg.lstsq - solution: extremely slow!!!
        (> 10 s for Pilatus 2M pixels ~1600*1400 )

        ptilde = np.empty((np.prod(shape), 2))
        for i in range(a.size):
            ptilde[i] = LA.lstsq(M[i],A[i], rcond=None)[0]

        # Direct SVD:
        (about 10 s for Pilatus 2M pixels ~1600*1400)
        M_pinv = LA.pinv(M)
        ptilde_svd = np.empty((np.prod(shape), 2))
        for i in range(a.size):
            ptilde_svd[i] = M_pinv[i] @ A[i]

        so have to do it explicitly and find failing solutions manually:
        (about 1 s for Pilatus 2M pixels ~1600*1400)
        """
        determ = b * c - a * d

        ptilde1 = (
            b * A2 - d * A1
        ) / determ  # this raises divide by 0 warnings, not sure how to suppress them
        ptilde2 = (
            c * A1 - a * A2
        ) / determ  # this raises divide by 0 warnings, not sure how to suppress them

        mask = np.isclose(chi, 0.0, atol=1e-12)  # solution fails at chi = 0
        # Special case chi = 0
        if mask.any():
            A = np.zeros((3, np.count_nonzero(mask)))
            A[1, :] = sintth[mask]
            A[2, :] = costth[mask]

            X = R.T @ A

            s = X[2] / self.dist

            ptilde1[mask] = X[0] / s
            ptilde2[mask] = X[1] / s

        p = (
            np.column_stack((ptilde1, ptilde2)) + np.array([self.poni1, self.poni2])
        ) / np.array([self.pixel1, self.pixel2])

        # ``calc_cartesian_positions`` places pixel index ``i`` at its centre,
        # ``(i + 0.5) * pitch``, so dividing a distance by the pitch above
        # yields a pixel edge coordinate. Shift onto the same pixel-centre
        # indexing the forward direction consumes, which makes this the exact
        # inverse of :meth:`primBeamPoints` rather than a half-pixel off it.
        p -= 0.5

        return p.reshape((*shape, 2))

    @property
    def Qmax(self):
        """in A-1"""
        qmin, qmax = self.Qrange
        return qmax

    @property
    def Qmin(self):
        """in A-1"""
        qmin, qmax = self.Qrange
        return qmin

    @property
    def Qrange(self):
        """Momentum transfer the detector reaches, at the calibrated arm.

        :returns: ``(Qmin, Qmax)`` in Angstrom^-1.
        :rtype: tuple[float, float]

        See :meth:`qRangeAtArm` for a detector whose arm has moved: a moving
        arm changes what the detector can reach, so this value goes stale.
        """
        return self.qRangeAtArm()

    def qRangeAtArm(self, gamma_arm=None, delta_arm=None):
        """Momentum transfer the detector reaches with its arm at a position.

        :param float gamma_arm, delta_arm:
            Detector arm position, as true scattering angles in rad; ``None``
            means the calibrated position.
        :returns: ``(Qmin, Qmax)`` in Angstrom^-1.
        :rtype: tuple[float, float]

        Swinging the arm out moves the whole accessible band to higher ``Q``,
        so a range taken at the calibrated position both under-reports what a
        moved detector reaches and wrongly keeps the low end at zero.
        """
        xx, yy = self._edge_pixcoord()
        param = np.asarray(self.paramAtArm(gamma_arm, delta_arm), dtype=float)
        # qFunction returns Q at pixel center
        Q = self.qFunction(xx, yy, param=param) / 10.0

        beam = self.pixelsPrimeBeam(0.0, 0.0, gamma_arm, delta_arm)[0]
        # beam on detector ?
        if (
            0 <= beam[1] <= self.detector.shape[1]
            and 0 <= beam[0] <= self.detector.shape[0]
        ):
            return 0.0, float(np.amax(Q))
        else:
            return float(np.amin(Q)), float(np.amax(Q))

    @property
    def rangegamdel_p(self):
        """Is not identical to rangegamdel_p_full_det, where each pixel is calculated, but close within 1e-5

        Not sure why...
        """  # noqa: E501
        xx, yy = self._edge_pixcoord()
        gamma_p, delta_p = self.primBeamPoints(xx, yy)
        delrange = np.amin(delta_p), np.amax(delta_p)
        gamrange = np.amin(gamma_p), np.amax(gamma_p)
        return gamrange, delrange

    def rangegamdel(self, alpha_i):
        """Is not identical to _rangegamdel_full_det, where each pixel is calculated, but close within 1e-5

        Not sure why...
        """  # noqa: E501
        xx, yy = self._edge_pixcoord()
        gamma, delta = self.surfaceAnglesPoint(xx, yy, alpha_i)
        delrange = np.amin(delta), np.amax(delta)
        gamrange = np.amin(gamma), np.amax(gamma)
        return gamrange, delrange

    @property
    def _rangegamdel_p_full_det(self):
        gamma_p, delta_p = self.primBeamAngles()
        delrange = np.amin(delta_p), np.amax(delta_p)
        gamrange = np.amin(gamma_p), np.amax(gamma_p)
        return gamrange, delrange

    def _rangegamdel_full_det(self, alpha_i):
        gamma_p, delta_p = self.surfaceAngles(alpha_i)
        delrange = np.amin(delta_p), np.amax(delta_p)
        gamrange = np.amin(gamma_p), np.amax(gamma_p)
        return gamrange, delrange

    def pixelsPrimeBeam(self, gamma_p, delta_p, gamma_arm=None, delta_arm=None):
        """Detector coordinates of laboratory-frame detector angles.

        :param gamma_p, delta_p: laboratory-frame detector angles in rad.
        :param float gamma_arm, delta_arm:
            Detector arm position of this frame, as true scattering angles
            ``gamma_p``/``delta_p`` in rad; ``None`` means the calibrated
            position.
        :returns: pixel coordinates as ``(dim1, dim2)``.
        :rtype: numpy.ndarray
        """
        gamma_p = np.atleast_1d(np.asarray(gamma_p))
        delta_p = np.atleast_1d(np.asarray(delta_p))
        assert gamma_p.shape == delta_p.shape

        singam = np.sin(gamma_p)
        cosgam = np.cos(gamma_p)

        azimuth = np.arctan2(np.sin(delta_p) * cosgam, singam)
        tth = np.arctan2(singam, np.cos(delta_p) * cosgam * np.cos(azimuth))

        chi = azimuth - self._deltaChi

        return self.pixelsTthChi(tth, chi, gamma_arm, delta_arm)

    def pixelsSurfaceAngles(
        self, gamma, delta, alpha_i, gamma_arm=None, delta_arm=None
    ):
        """Detector coordinates of surface-frame detector angles.

        With the arm supplied, this is the pixel a given pair of exit angles
        falls on *for that frame*. A rocking curve is the intensity at fixed
        ``(delta, gamma)`` versus sample rotation, so on a moving detector the
        pixel it has to be read from moves with the arm.

        :param gamma, delta: surface-frame detector angles in rad.
        :param alpha_i: incidence angle in rad.
        :param float gamma_arm, delta_arm:
            Detector arm position of this frame, as true scattering angles
            ``gamma_p``/``delta_p`` in rad; ``None`` means the calibrated
            position.
        :returns: pixel coordinates as ``(dim1, dim2)``.
        :rtype: numpy.ndarray
        """
        gamma_p = np.arcsin(
            np.cos(alpha_i) * np.sin(gamma)
            + np.sin(alpha_i) * np.cos(delta) * np.cos(gamma)
        )
        delta_p = np.arcsin((np.sin(delta) * np.cos(gamma)) / np.cos(gamma_p))
        return self.pixelsPrimeBeam(gamma_p, delta_p, gamma_arm, delta_arm)

    def pixelsCrystalAngles(
        self,
        gamma_cry,
        delta,
        alpha_i_cry,
        refraction_index,
        gamma_arm=None,
        delta_arm=None,
    ):
        gamma = vacAngles_singleArray(gamma_cry, refraction_index)
        alpha_i = vacAngles_singleArray(alpha_i_cry, refraction_index)
        return self.pixelsSurfaceAngles(gamma, delta, alpha_i, gamma_arm, delta_arm)

    # for now numerical solution, only iterables!
    # in rad
    def pixelsPrimPoint(self, gamma_p, delta_p, shape=None):
        warnings.warn(
            "pixelsPrimPoint is deprecated and will be removed in the future. Use pixelsPrimeBeam instead.",  # noqa: E501
            FutureWarning,
        )
        """
        gamma_p_det, delta_p_det = self.primBeamAngles(shape)
        ig = invert_geometry.InvertGeometry(gamma_p_det,delta_p_det)
        xy = []
        for gam_p, del_p in zip(gamma_p,delta_p):
            xy.append(ig(gam_p,del_p))
        return xy
        """
        return self.pixelsPrimeBeam(np.asarray(gamma_p), np.asarray(delta_p))

    def pixelsSurfacePoint(self, gamma, delta, alpha_i, shape=None):
        warnings.warn(
            "pixelsSurfacePoint is deprecated and will be removed in the future. Use pixelsSurfaceAngles instead.",  # noqa: E501
            FutureWarning,
        )
        """
        gamma_det, delta_det = self.surfaceAngles(alpha_i,shape)
        ig = invert_geometry.InvertGeometry(gamma_det,delta_det)
        xy = []
        for gam, delt in zip(gamma,delta):
            xy.append(ig(gam,delt))
        return xy
        """
        return self.pixelsSurfaceAngles(
            np.asarray(gamma), np.asarray(delta), np.asarray(alpha_i)
        )

    def correctionArray(self, shape=None):
        if shape is None:
            shape = self.get_shape()
        if self._cached_array.get("corrarr") is not None:
            return self._cached_array.get("corrarr")
        else:
            self._cached_array["corrarr"] = self.solidAngleArray(
                shape
            ) * self.polarizationArray(shape)
        return self._cached_array.get("corrarr")

    def _edge_pixcoord(self):
        edge_1_x = np.arange(self.detector.shape[0])  # pixel center
        edge_1_y = np.zeros(self.detector.shape[0])
        edge_2_x = edge_1_x
        edge_2_y = np.full(self.detector.shape[0], self.detector.shape[1] - 1)
        edge_3_x = np.zeros(self.detector.shape[1])
        edge_3_y = np.arange(self.detector.shape[1])
        edge_4_x = np.full(self.detector.shape[1], self.detector.shape[0] - 1)
        edge_4_y = edge_3_y
        #                   left edge right edge top edge  bottom edge
        xx = np.concatenate([edge_1_x, edge_2_x, edge_3_x, edge_4_x])
        yy = np.concatenate([edge_1_y, edge_2_y, edge_3_y, edge_4_y])
        return xx, yy

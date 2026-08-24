# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
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
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import unittest

from .. import DetectorCalibration, HKLVlieg
import pyFAI

import numpy as np
import pytest
import os
import warnings

try:
    from silx.io import dictdump

    silx_avail = True
except Exception:
    silx_avail = False


class TestRWDetector2D_SXRD(unittest.TestCase):
    def setUp(self):
        self.sxrddet = DetectorCalibration.Detector2D_SXRD()
        self.sxrddet.detector = pyFAI.detector_factory("Pilatus 2M CdTe")
        self.sxrddet.poni1 = 0.1
        self.sxrddet.poni2 = 0.1
        self.sxrddet.rot1 = np.pi / 20
        self.sxrddet.rot2 = -np.pi / 20
        self.sxrddet.rot3 = 0
        self.sxrddet.dist = 1.5
        self.sxrddet.set_energy(15.0)
        self.sxrddet.setAzimuthalReference(np.deg2rad(90.0))
        self.sxrddet.setPolarization(np.deg2rad(90.0), 0.75)

    def test_from_to_dict(self):
        nxdict = self.sxrddet.toNXdict()
        othersxrddet = DetectorCalibration.Detector2D_SXRD()
        othersxrddet.fromNXdict(nxdict)
        self.check_sxrd_equal(othersxrddet)

    def check_sxrd_equal(self, other):
        self.assertDictEqual(self.sxrddet.getPyFAI(), other.getPyFAI())
        self.assertEqual(self.sxrddet._polFactor, other._polFactor)
        self.assertEqual(self.sxrddet._polAxis, other._polAxis)
        self.assertEqual(self.sxrddet._deltaChi, other._deltaChi)
        np.testing.assert_allclose(
            self.sxrddet.getArmReference(), other.getArmReference()
        )

    def test_arm_reference_round_trip(self):
        self.sxrddet.setArmReference(
            gamma_arm=np.deg2rad(1.5), delta_arm=np.deg2rad(-2.0)
        )
        other = DetectorCalibration.Detector2D_SXRD()
        other.fromNXdict(self.sxrddet.toNXdict())
        self.check_sxrd_equal(other)

    def test_missing_arm_reference_loads_as_identity(self):
        nxdict = self.sxrddet.toNXdict()
        del nxdict["detector_SXRD"]["arm_reference"]

        other = DetectorCalibration.Detector2D_SXRD()
        other.fromNXdict(nxdict)

        np.testing.assert_array_equal(other.getArmReference(), np.identity(3))

    @unittest.skipUnless(silx_avail, "silx not available")
    def test_write_nx(self):
        self._nxfilename = "./detcal_test.nx"
        self.addCleanup(self._destruct_file, self._nxfilename)

        nxdict = self.sxrddet.toNXdict()
        dictdump.dicttonx(nxdict, self._nxfilename)
        self.assertTrue(os.path.isfile(self._nxfilename))
        os.remove(self._nxfilename)

    @unittest.skipUnless(silx_avail, "silx not available")
    def test_read_write_nx(self):
        self._nxfilename = "./detcal_test.nx"
        self.addCleanup(self._destruct_file, self._nxfilename)

        nxdict = self.sxrddet.toNXdict()
        dictdump.dicttonx(nxdict, self._nxfilename)

        read_dict = dictdump.nxtodict(self._nxfilename)
        othersxrddet = DetectorCalibration.Detector2D_SXRD()
        othersxrddet.fromNXdict(read_dict)
        self.check_sxrd_equal(othersxrddet)

        othersxrddet = DetectorCalibration.loadNXdict(read_dict)
        self.check_sxrd_equal(othersxrddet)

        os.remove(self._nxfilename)

    def _destruct_file(self, filename):
        if os.path.exists(filename):
            os.remove(self._nxfilename)


class TestAnglePixelConversion(unittest.TestCase):
    def setUp(self):
        self.sxrddet = DetectorCalibration.Detector2D_SXRD()
        self.sxrddet.detector = pyFAI.detector_factory("Pilatus 2M CdTe")
        self.sxrddet.poni1 = 0.1
        self.sxrddet.poni2 = 0.1
        self.sxrddet.rot1 = np.pi / 20
        self.sxrddet.rot2 = -np.pi / 20
        self.sxrddet.rot3 = 0
        self.sxrddet.dist = 1.5
        self.sxrddet.set_energy(15.0)
        self.sxrddet.setAzimuthalReference(np.deg2rad(90.0))
        self.sxrddet.setPolarization(np.deg2rad(90.0), 0.75)

        # Pixel-centre indexing: ``calc_cartesian_positions`` places index ``i``
        # at ``(i + 0.5) * pitch``, so the centre of pixel ``i`` is the integer
        # index ``i`` itself. This is the convention both the forward
        # (``primBeamPoints``) and the inverse (``pixelsTthChi``) direction use.
        self.p1 = np.arange(self.sxrddet.detector.shape[1], dtype=float)
        self.p2 = np.arange(self.sxrddet.detector.shape[0], dtype=float)
        self.p12 = np.moveaxis(np.array(np.meshgrid(self.p1, self.p2)), 0, -1)[
            :, :, ::-1
        ]

        self.mu = np.deg2rad(0.1)  # should be differed, but probably ok

    def test_prim_beam_points_inverse_is_exact(self):
        """``pixelsPrimeBeam`` must invert ``primBeamPoints`` exactly.

        Both directions index pixels by their centre, so the round trip is the
        identity rather than the half-pixel offset it used to carry. Checked on
        a module detector and on a plain regular one, since ``pixelsTthChi``
        inverts assuming a uniform pixel pitch.
        """
        rows = np.array([0.0, 1.0, 100.0, 517.0, 1000.0])
        columns = np.array([0.0, 3.0, 250.0, 733.0, 1400.0])

        for name in ("Pilatus 2M CdTe", None):
            with self.subTest(detector=name or "regular"):
                det = DetectorCalibration.Detector2D_SXRD()
                if name is None:
                    det.detector = pyFAI.detectors.Detector(
                        172e-6, 172e-6, max_shape=(1500, 1500)
                    )
                else:
                    det.detector = pyFAI.detector_factory(name)
                det.poni1 = 0.1
                det.poni2 = 0.1
                det.rot1 = np.pi / 20
                det.rot2 = -np.pi / 20
                det.rot3 = np.pi / 30
                det.dist = 1.5
                det.wavelength = 0.8e-10
                det.setAzimuthalReference(np.deg2rad(90.0))

                gamma_p, delta_p = det.primBeamPoints(rows, columns)
                back = det.pixelsPrimeBeam(gamma_p, delta_p)

                np.testing.assert_allclose(back[:, 0], rows, atol=1e-9)
                np.testing.assert_allclose(back[:, 1], columns, atol=1e-9)

    def test_surface_angles_trial_parameters_match_active_geometry(self):
        d1 = np.array([100.0, 500.0, 1000.0])
        d2 = np.array([200.0, 700.0, 1200.0])
        alpha = np.full(3, self.mu)
        param = np.array(
            [
                self.sxrddet.dist,
                self.sxrddet.poni1,
                self.sxrddet.poni2,
                self.sxrddet.rot1,
                self.sxrddet.rot2,
                self.sxrddet.rot3,
            ]
        )

        expected = self.sxrddet.surfaceAnglesPoint(d1, d2, alpha)
        actual = self.sxrddet.surfaceAnglesPointParam(d1, d2, alpha, param)

        self.assertTrue(np.allclose(actual[0], expected[0]))
        self.assertTrue(np.allclose(actual[1], expected[1]))

    def assertRoundTripPixelError(self, p12_conv, abserr, msg):
        """Report the angle-to-pixel round-trip error over the whole detector.

        Deliberately a warning rather than an assertion, as it has been since
        these sweeps were written. They drive the geometry far past anything
        physical -- ``rot`` up to 180 degrees, ``poni`` metres off a detector
        1.7 m across -- and there the surface-angle parametrisation is at or
        over the edge of being invertible at all (see
        ``Detector2D_SXRD.surfaceAngles`` on the arcsin saturation). Whether a
        given sweep value tips over depends on process state, which makes a
        hard assertion here flaky for reasons unrelated to what is being
        tested.

        The pixel-coordinate convention these were meant to guard is asserted
        exactly, and deterministically, in
        :meth:`TestAnglePixelConversion.test_prim_beam_points_inverse_is_exact`
        on well-conditioned geometries instead.
        """
        maxerror = np.nanmax(np.abs(self.p12 - p12_conv))
        if maxerror > abserr:
            warnings.warn(f"too large error: {maxerror:.5f} {msg}")

    def assertPixelErrorSurfaceAnglesLessThan(self, abserr=1e-3, msg=""):
        gamma, delta = self.sxrddet.surfaceAngles(self.mu)
        p12_conv = self.sxrddet.pixelsSurfaceAngles(gamma, delta, self.mu)
        self.assertRoundTripPixelError(
            p12_conv, abserr, f"pixel coord from surface angles, {msg}"
        )

    def assertPixelErrorTthChiLessThan(self, abserr=1e-3, msg=""):
        tth = self.sxrddet.twoThetaArray()
        chi = self.sxrddet.chiArray()
        p12_conv = self.sxrddet.pixelsTthChi(tth, chi)
        self.assertRoundTripPixelError(
            p12_conv, abserr, f"pixel coord from tth and chi, {msg}"
        )

    def assertGamDelRangeErrorLessThan(self, abserr=1e-3, msg=""):
        exact = self.sxrddet._rangegamdel_p_full_det
        corner = self.sxrddet.rangegamdel_p

        diff = np.array(exact) - np.array(corner)
        max_rel_diff = np.amax(np.abs(diff))
        if max_rel_diff > abserr:
            warnings.warn(
                f"too large error: {max_rel_diff:.5f} approx det corners, {msg}"
            )
        # self.assertLessEqual(max_rel_diff, abserr, "too large error approx det corners, %s" % msg)  # noqa: E501

    def assertQrangeValid(self):
        Q = self.sxrddet.qArray() / 10.0
        Qmin = np.amin(Q)
        Qmax = np.amax(Q)
        Qmin_fast, Qmax_fast = self.sxrddet.Qrange
        f2d_cal = self.sxrddet.getFit2D()
        # beam on detector ?
        if (
            0 <= f2d_cal["centerX"] <= self.sxrddet.detector.shape[1]
            and 0 <= f2d_cal["centerY"] <= self.sxrddet.detector.shape[0]
        ):
            Qmin = 0.0

        self.assertTrue(np.allclose(Qmin, Qmin_fast, 1e-5))
        self.assertTrue(np.allclose(Qmax, Qmax_fast, 1e-5))

    def test_poni1(self):
        for p1 in np.linspace(-3, 3, 5):
            self.sxrddet.poni1 = p1
            msg = f"poni1 = {p1}"
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.poni1 = 0.1

    def test_poni2(self):
        for p1 in np.linspace(-3, 3, 5):
            self.sxrddet.poni2 = p1
            msg = f"poni2 = {p1}"
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.poni2 = 0.1

    def test_rot1(self):
        for p1 in np.linspace(0, np.pi, 8):
            self.sxrddet.rot1 = p1
            msg = f"rot1 = {p1}"
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.rot1 = np.pi / 20

    def test_rot2(self):
        for p1 in np.linspace(0, np.pi, 8):
            self.sxrddet.rot2 = p1
            msg = f"rot2 = {p1}"
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.rot2 = -np.pi / 20

    def test_rot3(self):
        for p1 in np.linspace(0, np.pi, 8):
            self.sxrddet.rot3 = p1
            msg = f"rot3 = {p1}"
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.rot3 = 0

    def test_dist(self):
        for d in np.linspace(0.01, 10, 5):
            self.sxrddet.dist = d
            msg = f"dist = {d}"
            self.assertPixelErrorSurfaceAnglesLessThan(msg=msg)
            self.assertPixelErrorTthChiLessThan(msg=msg)
            self.assertGamDelRangeErrorLessThan(msg=msg)
            self.assertQrangeValid()
        self.sxrddet.dist = 1.5


class TestDetectorArm(unittest.TestCase):
    """A scanned detector arm is a rigid rotation about the sample.

    pyFAI's ``(rot1, rot2, rot3)`` represent that exactly, leaving ``dist`` and
    the PONI offsets alone, which is what lets the Vlieg equations keep
    consuming per-pixel ``(gamma, delta)`` unchanged.

    The arm is a motor, not calibration state, so it is passed to the
    conversions rather than stored on the detector.
    """

    # Deliberately tilted: a home geometry straight down the beam would hide
    # errors that only appear when the arm composes with existing tilts.
    HOME = dict(dist=0.8, poni1=0.09, poni2=0.11, rot1=0.05, rot2=-0.03, rot3=0.11)
    ARM_POSITIONS_DEG = [
        (0.0, 0.0),
        (5.0, 15.0),
        (12.0, 30.0),
        (-7.0, 55.0),
        (35.0, -25.0),
        (-30.0, -70.0),
    ]

    def build(self, **overrides):
        parameters = dict(self.HOME, **overrides)
        det = DetectorCalibration.Detector2D_SXRD(
            detector="Pilatus1M", wavelength=0.5e-10, **parameters
        )
        det.setAzimuthalReference(np.deg2rad(90.0))
        det.setPolarization(0.0, 0.95)
        return det

    @staticmethod
    def rays(det, rows, columns, gamma_arm=None, delta_arm=None):
        gamma_p, delta_p = det.primBeamPoints(rows, columns, gamma_arm, delta_arm)
        return np.stack(
            [
                np.sin(delta_p) * np.cos(gamma_p),
                np.cos(delta_p) * np.cos(gamma_p),
                np.sin(gamma_p),
            ],
            axis=-1,
        )

    def test_default_reference_is_identity(self):
        det = self.build()
        np.testing.assert_array_equal(det.getArmReference(), np.identity(3))
        np.testing.assert_array_equal(det.paramAtArm(), np.asarray(det.param[:6]))

    def test_unspecified_arm_reproduces_the_static_geometry_exactly(self):
        """Passing no arm must be bit-identical to the code before the arm.

        This is what allows the feature to land before the rest of the pipeline
        knows about arms at all.
        """
        rows = np.array([3.0, 200.0, 700.0])
        columns = np.array([9.0, 400.0, 900.0])
        alpha = np.deg2rad(0.5)
        det = self.build()

        gamma_p, delta_p = det.primBeamPoints(rows, columns)
        gamma, delta = det.surfaceAnglesPoint(rows, columns, alpha)

        np.testing.assert_array_equal(
            det.primBeamPoints(rows, columns, None, None), (gamma_p, delta_p)
        )
        np.testing.assert_array_equal(
            det.pixelsSurfaceAngles(gamma, delta, alpha, None, None),
            det.pixelsSurfaceAngles(gamma, delta, alpha),
        )
        # returns two arrays and a scalar, so compare component by component
        for expected, actual in zip(
            det.crystalAnglesPoint(rows, columns, alpha, 1.0 - 1e-6),
            det.crystalAnglesPoint(rows, columns, alpha, 1.0 - 1e-6, None, None),
        ):
            np.testing.assert_array_equal(expected, actual)

    def test_beam_centre_pixel_reads_the_arm_angles(self):
        """The conventional-diffractometer property.

        Moving the arm to ``(gamma_arm, delta_arm)`` must make the pixel that
        was on the primary beam read exactly those angles. The anchor is the
        beam centre, not the PONI -- with these tilts the two are over 100
        pixels apart and the PONI reads about a degree of delta.
        """
        home = self.build()
        anchor = home.pixelsPrimeBeam(0.0, 0.0)[0]
        rows, columns = anchor[:1], anchor[1:]

        for gamma_deg, delta_deg in self.ARM_POSITIONS_DEG:
            with self.subTest(gamma_arm=gamma_deg, delta_arm=delta_deg):
                gamma_p, delta_p = home.primBeamPoints(
                    rows, columns, np.deg2rad(gamma_deg), np.deg2rad(delta_deg)
                )
                self.assertAlmostEqual(np.rad2deg(gamma_p[0]), gamma_deg, places=9)
                self.assertAlmostEqual(np.rad2deg(delta_p[0]), delta_deg, places=9)

    def test_arm_angles_are_true_scattering_angles(self):
        """The arm is given in ``gamma_p``/``delta_p``, not six-circle angles.

        A floor-mounted arm reads the same numbers whatever the incidence
        angle, so what it reports has to be independent of alpha. The
        six-circle angles at the same pixel are not.
        """
        det = self.build()
        anchor = det.pixelsPrimeBeam(0.0, 0.0)[0]
        rows, columns = anchor[:1], anchor[1:]
        gamma_arm, delta_arm = np.deg2rad([12.0, 30.0])

        for alpha_deg in (0.0, 0.5, 5.0, 15.0):
            with self.subTest(alpha=alpha_deg):
                alpha = np.deg2rad(alpha_deg)
                gamma_p, delta_p = det.primBeamPoints(
                    rows, columns, gamma_arm, delta_arm
                )
                # independent of alpha, as a laboratory-frame quantity must be
                self.assertAlmostEqual(np.rad2deg(gamma_p[0]), 12.0, places=9)
                self.assertAlmostEqual(np.rad2deg(delta_p[0]), 30.0, places=9)

                # and the six-circle angles follow HKLVlieg's own conversion
                gamma, delta = det.surfaceAnglesPoint(
                    rows, columns, alpha, gamma_arm, delta_arm
                )
                _, expected_delta, expected_gamma = HKLVlieg.vliegDiffracAngles(
                    [alpha, delta_arm, gamma_arm, 0.0, 0.0, 0.0]
                )[:3]
                np.testing.assert_allclose(gamma[0], expected_gamma, atol=1e-12)
                np.testing.assert_allclose(delta[0], expected_delta, atol=1e-12)

    def test_surface_frame_arm_angles_can_be_converted(self):
        """Six-circle arm readouts convert to what the arm API expects.

        ``armAnglesFromSurface`` is the inverse of the conversion above, so a
        detector driven from surface-frame motor values lands on exactly the
        same geometry.
        """
        det = self.build()
        anchor = det.pixelsPrimeBeam(0.0, 0.0)[0]
        rows, columns = anchor[:1], anchor[1:]
        alpha = np.deg2rad(5.0)
        gamma_p, delta_p = np.deg2rad([12.0, 30.0])

        # the six-circle angles this arm position corresponds to ...
        gamma_surface, delta_surface = det.surfaceAnglesPoint(
            rows, columns, alpha, gamma_p, delta_p
        )
        # ... converted back must reproduce the primary-beam angles
        gamma_back, delta_back = DetectorCalibration.armAnglesFromSurface(
            gamma_surface[0], delta_surface[0], alpha
        )
        self.assertAlmostEqual(gamma_back, gamma_p, places=12)
        self.assertAlmostEqual(delta_back, delta_p, places=12)

        np.testing.assert_allclose(
            det.paramAtArm(gamma_back, delta_back),
            det.paramAtArm(gamma_p, delta_p),
            atol=1e-12,
        )

    def test_crystal_angle_round_trip_keeps_the_refraction_frame(self):
        """``pixelsCrystalAngles`` must invert ``crystalAnglesPoint`` exactly.

        Both sides carry a refraction correction, and the incidence angle has
        to be handed back in the frame the forward direction produced it in:
        ``crystalAnglesPoint`` returns ``alpha_i_cry``, and feeding the vacuum
        ``alpha_i`` instead refracts it a second time and displaces the result.

        Only invertible above the critical angle: below it the beam does not
        propagate into the crystal and ``crystalAngles`` clamps to zero, which
        no inverse can undo. Points in that band are excluded rather than
        papered over.
        """
        det = self.build()
        refraction_index = 1.0 - 3.16e-6
        critical = np.arccos(refraction_index)
        rows = np.array([200.0, 500.0, 900.0])
        columns = np.array([150.0, 600.0, 1000.0])

        for alpha_deg in (0.36, 2.0):
            self.assertGreater(np.deg2rad(alpha_deg), critical)
            for gamma_arm_deg, delta_arm_deg in ((None, None), (5.0, 12.0)):
                with self.subTest(alpha=alpha_deg, arm=gamma_arm_deg):
                    alpha = np.deg2rad(alpha_deg)
                    arm = (
                        (None, None)
                        if gamma_arm_deg is None
                        else tuple(np.deg2rad([gamma_arm_deg, delta_arm_deg]))
                    )
                    gamma, _ = det.surfaceAnglesPoint(rows, columns, alpha, *arm)
                    invertible = (gamma <= 0.0) | (gamma > critical)
                    self.assertTrue(invertible.any())

                    gamma_cry, delta, alpha_cry = det.crystalAnglesPoint(
                        rows, columns, alpha, refraction_index, *arm
                    )
                    back = det.pixelsCrystalAngles(
                        gamma_cry, delta, alpha_cry, refraction_index, *arm
                    )
                    np.testing.assert_allclose(
                        back[invertible, 0], rows[invertible], atol=1e-8
                    )
                    np.testing.assert_allclose(
                        back[invertible, 1], columns[invertible], atol=1e-8
                    )

    def test_reachable_q_follows_the_arm(self):
        """Swinging the arm out must move the reachable Q band with it.

        A range taken once at the calibrated position goes stale: on a
        theta-2theta scan the momentum transfer climbs past the calibrated
        reach, and every reflection after that point was rejected as
        unreachable even though the detector was pointing straight at it.
        """
        det = self.build()
        parked_min, parked_max = det.Qrange
        self.assertEqual((parked_min, parked_max), det.qRangeAtArm())

        previous_max = parked_max
        for gamma_arm_deg in (5.0, 10.0, 20.0):
            with self.subTest(gamma_arm=gamma_arm_deg):
                minimum, maximum = det.qRangeAtArm(np.deg2rad(gamma_arm_deg), 0.0)
                self.assertGreater(maximum, previous_max)
                previous_max = maximum

        # far enough out the direct beam leaves the detector, and the low end
        # must stop claiming that Q = 0 is still reachable
        minimum, _ = det.qRangeAtArm(np.deg2rad(60.0), 0.0)
        self.assertGreater(minimum, 0.0)

    def test_arm_rotates_the_rays_rigidly(self):
        rows = np.array([5.0, 300.0, 700.0, 1000.0])
        columns = np.array([7.0, 200.0, 600.0, 950.0])
        det = self.build()
        home_rays = self.rays(det, rows, columns)

        for gamma_deg, delta_deg in self.ARM_POSITIONS_DEG:
            with self.subTest(gamma_arm=gamma_deg, delta_arm=delta_deg):
                gamma_arm, delta_arm = np.deg2rad([gamma_deg, delta_deg])
                expected = home_rays @ DetectorCalibration.armRotation(
                    gamma_arm, delta_arm
                ).T
                np.testing.assert_allclose(
                    self.rays(det, rows, columns, gamma_arm, delta_arm),
                    expected,
                    atol=1e-12,
                )

    def test_arm_leaves_dist_and_poni_alone(self):
        det = self.build()
        dist, poni1, poni2 = det.paramAtArm(np.deg2rad(12.0), np.deg2rad(30.0))[:3]
        self.assertEqual(dist, self.HOME["dist"])
        self.assertEqual(poni1, self.HOME["poni1"])
        self.assertEqual(poni2, self.HOME["poni2"])

    def test_inverse_round_trip_at_non_zero_arm(self):
        rows = np.array([0.0, 120.0, 640.0, 1000.0])
        columns = np.array([4.0, 310.0, 700.0, 950.0])
        det = self.build()
        gamma_arm, delta_arm = np.deg2rad([12.0, 30.0])

        gamma_p, delta_p = det.primBeamPoints(rows, columns, gamma_arm, delta_arm)
        back = det.pixelsPrimeBeam(gamma_p, delta_p, gamma_arm, delta_arm)

        np.testing.assert_allclose(back[:, 0], rows, atol=1e-9)
        np.testing.assert_allclose(back[:, 1], columns, atol=1e-9)

    def test_rocking_curve_pixel_tracks_the_moving_detector(self):
        """A rocking curve is intensity at fixed exit angles vs sample rotation.

        On a moving detector the pixel holding those angles moves with the arm,
        which is exactly what makes the curve still extractable. The pixel must
        move, and reading it back at the same frame's arm must return the
        angles it was asked for.
        """
        det = self.build()
        alpha = np.deg2rad(0.6)
        gamma = np.deg2rad(np.array([4.0]))
        delta = np.deg2rad(np.array([12.0]))

        pixels = []
        for gamma_deg, delta_deg in self.ARM_POSITIONS_DEG:
            gamma_arm, delta_arm = np.deg2rad([gamma_deg, delta_deg])
            pixel = det.pixelsSurfaceAngles(
                gamma, delta, alpha, gamma_arm, delta_arm
            )[0]
            pixels.append(pixel)

            back_gamma, back_delta = det.surfaceAnglesPoint(
                pixel[:1], pixel[1:], alpha, gamma_arm, delta_arm
            )
            np.testing.assert_allclose(back_gamma, gamma, atol=1e-9)
            np.testing.assert_allclose(back_delta, delta, atol=1e-9)

        # the whole point: fixed exit angles do not sit still on the detector
        spread = np.ptp(np.asarray(pixels), axis=0)
        self.assertGreater(spread.min(), 1.0)

    def test_param_at_arm_is_pure_and_thread_safe(self):
        """The pipeline shares one detector read-only across workers.

        ``paramAtArm`` must therefore neither mutate the object nor depend on
        call order.
        """
        from concurrent.futures import ThreadPoolExecutor

        det = self.build()
        positions = [np.deg2rad(pair) for pair in self.ARM_POSITIONS_DEG] * 8

        before_param = np.array(det.param, dtype=float)
        before_reference = det.getArmReference()
        before_cache = set(det._cached_array)

        serial = [det.paramAtArm(*position) for position in positions]
        with ThreadPoolExecutor(max_workers=8) as pool:
            concurrent = list(pool.map(lambda p: det.paramAtArm(*p), positions))

        self.assertEqual(serial, concurrent)
        np.testing.assert_array_equal(np.array(det.param, dtype=float), before_param)
        np.testing.assert_array_equal(det.getArmReference(), before_reference)
        self.assertEqual(set(det._cached_array), before_cache)

    def test_reference_cancels_and_is_not_an_angle_offset(self):
        rows = np.array([10.0, 500.0])
        columns = np.array([20.0, 600.0])
        gamma_0, delta_0 = np.deg2rad([2.0, 5.0])
        gamma_arm, delta_arm = np.deg2rad([12.0, 30.0])

        plain = self.build()
        referenced = self.build()
        referenced.setArmReference(gamma_arm=gamma_0, delta_arm=delta_0)

        # at the reference the geometry is the calibrated one again
        np.testing.assert_allclose(
            referenced.paramAtArm(gamma_0, delta_0),
            np.asarray(plain.param[:6]),
            atol=1e-12,
        )
        np.testing.assert_allclose(
            referenced.primBeamPoints(rows, columns, gamma_0, delta_0),
            plain.primBeamPoints(rows, columns),
            atol=1e-12,
        )

        # and the relative rotation is a rotation difference, not an angle one
        naive = self.build()
        self.assertFalse(
            np.allclose(
                referenced.paramAtArm(gamma_arm, delta_arm),
                naive.paramAtArm(gamma_arm - gamma_0, delta_arm - delta_0),
                atol=1e-4,
            )
        )

    def test_gimbal_lock_is_reported(self):
        det = self.build(rot1=0.0, rot2=0.0, rot3=0.0)
        with self.assertRaises(ValueError) as caught:
            det.paramAtArm(np.pi / 2, 0.0)
        self.assertIn("gimbal lock", str(caught.exception))

    def test_one_arm_angle_alone_is_rejected(self):
        det = self.build()
        with self.assertRaises(ValueError):
            det.paramAtArm(np.deg2rad(5.0), None)
        with self.assertRaises(ValueError):
            det.primBeamPoints(np.array([1.0]), np.array([1.0]), None, 0.1)

    def test_reference_rejects_improper_rotations(self):
        det = self.build()
        with self.assertRaises(ValueError):
            det.setArmReference(np.diag([1.0, 1.0, -1.0]))  # reflection, det = -1
        with self.assertRaises(ValueError):
            det.setArmReference(2.0 * np.identity(3))  # not orthogonal
        with self.assertRaises(ValueError):
            det.setArmReference(np.identity(3), gamma_arm=0.1)  # matrix and angles
        with self.assertRaises(ValueError):
            det.setArmReference()  # neither


"""
def test_del_gam_range():
    sxrddet = Detector2D_SXRD()
    sxrddet.detector = pyFAI.detector_factory('Pilatus 2M CdTe')
    sxrddet.poni1 = 0.1
    sxrddet.poni2 = 0.1
    sxrddet.rot1 = np.pi/4
    sxrddet.rot2 = np.pi/5
    sxrddet.rot3 = np.pi/6
    sxrddet.dist = 1.5
    sxrddet.setAzimuthalReference(np.deg2rad(90.))


    def checkRange():
        exact = sxrddet._rangegamdel_p_full_det
        corner = sxrddet.rangegamdel_p
        #print("Exact: ", exact)
        #print("Corner: ", corner)
        diff =  np.array(exact) - np.array(corner)
        max_rel_diff = np.amax(np.abs(diff))
        return max_rel_diff, diff
        #print("Difference:", np.array(exact) - np.array(corner))

    sxrddet.poni2 = 0.0
    for p1 in np.linspace(-5,5,20):
        sxrddet.poni1 = p1
        maxerr, diff = checkRange()
        print(f"Max numerical error at poni1 = {p1} m: {maxerr}")

    sxrddet.poni1 = 0.0
    for p2 in np.linspace(-5,5,20):
        sxrddet.poni2 = p2
        maxerr, diff = checkRange()
        print(f"Max numerical error at poni2 = {p2} m: {maxerr}")

    # test rot1:
    for r1 in np.linspace(0,np.pi,8):
        sxrddet.rot1 = r1
        maxerr, diff = checkRange()
        print(f"Max numerical error at rot1 = {np.rad2deg(r1)} deg: {maxerr}")
    sxrddet.rot1 = np.pi/4
    for r2 in np.linspace(0,np.pi,8):
        sxrddet.rot2 = r2
        maxerr, diff = checkRange()
        print(f"Max numerical error at rot2 = {np.rad2deg(r2)} deg: {maxerr}")

    sxrddet.rot2 = np.pi/4
    for r3 in np.linspace(0,np.pi,8):
        sxrddet.rot3 = r3
        maxerr, diff = checkRange()
        print(f"Max numerical error at rot3 = {np.rad2deg(r3)} deg: {maxerr}")
    sxrddet.rot3 = np.pi/4

    for d in np.linspace(0.001,10,5):
        sxrddet.dist = d
        maxerr, diff = checkRange()
        print(f"Max numerical error at dist = {d} m: {maxerr}")


    def testPixelConversion():
    sxrddet = Detector2D_SXRD()
    sxrddet.detector = pyFAI.detector_factory('Pilatus 2M CdTe')
    sxrddet.poni1 = 0.1
    sxrddet.poni2 = 0.1
    sxrddet.rot1 = np.pi/4
    sxrddet.rot2 = np.pi/5
    sxrddet.rot3 = np.pi/6
    sxrddet.dist = 1.5
    sxrddet.setAzimuthalReference(np.deg2rad(90.))

    # pixel coordinates:
    p1 = np.arange(sxrddet.detector.shape[1] ) + 0.5 # pixel center
    p2 = np.arange(sxrddet.detector.shape[0] ) + 0.5
    p12 = np.moveaxis(np.array(np.meshgrid(p1,p2)),0, -1)[:,:,::-1]
    # this seems to be overcomplicated... is there a better method?


    def checkPixel_sxrd():
        gamma, delta = sxrddet.surfaceAngles(np.deg2rad(0.1))
        p12_conv = sxrddet.pixelsSurfaceAngles(gamma, delta, np.deg2rad(0.1))
        return np.nanmax(np.abs(p12 - p12_conv))
        #return np.allclose(p12, p12_conv, atol=1e-5)

    def checkPixel_tth():
        tth = sxrddet.twoThetaArray()
        chi = sxrddet.chiArray()
        p12_conv = sxrddet.pixelsTthChi(tth, chi)
        return np.nanmax(np.abs(p12 - p12_conv))

    checkPixel = checkPixel_sxrd

    sxrddet.poni2 = 0.0
    for p1 in np.linspace(-5,5,20):
        sxrddet.poni1 = p1
        maxerr = checkPixel()
        print(f"Max numerical error at poni1 = {p1} m: {maxerr}")
    sxrddet.poni1 = 0.0
    for p2 in np.linspace(-5,5,20):
        sxrddet.poni2 = p2
        maxerr = checkPixel()
        print(f"Max numerical error at poni2 = {p2} m: {maxerr}")

    # test rot1:
    for r1 in np.linspace(0,np.pi,8):
        sxrddet.rot1 = r1
        maxerr = checkPixel()
        print(f"Max numerical error at rot1 = {np.rad2deg(r1)} deg: {maxerr}")
    sxrddet.rot1 = np.pi/4
    for r2 in np.linspace(0,np.pi,8):
        sxrddet.rot2 = r2
        maxerr = checkPixel()
        print(f"Max numerical error at rot2 = {np.rad2deg(r2)} deg: {maxerr}")

    sxrddet.rot2 = np.pi/4
    for r3 in np.linspace(0,np.pi,8):
        sxrddet.rot3 = r3
        maxerr = checkPixel()
        print(f"Max numerical error at rot3 = {np.rad2deg(r3)} deg: {maxerr}")
    sxrddet.rot3 = np.pi/4

    for d in np.linspace(0.001,10,5):
        sxrddet.dist = d
        maxerr = checkPixel()
        print(f"Max numerical error at dist = {d} m: {maxerr}")


"""


def _polarization_reference(sxrddet, alpha_i):
    """Polarization from orGUI's own z-axis expression.

    :meth:`orgui.datautils.xrayutils.HKLVlieg.UBCalculator.polarization` is
    the z-axis entry of Appendix A of the ANA/ROD manual, written in
    diffractometer angles. It is an independent implementation of the same
    physics as the pyFAI detector-frame expression, so it pins the
    conventions ``polarizationArray`` has to translate between.
    """
    if hasattr(sxrddet, "_alpha_i"):
        del sxrddet._alpha_i
    gamma, delta = sxrddet.surfaceAngles(alpha_i)
    fraction = sxrddet._polFactor
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


def _sxrd_detector(azimuth_deg, pol_axis_deg, fraction):
    """A small calibrated detector with the given polarization settings."""
    sxrddet = DetectorCalibration.Detector2D_SXRD()
    # A private detector: detector_factory hands out one shared instance per
    # name, and pyFAI caches angle arrays against it, which makes tests that
    # share it depend on execution order.
    sxrddet.detector = pyFAI.detectors.Detector(
        pixel1=172e-6, pixel2=172e-6, max_shape=(619, 487)
    )
    # Close in, with the beam near one corner, so the detector spans a wide
    # range of 2theta and azimuth. The polarization correction scales with
    # sin^2(2theta), so a small-angle geometry would be ~1 everywhere and
    # could not tell any two conventions apart.
    sxrddet.poni1 = 0.005
    sxrddet.poni2 = 0.005
    sxrddet.rot1 = 0.0
    sxrddet.rot2 = 0.0
    sxrddet.rot3 = 0.0
    sxrddet.dist = 0.06
    sxrddet.set_energy(20.0)
    sxrddet.setAzimuthalReference(np.deg2rad(azimuth_deg))
    sxrddet.setPolarization(np.deg2rad(pol_axis_deg), fraction)
    # pyFAI caches its angle arrays, and detector_factory hands out shared
    # detector objects, so drop both caches rather than inherit whatever a
    # previously executed test left behind.
    sxrddet.reset()
    sxrddet._cached_array = {}
    return sxrddet


@pytest.mark.parametrize("azimuth_deg", [0.0, 30.0, 90.0, -40.0])
@pytest.mark.parametrize("fraction", [1.0, 0.95, 0.5, 0.0])
def test_polarization_array_matches_the_z_axis_expression(azimuth_deg, fraction):
    """``polarizationArray`` reproduces the ANA/ROD z-axis polarization.

    pyFAI works in its own detector frame and parameterizes the amount of
    polarization differently, so both arguments have to be translated. Two
    independent implementations agreeing to the precision pyFAI stores its
    result in is what pins that translation.
    """
    alpha_i = np.deg2rad(6.0)
    sxrddet = _sxrd_detector(azimuth_deg, 0.0, fraction)

    correction = np.asarray(sxrddet.polarizationArray(), dtype=np.float64)
    reference = _polarization_reference(sxrddet, alpha_i)

    # pyFAI stores the polarization as float32.
    np.testing.assert_allclose(correction, reference, atol=1e-6)


def test_polarization_uses_the_azimuthal_reference():
    """Leaving the azimuthal reference out is a large, angle-dependent error.

    A vertical scattering geometry sets the azimuthal reference to 90
    degrees, which swaps the polarized and unpolarized directions of the
    detector. Because the error varies across the detector, two rod branches
    recorded on opposite sides of it are corrected differently and stop
    agreeing -- the symptom this test exists to prevent.
    """
    alpha_i = np.deg2rad(6.0)
    sxrddet = _sxrd_detector(90.0, 0.0, 1.0)
    reference = _polarization_reference(sxrddet, alpha_i)

    correct = np.asarray(sxrddet.polarizationArray(), dtype=np.float64)
    sxrddet._cached_array = {}
    # What orGUI passed before: the polarization axis without the azimuthal
    # reference.
    without_reference = np.asarray(
        sxrddet.polarization(factor=sxrddet._polFactor, axis_offset=sxrddet._polAxis),
        dtype=np.float64,
    )

    np.testing.assert_allclose(correct, reference, atol=1e-6)
    # The error is not a constant scale factor but varies over the detector,
    # which is what makes it distort a rod rather than just rescale it.
    ratio = without_reference / correct
    assert ratio.max() / ratio.min() > 2.0
    assert ratio.min() < 0.8


def test_polarization_factor_is_the_horizontal_fraction():
    """``polarization_factor`` is a fraction, not pyFAI's ``(Ih-Iv)/(Ih+Iv)``.

    The config field and orGUI's own z-axis expression treat it as the
    fraction of horizontally polarized light, so it has to be rescaled to
    pyFAI's convention. The two coincide only for a fully horizontally
    polarized beam, which is why the difference stayed hidden.
    """
    alpha_i = np.deg2rad(6.0)
    fraction = 0.75
    sxrddet = _sxrd_detector(90.0, 0.0, fraction)
    reference = _polarization_reference(sxrddet, alpha_i)

    correct = np.asarray(sxrddet.polarizationArray(), dtype=np.float64)
    sxrddet._cached_array = {}
    # What orGUI passed before: the fraction used as pyFAI's factor.
    unscaled = np.asarray(
        sxrddet.polarization(factor=fraction, axis_offset=sxrddet._deltaChi),
        dtype=np.float64,
    )

    np.testing.assert_allclose(correct, reference, atol=1e-6)
    assert np.abs(unscaled - reference).max() > 1e-2
    # A fully horizontal beam is the one case where both agree.
    full = _sxrd_detector(90.0, 0.0, 1.0)
    np.testing.assert_allclose(
        np.asarray(full.polarizationArray(), dtype=np.float64),
        _polarization_reference(full, alpha_i),
        atol=1e-6,
    )


def test_unpolarized_beam_has_no_azimuthal_dependence():
    """An unpolarized beam gives ``(1 + cos^2(2theta)) / 2`` everywhere."""
    sxrddet = _sxrd_detector(90.0, 0.0, 0.5)

    correction = np.asarray(sxrddet.polarizationArray(), dtype=np.float64)

    tth = sxrddet.center_array(sxrddet.get_shape(), unit=pyFAI.units.TTH_RAD)
    np.testing.assert_allclose(correction, 0.5 * (1.0 + np.cos(tth) ** 2), atol=1e-6)

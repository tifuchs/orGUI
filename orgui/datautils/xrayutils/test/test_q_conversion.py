"""Equivalence of orGUI's QAlpha conversion with pyFAI's FiberIntegrator.

orGUI describes reciprocal space through the surface angles ``gamma`` and
``delta`` (:meth:`DetectorCalibration.Detector2D_SXRD.surfaceAngles`) combined
with the Vlieg diffraction equation
(:meth:`HKLVlieg.VliegAngles.QAlpha`). The "Q-plot" action of the GUI instead
rebins the displayed image onto a regular grid of in-plane and out-of-plane
momentum transfer using pyFAI's ``FiberIntegrator``.

Both routes have to describe the same reciprocal space. These tests pin that
down numerically: first per pixel, which is an exact identity, and then through
the full ``integrate2d_grazing_incidence`` call actually used by the GUI, where
the agreement is limited by the width of the rebinning grid.
"""

import unittest

import numpy as np
import pyFAI

from .. import DetectorCalibration, HKLVlieg

try:
    from pyFAI import units as pyFAI_units
    from pyFAI.integrator.fiber import FiberIntegrator

    fiber_available = True
except ImportError:  # pyFAI < 2025.1
    fiber_available = False

requires_fiber = unittest.skipUnless(
    fiber_available, "pyFAI >= 2025.1 with FiberIntegrator is required"
)

ENERGY = 78.0  # keV

# Azimuthal reference (in degrees) -> pyFAI EXIF sample orientation, as
# implemented by orgui.app.orGUI.orGUI._qConversionSampleOrientation.
#
# pyFAI pairs every orientation with an in-plane mirrored partner --
# (6, 7), (3, 4), (5, 8) and (1, 2). Partners share q_oop and differ only in
# the sign of q_ip, see TestSampleOrientationMirrorPairs. QAlpha returns the
# in-plane momentum transfer as the unsigned radial component
# sqrt(Qx**2 + Qy**2), so it cannot distinguish the partners and the tests
# below compare |q_ip|.
AZIMUTH_ORIENTATION = ((0.0, 7), (90.0, 4), (180.0, 8), (270.0, 1))

FLAT_DETECTOR = (0.0, 0.0, 0.0)
TILTED_DETECTOR = (np.deg2rad(1.7), np.deg2rad(-0.9), np.deg2rad(3.4))


def _geometry(azimuth_deg, rotations=FLAT_DETECTOR):
    """Return a Detector2D_SXRD with a fully specified geometry."""
    det = DetectorCalibration.Detector2D_SXRD()
    det.detector = pyFAI.detector_factory("Pilatus2m")
    det.wavelength = (12.39842 / ENERGY) * 1e-10
    det.dist = 0.729
    det.poni1 = 0.21
    det.poni2 = 0.14
    det.rot1, det.rot2, det.rot3 = rotations
    det.setAzimuthalReference(np.deg2rad(azimuth_deg))
    return det


def _vlieg_angles():
    lattice = HKLVlieg.Lattice(np.array([4.0, 4.0, 4.0]), np.array([90.0, 90.0, 90.0]))
    return HKLVlieg.VliegAngles(HKLVlieg.UBCalculator(lattice, ENERGY))


def _pixel_grid(det, stride=67):
    """Flattened pixel coordinates covering the whole detector."""
    shape = det.detector.shape
    rows, cols = np.meshgrid(
        np.arange(0, shape[0], stride, dtype=np.float64),
        np.arange(0, shape[1], stride, dtype=np.float64),
        indexing="ij",
    )
    return rows.ravel(), cols.ravel()


def _q_from_qalpha(det, angles, alpha_i, rows, cols):
    """In-plane and out-of-plane momentum transfer via surfaceAngles + QAlpha."""
    gamma, delta = det.surfaceAnglesPoint(rows, cols, alpha_i)
    Qxyz = angles.QAlpha(alpha_i, delta, gamma)
    q_ip = np.hypot(Qxyz[..., 0], Qxyz[..., 1])
    q_oop = Qxyz[..., 2]
    return q_ip, q_oop


def _q_from_pyfai(det, alpha_i, orientation, rows, cols):
    """Per-pixel q_ip and q_oop from the equations used by FiberIntegrator."""
    param = np.array(
        [det.dist, det.poni1, det.poni2, det.rot1, det.rot2, det.rot3],
        dtype=np.float64,
    )
    # calc_pos_zyx returns (along beam, slow/vertical, fast/horizontal), while
    # the pyFAI unit equations expect x=horizontal, y=vertical, z=beam.
    t3, t1, t2 = det.calc_pos_zyx(d1=rows, d2=cols, param=param, do_parallax=True)
    kwargs = {
        "x": t2,
        "y": t1,
        "z": t3,
        "wavelength": det.wavelength,
        "incident_angle": alpha_i,
        "tilt_angle": 0.0,
        "sample_orientation": orientation,
    }
    # pyFAI returns inverse nanometre, orGUI works in inverse Angstrom.
    return pyFAI_units.eq_qip(**kwargs) / 10.0, pyFAI_units.eq_qoop(**kwargs) / 10.0


@requires_fiber
class TestQAlphaMatchesPyFAIPerPixel(unittest.TestCase):
    """QAlpha and pyFAI agree pixel by pixel, before any rebinning."""

    def test_matches_for_every_gui_orientation(self):
        angles = _vlieg_angles()
        for azimuth_deg, orientation in AZIMUTH_ORIENTATION:
            for rotations in (FLAT_DETECTOR, TILTED_DETECTOR):
                for alpha_deg in (0.0, 0.08, 2.5):
                    with self.subTest(
                        azimuth=azimuth_deg,
                        orientation=orientation,
                        tilted=rotations is TILTED_DETECTOR,
                        alpha_i=alpha_deg,
                    ):
                        det = _geometry(azimuth_deg, rotations)
                        rows, cols = _pixel_grid(det)
                        alpha_i = np.deg2rad(alpha_deg)

                        q_ip, q_oop = _q_from_qalpha(det, angles, alpha_i, rows, cols)
                        ref_ip, ref_oop = _q_from_pyfai(
                            det, alpha_i, orientation, rows, cols
                        )

                        np.testing.assert_allclose(q_oop, ref_oop, atol=1e-9)
                        np.testing.assert_allclose(q_ip, np.abs(ref_ip), atol=1e-9)

    def test_total_momentum_transfer_is_orientation_independent(self):
        """|Q| is a physical invariant, so it must not depend on the labeling."""
        angles = _vlieg_angles()
        det = _geometry(90.0)
        rows, cols = _pixel_grid(det)
        alpha_i = np.deg2rad(0.6)

        q_ip, q_oop = _q_from_qalpha(det, angles, alpha_i, rows, cols)
        magnitude = np.hypot(q_ip, q_oop)

        for orientation in range(1, 9):
            with self.subTest(orientation=orientation):
                ref_ip, ref_oop = _q_from_pyfai(det, alpha_i, orientation, rows, cols)
                np.testing.assert_allclose(
                    magnitude, np.hypot(ref_ip, ref_oop), atol=1e-9
                )


@requires_fiber
class TestSampleOrientationMirrorPairs(unittest.TestCase):
    """Document why the comparisons above use |q_ip|."""

    def test_paired_orientations_share_qoop_and_mirror_qip(self):
        det = _geometry(90.0)
        rows, cols = _pixel_grid(det, stride=211)
        alpha_i = np.deg2rad(0.6)

        for first, second in ((6, 7), (3, 4), (5, 8), (1, 2)):
            with self.subTest(pair=(first, second)):
                ip_a, oop_a = _q_from_pyfai(det, alpha_i, first, rows, cols)
                ip_b, oop_b = _q_from_pyfai(det, alpha_i, second, rows, cols)

                np.testing.assert_allclose(oop_a, oop_b, atol=1e-12)
                np.testing.assert_allclose(ip_a, -ip_b, atol=1e-12)


@requires_fiber
class TestQAlphaMatchesFiberIntegrator(unittest.TestCase):
    """The rebinned FiberIntegrator map places intensity where QAlpha predicts."""

    def test_single_pixel_peak_lands_at_the_qalpha_position(self):
        angles = _vlieg_angles()
        alpha_i = np.deg2rad(0.6)

        for azimuth_deg, orientation in AZIMUTH_ORIENTATION:
            det = _geometry(azimuth_deg)
            integrator = FiberIntegrator(
                dist=det.dist,
                poni1=det.poni1,
                poni2=det.poni2,
                wavelength=det.wavelength,
                rot1=det.rot1,
                rot2=det.rot2,
                rot3=det.rot3,
                detector=det.detector,
            )

            for row, col in ((900, 700), (600, 1000)):
                with self.subTest(azimuth=azimuth_deg, pixel=(row, col)):
                    image = np.zeros(det.detector.shape, dtype=np.float64)
                    image[row, col] = 1000.0

                    # identical call to the one performed by the Q-plot action
                    result = integrator.integrate2d_grazing_incidence(
                        image,
                        sample_orientation=orientation,
                        incident_angle=alpha_i,
                        tilt_angle=0,
                        unit_oop="qoop_A^-1",
                        unit_ip="qip_A^-1",
                    )

                    peak = np.unravel_index(
                        np.argmax(result.intensity), result.intensity.shape
                    )
                    q_ip_binned = result.inplane[peak[1]]
                    q_oop_binned = result.outofplane[peak[0]]

                    q_ip, q_oop = _q_from_qalpha(
                        det,
                        angles,
                        alpha_i,
                        np.array([float(row)]),
                        np.array([float(col)]),
                    )

                    # the rebinning quantises the position to one grid cell
                    bin_ip = (
                        result.inplane.max() - result.inplane.min()
                    ) / result.inplane.size
                    bin_oop = (
                        result.outofplane.max() - result.outofplane.min()
                    ) / result.outofplane.size

                    self.assertLess(abs(abs(q_ip_binned) - q_ip[0]), bin_ip)
                    self.assertLess(abs(q_oop_binned - q_oop[0]), bin_oop)


if __name__ == "__main__":
    unittest.main()

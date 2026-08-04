"""Equivalence of the Q-plot conversion with orGUI's own QAlpha calculation.

orGUI describes reciprocal space through the surface angles ``gamma`` and
``delta`` (:meth:`DetectorCalibration.Detector2D_SXRD.surfaceAngles`) combined
with the Vlieg diffraction equation (:meth:`HKLVlieg.VliegAngles.QAlpha`). The
experimental Q-plot instead rebins a whole image onto a regular grid using
pyFAI's ``FiberIntegrator``, driven by :mod:`orgui.app.qconversion`.

The azimuthal reference is deliberately tested at values that are *not*
multiples of 90 degrees. pyFAI's ``sample_orientation`` flag can only express
quarter turns, so those cases are exactly the ones an implementation based on
that flag alone gets wrong.
"""

import unittest

import numpy as np
import pyFAI

from orgui.app import qconversion
from orgui.datautils.xrayutils import DetectorCalibration, HKLVlieg

requires_fiber = unittest.skipUnless(
    qconversion.HAS_FIBER_INTEGRATOR,
    "pyFAI >= 2025.1 with FiberIntegrator is required",
)

ENERGY = 78.0  # keV

FLAT_DETECTOR = (0.0, 0.0, 0.0)
TILTED_DETECTOR = (np.deg2rad(1.7), np.deg2rad(-0.9), np.deg2rad(3.4))

# not restricted to multiples of 90 degrees on purpose
AZIMUTHS = (0.0, 17.3, 45.0, 90.0, 135.0, 233.7, 341.9)
INCIDENCE_ANGLES = (0.0, 0.08, 0.6, 2.5)


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


def _ub_calculator():
    lattice = HKLVlieg.Lattice(np.array([4.0, 4.0, 4.0]), np.array([90.0, 90.0, 90.0]))
    ub = HKLVlieg.UBCalculator(lattice, ENERGY)
    ub.defaultU()
    return ub


def _pixel_grid(det, stride=251):
    """Flattened pixel coordinates covering the whole detector."""
    shape = det.detector.shape
    rows, cols = np.meshgrid(
        np.arange(0, shape[0], stride, dtype=np.float64),
        np.arange(0, shape[1], stride, dtype=np.float64),
        indexing="ij",
    )
    return rows.ravel(), cols.ravel()


def _pixel_positions(det, rows, cols):
    """pyFAI pixel positions as ``(x, y, z)`` = (horizontal, vertical, beam)."""
    param = np.array(
        [det.dist, det.poni1, det.poni2, det.rot1, det.rot2, det.rot3],
        dtype=np.float64,
    )
    t3, t1, t2 = det.calc_pos_zyx(d1=rows, d2=cols, param=param, do_parallax=True)
    return t2, t1, t3


def _q_from_conversion(det, rows, cols, alpha_i, frame="Q_alpha", **frame_angles):
    """Momentum transfer in inverse Angstrom from :mod:`qconversion`."""
    x, y, z = _pixel_positions(det, rows, cols)
    components = qconversion.qComponents(
        frame,
        x,
        y,
        z,
        det.wavelength,
        incident_angle=alpha_i,
        tilt_angle=qconversion.tiltAngleFromAzimuth(det.getAzimuthalReference()),
        sample_orientation=qconversion.SAMPLE_ORIENTATION,
        **frame_angles,
    )
    return np.stack([component / 10.0 for component in components], axis=-1)


@requires_fiber
class TestQAlphaEquivalence(unittest.TestCase):
    """The alpha frame must reproduce QAlpha component by component."""

    def test_matches_qalpha_for_arbitrary_azimuth(self):
        angles = HKLVlieg.VliegAngles(_ub_calculator())
        for azimuth_deg in AZIMUTHS:
            for rotations in (FLAT_DETECTOR, TILTED_DETECTOR):
                for alpha_deg in INCIDENCE_ANGLES:
                    with self.subTest(
                        azimuth=azimuth_deg,
                        tilted=rotations is TILTED_DETECTOR,
                        alpha_i=alpha_deg,
                    ):
                        det = _geometry(azimuth_deg, rotations)
                        rows, cols = _pixel_grid(det)
                        alpha_i = np.deg2rad(alpha_deg)

                        actual = _q_from_conversion(det, rows, cols, alpha_i)
                        gamma, delta = det.surfaceAnglesPoint(rows, cols, alpha_i)
                        expected = angles.QAlpha(alpha_i, delta, gamma)

                        np.testing.assert_allclose(actual, expected, atol=1e-9)

    def test_builtin_pyfai_units_cannot_express_a_general_azimuth(self):
        """Motivation for the custom units, see :mod:`qconversion`.

        ``sample_orientation`` only spans quarter turns, so no value of it can
        reproduce an azimuthal reference of 45 degrees.
        """
        from pyFAI import units as pyFAI_units

        angles = HKLVlieg.VliegAngles(_ub_calculator())
        det = _geometry(45.0)
        rows, cols = _pixel_grid(det)
        alpha_i = np.deg2rad(0.6)

        gamma, delta = det.surfaceAnglesPoint(rows, cols, alpha_i)
        expected = angles.QAlpha(alpha_i, delta, gamma)
        expected_oop = expected[..., 2]
        x, y, z = _pixel_positions(det, rows, cols)

        best = np.inf
        for orientation in range(1, 9):
            actual_oop = (
                pyFAI_units.eq_qoop(
                    x=x,
                    y=y,
                    z=z,
                    wavelength=det.wavelength,
                    incident_angle=alpha_i,
                    tilt_angle=0.0,
                    sample_orientation=orientation,
                )
                / 10.0
            )
            best = min(best, np.max(np.abs(expected_oop - actual_oop)))

        # the closest flag is still off by inverse Angstrom, not by rounding
        self.assertGreater(best, 1.0)


@requires_fiber
class TestFrames(unittest.TestCase):
    """The frame rotations must follow the Vlieg convention."""

    def test_q_phi_reproduces_hkl(self):
        ub = _ub_calculator()
        angles = HKLVlieg.VliegAngles(ub)
        det = _geometry(45.0)
        rows, cols = _pixel_grid(det, stride=401)
        alpha_i = np.deg2rad(0.6)
        omega, chi, phi = np.deg2rad(23.0), np.deg2rad(4.0), np.deg2rad(-7.0)

        q_phi = _q_from_conversion(
            det, rows, cols, alpha_i, frame="Q_phi", omega=omega, chi=chi, phi=phi
        )
        actual = np.einsum("ij,...j->...i", np.linalg.inv(ub.getUB()), q_phi)

        gamma, delta = det.surfaceAnglesPoint(rows, cols, alpha_i)
        h, k, ll = angles.anglesToHkl(alpha_i, delta, gamma, omega, chi, phi)

        np.testing.assert_allclose(actual[..., 0], h, atol=1e-9)
        np.testing.assert_allclose(actual[..., 1], k, atol=1e-9)
        np.testing.assert_allclose(actual[..., 2], ll, atol=1e-9)

    def test_q_cryst_is_b_times_hkl(self):
        ub = _ub_calculator()
        angles = HKLVlieg.VliegAngles(ub)
        det = _geometry(45.0)
        rows, cols = _pixel_grid(det, stride=401)
        alpha_i = np.deg2rad(0.6)
        omega, chi, phi = np.deg2rad(23.0), np.deg2rad(4.0), np.deg2rad(-7.0)

        actual = _q_from_conversion(
            det,
            rows,
            cols,
            alpha_i,
            frame="Q_cryst",
            omega=omega,
            chi=chi,
            phi=phi,
            U=ub.getU(),
        )

        gamma, delta = det.surfaceAnglesPoint(rows, cols, alpha_i)
        hkl = np.stack(
            angles.anglesToHkl(alpha_i, delta, gamma, omega, chi, phi), axis=-1
        )
        expected = np.einsum("ij,...j->...i", ub.lattice.B_mat, hkl)

        np.testing.assert_allclose(actual, expected, atol=1e-9)

    def test_frames_requiring_u_are_rejected_without_it(self):
        for frame in qconversion.FRAMES_REQUIRING_U:
            with self.subTest(frame=frame):
                with self.assertRaises(ValueError):
                    qconversion.frameMatrix(frame, alpha=0.1, omega=0.2)
                with self.assertRaises(ValueError):
                    qconversion.fiberUnits(frame)

    def test_every_frame_preserves_the_momentum_transfer(self):
        det = _geometry(45.0)
        rows, cols = _pixel_grid(det, stride=401)
        alpha_i = np.deg2rad(0.6)
        frame_angles = {
            "omega": np.deg2rad(23.0),
            "chi": np.deg2rad(4.0),
            "phi": np.deg2rad(-7.0),
            "U": _ub_calculator().getU(),
        }

        reference = np.linalg.norm(
            _q_from_conversion(det, rows, cols, alpha_i), axis=-1
        )
        for frame in qconversion.FRAMES:
            with self.subTest(frame=frame):
                q = _q_from_conversion(
                    det, rows, cols, alpha_i, frame=frame, **frame_angles
                )
                np.testing.assert_allclose(
                    np.linalg.norm(q, axis=-1), reference, atol=1e-9
                )

    def test_alpha_frame_matrix_is_the_identity(self):
        np.testing.assert_allclose(
            qconversion.frameMatrix("Q_alpha", alpha=0.3), np.identity(3), atol=1e-15
        )

    def test_unknown_frame_is_rejected(self):
        with self.assertRaises(ValueError):
            qconversion.frameMatrix("Q_nonsense")
        with self.assertRaises(ValueError):
            qconversion.fiberUnits("Q_nonsense")

    def test_frame_metadata_is_consistent(self):
        self.assertTrue(
            set(qconversion.FRAMES_REQUIRING_OMEGA) <= set(qconversion.FRAMES)
        )
        self.assertTrue(set(qconversion.FRAMES_REQUIRING_U) <= set(qconversion.FRAMES))
        # undoing U also requires the sample rotations to be undone first
        self.assertTrue(
            set(qconversion.FRAMES_REQUIRING_U)
            <= set(qconversion.FRAMES_REQUIRING_OMEGA)
        )
        self.assertNotIn(qconversion.DEFAULT_FRAME, qconversion.FRAMES_REQUIRING_OMEGA)
        self.assertEqual(set(qconversion.FRAME_LABELS), set(qconversion.FRAMES))


class TestInverseConversion(unittest.TestCase):
    """The Q-plot readout inverts the projection, so it must round-trip.

    The forward conversion keeps the out-of-plane component and collapses the
    two in-plane components into their radial distance. The inverse recovers
    the discarded direction from the Ewald condition, so a coordinate taken
    from the forward conversion has to lead back to the same detector angles.

    These tests do not need pyFAI's rebinning, only the geometry.
    """

    SAMPLE_ANGLES = dict(
        omega=np.deg2rad(23.0), chi=np.deg2rad(4.0), phi=np.deg2rad(-7.0)
    )

    def _forward(self, frame, alpha_i, delta, gamma, ub):
        """``(Q, q_ip, q_oop)`` of one detector angle pair, in the frame."""
        angles = HKLVlieg.VliegAngles(ub)
        q_alpha = angles.QAlpha(
            np.asarray(alpha_i), np.asarray(delta), np.asarray(gamma)
        )
        rotation = qconversion.frameMatrix(
            frame, alpha=alpha_i, U=ub.getU(), **self.SAMPLE_ANGLES
        )
        Q = rotation @ q_alpha
        # the signing convention of qIpOop
        return Q, np.copysign(np.hypot(Q[0], Q[1]), Q[0]), Q[2]

    def test_round_trip_recovers_the_detector_angles(self):
        ub = _ub_calculator()
        k = ub.getK()
        for frame in qconversion.FRAMES:
            for alpha_deg in INCIDENCE_ANGLES:
                for delta_deg, gamma_deg in ((12.0, 4.0), (-31.5, 22.0), (48.0, -3.0)):
                    with self.subTest(
                        frame=frame, alpha_i=alpha_deg, delta=delta_deg, gamma=gamma_deg
                    ):
                        alpha_i = np.deg2rad(alpha_deg)
                        delta, gamma = np.deg2rad([delta_deg, gamma_deg])
                        Q, q_ip, q_oop = self._forward(frame, alpha_i, delta, gamma, ub)

                        solutions = qconversion.qFromIpOop(
                            frame,
                            q_ip,
                            q_oop,
                            k,
                            alpha=alpha_i,
                            U=ub.getU(),
                            **self.SAMPLE_ANGLES,
                        )

                        self.assertTrue(solutions)
                        # the true momentum transfer is among the solutions
                        best = min(solutions, key=lambda s: np.linalg.norm(s - Q))
                        np.testing.assert_allclose(best, Q, atol=1e-9)

                        recovered = qconversion.detectorAngles(
                            best,
                            k,
                            frame=frame,
                            alpha=alpha_i,
                            U=ub.getU(),
                            **self.SAMPLE_ANGLES,
                        )
                        np.testing.assert_allclose(recovered, (delta, gamma), atol=1e-9)

    def test_the_alpha_frame_is_never_ambiguous(self):
        """The beam has no x component there, so the sign fixes the solution.

        This is why the default frame needs no disambiguation at all.
        """
        ub = _ub_calculator()
        rng = np.random.default_rng(20260803)
        for _ in range(200):
            alpha_i = np.deg2rad(rng.uniform(0.0, 3.0))
            delta = np.deg2rad(rng.uniform(-70.0, 70.0))
            gamma = np.deg2rad(rng.uniform(-20.0, 60.0))
            _, q_ip, q_oop = self._forward("Q_alpha", alpha_i, delta, gamma, ub)

            solutions = qconversion.qFromIpOop(
                "Q_alpha", q_ip, q_oop, ub.getK(), alpha=alpha_i
            )

            self.assertEqual(len(solutions), 1)

    def test_solutions_keep_the_sign_of_the_in_plane_coordinate(self):
        ub = _ub_calculator()
        for frame in qconversion.FRAMES:
            for delta_deg in (-40.0, -5.0, 5.0, 40.0):
                with self.subTest(frame=frame, delta=delta_deg):
                    alpha_i = np.deg2rad(0.6)
                    delta = np.deg2rad(delta_deg)
                    gamma = np.deg2rad(14.0)
                    _, q_ip, q_oop = self._forward(frame, alpha_i, delta, gamma, ub)

                    solutions = qconversion.qFromIpOop(
                        frame,
                        q_ip,
                        q_oop,
                        ub.getK(),
                        alpha=alpha_i,
                        U=ub.getU(),
                        **self.SAMPLE_ANGLES,
                    )

                    self.assertTrue(solutions)
                    for Q in solutions:
                        self.assertEqual(np.sign(Q[0]), np.sign(q_ip))

    def test_unreachable_coordinates_have_no_solution(self):
        """Beyond the Ewald sphere nothing can be reported."""
        ub = _ub_calculator()

        solutions = qconversion.qFromIpOop(
            "Q_alpha", 10.0 * ub.getK(), 0.0, ub.getK(), alpha=np.deg2rad(0.6)
        )

        self.assertEqual(solutions, ())

    def test_non_finite_coordinates_have_no_solution(self):
        ub = _ub_calculator()
        for q_ip, q_oop in ((np.nan, 0.0), (0.0, np.inf)):
            with self.subTest(q_ip=q_ip, q_oop=q_oop):
                self.assertEqual(
                    qconversion.qFromIpOop("Q_alpha", q_ip, q_oop, ub.getK()), ()
                )

    def test_beam_direction_matches_the_collapsed_conversion(self):
        """``beamDirection`` must be the ``c`` of the forward conversion."""
        for frame in qconversion.FRAMES:
            for alpha_deg in INCIDENCE_ANGLES:
                with self.subTest(frame=frame, alpha_i=alpha_deg):
                    ub = _ub_calculator()
                    alpha_i = np.deg2rad(alpha_deg)
                    _, c = qconversion.conversionCoefficients(
                        frame,
                        incident_angle=alpha_i,
                        tilt_angle=0.83,
                        U=ub.getU(),
                        **self.SAMPLE_ANGLES,
                    )

                    beam = qconversion.beamDirection(
                        frame, alpha=alpha_i, U=ub.getU(), **self.SAMPLE_ANGLES
                    )

                    np.testing.assert_allclose(beam, c, atol=1e-12)


@requires_fiber
class TestIntegrateImage(unittest.TestCase):
    """The rebinned map places intensity where QAlpha predicts."""

    def test_single_pixel_peak_lands_at_the_qalpha_position(self):
        angles = HKLVlieg.VliegAngles(_ub_calculator())
        alpha_i = np.deg2rad(0.6)

        for azimuth_deg in (45.0, 233.7):
            det = _geometry(azimuth_deg)
            for row, col in ((900, 700), (600, 1000)):
                with self.subTest(azimuth=azimuth_deg, pixel=(row, col)):
                    image = np.zeros(det.detector.shape, dtype=np.float64)
                    image[row, col] = 1000.0

                    result = qconversion.integrateImage(det, image, alpha_i)

                    peak = np.unravel_index(
                        np.argmax(result.intensity), result.intensity.shape
                    )
                    q_ip_binned = result.inplane[peak[1]]
                    q_oop_binned = result.outofplane[peak[0]]

                    gamma, delta = det.surfaceAnglesPoint(
                        np.array([float(row)]), np.array([float(col)]), alpha_i
                    )
                    expected = angles.QAlpha(alpha_i, delta, gamma)
                    expected_ip = np.hypot(expected[..., 0], expected[..., 1])[0]
                    expected_oop = expected[..., 2][0]

                    # the rebinning quantises the position to one grid cell
                    bin_ip = (
                        result.inplane.max() - result.inplane.min()
                    ) / result.inplane.size
                    bin_oop = (
                        result.outofplane.max() - result.outofplane.min()
                    ) / result.outofplane.size

                    self.assertLess(abs(abs(q_ip_binned) - expected_ip), bin_ip)
                    self.assertLess(abs(q_oop_binned - expected_oop), bin_oop)


@requires_fiber
class TestIntegratorCache(unittest.TestCase):
    """The pyFAI integrator is reused, but never at the cost of stale results.

    A fresh integrator resets on every call, which is slow and floods the log.
    Reusing one is only safe while the conversion is unchanged, because pyFAI
    caches the per-pixel unit arrays under the unit name.
    """

    def setUp(self):
        self.detector = _geometry(45.0)
        self.image = np.zeros(self.detector.detector.shape, dtype=np.float64)
        self.image[900, 700] = 1000.0
        self.alpha = np.deg2rad(0.6)
        qconversion._INTEGRATOR.clear()
        self.addCleanup(qconversion._INTEGRATOR.clear)
        self.addCleanup(qconversion._CONVERSION.clear)

    def _integrate(self, **keywords):
        qconversion._CONVERSION.clear()
        return qconversion.integrateImage(
            self.detector, self.image, self.alpha, npt=200, **keywords
        )

    def test_integrator_is_reused_when_nothing_changed(self):
        self._integrate()
        first = qconversion._INTEGRATOR._integrator
        self._integrate()

        self.assertIs(qconversion._INTEGRATOR._integrator, first)

    def test_integrator_is_rebuilt_when_the_conversion_changes(self):
        self._integrate(frame="Q_alpha")
        first = qconversion._INTEGRATOR._integrator
        self._integrate(frame="Q_chi", omega=np.deg2rad(10.0), chi=np.deg2rad(5.0))

        self.assertIsNot(qconversion._INTEGRATOR._integrator, first)

    def test_reuse_does_not_serve_stale_coordinates(self):
        """chi tilts the surface normal, so the axes really have to change."""

        def peak(chi_deg):
            result = self._integrate(
                frame="Q_chi", omega=np.deg2rad(10.0), chi=np.deg2rad(chi_deg)
            )
            index = np.unravel_index(
                np.argmax(result.intensity), result.intensity.shape
            )
            return result.inplane[index[1]], result.outofplane[index[0]]

        straight = peak(0.0)
        tilted = peak(9.0)
        again = peak(0.0)

        self.assertNotAlmostEqual(straight[1], tilted[1], places=3)
        np.testing.assert_allclose(straight, again, atol=1e-9)


class TestKernelBackends(unittest.TestCase):
    """The compiled kernel and the numpy fallback must agree."""

    def setUp(self):
        det = _geometry(45.0)
        rows, cols = _pixel_grid(det, stride=53)
        self.arguments = (
            "Q_alpha",
            *_pixel_positions(det, rows, cols),
            det.wavelength,
        )
        self.keywords = {
            "incident_angle": np.deg2rad(0.6),
            "tilt_angle": qconversion.tiltAngleFromAzimuth(np.deg2rad(45.0)),
        }
        self.addCleanup(qconversion._CONVERSION.clear)

    def _q_ip_oop(self):
        qconversion._CONVERSION.clear()
        return qconversion.qIpOop(*self.arguments, **self.keywords)

    def test_q_ip_oop_matches_the_reference_components(self):
        q_ip, q_oop = self._q_ip_oop()
        qx, qy, qz = qconversion.qComponents(*self.arguments, **self.keywords)

        np.testing.assert_allclose(q_ip, np.copysign(np.hypot(qx, qy), qx), atol=1e-9)
        np.testing.assert_allclose(q_oop, qz, atol=1e-9)

    def test_compiled_kernel_matches_the_numpy_fallback(self):
        if not qconversion.hasKernel():
            self.skipTest("the compiled conversion kernel is not available")

        compiled_ip, compiled_oop = self._q_ip_oop()

        kernel = qconversion._KERNEL
        qconversion._KERNEL = None  # force the numpy fallback
        try:
            fallback_ip, fallback_oop = self._q_ip_oop()
        finally:
            qconversion._KERNEL = kernel

        np.testing.assert_allclose(compiled_ip, fallback_ip, atol=1e-9)
        np.testing.assert_allclose(compiled_oop, fallback_oop, atol=1e-9)

    def test_repeated_call_is_served_from_the_cache(self):
        first = self._q_ip_oop()
        second = qconversion.qIpOop(*self.arguments, **self.keywords)
        # pyFAI evaluates the two units separately; the image is converted once
        self.assertIs(first[0], second[0])
        self.assertIs(first[1], second[1])


if __name__ == "__main__":
    unittest.main()

"""Regression tests for optional CTR resolution modeling."""

import unittest

import numpy as np

from .. import CTRplotutil, CTRresolution


def _angles(gamma):
    gamma = np.asarray(gamma, dtype=np.float64)
    zeros = np.zeros_like(gamma)
    return np.rec.fromarrays(
        [zeros, zeros, gamma, zeros, zeros, zeros],
        names="alpha,delta,gamma,omega,chi,phi",
    )


class QuadraticCrystal:
    """Crystal-like test double with intensity equal to L squared."""

    @staticmethod
    def F(h, k, l):  # noqa: N802,E741
        return np.asarray(l, dtype=np.complex128)


class TestResolutionFunctions(unittest.TestCase):
    def test_constant_and_gamma_dependent_widths(self):
        h = np.array([0.0, 1.0, 2.0])
        l = np.array([0.1, 0.2, 0.3])  # noqa: E741
        gamma = np.array([0.0, np.pi / 6.0, np.pi / 2.0])

        constant = CTRresolution.BoxResolution(0.2)
        varying = CTRresolution.GaussianResolution(0.2, 0.4)

        np.testing.assert_allclose(constant.width(h, h, l), 0.2)
        np.testing.assert_allclose(
            varying.width(h, h, l, _angles(gamma)),
            [0.2, 0.4, 0.6],
        )

    def test_box_full_width_and_gaussian_fwhm(self):
        box = CTRresolution.BoxResolution(0.4)
        gaussian = CTRresolution.GaussianResolution(0.4)

        np.testing.assert_array_equal(
            box.weights(np.array([-0.21, -0.2, 0.0, 0.2, 0.21]), 0.4),
            [0.0, 1.0, 1.0, 1.0, 0.0],
        )
        np.testing.assert_allclose(
            gaussian.weights(np.array([-0.2, 0.0, 0.2]), 0.4),
            [0.5, 1.0, 0.5],
        )

    def test_invalid_width_parameters_are_rejected(self):
        for value in (-0.1, np.inf, np.nan):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    CTRresolution.BoxResolution(value)
        with self.assertRaises(ValueError):
            CTRresolution.GaussianResolution(0.1, -0.2)

    def test_gamma_dependent_width_requires_angle_records(self):
        resolution = CTRresolution.BoxResolution(0.1, 0.2)
        with self.assertRaisesRegex(ValueError, "calcAnglesZmode"):
            resolution.width([0.0], [0.0], [1.0])


class TestFastConvolution(unittest.TestCase):
    def test_constant_intensity_is_preserved_on_unsorted_irregular_grid(self):
        l = np.array([3.0, 0.0, 1.0, 1.8])  # noqa: E741
        ctr = CTRplotutil.CTR((1.0, 0.0), l, np.full(l.shape, 4.0))
        ctrs = CTRplotutil.CTRCollection([ctr])

        result = CTRresolution.fast_convolve(
            ctrs, CTRresolution.GaussianResolution(1.2)
        )

        np.testing.assert_array_equal(result[0].l, l)
        np.testing.assert_allclose(result[0].sfI, 4.0)

    def test_box_uses_trapezoidal_weights_and_boundary_normalization(self):
        l = np.array([0.0, 1.0, 3.0])  # noqa: E741
        ctr = CTRplotutil.CTR((0.0, 0.0), l, np.array([1.0, 2.0, 4.0]))
        result = CTRresolution.fast_convolve(
            CTRplotutil.CTRCollection([ctr]),
            CTRresolution.BoxResolution(4.0),
        )[0]

        expected_intensity = np.array(
            [
                (0.5 * 1.0 + 1.5 * 4.0) / 2.0,
                (0.5 * 1.0 + 1.5 * 4.0 + 1.0 * 16.0) / 3.0,
                (1.5 * 4.0 + 1.0 * 16.0) / 2.5,
            ]
        )
        np.testing.assert_allclose(result.sfI, np.sqrt(expected_intensity))

    def test_zero_width_is_identity_and_output_metadata_is_clean(self):
        l = np.array([0.0, 0.4, 1.0])  # noqa: E741
        ctr = CTRplotutil.CTR(
            (1.0, 2.0),
            l,
            np.array([2.0, 3.0, 5.0]),
            err=np.ones(3),
            phi=np.array([10.0, 20.0, 30.0]),
            name="source rod",
        )
        ctr.angles = _angles([0.0, 0.1, 0.2])
        ctr.bgI = np.ones(3)
        ctrs = CTRplotutil.CTRCollection([ctr], name="source collection")
        ctrs.setPlotSettings("plot argument", color="red")

        result = CTRresolution.fast_convolve(
            ctrs, CTRresolution.GaussianResolution(0.0)
        )

        self.assertIsNot(result, ctrs)
        self.assertIsNot(result[0], ctr)
        self.assertEqual(result.name, ctrs.name)
        self.assertEqual(result.plotsett, ctrs.plotsett)
        self.assertEqual(result.plotkeyargs, ctrs.plotkeyargs)
        self.assertEqual(result[0].name, ctr.name)
        np.testing.assert_array_equal(result[0].sfI, ctr.sfI)
        np.testing.assert_array_equal(result[0].angles, ctr.angles)
        self.assertIsNone(result[0].err)
        self.assertIsNone(result[0].phi)
        self.assertFalse(hasattr(result[0], "bgI"))
        np.testing.assert_array_equal(ctr.err, np.ones(3))
        np.testing.assert_array_equal(ctr.phi, [10.0, 20.0, 30.0])

    def test_gamma_dependent_fast_convolution_requires_angles(self):
        ctr = CTRplotutil.CTR((0.0, 0.0), [0.0, 1.0], [1.0, 2.0])
        with self.assertRaisesRegex(ValueError, "calcAnglesZmode"):
            CTRresolution.fast_convolve(
                CTRplotutil.CTRCollection([ctr]),
                CTRresolution.BoxResolution(0.1, 0.2),
            )

    def test_invalid_ctr_data_are_rejected(self):
        cases = [
            CTRplotutil.CTR((0.0, 0.0), [0.0, 0.0], [1.0, 2.0]),
            CTRplotutil.CTR((0.0, 0.0), [0.0, np.nan], [1.0, 2.0]),
            CTRplotutil.CTR((0.0, 0.0), [0.0, 1.0], [1.0, -2.0]),
        ]
        for ctr in cases:
            with self.subTest(ctr=ctr):
                with self.assertRaises(ValueError):
                    CTRresolution.fast_convolve(
                        CTRplotutil.CTRCollection([ctr]),
                        CTRresolution.BoxResolution(0.2),
                    )


class TestStructureFactorSampling(unittest.TestCase):
    def test_box_quadrature_matches_quadratic_intensity_average(self):
        center = 2.0
        width = 0.6
        ctr = CTRplotutil.CTR((0.0, 0.0), [center], [0.0])
        result = CTRresolution.sample_structure_factor(
            CTRplotutil.CTRCollection([ctr]),
            QuadraticCrystal(),
            CTRresolution.BoxResolution(width),
            quadrature_order=5,
        )[0]

        expected_intensity = center**2 + width**2 / 12.0
        np.testing.assert_allclose(result.sfI, np.sqrt(expected_intensity))

    def test_gaussian_quadrature_matches_quadratic_intensity_average(self):
        center = 2.0
        fwhm = 0.6
        ctr = CTRplotutil.CTR((0.0, 0.0), [center], [0.0])
        result = CTRresolution.sample_structure_factor(
            CTRplotutil.CTRCollection([ctr]),
            QuadraticCrystal(),
            CTRresolution.GaussianResolution(fwhm),
            quadrature_order=5,
        )[0]

        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        expected_intensity = center**2 + sigma**2
        np.testing.assert_allclose(result.sfI, np.sqrt(expected_intensity))

    def test_zero_width_samples_the_central_structure_factor(self):
        ctr = CTRplotutil.CTR((1.0, 2.0), [0.5, 1.5], [99.0, 99.0])
        result = CTRresolution.sample_structure_factor(
            CTRplotutil.CTRCollection([ctr]),
            QuadraticCrystal(),
            CTRresolution.BoxResolution(0.0),
        )[0]
        np.testing.assert_allclose(result.sfI, [0.5, 1.5])

    def test_dense_irregular_fast_result_converges_to_slow_result(self):
        left = np.linspace(0.2, 1.0, 401) ** 1.3
        right = 1.0 + np.linspace(0.0, 0.8, 401)[1:] ** 1.7
        l = np.concatenate((left, right))  # noqa: E741
        ctr = CTRplotutil.CTR((0.0, 0.0), l, l)
        ctrs = CTRplotutil.CTRCollection([ctr])
        resolution = CTRresolution.GaussianResolution(0.3)

        fast = CTRresolution.fast_convolve(ctrs, resolution)[0]
        center_index = np.flatnonzero(l == 1.0)[0]
        template = CTRplotutil.CTR((0.0, 0.0), [1.0], [0.0])
        slow = CTRresolution.sample_structure_factor(
            CTRplotutil.CTRCollection([template]),
            QuadraticCrystal(),
            resolution,
        )[0]

        np.testing.assert_allclose(
            fast.sfI[center_index], slow.sfI[0], rtol=2e-4
        )

    def test_invalid_quadrature_order_is_rejected(self):
        ctrs = CTRplotutil.CTRCollection(
            [CTRplotutil.CTR((0.0, 0.0), [1.0], [0.0])]
        )
        for order in (0, 2, 2.5, True):
            with self.subTest(order=order):
                with self.assertRaises(ValueError):
                    CTRresolution.sample_structure_factor(
                        ctrs,
                        QuadraticCrystal(),
                        CTRresolution.BoxResolution(0.1),
                        quadrature_order=order,
                    )


class RecordingAngles:
    def __init__(self):
        self.calls = []

    def anglesZmode(self, hkl, fixedangle, **kwargs):
        self.calls.append((np.copy(hkl), fixedangle, kwargs))
        count = hkl.shape[1]
        values = np.zeros((count, 6), dtype=np.float64)
        values[:, 2] = np.arange(count) + 0.25
        return values


class TestCollectionAngleAPI(unittest.TestCase):
    def test_collection_forwards_all_z_mode_constraints(self):
        ctrs = CTRplotutil.CTRCollection(
            [
                CTRplotutil.CTR((0.0, 0.0), [0.1, 0.2], [1.0, 1.0]),
                CTRplotutil.CTR((1.0, 0.0), [0.3], [1.0]),
            ]
        )
        calculator = RecordingAngles()

        result = ctrs.calcAnglesZmode(
            calculator,
            fixedangle=0.12,
            fixed="out",
            chi=0.03,
            phi=-0.04,
            mirrorx=True,
        )

        self.assertEqual(len(result), 2)
        self.assertEqual(len(calculator.calls), 2)
        for call in calculator.calls:
            self.assertEqual(call[1], 0.12)
            self.assertEqual(
                call[2],
                {"fixed": "out", "chi": 0.03, "phi": -0.04, "mirrorx": True},
            )
        np.testing.assert_allclose(ctrs[0].angles["gamma"], [0.25, 1.25])
        np.testing.assert_allclose(ctrs[1].angles["gamma"], [0.25])


if __name__ == "__main__":
    unittest.main()

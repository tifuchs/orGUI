import numpy as np

from orgui.app import imagePeakFinder
from orgui.backend.scans import SimulationScan


def test_calc_image_range_subtracts_background_for_max_sum_and_rocking_curve():
    scan = SimulationScan((4, 4), 0.0, 2.0, 3)
    scan.images[:] = 10.0
    scan.images[1, 2, 2] = 50.0
    background = np.ones((4, 4), dtype=float) * 2.0
    background[2, 2] = 5.0

    images = imagePeakFinder.calc_image_range(
        scan,
        [0.0, 2.0],
        background_image=background,
        rocking_curve={"xy_start": [2.0, 2.0], "xsize": 2, "ysize": 2},
    )

    expected0 = scan.images[0] - background
    expected1 = scan.images[1] - background
    expected2 = scan.images[2] - background
    np.testing.assert_allclose(images["sum"], expected0 + expected1 + expected2)
    np.testing.assert_allclose(
        images["max"],
        np.maximum(np.maximum(expected0, expected1), expected2),
    )
    np.testing.assert_allclose(
        images["rocking_curve"],
        [
            np.nansum(expected0[1:3, 1:3]),
            np.nansum(expected1[1:3, 1:3]),
            np.nansum(expected2[1:3, 1:3]),
        ],
    )

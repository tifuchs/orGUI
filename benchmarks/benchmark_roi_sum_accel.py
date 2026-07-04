"""Benchmark ROI summing backends on synthetic Pilatus 6M images.

The benchmark compares three implementations of the same integration kernel:

- NumPy reference implementation
- C++ pybind11 extension, ``_roi_sum_cpp``
- Optional local Numba implementation, if numba is installed

The synthetic detector defaults to a Pilatus 6M ``(2527, 2463)`` image with
40 randomly placed center ROIs. Center ROI widths and heights are sampled from
10 to 100 pixels, and four background strips are generated around each ROI
with strip widths sampled from 5 to 50 pixels.

Run from the ``benchmarks/`` directory with orGUI already installed in the
active Python environment:

    python benchmark_roi_sum_accel.py
    python benchmark_roi_sum_accel.py --disk-images 16 --threads 4

The script imports ``orgui.app._roi_sum_cpp`` from the active Python
environment; if orGUI is not installed, the import fails.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import importlib
import json
import platform
import tempfile
import statistics
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np


PILATUS_6M_SHAPE = (2527, 2463)


def import_cpp_backend():
    """Import the C++ ROI backend from the active Python environment."""
    try:
        return importlib.import_module("orgui.app._roi_sum_cpp")
    except Exception as package_error:
        raise RuntimeError(
            "Could not import orgui.app._roi_sum_cpp from the active Python "
            "environment. Install orGUI with the C++ extension before running "
            "this benchmark."
        ) from package_error


def try_build_numba_backend():
    """Return numba-backed functions, or ``None`` if numba is unavailable."""
    try:
        from numba import njit
    except Exception:
        return None

    @njit(inline="always")
    def nansum_roi(data, roi, n):
        total = 0.0
        x0 = roi[n, 0, 0]
        x1 = roi[n, 0, 1]
        y0 = roi[n, 1, 0]
        y1 = roi[n, 1, 1]
        for y in range(y0, y1):
            for x in range(x0, x1):
                value = data[y, x]
                if not np.isnan(value):
                    total += value
        return total

    @njit(inline="always")
    def count_unmasked(mask, roi, n):
        total = 0.0
        x0 = roi[n, 0, 0]
        x1 = roi[n, 0, 1]
        y0 = roi[n, 1, 0]
        y1 = roi[n, 1, 1]
        for y in range(y0, y1):
            for x in range(x0, x1):
                if not mask[y, x]:
                    total += 1.0
        return total

    @njit(nogil=True, cache=False)
    def processImage_Carr(
        image,
        mask,
        correction,
        center,
        left,
        right,
        top,
        bottom,
        all_counters,
        correction_counters,
    ):
        for y in range(image.shape[0]):
            for x in range(image.shape[1]):
                if mask[y, x]:
                    image[y, x] = np.nan

        for i in range(center.shape[0]):
            all_counters[i, 0] = nansum_roi(image, center, i)
            all_counters[i, 2] = (
                nansum_roi(image, left, i)
                + nansum_roi(image, right, i)
                + nansum_roi(image, top, i)
                + nansum_roi(image, bottom, i)
            )
            correction_counters[i, 0] = nansum_roi(correction, center, i)
            correction_counters[i, 2] = (
                nansum_roi(correction, left, i)
                + nansum_roi(correction, right, i)
                + nansum_roi(correction, top, i)
                + nansum_roi(correction, bottom, i)
            )
            signal_pixels = count_unmasked(mask, center, i)
            background_pixels = (
                count_unmasked(mask, left, i)
                + count_unmasked(mask, right, i)
                + count_unmasked(mask, top, i)
                + count_unmasked(mask, bottom, i)
            )
            all_counters[i, 1] = signal_pixels
            all_counters[i, 3] = background_pixels
            correction_counters[i, 1] = signal_pixels
            correction_counters[i, 3] = background_pixels

    @njit(nogil=True, cache=False)
    def processImage_bg_Carr(
        image,
        background,
        mask,
        correction,
        center,
        left,
        right,
        top,
        bottom,
        all_counters,
        correction_counters,
        background_counters,
    ):
        processImage_Carr(
            image,
            mask,
            correction,
            center,
            left,
            right,
            top,
            bottom,
            all_counters,
            correction_counters,
        )
        for i in range(center.shape[0]):
            background_counters[i, 0] = nansum_roi(background, center, i)
            background_counters[i, 2] = (
                nansum_roi(background, left, i)
                + nansum_roi(background, right, i)
                + nansum_roi(background, top, i)
                + nansum_roi(background, bottom, i)
            )
            background_counters[i, 1] = all_counters[i, 1]
            background_counters[i, 3] = all_counters[i, 3]

    return {
        "processImage_Carr": processImage_Carr,
        "processImage_bg_Carr": processImage_bg_Carr,
    }


def make_rois(
    rng: np.random.Generator,
    shape: tuple[int, int],
    n_rois: int,
    roi_min: int,
    roi_max: int,
    bg_min: int,
    bg_max: int,
):
    """Return center and four background ROI arrays in accelerator format."""
    height, width = shape
    centers = []
    lefts = []
    rights = []
    tops = []
    bottoms = []
    for _ in range(n_rois):
        roi_w = int(rng.integers(roi_min, roi_max + 1))
        roi_h = int(rng.integers(roi_min, roi_max + 1))
        bg_w = int(rng.integers(bg_min, bg_max + 1))
        bg_h = int(rng.integers(bg_min, bg_max + 1))
        x0 = int(rng.integers(bg_w, width - roi_w - bg_w + 1))
        y0 = int(rng.integers(bg_h, height - roi_h - bg_h + 1))
        x1 = x0 + roi_w
        y1 = y0 + roi_h
        centers.append([[x0, x1], [y0, y1]])
        lefts.append([[x0 - bg_w, x0], [y0, y1]])
        rights.append([[x1, x1 + bg_w], [y0, y1]])
        tops.append([[x0, x1], [y0 - bg_h, y0]])
        bottoms.append([[x0, x1], [y1, y1 + bg_h]])
    return tuple(
        np.ascontiguousarray(np.asarray(values), dtype=np.int64)
        for values in (centers, lefts, rights, tops, bottoms)
    )


def make_inputs(args):
    """Create synthetic image data and ROI definitions."""
    rng = np.random.default_rng(args.seed)
    rois = make_rois(
        rng,
        args.shape,
        args.rois,
        args.roi_min,
        args.roi_max,
        args.bg_min,
        args.bg_max,
    )
    image = rng.poisson(args.signal_mean, size=args.shape).astype(np.float64)
    background = rng.uniform(
        args.background_min,
        args.background_max,
        size=args.shape,
    ).astype(np.float64)
    correction = rng.uniform(0.85, 1.15, size=args.shape).astype(np.float64)
    mask = rng.random(args.shape) < args.mask_fraction
    correction[mask] = np.nan
    background[mask] = np.nan
    return image, background, mask, correction, rois


def process_numpy(
    image,
    background,
    mask,
    correction,
    rois,
    all_counters,
    correction_counters,
    background_counters,
):
    """Reference NumPy implementation matching the accelerated API."""
    image[mask] = np.nan
    for i in range(rois[0].shape[0]):
        center = rois[0][i]
        all_counters[i, 0] = np.nansum(
            image[center[1, 0] : center[1, 1], center[0, 0] : center[0, 1]]
        )
        correction_counters[i, 0] = np.nansum(
            correction[
                center[1, 0] : center[1, 1],
                center[0, 0] : center[0, 1],
            ]
        )
        background_counters[i, 0] = np.nansum(
            background[
                center[1, 0] : center[1, 1],
                center[0, 0] : center[0, 1],
            ]
        )
        all_counters[i, 1] = correction_counters[i, 1] = background_counters[
            i,
            1,
        ] = np.sum(
            ~mask[
                center[1, 0] : center[1, 1],
                center[0, 0] : center[0, 1],
            ]
        )
        for roi_group in rois[1:]:
            roi = roi_group[i]
            ys = slice(roi[1, 0], roi[1, 1])
            xs = slice(roi[0, 0], roi[0, 1])
            all_counters[i, 2] += np.nansum(image[ys, xs])
            correction_counters[i, 2] += np.nansum(correction[ys, xs])
            background_counters[i, 2] += np.nansum(background[ys, xs])
            all_counters[i, 3] += np.sum(~mask[ys, xs])
        correction_counters[i, 3] = all_counters[i, 3]
        background_counters[i, 3] = all_counters[i, 3]


def run_one_backend(name, func, base_image, background, mask, correction, rois):
    """Run one backend once and return duration plus output counters."""
    image = base_image.copy(order="C")
    all_counters = np.zeros((rois[0].shape[0], 4), dtype=np.float64)
    correction_counters = np.zeros_like(all_counters)
    background_counters = np.zeros_like(all_counters)
    start = time.perf_counter()
    func(
        image,
        background,
        mask,
        correction,
        *rois,
        all_counters,
        correction_counters,
        background_counters,
    )
    duration = time.perf_counter() - start
    return duration, (all_counters, correction_counters, background_counters)


def integrate_image(func, image, background, mask, correction, rois):
    """Integrate one in-memory image with one backend function."""
    all_counters = np.zeros((rois[0].shape[0], 4), dtype=np.float64)
    correction_counters = np.zeros_like(all_counters)
    background_counters = np.zeros_like(all_counters)
    func(
        image,
        background,
        mask,
        correction,
        *rois,
        all_counters,
        correction_counters,
        background_counters,
    )
    return all_counters, correction_counters, background_counters


def benchmark_backend(
    name,
    func,
    base_image,
    background,
    mask,
    correction,
    rois,
    repeats,
):
    """Return median timing and last counters for one backend."""
    warmup_duration, counters = run_one_backend(
        name,
        func,
        base_image,
        background,
        mask,
        correction,
        rois,
    )
    timings = []
    for _ in range(repeats):
        duration, counters = run_one_backend(
            name,
            func,
            base_image,
            background,
            mask,
            correction,
            rois,
        )
        timings.append(duration)
    return {
        "name": name,
        "warmup_seconds": warmup_duration,
        "median_seconds": statistics.median(timings),
        "min_seconds": min(timings),
        "max_seconds": max(timings),
        "counters": counters,
    }


def make_disk_images(args, directory: Path):
    """Write random images to ``directory`` and return their paths."""
    rng = np.random.default_rng(args.seed + 1)
    paths = []
    for i in range(args.disk_images):
        image = rng.poisson(args.signal_mean, size=args.shape).astype(
            np.float64
        )
        path = directory / f"pilatus6m_random_{i:04d}.npy"
        np.save(path, image)
        paths.append(path)
    return paths


def process_disk_image(path, func, background, mask, correction, rois):
    """Read one image from disk and integrate it."""
    image = np.load(path)
    return integrate_image(func, image, background, mask, correction, rois)


def checksum_counters(counter_sets):
    """Return a deterministic checksum over image-level counter tuples."""
    return float(
        sum(
            np.nansum(counter)
            for counters in counter_sets
            for counter in counters
        )
    )


def benchmark_disk_backend(
    name,
    func,
    image_paths,
    background,
    mask,
    correction,
    rois,
    threads,
):
    """Benchmark reading images from disk and integrating them."""
    start = time.perf_counter()
    serial_counters = [
        process_disk_image(path, func, background, mask, correction, rois)
        for path in image_paths
    ]
    serial_seconds = time.perf_counter() - start

    start = time.perf_counter()
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                process_disk_image,
                path,
                func,
                background,
                mask,
                correction,
                rois,
            )
            for path in image_paths
        ]
        threaded_counters = [future.result() for future in futures]
    threaded_seconds = time.perf_counter() - start

    return {
        "name": name,
        "serial_seconds": serial_seconds,
        "threaded_seconds": threaded_seconds,
        "thread_speedup": serial_seconds / threaded_seconds,
        "serial_checksum": checksum_counters(serial_counters),
        "threaded_checksum": checksum_counters(threaded_counters),
        "serial_counters": serial_counters,
        "threaded_counters": threaded_counters,
    }


def wrap_numpy_backend():
    """Return a function with the C++/numba call signature."""

    def run(
        image,
        background,
        mask,
        correction,
        center,
        left,
        right,
        top,
        bottom,
        all_counters,
        correction_counters,
        background_counters,
    ):
        process_numpy(
            image,
            background,
            mask,
            correction,
            (center, left, right, top, bottom),
            all_counters,
            correction_counters,
            background_counters,
        )

    return run


def wrap_backend(module):
    """Return a function with the full background-image benchmark signature."""

    def run(
        image,
        background,
        mask,
        correction,
        center,
        left,
        right,
        top,
        bottom,
        all_counters,
        correction_counters,
        background_counters,
    ):
        module.processImage_bg_Carr(
            image,
            background,
            mask,
            correction,
            center,
            left,
            right,
            top,
            bottom,
            all_counters,
            correction_counters,
            background_counters,
        )

    return run


def check_results(reference, result, backend):
    """Raise if backend counters do not match the NumPy reference."""
    labels = ("image", "correction", "background")
    for label, expected, actual in zip(labels, reference, result, strict=True):
        np.testing.assert_allclose(
            actual,
            expected,
            rtol=1e-12,
            atol=1e-9,
            err_msg=f"{backend} {label} counters differ from NumPy",
        )


def print_results(results):
    """Print benchmark results as a compact table."""
    numpy_time = results["numpy"]["median_seconds"]
    print("\n=== Results ===")
    print(f"{'backend':>10s} {'median_s':>12s} {'min_s':>12s} {'speedup':>10s}")
    for name, result in results.items():
        speedup = numpy_time / result["median_seconds"]
        print(
            f"{name:>10s} "
            f"{result['median_seconds']:12.6g} "
            f"{result['min_seconds']:12.6g} "
            f"{speedup:10.2f}x"
        )


def print_disk_results(results):
    """Print disk-read benchmark results as a compact table."""
    print("\n=== Disk Read + Threaded Integration Results ===")
    print(
        f"{'backend':>10s} {'serial_s':>12s} "
        f"{'threaded_s':>12s} {'thread_gain':>12s}"
    )
    for name, result in results.items():
        print(
            f"{name:>10s} "
            f"{result['serial_seconds']:12.6g} "
            f"{result['threaded_seconds']:12.6g} "
            f"{result['thread_speedup']:12.2f}x"
        )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--rois", type=int, default=40)
    parser.add_argument("--shape", type=int, nargs=2, default=PILATUS_6M_SHAPE)
    parser.add_argument("--roi-min", type=int, default=10)
    parser.add_argument("--roi-max", type=int, default=100)
    parser.add_argument("--bg-min", type=int, default=5)
    parser.add_argument("--bg-max", type=int, default=50)
    parser.add_argument("--signal-mean", type=float, default=100.0)
    parser.add_argument("--background-min", type=float, default=5.0)
    parser.add_argument("--background-max", type=float, default=50.0)
    parser.add_argument("--mask-fraction", type=float, default=0.01)
    parser.add_argument(
        "--disk-images",
        type=int,
        default=0,
        help="If positive, write this many random .npy images and benchmark "
        "read + integrate in serial and with threads.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Thread count for --disk-images read + integrate benchmark.",
    )
    parser.add_argument(
        "--disk-dir",
        type=Path,
        help="Directory for temporary random .npy images. Defaults to a "
        "TemporaryDirectory that is deleted after the run.",
    )
    parser.add_argument("--json", type=Path, help="Optional JSON output path.")
    parser.add_argument(
        "--no-numba",
        action="store_true",
        help="Skip the optional numba comparison even if numba is installed.",
    )
    return parser.parse_args()


def main():
    """Run the ROI backend benchmark."""
    args = parse_args()
    args.shape = tuple(args.shape)
    image, background, mask, correction, rois = make_inputs(args)
    cpp_backend = import_cpp_backend()

    backends = {"numpy": wrap_numpy_backend(), "cpp": wrap_backend(cpp_backend)}
    if not args.no_numba:
        numba_backend = try_build_numba_backend()
        if numba_backend is not None:
            backends["numba"] = wrap_backend(SimpleNamespace(**numba_backend))
        else:
            print("numba unavailable; skipping numba backend")

    print("=== ROI Benchmark ===")
    print("detector Pilatus 6M")
    print("shape", args.shape)
    print("rois", args.rois)
    print("roi_size_px", (args.roi_min, args.roi_max))
    print("background_strip_px", (args.bg_min, args.bg_max))
    print("background_counts", (args.background_min, args.background_max))
    print("mask_fraction", args.mask_fraction)
    print("repeats", args.repeats)
    if args.disk_images:
        print("disk_images", args.disk_images)
        print("threads", args.threads)
    print("python", platform.python_version())
    print("numpy", np.__version__)
    print("cpp_backend", getattr(cpp_backend, "__file__", "installed module"))

    results = {}
    reference = None
    for name, func in backends.items():
        result = benchmark_backend(
            name,
            func,
            image,
            background,
            mask,
            correction,
            rois,
            args.repeats,
        )
        if name == "numpy":
            reference = result["counters"]
        else:
            check_results(reference, result["counters"], name)
        counters = result.pop("counters")
        result["counter_checksum"] = float(
            sum(np.nansum(counter) for counter in counters)
        )
        results[name] = result

    print_results(results)

    disk_results = {}
    if args.disk_images:
        if args.disk_dir:
            args.disk_dir.mkdir(parents=True, exist_ok=True)
            image_paths = make_disk_images(args, args.disk_dir)
            cleanup = None
            disk_dir = args.disk_dir
        else:
            cleanup = tempfile.TemporaryDirectory()
            disk_dir = Path(cleanup.name)
            image_paths = make_disk_images(args, disk_dir)
        print("\nDisk image directory", disk_dir)
        print("Disk image format .npy float64")
        for name, func in backends.items():
            disk_result = benchmark_disk_backend(
                name,
                func,
                image_paths,
                background,
                mask,
                correction,
                rois,
                args.threads,
            )
            for i, (serial, threaded) in enumerate(
                zip(
                    disk_result["serial_counters"],
                    disk_result["threaded_counters"],
                    strict=True,
                )
            ):
                check_results(serial, threaded, f"{name} threaded image {i}")
            del disk_result["serial_counters"]
            del disk_result["threaded_counters"]
            disk_results[name] = disk_result
        print_disk_results(disk_results)
        if cleanup is not None:
            cleanup.cleanup()

    if args.json:
        payload = {
            "config": {
                "detector": "Pilatus 6M",
                "shape": args.shape,
                "rois": args.rois,
                "roi_size_px": [args.roi_min, args.roi_max],
                "background_strip_px": [args.bg_min, args.bg_max],
                "background_counts": [
                    args.background_min,
                    args.background_max,
                ],
                "mask_fraction": args.mask_fraction,
                "repeats": args.repeats,
                "disk_images": args.disk_images,
                "threads": args.threads,
                "seed": args.seed,
            },
            "results": results,
            "disk_results": disk_results,
        }
        args.json.write_text(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

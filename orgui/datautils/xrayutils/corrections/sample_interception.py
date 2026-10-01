"""Flat sample/beam overlap: metres, radians, normalized densities [1/m].

The surface x axis follows projected beam propagation, x cross y is the
surface normal, and rotation is right-handed about that normal. Compute H
directly on the surface; the projection Jacobian gives f_hit = sin(alpha) H.
The incident beam is separable, even when the clipped overlap is not.
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _rotation(angle):
    return np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])


@dataclass(frozen=True)
class SampleShape:
    """Shape in surface coordinates [m], about an explicitly supplied origin.

    :param str kind: ``rectangle``, ``circle`` or ``polygon``.
    :param tuple dimensions: Rectangle length/width or circle diameter [m].
    :param tuple vertices: Polygon coordinates [m]; either winding is accepted.
    """

    kind: str
    dimensions: tuple = ()
    vertices: tuple = ()

    def __post_init__(self):
        if self.kind not in {"rectangle", "circle", "polygon"}:
            raise ValueError("sample shape must be rectangle, circle or polygon")
        dimensions = tuple(float(v) for v in self.dimensions)
        object.__setattr__(self, "dimensions", dimensions)
        if self.kind != "polygon":
            count = 2 if self.kind == "rectangle" else 1
            if len(dimensions) != count or not np.all(np.isfinite(dimensions)):
                raise ValueError("invalid sample dimensions, in metres")
            if min(dimensions) <= 0:
                raise ValueError("sample dimensions must be positive")
            return
        vertices = np.asarray(self.vertices, dtype=float)
        if vertices.ndim != 2 or vertices.shape[1] != 2:
            raise ValueError("polygon must contain (x, y) vertices")
        if len(vertices) > 1 and np.array_equal(vertices[0], vertices[-1]):
            vertices = vertices[:-1]
        if len(vertices) < 3 or not np.all(np.isfinite(vertices)):
            raise ValueError("polygon needs at least three finite vertices")
        scale = np.max(np.ptp(vertices, axis=0))
        if scale <= 0:
            raise ValueError("polygon has zero area")
        v = (vertices - vertices[0]) / scale
        if len(np.unique(v, axis=0)) != len(v):
            raise ValueError("polygon has repeated vertices")
        edges = np.roll(v, -1, axis=0) - v
        for previous, current in zip(np.roll(edges, 1, axis=0), edges):
            if abs(_cross(previous, current)) <= 1e-14 and previous @ current < 0:
                raise ValueError("polygon has overlapping adjacent edges")
        area = sum(_cross(a, b) for a, b in zip(v, np.roll(v, -1, axis=0)))
        if abs(area) <= 1e-14:
            raise ValueError("polygon has zero area")
        for i in range(len(v)):
            for j in range(i + 1, len(v)):
                if j == i + 1 or (i == 0 and j == len(v) - 1):
                    continue
                a, b, c, d = v[i], v[i] + edges[i], v[j], v[j] + edges[j]
                signs = [
                    _cross(b - a, c - a),
                    _cross(b - a, d - a),
                    _cross(d - c, a - c),
                    _cross(d - c, b - c),
                ]
                if (
                    signs[0] * signs[1] <= 0
                    and signs[2] * signs[3] <= 0
                    and np.all(
                        np.maximum(np.minimum(a, b), np.minimum(c, d))
                        <= np.minimum(np.maximum(a, b), np.maximum(c, d))
                    )
                ):
                    raise ValueError("polygon is self-intersecting")
        object.__setattr__(self, "vertices", tuple(map(tuple, vertices)))

    @property
    def area(self):
        """Geometric surface area, in square metres."""
        if self.kind == "circle":
            return np.pi * (self.dimensions[0] / 2) ** 2
        if self.kind == "rectangle":
            return np.prod(self.dimensions)
        v = np.asarray(self.vertices)
        return abs(sum(_cross(a, b) for a, b in zip(v, np.roll(v, -1, axis=0)))) / 2

    def outline(self, azimuth=0.0, offset=(0.0, 0.0), orientation=0.0):
        """Transformed outline [m]; azimuth/orientation are radians.

        The sample-fixed ``offset`` rotates with the shape. A circle outline
        is sampled for display only; integration uses exact circle sections.
        """
        if self.kind == "rectangle":
            length, width = self.dimensions
            v = (
                np.array(
                    [
                        [-length, -width],
                        [length, -width],
                        [length, width],
                        [-length, width],
                    ]
                )
                / 2
            )
        elif self.kind == "circle":
            angles = np.linspace(0, 2 * np.pi, 129)
            v = (
                self.dimensions[0]
                / 2
                * np.column_stack((np.cos(angles), np.sin(angles)))
            )
        else:
            v = np.asarray(self.vertices)
        return (v @ _rotation(orientation).T + offset) @ _rotation(azimuth).T


@dataclass(frozen=True)
class InterceptionResult:
    """Overlap arrays and quadrature error on H, with deterministic inputs."""

    illumination: np.ndarray
    fraction: np.ndarray
    error: np.ndarray
    effective_area: np.ndarray


def overlap(
    shape,
    vertical,
    horizontal,
    alpha,
    azimuth,
    *,
    offset=(0.0, 0.0),
    orientation=0.0,
    axis_vertical=0.0,
    axis_horizontal=0.0,
    rtol=1e-9,
    progress=None,
):
    """Integrate separable beam density over the actual transformed shape.

    :param SampleShape shape: Flat surface shape in metres.
    :param vertical: Normalized vertical profile, coordinates relative to its centre.
    :param horizontal: Normalized horizontal profile, same coordinate convention.
    :param alpha: Physical incidence angle(s) in radians, ``(0, pi/2]``.
    :param azimuth: Geometric rotation(s) in radians, broadcast with ``alpha``.
    :param offset: Sample-fixed displacement from the rotation axis [m].
    :param float orientation: Shape orientation at reference azimuth [rad].
    :param axis_vertical: Axis origin in vertical profile coordinates [m].
    :param axis_horizontal: Axis origin in horizontal profile coordinates [m].
    :param float rtol: Relative integration tolerance.
    :param progress: Optional ``(completed, total)`` callback; return ``False``
        to cancel between frames. No GUI objects enter this calculation.
    :returns: H, f_hit, estimated H error and peak-referenced effective area [m²].
    :rtype: InterceptionResult
    :raises ValueError: For invalid geometry/probability or unresolved quadrature.
    """
    if len(offset) != 2:
        raise ValueError("sample offset must have two coordinates, in metres")
    alpha, azimuth, zv, yh = np.broadcast_arrays(
        alpha, azimuth, axis_vertical, axis_horizontal
    )
    if (
        not np.all(
            np.isfinite(
                [
                    *np.ravel(alpha),
                    *np.ravel(azimuth),
                    *np.ravel(zv),
                    *np.ravel(yh),
                    *offset,
                    orientation,
                ]
            )
        )
        or np.any(alpha <= 0)
        or np.any(alpha > np.pi / 2)
    ):
        raise ValueError("finite angles in (0, pi/2] and finite placement are required")
    if not 1e-12 <= rtol <= 1e-3:
        raise ValueError("integration rtol must be between 1e-12 and 1e-3")
    for profile in (vertical, horizontal):
        if not np.isfinite(profile.peak_density) or profile.peak_density <= 0:
            raise ValueError("beam profile must have finite positive peak density")
        density = profile.density_at(profile.integration_points)
        if np.any(density < 0) or not np.all(np.isfinite(density)):
            raise ValueError("beam density must be finite and non-negative")
    illumination = np.empty(alpha.shape)
    error = np.empty(alpha.shape)
    # Deduplicate exact geometry only: no rounding of acquisition angles.
    cache = {}
    for completed, index in enumerate(np.ndindex(alpha.shape)):
        if progress is not None and progress(completed, alpha.size) is False:
            raise InterruptedError("sample interception cancelled")
        key = tuple(float(v[index]) for v in (alpha, azimuth, zv, yh))
        if key not in cache:
            cache[key] = _frame_overlap(
                shape, vertical, horizontal, *key, offset, orientation, rtol
            )
        illumination[index], error[index] = cache[key]
    fraction = np.sin(alpha) * illumination
    if np.any(fraction < -1e-12) or np.any(fraction > 1 + 1e-10):
        raise ValueError("overlap integration returned an invalid probability")
    return InterceptionResult(
        illumination,
        fraction,
        error,
        illumination / (vertical.peak_density * horizontal.peak_density),
    )


def _frame_overlap(
    shape, vertical, horizontal, alpha, azimuth, zv, yh, offset, orientation, rtol
):
    sine = np.sin(alpha)
    # An aligned rectangle separates exactly. Very narrow projected intervals
    # retain surface quadrature below, avoiding cancellation in generic CDFs.
    if (
        shape.kind == "rectangle"
        and (azimuth + orientation) % np.pi == 0
        and shape.dimensions[0] * sine * vertical.peak_density >= 1e-7
    ):
        center = _rotation(azimuth) @ np.asarray(offset)
        length, width = shape.dimensions
        zcenter = zv + center[0] * sine
        span = length * sine
        vertical_h = (
            vertical.interval_mass(zcenter - span / 2, zcenter + span / 2) / sine
        )
        horizontal_mass = horizontal.interval_mass(
            yh + center[1] - width / 2, yh + center[1] + width / 2
        )
        value = float(vertical_h * horizontal_mass)
        return value, abs(value) * np.finfo(float).eps * 8
    vertices = shape.outline(azimuth, offset, orientation)
    center = _rotation(azimuth) @ np.asarray(offset)
    lo, hi = np.min(vertices[:, 0]), np.max(vertices[:, 0])
    if shape.kind == "circle":
        radius = shape.dimensions[0] / 2
        lo, hi = center[0] - radius, center[0] + radius

        def intervals(x):
            half = np.sqrt(max(0.0, radius**2 - (x - center[0]) ** 2))
            return [(center[1] - half, center[1] + half)]
    else:
        ends = np.roll(vertices, -1, axis=0)
        edges = ends - vertices

        def intervals(x):
            hits = []
            for a, b, edge in zip(vertices, ends, edges):
                if min(a[0], b[0]) <= x < max(a[0], b[0]):
                    hits.append(a[1] + (x - a[0]) * edge[1] / edge[0])
            hits.sort()
            if len(hits) % 2:
                raise ValueError("polygon section has an odd number of intersections")
            return list(zip(hits[::2], hits[1::2]))

    peak = vertical.peak_density

    def integrand(t):
        x = lo + (hi - lo) * t
        mass = sum(
            float(horizontal.interval_mass(yh + a, yh + b)) for a, b in intervals(x)
        )
        return float(vertical.density_at(zv + x * sine)) / peak * mass

    # Break at shape vertices and vertical knots/scales. Horizontal transitions
    # across polygon edges are also explicit, so a narrow beam cannot be missed.
    points = list(vertices[:, 0]) + list((vertical.integration_points - zv) / sine)
    if shape.kind != "circle":
        for a, edge in zip(vertices, edges):
            if edge[1] != 0:
                for position in horizontal.integration_points:
                    t = (position - yh - a[1]) / edge[1]
                    if 0 < t < 1:
                        points.append(a[0] + t * edge[0])
    else:
        for position in horizontal.integration_points:
            height = position - yh - center[1]
            if abs(height) < radius:
                half = np.sqrt(radius**2 - height**2)
                points.extend([center[0] - half, center[0] + half])
    points = np.unique(
        [0.0, 1.0, *[(p - lo) / (hi - lo) for p in points if lo < p < hi]]
    )
    total = error = 0.0
    for a, b in zip(points[:-1], points[1:]):
        value = quad(
            integrand,
            a,
            b,
            epsabs=1e-12 * (b - a),
            epsrel=rtol,
            limit=200,
            full_output=1,
        )
        if len(value) != 3:
            raise ValueError("sample overlap quadrature did not converge: " + value[3])
        total += value[0]
        error += value[1]
    if error > max(1e-11, abs(total) * rtol * 5):
        raise ValueError("sample overlap quadrature exceeded its error tolerance")
    scale = (hi - lo) * peak
    return total * scale, error * scale

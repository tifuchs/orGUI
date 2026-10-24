"""Flat sample/beam overlap: metres, radians, normalized densities [1/m].

The surface x axis follows projected beam propagation, x cross y is the
surface normal, and rotation is right-handed about that normal. Compute H
directly on the surface; the projection Jacobian gives f_hit = sin(alpha) H.
The incident beam is separable, even when the clipped overlap is not.
"""

from dataclasses import dataclass

import numpy as np
from scipy import stats
from scipy.integrate import quad

from .beamprofile import (
    DistributionBeamProfile,
    GaussianBeamProfile,
    MeasuredBeamProfile,
    _SmoothedTopHat,
)

try:
    from .. import _sample_interception_cpp
except ImportError:
    _sample_interception_cpp = None


def _native_profile(profile):
    """Snapshot supported exact profiles in metres and normalized 1/m density.

    Never sample a distribution into an approximate lookup table. Custom
    profiles and overriding subclasses retain the Python implementation.
    """
    if type(profile) is GaussianBeamProfile:
        descriptor = dict(kind="normal", shift=0.0, scale=profile.sigma)
    elif type(profile) is MeasuredBeamProfile:
        descriptor = dict(kind="measured", shift=profile.sample_center,
                          positions=profile._z_raw, density=profile._p,
                          cumulative=profile._cum)
    elif type(profile) is DistributionBeamProfile:
        if profile._normal_distribution:
            descriptor = dict(kind="normal", scale=profile.rms_width,
                              shift=profile.sample_center - profile.centroid)
        elif profile._uniform_support is not None:
            descriptor = dict(kind="uniform", shift=profile.sample_center,
                              bounds=profile._uniform_support)
        else:
            descriptor = _native_distribution(profile)
            if descriptor is None:
                return None
    else:
        return None
    return dict(descriptor, peak=profile.peak_density,
                points=profile.integration_points)


def _native_distribution(profile):
    """Copy named distribution parameters; leave arbitrary SciPy families alone."""
    dist = profile._dist
    if type(dist) is _SmoothedTopHat:
        return dict(kind="smoothed", shift=profile.sample_center,
                    scale=dist.sigma, width=dist.width)
    family = getattr(dist, "dist", None)
    for name, parameters in (("trapezoid", ("c", "d")),
                             ("gennorm", ("beta",)), ("skewnorm", ("a",))):
        if type(family) is not type(getattr(stats, name)):
            continue
        # Frozen distributions allow shape, location and scale both as
        # positional arguments and keywords. Location is in metres.
        values = dict(loc=0.0, scale=1.0)
        values.update(zip((*parameters, "loc", "scale"), dist.args))
        values.update(dist.kwds)
        return dict(kind=name, shift=profile.sample_center - values["loc"],
                    scale=values["scale"],
                    parameters=[values[key] for key in parameters])
    return None


def _native_overlap(shape, vertical, horizontal, arrays, offset, orientation,
                    rtol, progress, workers):
    """Batch exact geometries; retry unresolved panels with adaptive Python."""
    if _sample_interception_cpp is None:
        return None
    profiles = (_native_profile(vertical), _native_profile(horizontal))
    if any(profile is None for profile in profiles):
        return None
    frames = np.column_stack([array.ravel() for array in arrays])
    unique, inverse = np.unique(frames, axis=0, return_inverse=True)
    inverse = inverse.reshape(-1)
    values = np.empty((len(unique), 2))
    shape_data = dict(kind=shape.kind, dimensions=shape.dimensions,
                      vertices=shape.outline() if shape.kind != "circle" else [])
    completed = 0
    # Count repeated acquisition geometries in the source-frame progress total.
    multiplicity = np.bincount(inverse, minlength=len(unique))
    for start in range(0, len(unique), 256):
        if progress is not None and progress(completed, len(frames)) is False:
            raise InterruptedError("sample interception cancelled")
        stop = min(start + 256, len(unique))
        value, error, resolved = _sample_interception_cpp.overlap_batch(
            shape_data, *profiles, unique[start:stop], offset, orientation,
            rtol, workers,
        )
        for index in np.flatnonzero(~resolved):
            if progress is not None and progress(completed, len(frames)) is False:
                raise InterruptedError("sample interception cancelled")
            value[index], error[index] = _frame_overlap(
                shape, vertical, horizontal, *unique[start + index],
                offset, orientation, rtol,
            )
        values[start:stop, 0], values[start:stop, 1] = value, error
        completed += int(multiplicity[start:stop].sum())
    if progress is not None and progress(completed, len(frames)) is False:
        raise InterruptedError("sample interception cancelled")
    return tuple(values[inverse, i].reshape(arrays[0].shape) for i in range(2))


# QUADPACK's embedded Gauss 10 / Kronrod 21 rule, mapped to [0, 1].
# Odd entries are the Gauss nodes. Keep the same rule/error estimator as
# quad's first panel, but evaluate profile arrays once instead of 21 times.
_GK_ABSCISSA = np.array([
    0.9956571630258080807, 0.9739065285171717201,
    0.9301574913557082260, 0.8650633666889845107,
    0.7808177265864168971, 0.6794095682990244062,
    0.5627571346686046833, 0.4333953941292471908,
    0.2943928627014601981, 0.1488743389816312109, 0.0,
])
_GK_NODES = (1 + np.r_[_GK_ABSCISSA, -_GK_ABSCISSA[-2::-1]]) / 2
_GK_HALF_WEIGHTS = np.array([
    0.0116946388673718743, 0.0325581623079647275,
    0.0547558965743519960, 0.0750396748109199528,
    0.0931254545836976055, 0.1093871588022976419,
    0.1234919762620658511, 0.1347092173114733259,
    0.1427759385770600808, 0.1477391049013384914,
    0.1494455540029169057,
])
_GK_WEIGHTS = np.r_[_GK_HALF_WEIGHTS, _GK_HALF_WEIGHTS[-2::-1]] / 2
_GAUSS_HALF_WEIGHTS = np.array([
    0.0666713443086881376, 0.1494513491505805931,
    0.2190863625159820440, 0.2692667193099963551,
    0.2955242247147528702,
])
_GAUSS_WEIGHTS = np.r_[_GAUSS_HALF_WEIGHTS, _GAUSS_HALF_WEIGHTS[::-1]] / 2


def _vector_quad(integrand, rtol):
    """Try one vectorized GK21 panel on [0, 1], else retain adaptive quad."""
    values = integrand(_GK_NODES)
    result = float(np.sum(_GK_WEIGHTS * values))
    gauss = float(np.sum(_GAUSS_WEIGHTS * values[1::2]))
    deviation = float(np.sum(_GK_WEIGHTS * np.abs(values - result)))
    absolute = float(np.sum(_GK_WEIGHTS * np.abs(values)))
    error = abs(result - gauss)
    if deviation != 0 and error != 0:
        error = deviation * min(1.0, (200 * error / deviation) ** 1.5)
    error = max(error, 50 * np.finfo(float).eps * absolute)
    # A saturated error estimate needs subdivision, as in QUADPACK. The
    # absolute tolerance remains per local panel, before its surface weight.
    if error <= max(1e-12, rtol * abs(result)) and (
        error != deviation or error == 0
    ):
        return result, error, None
    return quad(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=rtol,
                limit=200, full_output=1)


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
    workers=0,
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
        to cancel between bounded native batches or Python frames. No GUI
        objects enter this calculation.
    :param int workers: Native worker limit; 0 selects workers automatically.
        1 runs the native batch serially. Unsupported profiles and unresolved
        quadrature panels retain the Python/SciPy implementation.
    :returns: H, f_hit, estimated H error and peak-referenced effective area [m²].
    :rtype: InterceptionResult
    :raises ValueError: For invalid geometry/probability or unresolved quadrature.
    """
    if len(offset) != 2:
        raise ValueError("sample offset must have two coordinates, in metres")
    if not isinstance(workers, (int, np.integer)) or workers < 0:
        raise ValueError("workers must be a non-negative integer")
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
    native = _native_overlap(shape, vertical, horizontal, (alpha, azimuth, zv, yh),
                             offset, orientation, rtol, progress, int(workers))
    if native is None:
        illumination, error = _python_overlap(
            shape, vertical, horizontal, (alpha, azimuth, zv, yh),
            offset, orientation, rtol, progress,
        )
    else:
        illumination, error = native
    fraction = np.sin(alpha) * illumination
    if np.any(fraction < -1e-12) or np.any(fraction > 1 + 1e-10):
        raise ValueError("overlap integration returned an invalid probability")
    return InterceptionResult(
        illumination,
        fraction,
        error,
        illumination / (vertical.peak_density * horizontal.peak_density),
    )


def _python_overlap(shape, vertical, horizontal, arrays, offset, orientation,
                    rtol, progress):
    """Reference scalar loop, also supporting user-defined beam profiles."""
    alpha, azimuth, zv, yh = arrays
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
    return illumination, error


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
    center = _rotation(azimuth) @ np.asarray(offset)
    if shape.kind == "circle":
        radius = shape.dimensions[0] / 2
        lo, hi = center[0] - radius, center[0] + radius
        span = 2 * radius

        def _section(t):
            # x = cx + R sin(theta), theta = pi(t - 1/2). The exact
            # half-chord R cos(theta) avoids cancellation at circle edges.
            # dx/dt / span = (pi/2) cos(theta) preserves surface area [m²].
            theta = np.pi * (t - 0.5)
            cosine = np.cos(theta)
            x = center[0] + radius * np.sin(theta)
            half = radius * cosine
            return x, [(center[1] - half, center[1] + half)], np.pi / 2 * cosine
    else:
        vertices = shape.outline(azimuth, offset, orientation)
        lo, hi = np.min(vertices[:, 0]), np.max(vertices[:, 0])
        span = hi - lo
        ends = np.roll(vertices, -1, axis=0)
        edges = ends - vertices

        def _section_for_interval(a, b):
            # Vertex/beam breakpoints delimit fixed pairs of crossing edges.
            # Interpolate their heights in a local unit interval: recovering
            # an edge parameter from a global x at every node amplifies x
            # roundoff by 1/edge[0] for almost beam-perpendicular edges.
            hits = []
            for start, end, edge in zip(vertices, ends, edges):
                if (edge[0] != 0 and min(start[0], end[0]) <= a
                        and b <= max(start[0], end[0])):
                    y0 = start[1] + ((a - start[0]) / edge[0]) * edge[1]
                    dy = ((b - a) / edge[0]) * edge[1]
                    hits.append((y0, dy))
            hits.sort(key=lambda hit: hit[0] + hit[1] / 2)
            if len(hits) % 2:
                raise ValueError("polygon section has an odd number of intersections")

            def _local_section(t):
                heights = [y0 + dy * t for y0, dy in hits]
                return a + (b - a) * t, list(zip(
                    heights[::2], heights[1::2]
                )), 1.0

            return _local_section

    peak = vertical.peak_density
    # Custom profiles historically only needed scalar interval_mass. Do not
    # impose a new array contract on them (including overriding subclasses).
    vector_profiles = (DistributionBeamProfile, GaussianBeamProfile,
                       MeasuredBeamProfile)
    vectorized = (type(vertical) in vector_profiles
                  and type(horizontal) in vector_profiles)

    def integrand(t):
        x, intervals_y, jacobian = section(t)
        mass = sum(
            horizontal.interval_mass(yh + a, yh + b) for a, b in intervals_y
        )
        return vertical.density_at(zv + x * sine) / peak * mass * jacobian

    # Break at shape vertices and vertical knots/scales. Horizontal transitions
    # across polygon edges are also explicit, so a narrow beam cannot be missed.
    points = list((vertical.integration_points - zv) / sine)
    horizontal_points = horizontal.integration_points
    if shape.kind != "circle":
        points.extend(vertices[:, 0])
        for a, edge in zip(vertices, edges):
            if edge[1] != 0:
                for position in horizontal_points:
                    t = (position - yh - a[1]) / edge[1]
                    if 0 < t < 1:
                        points.append(a[0] + t * edge[0])
        # Keep physical x breakpoints; normalizing and reconstructing them
        # loses precision in the narrow corner sections near alignment.
        points = [lo, hi, *[p for p in points if lo < p < hi]]
    else:
        # Exact circular chords need no vertices from the display outline.
        # Only beam transitions supply interior quadrature breakpoints.
        points = [
            0.5 + np.arcsin((p - center[0]) / radius) / np.pi
            for p in points if abs(p - center[0]) < radius
        ]
        for position in horizontal_points:
            height = position - yh - center[1]
            if abs(height) < radius:
                theta = np.arccos(abs(height) / radius)
                points.extend([0.5 - theta / np.pi, 0.5 + theta / np.pi])
        points = [0.0, 1.0, *points]
    points = np.unique(points)
    total = error = 0.0
    for a, b in zip(points[:-1], points[1:]):
        if shape.kind == "circle":
            def section(t):
                return _section(a + (b - a) * t)

            weight = b - a
        else:
            section = _section_for_interval(a, b)
            weight = (b - a) / span
        if weight <= 8 * np.finfo(float).eps:
            # Rotated edges can produce knots only a few ULPs apart in the
            # unit interval. QUADPACK cannot resolve their interior nodes.
            # Density/peak and transverse mass are at most one; the largest
            # section Jacobian is pi/2 (circle). Retain the interval's area
            # and bound its entire contribution as error, in normalized units.
            total += integrand(0.5) * weight
            error += weight * (np.pi / 2 if shape.kind == "circle" else 1.0)
            continue
        if vectorized:
            value = _vector_quad(integrand, rtol)
        else:
            value = quad(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=rtol,
                         limit=200, full_output=1)
        if len(value) != 3:
            raise ValueError("sample overlap quadrature did not converge: " + value[3])
        total += value[0] * weight
        error += value[1] * weight
    if error > max(1e-11, abs(total) * rtol * 5):
        raise ValueError("sample overlap quadrature exceeded its error tolerance")
    scale = span * peak
    return total * scale, error * scale

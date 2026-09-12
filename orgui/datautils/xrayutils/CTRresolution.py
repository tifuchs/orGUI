# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
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
"""Optional one-dimensional resolution modeling for calculated CTRs.

Resolution functions in this module act on field intensity, ``abs(F)**2``,
along the CTR L direction. The low-level kernels return intensity directly;
the collection compatibility wrappers convert their results back to effective
amplitudes for :class:`~.CTRplotutil.CTRCollection`.
"""

from abc import ABC, abstractmethod
import copy
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

from .CTRplotutil import CTR, CTRCollection, MeasurementReduction


@dataclass(frozen=True)
class ResolutionFunction(ABC):
    """Base interface for an L-direction CTR resolution function.

    The effective full width is

    ``delta_l_0 + delta_l_1 * abs(sin(gamma)) +
    delta_l_2 * abs(sin(delta))``.

    The full central HKL coordinates and angle record are accepted by
    :meth:`width` so resolution models can depend on diffractometer angles.

    :param float delta_l_0:
        Constant L-width contribution in r.l.u.
    :param float delta_l_1:
        Gamma-dependent L-width contribution in r.l.u.
    :param float delta_l_2:
        Delta-dependent L-width contribution in r.l.u.
    """

    delta_l_0: float
    delta_l_1: float = 0.0
    delta_l_2: float = 0.0

    def __post_init__(self):
        for name, value in (
            ("delta_l_0", self.delta_l_0),
            ("delta_l_1", self.delta_l_1),
            ("delta_l_2", self.delta_l_2),
        ):
            if not np.isscalar(value) or not np.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be a finite, nonnegative scalar")

    def width(self, h, k, l, angles=None):  # noqa: E741
        """Return the effective L width at the supplied CTR points.

        :param numpy.ndarray h:
            H coordinates in r.l.u.
        :param numpy.ndarray k:
            K coordinates in r.l.u.
        :param numpy.ndarray l:
            L coordinates in r.l.u.
        :param angles:
            Optional structured angle array containing the enabled ``gamma``
            and/or ``delta`` fields in rad.
        :returns:
            Effective L widths in r.l.u., broadcast to the HKL shape.
        :rtype: numpy.ndarray
        :raises ValueError:
            If HKL shapes are incompatible or gamma is required but missing.
        """
        try:
            h_array, k_array, l_array = np.broadcast_arrays(
                np.asarray(h, dtype=np.float64),
                np.asarray(k, dtype=np.float64),
                np.asarray(l, dtype=np.float64),
            )
        except ValueError as exc:
            raise ValueError(
                "H, K, and L coordinates must be broadcast-compatible"
            ) from exc

        widths = np.full(l_array.shape, float(self.delta_l_0), dtype=np.float64)
        if self.delta_l_1 == 0.0 and self.delta_l_2 == 0.0:
            return widths

        if angles is None:
            raise ValueError(
                "Angle-dependent resolution requires CTR angle records; "
                "call CTRCollection.calcAnglesZmode first"
            )

        for coefficient, angle_name in (
            (self.delta_l_1, "gamma"),
            (self.delta_l_2, "delta"),
        ):
            if coefficient == 0.0:
                continue
            try:
                angle = np.asarray(angles[angle_name], dtype=np.float64)
            except (KeyError, TypeError, ValueError, IndexError) as exc:
                raise ValueError(
                    "Angle-dependent resolution requires an angle field named "
                    f"'{angle_name}'"
                ) from exc
            try:
                angle = np.broadcast_to(angle, l_array.shape)
            except ValueError as exc:
                raise ValueError(
                    f"{angle_name.capitalize()} angles must have the same shape "
                    "as the CTR points"
                ) from exc
            if not np.all(np.isfinite(angle)):
                raise ValueError(
                    f"{angle_name.capitalize()} angles must be finite and "
                    "expressed in radians"
                )
            widths += coefficient * np.abs(np.sin(angle))
        return widths

    @abstractmethod
    def weights(self, offsets, width):
        """Return unnormalized kernel weights for L offsets in r.l.u."""

    @abstractmethod
    def quadrature(self, width, order):
        """Return normalized L quadrature offsets and weights.

        :param numpy.ndarray width:
            Effective widths in r.l.u.
        :param int order:
            Quadrature order.
        :returns:
            Tuple of offsets shaped ``width.shape + (order,)`` and normalized
            one-dimensional quadrature weights.
        :rtype: tuple[numpy.ndarray, numpy.ndarray]
        """


class BoxResolution(ResolutionFunction):
    """Box resolution whose ``DeltaL`` is the full support width."""

    def weights(self, offsets, width):
        """Return one inside ``[-DeltaL / 2, DeltaL / 2]`` and zero outside."""
        offsets = np.asarray(offsets, dtype=np.float64)
        width = np.asarray(width, dtype=np.float64)
        return (np.abs(offsets) <= width / 2.0).astype(np.float64)

    def quadrature(self, width, order):
        """Return Gauss-Legendre quadrature for the normalized box."""
        nodes, weights = leggauss(order)
        width = np.asarray(width, dtype=np.float64)
        offsets = width[..., np.newaxis] * nodes / 2.0
        return offsets, weights / 2.0


class GaussianResolution(ResolutionFunction):
    """Gaussian resolution whose ``DeltaL`` is its FWHM."""

    def weights(self, offsets, width):
        """Return Gaussian weights with the configured full width at half maximum."""
        offsets, width = np.broadcast_arrays(
            np.asarray(offsets, dtype=np.float64),
            np.asarray(width, dtype=np.float64),
        )
        result = np.zeros(offsets.shape, dtype=np.float64)
        nonzero = width > 0.0
        result[nonzero] = np.exp(
            -4.0
            * np.log(2.0)
            * (offsets[nonzero] / width[nonzero]) ** 2
        )
        result[~nonzero] = offsets[~nonzero] == 0.0
        return result

    def quadrature(self, width, order):
        """Return Gauss-Hermite quadrature for the normalized Gaussian."""
        nodes, weights = hermgauss(order)
        width = np.asarray(width, dtype=np.float64)
        sigma = width / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        offsets = np.sqrt(2.0) * sigma[..., np.newaxis] * nodes
        return offsets, weights / np.sqrt(np.pi)


def _validate_ctr_arrays(ctr, require_amplitudes):
    l_array = np.asarray(ctr.l, dtype=np.float64)
    h_array = np.asarray(ctr.harr, dtype=np.float64)
    k_array = np.asarray(ctr.karr, dtype=np.float64)
    if (
        l_array.ndim != 1
        or h_array.shape != l_array.shape
        or k_array.shape != l_array.shape
    ):
        raise ValueError(
            f"{ctr!r} must contain one-dimensional, equal-length HKL arrays"
        )
    if (
        not np.all(np.isfinite(h_array))
        or not np.all(np.isfinite(k_array))
        or not np.all(np.isfinite(l_array))
    ):
        raise ValueError(f"{ctr!r} contains non-finite HKL coordinates")
    if np.unique(l_array).size != l_array.size:
        raise ValueError(f"{ctr!r} contains duplicate L coordinates")

    if not require_amplitudes:
        return h_array, k_array, l_array, None

    if not np.isrealobj(ctr.sfI):
        raise ValueError(f"{ctr!r} amplitudes must be real")
    amplitude = np.asarray(ctr.sfI, dtype=np.float64)
    if amplitude.shape != l_array.shape:
        raise ValueError(f"{ctr!r} amplitudes must have the same shape as L")
    if not np.all(np.isfinite(amplitude)) or np.any(amplitude < 0.0):
        raise ValueError(f"{ctr!r} amplitudes must be finite and nonnegative")
    return h_array, k_array, l_array, amplitude


def _validate_widths(widths, shape):
    widths = np.asarray(widths, dtype=np.float64)
    try:
        widths = np.broadcast_to(widths, shape)
    except ValueError as exc:
        raise ValueError("Resolution widths must match the CTR point shape") from exc
    if not np.all(np.isfinite(widths)) or np.any(widths < 0.0):
        raise ValueError("Effective resolution widths must be finite and nonnegative")
    return widths


def _validate_intensity_coordinates(h, k, l):  # noqa: E741
    try:
        h_array, k_array, l_array = np.broadcast_arrays(
            np.asarray(h, dtype=np.float64),
            np.asarray(k, dtype=np.float64),
            np.asarray(l, dtype=np.float64),
        )
    except ValueError as exc:
        raise ValueError(
            "H, K, and L coordinates must be broadcast-compatible"
        ) from exc
    if l_array.ndim == 0:
        h_array = h_array.reshape(1)
        k_array = k_array.reshape(1)
        l_array = l_array.reshape(1)
    if l_array.ndim != 1 or l_array.size == 0:
        raise ValueError(
            "Resolution coordinates must be nonempty one-dimensional arrays"
        )
    if not all(np.all(np.isfinite(array)) for array in (h_array, k_array, l_array)):
        raise ValueError("Resolution coordinates must be finite")
    if np.unique(l_array).size != l_array.size:
        raise ValueError("Resolution coordinates must not contain duplicate L values")
    return h_array, k_array, l_array


def _validate_intensity(intensity, shape, source="intensity"):
    if not np.isrealobj(intensity):
        raise ValueError(f"{source} must be real")
    intensity = np.asarray(intensity, dtype=np.float64)
    try:
        intensity = np.broadcast_to(intensity, shape)
    except ValueError as exc:
        raise ValueError(f"{source} returned an incompatible array shape") from exc
    if not np.all(np.isfinite(intensity)) or np.any(intensity < 0.0):
        raise ValueError(f"{source} must be finite and nonnegative")
    return intensity


def _validate_quadrature_order(quadrature_order):
    if (
        not isinstance(quadrature_order, int | np.integer)
        or isinstance(quadrature_order, bool | np.bool_)
        or quadrature_order <= 0
        or quadrature_order % 2 == 0
    ):
        raise ValueError("quadrature_order must be a positive odd integer")
    return int(quadrature_order)


def fast_convolve_intensity(h, k, l, intensity, resolution, angles=None):  # noqa: E741
    """Convolve field intensity on an existing irregular L grid.

    No interpolation or equidistant resampling is performed. Local composite
    trapezoidal weights account for unequal L spacing, and each kernel is
    normalized over the samples available at the rod boundaries. H and K are
    reference-frame coordinates in r.l.u.; ``angles`` contains optional
    central diffractometer records in rad.

    :param h: H coordinates in r.l.u., scalar or one-dimensional.
    :param k: K coordinates in r.l.u., scalar or one-dimensional.
    :param l: Unique L coordinates in r.l.u., scalar or one-dimensional.
    :param intensity: Finite, nonnegative field intensity at the central grid.
    :param ResolutionFunction resolution: L-resolution model.
    :param angles: Optional central structured angle records in rad.
    :returns: Convolved field intensity in the original point order.
    :rtype: numpy.ndarray
    """
    if not isinstance(resolution, ResolutionFunction):
        raise TypeError("resolution must be a ResolutionFunction")
    h_array, k_array, l_array = _validate_intensity_coordinates(h, k, l)
    intensity = _validate_intensity(intensity, l_array.shape)
    widths = _validate_widths(
        resolution.width(h_array, k_array, l_array, angles), l_array.shape
    )

    if l_array.size < 2:
        return np.ascontiguousarray(intensity)

    order = np.argsort(l_array)
    l_sorted = l_array[order]
    intensity_sorted = intensity[order]
    widths_sorted = widths[order]

    integration_weights = np.empty_like(l_sorted)
    integration_weights[0] = (l_sorted[1] - l_sorted[0]) / 2.0
    integration_weights[-1] = (l_sorted[-1] - l_sorted[-2]) / 2.0
    integration_weights[1:-1] = (l_sorted[2:] - l_sorted[:-2]) / 2.0

    convolved_sorted = np.empty_like(intensity_sorted)
    for index, (center, width) in enumerate(zip(l_sorted, widths_sorted)):
        if width == 0.0:
            convolved_sorted[index] = intensity_sorted[index]
            continue
        kernel_weights = resolution.weights(l_sorted - center, width)
        combined_weights = kernel_weights * integration_weights
        normalization = np.sum(combined_weights)
        if not np.isfinite(normalization) or normalization <= 0.0:
            raise ValueError(
                f"Resolution kernel has no support at L={center:g} r.l.u."
            )
        convolved_sorted[index] = (
            np.sum(combined_weights * intensity_sorted) / normalization
        )

    convolved = np.empty_like(convolved_sorted)
    convolved[order] = np.maximum(convolved_sorted, 0.0)
    return np.ascontiguousarray(convolved)


def sample_intensity(
    h,
    k,
    l,  # noqa: E741
    intensity,
    resolution,
    angles=None,
    quadrature_order=25,
):
    """Integrate a field-intensity callable around central CTR points.

    H and K remain fixed while L is displaced by the resolution quadrature.
    The callable receives flattened H, K, and L arrays in central-point-major,
    quadrature-point-minor order and must return broadcast-compatible finite,
    nonnegative field intensity. Coordinates are in r.l.u.; optional central
    diffractometer angle records are in rad.

    :param h: H coordinates in r.l.u., scalar or one-dimensional.
    :param k: K coordinates in r.l.u., scalar or one-dimensional.
    :param l: Unique L coordinates in r.l.u., scalar or one-dimensional.
    :param callable intensity: Callable ``intensity(h, k, l)``.
    :param ResolutionFunction resolution: L-resolution model.
    :param angles: Optional central structured angle records in rad.
    :param int quadrature_order: Positive odd quadrature order.
    :returns: Integrated field intensity in central-point order.
    :rtype: numpy.ndarray
    """
    if not callable(intensity):
        raise TypeError("intensity must be callable")
    if not isinstance(resolution, ResolutionFunction):
        raise TypeError("resolution must be a ResolutionFunction")
    quadrature_order = _validate_quadrature_order(quadrature_order)
    h_array, k_array, l_array = _validate_intensity_coordinates(h, k, l)
    widths = _validate_widths(
        resolution.width(h_array, k_array, l_array, angles), l_array.shape
    )
    offsets, integration_weights = resolution.quadrature(
        widths, quadrature_order
    )
    sample_shape = offsets.shape
    h_samples = np.broadcast_to(h_array[..., np.newaxis], sample_shape)
    k_samples = np.broadcast_to(k_array[..., np.newaxis], sample_shape)
    l_samples = l_array[..., np.newaxis] + offsets
    sampled = intensity(
        h_samples.ravel(), k_samples.ravel(), l_samples.ravel()
    )
    sampled = _validate_intensity(
        sampled, (l_samples.size,), "intensity callable"
    ).reshape(sample_shape)
    integrated = np.sum(sampled * integration_weights, axis=-1)
    return np.ascontiguousarray(np.maximum(integrated, 0.0))


def _new_collection(ctrs, amplitudes, *, preserve_measurement_metadata=False):
    result = CTRCollection(name=ctrs.name)
    result.plotsett = copy.deepcopy(ctrs.plotsett)
    result.plotkeyargs = copy.deepcopy(ctrs.plotkeyargs)
    for source, amplitude in zip(ctrs, amplitudes):
        target = CTR(
            source.hk,
            np.copy(source.l),
            np.ascontiguousarray(amplitude),
            err=None,
            phi=None,
            name=source.name,
            reduction=(
                copy.deepcopy(source.reduction)
                if preserve_measurement_metadata
                else MeasurementReduction()
            ),
            scan_geometry=(
                copy.deepcopy(source.scan_geometry)
                if preserve_measurement_metadata
                else None
            ),
        )
        target.harr = np.copy(source.harr)
        target.karr = np.copy(source.karr)
        target.weight = source.weight
        if hasattr(source, "angles"):
            target.angles = copy.deepcopy(source.angles)
        result.append(target)
    return result


def fast_convolve(ctrs, resolution):
    """Convolve CTR intensity using only the existing irregular L points.

    No interpolation or equidistant resampling is performed. Local composite
    trapezoidal weights account for unequal L spacing, and each kernel is
    renormalized over the samples available at the rod boundaries.

    :param CTRCollection ctrs:
        CTR amplitudes to convolve. Amplitudes must be finite and nonnegative.
    :param ResolutionFunction resolution:
        Box or Gaussian L-resolution function.
    :returns:
        A new collection containing effective amplitudes
        ``sqrt(convolved(abs(F)**2))``.
    :rtype: CTRCollection
    """
    if not isinstance(ctrs, CTRCollection):
        raise TypeError("ctrs must be a CTRCollection")
    if not isinstance(resolution, ResolutionFunction):
        raise TypeError("resolution must be a ResolutionFunction")
    ctrs._require_structure_factors("structure-factor resolution convolution")

    convolved = []
    for ctr in ctrs:
        h_array, k_array, l_array, amplitude = _validate_ctr_arrays(
            ctr, require_amplitudes=True
        )
        angles = getattr(ctr, "angles", None)
        convolved_intensity = fast_convolve_intensity(
            h_array,
            k_array,
            l_array,
            amplitude**2,
            resolution,
            angles,
        )
        convolved.append(np.sqrt(convolved_intensity))

    return _new_collection(
        ctrs, convolved, preserve_measurement_metadata=True
    )


def sample_structure_factor(ctrs, crystal, resolution, quadrature_order=25):
    """Sample crystal intensity around every CTR point along L.

    H and K remain fixed. Box functions use Gauss-Legendre quadrature and
    Gaussian functions use Gauss-Hermite quadrature. The effective width is
    evaluated from each central point and its central angle record.

    :param CTRCollection ctrs:
        Collection supplying the requested HKL points and optional angles.
    :param crystal:
        Forward model providing ``F2(h, k, l)``, or a crystal-like object
        providing only ``F(h, k, l)``. ``F2`` is preferred, because a mixed
        state has a squared structure factor but no unique complex amplitude.
    :param ResolutionFunction resolution:
        Box or Gaussian L-resolution function.
    :param int quadrature_order:
        Positive odd number of deterministic quadrature points. Defaults to
        25.
    :returns:
        A new collection containing effective amplitudes
        ``sqrt(integrated(F2))``.
    :rtype: CTRCollection
    """
    if not isinstance(ctrs, CTRCollection):
        raise TypeError("ctrs must be a CTRCollection")
    if not isinstance(resolution, ResolutionFunction):
        raise TypeError("resolution must be a ResolutionFunction")
    squared = getattr(crystal, "F2", None)
    if not callable(squared) and not (
        hasattr(crystal, "F") and callable(crystal.F)
    ):
        raise TypeError(
            "crystal must provide a callable F2(h, k, l) or F(h, k, l) method"
        )
    quadrature_order = _validate_quadrature_order(quadrature_order)

    sampled = []
    for ctr in ctrs:
        h_array, k_array, l_array, _ = _validate_ctr_arrays(
            ctr, require_amplitudes=False
        )
        angles = getattr(ctr, "angles", None)

        def structure_factor_intensity(h_samples, k_samples, l_samples):
            if callable(squared):
                values = np.asarray(squared(h_samples, k_samples, l_samples))
                source = "crystal.F2"
            else:
                values = (
                    np.abs(
                        np.asarray(
                            crystal.F(h_samples, k_samples, l_samples)
                        )
                    )
                    ** 2
                )
                source = "crystal.F"
            try:
                intensity = np.broadcast_to(values, h_samples.shape)
            except ValueError as exc:
                raise ValueError(
                    f"{source} returned an incompatible array shape"
                ) from exc
            if not np.all(np.isfinite(intensity)):
                raise ValueError(
                    f"{source} returned non-finite structure factors"
                )
            return intensity

        integrated = sample_intensity(
            h_array,
            k_array,
            l_array,
            structure_factor_intensity,
            resolution,
            angles,
            quadrature_order,
        )
        sampled.append(np.sqrt(integrated))

    return _new_collection(ctrs, sampled)

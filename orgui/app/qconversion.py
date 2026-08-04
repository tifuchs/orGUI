"""Reciprocal-space conversion of detector images for the Q-plot.

.. warning::

   **Experimental application code.** This module backs the experimental
   ``Q-plot`` action of the GUI. It is deliberately kept out of
   :mod:`orgui.datautils.xrayutils` so that it is not mistaken for the
   production reciprocal-space code. The authoritative per-pixel calculation is
   :meth:`~orgui.datautils.xrayutils.HKLVlieg.VliegAngles.QAlpha`; everything
   here only exists to display a whole image on a regular grid, and its
   conventions may still change.

The image is rebinned by pyFAI's ``FiberIntegrator``. pyFAI's built-in fiber
units cannot express orGUI's geometry, for two independent reasons:

* ``sample_orientation`` is the EXIF dihedral group acting on the image array,
  so it only spans quarter turns. orGUI's azimuthal reference is continuous and
  is therefore passed through ``tilt_angle``, which rotates about the incident
  beam. ``sample_orientation`` stays at the identity.
* pyFAI composes ``incident_angle`` and ``tilt_angle`` extrinsically as
  ``R_x(-tilt) @ R_y(incidence)``, i.e. about fixed axes. In orGUI the azimuth
  rotates the alpha rotation axis itself, so the incidence rotation has to be
  applied about the *already rotated* axis, ``R_y(incidence) @ R_x(-tilt)``.
  The two agree only when one of the angles vanishes.

This module therefore builds its own ``pyFAI.units.UnitFiber`` objects; pyFAI
still performs the rebinning.

Frames
------

``QAlpha`` returns the momentum transfer in the alpha frame, the frame of the
sample surface. The other frames are reached with the Vlieg rotation matrices,
following ``hkl = UB^-1 PHI^-1 CHI^-1 OMEGA^-1 Q_alpha``. ``Q_cryst``
additionally undoes the orientation matrix, so it holds the Cartesian
reciprocal-space vector of the crystal, ``B`` times ``hkl``. For every frame the
out-of-plane component is the ``z`` axis of that frame and the in-plane
component is the radial component in its ``xy`` plane.

Performance
-----------

Detector images are large, so the whole conversion is collapsed into a single
affine relation before touching any pixel data. With
``inv = 1 / |(x, y, z)|`` every component reduces to

.. code-block:: text

    q_j = k * ( (G[j, 0] * x + G[j, 1] * y + G[j, 2] * z) * inv - c[j] )

where ``G`` and ``c`` absorb the sample orientation map, the beam and incidence
rotations, the axis relabelling and the frame rotation. The compiled kernel in
``orgui/app/cpp/qconversion_cpp.cpp`` evaluates both displayed quantities in one
pass; a plain numpy fallback is used when the extension is unavailable.

``numexpr`` is deliberately **not** used here: version 2.11.0 returns wrong
results for a small fraction of multi-threaded evaluations of this kind of
expression, which is not acceptable for a coordinate transform.
"""

import importlib
import importlib.util
from pathlib import Path

import numpy as np

from ..datautils.xrayutils import HKLVlieg

try:
    from pyFAI import units as pyFAI_units
    from pyFAI.integrator.fiber import FiberIntegrator

    HAS_FIBER_INTEGRATOR = True
except ImportError:  # pyFAI < 2025.1
    HAS_FIBER_INTEGRATOR = False


# orGUI's azimuthal reference and pyFAI's fiber convention differ by a quarter
# turn about the incident beam.
AZIMUTH_OFFSET = np.pi / 2

# Kept at the identity on purpose, see the module docstring.
SAMPLE_ORIENTATION = 1

# Relative slack allowed when inverting the projection, so that a coordinate
# right on the horizon of the accessible region is not rejected by round-off.
_EWALD_TOLERANCE = 1e-9

FRAMES = ("Q_alpha", "Q_lab", "Q_omega", "Q_chi", "Q_phi", "Q_cryst")
DEFAULT_FRAME = "Q_alpha"

# Frames reached through the omega rotation. They are only defined for a single
# image, because a maximum or sum image combines many omega angles.
FRAMES_REQUIRING_OMEGA = ("Q_omega", "Q_chi", "Q_phi", "Q_cryst")

# Frames that additionally need the orientation matrix U.
FRAMES_REQUIRING_U = ("Q_cryst",)

FRAME_LABELS = {
    "Q_alpha": "Q alpha (surface)",
    "Q_lab": "Q lab",
    "Q_omega": "Q omega",
    "Q_chi": "Q chi",
    "Q_phi": "Q phi",
    "Q_cryst": "Q crystal",
}

# (beam, horz, vert) -> orGUI's (Qx, Qy, Qz)
_AXIS_RELABEL = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

_ORIENTATION_BASIS = {
    "x": (1.0, 0.0),
    "y": (0.0, 1.0),
    "(-x)": (-1.0, 0.0),
    "(-y)": (0.0, -1.0),
}

_KERNEL = None
_KERNEL_LOADED = False


def _import_kernel():
    """Import the packaged extension or an in-tree Meson build artifact."""
    try:
        return importlib.import_module("orgui.app._qconversion_cpp")
    except ModuleNotFoundError as package_error:
        repo_root = Path(__file__).resolve().parents[2]
        candidates = sorted(repo_root.glob("build/**/_qconversion_cpp*.pyd"))
        candidates += sorted(repo_root.glob("build/**/_qconversion_cpp*.so"))
        if not candidates:
            raise package_error
        spec = importlib.util.spec_from_file_location(
            "_qconversion_cpp", candidates[-1]
        )
        if spec is None or spec.loader is None:
            raise package_error
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module


def kernel():
    """Return the compiled kernel, or ``None`` when it is unavailable."""
    global _KERNEL, _KERNEL_LOADED
    if not _KERNEL_LOADED:
        try:
            _KERNEL = _import_kernel()
        except Exception:
            _KERNEL = None
        _KERNEL_LOADED = True
    return _KERNEL


def hasKernel():
    """Whether the compiled conversion kernel is available."""
    return kernel() is not None


def frameShortName(frame):
    """Return the frame name without the leading ``Q_``, for compact labels.

    :param str frame: one of :data:`FRAMES`.
    :returns: for example ``"alpha"`` for ``"Q_alpha"``.
    """
    if frame not in FRAMES:
        raise ValueError(f"Unknown frame {frame!r}, expected one of {FRAMES}")
    return frame[len("Q_") :]


def tiltAngleFromAzimuth(azimuth):
    """Return the pyFAI ``tilt_angle`` reproducing an orGUI azimuthal reference.

    :param float azimuth: azimuthal reference of the detector calibration, rad.
    :returns: tilt angle in rad.
    """
    return -(azimuth + AZIMUTH_OFFSET)


def frameMatrix(frame, alpha=0.0, omega=0.0, chi=0.0, phi=0.0, U=None):
    """Rotation taking a vector from the alpha frame into ``frame``.

    :param str frame: one of :data:`FRAMES`.
    :param float alpha: incidence angle in rad.
    :param float omega, chi, phi: sample circle angles in rad.
    :param U: orientation matrix, required for :data:`FRAMES_REQUIRING_U`.
    :returns: 3x3 rotation matrix.
    """
    if frame not in FRAMES:
        raise ValueError(f"Unknown frame {frame!r}, expected one of {FRAMES}")
    if frame in FRAMES_REQUIRING_U and U is None:
        raise ValueError(f"The {frame} frame requires the orientation matrix U")

    ALPHA, _, _, OMEGA, CHI, PHI = HKLVlieg.createVliegMatrices(
        [alpha, None, None, omega, chi, phi]
    )
    if frame == "Q_alpha":
        return np.identity(3)
    if frame == "Q_lab":
        return np.asarray(ALPHA)

    matrix = np.asarray(OMEGA).T
    if frame == "Q_omega":
        return matrix
    matrix = np.asarray(CHI).T @ matrix
    if frame == "Q_chi":
        return matrix
    matrix = np.asarray(PHI).T @ matrix
    if frame == "Q_phi":
        return matrix
    # Q_cryst: undo the orientation matrix as well, so that the result is the
    # Cartesian reciprocal-space vector of the crystal, i.e. B times hkl.
    return np.linalg.inv(np.asarray(U)) @ matrix


def _orientationMatrix(sample_orientation):
    """Embed the EXIF orientation map as a 3x3 matrix on ``(x, y, z)``."""
    mapping = pyFAI_units.MAPS_SAMPLE_ORIENTATION[sample_orientation]
    matrix = np.zeros((3, 3))
    matrix[0, 2] = 1.0  # the beam component picks z
    matrix[1, :2] = _ORIENTATION_BASIS[mapping["x"]]
    matrix[2, :2] = _ORIENTATION_BASIS[mapping["y"]]
    return matrix


def conversionCoefficients(
    frame,
    incident_angle=0.0,
    tilt_angle=0.0,
    sample_orientation=SAMPLE_ORIENTATION,
    omega=0.0,
    chi=0.0,
    phi=0.0,
    U=None,
):
    """Collapse the whole conversion into ``(G, c)``.

    The momentum transfer of a pixel at ``(x, y, z)`` is then
    ``k * (G @ (x, y, z) / |(x, y, z)| - c)`` in inverse nanometre, with
    ``k = 2 pi / wavelength``.

    :returns: ``(G, c)`` with shapes ``(3, 3)`` and ``(3,)``.
    """
    incidence = np.array(
        [
            [np.cos(incident_angle), 0.0, np.sin(incident_angle)],
            [0.0, 1.0, 0.0],
            [-np.sin(incident_angle), 0.0, np.cos(incident_angle)],
        ]
    )
    tilt = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(tilt_angle), np.sin(tilt_angle)],
            [0.0, -np.sin(tilt_angle), np.cos(tilt_angle)],
        ]
    )
    combined = (
        frameMatrix(frame, alpha=incident_angle, omega=omega, chi=chi, phi=phi, U=U)
        @ _AXIS_RELABEL
        @ incidence
        @ tilt
    )
    return (
        np.ascontiguousarray(combined @ _orientationMatrix(sample_orientation)),
        np.ascontiguousarray(combined[:, 0]),
    )


def _ipOopNumpy(x, y, z, G, c, k):
    """Fallback used when the compiled kernel is unavailable."""
    scratch = np.empty_like(x)
    accumulator = np.empty_like(x)

    def component(index):
        np.multiply(x, G[index, 0], out=accumulator)
        np.multiply(y, G[index, 1], out=scratch)
        np.add(accumulator, scratch, out=accumulator)
        np.multiply(z, G[index, 2], out=scratch)
        np.add(accumulator, scratch, out=accumulator)
        np.multiply(accumulator, inverse_norm, out=accumulator)
        np.subtract(accumulator, c[index], out=accumulator)
        np.multiply(accumulator, k, out=accumulator)
        return accumulator.copy()

    inverse_norm = np.empty_like(x)
    np.multiply(x, x, out=inverse_norm)
    np.multiply(y, y, out=scratch)
    np.add(inverse_norm, scratch, out=inverse_norm)
    np.multiply(z, z, out=scratch)
    np.add(inverse_norm, scratch, out=inverse_norm)
    np.sqrt(inverse_norm, out=inverse_norm)
    np.reciprocal(inverse_norm, out=inverse_norm)

    q_x = component(0)
    q_y = component(1)
    q_oop = component(2)
    np.multiply(q_x, q_x, out=scratch)
    np.multiply(q_y, q_y, out=q_y)
    np.add(scratch, q_y, out=scratch)
    np.sqrt(scratch, out=scratch)
    # signed by the in-plane direction of orGUI's delta angle
    return np.copysign(scratch, q_x), q_oop


class _ConversionCache:
    """Remember the last conversion.

    pyFAI evaluates the in-plane and the out-of-plane unit separately but hands
    both the very same position arrays, so the image is converted once instead
    of twice.
    """

    def __init__(self):
        self._arrays = None
        self._signature = None
        self._value = None

    def get(self, x, y, z, G, c, k):
        signature = (G.tobytes(), c.tobytes(), float(k))
        if self._arrays is not None and self._signature == signature:
            cached_x, cached_y, cached_z = self._arrays
            if cached_x is x and cached_y is y and cached_z is z:
                return self._value

        compiled = kernel()
        if compiled is not None:
            value = compiled.q_ip_oop(x, y, z, G, c, float(k))
        else:
            value = _ipOopNumpy(
                np.ascontiguousarray(x, dtype=np.float64),
                np.ascontiguousarray(y, dtype=np.float64),
                np.ascontiguousarray(z, dtype=np.float64),
                G,
                c,
                k,
            )
        self._arrays = (x, y, z)
        self._signature = signature
        self._value = value
        return value

    def clear(self):
        self._arrays = None
        self._signature = None
        self._value = None


_CONVERSION = _ConversionCache()


class _IntegratorCache:
    """Reuse the pyFAI integrator across conversions of the same geometry.

    A fresh ``FiberIntegrator`` starts with empty grazing-incidence parameters,
    so passing any non-zero incidence or tilt angle makes pyFAI reset and
    recompute all of its cached arrays. Keeping the integrator avoids that on
    every image change.

    The cache key contains the collapsed conversion coefficients, so a hit
    means the unit equations are byte for byte identical. That matters because
    pyFAI caches the per-pixel unit arrays under the unit *name*: reusing an
    integrator with a same-named but differently parameterised unit would
    silently return the previous image's coordinates.
    """

    def __init__(self):
        self._key = None
        self._detector = None
        self._integrator = None

    def get(self, detectorCal, G, c):
        key = (
            float(detectorCal.dist),
            float(detectorCal.poni1),
            float(detectorCal.poni2),
            float(detectorCal.rot1),
            float(detectorCal.rot2),
            float(detectorCal.rot3),
            float(detectorCal.wavelength),
            G.tobytes(),
            c.tobytes(),
        )
        if self._key == key and self._detector is detectorCal.detector:
            return self._integrator

        integrator = FiberIntegrator(
            dist=detectorCal.dist,
            poni1=detectorCal.poni1,
            poni2=detectorCal.poni2,
            wavelength=detectorCal.wavelength,
            rot1=detectorCal.rot1,
            rot2=detectorCal.rot2,
            rot3=detectorCal.rot3,
            detector=detectorCal.detector,
        )
        self._key = key
        self._detector = detectorCal.detector
        self._integrator = integrator
        return integrator

    def clear(self):
        self._key = None
        self._detector = None
        self._integrator = None


_INTEGRATOR = _IntegratorCache()


def qIpOop(
    frame,
    x,
    y,
    z,
    wavelength,
    incident_angle=0.0,
    tilt_angle=0.0,
    sample_orientation=SAMPLE_ORIENTATION,
    omega=0.0,
    chi=0.0,
    phi=0.0,
    U=None,
):
    """In-plane and out-of-plane momentum transfer, in inverse nanometre."""
    G, c = conversionCoefficients(
        frame,
        incident_angle=incident_angle,
        tilt_angle=tilt_angle,
        sample_orientation=sample_orientation,
        omega=omega,
        chi=chi,
        phi=phi,
        U=U,
    )
    return _CONVERSION.get(x, y, z, G, c, 2.0e-9 * np.pi / wavelength)


def qComponents(
    frame,
    x,
    y,
    z,
    wavelength,
    incident_angle=0.0,
    tilt_angle=0.0,
    sample_orientation=SAMPLE_ORIENTATION,
    omega=0.0,
    chi=0.0,
    phi=0.0,
    U=None,
):
    """All three Cartesian components expressed in ``frame``.

    This is the readable reference implementation; the Q-plot itself only needs
    the two quantities returned by :func:`qIpOop`.

    :returns: ``(Qx, Qy, Qz)`` of ``frame`` in inverse nanometre.
    """
    G, c = conversionCoefficients(
        frame,
        incident_angle=incident_angle,
        tilt_angle=tilt_angle,
        sample_orientation=sample_orientation,
        omega=omega,
        chi=chi,
        phi=phi,
        U=U,
    )
    k = 2.0e-9 * np.pi / wavelength
    inverse_norm = 1.0 / np.sqrt(x * x + y * y + z * z)
    return tuple(
        k * ((G[j, 0] * x + G[j, 1] * y + G[j, 2] * z) * inverse_norm - c[j])
        for j in range(3)
    )


def beamDirection(frame, alpha=0.0, omega=0.0, chi=0.0, phi=0.0, U=None):
    """Unit vector along the incident beam, expressed in ``frame``.

    In the alpha frame the beam is ``(0, cos alpha, -sin alpha)``, which is the
    vector :meth:`~orgui.datautils.xrayutils.HKLVlieg.VliegAngles.QAlpha`
    subtracts from the outgoing direction. It is the same vector as the ``c``
    returned by :func:`conversionCoefficients`, which is the beam column of the
    collapsed conversion.

    :param str frame: one of :data:`FRAMES`.
    :param float alpha: incidence angle in rad.
    :param float omega, chi, phi: sample circle angles in rad.
    :param U: orientation matrix, required for :data:`FRAMES_REQUIRING_U`.
    :returns: unit vector of shape ``(3,)``.
    """
    rotation = frameMatrix(frame, alpha=alpha, omega=omega, chi=chi, phi=phi, U=U)
    return rotation @ np.array([0.0, np.cos(alpha), -np.sin(alpha)])


def qFromIpOop(
    frame,
    q_ip,
    q_oop,
    k,
    alpha=0.0,
    omega=0.0,
    chi=0.0,
    phi=0.0,
    U=None,
):
    """Momentum transfer vectors that project onto one Q-plot coordinate.

    This inverts the projection performed by :func:`qIpOop`, which keeps the
    out-of-plane component and collapses the two in-plane components into their
    radial distance, signed by ``Qx``. The missing information is recovered from
    the Ewald condition: with ``q = Q / k`` and the beam direction ``c`` of
    :func:`beamDirection`, the outgoing direction ``q + c`` is a unit vector, so

    .. code-block:: text

        |q|^2 + 2 q.c = 0

    which is *linear* in ``qx`` and ``qy`` once ``qz`` and ``|q|`` are known.
    Intersecting that line with the circle ``qx^2 + qy^2 = (q_ip / k)^2`` leaves
    at most two solutions, of which only those with ``sign(Qx) == sign(q_ip)``
    reproduce the sign convention of :func:`qIpOop`.

    In the alpha frame the beam has no ``x`` component, the line is parallel to
    the ``x`` axis and the sign condition always leaves a single solution. In
    the rotated frames the beam generally does have one, and both solutions can
    carry the same sign; the caller has to decide between them, for example by
    keeping the one that falls onto the detector.

    ``k``, ``q_ip`` and ``q_oop`` must share the same unit, and the result is
    returned in that unit. Note that :func:`qIpOop` works in inverse nanometre
    while orGUI's own reciprocal space is in inverse Angstrom.

    :param str frame: one of :data:`FRAMES`.
    :param float q_ip: in-plane coordinate, signed as in :func:`qIpOop`.
    :param float q_oop: out-of-plane coordinate.
    :param float k: magnitude of the wave vector.
    :param float alpha: incidence angle in rad.
    :param float omega, chi, phi: sample circle angles in rad.
    :param U: orientation matrix, required for :data:`FRAMES_REQUIRING_U`.
    :returns:
        Tuple of zero, one or two momentum transfer vectors of shape ``(3,)``,
        expressed in ``frame``. Empty when the coordinate is outside the part
        of reciprocal space the Ewald sphere can reach.
    :rtype: tuple[numpy.ndarray, ...]
    """
    q_ip = float(q_ip)
    q_oop = float(q_oop)
    k = float(k)
    if not np.isfinite([q_ip, q_oop, k]).all() or k == 0.0:
        return ()

    c = beamDirection(frame, alpha=alpha, omega=omega, chi=chi, phi=phi, U=U)

    # everything below is in units of k, so that the Ewald sphere is the unit
    # sphere around the tip of the incident wave vector
    radius = abs(q_ip) / k
    q_z = q_oop / k
    squared_length = radius * radius + q_z * q_z

    # c[0] * q_x + c[1] * q_y = offset, from the Ewald condition
    offset = -0.5 * squared_length - c[2] * q_z
    normal_length = c[0] * c[0] + c[1] * c[1]
    if normal_length <= 0.0:
        # the beam is along z, so the condition no longer constrains the plane
        return ()

    # distance from the origin to the line, compared with the circle radius
    half_chord_squared = radius * radius - offset * offset / normal_length
    if half_chord_squared < 0.0:
        # tolerate round-off on the horizon of the accessible region
        if half_chord_squared < -_EWALD_TOLERANCE * max(radius * radius, 1.0):
            return ()
        half_chord_squared = 0.0
    half_chord = np.sqrt(half_chord_squared / normal_length)

    # foot of the perpendicular from the origin onto the line
    foot_x = c[0] * offset / normal_length
    foot_y = c[1] * offset / normal_length

    solutions = []
    for sign in (1.0, -1.0):
        q_x = foot_x - sign * c[1] * half_chord
        q_y = foot_y + sign * c[0] * half_chord
        # np.copysign in qIpOop signs the radius by Qx, so a solution is only
        # valid when its Qx carries the sign of the requested in-plane value
        if q_ip != 0.0 and np.sign(q_x) != np.sign(q_ip):
            continue
        solutions.append(k * np.array([q_x, q_y, q_z]))
        if half_chord == 0.0:
            # the line touches the circle, both signs give the same point
            break
    return tuple(solutions)


def detectorAngles(
    Q, k, frame=DEFAULT_FRAME, alpha=0.0, omega=0.0, chi=0.0, phi=0.0, U=None
):
    """Detector angles producing a momentum transfer expressed in ``frame``.

    The momentum transfer is rotated back into the alpha frame, where
    ``Q / k + c`` is the outgoing unit direction
    ``(sin delta cos gamma, cos delta cos gamma, sin gamma)`` of
    :meth:`~orgui.datautils.xrayutils.HKLVlieg.VliegAngles.QAlpha`.

    :param Q: momentum transfer of ``frame``, shape ``(3,)``.
    :param float k: magnitude of the wave vector, in the unit of ``Q``.
    :param str frame: one of :data:`FRAMES`.
    :param float alpha: incidence angle in rad.
    :param float omega, chi, phi: sample circle angles in rad.
    :param U: orientation matrix, required for :data:`FRAMES_REQUIRING_U`.
    :returns: ``(delta, gamma)`` in rad, or ``None`` when ``Q`` is not on the
        Ewald sphere and the direction cannot be normalised.
    :rtype: tuple[float, float] or None
    """
    rotation = frameMatrix(frame, alpha=alpha, omega=omega, chi=chi, phi=phi, U=U)
    q_alpha = np.linalg.inv(rotation) @ np.asarray(Q, dtype=np.float64)
    direction = q_alpha / float(k) + np.array([0.0, np.cos(alpha), -np.sin(alpha)])

    norm = np.linalg.norm(direction)
    if not np.isfinite(norm) or norm == 0.0:
        return None
    # the caller may have built Q from a rounded plot coordinate, so normalise
    # instead of insisting on an exactly unit direction
    direction = direction / norm

    gamma = np.arcsin(np.clip(direction[2], -1.0, 1.0))
    delta = np.arctan2(direction[0], direction[1])
    return float(delta), float(gamma)


def fiberUnits(frame=DEFAULT_FRAME, omega=0.0, chi=0.0, phi=0.0, U=None):
    """Build the in-plane and out-of-plane pyFAI units for ``frame``.

    :returns: ``(unit_ip, unit_oop)`` as ``pyFAI.units.UnitFiber``.
    """
    if not HAS_FIBER_INTEGRATOR:
        raise RuntimeError("pyFAI >= 2025.1 with FiberIntegrator is required")
    if frame not in FRAMES:
        raise ValueError(f"Unknown frame {frame!r}, expected one of {FRAMES}")
    if frame in FRAMES_REQUIRING_U and U is None:
        raise ValueError(f"The {frame} frame requires the orientation matrix U")

    def evaluate(index, x, y, z, wavelength, incident_angle, tilt, orientation):
        return qIpOop(
            frame,
            x,
            y,
            z,
            wavelength,
            incident_angle=incident_angle,
            tilt_angle=tilt,
            sample_orientation=orientation,
            omega=omega,
            chi=chi,
            phi=phi,
            U=U,
        )[index]

    def equation_ip(
        x,
        y,
        z,
        wavelength,
        incident_angle=0.0,
        tilt_angle=0.0,
        sample_orientation=SAMPLE_ORIENTATION,
    ):
        return evaluate(
            0, x, y, z, wavelength, incident_angle, tilt_angle, sample_orientation
        )

    def equation_oop(
        x,
        y,
        z,
        wavelength,
        incident_angle=0.0,
        tilt_angle=0.0,
        sample_orientation=SAMPLE_ORIENTATION,
    ):
        return evaluate(
            1, x, y, z, wavelength, incident_angle, tilt_angle, sample_orientation
        )

    label = FRAME_LABELS[frame]
    unit_ip = pyFAI_units.UnitFiber(
        f"qip_{frame}_A^-1",
        scale=0.1,
        label=rf"{label} $q_\parallel$ ($\AA^{{-1}}$)",
        equation=equation_ip,
        short_name=f"qip_{frame}",
        unit_symbol=r"\AA^{-1}",
        positive=False,
    )
    unit_oop = pyFAI_units.UnitFiber(
        f"qoop_{frame}_A^-1",
        scale=0.1,
        label=rf"{label} $q_\perp$ ($\AA^{{-1}}$)",
        equation=equation_oop,
        short_name=f"qoop_{frame}",
        unit_symbol=r"\AA^{-1}",
        positive=False,
    )
    return unit_ip, unit_oop


def outOfPlaneIncreasesWithRow(
    detectorCal, alpha, frame=DEFAULT_FRAME, omega=0.0, chi=0.0, phi=0.0, U=None
):
    """Whether the out-of-plane coordinate grows with the detector row index.

    The detector image is displayed with the row index increasing downwards,
    while the Q-plot draws the out-of-plane coordinate upwards. The two only
    agree visually when the out-of-plane coordinate *de*creases with the row
    index, which is the ordinary upward-scattering case.

    An azimuthal reference that points the surface normal the other way -- the
    inverted, downward-scattering geometry used at ID31 -- reverses this, and
    the reciprocal-space image then appears vertically mirrored with respect to
    the detector image. The caller can invert the Q-plot's ordinate to
    compensate.

    Only the sign matters here, so the wave vector and the beam offset drop
    out of the comparison and no wavelength is needed.

    :param detectorCal: :class:`DetectorCalibration.Detector2D_SXRD`.
    :param float alpha: incidence angle in rad.
    :param str frame: one of :data:`FRAMES`.
    :param float omega, chi, phi: sample circle angles in rad.
    :param U: orientation matrix, required for :data:`FRAMES_REQUIRING_U`.
    :rtype: bool
    """
    shape = detectorCal.detector.shape
    rows = np.array([0.0, float(shape[0] - 1)])
    cols = np.full(2, 0.5 * (shape[1] - 1))
    param = np.array(
        [
            detectorCal.dist,
            detectorCal.poni1,
            detectorCal.poni2,
            detectorCal.rot1,
            detectorCal.rot2,
            detectorCal.rot3,
        ],
        dtype=np.float64,
    )
    z, y, x = detectorCal.calc_pos_zyx(d1=rows, d2=cols, param=param, do_parallax=True)

    G, _ = conversionCoefficients(
        frame,
        incident_angle=alpha,
        tilt_angle=tiltAngleFromAzimuth(detectorCal.getAzimuthalReference()),
        omega=omega,
        chi=chi,
        phi=phi,
        U=U,
    )
    # the out-of-plane component up to the positive factor k and a constant
    # offset, both of which cancel in the comparison
    direction = np.stack([x, y, z], axis=-1)
    direction = direction / np.linalg.norm(direction, axis=-1, keepdims=True)
    out_of_plane = direction @ G[2]
    return bool(out_of_plane[1] > out_of_plane[0])


def integrateImage(
    detectorCal,
    image,
    alpha,
    frame=DEFAULT_FRAME,
    omega=0.0,
    chi=0.0,
    phi=0.0,
    U=None,
    npt=1000,
):
    """Rebin ``image`` onto a regular grid of momentum transfer. Experimental.

    :param detectorCal: :class:`DetectorCalibration.Detector2D_SXRD`.
    :param numpy.ndarray image: detector image.
    :param float alpha: incidence angle in rad.
    :param str frame: one of :data:`FRAMES`.
    :param float omega, chi, phi: sample circle angles in rad.
    :param U: orientation matrix, required for :data:`FRAMES_REQUIRING_U`.
    :param int npt: number of bins along each axis.
    :returns: the pyFAI 2D integration result.
    """
    if not HAS_FIBER_INTEGRATOR:
        raise RuntimeError("pyFAI >= 2025.1 with FiberIntegrator is required")

    unit_ip, unit_oop = fiberUnits(frame, omega=omega, chi=chi, phi=phi, U=U)
    tilt = tiltAngleFromAzimuth(detectorCal.getAzimuthalReference())
    # the coefficients identify the conversion completely, see _IntegratorCache
    G, c = conversionCoefficients(
        frame,
        incident_angle=alpha,
        tilt_angle=tilt,
        omega=omega,
        chi=chi,
        phi=phi,
        U=U,
    )
    integrator = _INTEGRATOR.get(detectorCal, G, c)
    return integrator.integrate2d_grazing_incidence(
        image,
        npt_ip=npt,
        npt_oop=npt,
        sample_orientation=SAMPLE_ORIENTATION,
        incident_angle=alpha,
        tilt_angle=tilt,
        unit_ip=unit_ip,
        unit_oop=unit_oop,
    )

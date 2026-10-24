"""Native batches versus the scalar reference, including thread boundaries."""

import threading
import sys

import numpy as np
import pytest
from scipy import stats

from orgui.datautils.xrayutils.corrections import beamprofile as bp
from orgui.datautils.xrayutils.corrections import sample_interception as si

native = pytest.importorskip("orgui.datautils.xrayutils._sample_interception_cpp")


def _inputs(kind="rectangle", profiles="normal", count=17):
    if kind == "polygon":
        shape = si.SampleShape("polygon", vertices=(
            (-.005, -.005), (.005, -.005), (.005, -.002), (-.002, -.002),
            (-.002, .002), (.005, .002), (.005, .005), (-.005, .005),
        ))
    else:
        shape = si.SampleShape(kind, (.01, .006) if kind == "rectangle" else (.01,))
    vertical = bp.gaussian_profile(160e-6, offset=20e-6)
    horizontal = bp.gaussian_profile(.0028, offset=-.0003)
    if profiles == "uniform":
        horizontal = bp.top_hat_profile(.004, offset=.0003)
    elif profiles == "measured":
        z, y = np.linspace(-.0004, .0004, 41), np.linspace(-.006, .006, 51)
        vertical = bp.MeasuredBeamProfile(z, vertical.density_at(z), offset=20e-6)
        horizontal = bp.MeasuredBeamProfile(y, horizontal.density_at(y), offset=-.0003)
    elif profiles in {"trapezoid", "triangle", "smoothed", "gennorm", "skewnorm"}:
        def make(width, offset):
            """Exercise both axes with shifted, non-Gaussian densities."""
            constructors = {
                "trapezoid": lambda: bp.trapezoid_profile(
                    width, width * .35, offset=offset),
                "triangle": lambda: bp.triangular_profile(width, offset=offset),
                "smoothed": lambda: bp.smoothed_top_hat_profile(
                    width, width * .08, offset=offset),
                "gennorm": lambda: bp.generalized_normal_profile(
                    width, 3.5, offset=offset),
                "skewnorm": lambda: bp.skew_normal_profile(
                    width, -2.5, center="median", offset=offset),
            }
            return constructors[profiles]()

        vertical, horizontal = make(160e-6, 20e-6), make(.0028, -.0003)
    return (shape, vertical, horizontal, np.linspace(.0001, .08, count),
            np.linspace(.01, 1.5, count))


@pytest.mark.parametrize("kind", ["rectangle", "circle", "polygon"])
@pytest.mark.parametrize("profiles", ["normal", "uniform", "measured", "trapezoid",
                                     "triangle", "smoothed", "gennorm", "skewnorm"])
def test_native_batches_match_reference_and_threads(kind, profiles, monkeypatch):
    """All shape types preserve shifted profiles, clipping, values and bounds."""
    args = _inputs(kind, profiles)
    kwargs = dict(offset=(.002, -.001), orientation=.13,
                  axis_vertical=10e-6, axis_horizontal=.0002, rtol=1e-10)
    serial = si.overlap(*args, workers=1, **kwargs)
    threaded = si.overlap(*args, workers=4, **kwargs)
    monkeypatch.setattr(si, "_sample_interception_cpp", None)
    reference = si.overlap(*args, **kwargs)
    for field in ("illumination", "fraction", "error", "effective_area"):
        np.testing.assert_array_equal(getattr(serial, field), getattr(threaded, field))
    np.testing.assert_allclose(serial.illumination, reference.illumination,
                               rtol=1e-10, atol=1e-13)
    assert np.all(np.abs(serial.illumination - reference.illumination)
                  <= serial.error + reference.error + 1e-13)


@pytest.mark.parametrize("profiles", ["trapezoid", "triangle", "smoothed",
                                     "gennorm", "skewnorm"])
def test_named_profiles_resolve_without_python(profiles, monkeypatch):
    """Both axes of each named family actually execute in the native kernel."""
    def forbidden(*args):
        """A successful result cannot be a silently scalar reference result."""
        pytest.fail("ordinary built-in profile batch unexpectedly used Python")

    monkeypatch.setattr(si, "_frame_overlap", forbidden)
    si.overlap(*_inputs("circle", profiles, count=5), workers=4)


@pytest.mark.parametrize("profile", [
    bp.DistributionBeamProfile(stats.trapezoid(.15, .8, .003, .004)),
    bp.DistributionBeamProfile(stats.gennorm(1.3, .003, .004)),
    bp.DistributionBeamProfile(stats.skewnorm(-4, .003, .004)),
    bp.DistributionBeamProfile(stats.trapezoid(c=0, d=1, loc=.003, scale=.004)),
])
def test_distribution_location_scale_and_interval_tails(profile, monkeypatch):
    """Positional/keyword frozen parameters preserve metres, shifts and tails."""
    shape = si.SampleShape("rectangle", (.01, .006))
    kwargs = dict(alpha=.04, azimuth=0, axis_horizontal=np.array([0, .01, -.01, 1e-5]))
    result = si.overlap(shape, bp.GaussianBeamProfile(160e-6), profile, **kwargs)
    monkeypatch.setattr(si, "_sample_interception_cpp", None)
    expected = si.overlap(shape, bp.GaussianBeamProfile(160e-6), profile, **kwargs)
    np.testing.assert_allclose(result.illumination, expected.illumination,
                               rtol=1e-10, atol=1e-13)


def test_arbitrary_scipy_distribution_keeps_fallback(monkeypatch):
    """Unsupported continuous families retain their complete SciPy behavior."""
    def forbidden(*args):
        """Custom distributions must never be serialized as a known family."""
        pytest.fail("unsupported distribution reached native code")

    monkeypatch.setattr(native, "overlap_batch", forbidden)
    shape, vertical, _, alpha, azimuth = _inputs(count=2)
    horizontal = bp.DistributionBeamProfile(stats.cauchy(scale=.003), center="median")
    assert si._native_profile(horizontal) is None
    result = si.overlap(shape, vertical, horizontal, alpha, azimuth)
    assert np.all(np.isfinite(result.illumination))


@pytest.mark.parametrize("exponent", [.4, .7, 1.3, 12.])
def test_native_subdivides_panels_before_scipy_fallback(exponent, monkeypatch):
    """Cusps/steep flanks require subdivision and resolve entirely in C++."""
    shape = si.SampleShape("circle", (.01,))
    args = (shape, bp.generalized_normal_profile(160e-6, exponent),
            bp.GaussianBeamProfile(.003), [.001, .04], [.3, .7])

    def forbidden(*args):
        """Prove that native adaptive quadrature resolved both geometries."""
        pytest.fail("native adaptive integration unexpectedly used Python")

    with monkeypatch.context() as context:
        context.setattr(si, "_frame_overlap", forbidden)
        actual = si.overlap(*args, rtol=1e-11)
    subdivisions = []
    original_quad = si.quad

    def adaptive(*args, **kwargs):
        """Count fixed panels that genuinely require adaptive subdivision."""
        subdivisions.append(1)
        return original_quad(*args, **kwargs)

    monkeypatch.setattr(si, "quad", adaptive)
    monkeypatch.setattr(si, "_sample_interception_cpp", None)
    expected = si.overlap(*args, rtol=1e-11)
    assert subdivisions
    np.testing.assert_allclose(actual.illumination, expected.illumination, rtol=1e-11)
    assert np.all(np.abs(actual.illumination - expected.illumination)
                  <= actual.error + expected.error + 1e-13)


@pytest.mark.parametrize("family", ["smoothed", "gennorm", "skewnorm"])
def test_native_profile_mass_retains_remote_tails(family, monkeypatch):
    """Huge integration intervals cannot hide the tails beyond quantile knots."""
    _, _, horizontal, *_ = _inputs(profiles=family, count=1)
    shape = si.SampleShape("rectangle", (.01, 1e7))
    vertical = bp.GaussianBeamProfile(160e-6)

    def forbidden(*args):
        """Remote tails must be resolved in native quadrature as well."""
        pytest.fail("native tail integration unexpectedly used Python")

    monkeypatch.setattr(si, "_frame_overlap", forbidden)
    actual = si.overlap(shape, vertical, horizontal, .04, 0.)
    half_span = .005 * np.sin(.04)
    expected = vertical.interval_mass(-half_span, half_span) / np.sin(.04)
    np.testing.assert_allclose(actual.illumination, expected, rtol=1e-12)


def test_native_array_shapes_deduplication_and_progress(monkeypatch):
    """Strides/broadcasting and duplicates survive bounded native dispatch."""
    args = _inputs(count=514)
    calls, updates = [], []
    caller = threading.get_ident()
    original = native.overlap_batch

    def batch(*args):
        """Count unique geometries passed to the compiled implementation."""
        calls.append(len(args[3]))
        return original(*args)

    def progress(completed, total):
        """Callbacks stay on the initiating Python thread."""
        assert threading.get_ident() == caller
        updates.append((completed, total))

    monkeypatch.setattr(native, "overlap_batch", batch)
    alpha = np.tile(args[3], (2, 1))[:, ::-1]
    azimuth = args[4][::-1]
    alpha.flags.writeable = False
    result = si.overlap(*args[:3], alpha, azimuth, progress=progress, workers=4)
    assert calls == [256, 256, 2]
    assert updates[0] == (0, 1028)
    assert updates[-1] == (1028, 1028)
    assert result.illumination.shape == (2, 514)
    np.testing.assert_array_equal(result.illumination[0], result.illumination[1])
    assert si.overlap(*args[:3], .01, .3).illumination.shape == ()
    assert si.overlap(*args[:3], [], []).illumination.shape == (0,)


def test_native_cancellation_stops_before_next_batch(monkeypatch):
    """Cancellation at a batch boundary cannot dispatch remaining frames."""
    calls = []
    original = native.overlap_batch

    def batch(*args):
        """Record actual native calls before cancellation is observed."""
        calls.append(len(args[3]))
        return original(*args)

    monkeypatch.setattr(native, "overlap_batch", batch)
    with pytest.raises(InterruptedError, match="cancelled"):
        si.overlap(*_inputs(count=514), progress=lambda done, total: done < 256)
    assert calls == [256]


def test_unresolved_native_panel_uses_reference(monkeypatch):
    """Only unresolved geometries take the scalar adaptive fallback."""
    original = native.overlap_batch
    scalar = si._frame_overlap
    retried = []

    def unresolved(*args):
        """Force one result to require adaptive calculation."""
        values, errors, resolved = original(*args)
        resolved[2] = False
        return values, errors, resolved

    def retry(*args):
        """Record each adaptive fallback geometry."""
        retried.append(args[3:7])
        return scalar(*args)

    args = _inputs(count=4)
    monkeypatch.setattr(native, "overlap_batch", unresolved)
    monkeypatch.setattr(si, "_frame_overlap", retry)
    result = si.overlap(*args)
    assert len(retried) == 1
    assert retried[0][0] == args[3][2]
    monkeypatch.setattr(si, "_sample_interception_cpp", None)
    expected = si.overlap(*args)
    np.testing.assert_allclose(result.illumination, expected.illumination, rtol=1e-11)


def test_custom_profile_never_enters_native(monkeypatch):
    """Arbitrary Python overrides are never approximated or called without GIL."""
    class Custom(bp.GaussianBeamProfile):
        def interval_mass(self, lo, hi):
            """Keep the documented scalar-only extension contract."""
            assert np.ndim(lo) == np.ndim(hi) == 0
            return super().interval_mass(lo, hi)

    def forbidden(*args):
        """Fail if a custom profile reaches native dispatch."""
        pytest.fail("custom profile must use Python")

    monkeypatch.setattr(native, "overlap_batch", forbidden)
    args = _inputs(count=2)
    si.overlap(args[0], args[1], Custom(.003), *args[3:])


def test_native_validation_rejects_malformed_inputs():
    """Invalid direct calls raise before worker threads or pointer access."""
    shape, vertical, horizontal, *_ = _inputs()
    descriptors = [dict(kind=shape.kind, dimensions=shape.dimensions,
                        vertices=shape.outline()),
                   si._native_profile(vertical), si._native_profile(horizontal)]
    for geometry in (np.ones(4), np.ones((3, 3)), np.full((1, 4), np.nan)):
        with pytest.raises(ValueError):
            native.overlap_batch(*descriptors, geometry, (0., 0.), 0., 1e-9, 4)
    with pytest.raises(ValueError, match="workers"):
        si.overlap(*_inputs(), workers=-1)


def test_native_call_releases_gil():
    """Another Python thread executes during a serial native batch."""
    shape, vertical, horizontal, alpha, azimuth = _inputs(count=1024)
    descriptors = [dict(kind=shape.kind, dimensions=shape.dimensions,
                        vertices=shape.outline()),
                   si._native_profile(vertical), si._native_profile(horizontal)]
    geometry = np.column_stack((alpha, azimuth, np.zeros((len(alpha), 2))))
    ready, start, observed = threading.Event(), threading.Event(), threading.Event()

    def observer():
        """Wait until the caller enters the region without Python scheduling."""
        ready.set()
        start.wait()
        observed.set()

    thread = threading.Thread(target=observer)
    thread.start()
    ready.wait()
    interval = sys.getswitchinterval()
    try:
        # Prevent an ordinary interpreter timeslice from satisfying the test.
        sys.setswitchinterval(10)
        start.set()
        native.overlap_batch(*descriptors, geometry, (0., 0.), 0., 1e-9, 1)
        released = observed.is_set()
    finally:
        sys.setswitchinterval(interval)
        thread.join()
    assert released

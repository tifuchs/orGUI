// Copyright (c) 2026 Timo Fuchs
// SPDX-License-Identifier: MIT
// Exact surface-chord geometry: metres, radians, normalized densities [1/m].
// No Python objects/callbacks are accessed while native frame workers run.
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace {
constexpr double pi = 3.1415926535897932384626433832795;
constexpr double eps = std::numeric_limits<double>::epsilon();
using Point = std::array<double, 2>;
using Frame = std::array<double, 4>;
using FloatArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

void finite(const std::vector<double>& values) {
    for (double value : values)
        if (!std::isfinite(value)) throw std::invalid_argument("finite values required");
}

double normal_cdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// A numerical failure is local to a geometry; Python/SciPy retries that row.
struct Unresolved {};

struct Profile {
    std::string kind;
    double shift, peak, scale = 1, left = 0, right = 1, parameter = 0, log_norm = 0;
    std::vector<double> points, positions, density, cumulative;

    explicit Profile(const py::dict& data)
        : kind(py::cast<std::string>(data["kind"])),
          shift(py::cast<double>(data["shift"])), peak(py::cast<double>(data["peak"])),
          points(py::cast<std::vector<double>>(data["points"])) {
        finite(points);
        finite({shift, peak});
        if (!(peak > 0)) throw std::invalid_argument("positive profile peak required");
        if (kind == "normal" || kind == "smoothed" || kind == "trapezoid"
            || kind == "gennorm" || kind == "skewnorm") {
            scale = py::cast<double>(data["scale"]);
            if (!(scale > 0) || !std::isfinite(scale))
                throw std::invalid_argument("positive finite profile scale required");
            if (kind == "smoothed") {
                parameter = py::cast<double>(data["width"]);
                if (!(parameter > 0) || !std::isfinite(parameter))
                    throw std::invalid_argument("positive finite smoothed width required");
                // Explicit edge transitions supplement quantiles when the
                // Gaussian blur is much narrower than the slit aperture.
                for (double edge : {-parameter / 2, parameter / 2})
                    for (double k : {-8., -4., -1., 0., 1., 4., 8.})
                        points.push_back(edge + k * scale - shift);
            } else if (kind != "normal") {
                auto params = py::cast<std::vector<double>>(data["parameters"]);
                finite(params);
                if (params.size() != (kind == "trapezoid" ? 2u : 1u))
                    throw std::invalid_argument("invalid distribution parameters");
                if (kind == "trapezoid") {
                    left = params[0]; right = params[1];
                    if (!(0 <= left && left <= right && right <= 1))
                        throw std::invalid_argument("invalid trapezoid shoulders");
                    for (double p : {0., left, right, 1.}) points.push_back(p * scale - shift);
                } else {
                    parameter = params[0];
                    if (kind == "gennorm") {
                        if (!(parameter > 0)) throw std::invalid_argument("positive exponent required");
                        log_norm = std::log(parameter / 2) - std::lgamma(1 / parameter);
                    } else if (parameter != 0) {
                        for (double k : {-8., -1., 1., 8.})
                            points.push_back(k * scale / std::abs(parameter) - shift);
                    }
                    points.push_back(-shift);
                }
            }
        } else if (kind == "uniform") {
            auto bounds = py::cast<Point>(data["bounds"]);
            left = bounds[0]; right = bounds[1];
            finite({left, right});
            if (!(right > left)) throw std::invalid_argument("ordered beam bounds required");
        } else if (kind == "measured") {
            positions = py::cast<std::vector<double>>(data["positions"]);
            density = py::cast<std::vector<double>>(data["density"]);
            cumulative = py::cast<std::vector<double>>(data["cumulative"]);
            finite(positions); finite(density); finite(cumulative);
            if (positions.size() < 2 || density.size() != positions.size()
                || cumulative.size() != positions.size())
                throw std::invalid_argument("matching measured profile arrays required");
            for (std::size_t i = 0; i < positions.size(); ++i) {
                if (density[i] < 0 || (i && positions[i] <= positions[i - 1]))
                    throw std::invalid_argument("ordered nonnegative measured profile required");
            }
        } else throw std::invalid_argument("unsupported native beam profile");
        std::sort(points.begin(), points.end());
        points.erase(std::unique(points.begin(), points.end()), points.end());
    }

    std::size_t segment(double x) const {
        auto it = std::upper_bound(positions.begin(), positions.end(), x);
        return std::min(positions.size() - 2,
                        static_cast<std::size_t>(it - positions.begin() - 1));
    }

    double pdf(double x) const {
        x += shift;
        if (kind == "normal") {
            double u = x / scale;
            return std::exp(-0.5 * u * u) / (std::sqrt(2 * pi) * scale);
        }
        if (kind == "uniform") return x >= left && x <= right ? 1 / (right - left) : 0;
        if (kind == "trapezoid") return trapezoid_pdf(x / scale) / scale;
        if (kind == "gennorm" || kind == "skewnorm" || kind == "smoothed")
            return unit_pdf(x / scale) / scale;
        if (x < positions.front() || x > positions.back()) return 0;
        auto i = segment(x);
        double t = (x - positions[i]) / (positions[i + 1] - positions[i]);
        return density[i] + t * (density[i + 1] - density[i]);
    }

    double unit_pdf(double u) const {
        if (kind == "gennorm")
            return std::exp(log_norm - std::pow(std::abs(u), parameter));
        if (kind == "skewnorm") {
            return 2 * std::exp(-u * u / 2) * normal_cdf(parameter * u)
                / std::sqrt(2 * pi);
        }
        double width = parameter / scale;
        double lo = u - width / 2, hi = u + width / 2;
        if (width * (1 + std::abs(u)) < 1e-7)
            return std::exp(-u * u / 2) / std::sqrt(2 * pi);
        return (lo > 0 ? normal_cdf(-lo) - normal_cdf(-hi)
                       : normal_cdf(hi) - normal_cdf(lo)) / width;
    }

    double measured_cdf(double x) const {
        x = std::clamp(x + shift, positions.front(), positions.back());
        auto i = segment(x);
        double d = positions[i + 1] - positions[i], t = (x - positions[i]) / d;
        return cumulative[i] + d * t * (density[i] + 0.5 * t * (density[i + 1] - density[i]));
    }

    double trapezoid_pdf(double x) const {
        if (x < 0 || x > 1) return 0;
        double height = 2 / (1 + right - left);
        if (x < left) return height * x / left;
        if (x > right) return height * (1 - x) / (1 - right);
        return height;
    }

    double numerical_mass(double lo, double hi, double* error) const;

    bool smoothed_mass(double lo, double hi, double& value, double& error) const {
        lo = (lo + shift) / scale; hi = (hi + shift) / scale;
        double width = parameter / scale;
        // Integrate Phi analytically. Reflect positive tails so a small
        // interval probability never comes from subtracting two near-ones.
        auto cdf = [&](double x, double& uncertainty) {
            auto primitive = [](double u, double& magnitude) {
                double a = u * normal_cdf(u), b = std::exp(-u * u / 2) / std::sqrt(2 * pi);
                magnitude = std::abs(a) + b;
                return a + b;
            };
            double a = 0, b = 0;
            double result = (primitive(x + width / 2, a) - primitive(x - width / 2, b)) / width;
            uncertainty = 16 * eps * (a + b) / width;
            return result;
        };
        if (lo >= 0) { double previous = lo; lo = -hi; hi = -previous; }
        double a = 0, b = 0;
        if (hi <= 0) {
            value = cdf(hi, a) - cdf(lo, b);
            error = a + b + 2 * eps * std::abs(value);
        } else {
            value = 1 - cdf(lo, a) - cdf(-hi, b);
            error = a + b + 4 * eps;
        }
        return std::isfinite(value) && std::isfinite(error) && value >= 0 && value <= 1
            && error <= std::max(2e-15, 5e-14 * value);
    }

    bool has_numerical_mass() const {
        return kind == "gennorm" || kind == "skewnorm" || kind == "smoothed";
    }

    double mass(double lo, double hi, double* error = nullptr) const {
        if (error) *error = 0;
        if (hi < lo) throw std::invalid_argument("reversed profile interval");
        if (lo == hi) return 0;
        if (kind == "measured") return measured_cdf(hi) - measured_cdf(lo);
        if (kind == "smoothed") {
            double value = 0, uncertainty = 0;
            if (smoothed_mass(lo, hi, value, uncertainty)) {
                if (error) *error = uncertainty;
                return value;
            }
        }
        if (has_numerical_mass())
            return numerical_mass(lo, hi, error);
        lo += shift; hi += shift;
        if (kind == "uniform")
            return std::max(0.0, std::min(hi, right) - std::max(lo, left)) / (right - left);
        lo /= scale; hi /= scale;
        if (kind == "trapezoid") {
            // Integrate each linear segment locally: no CDF cancellation in
            // tiny intervals or near the two support endpoints.
            double total = 0;
            std::array<double, 4> knots{0, left, right, 1};
            for (std::size_t i = 1; i < knots.size(); ++i) {
                double a = std::max(lo, knots[i - 1]), b = std::min(hi, knots[i]);
                if (b > a) total += (b - a) * (trapezoid_pdf(a) + trapezoid_pdf(b)) / 2;
            }
            return total;
        }
        double midpoint = (hi + lo) / 2;
        if ((hi - lo) * (1 + std::abs(midpoint)) < 1e-7)
            return (hi - lo) * std::exp(-midpoint * midpoint / 2) / std::sqrt(2 * pi);
        return lo > 0 ? normal_cdf(-lo) - normal_cdf(-hi) : normal_cdf(hi) - normal_cdf(lo);
    }
};

Point rotate(Point p, double angle) {
    double c = std::cos(angle), s = std::sin(angle);
    return {c * p[0] - s * p[1], s * p[0] + c * p[1]};
}

struct Shape {
    std::string kind;
    std::vector<double> dimensions;
    std::vector<Point> vertices;
    explicit Shape(const py::dict& data)
        : kind(py::cast<std::string>(data["kind"])),
          dimensions(py::cast<std::vector<double>>(data["dimensions"])) {
        finite(dimensions);
        if (kind != "rectangle" && kind != "circle" && kind != "polygon")
            throw std::invalid_argument("unsupported sample shape");
        if (kind == "circle" || kind == "rectangle") {
            if (dimensions.size() != (kind == "circle" ? 1u : 2u))
                throw std::invalid_argument("invalid sample dimensions");
            for (double d : dimensions)
                if (!(d > 0)) throw std::invalid_argument("positive sample dimensions required");
        }
        if (kind != "circle") {
            vertices = py::cast<std::vector<Point>>(data["vertices"]);
            if (vertices.size() < 3) throw std::invalid_argument("three vertices required");
            for (auto v : vertices) finite({v[0], v[1]});
        }
    }
};

struct Result { double value = 0, error = 0; bool resolved = false; };

// Same GK21 nodes, weights and conservative error estimator as _vector_quad.
constexpr double nodes[] = {
    .9956571630258080807, .9739065285171717201, .9301574913557082260,
    .8650633666889845107, .7808177265864168971, .6794095682990244062,
    .5627571346686046833, .4333953941292471908, .2943928627014601981,
    .1488743389816312109, 0
};
constexpr double weights[] = {
    .0116946388673718743, .0325581623079647275, .0547558965743519960,
    .0750396748109199528, .0931254545836976055, .1093871588022976419,
    .1234919762620658511, .1347092173114733259, .1427759385770600808,
    .1477391049013384914, .1494455540029169057
};
constexpr double gauss_weights[] = {
    .0666713443086881376, .1494513491505805931, .2190863625159820440,
    .2692667193099963551, .2955242247147528702
};

template<class Function> Result panel(const Function& integrand, double rtol, double atol) {
    std::array<double, 21> values;
    double total = 0, gauss = 0, deviation = 0, absolute = 0;
    for (int i = 0; i < 21; ++i) {
        int k = i <= 10 ? i : 20 - i;
        double node = (1 + (i <= 10 ? nodes[k] : -nodes[k])) / 2;
        values[i] = integrand(node);
        total += weights[k] * values[i] / 2;
        if (i % 2) gauss += gauss_weights[k / 2] * values[i] / 2;
    }
    for (int i = 0; i < 21; ++i) {
        int k = i <= 10 ? i : 20 - i;
        deviation += weights[k] * std::abs(values[i] - total) / 2;
        absolute += weights[k] * std::abs(values[i]) / 2;
    }
    double error = std::abs(total - gauss);
    if (deviation != 0 && error != 0)
        error = deviation * std::min(1.0, std::pow(200 * error / deviation, 1.5));
    error = std::max(error, 50 * eps * absolute);
    return {total, error, std::isfinite(total) && std::isfinite(error)
        && error <= std::max(atol, rtol * std::abs(total))
        && (error != deviation || error == 0)};
}

template<class Function>
Result adaptive(const Function& integrand, double rtol, double atol = 1e-12) {
    Result first = panel(integrand, rtol, atol);
    if (first.resolved || !std::isfinite(first.value) || !std::isfinite(first.error)) return first;
    struct Interval { double lo, hi; Result result; };
    std::vector<Interval> intervals{{0, 1, first}};
    double total = first.value, error = first.error;
    // Match the reference's 200-panel work bound. No extrapolation or
    // relaxed tolerance: roundoff/stagnation or exhaustion retries in SciPy.
    int stagnant = 0;
    while (intervals.size() < 200) {
        auto worst = std::max_element(intervals.begin(), intervals.end(),
            [](const Interval& a, const Interval& b) { return a.result.error < b.result.error; });
        Interval old = *worst;
        double mid = (old.lo + old.hi) / 2;
        if (mid == old.lo || mid == old.hi) return {};
        auto evaluate = [&](double lo, double hi) {
            return panel([&](double t) { return integrand(lo + (hi - lo) * t) * (hi - lo); },
                         rtol, atol);
        };
        Result a = evaluate(old.lo, mid), b = evaluate(mid, old.hi);
        if (!std::isfinite(a.value + b.value) || !std::isfinite(a.error + b.error)) return {};
        total += (a.value + b.value) - old.result.value;
        error = std::max(0.0, error + (a.error + b.error) - old.result.error);
        if (error <= std::max(atol, rtol * std::abs(total))) return {total, error, true};
        if (a.error + b.error >= .99 * old.result.error
            && std::abs(a.value + b.value - old.result.value) <= 1e-5 * std::abs(a.value + b.value))
            if (++stagnant >= 10) return {};
        *worst = {old.lo, mid, a};
        intervals.push_back({mid, old.hi, b});
    }
    return {};
}

double Profile::numerical_mass(double lo, double hi, double* error) const {
    // Integrate the exact density in dimensionless distribution coordinates,
    // splitting at quantiles and sharp transitions. Neither tails nor narrow
    // peaks are replaced by a lookup table. Inner errors enter the H bound.
    double a = (lo + shift) / scale, end = (hi + shift) / scale;
    if (!std::isfinite(a) || !std::isfinite(end)) throw Unresolved{};
    double total = 0, uncertainty = 0;
    auto integrate = [&](double b) {
        std::vector<Point> pending{{a, b}};
        int splits = 0;
        while (!pending.empty()) {
            auto bounds = pending.back(); pending.pop_back();
            double start = bounds[0], width = bounds[1] - start;
            // A very distant sample edge must not hide the residual tail
            // next to the last quantile from every Gauss-Kronrod node.
            if (width > 4 * std::max(1.0, std::min(std::abs(start), std::abs(bounds[1])))) {
                double middle = start + width / 2;
                if (++splits > 2048 || middle == start || middle == bounds[1]) throw Unresolved{};
                pending.push_back({middle, bounds[1]});
                pending.push_back({start, middle});
                continue;
            }
            Result part = adaptive([&](double t) { return unit_pdf(start + width * t) * width; },
                                   5e-14, 2e-15 / (points.size() + 1));
            if (!part.resolved) throw Unresolved{};
            total += part.value; uncertainty += part.error;
        }
        a = b;
    };
    for (double point : points) {
        double b = (point + shift) / scale;
        if (b > a && b < end) integrate(b);
    }
    if (end > a) integrate(end);
    if (!std::isfinite(total) || !std::isfinite(uncertainty)
        || uncertainty > std::max(2e-15, 5e-14 * std::abs(total))) throw Unresolved{};
    if (error) *error = uncertainty;
    return total;
}

Result frame_overlap(const Shape& shape, const Profile& vertical, const Profile& horizontal,
                     Frame f, Point offset, double orientation, double rtol) {
    double alpha = f[0], azimuth = f[1], zv = f[2], yh = f[3];
    double sine = std::sin(alpha);
    Point center = rotate(offset, azimuth);
    if (shape.kind == "rectangle" && std::fmod(azimuth + orientation, pi) == 0
        && shape.dimensions[0] * sine * vertical.peak >= 1e-7) {
        double zcenter = zv + center[0] * sine, span = shape.dimensions[0] * sine;
        double ve = 0, he = 0;
        double v = vertical.mass(zcenter - span / 2, zcenter + span / 2, &ve);
        double h = horizontal.mass(yh + center[1] - shape.dimensions[1] / 2,
                                   yh + center[1] + shape.dimensions[1] / 2, &he);
        double value = v * h / sine;
        double error = (ve * std::abs(h) + he * std::abs(v) + ve * he) / sine
            + std::abs(value) * eps * 8;
        return {value, error, std::isfinite(value) && std::isfinite(error)
            && error <= std::max(1e-11 * shape.dimensions[0] * vertical.peak, std::abs(value) * rtol * 5)};
    }
    bool circle = shape.kind == "circle";
    double radius = circle ? shape.dimensions[0] / 2 : 0;
    double lo = center[0] - radius, hi = center[0] + radius;
    std::vector<Point> vertices, edges;
    if (!circle) {
        lo = std::numeric_limits<double>::infinity(); hi = -lo;
        for (auto v : shape.vertices) {
            v = rotate(v, orientation);
            v = rotate({v[0] + offset[0], v[1] + offset[1]}, azimuth);
            vertices.push_back(v); lo = std::min(lo, v[0]); hi = std::max(hi, v[0]);
        }
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            auto a = vertices[i], b = vertices[(i + 1) % vertices.size()];
            edges.push_back({b[0] - a[0], b[1] - a[1]});
        }
    }
    double span = hi - lo;
    if (!(span > 0) || !std::isfinite(span)) return {};
    std::vector<double> points;
    for (double p : vertical.points) points.push_back((p - zv) / sine);
    if (!circle) {
        for (auto v : vertices) points.push_back(v[0]);
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            auto a = vertices[i], edge = edges[i];
            if (edge[1] != 0) for (double p : horizontal.points) {
                double t = (p - yh - a[1]) / edge[1];
                if (t > 0 && t < 1) points.push_back(a[0] + t * edge[0]);
            }
        }
        points.erase(std::remove_if(points.begin(), points.end(),
                     [=](double p) { return !(p > lo && p < hi); }), points.end());
        points.push_back(lo); points.push_back(hi);
    } else {
        std::vector<double> transformed{0, 1};
        for (double p : points) if (std::abs(p - center[0]) < radius)
            transformed.push_back(.5 + std::asin((p - center[0]) / radius) / pi);
        for (double p : horizontal.points) {
            double height = p - yh - center[1];
            if (std::abs(height) < radius) {
                double theta = std::acos(std::abs(height) / radius);
                transformed.push_back(.5 - theta / pi); transformed.push_back(.5 + theta / pi);
            }
        }
        points = std::move(transformed);
    }
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    double total = 0, error = 0;
    for (std::size_t i = 1; i < points.size(); ++i) {
        double a = points[i - 1], b = points[i], weight = (b - a) / (circle ? 1 : span);
        std::vector<Point> hits;
        if (!circle) {
            for (std::size_t j = 0; j < vertices.size(); ++j) {
                auto v = vertices[j], e = edges[j], end = vertices[(j + 1) % vertices.size()];
                if (e[0] != 0 && std::min(v[0], end[0]) <= a && b <= std::max(v[0], end[0]))
                    hits.push_back({v[1] + ((a - v[0]) / e[0]) * e[1], ((b - a) / e[0]) * e[1]});
            }
            std::stable_sort(hits.begin(), hits.end(), [](Point x, Point y) {
                return x[0] + x[1] / 2 < y[0] + y[1] / 2;
            });
            if (hits.size() % 2) return {};
        }
        double mass_error = 0;
        auto integrand = [&](double t) {
            double x, mass = 0, jacobian = 1, uncertainty = 0;
            if (circle) {
                double theta = pi * (a + (b - a) * t - .5), cosine = std::cos(theta);
                x = center[0] + radius * std::sin(theta);
                mass = horizontal.mass(yh + center[1] - radius * cosine,
                                       yh + center[1] + radius * cosine, &uncertainty);
                jacobian = pi / 2 * cosine;
            } else {
                x = a + (b - a) * t;
                for (std::size_t j = 0; j < hits.size(); j += 2) {
                    double local_error = 0;
                    mass += horizontal.mass(yh + hits[j][0] + hits[j][1] * t,
                                            yh + hits[j + 1][0] + hits[j + 1][1] * t, &local_error);
                    uncertainty += local_error;
                }
            }
            mass_error = std::max(mass_error, uncertainty);
            return vertical.pdf(zv + x * sine) / vertical.peak * mass * jacobian;
        };
        if (weight <= 8 * eps) {
            total += integrand(.5) * weight;
            error += weight * (circle ? pi / 2 : 1);
        } else {
            Result value = adaptive(integrand, rtol);
            if (!value.resolved) return {};
            total += value.value * weight;
            error += (value.error + mass_error * (circle ? pi / 2 : 1)) * weight;
        }
    }
    double scale = span * vertical.peak;
    return {total * scale, error * scale, std::isfinite(total) && std::isfinite(error)
        && error <= std::max(1e-11, std::abs(total) * rtol * 5)};
}

py::tuple overlap_batch(const py::dict& shape_data, const py::dict& vertical_data,
                       const py::dict& horizontal_data, FloatArray geometry,
                       Point offset, double orientation, double rtol, int workers) {
    // Snapshot all input arrays before releasing the GIL; external NumPy writes
    // cannot race native readers. The output arrays remain private until return.
    Shape shape(shape_data);
    Profile vertical(vertical_data), horizontal(horizontal_data);
    if (geometry.ndim() != 2 || geometry.shape(1) != 4)
        throw std::invalid_argument("geometry must have shape (frames, 4)");
    if (workers < 0 || !(rtol >= 1e-12 && rtol <= 1e-3))
        throw std::invalid_argument("invalid worker count or tolerance");
    finite({offset[0], offset[1], orientation});
    auto n = geometry.shape(0);
    std::vector<Frame> frames(static_cast<std::size_t>(n));
    auto input = geometry.unchecked<2>();
    for (py::ssize_t i = 0; i < n; ++i) {
        for (int j = 0; j < 4; ++j) frames[i][j] = input(i, j);
        finite({frames[i][0], frames[i][1], frames[i][2], frames[i][3]});
        if (!(frames[i][0] > 0 && frames[i][0] <= pi / 2))
            throw std::invalid_argument("incidence must be in (0, pi/2]");
    }
    py::array_t<double> values(n), errors(n);
    py::array_t<bool> resolved(n);
    double* out = values.mutable_data(); double* err = errors.mutable_data();
    bool* ok = resolved.mutable_data();
    unsigned hardware = std::max(1u, std::thread::hardware_concurrency());
    unsigned automatic = static_cast<unsigned>(std::min<py::ssize_t>(8, std::max<py::ssize_t>(1, n / 64)));
    // Two analytical interval masses cost less than launching threads;
    // numerically integrated profile masses can still benefit from workers.
    if (workers == 0 && shape.kind == "rectangle"
        && !vertical.has_numerical_mass() && !horizontal.has_numerical_mass()
        && std::all_of(frames.begin(), frames.end(), [&](Frame f) {
        return std::fmod(f[1] + orientation, pi) == 0
            && shape.dimensions[0] * std::sin(f[0]) * vertical.peak >= 1e-7;
    })) automatic = 1;
    unsigned limit = std::min(hardware, workers == 0 ? automatic : unsigned(workers));
    unsigned count = static_cast<unsigned>(std::min<py::ssize_t>(limit, n));
    {
        py::gil_scoped_release release;
        std::atomic<std::size_t> next{0};
        std::atomic<bool> failed{false};
        std::exception_ptr failure;
        std::mutex failure_mutex;
        auto work = [&]() {
            try {
                while (!failed.load()) {
                    auto i = next.fetch_add(1);
                    if (i >= frames.size()) break;
                    Result result;
                    try {
                        result = frame_overlap(shape, vertical, horizontal, frames[i], offset, orientation, rtol);
                    } catch (const Unresolved&) {
                        result = {};
                    } catch (const std::invalid_argument&) {
                        // Almost coincident polygon crossings can round in a
                        // different order. Let the reference decide that row.
                        result = {};
                    }
                    out[i] = result.value; err[i] = result.error; ok[i] = result.resolved;
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(failure_mutex);
                if (!failure) failure = std::current_exception();
                failed.store(true);
            }
        };
        std::vector<std::thread> threads;
        try {
            for (unsigned i = 1; i < count; ++i) threads.emplace_back(work);
            work();
        } catch (...) {
            failed.store(true);
            for (auto& thread : threads) thread.join();
            throw;
        }
        for (auto& thread : threads) thread.join();
        if (failure) std::rethrow_exception(failure);
    }
    return py::make_tuple(values, errors, resolved);
}
}  // namespace

PYBIND11_MODULE(_sample_interception_cpp, m) {
    m.doc() = "Batched exact surface-chord quadrature with GIL-free native workers.";
    m.def("overlap_batch", &overlap_batch,
          "Evaluate (alpha, azimuth, vertical axis, horizontal axis) rows [rad, rad, m, m].\n"
          "Returns H, its estimated error and a resolved mask; unresolved rows require adaptive fallback.",
          py::arg("shape"), py::arg("vertical"), py::arg("horizontal"), py::arg("geometry"),
          py::arg("offset"), py::arg("orientation"), py::arg("rtol"), py::arg("workers") = 0);
}

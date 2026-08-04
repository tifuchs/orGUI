// Reciprocal-space conversion kernel for the experimental Q-plot.
//
// The whole conversion of a detector image collapses into one affine relation,
// see orgui/app/qconversion.py. With inv = 1 / |(x, y, z)| every Cartesian
// component of the momentum transfer is
//
//     q_j = k * ( (G[j] . (x, y, z)) * inv - c[j] )
//
// so the in-plane and the out-of-plane component can be produced in a single
// pass over the pixel positions. That matters because the arrays are large:
// a Pilatus 6M has more than six million pixels.

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using DoubleArray = py::array_t<double, py::array::c_style | py::array::forcecast>;

// Returns the signed in-plane and the out-of-plane momentum transfer.
// The sign of the in-plane component follows orGUI's delta direction.
py::tuple q_ip_oop(
    const DoubleArray &x,
    const DoubleArray &y,
    const DoubleArray &z,
    const DoubleArray &conversion_matrix,
    const DoubleArray &offset,
    const double k)
{
    const py::buffer_info x_info = x.request();
    const py::buffer_info y_info = y.request();
    const py::buffer_info z_info = z.request();
    if (y_info.size != x_info.size || z_info.size != x_info.size) {
        throw std::invalid_argument(
            "x, y and z must have the same number of elements");
    }

    const py::buffer_info matrix_info = conversion_matrix.request();
    if (matrix_info.ndim != 2 || matrix_info.shape[0] != 3
        || matrix_info.shape[1] != 3) {
        throw std::invalid_argument("the conversion matrix must be 3x3");
    }
    const py::buffer_info offset_info = offset.request();
    if (offset_info.size != 3) {
        throw std::invalid_argument("the offset must have three elements");
    }

    const std::vector<py::ssize_t> shape(
        x_info.shape.begin(), x_info.shape.end());
    py::array_t<double> q_ip(shape);
    py::array_t<double> q_oop(shape);

    const double *const xp = static_cast<const double *>(x_info.ptr);
    const double *const yp = static_cast<const double *>(y_info.ptr);
    const double *const zp = static_cast<const double *>(z_info.ptr);
    const double *const gp = static_cast<const double *>(matrix_info.ptr);
    const double *const cp = static_cast<const double *>(offset_info.ptr);
    double *const ip = static_cast<double *>(q_ip.request().ptr);
    double *const oop = static_cast<double *>(q_oop.request().ptr);

    // Local copies keep the coefficients in registers and let the compiler
    // vectorise the loop.
    const double g00 = gp[0], g01 = gp[1], g02 = gp[2];
    const double g10 = gp[3], g11 = gp[4], g12 = gp[5];
    const double g20 = gp[6], g21 = gp[7], g22 = gp[8];
    const double c0 = cp[0], c1 = cp[1], c2 = cp[2];
    const py::ssize_t count = x_info.size;

    {
        py::gil_scoped_release release;
        for (py::ssize_t i = 0; i < count; ++i) {
            const double xi = xp[i];
            const double yi = yp[i];
            const double zi = zp[i];
            const double inv = 1.0 / std::sqrt(xi * xi + yi * yi + zi * zi);
            const double qx = k * ((g00 * xi + g01 * yi + g02 * zi) * inv - c0);
            const double qy = k * ((g10 * xi + g11 * yi + g12 * zi) * inv - c1);
            const double qz = k * ((g20 * xi + g21 * yi + g22 * zi) * inv - c2);
            const double radial = std::sqrt(qx * qx + qy * qy);
            ip[i] = qx >= 0.0 ? radial : -radial;
            oop[i] = qz;
        }
    }

    return py::make_tuple(q_ip, q_oop);
}

PYBIND11_MODULE(_qconversion_cpp, m)
{
    m.doc() = "Reciprocal-space conversion kernel for the experimental Q-plot.";
    m.def(
        "q_ip_oop",
        &q_ip_oop,
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        py::arg("conversion_matrix"),
        py::arg("offset"),
        py::arg("k"),
        "Return the signed in-plane and the out-of-plane momentum transfer.");
}

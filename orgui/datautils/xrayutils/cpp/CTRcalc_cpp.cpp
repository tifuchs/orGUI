#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using complex128 = std::complex<double>;

using Array1D = py::array_t<double, py::array::c_style>;
using Array2D = py::array_t<double, py::array::c_style>;
using Array3D = py::array_t<double, py::array::c_style>;

struct Matrix3 {
    double value[3][3];
};

struct Inputs {
    py::buffer_info h;
    py::buffer_info k;
    py::buffer_info l;
    py::buffer_info basis;
    py::buffer_info f_factors;
    py::buffer_info ref_hkl_transform;
    py::buffer_info b_mat;
    py::buffer_info r_mat;
    py::buffer_info r_mat_inv;
    py::buffer_info coherent_domain_matrix;
    py::buffer_info coherent_domain_occupancy;
};

inline const double *ptr(const py::buffer_info &info) {
    return static_cast<const double *>(info.ptr);
}

void require_shape(
    const py::buffer_info &info,
    const int ndim,
    const char *name
) {
    if (info.ndim != ndim) {
        throw py::value_error(
            std::string(name) + " has unexpected number of dimensions"
        );
    }
}

Inputs request_inputs(
    const Array1D &h,
    const Array1D &k,
    const Array1D &l,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &ref_hkl_transform,
    const Array2D &b_mat,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy
) {
    Inputs inputs{
        h.request(),
        k.request(),
        l.request(),
        basis.request(),
        f_factors.request(),
        ref_hkl_transform.request(),
        b_mat.request(),
        r_mat.request(),
        r_mat_inv.request(),
        coherent_domain_matrix.request(),
        coherent_domain_occupancy.request(),
    };

    require_shape(inputs.h, 1, "h");
    require_shape(inputs.k, 1, "k");
    require_shape(inputs.l, 1, "l");
    require_shape(inputs.basis, 2, "basis");
    require_shape(inputs.f_factors, 2, "f_factors");
    require_shape(inputs.ref_hkl_transform, 2, "refHKLTransform");
    require_shape(inputs.b_mat, 2, "B_mat");
    require_shape(inputs.r_mat, 2, "R_mat");
    require_shape(inputs.r_mat_inv, 2, "R_mat_inv");
    require_shape(inputs.coherent_domain_matrix, 3, "coherentDomainMatrix");
    require_shape(
        inputs.coherent_domain_occupancy,
        1,
        "coherentDomainOccupancy"
    );

    const auto n = inputs.h.shape[0];
    if (inputs.k.shape[0] != n || inputs.l.shape[0] != n) {
        throw py::value_error("h, k, and l must have the same length");
    }
    if (inputs.basis.shape[1] < 7) {
        throw py::value_error("basis must have at least 7 columns");
    }
    if (inputs.f_factors.shape[0] != inputs.basis.shape[0]
        || inputs.f_factors.shape[1] < 13) {
        throw py::value_error(
            "f_factors must have one row per basis atom and at least 13 columns"
        );
    }
    if (inputs.ref_hkl_transform.shape[0] != 3
        || inputs.ref_hkl_transform.shape[1] != 3
        || inputs.b_mat.shape[0] != 3
        || inputs.b_mat.shape[1] != 3
        || inputs.r_mat.shape[0] != 3
        || inputs.r_mat.shape[1] != 3
        || inputs.r_mat_inv.shape[0] != 3
        || inputs.r_mat_inv.shape[1] != 3) {
        throw py::value_error("transform matrices must be shape (3, 3)");
    }
    if (inputs.coherent_domain_matrix.shape[1] != 3
        || inputs.coherent_domain_matrix.shape[2] != 4) {
        throw py::value_error(
            "coherentDomainMatrix must have shape (N, 3, 4)"
        );
    }
    if (inputs.coherent_domain_occupancy.shape[0]
        != inputs.coherent_domain_matrix.shape[0]) {
        throw py::value_error(
            "coherentDomainOccupancy length must match coherentDomainMatrix"
        );
    }
    return inputs;
}

inline double get2(
    const double *data,
    const py::buffer_info &info,
    const py::ssize_t i,
    const py::ssize_t j
) {
    return data[i * info.shape[1] + j];
}

inline double get3(
    const double *data,
    const py::buffer_info &info,
    const py::ssize_t i,
    const py::ssize_t j,
    const py::ssize_t k
) {
    return data[(i * info.shape[1] + j) * info.shape[2] + k];
}

std::unique_ptr<Matrix3[]> effective_domain_matrices(const Inputs &inputs) {
    const auto n_domains = inputs.coherent_domain_matrix.shape[0];
    const double *domain = ptr(inputs.coherent_domain_matrix);
    const double *r_mat = ptr(inputs.r_mat);
    const double *r_mat_inv = ptr(inputs.r_mat_inv);
    auto matrices = std::make_unique<Matrix3[]>(n_domains);

    for (py::ssize_t d = 0; d < n_domains; ++d) {
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                double value = 0.0;
                for (int left = 0; left < 3; ++left) {
                    for (int right = 0; right < 3; ++right) {
                        value += (
                            r_mat_inv[row * 3 + left]
                            * get3(
                                domain,
                                inputs.coherent_domain_matrix,
                                d,
                                left,
                                right
                            )
                            * r_mat[right * 3 + col]
                        );
                    }
                }
                matrices[d].value[row][col] = value;
            }
        }
    }
    return matrices;
}

py::array_t<complex128> unitcell_core(
    const Array1D &h,
    const Array1D &k,
    const Array1D &l,
    const double atten,
    const bool apply_reference_transform,
    const bool apply_attenuation,
    const bool apply_bulk_lattice_sum,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &ref_hkl_transform,
    const Array2D &b_mat,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy
) {
    const Inputs inputs = request_inputs(
        h,
        k,
        l,
        basis,
        f_factors,
        ref_hkl_transform,
        b_mat,
        r_mat,
        r_mat_inv,
        coherent_domain_matrix,
        coherent_domain_occupancy
    );
    auto result = py::array_t<complex128>(inputs.h.shape[0]);
    auto result_info = result.request();
    auto *out = static_cast<complex128 *>(result_info.ptr);

    const double *h_data = ptr(inputs.h);
    const double *k_data = ptr(inputs.k);
    const double *l_data = ptr(inputs.l);
    const double *basis_data = ptr(inputs.basis);
    const double *ff_data = ptr(inputs.f_factors);
    const double *ref_data = ptr(inputs.ref_hkl_transform);
    const double *b_data = ptr(inputs.b_mat);
    const double *domain_data = ptr(inputs.coherent_domain_matrix);
    const double *occupancy_data = ptr(inputs.coherent_domain_occupancy);
    const auto domain_matrices = effective_domain_matrices(inputs);

    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    constexpr double dw_denominator = 16.0 * pi * pi;

    py::gil_scoped_release release;
    for (py::ssize_t p = 0; p < inputs.h.shape[0]; ++p) {
        double hh = h_data[p];
        double kk = k_data[p];
        double ll = l_data[p];
        if (apply_reference_transform) {
            const double h_in = hh;
            const double k_in = kk;
            const double l_in = ll;
            hh = ref_data[0] * h_in + ref_data[1] * k_in + ref_data[2] * l_in;
            kk = ref_data[3] * h_in + ref_data[4] * k_in + ref_data[5] * l_in;
            ll = ref_data[6] * h_in + ref_data[7] * k_in + ref_data[8] * l_in;
        }

        const double qx = b_data[0] * hh + b_data[1] * kk + b_data[2] * ll;
        const double qy = b_data[3] * hh + b_data[4] * kk + b_data[5] * ll;
        const double qz = b_data[6] * hh + b_data[7] * kk + b_data[8] * ll;
        const double q_para2 = qx * qx + qy * qy;
        const double q_perp2 = qz * qz;
        const double q2 = q_para2 + q_perp2;

        double amplitude_real = 0.0;
        double amplitude_imag = 0.0;
        for (py::ssize_t i = 0; i < inputs.basis.shape[0]; ++i) {
            double form_factor_real = (
                get2(ff_data, inputs.f_factors, i, 10)
                + get2(ff_data, inputs.f_factors, i, 11)
            );
            double form_factor_imag = get2(ff_data, inputs.f_factors, i, 12);
            for (int j = 0; j < 5; ++j) {
                form_factor_real += get2(ff_data, inputs.f_factors, i, j)
                    * std::exp(
                        -get2(ff_data, inputs.f_factors, i, j + 5) * q2
                    );
            }
            const double disorder_factor = std::exp(
                -(
                    get2(basis_data, inputs.basis, i, 4) * q_para2
                    + get2(basis_data, inputs.basis, i, 5) * q_perp2
                )
                / dw_denominator
            );
            const double occupancy = get2(basis_data, inputs.basis, i, 6);
            form_factor_real *= disorder_factor * occupancy;
            form_factor_imag *= disorder_factor * occupancy;

            const double x = get2(basis_data, inputs.basis, i, 1);
            const double y = get2(basis_data, inputs.basis, i, 2);
            const double z = get2(basis_data, inputs.basis, i, 3);
            for (py::ssize_t d = 0; d < inputs.coherent_domain_matrix.shape[0]; ++d) {
                const Matrix3 &mat = domain_matrices[d];
                const double x_rel = (
                    mat.value[0][0] * x
                    + mat.value[0][1] * y
                    + mat.value[0][2] * z
                    + get3(domain_data, inputs.coherent_domain_matrix, d, 0, 3)
                );
                const double y_rel = (
                    mat.value[1][0] * x
                    + mat.value[1][1] * y
                    + mat.value[1][2] * z
                    + get3(domain_data, inputs.coherent_domain_matrix, d, 1, 3)
                );
                const double z_rel = (
                    mat.value[2][0] * x
                    + mat.value[2][1] * y
                    + mat.value[2][2] * z
                    + get3(domain_data, inputs.coherent_domain_matrix, d, 2, 3)
                );
                const double phase = two_pi * (
                    hh * x_rel + kk * y_rel + ll * z_rel
                );
                double domain_real = occupancy_data[d] * std::cos(phase);
                double domain_imag = occupancy_data[d] * std::sin(phase);
                if (apply_attenuation) {
                    const double attenuation_factor = std::exp(atten * z_rel);
                    domain_real *= attenuation_factor;
                    domain_imag *= attenuation_factor;
                }
                amplitude_real += (
                    domain_real * form_factor_real
                    - domain_imag * form_factor_imag
                );
                amplitude_imag += (
                    domain_real * form_factor_imag
                    + domain_imag * form_factor_real
                );
            }
        }
        if (apply_bulk_lattice_sum) {
            const double denominator_phase = -two_pi * l_data[p];
            const double lattice_factor = std::exp(-atten);
            const double denominator_real = (
                1.0 - lattice_factor * std::cos(denominator_phase)
            );
            const double denominator_imag = (
                -lattice_factor * std::sin(denominator_phase)
            );
            const double denominator_norm = (
                denominator_real * denominator_real
                + denominator_imag * denominator_imag
            );
            const double divided_real = (
                amplitude_real * denominator_real
                + amplitude_imag * denominator_imag
            ) / denominator_norm;
            const double divided_imag = (
                amplitude_imag * denominator_real
                - amplitude_real * denominator_imag
            ) / denominator_norm;
            amplitude_real = divided_real;
            amplitude_imag = divided_imag;
        }
        out[p] = complex128{amplitude_real, amplitude_imag};
    }
    return result;
}

py::array_t<complex128> unitcell_F_uc_bulk(
    const Array1D &h,
    const Array1D &k,
    const Array1D &l,
    const double atten,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &ref_hkl_transform,
    const Array2D &b_mat,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy,
    const double /* uc_area */
) {
    return unitcell_core(
        h,
        k,
        l,
        atten,
        true,
        true,
        false,
        basis,
        f_factors,
        ref_hkl_transform,
        b_mat,
        r_mat,
        r_mat_inv,
        coherent_domain_matrix,
        coherent_domain_occupancy
    );
}

py::array_t<complex128> unitcell_F_uc_bulk_direct(
    const Array1D &h,
    const Array1D &k,
    const Array1D &l,
    const double atten,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &ref_hkl_transform,
    const Array2D &b_mat,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy,
    const double /* uc_area */
) {
    return unitcell_core(
        h,
        k,
        l,
        atten,
        false,
        true,
        false,
        basis,
        f_factors,
        ref_hkl_transform,
        b_mat,
        r_mat,
        r_mat_inv,
        coherent_domain_matrix,
        coherent_domain_occupancy
    );
}

py::array_t<complex128> unitcell_F_bulk(
    const Array1D &h,
    const Array1D &k,
    const Array1D &l,
    const double atten,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &ref_hkl_transform,
    const Array2D &b_mat,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy,
    const double /* uc_area */
) {
    return unitcell_core(
        h,
        k,
        l,
        atten,
        true,
        true,
        true,
        basis,
        f_factors,
        ref_hkl_transform,
        b_mat,
        r_mat,
        r_mat_inv,
        coherent_domain_matrix,
        coherent_domain_occupancy
    );
}

py::array_t<complex128> unitcell_F_uc(
    const Array1D &h,
    const Array1D &k,
    const Array1D &l,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &ref_hkl_transform,
    const Array2D &b_mat,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy,
    const double /* uc_area */
) {
    return unitcell_core(
        h,
        k,
        l,
        0.0,
        true,
        false,
        false,
        basis,
        f_factors,
        ref_hkl_transform,
        b_mat,
        r_mat,
        r_mat_inv,
        coherent_domain_matrix,
        coherent_domain_occupancy
    );
}

PYBIND11_MODULE(_CTRcalc_cpp, module) {
    module.doc() = "CPU C++ acceleration kernels for CTR calculations.";
    module.def("unitcell_F_uc_bulk", &unitcell_F_uc_bulk);
    module.def("unitcell_F_uc_bulk_direct", &unitcell_F_uc_bulk_direct);
    module.def("unitcell_F_bulk", &unitcell_F_bulk);
    module.def("unitcell_F_uc", &unitcell_F_uc);
}

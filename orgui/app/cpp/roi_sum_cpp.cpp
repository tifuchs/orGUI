#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using Array2D = py::array_t<double, py::array::c_style>;
using BoolArray2D = py::array_t<bool, py::array::c_style>;
using RoiArray3D = py::array_t<std::int64_t, py::array::c_style>;

void require_ndim(
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

void require_same_image_shape(
    const py::buffer_info &image,
    const py::buffer_info &other,
    const char *name
) {
    require_ndim(other, 2, name);
    if (other.shape[0] != image.shape[0] || other.shape[1] != image.shape[1]) {
        throw py::value_error(
            std::string(name) + " must have the same shape as image"
        );
    }
}

void require_roi_shape(const py::buffer_info &roi, const char *name) {
    require_ndim(roi, 3, name);
    if (roi.shape[1] != 2 || roi.shape[2] != 2) {
        throw py::value_error(
            std::string(name) + " must have shape (N, 2, 2)"
        );
    }
}

void require_counter_shape(
    const py::buffer_info &counter,
    const py::ssize_t n_rois,
    const char *name
) {
    require_ndim(counter, 2, name);
    if (counter.shape[0] != n_rois || counter.shape[1] != 4) {
        throw py::value_error(
            std::string(name) + " must have shape (N, 4)"
        );
    }
}

inline double *mutable_double_ptr(const py::buffer_info &info) {
    return static_cast<double *>(info.ptr);
}

inline const double *double_ptr(const py::buffer_info &info) {
    return static_cast<const double *>(info.ptr);
}

inline const bool *bool_ptr(const py::buffer_info &info) {
    return static_cast<const bool *>(info.ptr);
}

inline const std::int64_t *roi_ptr(const py::buffer_info &info) {
    return static_cast<const std::int64_t *>(info.ptr);
}

inline py::ssize_t index2(
    const py::buffer_info &info,
    const py::ssize_t row,
    const py::ssize_t col
) {
    return row * info.shape[1] + col;
}

inline double numpy_maximum(const double left, const double right) {
    if (std::isnan(left) || std::isnan(right)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return left > right ? left : right;
}

inline std::int64_t roi_value(
    const std::int64_t *roi,
    const py::buffer_info &info,
    const py::ssize_t n,
    const py::ssize_t axis,
    const py::ssize_t edge
) {
    return roi[(n * info.shape[1] + axis) * info.shape[2] + edge];
}

struct RoiSums {
    double image = 0.0;
    double correction = 0.0;
    double background = 0.0;
    double pixels = 0.0;
};

struct PolyFitResult {
    bool fitted = false;
    int order = -1;
    double sample_count = 0.0;
    double coefficients[6] = {0.0};
};

int polynomial_terms(const int order) {
    if (order <= 0) {
        return 1;
    }
    if (order == 1) {
        return 3;
    }
    return 6;
}

void polynomial_basis(
    const double x,
    const double y,
    const int order,
    double *basis
) {
    basis[0] = 1.0;
    if (order <= 0) {
        return;
    }
    basis[1] = x;
    basis[2] = y;
    if (order == 1) {
        return;
    }
    basis[3] = x * x;
    basis[4] = x * y;
    basis[5] = y * y;
}

bool solve_linear_system(double *matrix, double *rhs, const int n) {
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double pivot_abs = std::abs(matrix[col * n + col]);
        for (int row = col + 1; row < n; ++row) {
            const double candidate = std::abs(matrix[row * n + col]);
            if (candidate > pivot_abs) {
                pivot = row;
                pivot_abs = candidate;
            }
        }
        if (pivot_abs < 1e-12) {
            return false;
        }
        if (pivot != col) {
            for (int k = col; k < n; ++k) {
                std::swap(matrix[col * n + k], matrix[pivot * n + k]);
            }
            std::swap(rhs[col], rhs[pivot]);
        }
        const double diagonal = matrix[col * n + col];
        for (int row = col + 1; row < n; ++row) {
            const double factor = matrix[row * n + col] / diagonal;
            matrix[row * n + col] = 0.0;
            for (int k = col + 1; k < n; ++k) {
                matrix[row * n + k] -= factor * matrix[col * n + k];
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    for (int row = n - 1; row >= 0; --row) {
        double value = rhs[row];
        for (int col = row + 1; col < n; ++col) {
            value -= matrix[row * n + col] * rhs[col];
        }
        rhs[row] = value / matrix[row * n + row];
    }
    return true;
}

RoiSums sum_roi(
    const double *image_data,
    const double *correction_data,
    const double *background_data,
    const bool *mask_data,
    const py::ssize_t image_width,
    const std::int64_t *roi,
    const py::buffer_info &roi_info,
    const py::ssize_t n
) {
    const std::int64_t x0 = roi_value(roi, roi_info, n, 0, 0);
    const std::int64_t x1 = roi_value(roi, roi_info, n, 0, 1);
    const std::int64_t y0 = roi_value(roi, roi_info, n, 1, 0);
    const std::int64_t y1 = roi_value(roi, roi_info, n, 1, 1);
    RoiSums sums;
    for (std::int64_t y = y0; y < y1; ++y) {
        const py::ssize_t row_offset = y * image_width;
        for (std::int64_t x = x0; x < x1; ++x) {
            const py::ssize_t idx = row_offset + x;
            if (!mask_data[idx]) {
                sums.pixels += 1.0;
                const double image_value = image_data[idx];
                if (!std::isnan(image_value)) {
                    sums.image += image_value;
                }
            }

            const double correction_value = correction_data[idx];
            if (!std::isnan(correction_value)) {
                sums.correction += correction_value;
            }

            if (background_data != nullptr) {
                const double background_value = background_data[idx];
                if (!std::isnan(background_value)) {
                    sums.background += background_value;
                }
            }
        }
    }
    return sums;
}

void accumulate_fit_roi(
    const double *image_data,
    const bool *mask_data,
    const py::ssize_t image_width,
    const std::int64_t *roi,
    const py::buffer_info &roi_info,
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const int order,
    double *normal_matrix,
    double *normal_rhs,
    double &sample_count
) {
    const int n_terms = polynomial_terms(order);
    double basis[6];
    const std::int64_t x0 = roi_value(roi, roi_info, n, 0, 0);
    const std::int64_t x1 = roi_value(roi, roi_info, n, 0, 1);
    const std::int64_t y0 = roi_value(roi, roi_info, n, 1, 0);
    const std::int64_t y1 = roi_value(roi, roi_info, n, 1, 1);
    for (std::int64_t y = y0; y < y1; ++y) {
        const py::ssize_t row_offset = y * image_width;
        const double yn = (static_cast<double>(y) - y_origin) / scale;
        for (std::int64_t x = x0; x < x1; ++x) {
            const py::ssize_t idx = row_offset + x;
            if (mask_data[idx]) {
                continue;
            }
            const double value = image_data[idx];
            if (std::isnan(value)) {
                continue;
            }
            const double xn = (static_cast<double>(x) - x_origin) / scale;
            polynomial_basis(xn, yn, order, basis);
            for (int row = 0; row < n_terms; ++row) {
                normal_rhs[row] += basis[row] * value;
                for (int col = 0; col < n_terms; ++col) {
                    normal_matrix[row * n_terms + col] += (
                        basis[row] * basis[col]
                    );
                }
            }
            sample_count += 1.0;
        }
    }
}

double evaluate_fit_roi(
    const bool *mask_data,
    const py::ssize_t image_width,
    const std::int64_t *roi,
    const py::buffer_info &roi_info,
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const int order,
    const double *coefficients
) {
    const int n_terms = polynomial_terms(order);
    double basis[6];
    double total = 0.0;
    const std::int64_t x0 = roi_value(roi, roi_info, n, 0, 0);
    const std::int64_t x1 = roi_value(roi, roi_info, n, 0, 1);
    const std::int64_t y0 = roi_value(roi, roi_info, n, 1, 0);
    const std::int64_t y1 = roi_value(roi, roi_info, n, 1, 1);
    for (std::int64_t y = y0; y < y1; ++y) {
        const py::ssize_t row_offset = y * image_width;
        const double yn = (static_cast<double>(y) - y_origin) / scale;
        for (std::int64_t x = x0; x < x1; ++x) {
            const py::ssize_t idx = row_offset + x;
            if (mask_data[idx]) {
                continue;
            }
            const double xn = (static_cast<double>(x) - x_origin) / scale;
            polynomial_basis(xn, yn, order, basis);
            double value = 0.0;
            for (int term = 0; term < n_terms; ++term) {
                value += coefficients[term] * basis[term];
            }
            total += value;
        }
    }
    return total;
}

void fill_fit_roi_patch(
    const py::ssize_t image_width,
    const std::int64_t *roi,
    const py::buffer_info &roi_info,
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const int order,
    const double *coefficients,
    double *patch_data,
    const py::ssize_t patch_width
) {
    const int n_terms = polynomial_terms(order);
    double basis[6];
    const std::int64_t x0 = roi_value(roi, roi_info, n, 0, 0);
    const std::int64_t x1 = roi_value(roi, roi_info, n, 0, 1);
    const std::int64_t y0 = roi_value(roi, roi_info, n, 1, 0);
    const std::int64_t y1 = roi_value(roi, roi_info, n, 1, 1);
    for (std::int64_t y = y0; y < y1; ++y) {
        const py::ssize_t patch_row_offset = (y - y0) * patch_width;
        const double yn = (static_cast<double>(y) - y_origin) / scale;
        for (std::int64_t x = x0; x < x1; ++x) {
            const double xn = (static_cast<double>(x) - x_origin) / scale;
            polynomial_basis(xn, yn, order, basis);
            double value = 0.0;
            for (int term = 0; term < n_terms; ++term) {
                value += coefficients[term] * basis[term];
            }
            patch_data[patch_row_offset + (x - x0)] = value;
        }
    }
}

double residual_sse_fit_roi(
    const double *image_data,
    const bool *mask_data,
    const py::ssize_t image_width,
    const std::int64_t *roi,
    const py::buffer_info &roi_info,
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const int order,
    const double *coefficients
) {
    const int n_terms = polynomial_terms(order);
    double basis[6];
    double sse = 0.0;
    const std::int64_t x0 = roi_value(roi, roi_info, n, 0, 0);
    const std::int64_t x1 = roi_value(roi, roi_info, n, 0, 1);
    const std::int64_t y0 = roi_value(roi, roi_info, n, 1, 0);
    const std::int64_t y1 = roi_value(roi, roi_info, n, 1, 1);
    for (std::int64_t y = y0; y < y1; ++y) {
        const py::ssize_t row_offset = y * image_width;
        const double yn = (static_cast<double>(y) - y_origin) / scale;
        for (std::int64_t x = x0; x < x1; ++x) {
            const py::ssize_t idx = row_offset + x;
            if (mask_data[idx] || std::isnan(image_data[idx])) {
                continue;
            }
            const double xn = (static_cast<double>(x) - x_origin) / scale;
            polynomial_basis(xn, yn, order, basis);
            double fitted_value = 0.0;
            for (int term = 0; term < n_terms; ++term) {
                fitted_value += coefficients[term] * basis[term];
            }
            const double residual = image_data[idx] - fitted_value;
            sse += residual * residual;
        }
    }
    return sse;
}

PolyFitResult fit_local_polynomial_background(
    const double *image_data,
    const bool *mask_data,
    const py::ssize_t image_width,
    const std::int64_t *left_roi,
    const py::buffer_info &left_info,
    const std::int64_t *right_roi,
    const py::buffer_info &right_info,
    const std::int64_t *top_roi,
    const py::buffer_info &top_info,
    const std::int64_t *bottom_roi,
    const py::buffer_info &bottom_info,
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const int fit_order
) {
    PolyFitResult result;
    int order = fit_order;
    while (order >= 0 && !result.fitted) {
        const int n_terms = polynomial_terms(order);
        double normal_matrix[36] = {0.0};
        double normal_rhs[6] = {0.0};
        double sample_count = 0.0;
        accumulate_fit_roi(
            image_data,
            mask_data,
            image_width,
            left_roi,
            left_info,
            n,
            x_origin,
            y_origin,
            scale,
            order,
            normal_matrix,
            normal_rhs,
            sample_count
        );
        accumulate_fit_roi(
            image_data,
            mask_data,
            image_width,
            right_roi,
            right_info,
            n,
            x_origin,
            y_origin,
            scale,
            order,
            normal_matrix,
            normal_rhs,
            sample_count
        );
        accumulate_fit_roi(
            image_data,
            mask_data,
            image_width,
            top_roi,
            top_info,
            n,
            x_origin,
            y_origin,
            scale,
            order,
            normal_matrix,
            normal_rhs,
            sample_count
        );
        accumulate_fit_roi(
            image_data,
            mask_data,
            image_width,
            bottom_roi,
            bottom_info,
            n,
            x_origin,
            y_origin,
            scale,
            order,
            normal_matrix,
            normal_rhs,
            sample_count
        );
        if (sample_count >= n_terms
            && solve_linear_system(normal_matrix, normal_rhs, n_terms)) {
            result.fitted = true;
            result.order = order;
            result.sample_count = sample_count;
            for (int term = 0; term < n_terms; ++term) {
                result.coefficients[term] = normal_rhs[term];
            }
        } else {
            --order;
        }
    }
    return result;
}

void fill_image_and_correction_counters(
    const double *image_data,
    const py::buffer_info &image_info,
    const bool *mask_data,
    const py::buffer_info &mask_info,
    const double *correction_data,
    const py::buffer_info &correction_info,
    const double *background_data,
    const std::int64_t *center_roi,
    const py::buffer_info &center_info,
    const std::int64_t *left_roi,
    const py::buffer_info &left_info,
    const std::int64_t *right_roi,
    const py::buffer_info &right_info,
    const std::int64_t *top_roi,
    const py::buffer_info &top_info,
    const std::int64_t *bottom_roi,
    const py::buffer_info &bottom_info,
    double *all_counters,
    const py::buffer_info &all_info,
    double *correction_counters,
    const py::buffer_info &correction_counter_info,
    double *background_counters,
    const py::buffer_info *background_counter_info
) {
    const py::ssize_t n_rois = center_info.shape[0];
    const py::ssize_t image_width = image_info.shape[1];
    for (py::ssize_t i = 0; i < n_rois; ++i) {
        const RoiSums center = sum_roi(
            image_data,
            correction_data,
            background_data,
            mask_data,
            image_width,
            center_roi,
            center_info,
            i
        );
        const RoiSums left = sum_roi(
            image_data,
            correction_data,
            background_data,
            mask_data,
            image_width,
            left_roi,
            left_info,
            i
        );
        const RoiSums right = sum_roi(
            image_data,
            correction_data,
            background_data,
            mask_data,
            image_width,
            right_roi,
            right_info,
            i
        );
        const RoiSums top = sum_roi(
            image_data,
            correction_data,
            background_data,
            mask_data,
            image_width,
            top_roi,
            top_info,
            i
        );
        const RoiSums bottom = sum_roi(
            image_data,
            correction_data,
            background_data,
            mask_data,
            image_width,
            bottom_roi,
            bottom_info,
            i
        );
        const double image_background = (
            left.image + right.image + top.image + bottom.image
        );
        const double correction_background = (
            left.correction + right.correction + top.correction
            + bottom.correction
        );
        const double background_pixels = (
            left.pixels + right.pixels + top.pixels + bottom.pixels
        );

        all_counters[i * all_info.shape[1] + 0] = center.image;
        all_counters[i * all_info.shape[1] + 1] = center.pixels;
        all_counters[i * all_info.shape[1] + 2] = image_background;
        all_counters[i * all_info.shape[1] + 3] = background_pixels;
        correction_counters[i * correction_counter_info.shape[1] + 0] = (
            center.correction
        );
        correction_counters[i * correction_counter_info.shape[1] + 1] = (
            center.pixels
        );
        correction_counters[i * correction_counter_info.shape[1] + 2] = (
            correction_background
        );
        correction_counters[i * correction_counter_info.shape[1] + 3] = (
            background_pixels
        );
        if (background_data != nullptr) {
            background_counters[i * background_counter_info->shape[1] + 0] = (
                center.background
            );
            background_counters[i * background_counter_info->shape[1] + 1] = (
                center.pixels
            );
            background_counters[i * background_counter_info->shape[1] + 2] = (
                left.background + right.background + top.background
                + bottom.background
            );
            background_counters[i * background_counter_info->shape[1] + 3] = (
                background_pixels
            );
        }
    }
}

void processImage_Carr(
    Array2D image,
    const BoolArray2D mask,
    const Array2D correction,
    const RoiArray3D center_roi,
    const RoiArray3D left_roi,
    const RoiArray3D right_roi,
    const RoiArray3D top_roi,
    const RoiArray3D bottom_roi,
    Array2D all_counters,
    Array2D correction_counters
) {
    auto image_info = image.request();
    auto mask_info = mask.request();
    auto correction_info = correction.request();
    auto center_info = center_roi.request();
    auto left_info = left_roi.request();
    auto right_info = right_roi.request();
    auto top_info = top_roi.request();
    auto bottom_info = bottom_roi.request();
    auto all_info = all_counters.request();
    auto correction_counter_info = correction_counters.request();

    require_ndim(image_info, 2, "image");
    require_ndim(mask_info, 2, "mask");
    require_ndim(correction_info, 2, "C_corr");
    require_same_image_shape(image_info, mask_info, "mask");
    require_same_image_shape(image_info, correction_info, "C_corr");
    require_roi_shape(center_info, "croi");
    require_roi_shape(left_info, "leftroi");
    require_roi_shape(right_info, "rightroi");
    require_roi_shape(top_info, "toproi");
    require_roi_shape(bottom_info, "bottomroi");
    const py::ssize_t n_rois = center_info.shape[0];
    if (
        left_info.shape[0] != n_rois
        || right_info.shape[0] != n_rois
        || top_info.shape[0] != n_rois
        || bottom_info.shape[0] != n_rois
    ) {
        throw py::value_error("all ROI arrays must have the same length");
    }
    require_counter_shape(all_info, n_rois, "all_counters");
    require_counter_shape(correction_counter_info, n_rois, "all_Carr");

    double *image_data = mutable_double_ptr(image_info);
    const bool *mask_data = bool_ptr(mask_info);
    py::gil_scoped_release release;
    fill_image_and_correction_counters(
        image_data,
        image_info,
        mask_data,
        mask_info,
        double_ptr(correction_info),
        correction_info,
        nullptr,
        roi_ptr(center_info),
        center_info,
        roi_ptr(left_info),
        left_info,
        roi_ptr(right_info),
        right_info,
        roi_ptr(top_info),
        top_info,
        roi_ptr(bottom_info),
        bottom_info,
        mutable_double_ptr(all_info),
        all_info,
        mutable_double_ptr(correction_counter_info),
        correction_counter_info,
        nullptr,
        nullptr
    );
}

void processImage_bg_Carr(
    Array2D image,
    const Array2D background,
    const BoolArray2D mask,
    const Array2D correction,
    const RoiArray3D center_roi,
    const RoiArray3D left_roi,
    const RoiArray3D right_roi,
    const RoiArray3D top_roi,
    const RoiArray3D bottom_roi,
    Array2D all_counters,
    Array2D correction_counters,
    Array2D background_counters
) {
    auto image_info = image.request();
    auto background_info = background.request();
    auto mask_info = mask.request();
    auto correction_info = correction.request();
    auto center_info = center_roi.request();
    auto left_info = left_roi.request();
    auto right_info = right_roi.request();
    auto top_info = top_roi.request();
    auto bottom_info = bottom_roi.request();
    auto all_info = all_counters.request();
    auto correction_counter_info = correction_counters.request();
    auto background_counter_info = background_counters.request();

    require_ndim(image_info, 2, "image");
    require_ndim(background_info, 2, "bg");
    require_ndim(mask_info, 2, "mask");
    require_ndim(correction_info, 2, "C_corr");
    require_same_image_shape(image_info, background_info, "bg");
    require_same_image_shape(image_info, mask_info, "mask");
    require_same_image_shape(image_info, correction_info, "C_corr");
    require_roi_shape(center_info, "croi");
    require_roi_shape(left_info, "leftroi");
    require_roi_shape(right_info, "rightroi");
    require_roi_shape(top_info, "toproi");
    require_roi_shape(bottom_info, "bottomroi");
    const py::ssize_t n_rois = center_info.shape[0];
    if (
        left_info.shape[0] != n_rois
        || right_info.shape[0] != n_rois
        || top_info.shape[0] != n_rois
        || bottom_info.shape[0] != n_rois
    ) {
        throw py::value_error("all ROI arrays must have the same length");
    }
    require_counter_shape(all_info, n_rois, "all_counters");
    require_counter_shape(correction_counter_info, n_rois, "all_Carr");
    require_counter_shape(background_counter_info, n_rois, "all_Bgimg");

    double *image_data = mutable_double_ptr(image_info);
    const bool *mask_data = bool_ptr(mask_info);
    double *background_counter_data = mutable_double_ptr(
        background_counter_info
    );
    py::gil_scoped_release release;
    fill_image_and_correction_counters(
        image_data,
        image_info,
        mask_data,
        mask_info,
        double_ptr(correction_info),
        correction_info,
        double_ptr(background_info),
        roi_ptr(center_info),
        center_info,
        roi_ptr(left_info),
        left_info,
        roi_ptr(right_info),
        right_info,
        roi_ptr(top_info),
        top_info,
        roi_ptr(bottom_info),
        bottom_info,
        mutable_double_ptr(all_info),
        all_info,
        mutable_double_ptr(correction_counter_info),
        correction_counter_info,
        background_counter_data,
        &background_counter_info
    );
}

void processImage_polybg_Carr(
    const Array2D image,
    const BoolArray2D mask,
    const Array2D correction,
    const RoiArray3D center_roi,
    const RoiArray3D left_roi,
    const RoiArray3D right_roi,
    const RoiArray3D top_roi,
    const RoiArray3D bottom_roi,
    Array2D all_counters,
    Array2D correction_counters,
    int fit_order
) {
    if (fit_order < 0 || fit_order > 2) {
        throw py::value_error("fit_order must be 0, 1, or 2");
    }
    auto image_info = image.request();
    auto mask_info = mask.request();
    auto correction_info = correction.request();
    auto center_info = center_roi.request();
    auto left_info = left_roi.request();
    auto right_info = right_roi.request();
    auto top_info = top_roi.request();
    auto bottom_info = bottom_roi.request();
    auto all_info = all_counters.request();
    auto correction_counter_info = correction_counters.request();

    require_ndim(image_info, 2, "image");
    require_ndim(mask_info, 2, "mask");
    require_ndim(correction_info, 2, "C_corr");
    require_same_image_shape(image_info, mask_info, "mask");
    require_same_image_shape(image_info, correction_info, "C_corr");
    require_roi_shape(center_info, "croi");
    require_roi_shape(left_info, "leftroi");
    require_roi_shape(right_info, "rightroi");
    require_roi_shape(top_info, "toproi");
    require_roi_shape(bottom_info, "bottomroi");
    const py::ssize_t n_rois = center_info.shape[0];
    if (
        left_info.shape[0] != n_rois
        || right_info.shape[0] != n_rois
        || top_info.shape[0] != n_rois
        || bottom_info.shape[0] != n_rois
    ) {
        throw py::value_error("all ROI arrays must have the same length");
    }
    require_counter_shape(all_info, n_rois, "all_counters");
    require_counter_shape(correction_counter_info, n_rois, "all_Carr");

    const double *image_data = double_ptr(image_info);
    const double *correction_data = double_ptr(correction_info);
    const bool *mask_data = bool_ptr(mask_info);
    double *all_data = mutable_double_ptr(all_info);
    double *correction_counter_data = mutable_double_ptr(
        correction_counter_info
    );
    const py::ssize_t image_width = image_info.shape[1];
    const py::ssize_t counter_width = all_info.shape[1];
    const py::ssize_t correction_counter_width = correction_counter_info.shape[1];
    py::gil_scoped_release release;
    for (py::ssize_t i = 0; i < n_rois; ++i) {
        const RoiSums center = sum_roi(
            image_data,
            correction_data,
            nullptr,
            mask_data,
            image_width,
            roi_ptr(center_info),
            center_info,
            i
        );
        const RoiSums left = sum_roi(
            image_data,
            correction_data,
            nullptr,
            mask_data,
            image_width,
            roi_ptr(left_info),
            left_info,
            i
        );
        const RoiSums right = sum_roi(
            image_data,
            correction_data,
            nullptr,
            mask_data,
            image_width,
            roi_ptr(right_info),
            right_info,
            i
        );
        const RoiSums top = sum_roi(
            image_data,
            correction_data,
            nullptr,
            mask_data,
            image_width,
            roi_ptr(top_info),
            top_info,
            i
        );
        const RoiSums bottom = sum_roi(
            image_data,
            correction_data,
            nullptr,
            mask_data,
            image_width,
            roi_ptr(bottom_info),
            bottom_info,
            i
        );

        const double x0 = static_cast<double>(
            roi_value(roi_ptr(center_info), center_info, i, 0, 0)
        );
        const double x1 = static_cast<double>(
            roi_value(roi_ptr(center_info), center_info, i, 0, 1)
        );
        const double y0 = static_cast<double>(
            roi_value(roi_ptr(center_info), center_info, i, 1, 0)
        );
        const double y1 = static_cast<double>(
            roi_value(roi_ptr(center_info), center_info, i, 1, 1)
        );
        const double x_origin = 0.5 * (x0 + x1 - 1.0);
        const double y_origin = 0.5 * (y0 + y1 - 1.0);
        const double scale = std::max(
            1.0,
            std::max(std::abs(x1 - x0), std::abs(y1 - y0))
        );

        const PolyFitResult fit = fit_local_polynomial_background(
            image_data,
            mask_data,
            image_width,
            roi_ptr(left_info),
            left_info,
            roi_ptr(right_info),
            right_info,
            roi_ptr(top_info),
            top_info,
            roi_ptr(bottom_info),
            bottom_info,
            i,
            x_origin,
            y_origin,
            scale,
            fit_order
        );
        double fitted_background = 0.0;
        if (fit.fitted) {
            fitted_background = evaluate_fit_roi(
                mask_data,
                image_width,
                roi_ptr(center_info),
                center_info,
                i,
                x_origin,
                y_origin,
                scale,
                fit.order,
                fit.coefficients
            );
        }

        double stored_background = 0.0;
        if (fit.fitted && center.pixels > 0.0) {
            // Preserve the legacy background-subtraction formula:
            // cpixels / bgpixels * bgroi == fitted center background.
            stored_background = fitted_background * fit.sample_count / center.pixels;
        }
        all_data[i * counter_width + 0] = center.image;
        all_data[i * counter_width + 1] = center.pixels;
        all_data[i * counter_width + 2] = stored_background;
        all_data[i * counter_width + 3] = fit.fitted ? fit.sample_count : 0.0;
        correction_counter_data[i * correction_counter_width + 0] = (
            center.correction
        );
        correction_counter_data[i * correction_counter_width + 1] = (
            center.pixels
        );
        correction_counter_data[i * correction_counter_width + 2] = (
            left.correction + right.correction + top.correction
            + bottom.correction
        );
        correction_counter_data[i * correction_counter_width + 3] = (
            left.pixels + right.pixels + top.pixels + bottom.pixels
        );
    }
}

py::tuple interpolate_polybg_croi(
    const Array2D image,
    const BoolArray2D mask,
    const RoiArray3D center_roi,
    const RoiArray3D left_roi,
    const RoiArray3D right_roi,
    const RoiArray3D top_roi,
    const RoiArray3D bottom_roi,
    int fit_order
) {
    if (fit_order < 0 || fit_order > 2) {
        throw py::value_error("fit_order must be 0, 1, or 2");
    }
    auto image_info = image.request();
    auto mask_info = mask.request();
    auto center_info = center_roi.request();
    auto left_info = left_roi.request();
    auto right_info = right_roi.request();
    auto top_info = top_roi.request();
    auto bottom_info = bottom_roi.request();

    require_ndim(image_info, 2, "image");
    require_ndim(mask_info, 2, "mask");
    require_same_image_shape(image_info, mask_info, "mask");
    require_roi_shape(center_info, "croi");
    require_roi_shape(left_info, "leftroi");
    require_roi_shape(right_info, "rightroi");
    require_roi_shape(top_info, "toproi");
    require_roi_shape(bottom_info, "bottomroi");
    if (
        center_info.shape[0] < 1
        || left_info.shape[0] < 1
        || right_info.shape[0] < 1
        || top_info.shape[0] < 1
        || bottom_info.shape[0] < 1
    ) {
        throw py::value_error("all ROI arrays must contain at least one ROI");
    }

    const std::int64_t x0 = roi_value(roi_ptr(center_info), center_info, 0, 0, 0);
    const std::int64_t x1 = roi_value(roi_ptr(center_info), center_info, 0, 0, 1);
    const std::int64_t y0 = roi_value(roi_ptr(center_info), center_info, 0, 1, 0);
    const std::int64_t y1 = roi_value(roi_ptr(center_info), center_info, 0, 1, 1);
    if (x1 < x0 || y1 < y0) {
        throw py::value_error("center ROI edges must be sorted");
    }

    py::array_t<double> patch({y1 - y0, x1 - x0});
    auto patch_info = patch.request();
    double *patch_data = mutable_double_ptr(patch_info);
    for (py::ssize_t i = 0; i < patch_info.size; ++i) {
        patch_data[i] = std::numeric_limits<double>::quiet_NaN();
    }

    const double x_origin = 0.5 * (static_cast<double>(x0 + x1) - 1.0);
    const double y_origin = 0.5 * (static_cast<double>(y0 + y1) - 1.0);
    const double scale = std::max(
        1.0,
        std::max(std::abs(static_cast<double>(x1 - x0)),
                 std::abs(static_cast<double>(y1 - y0)))
    );
    const double *image_data = double_ptr(image_info);
    const bool *mask_data = bool_ptr(mask_info);
    const py::ssize_t image_width = image_info.shape[1];
    const py::ssize_t patch_width = patch_info.shape[1];

    bool success = false;
    int used_order = -1;
    double sample_count = 0.0;
    double rmse = std::numeric_limits<double>::quiet_NaN();
    {
        py::gil_scoped_release release;
        const PolyFitResult fit = fit_local_polynomial_background(
            image_data,
            mask_data,
            image_width,
            roi_ptr(left_info),
            left_info,
            roi_ptr(right_info),
            right_info,
            roi_ptr(top_info),
            top_info,
            roi_ptr(bottom_info),
            bottom_info,
            0,
            x_origin,
            y_origin,
            scale,
            fit_order
        );
        if (fit.fitted) {
            fill_fit_roi_patch(
                image_width,
                roi_ptr(center_info),
                center_info,
                0,
                x_origin,
                y_origin,
                scale,
                fit.order,
                fit.coefficients,
                patch_data,
                patch_width
            );
            const double sse = (
                residual_sse_fit_roi(
                    image_data,
                    mask_data,
                    image_width,
                    roi_ptr(left_info),
                    left_info,
                    0,
                    x_origin,
                    y_origin,
                    scale,
                    fit.order,
                    fit.coefficients
                )
                + residual_sse_fit_roi(
                    image_data,
                    mask_data,
                    image_width,
                    roi_ptr(right_info),
                    right_info,
                    0,
                    x_origin,
                    y_origin,
                    scale,
                    fit.order,
                    fit.coefficients
                )
                + residual_sse_fit_roi(
                    image_data,
                    mask_data,
                    image_width,
                    roi_ptr(top_info),
                    top_info,
                    0,
                    x_origin,
                    y_origin,
                    scale,
                    fit.order,
                    fit.coefficients
                )
                + residual_sse_fit_roi(
                    image_data,
                    mask_data,
                    image_width,
                    roi_ptr(bottom_info),
                    bottom_info,
                    0,
                    x_origin,
                    y_origin,
                    scale,
                    fit.order,
                    fit.coefficients
                )
            );
            const int n_terms = polynomial_terms(fit.order);
            const double dof = std::max(1.0, fit.sample_count - n_terms);
            rmse = std::sqrt(sse / dof);
            success = true;
            used_order = fit.order;
            sample_count = fit.sample_count;
        }
    }
    py::dict stats;
    stats["success"] = success;
    stats["order"] = used_order;
    stats["sample_count"] = sample_count;
    stats["rmse"] = rmse;
    return py::make_tuple(patch, stats);
}

void calcMaxSum(const Array2D image, Array2D sum_image, Array2D max_image) {
    auto image_info = image.request();
    auto sum_info = sum_image.request();
    auto max_info = max_image.request();
    require_ndim(image_info, 2, "image");
    require_same_image_shape(image_info, sum_info, "sumimg");
    require_same_image_shape(image_info, max_info, "maximg");
    const double *image_data = double_ptr(image_info);
    double *sum_data = mutable_double_ptr(sum_info);
    double *max_data = mutable_double_ptr(max_info);
    py::gil_scoped_release release;
    for (py::ssize_t row = 0; row < image_info.shape[0]; ++row) {
        for (py::ssize_t col = 0; col < image_info.shape[1]; ++col) {
            const auto idx = index2(image_info, row, col);
            sum_data[idx] += image_data[idx];
            max_data[idx] = numpy_maximum(image_data[idx], max_data[idx]);
        }
    }
}

void calcBgSub(Array2D image, const Array2D background) {
    auto image_info = image.request();
    auto background_info = background.request();
    require_ndim(image_info, 2, "image");
    require_same_image_shape(image_info, background_info, "bg");
    double *image_data = mutable_double_ptr(image_info);
    const double *background_data = double_ptr(background_info);
    py::gil_scoped_release release;
    for (py::ssize_t row = 0; row < image_info.shape[0]; ++row) {
        for (py::ssize_t col = 0; col < image_info.shape[1]; ++col) {
            const auto idx = index2(image_info, row, col);
            image_data[idx] -= background_data[idx];
        }
    }
}

void calcMaxSum_bg(
    const Array2D image,
    Array2D sum_image,
    Array2D max_image,
    const Array2D background
) {
    auto image_info = image.request();
    auto sum_info = sum_image.request();
    auto max_info = max_image.request();
    auto background_info = background.request();
    require_ndim(image_info, 2, "image");
    require_same_image_shape(image_info, sum_info, "sumimg");
    require_same_image_shape(image_info, max_info, "maximg");
    require_same_image_shape(image_info, background_info, "bg");
    const double *image_data = double_ptr(image_info);
    const double *background_data = double_ptr(background_info);
    double *sum_data = mutable_double_ptr(sum_info);
    double *max_data = mutable_double_ptr(max_info);
    py::gil_scoped_release release;
    for (py::ssize_t row = 0; row < image_info.shape[0]; ++row) {
        for (py::ssize_t col = 0; col < image_info.shape[1]; ++col) {
            const auto idx = index2(image_info, row, col);
            const double value = image_data[idx] - background_data[idx];
            sum_data[idx] += value;
            max_data[idx] = numpy_maximum(value, max_data[idx]);
        }
    }
}

PYBIND11_MODULE(_roi_sum_cpp, module) {
    module.doc() = "CPU C++ acceleration kernels for ROI image summing.";
    module.def("processImage_Carr", &processImage_Carr);
    module.def("processImage_bg_Carr", &processImage_bg_Carr);
    module.def(
        "processImage_polybg_Carr",
        &processImage_polybg_Carr,
        py::arg("image"),
        py::arg("mask"),
        py::arg("C_corr"),
        py::arg("croi"),
        py::arg("leftroi"),
        py::arg("rightroi"),
        py::arg("toproi"),
        py::arg("bottomroi"),
        py::arg("all_counters"),
        py::arg("all_Carr"),
        py::arg("fit_order") = 1
    );
    module.def(
        "interpolate_polybg_croi",
        &interpolate_polybg_croi,
        py::arg("image"),
        py::arg("mask"),
        py::arg("croi"),
        py::arg("leftroi"),
        py::arg("rightroi"),
        py::arg("toproi"),
        py::arg("bottomroi"),
        py::arg("fit_order") = 1
    );
    module.def("calcMaxSum", &calcMaxSum);
    module.def("calcBgSub", &calcBgSub);
    module.def("calcMaxSum_bg", &calcMaxSum_bg);
}

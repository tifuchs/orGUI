// Linear influence coefficients for unweighted local polynomial fits.
// Positions are pixels; input values are photon counts and moments counts^2.

using PolynomialWeights = std::vector<std::pair<py::ssize_t, double>>;

double polynomial_coordinate_scale(
    const py::buffer_info &center,
    const py::buffer_info *const background[4],
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const int order
) {
    double scale = std::max(1.0, std::max(
        static_cast<double>(roi_value(roi_ptr(center), center, n, 0, 1)
                            - roi_value(roi_ptr(center), center, n, 0, 0)),
        static_cast<double>(roi_value(roi_ptr(center), center, n, 1, 1)
                            - roi_value(roi_ptr(center), center, n, 1, 0))));
    if (order <= 2) {
        return scale;  // Preserve released degree 0/1/2 numerical coordinates.
    }
    for (const auto info : {&center, background[0], background[1],
                            background[2], background[3]}) {
        for (int axis = 0; axis < 2; ++axis) {
            const double origin = axis == 0 ? x_origin : y_origin;
            const double low = static_cast<double>(
                roi_value(roi_ptr(*info), *info, n, axis, 0));
            const double high = static_cast<double>(
                roi_value(roi_ptr(*info), *info, n, axis, 1));
            if (high > low) {
                scale = std::max(scale, std::max(std::abs(low - origin),
                                                std::abs(high - 1.0 - origin)));
            }
        }
    }
    return scale;
}

// Polynomial influence function for one occurrence of an ORIGINAL fit sample
// in the stored background counter B. Coordinates use the fit normalization.
bool polynomial_prediction_coefficients(
    const double *image_data,
    const bool *mask_data,
    const py::ssize_t width,
    const py::buffer_info &center,
    const py::buffer_info *const background[4],
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const PolyFitResult &fit,
    const double center_pixels,
    double *prediction,
    const bool *repaired_data = nullptr
) {
    if (!fit.fitted || center_pixels <= 0.0) {
        return false;
    }
    const int terms = polynomial_terms(fit.order);
    double basis[15];
    std::fill(prediction, prediction + terms, 0.0);
    const auto cx0 = roi_value(roi_ptr(center), center, n, 0, 0);
    const auto cx1 = roi_value(roi_ptr(center), center, n, 0, 1);
    const auto cy0 = roi_value(roi_ptr(center), center, n, 1, 0);
    const auto cy1 = roi_value(roi_ptr(center), center, n, 1, 1);
    for (auto y = cy0; y < cy1; ++y) {
        for (auto x = cx0; x < cx1; ++x) {
            const auto index = y * width + x;
            if (mask_data[index] && !(repaired_data && repaired_data[index])) {
                continue;
            }
            polynomial_basis((x - x_origin) / scale, (y - y_origin) / scale,
                             fit.order, basis);
            for (int j = 0; j < terms; ++j) {
                prediction[j] += basis[j];
            }
        }
    }
    double normal[225];
    std::copy(fit.normal_matrix, fit.normal_matrix + terms * terms, normal);
    if (!solve_polynomial_system(normal, prediction, fit.order)) {
        throw std::runtime_error("Cannot solve polynomial prediction weights");
    }
    const double stored_scale = fit.sample_count / center_pixels;
    for (int j = 0; j < terms; ++j) {
        prediction[j] *= stored_scale;
    }
    return true;
}

// Each returned coefficient multiplies one ORIGINAL detector pixel in the
// stored background counter B (not its center-normalized prediction).
// Repeated pixels in overlapping strips are combined before taking moments.
// A repaired mask changes only the prediction footprint, never the fit sources.
PolynomialWeights polynomial_background_weights(
    const double *image_data,
    const bool *mask_data,
    const py::ssize_t width,
    const py::buffer_info &center,
    const py::buffer_info *const background[4],
    const py::ssize_t n,
    const double x_origin,
    const double y_origin,
    const double scale,
    const PolyFitResult &fit,
    const double center_pixels,
    const bool *repaired_data = nullptr
) {
    PolynomialWeights weights;
    double prediction[15];
    if (!polynomial_prediction_coefficients(
            image_data, mask_data, width, center, background, n,
            x_origin, y_origin, scale, fit, center_pixels, prediction, repaired_data)) {
        return weights;
    }
    const int terms = polynomial_terms(fit.order);
    double basis[15];
    for (const auto info : {background[0], background[1],
                            background[2], background[3]}) {
        const auto x0 = roi_value(roi_ptr(*info), *info, n, 0, 0);
        const auto x1 = roi_value(roi_ptr(*info), *info, n, 0, 1);
        const auto y0 = roi_value(roi_ptr(*info), *info, n, 1, 0);
        const auto y1 = roi_value(roi_ptr(*info), *info, n, 1, 1);
        for (auto y = y0; y < y1; ++y) {
            for (auto x = x0; x < x1; ++x) {
                const auto index = y * width + x;
                if (mask_data[index]
                    || (image_data && std::isnan(image_data[index]))) {
                    continue;
                }
                polynomial_basis((x - x_origin) / scale, (y - y_origin) / scale,
                                 fit.order, basis);
                double weight = 0.0;
                for (int j = 0; j < terms; ++j) {
                    weight += basis[j] * prediction[j];
                }
                weights.emplace_back(index, weight);
            }
        }
    }
    std::sort(weights.begin(), weights.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    std::size_t written = 0;
    for (const auto &entry : weights) {
        if (written && weights[written - 1].first == entry.first) {
            weights[written - 1].second += entry.second;
        } else {
            weights[written++] = entry;
        }
    }
    weights.resize(written);
    return weights;
}

void polynomial_counter_moments(
    const double *image_data,
    const py::ssize_t width,
    const py::buffer_info &center,
    const py::ssize_t n,
    const PolynomialWeights &weights,
    const double original_center_count,
    const double center_pixels,
    const double background_pixels,
    double *moments
) {
    double background_variance = 0.0;
    double covariance = 0.0;
    const auto x0 = roi_value(roi_ptr(center), center, n, 0, 0);
    const auto x1 = roi_value(roi_ptr(center), center, n, 0, 1);
    const auto y0 = roi_value(roi_ptr(center), center, n, 1, 0);
    const auto y1 = roi_value(roi_ptr(center), center, n, 1, 1);
    for (const auto &entry : weights) {
        const double count = image_data[entry.first];
        background_variance += entry.second * entry.second * count;
        const auto y = entry.first / width;
        const auto x = entry.first % width;
        if (x >= x0 && x < x1 && y >= y0 && y < y1) {
            covariance += entry.second * count;
        }
    }
    const double ratio = background_pixels > 0.0
        ? center_pixels / background_pixels : 0.0;
    moments[0] = original_center_count;
    moments[1] = background_variance;
    moments[2] = covariance;
    moments[3] = original_center_count + ratio * ratio * background_variance
        - 2.0 * ratio * covariance;
}

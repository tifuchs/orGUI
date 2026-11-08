// Immutable geometry for repeated frames with a fixed detector mask and ROIs.
// All per-frame scratch is local, so frame workers may share one plan.
class PolynomialBackgroundPlan {
    struct Rectangle {
        std::int64_t x0, x1, y0, y1;
        bool contains(std::int64_t x, std::int64_t y) const {
            return x >= x0 && x < x1 && y >= y0 && y < y1;
        }
        bool overlaps(const Rectangle &other) const {
            return std::max(x0, other.x0) < std::min(x1, other.x1)
                && std::max(y0, other.y0) < std::min(y1, other.y1);
        }
    };

    struct Row {
        double influence[15] = {0.0};
        Rectangle center;
        Rectangle background[4];
        double x_origin = 0.0;
        double y_origin = 0.0;
        double inverse_scale = 1.0;
        double samples = 0.0;
        double center_pixels = 0.0;
        int order = -1;
        bool overlaps_center = false;
        bool overlapping_strips = false;
    };

    struct BackgroundSums {
        double image = 0.0;
        double variance = 0.0;
        double covariance = 0.0;
        double correction = 0.0;
        double pixels = 0.0;
        bool missing = false;
    };

public:
    PolynomialBackgroundPlan(
        const BoolArray2D mask,
        const RoiArray3D center,
        const RoiArray3D left,
        const RoiArray3D right,
        const RoiArray3D top,
        const RoiArray3D bottom,
        const int fit_order,
        const PixelRepairPlan* repair_plan = nullptr
    ) : order_(fit_order), repair_plan_(repair_plan) {
        if (order_ < 0 || order_ > 4) {
            throw py::value_error("fit_order must be between 0 and 4");
        }
        const auto mask_info = mask.request();
        require_ndim(mask_info, 2, "mask");
        height_ = mask_info.shape[0];
        width_ = mask_info.shape[1];
        mask_ = BoolArray2D({height_, width_});
        std::copy(bool_ptr(mask_info), bool_ptr(mask_info) + mask_info.size,
                  mask_.mutable_data());
        if (repair_plan_) {
            repair_plan_->require_shape(height_, width_);
            repaired_ = BoolArray2D({height_, width_});
            repair_plan_->repair_footprint(repaired_.mutable_data());
        }
        const RoiArray3D input[5] = {center, left, right, top, bottom};
        for (int j = 0; j < 5; ++j) {
            const auto info = input[j].request();
            require_roi_shape(info, "ROI");
            if (j && info.shape[0] != rois_[0].shape(0)) {
                throw py::value_error("all ROI arrays must have the same length");
            }
            rois_[j] = RoiArray3D(info.shape);
            std::copy(roi_ptr(info), roi_ptr(info) + info.size,
                      rois_[j].mutable_data());
            for (py::ssize_t n = 0; n < info.shape[0]; ++n) {
                for (int axis = 0; axis < 2; ++axis) {
                    const auto low = roi_value(roi_ptr(info), info, n, axis, 0);
                    const auto high = roi_value(roi_ptr(info), info, n, axis, 1);
                    const auto limit = axis == 0 ? width_ : height_;
                    if (low < 0 || high < low || high > limit) {
                        throw py::value_error("ROI bounds must lie inside the image");
                    }
                }
            }
        }
        const auto center_info = rois_[0].request();
        const auto left_info = rois_[1].request();
        const auto right_info = rois_[2].request();
        const auto top_info = rois_[3].request();
        const auto bottom_info = rois_[4].request();
        const py::buffer_info *background[4] = {
            &left_info, &right_info, &top_info, &bottom_info,
        };
        rows_.reserve(center_info.shape[0]);
        py::gil_scoped_release release;
        for (py::ssize_t n = 0; n < center_info.shape[0]; ++n) {
            rows_.push_back(prepare(nullptr, center_info, background, n,
                repair_plan_ ? repaired_.data() : nullptr));
        }
    }

    void process(
        const Array2D image,
        const Array2D correction,
        Array2D counters,
        Array2D correction_counters,
        py::object moments = py::none()
    ) const {
        const auto image_info = image.request();
        const auto correction_info = correction.request();
        const auto counter_info = counters.request();
        const auto correction_counter_info = correction_counters.request();
        const auto mask_info = mask_.request();
        require_same_image_shape(mask_info, image_info, "image");
        require_same_image_shape(mask_info, correction_info, "C_corr");
        require_counter_shape(counter_info, rows_.size(), "all_counters");
        require_counter_shape(correction_counter_info, rows_.size(), "all_Carr");
        Array2D moment_array;
        double *moment_data = nullptr;
        if (!moments.is_none()) {
            moment_array = moments.cast<Array2D>();
            const auto info = moment_array.request();
            require_counter_shape(info, rows_.size(), "moments");
            moment_data = mutable_double_ptr(info);
        }
        const auto center_info = rois_[0].request();
        const auto left_info = rois_[1].request();
        const auto right_info = rois_[2].request();
        const auto top_info = rois_[3].request();
        const auto bottom_info = rois_[4].request();
        const py::buffer_info *background[4] = {
            &left_info, &right_info, &top_info, &bottom_info,
        };
        const double *data = double_ptr(image_info);
        const double *correction_data = double_ptr(correction_info);
        const bool *mask_data = bool_ptr(mask_info);
        double *output = mutable_double_ptr(counter_info);
        double *corr_output = mutable_double_ptr(correction_counter_info);
        BoolArray2D repaired;
        const bool* repaired_data = nullptr;
        RoiRepairFrame repair_frame;
        if (repair_plan_) {
            repaired = BoolArray2D({height_, width_});
            std::fill(repaired.mutable_data(), repaired.mutable_data() + repaired.size(), false);
            repair_frame = repair_plan_->apply_roi_records(image, repaired, moment_data != nullptr);
            repaired_data = repaired.data();
        }
        py::gil_scoped_release release;
        for (py::ssize_t n = 0; n < static_cast<py::ssize_t>(rows_.size()); ++n) {
            const Row *row = &rows_[n];
            Row fallback;
            const RoiSums center = repaired_data ? sum_roi_with_mask(
                data, correction_data, nullptr, mask_data, repaired_data, true,
                width_, roi_ptr(center_info), center_info, n) : sum_roi(
                data, correction_data, nullptr, mask_data, width_,
                roi_ptr(center_info), center_info, n);
            if (repaired_data && center.pixels != row->center_pixels) {
                fallback = prepare(data, center_info, background, n, repaired_data);
                row = &fallback;
            }
            // Missing measured donors must be excluded, not replaced by zeros.
            // Rebuild only this ROI when its finite-source pattern changes.
            BackgroundSums sums = stream_background(
                *row, data, correction_data, mask_data, moment_data != nullptr);
            if (sums.missing) {
                fallback = prepare(data, center_info, background, n, repaired_data);
                row = &fallback;
                sums = stream_background(
                    *row, data, correction_data, mask_data, moment_data != nullptr);
            }
            output[4 * n] = center.image;
            output[4 * n + 1] = center.pixels;
            output[4 * n + 2] = sums.image;
            output[4 * n + 3] = row->samples;
            corr_output[4 * n] = center.correction;
            corr_output[4 * n + 1] = center.pixels;
            corr_output[4 * n + 2] = sums.correction;
            corr_output[4 * n + 3] = sums.pixels;
            if (moment_data) {
                const double ratio = row->samples > 0.0
                    ? center.pixels / row->samples : 0.0;
                moment_data[4 * n] = center.image;
                moment_data[4 * n + 1] = sums.variance;
                moment_data[4 * n + 2] = sums.covariance;
                moment_data[4 * n + 3] = center.image
                    + ratio * ratio * sums.variance - 2.0 * ratio * sums.covariance;
                if (repair_plan_) {
                    add_repair_covariance(repair_frame, repair_rectangle(center_info, n),
                        data, width_, repair_plan_->maximum_span(), center.pixels,
                        row->samples, [&](py::ssize_t index) {
                            const auto x = index % width_, y = index / width_;
                            int multiplicity = 0;
                            for (const auto& rectangle : row->background)
                                multiplicity += rectangle.contains(x, y);
                            if (!multiplicity || row->order < 0) return 0.0;
                            double basis[15];
                            polynomial_basis((x - row->x_origin) * row->inverse_scale,
                                (y - row->y_origin) * row->inverse_scale, row->order, basis);
                            double value = 0.0;
                            for (int j = 0; j < polynomial_terms(row->order); ++j)
                                value += basis[j] * row->influence[j];
                            return multiplicity * value;
                        }, moment_data + 4 * n, true);
                }
            }
        }
    }

    py::dict configuration() const {
        py::dict result;
        result["shape"] = py::make_tuple(height_, width_);
        result["fit_order"] = order_;
        result["rois"] = rows_.size();
        std::size_t storage = mask_.size() * sizeof(bool);
        for (const auto &roi : rois_) {
            storage += roi.size() * sizeof(std::int64_t);
        }
        storage += rows_.capacity() * sizeof(Row);
        if (repair_plan_) storage += repaired_.size() * sizeof(bool);
        result["storage_bytes"] = storage;
        return result;
    }

private:
    BackgroundSums stream_background(
        const Row &row, const double *data, const double *correction,
        const bool *mask, const bool moments
    ) const {
        BackgroundSums sums;
        const int degree = std::max(0, row.order);
        for (const auto &rectangle : row.background) {
            for (auto y = rectangle.y0; y < rectangle.y1; ++y) {
                const double yn = (y - row.y_origin) * row.inverse_scale;
                double powers[5] = {1.0};
                double coefficients[5] = {0.0};
                for (int j = 1; j <= degree; ++j) {
                    powers[j] = powers[j - 1] * yn;
                }
                for (int d = 0; d <= degree; ++d) {
                    for (int j = 0; j <= d; ++j) {
                        coefficients[d - j] +=
                            row.influence[d * (d + 1) / 2 + j] * powers[j];
                    }
                }
                double xn = (rectangle.x0 - row.x_origin) * row.inverse_scale;
                for (auto x = rectangle.x0; x < rectangle.x1;
                     ++x, xn += row.inverse_scale) {
                    const auto index = y * width_ + x;
                    if (!std::isnan(correction[index])) {
                        sums.correction += correction[index];
                    }
                    if (mask[index]) {
                        continue;
                    }
                    sums.pixels += 1.0;
                    const double value = data[index];
                    if (std::isnan(value)) {
                        sums.missing = true;
                        continue;
                    }
                    double weight = coefficients[degree];
                    for (int j = degree - 1; j >= 0; --j) {
                        weight = weight * xn + coefficients[j];
                    }
                    sums.image += weight * value;
                    if (moments) {
                        int multiplicity = 1;
                        if (row.overlapping_strips) {
                            multiplicity = 0;
                            for (const auto &other : row.background) {
                                multiplicity += other.contains(x, y);
                            }
                        }
                        // A repeated pixel appears m times in the stream;
                        // m*w^2 per occurrence gives (m*w)^2, not m*w^2.
                        sums.variance += multiplicity * weight * weight * value;
                        if (row.overlaps_center && row.center.contains(x, y)) {
                            sums.covariance += weight * value;
                        }
                    }
                }
            }
        }
        return sums;
    }

    Row prepare(const double *data, const py::buffer_info &center,
                const py::buffer_info *const background[4],
                const py::ssize_t n, const bool* repaired = nullptr) const {
        Row row;
        const auto x0 = roi_value(roi_ptr(center), center, n, 0, 0);
        const auto x1 = roi_value(roi_ptr(center), center, n, 0, 1);
        const auto y0 = roi_value(roi_ptr(center), center, n, 1, 0);
        const auto y1 = roi_value(roi_ptr(center), center, n, 1, 1);
        const double x_origin = 0.5 * (x0 + x1 - 1.0);
        const double y_origin = 0.5 * (y0 + y1 - 1.0);
        const double scale = polynomial_coordinate_scale(
            center, background, n, x_origin, y_origin, order_);
        row.center = {x0, x1, y0, y1};
        row.x_origin = x_origin;
        row.y_origin = y_origin;
        row.inverse_scale = 1.0 / scale;
        for (int j = 0; j < 4; ++j) {
            const auto info = background[j];
            row.background[j] = {
                roi_value(roi_ptr(*info), *info, n, 0, 0),
                roi_value(roi_ptr(*info), *info, n, 0, 1),
                roi_value(roi_ptr(*info), *info, n, 1, 0),
                roi_value(roi_ptr(*info), *info, n, 1, 1),
            };
            row.overlaps_center = row.overlaps_center
                || row.center.overlaps(row.background[j]);
            for (int k = 0; k < j; ++k) {
                row.overlapping_strips = row.overlapping_strips
                    || row.background[k].overlaps(row.background[j]);
            }
        }
        const bool *mask = mask_.data();
        for (auto y = y0; y < y1; ++y) {
            for (auto x = x0; x < x1; ++x) {
                const auto index = y * width_ + x;
                row.center_pixels += !mask[index] || (repaired && repaired[index]);
            }
        }
        const auto fit = fit_local_polynomial_background(
            data, mask, width_, roi_ptr(*background[0]), *background[0],
            roi_ptr(*background[1]), *background[1],
            roi_ptr(*background[2]), *background[2],
            roi_ptr(*background[3]), *background[3], n,
            x_origin, y_origin, scale, order_);
        row.samples = fit.fitted ? fit.sample_count : 0.0;
        row.order = fit.order;
        polynomial_prediction_coefficients(
            data, mask, width_, center, background, n,
            x_origin, y_origin, scale, fit, row.center_pixels, row.influence, repaired);
        return row;
    }

    BoolArray2D mask_;
    RoiArray3D rois_[5];
    py::ssize_t height_;
    py::ssize_t width_;
    int order_;
    std::vector<Row> rows_;
    const PixelRepairPlan* repair_plan_;
    BoolArray2D repaired_;
};

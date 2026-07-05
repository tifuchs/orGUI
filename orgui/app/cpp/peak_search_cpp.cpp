#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using Array2D = py::array_t<double, py::array::c_style | py::array::forcecast>;
using BoolArray2D = py::array_t<bool, py::array::c_style | py::array::forcecast>;

struct MaximumResult {
    bool valid = false;
    double value = std::numeric_limits<double>::quiet_NaN();
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
};

struct FrameMaximum {
    bool available = false;
    bool valid = false;
    double value = std::numeric_limits<double>::quiet_NaN();
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    double sharpness = std::numeric_limits<double>::quiet_NaN();
    double derivative_sharpness = std::numeric_limits<double>::quiet_NaN();
    double prominence = std::numeric_limits<double>::quiet_NaN();
};

struct Candidate {
    int index = -1;
    double value = std::numeric_limits<double>::quiet_NaN();
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    double sharpness = std::numeric_limits<double>::quiet_NaN();
    double derivative_sharpness = std::numeric_limits<double>::quiet_NaN();
    double prominence = std::numeric_limits<double>::quiet_NaN();
};

void require_ndim(const py::buffer_info &info, const int ndim, const char *name) {
    if (info.ndim != ndim) {
        throw py::value_error(std::string(name) + " has unexpected number of dimensions");
    }
}

void require_same_shape(
    const py::buffer_info &image,
    const py::buffer_info &other,
    const char *name
) {
    require_ndim(other, 2, name);
    if (image.shape[0] != other.shape[0] || image.shape[1] != other.shape[1]) {
        throw py::value_error(std::string(name) + " must have the same shape as image");
    }
}

inline py::ssize_t index2(
    const py::buffer_info &info,
    const py::ssize_t row,
    const py::ssize_t col
) {
    return row * info.shape[1] + col;
}

std::vector<double> finite_values(const std::vector<double> &values) {
    std::vector<double> finite;
    finite.reserve(values.size());
    for (const double value : values) {
        if (std::isfinite(value)) {
            finite.push_back(value);
        }
    }
    return finite;
}

double median(std::vector<double> values) {
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const auto mid = values.begin() + static_cast<std::ptrdiff_t>(values.size() / 2);
    std::nth_element(values.begin(), mid, values.end());
    const double upper = *mid;
    if (values.size() % 2 == 1) {
        return upper;
    }
    const auto lower_it = std::max_element(values.begin(), mid);
    return (*lower_it + upper) / 2.0;
}

double stddev(const std::vector<double> &values) {
    const std::vector<double> finite = finite_values(values);
    if (finite.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double mean = 0.0;
    for (const double value : finite) {
        mean += value;
    }
    mean /= static_cast<double>(finite.size());
    double variance = 0.0;
    for (const double value : finite) {
        const double diff = value - mean;
        variance += diff * diff;
    }
    variance /= static_cast<double>(finite.size());
    return std::sqrt(variance);
}

std::pair<double, double> robust_location_scale(const std::vector<double> &values) {
    std::vector<double> finite = finite_values(values);
    if (finite.empty()) {
        return {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        };
    }
    const double location = median(finite);
    std::vector<double> deviations;
    deviations.reserve(finite.size());
    for (const double value : finite) {
        deviations.push_back(std::abs(value - location));
    }
    return {location, 1.4826 * median(deviations)};
}

double robust_z(const double value, const std::vector<double> &values) {
    const auto [location, raw_scale] = robust_location_scale(values);
    if (!std::isfinite(location)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double scale = raw_scale;
    const double eps = std::numeric_limits<double>::epsilon();
    if (!std::isfinite(scale) || scale <= eps) {
        scale = stddev(values);
    }
    if (!std::isfinite(scale) || scale <= eps) {
        scale = eps;
    }
    return (value - location) / scale;
}

MaximumResult masked_maximum(
    const Array2D image,
    const py::object mask_obj = py::none(),
    const py::object background_obj = py::none(),
    const int mask_distance = 0
) {
    const auto image_info = image.request();
    require_ndim(image_info, 2, "image");
    py::buffer_info mask_info;
    py::buffer_info background_info;
    const bool has_mask = !mask_obj.is_none();
    const bool has_background = !background_obj.is_none();
    BoolArray2D mask_array;
    Array2D background_array;
    if (has_mask) {
        mask_array = py::cast<BoolArray2D>(mask_obj);
        mask_info = mask_array.request();
        require_same_shape(image_info, mask_info, "mask");
    }
    if (has_background) {
        background_array = py::cast<Array2D>(background_obj);
        background_info = background_array.request();
        require_same_shape(image_info, background_info, "background");
    }

    const double *image_data = static_cast<const double *>(image_info.ptr);
    const bool *mask_data = has_mask ? static_cast<const bool *>(mask_info.ptr) : nullptr;
    const double *background_data = has_background
        ? static_cast<const double *>(background_info.ptr)
        : nullptr;

    MaximumResult result;
    {
        py::gil_scoped_release release;
        for (py::ssize_t row = 0; row < image_info.shape[0]; ++row) {
            for (py::ssize_t col = 0; col < image_info.shape[1]; ++col) {
                const auto idx = index2(image_info, row, col);
                if (mask_data != nullptr && mask_data[idx]) {
                    continue;
                }
                if (mask_data != nullptr && mask_distance > 0) {
                    const py::ssize_t y0 = std::max<py::ssize_t>(0, row - mask_distance);
                    const py::ssize_t y1 = std::min<py::ssize_t>(
                        image_info.shape[0], row + mask_distance + 1
                    );
                    const py::ssize_t x0 = std::max<py::ssize_t>(0, col - mask_distance);
                    const py::ssize_t x1 = std::min<py::ssize_t>(
                        image_info.shape[1], col + mask_distance + 1
                    );
                    bool near_mask = false;
                    for (py::ssize_t yy = y0; yy < y1 && !near_mask; ++yy) {
                        for (py::ssize_t xx = x0; xx < x1; ++xx) {
                            if (mask_data[index2(mask_info, yy, xx)]) {
                                near_mask = true;
                                break;
                            }
                        }
                    }
                    if (near_mask) {
                        continue;
                    }
                }
                double value = image_data[idx];
                if (background_data != nullptr) {
                    value -= background_data[idx];
                }
                if (!std::isfinite(value)) {
                    continue;
                }
                if (!result.valid || value > result.value) {
                    result.valid = true;
                    result.value = value;
                    result.x = static_cast<double>(col);
                    result.y = static_cast<double>(row);
                }
            }
        }
    }
    return result;
}

double masked_roi_sum(
    const Array2D image,
    const py::object mask_obj,
    const py::object background_obj,
    const int y0,
    const int y1,
    const int x0,
    const int x1
) {
    const auto image_info = image.request();
    require_ndim(image_info, 2, "image");
    py::buffer_info mask_info;
    py::buffer_info background_info;
    const bool has_mask = !mask_obj.is_none();
    const bool has_background = !background_obj.is_none();
    BoolArray2D mask_array;
    Array2D background_array;
    if (has_mask) {
        mask_array = py::cast<BoolArray2D>(mask_obj);
        mask_info = mask_array.request();
        require_same_shape(image_info, mask_info, "mask");
    }
    if (has_background) {
        background_array = py::cast<Array2D>(background_obj);
        background_info = background_array.request();
        require_same_shape(image_info, background_info, "background");
    }
    const auto row0 = std::max<py::ssize_t>(0, y0);
    const auto row1 = std::min<py::ssize_t>(image_info.shape[0], y1);
    const auto col0 = std::max<py::ssize_t>(0, x0);
    const auto col1 = std::min<py::ssize_t>(image_info.shape[1], x1);
    const double *image_data = static_cast<const double *>(image_info.ptr);
    const bool *mask_data = has_mask ? static_cast<const bool *>(mask_info.ptr) : nullptr;
    const double *background_data = has_background
        ? static_cast<const double *>(background_info.ptr)
        : nullptr;

    double total = 0.0;
    {
        py::gil_scoped_release release;
        for (py::ssize_t row = row0; row < row1; ++row) {
            for (py::ssize_t col = col0; col < col1; ++col) {
                const auto idx = index2(image_info, row, col);
                if (mask_data != nullptr && mask_data[idx]) {
                    continue;
                }
                double value = image_data[idx];
                if (background_data != nullptr) {
                    value -= background_data[idx];
                }
                if (std::isfinite(value)) {
                    total += value;
                }
            }
        }
    }
    return total;
}

class PeakCandidateDetector {
public:
    PeakCandidateDetector(
        const int n_frames,
        const int burn_in,
        const int history,
        const int min_history,
        const double level_z,
        const double derivative_z,
        const int lookahead,
        const double min_prominence_z,
        const int refractory
    )
        : n_frames_(n_frames),
          burn_in_(burn_in),
          history_(history),
          min_history_(min_history),
          level_z_(level_z),
          derivative_z_(derivative_z),
          lookahead_(lookahead),
          min_prominence_z_(min_prominence_z),
          refractory_(refractory),
          next_allowed_(burn_in),
          frames_(static_cast<std::size_t>(n_frames))
    {
        if (n_frames < 0) {
            throw py::value_error("n_frames must be non-negative");
        }
    }

    std::vector<Candidate> push(
        const int index,
        const bool valid,
        const double value,
        const double x,
        const double y,
        const bool finish = false
    ) {
        if (index < 0 || index >= n_frames_) {
            throw py::index_error("frame index out of range");
        }
        FrameMaximum &frame = frames_[static_cast<std::size_t>(index)];
        frame.available = true;
        frame.valid = valid;
        frame.value = valid ? value : std::numeric_limits<double>::quiet_NaN();
        frame.x = valid ? x : std::numeric_limits<double>::quiet_NaN();
        frame.y = valid ? y : std::numeric_limits<double>::quiet_NaN();
        if (index > highest_seen_) {
            highest_seen_ = index;
        }
        return drain(finish);
    }

    std::vector<Candidate> finish() {
        return drain(true);
    }

private:
    bool ready_to_evaluate(const int index, const bool finish) const {
        if (index >= n_frames_) {
            return false;
        }
        if (!frames_[static_cast<std::size_t>(index)].available) {
            return false;
        }
        if (finish) {
            return true;
        }
        const int required = std::min(n_frames_ - 1, index + lookahead_);
        return highest_seen_ >= required;
    }

    std::vector<double> baseline_values(const int index) const {
        std::vector<double> baseline;
        const int start = std::max(0, index - history_);
        for (int i = start; i < index; ++i) {
            const FrameMaximum &frame = frames_[static_cast<std::size_t>(i)];
            if (frame.available && frame.valid && std::isfinite(frame.value)) {
                baseline.push_back(frame.value);
            }
        }
        return baseline;
    }

    std::vector<Candidate> drain(const bool finish) {
        std::vector<Candidate> candidates;
        py::gil_scoped_release release;
        while (ready_to_evaluate(eval_index_, finish)) {
            FrameMaximum &maximum = frames_[static_cast<std::size_t>(eval_index_)];
            if (!maximum.valid || !std::isfinite(maximum.value)) {
                ++eval_index_;
                continue;
            }

            const std::vector<double> baseline = baseline_values(eval_index_);
            if (eval_index_ < next_allowed_
                || static_cast<int>(baseline.size()) < min_history_) {
                ++eval_index_;
                continue;
            }

            std::vector<double> diffs;
            if (baseline.size() > 1) {
                diffs.reserve(baseline.size() - 1);
                for (std::size_t i = 1; i < baseline.size(); ++i) {
                    diffs.push_back(baseline[i] - baseline[i - 1]);
                }
            }
            const double previous = baseline.back();
            const double current_diff = maximum.value - previous;
            const double level_score = robust_z(maximum.value, baseline);
            const double diff_score = robust_z(current_diff, diffs);
            maximum.sharpness = level_score;
            maximum.derivative_sharpness = diff_score;

            if (!std::isfinite(level_score)
                || !std::isfinite(diff_score)
                || level_score < level_z_
                || diff_score < derivative_z_) {
                ++eval_index_;
                continue;
            }

            int candidate_index = eval_index_;
            double candidate_value = maximum.value;
            for (
                int look_index = eval_index_ + 1;
                look_index < std::min(n_frames_, eval_index_ + lookahead_ + 1);
                ++look_index
            ) {
                const FrameMaximum &look = frames_[static_cast<std::size_t>(look_index)];
                if (look.available && look.valid && look.value > candidate_value) {
                    candidate_index = look_index;
                    candidate_value = look.value;
                }
            }

            FrameMaximum &candidate_frame = frames_[static_cast<std::size_t>(candidate_index)];
            const auto [location, raw_scale] = robust_location_scale(baseline);
            double scale = raw_scale;
            const double eps = std::numeric_limits<double>::epsilon();
            if (!std::isfinite(scale) || scale <= eps) {
                scale = stddev(baseline);
            }
            if (!std::isfinite(scale) || scale <= eps) {
                scale = eps;
            }
            const double prominence = candidate_frame.value - location;
            candidate_frame.prominence = prominence;

            if (prominence / scale >= min_prominence_z_) {
                candidate_frame.sharpness = robust_z(candidate_frame.value, baseline);
                const std::vector<double> candidate_baseline = baseline_values(candidate_index);
                double prev_value = std::numeric_limits<double>::quiet_NaN();
                if (!candidate_baseline.empty()) {
                    prev_value = candidate_baseline.back();
                }
                candidate_frame.derivative_sharpness = robust_z(
                    candidate_frame.value - prev_value,
                    diffs
                );
                candidates.push_back(Candidate{
                    candidate_index,
                    candidate_frame.value,
                    candidate_frame.x,
                    candidate_frame.y,
                    candidate_frame.sharpness,
                    candidate_frame.derivative_sharpness,
                    candidate_frame.prominence
                });
                next_allowed_ = candidate_index + refractory_ + 1;
            }
            eval_index_ = std::max(eval_index_ + 1, candidate_index + 1);
        }
        return candidates;
    }

    int n_frames_;
    int burn_in_;
    int history_;
    int min_history_;
    double level_z_;
    double derivative_z_;
    int lookahead_;
    double min_prominence_z_;
    int refractory_;
    int next_allowed_;
    int eval_index_ = 0;
    int highest_seen_ = -1;
    std::vector<FrameMaximum> frames_;
};

PYBIND11_MODULE(_peak_search_cpp, module) {
    module.doc() = "C++ acceleration kernels for auto-Bragg peak search.";

    py::class_<MaximumResult>(module, "MaximumResult")
        .def_readonly("valid", &MaximumResult::valid)
        .def_readonly("value", &MaximumResult::value)
        .def_readonly("x", &MaximumResult::x)
        .def_readonly("y", &MaximumResult::y);

    py::class_<Candidate>(module, "Candidate")
        .def_readonly("index", &Candidate::index)
        .def_readonly("value", &Candidate::value)
        .def_readonly("x", &Candidate::x)
        .def_readonly("y", &Candidate::y)
        .def_readonly("sharpness", &Candidate::sharpness)
        .def_readonly("derivative_sharpness", &Candidate::derivative_sharpness)
        .def_readonly("prominence", &Candidate::prominence);

    py::class_<PeakCandidateDetector>(module, "PeakCandidateDetector")
        .def(
            py::init<int, int, int, int, double, double, int, double, int>(),
            py::arg("n_frames"),
            py::arg("burn_in"),
            py::arg("history"),
            py::arg("min_history"),
            py::arg("level_z"),
            py::arg("derivative_z"),
            py::arg("lookahead"),
            py::arg("min_prominence_z"),
            py::arg("refractory")
        )
        .def(
            "push",
            &PeakCandidateDetector::push,
            py::arg("index"),
            py::arg("valid"),
            py::arg("value"),
            py::arg("x"),
            py::arg("y"),
            py::arg("finish") = false
        )
        .def("finish", &PeakCandidateDetector::finish);

    module.def(
        "masked_maximum",
        &masked_maximum,
        py::arg("image"),
        py::arg("mask") = py::none(),
        py::arg("background") = py::none(),
        py::arg("mask_distance") = 0
    );
    module.def(
        "masked_roi_sum",
        &masked_roi_sum,
        py::arg("image"),
        py::arg("mask") = py::none(),
        py::arg("background") = py::none(),
        py::arg("y0"),
        py::arg("y1"),
        py::arg("x0"),
        py::arg("x1")
    );
}

// Copyright (c) 2026 Timo Fuchs
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <iterator>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <xxhash.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using FloatArray = py::array_t<double, py::array::c_style | py::array::forcecast>;
using BoolArray = py::array_t<bool, py::array::c_style | py::array::forcecast>;
using Int64Array = py::array_t<std::int64_t>;
using UInt32Array = py::array_t<std::uint32_t>;
using ContiguousUInt32Array = py::array_t<
    std::uint32_t,
    py::array::c_style | py::array::forcecast
>;

std::string xxh128_hex(const XXH128_hash_t hash) {
    XXH128_canonical_t canonical;
    XXH128_canonicalFromHash(&canonical, hash);
    constexpr char digits[] = "0123456789abcdef";
    std::string result(32, '0');
    for (std::size_t index = 0; index < sizeof(canonical.digest); ++index) {
        const unsigned char value = canonical.digest[index];
        result[2 * index] = digits[value >> 4];
        result[2 * index + 1] = digits[value & 0x0f];
    }
    return result;
}

py::buffer_info contiguous_buffer(const py::buffer &buffer) {
    py::buffer_info info = buffer.request();
    py::ssize_t expected_stride = info.itemsize;
    for (py::ssize_t axis = info.ndim - 1; axis >= 0; --axis) {
        if (info.shape[axis] > 1 && info.strides[axis] != expected_stride) {
            throw py::value_error("XXH3 input buffer must be C-contiguous");
        }
        expected_stride *= info.shape[axis];
    }
    return info;
}

std::string xxh3_128_buffer(const py::buffer &buffer) {
    const py::buffer_info info = contiguous_buffer(buffer);
    const std::size_t size = static_cast<std::size_t>(
        info.size * info.itemsize
    );
    XXH128_hash_t hash;
    {
        py::gil_scoped_release release;
        hash = XXH3_128bits(info.ptr, size);
    }
    return xxh128_hex(hash);
}

py::dict merge_sorted_batches(
    const ContiguousUInt32Array &left_chunk,
    const ContiguousUInt32Array &left_local,
    const FloatArray &left_intensity,
    const FloatArray &left_variance,
    const FloatArray &left_weight,
    const ContiguousUInt32Array &left_contributors,
    const ContiguousUInt32Array &right_chunk,
    const ContiguousUInt32Array &right_local,
    const FloatArray &right_intensity,
    const FloatArray &right_variance,
    const FloatArray &right_weight,
    const ContiguousUInt32Array &right_contributors
) {
    const auto lc = left_chunk.request();
    const auto ll = left_local.request();
    const auto li = left_intensity.request();
    const auto lv = left_variance.request();
    const auto lw = left_weight.request();
    const auto ln = left_contributors.request();
    const auto rc = right_chunk.request();
    const auto rl = right_local.request();
    const auto ri = right_intensity.request();
    const auto rv = right_variance.request();
    const auto rw = right_weight.request();
    const auto rn = right_contributors.request();
    const auto validate = [](const py::buffer_info &info, const char *name) {
        if (info.ndim != 1) {
            throw py::value_error(std::string(name) + " must be one-dimensional");
        }
    };
    validate(lc, "left_chunk_id");
    validate(ll, "left_local_voxel_id");
    validate(li, "left_weighted_intensity");
    validate(lv, "left_weighted_variance");
    validate(lw, "left_weight");
    validate(ln, "left_contributors");
    validate(rc, "right_chunk_id");
    validate(rl, "right_local_voxel_id");
    validate(ri, "right_weighted_intensity");
    validate(rv, "right_weighted_variance");
    validate(rw, "right_weight");
    validate(rn, "right_contributors");
    const py::ssize_t left_size = lc.size;
    const py::ssize_t right_size = rc.size;
    if (
        ll.size != left_size || li.size != left_size || lv.size != left_size
        || lw.size != left_size || ln.size != left_size
    ) {
        throw py::value_error("Left batch columns must have equal lengths");
    }
    if (
        rl.size != right_size || ri.size != right_size || rv.size != right_size
        || rw.size != right_size || rn.size != right_size
    ) {
        throw py::value_error("Right batch columns must have equal lengths");
    }
    const auto *lc_data = static_cast<const std::uint32_t *>(lc.ptr);
    const auto *ll_data = static_cast<const std::uint32_t *>(ll.ptr);
    const auto *li_data = static_cast<const double *>(li.ptr);
    const auto *lv_data = static_cast<const double *>(lv.ptr);
    const auto *lw_data = static_cast<const double *>(lw.ptr);
    const auto *ln_data = static_cast<const std::uint32_t *>(ln.ptr);
    const auto *rc_data = static_cast<const std::uint32_t *>(rc.ptr);
    const auto *rl_data = static_cast<const std::uint32_t *>(rl.ptr);
    const auto *ri_data = static_cast<const double *>(ri.ptr);
    const auto *rv_data = static_cast<const double *>(rv.ptr);
    const auto *rw_data = static_cast<const double *>(rw.ptr);
    const auto *rn_data = static_cast<const std::uint32_t *>(rn.ptr);
    const auto less = [](
        const std::uint32_t chunk_a,
        const std::uint32_t local_a,
        const std::uint32_t chunk_b,
        const std::uint32_t local_b
    ) {
        return chunk_a < chunk_b
            || (chunk_a == chunk_b && local_a < local_b);
    };
    py::ssize_t output_size = 0;
    {
        py::gil_scoped_release release;
        py::ssize_t left = 0;
        py::ssize_t right = 0;
        while (left < left_size && right < right_size) {
            if (
                lc_data[left] == rc_data[right]
                && ll_data[left] == rl_data[right]
            ) {
                ++left;
                ++right;
            } else if (less(
                lc_data[left],
                ll_data[left],
                rc_data[right],
                rl_data[right]
            )) {
                ++left;
            } else {
                ++right;
            }
            ++output_size;
        }
        output_size += left_size - left + right_size - right;
    }
    ContiguousUInt32Array output_chunk(output_size);
    ContiguousUInt32Array output_local(output_size);
    FloatArray output_intensity(output_size);
    FloatArray output_variance(output_size);
    FloatArray output_weight(output_size);
    ContiguousUInt32Array output_contributors(output_size);
    auto *oc_data = static_cast<std::uint32_t *>(output_chunk.request().ptr);
    auto *ol_data = static_cast<std::uint32_t *>(output_local.request().ptr);
    auto *oi_data = static_cast<double *>(output_intensity.request().ptr);
    auto *ov_data = static_cast<double *>(output_variance.request().ptr);
    auto *ow_data = static_cast<double *>(output_weight.request().ptr);
    auto *on_data = static_cast<std::uint32_t *>(
        output_contributors.request().ptr
    );
    {
        py::gil_scoped_release release;
        py::ssize_t left = 0;
        py::ssize_t right = 0;
        py::ssize_t output = 0;
        while (left < left_size || right < right_size) {
            const bool take_left = right >= right_size || (
                left < left_size && less(
                    lc_data[left],
                    ll_data[left],
                    rc_data[right],
                    rl_data[right]
                )
            );
            const bool equal = (
                left < left_size && right < right_size
                && lc_data[left] == rc_data[right]
                && ll_data[left] == rl_data[right]
            );
            if (equal) {
                oc_data[output] = lc_data[left];
                ol_data[output] = ll_data[left];
                oi_data[output] = li_data[left] + ri_data[right];
                ov_data[output] = lv_data[left] + rv_data[right];
                ow_data[output] = lw_data[left] + rw_data[right];
                on_data[output] = ln_data[left] + rn_data[right];
                ++left;
                ++right;
            } else if (take_left) {
                oc_data[output] = lc_data[left];
                ol_data[output] = ll_data[left];
                oi_data[output] = li_data[left];
                ov_data[output] = lv_data[left];
                ow_data[output] = lw_data[left];
                on_data[output] = ln_data[left];
                ++left;
            } else {
                oc_data[output] = rc_data[right];
                ol_data[output] = rl_data[right];
                oi_data[output] = ri_data[right];
                ov_data[output] = rv_data[right];
                ow_data[output] = rw_data[right];
                on_data[output] = rn_data[right];
                ++right;
            }
            ++output;
        }
    }
    py::dict result;
    result["chunk_id"] = std::move(output_chunk);
    result["local_voxel_id"] = std::move(output_local);
    result["weighted_intensity"] = std::move(output_intensity);
    result["weighted_variance"] = std::move(output_variance);
    result["weight"] = std::move(output_weight);
    result["contributors"] = std::move(output_contributors);
    return result;
}

struct Vec3 {
    double x;
    double y;
    double z;
};

struct Mat3 {
    std::array<double, 9> value{};
};

struct PixelRays {
    Vec3 base;
    Vec3 du;
    Vec3 dv;
    Vec3 duv;
};

struct FrameRotation {
    double sin_alpha;
    double cos_alpha;
    double sin_omega;
    double cos_omega;
    double sin_chi;
    double cos_chi;
    double sin_phi;
    double cos_phi;
};

struct CoordinateTransform {
    Mat3 matrix;
    Vec3 offset;
};

Vec3 apply(const Mat3 &matrix, const Vec3 &vector) {
    return {
        matrix.value[0] * vector.x + matrix.value[1] * vector.y
            + matrix.value[2] * vector.z,
        matrix.value[3] * vector.x + matrix.value[4] * vector.y
            + matrix.value[5] * vector.z,
        matrix.value[6] * vector.x + matrix.value[7] * vector.y
            + matrix.value[8] * vector.z,
    };
}

Vec3 rotate_x(const Vec3 &vector, const double angle) {
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    return {
        vector.x,
        cosine * vector.y - sine * vector.z,
        sine * vector.y + cosine * vector.z,
    };
}

Vec3 rotate_y(const Vec3 &vector, const double angle) {
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    return {
        cosine * vector.x + sine * vector.z,
        vector.y,
        -sine * vector.x + cosine * vector.z,
    };
}

Vec3 rotate_z(const Vec3 &vector, const double angle) {
    const double sine = std::sin(angle);
    const double cosine = std::cos(angle);
    return {
        cosine * vector.x - sine * vector.y,
        sine * vector.x + cosine * vector.y,
        vector.z,
    };
}

Vec3 rotate_x_sc(
    const Vec3 &vector,
    const double sine,
    const double cosine
) {
    return {
        vector.x,
        cosine * vector.y - sine * vector.z,
        sine * vector.y + cosine * vector.z,
    };
}

Vec3 rotate_y_sc(
    const Vec3 &vector,
    const double sine,
    const double cosine
) {
    return {
        cosine * vector.x + sine * vector.z,
        vector.y,
        -sine * vector.x + cosine * vector.z,
    };
}

Vec3 rotate_z_sc(
    const Vec3 &vector,
    const double sine,
    const double cosine
) {
    return {
        cosine * vector.x - sine * vector.y,
        sine * vector.x + cosine * vector.y,
        vector.z,
    };
}

enum class CoordinateFrame {
    Lab,
    Alpha,
    Omega,
    Chi,
    Phi,
    Crystal,
    Hkl,
};

CoordinateFrame parse_frame(const std::string &frame) {
    if (frame == "lab" || frame == "q_lab") {
        return CoordinateFrame::Lab;
    }
    if (frame == "alpha" || frame == "q_alpha") {
        return CoordinateFrame::Alpha;
    }
    if (frame == "omega" || frame == "q_omega") {
        return CoordinateFrame::Omega;
    }
    if (frame == "chi" || frame == "q_chi") {
        return CoordinateFrame::Chi;
    }
    if (frame == "phi" || frame == "q_phi") {
        return CoordinateFrame::Phi;
    }
    if (frame == "crystal" || frame == "q_crystal") {
        return CoordinateFrame::Crystal;
    }
    if (frame == "hkl") {
        return CoordinateFrame::Hkl;
    }
    throw py::value_error(
        "frame must be one of lab, alpha, omega, chi, phi, crystal, or hkl"
    );
}

Mat3 matrix_from_array(const FloatArray &array, const char *name) {
    const py::buffer_info info = array.request();
    if (info.ndim != 2 || info.shape[0] != 3 || info.shape[1] != 3) {
        throw py::value_error(std::string(name) + " must have shape (3, 3)");
    }
    const auto *data = static_cast<const double *>(info.ptr);
    Mat3 result;
    std::copy(data, data + 9, result.value.begin());
    return result;
}

std::array<double, 3> triple_from_array(
    const FloatArray &array,
    const char *name,
    const bool strictly_positive
) {
    const py::buffer_info info = array.request();
    if (info.ndim != 1 || info.shape[0] != 3) {
        throw py::value_error(std::string(name) + " must have shape (3,)");
    }
    const auto *data = static_cast<const double *>(info.ptr);
    std::array<double, 3> result{data[0], data[1], data[2]};
    for (const double value : result) {
        if (!std::isfinite(value) || (strictly_positive && value <= 0.0)) {
            throw py::value_error(std::string(name) + " contains an invalid value");
        }
    }
    return result;
}

std::array<std::int64_t, 3> shape_from_array(
    const Int64Array &array,
    const char *name
) {
    const py::buffer_info info = array.request();
    if (info.ndim != 1 || info.shape[0] != 3) {
        throw py::value_error(std::string(name) + " must have shape (3,)");
    }
    const auto *data = static_cast<const std::int64_t *>(info.ptr);
    std::array<std::int64_t, 3> result{data[0], data[1], data[2]};
    for (const std::int64_t value : result) {
        if (value <= 0) {
            throw py::value_error(std::string(name) + " values must be positive");
        }
    }
    return result;
}

struct Grid {
    std::array<double, 3> minimum;
    std::array<double, 3> step;
    std::array<std::int64_t, 3> shape;
    std::array<std::int64_t, 3> chunk_shape;
    std::array<std::uint64_t, 3> chunk_grid;
};

struct RecordKey {
    // Narrower than the flat voxel index they're decomposed from (see
    // record_key()): chunk is bounded by the grid's total chunk count and
    // local by one chunk's own voxel count, both far under 2^32 for any
    // chunk_shape a _GridSpec construction-time check (uint32-overflow
    // validation) allows through.
    std::uint32_t chunk;
    std::uint32_t local;
};

bool operator<(const RecordKey &left, const RecordKey &right) {
    return left.chunk < right.chunk
        || (left.chunk == right.chunk && left.local < right.local);
}

bool operator==(const RecordKey &left, const RecordKey &right) {
    return left.chunk == right.chunk && left.local == right.local;
}

struct Record {
    RecordKey key;
    double weighted_intensity;
    double weighted_variance;
    double weight;
    // Bounded by one tile's pixel count x max subdivision factor, far
    // under 2^32 in any realistic configuration -- see RecordKey's
    // comment for the analogous chunk/local bound.
    std::uint32_t contributors;
};

struct VoxelWeight {
    std::uint64_t voxel;
    double weight;
};

struct CachedVoxel {
    std::uint64_t voxel = 0;
    std::uint64_t generation = 0;
    bool valid = false;
};

struct LatticeVoxel {
    std::uint64_t voxel = 0;
    bool valid = false;
};

struct PixelCoordinateCache {
    std::size_t side = 0;
    std::uint64_t generation = 0;
    std::vector<CachedVoxel> values;

    void begin_pixel(const int max_depth) {
        // A dense dyadic cache is small through depth 3 (17^3 entries).
        // Higher depths retain the uncached path to keep memory bounded.
        if (max_depth > 3) {
            side = 0;
            return;
        }
        const std::size_t required_side =
            (static_cast<std::size_t>(1) << (max_depth + 1)) + 1;
        if (side != required_side) {
            side = required_side;
            values.assign(side * side * side, CachedVoxel{});
        }
        ++generation;
        if (generation == 0) {
            for (auto &value : values) {
                value.generation = 0;
            }
            generation = 1;
        }
    }
};

bool operator<(const VoxelWeight &left, const VoxelWeight &right) {
    return left.voxel < right.voxel;
}

struct Cell {
    double u0;
    double u1;
    double v0;
    double v1;
    double t0;
    double t1;
    double weight;
    int depth;
};

struct BlockProfile {
    std::uint64_t pixels_seen = 0;
    std::uint64_t valid_pixels = 0;
    std::uint64_t coordinate_evaluations = 0;
    std::uint64_t voxel_weights = 0;
    std::uint64_t maximum_weights_per_pixel = 0;
    std::uint64_t unreduced_records = 0;
    std::uint64_t reduced_records = 0;
    std::uint64_t mapping_nanoseconds = 0;
    std::uint64_t reduction_nanoseconds = 0;
};

class ReconstructionKernel {
public:
    ReconstructionKernel(
        const FloatArray &minimum,
        const FloatArray &step,
        const Int64Array &shape,
        const Int64Array &chunk_shape,
        const std::string &frame,
        const double wavevector,
        const FloatArray &ub_inverse,
        const FloatArray &u_inverse,
        const int max_depth = 2,
        int threads = 1,
        const std::size_t work_block_pixels = 4096,
        const std::size_t memory_budget_bytes = 512ULL * 1024ULL * 1024ULL
    )
        : frame_(parse_frame(frame)),
          wavevector_(wavevector),
          ub_inverse_(matrix_from_array(ub_inverse, "ub_inverse")),
          u_inverse_(matrix_from_array(u_inverse, "u_inverse")),
          max_depth_(max_depth),
          threads_(threads),
          work_block_pixels_(work_block_pixels),
          memory_budget_bytes_(memory_budget_bytes) {
        grid_.minimum = triple_from_array(minimum, "minimum", false);
        grid_.step = triple_from_array(step, "step", true);
        grid_.shape = shape_from_array(shape, "shape");
        grid_.chunk_shape = shape_from_array(chunk_shape, "chunk_shape");
        for (int axis = 0; axis < 3; ++axis) {
            grid_.chunk_grid[axis] = static_cast<std::uint64_t>(
                (grid_.shape[axis] + grid_.chunk_shape[axis] - 1)
                / grid_.chunk_shape[axis]
            );
        }
        if (!std::isfinite(wavevector_) || wavevector_ <= 0.0) {
            throw py::value_error("wavevector must be finite and positive");
        }
        if (max_depth_ < 0 || max_depth_ > 8) {
            throw py::value_error("max_depth must be between 0 and 8");
        }
        if (threads_ < 1) {
            threads_ = 1;
        }
        if (work_block_pixels_ < 1) {
            throw py::value_error("work_block_pixels must be positive");
        }
        if (memory_budget_bytes_ < 1024 * 1024) {
            throw py::value_error("memory_budget_bytes must be at least 1 MiB");
        }
    }

    py::dict accumulate(
        const FloatArray &intensity,
        const FloatArray &variance,
        const BoolArray &mask,
        const FloatArray &corner_rays,
        const FloatArray &angles_start,
        const FloatArray &angles_end,
        const bool profile = false
    ) const {
        const auto total_started = std::chrono::steady_clock::now();
        const py::buffer_info intensity_info = intensity.request();
        const py::buffer_info variance_info = variance.request();
        const py::buffer_info mask_info = mask.request();
        const py::buffer_info ray_info = corner_rays.request();
        const py::buffer_info start_info = angles_start.request();
        const py::buffer_info end_info = angles_end.request();
        validate_inputs(
            intensity_info,
            variance_info,
            mask_info,
            ray_info,
            start_info,
            end_info
        );

        const auto *intensity_data = static_cast<const double *>(intensity_info.ptr);
        const auto *variance_data = static_cast<const double *>(variance_info.ptr);
        const auto *mask_data = static_cast<const bool *>(mask_info.ptr);
        const auto *ray_data = static_cast<const double *>(ray_info.ptr);
        const auto *start_data = static_cast<const double *>(start_info.ptr);
        const auto *end_data = static_cast<const double *>(end_info.ptr);
        bool stationary = true;
        for (int index = 0; index < 4; ++index) {
            stationary = stationary && start_data[index] == end_data[index];
        }
        const std::vector<FrameRotation> rotations =
            frame_rotations(start_data, end_data);
        std::vector<CoordinateTransform> transforms;
        transforms.reserve(rotations.size());
        for (const FrameRotation &rotation : rotations) {
            transforms.push_back(coordinate_transform(rotation));
        }
        const std::size_t rows = static_cast<std::size_t>(intensity_info.shape[0]);
        const std::size_t cols = static_cast<std::size_t>(intensity_info.shape[1]);
        const std::size_t pixels = rows * cols;
        std::size_t worst_leaves = 1;
        const std::size_t subdivision_children = stationary ? 4 : 8;
        for (int depth = 0; depth < max_depth_; ++depth) {
            worst_leaves *= subdivision_children;
        }
        const std::size_t estimated_bytes_per_pixel =
            128 + 2 * worst_leaves * sizeof(Record);
        if (
            pixels > 0
            && estimated_bytes_per_pixel
                > memory_budget_bytes_ / pixels
        ) {
            throw py::value_error(
                "Detector tile exceeds the native memory budget at the configured "
                "adaptive depth; use a smaller detector tile or increase "
                "memory_budget_bytes"
            );
        }
        const std::size_t block_size = bounded_block_size();
        const std::size_t blocks = (pixels + block_size - 1) / block_size;
        std::vector<std::vector<Record>> block_results(blocks);
        std::vector<BlockProfile> block_profiles(
            profile ? blocks : 0
        );
        std::atomic<std::size_t> next_block{0};

        const auto blocks_started = std::chrono::steady_clock::now();
        {
            py::gil_scoped_release release;
            const int worker_count = static_cast<int>(
                std::min<std::size_t>(static_cast<std::size_t>(threads_), blocks)
            );
            std::vector<std::thread> workers;
            workers.reserve(static_cast<std::size_t>(worker_count));
            for (int worker = 0; worker < worker_count; ++worker) {
                workers.emplace_back([&, this]() {
                    while (true) {
                        const std::size_t block = next_block.fetch_add(1);
                        if (block >= blocks) {
                            break;
                        }
                        const std::size_t begin = block * block_size;
                        const std::size_t end = std::min(begin + block_size, pixels);
                        block_results[block] = accumulate_block(
                            begin,
                            end,
                            rows,
                            cols,
                            intensity_data,
                            variance_data,
                            mask_data,
                            ray_data,
                            transforms,
                            stationary,
                            profile ? &block_profiles[block] : nullptr
                        );
                    }
                });
            }
            for (auto &worker : workers) {
                worker.join();
            }
        }
        const auto blocks_finished = std::chrono::steady_clock::now();

        std::vector<Record> records;
        std::size_t record_count = 0;
        const auto concatenate_started = std::chrono::steady_clock::now();
        {
            py::gil_scoped_release release;
            for (const auto &block : block_results) {
                record_count += block.size();
            }
            records.reserve(record_count);
            for (auto &block : block_results) {
                records.insert(
                    records.end(),
                    std::make_move_iterator(block.begin()),
                    std::make_move_iterator(block.end())
                );
            }
        }
        const auto concatenate_finished = std::chrono::steady_clock::now();
        const auto final_reduce_started = concatenate_finished;
        {
            py::gil_scoped_release release;
            reduce_records(records);
        }
        const auto final_reduce_finished = std::chrono::steady_clock::now();
        const auto conversion_started = final_reduce_finished;
        py::dict result = records_to_python(records);
        const auto conversion_finished = std::chrono::steady_clock::now();
        if (profile) {
            BlockProfile combined;
            for (const BlockProfile &block : block_profiles) {
                combined.pixels_seen += block.pixels_seen;
                combined.valid_pixels += block.valid_pixels;
                combined.coordinate_evaluations += block.coordinate_evaluations;
                combined.voxel_weights += block.voxel_weights;
                combined.maximum_weights_per_pixel = std::max(
                    combined.maximum_weights_per_pixel,
                    block.maximum_weights_per_pixel
                );
                combined.unreduced_records += block.unreduced_records;
                combined.reduced_records += block.reduced_records;
                combined.mapping_nanoseconds += block.mapping_nanoseconds;
                combined.reduction_nanoseconds += block.reduction_nanoseconds;
            }
            const auto seconds = [](const auto start, const auto stop) {
                return std::chrono::duration<double>(stop - start).count();
            };
            py::dict details;
            details["pixels_seen"] = combined.pixels_seen;
            details["valid_pixels"] = combined.valid_pixels;
            details["coordinate_evaluations"] =
                combined.coordinate_evaluations;
            details["voxel_weights"] = combined.voxel_weights;
            details["maximum_weights_per_pixel"] =
                combined.maximum_weights_per_pixel;
            details["unreduced_block_records"] =
                combined.unreduced_records;
            details["reduced_block_records"] = combined.reduced_records;
            details["concatenated_records"] = record_count;
            details["final_records"] = records.size();
            details["block_mapping_cpu_seconds"] =
                static_cast<double>(combined.mapping_nanoseconds) * 1.0e-9;
            details["block_reduction_cpu_seconds"] =
                static_cast<double>(combined.reduction_nanoseconds) * 1.0e-9;
            details["block_wall_seconds"] =
                seconds(blocks_started, blocks_finished);
            details["concatenate_seconds"] =
                seconds(concatenate_started, concatenate_finished);
            details["final_reduce_seconds"] =
                seconds(final_reduce_started, final_reduce_finished);
            details["python_conversion_seconds"] =
                seconds(conversion_started, conversion_finished);
            details["total_seconds"] =
                seconds(total_started, conversion_finished);
            result["_profile"] = std::move(details);
        }
        return result;
    }

    py::dict configuration() const {
        py::dict result;
        result["max_depth"] = max_depth_;
        result["threads"] = threads_;
        result["work_block_pixels"] = work_block_pixels_;
        result["memory_budget_bytes"] = memory_budget_bytes_;
        return result;
    }

    py::array_t<double> coordinate(
        const FloatArray &corner_rays,
        const FloatArray &angles_start,
        const FloatArray &angles_end,
        const std::size_t row,
        const std::size_t column,
        const double u = 0.5,
        const double v = 0.5,
        const double t = 0.5
    ) const {
        const py::buffer_info ray_info = corner_rays.request();
        const py::buffer_info start_info = angles_start.request();
        const py::buffer_info end_info = angles_end.request();
        if (ray_info.ndim != 3 || ray_info.shape[2] != 3) {
            throw py::value_error("corner_rays must have shape (rows, columns, 3)");
        }
        if (
            ray_info.shape[0] < 2
            || ray_info.shape[1] < 2
            || row + 1 >= static_cast<std::size_t>(ray_info.shape[0])
            || column + 1 >= static_cast<std::size_t>(ray_info.shape[1])
        ) {
            throw py::value_error("Requested pixel is outside corner_rays");
        }
        if (start_info.ndim != 1 || start_info.shape[0] != 4) {
            throw py::value_error("angles_start must have shape (4,)");
        }
        if (end_info.ndim != 1 || end_info.shape[0] != 4) {
            throw py::value_error("angles_end must have shape (4,)");
        }
        for (const double value : {u, v, t}) {
            if (value < 0.0 || value > 1.0 || !std::isfinite(value)) {
                throw py::value_error("u, v, and t must be finite values in [0, 1]");
            }
        }
        const PixelRays prepared_rays = pixel_rays(
            row,
            column,
            static_cast<std::size_t>(ray_info.shape[1] - 1),
            static_cast<const double *>(ray_info.ptr)
        );
        const FrameRotation rotation = frame_rotation(
            t,
            static_cast<const double *>(start_info.ptr),
            static_cast<const double *>(end_info.ptr)
        );
        const Vec3 result = coordinate_at(
            prepared_rays,
            u,
            v,
            rotation
        );
        py::array_t<double> output(3);
        auto *data = static_cast<double *>(output.request().ptr);
        data[0] = result.x;
        data[1] = result.y;
        data[2] = result.z;
        return output;
    }

private:
    Grid grid_;
    CoordinateFrame frame_;
    double wavevector_;
    Mat3 ub_inverse_;
    Mat3 u_inverse_;
    int max_depth_;
    int threads_;
    std::size_t work_block_pixels_;
    std::size_t memory_budget_bytes_;

    std::size_t bounded_block_size() const {
        // A record is deliberately overestimated to retain headroom for the
        // adaptive cell stack and per-pixel temporary weights.
        constexpr std::size_t estimated_bytes_per_pixel = 512;
        const std::size_t per_worker = memory_budget_bytes_
            / static_cast<std::size_t>(std::max(threads_, 1));
        return std::max<std::size_t>(
            1,
            std::min(work_block_pixels_, per_worker / estimated_bytes_per_pixel)
        );
    }

    static void validate_inputs(
        const py::buffer_info &intensity,
        const py::buffer_info &variance,
        const py::buffer_info &mask,
        const py::buffer_info &rays,
        const py::buffer_info &start,
        const py::buffer_info &end
    ) {
        if (intensity.ndim != 2) {
            throw py::value_error("intensity must be a two-dimensional array");
        }
        if (
            variance.ndim != 2
            || variance.shape[0] != intensity.shape[0]
            || variance.shape[1] != intensity.shape[1]
        ) {
            throw py::value_error("variance must have the same shape as intensity");
        }
        if (
            mask.ndim != 2
            || mask.shape[0] != intensity.shape[0]
            || mask.shape[1] != intensity.shape[1]
        ) {
            throw py::value_error("mask must have the same shape as intensity");
        }
        if (
            rays.ndim != 3
            || rays.shape[0] != intensity.shape[0] + 1
            || rays.shape[1] != intensity.shape[1] + 1
            || rays.shape[2] != 3
        ) {
            throw py::value_error(
                "corner_rays must have shape (rows + 1, columns + 1, 3)"
            );
        }
        if (start.ndim != 1 || start.shape[0] != 4) {
            throw py::value_error(
                "angles_start must contain alpha, omega, chi, and phi"
            );
        }
        if (end.ndim != 1 || end.shape[0] != 4) {
            throw py::value_error(
                "angles_end must contain alpha, omega, chi, and phi"
            );
        }
    }

    PixelRays pixel_rays(
        const std::size_t row,
        const std::size_t column,
        const std::size_t columns,
        const double *rays
    ) const {
        const std::size_t stride = columns + 1;
        const auto value = [rays, stride](const std::size_t r, const std::size_t c) {
            const std::size_t offset = (r * stride + c) * 3;
            return Vec3{rays[offset], rays[offset + 1], rays[offset + 2]};
        };
        const Vec3 r00 = value(row, column);
        const Vec3 r10 = value(row + 1, column);
        const Vec3 r01 = value(row, column + 1);
        const Vec3 r11 = value(row + 1, column + 1);
        return {
            r00,
            {
                r10.x - r00.x,
                r10.y - r00.y,
                r10.z - r00.z,
            },
            {
                r01.x - r00.x,
                r01.y - r00.y,
                r01.z - r00.z,
            },
            {
                r11.x - r10.x - r01.x + r00.x,
                r11.y - r10.y - r01.y + r00.y,
                r11.z - r10.z - r01.z + r00.z,
            },
        };
    }

    Vec3 ray_at(
        const PixelRays &rays,
        const double u,
        const double v
    ) const {
        const double uv = u * v;
        Vec3 ray{
            rays.base.x + u * rays.du.x + v * rays.dv.x + uv * rays.duv.x,
            rays.base.y + u * rays.du.y + v * rays.dv.y + uv * rays.duv.y,
            rays.base.z + u * rays.du.z + v * rays.dv.z + uv * rays.duv.z,
        };
        const double norm = std::sqrt(ray.x * ray.x + ray.y * ray.y + ray.z * ray.z);
        if (norm > 0.0) {
            ray.x /= norm;
            ray.y /= norm;
            ray.z /= norm;
        }
        return ray;
    }

    FrameRotation frame_rotation(
        const double t,
        const double *angles_start,
        const double *angles_end
    ) const {
        std::array<double, 4> angles{};
        for (int index = 0; index < 4; ++index) {
            angles[index] = angles_start[index]
                + t * (angles_end[index] - angles_start[index]);
        }
        return {
            std::sin(-angles[0]),
            std::cos(-angles[0]),
            std::sin(angles[1]),
            std::cos(angles[1]),
            std::sin(-angles[2]),
            std::cos(-angles[2]),
            std::sin(-angles[3]),
            std::cos(-angles[3]),
        };
    }

    std::vector<FrameRotation> frame_rotations(
        const double *angles_start,
        const double *angles_end
    ) const {
        const std::size_t side =
            (static_cast<std::size_t>(1) << (max_depth_ + 1)) + 1;
        std::vector<FrameRotation> rotations;
        rotations.reserve(side);
        for (std::size_t index = 0; index < side; ++index) {
            rotations.push_back(
                frame_rotation(
                    static_cast<double>(index)
                        / static_cast<double>(side - 1),
                    angles_start,
                    angles_end
                )
            );
        }
        return rotations;
    }

    Vec3 apply_frame_rotation(
        Vec3 current,
        const FrameRotation &rotation
    ) const {
        if (frame_ == CoordinateFrame::Lab) {
            return current;
        }
        current = rotate_x_sc(
            current,
            rotation.sin_alpha,
            rotation.cos_alpha
        );
        if (frame_ == CoordinateFrame::Alpha) {
            return current;
        }
        current = rotate_z_sc(
            current,
            rotation.sin_omega,
            rotation.cos_omega
        );
        if (frame_ == CoordinateFrame::Omega) {
            return current;
        }
        current = rotate_y_sc(
            current,
            rotation.sin_chi,
            rotation.cos_chi
        );
        if (frame_ == CoordinateFrame::Chi) {
            return current;
        }
        current = rotate_x_sc(
            current,
            rotation.sin_phi,
            rotation.cos_phi
        );
        if (frame_ == CoordinateFrame::Phi) {
            return current;
        }
        if (frame_ == CoordinateFrame::Crystal) {
            return apply(u_inverse_, current);
        }
        return apply(ub_inverse_, current);
    }

    Vec3 coordinate_at(
        const PixelRays &rays,
        const double u,
        const double v,
        const FrameRotation &rotation
    ) const {
        const Vec3 ray = ray_at(rays, u, v);
        Vec3 current{
            wavevector_ * ray.x,
            wavevector_ * (ray.y - 1.0),
            wavevector_ * ray.z,
        };
        return apply_frame_rotation(current, rotation);
    }

    CoordinateTransform coordinate_transform(
        const FrameRotation &rotation
    ) const {
        const Vec3 x_axis = apply_frame_rotation(
            {wavevector_, 0.0, 0.0}, rotation
        );
        const Vec3 y_axis = apply_frame_rotation(
            {0.0, wavevector_, 0.0}, rotation
        );
        const Vec3 z_axis = apply_frame_rotation(
            {0.0, 0.0, wavevector_}, rotation
        );
        return {
            {
                {
                    x_axis.x, y_axis.x, z_axis.x,
                    x_axis.y, y_axis.y, z_axis.y,
                    x_axis.z, y_axis.z, z_axis.z,
                }
            },
            apply_frame_rotation({0.0, -wavevector_, 0.0}, rotation),
        };
    }

    Vec3 coordinate_from_unit_ray(
        const Vec3 &ray,
        const CoordinateTransform &transform
    ) const {
        const Vec3 rotated = apply(transform.matrix, ray);
        return {
            rotated.x + transform.offset.x,
            rotated.y + transform.offset.y,
            rotated.z + transform.offset.z,
        };
    }

    Vec3 coordinate_at(
        const PixelRays &rays,
        const double u,
        const double v,
        const CoordinateTransform &transform
    ) const {
        return coordinate_from_unit_ray(ray_at(rays, u, v), transform);
    }

    bool voxel_id(const Vec3 &coordinate, std::uint64_t &voxel) const {
        std::array<std::int64_t, 3> index{};
        const std::array<double, 3> values{
            coordinate.x,
            coordinate.y,
            coordinate.z,
        };
        for (int axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(values[axis])) {
                return false;
            }
            index[axis] = static_cast<std::int64_t>(
                std::floor((values[axis] - grid_.minimum[axis]) / grid_.step[axis])
            );
            if (index[axis] < 0 || index[axis] >= grid_.shape[axis]) {
                return false;
            }
        }
        voxel = (
            static_cast<std::uint64_t>(index[0])
                * static_cast<std::uint64_t>(grid_.shape[1])
            + static_cast<std::uint64_t>(index[1])
        ) * static_cast<std::uint64_t>(grid_.shape[2])
            + static_cast<std::uint64_t>(index[2]);
        return true;
    }

    RecordKey record_key(const std::uint64_t voxel) const {
        std::uint64_t remaining = voxel;
        const std::uint64_t shape_z =
            static_cast<std::uint64_t>(grid_.shape[2]);
        const std::uint64_t shape_y =
            static_cast<std::uint64_t>(grid_.shape[1]);
        const std::uint64_t index_z = remaining % shape_z;
        remaining /= shape_z;
        const std::uint64_t index_y = remaining % shape_y;
        const std::uint64_t index_x = remaining / shape_y;
        const std::uint64_t chunk_x = static_cast<std::uint64_t>(
            index_x / static_cast<std::uint64_t>(grid_.chunk_shape[0])
        );
        const std::uint64_t chunk_y = static_cast<std::uint64_t>(
            index_y / static_cast<std::uint64_t>(grid_.chunk_shape[1])
        );
        const std::uint64_t chunk_z = static_cast<std::uint64_t>(
            index_z / static_cast<std::uint64_t>(grid_.chunk_shape[2])
        );
        RecordKey key{};
        // Arithmetic stays uint64 throughout (chunk_grid/chunk_shape are
        // uint64 already); only the final, per-axis-bounded result
        // narrows to uint32, safe once the construction-time
        // uint32-overflow validation on chunk_shape/grid.shape holds.
        key.chunk = static_cast<std::uint32_t>((
            chunk_x * grid_.chunk_grid[1] + chunk_y
        ) * grid_.chunk_grid[2] + chunk_z);
        const std::uint64_t local_x = static_cast<std::uint64_t>(
            index_x % static_cast<std::uint64_t>(grid_.chunk_shape[0])
        );
        const std::uint64_t local_y = static_cast<std::uint64_t>(
            index_y % static_cast<std::uint64_t>(grid_.chunk_shape[1])
        );
        const std::uint64_t local_z = static_cast<std::uint64_t>(
            index_z % static_cast<std::uint64_t>(grid_.chunk_shape[2])
        );
        key.local = static_cast<std::uint32_t>((
            local_x * static_cast<std::uint64_t>(grid_.chunk_shape[1]) + local_y
        ) * static_cast<std::uint64_t>(grid_.chunk_shape[2]) + local_z);
        return key;
    }

    bool cached_voxel_key(
        const PixelRays &rays,
        const double u,
        const double v,
        const double t,
        const std::vector<CoordinateTransform> &transforms,
        PixelCoordinateCache &cache,
        std::uint64_t &voxel,
        BlockProfile *profile
    ) const {
        const double rotation_scale =
            static_cast<double>(transforms.size() - 1);
        const std::size_t it = static_cast<std::size_t>(
            std::llround(t * rotation_scale)
        );
        if (cache.side == 0) {
            if (profile != nullptr) {
                ++profile->coordinate_evaluations;
            }
            return voxel_id(
                coordinate_at(
                    rays,
                    u,
                    v,
                    transforms[it]
                ),
                voxel
            );
        }
        const double scale = static_cast<double>(cache.side - 1);
        const std::size_t iu = static_cast<std::size_t>(std::llround(u * scale));
        const std::size_t iv = static_cast<std::size_t>(std::llround(v * scale));
        CachedVoxel &cached = cache.values[
            (iu * cache.side + iv) * cache.side + it
        ];
        if (cached.generation != cache.generation) {
            if (profile != nullptr) {
                ++profile->coordinate_evaluations;
            }
            cached.valid = voxel_id(
                coordinate_at(
                    rays,
                    u,
                    v,
                    transforms[it]
                ),
                cached.voxel
            );
            cached.generation = cache.generation;
        }
        voxel = cached.voxel;
        return cached.valid;
    }

    LatticeVoxel stationary_lattice_voxel(
        const PixelRays &rays,
        const std::size_t u_index,
        const std::size_t v_index,
        const CoordinateTransform &transform,
        BlockProfile *profile
    ) const {
        LatticeVoxel result;
        if (profile != nullptr) {
            ++profile->coordinate_evaluations;
        }
        if (
            (u_index == 0 || u_index == 8)
            && (v_index == 0 || v_index == 8)
        ) {
            Vec3 ray = rays.base;
            if (u_index == 8) {
                ray.x += rays.du.x;
                ray.y += rays.du.y;
                ray.z += rays.du.z;
            }
            if (v_index == 8) {
                ray.x += rays.dv.x;
                ray.y += rays.dv.y;
                ray.z += rays.dv.z;
            }
            if (u_index == 8 && v_index == 8) {
                ray.x += rays.duv.x;
                ray.y += rays.duv.y;
                ray.z += rays.duv.z;
            }
            result.valid = voxel_id(
                coordinate_from_unit_ray(ray, transform),
                result.voxel
            );
            return result;
        }
        result.valid = voxel_id(
            coordinate_at(
                rays,
                static_cast<double>(u_index) * 0.125,
                static_cast<double>(v_index) * 0.125,
                transform
            ),
            result.voxel
        );
        return result;
    }

    static bool same_lattice_voxel(
        const LatticeVoxel &first,
        const LatticeVoxel &second,
        const LatticeVoxel &third,
        const LatticeVoxel &fourth
    ) {
        return first.valid
            && second.valid
            && third.valid
            && fourth.valid
            && first.voxel == second.voxel
            && first.voxel == third.voxel
            && first.voxel == fourth.voxel;
    }

    void split_pixel_stationary_depth2(
        const PixelRays &rays,
        const CoordinateTransform &transform,
        std::vector<VoxelWeight> &weights,
        BlockProfile *profile
    ) const {
        // Depth two uses exact eighth-pixel dyadics. Classify the root first,
        // then materialize only the lattice points needed by unresolved
        // children. This removes recursive stack/cache traffic without
        // evaluating points that a voxel-equality test has already accepted.
        std::array<LatticeVoxel, 25> corners;
        std::uint32_t evaluated = 0;
        const auto corner_index = [](const std::size_t u, const std::size_t v) {
            return (u / 2) * 5 + v / 2;
        };
        const auto evaluate_corner = [&](const std::size_t u, const std::size_t v) {
            const std::size_t index = corner_index(u, v);
            const std::uint32_t bit =
                std::uint32_t{1} << static_cast<unsigned int>(index);
            if ((evaluated & bit) == 0) {
                corners[index] =
                    stationary_lattice_voxel(rays, u, v, transform, profile);
                evaluated |= bit;
            }
        };
        evaluate_corner(0, 0);
        evaluate_corner(8, 0);
        evaluate_corner(0, 8);
        evaluate_corner(8, 8);
        if (
            same_lattice_voxel(
                corners[corner_index(0, 0)],
                corners[corner_index(8, 0)],
                corners[corner_index(0, 8)],
                corners[corner_index(8, 8)]
            )
        ) {
            weights.push_back({corners[corner_index(0, 0)].voxel, 1.0});
            return;
        }

        // Only the 3x3 depth-one lattice is needed to classify the four
        // children. The five additional quarter-pixel nodes inside a child
        // are evaluated below only if that child actually needs splitting.
        for (std::size_t u = 0; u <= 8; u += 4) {
            for (std::size_t v = 0; v <= 8; v += 4) {
                evaluate_corner(u, v);
            }
        }

        for (int child = 0; child < 4; ++child) {
            const std::size_t u0 = (child & 1) != 0 ? 4 : 0;
            const std::size_t v0 = (child & 2) != 0 ? 4 : 0;
            const std::size_t u1 = u0 + 4;
            const std::size_t v1 = v0 + 4;
            if (
                same_lattice_voxel(
                    corners[corner_index(u0, v0)],
                    corners[corner_index(u1, v0)],
                    corners[corner_index(u0, v1)],
                    corners[corner_index(u1, v1)]
                )
            ) {
                weights.push_back({
                    corners[corner_index(u0, v0)].voxel,
                    0.25,
                });
                continue;
            }
            // The four child corners are already available. These are the
            // five new nodes shared by its four grandchildren.
            evaluate_corner(u0 + 2, v0);
            evaluate_corner(u0, v0 + 2);
            evaluate_corner(u0 + 2, v0 + 2);
            evaluate_corner(u1, v0 + 2);
            evaluate_corner(u0 + 2, v1);
            for (int grandchild = 0; grandchild < 4; ++grandchild) {
                const std::size_t child_u0 =
                    u0 + ((grandchild & 1) != 0 ? 2 : 0);
                const std::size_t child_v0 =
                    v0 + ((grandchild & 2) != 0 ? 2 : 0);
                const std::size_t child_u1 = child_u0 + 2;
                const std::size_t child_v1 = child_v0 + 2;
                const LatticeVoxel &first =
                    corners[corner_index(child_u0, child_v0)];
                if (
                    same_lattice_voxel(
                        first,
                        corners[corner_index(child_u1, child_v0)],
                        corners[corner_index(child_u0, child_v1)],
                        corners[corner_index(child_u1, child_v1)]
                    )
                ) {
                    weights.push_back({first.voxel, 0.0625});
                    continue;
                }
                const LatticeVoxel centroid = stationary_lattice_voxel(
                    rays,
                    child_u0 + 1,
                    child_v0 + 1,
                    transform,
                    profile
                );
                if (centroid.valid) {
                    weights.push_back({centroid.voxel, 0.0625});
                }
            }
        }
    }

    void split_pixel(
        const PixelRays &rays,
        const std::vector<CoordinateTransform> &transforms,
        const bool stationary,
        std::vector<VoxelWeight> &weights,
        PixelCoordinateCache &cache,
        std::vector<Cell> &stack,
        BlockProfile *profile
    ) const {
        cache.begin_pixel(max_depth_);
        stack.clear();
        stack.push_back({
            0.0,
            1.0,
            0.0,
            1.0,
            stationary ? 0.5 : 0.0,
            stationary ? 0.5 : 1.0,
            1.0,
            0,
        });
        while (!stack.empty()) {
            const Cell cell = stack.back();
            stack.pop_back();
            std::uint64_t first_voxel = 0;
            bool first_valid = false;
            bool all_same = true;
            const int corner_count = stationary ? 4 : 8;
            for (int corner = 0; corner < corner_count; ++corner) {
                const double u = (corner & 1) != 0 ? cell.u1 : cell.u0;
                const double v = (corner & 2) != 0 ? cell.v1 : cell.v0;
                const double t = stationary
                    ? 0.5
                    : ((corner & 4) != 0 ? cell.t1 : cell.t0);
                std::uint64_t voxel = 0;
                const bool valid = cached_voxel_key(
                    rays,
                    u,
                    v,
                    t,
                    transforms,
                    cache,
                    voxel,
                    profile
                );
                if (corner == 0) {
                    first_valid = valid;
                    first_voxel = voxel;
                } else if (
                    valid != first_valid
                    || (valid && voxel != first_voxel)
                ) {
                    all_same = false;
                }
            }
            if (all_same && first_valid) {
                weights.push_back({first_voxel, cell.weight});
                continue;
            }
            if (cell.depth >= max_depth_) {
                std::uint64_t voxel = 0;
                const bool valid = cached_voxel_key(
                    rays,
                    0.5 * (cell.u0 + cell.u1),
                    0.5 * (cell.v0 + cell.v1),
                    0.5 * (cell.t0 + cell.t1),
                    transforms,
                    cache,
                    voxel,
                    profile
                );
                if (valid) {
                    weights.push_back({voxel, cell.weight});
                }
                continue;
            }
            const double um = 0.5 * (cell.u0 + cell.u1);
            const double vm = 0.5 * (cell.v0 + cell.v1);
            const double tm = 0.5 * (cell.t0 + cell.t1);
            const int child_count = stationary ? 4 : 8;
            const double child_weight =
                cell.weight / static_cast<double>(child_count);
            for (int child = child_count - 1; child >= 0; --child) {
                stack.push_back({
                    (child & 1) != 0 ? um : cell.u0,
                    (child & 1) != 0 ? cell.u1 : um,
                    (child & 2) != 0 ? vm : cell.v0,
                    (child & 2) != 0 ? cell.v1 : vm,
                    stationary
                        ? 0.5
                        : ((child & 4) != 0 ? tm : cell.t0),
                    stationary
                        ? 0.5
                        : ((child & 4) != 0 ? cell.t1 : tm),
                    child_weight,
                    cell.depth + 1,
                });
            }
        }
    }

    std::vector<Record> accumulate_block(
        const std::size_t begin,
        const std::size_t end,
        const std::size_t rows,
        const std::size_t columns,
        const double *intensity,
        const double *variance,
        const bool *mask,
        const double *rays,
        const std::vector<CoordinateTransform> &transforms,
        const bool stationary,
        BlockProfile *profile
    ) const {
        (void)rows;
        std::vector<Record> records;
        std::vector<VoxelWeight> weights;
        PixelCoordinateCache coordinate_cache;
        std::vector<Cell> stack;
        const std::size_t reserve_leaves = static_cast<std::size_t>(1)
            << std::min((stationary ? 2 : 3) * max_depth_, 9);
        records.reserve(2 * (end - begin));
        weights.reserve(reserve_leaves);
        stack.reserve(reserve_leaves);
        const auto mapping_started = std::chrono::steady_clock::now();
        for (std::size_t flat = begin; flat < end; ++flat) {
            if (profile != nullptr) {
                ++profile->pixels_seen;
            }
            if (
                mask[flat]
                || !std::isfinite(intensity[flat])
                || !std::isfinite(variance[flat])
                || variance[flat] < 0.0
            ) {
                continue;
            }
            if (profile != nullptr) {
                ++profile->valid_pixels;
            }
            const std::size_t row = flat / columns;
            const std::size_t column = flat % columns;
            const PixelRays prepared_rays = pixel_rays(
                row,
                column,
                columns,
                rays
            );
            weights.clear();
            if (max_depth_ == 0) {
                if (profile != nullptr) {
                    ++profile->coordinate_evaluations;
                }
                std::uint64_t voxel = 0;
                if (
                    voxel_id(
                        coordinate_at(
                            prepared_rays,
                            0.5,
                            0.5,
                            transforms[transforms.size() / 2]
                        ),
                        voxel
                    )
                ) {
                    weights.push_back({voxel, 1.0});
                }
            } else if (stationary && max_depth_ == 2) {
                split_pixel_stationary_depth2(
                    prepared_rays,
                    transforms[transforms.size() / 2],
                    weights,
                    profile
                );
            } else {
                split_pixel(
                    prepared_rays,
                    transforms,
                    stationary,
                    weights,
                    coordinate_cache,
                    stack,
                    profile
                );
            }
            if (profile != nullptr) {
                profile->voxel_weights += weights.size();
                profile->maximum_weights_per_pixel = std::max(
                    profile->maximum_weights_per_pixel,
                    static_cast<std::uint64_t>(weights.size())
                );
            }
            std::sort(weights.begin(), weights.end());
            std::size_t offset = 0;
            while (offset < weights.size()) {
                const std::uint64_t voxel = weights[offset].voxel;
                double weight = 0.0;
                do {
                    weight += weights[offset].weight;
                    ++offset;
                } while (
                    offset < weights.size()
                    && weights[offset].voxel == voxel
                );
                const RecordKey key = record_key(voxel);
                records.push_back({
                    key,
                    weight * intensity[flat],
                    weight * weight * variance[flat],
                    weight,
                    1,
                });
            }
        }
        const auto mapping_finished = std::chrono::steady_clock::now();
        if (profile != nullptr) {
            profile->unreduced_records = records.size();
            profile->mapping_nanoseconds = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    mapping_finished - mapping_started
                ).count()
            );
        }
        const auto reduction_started = mapping_finished;
        reduce_records(records);
        if (profile != nullptr) {
            profile->reduced_records = records.size();
            profile->reduction_nanoseconds = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - reduction_started
                ).count()
            );
        }
        return records;
    }

    static void reduce_records(std::vector<Record> &records) {
        std::stable_sort(
            records.begin(),
            records.end(),
            [](const Record &left, const Record &right) {
                return left.key < right.key;
            }
        );
        std::size_t write = 0;
        std::size_t read = 0;
        while (read < records.size()) {
            Record combined = records[read++];
            while (read < records.size() && records[read].key == combined.key) {
                combined.weighted_intensity += records[read].weighted_intensity;
                combined.weighted_variance += records[read].weighted_variance;
                combined.weight += records[read].weight;
                combined.contributors += records[read].contributors;
                ++read;
            }
            records[write++] = combined;
        }
        records.resize(write);
    }

    static py::dict records_to_python(const std::vector<Record> &records) {
        const py::ssize_t size = static_cast<py::ssize_t>(records.size());
        UInt32Array chunk_id(size);
        UInt32Array local_voxel_id(size);
        FloatArray weighted_intensity(size);
        FloatArray weighted_variance(size);
        FloatArray weight(size);
        UInt32Array contributors(size);
        auto *chunk_data = static_cast<std::uint32_t *>(chunk_id.request().ptr);
        auto *local_data = static_cast<std::uint32_t *>(local_voxel_id.request().ptr);
        auto *intensity_data = static_cast<double *>(weighted_intensity.request().ptr);
        auto *variance_data = static_cast<double *>(weighted_variance.request().ptr);
        auto *weight_data = static_cast<double *>(weight.request().ptr);
        auto *contributors_data = static_cast<std::uint32_t *>(
            contributors.request().ptr
        );
        {
            py::gil_scoped_release release;
            for (py::ssize_t index = 0; index < size; ++index) {
                const Record &record = records[static_cast<std::size_t>(index)];
                chunk_data[index] = record.key.chunk;
                local_data[index] = record.key.local;
                intensity_data[index] = record.weighted_intensity;
                variance_data[index] = record.weighted_variance;
                weight_data[index] = record.weight;
                contributors_data[index] = record.contributors;
            }
        }
        py::dict result;
        result["chunk_id"] = std::move(chunk_id);
        result["local_voxel_id"] = std::move(local_voxel_id);
        result["weighted_intensity"] = std::move(weighted_intensity);
        result["weighted_variance"] = std::move(weighted_variance);
        result["weight"] = std::move(weight);
        result["contributors"] = std::move(contributors);
        return result;
    }
};

PYBIND11_MODULE(_reciprocal_reconstruction_cpp, module) {
    module.doc() =
        "Native reciprocal-space coordinate conversion and footprint accumulation.";
    module.def(
        "xxh3_128",
        &xxh3_128_buffer,
        py::arg("buffer"),
        "Return the canonical XXH3-128 digest of a C-contiguous buffer."
    );
    module.def(
        "merge_sorted_batches",
        &merge_sorted_batches,
        py::arg("left_chunk_id"),
        py::arg("left_local_voxel_id"),
        py::arg("left_weighted_intensity"),
        py::arg("left_weighted_variance"),
        py::arg("left_weight"),
        py::arg("left_contributors"),
        py::arg("right_chunk_id"),
        py::arg("right_local_voxel_id"),
        py::arg("right_weighted_intensity"),
        py::arg("right_weighted_variance"),
        py::arg("right_weight"),
        py::arg("right_contributors"),
        "Linearly merge two sorted, already-reduced record batches."
    );
    py::class_<ReconstructionKernel>(module, "ReconstructionKernel")
        .def(
            py::init<
                const FloatArray &,
                const FloatArray &,
                const Int64Array &,
                const Int64Array &,
                const std::string &,
                double,
                const FloatArray &,
                const FloatArray &,
                int,
                int,
                std::size_t,
                std::size_t
            >(),
            py::arg("minimum"),
            py::arg("step"),
            py::arg("shape"),
            py::arg("chunk_shape"),
            py::arg("frame"),
            py::arg("wavevector"),
            py::arg("ub_inverse"),
            py::arg("u_inverse"),
            py::arg("max_depth") = 2,
            py::arg("threads") = 1,
            py::arg("work_block_pixels") = 4096,
            py::arg("memory_budget_bytes") = 512ULL * 1024ULL * 1024ULL
        )
        .def(
            "accumulate",
            &ReconstructionKernel::accumulate,
            py::arg("intensity"),
            py::arg("variance"),
            py::arg("mask"),
            py::arg("corner_rays"),
            py::arg("angles_start"),
            py::arg("angles_end"),
            py::arg("profile") = false
        )
        .def(
            "coordinate",
            &ReconstructionKernel::coordinate,
            py::arg("corner_rays"),
            py::arg("angles_start"),
            py::arg("angles_end"),
            py::arg("row"),
            py::arg("column"),
            py::arg("u") = 0.5,
            py::arg("v") = 0.5,
            py::arg("t") = 0.5
        )
        .def_property_readonly(
            "configuration",
            &ReconstructionKernel::configuration
        );
}

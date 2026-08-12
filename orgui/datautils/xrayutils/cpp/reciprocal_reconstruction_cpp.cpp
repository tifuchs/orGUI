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
#include <map>
#include <memory>
#include <memory_resource>
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
    // Single pass: allocate at the worst case (no duplicates -- output can
    // only be smaller than left_size + right_size when equal keys merge),
    // merge directly into it, then shrink to the true count. Replaces a
    // prior two-pass count-then-fill structure that walked the same
    // comparison logic twice per call.
    const py::ssize_t worst_case_size = left_size + right_size;
    ContiguousUInt32Array output_chunk(worst_case_size);
    ContiguousUInt32Array output_local(worst_case_size);
    FloatArray output_intensity(worst_case_size);
    FloatArray output_variance(worst_case_size);
    FloatArray output_weight(worst_case_size);
    ContiguousUInt32Array output_contributors(worst_case_size);
    auto *oc_data = static_cast<std::uint32_t *>(output_chunk.request().ptr);
    auto *ol_data = static_cast<std::uint32_t *>(output_local.request().ptr);
    auto *oi_data = static_cast<double *>(output_intensity.request().ptr);
    auto *ov_data = static_cast<double *>(output_variance.request().ptr);
    auto *ow_data = static_cast<double *>(output_weight.request().ptr);
    auto *on_data = static_cast<std::uint32_t *>(
        output_contributors.request().ptr
    );
    py::ssize_t output_size = 0;
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
        output_size = output;
    }
    if (output_size != worst_case_size) {
        // Freshly allocated above, not yet exposed to Python -- refcheck
        // would only ever see this one reference, but pass false directly
        // since that's the actual invariant being relied on.
        const std::vector<py::ssize_t> trimmed_shape{output_size};
        output_chunk.resize(trimmed_shape, false);
        output_local.resize(trimmed_shape, false);
        output_intensity.resize(trimmed_shape, false);
        output_variance.resize(trimmed_shape, false);
        output_weight.resize(trimmed_shape, false);
        output_contributors.resize(trimmed_shape, false);
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
    // shape as doubles, so voxel_id()'s per-axis range test compares in
    // floating point (before any integer conversion) rather than casting a
    // possibly out-of-range floor() result to int64 first.
    std::array<double, 3> shape_as_double;
    // Bit field layout of the internal flat voxel id (see voxel_id()): one
    // field per axis, wide enough for that axis' voxel count, packed
    // x:y:z from the most significant end. Bit-packing keeps the id
    // lexicographically ordered by (x, y, z) exactly as the row-major
    // encoding it replaces, so every ordering built on it is unchanged --
    // but record_key() can unpack it with shifts instead of divisions.
    std::array<int, 3> voxel_shift;
    std::array<std::uint64_t, 3> voxel_mask;
    // Per-axis chunk split. chunk_power_of_two selects shifts/masks over
    // divisions; every default chunk_shape (64, 64, 64) takes that path.
    std::array<int, 3> chunk_shift;
    std::array<std::uint64_t, 3> chunk_mask;
    bool chunk_power_of_two = false;
};

// Smallest `bits` with 2^bits >= size: the width of an index into `size`
// values, and, for a power-of-two `size`, exactly log2(size).
int ceil_log2(const std::int64_t size) {
    int bits = 0;
    while (
        bits < 63
        && (static_cast<std::uint64_t>(1) << bits)
            < static_cast<std::uint64_t>(size)
    ) {
        ++bits;
    }
    return bits;
}

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

struct RecordAccum {
    double weighted_intensity = 0.0;
    double weighted_variance = 0.0;
    double weight = 0.0;
    std::uint32_t contributors = 0;
};

// Records one pixel is assumed to leave behind, for every purpose that has
// to size memory before the data is seen: the per-worker arena, and the
// per-call budget precheck.
//
// Sizing by leaf count instead assumes every leaf of every pixel reaches a
// voxel no other leaf did. Real data does not behave that way at any depth:
// measured post-dedup density stays between 0.46 and 0.89 records per pixel
// from depth 0 through 8, and across both exposure models, because deeper
// subdivision's extra samples overwhelmingly land in voxels a neighbouring
// sample already reached. Leaves are transient -- they are sorted, merged
// and discarded per pixel -- while records are what stays resident.
//
// Four is therefore several times the density ever observed, and both users
// below take it as a ceiling on the leaf count rather than a replacement,
// so depth 0 stays exact: one leaf can only yield one record.
constexpr std::size_t reserved_records_per_pixel = 4;

struct VoxelWeight {
    std::uint64_t voxel;
    double weight;
};

struct LatticeVoxel {
    std::uint64_t voxel = 0;
    bool valid = false;
};

bool operator<(const VoxelWeight &left, const VoxelWeight &right) {
    return left.voxel < right.voxel;
}

// One work block of a frame group: a rectangle of the detector, taken
// across every frame of the group.
struct BrickExtent {
    std::size_t row_begin;
    std::size_t row_end;
    std::size_t column_begin;
    std::size_t column_end;
};

// A frame group's pixel data as the caller already holds it: one pointer
// per frame at the tile's first pixel, plus the distance between two
// rows of the caller's own array.
//
// The stride is what makes the tile a view rather than a copy. A tile cut
// out of whole corrected frames strides by the frame's row length and is
// read in place; a tile that has been gathered into its own
// ``(frames, rows, columns)`` buffer strides by the tile's own row
// length. The brick loop cannot tell the two apart, and reads the same
// values in the same order either way.
struct GroupPixels {
    std::vector<const double *> intensity;
    std::vector<const double *> variance;
    std::vector<const bool *> mask;
    std::size_t row_stride = 0;
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
        std::array<int, 3> voxel_bits{};
        int total_voxel_bits = 0;
        grid_.chunk_power_of_two = true;
        for (int axis = 0; axis < 3; ++axis) {
            grid_.chunk_grid[axis] = static_cast<std::uint64_t>(
                (grid_.shape[axis] + grid_.chunk_shape[axis] - 1)
                / grid_.chunk_shape[axis]
            );
            grid_.shape_as_double[axis] =
                static_cast<double>(grid_.shape[axis]);
            voxel_bits[axis] = std::max(1, ceil_log2(grid_.shape[axis]));
            total_voxel_bits += voxel_bits[axis];
            grid_.voxel_mask[axis] =
                (static_cast<std::uint64_t>(1) << voxel_bits[axis]) - 1;
            const std::int64_t chunk_size = grid_.chunk_shape[axis];
            if ((chunk_size & (chunk_size - 1)) != 0) {
                grid_.chunk_power_of_two = false;
            } else {
                grid_.chunk_shift[axis] = ceil_log2(chunk_size);
                grid_.chunk_mask[axis] =
                    static_cast<std::uint64_t>(chunk_size) - 1;
            }
        }
        grid_.voxel_shift[2] = 0;
        grid_.voxel_shift[1] = voxel_bits[2];
        grid_.voxel_shift[0] = voxel_bits[2] + voxel_bits[1];
        if (total_voxel_bits > 64) {
            // Only reachable for grids within a factor of eight of
            // overflowing a 64-bit voxel count outright, which no
            // physically meaningful reconstruction approaches.
            throw py::value_error(
                "grid shape needs more than 64 bits of voxel index; use a "
                "coarser step or a smaller grid"
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
        // Bounded by reserved_records_per_pixel for the reason given there.
        // Unbounded, this term is 4^depth (or 8^depth when the exposure
        // rotates) and rejects any useful detector tile from depth 3 up: at
        // depth 3 it claims 5248 bytes for a pixel that really costs about
        // 106, so a full Pilatus-6M frame "needs" 32.7 GB and the layout is
        // forced onto tiny tiles -- which then caps how many frames can be
        // in flight, since one worker's estimate divides into the budget.
        // Two resident copies of a pixel's records are still allowed for,
        // covering the per-block vectors and the merged output together.
        const std::size_t estimated_bytes_per_pixel =
            128
            + 2
                * std::min<std::size_t>(worst_leaves, reserved_records_per_pixel)
                * sizeof(Record);
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
        // One arena per worker thread, sized for one block's worst case
        // (every leaf a distinct voxel) and reused -- via
        // monotonic_buffer_resource::release(), not a fresh allocation --
        // across every block that thread processes. Falls back to the
        // default heap resource transparently if a pathological block
        // ever needs more than this, so undersizing here is a
        // performance risk, never a correctness one.
        const std::size_t bytes_per_node =
            sizeof(std::pair<const RecordKey, RecordAccum>) + 32;
        // Sized by reserved_records_per_pixel, as the precheck above is.
        // Undersizing here is a performance risk and never a correctness
        // one: the resource falls back to the heap transparently if a
        // block ever does exceed it.
        const std::size_t arena_bytes =
            block_size
                * std::min<std::size_t>(worst_leaves, reserved_records_per_pixel)
                * bytes_per_node
            + 4096;

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
                    std::unique_ptr<std::byte[]> arena_buffer(
                        new std::byte[arena_bytes]
                    );
                    std::pmr::monotonic_buffer_resource arena(
                        arena_buffer.get(), arena_bytes
                    );
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
                            profile ? &block_profiles[block] : nullptr,
                            arena
                        );
                    }
                });
            }
            for (auto &worker : workers) {
                worker.join();
            }
        }
        const auto blocks_finished = std::chrono::steady_clock::now();

        std::size_t record_count = 0;
        for (const auto &block : block_results) {
            record_count += block.size();
        }
        std::vector<Record> records;
        const auto merge_started = std::chrono::steady_clock::now();
        {
            py::gil_scoped_release release;
            // Each block's output is already sorted and deduped
            // (accumulate_block's pmr::map reduce) -- merge them directly
            // with the loser tree instead of concatenating into one
            // vector and blindly re-sorting it from scratch, which
            // ignores that pre-sortedness entirely.
            records = merge_sorted_blocks(block_results);
        }
        const auto merge_finished = std::chrono::steady_clock::now();
        py::dict result = records_to_python(records);
        const auto conversion_finished = std::chrono::steady_clock::now();
        if (profile) {
            result["_profile"] = summarize_profile(
                block_profiles,
                record_count,
                records.size(),
                total_started,
                blocks_started,
                blocks_finished,
                merge_started,
                merge_finished,
                conversion_finished
            );
        }
        return result;
    }

    py::dict accumulate_group(
        const FloatArray &intensity,
        const FloatArray &variance,
        const BoolArray &mask,
        const FloatArray &corner_rays,
        const FloatArray &angles_start,
        const FloatArray &angles_end,
        const bool profile = false
    ) const {
        // Several consecutive images mapped in one call, with the work
        // block a brick in (row, column, frame) rather than a run of one
        // image. On a rotation scan two adjacent frames are as close in
        // reciprocal space as two adjacent pixels are -- measured 0.71
        // voxels per frame against 0.49 and 0.80 for the next column and
        // row -- so a brick collapses redundancy the per-image block
        // cannot see, and it does so in the same cache-resident map.
        //
        // The caller chooses the group size by how many frames it passes;
        // one frame reproduces accumulate() exactly.
        //
        // This form takes one contiguous ``(frames, rows, columns)``
        // buffer per array. :meth:`accumulate_group_tile` takes the same
        // data as whole frames plus a rectangle, and is what the mapping
        // pipeline uses; both run the identical brick loop over the
        // identical values.
        const auto total_started = std::chrono::steady_clock::now();
        const py::buffer_info intensity_info = intensity.request();
        const py::buffer_info variance_info = variance.request();
        const py::buffer_info mask_info = mask.request();
        const py::buffer_info ray_info = corner_rays.request();
        const py::buffer_info start_info = angles_start.request();
        const py::buffer_info end_info = angles_end.request();
        validate_group_inputs(
            intensity_info,
            variance_info,
            mask_info,
            ray_info,
            start_info,
            end_info
        );

        const auto *intensity_data =
            static_cast<const double *>(intensity_info.ptr);
        const auto *variance_data =
            static_cast<const double *>(variance_info.ptr);
        const auto *mask_data = static_cast<const bool *>(mask_info.ptr);

        const std::size_t frames =
            static_cast<std::size_t>(intensity_info.shape[0]);
        const std::size_t rows =
            static_cast<std::size_t>(intensity_info.shape[1]);
        const std::size_t cols =
            static_cast<std::size_t>(intensity_info.shape[2]);

        GroupPixels pixels;
        pixels.row_stride = cols;
        pixels.intensity.reserve(frames);
        pixels.variance.reserve(frames);
        pixels.mask.reserve(frames);
        for (std::size_t frame = 0; frame < frames; ++frame) {
            const std::size_t offset = frame * rows * cols;
            pixels.intensity.push_back(intensity_data + offset);
            pixels.variance.push_back(variance_data + offset);
            pixels.mask.push_back(mask_data + offset);
        }
        return accumulate_group_impl(
            pixels,
            frames,
            rows,
            cols,
            static_cast<const double *>(ray_info.ptr),
            static_cast<const double *>(start_info.ptr),
            static_cast<const double *>(end_info.ptr),
            profile,
            total_started
        );
    }

    py::dict accumulate_group_tile(
        const py::sequence &intensity,
        const py::sequence &variance,
        const py::sequence &mask,
        const FloatArray &corner_rays,
        const FloatArray &angles_start,
        const FloatArray &angles_end,
        const std::size_t row_start,
        const std::size_t row_stop,
        const std::size_t column_start,
        const std::size_t column_stop,
        const bool profile = false
    ) const {
        // The same call as :meth:`accumulate_group`, taking whole frames
        // and the rectangle to map instead of a buffer already cut down
        // to it.
        //
        // Detector tiles partition the detector, so gathering each tile
        // into its own contiguous buffer copies every corrected frame
        // exactly once per group -- around 105 MB a frame, in Python,
        // holding the GIL, purely to give this call a shape it does not
        // need. The brick loop walks rows within a tile in either form,
        // so reading through the frame's own row stride costs it nothing.
        const auto total_started = std::chrono::steady_clock::now();
        const py::buffer_info ray_info = corner_rays.request();
        const py::buffer_info start_info = angles_start.request();
        const py::buffer_info end_info = angles_end.request();

        const std::size_t frames = static_cast<std::size_t>(py::len(intensity));
        if (frames < 1) {
            throw py::value_error("A frame group needs at least one frame");
        }
        if (
            static_cast<std::size_t>(py::len(variance)) != frames
            || static_cast<std::size_t>(py::len(mask)) != frames
        ) {
            throw py::value_error(
                "variance and mask must hold one frame each, as intensity does"
            );
        }
        if (row_stop <= row_start || column_stop <= column_start) {
            throw py::value_error("The tile rectangle must be non-empty");
        }
        const std::size_t rows = row_stop - row_start;
        const std::size_t cols = column_stop - column_start;

        // The casts hold references for as long as this call runs, which
        // is what keeps the frames alive while the GIL is released.
        std::vector<FloatArray> intensity_frames;
        std::vector<FloatArray> variance_frames;
        std::vector<BoolArray> mask_frames;
        intensity_frames.reserve(frames);
        variance_frames.reserve(frames);
        mask_frames.reserve(frames);
        GroupPixels pixels;
        pixels.intensity.reserve(frames);
        pixels.variance.reserve(frames);
        pixels.mask.reserve(frames);
        std::size_t frame_columns = 0;
        for (std::size_t frame = 0; frame < frames; ++frame) {
            intensity_frames.push_back(intensity[frame].cast<FloatArray>());
            variance_frames.push_back(variance[frame].cast<FloatArray>());
            mask_frames.push_back(mask[frame].cast<BoolArray>());
            const py::buffer_info intensity_info =
                intensity_frames.back().request();
            const py::buffer_info variance_info =
                variance_frames.back().request();
            const py::buffer_info mask_info = mask_frames.back().request();
            validate_tile_frame(
                intensity_info, variance_info, mask_info, row_stop, column_stop
            );
            const std::size_t columns_here =
                static_cast<std::size_t>(intensity_info.shape[1]);
            if (frame == 0) {
                frame_columns = columns_here;
            } else if (columns_here != frame_columns) {
                throw py::value_error(
                    "every frame in a group must have the same shape"
                );
            }
            const std::size_t offset = row_start * frame_columns + column_start;
            pixels.intensity.push_back(
                static_cast<const double *>(intensity_info.ptr) + offset
            );
            pixels.variance.push_back(
                static_cast<const double *>(variance_info.ptr) + offset
            );
            pixels.mask.push_back(
                static_cast<const bool *>(mask_info.ptr) + offset
            );
        }
        pixels.row_stride = frame_columns;
        validate_group_geometry(ray_info, start_info, end_info, frames, rows, cols);
        return accumulate_group_impl(
            pixels,
            frames,
            rows,
            cols,
            static_cast<const double *>(ray_info.ptr),
            static_cast<const double *>(start_info.ptr),
            static_cast<const double *>(end_info.ptr),
            profile,
            total_started
        );
    }

    py::dict accumulate_group_impl(
        const GroupPixels &pixels,
        const std::size_t frames,
        const std::size_t rows,
        const std::size_t cols,
        const double *ray_data,
        const double *start_data,
        const double *end_data,
        const bool profile,
        const std::chrono::steady_clock::time_point total_started
    ) const {
        const std::size_t frame_pixels = rows * cols;
        const std::size_t samples = frames * frame_pixels;

        // Exposure geometry is per frame: a group may mix stationary and
        // swept frames, and each needs its own rotation lattice.
        std::vector<std::vector<CoordinateTransform>> transforms(frames);
        std::vector<char> stationary(frames, 1);
        for (std::size_t frame = 0; frame < frames; ++frame) {
            const double *frame_start = start_data + frame * 4;
            const double *frame_end = end_data + frame * 4;
            bool still = true;
            for (int index = 0; index < 4; ++index) {
                still = still && frame_start[index] == frame_end[index];
            }
            stationary[frame] = still ? 1 : 0;
            for (const FrameRotation &rotation :
                 frame_rotations(frame_start, frame_end)) {
                transforms[frame].push_back(coordinate_transform(rotation));
            }
        }

        std::size_t worst_leaves = 1;
        for (int depth = 0; depth < max_depth_; ++depth) {
            worst_leaves *= 8;
        }
        const std::size_t estimated_bytes_per_pixel =
            128
            + 2
                * std::min<std::size_t>(worst_leaves, reserved_records_per_pixel)
                * sizeof(Record);
        if (
            samples > 0
            && estimated_bytes_per_pixel > memory_budget_bytes_ / samples
        ) {
            throw py::value_error(
                "Frame group exceeds the native memory budget at the "
                "configured adaptive depth; use fewer frames per group, a "
                "smaller detector tile, or increase memory_budget_bytes"
            );
        }

        // One brick holds the same sample count a single-image block would,
        // so the map's working set is unchanged; the frame extent is spent
        // out of the detector-plane extent rather than added to it.
        const std::size_t samples_per_brick = bounded_block_size();
        const std::size_t plane_pixels = std::max<std::size_t>(
            1, samples_per_brick / std::max<std::size_t>(frames, 1)
        );
        std::size_t brick_columns = std::max<std::size_t>(
            1,
            static_cast<std::size_t>(
                std::sqrt(static_cast<double>(plane_pixels))
            )
        );
        brick_columns = std::min(brick_columns, cols);
        std::size_t brick_rows =
            std::max<std::size_t>(1, plane_pixels / brick_columns);
        brick_rows = std::min(brick_rows, rows);

        std::vector<BrickExtent> bricks;
        for (std::size_t row = 0; row < rows; row += brick_rows) {
            for (std::size_t column = 0; column < cols; column += brick_columns) {
                bricks.push_back({
                    row,
                    std::min(row + brick_rows, rows),
                    column,
                    std::min(column + brick_columns, cols),
                });
            }
        }
        const std::size_t brick_count = bricks.size();
        std::vector<std::vector<Record>> block_results(brick_count);
        std::vector<BlockProfile> block_profiles(profile ? brick_count : 0);
        std::atomic<std::size_t> next_brick{0};

        const std::size_t bytes_per_node =
            sizeof(std::pair<const RecordKey, RecordAccum>) + 32;
        const std::size_t arena_bytes =
            brick_rows * brick_columns * frames
                * std::min<std::size_t>(worst_leaves, reserved_records_per_pixel)
                * bytes_per_node
            + 4096;

        const auto blocks_started = std::chrono::steady_clock::now();
        {
            py::gil_scoped_release release;
            const int worker_count = static_cast<int>(
                std::min<std::size_t>(
                    static_cast<std::size_t>(threads_), std::max<std::size_t>(brick_count, 1)
                )
            );
            std::vector<std::thread> workers;
            workers.reserve(static_cast<std::size_t>(worker_count));
            for (int worker = 0; worker < worker_count; ++worker) {
                workers.emplace_back([&, this]() {
                    std::unique_ptr<std::byte[]> arena_buffer(
                        new std::byte[arena_bytes]
                    );
                    std::pmr::monotonic_buffer_resource arena(
                        arena_buffer.get(), arena_bytes
                    );
                    std::vector<Vec3> centre_rays;
                    std::vector<VoxelWeight> weights;
                    weights.reserve(
                        static_cast<std::size_t>(1)
                        << std::min(3 * max_depth_, 9)
                    );
                    while (true) {
                        const std::size_t index = next_brick.fetch_add(1);
                        if (index >= brick_count) {
                            break;
                        }
                        block_results[index] = accumulate_brick(
                            bricks[index],
                            frames,
                            cols,
                            pixels,
                            ray_data,
                            transforms,
                            stationary,
                            centre_rays,
                            weights,
                            profile ? &block_profiles[index] : nullptr,
                            arena
                        );
                    }
                });
            }
            for (auto &worker : workers) {
                worker.join();
            }
        }
        const auto blocks_finished = std::chrono::steady_clock::now();

        std::size_t record_count = 0;
        for (const auto &block : block_results) {
            record_count += block.size();
        }
        std::vector<Record> records;
        const auto merge_started = std::chrono::steady_clock::now();
        {
            py::gil_scoped_release release;
            records = merge_sorted_blocks(block_results);
        }
        const auto merge_finished = std::chrono::steady_clock::now();
        py::dict result = records_to_python(records);
        const auto conversion_finished = std::chrono::steady_clock::now();
        if (profile) {
            py::dict details = summarize_profile(
                block_profiles,
                record_count,
                records.size(),
                total_started,
                blocks_started,
                blocks_finished,
                merge_started,
                merge_finished,
                conversion_finished
            );
            details["frames"] = frames;
            details["bricks"] = brick_count;
            details["brick_rows"] = brick_rows;
            details["brick_columns"] = brick_columns;
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

    static py::dict summarize_profile(
        const std::vector<BlockProfile> &block_profiles,
        const std::size_t pre_merge_records,
        const std::size_t final_records,
        const std::chrono::steady_clock::time_point total_started,
        const std::chrono::steady_clock::time_point blocks_started,
        const std::chrono::steady_clock::time_point blocks_finished,
        const std::chrono::steady_clock::time_point merge_started,
        const std::chrono::steady_clock::time_point merge_finished,
        const std::chrono::steady_clock::time_point conversion_finished
    ) {
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
        details["coordinate_evaluations"] = combined.coordinate_evaluations;
        details["voxel_weights"] = combined.voxel_weights;
        details["maximum_weights_per_pixel"] =
            combined.maximum_weights_per_pixel;
        details["unreduced_block_records"] = combined.unreduced_records;
        details["reduced_block_records"] = combined.reduced_records;
        details["pre_merge_records"] = pre_merge_records;
        details["final_records"] = final_records;
        details["block_mapping_cpu_seconds"] =
            static_cast<double>(combined.mapping_nanoseconds) * 1.0e-9;
        details["block_reduction_cpu_seconds"] =
            static_cast<double>(combined.reduction_nanoseconds) * 1.0e-9;
        details["block_wall_seconds"] = seconds(blocks_started, blocks_finished);
        details["merge_seconds"] = seconds(merge_started, merge_finished);
        details["python_conversion_seconds"] =
            seconds(merge_finished, conversion_finished);
        details["total_seconds"] = seconds(total_started, conversion_finished);
        return details;
    }

    static void validate_tile_frame(
        const py::buffer_info &intensity,
        const py::buffer_info &variance,
        const py::buffer_info &mask,
        const std::size_t row_stop,
        const std::size_t column_stop
    ) {
        // The frames are read in place, so their layout is part of the
        // contract rather than something a cast quietly repairs: one row
        // stride has to serve every frame and every array.
        if (intensity.ndim != 2) {
            throw py::value_error(
                "each frame's intensity must be a two-dimensional "
                "(rows, columns) array"
            );
        }
        for (const py::buffer_info *other : {&variance, &mask}) {
            if (
                other->ndim != 2
                || other->shape[0] != intensity.shape[0]
                || other->shape[1] != intensity.shape[1]
            ) {
                throw py::value_error(
                    "each frame's variance and mask must match its intensity "
                    "shape"
                );
            }
        }
        if (
            static_cast<std::size_t>(intensity.shape[0]) < row_stop
            || static_cast<std::size_t>(intensity.shape[1]) < column_stop
        ) {
            throw py::value_error(
                "the tile rectangle reaches outside the frame"
            );
        }
    }

    static void validate_group_geometry(
        const py::buffer_info &rays,
        const py::buffer_info &start,
        const py::buffer_info &end,
        const std::size_t frames,
        const std::size_t rows,
        const std::size_t cols
    ) {
        if (
            rays.ndim != 3
            || static_cast<std::size_t>(rays.shape[0]) != rows + 1
            || static_cast<std::size_t>(rays.shape[1]) != cols + 1
            || rays.shape[2] != 3
        ) {
            throw py::value_error(
                "corner_rays must have shape (rows + 1, columns + 1, 3)"
            );
        }
        for (const py::buffer_info *angles : {&start, &end}) {
            if (
                angles->ndim != 2
                || static_cast<std::size_t>(angles->shape[0]) != frames
                || angles->shape[1] != 4
            ) {
                throw py::value_error(
                    "angles_start and angles_end must have shape (frames, 4)"
                );
            }
        }
    }

    static void validate_group_inputs(
        const py::buffer_info &intensity,
        const py::buffer_info &variance,
        const py::buffer_info &mask,
        const py::buffer_info &rays,
        const py::buffer_info &start,
        const py::buffer_info &end
    ) {
        if (intensity.ndim != 3) {
            throw py::value_error(
                "intensity must be a three-dimensional (frames, rows, columns) "
                "array"
            );
        }
        for (const py::buffer_info *other : {&variance, &mask}) {
            if (
                other->ndim != 3
                || other->shape[0] != intensity.shape[0]
                || other->shape[1] != intensity.shape[1]
                || other->shape[2] != intensity.shape[2]
            ) {
                throw py::value_error(
                    "variance and mask must match the intensity shape"
                );
            }
        }
        if (intensity.shape[0] < 1) {
            throw py::value_error("A frame group needs at least one frame");
        }
        if (
            rays.ndim != 3
            || rays.shape[0] != intensity.shape[1] + 1
            || rays.shape[1] != intensity.shape[2] + 1
            || rays.shape[2] != 3
        ) {
            throw py::value_error(
                "corner_rays must have shape (rows + 1, columns + 1, 3)"
            );
        }
        for (const py::buffer_info *angles : {&start, &end}) {
            if (
                angles->ndim != 2
                || angles->shape[0] != intensity.shape[0]
                || angles->shape[1] != 4
            ) {
                throw py::value_error(
                    "angles_start and angles_end must have shape (frames, 4)"
                );
            }
        }
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
        // All three axes are evaluated unconditionally and combined at the
        // end: the early-exit form this replaces put a data-dependent
        // branch between every axis, and the per-axis std::isfinite() call
        // it needed is subsumed by the range test, since a NaN compares
        // false against both bounds. Out-of-range and non-finite input is
        // rejected exactly as before.
        const std::array<double, 3> values{
            coordinate.x,
            coordinate.y,
            coordinate.z,
        };
        std::array<std::uint64_t, 3> index{};
        int inside = 1;
        for (int axis = 0; axis < 3; ++axis) {
            const double scaled =
                (values[axis] - grid_.minimum[axis]) / grid_.step[axis];
            const double floored = std::floor(scaled);
            const bool axis_inside =
                floored >= 0.0 && floored < grid_.shape_as_double[axis];
            inside &= static_cast<int>(axis_inside);
            // Comparing before converting also keeps the conversion itself
            // in range, which the old cast-then-test order did not.
            index[axis] = static_cast<std::uint64_t>(
                static_cast<std::int64_t>(axis_inside ? floored : 0.0)
            );
        }
        voxel = (index[0] << grid_.voxel_shift[0])
            | (index[1] << grid_.voxel_shift[1])
            | index[2];
        return inside != 0;
    }

    RecordKey record_key(const std::uint64_t voxel) const {
        // Both halves of this decomposition used to be 64-bit integer
        // divisions -- four to recover the axis indices from a row-major
        // flat id, six more to split each index into chunk and local
        // parts. voxel_id()'s bit-packed id makes the first four shifts
        // and masks, and a power-of-two chunk_shape (every default one)
        // makes the rest shifts and masks too.
        const std::uint64_t index_x =
            (voxel >> grid_.voxel_shift[0]) & grid_.voxel_mask[0];
        const std::uint64_t index_y =
            (voxel >> grid_.voxel_shift[1]) & grid_.voxel_mask[1];
        const std::uint64_t index_z = voxel & grid_.voxel_mask[2];
        std::uint64_t chunk_x;
        std::uint64_t chunk_y;
        std::uint64_t chunk_z;
        std::uint64_t local_x;
        std::uint64_t local_y;
        std::uint64_t local_z;
        if (grid_.chunk_power_of_two) {
            chunk_x = index_x >> grid_.chunk_shift[0];
            chunk_y = index_y >> grid_.chunk_shift[1];
            chunk_z = index_z >> grid_.chunk_shift[2];
            local_x = index_x & grid_.chunk_mask[0];
            local_y = index_y & grid_.chunk_mask[1];
            local_z = index_z & grid_.chunk_mask[2];
        } else {
            chunk_x = index_x / static_cast<std::uint64_t>(grid_.chunk_shape[0]);
            chunk_y = index_y / static_cast<std::uint64_t>(grid_.chunk_shape[1]);
            chunk_z = index_z / static_cast<std::uint64_t>(grid_.chunk_shape[2]);
            local_x = index_x % static_cast<std::uint64_t>(grid_.chunk_shape[0]);
            local_y = index_y % static_cast<std::uint64_t>(grid_.chunk_shape[1]);
            local_z = index_z % static_cast<std::uint64_t>(grid_.chunk_shape[2]);
        }
        RecordKey key{};
        // Arithmetic stays uint64 throughout (chunk_grid/chunk_shape are
        // uint64 already); only the final, per-axis-bounded result
        // narrows to uint32, safe once the construction-time
        // uint32-overflow validation on chunk_shape/grid.shape holds.
        key.chunk = static_cast<std::uint32_t>((
            chunk_x * grid_.chunk_grid[1] + chunk_y
        ) * grid_.chunk_grid[2] + chunk_z);
        key.local = static_cast<std::uint32_t>((
            local_x * static_cast<std::uint64_t>(grid_.chunk_shape[1]) + local_y
        ) * static_cast<std::uint64_t>(grid_.chunk_shape[2]) + local_z);
        return key;
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

    // One node of the 3x3x3 dyadic lattice a moving cell's subdivision
    // spans: (iu, iv, it) each in {0, 1, 2}, low corner to high corner.
    static std::size_t lattice_index(
        const std::size_t iu,
        const std::size_t iv,
        const std::size_t it
    ) {
        return (iu * 3 + iv) * 3 + it;
    }

    LatticeVoxel evaluate_lattice_voxel(
        const PixelRays &rays,
        const double u,
        const double v,
        const double t,
        const std::vector<CoordinateTransform> &transforms,
        BlockProfile *profile
    ) const {
        // t selects one of the rotations frame_rotations() precomputed over
        // the dyadic lattice, by the same rounding the removed shared-corner
        // cache used, so coordinates here are bit-for-bit what it produced.
        LatticeVoxel result;
        if (profile != nullptr) {
            ++profile->coordinate_evaluations;
        }
        const double rotation_scale =
            static_cast<double>(transforms.size() - 1);
        const std::size_t it = static_cast<std::size_t>(
            std::llround(t * rotation_scale)
        );
        result.valid = voxel_id(
            coordinate_at(rays, u, v, transforms[it]),
            result.voxel
        );
        return result;
    }

    void subdivide_moving(
        const PixelRays &rays,
        const std::vector<CoordinateTransform> &transforms,
        const double u0,
        const double u1,
        const double v0,
        const double v1,
        const double t0,
        const double t1,
        const double weight,
        const int depth,
        const LatticeVoxel (&corners)[8],
        std::vector<VoxelWeight> &weights,
        BlockProfile *profile
    ) const {
        // Corner c sits at (u1 if c&1 else u0, v1 if c&2 else v0,
        // t1 if c&4 else t0) -- the same numbering the explicit-stack form
        // used, and the same order children are visited in, so leaves are
        // emitted in an unchanged sequence and their later summation
        // associates identically.
        bool all_same = true;
        for (int corner = 1; corner < 8; ++corner) {
            if (
                corners[corner].valid != corners[0].valid
                || (corners[0].valid && corners[corner].voxel != corners[0].voxel)
            ) {
                all_same = false;
                break;
            }
        }
        if (all_same && corners[0].valid) {
            weights.push_back({corners[0].voxel, weight});
            return;
        }
        const double um = 0.5 * (u0 + u1);
        const double vm = 0.5 * (v0 + v1);
        const double tm = 0.5 * (t0 + t1);
        if (depth >= max_depth_) {
            // A cell whose corners disagree at the deepest level assigns its
            // whole weight at its centre. Note an entirely out-of-grid cell
            // reaches here too (all_same holds, but no corner is valid),
            // exactly as before.
            const LatticeVoxel centre = evaluate_lattice_voxel(
                rays, um, vm, tm, transforms, profile
            );
            if (centre.valid) {
                weights.push_back({centre.voxel, weight});
            }
            return;
        }
        // Eight children share one 27-node lattice, of which the eight
        // corners are already in hand: only the nineteen nodes with a
        // midpoint index are new. That is what replaces the dense
        // per-pixel cache -- the sharing lives in the recursion instead of
        // in a side array that has to be indexed, generation-stamped, and
        // sized for a depth it mostly does not use.
        LatticeVoxel lattice[27];
        const double us[3] = {u0, um, u1};
        const double vs[3] = {v0, vm, v1};
        const double ts[3] = {t0, tm, t1};
        for (std::size_t iu = 0; iu < 3; ++iu) {
            for (std::size_t iv = 0; iv < 3; ++iv) {
                for (std::size_t it = 0; it < 3; ++it) {
                    LatticeVoxel &node = lattice[lattice_index(iu, iv, it)];
                    if (iu != 1 && iv != 1 && it != 1) {
                        const int corner = (iu == 2 ? 1 : 0)
                            | (iv == 2 ? 2 : 0)
                            | (it == 2 ? 4 : 0);
                        node = corners[corner];
                    } else {
                        node = evaluate_lattice_voxel(
                            rays, us[iu], vs[iv], ts[it], transforms, profile
                        );
                    }
                }
            }
        }
        const double child_weight = weight / 8.0;
        for (int child = 0; child < 8; ++child) {
            const std::size_t base_u = (child & 1) != 0 ? 1 : 0;
            const std::size_t base_v = (child & 2) != 0 ? 1 : 0;
            const std::size_t base_t = (child & 4) != 0 ? 1 : 0;
            LatticeVoxel child_corners[8];
            for (int corner = 0; corner < 8; ++corner) {
                child_corners[corner] = lattice[
                    lattice_index(
                        base_u + ((corner & 1) != 0 ? 1 : 0),
                        base_v + ((corner & 2) != 0 ? 1 : 0),
                        base_t + ((corner & 4) != 0 ? 1 : 0)
                    )
                ];
            }
            subdivide_moving(
                rays,
                transforms,
                (child & 1) != 0 ? um : u0,
                (child & 1) != 0 ? u1 : um,
                (child & 2) != 0 ? vm : v0,
                (child & 2) != 0 ? v1 : vm,
                (child & 4) != 0 ? tm : t0,
                (child & 4) != 0 ? t1 : tm,
                child_weight,
                depth + 1,
                child_corners,
                weights,
                profile
            );
        }
    }

    LatticeVoxel evaluate_stationary_voxel(
        const PixelRays &rays,
        const double u,
        const double v,
        const CoordinateTransform &transform,
        BlockProfile *profile
    ) const {
        // The stationary counterpart of evaluate_lattice_voxel: t is pinned,
        // so the rotation is chosen once by the caller instead of being
        // recovered from t on every node. transforms[size / 2] and
        // llround(0.5 * (size - 1)) are the same index for every depth --
        // both are 2^max_depth -- so this is bit-for-bit what the cached
        // path evaluated.
        LatticeVoxel result;
        if (profile != nullptr) {
            ++profile->coordinate_evaluations;
        }
        result.valid = voxel_id(
            coordinate_at(rays, u, v, transform),
            result.voxel
        );
        return result;
    }

    void subdivide_stationary(
        const PixelRays &rays,
        const CoordinateTransform &transform,
        const double u0,
        const double u1,
        const double v0,
        const double v1,
        const double weight,
        const int depth,
        const LatticeVoxel (&corners)[4],
        std::vector<VoxelWeight> &weights,
        BlockProfile *profile
    ) const {
        bool all_same = true;
        for (int corner = 1; corner < 4; ++corner) {
            if (
                corners[corner].valid != corners[0].valid
                || (corners[0].valid && corners[corner].voxel != corners[0].voxel)
            ) {
                all_same = false;
                break;
            }
        }
        if (all_same && corners[0].valid) {
            weights.push_back({corners[0].voxel, weight});
            return;
        }
        const double um = 0.5 * (u0 + u1);
        const double vm = 0.5 * (v0 + v1);
        if (depth >= max_depth_) {
            const LatticeVoxel centre = evaluate_stationary_voxel(
                rays, um, vm, transform, profile
            );
            if (centre.valid) {
                weights.push_back({centre.voxel, weight});
            }
            return;
        }
        // Four children share one 3x3 lattice whose corners are already in
        // hand: five new nodes per split, exactly what
        // split_pixel_stationary_depth2 materializes by hand at depth two.
        LatticeVoxel lattice[9];
        const double us[3] = {u0, um, u1};
        const double vs[3] = {v0, vm, v1};
        for (std::size_t iu = 0; iu < 3; ++iu) {
            for (std::size_t iv = 0; iv < 3; ++iv) {
                LatticeVoxel &node = lattice[iu * 3 + iv];
                if (iu != 1 && iv != 1) {
                    const int corner = (iu == 2 ? 1 : 0) | (iv == 2 ? 2 : 0);
                    node = corners[corner];
                } else {
                    node = evaluate_stationary_voxel(
                        rays, us[iu], vs[iv], transform, profile
                    );
                }
            }
        }
        const double child_weight = weight / 4.0;
        for (int child = 0; child < 4; ++child) {
            const std::size_t base_u = (child & 1) != 0 ? 1 : 0;
            const std::size_t base_v = (child & 2) != 0 ? 1 : 0;
            LatticeVoxel child_corners[4];
            for (int corner = 0; corner < 4; ++corner) {
                child_corners[corner] = lattice[
                    (base_u + ((corner & 1) != 0 ? 1 : 0)) * 3
                    + base_v + ((corner & 2) != 0 ? 1 : 0)
                ];
            }
            subdivide_stationary(
                rays,
                transform,
                (child & 1) != 0 ? um : u0,
                (child & 1) != 0 ? u1 : um,
                (child & 2) != 0 ? vm : v0,
                (child & 2) != 0 ? v1 : vm,
                child_weight,
                depth + 1,
                child_corners,
                weights,
                profile
            );
        }
    }

    void split_pixel_stationary(
        const PixelRays &rays,
        const CoordinateTransform &transform,
        std::vector<VoxelWeight> &weights,
        BlockProfile *profile
    ) const {
        LatticeVoxel corners[4];
        for (int corner = 0; corner < 4; ++corner) {
            corners[corner] = evaluate_stationary_voxel(
                rays,
                (corner & 1) != 0 ? 1.0 : 0.0,
                (corner & 2) != 0 ? 1.0 : 0.0,
                transform,
                profile
            );
        }
        subdivide_stationary(
            rays,
            transform,
            0.0,
            1.0,
            0.0,
            1.0,
            1.0,
            0,
            corners,
            weights,
            profile
        );
    }

    void split_pixel_moving(
        const PixelRays &rays,
        const std::vector<CoordinateTransform> &transforms,
        std::vector<VoxelWeight> &weights,
        BlockProfile *profile
    ) const {
        LatticeVoxel corners[8];
        for (int corner = 0; corner < 8; ++corner) {
            corners[corner] = evaluate_lattice_voxel(
                rays,
                (corner & 1) != 0 ? 1.0 : 0.0,
                (corner & 2) != 0 ? 1.0 : 0.0,
                (corner & 4) != 0 ? 1.0 : 0.0,
                transforms,
                profile
            );
        }
        subdivide_moving(
            rays,
            transforms,
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            1.0,
            1.0,
            0,
            corners,
            weights,
            profile
        );
    }

    std::vector<Record> accumulate_brick(
        const BrickExtent &brick,
        const std::size_t frames,
        const std::size_t columns,
        const GroupPixels &pixels,
        const double *rays,
        const std::vector<std::vector<CoordinateTransform>> &transforms,
        const std::vector<char> &stationary,
        std::vector<Vec3> &centre_rays,
        std::vector<VoxelWeight> &weights,
        BlockProfile *profile,
        std::pmr::monotonic_buffer_resource &arena
    ) const {
        arena.release();
        std::pmr::map<RecordKey, RecordAccum> tree(&arena);
        std::uint64_t unreduced_count = 0;
        RecordKey cached_key{};
        RecordAccum *cached_accumulator = nullptr;
        const auto accumulator_for = [&tree, &cached_key, &cached_accumulator](
            const RecordKey key
        ) -> RecordAccum & {
            if (cached_accumulator != nullptr && key == cached_key) {
                return *cached_accumulator;
            }
            RecordAccum &accumulator = tree[key];
            cached_key = key;
            cached_accumulator = &accumulator;
            return accumulator;
        };
        const std::size_t brick_rows = brick.row_end - brick.row_begin;
        const std::size_t brick_columns = brick.column_end - brick.column_begin;
        const std::size_t row_stride = pixels.row_stride;

        if (max_depth_ == 0) {
            // A pixel's centre ray is pure detector geometry -- the frame's
            // rotation enters only through the transform applied after it --
            // so it is the same for every frame in the group and is computed
            // once per brick instead of once per (pixel, frame). At 32x32
            // that buffer is 24 KiB.
            centre_rays.resize(brick_rows * brick_columns);
            for (std::size_t row = brick.row_begin; row < brick.row_end; ++row) {
                for (
                    std::size_t column = brick.column_begin;
                    column < brick.column_end;
                    ++column
                ) {
                    centre_rays[
                        (row - brick.row_begin) * brick_columns
                        + (column - brick.column_begin)
                    ] = ray_at(pixel_rays(row, column, columns, rays), 0.5, 0.5);
                }
            }
        }

        const auto mapping_started = std::chrono::steady_clock::now();
        for (std::size_t frame = 0; frame < frames; ++frame) {
            const std::vector<CoordinateTransform> &frame_transforms =
                transforms[frame];
            const CoordinateTransform &centre_transform =
                frame_transforms[frame_transforms.size() / 2];
            const double *const intensity = pixels.intensity[frame];
            const double *const variance = pixels.variance[frame];
            const bool *const mask = pixels.mask[frame];
            for (std::size_t row = brick.row_begin; row < brick.row_end; ++row) {
                const std::size_t row_offset = row * row_stride;
                for (
                    std::size_t column = brick.column_begin;
                    column < brick.column_end;
                    ++column
                ) {
                    const std::size_t flat = row_offset + column;
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
                    if (max_depth_ == 0) {
                        if (profile != nullptr) {
                            ++profile->coordinate_evaluations;
                        }
                        std::uint64_t voxel = 0;
                        if (
                            voxel_id(
                                coordinate_from_unit_ray(
                                    centre_rays[
                                        (row - brick.row_begin) * brick_columns
                                        + (column - brick.column_begin)
                                    ],
                                    centre_transform
                                ),
                                voxel
                            )
                        ) {
                            if (profile != nullptr) {
                                ++profile->voxel_weights;
                                profile->maximum_weights_per_pixel = std::max(
                                    profile->maximum_weights_per_pixel,
                                    static_cast<std::uint64_t>(1)
                                );
                            }
                            RecordAccum &accumulator =
                                accumulator_for(record_key(voxel));
                            accumulator.weighted_intensity += intensity[flat];
                            accumulator.weighted_variance += variance[flat];
                            accumulator.weight += 1.0;
                            accumulator.contributors += 1;
                            ++unreduced_count;
                        }
                        continue;
                    }
                    const PixelRays prepared_rays =
                        pixel_rays(row, column, columns, rays);
                    weights.clear();
                    if (stationary[frame] != 0 && max_depth_ == 2) {
                        split_pixel_stationary_depth2(
                            prepared_rays, centre_transform, weights, profile
                        );
                    } else if (stationary[frame] != 0) {
                        split_pixel_stationary(
                            prepared_rays, centre_transform, weights, profile
                        );
                    } else {
                        split_pixel_moving(
                            prepared_rays, frame_transforms, weights, profile
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
                        RecordAccum &accumulator =
                            accumulator_for(record_key(voxel));
                        accumulator.weighted_intensity +=
                            weight * intensity[flat];
                        accumulator.weighted_variance +=
                            weight * weight * variance[flat];
                        accumulator.weight += weight;
                        accumulator.contributors += 1;
                        ++unreduced_count;
                    }
                }
            }
        }
        const auto mapping_finished = std::chrono::steady_clock::now();
        if (profile != nullptr) {
            profile->unreduced_records = unreduced_count;
            profile->mapping_nanoseconds = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    mapping_finished - mapping_started
                ).count()
            );
        }
        std::vector<Record> records;
        records.reserve(tree.size());
        for (const auto &entry : tree) {
            records.push_back({
                entry.first,
                entry.second.weighted_intensity,
                entry.second.weighted_variance,
                entry.second.weight,
                entry.second.contributors,
            });
        }
        if (profile != nullptr) {
            profile->reduced_records = records.size();
            profile->reduction_nanoseconds = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - mapping_finished
                ).count()
            );
        }
        return records;
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
        BlockProfile *profile,
        std::pmr::monotonic_buffer_resource &arena
    ) const {
        (void)rows;
        // Reset the caller-owned arena to the start of its (already
        // allocated once, reused across every block this worker thread
        // processes) buffer -- no allocation or deallocation here, just a
        // cursor reset. This is what avoids the per-block malloc/free
        // churn a fresh arena-per-call would cause under heavy thread
        // contention.
        arena.release();
        std::pmr::map<RecordKey, RecordAccum> tree(&arena);
        std::vector<VoxelWeight> weights;
        // Subdivision now carries its cell state on the recursion stack --
        // nine LatticeVoxel per level stationary, twenty-seven moving, so
        // at most about four kilobytes at the deepest setting -- so the
        // only scratch a block still needs is the leaf list.
        weights.reserve(
            static_cast<std::size_t>(1)
            << std::min((stationary ? 2 : 3) * max_depth_, 9)
        );
        // Counts per-pixel-merged (voxel, weight) accumulation events --
        // the same quantity the old push_back-one-Record-per-event code
        // counted via records.size() before its later block-end reduce.
        // tree.size() alone cannot recover this: it's already the
        // post-dedup unique-key count by the time any pixel has been
        // processed, since accumulation happens directly into the map.
        std::uint64_t unreduced_count = 0;
        // Walked forward instead of recovered per pixel: the flat index
        // divided and taken modulo the row length was an integer division
        // and remainder on every pixel, and a block is a contiguous run.
        std::size_t row_cursor = begin / columns;
        std::size_t column_cursor = begin % columns;
        // Neighbouring pixels overwhelmingly land in the same voxel (only
        // about a third of adjacent pairs cross a voxel boundary on real
        // data), so remembering the last key's accumulator turns most
        // accumulations into two uint32 compares instead of a tree
        // descent. Map nodes never move, so the pointer stays valid for
        // the whole block.
        RecordKey cached_key{};
        RecordAccum *cached_accumulator = nullptr;
        const auto accumulator_for = [&tree, &cached_key, &cached_accumulator](
            const RecordKey key
        ) -> RecordAccum & {
            if (cached_accumulator != nullptr && key == cached_key) {
                return *cached_accumulator;
            }
            RecordAccum &accumulator = tree[key];
            cached_key = key;
            cached_accumulator = &accumulator;
            return accumulator;
        };
        const auto mapping_started = std::chrono::steady_clock::now();
        for (std::size_t flat = begin; flat < end; ++flat) {
            const std::size_t row = row_cursor;
            const std::size_t column = column_cursor;
            // Advanced here rather than at the bottom, so the skip paths
            // below cannot step over it.
            if (++column_cursor == columns) {
                column_cursor = 0;
                ++row_cursor;
            }
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
            const PixelRays prepared_rays = pixel_rays(
                row,
                column,
                columns,
                rays
            );
            if (max_depth_ == 0) {
                // One pixel centre, so exactly one (voxel, weight=1) pair:
                // the scratch vector, its sort, and the duplicate-merging
                // walk below all have nothing to do at this depth, and
                // multiplying by a weight of exactly one is exact.
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
                    if (profile != nullptr) {
                        ++profile->voxel_weights;
                        profile->maximum_weights_per_pixel = std::max(
                            profile->maximum_weights_per_pixel,
                            static_cast<std::uint64_t>(1)
                        );
                    }
                    RecordAccum &accumulator = accumulator_for(record_key(voxel));
                    accumulator.weighted_intensity += intensity[flat];
                    accumulator.weighted_variance += variance[flat];
                    accumulator.weight += 1.0;
                    accumulator.contributors += 1;
                    ++unreduced_count;
                }
                continue;
            }
            weights.clear();
            if (stationary && max_depth_ == 2) {
                // Kept in preference to the general recursion, which measured
                // 10% slower here: one flat 5x5 lattice per pixel shares
                // nodes between sibling subtrees, which corner-passing gives
                // up (22.02 against 20.82 evaluations per pixel). That
                // sharing is only affordable because a stationary depth-two
                // lattice is 25 nodes; see the design note on generalizing
                // it.
                split_pixel_stationary_depth2(
                    prepared_rays,
                    transforms[transforms.size() / 2],
                    weights,
                    profile
                );
            } else if (stationary) {
                split_pixel_stationary(
                    prepared_rays,
                    transforms[transforms.size() / 2],
                    weights,
                    profile
                );
            } else {
                split_pixel_moving(
                    prepared_rays,
                    transforms,
                    weights,
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
                // Accumulate directly into the tree instead of
                // push_back-ing a raw Record for a later sort+merge pass
                // -- the map already groups by key (and, iterated, is
                // already in sorted key order) as pixels are processed,
                // so no separate block-end reduce step is needed at all.
                RecordAccum &accumulator = accumulator_for(record_key(voxel));
                accumulator.weighted_intensity += weight * intensity[flat];
                accumulator.weighted_variance += weight * weight * variance[flat];
                accumulator.weight += weight;
                accumulator.contributors += 1;
                ++unreduced_count;
            }
        }
        const auto mapping_finished = std::chrono::steady_clock::now();
        if (profile != nullptr) {
            profile->unreduced_records = unreduced_count;
            profile->mapping_nanoseconds = static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    mapping_finished - mapping_started
                ).count()
            );
        }
        const auto reduction_started = mapping_finished;
        std::vector<Record> records;
        records.reserve(tree.size());
        for (const auto &entry : tree) {
            records.push_back({
                entry.first,
                entry.second.weighted_intensity,
                entry.second.weighted_variance,
                entry.second.weight,
                entry.second.contributors,
            });
        }
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

    // Loser tree (tournament tree) k-way merge of block_results -- each
    // block's output is already sorted and deduped (accumulate_block's
    // pmr::map reduce), so this replaces reduce_records' blind
    // concatenate-then-resort with a single O(N log k) pass that
    // exploits that pre-sortedness instead of ignoring it. One
    // comparison per tree level on replay, versus a priority-queue
    // heap's sift-down (~two comparisons per level).
    struct SortedRun {
        const Record *begin;
        const Record *end;
    };

    class LoserTree {
    public:
        explicit LoserTree(const std::vector<SortedRun> &runs) : runs_(runs) {
            std::size_t real_runs = runs.size();
            leaf_count_ = 1;
            while (leaf_count_ < real_runs) {
                leaf_count_ <<= 1;
            }
            if (leaf_count_ < 2) {
                leaf_count_ = 2;
            }
            run_index_.assign(leaf_count_, 0);
            current_key_.assign(leaf_count_, sentinel_key());
            for (std::size_t i = 0; i < real_runs; ++i) {
                if (runs[i].begin != runs[i].end) {
                    current_key_[i] = runs[i].begin[0].key;
                }
            }
            loser_.assign(leaf_count_, 0);
            build();
        }

        static RecordKey sentinel_key() {
            return RecordKey{
                std::numeric_limits<std::uint32_t>::max(),
                std::numeric_limits<std::uint32_t>::max(),
            };
        }

        bool exhausted() const {
            return current_key_[loser_[0]] == sentinel_key();
        }

        const Record &winner_record() const {
            return runs_[loser_[0]].begin[run_index_[loser_[0]]];
        }

        void advance_winner() {
            const std::uint32_t slot = loser_[0];
            ++run_index_[slot];
            const SortedRun &run = runs_[slot];
            current_key_[slot] = (run.begin + run_index_[slot] < run.end)
                ? run.begin[run_index_[slot]].key
                : sentinel_key();
            replay(slot);
        }

    private:
        void build() {
            std::vector<std::uint32_t> winner(2 * leaf_count_);
            for (std::uint32_t i = 0; i < leaf_count_; ++i) {
                winner[leaf_count_ + i] = i;
            }
            for (std::uint32_t node = leaf_count_ - 1; node >= 1; --node) {
                const std::uint32_t left = winner[2 * node];
                const std::uint32_t right = winner[2 * node + 1];
                if (
                    current_key_[left] < current_key_[right]
                    || current_key_[left] == current_key_[right]
                ) {
                    winner[node] = left;
                    loser_[node] = right;
                } else {
                    winner[node] = right;
                    loser_[node] = left;
                }
                if (node == 1) {
                    break;
                }
            }
            loser_[0] = winner[1];
        }

        void replay(const std::uint32_t slot) {
            std::uint32_t candidate = slot;
            std::uint32_t node = (leaf_count_ + slot) / 2;
            while (node >= 1) {
                if (!(current_key_[candidate] < current_key_[loser_[node]])) {
                    std::swap(candidate, loser_[node]);
                }
                if (node == 1) {
                    break;
                }
                node /= 2;
            }
            loser_[0] = candidate;
        }

        const std::vector<SortedRun> &runs_;
        std::uint32_t leaf_count_;
        std::vector<std::uint32_t> run_index_;
        std::vector<RecordKey> current_key_;
        std::vector<std::uint32_t> loser_;
    };

    static std::vector<Record> merge_sorted_blocks(
        const std::vector<std::vector<Record>> &block_results
    ) {
        std::vector<SortedRun> runs;
        runs.reserve(block_results.size());
        std::size_t total = 0;
        for (const auto &block : block_results) {
            runs.push_back({block.data(), block.data() + block.size()});
            total += block.size();
        }
        std::vector<Record> merged;
        merged.reserve(total);
        LoserTree tree(runs);
        while (!tree.exhausted()) {
            Record combined = tree.winner_record();
            const RecordKey winner_key = combined.key;
            tree.advance_winner();
            while (!tree.exhausted() && tree.winner_record().key == winner_key) {
                const Record &duplicate = tree.winner_record();
                combined.weighted_intensity += duplicate.weighted_intensity;
                combined.weighted_variance += duplicate.weighted_variance;
                combined.weight += duplicate.weight;
                combined.contributors += duplicate.contributors;
                tree.advance_winner();
            }
            merged.push_back(combined);
        }
        return merged;
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
            "accumulate_group",
            &ReconstructionKernel::accumulate_group,
            py::arg("intensity"),
            py::arg("variance"),
            py::arg("mask"),
            py::arg("corner_rays"),
            py::arg("angles_start"),
            py::arg("angles_end"),
            py::arg("profile") = false
        )
        .def(
            "accumulate_group_tile",
            &ReconstructionKernel::accumulate_group_tile,
            py::arg("intensity"),
            py::arg("variance"),
            py::arg("mask"),
            py::arg("corner_rays"),
            py::arg("angles_start"),
            py::arg("angles_end"),
            py::arg("row_start"),
            py::arg("row_stop"),
            py::arg("column_start"),
            py::arg("column_stop"),
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

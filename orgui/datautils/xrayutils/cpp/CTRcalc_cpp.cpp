#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstring>
#include <cstdint>
#include <list>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <xxhash.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

using complex128 = std::complex<double>;

using Array1D = py::array_t<double, py::array::c_style>;
using Array2D = py::array_t<double, py::array::c_style>;
using Array3D = py::array_t<double, py::array::c_style>;
using ComplexArray1D = py::array_t<complex128, py::array::c_style>;
using ComplexArray2D = py::array_t<complex128, py::array::c_style>;
using IntArray1D = py::array_t<std::int64_t, py::array::c_style>;

constexpr std::size_t form_factor_cache_default_budget = 256ULL * 1024 * 1024;

struct Hash128 {
    std::uint64_t low;
    std::uint64_t high;
};

Hash128 hash_values(const double *values, const std::size_t count) {
    const XXH128_hash_t hash = XXH3_128bits(values, count * sizeof(double));
    return {hash.low64, hash.high64};
}

bool same_hash(const Hash128 &left, const Hash128 &right) {
    return left.low == right.low && left.high == right.high;
}

struct CachedSpecies {
    std::array<double, 13> row;
    Hash128 hash;
    std::shared_ptr<const std::vector<double>> values;
    // Valid once inserted into FormFactorCache::lru_; lets touch() splice
    // this entry to the LRU tail in O(1) instead of scanning lru_.
    std::list<std::shared_ptr<CachedSpecies>>::iterator lru_it;
};

struct CachedQGrid {
    std::shared_ptr<const std::vector<double>> q2;
    Hash128 hash;
    std::vector<std::shared_ptr<CachedSpecies>> species;
};

class FormFactorCache {
public:
    std::shared_ptr<const std::vector<double>> find(
        const std::vector<double> &q2,
        const std::array<double, 13> &row
    ) {
        const Hash128 q_hash = hash_values(q2.data(), q2.size());
        const Hash128 row_hash = hash_values(row.data(), row.size());
        std::lock_guard<std::mutex> lock(mutex_);
        auto grid = find_grid(q2, q_hash);
        if (grid == grids_.end()) {
            ++misses_;
            return nullptr;
        }
        for (const auto &species : (*grid)->species) {
            if (same_hash(species->hash, row_hash)
                && std::memcmp(species->row.data(), row.data(), sizeof(row)) == 0) {
                touch(species);
                ++hits_;
                return species->values;
            }
        }
        ++misses_;
        return nullptr;
    }

    std::shared_ptr<const std::vector<double>> insert_or_get(
        const std::vector<double> &q2,
        const std::array<double, 13> &row,
        std::vector<double> values
    ) {
        const Hash128 q_hash = hash_values(q2.data(), q2.size());
        const Hash128 row_hash = hash_values(row.data(), row.size());
        std::lock_guard<std::mutex> lock(mutex_);
        auto grid = find_grid(q2, q_hash);
        if (grid != grids_.end()) {
            for (const auto &species : (*grid)->species) {
                if (same_hash(species->hash, row_hash)
                    && std::memcmp(species->row.data(), row.data(), sizeof(row)) == 0) {
                    touch(species);
                    ++hits_;
                    return species->values;
                }
            }
        }

        const std::size_t values_bytes = values.size() * sizeof(double);
        const std::size_t q_bytes = grid == grids_.end() ? q2.size() * sizeof(double) : 0;
        if (values_bytes + q_bytes > budget_bytes_) {
            return std::make_shared<const std::vector<double>>(std::move(values));
        }
        if (grid == grids_.end()) {
            auto new_grid = std::make_shared<CachedQGrid>();
            new_grid->q2 = std::make_shared<const std::vector<double>>(q2);
            new_grid->hash = q_hash;
            grids_.push_back(new_grid);
            grid = std::prev(grids_.end());
            grid_index_[q_hash.low].push_back(grid);
            resident_bytes_ += q_bytes;
        }
        auto species = std::make_shared<CachedSpecies>();
        species->row = row;
        species->hash = row_hash;
        species->values = std::make_shared<const std::vector<double>>(std::move(values));
        (*grid)->species.push_back(species);
        lru_.push_back(species);
        species->lru_it = std::prev(lru_.end());
        resident_bytes_ += values_bytes;
        evict_to_budget();
        return species->values;
    }

    void clear() {
        std::lock_guard<std::mutex> lock(mutex_);
        grids_.clear();
        grid_index_.clear();
        lru_.clear();
        resident_bytes_ = 0;
    }

    void reset_statistics() {
        std::lock_guard<std::mutex> lock(mutex_);
        hits_ = 0;
        misses_ = 0;
        evictions_ = 0;
    }

    void set_budget(const std::size_t bytes) {
        std::lock_guard<std::mutex> lock(mutex_);
        budget_bytes_ = bytes;
        evict_to_budget();
    }

    py::dict statistics() {
        std::lock_guard<std::mutex> lock(mutex_);
        py::dict result;
        result["hits"] = hits_;
        result["misses"] = misses_;
        result["evictions"] = evictions_;
        result["resident_bytes"] = resident_bytes_;
        result["budget_bytes"] = budget_bytes_;
        result["q_grids"] = grids_.size();
        result["species_entries"] = lru_.size();
        return result;
    }

private:
    using GridList = std::list<std::shared_ptr<CachedQGrid>>;

    // Bucketed by the low 64 bits of the Q-grid's hash; find_grid still
    // verifies the full 128-bit hash and the raw Q2 values within a bucket,
    // so a same-bucket collision can only cost extra comparisons, never a
    // wrong result.
    GridList::iterator find_grid(const std::vector<double> &q2, const Hash128 &hash) {
        const auto bucket = grid_index_.find(hash.low);
        if (bucket == grid_index_.end()) {
            return grids_.end();
        }
        for (const auto &grid : bucket->second) {
            const auto &stored = *(*grid)->q2;
            if (same_hash((*grid)->hash, hash) && stored.size() == q2.size()
                && std::memcmp(stored.data(), q2.data(), q2.size() * sizeof(double)) == 0) {
                return grid;
            }
        }
        return grids_.end();
    }

    void touch(const std::shared_ptr<CachedSpecies> &species) {
        lru_.splice(lru_.end(), lru_, species->lru_it);
    }

    void evict_to_budget() {
        while (resident_bytes_ > budget_bytes_ && !lru_.empty()) {
            const auto species = lru_.front();
            lru_.pop_front();
            for (auto grid = grids_.begin(); grid != grids_.end(); ++grid) {
                auto &entries = (*grid)->species;
                for (auto entry = entries.begin(); entry != entries.end(); ++entry) {
                    if (*entry == species) {
                        resident_bytes_ -= (*entry)->values->size() * sizeof(double);
                        entries.erase(entry);
                        ++evictions_;
                        if (entries.empty()) {
                            resident_bytes_ -= (*grid)->q2->size() * sizeof(double);
                            auto &bucket = grid_index_[(*grid)->hash.low];
                            bucket.erase(
                                std::remove(bucket.begin(), bucket.end(), grid),
                                bucket.end()
                            );
                            if (bucket.empty()) {
                                grid_index_.erase((*grid)->hash.low);
                            }
                            grids_.erase(grid);
                        }
                        goto evicted;
                    }
                }
            }
        evicted:;
        }
    }

    std::mutex mutex_;
    GridList grids_;
    std::unordered_map<std::uint64_t, std::vector<GridList::iterator>> grid_index_;
    std::list<std::shared_ptr<CachedSpecies>> lru_;
    std::size_t budget_bytes_ = form_factor_cache_default_budget;
    std::size_t resident_bytes_ = 0;
    std::uint64_t hits_ = 0;
    std::uint64_t misses_ = 0;
    std::uint64_t evictions_ = 0;
};

FormFactorCache &form_factor_cache() {
    static FormFactorCache cache;
    return cache;
}

py::dict form_factor_cache_stats() {
    return form_factor_cache().statistics();
}

void reset_form_factor_cache_stats() {
    form_factor_cache().reset_statistics();
}

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

std::unique_ptr<Matrix3[]> effective_domain_matrices(
    const py::buffer_info &coherent_domain_matrix,
    const py::buffer_info &r_mat_info,
    const py::buffer_info &r_mat_inv_info
) {
    const auto n_domains = coherent_domain_matrix.shape[0];
    const double *domain = ptr(coherent_domain_matrix);
    const double *r_mat = ptr(r_mat_info);
    const double *r_mat_inv = ptr(r_mat_inv_info);
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
                                coherent_domain_matrix,
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

std::unique_ptr<Matrix3[]> effective_domain_matrices(const Inputs &inputs) {
    return effective_domain_matrices(
        inputs.coherent_domain_matrix,
        inputs.r_mat,
        inputs.r_mat_inv
    );
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
            const double denominator_phase = -two_pi * ll;
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
    const Inputs inputs = request_inputs(
        h, k, l, basis, f_factors, ref_hkl_transform, b_mat, r_mat,
        r_mat_inv, coherent_domain_matrix, coherent_domain_occupancy
    );
    const auto n_points = inputs.h.shape[0];
    const auto n_atoms = inputs.basis.shape[0];
    const auto n_domains = inputs.coherent_domain_matrix.shape[0];
    auto result = py::array_t<complex128>(n_points);
    auto *out = static_cast<complex128 *>(result.request().ptr);

    const double *h_data = ptr(inputs.h);
    const double *k_data = ptr(inputs.k);
    const double *l_data = ptr(inputs.l);
    const double *basis_data = ptr(inputs.basis);
    const double *ff_data = ptr(inputs.f_factors);
    const double *ref_data = ptr(inputs.ref_hkl_transform);
    const double *b_data = ptr(inputs.b_mat);
    const double *domain_data = ptr(inputs.coherent_domain_matrix);
    const double *domain_occupancy = ptr(inputs.coherent_domain_occupancy);
    const auto domain_matrices = effective_domain_matrices(inputs);

    struct PointGeometry {
        double h;
        double k;
        double l;
        double q_para2;
        double q_perp2;
    };
    struct AtomDomainGeometry {
        double x;
        double y;
        double z;
        double weight;
    };
    struct SpeciesGroup {
        std::array<double, 13> row;
        std::vector<py::ssize_t> atoms;
    };

    std::vector<PointGeometry> points(n_points);
    std::vector<double> q2(n_points);
    std::vector<AtomDomainGeometry> atom_domains(n_atoms * n_domains);
    std::vector<SpeciesGroup> species;
    species.reserve(n_atoms);
    for (py::ssize_t i = 0; i < n_atoms; ++i) {
        std::array<double, 13> row;
        std::memcpy(row.data(), ff_data + i * inputs.f_factors.shape[1], sizeof(row));
        auto group = species.end();
        for (auto candidate = species.begin(); candidate != species.end(); ++candidate) {
            if (std::memcmp(candidate->row.data(), row.data(), sizeof(row)) == 0) {
                group = candidate;
                break;
            }
        }
        if (group == species.end()) {
            species.push_back({row, {}});
            group = std::prev(species.end());
        }
        group->atoms.push_back(i);
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    constexpr double dw_denominator = 16.0 * pi * pi;
    py::gil_scoped_release release;

    // These snapshots are deliberately made before hashing: exact numerical
    // reuse is defined by the computed float64 Q^2 grid, not source h/k/l.
    for (py::ssize_t p = 0; p < n_points; ++p) {
        const double h_in = h_data[p];
        const double k_in = k_data[p];
        const double l_in = l_data[p];
        const double hh = ref_data[0] * h_in + ref_data[1] * k_in + ref_data[2] * l_in;
        const double kk = ref_data[3] * h_in + ref_data[4] * k_in + ref_data[5] * l_in;
        const double ll = ref_data[6] * h_in + ref_data[7] * k_in + ref_data[8] * l_in;
        const double qx = b_data[0] * hh + b_data[1] * kk + b_data[2] * ll;
        const double qy = b_data[3] * hh + b_data[4] * kk + b_data[5] * ll;
        const double qz = b_data[6] * hh + b_data[7] * kk + b_data[8] * ll;
        points[p] = {hh, kk, ll, qx * qx + qy * qy, qz * qz};
        q2[p] = points[p].q_para2 + points[p].q_perp2;
    }
    for (py::ssize_t i = 0; i < n_atoms; ++i) {
        const double x = get2(basis_data, inputs.basis, i, 1);
        const double y = get2(basis_data, inputs.basis, i, 2);
        const double z = get2(basis_data, inputs.basis, i, 3);
        for (py::ssize_t d = 0; d < n_domains; ++d) {
            const Matrix3 &mat = domain_matrices[d];
            atom_domains[i * n_domains + d] = {
                mat.value[0][0] * x + mat.value[0][1] * y + mat.value[0][2] * z
                    + get3(domain_data, inputs.coherent_domain_matrix, d, 0, 3),
                mat.value[1][0] * x + mat.value[1][1] * y + mat.value[1][2] * z
                    + get3(domain_data, inputs.coherent_domain_matrix, d, 1, 3),
                mat.value[2][0] * x + mat.value[2][1] * y + mat.value[2][2] * z
                    + get3(domain_data, inputs.coherent_domain_matrix, d, 2, 3),
                domain_occupancy[d],
            };
        }
    }

    std::vector<std::shared_ptr<const std::vector<double>>> factors;
    factors.reserve(species.size());
    for (const auto &group : species) {
        auto values = form_factor_cache().find(q2, group.row);
        if (!values) {
            std::vector<double> calculated(n_points);
            for (py::ssize_t p = 0; p < n_points; ++p) {
                double value = group.row[10] + group.row[11];
                for (int term = 0; term < 5; ++term) {
                    value += group.row[term] * std::exp(-group.row[term + 5] * q2[p]);
                }
                calculated[p] = value;
            }
            values = form_factor_cache().insert_or_get(q2, group.row, std::move(calculated));
        }
        factors.push_back(std::move(values));
    }

    for (py::ssize_t p = 0; p < n_points; ++p) {
        double amplitude_real = 0.0;
        double amplitude_imag = 0.0;
        for (std::size_t g = 0; g < species.size(); ++g) {
            const double f_real = (*factors[g])[p];
            const double f_imag = species[g].row[12];
            for (const py::ssize_t i : species[g].atoms) {
                const double disorder = std::exp(-(
                    get2(basis_data, inputs.basis, i, 4) * points[p].q_para2
                    + get2(basis_data, inputs.basis, i, 5) * points[p].q_perp2
                ) / dw_denominator) * get2(basis_data, inputs.basis, i, 6);
                for (py::ssize_t d = 0; d < n_domains; ++d) {
                    const auto &atom = atom_domains[i * n_domains + d];
                    const double phase = two_pi * (
                        points[p].h * atom.x + points[p].k * atom.y + points[p].l * atom.z
                    );
                    const double domain_real = atom.weight * std::cos(phase);
                    const double domain_imag = atom.weight * std::sin(phase);
                    amplitude_real += disorder * (
                        domain_real * f_real - domain_imag * f_imag
                    );
                    amplitude_imag += disorder * (
                        domain_real * f_imag + domain_imag * f_real
                    );
                }
            }
        }
        out[p] = complex128{amplitude_real, amplitude_imag};
    }
    return result;
}

py::array_t<complex128> unitcell_zdensity_g(
    const Array1D &z,
    const double h,
    const double k,
    const double q_para2,
    const double c,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy,
    const double uc_area
) {
    const auto z_info = z.request();
    const auto basis_info = basis.request();
    const auto ff_info = f_factors.request();
    const auto r_info = r_mat.request();
    const auto r_inv_info = r_mat_inv.request();
    const auto domain_info = coherent_domain_matrix.request();
    const auto occupancy_info = coherent_domain_occupancy.request();

    require_shape(z_info, 1, "z");
    require_shape(basis_info, 2, "basis");
    require_shape(ff_info, 2, "f_factors");
    require_shape(r_info, 2, "R_mat");
    require_shape(r_inv_info, 2, "R_mat_inv");
    require_shape(domain_info, 3, "coherentDomainMatrix");
    require_shape(occupancy_info, 1, "coherentDomainOccupancy");
    if (basis_info.shape[1] < 7) {
        throw py::value_error("basis must have at least 7 columns");
    }
    if (ff_info.shape[0] != basis_info.shape[0] || ff_info.shape[1] < 13) {
        throw py::value_error(
            "f_factors must have one row per basis atom and at least 13 columns"
        );
    }
    if (r_info.shape[0] != 3 || r_info.shape[1] != 3
        || r_inv_info.shape[0] != 3 || r_inv_info.shape[1] != 3) {
        throw py::value_error("transform matrices must be shape (3, 3)");
    }
    if (domain_info.shape[1] != 3 || domain_info.shape[2] != 4) {
        throw py::value_error("coherentDomainMatrix must have shape (N, 3, 4)");
    }
    if (occupancy_info.shape[0] != domain_info.shape[0]) {
        throw py::value_error(
            "coherentDomainOccupancy length must match coherentDomainMatrix"
        );
    }

    const auto domain_matrices = effective_domain_matrices(
        domain_info,
        r_info,
        r_inv_info
    );
    const auto n_z = z_info.shape[0];
    const auto n_atoms = basis_info.shape[0];
    const auto n_domains = domain_info.shape[0];
    const auto *z_data = static_cast<const double *>(z_info.ptr);
    const auto *basis_data = static_cast<const double *>(basis_info.ptr);
    const auto *ff_data = static_cast<const double *>(ff_info.ptr);
    const auto *domain_data = static_cast<const double *>(domain_info.ptr);
    const auto *occupancy_data = static_cast<const double *>(occupancy_info.ptr);
    auto result = py::array_t<complex128>(n_z);
    auto *out = static_cast<complex128 *>(result.request().ptr);

    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double two_pi = 2.0 * pi;
    struct DensityAtomDomain {
        double delta_z2;
        double delta_para_q;
        double center_z;
        double base_real;
        double base_imag;
        std::array<double, 5> term_coefficients;
        std::array<double, 5> exp_dpara_q;
        std::array<double, 5> exp_dz;
        double phase_real;
        double phase_imag;
    };

    // zDensity_G evaluates the same atom/domain geometry for every profile
    // sample. Keep the z-independent work outside the profile loop; the
    // original Python implementation obtains the same benefit by operating
    // on a whole z vector for each atom/domain pair.
    std::vector<DensityAtomDomain> atom_domains;
    atom_domains.reserve(n_atoms * n_domains);
    for (py::ssize_t i = 0; i < n_atoms; ++i) {
        const double delta_z2 = get2(basis_data, basis_info, i, 5)
            / (8.0 * pi * pi);
        const double delta_para2 = get2(basis_data, basis_info, i, 4)
            / (8.0 * pi * pi);
        const double x = get2(basis_data, basis_info, i, 1);
        const double y = get2(basis_data, basis_info, i, 2);
        const double z_frac = get2(basis_data, basis_info, i, 3);
        for (py::ssize_t d = 0; d < n_domains; ++d) {
            const Matrix3 &mat = domain_matrices[d];
            const double x_rel = mat.value[0][0] * x + mat.value[0][1] * y
                + mat.value[0][2] * z_frac
                + get3(domain_data, domain_info, d, 0, 3);
            const double y_rel = mat.value[1][0] * x + mat.value[1][1] * y
                + mat.value[1][2] * z_frac
                + get3(domain_data, domain_info, d, 1, 3);
            const double z_rel = mat.value[2][0] * x + mat.value[2][1] * y
                + mat.value[2][2] * z_frac
                + get3(domain_data, domain_info, d, 2, 3);
            DensityAtomDomain entry{
                delta_z2,
                delta_para2 * q_para2,
                z_rel * c,
                (get2(ff_data, ff_info, i, 10) + get2(ff_data, ff_info, i, 11))
                    / std::sqrt(two_pi * delta_z2),
                get2(ff_data, ff_info, i, 12) / std::sqrt(two_pi * delta_z2),
                {}, {}, {},
                0.0, 0.0,
            };
            for (int term = 0; term < 5; ++term) {
                const double exponent = get2(ff_data, ff_info, i, term + 5);
                const double exp_dpara = exponent + 0.5 * delta_para2;
                entry.exp_dz[term] = exponent + 0.5 * delta_z2;
                entry.term_coefficients[term] = get2(ff_data, ff_info, i, term)
                    / std::sqrt(4.0 * pi * entry.exp_dz[term]);
                entry.exp_dpara_q[term] = exp_dpara * q_para2;
            }
            const double phase = -two_pi * (h * x_rel + k * y_rel);
            const double weight = get2(basis_data, basis_info, i, 6)
                * occupancy_data[d];
            entry.phase_real = weight * std::cos(phase);
            entry.phase_imag = weight * std::sin(phase);
            atom_domains.push_back(entry);
        }
    }

    py::gil_scoped_release release;
    std::fill(out, out + n_z, complex128{0.0, 0.0});
    for (const auto &entry : atom_domains) {
        for (py::ssize_t p = 0; p < n_z; ++p) {
            const double dz = z_data[p] - entry.center_z;
            const double dz_squared = dz * dz;
            const double core = std::exp(-0.5 * (
                entry.delta_para_q + dz_squared / entry.delta_z2
            ));
            double atom_real = entry.base_real * core;
            const double atom_imag = entry.base_imag * core;
            for (int term = 0; term < 5; ++term) {
                atom_real += entry.term_coefficients[term] * std::exp(
                    -entry.exp_dpara_q[term]
                    - dz_squared / (4.0 * entry.exp_dz[term])
                );
            }
            out[p] += complex128{
                entry.phase_real * atom_real - entry.phase_imag * atom_imag,
                entry.phase_real * atom_imag + entry.phase_imag * atom_real,
            };
        }
    }
    for (py::ssize_t p = 0; p < n_z; ++p) {
        out[p] /= uc_area;
    }
    return result;
}

void set_form_factor_cache_budget(const std::size_t bytes) {
    form_factor_cache().set_budget(bytes);
}

std::size_t form_factor_cache_expected_bytes(
    const std::size_t points,
    const std::size_t species
) {
    return points * (species + 1) * sizeof(double);
}

complex128 complex_sinc(const complex128 value) {
    if (value == complex128(0.0, 0.0)) {
        return complex128(1.0, 0.0);
    }
    return std::sin(value) / value;
}

complex128 complex_expm1(const complex128 value) {
    if (std::abs(value) >= 1.0e-5) {
        return std::exp(value) - complex128(1.0, 0.0);
    }
    complex128 sum = value;
    complex128 term = value;
    for (int order = 2; order <= 12; ++order) {
        term *= value / static_cast<double>(order);
        sum += term;
    }
    return sum;
}

complex128 optical_admittance(
    const complex128 refractive_index,
    const complex128 kz,
    const bool p_polarized
) {
    if (p_polarized) {
        return -refractive_index * refractive_index / kz;
    }
    return -kz;
}

void propagate_optical_state(
    complex128 &field,
    complex128 &derivative,
    const complex128 refractive_index,
    const complex128 kz,
    const double thickness,
    const bool p_polarized,
    const bool upward
) {
    const complex128 phase = kz * thickness;
    const complex128 cosine = std::cos(phase);
    const complex128 sinc = complex_sinc(phase);
    complex128 sin_over_admittance;
    complex128 admittance_sin;
    if (p_polarized) {
        sin_over_admittance = (
            kz * kz * thickness * sinc
            / (refractive_index * refractive_index)
        );
        admittance_sin = (
            refractive_index * refractive_index * thickness * sinc
        );
    } else {
        sin_over_admittance = thickness * sinc;
        admittance_sin = kz * kz * thickness * sinc;
    }
    const complex128 direction = upward
        ? complex128(0.0, 1.0)
        : complex128(0.0, -1.0);
    const complex128 next_field = (
        cosine * field + direction * sin_over_admittance * derivative
    );
    const complex128 next_derivative = (
        cosine * derivative + direction * admittance_sin * field
    );
    field = next_field;
    derivative = next_derivative;
}

py::tuple solve_layered_electric_field(
    const ComplexArray1D &refractive_indices,
    const Array1D &z_interfaces,
    const double k0,
    const Array1D &alpha,
    const std::string &polarization
) {
    const auto n_info = refractive_indices.request();
    const auto interface_info = z_interfaces.request();
    const auto alpha_info = alpha.request();
    require_shape(n_info, 1, "refractive_indices");
    require_shape(interface_info, 1, "z_interfaces");
    require_shape(alpha_info, 1, "alpha");
    if (n_info.shape[0] < 2) {
        throw py::value_error("at least two optical media are required");
    }
    if (interface_info.shape[0] != n_info.shape[0] - 1) {
        throw py::value_error("z_interfaces must have one fewer value than media");
    }
    if (!(k0 > 0.0) || !std::isfinite(k0)) {
        throw py::value_error("k0 must be finite and positive");
    }
    if (polarization != "s" && polarization != "p") {
        throw py::value_error("polarization must be 's' or 'p'");
    }

    const auto media = n_info.shape[0];
    const auto points = alpha_info.shape[0];
    auto kz = py::array_t<complex128>({media, points});
    auto A_plus = py::array_t<complex128>({media, points});
    auto A_minus = py::array_t<complex128>({media, points});
    auto A_plus_over_kz = py::array_t<complex128>({media, points});
    auto A_minus_over_kz = py::array_t<complex128>({media, points});
    auto reflection = py::array_t<complex128>(points);
    auto transmission = py::array_t<complex128>(points);
    auto z_reference = py::array_t<double>(media);

    const auto *n_data = static_cast<const complex128 *>(n_info.ptr);
    const auto *interface_data = static_cast<const double *>(interface_info.ptr);
    const auto *alpha_data = static_cast<const double *>(alpha_info.ptr);
    constexpr double half_pi = 1.570796326794896619231321691639751442;
    for (py::ssize_t point = 0; point < points; ++point) {
        if (!std::isfinite(alpha_data[point])
            || alpha_data[point] <= 0.0
            || alpha_data[point] > half_pi) {
            throw py::value_error(
                "alpha must contain values in (0, pi/2] radians"
            );
        }
    }
    auto *kz_data = static_cast<complex128 *>(kz.request().ptr);
    auto *plus_data = static_cast<complex128 *>(A_plus.request().ptr);
    auto *minus_data = static_cast<complex128 *>(A_minus.request().ptr);
    auto *plus_over_data = static_cast<complex128 *>(A_plus_over_kz.request().ptr);
    auto *minus_over_data = static_cast<complex128 *>(A_minus_over_kz.request().ptr);
    auto *reflection_data = static_cast<complex128 *>(reflection.request().ptr);
    auto *transmission_data = static_cast<complex128 *>(transmission.request().ptr);
    auto *reference_data = static_cast<double *>(z_reference.request().ptr);
    reference_data[0] = interface_data[0];
    for (py::ssize_t medium = 1; medium < media; ++medium) {
        reference_data[medium] = interface_data[medium - 1];
    }

    const bool p_polarized = polarization == "p";
    {
    py::gil_scoped_release release;
    std::vector<complex128> field_reference(media);
    std::vector<complex128> derivative_reference(media);
    for (py::ssize_t point = 0; point < points; ++point) {
        const double angle = alpha_data[point];
        const complex128 parallel_over_k0 = n_data[0] * std::cos(angle);
        for (py::ssize_t medium = 0; medium < media; ++medium) {
            complex128 value = -k0 * std::sqrt(
                n_data[medium] * n_data[medium]
                - parallel_over_k0 * parallel_over_k0
            );
            if (value.imag() < 0.0
                || (value.imag() == 0.0 && value.real() > 0.0)) {
                value = -value;
            }
            kz_data[medium * points + point] = value;
        }
        kz_data[point] = -k0 * n_data[0] * std::sin(angle);

        const complex128 terminal_kz = kz_data[(media - 1) * points + point];
        complex128 field_top;
        complex128 derivative_top;
        if (p_polarized) {
            field_top = -terminal_kz / (n_data[media - 1] * n_data[media - 1]);
            derivative_top = complex128(1.0, 0.0);
        } else {
            field_top = complex128(1.0, 0.0);
            derivative_top = -terminal_kz;
        }
        for (py::ssize_t medium = media - 2; medium > 0; --medium) {
            const double thickness = (
                interface_data[medium - 1] - interface_data[medium]
            );
            propagate_optical_state(
                field_top,
                derivative_top,
                n_data[medium],
                kz_data[medium * points + point],
                thickness,
                p_polarized,
                true
            );
        }

        const complex128 incident_admittance = optical_admittance(
            n_data[0], kz_data[point], p_polarized
        );
        const complex128 normalization = 2.0 / (
            field_top + derivative_top / incident_admittance
        );
        field_top *= normalization;
        derivative_top *= normalization;
        const complex128 tangential_reflection = 0.5 * (
            field_top - derivative_top / incident_admittance
        );
        // The transfer state stores tangential fields.  In the right-handed
        // Renaud p basis, an upward branch has the opposite tangential sign,
        // so convert at the public amplitude boundary and keep the stable
        // transfer propagation itself unchanged.
        reflection_data[point] = p_polarized
            ? -tangential_reflection : tangential_reflection;

        field_reference[0] = field_top;
        derivative_reference[0] = derivative_top;
        complex128 field_at_interface = field_top;
        complex128 derivative_at_interface = derivative_top;
        for (py::ssize_t interface = 0; interface < media - 1; ++interface) {
            field_reference[interface + 1] = field_at_interface;
            derivative_reference[interface + 1] = derivative_at_interface;
            if (interface + 1 < media - 1) {
                const double thickness = (
                    interface_data[interface] - interface_data[interface + 1]
                );
                propagate_optical_state(
                    field_at_interface,
                    derivative_at_interface,
                    n_data[interface + 1],
                    kz_data[(interface + 1) * points + point],
                    thickness,
                    p_polarized,
                    false
                );
            }
        }

        for (py::ssize_t medium = 0; medium < media; ++medium) {
            const complex128 medium_kz = kz_data[medium * points + point];
            complex128 difference(0.0, 0.0);
            if (medium_kz != complex128(0.0, 0.0)) {
                if (p_polarized) {
                    difference = (
                        derivative_reference[medium] * (-medium_kz)
                        / (n_data[medium] * n_data[medium])
                    );
                } else {
                    difference = derivative_reference[medium] / (-medium_kz);
                }
            }
            const complex128 plus = 0.5 * (
                field_reference[medium] + difference
            );
            const complex128 tangential_minus = 0.5 * (
                field_reference[medium] - difference
            );
            const complex128 minus = p_polarized
                ? -tangential_minus : tangential_minus;
            plus_data[medium * points + point] = plus;
            minus_data[medium * points + point] = minus;
            if (medium_kz != complex128(0.0, 0.0)) {
                plus_over_data[medium * points + point] = plus / medium_kz;
                minus_over_data[medium * points + point] = minus / medium_kz;
            } else {
                plus_over_data[medium * points + point] = complex128(0.0, 0.0);
                minus_over_data[medium * points + point] = complex128(0.0, 0.0);
            }
        }
        plus_data[point] = complex128(1.0, 0.0);
        minus_data[point] = reflection_data[point];
        plus_data[(media - 1) * points + point] = field_reference[media - 1];
        minus_data[(media - 1) * points + point] = complex128(0.0, 0.0);
        if (terminal_kz != complex128(0.0, 0.0)) {
            plus_over_data[(media - 1) * points + point] = (
                field_reference[media - 1] / terminal_kz
            );
        } else if (p_polarized) {
            plus_over_data[(media - 1) * points + point] = (
                -derivative_reference[media - 1]
                / (n_data[media - 1] * n_data[media - 1])
            );
        }
        minus_over_data[(media - 1) * points + point] = complex128(0.0, 0.0);
        transmission_data[point] = field_reference[media - 1];
    }
    }
    return py::make_tuple(
        kz,
        A_plus,
        A_minus,
        A_plus_over_kz,
        A_minus_over_kz,
        z_reference,
        reflection,
        transmission
    );
}

py::tuple sample_layered_electric_field(
    const ComplexArray2D &kz,
    const ComplexArray2D &A_plus,
    const ComplexArray2D &A_minus,
    const Array1D &z_reference,
    const Array1D &z_interfaces,
    const Array1D &z,
    const std::string &polarization
) {
    const auto kz_info = kz.request();
    const auto plus_info = A_plus.request();
    const auto minus_info = A_minus.request();
    const auto reference_info = z_reference.request();
    const auto interface_info = z_interfaces.request();
    const auto z_info = z.request();
    require_shape(kz_info, 2, "kz");
    require_shape(plus_info, 2, "A_plus");
    require_shape(minus_info, 2, "A_minus");
    require_shape(reference_info, 1, "z_reference");
    require_shape(interface_info, 1, "z_interfaces");
    require_shape(z_info, 1, "z");
    if (polarization != "s" && polarization != "p") {
        throw py::value_error("polarization must be 's' or 'p'");
    }
    if (plus_info.shape != kz_info.shape || minus_info.shape != kz_info.shape) {
        throw py::value_error("kz and amplitude arrays must have matching shapes");
    }
    const auto media = kz_info.shape[0];
    const auto points = kz_info.shape[1];
    if (reference_info.shape[0] != media || interface_info.shape[0] != media - 1) {
        throw py::value_error("field reference and interface shapes do not match media");
    }

    auto E = py::array_t<complex128>({z_info.shape[0], points});
    auto medium_index = py::array_t<std::int64_t>(z_info.shape[0]);
    const auto *kz_data = static_cast<const complex128 *>(kz_info.ptr);
    const auto *plus_data = static_cast<const complex128 *>(plus_info.ptr);
    const auto *minus_data = static_cast<const complex128 *>(minus_info.ptr);
    const auto *reference_data = static_cast<const double *>(reference_info.ptr);
    const auto *interface_data = static_cast<const double *>(interface_info.ptr);
    const auto *z_data = static_cast<const double *>(z_info.ptr);
    auto *E_data = static_cast<complex128 *>(E.request().ptr);
    auto *medium_data = static_cast<std::int64_t *>(medium_index.request().ptr);

    {
    py::gil_scoped_release release;
    for (py::ssize_t sample = 0; sample < z_info.shape[0]; ++sample) {
        py::ssize_t medium = 0;
        while (medium < media - 1 && z_data[sample] < interface_data[medium]) {
            ++medium;
        }
        medium_data[sample] = medium;
        const double depth = reference_data[medium] - z_data[sample];
        for (py::ssize_t point = 0; point < points; ++point) {
            const auto offset = medium * points + point;
            const complex128 phase = std::exp(
                complex128(0.0, 1.0) * kz_data[offset] * depth
            );
            E_data[sample * points + point] = plus_data[offset] * phase;
            if (polarization == "p") {
                E_data[sample * points + point] -= minus_data[offset] / phase;
            } else {
                E_data[sample * points + point] += minus_data[offset] / phase;
            }
        }
    }
    }
    return py::make_tuple(E, medium_index);
}

py::array_t<double> homogeneous_bulk_optical_constants(
    const Array2D &basis,
    const Array2D &f_factors,
    const double scale_over_volume
) {
    const auto basis_info = basis.request();
    const auto factors_info = f_factors.request();
    require_shape(basis_info, 2, "basis");
    require_shape(factors_info, 2, "f_factors");
    if (basis_info.shape[1] < 7
        || factors_info.shape[0] != basis_info.shape[0]
        || factors_info.shape[1] < 13) {
        throw py::value_error("bulk optical atom arrays are inconsistent");
    }
    if (!std::isfinite(scale_over_volume) || scale_over_volume <= 0.0) {
        throw py::value_error("scale_over_volume must be finite and positive");
    }
    auto result = py::array_t<double>(2);
    auto *out = static_cast<double *>(result.request().ptr);
    const auto *basis_data = static_cast<const double *>(basis_info.ptr);
    const auto *factor_data = static_cast<const double *>(factors_info.ptr);
    double real_sum = 0.0;
    double imaginary_sum = 0.0;
    {
    py::gil_scoped_release release;
    for (py::ssize_t atom = 0; atom < basis_info.shape[0]; ++atom) {
        double forward = get2(factor_data, factors_info, atom, 10);
        for (int term = 0; term < 5; ++term) {
            forward += get2(factor_data, factors_info, atom, term);
        }
        const double occupancy = get2(basis_data, basis_info, atom, 6);
        real_sum += occupancy * (
            forward + get2(factor_data, factors_info, atom, 11)
        );
        imaginary_sum += occupancy * get2(
            factor_data, factors_info, atom, 12
        );
    }
    }
    out[0] = scale_over_volume * real_sum;
    out[1] = scale_over_volume * imaginary_sum;
    return result;
}

#include "CTRdwba_cpp.hpp"

PYBIND11_MODULE(_CTRcalc_cpp, module) {
    module.doc() = "CPU C++ acceleration kernels for CTR calculations.";
    module.def("unitcell_F_uc_bulk", &unitcell_F_uc_bulk);
    module.def("unitcell_F_uc_bulk_direct", &unitcell_F_uc_bulk_direct);
    module.def("unitcell_F_bulk", &unitcell_F_bulk);
    module.def("unitcell_F_uc", &unitcell_F_uc);
    module.def("unitcell_zdensity_g", &unitcell_zdensity_g);
    module.def("form_factor_cache_stats", &form_factor_cache_stats);
    module.def("clear_form_factor_cache", []() { form_factor_cache().clear(); });
    module.def("reset_form_factor_cache_stats", &reset_form_factor_cache_stats);
    module.def("set_form_factor_cache_budget", &set_form_factor_cache_budget, py::call_guard<py::gil_scoped_release>());
    module.def("form_factor_cache_expected_bytes", &form_factor_cache_expected_bytes);
    module.def("solve_layered_electric_field", &solve_layered_electric_field);
    module.def("sample_layered_electric_field", &sample_layered_electric_field);
    module.def(
        "homogeneous_bulk_optical_constants",
        &homogeneous_bulk_optical_constants
    );
    module.def("unitcell_F_DWBA_bulk", &unitcell_F_DWBA_bulk);
    module.def("unitcell_F_DWBA_records", &unitcell_F_DWBA_records);
}

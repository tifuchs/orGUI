// Sparse provenance for repaired detector counts. All moments are counts^2.
// Median moments are conditional on one fitted homogeneous Poisson rate per
// connected donor neighborhood; they do not include interpolation/model bias.
// Linear-only neighborhoods retain independent per-donor count estimates.

using RoiMedianKey = std::tuple<int, int, int>;

// Tables depend only on donor counts/overlap, not on detector masks or rates.
// Share immutable tables across jobs, with a bounded process cache. The lock
// is used only during plan setup; frame processing reads plan-owned pointers.
struct RoiMedianTableCache {
    std::mutex mutex;
    std::map<int, std::shared_ptr<const MedianPoissonMomentsTable>> single;
    std::map<RoiMedianKey, std::shared_ptr<const MedianPoissonTable>> pairs;
    std::size_t bytes = 0;
    static constexpr std::size_t limit = 32 * 1024 * 1024;
};

RoiMedianTableCache& repair_median_table_cache() {
    static RoiMedianTableCache cache;
    return cache;
}

std::shared_ptr<const MedianPoissonMomentsTable> repair_median_table(int n) {
    auto& cache = repair_median_table_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    const auto found = cache.single.find(n);
    if (found != cache.single.end()) return found->second;
    auto table = std::make_shared<const MedianPoissonMomentsTable>(n);
    if (cache.bytes + table->storage_bytes() <= RoiMedianTableCache::limit) {
        cache.bytes += table->storage_bytes();
        cache.single[n] = table;
    }
    return table;
}

std::shared_ptr<const MedianPoissonTable> repair_median_pair_table(int n, int m, int k) {
    auto& cache = repair_median_table_cache();
    std::lock_guard<std::mutex> lock(cache.mutex);
    const RoiMedianKey key{n, m, k};
    const auto found = cache.pairs.find(key);
    if (found != cache.pairs.end()) return found->second;
    auto table = std::make_shared<const MedianPoissonTable>(n, m, k);
    if (cache.bytes + table->storage_bytes() <= RoiMedianTableCache::limit) {
        cache.bytes += table->storage_bytes();
        cache.pairs[key] = table;
    }
    return table;
}

struct RoiMedianTableBank {
    std::map<int, std::shared_ptr<const MedianPoissonMomentsTable>> single;
    std::map<RoiMedianKey, std::shared_ptr<const MedianPoissonTable>> pairs;

    MedianPoissonMoments moments(int n, double rate) const {
        if (!median_poisson_supported(n, n, n) || !std::isfinite(rate)
                || rate < 0.0 || rate > median_poisson_max_rate) {
            const double unavailable = std::numeric_limits<double>::quiet_NaN();
            return {unavailable, unavailable, unavailable};
        }
        const auto found = single.find(n);
        return found == single.end() ? median_poisson_interpolated_moments(n, rate)
            : found->second->moments(rate);
    }
    double covariance(int n, int m, int overlap, double rate) const {
        if (!median_poisson_supported(n, m, overlap) || !std::isfinite(rate)
                || rate < 0.0 || rate > median_poisson_max_rate)
            return std::numeric_limits<double>::quiet_NaN();
        if (n > m) std::swap(n, m);
        if (n == m && overlap == n) return moments(n, rate).variance;
        const auto found = pairs.find({n, m, overlap});
        return found == pairs.end() ? median_poisson_interpolated_covariance(n, m, overlap, rate)
            : found->second->covariance(rate);
    }
};

struct RoiRepairRecord {
    std::vector<py::ssize_t> targets;
    std::vector<std::pair<py::ssize_t, double>> donors;
    bool median = false;
    double value = 0.0;
    std::size_t cluster = 0;
    mutable double linear_coefficient = 0.0;
    mutable double residual_variance = 0.0;
    mutable bool moments_ready = false;
};

struct RoiRepairCluster {
    std::vector<py::ssize_t> donors;
    bool median = false;
    double rate = 0.0;
};

struct RoiRepairFrame {
    const RoiMedianTableBank* tables = nullptr;
    std::vector<RoiRepairRecord> records;
    std::vector<RoiRepairCluster> clusters;
    const std::vector<RoiRepairCluster>* cached_clusters = nullptr;
    std::vector<double> rates;

    const RoiRepairCluster& cluster(std::size_t index) const {
        return cached_clusters ? (*cached_clusters)[index] : clusters[index];
    }

    void estimate_rates(const double* image) {
        const auto& groups = cached_clusters ? *cached_clusters : clusters;
        rates.resize(groups.size());
        for (std::size_t i = 0; i < groups.size(); ++i) {
            double sum = 0.0;
            for (const auto donor : groups[i].donors) {
                if (!std::isfinite(image[donor]) || image[donor] < 0.0) {
                    sum = std::numeric_limits<double>::quiet_NaN();
                    break;
                }
                sum += image[donor];
            }
            rates[i] = sum / static_cast<double>(groups[i].donors.size());
        }
    }

    // Rebuild the sparse dependency graph from actual finite donor subsets.
    // This is necessary when NaNs break an otherwise connected neighborhood.
    void connect(const double* image) {
        std::vector<std::size_t> parent(records.size());
        for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = i;
        const auto root = [&](std::size_t i) {
            while (parent[i] != i) {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            return i;
        };
        std::unordered_map<py::ssize_t, std::size_t> owners;
        for (std::size_t i = 0; i < records.size(); ++i) {
            for (const auto& donor : records[i].donors) {
                const auto inserted = owners.emplace(donor.first, i);
                if (!inserted.second) parent[root(i)] = root(inserted.first->second);
            }
        }
        std::unordered_map<std::size_t, std::size_t> groups;
        for (std::size_t i = 0; i < records.size(); ++i) {
            const auto inserted = groups.emplace(root(i), clusters.size());
            if (inserted.second) clusters.emplace_back();
            auto& record = records[i];
            record.cluster = inserted.first->second;
            auto& cluster = clusters[record.cluster];
            cluster.median = cluster.median || record.median;
            for (const auto& donor : record.donors) cluster.donors.push_back(donor.first);
        }
        for (auto& cluster : clusters) {
            auto& donors = cluster.donors;
            std::sort(donors.begin(), donors.end());
            donors.erase(std::unique(donors.begin(), donors.end()), donors.end());
        }
        if (image) estimate_rates(image);
    }
};

struct RoiRepairRectangle {
    std::int64_t x0, x1, y0, y1;
    bool contains(py::ssize_t index, py::ssize_t width) const {
        const auto y = index / width;
        const auto x = index % width;
        return x >= x0 && x < x1 && y >= y0 && y < y1;
    }
};

RoiRepairRectangle repair_rectangle(const py::buffer_info& info, py::ssize_t i) {
    const auto p = roi_ptr(info);
    return {roi_value(p, info, i, 0, 0), roi_value(p, info, i, 0, 1),
            roi_value(p, info, i, 1, 0), roi_value(p, info, i, 1, 1)};
}

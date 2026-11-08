// Model-A median moments enter as an orthogonal linear projection on the
// original donor counts plus a residual. This is an algebraic covariance
// decomposition, not a replacement of the repaired median intensity.

void prepare_repair_moment(const RoiRepairFrame& frame, const RoiRepairRecord& record) {
    if (!record.median || record.moments_ready) return;
    record.moments_ready = true;
    const double rate = frame.rates[record.cluster];
    if (!std::isfinite(rate) || rate > median_poisson_max_rate
            || record.donors.size() > median_poisson_max_donors) {
        record.linear_coefficient = record.residual_variance =
            std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const auto moment = frame.tables
        ? frame.tables->moments(static_cast<int>(record.donors.size()), rate)
        : median_poisson_moments(static_cast<int>(record.donors.size()), rate);
    record.linear_coefficient = rate > 0.0 ? moment.donor_covariance / rate : 0.0;
    record.residual_variance = moment.variance
        - record.donors.size() * record.linear_coefficient * moment.donor_covariance;
}

double repair_overlap_sum(const RoiRepairRectangle& a, const RoiRepairRectangle& b,
                          const double* image, const bool* mask, py::ssize_t width) {
    double result = 0.0;
    for (auto y = std::max(a.y0, b.y0); y < std::min(a.y1, b.y1); ++y) {
        for (auto x = std::max(a.x0, b.x0); x < std::min(a.x1, b.x1); ++x) {
            const auto index = y * width + x;
            if (!mask[index] && !std::isnan(image[index])) result += image[index];
        }
    }
    return result;
}

template<class BackgroundWeight>
void add_repair_covariance(const RoiRepairFrame& frame,
                           const RoiRepairRectangle& center,
                           const double* image, py::ssize_t width, int max_span,
                           double center_pixels, double background_pixels,
                           BackgroundWeight background_weight, double* moments,
                           bool center_includes_repair = false) {
    struct Source {
        py::ssize_t index;
        double added;
        double variance;
    };
    std::vector<std::pair<const RoiRepairRecord*, int>> selected;
    std::vector<std::size_t> clusters;
    std::vector<Source> sources;
    std::vector<RoiRepairRecord>::const_iterator first = std::lower_bound(
        frame.records.begin(), frame.records.end(),
        center.y0 - max_span + 1,
        [width](const RoiRepairRecord& record, std::int64_t row) {
            return record.targets.front() / width < row;
        });
    for (; first != frame.records.end()
            && first->targets.front() / width < center.y1; ++first) {
        const RoiRepairRecord& record = *first;
        int count = 0;
        for (const auto target : record.targets) count += center.contains(target, width);
        if (!count) continue;
        prepare_repair_moment(frame, record);
        if (center_includes_repair) moments[0] -= count * record.value;
        selected.emplace_back(&record, count);
        clusters.push_back(record.cluster);
        const RoiRepairCluster& cluster = frame.cluster(record.cluster);
        for (const auto& donor : record.donors) {
            const double coefficient = record.median
                ? record.linear_coefficient : donor.second;
            sources.push_back({donor.first, count * coefficient,
                cluster.median || !std::isfinite(frame.rates[record.cluster])
                    ? frame.rates[record.cluster] : image[donor.first]});
        }
        if (record.median) moments[0] += count * count * record.residual_variance;
    }
    if (selected.empty()) return;
    std::sort(clusters.begin(), clusters.end());
    clusters.erase(std::unique(clusters.begin(), clusters.end()), clusters.end());
    // A single common rate is used consistently for every directly counted
    // donor in an involved median neighborhood, including fitted backgrounds.
    for (const auto group : clusters) {
        const auto& cluster = frame.cluster(group);
        if (cluster.median) for (const auto donor : cluster.donors) {
            sources.push_back({donor, 0.0, frame.rates[group]});
        }
    }
    std::sort(sources.begin(), sources.end(),
        [](const Source& a, const Source& b) { return a.index < b.index; });
    for (std::size_t i = 0; i < sources.size();) {
        const auto index = sources[i].index;
        const double variance = sources[i].variance;
        double added = 0.0;
        do { added += sources[i++].added; }
        while (i < sources.size() && sources[i].index == index);
        const double signal = center.contains(index, width) ? 1.0 : 0.0;
        const double background = background_weight(index);
        const double combined = signal + added;
        moments[0] += combined * combined * variance - signal * image[index];
        moments[1] += background * background * (variance - image[index]);
        moments[2] += background * (combined * variance - signal * image[index]);
    }
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const auto& a = *selected[i].first;
        if (!a.median) continue;
        for (std::size_t j = 0; j < i; ++j) {
            const auto& b = *selected[j].first;
            if (!b.median || a.cluster != b.cluster) continue;
            std::size_t x = 0, y = 0, overlap = 0;
            while (x < a.donors.size() && y < b.donors.size()) {
                if (a.donors[x].first < b.donors[y].first) ++x;
                else if (a.donors[x].first > b.donors[y].first) ++y;
                else { ++overlap; ++x; ++y; }
            }
            if (!overlap) continue;
            const double rate = frame.rates[a.cluster];
            if (!median_poisson_supported(static_cast<int>(a.donors.size()),
                    static_cast<int>(b.donors.size()), static_cast<int>(overlap))
                    || !std::isfinite(rate) || rate > median_poisson_max_rate) {
                moments[0] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            const double covariance = frame.tables
                ? frame.tables->covariance(static_cast<int>(a.donors.size()),
                    static_cast<int>(b.donors.size()), static_cast<int>(overlap), rate)
                : median_poisson_covariance(static_cast<int>(a.donors.size()),
                    static_cast<int>(b.donors.size()), static_cast<int>(overlap), rate);
            const double residual = covariance
                - overlap * a.linear_coefficient * b.linear_coefficient * rate;
            moments[0] += 2.0 * selected[i].second * selected[j].second * residual;
        }
    }
    const double ratio = background_pixels > 0.0 ? center_pixels / background_pixels : 0.0;
    moments[3] = moments[0] + ratio * ratio * moments[1] - 2.0 * ratio * moments[2];
}

    int maximum_span() const { return max_span_; }

    void require_shape(py::ssize_t height, py::ssize_t width) const {
        if (height != height_ || width != width_)
            throw py::value_error("image must match repair plan shape");
    }

    void repair_footprint(bool* output) const {
        std::fill(output, output + size_, false);
        for (const auto& component : components_)
            for (const auto target : component.targets) output[target] = true;
    }

    void prepare_moment_tables(const bool* region = nullptr) {
        moment_tables_ = std::make_unique<RoiMedianTableBank>();
        RoiRepairFrame geometry;
        for (const auto& component : components_) {
            RoiRepairRecord record;
            record.median = component.median;
            for (const auto& candidate : component.candidates)
                record.donors.emplace_back(candidate.index, 0.0);
            geometry.records.push_back(std::move(record));
        }
        geometry.connect(nullptr);
        moment_clusters_ = std::move(geometry.clusters);
        for (const auto& record : geometry.records) component_clusters_.push_back(record.cluster);
        std::unordered_map<py::ssize_t, std::vector<std::size_t>> owners;
        std::map<std::pair<std::size_t, std::size_t>, int> overlaps;
        for (std::size_t i = 0; i < components_.size(); ++i) {
            const auto& component = components_[i];
            if (!component.median) continue;
            if (region && std::none_of(component.targets.begin(), component.targets.end(),
                    [region](py::ssize_t target) { return region[target]; })) continue;
            const int n = static_cast<int>(component.candidates.size());
            if (!median_poisson_supported(n, n, n)) {
                ++unsupported_median_components_;
                continue;
            }
            if (!moment_tables_->single.count(n)) {
                moment_tables_->single[n] = repair_median_table(n);
            }
            for (const auto& candidate : component.candidates) {
                auto& previous = owners[candidate.index];
                for (const auto j : previous) ++overlaps[{j, i}];
                previous.push_back(i);
            }
        }
        for (const auto& overlap : overlaps) {
            int n = static_cast<int>(components_[overlap.first.first].candidates.size());
            int m = static_cast<int>(components_[overlap.first.second].candidates.size());
            if (n > m) std::swap(n, m);
            const RoiMedianKey key{n, m, overlap.second};
            if (n == m && overlap.second == n) continue;
            if (!moment_tables_->pairs.count(key)) {
                moment_tables_->pairs[key] = repair_median_pair_table(n, m, overlap.second);
            }
        }
    }

    // ROI-only provenance path. Reconstruction's released conditional-donor
    // variance API remains separate; the ROI consumer evaluates model moments.
    RoiRepairFrame apply_roi_records(Array2D image, BoolArray2D repaired,
                                     bool provenance = true) const {
        const auto image_info = image.request();
        const auto repaired_info = repaired.request();
        require_same_image_shape(image_info, repaired_info, "repaired");
        if (image_info.shape[0] != height_ || image_info.shape[1] != width_) {
            throw py::value_error("image must match repair plan shape");
        }
        double* data = mutable_double_ptr(image_info);
        bool* repaired_data = static_cast<bool*>(repaired_info.ptr);
        RoiRepairFrame frame;
        frame.tables = moment_tables_.get();
        bool complete_geometry = moment_tables_ != nullptr;
        py::gil_scoped_release release;
        std::vector<const PlannedRepairCandidate*> valid;
        std::vector<double> values;
        if (provenance) frame.records.reserve(components_.size());
        for (std::size_t component_index = 0; component_index < components_.size(); ++component_index) {
            const auto& component = components_[component_index];
            valid.clear();
            unsigned char sides = 0;
            for (const auto& candidate : component.candidates) {
                if (std::isnan(data[candidate.index])) {
                    complete_geometry = false;
                    continue;
                }
                valid.push_back(&candidate);
                sides |= candidate.sides;
            }
            if (static_cast<int>(valid.size()) < min_valid_neighbors_
                    || side_count(sides) < 2) {
                complete_geometry = false;
                continue;
            }
            RoiRepairRecord record;
            if (provenance) record.targets = component.targets;
            record.median = component.median;
            if (moment_tables_) record.cluster = component_clusters_[component_index];
            if (provenance) record.donors.reserve(valid.size());
            if (component.median) {
                values.clear();
                for (const auto candidate : valid) {
                    values.push_back(data[candidate->index]);
                    if (provenance) record.donors.emplace_back(candidate->index, 0.0);
                }
                std::sort(values.begin(), values.end());
                const auto mid = values.size() / 2;
                record.value = values.size() % 2 ? values[mid]
                    : 0.5 * (values[mid - 1] + values[mid]);
            } else {
                double total = 0.0;
                double numerator = 0.0;
                for (const auto candidate : valid) {
                    total += candidate->inverse_distance2;
                    numerator += data[candidate->index] * candidate->inverse_distance2;
                }
                record.value = numerator / total;
                for (const auto candidate : valid) {
                    const double weight = candidate->inverse_distance2 / total;
                    if (provenance) record.donors.emplace_back(candidate->index, weight);
                }
            }
            for (const auto target : component.targets) {
                data[target] = record.value;
                repaired_data[target] = true;
            }
            if (provenance) frame.records.push_back(std::move(record));
        }
        if (provenance) {
            if (complete_geometry) {
                frame.cached_clusters = &moment_clusters_;
                frame.estimate_rates(data);
            } else {
                frame.connect(data);
            }
        }
        return frame;
    }

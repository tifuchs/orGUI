// Native Renaud DWBA kernel, included after CTRcalc_cpp.cpp shared helpers.
py::array_t<complex128> unitcell_F_DWBA_bulk(
    const Array2D &Q_parallel,
    const IntArray1D &is_specular,
    const ComplexArray1D &reference_refractive_indices,
    const IntArray1D &field_i_index,
    const IntArray1D &field_f_index,
    const ComplexArray2D &kz_i,
    const ComplexArray2D &A_i_plus,
    const ComplexArray2D &A_i_minus,
    const ComplexArray2D &A_i_plus_over_kz,
    const ComplexArray2D &A_i_minus_over_kz,
    const ComplexArray2D &kz_f,
    const ComplexArray2D &A_f_plus,
    const ComplexArray2D &A_f_minus,
    const ComplexArray2D &A_f_plus_over_kz,
    const ComplexArray2D &A_f_minus_over_kz,
    const Array1D &z_reference_i,
    const Array1D &z_reference_f,
    const Array1D &alpha_i,
    const Array1D &alpha_f,
    const Array1D &cos_azimuth,
    const Array1D &sin_azimuth,
    const double k0,
    const std::string &polarization_i,
    const std::string &polarization_f,
    const IntArray1D &medium_index,
    const bool semi_infinite,
    const double attenuation,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &r_mat,
    const Array2D &r_mat_inv,
    const Array3D &coherent_domain_matrix,
    const Array1D &coherent_domain_occupancy
) {
    const auto Q_info = Q_parallel.request();
    const auto specular_info = is_specular.request();
    const auto reference_n_info = reference_refractive_indices.request();
    const auto field_i_index_info = field_i_index.request();
    const auto field_f_index_info = field_f_index.request();
    const auto kzi_info = kz_i.request();
    const auto Aip_info = A_i_plus.request();
    const auto Aim_info = A_i_minus.request();
    const auto Aipk_info = A_i_plus_over_kz.request();
    const auto Aimk_info = A_i_minus_over_kz.request();
    const auto kzf_info = kz_f.request();
    const auto Afp_info = A_f_plus.request();
    const auto Afm_info = A_f_minus.request();
    const auto Afpk_info = A_f_plus_over_kz.request();
    const auto Afmk_info = A_f_minus_over_kz.request();
    const auto zrefi_info = z_reference_i.request();
    const auto zreff_info = z_reference_f.request();
    const auto alphai_info = alpha_i.request();
    const auto alphaf_info = alpha_f.request();
    const auto cos_info = cos_azimuth.request();
    const auto sin_info = sin_azimuth.request();
    const auto medium_info = medium_index.request();
    const auto basis_info = basis.request();
    const auto factors_info = f_factors.request();
    const auto r_info = r_mat.request();
    const auto r_inv_info = r_mat_inv.request();
    const auto domain_info = coherent_domain_matrix.request();
    const auto occupancy_info = coherent_domain_occupancy.request();

    require_shape(Q_info, 2, "Q_parallel");
    require_shape(specular_info, 1, "is_specular");
    require_shape(reference_n_info, 1, "reference_refractive_indices");
    require_shape(field_i_index_info, 1, "field_i_index");
    require_shape(field_f_index_info, 1, "field_f_index");
    require_shape(kzi_info, 2, "kz_i");
    require_shape(Aip_info, 2, "A_i_plus");
    require_shape(Aim_info, 2, "A_i_minus");
    require_shape(Aipk_info, 2, "A_i_plus_over_kz");
    require_shape(Aimk_info, 2, "A_i_minus_over_kz");
    require_shape(kzf_info, 2, "kz_f");
    require_shape(Afp_info, 2, "A_f_plus");
    require_shape(Afm_info, 2, "A_f_minus");
    require_shape(Afpk_info, 2, "A_f_plus_over_kz");
    require_shape(Afmk_info, 2, "A_f_minus_over_kz");
    require_shape(zrefi_info, 1, "z_reference_i");
    require_shape(zreff_info, 1, "z_reference_f");
    require_shape(alphai_info, 1, "alpha_i");
    require_shape(alphaf_info, 1, "alpha_f");
    require_shape(cos_info, 1, "cos_azimuth");
    require_shape(sin_info, 1, "sin_azimuth");
    require_shape(medium_info, 1, "medium_index");
    require_shape(basis_info, 2, "basis");
    require_shape(factors_info, 2, "f_factors");
    require_shape(r_info, 2, "R_mat");
    require_shape(r_inv_info, 2, "R_mat_inv");
    require_shape(domain_info, 3, "coherentDomainMatrix");
    require_shape(occupancy_info, 1, "coherentDomainOccupancy");

    const auto points = Q_info.shape[0];
    const auto media = kzi_info.shape[0];
    const auto incident_fields = kzi_info.shape[1];
    const auto exit_fields = kzf_info.shape[1];
    const auto atoms = basis_info.shape[0];
    const auto domains = domain_info.shape[0];
    if (Q_info.shape[1] != 2) {
        throw py::value_error("Q_parallel must have shape (N, 2)");
    }
    if (specular_info.shape[0] != points) {
        throw py::value_error("is_specular must match the point count");
    }
    if (field_i_index_info.shape[0] != points
        || field_f_index_info.shape[0] != points) {
        throw py::value_error("field index arrays must match the point count");
    }
    if (reference_n_info.shape[0] != media) {
        throw py::value_error(
            "reference_refractive_indices must match the media count"
        );
    }
    if (kzf_info.shape[0] != media || incident_fields == 0
        || exit_fields == 0) {
        throw py::value_error(
            "incident and exit field tables must contain the same media"
        );
    }
    const auto same_incident_shape = [media, incident_fields](
        const py::buffer_info &info
    ) {
        return info.shape[0] == media && info.shape[1] == incident_fields;
    };
    const auto same_exit_shape = [media, exit_fields](
        const py::buffer_info &info
    ) {
        return info.shape[0] == media && info.shape[1] == exit_fields;
    };
    if (!same_incident_shape(Aip_info) || !same_incident_shape(Aim_info)
        || !same_incident_shape(Aipk_info) || !same_incident_shape(Aimk_info)
        || !same_exit_shape(Afp_info) || !same_exit_shape(Afm_info)
        || !same_exit_shape(Afpk_info) || !same_exit_shape(Afmk_info)) {
        throw py::value_error(
            "incident and exit field-table arrays have inconsistent shapes"
        );
    }
    if (zrefi_info.shape[0] != media || zreff_info.shape[0] != media) {
        throw py::value_error("z-reference arrays must match the media count");
    }
    if (alphai_info.shape[0] != points || alphaf_info.shape[0] != points
        || cos_info.shape[0] != points || sin_info.shape[0] != points) {
        throw py::value_error("geometry arrays must match the point count");
    }
    if (basis_info.shape[1] < 7 || factors_info.shape[0] != atoms
        || factors_info.shape[1] < 13 || medium_info.shape[0] != atoms) {
        throw py::value_error("atom, factor, and medium arrays are inconsistent");
    }
    if (r_info.shape[0] != 3 || r_info.shape[1] != 3
        || r_inv_info.shape[0] != 3 || r_inv_info.shape[1] != 3
        || domain_info.shape[1] != 3 || domain_info.shape[2] != 4
        || occupancy_info.shape[0] != domains) {
        throw py::value_error("domain and real-space transforms have invalid shapes");
    }
    if (polarization_i != "s" && polarization_i != "p") {
        throw py::value_error("polarization_i must be 's' or 'p'");
    }
    if (polarization_f != "s" && polarization_f != "p") {
        throw py::value_error("polarization_f must be 's' or 'p'");
    }
    if (!std::isfinite(attenuation) || attenuation < 0.0) {
        throw py::value_error("attenuation must be finite and non-negative");
    }
    const auto *reference_n_data = static_cast<const complex128 *>(
        reference_n_info.ptr
    );
    const auto *specular_data = static_cast<const std::int64_t *>(
        specular_info.ptr
    );
    const auto *field_i_index_data = static_cast<const std::int64_t *>(
        field_i_index_info.ptr
    );
    const auto *field_f_index_data = static_cast<const std::int64_t *>(
        field_f_index_info.ptr
    );
    for (py::ssize_t medium = 0; medium < media; ++medium) {
        if (!std::isfinite(reference_n_data[medium].real())
            || !std::isfinite(reference_n_data[medium].imag())) {
            throw py::value_error(
                "reference_refractive_indices must be finite"
            );
        }
    }
    for (py::ssize_t point = 0; point < points; ++point) {
        if (specular_data[point] != 0 && specular_data[point] != 1) {
            throw py::value_error("is_specular must contain only zero or one");
        }
        if (field_i_index_data[point] < 0
            || field_i_index_data[point] >= incident_fields
            || field_f_index_data[point] < 0
            || field_f_index_data[point] >= exit_fields) {
            throw py::value_error("field index contains an invalid column");
        }
    }
    const auto *medium_data = static_cast<const std::int64_t *>(medium_info.ptr);
    for (py::ssize_t atom = 0; atom < atoms; ++atom) {
        if (medium_data[atom] < 0 || medium_data[atom] >= media) {
            throw py::value_error("medium_index contains an invalid medium");
        }
    }

    auto result = py::array_t<complex128>(points);
    auto *out = static_cast<complex128 *>(result.request().ptr);
    const auto *Q_data = static_cast<const double *>(Q_info.ptr);
    const auto *kzi_data = static_cast<const complex128 *>(kzi_info.ptr);
    const auto *Aip_data = static_cast<const complex128 *>(Aip_info.ptr);
    const auto *Aim_data = static_cast<const complex128 *>(Aim_info.ptr);
    const auto *Aipk_data = static_cast<const complex128 *>(Aipk_info.ptr);
    const auto *Aimk_data = static_cast<const complex128 *>(Aimk_info.ptr);
    const auto *kzf_data = static_cast<const complex128 *>(kzf_info.ptr);
    const auto *Afp_data = static_cast<const complex128 *>(Afp_info.ptr);
    const auto *Afm_data = static_cast<const complex128 *>(Afm_info.ptr);
    const auto *Afpk_data = static_cast<const complex128 *>(Afpk_info.ptr);
    const auto *Afmk_data = static_cast<const complex128 *>(Afmk_info.ptr);
    const auto *zrefi_data = static_cast<const double *>(zrefi_info.ptr);
    const auto *zreff_data = static_cast<const double *>(zreff_info.ptr);
    const auto *alphai_data = static_cast<const double *>(alphai_info.ptr);
    const auto *alphaf_data = static_cast<const double *>(alphaf_info.ptr);
    const auto *cos_data = static_cast<const double *>(cos_info.ptr);
    const auto *sin_data = static_cast<const double *>(sin_info.ptr);
    const auto *basis_data = static_cast<const double *>(basis_info.ptr);
    const auto *factor_data = static_cast<const double *>(factors_info.ptr);
    const auto *r_data = static_cast<const double *>(r_info.ptr);
    const auto *occupancy_data = static_cast<const double *>(occupancy_info.ptr);
    const auto domain_matrices = effective_domain_matrices(
        domain_info, r_info, r_inv_info
    );

    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double classical_electron_radius = 2.8179403262e-5;
    constexpr double dw_denominator = 16.0 * pi * pi;
    constexpr int sigma_i[4] = {1, 1, -1, -1};
    constexpr int sigma_f[4] = {1, -1, 1, -1};
    const bool p_i = polarization_i == "p";
    const bool p_f = polarization_f == "p";
    const double c_vector[3] = {r_data[2], r_data[5], r_data[8]};
    const double a_vector[3] = {r_data[0], r_data[3], r_data[6]};
    const double b_vector[3] = {r_data[1], r_data[4], r_data[7]};
    const double area_cross[3] = {
        a_vector[1] * b_vector[2] - a_vector[2] * b_vector[1],
        a_vector[2] * b_vector[0] - a_vector[0] * b_vector[2],
        a_vector[0] * b_vector[1] - a_vector[1] * b_vector[0],
    };
    const double bulk_uc_area = std::sqrt(
        area_cross[0] * area_cross[0]
        + area_cross[1] * area_cross[1]
        + area_cross[2] * area_cross[2]
    );
    if (!(bulk_uc_area > 0.0) || !std::isfinite(bulk_uc_area)) {
        throw py::value_error("bulk lateral unit-cell area must be positive");
    }
    if (!(c_vector[2] > 0.0) || !std::isfinite(c_vector[2])) {
        throw py::value_error(
            "bulk surface-normal repeat must be finite and positive"
        );
    }
    const complex128 bulk_scattering_density = (
        k0 * k0 / (2.0 * pi * classical_electron_radius)
    ) * (complex128(1.0, 0.0) - reference_n_data[media - 1]);

    {
    py::gil_scoped_release release;
    const auto geometry_index = [points](
        const py::ssize_t medium,
        const int channel,
        const py::ssize_t point
    ) {
        return (medium * 4 + channel) * points + point;
    };
    const std::size_t geometry_size = static_cast<std::size_t>(
        media * 4 * points
    );
    std::vector<complex128> channel_Qz(geometry_size);
    std::vector<complex128> channel_Q2(geometry_size);
    std::vector<complex128> channel_coefficient(geometry_size);
    for (py::ssize_t medium = 0; medium < media; ++medium) {
        for (int channel = 0; channel < 4; ++channel) {
            for (py::ssize_t point = 0; point < points; ++point) {
                const auto incident_offset = medium * incident_fields
                    + field_i_index_data[point];
                const auto exit_offset = medium * exit_fields
                    + field_f_index_data[point];
                const auto index = geometry_index(medium, channel, point);
                const double Qx = Q_data[point * 2];
                const double Qy = Q_data[point * 2 + 1];
                const double Q_parallel_squared = Qx * Qx + Qy * Qy;
                const complex128 incident_amplitude = sigma_i[channel] > 0
                    ? Aip_data[incident_offset] : Aim_data[incident_offset];
                const complex128 exit_amplitude = sigma_f[channel] > 0
                    ? Afp_data[exit_offset] : Afm_data[exit_offset];
                const complex128 Qz = -(
                    static_cast<double>(sigma_i[channel])
                    * kzi_data[incident_offset]
                    + static_cast<double>(sigma_f[channel])
                    * kzf_data[exit_offset]
                );
                channel_Qz[index] = Qz;
                channel_Q2[index] = Q_parallel_squared + Qz * Qz;
                if (incident_amplitude == complex128(0.0, 0.0)
                    || exit_amplitude == complex128(0.0, 0.0)) {
                    channel_coefficient[index] = complex128(0.0, 0.0);
                    continue;
                }
                const double incident_sine = std::sin(alphai_data[point]);
                const double exit_sine = std::sin(alphaf_data[point]);
                complex128 vector_weight;
                if (!p_i && !p_f) {
                    vector_weight = (
                        incident_amplitude * exit_amplitude * cos_data[point]
                    );
                } else if (!p_i && p_f) {
                    vector_weight = (
                        -incident_amplitude
                        * static_cast<double>(sigma_f[channel])
                        * exit_amplitude
                        * exit_sine * sin_data[point]
                    );
                } else if (p_i && !p_f) {
                    vector_weight = (
                        -static_cast<double>(sigma_i[channel])
                        * incident_amplitude * incident_sine
                        * exit_amplitude * sin_data[point]
                    );
                } else {
                    const complex128 incident_over_kz = sigma_i[channel] > 0
                        ? Aipk_data[incident_offset] : Aimk_data[incident_offset];
                    const complex128 exit_over_kz = sigma_f[channel] > 0
                        ? Afpk_data[exit_offset] : Afmk_data[exit_offset];
                    const complex128 tangent = (
                        -static_cast<double>(sigma_i[channel])
                        * static_cast<double>(sigma_f[channel])
                        * incident_amplitude * exit_amplitude
                        * incident_sine * exit_sine * cos_data[point]
                    );
                    const complex128 incident_normal = (
                        -k0 * std::cos(alphai_data[point]) * incident_sine
                        * incident_over_kz
                    );
                    const complex128 exit_normal = (
                        -k0 * std::cos(alphaf_data[point]) * exit_sine
                        * exit_over_kz
                    );
                    vector_weight = tangent + incident_normal * exit_normal;
                }
                const complex128 reference_phase = std::exp(
                    complex128(0.0, 1.0) * (
                        static_cast<double>(sigma_i[channel])
                        * kzi_data[incident_offset] * zrefi_data[medium]
                        + static_cast<double>(sigma_f[channel])
                        * kzf_data[exit_offset] * zreff_data[medium]
                    )
                );
                channel_coefficient[index] = vector_weight * reference_phase;
            }
        }
    }

    struct AtomGroup {
        py::ssize_t medium;
        std::array<double, 13> factor_row;
        std::array<std::shared_ptr<const std::vector<double>>, 4> factors;
    };
    std::vector<AtomGroup> groups;
    std::vector<std::size_t> atom_group(static_cast<std::size_t>(atoms));
    for (py::ssize_t atom = 0; atom < atoms; ++atom) {
        std::array<double, 13> row;
        for (int column = 0; column < 13; ++column) {
            row[column] = get2(factor_data, factors_info, atom, column);
        }
        const auto medium = static_cast<py::ssize_t>(medium_data[atom]);
        std::size_t group_index = groups.size();
        for (std::size_t candidate = 0; candidate < groups.size(); ++candidate) {
            if (groups[candidate].medium == medium
                && std::memcmp(
                    groups[candidate].factor_row.data(),
                    row.data(),
                    sizeof(row)
                ) == 0) {
                group_index = candidate;
                break;
            }
        }
        if (group_index == groups.size()) {
            groups.push_back(AtomGroup{medium, row, {}});
        }
        atom_group[static_cast<std::size_t>(atom)] = group_index;
    }
    for (auto &group : groups) {
        for (int channel = 0; channel < 4; ++channel) {
            std::vector<double> complex_grid(
                1 + static_cast<std::size_t>(2 * points)
            );
            complex_grid[0] = std::numeric_limits<double>::quiet_NaN();
            for (py::ssize_t point = 0; point < points; ++point) {
                const complex128 value = channel_Q2[
                    geometry_index(group.medium, channel, point)
                ];
                complex_grid[1 + 2 * point] = value.real();
                complex_grid[2 + 2 * point] = value.imag();
            }
            auto factors = form_factor_cache().find(
                complex_grid, group.factor_row
            );
            if (!factors) {
                std::vector<double> calculated(
                    static_cast<std::size_t>(2 * points)
                );
                for (py::ssize_t point = 0; point < points; ++point) {
                    const complex128 Q_squared = channel_Q2[
                        geometry_index(group.medium, channel, point)
                    ];
                    complex128 value(
                        group.factor_row[10] + group.factor_row[11],
                        group.factor_row[12]
                    );
                    for (int term = 0; term < 5; ++term) {
                        value += group.factor_row[term] * std::exp(
                            -group.factor_row[term + 5] * Q_squared
                        );
                    }
                    calculated[2 * point] = value.real();
                    calculated[2 * point + 1] = value.imag();
                }
                factors = form_factor_cache().insert_or_get(
                    complex_grid,
                    group.factor_row,
                    std::move(calculated)
                );
            }
            group.factors[channel] = factors;
        }
    }

    const bool apply_empirical_attenuation = (
        semi_infinite && attenuation != 0.0
    );
    std::vector<double> atom_positions(
        static_cast<std::size_t>(atoms * domains * 3)
    );
    std::vector<double> atom_repeat_coordinates;
    if (apply_empirical_attenuation) {
        atom_repeat_coordinates.resize(
            static_cast<std::size_t>(atoms * domains)
        );
    }
    const auto *domain_data = static_cast<const double *>(domain_info.ptr);
    for (py::ssize_t atom = 0; atom < atoms; ++atom) {
        const double fractional[3] = {
            get2(basis_data, basis_info, atom, 1),
            get2(basis_data, basis_info, atom, 2),
            get2(basis_data, basis_info, atom, 3),
        };
        for (py::ssize_t domain = 0; domain < domains; ++domain) {
            double transformed[3] = {0.0, 0.0, 0.0};
            for (int row = 0; row < 3; ++row) {
                transformed[row] = get3(
                    domain_data, domain_info, domain, row, 3
                );
                for (int col = 0; col < 3; ++col) {
                    transformed[row] += (
                        domain_matrices[domain].value[row][col]
                        * fractional[col]
                    );
                }
            }
            for (int row = 0; row < 3; ++row) {
                double cartesian = 0.0;
                for (int col = 0; col < 3; ++col) {
                    cartesian += r_data[row * 3 + col] * transformed[col];
                }
                atom_positions[
                    static_cast<std::size_t>((atom * domains + domain) * 3 + row)
                ] = cartesian;
            }
            if (apply_empirical_attenuation) {
                atom_repeat_coordinates[
                    static_cast<std::size_t>(atom * domains + domain)
                ] = transformed[2];
            }
        }
    }

    std::vector<double> attenuated_domain_weights;
    if (apply_empirical_attenuation) {
        attenuated_domain_weights.resize(
            static_cast<std::size_t>(atoms * domains)
        );
        for (py::ssize_t atom = 0; atom < atoms; ++atom) {
            for (py::ssize_t domain = 0; domain < domains; ++domain) {
                const auto offset = static_cast<std::size_t>(
                    atom * domains + domain
                );
                attenuated_domain_weights[offset] = occupancy_data[domain]
                    * std::exp(
                        attenuation * atom_repeat_coordinates[offset]
                    );
            }
        }
    }

    for (py::ssize_t point = 0; point < points; ++point) {
        const double Qx = Q_data[point * 2];
        const double Qy = Q_data[point * 2 + 1];
        const double Q_parallel_squared = Qx * Qx + Qy * Qy;
        complex128 total(0.0, 0.0);

        for (int channel = 0; channel < 4; ++channel) {
            complex128 finite_atoms(0.0, 0.0);
            complex128 bulk_atoms(0.0, 0.0);
            for (py::ssize_t atom = 0; atom < atoms; ++atom) {
                const auto medium = static_cast<py::ssize_t>(medium_data[atom]);
                const auto index = geometry_index(medium, channel, point);
                const complex128 coefficient = channel_coefficient[index];
                if (coefficient == complex128(0.0, 0.0)) {
                    continue;
                }
                const complex128 Qz = channel_Qz[index];
                const auto &cached = *groups[
                    atom_group[static_cast<std::size_t>(atom)]
                ].factors[channel];
                complex128 atomic_factor(
                    cached[2 * point], cached[2 * point + 1]
                );
                atomic_factor *= std::exp(
                    -(
                        get2(basis_data, basis_info, atom, 4)
                        * Q_parallel_squared
                        + get2(basis_data, basis_info, atom, 5) * Qz * Qz
                    ) / dw_denominator
                );
                atomic_factor *= get2(basis_data, basis_info, atom, 6);

                complex128 domain_sum(0.0, 0.0);
                if (apply_empirical_attenuation && medium == media - 1) {
                    for (py::ssize_t domain = 0; domain < domains; ++domain) {
                        const auto domain_offset = static_cast<std::size_t>(
                            atom * domains + domain
                        );
                        const auto position_offset = 3 * domain_offset;
                        domain_sum += attenuated_domain_weights[domain_offset]
                            * std::exp(
                                complex128(0.0, 1.0) * (
                                    Qx * atom_positions[position_offset]
                                    + Qy * atom_positions[position_offset + 1]
                                    + Qz * atom_positions[position_offset + 2]
                                )
                            );
                    }
                } else {
                    for (py::ssize_t domain = 0; domain < domains; ++domain) {
                        const auto position_offset = static_cast<std::size_t>(
                            (atom * domains + domain) * 3
                        );
                        domain_sum += occupancy_data[domain] * std::exp(
                            complex128(0.0, 1.0) * (
                                Qx * atom_positions[position_offset]
                                + Qy * atom_positions[position_offset + 1]
                                + Qz * atom_positions[position_offset + 2]
                            )
                        );
                    }
                }
                const complex128 contribution = (
                    coefficient * atomic_factor * domain_sum
                );
                if (semi_infinite && medium == media - 1) {
                    bulk_atoms += contribution;
                } else {
                    finite_atoms += contribution;
                }
            }
            if (semi_infinite) {
                const py::ssize_t terminal_incident_offset =
                    (media - 1) * incident_fields + field_i_index_data[point];
                const py::ssize_t terminal_exit_offset =
                    (media - 1) * exit_fields + field_f_index_data[point];
                const complex128 Qz = -(
                    static_cast<double>(sigma_i[channel])
                    * kzi_data[terminal_incident_offset]
                    + static_cast<double>(sigma_f[channel])
                    * kzf_data[terminal_exit_offset]
                );
                const complex128 repeat_phase = (
                    Qx * c_vector[0] + Qy * c_vector[1] + Qz * c_vector[2]
                );
                const complex128 denominator = -complex_expm1(
                    -complex128(0.0, 1.0) * repeat_phase - attenuation
                );
                if (specular_data[point] != 0 && channel == 0) {
                    const auto terminal_index = geometry_index(
                        media - 1, channel, point
                    );
                    const complex128 mean_decay = (
                        complex128(0.0, 1.0) * Qz
                        + attenuation / c_vector[2]
                    );
                    complex128 mean_depth;
                    if (mean_decay == complex128(0.0, 0.0)) {
                        mean_depth = c_vector[2];
                    } else {
                        mean_depth = -complex_expm1(
                            -mean_decay * c_vector[2]
                        ) / mean_decay;
                    }
                    const double bulk_top = zrefi_data[media - 1];
                    const complex128 boundary_phase = std::exp(
                        complex128(0.0, 1.0) * Qz * bulk_top
                    );
                    bulk_atoms -= (
                        channel_coefficient[terminal_index]
                        * boundary_phase
                        * bulk_uc_area
                        * bulk_scattering_density
                        * mean_depth
                    );
                }
                if (bulk_atoms != complex128(0.0, 0.0)) {
                    const double pole_tolerance = 16.0
                        * std::numeric_limits<double>::epsilon()
                        * (1.0 + std::abs(repeat_phase));
                    if (attenuation == 0.0
                        && std::abs(denominator) <= pole_tolerance) {
                        throw std::domain_error(
                            "zero-attenuation semi-infinite bulk has an exact Bragg pole"
                        );
                    }
                    bulk_atoms /= denominator;
                }
            }
            total += finite_atoms + bulk_atoms;
        }
        out[point] = total;
    }
    }
    return result;
}

py::tuple unitcell_F_DWBA_records(
    const Array2D &Q_parallel,
    const IntArray1D &is_specular,
    const ComplexArray1D &reference_refractive_indices,
    const IntArray1D &field_i_index,
    const IntArray1D &field_f_index,
    const ComplexArray2D &kz_i,
    const ComplexArray2D &A_i_plus,
    const ComplexArray2D &A_i_minus,
    const ComplexArray2D &A_i_plus_over_kz,
    const ComplexArray2D &A_i_minus_over_kz,
    const ComplexArray2D &kz_f,
    const ComplexArray2D &A_f_plus,
    const ComplexArray2D &A_f_minus,
    const ComplexArray2D &A_f_plus_over_kz,
    const ComplexArray2D &A_f_minus_over_kz,
    const Array1D &z_reference_i,
    const Array1D &z_reference_f,
    const Array1D &z_interfaces,
    const Array1D &alpha_i,
    const Array1D &alpha_f,
    const Array1D &cos_azimuth,
    const Array1D &sin_azimuth,
    const double k0,
    const std::string &polarization_i,
    const std::string &polarization_f,
    const Array3D &reference_delta_beta,
    const double reference_area,
    const bool semi_infinite,
    const double attenuation,
    const Array2D &basis,
    const Array2D &f_factors,
    const Array2D &finite_positions,
    const IntArray1D &finite_atom_index,
    const IntArray1D &finite_record_index,
    const Array1D &finite_weight,
    const Array2D &bulk_positions,
    const IntArray1D &bulk_atom_index,
    const Array1D &bulk_weight,
    const Array1D &bulk_repeat_coordinate,
    const Array1D &bulk_repeat
) {
    const auto Q_info = Q_parallel.request();
    const auto specular_info = is_specular.request();
    const auto reference_n_info = reference_refractive_indices.request();
    const auto field_i_index_info = field_i_index.request();
    const auto field_f_index_info = field_f_index.request();
    const auto kzi_info = kz_i.request();
    const auto Aip_info = A_i_plus.request();
    const auto Aim_info = A_i_minus.request();
    const auto Aipk_info = A_i_plus_over_kz.request();
    const auto Aimk_info = A_i_minus_over_kz.request();
    const auto kzf_info = kz_f.request();
    const auto Afp_info = A_f_plus.request();
    const auto Afm_info = A_f_minus.request();
    const auto Afpk_info = A_f_plus_over_kz.request();
    const auto Afmk_info = A_f_minus_over_kz.request();
    const auto zrefi_info = z_reference_i.request();
    const auto zreff_info = z_reference_f.request();
    const auto interface_info = z_interfaces.request();
    const auto alphai_info = alpha_i.request();
    const auto alphaf_info = alpha_f.request();
    const auto cos_info = cos_azimuth.request();
    const auto sin_info = sin_azimuth.request();
    const auto reference_share_info = reference_delta_beta.request();
    const auto basis_info = basis.request();
    const auto factors_info = f_factors.request();
    const auto finite_position_info = finite_positions.request();
    const auto finite_atom_info = finite_atom_index.request();
    const auto finite_record_info = finite_record_index.request();
    const auto finite_weight_info = finite_weight.request();
    const auto bulk_position_info = bulk_positions.request();
    const auto bulk_atom_info = bulk_atom_index.request();
    const auto bulk_weight_info = bulk_weight.request();
    const auto bulk_coordinate_info = bulk_repeat_coordinate.request();
    const auto bulk_repeat_info = bulk_repeat.request();

    require_shape(Q_info, 2, "Q_parallel");
    require_shape(specular_info, 1, "is_specular");
    require_shape(reference_n_info, 1, "reference_refractive_indices");
    require_shape(field_i_index_info, 1, "field_i_index");
    require_shape(field_f_index_info, 1, "field_f_index");
    require_shape(kzi_info, 2, "kz_i");
    require_shape(Aip_info, 2, "A_i_plus");
    require_shape(Aim_info, 2, "A_i_minus");
    require_shape(Aipk_info, 2, "A_i_plus_over_kz");
    require_shape(Aimk_info, 2, "A_i_minus_over_kz");
    require_shape(kzf_info, 2, "kz_f");
    require_shape(Afp_info, 2, "A_f_plus");
    require_shape(Afm_info, 2, "A_f_minus");
    require_shape(Afpk_info, 2, "A_f_plus_over_kz");
    require_shape(Afmk_info, 2, "A_f_minus_over_kz");
    require_shape(zrefi_info, 1, "z_reference_i");
    require_shape(zreff_info, 1, "z_reference_f");
    require_shape(interface_info, 1, "z_interfaces");
    require_shape(alphai_info, 1, "alpha_i");
    require_shape(alphaf_info, 1, "alpha_f");
    require_shape(cos_info, 1, "cos_azimuth");
    require_shape(sin_info, 1, "sin_azimuth");
    require_shape(reference_share_info, 3, "reference_delta_beta");
    require_shape(basis_info, 2, "basis");
    require_shape(factors_info, 2, "f_factors");
    require_shape(finite_position_info, 2, "finite_positions");
    require_shape(finite_atom_info, 1, "finite_atom_index");
    require_shape(finite_record_info, 1, "finite_record_index");
    require_shape(finite_weight_info, 1, "finite_weight");
    require_shape(bulk_position_info, 2, "bulk_positions");
    require_shape(bulk_atom_info, 1, "bulk_atom_index");
    require_shape(bulk_weight_info, 1, "bulk_weight");
    require_shape(bulk_coordinate_info, 1, "bulk_repeat_coordinate");
    require_shape(bulk_repeat_info, 1, "bulk_repeat");

    const auto points = Q_info.shape[0];
    const auto media = kzi_info.shape[0];
    const auto incident_fields = kzi_info.shape[1];
    const auto exit_fields = kzf_info.shape[1];
    const auto records = reference_share_info.shape[0];
    const auto atoms = basis_info.shape[0];
    const auto finite_instances = finite_position_info.shape[0];
    const auto bulk_instances = bulk_position_info.shape[0];
    const bool homogeneous_bulk_only = records == 1;
    if (Q_info.shape[1] != 2 || records < 1 || media < 2) {
        throw py::value_error("DWBA geometry requires records and at least two media");
    }
    if (reference_n_info.shape[0] != media
        || reference_share_info.shape[1] != media
        || reference_share_info.shape[2] != 2
        || interface_info.shape[0] != media - 1
        || zrefi_info.shape[0] != media || zreff_info.shape[0] != media) {
        throw py::value_error("DWBA reference arrays have inconsistent shapes");
    }
    if (specular_info.shape[0] != points
        || field_i_index_info.shape[0] != points
        || field_f_index_info.shape[0] != points
        || alphai_info.shape[0] != points || alphaf_info.shape[0] != points
        || cos_info.shape[0] != points || sin_info.shape[0] != points) {
        throw py::value_error("DWBA point arrays have inconsistent shapes");
    }
    const auto same_incident_shape = [media, incident_fields](
        const py::buffer_info &info
    ) {
        return info.shape[0] == media && info.shape[1] == incident_fields;
    };
    const auto same_exit_shape = [media, exit_fields](
        const py::buffer_info &info
    ) {
        return info.shape[0] == media && info.shape[1] == exit_fields;
    };
    if (!same_incident_shape(Aip_info) || !same_incident_shape(Aim_info)
        || !same_incident_shape(Aipk_info) || !same_incident_shape(Aimk_info)
        || !same_exit_shape(Afp_info) || !same_exit_shape(Afm_info)
        || !same_exit_shape(Afpk_info) || !same_exit_shape(Afmk_info)) {
        throw py::value_error("DWBA field tables have inconsistent shapes");
    }
    if (basis_info.shape[1] < 7 || factors_info.shape[0] != atoms
        || factors_info.shape[1] < 13) {
        throw py::value_error("DWBA atom and factor arrays are inconsistent");
    }
    if (finite_position_info.shape[1] != 3
        || finite_atom_info.shape[0] != finite_instances
        || finite_record_info.shape[0] != finite_instances
        || finite_weight_info.shape[0] != finite_instances
        || bulk_position_info.shape[1] != 3
        || bulk_atom_info.shape[0] != bulk_instances
        || bulk_weight_info.shape[0] != bulk_instances
        || bulk_coordinate_info.shape[0] != bulk_instances
        || bulk_repeat_info.shape[0] != 3) {
        throw py::value_error("DWBA packed atomic arrays are inconsistent");
    }
    if (polarization_i != "s" && polarization_i != "p") {
        throw py::value_error("polarization_i must be 's' or 'p'");
    }
    if (polarization_f != "s" && polarization_f != "p") {
        throw py::value_error("polarization_f must be 's' or 'p'");
    }
    if (!std::isfinite(k0) || k0 <= 0.0
        || !std::isfinite(reference_area) || reference_area <= 0.0
        || !std::isfinite(attenuation) || attenuation < 0.0) {
        throw py::value_error("DWBA scales and attenuation must be finite and valid");
    }

    const auto *specular_data = static_cast<const std::int64_t *>(specular_info.ptr);
    const auto *field_i_index_data = static_cast<const std::int64_t *>(field_i_index_info.ptr);
    const auto *field_f_index_data = static_cast<const std::int64_t *>(field_f_index_info.ptr);
    const auto *finite_atom_data = static_cast<const std::int64_t *>(finite_atom_info.ptr);
    const auto *finite_record_data = static_cast<const std::int64_t *>(finite_record_info.ptr);
    const auto *bulk_atom_data = static_cast<const std::int64_t *>(bulk_atom_info.ptr);
    for (py::ssize_t point = 0; point < points; ++point) {
        if ((specular_data[point] != 0 && specular_data[point] != 1)
            || field_i_index_data[point] < 0
            || field_i_index_data[point] >= incident_fields
            || field_f_index_data[point] < 0
            || field_f_index_data[point] >= exit_fields) {
            throw py::value_error("DWBA point index arrays contain invalid values");
        }
    }
    for (py::ssize_t instance = 0; instance < finite_instances; ++instance) {
        if (finite_atom_data[instance] < 0 || finite_atom_data[instance] >= atoms
            || finite_record_data[instance] <= 0
            || finite_record_data[instance] >= records) {
            throw py::value_error("finite DWBA instance index is invalid");
        }
    }
    for (py::ssize_t instance = 0; instance < bulk_instances; ++instance) {
        if (bulk_atom_data[instance] < 0 || bulk_atom_data[instance] >= atoms) {
            throw py::value_error("bulk DWBA atom index is invalid");
        }
    }

    auto atomic = py::array_t<complex128>({records, points});
    auto reference = py::array_t<complex128>({records, points});
    auto *atomic_out = static_cast<complex128 *>(atomic.request().ptr);
    auto *reference_out = static_cast<complex128 *>(reference.request().ptr);
    std::fill(atomic_out, atomic_out + records * points, complex128(0.0, 0.0));
    std::fill(reference_out, reference_out + records * points, complex128(0.0, 0.0));

    const auto *Q_data = static_cast<const double *>(Q_info.ptr);
    const auto *reference_n_data = static_cast<const complex128 *>(reference_n_info.ptr);
    const auto *kzi_data = static_cast<const complex128 *>(kzi_info.ptr);
    const auto *Aip_data = static_cast<const complex128 *>(Aip_info.ptr);
    const auto *Aim_data = static_cast<const complex128 *>(Aim_info.ptr);
    const auto *Aipk_data = static_cast<const complex128 *>(Aipk_info.ptr);
    const auto *Aimk_data = static_cast<const complex128 *>(Aimk_info.ptr);
    const auto *kzf_data = static_cast<const complex128 *>(kzf_info.ptr);
    const auto *Afp_data = static_cast<const complex128 *>(Afp_info.ptr);
    const auto *Afm_data = static_cast<const complex128 *>(Afm_info.ptr);
    const auto *Afpk_data = static_cast<const complex128 *>(Afpk_info.ptr);
    const auto *Afmk_data = static_cast<const complex128 *>(Afmk_info.ptr);
    const auto *zrefi_data = static_cast<const double *>(zrefi_info.ptr);
    const auto *zreff_data = static_cast<const double *>(zreff_info.ptr);
    const auto *interface_data = static_cast<const double *>(interface_info.ptr);
    const auto *alphai_data = static_cast<const double *>(alphai_info.ptr);
    const auto *alphaf_data = static_cast<const double *>(alphaf_info.ptr);
    const auto *cos_data = static_cast<const double *>(cos_info.ptr);
    const auto *sin_data = static_cast<const double *>(sin_info.ptr);
    const auto *share_data = static_cast<const double *>(reference_share_info.ptr);
    const auto *basis_data = static_cast<const double *>(basis_info.ptr);
    const auto *factor_data = static_cast<const double *>(factors_info.ptr);
    const auto *finite_position_data = static_cast<const double *>(finite_position_info.ptr);
    const auto *finite_weight_data = static_cast<const double *>(finite_weight_info.ptr);
    const auto *bulk_position_data = static_cast<const double *>(bulk_position_info.ptr);
    const auto *bulk_weight_data = static_cast<const double *>(bulk_weight_info.ptr);
    const auto *bulk_coordinate_data = static_cast<const double *>(bulk_coordinate_info.ptr);
    const auto *repeat_data = static_cast<const double *>(bulk_repeat_info.ptr);
    for (py::ssize_t medium = 0; medium < media; ++medium) {
        if (!std::isfinite(reference_n_data[medium].real())
            || !std::isfinite(reference_n_data[medium].imag())) {
            throw py::value_error("reference refractive indices must be finite");
        }
    }
    if (!(repeat_data[2] > 0.0) || !std::isfinite(repeat_data[2])) {
        throw py::value_error("bulk repeat must point towards positive surface z");
    }

    constexpr double pi = 3.141592653589793238462643383279502884;
    constexpr double classical_electron_radius = 2.8179403262e-5;
    constexpr double dw_denominator = 16.0 * pi * pi;
    constexpr int sigma_i[4] = {1, 1, -1, -1};
    constexpr int sigma_f[4] = {1, -1, 1, -1};
    const bool p_i = polarization_i == "p";
    const bool p_f = polarization_f == "p";

    {
    py::gil_scoped_release release;
    const auto geometry_index = [points](
        const py::ssize_t medium,
        const int channel,
        const py::ssize_t point
    ) {
        return (medium * 4 + channel) * points + point;
    };
    const std::size_t geometry_size = static_cast<std::size_t>(media * 4 * points);
    std::vector<complex128> channel_Qz(geometry_size);
    std::vector<complex128> channel_Q2(geometry_size);
    std::vector<complex128> channel_coefficient(geometry_size);
    for (py::ssize_t medium = 0; medium < media; ++medium) {
        for (int channel = 0; channel < 4; ++channel) {
            for (py::ssize_t point = 0; point < points; ++point) {
                const auto incident_offset = medium * incident_fields
                    + field_i_index_data[point];
                const auto exit_offset = medium * exit_fields
                    + field_f_index_data[point];
                const auto index = geometry_index(medium, channel, point);
                const double Qx = Q_data[point * 2];
                const double Qy = Q_data[point * 2 + 1];
                const complex128 incident_amplitude = sigma_i[channel] > 0
                    ? Aip_data[incident_offset] : Aim_data[incident_offset];
                const complex128 exit_amplitude = sigma_f[channel] > 0
                    ? Afp_data[exit_offset] : Afm_data[exit_offset];
                const complex128 Qz = -(
                    static_cast<double>(sigma_i[channel]) * kzi_data[incident_offset]
                    + static_cast<double>(sigma_f[channel]) * kzf_data[exit_offset]
                );
                channel_Qz[index] = Qz;
                channel_Q2[index] = Qx * Qx + Qy * Qy + Qz * Qz;
                if (incident_amplitude == complex128(0.0, 0.0)
                    || exit_amplitude == complex128(0.0, 0.0)) {
                    channel_coefficient[index] = complex128(0.0, 0.0);
                    continue;
                }
                const double incident_sine = std::sin(alphai_data[point]);
                const double exit_sine = std::sin(alphaf_data[point]);
                complex128 vector_weight;
                if (!p_i && !p_f) {
                    vector_weight = incident_amplitude * exit_amplitude
                        * cos_data[point];
                } else if (!p_i && p_f) {
                    vector_weight = -incident_amplitude
                        * static_cast<double>(sigma_f[channel])
                        * exit_amplitude * exit_sine * sin_data[point];
                } else if (p_i && !p_f) {
                    vector_weight = -static_cast<double>(sigma_i[channel])
                        * incident_amplitude * incident_sine
                        * exit_amplitude * sin_data[point];
                } else {
                    const complex128 incident_over_kz = sigma_i[channel] > 0
                        ? Aipk_data[incident_offset] : Aimk_data[incident_offset];
                    const complex128 exit_over_kz = sigma_f[channel] > 0
                        ? Afpk_data[exit_offset] : Afmk_data[exit_offset];
                    const complex128 tangent = -static_cast<double>(sigma_i[channel])
                        * static_cast<double>(sigma_f[channel])
                        * incident_amplitude * exit_amplitude
                        * incident_sine * exit_sine * cos_data[point];
                    const complex128 incident_normal = -k0
                        * std::cos(alphai_data[point]) * incident_sine
                        * incident_over_kz;
                    const complex128 exit_normal = -k0
                        * std::cos(alphaf_data[point]) * exit_sine
                        * exit_over_kz;
                    vector_weight = tangent + incident_normal * exit_normal;
                }
                const complex128 reference_phase = std::exp(
                    complex128(0.0, 1.0) * (
                        static_cast<double>(sigma_i[channel])
                        * kzi_data[incident_offset] * zrefi_data[medium]
                        + static_cast<double>(sigma_f[channel])
                        * kzf_data[exit_offset] * zreff_data[medium]
                    )
                );
                channel_coefficient[index] = vector_weight * reference_phase;
            }
        }
    }

    struct FormFactorGroup {
        std::array<double, 13> factor_row;
        std::vector<std::shared_ptr<const std::vector<double>>> factors;
    };
    // Displacement and occupancy values are part of the exact state key.  This
    // lets the complete atomic factor be reused without merging refinable atoms.
    struct ScatteringState {
        std::size_t form_factor_group;
        std::array<double, 3> parameters;
    };
    std::vector<FormFactorGroup> form_factor_groups;
    std::vector<ScatteringState> scattering_states;
    std::vector<std::size_t> atom_scattering_state(static_cast<std::size_t>(atoms));
    for (py::ssize_t atom = 0; atom < atoms; ++atom) {
        std::array<double, 13> row;
        for (int column = 0; column < 13; ++column) {
            row[column] = get2(factor_data, factors_info, atom, column);
        }
        std::size_t form_factor_group = form_factor_groups.size();
        for (std::size_t candidate = 0; candidate < form_factor_groups.size(); ++candidate) {
            if (std::memcmp(
                form_factor_groups[candidate].factor_row.data(), row.data(), sizeof(row)
            ) == 0) {
                form_factor_group = candidate;
                break;
            }
        }
        if (form_factor_group == form_factor_groups.size()) {
            form_factor_groups.push_back(FormFactorGroup{
                row,
                std::vector<std::shared_ptr<const std::vector<double>>>(
                    static_cast<std::size_t>(media * 4)
                )
            });
        }
        const std::array<double, 3> parameters = {
            get2(basis_data, basis_info, atom, 4),
            get2(basis_data, basis_info, atom, 5),
            get2(basis_data, basis_info, atom, 6),
        };
        std::size_t scattering_state = scattering_states.size();
        for (std::size_t candidate = 0; candidate < scattering_states.size(); ++candidate) {
            if (scattering_states[candidate].form_factor_group == form_factor_group
                && std::memcmp(
                    scattering_states[candidate].parameters.data(),
                    parameters.data(),
                    sizeof(parameters)
                ) == 0) {
                scattering_state = candidate;
                break;
            }
        }
        if (scattering_state == scattering_states.size()) {
            scattering_states.push_back(ScatteringState{
                form_factor_group, parameters
            });
        }
        atom_scattering_state[static_cast<std::size_t>(atom)] = scattering_state;
    }
    for (auto &group : form_factor_groups) {
        for (py::ssize_t medium = 0; medium < media; ++medium) {
            for (int channel = 0; channel < 4; ++channel) {
                std::vector<double> complex_grid(1 + static_cast<std::size_t>(2 * points));
                complex_grid[0] = std::numeric_limits<double>::quiet_NaN();
                for (py::ssize_t point = 0; point < points; ++point) {
                    const complex128 value = channel_Q2[
                        geometry_index(medium, channel, point)
                    ];
                    complex_grid[1 + 2 * point] = value.real();
                    complex_grid[2 + 2 * point] = value.imag();
                }
                auto factors = form_factor_cache().find(complex_grid, group.factor_row);
                if (!factors) {
                    std::vector<double> calculated(static_cast<std::size_t>(2 * points));
                    for (py::ssize_t point = 0; point < points; ++point) {
                        const complex128 Q_squared = channel_Q2[
                            geometry_index(medium, channel, point)
                        ];
                        complex128 value(
                            group.factor_row[10] + group.factor_row[11],
                            group.factor_row[12]
                        );
                        for (int term = 0; term < 5; ++term) {
                            value += group.factor_row[term] * std::exp(
                                -group.factor_row[term + 5] * Q_squared
                            );
                        }
                        calculated[2 * point] = value.real();
                        calculated[2 * point + 1] = value.imag();
                    }
                    factors = form_factor_cache().insert_or_get(
                        complex_grid, group.factor_row, std::move(calculated)
                    );
                }
                group.factors[static_cast<std::size_t>(medium * 4 + channel)] = factors;
            }
        }
    }

    const auto medium_for_z = [media, interface_data](const double z) {
        py::ssize_t medium = 0;
        while (medium < media - 1 && z < interface_data[medium]) {
            ++medium;
        }
        return medium;
    };
    std::vector<py::ssize_t> finite_medium(static_cast<std::size_t>(finite_instances));
    for (py::ssize_t instance = 0; instance < finite_instances; ++instance) {
        finite_medium[static_cast<std::size_t>(instance)] = medium_for_z(
            finite_position_data[instance * 3 + 2]
        );
    }
    py::ssize_t terminal_repeat = 0;
    if (semi_infinite && bulk_instances > 0 && !homogeneous_bulk_only) {
        double maximum_z = -std::numeric_limits<double>::infinity();
        for (py::ssize_t instance = 0; instance < bulk_instances; ++instance) {
            maximum_z = std::max(maximum_z, bulk_position_data[instance * 3 + 2]);
        }
        const double terminal_boundary = interface_data[media - 2];
        if (maximum_z >= terminal_boundary) {
            terminal_repeat = static_cast<py::ssize_t>(
                std::floor((maximum_z - terminal_boundary) / repeat_data[2]) + 1.0
            );
        }
        if (terminal_repeat > 1000000) {
            throw std::domain_error("graded DWBA bulk requires too many explicit repeats");
        }
    }

    const py::ssize_t explicit_repeats = semi_infinite ? terminal_repeat : 1;
    constexpr std::size_t missing_pair = std::numeric_limits<std::size_t>::max();
    // Only state/medium pairs reached by a packed instance enter the scratch
    // table.  Its size is independent of the number of scan points and domains.
    struct ScatteringPair {
        std::size_t state;
        py::ssize_t medium;
    };
    std::vector<ScatteringPair> scattering_pairs;
    std::vector<std::size_t> pair_lookup(
        scattering_states.size() * static_cast<std::size_t>(media), missing_pair
    );
    const auto register_pair = [
        media, &pair_lookup, &scattering_pairs
    ](const std::size_t state, const py::ssize_t medium) {
        const auto key = state * static_cast<std::size_t>(media)
            + static_cast<std::size_t>(medium);
        auto &pair_index = pair_lookup[key];
        if (pair_index == missing_pair) {
            pair_index = scattering_pairs.size();
            scattering_pairs.push_back(ScatteringPair{state, medium});
        }
        return pair_index;
    };

    std::vector<std::size_t> finite_factor_index(
        static_cast<std::size_t>(finite_instances)
    );
    for (py::ssize_t instance = 0; instance < finite_instances; ++instance) {
        const auto atom = static_cast<std::size_t>(finite_atom_data[instance]);
        finite_factor_index[static_cast<std::size_t>(instance)] = register_pair(
            atom_scattering_state[atom],
            finite_medium[static_cast<std::size_t>(instance)]
        );
    }
    for (py::ssize_t repeat = 0; repeat < explicit_repeats; ++repeat) {
        for (py::ssize_t instance = 0; instance < bulk_instances; ++instance) {
            const double z = bulk_position_data[instance * 3 + 2]
                - repeat * repeat_data[2];
            const auto medium = homogeneous_bulk_only ? media - 1 : medium_for_z(z);
            const auto atom = static_cast<std::size_t>(bulk_atom_data[instance]);
            register_pair(atom_scattering_state[atom], medium);
        }
    }
    if (semi_infinite) {
        for (py::ssize_t instance = 0; instance < bulk_instances; ++instance) {
            const auto atom = static_cast<std::size_t>(bulk_atom_data[instance]);
            register_pair(atom_scattering_state[atom], media - 1);
        }
    }
    std::vector<complex128> atom_factors(
        scattering_pairs.size() * 4,
        complex128(0.0, 0.0)
    );

    for (py::ssize_t point = 0; point < points; ++point) {
        const double Qx = Q_data[point * 2];
        const double Qy = Q_data[point * 2 + 1];
        const double Q_parallel_squared = Qx * Qx + Qy * Qy;
        for (int channel = 0; channel < 4; ++channel) {
            for (
                std::size_t pair_index = 0;
                pair_index < scattering_pairs.size();
                ++pair_index
            ) {
                const auto pair = scattering_pairs[pair_index];
                const auto geometry = geometry_index(pair.medium, channel, point);
                if (channel_coefficient[geometry] == complex128(0.0, 0.0)) {
                    atom_factors[static_cast<std::size_t>(channel) * scattering_pairs.size()
                        + pair_index] = complex128(0.0, 0.0);
                    continue;
                }
                const auto &state = scattering_states[pair.state];
                const auto &cached = *form_factor_groups[
                    state.form_factor_group
                ].factors[static_cast<std::size_t>(pair.medium * 4 + channel)];
                complex128 factor(cached[2 * point], cached[2 * point + 1]);
                const complex128 Qz = channel_Qz[geometry];
                factor *= std::exp(-(
                    state.parameters[0] * Q_parallel_squared
                    + state.parameters[1] * Qz * Qz
                ) / dw_denominator);
                factor *= state.parameters[2];
                atom_factors[static_cast<std::size_t>(channel) * scattering_pairs.size()
                    + pair_index] = factor;
            }
            for (py::ssize_t instance = 0; instance < finite_instances; ++instance) {
                const auto medium = finite_medium[static_cast<std::size_t>(instance)];
                const auto geometry = geometry_index(medium, channel, point);
                const complex128 coefficient = channel_coefficient[geometry];
                if (coefficient == complex128(0.0, 0.0)) {
                    continue;
                }
                const complex128 Qz = channel_Qz[geometry];
                const double *position = finite_position_data + instance * 3;
                atomic_out[finite_record_data[instance] * points + point] +=
                    finite_weight_data[instance] * coefficient
                    * atom_factors[
                        static_cast<std::size_t>(channel) * scattering_pairs.size()
                        + finite_factor_index[static_cast<std::size_t>(instance)]
                    ]
                    * std::exp(complex128(0.0, 1.0) * (
                        Qx * position[0] + Qy * position[1] + Qz * position[2]
                    ));
            }
        }
        for (py::ssize_t repeat = 0; repeat < explicit_repeats; ++repeat) {
            for (py::ssize_t instance = 0; instance < bulk_instances; ++instance) {
                const double position[3] = {
                    bulk_position_data[instance * 3] - repeat * repeat_data[0],
                    bulk_position_data[instance * 3 + 1] - repeat * repeat_data[1],
                    bulk_position_data[instance * 3 + 2] - repeat * repeat_data[2],
                };
                const auto medium = homogeneous_bulk_only
                    ? media - 1 : medium_for_z(position[2]);
                const auto atom = static_cast<std::size_t>(bulk_atom_data[instance]);
                const auto factor_index = pair_lookup[
                    atom_scattering_state[atom] * static_cast<std::size_t>(media)
                    + static_cast<std::size_t>(medium)
                ];
                const double empirical = semi_infinite
                    ? std::exp(attenuation * (bulk_coordinate_data[instance] - repeat))
                    : 1.0;
                for (int channel = 0; channel < 4; ++channel) {
                    const auto geometry = geometry_index(medium, channel, point);
                    const complex128 coefficient = channel_coefficient[geometry];
                    if (coefficient == complex128(0.0, 0.0)) {
                        continue;
                    }
                    const complex128 Qz = channel_Qz[geometry];
                    atomic_out[point] += bulk_weight_data[instance] * empirical
                        * coefficient * atom_factors[
                            static_cast<std::size_t>(channel) * scattering_pairs.size()
                            + factor_index
                        ]
                        * std::exp(complex128(0.0, 1.0) * (
                            Qx * position[0] + Qy * position[1] + Qz * position[2]
                        ));
                }
            }
        }
        if (semi_infinite && bulk_instances > 0) {
            const py::ssize_t medium = media - 1;
            const int channel = 0;
            const auto geometry = geometry_index(medium, channel, point);
            complex128 numerator(0.0, 0.0);
            const complex128 Qz = channel_Qz[geometry];
            for (py::ssize_t instance = 0; instance < bulk_instances; ++instance) {
                const double position[3] = {
                    bulk_position_data[instance * 3]
                        - terminal_repeat * repeat_data[0],
                    bulk_position_data[instance * 3 + 1]
                        - terminal_repeat * repeat_data[1],
                    bulk_position_data[instance * 3 + 2]
                        - terminal_repeat * repeat_data[2],
                };
                const auto atom = static_cast<std::size_t>(bulk_atom_data[instance]);
                const auto factor_index = pair_lookup[
                    atom_scattering_state[atom] * static_cast<std::size_t>(media)
                    + static_cast<std::size_t>(medium)
                ];
                const double empirical = std::exp(
                    attenuation * (bulk_coordinate_data[instance] - terminal_repeat)
                );
                numerator += bulk_weight_data[instance] * empirical
                    * channel_coefficient[geometry]
                    * atom_factors[factor_index]
                    * std::exp(complex128(0.0, 1.0) * (
                        Qx * position[0] + Qy * position[1] + Qz * position[2]
                    ));
            }
            if (numerator != complex128(0.0, 0.0)) {
                const complex128 repeat_phase = Qx * repeat_data[0]
                    + Qy * repeat_data[1] + Qz * repeat_data[2];
                const complex128 denominator = -complex_expm1(
                    -complex128(0.0, 1.0) * repeat_phase - attenuation
                );
                const double pole_tolerance = 16.0
                    * std::numeric_limits<double>::epsilon()
                    * (1.0 + std::abs(repeat_phase));
                if (attenuation == 0.0 && std::abs(denominator) <= pole_tolerance) {
                    throw std::domain_error(
                        "zero-attenuation semi-infinite bulk has an exact Bragg pole"
                    );
                }
                atomic_out[point] += numerator / denominator;
            }
        }

        if (semi_infinite && specular_data[point] != 0) {
            for (py::ssize_t record = 0; record < records; ++record) {
                for (py::ssize_t medium = 1; medium < media; ++medium) {
                    const auto share_offset = (record * media + medium) * 2;
                    const complex128 density = (
                        k0 * k0 / (2.0 * pi * classical_electron_radius)
                    ) * complex128(
                        share_data[share_offset], share_data[share_offset + 1]
                    );
                    if (density == complex128(0.0, 0.0)) {
                        continue;
                    }
                    for (int channel = 0; channel < 4; ++channel) {
                        const auto geometry = geometry_index(medium, channel, point);
                        const complex128 coefficient = channel_coefficient[geometry];
                        if (coefficient == complex128(0.0, 0.0)) {
                            continue;
                        }
                        complex128 q = channel_Qz[geometry];
                        if (record == 0 && attenuation != 0.0) {
                            q -= complex128(0.0, attenuation / repeat_data[2]);
                        }
                        const complex128 iq = complex128(0.0, 1.0) * q;
                        complex128 phi;
                        if (medium == media - 1) {
                            if (iq == complex128(0.0, 0.0)) {
                                throw std::domain_error(
                                    "terminal reference integral does not converge"
                                );
                            }
                            const double upper = interface_data[media - 2];
                            phi = std::exp(iq * upper) / iq;
                        } else {
                            const double lower = interface_data[medium];
                            const double upper = interface_data[medium - 1];
                            const double thickness = upper - lower;
                            if (iq == complex128(0.0, 0.0)) {
                                phi = complex128(thickness, 0.0);
                            } else {
                                phi = std::exp(iq * lower)
                                    * complex_expm1(iq * thickness) / iq;
                            }
                        }
                        reference_out[record * points + point] +=
                            reference_area * density * coefficient * phi;
                    }
                }
            }
        }
    }
    }
    return py::make_tuple(atomic, reference);
}

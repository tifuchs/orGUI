#pragma once

// Model-conditioned counting moments, in raw detector counts squared.
// All donors in a connected repair neighborhood have ONE independent Poisson
// rate. This is not a model for spatial gradients or interpolation bias.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

struct MedianPoissonMoments {
    double mean;
    double variance;
    double donor_covariance;
};

// The finite-difference expansion is validated through 32 donors, including
// complete overlap. Larger nearly identical samples suffer cancellation in
// double precision; rejecting them is preferable to false counting precision.
inline constexpr int median_poisson_max_donors = 32;
inline constexpr double median_poisson_max_rate = 1.e8;
inline constexpr double median_poisson_table_max_rate = 1048576.;

inline bool median_poisson_supported(int n, int m, int shared) {
    return n >= 1 && m >= 1
        && n <= median_poisson_max_donors && m <= median_poisson_max_donors
        && shared >= 0 && shared <= std::min(n,m);
}

namespace roi_median_poisson_detail {

inline double choose(int n, int k) {
    if (k < 0 || k > n) return 0.;
    k = std::min(k, n-k);
    double value = 1.;
    for (int j = 1; j <= k; ++j) value *= double(n-k+j)/j;
    return value;
}

inline void require_arguments(int n, double rate) {
    if (n < 1 || n > median_poisson_max_donors)
        throw std::range_error("Poisson median moments support 1..32 donors");
    if (!std::isfinite(rate) || rate < 0. || rate > median_poisson_max_rate)
        throw std::range_error("Poisson median rate must be finite and in [0, 1e8]");
}

struct PoissonSupport {
    int first;
    std::vector<double> mass;
    std::vector<double> cdf;
    std::vector<double> survival;

    explicit PoissonSupport(double rate) {
        const double sigma = std::sqrt(rate);
        first = std::max(0, static_cast<int>(std::floor(rate-12.*sigma-32.)));
        const int last = static_cast<int>(std::ceil(rate+12.*sigma+64.));
        const int mode = static_cast<int>(std::floor(rate));
        const int length = last-first+1;
        mass.resize(length);
        cdf.resize(length);
        survival.resize(length);
        // Relative PMFs avoid cancellation in -lambda+k*log(lambda)-lgamma(k).
        // Normalization drops only the bounded far tails of this support.
        mass[mode-first] = 1.;
        for (int k = mode; k > first; --k)
            mass[k-1-first] = mass[k-first] * k/rate;
        for (int k = mode+1; k <= last; ++k)
            mass[k-first] = mass[k-1-first] * rate/k;
        double total = 0.;
        for (double value : mass) total += value;
        double partial = 0.;
        for (int j = 0; j < length; ++j) {
            mass[j] /= total;
            partial += mass[j];
            cdf[j] = std::min(1., partial);
        }
        partial = 0.;
        for (int j = length-1; j >= 0; --j) {
            survival[j] = partial;
            partial += mass[j];
        }
    }
};

struct DifferenceTerm {
    int left_power;
    int right_power;
    double coefficient;
};

inline std::vector<std::vector<DifferenceTerm>> difference_terms(int n, int shared) {
    // E[Delta^j 1{Binomial(n-j,p) <= rank-1}]. Even medians average
    // both middle order statistics, preserving Poisson ties exactly.
    std::vector<std::vector<DifferenceTerm>> result(shared+1);
    const int low = (n+1)/2;
    for (int j = 1; j <= shared; ++j) {
        std::vector<double> coefficients(n-j+1, 0.);
        for (int a = 0; a < j; ++a) {
            const double sign = (j-a)%2 ? -1. : 1.;
            const double factor = sign*choose(j-1,a)*(n%2 ? 1. : .5);
            for (int r = low; r <= low+(n%2 ? 0 : 1); ++r) {
                const int index = r-1-a;
                if (index >= 0 && index <= n-j)
                    coefficients[index] += factor*choose(n-j,index);
            }
        }
        for (int index = 0; index <= n-j; ++index)
            if (coefficients[index])
                result[j].push_back({index,n-j-index,coefficients[index]});
    }
    return result;
}

inline double evaluate_difference(const std::vector<DifferenceTerm>& terms,
                                 const std::vector<double>& p,
                                 const std::vector<double>& q) {
    double result = 0.;
    for (const auto& term : terms)
        result += term.coefficient*p[term.left_power]*q[term.right_power];
    return result;
}

}  // namespace roi_median_poisson_detail

// Deterministic evaluation, with support tails beyond 12 sigma plus padding.
// No fitted Gaussian/continuous-median approximation or simulation is used.
inline MedianPoissonMoments median_poisson_moments(int n, double rate) {
    using namespace roi_median_poisson_detail;
    require_arguments(n, rate);
    if (!rate) return {0., 0., 0.};
    if (n <= 2) return {rate, rate/n, rate/n};
    const PoissonSupport support(rate);
    double mean = 0., second = 0., donor = 0.;
    double gap_second = 0., gap_prefix = 0.;
    const int low = (n+1)/2;
    const double gap_factor = n%2 ? 0. : choose(n, n/2);
    std::vector<double> p_power(n+1), q_power(n+1), binomial(n+1);
    for (int j = 0; j <= n; ++j) binomial[j] = choose(n,j);
    const double donor_low = choose(n-1,low-1), donor_high = choose(n-1,low);
    for (std::size_t i = 0; i < support.mass.size(); ++i) {
        const double p = support.cdf[i], q = support.survival[i];
        p_power[0] = q_power[0] = 1.;
        for (int j = 1; j <= n; ++j) {
            p_power[j] = p_power[j-1]*p;
            q_power[j] = q_power[j-1]*q;
        }
        double tail = 0.;
        for (int j = 0; j < low; ++j)
            tail += binomial[j]*p_power[j]*q_power[n-j];
        if (!(n%2)) tail += .5*binomial[low]*p_power[low]*q_power[n-low];
        mean += tail;
        second += (2.*i+1.)*tail;
        double derivative = donor_low*p_power[low-1]*q_power[n-low];
        if (!(n%2)) derivative = .5*(derivative
            + donor_high*p_power[low]*q_power[n-low-1]);
        donor += support.mass[i]*derivative;
        if (!(n%2)) {
            const double left = p_power[n/2];
            const double right = q_power[n/2];
            gap_second += (2.*gap_prefix+left)*right;
            gap_prefix += left;
        }
    }
    const double variance = second-mean*mean-.25*gap_factor*gap_second;
    if (variance < -1.e-10*std::max(1., rate))
        throw std::runtime_error("Poisson median variance lost numerical accuracy");
    return {support.first+mean, std::max(0., variance), rate*donor};
}

inline double median_poisson_covariance(int n, int m, int shared, double rate) {
    using namespace roi_median_poisson_detail;
    require_arguments(n, rate);
    require_arguments(m, rate);
    if (shared < 0 || shared > std::min(n,m))
        throw std::invalid_argument("invalid shared median donor count");
    if (!shared || !rate) return 0.;
    if (n == m && shared == n) return median_poisson_moments(n, rate).variance;
    if (n <= 2) return double(shared)/n * median_poisson_moments(m, rate).donor_covariance;
    if (m <= 2) return double(shared)/m * median_poisson_moments(n, rate).donor_covariance;
    const PoissonSupport support(rate);
    std::vector<double> prefix_n(shared+1, 0.), prefix_m(shared+1, 0.);
    const auto terms_n = difference_terms(n,shared);
    const auto terms_m = n == m ? terms_n : difference_terms(m,shared);
    const int largest = std::max(n,m);
    std::vector<double> p_power(largest+1), q_power(largest+1), factors(shared+1);
    for (int j = 1; j <= shared; ++j) factors[j] = choose(shared,j);
    double covariance = 0.;
    for (std::size_t i = 0; i < support.mass.size(); ++i) {
        const double p = support.cdf[i], q = support.survival[i];
        p_power[0] = q_power[0] = 1.;
        for (int j = 1; j <= largest; ++j) {
            p_power[j] = p_power[j-1]*p;
            q_power[j] = q_power[j-1]*q;
        }
        for (int j = 1; j <= shared; ++j) {
            const double dn = evaluate_difference(terms_n[j],p_power,q_power);
            const double dm = n == m ? dn : evaluate_difference(terms_m[j],p_power,q_power);
            const double an = p_power[j]*dn, am = p_power[j]*dm;
            const double bn = q_power[j]*dn, bm = q_power[j]*dm;
            // Finite Hoeffding expansion of shared Bernoulli threshold
            // indicators. For x<y each term separates into F(x)^j D_nj(F(x))
            // and (1-F(y))^j D_mj(F(y)), allowing exact prefix summation.
            covariance += factors[j]*(prefix_n[j]*bm
                + prefix_m[j]*bn + an*bm);
            prefix_n[j] += an;
            prefix_m[j] += am;
        }
    }
    if (!std::isfinite(covariance) || covariance < -1.e-10*std::max(1.,rate))
        throw std::runtime_error("Poisson median covariance lost numerical accuracy");
    return std::max(0., covariance);
}

// Every moment type uses the SAME lambda grid. Linear interpolation therefore
// mixes whole covariance matrices and preserves positive semidefiniteness;
// directly counted donor variance lambda has the same exact interpolation.
inline const std::vector<double>& median_poisson_rate_grid() {
    static const std::vector<double> grid = [] {
        std::vector<double> result{0.};
        // Relative spacing resolves the small-rate lambda^ceil(n/2) onset;
        // a uniform low-rate grid badly overestimates rare-event variances.
        for (double rate = 1.e-8; rate < 4.; rate *= 1.005)
            result.push_back(rate);
        result.push_back(4.);
        for (double step = .03125, end = 16.; end <= 4096.; step *= 4., end *= 4.)
            for (double rate = result.back()+step; rate <= end; rate += step)
                result.push_back(rate);
        for (double rate = 8192.; rate <= median_poisson_table_max_rate; rate *= 2.)
            result.push_back(rate);
        return result;
    }();
    return grid;
}

// A previously unseen NaN donor subset must use these same-grid fallbacks,
// rather than mix exact-at-lambda entries with interpolated entries from other
// repair types in one joint covariance matrix. No persistent cache is mutated.
inline MedianPoissonMoments median_poisson_interpolated_moments(int n, double rate) {
    const auto& grid = median_poisson_rate_grid();
    if (rate < grid[1] || !std::isfinite(rate) || rate > grid.back())
        return median_poisson_moments(n,rate);
    const auto upper = std::lower_bound(grid.begin(),grid.end(),rate);
    if (*upper == rate) return median_poisson_moments(n,rate);
    const double low = *(upper-1), high = *upper;
    const double weight = (rate-low)/(high-low);
    const auto a = median_poisson_moments(n,low);
    const auto b = median_poisson_moments(n,high);
    return {a.mean+weight*(b.mean-a.mean),
        a.variance+weight*(b.variance-a.variance),
        a.donor_covariance+weight*(b.donor_covariance-a.donor_covariance)};
}

inline double median_poisson_interpolated_covariance(int n, int m, int shared,
                                                    double rate) {
    const auto& grid = median_poisson_rate_grid();
    if (rate < grid[1] || !std::isfinite(rate) || rate > grid.back())
        return median_poisson_covariance(n,m,shared,rate);
    const auto upper = std::lower_bound(grid.begin(),grid.end(),rate);
    if (*upper == rate) return median_poisson_covariance(n,m,shared,rate);
    const double low = *(upper-1), high = *upper;
    const double weight = (rate-low)/(high-low);
    const double a = median_poisson_covariance(n,m,shared,low);
    const double b = median_poisson_covariance(n,m,shared,high);
    return a+weight*(b-a);
}

class MedianPoissonMomentsTable {
public:
    explicit MedianPoissonMomentsTable(int n) : n_(n) {
        values_.reserve(median_poisson_rate_grid().size());
        for (double rate : median_poisson_rate_grid())
            values_.push_back(median_poisson_moments(n, rate));
    }

    MedianPoissonMoments moments(double rate) const {
        const auto& grid = median_poisson_rate_grid();
        if (rate < 0. || !std::isfinite(rate) || rate > grid.back()
            || (rate > 0. && rate < grid[1]))
            return median_poisson_moments(n_, rate);
        const auto upper = std::lower_bound(grid.begin(), grid.end(), rate);
        const std::size_t hi = static_cast<std::size_t>(upper-grid.begin());
        if (!hi || *upper == rate) return values_[hi];
        const std::size_t lo = hi-1;
        const double weight = (rate-grid[lo])/(grid[hi]-grid[lo]);
        const auto& a = values_[lo];
        const auto& b = values_[hi];
        return {a.mean+weight*(b.mean-a.mean),
            a.variance+weight*(b.variance-a.variance),
            a.donor_covariance+weight*(b.donor_covariance-a.donor_covariance)};
    }

    std::size_t storage_bytes() const { return sizeof(*this)+values_.capacity()*sizeof(MedianPoissonMoments); }
    int donors() const { return n_; }

private:
    int n_;
    std::vector<MedianPoissonMoments> values_;
};

class MedianPoissonTable {
public:
    MedianPoissonTable(int n, int m, int shared) : n_(n), m_(m), shared_(shared) {
        values_.reserve(median_poisson_rate_grid().size());
        for (double rate : median_poisson_rate_grid())
            values_.push_back(median_poisson_covariance(n,m,shared,rate));
    }

    double covariance(double rate) const {
        const auto& grid = median_poisson_rate_grid();
        if (rate < 0. || !std::isfinite(rate) || rate > grid.back()
            || (rate > 0. && rate < grid[1]))
            return median_poisson_covariance(n_,m_,shared_,rate);
        const auto upper = std::lower_bound(grid.begin(),grid.end(),rate);
        const std::size_t hi = static_cast<std::size_t>(upper-grid.begin());
        if (!hi || *upper == rate) return values_[hi];
        const std::size_t lo = hi-1;
        const double weight = (rate-grid[lo])/(grid[hi]-grid[lo]);
        return values_[lo]+weight*(values_[hi]-values_[lo]);
    }

    std::size_t storage_bytes() const { return sizeof(*this)+values_.capacity()*sizeof(double); }

private:
    int n_, m_, shared_;
    std::vector<double> values_;
};

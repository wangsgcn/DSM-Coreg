#include "dsm_coreg/alignment.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dsm_coreg {
namespace {

inline bool finite(double value) {
    return std::isfinite(value);
}

// Compact 3x3 matrix used by the normal-equation solver.  Keeping this local
// avoids introducing a large linear-algebra dependency for a problem whose
// parameter vector always has exactly three elements: tx, ty, and tz.
struct Mat3 {
    double a[9] = {0.0};

    double& operator()(int r, int c) { return a[r * 3 + c]; }
    double operator()(int r, int c) const { return a[r * 3 + c]; }
};

Mat3 zero_mat3() {
    return Mat3{};
}

// Gaussian elimination with row pivoting.  For a 3x3 system this is fast and
// sufficient, while still protecting against obviously singular/ill-conditioned
// normal equations.  A more general solver (Eigen/LAPACK) can be substituted
// later without changing the alignment interface.
bool solve_3x3(Mat3 matrix, const double rhs_input[3], double solution[3]) {
    double rhs[3] = {rhs_input[0], rhs_input[1], rhs_input[2]};
    int pivot[3] = {0, 1, 2};

    for (int k = 0; k < 3; ++k) {
        int best = k;
        double best_abs = std::abs(matrix(pivot[k], k));
        for (int r = k + 1; r < 3; ++r) {
            const double candidate = std::abs(matrix(pivot[r], k));
            if (candidate > best_abs) {
                best_abs = candidate;
                best = r;
            }
        }

        if (best_abs < 1.0e-18) return false;
        std::swap(pivot[k], pivot[best]);

        const int pk = pivot[k];
        const double diagonal = matrix(pk, k);

        for (int r = k + 1; r < 3; ++r) {
            const int pr = pivot[r];
            const double factor = matrix(pr, k) / diagonal;
            if (factor == 0.0) continue;

            for (int c = k; c < 3; ++c) {
                matrix(pr, c) -= factor * matrix(pk, c);
            }
            rhs[pr] -= factor * rhs[pk];
        }
    }

    for (int k = 2; k >= 0; --k) {
        const int pk = pivot[k];
        double sum = rhs[pk];
        for (int c = k + 1; c < 3; ++c) {
            sum -= matrix(pk, c) * solution[c];
        }

        const double diagonal = matrix(pk, k);
        if (std::abs(diagonal) < 1.0e-18) return false;
        solution[k] = sum / diagonal;
    }

    return true;
}

inline double huber_weight(double residual, double delta) {
    if (delta <= 0.0) return 1.0;
    const double abs_r = std::abs(residual);
    return (abs_r <= delta) ? 1.0 : (delta / abs_r);
}

// nth_element computes a median in linear expected time and avoids sorting the
// entire residual vector for every coarse-search candidate.
double median_inplace(std::vector<double>& values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();

    const std::size_t n = values.size();
    const std::size_t mid = n / 2;
    std::nth_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid), values.end());
    double median = values[mid];

    if (n % 2 == 0) {
        const auto max_lower = std::max_element(values.begin(), values.begin() + static_cast<std::ptrdiff_t>(mid));
        median = 0.5 * (median + *max_lower);
    }
    return median;
}

struct FitDiagnostics {
    double weighted_rmse = std::numeric_limits<double>::quiet_NaN();
    int used = 0;
};

FitDiagnostics evaluate_fit(const std::vector<Sample>& samples,
                            const std::vector<std::vector<int>>& block_to_sample_indices,
                            const std::vector<int>& block_counts,
                            const DemRaster& target,
                            const Translation& translation,
                            double huber_delta_m) {
    double weighted_rss = 0.0;
    double sum_weight = 0.0;
    int used = 0;

    auto block_count = [&](int b) -> int {
        return block_counts.empty() ? 1 : block_counts[static_cast<std::size_t>(b)];
    };

    for (int b = 0; b < static_cast<int>(block_to_sample_indices.size()); ++b) {
        const int multiplicity = block_count(b);
        if (multiplicity <= 0) continue;

        for (const int index : block_to_sample_indices[static_cast<std::size_t>(b)]) {
            const Sample& sample = samples[static_cast<std::size_t>(index)];
            const double x_sample = sample.x - translation.tx;
            const double y_sample = sample.y - translation.ty;

            double z_target = 0.0;
            if (!target.bilinear_sample(x_sample, y_sample, z_target)) continue;

            const double residual = sample.z_reference - (z_target + translation.tz);
            if (!finite(residual)) continue;

            double weight = huber_weight(residual, huber_delta_m);
            weight *= static_cast<double>(multiplicity);

            weighted_rss += weight * residual * residual;
            sum_weight += weight;
            ++used;
        }
    }

    FitDiagnostics out;
    out.used = used;
    if (sum_weight > 0.0) {
        out.weighted_rmse = std::sqrt(weighted_rss / sum_weight);
    }
    return out;
}

}  // namespace

SampleSet build_samples(const DemRaster& reference,
                        const DemRaster& target,
                        const Translation& current,
                        int border_px,
                        double max_abs_dz_m) {
    reference.validate_north_up("reference");
    target.validate_north_up("target");

    SampleSet output;
    const int border = std::max(0, border_px);
    const int inner_width = std::max(0, reference.width - 2 * border);
    const int inner_height = std::max(0, reference.height - 2 * border);
    output.samples.reserve(static_cast<std::size_t>(inner_width) * static_cast<std::size_t>(inner_height));

    for (int r = border; r < reference.height - border; ++r) {
        for (int c = border; c < reference.width - border; ++c) {
            const double z_reference = reference.at(c, r);
            if (!reference.valid_z(z_reference)) {
                ++output.skipped_reference_invalid;
                continue;
            }

            double x = 0.0;
            double y = 0.0;
            reference.pixel_center_to_world(c, r, x, y);

            // This is the important change relative to the legacy code: the
            // target is sampled using the current alignment estimate, not at
            // zero horizontal shift.  The pre-filter therefore removes genuine
            // surface disagreement rather than accidentally rejecting samples
            // solely because the DSMs have not been aligned yet.
            double z_target = 0.0;
            if (!target.bilinear_sample(x - current.tx, y - current.ty, z_target)) {
                ++output.skipped_target_invalid;
                continue;
            }

            if (max_abs_dz_m > 0.0) {
                const double residual = z_reference - (z_target + current.tz);
                if (!finite(residual) || std::abs(residual) > max_abs_dz_m) {
                    ++output.skipped_dz_filter;
                    continue;
                }
            }

            output.samples.push_back(Sample{x, y, z_reference, c, r});
        }
    }

    return output;
}

CoarseSearchResult coarse_xy_search(const DemRaster& reference,
                                    const DemRaster& target,
                                    const Translation& initial,
                                    const CoarseSearchConfig& cfg,
                                    Logger& logger) {
    CoarseSearchResult best;
    best.robust_score_m = std::numeric_limits<double>::infinity();

    // Candidate count is determined explicitly rather than through floating
    // point loop termination so repeated runs evaluate the same set of shifts.
    const int nx = static_cast<int>(std::floor((2.0 * cfg.range_x_m) / cfg.step_m + 0.5)) + 1;
    const int ny = static_cast<int>(std::floor((2.0 * cfg.range_y_m) / cfg.step_m + 0.5)) + 1;
    const int total_candidates = nx * ny;

    const double start_tx = initial.tx - cfg.range_x_m;
    const double start_ty = initial.ty - cfg.range_y_m;

    // Each candidate is independent, which makes this loop safe to parallelize.
    // The expensive bootstrap still dominates overall runtime for large jobs;
    // this search operates only on the coarsest DSM and a strided sample set.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int candidate_index = 0; candidate_index < total_candidates; ++candidate_index) {
        const int ix = candidate_index % nx;
        const int iy = candidate_index / nx;
        const double tx = start_tx + static_cast<double>(ix) * cfg.step_m;
        const double ty = start_ty + static_cast<double>(iy) * cfg.step_m;

        std::vector<double> dz;
        const std::size_t rough_capacity =
            static_cast<std::size_t>(std::max(1, reference.width / cfg.sample_stride)) *
            static_cast<std::size_t>(std::max(1, reference.height / cfg.sample_stride));
        dz.reserve(rough_capacity);

        for (int r = 0; r < reference.height; r += cfg.sample_stride) {
            for (int c = 0; c < reference.width; c += cfg.sample_stride) {
                const double z_reference = reference.at(c, r);
                if (!reference.valid_z(z_reference)) continue;

                double x = 0.0;
                double y = 0.0;
                reference.pixel_center_to_world(c, r, x, y);

                double z_target = 0.0;
                if (!target.bilinear_sample(x - tx, y - ty, z_target)) continue;

                const double difference = z_reference - z_target;
                if (finite(difference)) dz.push_back(difference);
            }
        }

        if (static_cast<int>(dz.size()) < cfg.min_samples) continue;

        // Given tx and ty, the robust vertical nuisance parameter is the median
        // elevation difference.  This prevents an unknown datum/height bias from
        // distorting the XY grid-search score.
        std::vector<double> dz_for_median = dz;
        const double tz = median_inplace(dz_for_median);

        std::vector<double> absolute_residuals;
        absolute_residuals.reserve(dz.size());
        for (const double value : dz) {
            absolute_residuals.push_back(std::abs(value - tz));
        }
        const double score = median_inplace(absolute_residuals);
        if (!finite(score)) continue;

#ifdef _OPENMP
#pragma omp critical(dsm_coreg_coarse_best)
#endif
        {
            if (score < best.robust_score_m) {
                best.translation = Translation{tx, ty, tz};
                best.robust_score_m = score;
                best.samples_used = static_cast<int>(dz.size());
                best.ok = true;
            }
        }
    }

    best.candidates_evaluated = total_candidates;

    if (best.ok) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(4)
            << "Coarse XY search selected tx=" << best.translation.tx
            << " m, ty=" << best.translation.ty
            << " m, tz=" << best.translation.tz
            << " m; median absolute residual=" << best.robust_score_m
            << " m; samples=" << best.samples_used
            << "; candidates=" << best.candidates_evaluated;
        logger.info(oss.str());
    } else {
        logger.warning("Coarse XY search did not find a candidate with sufficient overlap.");
    }

    return best;
}

std::vector<std::vector<int>> make_single_block(const std::vector<Sample>& samples) {
    std::vector<std::vector<int>> blocks(1);
    blocks[0].resize(samples.size());
    std::iota(blocks[0].begin(), blocks[0].end(), 0);
    return blocks;
}

FitResult gauss_newton_fit(
    const std::vector<Sample>& samples,
    const std::vector<std::vector<int>>& block_to_sample_indices,
    const std::vector<int>& block_counts,
    const DemRaster& target,
    const GradRaster& gradient,
    const Translation& initial,
    const AlignmentConfig& cfg,
    bool debug_iterations,
    Logger* logger) {

    FitResult output;
    Translation current = initial;

    auto block_count = [&](int b) -> int {
        return block_counts.empty() ? 1 : block_counts[static_cast<std::size_t>(b)];
    };

    for (int iteration = 0; iteration < cfg.max_iterations; ++iteration) {
        Mat3 hessian = zero_mat3();
        double gradient_vector[3] = {0.0, 0.0, 0.0};
        int sample_used = 0;

        for (int b = 0; b < static_cast<int>(block_to_sample_indices.size()); ++b) {
            const int multiplicity = block_count(b);
            if (multiplicity <= 0) continue;

            for (const int index : block_to_sample_indices[static_cast<std::size_t>(b)]) {
                const Sample& sample = samples[static_cast<std::size_t>(index)];
                const double x_sample = sample.x - current.tx;
                const double y_sample = sample.y - current.ty;

                double z_target = 0.0;
                double dz_dx = 0.0;
                double dz_dy = 0.0;
                if (!target.bilinear_sample(x_sample, y_sample, z_target)) continue;
                if (!gradient.bilinear_sample(x_sample, y_sample, dz_dx, dz_dy)) continue;

                const double residual = sample.z_reference - (z_target + current.tz);
                if (!finite(residual) || !finite(dz_dx) || !finite(dz_dy)) continue;

                // Huber IRLS weight.  Multiplying by the bootstrap block
                // multiplicity is equivalent to repeating all samples in the
                // block the selected number of times, without physically
                // duplicating them in memory.
                double weight = huber_weight(residual, cfg.huber_delta_m);
                weight *= static_cast<double>(multiplicity);

                // For r = z_ref - [z_tgt(x-tx,y-ty) + tz], the Jacobian is
                // [ +dz/dx, +dz/dy, -1 ].  The positive signs for tx/ty follow
                // from the chain rule applied to (x - tx, y - ty).
                const double jacobian[3] = {dz_dx, dz_dy, -1.0};

                for (int row = 0; row < 3; ++row) {
                    for (int col = 0; col < 3; ++col) {
                        hessian(row, col) += weight * jacobian[row] * jacobian[col];
                    }
                    gradient_vector[row] += weight * jacobian[row] * residual;
                }

                ++sample_used;
            }
        }

        if (sample_used < cfg.min_samples) {
            output.ok = false;
            return output;
        }

        if (cfg.lambda > 0.0) {
            hessian(0, 0) += cfg.lambda;
            hessian(1, 1) += cfg.lambda;
            hessian(2, 2) += cfg.lambda;
        }

        const double rhs[3] = {-gradient_vector[0], -gradient_vector[1], -gradient_vector[2]};
        double update[3] = {0.0, 0.0, 0.0};
        if (!solve_3x3(hessian, rhs, update)) {
            output.ok = false;
            return output;
        }

        current.tx += update[0];
        current.ty += update[1];
        current.tz += update[2];

        const double step = std::sqrt(update[0] * update[0] +
                                      update[1] * update[1] +
                                      update[2] * update[2]);

        if (debug_iterations && logger != nullptr) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(6)
                << "GN iter=" << iteration
                << " used=" << sample_used
                << " tx=" << current.tx
                << " ty=" << current.ty
                << " tz=" << current.tz
                << " step=" << step;
            logger->debug(oss.str());
        }

        output.iterations = iteration + 1;
        if (step < cfg.step_tolerance_m) {
            output.converged = true;
            break;
        }
    }

    // Recompute diagnostics at the FINAL parameter vector.  The legacy program
    // reported an RMSE accumulated before the final parameter update; the
    // difference is usually tiny near convergence, but recomputation is cleaner
    // and is important when termination occurs at the maximum iteration count.
    const FitDiagnostics diagnostics = evaluate_fit(samples,
                                                    block_to_sample_indices,
                                                    block_counts,
                                                    target,
                                                    current,
                                                    cfg.huber_delta_m);

    output.translation = current;
    output.weighted_rmse_m = diagnostics.weighted_rmse;
    output.samples_used = diagnostics.used;
    output.ok = diagnostics.used >= cfg.min_samples &&
                finite(current.tx) && finite(current.ty) && finite(current.tz);

    // Reaching max_iterations without satisfying the strict step tolerance does
    // not automatically invalidate a numerically stable result.  The separate
    // 'converged' flag is written to JSON so downstream users can distinguish
    // the two cases.
    return output;
}

}  // namespace dsm_coreg

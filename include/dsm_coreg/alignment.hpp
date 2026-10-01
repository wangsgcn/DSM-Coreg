#pragma once

#include "dsm_coreg/config.hpp"
#include "dsm_coreg/logging.hpp"
#include "dsm_coreg/raster.hpp"

#include <cstddef>
#include <vector>

namespace dsm_coreg {

struct Translation {
    double tx = 0.0;
    double ty = 0.0;
    double tz = 0.0;
};

struct Sample {
    double x = 0.0;
    double y = 0.0;
    double z_reference = 0.0;
    int ref_c = -1;
    int ref_r = -1;
};

struct SampleSet {
    std::vector<Sample> samples;
    std::size_t skipped_reference_invalid = 0;
    std::size_t skipped_target_invalid = 0;
    std::size_t skipped_dz_filter = 0;
};

struct FitResult {
    Translation translation;
    double weighted_rmse_m = 0.0;
    int samples_used = 0;
    int iterations = 0;
    bool converged = false;
    bool ok = false;
};

struct CoarseSearchResult {
    Translation translation;
    double robust_score_m = 0.0;
    int samples_used = 0;
    int candidates_evaluated = 0;
    bool ok = false;
};

struct AlignmentLevelRecord {
    int pyramid_factor = 1;
    double reference_pixel_size_x_m = 0.0;
    double reference_pixel_size_y_m = 0.0;
    std::size_t candidate_samples = 0;
    FitResult fit;
};

// Build reference samples after evaluating the target DSM at the CURRENT
// translation estimate.  This is a deliberate improvement over the legacy
// implementation, which applied max_abs_dz at zero shift and could therefore
// discard exactly the building-edge samples needed to recover a larger offset.
SampleSet build_samples(const DemRaster& reference,
                        const DemRaster& target,
                        const Translation& current,
                        int border_px,
                        double max_abs_dz_m);

// A robust global XY search used only when explicitly enabled.  For each XY
// candidate, tz is estimated by the median reference-target elevation
// difference and the candidate is scored by median absolute residual.  The
// method is intentionally simple and robust to construction/vegetation changes.
CoarseSearchResult coarse_xy_search(const DemRaster& reference,
                                    const DemRaster& target,
                                    const Translation& initial,
                                    const CoarseSearchConfig& cfg,
                                    Logger& logger);

// Gauss-Newton/Huber solver.  block_to_sample_indices and block_counts support
// the spatial bootstrap efficiently: each bootstrap replicate changes only the
// multiplicity of whole blocks, not the underlying sample vector.
FitResult gauss_newton_fit(
    const std::vector<Sample>& samples,
    const std::vector<std::vector<int>>& block_to_sample_indices,
    const std::vector<int>& block_counts,
    const DemRaster& target,
    const GradRaster& gradient,
    const Translation& initial,
    const AlignmentConfig& cfg,
    bool debug_iterations,
    Logger* logger);

// Utility for ordinary (non-bootstrap) fitting: place all samples in one block.
std::vector<std::vector<int>> make_single_block(const std::vector<Sample>& samples);

}  // namespace dsm_coreg

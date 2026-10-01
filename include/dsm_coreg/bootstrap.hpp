#pragma once

#include "dsm_coreg/alignment.hpp"
#include "dsm_coreg/config.hpp"
#include "dsm_coreg/logging.hpp"
#include "dsm_coreg/raster.hpp"
#include "dsm_coreg/timing.hpp"

#include <cstdint>
#include <filesystem>
#include <map>

namespace dsm_coreg {

struct Stats3 {
    double mean[3] = {0.0, 0.0, 0.0};
    double sd[3] = {0.0, 0.0, 0.0};
    double p05[3] = {0.0, 0.0, 0.0};
    double p95[3] = {0.0, 0.0, 0.0};
};

struct BootstrapResult {
    Stats3 stats;
    int requested_replicates = 0;
    int successful_replicates = 0;
    int failed_replicates = 0;

    // block size -> number of bootstrap replicates that used that size
    std::map<int, std::uint32_t> block_size_usage;

    // Optional audit/debug artifact.  When enabled, the GeoPackage contains
    // one polygon layer per bootstrap replicate plus a non-spatial replicate
    // summary table.  The translations in that summary are global replicate
    // estimates, not local per-block fits.
    bool geopackage_written = false;
    std::filesystem::path geopackage_path;

    bool ok = false;
};

// Run the spatial block bootstrap on the final/native-resolution sample set.
//
// reference is required only for optional GeoPackage audit output: its native
// grid, geotransform, and CRS define the polygons representing selected blocks.
// It is not used by the statistical fit itself.
//
// timing receives fine-grained timing for partition construction, replicate
// fitting, and optional GeoPackage serialization in addition to the enclosing
// spatial_block_bootstrap timer maintained by main.cpp.
BootstrapResult run_block_bootstrap(const std::vector<Sample>& samples,
                                    const DemRaster& reference,
                                    const DemRaster& target,
                                    const GradRaster& gradient,
                                    const FitResult& base_fit,
                                    const AlignmentConfig& alignment_cfg,
                                    const BootstrapConfig& bootstrap_cfg,
                                    bool debug_iterations,
                                    Logger& logger,
                                    TimingRegistry& timing);

}  // namespace dsm_coreg

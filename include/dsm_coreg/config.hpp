#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace dsm_coreg {

// -----------------------------------------------------------------------------
// Configuration model
// -----------------------------------------------------------------------------
//
// The application uses TOML as its human-edited input configuration format.
// TOML is a better fit than JSON for this purpose because experiment files can
// contain comments and are easy to scan by eye.  The program still writes its
// scientific results as JSON because JSON is convenient for downstream Python,
// C++, reporting, and batch-processing workflows.
//
// The configuration file is intentionally the authoritative experiment record.
// The command line remains small (normally only --config <file.toml>) so that a
// run can be reproduced later without reconstructing a long sequence of CLI
// switches.
// -----------------------------------------------------------------------------

struct InputConfig {
    std::filesystem::path reference_dsm;
    std::filesystem::path target_dsm;
};

struct AlignmentConfig {
    // Initial translation estimate, in map/elevation units (normally meters).
    // Positive tx/ty means the target DSM is translated in the positive map-x
    // / map-y direction to align it with the reference DSM.
    double initial_tx_m = 0.0;
    double initial_ty_m = 0.0;
    double initial_tz_m = 0.0;

    // Number of reference pixels skipped along each border. At coarse pyramid
    // levels this value is scaled down so that the skipped physical distance is
    // approximately comparable to the native-resolution setting.
    int border_px = 5;

    // Optional robust pre-filter. Samples are retained when the elevation
    // difference evaluated AT THE CURRENT ALIGNMENT ESTIMATE is no larger than
    // this value. A non-positive value disables the filter.
    double max_abs_dz_m = 50.0;

    // Huber tuning threshold (meters). Set <= 0 to recover ordinary least
    // squares. In heterogeneous urban DSMs this should normally remain > 0.
    double huber_delta_m = 1.0;

    int max_iterations = 30;

    // Convergence tolerance on the Euclidean norm of [dtx, dty, dtz].
    double step_tolerance_m = 1.0e-3;

    // Small diagonal damping term. This behaves like a simple
    // Levenberg-Marquardt stabilization of the normal equations.
    double lambda = 1.0e-3;

    // Minimum number of samples that must participate in a fit.
    int min_samples = 200;
};

struct PyramidConfig {
    bool enabled = true;

    // Downsampling factors, processed from coarse to fine. A value of 1 means
    // native resolution. The default 8-4-2-1 sequence increases the capture
    // range of Gauss-Newton without forcing the optional global search.
    std::vector<int> factors{8, 4, 2, 1};
};

struct CoarseSearchConfig {
    // Optional because the existing Gauss-Newton/Huber solver already converges
    // reliably for the user's observed several-meter offsets. Enable this only
    // when the initial geolocation may be substantially worse.
    bool enabled = false;

    // Search interval centered on the configured initial tx/ty.
    double range_x_m = 20.0;
    double range_y_m = 20.0;
    double step_m = 2.0;

    // Evaluate only every N-th reference pixel at the coarsest pyramid level.
    // This keeps the optional global search inexpensive relative to the full-
    // resolution optimization and bootstrap.
    int sample_stride = 4;

    // Reject candidate shifts with too little overlap.
    int min_samples = 200;
};

struct BootstrapOutputConfig {
    // Optional spatial audit product.  When enabled, dsm_coreg writes a
    // GeoPackage containing exactly one polygon feature layer per requested
    // bootstrap replicate.  Each layer shows which blocks were selected and
    // how many times they were drawn.
    bool write_geopackage = false;

    // Relative paths are resolved against output.directory after the complete
    // TOML document has been parsed.  An absolute path is also accepted.
    std::filesystem::path geopackage = "bootstrap_replicates.gpkg";
};

struct BootstrapConfig {
    bool enabled = true;
    int replicates = 200;

    // "fixed" or "multiscale".
    std::string block_mode = "multiscale";

    // Used when block_mode == "fixed".
    int block_px = 64;

    // Used when block_mode == "multiscale". One block size is selected for
    // each bootstrap replicate. Results are pooled, and usage counts are
    // written to the results JSON for transparency.
    int block_min_px = 32;
    int block_max_px = 128;
    int block_step_px = 16;

    std::uint64_t seed = 12345;

    // If <= 0, the program automatically requires at least max(20, N/5)
    // successful replicates.
    int min_successful_replicates = 0;

    BootstrapOutputConfig output;
};

struct RuntimeConfig {
    // 0 means use the OpenMP runtime default.
    int threads = 0;

    // By default the application writes diagnostics to the log file and keeps
    // stdout quiet. Setting this true prints only a compact final summary.
    bool console_summary = false;

    // When true, iteration-level diagnostics are written to the log file.
    bool debug_iterations = false;
};

struct OutputConfig {
    std::filesystem::path directory = "alignment_results";
    std::string results_json = "alignment_results.json";
    std::string log_file = "alignment.log";

    // The effective configuration is written as TOML so that it can be reused
    // directly as the starting point for a later experiment.
    std::string effective_config_toml = "effective_config.toml";

    bool write_aligned_dsm = true;
    std::string aligned_dsm = "target_aligned.tif";
};

struct ToolConfig {
    int config_version = 1;
    InputConfig input;
    AlignmentConfig alignment;
    PyramidConfig pyramid;
    CoarseSearchConfig coarse_search;
    BootstrapConfig bootstrap;
    RuntimeConfig runtime;
    OutputConfig output;

    // Directory containing the configuration file. Relative input/output paths
    // are resolved against this directory, not against the process working
    // directory. This makes archived experiment folders portable.
    std::filesystem::path config_directory;
};

// Parse, type-check, and semantically validate a TOML configuration file.
// Unknown keys are rejected deliberately: a misspelled scientific parameter
// should fail loudly instead of silently falling back to a default.
ToolConfig load_config(const std::filesystem::path& config_path);

// Write the fully resolved/effective configuration, including default values,
// so every run has a complete provenance record next to its result JSON.
void write_effective_config_toml(const ToolConfig& cfg,
                                 const std::filesystem::path& output_path);

}  // namespace dsm_coreg

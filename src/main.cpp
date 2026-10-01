#include "dsm_coreg/alignment.hpp"
#include "dsm_coreg/bootstrap.hpp"
#include "dsm_coreg/config.hpp"
#include "dsm_coreg/logging.hpp"
#include "dsm_coreg/raster.hpp"
#include "dsm_coreg/results.hpp"
#include "dsm_coreg/timing.hpp"

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef DSM_COREG_VERSION
#define DSM_COREG_VERSION "unknown"
#endif

namespace {

struct CliOptions {
    std::filesystem::path config_path;
    bool validate_only = false;
};

void print_usage(const char* executable) {
    std::cout
        << "Robust DSM co-registration with spatial block-bootstrap uncertainty\n\n"
        << "Usage:\n"
        << "  " << executable << " --config <config.toml> [--validate-only]\n"
        << "  " << executable << " --version\n\n"
        << "Options:\n"
        << "  --config FILE       TOML configuration file.\n"
        << "  --validate-only     Validate configuration, input metadata, and horizontal CRS;\n"
        << "                      do not run alignment or bootstrap.\n"
        << "  --version           Print program version.\n"
        << "  --help              Show this help text.\n";
}

CliOptions parse_cli(int argc, char** argv) {
    CliOptions options;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        }
        if (arg == "--version") {
            std::cout << "dsm_coreg " << DSM_COREG_VERSION << '\n';
            std::exit(0);
        }
        if (arg == "--validate-only") {
            options.validate_only = true;
            continue;
        }
        if (arg == "--config") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--config requires a file path.");
            }
            options.config_path = argv[++i];
            continue;
        }

        throw std::runtime_error("Unknown command-line argument: " + arg);
    }

    if (options.config_path.empty()) {
        throw std::runtime_error("A configuration file is required. Use --config <config.toml>.");
    }
    return options;
}

void validate_north_up_metadata(const dsm_coreg::DemMetadata& metadata,
                                const std::string& name) {
    if (std::abs(metadata.gt[2]) > 1.0e-12 || std::abs(metadata.gt[4]) > 1.0e-12) {
        throw std::runtime_error(name +
                                 " uses a rotated geotransform. Warp it to a north-up grid before alignment.");
    }
    if (metadata.gt[1] == 0.0 || metadata.gt[5] == 0.0) {
        throw std::runtime_error(name + " has an invalid zero pixel size.");
    }
}

std::string translation_string(const dsm_coreg::Translation& t) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4)
        << "tx=" << t.tx << " m, ty=" << t.ty << " m, tz=" << t.tz
        << " m, horizontal=" << std::hypot(t.tx, t.ty) << " m";
    return oss.str();
}

}  // namespace

int main(int argc, char** argv) {
    using namespace dsm_coreg;
    using Clock = std::chrono::steady_clock;

    const auto program_start = Clock::now();
    TimingRegistry timing;
    std::unique_ptr<Logger> logger;

    try {
        const CliOptions cli = parse_cli(argc, argv);

        ToolConfig cfg;
        {
            ScopedTimer timer(timing, "config_load");
            cfg = load_config(cli.config_path);
        }

        std::filesystem::create_directories(cfg.output.directory);
        logger = std::make_unique<Logger>(cfg.output.directory / cfg.output.log_file);
        logger->info(std::string("dsm_coreg version ") + DSM_COREG_VERSION + " started.");
        logger->info("Configuration: " + std::filesystem::absolute(cli.config_path).string());

        GDALAllRegister();
        OGRRegisterAll();

#ifdef _OPENMP
        if (cfg.runtime.threads > 0) omp_set_num_threads(cfg.runtime.threads);
        logger->info("OpenMP enabled; requested threads=" + std::to_string(cfg.runtime.threads));
#else
        if (cfg.runtime.threads > 0) {
            logger->warning("runtime.threads was set, but this binary was built without OpenMP.");
        }
        logger->info("OpenMP disabled at compile time.");
#endif

        // Save the resolved configuration immediately.  If a long alignment is
        // interrupted later, the output directory still records exactly which
        // inputs and parameters were intended for the run.
        write_effective_config_toml(cfg, cfg.output.directory / cfg.output.effective_config_toml);

        DemMetadata reference_metadata;
        DemMetadata target_metadata;
        {
            ScopedTimer timer(timing, "metadata_and_crs_validation");
            reference_metadata = read_dem_metadata(cfg.input.reference_dsm);
            target_metadata = read_dem_metadata(cfg.input.target_dsm);
            validate_north_up_metadata(reference_metadata, "Reference DSM");
            validate_north_up_metadata(target_metadata, "Target DSM");
            validate_horizontal_crs_compatible(reference_metadata.projection_wkt,
                                               target_metadata.projection_wkt);
        }
        logger->info("Input metadata and horizontal CRS validation passed.");

        if (cli.validate_only) {
            logger->info("Validation-only mode completed successfully; no alignment was run.");
            if (cfg.runtime.console_summary) {
                std::cout << "Configuration and DSM metadata validation succeeded.\n"
                          << "Log: " << (cfg.output.directory / cfg.output.log_file) << '\n';
            }
            return 0;
        }

        DemRaster reference;
        DemRaster target;
        {
            ScopedTimer timer(timing, "read_reference_dsm");
            reference = read_dem_geotiff(cfg.input.reference_dsm);
        }
        {
            ScopedTimer timer(timing, "read_target_dsm");
            target = read_dem_geotiff(cfg.input.target_dsm);
        }

        {
            std::ostringstream oss;
            oss << "Reference DSM: " << reference.width << 'x' << reference.height
                << ", pixel=" << std::abs(reference.gt[1]) << " x "
                << std::abs(reference.gt[5]) << " m";
            logger->info(oss.str());
        }
        {
            std::ostringstream oss;
            oss << "Target DSM: " << target.width << 'x' << target.height
                << ", pixel=" << std::abs(target.gt[1]) << " x "
                << std::abs(target.gt[5]) << " m";
            logger->info(oss.str());
        }

        std::vector<int> factors = cfg.pyramid.enabled ? cfg.pyramid.factors
                                                       : std::vector<int>{1};

        Translation current{cfg.alignment.initial_tx_m,
                            cfg.alignment.initial_ty_m,
                            cfg.alignment.initial_tz_m};
        RunResults run_results;

        // These are retained only for the native-resolution level because the
        // bootstrap must operate on the final sample set/gradient field.  Coarse
        // pyramid products are released after each level to limit memory use.
        std::optional<SampleSet> native_samples;
        std::optional<GradRaster> native_gradient;
        FitResult base_fit;

        {
            ScopedTimer alignment_timer(timing, "alignment_total");

            for (std::size_t level_index = 0; level_index < factors.size(); ++level_index) {
                const int factor = factors[level_index];

                std::optional<DemRaster> coarse_reference;
                std::optional<DemRaster> coarse_target;
                const DemRaster* reference_level = &reference;
                const DemRaster* target_level = &target;

                if (factor > 1) {
                    {
                        ScopedTimer timer(timing, "pyramid_build_factor_" + std::to_string(factor));
                        coarse_reference = downsample_average(reference, factor);
                        coarse_target = downsample_average(target, factor);
                    }
                    reference_level = &(*coarse_reference);
                    target_level = &(*coarse_target);
                }

                // Optional global initialization is performed only once at the
                // coarsest level.  Subsequent levels inherit the previous
                // estimate and therefore need only local Gauss-Newton refinement.
                if (level_index == 0 && cfg.coarse_search.enabled) {
                    ScopedTimer timer(timing, "coarse_xy_search");
                    CoarseSearchResult coarse = coarse_xy_search(*reference_level,
                                                                 *target_level,
                                                                 current,
                                                                 cfg.coarse_search,
                                                                 *logger);
                    run_results.coarse_search = coarse;
                    if (!coarse.ok) {
                        throw std::runtime_error("Optional coarse XY search failed to find a valid initialization.");
                    }
                    current = coarse.translation;
                }

                GradRaster gradient;
                {
                    ScopedTimer timer(timing, "gradient_factor_" + std::to_string(factor));
                    gradient = compute_target_gradients(*target_level);
                }

                // Scale the native-pixel border to the pyramid level while
                // preserving a similar physical exclusion width.  ceil() avoids
                // dropping a requested non-zero border to zero at coarse scales.
                const int level_border =
                    (cfg.alignment.border_px == 0)
                        ? 0
                        : std::max(1, static_cast<int>(std::ceil(
                                          static_cast<double>(cfg.alignment.border_px) /
                                          static_cast<double>(factor))));

                SampleSet sample_set;
                {
                    ScopedTimer timer(timing, "sample_build_factor_" + std::to_string(factor));
                    sample_set = build_samples(*reference_level,
                                               *target_level,
                                               current,
                                               level_border,
                                               cfg.alignment.max_abs_dz_m);
                }

                if (sample_set.samples.size() <
                    static_cast<std::size_t>(cfg.alignment.min_samples)) {
                    throw std::runtime_error(
                        "Too few samples at pyramid factor " + std::to_string(factor) +
                        ": " + std::to_string(sample_set.samples.size()) +
                        ". Consider increasing alignment.max_abs_dz_m or reviewing overlap/CRS.");
                }

                const std::vector<std::vector<int>> single_block =
                    make_single_block(sample_set.samples);
                const std::vector<int> empty_counts;

                FitResult fit;
                {
                    ScopedTimer timer(timing, "alignment_factor_" + std::to_string(factor));
                    fit = gauss_newton_fit(sample_set.samples,
                                           single_block,
                                           empty_counts,
                                           *target_level,
                                           gradient,
                                           current,
                                           cfg.alignment,
                                           cfg.runtime.debug_iterations,
                                           logger.get());
                }

                if (!fit.ok) {
                    throw std::runtime_error("Gauss-Newton fit failed at pyramid factor " +
                                             std::to_string(factor) + ".");
                }

                current = fit.translation;
                base_fit = fit;

                AlignmentLevelRecord record;
                record.pyramid_factor = factor;
                record.reference_pixel_size_x_m = std::abs(reference_level->gt[1]);
                record.reference_pixel_size_y_m = std::abs(reference_level->gt[5]);
                record.candidate_samples = sample_set.samples.size();
                record.fit = fit;
                run_results.levels.push_back(record);

                logger->info("Alignment factor " + std::to_string(factor) +
                             " completed: " + translation_string(fit.translation) +
                             ", weighted RMSE=" + std::to_string(fit.weighted_rmse_m) +
                             " m, samples=" + std::to_string(fit.samples_used));

                if (factor == 1) {
                    native_samples = std::move(sample_set);
                    native_gradient = std::move(gradient);
                }
            }
        }

        if (!base_fit.ok || !native_samples.has_value() || !native_gradient.has_value()) {
            throw std::runtime_error("Native-resolution alignment did not complete successfully.");
        }
        run_results.base_fit = base_fit;

        if (cfg.bootstrap.enabled) {
            ScopedTimer timer(timing, "spatial_block_bootstrap");
            run_results.bootstrap = run_block_bootstrap(native_samples->samples,
                                                        reference,
                                                        target,
                                                        *native_gradient,
                                                        base_fit,
                                                        cfg.alignment,
                                                        cfg.bootstrap,
                                                        cfg.runtime.debug_iterations,
                                                        *logger,
                                                        timing);
        }

        if (cfg.output.write_aligned_dsm) {
            ScopedTimer timer(timing, "write_aligned_dsm");
            write_aligned_dsm(cfg.output.directory / cfg.output.aligned_dsm,
                              target,
                              base_fit.translation.tx,
                              base_fit.translation.ty,
                              base_fit.translation.tz);
        }

        const auto before_results = Clock::now();
        const std::chrono::duration<double> total_elapsed = before_results - program_start;
        timing.add_seconds("total_before_results_json", total_elapsed.count());

        const std::filesystem::path results_path =
            cfg.output.directory / cfg.output.results_json;
        write_results_json(cfg, run_results, timing, results_path);

        logger->info("Alignment completed successfully: " + translation_string(base_fit.translation));
        logger->info("Results JSON: " + results_path.string());

        if (cfg.runtime.console_summary) {
            std::cout << std::fixed << std::setprecision(4)
                      << "DSM alignment completed successfully.\n"
                      << "tx = " << base_fit.translation.tx << " m\n"
                      << "ty = " << base_fit.translation.ty << " m\n"
                      << "tz = " << base_fit.translation.tz << " m\n"
                      << "horizontal shift = "
                      << std::hypot(base_fit.translation.tx, base_fit.translation.ty) << " m\n"
                      << "Results: " << results_path << '\n'
                      << "Log: " << (cfg.output.directory / cfg.output.log_file) << '\n';
        }

        return 0;
    } catch (const std::exception& e) {
        if (logger) logger->error(e.what());
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }
}

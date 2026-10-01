#include "dsm_coreg/results.hpp"

#include <cpl_json.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#ifndef DSM_COREG_VERSION
#define DSM_COREG_VERSION "unknown"
#endif

namespace dsm_coreg {
namespace {

CPLJSONObject translation_json(const Translation& t) {
    CPLJSONObject obj;
    obj.Add("tx_m", t.tx);
    obj.Add("ty_m", t.ty);
    obj.Add("tz_m", t.tz);
    obj.Add("horizontal_shift_m", std::hypot(t.tx, t.ty));
    return obj;
}

CPLJSONObject parameter_stats_json(const Stats3& stats, int index) {
    CPLJSONObject obj;
    obj.Add("mean_m", stats.mean[index]);
    obj.Add("std_m", stats.sd[index]);
    obj.Add("p05_m", stats.p05[index]);
    obj.Add("p95_m", stats.p95[index]);
    return obj;
}

}  // namespace

void write_results_json(const ToolConfig& cfg,
                        const RunResults& results,
                        const TimingRegistry& timing,
                        const std::filesystem::path& output_path) {
    CPLJSONObject root;
    root.Add("tool", "dsm_coreg");
    root.Add("version", DSM_COREG_VERSION);

    CPLJSONObject input;
    input.Add("reference_dsm", cfg.input.reference_dsm.string());
    input.Add("target_dsm", cfg.input.target_dsm.string());
    root.Add("input", input);

    CPLJSONObject alignment = translation_json(results.base_fit.translation);
    alignment.Add("weighted_rmse_m", results.base_fit.weighted_rmse_m);
    alignment.Add("samples_used", results.base_fit.samples_used);
    alignment.Add("iterations", results.base_fit.iterations);
    alignment.Add("converged", results.base_fit.converged);
    alignment.Add("ok", results.base_fit.ok);
    root.Add("alignment", alignment);

    CPLJSONArray levels;
    for (const AlignmentLevelRecord& record : results.levels) {
        CPLJSONObject item;
        item.Add("pyramid_factor", record.pyramid_factor);
        item.Add("reference_pixel_size_x_m", record.reference_pixel_size_x_m);
        item.Add("reference_pixel_size_y_m", record.reference_pixel_size_y_m);
        item.Add("candidate_samples", static_cast<GInt64>(record.candidate_samples));
        item.Add("translation", translation_json(record.fit.translation));
        item.Add("weighted_rmse_m", record.fit.weighted_rmse_m);
        item.Add("samples_used", record.fit.samples_used);
        item.Add("iterations", record.fit.iterations);
        item.Add("converged", record.fit.converged);
        item.Add("ok", record.fit.ok);
        levels.Add(item);
    }
    root.Add("alignment_levels", levels);

    if (results.coarse_search.has_value()) {
        const CoarseSearchResult& coarse = *results.coarse_search;
        CPLJSONObject obj;
        obj.Add("enabled", true);
        obj.Add("ok", coarse.ok);
        obj.Add("translation", translation_json(coarse.translation));
        obj.Add("robust_score_m", coarse.robust_score_m);
        obj.Add("samples_used", coarse.samples_used);
        obj.Add("candidates_evaluated", coarse.candidates_evaluated);
        root.Add("coarse_search", obj);
    } else {
        CPLJSONObject obj;
        obj.Add("enabled", false);
        root.Add("coarse_search", obj);
    }

    if (results.bootstrap.has_value()) {
        const BootstrapResult& bootstrap = *results.bootstrap;
        CPLJSONObject obj;
        obj.Add("enabled", true);
        obj.Add("method", "spatial_block_bootstrap");
        obj.Add("requested_replicates", bootstrap.requested_replicates);
        obj.Add("successful_replicates", bootstrap.successful_replicates);
        obj.Add("failed_replicates", bootstrap.failed_replicates);
        obj.Add("tx", parameter_stats_json(bootstrap.stats, 0));
        obj.Add("ty", parameter_stats_json(bootstrap.stats, 1));
        obj.Add("tz", parameter_stats_json(bootstrap.stats, 2));

        CPLJSONObject usage;
        for (const auto& [block_px, count] : bootstrap.block_size_usage) {
            usage.Add(std::to_string(block_px), static_cast<GInt64>(count));
        }
        obj.Add("block_size_usage_replicates", usage);
        obj.Add("bootstrap_geopackage_written", bootstrap.geopackage_written);
        if (bootstrap.geopackage_written) {
            obj.Add("bootstrap_geopackage", bootstrap.geopackage_path.string());
        }
        root.Add("uncertainty", obj);
    } else {
        CPLJSONObject obj;
        obj.Add("enabled", false);
        root.Add("uncertainty", obj);
    }

    CPLJSONObject timing_json;
    for (const auto& [name, seconds] : timing.snapshot()) {
        timing_json.Add(name, seconds);
    }
    root.Add("timing_seconds", timing_json);

    CPLJSONObject outputs;
    outputs.Add("results_json", output_path.string());
    outputs.Add("log_file", (cfg.output.directory / cfg.output.log_file).string());
    outputs.Add("effective_config_toml",
                (cfg.output.directory / cfg.output.effective_config_toml).string());
    if (cfg.output.write_aligned_dsm) {
        outputs.Add("aligned_dsm", (cfg.output.directory / cfg.output.aligned_dsm).string());
    }
    if (cfg.bootstrap.enabled && cfg.bootstrap.output.write_geopackage) {
        outputs.Add("bootstrap_geopackage", cfg.bootstrap.output.geopackage.string());
    }
    root.Add("outputs", outputs);

    std::ofstream stream(output_path, std::ios::out | std::ios::trunc);
    if (!stream) {
        throw std::runtime_error("Failed to open results JSON for writing: " + output_path.string());
    }
    stream << root.Format(CPLJSONObject::PrettyFormat::Pretty) << '\n';
    if (!stream) {
        throw std::runtime_error("Failed while writing results JSON: " + output_path.string());
    }
}

}  // namespace dsm_coreg

#pragma once

#include "dsm_coreg/alignment.hpp"
#include "dsm_coreg/bootstrap.hpp"
#include "dsm_coreg/config.hpp"
#include "dsm_coreg/timing.hpp"

#include <filesystem>
#include <optional>
#include <vector>

namespace dsm_coreg {

struct RunResults {
    FitResult base_fit;
    std::vector<AlignmentLevelRecord> levels;
    std::optional<CoarseSearchResult> coarse_search;
    std::optional<BootstrapResult> bootstrap;
};

void write_results_json(const ToolConfig& cfg,
                        const RunResults& results,
                        const TimingRegistry& timing,
                        const std::filesystem::path& output_path);

}  // namespace dsm_coreg

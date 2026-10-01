#include "dsm_coreg/config.hpp"

#include <toml++/toml.hpp>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace dsm_coreg {
namespace {

// Convert a TOML key to std::string once so error messages are easy to build.
std::string key_string(const toml::key& key) {
    return std::string(key.str());
}

// TOML is intentionally treated as a strict experiment specification.  A typo
// such as "huber_dleta_m" should not be ignored, because silently falling back
// to a default could invalidate a long scientific run.  Every section therefore
// checks its keys against an explicit allow-list.
void validate_allowed_keys(const toml::table& table,
                           const std::set<std::string>& allowed,
                           const std::string& context) {
    for (const auto& [key, node] : table) {
        (void)node;
        const std::string name = key_string(key);
        if (!allowed.contains(name)) {
            throw std::runtime_error(
                "Unknown configuration key '" + context + "." + name +
                "'. Unknown keys are rejected to prevent silent typos.");
        }
    }
}

// Require a TOML table such as [input].  The input section is mandatory; other
// sections may be absent and then use the defaults declared in config.hpp.
const toml::table& required_table(const toml::table& parent,
                                  std::string_view name) {
    const toml::node* node = parent.get(name);
    if (node == nullptr || !node->is_table()) {
        throw std::runtime_error(
            "Configuration section '[" + std::string(name) + "]' is required.");
    }
    return *node->as_table();
}

const toml::table* optional_table(const toml::table& parent,
                                  std::string_view name) {
    const toml::node* node = parent.get(name);
    if (node == nullptr) return nullptr;
    if (!node->is_table()) {
        throw std::runtime_error(
            "Configuration entry '" + std::string(name) + "' must be a TOML table.");
    }
    return node->as_table();
}

std::string required_string(const toml::table& table,
                            std::string_view key,
                            const std::string& context) {
    const toml::node* node = table.get(key);
    if (node == nullptr) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) + "' is required.");
    }
    const auto value = node->value<std::string>();
    if (!value.has_value()) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' must be a string.");
    }
    return *value;
}

std::string optional_string(const toml::table& table,
                            std::string_view key,
                            const std::string& default_value,
                            const std::string& context) {
    const toml::node* node = table.get(key);
    if (node == nullptr) return default_value;
    const auto value = node->value<std::string>();
    if (!value.has_value()) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' must be a string.");
    }
    return *value;
}

bool optional_bool(const toml::table& table,
                   std::string_view key,
                   bool default_value,
                   const std::string& context) {
    const toml::node* node = table.get(key);
    if (node == nullptr) return default_value;
    const auto value = node->value<bool>();
    if (!value.has_value()) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' must be true or false.");
    }
    return *value;
}

std::int64_t optional_int64(const toml::table& table,
                            std::string_view key,
                            std::int64_t default_value,
                            const std::string& context) {
    const toml::node* node = table.get(key);
    if (node == nullptr) return default_value;
    const auto value = node->value<std::int64_t>();
    if (!value.has_value()) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' must be an integer.");
    }
    return *value;
}

int optional_int(const toml::table& table,
                 std::string_view key,
                 int default_value,
                 const std::string& context) {
    const std::int64_t value = optional_int64(table, key, default_value, context);
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' is outside the range supported by int.");
    }
    return static_cast<int>(value);
}

std::uint64_t optional_uint64(const toml::table& table,
                              std::string_view key,
                              std::uint64_t default_value,
                              const std::string& context) {
    // TOML integer values are signed 64-bit.  Restrict the seed to the positive
    // signed range when it is supplied in the file, then convert explicitly.
    const toml::node* node = table.get(key);
    if (node == nullptr) return default_value;
    const auto value = node->value<std::int64_t>();
    if (!value.has_value() || *value < 0) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' must be a non-negative integer.");
    }
    return static_cast<std::uint64_t>(*value);
}

double optional_double(const toml::table& table,
                       std::string_view key,
                       double default_value,
                       const std::string& context) {
    const toml::node* node = table.get(key);
    if (node == nullptr) return default_value;

    // TOML distinguishes integer and floating-point nodes.  Scientific users
    // reasonably expect both "20" and "20.0" to work for a quantity in meters,
    // so accept either representation and convert integers to double.
    if (const auto value = node->value<double>(); value.has_value()) {
        return *value;
    }
    if (const auto value = node->value<std::int64_t>(); value.has_value()) {
        return static_cast<double>(*value);
    }

    throw std::runtime_error(
        "Configuration key '" + context + "." + std::string(key) +
        "' must be numeric.");
}

std::vector<int> optional_int_array(const toml::table& table,
                                    std::string_view key,
                                    const std::vector<int>& default_value,
                                    const std::string& context) {
    const toml::node* node = table.get(key);
    if (node == nullptr) return default_value;
    if (!node->is_array()) {
        throw std::runtime_error(
            "Configuration key '" + context + "." + std::string(key) +
            "' must be an integer array.");
    }

    const toml::array& array = *node->as_array();
    std::vector<int> values;
    values.reserve(array.size());

    for (std::size_t i = 0; i < array.size(); ++i) {
        const toml::node& item = array[i];
        const auto value = item.value<std::int64_t>();
        if (!value.has_value() ||
            *value < std::numeric_limits<int>::min() ||
            *value > std::numeric_limits<int>::max()) {
            throw std::runtime_error(
                "Every entry in '" + context + "." + std::string(key) +
                "' must be an integer representable by int.");
        }
        values.push_back(static_cast<int>(*value));
    }
    return values;
}

std::filesystem::path resolve_path(const std::filesystem::path& base,
                                   const std::filesystem::path& value) {
    if (value.empty()) return value;
    if (value.is_absolute()) return value.lexically_normal();
    return (base / value).lexically_normal();
}

void validate_semantics(const ToolConfig& cfg) {
    if (cfg.config_version != 1) {
        throw std::runtime_error(
            "Unsupported config_version=" + std::to_string(cfg.config_version) +
            ". This build supports config_version=1.");
    }

    if (cfg.input.reference_dsm.empty() || cfg.input.target_dsm.empty()) {
        throw std::runtime_error(
            "Both input.reference_dsm and input.target_dsm are required.");
    }
    if (!std::filesystem::exists(cfg.input.reference_dsm)) {
        throw std::runtime_error(
            "Reference DSM does not exist: " + cfg.input.reference_dsm.string());
    }
    if (!std::filesystem::exists(cfg.input.target_dsm)) {
        throw std::runtime_error(
            "Target DSM does not exist: " + cfg.input.target_dsm.string());
    }

    if (cfg.alignment.border_px < 0) {
        throw std::runtime_error("alignment.border_px must be >= 0.");
    }
    if (cfg.alignment.max_iterations < 1) {
        throw std::runtime_error("alignment.max_iterations must be >= 1.");
    }
    if (cfg.alignment.step_tolerance_m <= 0.0) {
        throw std::runtime_error("alignment.step_tolerance_m must be > 0.");
    }
    if (cfg.alignment.lambda < 0.0) {
        throw std::runtime_error("alignment.lambda must be >= 0.");
    }
    if (cfg.alignment.min_samples < 3) {
        throw std::runtime_error("alignment.min_samples must be >= 3.");
    }

    if (cfg.pyramid.enabled) {
        if (cfg.pyramid.factors.empty()) {
            throw std::runtime_error(
                "pyramid.factors cannot be empty when pyramid.enabled=true.");
        }
        for (const int factor : cfg.pyramid.factors) {
            if (factor < 1) {
                throw std::runtime_error("All pyramid factors must be >= 1.");
            }
        }
        for (std::size_t i = 1; i < cfg.pyramid.factors.size(); ++i) {
            if (cfg.pyramid.factors[i] >= cfg.pyramid.factors[i - 1]) {
                throw std::runtime_error(
                    "pyramid.factors must be strictly descending, e.g. [8, 4, 2, 1].");
            }
        }
        if (cfg.pyramid.factors.back() != 1) {
            throw std::runtime_error(
                "pyramid.factors must end with 1 so the final fit uses native resolution.");
        }
    }

    if (cfg.coarse_search.enabled) {
        if (cfg.coarse_search.range_x_m < 0.0 || cfg.coarse_search.range_y_m < 0.0) {
            throw std::runtime_error("coarse_search ranges must be >= 0.");
        }
        if (cfg.coarse_search.step_m <= 0.0) {
            throw std::runtime_error("coarse_search.step_m must be > 0.");
        }
        if (cfg.coarse_search.sample_stride < 1) {
            throw std::runtime_error("coarse_search.sample_stride must be >= 1.");
        }
        if (cfg.coarse_search.min_samples < 3) {
            throw std::runtime_error("coarse_search.min_samples must be >= 3.");
        }
    }

    if (cfg.bootstrap.enabled) {
        if (cfg.bootstrap.replicates < 1) {
            throw std::runtime_error("bootstrap.replicates must be >= 1.");
        }
        if (cfg.bootstrap.block_mode != "fixed" &&
            cfg.bootstrap.block_mode != "multiscale") {
            throw std::runtime_error(
                "bootstrap.block_mode must be 'fixed' or 'multiscale'.");
        }
        if (cfg.bootstrap.block_px < 4) {
            throw std::runtime_error("bootstrap.block_px must be >= 4.");
        }
        if (cfg.bootstrap.block_min_px < 4 || cfg.bootstrap.block_max_px < 4) {
            throw std::runtime_error(
                "bootstrap.block_min_px and block_max_px must be >= 4.");
        }
        if (cfg.bootstrap.block_max_px < cfg.bootstrap.block_min_px) {
            throw std::runtime_error(
                "bootstrap.block_max_px must be >= block_min_px.");
        }
        if (cfg.bootstrap.block_step_px < 1) {
            throw std::runtime_error("bootstrap.block_step_px must be >= 1.");
        }
        if (cfg.bootstrap.min_successful_replicates < 0) {
            throw std::runtime_error(
                "bootstrap.min_successful_replicates must be >= 0.");
        }
        if (cfg.bootstrap.output.write_geopackage &&
            cfg.bootstrap.output.geopackage.empty()) {
            throw std::runtime_error(
                "bootstrap.output.geopackage cannot be empty when "
                "write_geopackage=true.");
        }
    }

    if (cfg.runtime.threads < 0) {
        throw std::runtime_error("runtime.threads must be >= 0.");
    }

    if (cfg.output.directory.empty()) {
        throw std::runtime_error("output.directory cannot be empty.");
    }
    if (cfg.output.results_json.empty() ||
        cfg.output.log_file.empty() ||
        cfg.output.effective_config_toml.empty()) {
        throw std::runtime_error(
            "output.results_json, output.log_file, and output.effective_config_toml "
            "cannot be empty.");
    }
    if (cfg.output.write_aligned_dsm && cfg.output.aligned_dsm.empty()) {
        throw std::runtime_error(
            "output.aligned_dsm cannot be empty when write_aligned_dsm=true.");
    }
}

// Helper used only by write_effective_config_toml().  Keeping TOML construction
// here centralizes the output schema and ensures every default is recorded.
toml::array make_int_array(const std::vector<int>& values) {
    toml::array array;
    for (const int value : values) array.push_back(value);
    return array;
}

}  // namespace

ToolConfig load_config(const std::filesystem::path& config_path) {
    if (config_path.empty()) {
        throw std::runtime_error("Configuration path is empty.");
    }
    if (!std::filesystem::exists(config_path)) {
        throw std::runtime_error(
            "Configuration file does not exist: " + config_path.string());
    }

    toml::table root;
    try {
        root = toml::parse_file(config_path.string());
    } catch (const toml::parse_error& error) {
        std::ostringstream oss;
        oss << error;
        throw std::runtime_error(
            "Failed to parse TOML configuration '" + config_path.string() + "': " + oss.str());
    }

    validate_allowed_keys(
        root,
        {"config_version", "input", "alignment", "pyramid", "coarse_search",
         "bootstrap", "runtime", "output"},
        "root");

    ToolConfig cfg;
    cfg.config_directory = std::filesystem::absolute(config_path).parent_path();

    // ------------------------------- root -----------------------------------
    if (const toml::node* node = root.get("config_version"); node != nullptr) {
        const auto value = node->value<std::int64_t>();
        if (!value.has_value()) {
            throw std::runtime_error("config_version must be an integer.");
        }
        cfg.config_version = static_cast<int>(*value);
    }

    // ------------------------------ [input] ---------------------------------
    const toml::table& input = required_table(root, "input");
    validate_allowed_keys(input, {"reference_dsm", "target_dsm"}, "input");
    cfg.input.reference_dsm = resolve_path(
        cfg.config_directory,
        required_string(input, "reference_dsm", "input"));
    cfg.input.target_dsm = resolve_path(
        cfg.config_directory,
        required_string(input, "target_dsm", "input"));

    // ---------------------------- [alignment] -------------------------------
    if (const toml::table* alignment = optional_table(root, "alignment")) {
        validate_allowed_keys(
            *alignment,
            {"initial_tx_m", "initial_ty_m", "initial_tz_m", "border_px",
             "max_abs_dz_m", "huber_delta_m", "max_iterations",
             "step_tolerance_m", "lambda", "min_samples"},
            "alignment");

        cfg.alignment.initial_tx_m = optional_double(
            *alignment, "initial_tx_m", cfg.alignment.initial_tx_m, "alignment");
        cfg.alignment.initial_ty_m = optional_double(
            *alignment, "initial_ty_m", cfg.alignment.initial_ty_m, "alignment");
        cfg.alignment.initial_tz_m = optional_double(
            *alignment, "initial_tz_m", cfg.alignment.initial_tz_m, "alignment");
        cfg.alignment.border_px = optional_int(
            *alignment, "border_px", cfg.alignment.border_px, "alignment");
        cfg.alignment.max_abs_dz_m = optional_double(
            *alignment, "max_abs_dz_m", cfg.alignment.max_abs_dz_m, "alignment");
        cfg.alignment.huber_delta_m = optional_double(
            *alignment, "huber_delta_m", cfg.alignment.huber_delta_m, "alignment");
        cfg.alignment.max_iterations = optional_int(
            *alignment, "max_iterations", cfg.alignment.max_iterations, "alignment");
        cfg.alignment.step_tolerance_m = optional_double(
            *alignment, "step_tolerance_m", cfg.alignment.step_tolerance_m, "alignment");
        cfg.alignment.lambda = optional_double(
            *alignment, "lambda", cfg.alignment.lambda, "alignment");
        cfg.alignment.min_samples = optional_int(
            *alignment, "min_samples", cfg.alignment.min_samples, "alignment");
    }

    // ----------------------------- [pyramid] --------------------------------
    if (const toml::table* pyramid = optional_table(root, "pyramid")) {
        validate_allowed_keys(*pyramid, {"enabled", "factors"}, "pyramid");
        cfg.pyramid.enabled = optional_bool(
            *pyramid, "enabled", cfg.pyramid.enabled, "pyramid");
        cfg.pyramid.factors = optional_int_array(
            *pyramid, "factors", cfg.pyramid.factors, "pyramid");
    }

    // -------------------------- [coarse_search] -----------------------------
    if (const toml::table* coarse = optional_table(root, "coarse_search")) {
        validate_allowed_keys(
            *coarse,
            {"enabled", "range_x_m", "range_y_m", "step_m", "sample_stride",
             "min_samples"},
            "coarse_search");

        cfg.coarse_search.enabled = optional_bool(
            *coarse, "enabled", cfg.coarse_search.enabled, "coarse_search");
        cfg.coarse_search.range_x_m = optional_double(
            *coarse, "range_x_m", cfg.coarse_search.range_x_m, "coarse_search");
        cfg.coarse_search.range_y_m = optional_double(
            *coarse, "range_y_m", cfg.coarse_search.range_y_m, "coarse_search");
        cfg.coarse_search.step_m = optional_double(
            *coarse, "step_m", cfg.coarse_search.step_m, "coarse_search");
        cfg.coarse_search.sample_stride = optional_int(
            *coarse, "sample_stride", cfg.coarse_search.sample_stride, "coarse_search");
        cfg.coarse_search.min_samples = optional_int(
            *coarse, "min_samples", cfg.coarse_search.min_samples, "coarse_search");
    }

    // ---------------------------- [bootstrap] -------------------------------
    if (const toml::table* bootstrap = optional_table(root, "bootstrap")) {
        validate_allowed_keys(
            *bootstrap,
            {"enabled", "replicates", "block_mode", "block_px", "block_min_px",
             "block_max_px", "block_step_px", "seed",
             "min_successful_replicates", "output"},
            "bootstrap");

        cfg.bootstrap.enabled = optional_bool(
            *bootstrap, "enabled", cfg.bootstrap.enabled, "bootstrap");
        cfg.bootstrap.replicates = optional_int(
            *bootstrap, "replicates", cfg.bootstrap.replicates, "bootstrap");
        cfg.bootstrap.block_mode = optional_string(
            *bootstrap, "block_mode", cfg.bootstrap.block_mode, "bootstrap");
        cfg.bootstrap.block_px = optional_int(
            *bootstrap, "block_px", cfg.bootstrap.block_px, "bootstrap");
        cfg.bootstrap.block_min_px = optional_int(
            *bootstrap, "block_min_px", cfg.bootstrap.block_min_px, "bootstrap");
        cfg.bootstrap.block_max_px = optional_int(
            *bootstrap, "block_max_px", cfg.bootstrap.block_max_px, "bootstrap");
        cfg.bootstrap.block_step_px = optional_int(
            *bootstrap, "block_step_px", cfg.bootstrap.block_step_px, "bootstrap");
        cfg.bootstrap.seed = optional_uint64(
            *bootstrap, "seed", cfg.bootstrap.seed, "bootstrap");
        cfg.bootstrap.min_successful_replicates = optional_int(
            *bootstrap,
            "min_successful_replicates",
            cfg.bootstrap.min_successful_replicates,
            "bootstrap");

        // [bootstrap.output] contains optional audit/debug artifacts associated
        // specifically with the bootstrap.  Keeping these settings nested under
        // bootstrap avoids mixing statistical output with the primary aligned
        // DSM/result files in the top-level [output] section.
        if (const toml::table* bootstrap_output = optional_table(*bootstrap, "output")) {
            validate_allowed_keys(
                *bootstrap_output,
                {"write_geopackage", "geopackage"},
                "bootstrap.output");

            cfg.bootstrap.output.write_geopackage = optional_bool(
                *bootstrap_output,
                "write_geopackage",
                cfg.bootstrap.output.write_geopackage,
                "bootstrap.output");
            cfg.bootstrap.output.geopackage = optional_string(
                *bootstrap_output,
                "geopackage",
                cfg.bootstrap.output.geopackage.string(),
                "bootstrap.output");
        }
    }

    // ----------------------------- [runtime] --------------------------------
    if (const toml::table* runtime = optional_table(root, "runtime")) {
        validate_allowed_keys(
            *runtime,
            {"threads", "console_summary", "debug_iterations"},
            "runtime");
        cfg.runtime.threads = optional_int(
            *runtime, "threads", cfg.runtime.threads, "runtime");
        cfg.runtime.console_summary = optional_bool(
            *runtime, "console_summary", cfg.runtime.console_summary, "runtime");
        cfg.runtime.debug_iterations = optional_bool(
            *runtime, "debug_iterations", cfg.runtime.debug_iterations, "runtime");
    }

    // ------------------------------ [output] --------------------------------
    if (const toml::table* output = optional_table(root, "output")) {
        validate_allowed_keys(
            *output,
            {"directory", "results_json", "log_file", "effective_config_toml",
             "write_aligned_dsm", "aligned_dsm"},
            "output");

        cfg.output.directory = optional_string(
            *output, "directory", cfg.output.directory.string(), "output");
        cfg.output.results_json = optional_string(
            *output, "results_json", cfg.output.results_json, "output");
        cfg.output.log_file = optional_string(
            *output, "log_file", cfg.output.log_file, "output");
        cfg.output.effective_config_toml = optional_string(
            *output,
            "effective_config_toml",
            cfg.output.effective_config_toml,
            "output");
        cfg.output.write_aligned_dsm = optional_bool(
            *output, "write_aligned_dsm", cfg.output.write_aligned_dsm, "output");
        cfg.output.aligned_dsm = optional_string(
            *output, "aligned_dsm", cfg.output.aligned_dsm, "output");
    }

    // Resolve the output directory after parsing so all generated files remain
    // grouped with the experiment even when dsm_coreg is launched from a
    // different working directory.
    cfg.output.directory = resolve_path(cfg.config_directory, cfg.output.directory);

    // A relative bootstrap GeoPackage name is interpreted relative to the main
    // output directory rather than the TOML file.  This keeps all products from
    // one run together while still allowing an explicit absolute path when
    // desired.
    if (!cfg.bootstrap.output.geopackage.is_absolute()) {
        cfg.bootstrap.output.geopackage =
            (cfg.output.directory / cfg.bootstrap.output.geopackage).lexically_normal();
    } else {
        cfg.bootstrap.output.geopackage =
            cfg.bootstrap.output.geopackage.lexically_normal();
    }

    validate_semantics(cfg);
    return cfg;
}

void write_effective_config_toml(const ToolConfig& cfg,
                                 const std::filesystem::path& output_path) {
    // Build a fresh TOML document from the strongly typed configuration rather
    // than copying the user's input text.  This records every default value and
    // every resolved path actually used by the program.
    toml::table root;
    root.insert("config_version", cfg.config_version);

    toml::table input;
    input.insert("reference_dsm", cfg.input.reference_dsm.string());
    input.insert("target_dsm", cfg.input.target_dsm.string());
    root.insert("input", std::move(input));

    toml::table alignment;
    alignment.insert("initial_tx_m", cfg.alignment.initial_tx_m);
    alignment.insert("initial_ty_m", cfg.alignment.initial_ty_m);
    alignment.insert("initial_tz_m", cfg.alignment.initial_tz_m);
    alignment.insert("border_px", cfg.alignment.border_px);
    alignment.insert("max_abs_dz_m", cfg.alignment.max_abs_dz_m);
    alignment.insert("huber_delta_m", cfg.alignment.huber_delta_m);
    alignment.insert("max_iterations", cfg.alignment.max_iterations);
    alignment.insert("step_tolerance_m", cfg.alignment.step_tolerance_m);
    alignment.insert("lambda", cfg.alignment.lambda);
    alignment.insert("min_samples", cfg.alignment.min_samples);
    root.insert("alignment", std::move(alignment));

    toml::table pyramid;
    pyramid.insert("enabled", cfg.pyramid.enabled);
    pyramid.insert("factors", make_int_array(cfg.pyramid.factors));
    root.insert("pyramid", std::move(pyramid));

    toml::table coarse;
    coarse.insert("enabled", cfg.coarse_search.enabled);
    coarse.insert("range_x_m", cfg.coarse_search.range_x_m);
    coarse.insert("range_y_m", cfg.coarse_search.range_y_m);
    coarse.insert("step_m", cfg.coarse_search.step_m);
    coarse.insert("sample_stride", cfg.coarse_search.sample_stride);
    coarse.insert("min_samples", cfg.coarse_search.min_samples);
    root.insert("coarse_search", std::move(coarse));

    toml::table bootstrap;
    bootstrap.insert("enabled", cfg.bootstrap.enabled);
    bootstrap.insert("replicates", cfg.bootstrap.replicates);
    bootstrap.insert("block_mode", cfg.bootstrap.block_mode);
    bootstrap.insert("block_px", cfg.bootstrap.block_px);
    bootstrap.insert("block_min_px", cfg.bootstrap.block_min_px);
    bootstrap.insert("block_max_px", cfg.bootstrap.block_max_px);
    bootstrap.insert("block_step_px", cfg.bootstrap.block_step_px);
    bootstrap.insert("seed", static_cast<std::int64_t>(cfg.bootstrap.seed));
    bootstrap.insert("min_successful_replicates",
                     cfg.bootstrap.min_successful_replicates);

    toml::table bootstrap_output;
    bootstrap_output.insert("write_geopackage",
                            cfg.bootstrap.output.write_geopackage);
    bootstrap_output.insert("geopackage",
                            cfg.bootstrap.output.geopackage.string());
    bootstrap.insert("output", std::move(bootstrap_output));

    root.insert("bootstrap", std::move(bootstrap));

    toml::table runtime;
    runtime.insert("threads", cfg.runtime.threads);
    runtime.insert("console_summary", cfg.runtime.console_summary);
    runtime.insert("debug_iterations", cfg.runtime.debug_iterations);
    root.insert("runtime", std::move(runtime));

    toml::table output;
    output.insert("directory", cfg.output.directory.string());
    output.insert("results_json", cfg.output.results_json);
    output.insert("log_file", cfg.output.log_file);
    output.insert("effective_config_toml", cfg.output.effective_config_toml);
    output.insert("write_aligned_dsm", cfg.output.write_aligned_dsm);
    output.insert("aligned_dsm", cfg.output.aligned_dsm);
    root.insert("output", std::move(output));

    std::ofstream stream(output_path, std::ios::out | std::ios::trunc);
    if (!stream) {
        throw std::runtime_error(
            "Failed to open effective TOML configuration for writing: " +
            output_path.string());
    }
    stream << root << '\n';
    if (!stream) {
        throw std::runtime_error(
            "Failed while writing effective TOML configuration: " +
            output_path.string());
    }
}

}  // namespace dsm_coreg

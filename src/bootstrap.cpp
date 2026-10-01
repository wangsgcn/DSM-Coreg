#include "dsm_coreg/bootstrap.hpp"
#include "dsm_coreg/bootstrap_geopackage.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace dsm_coreg {
namespace {

// -----------------------------------------------------------------------------
// Spatial block partition
// -----------------------------------------------------------------------------
//
// A BlockPartition maps each populated spatial block to the indices of all
// samples that fall inside that block.  It also retains the block's integer row
// and column on the reference grid.  The row/column information is not needed by
// Gauss-Newton itself, but it is essential for reconstructing block polygons in
// the optional bootstrap GeoPackage audit output.
// -----------------------------------------------------------------------------
struct BlockPartition {
    int block_px = 0;
    int num_blocks = 0;

    // block_to_indices[b] lists the indices into the global `samples` vector for
    // every sample belonging to block b.
    std::vector<std::vector<int>> block_to_indices;

    // Spatial location of each compact block id.  These arrays have length
    // num_blocks and are indexed by the same block id as block_to_indices.
    std::vector<int> block_rows;
    std::vector<int> block_cols;
};

// Pack two signed 32-bit integers into a unique 64-bit hash key.  block_row and
// block_col are non-negative in normal use, but the explicit uint32 conversion
// keeps the packing operation well defined.
std::uint64_t block_key(int block_row, int block_col) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(block_row)) << 32U) |
           static_cast<std::uint32_t>(block_col);
}

// Partition all valid reference samples into non-overlapping block_px x block_px
// cells.  Only populated blocks receive a compact block id, which avoids storing
// large empty grids over NoData regions or irregular overlap footprints.
BlockPartition build_partition(const std::vector<Sample>& samples, int block_px) {
    BlockPartition partition;
    partition.block_px = block_px;

    std::unordered_map<std::uint64_t, int> id_by_cell;
    id_by_cell.reserve(4096);

    std::vector<int> sample_block(samples.size(), -1);
    int next_id = 0;

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const Sample& sample = samples[i];
        if (sample.ref_c < 0 || sample.ref_r < 0) continue;

        const int block_col = sample.ref_c / block_px;
        const int block_row = sample.ref_r / block_px;
        const std::uint64_t key = block_key(block_row, block_col);

        // Try to insert this spatial cell.  For a newly encountered cell, the
        // proposed value next_id becomes its permanent compact block id.
        auto [it, inserted] = id_by_cell.emplace(key, next_id);
        if (inserted) {
            partition.block_rows.push_back(block_row);
            partition.block_cols.push_back(block_col);
            ++next_id;
        }

        sample_block[i] = it->second;
    }

    partition.num_blocks = next_id;
    partition.block_to_indices.assign(
        static_cast<std::size_t>(partition.num_blocks), {});

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const int block_id = sample_block[i];
        if (block_id >= 0) {
            partition.block_to_indices[static_cast<std::size_t>(block_id)]
                .push_back(static_cast<int>(i));
        }
    }

    return partition;
}

// Linear interpolation between adjacent order statistics gives a stable
// percentile estimate without introducing another statistics dependency.
double percentile(std::vector<double> values, double q) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());

    const double position = q * static_cast<double>(values.size() - 1);
    const std::size_t i0 = static_cast<std::size_t>(std::floor(position));
    const std::size_t i1 = std::min(values.size() - 1, i0 + 1);
    const double fraction = position - static_cast<double>(i0);
    return (1.0 - fraction) * values[i0] + fraction * values[i1];
}

// Compute the marginal bootstrap summary statistics for tx, ty, and tz.  The
// standard deviation here is the sample standard deviation across successful
// GLOBAL replicate estimates.  It should be interpreted as the estimated
// sampling variability/spatial stability of the global registration estimator,
// not as absolute geolocation accuracy of either input DSM.
Stats3 compute_stats(const std::vector<std::array<double, 3>>& estimates) {
    Stats3 stats;
    if (estimates.empty()) return stats;

    const double n = static_cast<double>(estimates.size());
    for (const auto& e : estimates) {
        stats.mean[0] += e[0];
        stats.mean[1] += e[1];
        stats.mean[2] += e[2];
    }
    stats.mean[0] /= n;
    stats.mean[1] /= n;
    stats.mean[2] /= n;

    if (estimates.size() > 1) {
        double sumsq[3] = {0.0, 0.0, 0.0};
        for (const auto& e : estimates) {
            for (int k = 0; k < 3; ++k) {
                const double d = e[static_cast<std::size_t>(k)] - stats.mean[k];
                sumsq[k] += d * d;
            }
        }
        for (int k = 0; k < 3; ++k) {
            stats.sd[k] = std::sqrt(
                sumsq[k] / static_cast<double>(estimates.size() - 1));
        }
    }

    std::vector<double> tx;
    std::vector<double> ty;
    std::vector<double> tz;
    tx.reserve(estimates.size());
    ty.reserve(estimates.size());
    tz.reserve(estimates.size());
    for (const auto& e : estimates) {
        tx.push_back(e[0]);
        ty.push_back(e[1]);
        tz.push_back(e[2]);
    }

    stats.p05[0] = percentile(tx, 0.05);
    stats.p95[0] = percentile(tx, 0.95);
    stats.p05[1] = percentile(ty, 0.05);
    stats.p95[1] = percentile(ty, 0.95);
    stats.p05[2] = percentile(tz, 0.05);
    stats.p95[2] = percentile(tz, 0.95);

    return stats;
}

// SplitMix64-style mixing.  Each bootstrap replicate gets an RNG seed derived
// from its replicate index, not from the OpenMP thread that happened to execute
// it.  Therefore replicate #37 receives the same random block selection whether
// the program uses 1, 4, or 32 threads.
std::uint64_t mix_seed(std::uint64_t seed, std::uint64_t replicate_index) {
    std::uint64_t z = seed + 0x9e3779b97f4a7c15ULL * (replicate_index + 1ULL);
    z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31U);
}

// One slot is preallocated per bootstrap replicate.  OpenMP worker threads only
// modify their own slot, so no mutex is required during the expensive parallel
// fitting stage.
struct ReplicateSlot {
    int replicate_index = -1;
    int block_px = 0;
    FitResult fit;

    // Filled only when bootstrap.output.write_geopackage=true.  Avoiding this
    // vector during ordinary production runs prevents unnecessary memory use.
    std::vector<BootstrapBlockSelection> selected_blocks;
};

}  // namespace

BootstrapResult run_block_bootstrap(const std::vector<Sample>& samples,
                                    const DemRaster& reference,
                                    const DemRaster& target,
                                    const GradRaster& gradient,
                                    const FitResult& base_fit,
                                    const AlignmentConfig& alignment_cfg,
                                    const BootstrapConfig& bootstrap_cfg,
                                    bool debug_iterations,
                                    Logger& logger,
                                    TimingRegistry& timing) {
    BootstrapResult output;
    output.requested_replicates = bootstrap_cfg.replicates;

    if (!bootstrap_cfg.enabled) {
        output.ok = true;
        return output;
    }
    if (!base_fit.ok) {
        throw std::runtime_error(
            "Bootstrap cannot run because the base alignment fit is invalid.");
    }

    // -------------------------------------------------------------------------
    // Step 1: Build the spatial partitions once.
    // -------------------------------------------------------------------------
    // In fixed mode there is one partition.  In multiscale mode we build one
    // partition for each candidate block size, then each replicate randomly
    // chooses one of these partitions.  Reusing the partitions is important:
    // assigning millions of pixels to blocks inside every replicate would be
    // needlessly expensive.
    std::vector<int> block_sizes;
    std::vector<BlockPartition> partitions;
    {
        ScopedTimer timer(timing, "bootstrap_partition_build");

        if (bootstrap_cfg.block_mode == "fixed") {
            block_sizes.push_back(bootstrap_cfg.block_px);
        } else {
            for (int size = bootstrap_cfg.block_min_px;
                 size <= bootstrap_cfg.block_max_px;
                 size += bootstrap_cfg.block_step_px) {
                block_sizes.push_back(size);
            }
        }

        partitions.reserve(block_sizes.size());
        for (const int block_px : block_sizes) {
            BlockPartition partition = build_partition(samples, block_px);
            if (partition.num_blocks < 2) {
                throw std::runtime_error(
                    "Bootstrap block size " + std::to_string(block_px) +
                    " px produces fewer than two populated blocks.");
            }
            partitions.push_back(std::move(partition));
        }
    }

    {
        std::ostringstream oss;
        oss << "Starting spatial block bootstrap: replicates="
            << bootstrap_cfg.replicates
            << ", block_mode=" << bootstrap_cfg.block_mode
            << ", block_sizes_px=";
        for (std::size_t i = 0; i < block_sizes.size(); ++i) {
            if (i > 0) oss << ',';
            oss << block_sizes[i];
        }
        if (bootstrap_cfg.output.write_geopackage) {
            oss << ", GeoPackage audit output="
                << bootstrap_cfg.output.geopackage.string();
        }
        logger.info(oss.str());
    }

    // Pre-allocate one independent slot per requested replicate.  This also
    // preserves deterministic replicate ordering in the later GeoPackage even
    // though OpenMP executes replicates in a dynamic order.
    std::vector<ReplicateSlot> slots(
        static_cast<std::size_t>(bootstrap_cfg.replicates));

    // -------------------------------------------------------------------------
    // Step 2: Generate each bootstrap sample and solve one GLOBAL translation.
    // -------------------------------------------------------------------------
    {
        ScopedTimer timer(timing, "bootstrap_replicate_fits");

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
        for (int replicate = 0; replicate < bootstrap_cfg.replicates; ++replicate) {
            std::mt19937_64 rng(
                mix_seed(bootstrap_cfg.seed,
                         static_cast<std::uint64_t>(replicate)));

            // In multiscale mode one block size is chosen for the ENTIRE
            // replicate.  Therefore every selected block in a particular
            // bootstrap_XXXXXX layer has the same block_px value.
            std::uniform_int_distribution<int> partition_distribution(
                0, static_cast<int>(partitions.size()) - 1);
            const int partition_index =
                (partitions.size() == 1)
                    ? 0
                    : partition_distribution(rng);
            const BlockPartition& partition =
                partitions[static_cast<std::size_t>(partition_index)];

            const int K = partition.num_blocks;
            std::uniform_int_distribution<int> block_distribution(0, K - 1);
            std::vector<int> block_counts(static_cast<std::size_t>(K), 0);

            // Classical nonparametric block bootstrap:
            //
            //   - There are K available spatial blocks.
            //   - Draw K times WITH replacement.
            //   - block_counts[j] is therefore 0, 1, 2, 3, ... .
            //
            // gauss_newton_fit() uses these multiplicities as weights, which is
            // equivalent to physically copying a block that was drawn multiple
            // times but avoids copying millions of samples.
            for (int draw = 0; draw < K; ++draw) {
                const int selected = block_distribution(rng);
                ++block_counts[static_cast<std::size_t>(selected)];
            }

            // Bootstrap fits deliberately start from the converged full-data
            // solution.  The pyramid/coarse global initialization is NOT rerun
            // hundreds of times because the bootstrap is estimating local
            // uncertainty of the already identified global optimum.
            const bool log_this_replicate =
                debug_iterations && bootstrap_cfg.replicates <= 20;

            FitResult fit = gauss_newton_fit(
                samples,
                partition.block_to_indices,
                block_counts,
                target,
                gradient,
                base_fit.translation,
                alignment_cfg,
                log_this_replicate,
                log_this_replicate ? &logger : nullptr);

            // Treat a numerically non-finite translation as a failed fit even if
            // the optimizer returned ok=true.  This prevents NaNs from entering
            // the bootstrap statistics or the GeoPackage summary table.
            if (!(fit.ok &&
                  std::isfinite(fit.translation.tx) &&
                  std::isfinite(fit.translation.ty) &&
                  std::isfinite(fit.translation.tz))) {
                fit.ok = false;
            }

            ReplicateSlot& slot = slots[static_cast<std::size_t>(replicate)];
            slot.replicate_index = replicate;
            slot.block_px = partition.block_px;
            slot.fit = fit;

            // Optional audit artifact: record every UNIQUE selected block once,
            // together with multiplicity.  This is much smaller and clearer than
            // writing the same polygon multiple times when it was drawn twice or
            // three times.
            if (bootstrap_cfg.output.write_geopackage) {
                // Approximately 63.2% of K categories are expected to be unique
                // when drawing K times with replacement.  Reserving 2K/3 is a
                // useful approximation that reduces reallocations without
                // requiring an additional counting pass.
                slot.selected_blocks.reserve(
                    static_cast<std::size_t>((2 * K) / 3 + 1));

                for (int block_id = 0; block_id < K; ++block_id) {
                    const int multiplicity =
                        block_counts[static_cast<std::size_t>(block_id)];
                    if (multiplicity <= 0) continue;

                    BootstrapBlockSelection selected;
                    selected.block_id = block_id;
                    selected.block_row =
                        partition.block_rows[static_cast<std::size_t>(block_id)];
                    selected.block_col =
                        partition.block_cols[static_cast<std::size_t>(block_id)];
                    selected.block_px = partition.block_px;
                    selected.sample_count = static_cast<int>(
                        partition.block_to_indices[
                            static_cast<std::size_t>(block_id)].size());
                    selected.multiplicity = multiplicity;
                    slot.selected_blocks.push_back(selected);
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Step 3: Collect successful GLOBAL replicate estimates.
    // -------------------------------------------------------------------------
    std::vector<std::array<double, 3>> successful;
    successful.reserve(slots.size());

    for (const ReplicateSlot& slot : slots) {
        ++output.block_size_usage[slot.block_px];
        if (slot.fit.ok) {
            successful.push_back({slot.fit.translation.tx,
                                  slot.fit.translation.ty,
                                  slot.fit.translation.tz});
        }
    }

    output.successful_replicates = static_cast<int>(successful.size());
    output.failed_replicates =
        output.requested_replicates - output.successful_replicates;

    // -------------------------------------------------------------------------
    // Step 4: Optional GeoPackage audit output.
    // -------------------------------------------------------------------------
    // Write this BEFORE enforcing the minimum-success threshold.  If a bootstrap
    // run fails badly, the block-selection layers are especially valuable for
    // diagnosing the failure, so preserving them is preferable to throwing away
    // the intermediate artifact.
    if (bootstrap_cfg.output.write_geopackage) {
        ScopedTimer timer(timing, "bootstrap_geopackage_write");

        std::vector<BootstrapReplicateArtifact> artifacts;
        artifacts.reserve(slots.size());

        for (ReplicateSlot& slot : slots) {
            BootstrapReplicateArtifact artifact;
            artifact.replicate_index = slot.replicate_index;
            artifact.block_px = slot.block_px;
            artifact.fit = slot.fit;
            artifact.selected_blocks = std::move(slot.selected_blocks);
            artifacts.push_back(std::move(artifact));
        }

        write_bootstrap_geopackage(
            bootstrap_cfg.output.geopackage,
            reference,
            artifacts);

        output.geopackage_written = true;
        output.geopackage_path = bootstrap_cfg.output.geopackage;
        logger.info(
            "Bootstrap GeoPackage written: " +
            output.geopackage_path.string() +
            " (one polygon layer per requested replicate plus "
            "bootstrap_replicates summary table).");
    }

    // -------------------------------------------------------------------------
    // Step 5: Validate replicate count and compute uncertainty statistics.
    // -------------------------------------------------------------------------
    const int minimum_success =
        (bootstrap_cfg.min_successful_replicates > 0)
            ? bootstrap_cfg.min_successful_replicates
            : std::max(20, bootstrap_cfg.replicates / 5);

    if (output.successful_replicates < minimum_success) {
        throw std::runtime_error(
            "Too few successful bootstrap fits: " +
            std::to_string(output.successful_replicates) +
            " of " + std::to_string(output.requested_replicates) +
            (output.geopackage_written
                 ? ". The bootstrap GeoPackage was still written for diagnosis."
                 : "."));
    }

    output.stats = compute_stats(successful);
    output.ok = true;

    {
        std::ostringstream oss;
        oss << "Bootstrap complete: successful="
            << output.successful_replicates
            << ", failed=" << output.failed_replicates
            << ", sd(tx)=" << output.stats.sd[0]
            << " m, sd(ty)=" << output.stats.sd[1]
            << " m, sd(tz)=" << output.stats.sd[2] << " m";
        logger.info(oss.str());
    }

    return output;
}

}  // namespace dsm_coreg

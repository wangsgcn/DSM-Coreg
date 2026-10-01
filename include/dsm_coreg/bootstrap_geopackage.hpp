#pragma once

#include "dsm_coreg/alignment.hpp"
#include "dsm_coreg/raster.hpp"

#include <filesystem>
#include <vector>

namespace dsm_coreg {

// -----------------------------------------------------------------------------
// Bootstrap GeoPackage audit structures
// -----------------------------------------------------------------------------
//
// The statistical bootstrap produces ONE global translation estimate for each
// replicate.  The individual blocks selected within a replicate do *not* have
// separate tx/ty/tz values.  These structures preserve that distinction:
//
//   BootstrapReplicateArtifact::fit
//       = the global tx/ty/tz solution for the complete bootstrap replicate.
//
//   BootstrapBlockSelection
//       = one spatial block that participated in that replicate, together with
//         the number of times it was drawn (its bootstrap multiplicity).
//
// Keeping the two concepts separate avoids a common interpretation error where
// a replicate-level parameter estimate is accidentally treated as a local
// per-block estimate.
// -----------------------------------------------------------------------------

struct BootstrapBlockSelection {
    int block_id = -1;
    int block_row = -1;
    int block_col = -1;
    int block_px = 0;

    // Number of reference samples/pixels that belong to this populated block.
    // This is useful when a block intersects NoData or the filtered sample mask.
    int sample_count = 0;

    // Number of times this block was selected in the replicate.  The bootstrap
    // draws K blocks with replacement from K available blocks, so multiplicity
    // can be 1, 2, 3, ... .  Blocks with multiplicity 0 are intentionally not
    // stored in the per-replicate polygon layer.
    int multiplicity = 0;
};

struct BootstrapReplicateArtifact {
    // Zero-based replicate index used internally by the computation.  The
    // GeoPackage layer name and summary table expose replicate_id as one-based
    // (1, 2, ..., B), which is more natural when browsing the file manually.
    int replicate_index = -1;

    // Block size used by this replicate.  In fixed mode this is the same for
    // every replicate; in multiscale mode it can differ between replicates.
    int block_px = 0;

    // The GLOBAL alignment fit obtained from all selected blocks in this
    // replicate.  tx/ty/tz therefore belong to the replicate, not to any one
    // polygon below.
    FitResult fit;

    // Unique selected blocks.  A block selected multiple times appears only
    // once here; multiplicity records the number of draws.  This representation
    // is much smaller and clearer than duplicating identical polygons.
    std::vector<BootstrapBlockSelection> selected_blocks;
};

// Write a GeoPackage designed for visual inspection and reproducibility of the
// spatial block bootstrap.
//
// The output contains:
//
//   1. Exactly one POLYGON feature layer per requested bootstrap replicate:
//        bootstrap_000001
//        bootstrap_000002
//        ...
//
//      Each polygon represents a selected spatial block.  Polygon attributes
//      describe block identity, location, sample count, and multiplicity.  The
//      replicate-level tx/ty/tz are NOT repeated on every polygon.
//
//   2. One non-spatial table named "bootstrap_replicates" with one row per
//      replicate.  This table stores the GLOBAL tx/ty/tz, fit diagnostics, and
//      block-selection counts.  It is the authoritative tabular location for
//      replicate-level translation values.
//
//   3. Replicate-level layer metadata (tx, ty, tz, success, etc.) is also added
//      to each polygon layer as a convenience for tools that expose OGR layer
//      metadata.  The summary table remains the portable/reliable representation.
//
// The reference DSM supplies the geotransform, dimensions, and CRS used to turn
// block row/column indices into map-coordinate polygons.
void write_bootstrap_geopackage(
    const std::filesystem::path& output_path,
    const DemRaster& reference,
    const std::vector<BootstrapReplicateArtifact>& replicates);

}  // namespace dsm_coreg

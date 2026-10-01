#include "dsm_coreg/bootstrap_geopackage.hpp"

#include <gdal_priv.h>
#include <ogrsf_frmts.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace dsm_coreg {
namespace {

// GDAL datasets are C-style objects that must be closed with GDALClose().
// Wrapping the raw pointer in a unique_ptr makes cleanup automatic even when an
// exception is thrown while creating one of hundreds of replicate layers.
struct GDALDatasetCloser {
    void operator()(GDALDataset* dataset) const {
        if (dataset != nullptr) GDALClose(dataset);
    }
};
using DatasetPtr = std::unique_ptr<GDALDataset, GDALDatasetCloser>;

// Create an OGR field and fail loudly if the GeoPackage driver rejects it.
// Scientific audit output should never silently omit a requested attribute.
void add_field(OGRLayer* layer,
               const char* name,
               OGRFieldType type,
               int width = 0,
               int precision = 0) {
    OGRFieldDefn field(name, type);
    if (width > 0) field.SetWidth(width);
    if (precision > 0) field.SetPrecision(precision);

    if (layer->CreateField(&field) != OGRERR_NONE) {
        throw std::runtime_error(
            "Failed to create GeoPackage field '" + std::string(name) +
            "' in layer '" + layer->GetName() + "'.");
    }
}

// Convert an integer replicate index to a stable, lexically sortable layer
// name.  Six digits comfortably supports hundreds of thousands of replicates
// while keeping ordinary 200- or 1000-replicate files easy to browse.
std::string replicate_layer_name(int replicate_index) {
    std::ostringstream oss;
    oss << "bootstrap_" << std::setw(6) << std::setfill('0')
        << (replicate_index + 1);
    return oss.str();
}

// Convert a value to a high-precision string for layer metadata.  The summary
// table stores the same values as native floating-point fields; metadata is only
// a convenience copy and therefore uses text, as required by GDAL metadata APIs.
std::string double_string(double value) {
    std::ostringstream oss;
    oss << std::setprecision(17) << value;
    return oss.str();
}

// Compute one block polygon from its integer position in the reference raster.
//
// A bootstrap block is defined on the reference pixel grid.  block_row and
// block_col identify the upper-left block cell in units of block_px pixels.  The
// block footprint is clipped at the raster edges, then its four pixel-corner
// coordinates are transformed by the GDAL geotransform.
//
// Although dsm_coreg currently requires north-up rasters, the complete affine
// GDAL transform is used here (including gt[2] and gt[4]) so this helper remains
// mathematically correct if that restriction is relaxed later.
OGRPolygon make_block_polygon(const DemRaster& reference,
                              int block_row,
                              int block_col,
                              int block_px) {
    const int c0 = std::clamp(block_col * block_px, 0, reference.width);
    const int r0 = std::clamp(block_row * block_px, 0, reference.height);
    const int c1 = std::clamp(c0 + block_px, 0, reference.width);
    const int r1 = std::clamp(r0 + block_px, 0, reference.height);

    auto corner_to_world = [&](int c, int r, double& x, double& y) {
        x = reference.gt[0] + static_cast<double>(c) * reference.gt[1] +
            static_cast<double>(r) * reference.gt[2];
        y = reference.gt[3] + static_cast<double>(c) * reference.gt[4] +
            static_cast<double>(r) * reference.gt[5];
    };

    double x00 = 0.0, y00 = 0.0;
    double x10 = 0.0, y10 = 0.0;
    double x11 = 0.0, y11 = 0.0;
    double x01 = 0.0, y01 = 0.0;
    corner_to_world(c0, r0, x00, y00);
    corner_to_world(c1, r0, x10, y10);
    corner_to_world(c1, r1, x11, y11);
    corner_to_world(c0, r1, x01, y01);

    OGRLinearRing ring;
    ring.addPoint(x00, y00);
    ring.addPoint(x10, y10);
    ring.addPoint(x11, y11);
    ring.addPoint(x01, y01);
    ring.addPoint(x00, y00);  // Explicitly close the ring.

    OGRPolygon polygon;
    polygon.addRing(&ring);
    return polygon;
}

// Construct the GeoPackage's non-spatial replicate summary table.  This table
// stores each replicate's global tx/ty/tz exactly once, which is statistically
// cleaner than duplicating the same values across every block polygon.
OGRLayer* create_summary_table(GDALDataset* dataset) {
    OGRLayer* layer = dataset->CreateLayer(
        "bootstrap_replicates", nullptr, wkbNone, nullptr);
    if (layer == nullptr) {
        throw std::runtime_error(
            "Failed to create non-spatial table 'bootstrap_replicates'.");
    }

    add_field(layer, "replicate_id", OFTInteger);
    add_field(layer, "block_px", OFTInteger);
    add_field(layer, "success", OFTInteger);
    add_field(layer, "converged", OFTInteger);
    add_field(layer, "tx_m", OFTReal, 0, 10);
    add_field(layer, "ty_m", OFTReal, 0, 10);
    add_field(layer, "tz_m", OFTReal, 0, 10);
    add_field(layer, "shift_xy_m", OFTReal, 0, 10);
    add_field(layer, "wrmse_m", OFTReal, 0, 10);
    add_field(layer, "samples_used", OFTInteger);
    add_field(layer, "iterations", OFTInteger);
    add_field(layer, "unique_blocks", OFTInteger);
    add_field(layer, "total_draws", OFTInteger);

    return layer;
}

// Append exactly one summary record for a replicate.  Failed replicates are
// intentionally retained: success=0 and the floating-point fit fields remain
// NULL.  Keeping failures makes the GeoPackage a complete audit trail rather
// than a filtered view containing only successful bootstrap samples.
void write_summary_feature(OGRLayer* summary,
                           const BootstrapReplicateArtifact& replicate) {
    std::unique_ptr<OGRFeature, decltype(&OGRFeature::DestroyFeature)> feature(
        OGRFeature::CreateFeature(summary->GetLayerDefn()),
        &OGRFeature::DestroyFeature);

    if (!feature) {
        throw std::runtime_error("Failed to allocate bootstrap summary feature.");
    }

    const int replicate_id = replicate.replicate_index + 1;
    feature->SetField("replicate_id", replicate_id);
    feature->SetField("block_px", replicate.block_px);
    feature->SetField("success", replicate.fit.ok ? 1 : 0);
    feature->SetField("converged", replicate.fit.converged ? 1 : 0);
    feature->SetField("samples_used", replicate.fit.samples_used);
    feature->SetField("iterations", replicate.fit.iterations);
    feature->SetField("unique_blocks",
                      static_cast<int>(replicate.selected_blocks.size()));

    int total_draws = 0;
    for (const BootstrapBlockSelection& block : replicate.selected_blocks) {
        total_draws += block.multiplicity;
    }
    feature->SetField("total_draws", total_draws);

    // A failed fit has no meaningful global translation.  Leaving these fields
    // NULL is preferable to writing zero, because zero could be mistaken for a
    // valid estimate of no translation.
    if (replicate.fit.ok &&
        std::isfinite(replicate.fit.translation.tx) &&
        std::isfinite(replicate.fit.translation.ty) &&
        std::isfinite(replicate.fit.translation.tz)) {
        feature->SetField("tx_m", replicate.fit.translation.tx);
        feature->SetField("ty_m", replicate.fit.translation.ty);
        feature->SetField("tz_m", replicate.fit.translation.tz);
        feature->SetField(
            "shift_xy_m",
            std::hypot(replicate.fit.translation.tx, replicate.fit.translation.ty));
        if (std::isfinite(replicate.fit.weighted_rmse_m)) {
            feature->SetField("wrmse_m", replicate.fit.weighted_rmse_m);
        }
    }

    if (summary->CreateFeature(feature.get()) != OGRERR_NONE) {
        throw std::runtime_error(
            "Failed to write a record to 'bootstrap_replicates'.");
    }
}

// Create the polygon layer for one replicate and write all uniquely selected
// blocks.  Replicate-level translations are stored as layer metadata and in the
// summary table, never as per-polygon attributes.
void write_replicate_layer(GDALDataset* dataset,
                           OGRSpatialReference* spatial_reference,
                           const DemRaster& reference,
                           const BootstrapReplicateArtifact& replicate) {
    const std::string layer_name = replicate_layer_name(replicate.replicate_index);

    OGRLayer* layer = dataset->CreateLayer(
        layer_name.c_str(), spatial_reference, wkbPolygon, nullptr);
    if (layer == nullptr) {
        throw std::runtime_error(
            "Failed to create bootstrap replicate layer '" + layer_name + "'.");
    }

    // Polygon-level fields describe only the block selection itself.
    add_field(layer, "block_id", OFTInteger);
    add_field(layer, "block_row", OFTInteger);
    add_field(layer, "block_col", OFTInteger);
    add_field(layer, "block_px", OFTInteger);
    add_field(layer, "multiplicity", OFTInteger);
    add_field(layer, "n_samples", OFTInteger);

    // Add global replicate information as layer metadata.  Many GDAL-aware
    // applications can inspect this with ogrinfo.  QGIS support for arbitrary
    // layer metadata varies, which is why bootstrap_replicates also stores the
    // same information in a portable tabular form.
    const std::string replicate_id = std::to_string(replicate.replicate_index + 1);
    const std::string block_px = std::to_string(replicate.block_px);
    const std::string success = replicate.fit.ok ? "1" : "0";
    const std::string converged = replicate.fit.converged ? "1" : "0";
    layer->SetMetadataItem("replicate_id", replicate_id.c_str());
    layer->SetMetadataItem("block_px", block_px.c_str());
    layer->SetMetadataItem("fit_success", success.c_str());
    layer->SetMetadataItem("fit_converged", converged.c_str());

    if (replicate.fit.ok &&
        std::isfinite(replicate.fit.translation.tx) &&
        std::isfinite(replicate.fit.translation.ty) &&
        std::isfinite(replicate.fit.translation.tz)) {
        const std::string tx = double_string(replicate.fit.translation.tx);
        const std::string ty = double_string(replicate.fit.translation.ty);
        const std::string tz = double_string(replicate.fit.translation.tz);
        const std::string shift = double_string(
            std::hypot(replicate.fit.translation.tx, replicate.fit.translation.ty));
        const std::string rmse = double_string(replicate.fit.weighted_rmse_m);
        layer->SetMetadataItem("tx_m", tx.c_str());
        layer->SetMetadataItem("ty_m", ty.c_str());
        layer->SetMetadataItem("tz_m", tz.c_str());
        layer->SetMetadataItem("horizontal_shift_m", shift.c_str());
        layer->SetMetadataItem("weighted_rmse_m", rmse.c_str());
    }

    // A transaction is important for performance.  Without it, SQLite (the
    // storage engine underneath GeoPackage) may sync to disk for every polygon.
    // With hundreds of replicate layers that can become dramatically slower.
    if (layer->StartTransaction() != OGRERR_NONE) {
        throw std::runtime_error(
            "Failed to start transaction for layer '" + layer_name + "'.");
    }

    bool committed = false;
    try {
        for (const BootstrapBlockSelection& block : replicate.selected_blocks) {
            std::unique_ptr<OGRFeature, decltype(&OGRFeature::DestroyFeature)> feature(
                OGRFeature::CreateFeature(layer->GetLayerDefn()),
                &OGRFeature::DestroyFeature);
            if (!feature) {
                throw std::runtime_error(
                    "Failed to allocate polygon feature for layer '" + layer_name + "'.");
            }

            feature->SetField("block_id", block.block_id);
            feature->SetField("block_row", block.block_row);
            feature->SetField("block_col", block.block_col);
            feature->SetField("block_px", block.block_px);
            feature->SetField("multiplicity", block.multiplicity);
            feature->SetField("n_samples", block.sample_count);

            OGRPolygon polygon = make_block_polygon(
                reference, block.block_row, block.block_col, block.block_px);
            if (feature->SetGeometry(&polygon) != OGRERR_NONE) {
                throw std::runtime_error(
                    "Failed to attach block geometry in layer '" + layer_name + "'.");
            }

            if (layer->CreateFeature(feature.get()) != OGRERR_NONE) {
                throw std::runtime_error(
                    "Failed to write block polygon in layer '" + layer_name + "'.");
            }
        }

        if (layer->CommitTransaction() != OGRERR_NONE) {
            throw std::runtime_error(
                "Failed to commit transaction for layer '" + layer_name + "'.");
        }
        committed = true;
    } catch (...) {
        if (!committed) layer->RollbackTransaction();
        throw;
    }
}

}  // namespace

void write_bootstrap_geopackage(
    const std::filesystem::path& output_path,
    const DemRaster& reference,
    const std::vector<BootstrapReplicateArtifact>& replicates) {
    if (replicates.empty()) {
        throw std::runtime_error(
            "Cannot write bootstrap GeoPackage because no replicate artifacts were supplied.");
    }

    reference.validate_north_up("reference");

    // Ensure the parent directory exists even when the user chose an absolute
    // GeoPackage path outside output.directory.
    if (!output_path.parent_path().empty()) {
        std::filesystem::create_directories(output_path.parent_path());
    }

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GPKG");
    if (driver == nullptr) {
        throw std::runtime_error(
            "GDAL GeoPackage (GPKG) driver is unavailable in this build.");
    }

    // Reproducibility is easier when each run creates a clean file.  Otherwise
    // stale replicate layers from a previous run with more replicates could
    // remain in the GeoPackage and be mistaken for current output.
    if (std::filesystem::exists(output_path)) {
        if (driver->Delete(output_path.string().c_str()) != CE_None) {
            throw std::runtime_error(
                "Failed to remove existing bootstrap GeoPackage: " +
                output_path.string());
        }
    }

    DatasetPtr dataset(driver->Create(
        output_path.string().c_str(), 0, 0, 0, GDT_Unknown, nullptr));
    if (!dataset) {
        throw std::runtime_error(
            "Failed to create bootstrap GeoPackage: " + output_path.string());
    }

    // Build the CRS from the reference DSM WKT.  All block polygons use the
    // reference raster grid, so the reference CRS is the correct layer CRS.
    OGRSpatialReference spatial_reference;
    if (!reference.projection_wkt.empty()) {
        if (spatial_reference.SetFromUserInput(reference.projection_wkt.c_str()) != OGRERR_NONE) {
            throw std::runtime_error(
                "Failed to construct GeoPackage spatial reference from reference DSM CRS.");
        }
        spatial_reference.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    }
    OGRSpatialReference* srs_ptr = reference.projection_wkt.empty()
                                       ? nullptr
                                       : &spatial_reference;

    OGRLayer* summary = create_summary_table(dataset.get());

    // Write summary rows in one transaction.  This table is small, but a single
    // transaction is still cleaner and faster than one SQLite transaction per
    // replicate.
    if (summary->StartTransaction() != OGRERR_NONE) {
        throw std::runtime_error(
            "Failed to start transaction for bootstrap_replicates table.");
    }
    bool summary_committed = false;
    try {
        for (const BootstrapReplicateArtifact& replicate : replicates) {
            write_summary_feature(summary, replicate);
        }
        if (summary->CommitTransaction() != OGRERR_NONE) {
            throw std::runtime_error(
                "Failed to commit bootstrap_replicates summary table.");
        }
        summary_committed = true;
    } catch (...) {
        if (!summary_committed) summary->RollbackTransaction();
        throw;
    }

    // Exactly one spatial feature layer is created for every requested
    // replicate, including failed fits.  A failed replicate still has a valid
    // block selection that is useful for diagnosing why the fit failed.
    for (const BootstrapReplicateArtifact& replicate : replicates) {
        write_replicate_layer(dataset.get(), srs_ptr, reference, replicate);
    }

    // GDALClose(), invoked automatically by DatasetPtr, flushes GeoPackage
    // metadata and SQLite pages to disk.
}

}  // namespace dsm_coreg

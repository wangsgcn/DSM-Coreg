#include "dsm_coreg/raster.hpp"

#include <cpl_conv.h>
#include <gdal_priv.h>
#include <ogr_spatialref.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace dsm_coreg {
namespace {

inline bool finite(double value) {
    return std::isfinite(value);
}

void copy_geotransform(const double src[6], double dst[6]) {
    std::copy(src, src + 6, dst);
}

}  // namespace

void DemRaster::validate_north_up(const std::string& name) const {
    // The alignment implementation converts between map coordinates and pixel
    // coordinates using the simplified north-up affine form.  Rotated rasters
    // should first be warped to a north-up projected grid with GDAL.
    if (std::abs(gt[2]) > 1.0e-12 || std::abs(gt[4]) > 1.0e-12) {
        throw std::runtime_error(name +
                                 ": rotated geotransform is not supported (gt[2] or gt[4] != 0). "
                                 "Warp the raster to a north-up grid before alignment.");
    }
    if (gt[1] == 0.0 || gt[5] == 0.0) {
        throw std::runtime_error(name + ": invalid zero raster pixel size.");
    }
}

bool DemRaster::in_bounds(int c, int r) const {
    return c >= 0 && c < width && r >= 0 && r < height;
}

double DemRaster::at(int c, int r) const {
    return z[static_cast<std::size_t>(r) * static_cast<std::size_t>(width) +
             static_cast<std::size_t>(c)];
}

bool DemRaster::valid_z(double value) const {
    if (!finite(value)) return false;
    if (has_no_data && finite(no_data) && std::abs(value - no_data) <= 1.0e-12) return false;
    return true;
}

void DemRaster::pixel_center_to_world(int c, int r, double& x, double& y) const {
    x = gt[0] + (static_cast<double>(c) + 0.5) * gt[1] +
        (static_cast<double>(r) + 0.5) * gt[2];
    y = gt[3] + (static_cast<double>(c) + 0.5) * gt[4] +
        (static_cast<double>(r) + 0.5) * gt[5];
}

void DemRaster::world_to_pixel_center(double x, double y, double& c, double& r) const {
    c = (x - gt[0]) / gt[1] - 0.5;
    r = (y - gt[3]) / gt[5] - 0.5;  // gt[5] is normally negative.
}

bool DemRaster::bilinear_sample(double x, double y, double& value) const {
    double cf = 0.0;
    double rf = 0.0;
    world_to_pixel_center(x, y, cf, rf);

    const int c0 = static_cast<int>(std::floor(cf));
    const int r0 = static_cast<int>(std::floor(rf));
    const int c1 = c0 + 1;
    const int r1 = r0 + 1;

    if (!in_bounds(c0, r0) || !in_bounds(c1, r1)) return false;

    const double z00 = at(c0, r0);
    const double z10 = at(c1, r0);
    const double z01 = at(c0, r1);
    const double z11 = at(c1, r1);
    if (!valid_z(z00) || !valid_z(z10) || !valid_z(z01) || !valid_z(z11)) return false;

    const double dx = cf - static_cast<double>(c0);
    const double dy = rf - static_cast<double>(r0);

    value = (1.0 - dx) * (1.0 - dy) * z00 +
            dx * (1.0 - dy) * z10 +
            (1.0 - dx) * dy * z01 +
            dx * dy * z11;
    return finite(value);
}

bool GradRaster::in_bounds(int c, int r) const {
    return c >= 0 && c < width && r >= 0 && r < height;
}

bool GradRaster::valid(double value) const {
    return finite(value);
}

double GradRaster::at_gx(int c, int r) const {
    return gx[static_cast<std::size_t>(r) * static_cast<std::size_t>(width) +
              static_cast<std::size_t>(c)];
}

double GradRaster::at_gy(int c, int r) const {
    return gy[static_cast<std::size_t>(r) * static_cast<std::size_t>(width) +
              static_cast<std::size_t>(c)];
}

void GradRaster::world_to_pixel_center(double x, double y, double& c, double& r) const {
    c = (x - gt[0]) / gt[1] - 0.5;
    r = (y - gt[3]) / gt[5] - 0.5;
}

bool GradRaster::bilinear_sample(double x, double y, double& out_gx, double& out_gy) const {
    double cf = 0.0;
    double rf = 0.0;
    world_to_pixel_center(x, y, cf, rf);

    const int c0 = static_cast<int>(std::floor(cf));
    const int r0 = static_cast<int>(std::floor(rf));
    const int c1 = c0 + 1;
    const int r1 = r0 + 1;
    if (!in_bounds(c0, r0) || !in_bounds(c1, r1)) return false;

    const double gx00 = at_gx(c0, r0);
    const double gx10 = at_gx(c1, r0);
    const double gx01 = at_gx(c0, r1);
    const double gx11 = at_gx(c1, r1);
    const double gy00 = at_gy(c0, r0);
    const double gy10 = at_gy(c1, r0);
    const double gy01 = at_gy(c0, r1);
    const double gy11 = at_gy(c1, r1);

    if (!valid(gx00) || !valid(gx10) || !valid(gx01) || !valid(gx11) ||
        !valid(gy00) || !valid(gy10) || !valid(gy01) || !valid(gy11)) {
        return false;
    }

    const double dx = cf - static_cast<double>(c0);
    const double dy = rf - static_cast<double>(r0);

    out_gx = (1.0 - dx) * (1.0 - dy) * gx00 +
             dx * (1.0 - dy) * gx10 +
             (1.0 - dx) * dy * gx01 +
             dx * dy * gx11;

    out_gy = (1.0 - dx) * (1.0 - dy) * gy00 +
             dx * (1.0 - dy) * gy10 +
             (1.0 - dx) * dy * gy01 +
             dx * dy * gy11;

    return finite(out_gx) && finite(out_gy);
}

DemMetadata read_dem_metadata(const std::filesystem::path& path) {
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.string().c_str(), GA_ReadOnly));
    if (ds == nullptr) {
        throw std::runtime_error("Failed to open raster: " + path.string());
    }

    DemMetadata metadata;
    metadata.width = ds->GetRasterXSize();
    metadata.height = ds->GetRasterYSize();
    if (ds->GetGeoTransform(metadata.gt) != CE_None) {
        GDALClose(ds);
        throw std::runtime_error("Failed to read geotransform: " + path.string());
    }
    if (const char* wkt = ds->GetProjectionRef(); wkt != nullptr) {
        metadata.projection_wkt = wkt;
    }
    GDALClose(ds);
    return metadata;
}

DemRaster read_dem_geotiff(const std::filesystem::path& path) {
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(path.string().c_str(), GA_ReadOnly));
    if (ds == nullptr) {
        throw std::runtime_error("Failed to open raster: " + path.string());
    }

    GDALRasterBand* band = ds->GetRasterBand(1);
    if (band == nullptr) {
        GDALClose(ds);
        throw std::runtime_error("Raster has no band 1: " + path.string());
    }

    DemRaster dem;
    dem.width = ds->GetRasterXSize();
    dem.height = ds->GetRasterYSize();
    if (ds->GetGeoTransform(dem.gt) != CE_None) {
        GDALClose(ds);
        throw std::runtime_error("Failed to read geotransform: " + path.string());
    }
    if (const char* wkt = ds->GetProjectionRef(); wkt != nullptr) {
        dem.projection_wkt = wkt;
    }

    int success = 0;
    dem.no_data = band->GetNoDataValue(&success);
    dem.has_no_data = (success != 0);

    dem.z.resize(static_cast<std::size_t>(dem.width) * static_cast<std::size_t>(dem.height));
    const CPLErr err = band->RasterIO(GF_Read,
                                      0,
                                      0,
                                      dem.width,
                                      dem.height,
                                      dem.z.data(),
                                      dem.width,
                                      dem.height,
                                      GDT_Float64,
                                      0,
                                      0);
    GDALClose(ds);
    if (err != CE_None) {
        throw std::runtime_error("RasterIO read failed: " + path.string());
    }

    dem.validate_north_up(path.filename().string());
    return dem;
}

void validate_horizontal_crs_compatible(const std::string& reference_wkt,
                                        const std::string& target_wkt) {
    if (reference_wkt.empty() || target_wkt.empty()) {
        throw std::runtime_error("Both DSMs must contain a CRS definition.");
    }

    OGRSpatialReference reference_srs;
    OGRSpatialReference target_srs;
    reference_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);
    target_srs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    if (reference_srs.SetFromUserInput(reference_wkt.c_str()) != OGRERR_NONE ||
        target_srs.SetFromUserInput(target_wkt.c_str()) != OGRERR_NONE) {
        throw std::runtime_error("Failed to parse one or both raster CRS definitions.");
    }

    // Remove vertical components before comparison.  The alignment model allows
    // a global tz offset, so differing vertical CRS realizations should not make
    // the horizontal co-registration tool reject otherwise compatible rasters.
    if (reference_srs.StripVertical() != OGRERR_NONE ||
        target_srs.StripVertical() != OGRERR_NONE) {
        throw std::runtime_error("Failed to isolate horizontal CRS components.");
    }

    if (!reference_srs.IsProjected() || !target_srs.IsProjected()) {
        throw std::runtime_error("DSM alignment requires projected horizontal coordinate systems.");
    }

    const double reference_units = reference_srs.GetLinearUnits(nullptr);
    const double target_units = target_srs.GetLinearUnits(nullptr);
    if (std::abs(reference_units - 1.0) > 1.0e-9 ||
        std::abs(target_units - 1.0) > 1.0e-9) {
        throw std::runtime_error("DSM alignment currently requires projected coordinates in meters.");
    }

    if (!reference_srs.IsSame(&target_srs)) {
        throw std::runtime_error(
            "Reference and target DSM horizontal CRS definitions do not match. "
            "Reproject both rasters to the same projected CRS before alignment.");
    }
}

DemRaster downsample_average(const DemRaster& input, int factor) {
    if (factor <= 1) return input;
    input.validate_north_up("downsample input");

    const int out_width = input.width / factor;
    const int out_height = input.height / factor;
    if (out_width < 2 || out_height < 2) {
        throw std::runtime_error("Pyramid factor " + std::to_string(factor) +
                                 " is too large for the raster dimensions.");
    }

    DemRaster output;
    output.width = out_width;
    output.height = out_height;
    output.projection_wkt = input.projection_wkt;
    output.has_no_data = true;
    output.no_data = std::numeric_limits<double>::quiet_NaN();
    copy_geotransform(input.gt, output.gt);
    output.gt[1] *= static_cast<double>(factor);
    output.gt[2] *= static_cast<double>(factor);
    output.gt[4] *= static_cast<double>(factor);
    output.gt[5] *= static_cast<double>(factor);
    output.z.assign(static_cast<std::size_t>(out_width) * static_cast<std::size_t>(out_height),
                    std::numeric_limits<double>::quiet_NaN());

    for (int r_out = 0; r_out < out_height; ++r_out) {
        for (int c_out = 0; c_out < out_width; ++c_out) {
            double sum = 0.0;
            int count = 0;
            const int r0 = r_out * factor;
            const int c0 = c_out * factor;

            for (int dr = 0; dr < factor; ++dr) {
                for (int dc = 0; dc < factor; ++dc) {
                    const double value = input.at(c0 + dc, r0 + dr);
                    if (input.valid_z(value)) {
                        sum += value;
                        ++count;
                    }
                }
            }

            if (count > 0) {
                output.z[static_cast<std::size_t>(r_out) * static_cast<std::size_t>(out_width) +
                         static_cast<std::size_t>(c_out)] = sum / static_cast<double>(count);
            }
        }
    }

    return output;
}

GradRaster compute_target_gradients(const DemRaster& target) {
    target.validate_north_up("target");

    GradRaster gradient;
    gradient.width = target.width;
    gradient.height = target.height;
    copy_geotransform(target.gt, gradient.gt);
    gradient.gx.assign(static_cast<std::size_t>(gradient.width) * gradient.height,
                       std::numeric_limits<double>::quiet_NaN());
    gradient.gy.assign(static_cast<std::size_t>(gradient.width) * gradient.height,
                       std::numeric_limits<double>::quiet_NaN());

    // IMPORTANT: use signed pixel sizes.  gt[5] is normally negative for a
    // north-up raster, and retaining that sign converts image-row differences
    // into the correct map-coordinate derivative dz/dy.
    const double dx = target.gt[1];
    const double dy = target.gt[5];

    for (int r = 0; r < gradient.height; ++r) {
        for (int c = 0; c < gradient.width; ++c) {
            if (!target.valid_z(target.at(c, r))) continue;

            const int c_left = std::max(0, c - 1);
            const int c_right = std::min(gradient.width - 1, c + 1);
            const double z_left = target.at(c_left, r);
            const double z_right = target.at(c_right, r);
            if (c_left != c_right && target.valid_z(z_left) && target.valid_z(z_right)) {
                gradient.gx[static_cast<std::size_t>(r) * gradient.width + c] =
                    (z_right - z_left) / (static_cast<double>(c_right - c_left) * dx);
            }

            const int r_up = std::max(0, r - 1);
            const int r_down = std::min(gradient.height - 1, r + 1);
            const double z_up = target.at(c, r_up);
            const double z_down = target.at(c, r_down);
            if (r_up != r_down && target.valid_z(z_up) && target.valid_z(z_down)) {
                gradient.gy[static_cast<std::size_t>(r) * gradient.width + c] =
                    (z_down - z_up) / (static_cast<double>(r_down - r_up) * dy);
            }
        }
    }

    return gradient;
}

void write_aligned_dsm(const std::filesystem::path& path,
                       const DemRaster& target,
                       double tx_m,
                       double ty_m,
                       double tz_m) {
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("GTiff");
    if (driver == nullptr) {
        throw std::runtime_error("GDAL GTiff driver is unavailable.");
    }

    char** options = nullptr;
    options = CSLSetNameValue(options, "TILED", "YES");
    options = CSLSetNameValue(options, "COMPRESS", "DEFLATE");
    options = CSLSetNameValue(options, "PREDICTOR", "3");
    options = CSLSetNameValue(options, "BIGTIFF", "IF_SAFER");

    GDALDataset* out = driver->Create(path.string().c_str(),
                                      target.width,
                                      target.height,
                                      1,
                                      GDT_Float32,
                                      options);
    CSLDestroy(options);
    if (out == nullptr) {
        throw std::runtime_error("Failed to create aligned DSM: " + path.string());
    }

    double shifted_gt[6];
    copy_geotransform(target.gt, shifted_gt);
    shifted_gt[0] += tx_m;
    shifted_gt[3] += ty_m;

    if (out->SetGeoTransform(shifted_gt) != CE_None ||
        out->SetProjection(target.projection_wkt.c_str()) != CE_None) {
        GDALClose(out);
        throw std::runtime_error("Failed to write aligned DSM georeferencing.");
    }

    GDALRasterBand* band = out->GetRasterBand(1);
    const float output_nodata = -9999.0F;
    band->SetNoDataValue(output_nodata);

    std::vector<float> row(static_cast<std::size_t>(target.width));
    for (int r = 0; r < target.height; ++r) {
        for (int c = 0; c < target.width; ++c) {
            const double value = target.at(c, r);
            row[static_cast<std::size_t>(c)] = target.valid_z(value)
                                                  ? static_cast<float>(value + tz_m)
                                                  : output_nodata;
        }

        if (band->RasterIO(GF_Write,
                           0,
                           r,
                           target.width,
                           1,
                           row.data(),
                           target.width,
                           1,
                           GDT_Float32,
                           0,
                           0) != CE_None) {
            GDALClose(out);
            throw std::runtime_error("Failed while writing aligned DSM rows.");
        }
    }

    band->SetDescription("Aligned target DSM");
    out->SetMetadataItem("DSM_COREG_TX_M", std::to_string(tx_m).c_str());
    out->SetMetadataItem("DSM_COREG_TY_M", std::to_string(ty_m).c_str());
    out->SetMetadataItem("DSM_COREG_TZ_M", std::to_string(tz_m).c_str());
    GDALClose(out);
}

}  // namespace dsm_coreg

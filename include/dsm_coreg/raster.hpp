#pragma once

#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace dsm_coreg {

struct DemMetadata {
    int width = 0;
    int height = 0;
    double gt[6] = {0, 1, 0, 0, 0, -1};
    std::string projection_wkt;
};

// In-memory DSM representation.  Elevations are stored as double precision so
// interpolation and optimization are not affected by repeated float rounding.
struct DemRaster {
    int width = 0;
    int height = 0;
    std::vector<double> z;
    double gt[6] = {0, 1, 0, 0, 0, -1};
    std::string projection_wkt;
    bool has_no_data = false;
    double no_data = std::numeric_limits<double>::quiet_NaN();

    void validate_north_up(const std::string& name) const;

    [[nodiscard]] bool in_bounds(int c, int r) const;
    [[nodiscard]] double at(int c, int r) const;
    [[nodiscard]] bool valid_z(double value) const;

    void pixel_center_to_world(int c, int r, double& x, double& y) const;
    void world_to_pixel_center(double x, double y, double& c, double& r) const;

    // Bilinear sampling in world coordinates.  The function requires all four
    // neighbors to be valid; this intentionally prevents NoData from being
    // silently interpolated into an apparently valid elevation.
    bool bilinear_sample(double x, double y, double& value) const;
};

struct GradRaster {
    int width = 0;
    int height = 0;
    std::vector<double> gx;  // dz / dx in map coordinates
    std::vector<double> gy;  // dz / dy in map coordinates
    double gt[6] = {0, 1, 0, 0, 0, -1};

    [[nodiscard]] bool in_bounds(int c, int r) const;
    [[nodiscard]] bool valid(double value) const;
    [[nodiscard]] double at_gx(int c, int r) const;
    [[nodiscard]] double at_gy(int c, int r) const;

    void world_to_pixel_center(double x, double y, double& c, double& r) const;
    bool bilinear_sample(double x, double y, double& out_gx, double& out_gy) const;
};

DemMetadata read_dem_metadata(const std::filesystem::path& path);
DemRaster read_dem_geotiff(const std::filesystem::path& path);

// Validate that both rasters use the same projected horizontal CRS and meter
// units.  Vertical CRS differences are intentionally ignored because tz is
// jointly estimated as a nuisance parameter.
void validate_horizontal_crs_compatible(const std::string& reference_wkt,
                                        const std::string& target_wkt);

// Average valid native pixels inside non-overlapping factor x factor cells.
// Partial cells at the far right/bottom edge are omitted so every output pixel
// represents the same physical footprint and the coarse-grid geotransform is
// exact.
DemRaster downsample_average(const DemRaster& input, int factor);

GradRaster compute_target_gradients(const DemRaster& target);

// Write a georeferenced target DSM after applying the estimated rigid
// translation.  Horizontal alignment is represented by shifting the raster
// geotransform; vertical alignment is applied by adding tz to valid elevations.
void write_aligned_dsm(const std::filesystem::path& path,
                       const DemRaster& target,
                       double tx_m,
                       double ty_m,
                       double tz_m);

}  // namespace dsm_coreg

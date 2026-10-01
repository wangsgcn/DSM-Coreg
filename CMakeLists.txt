cmake_minimum_required(VERSION 3.24)

project(dsm_coreg
    VERSION 0.4.0
    DESCRIPTION "Robust DSM co-registration with spatial block-bootstrap uncertainty"
    LANGUAGES CXX)

# -----------------------------------------------------------------------------
# Language and compiler requirements
# -----------------------------------------------------------------------------
#
# The code intentionally targets C++20.  GCC 14 is the recommended compiler for
# this project, but the compiler path is not hard-coded.  Select GCC 14 at CMake
# configure time, for example:
#
#   cmake -S . -B build \
#       -DCMAKE_CXX_COMPILER=g++-14 \
#       -DCMAKE_BUILD_TYPE=Release
#
# Keeping compiler selection outside the project makes the source portable to a
# workstation, Docker container, CI system, or another Linux distribution.
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# Required dependency: GDAL
# -----------------------------------------------------------------------------
# GDAL is used for GeoTIFF I/O, geotransforms, CRS validation, warping-related
# metadata, and JSON result serialization through CPLJSON.
find_package(GDAL REQUIRED)

# -----------------------------------------------------------------------------
# Required dependency: toml++
# -----------------------------------------------------------------------------
# toml++ is NOT part of the C++ standard library.  It is a third-party, modern
# C++ TOML parser.  The project first looks for a system/package-manager install
# that exports the CMake target tomlplusplus::tomlplusplus.
#
# If it is not installed, the default behavior is to fetch a pinned release with
# CMake FetchContent.  Pinning the version makes builds reproducible instead of
# silently following the newest upstream commit.
#
# For an offline/air-gapped environment, install toml++ beforehand and configure
# with:
#
#   -DDSM_COREG_FETCH_TOMLPLUSPLUS=OFF
#
# CMake will then fail clearly if no local toml++ package is available.
option(DSM_COREG_FETCH_TOMLPLUSPLUS
       "Fetch toml++ with CMake FetchContent when no local package is found"
       ON)

find_package(tomlplusplus CONFIG QUIET)

if(NOT TARGET tomlplusplus::tomlplusplus)
    if(DSM_COREG_FETCH_TOMLPLUSPLUS)
        include(FetchContent)
        FetchContent_Declare(
            tomlplusplus
            GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
            GIT_TAG        v3.4.0
            GIT_SHALLOW    TRUE
        )
        FetchContent_MakeAvailable(tomlplusplus)
    else()
        message(FATAL_ERROR
            "toml++ was not found. Install toml++ so CMake can find "
            "tomlplusplus::tomlplusplus, or reconfigure with "
            "-DDSM_COREG_FETCH_TOMLPLUSPLUS=ON on a machine with network access.")
    endif()
endif()

# -----------------------------------------------------------------------------
# Optional dependency: OpenMP
# -----------------------------------------------------------------------------
# If available, OpenMP parallelizes computationally expensive bootstrap work and
# the optional coarse-search stage.  The program still builds and runs serially
# when OpenMP is unavailable.
find_package(OpenMP QUIET)

# -----------------------------------------------------------------------------
# Executable
# -----------------------------------------------------------------------------
add_executable(dsm_coreg
    src/main.cpp
    src/config.cpp
    src/logging.cpp
    src/raster.cpp
    src/alignment.cpp
    src/bootstrap.cpp
    src/bootstrap_geopackage.cpp
    src/results.cpp
)

target_include_directories(dsm_coreg
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(dsm_coreg
    PRIVATE
        GDAL::GDAL
        tomlplusplus::tomlplusplus
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(dsm_coreg PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(dsm_coreg PRIVATE DSM_COREG_HAS_OPENMP=1)
endif()

# Make the program version available to --version and to the machine-readable
# result JSON without duplicating the version string in C++ source files.
target_compile_definitions(dsm_coreg
    PRIVATE
        DSM_COREG_VERSION="${PROJECT_VERSION}"
)

# Strong warnings are useful for scientific/research software.  They are not
# promoted to errors so a harmless compiler-version-specific warning does not
# prevent a valid build on another system.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(dsm_coreg PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Release builds are strongly recommended for large rasters and hundreds of
# bootstrap replicates.  Avoid forcing -march=native here because it would make
# the resulting executable less portable between machines.
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "GDAL version: ${GDAL_VERSION}")
message(STATUS "toml++ target: tomlplusplus::tomlplusplus")
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP: enabled")
else()
    message(STATUS "OpenMP: not found; building serial implementation")
endif()

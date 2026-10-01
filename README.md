# DSM-Coreg

**Robust co-registration and spatial block-bootstrap uncertainty estimation for high-resolution digital surface models**

DSM-Coreg is an open-source C++20 research software tool for estimating the relative translation between two already georeferenced digital surface models (DSMs).

The `dsm_coreg` command-line program jointly estimates horizontal and vertical offsets `(tx, ty, tz)` using sub-pixel DSM sampling, Huber-weighted Gauss-Newton optimization, and a coarse-to-fine DSM pyramid.

The tool is designed for cross-sensor and multi-temporal DSM pairs in which a global geolocation offset may coexist with large local surface differences caused by physical change, vegetation, occlusion, interpolation, reconstruction artifacts, or sensor-dependent surface representation.

Spatial block-bootstrap resampling is provided to quantify the precision and spatial stability of the recovered scene-wide translation while retaining local spatial dependence.

---

## Key Features

- Joint estimation of horizontal and vertical translation `(tx, ty, tz)`
- Sub-pixel bilinear sampling of the target DSM
- Huber robust loss for large local surface discrepancies
- Coarse-to-fine multi-resolution optimization
- Optional coarse horizontal XY search
- Spatial block-bootstrap uncertainty estimation
- Fixed and multiscale bootstrap block sizes
- Bootstrap standard deviations and percentile intervals
- Optional bootstrap GeoPackage audit output
- GeoTIFF input and aligned DSM output
- TOML configuration
- Machine-readable JSON results
- Detailed logging and runtime diagnostics
- Explicit CRS and raster validation
- OpenMP-based multi-threaded processing

---

## Intended Use

DSM-Coreg is intended for relative co-registration of already georeferenced DSMs whose dominant remaining geometric discrepancy can be represented by a global translation.

Example applications include:

- lidar-to-photogrammetric DSM alignment
- multi-temporal DSM comparison
- change-detection preprocessing
- DSM validation
- DSM fusion
- cross-sensor surface comparison

The current model estimates only a global translation. It does **not** estimate absolute geolocation accuracy, rotation, scale, shear, or spatially varying warping.

---

## Quick Start

A normal run is:

```bash
./build/dsm_coreg --config config/dsm_coreg.example.toml
```

Validate the configuration, input paths, raster metadata, projected meter units, and CRS without running the alignment:

```bash
./build/dsm_coreg \
  --config config/dsm_coreg.example.toml \
  --validate-only
```

Show the program version:

```bash
./build/dsm_coreg --version
```

---

## Source Layout

```text
CMakeLists.txt

config/
  dsm_coreg.example.toml

include/dsm_coreg/
  alignment.hpp
  bootstrap.hpp
  bootstrap_geopackage.hpp
  config.hpp
  logging.hpp
  raster.hpp
  results.hpp
  timing.hpp

src/
  alignment.cpp
  bootstrap.cpp
  bootstrap_geopackage.cpp
  config.cpp
  logging.cpp
  main.cpp
  raster.cpp
  results.cpp
```

---

## Requirements

Required:

- C++20 compiler
- CMake
- GDAL development package
- toml++

Optional:

- OpenMP for multi-threaded processing

GCC 14 is recommended for the current development configuration.

---

## Build

### Standard Build

If toml++ is already installed, CMake uses the installed package. Otherwise, the project can fetch a pinned toml++ release automatically.

```bash
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=g++-14 \
  -DCMAKE_BUILD_TYPE=Release

cmake --build build -j
```

The executable is:

```text
build/dsm_coreg
```

### Offline / Air-Gapped Build

Install toml++ locally and disable automatic downloading:

```bash
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=g++-14 \
  -DCMAKE_BUILD_TYPE=Release \
  -DDSM_COREG_FETCH_TOMLPLUSPLUS=OFF
```

If toml++ is installed in a non-standard prefix:

```bash
cmake -S . -B build \
  -DCMAKE_CXX_COMPILER=g++-14 \
  -DCMAKE_BUILD_TYPE=Release \
  -DDSM_COREG_FETCH_TOMLPLUSPLUS=OFF \
  -DCMAKE_PREFIX_PATH=/opt/tomlplusplus
```

---

## Configuration

DSM-Coreg uses TOML for human-editable configuration and JSON for machine-readable results.

A typical run is:

```bash
./build/dsm_coreg --config config/dsm_coreg.example.toml
```

Relative file paths in the TOML configuration are resolved relative to the configuration file itself rather than the shell working directory.

Unknown configuration keys are rejected deliberately to reduce the risk of silent mistakes in long scientific runs.

For example, a misspelled key such as:

```toml
huber_dleta_m = 1.0
```

produces an error instead of silently falling back to a default value.

### `[input]`

```toml
[input]
reference_dsm = "reference.tif"
target_dsm = "target.tif"
```

- `reference_dsm`: designated reference DSM
- `target_dsm`: DSM to be translated into alignment with the reference

### `[alignment]`

Important parameters include:

```toml
[alignment]
initial_tx_m = 0.0
initial_ty_m = 0.0
initial_tz_m = 0.0

border_px = 5
max_abs_dz_m = 50.0
huber_delta_m = 1.0

max_iterations = 30
step_tolerance_m = 0.001
lambda = 0.001
min_samples = 200
```

- `initial_tx_m`, `initial_ty_m`, `initial_tz_m`: initial translation estimate
- `border_px`: reference-image border excluded from fitting
- `max_abs_dz_m`: optional current-alignment residual screen; set `<= 0` to disable
- `huber_delta_m`: Huber threshold; set `<= 0` for ordinary least squares
- `max_iterations`: maximum Gauss-Newton iterations per pyramid level
- `step_tolerance_m`: convergence threshold on the 3-D translation update
- `lambda`: diagonal damping term
- `min_samples`: minimum number of participating samples

### `[pyramid]`

```toml
[pyramid]
enabled = true
factors = [8, 4, 2, 1]
```

The pyramid is processed from coarse to fine.

Translations are represented in map units, so the same `(tx, ty, tz)` estimate passes directly from one resolution level to the next.

### `[coarse_search]`

The optional coarse XY search can initialize the nonlinear optimizer when the initial horizontal displacement is too large for reliable local convergence.

```toml
[coarse_search]
enabled = false
range_x_m = 20.0
range_y_m = 20.0
step_m = 2.0
sample_stride = 4
min_samples = 200
```

The coarse search is used only for initialization; the final translation is estimated by the robust nonlinear optimizer.

---

## Spatial Block Bootstrap

The bootstrap is applied after the final native-resolution alignment.

For each replicate:

1. the valid reference-domain samples are partitioned into spatial blocks;
2. blocks are sampled with replacement;
3. the complete scene-wide translation `(tx, ty, tz)` is re-estimated;
4. one global translation vector is recorded for that replicate.

The bootstrap therefore estimates the sampling variability of the **global co-registration solution**, not independent local translations.

### Fixed Block Mode

```toml
[bootstrap]
enabled = true
replicates = 200
block_mode = "fixed"
block_px = 64
seed = 12345
```

`block_px` is expressed in **reference-DSM pixels**.

For example:

- a 64-pixel block on a 0.25 m reference grid corresponds to 16 m;
- a 32-pixel block on a 0.5 m reference grid also corresponds to 16 m.

### Multiscale Block Mode

```toml
[bootstrap]
enabled = true
replicates = 200
block_mode = "multiscale"

block_min_px = 32
block_max_px = 128
block_step_px = 16

seed = 12345
```

Multiscale mode evaluates robustness across a range of spatial resampling scales.

Because the ensemble combines both block-selection variability and block-size variability, it should be interpreted as a robustness analysis rather than as a direct substitute for uncertainty at one fixed spatial scale.

### Bootstrap Interpretation

Bootstrap standard deviations and percentile intervals characterize the precision and spatial stability of the relative scene-wide translation under spatial resampling.

They do **not** represent the absolute geolocation accuracy of either input DSM.

---

## Optional Bootstrap GeoPackage

For auditing or debugging, the selected bootstrap blocks can be written to a GeoPackage:

```toml
[bootstrap.output]
write_geopackage = true
geopackage = "bootstrap_replicates.gpkg"
```

The GeoPackage contains one polygon layer per requested replicate:

```text
bootstrap_000001
bootstrap_000002
...
```

Each layer contains the unique blocks selected for that replicate, with attributes such as:

- `block_id`
- `block_row`
- `block_col`
- `block_px`
- `multiplicity`
- `n_samples`

A block selected multiple times appears once with an increased `multiplicity`.

The GeoPackage also contains a non-spatial `bootstrap_replicates` table with one row per replicate containing the global translation estimate, horizontal shift, weighted RMSE, convergence state, and sampling counts.

This output is optional because writing hundreds of layers can add substantial I/O time and disk usage.

---

## Sign Convention

For a reference sample at `(x, y)`, DSM-Coreg uses the residual

```text
r = z_ref(x, y) - [ z_target(x - tx, y - ty) + tz ]
```

Therefore:

- positive `tx` shifts the target DSM in the positive map-x direction;
- positive `ty` shifts the target DSM in the positive map-y direction;
- `tz` is added to valid target elevations.

The aligned target surface is

```text
z_aligned(x, y) = z_target(x - tx, y - ty) + tz
```

---

## Outputs

A typical output directory contains:

```text
alignment_results.json
alignment.log
effective_config.toml
target_aligned.tif
```

### `alignment_results.json`

Contains machine-readable results such as:

- final `(tx, ty, tz)`
- horizontal translation magnitude
- weighted RMSE
- number of participating samples
- convergence state
- per-pyramid-level results
- bootstrap means
- bootstrap standard deviations
- percentile intervals
- timing information

### `effective_config.toml`

Records the complete resolved configuration used for the run, including defaults and resolved file paths.

It can be archived with an experiment or reused as the starting point for a later run.

### `alignment.log`

Contains processing progress, warnings, convergence diagnostics, and timing information.

---

## Model Scope and Limitations

DSM-Coreg assumes that the dominant relative misregistration between the two input DSMs can be represented by a constant 3-D translation.

This assumption is appropriate when both DSMs are already georeferenced in a common projected coordinate system and the remaining discrepancy is dominated by translation.

The current version does not estimate:

- rotation
- scale
- shear
- nonrigid deformation
- spatially varying geolocation error

Large local DSM differences caused by construction, demolition, vegetation, occlusion, interpolation, or reconstruction differences are not expected to disappear after alignment. The purpose of the robust estimator is to recover the dominant shared geometric translation without forcing all local surfaces into agreement.

---

## Research Background

DSM-Coreg was developed to support research on robust co-registration of high-resolution lidar-derived and satellite-photogrammetric DSMs.

The method combines:

- scene-wide translation estimation;
- Huber-weighted robust optimization;
- coarse-to-fine alignment;
- spatial block-bootstrap uncertainty estimation;
- bidirectional consistency analysis; and
- block-size sensitivity analysis.

A manuscript describing the method and experimental evaluation is in preparation.

### Citation

If you use DSM-Coreg in published research, please cite the associated paper once a formal citation is available.

A BibTeX entry will be added here after publication.

---

## License

DSM-Coreg is released under the **MIT License**.

You are free to use, copy, modify, merge, publish, distribute, sublicense, and
sell copies of the software, subject to the conditions of the MIT License.

Commercial use is permitted.

See the repository's `LICENSE` file for the complete license text.

---

## Contributing

Bug reports, reproducible test cases, documentation improvements, and
contributions are welcome.

Before submitting substantial changes, please open an issue to discuss the
proposed work.

---

## Disclaimer

DSM-Coreg is provided under the MIT License without warranty of any kind.
Users are responsible for validating the suitability of the software and its
outputs for their own applications.

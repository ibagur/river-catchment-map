# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Pipeline

```r
# Full pipeline (from R or RStudio, working directory = repo root)
source("run_analysis.R")

# From terminal
Rscript run_analysis.R

# Regenerate only the map (fastest for styling iterations — already the default)
# In run_analysis.R, force_rerun["map"] is set to TRUE by default.
# Set other flags to TRUE only when upstream data changes.
```

## Global Parameters (in `run_analysis.R`)

All tunable parameters live at the top of `run_analysis.R`:

| Parameter | Default | Purpose |
|---|---|---|
| `EPSG_TARGET` | `32737` | UTM Zone 37S — all spatial data is projected to this |
| `CONTOUR_INT` | `50` | Contour interval in metres |
| `FLOW_THRESH` | `5000` | Min flow accumulation cells for pour point detection |
| `POUR_SNAP_DIST` | `100` | Metres to snap pour points to stream network |
| `DISTRICT_NAME` | `"Mecufi"` | Admin boundary filter string (case-insensitive) |

Each step can be individually re-triggered via `force_rerun` flags.

## Architecture

Six sourced scripts execute sequentially from `run_analysis.R`:

```
run_analysis.R
├── R/00_setup.R        — packages, WhiteboxTools init, dirs, shared helpers
├── R/01_load_data.R    — admin boundary, DEM tiles, waterways, roads → gis/derived/
├── R/02_dem_processing.R — reproject, sink-fill, hillshade, slope, aspect
├── R/03_hydrology.R    — D8 flow direction/accumulation, pour points, watersheds
├── R/04_contours.R     — contour generation, smoothing, major/minor classification
└── R/05_map.R          — tmap v4 cartographic composition → output/
```

### Caching pattern
Every script checks for its output files in `gis/derived/` before running. If they exist and the corresponding `force_rerun` flag is `FALSE`, the script loads from cache and skips processing. This means DEM processing (~2 min) and hydrology (~3 min) only re-run when needed.

### Key data flow variables
Scripts communicate through shared R objects (no explicit `return()`). After each script, the pipeline environment holds:

- `mecufi_utm` — district boundary sf (UTM 37S)
- `dem_filled` — sink-filled COP30 DEM SpatRaster (UTM 37S)
- `flow_acc` — D8 flow accumulation SpatRaster
- `watersheds_utm` — delineated catchment polygons sf
- `contours_utm` — smoothed contour lines sf with `is_major` column
- `waterways_lines_utm`, `waterways_polys_utm`, `roads_utm` — HOT OSM features sf

### Helper functions (`R/00_setup.R`)
- `assert_crs(x, epsg)` — errors if CRS doesn't match
- `assert_overlap(raster, vector)` — errors if extents don't overlap (handles CRS reprojection)
- `assert_rows(x, n_min)` — errors if object has too few features
- `wbt_path(...)` — wraps `normalizePath()` for WhiteboxTools (requires absolute paths)

## Input Data Layout

```
gis/
├── admin/               moz_admin_boundaries.shp.zip          (OCHA; ADM2 layer = districts)
├── raster/
│   ├── elevation/       rasters_COP30.tar.gz                  (Copernicus GLO-30 DEM tiles)
│   │                    rasters_AW3D30.tar.gz                 (ALOS AW3D30 DEM — unused)
│   └── population/      MOZ_population_v1_1_gridded.tif       (Maxar building-footprint pop.)
│                        moz_pop_2026_CN_100m_R2025A_v1.tif    (WorldPop 2026, 100 m — unused)
├── waterways/           hotosm_moz_waterways_lines_shp.zip
│                        hotosm_moz_waterways_polygons_shp.zip
└── roads/               moz_roads_shp.zip
```

All archives are read directly via `/vsizip/` (shapefiles) or `untar()` (rasters) — no manual extraction needed unless `gis/derived/` cache is empty.

## Output

```
output/
├── mecufi_catchment_map.png   — 300 dpi, 250×200 mm (primary deliverable)
├── mecufi_contours.shp        — smoothed contour lines
├── mecufi_watersheds.shp      — delineated catchment polygons
└── mecufi_flow_acc.tif        — D8 flow accumulation (DEFLATE compressed)
```

## Key Dependency Notes

- **tmap v4 required** — `00_setup.R` throws an error on v3. Install: `install.packages("tmap")`
- **WhiteboxTools** — auto-installed via `whitebox::install_whitebox()` if missing
- **WhiteboxTools requires absolute paths** — always use `wbt_path()` wrapper, never relative paths
- **D8 flow accumulation input** — `wbt_d8_flow_accumulation()` takes the filled DEM directly (NOT the flow direction pointer), with `out_type = "cells"`

## Map Visualization Details

The river glow effect in `05_map.R` uses:
1. Channel mask at `flow_thresh_display = 1500` cells
2. `terra::distance()` to compute metres-to-nearest-channel per cell
3. Exponential decay with `decay_radius_m = 250` m
4. Log-normalised flow accumulation weighting (70% weight) to distinguish main stems from tributaries

Hypsometric breaks are hardcoded to 0–375 m (25 m bands) matching the Mecufi DEM range. Update `hyps_breaks` and `hyps_colors` vectors in `05_map.R` if running on a different extent.

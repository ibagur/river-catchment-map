# River Catchment Map — Mecufi District, Cabo Delgado, Mozambique

An R pipeline that produces publication-quality river catchment maps for Mecufi District, combining hydrological analysis (D8 flow routing, watershed delineation) with cartographic composition in tmap v4. A flood exposure layer overlays population density weighted by drainage proximity, producing two map variants that compare different population datasets.

**Primary outputs** (300 dpi, 250 × 200 mm):
- `output/mecufi_catchment_map_real.png` — flood exposure using Maxar building-footprint counts
- `output/mecufi_catchment_map_estimated.png` — flood exposure using WorldPop 2026 RF model

---

## Quick Start

```r
# From R or RStudio (working directory = repo root)
source("run_analysis.R")
```

```bash
# From the terminal
Rscript run_analysis.R
```

For map styling iterations only (fastest — skips hydrology and DEM steps):

```r
# In run_analysis.R, force_rerun["map"] is TRUE by default.
# All other flags default to FALSE, so only the map step reruns.
source("run_analysis.R")
```

---

## Requirements

### R packages

```r
install.packages(c(
  "sf", "terra", "whitebox", "tidyterra",
  "tmap",      # >= 4.0 required
  "smoothr", "dplyr", "stringr", "scales"
))
```

**tmap v4 is required.** The setup script will error on v3.

### WhiteboxTools

Installed automatically on first run via `whitebox::install_whitebox()`. No manual setup needed.

---

## Input Data

Place source files in `gis/` before running. The pipeline reads archives directly — no manual extraction needed.

| Folder | File | Source |
|--------|------|--------|
| `gis/admin/` | `moz_admin_boundaries.shp.zip` | OCHA — ADM2 district boundaries |
| `gis/raster/elevation/` | `rasters_COP30.tar.gz` | Copernicus GLO-30 DEM tiles |
| `gis/raster/population/` | `MOZ_population_v1_1_gridded.tif` | Maxar — building-footprint constrained counts |
| `gis/raster/population/` | `moz_pop_2026_CN_100m_R2025A_v1.tif` | WorldPop 2026 RF model, 100 m |
| `gis/waterways/` | `hotosm_moz_waterways_lines_shp.zip` | HOT OSM |
| `gis/waterways/` | `hotosm_moz_waterways_polygons_shp.zip` | HOT OSM |
| `gis/roads/` | `hotosm_moz_roads_lines_shp.zip` | HOT OSM |
| `gis/places/` | `hotosm_moz_populated_places_points_shp.zip` | HOT OSM |

---

## Pipeline Architecture

Six scripts execute sequentially from `run_analysis.R`:

```
run_analysis.R
├── R/00_setup.R          packages, WhiteboxTools init, dirs, shared helpers
├── R/01_load_data.R      admin boundary, DEM tiles, waterways, roads → gis/derived/
├── R/02_dem_processing.R reproject, sink-fill, hillshade, slope, aspect
├── R/03_hydrology.R      D8 flow direction/accumulation, pour points, watersheds
├── R/04_contours.R       contour generation, smoothing, major/minor classification
└── R/05_map.R            tmap v4 cartographic composition → output/
```

### Caching

Every script checks for its outputs in `gis/derived/` before running. If they exist and the corresponding `force_rerun` flag is `FALSE`, the step is skipped. DEM processing (~2 min) and hydrology (~3 min) only re-run when needed.

```r
force_rerun <- c(
  load_data = FALSE,
  dem_proc  = FALSE,
  hydrology = FALSE,
  contours  = FALSE,
  map       = TRUE    # always regenerate for styling iterations
)
```

---

## Global Parameters

All tunable parameters live at the top of `run_analysis.R`:

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `EPSG_TARGET` | `32737` | UTM Zone 37S — all spatial data projected to this |
| `CONTOUR_INT` | `50` | Contour interval in metres |
| `FLOW_THRESH` | `5000` | Min flow accumulation cells for pour point detection |
| `POUR_SNAP_DIST` | `100` | Metres to snap pour points to stream network |
| `DISTRICT_NAME` | `"Mecufi"` | Admin boundary filter string (case-insensitive) |

---

## Outputs

| File | Description |
|------|-------------|
| `output/mecufi_catchment_map_real.png` | Map with Maxar population exposure — 300 dpi, 250 × 200 mm |
| `output/mecufi_catchment_map_estimated.png` | Map with WorldPop 2026 exposure — 300 dpi, 250 × 200 mm |
| `output/mecufi_contours.shp` | Smoothed contour lines with `is_major` column |
| `output/mecufi_watersheds.shp` | Delineated catchment polygons |
| `output/mecufi_flow_acc.tif` | D8 flow accumulation (DEFLATE compressed) |

Intermediate processing files are cached in `gis/derived/` (GeoPackage and GeoTIFF format).

---

## Map Composition

The map layers, rendered bottom to top:

1. **Hillshade** — Copernicus GLO-30, single light source
2. **Hypsometric DEM** — 15 discrete 25 m elevation bands (0–375 m), semi-transparent
3. **River glow** — distance-decay halo centred on channel pixels; main stems glow wider than tributaries via log-normalised flow accumulation weighting
4. **Contours** — minor (50 m, light grey) and major (250 m, darker grey)
5. **Watershed polygons** — semi-transparent blue fills with steel-blue borders
6. **Flood exposure** — population × drainage proximity score (YlOrRd, four discrete classes); varies between map variants
7. **Waterway polygons** — solid water bodies (lakes, wide rivers)
8. **Waterway lines** — stream and river network
9. **Roads** — HOT OSM road network
10. **Populated places** — point markers and name labels
11. **District boundary** — dark grey border, always rendered on top

### Dual-output architecture

`05_map.R` wraps the full layer stack in a `build_map()` helper parameterised by the exposure raster, legend label, and population credits suffix. Both variants are produced in a single pipeline run; all other layers are shared.

### Flood exposure methodology

```r
# 1. Proximity: exponential decay from drainage channels (e-folding = 400 m)
flood_channels <- terra::ifel(flow_acc >= 500L, 1L, NA)
dist_to_flood  <- terra::distance(flood_channels)
flood_prox     <- exp(-dist_to_flood / 400)

# 2. Raw exposure: population × proximity
exposure_raw <- pop_resampled * flood_prox

# 3. Log-normalise to [0, 1] relative to district maximum
exposure_norm <- log1p(exposure_raw) / max(log1p(exposure_raw))
```

The Maxar variant masks to cells with confirmed building footprints; the WorldPop variant applies a population threshold (>= 1 person per cell) plus a bottom-quartile cut to suppress background noise from the continuous RF model.

### River glow technique

```r
# Threshold flow accumulation → binary channel mask
channels_mask <- terra::ifel(flow_acc >= 1500L, 1L, NA)

# Distance to nearest channel pixel (metres, UTM)
dist_to_channel <- terra::distance(channels_mask)

# Exponential spatial decay (e-folding radius = 250 m)
glow_spatial <- exp(-dist_to_channel / 250)

# Log-normalised flow magnitude weight → main stems glow brighter
flow_norm <- log1p(flow_acc) / max(log1p(flow_acc))

# Blend: 30% flat glow + 70% flow-weighted
river_glow <- glow_spatial * (0.30 + 0.70 * flow_norm)
```

### Hypsometric colour scheme

Hardcoded to the Mecufi DEM range (0–375 m). Update `hyps_breaks` and `hyps_colors` in `R/05_map.R` if running on a different area.

---

## Helper Functions (`R/00_setup.R`)

| Function | Purpose |
|----------|---------|
| `assert_crs(x, epsg)` | Errors if CRS does not match the expected EPSG code |
| `assert_overlap(raster, vector)` | Errors if extents do not overlap (handles CRS reprojection) |
| `assert_rows(x, n_min)` | Errors if object has fewer than `n_min` features |
| `wbt_path(...)` | Wraps `normalizePath()` — required because WhiteboxTools needs absolute paths |

---

## Data Sources

| Dataset | Provider | Licence |
|---------|----------|---------|
| Copernicus GLO-30 DEM | ESA / Copernicus | Free for non-commercial use |
| Admin boundaries (ADM2) | OCHA | HDX Open Data |
| Waterways, roads, places | HOT OpenStreetMap | ODbL |
| Population v1.1 (building-footprint) | Maxar / GRID3 | CC BY 4.0 |
| Population 2026 RF model | WorldPop | CC BY 4.0 |

---

## Notes

- **WhiteboxTools requires absolute paths** — always use the `wbt_path()` wrapper, never relative paths.
- **D8 flow accumulation input** — `wbt_d8_flow_accumulation()` takes the filled DEM directly (not the flow direction pointer), with `out_type = "cells"`.
- The pipeline targets Mecufi District by filtering the OCHA ADM2 boundary on `DISTRICT_NAME`. To adapt for another district, change that parameter and supply a matching DEM covering the new extent.

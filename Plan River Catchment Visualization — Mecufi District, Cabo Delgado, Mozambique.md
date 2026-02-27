# Updated Plan: River Catchment Visualization — Mecufi District, Cabo Delgado, Mozambique

## Input Files (from your data)

| Category                   | Files to Use                            | Format Choice                                 |
| -------------------------- | --------------------------------------- | --------------------------------------------- |
| Admin boundary             | `moz_admin_boundaries.shp.zip`          | Shapefile (simplest for R)                    |
| DEM (primary)              | `rasters_COP30.tar.gz`                  | Copernicus GLO-30 (better quality for Africa) |
| DEM (secondary/validation) | `rasters_AW3D30.tar.gz`                 | JAXA — use as cross-check if needed           |
| Waterways (lines)          | `hotosm_moz_waterways_lines_shp.zip`    | Shapefile                                     |
| Waterways (polygons)       | `hotosm_moz_waterways_polygons_shp.zip` | Shapefile (for larger water bodies)           |

> **COP30 preferred over AW3D30** for this region — Copernicus GLO-30 has better void-filling and is more consistent over coastal Mozambique terrain.

------

## Revised Step-by-Step Workflow

### Step 0 — Environment Setup & File Extraction

- Extract all `.zip` and `.tar.gz` archives
- Install/load: `sf`, `terra`, `whitebox`, `tidyterra`, `tmap`, `smoothr`
- Initialize WhiteboxTools: `whitebox::wbt_init()`

### Step 1 — Load and Filter to Mecufi District

- Load `moz_admin_boundaries.shp` with `sf::st_read()`
- **Filter** to Mecufi district: `dplyr::filter(boundaries, district == "Mecufi")` — verify the exact attribute field name
- Load COP30 DEM raster with `terra::rast()`
- Reproject all layers to **UTM Zone 37S (EPSG:32737)** — appropriate projected CRS for Cabo Delgado
- Clip DEM to Mecufi boundary: `terra::crop()` + `terra::mask()`
- Load and clip waterway lines and polygons to Mecufi extent

### Step 2 — DEM Preprocessing & Terrain Derivatives

- **Sink filling**: `whitebox::wbt_fill_depressions()` on the clipped COP30 DEM — essential before any hydrological analysis
- **Hillshade**: `terra::shade(slope, aspect, angle=45, direction=315)` — produces the grey relief background
- **Slope**: `terra::terrain(dem, v="slope")`
- **Aspect**: `terra::terrain(dem, v="aspect")`
- **Hypsometric color ramp**: classify DEM values into elevation bands with a brown→yellow→green color ramp (matching the example image style)

### Step 3 — Contour Lines

- Generate contours from filled DEM: `terra::as.contour(dem_filled, nlevels=...)`
- Suggested interval: **50m or 100m** given Mecufi's coastal-to-inland relief (likely modest elevation range ~0–600m)
- Convert to `sf`, clip to district boundary
- Optionally smooth with `smoothr::smooth(method="ksmooth")`

### Step 4 — Hydrological Analysis & Catchment Delineation

Using **WhiteboxTools** pipeline on the sink-filled DEM:

1. `wbt_d8_pointer()` → flow direction raster
2. `wbt_d8_flow_accumulation()` → upstream contributing area per cell
3. Extract **pour points**: use waterway line endpoints where rivers exit the district, or manually define outlet on the coast
4. `wbt_snap_pour_points()` → snap outlets to highest accumulation cell within tolerance
5. `wbt_watershed()` → delineate catchment raster per outlet
6. Convert catchment raster → vector polygon: `terra::as.polygons()` → `sf::st_as_sf()`

> For Mecufi specifically: rivers likely drain east toward the Indian Ocean coast — confirm outlet direction before setting pour points.

### Step 5 — River Shading / Catchment Influence Layer

Two complementary approaches (combine both for best visual result):

**A) Flow accumulation shading** (continuous gradient effect):

- Log-transform: `log1p(flow_accumulation)`
- Map to a blue color ramp, overlay semi-transparently on hillshade
- This naturally creates wider shading along major rivers, thinning toward headwaters — matching the example

**B) Buffered river zones** (explicit catchment shading):

- Multi-ring buffers around waterway lines at e.g. 500m, 1km, 2km with decreasing opacity
- `sf::st_buffer()` + alpha transparency in `tmap`

### Step 6 — Final Map Composition with `tmap`

```
Layer stack (bottom → top):
1. Hillshade raster          — grey relief base
2. Hypsometric DEM           — elevation color ramp (brown/yellow/green)
3. Flow accumulation shading — blue gradient along rivers
4. Contour lines             — thin grey/white, semi-transparent
5. Catchment polygons        — outlined, lightly filled per sub-basin
6. Waterway polygons         — filled blue (lakes, wide rivers)
7. Waterway lines            — blue lines, width scaled to stream order
8. Mecufi district boundary  — bold red/dark outline
9. Map furniture             — north arrow, scale bar, graticule, legend
```

------

## Key Parameters to Confirm Before Running

| Parameter                   | Suggestion                 | To Verify                                  |
| --------------------------- | -------------------------- | ------------------------------------------ |
| Contour interval            | 50m                        | Check actual elevation range of Mecufi DEM |
| UTM zone                    | EPSG:32737 (37S)           | Confirm district falls within this zone    |
| Admin filter field          | `"district"` or `"NAME_2"` | Inspect shapefile attribute table          |
| Pour point location         | Coast / river mouths       | Inspect waterway line endpoints            |
| Flow accumulation threshold | 1000–5000 cells            | Tune to match visible stream density       |
| DEM choice                  | COP30 primary              | Compare with AW3D30 if artifacts appear    |

------

## Output Files to Generate

- `mecufi_catchment_map.png` — final high-res map
- `mecufi_contours.shp` — contour lines
- `mecufi_watersheds.shp` — delineated catchment polygons
- `mecufi_flow_acc.tif` — flow accumulation raster (useful for future analysis)
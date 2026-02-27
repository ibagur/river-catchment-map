# ______________________________________________________________________________
# 03_hydrology.R — D8 flow routing, pour point detection, watershed delineation
# ______________________________________________________________________________
# Purpose : Compute D8 flow direction and accumulation, automatically identify
#           pour points at the coastal district boundary, snap them to the
#           stream network, and delineate watersheds.
# Inputs  : gis/derived/mecufi_dem_filled.tif
#           gis/derived/mecufi_boundary_utm.gpkg
#           gis/derived/mecufi_waterways_lines_utm.gpkg  (fallback pour points)
# Outputs : gis/derived/mecufi_flow_dir.tif
#           gis/derived/mecufi_flow_acc.tif
#           gis/derived/mecufi_pour_points.gpkg
#           gis/derived/mecufi_pour_points_snapped.shp   (WBT requires .shp)
#           gis/derived/mecufi_watersheds_raw.tif
#           gis/derived/mecufi_watersheds.gpkg
# ______________________________________________________________________________

# ______________________________________________________________________________
# CACHE CHECK ----
# ______________________________________________________________________________

derived_files_03 <- c(
  "gis/derived/mecufi_flow_acc.tif",
  "gis/derived/mecufi_watersheds.gpkg"
)

if (all(file.exists(derived_files_03)) && !force_rerun["hydrology"]) {
  message("03_hydrology: cache hit — loading from gis/derived/")

  flow_acc        <- terra::rast("gis/derived/mecufi_flow_acc.tif")
  watersheds_utm  <- sf::st_read("gis/derived/mecufi_watersheds.gpkg", quiet = TRUE)

} else {

# ______________________________________________________________________________
# LOAD INPUTS ----
# ______________________________________________________________________________

dem_filled <- terra::rast("gis/derived/mecufi_dem_filled.tif")
mecufi_utm <- sf::st_read("gis/derived/mecufi_boundary_utm.gpkg", quiet = TRUE)

filled_path <- wbt_path("gis/derived/mecufi_dem_filled.tif")

# ______________________________________________________________________________
# D8 FLOW DIRECTION ----
# ______________________________________________________________________________

message("03_hydrology: computing D8 flow direction...")

flow_dir_path <- wbt_path("gis/derived/mecufi_flow_dir.tif")

whitebox::wbt_d8_pointer(
  dem    = filled_path,
  output = flow_dir_path
)

# ______________________________________________________________________________
# D8 FLOW ACCUMULATION ----
# ______________________________________________________________________________
# CRITICAL: input is the filled DEM, NOT the flow direction pointer.

message("03_hydrology: computing D8 flow accumulation...")

flow_acc_path <- wbt_path("gis/derived/mecufi_flow_acc.tif")

whitebox::wbt_d8_flow_accumulation(
  input    = filled_path,
  output   = flow_acc_path,
  out_type = "cells"
)

flow_acc <- terra::rast("gis/derived/mecufi_flow_acc.tif")

# ______________________________________________________________________________
# POUR POINT IDENTIFICATION (automated) ----
# ______________________________________________________________________________

message("03_hydrology: identifying pour points at district boundary...")

## ---- Create 500 m interior ring ----
# Negative buffer creates an inner ring; difference gives the boundary strip
mecufi_inner   <- sf::st_buffer(mecufi_utm, -500)
boundary_strip <- sf::st_difference(mecufi_utm, mecufi_inner)
boundary_vect  <- terra::vect(boundary_strip)

## ---- Mask flow accumulation to boundary strip ----
flow_acc_boundary <- terra::mask(flow_acc, boundary_vect)

## ---- Apply FLOW_THRESH as absolute minimum, then thin to major outlets ----
boundary_vals <- terra::values(flow_acc_boundary, na.rm = TRUE)

if (length(boundary_vals) == 0) {
  warning("No flow accumulation values found in boundary strip — using whole district")
  boundary_vals <- terra::values(flow_acc, na.rm = TRUE)
  flow_acc_boundary <- flow_acc
}

# Step 1: apply absolute minimum accumulation threshold (FLOW_THRESH)
# This filters to cells that represent meaningful river channels
high_flow <- flow_acc_boundary >= FLOW_THRESH

n_above_thresh <- sum(terra::values(high_flow, na.rm = TRUE), na.rm = TRUE)
message(sprintf("  %d boundary cells above FLOW_THRESH (%d)", n_above_thresh, FLOW_THRESH))

# If threshold too high, fall back to 95th percentile
if (n_above_thresh < 3) {
  thresh_95 <- quantile(boundary_vals, 0.95, na.rm = TRUE)
  high_flow  <- flow_acc_boundary >= thresh_95
  message(sprintf("  FLOW_THRESH too high — using 95th percentile (%.0f cells)", thresh_95))
}

## ---- Convert high-flow cells to points ----
pour_pts_raw <- terra::as.points(terra::mask(flow_acc_boundary, high_flow),
                                  values = TRUE, na.rm = TRUE)

if (length(pour_pts_raw) == 0 || nrow(pour_pts_raw) < 1) {
  message("Automated pour point detection failed — using waterway line endpoints as fallback")

  ## ---- FALLBACK: downstream endpoints of OSM waterway lines ----
  ww_lines <- sf::st_read("gis/derived/mecufi_waterways_lines_utm.gpkg", quiet = TRUE)

  endpoints <- ww_lines %>%
    dplyr::rowwise() %>%
    dplyr::mutate(geometry = sf::st_point(
      tail(sf::st_coordinates(geometry)[, 1:2], 1)
    ) %>% sf::st_sfc(crs = EPSG_TARGET)) %>%
    sf::st_as_sf() %>%
    sf::st_filter(boundary_strip)

  pour_pts_sf <- endpoints %>%
    dplyr::select(geometry) %>%
    dplyr::slice_head(n = 10)

} else {

  ## ---- Spatially thin: one point per 5 km grid cell, keep top 10 by flow ----
  # 5 km thinning + hard cap prevents hundreds of small catchments
  thin_km     <- 5000   # 5 km in UTM metres
  max_pts     <- 15L    # hard cap for the final number of pour points

  pour_pts_sf <- pour_pts_raw %>%
    sf::st_as_sf() %>%
    dplyr::rename(flow_acc_val = 1) %>%
    dplyr::mutate(
      grid_x = floor(sf::st_coordinates(.)[, 1] / thin_km),
      grid_y = floor(sf::st_coordinates(.)[, 2] / thin_km)
    ) %>%
    dplyr::group_by(grid_x, grid_y) %>%
    dplyr::slice_max(flow_acc_val, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::slice_max(flow_acc_val, n = max_pts, with_ties = FALSE) %>%
    dplyr::select(-grid_x, -grid_y)

  message(sprintf("Identified %d pour points after 5 km thinning + top-%d filter",
                  nrow(pour_pts_sf), max_pts))
}

## ---- Validate: pour points must be inside (or on edge of) district ----
# Use a small buffer because raster cell centroids from terra::as.points() can
# land microscopically outside the vector polygon boundary due to grid precision.
pour_pts_inside <- sf::st_filter(pour_pts_sf, sf::st_buffer(mecufi_utm, 300))
if (nrow(pour_pts_inside) == 0) {
  warning("Pour points are outside the district boundary after detection — check geometry")
  pour_pts_inside <- pour_pts_sf   # proceed anyway
}
message(sprintf("  %d pour points retained after boundary filter (of %d candidates)",
                nrow(pour_pts_inside), nrow(pour_pts_sf)))

sf::st_write(pour_pts_inside, "gis/derived/mecufi_pour_points.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

# ______________________________________________________________________________
# SNAP POUR POINTS TO STREAM NETWORK ----
# ______________________________________________________________________________

message("03_hydrology: snapping pour points to stream network...")

# WBT snap_pour_points requires shapefiles (not gpkg)
pour_pts_shp  <- wbt_path("gis/derived/mecufi_pour_points_snapped.shp")
pour_pts_input <- wbt_path("gis/derived/mecufi_pour_points_raw.shp")

# Write unsnapped pour points as shapefile
sf::st_write(pour_pts_inside, gsub("\\.gpkg$", "_raw.shp",
             "gis/derived/mecufi_pour_points.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

pour_pts_input <- wbt_path("gis/derived/mecufi_pour_points_raw.shp")

whitebox::wbt_snap_pour_points(
  pour_pts  = pour_pts_input,
  flow_accum = flow_acc_path,
  output    = pour_pts_shp,
  snap_dist = POUR_SNAP_DIST
)

## ---- Validate snapped points inside district ----
snapped_pts <- sf::st_read("gis/derived/mecufi_pour_points_snapped.shp", quiet = TRUE)
snapped_inside <- sf::st_filter(snapped_pts,
                                sf::st_buffer(mecufi_utm, POUR_SNAP_DIST + 50))

if (nrow(snapped_inside) == 0) {
  warning(sprintf(
    "All snapped pour points outside district + %dm buffer — reducing snap distance",
    POUR_SNAP_DIST
  ))
}

message(sprintf("%d of %d pour points remain inside district after snapping",
                nrow(snapped_inside), nrow(snapped_pts)))

# ______________________________________________________________________________
# WATERSHED DELINEATION ----
# ______________________________________________________________________________

message("03_hydrology: delineating watersheds...")

watersheds_path <- wbt_path("gis/derived/mecufi_watersheds_raw.tif")

whitebox::wbt_watershed(
  d8_pntr   = flow_dir_path,
  pour_pts  = pour_pts_shp,
  output    = watersheds_path
)

## ---- Convert raster watersheds to polygons ----
watersheds_r <- terra::rast("gis/derived/mecufi_watersheds_raw.tif")

watersheds_utm <- watersheds_r %>%
  terra::as.polygons(dissolve = TRUE) %>%
  sf::st_as_sf() %>%
  sf::st_intersection(mecufi_utm) %>%
  dplyr::rename(watershed_id = 1) %>%
  dplyr::filter(!sf::st_is_empty(geometry))

assert_rows(watersheds_utm, 1, "watersheds_utm")

sf::st_write(watersheds_utm, "gis/derived/mecufi_watersheds.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

message(sprintf("03_hydrology complete: %d watershed polygons delineated",
                nrow(watersheds_utm)))

} # end cache-miss block

# ______________________________________________________________________________
# 01_load_data.R — Load and filter all spatial layers to Mecufi District
# ______________________________________________________________________________
# Purpose : Read admin boundaries, DEM tiles, waterways, and roads from their
#           compressed source archives. Filter all layers to Mecufi District,
#           reproject to UTM 37S, and cache to gis/derived/.
# Inputs  : gis/admin/moz_admin_boundaries.shp.zip
#           gis/raster/elevation/rasters_COP30.tar.gz
#           gis/raster/population/MOZ_population_v1_1_gridded.tif  (Maxar)
#           gis/raster/population/moz_pop_2026_CN_100m_R2025A_v1.tif (WorldPop)
#           gis/roads/moz_roads_shp.zip
#           gis/waterways/hotosm_moz_waterways_lines_shp.zip
#           gis/waterways/hotosm_moz_waterways_polygons_shp.zip
# Outputs : gis/derived/mecufi_boundary_utm.gpkg
#           gis/derived/mecufi_dem_raw.tif       (mosaiced COP30 tiles, WGS84)
#           gis/derived/mecufi_waterways_lines_utm.gpkg
#           gis/derived/mecufi_waterways_polys_utm.gpkg
#           gis/derived/mecufi_roads_utm.gpkg
#           gis/derived/mecufi_pop_real_utm.tif  (Maxar building-footprint pop.)
#           gis/derived/mecufi_pop_estimated_utm.tif (WorldPop 2026 RF 100 m)
# ______________________________________________________________________________

# ______________________________________________________________________________
# CACHE CHECK ----
# ______________________________________________________________________________

derived_files_01 <- c(
  "gis/derived/mecufi_boundary_utm.gpkg",
  "gis/derived/mecufi_dem_raw.tif",
  "gis/derived/mecufi_waterways_lines_utm.gpkg",
  "gis/derived/mecufi_waterways_polys_utm.gpkg",
  "gis/derived/mecufi_roads_utm.gpkg",
  "gis/derived/mecufi_places_utm.gpkg",
  "gis/derived/mecufi_pop_real_utm.tif",
  "gis/derived/mecufi_pop_estimated_utm.tif"
)

if (all(file.exists(derived_files_01)) && !force_rerun["load_data"]) {
  message("01_load_data: cache hit — loading from gis/derived/")

  mecufi_utm          <- sf::st_read("gis/derived/mecufi_boundary_utm.gpkg",    quiet = TRUE)
  dem_raw             <- terra::rast("gis/derived/mecufi_dem_raw.tif")
  waterways_lines_utm <- sf::st_read("gis/derived/mecufi_waterways_lines_utm.gpkg", quiet = TRUE)
  waterways_polys_utm <- sf::st_read("gis/derived/mecufi_waterways_polys_utm.gpkg", quiet = TRUE)
  roads_utm           <- sf::st_read("gis/derived/mecufi_roads_utm.gpkg",       quiet = TRUE)
  places_utm          <- sf::st_read("gis/derived/mecufi_places_utm.gpkg",      quiet = TRUE)
  pop_real_utm        <- terra::rast("gis/derived/mecufi_pop_real_utm.tif")
  pop_estimated_utm   <- terra::rast("gis/derived/mecufi_pop_estimated_utm.tif")

} else {

# ______________________________________________________________________________
# ADMIN BOUNDARY ----
# ______________________________________________________________________________

## ---- Read district-level (ADM2) shapefile from zip ----
# The zip contains admin0–admin3 layers; moz_admin2.shp = district level
admin_zip  <- normalizePath("gis/admin/moz_admin_boundaries.shp.zip")
zip_contents <- unzip(admin_zip, list = TRUE)

# Prefer ADM2; fall back to first .shp if not found
shp_candidates <- zip_contents$Name[grepl("\\.shp$", zip_contents$Name)]
shp_name <- if (any(grepl("admin2\\.shp$", shp_candidates))) {
  shp_candidates[grepl("admin2\\.shp$", shp_candidates)][1]
} else {
  shp_candidates[1]
}

admin_raw <- sf::st_read(
  paste0("/vsizip/", admin_zip, "/", shp_name),
  quiet = TRUE
)

## ---- Diagnostic: reveal field names and Mecufi rows ----
message("\n--- Admin boundary field names ---")
message(paste(names(admin_raw), collapse = ", "))

# Search all character columns for rows matching DISTRICT_NAME
# Use st_drop_geometry() before names() to avoid geometry column leaking in
adm_df    <- sf::st_drop_geometry(admin_raw)
char_cols <- names(adm_df)[sapply(adm_df, is.character)]
mecufi_matches <- adm_df %>%
  dplyr::filter(dplyr::if_any(
    dplyr::all_of(char_cols),
    ~ stringr::str_detect(., stringr::regex(DISTRICT_NAME, ignore_case = TRUE))
  ))

message(sprintf("\n--- Rows matching '%s' ---", DISTRICT_NAME))
print(mecufi_matches[, char_cols])

## ---- Identify the filter column (first column containing DISTRICT_NAME) ----
filter_col <- char_cols[sapply(char_cols, function(col) {
  any(stringr::str_detect(
    admin_raw[[col]],
    stringr::regex(DISTRICT_NAME, ignore_case = TRUE)
  ), na.rm = TRUE)
})][1]

if (is.na(filter_col)) {
  stop(sprintf(
    "Could not find '%s' in any character column of admin boundary.\nColumns searched: %s",
    DISTRICT_NAME, paste(char_cols, collapse = ", ")
  ))
}
message(sprintf("Filtering on column: %s", filter_col))

## ---- Filter to Mecufi and reproject ----
mecufi_utm <- admin_raw %>%
  dplyr::filter(stringr::str_detect(
    .data[[filter_col]],
    stringr::regex(DISTRICT_NAME, ignore_case = TRUE)
  )) %>%
  sf::st_transform(EPSG_TARGET)

assert_rows(mecufi_utm, 1, "mecufi_utm")
sf::st_write(mecufi_utm, "gis/derived/mecufi_boundary_utm.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

# ______________________________________________________________________________
# DEM (COP30) ----
# ______________________________________________________________________________

## ---- Extract tar.gz into dedicated subdirectory to avoid mixing with other rasters ----
raster_tar <- "gis/raster/elevation/rasters_COP30.tar.gz"
raster_dir  <- "gis/raster/elevation/cop30"

dir.create(raster_dir, showWarnings = FALSE, recursive = TRUE)

if (length(list.files(raster_dir, pattern = "\\.tif$", recursive = TRUE)) == 0) {
  message("Extracting COP30 raster archive...")
  untar(raster_tar, exdir = raster_dir)
}

## ---- Find all DEM tiles ----
dem_files <- list.files(
  raster_dir,
  pattern    = "\\.tif$",
  full.names = TRUE,
  recursive  = TRUE
)

if (length(dem_files) == 0) {
  stop("No .tif files found in gis/raster/ after extraction. Check archive contents.")
}
message(sprintf("Found %d DEM tile(s): %s", length(dem_files),
                paste(basename(dem_files), collapse = ", ")))

## ---- Mosaic multiple tiles (or use single tile directly) ----
if (length(dem_files) == 1) {
  dem_raw <- terra::rast(dem_files[1])
} else {
  tile_list <- lapply(dem_files, terra::rast)
  dem_raw   <- terra::mosaic(terra::sprc(tile_list), fun = "mean")
}

## ---- Crop to Mecufi with 5 km buffer (in WGS84 to match DEM CRS) ----
mecufi_wgs84   <- sf::st_transform(mecufi_utm, 4326)
mecufi_buf5km  <- sf::st_buffer(sf::st_transform(mecufi_utm, 4326),
                                dist = 0.05)  # ~5 km in degrees at this latitude
mecufi_buf_vect <- terra::vect(mecufi_buf5km)

dem_raw <- terra::crop(dem_raw, mecufi_buf_vect)

assert_overlap(dem_raw, mecufi_utm)
terra::writeRaster(dem_raw, "gis/derived/mecufi_dem_raw.tif",
                   overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

# ______________________________________________________________________________
# WATERWAYS ----
# ______________________________________________________________________________

## ---- Create 500 m buffer in UTM for spatial filtering ----
mecufi_buf500 <- sf::st_buffer(mecufi_utm, 500)

## ---- Waterway lines ----
ww_lines_zip  <- normalizePath("gis/waterways/hotosm_moz_waterways_lines_shp.zip")
ww_lines_contents <- unzip(ww_lines_zip, list = TRUE)
ww_lines_shp  <- ww_lines_contents$Name[grepl("\\.shp$", ww_lines_contents$Name)][1]

waterways_lines_utm <- sf::st_read(
  paste0("/vsizip/", ww_lines_zip, "/", ww_lines_shp),
  quiet = TRUE
) %>%
  sf::st_transform(EPSG_TARGET) %>%
  sf::st_filter(mecufi_buf500)

sf::st_write(waterways_lines_utm, "gis/derived/mecufi_waterways_lines_utm.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

## ---- Waterway polygons ----
ww_polys_zip  <- normalizePath("gis/waterways/hotosm_moz_waterways_polygons_shp.zip")
ww_polys_contents <- unzip(ww_polys_zip, list = TRUE)
ww_polys_shp  <- ww_polys_contents$Name[grepl("\\.shp$", ww_polys_contents$Name)][1]

waterways_polys_utm <- sf::st_read(
  paste0("/vsizip/", ww_polys_zip, "/", ww_polys_shp),
  quiet = TRUE
) %>%
  sf::st_transform(EPSG_TARGET) %>%
  sf::st_filter(mecufi_buf500)

sf::st_write(waterways_polys_utm, "gis/derived/mecufi_waterways_polys_utm.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

# ______________________________________________________________________________
# ROADS ----
# ______________________________________________________________________________

roads_zip  <- normalizePath("gis/roads/hotosm_moz_roads_lines_shp.zip")
roads_contents <- unzip(roads_zip, list = TRUE)
roads_shp  <- roads_contents$Name[grepl("[.]shp$", roads_contents$Name)][1]

## ---- 2 km buffer for roads (show context near district boundary) ----
mecufi_buf2km <- sf::st_buffer(mecufi_utm, 2000)

# Keep motorway → tertiary (named/classified roads) plus tracks and paths;
# exclude service, living_street, steps, construction, etc.
osm_road_types <- c(
  "motorway", "motorway_link",
  "trunk", "trunk_link",
  "primary", "primary_link",
  "secondary", "secondary_link",
  "tertiary", "tertiary_link",
  "unclassified", "residential",
  "track", "path", "footway", "bridleway"
)

roads_utm <- sf::st_read(
  paste0("/vsizip/", roads_zip, "/", roads_shp),
  quiet = TRUE
) %>%
  dplyr::filter(highway %in% osm_road_types) %>%
  sf::st_transform(EPSG_TARGET) %>%
  sf::st_filter(mecufi_buf2km)

sf::st_write(roads_utm, "gis/derived/mecufi_roads_utm.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

# ______________________________________________________________________________
# POPULATED PLACES ----
# ______________________________________________________________________________

places_zip  <- normalizePath("gis/places/hotosm_moz_populated_places_points_shp.zip")
places_shp  <- "hotosm_moz_populated_places_points_shp.shp"

places_utm <- sf::st_read(
  paste0("/vsizip/", places_zip, "/", places_shp),
  quiet = TRUE
) %>%
  dplyr::filter(place %in% c("city", "town", "village", "hamlet")) %>%
  sf::st_transform(EPSG_TARGET) %>%
  sf::st_filter(mecufi_utm)

sf::st_write(places_utm, "gis/derived/mecufi_places_utm.gpkg",
             delete_dsn = TRUE, quiet = TRUE)

# ______________________________________________________________________________
# POPULATION RASTERS ----
# ______________________________________________________________________________

## ---- Population raster (Maxar / MOZ gridded v1.1) ----
# Building-footprint-constrained population counts, ~93 m native resolution.
# Cropped to 500 m buffer, projected to UTM 37S, masked to district boundary.
pop_real_src <- "gis/raster/population/MOZ_population_v1_1_gridded.tif"

# Raster is WGS84 — buffer and crop in WGS84, then project to UTM
mecufi_buf_wgs84 <- sf::st_transform(sf::st_buffer(mecufi_utm, 500), 4326)

pop_real_utm <- terra::rast(pop_real_src) %>%
  terra::crop(terra::vect(mecufi_buf_wgs84)) %>%
  terra::project(paste0("EPSG:", EPSG_TARGET), method = "bilinear") %>%
  terra::mask(terra::vect(mecufi_utm))

assert_overlap(pop_real_utm, mecufi_utm)

terra::writeRaster(
  pop_real_utm,
  "gis/derived/mecufi_pop_real_utm.tif",
  overwrite = TRUE,
  gdal      = "COMPRESS=DEFLATE"
)

## ---- Population raster (WorldPop 2026 RF / 100 m) ----
# Random-forest modelled counts constrained to census totals, 100 m resolution.
# Continuous grid-wide estimates — no building mask.
pop_est_src <- "gis/raster/population/moz_pop_2026_CN_100m_R2025A_v1.tif"
mecufi_buf_wgs84_est <- sf::st_transform(sf::st_buffer(mecufi_utm, 500), 4326)

pop_estimated_utm <- terra::rast(pop_est_src) %>%
  terra::crop(terra::vect(mecufi_buf_wgs84_est)) %>%
  terra::project(paste0("EPSG:", EPSG_TARGET), method = "bilinear") %>%
  terra::mask(terra::vect(mecufi_utm))

assert_overlap(pop_estimated_utm, mecufi_utm)

terra::writeRaster(
  pop_estimated_utm,
  "gis/derived/mecufi_pop_estimated_utm.tif",
  overwrite = TRUE,
  gdal      = "COMPRESS=DEFLATE"
)

message(sprintf(
  "01_load_data complete: %d waterway lines, %d waterway polys, %d road features, %d places",
  nrow(waterways_lines_utm), nrow(waterways_polys_utm), nrow(roads_utm), nrow(places_utm)
))

} # end cache-miss block

# ______________________________________________________________________________
# 02_dem_processing.R — DEM reprojection, sink fill, and terrain derivatives
# ______________________________________________________________________________
# Purpose : Reproject the raw COP30 DEM to UTM 37S, crop and mask to Mecufi,
#           fill sinks, and compute hillshade, slope, and aspect for display
#           and hydrological analysis.
# Inputs  : gis/derived/mecufi_dem_raw.tif   (WGS84)
#           gis/derived/mecufi_boundary_utm.gpkg
# Outputs : gis/derived/mecufi_dem_utm.tif       (UTM, not filled)
#           gis/derived/mecufi_dem_filled.tif     (UTM, sink-filled)
#           gis/derived/mecufi_slope_deg.tif      (degrees)
#           gis/derived/mecufi_aspect_deg.tif     (degrees)
#           gis/derived/mecufi_hillshade.tif
# ______________________________________________________________________________

# ______________________________________________________________________________
# CACHE CHECK ----
# ______________________________________________________________________________

derived_files_02 <- c(
  "gis/derived/mecufi_dem_filled.tif",
  "gis/derived/mecufi_hillshade.tif",
  "gis/derived/mecufi_slope_deg.tif",
  "gis/derived/mecufi_aspect_deg.tif"
)

if (all(file.exists(derived_files_02)) && !force_rerun["dem_proc"]) {
  message("02_dem_processing: cache hit — loading from gis/derived/")

  dem_filled <- terra::rast("gis/derived/mecufi_dem_filled.tif")
  hillshade  <- terra::rast("gis/derived/mecufi_hillshade.tif")
  slope_deg  <- terra::rast("gis/derived/mecufi_slope_deg.tif")
  aspect_deg <- terra::rast("gis/derived/mecufi_aspect_deg.tif")

} else {

# ______________________________________________________________________________
# REPROJECT DEM ----
# ______________________________________________________________________________

## ---- Load inputs ----
dem_raw    <- terra::rast("gis/derived/mecufi_dem_raw.tif")
mecufi_utm <- sf::st_read("gis/derived/mecufi_boundary_utm.gpkg", quiet = TRUE)

## ---- Reproject to UTM 37S ----
message("02_dem_processing: reprojecting DEM to UTM 37S...")
dem_utm <- terra::project(dem_raw,
                          y      = paste0("EPSG:", EPSG_TARGET),
                          method = "bilinear")

# Crop and mask to Mecufi boundary
mecufi_vect <- terra::vect(mecufi_utm)
dem_utm     <- terra::crop(dem_utm, mecufi_vect)
dem_utm     <- terra::mask(dem_utm, mecufi_vect)

terra::writeRaster(dem_utm, "gis/derived/mecufi_dem_utm.tif",
                   overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

## ---- Report elevation range (use to tune CONTOUR_INT and hypsometric breaks) ----
elev_vals <- terra::values(dem_utm, na.rm = TRUE)
message(sprintf(
  "Elevation range: %.0f – %.0f m  |  mean: %.0f m  |  SD: %.0f m",
  min(elev_vals), max(elev_vals), mean(elev_vals), sd(elev_vals)
))

# ______________________________________________________________________________
# SINK FILLING ----
# ______________________________________________________________________________

message("02_dem_processing: filling sinks with WhiteboxTools...")

dem_utm_path    <- wbt_path("gis/derived/mecufi_dem_utm.tif")
filled_path     <- wbt_path("gis/derived/mecufi_dem_filled.tif")

whitebox::wbt_fill_depressions(
  dem  = dem_utm_path,
  output = filled_path
)

dem_filled <- terra::rast("gis/derived/mecufi_dem_filled.tif")

# ______________________________________________________________________________
# SLOPE AND ASPECT ----
# ______________________________________________________________________________

message("02_dem_processing: computing slope and aspect...")

slope_path  <- wbt_path("gis/derived/mecufi_slope_deg.tif")
aspect_path <- wbt_path("gis/derived/mecufi_aspect_deg.tif")

whitebox::wbt_slope(
  dem    = filled_path,
  output = slope_path,
  units  = "degrees"
)

whitebox::wbt_aspect(
  dem    = filled_path,
  output = aspect_path
)

slope_deg  <- terra::rast("gis/derived/mecufi_slope_deg.tif")
aspect_deg <- terra::rast("gis/derived/mecufi_aspect_deg.tif")

# ______________________________________________________________________________
# HILLSHADE ----
# ______________________________________________________________________________
# NOTE: terra::shade() requires radians. WBT slope/aspect are in degrees.
# Convert: multiply by pi/180.

message("02_dem_processing: computing hillshade...")

hillshade <- terra::shade(
  slope     = slope_deg * pi / 180,
  aspect    = aspect_deg * pi / 180,
  angle     = 45,
  direction = 315
)

terra::writeRaster(hillshade, "gis/derived/mecufi_hillshade.tif",
                   overwrite = TRUE, gdal = "COMPRESS=DEFLATE")

assert_crs(dem_filled, EPSG_TARGET)
assert_overlap(dem_filled, mecufi_utm)

message("02_dem_processing complete.")

} # end cache-miss block

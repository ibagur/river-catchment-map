# ______________________________________________________________________________
# 04_contours.R — Contour generation, smoothing, and classification
# ______________________________________________________________________________
# Purpose : Generate elevation contours from the filled DEM at CONTOUR_INT
#           intervals, smooth them to remove pixel jaggedness, classify into
#           major/minor lines, and clip to the district boundary.
# Inputs  : gis/derived/mecufi_dem_filled.tif
#           gis/derived/mecufi_boundary_utm.gpkg
# Outputs : gis/derived/mecufi_contours.gpkg
# ______________________________________________________________________________

# ______________________________________________________________________________
# CACHE CHECK ----
# ______________________________________________________________________________

derived_file_04 <- "gis/derived/mecufi_contours.gpkg"

if (file.exists(derived_file_04) && !force_rerun["contours"]) {
  message("04_contours: cache hit — loading from gis/derived/")

  contours_utm <- sf::st_read(derived_file_04, quiet = TRUE)

} else {

# ______________________________________________________________________________
# LOAD INPUTS ----
# ______________________________________________________________________________

dem_filled <- terra::rast("gis/derived/mecufi_dem_filled.tif")
mecufi_utm <- sf::st_read("gis/derived/mecufi_boundary_utm.gpkg", quiet = TRUE)

# ______________________________________________________________________________
# CONTOUR GENERATION ----
# ______________________________________________________________________________

## ---- Compute actual elevation range for contour levels ----
elev_vals <- terra::values(dem_filled, na.rm = TRUE)
elev_min  <- floor(min(elev_vals) / CONTOUR_INT) * CONTOUR_INT
elev_max  <- ceiling(max(elev_vals) / CONTOUR_INT) * CONTOUR_INT

contour_levels <- seq(elev_min, elev_max, by = CONTOUR_INT)
message(sprintf("04_contours: generating %d contour levels (%d – %d m at %d m intervals)",
                length(contour_levels), elev_min, elev_max, CONTOUR_INT))

## ---- Generate contours ----
contours_raw <- terra::as.contour(dem_filled, levels = contour_levels) %>%
  sf::st_as_sf()

assert_rows(contours_raw, 1, "contours_raw")

# ______________________________________________________________________________
# CLASSIFY MAJOR / MINOR ----
# ______________________________________________________________________________
# Major contours: every 2nd interval (i.e., level divisible by CONTOUR_INT * 2)

contours_classified <- contours_raw %>%
  dplyr::rename(level = 1) %>%
  dplyr::mutate(
    level    = as.numeric(level),
    is_major = (level %% (CONTOUR_INT * 2L)) == 0
  )

# ______________________________________________________________________________
# SMOOTH AND CLIP ----
# ______________________________________________________________________________

message("04_contours: smoothing contours...")

contours_smooth <- contours_classified %>%
  smoothr::smooth(method = "ksmooth", smoothness = 2)

## ---- Clip to district boundary ----
contours_utm <- contours_smooth %>%
  sf::st_intersection(mecufi_utm) %>%
  dplyr::filter(!sf::st_is_empty(geometry))

assert_rows(contours_utm, 1, "contours_utm")

sf::st_write(contours_utm, derived_file_04, delete_dsn = TRUE, quiet = TRUE)

message(sprintf("04_contours complete: %d contour lines (%d major, %d minor)",
                nrow(contours_utm),
                sum(contours_utm$is_major),
                sum(!contours_utm$is_major)))

} # end cache-miss block

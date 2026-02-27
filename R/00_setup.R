# ______________________________________________________________________________
# 00_setup.R — Packages, WhiteboxTools, directories, and helper functions
# ______________________________________________________________________________
# Purpose : Load all required packages, initialise WhiteboxTools, create output
#           directories, and define shared helper functions used throughout the
#           pipeline.
# Called  : source("R/00_setup.R") from run_analysis.R
# ______________________________________________________________________________

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(whitebox)
  library(tidyterra)
  library(tmap)
  library(smoothr)
  library(dplyr)
  library(stringr)
  library(scales)
})

# ______________________________________________________________________________
# GUARDS ----
# ______________________________________________________________________________

## ---- tmap version guard ----
stopifnot(
  "tmap >= 4 is required — run: install.packages('tmap')" =
    as.integer(substr(packageVersion("tmap"), 1, 1)) >= 4L
)

## ---- WhiteboxTools initialisation ----
if (!check_whitebox_binary()) {
  message("WhiteboxTools not found — installing...")
  install_whitebox()
}
wbt_init()
wbt_verbose(FALSE)   # suppress per-tool console chatter

# ______________________________________________________________________________
# OUTPUT DIRECTORIES ----
# ______________________________________________________________________________

dir.create("gis/derived", showWarnings = FALSE, recursive = TRUE)
dir.create("output",      showWarnings = FALSE, recursive = TRUE)

# ______________________________________________________________________________
# HELPER FUNCTIONS ----
# ______________________________________________________________________________

## ---- Helper: assert_crs ----
# Errors if the CRS of x does not match the expected EPSG code.
assert_crs <- function(x, epsg) {
  actual <- if (inherits(x, "SpatRaster") || inherits(x, "SpatVector")) {
    as.integer(terra::crs(x, describe = TRUE)$code)
  } else {
    sf::st_crs(x)$epsg
  }
  if (!isTRUE(actual == as.integer(epsg))) {
    stop(sprintf(
      "CRS mismatch: expected EPSG:%d, got EPSG:%s\n  Object class: %s",
      as.integer(epsg), actual, class(x)[1]
    ))
  }
  invisible(x)
}

## ---- Helper: assert_overlap ----
# Errors if the extents of a SpatRaster and an sf/SpatVector do not overlap.
# Reprojects the vector to the raster CRS when they differ.
assert_overlap <- function(raster, vector) {
  r_crs <- terra::crs(raster)
  r_ext <- as.vector(terra::ext(raster))

  # Reproject vector to raster CRS for comparison
  if (inherits(vector, "sf")) {
    v_crs <- sf::st_crs(vector)
    r_crs_sf <- sf::st_crs(raster)
    if (!isTRUE(v_crs == r_crs_sf)) {
      vector <- sf::st_transform(vector, r_crs_sf)
    }
    v_bb  <- as.vector(sf::st_bbox(vector))
    v_ext <- c(v_bb[1], v_bb[3], v_bb[2], v_bb[4])
  } else {
    vect_reproj <- if (!terra::same.crs(raster, vector)) {
      terra::project(vector, raster)
    } else {
      vector
    }
    v_ext <- as.vector(terra::ext(vect_reproj))
  }

  x_overlap <- r_ext[1] < v_ext[2] && r_ext[2] > v_ext[1]
  y_overlap <- r_ext[3] < v_ext[4] && r_ext[4] > v_ext[3]
  if (!x_overlap || !y_overlap) {
    stop(sprintf(
      "No spatial overlap between raster and vector.\n  Raster ext: %s\n  Vector ext (reprojected): %s",
      paste(round(r_ext, 4), collapse = ", "),
      paste(round(v_ext, 4), collapse = ", ")
    ))
  }
  invisible(raster)
}

## ---- Helper: assert_rows ----
# Errors if x has fewer than n_min rows.
assert_rows <- function(x, n_min = 1L, label = deparse(substitute(x))) {
  n <- if (inherits(x, "sf")) nrow(x) else nrow(as.data.frame(x))
  if (n < n_min) {
    stop(sprintf("Expected >= %d rows in '%s', got %d", n_min, label, n))
  }
  invisible(x)
}

## ---- Helper: wbt_path ----
# Wraps normalizePath() with mustWork = FALSE.
# WhiteboxTools REQUIRES absolute paths — always use this wrapper.
wbt_path <- function(...) {
  normalizePath(file.path(...), mustWork = FALSE)
}

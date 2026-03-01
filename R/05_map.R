# ______________________________________________________________________________
# 05_map.R — tmap v4 cartographic composition and export
# ______________________________________________________________________________
# Purpose : Compose publication-quality maps of Mecufi District river
#           catchments, blending hillshade, hypsometric DEM, flow accumulation,
#           contours, catchment polygons, waterways, roads, and map furniture.
#           Produces two variants differing only in the population exposure layer:
#             _real      — Maxar building-footprint counts (sparse, validated)
#             _estimated — WorldPop 2026 RF model (continuous, grid-wide)
# Inputs  : gis/derived/mecufi_dem_filled.tif
#           gis/derived/mecufi_hillshade.tif
#           gis/derived/mecufi_flow_acc.tif
#           gis/derived/mecufi_contours.gpkg
#           gis/derived/mecufi_watersheds.gpkg
#           gis/derived/mecufi_waterways_lines_utm.gpkg
#           gis/derived/mecufi_waterways_polys_utm.gpkg
#           gis/derived/mecufi_roads_utm.gpkg
#           gis/derived/mecufi_boundary_utm.gpkg
#           gis/derived/mecufi_pop_real_utm.tif
#           gis/derived/mecufi_pop_estimated_utm.tif
# Outputs : output/mecufi_catchment_map_real.png      (300 dpi)
#           output/mecufi_catchment_map_estimated.png  (300 dpi)
#           output/mecufi_contours.shp
#           output/mecufi_watersheds.shp
#           output/mecufi_flow_acc.tif
# ______________________________________________________________________________

if (!force_rerun["map"]) {
  if (all(file.exists(c("output/mecufi_catchment_map_real.png",
                        "output/mecufi_catchment_map_estimated.png")))) {
    message("05_map: cache hit — both maps exist (set force_rerun['map'] = TRUE to regenerate)")
    stop("Maps exist and force_rerun['map'] is FALSE — stopping gracefully.", call. = FALSE)
  }
}

# ______________________________________________________________________________
# LOAD ALL LAYERS ----
# ______________________________________________________________________________

message("05_map: loading layers...")

dem_filled          <- terra::rast("gis/derived/mecufi_dem_filled.tif")
hillshade_r         <- terra::rast("gis/derived/mecufi_hillshade.tif")
flow_acc            <- terra::rast("gis/derived/mecufi_flow_acc.tif")
contours_utm        <- sf::st_read("gis/derived/mecufi_contours.gpkg",           quiet = TRUE)
watersheds_utm      <- sf::st_read("gis/derived/mecufi_watersheds.gpkg",          quiet = TRUE)
waterways_lines_utm <- sf::st_read("gis/derived/mecufi_waterways_lines_utm.gpkg", quiet = TRUE)
waterways_polys_utm <- sf::st_read("gis/derived/mecufi_waterways_polys_utm.gpkg", quiet = TRUE)
roads_utm           <- sf::st_read("gis/derived/mecufi_roads_utm.gpkg",           quiet = TRUE)
places_utm          <- sf::st_read("gis/derived/mecufi_places_utm.gpkg",          quiet = TRUE)
mecufi_utm          <- sf::st_read("gis/derived/mecufi_boundary_utm.gpkg",        quiet = TRUE)
pop_real_utm        <- terra::rast("gis/derived/mecufi_pop_real_utm.tif")
pop_estimated_utm   <- terra::rast("gis/derived/mecufi_pop_estimated_utm.tif")

# ______________________________________________________________________________
# PREPARE RASTERS ----
# ______________________________________________________________________________

## ---- River glow: distance-decay halo centred on channel pixels ----
# Technique:
#   1. Threshold flow_acc to get binary channel mask (NA = off-channel)
#   2. terra::distance() gives metres to nearest channel for every cell
#   3. Exponential decay: e^(-d / r) → 1 at channel, fades to ~0 at radius r
#   4. Modulate intensity by log(flow_acc) so main stems glow brighter/wider
#      than tributaries even at the same spatial distance
flow_thresh_display <- 1500L   # ~1.4 km² upstream area → medium+ streams only
decay_radius_m      <- 250     # e-folding glow width in UTM metres (~8 pixels at 30 m)

# Step 1: channel mask — non-NA pixels are the river network
channels_mask <- terra::ifel(flow_acc >= flow_thresh_display, 1L, NA)

# Step 2: distance to nearest channel pixel (metres; UTM so units are correct)
message("05_map: computing distance-to-channel for river glow...")
dist_to_channel <- terra::distance(channels_mask)

# Step 3: exponential spatial decay
glow_spatial <- exp(-dist_to_channel / decay_radius_m)

# Step 4: flow-magnitude weight (log-normalised to [0,1])
flow_log  <- log1p(terra::subst(flow_acc, NA, 0))
flow_norm <- flow_log / max(terra::values(flow_log, na.rm = TRUE), na.rm = TRUE)

# Blend: 30 % flat glow + 70 % flow-weighted → tributaries show but main
# stems are distinctly stronger
river_glow <- glow_spatial * (0.30 + 0.70 * flow_norm)

# Clip near-zero values to NA so off-channel terrain stays fully transparent
river_glow_display <- terra::ifel(river_glow < 0.05, NA, river_glow)

## ---- Flood proximity: shared base for both exposure layers ----
# Lower threshold (500 cells) than display (1500) to capture more of the
# network — representing the broader floodplain zone rather than visual halo.
flood_channels <- terra::ifel(flow_acc >= 500L, 1L, NA)

message("05_map: computing distance-to-flood-channel for exposure layers...")
dist_to_flood  <- terra::distance(flood_channels)
flood_prox     <- exp(-dist_to_flood / 400)   # 400 m e-folding decay

## ---- Exposure: Real (Maxar) — sparse, building-masked ----
pop_real_resampled <- terra::resample(pop_real_utm, flow_acc, method = "bilinear")
pop_real_resampled <- terra::subst(pop_real_resampled, NA, 0)
exposure_raw_real  <- pop_real_resampled * flood_prox
exposure_log_real  <- log1p(exposure_raw_real)
exposure_norm_real <- exposure_log_real / terra::global(exposure_log_real, "max", na.rm = TRUE)[[1]]

# Mask to cells where Maxar confirms buildings (nearest-neighbour preserves binary signal)
pop_mask          <- terra::resample(pop_real_utm, flow_acc, method = "near")
exposure_disp_real <- terra::ifel(is.na(pop_mask) | pop_mask <= 0, NA, exposure_norm_real)

# Drop noisiest bottom decile to reduce speckle
low_cut_real       <- terra::global(exposure_disp_real,
                                    function(x) quantile(x, 0.10, na.rm = TRUE))[[1]]
exposure_disp_real <- terra::ifel(exposure_disp_real < low_cut_real, NA, exposure_disp_real)

## ---- Exposure: Estimated (WorldPop) — population-masked, no building constraint ----
pop_est_resampled   <- terra::resample(pop_estimated_utm, flow_acc, method = "bilinear")
pop_est_resampled   <- terra::subst(pop_est_resampled, NA, 0)
exposure_raw_est    <- pop_est_resampled * flood_prox
exposure_log_est    <- log1p(exposure_raw_est)
exposure_norm_est   <- exposure_log_est / terra::global(exposure_log_est, "max", na.rm = TRUE)[[1]]

# WorldPop assigns near-zero values to every land cell, creating a low yellow
# wash across the whole map. Two-stage filter:
#   1. Mask out cells where WorldPop population < 1 person per ~30 m cell
#      (after bilinear resampling from 100 m) — removes background noise
#   2. Drop bottom quartile of surviving values to sharpen the visual contrast
pop_est_mask        <- terra::ifel(pop_est_resampled >= 1, 1, NA)
exposure_disp_est   <- terra::ifel(is.na(pop_est_mask), NA, exposure_norm_est)
low_cut_est         <- terra::global(exposure_disp_est,
                                     function(x) quantile(x, 0.25, na.rm = TRUE))[[1]]
exposure_disp_est   <- terra::ifel(exposure_disp_est < low_cut_est, NA, exposure_disp_est)

## ---- Hypsometric colour intervals (25 m bands) ----
# 15 discrete bands give ~25 m resolution across the 0-366 m range.
# Breaks and colours are hardcoded to actual data range (0-366 m).
# Adjust breaks if running on a different DEM extent.
elev_vals  <- terra::values(dem_filled, na.rm = TRUE)
elev_min_r <- min(elev_vals)
elev_max_r <- max(elev_vals)

hyps_breaks <- c(0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)

# 15 colours: dark brown at sea level → warm tan lowlands →
# pale yellow mid-slopes → yellow-green → progressively darker greens highlands
hyps_colors <- c(
  "#6B3A1F",  #   0– 25 m  deep brown (coastal floodplain)
  "#8B5530",  #  25– 50 m  dark reddish-brown
  "#A67245",  #  50– 75 m  warm brown
  "#BE9060",  #  75–100 m  orange-brown
  "#CFA978",  # 100–125 m  tan
  "#DCC090",  # 125–150 m  sandy tan
  "#E8D4A8",  # 150–175 m  pale sand
  "#EDE0BB",  # 175–200 m  very pale yellow
  "#DCE4A0",  # 200–225 m  light yellow-green
  "#C4D880",  # 225–250 m  yellow-green
  "#A8C860",  # 250–275 m  medium yellow-green
  "#88B848",  # 275–300 m  medium green
  "#6AA038",  # 300–325 m  grass green
  "#508830",  # 325–350 m  darker green
  "#386820"   # 350–375 m  deep forest green
)

## ---- Split contours into major/minor subsets ----
contours_minor <- contours_utm %>% dplyr::filter(!is_major)
contours_major <- contours_utm %>% dplyr::filter( is_major)

# ______________________________________________________________________________
# MAP COMPOSITION ----
# ______________________________________________________________________________
# build_map() assembles the full layer stack, parameterised only by the
# exposure raster, its legend label, and the population credits suffix.
# All other objects are read from the shared pipeline environment.

## ---- Helper: assemble full tmap layer stack ----
build_map <- function(exposure_disp, pop_label, credits_suffix) {

  m <-
    ## Layer 1: Hillshade base
    tm_shape(hillshade_r) +
      tm_raster(
        col.scale  = tm_scale_continuous(values = "greys"),
        col_alpha  = 1.0,
        col.legend = tm_legend_hide()
      ) +

    ## Layer 2: Hypsometric DEM — 25 m discrete intervals
    tm_shape(dem_filled) +
      tm_raster(
        col.scale  = tm_scale_intervals(
          breaks = hyps_breaks,
          values = hyps_colors,
          labels = paste0(head(hyps_breaks, -1), "–", tail(hyps_breaks, -1), " m")
        ),
        col_alpha  = 0.72,
        col.legend = tm_legend(
          title    = "Elevation (m)",
          position = tm_pos_out("right", "center")
        )
      ) +

    ## Layer 3: River glow — distance-decay halo, light cyan → deep navy
    # Values near 1 = on-channel (dark); values near 0.015 = glow edge (light).
    tm_shape(river_glow_display) +
      tm_raster(
        col.scale  = tm_scale_continuous(
          values = c("#C8E6F5", "#56B4E9", "#0077B6", "#023E8A", "#03045E"),
          limits = c(0, 1)
        ),
        col_alpha  = 0.80,
        col.legend = tm_legend_hide()
      ) +

    ## Layer 4: Minor contours
    tm_shape(contours_minor) +
      tm_lines(
        col       = "#B0B0B0",
        lwd       = 0.2,
        col_alpha = 0.45
      ) +

    ## Layer 5: Major contours
    tm_shape(contours_major) +
      tm_lines(
        col       = "#888888",
        lwd       = 0.55,
        col_alpha = 0.65
      ) +

    ## Layer 6: Watershed polygons
    tm_shape(watersheds_utm) +
      tm_polygons(
        fill        = tm_const(),
        fill.scale  = tm_scale(values = "#2C5F8A"),
        fill_alpha  = 0.20,
        col         = "#2C5F8A",
        lwd         = 1.2,
        col.legend  = tm_legend_hide(),
        fill.legend = tm_legend_hide()
      ) +

    ## Layer 7: Flood exposure (population × drainage proximity)
    # Semi-transparent YlOrRd raster: yellow = low exposure, red = high exposure.
    tm_shape(exposure_disp) +
      tm_raster(
        col.scale  = tm_scale_intervals(
          breaks   = c(0, 0.33, 0.60, 0.80, 1.0),
          values   = "brewer.yl_or_rd",
          labels   = c("Low", "Medium", "High", "Very high"),
          value.na = NA
        ),
        col_alpha  = 0.70,
        col.legend = tm_legend(
          title       = pop_label,
          position    = tm_pos_in("left", "bottom"),
          orientation = "landscape",
          text.size   = 0.40
        )
      ) +

    ## Map furniture
    tm_graticules(
      col         = "grey70",
      alpha       = 0.4,
      lwd         = 0.3,
      labels.size = 0.45
    ) +
    tm_compass(
      type     = "8star",
      position = c("right", "top"),
      size     = 2.0
    ) +
    tm_scalebar(
      breaks    = c(0, 10, 20),
      position  = tm_pos_in("left", "bottom"),
      text.size = 0.5
    ) +
    tm_credits(
      text     = paste0(
        "Elevation: Copernicus GLO-30 DEM | Admin: OCHA | ",
        "Waterways/Roads: HOT OSM | Projection: UTM Zone 37S (EPSG:32737) | ",
        credits_suffix
      ),
      position = tm_pos_out("center", "bottom"),
      size     = 0.45
    ) +
    tm_title(
      text     = "River Catchments — Mecufi District, Cabo Delgado",
      size     = 1.0,
      fontface = "bold"
    ) +
    tm_layout(
      bg.color      = "#D6EAF8",
      outer.margins = 0.02,
      inner.margins = 0.04,
      frame         = TRUE,
      frame.lwd     = 1.5
    )

  ## Conditionally append waterway polygon layer
  if (nrow(waterways_polys_utm) > 0) {
    m <- m +
      tm_shape(waterways_polys_utm) +
        tm_polygons(
          fill        = tm_const(),
          fill.scale  = tm_scale(values = "#5BA3CE"),
          fill_alpha  = 0.80,
          col         = "#3A7BA8",
          lwd         = 0.5,
          col.legend  = tm_legend_hide(),
          fill.legend = tm_legend_hide()
        )
  }

  ## Conditionally append waterway line layer
  if (nrow(waterways_lines_utm) > 0) {
    m <- m +
      tm_shape(waterways_lines_utm) +
        tm_lines(
          col       = "#3A7BA8",
          lwd       = 0.8,
          col_alpha = 1.0
        )
  }

  ## Conditionally append roads layer
  if (nrow(roads_utm) > 0) {
    m <- m +
      tm_shape(roads_utm) +
        tm_lines(
          col        = "#CC2200",
          lwd        = 0.8,
          col_alpha  = 0.70,
          col.legend = tm_legend_hide()
        )
  }

  ## Conditionally append populated places layer
  if (nrow(places_utm) > 0) {
    m <- m +
      tm_shape(places_utm) +
        tm_dots(
          size      = 0.25,
          fill      = "#CC2200",
          col       = "#CC2200",
          lwd       = 0.8
        ) +
      tm_shape(places_utm) +
        tm_text(
          text      = "name",
          size      = 0.45,
          col       = "#444444",
          fontface  = "bold",
          options   = opt_tm_text(shadow = FALSE),
          xmod      = 0.4,
          ymod      = 0.4
        )
  }

  ## District boundary — always on top of optional layers
  m <- m +
    tm_shape(mecufi_utm) +
      tm_borders(
        col = "#444444",
        lwd = 2.0
      )

  m
}

message("05_map: composing maps...")
tmap_mode("plot")

map_real <- build_map(
  exposure_disp  = exposure_disp_real,
  pop_label      = "Flood exposure (Maxar)",
  credits_suffix = "Population: Maxar v1.1"
)

map_estimated <- build_map(
  exposure_disp  = exposure_disp_est,
  pop_label      = "Flood exposure (WorldPop 2026)",
  credits_suffix = "Population: WorldPop 2026 RF 100m"
)

# ______________________________________________________________________________
# EXPORT ----
# ______________________________________________________________________________

message("05_map: exporting output/mecufi_catchment_map_real.png ...")
tmap_save(map_real,
          filename = "output/mecufi_catchment_map_real.png",
          width    = 250,
          height   = 200,
          units    = "mm",
          dpi      = 300)

message("05_map: exporting output/mecufi_catchment_map_estimated.png ...")
tmap_save(map_estimated,
          filename = "output/mecufi_catchment_map_estimated.png",
          width    = 250,
          height   = 200,
          units    = "mm",
          dpi      = 300)

## ---- Export auxiliary GIS outputs ----
sf::st_write(contours_utm,
             "output/mecufi_contours.shp",
             delete_dsn = TRUE, quiet = TRUE)

sf::st_write(watersheds_utm,
             "output/mecufi_watersheds.shp",
             delete_dsn = TRUE, quiet = TRUE)

terra::writeRaster(flow_acc,
                   "output/mecufi_flow_acc.tif",
                   overwrite = TRUE,
                   gdal      = "COMPRESS=DEFLATE")

message("05_map complete.")
message("Outputs: output/mecufi_catchment_map_real.png + output/mecufi_catchment_map_estimated.png (300 dpi each)")

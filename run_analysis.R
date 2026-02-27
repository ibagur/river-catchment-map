# ______________________________________________________________________________
# River Catchment Map — Mecufi District, Cabo Delgado, Mozambique
# ______________________________________________________________________________
# Purpose : Orchestrate the full river catchment map pipeline.
#           Source each module in order; intermediate outputs are cached in
#           gis/derived/ so individual steps can be rerun without repeating
#           the full pipeline.
# Output  : output/mecufi_catchment_map.png (300 dpi)
# Usage   : source("run_analysis.R")  or  Rscript run_analysis.R
# ______________________________________________________________________________

# ______________________________________________________________________________
# GLOBAL PARAMETERS ----
# ______________________________________________________________________________

EPSG_TARGET    <- 32737L     # UTM Zone 37S
CONTOUR_INT    <- 50L        # Contour interval (metres); tune after Step 4 prints elev range
FLOW_THRESH    <- 5000L      # Flow accumulation cell count threshold for pour points
POUR_SNAP_DIST <- 100L       # Pour point snap distance (UTM metres)
DISTRICT_NAME  <- "Mecufi"

# ______________________________________________________________________________
# FORCE-RERUN FLAGS ----
# ______________________________________________________________________________
# Set TRUE to regenerate a cached step (and all downstream steps).
# Hydrology and DEM processing are the slowest; leave FALSE for styling runs.

force_rerun <- c(
  load_data = FALSE,
  dem_proc  = FALSE,
  hydrology = FALSE,
  contours  = FALSE,
  map       = TRUE           # Always regenerate map for styling iteration
)

# ______________________________________________________________________________
# PIPELINE ----
# ______________________________________________________________________________

source("R/00_setup.R")
source("R/01_load_data.R")
source("R/02_dem_processing.R")
source("R/03_hydrology.R")
source("R/04_contours.R")
source("R/05_map.R")

message("Pipeline complete. Check output/mecufi_catchment_map.png")

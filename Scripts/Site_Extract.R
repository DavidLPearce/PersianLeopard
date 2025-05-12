
# ------------------------------------------------------------------------------
#
#                         Load/Install Packages
#
# ------------------------------------------------------------------------------

library(terra)
library(sp)
library(sf)
library(raster)
library(spatialEco)
library(progress)

set.seed(123)
options(scipen = 9999)
setwd(".")

# ------------------------------------------------------------------------------
#
#                         Read in data
#
# ------------------------------------------------------------------------------


# Camera sites
sites <- vect("D:/Persian_Leopard_GIS/CameraLocationsNew/CameraLocationsNew.shp")

# Elevation
DEM <- rast("D:/Persian_Leopard_GIS/DEM20kmBuffer/DEM20kmbuffer.tif")

# Rivers
rivers <- vect("D:/Persian_Leopard_GIS/Rivers/riverclip.shp")

# Pastures
pasture <- vect("D:/Persian_Leopard_GIS/PastureLocations/pastureloc.shp")

# Villages
villages <- vect("D:/Persian_Leopard_GIS/VillageLocations/vilage.shp")

# ------------------------------------------------------------------------------
#
#                         Data Wrangling
#
# ------------------------------------------------------------------------------


# Check the CRS of each layer
print(sites)
print(DEM)
print(rivers)
print(pasture)
print(villages)

# Visually inspect
plot(sites)
plot(DEM)
plot(rivers)
plot(pasture)
plot(villages)

# Vector ruggedness measure to raster
VRM <- spatialEco::vrm(dem, s = 3) # Extract VRM # 3x3 window
plot(VRM)
VRM_rast <- VRM > 0.05 & VRM <= 0.20# Classify VRM values into raster
plot(VRM_rast)

# Dataframe for covariates
cov_df <- as.data.frame(sites)
cov_df <- cov_df[,-c(4,5)]
colnames(cov_df)[1] <- "SiteID"
print(cov_df)

# ------------------------------------------------------------------------------
#
#                         Extract Landscape Metrics
#
# ------------------------------------------------------------------------------

# Sequence of buffers
buffer_vec <- seq(from = 100, to = 15000, by = 100)

# Initialize the progress bar
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent | Site: :site | Elapsed: :elapsed",
  total = NROW(sites),
  clear = FALSE,
  width = 100
)

# row = 1
# buffer = 1

# Loop
for (row in 1:NROW(sites)) {

  # Subset site
  site_sub <- sites[row]
  
  # ---------------------
  # Distance Metrics
  # ---------------------
  
  # Nearest Pasture Distance
  past_dist <- distance(site_sub, pasture)
  nearest_past_dist <- apply(past_dist, 1, min)
  nearest_past_dist_col <- paste0("nrstPastureDist")
  cov_df[row, nearest_past_dist_col] <- nearest_past_dist
  
  # Nearest Village Distance
  vil_dist <- distance(site_sub, villages)
  nearest_vil_dist <- apply(vil_dist, 1, min)
  nearest_vil_dist_col <- paste0("nrstVillageDist")
  cov_df[row, nearest_vil_dist_col] <- nearest_vil_dist
  
  # Nearest River
  riv_dist <- distance(site_sub, rivers)
  nearest_riv_dist <- apply(riv_dist, 1, min)
  nearest_riv_dist_col <- paste0("nrstRiverDist")
  cov_df[row, nearest_riv_dist_col] <- nearest_riv_dist
  
  # ---------------------
  # Buffer Metrics 
  # ---------------------
  
  # Progress bar for buffers (re-initialize each site)
  buffer_pb <- progress_bar$new(
    format = paste0(sites$Field1[row], 
                    " | Buffer [:bar] :percent | Elapsed: :elapsed"),
    total = length(buffer_vec),
    clear = FALSE,
    width = 80
  )
  
  for (buffer in seq_along(buffer_vec)){
    
    # buffer size
    buffer_size <- buffer_vec[buffer]
  
    # Create a buffer around the site
    site_buffer <- terra::buffer(site_sub, width = buffer_size)
  
    # Subset DEM
    DEM_subset <- terra::crop(DEM, site_buffer)
    # plot(DEM_subset)
    
    # Subset VRM raster
    VRM_subset <- terra::crop(VRM_rast, site_buffer)
    # plot(VRM_subset)
    
    # Extract Elevation
    elev_values <- terra::extract(DEM_subset, site_buffer)
    elev_mn <- mean(elev_values[, 2])
    elev_mn_col <- paste0("elev_mn_", buffer_size, "m")
    cov_df[row, elev_mn_col] <- elev_mn

    # Extract Slope
    slope_values <- terra::terrain(x = DEM_subset, v = "slope", unit = "degrees")
    # plot(slope_values)
    slope_df <- as.data.frame(slope_values)
    slope_mn <- mean(slope_df[,1])
    slope_mn_col <- paste0("slope_mn_", buffer_size, "m")
    cov_df[row, slope_mn_col] <- slope_mn

    # Extract Aspect
    aspect_values <- terra::terrain(x = DEM_subset, v = "aspect")
    # plot(aspect_values)
    aspect_df <- as.data.frame(aspect_values)
    aspect_mn <- mean(aspect_df[,1])
    slope_mn_col <- paste0("slope_mn_", buffer_size, "m")
    cov_df[row, slope_mn_col] <- slope_mn

    # Extract Vector Ruggedness Measure
    VRM_values <- spatialEco::vrm(DEM_subset, s = 3) # 3x3 window
    # plot(VRM_values)
    VRM_df <- as.data.frame(VRM_values)
    VRM_mn <- mean(VRM_df[,1])
    VRM_mn_col <- paste0("VRM_mn_", buffer_size, "m")
    cov_df[row, VRM_mn_col] <- VRM_mn

    # Extract VRM Aggregation Index
    # lsm_c_ai: equals the number of like adjacencies divided by the theoretical maximum possible number of like adjacencies for that class
    # Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
    VRM_ai <- landscapemetrics::calculate_lsm(VRM_subset, what = "lsm_c_ai")
    VRM_ai <- VRM_ai[which(VRM_ai$class == 1),] # VRM
    VRM_ai_col <- paste0("VRM_ai_", buffer_size, "m")
    if (nrow(VRM_ai) > 0) {
      cov_df[row, VRM_ai_col] <- VRM_ai$value
    } else {
      cov_df[row, VRM_ai_col] <- NA
    }
    
    # Extract VRM Patch Area
    # lsm_p_area: Calculates the area of each patch.(m^2)
    VRM_Parea <- landscapemetrics::calculate_lsm(VRM_subset, what = "lsm_p_area")
    VRM_Parea <- VRM_Parea[which(VRM_Parea$class == 1),] # VRM
    VRM_mnParea <- mean(VRM_Parea$value)
    VRM_mnParea_col <- paste0("VRM_mnParea_", buffer_size, "m")
    if (nrow(VRM_ai) > 0) {
      cov_df[row, VRM_mnParea_col] <- VRM_mnParea
    } else {
      cov_df[row, VRM_mnParea_col] <- NA
    }

    # Pastures within buffer
    past_in_buf <- terra::intersect(pasture, site_buffer)
    # plot(site_buffer) # visually inspect
    # plot(past_in_buf, add = TRUE)
    Npastures <- nrow(past_in_buf)
    Npastures_col <- paste0("Npastures_", buffer_size, "m")
    cov_df[row, Npastures_col] <- Npastures

    # Villages within buffer
    vil_in_buf <- terra::intersect(villages, site_buffer)
    # plot(site_buffer) # visually inspect
    # plot(vil_in_buf, add = TRUE)
    Nvillages <- nrow(vil_in_buf)
    Nvillages_col <- paste0("Nvillages_", buffer_size, "m")
    cov_df[row, Nvillages_col] <- Nvillages
    
    # Rivers within buffer
    riv_in_buf <- terra::intersect(rivers, site_buffer)
    # plot(site_buffer) # visually inspect
    # plot(riv_in_buf, add = TRUE)
    Nriver <- nrow(riv_in_buf)
    Nriver_col <- paste0("Nriver_", buffer_size, "m")
    cov_df[row, Nriver_col] <- Nriver
    
    # Update Buffer Progress Bar
    buffer_pb$tick()

  } # end buffer
  
  # Update Site Progress Bar
  pb$tick
  
} # end site
# --------------------------- End Extraction Loop -----------------------------

# Take a look
head(cov_df)
View(cov_df)


# Export
saveRDS(cov_df, "Data/SiteCovariates.rds")

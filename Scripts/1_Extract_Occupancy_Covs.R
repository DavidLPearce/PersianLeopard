
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
VRM <- spatialEco::vrm(DEM, s = 3) # Extract VRM # 3x3 window
plot(VRM)
VRM_rast <- VRM > 0.05 & VRM <= 0.20# Classify VRM values into raster
plot(VRM_rast)



# ------------------------------------------------------------------------------
#
#                         Extract Landscape Metrics
#
# ------------------------------------------------------------------------------

# Sequence of buffers
buffer_vec <-  seq(250, 15000, by = 250) 
sigma_vec <- seq(250, 15000, by = 250) 

# Dataframe for covariates
dist_df <- as.data.frame(sites)
dist_df <- dist_df[,-c(4:5)]
colnames(dist_df)[1] <- "SiteID"
print(dist_df)

elev_df <- cov_df
slope_df <- cov_df
aspect_df <- cov_df
VRM_df <- cov_df
VRM_AI_df <- cov_df
Vrm_Patch_df <- cov_df
Npast_df <- cov_df
Nvil_df <- cov_df
Nriver_df <- cov_df

# Distance-weighted function
weighted_metric <- function(site_coords, raster_df, value_col, sigma) {
  if (nrow(raster_df) == 0) return(NA)
  dists <- sqrt((raster_df$x - site_coords[1])^2 + (raster_df$y - site_coords[2])^2)
  weights <- exp(-(dists^2) / (2 * sigma^2))
  weighted.mean(raster_df[[value_col]], weights, na.rm = TRUE)
}

# Initialize the progress bar
total_steps <- NROW(sites) * length(buffer_vec) * length(sigma_vec)
pb <- progress_bar$new(
  format = "Processing [:bar] :percent | :message | Elapsed: :elapsed",
  total = total_steps,
  clear = FALSE,
  width = 120
)

# row = 1
# buffer = 1
# sigma = 1

# Loop
for (row in 1:NROW(sites)) {

  # Subset site
  site_sub <- sites[row]
  
  # Extract site coordinates
  site_coords <- as.data.frame(site_sub)
  site_coords <- c(site_coords$x[1], site_coords$y[1])
  
  # ---------------------
  # Distance Metrics
  # ---------------------
  
  # Nearest Pasture Distance
  past_dist <- distance(site_sub, pasture)
  nearest_past_dist <- apply(past_dist, 1, min)
  nearest_past_dist_col <- paste0("nrstPastureDist")
  dist_df[row, nearest_past_dist_col] <- nearest_past_dist
  
  # Nearest Village Distance
  vil_dist <- distance(site_sub, villages)
  nearest_vil_dist <- apply(vil_dist, 1, min)
  nearest_vil_dist_col <- paste0("nrstVillageDist")
  dist_df[row, nearest_vil_dist_col] <- nearest_vil_dist
  
  # Nearest River
  riv_dist <- distance(site_sub, rivers)
  nearest_riv_dist <- apply(riv_dist, 1, min)
  nearest_riv_dist_col <- paste0("nrstRiverDist")
  dist_df[row, nearest_riv_dist_col] <- nearest_riv_dist
  
  # ---------------------
  # Buffer Metrics 
  # ---------------------

  for (buffer_idx in seq_along(buffer_vec)){
    
    # buffer size
    buffer_size <- buffer_vec[buffer_idx]
    
    # Print progress message
    if (!pb$finished) {
      pb$message(sprintf("Site: %s | Buffer: %sm", 
                         sites$Field1[row], buffer_size))
      pb$tick()
    }

    # Create a buffer around the site
    site_buffer <- terra::buffer(site_sub, width = buffer_size)
  
    # Subset DEM
    DEM_subset <- terra::crop(DEM, site_buffer)
    # plot(DEM_subset)
    
    # Subset VRM raster
    VRM_subset <- terra::crop(VRM_rast, site_buffer)
    # plot(VRM_subset)
    
    
    # Extract VRM Aggregation Index
    # lsm_c_ai: equals the number of like adjacencies divided by the theoretical maximum possible number of like adjacencies for that class
    # Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
    VRM_ai <- landscapemetrics::calculate_lsm(VRM_subset, what = "lsm_c_ai")
    VRM_ai <- VRM_ai[which(VRM_ai$class == 1),] # VRM
    VRM_ai_col <- paste0("VRM_ai_", buffer_size, "m")
    if (nrow(VRM_ai) > 0) {
      VRM_AI_df[row, VRM_ai_col] <- VRM_ai$value
    } else {
      VRM_AI_df[row, VRM_ai_col] <- NA
    }
    
    # Extract VRM Patch Area
    # lsm_p_area: Calculates the area of each patch.(m^2)
    VRM_Parea <- landscapemetrics::calculate_lsm(VRM_subset, what = "lsm_p_area")
    VRM_Parea <- VRM_Parea[which(VRM_Parea$class == 1),] # VRM
    VRM_Parea <- mean(VRM_Parea$value)
    VRM_Parea_col <- paste0("VRM_Parea_", buffer_size, "m")
    if (nrow(VRM_ai) > 0) {
      Vrm_Patch_df[row, VRM_Parea_col] <- VRM_Parea
    } else {
      Vrm_Patch_df[row, VRM_Parea_col] <- NA
    }
    
    # Pastures within buffer
    past_in_buf <- terra::intersect(pasture, site_buffer)
    # plot(site_buffer) # visually inspect
    # plot(past_in_buf, add = TRUE)
    Npastures <- nrow(past_in_buf)
    Npastures_col <- paste0("Npastures_", buffer_size, "m")
    Npast_df[row, Npastures_col] <- Npastures

    # Villages within buffer
    vil_in_buf <- terra::intersect(villages, site_buffer)
    # plot(site_buffer) # visually inspect
    # plot(vil_in_buf, add = TRUE)
    Nvillages <- nrow(vil_in_buf)
    Nvillages_col <- paste0("Nvillages_", buffer_size, "m")
    Nvil_df[row, Nvillages_col] <- Nvillages
    
    # Rivers within buffer
    riv_in_buf <- terra::intersect(rivers, site_buffer)
    # plot(site_buffer) # visually inspect
    # plot(riv_in_buf, add = TRUE)
    Nriver <- nrow(riv_in_buf)
    Nriver_col <- paste0("Nriver_", buffer_size, "m")
    Nriver_df[row, Nriver_col] <- Nriver
  
    # ---------------------
    # Distance-weighted 
    # ---------------------
    for (sigma_idx  in seq_along(sigma_vec)){
    
      # sigma size
      sigma <- sigma_vec[sigma_idx ]  
      

      # Extract Elevation
      elev_values <- terra::extract(DEM_subset, site_buffer)
      elev_tmp <- terra::extract(DEM_subset, site_buffer, xy = TRUE)
      elev <- weighted_metric(site_coords, elev_tmp, "DEM20kmbuffer", sigma = sigma)
      elev_col <- paste0("elev_", buffer_size, "m_sigma", sigma)
      elev_df[row, elev_col] <- elev
  
      # Extract Slope
      slope_values <- terra::terrain(x = DEM_subset, v = "slope", unit = "degrees")
      # plot(slope_values)
      slope_tmp <- terra::extract(slope_values, site_buffer, xy = TRUE)
      slope <- weighted_metric(site_coords, slope_tmp, "slope", sigma = sigma)
      slope_col <- paste0("slope_", buffer_size, "m_sigma", sigma)
      slope_df[row, slope_col] <- slope
  
      # Extract Aspect
      aspect_values <- terra::terrain(x = DEM_subset, v = "aspect")
      # plot(aspect_values)
      aspect_tmp <- terra::extract(aspect_values, site_buffer, xy = TRUE)
      aspect <- weighted_metric(site_coords, aspect_tmp, "aspect", sigma = sigma)
      aspect_col <- paste0("aspect_", buffer_size, "m_sigma", sigma)
      aspect_df[row, aspect_col] <- aspect
  
      # Extract Vector Ruggedness Measure
      VRM_values <- spatialEco::vrm(DEM_subset, s = 3) # 3x3 window
      # plot(VRM_values)
      VRM_tmp <- terra::extract(VRM_values, site_buffer, xy = TRUE)
      VRM <- weighted_metric(site_coords, VRM_df, "lyr1", sigma = sigma)
      VRM_col <- paste0("VRM_", VRM_tmp, "m_sigma", sigma)
      VRM_df[row, VRM_col] <- VRM

    } # end sigma
  } # end buffer
} # end site
# --------------------------- End Extraction Loop -----------------------------

# Take a look
head(cov_df)




# Export
saveRDS(cov_df, "Data/OccuCovData/cov_df.rds")



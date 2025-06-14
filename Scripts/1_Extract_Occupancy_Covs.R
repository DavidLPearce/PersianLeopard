
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
library(landscapemetrics)
library(progress)
packageVersion("landscapemetrics")
set.seed(123)
options(scipen = 9999)
setwd(".")

# ------------------------------------------------------------------------------
#
#                         Read in data
#
# ------------------------------------------------------------------------------


# Camera sites
sites <- vect("./Data/SpatialData/CameraLocationsNew/CameraLocationsNew.shp")

# Elevation
DEM <- rast("./Data/SpatialData/DEM20kmBuffer/DEM20kmbuffer.tif")

# # Soil Adjusted Vegetation Index
# SAVI <- rast("./Data/SpatialData/SoilAdjustedVegitationIndex/SAVI_Raster.tif")

# Rivers
rivers <- vect("./Data/SpatialData/Rivers/riverclip.shp")

# Pastures
pasture <- vect("./Data/SpatialData/PastureLocations/pastureloc.shp")

# Villages
villages <- vect("./Data/SpatialData/VillageLocations/vilage.shp")

# Study area boundary
study_area <- vect("D:/Persian_Leopard_GIS/StudyAreaBoundary/marz1.shp")


# ------------------------------------------------------------------------------
#
#                         Data Wrangling
#
# ------------------------------------------------------------------------------


# Check the CRS of each layer
print(sites)
print(DEM)
# print(SAVI)
print(rivers)
print(pasture)
print(villages)

# Reproject to match DEM's CRS
# SAVI <- project(SAVI, crs(DEM)) 
study_area <- project(study_area, crs(DEM))

# Visually inspect
plot(DEM)
# plot(SAVI)
plot(sites, add = TRUE, col = "black")
plot(rivers, add = TRUE, col = "blue")
plot(pasture, add = TRUE, col = "green")
plot(villages, add = TRUE, col = "red")
plot(study_area, add = TRUE, border  = "yellow")
dev.off() # Clear

# # Vector Ruggedness Measure Local
# # --------------------------------------------------
# # See Dilts et al. 2022
# # https://link.springer.com/article/10.1007/s10980-023-01646-6
# 
# # Note: VRM is scale dependent and needs to be calculated at each scale i.e., each buffer
# # But as an example:
# 
# # # Smooth DEM and create a raster of surface variation
# w <- matrix(1, nrow = 3, ncol = 3)
# DEM_smooth <- focal(DEM, w = w, fun = mean, na.policy = "omit", na.rm = TRUE)
# surf_vari_rast <- DEM - DEM_smooth
# plot(surf_vari_rast)
# 
# # Calculate Vector Ruggedness Measure using surface variation
# VRML <- spatialEco::vrm(surf_vari_rast, s=3) # Extract VRM # 3x3 window
# plot(VRML)
# dev.off() # Clear
# # writeRaster(VRML, "./Data/SpatialData/VRM/VRML.tif") # Export
# 
# # Top 10% of cells as rugged
# vals <- values(VRML, mat = FALSE, na.rm = TRUE)
# threshold <- quantile(vals, probs = 0.9, na.rm = TRUE)
# VRML_rugged <- VRML > threshold
# plot(VRML_rugged)
# dev.off() # Clear
# # writeRaster(VRML_rugged, "./Data/SpatialData/VRM/VRML_rugged.tif") # Export
# # --------------------------------------------------

# ------------------------------------------------------------------------------
#
#                         Extract Landscape Metrics
#
# ------------------------------------------------------------------------------

# Dataframe for covariates
# Separate since each covariate will be fit in 
# a separate occupancy model to determine scale
dist_df <- as.data.frame(sites)
dist_df <- dist_df[,-c(4:5)]
colnames(dist_df)[1] <- "SiteID"
print(dist_df)

elev_MN_df <- dist_df
elev_SD_df <- dist_df
slope_MN_df <- dist_df
slope_SD_df <- dist_df
aspect_SINE_df <- dist_df
aspect_COS_df <- dist_df
VRML_MN_df <- dist_df
VRML_PRP_df <- dist_df
VRML_PD_df <- dist_df
VRML_PA_df <- dist_df
VRML_COH_df<- dist_df
Npast_df <- dist_df
Nvil_df <- dist_df
Nriver_df <- dist_df
# SAVI_MN_df <- dist_df
# SAVI_SD_df <- dist_df


# Define buffers
buffer_vec <- seq(250, 15000, by = 250)

# Function for smooth filter
w <- matrix(1, nrow = 3, ncol = 3)

# Progress bar
total_steps <- NROW(sites) * length(buffer_vec)
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent | :site | Buffer: :buffer | Elapsed: :elapsed| ETA: :eta",
  total = total_steps,
  clear = FALSE,
  width = 120
)

# row = 1
# buffer_idx = 1

# ---------------------
# Extract Metrics Loop
# ---------------------
for (row in 1:NROW(sites)) {
  
  # Subset site
  site_sub <- sites[row]
  
  site_id <- as.data.frame(site_sub)
  site_id <- site_id[, "Field1"]
  
  # ---------------------
  # Distance Metrics
  # ---------------------
  
  # ---- Nearest Pasture Distance
  past_dist <- distance(site_sub, pasture)
  nearest_past_dist <- apply(past_dist, 1, min)
  nearest_past_dist_col <- paste0("nrstPastureDist")
  dist_df[row, nearest_past_dist_col] <- nearest_past_dist

  # ---- Nearest Village Distance
  vil_dist <- distance(site_sub, villages)
  nearest_vil_dist <- apply(vil_dist, 1, min)
  nearest_vil_dist_col <- paste0("nrstVillageDist")
  dist_df[row, nearest_vil_dist_col] <- nearest_vil_dist

  # ---- Nearest River
  riv_dist <- distance(site_sub, rivers)
  nearest_riv_dist <- apply(riv_dist, 1, min)
  nearest_riv_dist_col <- paste0("nrstRiverDist")
  dist_df[row, nearest_riv_dist_col] <- nearest_riv_dist
  
  # ---------------------
  # Buffer Metrics 
  # ---------------------
  for (buffer_idx in seq_along(buffer_vec)) {
    
    # Get Buffer size
    buffer_size <- buffer_vec[buffer_idx]
      
    # Progress bar message
    pb$tick(tokens = list(site = sites$Field1[row], buffer = buffer_size))

    # Create a buffer around the site
    site_buffer <- terra::buffer(site_sub, width = buffer_size)

    # Subset DEM
    DEM_subset <- terra::crop(DEM, site_buffer)
    # plot(DEM_subset)
    
    # Smoothed DEM
    DEM_subset_smooth <- focal(DEM_subset, w = w, fun = mean, na.policy = "omit", na.rm = TRUE)
    surf_vari_subset <- DEM_subset - DEM_subset_smooth
    # plot(surf_vari_subset)
    
    # Calculate VRM
    VRML_subset <- spatialEco::vrm(surf_vari_subset, s=3)
    # plot(VRML_subset)
    
    # Ruggedness raster
    VRML_vals <- values(VRML_subset, mat = FALSE, na.rm = TRUE)
    VRML_threshold <- quantile(VRML_vals, probs = 0.9, na.rm = TRUE)
    VRML_rugged_subset <- VRML_subset > VRML_threshold
    # plot(VRML_rugged_subset)
    
    # # Subset SAVI
    # SAVI_subset <- terra::crop(SAVI, site_buffer)
    # # plot(SAVI_subset)

    # ---- Extract Elevation
    elev_MN_value <- terra::extract(DEM_subset, site_buffer, fun = mean, na.rm = TRUE) # Mean
    elev_MN_col <- paste0("elev_MN_buffer", buffer_size)
    elev_MN_df[row, elev_MN_col] <- elev_MN_value[,"DEM20kmbuffer"]

    elev_SD_value <- terra::extract(DEM_subset, site_buffer, fun = sd, na.rm = TRUE) # Standard Deviation
    elev_SD_col <- paste0("elev_SD_buffer", buffer_size)
    elev_SD_df[row, elev_SD_col] <- elev_SD_value[,"DEM20kmbuffer"]

    # ---- Extract Slope
    slope_subset<- terrain(DEM_subset, v = "slope", unit = "degrees")

    slope_MN_value <- terra::extract(slope_subset, site_buffer, fun = mean, na.rm = TRUE)# Mean
    slope_MN_col <- paste0("slope_MN_buffer", buffer_size)
    slope_MN_df[row, slope_MN_col] <- slope_MN_value[,"slope"]

    slope_SD_value <- terra::extract(slope_subset, site_buffer, fun = sd, na.rm = TRUE) # Standard Deviation
    slope_SD_col <- paste0("slope_SD_buffer", buffer_size)
    slope_SD_df[row, slope_SD_col] <- slope_SD_value[,"slope"]

    # ---- Extract Aspect
    aspect_subset <- terra::terrain(x = DEM_subset, v = "aspect")
    aspect_values <- as.data.frame(aspect_subset)
    aspect_values_rad <- aspect_values * pi / 180 # radians

    aspect_sine <- mean(sin(aspect_values_rad[,1])) # Sine (Eastness)
    aspect_SINE_col <- paste0("aspect_SINE_buffer", buffer_size)
    aspect_SINE_df[row, aspect_SINE_col] <- aspect_sine

    aspect_cos <- mean(cos(aspect_values_rad[,1]))# Cosine (Northness)
    aspect_COS_col <- paste0("aspect_COS_buffer", buffer_size)
    aspect_COS_df[row, aspect_COS_col] <- aspect_cos

    # ---- Extract Vector Ruggedness Measure Local (VRML)
    VRML_MN_value <- terra::extract(VRML_subset, site_buffer, fun = mean, na.rm = TRUE)
    VRML_MN_col <- paste0("VRML_MN_buffer", buffer_size)
    VRML_MN_df[row, VRML_MN_col] <- VRML_MN_value[, "lyr1"]
    
    # ---- Extract VRML Proportion
    VRML_PRP_value <- terra::extract(VRML_rugged_subset, site_buffer, fun = mean, na.rm = TRUE)
    VRML_PRP_col <- paste0("VRML_PRP_buffer", buffer_size)
    VRML_PRP_df[row, VRML_PRP_col] <- VRML_PRP_value[, "lyr1"]

    # ---- Extract VRML Patch Density
    VRML_PD_value <- landscapemetrics::calculate_lsm(VRML_rugged_subset, what = "lsm_c_pd")
    VRML_PD_value <- VRML_PD_value[which(VRML_PD_value$class == 1),] # 1 = rugged
    VRML_PD_col <- paste0("VRML_PD_buffer", buffer_size)
    if (nrow(VRML_PD_value) > 0 && any(VRML_PD_value[["value"]] > 0, na.rm = TRUE)) {
      VRML_PD_df[row, VRML_PD_col] <- VRML_PD_value[["value"]]
    } else {
      VRML_PD_df[row, VRML_PD_col] <- NA
    }

    # ---- Extract VRML Patch Area
    VRML_PA_value <- landscapemetrics::calculate_lsm(VRML_rugged_subset, what = "lsm_p_area")
    VRML_PA_value <- VRML_PA_value[which(VRML_PA_value$class == 1),] # 1 = rugged
    VRML_PA_value <- mean(VRML_PA_value$value, na.rm = TRUE)
    VRML_PA_col <- paste0("VRML_PA_buffer", buffer_size)
    if (!is.na(VRML_PA_value) && VRML_PA_value > 0) {
      VRML_PA_df[row, VRML_PA_col] <- VRML_PA_value
    } else {
      VRML_PA_df[row, VRML_PA_col] <- NA
    }

    # ---- Extract VRML Cohesion Index
    VRML_COH_value <- landscapemetrics::calculate_lsm(VRML_rugged_subset, what = "lsm_c_cohesion")
    VRML_COH_value <- VRML_COH_value[which(VRML_COH_value$class == 1),] # 1 = rugged
    VRML_COH_col <- paste0("VRML_COH_buffer", buffer_size)
    if (VRML_COH_value[, "value"] > 0) {
      VRML_COH_df[row, VRML_COH_col] <- VRML_COH_value[,"value"]
    } else {
      VRML_COH_df[row, VRML_COH_col] <- NA
    }

    # ---- Pastures within buffer
    past_in_buf <- terra::intersect(pasture, site_buffer)
    Npastures <- nrow(past_in_buf)
    Npastures_col <- paste0("Npastures_buffer", buffer_size)
    Npast_df[row, Npastures_col] <- Npastures

    # ---- Villages within buffer
    vil_in_buf <- terra::intersect(villages, site_buffer)
    Nvillages <- nrow(vil_in_buf)
    Nvillages_col <- paste0("Nvillages_buffer", buffer_size)
    Nvil_df[row, Nvillages_col] <- Nvillages

    # ---- Rivers within buffer
    riv_in_buf <- terra::intersect(rivers, site_buffer)
    Nriver <- nrow(riv_in_buf)
    Nriver_col <- paste0("Nriver_buffer", buffer_size)
    Nriver_df[row, Nriver_col] <- Nriver

    # # ---- Soil Adjusted Vegetation Index
    # SAVI_MN_value <- extract(SAVI_subset, site_buffer, fun = mean, na.rm = TRUE) # Mean
    # SAVI_MN_col <- paste0("SAVI_MN_buffer", buffer_size)
    # SAVI_MN_df[row, SAVI_MN_col] <- SAVI_MN_value[,"SAVI"]
    # 
    # SAVI_SD_value <- extract(SAVI_subset, site_buffer, fun = sd, na.rm = TRUE) # Standard Deviation
    # SAVI_SD_col <- paste0("SAVI_SD_buffer", buffer_size)
    # SAVI_SD_df[row, SAVI_SD_col] <- SAVI_SD_value[,"SAVI"]

  } # End Buffer Loop
} # End Site Loop
# ----- End Loop -----  

# Take a look
head(dist_df, 5)
head(elev_MN_df, 5)
head(elev_SD_df, 5)
head(slope_MN_df, 5)
head(slope_SD_df, 5)
head(aspect_SINE_df, 5)
head(aspect_COS_df, 5)
head(VRML_MN_df, 5)
head(VRML_PRP_df, 5)
head(VRML_PD_df, 5)
head(VRML_PA_df, 5)
head(VRML_COH_df, 5)
head(Npast_df, 5)
head(Nvil_df, 5)
head(Nriver_df, 5)
# head(SAVI_MN_df, 5)
# head(SAVI_SD_df, 5)

# Export
saveRDS(dist_df, "Data/OccuCovData/dist_df.rds")
saveRDS(elev_MN_df, "Data/OccuCovData/elev_MN_df.rds")
saveRDS(elev_SD_df, "Data/OccuCovData/elev_SD_df.rds")
saveRDS(slope_MN_df, "Data/OccuCovData/slope_MN_df.rds")
saveRDS(slope_SD_df, "Data/OccuCovData/slope_SD_df.rds")
saveRDS(aspect_SINE_df, "Data/OccuCovData/aspect_SINE_df.rds")
saveRDS(aspect_COS_df, "Data/OccuCovData/aspect_COS_df.rds")
saveRDS(VRML_MN_df, "Data/OccuCovData/VRML_MN_df.rds")
saveRDS(VRML_PRP_df, "Data/OccuCovData/VRML_PRP_df.rds")
saveRDS(VRML_PD_df, "Data/OccuCovData/VRML_PD_df.rds")
saveRDS(VRML_PA_df, "Data/OccuCovData/VRML_PA_df.rds")
saveRDS(VRML_COH_df, "Data/OccuCovData/VRML_COH_df.rds")
saveRDS(Npast_df, "Data/OccuCovData/Npast_df.rds")
saveRDS(Nvil_df, "Data/OccuCovData/Nvil_df.rds")
saveRDS(Nriver_df, "Data/OccuCovData/Nriver_df.rds")
saveRDS(SAVI_MN_df, "Data/OccuCovData/SAVI_MN_df.rds")
saveRDS(SAVI_SD_df, "Data/OccuCovData/SAVI_SD_df.rds")

# ------------------- End of Script -------------------
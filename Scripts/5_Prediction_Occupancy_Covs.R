
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

# Elevation
DEM <- rast("./Data/SpatialData/DEM20kmBuffer/DEM20kmbuffer.tif")

# Study area boundary
study_area <- vect("E:/Persian_Leopard_GIS/StudyAreaBoundary/marz1.shp")

# ------------------------------------------------------------------------------
#
#                         Data Wrangling
#
# ------------------------------------------------------------------------------


# Check the CRS of each layer
print(DEM)
print(study_area)

# Reproject study_area to match DEM's CRS
study_area <- project(study_area, crs(DEM))

# Visually inspect
plot(DEM)
plot(study_area, add = TRUE, border  = "yellow")

# Buffer by 10,000 meters (10 km)
study_area_buff <- buffer(study_area, width = 10000)

# Plot to verify
plot(study_area_buff, add = TRUE, border = "red")

study_area <- study_area_buff

# Create prediction sites sites
grid_res <- 1000  # m spacing
template <- rast(ext(DEM), resolution = grid_res, crs = crs(DEM))

# Convert to grid points
pred_sites <- as.points(template)

# Subset points to  study area
pred_sites <- pred_sites[study_area, ]

# Plot
plot(DEM)
plot(pred_sites, pch = 16, cex = 0.3, col = "red", add = TRUE)

# Add point ID 
pred_sites$PointID <- 1:nrow(pred_sites)

# Creating a dataframe to hold covariate values
coords <- crds(pred_sites, df = TRUE)   
PointID   <- as.data.frame(pred_sites)      
pred_cov_df <- cbind(PointID, coords)
head(pred_cov_df, 25)
nrow(pred_cov_df)

# ------------------------------------------------------------------------------
#
#                         Extract Landscape Metrics
#
# ------------------------------------------------------------------------------
 
# Best model was 
# psi ~ aspect_SINE_buffer4000 + aspect_COS_buffer5000 + VRML_COH_buffer750

# Define scale
sine_buffer <- 4000
cos_buffer <- 5000
vrml_buffer <- 750

# Function for smooth filter
w <- matrix(1, nrow = 3, ncol = 3)

# Progress bar
pb <- progress_bar$new(
  format = "Site :current of :total | [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
  total = NROW(pred_cov_df),
  clear = FALSE,
  width = 100
)

# ---------------------
# Extract Metrics Loop
# ---------------------
for (row in 1:NROW(pred_cov_df)) {
  
  # Subset site
  site_sub <- pred_cov_df[row,]
  
  # Convert data.frame to a SpatVector  
  site_sub <- vect(site_sub, geom = c("x", "y"), crs = crs(DEM))

  # Create a buffer around the site
  aspect_sine_scale <- terra::buffer(site_sub, width = sine_buffer) # Aspect Sine
  aspect_cos_scale <- terra::buffer(site_sub, width = cos_buffer) # Aspect Cosine
  vrml_coh_scale <- terra::buffer(site_sub, width = vrml_buffer)# VRML Cohesion Index
  
  # Subset DEM
  DEM_aspect_sine <- terra::crop(DEM, aspect_sine_scale)
  DEM_aspect_cos <- terra::crop(DEM, aspect_sine_scale)
  DEM_vrml_coh <- terra::crop(DEM, aspect_sine_scale)

  # Smoothed DEM
  DEM_subset_smooth <- focal(DEM_vrml_coh, w = w, fun = mean, na.policy = "omit", na.rm = TRUE)
  surf_vari_subset <- DEM_vrml_coh - DEM_subset_smooth
  # plot(surf_vari_subset)
  
  # Calculate VRM
  VRML_subset <- spatialEco::vrm(surf_vari_subset, s=3)
  # plot(VRML_subset)
  
  # Ruggedness raster
  VRML_vals <- values(VRML_subset, mat = FALSE, na.rm = TRUE)
  VRML_threshold <- quantile(VRML_vals, probs = 0.9, na.rm = TRUE)
  VRML_rugged_subset <- VRML_subset > VRML_threshold
  # plot(VRML_rugged_subset)

  # ---- Extract Aspect Sine
  aspect_sine_subset <- terra::terrain(x = DEM_aspect_sine, v = "aspect")
  aspect_sine_values <- as.data.frame(aspect_sine_subset)
  aspect_sine_values_rad <- aspect_sine_values * pi / 180 # radians
  
  aspect_sine <- mean(sin(aspect_sine_values_rad[,1])) # Sine (Eastness)
  aspect_SINE_col <- paste0("aspect_SINE_buffer4000")
  pred_cov_df[row, aspect_SINE_col] <- aspect_sine
  
  # ---- Extract Aspect Cosine
  aspect_cos_subset <- terra::terrain(x = DEM_aspect_cos, v = "aspect")
  aspect_cos_values <- as.data.frame(aspect_cos_subset)
  aspect_cos_values_rad <- aspect_cos_values * pi / 180 # radians
  
  aspect_cos <- mean(cos(aspect_cos_values_rad[,1]))# Cosine (Northness)
  aspect_COS_col <- paste0("aspect_COS_buffer5000")
  pred_cov_df[row, aspect_COS_col] <- aspect_cos
  
  # ---- Extract VRML Cohesion Index
  VRML_COH_value <- landscapemetrics::calculate_lsm(VRML_rugged_subset, what = "lsm_c_cohesion")
  VRML_COH_value <- VRML_COH_value[which(VRML_COH_value$class == 1),] # 1 = rugged
  VRML_COH_col <- paste0("VRML_COH_buffer750")
  if (VRML_COH_value[, "value"] > 0) {
    pred_cov_df[row, VRML_COH_col] <- VRML_COH_value[,"value"]
  } else {
    pred_cov_df[row, VRML_COH_col] <- NA
  }
  
  # Update progress bar
  pb$tick(tokens = list(current = row))
    
} # End Site Loop
# ----- End Loop -----   


# Take a look
head(pred_cov_df, 25)


# Save predictions
saveRDS(pred_cov_df, "Data/OccuCovData/pred_cov_df.rds")

# ------------------- End of Script -------------------

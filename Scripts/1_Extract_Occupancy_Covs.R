
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
VRM_rast <- VRM > 0.03 & VRM <= 1 # Classify VRM values into raster
plot(VRM_rast)



# ------------------------------------------------------------------------------
#
#                         Extract Landscape Metrics
#
# ------------------------------------------------------------------------------

# Dataframe for covariates
dist_df <- as.data.frame(sites)
dist_df <- dist_df[,-c(4:5)]
colnames(dist_df)[1] <- "SiteID"
print(dist_df)

elev_df <- dist_df
slope_df <- dist_df
aspect_df <- dist_df
VRM_df <- dist_df
VRM_AI_df <- dist_df
VRM_Patch_df <- dist_df
Npast_df <- dist_df
Nvil_df <- dist_df
Nriver_df <- dist_df

# Define concentric buffers aka rings
ring_vec <- seq(250, 15000, by = 250)

# Initialize the progress bar
total_steps <- NROW(sites) * (length(ring_vec) - 1)
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent | Site :site Ring :ring",
  total = total_steps, clear = FALSE, width = 80
)

row = 1
ring_idx = 2

# ---------------------
# Extract Metrics
# ---------------------
for (row in 1:NROW(sites)) {
  
  # Subset site
  site_sub <- sites[row]
  
  site_id <- as.data.frame(site_sub)
  site_id <- site_id[, "Field1"]
  
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
  # Ring Metrics 
  # ---------------------
  for (ring_idx in 2:length(ring_vec)) {
    
    # Extract ring limits
    inner_r <- ring_vec[ring_idx - 1]
    outer_r <- ring_vec[ring_idx]
        
    # Progress Bar Update
    pb$tick(tokens = list(site = site_id, 
                          ring = paste0(inner_r, "-", outer_r, "m")))
    
    # Create ring
    outer_buffer <- terra::buffer(site_sub, width = outer_r)
    inner_buffer <- terra::buffer(site_sub, width = inner_r)
    ring_polygon <- erase(outer_buffer, inner_buffer)
    # plot(ring_polygon)

    # Subset DEM
    DEM_subset <- terra::crop(DEM, outer_buffer)
    # plot(DEM_subset)
    
    # Subset VRM raster
    VRM_subset <- terra::crop(VRM_rast, outer_buffer)
    # plot(VRM_subset)
    VRM_ring <- mask(crop(VRM_subset, ring_polygon), ring_polygon)
    # plot(VRM_ring)
    
    # Extract Elevation
    elev_value <- extract(DEM_subset, ring_polygon, fun = mean, na.rm = TRUE)
    elev_col <- paste0("elev_", outer_r)
    elev_df[row, elev_col] <- elev_value[,"DEM20kmbuffer"]
    
    # Extract Slope
    slope <- terrain(DEM_subset, v = "slope", unit = "degrees")
    slope_value <- extract(slope, ring_polygon, fun = mean, na.rm = TRUE)
    slope_col <- paste0("slope_", outer_r)
    slope_df[row, slope_col] <- slope_value[,"slope"]
    
    
    # Extract Aspect
    aspect <- terra::terrain(x = DEM_subset, v = "aspect")
    aspect_value <- extract(aspect, ring_polygon, fun = mean, na.rm = TRUE)
    aspect_col <- paste0("aspect_", outer_r)
    aspect_df[row, aspect_col] <- aspect_value[,"aspect"]
    
    # Extract Vector Ruggedness Measure
    VRM <- spatialEco::vrm(DEM_subset, s = 3) # 3x3 window
    # plot(VRM)
    VRM_value <- extract(VRM, ring_polygon, fun = mean, na.rm = TRUE)
    VRM_col <- paste0("VRM_", outer_r)
    VRM_df[row, VRM_col] <- VRM_value[, "lyr1"]
    
    # Extract VRM Aggregation Index
    # lsm_c_ai: equals the number of like adjacencies divided by the theoretical maximum possible number of like adjacencies for that class
    # Equals 0 for maximally disaggregated and 100 for maximally aggregated classes.
    VRM_ai <- landscapemetrics::calculate_lsm(VRM_ring, what = "lsm_c_ai")
    VRM_ai <- VRM_ai[which(VRM_ai$class == 1),] # VRM
    VRM_ai_col <- paste0("VRM_AI_", outer_r)
    if (nrow(VRM_ai) > 0) {
      VRM_AI_df[row, VRM_ai_col] <- VRM_ai[, "value"]
    } else {
      VRM_AI_df[row, VRM_ai_col] <- NA
    }
    
    # Extract VRM Patch Area
    # lsm_p_area: Calculates the area of each patch.(m^2)
    VRM_Parea <- landscapemetrics::calculate_lsm(VRM_ring, what = "lsm_p_area")
    VRM_Parea <- VRM_Parea[which(VRM_Parea$class == 1),] # VRM
    VRM_Parea <- mean(VRM_Parea$value)
    VRM_Parea_col <- paste0("VRM_Parea_", outer_r)
    if (nrow(VRM_ai) > 0) {
      VRM_Patch_df[row, VRM_Parea_col] <- VRM_Parea
    } else {
      VRM_Patch_df[row, VRM_Parea_col] <- NA
    }
    
    # Pastures within buffer
    past_in_buf <- terra::intersect(pasture, ring_polygon)
    # plot(ring_polygon) # visually inspect
    # plot(past_in_buf, add = TRUE)
    Npastures <- nrow(past_in_buf)
    Npastures_col <- paste0("Npastures_", outer_r)
    Npast_df[row, Npastures_col] <- Npastures
    
    # Villages within buffer
    vil_in_buf <- terra::intersect(villages, ring_polygon)
    # plot(ring_polygon) # visually inspect
    # plot(vil_in_buf, add = TRUE)
    Nvillages <- nrow(vil_in_buf)
    Nvillages_col <- paste0("Nvillages_", outer_r)
    Nvil_df[row, Nvillages_col] <- Nvillages
    
    # Rivers within buffer
    riv_in_buf <- terra::intersect(rivers, ring_polygon)
    # plot(ring_polygon) # visually inspect
    # plot(riv_in_buf, add = TRUE)
    Nriver <- nrow(riv_in_buf)
    Nriver_col <- paste0("Nriver_", outer_r)
    Nriver_df[row, Nriver_col] <- Nriver

  } # end ring loop
} # end site loop

# --------------------------- End Extraction Loop -----------------------------

# Take a look
head(dist_df)
head(elev_df)
head(slope_df)
head(aspect_df)
head(VRM_df)
head(VRM_AI_df)
head(VRM_Patch_df)
head(Npast_df)
head(Nvil_df)
head(Nriver_df)



# Export
saveRDS(dist_df, "Data/OccuCovData/dist_df.rds")
saveRDS(elev_df, "Data/OccuCovData/elev_df.rds")
saveRDS(slope_df, "Data/OccuCovData/slope_df.rds")
saveRDS(aspect_df, "Data/OccuCovData/aspect_df.rds")
saveRDS(VRM_df, "Data/OccuCovData/VRM_df.rds")
saveRDS(VRM_AI_df, "Data/OccuCovData/VRM_AI_df.rds")
saveRDS(VRM_Patch_df, "Data/OccuCovData/VRM_Patch_df.rds")
saveRDS(Npast_df, "Data/OccuCovData/Npast_df.rds")
saveRDS(Nvil_df, "Data/OccuCovData/Nvil_df.rds")
saveRDS(Nriver_df, "Data/OccuCovData/Nriver_df.rds")


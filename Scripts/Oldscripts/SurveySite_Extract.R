
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
cam_dat <- read.csv("./KansasCamera_data.csv")

# Landcover
lndcvr <- rast("./kansasraster.tif")

# Elevation
dem <- rast("./KS_DEM.tif")

# Streams
stream <- vect("./Streams/AllStreams.shp")


# ------------------------------------------------------------------------------
#
#                         Extract data by site
#
# ------------------------------------------------------------------------------

# Getting site coordinates
site_coords <- unique(cam_dat[, c("Site", "Latitude", "Longitude")])

# Convert siteXY to spatial points
site_coords_sp <- SpatialPoints(coords = site_coords[, c("Longitude", "Latitude")],
                             proj4string = CRS("+proj=longlat +datum=WGS84"))

# Convert the coordinates to a SpatVector object
site_pointvec <- vect(site_coords_sp, type = "points", crs = crs(dem))


# continuous to discrete
lndcvr <- as.factor(lndcvr)

# reproject landcover raster to lat and long
# for extent
lndcvr_reproj <- project(lndcvr, crs(dem))


# Extracting elevation
elev_values <- terra::extract(dem, site_pointvec)
elev_values <- as.data.frame(elev_values)
elev_values$Site <- site_coords$Site
head(elev_values)

# Adding elevation data to site_covs
site_covs <- site_coords
site_covs$Elev <- elev_values$KS_DEM





site = 293
############
# Loop
############

# Initialize the progress bar
pb <- progress_bar$new(
  format = "  Progress [:bar] :percent in :elapsed",
  total = NROW(site_coords$Site),    # Total number of iterations
  width = 60               # Width of the progress bar
)


# Loop to extract slope, aspect, distances
for (site in 1:NROW(site_coords[,'Site'])) {

  # Update the progress bar
  pb$tick()

  # subset the site
  site_sub <- site_coords[which(site_coords[,'Site'] == site),]

  # setting projection
  site_sub_coords <- SpatialPoints(coords = site_sub[, c("Longitude", "Latitude")],
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))

  # Convert the coordinates to a SpatVector object
  site_sub_vec <- vect(site_sub_coords, type = "points", crs = crs(dem))

  # Transform the CRS of site_pointvec to match the CRS of stream
  site_sub_vec_stream <- terra::project(site_sub_vec, crs(stream))

  # buffer the site
  site_buffer <- terra::buffer(site_sub_coords, width = 3000) # adjust as needed

  # subset DEM
  dem_subset <- terra::crop(dem, site_buffer)

  # get slope
  slope_sub <- terra::terrain(x = dem_subset, v = "slope", unit = "degrees")

  # get aspect
  aspect_sub <- terra::terrain(x = dem_subset, v = "aspect")

  # get vector ruggedness
  vrm_sub <- spatialEco::vrm(dem_subset, s = 3)

  # get distances to stream
  stream_distances <- terra::distance(site_sub_vec_stream, stream)

  # extract slope
  slope_values <- terra::extract(slope_sub, site_sub_vec)

  # extract aspect
  aspect_values <- terra::extract(aspect_sub, site_sub_vec)

  # extract vector ruggedness values
  vrm_values <- terra::extract(vrm_sub, site_pointvec)

  # extract min stream distance
  stream_distances_min <- min(stream_distances)

  # adding  values to site_covs
  site_covs[site,'Slope'] <- slope_values[,'slope']
  site_covs[site,'Aspect'] <- aspect_values[,'aspect']
  site_covs[site,'VRM'] <- vrm_values[,'lyr1']
  site_covs[site,'StreamsDist'] <- stream_distances_min
}

# take a look
head(site_covs)
tail(site_covs)


View(site_covs)

library(terra)

# Load DEM

# Elevation
dem <- rast("./Data/SpatialData/DEM20kmBuffer/DEM20kmbuffer.tif")

# Compute curvature
flow <- terrain(dem, v = "flowdir", unit = "radians")

# plot
plot(dem, main = "DEM")
plot(flow, main = "flowdir")





library(sf)

# Replace with your actual file path (just the .shp file is fine)
locations <- st_read("D:/Persian_Leopard_GIS/CameraLocations/Camerloc.shp")

# View the data
print(locations)

# Assign CRS: WGS84 / UTM Zone 40N (EPSG:32640)
st_crs(locations) <- 32640
print(locations)


plot(locations["geometry"])

# Save as ESRI Shapefile
st_write(locations, "D:/Persian_Leopard_GIS/CameraLocationsNew/CameraLocationsNew.shp", delete_layer = TRUE)

# Libraries
library(raster)
library(rgdal)
library(here)

# Load sediment data (Reprojected in ArcGIS to WGS84, reduced raster size to northern region)
grain_size <- brick(here("data/Environmental_data/sediment", "Reproject_Clip_Mean_Grain.tif"))

# Project the data
project_grain <- rasterToPoints(grain_size, spatial = T)
proj4string(project_grain)
projection <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
new_grain <- spTransform(project_grain, CRS(projection))
proj4string(new_grain)
new_grain@data <- data.frame(new_grain@data, long = coordinates(new_grain)[,1],
                             lat = coordinates(new_grain)[,2])
grain_df <- as.data.frame(new_grain)
grain_df <- grain_df[-c(4,5)]

save(grain_df, file = "grain_df.Rdata")

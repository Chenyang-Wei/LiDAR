# This script is based on the NEON tutorial series titled 
#   "Introduction to Light Detection and Ranging (LiDAR) – 
#   Explore Point Clouds and Work with LiDAR Raster Data in R"
#   (https://www.neonscience.org/resources/learning-hub/tutorials/
#   introduction-light-detection-and-ranging-lidar-explore-point).
# 
# Last updated: 1/26/2024.


# Create a lidar-derived Canopy Height Model (CHM) ------------------------

# Load needed packages.
library(raster)
library(rgdal)

# Set working directory.
wd <- "C:/Research_Projects/LiDAR"
setwd(wd)

# Assign raster to object.
dsm <- raster(
  file.path(
    wd, "NEON Tutorials", "Data",
    "NEON-DS-Field-Site-Spatial-Data/SJER/DigitalSurfaceModel/SJER2013_DSM.tif"))

# View info about the raster.
dsm

# Plot the DSM.
plot(dsm, 
     main = "Lidar Digital Surface Model\nSJER, California")

# Import the digital terrain model.
dtm <- raster(
  file.path(
    wd, "NEON Tutorials", "Data",
    "NEON-DS-Field-Site-Spatial-Data/SJER/DigitalTerrainModel/SJER2013_DTM.tif"))
dtm

plot(dtm, 
     main = "Lidar Digital Terrain Model\nSJER, California")

# Use raster math to create CHM.
chm <- dsm - dtm

# View CHM attributes.
chm

plot(chm, 
     main = "Lidar Canopy Height Model\nSJER, California")

# Convert meters to feet.
metersToFeet <- function(meters) {
  return(meters * 3.28084)
}
chm_inFeet <- metersToFeet(chm)

plot(chm_inFeet, 
     main = "Lidar Canopy Height Model\nSJER, California (in feet)")

# Create a function that subtracts one raster from another.
canopyCalc <- function(DSM, DTM) {
  return(DSM - DTM)
}

# Use the function to create the final CHM.
chm2 <- canopyCalc(dsm, dtm)
chm2

# Or use the "overlay" function.
chm3 <- overlay(dsm, dtm, fun = canopyCalc)
chm3

# Write out the CHM in tiff format.
dir.create(file.path(
  wd, "NEON Tutorials", "Results"))
writeRaster(chm,
            file.path(
              wd, "NEON Tutorials", "Results",
              "chm_SJER.tif"),
            "GTiff")


# Extract Values from a Raster in R ---------------------------------------



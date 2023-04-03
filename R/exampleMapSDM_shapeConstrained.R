########################################################################################
# Importing and mapping shape-constrained SDM predictions
# RW: Modified from code by Barb Muhling downloaded from FSII Google Drive 3/28/2023
########################################################################################

library(reshape2)
library(ggplot2)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(dplyr)
library(maps)
library(sf)
library(stars)

# Open the netcdf containing the predictions you want
# Filename structure is species_m_yyyy_SDM.nc
testnc <- nc_open("C:/Users/r.wildermuth/Downloads/sard_4_2015_GAMBOOST.nc")
# Look at the file structure
print(testnc)
# Extract the spatiotemporal variables
lon <- ncvar_get(testnc, "lon")
lat <- ncvar_get(testnc, "lat")
tim <- ncvar_get(testnc, "time")
# Get the SDM predictions
predSDM <- ncvar_get(testnc, "predMBOOST")
# Close the netcdf
nc_close(testnc)

# Reshape the 3D array so we can map it, change the time field to be date
dimnames(predSDM) <- list(lon = lon, lat = lat, tim = tim)
sdmMelt <- reshape2::melt(predSDM, value.name = "predSDM")
sdmMelt$dt <- as.Date("1900-01-01") + days(sdmMelt$tim)

# # Optional (but recommended): trim predictions to within 300-500km of the coast
# if(!exists("distLand")) {
#   distLand <- (read.csv("DistLandROMSPoints.csv", head=TRUE, sep=","))[c("lon","lat","distLand")]
# }
# sdmMelt <- dplyr::full_join(sdmMelt, distLand, by = c("lon", "lat"))
# sdmMelt <- subset(sdmMelt, sdmMelt$distLand < 500000)

# Use CalCurrent Atlantis extent instead
ccAtl <- st_read("C:/Users/r.wildermuth/Documents/FutureSeas/MapFiles/emocc_whole_domain.shp")
ccAtl <- st_transform(ccAtl, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# Map prep
pac.coast <- borders("world", colour="gray50", fill="gray50", xlim = c(-140, -100), ylim = c(20, 60))
mycols <- colors()[c(473,562,71,610,655,653,621,34)]
mypalette <- colorRampPalette(mycols)(255)

# Select a date to map. This example just uses the median value
myDate <- round(median(sdmMelt$dt), 0) 
p1 <- ggplot(data = subset(sdmMelt, sdmMelt$dt == myDate)) +
  geom_tile(aes(x = lon, y = lat, fill = predSDM)) +
  scale_fill_gradientn(colours = mypalette, limits = c(0, 0.8), na.value = NA) +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  pac.coast + 
  geom_sf(data = ccAtl, color = "white", size = 1.5, fill = NA) +
  coord_sf(xlim = c(-140, -110), ylim = c(25, 53)) 
p1

# Subset to CalCurrent Atlantis extent
# from https://luisdva.github.io/rstats/GIS-with-R/
sdmRaster <- read_stars("C:/Users/r.wildermuth/Downloads/sard_4_2015_GAMBOOST.nc")
st_crs(sdmRaster) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
sdmRaster <- st_transform(sdmRaster, "+proj=longlat +ellps=WGS84 +datum=WGS84")

croppedSDM <- st_crop(x = sdmRaster,
                      y = ccAtl)
plot(croppedSDM)

ggplot() +
  geom_stars(data = croppedSDM)

sdmSlice15 <- slice(sdmRaster, index = 15, along = "time")

croppedSlice15 <- st_crop(x = sdmSlice15,
                          y = ccAtl)
plot(croppedSlice15)
sdmSlice10 <- slice(sdmRaster, index = 10, along = "time")

croppedSlice10 <- st_crop(x = sdmSlice10,
                          y = ccAtl)
plot(croppedSlice10)

ggplot() +
  geom_stars(data = croppedSlice) +
  scale_fill_gradientn(colours = mypalette, limits = c(0, 0.8), na.value = NA) +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  pac.coast + 
  geom_sf(data = ccAtl, color = "white", size = 1.5, fill = NA) +
  coord_sf(xlim = c(-140, -110), ylim = c(25, 53))

# Find upper 95% quantile
upper95 <- quantile(croppedSlice15$sard_4_2015_GAMBOOST.nc, probs = 0.95, na.rm = TRUE)
probMean <- mean(croppedSlice15$sard_4_2015_GAMBOOST.nc, na.rm = TRUE)
probSD <- sd(croppedSlice15$sard_4_2015_GAMBOOST.nc, na.rm = TRUE)
probMean + 1.96*probSD

# Try masking low probability cells
topProbSlice10 <- croppedSlice10
topProbSlice10[croppedSlice10 < upper95] <- NA

ggplot() +
  geom_stars(data = topProbSlice10) +
  scale_fill_gradientn(colours = mypalette, limits = c(0, 0.8), na.value = NA) +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  pac.coast + 
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  coord_sf(xlim = c(-140, -110), ylim = c(25, 53))

topProbSlice15 <- croppedSlice15
topProbSlice15[croppedSlice15 < upper95] <- NA

sum(topProbSlice10[[1]], na.rm = TRUE)
sum(topProbSlice15[[1]], na.rm = TRUE)

topProbSDM <- croppedSDM
topProbSDM[croppedSDM < upper95] <- NA
plot(topProbSDM)

# Code to create study area figure for DFA pub
# Created: 12/11/2024, Robert Wildermuth

library(tidyverse)
library(lubridate)
library(sf)
library(stars)
library(cubelyr)
library(gridExtra)

# Sardine spawning grounds ------------------------------------------------

# code from ProcessData.Rmd

# data file path
datPath <- "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/"

# sardine and anchovy have different breeding seasons
sardFiles <- expand.grid("sard", 3:7, 1998:2022, "mboost_roms.nc") %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, Var4, sep = "_")) %>%
  pull(sdmFile)

sardFiles <- paste0(datPath, "SDMoutput/sardine_noSSB/", sardFiles)
sardSDMs <- read_stars(sardFiles, proxy = FALSE, quiet = TRUE)
# Time fix from Barb
wrongTimes <- as.Date(st_get_dimension_values(sardSDMs, "time"))
wrongTimes[c(1, length(wrongTimes))] # Can see is wrong baseline
rightTimes <- wrongTimes + 693579 # 693579 is the number of days since 0000-00-00
sardSDMs <- st_set_dimensions(sardSDMs, 3, values = rightTimes, names = "time")
st_get_dimension_values(sardSDMs, "time")[c(1, length(wrongTimes))] # Check times look ok now

st_crs(sardSDMs) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
sardSDMs <- st_transform(sardSDMs, "+proj=longlat +ellps=WGS84 +datum=WGS84")
#st_get_dimension_values(sardSDMs, "time")

# Use CalCurrent Atlantis extent to limit sample frame
ccAtl <- st_read("C:/Users/r.wildermuth/Documents/FutureSeas/MapFiles/emocc_whole_domain.shp")
ccAtl <- st_transform(ccAtl, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# restrict to area south of 40deg N
newAtlbbox <- st_bbox(ccAtl)
newAtlbbox$ymax <- 40.0
newAtlbbox <- st_bbox(unlist(newAtlbbox))
ccAtl <- st_crop(ccAtl, newAtlbbox)

# crop to Atlantis extent
sardSDMs <- st_crop(x = sardSDMs,
                    y = ccAtl)

# show first 5 days of SDM
sardSlice <- slice(sardSDMs, index = 1:5, along = "time")
ggplot() + geom_stars(data = sardSlice) + facet_wrap(~time)

# Calculate high probability habitat

# discrimination threshold from Barb Muhling:
sardHabThresh <- 0.45

# Try masking low probability cells
sardSpawn <- sardSDMs
sardSpawn[sardSDMs < sardHabThresh] <- NA
# Don't allow spawning habitat above 40deg north
#sardSpawn2 <- sardSpawn %>% filter(y < 40)

# Map prep
pac.coast <- borders("world", colour="grey", fill="grey", xlim = c(-140, -100), ylim = c(20, 60))
mycols <- RColorBrewer::brewer.pal(9, "Greens")#colors()[c(473,562,71,610,655,653,621,34)]
mypalette <- colorRampPalette(mycols)(255)

ggplot() +
  geom_stars(data = sardSpawnSlice) +
  scale_fill_gradientn(colours = mypalette, limits = c(0, 0.8), na.value = NA) +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  pac.coast + 
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  coord_sf(xlim = c(-130, -113), ylim = c(27, 42)) +
  facet_wrap(~time)

sardSpawnYr <- aggregate(sardSpawn, by = "years", FUN = sum, na.rm = TRUE)
st_get_dimension_values(sardSpawnYr, "time")

# reorder and then make 0s NAs to turn them grey
sardSpawnYr <- aperm(sardSpawnYr, c(2,3,1))
sardSpawnYr[sardSpawnYr == 0] <- NA

# aggregate over time
by_t <- as.Date(c("1998-01-01", "2022-01-01"))
sardSpawnGrounds <- aggregate(sardSpawnYr, by = by_t, FUN = sum, 
                              na.rm = TRUE, rightmost.closed = TRUE)
# rearrange the dimensions and get only the aggregated time slice
sardSpawnGrounds <- aperm(sardSpawnGrounds, c(2,3,1))
sardSpawnGrounds <- sardSpawnGrounds[, , , 1]

sardSpawnGroundsCont <- sardSpawnGrounds
sardSpawnGroundsCont[sardSpawnGrounds == 0] <- NA

# Get EEZ
eezs <- st_read("C:/Users/r.wildermuth/Documents/FutureSeas/MapFiles/World_EEZ_v11_20191118/eez_boundaries_v11.shp")
eezs <- st_transform(eezs, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# Get CalCOFI core area
coreCalCOFI <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/MapFiles/CalCOFI core stations.csv") %>% 
                rename(lat = `Lat (dec)`,
                       lon = `Lon (dec)`)
coreCalCOFI <- coreCalCOFI %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull()
st_crs(coreCalCOFI) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# NEMURO ZM average areas
nemuroZM <- data.frame(lat = c(31, 34.5, 34.5, 31, 40.5, 44.5, 44.5, 40.5),
                       lon = c(-116, -116, -124, -124, -124, -124, -130, -130),
                       NorS = c(rep("SoCal", 4), rep("NorCal", 4)))
nemuroZM <- nemuroZM %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  group_by(NorS) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull()
st_crs(nemuroZM) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# Reference lines for upwelling and poleward transport
transpLines <- data.frame(lat = c(44.6, 44.6,
                                  32, 32, 32,
                                  33, 33,
                                  36, 36,
                                  39, 39),
                          lon = c(-124, -125,
                                  -116, -120, -124,
                                  -118, -130,
                                  -122, -130,
                                  -124, -130),
                          ID = c(rep("Newport", 2),
                                 rep("Poleward", 3),
                                 rep("upwell33", 2),
                                 rep("upwell36", 2),
                                 rep("upwell39", 2)))
transpLines <- transpLines %>% 
                  st_as_sf(coords = c("lon", "lat"), 
                           crs = "+proj=longlat +ellps=WGS84 +datum=WGS84") %>% 
                  group_by(ID) %>% 
                  summarize() %>%
                  st_cast("MULTILINESTRING")

sardPlot <- ggplot() +
  geom_stars(data = sardSpawnGroundsCont) +
  pac.coast +
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  geom_sf(data = eezs, color = "steelblue", size = 1.5, fill = NA) +
  geom_sf(data = coreCalCOFI, color = "orangered", size = 1.5, fill = NA) +
  geom_sf(data = nemuroZM, color = "darkgreen", size = 1.5, fill = NA) +
  geom_sf(data = transpLines, color = "purple", size = 3) +
  scale_fill_viridis_c(option = "mako",
                       direction = -1,
                       na.value = NA) +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  
  coord_sf(xlim = c(-130, -115), ylim = c(28, 48)) +
  facet_wrap(~time) +
  theme_minimal() +
  labs(fill = "Cumulative \nsardine habitat")

# Change to categorical (1 = in spawning grounds)
sardSpawnGrounds[sardSpawnGrounds > 0] <- 1

ggplot() +
  geom_stars(data = sardSpawnGrounds) +
  pac.coast +
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  # 
  coord_sf(xlim = c(-130, -115), ylim = c(28, 42)) +
  facet_wrap(~time)

# Anchovy spawning grounds -----------------------------------------------------------------

anchFiles <- expand.grid("anch", 2:4, 1998:2022, "mboost_roms.nc") %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, Var4, sep = "_")) %>%
  pull(sdmFile)

anchFiles <- paste0(datPath, "SDMoutput/anchovy_noSSB/", anchFiles)
anchSDMs <- read_stars(anchFiles, proxy = FALSE, quiet = TRUE)
# Time fix from Barb
wrongTimes <- as.Date(st_get_dimension_values(anchSDMs, "time"))
wrongTimes[c(1, length(wrongTimes))] # Can see is wrong baseline
rightTimes <- wrongTimes + 693579 # 693579 is the number of days since 0000-00-00
anchSDMs <- st_set_dimensions(anchSDMs, 3, values = rightTimes, names = "time")
st_get_dimension_values(anchSDMs, "time")[c(1, length(wrongTimes))] # Check times look ok now

st_crs(anchSDMs) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
anchSDMs <- st_transform(anchSDMs, "+proj=longlat +ellps=WGS84 +datum=WGS84")
#st_get_dimension_values(anchSDMs, "time")

# crop to Atlantis extent
anchSDMs <- st_crop(x = anchSDMs,
                    y = ccAtl)

# Calculate high probability habitat
# discrimination threshold from Barb Muhling:
anchHabThresh <- 0.29

# Try masking low probability cells
anchSpawn <- anchSDMs
anchSpawn[anchSDMs < anchHabThresh] <- NA

anchSpawnYr <- aggregate(anchSpawn, by = "years", FUN = sum, na.rm = TRUE)
st_get_dimension_values(anchSpawnYr, "time")

anchSpawnGrounds <- aggregate(anchSpawnYr, by = by_t, FUN = sum, 
                              na.rm = TRUE, rightmost.closed = TRUE)
# Change to categorical (1 = in spawning grounds)
# anchSpawnGrounds[anchSpawnGrounds > 0] <- 1
# rearrange the dimensions and get only the aggregated time slice
anchSpawnGrounds <- aperm(anchSpawnGrounds, c(2,3,1))
anchSpawnGrounds <- anchSpawnGrounds[, , , 1]

anchSpawnGroundsCont <- anchSpawnGrounds
anchSpawnGroundsCont[anchSpawnGrounds == 0] <- NA

anchPlot <- ggplot() +
  geom_stars(data = anchSpawnGroundsCont) +
  pac.coast +
  geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
  geom_sf(data = eezs, color = "steelblue", size = 1.5, fill = NA) +
  geom_sf(data = coreCalCOFI, color = "orangered", size = 1.5, fill = NA) +
  geom_sf(data = nemuroZM, color = "darkgreen", size = 1.5, fill = NA) +
  geom_sf(data = transpLines, color = "purple", size = 3) +
  scale_fill_viridis_c(option = "mako",
                       direction = -1,
                       na.value = NA) +
  guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
  
  coord_sf(xlim = c(-130, -115), ylim = c(28, 48)) +
  facet_wrap(~time) +
  theme_minimal()+
  labs(fill = "Cumulative \nanchovy habitat")

grid.arrange(sardPlot, anchPlot)

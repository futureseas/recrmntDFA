# Projection data

library(tidyverse)
library(lubridate)
library(sf)
library(stars)
library(cubelyr)
library(units)
library(tidync)

# data file path
datPath <- "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/"


# SDM-related indicators --------------------------------------------------

# Species distribution model outputs
# Spawning habitat area and duration are calculated from ROMS-based species 
# distribution models produced by Barb Muhling and downloaded 3/29/2023. Habitat 
# predictions are derived from SDMs without SSB included as a predictor. For now 
# using only GAM estimates of SDM.

# From conversation with Barb on 3/29/23:
#   - ROMS-based estimates have more predictors, thus better predictions of habitat
# - ROMS sample frame likely represents northern Pacific sardine stock better
# - ROMS predictions need updates before final version of analysis

## Spawning habitat area
# A proxy of area of spawning habitat is determined by compiling daily habitat 
# estimates across years within the traditional breeding season (for sardine: 
# March - July). Breeding habitat is defined as SDM cells with probabilities of 
# occurrence at or above the upper 95% confidence interval. The annual index is 
# calculated as the sum of probabilities of occurrence in breeding habitat over 
# all cells and days within the traditional breeding season for each year.


## Fxn to process annual habitat values
CalcAnnualHab <- function(sdmFiles, # character vector of filepaths to SDM output
                          threshold, # discrimination threshold from the SDM
                          shp, # shapefile of study are to clip map to
                          dailyHab = FALSE # logical whether to include daily high quality habitat stars object
){
  sardSDMs <- read_stars(sdmFiles, proxy = FALSE, quiet = TRUE)
  
  # Time fix from Barb
  wrongTimes <- as.Date(st_get_dimension_values(sardSDMs, "time"))
  # wrongTimes[c(1, length(wrongTimes))] # Can see is wrong baseline
  rightTimes <- wrongTimes + 693579 # 693579 is the number of days since 0000-00-00
  sardSDMs <- st_set_dimensions(sardSDMs, 3, values = rightTimes, names = "time")
  # st_get_dimension_values(sardSDMs, "time")[c(1, length(wrongTimes))] # Check times look ok now
  
  st_crs(sardSDMs) <- "+proj=longlat +ellps=WGS84 +datum=WGS84"
  sardSDMs <- st_transform(sardSDMs, "+proj=longlat +ellps=WGS84 +datum=WGS84")
  #st_get_dimension_values(sardSDMs, "time")
  
  # crop to Atlantis extent
  sardSDMs <- st_crop(x = sardSDMs,
                      y = shp)
  
  # Try masking low probability cells
  sardSpawn <- sardSDMs
  sardSpawn[sardSDMs < threshold] <- NA
  
  sardSpawnYr <- aggregate(sardSpawn, by = "years", FUN = sum, na.rm = TRUE)
  # st_get_dimension_values(sardSpawnYr, "time")
  
  # reorder and then make 0s NAs to turn them grey
  sardSpawnYr <- aperm(sardSpawnYr, c(2,3,1))
  sardSpawnYr[sardSpawnYr == 0] <- NA
  
  sardSpawnHab <- aggregate(sardSpawnYr, shp, FUN = sum, na.rm = TRUE)
  
  if(dailyHab){
    return(list(sardSpawnHab, sardSpawn))
  } else {
    return(list(sardSpawnHab))
  }
  
}

# Use CalCurrent Atlantis extent to limit sample frame
ccAtl <- st_read("C:/Users/r.wildermuth/Documents/FutureSeas/MapFiles/emocc_whole_domain.shp")
ccAtl <- st_transform(ccAtl, "+proj=longlat +ellps=WGS84 +datum=WGS84")

# restrict to area south of 40deg N
newAtlbbox <- st_bbox(ccAtl)
newAtlbbox$ymax <- 40.0
newAtlbbox <- st_bbox(unlist(newAtlbbox))
ccAtl <- st_crop(ccAtl, newAtlbbox)

## Nursury habitat area

# discrimination threshold from Barb Muhling:
sardLarvHabThresh <- 0.17

# sardine and anchovy have different breeding seasons
sardFiles19802040 <- expand.grid("sardLarvaeSDMs_proj", 3:7, 1980:2040) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

sardFiles19802040 <- paste0(datPath, "SDMoutput/proj_sardLarv/", sardFiles19802040, ".nc")

sardLrvHab19802040 <- CalcAnnualHab(sdmFiles = sardFiles19802040, threshold = sardLarvHabThresh, 
                                    shp = ccAtl)

sardFiles20412100 <- expand.grid("sardLarvaeSDMs_proj", 3:7, 2041:2100) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

sardFiles20412100 <- paste0(datPath, "SDMoutput/proj_sardLarv/", sardFiles20412100, ".nc")

sardLrvHab20412100 <- CalcAnnualHab(sdmFiles = sardFiles20412100, threshold = sardLarvHabThresh, 
                                    shp = ccAtl)


#  Anchovy 

# discrimination threshold from Barb Muhling:
anchLarvHabThresh <- 0.42

anchFiles19802040 <- expand.grid("anchLarvaeSDMs_proj", 2:4, 1980:2040) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

anchFiles19802040 <- paste0(datPath, "SDMoutput/proj_anchLarv/", anchFiles19802040, ".nc")

anchLrvHab19802040 <- CalcAnnualHab(sdmFiles = anchFiles19802040, threshold = anchLarvHabThresh, 
                                    shp = ccAtl)

anchFiles20412100 <- expand.grid("anchLarvaeSDMs_proj", 2:4, 2041:2100) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

anchFiles20412100 <- paste0(datPath, "SDMoutput/proj_anchLarv/", anchFiles20412100, ".nc")

anchLrvHab20412100 <- CalcAnnualHab(sdmFiles = anchFiles20412100, threshold = anchLarvHabThresh, 
                                    shp = ccAtl)

nurseTS <- rbind(data.frame(year = 1980:2040,
                      sardNurseHab_GFDL = t(sardLrvHab19802040[[1]]$predGAMBOOST_GFDL),
                      sardNurseHab_IPSL = t(sardLrvHab19802040[[1]]$predGAMBOOST_IPSL),
                      sardNurseHab_HAD = t(sardLrvHab19802040[[1]]$predGAMBOOST_HADL),
                      anchNurseHab_GFDL = t(anchLrvHab19802040[[1]]$predGAMBOOST_GFDL),
                      anchNurseHab_IPSL = t(anchLrvHab19802040[[1]]$predGAMBOOST_IPSL),
                      anchNurseHab_HAD = t(anchLrvHab19802040[[1]]$predGAMBOOST_HADL)),
                data.frame(year = 2041:2100,
                           sardNurseHab_GFDL = t(sardLrvHab20412100[[1]]$predGAMBOOST_GFDL),
                           sardNurseHab_IPSL = t(sardLrvHab20412100[[1]]$predGAMBOOST_IPSL),
                           sardNurseHab_HAD = t(sardLrvHab20412100[[1]]$predGAMBOOST_HADL),
                           anchNurseHab_GFDL = t(anchLrvHab20412100[[1]]$predGAMBOOST_GFDL),
                           anchNurseHab_IPSL = t(anchLrvHab20412100[[1]]$predGAMBOOST_IPSL),
                           anchNurseHab_HAD = t(anchLrvHab20412100[[1]]$predGAMBOOST_HADL)))

nurseTS %>% pivot_longer(cols = -year, names_to = c("spp", "ESM"), names_sep = "_") %>%
  ggplot() +
  geom_line(aes(x = year, y = value, col = ESM)) +
  facet_wrap(~spp)


## Spawning Habitat

# discrimination threshold from Barb Muhling:
sardHabThresh <- 0.45

# sardine and anchovy have different breeding seasons
sardFiles19802040 <- expand.grid("sardSDMs_proj", 3:7, 1980:2040) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

sardFiles19802040 <- paste0(datPath, "SDMoutput/proj_sardine/", sardFiles19802040, ".nc")

sardSpHab19802040 <- CalcAnnualHab(sdmFiles = sardFiles19802040, threshold = sardHabThresh, 
                           shp = ccAtl, dailyHab = TRUE)

sardFiles20412100 <- expand.grid("sardSDMs_proj", 3:7, 2041:2100) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

sardFiles20412100 <- paste0(datPath, "SDMoutput/proj_sardine/", sardFiles20412100, ".nc")

sardSpHab20412100 <- CalcAnnualHab(sdmFiles = sardFiles20412100, threshold = sardHabThresh, 
                           shp = ccAtl, dailyHab = TRUE)
# check space use:
sort( sapply(ls(),function(x){object.size(get(x))})) 

# # show first 5 days of SDM
# sardSlice <- slice(sardSDMs, index = 1:5, along = "time")
# ggplot() + geom_stars(data = sardSlice) + facet_wrap(~time)
# 
# # Calculate high probability habitat
# # probMean <- mean(sardSDMs$sard_3_1998_mboost_roms.nc, na.rm = TRUE)
# # probSD <- sd(sardSDMs$sard_3_1998_mboost_roms.nc, na.rm = TRUE)
# # upper95 <- probMean + 1.96*probSD
# 
# 
# 
# # Don't allow spawning habitat above 40deg north
# #sardSpawn2 <- sardSpawn %>% filter(y < 40)
# 
# sardSpawnSlice <- slice(sardSpawn, index = 15:25, along = "time")
# # ggplot() + geom_stars(data = sardSpawnSlice) + facet_wrap(~time)
# 
# # Map prep
# pac.coast <- borders("world", colour="brown", fill="brown", xlim = c(-140, -100), ylim = c(20, 60))
# mycols <- RColorBrewer::brewer.pal(9, "Greens")#colors()[c(473,562,71,610,655,653,621,34)]
# mypalette <- colorRampPalette(mycols)(255)
# 
# ggplot() +
#   geom_stars(data = sardSpawnSlice) +
#   scale_fill_gradientn(colours = mypalette, limits = c(0, 0.8), na.value = NA) +
#   guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
#   pac.coast + 
#   geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
#   coord_sf(xlim = c(-130, -113), ylim = c(27, 42)) +
#   facet_wrap(~time)
# 
# 
# ggplot() +
#   geom_stars(data = sardSpawnYr["predGAMBOOST_GFDL"]) +
#   pac.coast +
#   geom_sf(data = ccAtl, color = "black", size = 1.5, fill = NA) +
#   scale_fill_gradientn(colours = mycols[3:9],
#                        limits = c(0, max(unlist(sardSpawnYr), na.rm = TRUE)),
#                        na.value = NA) +
#   guides(fill = guide_colorbar(barwidth=0.5, barheight=5)) +
#   
#   coord_sf(xlim = c(-130, -115), ylim = c(28, 42)) +
#   facet_wrap(~time) +
#   theme_minimal()
# 

#  Anchovy 

# discrimination threshold from Barb Muhling:
anchHabThresh <- 0.29

anchFiles19802040 <- expand.grid("anchSDMs_proj", 2:4, 1980:2040) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

anchFiles19802040 <- paste0(datPath, "SDMoutput/proj_anchovy/", anchFiles19802040, ".nc")

anchSpHab19802040 <- CalcAnnualHab(sdmFiles = anchFiles19802040, threshold = anchHabThresh, 
                           shp = ccAtl, dailyHab = TRUE)

anchFiles20412100 <- expand.grid("anchSDMs_proj", 2:4, 2041:2100) %>% 
  mutate(sdmFile = paste(Var1, Var2, Var3, sep = "_")) %>%
  pull(sdmFile)

anchFiles20412100 <- paste0(datPath, "SDMoutput/proj_anchovy/", anchFiles20412100, ".nc")

anchSpHab20412100 <- CalcAnnualHab(sdmFiles = anchFiles20412100, threshold = anchHabThresh, 
                           shp = ccAtl, dailyHab = TRUE)

spawnTS <- rbind(data.frame(year = 1980:2040,
                            sardSpawnHab_GFDL = t(sardSpHab19802040[[1]]$predGAMBOOST_GFDL),
                            sardSpawnHab_IPSL = t(sardSpHab19802040[[1]]$predGAMBOOST_IPSL),
                            sardSpawnHab_HAD = t(sardSpHab19802040[[1]]$predGAMBOOST_HADL),
                            anchSpawnHab_GFDL = t(anchSpHab19802040[[1]]$predGAMBOOST_GFDL),
                            anchSpawnHab_IPSL = t(anchSpHab19802040[[1]]$predGAMBOOST_IPSL),
                            anchSpawnHab_HAD = t(anchSpHab19802040[[1]]$predGAMBOOST_HADL)),
                 data.frame(year = 2041:2100,
                            sardSpawnHab_GFDL = t(sardSpHab20412100[[1]]$predGAMBOOST_GFDL),
                            sardSpawnHab_IPSL = t(sardSpHab20412100[[1]]$predGAMBOOST_IPSL),
                            sardSpawnHab_HAD = t(sardSpHab20412100[[1]]$predGAMBOOST_HADL),
                            anchSpawnHab_GFDL = t(anchSpHab20412100[[1]]$predGAMBOOST_GFDL),
                            anchSpawnHab_IPSL = t(anchSpHab20412100[[1]]$predGAMBOOST_IPSL),
                            anchSpawnHab_HAD = t(anchSpHab20412100[[1]]$predGAMBOOST_HADL)))

spawnTS %>% pivot_longer(cols = -year, names_to = c("spp", "ESM"), names_sep = "_") %>%
  ggplot() +
  geom_line(aes(x = year, y = value, col = ESM)) +
  facet_wrap(~spp)




## Spawning season timing and duration
load(file = paste0(datPath, "histSardSpawnGr.RData"))
load(file = paste0(datPath, "histAnchSpawnGr.RData"))

#Find proportion of breeding habitat in spawning ground
# convert spawning days to categorical
sardSpawn1NA19802040 <- sardSpHab19802040[[2]]
# find number of high probability cells at each time
# sardAggSpwn19802040 <- aggregate(sardSpawn1NA19802040, by = ccAtl, FUN = sum, na.rm = TRUE)
# sardAggSpwn <- t(sardAggSpwn$sard_3_1998_mboost_roms.nc)
sardSpawn1NA19802040[sardSpHab19802040[[2]] > sardHabThresh] <- 1 # should only be NAs and 1s

# find number of cells in spawning grounds
sardSpGrdDenom <- aggregate(sardSpawnGrounds, by = ccAtl, FUN = sum, na.rm = TRUE)
# sardSpGrdDenom <- as.numeric(sardSpGrdDenom$sard_3_1998_mboost_roms.nc)

# find number of high probability cells at each time
sardSpGrdNum19802040 <- aggregate(sardSpawn1NA19802040, by = ccAtl, FUN = sum, na.rm = TRUE)
# sardSpGrdNum <- t(sardSpGrdNum$sard_3_1998_mboost_roms.nc) 


# convert spawning days to categorical
sardSpawn1NA20412100 <- sardSpHab20412100[[2]]
# find number of high probability cells at each time
# sardAggSpwn20412100 <- aggregate(sardSpawn1NA20412100, by = ccAtl, FUN = sum, na.rm = TRUE)
# sardAggSpwn <- t(sardAggSpwn$sard_3_1998_mboost_roms.nc)
sardSpawn1NA20412100[sardSpHab20412100[[2]] > sardHabThresh] <- 1 # should only be NAs and 1s

# find number of high probability cells at each time
sardSpGrdNum20412100 <- aggregate(sardSpawn1NA20412100, by = ccAtl, FUN = sum, na.rm = TRUE)

sardSpDates <- rbind(data.frame(date = st_get_dimension_values(sardSpHab19802040[[2]], which = "time"),
                          # propSpawn = sardSpGrdNum/sardSpGrdDenom,
                          # aggSpawn = sardAggSpwn,
                          propSpawn_GFDL = t(sardSpGrdNum19802040$predGAMBOOST_GFDL),
                          propSpawn_IPSL = t(sardSpGrdNum19802040$predGAMBOOST_IPSL),
                          propSpawn_HAD = t(sardSpGrdNum19802040$predGAMBOOST_HADL)),
                data.frame(date = st_get_dimension_values(sardSpHab20412100[[2]], which = "time"),
                          # propSpawn = sardSpGrdNum/sardSpGrdDenom,
                          # aggSpawn = sardAggSpwn,
                          propSpawn_GFDL = t(sardSpGrdNum20412100$predGAMBOOST_GFDL),
                          propSpawn_IPSL = t(sardSpGrdNum20412100$predGAMBOOST_IPSL),
                          propSpawn_HAD = t(sardSpGrdNum20412100$predGAMBOOST_HADL))) %>%
  mutate(year = year(ymd(date)),
         dayofYr = yday(date),
         propSpawn_GFDL = propSpawn_GFDL/as.vector(sardSpGrdDenom$sard_3_1998_mboost_roms.nc),
         propSpawn_IPSL = propSpawn_IPSL/as.vector(sardSpGrdDenom$sard_3_1998_mboost_roms.nc),
         propSpawn_HAD = propSpawn_HAD/as.vector(sardSpGrdDenom$sard_3_1998_mboost_roms.nc))

# summary(sardSpGrdNum/sardSpGrdDenom)

# sardSpDates %>% #filter(year == 2015) %>%
#   ggplot(aes(x = dayofYr, y = propSpawn)) +
#   geom_line() +
#   geom_hline(yintercept = quantile(sardSpGrdNum/sardSpGrdDenom, probs = c(0.1, 0.2, 0.5))) +
#   facet_wrap(~year, scales = "free_y")
# 
# sardSpDates %>% #filter(year == 2015) %>%
#   ggplot(aes(x = dayofYr, y = aggSpawn)) +
#   geom_line() +
#   geom_hline(yintercept = quantile(sardAggSpwn, probs = c(0.1, 0.2, 0.5))) +
#   facet_wrap(~year, scales = "free_y")

# Calculate number of days 5% of spawning grounds contain high quality habitat
sardSpawnDays <- sardSpDates %>% pivot_longer(c(propSpawn_GFDL, 
                                                propSpawn_IPSL, 
                                                propSpawn_HAD), 
                                              names_to = "ESM",
                                              values_to = "propSpawn") %>%
  filter(propSpawn >= 0.05) %>% 
  group_by(year, ESM) %>%
  summarize(daysAbove5pct = n()) %>% 
  mutate(ESM = sub("propSpawn_", "", ESM))

# convert spawning days to categorical
anchSpawn1NA19802040 <- anchSpHab19802040[[2]]
anchSpawn1NA19802040[anchSpHab19802040[[2]] > anchHabThresh] <- 1 # should only be NAs and 1s

# find number of cells in spawning grounds
anchSpGrdDenom <- aggregate(anchSpawnGrounds, by = ccAtl, FUN = sum, na.rm = TRUE)
# anchSpGrdDenom <- as.numeric(anchSpGrdDenom$anch_2_1998_mboost_roms.nc)

# find number of high probability cells at each time
anchSpGrdNum19802040 <- aggregate(anchSpawn1NA19802040, by = ccAtl, FUN = sum, na.rm = TRUE)

anchSpawn1NA20412100 <- anchSpHab20412100[[2]]
anchSpawn1NA20412100[anchSpHab20412100[[2]] > anchHabThresh] <- 1 # should only be NAs and 1s

# find number of high probability cells at each time
anchSpGrdNum20412100 <- aggregate(anchSpawn1NA20412100, by = ccAtl, FUN = sum, na.rm = TRUE)

anchSpDates <- rbind(data.frame(date = st_get_dimension_values(anchSpHab19802040[[2]], which = "time"),
                                propSpawn_GFDL = t(anchSpGrdNum19802040$predGAMBOOST_GFDL),
                                propSpawn_IPSL = t(anchSpGrdNum19802040$predGAMBOOST_IPSL),
                                propSpawn_HAD = t(anchSpGrdNum19802040$predGAMBOOST_HADL)),
                     data.frame(date = st_get_dimension_values(anchSpHab20412100[[2]], which = "time"),
                                propSpawn_GFDL = t(anchSpGrdNum20412100$predGAMBOOST_GFDL),
                                propSpawn_IPSL = t(anchSpGrdNum20412100$predGAMBOOST_IPSL),
                                propSpawn_HAD = t(anchSpGrdNum20412100$predGAMBOOST_HADL))) %>%
  mutate(year = year(ymd(date)),
         dayofYr = yday(date),
         propSpawn_GFDL = propSpawn_GFDL/as.vector(anchSpGrdDenom$anch_2_1998_mboost_roms.nc),
         propSpawn_IPSL = propSpawn_IPSL/as.vector(anchSpGrdDenom$anch_2_1998_mboost_roms.nc),
         propSpawn_HAD = propSpawn_HAD/as.vector(anchSpGrdDenom$anch_2_1998_mboost_roms.nc))

# summary(anchSpGrdNum/anchSpGrdDenom)
# 
# anchSpDates %>% #filter(year == 2015) %>%
#   ggplot(aes(x = dayofYr, y = propSpawn)) +
#   geom_line() +
#   geom_hline(yintercept = quantile(anchSpGrdNum/anchSpGrdDenom, probs = c(0.1, 0.2, 0.5))) +
#   facet_wrap(~year, scales = "free_y")

# Calculate number of days 40% of spawning grounds contain high quality habitat
anchSpawnDays <- anchSpDates %>% pivot_longer(c(propSpawn_GFDL, 
                                                propSpawn_IPSL, 
                                                propSpawn_HAD), 
                                              names_to = "ESM",
                                              values_to = "propSpawn") %>% 
  filter(propSpawn >= 0.4) %>% 
  group_by(year, ESM) %>%
  summarize(daysAbove40pct = n()) %>% 
  mutate(ESM = sub("propSpawn_", "", ESM))


# SST ---------------------------------------------------------------------

# From Mer Pozo Buil (6/4/2023): "monthly SST averaged for the diagonal CalCOFI region for the 3 datasets:
# ROMS-ESMS: includes the 3 projections (ROMS-GFDL, ROMS-IPSL, and ROMS-HAD) from 1980-2100
# WCRNT: includes ROMS real-near time from 2011-Oct2022
# WCRA: includes ROMS reanalysis from 1980-2010"
sstProj <- read_csv(paste0(datPath, "ROMS-ESMS_CalCOFI_monave_SST_198001-210012.csv"))

sstProj <- sstProj %>% mutate(seas = case_when(month %in% 1:6 ~ "springSST",
                                               month %in% 7:12 ~ "summerSST")) %>%
  group_by(year, seas) %>%
  summarize(sstGFDL = mean(sst_roms_gfdl),
            sstIPSL = mean(sst_roms_ipsl),
            sstHAD = mean(sst_roms_had)) %>% 
  pivot_longer(cols = c(sstGFDL, sstIPSL, sstHAD),
               names_to = "ESM",
               names_prefix = "sst",
               values_to = "SST")

sstProj %>% ggplot(aes(x = year, y = SST, color = ESM)) +
  geom_line() +
  facet_wrap(~seas, ncol = 1)+
  theme_classic()

sstProj <- sstProj %>% pivot_wider(names_from = seas, values_from = SST)

# Poleward transport ------------------------------------------------------

# Northward transport across 32N
# Downloaded from Mike Jacox (8/18/2023)
offshoreGFDL <- read_ncdf(paste0(datPath, "monthly_northward_velocity_gfdl_0-50m_124.0-120.0W_31.9-32.1N_1980-2100.nc"))
offshoreIPSL <- read_ncdf(paste0(datPath, "monthly_northward_velocity_ipsl_0-50m_124.0-120.0W_31.9-32.1N_1980-2100.nc"))
offshoreHAD <- read_ncdf(paste0(datPath, "monthly_northward_velocity_had_0-50m_124.0-120.0W_31.9-32.1N_1980-2100.nc"))

offshoreTrans <- rbind(data.frame(year = offshoreGFDL$year,
                                  month = offshoreGFDL$month,
                                  offshTrans = offshoreGFDL$v,
                                  ESM = "GFDL"),
                       data.frame(year = offshoreIPSL$year,
                                  month = offshoreIPSL$month,
                                  offshTrans = offshoreIPSL$v,
                                  ESM = "IPSL"),
                       data.frame(year = offshoreHAD$year,
                                  month = offshoreHAD$month,
                                  offshTrans = offshoreHAD$v,
                                  ESM = "HAD"))

transport <- offshoreTrans %>% mutate(seas = case_when(month %in% 1:6 ~ "spring",
                                                       month %in% 6:12 ~ "summer"),
                                      season = case_when(month %in% 3:5 ~ "spring",
                                                         month %in% 6:8 ~ "summer",
                                                         month %in% 9:11 ~ "fall",
                                                         month %in% c(12,1,2) ~ "winter")) %>%
  group_by(year, seas, ESM) %>%
  summarize(avgOffshTransp = mean(offshTrans)) #%>% 
# ggplot(aes(x = year, y = avgOffshTran, color = season)) +
# geom_line()

# pivot_wider(names_from = seas, names_prefix = "avgOffTrans",
#             values_from = avgOffshTran)

nearshoreGFDL <- read_ncdf(paste0(datPath, "monthly_northward_velocity_gfdl_0-50m_120.0-115.0W_31.9-32.1N_1980-2100.nc"))
nearshoreIPSL <- read_ncdf(paste0(datPath, "monthly_northward_velocity_ipsl_0-50m_120.0-115.0W_31.9-32.1N_1980-2100.nc"))
nearshoreHAD <- read_ncdf(paste0(datPath, "monthly_northward_velocity_had_0-50m_120.0-115.0W_31.9-32.1N_1980-2100.nc"))

nearshoreTrans <- rbind(data.frame(year = nearshoreGFDL$year,
                                   month = nearshoreGFDL$month,
                                   nearshTrans = nearshoreGFDL$v,
                                   ESM = "GFDL"),
                        data.frame(year = nearshoreIPSL$year,
                                   month = nearshoreIPSL$month,
                                   nearshTrans = nearshoreIPSL$v,
                                   ESM = "IPSL"),
                        data.frame(year = nearshoreHAD$year,
                                   month = nearshoreHAD$month,
                                   nearshTrans = nearshoreHAD$v,
                                   ESM = "HAD"))


transport <- nearshoreTrans %>% mutate(seas = case_when(month %in% 1:6 ~ "spring",
                                                        month %in% 6:12 ~ "summer"),
                                       season = case_when(month %in% 3:5 ~ "spring",
                                                          month %in% 6:8 ~ "summer",
                                                          month %in% 9:11 ~ "fall",
                                                          month %in% c(12,1,2) ~ "winter")) %>%
  group_by(year, seas, ESM) %>%
  summarize(avgNearshTransp = mean(nearshTrans)) %>%
  full_join(y = transport, by = c("year", "seas", "ESM"))

transport %>% 
  ggplot(aes(x = year, y = avgNearshTransp, color = ESM)) +
  geom_line() + 
  facet_grid(rows = vars(ESM), cols = vars(seas))+
  theme_classic()

transport %>% 
  ggplot(aes(x = year, y = avgOffshTransp, color = ESM)) +
  geom_line() + 
  facet_grid(rows = vars(ESM), cols = vars(seas))+
  theme_classic()

transport <- transport %>%  pivot_wider(names_from = seas, #names_prefix = "avgNearTrans",
                                        values_from = c(avgNearshTransp, avgOffshTransp)) #%>%
# full_join(y = transport, by = "year")

# NEMURO Zooplk'n ---------------------------------------------------------

# NEMURO plankton
nemuroZ <- read_csv(paste0(datPath, "ROMS-NEMURO spring plankton.csv"))

nemuroZ %>% filter(GCM != "HIST") %>%
  ggplot(aes(x = year, y = ZM, color = GCM)) +
  geom_line() +
  facet_wrap(~zone, ncol = 1)+
  theme_classic()

nemuroZ <- nemuroZ %>% filter(GCM != "HIST") %>%           
  pivot_wider(names_from = zone, values_from = c(PS, PL, ZS, ZM, ZL)) %>%
  rename(ESM = GCM)


# ROMS upwelling indicators -----------------------------------------------
# library(ncdf4)

# functionalize getting stuff from an .nc file
GetUpwelling.nc <- function(path,
                            #variable,
                            esm
){
  # uses library(ncdf4)
  # ncBEUTI <- nc_open(path)
  # # Extract the spatiotemporal variables
  # lat <- ncvar_get(ncBEUTI, "latitude")
  # tim <- ncvar_get(ncBEUTI, "time")
  # # Get the SDM predictions
  # beutiGFDL <- ncvar_get(ncBEUTI, variable)
  # # Close the netcdf
  # nc_close(ncBEUTI)
  # 
  # # Reshape the 3D array so we can map it, change the time field to be date
  # dimnames(beutiGFDL) <- list(lat = lat, time = tim)
  # beutiGFDL <- reshape2::melt(beutiGFDL, value.name = variable)
  # beutiGFDL$dt <- as.Date("1970-01-01") + days(beutiGFDL$time)
  # beutiGFDL <- beutiGFDL %>% mutate(year = year(dt),
  #                                   month = month(dt),
  #                                   ESM = esm)
  
  # Use library(tidync)
  ncObj <- tidync(path)
  varTbl <- ncObj %>% hyper_tibble()
  timeTbl <- ncObj %>% activate("D0") %>% hyper_tibble()
  
  varTbl <- varTbl %>% full_join(y = timeTbl, by = "time") %>%
              mutate(ESM = esm)
  
  return(varTbl)
}

beutiGFDL <- GetUpwelling.nc(path = paste0(datPath, "beuti_gfdl_monthly_1980-2100.nc"),
                             # variable = "BEUTI", 
                             esm = "GFDL")
# beutiGFDL <- tidync(paste0(datPath, "beuti_gfdl_monthly_1980-2100.nc"))
# beutiGFDL %>% hyper_tibble()
# beutiGFDL %>% activate("D0") %>% hyper_tibble()

beutiHAD <- GetUpwelling.nc(path = paste0(datPath, "beuti_had_monthly_1980-2100.nc"),
                             # variable = "BEUTI", 
                            esm = "HAD")
beutiIPSL <- GetUpwelling.nc(path = paste0(datPath, "beuti_ipsl_monthly_1980-2100.nc"),
                             # variable = "BEUTI", 
                             esm = "IPSL")
  
upwelling <- rbind(beutiGFDL, beutiHAD, beutiIPSL)

cutiGFDL <- GetUpwelling.nc(path = paste0(datPath, "cuti_gfdl_monthly_1980-2100.nc"),
                             # variable = "CUTI", 
                            esm = "GFDL")
cutiHAD <- GetUpwelling.nc(path = paste0(datPath, "cuti_had_monthly_1980-2100.nc"),
                           # variable = "CUTI", 
                           esm = "HAD")
cutiIPSL <- GetUpwelling.nc(path = paste0(datPath, "cuti_ipsl_monthly_1980-2100.nc"),
                            # variable = "CUTI", 
                            esm = "IPSL")

upwelling <- upwelling %>% full_join(y = rbind(cutiGFDL, cutiHAD, cutiIPSL),
                                     by = c("latitude", "time", "year", "month", "ESM"))

# STI and LUSI
# ncPhen <- nc_open(paste0(datPath, "cuti_phenology_gfdl_1980-2100.nc"))
# # Extract the spatiotemporal variables
# lat <- ncvar_get(ncPhen, "latitude")
# tim <- ncvar_get(ncPhen, "time")
# # Get the SDM predictions
# stiGFDL <- ncvar_get(ncPhen, "STI")
# lusiGFDL <- ncvar_get(ncPhen, "LUSI")
# # Close the netcdf
# nc_close(ncPhen)
# 
# # Reshape the 3D array so we can map it, change the time field to be date
# dimnames(stiGFDL) <- list(lat = lat)#, time = tim)
# stiGFDL <- reshape2::melt(stiGFDL, value.name = "STI")
# dimnames(lusiGFDL) <- list(lat = lat)#, time = tim)
# lusiGFDL <- reshape2::melt(lusiGFDL, value.name = "LUSI")
# # beutiGFDL$dt <- as.Date("1970-01-01") + days(beutiGFDL$time) !RW: not in date format
# phenolGFDL <- stiGFDL %>% full_join(y = lusiGFDL, by = c("lat", "Var2")) %>%
#                 mutate(#year = year(dt),
#                        #month = month(dt),
#                        ESM = "GFDL")


# gfdlPhen <- tidync(paste0(datPath, "cuti_phenology_gfdl_1980-2100.nc"))
# varTbl <- gfdlPhen %>% hyper_tibble()
# timeTbl <- gfdlPhen %>% activate("D0") %>% hyper_tibble()
# varTbl <- varTbl %>% full_join(y = timeTbl, by = "time") %>%
#   mutate(ESM = esm)

gfdlPhen <- GetUpwelling.nc(path = paste0(datPath, "cuti_phenology_gfdl_1980-2100.nc"),
                            esm = "GFDL")
ipslPhen <- GetUpwelling.nc(path = paste0(datPath, "cuti_phenology_ipsl_1980-2100.nc"),
                            esm = "IPSL")
hadPhen <- GetUpwelling.nc(path = paste0(datPath, "cuti_phenology_had_1980-2100.nc"),
                            esm = "HAD")

# ncPhen <- nc_open(paste0(datPath, "cuti_phenology_ipsl_1980-2100.nc"))
# # Extract the spatiotemporal variables
# lat <- ncvar_get(ncPhen, "latitude")
# yr <- ncvar_get(ncPhen, "year")
# # Get the SDM predictions
# stiIPSL <- ncvar_get(ncPhen, "STI")
# lusiIPSL <- ncvar_get(ncPhen, "LUSI")
# # Close the netcdf
# nc_close(ncPhen)
# 
# # Reshape the 3D array so we can map it, change the time field to be date
# dimnames(stiIPSL) <- list(lat = lat)#, time = tim)
# stiIPSL <- reshape2::melt(stiIPSL, value.name = "STI")
# dimnames(lusiIPSL) <- list(lat = lat)#, time = tim)
# lusiIPSL <- reshape2::melt(lusiIPSL, value.name = "LUSI")
# dimnames(yr) <- list(lat = lat)#, time = tim)
# lusiIPSL <- reshape2::melt(lusiIPSL, value.name = "LUSI")
# # beutiGFDL$dt <- as.Date("1970-01-01") + days(beutiGFDL$time) !RW: not in date format
# phenolIPSL <- stiIPSL %>% full_join(y = lusiIPSL, by = c("lat", "Var2")) %>%
#                 mutate(#year = year(dt),
#                        #month = month(dt),
#                        ESM = "IPSL")
# 
# ncPhen <- nc_open(paste0(datPath, "cuti_phenology_had_1980-2100.nc"))
# # Extract the spatiotemporal variables
# lat <- ncvar_get(ncPhen, "latitude")
# tim <- ncvar_get(ncPhen, "time")
# # Get the SDM predictions
# stiHAD <- ncvar_get(ncPhen, "STI")
# lusiHAD <- ncvar_get(ncPhen, "LUSI")
# # Close the netcdf
# nc_close(ncPhen)
# 
# # Reshape the 3D array so we can map it, change the time field to be date
# dimnames(stiHAD) <- list(lat = lat)#, time = tim)
# stiHAD <- reshape2::melt(stiHAD, value.name = "STI")
# dimnames(lusiHAD) <- list(lat = lat)#, time = tim)
# lusiHAD <- reshape2::melt(lusiHAD, value.name = "LUSI")
# # beutiGFDL$dt <- as.Date("1970-01-01") + days(beutiGFDL$time) !RW: not in date format
# phenolHAD <- stiHAD %>% full_join(y = lusiHAD, by = c("lat", "Var2")) %>%
#                 mutate(#year = year(dt),
#                        #month = month(dt),
#                        ESM = "HAD")

upwelling <- upwelling %>% # need to aggregate to season first
                filter(latitude %in% c(33, 39), month %in% 5:7) %>%
                group_by(latitude, year, ESM) %>%
                summarize(BEUTI = mean(BEUTI),
                          CUTI = mean(CUTI))
 


  
upwelling <- upwelling %>% full_join(y = rbind(gfdlPhen, ipslPhen, hadPhen),
                                     by = c("latitude", "ESM", "year")) %>%
                filter(latitude %in% c(33, 36, 39)) %>% 
                pivot_wider(id_cols = c(year, ESM), names_from = latitude, 
                            values_from = c(BEUTI, CUTI, STI, LUSI))

# HCI ---------------------------------------------------------------------
hci <- read_csv(paste0(datPath, "HCI_rolling_fixed.csv"))

hci <- hci %>% filter(month %in% 2:7) %>% 
          group_by(region, year, projection) %>%
          summarize(HCI = mean(hci.rolling)) %>%
          rename(ESM = projection) %>%
          mutate(ESM = toupper(ESM),
                 region = sub(pattern = "egion ", replacement = "", x = region)) %>%
          pivot_wider(id_cols = c(year, ESM), names_from = region, 
                      names_prefix = "HCI_", values_from = HCI)

# Join together and calculate zscores 

# some rearranging for joining data
nurseTS <- nurseTS %>% pivot_longer(cols = -year, names_to = c("spp", "ESM"), names_sep = "_") %>%
              pivot_wider(id_cols = c(year, ESM), names_from = spp, values_from = value)
spawnTS <- spawnTS %>% pivot_longer(cols = -year, names_to = c("spp", "ESM"), names_sep = "_") %>%
              pivot_wider(id_cols = c(year, ESM), names_from = spp, values_from = value)

projDat <- sstProj %>% full_join(y= transport, by = c("year", "ESM")) %>%
              full_join(y= nemuroZ, by = c("year", "ESM")) %>%
              full_join(y = upwelling, by = c("year", "ESM")) %>%
              full_join(y = hci, by = c("year", "ESM")) %>%
              full_join(y = nurseTS, by = c("year", "ESM")) %>%
              full_join(y = spawnTS, by = c("year", "ESM")) %>%
              full_join(sardSpawnDays, by = c("year", "ESM")) %>%
              full_join(anchSpawnDays, by = c("year", "ESM")) 

projDat <- projDat %>% mutate(avgNearTransspring = as.numeric(avgNearshTransp_spring),
                              avgNearTranssummer = as.numeric(avgNearshTransp_summer),
                              avgOffTransspring = as.numeric(avgOffshTransp_spring),
                              avgOffTranssummer = as.numeric(avgOffshTransp_summer))

# read training dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

# select only projectable variables 
allDat <- datDFA %>% filter(year %in% 1980:2019) %>%
            select(year, springSST, summerSST, avgNearTransspring, avgNearTranssummer,
                   avgOffTransspring, avgOffTranssummer, ZM_NorCal, ZM_SoCal,
                   HCI_R3, HCI_R4, BEUTI_33N, BEUTI_39N, CUTI_33N, CUTI_39N,
                   OC_LUSI_33N, OC_LUSI_36N, OC_LUSI_39N, OC_STI_33N, OC_STI_36N, 
                   OC_STI_39N, sardSpawnHab, anchSpawnHab, daysAbove5pct, 
                   daysAbove40pct, sardNurseHab, anchNurseHab)

# check mean and sd in historical period are same
allDat %>% pivot_longer(-year, names_to = "var", values_to = "vals") %>%
  group_by(var) %>%
  summarize(histMean = mean(vals, na.rm = TRUE),
            histSD = sd(vals, na.rm = TRUE))

projHistSmry <- projDat %>% select(year, ESM, springSST, summerSST, avgNearTransspring, avgNearTranssummer,
                                   avgOffTransspring, avgOffTranssummer, ZM_NorCal, ZM_SoCal,
                                   HCI_R3, HCI_R4, BEUTI_33, BEUTI_39, CUTI_33, CUTI_39,
                                   LUSI_33, LUSI_36, LUSI_39, STI_33, STI_36, 
                                   STI_39, sardSpawnHab, anchSpawnHab, daysAbove5pct, 
                                   daysAbove40pct, sardNurseHab, anchNurseHab) %>%
                  pivot_longer(-c(year, ESM), names_to = "var", values_to = "vals") %>%
                  filter(year < 2020) %>%
                  group_by(var, ESM) %>%
                  summarize(histMean = mean(vals, na.rm = TRUE),
                            histSD = sd(vals, na.rm = TRUE))
## RW!: mean and SD are pretty close but not spot on. Which should be used for zscoring?
##      Do we need to bias correct?

# For now, zscore using mean of historical period from projection set
projDat <- projDat %>% select(year, ESM, springSST, summerSST, avgNearTransspring, avgNearTranssummer,
                              avgOffTransspring, avgOffTranssummer, ZM_NorCal, ZM_SoCal,
                              HCI_R3, HCI_R4, BEUTI_33, BEUTI_39, CUTI_33, CUTI_39,
                              LUSI_33, LUSI_36, LUSI_39, STI_33, STI_36, 
                              STI_39, sardSpawnHab, anchSpawnHab, daysAbove5pct, 
                              daysAbove40pct, sardNurseHab, anchNurseHab) %>%
  pivot_longer(-c(year, ESM), names_to = "var", values_to = "vals") %>%
  full_join(y = projHistSmry, by = c("var", "ESM")) %>%
  mutate(scaled = (vals - histMean)/histSD) %>%
  select(year, ESM, var, scaled) %>%
  pivot_wider(names_from = var, values_from = scaled)

# add other empty variables
projDat[c("NCOPspring", "NCOPsummer", "SCOPspring", "SCOPsummer",  
          "swfscRockfishSurv_Myctophids", "avgSSWIspring", "avgSSWIsummer", 
          "sardLarv", "anchLarv", "mesopelLarv", "anchYoY", "age1SprSardmeanWAA", 
          "meanSSBwt", "copBio", "naupBio", 
          "anchRec", "sardRec", "albacore", "hake")] <- NA

# select variables in same order as provided to model
projDat <- projDat %>% filter(year >= 2020) %>% 
  select(year, ESM, "HCI_R3", "HCI_R4", "NCOPspring", "NCOPsummer", "SCOPspring",
         "SCOPsummer", "BEUTI_33", "BEUTI_39", "CUTI_33", "CUTI_39",
         "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
         "STI_36", "STI_39", "swfscRockfishSurv_Myctophids",
         "avgSSWIspring", "avgSSWIsummer", "sardLarv", "anchLarv",
         "mesopelLarv", "anchYoY", "age1SprSardmeanWAA", "meanSSBwt",
         "copBio", "naupBio", "ZM_NorCal", "ZM_SoCal", "sardSpawnHab",
         "anchSpawnHab", "daysAbove5pct", "daysAbove40pct", 
         "sardNurseHab", "anchNurseHab", "anchRec", "sardRec", 
         "springSST", "summerSST", "avgNearTransspring", 
         "avgNearTranssummer", "avgOffTransspring", "avgOffTranssummer",
         "albacore", "hake")

# write_csv(projDat, file = "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/formattedDFAprojDat.csv")

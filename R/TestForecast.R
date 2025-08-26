# test forecast performance based on different projection horizons
# Created: 8/22/2025, Robert Wildermuth

library(tidyverse)
library(MARSS)
source("R/OSAResids.R")
source("R/LFOXV.R")

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

# Assess LFOIC for last 10 years of data
peels <- 10

# keep observation data indicators, remove associated model-derived indicators
localModel <- c("HCI_R4", "NCOPspring", "NCOPsummerlag1", "SCOPsummerlag1", 
                "BEUTI_33N", "BEUTI_39N", "CUTI_33N", "OC_LUSI_33N", "OC_LUSI_39N",
                "OC_STI_33N", "OC_STI_39N", "swfscRockfishSurv_Myctophids",
                "avgSSWIspring", "avgSSWIsummer", "sardLarv", "anchLarv", 
                "mesopelLarv", "anchYoY", "age1SprSardmeanWAA", "meanSSBwt",
                "copBio", "naupBio", "sardSpawnHab", "anchSpawnHab", 
                "daysAbove5pct", "daysAbove40pct", "sardNurseHab", "anchNurseHab",
                "anchRec", "sardRec", "summerSST", "albacore", "hake")

# remove observation indicators, keep associated model-derived indicators
projectModel <- c("HCI_R4", "BEUTI_33N", "BEUTI_39N", "CUTI_33N", "OC_LUSI_33N",
                  "OC_LUSI_39N", "OC_STI_33N", "OC_STI_39N", "ZM_NorCal",
                  "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", "daysAbove5pct",
                  "daysAbove40pct", "sardNurseHab", "anchNurseHab", "anchRec", 
                  "sardRec", "summerSST", "avgNearTransspring", "avgNearTranssummer", 
                  "avgOffTransspring", "avgOffTranssummer")

# Data for full historical dataset ---------------------------------------

locDat <- datDFA%>% select(year, all_of(localModel), all_of(projectModel))

datNames <- names(locDat)[-1]

# subset  1990 to 2019
locDat <- locDat %>% filter(year %in% 1990:2019)

# transpose for MARSS formatting
itDat <- locDat %>% select(-year) %>% t()

locEqRMSE1 <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = 2,
                  peels = peels,
                  horizon = 1) 
locEqRMSE <- locEqRMSE1 %>% filter(variable == "totRMSE", 
                                    resType == "resid.Inf",
                                    predHoriz == 1)
locEqRMSE2 <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = 2,
                  peels = peels,
                  horizon = 2)
locEqRMSE <- locEqRMSE2 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 2) %>%
                bind_rows(locEqRMSE)
locEqRMSE3 <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = 2,
                  peels = peels,
                  horizon = 3)
locEqRMSE <- locEqRMSE3 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 3) %>%
                bind_rows(locEqRMSE)
locEqRMSE4 <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = 2,
                  peels = peels,
                  horizon = 4)
locEqRMSE <- locEqRMSE4 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 4) %>%
                bind_rows(locEqRMSE)
locEqRMSE5 <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = 2,
                  peels = peels,
                  horizon = 5)
locEqRMSE <- locEqRMSE5 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 5) %>%
                bind_rows(locEqRMSE)

locEqRMSE$model <- "Local"

# Data for historical projection dataset ---------------------------------------

projDat <- datDFA%>% select(year, all_of(projectModel))

datNames <- names(projDat)[-1]

# subset  1990 to 2019
projDat <- projDat %>% filter(year %in% 1990:2019)

# transpose for MARSS formatting
itDat <- projDat %>% select(-year) %>% t()

projEqRMSE1 <- LFOXV(dfaDat = itDat,
                  Rstructure = "diagonal and equal",
                  mTrends = 3,
                  peels = peels,
                  horizon = 1)
projEqRMSE <- projEqRMSE1 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 1)
projEqRMSE2 <- LFOXV(dfaDat = itDat,
                    Rstructure = "diagonal and equal",
                    mTrends = 3,
                    peels = peels,
                    horizon = 2)
projEqRMSE <- projEqRMSE2 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 2) %>%
                bind_rows(projEqRMSE)
projEqRMSE3 <- LFOXV(dfaDat = itDat,
                    Rstructure = "diagonal and equal",
                    mTrends = 3,
                    peels = peels,
                    horizon = 3)
projEqRMSE <- projEqRMSE3 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 3) %>%
                bind_rows(projEqRMSE)
projEqRMSE4 <- LFOXV(dfaDat = itDat,
                    Rstructure = "diagonal and equal",
                    mTrends = 3,
                    peels = peels,
                    horizon = 4)
projEqRMSE <- projEqRMSE4 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 4) %>%
                bind_rows(projEqRMSE)
projEqRMSE5 <- LFOXV(dfaDat = itDat,
                    Rstructure = "diagonal and equal",
                    mTrends = 3,
                    peels = peels,
                    horizon = 5)
projEqRMSE <- projEqRMSE5 %>% filter(variable == "totRMSE", 
                                   resType == "resid.Inf",
                                   predHoriz == 5) %>%
                bind_rows(projEqRMSE)
projEqRMSE$model <- "Projection"
# Calculate Persistence Prediction RMSE -----------------------------------

datLen <- length(1990:2019)

perstResids <- datDFA %>% select(year, sardRec, anchRec) %>% 
                  filter(year >= 1990, year < 2020) %>% 
                  mutate(zscoreSardRec = zscore(sardRec),
                         zscoreAnchRec = zscore(anchRec),
                         perst1Sard = c(NA, zscoreSardRec[1:(datLen-1)]),
                         perst1Anch = c(NA, zscoreAnchRec[1:(datLen-1)]),
                         resid1Sard = zscoreSardRec - perst1Sard,
                         resid1Anch = zscoreAnchRec - perst1Anch,
                         perst2Sard = c(NA,NA, zscoreSardRec[1:(datLen-2)]),
                         perst2Anch = c(NA,NA, zscoreAnchRec[1:(datLen-2)]),
                         resid2Sard = zscoreSardRec - perst2Sard,
                         resid2Anch = zscoreAnchRec - perst2Anch,
                         perst3Sard = c(NA,NA,NA, zscoreSardRec[1:(datLen-3)]),
                         perst3Anch = c(NA,NA,NA, zscoreAnchRec[1:(datLen-3)]),
                         resid3Sard = zscoreSardRec - perst3Sard,
                         resid3Anch = zscoreAnchRec - perst3Anch,
                         perst4Sard = c(NA,NA,NA,NA, zscoreSardRec[1:(datLen-4)]),
                         perst4Anch = c(NA,NA,NA,NA, zscoreAnchRec[1:(datLen-4)]),
                         resid4Sard = zscoreSardRec - perst4Sard,
                         resid4Anch = zscoreAnchRec - perst4Anch,
                         perst5Sard = c(NA,NA,NA,NA,NA, zscoreSardRec[1:(datLen-5)]),
                         perst5Anch = c(NA,NA,NA,NA,NA, zscoreAnchRec[1:(datLen-5)]),
                         resid5Sard = zscoreSardRec - perst5Sard,
                         resid5Anch = zscoreAnchRec - perst5Anch,
                         resid0SardMean = zscoreSardRec - 0, # residual from mean S-R value is just the rec dev
                         resid0AnchMean = zscoreAnchRec - 0) %>% # residual from mean S-R value is just the rec dev
                  select(year, resid1Sard, resid2Sard, resid3Sard, resid4Sard, 
                         resid5Sard, resid0SardMean, resid1Anch, resid2Anch, 
                         resid3Anch, resid4Anch, resid5Anch, resid0AnchMean) %>% 
                  pivot_longer(cols = -year, names_prefix = "resid", names_sep = 1, 
                               names_to = c("predHoriz", "variable"), 
                               values_to = "perstResid") %>% 
                  filter(year %in% 2010:2019) %>%
                  group_by(variable, predHoriz) %>%
                  summarize(sosRes = sum(perstResid^2, na.rm = TRUE),
                            peels = sum(!is.na(perstResid))) %>%
                  mutate(RMSE = sqrt(sosRes/peels),
                         resType = "resid.Perst",
                         predHoriz = as.numeric(predHoriz)) %>%
                  select(-sosRes)

perstResids <- perstResids %>% group_by(resType, predHoriz) %>%
                  summarize(totRMSE = sum(RMSE))

testFcast <- bind_rows(locEqRMSE, projEqRMSE) #, perstResids)
# plot out best performing model structures over prediction horizons
testFcast %>% filter(resType %in% "resid.Inf", variable == "totRMSE") %>%
  ggplot(aes(x = predHoriz, y = RMSE)) +
  geom_line(aes(color = as.character(model)), linewidth = 1) + #paste(mTrends, Rstructure, sep = "-"))) +
  geom_point(data = perstResids %>% filter(predHoriz != 0), aes(y = totRMSE)) +
  geom_hline(yintercept = perstResids %>% filter(predHoriz == 0) %>% pull(totRMSE)) +
  # scale_color_viridis_d() +
  labs(color = "Model", x = "Prediction Horizon") +
  theme_classic()

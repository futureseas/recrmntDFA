# Model selection of recruitment DFAs based on LFOIC
# Created: 12/14/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)
source("R/LFOXV.R")
source("R/OSAResids.R")

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

# allDat <- datDFA%>% select(year, all_of(localModel), all_of(projectModel))
allDat <- datDFA%>% select(year, all_of(localModel), all_of(projectModel),
                           ENSO, NPGO, PDOspring, PDOsummer)

datNames <- names(allDat)[-1]

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- matrix(list(0),length(datNames),length(datNames))
diag(Rcustom) <- c("HCI", 
                   "COP", "COP", "COP",
                   "BEUTI", "BEUTI",
                   "CUTI", 
                   "LUSI", "LUSI",
                   "STI", "STI", 
                   "RREAS",
                   "SSWI", "SSWI",
                   "CalCOFI", "CalCOFI", "CalCOFI",
                   "RREAS",
                   "WAA", "WAA",
                   "PRPOOS", "PRPOOS", 
                   "sardSDM", "anchSDM", "sardSDM", "anchSDM",
                   "sardlarvSDM", "anchlarvSDM",
                   "anchRec",
                   "sardRec",
                   "SST",
                   "Alb", "Hake",
                   "NEMURO", "NEMURO",
                   "Transp", "Transp", "Transp", "Transp", #)
                   "Basin", "Basin", "Basin", "Basin")

# Data for historical projection dataset ---------------------------------------

# # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
# allDat <- datDFA%>% select(year, all_of(projectModel))
# 
# datNames <- names(allDat)[-1]
# 
# # Create a custom R obs error matrix assuming each data source has it's own common error
# Rcustom <- matrix(list(0),length(datNames),length(datNames))
# diag(Rcustom) <- c("HCI",
#                    "BEUTI", "BEUTI",
#                    "CUTI",
#                    "LUSI", "LUSI",
#                    "STI", "STI",
#                    "NEMURO", "NEMURO",
#                    "sardSDM", "anchSDM", "sardSDM", "anchSDM",
#                    "sardlarvSDM", "anchlarvSDM",
#                    "sardRec", "anchRec",
#                    "SST",
#                    "Transp", "Transp", "Transp", "Transp")


# LFOIC  ------------------------------------------------------------------

# Set up table of model structures to hold LFOIC vals
xvModSel <- tibble(initYr = 0, 
                   Rstructure = "", 
                   mTrends = 0)[0,]

# loop over initial dates
for(y in c(#1980, 1985, 
           1990)){
  cat("\n")
  print(y)
  
  # subset from 1980, 1985, or 1990 to 2019
  initDat <- allDat %>% filter(year %in% y:2019)
  
  # transpose for MARSS formatting
  itDat <- initDat %>% select(-year) %>% t()
  
  # loop over number of trends
  for(m in 1:7){
    cat("\n Trends: ", m)
    cat("\n Diagonal and equal R matrix")
    itEqRMSE <- LFOXV(dfaDat = itDat,
                       Rstructure = "diagonal and equal",
                       mTrends = m,
                       peels = peels)

    itEqRMSE <- itEqRMSE %>% mutate(initYr = y,
                                    mTrends = m,
                                    Rstructure = "diag & equal")

    itEqRMSE %>% filter(resType == "resid.Inf", variable == "totRMSE") %>% print()

    cat("\n Diagonal and unequal R matrix")

    itUneqRMSE <- LFOXV(dfaDat = itDat,
                      Rstructure = "diagonal and unequal",
                      mTrends = m,
                      peels = peels)

    itUneqRMSE <- itUneqRMSE %>% mutate(initYr = y,
                                    mTrends = m,
                                    Rstructure = "diag & unequal")

    itUneqRMSE %>% filter(resType == "resid.Inf", variable == "totRMSE") %>% print()

    cat("\n Custom R matrix")

    itCustRMSE <- LFOXV(dfaDat = itDat,
                        Rstructure = Rcustom,
                        mTrends = m,
                        peels = peels)

    itCustRMSE <- itCustRMSE %>% mutate(initYr = y,
                                    mTrends = m,
                                    Rstructure = "custom R by SDM")
    itCustRMSE %>% filter(resType == "resid.Inf", variable == "totRMSE") %>% print()

    xvModSel <- xvModSel %>% bind_rows(itUneqRMSE, itCustRMSE,
                                       itEqRMSE)
    # xvModSel <- xvModSel %>% bind_rows(itUneqRMSE, itCustRMSE)
  } # end trends loop
} # end year loop 

xvModSel <- xvModSel %>% mutate(nIndices = length(datNames),
                                peels = peels)
  
# write_csv(xvModSel, file = "fullHistoricalModelSelection_LocalandProj.csv")
write_csv(xvModSel, file = "oceanBasinHistoricalModelSelection_LocalandProj.csv")
# write_csv(xvModSel, file = "histProjectionModelSelection.csv")


# Calculate Persistence Prediction RMSE -----------------------------------

allDat %>% select(year, sardRec, anchRec) %>% 
  filter(year >= 1990, year < 2020) %>% 
  mutate(zscoreSardRec = zscore(sardRec),
         zscoreAnchRec = zscore(anchRec),
         perstSard = c(NA, zscoreSardRec[1:29]),
         perstAnch = c(NA, zscoreAnchRec[1:29]),
         residSard = zscoreSardRec - perstSard,
         residAnch = zscoreAnchRec - perstAnch) %>% 
  select(year, residAnch, residSard) %>% 
  pivot_longer(cols = -year, names_prefix = "resid", names_to = "VoI", 
               values_to = "perstResid") %>% 
  filter(year %in% 2010:2019) %>%
  group_by(VoI) %>%
  summarize(sosRes = sum(perstResid^2, na.rm = TRUE),
            nObs = sum(!is.na(perstResid))) %>%
  mutate(RMSE = sqrt(sosRes/nObs)) %>% pull(RMSE) %>% sum() #%>% summarize(totRMSE = sum(RMSE))


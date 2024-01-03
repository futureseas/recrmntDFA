# Model selection of recruitment DFAs based on LFOIC
# Created: 12/14/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)
source("R/LFOXV.R")

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

# Assess LFOIC for last 10 years of data
peels <- 10


# LFOIC for full historical dataset ---------------------------------------

# remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
allDat <- datDFA%>% select(-c(NOI,
                              ENSO,
                              NPGO,
                              #All_Copepods, # ~same as calanoid copepods
                              #euphausiids, # too large for larval mouth gape
                              anchBioSmrySeas1,
                              sardBioSmrySeas1,
                              anchBioSmrySeas2, # leave out biomass since not fit well
                              sardBioSmrySeas2,
                              # NCOPspring,
                              # SCOPspring, # summer copepod index had higher loadings
                              PDOsummer, # lower loading than spring, may want to try a lag
                              PDOspring,
                              NCOPsummer,
                              SCOPsummer,
                              # BEUTI_33N, # oceanography at 39N had highest loadings
                              # OC_LUSI_33N,
                              # OC_LUSI_36N,
                              # OC_STI_33N,
                              # OC_STI_36N,
                              copMeanSize,
                              copPropBioSize,
                              naupMeanSize,
                              naupPropBioSize,
                              ZL_NorCal,
                              ZL_SoCal,
                              age1SprAnchmeanWAA)) # not enough data in time windowselect(-c(#sprCalCOFILarvalSardine,

datNames <- names(allDat)[-1]

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- matrix(list(0),length(datNames),length(datNames)) 
diag(Rcustom) <- c("HCI", "HCI",
                   "COP", "COP", "COP", "COP",
                   "BEUTI", "BEUTI", 
                   "CUTI", "CUTI",
                   "LUSI", "LUSI", "LUSI", 
                   "STI", "STI", "STI",
                   "RREAS",
                   "SSWI", "SSWI",
                   "CalCOFI", "CalCOFI", "CalCOFI",
                   "RREAS",
                   "WAA", "WAA",
                   "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS",
                   "NEMURO", "NEMURO",
                   "SDM", "SDM", "TIME", "TIME", "SDM", "SDM",
                   "anchRec", 
                   "sardRec",
                   # "anchBio",
                   # "sardBio",
                   "SST", "SST",
                   "Transp", "Transp", "Transp", "Transp",
                   "Alb", "Hake")

# Set up table of model structures to hold LFOIC vals
xvModSel <- expand_grid(initYr = c(1980, 1985, 1990), 
                        Rstructure = c("Rcustom", "diag & equal", "diag & unequal"), 
                        mTrends = 1:5) %>%
              mutate(nIndices = length(datNames),
                     peels = peels,
                     LFOIC = NA)

# loop over initial dates
for(y in c(1980, 1985, 1990)){
  cat("\n")
  print(y)
  
  # subset from 1980, 1985, or 1990 to 2019
  initDat <- allDat %>% filter(year %in% y:2019)
  
  # transpose for MARSS formatting
  itDat <- initDat %>% select(-year) %>% t()
  
  # loop over number of trends
  for(m in 1:5){
    cat("\n Trends: ", m)
    cat("\n Diagonal and equal R matrix")
    itEqRMSE <- LFOXV(dfaDat = itDat, 
                       Rstructure = "diagonal and equal", 
                       mTrends = m, 
                       peels = peels)
    
    xvModSel[xvModSel$initYr == y &
               xvModSel$mTrends == m &
               xvModSel$Rstructure == "diag & equal", "LFOIC"] <- itEqRMSE
    
    cat("\n Diagonal and unequal R matrix")
    
    itUneqRMSE <- LFOXV(dfaDat = itDat, 
                      Rstructure = "diagonal and unequal", 
                      mTrends = m, 
                      peels = peels)
    
    xvModSel[xvModSel$initYr == y &
               xvModSel$mTrends == m &
               xvModSel$Rstructure == "diag & unequal", "LFOIC"] <- itUneqRMSE
    
    cat("\n Custom R matrix")
    
    itCustRMSE <- LFOXV(dfaDat = itDat, 
                      Rstructure = Rcustom, 
                      mTrends = m, 
                      peels = peels)
    
    xvModSel[xvModSel$initYr == y &
               xvModSel$mTrends == m &
               xvModSel$Rstructure == "Rcustom", "LFOIC"] <- itCustRMSE
  }
}
 
  
write_csv(xvModSel, file = "fullHistoricalModelSelection.csv")



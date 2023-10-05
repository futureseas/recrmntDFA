# Fit Bayesian DFA for recruitment indicators
# Created: 9/29/2023, Robert Wildermuth

library(tidyverse)
library(bayesdfa)

datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

allDat <- datDFA %>% filter(year %in% 1980:2019) %>%
            # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
            select(-c(NOI,
                      ENSO,
                      NPGO,
                      #All_Copepods, # ~same as calanoid copepods
                      #euphausiids, # too large for larval mouth gape
                      anchBioSmrySeas1,
                      sardBioSmrySeas1,
                      # anchBioSmrySeas2, # leave out biomass since not fit well
                      # sardBioSmrySeas2,
                      # NCOPspring,
                      # SCOPspring, # summer copepod index had higher loadings
                      PDOsummer, # lower loading than spring, may want to try a lag
                      PDOspring,
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
                      # age1SprAnchmeanWAA # not enough data in time window
                      )) 

datNames <- names(allDat)[-1]

# summarize data density by variable and year
datDense <- allDat %>% pivot_longer(cols = -year) %>%
              mutate(datAvail = !is.na(value)) 
with(datDense ,
  table(name, datAvail))
with(datDense ,
     table(year, datAvail)) 
datDense %>% group_by(year) %>%
  summarize(nInd = n(),
            nMissing = nInd - sum(datAvail)) %>%
  ggplot(aes(x = year, y = nMissing)) +
  geom_col()

# transpose for MARSS formatting
allDat <- allDat %>% select(-year) %>% t()

chains <- 3
iter <- 10

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- c(1, #"HCI",
             rep(2, 4), #"COP", "COP", "COP", "COP",
             3, 3, #"BEUTI", "BEUTI", 
             4, 4, #"CUTI", "CUTI",
             5, 5, 5, #"LUSI", "LUSI", "LUSI", 
             6, 6, 6, #"STI", "STI", "STI",
             7, #"RREAS",
             8, 8, #"SSWI", "SSWI",
             9, 9, 9, #"CalCOFI", "CalCOFI", "CalCOFI",
             7, #"RREAS",
             10, 10, 10, #"WAA", "WAA",
             11, 11, #"PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS",
             12, 12, #"NEMURO", "NEMURO",
             13, 13, #"SDM", "SDM", 
             14, 14, #"TIME", "TIME", 
             13, 13, #"SDM", "SDM",
             15, #"anchRec", 
             16, #"sardRec",
             17, #"anchBio",
             18, #"sardBio",
             19, 19, #"SST", "SST",
             rep(20, 4), #"Transp", "Transp", "Transp", "Transp",
             21, #"Alb",
             22) #"Hake"

fit1trend <- fit_dfa(
                y = allDat, num_trends = 1, scale="zscore",
                varIndx = rep(1, nrow(allDat)), #rep_len(1:2, nrow(allDat)),
                # varIndx = 1:47, # mapping of unique R matrix variances to indicators
                # varIndx = Rcustom,
                iter = iter, chains = chains, thin = 1)

# check convergence 
is_converged(fit1trend)

# Rotate trends
trendsRot <- rotate_trends(fit1trend)

plot_trends(trendsRot)

# plot the fitted values
plot_fitted(fit1trend)

# plot the loadings
plot_loadings(trendsRot) + 
  theme_classic() + coord_flip()

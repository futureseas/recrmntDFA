# Fit Bayesian DFA for recruitment indicators
# Created: 9/29/2023, Robert Wildermuth

library(tidyverse)
library(bayesdfa)

datDFA <- read_csv("recrDFAdat.csv")

allDat <- datDFA %>% filter(year %in% 1980:2019) %>%
            # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
            select(-c(NOI,
                      ENSO,
                      NPGO,
                      anchBioSmrySeas1,
                      sardBioSmrySeas1,
                      anchBioSmrySeas2, # accounted for in S-R relationship
                      sardBioSmrySeas2,
                      # NCOPspring,
                      # SCOPspring, # summer copepod index had higher loadings
                      NCOPsummer, # only have NEMURO pojections of spring
                      SCOPsummer,
                      PDOsummer, 
                      PDOspring,
                      BEUTI_39N, # keep more local variables
                      CUTI_39N,
                      OC_LUSI_39N,
                      OC_STI_39N,
                      # avgNearTransspring, # spring transport had lower MARSS loadings       
                      # avgOffTransspring,
                      # naupBio, # had low MARSS loadings
                      copMeanSize,
                      copPropBioSize,
                      naupMeanSize,
                      naupPropBioSize,
                      ZL_NorCal,
                      ZL_SoCal
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

# overall proportion of missing data
sum(is.na(allDat))/(dim(allDat)[1]*dim(allDat)[2])

# transpose for MARSS formatting
allDat <- allDat %>% select(-year) %>% t()

chains <- 4
iter <- 5000

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- c(1, #"HCI",
             2, 2, #rep(2, 4), #"COP", "COP", "COP", "COP",
             3, #3, #"BEUTI", "BEUTI", 
             4, #4, #"CUTI", "CUTI",
             5, 5, #5, #"LUSI", "LUSI", "LUSI", 
             6, 6, #6, #"STI", "STI", "STI",
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
             #17, #"anchBio",
             #18, #"sardBio",
             17, 17, #19, 19, #"SST", "SST",
             rep(18, 4), #"Transp", "Transp", "Transp", "Transp",
             19, #"Alb",
             20) #"Hake"

fit4trends <- fit_dfa(
                y = allDat, num_trends = 1, scale="zscore",
                varIndx = rep(1, nrow(allDat)), #rep_len(1:2, nrow(allDat)),
                # varIndx = 1:nrow(allDat), # mapping of unique R matrix variances to indicators
                # varIndx = Rcustom,
                estimate_trend_ar = TRUE,
                iter = iter, chains = chains, thin = 1,
                cores = parallel::detectCores()-2)
# RW: will not fit models with more than ~40 index time series


# check convergence 
is_converged(fit4trends)

# need convergence diagnostics
# thinning?
shinystan::launch_shinystan(fit4trends$model)

# Rotate trends
trendsRot <- rotate_trends(fit4trends)

plot_trends(trendsRot)

# plot the fitted values
plot_fitted(fit4trends)

# plot the loadings
plot_loadings(trendsRot) + 
  theme_classic() + coord_flip()

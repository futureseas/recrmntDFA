# Fit Bayesian DFA for recruitment indicators
# Created: 9/29/2023, Robert Wildermuth

library(tidyverse)
library(bayesdfa)
library(rstan)

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

# long format data_shape
longDat <- allDat %>% #select(year, HCI, OC_LUSI_33N, ZM_SoCal) %>%
              pivot_longer(cols = -year, names_to = "VoI", values_to = "obs") %>%
              rename(time = year) %>%
              mutate(ts = as.numeric(factor(VoI, levels = datNames)),
                     time = time- 1979)

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

fit1trend <- fit_dfa(
                y = longDat[, c("obs", "ts", "time")], # allDat,#
                num_trends = 1, scale="zscore",
                # varIndx = rep(1, nrow(allDat)), #rep_len(1:2, nrow(allDat)),
                # varIndx = 1:nrow(allDat), # mapping of unique R matrix variances to indicators
                varIndx = Rcustom,
                est_correlation = FALSE,
                estimate_trend_ar = TRUE,
                data_shape = "long",
                iter = iter, chains = chains, thin = 3,
                cores = parallel::detectCores()-2)

fit2trends <- fit_dfa(
                y = longDat, #allDat, 
                num_trends = 2, scale="zscore",
                # varIndx = rep(1, nrow(allDat)), #rep_len(1:2, nrow(allDat)),
                # varIndx = 1:nrow(allDat), # mapping of unique R matrix variances to indicators
                varIndx = Rcustom,
                est_correlation = FALSE,
                estimate_trend_ar = TRUE,
                data_shape = "long",
                iter = iter, chains = chains, thin = 3,
                cores = parallel::detectCores()-2)
# RW: will not fit models with more than ~40 index time series
# save(fit2trends, file = "bayesFit_1980to2019_noBio_strend_Rcustom.RData")


fit3trends <- fit_dfa(
  y = allDat, num_trends = 3, scale="zscore",
  # varIndx = rep(1, nrow(allDat)), #rep_len(1:2, nrow(allDat)),
  # varIndx = 1:nrow(allDat), # mapping of unique R matrix variances to indicators
  varIndx = Rcustom,
  est_correlation = FALSE,
  estimate_trend_ar = TRUE,
  iter = iter, chains = chains, thin = 1,
  cores = parallel::detectCores()-2)

fit4trends <- fit_dfa(
  y = allDat, num_trends = 4, scale="zscore",
  # varIndx = rep(1, nrow(allDat)), #rep_len(1:2, nrow(allDat)),
  # varIndx = 1:nrow(allDat), # mapping of unique R matrix variances to indicators
  varIndx = Rcustom,
  est_correlation = TRUE,
  estimate_trend_ar = FALSE,
  iter = iter, chains = chains, thin = 1,
  cores = parallel::detectCores()-2)

# check convergence 
is_converged(fit2trends)

# need convergence diagnostics
# thinning?
shinystan::launch_shinystan(fit1trend$model)
shinystan::launch_shinystan(fit2trends$model)
shinystan::launch_shinystan(fit3trends$model)
# For now just look at monitor to evaluate convergence
npars <- dim(fit1trend$monitor)[1]
fit1trend$monitor[(npars-25):npars,1:22]
fit2trends$monitor[(npars-25):npars,1:22]
fit3trends$monitor[(npars-25):npars,1:22]
traceplot(fit1trend$model, pars = "xstar") + facet_wrap(~chain)
traceplot(fit2trends$model, pars = "phi") 
traceplot(fit3trends$model, pars = "xstar") 
stan_ac(fit1trend$model, pars = "Z")
stan_ac(fit2trends$model, pars = "Z")
stan_ac(fit3trends$model, pars = "sigma")

fit1trend$monitor %>% filter(Rhat > 1.05)

fit2trends$monitor %>% filter(Rhat > 1.05)
fit3trends$monitor %>% filter(Rhat > 1.05) %>% as.data.frame %>% select(mean, n_eff, Rhat)
fit3trends$monitor %>% filter(n_eff < 100) %>% as.data.frame %>% select(mean, n_eff, Rhat)

# see how chain flipping works
# flipped_chains = bayesdfa::invert_chains(fit1trends$model, trends = 1)
# head(fit4trends$model@sim$samples[[3]]$`x[1,1]`)
# head(flipped_chains$model@sim$samples[[3]]$`x[1,1]`)

# Rotate trends
trendsRot <- rotate_trends(fit2trends)

plot_trends(trendsRot) + geom_hline(yintercept = 0)

# plot the fitted values
plot_fitted(fit2trends,
            names = datNames)
dfa_fitted(fit2trends, names = datNames) %>% filter(ID %in% c("sardRec", "sardLarv",
                                                              "anchRec", "anchLarv",
                                                              "anchYoY"))

# plot the loadings
plot_loadings(trendsRot, names = datNames
              ,threshold = 0.95 # 95% of posterior density is +/- 0
              ) + 
  theme_classic() 
dfa_loadings(trendsRot, names = datNames) %>% filter(#trend == "Trend 2",
                                                     name %in% c("sardRec", "sardLarv",
                                                                 "anchRec", "anchLarv",
                                                                 "anchYoY"))#,
                                                     prob_diff0 > 0.5)

dfa_trends(trendsRot)

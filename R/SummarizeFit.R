# Summary information for MLE and Bayesian DFA loadings and fits to data
# Created: 11/8/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# load(file = "marssFit_1980to2019_noBio_3trend_Rcustom.RData")
# load(file = "marssFit_1980to2019_ProjDFA_3trend_Rcustom.RData")
load(file = "marssFit_1990to2019_noBio_5trend_DiagEql.RData")

smryDFA <- overallDFA

# Estimated fit of MLE DFA ------------------------------------------------

# Look at factor loadings
# get the inverse of the rotation matrix 
Z.est <- coef(smryDFA, type = "matrix")$Z 
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(smryDFA, type = "matrix")$Z)$rotmat
} 

# rotate factor loadings 
Z.rot <- Z.est %*% H.inv 
# rotate trends 
trends.rot <- solve(H.inv) %*% smryDFA$states

# Add CIs to marssMLE object 
smryDFA <- MARSSparamCIs(smryDFA) 
# Use coef() to get the upper and lower CIs 
Z.low <- coef(smryDFA, type = "Z", what = "par.lowCI") 
Z.up <- coef(smryDFA, type = "Z", what = "par.upCI") 
Z.rot.up <- Z.up %*% H.inv 
Z.rot.low <- Z.low %*% H.inv 
df <- data.frame(ind = rownames(Z.rot),
                 trend = rep(1:ncol(Z.rot), each = nrow(Z.rot)),
                 est = as.vector(Z.rot), 
                 conf.up = as.vector(Z.rot.up), 
                 conf.low = as.vector(Z.rot.low))

# Review significance of loadings and trends
df <- df %>% mutate(diff0 = case_when(est > 0 & conf.low > 0 ~ TRUE, # RW!: low CI not always < high CI
                                      est < 0 & conf.up < 0 ~ TRUE,
                                      .default = FALSE))

df %>% filter(ind %in% c("sardRec", "anchRec", "sardLarv", 
                         "anchLarv", "anchYoY")) %>%
  arrange(ind)

# significance of loadings for projectable indicators
df %>% 
  # filter(ind %in% c("HCI_R3", "HCI_R4", "BEUTI_33N", "BEUTI_39N", "CUTI_33N",
  #                        "CUTI_39N",
  #                        "OC_LUSI_33N", "OC_LUSI_36N", "OC_LUSI_39N", "OC_STI_33N",
  #                        "OC_STI_36N", "OC_STI_39N", "ZM_NorCal", "ZM_SoCal",
  #                        "sardSpawnHab",
  #                        "anchSpawnHab", "daysAbove5pct", "daysAbove40pct",
  #                        "sardNurseHab", "anchNurseHab",
  #                        "springSST", "summerSST", "avgNearTransspring",
  #                        "avgNearTranssummer", "avgOffTransspring", "avgOffTranssummer")) %>%
  arrange(ind) %>% 
  # group_by(ind) %>% summarize(anyDiff0 = sum(diff0)) %>% filter(anyDiff0 == 0)
  # filter(diff0 == 1, abs(est) < 0.05)
  # filter(diff0 == 1) %>% arrange(abs(est))
  # filter(diff0 == 1, abs(est) > 0.25) %>% pull(ind) %>% unique() #%>% length()
  group_by(ind) %>% summarize(sumLoadings = sum(abs(est))) %>% arrange(sumLoadings) %>% print(n=45)
  
# CUTI_33N has lowest overall trend loadings
df %>% filter(ind == "CUTI_33N")

# indicators with no significant loadings
df %>% group_by(ind) %>% summarize(anyDiff0 = sum(diff0)) %>% filter(anyDiff0 == 0)

# Find variables with highest, most precise loadings and lowest, least precise loadings
df <- df %>% mutate(interval = abs(conf.up-conf.low),
                    magnLoading = abs(est))
df %>% arrange(trend, interval, magnLoading)
hiLoading <- df %>% group_by(trend) %>% slice_max(magnLoading, n = 10) %>% print(n = 40)
mostPrecise <- df %>% group_by(trend) %>% slice_min(interval, n = 10) %>% print(n = 40)
loLoading <- df %>% group_by(trend) %>% slice_min(magnLoading, n = 10) %>% print(n = 40)
leastPrecise <- df %>% group_by(trend) %>% slice_max(interval, n = 10) %>% print(n = 40)

hiLoading <- hiLoading %>% group_by(ind) %>% summarize(counts = n()) %>% arrange(counts) %>% print(n = 30)
mostPrecise <- mostPrecise %>% group_by(ind) %>% summarize(counts = n()) %>% arrange(counts) %>% print(n = 30)
hiLoading %>% inner_join(y = mostPrecise, by = "ind") %>% print(n = 30)

loLoading <- loLoading %>% group_by(ind) %>% summarize(counts = n()) %>% arrange(counts) %>% print(n = 30)
leastPrecise <- leastPrecise %>% group_by(ind) %>% summarize(counts = n()) %>% arrange(counts) %>% print(n = 30)
loLoading %>% inner_join(y = leastPrecise, by = "ind") %>% print(n = 30)

# Calculate RMSE for each indicator on estimates conditioned on all the data (tT)
resids <- residuals(smryDFA, type = "tT") %>% filter(name == "model")

resids %>% group_by(.rownames) %>% summarize(sosRes = sum(.resids^2, na.rm = TRUE),
                                             nObs = n() - sum(is.na(value))) %>%
  mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) %>%
  # filter(.rownames %in% c("sardRec", "anchRec", "sardLarv",
  #                   "anchLarv", "anchYoY")) %>%
  summarize(totRMSE = sum(RMSE))



# Estimated fit of Bayesian DFA -------------------------------------------

# Rotate trends
trendsRot <- rotate_trends(fit1trend)

# Review significance of loadings and trends
dfa_loadings(trendsRot, names = datNames) %>% filter(#trend == "Trend 2",
                                                      # name %in% c("sardRec", "sardLarv",
                                                      #             "anchRec", "anchLarv",
                                                      #             "anchYoY"))#,
                                                    prob_diff0 <= 0.65)




# Calculate RMSE for each indicator
dfa_fitted(fit1trend, names = datNames) %>% 
  group_by(ID) %>% summarize(sosRes = sum((y - estimate)^2, na.rm = TRUE),
                             nObs = n() - sum(is.na(y))) %>%
  mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) %>%
  filter(ID %in% c("sardRec", "sardLarv", "anchRec", "anchLarv", "anchYoY")) %>%
  summarize(totRMSE = sum(RMSE))

# Calculate Leave-Future-Out Cross-validation scores

lfoxv1trend <- dfa_cv(fit1trend, cv_method = "lfocv", 
                      iter = fit1trend$sampling_args$iter,
                      chains = fit1trend$sampling_args$chains,
                      thin = fit1trend$sampling_args$thin,
                      cores = parallel::detectCores()-2)

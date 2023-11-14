# Summary information for MLE and Bayesian DFA loadings and fits to data
# Created: 11/8/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# Estimated fit of MLE DFA ------------------------------------------------

# Look at factor loadings
# get the inverse of the rotation matrix 
Z.est <- coef(overallDFA, type = "matrix")$Z 
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(overallDFA, type = "matrix")$Z)$rotmat
} 

# rotate factor loadings 
Z.rot <- Z.est %*% H.inv 
# rotate trends 
trends.rot <- solve(H.inv) %*% overallDFA$states

# Add CIs to marssMLE object 
overallDFA <- MARSSparamCIs(overallDFA) 
# Use coef() to get the upper and lower CIs 
Z.low <- coef(overallDFA, type = "Z", what = "par.lowCI") 
Z.up <- coef(overallDFA, type = "Z", what = "par.upCI") 
Z.rot.up <- Z.up %*% H.inv 
Z.rot.low <- Z.low %*% H.inv 
df <- data.frame(ind = rownames(Z.rot),
                 trend = rep(1:ncol(Z.rot), each = nrow(Z.rot)),
                 est = as.vector(Z.rot), 
                 conf.up = as.vector(Z.rot.up), 
                 conf.low = as.vector(Z.rot.low))

# Review significance of loadings and trends
df <- df %>% mutate(diff0 = case_when(est > 0 & conf.low > 0 ~ TRUE,
                                      est < 0 & conf.up < 0 ~ TRUE,
                                      .default = FALSE))

df %>% filter(ind %in% c("sardRec", "anchRec", "sardLarv", 
                         "anchLarv", "anchYoY"))

# indicators with no significant loadings
df %>% group_by(ind) %>% summarize(anyDiff0 = sum(diff0)) %>% filter(anyDiff0 == 0)

# Calculate RMSE for each indicator on estimates conditioned on all the data (tT)
resids <- residuals(overallDFA, type = "tT") %>% filter(name == "model")

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

# Test ability to fit many data variables

library(bayesdfa)
chains <- 3
iter <- 10

simDat20 <- sim_dfa(
              num_trends = 1,
              num_years = 15,
              num_ts = 20
            )
dim(simDat20$y_sim)

fit1trend <- fit_dfa(
  y = simDat20$y_sim, num_trends = 1, scale="zscore",
  varIndx = 1:nrow(simDat20$y_sim), # mapping of unique R matrix variances to indicators
  iter = iter, chains = chains, thin = 1)

simDat35 <- sim_dfa(
  num_trends = 1,
  num_years = 15,
  num_ts = 35
)
dim(simDat35$y_sim)

fit1trend <- fit_dfa(
  y = simDat35$y_sim, num_trends = 1, scale="zscore",
  varIndx = 1:nrow(simDat35$y_sim), # mapping of unique R matrix variances to indicators
  iter = iter, chains = chains, thin = 1)

simDat40 <- sim_dfa(
  num_trends = 1,
  num_years = 40, # initialization eventually worked with this sample size (num_years = 10000)
  num_ts = 40
)
dim(simDat40$y_sim)

fit1trend <- fit_dfa(
  y = simDat40$y_sim, num_trends = 1, scale="zscore",
  varIndx = 1:nrow(simDat40$y_sim), # mapping of unique R matrix variances to indicators
  iter = iter, chains = chains, thin = 1)

simDat41 <- sim_dfa(
  num_trends = 1,
  num_years = 100,
  num_ts = 41
)
dim(simDat41$y_sim)

fit1trend <- fit_dfa(
  y = simDat41$y_sim, num_trends = 1, scale="zscore",
  varIndx = 1:nrow(simDat41$y_sim), # mapping of unique R matrix variances to indicators
  iter = iter, chains = chains, thin = 1)

simDat50 <- sim_dfa(
  num_trends = 1,
  num_years = 100,
  num_ts = 50
)
dim(simDat50$y_sim)

fit1trend <- fit_dfa(
  y = simDat50$y_sim, num_trends = 1, scale="zscore",
  varIndx = 1:nrow(simDat50$y_sim), # mapping of unique R matrix variances to indicators
  iter = iter, chains = chains, thin = 1)

bs1trend <- fit_dfa(
  y = simDat50$y_sim, num_trends = 1, scale="zscore",
  trend_model = "bs", n_knots = 7,
  # varIndx = 1:nrow(simDat50$y_sim), # mapping of unique R matrix variances to indicators
  iter = iter, chains = chains, thin = 1)

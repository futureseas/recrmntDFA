# Leave-Future-Out Crossvalidation for MARSS DFA model selection
# Created: 12/13/2023, Robert Wildermuth

LFOXV <- function(dfaDat, # data matrix formatted for MARSS input (variables in rows, time in columns)
                  Rstructure, # either MARSS model 'R' matrix character arguments or custom matrix format
                  mTrends, # number of trends as MARSS model 'm' argument
                  peels, # number of years to peel from end of time series
                  colsRMSE = c("sardRec", "anchRec") # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
){
  LFOIC <- 0
  
  for(i in 0:peels){
    
    cat("\nFitting peel: ", i)
    
    peelDFA <- MARSS(y = dfaDat[,1:(ncol(dfaDat)-i)], 
                     form = "dfa",
                     method = "BFGS",
                     inits = list(x0 = matrix(1, 1, 1)),
                     z.score = TRUE,
                     model = list(R = Rstructure,
                                  m = mTrends),
                     silent = TRUE)
    
    # check convergence
    if(peelDFA$convergence != 0){
      return(-999) # default RMSE for failed convergence
    }
    
    # calculate RMSE
    resids <- residuals(peelDFA, type = "tT") %>% filter(name == "model")
    
    peelRMSE <- resids %>% group_by(.rownames) %>% 
                  summarize(sosRes = sum(.resids^2, na.rm = TRUE),
                            nObs = n() - sum(is.na(value))) %>%
                  mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) %>%
                  filter(.rownames %in% colsRMSE) %>%
                  summarize(totRMSE = sum(RMSE))
    
    LFOIC <- LFOIC + peelRMSE$totRMSE
  }
  
  return(LFOIC)
}
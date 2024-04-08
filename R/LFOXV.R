# Leave-Future-Out Crossvalidation for MARSS DFA model selection
# Created: 12/13/2023, Robert Wildermuth

# Informed Innovation Residuals -------------------------------------------

# fxn for calculating informed innovation residuals
infInnoResid <- function(objMARSS, # MARSS model fit
                         colsRMSE = c("sardRec", "anchRec") # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
                         ){
  
  
  origDat <- objMARSS$call$data
  # need to scale data for residual calcs and 'newdata' input
  origDat <- apply(origDat, 1, scale, simplify = TRUE) %>% t() 
  
  # create 'yStar' by removing last observation of 'colsRMSE' variables
  yStar <- origDat
  yStar[colsRMSE, ncol(yStar)] <- NA
  
  yExpct <- predict(object = objMARSS, n.ahead = 0, interval = "none", 
                    newdata = list(t = 1:ncol(objMARSS$call$data),
                                   y = yStar),
                    type = 'ytT', x0 = "use.model")
  
  # check
  # yExpct$pred %>% filter(.rownames == "anchRec") %>% 
  #   mutate(fullDat = origDat["anchRec",], newdat = yStar["anchRec",])
  
  # transpose for column addition
  trCol <- origDat %>% t() 
  trCol <- as.data.frame(trCol) %>% mutate(t = 1:nrow(trCol)) %>% 
              pivot_longer(cols = -t, names_to = ".rownames", 
                           values_to = "origDat")
  
  infInnoResids <- yExpct$pred %>% full_join(y = trCol, by = c("t", ".rownames"))
  
  # infInnoResids %>% filter(y != origDat | (is.na(y) & !is.na(origDat)))
  # infInnoResids %>% filter(.rownames %in% colsRMSE)
  
  infInnoResids <- infInnoResids %>% mutate(infInnoResid = origDat - estimate)
  
  return(infInnoResids)
}

# Projection Residuals -------------------------------------------

# fxn for calculating projection residuals
ProjectResid <- function(objMARSS, # MARSS model fit
                         colsRMSE = c("sardRec", "anchRec") # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
){
  
  
  origDat <- objMARSS$call$data
  # need to scale data for residual calcs and 'newdata' input
  origDat <- apply(origDat, 1, scale, simplify = TRUE) %>% t() 
  
  # create 'yStar' by removing all observations of 'colsRMSE' variables
  yStar <- origDat
  yStar[colsRMSE, ] <- NA
  
  yExpct <- predict(object = objMARSS, n.ahead = 0, interval = "none", 
                    newdata = list(t = 1:ncol(objMARSS$call$data),
                                   y = yStar),
                    type = 'ytT', x0 = "use.model")
  
  # check
  # yExpct$pred %>% filter(.rownames == "anchRec") %>% 
  #   mutate(fullDat = origDat["anchRec",], newdat = yStar["anchRec",])
  
  # transpose for column addition
  trCol <- origDat %>% t() 
  trCol <- as.data.frame(trCol) %>% mutate(t = 1:nrow(trCol)) %>% 
    pivot_longer(cols = -t, names_to = ".rownames", 
                 values_to = "origDat")
  
  projResids <- yExpct$pred %>% full_join(y = trCol, by = c("t", ".rownames"))
  
  # projResids %>% filter(y != origDat | (is.na(y) & !is.na(origDat)))
  # projResids %>% filter(.rownames %in% colsRMSE)
  
  projResids <- projResids %>% mutate(projResid = origDat - estimate)
  
  return(projResids)
}


# Leave-Future-Out Cross-Validation ---------------------------------------

LFOXV <- function(dfaDat, # data matrix formatted for MARSS input (variables in rows, time in columns)
                  Rstructure, # either MARSS model 'R' matrix character arguments or custom matrix format
                  mTrends, # number of trends as MARSS model 'm' argument
                  peels, # number of years to peel from end of time series
                  colsRMSE = c("sardRec", "anchRec") # column names to calculate RMSE for. Default: c("sardRec", "anchRec")
){
  
  peelRMSE <- tibble(.rownames = '', 
                     t = 0)[0,]
  
  
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
    smResids <- residuals(peelDFA, type = "tT") %>% filter(name == "model") # smoothation resids (Y|data from 1:T)
    innoResids <- residuals(peelDFA, type = "tt1") %>% filter(name == "model") # innovation resids (Y|data from 1:(t-1)) - one-step-ahead resids
    contResids <- residuals(peelDFA, type = "tt") %>% filter(name == "model") # contemporaneous resids (Y|data from 1:t)
    
    infInnoRes <- infInnoResid(objMARSS = peelDFA, colsRMSE = colsRMSE)
    projRes <- ProjectResid(objMARSS = peelDFA, colsRMSE = colsRMSE)
    
    # combine and compare values
    smResids <- smResids %>% select(.rownames, t, value, .fitted, .resids) %>% 
                  rename_with(~paste0( .x, "_sm"), c(value, .fitted, .resids))
    innoResids <- innoResids %>% select(.rownames, t, value, .fitted, .resids) %>% 
                    rename_with(~paste0( .x, "_inno"), c(value, .fitted, .resids))
    contResids <- contResids %>% select(.rownames, t, value, .fitted, .resids) %>% 
                    rename_with(~paste0( .x, "_cont"), c(value, .fitted, .resids))
    
    resids <- infInnoRes %>% full_join(y = smResids, by = c("t", ".rownames")) %>%
                full_join(y = innoResids, by = c("t", ".rownames")) %>%
                full_join(y = contResids, by = c("t", ".rownames")) %>%
                full_join(y = projRes, by = c("t", ".rownames"))
      
    # # Checks
    # resids %>% select(.rownames, t, y, origDat, value_sm, value_inno, value_cont) %>%
    #   mutate(checkInfInno = round(y, digits = 10) == round(origDat, digits = 10),
    #          checkInno = round(origDat, digits = 10) == round(value_inno, digits = 10),
    #          checkSm = round(origDat, digits = 10) == round(value_sm, digits = 10),
    #          checkCont = round(origDat, digits = 10) == round(value_cont, digits = 10)) %>%
    #   # filter(isFALSE(checkInfInno)| is.na(checkCont))
    #   #        #isFALSE(checkInno))
    #   #        #isFALSE(checkSm))
    #   #        isFALSE(checkCont) | is.na(checkCont))
    #   filter(is.na(y))
    # resids %>% select(-c(y, origDat, value_sm, value_inno, value_cont)) %>%
    #   filter(round(estimate, digits = 5) != round(.fitted_sm, digits = 5)) # many estimates not equal across time series
    # resids %>% select(.rownames, t, infInnoResid, .resids_sm, .resids_inno, .resids_cont) %>%
    #   #filter(.rownames %in% colsRMSE) %>% 
    #   pivot_longer(cols = c(infInnoResid, .resids_sm, .resids_inno, .resids_cont)) %>%
    #   ggplot(aes(x = t, y = value, colour = name)) +
    #   geom_line() +
    #   facet_wrap(~.rownames)
    # # diffs between informed innovation and naive innovation resids greater than 
    # # btwn informed innovation and contemporaneous resids
    # # informed innovations and smoothations basically the same (within a peel)
    # resids %>% select(.rownames, t, infInnoResid, .resids_sm, .resids_inno, .resids_cont) %>%
    #   filter(t == ncol(dfaDat))
    
    # collect residuals for datum of interest
    peelRMSE <- resids %>% select(.rownames, t, infInnoResid, .resids_sm, 
                                  .resids_inno, .resids_cont, projResid) %>% 
                  filter(.rownames %in% colsRMSE, t == max(t)) %>%
                  mutate(peel = i) %>%
                  bind_rows(peelRMSE)
  } # end peel for-loop
    
  overPeelthenVar <- peelRMSE %>% 
                  pivot_longer(cols = c(infInnoResid, .resids_sm, projResid,
                                        .resids_inno, .resids_cont),
                               values_to = "resids", names_to = "resType") %>%
                  # Calc RMSE over all time points
                  group_by(.rownames, resType) %>%
                  summarize(sosRes = sum(resids^2, na.rm = TRUE),
                            nObs = sum(!is.na(resids)),
                            avgAbsRes = mean(abs(resids), na.rm = TRUE)) %>%
                  mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) %>%
                  group_by(resType) %>%
                  summarize(totRMSE = sum(RMSE),
                            meanRMSE = mean(RMSE))
  # Over variables then peels
  overVarthenPeel <- peelRMSE %>% 
                      pivot_longer(cols = c(infInnoResid, .resids_sm, projResid,
                                            .resids_inno, .resids_cont),
                                   values_to = "resids", names_to = "resType") %>%
                      # Calc RMSE over all time points
                      group_by(peel, resType) %>%
                      summarize(sosRes = sum(resids^2, na.rm = TRUE),
                                nObs = sum(!is.na(resids)),
                                avgAbsRes = mean(abs(resids), na.rm = TRUE)) %>%
                      mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) %>%
                      group_by(resType) %>%
                      summarize(totRMSE = sum(RMSE),
                                meanRMSE = mean(RMSE))
   # OR can calculate over all vars of interest AND peel
  overPeelandVar <- peelRMSE %>% 
                      pivot_longer(cols = c(infInnoResid, .resids_sm, projResid,
                                            .resids_inno, .resids_cont),
                                   values_to = "resids", names_to = "resType") %>%
                      # Calc RMSE over all time points
                      group_by(resType) %>%
                      summarize(sosRes = sum(resids^2, na.rm = TRUE),
                                nObs = sum(!is.na(resids)),
                                avgAbsRes = mean(abs(resids), na.rm = TRUE)) %>%
                      mutate(RMSE = sqrt(sosRes/nObs)) %>% arrange(RMSE) 
    
  LFOIC <- overPeelandVar %>% filter(resType == "infInnoResid") %>% pull(RMSE)

  
  return(LFOIC)
}
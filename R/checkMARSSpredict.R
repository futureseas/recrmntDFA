# Check how predict() function works for projection with forced variables
# Created: 3/25/2024, Robert Wildermuth

library(tidyverse)
library(MARSS)
source("R/LFOXV.R")

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

# LFOIC for historical projection dataset ---------------------------------------

# remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
allDat <- datDFA%>% select(c("year", "HCI_R3", "HCI_R4", "BEUTI_33N", "BEUTI_39N", 
                             "CUTI_33N", "CUTI_39N",
                             "OC_LUSI_33N", "OC_LUSI_36N", "OC_LUSI_39N", "OC_STI_33N",
                             "OC_STI_36N", "OC_STI_39N", "ZM_NorCal", "ZM_SoCal",
                             "sardSpawnHab",
                             "anchSpawnHab", "daysAbove5pct", "daysAbove40pct",
                             "sardNurseHab", "anchNurseHab",
                             "springSST", "summerSST", "avgNearTransspring",
                             "avgNearTranssummer", "avgOffTransspring", "avgOffTranssummer",
                             # Variables of interest
                             "sardRec", "anchRec")) 

datNames <- names(allDat)[-1]

# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- matrix(list(0),length(datNames),length(datNames)) 
diag(Rcustom) <- c("HCI", "HCI",
                   "BEUTI", "BEUTI", 
                   "CUTI", "CUTI",
                   "LUSI", "LUSI", "LUSI", 
                   "STI", "STI", "STI",
                   "NEMURO", "NEMURO",
                   "SDM", "SDM", "TIME", "TIME", "SDM", "SDM",
                   "SST", "SST",
                   "Transp", "Transp", "Transp", "Transp",
                   "sardRec", "anchRec")

Rstructure <- "diagonal and equal" 
mTrends <- 3

# transpose for MARSS formatting
chkDat <- allDat %>% filter(year %in% 1990:2019) %>% select(-year) %>% t()
# need to scale data for residual calcs and 'newdata' input
sclDat <- apply(chkDat, 1, scale, simplify = TRUE) %>% t()

fitDFA <- MARSS(y = chkDat, 
                 form = "dfa",
                 method = "BFGS",
                 inits = list(x0 = matrix(1, 1, 1)),
                 z.score = TRUE,
                 model = list(R = Rstructure,
                              m = mTrends),
                 silent = TRUE)

# Try with last 5 years removed

last5Dat <- chkDat
last5Dat[,26:30] <- NA

last5DFA <- MARSS(y = last5Dat, 
                  form = "dfa",
                  method = "BFGS",
                  inits = list(x0 = matrix(1, 1, 1)),
                  z.score = TRUE,
                  model = list(R = Rstructure,
                               m = mTrends),
                  silent = TRUE)

noEmptyDFA <- MARSS(y = last5Dat[,1:25], 
                  form = "dfa",
                  method = "BFGS",
                  inits = list(x0 = matrix(1, 1, 1)),
                  z.score = TRUE,
                  model = list(R = Rstructure,
                               m = mTrends),
                  silent = TRUE)

# compare Z and x estimates
zs <- data.frame(fullDatZs = fitDFA$par$Z,
                 last5DatZs = last5DFA$par$Z,
                 noEmptyDatZs = noEmptyDFA$par$Z)

fullX <- bind_cols(t(fitDFA$states), t = 1:nrow(t(fitDFA$states)))
last5X <- bind_cols(t(last5DFA$states), t = 1:nrow(t(last5DFA$states)))
noEmptyX <- bind_cols(t(noEmptyDFA$states), t = 1:nrow(t(noEmptyDFA$states)))
xs <- full_join(fullX, last5X, by = "t") %>% full_join(noEmptyX, by = "t")
tail(xs, n = 10)

test1 <- fitted(noEmptyDFA, type = "ytT", interval = "none")
test2 <- fitted(last5DFA, type = "ytT", interval = "none")
yFitted <- full_join(test1, test2, by = c(".rownames", "t"))
yFitted %>% filter(t >24)

test3 <- fitted(noEmptyDFA, type = "xtT", interval = "none")
test4 <- fitted(last5DFA, type = "xtT", interval = "none")
xFitted <- full_join(test3, test4, by = c(".rownames", "t"))
xFitted %>% filter(t >24)

test5 <- predict(noEmptyDFA, type = "xtT", interval = "none",
                 x0 = "use.model", n.ahead = 0)
test6 <- predict(last5DFA, type = "xtT", interval = "none",
                 x0 = "use.model", n.ahead = 0)
xPredict <- full_join(test5$pred, test6$pred, by = c(".rownames", "t"))
xPredict %>% filter(t >24) # same as fitted() output

# try predict() with the full time series as newdata
fullPredict <- predict(last5DFA, type = "xtT", interval = "none",
                       x0 = "use.model", n.ahead = 0,
                       newdata = list(t = 1:ncol(sclDat),
                                      y = sclDat))
zs$fullPredZs <- fullPredict$model$par$Z

fullPredXs <- fullPredict$pred %>% pivot_wider(names_from = ".rownames",
                                               values_from = c(".x", "estimate"))
xs <- full_join(xs, fullPredXs, by = "t") # '.x' from 'fullPredict' within ~same order of magnitude?

fullPredYs <- predict(last5DFA, type = "ytT", interval = "none",
                       x0 = "use.model", n.ahead = 0,
                       newdata = list(t = 1:ncol(sclDat),
                                      y = sclDat))
fullPredYs
# transpose for column addition
trCol <- sclDat %>% t() 
trCol <- as.data.frame(trCol) %>% mutate(t = 1:nrow(trCol)) %>% 
  pivot_longer(cols = -t, names_to = ".rownames", 
               values_to = "fullDat")
yPredict <- fullPredYs$pred %>% full_join(trCol, by = c(".rownames", "t"))

# try forecast() with the full time series as newdata
fullForecast <- forecast(noEmptyDFA, type = "xtT", interval = "none",
                         h = 5, # number of intervals to forecast
                         newdata = list(y = sclDat[, 26:30]))
xPredFore <- full_join(fullPredict$pred, fullForecast$pred, by = c(".rownames", "t"))

fullForeY <- forecast(noEmptyDFA, type = "ytT", interval = "none",
                         h = 5, # number of intervals to forecast
                         newdata = list(y = sclDat[, 26:30]))
fullForeY$pred %>% filter(.rownames == "avgOffTranssummer")
sclDat[26,]
noEmptyDFA$marss$data[26,]
chkDat[26,]
###################################################
# compare last time interval
colsRMSE <- c("sardRec", "anchRec")
infInnoRes <- infInnoResid(objMARSS = fitDFA, colsRMSE = colsRMSE)
# expected value 'estimate' calculated from predict() 
infInnoRes %>% filter(.rownames %in% colsRMSE, t == max(t))

# reconstruct Z matrix
Z <- matrix(0, nrow = 28, ncol = 3)
Z[,1] <- fitDFA$par$Z[1:28]
Z[2:28,2] <- fitDFA$par$Z[29:(28+27)]
Z[3:28,3] <- fitDFA$par$Z[(28+27+1):(28+27+26)]

y29 <- sclDat[, 30]
y29[is.na(y29)] <- 0 # expect missing data is 0
y29 %*% Z

# compare future intervals

# read prepped dataset
projDat <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/formattedDFAprojDat.csv")

datNames <- names(projDat)[-(1:2)]

# transpose for MARSS formatting
projDatIPSL <- projDat %>% filter(ESM == "IPSL") %>%
  select(-year, -ESM) %>%
  select(rownames(projectDFA$model$data)) %>%
  t()

forecastIPSL <- predict(object = projectDFA,
                        # n.ahead = 81, # alone, this sets value of X as last state of fitted model
                        # newdata = list(t = 41:121,
                        #                y = testProj),
                        newdata = list(t = 41:121,
                                       y = projDatIPSL), # data are zscored in FormatProjectionData.R
                        x0 = "use.model",
                        type = "ytt")
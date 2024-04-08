# Project DFA
# Created: 11/13/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# read prepped dataset
projDat <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/formattedDFAprojDat.csv")

datNames <- names(projDat)[-(1:2)]

# load fitted MARSS output
# load(file = "marssFit_1980to2019_noBio_3trend_Rcustom.RData")
load(file = "marssFit_1980to2019_ProjDFA_4trend_Rcustom.RData")

# transpose for MARSS formatting
projDatIPSL <- projDat %>% filter(ESM == "IPSL") %>%
                  select(-year, -ESM) %>%
                  select(c("BEUTI_33", "BEUTI_39", "CUTI_33", "CUTI_39", 
                           "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
                           "STI_36", "STI_39", "avgSSWIspring", "avgSSWIsummer", 
                           "sardLarv", "anchLarv", "anchYoY", "ZM_NorCal",
                           "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", 
                           "daysAbove5pct", "daysAbove40pct", "sardNurseHab", 
                           "anchNurseHab", "anchRec", "sardRec", "springSST", 
                           "summerSST", "avgNearTransspring", "avgNearTranssummer",
                           "avgOffTransspring", "avgOffTranssummer")) %>%
                  t()

# Testing what predict() does for estimating latent states
# testProj <- matrix(0, nrow = nrow(projDatIPSL), ncol = 81) # makes projected values of X=0 also
# testProj <- matrix(NA, nrow = nrow(projDatIPSL), ncol = 81) # makes projected values of X=0 also
# testProj[which(datNames %in% c("HCI", "anchRec")), ] <- 0 # sets most precise, highest loading variable to mean
# testProj[which(datNames %in% c("meanSSBwt", "naupBio")), ] <- 0 # sets least precise, lowest loading variable to mean

forecastIPSL <- forecast(object = projectDFA,
                        # n.ahead = 81, # alone, this sets value of X as last state of fitted model
                        # newdata = list(t = 41:121,
                        #                y = testProj),
                        newdata = list(y = projDatIPSL), # data are zscored in FormatProjectionData.R
                        h = 81,
                        interval = "none",
                        type = "ytt")
                        # type = "xtT") #
                      # RW!: do we want projections with full knowledge (ytT), or just to "present" (ytt)?
plot(forecastIPSL)


projDatGFDL <- projDat %>% filter(ESM == "GFDL") %>%
                  select(-year, -ESM) %>%
                  select(c("BEUTI_33", "BEUTI_39", "CUTI_33", "CUTI_39", 
                           "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
                           "STI_36", "STI_39", "avgSSWIspring", "avgSSWIsummer", 
                           "sardLarv", "anchLarv", "anchYoY", "ZM_NorCal",
                           "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", 
                           "daysAbove5pct", "daysAbove40pct", "sardNurseHab", 
                           "anchNurseHab", "anchRec", "sardRec", "springSST", 
                           "summerSST", "avgNearTransspring", "avgNearTranssummer",
                           "avgOffTransspring", "avgOffTranssummer")) %>% t()

forecastGFDL <- forecast(object = projectDFA,
                        newdata = list(y = projDatGFDL), # data are zscored in FormatProjectionData.R
                        h = 81,
                        interval = "confidence",
                        type = "ytt")
plot(forecastGFDL)


projDatHAD <- projDat %>% filter(ESM == "HAD") %>%
                select(-year, -ESM) %>% 
                select(c("BEUTI_33", "BEUTI_39", "CUTI_33", "CUTI_39", 
                         "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
                         "STI_36", "STI_39", "avgSSWIspring", "avgSSWIsummer", 
                         "sardLarv", "anchLarv", "anchYoY", "ZM_NorCal",
                         "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", 
                         "daysAbove5pct", "daysAbove40pct", "sardNurseHab", 
                         "anchNurseHab", "anchRec", "sardRec", "springSST", 
                         "summerSST", "avgNearTransspring", "avgNearTranssummer",
                         "avgOffTransspring", "avgOffTranssummer")) %>% t()

forecastHAD <- forecast(object = projectDFA,
                       newdata = list(y = projDatHAD), # data are zscored in FormatProjectionData.R
                       h = 81,
                       interval = "none",
                       type = "ytt")
plot(forecastHAD)

# forecastHAD
# forecastGFDL
forecastIPSL$pred %>% #filter(!is.na(y)) %>%
  ggplot(aes(x = t)) +
  geom_vline(xintercept = 40, color = "grey") +
  geom_point(aes(y = y), color = "darkblue") +
  geom_line(aes(y = estimate), linewidth = 1) +
  geom_hline(yintercept = 0) +
  facet_wrap(~.rownames, scales = "free_y") +
  theme_classic()

forecastHAD$pred %>% filter(.rownames %in% c("sardRec", "anchRec", 
                                             "sardLarv", "anchLarv")) %>%
  ggplot(aes(x = t)) +
  geom_point(aes(y = y), color = "darkblue") +
  geom_line(aes(y = estimate), linewidth = 1) +
  geom_hline(yintercept = 0) +
  facet_wrap(~.rownames, scales = "free_y") +
  theme_classic()

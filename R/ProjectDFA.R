# Project DFA
# Created: 11/13/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# read prepped dataset
projDat <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/formattedDFAprojDat.csv")

datNames <- names(projDat)[-(1:2)]

# load fitted MARSS output
load(file = "marssFit_1980to2019_noBio_3trend_Rcustom.RData")

# transpose for MARSS formatting
projDatIPSL <- projDat %>% filter(ESM == "IPSL") %>%
              select(-year, -ESM) %>% t()

forecastIPSL <- predict(object = overallDFA,
                        newdata = list(t = 41:121,
                                       y = projDatIPSL), # data are zscored in FormatProjectionData.R
                        x0 = "use.model",
                        type = "ytt")# type = "xtT") # 
                      # RW!: do we want projections with full knowledge (ytT), or just to "present" (ytt)?
plot(forecastIPSL)


projDatGFDL <- projDat %>% filter(ESM == "GFDL") %>%
  select(-year, -ESM) %>% t()

forecastGFDL <- predict(object = overallDFA,
                        newdata = list(t = 41:121,
                                       y = projDatGFDL), # data are zscored in FormatProjectionData.R
                        x0 = "use.model",
                        type = "ytt")# type = "xtT") #
plot(forecastGFDL)


projDatHAD <- projDat %>% filter(ESM == "HAD") %>%
  select(-year, -ESM) %>% t()

forecastHAD <- predict(object = overallDFA,
                        newdata = list(t = 41:121,
                                       y = projDatHAD), # data are zscored in FormatProjectionData.R
                        x0 = "use.model",
                        type = "ytt")# type = "xtT") #
plot(forecastHAD)
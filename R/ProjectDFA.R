# Project DFA
# Created: 11/13/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# read prepped dataset
projDat <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/formattedDFAprojDat.csv")

datNames <- names(projDat)[-(1:2)]

# load fitted MARSS output
load(file = "marssFit_1990to2019_ProjDFA_5trend_Rcustom.RData")

# transpose for MARSS formatting
projDatIPSL <- projDat %>% filter(ESM == "IPSL") %>%
                  select(-year, -ESM) %>%
                  select(c("HCI_R3", "HCI_R4", "BEUTI_33", "BEUTI_39", "CUTI_33", 
                           "CUTI_39", "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
                           "STI_36", "STI_39", "sardLarv", "anchLarv", "anchYoY", 
                           "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", 
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
                        interval = "confidence",
                        type = "ytt")
fcastXsIPSL <- forecast(object = projectDFA,
                        newdata = list(y = projDatIPSL), # data are zscored in FormatProjectionData.R
                        h = 81,
                        interval = "confidence",
                        type = "xtt")
                        # type = "xtT") #
                      # RW!: do we want projections with full knowledge (ytT), or just to "present" (ytt)?
plot(forecastIPSL)


projDatGFDL <- projDat %>% filter(ESM == "GFDL") %>%
                  select(-year, -ESM) %>%
                  select(c("HCI_R3", "HCI_R4", "BEUTI_33", "BEUTI_39", "CUTI_33", 
                           "CUTI_39", "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
                           "STI_36", "STI_39", "sardLarv", "anchLarv", "anchYoY", 
                           "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", 
                           "daysAbove5pct", "daysAbove40pct", "sardNurseHab", 
                           "anchNurseHab", "anchRec", "sardRec", "springSST", 
                           "summerSST", "avgNearTransspring", "avgNearTranssummer",
                           "avgOffTransspring", "avgOffTranssummer")) %>% t()

forecastGFDL <- forecast(object = projectDFA,
                        newdata = list(y = projDatGFDL), # data are zscored in FormatProjectionData.R
                        h = 81,
                        interval = "confidence",
                        type = "ytt")
fcastXsGFDL <- forecast(object = projectDFA,
                        newdata = list(y = projDatGFDL), # data are zscored in FormatProjectionData.R
                        h = 81,
                        interval = "confidence",
                        type = "xtt")
plot(forecastGFDL)


projDatHAD <- projDat %>% filter(ESM == "HAD") %>%
                select(-year, -ESM) %>% 
                select(c("HCI_R3", "HCI_R4", "BEUTI_33", "BEUTI_39", "CUTI_33", 
                         "CUTI_39", "LUSI_33", "LUSI_36", "LUSI_39", "STI_33", 
                         "STI_36", "STI_39", "sardLarv", "anchLarv", "anchYoY", 
                         "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", 
                         "daysAbove5pct", "daysAbove40pct", "sardNurseHab", 
                         "anchNurseHab", "anchRec", "sardRec", "springSST", 
                         "summerSST", "avgNearTransspring", "avgNearTranssummer",
                         "avgOffTransspring", "avgOffTranssummer")) %>% t()

forecastHAD <- forecast(object = projectDFA,
                       newdata = list(y = projDatHAD), # data are zscored in FormatProjectionData.R
                       h = 81,
                       interval = "confidence",
                       type = "ytt")
fcastXsHAD <- forecast(object = projectDFA,
                       newdata = list(y = projDatHAD), # data are zscored in FormatProjectionData.R
                       h = 81,
                       interval = "confidence",
                       type = "xtt")
plot(forecastHAD)

# forecastHAD
# forecastGFDL
forecastIPSL$pred %>% #filter(!is.na(y)) %>%
  ggplot(aes(x = t)) +
  geom_vline(xintercept = 30, color = "grey") +
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


# Plots for MS ------------------------------------------------------------

# projection time series
projTSGFDL <- bind_rows(forecastGFDL$pred, fcastXsGFDL$pred) %>% mutate(ESM = "GFDL")
names(projTSGFDL) <- make.names(names(projTSGFDL))
projTSGFDL %>% 
  # filter(.rownames %in% c("HCI_R4", "BEUTI_39N", "OC_STI_33N", 
  #                        "ZM_NorCal", "anchSpawnHab", "daysAbove40pct",
  #                        "springSST", "avgNearTransspring",
  #                        "X2", "X3", "X4", "X5",
  #                        "sardRec", "anchRec")) %>%
  # mutate(.rownames = factor(.rownames, levels = c("anchSpawnHab", "daysAbove40pct",
  #                                                 "HCI_R4", "BEUTI_39N", 
  #                                                 "OC_STI_33N", "avgNearTransspring",
  #                                                 "springSST", "ZM_NorCal", 
  #                                                 "X2", "X3", "X4", "X5",
  #                                                 "anchRec", "sardRec"))) %>%

  # filter(.rownames %in% paste0("X", 1:5)) %>%

  filter(.rownames %in% c("HCI_R3", "HCI_R4", "BEUTI_33N", "BEUTI_39N", "CUTI_33N",
                         "CUTI_39N", "OC_LUSI_33N", "OC_LUSI_36N", "OC_LUSI_39N",
                         "OC_STI_33N", "OC_STI_36N", "OC_STI_39N",
                         "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab",
                         "daysAbove5pct", "daysAbove40pct", "sardNurseHab",
                         "anchNurseHab",  "springSST",
                         "summerSST", "avgNearTransspring", "avgNearTranssummer",
                         "avgOffTransspring", "avgOffTranssummer")) %>%

  ggplot(aes(x = t, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lo.95, ymax = Hi.95), alpha = 0.3) +
  facet_wrap(~.rownames, scales = "free", ncol = 4) +
  geom_point(aes(y = y), color = "steelblue", size = 0.2) +
  geom_vline(xintercept = 30.5) + geom_hline(yintercept = 0) +
  theme_classic() +
  labs(title = "GFDL")

projTSGFDL %>% 
  filter(.rownames %in% c("sardLarv", "anchLarv", "anchYoY", "anchRec", "sardRec" )) %>%
  ggplot(aes(x = t, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lo.95, ymax = Hi.95), alpha = 0.3) +
  facet_wrap(~.rownames, scales = "free") +
  geom_point(aes(y = y), color = "steelblue", shape = 1, size = 1) +
  geom_vline(xintercept = 30.5) +
  theme_classic() +
  labs(title = "GFDL")

# Histogram of states across all ESMs
projTSIPSL <- bind_rows(forecastIPSL$pred, fcastXsIPSL$pred) %>% mutate(ESM = "IPSL")
names(projTSIPSL) <- make.names(names(projTSIPSL))

projTSHAD <- bind_rows(forecastHAD$pred, fcastXsHAD$pred) %>% mutate(ESM = "HAD")
names(projTSHAD) <- make.names(names(projTSHAD))

# projTSIPSL %>% 
projTSHAD %>% 
  # filter(.rownames %in% paste0("X", 1:5)) %>%
  
  filter(.rownames %in% c("HCI_R3", "HCI_R4", "BEUTI_33N", "BEUTI_39N", "CUTI_33N",
                          "CUTI_39N", "OC_LUSI_33N", "OC_LUSI_36N", "OC_LUSI_39N",
                          "OC_STI_33N", "OC_STI_36N", "OC_STI_39N",
                          "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab",
                          "daysAbove5pct", "daysAbove40pct", "sardNurseHab",
                          "anchNurseHab",  "springSST",
                          "summerSST", "avgNearTransspring", "avgNearTranssummer",
                          "avgOffTransspring", "avgOffTranssummer")) %>%
  #   
  ggplot(aes(x = t, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lo.95, ymax = Hi.95), alpha = 0.3) +
  facet_wrap(~.rownames, scales = "free", ncol = 4) +
  geom_point(aes(y = y), color = "steelblue", size = 0.2) +
  geom_vline(xintercept = 30.5) + geom_hline(yintercept = 0) +
  theme_classic() +
  labs(title = "HAD")

# projTSIPSL %>% 
projTSHAD %>% 
  filter(.rownames %in% c("sardLarv", "anchLarv", "anchYoY", "anchRec", "sardRec" )) %>%
  ggplot(aes(x = t, y = estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lo.95, ymax = Hi.95), alpha = 0.3) +
  facet_wrap(~.rownames, scales = "free") +
  geom_point(aes(y = y), color = "steelblue", shape = 1, size = 1) +
  geom_vline(xintercept = 30.5) + geom_hline(yintercept = 0) +
  theme_classic() +
  labs(title = "HAD")


allProj <- bind_rows(projTSHAD, projTSGFDL, projTSIPSL) %>% 
            mutate(Year = t + 1989,
                   periods = case_when(Year <= 2019 ~ "Historical",
                                       Year <= 2039 & Year > 2019 ~ "2020-2039",
                                       Year <= 2059 & Year > 2039 ~ "2040-2059",
                                       Year <= 2079 & Year > 2059 ~ "2060-2079",
                                       Year > 2079 ~ "2080-2100"),
                   periods = factor(periods, levels = c("Historical","2020-2039",
                                                        "2040-2059", "2060-2079", 
                                                        "2080-2100")))
# rotate factor loadings
Z.est <- coef(projectDFA, type = "matrix")$Z
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(projectDFA, type = "matrix")$Z)$rotmat
} 

Z.rot <- Z.est %*% H.inv
plotZs <- as.data.frame(Z.rot)
plotZs$.rownames <- row.names(Z.rot)
plotZs <- plotZs %>%  filter(.rownames %in% c("avgNearTranssummer", "sardNurseHab", "HCI_R4", # (nearly) significant
                                              "OC_STI_33N", "sardSpawnHab", "daysAbove5pct", # not sig but strong
                                              "X4", "sardRec")) %>% 
            mutate(.rownames = factor(.rownames, 
                                      levels = c("avgNearTranssummer", "sardNurseHab", "HCI_R4", # (nearly) significant
                                                 "OC_STI_33N", "sardSpawnHab", "daysAbove5pct", # not sig but strong
                                                 "X4", "sardRec")))

# Example plot of forcing variables, primary latent trend, and sardine recruitment projection
allProj %>%
  # full_join(y = Z.rot, by = ".rownames") %>%
  mutate(t = t+1989) %>%
  filter(.rownames %in% c("avgNearTranssummer", "sardNurseHab", "HCI_R4", # (nearly) significant
                          "OC_STI_33N", "sardSpawnHab", "daysAbove5pct", # not sig but strong
                          "X4", "sardRec")) %>%
  mutate(.rownames = factor(.rownames, levels = c("avgNearTranssummer", "sardNurseHab", "HCI_R4", # (nearly) significant
                                                  "OC_STI_33N", "sardSpawnHab", "daysAbove5pct", # not sig but strong
                                                  "X4", "sardRec"))) %>%
  ggplot(aes(x = t, y = estimate, color = ESM, fill = ESM)) +
  geom_line() +
  geom_ribbon(aes(ymin = Lo.95, ymax = Hi.95), alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~.rownames, scales = "free", nrow = 3) +
  geom_point(aes(y = y), color = "black", shape = 1, size = 1) +
  geom_vline(xintercept = 2019.5) + geom_hline(yintercept = 0) +
  geom_text(data = plotZs,
            mapping = aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
                          label = round(V4, digits = 3)), 
            inherit.aes = FALSE,
            color = "black") +
  theme_classic()



allProjSmry <- allProj %>% group_by(.rownames, ESM, periods) %>%
                  summarize(varMean = mean(estimate),
                            varMed = median(estimate),
                            varSD = sd(estimate)) %>% 
                  filter(periods != "Historical")

allProjSmry <- allProj %>% summarize(.by = c(.rownames, periods),
                                      varMean = mean(estimate),
                                      varMed = median(estimate),
                                      varSD = sd(estimate),
                                      ESM = "Total") %>%
                  bind_rows(allProjSmry)

histHisto <- allProj %>% filter(periods == "Historical", ESM == "GFDL") %>%
                mutate(ESM = "Observations",
                       val = y)

projHisto <- allProj %>% filter(periods != "Historical") %>%
                mutate(val = estimate) %>% 
                bind_rows(histHisto)

projHisto %>% 
  filter(.rownames %in% c(#"sardLarv", "anchLarv", "anchYoY", 
                          "anchRec", "sardRec" )) %>%
  ggplot() + geom_vline(xintercept = 0, color = "grey") +
  geom_histogram(aes(x = val, fill = ESM)) +
  facet_grid(rows = vars(.rownames), cols = vars(periods)) +
  theme_bw() +
  geom_vline(data = allProjSmry %>% 
                      filter(.rownames %in% c(#"sardLarv", "anchLarv", "anchYoY", 
                                              "anchRec", "sardRec"),
                             periods != "Historical"), 
             aes(xintercept = varMean, color = ESM), linewidth = 0.5) +
  xlab("Standardized Index Estimate")

projProp <- allProj %>% filter(periods == "Historical") %>%
              group_by(.rownames, ESM, periods) %>%
              summarize(varMin = min(y, na.rm = TRUE),
                        varMax = max(y, na.rm = TRUE),
                        varMean = mean(y, na.rm = TRUE)) %>%
              select(-periods) %>%
              full_join(y = allProj, by = c(".rownames", "ESM")) %>%
              mutate(inPropCat = case_when(y >= varMax ~ "histMax",
                                           y <= varMin ~ "histMin",
                                           y > 1 ~ ">1 SD",
                                           y < -1 ~ "<1 SD",
                                           y < 1 & y > -1 ~ "average"),
                     expPropCat = case_when(estimate >= varMax ~ "histMax",
                                            estimate <= varMin ~ "histMin",
                                            estimate > 1 ~ ">1 SD",
                                            estimate < -1 ~ "<1 SD",
                                            estimate < 1 & estimate > -1 ~ "average"),
                     finalCat = case_when(periods == "Historical" ~ inPropCat,
                                          periods != "Historical" ~ expPropCat),
                     finalCat = factor(finalCat, levels = c("histMax", ">1 SD",
                                                                "average", "<1 SD", 
                                                                "histMin")))


projProp %>% filter(.rownames %in% c("anchRec", "sardRec")) %>%
  select(.rownames, ESM, t, periods, finalCat) %>%
  na.omit() %>%
  ggplot(aes(x = periods, y = t, fill = finalCat)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(rows = vars(.rownames), cols = vars(ESM)) +
  theme_bw() +
  labs(y = "Recruitment Deviation Proportion", fill = "Category")

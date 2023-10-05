# Projection data

library(tidyverse)
library(lubridate)
library(sf)
library(stars)
library(cubelyr)
library(units)

# data file path
datPath <- "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/"

# From Mer Pozo Buil (6/4/2023): "monthly SST averaged for the diagonal CalCOFI region for the 3 datasets:
# ROMS-ESMS: includes the 3 projections (ROMS-GFDL, ROMS-IPSL, and ROMS-HAD) from 1980-2100
# WCRNT: includes ROMS real-near time from 2011-Oct2022
# WCRA: includes ROMS reanalysis from 1980-2010"
sstProj <- read_csv(paste0(datPath, "ROMS-ESMS_CalCOFI_monave_SST_198001-210012.csv"))

sstProj <- sstProj %>% mutate(seas = case_when(month %in% 1:6 ~ "spring",
                                               month %in% 7:12 ~ "summer")) %>%
  group_by(year, seas) %>%
  summarize(sstGFDL = mean(sst_roms_gfdl),
            sstIPSL = mean(sst_roms_ipsl),
            sstHAD = mean(sst_roms_had)) %>% 
  pivot_longer(cols = c(sstGFDL, sstIPSL, sstHAD),
               names_to = "ESM",
               names_prefix = "sst",
               values_to = "SST")

sstProj %>% ggplot(aes(x = year, y = SST, color = ESM)) +
  geom_line() +
  facet_wrap(~seas, ncol = 1)+
  theme_classic()

# Northward transport across 32N
# Downloaded from Mike Jacox (8/18/2023)
transportOffshore <- read_ncdf(paste0(datPath, "wcra31_monthly_northward_velocity_0-50m_124.0-120.0W_31.9-32.1N_1981-2010.nc"))

offshoreTrans <- data.frame(year = transportOffshore$year,
                            month = transportOffshore$month,
                            offshTrans = transportOffshore$v)

transport <- offshoreTrans %>% mutate(seas = case_when(month %in% 1:6 ~ "spring",
                                                       month %in% 6:12 ~ "summer"),
                                      season = case_when(month %in% 3:5 ~ "spring",
                                                         month %in% 6:8 ~ "summer",
                                                         month %in% 9:11 ~ "fall",
                                                         month %in% c(12,1,2) ~ "winter")) %>%
  group_by(year, seas) %>%
  summarize(avgOffshTran = mean(offshTrans)) %>% 
  # ggplot(aes(x = year, y = avgOffshTran, color = season)) +
  # geom_line()
  
  pivot_wider(names_from = seas, names_prefix = "avgOffTrans",
              values_from = avgOffshTran)

nearshoreGFDL <- read_ncdf(paste0(datPath, "monthly_northward_velocity_gfdl_0-50m_120.0-115.0W_31.9-32.1N_1980-2100.nc"))
nearshoreIPSL <- read_ncdf(paste0(datPath, "monthly_northward_velocity_ipsl_0-50m_120.0-115.0W_31.9-32.1N_1980-2100.nc"))
nearshoreHAD <- read_ncdf(paste0(datPath, "monthly_northward_velocity_had_0-50m_120.0-115.0W_31.9-32.1N_1980-2100.nc"))

nearshoreTrans <- rbind(data.frame(year = nearshoreGFDL$year,
                             month = nearshoreGFDL$month,
                             nearshTrans = nearshoreGFDL$v,
                             ESM = "GFDL"),
                        data.frame(year = nearshoreIPSL$year,
                                   month = nearshoreIPSL$month,
                                   nearshTrans = nearshoreIPSL$v,
                                   ESM = "IPSL"),
                        data.frame(year = nearshoreHAD$year,
                                   month = nearshoreHAD$month,
                                   nearshTrans = nearshoreHAD$v,
                                   ESM = "HAD"))


transport <- nearshoreTrans %>% mutate(seas = case_when(month %in% 1:6 ~ "spring",
                                                        month %in% 6:12 ~ "summer"),
                                       season = case_when(month %in% 3:5 ~ "spring",
                                                          month %in% 6:8 ~ "summer",
                                                          month %in% 9:11 ~ "fall",
                                                          month %in% c(12,1,2) ~ "winter")) %>%
  group_by(year, seas, ESM) %>%
  summarize(avgNearshTransp = mean(nearshTrans)) 

transport %>% 
  ggplot(aes(x = year, y = avgNearshTransp, color = ESM)) +
  geom_line() + 
  facet_grid(rows = vars(ESM), cols = vars(seas))+
  theme_classic()

  # pivot_wider(names_from = seas, names_prefix = "avgNearTrans",
  #             values_from = avgNearshTran) %>%
  # full_join(y = transport, by = "year")

# NEMURO plankton
nemuroZ <- read_csv(paste0(datPath, "ROMS-NEMURO spring plankton.csv"))

nemuroZ %>% filter(GCM != "HIST") %>%
              ggplot(aes(x = year, y = ZM, color = GCM)) +
              geom_line() +
              facet_wrap(~zone, ncol = 1)+
              theme_classic()

nemuroZ <- nemuroZ %>% filter(GCM != "HIST") %>%           
              pivot_wider(names_from = zone, values_from = c(PS, PL, ZS, ZM, ZL))


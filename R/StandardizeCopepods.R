# GAM standardization of PRPOOS copepod survey indices
# Created: 8/14/2023, Robert Wildermuth

library(tidyverse)
library(mgcv)

# data file path
datPath <- "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/"

# Copepods ----------------------------------------------------------------

#Size fraction data for copepods
copepodSize <- read_csv(file= paste0(datPath, "SpringPRPOOSdata_Copepods_sizeFractions.csv"),
                        na = "NaN")
copepodSize <- copepodSize %>% mutate(Time_PST = dmy_hms(Time_PST),
                                      year = year(Time_PST), month = month(Time_PST)) %>%
                 filter(!is.na(year),
                        Night0Day1 == 1) # larvae only forage during day

# need to get total biomass C per sample
# apply maximum size of 0.95 mm based on max width of larval forage ~400 micrometers
# in larvae <20 mm and length:width ratio of 2.306 (Arthur 1976, 1977)
copTotBio <- copepodSize %>% mutate(size = str_match(string = `Size Range (Major Axis)`, 
                                                     pattern = "-\\s*(.*?)\\s*mm")[,2],
                                    size = as.numeric(size)) %>%
                filter(size <= 0.95) %>%
                group_by(year, month, Cruise, Line, Station, Time_PST, Latitude, Longitude) %>%
                summarize(totBio = sum(Carbon_mg_sum)) %>%
                # get Julien day for standardization
                mutate(dayYr = yday(Time_PST))

# convert year to factor
copTotBio$year <- as.factor(copTotBio$year)

copTotBio$station <- as.factor(copTotBio$Station)

# sample sizes by year/line
with(copTotBio, table(year,Station, Line))

########## fit a Tweedie GAMM using mgcv
# Set this up like Hunsicker et al. 2022 Appendix
copFit <- gam(totBio ~ year + te(Latitude, Longitude, bs = "fs"),
              family = Tweedie(p = 1.25, link = "log"),
              data = copTotBio, method = "REML",
              control = gam.control(nthreads = 3))
gam.check(copFit)

# save(copFit,
#      file = "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/copepodIndex.RData")

# compare to day-of-year smoother
copDoYFit <- gam(totBio ~ year + s(dayYr) + te(Latitude, Longitude, bs = "fs"),
              family = Tweedie(p = 1.25, link = "log"),
              data = copTotBio, method = "REML",
              control = gam.control(nthreads = 3))
gam.check(copDoYFit)


# Nauplii -----------------------------------------------------------------

#Size fraction data for nauplii 
naupliiSize <- read_csv(file= paste0(datPath, "SpringPRPOOSdata_Nauplii_sizeFractions.csv"),
                        na = "NaN")
naupliiSize <- naupliiSize %>% mutate(Time_PST = dmy_hms(Time_PST),
                                      year = year(Time_PST), month = month(Time_PST)) %>%
                  filter(!is.na(year),
                         Night0Day1 == 1) # larvae only forage during day

# need to get total biomass C per sample
# apply maximum size of 0.95 mm based on max width of larval forage ~400 micrometers
# in larvae <20 mm and length:width ratio of 2.306 (Arthur 1976, 1977)
naupTotBio <- naupliiSize %>% mutate(size = str_match(string = `Size Range (Major Axis)`, 
                                                     pattern = "-\\s*(.*?)\\s*mm")[,2],
                                     size = as.numeric(size)) %>%
                filter(size <= 0.95) %>%
                group_by(year, month, Cruise, Line, Station, Time_PST, Latitude, Longitude) %>%
                summarize(totBio = sum(Carbon_mg_sum)) %>%
                # get Julien day for standardization
                mutate(dayYr = yday(Time_PST))

# convert year to factor
naupTotBio$year <- as.factor(naupTotBio$year)

########## fit a Tweedie GAMM using mgcv
# Set this up like Hunsicker et al. 2022 Appendix
naupFit <- gam(totBio ~ year + te(Latitude, Longitude, bs = "fs"),
              family = Tweedie(p = 1.25, link = "log"),
              data = naupTotBio, method = "REML",
              control = gam.control(nthreads = 3))
gam.check(naupFit)

# save(naupFit,
#      file = "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/naupliiIndex.RData")

# compare to day-of-year smoother
naupDoYFit <- gam(totBio ~ year + s(dayYr) + te(Latitude, Longitude, bs = "fs"),
                 family = Tweedie(p = 1.25, link = "log"),
                 data = naupTotBio, method = "REML",
                 control = gam.control(nthreads = 3))
gam.check(naupDoYFit)

# delta GLMM standardization of CalCOFI mesopelagics
# based on code from Charlie Hinchliffe
# Created: 6/5/2023, Charlie Hinchliffe, Robert Wildermuth

library(tidyverse)
library(rstanarm)
library(mgcv)

# CalCOFI larval data provided by Andrew Thompson 5/23/2023 for 66 regularly sampled core stations
mydata <- read.csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/CalCOFI_ichthioplankton_20230523.csv")

mesopelDat <- mydata %>% 
                select(unique.code, S_C, S_SC, S_L, S_S, longitude, latitude, year, season,
                       Protomyctophum.crockeri, Nannobrachium1, Triphoturus.mexicanus, 
                       Symbolophorus.californiensis, Bathylagoides.wesethi, 
                       Diogenichthys.atlanticus, Chauliodus.macouni, Gonostomatidae1, 
                       Vinciguerria1) %>%
                # mutate(std.egg.count = (raw.egg.count/percent.sorted)*standard.haul.factor)  %>% 
                mutate(cpue = rowSums(across(c(Protomyctophum.crockeri, Nannobrachium1, Triphoturus.mexicanus, 
                                  Symbolophorus.californiensis, Bathylagoides.wesethi, 
                                  Diogenichthys.atlanticus, Chauliodus.macouni, 
                                  Gonostomatidae1, Vinciguerria1)))) %>%
                filter(season == "spring") %>% 
                rename(station = S_S) 

summary(is.finite(mesopelDat$cpue))

#Generate distance blocks from station numbers
mesopelDat$dist_block <- NA_character_ 
mesopelDat$dist_block[mesopelDat$station < 55] <- "Block1"
mesopelDat$dist_block[mesopelDat$station >= 55 & mesopelDat$station < 70] <- "Block2"
mesopelDat$dist_block[mesopelDat$station >= 70 & mesopelDat$station < 90] <- "Block3"
mesopelDat$dist_block[mesopelDat$station >= 90] <- "Block4"

# see about latitude blocks
mesopelDat <- mesopelDat %>%
                mutate(latBlock = case_when(latitude < 32.5 ~ "South",
                                            latitude >= 32.5 ~ "North"))
with(mesopelDat, table(dist_block, latBlock))
with(mesopelDat, table(year, dist_block, latBlock))

# add column with response as binary
mesopelDat$bin <- as.numeric(mesopelDat$cpue>0)

# get vector of years for plots (later)
yr.vec <- unique(mesopelDat$year)

# convert year to factor
mesopelDat$year <- as.factor(mesopelDat$year)

mesopelDat$station <- as.factor(mesopelDat$station)

# sample sizes by year/line
with(mesopelDat, table(year,dist_block))


# try removing intercept
mesopelDat$int <- 1

########## fit a Tweedie GAMM using mgcv
# Set this up like Hunsicker et al. 2022 Appendix
gamFit <- gam(cpue ~ year + te(latitude, longitude, bs = "fs"),
              family = Tweedie(p = 1.25, link = "log"),
              data = mesopelDat, method = "REML",
              control = gam.control(nthreads = 3))
gam.check(gamFit)
mesopelLarv <- gamFit
# save(mesopelLarv,
#      file = "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/soMesopelLarvalIndex.RData")

# For checking with other methods
########## fit a binomial GLM using rstanarm
# using default priors for example only... better to think more about this
binfit <- stan_glm(bin ~ year + latitude + longitude - 1,
                   family = binomial,
                   data = mesopelDat,
                   iter = 8000,
                   cores = parallel::detectCores()-2,
                   thin = 3)

summary(binfit, digits=4)

prior_summary(binfit)
# diagnostics using shiny interface:
# some parameters have effective sample sizes that are <10% of total sample size
# uncomment next line to launch diagnostic plots, etc. in web browser
# launch_shinystan(binfit)

# compare stan point estimates to glm()
binfit.mle <- glm(bin ~ year + latitude + longitude -1,
                  family=binomial,
                  data=mesopelDat)
summary(binfit.mle)
# default priors aren't that informative in this example
plot(coef(binfit), coef(binfit.mle), xlab="posterior medians from STAN", ylab="MLEs from glm()"); abline(0,1)

# get dataframe of posterior draws
# NOTE THAT R USES 'TREATMENT' CONTRASTS BY DEFAULT,
# SO YEAR EFFECTS ARE OFFSETS FROM THE REFERENCE LEVEL
bin.draws <- as.data.frame(binfit)
summary(bin.draws); dim(bin.draws)

# get logit scale prob's for year effects
# these are for the baseline 'line' level, but it's a relative index, so not an issue
# bin.yrs.logit <- cbind.data.frame(year1951=bin.draws[,1],
#                                   bin.draws[,1] + bin.draws[,2:length(yr.vec)])
bin.yrs.logit <- bin.draws %>% select(grep("year", names(bin.draws), 
                                           fixed = TRUE, value = TRUE))

summary(bin.yrs.logit)
# convert to proportions
bin.yrs.prop <- exp(bin.yrs.logit) / (1 + exp(bin.yrs.logit))
summary(bin.yrs.prop)
# plot annual proportion positive
plot(rev(yr.vec), apply(bin.yrs.prop, 2, mean), type='o', xlab='year', ylab='proportion positive')

########## fit the model for the conditional mean (given positive)
# I'm using a gamma, but you can do whatever
gamfit <- stan_glm(cpue ~ year + latitude + longitude -1,
                   family = Gamma(link='log'), # gaussian(),
                   data = subset(mesopelDat, cpue>0),
                   iter = 8000,
                   cores = parallel::detectCores()-2,
                   thin = 3)

summary(gamfit, digits=4)

prior_summary(gamfit)
# uncomment next line to launch diagnostic plots, etc. in web browser
# launch_shinystan(gamfit)

# compare stan point estimates to glm()
gamfit.mle <- glm(cpue ~ year + latitude + longitude -1,
                  family = Gamma(link='log'),
                  data = subset(mesopelDat, cpue>0))

summary(gamfit.mle)

# default priors aren't that informative in this example
plot(coef(gamfit), coef(gamfit.mle), xlab="posterior medians from STAN", ylab="MLEs from glm()"); abline(0,1)

# get dataframe of posterior draws
# NOTE THAT R USES 'TREATMENT' CONTRASTS BY DEFAULT,
# SO YEAR EFFECTS ARE OFFSETS FROM THE REFERENCE LEVEL
gamma.draws <- as.data.frame(gamfit)
summary(gamma.draws); dim(gamma.draws)

# get log scale prob's for year effects
# these are for the baseline 'line' level, but it's a relative index
# the 'old' delta_glm function calculated "least squares means" (added the mean of the other covariates' levels)
# gamma.yrs.log <- cbind.data.frame(year1951=gamma.draws[,1],
#                                   gamma.draws[,1] + gamma.draws[,2:length(yr.vec)])
gamma.yrs.log <- gamma.draws %>% select(grep("year", names(gamma.draws), 
                                           fixed = TRUE, value = TRUE))

summary(gamma.yrs.log)

# convert to arithmetic scale
gamma.yrs <- exp(gamma.yrs.log)

summary(gamma.yrs)

# plot annual conditional means
plot(rev(yr.vec), apply(gamma.yrs, 2, mean), type='o', xlab='year', ylab='conditional mean CPUE')

########## use the two GLMs to compute the index
# multiply the proportion positives by the conditional means
index.draws <- bin.yrs.prop * gamma.yrs

# get mean index and 95% interval
index <- apply(index.draws, 2, median)
logSE <- apply(index.draws, 2, function(x) { sd(log(x)) }) # this is a common input to our assessment models
int95 <- apply(index.draws, 2, quantile, probs=c(0.025,0.975))
index.df <- cbind.data.frame(index, logSE, t(int95))

index.df 


# formate for SS adding season and survey numbers
years <- as.data.frame(yr.vec) %>% arrange(yr.vec)

ss_data <- cbind(years, index.df) 

ss_data <- ss_data %>% 
  mutate(seas = 1, survey = 6) %>% 
  select(yr.vec, seas, survey, index, logSE) %>% 
  rename(year = yr.vec, obs = index, se_log = logSE) %>% 
  rename(index = survey)

#write.csv(ss_data, "data/CalCOFI_index_year_distblock_gamma.csv", row.names = FALSE)

# plot index and 95% interval
ss_data2 <- cbind(years, index.df) %>% 
  ungroup()

ss_data2 <- ss_data2 %>% 
  mutate(seas = 1, survey = 6, yr.vex = as.numeric(yr.vec))  %>% 
  rename(year = yr.vec, obs = index, se_log = logSE, low = "2.5%", high = "97.5%") %>% 
  mutate(year = as.numeric(as.character(year)) ,low = as.numeric(low), high = as.numeric(high)) %>% 
  rename(index = survey) 

mesopelPlot <- ggplot(ss_data2, aes(year, obs)) +
                geom_line() +
                geom_point() +
                geom_line(data = ss_data2, aes(x = year, y = low, col = "red")) +
                geom_line(data = ss_data2, aes(x = year, y = high, col = "red")) +
                theme_bw() + 
                theme(legend.position = "none")

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

# Compare to abundance index
indxAnom <- data.frame(year = years,
                       abund = apply(gamma.yrs, 2, median)) %>% 
              mutate(indx = scale(log(abund), scale = FALSE)) %>%
              rename(year = yr.vec)
ggplot(indxAnom, aes(x = year, y = indx)) +
  geom_line() +
  geom_point() + 
  geom_line(data = datDFA %>% filter(!is.na(sprCalCOFISouthernMesopels)), 
            aes(x = year, y = scale(sprCalCOFISouthernMesopels, scale = FALSE)),
            col = "darkblue")

# Pull GAM index
indxGAM <- data.frame(year = rev(yr.vec),
                      coefGAM = as.vector(coef(gamFit))[1:length(yr.vec)])
indxGAM$index <- indxGAM$coefGAM
indxGAM <- indxGAM %>% mutate(index = index + c(0, coefGAM[2:length(yr.vec)]))

ggplot(ss_data2, aes(year, scale(log(obs), scale = FALSE))) +
  geom_line() +
  geom_point() +
  # geom_line(data = ss_data2, aes(x = year, y = low, col = "red")) +
  # geom_line(data = ss_data2, aes(x = year, y = high, col = "red")) +
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_line(data = datDFA %>% filter(!is.na(sprCalCOFISouthernMesopels)), 
            aes(x = year, y = scale(sprCalCOFISouthernMesopels, scale = FALSE)),
            col = "darkblue") +
  geom_line(data = indxGAM, 
            aes(x = year, y = scale(index, center = FALSE)),
            col = "orangered")

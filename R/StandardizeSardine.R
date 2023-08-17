# delta GLMM standardization of CalCOFI mesopelagics
# based on code from Charlie Hinchliffe
# Created: 6/5/2023, Charlie Hinchliffe, Robert Wildermuth

library(tidyverse)
library(rstanarm)
library(mgcv)

# CalCOFI larval data provided by Andrew Thompson 5/23/2023 for 66 regularly sampled core stations
# data is scaled for proportion of sample sorted and unit is No./10m^2 (cylindrical)
mydata <- read.csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/CalCOFI_ichthioplankton_20230523.csv")

sardDat <- mydata %>% 
                select(unique.code, S_C, S_SC, S_L, S_S, longitude, latitude, year, season,
                       Sardinops.sagax) %>%
                rename(cpue = Sardinops.sagax,
                       station = S_S) %>%
                filter(season == "spring") 

summary(is.finite(sardDat$cpue))

#Generate distance blocks from station numbers
sardDat$dist_block <- NA_character_ 
sardDat$dist_block[sardDat$station < 55] <- "Block1"
sardDat$dist_block[sardDat$station >= 55 & sardDat$station < 70] <- "Block2"
sardDat$dist_block[sardDat$station >= 70 & sardDat$station < 90] <- "Block3"
sardDat$dist_block[sardDat$station >= 90] <- "Block4"

pac.coast <- borders("world", colour="gray50", fill="gray50", 
                     xlim = c(-126, -116), ylim = c(29, 36))
sardDat %>% ggplot(aes(x = longitude, y = latitude, col = dist_block)) +
  geom_point() +
  labs(title = "Sample locations w/in Blocks") +
  pac.coast +
  coord_sf(xlim = c(-126, -116), ylim = c(29, 36)) 

# see about latitude blocks
sardDat <- sardDat %>%
                mutate(latBlock = case_when(latitude < 32.5 ~ "South",
                                            latitude >= 32.5 ~ "North"))
with(sardDat, table(dist_block, latBlock))
with(sardDat, table(year, dist_block, latBlock))

# add column with response as binary
sardDat$bin <- as.numeric(sardDat$cpue>0)

# plot confirmed observations
sardDat %>% 
  filter(bin == 1) %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_point(aes(size = cpue, col = dist_block, alpha = 0.1)) +
  pac.coast +
  coord_sf(xlim = c(-126, -116), ylim = c(29, 36)) +
  labs(title = "Sardine observations w/in Blocks") +
  geom_segment(data = data.frame(x = c(-116, -117, -118, -120, -122, -124),
                                 y = c(31, 33, 33, 33, 33, 33),
                                 xend = c(-124, -124, -118, -120, -122, -124),
                                 yend = c(31, 33, 31, 31, 31, 31)),
               mapping = aes(x = x, y = y, xend = xend, yend = yend))

# get vector of years for plots (later)
yr.vec <- unique(sardDat$year)

# convert year to factor
sardDat$year <- as.factor(sardDat$year)

sardDat$station <- as.factor(sardDat$station)

# sample sizes by year/line
with(sardDat, table(year,dist_block))


# try removing intercept
sardDat$int <- 1

########## fit a Tweedie GAMM using mgcv
# Set this up like Hunsicker et al. 2022 Appendix
gamFit <- gam(cpue ~ year + te(latitude, longitude, bs = "fs"),
              family = Tweedie(p = 1.25, link = "log"),
              data = sardDat, method = "REML",
              control = gam.control(nthreads = 3))
gam.check(gamFit)

sardLarv <- gamFit
# save(sardLarv,
#      file = "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/sardineLarvalIndex.RData")

########## fit a binomial GLM using rstanarm
# using default priors for example only... better to think more about this
# binfit <- stan_glm(bin ~ year + dist_block - 1,
#                    family = binomial,
#                    data = sardDat,
#                    iter=4000,
#                    cores = parallel::detectCores()-2,
#                    thin = 3)

# try with location coefficients instead
binfit <- stan_glm(bin ~ year + latitude + longitude - 1,
                   family = binomial,
                   data = sardDat,
                   iter=4000,
                   cores = parallel::detectCores()-2,
                   thin = 1)

summary(binfit, digits=4)

prior_summary(binfit)
# diagnostics using shiny interface:
# some parameters have effective sample sizes that are <10% of total sample size
# uncomment next line to launch diagnostic plots, etc. in web browser
# launch_shinystan(binfit)

# compare stan point estimates to glm()
# binfit.mle <- glm(bin ~ year + dist_block -1,
binfit.mle <- glm(bin ~ year + latitude + longitude -1,
                  family=binomial,
                  data=sardDat)
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
# gamfit <- stan_glm(cpue ~ year + dist_block - 1,
#                    family = Gamma(link='log'), # gaussian(),
#                    data = subset(sardDat, cpue>0),
#                    iter = 4000,
#                    cores = parallel::detectCores()-2,
#                    thin = 3)

# try with location coefficients instead
gamfit <- stan_glm(cpue ~ year + latitude + longitude - 1,
                   family = Gamma(link='log'), # gaussian(),
                   data = subset(sardDat, cpue>0),
                   iter = 4000,
                   cores = parallel::detectCores()-2,
                   thin = 1)

summary(gamfit, digits=4)

prior_summary(gamfit)
# uncomment next line to launch diagnostic plots, etc. in web browser
# launch_shinystan(gamfit)

# compare stan point estimates to glm()
# gamfit.mle <- glm(cpue ~ year + dist_block -1,
gamfit.mle <- glm(cpue ~ year + latitude + longitude -1,
                  family = Gamma(link='log'),
                  data = subset(sardDat, cpue>0))

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

# some years with no positive cpue
posYrs <- as.numeric(sub(pattern = "year", replacement = "", x = names(gamma.yrs), fixed = TRUE))

# plot annual conditional means
plot(posYrs, apply(gamma.yrs, 2, mean), type='o', xlab='year', ylab='conditional mean CPUE')

########## use the two GLMs to compute the index
# multiply the proportion positives by the conditional means
index.draws <- bin.yrs.prop[, names(gamma.yrs)] * gamma.yrs

# get mean index and 95% interval
index <- apply(index.draws, 2, median)
logSE <- apply(index.draws, 2, function(x) { sd(log(x)) }) # this is a common input to our assessment models
int95 <- apply(index.draws, 2, quantile, probs=c(0.025,0.975))
index.df <- cbind.data.frame(index, logSE, t(int95))

index.df 


# formate for SS adding season and survey numbers
# years <- as.data.frame(yr.vec) %>% arrange(yr.vec)

ss_data <- cbind(years = posYrs, index.df) 

ss_data <- ss_data %>% 
  mutate(seas = 1, survey = 6) %>% 
  select(yr.vec, seas, survey, index, logSE) %>% 
  rename(year = yr.vec, obs = index, se_log = logSE) %>% 
  rename(index = survey)

#write.csv(ss_data, "data/CalCOFI_index_year_distblock_gamma.csv", row.names = FALSE)

# plot index and 95% interval
ss_data2 <- cbind(years = posYrs, index.df) %>% 
  ungroup()

ss_data2 <- ss_data2 %>% 
  #mutate(seas = 1, survey = 6, yr.vex = as.numeric(yr.vec))  %>% 
  rename(year = years, obs = index, se_log = logSE, low = "2.5%", high = "97.5%") %>% 
  mutate(year = as.numeric(as.character(year)) ,low = as.numeric(low), high = as.numeric(high)) #%>% 
  #rename(index = survey) 

sardPlot <- ggplot(ss_data2, aes(year, obs)) +
                geom_line() +
                geom_point() +
                geom_line(data = ss_data2, aes(x = year, y = low, col = "red")) +
                geom_line(data = ss_data2, aes(x = year, y = high, col = "red")) +
                theme_bw() + 
                theme(legend.position = "none")

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")

# Compare to abundance index
indxAnom <- data.frame(year = posYrs,
                       abund = apply(gamma.yrs, 2, median)) %>% 
              mutate(indx = scale(log(abund), scale = FALSE)) 
ggplot(indxAnom, aes(x = year, y = indx)) +
  geom_line() +
  geom_point() + 
  geom_line(data = datDFA %>% filter(!is.na(sprCalCOFILarvalSardine)), 
            aes(x = year, y = scale(sprCalCOFILarvalSardine, scale = FALSE)),
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
  geom_line(data = datDFA %>% filter(!is.na(sprCalCOFILarvalSardine)), 
            aes(x = year, y = scale(sprCalCOFILarvalSardine, scale = FALSE)),
            col = "darkblue") +
  geom_line(data = indxGAM, 
            aes(x = year, y = index),
            col = "orangered")

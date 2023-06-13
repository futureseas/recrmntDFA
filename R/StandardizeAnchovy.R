# delta GLMM standardization of CalCOFI mesopelagics
# based on code from Charlie Hinchliffe
# Created: 6/5/2023, Charlie Hinchliffe, Robert Wildermuth

library(tidyverse)
library(rstanarm)

# CalCOFI larval data provided by Andrew Thompson 5/23/2023 for 66 regularly sampled core stations
mydata <- read.csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/CalCOFI_ichthioplankton_20230523.csv")

anchDat <- mydata %>% 
                select(unique.code, S_C, S_SC, S_L, S_S, longitude, latitude, year, season,
                       Engraulis.mordax) %>%
                rename(cpue = Engraulis.mordax,
                       station = S_S) %>%
                filter(season == "spring") 

summary(is.finite(anchDat$cpue))

#Generate distance blocks from station numbers
anchDat$dist_block <- NA_character_ 
anchDat$dist_block[anchDat$station < 55] <- "Block1"
anchDat$dist_block[anchDat$station >= 55 & anchDat$station < 70] <- "Block2"
anchDat$dist_block[anchDat$station >= 70 & anchDat$station < 90] <- "Block3"
anchDat$dist_block[anchDat$station >= 90] <- "Block4"

# see about latitude blocks
anchDat <- anchDat %>%
                mutate(latBlock = case_when(latitude < 32.5 ~ "South",
                                            latitude >= 32.5 ~ "North"))
with(anchDat, table(dist_block, latBlock))
with(anchDat, table(year, dist_block, latBlock))

# add column with response as binary
anchDat$bin <- as.numeric(anchDat$cpue>0)

# plot confirmed observations
pac.coast <- borders("world", colour="gray50", fill="gray50", 
                     xlim = c(-126, -116), ylim = c(29, 36))
anchDat %>% 
  filter(bin == 1) %>%
  ggplot(aes(x = longitude, y = latitude)) +
  geom_point(aes(size = cpue, col = dist_block, alpha = 0.1)) +
  pac.coast +
  coord_sf(xlim = c(-126, -116), ylim = c(29, 36)) +
  labs(title = "Anchovy observations w/in Blocks") +
  geom_segment(data = data.frame(x = -116, 
                                 y = 31, 
                                 xend = -124,
                                 yend = 31),
               mapping = aes(x = x, y = y, xend = xend, yend = yend))

# get vector of years for plots (later)
yr.vec <- unique(anchDat$year)

# convert year to factor
anchDat$year <- as.factor(anchDat$year)

anchDat$station <- as.factor(anchDat$station)

# sample sizes by year/line
with(anchDat, table(year,dist_block))


# try removing intercept
anchDat$int <- 1

########## fit a binomial GLM using rstanarm
# using default priors for example only... better to think more about this
binfit <- stan_glm(bin ~ year + dist_block - 1,
                   family = binomial,
                   data = anchDat,
                   iter=4000,
                   cores = parallel::detectCores()-2,
                   thin = 3)

summary(binfit, digits=4)

prior_summary(binfit)
# diagnostics using shiny interface:
# some parameters have effective sample sizes that are <10% of total sample size
# uncomment next line to launch diagnostic plots, etc. in web browser
# launch_shinystan(binfit)

# compare stan point estimates to glm()
binfit.mle <- glm(bin ~ year + dist_block -1,
                  family=binomial,
                  data=anchDat)
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
gamfit <- stan_glm(cpue ~ year + dist_block -1,
                   family = Gamma(link='log'), # gaussian(),
                   data = subset(anchDat, cpue>0),
                   iter = 4000,
                   cores = parallel::detectCores()-2,
                   thin = 3)

summary(gamfit, digits=4)

prior_summary(gamfit)
# uncomment next line to launch diagnostic plots, etc. in web browser
# launch_shinystan(gamfit)

# compare stan point estimates to glm()
gamfit.mle <- glm(cpue ~ year + dist_block -1,
                  family = Gamma(link='log'),
                  data = subset(anchDat, cpue>0))

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

anchPlot <- ggplot(ss_data2, aes(year, obs)) +
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
  geom_line(data = datDFA %>% filter(!is.na(sprCalCOFILarvalAnchovy)), 
            aes(x = year, y = scale(sprCalCOFILarvalAnchovy, scale = FALSE)),
            col = "darkblue")


ggplot(ss_data2, aes(year, scale(log(obs), scale = FALSE))) +
  geom_line() +
  geom_point() +
  # geom_line(data = ss_data2, aes(x = year, y = low, col = "red")) +
  # geom_line(data = ss_data2, aes(x = year, y = high, col = "red")) +
  theme_bw() + 
  theme(legend.position = "none") + 
  geom_line(data = datDFA %>% filter(!is.na(sprCalCOFILarvalAnchovy)), 
            aes(x = year, y = scale(sprCalCOFILarvalAnchovy, scale = FALSE)),
            col = "darkblue")

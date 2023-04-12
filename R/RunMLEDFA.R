# MLE dynamic factor analysis of sardine recruitment using MARSS package 
# Created: 4/4/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")


# Sardine model -----------------------------------------------------------

# subset for sardine DFA from 1990 to 2019
sardDat <- datDFA %>% filter(year %in% 1990:2019) %>%
              # remove anchovy-related timeseries
              select(-c(sprCalCOFILarvalAnchovy,
                        age1SprAnchmeanWAA,
                        anchSpawnHab,
                        anchBioSmrySeas1, 
                        anchBioSmrySeas2,
                        anchRec))
# remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
sardDat <- sardDat %>% select(-c(sardBioSmrySeas1, All_Copepods))

datNames <- names(sardDat)[-1]

# transpose for MARSS formatting
sardDat <- sardDat %>% select(-year) %>% t()

# do simple DFA fit

simpleDFA <- MARSS(y = sardDat, 
                   form = "dfa",
                   control = list(maxit = 1000,
                                  allow.degen = TRUE),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(#R = "diagonal and equal", # observation errors are the same
                                R = "diagonal and unequal", # observation errors independent
                                #R = "equalvarcov", # observation errors equal and covars equal
                                #R = "unconstrained", # all observation errors independent
                                m = 1) # one latent process
)
# Second trend is mostly explained by missing data in the zooplankton dataset. Stick to one trend

# plot(c(1990:2019), y= simpleDFA$states, type = "l")
# abline(h = 0)
# 
# # plot loadings
# # new df with coordinates
# loadingsDF <- data.frame(index = datNames,
#                          vals = simpleDFA$par$Z, 
#                          dummy0 = 0)
# 
# loadingsDF %>% ggplot(aes(y = index)) +
#   geom_segment(aes(x = dummy0,
#                    yend = index,
#                    xend = vals)) +
#   labs(x = "Loadings", y = "Index") +
#   geom_vline(xintercept = 0, color = "grey") +
#   theme_classic()

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(simpleDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 
ggplot(data = subset(d, name=="model")) + 
  geom_point(aes(t, value)) + 
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(t, .fitted), col="blue") + 
  facet_wrap(~.rownames) + xlab("Time Step") + ylab("Count")

# Look at factor loadings
# get the inverse of the rotation matrix 
Z.est <- coef(simpleDFA, type = "matrix")$Z 
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(simpleDFA, type = "matrix")$Z)$rotmat
} 
  
# rotate factor loadings 
Z.rot <- Z.est %*% H.inv 
# rotate trends 
trends.rot <- solve(H.inv) %*% simpleDFA$states

# Add CIs to marssMLE object 
simpleDFA <- MARSSparamCIs(simpleDFA) 
# Use coef() to get the upper and lower CIs 
Z.low <- coef(simpleDFA, type = "Z", what = "par.lowCI") 
Z.up <- coef(simpleDFA, type = "Z", what = "par.upCI") 
Z.rot.up <- Z.up %*% H.inv 
Z.rot.low <- Z.low %*% H.inv 
df <- data.frame(est = as.vector(Z.rot), 
                 conf.up = as.vector(Z.rot.up), 
                 conf.low = as.vector(Z.rot.low) )

# plot the factor loadings 
spp <- rownames(sardDat) 
N.ts <- nrow(sardDat)
minZ <- 0.05 
m <- dim(trends.rot)[1] 
ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot))) 
par(mfrow = c(ceiling(m / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1)) 
for (i in 1:m) { 
  plot(c(1:N.ts)[abs(Z.rot[, i]) > minZ], as.vector(Z.rot[abs(Z.rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1) ) 
  for (j in 1:N.ts) { 
    if (Z.rot[j, i] > minZ) { 
      text(j,-0.05, spp[j], srt = 90, adj = 1, cex = 0.9) 
      } 
    if (Z.rot[j, i] <-minZ) { 
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9) 
      } 
    abline(h = 0, lwd = 1, col = "gray") 
    } # end j loop 
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5) 
  } # end i loop

# get ts of trends
ts.trends <- t(trends.rot)
# par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
# loop over each trend
for (i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[, i],
       ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
       type = "n", lwd = 2, bty = "L",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  # draw zero-line
  abline(h = 0, col = "gray")
  # plot trend line
  par(new = TRUE)
  plot(ts.trends[, i],
       ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
       type = "l", lwd = 2, bty = "L",
       xlab = "", ylab = "", xaxt = "n"
  )
  # add panel labels
  mtext(paste("Trend", i, sep = " "), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(sardDat)[2]) + 1, 1980 + 0:dim(sardDat)[2])
} # end i loop (trends)


# Anchovy model -----------------------------------------------------------


# subset for anchovy DFA from 1990 to 2019
anchDat <- datDFA %>% filter(year %in% 1990:2019) %>%
  # remove sardine-related timeseries
  select(-c(sprCalCOFILarvalSardine,
            age1SprSardmeanWAA,
            sardSpawnHab,
            sardBioSmrySeas1, 
            sardBioSmrySeas2,
            sardRec))
# remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
anchDat <- anchDat %>% select(-c(All_Copepods,
                                 anchBioSmrySeas1,
                                 anchBioSmrySeas2,
                                 anchRec,
                                 age1SprAnchmeanWAA)) # not enough data in time window

datNames <- names(anchDat)[-1]

# transpose for MARSS formatting
anchDat <- anchDat %>% select(-year) %>% t()

# do simple DFA fit

anchDFA <- MARSS(y = anchDat, 
                   form = "dfa",
                   control = list(maxit = 1000,
                                  allow.degen = TRUE),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(#R = "diagonal and equal", # observation errors are the same
                     R = "diagonal and unequal", # observation errors independent
                     #R = "equalvarcov", # observation errors equal and covars equal
                     #R = "unconstrained", # all observation errors independent
                     m = 1) # one latent process
)
# Second trend is mostly explained by missing data in the zooplankton dataset. Stick to one trend

# plot(c(1990:2019), y= anchDFA$states, type = "l")
# abline(h = 0)
# 
# # plot loadings
# # new df with coordinates
# loadingsDF <- data.frame(index = datNames,
#                          vals = anchDFA$par$Z, 
#                          dummy0 = 0)
# 
# loadingsDF %>% ggplot(aes(y = index)) +
#   geom_segment(aes(x = dummy0,
#                    yend = index,
#                    xend = vals)) +
#   labs(x = "Loadings", y = "Index") +
#   geom_vline(xintercept = 0, color = "grey") +
#   theme_classic()

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(anchDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 
ggplot(data = subset(d, name=="model")) + 
  geom_point(aes(t, value)) + 
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(t, .fitted), col="blue") + 
  facet_wrap(~.rownames) + xlab("Time Step") + ylab("Count")

# Look at factor loadings
# get the inverse of the rotation matrix 
Z.est <- coef(anchDFA, type = "matrix")$Z 
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(anchDFA, type = "matrix")$Z)$rotmat
} 

# rotate factor loadings 
Z.rot <- Z.est %*% H.inv 
# rotate trends 
trends.rot <- solve(H.inv) %*% anchDFA$states

# Add CIs to marssMLE object 
anchDFA <- MARSSparamCIs(anchDFA) 
# Use coef() to get the upper and lower CIs 
Z.low <- coef(anchDFA, type = "Z", what = "par.lowCI") 
Z.up <- coef(anchDFA, type = "Z", what = "par.upCI") 
Z.rot.up <- Z.up %*% H.inv 
Z.rot.low <- Z.low %*% H.inv 
df <- data.frame(est = as.vector(Z.rot), 
                 conf.up = as.vector(Z.rot.up), 
                 conf.low = as.vector(Z.rot.low) )

# plot the factor loadings 
spp <- rownames(anchDat) 
N.ts <- nrow(anchDat)
minZ <- 0.05 
m <- dim(trends.rot)[1] 
ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot))) 
par(mfrow = c(ceiling(m / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1)) 
for (i in 1:m) { 
  plot(c(1:N.ts)[abs(Z.rot[, i]) > minZ], as.vector(Z.rot[abs(Z.rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1) ) 
  for (j in 1:N.ts) { 
    if (Z.rot[j, i] > minZ) { 
      text(j,-0.05, spp[j], srt = 90, adj = 1, cex = 0.9) 
    } 
    if (Z.rot[j, i] <-minZ) { 
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9) 
    } 
    abline(h = 0, lwd = 1, col = "gray") 
  } # end j loop 
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5) 
} # end i loop

# get ts of trends
ts.trends <- t(trends.rot)
# par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
# loop over each trend
for (i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[, i],
       ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
       type = "n", lwd = 2, bty = "L",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  # draw zero-line
  abline(h = 0, col = "gray")
  # plot trend line
  par(new = TRUE)
  plot(ts.trends[, i],
       ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
       type = "l", lwd = 2, bty = "L",
       xlab = "", ylab = "", xaxt = "n"
  )
  # add panel labels
  mtext(paste("Trend", i, sep = " "), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(anchDat)[2]) + 1, 1980 + 0:dim(anchDat)[2])
} # end i loop (trends)

# All together -----------------------------------------------------------


# subset for anchovy DFA from 1990 to 2019
allDat <- datDFA %>% filter(year %in% 1990:2019) %>%
            # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
            select(-c(All_Copepods,
                      anchBioSmrySeas1, 
                      sardBioSmrySeas1,
                      anchBioSmrySeas2, # leave out biomass since not fit well
                      sardBioSmrySeas2,
                      NCOPspring,
                      SCOPspring, # summer copepod index had higher loadings
                      PDOsummer, # lower loading than spring, may want to try a lag
                      BEUTI_33N, # oceanography at 39N had highest loadings
                      OC_LUSI_33N,
                      OC_LUSI_36N,
                      OC_STI_33N,
                      OC_STI_36N,
                      age1SprAnchmeanWAA)) # not enough data in time windowselect(-c(#sprCalCOFILarvalSardine,

datNames <- names(allDat)[-1]

# transpose for MARSS formatting
allDat <- allDat %>% select(-year) %>% t()

# do simple DFA fit

overallDFA <- MARSS(y = allDat, 
                 form = "dfa",
                 control = list(maxit = 2000,
                                allow.degen = TRUE),
                 inits = list(x0 = matrix(1, 1, 1)),
                 z.score = TRUE,
                 model = list(#R = "diagonal and equal", # observation errors are the same
                               R = "diagonal and unequal", # observation errors independent
                               #R = "equalvarcov", # observation errors equal and covars equal
                               #R = "unconstrained", # all observation errors independent
                               m = 1) # one latent process
)
# Second trend is mostly explained by missing data in the zooplankton dataset. Stick to one trend

# plot(c(1990:2019), y= overallDFA$states, type = "l")
# abline(h = 0)
# 
# # plot loadings
# # new df with coordinates
# loadingsDF <- data.frame(index = datNames,
#                          vals = overallDFA$par$Z,
#                          dummy0 = 0)
# 
# loadingsDF %>% ggplot(aes(y = index)) +
#   geom_segment(aes(x = dummy0,
#                    yend = index,
#                    xend = vals)) +
#   labs(x = "Loadings", y = "Index") +
#   geom_vline(xintercept = 0, color = "grey") +
#   theme_classic()

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(overallDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 
ggplot(data = subset(d, name=="model")) + 
  geom_point(aes(t, value)) + 
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(t, .fitted), col="blue") + 
  facet_wrap(~.rownames) + xlab("Time Step") + ylab("Count")

# Look at factor loadings
# get the inverse of the rotation matrix 
Z.est <- coef(overallDFA, type = "matrix")$Z 
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(overallDFA, type = "matrix")$Z)$rotmat
} 

# rotate factor loadings 
Z.rot <- Z.est %*% H.inv 
# rotate trends 
trends.rot <- solve(H.inv) %*% overallDFA$states

# Add CIs to marssMLE object 
overallDFA <- MARSSparamCIs(overallDFA) 
# Use coef() to get the upper and lower CIs 
Z.low <- coef(overallDFA, type = "Z", what = "par.lowCI") 
Z.up <- coef(overallDFA, type = "Z", what = "par.upCI") 
Z.rot.up <- Z.up %*% H.inv 
Z.rot.low <- Z.low %*% H.inv 
df <- data.frame(est = as.vector(Z.rot), 
                 conf.up = as.vector(Z.rot.up), 
                 conf.low = as.vector(Z.rot.low) )

# plot the factor loadings 
spp <- rownames(allDat) 
N.ts <- nrow(allDat)
minZ <- 0.05 
m <- dim(trends.rot)[1] 
ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot))) 
par(mfrow = c(ceiling(m / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1)) 
for (i in 1:m) { 
  plot(c(1:N.ts)[abs(Z.rot[, i]) > minZ], as.vector(Z.rot[abs(Z.rot[, i]) > minZ, i]), 
       type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1) ) 
  for (j in 1:N.ts) { 
    if (Z.rot[j, i] > minZ) { 
      text(j,-0.05, spp[j], srt = 90, adj = 1, cex = 0.9) 
    } 
    if (Z.rot[j, i] <-minZ) { 
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9) 
    } 
    abline(h = 0, lwd = 1, col = "gray") 
  } # end j loop 
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5) 
} # end i loop

# get ts of trends
ts.trends <- t(trends.rot)
# par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
# loop over each trend
for (i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[, i],
       ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
       type = "n", lwd = 2, bty = "L",
       xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  # draw zero-line
  abline(h = 0, col = "gray")
  # plot trend line
  par(new = TRUE)
  plot(ts.trends[, i],
       ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
       type = "l", lwd = 2, bty = "L",
       xlab = "", ylab = "", xaxt = "n"
  )
  # add panel labels
  mtext(paste("Trend", i, sep = " "), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(allDat)[2]) + 1, 1980 + 0:dim(allDat)[2])
} # end i loop (trends)


# Plots for ECCWO poster --------------------------------------------------

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(overallDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 

# combine trend states with variables for plotting with model fits
ts.trends <- data.frame(trendStates = ts.trends,
                        t = 1990:2019)

d %>% filter(name=="model",
             .rownames %in% c("NPGO", "anchSpawnHab",
                              "SCOPsummer", "sprCalCOFISouthernMesopels",
                              "anchRec", "sprCalCOFILarvalAnchovy")) %>%
  mutate(t = t+1989) %>%
  full_join(y = ts.trends, by = "t") %>%
  ggplot(aes(x = t)) + 
  geom_point(aes(y = value)) + 
  geom_ribbon(aes(ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(y = .fitted), col="blue") + 
  geom_line(aes(y = trendStates), col = "darkgreen") +
  facet_wrap(~.rownames) + xlab("Year") + ylab("Variable Anomaly") +
  theme_minimal() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

# plot loadings
# new df with coordinates
loadingsDF <- data.frame(index = rownames(allDat),
                         vals = overallDFA$par$Z,
                         dummy0 = 0)

# order loadings by magnitude and arrange for plotting
varArrang <- loadingsDF %>% arrange(vals) %>% pull(index)
# leave out response variables
varArrang <- varArrang[-which(varArrang %in% c("sprCalCOFILarvalAnchovy",
                                               "sardRec",
                                               "anchRec",
                                               "sprCalCOFILarvalSardine"))]
loadingsDF <- loadingsDF %>% mutate(index = factor(index, 
                                                   level = c("anchRec",
                                                      "sprCalCOFILarvalAnchovy",
                                                      "sardRec",
                                                      "sprCalCOFILarvalSardine",
                                                      varArrang)))

loadingsDF %>% ggplot(aes(y = index)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = vals,
                   linewidth = 4)) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 4.5, color = "black") +
  theme_classic()
# MLE dynamic factor analysis of sardine recruitment using MARSS package 
# Created: 4/4/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)

# read prepped dataset
datDFA <- read_csv("C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/recrmntDFA/recrDFAdat.csv")


# Sardine model -----------------------------------------------------------

# subset for sardine DFA from 1990 to 2019
sardDat <- datDFA %>% filter(year %in% 1950:2019) %>%
              # remove anchovy-related timeseries
              select(-c(anchLarv,
                        age1SprAnchmeanWAA,
                        anchSpawnHab,
                        daysAbove40pct,
                        anchNurseHab,
                        anchBioSmrySeas1, 
                        anchBioSmrySeas2,
                        anchRec))
# remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
sardDat <- sardDat %>% select(-c(sardBioSmrySeas1))

datNames <- names(sardDat)[-1]

# transpose for MARSS formatting
sardDat <- sardDat %>% select(-year) %>% t()

# do simple DFA fit

simpleDFA <- MARSS(y = sardDat, 
                   form = "dfa",
                   control = list(maxit = 1000,
                                  allow.degen = TRUE,
                                  conv.test.slope.tol = 0.1),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(#R = "diagonal and equal", # observation errors are the same
                                R = "diagonal and unequal", # observation errors independent
                                #R = "equalvarcov", # observation errors equal and covars equal
                                #R = "unconstrained", # all observation errors independent
                                m = 2) # one latent process
)
# Second trend is mostly explained by missing data in the zooplankton dataset. Stick to one trend

# plot(c(1990:2019), y= simpleDFA$states, type = "l")
# abline(h = 0)

# plot loadings
# new df with coordinates
loadingsDF <- data.frame(index = datNames,
                         vals = simpleDFA$par$Z,
                         dummy0 = 0)

loadingsDF %>% ggplot(aes(y = index)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = vals)) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  theme_classic()

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
  axis(1, 12 * (0:dim(sardDat)[2]) + 1, 1950 + 0:dim(sardDat)[2])
} # end i loop (trends)


# Anchovy model -----------------------------------------------------------


# subset for anchovy DFA from 1990 to 2019
anchDat <- datDFA %>% filter(year %in% 1950:2019) %>%
  # remove sardine-related timeseries
  select(-c(sardLarv,
            age1SprSardmeanWAA,
            meanSSBwt,
            sardSpawnHab,
            daysAbove5pct,
            sardNurseHab,
            sardBioSmrySeas1, 
            sardBioSmrySeas2,
            sardRec))
# remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
anchDat <- anchDat %>% select(-c(anchBioSmrySeas1,
                                 # anchBioSmrySeas2,
                                 # anchRec,
                                 age1SprAnchmeanWAA)) # not enough data in time window

datNames <- names(anchDat)[-1]

# transpose for MARSS formatting
anchDat <- anchDat %>% select(-year) %>% t()

# do simple DFA fit

anchDFA <- MARSS(y = anchDat, 
                   form = "dfa",
                   control = list(maxit = 1000,
                                  allow.degen = TRUE,
                                  conv.test.slope.tol = 0.1),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(#R = "diagonal and equal", # observation errors are the same
                     R = "diagonal and unequal", # observation errors independent
                     #R = "equalvarcov", # observation errors equal and covars equal
                     #R = "unconstrained", # all observation errors independent
                     m = 1) # one latent process
)

# plot(c(1990:2019), y= anchDFA$states, type = "l")
# abline(h = 0)

# plot loadings
# new df with coordinates
loadingsDF <- data.frame(index = datNames,
                         vals = anchDFA$par$Z,
                         dummy0 = 0)

loadingsDF %>% ggplot(aes(y = index)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = vals)) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  theme_classic()

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


# subset from 1980 to 2019
allDat <- datDFA %>% filter(year %in% 1980:2019) %>%
            # remove contemporary adult biomass with recruits, should be S2 biomass -> S1 recs
            select(-c(NOI,
                      ENSO,
                      NPGO,
                      #All_Copepods, # ~same as calanoid copepods
                      #euphausiids, # too large for larval mouth gape
                      anchBioSmrySeas1,
                      sardBioSmrySeas1,
                      anchBioSmrySeas2, # leave out biomass since not fit well
                      sardBioSmrySeas2,
                      # NCOPspring,
                      # SCOPspring, # summer copepod index had higher loadings
                      PDOsummer, # lower loading than spring, may want to try a lag
                      PDOspring,
                      # BEUTI_33N, # oceanography at 39N had highest loadings
                      # OC_LUSI_33N,
                      # OC_LUSI_36N,
                      # OC_STI_33N,
                      # OC_STI_36N,
                      copMeanSize,
                      copPropBioSize,
                      naupMeanSize,
                      naupPropBioSize,
                      ZL_NorCal,
                      ZL_SoCal,
                      age1SprAnchmeanWAA)) # not enough data in time windowselect(-c(#sprCalCOFILarvalSardine,

datNames <- names(allDat)[-1]

# transpose for MARSS formatting
allDat <- allDat %>% select(-year) %>% t()


# Create a custom R obs error matrix assuming each data source has it's own common error
Rcustom <- matrix(list(0),length(datNames),length(datNames)) 
diag(Rcustom) <- c("HCI",
                   "COP", "COP", "COP", "COP",
                   "BEUTI", "BEUTI", 
                   "CUTI", "CUTI",
                   "LUSI", "LUSI", "LUSI", 
                   "STI", "STI", "STI",
                   "RREAS",
                   "SSWI", "SSWI",
                   "CalCOFI", "CalCOFI", "CalCOFI",
                   "RREAS",
                   "WAA", "WAA",
                   "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS", "PRPOOS", #"PRPOOS", "PRPOOS",
                   "NEMURO", "NEMURO",
                   "SDM", "SDM", "TIME", "TIME", "SDM", "SDM",
                   "anchRec", 
                   "sardRec",
                   # "anchBio",
                   # "sardBio",
                   "SST", "SST",
                   "Transp", "Transp", "Transp", "Transp",
                   "Alb", "Hake")

overallDFA <- MARSS(y = allDat, 
                    form = "dfa",
                    method = "BFGS",
                 # control = list(maxit = 10000,
                 #                conv.test.slope.tol = 0.1,
                 #                allow.degen = TRUE),
                 inits = list(x0 = matrix(1, 1, 1)),
                 z.score = TRUE,
                 model = list(#R = "diagonal and equal", # observation errors are the same
                               # R = "diagonal and unequal", # observation errors independent
                               # R = "equalvarcov", # observation errors equal and covars equal
                               # R = "unconstrained", # all observation errors independent
                               R = Rcustom,
                               m = 4) # number of latent processes
)

# save(overallDFA, file = "marssFit_1980to2019_noBio_3trend_Rcustom.RData")
# save(overallDFA, file = "marssFit_1980to2019_noBio_4trend_Rcustom.RData")

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

# plot(c(1990:2019), y= overallDFA$states, type = "l")
# abline(h = 0)
# 
# plot loadings
# new df with coordinates
loadingsDF <- data.frame(index = datNames,
                         trend1 = Z.rot[, 1],
                         trend2 = Z.rot[, 2],
                         trend3 = Z.rot[, 3],
                         trend4 = Z.rot[, 4],
                         dummy0 = 0) %>%
                pivot_longer(cols = c(trend1, trend2, trend3, trend4), 
                             names_to = "Trend")

loadingsDF %>% ggplot(aes(y = index)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = value)) +
  facet_wrap(~Trend) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  theme_classic()

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(overallDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 
ggplot(data = subset(d, name=="model" & 
                       .rownames %in% c("sardRec", "anchRec", "sardLarv", 
                                        "anchLarv", "anchBioSmrySeas2", 
                                        "sardBioSmrySeas2"))) + 
  geom_point(aes(t, value)) + 
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(t, .fitted), col="blue") + 
  facet_wrap(~.rownames) + xlab("Time Step") + ylab("Count")



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
par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
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
  axis(1,  (0:dim(allDat)[2]) , 1980 + 0:dim(allDat)[2])
} # end i loop (trends)

# Projection Model -----------------------------------------------------------

# subset from 1980 to 2019
projDat <- datDFA %>% filter(year %in% 1980:2019) %>%
  # select only projectable vars and vars of interest
  select("year", "BEUTI_33N", "BEUTI_39N", "CUTI_33N", "CUTI_39N", "OC_LUSI_33N",
         "OC_LUSI_36N", "OC_LUSI_39N", "OC_STI_33N", "OC_STI_36N", "OC_STI_39N",
         "avgSSWIspring", "avgSSWIsummer", "sardLarv", "anchLarv", "anchYoY",
         "ZM_NorCal", "ZM_SoCal", "sardSpawnHab", "anchSpawnHab", "daysAbove5pct",
         "daysAbove40pct", "sardNurseHab", "anchNurseHab", "anchRec", "sardRec",
         "springSST", "summerSST", "avgNearTransspring", "avgNearTranssummer",
         "avgOffTransspring", "avgOffTranssummer") 

datNames <- names(projDat)[-1]

# transpose for MARSS formatting
projDat <- projDat %>% select(-year) %>% t()

# Create a custom R obs error matrix assuming each data source has it's own common error
RcustProj <- matrix(list(0),length(datNames),length(datNames)) 
diag(RcustProj) <- c("BEUTI", "BEUTI", 
                     "CUTI", "CUTI",
                     "LUSI", "LUSI", "LUSI", 
                     "STI", "STI", "STI",
                     "SSWI", "SSWI",
                     "CalCOFI", "CalCOFI",
                     "RREAS",
                     "NEMURO", "NEMURO",
                     "SDM", "SDM", "TIME", "TIME", "SDM", "SDM",
                     "anchRec",
                     "sardRec",
                     "SST", "SST",
                     "Transp", "Transp", "Transp", "Transp")

projectDFA <- MARSS(y = projDat, 
                    form = "dfa",
                    method = "BFGS",
                    # control = list(maxit = 10000,
                    #                conv.test.slope.tol = 0.1,
                    #                allow.degen = TRUE),
                    inits = list(x0 = matrix(1, 1, 1)),
                    z.score = TRUE,
                    model = list(#R = "diagonal and equal", # observation errors are the same
                      # R = "diagonal and unequal", # observation errors independent
                      # R = "equalvarcov", # observation errors equal and covars equal
                      # R = "unconstrained", # all observation errors independent
                      R = RcustProj,
                      m = 4) # number of latent processes
)

# save(projectDFA, file = "marssFit_1980to2019_ProjDFA_3trend_Rcustom.RData")
# save(projectDFA, file = "marssFit_1980to2019_ProjDFA_4trend_Rcustom.RData")

# Look at factor loadings
# get the inverse of the rotation matrix 
Z.est <- coef(projectDFA, type = "matrix")$Z 
H.inv <- 1 
if (ncol(Z.est) > 1){
  H.inv <- varimax(coef(projectDFA, type = "matrix")$Z)$rotmat
} 

# rotate factor loadings 
Z.rot <- Z.est %*% H.inv 
# rotate trends 
trends.rot <- solve(H.inv) %*% projectDFA$states

# Add CIs to marssMLE object 
projectDFA <- MARSSparamCIs(projectDFA) 
# Use coef() to get the upper and lower CIs 
Z.low <- coef(projectDFA, type = "Z", what = "par.lowCI") 
Z.up <- coef(projectDFA, type = "Z", what = "par.upCI") 
Z.rot.up <- Z.up %*% H.inv 
Z.rot.low <- Z.low %*% H.inv 
df <- data.frame(est = as.vector(Z.rot), 
                 conf.up = as.vector(Z.rot.up), 
                 conf.low = as.vector(Z.rot.low) )

# plot(c(1990:2019), y= projectDFA$states, type = "l")
# abline(h = 0)
# 
# plot loadings
# new df with coordinates
loadingsDF <- data.frame(index = datNames,
                         trend1 = Z.rot[, 1],
                         trend2 = Z.rot[, 2],
                         trend3 = Z.rot[, 3],
                         trend4 = Z.rot[, 4],
                         dummy0 = 0) %>%
  pivot_longer(cols = c(trend1, trend2, trend3, trend4), 
               names_to = "Trend")

loadingsDF %>% ggplot(aes(y = index)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = value)) +
  facet_wrap(~Trend) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  theme_classic()

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(projectDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 
ggplot(data = subset(d, name=="model" & 
                       .rownames %in% c("sardRec", "anchRec", "sardLarv", 
                                        "anchLarv", "anchBioSmrySeas2", 
                                        "sardBioSmrySeas2"))) + 
  geom_point(aes(t, value)) + 
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(t, .fitted), col="blue") + 
  facet_wrap(~.rownames) + xlab("Time Step") + ylab("Count")



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
par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
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
  axis(1,  (0:dim(allDat)[2]) , 1980 + 0:dim(allDat)[2])
} # end i loop (trends)

# Plots for ECCWO poster --------------------------------------------------

# fits to data from pg 137 in MARSS User Guide
alpha <- 0.05 
d <- residuals(overallDFA, type = "tT") 
d$up <- qnorm(1- alpha / 2) * d$.sigma + d$.fitted 
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted 

# combine trend states with variables for plotting with model fits
ts.trends <- data.frame(trendStates = ts.trends,
                        t = 1980:2019)

d %>% filter(name=="model",
             # .rownames %in% c("NPGO", "anchSpawnHab",
             #                  "SCOPsummer", "sprCalCOFISouthernMesopels",
             #                  "anchRec", "sprCalCOFILarvalAnchovy")) %>%
             .rownames %in% c("sardRec", "anchRec", "sardLarv", "anchLarv")) %>%
  mutate(t = t+1979) %>%
  full_join(y = ts.trends, by = "t") %>%
  ggplot(aes(x = t)) + 
  geom_point(aes(y = value)) + 
  geom_ribbon(aes(ymin = lo, ymax = up), linetype = 2, alpha = 0.2) + 
  geom_line(aes(y = .fitted), col="blue") + 
  # geom_line(aes(y = trendStates), col = "darkgreen") +
  facet_wrap(~.rownames, ncol = 1) + xlab("Year") + ylab("Variable Anomaly") +
  theme_minimal() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

# plot loadings
# new df with coordinates
loadingsDF <- data.frame(vals = Z.rot, 
                         index = rownames(Z.rot),
                         dummy0 = 0)

# order loadings by magnitude and arrange for plotting
varArrang <- loadingsDF %>% arrange(vals.1) %>% pull(index)
# leave out response variables
varArrang <- varArrang[-which(varArrang %in% c("sardRec", "anchRec", 
                                               "sardLarv", "anchLarv"))]
loadingsDF <- loadingsDF %>% mutate(index = factor(index, 
                                                   level = c("sardRec", "anchRec", 
                                                             "sardLarv", "anchLarv",
                                                      varArrang)))

myCols <- c("#F8766D", "black","#FFB000",   "#619CFF", "#00BA38")
names(myCols) <- levels(c("Foraging", "Interest Var", "Preconditioning",  "Predation","Temperature"
                          ))

test1 <- loadingsDF %>% pivot_longer(cols = grep("vals.", names(loadingsDF), value = TRUE),
                            names_to = "trend",
                            names_prefix = "vals.",
                            values_to = "vals") %>%
  mutate(vals = case_when(abs(vals) < 0.05 ~ 0,
                          TRUE ~ vals),
         colCode = case_when(index %in% c("age1SprSardmeanWAA", "meanSSBwt",
                                          "NCOPspring", "NCOPsummer",                  
                                          "SCOPspring", "SCOPsummer",
                                          "ZM_NorCal", "ZM_SoCal") ~"Preconditioning",#  "#FFB000",
                             index %in% c("HCI", "sardSpawnHab", "anchSpawnHab",                
                                          "daysAbove5pct", "daysAbove40pct",
                                          "springSST", "summerSST") ~ "Temperature", #"#00BA38",
                             index %in% c("BEUTI_33N", "BEUTI_39N", "CUTI_33N",
                                          "CUTI_39N", "OC_LUSI_33N","OC_LUSI_36N",
                                          "OC_LUSI_39N", "OC_STI_33N", "OC_STI_36N",
                                          "OC_STI_39N", "swfscRockfishSurv_Myctophids",
                                          "avgSSWIspring", "avgSSWIsummer",
                                          "mesopelLarv", "copBio", "naupBio",
                                          "sardNurseHab", "anchNurseHab",
                                          "avgNearTransspring", "avgNearTranssummer",
                                          "avgOffTransspring", "avgOffTranssummer") ~ "Foraging",#"#F8766D",
                             index %in% c("albacore", "hake") ~ "Predation",#"#619CFF",
                             TRUE ~ "Interest Var" ),
         colCode = as.factor(colCode)) #%>% #"black")) %>%
  #filter(abs(vals) > 0.1) 
test1 %>%
  ggplot(aes(y = index, color = colCode)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = vals,
                   linewidth = 4)) +
  scale_color_manual(values = myCols) +
  # scale_color_manual(values = c("#FFB000", "#00BA38", "#F8766D", "#619CFF", "black"),
  #                    labels = c("Preconditioning", "Temperature",
  #                               "Foraging", "Predation", "Interest Var")) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  geom_hline(yintercept = 4.5, color = "black") +
  theme_classic() +
  facet_wrap(~trend, nrow = 1) #+
  # theme(legend.position = "none")

loadingsDF %>% filter(index %in% c("sardRec", "anchRec", 
                                   "sardLarv", "anchLarv"))

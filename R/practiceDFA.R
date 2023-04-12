# Practice DFA using CalCurrent IEA data
# Created: 1/10/2023, Robert Wildermuth

library(tidyverse)
library(MARSS)
library(lubridate)


# get data files to read in and collate

# data file path
datPath <- "C:/Users/r.wildermuth/Documents/FutureSeas/RecruitmentIndex/DFA_data/"
datFiles <- list.files(datPath)

datFiles <- datFiles[-which(datFiles %in% c("Anchovy_Sardine_len_wt_age.csv", 
                                            "AtlantisModel_C_biomass_data_31Mar2021.csv",
                                            "SDMoutput",
                                            "anchovySSfiles"))]

# This won't collate on year-month before joining dfs
# df <-
#   list.files(path = datPath, pattern = "*.csv") %>% 
#   map_dfr(~read_csv(paste0(datPath,.),
#                     col_types = c("T","d"),
#                     skip = 1))
# df

# df names
datNames <- sub("_dl20230109.csv", "", datFiles)
datNames <- sub("cciea_", "", datNames)

for(i in 1:length(datFiles)) {
  assign(datNames[i],
         read_csv(paste0(datPath,
                          datFiles[i]),
                  col_names = c("time", datNames[i]),
                  col_types = c("T","d"),
                  skip = 2))
}

recrDat <- lapply(mget(datNames), 
                  function(x) {x %>% mutate(year = year(time),
                                            month = month(time)) %>%
                                  select(-time)})
recrDat <- recrDat %>% reduce(full_join, by = c("year", "month"))

# Note: going to need to figure out how to deal with annual vs monthly datasets

# for now just use January (1) values for practice
# also narrow to where we have larval data
samplDat <- recrDat %>% filter(month == 1, !is.na(sprCalCOFILarvalAnchovy)) %>% print(n = 23)
# mostly complete dataset

# transpose for MARSS formatting
samplDat <- samplDat %>% select(-year, -month) %>% t()

# do simple DFA fit

simpleDFA <- MARSS(y = samplDat, 
                   form = "dfa",
                   control = list(maxit = 1000,
                                  allow.degen = TRUE),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(R = "diagonal and unequal", # observation errors independent
                                m = 1) # one latent process
                   )

plot(c(1998:2019,2021), y= simpleDFA$states, type = "l")
abline(h = 0)

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

# look at only non-larval data --------------------------------------------

# for now just use January (1) values for practice
# also narrow to where we have larval data
envtlDat <- recrDat %>% filter(month == 1, year >= 1967) %>% print(n = 125)
# mostly complete dataset

# transpose for MARSS formatting
envtlDat <- envtlDat %>% select(-year, -month, 
                                -sprCalCOFILarvalAnchovy, 
                                -sprCalCOFILarvalSardine) %>% t()

# fit DFA 
envtlDFA <- MARSS(y = envtlDat, 
                   form = "dfa",
                   control = list(maxit = 1000,
                                  allow.degen = TRUE),
                   inits = list(x0 = matrix(1, 1, 1)),
                   z.score = TRUE,
                   model = list(R = "diagonal and unequal", # observation errors independent
                                m = 1) # one latent process
)

plot(c(1967:2022), y= envtlDFA$states, type = "l")
abline(h = 0)

# plot loadings
# new df with coordinates
envtlLoads <- data.frame(index = datNames[-(16:17)],
                         vals = envtlDFA$par$Z, 
                         dummy0 = 0)

envtlLoads %>% ggplot(aes(y = index)) +
  geom_segment(aes(x = dummy0,
                   yend = index,
                   xend = vals)) +
  labs(x = "Loadings", y = "Index") +
  geom_vline(xintercept = 0, color = "grey") +
  theme_classic()
# RW: Pretty similar patterns for loadings compared to shorter, larval dataset w/ different magnitudes
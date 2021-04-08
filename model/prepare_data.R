#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(lubridate)

# import population estimate
pop_est <- read_csv("data/pop_est.csv")

# import August 2016 to combine
august_2016_counts <- read_csv("data/august_2016_counts.csv")

# import betas
betas <- read_csv("data/betas.csv")





#==========
#========== Calculate mean detection probabilites
#==========

# average over day / night categories
base_beta <- {betas %>%
  filter(name %in% c("(Intercept)","day_nightN","size_classlarge:day_nightN")) %>%
  summarize(base_beta = mean(mean))}$base_beta

# beta by station
beta_station <- betas %>%
  filter(!(name %in% c("(Intercept)","day_nightN","size_classlarge:day_nightN",
                       "size_classlarge"))) %>%
  mutate(station = strsplit(name, "station") %>% map_chr(~.x[2])) %>%
  select(station, mean) %>%
  rename(station_beta = mean) %>%
  bind_rows(tibble(station = "124",
                   station_beta = 0))

# beta by size_class
beta_size <- betas %>%
  filter(name %in% c("size_classlarge")) %>%
  select(mean) %>%
  rename(size_beta = mean) %>%
  mutate(size_class = "large") %>%
  bind_rows(tibble(size_class = "small",
                   size_beta = 0))

# 
august_2016_est <- august_2016_counts %>%
  full_join(beta_station) %>%
  full_join(beta_size) %>%
  mutate(alpha = base_beta + station_beta + size_beta,
         phi = exp(alpha) / (1 + exp(alpha)),
         est = (count / 5) / phi,
         date = as_date(paste(year, month, "01", sep = "-"))) %>%
  group_by(date, basin, size_class) %>%
  summarize(mean = mean(est)) %>%
  ungroup()

  





#==========
#========== Prepate data
#==========

# prepare data
data <- pop_est  %>%
  bind_rows(august_2016_est) %>%
  mutate(basin = factor(basin, levels = c("south", "north")),
         stage = factor(size_class, levels = c("small","large")),
         state = paste(basin, stage, sep = "_"),
         time = (as.numeric(date) - min(as.numeric(date)))/365,
         scale = ifelse(basin == "south", 2, 1),
         season = ifelse(lubridate::month(date) %in% c(6, 7), 
                         "summer",
                         ifelse(lubridate::month(date) %in% c(8, 9),
                                "winter",
                                NA)),
         mean_scale = scale*mean,
         sd = scale * sd) %>%
  select(state, stage, basin, date, season, time, mean, mean_scale) %>%
  arrange(date, stage) 

# mean densities and abundance
data %>%
  group_by(basin, stage) %>%
  summarize(density = median(mean),
            abundance = median(mean_scale))

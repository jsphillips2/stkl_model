#==========
#========== Load packages
#==========

library(tidyverse)
library(cowplot)
library(lemon)
library(WaveletComp)




#==========
#========== Plot defaults
#==========

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,-4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,10,0,0)),
                  axis.title.x = element_text(margin = margin(10,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

# define basin colors
basin_colors <- c("dodgerblue","gray20")
names(basin_colors) <- c("south","north")

# define datebreaks / limits
date_breaks <- lubridate::as_date(c("1990-06-01",
                                    "2000-06-01",
                                    "2010-06-01",
                                    "2020-06-01"))
date_limits <-  lubridate::as_date(c("1990-08-01",
                                     "2020-08-01"))


# define datebreaks / limits
year_breaks <- c(1990, 2000, 2010, 2020)
year_limits = c(1990, 2020)



#==========
#========== Import data
#==========

# main output
out_in <- readRDS("output/fit_full.rds")

# annual projections
# proj_sum <-read_rds(paste0("analysis/demographic_model/model/output/",name,"/proj_sum.rds"))
# n_sum <- proj_sum$n_sum %>%
#   mutate(year = (1991:2020)[time])
# sens_sum <- proj_sum$sens_sum %>%
#   mutate(year = (1991:2020)[time])
# aa_sum <- proj_sum$aa_sum %>%ungroup()
# aa_cont_sum <- proj_sum$aa_cont_sum %>% ungroup()




#==========
#========== Additional preparatinos
#==========

# extract data 
data_prep <- out_in$data_prep %>%
  mutate(stage = factor(stage, 
                        levels = c("small","large"),
                        labels = c("juvenile","adult")))

# extract data passed to stan
data_list <- out_in$data_list

# extract specifications
iter <- out_in$mcmc_specs$iter
chains <- out_in$mcmc_specs$chains

# extract full fit
fit <- out_in$fit

# extract posterior summary
fit_summary <- out_in$fit_summary

# data frame for matching stages and basins to id's
state_match <- data_prep %>% 
  tidyr::expand(nesting(state, basin, stage))  %>%
  arrange(basin, stage) %>%
  mutate(st = row_number(),
         state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "adult.south",
                                   "juvenile.north",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "adult\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nnorth")))

# vector of dates for matching with time id's
date_match <- data_prep$date %>% unique()

annual_output <-read_rds("analysis/annual_output.rds")

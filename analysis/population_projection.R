#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(Matrix)
library(matrixcalc)
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=18),
                  legend.text = element_text(size=18),
                  axis.text=element_text(size=18, color="black"),
                  axis.title = element_text(size = 18),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))



# import model
out_in <- read_rds("output/fit_full.rds")

# extract data
data_prep <- out_in$data_prep
data_list <- out_in$data_list
fit <- out_in$fit 
fit_summary <- out_in$fit_summary 
pj <- data_list$pj
nt <- data_list$nt
n <- data_list$n
b <- data_list$b
dates <- unique(data_prep$date)
years_all <- lubridate::year(dates)
years <- years_all[1:(length(years_all) - 2)]

# paramter names for extraction
vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rt") |
           str_detect(fit_summary$var, "rf") |
           str_detect(fit_summary$var, "qt") |
           str_detect(fit_summary$var, "qf") |
           str_detect(fit_summary$var, "gt") |
           str_detect(fit_summary$var, "gf") |
           str_detect(fit_summary$var, "dt") |
           str_detect(fit_summary$var, "df") | 
           str_detect(fit_summary$var, "x0"))}$var

# extract fit
extract_full <-  rstan::extract(fit, pars = vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(vars) %>%
  sample_n(500) %>%
  mutate(id = row_number()) %>%
  gather(var, val, -id) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])))




#==========
#========== Iterate
#==========

# iterations of MCMC
ids <- unique(extract_full$id)

setup = list(pj = pj,
             nt = nt,
             n = n,
             b = b,
             dates = dates,
             years_all = years_all,
             years = years,
             ids = ids,
             theta_names = theta_names)

start_time <- Sys.time()
annual_output <- parallel::mclapply(ids, function(id_){
  
  # extract data corresponding to iteration
  extract_ <- extract_full %>%
    filter(id == id_)
  
  # project annual dynamics
  annual_proj <- annual_proj_fn(extract_ = extract_, 
                                pj_ = pj, 
                                nt_ = nt, 
                                n_ = n, 
                                b_ = b, 
                                years_ = years, 
                                theta_names_ = theta_names) 
  
  
  # project population size derivatives over 1 time step
  dX <- dX_fn(annual_proj_ = annual_proj, 
              ip_ = 1)
  
  # project population size derivatives over many time steps (approx. asymptotic)
  # dX_asym <- dX_fn(annual_proj_ = annual_proj, 
  #             ip_ = 100)
  
  # caluclate transient growth rate, sensitivity, and elasticity
  sens <- sens_fn(dX)
  
  # caluclate asymptotic growth rate, sensitivity, and elasticity 
  # sens_asym <- sens_fn(dX_asym)
  
  return(list(setup = setup,
              annual_proj = annual_proj,
              dX = dX,
              # dX_asym  = dX_asym,
              sens = sens
              # sens_asym = sens_asym
              ))
  
})
end_time <- Sys.time()
end_time - start_time

# write_rds(annual_output, "analysis/annual_output.rds")

#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(rstan)
library(Matrix)
library(lubridate)

# rstan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# import data
data <- read_csv("data/hornsili_cpue_clean.csv") %>%
  mutate(basin = factor(basin, 
                        levels = c("south", 
                                   "north")),
         stage = factor(stage, 
                        levels = c("small",
                                   "large")),
         state = factor(state,
                        levels = c("south_small",
                                   "south_large",
                                   "north_small",
                                   "north_large")),
         station = factor(station,
                          levels = c("23",
                                     "27",
                                     "41",
                                     "44",
                                     "135",
                                     "DN",
                                     "128",
                                     "124")))

#=========================================================================================





#=========================================================================================
#========== Define demographic matrices
#=========================================================================================

# number of stages
stages <- data$stage %>% unique() %>% length()

# number of sites
sites <- data$basin %>% unique() %>% length()

# function for block-diagonal matrix from submatrix
block_fn <- function(v_, r_, c_) {
  as.matrix(Reduce("bdiag", 
                   lapply(1:c_, function(x){
                     matrix(v_, nrow = r_, ncol = r_, byrow = T)
                   })))
}

# site-specific stage transitions (time-varyings)
qt_block <- block_fn(v_ = c(1, 0,
                            0, 1),
                     c_ = sites, r_ = stages)

# site-specific stage transitions (fixed)
qf_block <- block_fn(v_ = c(0, 0,
                            0, 0),
                     c_ = sites, r_ = stages)

# site-specific stage transitions (time-varyings)
gt_block <- block_fn(v_ = c(0, 0,
                            0, 0),
                     c_ = sites, r_ = stages)

# site-specific stage transitions (fixed)
gf_block <- block_fn(v_ = c(0, 0,
                            1, 0),
                     c_ = sites, r_ = stages)

# stage-specific site transitions (time-varying)
# diagonal should be 0
dt_block <- matrix(c(0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0),
                   nrow = stages * sites, ncol = sites * sites)


# stage-specific site transitions (fixed)
# diagonal should be 0
df_block <- matrix(c(0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0),
                   nrow = stages * sites, ncol = sites * sites)

# site-specific recruitment (time-varying)
rt_block <- block_fn(v_ = c(0, 1,
                            0, 0),
                     c_ = sites, r_ = stages)

# site-specific recruitment (fixed)
rf_block <- block_fn(v_ = c(0, 0,
                            0, 0),
                     c_ = sites, r_ = stages)

#=========================================================================================





#=========================================================================================
#========== Package data
#=========================================================================================

# set indeces
nc = stages
nt <- length(unique(data$date))
ns <- length(unique(data$station))
n <- length(unique(data$state))

# define positions in for filling matrices
pja <- lapply(list(rt_block, rf_block, 
                   qt_block, qf_block,
                   gt_block, gf_block, 
                   dt_block, df_block), 
              function(x_) {
                x1  = cbind(c(x_), 
                            rep(1:nrow(x_), ncol(x_)), 
                            rep(1:ncol(x_), each = nrow(x_)))
                x1[,1] = cumsum(x1[, 1])*x1[, 1] + 1
                return(x1)
              }) 
pj <- array(unlist(pja), dim = c(n * n, 3, length(pja)))

# extract dimensions for demographic rates
j <- lapply(pja, function(x){max(x[,1] - 1)}) %>% unlist()

# adjust so that development rate is constant between basins
pj[12,1,6] <- 2
j[6] <- 1

# extract CPUE
# extract CPUE
data_sort <- data %>% arrange(date, station, stage)
yy <- array(data_sort$cpue,
            dim = c(nc, 
                    ns,
                    nt))


# check
plot(yy[1, 1, ], type = "l", ylim = c(0, 400))
for (k in 2:5) {
  lines(yy[1, k, ] , type = "l", col = k)
}
ggplot(data = data_sort %>% filter(basin == "south", stage == "small") %>% 
         group_by(station) %>%
         arrange(date) %>%
         mutate(index = row_number()),
       aes(x = index,
           y = cpue,
           color = station))+
  geom_line()+
  ylim(c(0, 400))+
  theme_bw()

plot(yy[1, 6, ], type = "l", ylim = c(0, 400))
for (k in 7:8) {
  lines(yy[1, k, ] , type = "l", col = k)
}
ggplot(data = data_sort %>% 
         filter(basin == "north", stage == "small") %>% 
         group_by(station) %>%
         arrange(date) %>%
         mutate(index = row_number()),
       aes(x = index,
           y = cpue,
           color = station))+
  geom_line()+
  ylim(c(0, 400))+
  theme_bw()

plot(yy[2, 1, ], type = "l", ylim = c(0, 400))
for (k in 2:5) {
  lines(yy[2, k, ] , type = "l", col = k)
}
ggplot(data = data_sort %>% filter(basin == "south", stage == "large") %>% 
         group_by(station) %>%
         arrange(date) %>%
         mutate(index = row_number()),
       aes(x = index,
           y = cpue,
           color = station))+
  geom_line()+
  ylim(c(0, 400))+
  theme_bw()

plot(yy[2, 6, ], type = "l", ylim = c(0, 400))
for (k in 7:8) {
  lines(yy[2, k, ] , type = "l", col = k)
}
ggplot(data = data_sort %>% filter(basin == "north", stage == "large") %>% 
         group_by(station) %>%
         arrange(date) %>%
         mutate(index = row_number()),
       aes(x = index,
           y = cpue,
           color = station))+
  geom_line()+
  ylim(c(0, 400))+
  theme_bw()

# scale CPUE by global meaan
y <- yy / mean(yy)

# define mapping of stages and classes to states
m <- data_sort %>% 
  select(basin, stage, state, station) %>%
  unique() %>%
  arrange(stage, basin, station) %>%
  mutate(val = as.numeric(state)) %>%
  select(-basin, -state) %>%
  pivot_wider(names_from = stage,
              values_from = val) %>%
  select(small, large) %>%
  as.matrix()

# standardize season to 6 months (0.5 years)
b <- rep(0.5, nt - 1)

# priors
p <- c(2, 2, 2, 2)

# scaling factor for relative basin area
s <- c(2, 2, 1, 1)

# package data
data_list <- list(n = n,
                  nc = nc,
                  nt = nt,
                  ns = ns,
                  m = m,
                  y = y,
                  j = j,
                  pj = pj,
                  b = b,
                  s = s,
                  p = p)

#=========================================================================================





#=========================================================================================
#========== Fit model
#=========================================================================================

# MCMC specifications
# chains <- 4
# iter <- 4000
# adapt_delta <- 0.9
# max_treedepth <- 11

# fit model
# start_time <- Sys.time()
# fit <- stan(file = "model/multi_state.stan",
#             data = data_list,
#             seed=2e3,
#             chains = chains,
#             iter = iter,
#             control = list(adapt_delta = adapt_delta, 
#                            max_treedepth = max_treedepth))
# end_time <- Sys.time()
# end_time - start_time

# summarize fit
# fit_summary = rstan::summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
#   {as_tibble(.) %>%
#       mutate(var = rownames(rstan::summary(fit)$summary))}

# package for export
# out <- list(data_sort = data_sort,
#             mat_list = list(rt_block, rf_block, qt_block, qf_block, dt_block, df_block),
#             data_list = data_list,
#             mcmc_specs = list(chains = chains,
#                               iter = iter,
#                               adapt_delta = adapt_delta,
#                               max_treedepth = max_treedepth),
#             fit = fit,
#             fit_summary = fit_summary)
# # export
# write_rds(out,  "output/fit_no_move.rds")

#=========================================================================================





#=========================================================================================
#========== Examine overall fit
#=========================================================================================

# data frame for matching states
state_match <- data %>%
  select(basin, stage, state) %>%
  unique() %>%
  mutate(st = as.numeric(state))

# vectors matching dates
date_match <- data$date %>% unique()

# estimated relative population size
x_clean <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "x"), !str_detect(fit_summary$var, "x0")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) %>%
  full_join(state_match) %>%
  select(basin, state, stage, date, name, lo, mi, hi)

# plot
x_clean %>%
  ggplot(aes(x = date,
             y = mi,
             color = basin, fill = basin))+
  facet_grid(stage~basin)+
  geom_line()+
  geom_ribbon(data = x_clean, 
              aes(x = date, 
                  ymin = lo, 
                  ymax = hi),
              linetype = 0, alpha = 0.2)+
  scale_color_manual("",values = c("dodgerblue","gray35"), guide = F)+
  scale_fill_manual("",values = c("dodgerblue","gray35"), guide = F)+
  scale_x_date("Date",
               breaks = lubridate::as_date(c("1994-06-01","2003-06-01","2012-06-01")),
               labels = c("1994","2003","2012"))

#=========================================================================================






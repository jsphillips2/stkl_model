#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(rstan)
library(Matrix)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# import and process data
source("model/prepare_data.R")

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

# prepare data
data_prep <- data %>%
  arrange(date, basin, stage)


# set indeces
n = stages * sites
nt <- length(unique(data_prep$date))

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

# extract abundance
yy <- matrix(data_prep$mean_scale,
            nrow = length(unique(data_prep$state)),
            ncol = length(unique(data_prep$date)))
y_scale <- mean(yy)
y <- yy / y_scale

# extract season
b <- {data_prep %>% 
  filter(stage == "small",basin == "north") %>%
  mutate(b = c(diff(time), 0)) %>%
  filter(time < max(time))}$b

# UPDATE: standardize season to 6 months (0.5 years)
b <- rep(0.5, length(b))

# priors
p <- c(2, 2, 2, 2)

# package data
data_list <- list(n = n,
                  nt = nt,
                  y = y,
                  j = j,
                  pj = pj,
                  b = b,
                  p = p)

#=========================================================================================





#=========================================================================================
#========== Fit model
#=========================================================================================

# models
models <- c("multi_state")
model <- models[1]

# model
model_path <- paste0("model/",model,".stan")


# MCMC specifications (for testing)
# chains <- 4
# iter <- 100
# adapt_delta <- 0.8
# max_treedepth <- 10

# MCMC specifications
# chains <- 4
# iter <- 8000
# adapt_delta <- 0.9
# max_treedepth <- 11

# fit model
# start_time <- Sys.time()
# fit <- stan(file = model_path, data = data_list, seed=2e3,
#             chains = chains, iter = iter,
#             control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))
# end_time <- Sys.time()
# end_time - start_time
# 
# fit_summary = rstan::summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
#   {as_tibble(.) %>%
#       mutate(var = rownames(rstan::summary(fit)$summary))}
# 
# # data frame for matching
# state_match <- data_prep %>% 
#   tidyr::expand(nesting(state, basin, stage))  %>%
#   arrange(basin, stage) %>%
#   mutate(st = row_number(),
#          state = factor(interaction(stage, basin),
#                         levels = c("juvenile.south",
#                                    "adult.south",
#                                    "juvenile.north",
#                                    "adult.north"),
#                         labels = c("juvenile\nsouth",
#                                    "adult\nsouth",
#                                    "juvenile\nnorth",
#                                    "adult\nnorth")))
# 
# # year match
# date_match <- data_prep$date %>% unique()
# 
# # X
# x_clean <- fit_summary %>%
#   select(var, `16%`, `50%`, `84%`) %>%
#   rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
#   filter(str_detect(fit_summary$var, "x"), !str_detect(fit_summary$var, "x0")) %>%
#   mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
#          st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#          time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
#          date = date_match[time]) %>%
#   full_join(state_match) %>%
#   select(basin, state, stage, date, name, lo, mi, hi)
# 
# # plot 
# data_prep %>%
#   ggplot(aes(color = basin, fill = basin))+
#   facet_grid(stage~basin)+
#   geom_point(aes(x = date, mean_scale / mean(yy)), size = 0.7, shape = 1, stroke = 0.3)+
#   geom_ribbon(data = x_clean, aes(x = date, ymin = lo, ymax = hi),
#               linetype = 0, alpha = 0.2)+
#   geom_line(aes(x = date, mean_scale / mean(yy)), size = 0.1)+
#   geom_line(data = x_clean, aes(x = date, y = mi), size = 0.4)+
#   scale_color_manual("",values = c("dodgerblue","gray35"), guide = F)+
#   scale_fill_manual("",values = c("dodgerblue","gray35"), guide = F)+
#   scale_x_date("Date", 
#                breaks = lubridate::as_date(c("1994-06-01","2003-06-01","2012-06-01")),
#                labels = c("1994","2003","2012"))
# 
# # package for export
# out <- list(data_prep = data_prep,
#             mat_list = list(rt_block, rf_block, qt_block, qf_block, dt_block, df_block),
#             y_scale = y_scale,
#             data_list = data_list,
#             mcmc_specs = list(chains = chains,
#                               iter = iter,
#                               adapt_delta = adapt_delta,
#                               max_treedepth = max_treedepth),
#             fit = fit,
#             fit_summary = fit_summary)
# 
# # export
# saveRDS(out,  "output/fit_no_move.rds")

#=========================================================================================
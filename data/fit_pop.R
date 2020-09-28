#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(runjags)
library(lubridate)

# set plot theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=16),
                  legend.text = element_text(size=16),
                  axis.text=element_text(size=16, color="black"),
                  axis.title = element_text(size = 16),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))






#==========
#========== Prepare data
#==========



### Primary data

# import abundance data
abundance <- read_csv("data/clean_data/counts89_20.csv") %>%
  filter(station != "BV") %>%
  mutate(basin = ifelse(basin == "east", "south", basin),
         basin = factor(basin, levels = c("north","south")),
         size_class = factor(size_class, levels = c("small","large")),
         day_night = factor(day_night, levels = c("D","N")),
         count = ifelse(year==2016 & is.na(count), 0, count),
         date = ymd(paste(year, month, 1, sep = "_")),
         date_z = (as.numeric(date) - min(as.numeric(date)))/(2*sd(as.numeric(date)))
  ) %>%
  arrange(basin, station, size_class, date, day_night)

# extract observations
y <- abundance$count %>% round(0)



### Define model matrix for detection probabilites
x <- model.matrix(~ station + size_class * day_night, 
                  abundance %>% select(station, size_class, day_night) %>% unique())

# indeces for model matrix
n_k = nrow(x)
n_j <- ncol(x)



### Define mappings

data_prep <- abundance %>%
  mutate(
    l_group = as.numeric(factor(paste(basin, size_class, date))),
    n_group = as.numeric(factor(paste(station, size_class, date))),
    p_group = as.numeric(factor(paste(station, 
                                      size_class, 
                                      day_night 
    ))))

# basin abundance (lambda) to station abundance (nu) 
lam_nu <- {data_prep %>% 
  group_by(n_group) %>%
  summarize(l_group = unique(l_group)) %>% 
  ungroup() %>%
  arrange(n_group)}$l_group
n_nu <- length(lam_nu)
n_lam <- max(lam_nu)

# station abundance (nu) to trap observations (y)
nu_y <- {data_prep %>% 
  select(n_group)}$n_group

# detection probability (phi) to trap observations (y)
phi_y <- {data_prep %>% 
    select(p_group)}$p_group



### Priors
u_lam <- 1000 # basin abundance mean
u_beta <- 0 # coefficient mean (separate value for intercept) 
s_beta <- 2 # coefficient sd (separate value for intercept)



### Package data
data_list = list(y = y, # observed trap counts
                 x = x, # model matrix for detection probabilities 
                 n_k = n_k, # number of rows in model matrix 
                 n_j = n_j, # number of coefficients 
                 n_nu = n_nu, # number of station estimates 
                 n_lam = n_lam, # number of basin estimates 
                 lam_nu = lam_nu, # mapping basin to station abundance 
                 nu_y = nu_y, # mapping station abundance to observed trap counts
                 phi_y = phi_y, # mapping of detection probability to observed trap counts
                 u_lam = u_lam, # basin abundance mean
                 u_beta = u_beta, # coefficient mean
                 s_beta = s_beta # coefficient sd
                 )





#==========
#========== Fit model
#==========



### Initial values

# basin abundance (mean across stations)
lp <- {data_prep %>%
    group_by(l_group) %>%
    summarize(count = mean(count))}$count

# station abundance (max across traps)
nmax <- {data_prep %>%
    group_by(n_group) %>%
    summarize(nmax = max(count))}$nmax

# set seed
seed = 2e4

# function for initial values
init_fn <- function(){
  list(
       beta = runif(n_j, u_beta - s_beta, u_beta + s_beta),
       lam = runif(lp, 
                   (lp - 0.5*lp)/(exp(-2)/(1+exp(-2))), 
                   (lp + 0.5*lp)/(exp(-2)/(1+exp(-2)))),
       nu = as.integer(round(runif(length(nmax), 
                             nmax/(exp(-2)/(1+exp(-2))), 
                             1.5*(nmax + 1)/(exp(-2)/(1+exp(-2)))))),
       .RNG.name="base::Wichmann-Hill", .RNG.seed=seed)
}



### Fit model

# file path for model
model_path <- "analysis/population_estimate/pop_est.txt"

# variables to monitor
# monitor <- c("u_lam","s_lam","lam","nu","beta","phi")
monitor <- c("lam","nu","beta","phi")


# MCMC specifications (for testing)
# n.chains <- 1
# adapt <- 10
# burnin <- 10
# sample <- 10

# MCMC specifications
# n.chains <- 4
# adapt <- 5000
# burnin <- 5000
# sample <- 5000

# fit
# start_time <- Sys.time()
# fit <- run.jags(model = model_path, data = data_list, inits = init_fn,
#                 n.chains = n.chains, burnin = burnin, sample = sample,
#                 adapt = adapt, monitor = monitor)
# end_time <- Sys.time()
# end_time - start_time

# export fit
# saveRDS(fit, paste0("analysis/population_estimate/fit.rds"))

# check diagnostic
psrf_check <- coda::gelman.diag(fit$mcmc)
psrf_clean <- psrf_check$psrf %>%
  as_tibble() %>%
  mutate(rowname =  rownames(psrf_check$psrf)) %>%
  rename(est  = `Point est.`,
         upp = `Upper C.I.`) %>%
  select(rowname, est, upp) %>%
  arrange(-est)

# summarize
options(mc.cores = parallel::detectCores()-2)
fit_sum <- fit$mcmc %>%
  parallel::mclapply(function(x_){
    d_ = x_ %>%
      as_tibble()
    return(d_)
  }) %>%
  bind_rows() %>%
  gather(rowname, val) %>%
  group_by(rowname) %>%
  summarize(mean = mean(val),
            sd = sd(val),
            lower95 = quantile(val, probs = 0.025),
            lower68 = quantile(val, probs = 0.16),
            median = quantile(val, probs = 0.5),
            upper68 = quantile(val, probs = 0.84),
            upper95 = quantile(val, probs = 0.975)) 

# export summary
# write_csv(fit_sum, paste0("analysis/population_estimate/fit.csv"))





#==========
#========== Examine fit
#==========



### Import data
# fit <- readRDS(paste0("analysis/population_estimate/fit.rds"))
# fit_sum <- read_csv("analysis/population_estimate/fit.csv")

### MCMC chains
bayesplot::mcmc_trace(x = fit$mcmc,
                      pars = fit_sum$rowname[str_detect(fit_sum$rowname, "beta\\[")])
# bayesplot::mcmc_trace(x = fit$mcmc, 
#                       pars = fit_sum$rowname[str_detect(fit_sum$rowname, "phi\\[")])
# bayesplot::mcmc_trace(x = fit$mcmc, pars = "beta[8]")



### Population estimates

# summarize
l_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "lam\\[")) %>%
  mutate(l_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower95, lower68, median, upper68, upper95, l_group) %>%
  full_join(data_prep %>% 
              expand(nesting(date, basin, size_class, l_group))) %>%
  arrange(size_class, basin, date)

# plot
l_fit %>%
  mutate(median = ifelse(basin == "south", 2 * median, median)) %>%
  ggplot(aes(date, median, color = size_class))+
  facet_wrap(~basin, nrow = 2)+
  geom_hline(yintercept = l_fit$median %>% mean, color = "gray50", size = 0.3)+
  geom_line(size = 0.5)+
  geom_point(size = 1)+
  scale_x_date("Date",breaks = c("1995-06-01","2005-06-01","2015-06-01") %>% as.Date, 
               labels = c(1995,2005,2015))+
  scale_y_continuous("Population density")+
  scale_color_manual("",values = c("gray30","dodgerblue"))+
  scale_fill_manual("",values = c("gray30","dodgerblue"))

# plot (log scale)
l_fit %>%
  ggplot(aes(date, median, color = size_class))+
  facet_wrap(~basin, nrow = 2)+
  geom_hline(yintercept = l_fit$median %>% mean, color = "gray50", size = 0.3)+
  geom_line(size = 0.5)+
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = size_class), alpha = 0.5, linetype = 0)+
  geom_point(size = 1)+
  scale_x_date("Date",breaks = c("1995-06-01","2005-06-01","2015-06-01") %>% as.Date, 
               labels = c(1995,2005,2015))+
  scale_y_continuous("Population density", trans = "log", breaks = c(1,10,100,1000))+
  scale_color_manual("",values = c("gray30","dodgerblue"))+
  scale_fill_manual("",values = c("gray30","dodgerblue"))

# mean detection probability
u_phi <- fit$mcmc %>%
  parallel::mclapply(function(x_){
    d_ = x_ %>%
      as_tibble() %>%
      select(colnames(x_)[str_detect(colnames(x_), "phi\\[")])
    return(d_)
  }) %>%
  bind_rows()  %>% 
  apply(1, mean) %>%
  median()

# combine estimates with observed means
comp <- data_prep %>%
  group_by(month, year, date, basin, size_class) %>%
  summarize(count = mean(count)/u_phi) %>%
  ungroup()  %>%
  full_join(l_fit %>%
              select(date, basin, size_class, median)) %>%
  gather(var, val, count, median)

# plot with observed data (scaled by detetction probability)
data_prep %>%
  mutate(count = count) %>%
  ggplot(aes(date, count))+
  facet_grid(basin~size_class)+
  geom_point(aes(shape = station), alpha = 0.5, size = 0.8, color = "dodgerblue")+
  geom_line(data = comp,
            aes(y = u_phi * val, color = var),
            size = 0.8)+
  scale_color_manual("", values = c("firebrick","black"), labels = c("mean",
                                                                     "estimate"))+
  scale_x_date("Date", breaks = as.Date(c("1995-01-01","2005-01-01","2015-01-01")),
               labels = c("1995","2005","2015"))+
  scale_y_continuous("Catch & Density", 
                     limits = c(0, 650), breaks = c(100, 300, 500))+
  scale_shape_manual("", values = 1:8)

# export estimate
# write_csv(l_fit, "analysis/population_estimate/pop_est.csv")



### betas (effects on logit scale)

# extract betas
betas <- fit_sum %>%
  filter(rowname %in% fit_sum$rowname[str_detect(fit_sum$rowname, "beta\\[")]) %>%
  mutate(id = strsplit(rowname, "\\[|\\]") %>% map_int(~as.integer(.x[2])),
         name = colnames(x)[id])
# plot
betas %>%
  filter(name != "(Intercept)") %>%
  ggplot(aes(name, median))+
  geom_point()+
  geom_hline(yintercept = 0, color = "gray50", size = 0.3)+
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0)+
  scale_x_discrete("")+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1),
        )




#### phis 
phis <- fit_sum %>%
  filter(rowname %in% fit_sum$rowname[str_detect(fit_sum$rowname, "phi\\[")]) %>%
  mutate(p_group = strsplit(rowname, "\\[|\\]") %>% map_int(~as.integer(.x[2]))) %>%
  full_join(data_prep %>%
              select(station, day_night, size_class, p_group) %>%
              unique())

phis %>%
  mutate(day_night = factor(day_night, levels = c("D","N"), labels = c("day","night")),
         station = factor(station, 
                          levels = c("23","27","41","44",
                                     "135","DN","128","124"))) %>%
  ggplot(aes(station, median, color = size_class))+
  facet_wrap(~day_night)+
  geom_hline(yintercept = 0.2, color= "gray70")+
  geom_point()+
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0)+
  scale_color_manual("",values = c("black","dodgerblue"))+
  scale_y_continuous("Detection probability \n(relative to basin mean)",
                     limits = c(0.05, 0.35))+
  scale_x_discrete("")+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1))


betas %>%
  filter(name != "(Intercept)") %>%
  ggplot(aes(name, median))+
  geom_point()+
  geom_hline(yintercept = 0, color = "gray50", size = 0.3)+
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0)+
  scale_x_discrete("")
# plot
betas %>%
  filter(val > 0) %>%
  ggplot(aes(val, mean))+
  facet_wrap(~var, scales = "free_x")+
  geom_hline(yintercept = 0, color = "gray50", size = 0.3)+
  geom_point(size = 0.5)+
  geom_errorbar(aes(ymin = lower95, ymax = upper95), width = 0)+
  scale_x_continuous("Predictor level", breaks = NULL)+
  scale_y_continuous("Effect (logit scale)", 
                     limits = c(-1.6, 1.6), breaks = c(-1,0,1))



### Population by station

# summarize
nu_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "nu\\[")) %>%
  mutate(n_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower95, lower68, median, upper68, upper95, n_group) %>%
  full_join(data_prep %>% 
              expand(nesting(date, station, basin, size_class, n_group))) %>%
  arrange(size_class, basin, station, date)

# plot
nu_fit %>%
  ggplot(aes(date, median, color = station))+
  facet_grid(basin~size_class)+
  geom_hline(yintercept = nu_fit$median %>% mean, color = "gray50", size = 0.3)+
  geom_line(size = 0.5)+
  geom_ribbon(aes(ymin = lower95, ymax = upper95, fill = station), alpha = 0.5, linetype = 0)+
  scale_x_date("Date",breaks = c("1995-06-01","2005-06-01","2015-06-01") %>% as.Date, 
               labels = c(1995,2005,2015))+
  scale_y_continuous("Population density")+
  scale_color_viridis_d("")+
  scale_fill_viridis_d("")


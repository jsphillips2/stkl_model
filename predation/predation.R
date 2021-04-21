#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(cowplot)
library(lemon)
library(rstan)
library(loo)
library(WaveletComp)
library(lubridate)
library(nlme)


# set cores
options(mc.cores = parallel::detectCores()-4)

# import trapping data and model population estimates



# import model objects
fit_full <- read_rds("output/fit_full.rds")

# define datebreaks / limits
date_breaks <- lubridate::as_date(c("1995-06-01",
                                    "2005-06-01",
                                    "2015-06-01"))
date_limits <-  lubridate::as_date(c("1990-08-01",
                                     "2020-08-01"))


# define datebreaks / limits
year_breaks <- c(1995, 2005, 2015)
year_limits = c(1990, 2020)

# extract data 
data_prep <- fit_full$data_prep %>%
  mutate(stage = factor(stage, 
                        levels = c("small","large"),
                        labels = c("juvenile","adult"))) %>%
  mutate(state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "juvenile.north",
                                   "adult.south",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nsouth",
                                   "adult\nnorth")))

# data frame for matching stages and basins to id's
state_match <- data_prep %>% 
  tidyr::expand(nesting(state, basin, stage))  %>%
  arrange(basin, stage) %>%
  mutate(st = row_number(),
         state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "juvenile.north",
                                   "adult.south",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nsouth",
                                   "adult\nnorth")))

# vector of dates for matching with time id's
date_match <- data_prep$date %>% unique()

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(t = 1,
                                       r = 1,
                                       b = 1,
                                       l = 1),
                  legend.margin = margin(t = 0,
                                         r = 0,
                                         b = 0,
                                         l = -4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black",family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,5,0,0)),
                  axis.title.x = element_text(margin = margin(5,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

#=========================================================================================




#=========================================================================================
#========== Lambdas
#=========================================================================================


stkl_lam <- read_csv("output/lam_sum.csv") %>%
  filter(var == "l_asym")
charr_lam <- read_csv("predation/charr/lambda_full.csv") %>%
  group_by(year) %>%
    summarize(lo = quantile(lam, probs = c(0.16), na.rm = T),
              mi = median(lam, na.rm = T),
              hi = quantile(lam, probs = c(0.84), na.rm = T))


ggplot(data  = stkl_lam,
       aes(x = year,
           y = mi))+
  geom_hline(yintercept = 1,
             linetype =2,
             size = 0.2,
             color = "gray50")+
  geom_line(color = "dodgerblue",
            linetype = 1)+
  geom_line(data = charr_lam %>% filter(year > 1990, year < 2019),
            color = "firebrick")+
  scale_y_continuous(trans = "log",
                     breaks = c(1/9, 1/3, 1, 3),
                     labels = c("1/9", "1/3", "1", "3"),
                     limits = c(1/9, 4))

  


wave_prep <- stkl_lam %>%
  select(year, mi) %>%
  rename(stkl = mi) %>%
  full_join(charr_lam %>%
              select(year, mi) %>%
              rename(charr = mi)) %>%
  na.omit() %>%
  mutate(stkl = log(stkl),
         charr = log(charr))

# wavelet transform
set.seed(34)
wave_r <- analyze.coherency(wave_prep,
                          my.pair = c("stkl","charr"),
                          loess.span = 0,
                          lowerPeriod = 2,
                          upperPeriod = 29,
                          dj = 1 / 20,
                          dt = 1,
                          make.pval = T,
                          n.sim = 2000)

wc.image(wave_r, 
         n.levels = 250,
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "", 
         periodlab = "period (days)")

wt.image(wave_r, 
         my.series = "stkl",
         n.levels = 250,
         timelab = "", 
         periodlab = "period (days)")

wt.image(wave_r, 
         my.series = "charr",
         n.levels = 250,
         timelab = "", 
         periodlab = "period (days)")

#=========================================================================================





#=========================================================================================
#========== Stkl lambda and charr abundance
#=========================================================================================

stkl_lam <- read_csv("output/lam_sum.csv") %>%
  filter(var == "l_asym")

charr_sites <- read_csv("predation/charr/site_data.csv") 

fit_sum <- read_csv("predation/charr/fit_sum.csv") 

data_list <- read_rds("predation/charr/data_list.rds") 

data <- read_csv("predation/charr/myvatn_char_clean_1986_2020.csv")


# extract scaling parameter
k <- data_list$k

# extract density estimate
x_fit <- fit_sum %>%
  filter(str_detect(.$var, "x\\[")) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) %>%
  # scale by k
  mutate(lo = k * lo,
         mi = k * mi,
         hi = k * hi) %>%
  select(year, stage, lo, mi, hi)



test  = stkl_lam %>%
  select(year, mi) %>%
  left_join(x_fit %>%
              # filter(stage %in% c("age 2", "age 4")) %>%
              group_by(year) %>%
              summarize(charr = sum(mi)) %>%
              ungroup() %>%
              mutate(charr = charr / max(charr))) %>%
  pivot_longer(cols = c(mi, charr)) %>%
  group_by(name) %>%
  mutate(value = log(value),
         z = (value - mean(value)) / sd(value))

ggplot(data  = test,
       aes(x = year,
           y = z,
           color = name))+
  geom_hline(yintercept = 1,
             linetype =2,
             size = 0.2,
             color = "gray50")+
  geom_line(linetype = 1)

ggplot(data  = test %>%
         select(-value) %>%
         pivot_wider(names_from = name,
                     values_from = z) %>%
         arrange(year) %>%
         mutate(lag = dplyr::lag(charr, n = 3L)),
       aes(x = lag,
           y = mi))+
  geom_point()+
  geom_smooth(method = "lm", se = F)
# +
#   scale_y_continuous(trans = "log",
#                      breaks = c(1/9, 1/3, 1, 3),
#                      labels = c("1/9", "1/3", "1", "3"),
#                      limits = c(1/9, 4))


#=========================================================================================





#=========================================================================================
#========== Stkl surv and charr abundance
#=========================================================================================

surv_sum <- read_csv("output/surv_sum_write.csv") %>%
  full_join(state_match) %>%
  select(-st)

charr_sites <- read_csv("predation/charr/site_data.csv") 

fit_sum <- read_csv("predation/charr/fit_sum.csv") 

data_list <- read_rds("predation/charr/data_list.rds") 

data <- read_csv("predation/charr/myvatn_char_clean_1986_2020.csv")


# extract scaling parameter
k <- data_list$k

# extract density estimate
x_fit <- fit_sum %>%
  filter(str_detect(.$var, "x\\[")) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) %>%
  # scale by k
  mutate(lo = k * lo,
         mi = k * mi,
         hi = k * hi) %>%
  select(year, stage, lo, mi, hi)


charr_dens <- x_fit %>%
  filter(stage %in% c("age 4+")) %>%
  group_by(year) %>%
  summarize(charr = sum(mi)) %>%
  ungroup() %>%
  mutate(charr = charr / max(charr))

test  = surv_sum %>%
  mutate(year = year(date)) %>%
  left_join(charr_dens %>%
              arrange(year) %>%
              mutate(charr_lag = lag(charr))) %>%
  mutate(z = (mi - mean(mi)) / sd(mi),
         charr_z = log(charr),
         charr_z = (charr_z - mean(charr_z)) / sd(charr_z),
         charr_lag_z = log(charr_lag),
         charr_lag_z = (charr_lag_z - mean(charr_lag_z)) / sd(charr_lag_z))


test_fit <- test %>%
  mutate(mi = log(mi / (1 - mi)),
         z = (mi - mean(mi)) / sd(mi),
         year_z = (year - mean(year)) / sd(year)) %>%
  select(date, year_z, state,stage,basin, charr_z, mi, z)

ggplot(data  = test_fit,
       aes(x = date,
           y = z))+
  facet_wrap(~state)+
  geom_line()+
  geom_line(aes(y = charr_z),
            color = "red")

ggplot(data  = test_fit,
       aes(x = charr_z,
           y = z))+
  facet_wrap(~state)+
  geom_point()+
  geom_smooth(method = "lm",
              se = F)

m <- gls(z ~ charr_z * state + year_z * state,
          correlation = corCAR1(form = ~ date | state),
          data = test_fit)
ttab <- summary(m, type = "marginal")$tTable %>% round(2)

coef_pos <- c(2 ,7, 8, 9)
summary(m, type = "marginal")$tTable[coef_pos, 1:2]

mm <- model.matrix(~charr_z + charr_z:state, 
                   data = test_fit %>%
                     expand(charr_z = 1,
                            state))[,c(2:5)]
vc <- as.matrix(vcov(m)[coef_pos, coef_pos])

apply(mm, 1, function(x){sqrt(t(x) %*% vc %*% x)})



charr_effect <- tibble(state = levels(test_fit$state),
                       est = c(mm %*% summary(m, type = "marginal")$tTable[coef_pos, 1]),
                       se = c(apply(mm, 1, function(x){sqrt(t(x) %*% vc %*% x)})))

ggplot(data = charr_effect,
       aes(x = state, 
           y = est))+
  geom_hline(yintercept = 0,
             linetype = 2)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = est -  se,
                    ymax = est +  se),
                width = 0,
                size = 1)

#=========================================================================================
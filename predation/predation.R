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
library(car)
library(AICcmodavg)


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
#========== Stkl surv and charr abundance
#=========================================================================================


# load data
surv_sum <- read_csv("output/surv_sum_write.csv") %>%
  full_join(state_match) %>%
  select(-st)

charr_sites <- read_csv("predation/charr/site_data.csv") 

fit_sum <- read_csv("predation/charr/fit_sum.csv") 

data_list <- read_rds("predation/charr/data_list.rds") 

data <- read_csv("predation/charr/myvatn_char_clean_1986_2020.csv")


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
  select(year, stage, mi)

# create charr index
charr_index <- x_fit %>%
  group_by(stage) %>%
  mutate(mi = log(mi),
         mi = (mi - mean(mi)) / sd(mi)) %>%
  pivot_wider(names_from = stage,
              values_from = mi) %>%
  rename(charr_1 = `age 1`,
         charr_2 = `age 2`,
         charr_3 = `age 3`,
         charr_4 = `age 4+`)


# combine stickleback survival with charr index
comb <- surv_sum %>%
  mutate(year = year(date),
         l_mi = log(mi / (1 - mi))) %>%
  group_by(year, state, basin, stage) %>%
  ungroup() %>%
  left_join(charr_index) %>%
  mutate(year_z = (year - mean(year)) / sd(year))


# fit model for juveniles
mod_juv <- gls(l_mi ~ year_z  * basin + 
                      charr_1 * basin + 
                      charr_2 * basin + 
                      charr_3 * basin + 
                      charr_4 * basin ,
              correlation = corCAR1(form = ~ date | basin),
              data = comb %>% filter(stage == "juvenile"))
summary(mod_juv)$tTable %>% round(3)
summary(update(mod_juv, .~. - charr_1 : basin - 
                              charr_2 : basin - 
                              charr_3 : basin - 
                              charr_4 : basin))$tTable %>% round(3)


# fit model for adults
mod_adl <- gls(l_mi ~ year_z  * basin + 
                      charr_1 * basin + 
                      charr_2 * basin + 
                      charr_3 * basin + 
                      charr_4 * basin,
               correlation = corCAR1(form = ~ date | basin),
               data = comb %>% filter(stage == "adult"))
summary(mod_adl)$tTable %>% round(3)
summary(update(mod_adl, .~. - charr_1 : basin - 
                              charr_2 : basin - 
                              charr_3 : basin - 
                              charr_4 : basin))$tTable %>% round(3)


# plot
pp <- ggplot(data = comb %>%
               select(date, l_mi, state, charr_4) %>%
               pivot_longer(cols = c(l_mi, charr_4)),
             aes(x = date,
                 y = value,
                 color = name))+
  facet_rep_wrap(~state)+
  geom_line()+
  scale_y_continuous("Value",
                     limits = c(-5.5, 5.5),
                     breaks = c(-4, 0, 4))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  scale_color_manual("",
                     values = c("blue", "black"),
                     labels = c("charr age 4+",
                                "stkl surv prob"))+
  theme(legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = "top",
        legend.text = element_text(margin = margin(l = -10)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.9, "lines"),
        plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "top", 
                    bottom='none',
                    gap = 0)

pp
# cairo_pdf(file = "predation/p_4.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# pp
# dev.off()


comb_plot <- comb %>%
  pivot_longer(cols = c(charr_1, charr_2, charr_3, charr_4)) %>%
  mutate(name = factor(name, 
                       levels = c("charr_1",
                                  "charr_2",
                                  "charr_3",
                                  "charr_4"),
                       labels = c("charr age 1",
                                  "charr age 2",
                                  "charr age 3",
                                  "charr age 4+")))

p_juv <- comb %>%
  expand(basin, 
         year_z = 0, 
         charr_1 = seq(from = -2, to = 2, length.out = 100),
         charr_2 = 0,
         charr_3 = 0,
         charr_4 = 0) %>%
  bind_rows(comb %>%
              expand(basin, 
                     year_z = 0, 
                     charr_1 = 0,
                     charr_2 = seq(from = -2, to = 2, length.out = 100),
                     charr_3 = 0,
                     charr_4 = 0)) %>%
  bind_rows(comb %>% expand(basin, 
                     year_z = 0,
                     charr_1 = 0,
                     charr_2 = 0,
                     charr_3 = seq(from = -2, to = 2, length.out = 100),
                     charr_4 = 0)) %>%
  bind_rows(comb %>%
              expand(basin, 
                     year_z = 0, 
                     charr_1 = 0,
                     charr_2 = 0,
                     charr_3 = 0,
                     charr_4 = seq(from = -2, to = 2, length.out = 100))) %>%
  filter(!(charr_1 == 0 && charr_2 == 0 && charr_3 == 0 && charr_4 == 0))
            
p_juv_pred <- predictSE.gls(mod_juv, newdata = p_juv, print.matrix = T)
p_juv$fit <- p_juv_pred[,"fit"]
p_juv$se <- p_juv_pred[,"se.fit"]

pp_juv <- p_juv %>%
  pivot_longer(cols = c(charr_1, charr_2, charr_3, charr_4)) %>%
  filter(value != 0) %>%
  mutate(name = factor(name, 
                       levels = c("charr_1",
                                  "charr_2",
                                  "charr_3",
                                  "charr_4"),
                       labels = c("charr age 1",
                                  "charr age 2",
                                  "charr age 3",
                                  "charr age 4+")))

pp2 <- ggplot(data = comb_plot %>% filter(stage == "juvenile"),
              aes(x = value,
                  y = l_mi,
                  color = basin))+
  facet_rep_wrap(~name)+
  geom_point(alpha = 0.5,
             size = 1)+
  geom_ribbon(data = pp_juv,
              aes(x = value,
                  ymin = fit - se,
                  ymax = fit + se,
                  fill = basin),
              alpha = 0.2,
              inherit.aes = F)+
  geom_line(data = pp_juv,
            aes(x = value,
                y = fit))+
  scale_x_continuous("Charr abundance index",
                     limits = c(-2.5, 2.5),
                     breaks = c(-2, 0, 2))+
  scale_y_continuous("Survival probability (logit)",
                     limits = c(-5.5, 5.5),
                     breaks = c(-4, 0, 4))+
  scale_color_manual("",values = c("dodgerblue","firebrick"))+
  scale_fill_manual("",values = c("dodgerblue","firebrick"))+
  theme(legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = "top",
        legend.text = element_text(margin = margin(l = -10)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.9, "lines"),
        plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "top", 
                    bottom='none',
                    gap = 0)
pp2
# cairo_pdf(file = "predation/p_juv.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# pp2
# dev.off()




p_adl <- comb %>%
  expand(basin, 
         year_z = 0, 
         charr_1 = seq(from = -2, to = 2, length.out = 100),
         charr_2 = 0,
         charr_3 = 0,
         charr_4 = 0) %>%
  bind_rows(comb %>%
              expand(basin, 
                     year_z = 0, 
                     charr_1 = 0,
                     charr_2 = seq(from = -2, to = 2, length.out = 100),
                     charr_3 = 0,
                     charr_4 = 0)) %>%
  bind_rows(comb %>% expand(basin, 
                            year_z = 0,
                            charr_1 = 0,
                            charr_2 = 0,
                            charr_3 = seq(from = -2, to = 2, length.out = 100),
                            charr_4 = 0)) %>%
  bind_rows(comb %>%
              expand(basin, 
                     year_z = 0, 
                     charr_1 = 0,
                     charr_2 = 0,
                     charr_3 = 0,
                     charr_4 = seq(from = -2, to = 2, length.out = 100))) %>%
  filter(!(charr_1 == 0 && charr_2 == 0 && charr_3 == 0 && charr_4 == 0))

p_adl_pred <- predictSE.gls(mod_adl, newdata = p_adl, print.matrix = T)
p_adl$fit <- p_adl_pred[,"fit"]
p_adl$se <- p_adl_pred[,"se.fit"]

pp_adl <- p_adl %>%
  pivot_longer(cols = c(charr_1, charr_2, charr_3, charr_4)) %>%
  filter(value != 0) %>%
  mutate(name = factor(name, 
                       levels = c("charr_1",
                                  "charr_2",
                                  "charr_3",
                                  "charr_4"),
                       labels = c("charr age 1",
                                  "charr age 2",
                                  "charr age 3",
                                  "charr age 4+")))

pp3 <- ggplot(data = comb_plot %>% filter(stage == "adult"),
              aes(x = value,
                  y = l_mi,
                  color = basin))+
  facet_rep_wrap(~name)+
  geom_point(alpha = 0.5,
             size = 1)+
  geom_ribbon(data = pp_adl,
              aes(x = value,
                  ymin = fit - se,
                  ymax = fit + se,
                  fill = basin),
              alpha = 0.2,
              inherit.aes = F)+
  geom_line(data = pp_adl,
            aes(x = value,
                y = fit))+
  scale_x_continuous("Charr abundance index",
                     limits = c(-2.5, 2.5),
                     breaks = c(-2, 0, 2))+
  scale_y_continuous("Survival probability (logit)",
                     limits = c(-3.5, 3.5),
                     breaks = c(-3, 0, 3))+
  scale_color_manual("",values = c("dodgerblue","firebrick"))+
  scale_fill_manual("",values = c("dodgerblue","firebrick"))+
  theme(legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = "top",
        legend.text = element_text(margin = margin(l = -10)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.9, "lines"),
        plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "top", 
                    bottom='none',
                    gap = 0)

pp3
# cairo_pdf(file = "predation/p_adl.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# pp3
# dev.off()


#=========================================================================================

#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(cowplot)


source("analysis/demographic_model/model/figures/fig_setup.R")


# extract posterior summary
fit_summary <- out_in$fit_summary

# # extract data
# data_prep <- out_in$data_prep %>%
#   mutate(stage = factor(stage, 
#                         levels = c("small","large"),
#                         labels = c("juvenile","adult")))

data_list <- out_in$data_list

# time match
# date_match <- data_prep$date %>% unique()






#==========
#========== Fig: survival
#==========



trans_season <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "qte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(st, mi) %>%
  rename(se = mi)



trans_rates <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "qt"),
         !str_detect(fit_summary$var, "qte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b1 = ifelse(lubridate::month(date) == 6, 0, 1)) %>%
  full_join(trans_season) %>%
  full_join(state_match)

surv_probs <- trans_rates %>%
  mutate(p = exp(-(1 / 12) * mi * exp(se * 0.5)),
         p_se = exp(-(1 / 12) * mi * exp(se * b1)))






p1 <- surv_probs %>%
  ggplot(aes(date, p, color = basin))+
  facet_wrap(~stage, nrow = 2)+
  geom_hline(yintercept = 0.5,size = 0.2, color = "gray50", linetype = 2)+
  geom_line(aes(y = p_se), size = 0.1)+
  geom_line(size = 0.4)+
  scale_color_manual("",values = c("dodgerblue","gray20"))+
  scale_y_continuous("Survival probability (monthly)", breaks = c(0.1, 0.5, 0.9), 
                     limits = c(0, 1))+
  scale_x_date("Year", 
               breaks = lubridate::as_date(c("1990-06-01",
                                             "2000-06-01",
                                             "2010-06-01",
                                             "2020-06-01")),
               labels = c("1990","2000","2010","2020"),
               limits = lubridate::as_date(c("1990-08-01",
                                             "2020-08-01")))+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", bottom='both')
p1

# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fit_surv.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p1
# dev.off()



#==========
#========== Fig: recruitment
#==========


rec_season <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "rte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(st, mi) %>%
  rename(se = mi)


rec <- fit_summary %>% select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "rt"),
         !str_detect(fit_summary$var, "rte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         basin = levels(state_match$basin)[
           str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))],
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b2 = data_list$b[time,2]) %>%
  full_join(rec_season %>%
              mutate(basin = levels(state_match$basin)[st])) %>%
  mutate(b1 = ifelse(lubridate::month(date) == 6, 0, 1),
         rec = b2 * mi * exp(se * b1),
         year = lubridate::year(date)) %>%
  group_by(year, basin) %>%
  summarize(rec = sum(rec))

p2 <- rec %>%
  ggplot(aes(year, rec, color = basin))+
  geom_line(size = 0.4)+
  scale_color_manual("",values = c("dodgerblue","gray20"))+
  scale_y_continuous("Per capita births (annual)",
                     limits = c(0, 13),
                     breaks = c(0, 4, 8, 12))+
  scale_x_continuous("Year",
                     limits = c(1990, 2020),
                     breaks = c(1990, 2000, 2010, 2020))+
  theme(strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", bottom='both')
p2

# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fit_rec.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p2
# dev.off()


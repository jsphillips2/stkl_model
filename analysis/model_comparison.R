#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(loo)


# import model
out_full <- read_rds("output/fit_full.rds")
out_null <- read_rds("output/fit_null.rds")
fit_full <- out_full$fit
fit_null <- out_null$fit
fit_summary_full <- out_full$fit_summary
fit_summary_null <- out_null$fit_summary






#==========
#========== LOO & Log0likd
#==========

# loo
loos <- lapply(c(fit_full,fit_null), function(x){
  
  fit_ <- x
  
  log_lik <- extract_log_lik(fit_, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik))
  
  loo <- loo(log_lik, r_eff = r_eff, cores = 10)
}) %>%
  set_names(c("full","null"))

# looic
comp = {loo_compare(loos[[1]],
                    loos[[2]]) %>%
    as_tibble() %>%
    as.matrix() %>%
    as_tibble() %>%
    mutate(model = c("full","null")[as.numeric(rownames(.))])} %>%
  select(model, looic) %>%
  mutate(looic = round(looic, 1),
         loo_dev = round(looic - looic[1], 1)) %>%
  full_join(tibble(model = c("full",
                              "null"),
                   loglik = c({fit_summary_full %>%
                                filter(var == "log_lik_sum")}$`50%`,
                              {fit_summary_null %>%
                                  filter(var == "log_lik_sum")}$`50%`)))







#==========
#========== Plot
#==========



# extract scaling parameter
y_scale <- out_null$y_scale

# extract fit
x_clean <- fit_summary_null %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary_null$var, "x"), !str_detect(fit_summary_null$var, "x0")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) %>%
  full_join(state_match) %>%
  select(basin, state, stage, date, name, lo, mi, hi)





#==========
#========== Plot null
#==========

# plot labels
labs <- x_clean %>%
  tidyr::expand(stage) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 10)

# plot 
p1 <-  ggplot(data = x_clean,
              aes(x = date, 
                  y = mi, 
                  color = basin, 
                  fill = basin))+
  facet_rep_wrap(~stage, 
                 nrow = 2)+
  geom_text(data = labs,
            aes(label = stage, 
                x = x, 
                y = y),
            color = "black", 
            size = 3,
            inherit.aes = F)+
  geom_ribbon(aes(ymin = lo, 
                  ymax = hi),
              alpha = 0.2,
              linetype = 0)+
  geom_point(data = data_prep,
             aes(y = mean_scale / y_scale),
             shape = 1,
             size = 0.8,
             stroke = 0.35)+
  geom_line(size = 0.35)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_fill_manual(name = "",
                    values = basin_colors,
                    guide = F)+
  scale_y_continuous(name = Relative~abundance,
                     limits = c(0, 12),
                     breaks = c(0, 3, 6, 9))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(-1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.15, 0.42),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')+
  guides(fill = guide_legend(override.aes = list(fill = NA)),
         color = guide_legend(override.aes = list(size = 0.5,
                                                  shape = NA)))

# examine plot
p1

# export
# cairo_pdf(file = "figures/figs/fig_fit_null.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# p1
# dev.off()
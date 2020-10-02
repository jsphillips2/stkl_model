#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")





#==========
#========== Calculate monthly survival
#==========

# extract variable names
dt_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "dt"))}$var

# extract survival
dt_full <- rstan::extract(fit, pars = dt_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(dt_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  sample_n(4000) %>%
  gather(var, dt, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b = data_list$b[time],
         basin = levels(state_match$basin)[st])

# calculate monthly survival
move_month <- dt_full %>%
  mutate(p = 1 - exp(-mean(b) * dt),
         p_se = 1 - exp(-b * dt),
         month = lubridate::month(date),
         season = factor(ifelse(month < 8, 0, 1),
                         levels = c(0, 1),
                         labels = c("summer", "winter"))) %>%
  select(chain_step, date, season, basin, p, p_se) %>%
  gather(var, val, p, p_se) %>%
  group_by(var, basin,season, date) %>%
  summarize(lo = quantile(val, prob = 0.16),
            mi = quantile(val, prob = 0.5),
            hi = quantile(val, prob = 0.84))




#==========
#========== Plot
#==========

# plot
p1 <- ggplot(data = move_month %>% 
               filter(var == "p"),
             aes(x = date,
                 y = mi,
                 color = basin))+
  geom_hline(yintercept = 0.5,
             size = 0.2,
             color = "gray50",
             linetype = 2)+
  geom_line(data = move_month %>% 
              filter(var == "p_se"),
            size = 0.2,
            alpha = 0.5)+
  geom_line(size = 0.4)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Dispersal~probability, 
                     breaks = c(0, 0.5, 1), 
                     limits = c(0, 1))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.15, 0.45),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')

# examine plot
p1

# export
# cairo_pdf(file = "figures/figs/fig_surv.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# p1
# dev.off()
  
  
  
  

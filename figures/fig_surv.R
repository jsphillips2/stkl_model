#==========
#========== Preliminaries
#==========

source("analysis/demographic_model/model/figures/fig_setup.R")





#==========
#========== Calculate monthly survival
#==========

# extract variable names
qs_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "qs"))}$var
qte_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "qte"))}$var
qt_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "qt"),
           !str_detect(fit_summary$var, "qte"))}$var

# extract scaling parameter
qs_full <- rstan::extract(fit, pars = qs_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(qs_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  gather(var, qs, -chain, -step, -chain_step) 

# extract season effects
qte_full <- rstan::extract(fit, pars = qte_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(qte_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  gather(var, qte, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  full_join(state_match) %>%
  full_join(qs_full %>%
              select(chain_step, qs))

# extract survival
qt_full <- rstan::extract(fit, pars = qt_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(qt_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  gather(var, qt, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b1 = data_list$b[time,1],
         b2 = data_list$b[time,2]) %>%
  full_join(state_match) %>%
  full_join(qte_full %>% select(chain_step, basin, stage, state, qte, qs))

# calculate monthly survival
surv_month <- qt_full %>%
  mutate(p = exp(-(1 / 12) * qs * qt * exp(qte * 0.5)),
         p_se = exp(-(1 / 12) * qs * qt * exp(qte * b1)),
         month = lubridate::month(date),
         season = factor(ifelse(month < 8, 0, 1),
                         levels = c(0, 1),
                         labels = c("summer", "winter"))) %>%
  select(chain_step, date, season, basin, stage, state, p, p_se) %>%
  gather(var, val, p, p_se) %>%
  group_by(var, basin, stage, state, season, date) %>%
  summarize(lo = quantile(val, prob = 0.16),
            mi = quantile(val, prob = 0.5),
            hi = quantile(val, prob = 0.84))




#==========
#========== Plot
#==========

# plot
p1 <- ggplot(data = surv_month %>% 
               filter(var == "p"),
             aes(x = date,
                 y = mi,
                 color = basin))+
  facet_rep_wrap(~stage,
                 nrow = 2)+
  geom_hline(yintercept = 0.5,
             size = 0.2,
             color = "gray50",
             linetype = 2)+
  geom_line(data = surv_month %>% 
              filter(var == "p_se"),
            size = 0.2,
            alpha = 0.5)+
  geom_line(size = 0.4)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Survival~probability~(monthly), 
                     breaks = c(0.1, 0.5, 0.9), 
                     limits = c(0, 1))+
  scale_x_date(name = "Year",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.15, 0.675),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')

# examine plot
p1

# export
# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fig_surv.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# p1
# dev.off()
  
  
  
  

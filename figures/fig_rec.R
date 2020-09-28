#==========
#========== Preliminaries
#==========

source("analysis/demographic_model/model/figures/fig_setup.R")





#==========
#========== Calculate annual recruitment
#==========

# extract variable names
rs_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rs"))}$var
rte_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rte"))}$var
rt_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rt"),
           !str_detect(fit_summary$var, "rte"))}$var

# extract scaling parameter
rs_full <- rstan::extract(fit, pars = rs_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(rs_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  gather(var, rs, -chain, -step, -chain_step) 

# extract season effects
rte_full <- rstan::extract(fit, pars = rte_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(rte_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  gather(var, rte, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         basin = levels(state_match$basin)[
           str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))]) %>%
  full_join(rs_full %>%
              select(chain_step, rs))

# extract recruitment
rt_full <- rstan::extract(fit, pars = rt_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(rt_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  gather(var, rt, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         basin = levels(state_match$basin)[
           str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))],
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b1 = data_list$b[time,1],
         b2 = data_list$b[time,2]) %>%
  full_join(rte_full %>% select(chain_step, basin, rte, rs))

# calculate annual recruitment
rec_annual <- rt_full %>%
  mutate(rec = rs * rt * b2* exp(rte * b1),
         year = lubridate::year(date)) %>%
  group_by(chain_step, basin, year) %>%
  summarize(rec = sum(rec)) %>%
  group_by(basin, year) %>%
  summarize(lo = quantile(rec, prob = 0.16),
            mi = quantile(rec, prob = 0.5),
            hi = quantile(rec, prob = 0.84))




#==========
#========== Plot
#==========

p1 <- ggplot(data = rec_annual,
             aes(x = year, 
                 y = mi, 
                 color = basin))+
  geom_line(size = 0.4)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Fecundity~(annual),
                     breaks = c(2, 8, 14),
                     limits = c(0, 16))+
  scale_x_continuous(name = "Year",
               limits = year_limits,
               breaks = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.2, 0.85),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')
  
p1

# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fit_fec.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p1
# dev.off()
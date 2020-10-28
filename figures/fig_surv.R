#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- annual_output[[1]]$setup$years
ids <- annual_output[[1]]$setup$ids





#==========
#========== Calculate monthly survival
#==========

# extract
surv_full <- parallel::mclapply(ids, function(i_){
  annual_proj_ = annual_output[[i_]]$annual_proj
  QQ_ = annual_proj_$QQ
  QQs_ = annual_proj_$QQs
  x_ = annual_proj_$x
  lapply(1:dim(QQ_)[3], function(t_){
    tibble(id = i_,
           surv = diag(QQ_[,,t_]),
           surv_s = diag(QQs_[,,t_]),
           x = x_[, t_],
           date = date_match[t_],
           st = 1:4
    ) 
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  full_join(state_match)
  ungroup()

# summarize
surv_sum <- surv_full %>%
  gather(var, val, surv, surv_s) %>%
  group_by(date, basin, stage, var) %>%
  summarize(lo = quantile(val, probs = c(0.16)),
            mi = median(val),
            hi = quantile(val, probs = c(0.84))) %>%
  ungroup()





#==========
#========== Plot dispersal
#==========

# plot labels
labs <- surv_sum %>%
  tidyr::expand(stage) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 1.075)

# plot
p1 <- surv_sum %>%
  filter(var %in% c("surv", "surv_s")) %>%
  select(date, basin, stage, var, mi) %>%
  spread(var, mi) %>%
  ggplot(aes(x = date, 
             y = surv_s, 
             color = basin))+
  facet_rep_wrap(~stage,
             nrow = 2)+
  geom_text(data = labs,
            aes(label = stage, 
                x = x, 
                y = y),
            color = "black", 
            size = 3,
            inherit.aes = F)+
  geom_hline(yintercept = 0.5,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_line(aes(y = surv),
            size = 0.2,
            alpha = 0.5)+
  geom_line(size = 0.4)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Survival~probability, 
                     breaks = c(0, 0.5, 1), 
                     labels = c("0","0.5","1"),
                     limits = c(0, 1.2))+
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
                    bottom='both')

# examine plot
p1

# export
# cairo_pdf(file = "figures/figs/fig_surv.pdf",
#           width = 3.5, height = 4, family = "Arial")
# p1
# dev.off()






#==========
#========== Calculate development probability
#==========

grow <- parallel::mclapply(ids, function(i_){
  annual_proj_ = annual_output[[i_]]$annual_proj
  GG_ = annual_proj_$GGs
  lapply(1:dim(GG_)[3], function(t_){
    tibble(id = i_,
           g = 1 - diag(GG_[,,t_]),
           t = t_
    ) %>%
      filter(g > 0) 
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  ungroup()


grow %>%
  group_by(t) %>%
  summarize(lo = quantile(g, probs = c(0.16), na.rm = T),
            mi = median(g, na.rm = T),
            hi = quantile(g, probs = c(0.84), na.rm = T)) %>%
  select(-t) %>%
  unique()

#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- proj_output[[1]]$setup$years
ids <- proj_output[[1]]$setup$ids





#==========
#========== Calculate monthly survival
#==========

# extract
surv_full <- parallel::mclapply(ids, function(i_){
  proj_ = proj_output[[i_]]$proj
  QQ_ = proj_$QQ
  QQs_ = proj_$QQs
  x_ = proj_$x
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
  full_join(state_match) %>%
  ungroup()

# summarize
surv_sum <- surv_full %>%
  gather(var, val, surv, surv_s) %>%
  group_by(date, basin, stage, var) %>%
  summarize(lo = quantile(val, probs = c(0.16)),
            mi = median(val),
            hi = quantile(val, probs = c(0.84))) %>%
  ungroup()

# juvenile covariance matrix
juv_covmat <- surv_sum %>%
  ungroup() %>%
  filter(stage == "juvenile",
         var == "surv_s") %>%
  select(date, basin,  mi) %>%
  spread(basin, mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov() %>%
  round(2)

# adult coveriance matrix
adult_covmat <- surv_sum %>%
  ungroup() %>%
  filter(stage == "adult",
         var == "surv_s") %>%
  select(date, basin,  mi) %>%
  spread(basin, mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov() %>%
  round(3)

# combine
covmat <- tibble(stage = factor(c("juvenile","adult")),
                 cov = c(format(juv_covmat[2,1], digits = 2, nsmall = 3), 
                         format(adult_covmat[2,1], digits = 2, nsmall = 2)),
                 var = c(format(mean(diag(juv_covmat)), digits = 2, nsmall = 2),
                         format(mean(diag(adult_covmat)), digits = 2, nsmall = 2)))





#==========
#========== Plot survival
#==========

# plot labels
labs <- surv_sum %>%
  tidyr::expand(stage) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 1.08)

covmat <- covmat %>%
  mutate(x = lubridate::as_date("2017-06-01"),
         y_cov = 1.03,
         y_var = 0.92,
         lab_cov = paste0("cov*~`=`*~`", cov,"`"),
         lab_var = paste0("~bar(var)*~`=`*~`", var,"`"))

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
            size = 3.5,
            inherit.aes = F)+
  geom_text(data = covmat,
            aes(x = x, 
                y = y_cov,
                label = lab_cov),
            color = "black", 
            size = 2.8,
            inherit.aes = F,
            parse = T)+
  geom_text(data = covmat,
            aes(x = x, 
                y = y_var,
                label = lab_var),
            color = "black", 
            size = 2.8,
            inherit.aes = F,
            parse = T)+
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
  proj_ = proj_output[[i_]]$proj
  GG_ = proj_$GGs
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

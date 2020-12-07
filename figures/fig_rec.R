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
#========== Extract and summarize dispersal
#==========

# extract
rec_full <- parallel::mclapply(ids, function(i_){
  proj_ = proj_output[[i_]]$proj
  RR_ = proj_$RR
  x_ = proj_$x
  lapply(1:dim(RR_)[3], function(t_){
    tibble(id = i_,
           rec = c(RR_[c(1,3),c(2,4),t_])[c(1,4)],
           x = x_[c(2,4), t_],
           date = date_match[t_]
    ) %>%
      filter(rec > 0) %>%
      group_by(date) %>%
      mutate(basin = 1:2,
             basin = factor(basin,
                            levels = c(1,2),
                            labels = c("south","north")),
             rec_tot = rec * x)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  ungroup()

# summarize
rec_sum <- rec_full %>%
  gather(var, val, rec, rec_tot, x) %>%
  group_by(date, basin,var) %>%
  summarize(lo = quantile(val, probs = c(0.16)),
            mi = median(val),
            hi = quantile(val, probs = c(0.84)))

# covariance matirx
rec_covmat <- rec_sum %>%
  ungroup() %>%
  filter(var == "rec") %>%
  select(date, basin, mi) %>%
  spread(basin, mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()

tot_covmat <- rec_sum %>%
  ungroup() %>%
  filter(var == "rec_tot") %>%
  select(date, basin, mi) %>%
  spread(basin, mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()

covmat <- tibble(type = c("rec","tot"),
                 cov = c(format(rec_covmat[2,1], digits = 2, nsmall = 2), 
                         format(tot_covmat[2,1], digits = 2, nsmall = 2)),
                 var = c(format(mean(diag(rec_covmat)), digits = 2, nsmall = 2),
                         format(mean(diag(tot_covmat)), digits = 2, nsmall = 2)))



#==========
#========== Plot per capita recruitent
#==========

# label
covmat <- covmat %>%
  filter(type == "rec") %>%
  mutate(x = lubridate::as_date("2018-01-01"),
         y_cov = 6.2,
         y_var = 5.6,
         lab_cov = paste0("cov*~`=`*~`", cov,"`"),
         lab_var = paste0("~bar(var)*~`=`*~`", var,"`"))

# plot
p1 <- ggplot(data = rec_sum %>%
               filter(var == "rec"),
             aes(x = date, 
                 y = mi, 
                 color = basin))+
  geom_line(size = 0.4)+
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
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Recruitment~(capita^{-1}),
                     breaks = c(0, 2, 4, 6),
                     limits = c(0, 6.5))+
  scale_x_date(name = "",
               limits = date_limits,
               breaks = date_breaks,
               labels = NULL)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(l = 0, r = 0, t = 0, b = -10),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.22, 0.92),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.title.y = element_text(margin = margin(t = 0, r = 16, b = 0, l = 0)))+
  coord_capped_cart(left = "both", 
                    bottom='both')

p1





#==========
#========== Plot total recruitent
#==========

# label
covmat <- covmat %>%
  filter(type == "tot") %>%
  mutate(x = lubridate::as_date("2018-01-01"),
         y_cov = 10,
         y_var = 8.8,
         lab_cov = paste0("cov*~`=`*~`", cov,"`"),
         lab_var = paste0("~bar(var)*~`=`*~`", var,"`"))

# plot
p2 <- ggplot(data = rec_sum %>%
               filter(var == "rec_tot"),
             aes(x = date, 
                 y = mi, 
                 color = basin))+
  geom_line(size = 0.4)+
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
  scale_color_manual(name = "",
                     values = basin_colors,
                     guide = F)+
  scale_y_continuous(name = Recruitment~(total),
                     breaks = c(0, 2, 4, 6, 8, 10),
                     limits = c(0, 10.5))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')

p2





#==========
#========== Combine plots
#==========

# combine
p3 <- plot_grid(p1, p2,
                ncol = 1,
                align = "v",
                axis = "tblr",
                label_size = 12)

# export
p3

# export
# cairo_pdf(file = "figures/figs/fig_rec.pdf",
#           width = 3.5, height = 4, family = "Arial")
# p3
# dev.off()
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
#========== Extract and summarize dispersal
#==========

# extract
rec_full <- parallel::mclapply(ids, function(i_){
  annual_proj_ = annual_output[[i_]]$annual_proj
  RR_ = annual_proj_$RR
  x_ = annual_proj_$x
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








#==========
#========== Plot per capita recruitent
#==========

# plot
p1 <- ggplot(data = rec_sum %>%
               filter(var == "rec"),
             aes(x = date, 
                 y = mi, 
                 color = basin))+
  geom_line(size = 0.4)+
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

# plot
p2 <- ggplot(data = rec_sum %>%
               filter(var == "rec_tot"),
             aes(x = date, 
                 y = mi, 
                 color = basin))+
  geom_line(size = 0.4)+
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
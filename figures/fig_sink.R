#==========
#========== Preliminaries
#==========

# source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- proj_output[[1]]$setup$years
ids <- proj_output[[1]]$setup$ids





#==========
#========== Extract and summarize dynamics
#==========
sink_full <- parallel::mclapply(ids, function(i_){
  proj_output[[i_]]$proj$X_noD_proj %>% 
    t() %>% 
    as_tibble(.name_repair = "universal") %>%
    set_names(c("1","2","3","4")) %>%
    mutate(year = 1990 + row_number()) %>%
    gather(state, val, -year) %>%
    mutate(stage = ifelse(state %in% c(1, 3), "juvenile", "adult") %>%
             factor(levels = c("juvenile", "adult")),
           basin = ifelse(state < 3, "south", "north") %>%
             factor(levels = c("south", "north")),
           id = i_)
}) %>%
  bind_rows() 

sink_sum <- sink_full %>%
  mutate(val = ifelse(val < 0.01, 0.01, val)) %>%
  group_by(year, stage, basin) %>%
  summarize(lo = quantile(val, probs = c(0.16)),
            mi = median(val),
            hi = quantile(val, probs = c(0.84))) %>%
  ungroup()




#==========
#========== Plot 
#==========

# plot labels
labs <- sink_sum %>%
  tidyr::expand(stage) %>%
  mutate(x = 2005,
         y = 2.5e4)

# plot 
p1 <-  ggplot(data = sink_sum,
              aes(x = year, 
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
  geom_line(size = 0.35)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_fill_manual(name = "",
                    values = basin_colors,
                    guide = F)+
  geom_ribbon(aes(ymin = lo, 
                  ymax = hi),
              alpha = 0.2,
              linetype = 0)+
  scale_y_continuous(name = Projected~abundance~(no~dispersal),
                     trans = "log",
                     limits = c(0.01, 10e4),
                     breaks = c(0.01, 10, 10000),
                     labels = c("0", "10", expression(10^3))
                     )+
  scale_x_continuous(name = "Year",
                     limits = year_limits,
                     breaks = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(-1, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.15, 0.4),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        # axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        plot.margin = margin(t = 1, r = 3, b = 1, l = 1))+
  coord_capped_cart(left = "both", 
                    bottom='both')+
  guides(fill = guide_legend(override.aes = list(fill = NA)),
         color = guide_legend(override.aes = list(size = 0.5,
                                                  shape = NA)))

# examine plot
p1

# export
# cairo_pdf(file = "figures/figs/fig_sink.pdf",
#           width = 3.5, height = 4, family = "Arial")
# p1
# dev.off()

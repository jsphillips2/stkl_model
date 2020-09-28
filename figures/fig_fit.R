#==========
#========== Preliminaries
#==========

source("analysis/demographic_model/model/figures/fig_setup.R")




#==========
#========== Plot set up
#==========

# extract scaling parameter
y_scale <- out_in$y_scale

# extract fit
x_clean <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "x"), !str_detect(fit_summary$var, "x0")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) %>%
  full_join(state_match) %>%
  select(basin, state, stage, date, name, lo, mi, hi)





#==========
#========== Plot 
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
                     breaks = c(2,6,10))+
  scale_x_date(name = "Year",
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
# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fig_fit.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# p1
# dev.off()
  
  
  
  
  
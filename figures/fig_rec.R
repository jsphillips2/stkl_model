#==========
#========== Preliminaries
#==========

source("analysis/demographic_model/model/figures/fig_setup.R")





#==========
#========== Calculate annual recruitment
#==========

# extract variable names
rt_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rt"))}$var

# extract recruitment
rt_full <- fit_summary %>%
  rename(lo = `16%`,
         mi = `50%`,
         hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "rt")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         basin = levels(state_match$basin)[
           str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))],
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b = data_list$b[time])
 


#==========
#========== Plot
#==========

p1 <- ggplot(data = rt_full,
             aes(x = date, 
                 y = mi, 
                 color = basin))+
  geom_line(size = 0.4)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Recruitment~(adult^{-1}),
                     breaks = c(0, 2, 4, 6),
                     limits = c(-0.1, 6.1))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.2, 1),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')
  
p1

# cairo_pdf(file = "figures/figs/fig_rec.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p1
# dev.off()
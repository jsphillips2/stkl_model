#==========
#========== Preliminaries
#==========

source("analysis/demographic_model/model/figures/fig_setup.R")






#==========
#========== Prepare data
#==========

# filter out movement
move_annual <- aa_cont_sum %>%
  mutate(year = unique(lubridate::year(date_match))[time]) %>%
  full_join(state_match %>%
              rename(recipient_state = state,
                     recipient_basin = basin,
                     recipient_stage = stage,
                     row = st)) %>%
  full_join(state_match %>%
              rename(donor_state = state,
                     donor_basin = basin,
                     donor_stage = stage,
                     col = st)) %>%
  filter(recipient_stage == donor_stage,
         recipient_basin != donor_basin)

# calculate net movement
net_move <- move_annual %>%
  group_by(year, donor_basin, recipient_basin, recipient_stage) %>%
  summarize(move = sum(value)) %>%
  ungroup() %>%
  select(year, recipient_basin, recipient_stage, move) %>%
  spread(recipient_basin, move) %>%
  mutate(net_move =  north - south) 

# linearly interpolate by basin (for changing line colors)
interp_juv <- approxfun(y = {net_move %>% filter(recipient_stage == "juvenile")}$net_move,
                 x = {net_move %>% filter(recipient_stage == "juvenile")}$year)
interp_adult <- approxfun(y = {net_move %>% filter(recipient_stage == "adult")}$net_move,
                 x = {net_move %>% filter(recipient_stage == "adult")}$year)

# combine linear interpolations
net_move_clean <- net_move  %>%
  tidyr::expand(recipient_stage, 
                year = seq(min(year),
                           max(year), 
                           by = 0.001)) %>%
  mutate(net_move = ifelse(recipient_stage == "juvenile", 
                                   interp_juv(year), 
                                   interp_adult(year)),
         recipient_basin = factor(ifelse(net_move < 0, 
                                         "southward",
                                         "northward"))) 




#==========
#========== Plot
#==========

p1 <- ggplot(data = net_move_clean,
             aes(x = year, 
                 y = net_move, 
                 color = recipient_basin))+
  facet_rep_wrap(~recipient_stage, nrow = 2)+
  geom_line()+
  geom_ribbon(data = net_move_clean %>% 
                filter(recipient_basin == "northward"),
              aes(ymin = 0, 
                  ymax = net_move), 
              color = "gray20",
              fill = "white",
              outline.type = "upper",
              size = 0.4)+
  geom_ribbon(data = net_move_clean %>% 
                filter(recipient_basin == "southward"),
              aes(ymin = net_move, 
                  ymax = 0), 
              color = "dodgerblue",
              fill = "white",
              outline.type = "lower",
              size = 0.4)+
  geom_ribbon(data = net_move_clean %>% 
                filter(recipient_basin == "southward"),
              aes(ymin = -0.06, 
                  ymax = 0.06), 
              color = "white",
              fill = "white",
              linetype = 0)+
  geom_hline(yintercept = 0,size = 0.2, color = "gray50", linetype = 2)+
  scale_color_manual(name = "",
                     values = basin_colors %>% 
                       set_names(c("southward","northward")))+
  scale_y_continuous("Net movment (annual)",
                     limits = c(-4.7, 2.7),
                     breaks = c(-4, -2, 0, 2))+
  scale_x_continuous(name = "Year",
                     limits = year_limits,
                     breaks = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.2, 0.675),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')

# examine plot
p1
  
# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fit_move.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# p1
# dev.off()


  








  
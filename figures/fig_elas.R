#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- annual_output[[1]]$setup$years
theta_names <- annual_output[[1]]$setup$theta_names
ids <- annual_output[[1]]$setup$ids





#==========
#========== Extract and summarize elasticities
#==========

# extract
elast <- parallel::mclapply(ids, function(i_){
  sens_ = annual_output[[i_]]$sens
  lapply(1:length(sens_), function(y_){
    tibble(id = i_,
           year = unique(years)[y_],
           elas = sens_[[y_]]$l_elas_p
    ) %>%
      mutate(theta = theta_names[row_number()])
  }) %>%
    bind_rows()
}) %>%
  bind_rows() 


# summarize
elast_sum <- elast %>%
  mutate(name = str_split(theta, "") %>% map_chr(~as.character(.x[[1]])),
         index = str_split(theta, "") %>% map_chr(~as.character(.x[[2]])),
         name = paste0(name, index),
         season = str_split(theta, "") %>% map_chr(~as.character(.x[[3]]))) %>%
  group_by(id, year, name, index) %>%
  summarize(elas = sum(elas)) %>%
  group_by(year, name, index) %>%
  summarize(lo = quantile(elas, probs = c(0.16), na.rm = T),
            mi = median(elas, na.rm = T),
            hi = quantile(elas, probs = c(0.84), na.rm = T)) %>%
  ungroup()


# mean across time
elast_sum_mean <- elast_sum %>%
  group_by(name) %>%
  summarize(mi = mean(mi)) %>%
  arrange(-mi) %>%
  ungroup()

# set parameter order
theta_order <- elast_sum_mean$name 
theta_labs <- c(expression(phi["j,n"]),
                expression(rho[n]),
                expression(phi["a,n"]),
                expression(gamma),
                expression(phi["a,s"]),
                expression(delta[s%->%n]),
                expression(rho[s]),
                expression(phi["j,s"]),
                expression(delta[n%->%s]))



#==========
#========== Plot
#==========

p1 <- elast_sum %>%
  ungroup() %>%
  mutate(name = factor(name, levels = theta_order),
         pos = as.numeric(name) - 
           0.25 * (year - mean(year)) / max((year - mean(year)))) %>%
  ggplot(aes(pos, mi))+
  geom_hline(yintercept = 0,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_errorbar(aes(ymin = lo,
                    ymax = hi,
                    color = year),
                width = 0,
                size = 0.2)+
  geom_point(aes(fill = year,
                 color = year),
             size = 1,
             shape = 21, 
             stroke = 0.2)+
  geom_segment(data = elast_sum_mean %>%
                 mutate(name = factor(name, levels = theta_order),
                        pos = as.numeric(name),
                        xmin = pos - 0.3,
                        xmax = pos + 0.3),
               aes(x = xmin,
                   y = mi,
                   xend = xmax,
                   yend = mi),
               size = 0.5)+
  scale_x_reverse(Demographic~rate,
                     breaks = 1:9,
                     labels = theta_labs)+
  scale_y_continuous(Proportional~sensitivity~of~lambda,
                     breaks = c(-1.5, 0, 1.5),
                     labels = c("-1.5","0","1.5"),
                     limits = c(-1.61, 1.61))+
  scale_fill_gradient2("",
                       low ="dodgerblue", 
                       mid="gray70", 
                       high="firebrick", 
                       midpoint = 2005,
                       breaks = c(1995,2015),
                       guide = guide_colourbar(ticks = F))+
  scale_color_gradient2(low ="dodgerblue",
                        mid="gray70",
                        high="firebrick",
                        midpoint = 2005,
                        guide = F)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_flip(left = "both", 
                    bottom='both')

# examine plot
p1

# export
# cairo_pdf(file = "figures/figs/fig_elas.pdf",
#           width = 3.5, height = 4.5, family = "Arial")
# p1
# dev.off()



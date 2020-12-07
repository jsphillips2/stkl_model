#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- proj_output[[1]]$setup$years
ids <- proj_output[[1]]$setup$ids
year_limits_spec <- c(1990, 2020) 
year_breaks_spec <- c(1995, 2005, 2015) 




#==========
#========== Extract and summarize elasticities
#==========

# extract
lam_full <- parallel::mclapply(ids, function(i_){
  sens_ = proj_output[[i_]]$sens
  lapply(1:length(sens_), function(y_){
    tibble(id = i_,
           year = unique(years)[y_],
           l = sens_[[y_]]$l,
           l_asym = sens_[[y_]]$l_asym,
           r = sens_[[y_]]$r,
           r_asym = sens_[[y_]]$r_asym 
    )
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  gather(var, val, l, l_asym, r, r_asym) 

# summarize
lam <- lam_full %>%
  group_by(year, var) %>%
  summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
            mi = median(val, na.rm = T),
            hi = quantile(val, probs = c(0.84), na.rm = T)) %>%
  ungroup()

# summarize median lambda
lam_full %>%
  filter(var %in% c("l","l_asym")) %>%
  group_by(id, var) %>%
  summarize(geom = prod(val)^(1/length(val)),
            arith = mean(val)) %>%
  ungroup() %>%
  gather(type, val, geom, arith) %>%
  group_by(var, type) %>%
  summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
            mi = median(val, na.rm = T),
            hi = quantile(val, probs = c(0.84), na.rm = T))




#==========
#========== Wavelet analysis
#==========

# extract relevant data
wave_prep <- lam %>%
  filter(var %in% c("r", "r_asym")) %>%
  select(year, var, mi) %>%
  spread(var, mi)

# wavelet transform
wave_r <- analyze.wavelet(wave_prep,
                          "r",
                          loess.span = 0,
                          lowerPeriod = 2,
                          upperPeriod = 29,
                          dj = 1 / 20,
                          dt = 1)
wave_asym <- analyze.wavelet(wave_prep,
                             "r_asym",
                             loess.span = 0,
                             lowerPeriod = 2,
                             upperPeriod = 29,
                             dj = 1 / 20,
                             dt = 1)

# extract values from wavelet transform
wave_d <- as_tibble(wave_r$Power) %>%
  mutate(period = wave_r$Period) %>%
  gather(id, power, -period) %>%
  mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
         year = 1990 + id,
         type = "Transient") %>%
  bind_rows(as_tibble(wave_asym$Power) %>%
              mutate(period = wave_asym$Period) %>%
              gather(id, power, -period) %>%
              mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
                     year = 1990 + id,
                     type = "Asymptotic")) %>%
  mutate(type = factor(type, levels = c("Transient","Asymptotic")))

# define power quantiles
qt <- quantile(wave_d$power, prob = seq(0,1, length.out = 100))

# define power levles
wave_d$power_quant <- qt[cut(wave_d$power, breaks = qt)]
wave_d$power_quant[is.na(wave_d$power_quant)] <- 0





#==========
#========== Plot lambda
#==========

p1 <- lam %>%
  filter(var %in% c("l","l_asym")) %>%
  mutate(type = factor(var,
                       levels = c("l","l_asym"),
                       labels = c("Transient"," Asymptotic"))) %>%
  ggplot(aes(year, mi))+
  facet_wrap(~type)+
  geom_hline(yintercept = 1, 
             size = 0.2, 
             color = "black", 
             linetype = 2)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2,
              linetype = 0)+
  geom_line()+
  scale_y_continuous(Population~growth~rate~(lambda),
                     trans = "log",
                     breaks = c(1/9, 1/3, 1, 3),
                     labels = c("1/9", "1/3", "1", "3"))+
  scale_x_continuous(name = "",
                     limits = year_limits_spec,
                     breaks = year_breaks_spec,
                     labels = NULL)+
  theme(panel.border = element_blank(),
        panel.spacing.x = unit(-1, "lines"),
        plot.margin = margin(l = 0, r = 0, t = 0, b = -10),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", bottom='both')

# examine plot
p1

  



#==========
#========== Plot wavelet
#==========

p2 <- ggplot(data = wave_d,
       aes(x = year, 
           y= period))+
  facet_wrap(~type)+
  geom_tile(aes(fill = power_quant)) +
  geom_contour(aes(z = power_quant),
               bins = 4,
               color = "black",
               size = 0.3)+
  scale_y_continuous("Period (years)",
                     trans = "log",
                     breaks = 2 * c(1:5),
                     limits = c(1.9, 10.9))+
  scale_x_continuous(name = "Year",
                     limits = year_limits_spec,
                     breaks = year_breaks_spec)+
  scale_fill_gradient2(name = "",
                       low ="dodgerblue", 
                       mid="gray90", 
                       high="firebrick", 
                       midpoint = 0.1,
                       limits = c(0, 1.3),
                       breaks = c(0.1, 1.2),
                       labels = c("Low power",
                                  "High power"),
                       guide = F)+
  theme(strip.text = element_blank(),
        legend.key.width = unit(0.9, "lines"),
        legend.key.height= unit(0.5, "lines"),
        legend.direction = "vertical",
        legend.spacing.x = unit(0.2, "lines"), 
        panel.border = element_blank(),
        panel.spacing.x = unit(-1, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  coord_capped_cart(left = "both", bottom='both')

# examine plot
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
# cairo_pdf(file = "figures/figs/fig_lam.pdf",
#           width = 5, height = 4, family = "Arial")
# p3
# dev.off()













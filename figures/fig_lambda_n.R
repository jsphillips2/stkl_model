#==========
#========== Preliminaries
#==========

source("analysis/demographic_model/model/figures/fig_setup.R")




#==========
#========== Prepare data
#==========

# extract lambda
lambda <- sens_sum %>%
  filter(var == "l") %>%
  select(var, year, value) %>%
  unique() 

# extract r
r_d <-sens_sum %>%
  filter(var == "r") %>%
  select(var, year, value) %>%
  unique() 





#==========
#========== Wavelet analysis
#==========

# wavelet transform
wave <- analyze.wavelet(r_d,
                        "value",
                        loess.span = 0,
                        lowerPeriod = 2,
                        upperPeriod = 29,
                        dj = 1 / 20,
                        dt = 1,
                        make.pval = T,
                        n.sim = 4000,
                        method = "white.noise")

# extract values from wavelet transform
wave_d <- as_tibble(wave$Power) %>%
  mutate(period = wave$Period) %>%
  gather(id, power, -period) %>%
  mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
         year = 1990 + id) %>%
  full_join(as_tibble(wave$Power.pval) %>%
              mutate(period = wave$Period) %>%
              gather(id, pval, -period) %>%
              mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
                     year = 1990 + id,
                     bin = ifelse(pval < 0.05, 1, 0)))

# define power quantiles
qt <- quantile(wave_d$power, prob = seq(0,1, length.out = 100))

# define power levles
wave_d$power_quant <- qt[cut(wave_d$power, breaks = qt)]
wave_d$power_quant[is.na(wave_d$power_quant)] <- 0





#==========
#========== Density dependence
#==========

# calculate total populatino size
n_real <- ((c(0, {sens_sum %>%
    filter(var == "r") %>%
    select(time, value) %>%
    unique()}$value) %>% cumsum()) + log(sum({n_sum %>% filter(time == 1)}$value))) %>%
  exp()

# combine with growth rate
n_r_clean <- sens_sum %>%
  filter(var == "r") %>%
  select(year, value) %>%
  unique() %>%
  mutate(n = n_real[1:29])

# fit model
# m <- rstanarm::stan_glm(value ~ n,
#                         data = n_r_clean,
#                         family = gaussian(link = "identity"),
#                         seed = 2e3)



#==========
#========== Prepare panel A: Lambda
#==========

p1 <- ggplot(data = lambda,
             aes(x = year, y =value))+
  geom_hline(yintercept = 1,
             size = 0.2, 
             color = "gray50", 
             linetype = 2)+
  geom_line(size = 0.4)+
  scale_y_continuous(name = "Growth rate ("~lambda~")",
                     limits = c(1 / 3, 3),
                     breaks = c(0.5, 1, 2),
                     trans = "log")+
  scale_x_continuous(name = "",
                     limits = year_limits,
                     breaks = year_breaks,
                     labels = NULL)+
  theme(legend.position = c(0.5,0.95),
        legend.title = element_blank(),
        legend.key.width = unit(1.15, "lines"),
        legend.key.height= unit(0.9, "lines"),
        legend.text = element_text(size = 8),
        plot.margin = margin(0,0,-10,0),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.line.y = element_line(size = 0.25),
        axis.line.x = element_line(size = 0.25),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))
        )+
  coord_capped_cart(left = "both",
                    bottom='both') 
p1




#==========
#========== Prepare panel B: Wavelet
#==========



p2_base <- ggplot(data = wave_d,
                  aes(x = year, 
                      y= period))+
  geom_tile(aes(fill = power_quant)) +
  geom_contour(aes(z = bin),
               bins = 2,
               color = "black",
               size = 0.4)+
  scale_y_continuous("Period (years)",
                     trans = "log",
                     breaks = c(3, 6, 12, 24))+
  scale_x_continuous(name = "",
                     limits = year_limits,
                     breaks = year_breaks)+
  scale_fill_gradient2(name = "",
                       low ="dodgerblue", 
                       mid="gray90", 
                       high="firebrick", 
                       midpoint = 0.1,
                       limits = c(0, 1.1),
                       breaks = c(0.095, 1.05),
                       labels = c("Low power",
                                  "High power"),
                       guide = guide_colourbar(ticks = F))+
  theme(legend.key.width = unit(0.9, "lines"),
        legend.key.height= unit(0.5, "lines"),
        legend.direction = "vertical",
        legend.spacing.x = unit(0.2, "lines"), 
        plot.margin = margin(-10,0,0,0),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  coord_capped_cart(left = "both", bottom='both')
p2 <- p2_base + theme(legend.position = "none")
p2
p2_leg <- get_legend(p2_base)



#==========
#========== Prepare panel C: Deviation
#==========


p3 <- ggplot(data = n_r_clean,
             aes(x = n, 
                 y = value))+
  geom_hline(yintercept = 0, 
             size = 0.2, 
             color = "gray50", 
             linetype = 2)+
  geom_vline(xintercept = 2285.691,
             size = 0.2, 
             color = "gray50", 
             linetype = 2)+
  geom_smooth(method = "lm", 
              se = F, 
              color = "black", 
              size = 0.4)+
  geom_point(size = 0.8,
             shape = 1)+
  scale_y_continuous(name = "",
                     limits = log(c(1 / 3, 3)),
                     breaks = NULL)+
  scale_x_continuous(name = "Total abundance",
                     limits = c(0, 31),
                     breaks = c(4, 16, 28))+
  theme(plot.margin = margin(0,0,0,-10),
        panel.border = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)))+
  coord_capped_cart(bottom = "both")

p3




#==========
#========== Combine
#==========

p4 <- plot_grid(p1, p2,
                ncol = 1,
                align = "v",
                axis = "tblr",
                label_size = 12)
p4

p5 <- plot_grid(p2_leg, NULL,
                ncol = 2,
                rel_widths = c(0.5, 0.5),
                label_size = 12)
p5

p6 <- plot_grid(NULL,p3, p5, NULL,
                ncol = 1,
                rel_heights = c(0.025, 1, 0.5, 0.305),
                label_size = 12)
p6


p7 <- plot_grid(p4, p6,
                ncol = 2,
                rel_heights = c(1, 1),
                align = "h",
                axis = "tblr",
                label_size = 12)
p7

# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fig_lam.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p7
# dev.off()

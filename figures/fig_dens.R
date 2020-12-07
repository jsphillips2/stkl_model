#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")
library(nlme)

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- proj_output[[1]]$setup$years
ids <- proj_output[[1]]$setup$ids
year_limits_spec <- c(1990, 2020) 
year_breaks_spec <- c(1995, 2005, 2015) 




#==========
#========== Extract and summarize 
#==========

# extract
lam_x_full <- parallel::mclapply(ids, function(i_){
  sens_ = proj_output[[i_]]$sens
  proj_ = proj_output[[i_]]$proj$X_proj
  out_ = lapply(1:length(sens_), function(y_){
    tibble(id = i_,
           year = unique(years)[y_],
           l_asym = sens_[[y_]]$l_asym,
           r_asym = sens_[[y_]]$r_asym ,
           x = sum(proj_[,y_])
    )
  }) %>%
    bind_rows()
  m_ <- gls(r_asym ~ x, data = out_, correlation = corCAR1(form = ~year))
  out_ = out_ %>%
    mutate(int = coef(m_)[1],
           slope = coef(m_)[2],
           phi_logit = as.numeric(m_$modelStruct$corStruct),
           phi = exp(phi_logit) / (1 + exp(phi_logit)))
  return(out_)
    
}) %>%
  bind_rows() %>%
  select(-phi_logit) %>%
  gather(var, val, l_asym, r_asym, x, int, slope, phi) 

# summarize
lam_x <- lam_x_full %>%
  group_by(year, var) %>%
  summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
            mi = median(val, na.rm = T),
            hi = quantile(val, probs = c(0.84), na.rm = T)) %>%
  ungroup()

# restructure for examining DD
lam_wide <- lam_x %>%
  filter(var == "x") %>%
  rename(lo_x = lo,
         mi_x = mi,
         hi_x = hi) %>%
  select(-var) %>%
  full_join(lam_x %>%
              filter(var %in% c("r_asym",
                                "l_asym")))



# predicted values
preds <- lam_x_full %>%
  filter(var %in% c("int","slope")) %>%
  select(-year) %>%
  unique() %>%
  spread(var, val) %>%
  expand(nesting(id, int, slope),
         mi_x = seq(min(lam_wide$mi_x),
                    max(lam_wide$mi_x),
                    length.out = 100)) %>%
  mutate(val = int + slope * mi_x) %>%
  group_by(mi_x) %>%
  summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
            mi = median(val, na.rm = T),
            hi = quantile(val, probs = c(0.84), na.rm = T)) %>%
  ungroup()


# density dependence summaries
dens_sum <- lam_x_full %>%
  filter(var %in% c("int","slope","phi")) %>%
  select(-year) %>%
  unique() %>%
  spread(var, val) %>%
  mutate(k = int / abs(slope),
         lam0 = exp(int)) %>%
  gather(var, val, -id) %>%
  group_by(var) %>%
  summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
            mi = median(val, na.rm = T),
            hi = quantile(val, probs = c(0.84), na.rm = T)) %>%
  ungroup()



#==========
#========== Extract and summarize 
#==========

p1 <- ggplot(data = lam_wide %>%
               filter(var == "l_asym"),
             aes(x = mi_x, 
                 y = mi, 
                 color = year))+
  geom_hline(yintercept = 1, 
             size = 0.2, 
             color = "black",
             linetype = 2)+
  geom_vline(xintercept = {dens_sum %>% filter(var == "k")}$mi, 
             size = 0.2, 
             color = "black",
             linetype = 2)+
  geom_errorbarh(aes(xmin = lo_x, 
                     xmax  = hi_x), 
                 size = 0.2)+
  geom_errorbar(aes(ymin = lo, 
                    ymax  = hi), 
                size = 0.2)+
  geom_point(size = 1)+
  geom_line(data = preds,
            aes(x = mi_x,
                y = exp(mi)),
            inherit.aes = F)+
  geom_ribbon(data = preds,
            aes(x = mi_x,
                ymin = exp(lo),
                ymax = exp(hi)),
            inherit.aes = F,
            linetype = 0,
            alpha = 0.2)+
  scale_y_continuous(Asymptotic~growth~rate~(lambda),
                     trans = "log",
                     breaks = c(1/4, 1, 4),
                     labels = c("1/4", "1", "4"),
                     limits = c(1/4.5, 4))+
  scale_x_continuous(name = Total~relative~abundance,
                     limits = c(0, 16),
                     breaks = c(0, 4, 8, 12, 16))+
  scale_color_gradient2("",
                        low ="dodgerblue", 
                        mid="gray70", 
                        high="firebrick", 
                        midpoint = 2005,
                        breaks = c(1995,2015),
                        guide = guide_colourbar(ticks = F))+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.75,0.9),
        legend.direction = "horizontal",
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')

# exmine plot
p1

# export
# cairo_pdf(file = "figures/figs/fig_dens.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p1
# dev.off()

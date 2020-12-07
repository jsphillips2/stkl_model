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
disp_full <- parallel::mclapply(ids, function(i_){
  proj_ = proj_output[[i_]]$proj
  DD_ = proj_$DD
  DDs_ = proj_$DDs
  x_ = proj_$x
  lapply(1:dim(DD_)[3], function(t_){
    tibble(id = i_,
           disp = 1 - diag(DD_[,,t_]),
           disp_s = 1 - diag(DDs_[,,t_]),
           x = x_[, t_],
           date = date_match[t_]
    ) %>%
      filter(disp > 0) %>%
      group_by(date) %>%
      mutate(basin = 1:2,
             basin = factor(basin,
                            levels = c(1,2),
                            labels = c("south","north")),
             move = disp * x,
             move_s = disp_s * x)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  ungroup()

# summarize
disp_sum <- disp_full %>%
  gather(var, val, disp, disp_s, move, move_s, x) %>%
  group_by(date, basin,var) %>%
  summarize(lo = quantile(val, probs = c(0.16)),
            mi = median(val),
            hi = quantile(val, probs = c(0.84)))

# calculate net and cumulative movement 
disp_net <- disp_full %>%
  select(id, date, basin, move) %>%
  spread(basin, move) %>%
  group_by(id) %>%
  arrange(id, date) %>%
  mutate(net = south - north,
         cumsum = cumsum(net)) %>%
  select(-south, -north) %>%
  gather(var, val, net, cumsum) %>%
  group_by(date, var) %>%
  summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
            mi = median(val, na.rm = T),
            hi = quantile(val, probs = c(0.84), na.rm = T)) %>%
  ungroup()

# covariance matirx
disp_covmat <- disp_sum %>%
  ungroup() %>%
  filter(var == "disp_s") %>%
  select(date, basin, mi) %>%
  spread(basin, mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()

covmat <- tibble(cov = format(disp_covmat[2,1], digits = 2, nsmall = 2),
                 var = format(mean(diag(disp_covmat)), digits = 2, nsmall = 2))







#==========
#========== Plot dispersal
#==========

# label
covmat <- covmat %>%
  mutate(x = lubridate::as_date("2018-01-01"),
         y_cov = 1.03,
         y_var = 0.93,
         lab_cov = paste0("cov*~`=`*~`", cov,"`"),
         lab_var = paste0("bar(var)*~`=`*~`", var,"`"))

# plot
p1 <- disp_sum %>%
  filter(var %in% c("disp", "disp_s")) %>%
  select(date, basin, var, mi) %>%
  spread(var, mi) %>%
  ggplot(aes(x = date, 
             y = disp_s, 
             color = basin))+
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
  geom_hline(yintercept = 0.5,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_line(aes(y = disp),
            size = 0.2,
            alpha = 0.5)+
  geom_line(size = 0.4)+
  scale_color_manual(name = "",
                     values = basin_colors)+
  scale_y_continuous(name = Dispersal~probability, 
                     breaks = c(0, 0.5, 1), 
                     labels = c("0","0.5","1"),
                     limits = c(0, 1.1))+
  scale_x_date(name = "",
               limits = date_limits,
               breaks = date_breaks,
               labels = NULL)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(l = 0, r = 0, t = 0, b = -10),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.15, 0.9),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom='both')

# examine plot
p1




#==========
#========== Plot net movement
#==========

# plot labels
labs <- tibble(x = lubridate::as_date("1993-07-01"),
               y = c(-2.9, 0.9),
               label = c("southward",
                          "northward"))
labs2 <- tibble(x = lubridate::as_date("2012-07-01"),
                y = -1.5,
                label = "cumulative")

labs3 <- tibble(x = lubridate::as_date("2012-07-01"),
                y = 0.3,
                label = "incremental")

# plot
p2 <- disp_net %>%
  ggplot(aes(x = date, 
             y = mi, 
             linetype = var))+
  geom_hline(yintercept = 0,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_line(size = 0.4)+
  geom_text(data = labs,
            aes(label = label, 
                x = x, 
                y = y),
            color = "black", 
            size = 2.8,
            inherit.aes = F)+
  geom_text(data = labs2,
            aes(label = label, 
                x = x, 
                y = y),
            color = "black", 
            size = 3,
            inherit.aes = F)+
  geom_text(data = labs3,
            aes(label = label, 
                x = x, 
                y = y),
            color = "black", 
            size = 3,
            inherit.aes = F)+
  scale_linetype_manual("",
                        values = c(3, 1), 
                        guide = F)+
  scale_y_continuous(name = Net~movement,
                     breaks = c(-3, -2, -1, 0, 1),
                     limits = c(-3, 1.2))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.15, 0.95),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.title.y = element_text(margin = margin(t = 0, r = 17, b = 0, l = 0)))+
  coord_capped_cart(left = "both", 
                    bottom='both')

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
# cairo_pdf(file = "figures/figs/fig_disp.pdf",
#           width = 3.5, height = 4, family = "Arial")
# p3
# dev.off()
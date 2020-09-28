#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(cowplot)


options(mc.cores = parallel::detectCores()-2)

# import data
name <- "full"
proj_sum <-read_rds(paste0("analysis/demographic_model/model/output/",name,"/proj_sum.rds"))
sens_sum <- proj_sum$sens_sum %>%
  mutate(year = (1991:2020)[time])
var_sum <- proj_sum$var_sum
pars_cont <- var_sum$pars_cont

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,-4),
                  strip.text = element_text(size = 8),
                  strip.text.y = element_text(angle = -90, margin=margin(0,0,0,2)),
                  strip.text.x = element_text(margin=margin(0,0,2,0)),
                  legend.text = element_text(size = 10),
                  axis.text = element_text(size = 10, color="black"),
                  axis.title.y = element_text(angle = 90, margin=margin(0,10,0,0)),
                  axis.title.x = element_text(margin = margin(10,0,0,0)),
                  panel.spacing = unit(0.1, "lines")))





#==========
#========== Prepare panel A: Sensitivity
#==========

sens_sum %>%
  filter(var == "l_sens")




# annotation
labs_a <- sens_sum %>%
  tidyr::expand(recipient, donor) %>%
  mutate(x = 2003.5,
         y = 1.25 * 0.98,
         labels = interaction(recipient, donor) %>%
           factor(levels = c("juvenile.juvenile",
                             "juvenile.adult",
                             "adult.juvenile",
                             "adult.adult"),
                  labels = c("juvenile to juvenile",
                             "adult to juvenile",
                             "juvenile to adult",
                             "adult to adult")))
# plot
p1a <- sens_sum %>%
  filter(var == "l_sens") %>%
  ggplot(aes(year, value))+
  facet_grid(recipient ~ donor)+
  geom_text(data = labs_a,
            aes(label = labels, x = x, y = y),
            color = "black", size = 3)+
  geom_line(aes(color = group, linetype = group), size = 0.4)+
  scale_color_manual("",
                     values = c("dodgerblue","dodgerblue","gray35","gray35"))+
  scale_linetype_manual("", values = c(1,2,1,2))+
  # scale_y_continuous("Sensitivity of "~lambda,
  #                    limits = c(-0.05, 1.25),
  #                    breaks = c(0.2, 0.6, 1.0))+
  # scale_x_continuous("",breaks = NULL,limits = c(1990, 2016))+
  theme(legend.key.width = unit(0.26,"in"),
        legend.key.height = unit(0.4,"in"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text.align = 0.35,
        legend.text = element_text(size = 8),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(5,1,7,1))+
  guides(linetype = guide_legend(override.aes = list(size = 0.6)))
p1_leg <- get_legend(p1a)
p1 <- p1a + theme(legend.position = "none")




#==========
#========== Prepare panel B: Contributions
#==========

# create matrix of parameter identies
id_mat <- crossing(row = 1:4, col = 1:4) %>%
  arrange(col, row) %>%
  mutate(id = row_number(),
         donor_basin = ifelse(col %in% c(1,2), "south", "north"),
         recipient_basin = ifelse(row %in% c(1,2), "south", "north"),
         donor_stage = ifelse(col %in% c(1, 3), "juvenile", "adult") %>%
           factor(levels = c("juvenile","adult")),
         recipient_stage = ifelse(row %in% c(1, 3), "juvenile", "adult") %>%
           factor(levels = c("juvenile","adult")),
         type = ifelse( (col %in% c(1,2) & row %in% c(3,4)) | 
                          (col %in% c(3,4) & row %in% c(1,2)),
                        "between", "within"),
         group = interaction(donor_basin, type) %>%
           factor(levels = c("south.within",
                             "north.between",
                             "north.within",
                             "south.between"),
                  labels = c("south \nto south",
                             "north \nto south",
                             "north \nto north",
                             "south \nto north")))

# fill matrix with values
pars_cont_mat <- matrix(pars_cont$value, nrow = 29, ncol = 20)

# create tibble
pars_cont_tibble <- as_tibble(pars_cont_mat) %>%
  mutate(time = row_number(),
         year = (1991:2016)[time]) %>%
  gather(var, value, -time, -year) %>%
  mutate(var = str_split(var, "V") %>% map_int(~as.integer(.x[2]))) %>%
  filter(var < 17) %>%
  mutate(row = rep(1:4, 4)[var],
         col = rep(1:4, each = 4)[var]) %>%
  full_join(id_mat)

# annotation
labs_b <- pars_cont_tibble %>%
  tidyr::expand(recipient_stage, donor_stage) %>%
  mutate(x = 2003.5,
         y = 0.26 * 0.98,
         labels = interaction(recipient_stage, donor_stage) %>%
           factor(levels = c("juvenile.juvenile",
                             "juvenile.adult",
                             "adult.juvenile",
                             "adult.adult"),
                  labels = c("juvenile to juvenile",
                             "adult to juvenile",
                             "juvenile to adult",
                             "adult to adult")))
# plot
p2 <- pars_cont_tibble %>%
  ggplot(aes(year, value))+
  facet_grid(recipient_stage ~ donor_stage)+
  geom_text(data = labs_b,
            aes(label = labels, x = x, y = y),
            color = "black", size = 3)+
  geom_hline(yintercept = 0, alpha = 0.5, color = "gray50", size = 0.3)+
  geom_line(aes(color = group, linetype = group), size = 0.4)+
  scale_color_manual("",
                     values = c("dodgerblue","dodgerblue","gray35","gray35"),
                     guide = F)+
  scale_linetype_manual("", values = c(1,2,1,2),
                        guide = F)+
  # scale_y_continuous("Contribution to "~Delta~lambda,
  #                    limits = c(-0.26, 0.26),
  #                    breaks = c(-0.2, 0, 0.2))+
  # scale_x_continuous("Year",
  #                    limits = c(1990, 2016),
  #                    breaks = c(1994, 2003, 2012),
  #                    labels = c(1994, 2003, 2012))+
  theme(legend.key.width = unit(0.22,"in"),
        legend.key.height = unit(0.45,"in"),
        legend.text.align = 0.35,
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = margin(5,1,7,1))
p2




#==========
#========== Combine
#==========

p3 <- plot_grid(p1, p2,
                ncol = 1,
                align = "v",
                axis = "tblr",
                vjust = -0.75,
                labels = c("a","b"),
                label_size = 12)
p4 <- plot_grid(p1_leg, NULL, p3,ncol = 1, 
                rel_heights = c(1,0.25, 10))
p4

# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fig_sens.pdf",
#           width = 3.5, height = 6.5, family = "Arial")
# p4
# dev.off()







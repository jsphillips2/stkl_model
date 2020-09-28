#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(cowplot)

# import data
name <- "full"
proj_sum <-read_rds(paste0("analysis/demographic_model/model/output/",name,"/proj_sum.rds"))
var_sum <- proj_sum$var_sum




#==========
#========== Prepare panel A: Covariance matrix 
#==========

# fill covariance matrix
covar <- var_sum$covar
covar_mat <- matrix(covar$value, nrow = 20, ncol = 20)

# check diagaonl
diag(covar_mat)[1:16]

# these values are ordered first by row, then by colum. 
# the two largest values at positioins 5 and 15 correpsond to positions [1,2] and [3,4]

# create matrix of parameter identies
id_mat <- crossing(row = 1:4, col = 1:4) %>%
  arrange(col, row) %>%
  mutate(par_row = factor(paste(row, col, sep = ","),
                          levels = paste(row, col, sep = ",")),
         par_col = par_row,
         donor_basin = ifelse(col %in% c(1,2), "south", "north"),
         recipient_basin = ifelse(row %in% c(1,2), "south", "north"),
         donor_stage = ifelse(col %in% c(1, 3), "juvenile", "adult") %>%
           factor(levels = c("juvenile","adult")),
         recipient_stage = ifelse(row %in% c(1, 3), "juvenile", "adult") %>%
           factor(levels = c("juvenile","adult")),
         basin_group = interaction(donor_basin, recipient_basin) %>% 
           factor(levels = c("south.south",
                             "north.south",
                             "south.north",
                             "north.north"),
                  labels = c("south",
                             "north → south",
                             "south → north",
                             "north")),
         stage_group = interaction(donor_stage, recipient_stage) %>%
           factor(levels = c(sprintf("juvenile.juvenile"),
                             sprintf("juvenile.adult"),
                             sprintf("adult.juvenile"),
                             sprintf("adult.adult")),
                  labels = c(sprintf("j \u2192 j"),
                             sprintf("j \u2192 a"),
                             sprintf("a \u2192 j"),
                             sprintf("a \u2192 a"))))


# extract values from covar_mat corresponding to projection matrix elements
covar_mat_proj <- covar_mat[1:16, 1:16]

# add column and row names
colnames(covar_mat_proj) <- levels(id_mat$par_row)
rownames(covar_mat_proj) <- levels(id_mat$par_row)

# check names on diagonal
diag(covar_mat_proj)

# set naming order
order_names <- {id_mat %>%
    mutate(par_row = as.character(par_row)) %>%
    arrange(basin_group, stage_group)}$par_row

# reorder variance covariance matrix
covar_proj_order <- covar_mat_proj[order_names,order_names]

# check names on diagonal
diag(covar_proj_order)

# create tibble of covariance matrix elements and combine with names
covar_tibble <- tibble(z = c(covar_proj_order)) %>%
  mutate(position = row_number(),
         plot_row = rep(1:16, 16),
         plot_col = rep(1:16, each = 16))

# annotation
stage_labs <- id_mat %>%
  arrange(basin_group, stage_group) %>%
  mutate(plot_row = row_number()) %>%
  full_join(covar_tibble %>%
              filter(plot_row == plot_col)) %>%
  select(plot_row, plot_col, stage_group) %>%
  mutate(plot_row = plot_row - 1,
         plot_col = plot_col + 0.25)
basin_labs <- tibble(id = 1:4,
                     x = -1.5,
                     y = c(1,5,9,13) - 0.25,
                     yend = c(4,8,12,16) + 0.25,
                     ymid = c(2.5,6.5,10.5,14.5),
                     labels = c("south",
                                "north \u2192 south",
                                "south \u2192 north",
                                "north"))

# plot
p1 <- covar_tibble %>%
  filter(plot_row >= plot_col) %>%
  ggplot(aes(plot_col, -plot_row))+
  # covariance tiles
  geom_tile(aes(fill = sqrt(abs(z)) * z/abs(z)))+
  # stage labels - rows
  geom_text(data = stage_labs,
            aes(x = -0.5, y = -plot_row - 1,
                label = stage_group), size = 3)+
  # basin segment - rows
  geom_segment(data = basin_labs,
               aes(x = x, xend = x,
                   y = -y, yend = -yend,
                   group = id),
               size = 0.25)+
  # basin labels - rows
  geom_text(data = basin_labs,
            aes(x = -2,
                y = -ymid,
                label = labels),
            angle = 90, size = 3)+
  # stage labels - col
  geom_text(data = stage_labs,
            aes(y = 0.5 - 18, x = plot_row + 1,
                label = stage_group), 
            angle = 90, size = 3)+
  # basin segment - col
  geom_segment(data = basin_labs,
               aes(y = -x - 20, yend = -x - 20,
                   x = y, xend = yend,
                   group = id),
               size = 0.25)+
  # basin labels - col
  geom_text(data = basin_labs,
            aes(y = 2 - 21,
                x = ymid,
                label = labels),
            size = 3)+
  # fill scale
  scale_fill_gradient2("Covariance",
                       low = "navy",high = "firebrick", mid = "gray90",
                       midpoint = 0,
                       breaks = c(-0.22, 0, 0.22, 0.5),
                       labels = c(-0.05, 0, 0.05, 0.25),
                       limits = c(-0.24,0.52))+
  scale_x_continuous(limits = c(-2,17.5))+
  coord_equal()+
  theme_void()+
  theme(plot.margin = margin(c(0, 0, 5, 0)),
        legend.position = c(0.86,0.8),
        legend.title = element_text(size = 10))

p1



#==========
#========== Panel B: Variance partitiong
#==========

# extract contribution matrix and define tibble
cont_mat <- var_sum$cont %>%
  mutate(value = value/sum(value)) %>%
  filter(id < 17) %>%
  full_join(id_mat %>%
              select(donor_basin, recipient_basin, basin_group, stage_group) %>%
              mutate(id = row_number())) %>%
  full_join(crossing(plot_col = 1:4,
                     plot_row = 1:4) %>%
              mutate(id = row_number()))

# summarize by basin
cont_mat_sum <- cont_mat %>%
  group_by(donor_basin, recipient_basin) %>%
  summarize(value = sum(value),
            plot_col = mean(plot_col) + 0.5,
            plot_row = mean(plot_row) + 0.5) %>%
  mutate(label = paste0(100*round(value,2),"%"))

# annotation
stage_labs_b <- tibble(x = c(0.25, 0.25, 0.25, 0.25,
                             1, 2, 3, 4),
                       y = c(-1, -2, -3, -4,
                             -0.25, -0.25, -0.25, -0.25),
                       labels = rep(c("j","a"), 4))
basin_labs_b <- tibble(x = c(-0.25, -0.25,
                             1.5, 3.5),
                       y = c(-1.5, -3.5,
                             0.25, 0.25),
                       angle = c(90, 90, 0, 0),
                       labels = rep(c("south","north"), 2))
basin_segs_b <- tibble(x = c(0, 0, 0.9, 2.9),
                           xend = c(0, 0, 2.1, 4.1),
                           y = c(-0.9, -2.9, 0, 0),
                           yend = c(-2.1, -4.1, 0, 0),
                           id = c(1,2,3,4))
donor_labs_b <- tibble(x = c(-0.75, 2.5),
                       y = c(-2.5, 0.75),
                       labels = c("recipient", "donor"),
                       angle = c(90, 0))
donor_segs_b <- tibble(x = c(-0.5, 1.5),
                           xend = c(-0.5, 3.5),
                           y = c(-1.5, 0.5),
                           yend = c(-3.5, 0.5),
                           id = c(1,2))

# define breaks
bl <- c(log(2/100), log(6/100), log(18/100)) 

# plot
p2 <- cont_mat %>%
  ggplot(aes(plot_col, -plot_row))+
  # contribution tiles
  geom_tile(aes(fill = log(value)))+
  # stage labels
  geom_text(data = stage_labs_b,
            aes(x = x, y = y, label = labels),
            size = 3)+
  # basin labels
  geom_text(data = basin_labs_b,
            aes(x = x, y = y, label = labels, angle  = angle),
            size = 3)+
  # donor labels
  geom_text(data = donor_labs_b,
            aes(x = x, y = y, label = labels, angle  = angle),
            size = 3)+
  # basin segments
  geom_segment(data = basin_segs_b,
               aes(x = x, xend = xend, y = y, yend = yend),
               size = 0.25)+
  # donor segments
  geom_segment(data = donor_segs_b,
               aes(x = x, xend = xend, y = y, yend = yend),
               size = 0.25)+
  # total contributions
  geom_text(data  = cont_mat_sum, 
            aes(x = plot_col - 0.5,
                y = -plot_row + 0.5,
                label = label, fontface = "bold"),
            size = 3)+
  # central segments
  geom_segment(aes(x = 0.5, xend = 4.5,
                   y = -2.5, yend = -2.5),
               size = 0.25)+
  geom_segment(aes(x = 2.5, xend = 2.5,
                   y = -0.5, yend = -4.5),
               size = 0.25)+
  # fill scale
  scale_fill_gradient2("Contribution",
                       low = "navy",high = "firebrick", mid = "gray100",
                       midpoint = -4,
                       limits = c(-4.11, -1.51),
                       breaks = bl,
                       labels = paste0(100*round(exp(bl), 2), "%"))+
  coord_equal()+
  theme_void()+
  theme(plot.margin = margin(c(5, 0, 0, 0)),
        legend.position = c(1.25,0.4),
        legend.title = element_text(size = 10))
p2




#==========
#========== Assemble panels and save
#==========

p3 <- plot_grid(p1, p2, nrow = 2, 
          rel_heights = c(2, 1), 
          labels = c("a",
                     "b"),
          align = "hv",
          axis = "tblr",
          label_size = 12,
          hjust = c(-5, -5),
          vjust = c(1.1,1)
          )
p3


# cairo_pdf(file = "analysis/demographic_model/model/figures/figs/fig_cov.pdf",
#           width = 6, height = 8, family = "Arial")
# p3
# dev.off()

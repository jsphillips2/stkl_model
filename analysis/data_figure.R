#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(lubridate)
library(cowplot)
library(lemon)
library(ggsn)

# import abundance data
abundance <- read_csv("data/counts89_20.csv") %>%
  filter(station != "BV") %>%
  mutate(basin = ifelse(basin == "east", "south", basin),
         basin = factor(basin, levels = c("north","south")),
         size_class = factor(size_class, levels = c("small","large")),
         day_night = factor(day_night, levels = c("D","N")),
         count = ifelse(year==2016 & is.na(count), 0, count),
         date = ymd(paste(year, month, 1, sep = "_")),
         l_group = as.numeric(factor(paste(basin, size_class, date)))
  ) %>%
  mutate(stage = factor(size_class,
                        levels = c("small",
                                   "large"),
                        labels = c("juvenile",
                                   "adult")),
         state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "juvenile.north",
                                   "adult.south",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nsouth",
                                   "adult\nnorth")),
         station = factor(station,
                          levels = c("23","27","41","44","135","124","128","DN")))

# import fits
fit_sum <- read_csv("data/fit.csv")

# import betas
betas <- read_csv("data/betas.csv")

# define datebreaks / limits
date_breaks <- lubridate::as_date(c("1995-06-01",
                                    "2005-06-01",
                                    "2015-06-01"))
date_limits <-  lubridate::as_date(c("1990-08-01",
                                     "2020-08-01"))


# define datebreaks / limits
year_breaks <- c(1995, 2005, 2015)
year_limits = c(1990, 2020)

# vector of dates for matching with time id's
date_match <- abundance$date %>% unique()

# define station colors
station_colors <- c("skyblue3","firebrick3","magenta4","darkorange","goldenrod",
                    "purple","darkgreen","dodgerblue2")
names(station_colors) <- c("23","27","41","44","135","124","128","DN")

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(t = 1,
                                       r = 1,
                                       b = 1,
                                       l = 1),
                  legend.margin = margin(t = 0,
                                         r = 0,
                                         b = 0,
                                         l = -4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black",family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,5,0,0)),
                  axis.title.x = element_text(margin = margin(5,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

#=========================================================================================





#=========================================================================================
#========== Detection-corrected fit
#=========================================================================================

# beta intercept
beta_int <- {betas %>%
    filter(name %in% c("(Intercept)"))}$median

# beta night
base_night <- {betas %>%
    filter(name %in% c("day_nightN","size_classlarge:day_nightN")) %>%
    summarize(val = mean(median))}$val

# beta basin (average over station)
beta_basin <- betas %>%
  filter(!(name %in% c("(Intercept)","day_nightN","size_classlarge:day_nightN",
                       "size_classlarge"))) %>%
  mutate(station = strsplit(name, "station") %>% map_chr(~.x[2])) %>%
  select(station, median) %>%
  left_join(abundance %>%
              tidyr::expand(nesting(station, basin))) %>%
  group_by(basin) %>%
  summarize(beta_basin = mean(median))

# beta by size_class
beta_size <- betas %>%
  filter(name %in% c("size_classlarge")) %>%
  select(mean) %>%
  rename(beta_size = mean) %>%
  mutate(size_class = "large") %>%
  bind_rows(tibble(size_class = "small",
                   beta_size = 0))

# extract fit
l_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "lam\\[")) %>%
  mutate(l_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower95, lower68, median, upper68, upper95, l_group) %>%
  full_join(abundance %>%
              expand(nesting(date, basin, size_class, l_group))) %>%
  arrange(size_class, basin, date)


# calculate expected catch
trap_fit <- l_fit %>%
  select(date, basin, size_class, median) %>%
  left_join(beta_basin) %>%
  left_join(beta_size) %>%
  mutate(alpha = beta_int + base_night + beta_basin + beta_size,
         phi = exp(alpha) / (1 + exp(alpha)),
         est = phi * median) %>%
  select(date, basin, size_class, phi, est) %>%
  mutate(stage = factor(size_class,
                        levels = c("small",
                                   "large"),
                        labels = c("juvenile",
                                   "adult")),
         state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "juvenile.north",
                                   "adult.south",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nsouth",
                                   "adult\nnorth")))


#=========================================================================================




#=========================================================================================
#========== Plot map
#=========================================================================================

# import files
source("analysis/map_prep.R")
stations <- read_csv("data/stations.csv") %>%
  mutate(station = factor(station,
                          levels = c("23","27","41","44","135","124","128","DN")))

# plto annotation setup
anchor_scale = c(0403100, 7281500)
anchor_north = c(0405500, 7281000)
names(anchor_scale) = c("x","y")
names(anchor_north) = c("x","y")

# plot
p_map <- myv_df %>%
  ggplot(aes(long,lat))+
  geom_polygon(aes(fill = piece), size = 0.3, color = "black")+
  coord_equal()+
  theme_void() + 
  theme(plot.margin = margin(0,0,0,0))+
  geom_text(inherit.aes = F,
            data = stations %>%
              mutate(long = ifelse(station == "DN", long + 250, long),
                     lat = ifelse(station == "DN", lat + 100, lat),
                     lat = ifelse(station == "128", lat + 100, lat)),
            fontface = "bold",
            aes(x = long,
                y = lat,
                label = station,
                color = station),
            size = 2.5)+
  scale_fill_manual("",values = c("gray95",rep("white", 18)), guide = F)+
  scale_color_manual("",
                     values = station_colors,
                     guide = F)+
  scalebar(myv_df, dist = 2, dist_unit = "km",st.bottom = F,
           location = "topleft",st.size = 2.5, st.dist = 0.03,
           border.size = 0.2,anchor = anchor_scale,
           transform = F, model = "WGS84")+
  north(myv_df,symbol = 3,
        anchor = anchor_north)+
  theme(plot.margin = margin(t = 0))
p_map



#=========================================================================================





#=========================================================================================
#========== Plot trapping data
#=========================================================================================

# plot labels
labs <- abundance %>%
  tidyr::expand(state) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 570)

p_trap <- ggplot(data = abundance,
                 aes(x = date,
                     y = count))+
  facet_rep_wrap(~state)+
  geom_text(data = labs,
            aes(label = state,
                x = x,
                y = y),
            color = "black",
            size = 3.2,
            inherit.aes = F)+
  geom_point(aes(color = station),
             size = 0.6,
             alpha = 0.5,
             shape = 16)+
  geom_line(data = trap_fit,
            aes(y = est),
            size = 0.5)+
  scale_y_continuous(name = Trap~catch,
                     limits = c(0, 600),
                     breaks = c(0, 300, 600))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks[c(1,3)],
               labels = year_breaks[c(1,3)])+
  scale_color_manual("",
                     values = station_colors,
                     guide = F)+
  theme(legend.key.height = unit(1, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.text = element_text(margin = margin(l = -10)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.9, "lines"),
        plot.margin = margin(r = 5),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "top", 
                    bottom='none',
                    gap = 0)

# examine
p_trap



#=========================================================================================





#=========================================================================================
#========== Combine
#=========================================================================================

p_data <- plot_grid(p_map, NULL,p_trap,
                  nrow = 3,
                  rel_heights = c(0.8, 0.1, 1),
                  ncol = 1,
                  labels = c("a",
                             "",
                             "b"),
                  label_size = 12,
                  label_fontface = "plain",
                  hjust = c(-0.5, 0, -0.5),
                  vjust = c(1,0,-2)
)


# examine
p_data


# export
# cairo_pdf(file = "analysis/figures/fig_data.pdf",
#           width = 3.5, height = 7, family = "Arial")
# p_data
# dev.off()


#=========================================================================================
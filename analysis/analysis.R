#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(cowplot)
library(lemon)
library(rstan)
library(loo)
library(WaveletComp)

# set cores
options(mc.cores = parallel::detectCores()-4)

# import trapping data and model population estimates



# import model objects
fit_full <- read_rds("output/fit_full.rds")
fit_no_juv_move <- read_rds("output/fit_no_juv_move.rds")
fit_no_move <- read_rds("output/fit_no_move.rds")
fit_null <- read_rds("output/fit_null.rds")

# import demographic projection
proj_no_juv_move <- read_rds("output/proj_no_juv_move.rds")

# define named list with model
models <- list(fit_full,
               fit_no_juv_move,
               fit_no_move,
               fit_null)
model_names <- c("fit_full",
                 "fit_no_juv_move",
                 "fit_no_move",
                 "fit_null")
names(models) <- model_names

# define datebreaks / limits
date_breaks <- lubridate::as_date(c("1995-06-01",
                                    "2005-06-01",
                                    "2015-06-01"))
date_limits <-  lubridate::as_date(c("1990-08-01",
                                     "2020-08-01"))


# define datebreaks / limits
year_breaks <- c(1995, 2005, 2015)
year_limits = c(1990, 2020)

# extract data 
data_prep <- fit_full$data_prep %>%
  mutate(stage = factor(stage, 
                        levels = c("small","large"),
                        labels = c("juvenile","adult"))) %>%
  mutate(state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "juvenile.north",
                                   "adult.south",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nsouth",
                                   "adult\nnorth")))

# data frame for matching stages and basins to id's
state_match <- data_prep %>% 
  tidyr::expand(nesting(state, basin, stage))  %>%
  arrange(basin, stage) %>%
  mutate(st = row_number(),
         state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "juvenile.north",
                                   "adult.south",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nsouth",
                                   "adult\nnorth")))

# vector of dates for matching with time id's
date_match <- data_prep$date %>% unique()

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
#========== Model comparison
#=========================================================================================

# # extract LOOIC and posterior log likelihood
# model_compare <- model_names %>%
#   lapply(function(x){
#     xx = models[x]
#     tibble(model = x,
#            looic = loo(xx[[1]]$fit)$estimates["looic","Estimate"],
#            log_lik = {xx[[1]]$fit_summary %>% 
#                filter(str_detect(var, "log_lik_sum"))}$`50%`
#     )
#   }) %>%
#   bind_rows()
# 
# write_csv(model_compare, "output/model_compare.csv")

# read survival
# model_compare <- read_csv("output/model_compare.csv")

#=========================================================================================





#=========================================================================================
#========== Stage-transition
#=========================================================================================

# # extract
ids <- proj_no_juv_move [[1]]$setup$ids
# extract
trans_full <- parallel::mclapply(ids, function(i_){
  proj_ = proj_no_juv_move [[i_]]$proj
  tibble(gg = proj_$GG[2,1,1]) 
}) %>%
  bind_rows()
  
# summarize
trans_sum_write <- trans_full %>%
  summarize(lo = quantile(gg, probs = c(0.16)),
            mi = median(gg),
            hi = quantile(gg, probs = c(0.84))) %>%
  ungroup()

trans_sum_write

#=========================================================================================





#=========================================================================================
#========== Model fits
#=========================================================================================

# extract scaling parameter
y_scale <- fit_full$y_scale

# extract fits
# fits_write <- model_names %>%
#   lapply(function(x){
#     xx = models[x]
#     xxx = xx[[1]]$fit_summary
#     xxx %>%
#       select(var, `16%`, `50%`, `84%`) %>%
#       rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
#       filter(str_detect(xxx$var, "x"), !str_detect(xxx$var, "x0")) %>%
#       mutate(model = x,
#              name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
#              st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#              time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
#              date = date_match[time]) 
#   }) %>%
#   bind_rows()
# 
# write_csv(fits_write, "output/fits_write.csv")

# read survival
fits <- read_csv("output/fits_write.csv") %>%
  mutate(model = factor(model,
                        levels = c("fit_full",
                                   "fit_no_juv_move",
                                   "fit_no_move",
                                   "fit_null"),
                        labels = c("full model",
                                   "adult dispersal",
                                   "no dispersal",
                                   "fixed rates"))) %>%
  full_join(state_match) %>%
  select(model, basin, state, stage, date, name, lo, mi, hi)

# plot labels
labs <- fits %>%
  tidyr::expand(state) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 11)

# plot
p_fit <- ggplot(data = fits,
             aes(x = date,
                 y = mi))+
  facet_rep_wrap(~state)+
  geom_text(data = labs,
            aes(label = state,
                x = x,
                y = y),
            color = "black",
            size = 3.2,
            inherit.aes = F)+
  geom_point(data = data_prep,
             aes(y = mean_scale / y_scale),
             shape = 21,
             size = 0.8,
             alpha = 0.4,
             fill = "black",
             stroke = 0.2)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi,
                  fill = model),
              alpha = 0.2)+
  geom_line(aes(color = model),
            size = 0.3,
            alpha = 1)+
  scale_y_continuous(name = Relative~abundance,
                     limits = c(0, 12),
                     breaks = c(0, 4, 8, 12))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  scale_color_manual("",
                     values = c("black", "firebrick2","royalblue3","orange1"))+
  scale_fill_manual("",
                     values = c("gray70", "firebrick2","royalblue3","orange1"))+
    
  theme(legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = "top",
        legend.text = element_text(margin = margin(l = -10)),
        legend.key.size = unit(0.7, "lines"),
        legend.spacing.x = unit(0.9, "lines"),
        plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
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
p_fit

# export
# cairo_pdf(file = "analysis/figures/fig_fit.pdf",
#           width = 3.5, height = 3.5, family = "Arial")
# p_fit
# dev.off()


#=========================================================================================





#=========================================================================================
#========== Recruitment
#=========================================================================================

# # extract
# ids <- proj_no_juv_move [[1]]$setup$ids
# # extract
# rec_full <- parallel::mclapply(ids, function(i_){
#   proj_ = proj_no_juv_move [[i_]]$proj
#   RR_ = proj_$RR
#   x_ = proj_$x
#   lapply(1:dim(RR_)[3], function(t_){
#     tibble(id = i_,
#            rec = c(RR_[c(1,3),c(2,4),t_])[c(1,4)],
#            x = x_[c(2,4), t_],
#            date = date_match[t_]
#     ) %>%
#       filter(rec > 0) %>%
#       group_by(date) %>%
#       mutate(basin = 1:2,
#              basin = factor(basin,
#                             levels = c(1,2),
#                             labels = c("south","north")))
#   }) %>%
#     bind_rows()
# }) %>%
#   bind_rows() %>%
#   ungroup()
# 
# # summarize
# rec_sum_write <- rec_full %>%
#   group_by(date,  basin) %>%
#   summarize(lo = quantile(rec, probs = c(0.16)),
#             mi = median(rec),
#             hi = quantile(rec, probs = c(0.84))) %>%
#   ungroup()
# 
# write_csv(rec_sum_write, "output/rec_sum_write.csv")

# read recruitment
rec_sum <- read_csv("output/rec_sum_write.csv") %>%
  mutate(basin = factor(basin,
                        levels = c("south", "north"),
                        labels = c("south","north")))

# plot labels
labs <- rec_sum %>%
  tidyr::expand(basin) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 9)

# plot
p_rec_a <- ggplot(data = rec_sum,
                 aes(x = date,
                     y = mi))+
  facet_rep_wrap(~basin,
                 nrow = 1)+
  geom_text(data = labs,
            aes(label = basin,
                x = x,
                y = y),
            color = "black",
            size = 3.2,
            inherit.aes = F)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)+
  geom_line(size = 0.4)+
  scale_y_continuous(name = Recruitment~catpia^{-1},
                     limits = c(0, 10),
                     breaks = c(0, 4, 8))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
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
p_rec_a


#  covariance matrix
r_covmat <- rec_sum %>%
  select(date, basin,  mi) %>%
  pivot_wider(names_from =  basin, values_from = mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()
r_covmat <- (r_covmat / mean(r_covmat)) 
r_covmat_d <- r_covmat %>% 
  as_tibble() %>%
  mutate(row = c("south","north")) %>%
  pivot_longer(cols = c(south,north)) %>%
  mutate(x = as.numeric(factor(row,levels = c("south","north"))),
         y = as.numeric(factor(name,levels = c("south","north")))) %>%
  filter(!(x == 2 & y == 1))

# plot covariance
p_rec_b <- ggplot(data = r_covmat_d,
                    aes(x = x,
                        y = -y,
                        fill = value))+
  geom_tile()+
  geom_text(aes(x = x,
                y = -y,
                label = round(value, 1)),
            size = 3.5,
            fontface = "plain")+
  scale_x_continuous("",
                     breaks = c(1, 2),
                     labels = c("south",
                                "north"))+
  scale_y_continuous("",
                     breaks = c(-1, -2),
                     labels = c("south",
                                "north"))+
  scale_fill_gradient(name = "",
                       low="gray90", 
                       high="firebrick", 
                       limits = c(0.5, 1.7),
                       guide = F)+
  coord_equal()+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine
p_rec_b

# combine
p_rec <- plot_grid(NULL, p_rec_a, NULL, p_rec_b,
                   ncol = 1,
                   rel_heights = c(0.05, 1.3, 0.15, 1),
                   align = "v",
                   axis = "tblr",
                   labels = c("",
                              "a",
                              "",
                              "b"),
                   label_size = 12,
                   label_fontface = "plain",
                   hjust = c(0, 0, 0),
                   vjust = c(0,0.1,0,0))

# examine
p_rec

# export
# cairo_pdf(file = "analysis/figures/fig_rec.pdf",
#           width = 3.5, height = 4.3, family = "Arial")
# p_rec
# dev.off()


#=========================================================================================





#=========================================================================================
#========== Survival
#=========================================================================================

# # extract
# ids <- proj_no_juv_move [[1]]$setup$ids
# surv_full <- parallel::mclapply(ids, function(i_){
#   proj_ = proj_no_juv_move [[i_]]$proj
#   QQ_ = proj_$QQ
#   QQs_ = proj_$QQs
#   x_ = proj_$x
#   lapply(1:dim(QQ_)[3], function(t_){
#     tibble(id = i_,
#            surv = diag(QQ_[,,t_]),
#            x = x_[, t_],
#            date = date_match[t_],
#            st = 1:4
#     ) 
#   }) %>%
#     bind_rows()
# }) %>%
#   bind_rows() 
# 
# # summarize
# surv_sum_write <- surv_full %>%
#   group_by(date, st) %>%
#   summarize(lo = quantile(surv, probs = c(0.16)),
#             mi = median(surv),
#             hi = quantile(surv, probs = c(0.84))) %>%
#   ungroup()
# 
# write_csv(surv_sum_write, "output/surv_sum_write.csv")

# read survival
surv_sum <- read_csv("output/surv_sum_write.csv") %>%
  full_join(state_match) %>%
  select(-st)
  

# plot labels
labs <- surv_sum %>%
  tidyr::expand(state) %>%
  mutate(x = lubridate::as_date("2005-07-01"),
         y = 1.2)

# plot
p_surv_a <- ggplot(data = surv_sum,
                aes(x = date,
                    y = mi))+
  facet_rep_wrap(~state)+
  geom_hline(yintercept = 0.5,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_text(data = labs,
            aes(label = state,
                x = x,
                y = y),
            color = "black",
            size = 3.2,
            inherit.aes = F)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)+
  geom_line(size = 0.4)+
  scale_y_continuous(name = Survival~probability,
                     limits = c(0, 1.3),
                     breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1"))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "top", 
                    bottom='none',
                    gap = 0)


# examine
p_surv_a

#  covariance matrix
s_covmat <- surv_sum %>%
  select(date, state,  mi) %>%
  pivot_wider(names_from =  state, values_from = mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()
s_covmat <- (s_covmat / mean(s_covmat)) 
s_covmat <- s_covmat %>% 
  as_tibble() %>%
  mutate(row = c("juvenile\nsouth",
                 "adult\nsouth",
                 "juvenile\nnorth",
                 "adult\nnorth")) %>%
  pivot_longer(cols = c("juvenile\nsouth",
                        "adult\nsouth",
                        "juvenile\nnorth",
                        "adult\nnorth")) %>%
  mutate(x = as.numeric(factor(row,levels = c("juvenile\nsouth",
                                              "adult\nsouth",
                                              "juvenile\nnorth",
                                              "adult\nnorth"))),
         y = as.numeric(factor(name,levels = c("juvenile\nsouth",
                                               "adult\nsouth",
                                               "juvenile\nnorth",
                                               "adult\nnorth")))) %>%
  filter(!(x > y))

# plot covariance
p_surv_b <- ggplot(data = s_covmat,
                    aes(x = x,
                        y = -y,
                        fill = value))+
  geom_tile()+
  geom_text(aes(x = x,
                y = -y,
                label = round(value, 1)),
            size = 3.5,
            fontface = "plain")+
  scale_x_continuous("",
                     breaks = c(1:4),
                     labels = c("juvenile\nsouth",
                                "adult\nsouth",
                                "juvenile\nnorth",
                                "adult\nnorth"))+
  scale_y_continuous("",
                     breaks = -c(1:4),
                     labels = c("juvenile\nsouth",
                                "adult\nsouth",
                                "juvenile\nnorth",
                                "adult\nnorth"))+
  scale_fill_gradient(name = "",
                      low="gray90", 
                      high="firebrick", 
                      limits = c(0.5, 1.65),
                      guide = F)+
  coord_equal()+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine
p_surv_b


# combine
p_surv <- plot_grid(NULL, p_surv_a, NULL, p_surv_b,
                   ncol = 1,
                   rel_heights = c(0.05, 1.3, 0.15, 1),
                   align = "v",
                   axis = "tblr",
                   labels = c("",
                              "a",
                              "",
                              "b"),
                   label_size = 12,
                   label_fontface = "plain",
                   hjust = c(0, 0, 0),
                   vjust = c(0,0,0,0))

# examine
p_surv

# export
# cairo_pdf(file = "analysis/figures/fig_surv.pdf",
#           width = 3.5, height = 6.75, family = "Arial")
# p_surv
# dev.off()



#=========================================================================================




#=========================================================================================
#========== Dispersal
#=========================================================================================

# # extract
# ids <- proj_no_juv_move [[1]]$setup$ids
# disp_full <- parallel::mclapply(ids, function(i_){
#   proj_ = proj_no_juv_move[[i_]]$proj
#   QQ_ = proj_$QQ
#   GG_ = proj_$GG
#   DD_ = proj_$DD
#   x_ = proj_$x
#   lapply(1:dim(DD_)[3], function(t_){
#     P_ = DD_[,,t_] %*% GG_[,,t_] %*%  QQ_[,,t_]
#     tibble(id = i_,
#            basin = 1:2,
#            disp = 1 - diag(DD_[,,t_])[c(2,4)],
#            move = c(sum(P_[3:4, 1:2] %*% x_[1:2, t_]),
#                     sum(P_[1:2, 3:4] %*% x_[3:4, t_])),
#            date = date_match[t_]
#     ) 
#   }) %>%
#     bind_rows()
# }) %>%
#   bind_rows() %>%
#   ungroup() %>%
#   mutate(basin = factor(basin,
#                         levels = c(1,2),
#                         labels = c("south","north")))
# 
# # summarize
# disp_sum <- disp_full %>%
#   gather(var, val, disp, move) %>%
#   group_by(date, basin,var) %>%
#   summarize(lo = quantile(val, probs = c(0.16)),
#             mi = median(val),
#             hi = quantile(val, probs = c(0.84)))
# 
# # calculate net and cumulative movement
# disp_net <- disp_full %>%
#   select(id, date, basin, move) %>%
#   spread(basin, move) %>%
#   group_by(id) %>%
#   arrange(id, date) %>%
#   mutate(net = south - north) %>%
#   select(-south, -north) %>%
#   group_by(date) %>%
#   summarize(lo = quantile(net, probs = c(0.16), na.rm = T),
#             mi = median(net, na.rm = T),
#             hi = quantile(net, probs = c(0.84), na.rm = T)) %>%
#   ungroup()
# 
# write_csv(disp_sum, "output/disp_sum.csv")
# write_csv(disp_net, "output/disp_net.csv")

# read dispersal
# disp_sum <- read_csv("output/disp_sum.csv") %>%
#   mutate(basin = factor(basin,
#                         levels = c("south","north"),
#                         labels = c("south","north")))
# disp_net <- read_csv("output/disp_net.csv")

# plot labels
labs <- tibble(x = lubridate::as_date("1993-07-01"),
               y = c(-3.3, 1),
               label = c("southward",
                         "northward"))

# plot
p_disp <- ggplot(data = disp_net,
                 aes(x = date, 
                     y = mi))+
  geom_hline(yintercept = 0,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)+
  geom_line(size = 0.4)+
  geom_text(data = labs,
            aes(label = label,
                x = x,
                y = y),
            color = "black",
            size = 3.2,
            inherit.aes = F)+
  scale_y_continuous(name = Net~dispersal,
                     breaks = c(-3, -2, -1, 0, 1),
                     limits = c(-3.35, 1.5))+
  scale_x_date(name = "Date",
               limits = date_limits,
               breaks = date_breaks,
               labels = year_breaks)+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
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
p_disp

# export
# cairo_pdf(file = "analysis/figures/fig_disp.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p_disp
# dev.off()


#=========================================================================================





#=========================================================================================
#========== Lambda
#=========================================================================================

# extract
# ids <- proj_no_juv_move [[1]]$setup$ids
# years <- proj_no_juv_move [[1]]$setup$years
# lam_full <- parallel::mclapply(ids, function(i_){
#   sens_ = proj_no_juv_move[[i_]]$sens
#   lapply(1:length(sens_), function(y_){
#     tibble(id = i_,
#            year = unique(years)[y_],
#            l = sens_[[y_]]$l,
#            l_asym = sens_[[y_]]$l_asym,
#            r = sens_[[y_]]$r,
#            r_asym = sens_[[y_]]$r_asym 
#     )
#   }) %>%
#     bind_rows()
# }) %>%
#   bind_rows() %>%
#   gather(var, val, l, l_asym, r, r_asym) 
# 
# # summarize
# lam_sum <- lam_full %>%
#   group_by(year, var) %>%
#   summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
#             mi = median(val, na.rm = T),
#             hi = quantile(val, probs = c(0.84), na.rm = T)) %>%
#   ungroup()
# 
# # summarize median lambda
# lam_sum_all <- lam_full %>%
#   filter(var %in% c("l","l_asym")) %>%
#   group_by(id, var) %>%
#   summarize(geom = prod(val)^(1/length(val)),
#             arith = mean(val)) %>%
#   ungroup() %>%
#   gather(type, val, geom, arith) %>%
#   group_by(var, type) %>%
#   summarize(lo = quantile(val, probs = c(0.16), na.rm = T),
#             mi = median(val, na.rm = T),
#             hi = quantile(val, probs = c(0.84), na.rm = T))
# 
# write_csv(lam_sum, "output/lam_sum.csv")
# write_csv(lam_sum_all, "output/lam_sum_all.csv")

lam_sum <- read_csv("output/lam_sum.csv")
lam_sum_all <- read_csv( "output/lam_sum_all.csv")

# extract relevant data
wave_prep <- lam_sum %>%
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
wave_d <- as_tibble(as.data.frame(wave_r$Power)) %>%
  mutate(period = wave_r$Period) %>%
  gather(id, power, -period) %>%
  mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
         year = 1990 + id,
         type = "Transient") %>%
  bind_rows(as.data.frame(as_tibble(wave_asym$Power)) %>%
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


# plot panel a
p_lam_a <- ggplot(data = lam_sum %>%
                    filter(var %in% c("l","l_asym")) %>%
                    mutate(type = factor(var,
                                         levels = c("l","l_asym"),
                                         labels = c("Transient"," Asymptotic"))),
                  aes(x = year, 
                      y = mi))+
  facet_wrap(~type)+
  geom_hline(yintercept = 1,
             size = 0.2,
             color = "black",
             linetype = 2)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)+
  geom_line(size = 0.4)+
  scale_y_continuous(Population~growth~rate~(lambda),
                     trans = "log",
                     breaks = c(1/9, 1/3, 1, 3),
                     labels = c("1/9", "1/3", "1", "3"))+
  scale_x_continuous(name = "",
               limits = year_limits,
               breaks = year_breaks,
               labels = NULL)+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(0.25, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine
p_lam_a
  
# plot panel b
p_lam_b <- ggplot(data = wave_d,
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
                     limits = year_limits,
                     breaks = year_breaks)+
  scale_fill_gradient2(name = "",
                       low ="dodgerblue", 
                       mid="gray90", 
                       high="firebrick", 
                       midpoint = 0.1,
                       limits = c(0, 1.3),
                       breaks = c(0.1, 1.2),
                       guide = F)+
  theme(strip.text = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0.25, "lines"),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  coord_capped_cart(left = "top", 
                    bottom='none',
                    gap = 0)

# examine plot
p_lam_b

# combine
p_lam <- plot_grid(NULL, p_lam_a, NULL, p_lam_b,
                ncol = 1,
                rel_heights = c(0.05, 1, 0.075, 1),
                align = "v",
                axis = "tblr",
                labels = c("",
                           "a",
                           "",
                           "b"),
                label_size = 12,
                label_fontface = "plain",
                hjust = c(0, 0, 0),
                vjust = c(0,0,0,0))

# examine
p_lam

# export
# cairo_pdf(file = "analysis/figures/fig_lam.pdf",
#           width = 5, height = 4.25, family = "Arial")
# p_lam
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Elasticity
#=========================================================================================

# extract
# years <- proj_no_juv_move [[1]]$setup$years
# ids <- proj_no_juv_move [[1]]$setup$ids
# theta_names <- proj_no_juv_move[[1]]$setup$theta_names
# 
# elast <- parallel::mclapply(ids, function(i_){
#   sens_ = proj_no_juv_move[[i_]]$sens
#   lapply(1:length(sens_), function(y_){
#     tibble(id = i_,
#            year = unique(years)[y_],
#            elas = sens_[[y_]]$l_elas_p
#     ) %>%
#       mutate(theta = theta_names[row_number()])
#   }) %>%
#     bind_rows()
# }) %>%
#   bind_rows() 
# 
# 
# # summarize
# elast_sum <- elast %>%
#   mutate(name = str_split(theta, "") %>% map_chr(~as.character(.x[[1]])),
#          index = str_split(theta, "") %>% map_chr(~as.character(.x[[2]])),
#          name = paste0(name, index),
#          season = str_split(theta, "") %>% map_chr(~as.character(.x[[3]]))) %>%
#   group_by(id, year, name, index) %>%
#   summarize(elas = sum(elas)) %>%
#   group_by(year, name, index) %>%
#   summarize(lo = quantile(elas, probs = c(0.16), na.rm = T),
#             mi = median(elas, na.rm = T),
#             hi = quantile(elas, probs = c(0.84), na.rm = T)) %>%
#   ungroup()
# 
# 
# # mean across time
# elast_sum_mean <- elast_sum %>%
#   group_by(name) %>%
#   summarize(mi = mean(mi)) %>%
#   arrange(-mi) %>%
#   ungroup()
# 
# write_csv(elast_sum, "output/elast_sum.csv")
# write_csv(elast_sum_mean, "output/elast_sum_mean.csv")

elast_sum <- read_csv("output/elast_sum.csv") %>%
  filter(name != "g1")
elast_sum_mean <- read_csv("output/elast_sum_mean.csv") %>%
  filter(name != "g1")

# set parameter order
theta_order <- elast_sum_mean$name 
theta_labs <- c(expression("surv"["j,n"]),
                expression("recr"[n]),
                expression("surv"["a,n"]),
                # expression("trans"),
                expression("surv"["a,s"]),
                expression("disp"[s%->%n]),
                expression("recr"[s]),
                expression("surv"["j,s"]),
                expression("disp"[n%->%s]))

# plot
p_elas <- ggplot(data = elast_sum %>%
                   ungroup() %>%
                   mutate(name = factor(name, levels = theta_order),
                          pos = as.numeric(name) - 
                            0.25 * (year - mean(year)) / max((year - mean(year)))) ,
                 aes(pos, mi))+
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
                  breaks = 1:8,
                  labels = theta_labs)+
  scale_y_continuous(Elasticity~of~lambda~(transient),
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
p_elas

# export
# cairo_pdf(file = "analysis/figures/fig_elas.pdf",
#           width = 3.5, height = 4, family = "Arial")
# p_elas
# dev.off()

#=========================================================================================
# load packages
library(tidyverse)
library(Matrix)
library(matrixcalc)

options(mc.cores = parallel::detectCores()-4)


# import model
out_in <- read_rds("output/fit_full.rds")

# extract data
data_prep <- out_in$data_prep
data_list <- out_in$data_list
fit <- out_in$fit 
fit_summary <- out_in$fit_summary 
pj <- data_list$pj
nt <- data_list$nt
n <- data_list$n
b <- data_list$b
dates <- unique(data_prep$date)
years_all <- lubridate::year(dates)
years <- years_all[1:(length(years_all) - 2)]

# paramter names for extraction
vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rt") |
           str_detect(fit_summary$var, "rf") |
           str_detect(fit_summary$var, "qt") |
           str_detect(fit_summary$var, "qf") |
           str_detect(fit_summary$var, "gt") |
           str_detect(fit_summary$var, "gf") |
           str_detect(fit_summary$var, "dt") |
           str_detect(fit_summary$var, "df") | 
           str_detect(fit_summary$var, "x0"))}$var

# extract fit
extract_full <-  rstan::extract(fit, pars = vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(vars) %>%
  sample_n(1000) %>%
  mutate(id = row_number()) %>%
  gather(var, val, -id) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])))

# select row
id_ <- 50

extract <- extract_full %>%
  filter(id == id_)







plot(colSums(x), type = "l", ylim = c(0, 17))
points(colSums(data_list$y))

years_unique <- unique(years)
ny <- length(years_unique)


aa <- lapply(1 : (ny - 1), function(i_){
  ts <- sort(which(years == years_unique[i_]), decreasing = T)
  
  Reduce("%*%", 
         lapply(ts, function(t_){AA[,,t_]}))
  
})

aa_array <- aa %>% unlist() %>% array(dim = c(n, n, (ny - 1)))

n_proj <- array(0, dim = c(n, ny - 1))

n_proj[, 1] <- x[, 1]

for (y in 2 : (ny - 1)) {
  
  n_proj[, y] <- aa_array[, , y - 1] %*% n_proj[, y - 1]
  
}

plot(colSums(n_proj), type = "l")

# AA
# 
# 
# 
# surv <- lapply(2:nt, function(t_){
#   return(tibble(t = t_ - 1,
#                 row = rep(1:n, n),
#                 col = rep(1:n, each = n),
#                 val = c(QQ[,,t_ - 1])))
# }) %>%
#   bind_rows() %>%
#   filter(row == col) %>%
#   mutate(basin = factor(ifelse(col < 3, "south", "north"),
#                         levels = c("south","north")),
#          stage = factor(ifelse(col %in% c(1, 3), "juvenile", "adult"),
#                         levels = c("juvenile","adult")))
# 
# 
# surv %>%
#   ggplot(aes(t, val, color = basin))+
#   facet_wrap(~stage, nrow = 2)+
#   geom_line()+
#   scale_color_manual(values = c("dodgerblue","gray20"))+
#   theme_bw()
# 
# 
# disp <- lapply(2:nt, function(t_){
#   return(tibble(t = t_ - 1,
#                 row = rep(1:n, n),
#                 col = rep(1:n, each = n),
#                 val = c(DD[,,t_ - 1])))
# }) %>%
#   bind_rows() %>%
#   filter(row %in% c(2, 4),
#          col %in% c(2, 4),
#          row != col) %>%
#   mutate(basin = factor(ifelse(col < 3, "south", "north"),
#                         levels = c("south","north")))
# 
# 
# disp %>%
#   ggplot(aes(t, val, color = basin))+
#   geom_line()+
#   scale_color_manual(values = c("dodgerblue","gray20"))+
#   theme_bw()+
#   scale_y_continuous(limits = c(0, 1))

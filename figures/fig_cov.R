#==========
#========== Preliminaries
#==========

source("figures/fig_setup.R")
source("analysis/population_projection_functions.R")

options(mc.cores = parallel::detectCores()-6)

# extract setup values
years <- annual_output[[1]]$setup$years
ids <- annual_output[[1]]$setup$ids
years_unique <- unique(years)
ny = length(years_unique)






#==========
#========== Extract and summarize dispersal
#==========

# set theta labels
theta_sub <- matrix(c("r1s","r1w",
                      "r2s","r2w",
                      "s1s","s1w",
                      "s2s","s2w",
                      "s3s","s3w",
                      "s4s","s4w",
                      "g1s","g1w",
                      "d1s","d1w",
                      "d2s","d2w"),
                    nrow = 9,
                    byrow = T)

theta_order <- c("d2s",
                 "d1s",
                 "s4s",
                 "s3s",
                 "s2s",
                 "s1s",
                 "r2s",
                 "r1s")

theta_labs <- c(expression(delta[n%->%s]),
                expression(delta[s%->%n]),
                expression(phi["a,n"]),
                expression(phi["j,n"]),
                expression(phi["a,s"]),
                expression(phi["j,s"]),
                expression(rho[n]),
                expression(rho[s]))


# extract theta_s
theta_full <- parallel::mclapply(ids, function(i_){
  
  xx_ <- annual_output[[i_]]$annual_proj
  
  theta_ <- sapply(1:ny ,function(y){
    ts_ <- which(years == years_unique[y])
    c(xx_$RR[1,2,ts_[1]], xx_$RR[3,4,ts_[1]],
      xx_$QQs[1,1,ts_[1]], xx_$QQs[2,2,ts_[1]], xx_$QQs[3,3,ts_[1]], xx_$QQs[4,4,ts_[1]],
      xx_$GGs[2,1,ts_[1]], xx_$DDs[4,2,ts_[1]], xx_$DDs[2,4,ts_[1]],
      xx_$RR[1,2,ts_[2]], xx_$RR[3,4,ts_[2]],
      xx_$QQs[1,1,ts_[2]], xx_$QQs[2,2,ts_[2]], xx_$QQs[3,3,ts_[2]], xx_$QQs[4,4,ts_[2]],
      xx_$GGs[2,1,ts_[2]], xx_$DDs[4,2,ts_[2]], xx_$DDs[2,4,ts_[2]])
  }) %>% t()
  colnames(theta_) <- theta_names
  
  xx_m <- sapply(1:nrow(theta_sub),
                 function(i_){
                   (sapply(1:ny,
                           function(j_){theta_[j_, theta_sub[i_,]]}))
                 })
  
  out <- as_tibble(xx_m) %>%
    set_names(theta_sub[,1]) %>%
    mutate(time = row_number(),
           i = i_) 
  
  return(out)
}) %>%
  bind_rows()

# summarize theta_s
theta_sum <- theta_full %>%
  gather(var, val, -time, -i) %>%
  group_by(time, var) %>%
  summarize(mi = median(val, na.rm = T)) %>%
  ungroup()

theta_sum %>%
  filter(var %in% c("r1s","r2s")) %>%
  ggplot(aes(time, mi, color = var))+
  geom_line()


# correlation matrix
cor_theta <- theta_sum %>%
  filter(var != "g1s") %>%
  spread(var, mi) %>%
  select(-time) %>%
  cor(method = "spearman") %>%
  as_tibble() %>%
  {mutate(., row = names(.))} %>%
  gather(col, val, -row) %>%
  mutate(row = factor(row, 
                      levels = theta_order),
         col = factor(col, 
                      levels = theta_order),
         col = factor(col),
         x = as.numeric(row),
         y = as.numeric(col))




ggplot(data = cor_theta %>%
         filter(x >= y),
       aes(x = x, 
           y= y))+
  geom_tile(aes(fill = val)) +
  scale_fill_gradient2(name = "Correlation",
                       low ="dodgerblue", 
                       mid="white", 
                       high="firebrick", 
                       midpoint = 0,
                       limits = c(-1, 1),
                       breaks = c(-1, 0, 1))+
  scale_x_reverse("",
                  breaks = 1:length(theta_labs),
                  labels = theta_labs)+
  scale_y_continuous("",
                     breaks = 1:length(theta_labs),
                     labels = theta_labs)+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.5,0.9),)+
  coord_flex_fixed(left= capped_vertical('both'), 
                  bottom = capped_horizontal('both'),
                  ratio = 1)






# correlation matrix
cor_theta <- theta_sum %>%
  filter(var != "g1s") %>%
  group_by(var) %>%
  mutate(mi = (mi - mean(mi))/sd(mi)) %>%
  spread(var, mi) %>%
  select(-time) %>%
  cov() %>%
  as_tibble() %>%
  {mutate(., row = names(.))} %>%
  gather(col, val, -row) %>%
  mutate(row = factor(row, 
                      levels = theta_order),
         col = factor(col, 
                      levels = theta_order),
         col = factor(col),
         x = as.numeric(row),
         y = as.numeric(col))




ggplot(data = cor_theta %>%
         filter(x >= y),
       aes(x = x, 
           y= y))+
  geom_tile(aes(fill = val)) +
  scale_fill_gradient2(name = "Correlation",
                       low ="dodgerblue", 
                       mid="white", 
                       high="firebrick", 
                       midpoint = 0)+
  scale_x_reverse("",
                  breaks = 1:length(theta_labs),
                  labels = theta_labs)+
  scale_y_continuous("",
                     breaks = 1:length(theta_labs),
                     labels = theta_labs)+
  theme(panel.border = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.7, "lines"),
        legend.position = c(0.5,0.9),)+
  coord_flex_fixed(left= capped_vertical('both'), 
                   bottom = capped_horizontal('both'),
                   ratio = 1)









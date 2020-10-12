### define elements of annual projection matrix
a_list <- list(a11 = expression(s1w*(1-g1w)*s1s*(1-g1s)+r1w*s1s*g1s*(1-d1s)),
               a21 = expression(s1w*g1w*(1-d1w)*s1s*(1-g1s)+
                                  s2w*(1-d1w)*s1s*g1s*(1-d1s)+s4w*d2w*s1s*g1s*d1s),
               a31 = expression(r2w*s1s*g1s*d1s),
               a41 = expression(s1w*g1w*d1w*s1s*(1-g1s)+s2w*d1w*s1s*g1s*(1-d1s)+
                                  s4w*(1-d2w)*s1s*g1s*d1s),
               a12 = expression(s1w*(1-g1w)*r1s+r1w*s2s*(1-d1s)),
               a22 = expression(s1w*g1w*(1-d1w)*r1s+s2w*(1-d1w)*s2s*(1-d1s)+
                                  s4w*d2w*s2s*d1s),
               a32 = expression(r2w*s2s*d1s),
               a42 = expression(s1w*g1w*d1w*r1s+s2w*d1w*s2s*(1-d1s)+
                                  s4w*(1-d2w)*s2s*d1s),
               a13 = expression(r1w*s3s*g1s*d2s),
               a23 = expression(s2w*(1-d1w)*s3s*g1s*d2s+s3w*g1w*d2w*s3s*(1-g1s)+
                                  s4w*d2w*s3s*g1s*(1-d2s)),
               a33 = expression(s3w*(1-g1w)*s3s*(1-g1s)+r2w*s3s*g1s*(1-d2s)),
               a43 = expression(s2w*d1w*s3s*g1s*d2s+s3w*g1w*(1-d2w)*s3s*(1-g1s)+
                                  s4w*(1-d2w)*s3s*g1s*(1-d2s)),
               a14 = expression(r1w*s4s*d2s),
               a24 = expression(s2w*(1-d1w)*s4s*d2s+s3w*g1w*d2w*r2s+s4w*d2w*s4s*(1-d2s)),
               a34 = expression(s3w*(1-g1w)*r2s+r2w*s4s*(1-d2s)),
               a44 = expression(s2w*d1w*s4s*d2s+s3w*g1w*(1-d2w)*r2s+
                                 s4w*(1-d2w)*s4s*(1-d2s)))



### define function for evaluating matrix elements 
m_fn <- function(v_){
  
  r1s <- v_["r1s"]
  r2s <- v_["r2s"]
  s1s <- v_["s1s"]
  s2s <- v_["s2s"]
  s3s <- v_["s3s"]
  s4s <- v_["s4s"]
  g1s <- v_["g1s"]
  d1s <- v_["d1s"]
  d2s <- v_["d2s"]
  r1w <- v_["r1w"]
  r2w <- v_["r2w"]
  s1w <- v_["s1w"]
  s2w <- v_["s2w"]
  s3w <- v_["s3w"]
  s4w <- v_["s4w"]
  g1w <- v_["g1w"]
  d1w <- v_["d1w"]
  d2w <- v_["d2w"]
  
  a11 = eval(a_list[[1]])
  a21 = eval(a_list[[2]])
  a31 = eval(a_list[[3]])
  a41 = eval(a_list[[4]])
  
  a12 = eval(a_list[[5]])
  a22 = eval(a_list[[6]])
  a32 = eval(a_list[[7]])
  a42 = eval(a_list[[8]])
  
  a13 = eval(a_list[[9]])
  a23 = eval(a_list[[10]])
  a33 = eval(a_list[[11]])
  a43 = eval(a_list[[12]])
  
  a14 = eval(a_list[[13]])
  a24 = eval(a_list[[14]])
  a34 = eval(a_list[[15]])
  a44 = eval(a_list[[16]])
  
  matrix(c(a11, a12, a13, a14,
           a21, a22, a23, a24,
           a31, a32, a33, a34,
           a41, a42, a43, a44),
         nrow = 4,
         ncol = 4,
         byrow = T)
  
}





### define variables to differential with respect to 
par_names <- c("r1s","r2s",
               "s1s","s2s","s3s","s4s",
               "g1s","d1s","d2s",
               "r1w","r2w",
               "s1w","s2w","s3w","s4w",
               "g1w","d1w","d2w")




### extract parameter values
pars <- sapply(1:(ny - 1) ,function(y){
  ts_ <- which(years == years_unique[y])
  c(RR[1,2,ts_[1]],RR[3,4,ts_[1]],
    QQ[1,1,ts_[1]], QQ[2,2,ts_[1]], QQ[3,3,ts_[1]], QQ[4,4,ts_[1]],
    GG[2,1,ts_[1]], DD[4,2,ts_[1]], DD[2,4,ts_[1]],
    RR[1,2,ts_[2]],RR[3,4,ts_[2]],
    QQ[1,1,ts_[2]], QQ[2,2,ts_[2]], QQ[3,3,ts_[2]], QQ[4,4,ts_[2]],
    GG[2,1,ts_[2]], DD[4,2,ts_[2]], DD[2,4,ts_[2]])
}) %>% t()
colnames(pars) <- par_names





### define function for derivatives
d_list <- lapply(a_list,
                 function(x_){
                   xx_ <- x_
                   deriv(xx_, 
                         par_names,
                         function(r1s,r2s,s1s,s2s,s3s,s4s,g1s,d1s,d2s,
                                  r1w,r2w,s1w,s2w,s3w,s4w,g1w,d1w,d2w){})})





### define function for evaluating derivatives
eval_fn <- function(v_, i_){
  
  r1s <- v_["r1s"]
  r2s <- v_["r2s"]
  s1s <- v_["s1s"]
  s2s <- v_["s2s"]
  s3s <- v_["s3s"]
  s4s <- v_["s4s"]
  g1s <- v_["g1s"]
  d1s <- v_["d1s"]
  d2s <- v_["d2s"]
  r1w <- v_["r1w"]
  r2w <- v_["r2w"]
  s1w <- v_["s1w"]
  s2w <- v_["s2w"]
  s3w <- v_["s3w"]
  s4w <- v_["s4w"]
  g1w <- v_["g1w"]
  d1w <- v_["d1w"]
  d2w <- v_["d2w"]
  
  d_list[[i_]](r1s,r2s,s1s,s2s,s3s,s4s,g1s,d1s,d2s,
              r1w,r2w,s1w,s2w,s3w,s4w,g1w,d1w,d2w)
}




#### 
np <- ncol(pars)

Is_ <- diag(nrow = n, ncol = n) 

# initial state distribution
n_ <- x0

# initial derivatives of population size with respect to matrix elemetns
# set to 0
dn_ <- matrix(0, nrow = n, ncol = n^2) 
dnp_ <- matrix(0, nrow = n, ncol = np) 

# initialize lists for storing results
n_out_ <- list()
dn_out_ <- list()
dnp_out_ <- list()

ip_ <- 100


# loop over time
for (t_ in 1:(ny - 1)) {
  
  # Store results
  n_out_[[t_]] = n_
  dn_out_[[t_]] = dn_
  dnp_out_[[t_]] = dnp_ 
  
  dadp_ <- matrix(sapply(1:(n * n), function(x_){attributes(eval_fn(pars[t_,], x_))$gradient}),
                 ncol = ncol(pars),
                 byrow = T)
  
  # loop over projection intervals
  for (pt_ in 1:ip_) {
    
    # project derivatives
    dn_ = aa_array[, , t_] %*% dn_ + kronecker(t(n_), Is_) 
    
    # project derivatives
    dnp_ = aa_array[, , t_] %*% dnp_ + kronecker(t(n_), Is_) %*% dadp_ %*% diag(pars[t_,])
    
    # Project abundances
    n_ = aa_array[, , t_] %*% n_
  }
  
}
           

plot(colSums(n_proj), type = "l")
plot(sapply(n_out_, sum), type = "l")



# vector for summing across states
c_ <- rep(1, n)



# Loop over time
l_out_ <- lapply(2:(ny - 1), function(x_){
  
  # proportional rate of increase (lambda)
  l_ = (sum(n_out_[[x_]])/sum(n_out_[[x_ - 1]]))^(1/ip_) # realized 
  l_asym_ = c(demogR::eigen.analysis(aa_array[, , x_ - 1])$lambda) # asymptotic
  
  # intrinsic growth rate
  r_ = log(l_) # realized 
  r_asym_ = log(l_asym_) # asymptotic
  
  
  # sensitivity of realized r and lambda 
  r_sens_ = t(((t(c_) %*% dn_out_[[x_]])/sum(n_out_[[x_]]) - 
                 (t(c_) %*% dn_out_[[x_ - 1]])/sum(n_out_[[x_ - 1]])))
  r_sens_ = r_sens_/ip_
  l_sens_ = l_ * r_sens_
  
  # sensitivity of asymptotic r and lambda 
  l_sens_asym_ = c(demogR::eigen.analysis(aa_array[, , x_ - 1])$sensitivities)
  r_sens_asym_ =  l_sens_asym_/l_asym_
  
  
  # sensitivity of realized r and lambda with respect to pars
  l_elas_p_ = t(((t(c_) %*% dnp_out_[[x_]])/sum(n_out_[[x_]]) - 
                 (t(c_) %*% dnp_out_[[x_ - 1]])/sum(n_out_[[x_ - 1]])))
  l_elas_p_ = l_elas_p_ / ip_
  r_elas_p_ = l_elas_p_ / l_

  
  # dataframe for export
  list(time = x_ - 1,
         row = rep(1:n, n),
         col = rep(1:n, each = n),
         l = l_,
         l_asym = l_asym_,
         r = r_,
         r_asym = r_asym_,
         l_sens = l_sens_,
         l_sens_asym = l_sens_asym_,
         r_sens = r_sens_,
         r_sens_asym = r_sens_asym_,
         l_elas_p = l_elas_p_)
})

l_elas_p <- sapply(l_out_, function(x_){
  x_$l_elas_p 
})

l_elas_p_clean <- as_tibble(t(l_elas_p)) %>%
  set_names(par_names) %>%
  mutate(time = row_number()) %>%
  gather(var, val, -time) %>%
  mutate(name = str_split(var, "") %>% map_chr(~as.character(.x[[1]])),
         index = str_split(var, "") %>% map_chr(~as.character(.x[[2]])),
         name = paste0(name, index),
         season = str_split(var, "") %>% map_chr(~as.character(.x[[3]]))) %>%
  group_by(name) %>%
  arrange(name, time, season) %>%
  mutate(season = ifelse(season == "s", 0, 1),
         new_time = time + cumsum(season))


l_elas_p_clean %>%
  ggplot(aes(new_time, val))+
  facet_wrap(~name)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(size = 1)+
  theme_classic()

l_elas_p_clean %>%
  ggplot(aes(name, val, color = new_time))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_jitter(size = 2, height = 0, width = 0.1)+
  theme_classic()


l_elas_p_season <- l_elas_p_clean %>%
  group_by(name, time) %>%
  summarize(val = sum(val))

l_elas_p_season %>%
  ggplot(aes(time, val))+
  facet_wrap(~name)+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(size = 1)+
  theme_classic()

l_elas_p_season %>%
  ggplot(aes(name, val, color = time))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_violin(linetype = 0, fill = "gray70", width = 1.5)+
  geom_jitter(size = 2, height = 0, width = 0.1)+
  theme_classic()




l_elas_mean <- l_elas_p_season %>%
    group_by(name) %>%
    summarize(val = mean(val)) %>%
    arrange(-val)

l_elas_p_season %>%
  ungroup() %>%
  mutate(name = factor(name, 
                       levels = l_elas_mean$name),
         pos = as.numeric(name)) %>%
  ggplot(aes(pos, val))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_jitter(aes(fill = time),
              size = 2, height = 0, width = 0.1,
              shape = 21, stroke = 0.5)+
  geom_segment(data = l_elas_mean %>% 
                 mutate(name = factor(name, 
                                      levels = l_elas_mean$name),
                        pos = as.numeric(name),
                        xmin = pos - 0.2,
                        xmax = pos + 0.2),
               aes(x = xmin,
                   y = val,
                   xend = xmax,
                   yend = val),
               size = 2)+
  scale_x_continuous("Demographic rate",
                     breaks = 1:9,
                     labels = l_elas_mean$name)+
  scale_y_continuous(Proportional~sensitivity~of~lambda)+
  scale_fill_gradient2(low ="dodgerblue", 
                       mid="gray90", 
                       high="firebrick", 
                       midpoint = 15)+
  coord_flip()+
  theme_classic()






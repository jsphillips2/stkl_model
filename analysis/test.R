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


par_names <- c("r1s","r2s",
               "s1s","s2s","s3s","s4s",
               "g1s","d1s","d2s",
               "r1w","r2w",
               "s1w","s2w","s3w","s4w",
               "g1w","d1w","d2w")


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


d_list <- lapply(a_list,
                 function(x_){
                   xx_ <- x_
                   deriv(xx_, 
                         par_names,
                         function(r1s,r2s,s1s,s2s,s3s,s4s,g1s,d1s,d2s,
                                  r1w,r2w,s1w,s2w,s3w,s4w,g1w,d1w,d2w){})})

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

ip_ <- 10


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
    dn_ = aa_array[, , t_] %*% dn_ + kronecker(t(n_), Is_) %*% diag(c(aa_array[, , t_]))
    
    
    # project derivatives
    dnp_ = aa_array[, , t_] %*% dnp_ + kronecker(t(n_), Is_) %*% dadp_ %*% diag(pars[t_,])
    
    # Project abundances
    n_ = aa_array[, , t_] %*% n_
  }
  
}


# vector for summing across states
c_ <- rep(1, n)



# Loop over time
l_out_2 <- lapply(2:(ny - 1), function(x_){
  
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
  r_sens_p_ = t(((t(c_) %*% dnp_out_[[x_]])/sum(n_out_[[x_]]) - 
                   (t(c_) %*% dnp_out_[[x_ - 1]])/sum(n_out_[[x_ - 1]])))
  r_sens_p_ = r_sens_p_/ip_
  l_sens_p_ = l_ * r_sens_p_
  
  # elasticity
  l_elas_asym_ = c(demogR::eigen.analysis(aa_array[, , x_ - 1])$elasticities)
  r_elas_asym_ =  l_elas_asym_ / l_asym_
  l_elas_ = t(((t(c_) %*% dn_out_[[x_]])/sum(n_out_[[x_]]) - 
                 (t(c_) %*% dn_out_[[x_ - 1]])/sum(n_out_[[x_ - 1]])))
  l_elas_ = l_elas_/ip_
  r_elas_ = l_elas_ / l_
  
  
  # dataframe for export
  out = list(time = x_ - 1,
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
       r_sens_p = r_sens_p_,
       l_sens_p = l_sens_p_,
       l_elas_asym = l_elas_asym_,
       r_elas_asym = r_elas_asym_,
       l_elas = l_elas_,
       r_elas = r_elas_
  )
  
  return(out)
})


l_elas_asym_out <- sapply(l_out_2, function(x_){
  x_$l_elas_asym 
})

r_elas_asym_out <- sapply(l_out_2, function(x_){
  x_$r_elas_asym 
})

l_elas_out <- sapply(l_out_2, function(x_){
  x_$l_elas
})

r_elas_out <- sapply(l_out_2, function(x_){
  x_$r_elas
})

plot(l_elas_asym_out[1, ], col = "blue", type = "l")
lines(l_elas_out[1, ], col = "red")

plot(r_elas_asym_out[1, ], col = "blue", type = "l")
lines(r_elas_out[1, ], col = "red")

colSums(r_elas_asym_out)
colSums(r_elas_out)

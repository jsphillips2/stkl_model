
source("analysis/matrix_functions.R")


#========================================================================================= 
# Seasonal and annual population projection
#=========================================================================================

proj_fn <- function(extract_, pj_, nt_, n_, b_, years_, theta_names_, nod_) {
  
  ### extract parameters
  x0 <- extract_ %>%
    filter(name == "x0") %>%
    {array(.$val, dim = n)}
  
  
  rt <- extract_ %>%
    filter(name == "rt") %>%
    {array(.$val, dim = c(nt - 1, max(pj_[, 1, 1]) - 1))} %>%
    t()
  
  rf <- extract_ %>%
    filter(name == "rf") %>%
    {array(.$val, dim = c(max(pj_[, 1 , 2]) - 1))}
  
  qt <- extract_ %>%
    filter(name == "qt") %>%
    {array(.$val, dim = c(nt_ - 1, max(pj_[, 1 ,3]) - 1))} %>%
    t()
  
  qf <- extract_ %>%
    filter(name == "qf") %>%
    {array(.$val, dim = c(max(pj_[, 1 , 4]) - 1))}
  
  gt <- extract_ %>%
    filter(name == "gt") %>%
    {array(.$val, dim = c(nt_ - 1, max(pj_[, 1 ,5]) - 1))} %>%
    t()
  
  gf <- extract_ %>%
    filter(name == "gf") %>%
    {array(.$val, dim = c(max(pj_[, 1 , 6]) - 1))}
  
  dt <- extract_ %>%
    filter(name == "dt") %>%
    {array(.$val, dim = c(nt_ - 1, max(pj_[, 1 ,7]) - 1))} %>%
    t()
  
  df <- extract_ %>%
    filter(name == "df") %>%
    {array(.$val, dim = c(max(pj_[, 1 , 8]) - 1))}
  
  
  
  
  
  ### define matrices for filling  
  Rt <- matrix(0, nrow = n_, ncol = n_)
  Rf <- matrix(0, nrow = n_, ncol = n_)
  Qt <- matrix(0, nrow = n_, ncol = n_)
  Qf <- matrix(0, nrow = n_, ncol = n_)
  Gt <- matrix(0, nrow = n_, ncol = n_)
  Gf <- matrix(0, nrow = n_, ncol = n_)
  Dt <- matrix(0, nrow = n_, ncol = n_)
  Df <- matrix(0, nrow = n_, ncol = n_)
  no_D <- matrix(0, nrow = n_, ncol = n_)
  
  ### define arrays for filling 
  RR <- array(0, dim = c(n_, n_, nt_ - 1))
  QQ <- array(0, dim = c(n_, n_, nt_ - 1))
  GG <- array(0, dim = c(n_, n_, nt_ - 1))
  DD <- array(0, dim = c(n_, n_, nt_ - 1))
  AA <- array(0, dim = c(n_, n_, nt_ - 1))
  x <- array(0, dim = c(n_, nt_))
  noDD <- array(0, dim = c(n_, n_, nt_ - 1))   # matrix for no dispersal
  AA_noD <- array(0, dim = c(n_, n_, nt_ - 1)) # matrix for no dispersal
  x_noD <- array(0, dim = c(n_, nt_))          # matrix for no dispersal
  QQs <- array(0, dim = c(n_, n_, nt_ - 1)) # standardized time steps
  GGs <- array(0, dim = c(n_, n_, nt_ - 1)) # standardized time steps
  DDs <- array(0, dim = c(n_, n_, nt_ - 1)) # standardized time steps
  
  
  
  #### define function transition probability matrix
  trans_fn <- function(m0_){
    m1_ = m0_
    diag(m1_) <- 0
    v_ <- colSums(m0_)
    m2_ <- diag(v_)
    m3_ <- m1_ - m2_
    m4_ <- expm(m3_)
    return(as.matrix(m4_))
  }
  
  
  
  
  ### set initial population size
  x[,1] <- x0
  x_noD[,1] <- x0
  
  
  
  
  
  #### project over time
  for (t in 2:nt_) {
    for (i in 1 : (n_ * n_)) {
      
      # fill matrices
      Rt[pj_[i, 2, 1], pj_[i, 3, 1]] <- c(0, rt[, t - 1])[pj_[i, 1, 1]]
      Rf[pj_[i, 2, 2], pj_[i, 3, 2]] <- c(0, rf)[pj_[i, 1, 2]]
      Qt[pj_[i, 2, 3], pj_[i, 3, 3]] <- c(0, qt[, t - 1])[pj_[i, 1, 3]]
      Qf[pj_[i, 2, 2], pj_[i, 3, 4]] <- c(0, qf)[pj_[i, 1, 4]]
      Gt[pj_[i, 2, 5], pj_[i, 3, 5]] <- c(0, gt[, t - 1])[pj_[i, 1, 5]]
      Gf[pj_[i, 2, 6], pj_[i, 3, 6]] <- c(0, gf)[pj_[i, 1, 6]]
      Dt[pj_[i, 2, 7], pj_[i, 3, 7]] <- c(0, dt[, t - 1])[pj_[i, 1, 7]]
      Df[pj_[i, 2, 8], pj_[i, 3, 8]] <- c(0, df)[pj_[i, 1, 8]]
      
      # calculate projection matricees
      R <- Rt + Rf
      Q <- trans_fn(b_[t - 1] * (Qt + Qf))
      G <- trans_fn(b_[t - 1] * (Gt + Gf))
      D <- trans_fn(b_[t - 1] * (Dt + Df))
      A <-  R + (D %*% G %*% Q)
      diag(no_D) <- nod_ # no dispersal (retain loss on diagonal)        
      A_noD <- R + (no_D %*% G %*% Q) # no dispersal (retain loss on diagonal)   
      
      # project dynamics
      x[, t] <- A %*% x[, t - 1]
      x_noD[, t] <- A_noD %*% x_noD[, t - 1]
      
      # store projection matrices for export
      RR[, , t - 1] <- R
      QQ[, , t - 1] <- Q
      GG[, , t - 1] <- G
      DD[, , t - 1] <- D
      AA[, , t - 1] <- A
      noDD[, , t - 1] <- no_D
      AA_noD[, , t - 1] <- A_noD
      
      # store projection matrices with standardized time steps
      QQs[, , t - 1] <- trans_fn(mean(b_) * (Qt + Qf))
      GGs[, , t - 1] <- trans_fn(mean(b_) * (Gt + Gf))
      DDs[, , t - 1] <- trans_fn(mean(b_) * (Dt + Df))
      
    }
  }
  
  
  
  
  
  # define unique years
  years_unique <- unique(years_)
  ny <- length(years_unique)
  
  
  
  
  
  ### annual projection matrix
  aa <- lapply(1 : ny, function(i_){
    ts <- sort(which(years_ == years_unique[i_]), decreasing = T)
    
    Reduce("%*%", 
           lapply(ts, function(t_){AA[,,t_]}))
    
  })
  aa_array <- aa %>% unlist() %>% array(dim = c(n_, n_, ny))
  
  
  
  
  ### annual projection matrix (no dispersal)
  aa_noD <- lapply(1 : ny, function(i_){
    ts <- sort(which(years_ == years_unique[i_]), decreasing = T)
    
    Reduce("%*%", 
           lapply(ts, function(t_){AA_noD[,,t_]}))
    
  })
  aa_noD_array <- aa_noD %>% unlist() %>% array(dim = c(n_, n_, ny))
  
  
  
  
  #### project annual population size 
  X_proj <- array(0, dim = c(n_, ny))
  X_proj[, 1] <- x[, 1]
  for (y in 2 : ny) {
    
    X_proj[, y] <- aa_array[, , y - 1] %*% X_proj[, y - 1]
    
  }
  
  #### project annual population size (no dispersal)
  X_noD_proj <- array(0, dim = c(n_, ny))
  X_noD_proj[, 1] <- x[, 1]
  for (y in 2 : ny) {
    
    X_noD_proj[, y] <- aa_noD_array[, , y - 1] %*% X_noD_proj[, y - 1]
    
  }
  
  
  
  
  
  ### extract parameter values
  theta <- sapply(1:ny ,function(y){
    ts_ <- which(years == years_unique[y])
    c(RR[1,2,ts_[1]],RR[3,4,ts_[1]],
      QQ[1,1,ts_[1]], QQ[2,2,ts_[1]], QQ[3,3,ts_[1]], QQ[4,4,ts_[1]],
      GG[2,1,ts_[1]], DD[4,2,ts_[1]], DD[2,4,ts_[1]],
      RR[1,2,ts_[2]], RR[3,4,ts_[2]],
      QQ[1,1,ts_[2]], QQ[2,2,ts_[2]], QQ[3,3,ts_[2]], QQ[4,4,ts_[2]],
      GG[2,1,ts_[2]], DD[4,2,ts_[2]], DD[2,4,ts_[2]])
  }) %>% t()
  colnames(theta) <- theta_names_
  
  return(list(RR = RR,
              QQ = QQ,
              GG = GG,
              DD = DD,
              AA = AA,
              noDD = noDD,
              AA_noD = AA_noD,
              QQs = QQs,
              GGs = GGs,
              DDs = DDs,
              x = x,
              x_noD = x_noD,
              aa = aa,
              aa_noD = aa_noD,
              aa_array = aa_array,
              aa_noD_array = aa_noD_array,
              X_proj = X_proj,
              X_noD_proj = X_noD_proj,
              theta = theta,
              n = n_, 
              years = years_))

}





#========================================================================================= 
# Derivatives in total abundance relative to matrix elements and parameters (theta) 
#=========================================================================================

dX_fn <- function(annual_proj_, ip_) {
  
  # define values from annual_proj
  n = annual_proj_$n
  years = annual_proj_$years
  aa_array = annual_proj_$aa_array
  X_proj = annual_proj_$X_proj
  theta = annual_proj_$theta
  
  # define unique years
  years_unique <- unique(years)
  ny <- length(years_unique)
  
  # number of parameters
  np <- ncol(theta)
  
  # identity matrix  
  Is <- diag(nrow = n, ncol = n) 
  
  # initial state distribution
  X <- X_proj[, 1]
  X_prop <- X
  
  # initial derivatives of population size with respect to matrix elemetns
  # set to 0
  dXdA <- matrix(0, nrow = n, ncol = n^2) 
  dXdT <- matrix(0, nrow = n, ncol = np) 
  
  # initialize lists for storing results
  dXdA_out <- list()
  dXdT_out <- list()
  X_prop_out <- list()
  X_out <- list()
  
  # loop over time
  for (y in 1:ny) {
    
    # Store results
    dXdA_out[[y]] = dXdA
    dXdT_out[[y]] = dXdT 
    X_prop_out[[y]] = X_prop
    X_out[[y]] = X
    
    
    dAdT <- matrix(sapply(1:(n* n), 
                          function(j_){attributes(a_deriv_eval_fn(theta[y,], 
                                                                  j_))$gradient}),
                   ncol = ncol(theta),
                   byrow = T)
    
    # loop over projection intervals
    for (ipy in 1:ip_) {
      
      # project derivatives
      dXdA = aa_array[, , y] %*% dXdA + kronecker(t(X_prop), Is) 
      
      # project derivatives
      dXdT = aa_array[, , y] %*% dXdT + kronecker(t(X_prop), Is) %*% dAdT %*% diag(theta[y,])
      
      # Project abundances
      X_prop = aa_array[, , y] %*% X_prop
    }
    
    X = aa_array[, , y] %*% X
    
  }
  
  return(list(n = n,
              years = years,
              ip = ip_,
              aa_array = aa_array,
              dXdA_out = dXdA_out,
              dXdT_out = dXdT_out,
              X_prop_out = X_prop_out,
              X_out = X_out))
  
}





#========================================================================================= 
# Sensitivity analysis
#=========================================================================================

sens_fn <- function(dX_) {
  
  # define values from dX
  n = dX_$n
  years = dX_$years
  aa_array = dX_$aa_array
  ip = dX_$ip
  dXdA_out = dX_$dXdA_out
  dXdT_out = dX_$dXdT_out
  X_prop_out = dX_$X_prop_out
  X_out = dX_$X_out
  
  # define unique years
  years_unique <- unique(years)
  ny <- length(years_unique)
  
  # vector for summing across states
  c <- rep(1, n)
  
  # Loop over time
  l_out <- lapply(2:ny, function(y){
    
    # proportional rate of increase (lambda)
    l = (sum(X_prop_out[[y]])/sum(X_prop_out[[y - 1]]))^(1/ip) # realized 
    l_asym = c(demogR::eigen.analysis(aa_array[, , y - 1])$lambda) # asymptotic
    
    # intrinsic growth rate
    r = log(l) # realized 
    r_asym = log(l_asym) # asymptotic
    
    
    # sensitivity of realized r and lambda 
    r_sens = t(((t(c) %*% dXdA_out[[y]])/sum(X_prop_out[[y]]) - 
                   (t(c) %*% dXdA_out[[y - 1]])/sum(X_prop_out[[y - 1]])))
    r_sens = r_sens/ip
    l_sens = l * r_sens
    
    # sensitivity of asymptotic r and lambda 
    l_sens_asym = c(demogR::eigen.analysis(aa_array[, , y - 1])$sensitivities)
    r_sens_asym =  l_sens_asym/l_asym
    
    
    # sensitivity of realized r and lambda with respect to pars
    l_elas_p = t(((t(c) %*% dXdT_out[[y]])/sum(X_prop_out[[y]]) - 
                     (t(c) %*% dXdT_out[[y - 1]])/sum(X_prop_out[[y - 1]])))
    l_elas_p = l_elas_p / ip
    r_elas_p = l_elas_p / l
    
    
    # dataframe for export
    list(l = l,
         l_asym = l_asym,
         r = r,
         r_asym = r_asym,
         l_sens = l_sens,
         l_sens_asym = l_sens_asym,
         r_sens = r_sens,
         r_sens_asym = r_sens_asym,
         l_elas_p = l_elas_p)
  })
  
  names(l_out) = years_unique[1 : (ny - 1)]
  return(l_out)
}



#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(Matrix)
library(matrixcalc)

options(mc.cores = parallel::detectCores()-4)


# import model
out_in <- readRDS(paste0("output/fit_full.rds"))


#==========
#========== Prepare data
#==========

# extract data
data_prep <- out_in$data_prep %>%
  mutate(stage = factor(stage, 
                        levels = c("small","large"),
                        labels = c("juvenile","adult")))

# extract fit
fit <- out_in$fit

# summarize fit
fit_summary <- rstan::summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_tibble(.) %>%
    mutate(var = rownames(rstan::summary(fit)$summary))}

# MCMC specs
iter <- out_in$mcmc_specs$iter
chains <- out_in$mcmc_specs$chains

# data passed to stan
data_list <- out_in$data_list

# data frame for matching
state_match <- data_prep %>% 
  tidyr::expand(nesting(state, basin, stage))  %>%
  arrange(basin, stage) %>%
  mutate(st = row_number(),
         state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "adult.south",
                                   "juvenile.north",
                                   "adult.north"),
                        labels = c("juvenile\nsouth",
                                   "adult\nsouth",
                                   "juvenile\nnorth",
                                   "adult\nnorth")))

# year match
date_match <- data_prep$date %>% unique()

# extract projection matrix
A_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "AA"))}$var
A <-  rstan::extract(fit, pars = A_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(A_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  sample_n(4000) %>%
  gather(var, value, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4])),
         date = date_match[time],
         year = lubridate::year(date),
         stage = rep(c("juvenile","adult"), 2)[col] %>%
           factor(levels = c("juvenile","adult")),
         basin = rep(c("south","north"), each = 2)[col] %>%
           factor(levels = c("south","north"))) %>%
  arrange(time, col, row)

# extract population size vector
n0_vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "x0"))}$var
n0 <- rstan::extract(fit, pars = n0_vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(n0_vars) %>%
  mutate(chain = rep(1:out_in$mcmc_specs$chains, each = iter/2), 
         step = rep(c(1:(out_in$mcmc_specs$iter/2)), chains),
         chain_step = paste(chain,step, sep="_")) %>%
  filter(chain_step %in% unique(A$chain_step)) %>%
  gather(var, value, -chain, -step, -chain_step) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  full_join(state_match) 





#==========
#========== Functions for projection analysis
#==========

##### Function for annual projection matrix
aa_fun <- function(x_) {
  
  x_ %>%
    
    # extract relevant variables
    select(year, date, row, col, value) %>% 
    
    # split by year
    split(.$year) %>%
    
    # apply matrix multiplication to all matrices within a given year
    lapply(function(x2_){
      Reduce("%*%",
             
             # split by date
             x2_ %>% 
               split(.$date) %>%
               lapply(function(x3_){
                 
                 # arrange by columns, then rows
                 x3_ = x3_ %>% arrange(col, row)
                 
                 # fill matrix
                 matrix(x3_$value, nrow = max(x3_$row), ncol = max(x3_$col), byrow = F)
          }))
    })
}



##### Function for annual projections and associated analyses
# Following Caswell 2007; Koons et al 2016; Koons et al. 2017
# Default to one  projection interval per annual matrix; can use to check again asympt.
proj_fun <- function(aa_, n0_, ip_ = 1) {
  
  
  
  ##### Project abundance & associated derivatives for sensitivty analysis
  
  # extract matrix dimensions
  s_ <- nrow(aa_[[1]])
  
  # number of time steps
  tmax_ <- length(aa_)
  
  # convert list of projection matrices to array
  aa_array_ <- aa_ %>% unlist() %>% array(dim = c(s_, s_, tmax_))
  
  # indentity matrix
  Is_ <- diag(nrow = s_, ncol = s_) 
  
  # initial state distribution
  n_ <- n0_ 
  
  # initial derivatives of population size with respect to matrix elemetns
  # set to 0
  dn_ <- matrix(0, nrow = s_, ncol = s_^2) 
  
  # initialize lists for storing results
  n_out_ <- list()
  dn_out_ <- list()
  
  # loop over time
  for (t_ in 1:(tmax_)) {
    
    # Store results
    n_out_[[t_]] = n_
    dn_out_[[t_]] = dn_
    
    # loop over projection intervals
    for (pt_ in 1:ip_) {
      
      # project derivatives
      dn_ = aa_array_[, , t_] %*% dn_ + kronecker(t(n_), Is_) 
      
      # Project abundances
      n_ = aa_array_[, , t_] %*% n_
    }
    
  }
  
  
  
  ##### Caculate annual state contributions
  aa_cont_ <- lapply(1:(tmax_ - 1), function(x_){
    hadamard.prod(aa_[[x_]], 
                  matrix(rep(n_out_[[x_]],
                             each = s_),
                         nrow=s_))
    }
  )
  
  
  
  ##### Population growth rate & sensitivty analysis
  
  # vector for summing across states
  c_ <- rep(1, s_)
  
  # Loop over time
  l_out_ <- lapply(2:tmax_, function(x_){
    
    # proportional rate of increase (lambda)
    l_ = (sum(n_out_[[x_]])/sum(n_out_[[x_ - 1]]))^(1/ip_) # realized 
    l_asym_ = c(demogR::eigen.analysis(aa_array_[, , x_ - 1])$lambda) # asymptotic
    
    # intrinsic growth rate
    r_ = log(l_) # realized 
    r_asym_ = log(l_asym_) # asymptotic
    
    
    # sensitivity of realized r and lambda 
    r_sens_ = t(((t(c_) %*% dn_out_[[x_]])/sum(n_out_[[x_]]) - 
                  (t(c_) %*% dn_out_[[x_ - 1]])/sum(n_out_[[x_ - 1]])))
    r_sens_ = r_sens_/ip_
    l_sens_ = l_ * r_sens_
    
    # sensitivity of asymptotic r and lambda 
    l_sens_asym_ = c(demogR::eigen.analysis(aa_array_[, , x_ - 1])$sensitivities)
    r_sens_asym_ =  l_sens_asym_/l_asym_
    
    # dataframe for export
    tibble(time = x_ - 1,
           row = rep(1:s_, s_),
           col = rep(1:s_, each = s_),
           l = l_,
           l_asym = l_asym_,
           r = r_,
           r_asym = r_asym_,
           l_sens = l_sens_,
           l_sens_asym = l_sens_asym_,
           r_sens = r_sens_,
           r_sens_asym = r_sens_asym_
    )
  }) %>%
    bind_rows()
  
  
  
  ##### Variance partitioning of population growth rate
  
  # convert list of stage distributions to array
  n_array_ <- n_out_ %>% 
    unlist() %>% 
    array(dim = c(s_, tmax_))
  
  # normalize stage distributions by magnitude
  n_array_norm_ <- n_array_ %>%
    apply(2, function(x_){x_/sum(x_)})
  
  # sensitivity of projection matrix elements
  a_sens_ <- rep(rowSums(n_array_norm_)/tmax_, each = s_)
  
  # sensitivity of abundances
  n_sens_ <- (aa_ %>% lapply(colSums) %>% 
               unlist() %>% 
               array(dim = c(s_, tmax_)) %>% 
               rowSums)/tmax_
  
  # combine sensitivities
  sens_ <- c(a_sens_, n_sens_)
  
  # time x parameter matrix
  pars_ <- cbind(t(aa_ %>% 
                     sapply(function(x_){t(c(x_))})), 
                 t(n_array_norm_))
  
  # variance-covariance matrix of parameters
  covar_ <- cov(pars_)
  
  # contribution of parameter combinations to variance in lambda
  cont_mat_ <- matrix(0, nrow = nrow(covar_), ncol = ncol(covar_))
  for (i in 1:nrow(covar_)) {
    for(j in 1:ncol(covar_)) {
      cont_mat_[i, j] = covar_[i, j] * sens_[i] * sens_[j]
    }
  }
  
  # total contribution of each parameter to variance in lambda
  cont_ <- rowSums(cont_mat_)
  
  # approximated variance in lambda
  total_cont_ <- sum(cont_)
  
  # error and log error in approximation
  approx_err_ <- var(l_out_$l)/total_cont_

  # deviations in parameters
  pars_dev_ <- pars_ %>% apply(2, diff) 
  
  # contribution of parameter deviation to change in lambda
  pars_cont <- pars_dev_ %>% 
    apply(1, function(x){x * sens_}) %>% 
    matrix(nrow = ncol(pars_), ncol = nrow(pars_) - 1) %>% t()
  
  # bundle variance partitioning 
  var_part_ <- list(sens = sens_,
                   covar = covar_,
                   cont_mat = cont_mat_,
                   cont = cont_,
                   approx_err = approx_err_,
                   pars_cont = pars_cont)
  
  # return output
  return(list(n_array_ = n_array_, 
              aa_cont_array_ = aa_cont_,
              l_out = l_out_,
              var_part = var_part_
              ))
  
}





#==========
#========== Projection analysis
#==========

###### Apply projection analysis to posterior draws
# extract iterations
chain_steps <- unique(A$chain_step)
proj_full <- parallel::mclapply(chain_steps, function(x_){
  aa_ <- aa_fun(A %>% filter(chain_step == x_))
  n0_ <- {n0 %>% filter(chain_step == x_) %>% arrange(st)}$value
  proj_ <- proj_fun(aa_, n0_)
  return(list(aa = aa_, proj = proj_))
})
names(proj_full) = chain_steps



##### Summarize abundance projection
n_sum <- parallel::mclapply(chain_steps, function(x_){
  x2_ <- proj_full[[x_]]
  n_ <- t(x2_$proj$n_array_)
  colnames(n_) = 1:4
  n_ <- as_tibble(n_) %>%
    mutate(time = row_number()) %>%
    gather(st, value, -time) %>%
    mutate(iter = x_,
           st = as.numeric(factor(st, levels = c(1:4))))
}) %>%
  bind_rows() %>%
  group_by(time, st) %>%
  summarize(value_lo = quantile(value, probs = c(0.16)),
            value_hi = quantile(value, probs = c(0.84)),
            value = median(value)) %>%
  full_join(state_match) %>%
  ungroup()


##### Summarize per capita contributions
aa_sum <- parallel::mclapply(chain_steps, function(x_){
  x2_ <- proj_full[[x_]]
  aa <- x2_$aa
  lapply(1:length(aa), function(x3_) {
      tibble(iter = x_,
             time = x3_,
             row = rep(1:nrow(aa[[1]]), nrow(aa[[1]])),
             col = rep(1:nrow(aa[[1]]), each = nrow(aa[[1]])),
             value = c(aa[[x3_]]))
    }) %>% bind_rows()
}) %>% bind_rows()  %>%
  group_by(time, row, col) %>%
  summarize(value_lo = quantile(value, probs = c(0.16)),
            value_hi = quantile(value, probs = c(0.84)),
            value = median(value)) %>%
  ungroup()



##### Summarize total contributions
aa_cont_sum <- parallel::mclapply(chain_steps, function(x_){
  x2_ <- proj_full[[x_]]
  aa <- x2_$proj$aa_cont_
  lapply(1:length(aa), function(x3_) {
    tibble(iter = x_,
           time = x3_,
           row = rep(1:nrow(aa[[1]]), nrow(aa[[1]])),
           col = rep(1:nrow(aa[[1]]), each = nrow(aa[[1]])),
           value = c(aa[[x3_]]))
  }) %>% bind_rows()
}) %>% bind_rows()  %>%
  group_by(time, row, col) %>%
  summarize(value_lo = quantile(value, probs = c(0.16)),
            value_hi = quantile(value, probs = c(0.84)),
            value = median(value)) %>%
  ungroup()




##### Summarize sensitivity analysis
sens_sum <- parallel::mclapply(chain_steps, function(x_){
  x2_ <- proj_full[[x_]]$proj$l_out  %>% # need to update lout
    mutate(iter = x_)
}) %>%
  bind_rows() %>%
  gather(var, value, -time, -row, -col, - iter) %>%
  group_by(var, time, row, col) %>%
  summarize(value_lo = quantile(value, probs = c(0.16)),
            value_hi = quantile(value, probs = c(0.84)),
            value = median(value)) %>%
  ungroup() %>%
  mutate(basin = ifelse(col < 3, "south", "north") %>%
           factor(levels = c("south","north")),
         stage = ifelse(col %in% c(1, 3), "juvenile","adult") %>%
           factor(levels = c("juvenile","adult")),
         recipient = ifelse(row %in% c(1, 3), "juvenile", "adult") %>%
           factor(levels = c("juvenile","adult")),
         donor = ifelse(col %in% c(1, 3), "juvenile", "adult") %>%
           factor(levels = c("juvenile","adult")),
         type = ifelse( (col %in% c(1,2) & row %in% c(3,4)) | 
                          (col %in% c(3,4) & row %in% c(1,2)),
                        "between", "within"),
         group = interaction(basin, type) %>%
           factor(levels = c("south.within",
                             "north.between",
                             "north.within",
                             "south.between"),
                  labels = c("south \nto south",
                             "north \nto south",
                             "north \nto north",
                             "south \nto north")))



##### Summarize variance partitioning
var_sum_collapsed <- parallel::mclapply(chain_steps, function(x_){
  x2_ <- proj_full[[x_]]$proj$var_part  
  x3_ <- tibble(iter = x_,
                var = names(unlist(x2_)),
                value = unlist(x2_))
})  %>%
  bind_rows() %>%
  mutate(order = row_number()) %>%
  group_by(var) %>%
  summarize(order = min(order),
            value_lo = quantile(value, probs = c(0.16)),
            value_hi = quantile(value, probs = c(0.84)),
            value = median(value)) %>%
  ungroup() %>%
  arrange(order)
var_sum <- list(sens = var_sum_collapsed %>%
                  filter(str_detect(var_sum_collapsed$var, "sens")) %>%
                  mutate(id = str_split(var, "sens") %>% map_int(~as.integer(.x[2]))) %>%
                  arrange(id),
                covar = var_sum_collapsed %>%
                  filter(str_detect(var_sum_collapsed$var, "covar")) %>%
                  mutate(id = str_split(var, "covar") %>% map_int(~as.integer(.x[2]))) %>%
                  arrange(id),
                cont_mat = var_sum_collapsed %>%
                  filter(str_detect(var_sum_collapsed$var, "cont_mat")) %>%
                  mutate(id = str_split(var, "cont_mat") %>% map_int(~as.integer(.x[2]))) %>%
                  arrange(id),
                cont = var_sum_collapsed %>%
                  filter(str_detect(var_sum_collapsed$var, "cont"),
                         !str_detect(var_sum_collapsed$var, "cont_mat"),
                         !str_detect(var_sum_collapsed$var, "pars_cont")) %>%
                  mutate(id = str_split(var, "cont") %>% map_int(~as.integer(.x[2]))) %>%
                  arrange(id),
                pars_cont = var_sum_collapsed %>%
                  filter(str_detect(var_sum_collapsed$var, "pars_cont")) %>%
                  mutate(id = str_split(var, "pars_cont") %>% map_int(~as.integer(.x[2]))) %>%
                  arrange(id),
                appox_err = var_sum_collapsed %>%
                  filter(str_detect(var_sum_collapsed$var, "approx_err")) 
                )


# Package and export
proj_sum <- list(n_sum = n_sum,
                 aa_sum = aa_sum,
                 aa_cont_sum = aa_cont_sum,
                 sens_sum = sens_sum,
                 var_sum = var_sum)
# write_rds(proj_full,
#           paste0("analysis/proj.rds"))
# write_rds(proj_sum,
#           paste0("analysis/proj_sum.rds"))


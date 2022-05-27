#=========================================================================================
#========== Define functions
#=========================================================================================

# function for block-diagonal matrix from submatrix
block_fn <- function(v_, r_, c_) {
  as.matrix(Reduce("bdiag", 
                   lapply(1:c_, function(x){
                     matrix(v_, nrow = r_, ncol = r_, byrow = T)
                   })))
}

#=========================================================================================





#=========================================================================================
#========== Full model
#=========================================================================================

if(model == "full") {
  
  # site-specific stage transitions (time-varying)
  qt_block <- block_fn(v_ = c(1, 0,
                              0, 1),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  qf_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (time-varying)
  gt_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  gf_block <- block_fn(v_ = c(0, 0,
                              1, 0),
                       c_ = sites, r_ = stages)
  
  # stage-specific site transitions (time-varying)
  # diagonal should be 0
  dt_block <- matrix(c(0, 0, 1, 0,
                       0, 0, 0, 1,
                       1, 0, 0, 0,
                       0, 1, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  
  # stage-specific site transitions (fixed)
  # diagonal should be 0
  df_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  # site-specific recruitment (time-varying)
  rt_block <- block_fn(v_ = c(0, 1,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific recruitment (fixed)
  rf_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
}

#=========================================================================================





#=========================================================================================
#========== No juvenile movement
#=========================================================================================

if(model == "no_juv_move") {
  
  # site-specific stage transitions (time-varyings)
  qt_block <- block_fn(v_ = c(1, 0,
                              0, 1),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  qf_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (time-varyings)
  gt_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  gf_block <- block_fn(v_ = c(0, 0,
                              1, 0),
                       c_ = sites, r_ = stages)
  
  # stage-specific site transitions (time-varying)
  # diagonal should be 0
  dt_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 1,
                       0, 0, 0, 0,
                       0, 1, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  
  # stage-specific site transitions (fixed)
  # diagonal should be 0
  df_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  # site-specific recruitment (time-varying)
  rt_block <- block_fn(v_ = c(0, 1,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific recruitment (fixed)
  rf_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
}

#=========================================================================================





#=========================================================================================
#========== No movement (juveniles or adults)
#=========================================================================================

if(model == "no_move") {
  
  # site-specific stage transitions (time-varyings)
  qt_block <- block_fn(v_ = c(1, 0,
                              0, 1),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  qf_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (time-varyings)
  gt_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  gf_block <- block_fn(v_ = c(0, 0,
                              1, 0),
                       c_ = sites, r_ = stages)
  
  # stage-specific site transitions (time-varying)
  # diagonal should be 0
  dt_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  
  # stage-specific site transitions (fixed)
  # diagonal should be 0
  df_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  # site-specific recruitment (time-varying)
  rt_block <- block_fn(v_ = c(0, 1,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific recruitment (fixed)
  rf_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
}

#=========================================================================================





#=========================================================================================
#========== No movement (juveniles or adults)
#=========================================================================================

if(model == "null") {
  
  # site-specific stage transitions (time-varyings)
  qt_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  qf_block <- block_fn(v_ = c(1, 0,
                              0, 1),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (time-varyings)
  gt_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific stage transitions (fixed)
  gf_block <- block_fn(v_ = c(0, 0,
                              1, 0),
                       c_ = sites, r_ = stages)
  
  # stage-specific site transitions (time-varying)
  # diagonal should be 0
  dt_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0,
                       0, 0, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  
  # stage-specific site transitions (fixed)
  # diagonal should be 0
  df_block <- matrix(c(0, 0, 0, 0,
                       0, 0, 0, 1,
                       0, 0, 0, 0,
                       0, 1, 0, 0),
                     nrow = stages * sites, ncol = sites * sites)
  
  # site-specific recruitment (time-varying)
  rt_block <- block_fn(v_ = c(0, 0,
                              0, 0),
                       c_ = sites, r_ = stages)
  
  # site-specific recruitment (fixed)
  rf_block <- block_fn(v_ = c(0, 1,
                              0, 0),
                       c_ = sites, r_ = stages)
}

#=========================================================================================
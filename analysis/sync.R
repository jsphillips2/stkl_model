covmat <- surv_sum %>%
  select(date, state,  mi) %>%
  pivot_wider(names_from =  state, values_from = mi) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()
covmat <- (covmat / mean(covmat))  %>% round(1)
covmat



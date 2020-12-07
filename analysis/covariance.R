test <- surv_sum %>%
  ungroup() %>%
  filter(var == "surv_s") %>%
  mutate(state = factor(interaction(stage, basin),
                        levels = c("juvenile.south",
                                   "adult.south",
                                   "juvenile.north",
                                   "adult.north"),
                        labels = c("surv_js",
                                   "surv_as",
                                   "surv_jn",
                                   "surv_an")))%>%
  select(date, state, mi) %>%
  spread(state, mi) %>%
  full_join(disp_sum %>%
              ungroup() %>%
              filter(var == "disp_s") %>%
              mutate(basin = factor(basin,
                                    levels = c("north","south"),
                                    labels = c("disp_n",
                                               "disp_s"))) %>%
              select(date, basin, mi) %>%
              spread(basin, mi)) %>%
  full_join(rec_sum %>%
              ungroup() %>%
              filter(var == "rec") %>%
              mutate(basin = factor(basin,
                                    levels = c("north","south"),
                                    labels = c("rec_n",
                                               "rec_s"))) %>%
              select(date, basin, mi) %>%
              spread(basin, mi)
  )


covm <- test %>%
  gather(var, val, -date) %>%
  group_by(var) %>%
  mutate(val = (val - mean(val))/sd(val)) %>%
  spread(var, val) %>%
  select(-date) %>%
  as.matrix() %>%
  cov()

corm <- test %>%
  select(-date) %>%
  as.matrix() %>%
  cor(method = "spearman")

corm_clean <- round(corm[rownames(covm), colnames(covm)], 2)
corm_clean[upper.tri(corm_clean)] <- 0
corm_clean


covm_clean <- round(covm, 2)
covm_clean[upper.tri(covm_clean)] <- 0
covm_clean

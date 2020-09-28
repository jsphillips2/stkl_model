# import model
name <- "full"
out_full <- readRDS(
  paste0("analysis/demographic_model/model/output/",name,"/fit_",name,".rds"))

name <- "no_move"
out_no_move <- readRDS(
  paste0("analysis/demographic_model/model/output/",name,"/fit_",name,".rds"))



# extract fit
fit <- out_full$fit

# summarize fit
fit_summary <- rstan::summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_tibble(.) %>%
    mutate(var = rownames(rstan::summary(fit)$summary))}



# extract fit
fit_nm <- out_no_move$fit

# summarize fit
fit_summary_nm <- rstan::summary(fit_nm, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_tibble(.) %>%
    mutate(var = rownames(rstan::summary(fit_nm)$summary))}

fit_summary %>%
  filter(var == "log_lik_sum")


fit_summary_nm %>%
  filter(var == "log_lik_sum")



# loo
loos <- lapply(c(fit,fit_nm), function(x){
  
  fit_ <- x
  
  log_lik <- extract_log_lik(fit_, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik))
  
  loo <- loo(log_lik, r_eff = r_eff, cores = 10)
}) %>%
  set_names(c("fit","fit_nm"))

# looic
elpd <- compare(loos[[1]],
                 loos[[2]])





x_clean_nm <- fit_summary_nm %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary_nm$var, "x"), !str_detect(fit_summary_nm$var, "x0")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) %>%
  full_join(state_match) %>%
  select(basin, stage, date, name, lo, mi, hi)

x_clean %>%
  mutate(model = "full") %>%
  bind_rows(x_clean_nm %>%
              mutate(model = "no_move")) %>%
  ggplot(aes(date, mi))+
  facet_grid(stage~basin)+
  geom_line(data = data_prep, aes(x = date, mean_scale), size = 0.1)+
  geom_point(data = data_prep,
             aes(x = date, mean_scale), size = 0.7, shape = 1, stroke = 0.3)+
  geom_line(aes(x = date, y = mi,color = model), size = 0.4)+
  scale_y_continuous("Abundance", breaks = c(300,900,1500))+
  scale_x_date("Date", 
               breaks = lubridate::as_date(c("1994-06-01","2004-06-01","2014-06-01")),
               labels = c("1994","2004","2014"))

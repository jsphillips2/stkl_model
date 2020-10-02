out_in <- read_rds("output/fit_full.rds")
fit_summary <- out_in$fit_summary

qt <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "qt"),
         !str_detect(fit_summary$var, "qte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         b = b[time]) 


qt %>%
  full_join(state_match) %>%
  mutate(mi = exp(-(2/12) * mi),
         lo = exp(-(2/12) * lo),
         hi = exp(-(2/12) * hi)) %>%
  ggplot(aes(date, mi, color = basin))+
  facet_wrap(~stage, nrow = 2)+
  geom_hline(yintercept = 0)+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi,
                  fill = basin),
              alpha = 0.1,
              linetype = 0)+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Survival probability",
                     breaks = c(0, 0.2, 0.5, 0.8),
                     limits = c(0, 1))




dt <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "dt"),
         !str_detect(fit_summary$var, "dte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time],
         basin = levels(state_match$basin)[st],
         b = b[time]) 


test <- dt %>%
  mutate(prob = 1 - exp(- b * mi),
         prob_s= 1 - exp(- 0.5 * mi)) %>%
  select(date, b, basin, mi, prob, prob_s)

test %>%
  ggplot(aes(date, prob, color = basin))+
  geom_line(size = 0.2)+  
  geom_line(aes(y = prob_s), size = 0.5)+  
  scale_color_manual(values = c("gray20","dodgerblue"))+
  scale_fill_manual(values = c("gray20","dodgerblue"))
  


gt <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "gt"),
         !str_detect(fit_summary$var, "gte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) %>%
  full_join(state_match)


gt %>%
  # mutate(mi = 1 - exp(- (1 / 12) * mi)) %>%
  ggplot(aes(date, mi, color = basin))+
  facet_wrap(~stage)+
  geom_line()




rt <- fit_summary %>% select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "rt"),
         !str_detect(fit_summary$var, "rte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         basin = levels(state_match$basin)[
           str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))],
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time])


rt %>%
  ggplot(aes(date, mi, color = basin))+
  geom_line()+
  scale_color_manual(values = c("dodgerblue","gray20"))



fit_summary %>% filter(var == "ps")

library(truncnorm)
hist(rgamma(1e6, 1.5, 1.5 / 2), breaks = 1e3, col = "gray90", lty = 0)
hist(rtruncnorm(1e6, a = 0, b = Inf, mean = 1.62, sd =0.298), breaks = 1e3, add = T, col = "blue", lty = 0)
hist(rtruncnorm(1e6, a = 0, b = Inf,mean = 0.965, sd = 0.261), breaks = 1e3, add = T, col = "darkgreen", lty = 0)
hist(rtruncnorm(1e6, a = 0, b = Inf,mean = 0.636, sd = 0.0555), breaks = 1e3, add = T, col = "red", lty = 0)


# plot 
data_prep %>%
  ggplot(aes(color = basin, fill = basin))+
  facet_grid(stage~basin)+
  geom_point(aes(x = date, mean_scale / mean(yy)), size = 2, shape = 16, alpha =0.7)+
  geom_ribbon(data = x_clean, aes(x = date, ymin = lo, ymax = hi),
              linetype = 0, alpha = 0.2)+
  geom_line(data = x_clean, aes(x = date, y = mi), size = 0.8)+
  scale_color_manual("",values = c("dodgerblue","gray35"), guide = F)+
  scale_fill_manual("",values = c("dodgerblue","gray35"), guide = F)+
  scale_x_date("Date", 
               breaks = lubridate::as_date(c("1994-06-01","2003-06-01","2012-06-01")),
               labels = c("1994","2003","2012"))
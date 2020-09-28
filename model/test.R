out_in <- read_rds("output/fit_adult_b.rds")
fit_summary <- out_in$fit_summary

qt <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "qt"),
         !str_detect(fit_summary$var, "qte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) 


qt %>%
  full_join(state_match) %>%
  mutate(mi = exp(-0.5 * mi),
         lo = exp(-0.5 * lo),
         hi = exp(-0.5 * hi)) %>%
  ggplot(aes(date, mi, color = basin))+
  facet_wrap(~stage)+
  geom_line()+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Survival probability",
                     breaks = c(0.2, 0.5, 0.8))




dt <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "dt"),
         !str_detect(fit_summary$var, "dte")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         st = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         date = date_match[time]) %>%
  full_join(state_match)


dt %>%
  # mutate(mi = 1 - exp(- (1 / 12) * mi)) %>%
  ggplot(aes(date, mi, color = basin))+
  facet_wrap(~stage)+
  geom_line()


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
  geom_line()




dd <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`, Rhat) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "DD")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))) %>%
  arrange(time, col, row)

dd %>%
  filter() %>%
  ggplot(aes(time, mi))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_grid(row~col)+
  geom_line()

dd %>%
  filter(mi > 0,
         row != col) %>%
  mutate(stage = factor(ifelse(col %in% c(1,3),"j","a"),
                        levels = c("j","a")),
         basin = factor(ifelse(col %in% c(1,2),"s","n"),
                        levels = c("s","n"))) %>%
  ggplot(aes(time, mi,color = basin))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_wrap(~stage)+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi, 
                  fill = basin),
              alpha = 0.2,
              linetype = 0)+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Dispersal probability",
                     breaks = c(0.2, 0.5, 0.8))



qq <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`, Rhat) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "QQ")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))) %>%
  arrange(time, col, row)

qq %>%
  filter() %>%
  ggplot(aes(time, mi))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_grid(row~col)+
  geom_line()

qq %>%
  filter(mi > 0) %>%
  mutate(donor_stage = factor(ifelse(col %in% c(1,3),"j","a"),
                        levels = c("j","a")),
         rec_stage = factor(ifelse(row %in% c(1,3),"j","a"),
                              levels = c("j","a")),
         basin = factor(ifelse(col %in% c(1,2),"s","n"),
                        levels = c("s","n"))) %>%
  ggplot(aes(time, mi,color = basin))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_grid(rec_stage~donor_stage)+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi, 
                  fill = basin),
              alpha = 0.2,
              linetype = 0)+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Transition probability",
                     breaks = c(0.2, 0.5, 0.8))


rr <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`, Rhat) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "RR")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))) %>%
  arrange(time, col, row)

rr %>%
  filter(mi > 0) %>%
  mutate(basin = factor(ifelse(col %in% c(1,2),"s","n"),
                        levels = c("s","n"))) %>%
  ggplot(aes(time, mi,color = basin))+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi, 
                  fill = basin),
              alpha = 0.2,
              linetype = 0)+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Fecundity")




gg <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`, Rhat) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "GG")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))) %>%
  arrange(time, col, row)


gg %>%
  filter(mi > 0) %>%
  mutate(donor_stage = factor(ifelse(col %in% c(1,3),"j","a"),
                              levels = c("j","a")),
         rec_stage = factor(ifelse(row %in% c(1,3),"j","a"),
                            levels = c("j","a")),
         basin = factor(ifelse(col %in% c(1,2),"s","n"),
                        levels = c("s","n"))) %>%
  ggplot(aes(time, mi,color = basin))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_grid(rec_stage~donor_stage)+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi, 
                  fill = basin),
              alpha = 0.2,
              linetype = 0)+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Transition probability",
                     breaks = c(0.2, 0.5, 0.8))






qq_y <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`, Rhat) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "QQ")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4])),
         date = date_match[time],
         year = lubridate::year(date)) %>%
  group_by(year, row, col) %>%
  summarize(p = prod(mi))


qq_y %>%
  filter(p > 0,
         year < 2020) %>%
  mutate(donor_stage = factor(ifelse(col %in% c(1,3),"j","a"),
                              levels = c("j","a")),
         rec_stage = factor(ifelse(row %in% c(1,3),"j","a"),
                            levels = c("j","a")),
         basin = factor(ifelse(col %in% c(1,2),"s","n"),
                        levels = c("s","n"))) %>%
  ggplot(aes(year, p,color = basin))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_grid(rec_stage~donor_stage)+
  geom_line()+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Transition probability",
                     breaks = c(0.2, 0.5, 0.8))



aa <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`, Rhat) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "AA")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))) %>%
  arrange(time, col, row)


aa %>%
  filter(mi > 0) %>%
  mutate(basin = factor(ifelse(col %in% c(1,2),"s","n"),
                        levels = c("s","n")),
         mi = ifelse((row ==1 & col == 2) | (row == 3 & col == 4), mi/ 10, mi),
         lo = ifelse((row ==1 & col == 2) | (row == 3 & col == 4), lo/ 10, lo),
         hi = ifelse((row ==1 & col == 2) | (row == 3 & col == 4), hi/ 10, hi)) %>%
  ggplot(aes(time, mi,color = basin))+
  geom_hline(yintercept = 0.5, size = 0.2)+
  facet_grid(row~col)+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi, 
                  fill = basin),
              alpha = 0.2,
              linetype = 0)+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  scale_fill_manual(values = c("dodgerblue","gray20"))+
  scale_y_continuous("Per capita contribution",
                     breaks = c(0.2, 0.5, 0.8))

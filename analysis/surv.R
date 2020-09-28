test <- fit_summary %>%
  filter(str_detect(fit_summary$var, "AA"))  %>%
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
  arrange(time, col, row) %>%
  select(year, stage, basin, row, col, `50%`)


test %>% 
  filter(col == 1)

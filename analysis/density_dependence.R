library(nlme)

td <- sens_sum %>%
  filter(var == "l") %>%
  select(year, value) %>%
  unique() %>%
  mutate(n = n_real[1:29])
  
td %>%
  ggplot(aes(n, log(value)))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_vline(xintercept = 2285.691, linetype = 2)+
  geom_smooth(method = "lm", se = F, color = "black")+
  geom_point(size = 2)

summary(gls(log(value) ~ n, data = td, correlation = corAR1(form = ~year)))



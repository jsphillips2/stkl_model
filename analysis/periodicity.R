#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(cowplot)
library(WaveletComp)



options(mc.cores = parallel::detectCores()-2)

# import data
name <- "full"
proj_sum <-read_rds(paste0("analysis/demographic_model/model/output/",name,
                           "/proj_sum_","fit_full_update",".rds"))
n_sum <- proj_sum$n_sum %>%
  mutate(year = (1991:2020)[time])
sens_sum <- proj_sum$sens_sum %>%
  mutate(year = (1991:2020)[time])

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size = 10),
                  strip.text.y = element_text(angle = -90, margin=margin(0,0,0,2)),
                  strip.text.x = element_text(margin=margin(0,0,2,0)),
                  legend.text = element_text(size = 10),
                  axis.text = element_text(size = 10, color="black"),
                  axis.title.y = element_text(angle = 90, margin=margin(0,10,0,0)),
                  axis.title.x = element_text(margin = margin(10,0,0,0))))







#==========
#========== Extract lambda
#==========

x <- {sens_sum %>%
    filter(var == "l" | var == "l_asym") %>%
    mutate(var = factor(var, 
                        levels = c("l","l_asym"),
                        labels = c("Realized", "Asymptotic"))) %>%
    select(var, year, value) %>%
    unique() %>%
    filter(var == "Asymptotic")}$value %>% log()

x_n <- length(x)

plot(x, type = "l")



#==========
#========== Autocorrelation
#==========

acf(x)





#==========
#========== Periodogram (Discrrete Fourier Transform)
#==========

fourier_data <- tibble(fourier = fft(x)[1 : (x_n / 2 + 1)]) %>%
  mutate(strength = Mod(fourier),
         frequency = 0 : (x_n / 2),
         period = x_n / frequency)

fourier_data %>%
  filter(frequency > 0) %>%
  ggplot(aes(frequency, strength))+
  geom_point(size = 3)+
  geom_segment( aes(x = frequency, xend = frequency, 
                    y = 0, yend=strength))+
  scale_y_continuous("Strength")+
  scale_x_continuous(Frequency~(year^-1),
                     breaks = seq(1, 12, by = 3),
                     labels = paste(rep(1, 4),
                                     rep("/", 4),
                                     round(x_n / seq(1, 12, by = 3), 2)))
  




#==========
#========== Wavelet analysis
#==========

df <- data.frame(x = x)
  
wave <- analyze.wavelet(df,
                        "x",
                        loess.span = 0,
                        lowerPeriod = 2,
                        upperPeriod = 29,
                        dj = 1 / 20,
                        dt = 1,
                        make.pval = T,
                        n.sim = 4000,
                        method = "white.noise")

wave_d <- as_tibble(wave$Power) %>%
  mutate(period = wave$Period) %>%
  gather(id, power, -period) %>%
  mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
         year = 1990 + id) %>%
  full_join(as_tibble(wave$Power.pval) %>%
              mutate(period = wave$Period) %>%
              gather(id, pval, -period) %>%
              mutate(id = str_split(id, "V") %>% map_int(~as.integer(.x[2])),
                     year = 1990 + id,
                     bin = ifelse(pval < 0.05, 1, 0)))

qt <- quantile(wave_d$power, prob = seq(0,1, length.out = 100))

wave_d$power_quant <- qt[cut(wave_d$power, breaks = qt)]

wave_d$power_quant[is.na(wave_d$power_quant)] <- 0

wave_d %>%
  ggplot(aes(year, period))+
  geom_tile(aes(fill = power_quant)) +
  geom_contour(aes(z = bin),
               bins = 2,
               color = "black")+
  scale_y_continuous("Period (years)",
                     trans = "log",
                     breaks = c(2,4,8,16,32))+
  scale_x_continuous("Year")+
  scale_fill_gradient2("Power (quantile)",
                       low="dodgerblue", mid="gray90", high="firebrick", 
                       midpoint = 0.1) 




wave_d %>%
  mutate(bin = ifelse(pval < 0.01, 1, 0))



which(wave$Power.pval < 0.01)

wt.image(wave, n.levels = 250,legend.params = list(lab = "wavelet power levels"))

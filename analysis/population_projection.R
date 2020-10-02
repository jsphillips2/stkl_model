# load packages
library(tidyverse)
library(Matrix)
library(matrixcalc)

options(mc.cores = parallel::detectCores()-4)


# import model
out_in <- read_rds("output/fit_full.rds")

# extract data
data_list <- out_in$data_list
fit <- out_in$fit 
fit_summary <- out_in$fit_summary 
pj <- data_list$pj
nt <- data_list$nt
n <- data_list$n
b <- data_list$b

# select row
id <- 50

# paramter names for extraction
vars <- {fit_summary %>%
    filter(str_detect(fit_summary$var, "rt") |
           str_detect(fit_summary$var, "rf") |
           str_detect(fit_summary$var, "qt") |
           str_detect(fit_summary$var, "qf") |
           str_detect(fit_summary$var, "gt") |
           str_detect(fit_summary$var, "gf") |
           str_detect(fit_summary$var, "dt") |
           str_detect(fit_summary$var, "df") | 
           str_detect(fit_summary$var, "x0"))}$var

# extract fit
extract <-  rstan::extract(fit, pars = vars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(vars) %>%
  filter(row_number() == id) %>%
  gather(var, val) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])))


x0 <- extract %>%
  filter(name == "x0") %>%
  {array(.$val, dim = n)}


rt <- extract %>%
  filter(name == "rt") %>%
  {array(.$val, dim = c(nt - 1, max(pj[, 1, 1]) - 1))} %>%
  t()

rf <- extract %>%
  filter(name == "rf") %>%
  {array(.$val, dim = c(max(pj[, 1 , 2]) - 1))}

qt <- extract %>%
  filter(name == "qt") %>%
  {array(.$val, dim = c(nt - 1, max(pj[, 1 ,3]) - 1))} %>%
  t()

qf <- extract %>%
  filter(name == "qf") %>%
  {array(.$val, dim = c(max(pj[, 1 , 4]) - 1))}

gt <- extract %>%
  filter(name == "gt") %>%
  {array(.$val, dim = c(nt - 1, max(pj[, 1 ,5]) - 1))} %>%
  t()

gf <- extract %>%
  filter(name == "gf") %>%
  {array(.$val, dim = c(max(pj[, 1 , 6]) - 1))}

dt <- extract %>%
  filter(name == "dt") %>%
  {array(.$val, dim = c(nt - 1, max(pj[, 1 ,7]) - 1))} %>%
  t()

df <- extract %>%
  filter(name == "df") %>%
  {array(.$val, dim = c(max(pj[, 1 , 8]) - 1))}
  
  
Rt <- matrix(0, nrow = n, ncol = n)
Rf <- matrix(0, nrow = n, ncol = n)
Qt <- matrix(0, nrow = n, ncol = n)
Qf <- matrix(0, nrow = n, ncol = n)
Gt <- matrix(0, nrow = n, ncol = n)
Gf <- matrix(0, nrow = n, ncol = n)
Dt <- matrix(0, nrow = n, ncol = n)
Df <- matrix(0, nrow = n, ncol = n)

RR <- array(0, dim = c(n, n, nt - 1))
QQ <- array(0, dim = c(n, n, nt - 1))
GG <- array(0, dim = c(n, n, nt - 1))
DD <- array(0, dim = c(n, n, nt - 1))
AA <- array(0, dim = c(n, n, nt - 1))
x <- array(0, dim = c(n, nt))

trans_fn <- function(m0_){
  m1_ = m0_
  diag(m1_) <- 0
  v_ <- colSums(m0_)
  m2_ <- diag(v_)
  m3_ <- m1_ - m2_
  m4_ <- expm(m3_)
  return(as.matrix(m4_))
}

x[,1] <- x0

for (t in 2:nt) {
  for (i in 1 : (n * n)) {
    Rt[pj[i, 2, 1], pj[i, 3, 1]] <- c(0, rt[, t - 1])[pj[i, 1, 1]]
    Rf[pj[i, 2, 2], pj[i, 3, 2]] <- c(0, rf)[pj[i, 1, 2]]
    Qt[pj[i, 2, 3], pj[i, 3, 3]] <- c(0, qt[, t - 1])[pj[i, 1, 3]]
    Qf[pj[i, 2, 2], pj[i, 3, 4]] <- c(0, qf)[pj[i, 1, 4]]
    Gt[pj[i, 2, 5], pj[i, 3, 5]] <- c(0, gt[, t - 1])[pj[i, 1, 5]]
    Gf[pj[i, 2, 6], pj[i, 3, 6]] <- c(0, gf)[pj[i, 1, 6]]
    Dt[pj[i, 2, 7], pj[i, 3, 7]] <- c(0, dt[, t - 1])[pj[i, 1, 7]]
    Df[pj[i, 2, 8], pj[i, 3, 8]] <- c(0, df)[pj[i, 1, 8]]
    
    R <- Rt + Rf
    Q <- trans_fn(b[t - 1] * (Qt + Qf))
    G <- trans_fn(b[t - 1] * (Gt + Gf))
    D <- trans_fn(b[t - 1] * (Dt + Df))
    
    A =  R + (D %*% G %*% Q)
    
    x[, t] <- A %*% x[, t - 1]
    
    RR[, , t - 1] <- R
    QQ[, , t - 1] <- Q
    GG[, , t - 1] <- G
    DD[, , t - 1] <- D
    AA[, , t - 1] <- A
    
  }
}

AA



surv <- lapply(2:nt, function(t_){
  return(tibble(t = t_ - 1,
                row = rep(1:n, n),
                col = rep(1:n, each = n),
                val = c(QQ[,,t_ - 1])))
}) %>%
  bind_rows() %>%
  filter(row == col) %>%
  mutate(basin = factor(ifelse(col < 3, "south", "north"),
                        levels = c("south","north")),
         stage = factor(ifelse(col %in% c(1, 3), "juvenile", "adult"),
                        levels = c("juvenile","adult")))
  

surv %>%
  ggplot(aes(t, val, color = basin))+
  facet_wrap(~stage, nrow = 2)+
  geom_line()+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  theme_bw()


disp <- lapply(2:nt, function(t_){
  return(tibble(t = t_ - 1,
                row = rep(1:n, n),
                col = rep(1:n, each = n),
                val = c(DD[,,t_ - 1])))
}) %>%
  bind_rows() %>%
  filter(row %in% c(2, 4),
         col %in% c(2, 4),
         row != col) %>%
  mutate(basin = factor(ifelse(col < 3, "south", "north"),
                        levels = c("south","north")))


disp %>%
  ggplot(aes(t, val, color = basin))+
  geom_line()+
  scale_color_manual(values = c("dodgerblue","gray20"))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1))

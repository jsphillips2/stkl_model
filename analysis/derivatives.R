library(Ryacas)

S <- matrix(c("s1s*(1-g1s)",       "r1s",               "0",                 "0",
              "s1s*g1s*(1-d1s)",  "s2s*(1-d1s)",   "s3s*g1s*d2s",    "s4s*d2s",
              "0",                   "0",               "s3s*(1-g1s)",     "r2s",
              "s1s*g1s*d1s",      "s2s*d1s",       "s3s*g1s*(1-d2s)", "s4s*(1-d2s)"), 
            nrow = 4, 
            ncol = 4,
            byrow = T)
W <- matrix(c("s1w*(1-g1w)",       "r1w",               "0",                 "0",
              "s1w*g1w*(1-d1w)",  "s2w*(1-d1w)",   "s3w*g1w*d2w",    "s4w*d2w",
              "0",                   "0",               "s3w*(1-g1w)",     "r2w",
              "s1w*g1w*d1w",      "s2w*d1w",       "s3w*g1w*(1-d2w)", "s4w*(1-d2w)"), 
            nrow = 4, 
            ncol = 4,
            byrow = T)


SS <- ysym(S)
WW <- ysym(W)
NN <- ysym(N)

SSWW <- WW %*% SS

pars <- c("s1s",
          "s2s",
          "s3s",
          "s4s",
          "g1s",
          "d1s",
          "d2s",
          "s1w",
          "s2w",
          "s3w",
          "s4w",
          "g1w",
          "d1w",
          "d2w")



derivs <- deriv(SSWW,
                pars)


N <- c("x1","x2","x3","x4")


L <- sum(WW %*% SS * N) / sum(NN)

L_derivs <- deriv(L,
                c(pars, N))

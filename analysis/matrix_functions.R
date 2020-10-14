### define variables to differentiate with respect to 
theta_names <- c("r1s","r2s",
               "s1s","s2s","s3s","s4s",
               "g1s","d1s","d2s",
               "r1w","r2w",
               "s1w","s2w","s3w","s4w",
               "g1w","d1w","d2w")





### define elements of annual projection matrix
a_list <- list(a11 = expression(s1w*(1-g1w)*s1s*(1-g1s)+r1w*s1s*g1s*(1-d1s)),
               a21 = expression(s1w*g1w*(1-d1w)*s1s*(1-g1s)+
                                  s2w*(1-d1w)*s1s*g1s*(1-d1s)+s4w*d2w*s1s*g1s*d1s),
               a31 = expression(r2w*s1s*g1s*d1s),
               a41 = expression(s1w*g1w*d1w*s1s*(1-g1s)+s2w*d1w*s1s*g1s*(1-d1s)+
                                  s4w*(1-d2w)*s1s*g1s*d1s),
               a12 = expression(s1w*(1-g1w)*r1s+r1w*s2s*(1-d1s)),
               a22 = expression(s1w*g1w*(1-d1w)*r1s+s2w*(1-d1w)*s2s*(1-d1s)+
                                  s4w*d2w*s2s*d1s),
               a32 = expression(r2w*s2s*d1s),
               a42 = expression(s1w*g1w*d1w*r1s+s2w*d1w*s2s*(1-d1s)+
                                  s4w*(1-d2w)*s2s*d1s),
               a13 = expression(r1w*s3s*g1s*d2s),
               a23 = expression(s2w*(1-d1w)*s3s*g1s*d2s+s3w*g1w*d2w*s3s*(1-g1s)+
                                  s4w*d2w*s3s*g1s*(1-d2s)),
               a33 = expression(s3w*(1-g1w)*s3s*(1-g1s)+r2w*s3s*g1s*(1-d2s)),
               a43 = expression(s2w*d1w*s3s*g1s*d2s+s3w*g1w*(1-d2w)*s3s*(1-g1s)+
                                  s4w*(1-d2w)*s3s*g1s*(1-d2s)),
               a14 = expression(r1w*s4s*d2s),
               a24 = expression(s2w*(1-d1w)*s4s*d2s+s3w*g1w*d2w*r2s+s4w*d2w*s4s*(1-d2s)),
               a34 = expression(s3w*(1-g1w)*r2s+r2w*s4s*(1-d2s)),
               a44 = expression(s2w*d1w*s4s*d2s+s3w*g1w*(1-d2w)*r2s+
                                  s4w*(1-d2w)*s4s*(1-d2s)))





### define function for evaluating matrix elements 
a_mat_fn <- function(v_){
  
  r1s <- v_["r1s"]; r2s <- v_["r2s"]; 
  s1s <- v_["s1s"]; s2s <- v_["s2s"]; s3s <- v_["s3s"]; s4s <- v_["s4s"]
  g1s <- v_["g1s"]; d1s <- v_["d1s"]; d2s <- v_["d2s"]
  r1w <- v_["r1w"]; r2w <- v_["r2w"]
  s1w <- v_["s1w"]; s2w <- v_["s2w"]; s3w <- v_["s3w"]; s4w <- v_["s4w"]
  g1w <- v_["g1w"]; d1w <- v_["d1w"]; d2w <- v_["d2w"]
  
  a11 = eval(a_list[[1]]); a21 = eval(a_list[[2]])
  a31 = eval(a_list[[3]]); a41 = eval(a_list[[4]])
  
  a12 = eval(a_list[[5]]); a22 = eval(a_list[[6]])
  a32 = eval(a_list[[7]]); a42 = eval(a_list[[8]])
  
  a13 = eval(a_list[[9]]); a23 = eval(a_list[[10]])
  a33 = eval(a_list[[11]]); a43 = eval(a_list[[12]])
  
  a14 = eval(a_list[[13]]); a24 = eval(a_list[[14]])
  a34 = eval(a_list[[15]]); a44 = eval(a_list[[16]])
  
  matrix(c(a11, a12, a13, a14,
           a21, a22, a23, a24,
           a31, a32, a33, a34,
           a41, a42, a43, a44),
         nrow = 4,
         ncol = 4,
         byrow = T)
  
}





### define function for derivatives
a_deriv_fn <- lapply(a_list,
                     function(x_){
                       xx_ <- x_
                       deriv(xx_, 
                             theta_names,
                             function(r1s,r2s,s1s,s2s,s3s,s4s,g1s,d1s,d2s,
                                      r1w,r2w,s1w,s2w,s3w,s4w,g1w,d1w,d2w){})})





### define function for evaluating derivatives
a_deriv_eval_fn <- function(v_, i_){
  
  r1s <- v_["r1s"]; r2s <- v_["r2s"]; 
  s1s <- v_["s1s"]; s2s <- v_["s2s"]; s3s <- v_["s3s"]; s4s <- v_["s4s"]
  g1s <- v_["g1s"]; d1s <- v_["d1s"]; d2s <- v_["d2s"]
  r1w <- v_["r1w"]; r2w <- v_["r2w"]
  s1w <- v_["s1w"]; s2w <- v_["s2w"]; s3w <- v_["s3w"]; s4w <- v_["s4w"]
  g1w <- v_["g1w"]; d1w <- v_["d1w"]; d2w <- v_["d2w"]
  
  a_deriv_fn[[i_]](r1s,r2s,s1s,s2s,s3s,s4s,g1s,d1s,d2s,
               r1w,r2w,s1w,s2w,s3w,s4w,g1w,d1w,d2w)
}


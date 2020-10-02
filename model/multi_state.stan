//=======================================================================================

// MULTISTATE MATRIX POPULATION MODEL WITH TIME-VARYING PARAMETERS

//=======================================================================================




//=======================================================================================


functions {
  
  // convert transition rate matrix to transition probability matrix
  matrix trans_prob(matrix m0) {
    
    // declare variables
    int k = cols(m0);
    vector [k] v;
    matrix [k, k] m1;
    matrix [k, k] m2;
    matrix [k, k] m3;
    matrix [k, k] m4;
    
    // copy original transition rate matrix
    m1 = m0; 
    
    // replace diagonal elements of transition rate matrix with 0s
    for (i in 1:k) {
      m1[i, i] = 0; 
    }
    
    // sum rates for each starting state
    for (i in 1:k) {
      v[i] = sum(m0[,i]); 
    } 
    
    // diagonal matrix with summed rates
    m2 = diag_matrix(v); 
    
    // fill diagonal of transition rate matrix with negative summed rates
    m3 = m1 - m2; 
    
    // exponentiate transition matrix (solution to ODE for markov process)
    m4 = matrix_exp(m3); 
    
    // Return
    return m4;
  }
  
}


//=======================================================================================


data {
  
  // declare variables
  int <lower = 1> n; // number of states
  int <lower = 1> nt; // number of time-steps
  matrix <lower = 0> [n, nt] y; // observed abundance
  int <lower = 0> j[8]; // number of positions for each demographic parameter
  int <lower = 0> pj [n * n, 3, 8]; // position array for rates
  vector [nt - 1] b; // season index and tim step sizes 
  real p[4]; // values for priors

}


//=======================================================================================


parameters {
  
  // declare variables
  vector <lower = 0> [n] x0; // initial abundance
  matrix <lower = 0> [j[1], nt - 1] rt; // recruitment
  vector <lower = 0> [j[2]] rf; // recruitment (fixed)
  matrix <lower = 0> [j[3], nt - 1] qt; // transition rate
  vector <lower = 0> [j[4]] qf; // transition rate (fixed)
  matrix <lower = 0> [j[5], nt - 1] gt; // transition rate
  vector <lower = 0> [j[6]] gf; // transition rate (fixed) 
  matrix <lower = 0> [j[7], nt - 1] dt; // movement rate
  vector <lower = 0> [j[8]] df; // movement rate (fixed)
  real <lower = 0> ys; // standard deviation for likelihood
  real <lower = 0> rs; // scale for recruitment random walk
  real <lower = 0> ps; // scale for transition / movement rate random walk

}


//=======================================================================================


transformed parameters {
  
  // declare variables
  matrix [n, nt] x; // abundance
  
  
  
  // initial abundance
  x[, 1] = x0; 
  
   
  
  // generate demographic matrices and project dynamics
  {
    matrix [n, n] Rt; // recruitment matrix
    matrix [n, n] Rf; // recruitment matrix (fixed)
    matrix [n, n] Qt; // transition rate matrix
    matrix [n, n] Qf; // transition rate matrix (fixed)
    matrix [n, n] Gt; // transition rate matrix
    matrix [n, n] Gf; // transition rate matrix (fixed)
    matrix [n, n] Dt; // movement rate matrix
    matrix [n, n] Df; // movement rate matrix (fixed)
    matrix [n, n] Q; // transition probability matrix
    matrix [n, n] G; // transition probability matrix
    matrix [n, n] D; // transition probability matrix
    matrix [n, n] R; // transition probability matrix
    matrix [n, n] A; // full projection matrix
    
    
    
    // initialize matrices with 0's
    Rt = rep_matrix(0, n, n); 
    Rf = rep_matrix(0, n, n);
    Qt = rep_matrix(0, n, n); 
    Qf = rep_matrix(0, n, n);
    Gt = rep_matrix(0, n, n); 
    Gf = rep_matrix(0, n, n);
    Dt = rep_matrix(0, n, n); 
    Df = rep_matrix(0, n, n);
    
    
    
    // loop over time
    for (t in 2:nt) { 
      
      // fill demographic matrices
        if(j[1] > 0 || j[2] > 0 || j[3] > 0 || j[4] > 0) {
          
          for (i in 1 : (n * n)) {
            // time-varying recruitment
            if (j[1] > 0) {
                Rt[pj[i, 2, 1], pj[i, 3, 1]] = 
                  append_row(0, rt[, t - 1])[pj[i, 1, 1]];
            }
            
            // fixed recruitment
            if (j[2] > 0) { 
                Rf[pj[i, 2, 2], pj[i, 3, 2]] = 
                 append_row(0, rf)[pj[i, 1, 2]];
            }
            
            // time-varying transition rate
            if (j[3] > 0) { 
                Qt[pj[i, 2, 3], pj[i, 3, 3]] = 
                  append_row(0, qt[, t - 1])[pj[i, 1, 3]];
            }
            
            // fixed transition
            if (j[4] > 0) { 
                Qf[pj[i, 2, 4], pj[i, 3, 4]] = 
                  append_row(0, qf)[pj[i, 1, 4]];
            }
            
            // time-varying transition rate
            if (j[5] > 0) { 
                Gt[pj[i, 2, 5], pj[i, 3, 5]] = 
                  append_row(0, gt[, t - 1])[pj[i, 1, 5]];
            }
            
            // fixed transition
            if (j[6] > 0) { 
                Gf[pj[i, 2, 6], pj[i, 3, 6]] = 
                  append_row(0, gf)[pj[i, 1, 6]];
            }
            
            // time-varying movement rate
            if (j[7] > 0) { 
                Dt[pj[i, 2, 7], pj[i, 3, 7]] = 
                  append_row(0, dt[, t - 1])[pj[i, 1, 7]];
            }
            
            // fixed movement
            if (j[8] > 0) { 
                Df[pj[i, 2, 8], pj[i, 3, 8]] = 
                  append_row(0, df)[pj[i, 1, 8]];
            }    
            
          } // n
        }
      
      
      
      // calculate transition probability matrix from transition rate matrix
      Q = trans_prob(b[t - 1] * (Qt + Qf)); 
      G = trans_prob(b[t - 1] * (Gt + Gf)); 
      D = trans_prob(b[t - 1] * (Dt + Df));
      R = Rt + Rf;
      
      // projection matrix
      A =  R + (D * G * Q);

      
      // project dynamics
      x[, t] = A * x[, t - 1]; 
      
    } // t
    
  }
}


//=======================================================================================


model {
  
  // scales
  ys ~ gamma(1.5, 1.5 / p[1]); 
  rs ~ gamma(1.5, 1.5 / p[2]); 
  ps ~ gamma(1.5, 1.5 / p[3]);
  
  // fixed rates
  if(j[2] > 0) {
    for (i in 1:j[2]) {
      rf[i] ~ gamma(1.5, 1.5 / p[2]); 
    }
  }
  if(j[4] > 0) {
    for (i in 1:j[4]) {
      qf[i] ~ gamma(1.5, 1.5 / p[3]); 
    }
  }
  if(j[6] > 0) {
    for (i in 1:j[6]) {
      gf[i] ~ gamma(1.5, 1.5 / p[3]); 
    }
  }  
  if(j[8] > 0) {
    for (i in 1:j[8]) {
      df[i] ~ gamma(1.5, 1.5 / p[3]); 
    }
  }  
  

  // time-varying rates
    if (j[1] > 0 || j[3] > 0 || j[5] > 0 || j[7] > 0) {
     for (t in 1:nt) {
        if (j[1] > 0) {
          for (i in 1:j[1]) {
            if (t == 1) {
              rt[i, t] ~ gamma(1.5, 1.5 / p[2]); 
            }
            if (t > 1 && t < nt) {
              rt[i, t] ~ normal(rt[i, t - 1], rs) T[0, ]; 
            }
          } // i
        }
        if (j[3] > 0) {
          for (i in 1:j[3]) {
            if (t == 1) {
              qt[i, t] ~ gamma(1.5, 1.5 / p[3]); 
            }
            if (t > 1 && t < nt) {
              qt[i, t] ~ normal(qt[i, t - 1], ps) T[0, ]; 
            }
          } // i
        }
        if (j[5] > 0) {
          for (i in 1:j[5]) {
            if (t == 1) {
              gt[i, t] ~ gamma(1.5, 1.5 / p[3]); 
            }
            if (t > 1 && t < nt) {
              gt[i, t] ~ normal(gt[i, t - 1], ps) T[0, ]; 
            }
          } // i
        }
        if (j[7] > 0) {
          for (i in 1:j[7]) {
            if (t == 1) {
              dt[i, t] ~ gamma(1.5, 1.5 / p[3]); 
            }
            if (t > 1 && t < nt) {
              dt[i, t] ~ normal(dt[i, t - 1], ps) T[0, ]; 
            }
          } // i
        }  
      } // t 
    }
    
    
  
  // initial abundance
  for (i in 1:n) {
    x0[i] ~ gamma(1.5, 1.5 / p[4]);
  } // i
  
  
  
  // likelihood
  for (i in 1:n){
    y[i, 1:nt] ~ normal(x[i, 1:nt], ys); 
  } // i

  
}


//=======================================================================================


generated quantities {
  
  // declare variables
  real log_lik [n * nt]; // pointwise log-likelihood
  real log_lik_sum; // total log-likelihood
  
  {
    
    // convert matrices to vectors
    vector [n * nt] yy;
    vector [n * nt] xx;
    yy = to_vector(y);
    xx = to_vector(x);
    
    // pointwise log-likelihood
    for (i in 1 : (n * nt)) {
      log_lik[i] = normal_lpdf(yy[i] | xx[i], ys);
    }
    
  }  
  
  // total log-likelihood
  log_lik_sum = sum(log_lik);
  
}


//=======================================================================================


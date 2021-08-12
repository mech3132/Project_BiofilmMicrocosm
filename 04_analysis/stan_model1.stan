//Stan model for simple linear regression
data {
  int <lower = 0> N; // Sample size
  int <lower = 1> D; // number predictors
  int <lower = 1> L; // number of Levels
  int <lower = 1, upper = L> ll[N]
  row_vector[D] x[N];
  // matrix[N,K] X; // data matrix
  vector[N] y; // Outcome
}
parameters {
  real a; //Intercept
  vector[K] b; // vector of coefficients
  //real bCV; // CV_log10
  //real bR; // B_obs_rich
  real < lower = 0 > sigma; //Error SD
}
model{
  y ~ normal(a + X*b, sigma);
}
generated quantities {
  real y_rep[N];
  
  for (n in 1:N) 
    y_rep[n] = normal_rng(a + X[n]*b, sigma);
    
  //real y_rep[N] = normal_rng(a + xCV * bCV + xR * bR, sigma);
} // posterior predictive distribution

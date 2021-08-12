// Testing out model
  data {
  int<lower = 1> N;               // Number of observations.
  //int<lower = 1> K;               // Number of predictors
  vector[N] y;                    // Vector of observations.
  //matrix[N,K] x;                  // predictor matrix
  
  vector[N] xC;                   //CV
  vector[N] xR;                   //Richness
  }
  parameters {
  real a;                         // Intercept
  real<lower = 0> sigma;          // Variance of the likelihood.
  real <lower=1> iC;                        // Effect of richness on C (invisible)
  real <lower=0> sigC;            //sigma of CV
  real bC;                        // predictor for C
  real bR;                        // predictor for R
  }
  model{
  // Priors.
  a ~ normal(0,5);
  bC ~ normal(0,5);
  bR ~ normal(0,5);
  sigma ~ exponential(0.88);
  
  //Latent priors
  iC ~ exponential(1);
  sigC ~ exponential(1);
  
  // latent variables
  xC ~ normal(xR*iC, sigC);
  
  // Model
  y ~ normal(a + xC*bC + xR*bR, sigma);
  }
  generated quantities {
  //vector[N] y_pred;
  //for (n in 1:N) {
  //  y_pred[n] = normal_rng(a + xC*bC + xR*bR, sigma);
  //}
  }
  

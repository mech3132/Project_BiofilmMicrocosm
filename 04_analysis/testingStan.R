### Random code tidbits


### Trying custom stan 
dat_for_stan <- dat_forBayes %>% select(CV_log10, B_obsRich,Micro_Bd_density_log10, DateStart) %>% drop_na()
ll <- match(dat_for_stan$DateStart, levels(dat_for_stan$DateStart))
# stan_data <- list(N = nrow(dat_for_stan)
#                   , xCV = dat_for_stan$CV_log10
#                   , xR = dat_for_stan$B_obsRich
#                   , lD = dat_for_stan$DateStart
#                   , y = dat_for_stan$Micro_Bd_density_log10
#                   )
stan_data <- list(N = nrow(dat_for_stan)
                  , L = length(levels(dat_for_stan$DateStart))
                  , ll = ll
                  , X = as.matrix(cbind(xC=dat_for_stan$CV_log10, xR=dat_for_stan$B_obsRich))
                  , K = ncol(as.matrix(cbind(xC=dat_for_stan$CV_log10, xR=dat_for_stan$B_obsRich)))
                  , y = dat_for_stan$Micro_Bd_density_log10
)

lm_model <- lm(stan_data$y ~ stan_data$xCV + stan_data$xR)
lm_coef <- coef(lm_model)
lm_coef <- t(lm_coef) %>% as_tibble() %>% rename(Intercept = "(Intercept)", CV = "stan_data$xCV", Rich = "stan_data$xR")
ggplot(dat_for_stan, aes(x=CV_log10, y=Micro_Bd_density_log10)) + geom_point() +
  geom_abline(data = lm_coef, aes(intercept=Intercept, slope=CV))

plot(lm_model$fitted.values ~ dat_for_stan$Micro_Bd_density_log10)

write("//Stan model for simple linear regression
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
} // posterior predictive distribution", 
"stan_model1.stan")

write(rstan::get_stancode(stanlmer_micro$stanfit), "test_stan.stan")
# Check
stanc("stan_model1.stan")
# Save filepath
stan_model1 <- "stan_model1.stan"

fit <- stan(file = stan_model1, data = stan_data, warmup = 500, iter = 1000, chains = 4, cores = 1, thin = 1)

posterior <- rstan::extract(fit)
posterior$b
stan_coef <- data.frame(Intercept = mean(posterior$a), bCV = mean(posterior$bCV), bR = mean(posterior$bR))
stan_pred <- data.frame(y_obs = dat_for_stan$Micro_Bd_density_log10
                        , y_pred = colMeans(posterior$y_rep)
                        , y_pred_lwr95 = apply(posterior$y_rep, MARGIN=2, function(x) quantile(x, 0.025))
                        , y_pred_upr95 = apply(posterior$y_rep, MARGIN=2, function(x) quantile(x, 0.975)))

ggplot(dat_for_stan, aes(x=CV_log10, y=Micro_Bd_density_log10)) + geom_point() +
  geom_abline(data = lm_coef, aes(intercept=Intercept, slope=CV)) +
  geom_abline(data = stan_coef, aes(intercept=Intercept, slope=bCV), col = "red", lty=2)

ggplot(dat_for_stan, aes(x=B_obsRich, y=Micro_Bd_density_log10)) + geom_point() +
  geom_abline(data = lm_coef, aes(intercept=Intercept, slope=Rich)) +
  geom_abline(data = stan_coef, aes(intercept=Intercept, slope=bR), col = "red", lty=2)

ggplot(stan_pred,aes(x=y_obs, y=y_pred)) +geom_point() #+ geom_segment(aes(x = y_obs, xend = y_obs, y = y_pred_lwr95, yend = y_pred_upr95 ), col="red") 




### Try just predictors

stanlmer_CV <- stan_lmer(CV_log10 ~ B_obsRich*Inhibitory + (1|DateStart), data = dat_forBayes
                         , prior = normal(0,5,autoscale=TRUE)
                         # , prior_intercept = exponential(autoscale=TRUE)
                         , adapt_delta = 0.999)
summary(stanlmer_CV)
# 
# stanlmer_micro_samps <- rstan::extract(stanlmer_micro$stanfit)
# all_samps_stanlmer <- stanlmer_micro_samps$alpha %>% as_tibble() %>% dplyr::rename(Intercept="V1") %>%
#   cbind(stanlmer_micro_samps$beta) %>% as_tibble() %>%
#   dplyr::rename(CV_log10=`1`, B_obsRich=`2`, qPCR_log10_Bd=`3`, Inhibitory = `4`, B_obsRich_Inhibitory = `5`) %>%
#   cbind(stanlmer_micro_samps$b[,1:6]) %>% as_tibble() %>%
#   dplyr::rename_at(vars(c(`1`,`2`, `3`,`4`,`5`,`6`)), ~levels(dat_forBayes$DateStart)) 
# # Get quantiles
# all_samps_stanlmer %>% pivot_longer(everything(), names_to = "Predictor", values_to = "Sample") %>%
#   group_by(Predictor) %>% mutate(lwr95 = quantile(Sample, 0.025), upr95 = quantile(Sample, 0.975)) %>% ungroup() %>% 
#   ggplot(aes(x=Sample)) + geom_histogram() + facet_wrap(Predictor~., scales="free") +
#   geom_vline(aes(xintercept=upr95), col="red", lty=2) + geom_vline(aes(xintercept=lwr95), col="red", lty=2)



stanlmer_qPCR <- stan_lmer(qPCR_log10_Bd ~ B_obsRich*Inhibitory + (1|DateStart), data = dat_forBayes
                           , prior = normal(0,5,autoscale=TRUE)
                           # , prior_intercept = exponential(autoscale=TRUE)
                           , adapt_delta = 0.999)
summary(stanlmer_qPCR)
# 
# stanlmer_micro_samps <- rstan::extract(stanlmer_micro$stanfit)
# all_samps_stanlmer <- stanlmer_micro_samps$alpha %>% as_tibble() %>% dplyr::rename(Intercept="V1") %>%
#   cbind(stanlmer_micro_samps$beta) %>% as_tibble() %>%
#   dplyr::rename(CV_log10=`1`, B_obsRich=`2`, qPCR_log10_Bd=`3`, Inhibitory = `4`, B_obsRich_Inhibitory = `5`) %>%
#   cbind(stanlmer_micro_samps$b[,1:6]) %>% as_tibble() %>%
#   dplyr::rename_at(vars(c(`1`,`2`, `3`,`4`,`5`,`6`)), ~levels(dat_forBayes$DateStart)) 
# # Get quantiles
# all_samps_stanlmer %>% pivot_longer(everything(), names_to = "Predictor", values_to = "Sample") %>%
#   group_by(Predictor) %>% mutate(lwr95 = quantile(Sample, 0.025), upr95 = quantile(Sample, 0.975)) %>% ungroup() %>% 
#   ggplot(aes(x=Sample)) + geom_histogram() + facet_wrap(Predictor~., scales="free") +
#   geom_vline(aes(xintercept=upr95), col="red", lty=2) + geom_vline(aes(xintercept=lwr95), col="red", lty=2)





##### DATA #####
dat_for_stan <- dat %>% filter(Rep !="Rep13") %>% 
  select(Micro_Bd_density_log10,DateStart, CV_log10, B_obsRich, qPCR_log10_Bd, qPCR_log10_bact) %>% drop_na()
# gg <- match(dat_for_stan$DateStart, levels(dat_for_stan$DateStart))
stan_data <- list(N = nrow(dat_for_stan)
                  # , x = cbind(CV=dat_for_stan$CV_log10, R = dat_for_stan$B_obsRich)
                  , xC = dat_for_stan$CV_log10
                  , xR = dat_for_stan$B_obsRich
                  # , K = 2
                  # , g = gg
                  # , x = as.matrix(cbind(xC=dat_for_stan$CV_log10))
                  # , K = ncol(as.matrix(cbind(xC=dat_for_stan$CV_log10)))
                  , y = dat_for_stan$Micro_Bd_density_log10
)

##### STAN MODEL #####

write("// Testing out model
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
  iC ~ normal(0,5);
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
  ",
      "basic_stan_model.stan")
stanc("basic_stan_model.stan")
basic_stan_model <- "basic_stan_model.stan"

##### RUN AUTO RSTANARM #####
fit_auto <- stan_lmer(Micro_Bd_density_log10 ~ CV_log10 + B_obsRich + (1|DateStart)
                    , data=dat_for_stan
                    , prior = normal(0,5, autoscale = TRUE))
prior_summary(fit_auto)
##### RUN CUSTOM RSTANARM #####
fit_custom <- stan(file = basic_stan_model, data = stan_data
                   , warmup = 500, iter = 1000, chains = 4, cores = 1, thin = 1
                   ,  control = list(adapt_delta = 0.9999))

#### COMPARE OUTPUTS #####
posterior <- rstan::extract(fit_custom)
fit_auto

hist(posterior$a)
hist(posterior$iC)
hist(posterior$sigC)
hist(posterior$bC)
hist(posterior$bR)
hist(posterior$sigma)
posterior$alpha

summary(fit_auto)

stan()

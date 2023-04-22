functions {
  real[] sir2(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real beta = theta[1];
      real gamma = theta[2];
      real initS = theta[3];
      real initI = theta[4];

      real S = y[1] + initS;
      real I = y[2] + initI;
      real R = y[3];
      real N = S + I + R;
      real dy_dt[3];

      dy_dt[1] = -beta * I * S / N;
      dy_dt[2] =  beta * I * S / N - gamma * I;
      dy_dt[3] =  gamma * I;
      return(dy_dt);
  }
  real[,] simulate_sir2(real[] t, real initS, real initI, real beta, real gamma) {
    int N_t = num_elements(t);
    real t0 = -0.1;
    real theta[4] = {beta, gamma, initS, initI};
    real y[N_t, 3];
    y = integrate_ode_rk45(sir2, rep_array(0.0, 3), t0, t, theta, rep_array(0.0, 0), rep_array(0, 0));
    return(y);
  }
  
  int[] simulate_cases2_rng(real[] t, real initS, real initI, real beta, real gamma, real phi) {
    int N_t = num_elements(t);
    real y[N_t, 3] = simulate_sir2(t, initS, initI, beta, gamma);
    int cases[N_t] = neg_binomial_2_rng(y[,2], phi);
    return(cases);
  }
}
data {
  int N_t;
  real t[N_t];
  int cases[N_t];
}
parameters {
  real<lower = 0> beta;
  real<lower = 0> gamma;
  real<lower = 0> phi_inv;
  real log_initS_raw;
  real<lower = 0> initI;
}
transformed parameters {
  real phi = 1. / phi_inv;
  real<lower = 0> initS = 1000 + exp(log_initS_raw);
  real ys[N_t, 3] = simulate_sir2(t, initS, initI, beta, gamma);
}
model {
  // priors
  phi_inv ~ exponential(5);
  beta ~ normal(0, 5);
  gamma ~ normal(0.4, 0.5);
  log_initS_raw ~ normal(5, 10);
  initI ~ std_normal();
  
  // likelihood
  cases ~ neg_binomial_2(ys[,2], phi);
}
generated quantities {
  int yrep[N_t] = neg_binomial_2_rng(ys[,2], phi);
  real log_lik[N_t];
  for (i in 1:N_t)
    log_lik[i] = neg_binomial_2_lpmf(cases[i] | ys[i,2], phi);
}

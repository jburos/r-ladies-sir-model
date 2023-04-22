functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {

      real beta = theta[1];
      real gamma = theta[2];

      real S = y[1];
      real I = y[2];
      real R = y[3];
      real N = S + I + R;
      real dy_dt[3];

      dy_dt[1] = -beta * I * S / N;
      dy_dt[2] =  beta * I * S / N - gamma * I;
      dy_dt[3] =  gamma * I;
      return(dy_dt);
  }
  real[,] simulate_sir(real[] t, data int N, data int i0, real beta, real gamma) {
    int N_t = num_elements(t);
    real y0[3] = {N - i0 + 0.0, i0 + 0.0, 0.0};
    real t0 = -0.1;
    real theta[2] = {beta, gamma};
    real y[N_t, 3];
    y = integrate_ode_rk45(sir, y0, t0, t, theta, rep_array(0.0, 0), rep_array(0, 0));
    return(y);
  }
  
  int[] simulate_cases_rng(real[] t, data int N, data int i0, real beta, real gamma, real phi) {
    int N_t = num_elements(t);
    real y[N_t, 3] = simulate_sir(t, N, i0, beta, gamma);
    int cases[N_t] = neg_binomial_2_rng(y[,2], phi);
    return(cases);
  }
}
model {
}

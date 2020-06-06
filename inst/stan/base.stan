functions {
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan

  vector test_csr_matrix_times_vector(int m, int n, vector w,
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b);
  }
}

data {
  int <lower=1> M; // number of countries
  int <lower=1> N0; // number of days for which to impute infections
  int <lower=1> starts[M]; // the start index of each group
  int<lower=1> NC[M]; // days of observed data for each group.
  int<lower=1> N2; // total period for the simulation
  int<lower=1> NS; // maximum number of simulation days for any given group
  int<lower=0> N_obs; // total size of the observation vector
  int<lower=0> R; // number of different observation types
  real<lower=0> noise_scale[R];
  int obs[N_obs]; // vector of observations
  int obs_group[N_obs]; // group (1 to M) to which each observation belongs
  int obs_type[N_obs]; // type of observation (1 to r). 
  int obs_date[N_obs]; // observation date (1 to N2)
  matrix[NS, R] P; // the 'pvec' for each type of observation.
  matrix[M, R] means; // mean values for observed events / total cases (for example, IFR)
  real pop[M];
  real SI[NS]; // fixed SI using empirical data
#include /data/NKX.stan
#include /data/data_glm.stan
#include /data/hyperparameters.stan
#include /data/glmer_stuff.stan
#include /data/glmer_stuff2.stan
}

transformed data {
  real aux = not_a_number();
  vector[NS] SI_rev; // SI in reverse order
  vector[NS] P_rev[M]; // P in reversed order
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
#include /tdata/tdata_glm.stan

  for(i in 1:NS)
    SI_rev[i] = SI[NS-i+1];

  for(m in 1:M){
    for(i in 1:NS) {
      P_rev[m, i] = P[N2-i+1,m];
    }
  }
}

parameters {
  real gamma[has_intercept];
#include /parameters/parameters_glm.stan
  real<lower=0> mu[M];
  real<lower=0> kappa;
  real<lower=0> y[M];
  real<lower=0> phi;
  real<lower=0> tau2;
  real<lower=0> noise[M, R];
}

transformed parameters {
  real E_obs[N_obs]; // expected values of the observations 
  vector[N] R0_vec;
  vector[N] Rt_vec;
  matrix[N2, M] Rt = rep_matrix(0,N2,M);
  matrix[N2, M] prediction = rep_matrix(0,N2,M);
  vector[N] eta;  // linear predictor
  
#include /tparameters/tparameters_glm.stan
#include /model/make_eta.stan


  { // generate vector of R0s from group means
    int idx = NC[1]+1;
    R0_vec[1:NC[1]] = rep_vector(mu[1], NC[1]);
    for (m in 2:M) {
      R0_vec[idx:(idx+NC[m]-1)] = rep_vector(mu[m], NC[m]);
      idx += NC[m];
    }
  }

  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p,
                             1.0, tau, scale, zeta, rho, z_T);
      b = make_b(z_b, theta_L, p, l);
    }
  }

  if (t > 0) {
#include /model/eta_add_Zb.stan
  }
  if (has_intercept == 1) {
    eta += gamma[1];
  }
  else {
#include /model/eta_no_intercept.stan
  }

  // todo: Add branching logic for different link functions.
  Rt_vec = R0_vec * 2 .* inv_logit(-eta);

  { // predict cases over time
    int idx=1, n0, n1, n2;
    matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
    for (m in 1:M){
      // time indices for group m: start date, final seed date, final date
      n0 = starts[m];
      n1 = n0 + N0 - 1;
      n2 = n0 + NC[m] - 1;

      // impute unadjusted Rt from the linear predictor
      Rt[n0:n2,m] = Rt_vec[idx:(idx+NC[m]-1)];
      idx += NC[m];

      prediction[n0:n1,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
      cumm_sum[n0:n1,m] = cumulative_sum(prediction[n0:n1,m]);

      for(i in (n1+1):n2) {
        real convolution = dot_product(sub_col(prediction, n0, m, i-n0), tail(SI_rev, i-n0));
        prediction[i,m] = (pop[m] - cumm_sum[i-1,m]) * (1 - exp(-Rt[i,m] * convolution / pop[m])) 
        cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i,m];
      }
    }
  }

  {
    // compute expected values of the observations
    int m, dt, tp, n0;
    for (i in 1:N_obs) {
      m = obs_group[i];
      dt = obs_date[i];
      tp = obs_type[i];
      n0 = starts[m];
      if (dt == 1)
        E_obs[i] = 1e-15 * prediction[1,m];
      else
        E_obs[i] = noise[m, tp] * means[m, tp] * dot_product(sub_col(prediction, n0, m, dt-n0), tail(P_rev[tp], dt-n0));
    }
  }
  
}

model {
  tau2 ~ exponential(0.03);
  for (m in 1:M) {
    y[m] ~ exponential(1/tau2);
  }
  phi ~ normal(0,5);
  kappa ~ normal(0,0.5);
  mu ~ normal(3.28, kappa);

  for (r in 1:R)
    noise[,r] ~ normal(1, noise_scale[r]);

#include /model/priors_glm.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau,
                          regularization, delta, shape, t, p);
  }

  if (prior_PD == 0) {
    for (i in 1:N_obs)
      obs[i] ~ neg_binomial_2(E_obs[i], phi);
  }
}

generated quantities {
  real alpha[has_intercept];

  matrix[N2,M] prediction_out = prediction;
  matrix[N2,M] Rt_adj_out = Rt_adj;
  
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  }
  
}

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
  int<lower=1> NC[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int deaths[N2, M];
  matrix[N2, M] f; // ifr
  real pop[M];
  real SI[N2]; // fixed SI using empirical data
#include /data/NKX.stan
#include /data/data_glm.stan
#include /data/hyperparameters.stan
#include /data/glmer_stuff.stan
#include /data/glmer_stuff2.stan
}

transformed data {
  real aux = not_a_number();
  vector[N2] SI_rev; // SI in reverse order
  vector[N2] f_rev[M]; // f in reversed order
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
#include /tdata/tdata_glm.stan

  for(i in 1:N2)
    SI_rev[i] = SI[N2-i+1];

  for(m in 1:M){
    for(i in 1:N2) {
      f_rev[m, i] = f[N2-i+1,m];
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
  real <lower=0> ifr_noise[M];
}

transformed parameters {
  vector[N] Rt_vec;
  matrix[N2, M] prediction = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
  matrix[N2, M] Rt = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_adj = Rt;
  vector[N] eta;  // linear predictor
  vector[N] R0_vec;
  
#include /tparameters/tparameters_glm.stan
#include /model/make_eta.stan


  {
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

  # Todo: Add branching logic for different link functions.
  Rt_vec = R0_vec * 2 .* inv_logit(-eta);
  {
    int idx = 1;
    for (m in 1:M) {
      Rt[1:NC[m],m] = Rt_vec[idx:(idx+NC[m]-1)];
      idx += NC[m];
    }
  }

  {
    matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
    for (m in 1:M){

      prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
      cumm_sum[2:N0,m] = cumulative_sum(prediction[2:N0,m]);
      Rt_adj[1:N0,m] = Rt[1:N0,m];

      for (i in (N0+1):N2) {
        real convolution = dot_product(sub_col(prediction, 1, m, i-1), tail(SI_rev, i-1));
        cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i-1,m];
        Rt_adj[i,m] = ((pop[m]-cumm_sum[i,m]) / pop[m]) * Rt[i,m];
        prediction[i, m] = fmin(pop[m] - cumm_sum[i,m], prediction[i, m] + Rt_adj[i,m] * convolution);
      }

      E_deaths[1, m]= 1e-15 * prediction[1,m];
      for (i in 2:N2){
        E_deaths[i,m] = ifr_noise[m] * dot_product(sub_col(prediction, 1, m, i-1), tail(f_rev[m], i-1));
      }
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
  ifr_noise ~ normal(1,0.1);

#include /model/priors_glm.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau,
                          regularization, delta, shape, t, p);
  }

  if (prior_PD == 0) {
    for(m in 1:M){
      for (i in 1:NC[m]) {
        if (deaths[i,m] != -1) {
          deaths[i,m] ~ neg_binomial_2(E_deaths[i,m], phi);
        }
      }
    }
  }
}

generated quantities {
  real alpha[has_intercept];
  
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  }
  
}

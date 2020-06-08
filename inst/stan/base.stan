functions {
#include /functions/reverse.stan
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan

  vector test_csr_matrix_times_vector(int m, int n, vector w,
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b);
  }
}

data {
#include /data/data_obs.stan
#include /data/data_indices.stan
#include /data/data_model.stan
matrix<lower=0,upper=1>[M, R] means; // mean values for observed events / total cases (for example, IFR)
vector<lower=0>[R] noise_scales;
#include /data/NKX.stan
#include /data/data_glm.stan
#include /data/hyperparameters.stan
#include /data/glmer_stuff.stan
#include /data/glmer_stuff2.stan
}

transformed data {
  real aux = not_a_number();
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
#include /tdata/tdata_reverse.stan
#include /tdata/tdata_glm.stan

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);
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
  vector[N_obs] E_obs; // expected values of the observations 
  vector[N] eta;  // linear predictor

#include /tparameters/infections_rt.stan
#include /tparameters/tparameters_glm.stan
#include /model/make_eta.stan

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

#include /tparameters/gen_infections.stan

  {  // compute expected values of the observations
    for (i in 1:N_obs) {
      int m = obs_group[i];
      int dt = obs_date[i];
      int tp = obs_type[i];
      int n0 = starts[m];
      if (dt == 1)
        E_obs[i] = 1e-15 * infections[1,m];
      else
        E_obs[i] = noise[m, tp] * means[m, tp] * dot_product(sub_col(infections, n0, m, dt-n0), tail(pvecs_rev[tp], dt-n0));
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
    noise[,r] ~ normal(1, noise_scales[r]);

#include /model/priors_glm.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau,
                          regularization, delta, shape, t, p);
  }

  if (prior_PD == 0) {
    for (i in 1:N_obs)
      obs[i] ~ neg_binomial_2(E_obs[i] + 1e-15, phi);
  }
}

generated quantities {
  real alpha[has_intercept];
  
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  }
}

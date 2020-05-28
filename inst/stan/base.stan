functions {
   /*
   * Return matrix formed by data vector.
   *
   * @param NC Number of observations in each group
   * @param N2 Total number of observed days + # of days to forecast
   * @param M Number of groups
   * @return A real-valued matrix parsed from `vec'
   */
  matrix vec_to_mat(int[] NC, int N2, int M, vector vec) {
    matrix[N2, M] mat = rep_matrix(0,N2,M);
    mat[1:NC[1],1] = vec[1:NC[1]];
    for (m in 2:M) {
       mat[m:NC[m],m] = vec[(NC[m-1]+1):NC[m]];
    }
    return mat;
  }

   /*
   * Returns vector of R0 for each observation
   *
   * @param NC Number of observations in each group
   * @param len Total number of observations
   * @param M Number of groups
   * @param mu Current estimate of R0 for each group
   * @return A real-valued vector
   */
  vector r0(int[] NC, int len, int M, real[] mu) {
    vector[len] vec;
      vec[1:NC[1]] = rep_vector(mu[1],NC[1]);
    for (m in 2:M) {
      vec[(NC[m-1]+1):NC[m]] = rep_vector(mu[m], NC[m] - NC[m-1]);
    }
    return vec;
  }

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
  int EpidemicStart[M];
  real pop[M];
  real SI[N2]; // fixed SI using empirical data
#include /data/NKX.stan
#include /data/data_glm.stan
#include /data/weights_offset.stan
#include /data/hyperparameters.stan
#include /data/glmer_stuff.stan
#include /data/glmer_stuff2.stan
}

transformed data {
  int<lower=1> len = sum(NC); //total data observations
  vector[N2] SI_rev; // SI in reverse order
  vector[N2] f_rev[M]; // f in reversed order
  int<lower=1> V[special_case ? t : 0, len] = make_V(len, special_case ? t : 0, v);
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
  real<lower=0> aux_unscaled;
  real<lower=0> mu[M];
  real<lower=0> kappa;
  real<lower=0> y[M];
  real<lower=0> phi;
  real<lower=0> tau2;
  real <lower=0> ifr_noise[M];
}

transformed parameters {
  vector[len] Rt_vec = rep_vector(0, len);
  matrix[N2, M] prediction = rep_matrix(0,N2,M);
  matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
  matrix[N2, M] Rt = rep_matrix(0,N2,M);
  matrix[N2, M] Rt_adj = Rt;
  vector[N] eta;  // linear predictor

  // aux has to be defined first in the hs case
  real aux = prior_dist_for_aux == 0 ? aux_unscaled : (prior_dist_for_aux <= 2 ?
             prior_scale_for_aux * aux_unscaled + prior_mean_for_aux :
             prior_scale_for_aux * aux_unscaled);

#include /tparameters/tparameters_glm.stan
#include /model/make_eta.stan

  if (prior_dist_for_aux == 0) // none
    aux = aux_unscaled;
  else {
    aux = prior_scale_for_aux * aux_unscaled;
    if (prior_dist_for_aux <= 2) // normal or student_t
      aux += prior_mean_for_aux;
  }

  if (t > 0) {
    if (special_case == 1) {
      int start = 1;
      theta_L = scale .* tau * aux;
      if (t == 1) b = theta_L[1] * z_b;
      else for (i in 1:t) {
        int end = start + l[i] - 1;
        b[start:end] = theta_L[i] * z_b[start:end];
        start = end + 1;
      }
    }
    else {
      theta_L = make_theta_L(len_theta_L, p,
                             aux, tau, scale, zeta, rho, z_T);
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
  # Todo: Add branching logic for different weights.
  Rt_vec = r0(NC, len, M, mu) * 2 .* inv_logit(-eta);
  Rt = vec_to_mat(NC, N2, M, Rt_vec);

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
        prediction[i, m] = prediction[i, m] + Rt_adj[i,m] * convolution;
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

  // Log-priors
  if (prior_dist_for_aux > 0 && prior_scale_for_aux > 0) {
    real log_half = -0.693147180559945286;
    if (prior_dist_for_aux == 1)
      target += normal_lpdf(aux_unscaled | 0, 1) - log_half;
    else if (prior_dist_for_aux == 2)
      target += student_t_lpdf(aux_unscaled | prior_df_for_aux, 0, 1) - log_half;
    else
     target += exponential_lpdf(aux_unscaled | 1);
  }

#include /model/priors_glm.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau,
                          regularization, delta, shape, t, p);
  }

  if (prior_PD == 0) {
    for(m in 1:M){
      deaths[EpidemicStart[m]:NC[m], m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:NC[m], m], phi);
    }
  }
}
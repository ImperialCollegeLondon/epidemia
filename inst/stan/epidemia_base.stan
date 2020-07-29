functions {
#include /functions/reverse.stan
#include /functions/common_functions.stan
#include /functions/continuous_likelihoods.stan
#include /functions/linkinv.stan

  vector test_csr_matrix_times_vector(int m, int n, vector w,
                                      int[] v, int[] u, vector b) {
    return csr_matrix_times_vector(m, n, w, v, u, b);
  }
}

data {
#include /data/data_indices.stan
#include /data/data_obs.stan
int obs[N_obs]; // vector of observations
#include /data/data_obs_mm.stan
#include /data/data_model.stan
#include /data/NKX.stan
#include /data/data_glm.stan
#include /data/data_ac.stan
#include /data/hyperparameters.stan
#include /data/glmer_stuff.stan
#include /data/glmer_stuff2.stan
}

transformed data {
  real aux = not_a_number();
  int<lower=1> V[special_case ? t : 0, N] = make_V(N, special_case ? t : 0, v);
  int<lower=1> ac_V[ac_nterms, N] = make_V(N, ac_nterms, ac_v);
#include /tdata/tdata_reverse.stan
#include /tdata/tdata_glm.stan

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);

}

parameters {
  vector[num_ointercepts] ogamma;
  real gamma[has_intercept];
  vector<lower=0>[num_oaux] oaux_raw;
#include /parameters/parameters_glm.stan
#include /parameters/parameters_ac.stan
#include /parameters/parameters_obs.stan
  vector<lower=0>[M] y_raw;
  real<lower=0> tau_raw;
}

transformed parameters {
  vector[N_obs] oeta;
  vector[N_obs] E_obs; // expected values of the observations 
  vector[N] eta;  // linear predictor
  real<lower=0> tau2 = prior_scale_for_tau * tau_raw;
  vector<lower=0>[M] y = tau2 * y_raw;
  vector<lower=0>[num_oaux] oaux = oaux_raw;

#include /tparameters/infections_rt.stan
#include /tparameters/tparameters_ac.stan
#include /tparameters/tparameters_obs.stan
#include /tparameters/tparameters_glm.stan

// transform auxiliary parameters
  for (i in 1:num_oaux) {
    if (prior_dist_for_oaux[i] > 0) {
      if (prior_scale_for_oaux[i] > 0) {
        oaux[i] *= prior_scale_for_oaux[i];
      }
      if (prior_dist_for_oaux[i] <= 2) {
        oaux[i] += prior_mean_for_oaux[i];
      }
    }
  }

  {
    int i = 1;
    for (proc in 1:ac_nproc) { // this treats ac terms as random walks for now (to be extended to AR(p))
        ac_beta[i:(i+ac_ntime[proc]-1)] = cumulative_sum(ac_noise[i:(i+ac_ntime[proc]-1)]);
        i += ac_ntime[proc];
    }
  }
  
#include /tparameters/make_eta.stan
#include /tparameters/make_oeta.stan
#include /tparameters/gen_infections.stan
#include /tparameters/gen_eobs.stan
}

model {
  target += exponential_lpdf(tau_raw | 1);
  target += exponential_lpdf(y_raw | 1);

#include /model/priors_glm.stan
#include /model/priors_ac.stan
#include /model/priors_obs.stan
  if (t > 0) {
    real dummy = decov_lp(z_b, z_T, rho, zeta, tau,
                          regularization, delta, shape, t, p);
  }

  if (prior_PD == 0) {
    int i = 1;
    for (r in 1:R) {
      if (ofamily[r] == 1) { // poisson
        target += poisson_lpmf(segment(obs, i, oN[r]) | segment(E_obs, i, oN[r]));
      }
      else { // neg binom
        target += neg_binomial_2_lpmf(segment(obs, i, oN[r]) | 
        segment(E_obs, i, oN[r]), oaux[has_oaux[r]]);
      }
      i += oN[r];
    }
  }

}

generated quantities {
  real alpha[has_intercept];
  
  if (has_intercept == 1) {
    alpha[1] = gamma[1] - dot_product(xbar, beta);
  }
}

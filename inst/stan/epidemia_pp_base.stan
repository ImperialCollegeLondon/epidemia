functions {
#include functions/reverse.stan
#include functions/linkinv.stan
}

data {
#include data/data_indices.stan
#include data/data_obs.stan
#include data/data_model.stan
}

transformed data {
#include tdata/tdata_reverse.stan

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);
}

parameters {
  vector<lower=0>[M+1] y;
  vector<lower=0>[num_oaux+1] oaux;
  vector[N] eta;
  vector[N_obs] oeta;
}

generated quantities {
  vector[N_obs] E_obs;
  int obs[N_obs];
#include /tparameters/infections_rt.stan
#include /tparameters/gen_infections.stan
#include /tparameters/gen_eobs.stan

  {
    int i = 1;
    for (r in 1:R) {
      if (ofamily[r] == 1) { // poisson
        obs[i:(i+oN[r]-1)] = poisson_rng(linkinv(segment(E_obs, i, oN[r]) + 1e-15, olink[r]));
      }
      else { // neg binom
        obs[i:(i+oN[r]-1)] = neg_binomial_2_rng(linkinv(segment(E_obs, i, oN[r]) + 1e-15, olink[r]), 
          oaux[has_oaux[r]]);
      }
      i += oN[r];
    }
  }
}


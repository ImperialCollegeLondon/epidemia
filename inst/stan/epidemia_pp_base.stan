functions {
#include functions/reverse.stan
#include functions/linkinv.stan
}

data {
#include data/data_indices.stan
#include data/data_obs.stan
#include data/data_model.stan
#include /data/data_inf.stan
}

transformed data {
#include tdata/tdata_reverse.stan

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);
}

parameters {
  vector<lower=0>[M+2] y;
  vector<lower=0>[num_oaux+2] oaux;
  vector[N+2] eta;
  vector[N_obs+2] oeta;
  vector[latent ? N - M * N0 + 2 : 2] inf_noise;
  vector<lower=0>[latent ? 3 : 2] inf_aux;
}

generated quantities {
  vector[N_obs] E_obs;
#include /tparameters/infections_rt.stan
#include /tparameters/gen_infections.stan
#include /tparameters/gen_eobs.stan
}


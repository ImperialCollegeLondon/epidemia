functions {
#include functions/reverse.stan
}

data {
#include data/data_indices.stan
int<lower=1> N; 
int<lower=0> N_obs; // total size of the observation vector
int obs_group[N_obs]; // group (1 to M) to which each observation belongs
int obs_type[N_obs]; // type of observation (1 to r). 
int obs_date[N_obs]; // observation date (1 to N2)
real<lower=0> r0; // the prior expected value for r0
vector<lower=0,upper=1>[NS] pvecs[R]; // the 'pvec' for each type of observation.
vector<lower=0>[M] pop;
simplex[NS] si; // fixed serial interval using empirical data
}

transformed data {
#include tdata/tdata_reverse.stan

for(r in 1:R)
      pvecs_rev[r] = reverse(pvecs[r]);
}

parameters {
  vector<lower=0>[M+1] y;
  real<lower=0> phi[R+1];
  vector[N] eta;
  vector[N_obs] oeta;
}

generated quantities {
  matrix[N2, M] pred[R];
  vector[N_obs] E_obs;
  vector[N_obs] obs;
#include /tparameters/infections_rt.stan

  // initialise to 0
  for (r in 1:R)
    pred[r] = rep_matrix(0, N2, M);

#include /tparameters/gen_infections.stan
#include /tparameters/gen_eobs.stan

    obs = neg_binomial_2_rng(E_obs + 1e-15, phi[obs_type]);
}


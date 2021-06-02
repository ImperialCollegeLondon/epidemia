functions {
#include functions/reverse.stan
#include functions/linkinv.stan
#include functions/trunc_normal.stan
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
      pvecs_rev[r] = reverse2(pvecs[r]);
}

parameters {
  vector<lower=0>[M+2] seeds;
  vector<lower=0>[num_oaux+2] oaux;
  vector[N+2] eta;
  vector[N_obs+2] oeta;
  vector<lower=-1>[latent ? N - M * N0 + 2 : 2] infections_raw;
  vector<lower=0>[latent ? 3 : 2] inf_aux;
  vector<lower=0,upper=1>[S0_fixed ? 2 : M+2] S0;
  vector<lower=0,upper=1>[veps_fixed ? 2 : M+2] veps;
}
generated quantities {
  vector[N_obs] E_obs;
  matrix<lower=0>[N2, M] Rt = rep_matrix(0,N2,M);
  matrix<lower=0>[N2, M] infectiousness;
#include /tparameters/infections_rt.stan
#include /generated_quantities/gen_infections_pp.stan
#include /tparameters/gen_eobs.stan

  if (pop_adjust) {
    for (m in 1:M) {
      Rt[1:N2, m] = susc[1:N2, m] ./ pops[m] .* Rt_unadj[1:N2,m];
    }
  } else {
    Rt = Rt_unadj;
  }

  infectiousness = load / max(gen);
}


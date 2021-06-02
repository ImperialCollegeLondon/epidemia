
vector<lower=0>[hseeds] seeds_aux_raw;
vector<lower=0>[M] seeds_raw;
vector<lower=0>[latent ? N - M * N0 : 0] infections_raw;
vector<lower=0,upper=1>[S0_fixed ? 0 : M] S0;
vector<lower=0, upper=1>[veps_fixed ? 0 : M] veps;
